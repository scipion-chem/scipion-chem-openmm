#Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# -*- coding: utf-8 -*-
"""
Sets up an Open Free Energy relative binding free energy (RBFE) campaign: assigns partial
charges once per ligand, plans a minimal-spanning alchemical network with an atom mapper, and
writes one transformation JSON per (edge, leg) for `openfe quickrun` to execute.

Follows openfe's own RBFE python tutorial step for step
(https://docs.openfree.energy/en/stable/tutorials/rbfe_python_tutorial.html). Everything not
set from the params file is left at RelativeHybridTopologyProtocol.default_settings().
"""
import os, sys

import openfe
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
from openff.units import unit
from rdkit import Chem

from utils import parseParams, ensureEnoughSamples

SOLVENT, COMPLEX = 'solvent', 'complex'


def loadChargedLigands(sdfFile, chargeMethod):
  supp = Chem.SDMolSupplier(sdfFile, removeHs=False)
  ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp if mol is not None]
  if len(ligands) < 2:
    raise ValueError(f'Need at least 2 readable ligands in {sdfFile}, got {len(ligands)}')

  # Charge once per ligand and reuse everywhere (the paper's approach: avoids the
  # conformer-dependent charge irreproducibility of charging per-transformation).
  backend = 'openeye' if chargeMethod == 'am1bccelf10' else 'ambertools'
  cSettings = OpenFFPartialChargeSettings(partial_charge_method=chargeMethod,
                                         off_toolkit_backend=backend)
  return bulk_assign_partial_charges(
    molecules=ligands, overwrite=False,
    method=cSettings.partial_charge_method,
    toolkit_backend=cSettings.off_toolkit_backend,
    generate_n_conformers=cSettings.number_of_conformers,
    nagl_model=cSettings.nagl_model, processors=1)


def buildMapper(pDic):
  elementChange = eval(pDic['elementChange'])
  if pDic['mapper'] == 'Kartograf':
    from kartograf import KartografAtomMapper
    return KartografAtomMapper(atom_max_distance=float(pDic['max3d']),
                               map_hydrogens_on_hydrogens_only=not elementChange)
  return openfe.LomapAtomMapper(max3d=float(pDic['max3d']), element_change=elementChange)


def buildSettings(pDic):
  """Apply the handful of exposed settings onto RelativeHybridTopologyProtocol's defaults.

  Attributes are set DIRECTLY, never via hasattr/getattr guards: a guard silently skips a
  misspelled field and leaves the default in place, which is how an earlier version of this
  function shipped a wrong lambda schedule ('lambda_elec_windows' does not exist on the RFE
  LambdaSettings - the field is 'lambda_windows'). Field names below were confirmed by
  introspecting the installed openfe, not read off the docs."""
  settings = RelativeHybridTopologyProtocol.default_settings()
  settings.protocol_repeats = int(pDic['protocolRepeats'])
  settings.forcefield_settings.small_molecule_forcefield = pDic['smallMolFF']
  settings.partial_charge_settings.partial_charge_method = pDic['chargeMethod']
  settings.thermo_settings.temperature = float(pDic['temperature']) * unit.kelvin
  settings.solvation_settings.solvent_padding = float(pDic['solventPadding']) * unit.nanometer
  settings.simulation_settings.minimization_steps = int(pDic['minimizationSteps'])
  settings.simulation_settings.equilibration_length = float(pDic['equilLength']) * unit.nanosecond
  settings.simulation_settings.production_length = float(pDic['productionLength']) * unit.nanosecond
  if 'computePlatform' in pDic:
    settings.engine_settings.compute_platform = pDic['computePlatform']
  if 'gpuIndex' in pDic:
    # Optional[list[int]] - a bare int is rejected by openfe's pydantic validation.
    settings.engine_settings.gpu_device_index = [int(i) for i in pDic['gpuIndex'].split()]

  # These two MUST be set together: openfe's own _validate rejects the protocol outright with
  # "Number of replicas in simulation_settings: N must equal the number of lambda windows in
  # lambda_settings: M" if they disagree.
  nWindows = int(pDic['nReplicas'])
  settings.lambda_settings.lambda_windows = nWindows
  settings.simulation_settings.n_replicas = nWindows

  # Keep enough MBAR samples even for very short productions - without this, pymbar hangs
  # rather than merely losing precision (see ensureEnoughSamples).
  ensureEnoughSamples(settings.simulation_settings)
  return settings


if __name__ == '__main__':
    pDic = parseParams(sys.argv[1], sep='::')
    print(f'openfe version: {openfe.__version__}')

    ligands = loadChargedLigands(pDic['ligandsSdf'], pDic['chargeMethod'])
    protein = openfe.ProteinComponent.from_pdb_file(pDic['proteinPdb'])
    solvent = openfe.SolventComponent()

    network = openfe.ligand_network_planning.generate_minimal_spanning_network(
      ligands=ligands, mappers=[buildMapper(pDic)],
      scorer=openfe.lomap_scorers.default_lomap_score)

    with open(pDic['networkFile'], 'w') as f:
      f.write(network.to_graphml())

    protocol = RelativeHybridTopologyProtocol(buildSettings(pDic))
    transformDir = pDic['transformDir']
    os.makedirs(transformDir, exist_ok=True)

    edgeLines = []
    for i, mapping in enumerate(network.edges):
      nameA, nameB = mapping.componentA.name, mapping.componentB.name
      # Scipion-side file-name-safe edge id; the real ligand names are recorded in edgesFile
      # so the protocol can label results without re-deriving the network.
      edgeName = f'edge{i}'
      edgeLines.append(f'{edgeName} :: {nameA} :: {nameB}')

      for leg in (SOLVENT, COMPLEX):
        sysADict = {'ligand': mapping.componentA, 'solvent': solvent}
        sysBDict = {'ligand': mapping.componentB, 'solvent': solvent}
        if leg == COMPLEX:
          sysADict['protein'] = protein
          sysBDict['protein'] = protein

        sysA = openfe.ChemicalSystem(sysADict, name=f'{nameA}_{leg}')
        sysB = openfe.ChemicalSystem(sysBDict, name=f'{nameB}_{leg}')
        transformation = openfe.Transformation(stateA=sysA, stateB=sysB, mapping=mapping,
                                               protocol=protocol, name=f'{edgeName}_{leg}')
        transformation.dump(os.path.join(transformDir, f'{edgeName}_{leg}.json'))
        print(f'wrote {edgeName}_{leg}.json  ({nameA} -> {nameB}, {leg})')

    with open(pDic['edgesFile'], 'w') as f:
      f.write('\n'.join(edgeLines) + '\n')
    print(f'planned {len(edgeLines)} edge(s), {2 * len(edgeLines)} transformation(s)')
