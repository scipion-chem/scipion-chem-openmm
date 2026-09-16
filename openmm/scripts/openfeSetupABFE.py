#Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# -*- coding: utf-8 -*-
"""
Sets up an Open Free Energy absolute binding free energy (ABFE) calculation: charges the input
ligands and writes one transformation JSON per ligand for `openfe quickrun` to execute.

Follows openfe's own ABFE tutorial step for step
(https://docs.openfree.energy/en/stable/tutorials/abfe_tutorial.html): state A is the solvated
complex with the ligand, state B is the same system WITHOUT the ligand - the Protocol derives
the solvent leg from those itself, so no separate solvent state is defined here. Everything not
set from the params file is left at AbsoluteBindingProtocol.default_settings().
"""
import os, sys

import openfe
from openfe.protocols.openmm_afe import AbsoluteBindingProtocol
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
from openff.units import unit
from rdkit import Chem

from utils import parseParams, ensureEnoughSamples, ensureTrajectoryFrames


def loadChargedLigands(sdfFile, chargeMethod):
  supp = Chem.SDMolSupplier(sdfFile, removeHs=False)
  ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp if mol is not None]
  if not ligands:
    raise ValueError(f'No readable ligands found in {sdfFile}')

  backend = 'openeye' if chargeMethod == 'am1bccelf10' else 'ambertools'
  cSettings = OpenFFPartialChargeSettings(partial_charge_method=chargeMethod,
                                         off_toolkit_backend=backend)
  return bulk_assign_partial_charges(
    molecules=ligands, overwrite=False,
    method=cSettings.partial_charge_method,
    toolkit_backend=cSettings.off_toolkit_backend,
    generate_n_conformers=cSettings.number_of_conformers,
    nagl_model=cSettings.nagl_model, processors=1)


def buildSettings(pDic):
  """Apply the exposed settings onto AbsoluteBindingProtocol's defaults.

  This protocol keeps SEPARATE blocks per leg (complex_* / solvent_*), so each gets its own
  padding/equilibration/production. Fields are accessed directly, not through getattr guards: a
  guard would silently keep the default if a name changed upstream.

  The lambda schedules are deliberately untouched - they are parallel LISTS whose length must
  match the leg's n_replicas, so windows are not a safe single-number knob as in RBFE.

  One deviation from openfe: it pads the legs differently (1.0 nm complex, 1.5 nm solvent) while
  the form exposes one value for both. Safe (the floor is 1.3 nm anyway), but it does solvate
  the complex leg more generously, i.e. more water over 30 windows."""
  settings = AbsoluteBindingProtocol.default_settings()
  settings.protocol_repeats = int(pDic['protocolRepeats'])
  settings.forcefield_settings.small_molecule_forcefield = pDic['smallMolFF']
  settings.partial_charge_settings.partial_charge_method = pDic['chargeMethod']
  settings.thermo_settings.temperature = float(pDic['temperature']) * unit.kelvin
  settings.restraint_settings.host_min_distance = float(pDic['hostMinDistance']) * unit.nanometer
  settings.restraint_settings.host_max_distance = float(pDic['hostMaxDistance']) * unit.nanometer
  if 'computePlatform' in pDic:
    settings.engine_settings.compute_platform = pDic['computePlatform']
  if 'gpuIndex' in pDic:
    # Optional[list[int]] - a bare int is rejected by openfe's pydantic validation.
    settings.engine_settings.gpu_device_index = [int(i) for i in pDic['gpuIndex'].split()]

  padding = float(pDic['solventPadding']) * unit.nanometer
  equil = float(pDic['equilLength']) * unit.nanosecond
  prod = float(pDic['productionLength']) * unit.nanosecond
  minSteps = int(pDic['minimizationSteps'])
  # 0/absent means "leave openfe's own per-leg pre-equilibration alone" - see the form help.
  preEquil = float(pDic.get('preEquilLength', 0) or 0)

  for prefix in ('complex', 'solvent'):
    getattr(settings, f'{prefix}_solvation_settings').solvent_padding = padding
    simSettings = getattr(settings, f'{prefix}_simulation_settings')
    simSettings.equilibration_length = equil
    simSettings.production_length = prod
    simSettings.minimization_steps = minSteps
    # See ensureEnoughSamples: too few MBAR samples hangs pymbar rather than just being noisy.
    ensureEnoughSamples(simSettings)

    # Plain-MD pre-equilibration, before any alchemical window. openfe's defaults total 6.55 ns
    # per repeat and no window setting touches them, so a short run must shrink this explicitly.
    equilSettings = getattr(settings, f'{prefix}_equil_simulation_settings')
    equilSettings.minimization_steps = minSteps
    if preEquil > 0:
      equilSettings.equilibration_length_nvt = preEquil * unit.nanosecond
      equilSettings.equilibration_length = preEquil * unit.nanosecond
      equilSettings.production_length = preEquil * unit.nanosecond
      # Boresch anchors come from the RMSF over this trajectory, so a shortened pre-equilibration
      # must write frames more often or it writes none at all.
      ensureTrajectoryFrames(equilSettings, getattr(settings, f'{prefix}_equil_output_settings'))
  return settings


def safeName(name, index):
  """File-name-safe transformation id: docking-derived names may contain spaces or slashes."""
  cleaned = ''.join(c if c.isalnum() or c in '-_' else '_' for c in name)
  return f'lig{index}_{cleaned}' if cleaned else f'lig{index}'


if __name__ == '__main__':
    pDic = parseParams(sys.argv[1], sep='::')
    print(f'openfe version: {openfe.__version__}')

    ligands = loadChargedLigands(pDic['ligandsSdf'], pDic['chargeMethod'])
    protein = openfe.ProteinComponent.from_pdb_file(pDic['proteinPdb'])
    solvent = openfe.SolventComponent()

    protocol = AbsoluteBindingProtocol(settings=buildSettings(pDic))
    transformDir = pDic['transformDir']
    os.makedirs(transformDir, exist_ok=True)

    nameLines = []
    for i, ligand in enumerate(ligands):
      transName = safeName(ligand.name, i)

      # State A: ligand bound in the solvated complex. State B: the same system with the
      # ligand removed - openfe builds the solvent leg from these itself.
      systemA = openfe.ChemicalSystem({'ligand': ligand, 'protein': protein, 'solvent': solvent},
                                      name=ligand.name)
      systemB = openfe.ChemicalSystem({'protein': protein, 'solvent': solvent})

      transformation = openfe.Transformation(stateA=systemA, stateB=systemB, mapping=None,
                                             protocol=protocol, name=transName)
      transformation.dump(os.path.join(transformDir, f'{transName}.json'))
      nameLines.append(f'{transName} :: {ligand.name}')
      print(f'wrote {transName}.json  ({ligand.name})')

    with open(pDic['ligandsFile'], 'w') as f:
      f.write('\n'.join(nameLines) + '\n')
    print(f'prepared {len(nameLines)} ABFE transformation(s)')
