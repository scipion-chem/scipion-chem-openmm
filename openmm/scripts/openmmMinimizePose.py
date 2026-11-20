# **************************************************************************
# *
# * Authors:  Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
Minimize ONE ligand pose together with the side chains of the surrounding pocket.

 - Receptor: protein only (Amber force field), hydrogens added, implicit solvent.
 - Ligand: parametrized from the SDF (GAFF / SMIRNOFF / Espaloma through openmmforcefields).
 - Only the ligand and the side chains of residues within `pocketRadius` of it are free to move.
   The rest of the protein is frozen (mass 0).
 - The ligand heavy atoms are harmonically restrained to their initial positions with a
   force constant that is released stage by stage (`stages`, in kcal/mol/A^2).
 - Optionally the ligand coordinates are randomly perturbed (`perturbSigma`, A) before minimizing,
   e.g. to let perfectly flat saturated rings pucker.
 - No validation is done here: check the result with a PoseBusters protocol and run again if needed.

Usage: python openmmMinimizePose.py params.txt      (key :: value lines)
Writes (in the working dir): <name>_min.sdf, <name>_min_receptor.pdb, <name>_min_complex.pdb, result.json
"""

import sys, os, json, shutil, traceback
import numpy as np

from openmm.app import *
from openmm import *
from openmm.unit import nanometer, angstrom, dalton, picoseconds, kilojoules_per_mole

from utils import parseParams, getGenerator

PROT_RESNAMES = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL '
                    'HID HIE HIP HIN CYX CYM ASH GLH LYN ACE NME NHE'.split())
BACKBONE = {'N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'HA2', 'HA3', 'C', 'O', 'OXT', 'HXT'}
KCAL_A2_TO_KJ_NM2 = 418.4


########################################################################################
# Loading
########################################################################################
def loadReceptor(recFile):
  """Protein-only Modeller (no waters/ligands/hetero), missing heavy atoms added if pdbfixer is available"""
  try:
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=recFile)
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # never build loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    top, pos = fixer.topology, fixer.positions
  except ImportError:
    pdb = PDBFile(recFile)
    top, pos = pdb.topology, pdb.positions

  modeller = Modeller(top, pos)
  modeller.delete([r for r in modeller.topology.residues() if r.name not in PROT_RESNAMES])
  return modeller


def loadLigand(ligFile, charges='default'):
  """Returns the RDKit molecule (with explicit Hs and the input coordinates) and the openff Molecule"""
  from rdkit import Chem
  from openff.toolkit.topology import Molecule

  suppl = Chem.SDMolSupplier(ligFile, removeHs=False)
  mol = suppl[0] if len(suppl) else None
  if mol is None:
    suppl = Chem.SDMolSupplier(ligFile, removeHs=False, sanitize=False)
    mol = suppl[0]
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, catchErrors=True)
  mol = Chem.AddHs(mol, addCoords=True)  # only adds the missing hydrogens

  offMol = Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
  if charges.lower() == 'gasteiger':
    offMol.assign_partial_charges('gasteiger')
  return mol, offMol


def addLigandToModeller(modeller, offMol):
  ligTop = offMol.to_topology().to_openmm()
  for res in ligTop.residues():
    res.name = 'LIG'
  ligPos = offMol.conformers[0].m_as('nanometer') * nanometer
  modeller.add(ligTop, ligPos)
  return modeller


########################################################################################
# Minimization
########################################################################################
def getSelections(top, pos, radiusNm):
  """Ligand atoms, ligand heavy atoms and the free (mobile) atoms: ligand + side chains of pocket residues"""
  atoms = list(top.atoms())
  ligAtoms = [a.index for a in atoms if a.residue.name == 'LIG']
  ligHeavy = [i for i in ligAtoms if atoms[i].element is not None and atoms[i].element.symbol != 'H']
  protHeavy = [a.index for a in atoms if a.residue.name != 'LIG' and a.element is not None
               and a.element.symbol != 'H']

  d = np.linalg.norm(pos[protHeavy][:, None, :] - pos[ligHeavy][None, :, :], axis=2).min(axis=1)
  pocketRes = {atoms[protHeavy[k]].residue.index for k in np.where(d <= radiusNm)[0]}

  mobile = set(ligAtoms)
  for a in atoms:
    if a.residue.index in pocketRes and a.residue.name != 'LIG' and a.name not in BACKBONE:
      mobile.add(a.index)
  return ligAtoms, ligHeavy, sorted(mobile), pocketRes


def minimizeAttempt(modeller, forcefield, pDic, radiusA, sigmaA, seed):
  top = modeller.topology
  pos0 = np.array(modeller.positions.value_in_unit(nanometer))
  nAtoms = len(pos0)

  ligAtoms, ligHeavy, mobile, pocketRes = getSelections(top, pos0, radiusA / 10.)
  system = forcefield.createSystem(top, nonbondedMethod=NoCutoff, constraints=None, removeCMMotion=False)

  mobileSet = set(mobile)
  for i in range(nAtoms):
    if i not in mobileSet:
      system.setParticleMass(i, 0 * dalton)

  rest = CustomExternalForce('0.5*kLig*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
  rest.addGlobalParameter('kLig', 0.0)
  for p in ('x0', 'y0', 'z0'):
    rest.addPerParticleParameter(p)
  for i in ligHeavy:
    rest.addParticle(i, pos0[i])
  system.addForce(rest)

  integrator = VerletIntegrator(0.001 * picoseconds)
  platform = None
  try:
    platform = max((Platform.getPlatform(i) for i in range(Platform.getNumPlatforms())), key=lambda p: p.getSpeed())
    if platform.getName() == 'CPU':
      platform.setPropertyDefaultValue('Threads', str(pDic.get('cpuThreads', 1)))
  except Exception:
    platform = None
  context = Context(system, integrator, platform) if platform else Context(system, integrator)

  pos = pos0.copy()
  rng = np.random.RandomState(seed)
  if sigmaA > 0:  # restraints reference stays at the original pose
    pos[ligAtoms] += rng.normal(0, sigmaA / 10., size=(len(ligAtoms), 3))
  context.setPositions(pos)

  e0 = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
  tol = float(pDic.get('tolerance', 10)) * kilojoules_per_mole / nanometer
  maxIter = int(pDic.get('maxIter', 2000))
  for k in [float(x) for x in pDic.get('stages', '10 2 0.5 0').replace(',', ' ').split()]:
    context.setParameter('kLig', k * KCAL_A2_TO_KJ_NM2)
    LocalEnergyMinimizer.minimize(context, tol, maxIter)

  state = context.getState(getEnergy=True, getPositions=True)
  posF = np.array(state.getPositions(asNumpy=True).value_in_unit(nanometer))
  e1 = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
  rmsd = float(np.sqrt(((posF[ligHeavy] - pos0[ligHeavy]) ** 2).sum(axis=1).mean()) * 10.)  # Angstrom
  return {'pos': posF, 'rmsd': rmsd, 'e0': e0, 'e1': e1, 'nMobile': len(mobile), 'nPocketRes': len(pocketRes),
          'radius': radiusA, 'sigma': sigmaA}


def getRDKitConnectivity(mol):
  """Return indexed RDKit bonds, including bond order and aromaticity."""
  return sorted(
    (
      min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
      max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
      str(b.GetBondType()),
      bool(b.GetIsAromatic())
    )
    for b in mol.GetBonds()
  )


def getRDKitAdjacency(mol):
  """Return the atom-index pairs connected by bonds in an RDKit molecule."""
  return {
    (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
     max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
    for b in mol.GetBonds()
  }


def getOpenMMLigandAdjacency(topology, ligAtoms):
  """Return ligand bonds using ligand-local indices in OpenMM atom order."""
  globalToLocal = {globalIdx: localIdx for localIdx, globalIdx in enumerate(ligAtoms)}
  adjacency = set()

  for bond in topology.bonds():
    i, j = bond.atom1.index, bond.atom2.index
    if i in globalToLocal and j in globalToLocal:
      a, b = globalToLocal[i], globalToLocal[j]
      adjacency.add((min(a, b), max(a, b)))

  return adjacency


def getCIPAssignments(mol):
  """Return assigned and unassigned tetrahedral stereocentres by atom index."""
  from rdkit import Chem
  Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
  return tuple(Chem.FindMolChiralCenters(
    mol, force=True, includeUnassigned=True, includeCIP=True
  ))


def validateLigandTopology(rdMol, modeller, ligAtoms):
  """Check atom count, atom identity/order and connectivity before minimization output."""
  topAtoms = list(modeller.topology.atoms())

  if len(ligAtoms) != rdMol.GetNumAtoms():
    raise ValueError(
      f"Atom count mismatch: OpenMM={len(ligAtoms)}, RDKit={rdMol.GetNumAtoms()}"
    )

  for k, globalIdx in enumerate(ligAtoms):
    rdAtom = rdMol.GetAtomWithIdx(k)
    ommAtom = topAtoms[globalIdx]

    if ommAtom.element is None:
      raise ValueError(f"OpenMM atom {globalIdx} has no element assigned")

    rdElement = rdAtom.GetSymbol()
    ommElement = ommAtom.element.symbol
    if rdElement != ommElement:
      raise ValueError(
        f"Atom mapping mismatch at ligand atom {k}: "
        f"RDKit={rdElement}, OpenMM={ommElement} (global index {globalIdx})"
      )

  rdAdjacency = getRDKitAdjacency(rdMol)
  ommAdjacency = getOpenMMLigandAdjacency(modeller.topology, ligAtoms)

  if rdAdjacency != ommAdjacency:
    missing = sorted(rdAdjacency - ommAdjacency)
    extra = sorted(ommAdjacency - rdAdjacency)
    raise ValueError(
      "RDKit/OpenMM ligand connectivity mismatch. "
      f"Bonds missing from OpenMM (RDKit local atom indices): {missing}; "
      f"extra OpenMM bonds: {extra}"
    )

  # Explicitly report hydrogen-neighbour counts as a more readable diagnostic.
  for k, rdAtom in enumerate(rdMol.GetAtoms()):
    rdH = sum(n.GetSymbol() == 'H' for n in rdAtom.GetNeighbors())
    ommH = sum(
      1 for a, b in ommAdjacency
      if (a == k and topAtoms[ligAtoms[b]].element.symbol == 'H') or
         (b == k and topAtoms[ligAtoms[a]].element.symbol == 'H')
    )
    if rdH != ommH:
      raise ValueError(
        f"Hydrogen-neighbour mismatch at ligand atom {k}: "
        f"RDKit={rdH}, OpenMM={ommH}"
      )

  print(
    f"Ligand topology validation passed: {rdMol.GetNumAtoms()} atoms, "
    f"{len(rdAdjacency)} bonds; RDKit and OpenMM atom order/connectivity agree.",
    flush=True
  )


def validateWrittenLigand(inputMol, outputFile):
  """Check that the written SDF retains the input RDKit graph and stereochemistry."""
  from rdkit import Chem

  suppl = Chem.SDMolSupplier(outputFile, removeHs=False, sanitize=True)
  outputMol = suppl[0] if len(suppl) else None
  if outputMol is None:
    raise ValueError(f"Could not read minimized ligand SDF: {outputFile}")

  if inputMol.GetNumAtoms() != outputMol.GetNumAtoms():
    raise ValueError(
      f"Output SDF atom count changed: input={inputMol.GetNumAtoms()}, "
      f"output={outputMol.GetNumAtoms()}"
    )

  inputConnectivity = getRDKitConnectivity(inputMol)
  outputConnectivity = getRDKitConnectivity(outputMol)
  if inputConnectivity != outputConnectivity:
    raise ValueError(
      "Ligand connectivity changed between input and output SDF. "
      f"Input bonds: {inputConnectivity}; output bonds: {outputConnectivity}"
    )

  inputStereo = getCIPAssignments(Chem.Mol(inputMol))
  outputStereo = getCIPAssignments(Chem.Mol(outputMol))
  if inputStereo != outputStereo:
    raise ValueError(
      "Ligand stereochemistry changed between input and output SDF. "
      f"Input stereocentres: {inputStereo}; output stereocentres: {outputStereo}"
    )

  print(
    f"Output SDF validation passed: {outputMol.GetNumAtoms()} atoms, "
    f"{len(outputConnectivity)} bonds; bond orders, aromaticity and CIP "
    "assignments are preserved.",
    flush=True
  )


########################################################################################
# Writing and checking
########################################################################################
def writeLigandSdf(rdMol, ligAtoms, pos, outFile):
  from rdkit import Chem
  from rdkit.Geometry import Point3D
  m = Chem.Mol(rdMol)
  conf = m.GetConformer()
  for k, i in enumerate(ligAtoms):
    x, y, z = pos[i] * 10.
    conf.SetAtomPosition(k, Point3D(float(x), float(y), float(z)))
  w = Chem.SDWriter(outFile)
  w.write(m)
  w.close()


def writeProtein(modeller, pos, outFile, withLigand=False):
  m = Modeller(modeller.topology, pos * nanometer)
  if not withLigand:
    m.delete([r for r in m.topology.residues() if r.name == 'LIG'])
  with open(outFile, 'w') as f:
    PDBFile.writeFile(m.topology, m.positions, f, keepIds=True)


########################################################################################
if __name__ == "__main__":
  pDic = parseParams(sys.argv[1], sep='::')
  name = pDic['molName']
  recFile, ligFile = os.path.abspath(pDic['receptorFile']), os.path.abspath(pDic['ligandFile'])
  outLig, outRec, outCplx = f'{name}_min.sdf', f'{name}_min_receptor.pdb', f'{name}_min_complex.pdb'
  result = {'name': name, 'status': 'error', 'message': ''}

  try:
    forcefield = ForceField(pDic['mFF'], *([pDic['implicitFF']] if pDic.get('implicitFF', 'None') != 'None' else []))
    modeller = loadReceptor(recFile)
    modeller.addHydrogens(forcefield, pH=float(pDic.get('hPH', 7.0)))

    rdMol, offMol = loadLigand(ligFile, pDic.get('ligCharges', 'default'))
    gen = getGenerator(pDic['ligandFF'])(molecules=offMol, forcefield=pDic['ligandFF'], cache='molecules_ff.json')
    forcefield.registerTemplateGenerator(gen.generator)
    modeller = addLigandToModeller(modeller, offMol)
    ligAtoms = [a.index for a in modeller.topology.atoms() if a.residue.name == 'LIG']

    res = minimizeAttempt(modeller, forcefield, pDic, float(pDic.get('pocketRadius', 8)),
                          float(pDic.get('perturbSigma', 0)), seed=int(pDic.get('seed', 0)))
    print(f"E {res['e0']:.0f} -> {res['e1']:.0f} kJ/mol, ligand RMSD {res['rmsd']:.2f} A, "
          f"{res['nPocketRes']} pocket residues", flush=True)

    # Validate that OpenMM retained the exact atom ordering and bond graph
    # expected by the original RDKit molecule before transferring coordinates.
    validateLigandTopology(rdMol, modeller, ligAtoms)

    writeLigandSdf(rdMol, ligAtoms, res['pos'], outLig)
    validateWrittenLigand(rdMol, outLig)

    writeProtein(modeller, res['pos'], outRec)
    writeProtein(modeller, res['pos'], outCplx, withLigand=True)
    result.update({'status': 'minimized', 'rmsd': res['rmsd'], 'e0': res['e0'], 'e1': res['e1'],
                   'nPocketRes': res['nPocketRes']})
  except Exception as e:
    result.update({'status': 'error', 'message': f'{e}'})
    traceback.print_exc()
    try:  # keep the input so the output set is complete
      shutil.copy(ligFile, outLig)
      shutil.copy(recFile, outRec)
    except Exception:
      pass

  with open('result.json', 'w') as f:
    json.dump(result, f, indent=1)
  print(json.dumps(result), flush=True)
