# Utils for the scripts
import os

def parseParams(paramsFile, listParams=[], sep=':'):
  paramsDic = {}
  with open(paramsFile) as f:
    for line in f:
      key, value = line.strip().split(sep)
      if key in listParams:
        paramsDic[key.strip()] = value.strip().split()
      else:
        paramsDic[key.strip()] = value.strip()
  return paramsDic

def parseMoleculeFile(molFile, sanitize=True):
  from rdkit import Chem
  if molFile.endswith('.mol2'):
    mol = Chem.MolFromMol2File(molFile)
  elif molFile.endswith('.mol'):
    mol = Chem.MolFromMolFile(molFile)
  elif molFile.endswith('.pdb'):
    mol = Chem.MolFromPDBFile(molFile)
  elif molFile.endswith('.smi'):
    with open(molFile, "r") as f:
      line = f.readline()
      if line.startswith('SMILES'):
        line = f.readline()
      mol = Chem.MolFromSmiles(line)
  elif molFile.endswith('.sdf'):
    suppl = Chem.SDMolSupplier(molFile, sanitize=sanitize)
    for mol in suppl:
      break
  else:
    mol = Chem.MolFromSmiles(molFile)

  return mol

def getMolFilesDic(molFiles, sanitize=True):
  molsDict = {}
  for molFile in molFiles:
    m = parseMoleculeFile(molFile, sanitize=sanitize)
    if m:
      molsDict[m] = molFile

  mols = list(molsDict.keys())
  return molsDict, mols

def writeMol(mol, outFile, cid=-1, setName=False):
  from rdkit import Chem
  w = Chem.SDWriter(outFile)
  molName = os.path.split(os.path.splitext(outFile)[0])[-1]
  if setName:
    mol.SetProp('_Name', molName)
  w.write(mol, cid)
  w.close()

def getBaseName(file):
  return os.path.splitext(os.path.basename(file.strip()))[0]

MIN_MBAR_SAMPLES = 100

def ensureEnoughSamples(simSettings, minSamples=MIN_MBAR_SAMPLES):
  '''Shorten an openfe MultiStateSimulationSettings' time_per_iteration so the run still yields
  at least `minSamples` energy samples per lambda state.

  MBAR needs a decent number of samples; with too few it does not merely give a poor estimate,
  it HANGS. Confirmed the hard way on a real run: production_length=0.01 ns against openfe's
  default time_per_iteration of 2.5 ps gives 4 samples per state, and pymbar then span at 100%
  CPU on one core for 35+ minutes, emitting "Did not converge to within specified tolerance,
  max_delta = 0.000000e+00, iterations completed = 9999" 118 times over - once per bootstrap
  resample - because a degenerate solve never satisfies the 1e-12 tolerance.

  Only ever shortens: at openfe's production defaults (5 ns / 2.5 ps = 2000 samples) this is a
  no-op, so it costs nothing in real runs and only rescues short ones.'''
  from openff.units import unit
  prodPs = simSettings.production_length.m_as(unit.picosecond)
  perIterPs = simSettings.time_per_iteration.m_as(unit.picosecond)
  if perIterPs > 0 and prodPs / perIterPs < minSamples:
    simSettings.time_per_iteration = (prodPs / minSamples) * unit.picosecond
  return simSettings

MIN_EQUIL_TRAJ_FRAMES = 10

def ensureTrajectoryFrames(simSettings, outputSettings, minFrames=MIN_EQUIL_TRAJ_FRAMES):
  '''Shorten an openfe pre-equilibration output block's trajectory_write_interval so the run
  actually writes at least `minFrames` frames.

  openfe's ABFE complex leg picks its Boresch restraint anchor atoms from the RMSF measured over
  this pre-equilibration trajectory, so an empty trajectory is a hard failure, not a cosmetic
  one. Confirmed the hard way: with a 10 ps pre-equilibration production against the default
  20 ps trajectory_write_interval, NO frames are written - production_equil.xtc ends up 0 bytes -
  and the complex Setup unit dies with "OSError: XDR read error = endoffile" (retried 3 times,
  identically). The solvent leg survives the same empty file because it never reads it back.

  Only ever shortens the interval: at openfe's own pre-equilibration lengths (5 ns complex /
  0.5 ns solvent production vs a 20 ps interval = 250 / 25 frames) this is a no-op.'''
  from openff.units import unit
  prodPs = simSettings.production_length.m_as(unit.picosecond)
  intervalPs = outputSettings.trajectory_write_interval.m_as(unit.picosecond)
  if intervalPs > 0 and prodPs / intervalPs < minFrames:
    outputSettings.trajectory_write_interval = (prodPs / minFrames) * unit.picosecond
  return outputSettings

def getGenerator(ligFF):
  from openmmforcefields.generators import EspalomaTemplateGenerator, GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
  if 'espaloma' in ligFF.lower():
    gen = EspalomaTemplateGenerator
  elif 'gaff' in ligFF.lower():
    gen = GAFFTemplateGenerator
  elif 'smirnoff' in ligFF.lower() or 'openff' in ligFF.lower():
    gen = SMIRNOFFTemplateGenerator
  return gen

def addMoleculesFF(forcefield, ligFile, ligFF):
  '''Update forcefiled with Espaloma parameters for ligand'''
  from openff.toolkit.topology import Molecule
  molecule = Molecule.from_file(ligFile)
  generator = getGenerator(ligFF)
  tempGenerator = generator(molecules=molecule, forcefield=ligFF, cache="molecules_ff.json")
  forcefield.registerTemplateGenerator(tempGenerator.generator)
  return forcefield
