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

# openfe minimums enforced here so that a short run errors out neither silently nor fatally.
MIN_MBAR_SAMPLES = 100        # fewer samples/state and pymbar spins to its 10000-iteration ceiling
MIN_EQUIL_TRAJ_FRAMES = 10    # fewer frames and the ABFE complex leg dies on an empty .xtc

def shortenInterval(obj, attr, prodLength, minCount):
  """Shorten obj.attr so that prodLength still yields at least minCount of whatever it counts.
  Only ever shortens, so it is a no-op at openfe's own defaults."""
  from openff.units import unit
  prodPs = prodLength.m_as(unit.picosecond)
  intervalPs = getattr(obj, attr).m_as(unit.picosecond)
  if intervalPs > 0 and prodPs / intervalPs < minCount:
    setattr(obj, attr, (prodPs / minCount) * unit.picosecond)

def ensureEnoughSamples(simSettings, minSamples=MIN_MBAR_SAMPLES):
  """Too few MBAR samples does not merely cost precision: pymbar never converges and burns CPU
  for tens of minutes, once per bootstrap resample."""
  shortenInterval(simSettings, 'time_per_iteration', simSettings.production_length, minSamples)
  return simSettings

def ensureTrajectoryFrames(simSettings, outputSettings, minFrames=MIN_EQUIL_TRAJ_FRAMES):
  """The ABFE complex leg picks its Boresch anchors from the RMSF over the pre-equilibration
  trajectory, so a 0-frame .xtc is fatal ("XDR read error = endoffile"), not cosmetic."""
  shortenInterval(outputSettings, 'trajectory_write_interval',
                  simSettings.production_length, minFrames)
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
