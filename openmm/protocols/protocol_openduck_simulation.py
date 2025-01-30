# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
This file will run a OpenDuck undocking simulation using OpenMM
"""
import os, json, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC
from pwchem.utils import RESIDUES1TO3, convertToSdf

from .. import Plugin
from ..objects import OpenMMSystem
from ..constants import OPENMM_DIC

program = 'openduck openmm-full-protocol'
report = 'openduck report'
scriptName = 'rdkit_addHydrogens.py'
buildSystem = 'openmmBuildSystem.py'

class ProtOpenDuckSimulation(EMProtocol):
    """
    This protocol will start a undocking simulation using OpenMM
    """
    _label = 'openduck undocking simulation'


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution: ",
                       help="This protocol has both CPU and GPU implementation.\
                                                 Select the one you want to use.")
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

        form.addSection(label=Message.LABEL_INPUT)
        iGroup = form.addGroup('Input')
        iGroup.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input set of molecules:', allowsNull=False,
                        help='Input set of docked molecules. One of them will be prepared together with its target')
        iGroup.addParam('inputLigand', params.StringParam,
                        label='Ligand to prepare: ',
                        help='Specific ligand to prepare in the system')

        iGroup.addParam('doHydrogens', params.BooleanParam, default=True,
                        label='ReAssign hydrogens: ',
                        help='Hydrogens are reomved (if present) and readded')
        iGroup.addParam('doGasteiger', params.BooleanParam, default=True,
                        label='Recalculate gasteiger charges: ',
                        help='Charges of the molecule are recomputed using the gasteiger method')

        intGroup = form.addGroup('Interaction')
        intGroup.addParam('intChain', params.StringParam, allowsNull=False,
                          label='Chain of interaction:',
                          help='Specify the chain of the structure interacting with the ligand with a hydrogen bond')
        intGroup.addParam('intResidue', params.StringParam, allowsNull=False,
                          label='Residue of interaction:',
                          help='Specify the residue of the structure interacting with the ligand with a hydrogen bond')
        intGroup.addParam('intAtom', params.StringParam, allowsNull=False,
                          label='Atom of interaction:',
                          help='Specify the atom of the structure interacting with the ligand with a hydrogen bond')

        cGroup = form.addGroup('Chunking')
        cGroup.addParam('doChunk', params.BooleanParam, default=True, label='Perform chunking: ',
                        help='Chunk the protein structure to only take into account the surroundings of the '
                             'defined interaction in the MD simulation')
        cGroup.addParam('cutoff', params.FloatParam, default=9, label="Chunking cutoff (A): ", condition='doChunk',
                        help='Cutoff distance to define chunking (in Angstroms)')

        form.addSection(label='Parametrization')
        pGroup = form.addGroup('Parametrization')
        pGroup.addParam('smallFF', params.EnumParam, default=0, label="Small molecule forcefield: ",
                        choices=['SMIRNOFF', 'GAFF', 'ESPALOMA'],
                        help='Small molecule forcefield to use for parameterization.')
        pGroup.addParam('waterFF', params.EnumParam, default=0, label="Water forcefield: ",
                        choices=['tip3p', 'spce'],
                        help='Water model to parameterize the solvent.')
        pGroup.addParam('proteinFF', params.EnumParam, default=0, label="Protein forcefield: ",
                        choices=['amber99sb', 'amber14-all'],
                        help='Protein forcefield to parameterize the chunked protein.')
        pGroup.addParam('ion', params.FloatParam, default=0.1, label="Ionic strength (M): ",
                        help='Ionic strength (concentration) of the counter ion salts (Na+/Cl-)')
        pGroup.addParam('buffer', params.FloatParam, default=10, label="Solvent buffer (A): ",
                        help='Buffer distance between the periodic box and the protein (in Angstroms)')
        pGroup.addParam('doHMR', params.BooleanParam, default=False, label='Perform HMR: ',
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Perform Hydrogen Mass Repartition on the topology and run at dt=0.04 ps.')

        sGroup = form.addGroup('Simulation')
        sGroup.addParam('fConst', params.FloatParam, default=1.0, label="Force constant for equilibration: ",
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Force constant for equilibration')
        sGroup.addParam('nCycles', params.IntParam, default=20, label="Number of MD/SMD cycles: ",
                        help='Number of MD/SMD cycles to perform')
        sGroup.addParam('mdLength', params.FloatParam, default=0.5, label="MD length (ns): ",
                        help='Length of MD sampling between SMD runs in ns')
        sGroup.addParam('iniVel', params.FloatParam, default=0.00001, label="Initial velocities: ",
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Set initial velocities when heating.')
        sGroup.addParam('iniDist', params.FloatParam, default=2.5, label="Initial hydrogen bond distance (A): ",
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Set initial hydrogen bond distance for SMD in Angstroms')

    def _insertAllSteps(self):
      self._insertFunctionStep('prepareStep')
      self._insertFunctionStep('simulateStep')
      self._insertFunctionStep('createOutputStep')


    def prepareStep(self):
      # pdbfixer in receptor, addHs in ligand
      args = f'{self.getInputReceptorFile()} --output {self.getPreparedReceptorFile()}'
      pwchemPlugin.runOPENBABEL(self, 'pdbfixer', args=args, cwd=self._getExtraPath())

      self.addLigandHydrogens()

    def simulateStep(self):
      os.mkdir(self.getOutputDir())
      paramsFile = self.writeSimParamsFile()
      Plugin.runOpenMM(self, program, args=f'-y {paramsFile}', cwd=self.getOutputDir())

    def createOutputStep(self):
      args = f'-p {self.getOutputDir()} -f openmm --plot -of csv -o openDuckW_min.csv'
      Plugin.runOpenMM(self, report, args=args, cwd=self._getPath())

      args = f'-p {self.getOutputDir()} -f openmm --plot -d jarzynski -of csv -o openDuckW_jarzynski.csv'
      Plugin.runOpenMM(self, report, args=args, cwd=self._getPath())

      outPDB, outDCD = self.copyOutputSimulations()
      args = self.writeBuildParamsFile(os.path.abspath(outPDB))
      pwchemPlugin.runScript(self, buildSystem, args, env=OPENMM_DIC, cwd=self._getPath())
      Plugin.runOpenMM(self, buildSystem, args=args, cwd=self._getPath())

      systemFile = outDCD.replace('.dcd', '_system.xml')
      mFF, wFF = self.getFFFiles()

      nFrames = self.nSteps.get() // self.nTraj.get()
      nTime = nFrames * self.stepSize.get()
      outSystem = OpenMMSystem(filename=outPDB, serieFile=systemFile,
                               ff=mFF, wff=wFF)
      outSystem.setOriStructFile(outPDB)
      outSystem.setTrajectoryFile(outDCD)

      self._defineOutputs(outputSystem=outSystem)


    def _warnings(self):
      ws = []
      return ws

    def _summary(self):
      summ = []
      if os.path.exists(self.getOutWorkFile(False)):
        wqb = self.parseOutputCSV(False)
        if wqb is not None:
          summ += [f'Min Wqb: {wqb}']
      if os.path.exists(self.getOutWorkFile(True)):
        wqb = self.parseOutputCSV(True)
        if wqb is not None:
          summ += [f'Jarzynski Wqb: {wqb}\n']
      return summ

    ##################### UTILS FUNCTIONS ##################################

    def getOutWorkFile(self, jar=True):
      jarStr = 'jarzynski' if jar else 'min'
      return self._getPath(f'openDuckW_{jarStr}.csv')

    def parseOutputCSV(self, jar=True):
      with open(self.getOutWorkFile(jar)) as f:
        f.readline()
        score = float(f.readline().strip().split(',')[1])
      return score

    def getInputReceptorFile(self):
      return os.path.abspath(self.inputSetOfMols.get().getProteinFile())

    def getPreparedReceptorFile(self):
      recFile = self.getInputReceptorFile()
      return os.path.abspath(self._getExtraPath(recFile.split('/')[-1]))

    def getInteractionStr(self):
      chainID, atomID = json.loads(self.intChain.get())['chain'], json.loads(self.intAtom.get())['atom']
      resID = json.loads(self.intResidue.get())['index'].split('-')[0]
      resType = RESIDUES1TO3[json.loads(self.intResidue.get())['residues']]
      return f'{chainID}_{resType}_{resID}_{atomID}'

    def getSpecifiedMolFile(self):
        myMol = None
        for mol in self.inputSetOfMols.get():
          if mol.__str__() == self.inputLigand.get():
            myMol = mol.clone()
            break
        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            molFile = myMol.getPoseFile()
            return convertToSdf(self, molFile)

    def addLigandHydrogens(self):
      os.mkdir(self.getLigandFileDir())
      paramFile = self.writePrepParamsFile()
      pwchemPlugin.runScript(self, scriptName, paramFile, env=RDKIT_DIC, cwd=self._getPath())
      os.rename(self.getPreparedLigandFile(), self.getPreparedLigandFile().replace('.sdf', '.mol'))

    def getLigandFileDir(self):
      return os.path.abspath(self._getExtraPath('ligand'))

    def getPreparedLigandFile(self):
      ligFile = os.listdir(self.getLigandFileDir())[0]
      return os.path.join(self.getLigandFileDir(), ligFile)

    def writePrepParamsFile(self):
        molFiles = [self.getSpecifiedMolFile()]
        paramsFile = self.getLigParamFile()
        with open(paramsFile, 'w') as f:
            molFs = [os.path.abspath(molFile) for molFile in molFiles]
            f.write(f"ligandFiles: {' '.join(molFs)}\n")

            f.write(f'outputDir: {self.getLigandFileDir()}\n')
            f.write(f'doHydrogens: {self.doHydrogens.get()}\n')
            f.write(f'doGasteiger: {self.doGasteiger.get()}\n')
        return paramsFile

    def writeSimParamsFile(self):
      intStr = self.getInteractionStr()
      with open(self.getSimParamsFile(), 'w') as f:
        f.write('# Main Arguments\n')
        f.write(f'interaction : {intStr}\n')
        f.write(f'receptor_pdb : {self.getPreparedReceptorFile()}\n')
        f.write(f'ligand_mol : {self.getPreparedLigandFile()}\n')
        if getattr(self, params.USE_GPU).get():
          f.write(f'gpu_id : {getattr(self, params.GPU_LIST).get()}\n')

        f.write('\n# Chunking Arguments\n')
        f.write(f'do_chunk : {self.doChunk.get()}\n')
        if self.doChunk.get():
          f.write(f'cutoff : {self.cutoff.get()}\n')

        f.write('\n# Preparation Arguments\n')
        f.write(f'small_molecule_forcefield : {self.getEnumText("smallFF").lower()}\n')
        f.write(f'protein_forcefield : {self.getEnumText("proteinFF").lower()}\n')
        f.write(f'water_model : {self.getEnumText("waterFF").lower()}\n')

        f.write(f'ionic_strength : {self.ion.get()}\n')
        f.write(f'solvent_buffer_distance : {self.buffer.get()}\n')
        f.write(f'force_constant_eq : {self.fConst.get()}\n')

        f.write('\n# Production Arguments\n')
        f.write(f'smd_cycles : {self.nCycles.get()}\n')
        f.write(f'md_length : {self.mdLength.get()}\n')
        f.write(f'init_velocities : {self.iniVel.get()}\n')
        f.write(f'init_distance : {self.iniDist.get()}\n')

      return self.getSimParamsFile()

    def writeBuildParamsFile(self, outPDB):
        paramsFile = self.getLigParamFile()
        with open(paramsFile, 'w') as f:
            # receptofile is actually the complex file
            f.write(f"receptorFile :: {outPDB}\n")
            f.write(f"ligandFile :: {os.path.abspath(self.getSpecifiedMolFile())}\n")

            mFF, wFF = self.getFFFiles()
            f.write(f'mFF :: {mFF}\n')
            f.write(f'wFF :: {wFF}\n')
            f.write(f'ligandFF :: {self.getEnumText("smallFF").lower()}\n')

            # Parameters used by OpenDuck
            f.write(f'nonbondedMethod :: PME\n')
            f.write(f'nonbondedCutoff :: 0.9\n')
            f.write(f'constraints :: HBonds\n')

        return paramsFile

    def getLigParamFile(self):
      return os.path.abspath(self._getExtraPath('addHydrogens.txt'))

    def getSimParamsFile(self):
      return os.path.abspath(self._getExtraPath('simulationParams.yaml'))

    def getOutputDir(self, path=''):
      return os.path.abspath(os.path.join(self._getExtraPath('simulation'), path))

    def getFFFiles(self):
      mFF = '{}.xml'.format(self.getEnumText('proteinFF'))
      wFF = '{}.xml'.format(self.getEnumText('waterFF'))

      return mFF, wFF

    def copyOutputSimulations(self):
      outPDB = self.getOutputDir('duck_runs/smd_0_300.pdb')
      nOutPDB = self._getPath(f'{self.inputLigand.get()}_openduck.pdb')
      shutil.copy(outPDB, nOutPDB)

      outDCD = self.getOutputDir('duck_runs/smd_0_300.dcd')
      nOutDCD = self._getPath(f'{self.inputLigand.get()}_openduck.dcd')
      shutil.copy(outDCD, nOutDCD)

      return nOutPDB, nOutDCD



