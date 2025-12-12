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
This module will prepare the system for the simulation
"""
import os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem import Plugin as pwchemPlugin

from .. import Plugin
from ..constants import OPENMM_DIC
from ..objects import OpenMMSystem
from pwem.convert import cifToPdb


class ProtOpenMMSystemSimulation(EMProtocol):
    """
    This protocol will start a Molecular Dynamics simulation.
    """
    _label = 'system simulation'
    stepsExecutionMode = params.STEPS_PARALLEL


    # -------------------------- DEFINE param functions ----------------------
    def _defineMinimization(self, form):
        form.addParam('addMinimization', params.BooleanParam, default=True, label="Add minimization: ",
                      help='Add energy minimization')
        form.addParam('minimTol', params.FloatParam, default=10, label="Minimization tolerance (kJ/mol): ",
                      condition='addMinimization',
                      help='This specifies how precisely the energy minimum must be located.  Minimization is halted '
                           'once the root-mean-square value of all force components reaches this tolerance.')
        form.addParam('maxIter', params.IntParam, default=10000, label="Maximum iterations: ",
                      condition='addMinimization',
                      help='The maximum number of iterations to perform.  If this is 0, minimization is continued until'
                           ' the results converge without regard to how many iterations it takes.')
        return form

    def _defineIntegrator(self, form):
        form.addParam('integrator', params.EnumParam, default=1, label="Simulation integrator: ",
                      choices=['Verlet', 'Langevin', 'LangevinMiddle', 'NoseHoover', 'Brownian', 'VariableVerlet',
                               'VariableLangevin'],
                      help='http://docs.openmm.org/latest/userguide/theory/04_integrators.html')

        form.addParam('stepSize', params.FloatParam, default=0.004, label="Step size for integration (ps): ",
                      condition='not integrator in [5, 6]',
                      help='The step size with which to integrate the system (in picoseconds)')
        form.addParam('fricCoef', params.FloatParam, default=1, label="Friction coefficient (1/ps): ",
                      condition='integrator in [1, 2, 4, 6]',
                      help='The friction coefficient which couples the system to the heat bath (in inverse picoseconds)')
        form.addParam('temperature', params.FloatParam, default=300, label="Simulation temperature (K): ",
                      condition='integrator in [1, 2, 3, 4, 6]', help='Temperature for the simulation')
        form.addParam('colFreq', params.FloatParam, default=1, label="Collision frequency (1/ps): ",
                      condition='integrator in [3]',
                      help='The friction coefficient which couples the system to the heat bath (in inverse picoseconds)')

        form.addParam('errTol', params.FloatParam, default=0.001, label="Error tolerance: ",
                      condition='integrator in [5, 6]', help='The error tolerance')
        return form

    def _defineBarostat(self, form):
        form.addParam('addBarostat', params.BooleanParam, default=False, label="Add barostat: ",
                      help='Add MonteCarlo Barostat to run a NPT simulation')
        form.addParam('pressure', params.FloatParam, default=1, label="Pressure (bar): ", condition='addBarostat',
                      help='The default pressure acting on the system (in bar)')
        form.addParam('barFreq', params.IntParam, default=25, label="Barostat frequency: ",
                      condition='addBarostat',
                      help='The frequency at which Monte Carlo pressure changes should be attempted (in time steps)')
        return form

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
        form.addParam('inputSystem', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='OpenMMSystem', help='OpenMMSystem to execute the simulation over')
        form.addParam('nSteps', params.IntParam, default=10000, label="Number of simulation steps: ",
                      help='Number of steps for simulation')

        form.addParam('cph', params.BooleanParam, default=False, label="Simulate at non-neutral ph: ",
                        help='Whether to simulate at non-neutral ph. If True, it will use OpenMM Constant pH (https://github.com/openmm/openmm-cph).')

        phGroup = form.addGroup('pH', condition='cph')
        phGroup.addParam('singlePH', params.BooleanParam, default=True, label="Use one single pH value: ",
                         help='Choose whether to use one single pH value or an array of them.')
        phGroup.addParam('onePH', params.FloatParam, default=7.5, label="pH value: ", condition='singlePH',
                         help='The pH value to use.')
        phGroup.addParam('manyPH', params.StringParam, default='6.5, 7.0, 7.5, 8.0, 8.5', label="pH values: ",
                         condition='not singlePH',
                         help='The pH values to use, separated with commas.')
        titrGroup = form.addGroup('Titration', condition='cph')
        titrGroup.addParam('residuesToTitrate', params.StringParam, default='ASP, GLU, CYS, HIS, LYS',
                           label="Residues to titrate",
                           help='Residues that will be considered for constant pH titration')
        ffGroup = form.addGroup('Constant pH parameters', condition='cph')
        ffGroup.addParam('implicitModel', params.EnumParam, default=1,
                         choices=['OBC1', 'OBC2', 'GBn', 'GBn2'],
                         label='Implicit solvent model: ',
                         help='Choose the Generalized Born implicit solvent model (only for Amber force fields)')
        ffGroup.addParam('constraintsImp', params.EnumParam, default=1, label="Implicit force field constraints: ",
                         choices=['None', 'HBonds', 'AllBonds', 'HAngles'], condition='cph',
                         help='You can optionally tell OpenMM to constrain certain bond lengths and angles.'
                              'https://docs.openmm.org/latest/userguide/application/02_running_sims.html#constraints')
        ffGroup.addParam('nonbondedMethodImp', params.EnumParam, default=1, condition='cph',
                         choices=['NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'Ewald', 'PME', 'LJPME'],
                         label="Implicit non bonded method: ",
                         help='Non bonded method to simulate the non bonded atom interactions')
        ffGroup.addParam('implicitCutoff', params.FloatParam, default=2.0, label="Implicit cutoff (nm)",
                         expertLevel=params.LEVEL_ADVANCED, condition='cph',
                         help='Cutoff distance for nonbonded interactions in implicit solvent')
        simGroup = form.addGroup('Constant pH simulation', condition='cph', expertLevel=params.LEVEL_ADVANCED,)
        simGroup.addParam('relaxSteps', params.IntParam, default=500, label="Relaxation steps",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of steps for initial relaxation of the system')
        simGroup.addParam('equilSteps', params.IntParam, default=100, label="Equilibration steps",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of equilibration cycles')
        simGroup.addParam('stepEquil', params.IntParam, default=1, label="Steps per equilibration cycle",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of MD steps per equilibration cycle')
        simGroup.addParam('stepProd', params.IntParam, default=1, label="Steps per production cycle",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of MD steps per production cycle')


        tGroup = form.addGroup('Trajectory')
        tGroup.addParam('nTraj', params.IntParam, default=100, label="Steps interval: ",
                        help='Save the state of the system each x steps for the trajectory')
        tGroup.addParam('useOpenmmdl', params.BooleanParam, default=True, label="Analyze trajectory with OpenMMDL: ",
                        help='Whether to analyze the trajectory with OpenMMDL')

        mGroup = form.addGroup('Minimization')
        self._defineMinimization(mGroup)

        iGroup = form.addGroup('Integrator')
        self._defineIntegrator(iGroup)

        bGroup = form.addGroup('Barostat')
        self._defineBarostat(bGroup)

        form.addParallelSection(threads=4, mpi=1)


    def _insertAllSteps(self):
      if (self.cph.get()):
        self._insertFunctionStep(self.createParamsFileStep)
        self._insertFunctionStep(self.cphSimulateStep)
      else:
        self._insertFunctionStep(self.simulateStep)
      if self.useOpenmmdl.get() and self.inputSystem.get().getLigTopologyFile():
        self._insertFunctionStep(self.analyzeStep)
      self._insertFunctionStep(self.createOutputStep)

    def createParamsFileStep(self):
        """Write simulation parameters to TXT file."""
        paramsFile = self.getParamsFile()

        recFile = self.getStructureFile() #cif file
        #pdbFile = self._getExtraPath(f'{self.getSystemName()}.pdb')
        #cifToPdb(os.path.abspath(recFile), (pdbFile))
        molFile = self.inputSystem.get().getLigTopologyFile()

        txtFilePath = os.path.join(os.path.dirname(recFile), 'extra')
        txtFile = os.path.join(txtFilePath, 'solvationParams.txt')
        pdbFile = os.path.join(txtFilePath, f'{self.getSystemName().split("_")[0]}.pdb')

        solvParams = self.readSolvParams(txtFile)
        with open(paramsFile, 'w') as f:
            # Imports in script
            home = Plugin.getVar(OPENMM_DIC['home'])
            scriptsDir = os.path.join(home, 'openmm-cph')

            f.write(f"constantPHScript = {os.path.abspath(os.path.join(scriptsDir, 'constantph.py'))}\n")
            f.write(f"referenceEnergyScript = {os.path.abspath(os.path.join(scriptsDir, 'reference_energy.py'))}\n")

            # Input system
            f.write(f"inputPdb = {os.path.abspath(pdbFile)}\n")

            # Force fields
            mff, wff = self.getFFFiles()
            f.write(f"explicitFF = {mff}\n")
            f.write(f"explicitSolvent = {wff}\n")
            f.write(f"constraintsExp = {solvParams.get('constraints')}\n")
            mffImp, wffImp = self.getImplicitFF()
            f.write(f"implicitFF = {mffImp}\n")
            f.write(f"implicitSolvent = {wffImp}\n")
            f.write(f"constraintsImp = {self.getEnumText('constraintsImp')}\n")
            if molFile:
                f.write(f"ligandFile = {os.path.abspath(molFile)}\n")
                f.write(f"ligandFF = {solvParams.get('ligandFF')}\n")

            # Cutoffs and hydrogen mass
            f.write(f"nonBondedMethodExp = {solvParams.get('nonbondedMethod')}\n")
            f.write(f"nonBondedMethodImp = {self.getEnumText('nonbondedMethodImp')}\n")
            f.write(f"explicitCutoff = {solvParams.get('nonbondedCutoff')}\n")
            f.write(f"implicitCutoff = {self.implicitCutoff.get()}\n")

            # Simulation steps
            f.write(f"nSteps = {self.nSteps.get()}\n")
            f.write(f"relaxSteps = {self.relaxSteps.get()}\n")
            f.write(f"equilSteps = {self.equilSteps.get()}\n")
            f.write(f"stepEquil = {self.stepEquil.get()}\n")
            f.write(f"stepProd = {self.stepProd.get()}\n")
            f.write(f"reportEvery = {self.nTraj.get()}\n")

            # Residues to titrate
            f.write(f"residuesToTitrate = {self.residuesToTitrate.get()}\n")

            home = Plugin.getVar(OPENMM_DIC['home'])
            aspPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/ASP.pdb')
            gluPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/GLU.pdb')
            cysPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/CYS.pdb')
            hisPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/HIS.pdb')
            lysPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/LYS.pdb')
            # Reference models
            f.write(f"aspModel = {aspPDB}\n")
            f.write(f"gluModel = {gluPDB}\n")
            f.write(f"hisModel = {hisPDB}\n")
            f.write(f"cysModel = {cysPDB}\n")
            f.write(f"lysModel = {lysPDB}\n")

            # pH values
            f.write(f"singlePH = {self.singlePH.get()}\n")
            f.write(f"onePH = {self.onePH.get()}\n")
            f.write(f"manyPH = {self.manyPH.get()}\n")

            # Minimization
            f.write(f"addMinimization = {str(self.addMinimization.get())}\n")
            f.write(f"minimTol = {self.minimTol.get()}\n")
            f.write(f"maxIter = {self.maxIter.get()}\n")

            # Add hydrogens
            f.write(f"addHydrogens = {solvParams.get('addH')}\n")
            f.write(f"hPH = {solvParams.get('hPH')}\n")

            # Barostat
            f.write(f"addBarostat = {str(self.addBarostat.get())}\n")
            if self.addBarostat.get():
                f.write(f"pressure = {self.pressure.get()}\n")

            # Integrator
            f.write(f"integrator = {self.getEnumText('integrator')}\n")
            f.write(f"stepSize = {self.stepSize.get()}\n")
            f.write(f"fricCoef = {self.fricCoef.get()}\n")
            f.write(f"temperature = {self.temperature.get()}\n")
            f.write(f"colFreq = {self.colFreq.get()}\n")
            f.write(f"errTol = {self.errTol.get()}\n")

            #Output paths
            sysName = self.getSystemName()
            trajFile = self._getPath(f"{sysName}.dcd")
            f.write(f"trajFile = {os.path.abspath(trajFile)}\n")
            logFile = self._getPath("md_log.txt")
            f.write(f"logFile = {os.path.abspath(logFile)}\n")
            finalPdb = self._getPath(f"{sysName}.pdb")
            f.write(f"finalPdb = {os.path.abspath(finalPdb)}\n")
            finalCif = self._getPath(f"{sysName}.cif")
            f.write(f"finalCif = {os.path.abspath(finalCif)}\n")
            f.write(f'systemXml = {self.getSystemFile()}\n')

            # Solvation box etc
            if solvParams.get('boxSize'):
                f.write(f"boxSize = {solvParams.get('boxSize')}\n")
            else:
                f.write(f"padding = {solvParams.get('padDist')}\n")

            f.write(f"saltConc = {solvParams.get('saltConc')}\n")
            f.write(f"neutralize = {solvParams.get('neutralize')}\n")
            f.write(f"cationType = {solvParams.get('cationType')}\n")
            f.write(f"anionType = {solvParams.get('anionType')}\n")

        print(f"Parameters file created at: {paramsFile}")

    def cphSimulateStep(self):
        paramsFile = self.getParamsFile()
        Plugin.runScript(self, 'openmmConstantpH.py', args=f'--params {paramsFile}', env=OPENMM_DIC,
                         cwd=self._getPath())

    def simulateStep(self):
      sysFile, structFile = self.getSystemFile(), self.getStructureFile()

      with open(self.getParamsFile(), 'w') as f:
        f.write(f'systemFile :: {sysFile}\n')
        f.write(f'structureFile :: {structFile}\n')
        f.write(f'nSteps :: {self.nSteps.get()}\n')

        integrator = self.getEnumText('integrator')
        f.write(f'integrator :: {integrator}\n')
        if self.integrator.get() not in [0, 5]:
          f.write(f'temperature :: {self.temperature.get()}\n')

        if self.integrator.get() not in [5, 6]:
          f.write(f'stepSize :: {self.stepSize.get()}\n')

        if self.integrator.get() not in [0, 3, 5]:
          f.write(f'fricCoef :: {self.fricCoef.get()}\n')

        f.write(f'addMinimization :: {self.addMinimization.get()}\n')
        if self.addMinimization:
          f.write(f'minimTol :: {self.minimTol.get()}\n')
          f.write(f'maxIter :: {self.maxIter.get()}\n')

        f.write(f'addBarostat :: {self.addBarostat.get()}\n')
        if self.addBarostat:
          f.write(f'pressure :: {self.pressure.get()}\n')
          f.write(f'temperature :: {self.temperature.get()}\n')

        f.write(f'nTraj :: {self.nTraj.get()}\n')
        if getattr(self, params.USE_GPU).get():
          f.write(f'gpus :: {getattr(self, params.GPU_LIST)}\n')

      Plugin.runScript(self, 'openmmSimulateSystem.py', args=self.getParamsFile(), env=OPENMM_DIC,
                             cwd=self._getPath())

    def getNFrames(self):
      nFrames = self.nSteps.get() // self.nTraj.get()
      return nFrames

    def analyzeStep(self):
        '''Run OpenMMDL analysis'''
        oDir = self._getExtraPath('OpenMMDL')
        if not os.path.exists(oDir):
          os.mkdir(oDir)
        systemName = self.getSystemName()
        outTopFile, outDcdFile = os.path.abspath(self._getPath(f'{systemName}.pdb')), \
                                 os.path.abspath(self._getPath(f'{systemName}.dcd'))

        args = f'-t {outTopFile} -d {outDcdFile} -n LIG -c {self.numberOfThreads.get()}'
        pwchemPlugin.runCondaCommand(self, args, OPENMM_DIC, 'openmmdl_analysis', cwd=oDir)


    def createOutputStep(self):
      systemName = self.getSystemName()
      systemFile = os.path.relpath(self.getSystemFile())
      outTopFile, outDcdFile = self._getPath(f'{systemName}.pdb'), self._getPath(f'{systemName}_wrapped.dcd')
      outCifFile = self._getPath(f'{systemName}.cif')

      mFF, wFF = self.getFFFiles()
      nFrames = self.getNFrames()
      nTime = nFrames * self.stepSize.get()
      outSystem = OpenMMSystem(filename=outTopFile, serieFile=systemFile, cifFile=outCifFile,
                               repFile=self._getPath('md_log.txt'),
                               ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)
      outSystem.setTrajectoryFile(outDcdFile)
      print(outDcdFile)

      ligFile = self.inputSystem.get().getLigTopologyFile()
      if ligFile:
        outSystem.setLigTopologyFile(ligFile)

      anaDir = self._getExtraPath('OpenMMDL')
      if os.path.exists(anaDir):
        outSystem.setOpenmmdlDir(anaDir)

      self._defineOutputs(outputSystem=outSystem)


    def _warnings(self):
      ws = []
      if not self.addMinimization.get():
        ws.append('Running the simulation without a prior minimization might lead to errors in the simulation.\n')
      return ws


    def getWaterModel(self, wFF):
      model = 'tip3p'
      if 'spce' in wFF:
        model = 'spce'
      elif 'tip4p' in wFF:
        model = 'tip4pew'
      elif 'tip5p' in wFF:
        model = 'tip5p'
      return model

    def getFFFiles(self):
      system = self.inputSystem.get()
      return system.getForceField(), system.getWaterForceField()

    def getNBParams(self):
      system = self.inputSystem.get()
      return system._nbMethod.get(), system._nbCutoff.get()

    def getParamsFile(self):
      return os.path.abspath(self._getExtraPath('simulationParams.txt'))

    def getStructureFile(self):
      return os.path.abspath(self.inputSystem.get().getCifFile())

    def getSystemFile(self):
      return os.path.abspath(self.inputSystem.get().getSerieFile())

    def getSystemName(self):
      return self.inputSystem.get().getSystemName()

    def getImplicitFF(self):
        mFF, _ = self.getFFFiles()
        model = self.implicitModel.get()

        implicitDict = {
            'OBC1': 'implicit/obc1.xml',
            'OBC2': 'implicit/obc2.xml',
            'GBn': 'implicit/gbn.xml',
            'GBn2': 'implicit/gbn2.xml'
        }

        wFF = implicitDict.get(model, 'implicit/gbn.xml')
        return mFF, wFF

    def readSolvParams(self, txtFile):
        """
        Reads a solvationParams.txt file and returns a dictionary
        of parameter names and their values.
        """
        if not os.path.exists(txtFile):
            raise FileNotFoundError(f"File not found: {txtFile}")

        paramsDict = {}
        with open(txtFile, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if '::' in line:
                    key, value = line.split('::', 1)
                    key = key.strip()
                    value = value.strip()
                    if value.lower() in ['true', 'false']:
                        value = value.lower() == 'true'
                    else:
                        try:
                            if '.' in value:
                                value = float(value)
                            else:
                                value = int(value)
                        except ValueError:
                            pass
                    paramsDict[key] = value

        return paramsDict
