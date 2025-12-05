# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
This module will run the simulation for a constant pH
"""
import os, sys

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem import Plugin as pwchemPlugin

from .. import Plugin
from pwchem.utils import getBaseName, convertToSdf
from ..constants import OPENMM_DIC
from ..objects import OpenMMSystem

from pwem.convert import cifToPdb


class ProtOpenMMSystemSimulationConstantPH(EMProtocol):
    """
    This protocol will start a Molecular Dynamics simulation with constant pH specified by user.
    """
    _label = 'constant pH system simulation'

    stepsExecutionMode = params.STEPS_PARALLEL


    IMPLICIT_SOLVENT_MAP = {
        'obc1': 'implicit/amber99_obc.xml',
        'obc2': 'implicit/amber99_obc2.xml',
        'gbn': 'implicit/gbn.xml',
        'gbn2': 'implicit/gbn2.xml',
        'hct': 'implicit/hct.xml'
    }

    # -------------------------- DEFINE param functions ----------------------
    def _defineFFParams(self, form, ligandCondition='True'):
        form.addParam('ffType', params.EnumParam, default=0, choices=['Amber14', 'CHARMM36', 'Old'],
                      label="Main atomic force field: ", help='Main force field to use')
        form.addParam('ffAmberType', params.EnumParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                      condition='ffType==0', label="Amber atomic force field: ",
                      choices=['All', 'protein.ff14SB', 'protein.ff15ipq', 'DNA.OL15', 'DNA.bsc1', 'RNA.OL3', 'lipid17'],
                      help='Amber main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')
        form.addParam('ffAmberWaterType', params.EnumParam, default=3, condition='ffType==0',
                      label="Amber water force field: ",
                      choices=['SPCE', 'OPC', 'OPC3', 'tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb'],
                      help='Water amber force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')

        form.addParam('ffCHARMMWaterType', params.EnumParam, default=0, condition='ffType==1',
                      label="CHARMM water force field: ", expertLevel=params.LEVEL_ADVANCED,
                      choices=['Water', 'SPCE', 'tip3p-pme-b', 'tip3p-pme-f', 'tip4pew', 'tip4p2005', 'tip5p', 'tip5pew'],
                      help='Water CHARMM force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#charmm36')

        form.addParam('ffOldType', params.EnumParam, default=0,
                      choices=['amber96', 'amber99sb', 'amber99sbildn', 'amber99sbnmr', 'amber03', 'amber10', 'charmm_polar_2013'],
                      condition='ffType==2', label="Older force field: ",
                      help='Select an older main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#older-force-fields')

        form.addParam('ffWaterType', params.EnumParam, default=0,
                      choices=['tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb', 'tip5p', 'spce', 'swm4ndp', 'opc', 'opc3'],
                      condition='ffType==2', label="Water force field: ",
                      help='Select an water force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#water-models')

        form.addParam('constraints', params.EnumParam, default=1, label="Forcefield constraints: ",
                      choices=['None', 'HBonds', 'AllBonds', 'HAngles'],
                      help='You can optionally tell OpenMM to constrain certain bond lengths and angles.'
                           'http://docs.openmm.org/latest/userguide/application/02_running_sims.html#constraints')
        return form

    def _defineFFImpParams(self, form, ligandCondition='True'):
        form.addParam('ffTypeImp', params.EnumParam, default=0, choices=['Amber14', 'CHARMM36', 'Old'],
                      label="Implicit atomic force field: ", help='Implicit force field to use. It will be used when deciding whether to change the protonation states of residues.')
        form.addParam('ffAmberTypeImp', params.EnumParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                      condition='ffTypeImp==0', label="Amber atomic force field: ",
                      choices=['All', 'protein.ff14SB', 'protein.ff15ipq', 'DNA.OL15', 'DNA.bsc1', 'RNA.OL3', 'lipid17'],
                      help='Amber main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')
        form.addParam('ffAmberWaterTypeImp', params.EnumParam, default=3, condition='ffTypeImp==0',
                      label="Amber water force field: ",
                      choices=['SPCE', 'OPC', 'OPC3', 'tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb'],
                      help='Water amber force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')

        form.addParam('ffCHARMMWaterTypeImp', params.EnumParam, default=0, condition='ffTypeImp==1',
                      label="CHARMM water force field: ", expertLevel=params.LEVEL_ADVANCED,
                      choices=['Water', 'SPCE', 'tip3p-pme-b', 'tip3p-pme-f', 'tip4pew', 'tip4p2005', 'tip5p', 'tip5pew'],
                      help='Water CHARMM force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#charmm36')
        form.addParam('ffOldTypeImp', params.EnumParam, default=0,
                      choices=['amber96', 'amber99sb', 'amber99sbildn', 'amber99sbnmr', 'amber03', 'amber10', 'charmm_polar_2013'],
                      condition='ffTypeImp==2', label="Older force field: ",
                      help='Select an older main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#older-force-fields')

        form.addParam('ffWaterTypeImp', params.EnumParam, default=0,
                      choices=['tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb', 'tip5p', 'spce', 'swm4ndp', 'opc', 'opc3'],
                      condition='ffTypeImp==2', label="Implicit water force field: ",
                      help='Select an water force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#water-models')

        form.addParam('constraintsImp', params.EnumParam, default=1, label="Forcefield constraints: ",
                      choices=['None', 'HBonds', 'AllBonds', 'HAngles'],
                      help='You can optionally tell OpenMM to constrain certain bond lengths and angles.'
                           'http://docs.openmm.org/latest/userguide/application/02_running_sims.html#constraints')
        return form

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

        form.addParam('stepSize', params.FloatParam, default=0.002, label="Step size for integration (ps): ",
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

        form.addParam('inputStructure', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='AtomStruct',
                      help='Atomic structure to be prepared for MD by solvation, ions addition etc.')

        phGroup = form.addGroup('pH ')
        phGroup.addParam('singlePH', params.BooleanParam, default=True, label="Use one single pH value: ",
                      help='Choose whether to use one single pH value or an array of them.')
        phGroup.addParam('onePH', params.FloatParam, default=7.5, label="pH value: ", condition='singlePH',
                      help='The pH value to use.')
        phGroup.addParam('manyPH', params.StringParam, default='6.5, 7.0, 7.5, 8.0, 8.5', label="pH values: ", condition='not singlePH',
                      help='The pH values to use, separated with commas.')

        ffGroup = form.addGroup('Main Force Field Parameters')
        self._defineFFParams(ffGroup, False)
        ffGroup = form.addGroup('Implicit Force Field Parameters')
        self._defineFFImpParams(ffGroup, False)

        simGroup = form.addGroup('Simulation Steps')
        simGroup.addParam('saveInterval', params.IntParam, default=100, label="Trajectory save interval: ",
                          help="Number of steps between saving frames in the trajectory")
        simGroup.addParam('relaxSteps', params.IntParam, default=500, label="Relaxation steps",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of steps for initial relaxation of the system')
        simGroup.addParam('equilSteps', params.IntParam, default=100, label="Equilibration steps",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of equilibration cycles')
        simGroup.addParam('stepEquil', params.IntParam, default=1, label="Steps per equilibration cycle",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of MD steps per equilibration cycle')
        simGroup.addParam('prodSteps', params.IntParam, default=10000, label="Production steps",
                          help='Number of production cycles')
        simGroup.addParam('stepProd', params.IntParam, default=1, label="Steps per production cycle",
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Number of MD steps per production cycle')
        simGroup.addParam('useOpenmmdl', params.BooleanParam, default=True, label="Analyze trajectory with OpenMMDL: ",
                        help='Whether to analyze the trajectory with OpenMMDL')

        titrGroup = form.addGroup('Titration')
        titrGroup.addParam('residuesToTitrate', params.StringParam, default='ASP, GLU, CYS, HIS, LYS',
                           label="Residues to titrate",
                           help='Residues that will be considered for constant pH titration')

        fform = form.addGroup('Non bonded interactions')
        fform.addParam('hydrogenMass', params.FloatParam, default=1.5, label="Hydrogen mass (amu)",
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Mass of hydrogens to use in simulations (can accelerate integration)')
        fform.addParam('nonbondedMethodExp', params.EnumParam, default=2,
                       choices=['NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'Ewald', 'PME', 'LJPME'],
                       label="Explicit non bonded method: ",
                       help='Non bonded method to simulate the non bonded atom interactions')
        fform.addParam('explicitCutoff', params.FloatParam, default=0.9, label="Explicit cutoff (nm)",
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Cutoff distance for nonbonded interactions in explicit solvent')
        fform.addParam('nonbondedMethodImp', params.EnumParam, default=1,
                       choices=['NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'Ewald', 'PME', 'LJPME'],
                       label="Implicit non bonded method: ",
                       help='Non bonded method to simulate the non bonded atom interactions')
        fform.addParam('implicitCutoff', params.FloatParam, default=2.0, label="Implicit cutoff (nm)",
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Cutoff distance for nonbonded interactions in implicit solvent')
        formH = form.addGroup('Hydrogens')
        formH.addParam('addH', params.BooleanParam, default=False,
                       label='Add hydrogens to the system: ', help='Add hydrogens to the system')
        formH.addParam('hPH', params.FloatParam, default=7.5,
                       label='PH for hydrogen addition: ', help='The pH based on which to select variants')

        mGroup = form.addGroup('Minimization')
        self._defineMinimization(mGroup)

        iGroup = form.addGroup('Integrator')
        self._defineIntegrator(iGroup)

        bGroup = form.addGroup('Barostat')
        self._defineBarostat(bGroup)


        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
      self._insertFunctionStep(self.createParamsFileStep)
      self._insertFunctionStep(self.productionRunStep)
      #if self.useOpenmmdl.get() and self.inputSystem.get().getLigTopologyFile():
      #    self._insertFunctionStep(self.analyzeStep)
      self._insertFunctionStep(self.createOutputStep)

    def createParamsFileStep(self):
        """Write simulation parameters to TXT file."""
        paramsFile = self.getParamsFile()
        with open(paramsFile, 'w') as f:
            # Imports in script
            home = Plugin.getVar(OPENMM_DIC['home'])
            scripts_dir = os.path.join(home, 'openmm-cph')

            f.write(f"constantPHScript = {os.path.abspath(os.path.join(scripts_dir, 'constantph.py'))}\n")
            f.write(f"referenceEnergyScript = {os.path.abspath(os.path.join(scripts_dir, 'reference_energy.py'))}\n")

            # Input system
            f.write(f"inputPdb = {self.getStructureFile()}\n")

            # Force fields
            mff, wff = self.getFFFiles()
            f.write(f"explicitFF = {mff}\n")
            f.write(f"explicitSolvent = {wff}\n")
            f.write(f"constraintsExp = {self.getEnumText('constraints')}\n")
            mffImp, wffImp = self.getFFFilesImp()
            f.write(f"implicitFF = {mffImp}\n")
            f.write(f"implicitSolvent = {wffImp}\n")
            f.write(f"constraintsImp = {self.getEnumText('constraintsImp')}\n")
            # f.write(f"ligandFile = {os.path.abspath(self.inputSystem.get().getLigTopologyFile())}\n")
            # f.write(f"ligandFF = {os.path.abspath(self._getExtraPath('ligandFF.xml'))}\n")

            # Cutoffs and hydrogen mass
            f.write(f"nonBondedMethodExp = {self.getEnumText('nonbondedMethodExp')}\n")
            f.write(f"nonBondedMethodImp = {self.getEnumText('nonbondedMethodImp')}\n")
            f.write(f"explicitCutoff = {self.explicitCutoff.get()}\n")
            f.write(f"implicitCutoff = {self.implicitCutoff.get()}\n")
            f.write(f"hydrogenMass = {self.hydrogenMass.get()}\n")

            # Simulation steps
            f.write(f"relaxSteps = {self.relaxSteps.get()}\n")
            f.write(f"equilSteps = {self.equilSteps.get()}\n")
            f.write(f"stepEquil = {self.stepEquil.get()}\n")
            f.write(f"prodSteps = {self.prodSteps.get()}\n")
            f.write(f"stepProd = {self.stepProd.get()}\n")
            f.write(f"reportEvery = {self.saveInterval.get()}\n")

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
            f.write(f"addHydrogens = {self.addH.get()}\n")
            f.write(f"hPH = {self.hPH.get()}\n")

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


        print(f"Parameters file created at: {paramsFile}")

    def productionRunStep(self):
        paramsFile = self.getParamsFile()
        Plugin.runScript(self, 'openmmConstantpH.py', args=f'--params {paramsFile}', env=OPENMM_DIC, cwd=self._getPath())


    def createOutputStep(self):
      systemName = self.getSystemName()
      systemFile = os.path.relpath(self.getSystemFile())
      systemFile = self._getPath(f'{systemName}.xml')
      outTopFile, outDcdFile = self._getPath(f'{systemName}.pdb'), self._getPath(f'{systemName}.dcd')
      outCifFile = self._getPath(f'{systemName}.cif')

      mFF, wFF = self.getFFFiles()
      nFrames = self.getNFrames()
      nTime = nFrames * self.stepSize.get()
      outSystem = OpenMMSystem(filename=outTopFile, serieFile=systemFile, cifFile=outCifFile,
                               repFile=self._getPath('md_log.txt'),
                               ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)
      outSystem.setTrajectoryFile(outDcdFile)

      #ligFile = self.inputSystem.get().getLigTopologyFile()
      #if ligFile:
      #  outSystem.setLigTopologyFile(ligFile)

      anaDir = self._getExtraPath('OpenMMDL')
      if os.path.exists(anaDir):
        outSystem.setOpenmmdlDir(anaDir)

      self._defineOutputs(outputSystem=outSystem)

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

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
      summary = []
      return summary

    def _methods(self):
      methods = []
      return methods

    def _warnings(self):
      ws = []
      if not self.addMinimization.get():
        ws.append('Running the simulation without a prior minimization might lead to errors in the simulation.\n')
      return ws

    # --------------------------- UTILS functions -----------------------------------
    def getFFFiles(self):
        if self.ffType.get() == 0:
            mFF = 'amber14-all.xml' if self.ffAmberType.get() == 0 \
                else 'amber14/{}.xml'.format(self.getEnumText('ffAmberType'))
            wFF = 'amber14/{}.xml'.format(self.getEnumText('ffAmberWaterType').lower())

        elif self.ffType.get() == 1:
            mFF = 'charmm36.xml'
            wFF = 'charmm36/{}.xml'.format(self.getEnumText('ffCHARMMWaterType').lower())

        elif self.ffType.get() == 2:
            mFF = '{}.xml'.format(self.getEnumText('ffOldType'))
            wFF = '{}.xml'.format(self.getEnumText('ffWaterType'))

        return mFF, wFF

    def getFFFilesImp(self):
        if self.ffTypeImp.get() == 0:
            mFF = 'amber14-all.xml' if self.ffAmberTypeImp.get() == 0 \
                else 'amber14/{}.xml'.format(self.getEnumText('ffAmberType'))
            wFF = 'amber14/{}.xml'.format(self.getEnumText('ffAmberWaterType').lower())

        elif self.ffTypeImp.get() == 1:
            mFF = 'charmm36.xml'
            wFF = 'charmm36/{}.xml'.format(self.getEnumText('ffCHARMMWaterType').lower())

        elif self.ffTypeImp.get() == 2:
            mFF = '{}.xml'.format(self.getEnumText('ffOldType'))
            wFF = '{}.xml'.format(self.getEnumText('ffWaterType'))

        return mFF, wFF
    def getNBParams(self):
      system = self.inputSystem.get()
      return system._nbMethod.get(), system._nbCutoff.get()

    def getParamsFile(self):
      return os.path.abspath(self._getExtraPath('simulationParams.txt'))

    def getStructureFile(self):
        proteinFile = self.inputStructure.get().getFileName()
        name = os.path.splitext(os.path.basename(proteinFile))[0]
        pdbFile = self._getExtraPath(f'{name}_system.pdb')
        cifToPdb(os.path.abspath(proteinFile), (pdbFile))
        return os.path.abspath(pdbFile)

    def getSystemFile(self): #we will need to change this
        proteinFile = self.inputStructure.get().getFileName()
        name = os.path.splitext(os.path.basename(proteinFile))[0]
        systemName = self._getPath(f'{name}_system.xml')
        return os.path.abspath(systemName)

    def getSystemName(self):
        return getBaseName(self.getStructureFile())

    def getNFrames(self):
        totalSteps = self.prodSteps.get() * self.stepProd.get()
        saveInterval = self.saveInterval.get()
        nFrames = totalSteps // saveInterval
        return nFrames
