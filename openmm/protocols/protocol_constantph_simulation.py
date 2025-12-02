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
        form.addParam('inputSystem', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='OpenMMSystem', help='OpenMMSystem to execute the simulation over')

        phGroup = form.addGroup('pH ')
        phGroup.addParam('singlePH', params.BooleanParam, default=True, label="Use one single pH value: ",
                      help='Choose whether to use one single pH value or an array of them.')
        phGroup.addParam('onePH', params.FloatParam, default=7.5, label="pH value: ", condition='singlePH',
                      help='The pH value to use.')
        phGroup.addParam('manyPH', params.StringParam, default='6.5, 7.0, 7.5, 8.0, 8.5', label="pH values: ", condition='not singlePH',
                      help='The pH values to use, separated with commas.')

        ffGroup = form.addGroup('Force Field Parameters')
        ffGroup.addParam('implicitSolvent', params.EnumParam, default=0,
                         choices=['obc1', 'obc2', 'gbn', 'gbn2', 'hct'],
                         label='Implicit solvent model: ',
                         help='Implicit solvent to use for constant pH simulation.')
        ffGroup.addParam('explicitCutoff', params.FloatParam, default=0.9, label="Explicit cutoff (nm)",
                         expertLevel=params.LEVEL_ADVANCED,
                         help='Cutoff distance for nonbonded interactions in explicit solvent')
        ffGroup.addParam('implicitCutoff', params.FloatParam, default=2.0, label="Implicit cutoff (nm)",
                         expertLevel=params.LEVEL_ADVANCED,
                         help='Cutoff distance for nonbonded interactions in implicit solvent')
        ffGroup.addParam('hydrogenMass', params.FloatParam, default=1.5, label="Hydrogen mass (amu)",
                         expertLevel=params.LEVEL_ADVANCED,
                         help='Mass of hydrogens to use in simulations (can accelerate integration)')
        ffGroup.addParam('constraints', params.EnumParam, default=1,
                         label="Forcefield constraints",
                         choices=['None', 'HBonds', 'AllBonds', 'HAngles'],
                         help='Optional bond/angle constraints for OpenMM. '
                              'http://docs.openmm.org/latest/userguide/application/02_running_sims.html#constraints')

        simGroup = form.addGroup('Simulation Steps')
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

        titrGroup = form.addGroup('Titration')
        titrGroup.addParam('residuesToTitrate', params.StringParam, default='ASP, GLU, CYS, HIS, LYS',
                           label="Residues to titrate",
                           help='Residues that will be considered for constant pH titration')

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

            # Force fields (assuming methods exist to get these)
            f.write(f"explicitFF = {self.inputSystem.get().getForceField()}\n")
            f.write(f"explicitSolvent = {self.inputSystem.get().getWaterForceField()}\n")
            f.write(f"implicitFF = {self.inputSystem.get().getForceField()}\n")
            solventKey = self.getEnumText("implicitSolvent")
            implicitFF_file = self.IMPLICIT_SOLVENT_MAP[solventKey]
            f.write(f"implicitSolvent = {implicitFF_file}\n")

            f.write(f"constraints = {self.constraints.get()}\n")

            # Cutoffs and hydrogen mass
            f.write(f"explicitCutoff = {self.explicitCutoff.get()}\n")
            f.write(f"implicitCutoff = {self.implicitCutoff.get()}\n")
            f.write(f"hydrogenMass = {self.hydrogenMass.get()}\n")

            # Simulation steps
            f.write(f"relaxSteps = {self.relaxSteps.get()}\n")
            f.write(f"equilSteps = {self.equilSteps.get()}\n")
            f.write(f"stepEquil = {self.stepEquil.get()}\n")
            f.write(f"prodSteps = {self.prodSteps.get()}\n")
            f.write(f"stepProd = {self.stepProd.get()}\n")

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
            if (self.singlePH.get() ):
                f.write(f"onePH = {self.onePH.get()}\n")
            else:
                f.write(f"manyPH = {self.manyPH.get()}\n")

            # Minimization
            f.write(f"addMinimization = {str(self.addMinimization.get())}\n")
            f.write(f"minimTol = {self.minimTol.get()}\n")
            f.write(f"maxIter = {self.maxIter.get()}\n")

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

        print(f"Parameters file created at: {paramsFile}")

    def productionRunStep(self):
        paramsFile = self.getParamsFile()
        Plugin.runScript(self, 'openmmConstantpH.py', args=f'--params {paramsFile}', env=OPENMM_DIC, cwd=self._getPath())


    def createOutputStep(self): #todo this when i see how and if it works
      systemName = self.getSystemName()
      systemFile = os.path.relpath(self.getSystemFile())
      outTopFile, outDcdFile = self._getPath(f'{systemName}.pdb'), self._getPath(f'{systemName}.dcd')
      outCifFile = self._getPath(f'{systemName}.cif')

      mFF, wFF = self.getFFFiles()
      nFrames = self.getNFrames()
      nTime = nFrames * self.stepSize.get()
      outSystem = OpenMMSystem(filename=outTopFile, serieFile=systemFile, cifFile=outCifFile,
                               repFile=self._getPath('md_log.txt'),
                               ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)
      outSystem.setTrajectoryFile(outDcdFile)

      ligFile = self.inputSystem.get().getLigTopologyFile()
      if ligFile:
        outSystem.setLigTopologyFile(ligFile)

      anaDir = self._getExtraPath('OpenMMDL')
      if os.path.exists(anaDir):
        outSystem.setOpenmmdlDir(anaDir)

      self._defineOutputs(outputSystem=outSystem)

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
      system = self.inputSystem.get()
      return system.getForceField(), system.getWaterForceField()

    def getNBParams(self):
      system = self.inputSystem.get()
      return system._nbMethod.get(), system._nbCutoff.get()

    def getParamsFile(self):
      return os.path.abspath(self._getExtraPath('simulationParams.txt'))

    def getStructureFile(self):
      name = os.path.splitext(os.path.basename(self.inputSystem.get().getCifFile()))[0]
      pdbFile = self._getExtraPath(f'{name}.pdb')
      cifToPdb(os.path.abspath(self.inputSystem.get().getCifFile()), (pdbFile))
      return os.path.abspath(pdbFile)

    def getSystemFile(self):
      return os.path.abspath(self.inputSystem.get().getSerieFile())

    def getSystemName(self):
      return self.inputSystem.get().getSystemName()
