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
      self._insertFunctionStep(self.simulateStep)
      if self.useOpenmmdl.get() and self.inputSystem.get().getLigTopologyFile():
        self._insertFunctionStep(self.analyzeStep)
      self._insertFunctionStep(self.createOutputStep)


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
        pwchemPlugin.runCondaCommand(self, args, OPENMM_DIC, 'openmmdl analysis', cwd=oDir)


    def createOutputStep(self):
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
