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

from .. import Plugin
from ..constants import OPENMM_DIC
from ..objects import OpenMMSystem


class ProtOpenMMSystemSimulation(EMProtocol):
    """
    This protocol will start a Molecular Dynamics simulation.
    """
    _label = 'system simulation'


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
        form.addParam('inputSystem', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='OpenMMSystem', help='OpenMMSystem to execute the simulation over')
        form.addParam('nSteps', params.IntParam, default=10000, label="Number of simulation steps: ",
                      help='Number of steps for simulation')

        tGroup = form.addGroup('Trajectory')
        tGroup.addParam('nTraj', params.IntParam, default=100, label="Steps interval: ",
                        help='Save the state of the system each x steps for the trajectory')

        mGroup = form.addGroup('Minimization')
        mGroup.addParam('addMinimization', params.BooleanParam, default=True, label="Add minimization: ",
                      help='Add energy minimization')
        mGroup.addParam('minimTol', params.FloatParam, default=10, label="Minimization tolerance (kJ/mol): ",
                      condition='addMinimization',
                      help='This specifies how precisely the energy minimum must be located.  Minimization is halted '
                           'once the root-mean-square value of all force components reaches this tolerance.')
        mGroup.addParam('maxIter', params.IntParam, default=10000, label="Maximum iterations: ",
                      condition='addMinimization',
                      help='The maximum number of iterations to perform.  If this is 0, minimization is continued until'
                           ' the results converge without regard to how many iterations it takes.')

        iGroup = form.addGroup('Integrator')
        iGroup.addParam('integrator', params.EnumParam, default=1, label="Simulation integrator: ",
                      choices=['Verlet', 'Langevin', 'LangevinMiddle', 'NoseHoover', 'Brownian', 'VariableVerlet',
                               'VariableLangevin'],
                      help='http://docs.openmm.org/latest/userguide/theory/04_integrators.html')

        iGroup.addParam('stepSize', params.FloatParam, default=0.004, label="Step size for integration (ps): ",
                      condition='not integrator in [5, 6]',
                      help='The step size with which to integrate the system (in picoseconds)')
        iGroup.addParam('fricCoef', params.FloatParam, default=1, label="Friction coefficient (1/ps): ",
                      condition='integrator in [1, 2, 4, 6]',
                      help='The friction coefficient which couples the system to the heat bath (in inverse picoseconds)')
        iGroup.addParam('temperature', params.FloatParam, default=300, label="Simulation temperature (K): ",
                      condition='integrator in [1, 2, 3, 4, 6]',
                      help='Temperature for the simulation')
        iGroup.addParam('colFreq', params.FloatParam, default=1, label="Collision frequency (1/ps): ",
                      condition='integrator in [3]',
                      help='The friction coefficient which couples the system to the heat bath (in inverse picoseconds)')

        iGroup.addParam('errTol', params.FloatParam, default=0.001, label="Error tolerance: ",
                      condition='integrator in [5, 6]',
                      help='The error tolerance')

        bGroup = form.addGroup('Barostat')
        bGroup.addParam('addBarostat', params.BooleanParam, default=False, label="Add barostat: ",
                      help='Add MonteCarlo Barostat to run a NPT simulation')
        bGroup.addParam('pressure', params.FloatParam, default=1, label="Pressure (bar): ", condition='addBarostat',
                      help='The default pressure acting on the system (in bar)')
        bGroup.addParam('barFreq', params.IntParam, default=25, label="Barostat frequency: ",
                      condition='addBarostat',
                      help='The frequency at which Monte Carlo pressure changes should be attempted (in time steps)')

    def _insertAllSteps(self):
      self._insertFunctionStep('simulateStep')
      self._insertFunctionStep('createOutputStep')


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


    def createOutputStep(self):
      systemName = self.getSystemName()
      oriStructFile, systemFile = self.getStructureFile(), self.getSystemFile()
      outPdbFile, outDcdFile = self._getPath(f'{systemName}.pdb'), self._getPath(f'{systemName}.dcd')

      mFF, wFF = self.getFFFiles()
      nFrames = self.nSteps.get() // self.nTraj.get()
      nTime = nFrames * self.stepSize.get()
      outSystem = OpenMMSystem(filename=outPdbFile, serieFile=systemFile,
                               repFile=self._getPath('md_log.txt'),
                               ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)
      outSystem.setOriStructFile(oriStructFile)
      outSystem.setTrajectoryFile(outDcdFile)

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
      return os.path.abspath(self.inputSystem.get().getFileName())

    def getSystemFile(self):
      return os.path.abspath(self.inputSystem.get().getSerieFile())

    def getSystemName(self):
      return self.inputSystem.get().getSystemName()
