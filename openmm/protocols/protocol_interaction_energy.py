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
This protocol will calculate the Coulomb and LJ energies of interaction between a protein receptor
and a docked ligand. It can be used to get the average interaction energy in a trajectory.
"""
import os
import numpy as np

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from openmm import Plugin
from openmm.constants import OPENMM_DIC

class ProtOpenMMInteractionEnergy(EMProtocol):
    """
    This protocol will calculate the interaction energy of protein and ligand in a system
    """
    _label = 'system interaction energy'


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSystem', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='OpenMMSystem', help='OpenMMSystem to execute the calculation over')

        form.addParam('nTraj', params.IntParam, default=10, label="Trajectory sampling frequency: ",
                      condition='not inputSystem or inputSystem.hasTrajectory()',
                      help='The interaction energy will be calculated over the frames of the trajectory, taking a '
                           'frame for each of this number of steps.')

        mGroup = form.addGroup('Minimization', condition='not inputSystem or not inputSystem.hasTrajectory()')
        mGroup.addParam('addMinimization', params.BooleanParam, default=True, label="Add minimization: ",
                        help='Add energy minimization to the original system if there is no trajectory')
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

    def _insertAllSteps(self):
      self._insertFunctionStep('simulateStep')
      self._insertFunctionStep('createOutputStep')


    def simulateStep(self):
      sysFile, structFile = self.getSystemFile(), self.getStructureFile()
      trajFile = self.getSystemTrajFile()

      with open(self.getParamsFile(), 'w') as f:
        f.write(f'systemFile :: {sysFile}\n')
        f.write(f'structureFile :: {structFile}\n')
        if trajFile:
            f.write(f'trajFile :: {os.path.abspath(trajFile)}\n')
            f.write(f'nTraj :: {self.nTraj.get()}\n')
        else:
            f.write(f'addMin :: {self.addMinimization.get()}\n')
            if self.addMinimization.get():
              f.write(f'minimTol :: {self.minimTol.get()}\n')
              f.write(f'maxIter :: {self.maxIter.get()}\n')

        f.write(f'mFF :: {self.getSystemFF()}\n')

        integrator = self.getEnumText('integrator')
        f.write(f'integrator :: {integrator}\n')
        if self.integrator.get() not in [0, 5]:
            f.write(f'temperature :: {self.temperature.get()}\n')

        if self.integrator.get() not in [5, 6]:
            f.write(f'stepSize :: {self.stepSize.get()}\n')

        if self.integrator.get() not in [0, 3, 5]:
            f.write(f'fricCoef :: {self.fricCoef.get()}\n')

      Plugin.runScript(self, 'openmmInteractionEnergy.py', args=self.getParamsFile(), env=OPENMM_DIC,
                             cwd=self._getPath())

    def createOutputStep(self):
      outSystem = self.inputSystem.get().clone()
      if not outSystem.hasTrajectory() and self.addMinimization.get():
        outSystem.setFileName(self._getPath(f'{self.getSystemName()}.pdb'))
      elif outSystem.hasTrajectory():
          repFile = outSystem.getReportFile()

          data = np.loadtxt(repFile, delimiter=',')
          cEs, ljEs = self.parseEnergies(self.getOutputFile())
          cEs, ljEs = np.array(cEs).reshape(-1, 1), np.array(ljEs).reshape(-1, 1)
          data = np.hstack((data, cEs, ljEs))

          newRepFile = self._getPath('md_log.txt')
          headerStr = self.getHeaderStr(repFile) + f',"Coulomb Energy (KJ/mol)","LJ Energy (KJ/mol)"'
          np.savetxt(newRepFile, data, delimiter=",", comments="", fmt="%f", header=headerStr)
          outSystem.setReportFile(newRepFile)

      self._defineOutputs(outputSystem=outSystem)

####################### UTILS FUNCTIONS ############################

    def parseEnergies(self, resFile):
      with open(resFile) as f:
        coulombEnergies = [float(energy) for energy in f.readline().split(':')[1].split()]
        ljEnergies = [float(energy) for energy in f.readline().split(':')[1].split()]
      return coulombEnergies, ljEnergies

    def getParamsFile(self):
      return os.path.abspath(self._getExtraPath('simulationParams.txt'))

    def getStructureFile(self):
      return os.path.abspath(self.inputSystem.get().getFileName())

    def getSystemFile(self):
      return os.path.abspath(self.inputSystem.get().getSerieFile())

    def getSystemName(self):
      return self.inputSystem.get().getSystemName()

    def getSystemTrajFile(self):
      return self.inputSystem.get().getTrajectoryFile()

    def getSystemFF(self):
      return self.inputSystem.get().getForceField()

    def getHeaderStr(self, repFile):
      with open(repFile) as f:
        headerLine = f.readline().strip()
      return headerLine

####################### SUMMARY FUNCTIONS ############################

    def getSummaryStr(self, coulombEnergies, ljEnergies):
      avg_co, avg_lj = np.mean(coulombEnergies), np.mean(ljEnergies)
      if len(ljEnergies) > 1:
        std_co, std_lj = np.std(coulombEnergies), np.std(ljEnergies)
      else:
        std_co, std_lj = 0, 0

      energyStr = f'Average Coulomb energy:\t\t{avg_co:.4f} ± {std_co:.4f} kJ/mol\n' \
                  f'Average LJ energy:\t\t{avg_lj:.4f} ± {std_lj:.4f} kJ/mol"\n'
      return energyStr

    def getOutputFile(self):
      return self._getPath('energy_results.tsv')

    def _summary(self):
      s = []
      resFile = self.getOutputFile()
      if os.path.exists(resFile):
          s = self.getSummaryStr(*self.parseEnergies(resFile))
      return s
