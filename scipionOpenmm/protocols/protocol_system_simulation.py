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
from pwem.convert import AtomicStructHandler

from pwchem.utils import runOpenBabel, getBaseName
from pwchem import Plugin as pwchemPlugin

from scipionOpenmm import Plugin as openmmPlugin
from scipionOpenmm.constants import OPENMM_DIC
from scipionOpenmm.objects import OpenMMSystem


class ProtOpenMMSystemSimulation(EMProtocol):
    """
    This protocol will start a Molecular Dynamics simulation.
    """
    _label = 'system simulation'


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSystem', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='OpenMMSystem', help='OpenMMSystem to execute the simulation over')
        form.addParam('nSteps', params.IntParam, default=10000, label="Number of simualtion steps: ",
                      help='Number of steps for simulation')

        cGroup = form.addGroup('Constraints')
        cGroup.addParam('constraints', params.EnumParam, default=0, label="Constraints: ",
                        choices=['None', 'HBonds', 'AllBonds', 'HAngles'],
                        help='http://docs.openmm.org/latest/userguide/application/02_running_sims.html#constraints')

        mGroup = form.addGroup('Minimization')
        mGroup.addParam('addMinimization', params.BooleanParam, default=False, label="Add minimization: ",
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
      outDir = os.path.abspath(self._getExtraPath())
      inFile = self.getSystemFilename()
      systemBasename = self.getSystemName()

      with open(self.getParamsFile(), 'w') as f:
        f.write('inputFile :: {}\n'.format(inFile))
        mFF, wFF = self.getFFFiles()
        f.write('mFF :: {}\nwFF :: {}\n'.format(mFF, wFF))
        f.write('nSteps :: {}\n'.format(self.nSteps.get()))

        f.write('constraints :: {}\n'.format(self.getEnumText('constraints')))

        nbMethod, nbCutOff = self.getNBParams()
        f.write('nbMethod :: {}\nnbCutoff :: {}\n'.format(nbMethod, nbCutOff))

        integrator = self.getEnumText('integrator')
        f.write('integrator :: {}\n'.format(integrator))
        if not self.integrator.get() in [0, 5]:
          f.write('temperature :: {}\n'.format(self.temperature.get()))

        if not self.integrator.get() in [5, 6]:
          f.write('stepSize :: {}\n'.format(self.stepSize.get()))

        if not self.integrator.get() in [0, 3, 5]:
          f.write('fricCoef :: {}\n'.format(self.fricCoef.get()))

        f.write('addMinimization :: {}\n'.format(self.addMinimization.get()))
        if self.addMinimization:
          f.write('minimTol :: {}\n'.format(self.minimTol.get()))
          f.write('maxIter :: {}\n'.format(self.maxIter.get()))

        f.write('addBarostat :: {}\n'.format(self.addBarostat.get()))
        if self.addBarostat:
          f.write('pressure :: {}\n'.format(self.pressure.get()))
          f.write('temperature :: {}\n'.format(self.temperature.get()))

      openmmPlugin.runScript(self, 'openmmSimulateSystem.py', args=self.getParamsFile(), env=OPENMM_DIC,
                             cwd=self._getPath())


    def createOutputStep(self):
      systemBasename = self.getSystemName()
      outSystemFile = self._getPath('{}_system.pdb'.format(systemBasename))

      mFF, wFF = self.getFFFiles()
      outSystem = OpenMMSystem(filename=outSystemFile, ff=mFF, wff=wFF,
                               nonbondedMethod=self.getEnumText('nonbondedMethod'),
                               nonbondedCutoff=self.nonbondedCutoff.get())

      self._defineOutputs(outputSystem=outSystem)
      self._defineSourceRelation(self.inputStructure, outSystem)


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

    def getSystemFilename(self):
      return os.path.abspath(self.inputSystem.get().getFileName())

    def getSystemName(self):
      return getBaseName(self.getSystemFilename())
