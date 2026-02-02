# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pwem.protocols import ProtImportPdb
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwchem.protocols import ProtChemPrepareReceptor
from pwchem.tests import TestPrepareReceptor, TestExtractLigand

from ..protocols import ProtOpenMMSystemPrep, ProtOpenMMSystemSimulation, ProtOpenMMInteractionEnergy

STRUCTURE, LIGAND = 0, 1

class TestOpenMMPrepareSystem(TestPrepareReceptor, TestExtractLigand):
    @classmethod
    def _runImportPDB(cls):
      protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=0, pdbId='4erf')
      cls.launchProtocol(protImportPDB)
      cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareSystem(cls, protPrepare, inputFrom=STRUCTURE):
        protPrepareS = cls.newProtocol(
            ProtOpenMMSystemPrep, inputFrom=inputFrom)

        if inputFrom == STRUCTURE:
            protPrepareS.inputStructure.set(protPrepare)
            protPrepareS.inputStructure.setExtended('outputStructure')
        else:
            protPrepareS.inputSetOfMols.set(protPrepare)
            protPrepareS.inputSetOfMols.setExtended('outputSmallMolecules')
            protPrepareS.inputLigand.set('SmallMolecule (g1_4erf_0R3-1_1 molecule)')

        cls.launchProtocol(protPrepareS)
        return protPrepareS

    def test(self):
        self._runPrepareReceptor()
        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)

        protPrepare = self._runPrepareSystem(self.protPrepareReceptor)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))

    def test2(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))


class TestOpenMMSimulation(TestOpenMMPrepareSystem):
  @classmethod
  def _runSimulation(cls, protPrepareS):
      protSim = cls.newProtocol(
        ProtOpenMMSystemSimulation,
        inputSystem=protPrepareS.outputSystem, stepSize=0.002,
        maxIter=20, nSteps=10, nTraj=5)

      cls.launchProtocol(protSim)
      return protSim

  @classmethod
  def _runSimulationCPH(cls, protPrepareS):
      protSim = cls.newProtocol(
          ProtOpenMMSystemSimulation,
          inputSystem=protPrepareS.outputSystem,
          cph=True, singlePH=True, onePH=3.0, stepSize=0.002,
          maxIter=20, nSteps=10, nTraj=5)

      cls.launchProtocol(protSim)
      return protSim

  def test(self):
      self._runPrepareReceptor()
      self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
      protPrepare = self._runPrepareSystem(self.protPrepareReceptor)
      self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

      protSim = self._runSimulation(protPrepare)
      self._waitOutput(protSim, 'outputSystem', sleepTime=10)
      self.assertIsNotNone(getattr(protSim, 'outputSystem', None))

  def test2(self):
      protExtract = self._runExtractLigand(self.protImportPDB)
      self._waitOutput(protExtract, 'outputSmallMolecules')

      protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
      self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

      protSim = self._runSimulation(protPrepare)
      self._waitOutput(protSim, 'outputSystem', sleepTime=10)
      self.assertIsNotNone(getattr(protSim, 'outputSystem', None))

  def test_cph(self):
      self._runPrepareReceptor()
      self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
      protPrepare = self._runPrepareSystem(self.protPrepareReceptor)
      self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

      protSim = self._runSimulationCPH(protPrepare)
      self._waitOutput(protSim, 'outputSystem', sleepTime=10)
      self.assertIsNotNone(getattr(protSim, 'outputSystem', None))

  def test2_cph(self):
      protExtract = self._runExtractLigand(self.protImportPDB)
      self._waitOutput(protExtract, 'outputSmallMolecules')

      protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
      self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

      protSim = self._runSimulationCPH(protPrepare)
      self._waitOutput(protSim, 'outputSystem', sleepTime=10)
      self.assertIsNotNone(getattr(protSim, 'outputSystem', None))

class TestOpenMMInteractions(TestOpenMMSimulation):
    @classmethod
    def _runInteractions(cls, protIn, inputFrom=STRUCTURE):
        protInt = cls.newProtocol(
          ProtOpenMMInteractionEnergy, inputFrom=inputFrom, maxIter=500)

        if inputFrom == STRUCTURE:
            protInt.inputSystem.set(protIn)
            protInt.inputSystem.setExtended('outputSystem')
        else:
            protInt.inputSetOfMols.set(protIn)
            protInt.inputSetOfMols.setExtended('outputSmallMolecules')

        cls.launchProtocol(protInt)
        return protInt

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)

        protInt = self._runInteractions(protSim)
        self._waitOutput(protInt, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protInt, 'outputSystem', None))

    def test2(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protInt = self._runInteractions(protExtract, inputFrom=LIGAND)
        self._waitOutput(protInt, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protInt, 'outputSmallMolecules', None))
