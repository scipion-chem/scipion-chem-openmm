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

import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

from ..protocols import ProtOpenMMReceptorPrep, ProtOpenMMSystemPrep, ProtOpenMMSystemSimulation

class TestOpenMMPrepareReceptor(BaseTest):
  @classmethod
  def setUpClass(cls):
    cls.ds = DataSet.getDataSet('model_building_tutorial')
    setupTestProject(cls)
    cls._runImportPDB()
    cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

  @classmethod
  def _runImportPDB(cls):
    cls.protImportPDB = cls.newProtocol(
      ProtImportPdb,
      inputPdbData=1,
      pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_mut1.pdb'))
    cls.proj.launchProtocol(cls.protImportPDB, wait=False)

  @classmethod
  def _runPrepareReceptor(cls):
    protPrepare = cls.newProtocol(
      ProtOpenMMReceptorPrep,
      inputAtomStruct=cls.protImportPDB.outputPdb,
      addRes=True)

    cls.launchProtocol(protPrepare)
    return protPrepare

  def test(self):
    protPrepare = self._runPrepareReceptor()
    self._waitOutput(protPrepare, 'outputStructure', sleepTime=10)
    self.assertIsNotNone(getattr(protPrepare, 'outputStructure', None))


class TestOpenMMPrepareSystem(TestOpenMMPrepareReceptor):
    @classmethod
    def _runPrepareSystem(cls, protPrepare):
        protPrepareS = cls.newProtocol(
            ProtOpenMMSystemPrep,
            inputStructure=protPrepare.outputStructure)

        cls.launchProtocol(protPrepareS)
        return protPrepareS

    def test(self):
        protPrepareRec = self._runPrepareReceptor()
        self._waitOutput(protPrepareRec, 'outputStructure', sleepTime=10)
        protPrepare = self._runPrepareSystem(protPrepareRec)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))


class TestOpenMMSimulation(TestOpenMMPrepareSystem):
  @classmethod
  def _runSimulation(cls, protPrepareS):
    protSim = cls.newProtocol(
      ProtOpenMMSystemSimulation,
      inputSystem=protPrepareS.outputSystem,
      maxIter=1000, nSteps=1000)

    cls.launchProtocol(protSim)
    return protSim

  def test(self):
    protPrepareRec = self._runPrepareReceptor()
    self._waitOutput(protPrepareRec, 'outputStructure', sleepTime=10)
    protPrepare = self._runPrepareSystem(protPrepareRec)
    self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

    protSim = self._runSimulation(protPrepare)
    self._waitOutput(protSim, 'outputSystem', sleepTime=10)
    self.assertIsNotNone(getattr(protSim, 'outputSystem', None))
