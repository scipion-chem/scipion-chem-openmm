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

from pwchem.tests import TestPrepareReceptor

from ..protocols import ProtOpenMMSystemPrep, ProtOpenMMSystemSimulation

class TestOpenMMPrepareSystem(TestPrepareReceptor):
    @classmethod
    def _runPrepareSystem(cls, protPrepare):
        protPrepareS = cls.newProtocol(
            ProtOpenMMSystemPrep,
            inputStructure=protPrepare.outputStructure)

        cls.launchProtocol(protPrepareS)
        return protPrepareS

    def test(self):
        self._runPrepareReceptor()
        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
        protPrepare = self._runPrepareSystem(self.protPrepareReceptor)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))


class TestOpenMMSimulation(TestOpenMMPrepareSystem):
  @classmethod
  def _runSimulation(cls, protPrepareS):
    protSim = cls.newProtocol(
      ProtOpenMMSystemSimulation,
      inputSystem=protPrepareS.outputSystem,
      maxIter=50, nSteps=100)

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
