# **************************************************************************
# *
# * Authors:     Joaquin Algorta (joaquin.algorta@cnb.csic.es)
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

from pwem.protocols import ProtImportPdb
from pyworkflow.tests import DataSet
from pwchem.tests import TestExtractLigand

from ..protocols import ProtOpenFERBFE, ProtOpenFEABFE

# 1uaz chain A (retinal), as in TestGromacsPmxABFE.
chainStr = '{"model": 0, "chain": "A", "residues": 236}'
# JNK1 chain A, residues 1-358 - real residue range of the FEP/ PDBs below.
jnk1ChainStr = '{"model": 0, "chain": "A", "residues": 358}'


class TestOpenFERBFE(TestExtractLigand):
    """Smoke test for the OpenFE relative binding free energy (RBFE) protocol.

    `test` runs one published benchmark edge, JNK1 18625-1 -> 18626-1 from pmx's
    protLig_benchmark: congeneric enough for LOMAP to map (one chlorine moved). One edge and one
    repeat, so it cannot assert on RBFE_dG.

    `test2` is the only test producing a per-ligand DG, so it must satisfy both openfe minimums
    for `gather --report dg`: at least 3 edges (hence 4 ligands, N-1 edges) and at least 2
    repeats. Missing either leaves the RBFE_dG column empty."""

    @classmethod
    def _runImportPDB(cls):
        """setUpClass calls this, so import the 2-ligand JNK1 PDB `test` needs rather than the
        unrelated structure the base class would fetch."""
        cls.protImportPDB = cls._runImportJNK1PDB()

    @classmethod
    def _runImportJNK1PDB(cls, quad=False):
        pdbName = ('FEP/jnk1_18624-1_18625-1_18626-1_18627-1.pdb' if quad
                   else 'FEP/jnk1_18625-1_18626-1.pdb')
        dsLig = DataSet.getDataSet('smallMolecules')
        protImport = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                     pdbFile=dsLig.getFile(pdbName))
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def _runOpenFERBFE(cls, protExtract, label, repeats=1):
        protRBFE = cls.newProtocol(
            ProtOpenFERBFE,
            protocolRepeats=repeats, nReplicas=11, productionLength=0.01, equilLength=0.005)
        protRBFE.inputSetOfMols.set(protExtract)
        protRBFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protRBFE.setObjLabel(label)

        cls.launchProtocol(protRBFE)
        return protRBFE

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB, jnk1ChainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protRBFE = self._runOpenFERBFE(protExtract, 'openfe - RBFE (JNK1 18625-1 -> 18626-1)')
        self._waitOutput(protRBFE, 'outputSmallMolecules', sleepTime=10)
        outMols = getattr(protRBFE, 'outputSmallMolecules', None)
        self.assertIsNotNone(outMols)
        self.assertEqual(len(outMols), 2)
        # One edge, so no DG table exists - this proves the edge ran and gave a relative dG.
        self.assertEqual(len(protRBFE.parseEdges()), 1)
        self.assertTrue(os.path.exists(protRBFE.getGatherFile('ddg')))

    def test2(self):
        """4 ligands -> a 3-edge minimal spanning network, the smallest that yields DG values."""
        protImportQuad = self._runImportJNK1PDB(quad=True)
        self._waitOutput(protImportQuad, 'outputPdb')

        protExtract = self._runExtractLigand(protImportQuad, jnk1ChainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        # repeats=2, not the 1 every other test uses: gather refuses a single repeat, so this is
        # the only way to exercise the DG column. It roughly doubles the runtime.
        protRBFE = self._runOpenFERBFE(protExtract, 'openfe - RBFE (4-ligand network, JNK1)',
                                       repeats=2)
        self._waitOutput(protRBFE, 'outputSmallMolecules', sleepTime=10)
        outMols = getattr(protRBFE, 'outputSmallMolecules', None)
        self.assertIsNotNone(outMols)
        self.assertEqual(len(outMols), 4)
        self.assertEqual(len(protRBFE.parseEdges()), 3)          # N ligands -> N-1 edges
        # 3 edges clears openfe's minimum, so every ligand gets a DG.
        self.assertTrue(os.path.exists(protRBFE.getGatherFile('dg')))
        for mol in outMols:
            self.assertIsNotNone(mol.RBFE_dG.get())


class TestOpenFEABFE(TestExtractLigand):
    """Smoke test for the OpenFE absolute binding free energy (ABFE) protocol.

    `test` runs retinal in 1uaz, the same single-ligand system TestGromacsPmxABFE uses. The
    protocol takes ONE ligand, wizard-picked, so there is no multi-ligand test.

    preEquilLength is essential here: openfe's per-leg pre-equilibration defaults total 6.55 ns
    per repeat and no window setting touches them, so without shrinking it the run spends hours
    before the first alchemical window."""

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(ProtImportPdb, inputPdbData=0, pdbId='1uaz')
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runOpenFEABFE(cls, protExtract, ligandName, label):
        protABFE = cls.newProtocol(
            ProtOpenFEABFE,
            protocolRepeats=1, productionLength=0.01, equilLength=0.005, preEquilLength=0.01)
        protABFE.inputSetOfMols.set(protExtract)
        protABFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protABFE.inputLigand.set(ligandName)
        protABFE.setObjLabel(label)

        cls.launchProtocol(protABFE)
        return protABFE

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB, chainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        ligandName = str(next(iter(protExtract.outputSmallMolecules)))
        protABFE = self._runOpenFEABFE(protExtract, ligandName, 'openfe - ABFE (1uaz RET)')
        self._waitOutput(protABFE, 'outputSmallMolecules', sleepTime=10)
        outMols = getattr(protABFE, 'outputSmallMolecules', None)
        self.assertIsNotNone(outMols)
        self.assertEqual(len(outMols), 1)        # ABFE runs exactly one ligand
        self.assertIsNotNone(next(iter(outMols)).ABFE_dG.get())
