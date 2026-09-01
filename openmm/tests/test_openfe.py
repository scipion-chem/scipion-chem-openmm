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

    `test` runs a real published RBFE benchmark edge: ligands 18625-1 -> 18626-1 against JNK1,
    from pmx's protLig_benchmark. The two are genuinely congeneric (RDKit Morgan Tanimoto 0.77 -
    one chlorine moved to a different ring position), which is what LOMAP needs to map them at
    all. Experimental ddG is -3.22 kJ/mol (-0.77 kcal/mol); not asserted, since the settings
    below are tiny.

    `test2` gives the protocol three ligands so openfe's own generate_minimal_spanning_network
    plans a real 2-edge network: both edges land in one results.txt and no single scalar ddG is
    published."""

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(ProtImportPdb, inputPdbData=0, pdbId='1uaz')
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runImportJNK1PDB(cls, trio=False):
        pdbName = ('FEP/jnk1_18624-1_18625-1_18626-1.pdb' if trio
                   else 'FEP/jnk1_18625-1_18626-1.pdb')
        dsLig = DataSet.getDataSet('smallMolecules')
        protImport = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                     pdbFile=dsLig.getFile(pdbName))
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def _runOpenFERBFE(cls, protExtract, label):
        protRBFE = cls.newProtocol(
            ProtOpenFERBFE,
            protocolRepeats=1, nReplicas=11, productionLength=0.01, equilLength=0.005)
        protRBFE.inputSetOfMols.set(protExtract)
        protRBFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protRBFE.setObjLabel(label)

        cls.launchProtocol(protRBFE)
        return protRBFE

    def test(self):
        protImportJNK1 = self._runImportJNK1PDB()
        self._waitOutput(protImportJNK1, 'outputPdb')

        protExtract = self._runExtractLigand(protImportJNK1, jnk1ChainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protRBFE = self._runOpenFERBFE(protExtract, 'openfe - RBFE (JNK1 18625-1 -> 18626-1)')
        self._waitOutput(protRBFE, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protRBFE, 'outputSystem', None))
        # Two ligands -> one network edge -> the scalar ddG is published on the output.
        self.assertIsNotNone(protRBFE.outputSystem.getFreeEnergy())

    def test2(self):
        """3 ligands -> openfe plans a 2-edge minimal spanning network."""
        protImportTrio = self._runImportJNK1PDB(trio=True)
        self._waitOutput(protImportTrio, 'outputPdb')

        protExtract = self._runExtractLigand(protImportTrio, jnk1ChainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protRBFE = self._runOpenFERBFE(protExtract, 'openfe - RBFE (3-ligand network, JNK1)')
        self._waitOutput(protRBFE, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protRBFE, 'outputSystem', None))
        self.assertEqual(len(protRBFE.parseEdges()), 2)
        # A multi-edge network reports per-edge values in results.txt and no single scalar.
        self.assertIsNone(protRBFE.outputSystem.getFreeEnergy())


class TestOpenFEABFE(TestExtractLigand):
    """Smoke test for the OpenFE absolute binding free energy (ABFE) protocol.

    `test` runs retinal in 1uaz - the same single-ligand system TestGromacsPmxABFE uses. openfe
    handles the whole double-decoupling cycle itself (automatic Boresch restraints, both legs).

    `test2` runs two ligands together: ABFE needs no atom mapping and treats each ligand
    independently, so a set yields one dG_bind per ligand and no single scalar on the output.

    preEquilLength is essential to these tests: openfe's per-leg pre-equilibration defaults total
    6.55 ns per repeat (complex 0.25+0.5+5.0, solvent 0.1+0.2+0.5) and are untouched by any
    window setting, so without shrinking it the run spends hours before the first alchemical
    window. Shrinking it also shortens the trajectory write interval automatically - the complex
    leg picks its Boresch anchors from the RMSF over that trajectory, and a 10 ps
    pre-equilibration against the default 20 ps interval writes zero frames and dies with
    "XDR read error = endoffile"."""

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(ProtImportPdb, inputPdbData=0, pdbId='1uaz')
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runImportJNK1PDB(cls):
        dsLig = DataSet.getDataSet('smallMolecules')
        protImport = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                     pdbFile=dsLig.getFile('FEP/jnk1_18625-1_18626-1.pdb'))
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def _runOpenFEABFE(cls, protExtract, label):
        protABFE = cls.newProtocol(
            ProtOpenFEABFE,
            protocolRepeats=1, productionLength=0.01, equilLength=0.005, preEquilLength=0.01)
        protABFE.inputSetOfMols.set(protExtract)
        protABFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protABFE.setObjLabel(label)

        cls.launchProtocol(protABFE)
        return protABFE

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB, chainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protABFE = self._runOpenFEABFE(protExtract, 'openfe - ABFE (1uaz RET)')
        self._waitOutput(protABFE, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protABFE, 'outputSystem', None))
        # One ligand -> the scalar dG_bind is published on the output.
        self.assertIsNotNone(protABFE.outputSystem.getFreeEnergy())

    def test2(self):
        """Two ligands -> one dG_bind each, no single scalar on the output."""
        protImportJNK1 = self._runImportJNK1PDB()
        self._waitOutput(protImportJNK1, 'outputPdb')

        protExtract = self._runExtractLigand(protImportJNK1, jnk1ChainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protABFE = self._runOpenFEABFE(protExtract, 'openfe - ABFE (JNK1, 2 ligands)')
        self._waitOutput(protABFE, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protABFE, 'outputSystem', None))
        self.assertEqual(len(protABFE.parseLigandNames()), 2)
        self.assertIsNone(protABFE.outputSystem.getFreeEnergy())
