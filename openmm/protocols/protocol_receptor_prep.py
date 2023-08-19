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
This module will prepare a PDB receptor for OpenMM simulations
"""
import os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.objects import AtomStruct

from pwchem.protocols import ProtChemPrepareReceptor

from .. import Plugin


class ProtOpenMMReceptorPrep(ProtChemPrepareReceptor):
    """
    This protocol uses PDBFixer for receptor preparation (https://github.com/openmm/pdbfixer).
    It will fix the PDB file so it can latterly be used for OpenMM simulations
    """
    _label = 'receptor preparation'

    # -------------------------- DEFINE constants ----------------------------


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputAtomStruct', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='AtomStruct', help='Atom structure to convert to OpenMM system')
        form.addParam('addAtoms', params.EnumParam, default=0,
                      label="Add missing atoms: ", choices=['All', 'Heavy', 'Hydrogen', 'None'],
                      help='Use PDBFixer to add the missing atoms specified in the PDB atomic structure')
        form.addParam('addRes', params.BooleanParam, default=False, label="Add missing residues: ",
                      help='Use PDBFixer to add missing residues')
        form.addParam('repNonStd', params.BooleanParam, default=False, label="Replace non-standard residues: ",
                      help='Use PDBFixer to replace nonstandard residues with standard equivalents')

        self.defineCleanParams(form)

    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('pdbFixerStep')
        self._insertFunctionStep('createOutput')

    def pdbFixerStep(self):
        pdbFile = os.path.abspath(self.getPreparedFile())
        addResStr = ' --add-residues' if self.addRes else ''
        repNStdStr = ' --replace-nonstandard' if self.repNonStd else ''
        Plugin.runOpenMM(self, 'pdbfixer', args='{} --add-atoms={}{}{} --output {}'.
                               format(pdbFile, self.getEnumText('addAtoms').lower(), addResStr, repNStdStr, pdbFile))

    def createOutput(self):
        fnOut = self.getPreparedFile()
        if os.path.exists(fnOut):
            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputAtomStruct, target)

