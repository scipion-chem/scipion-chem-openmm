# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
Protocol to strip water from an OpenMMSystem produced by ProtOpenMMSystemSimulation.
"""
import os
from pyworkflow.protocol import params
from pyworkflow.utils import Message
import subprocess
from pwem.protocols import EMProtocol

from .. import Plugin
from ..constants import OPENMM_DIC
from ..objects import OpenMMSystem

class ProtStripWater(EMProtocol):
    _label = 'strip water from OpenMM system'

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        g = form.addGroup('Input')
        g.addParam('inputSystem', params.PointerParam, label='Input system: ', pointerClass='OpenMMSystem',
                   help='OpenMMSystem produced by the simulation protocol.', allowsNull=False, important=True)

        g = form.addGroup('Options')
        g.addParam('keepIons', params.BooleanParam, default=True,
                   label='Keep ions: ', help='If True, ions are preserved; otherwise they are removed with water.')

    def _insertAllSteps(self):
        self._insertFunctionStep("stripStep")
        self._insertFunctionStep("createOutputStep")

    def stripStep(self):
        pdbFile = self.inputSystem.get().getFileName()
        name = os.path.splitext(os.path.basename(pdbFile))[0]
        systemPath = os.path.dirname(pdbFile)
        pdbIn = pdbFile
        dcdIn = os.path.join(systemPath, f"{name}.dcd")

        paramsFile = self._getExtraPath('strip_params.txt')
        outPrefix = self._getPath(f'{name}_clean')
        with open(paramsFile, 'w') as f:
            f.write(f'pdbIn :: {os.path.abspath(pdbIn)}\n')
            f.write(f'dcdIn :: {os.path.abspath(dcdIn)}\n')
            f.write(f'outPrefix :: {os.path.abspath(outPrefix)}\n')
            f.write(f'keepIons :: {self.keepIons.get()}\n')

        Plugin.runScript(self, 'stripWater.py', args=os.path.abspath(paramsFile), env=OPENMM_DIC, cwd=self._getPath())

    def createOutputStep(self):
        pdbFile = self.inputSystem.get().getFileName()
        name = os.path.splitext(os.path.basename(pdbFile))[0]
        suffix = '_clean'
        outPdb = self._getPath(f'{name}{suffix}.pdb')
        outDcd = self._getPath(f'{name}{suffix}.dcd')

        outCif = self._getPath(f'{name}{suffix}.cif')
        subprocess.run(["obabel", "-ipdb", outPdb, "-ocif", "-O", outCif], check=True)

        mFF, wFF = self.getFFFiles()
        nFrames = self.inputSystem.get().getNFrames()
        nTime = self.inputSystem.get().getNTime()

        outSystem = OpenMMSystem(filename=outPdb, serieFile=outDcd if os.path.exists(outDcd) else '',
                                 cifFile=outCif, ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)

        if os.path.exists(outDcd):
            outSystem.setTrajectoryFile(outDcd)

        self._defineOutputs(outputSystem=outSystem)

    def getFFFiles(self):
        system = self.inputSystem.get()
        return system.getForceField(), system.getWaterForceField()