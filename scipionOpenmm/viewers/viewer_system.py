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

import os, glob
import pyworkflow.protocol.params as params

from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.viewers import ChimeraViewer, EmPlotter

from pwchem.viewers import VmdViewPopen, MDSystemViewer, MDSystemPViewer
from pwchem.utils import natural_sort
from pwchem.constants import TCL_MD_STR

from scipionOpenmm import Plugin as openmmPlugin
from scipionOpenmm.objects import OpenMMSystem
from scipionOpenmm.protocols import ProtOpenMMSystemPrep

class OpenMMSystemPViewer(MDSystemPViewer):
    """ Visualize the output of OpenMM simulation """
    _label = 'Viewer OpenMM System'
    _targets = [OpenMMSystem]

    def __init__(self, **args):
      super().__init__(**args)

    def getMDSystem(self, objType=OpenMMSystem):
        if isinstance(self.protocol, objType):
            return self.protocol
        else:
            return self.protocol.outputSystem

    def _showMdVMD(self, paramName=None):
      system = self.getMDSystem()

      outTcl = os.path.join(os.path.dirname(system.getTrajectoryFile()), 'vmdSimulation.tcl')
      systExt = os.path.splitext(system.getOriStructFile())[1][1:]
      trjExt = os.path.splitext(system.getTrajectoryFile())[1][1:]
      with open(outTcl, 'w') as f:
        f.write(TCL_MD_STR % (system.getSystemFile(), systExt, system.getTrajectoryFile(), trjExt))
      args = '-e {}'.format(outTcl)

      return [VmdViewPopen(args)]
