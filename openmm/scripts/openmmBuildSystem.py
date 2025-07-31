#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# # -*- coding: utf-8 -*-
# # **************************************************************************
# # *
# # * Authors: Daniel Del Hoyo Gómez (ddelhoyo@cnb.csic.es)
# # *
# # *
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import sys
from openmm.app import *
from openmm import *
from openmm.unit import *

from utils import parseParams, addMoleculesFF

if __name__ == "__main__":
    pDic = parseParams(sys.argv[1], sep='::')
    sysName = os.path.splitext(os.path.basename(pDic['receptorFile']))[0]

    pdb = PDBFile(pDic['receptorFile'])
    forcefield = ForceField(pDic['mFF'], pDic['wFF'])

    ligFile = pDic['ligandFile'] if 'ligandFile' in pDic else None
    if ligFile:
      forcefield = addMoleculesFF(forcefield, ligFile, pDic['ligandFF'])

    sysKwargs = {"nonbondedMethod": eval(pDic['nonbondedMethod']),
                 "nonbondedCutoff": eval(pDic['nonbondedCutoff']) * nanometers,
                 "constraints": eval(pDic['constraints'])
                 }
    system = forcefield.createSystem(pdb.topology, **sysKwargs)

    with open(f'{sysName}_system.xml', 'w') as output:
      output.write(XmlSerializer.serialize(system))






