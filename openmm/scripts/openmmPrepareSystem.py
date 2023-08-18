#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Daniel Del Hoyo Gómez (ddelhoyo@cnb.csic.es)
# # # *
# # # *
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

import sys, os
from openmm.app import *
from openmm import *
from openmm.unit import *

def parseParams(paramsFile):
  paramsDic = {}
  with open(paramsFile) as f:
    for line in f:
      key, value = line.strip().split('::')
      paramsDic[key.strip()] = value.strip()
  return paramsDic


if __name__ == "__main__":
    pDic = parseParams(sys.argv[1])
    sysName = os.path.splitext(os.path.basename(pDic['inputFile']))[0]

    pdb = PDBFile(pDic['inputFile'])
    forcefield = ForceField(pDic['mFF'], pDic['wFF'])

    modeller = Modeller(pdb.topology, pdb.positions)
    if eval(pDic['addH']):
      modeller.addHydrogens(forcefield, pH=float(pDic['hPH']))

    # todo: infer model arg from wFF
    if 'boxSize' in pDic:
      bSize = list(map(float, pDic['boxSize'].split(',')))
      kwargs = {"boxSize": Vec3(bSize[0], bSize[1], bSize[2])*nanometers}
    else:
      kwargs = {"padding": float(pDic['padDist'])}

    kwargs.update({"ionicStrength": float(pDic['saltConc'])*molar, "neutralize": eval(pDic['neutralize']),
                   "positiveIon": pDic['cationType'], "negativeIon": pDic['anionType']})

    modeller.addSolvent(forcefield, model=pDic['wModel'], **kwargs)

    PDBFile.writeFile(modeller.topology, modeller.positions,
                      open('{}_system.pdb'.format(sysName), 'w'))




