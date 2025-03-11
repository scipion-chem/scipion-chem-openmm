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

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import EspalomaTemplateGenerator, GAFFTemplateGenerator, SMIRNOFFTemplateGenerator

from utils import parseParams

def getGenerator(ligFF):
  if 'espaloma' in ligFF.lower():
    gen = EspalomaTemplateGenerator
  elif 'gaff' in ligFF.lower():
    gen = GAFFTemplateGenerator
  elif 'smirnoff' in ligFF.lower() or 'openff' in ligFF.lower():
    gen = SMIRNOFFTemplateGenerator
  return gen

def addLigand(modeller, ligFile):
  '''Update modeller object of receptor with the ligand topology and positions'''
  molecule = Molecule.from_file(ligFile)
  ligTop = molecule.to_topology().to_openmm()
  for residue in ligTop.residues():
    residue.name = "LIG"
  positions = molecule.conformers[0]

  ligPos = positions.m_as("nanometer") * nanometers

  modeller.add(ligTop, ligPos)
  return modeller

def addMoleculesFF(forcefield, ligFile, ligFF):
  '''Update forcefiled with Espaloma parameters for ligand'''
  molecule = Molecule.from_file(ligFile)
  generator = getGenerator(ligFF)
  tempGenerator = generator(molecules=molecule, forcefield=ligFF, cache="molecules_ff.json")
  forcefield.registerTemplateGenerator(tempGenerator.generator)
  return forcefield


if __name__ == "__main__":
    pDic = parseParams(sys.argv[1], sep='::')
    sysName = os.path.splitext(os.path.basename(pDic['receptorFile']))[0]

    pdb = PDBFile(pDic['receptorFile'])
    forcefield = ForceField(pDic['mFF'], pDic['wFF'])

    modeller = Modeller(pdb.topology, pdb.positions)
    ligFile = pDic['ligandFile'] if 'ligandFile' in pDic else None
    if ligFile:
      ligFF = pDic['ligandFF']
      modeller = addLigand(modeller, ligFile)
      forcefield = addMoleculesFF(forcefield, ligFile, ligFF)

    if eval(pDic['addH']):
      modeller.addHydrogens(forcefield, pH=float(pDic['hPH']))

    if 'boxSize' in pDic:
      bSize = list(map(float, pDic['boxSize'].split(',')))
      kwargs = {"boxSize": Vec3(bSize[0], bSize[1], bSize[2]) * nanometers}
    else:
      kwargs = {"padding": float(pDic['padDist'])}

    kwargs.update({"ionicStrength": float(pDic['saltConc']) * molar, "neutralize": eval(pDic['neutralize']),
                   "positiveIon": pDic['cationType'], "negativeIon": pDic['anionType']})

    modeller.addSolvent(forcefield, model=pDic['wModel'], **kwargs)

    # Save PDB for visualization
    PDBFile.writeFile(modeller.topology, modeller.positions,
                      open(f'{sysName}_system.pdb', 'w'))

    sysKwargs = {"nonbondedMethod": eval(pDic['nonbondedMethod'])}
    sysKwargs.update({"nonbondedCutoff": float(pDic['nonbondedCutoff']) * nanometer})
    sysKwargs.update({"constraints": eval(pDic['constraints'])})
    system = forcefield.createSystem(modeller.topology, **sysKwargs)

    with open(f'{sysName}_system.xml', 'w') as output:
      output.write(XmlSerializer.serialize(system))






