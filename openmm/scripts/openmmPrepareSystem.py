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

from openff.toolkit.topology import Molecule

from utils import parseParams, addMoleculesFF

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

if __name__ == "__main__":
    pDic = parseParams(sys.argv[1], sep='::')
    sysName = os.path.splitext(os.path.basename(pDic['receptorFile']))[0]
    
    recFile = pDic['receptorFile']
    parser = PDBFile if recFile.endswith('.pdb') else PDBxFile
    pdb = parser(recFile)
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
    PDBFile.writeFile(modeller.topology, modeller.positions, open(f'{sysName}_system.pdb', 'w'))
    PDBxFile.writeFile(modeller.topology, modeller.positions, open(f'{sysName}_system.cif', 'w'))

    sysKwargs = {"nonbondedMethod": eval(pDic['nonbondedMethod'])}
    sysKwargs.update({"nonbondedCutoff": float(pDic['nonbondedCutoff']) * nanometer})
    sysKwargs.update({"constraints": eval(pDic['constraints'])})
    system = forcefield.createSystem(modeller.topology, **sysKwargs)

    with open(f'{sysName}_system.xml', 'w') as output:
      output.write(XmlSerializer.serialize(system))

    with open(f'{sysName}_topology.pdb', 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)

        # Write bonds belonging to water molecules to make the loading in post-processing protocols possible
        water_names = {'HOH', 'WAT', 'TIP3P', 'TIP4P', 'SOL', 'SPC'}
        conect_dict = {}

        for bond in modeller.topology.bonds():
            a1, a2 = bond[0], bond[1]
            if a1.residue.name in water_names or a2.residue.name in water_names:
                i, j = a1.index + 1, a2.index + 1
                conect_dict.setdefault(i, []).append(j)
                conect_dict.setdefault(j, []).append(i)
        for atom_idx, connections in sorted(conect_dict.items()):
            conns = sorted(set(connections))

            for i in range(0, len(conns), 4):
                chunk = conns[i:i + 4]
                conect_line = f"CONECT{atom_idx:5d}" + "".join(f"{bonded_atom:5d}" for bonded_atom in chunk)
                f.write(conect_line + "\n")
