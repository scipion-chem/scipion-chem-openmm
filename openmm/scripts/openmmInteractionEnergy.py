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

# General imports
import sys, os
import numpy as np
import mdtraj as md

# Openmm imports
from openmm.app import *
from openmm import *
from openmm.unit import bar, kelvin, picosecond, kilojoules_per_mole, nanometer

from utils import parseParams

CATION_NAMES = ['Cs+', 'K+', 'Li+', 'Na+', 'Rb+']
ANION_NAMES = ['Cl-', 'Br-', 'F-', 'I-']
ION_NAMES = [ion[:-1] for ion in zip(CATION_NAMES, ANION_NAMES)]

def setCharmmForces(system):
		for force in system.getForces():
			if isinstance(force, NonbondedForce):
				force.setForceGroup(0)
				force.addGlobalParameter("protein_scale", 1)
				force.addGlobalParameter("ligand_scale", 1)
				for i in range(force.getNumParticles()):
					charge, sigma, epsilon = force.getParticleParameters(i)

					param = "protein_scale" if i in protein else "ligand_scale"
					force.setParticleParameters(i, 0, 0, 0)
					force.addParticleParameterOffset(param, i, charge, sigma, epsilon)
				for i in range(force.getNumExceptions()):
					p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
					force.setExceptionParameters(i, p1, p2, 0, 0, 0)
			elif isinstance(force, CustomNonbondedForce):
				force.setForceGroup(1)
				force.addInteractionGroup(protein, ligand)
			else:
				force.setForceGroup(2)
		return system

def setAmberForces(system):
		for force in system.getForces():
			if isinstance(force, NonbondedForce):
				force.setForceGroup(0)
				force.addGlobalParameter("receptor_coulomb_scale", 1)
				force.addGlobalParameter("receptor_lj_scale", 1)
				force.addGlobalParameter("ligand_coulomb_scale", 1)
				force.addGlobalParameter("ligand_lj_scale", 1)
				for i in range(force.getNumParticles()):
					charge, sigma, epsilon = force.getParticleParameters(i)
					force.setParticleParameters(i, 0, 0, 0)
					if i in protein:
						force.addParticleParameterOffset("receptor_coulomb_scale", i, charge, 0, 0)
						force.addParticleParameterOffset("receptor_lj_scale", i, 0, sigma, epsilon)
					elif i in ligand:
						force.addParticleParameterOffset("ligand_coulomb_scale", i, charge, 0, 0)
						force.addParticleParameterOffset("ligand_lj_scale", i, 0, sigma, epsilon)
				for i in range(force.getNumExceptions()):
					p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
					force.setExceptionParameters(i, p1, p2, 0, 0, 0)
			else:
				force.setForceGroup(2)
		return system

def setSystemForces(system, ff):
		if 'charmm' in ff:
				system = setCharmmForces(system)
		else:
				system = setAmberForces(system)
		return system

def coulomb_energy(context, protein_scale, ligand_scale):
		context.setParameter("protein_scale", protein_scale)
		context.setParameter("ligand_scale", ligand_scale)
		return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()

def energy(context, receptor_coulomb_scale, receptor_lj_scale, ligand_coulomb_scale, ligand_lj_scale):
		context.setParameter("receptor_coulomb_scale", receptor_coulomb_scale)
		context.setParameter("receptor_lj_scale", receptor_lj_scale)
		context.setParameter("ligand_coulomb_scale", ligand_coulomb_scale)
		context.setParameter("ligand_lj_scale", ligand_lj_scale)
		return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()

def getCharmmEnergies(context):
		total_coulomb = coulomb_energy(context, 1, 1)
		protein_coulomb = coulomb_energy(context, 1, 0)
		ligand_coulomb = coulomb_energy(context, 0, 1)

		coulomb = total_coulomb - protein_coulomb - ligand_coulomb
		lj = context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
		return coulomb, lj

def getAmberEnergies(context):
		total_coulomb = energy(context, 1, 0, 1, 0)
		receptor_coulomb = energy(context, 1, 0, 0, 0)
		ligand_coulomb = energy(context, 0, 0, 1, 0)
		total_lj = energy(context, 0, 1, 0, 1)
		receptor_lj = energy(context, 0, 1, 0, 0)
		ligand_lj = energy(context, 0, 0, 0, 1)

		coulomb = total_coulomb - receptor_coulomb - ligand_coulomb
		lj = total_lj - receptor_lj - ligand_lj
		return coulomb, lj

def getInteractionEnergies(context, ff):
		if 'charmm' in ff:
			coulomb, lj = getCharmmEnergies(context)
		else:
			coulomb, lj = getAmberEnergies(context)
		return  coulomb, lj

def buildIntegrator(pDic):
	intArgs = []
	intClass = eval('{}Integrator'.format(pDic['integrator']))
	if pDic['integrator'] in ['Langevin', 'LangevinMiddle', 'NoseHoover', 'Brownian', 'VariableLangevin']:
		intArgs.append(float(pDic['temperature']) * kelvin)

	if pDic['integrator'] in ['Langevin', 'LangevinMiddle', 'Brownian', 'VariableLangevin']:
		intArgs.append(float(pDic['fricCoef']) / picosecond)

	if pDic['integrator'] not in ['VariableVerlet', 'VariableLangevin']:
		intArgs.append(float(pDic['stepSize']) * picosecond)

	integrator = intClass(*intArgs)
	return integrator

def saveEnergies(coulombEnergies, ljEnergies):
	with open('energy_results.tsv', 'w') as f:
		coulombEnergies = [str(c._value) for c in coulombEnergies]
		ljEnergies = [str(c._value) for c in ljEnergies]

		f.write('Coulomb energies:\t' + '\t'.join(coulombEnergies) + '\n')
		f.write('LJ energies:\t' + '\t'.join(ljEnergies) + '\n')

if __name__ == "__main__":
	pDic = parseParams(sys.argv[1], sep='::')
	sysFile, pdbFile = pDic['systemFile'], pDic['structureFile']
	sysName = os.path.splitext(os.path.basename(sysFile))[0]
	pdb = PDBFile(pdbFile)
	with open(sysFile) as input:
		system = XmlSerializer.deserialize(input.read())

	if 'trajFile' in pDic:
			trajFile = pDic['trajFile']
			traj = md.load(trajFile, top=pdbFile)  # MDTraj trajectory
	elif eval(pDic['addMin']):
			integrator = buildIntegrator(pDic)
			simulation = Simulation(pdb.topology, system, integrator)
			simulation.context.setPositions(pdb.positions)

			simulation.minimizeEnergy(tolerance=float(pDic['minimTol']) * kilojoules_per_mole / nanometer,
																maxIterations=int(pDic['maxIter']))
			positions = simulation.context.getState(getPositions=True).getPositions()
			PDBFile.writeFile(simulation.topology, positions, open(f'{sysName}.pdb', 'w'))
	else:
			positions = pdb.positions

	solventNames = ['HOH'] + ION_NAMES
	solvent = set([a.index for a in pdb.topology.atoms() if a.residue.name in solventNames])
	ligand = set([a.index for a in pdb.topology.atoms() if a.residue.name in ('LIG')])
	protein = set([a.index for a in pdb.topology.atoms() if a.index not in solvent and a.index not in ligand])

	ff = pDic['mFF']
	system = setSystemForces(system, ff)

	if 'trajFile' in pDic:
			frames = [t.xyz[0] for t in traj[:]]
	else:
			frames = [positions]

	coEs, ljEs = [], []
	for i, positions in enumerate(frames):  # Subsample frames
		if len(frames) > 10 and (i+1) % (len(frames)//10) == 0:
			print(f'Iteration: {i+1}/{len(frames)}')
			sys.stdout.flush()
		integrator = buildIntegrator(pDic)
		context = Context(system, integrator)
		context.setPositions(positions)

		coulomb, lj = getInteractionEnergies(context, ff)
		coEs.append(coulomb), ljEs.append(lj)


	saveEnergies(coEs, ljEs)
	avg_co, avg_lj = np.mean(coEs), np.mean(ljEs)
	if len(frames) > 1:
		std_co, std_lj = np.std(coEs), np.std(ljEs)
	else:
		std_co, std_lj = 0, 0


