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
import sys

# Openmm imports
from openmm.app import PDBFile, ForceField, Simulation, StateDataReporter,\
	DCDReporter, NoCutoff, HBonds, PDBxFile, CheckpointReporter
from openmm import *
from openmm.unit import *

from utils import parseParams

def getDoneSteps(statusFile):
	with open(statusFile, 'r') as file:
		lines = file.readlines()
		doneSteps = int(lines[-1].split(',')[0].strip())
	return doneSteps

if __name__ == "__main__":
	pDic = parseParams(sys.argv[1], sep='::')
	sysFile, recFile = pDic['systemFile'], pDic['structureFile']

	parser = PDBFile if recFile.endswith('.pdb') else PDBxFile
	pdb = parser(recFile)
	with open(sysFile) as input:
		system = XmlSerializer.deserialize(input.read())

	sysName = os.path.splitext(os.path.basename(sysFile))[0]
	nTraj = int(pDic['nTraj'])

	if eval(pDic['addBarostat']):
		system.addForce(MonteCarloBarostat(float(pDic['pressure']) * bar, float(pDic['temperature']) * kelvin))

	intArgs = []
	intClass = eval('{}Integrator'.format(pDic['integrator']))
	if pDic['integrator'] in ['Langevin', 'LangevinMiddle', 'NoseHoover', 'Brownian', 'VariableLangevin']:
		intArgs.append(float(pDic['temperature']) * kelvin)

	if pDic['integrator'] in ['Langevin', 'LangevinMiddle', 'Brownian', 'VariableLangevin']:
		intArgs.append(float(pDic['fricCoef']) / picosecond)

	if pDic['integrator'] not in ['VariableVerlet', 'VariableLangevin']:
		intArgs.append(float(pDic['stepSize']) * picosecond)

	integrator = intClass(*intArgs)

	kwargs = {}
	if 'gpus' in pDic:
		kwargs['platform'] = Platform.getPlatformByName('CUDA')
		kwargs['platformProperties'] = {'DeviceIndex': pDic['gpus'].strip()}

	simulation = Simulation(pdb.topology, system, integrator, **kwargs)
	simulation.context.setPositions(pdb.positions)
	statusFile = "md_log.txt"

	doneSteps = 0
	if 'chkFile' in pDic:
		chkFile = pDic['chkFile']
		if os.path.exists(chkFile):
			print("Loading checkpoint...")
			doneSteps = getDoneSteps(statusFile)
			with open(chkFile, 'rb') as f:
				simulation.loadCheckpoint(f)

	if eval(pDic['addMinimization']) and not os.path.exists(chkFile):
		print('Running {} minimization steps or until <= {} kJ/mol'.format(pDic['maxIter'], pDic['minimTol']))
		sys.stdout.flush()
		simulation.reporters.append(StateDataReporter(sys.stdout, nTraj, step=True,
													  potentialEnergy=True, temperature=True, volume=True))
		simulation.reporters.append(StateDataReporter("min_log.txt", nTraj, step=True,
													  potentialEnergy=True, temperature=True, volume=True))
		simulation.minimizeEnergy(tolerance=float(pDic['minimTol'])*kilojoules_per_mole/nanometer,
								  maxIterations=int(pDic['maxIter']))

		minPositions = simulation.context.getState(getPositions=True).getPositions()
		PDBFile.writeFile(simulation.topology, minPositions, open(f'{sysName}_minimized.pdb', 'w'))

	# Set up the reporters to report energies every 1000 steps.
	trjFile = f'{sysName}.dcd'
	appe = False
	if os.path.exists(trjFile):
		appe = True
	simulation.reporters.append(DCDReporter(trjFile, nTraj, append=appe))
	simulation.reporters.append(StateDataReporter(statusFile, nTraj, step=True, append=appe,
												  potentialEnergy=True, temperature=True, volume=True))
	simulation.reporters.append(CheckpointReporter(chkFile, nTraj))

	# run simulation
	todoSteps = int(pDic['nSteps']) - doneSteps
	print(f'Running {todoSteps} steps simulation')
	sys.stdout.flush()
	simulation.step(todoSteps)

	positions = simulation.context.getState(getPositions=True).getPositions()
	PDBFile.writeFile(simulation.topology, positions, open(f'{sysName}.pdb', 'w'))
	PDBxFile.writeFile(simulation.topology, positions, open(f'{sysName}.cif', 'w'))
