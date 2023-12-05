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
    nTraj = int(pDic['nTraj'])

    pdb = PDBFile(pDic['inputFile'])
    forcefield = ForceField(pDic['mFF'], pDic['wFF'])

    sysKwargs = {"nonbondedMethod": eval(pDic['nbMethod'])}
    sysKwargs.update({"nonbondedCutoff": float(pDic['nbCutoff']) * nanometer})
    sysKwargs.update({"constraints": eval(pDic['constraints'])})
    system = forcefield.createSystem(pdb.topology, **sysKwargs)

    if eval(pDic['addBarostat']):
      system.addForce(MonteCarloBarostat(float(pDic['pressure']) * bar, float(pDic['temperature']) * kelvin))

    intArgs = []
    intClass = eval('{}Integrator'.format(pDic['integrator']))
    if pDic['integrator'] in ['Langevin', 'LangevinMiddle', 'NoseHoover', 'Brownian', 'VariableLangevin']:
      intArgs.append(float(pDic['temperature']) * kelvin)

    if pDic['integrator'] in ['Langevin', 'LangevinMiddle', 'Brownian', 'VariableLangevin']:
      intArgs.append(float(pDic['fricCoef']) / picosecond)

    if pDic['integrator'] not in ['VariableVerlet', 'VariableLangevin']:
      intArgs.append(float(pDic['stepSize']) * picoseconds)

    integrator = intClass(*intArgs)

    properties = {}
    if 'gpus' in pDic:
      properties.update({'DeviceIndex': pDic['gpus'].strip()})
    simulation = Simulation(pdb.topology, system, integrator, platformProperties=properties)
    simulation.context.setPositions(pdb.positions)

    if eval(pDic['addMinimization']):
      print('Running {} minimization steps or until <= {} kJ/mol'.format(pDic['maxIter'], pDic['minimTol']))
      sys.stdout.flush()
      simulation.reporters.append(StateDataReporter(sys.stdout, nTraj, step=True,
                                                    potentialEnergy=True, temperature=True, volume=True))
      simulation.reporters.append(StateDataReporter("min_log.txt", nTraj, step=True,
                                                    potentialEnergy=True, temperature=True, volume=True))
      simulation.minimizeEnergy(tolerance=float(pDic['minimTol'])*kilojoules_per_mole/nanometer,
                                maxIterations=int(pDic['maxIter']))

    # Set up the reporters to report energies every 1000 steps.
    simulation.reporters.append(DCDReporter(f'{sysName}.dcd', nTraj))
    simulation.reporters.append(StateDataReporter("md_log.txt", nTraj, step=True,
                                                  potentialEnergy=True, temperature=True, volume=True))
    # run simulation
    print('Running {} steps simulation'.format(pDic['nSteps']))
    sys.stdout.flush()
    simulation.step(int(pDic['nSteps']))

    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'{sysName}.pdb', 'w'))



