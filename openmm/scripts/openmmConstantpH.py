#!/usr/bin/env python3
"""
Parametrizable Constant pH Simulation in OpenMM (TXT config)
Supports ASP, GLU, CYS, HIS (3-state) and LYS

This version has extra prints for progress/debugging.
"""

import argparse
import time
import math

from openmm import *
from openmm.app import *
from openmm.unit import *
import os, sys

from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField as offForceField
from openff.toolkit.utils import get_data_file_path

# ---------------------------
# Helpers
# ---------------------------

def parsePhValues(singlePH, onePH, manyPH):
    if singlePH:
        return [onePH]
    else:
        return [float(x.strip()) for x in manyPH.split(',')]


def parseTxtConfig(filename):
    params = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            k, v = [x.strip() for x in line.split('=', 1)]

            if v.lower() in ['true', 'false']:
                params[k] = v.lower() == 'true'
            else:
                try:
                    # try numeric
                    params[k] = float(v) if '.' in v else int(v)
                except ValueError:
                    # list or string
                    params[k] = [x.strip() for x in v.split(',')] if ',' in v else v
    return params

def generateLigandFF(ligFile, ligFF):
    sdf = get_data_file_path(ligFile)
    ligand = Molecule(sdf)

    if ligand.n_conformers == 0:
        print("[INFO] No conformers found, generating 3D conformer...")
        ligand.generate_conformers(n_conformers=1)

    offFF = offForceField("openff-2.1.0.offxml")
    system = offFF.create_openmm_system(ligand.to_topology())

    with open(ligFF, "w") as f:
        f.write(XmlSerializer.serialize(system))

    print(f"[INFO] Ligand XML force field saved as: {ligFF}")

# ---------------------------
# Integrator factory
# ---------------------------
def createIntegrator(params, temperature):
    stepSize = params.get('stepSize', 0.004) * picoseconds
    fric = params.get('fricCoef', 1.0) / picosecond
    colFreq = params.get('colFreq', 1.0) / picosecond
    errTol = params.get('errTol', 0.001)

    integratorName = params.get('integrator', 'Langevin')

    print(f"[integrator] Creating integrator '{integratorName}' stepSize={stepSize}, temp={temperature}")
    if integratorName == 'Verlet':
        return VerletIntegrator(stepSize)
    elif integratorName == 'Langevin':
        return LangevinIntegrator(temperature, fric, stepSize)
    elif integratorName == 'LangevinMiddle':
        return LangevinMiddleIntegrator(temperature, fric, stepSize)
    elif integratorName == 'NoseHoover':
        return NoseHooverIntegrator(temperature, 1.0 / picosecond, stepSize)
    elif integratorName == 'Brownian':
        return BrownianIntegrator(temperature, stepSize)
    elif integratorName == 'VariableVerlet':
        return VariableVerletIntegrator(stepSize, errTol)
    elif integratorName == 'VariableLangevin':
        return VariableLangevinIntegrator(temperature, fric, stepSize, errTol)
    else:
        raise ValueError(f"Unknown integrator: {integratorName}")


def computeRef(modelFile, variantsDict, targetPKa, params,
                explicitFF, implicitFF, explicitParams, implicitParams,
                integrator, relaxationIntegrator):
    """
    Compute reference energies for a model residue (ASP, GLU, etc.) for constant pH simulation.

    Returns a dict: {index: [energy_state0, energy_state1, ...]}
    """
    pdb = PDBFile(modelFile)

    cph = ConstantPH(
        pdb.topology,
        pdb.positions,
        [7.0],
        explicitFF,
        implicitFF,
        variantsDict,
        {index: [0.0] * len(states) for index, states in variantsDict.items()},
        250,
        explicitParams,
        implicitParams,
        integrator,
        relaxationIntegrator
    )

    # Compute reference energies
    finder = ReferenceEnergyFinder(cph, targetPKa, params.get('temperature', 300) * kelvin)

    totalIterations = params['equilSteps'] * params['stepEquil']
    chunk = params['relaxSteps']

    startTime = time.time()
    for start in range(0, totalIterations, chunk):
        finder.findReferenceEnergies(iterations=chunk, substeps=10)
        try:
            pos = cph.simulation.context.getState(getPositions=True).getPositions()
            for p in pos[:10]:
                if math.isnan(p.x) or math.isnan(p.y) or math.isnan(p.z):
                    raise ValueError("NaN in positions during reference energy computation")
        except Exception as e:
            print(f"[computeRef] ERROR while checking positions: {e}")
            raise

    elapsed = time.time() - startTime

    # Extract reference energies
    refenergies = {index: cph.titrations[index].referenceEnergies for index in variantsDict}
    return refenergies


# ---------------------------
# Main simulation function
# ---------------------------

def runConstantPhSimulation(params):

    print("\n--- Loading system ---")
    print(f"[params] inputPdb: {params.get('inputPdb')}")
    pdb = PDBFile(params['inputPdb'])

    print("[run] Creating force fields...")
    #explicitFF = ForceField(params['explicitFF'], params['explicitSolvent'])
    #todo this includes ligand ff
    explicitFF = ForceField(params['explicitFF'], params['explicitSolvent'], params['ligandFF'])
    print("  explicit ForceField created.")
    implicitFF = ForceField(params['implicitFF'], params['implicitSolvent'])
    print("  implicit ForceField created.")

    explicitParams = dict( #todo this was PME why tf does it not work!?
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=params['explicitCutoff'] * nanometers,
        constraints=params['constraints'],
        hydrogenMass=params['hydrogenMass'] * amu
    )

    implicitParams = dict(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=params['implicitCutoff'] * nanometers,
        constraints=params['constraints']
    )

    print(f"[run] NB params explicit_cutoff={explicitParams['nonbondedCutoff']}, implicit_cutoff={implicitParams['nonbondedCutoff']}")

    temperature = params['temperature'] * kelvin

    # ---------------------------
    # Create integrators
    # ---------------------------
    integrator = createIntegrator(params, temperature)
    relaxationIntegrator = createIntegrator(params, temperature)  # Can have different params if needed
    print("[run] Integrators created.")

    # -----------------------------------
    # Reference energies
    # -----------------------------------

    print("\n--- Computing reference energies ---")
    refenergies = {}
    variantsDict = {}

    # ASP
    if 'ASP' in params['residuesToTitrate']:
        print("[run] Computing ASP reference energy...")
        refenergies['ASP'] = computeRef(
            params['aspModel'],
            {1: ['ASP', 'ASH']},
            3.9,
            params,
            explicitFF, implicitFF,
            explicitParams, implicitParams,
            integrator, relaxationIntegrator
        )[1]
        variantsDict['ASP'] = ['ASP', 'ASH']

    # GLU
    if 'GLU' in params['residuesToTitrate']:
        print("[run] Computing GLU reference energy...")
        refenergies['GLU'] = computeRef(
            params['gluModel'],
            {1: ['GLU', 'GLH']},
            4.2,
            params,
            explicitFF, implicitFF,
            explicitParams, implicitParams,
            integrator, relaxationIntegrator
        )[1]
        variantsDict['GLU'] = ['GLU', 'GLH']

    # CYS
    if 'CYS' in params['residuesToTitrate']:
        print("[run] Computing CYS reference energy...")
        refenergies['CYS'] = computeRef(
            params['cysModel'],
            {1: ['CYS', 'CYX']},
            7.1,
            params,
            explicitFF, implicitFF,
            explicitParams, implicitParams,
            integrator, relaxationIntegrator
        )[1]
        variantsDict['CYS'] = ['CYS', 'CYX']

    # HIS (3 states)
    if 'HIS' in params['residuesToTitrate']:
        print("[run] Computing HIS reference energies (HID/HIE)...")
        hid = computeRef(
            params['hisModel'],
            {1: ['HIP', 'HID']},
            7.1,
            params,
            explicitFF, implicitFF,
            explicitParams, implicitParams,
            integrator, relaxationIntegrator
        )[1]

        hie = computeRef(
            params['hisModel'],
            {1: ['HIP', 'HIE']},
            6.5,
            params,
            explicitFF, implicitFF,
            explicitParams, implicitParams,
            integrator, relaxationIntegrator
        )[1]

        refenergies['HIS'] = [
            0.0 * kilojoules_per_mole,
            hid[1],
            hie[1]
        ]
        variantsDict['HIS'] = ['HIP', 'HID', 'HIE']

    # LYS
    if 'LYS' in params['residuesToTitrate']:
        print("[run] Computing LYS reference energy...")
        refenergies['LYS'] = computeRef(
            params['lysModel'],
            {1: ['LYS', 'LYN']},
            10.5,
            params,
            explicitFF, implicitFF,
            explicitParams, implicitParams,
            integrator, relaxationIntegrator
        )[1]
        variantsDict['LYS'] = ['LYS', 'LYN']

    # -----------------------------------
    # Assign residues
    # -----------------------------------
    phValues = parsePhValues(
        params['singlePH'], params['onePH'], params['manyPH']
    )
    print(f"[run] pH values to run: {phValues}")

    simVariants = {}
    simRefenergies = {}

    for residue in pdb.topology.residues():
        # Ignorar ligandos LIG
        if residue.name == 'LIG':
            print(f"[run] Ignoring ligand residue {residue.name} at index {residue.index}")
            continue
        if residue.name in variantsDict:
            simVariants[residue.index] = variantsDict[residue.name]
            simRefenergies[residue.index] = refenergies[residue.name]

    print("Titrated residues:")
    for k, v in simVariants.items():
        print("  Residue", k, "->", v)

    # -----------------------------------
    # Build ConstantPH Simulation
    # -----------------------------------
    #todo to test, ideally i want to keep ligands
    modeller = Modeller(pdb.topology, pdb.positions)

    # Remove ligands
    ligands = [res for res in modeller.topology.residues() if res.name == 'LIG']
    if ligands:
        print(f"[run] Removing {len(ligands)} LIG residues")
        modeller.delete(ligands)

    # Use the filtered topology and positions in ConstantPH
    filteredTopology = modeller.topology
    filteredPositions = modeller.positions
    #todo try to see if it works with normal pdb
    print("\n--- Creating simulation ---")
    cph = ConstantPH(
        pdb.topology, pdb.positions, phValues,
        explicitFF, implicitFF,
        simVariants, simRefenergies,
        params['relaxSteps'],
        explicitParams, implicitParams,
        integrator, relaxationIntegrator
    )
    print("[run] ConstantPH object created.")

    trajFile = params['trajFile']
    logFile = params['logFile']
    finalPdb = params['finalPdb']
    reportEvery = 1000

    cph.simulation.reporters.append(DCDReporter(trajFile, reportEvery))
    totalSteps = params['equilSteps'] * params['stepEquil'] + params['prodSteps'] * params['stepProd']
    cph.simulation.reporters.append(
        StateDataReporter(
            logFile,
            reportEvery,
            step=True,
            temperature=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            volume=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=totalSteps,
            separator='\t'
        )
    )

    if params.get('addBarostat', False):
        print("Adding barostat...")
        cph.simulation.system.addForce(
            MonteCarloBarostat(params['pressure'] * bar, temperature)
        )
        cph.simulation.context.reinitialize(preserveState=True)
        print("[run] Barostat added and context reinitialized.")

    if params.get('addMinimization', True):
        print("Minimizing...")
        cph.simulation.minimizeEnergy(
            tolerance=params['minimTol'] * kilojoules_per_mole / nanometer ,
            maxIterations=params['maxIter']
        )
        print("[run] Minimization finished.")

    # -----------------------------------
    # Equilibration
    # -----------------------------------

    print("\n--- Equilibration ---")
    for idx in range(params['equilSteps']):
        cph.simulation.step(params['stepEquil'])
        cph.attemptMCStep(temperature)
        if idx % max(1, params.get('equilSteps') // 10) == 0:
            print(f"[equil] Completed equilibration cycle {idx+1}/{params['equilSteps']}")

    # -----------------------------------
    # Production
    # -----------------------------------

    print("\n--- Production ---")
    for prodIdx in range(params['prodSteps']):
        cph.simulation.step(params['stepProd'])
        cph.attemptMCStep(temperature)

        try:
            state = cph.simulation.context.getState(getPositions=True)
            positions = state.getPositions()
            for i, p in enumerate(positions[:20]):
                if math.isnan(p.x) or math.isnan(p.y) or math.isnan(p.z):
                    raise ValueError(f"NaN detected in position at prod step {prodIdx}, atom {i}")
        except Exception as e:
            print(f"[production] ERROR: {e}")
            raise

        if prodIdx % max(1, params.get('prodSteps') // 10) == 0:
            states = [simVariants[i][cph.titrations[i].currentIndex] for i in simVariants]
            print(f"[production] step {prodIdx+1}/{params['prodSteps']} pH: {cph.pH[cph.currentPHIndex]} states: {states}")

    print("[run] Production finished successfully.")

    state = cph.simulation.context.getState(getPositions=True)
    with open(finalPdb, "w") as f:
        PDBFile.writeFile(cph.simulation.topology, state.getPositions(), f)

    finalCif = params['finalCif']
    with open(finalCif, "w") as f:
        PDBxFile.writeFile(cph.simulation.topology, state.getPositions(), f)

    print("[run] Final snapshot written.")


# ---------------------------
# Entry point
# ---------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True,
                        help="TXT configuration file")

    args = parser.parse_args()

    params = parseTxtConfig(args.params)

    constantphPath = params.get('constantPHScript')
    referenceEnergyPath = params.get('referenceEnergyScript')

    if constantphPath:
        scriptsDir = os.path.dirname(constantphPath)
        if scriptsDir and scriptsDir not in sys.path:
            sys.path.insert(0, scriptsDir)

    try:
        from constantph import ConstantPH
        from reference_energy import ReferenceEnergyFinder
    except Exception as e:
        print(f"[startup] Could not import constantph/reference_energy: {e}")
        print("[startup] Make sure constantph.py and reference_energy.py are on sys.path or pass 'constantPHScript' in params.")
        raise

    generateLigandFF(params["ligandFile"], params["ligandFF"])

    runConstantPhSimulation(params)