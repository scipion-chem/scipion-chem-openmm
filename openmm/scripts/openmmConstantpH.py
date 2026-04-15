#!/usr/bin/env python3
"""
Runs constant pH simulations
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
from openmm import XmlSerializer
from utils import parseParams, addMoleculesFF

# ---------------------------
# Helpers
# ---------------------------

def parsePhValues(singlePH, onePH, manyPH):
    if singlePH:
        return [float(onePH)]

    if isinstance(manyPH, (list, tuple)):
        return [float(v) for v in manyPH]

    values = []
    for item in str(manyPH).split(','):
        item = item.strip().strip('"').strip("'")
        if not item:
            continue
        try:
            values.append(float(item))
        except ValueError:
            raise ValueError(f"Invalid pH value '{item}' in '{manyPH}'")

    if not values:
        raise ValueError(f"No valid pH values parsed from '{manyPH}'")

    return values


def parseValue(v):
    valueLower = v.lower()

    if valueLower in ("true", "false"):
        return valueLower == "true"

    try:
        return float(v) if "." in v else int(v)
    except ValueError:
        return [x.strip() for x in v.split(",")] if "," in v else v


def parseTxtConfig(filename):
    params = {}

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            k, v = (x.strip() for x in line.split("=", 1))
            params[k] = parseValue(v)

    return params

# ---------------------------
# Integrator factory
# ---------------------------
def createIntegrator(params, temperature):
    stepSize = params.get('stepSize', 0.004) * picoseconds
    fric = params.get('fricCoef', 1.0) / picosecond
    errTol = params.get('errTol', 0.001)

    integratorName = params.get('integrator', 'Langevin')

    if integratorName == 'Verlet':
        integrator = VerletIntegrator(stepSize)
    elif integratorName == 'Langevin':
        integrator = LangevinIntegrator(temperature, fric, stepSize)
    elif integratorName == 'LangevinMiddle':
        integrator = LangevinMiddleIntegrator(temperature, fric, stepSize)
    elif integratorName == 'NoseHoover':
        integrator = NoseHooverIntegrator(temperature, 1.0 / picosecond, stepSize)
    elif integratorName == 'Brownian':
        integrator = BrownianIntegrator(temperature, stepSize)
    elif integratorName == 'VariableVerlet':
        integrator = VariableVerletIntegrator(stepSize, errTol)
    elif integratorName == 'VariableLangevin':
        integrator = VariableLangevinIntegrator(temperature, fric, stepSize, errTol)
    else:
        raise ValueError(f"Unknown integrator: {integratorName}")

    return integrator


def computeRef(modelFile, variantsDict, targetPKa, params,
                explicitFF, implicitFF, explicitParams, implicitParams,
                integrator, relaxationIntegrator):
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

    for _ in range(0, totalIterations, chunk):
        finder.findReferenceEnergies(iterations=chunk, substeps=10)
        try:
            pos = cph.simulation.context.getState(getPositions=True).getPositions()
            for p in pos[:10]:
                if math.isnan(p.x) or math.isnan(p.y) or math.isnan(p.z):
                    raise ValueError("NaN in positions during reference energy computation")
        except Exception as e:
            print(f"[computeRef] ERROR while checking positions: {e}")
            raise

    refenergies = {index: cph.titrations[index].referenceEnergies for index in variantsDict}
    return refenergies

def addLigand(modeller, ligFile):
    '''Add ligand topology and coordinates to the modeller'''
    molecule = Molecule.from_file(ligFile)
    ligTop = molecule.to_topology().to_openmm()
    for residue in ligTop.residues():
        residue.name = "LIG"
    positions = molecule.conformers[0]

    ligPos = positions.m_as("nanometer") * nanometers

    modeller.add(ligTop, ligPos)
    return modeller

# ---------------------------
# Main simulation function
# ---------------------------

def runConstantPhSimulation(params):
    pdb = PDBFile(params['inputPdb'])

    explicitFF = ForceField(params['explicitFF'], params['explicitSolvent'])
    print("  explicit ForceField created.")
    implicitFF = ForceField(params['implicitFF'], params['implicitSolvent'])
    print("  implicit ForceField created.")

    nonbondedMethods = {
        'NoCutoff': NoCutoff,
        'CutoffNonPeriodic': CutoffNonPeriodic,
        'CutoffPeriodic': CutoffPeriodic,
        'Ewald': Ewald,
        'PME': PME,
        'LJPME': LJPME
    }
    nonBondedMethodExp = nonbondedMethods[params['nonBondedMethodExp']]
    nonBondedMethodImp = nonbondedMethods[params['nonBondedMethodImp']]

    explicitParams = dict(
        nonbondedMethod=nonBondedMethodExp,
        nonbondedCutoff=params['explicitCutoff'] * nanometers,
        constraints=params['constraintsExp'],
        hydrogenMass=1.5 * amu
    )

    implicitParams = dict(
        nonbondedMethod=nonBondedMethodImp,
        nonbondedCutoff=params['implicitCutoff'] * nanometers,
        constraints=params['constraintsImp']
    )

    temperature = params['temperature'] * kelvin

    # ---------------------------
    # Create integrators
    # ---------------------------
    integrator = createIntegrator(params, temperature)
    relaxationIntegrator = createIntegrator(params, temperature)
    print("[run] Integrators created.")

    # -----------------------------------
    # Reference energies
    # -----------------------------------
    refenergies = {}
    variantsDict = {}

    modeller = Modeller(pdb.topology, pdb.positions)

    if params.get('ligandFile'):
        print("[run] Adding ligand coordinates...")
        modeller = addLigand(modeller, params['ligandFile'])

        print("[run] Adding ligand template to ForceField...")
        explicitFF = addMoleculesFF(explicitFF, params['ligandFile'], params['ligandFF'])
        implicitFF = addMoleculesFF(implicitFF, params['ligandFile'], params['ligandFF'])

    if params['addHydrogens']:
        print("[run] Adding hydrogens...")
        modeller.addHydrogens(explicitFF, pH=float(params['hPH']))


    boxSize = params.get('boxSize', None)
    if boxSize:
        bSize = list(map(float, params['boxSize'].split(',')))
        kwargs = {"boxSize": Vec3(bSize[0], bSize[1], bSize[2]) * nanometers}
    else:
        kwargs = {"padding": float(params['padding'])}

    kwargs.update({"ionicStrength": float(params['saltConc']) * molar, "neutralize": (params['neutralize']),
                   "positiveIon": params['cationType'], "negativeIon": params['anionType']})

    waterModel = params['explicitSolvent']
    wModel = os.path.splitext(os.path.basename(waterModel))[0]
    modeller.addSolvent(explicitFF, model=wModel, **kwargs)

    resInfo = {
        'ASP': {'model': 'aspModel', 'states': [(['ASP', 'ASH'], 3.9)]},
        'GLU': {'model': 'gluModel', 'states': [(['GLU', 'GLH'], 4.2)]},
        'CYS': {'model': 'cysModel', 'states': [(['CYS', 'CYX'], 7.1)]},
        'LYS': {'model': 'lysModel', 'states': [(['LYS', 'LYN'], 10.5)]},
        'HIS': {'model': 'hisModel', 'states': [(['HIP', 'HID'], 7.1), (['HIP', 'HIE'], 6.5)]},
    }

    for res, info in resInfo.items():
        if res not in params['residuesToTitrate']:
            continue

        print(f"[run] Computing {res} reference energy{'s' if len(info['states']) > 1 else ''}...")

        variantsDict[res] = []
        refenergies[res] = []

        if res != 'HIS':
            # Standard 2-state residues
            stateVariants, pKa = info['states'][0]
            ref = computeRef(
                params[info['model']],
                {1: stateVariants},
                pKa,
                params,
                explicitFF, implicitFF,
                explicitParams, implicitParams,
                integrator, relaxationIntegrator
            )[1]

            variantsDict[res] = stateVariants
            refenergies[res] = ref
        else:
            # HIS: 3 states
            hidVariants, hidPka = info['states'][0]
            hieVariants, hiePka = info['states'][1]

            hid = computeRef(
                params[info['model']],
                {1: hidVariants},
                hidPka,
                params,
                explicitFF, implicitFF,
                explicitParams, implicitParams,
                integrator, relaxationIntegrator
            )[1]

            hie = computeRef(
                params[info['model']],
                {1: hieVariants},
                hiePka,
                params,
                explicitFF, implicitFF,
                explicitParams, implicitParams,
                integrator, relaxationIntegrator
            )[1]

            # OpenMM expects referenceEnergies length == number of variants (3 for HIS)
            refenergies['HIS'] = [0.0 * kilojoules_per_mole, hid[1], hie[1]]
            variantsDict['HIS'] = ['HIP', 'HID', 'HIE']

    # -----------------------------------
    # Assign residues
    # -----------------------------------
    phValues = parsePhValues(
        params['singlePH'], params['onePH'], params['manyPH']
    )

    simVariants = {}
    simRefenergies = {}

    for residue in pdb.topology.residues():
        if residue.name == 'LIG':
            continue
        if residue.name in variantsDict:
            simVariants[residue.index] = variantsDict[residue.name]
            simRefenergies[residue.index] = refenergies[residue.name]

    # -----------------------------------
    # Build ConstantPH Simulation
    # -----------------------------------

    filteredTopology = modeller.topology
    filteredPositions = modeller.positions

    cph = ConstantPH(
        filteredTopology, filteredPositions, phValues,
        explicitFF, implicitFF,
        simVariants, simRefenergies,
        params['relaxSteps'],
        explicitParams, implicitParams,
        integrator, relaxationIntegrator
    )
    print("[run] ConstantPH object created.")

    systemXml = params['systemXml']
    with open(systemXml, "w") as f:
        f.write(XmlSerializer.serialize(cph.simulation.system))

    trajFile = params['trajFile']
    logFile = params['logFile']
    finalPdb = params['finalPdb']
    reportEvery = params['reportEvery']

    #cph.simulation.reporters.append(DCDReporter(trajFile, reportEvery))
    cph.simulation.reporters.append(
        DCDReporter(trajFile, reportEvery, enforcePeriodicBox=True)
    )
    totalSteps = params['nSteps']
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
            progress=False,
            remainingTime=False,
            speed=True,
            totalSteps=totalSteps,
            separator=','
        )
    )

    if params.get('addBarostat', False):
        cph.simulation.system.addForce(
            MonteCarloBarostat(params['pressure'] * bar, temperature)
        )
        cph.simulation.context.reinitialize(preserveState=True)
        print("[run] Barostat added and context reinitialized.")

    if params.get('addMinimization', True):
        cph.simulation.minimizeEnergy(
            tolerance=params['minimTol'] * kilojoules_per_mole / nanometer ,
            maxIterations=params['maxIter']
        )
        print("[run] Minimization finished.")

    nSteps = params['nSteps']
    equilSteps = params.get('equilSteps', nSteps // 10)
    prodSteps = nSteps - equilSteps
    stepEquil = params.get('stepEquil', 1)
    stepProd = params.get('stepProd', 1)
    # -----------------------------------
    # Equilibration
    # -----------------------------------

    for idx in range(equilSteps):
        cph.simulation.step(stepEquil)
        cph.attemptMCStep(temperature)

        if idx % max(1, equilSteps // 10) == 0:
            print(f"[equil] Completed equilibration cycle {idx + 1}/{equilSteps}")

    # -----------------------------------
    # Production
    # -----------------------------------
    for prodIdx in range(prodSteps):
        cph.simulation.step(stepProd)
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

        if prodIdx % max(1, prodSteps // 10) == 0:
            states = [simVariants[i][cph.titrations[i].currentIndex] for i in simVariants]
            print(f"[production] step {prodIdx + 1}/{prodSteps} pH: {cph.pH[cph.currentPHIndex]} states: {states}")

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

    runConstantPhSimulation(params)