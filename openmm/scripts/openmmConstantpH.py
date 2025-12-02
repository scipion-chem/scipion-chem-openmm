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

# ---------------------------
# Helpers
# ---------------------------

def parse_ph_values(singlePH, onePH, manyPH):
    if singlePH:
        return [onePH]
    else:
        return [float(x.strip()) for x in manyPH.split(',')]


def parse_txt_config(filename):
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


# ---------------------------
# Integrator factory
# ---------------------------
def create_integrator(params, temperature):
    step_size = params.get('stepSize', 0.004) * picoseconds
    fric = params.get('fricCoef', 1.0) / picosecond
    col_freq = params.get('colFreq', 1.0) / picosecond
    err_tol = params.get('errTol', 0.001)

    integrator_name = params.get('integrator', 'Langevin')

    print(f"[integrator] Creating integrator '{integrator_name}' stepSize={step_size}, temp={temperature}")
    if integrator_name == 'Verlet':
        return VerletIntegrator(step_size)
    elif integrator_name == 'Langevin':
        return LangevinIntegrator(temperature, fric, step_size)
    elif integrator_name == 'LangevinMiddle':
        return LangevinMiddleIntegrator(temperature, fric, step_size)
    elif integrator_name == 'NoseHoover':
        return NoseHooverIntegrator(temperature, 1.0 / picosecond, step_size)
    elif integrator_name == 'Brownian':
        return BrownianIntegrator(temperature, step_size)
    elif integrator_name == 'VariableVerlet':
        return VariableVerletIntegrator(step_size, err_tol)
    elif integrator_name == 'VariableLangevin':
        return VariableLangevinIntegrator(temperature, fric, step_size, err_tol)
    else:
        raise ValueError(f"Unknown integrator: {integrator_name}")


def compute_ref(model_file, variants_dict, target_pKa, params,
                explicitFF, implicitFF, explicit_params, implicit_params,
                integrator, relaxation_integrator):
    """
    Compute reference energies for a model residue (ASP, GLU, etc.) for constant pH simulation.

    Returns a dict: {index: [energy_state0, energy_state1, ...]}
    """
    print(f"[compute_ref] Loading model PDB: {model_file}")
    pdb = PDBFile(model_file)

    print("[compute_ref] Initializing ConstantPH for reference computation...")
    cph = ConstantPH(
        pdb.topology,
        pdb.positions,
        [7.0],  # arbitrary pH for reference computation
        explicitFF,
        implicitFF,
        variants_dict,
        {index: [0.0] * len(states) for index, states in variants_dict.items()},
        250,  # short "relaxation" steps
        explicit_params,
        implicit_params,
        integrator,
        relaxation_integrator
    )

    # quick positions print (first 10 atoms)
    try:
        positions = cph.simulation.context.getState(getPositions=True).getPositions()
        print(f"[compute_ref] First positions (nm): {positions[:10]}")
    except Exception as e:
        print(f"[compute_ref] Warning: couldn't read positions: {e}")

    # Compute reference energies
    print(f"[compute_ref] Starting ReferenceEnergyFinder for target pKa={target_pKa}")
    finder = ReferenceEnergyFinder(cph, target_pKa, params.get('temperature', 300) * kelvin)

    total_iterations = params.get('ref_total_iterations', 20000)
    chunk = params.get('ref_chunk', 200)
    start_time = time.time()
    for start in range(0, total_iterations, chunk):
        print(f"[compute_ref] Running finder iterations {start}..{start+chunk-1}")
        finder.findReferenceEnergies(iterations=chunk, substeps=10)
        # small sanity check: are positions finite?
        try:
            pos = cph.simulation.context.getState(getPositions=True).getPositions()
            for p in pos[:10]:
                if math.isnan(p.x) or math.isnan(p.y) or math.isnan(p.z):
                    raise ValueError("NaN in positions during reference energy computation")
        except Exception as e:
            print(f"[compute_ref] ERROR while checking positions: {e}")
            raise

    elapsed = time.time() - start_time
    print(f"[compute_ref] Finished reference finder in {elapsed:.1f}s")

    # Extract reference energies
    ref_energies = {index: cph.titrations[index].referenceEnergies for index in variants_dict}
    print(f"[compute_ref] Reference energies keys: {list(ref_energies.keys())}")
    return ref_energies


# ---------------------------
# Main simulation function
# ---------------------------

def run_constant_ph_simulation(params):

    print("\n--- Loading system ---")
    print(f"[params] inputPdb: {params.get('inputPdb')}")
    pdb = PDBFile(params['inputPdb'])
    print("[run] PDB loaded. Topology residues:", sum(1 for _ in pdb.topology.residues()))

    print("[run] Creating force fields...")
    print(f"  explicitFF: {params.get('explicitFF')}  explicitSolvent: {params.get('explicitSolvent')}")
    explicitFF = ForceField(params['explicitFF'], params['explicitSolvent'])
    print("  explicit ForceField created.")
    print(f"  implicitFF: {params.get('implicitFF')}  implicitSolvent: {params.get('implicitSolvent')}")
    implicitFF = ForceField(params['implicitFF'], params['implicitSolvent'])
    print("  implicit ForceField created.")

    explicit_params = dict(
        nonbondedMethod=PME,
        nonbondedCutoff=params['explicitCutoff'] * nanometers,
        constraints=params['constraints'],
        hydrogenMass=params['hydrogenMass'] * amu
    )

    implicit_params = dict(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=params['implicitCutoff'] * nanometers,
        constraints=params['constraints']
    )

    print(f"[run] NB params explicit_cutoff={explicit_params['nonbondedCutoff']}, implicit_cutoff={implicit_params['nonbondedCutoff']}")

    temperature = params['temperature'] * kelvin
    print(f"[run] Simulation temperature: {temperature}")

    # ---------------------------
    # Create integrators
    # ---------------------------
    integrator = create_integrator(params, temperature)
    relaxation_integrator = create_integrator(params, temperature)  # Can have different params if needed
    print("[run] Integrators created.")

    # -----------------------------------
    # Reference energies
    # -----------------------------------

    print("\n--- Computing reference energies ---")
    ref_energies = {}
    variants_dict = {}

    # ASP
    if 'ASP' in params['residuesToTitrate']:
        print("[run] Computing ASP reference energy...")
        ref_energies['ASP'] = compute_ref(
            params['aspModel'],
            {1: ['ASP', 'ASH']},
            3.9,
            params,
            explicitFF, implicitFF,
            explicit_params, implicit_params,
            integrator, relaxation_integrator
        )[1]
        variants_dict['ASP'] = ['ASP', 'ASH']

    # GLU
    if 'GLU' in params['residuesToTitrate']:
        print("[run] Computing GLU reference energy...")
        ref_energies['GLU'] = compute_ref(
            params['gluModel'],
            {1: ['GLU', 'GLH']},
            4.2,
            params,
            explicitFF, implicitFF,
            explicit_params, implicit_params,
            integrator, relaxation_integrator
        )[1]
        variants_dict['GLU'] = ['GLU', 'GLH']

    # CYS
    if 'CYS' in params['residuesToTitrate']:
        print("[run] Computing CYS reference energy...")
        ref_energies['CYS'] = compute_ref(
            params['cysModel'],
            {1: ['CYS', 'CYX']},
            7.1,
            params,
            explicitFF, implicitFF,
            explicit_params, implicit_params,
            integrator, relaxation_integrator
        )[1]
        variants_dict['CYS'] = ['CYS', 'CYX']

    # HIS (3 states)
    if 'HIS' in params['residuesToTitrate']:
        print("[run] Computing HIS reference energies (HID/HIE)...")
        hid = compute_ref(
            params['hisModel'],
            {1: ['HIP', 'HID']},
            7.1,
            params,
            explicitFF, implicitFF,
            explicit_params, implicit_params,
            integrator, relaxation_integrator
        )[1]

        hie = compute_ref(
            params['hisModel'],
            {1: ['HIP', 'HIE']},
            6.5,
            params,
            explicitFF, implicitFF,
            explicit_params, implicit_params,
            integrator, relaxation_integrator
        )[1]

        ref_energies['HIS'] = [
            0.0 * kilojoules_per_mole,
            hid[1],
            hie[1]
        ]
        variants_dict['HIS'] = ['HIP', 'HID', 'HIE']

    # LYS
    if 'LYS' in params['residuesToTitrate']:
        print("[run] Computing LYS reference energy...")
        ref_energies['LYS'] = compute_ref(
            params['lysModel'],
            {1: ['LYS', 'LYN']},
            10.5,
            params,
            explicitFF, implicitFF,
            explicit_params, implicit_params,
            integrator, relaxation_integrator
        )[1]
        variants_dict['LYS'] = ['LYS', 'LYN']

    # -----------------------------------
    # Assign residues
    # -----------------------------------
    ph_values = parse_ph_values(
        params['singlePH'], params['onePH'], params['manyPH']
    )
    print(f"[run] pH values to run: {ph_values}")

    sim_variants = {}
    sim_ref_energies = {}

    for residue in pdb.topology.residues():
        if residue.name in variants_dict:
            sim_variants[residue.index] = variants_dict[residue.name]
            sim_ref_energies[residue.index] = ref_energies[residue.name]

    print("Titrated residues:")
    for k, v in sim_variants.items():
        print("  Residue", k, "->", v)

    # -----------------------------------
    # Build ConstantPH Simulation
    # -----------------------------------

    print("\n--- Creating simulation ---")
    cph = ConstantPH(
        pdb.topology, pdb.positions, ph_values,
        explicitFF, implicitFF,
        sim_variants, sim_ref_energies,
        params['relaxSteps'],
        explicit_params, implicit_params,
        integrator, relaxation_integrator
    )
    print("[run] ConstantPH object created.")

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
            tolerance=params['minimTol'] * kilojoules_per_mole,
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
    for prod_idx in range(params['prodSteps']):
        cph.simulation.step(params['stepProd'])
        cph.attemptMCStep(temperature)

        # check positions quickly (first atoms)
        try:
            state = cph.simulation.context.getState(getPositions=True)
            positions = state.getPositions()
            # check for NaNs in first 20 atoms
            for i, p in enumerate(positions[:20]):
                if math.isnan(p.x) or math.isnan(p.y) or math.isnan(p.z):
                    raise ValueError(f"NaN detected in position at prod step {prod_idx}, atom {i}")
        except Exception as e:
            print(f"[production] ERROR: {e}")
            # re-raise to preserve stack trace if you want script to stop
            raise

        if prod_idx % max(1, params.get('prodSteps') // 10) == 0:
            states = [sim_variants[i][cph.titrations[i].currentIndex] for i in sim_variants]
            print(f"[production] step {prod_idx+1}/{params['prodSteps']} pH: {cph.pH[cph.currentPHIndex]} states: {states}")

    print("[run] Production finished successfully.")


# ---------------------------
# Entry point
# ---------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True,
                        help="TXT configuration file")

    args = parser.parse_args()

    params = parse_txt_config(args.params)

    # optional: paths to constantph/reference_energy scripts (if you pass them in params)
    constantph_path = params.get('constantPHScript')
    reference_energy_path = params.get('referenceEnergyScript')

    if constantph_path:
        scripts_dir = os.path.dirname(constantph_path)
        if scripts_dir and scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)

    # Import dynamically if needed (these files must be on sys.path)
    try:
        from constantph import ConstantPH
        from reference_energy import ReferenceEnergyFinder
    except Exception as e:
        print(f"[startup] Could not import constantph/reference_energy: {e}")
        print("[startup] Make sure constantph.py and reference_energy.py are on sys.path or pass 'constantPHScript' in params.")
        raise

    run_constant_ph_simulation(params)