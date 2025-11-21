#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
strip_water.py
Reads a params file with keys:
  pdbIn :: /path/to/pdb
  dcdIn :: /path/to/dcd   (may be empty)
  outPrefix :: /path/prefix
  keepIons :: True/False

Produces:
  {outPrefix}.pdb
  {outPrefix}.dcd   (only if dcdIn provided)
"""
import sys
import os
import mdtraj as md

def read_params(path):
    d = {}
    with open(path) as f:
        for line in f:
            if '::' in line:
                k,v = line.split('::',1)
                d[k.strip()] = v.strip()
    return d

def main(params_path):
    params = read_params(params_path)
    pdbIn = params.get('pdbIn','')
    dcdIn = params.get('dcdIn','').strip()
    outPrefix = params.get('outPrefix','systemNowater')
    keepIons = params.get('keepIons','True').lower() in ('1','true','yes')

    if not os.path.exists(pdbIn):
        print("ERROR: pdbIn not found:", pdbIn)
        sys.exit(1)

    # load pdb
    topo = md.load_pdb(pdbIn)

    # build selection string
    if keepIons:
        sel = "not water"
    else:
        # remove water and ions (ions identified as elements or common resnames)
        # keep only protein, nucleic acids, ligands, cofactors
        sel = "not water and not (resname NA CL K MG CA)".lower()
        # mdtraj expects lower-case selections, but 'resname' is matched case-insensitively in mdtraj
        # fallback: select by "not solvent" if mdtraj supports:
        # sel = "not water and not solvent"  # but "solvent" is not always supported

    # do selection on PDB topology
    notWater = topo.topology.select("not water") if keepIons else topo.topology.select("not water and not element Na and not element Cl and not element K and not element Mg and not element Ca")

    # save new pdb (first frame)
    pdbOut = f"{outPrefix}.pdb"
    topo.atom_slice(notWater)[0].save_pdb(pdbOut)
    print("Saved:", pdbOut)

    # if DCD provided, process it
    if dcdIn and os.path.exists(dcdIn):
        print("Loading DCD:", dcdIn)
        traj = md.load_dcd(dcdIn, top=pdbIn)  # use original pdb top to map indices
        print("Slicing trajectory (removing water)...")
        trajNow = traj.atom_slice(notWater)
        dcdOut = f"{outPrefix}.dcd"
        trajNow.save_dcd(dcdOut)
        print("Saved:", dcdOut)
    else:
        print("No DCD input provided or DCD missing; only generated PDB.")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: strip_water.py params.txt")
        sys.exit(1)
    main(sys.argv[1])