#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
strip_water.py
Reads a params file with keys:
  cifIn :: /path/to/pdb
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
import MDAnalysis as mda
from openmm.app import PDBxFile

def readParams(path):
    d = {}
    with open(path) as f:
        for line in f:
            if '::' in line:
                k,v = line.split('::',1)
                d[k.strip()] = v.strip()
    return d

def main(paramsPath):
    params = readParams(paramsPath)
    cifIn = params.get('cifIn','')
    dcdIn = params.get('dcdIn','').strip()
    outPrefix = params.get('outPrefix','systemNowater')
    keepIons = params.get('keepIons','True').lower() in ('1','true','yes')

    if not os.path.exists(cifIn):
        print("ERROR: cifIn not found:", cifIn)
        sys.exit(1)

    # load pdb
    cifObj = PDBxFile(cifIn)
    u = mda.Universe(cifObj.topology, dcdIn, topology_format='OPENMMTOPOLOGY')
    solvent = {'HOH', 'WAT', 'TIP3', 'TIP4', 'TIP5', 'SOL', 'SPC', 'SPCE'}
    if not keepIons:
        solvent = solvent.union({'NA', 'CL', 'K', 'MG', 'CA'})
    not_solvent = u.select_atoms(f"not resname {' '.join(solvent)}")
    print(f"Kept {not_solvent.n_atoms} atoms (removed ~{(u.atoms.n_atoms - not_solvent.n_atoms)} solvent atoms)")

    with mda.Writer(f"{outPrefix}.dcd", not_solvent.n_atoms) as W:
        for ts in u.trajectory:
            W.write(not_solvent)

    not_solvent.write(f"{outPrefix}.pdb")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: strip_water.py params.txt")
        sys.exit(1)
    main(sys.argv[1])