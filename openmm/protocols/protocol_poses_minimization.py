# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This protocol minimizes a set of docked / cofolded poses together with the side chains of their binding
pocket using OpenMM (implicit solvent). It is meant to clean poses with steric clashes or bad ligand geometry
(bond lengths, flat saturated rings...) while keeping them close to the original pose.
"""
import os, json

from pyworkflow.protocol import params
from pyworkflow.utils import Message
import pyworkflow.object as pwobj

from pwchem.utils import getBaseName, convertToSdf

from .. import Plugin
from ..constants import OPENMM_DIC
from .protocol_system_prep import ProtOpenMMSystemPrep

IMPLICIT_FILES = ['implicit/obc2.xml', 'implicit/gbn2.xml', 'implicit/obc1.xml', 'implicit/hct.xml', 'None']
IMPLICIT_NAMES = ['OBC2', 'GBn2', 'OBC1', 'HCT', 'None (vacuum)']
CHARGES = ['Force field default', 'Gasteiger']
SCRIPT_NAME = 'openmmMinimizePose.py'


class ProtOpenMMPosesMinimization(ProtOpenMMSystemPrep):
    """
    Minimize the poses of a set of small molecules together with the side chains of the surrounding pocket.

    AI Generated:

        ProtOpenMMPosesMinimization - User Manual

        Overview
        --------
        For each pose of the input set, the receptor (protein only) and the ligand are parametrized
        (Amber protein force field + GAFF / SMIRNOFF / Espaloma for the ligand) and minimized in implicit
        solvent. Only the ligand and the side chains of the residues within a given radius of it can move;
        the rest of the receptor is frozen. The ligand heavy atoms are restrained to the input pose, and the
        restraint is released progressively (stages), so clashes are solved mostly by the side chains.

        The protocol only minimizes: no validation is done. Use a PoseBusters protocol over the output and, if
        needed, run this protocol again over the poses that still fail (changing the number of stages,
        iterations, pocket radius or adding a random perturbation of the ligand, which helps flat saturated rings).

        Inputs
        ------
        - **inputSetOfMols**: SetOfSmallMolecules with poses (and the receptor of each pose, or of the set).

        Outputs
        -------
        - **outputSmallMolecules**: a copy of the input set where the pose and receptor files are the minimized
          ones (per molecule), with the attributes minStatus (minimized / error), minRMSD (ligand heavy atoms vs
          the input pose), minEnergyInit and minEnergyFinal (kJ/mol). Molecules that failed keep the input files.
          The minimized complex is also written in extra/<molecule>/<molecule>_min_complex.pdb.
    """
    _label = 'poses minimization'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, label="Input poses: ", important=True,
                      pointerClass='SetOfSmallMolecules',
                      help='Set of small molecules with poses (and receptors) to minimize')

        form.addSection('Minimization')
        form.addParam('pocketRadius', params.FloatParam, default=8.0, label='Pocket radius (A): ',
                      help='Side chains of the residues with any heavy atom closer than this distance to the '
                           'ligand are free to move. The rest of the protein is frozen. In the retries the radius '
                           'is increased 2 A at a time.')
        form.addParam('stages', params.StringParam, default='10 2 0.5 0', label='Ligand restraint stages: ',
                      help='Force constants (kcal/mol/A^2) of the harmonic restraint on the ligand heavy atoms '
                           'with respect to the input pose. One minimization is done per value, in order. '
                           '0 means free ligand.')
        form.addParam('maxIter', params.IntParam, default=2000, expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum iterations per stage: ', help='Maximum minimizer iterations in each stage')
        form.addParam('tolerance', params.FloatParam, default=10.0, expertLevel=params.LEVEL_ADVANCED,
                      label='Tolerance (kJ/mol/nm): ', help='Minimizer convergence threshold on the forces')
        form.addParam('perturbSigma', params.FloatParam, default=0.0, label='Initial ligand perturbation (A): ',
                      help='Standard deviation of a random displacement applied to the ligand atoms before '
                           'minimizing (the restraint keeps the input pose as reference). 0 = none. A small value '
                           '(0.1) lets perfectly flat saturated rings pucker.')
        form.addParam('seed', params.IntParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                      condition='perturbSigma>0', label='Random seed: ', help='Seed of the random perturbation')
        form.addParam('implicitSolvent', params.EnumParam, default=0, choices=IMPLICIT_NAMES,
                      label='Implicit solvent: ', help='Generalized Born model used during the minimization')
        form.addParam('hPH', params.FloatParam, default=7.0, label='pH for hydrogens: ',
                      help='pH used to add the hydrogens of the receptor')

        form.addSection('Force fields')
        self._defineFFParams(form)
        form.addParam('ligCharges', params.EnumParam, default=0, choices=CHARGES, expertLevel=params.LEVEL_ADVANCED,
                      label='Ligand partial charges: ',
                      help='Force field default (e.g. AM1-BCC) or fast Gasteiger charges, useful if AM1-BCC '
                           '(AmberTools) is not available')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        convStep = self._insertFunctionStep(self.convertInputStep, prerequisites=[])
        minSteps = []
        for molName, _, _ in self.getMolsInfo():
            minSteps.append(self._insertFunctionStep(self.minimizeStep, molName, prerequisites=[convStep]))
        self._insertFunctionStep(self.createOutputStep, prerequisites=minSteps)

    def convertInputStep(self):
        for molName, molFile, _ in self.getMolsInfo():
            oDir = self.getMolDir(molName)
            os.makedirs(oDir, exist_ok=True)
            sdfFile = os.path.join(oDir, f'{molName}_input.sdf')
            convertToSdf(self, molFile, sdfFile=sdfFile, overWrite=True)

    def minimizeStep(self, molName):
        oDir = self.getMolDir(molName)
        recFile = [r for n, _, r in self.getMolsInfo() if n == molName][0]
        paramsFile = os.path.join(oDir, 'minimizationParams.txt')
        with open(paramsFile, 'w') as f:
            f.write(f'molName :: {molName}\n')
            f.write(f'receptorFile :: {recFile}\n')
            f.write(f'ligandFile :: {os.path.join(oDir, molName + "_input.sdf")}\n')
            f.write(f'mFF :: {self.getFFFiles()[0]}\n')
            f.write(f'implicitFF :: {IMPLICIT_FILES[self.implicitSolvent.get()]}\n')
            f.write(f'ligandFF :: {self.getLigandFFVersion()}\n')
            f.write(f'ligCharges :: {CHARGES[self.ligCharges.get()].split()[0].lower()}\n')
            f.write(f'hPH :: {self.hPH.get()}\n')
            f.write(f'pocketRadius :: {self.pocketRadius.get()}\n')
            f.write(f'stages :: {self.stages.get().strip()}\n')
            f.write(f'maxIter :: {self.maxIter.get()}\n')
            f.write(f'tolerance :: {self.tolerance.get()}\n')
            f.write(f'perturbSigma :: {self.perturbSigma.get()}\n')
            f.write(f'seed :: {self.seed.get() if self.perturbSigma.get() > 0 else 0}\n')
            f.write('cpuThreads :: 1\n')

        Plugin.runScript(self, SCRIPT_NAME, args=paramsFile, env=OPENMM_DIC, cwd=oDir)

    def createOutputStep(self):
        inMols = self.inputSetOfMols.get()
        outputSet = inMols.createCopy(self._getPath(), copyInfo=True)
        names = [n for n, _, _ in self.getMolsInfo()]

        for molName, mol in zip(names, inMols):
            res = self.getResult(molName)
            oDir = self.getMolDir(molName)
            outLig = os.path.join(oDir, f'{molName}_min.sdf')
            outRec = os.path.join(oDir, f'{molName}_min_receptor.pdb')

            if os.path.exists(outLig) and os.path.exists(outRec):
                mol.setPoseFile(outLig)
                mol.setProteinFile(outRec)
            mol.minStatus = pwobj.String(res.get('status', 'error'))
            mol.minRMSD = pwobj.Float(res.get('rmsd'))
            mol.minEnergyInit = pwobj.Float(res.get('e0'))
            mol.minEnergyFinal = pwobj.Float(res.get('e1'))
            outputSet.append(mol)

        self._defineOutputs(outputSmallMolecules=outputSet)
        self._defineSourceRelation(self.inputSetOfMols, outputSet)

    ####################### UTILS FUNCTIONS ############################
    def getMolsInfo(self):
        """[(unique molecule name, pose file, receptor file)] in set order"""
        info, used = [], set()
        setRec = self.inputSetOfMols.get().getProteinFile()
        for i, mol in enumerate(self.inputSetOfMols.get()):
            name = mol.getUniqueName().replace('/', '_').replace(' ', '_')
            if name in used:
                name = f'{name}_{i}'
            used.add(name)
            rec = mol.getProteinFile() or setRec
            info.append((name, os.path.abspath(mol.getPoseFile()), os.path.abspath(rec) if rec else None))
        return info

    def getMolDir(self, molName):
        return os.path.abspath(self._getExtraPath(molName))

    def getResult(self, molName):
        resFile = os.path.join(self.getMolDir(molName), 'result.json')
        if os.path.exists(resFile):
            with open(resFile) as f:
                return json.load(f)
        return {'status': 'error'}

    def _summary(self):
        s = []
        if self.isFinished() and hasattr(self, 'outputSmallMolecules'):
            rmsds, nErr = [], 0
            for molName, _, _ in self.getMolsInfo():
                res = self.getResult(molName)
                if res.get('status') == 'minimized':
                    rmsds.append(res['rmsd'])
                else:
                    nErr += 1
            s.append(f'Minimized poses: {len(rmsds)}, failed (input kept): {nErr}')
            if rmsds:
                s.append(f'Ligand RMSD to the input pose: mean {sum(rmsds) / len(rmsds):.2f} A, max {max(rmsds):.2f} A')
        return s

    def _warnings(self):
        ws = []
        if self.inputSetOfMols.get() and len(self.inputSetOfMols.get()) > 100:
            ws.append(f'{len(self.inputSetOfMols.get())} poses will be minimized (about 0.5-1.5 min each on CPU). '
                      f'Consider running it only over the poses that need it.')
        return ws

    def _validate(self):
        vs = []
        inSet = self.inputSetOfMols.get()
        if inSet:
            if not inSet.getProteinFile() and not all(m.getProteinFile() for m in inSet):
                vs.append('The input set must have a receptor file (set or molecule level)')
            if self.ffType.get() != 0:
                vs.append('Only the Amber force field family is supported for the minimization of poses')
            try:
                float(self.stages.get().replace(',', ' ').split()[0])
            except Exception:
                vs.append('The restraint stages must be a list of numbers, e.g.: 10 2 0.5 0')
        return vs
