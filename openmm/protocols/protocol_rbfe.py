# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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
Relative binding free energy (RBFE) for a set of congeneric ligands in a shared receptor,
wrapping Open Free Energy's `RelativeHybridTopologyProtocol` - the hybrid-topology, HREX,
OpenMM-engine protocol benchmarked on >1700 ligands in Baumann et al., J. Chem. Inf. Model.
2026, 66, 6429-6452.

Follows openfe's own documented three-stage workflow (see
claude/decisions/openmm/OpenFE_FEP.md):

  1. prepareInputsStep - receptor PDB (pdbfixer) + one multi-molecule ligand SDF, with
                     hydrogens added.
  2. planStep      - openfeSetupRBFE.py: assign partial charges once per ligand, build a
                     LigandNetwork with an atom mapper (LOMAP or Kartograf) + LOMAP scorer +
                     minimal-spanning-network planner, then emit one transformation JSON per
                     (edge, leg) - two legs per edge, 'solvent' and 'complex'.
  3. runAllTransformationsStep - `openfe quickrun` per transformation JSON.
  4. createOutputStep - ddG_bind for each edge from the thermodynamic cycle
                     ddG_bind(A->B) = dG_complex(A->B) - dG_solvent(A->B).

Everything not on the form stays at `RelativeHybridTopologyProtocol.default_settings()`, which
is the paper's unmodified default protocol - so the defaults here are the benchmark protocol.

This module is deliberately self-contained (it duplicates a little input-prep/result-parsing
code with protocol_abfe.py) so that one protocol is one file.
"""

import os
import json
import math

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import convertToSdf, mergeSDFs, getBaseName, addHydrogensToMol

from .. import Plugin
from ..constants import (OPENMM_DIC, DEFAULT_TEMPERATURE, DEFAULT_SOLVENT_PADDING,
                          DEFAULT_PROTOCOL_REPEATS, DEFAULT_EQUIL_LENGTH,
                          DEFAULT_RBFE_PRODUCTION, DEFAULT_RBFE_N_REPLICAS,
                          DEFAULT_MINIMIZATION_STEPS, MIN_SOLVENT_PADDING)

SMALL_MOL_FFS = ['openff-2.2.0', 'openff-2.1.0', 'openff-2.0.0', 'gaff-2.11', 'espaloma-0.3.2']
CHARGE_METHODS = ['am1bcc', 'am1bccelf10', 'nagl', 'espaloma']
MAPPERS = ['LOMAP', 'Kartograf']
COMPLEX, SOLVENT = 'complex', 'solvent'


class ProtOpenFERBFE(EMProtocol):
    """
    Relative binding free energy between congeneric ligands, via Open Free Energy's
    hybrid-topology RBFE protocol (LOMAP/Kartograf atom mapping, alchemical network planning,
    Hamiltonian replica exchange, MBAR) on the OpenMM GPU engine.

    Give it a docked set of ligands sharing one receptor: it plans a minimal-spanning
    alchemical network over them and computes ddG_bind for every edge. Defaults reproduce the
    protocol benchmarked in Baumann et al. 2026 (openff-2.2.0/Sage, Amber14SB, TIP3P, 298.15 K,
    11 lambda windows, 5 ns/window, 3 repeats).
    """
    _label = 'RBFE (OpenFE relative binding free energy)'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label='Use GPU for execution: ',
                       help='Sets openfe "engine_settings.compute_platform" to CUDA. Alchemical free '
                            'energy is impractical without a GPU.')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label='Choose GPU IDs')
        form.addParallelSection(threads=4, mpi=1)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False, important=True,
                      help='Set of docked/aligned ligands sharing a common receptor. The receptor '
                           'is taken from the set itself, so the ligands must already be posed in '
                           'its binding site (e.g. the output of a docking protocol). At least 2 '
                           'ligands are needed - use the ABFE protocol for a single ligand.')

        form.addSection(label='RBFE settings')
        nGroup = form.addGroup('Alchemical network')
        nGroup.addParam('atomMapper', params.EnumParam, default=0, choices=MAPPERS,
                        display=params.EnumParam.DISPLAY_HLIST, label='Atom mapper: ',
                        help='LOMAP: maximum-common-substructure mapping (the paper\'s default). '
                             'Kartograf: 3D-geometry-based mapping - can map atoms LOMAP will not, '
                             'but is sensitive to how well the input poses are aligned (a documented '
                             'source of large errors in the paper when poses are poorly aligned).')
        nGroup.addParam('max3d', params.FloatParam, default=1.0, expertLevel=params.LEVEL_ADVANCED,
                        label='Max 3D distance (A): ',
                        help='Mapped atoms further apart than this in the input poses are not mapped '
                             '(LomapAtomMapper max3d / Kartograf atom_max_distance).')
        nGroup.addParam('elementChange', params.BooleanParam, default=False,
                        expertLevel=params.LEVEL_ADVANCED, label='Allow element changes: ',
                        help='Whether the mapper may map atoms of different elements onto each other.')

        sGroup = form.addGroup('Sampling')
        sGroup.addParam('nReplicas', params.IntParam, default=DEFAULT_RBFE_N_REPLICAS,
                        label='Lambda windows: ',
                        help='openfe "lambda_settings.lambda_windows" (and the matching '
                             '"simulation_settings.n_replicas"), sampled with Hamiltonian replica '
                             'exchange. 11 is the paper\'s value for neutral transformations; it used '
                             '22 for charge-changing ones (and 20 ns/window instead of 5). This is a '
                             'CORRECTNESS parameter, not a cost knob - see this protocol\'s warnings '
                             'before lowering it.')
        sGroup.addParam('productionLength', params.FloatParam, default=DEFAULT_RBFE_PRODUCTION,
                        label='Production per window (ns): ',
                        help='openfe "simulation_settings.production_length". 5 ns/window is the '
                             'paper\'s default for neutral transformations; it reports results are '
                             'already ~converged by 80% of that. This IS the cost knob - lower it '
                             'rather than the window count.')
        sGroup.addParam('minimizationSteps', params.IntParam, default=DEFAULT_MINIMIZATION_STEPS,
                        expertLevel=params.LEVEL_ADVANCED, label='Minimization steps per window: ',
                        help='openfe "simulation_settings.minimization_steps". Applied to EVERY '
                             'lambda window, so at short sampling times it can dominate the run: '
                             '5000 steps x 11 windows is more integration work than 10 ps of MD per '
                             'window. BUT: Lowering this is NOT a safe way to shorten a run: measured at 100 steps, a freshly-solvated system still has clashes and the first MD blows up with "OpenMMException: Particle coordinate is NaN". OpenMM also stops early once converged, so 5000 is an upper bound rather than a fixed cost. Cut sampling time instead.')

        gGroup = form.addGroup('Force field and thermodynamics')
        gGroup.addParam('smallMolFF', params.EnumParam, default=0, choices=SMALL_MOL_FFS,
                        label='Small molecule force field: ',
                        help='openfe "small_molecule_forcefield". Default openff-2.2.0 is the Open '
                             'Force Field Sage 2.2.0 used throughout the benchmark paper. Protein and '
                             'water force fields are left at openfe defaults (Amber14SB + TIP3P).')
        gGroup.addParam('chargeMethod', params.EnumParam, default=0, choices=CHARGE_METHODS,
                        label='Ligand partial charges: ',
                        help='Charges are assigned ONCE per unique ligand and reused across all '
                             'transformations and repeats (openfe bulk_assign_partial_charges), which '
                             'is what the paper does to avoid conformer-dependent charge '
                             'irreproducibility. am1bcc (AmberTools/Antechamber) is the paper\'s choice.')
        gGroup.addParam('temperature', params.FloatParam, default=DEFAULT_TEMPERATURE,
                        label='Temperature (K): ',
                        help='openfe "thermo_settings.temperature". 298.15 K in the paper.')
        gGroup.addParam('solventPadding', params.FloatParam, default=DEFAULT_SOLVENT_PADDING,
                        label='Solvent padding (nm): ',
                        help='openfe "solvation_settings.solvent_padding" - distance from the solute to '
                             'the edge of the DODECAHEDRAL water box. Do not lower below ~1.3 nm: that '
                             'box\'s c.z component is only 0.7071 of its length, so a smaller pad makes '
                             'the half-box shorter than OpenMM\'s 0.9 nm nonbonded cutoff and the '
                             'solvent leg dies outright. The paper quotes 1.2 nm, which is not safe '
                             'with openfe\'s dodecahedral-box + 0.9 nm-cutoff defaults.')
        gGroup.addParam('protocolRepeats', params.IntParam, default=DEFAULT_PROTOCOL_REPEATS,
                        label='Independent repeats: ',
                        help='openfe "protocol_repeats" - independent replicas of each transformation, '
                             'differing only in initial velocities, averaged into the final estimate. '
                             '3 is both openfe\'s default and the paper\'s protocol. With 1 repeat '
                             'openfe reports an uncertainty of exactly 0, meaning "no estimate".')
        gGroup.addParam('equilLength', params.FloatParam, default=DEFAULT_EQUIL_LENGTH,
                        expertLevel=params.LEVEL_ADVANCED, label='Equilibration per window (ns): ',
                        help='openfe "simulation_settings.equilibration_length". 1 ns in the paper.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        prepStep = self._insertFunctionStep(self.prepareInputsStep)
        planStep = self._insertFunctionStep(self.planStep, prerequisites=[prepStep])
        # The transformation list only exists once planStep has run, so the per-transformation
        # quickruns are fanned out from a single dispatching step rather than inserted here.
        runStep = self._insertFunctionStep(self.runAllTransformationsStep, prerequisites=[planStep])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[runStep])

    # -- input preparation --------------------------------------------------------
    def getReceptorFile(self):
        return os.path.abspath(self.inputSetOfMols.get().getProteinFile())

    def getSystemName(self):
        return getBaseName(self.getReceptorFile())

    def getReceptorPDB(self):
        """pdbfixer-cleaned receptor, the same way ProtOpenMMSystemPrep prepares one.
        openfe's ProteinComponent.from_pdb_file needs a complete, sane PDB."""
        recPDB = os.path.abspath(self._getExtraPath(f'{self.getSystemName()}_receptor.pdb'))
        if not os.path.exists(recPDB):
            args = f'{self.getReceptorFile()} --output {recPDB}'
            pwchemPlugin.runOPENBABEL(self, 'pdbfixer', args=args, cwd=self._getExtraPath())
        return recPDB

    def getLigandsSDF(self):
        return os.path.abspath(self._getExtraPath('ligands.sdf'))

    def prepareInputsStep(self):
        """Write the receptor PDB and a single multi-molecule SDF of the input ligands - the two
        inputs the setup script consumes (Chem.SDMolSupplier + ProteinComponent.from_pdb_file).

        Hydrogens are (re)added to every ligand, which is NOT optional for openfe: the OpenFF
        toolkit builds a full valence model and rejects anything it reads as a radical. A ligand
        extracted straight from a crystallographic PDB has no hydrogens at all, so every carbon
        comes back bare and parametrisation dies with "RadicalsNotSupportedError: ... Found 2
        radical electrons on molecule [C][C]/C([C])=[C]/...". Ligands that already carry explicit
        hydrogens are unaffected - rdkit only fills in missing ones."""
        self.getReceptorPDB()

        sdfFiles = []
        for mol in self.inputSetOfMols.get():
            mol = mol.clone()
            molFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
            sdfFile = os.path.abspath(convertToSdf(self, os.path.abspath(molFile)))
            # Not convertToSdf(addHydrogens=True): that flag is skipped entirely when the input is
            # already an .sdf (the function returns before reaching it), which would silently leave
            # an unprotonated ligand for exactly the inputs most likely to need it.
            addHydrogensToMol(self, os.path.dirname(sdfFile), sdfFile)
            sdfFiles.append(sdfFile)
        mergeSDFs(sdfFiles, self.getLigandsSDF())

    # -- planning -----------------------------------------------------------------
    def getTransformDir(self):
        return os.path.abspath(self._getExtraPath('transformations'))

    def getEdgesFile(self):
        return os.path.abspath(self._getExtraPath('edges.txt'))

    def planStep(self):
        """Build the ligand network and dump one transformation JSON per (edge, leg)."""
        os.makedirs(self.getTransformDir(), exist_ok=True)
        paramsFile = os.path.abspath(self._getExtraPath('rbfeSetup.txt'))
        with open(paramsFile, 'w') as f:
            f.write(f'ligandsSdf :: {self.getLigandsSDF()}\n')
            f.write(f'proteinPdb :: {self.getReceptorPDB()}\n')
            f.write(f'smallMolFF :: {self.getEnumText("smallMolFF")}\n')
            f.write(f'chargeMethod :: {self.getEnumText("chargeMethod")}\n')
            f.write(f'temperature :: {self.temperature.get()}\n')
            f.write(f'solventPadding :: {self.solventPadding.get()}\n')
            f.write(f'protocolRepeats :: {self.protocolRepeats.get()}\n')
            f.write(f'equilLength :: {self.equilLength.get()}\n')
            f.write(f'minimizationSteps :: {self.minimizationSteps.get()}\n')
            f.write(f'transformDir :: {self.getTransformDir()}\n')
            f.write(f'mapper :: {self.getEnumText("atomMapper")}\n')
            f.write(f'max3d :: {self.max3d.get()}\n')
            f.write(f'elementChange :: {self.elementChange.get()}\n')
            f.write(f'nReplicas :: {self.nReplicas.get()}\n')
            f.write(f'productionLength :: {self.productionLength.get()}\n')
            f.write(f'networkFile :: {os.path.abspath(self._getExtraPath("ligand_network.graphml"))}\n')
            f.write(f'edgesFile :: {self.getEdgesFile()}\n')
            if getattr(self, params.USE_GPU).get():
                f.write('computePlatform :: cuda\n')  # openfe's own default spelling is lowercase
                # Without this openfe leaves engine_settings.gpu_device_index at None and OpenMM
                # picks device 0, silently ignoring the GPU chosen on the form.
                gpuIds = getattr(self, params.GPU_LIST).get().replace(',', ' ').split()
                if gpuIds:
                    f.write(f'gpuIndex :: {" ".join(gpuIds)}\n')

        Plugin.runScript(self, 'openfeSetupRBFE.py', args=paramsFile, env=OPENMM_DIC,
                         cwd=self._getExtraPath())

    # -- execution ----------------------------------------------------------------
    def getResultsDir(self):
        return os.path.abspath(self._getExtraPath('results'))

    def getResultsFile(self):
        return self._getExtraPath('results.txt')

    def getTransformationNames(self):
        """Names of the transformation JSONs the setup step produced. Read from disk (not
        precomputed in _insertAllSteps) because the network planner decides how many edges
        there are."""
        tDir = self.getTransformDir()
        if not os.path.exists(tDir):
            return []
        return sorted(getBaseName(f) for f in os.listdir(tDir) if f.endswith('.json'))

    def quickrunStep(self, name):
        """One `openfe quickrun` per transformation - openfe's own recommended unit of work."""
        os.makedirs(self.getResultsDir(), exist_ok=True)
        workDir = os.path.join(self.getResultsDir(), name)
        os.makedirs(workDir, exist_ok=True)

        transFile = os.path.join(self.getTransformDir(), f'{name}.json')
        outFile = os.path.join(self.getResultsDir(), f'{name}.json')
        args = f'{transFile} -o {outFile} -d {workDir}'
        # runOpenMM, not a dedicated openfe runner: openfe is installed into the same conda env as
        # OpenMM (see Plugin.addOPENMMPackage), and runOpenMM is exactly "activate that env, then
        # run this program" - the same helper ProtOpenDuckSimulation uses for the OpenDuck CLI.
        Plugin.runOpenMM(self, 'openfe quickrun', args, cwd=self.getResultsDir())

    def runAllTransformationsStep(self):
        """Run `openfe quickrun` on every planned transformation. Sequential within this step:
        each individual quickrun already saturates a GPU, so running several at once on one
        device would only contend for it."""
        names = self.getTransformationNames()
        if not names:
            raise RuntimeError('The planning step produced no transformations - check that the '
                               'input ligands are congeneric enough for the atom mapper to connect '
                               'them (see the ligand_network.graphml / setup log in extra/).')
        for i, name in enumerate(names, start=1):
            self.info(f'openfe quickrun {i}/{len(names)}: {name}')
            self.quickrunStep(name)

    # -- results ------------------------------------------------------------------
    def parseResult(self, name):
        """(estimate, uncertainty) in kcal/mol from a quickrun result JSON, or (None, None).

        Never raises - a single failed/missing transformation must not take down the whole
        protocol's output step. openfe serialises quantities as {'magnitude': x, 'unit': ...} in
        some versions and as a plain number in others, so both are handled."""
        resFile = os.path.join(self.getResultsDir(), f'{name}.json')
        if not os.path.exists(resFile):
            return None, None
        try:
            with open(resFile) as f:
                data = json.load(f)
            return self._asFloat(data.get('estimate')), self._asFloat(data.get('uncertainty'))
        except Exception as e:
            self.warning(f'Could not parse the openfe result for "{name}": {e}')
            return None, None

    @staticmethod
    def _asFloat(value):
        if value is None:
            return None
        if isinstance(value, dict):
            value = value.get('magnitude')
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def parseEdges(self):
        """(edgeName, ligandA, ligandB) per planned edge, written by the setup script so the leg
        JSONs can be paired back up without re-deriving the network here."""
        edges = []
        if not os.path.exists(self.getEdgesFile()):
            return edges
        with open(self.getEdgesFile()) as f:
            for line in f:
                parts = [p.strip() for p in line.strip().split('::')]
                if len(parts) == 3:
                    edges.append(tuple(parts))
        return edges

    def uncertaintyNote(self):
        """openfe derives a transformation's uncertainty from the SPREAD ACROSS REPEATS, so with a
        single repeat it reports exactly 0.0 - which in a results table reads like a perfectly
        precise number when it actually means "no uncertainty estimate exists"."""
        if self.protocolRepeats.get() > 1:
            return ''
        return ('# NOTE: run with a single repeat, so every "+/- 0.00" below means NO uncertainty\n'
                '# estimate was available, not a precise result - openfe derives it from the\n'
                '# spread across repeats. Use 3 repeats (the default) for a real error bar.\n')

    def writeResultsFile(self):
        """Combine both legs of every edge and write extra/results.txt. Returns
        (edgeResults, nPlannedEdges). Kept separate from createOutputStep so the free-energy
        bookkeeping can be tested without the project/mapper machinery that registering an
        output object requires."""
        lines, edgeResults = [], {}
        edges = self.parseEdges()
        for edgeName, ligA, ligB in edges:
            dgC, errC = self.parseResult(f'{edgeName}_{COMPLEX}')
            dgS, errS = self.parseResult(f'{edgeName}_{SOLVENT}')
            if dgC is None or dgS is None:
                lines.append(f'{ligA} -> {ligB}: incomplete (complex={dgC}, solvent={dgS})')
                continue
            ddg = dgC - dgS
            err = math.sqrt((errC or 0.0) ** 2 + (errS or 0.0) ** 2)
            edgeResults[(ligA, ligB)] = (ddg, err)
            lines.append(f'{ligA} -> {ligB}: ddG_bind = {ddg:+.2f} +/- {err:.2f} kcal/mol '
                        f'(complex {dgC:+.2f}, solvent {dgS:+.2f})')

        header = ('# Relative binding free energies (OpenFE hybrid topology RBFE)\n'
                  '# ddG_bind(A->B) = dG_complex(A->B) - dG_solvent(A->B), in kcal/mol\n'
                  '# A negative value means ligand B binds more strongly than ligand A.\n')
        with open(self.getResultsFile(), 'w') as f:
            f.write(header + self.uncertaintyNote() + '\n'.join(lines) + '\n')
        return edgeResults, len(edges)

    def createOutputStep(self):
        from ..objects import OpenMMSystem
        edgeResults, nEdges = self.writeResultsFile()

        outSystem = OpenMMSystem(filename=self.getReceptorPDB())
        outSystem.setFreeEnergyFile(self.getResultsFile())
        # Only meaningful as a single number when the network genuinely HAS one edge; with a whole
        # network the per-edge table in results.txt is the real result. Keyed on the planned edge
        # count, not on how many succeeded - otherwise a 5-edge network in which 4 edges failed
        # would publish the surviving one as if it were the whole answer.
        if nEdges == 1 and edgeResults:
            outSystem.setFreeEnergy(list(edgeResults.values())[0][0])
        self._defineOutputs(outputSystem=outSystem)
        self._defineSourceRelation(self.inputSetOfMols, outSystem)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        mols = self.inputSetOfMols.get()
        if mols is not None:
            if not mols.isDocked():
                errors.append('The input set of molecules must be docked: binding free energy needs '
                              'ligands posed in the receptor\'s binding site.')
            if len(mols) < 2:
                errors.append('RBFE needs at least 2 ligands to form an alchemical transformation. '
                              'Use the ABFE protocol for a single ligand.')
        if self.protocolRepeats.get() < 1:
            errors.append('Independent repeats must be at least 1.')
        if self.nReplicas.get() < 2:
            errors.append('At least 2 lambda windows are needed.')
        if self.solventPadding.get() < MIN_SOLVENT_PADDING:
            errors.append(
                f'Solvent padding must be at least {MIN_SOLVENT_PADDING} nm. openfe solvates into a '
                f'dodecahedral box whose c.z component is only 0.7071 of the box length, and OpenMM '
                f'requires the nonbonded cutoff (0.9 nm) to be no more than half of it - so for a '
                f'small ligand anything below ~1.27 nm makes the solvent leg die with '
                f'"NonbondedForce: The cutoff distance cannot be greater than half the periodic box '
                f'size". Measured: 1.2 nm padding gives a half-box of 0.861 nm.')
        return errors

    def _warnings(self):
        ws = []
        if not getattr(self, params.USE_GPU).get():
            ws.append('Running without a GPU: alchemical free energy calculations are impractically '
                     'slow on CPU (the paper reports ~5.8 h per neutral edge on an NVIDIA L40S).')
        if self.protocolRepeats.get() == 1:
            ws.append('With a single repeat there is no inter-repeat reproducibility estimate, and '
                     'openfe reports an uncertainty of exactly 0 - which means "no estimate", not '
                     '"exact". The paper\'s protocol uses 3 repeats.')
        if self.nReplicas.get() < 8:
            ws.append(f'{self.nReplicas.get()} lambda windows is very few for a hybrid-topology '
                     'transformation and is likely to be counterproductive: with too few windows '
                     'consecutive states barely share phase space, so MBAR cannot converge. '
                     'Measured with 2 windows on a real edge: off-diagonal overlap 0.00035 (the '
                     'benchmark paper treats <0.03 as unconverged), a meaningless estimate, AND '
                     'roughly 25 minutes of wasted CPU because every bootstrap resample runs MBAR '
                     'to its 10000-iteration limit. Cut "Production per window" instead - that '
                     'costs precision without destroying the alchemical path.')

        mols = self.inputSetOfMols.get()
        nLig = len(mols) if mols is not None else 0
        if nLig:
            # A minimal spanning network over N ligands has N-1 edges, each run as 2 legs.
            nEdges = max(nLig - 1, 1)
            nSims = 2 * nEdges * self.protocolRepeats.get()
            ns = nSims * self.nReplicas.get() * self.productionLength.get()
            ws.append(f'{nLig} ligands -> ~{nEdges} network edge(s) x 2 legs x '
                     f'{self.protocolRepeats.get()} repeat(s) = ~{nSims} replica-exchange '
                     f'simulations, ~{ns:.0f} ns of GPU MD in total. The paper reports ~5.8 h per '
                     'neutral edge on an NVIDIA L40S, so budget accordingly.')
        ws.append('Charge-changing edges need a larger schedule than the neutral default (the paper '
                 'used 22 windows and 20 ns/window for those); this protocol applies one schedule to '
                 'every edge, so raise the window count/production length if your set mixes charge '
                 'states.')
        return ws

    def _summary(self):
        if self.isFinished() and os.path.exists(self.getResultsFile()):
            with open(self.getResultsFile()) as f:
                return [f.read()]
        return ['The protocol has not finished.']

    def _methods(self):
        if not self.isFinished():
            return []
        return ['Relative binding free energies were computed with the Open Free Energy (OpenFE) '
                'hybrid-topology RBFE protocol: an alchemical network was planned over the input '
                f'ligands with the {self.getEnumText("atomMapper")} atom mapper and the LOMAP scorer, '
                'and each edge was run as solvent and complex legs of Hamiltonian replica-exchange '
                'alchemical sampling on the OpenMM engine, analysed with MBAR.']
