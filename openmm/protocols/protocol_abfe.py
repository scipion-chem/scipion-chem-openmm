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
Absolute binding free energy (ABFE) for ONE ligand in a receptor, wrapping Open Free Energy's
`AbsoluteBindingProtocol`. The ligand is picked from the input set with a wizard - one ABFE
protocol per ligand, since each is a full 44-window double-decoupling calculation.

openfe handles the whole double-decoupling cycle internally: state A is the ligand bound in the
solvated complex, state B is the same system with the ligand removed, and the Protocol itself
builds both the complex and solvent legs from those two states. It applies automatic Boresch
orientational restraints (1 bond + 2 angles + 3 dihedrals over 3 protein and 3 ligand atoms,
anchors picked from RMSF/secondary-structure/distance heuristics), fully annihilates ligand
electrostatics, decouples Lennard-Jones, samples with Hamiltonian replica exchange on the OpenMM
engine and estimates free energies with MBAR - so this protocol contributes only input
preparation, settings, execution and result collection.

  1. prepareInputsStep - receptor PDB (pdbfixer) + an SDF holding the selected ligand, with
                     hydrogens added.
  2. setupStep     - openfeSetupABFE.py: charge the ligand and build its Transformation JSON.
  3. runAllTransformationsStep - `openfe quickrun` on that transformation.
  4. createOutputStep - read the result's estimate/uncertainty (dG_bind, kcal/mol).

Everything not on the form stays at `AbsoluteBindingProtocol.default_settings()` (3 repeats,
30 complex / 14 solvent lambda windows, 10 ns production per replica).

Note ABFE is substantially more expensive than RBFE and, unlike the RBFE protocol, was NOT part
of the Baumann et al. 2026 benchmark - so its accuracy is not backed by that paper's statistics.
See claude/decisions/openmm/OpenFE_FEP.md SS3.1.

This module is deliberately self-contained (it duplicates a little input-prep/result-parsing
code with protocol_rbfe.py) so that one protocol is one file.
"""

import os
import json

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import convertToSdf, mergeSDFs, getBaseName, addHydrogensToMol

from .. import Plugin
from ..constants import (OPENMM_DIC, DEFAULT_TEMPERATURE, DEFAULT_SOLVENT_PADDING,
                          DEFAULT_PROTOCOL_REPEATS, DEFAULT_EQUIL_LENGTH,
                          DEFAULT_MINIMIZATION_STEPS, MIN_SOLVENT_PADDING)

SMALL_MOL_FFS = ['openff-2.2.0', 'openff-2.1.0', 'openff-2.0.0', 'gaff-2.11', 'espaloma-0.3.2']
CHARGE_METHODS = ['am1bcc', 'am1bccelf10', 'nagl', 'espaloma']


class ProtOpenFEABFE(EMProtocol):
    """
    Absolute binding free energy via Open Free Energy's AbsoluteBindingProtocol: automatic
    Boresch restraints, electrostatics annihilation + Lennard-Jones decoupling, Hamiltonian
    replica exchange on OpenMM, MBAR analysis.

    Unlike the RBFE protocol this needs no congeneric series and no atom mapping - a ligand is
    evaluated on its own, so structurally unrelated ligands can be compared by running one ABFE
    protocol each. It is however considerably more expensive per ligand.
    """
    _label = 'ABFE (OpenFE absolute binding free energy)'
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
                           'is taken from the set itself, so the ligand must already be posed in '
                           'its binding site (e.g. the output of a docking protocol).')
        form.addParam('inputLigand', params.StringParam, label='Ligand: ',
                      help='The single ligand to compute the absolute binding free energy for, '
                           'picked from the set above with the wizard. ABFE is run for ONE ligand '
                           'at a time: each one is a full double-decoupling calculation over 44 '
                           'lambda windows plus per-leg pre-equilibration, so running a whole set '
                           'in a single protocol would be many GPU-hours with no way to inspect or '
                           'restart an individual ligand. Add one ABFE protocol per ligand instead.')

        form.addSection(label='ABFE settings')
        rGroup = form.addGroup('Boresch restraints')
        rGroup.addParam('hostMinDistance', params.FloatParam, default=0.5,
                        label='Host anchor min distance (nm): ',
                        help='openfe "restraint_settings.host_min_distance" - the closest a receptor '
                             'atom may be to the ligand to be considered as a Boresch restraint '
                             'anchor. openfe picks the actual anchor atoms itself within this shell, '
                             'using RMSF, secondary structure and distance heuristics.')
        rGroup.addParam('hostMaxDistance', params.FloatParam, default=1.5,
                        label='Host anchor max distance (nm): ',
                        help='openfe "restraint_settings.host_max_distance" - the furthest a receptor '
                             'atom may be to be considered as a restraint anchor.')

        sGroup = form.addGroup('Sampling')
        sGroup.addParam('productionLength', params.FloatParam, default=10.0,
                        label='Production per window (ns): ',
                        help='openfe "*_simulation_settings.production_length" for both legs. '
                             'openfe\'s ABFE default is 10 ns per replica, across 30 complex + 14 '
                             'solvent lambda windows. ABFE convergence is dominated by the '
                             'near-fully-decoupled windows, so reduce with care.')
        sGroup.addParam('minimizationSteps', params.IntParam, default=DEFAULT_MINIMIZATION_STEPS,
                        expertLevel=params.LEVEL_ADVANCED, label='Minimization steps per window: ',
                        help='openfe "*_simulation_settings.minimization_steps". Applied to EVERY '
                             'one of the 44 lambda windows, so at short sampling times it dominates '
                             'the run. BUT: Lowering this is NOT a safe way to shorten a run: measured at 100 steps, a freshly-solvated system still has clashes and the first MD blows up with "OpenMMException: Particle coordinate is NaN". OpenMM also stops early once converged, so 5000 is an upper bound rather than a fixed cost. Cut sampling time instead.')
        sGroup.addParam('preEquilLength', params.FloatParam, default=0.0,
                        expertLevel=params.LEVEL_ADVANCED,
                        label='Pre-equilibration per leg (ns, 0 = openfe default): ',
                        help='Before the alchemical windows, openfe runs a plain MD pre-equilibration '
                             'of each leg (NVT + NPT equilibration + production). Its defaults are '
                             'large and asymmetric - 0.25+0.5+5.0 ns for the complex leg and '
                             '0.1+0.2+0.5 ns for the solvent leg, i.e. 6.55 ns per repeat before any '
                             'free energy sampling starts, which is hours on a workstation GPU. '
                             'Leave at 0 to keep openfe\'s own values; set a value to use it for all '
                             'three phases of BOTH legs (mainly useful for quick smoke runs).')

        gGroup = form.addGroup('Force field and thermodynamics')
        gGroup.addParam('smallMolFF', params.EnumParam, default=0, choices=SMALL_MOL_FFS,
                        label='Small molecule force field: ',
                        help='openfe "small_molecule_forcefield". Default openff-2.2.0 is the Open '
                             'Force Field Sage 2.2.0. Protein and water force fields are left at '
                             'openfe defaults (Amber14SB + TIP3P).')
        gGroup.addParam('chargeMethod', params.EnumParam, default=0, choices=CHARGE_METHODS,
                        label='Ligand partial charges: ',
                        help='Charges are assigned ONCE per unique ligand and reused across all '
                             'transformations and repeats (openfe bulk_assign_partial_charges), '
                             'avoiding conformer-dependent charge irreproducibility. am1bcc uses '
                             'AmberTools/Antechamber.')
        gGroup.addParam('temperature', params.FloatParam, default=DEFAULT_TEMPERATURE,
                        label='Temperature (K): ',
                        help='openfe "thermo_settings.temperature".')
        gGroup.addParam('solventPadding', params.FloatParam, default=DEFAULT_SOLVENT_PADDING,
                        label='Solvent padding (nm): ',
                        help='openfe "*_solvation_settings.solvent_padding" - distance from the solute '
                             'to the edge of the DODECAHEDRAL water box. Do not lower below ~1.3 nm: '
                             'that box\'s c.z component is only 0.7071 of its length, so a smaller pad '
                             'makes the half-box shorter than OpenMM\'s 0.9 nm nonbonded cutoff and '
                             'the solvent leg dies outright. Note openfe itself pads the two legs '
                             'differently (1.0 nm complex / 1.5 nm solvent); this single value is '
                             'applied to both.')
        gGroup.addParam('protocolRepeats', params.IntParam, default=DEFAULT_PROTOCOL_REPEATS,
                        label='Independent repeats: ',
                        help='openfe "protocol_repeats" - independent replicas of each transformation, '
                             'averaged into the final estimate. With 1 repeat openfe reports an '
                             'uncertainty of exactly 0, meaning "no estimate".')
        gGroup.addParam('equilLength', params.FloatParam, default=DEFAULT_EQUIL_LENGTH,
                        expertLevel=params.LEVEL_ADVANCED, label='Equilibration per window (ns): ',
                        help='openfe "*_simulation_settings.equilibration_length", per lambda window '
                             '(distinct from the per-leg pre-equilibration above).')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        prepStep = self._insertFunctionStep(self.prepareInputsStep)
        setupStep = self._insertFunctionStep(self.setupStep, prerequisites=[prepStep])
        runStep = self._insertFunctionStep(self.runAllTransformationsStep, prerequisites=[setupStep])
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

    def getSpecifiedMol(self):
        """The one ligand picked on the form, matched by its string representation - the same
        lookup ProtOpenMMSystemPrep.getSpecifiedMolFile uses."""
        for mol in self.inputSetOfMols.get():
            if mol.__str__() == self.inputLigand.get():
                return mol.clone()
        raise ValueError(f'Ligand "{self.inputLigand.get()}" not found in the input set of '
                         f'molecules')

    def prepareInputsStep(self):
        """Write the receptor PDB and an SDF holding the ONE selected ligand - the two inputs the
        setup script consumes (Chem.SDMolSupplier + ProteinComponent.from_pdb_file).

        Hydrogens are (re)added, which is NOT optional for openfe: the OpenFF toolkit builds a
        full valence model and rejects anything it reads as a radical. A ligand extracted straight
        from a crystallographic PDB has no hydrogens at all, so every carbon comes back bare and
        parametrisation dies with "RadicalsNotSupportedError: ... Found 2 radical electrons on
        molecule [C][C]/C([C])=[C]/..." (confirmed on retinal from 1uaz). A ligand that already
        carries explicit hydrogens is unaffected - rdkit only fills in missing ones."""
        self.getReceptorPDB()

        mol = self.getSpecifiedMol()
        molFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
        sdfFile = os.path.abspath(convertToSdf(self, os.path.abspath(molFile)))
        # Not convertToSdf(addHydrogens=True): that flag is skipped entirely when the input is
        # already an .sdf (the function returns before reaching it), which would silently leave
        # an unprotonated ligand for exactly the inputs most likely to need it.
        addHydrogensToMol(self, os.path.dirname(sdfFile), sdfFile)
        mergeSDFs([sdfFile], self.getLigandsSDF())

    # -- setup --------------------------------------------------------------------
    def getTransformDir(self):
        return os.path.abspath(self._getExtraPath('transformations'))

    def getLigandNamesFile(self):
        return os.path.abspath(self._getExtraPath('ligands.txt'))

    def setupStep(self):
        """One openfe Transformation (and JSON) for the selected ligand."""
        os.makedirs(self.getTransformDir(), exist_ok=True)
        paramsFile = os.path.abspath(self._getExtraPath('abfeSetup.txt'))
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
            f.write(f'hostMinDistance :: {self.hostMinDistance.get()}\n')
            f.write(f'hostMaxDistance :: {self.hostMaxDistance.get()}\n')
            f.write(f'productionLength :: {self.productionLength.get()}\n')
            f.write(f'ligandsFile :: {self.getLigandNamesFile()}\n')
            if self.preEquilLength.get() > 0:
                f.write(f'preEquilLength :: {self.preEquilLength.get()}\n')
            if getattr(self, params.USE_GPU).get():
                f.write('computePlatform :: cuda\n')  # openfe's own default spelling is lowercase
                # Without this openfe leaves engine_settings.gpu_device_index at None and OpenMM
                # picks device 0, silently ignoring the GPU chosen on the form.
                gpuIds = getattr(self, params.GPU_LIST).get().replace(',', ' ').split()
                if gpuIds:
                    f.write(f'gpuIndex :: {" ".join(gpuIds)}\n')

        Plugin.runScript(self, 'openfeSetupABFE.py', args=paramsFile, env=OPENMM_DIC,
                         cwd=self._getExtraPath())

    # -- execution ----------------------------------------------------------------
    def getResultsDir(self):
        return os.path.abspath(self._getExtraPath('results'))

    def getResultsFile(self):
        return self._getExtraPath('results.txt')

    def getTransformationNames(self):
        """Names of the transformation JSONs the setup step produced, read from disk."""
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
        """`openfe quickrun` for the selected ligand's transformation."""
        names = self.getTransformationNames()
        if not names:
            raise RuntimeError('The setup step produced no transformations - check the setup log in '
                               'extra/ (a ligand may have failed charge assignment).')
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

    def parseLigandNames(self):
        """(transformationName, ligandName) pairs written by the setup script, so results can be
        labelled with the real ligand name rather than the sanitized file name."""
        pairs = []
        if not os.path.exists(self.getLigandNamesFile()):
            return pairs
        with open(self.getLigandNamesFile()) as f:
            for line in f:
                parts = [p.strip() for p in line.strip().split('::')]
                if len(parts) == 2:
                    pairs.append(tuple(parts))
        return pairs

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
        """Collect every ligand's dG_bind and write extra/results.txt. Returns
        (results, nRequestedLigands). Kept separate from createOutputStep so the result
        bookkeeping can be tested without the project/mapper machinery that registering an
        output object requires."""
        lines, results = [], {}
        ligands = self.parseLigandNames()
        for transName, ligName in ligands:
            dg, err = self.parseResult(transName)
            if dg is None:
                lines.append(f'{ligName}: incomplete (no result)')
                continue
            results[ligName] = (dg, err)
            lines.append(f'{ligName}: dG_bind = {dg:+.2f} +/- {err or 0.0:.2f} kcal/mol')

        header = ('# Absolute binding free energies (OpenFE AbsoluteBindingProtocol)\n'
                  '# dG_bind in kcal/mol; more negative means stronger binding.\n')
        with open(self.getResultsFile(), 'w') as f:
            f.write(header + self.uncertaintyNote() + '\n'.join(lines) + '\n')
        return results, len(ligands)

    def createOutputStep(self):
        from ..objects import OpenMMSystem
        results, nLigands = self.writeResultsFile()

        outSystem = OpenMMSystem(filename=self.getReceptorPDB())
        outSystem.setFreeEnergyFile(self.getResultsFile())
        # Exactly one ligand by construction, so publish its dG_bind as the scalar - but only if
        # it actually produced a result (a failed quickrun leaves nothing to report).
        if nLigands == 1 and results:
            outSystem.setFreeEnergy(list(results.values())[0][0])
        self._defineOutputs(outputSystem=outSystem)
        self._defineSourceRelation(self.inputSetOfMols, outSystem)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        mols = self.inputSetOfMols.get()
        if mols is not None and not mols.isDocked():
            errors.append('The input set of molecules must be docked: binding free energy needs '
                          'ligands posed in the receptor\'s binding site.')
        if self.protocolRepeats.get() < 1:
            errors.append('Independent repeats must be at least 1.')
        if mols is not None and self.inputLigand.get() not in [m.__str__() for m in mols]:
            errors.append(f'Ligand "{self.inputLigand.get()}" not found in the input set of '
                          f'molecules - pick one with the wizard.')
        if self.hostMinDistance.get() >= self.hostMaxDistance.get():
            errors.append('The host anchor min distance must be smaller than the max distance.')
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
                     'slow on CPU.')
        if self.protocolRepeats.get() == 1:
            ws.append('With a single repeat there is no inter-repeat reproducibility estimate, and '
                     'openfe reports an uncertainty of exactly 0 - which means "no estimate", not '
                     '"exact".')

        preEquil = (3 * self.preEquilLength.get() * 2 if self.preEquilLength.get() > 0
                    else 6.55)
        ns = self.protocolRepeats.get() * (44 * self.productionLength.get() + preEquil)
        ws.append(f'One ligand x {self.protocolRepeats.get()} repeat(s): a full double-decoupling '
                 f'calculation over 30 complex + 14 solvent lambda windows, plus ~{preEquil:.2f} ns '
                 f'of per-leg pre-equilibration - roughly {ns:.0f} ns of GPU MD. Absolute binding '
                 f'free energy is far more expensive than the relative (RBFE) protocol; if your '
                 f'ligands share a scaffold, prefer RBFE.')
        ws.append('The OpenFE ABFE protocol was not part of the Baumann et al. 2026 benchmark (that '
                 'study evaluated the RBFE protocol only), so there are no published accuracy '
                 'statistics backing it the way there are for RBFE.')
        return ws

    def _summary(self):
        if self.isFinished() and os.path.exists(self.getResultsFile()):
            with open(self.getResultsFile()) as f:
                return [f.read()]
        return ['The protocol has not finished.']

    def _methods(self):
        if not self.isFinished():
            return []
        return ['Absolute binding free energies were computed with the Open Free Energy (OpenFE) '
                'AbsoluteBindingProtocol: automatic Boresch orientational restraints, annihilation '
                'of ligand electrostatics and decoupling of Lennard-Jones interactions across '
                'complex and solvent legs, sampled with Hamiltonian replica exchange on the OpenMM '
                'engine and analysed with MBAR.']
