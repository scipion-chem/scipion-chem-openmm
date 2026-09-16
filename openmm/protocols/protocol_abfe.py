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
    Absolute binding free energy (ABFE) for ONE ligand in a receptor, wrapping Open Free Energy's
    `AbsoluteBindingProtocol`. The ligand is picked from the input set with a wizard - one ABFE
    protocol per ligand, since each is a full 44-window double-decoupling calculation.

    openfe runs the whole double-decoupling cycle itself - state A is the bound ligand, state B
    the same system without it - deriving both legs, applying automatic Boresch restraints,
    annihilating electrostatics, decoupling LJ, sampling with HREX and estimating with MBAR. This
    protocol contributes input preparation, settings, execution and result collection.

      prepareInputsStep         - receptor PDB + an SDF with the selected ligand
      setupStep                 - openfeSetupABFE.py: charge it, build its Transformation JSON
      runAllTransformationsStep - `openfe quickrun` on that transformation
      gatherStep                - `openfe gather-abfe` as dg / raw, into extra/*.tsv
      createOutputStep          - one-molecule SetOfSmallMolecules with ABFE_dG / ABFE_err
    """
    _label = 'ABFE (OpenFE absolute binding free energy)'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label='Use GPU for execution: ')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label='Choose GPU IDs')
        form.addParallelSection(threads=4, mpi=1)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False, important=True,
                      help='Set of docked molecules.')
        form.addParam('inputLigand', params.StringParam, label='Ligand: ',
                      help='The single molecule to compute the absolute binding free energy for, '
                           'picked from the set above with the wizard.')

        form.addSection(label='ABFE settings')
        rGroup = form.addGroup('Boresch restraints')
        rGroup.addParam('hostMinDistance', params.FloatParam, default=0.5,
                        label='Host anchor min distance (nm): ',
                        help='The closest a receptor atom may be to the ligand to be considered as a Boresch restraint '
                             'anchor. openfe picks the actual anchor atoms itself within this shell, '
                             'using RMSF, secondary structure and distance heuristics.')
        rGroup.addParam('hostMaxDistance', params.FloatParam, default=1.5,
                        label='Host anchor max distance (nm): ',
                        help='The furthest a receptor atom may be to be considered as a restraint anchor.')

        sGroup = form.addGroup('Sampling')
        sGroup.addParam('productionLength', params.FloatParam, default=10.0,
                        label='Production per window (ns): ',
                        help='Production length for both legs. '
                             'openfe\'s ABFE default is 10 ns per replica, across 30 complex + 14 '
                             'solvent lambda windows. ABFE convergence is dominated by the '
                             'near-fully-decoupled windows, so reduce with care.')
        sGroup.addParam('minimizationSteps', params.IntParam, default=DEFAULT_MINIMIZATION_STEPS,
                        expertLevel=params.LEVEL_ADVANCED, label='Minimization steps per window: ',
                        help='Minimization steps applied to every '
                             'one of the 44 windows. NOT a safe way to shorten a run: at 100 steps '
                             'a freshly-solvated system still clashes and the first MD dies. Cut sampling time instead.')
        sGroup.addParam('preEquilLength', params.FloatParam, default=0.0,
                        expertLevel=params.LEVEL_ADVANCED,
                        label='Pre-equilibration per leg: ',
                        help='Before the alchemical windows, openfe runs a plain MD pre-equilibration '
                             'of each leg (NVT + NPT equilibration). Its defaults are '
                             'large and asymmetric - 0.25+0.5+5.0 ns for the complex leg and '
                             '0.1+0.2+0.5 ns for the solvent leg.'
                             'Leave at 0 to keep openfe\'s own values; set a value to use it for all '
                             'three phases of BOTH legs.')

        gGroup = form.addGroup('Force field and thermodynamics')
        gGroup.addParam('smallMolFF', params.EnumParam, default=0, choices=SMALL_MOL_FFS,
                        label='Small molecule force field: ',
                        help='Default openff-2.2.0 is the Open '
                             'Force Field Sage 2.2.0. Protein and water force fields are left at '
                             'openfe defaults (Amber14SB + TIP3P).')
        gGroup.addParam('chargeMethod', params.EnumParam, default=0, choices=CHARGE_METHODS,
                        label='Ligand partial charges: ',
                        help='Charges are assigned ONCE per unique ligand and reused across all '
                             'transformations and repeats, avoiding conformer-dependent charge irreproducibility. '
                             'am1bcc uses AmberTools/Antechamber.')
        gGroup.addParam('temperature', params.FloatParam, default=DEFAULT_TEMPERATURE,
                        label='Temperature (K): ')
        gGroup.addParam('solventPadding', params.FloatParam, default=DEFAULT_SOLVENT_PADDING,
                        label='Solvent padding (nm): ',
                        help='Distance from the solute '
                             'to the edge of the DODECAHEDRAL water box. Do not lower below ~1.3 nm: '
                             'that box\'s c.z component is only 0.7071 of its length, so a smaller pad '
                             'makes the half-box shorter than OpenMM\'s 0.9 nm nonbonded cutoff and '
                             'the solvent leg dies outright. Note openfe itself pads the two legs '
                             'differently (1.0 nm complex / 1.5 nm solvent); this single value is '
                             'applied to both.')
        gGroup.addParam('protocolRepeats', params.IntParam, default=DEFAULT_PROTOCOL_REPEATS,
                        label='Independent repeats: ',
                        help='Independent replicas of each transformation, '
                             'averaged into the final estimate. With 1 repeat openfe reports an '
                             'uncertainty of exactly 0, meaning "no estimate".')
        gGroup.addParam('equilLength', params.FloatParam, default=DEFAULT_EQUIL_LENGTH,
                        expertLevel=params.LEVEL_ADVANCED, label='Equilibration per window (ns): ',
                        help='Equilibration length, per lambda window '
                             '(distinct from the per-leg pre-equilibration above).')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        prepStep = self._insertFunctionStep(self.prepareInputsStep)
        setupStep = self._insertFunctionStep(self.setupStep, prerequisites=[prepStep])
        runStep = self._insertFunctionStep(self.runAllTransformationsStep, prerequisites=[setupStep])
        gatherStep = self._insertFunctionStep(self.gatherStep, prerequisites=[runStep])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[gatherStep])

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
        """Receptor PDB + an SDF with the ONE selected ligand, the two inputs the setup script
        consumes."""
        self.getReceptorPDB()

        mol = self.getSpecifiedMol()
        molFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
        sdfFile = os.path.abspath(convertToSdf(self, os.path.abspath(molFile)))
        # Not convertToSdf(addHydrogens=True): that flag is skipped for .sdf inputs, i.e. exactly
        # the ones most likely to need it.
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

    # -- analysis -----------------------------------------------------------------
    def getGatherFile(self, report):
        return self._getExtraPath(f'gather_{report}.tsv')

    def gatherStep(self):
        """`openfe gather-abfe` once per report type. The `dg` report is where ABFE_dG comes from,
        because it picks its error column by repeat count: MBAR uncertainty for one repeat, the
        standard deviation for several.
        ABFE gathering is still experimental, so createOutputStep falls back to the JSON."""
        resultFiles = [os.path.join(self.getResultsDir(), f'{name}.json')
                       for name in self.getTransformationNames()]
        resultFiles = [f for f in resultFiles if os.path.exists(f)]
        if not resultFiles:
            raise RuntimeError('No openfe result JSON was produced - the transformation failed. '
                               'Check the quickrun log in extra/results/.')

        for report in ('dg', 'raw'):
            # No `ddg` report here - nothing is relative - but --allow-partial still matters:
            # if one leg of the cycle failed, the other is still tabulated.
            args = (f'{" ".join(resultFiles)} --report {report} --allow-partial '
                    f'-o {os.path.abspath(self.getGatherFile(report))}')
            try:
                Plugin.runOpenMM(self, 'openfe gather-abfe', args, cwd=self._getExtraPath())
            except Exception as e:
                self.warning(f'`openfe gather-abfe --report {report}` failed ({e}). '
                             f'extra/gather_{report}.tsv will be missing.')

    def parseGatherDG(self):
        """{ligandName: (dG_bind, uncertainty)} from the `dg` report. Read positionally, not by
        header: the uncertainty column is named "MBAR uncertainty" for one repeat and "std dev
        uncertainty" for several, but its position never moves."""
        values = {}
        gatherFile = self.getGatherFile('dg')
        if not os.path.exists(gatherFile):
            return values
        with open(gatherFile) as f:
            next(f, None)
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 3 and parts[0].strip():
                    values[parts[0].strip()] = (self._asFloat(parts[1]), self._asFloat(parts[2]))
        return values

    # -- results ------------------------------------------------------------------
    def parseResult(self, name):
        """(estimate, uncertainty) in kcal/mol from a quickrun result JSON, or (None, None).
        Never raises: one failed transformation must not take down the output step. Quantities
        come as {'magnitude': x, ...} or as plain numbers depending on the version."""
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
        """(transformationName, ligandName) pairs written by the setup script, so results carry the
        real ligand name rather than the sanitized file name."""
        if not os.path.exists(self.getLigandNamesFile()):
            return []
        with open(self.getLigandNamesFile()) as f:
            rows = [tuple(p.strip() for p in line.strip().split('::')) for line in f]
        return [r for r in rows if len(r) == 2]

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
        """Collect every ligand's dG_bind into extra/results.txt -> (results, nLigands).
        Separate from createOutputStep so the bookkeeping is testable without a project."""
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

    def publishedFreeEnergy(self, jsonResults):
        """The (dG_bind, uncertainty) published on the output molecule, in kcal/mol.

        Prefers the gather `dg` report, falls back to the result JSON. Same dG_bind either way;
        gather is preferred for its ERROR, which is the MBAR uncertainty rather than the JSON's
        0.0 for a single repeat. The fallback keeps a finished simulation from losing its result
        to an experimental reporter."""
        gathered = self.parseGatherDG()
        # One ligand by construction, so any row gather produced is this one - no name matching.
        for dg, err in gathered.values():
            if dg is not None:
                return dg, err

        if gathered:
            self.warning('The gather-abfe "dg" report had no usable value; falling back to the '
                         'result JSON.')
        elif os.path.exists(self.getResultsDir()):
            self.info('No gather-abfe "dg" report; taking dG_bind from the result JSON instead. '
                      'Note its uncertainty is 0.0 for a single-repeat run, meaning "no estimate".')
        # A failed quickrun leaves nothing at all, in which case the columns stay empty.
        return list(jsonResults.values())[0] if jsonResults else (None, None)

    def createOutputStep(self):
        import pyworkflow.object as pwobj
        from pwchem.objects import SetOfSmallMolecules

        results, _ = self.writeResultsFile()
        dg, err = self.publishedFreeEnergy(results)

        outMols = SetOfSmallMolecules.createCopy(self.inputSetOfMols.get(), self._getPath(),
                                                 copyInfo=True)
        mol = self.getSpecifiedMol()

        setattr(mol, 'ABFE_dG', pwobj.Float(dg))
        setattr(mol, 'ABFE_err', pwobj.Float(err))
        outMols.append(mol)

        setattr(outMols, '_freeEnergyFile', pwobj.String(self.getResultsFile()))

        self._defineOutputs(outputSmallMolecules=outMols)
        self._defineSourceRelation(self.inputSetOfMols, outMols)

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
                f'Solvent padding must be at least {MIN_SOLVENT_PADDING} nm.')
        return errors

    def _warnings(self):
        ws = []
        if not getattr(self, params.USE_GPU).get():
            ws.append('Running without a GPU: alchemical free energy calculations are VERY '
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
        return ws

    def _summary(self):
        if not (self.isFinished() and os.path.exists(self.getResultsFile())):
            return ['The protocol has not finished.']
        with open(self.getResultsFile()) as f:
            summary = [f.read()]
        reports = [r for r in ('dg', 'raw') if os.path.exists(self.getGatherFile(r))]
        if reports:
            summary.append('openfe gather-abfe reports in extra/: '
                           + ', '.join(f'gather_{r}.tsv' for r in reports))
        return summary

    def _methods(self):
        if not self.isFinished():
            return []
        return ['Absolute binding free energies were computed with the Open Free Energy (OpenFE) '
                'AbsoluteBindingProtocol: automatic Boresch orientational restraints, annihilation '
                'of ligand electrostatics and decoupling of Lennard-Jones interactions across '
                'complex and solvent legs, sampled with Hamiltonian replica exchange on the OpenMM '
                'engine and analysed with MBAR.']
