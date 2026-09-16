
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
                          DEFAULT_MINIMIZATION_STEPS, MIN_SOLVENT_PADDING, MIN_MLE_EDGES,
                          MIN_MLE_REPEATS)

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
                       help='Alchemical free energy is impractical without a GPU.')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label='Choose GPU IDs')
        form.addParallelSection(threads=4, mpi=1)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False, important=True,
                      help='Set of docked molecules sharing a common receptor, which is taken '
                           'from the set itself. At least 2 are needed - use ABFE for one ligand.')

        form.addSection(label='RBFE settings')
        nGroup = form.addGroup('Alchemical network')
        nGroup.addParam('atomMapper', params.EnumParam, default=0, choices=MAPPERS,
                        display=params.EnumParam.DISPLAY_HLIST, label='Atom mapper: ',
                        help='LOMAP: maximum-common-substructure mapping. Kartograf: '
                             '3D-geometry-based - maps atoms LOMAP will not, but is sensitive to '
                             'how well the input poses are aligned.')
        nGroup.addParam('max3d', params.FloatParam, default=1.0, expertLevel=params.LEVEL_ADVANCED,
                        label='Max 3D distance (A): ',
                        help='Atoms further apart than this in the input poses are not mapped.')
        nGroup.addParam('elementChange', params.BooleanParam, default=False,
                        expertLevel=params.LEVEL_ADVANCED, label='Allow element changes: ',
                        help='Whether the mapper may map atoms of different elements onto each other.')

        sGroup = form.addGroup('Sampling')
        sGroup.addParam('nReplicas', params.IntParam, default=DEFAULT_RBFE_N_REPLICAS,
                        label='Lambda windows: ',
                        help='Lambda windows, sampled with Hamiltonian replica exchange. 11 for '
                             'neutral transformations, 22 for charge-changing ones. A CORRECTNESS '
                             'parameter, not a cost knob - see this protocol\'s warnings before '
                             'lowering it.')
        sGroup.addParam('productionLength', params.FloatParam, default=DEFAULT_RBFE_PRODUCTION,
                        label='Production per window (ns): ',
                        help='Production length per window. 5 ns is the default for neutral '
                             'transformations, and results are ~converged by 80% of that. This IS '
                             'the cost knob - lower it rather than the window count.')
        sGroup.addParam('minimizationSteps', params.IntParam, default=DEFAULT_MINIMIZATION_STEPS,
                        expertLevel=params.LEVEL_ADVANCED, label='Minimization steps per window: ',
                        help='Minimization steps applied to every lambda window. NOT a safe way '
                             'to shorten a run: at 100 steps a freshly-solvated system still '
                             'clashes and the first MD dies. Cut sampling time instead.')

        gGroup = form.addGroup('Force field and thermodynamics')
        gGroup.addParam('smallMolFF', params.EnumParam, default=0, choices=SMALL_MOL_FFS,
                        label='Small molecule force field: ',
                        help='Default openff-2.2.0 is the Open Force Field Sage 2.2.0. Protein '
                             'and water force fields are left at openfe defaults (Amber14SB + '
                             'TIP3P).')
        gGroup.addParam('chargeMethod', params.EnumParam, default=0, choices=CHARGE_METHODS,
                        label='Ligand partial charges: ',
                        help='Charges are assigned ONCE per unique ligand and reused across all '
                             'transformations and repeats, avoiding conformer-dependent charge '
                             'irreproducibility. am1bcc uses AmberTools/Antechamber.')
        gGroup.addParam('temperature', params.FloatParam, default=DEFAULT_TEMPERATURE,
                        label='Temperature (K): ')
        gGroup.addParam('solventPadding', params.FloatParam, default=DEFAULT_SOLVENT_PADDING,
                        label='Solvent padding (nm): ',
                        help='Distance from the solute to the edge of the DODECAHEDRAL water '
                             'box. Do not lower below ~1.3 nm: that box\'s c.z component is only '
                             '0.7071 of its length, so a smaller pad makes the half-box shorter '
                             'than OpenMM\'s 0.9 nm nonbonded cutoff and the solvent leg dies.')
        gGroup.addParam('protocolRepeats', params.IntParam, default=DEFAULT_PROTOCOL_REPEATS,
                        label='Independent repeats: ',
                        help='Independent replicas of each transformation, differing only in '
                             'initial velocities, averaged into the final estimate. With 1 repeat '
                             'openfe reports an uncertainty of exactly 0, meaning "no estimate".')
        gGroup.addParam('equilLength', params.FloatParam, default=DEFAULT_EQUIL_LENGTH,
                        expertLevel=params.LEVEL_ADVANCED, label='Equilibration per window (ns): ',
                        help='Equilibration length, per lambda window.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        prepStep = self._insertFunctionStep(self.prepareInputsStep)
        planStep = self._insertFunctionStep(self.planStep, prerequisites=[prepStep])
        runStep = self._insertFunctionStep(self.runAllTransformationsStep, prerequisites=[planStep])
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

    def prepareInputsStep(self):
        """Receptor PDB + one multi-molecule SDF, the two inputs the setup script consumes.

        Adding hydrogens is not optional: OpenFF builds a full valence model, so a ligand taken
        straight from a crystallographic PDB dies with RadicalsNotSupportedError."""
        self.getReceptorPDB()

        sdfFiles = []
        for mol in self.inputSetOfMols.get():
            mol = mol.clone()
            molFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
            sdfFile = os.path.abspath(convertToSdf(self, os.path.abspath(molFile)))
            # Not convertToSdf(addHydrogens=True): that flag is skipped for .sdf inputs, i.e.
            # exactly the ones most likely to need it.
            addHydrogensToMol(self, os.path.dirname(sdfFile), sdfFile)
            sdfFiles.append(sdfFile)
        mergeSDFs(sdfFiles, self.getLigandsSDF())

    # -- planning -----------------------------------------------------------------
    def getTransformDir(self):
        return os.path.abspath(self._getExtraPath('transformations'))

    def getEdgesFile(self):
        return os.path.abspath(self._getExtraPath('edges.txt'))

    def getLigandNamesFile(self):
        return os.path.abspath(self._getExtraPath('ligands.txt'))

    def getNetworkFile(self):
        """The planned LigandNetwork, openfe's own graphml serialization. Kept project-relative:
        it is stored on the output set and read back by the viewer, which launches
        `openfe view-ligand-network` on it."""
        return self._getExtraPath('ligand_network.graphml')

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
            f.write(f'networkFile :: {os.path.abspath(self.getNetworkFile())}\n')
            f.write(f'edgesFile :: {self.getEdgesFile()}\n')
            f.write(f'ligandsFile :: {self.getLigandNamesFile()}\n')
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
        # runOpenMM, not a dedicated runner: openfe lives in the OpenMM conda env.
        Plugin.runOpenMM(self, 'openfe quickrun', args, cwd=self.getResultsDir())

    def runAllTransformationsStep(self):
        """Sequential: each quickrun already saturates a GPU, so parallel ones only contend."""
        names = self.getTransformationNames()
        if not names:
            raise RuntimeError('The planning step produced no transformations - check that the '
                               'input ligands are congeneric enough for the atom mapper to connect '
                               'them (see the ligand_network.graphml / setup log in extra/).')
        for i, name in enumerate(names, start=1):
            self.info(f'openfe quickrun {i}/{len(names)}: {name}')
            self.quickrunStep(name)

    # -- analysis -----------------------------------------------------------------
    def getGatherFile(self, report):
        return self._getExtraPath(f'gather_{report}.tsv')

    def mleBlockers(self, nEdges=None):
        """Reasons `gather --report dg` cannot run (empty if it can). Both are openfe minimums
        enforced with sys.exit(1); checked here so _warnings() can say so before the GPU time is
        spent. nEdges defaults to the planned network, which only exists after planStep."""
        reasons = []
        if nEdges is None:
            nEdges = len(self.parseEdges())
        if nEdges < MIN_MLE_EDGES:
            reasons.append(f'it needs at least {MIN_MLE_EDGES} network edges and this run has '
                           f'{nEdges}, i.e. {MIN_MLE_EDGES + 1} ligands')
        if self.protocolRepeats.get() < MIN_MLE_REPEATS:
            reasons.append(f'it needs at least {MIN_MLE_REPEATS} repeats per edge and this run '
                           f'has {self.protocolRepeats.get()}, so there is no uncertainty for the '
                           f'maximum-likelihood solve to weight each edge by')
        return reasons

    def gatherStep(self):
        """`openfe gather` once per report type. The `dg` report is more than a re-read of the
        JSONs: it solves the whole network by maximum likelihood, so a ligand's value is informed
        by every path reaching it. That is what the output molecules carry.
        JSONs are listed explicitly because gather walks directories recursively."""
        resultFiles = [os.path.join(self.getResultsDir(), f'{name}.json')
                       for name in self.getTransformationNames()]
        resultFiles = [f for f in resultFiles if os.path.exists(f)]
        if not resultFiles:
            raise RuntimeError('No openfe result JSON was produced - every transformation failed. '
                               'Check the quickrun logs in extra/results/.')

        for report in ('dg', 'ddg', 'raw'):
            if report == 'dg' and self.mleBlockers():
                # Not an error, and not worth running a command that can only exit 1.
                self.info(f'Skipping the "dg" report: {" and ".join(self.mleBlockers())}. The '
                          f'per-ligand RBFE_dG column will be empty; the per-edge ddG values are '
                          f'in extra/results.txt and gather_ddg.tsv.')
                continue
            # --allow-partial: a network in which some edge failed should still yield a table for
            # the edges that worked, rather than aborting the whole analysis.
            args = (f'{" ".join(resultFiles)} --report {report} --allow-partial '
                    f'-o {os.path.abspath(self.getGatherFile(report))}')
            try:
                Plugin.runOpenMM(self, 'openfe gather', args, cwd=self._getExtraPath())
            except Exception as e:
                # Never fatal: results.txt comes from the JSONs and does not depend on gather.
                self.warning(f'`openfe gather --report {report}` failed ({e}). '
                             f'extra/gather_{report}.tsv will be missing.')

    def parseGatherDG(self):
        """{ligandName: (dG, uncertainty)} from the `dg` report. Maximum-likelihood estimates
        centred on zero, so only DIFFERENCES between ligands mean anything."""
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

    @staticmethod
    def parseSetupFile(path, nFields):
        """Rows of a "::"-separated file written by the setup script, as nFields-tuples."""
        if not os.path.exists(path):
            return []
        with open(path) as f:
            rows = [tuple(p.strip() for p in line.strip().split('::')) for line in f]
        return [r for r in rows if len(r) == nFields]

    def parseLigandNames(self):
        """The openfe ligand names in input-set order, written by the setup script."""
        return [name for _, name in self.parseSetupFile(self.getLigandNamesFile(), 2)]

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

    def parseEdges(self):
        """(edgeName, ligandA, ligandB) per planned edge, so the leg JSONs can be paired back up
        without re-deriving the network here."""
        return self.parseSetupFile(self.getEdgesFile(), 3)

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
        """Combine both legs of every edge into extra/results.txt -> (edgeResults, nEdges).
        Separate from createOutputStep so the bookkeeping is testable without a project."""
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

    def matchedLigandNames(self):
        """The openfe ligand name per input molecule, in the set's own order - a positional
        match, since the SDF is written and read in that order. The fallback only matters when a
        ligand was dropped on the way, which shifts every later index."""
        mols = list(self.inputSetOfMols.get())
        names = self.parseLigandNames()
        if len(names) == len(mols):
            return names

        self.warning(f'{len(names)} ligand(s) reached openfe but the input set has {len(mols)} - '
                     f'some were dropped during preparation. Matching results by name instead of '
                     f'by position; molecules that cannot be matched will have no free energy.')
        matched = []
        for mol in mols:
            molName = mol.getMolName() or ''
            hits = [n for n in names if n == molName or molName in n or n in molName]
            matched.append(hits[0] if len(hits) == 1 else None)
        return matched

    def createOutputStep(self):
        import pyworkflow.object as pwobj
        from pwchem.objects import SetOfSmallMolecules

        self.writeResultsFile()
        dgs = self.parseGatherDG()
        names = self.matchedLigandNames()

        outMols = SetOfSmallMolecules.createCopy(self.inputSetOfMols.get(), self._getPath(),
                                                 copyInfo=True)
        for mol, name in zip(self.inputSetOfMols.get(), names):
            mol = mol.clone()
            dg, err = dgs.get(name, (None, None))
            setattr(mol, 'RBFE_dG', pwobj.Float(dg))
            setattr(mol, 'RBFE_err', pwobj.Float(err))
            outMols.append(mol)

        # Set-level metadata: not per-molecule columns, but the files the viewer needs.
        setattr(outMols, '_networkFile', pwobj.String(self.getNetworkFile()))
        setattr(outMols, '_freeEnergyFile', pwobj.String(self.getResultsFile()))

        self._defineOutputs(outputSmallMolecules=outMols)
        self._defineSourceRelation(self.inputSetOfMols, outMols)

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
            errors.append(f'Solvent padding must be at least {MIN_SOLVENT_PADDING} nm.')
        return errors

    def _warnings(self):
        ws = []
        if not getattr(self, params.USE_GPU).get():
            ws.append('Running without a GPU: alchemical free energy calculations are VERY slow '
                     'on CPU.')
        if self.protocolRepeats.get() == 1:
            ws.append('With a single repeat there is no inter-repeat reproducibility estimate, and '
                     'openfe reports an uncertainty of exactly 0 - which means "no estimate", not '
                     '"exact".')
        if self.nReplicas.get() < 8:
            ws.append(f'{self.nReplicas.get()} lambda windows is very few and likely '
                     'counterproductive: consecutive states barely share phase space, so MBAR '
                     'cannot converge and the estimate is meaningless however long it runs. Cut '
                     '"Production per window" instead - that costs precision without destroying '
                     'the alchemical path.')

        mols = self.inputSetOfMols.get()
        nLig = len(mols) if mols is not None else 0
        if nLig:
            # A minimal spanning network over N ligands has N-1 edges, each run as 2 legs.
            nEdges = max(nLig - 1, 1)
            nSims = 2 * nEdges * self.protocolRepeats.get()
            ns = nSims * self.nReplicas.get() * self.productionLength.get()
            ws.append(f'{nLig} ligands -> ~{nEdges} network edge(s) x 2 legs x '
                     f'{self.protocolRepeats.get()} repeat(s) = ~{nSims} replica-exchange '
                     f'simulations, ~{ns:.0f} ns of GPU MD in total.')
            blockers = self.mleBlockers(nEdges=nLig - 1)
            if blockers:
                ws.append('The output molecules will have an EMPTY RBFE_dG column, because openfe '
                         'cannot solve the network for per-ligand free energies: '
                         + '; '.join(blockers) + '. The relative ddG of each edge is still '
                         'computed and written to extra/results.txt and gather_ddg.tsv.')
        ws.append('Charge-changing edges need a larger schedule than the neutral default (22 '
                 'windows, 20 ns/window); this protocol applies one schedule to every edge, so '
                 'raise the window count/production length if your set mixes charge states.')
        return ws

    def _summary(self):
        if not (self.isFinished() and os.path.exists(self.getResultsFile())):
            return ['The protocol has not finished.']
        with open(self.getResultsFile()) as f:
            summary = [f.read()]
        reports = [r for r in ('dg', 'ddg', 'raw') if os.path.exists(self.getGatherFile(r))]
        if reports:
            summary.append('openfe gather reports in extra/: '
                           + ', '.join(f'gather_{r}.tsv' for r in reports))
        return summary

    def _methods(self):
        if not self.isFinished():
            return []
        return ['Relative binding free energies were computed with the Open Free Energy (OpenFE) '
                'hybrid-topology RBFE protocol: an alchemical network was planned over the input '
                f'ligands with the {self.getEnumText("atomMapper")} atom mapper and the LOMAP scorer, '
                'and each edge was run as solvent and complex legs of Hamiltonian replica-exchange '
                'alchemical sampling on the OpenMM engine, analysed with MBAR.']
