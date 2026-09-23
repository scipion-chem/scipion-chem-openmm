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
import math

from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem.utils import convertToSdf, mergeSDFs, addHydrogensToMol

from .. import Plugin
from .protocol_openfe_base import ProtOpenFEBase
from ..constants import (OPENMM_DIC, DEFAULT_TEMPERATURE, DEFAULT_SOLVENT_PADDING,
                          DEFAULT_PROTOCOL_REPEATS, DEFAULT_EQUIL_LENGTH,
                          DEFAULT_RBFE_PRODUCTION, DEFAULT_RBFE_N_REPLICAS,
                          DEFAULT_MINIMIZATION_STEPS, MIN_MLE_EDGES,
                          MIN_MLE_REPEATS)

SMALL_MOL_FFS = ['openff-2.2.0', 'openff-2.1.0', 'openff-2.0.0', 'gaff-2.11', 'espaloma-0.3.2']
CHARGE_METHODS = ['am1bcc', 'am1bccelf10', 'nagl', 'espaloma']
MAPPERS = ['LOMAP', 'Kartograf']
COMPLEX, SOLVENT = 'complex', 'solvent'


class ProtOpenFERBFE(ProtOpenFEBase):
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
        self._defineGpuParams(form)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False, important=True,
                      help='Set of docked molecules sharing a common receptor, which is taken '
                           'from the set itself. At least 2 are needed - use ABFE for one ligand.')
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
        sGroup.addParam('protocolRepeats', params.IntParam, default=DEFAULT_PROTOCOL_REPEATS,
                        label='Independent repeats: ',
                        help='Independent replicas of each transformation, differing only in '
                             'initial velocities, averaged into the final estimate. With 1 repeat '
                             'openfe reports an uncertainty of exactly 0, meaning "no estimate".')
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
        sGroup.addParam('equilLength', params.FloatParam, default=DEFAULT_EQUIL_LENGTH,
                        expertLevel=params.LEVEL_ADVANCED, label='Equilibration per window (ns): ',
                        help='Equilibration length, per lambda window.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        prepStep = self._insertFunctionStep(self.prepareInputsStep)
        planStep = self._insertFunctionStep(self.planStep, prerequisites=[prepStep])
        runStep = self._insertFunctionStep(self.runAllTransformationsStep, prerequisites=[planStep])
        gatherStep = self._insertFunctionStep(self.gatherStep, prerequisites=[runStep])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[gatherStep])

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
    def getEdgesFile(self):
        return os.path.abspath(self._getExtraPath('edges.txt'))

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
            self.writeCommonSetupParams(f)
            f.write(f'mapper :: {self.getEnumText("atomMapper")}\n')
            f.write(f'max3d :: {self.max3d.get()}\n')
            f.write(f'elementChange :: {self.elementChange.get()}\n')
            f.write(f'nReplicas :: {self.nReplicas.get()}\n')
            f.write(f'networkFile :: {os.path.abspath(self.getNetworkFile())}\n')
            f.write(f'edgesFile :: {self.getEdgesFile()}\n')
            self.writeGpuParams(f)

        Plugin.runScript(self, 'openfeSetupRBFE.py', args=paramsFile, env=OPENMM_DIC,
                         cwd=self._getExtraPath())

    # -- execution ----------------------------------------------------------------
    def runAllTransformationsStep(self):
        super().runAllTransformationsStep(
            noneMsg='The planning step produced no transformations - check that the input ligands '
                    'are congeneric enough for the atom mapper to connect them (see the '
                    'ligand_network.graphml / setup log in extra/).')

    # -- analysis -----------------------------------------------------------------
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
        resultFiles = self.existingResultFiles()
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

    def parseLigandNames(self):
        """The openfe ligand names in input-set order, written by the setup script."""
        return [name for _, name in self.parseSetupFile(self.getLigandNamesFile(), 2)]

    # -- results ------------------------------------------------------------------
    def parseEdges(self):
        """(edgeName, ligandA, ligandB) per planned edge, so the leg JSONs can be paired back up
        without re-deriving the network here."""
        return self.parseSetupFile(self.getEdgesFile(), 3)

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
        errors = self.validateCommon()
        mols = self.inputSetOfMols.get()
        if mols is not None and len(mols) < 2:
            errors.append('RBFE needs at least 2 ligands to form an alchemical transformation. '
                          'Use the ABFE protocol for a single ligand.')
        if self.nReplicas.get() < 2:
            errors.append('At least 2 lambda windows are needed.')
        return errors

    def _warnings(self):
        ws = self.commonWarnings()
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
        return self.summaryWithReports(('dg', 'ddg', 'raw'))

    def _methods(self):
        if not self.isFinished():
            return []
        return ['Relative binding free energies were computed with the Open Free Energy (OpenFE) '
                'hybrid-topology RBFE protocol: an alchemical network was planned over the input '
                f'ligands with the {self.getEnumText("atomMapper")} atom mapper and the LOMAP scorer, '
                'and each edge was run as solvent and complex legs of Hamiltonian replica-exchange '
                'alchemical sampling on the OpenMM engine, analysed with MBAR.']
