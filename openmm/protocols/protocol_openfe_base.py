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
Not a runnable protocol: it has no _label and is absent from protocols.conf, so it never appears
in the GUI tree. Just provide shared methods for ABFE and RBFE protocols
"""

import os
import json

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import getBaseName

from .. import Plugin
from ..constants import MIN_SOLVENT_PADDING


class ProtOpenFEBase(EMProtocol):
    """Shared execution and result-parsing helpers for the OpenFE free energy protocols."""

    # -- form helpers -------------------------------------------------------------
    def _defineGpuParams(self, form):
        """The hidden GPU block every openfe protocol needs; the visible form stays in the
        protocol itself."""
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label='Use GPU for execution: ',
                       help='Alchemical free energy is impractical without a GPU.')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label='Choose GPU IDs')
        form.addParallelSection(threads=4, mpi=1)

    def writeCommonSetupParams(self, f):
        """The params-file keys both setup scripts read. Protocol-specific keys are written by
        the protocol around this call."""
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
        f.write(f'productionLength :: {self.productionLength.get()}\n')
        f.write(f'ligandsFile :: {self.getLigandNamesFile()}\n')

    def writeGpuParams(self, f):
        """GPU lines of a setup-script params file."""
        if getattr(self, params.USE_GPU).get():
            f.write('computePlatform :: cuda\n')  # openfe's own default spelling is lowercase
            # Without this openfe leaves engine_settings.gpu_device_index at None and OpenMM
            # picks device 0, silently ignoring the GPU chosen on the form.
            gpuIds = getattr(self, params.GPU_LIST).get().replace(',', ' ').split()
            if gpuIds:
                f.write(f'gpuIndex :: {" ".join(gpuIds)}\n')

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

    # -- paths --------------------------------------------------------------------
    def getTransformDir(self):
        return os.path.abspath(self._getExtraPath('transformations'))

    def getLigandNamesFile(self):
        return os.path.abspath(self._getExtraPath('ligands.txt'))

    def getResultsDir(self):
        return os.path.abspath(self._getExtraPath('results'))

    def getResultsFile(self):
        return self._getExtraPath('results.txt')

    def getGatherFile(self, report):
        return self._getExtraPath(f'gather_{report}.tsv')

    # -- execution ----------------------------------------------------------------
    def getTransformationNames(self):
        """Names of the transformation JSONs the setup step produced. Read from disk, not
        precomputed in _insertAllSteps, because the setup step decides how many there are."""
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

    def runAllTransformationsStep(self, noneMsg='The setup step produced no transformations.'):
        """Sequential: each quickrun already saturates a GPU, so parallel ones only contend."""
        names = self.getTransformationNames()
        if not names:
            raise RuntimeError(noneMsg)
        for i, name in enumerate(names, start=1):
            self.info(f'openfe quickrun {i}/{len(names)}: {name}')
            self.quickrunStep(name)

    def existingResultFiles(self):
        """The quickrun result JSONs that exist, listed explicitly rather than handing gather the
        results directory: it walks recursively, and the working directories hold thousands of
        files."""
        files = [os.path.join(self.getResultsDir(), f'{name}.json')
                 for name in self.getTransformationNames()]
        return [f for f in files if os.path.exists(f)]

    # -- result parsing -----------------------------------------------------------
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

    def parseGatherDG(self):
        """{ligandName: (dG, uncertainty)} from the `dg` report. Read positionally, not by header:
        the uncertainty column is named differently depending on the repeat count, but its
        position never moves."""
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

    def uncertaintyNote(self):
        """openfe derives a transformation's uncertainty from the SPREAD ACROSS REPEATS, so with a
        single repeat it reports exactly 0.0 - which reads like a precise number when it really
        means "no uncertainty estimate exists"."""
        if self.protocolRepeats.get() > 1:
            return ''
        return ('# NOTE: run with a single repeat, so every "+/- 0.00" below means NO uncertainty\n'
                '# estimate was available, not a precise result - openfe derives it from the\n'
                '# spread across repeats. Use 3 repeats (the default) for a real error bar.\n')

    # -- info ---------------------------------------------------------------------
    def validateCommon(self):
        """Checks both protocols make; each appends its own on top."""
        errors = []
        mols = self.inputSetOfMols.get()
        if mols is not None and not mols.isDocked():
            errors.append('The input set of molecules must be docked: binding free energy needs '
                          'ligands posed in the receptor\'s binding site.')
        if self.protocolRepeats.get() < 1:
            errors.append('Independent repeats must be at least 1.')
        if self.solventPadding.get() < MIN_SOLVENT_PADDING:
            errors.append(f'Solvent padding must be at least {MIN_SOLVENT_PADDING} nm.')
        return errors

    def commonWarnings(self):
        """Warnings both protocols raise; each appends its own cost estimate on top."""
        ws = []
        if not getattr(self, params.USE_GPU).get():
            ws.append('Running without a GPU: alchemical free energy calculations are VERY '
                      'slow on CPU.')
        if self.protocolRepeats.get() == 1:
            ws.append('With a single repeat there is no inter-repeat reproducibility estimate, and '
                      'openfe reports an uncertainty of exactly 0 - which means "no estimate", not '
                      '"exact".')
        return ws

    def summaryWithReports(self, reports, gatherName='openfe gather'):
        """results.txt plus whichever gather tables were produced."""
        if not (self.isFinished() and os.path.exists(self.getResultsFile())):
            return ['The protocol has not finished.']
        with open(self.getResultsFile()) as f:
            summary = [f.read()]
        found = [r for r in reports if os.path.exists(self.getGatherFile(r))]
        if found:
            summary.append(f'{gatherName} reports in extra/: '
                           + ', '.join(f'gather_{r}.tsv' for r in found))
        return summary
