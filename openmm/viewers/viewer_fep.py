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
import subprocess

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.gui.text import openTextFileEditor

from .. import Plugin as openmmPlugin
from ..constants import OPENMM_DIC
from ..protocols import ProtOpenFERBFE, ProtOpenFEABFE

RBFE_REPORTS = ['dg', 'ddg', 'raw']
ABFE_REPORTS = ['dg', 'raw']


class OpenFEFreeEnergyViewer(pwviewer.ProtocolViewer):
  """ Visualize the ligand network and the analysis reports of an OpenFE free energy run """
  _label = 'Viewer OpenFE free energy'
  _targets = [ProtOpenFERBFE, ProtOpenFEABFE]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def isRBFE(self):
    return isinstance(self.protocol, ProtOpenFERBFE)

  def getReports(self):
    return RBFE_REPORTS if self.isRBFE() else ABFE_REPORTS

  def _defineParams(self, form):
    form.addSection(label='Free energy results')

    if self.isRBFE():
      # Only RBFE has a network: ABFE treats each ligand independently
      group = form.addGroup('Ligand network')
      group.addParam('displayNetwork', LabelParam, label='View ligand network: ',
                     help='Opens openfe\'s own interactive network window on the planned '
                          '.graphml. Nodes are the ligands; clicking an edge draws both '
                          'structures with the mapped atoms colour-matched, which is the only '
                          'way to check the atom mapping that the whole calculation rests on.')

    group = form.addGroup('Analysis reports')
    group.addParam('report', EnumParam, label='Report: ', default=0,
                   choices=self.getReports(), display=EnumParam.DISPLAY_HLIST,
                   help='The `openfe gather` tables written to extra/.\n'
                        'dg: one row per ligand, maximum-likelihood free energy over the whole '
                        'network with its uncertainty. Centred on zero, so only differences '
                        'between ligands mean anything. This is the column on the output '
                        'molecules.\n'
                        'ddg: one row per edge, the relative free energy of that transformation.\n'
                        'raw: every repeat and every leg separately, uncombined - use this to see '
                        'whether the repeats actually agree, and how large each repeat\'s own '
                        'MBAR uncertainty is.')
    group.addParam('displayReport', LabelParam, label='Open report: ',
                   help='Opens the selected gather report as a tab-separated table')
    group.addParam('displaySummary', LabelParam, label='Open results summary: ',
                   help='Opens extra/results.txt, the free energies computed directly from the '
                        'result JSONs rather than through gather')

  def _getVisualizeDict(self):
    return {'displayNetwork': self._viewNetwork,
            'displayReport': self._viewReport,
            'displaySummary': self._viewSummary}

  # ------------------------------- views ---------------------------------------
  def _viewNetwork(self, paramName=None):
    netFile = os.path.abspath(self.protocol.getNetworkFile())
    if not os.path.exists(netFile):
      return [self.errorMessage(f'No ligand network found at {netFile}. The planning step must '
                                f'have failed before writing it.', title='Network not found')]

    # openfe's viewer is an interactive matplotlib window.
    cmd = (f'{openmmPlugin.getEnvActivationCommand(OPENMM_DIC)} && '
           f'openfe view-ligand-network {netFile}')
    subprocess.Popen(cmd, shell=True, cwd=os.path.dirname(netFile))
    return []

  def _viewReport(self, paramName=None):
    report = self.getReports()[self.report.get()]
    reportFile = os.path.abspath(self.protocol.getGatherFile(report))
    if not os.path.exists(reportFile):
      return [self.errorMessage(f'No "{report}" report at {reportFile}. `openfe gather` may have '
                                f'failed - check the protocol log.', title='Report not found')]
    openTextFileEditor(reportFile)
    return []

  def _viewSummary(self, paramName=None):
    resultsFile = os.path.abspath(self.protocol.getResultsFile())
    if not os.path.exists(resultsFile):
      return [self.errorMessage(f'No results file at {resultsFile}.', title='Results not found')]
    openTextFileEditor(resultsFile)
    return []
