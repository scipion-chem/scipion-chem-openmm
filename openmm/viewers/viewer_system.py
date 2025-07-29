# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import numpy as np
import matplotlib.pyplot as plt

import pyworkflow.protocol.params as params

from pwchem.viewers import VmdViewPopen, MDSystemPViewer

from ..objects import OpenMMSystem

PENERGY, TEMP, VOL = 'Potential Energy (kJ/mole)', "Temperature (K)", "Box Volume (nm^3)"

class OpenMMSystemPViewer(MDSystemPViewer):
    """ Visualize the output of OpenMM simulation """
    _label = 'Viewer OpenMM System'
    _targets = [OpenMMSystem]

    def __init__(self, **args):
      super().__init__(**args)

    def _defineReportParams(self, form):
      group = form.addGroup('OpenMM reporter analysis')
      group.addParam('repFeature', params.StringParam, label='Display reporter feature: ', default='',
                     help='Which feature of the reporter to plot'
                     )
      group.addParam('displayReporter', params.LabelParam,
                     label='Plot reporter trajectory analysis: ',
                     help='Plots a graph with the reporter feature chosen over the trajectory')

    def _defineParams(self, form):
      super()._defineParams(form)

      if self.getMDSystem().hasTrajectory():
          self._defineReportParams(form)

    def _getVisualizeDict(self):
      dispDic = super()._getVisualizeDict()
      dispDic.update({'displayReporter': self._showReportParameter})
      return dispDic

    def getMDSystem(self, objType=OpenMMSystem):
        if isinstance(self.protocol, objType):
            return self.protocol
        else:
            return self.protocol.outputSystem

    def getMDFeatures(self):
        mdSystem = self.getMDSystem()
        repFile = mdSystem.getReportFile()
        with open(repFile) as f:
          headers = f.readline().strip()[1:].replace('"', '').split(',')

        idxHeaders = [(i+1, h) for i, h in enumerate(headers[1:])]
        return idxHeaders


    def _showReportParameter(self, paramName=None):
      system = self.getMDSystem()
      repFile = system.getReportFile()

      data = np.loadtxt(repFile, delimiter=',')
      step = data[:, 0]

      valIdx, valName = eval(self.repFeature.get())
      values = data[:, valIdx]
      plt.plot(step, values)
      plt.title(f'"{system.getSystemName()}" trajectory "{valName}"')
      plt.xlabel("Step")
      plt.ylabel(f"{valName}")
      plt.show()


