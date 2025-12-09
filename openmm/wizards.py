# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
This wizard will show the structure of the pdb using a matplotlib viewer
to select the radius of the sphere that contains the protein or a desired zone.
"""

import os

from pwchem.wizards import SelectMultiChainWizard, SelectElementWizard, \
  SelectChainWizardQT, SelectResidueWizardQT, SelectAtomWizardQT, VariableWizard
from pwchem.utils import getBaseName
from pwchem.viewers import PyMolViewer

from openmm.protocols import ProtOpenMMSystemPrep, ProtOpenDuckSimulation, ProtOpenMMSystemSimulationConstantPH
from openmm.viewers import OpenMMSystemPViewer

SelectElementWizard().addTarget(protocol=ProtOpenMMSystemPrep,
                               targets=['inputLigand'],
                               inputs=['inputSetOfMols'],
                               outputs=['inputLigand'])

SelectElementWizard().addTarget(protocol=ProtOpenMMSystemSimulationConstantPH,
                               targets=['inputLigand'],
                               inputs=['inputSetOfMols'],
                               outputs=['inputLigand'])

SelectElementWizard().addTarget(protocol=ProtOpenDuckSimulation,
                               targets=['inputLigand'],
                               inputs=['inputSetOfMols'],
                               outputs=['inputLigand'])

SelectChainWizardQT().addTarget(protocol=ProtOpenDuckSimulation,
                              targets=['intChain'],
                              inputs=['inputSetOfMols'],
                              outputs=['intChain'])

SelectResidueWizardQT().addTarget(protocol=ProtOpenDuckSimulation,
                                targets=['intResidue'],
                                inputs=['inputSetOfMols', 'intChain'],
                                outputs=['intResidue'])

SelectAtomWizardQT().addTarget(protocol=ProtOpenDuckSimulation,
                               targets=['intAtom'],
                               inputs=['inputSetOfMols', 'intChain', 'intResidue'],
                               outputs=['intAtom'])

SelectElementWizard().addTarget(protocol=OpenMMSystemPViewer,
                                targets=['repFeature'],
                                inputs=['getMDFeatures'],
                                outputs=['repFeature'])

class ViewInputComplexWizard(VariableWizard):
  """Visualize the chosen ligand with the correspondant labels"""
  _targets, _inputs, _outputs = [], {}, {}

  def getMol(self, inSet, molName):
    myMol = None
    for mol in inSet:
      if mol.__str__() == molName:
        myMol = mol.clone()
        break
    if myMol == None:
      print('The input ligand is not found')
      return None
    else:
      return myMol

  def writePmlFile(self, pmlFile, topoFile, sysName):
    pmlStr = ''

    topoFile = os.path.abspath(topoFile)
    pmlStr += f'load {topoFile}, {sysName}\nhide spheres, {sysName}\nshow sticks, {sysName}\nhide everything, resn HOH or resn WAT\n'
    pmlStr += 'zoom resn LIG\nselect ligand, resn LIG\nlabel ligand and name CA, "%-s" % (ID)\n' \
              'label ligand and not name CA, "%-s" % (ID)\nhide everything, (resn LIG) and (elem H)\ncolor gray70, not resn LIG'

    with open(pmlFile, 'w') as f:
      f.write(pmlStr)

  def show(self, form, *params):
    protocol = form.protocol
    project = protocol.getProject()

    system = protocol.getMDSystem()
    topoFile = system.getFileName()
    sysName = getBaseName(topoFile)

    pmlsDir = project.getTmpPath()
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(sysName))
    self.writePmlFile(pmlFile, topoFile, sysName)

    pymolV = PyMolViewer(project=project)
    view = pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))[0]
    view.show()

ViewInputComplexWizard().addTarget(protocol=OpenMMSystemPViewer,
                              targets=['viewLigand'],
                              inputs=[],
                              outputs=[])
