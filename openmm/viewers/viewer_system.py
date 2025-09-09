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

import os, csv, subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

import pyworkflow.protocol.params as params

from pwchem.viewers import VmdViewPopen, MDSystemPViewer
from pwchem.constants import RDKIT_DIC

from .. import Plugin as openmmPlugin
from ..objects import OpenMMSystem

PENERGY, TEMP, VOL = 'Potential Energy (kJ/mole)', "Temperature (K)", "Box Volume (nm^3)"


def parseResidueLine(sline, outDic, nFrames):
    frame = int(sline[12])

    resNr, resType, resChain = sline[1:4]
    resId = f'{resChain}:{resType}_{resNr}'

    intType = sline[13]

    if resId not in outDic:
      outDic[resId] = {i + 1: [] for i in range(nFrames)}
    outDic[resId][frame].append(intType)
    return outDic


def parseLigandLine(sline, outDic, nFrames):
    frame = int(sline[12])
    intType = sline[13]

    if intType == 'hbond':
      protIsDon = sline[18]
      ligIds = [sline[21] if protIsDon.upper() == 'TRUE' else sline[19]]
    elif intType == 'waterbridge':
      protIsDon = sline[18]
      ligIds = [sline[27] if protIsDon.upper() == 'TRUE' else sline[26]]
    elif intType in ['saltbridge', 'pistacking', 'pication']:
      ligIds = sline[33].split(',')
    elif intType in ['halogen']:
      ligIds = sline[40].split(',')
    else:
      ligIds = [sline[8]]

    for ligId in ligIds:
      if type(ligId) == str:
        ligId = eval(ligId)
      ligId = int(ligId)
      if ligId not in outDic:
        outDic[ligId] = {i + 1: [] for i in range(nFrames)}
      outDic[ligId][frame].append(intType)

    return outDic

def parseInteractionsLine(sline, outDic, nFrames):
    frame = int(sline[12])
    intType = sline[13]

    resNr, resType, resChain = sline[1:4]
    resId = f'{resChain}:{resType}_{resNr}'

    if intType == 'hbond':
      protIsDon = sline[18]
      ligIds = [sline[21] if protIsDon.upper() == 'TRUE' else sline[19]]
    elif intType == 'waterbridge':
      protIsDon = sline[18]
      ligIds = [sline[27] if protIsDon.upper() == 'TRUE' else sline[26]]
    elif intType in ['saltbridge', 'pistacking', 'pication']:
      ligIds = sline[33].split(',')
    elif intType in ['halogen']:
      ligIds = sline[40].split(',')
    else:
      ligIds = [sline[8]]

    if intType in ['pistacking', 'pication']:
      ligIds = ['_'.join(ligIds)]

    for ligId in ligIds:
      pairId = (resId, ligId)
      if pairId not in outDic:
        outDic[pairId] = {}

      if not intType in outDic[pairId]:
        outDic[pairId][intType] = [frame]

      if not frame in outDic[pairId][intType]:
        outDic[pairId][intType].append(frame)

    return outDic


def cleanChain(d, pair=False):
  if not pair:
      chains = [rid.split(':')[0] for rid in d]
      if len(set(chains)) == 1:
        d = {rid.split(':')[1]: v for rid, v in d.items()}
  else:
      chains = [rid[0].split(':')[0] for rid in d]
      if len(set(chains)) == 1:
        d = {(rid[0].split(':')[1], rid[1]): v for rid, v in d.items()}
  return d

def heatmap(data, rowLabels, colLabels, totalTime):
  """
  Create a heatmap interactions heatmap for residues or ligand ids
  """

  # Get the range of integer values
  vmin, vmax = np.min(data), np.max(data)
  n_values = vmax - vmin + 1

  # Create discrete colormap from a continuous colormap
  continuous_cmap = plt.cm.YlGn  # or 'plasma', 'inferno', 'magma', 'coolwarm', etc.
  discrete_cmap = ListedColormap(continuous_cmap(np.linspace(0, 1, n_values)))

  # Create boundaries at half-integers
  bounds = np.arange(vmin - 0.5, vmax + 1.5, 1)
  norm = BoundaryNorm(bounds, n_values)

  fig, ax = plt.subplots(figsize=(10, 8))
  im = ax.imshow(data, cmap=discrete_cmap, norm=norm, interpolation='none', aspect='auto')

  # Create discrete colorbar
  cbar = plt.colorbar(im, ax=ax, ticks=np.arange(vmin, vmax + 1))
  cbar.set_ticklabels([str(i) for i in range(vmin, vmax + 1)])

  # Show all ticks and label them with the respective list entries.
  ax.set_yticks(np.arange(data.shape[0]), labels=rowLabels)

  # Calculate time values and create 4 equally spaced y-ticks
  nFrames = len(colLabels)
  stepSize = totalTime / nFrames
  x_tick_positions = np.linspace(0, nFrames - 1, 4)  # 4 equally spaced frame positions
  x_tick_times = x_tick_positions * stepSize  # Convert to time

  # Set y-ticks with time labels
  ax.set_xticks(x_tick_positions)
  ax.set_xticklabels([f'{time:.3f} ps' for time in x_tick_times])

  # Let the horizontal axes labeling appear on top.
  ax.tick_params(top=True, bottom=False,
                 labeltop=True, labelbottom=False)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
           rotation_mode="anchor")

  return im, cbar

def histogram(resContacts, nFrames):
  # Extract residue IDs and interaction types
  residues = list(resContacts.keys())
  interaction_types = sorted(set().union(*[set(interactions.keys()) for interactions in resContacts.values()]))

  # Prepare data for stacking
  stack_data = {itype: [] for itype in interaction_types}
  for residue in residues:
    for itype in interaction_types:
      stack_data[itype].append(resContacts[residue].get(itype, 0))

  # Convert to numpy arrays for easier handling
  bottom = np.zeros(len(residues))
  fig, ax = plt.subplots(figsize=(12, 8))

  # Create stacked bars
  colors = plt.cm.Set1(np.linspace(0, 1, len(interaction_types)))
  bars = []

  for i, (itype, values) in enumerate(stack_data.items()):
    values = [v / nFrames for v in values]
    bar = ax.bar(residues, values, bottom=bottom, label=itype, color=colors[i])
    bars.append(bar)
    bottom += np.array(values)

  plt.xlabel('Residue ID')
  plt.ylabel('Interaction Count')
  plt.title('Residue Interaction Counts by Type')
  plt.legend()
  plt.xticks(rotation=45)
  plt.tight_layout()
  plt.show()


class OpenMMSystemPViewer(MDSystemPViewer):
    """ Visualize the output of OpenMM simulation """
    _label = 'Viewer OpenMM System'
    _targets = [OpenMMSystem]

    def __init__(self, **args):
      super().__init__(**args)

    def _defineReportParams(self, form):
      group = form.addGroup('OpenMM reporter analysis')
      group.addParam('repFeature', params.StringParam, label='Display reporter feature: ', default='',
                     help='Which feature of the reporter to plot')
      group.addParam('displayReporter', params.LabelParam, label='Plot reporter trajectory analysis: ',
                     help='Plots a graph with the reporter feature chosen over the trajectory')

      form.addSection('Receptor-ligand interactions')
      group = form.addGroup('OpenMMDL analysis')
      group.addParam('openmmdlAnalysis', params.EnumParam, label='OpenMMDL analysis: ', default=0,
                     choices=['RMSD', 'Barcodes', 'Binding Modes Markov States'],
                     help='Show the OpenMMDL interaction Markov States generated')
      group.addParam('barcodeType', params.EnumParam, label='Which barcode to display: ',
                     condition='openmmdlAnalysis==1', choices=self.getBarcodeTypes(),
                     help='Which feature of the barcodes to plot')
      group.addParam('displayOpenMMDL', params.LabelParam, label='Display OpenMMDL analysis: ',
                     help='Show the OpenMMDL barcodes, RMSD or interaction Markov states generated')

      group = form.addGroup('Receptor-ligand interactions')
      group.addParam('threshold', params.FloatParam, label='Interaction threshold: ', default=0.1,
                     help='Proportion of time through the simulation that a interaction must appear to be considered')

      group.addParam('displayInteractions', params.LabelParam, label='Display interactions diagram: ',
                     help='Displays the ligand in 2D with the receptor residues interactions')

      group.addParam('target', params.EnumParam, label='Show interactions of: ', default=0,
                     choices=['Receptor', 'Ligand'], display=params.EnumParam.DISPLAY_HLIST,
                     help='Whether to show the interactions of the receptor residues or the ligand atoms')

      group.addParam('viewLigand', params.LabelParam, label='View ligand: ', condition='target==1',
                     help='Displays the ligand in PyMol with the annotated atom numbers')
      group.addParam('displayHeatmap', params.LabelParam, label='Display interactions over time: ',
                     help='Display the number of interactions of each residue/atom over the simulation time')
      group.addParam('displayHistogram', params.LabelParam, label='Display interactions type: ',
                     help='Display the type of interactions of each residue/atom through the simulation')


    def _defineParams(self, form):
      super()._defineParams(form)

      if self.getMDSystem().hasTrajectory():
          self._defineReportParams(form)

    def _getVisualizeDict(self):
      dispDic = super()._getVisualizeDict()
      dispDic.update(
        {'displayReporter': self._showReportParameter,
         'displayOpenMMDL': self.showOpenMMDL,
         'displayHeatmap': self.showHeatmap,
         'displayHistogram': self.showHistogram,
         'displayInteractions': self.showInteractionsDiagram
         })
      return dispDic

    def _showReportParameter(self, paramName=None):
      system = self.getMDSystem()
      repFile = system.getReportFile()

      data = np.loadtxt(repFile, delimiter=',', ndmin=2)
      step = data[:, 0]

      valIdx, valName = eval(self.repFeature.get())
      values = data[:, valIdx]
      plt.plot(step, values, 'o-')
      plt.title(f'"{system.getSystemName()}" trajectory "{valName}"')
      plt.xlabel("Step")
      plt.ylabel(f"{valName}")
      plt.show()

    def showOpenMMDL(self, paramName=None):
      system = self.getMDSystem()
      anaDir = system.getOpenmmdlDir()

      option = self.openmmdlAnalysis.get()
      if option == 0:
        imgFile = os.path.join(anaDir, 'RMSD/RMSD_over_time.png')
        title = 'RMSD over time'
      elif option == 1:
        barType = self.getEnumText("barcodeType")
        imgFile = os.path.join(anaDir, f'Barcodes/{barType}_interactions.png')
        title = f'{barType} barcodes'
      elif option == 2:
        imgFile = os.path.join(anaDir, 'Binding_Modes_Markov_States/all_binding_modes_arranged.png')
        title = 'Binding_Modes_Markov_States'
      self.displayImage(imgFile, title=title)

    def showHeatmap(self, paramName=None):
      outLabel = 'Residue' if self.target.get() == 0 else 'Ligand'
      outDic = self.parseOpenMMDL(which=outLabel)

      system = self.getMDSystem()
      tTime = system.getNTime()
      self.makeHeatmapPlot(outDic, tTime, outLabel, self.threshold.get())

    def showHistogram(self, paramName=None):
      outLabel = 'Residue' if self.target.get() == 0 else 'Ligand'
      outDic = self.parseOpenMMDL(which=outLabel)

      self.makeHistogramPlot(outDic, self.threshold.get())

    def showInteractionsDiagram(self, paramName=None):
      outDic = self.parseOpenMMDL(which='interactions')
      self.makeInteractionsPlot(outDic)

    ############ UTILS FUNCTIONS ##################

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

    def displayImage(self, imageFile, title=''):
      import matplotlib.image as mpimg
      img = mpimg.imread(imageFile)
      plt.imshow(img)
      plt.axis('off')  # Hide axes
      plt.title(title)
      plt.show()

    def getBarcodeTypes(self):
      system = self.getMDSystem()
      anaDir = system.getOpenmmdlDir()
      if not os.path.exists(anaDir):
        return []
      else:
        anaDir = os.path.join(anaDir, 'Barcodes')
        types = []
        for file in os.listdir(anaDir):
          if '_interactions.png' in file:
            types.append(file.replace('_interactions.png', ''))
        return types

    def filterInteractions(self, outDic, nFrames):
      th = self.threshold.get()

      newDic = {}
      for pairId, intDic in outDic.items():
        for intType, frames in intDic.items():
          prop = len(frames) / nFrames
          if prop > th:
            if not pairId in newDic:
              newDic[pairId] = {}
            newDic[pairId][intType] = prop
      return newDic

    def parseOpenMMDL(self, which='residue'):
      '''Parse the OpenMMDL output and return a dictionary with the info of "which" group.
      Which in ['residue', 'ligand', 'interaction']
      '''
      which = which.lower()
      system = self.getMDSystem()
      anaDir = system.getOpenmmdlDir()

      inFile = os.path.join(anaDir, 'interactions_gathered.csv')
      nFrames = system.getNFrames()
      parseFuncs = {'residue': parseResidueLine, 'ligand': parseLigandLine, 'interactions': parseInteractionsLine}
      parseFunc = parseFuncs[which]

      outDic = {}
      with open(inFile) as f:
        csv_reader = csv.reader(f)
        next(csv_reader)
        for sline in csv_reader:
          outDic = parseFunc(sline, outDic, nFrames)

      if which == 'residue':
        outDic = cleanChain(outDic)

      elif which == 'interactions':
        outDic = cleanChain(outDic, pair=True)
        outDic = self.filterInteractions(outDic, nFrames)

      outDic = {k: outDic[k] for k in sorted(outDic)}
      return outDic


    def makeHeatmapPlot(self, d, totalTime, title='', th=0.1):
      df, filtIds = [], []
      for resId, frameDic in d.items():
        newContact = []
        for frame, ints in frameDic.items():
          newContact.append(len(ints))

        if sum(newContact) / len(newContact) > th:
          df.append(newContact)
          filtIds.append(resId)

      df = np.array(df)
      frames = list(range(1, len(df[0])+1))

      im, _ = heatmap(df, filtIds, frames, totalTime)
      if title:
        plt.savefig(f'{title}.png')
      plt.tight_layout()
      plt.show()

    def makeHistogramPlot(self, d, th=0.1):
      resContacts = {}
      for resId, frameDic in d.items():
        resDic, resCount = {}, 0
        for frame, ints in frameDic.items():
          for intType in ints:
            intType = intType.capitalize()
            if not intType in resDic:
              resDic[intType] = 0
            resDic[intType] += 1
            resCount += 1

        if resCount / len(frameDic) > th:
          resContacts[str(resId)] = resDic

      histogram(resContacts, len(frameDic))


    def makeInteractionsPlot(self, d):
      system = self.getMDSystem()
      topFile = system.getFileName()
      inDir = os.path.dirname(topFile)
      ligTopFile = os.path.abspath(os.path.join(inDir, f'{system.getSystemName()}_ligand.pdb'))
      subprocess.check_call(f'grep " {system.getLigandID()} " {topFile} > {ligTopFile}', shell=True)

      anaDir = system.getOpenmmdlDir()
      paramsFile = os.path.abspath(os.path.join(anaDir, 'drawInteractionsParams.txt'))
      outFile = os.path.abspath(os.path.join(anaDir, f'{system.getSystemName()}_interactions.png'))
      with open(paramsFile, 'w') as f:
        f.write(f'molFile :: {ligTopFile}\n')
        f.write(f'intDic :: {d}\n')
        f.write(f'outFile :: {outFile}\n')

      openmmPlugin.runScript(self, 'openmmDrawInteractions.py', args=paramsFile, env=RDKIT_DIC,
                             popen=True, cwd=self._getPath())




