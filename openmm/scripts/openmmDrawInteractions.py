#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Daniel Del Hoyo Gómez (ddelhoyo@cnb.csic.es)
# # # *
# # # *
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import math, sys
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.figure import Figure

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from utils import parseParams, getMolFilesDic

RES_COLORS = {'Apolar': 'green', 'Negative': 'blue', 'Positive': 'red', 'Polar': 'cyan', 'Aromatic': 'purple'}
RES_GROUPS = {'Apolar': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET'], 'Aromatic': ['PHE', 'TYR', 'TRP'],
							'Polar': ['SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN'], 'Negative': ['ASP', 'GLU'],
							'Positive': ['LYS', 'ARG', 'HIS']}

COLOR_DIC = {res: RES_COLORS[group] for group, resGroup in RES_GROUPS.items() for res in resGroup}

COLOR_INT_DIC = {'halogen': 'red', 'hydrophobic': 'green', 'pication': 'orange', 'pistacking': 'brown',
								 'saltbridge': 'gray', 'waterbridge': 'blue', 'hbond': 'cyan'}

DEF_SAVE_FILE = 'moleculeInteractions.png'


def getExternalPosition(atomPoint, molCenter, distance=80):
	"""Returns a position afar from the molecule to draw a residue circle"""
	dx = atomPoint[0] - molCenter[0]
	dy = atomPoint[1] - molCenter[1]

	length = math.sqrt(dx ** 2 + dy ** 2)
	if length > 0:
		dx /= length
		dy /= length

	newX = atomPoint[0] + dx * distance
	newY = atomPoint[1] + dy * distance

	return newX, newY


def getMoleculeCenter(drawer, mol):
	"""Get the center of the drawn molecules"""
	coordsX, coordsY = [], []

	for i in range(mol.GetNumAtoms()):
		try:
			atomPoint = drawer.GetDrawCoords(i)
			coordsX.append(atomPoint.x)
			coordsY.append(atomPoint.y)
		except:
			continue

	if coordsX and coordsY:
		centroX, centroY = np.mean(coordsX), np.mean(coordsY)
		return centroX, centroY
	else:
		return None


def setConnectionLineCoords(c1, c2, displac, margin=50):
	'''Draw the connection lines based on the extreme coordinates, with some margin. 
	Also, the lines can be perpendiculary displaced to accommodate several interactions that do not overlap'''
	x1, x2 = c1[0], c2[0]
	y1, y2 = c1[1], c2[1]
	dx, dy = x2 - x1, y2 - y1
	length = np.hypot(dx, dy)
	displac = displac * 0.01	

	ux, uy = dx / length, dy / length
	x1New, x2New = x1 + ux * margin/2 + displac * dx, x2 - ux * margin + displac * dx
	y1New, y2New = y1 + uy * margin/2 - displac * dy, y2 - uy * margin - displac * dy
	return (x1New, y1New), (x2New, y2New)


def getBoundingBox(allPoints, margin=200):
	'''Get the bounding box for the figure given all the points drawn, with a given margin to accommodate them'''
	allX = [p[0] for p in allPoints]
	allY = [p[1] for p in allPoints]

	minX, maxX = min(allX), max(allX)
	minY, maxY = min(allY), max(allY)

	bboxX = [minX - margin, maxX + margin]
	bboxY = [minY - margin, maxY + margin]
	return bboxX, bboxY

def parsePDBAtomNumbers(ligFile):
		'''Parse the PDB to transform PDB atom numbers in RDKit atom numbers'''
		d = {}
		with open(ligFile) as f:
				for i, line in enumerate(f):
						d[int(line[7:11])] = i
		return d

def removeNumbers(s):
	'''Remove digits form string'''
	return ''.join(c for c in s if not c.isdigit())


def getCirclesRepulsion(circle, allCircles, interactDistance, molCenter, kRep):
		repulsionForce = np.zeros(2)
		for otherCircle in allCircles:
			if circle != otherCircle:
				direction = circle - otherCircle
				distance = np.linalg.norm(direction)
				if distance < interactDistance:
					if distance == 0:
						direction = circle - molCenter
						direction[0] = -direction[0]
						distance = 100

					direction /= distance
					repF = direction * kRep / (distance)
					repulsionForce += repF
		return repulsionForce

def getAtomsRepulsion(circle, atomPositions, atomInteractRadius, kRep):
	repulsionForce = np.zeros(2)
	for atomPos in atomPositions:
		direction = circle - atomPos
		distance = np.linalg.norm(direction)
		if distance < atomInteractRadius:
			direction /= distance
			repForce = direction * kRep / (distance)
			repulsionForce += repForce
	return repulsionForce

def optimizeCirclePositions(circlePositions, atomPositions, molCenter,
														interactDistance, atomInteractRadius, iterations=100):
	"""Optimize the initial positions of the resdiue circles to avoid overlaps using repulsion forces between circles and 
  other circles and atoms. 
  Originally, the circle positions are next to the interactions atoms.
  """
	circles = np.array(circlePositions)
	kRep = 1000

	for _ in range(iterations):
		repForces = []
		for i, circle in enumerate(circles):
			# Repulsion against other circles
			repulsionForce = getCirclesRepulsion(circle, circles, interactDistance, molCenter, kRep)
			# Repulsion against closeby atoms
			repulsionForce += getAtomsRepulsion(circle, atomPositions, atomInteractRadius, kRep)

			# Update position
			circles[i] += repulsionForce
			repForces.append(sum(repulsionForce))

		# If the positions are not updated much anymore, exit
		if max(repForces) < 0.1:
			break
	return circles

def getDisplacements(n):
	a = list(range(n))
	meanA = sum(a) / len(a)
	return [i-meanA for i in a]


def normalizeAngle(angle, ofx, ofy):
	"""
  Normalize the angle in range[-90, 90] degrees so the text is always upside down
  """
	angle = angle % 360
	if angle > 180:
		angle -= 360

	if angle > 90:
		angle -= 180
		ofx, ofy = -ofx, -ofy
	elif angle < -90:
		angle += 180
		ofx, ofy = -ofx, -ofy

	return angle, ofx, ofy


def addParallelLabel(ax, x, y, texto, color, distancia=5, offSign=1, **kwargs):
	"""
  Add a label parallel to the interaction line, with the desired color at a desired perpedicular distance to the line.
  Offsign can be used to control whether the label will be placed in one or the other side of the line

  """
	dx = x[-1] - x[0]
	dy = y[-1] - y[0]
	angle = np.degrees(np.arctan2(dy, dx))

	length = np.sqrt(dx ** 2 + dy ** 2)
	if length > 0:
		ux = -dy / length  
		uy = dx / length  
	else:
		ux, uy = 0, 1

	midX = np.mean(x)
	midY = np.mean(y)
	offsetX = distancia * ux * offSign
	offsetY = distancia * uy * offSign
	angle, offsetX, offsetY = normalizeAngle(angle, offsetX, offsetY)

	label = ax.annotate(texto, xy=(midX, midY), xytext=(offsetX, offsetY), color=color, textcoords='offset points',
							rotation=angle, ha='center', va='center', **kwargs)
	return label


class MoleculeInteractions:
	'''Main object to draw a molecule and its interacting residues, obtained from OpenMMDL.
	The residues are placed surrounding the moelcule, next to their interacting atoms, but their positions can be 
	dynamically changed.
	'''
	def __init__(self, molFile, inDic, outFile):
		self.canvasSize = 2000
		self.circleRadius = 40
		self.interactDistance = 200
		self.outFile = outFile
		
		self.circles, self.resCircles = {}, {}
		self.cirLabels = {}
		self.selectedCircle = None
		
		molsDict, _ = getMolFilesDic([molFile], True)
		self.mol = list(molsDict.keys())[0]
		self.atomNumDic = parsePDBAtomNumbers(molFile)

		self.nDic = self.processInputDic(inDic)
		self.drawMolecule()
		self.cirDic, allPoints = self.getDrawnPoints(self.nDic)

		self.getBoundingBoxes(allPoints)
		self.drawAminoacids()
		self.drawInteractions()
		self.addLegends()

	########## Actions tools ########################
	def onPress(self, event):
		if event.inaxes != self.ax:
			return

		for circle in self.circles:
			contains, _ = circle.contains(event)
			if contains:
				self.selectedCircle = circle
				break

	def onRelease(self, event):
		self.selectedCircle = None

	def onMotion(self, event):
		if self.selectedCircle is None or event.inaxes != self.ax:
			return

		#  Moving residue Circle and label
		self.selectedCircle.center = (event.xdata, event.ydata)
		self.cirLabels[self.selectedCircle].set_position((event.xdata, event.ydata))

		# Redoing interactions in new position
		self.removeCircleInteractions(self.selectedCircle)
		self.drawCircleInteractions(self.selectedCircle)

		self.fig.canvas.draw()

	def on_close(self, event):
		event.canvas.figure.savefig(self.outFile, dpi=300, bbox_inches='tight')
		plt.close('all')

	########## Main tools ########################
	def drawMolecule(self):
			'''Draw the moelcule using RDKit'''
			rdDepictor.Compute2DCoords(self.mol)
			Chem.RemoveStereochemistry(self.mol)
			self.mol = rdMolDraw2D.PrepareMolForDrawing(self.mol)

			self.drawer = rdMolDraw2D.MolDraw2DCairo(int(self.canvasSize * 0.5), int(self.canvasSize * 0.5))
			self.drawer.DrawMolecule(self.mol)
			self.drawer.FinishDrawing()

	def getDrawnPoints(self, nDic):
		'''Places the residues given the molecule atoms positions and returns that information
		'''
		molPoints = self.getMolPoints()
		cirDic = self.getCirclesDic(nDic)

		centerX, centerY = getMoleculeCenter(self.drawer, self.mol)
		newCircles = optimizeCirclePositions(list(cirDic.values()), molPoints, [centerX, centerY],
																				 self.interactDistance, self.interactDistance)

		for i, resId in enumerate(cirDic):
			newCircle = newCircles[i]
			molPoints.append(newCircle)
			cirDic[resId] = newCircle

		return cirDic, molPoints

	def getBoundingBoxes(self, allPoints):
			'''Sets the bounding boxes based on the drawn positions'''
			bboxX, bboxY = getBoundingBox(allPoints)
			bboxWidth, bboxHeight = bboxX[1] - bboxX[0], bboxY[1] - bboxY[0]

			self.drawer.WriteDrawingText('mol_temp_full.png')
			imgFull = mpimg.imread('mol_temp_full.png')

			self.fig, self.ax = plt.subplots(figsize=(10, 10 * bboxHeight / bboxWidth))

			self.ax.imshow(imgFull)
			self.ax.set_xlim(bboxX[0], bboxX[1])
			self.ax.set_ylim(bboxY[0], bboxY[1])
			self.ax.set_axis_off()

			self.fig.canvas.mpl_connect('button_press_event', self.onPress)
			self.fig.canvas.mpl_connect('button_release_event', self.onRelease)
			self.fig.canvas.mpl_connect('motion_notify_event', self.onMotion)
			self.fig.canvas.mpl_connect('close_event', self.on_close)
	
	def drawAminoacids(self):
		'''Draw all the interacting resiudes and their labels'''
		for resId, (circleX, circleY) in self.cirDic.items():
			circleCoords = (circleX, circleY)
			circle = Circle(circleCoords, self.circleRadius, facecolor=COLOR_DIC[removeNumbers(resId.upper())],
											linewidth=2, alpha=0.2, transform=self.ax.transData)
			self.ax.add_patch(circle)
			cirLabel = self.ax.text(circleX, circleY, resId, color='black', weight='bold', fontsize=8,
															ha='center', va='center', transform=self.ax.transData)

			self.resCircles[resId] = circle
			self.circles[circle] = resId
			self.cirLabels[circle] = cirLabel

	def drawInteractions(self):
		'''Draw all the interactions between residues and moelcule atoms'''
		self.circleInteractions = {}
		for resId in self.cirDic:
			circle = self.resCircles[resId]
			self.drawCircleInteractions(circle)

	def addLegends(self):
		'''Add residues and interaction legends to the graph'''
		residuePatches = []
		for resGroup, color in RES_COLORS.items():
			element = Patch(facecolor=color, alpha=0.4, edgecolor='black', label=resGroup, linewidth=1)
			residuePatches.append(element)

		residueLegend = self.ax.legend(handles=residuePatches, title="Residues type", frameon=True,
																		loc='upper left')
		self.ax.add_artist(residueLegend)

		interactionPatches = []
		for intType, color in COLOR_INT_DIC.items():
			element = Line2D([0], [0], color=color, linewidth=2.5,
												linestyle='--', label=intType)
			interactionPatches.append(element)

		interactionLegend = self.ax.legend(handles=interactionPatches, title="Interactions type", frameon=True,
																		 loc='upper right')
		self.ax.add_artist(interactionLegend)


	def display(self):
		'''Main function to display the generated graph'''
		plt.tight_layout(pad=5)
		plt.show()


	########## HELP tools ########################
	def processInputDic(self, intDic):
		'''Processes the input dictionary containing the interacting information'''
		nDic = {}
		for pairId, iDic in intDic.items():
			resId, atomIds = pairId[0].replace('_', ''), pairId[1]
			if '_' in atomIds:
				nAtomIds = '_'.join([str(self.atomNumDic[int(atomId)]) for atomId in atomIds.split('_')])
			else:
				nAtomIds = self.atomNumDic[int(eval(atomIds))]

			if resId not in nDic:
				nDic[resId] = {}
			nDic[resId][nAtomIds] = {}
			for intType, n in iDic.items():
				nDic[resId][nAtomIds][intType] = n
		return nDic

	def getCirclesDic(self, nDic):
		'''Gets the original positions of the circles, close to their interacting atoms'''
		cirDic = {}
		centerX, centerY = getMoleculeCenter(self.drawer, self.mol)
		for resId, atomDic in nDic.items():
			atomIdx = list(atomDic.keys())[0]
			if isinstance(atomIdx, str):
				atomIdx = int(atomIdx.split('_')[0])
				iDist = self.interactDistance * 2
			else:
				iDist = self.interactDistance

			atomPoint = self.drawer.GetDrawCoords(atomIdx)
			atomPoint = [atomPoint.x, atomPoint.y]
			circleX, circleY = getExternalPosition(atomPoint, (centerX, centerY), iDist)
			cirDic[resId] = [circleX, circleY]
		return cirDic
	
	def getMolPoints(self):
		'''Returns the molecule atom coordinates'''
		allPoints = []
		for i in range(self.mol.GetNumAtoms()):
			atomPoint = self.drawer.GetDrawCoords(i)
			allPoints.append((atomPoint.x, atomPoint.y))
		return allPoints

	def addInteractLine(self, atomCoords, circleCoords, disp, intColor, n, offSign):
		'''Adds a interaction line and its label given the circle and atom coordinates
		'''
		c1, c2 = setConnectionLineCoords(atomCoords, circleCoords, disp)
		xs, ys = [c1[0], c2[0]], [c1[1], c2[1]]
		line = Line2D(xs, ys, color=intColor, alpha=0.6, linewidth=1.5, linestyle='--')
		self.ax.add_line(line)
		label = addParallelLabel(self.ax, xs, ys, n, intColor, offSign=offSign)
		return line, label

	def drawCircleInteractions(self, circle):
		'''Functions that draws all the interactions of a given circle'''
		resId = self.circles[circle]
		circleCoords = circle.get_center()

		circleLines = []
		for atomIdx, intsDic in self.nDic[resId].items():
			if isinstance(atomIdx, str):
				piCoords = []
				for atomAro in atomIdx.split('_'):
					atomPoint = self.drawer.GetDrawCoords(int(atomAro))
					piCoords.append([atomPoint.x, atomPoint.y])

				piPoint = np.mean(np.array(piCoords), axis=0)
				atomCoords = (piPoint[0], piPoint[1])

			else:
				atomPoint = self.drawer.GetDrawCoords(atomIdx)
				atomCoords = (atomPoint.x, atomPoint.y)

			nInts = len(intsDic)
			displacs = getDisplacements(nInts)
			for i, (intType, n) in enumerate(intsDic.items()):
				offSign = 1 if i == 0 else -1
				intColor = COLOR_INT_DIC[intType.lower()]
				linea, label = self.addInteractLine(atomCoords, circleCoords, displacs[i], intColor, n, offSign)
				circleLines.append([linea, label])

		self.circleInteractions[circle] = circleLines

	def removeCircleInteractions(self, circle):
		'''Remove the existing interactions of a circle'''
		for line, label in self.circleInteractions[circle]:
			line.remove(), label.remove()


if __name__ == "__main__":
		pDic = parseParams(sys.argv[1], sep='::')
		molFile, intDic = pDic['molFile'], eval(pDic['intDic'])
		outFile = pDic['outFile'] if 'outFile' in pDic else DEF_SAVE_FILE

		app = MoleculeInteractions(molFile, intDic, outFile)
		app.display()

