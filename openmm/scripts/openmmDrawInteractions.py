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

import math, sys, random
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.image as mpimg

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from utils import parseParams, getMolFilesDic

RES_COLORS = {'Apolar': 'green', 'Negative': 'blue', 'Positive': 'red', 'Polar': 'cyan', 'Aromatic': 'purple'}
RES_GROUPS = {'Apolar': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET'], 'Aromatic': ['PHE', 'TYR', 'TRP'],
							'Polar': ['SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN'], 'Negative': ['ASP', 'GLU'],
							'Positive': ['LYS', 'ARG', 'HIS']}

COLOR_DIC = {res: RES_COLORS[group] for group, resGroup in RES_GROUPS.items() for res in resGroup}

COLOR_INT_DIC = {'halogen': 'red', 'hydrophobic': 'green', 'pication': 'orange', 'pistacking': 'brown', 'saltbridge': 'gray'}


def getExternalPosition(atomPoint, molCenter, distance=80):
	"""Returns a position afar from the molecule to draw a residue circle"""
	dx = atomPoint.x - molCenter[0]
	dy = atomPoint.y - molCenter[1]

	length = math.sqrt(dx ** 2 + dy ** 2)
	if length > 0:
		dx /= length
		dy /= length

	new_x = atomPoint.x + dx * distance
	new_y = atomPoint.y + dy * distance

	return new_x, new_y


def getMoleculeCenter(drawer, mol):
	"""Get the center of the drawn molecules"""
	coords_x = []
	coords_y = []

	for i in range(mol.GetNumAtoms()):
		try:
			atom_point = drawer.GetDrawCoords(i)
			coords_x.append(atom_point.x)
			coords_y.append(atom_point.y)
		except:
			continue

	if coords_x and coords_y:
		centro_x = np.mean(coords_x)
		centro_y = np.mean(coords_y)
		return centro_x, centro_y
	else:
		return None


def setConnectionLineCoords(c1, c2, displac, margin=50):
	x1, x2 = c1[0], c2[0]
	y1, y2 = c1[1], c2[1]
	dx, dy = x2 - x1, y2 - y1
	length = np.hypot(dx, dy)
	displac = displac * 0.02

	ux, uy = dx / length, dy / length
	x1_new, x2_new = x1 + ux * margin/2 + displac * dx, x2 - ux * margin + displac * dx
	y1_new, y2_new = y1 + uy * margin/2 - displac * dy, y2 - uy * margin - displac * dy
	return (x1_new, y1_new), (x2_new, y2_new)


def getBoundingBox(allPoints, margin=200):
	all_x = [p[0] for p in allPoints]
	all_y = [p[1] for p in allPoints]

	min_x, max_x = min(all_x), max(all_x)
	min_y, max_y = min(all_y), max(all_y)

	bbox_x = [min_x - margin, max_x + margin]
	bbox_y = [min_y - margin, max_y + margin]
	return bbox_x, bbox_y

def parsePDBAtomNumbers(ligFile):
		d = {}
		with open(ligFile) as f:
				for i, line in enumerate(f):
						d[int(line[7:11])] = i
		return d

def remove_numbers(s):
	return ''.join(c for c in s if not c.isdigit())


def optimize_circle_positions(circle_positions, atomPositions, molCenter, atomInteractRadius, iterations=100):
	"""
  Optimiza las posiciones de los círculos para minimizar solapamientos
  """
	circles = np.array(circle_positions)
	k_repulsion_c, k_repulsion_a = 1000, 1000  # Constante de repulsión

	for _ in range(iterations):
		repForces = []
		for i, circle in enumerate(circles):
			# Fuerza de repulsión con otros círculos
			repulsion_force = np.zeros(2)
			for j, other_circle in enumerate(circles):
				if i != j:
					direction = circle - other_circle
					distance = np.linalg.norm(direction)
					if distance < interactDistance:
						if distance == 0:
							direction = circle - molCenter
							direction[0] = -direction[0]
							distance = 100

						direction /= distance
						# Repulsión inversamente proporcional a la distancia
						repF = direction * k_repulsion_c / (distance)
						repulsion_force += repF

			for atomPos in atomPositions:
				direction = circle - atomPos
				distance = np.linalg.norm(direction)
				if distance < atomInteractRadius:
					direction /= distance
					# Repulsión inversamente proporcional a la distancia
					repForce = direction * k_repulsion_a / (distance)
					repulsion_force += repForce
					# print('distance: ', distance, atomInteractRadius, repForce)

			# Actualizar posición
			circles[i] += repulsion_force
			repForces.append(sum(repulsion_force))

		# If the positions are not updated much anymore, exit
		if max(repForces) < 0.1:
			break
	return circles

def getDisplacements(n):
	a = list(range(n))
	meanA = sum(a) / len(a)
	return [i-meanA for i in a]

if __name__ == "__main__":
		# Drawing parameters
		canvas_size = 2000
		circleRadius = 40
		interactDistance = 200

		pDic = parseParams(sys.argv[1], sep='::')
		sanitize = True
		molsDict, _ = getMolFilesDic([pDic['molFile']], sanitize)

		mol = list(molsDict.keys())[0]
		rdDepictor.Compute2DCoords(mol)
		Chem.RemoveStereochemistry(mol)
		mol = rdMolDraw2D.PrepareMolForDrawing(mol)

		drawer = rdMolDraw2D.MolDraw2DCairo(int(canvas_size * 0.5), int(canvas_size * 0.5))
		drawer.DrawMolecule(mol)
		drawer.FinishDrawing()

		intDic = eval(pDic['intDic'])
		atomNumDic = parsePDBAtomNumbers(pDic['molFile'])
		print('atomNumDic: ', atomNumDic)

		nDic = {}
		for pairId, iDic in intDic.items():
				resId, atomId = pairId[0].replace('_', ''), atomNumDic[pairId[1]]
				if not resId in nDic:
					nDic[resId] = {}
				nDic[resId][atomId] = {}
				for intType, n in iDic.items():
					nDic[resId][atomId][intType] = n

		allPoints = []
		for i in range(mol.GetNumAtoms()):
			atom_point = drawer.GetDrawCoords(i)
			allPoints.append((atom_point.x, atom_point.y))

		cirDic = {}
		centerX, centerY = getMoleculeCenter(drawer, mol)
		for resId, atomDic in nDic.items():
			atom_idx = list(atomDic.keys())[0]
			atom_point = drawer.GetDrawCoords(atom_idx)
			circle_x, circle_y = getExternalPosition(atom_point, (centerX, centerY), interactDistance)
			cirDic[resId] = [circle_x, circle_y]

		# todo: Aromatic ints connect to center of cycle. Type interactions
		# todo: dynamic???
		newCircles = optimize_circle_positions(list(cirDic.values()), allPoints, [centerX, centerY], interactDistance)

		for i, resId in enumerate(cirDic):
			newCircle = newCircles[i]
			allPoints.append(newCircle)
			cirDic[resId] = newCircle

		# 4. CALCULAR BOUNDING BOX DE TODOS LOS ELEMENTOS
		bbox_x, bbox_y = getBoundingBox(allPoints)
		bbox_width = bbox_x[1] - bbox_x[0]
		bbox_height = bbox_y[1] - bbox_y[0]

		# 5. GUARDAR IMAGEN COMPLETA Y LUEGO RECORTAR
		drawer.WriteDrawingText('mol_temp_full.png')
		img_full = mpimg.imread('mol_temp_full.png')

		# 6. CREAR FIGURA Y AJUSTAR VISTA AL BOUNDING BOX
		fig, ax = plt.subplots(figsize=(10, 10 * bbox_height / bbox_width))  # Proporciones correctas

		# Mostrar solo la región del bounding box
		ax.imshow(img_full)
		ax.set_xlim(bbox_x[0], bbox_x[1])
		ax.set_ylim(bbox_y[0], bbox_y[1])  # Recordar: eje Y invertido
		ax.set_axis_off()

		# 7. DIBUJAR CÍRCULOS Y CONECTORES (usando las coordenadas ya calculadas)
		for resId, (circle_x, circle_y) in cirDic.items():
			circleCoords = (circle_x, circle_y)
			# Dibujar círculo
			circle = Circle(circleCoords, circleRadius,
				facecolor=COLOR_DIC[remove_numbers(resId.upper())], linewidth=3, alpha=0.2, transform=ax.transData)
			ax.add_patch(circle)

			# Etiqueta
			ax.text(circle_x, circle_y, resId,
							color='black', weight='bold', fontsize=10, ha='center', va='center', transform=ax.transData)

			piStackAtoms, piCatAtoms = [], []
			for atom_idx, intsDic in nDic[resId].items():
					atom_point = drawer.GetDrawCoords(atom_idx)
					atomCoords = (atom_point.x, atom_point.y)

					nInts = len(intsDic)
					displacs = getDisplacements(nInts)
					for i, (intType, n) in enumerate(intsDic.items()):
							if intType.lower() == 'pistacking':
									piStackAtoms.append(atomCoords)
							elif intType.lower() == 'pication':
									piCatAtoms.append(atomCoords)
							else:
									(x1_new, y1_new), (x2_new, y2_new) = setConnectionLineCoords(atomCoords, circleCoords, displacs[i])
									ax.plot([x1_new, x2_new], [y1_new, y2_new], COLOR_INT_DIC[intType.lower()],
													linestyle='--', alpha=0.6, linewidth=1.5)


			if piStackAtoms:
					aroCenter = np.mean(np.array(piStackAtoms), axis=0)
					(x1_new, y1_new), (x2_new, y2_new) = setConnectionLineCoords(aroCenter, circleCoords, -0.5)
					ax.plot([x1_new, x2_new], [y1_new, y2_new], COLOR_INT_DIC['pistacking'],
									linestyle='--', alpha=0.6, linewidth=1.5)
			if piCatAtoms:
					aroCenter = np.mean(np.array(piCatAtoms), axis=0)
					(x1_new, y1_new), (x2_new, y2_new) = setConnectionLineCoords(aroCenter, circleCoords, 0.5)
					ax.plot([x1_new, x2_new], [y1_new, y2_new], COLOR_INT_DIC['pication'],
									linestyle='--', alpha=0.6, linewidth=1.5)


		# 8. GUARDAR RESULTADO FINAL
		plt.tight_layout(pad=0)
		plt.savefig('molecula_ajustada_automaticamente.png', dpi=300, bbox_inches='tight')
		plt.show()
