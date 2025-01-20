# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
This module will prepare the system for the simulation
"""
import os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.utils import getBaseName, convertToSdf

from .. import Plugin
from ..constants import OPENMM_DIC
from ..objects import OpenMMSystem

STRUCTURE, LIGAND = 0, 1
GAFF_Vs = ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11']
SMIRNOFF_Vs = ['openff-1.0.1', 'openff-1.1.1', 'openff-1.0.0-RC1', 'openff-1.2.0', 'openff-1.1.0', 'openff-1.0.0', 'openff-1.0.0-RC2', 'smirnoff99Frosst-1.0.2', 'smirnoff99Frosst-1.0.0', 'smirnoff99Frosst-1.1.0', 'smirnoff99Frosst-1.0.4', 'smirnoff99Frosst-1.0.8', 'smirnoff99Frosst-1.0.6', 'smirnoff99Frosst-1.0.3', 'smirnoff99Frosst-1.0.1', 'smirnoff99Frosst-1.0.5', 'smirnoff99Frosst-1.0.9', 'smirnoff99Frosst-1.0.7']
SMIRNOFF_Vs.sort()
ESPALOMA_Vs = ['espaloma-0.3.2']

class ProtOpenMMSystemPrep(EMProtocol):
    """
    This protocol will start a Molecular Dynamics preparation. It will create the system
    and the topology, structure, and position restriction files

    It is necessary to insert a cleaned PDB structure from Protocol Import Atomic Structure
    or other similar protocols.
    """
    _label = 'system preparation'
    _cations = ['Cs+', 'K+', 'Li+', 'Na+', 'Rb+']
    _anions = ['Cl-', 'Br-', 'F-', 'I-']

    # -------------------------- DEFINE constants ----------------------------


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputFrom', params.EnumParam, default=STRUCTURE,
                      label='Input from: ', choices=['AtomStruct', 'SetOfSmallMolecules'],
                      help='Type of input you want to use')
        form.addParam('inputStructure', params.PointerParam, pointerClass='SchrodingerAtomStruct, AtomStruct',
                      label='Input structure to be prepared for MD:', allowsNull=False, condition='inputFrom==0',
                      help='Atomic structure to be prepared for MD by solvation, ions addition etc')
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Input set of molecules:', allowsNull=False, condition='inputFrom==1',
                      help='Input set of docked molecules. One of them will be prepared together with its target')
        form.addParam('inputLigand', params.StringParam, condition='inputFrom==1',
                      label='Ligand to prepare: ',
                      help='Specific ligand to prepare in the system')

        ffGroup = form.addGroup('System force fields')
        ffGroup.addParam('ffType', params.EnumParam, default=0, choices=['Amber14', 'CHARMM36', 'Old'],
                         label="Main atomic force field: ", help='Main force field to use')
        ffGroup.addParam('ffAmberType', params.EnumParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                         condition='ffType==0', label="Amber atomic force field: ",
                         choices=['All', 'protein.ff14SB', 'protein.ff15ipq', 'DNA.OL15', 'DNA.bsc1', 'RNA.OL3', 'lipid17'],
                         help='Amber main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')
        ffGroup.addParam('ffAmberWaterType', params.EnumParam, default=3, condition='ffType==0', label="Amber water force field: ",
                         choices=['SPCE', 'OPC', 'OPC3', 'tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb'],
                         help='Water amber force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')

        ffGroup.addParam('ffCHARMMWaterType', params.EnumParam, default=0, condition='ffType==1',
                         label="CHARMM water force field: ", expertLevel=params.LEVEL_ADVANCED,
                         choices=['Water', 'SPCE', 'tip3p-pme-b', 'tip3p-pme-f', 'tip4pew', 'tip4p2005', 'tip5p', 'tip5pew'],
                         help='Water CHARMM force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#charmm36')

        # ffGroup.addParam('ffAMOEBAType', params.EnumParam, default=0, expertLevel=params.LEVEL_ADVANCED,
        #                  choices=['2018', '2013', '2009'], condition='ffType==2', label="AMOEBA atomic force field: ",
        #                  help='AMOEBA main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amoeba')
        # ffGroup.addParam('useAMOEBAImplicit', params.BooleanParam, default=False,
        #                  label='Use AMOEBA implicit solvent: ', condition='ffType==2',
        #                  help='Whether to use the implicit or explicit AMOEBA solvent model')

        ffGroup.addParam('ffOldType', params.EnumParam, default=0,
                         choices=['amber96', 'amber99sb', 'amber99sbildn', 'amber99sbnmr', 'amber03', 'amber10', 'charmm_polar_2013'],
                         condition='ffType==2', label="Older force field: ",
                         help='Select an older main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#older-force-fields')
        # ffGroup.addParam('useOldImplicit', params.BooleanParam, default=False,
        #                  label='Use implicit solvent: ', condition='ffType==3',
        #                  help='Whether to use the implicit or explicit solvent model for Amber old force fields')

        ffGroup.addParam('ffWaterType', params.EnumParam, default=0,
                         choices=['tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb', 'tip5p', 'spce', 'swm4ndp', 'opc', 'opc3'],
                         condition='ffType==2', label="Water force field: ",
                         help='Select an water force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#water-models')

        ffGroup.addParam('ffSmallType', params.EnumParam, default=2, choices=['GAFF', 'SMIRNOFF', 'ESPALOMA'],
                         condition='inputFrom==1', label="Small molecules force field: ",
                         help='Small molecules force field to use')
        ffGroup.addParam('gaffVersion', params.EnumParam, default=4, choices=GAFF_Vs, expertLevel=params.LEVEL_ADVANCED,
                         condition='inputFrom==1 and ffSmallType==0', label="GAFF force field: ",
                         help='GAFF force field to use')
        ffGroup.addParam('smirnoffVersion', params.EnumParam, default=6, choices=SMIRNOFF_Vs,
                         expertLevel=params.LEVEL_ADVANCED,
                         condition='inputFrom==1 and ffSmallType==1', label="SMIRNOFF force field: ",
                         help='SMIRNOFF force field to use')
        ffGroup.addParam('espalomaVersion', params.EnumParam, default=0, choices=ESPALOMA_Vs,
                         expertLevel=params.LEVEL_ADVANCED,
                         condition='inputFrom==1 and ffSmallType==2', label="ESPALOMA force field: ",
                         help='ESPALOMA force field to use')

        ffGroup.addParam('constraints', params.EnumParam, default=1, label="Forcefield constraints: ",
                         choices=['None', 'HBonds', 'AllBonds', 'HAngles'],
                         help='http://docs.openmm.org/latest/userguide/application/02_running_sims.html#constraints')

        ffGroup = form.addGroup('Non bonded interactions')
        ffGroup.addParam('nonbondedMethod', params.EnumParam, default=0,
                         choices=['NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'Ewald', 'PME', 'LJPME'],
                         label="Non bonded method: ",
                         help='Non bonded method to simulate the non bonded atom interactions')
        ffGroup.addParam('nonbondedCutoff', params.FloatParam, default=1.0, expertLevel=params.LEVEL_ADVANCED,
                         label='Distance cutoff for non bonded interactions (nm): ', condition='nonbondedMethod!=0',
                         help='TThe cutoff distance to use for nonbonded interactions')

        ffGroup = form.addGroup('Hydrogens')
        ffGroup.addParam('addH', params.BooleanParam, default=False,
                         label='Add hydrogens to the system: ', help='Add hydrogens to the system')
        ffGroup.addParam('hPH', params.FloatParam, default=7.0, expertLevel=params.LEVEL_ADVANCED,
                         label='PH for hydrogen addition: ', help='The pH based on which to select variants')
        # todo: allow the use of variants

        form.addSection(label='Solvent box')
        sGroup = form.addGroup('Boundary box')
        sGroup.addParam('sizeType', params.EnumParam, label="System size type: ", default=1,
                        choices=['Absolute', 'Padding'], display=params.EnumParam.DISPLAY_HLIST,
                        help='Absolute: absolute size of the box (diameter)\n'
                             'Buffer: distance from the solute to the edge of the box\n')
        line = sGroup.addLine('Box size (nm):', condition='sizeType == 0',
                              help='Distances of the bounding box (nm).\nIf BSS, then it will be the value of the '
                                   'image distance')
        line.addParam('distA', params.FloatParam, default=5.0, label='a: ')
        line.addParam('distB', params.FloatParam, default=5.0, label='b: ')
        line.addParam('distC', params.FloatParam, default=5.0, label='c: ')
        sGroup.addParam('padDist', params.FloatParam, condition='sizeType == 1',
                        default=1.0, label='Padding distance: ',
                        help='Distance (nm) from the solute to the edge of the box.')

        iGroup = form.addGroup('Ions')
        iGroup.addParam('saltConc', params.FloatParam, default=0, label='Salt concentration (M): ',
                        help='Ionic strength to prepare the system')

        iGroup.addParam('neutralize', params.BooleanParam, default=True, label='Neutralize system: ',
                        help='Whether to add ions to the system until neutralize.')

        iGroup.addParam('cationType', params.EnumParam,
                      label='Cation to add: ', choices=self._cations, default=3,
                      help='Which cation to add in the system')

        iGroup.addParam('anionType', params.EnumParam,
                      label='Anions to add: ', choices=self._anions, default=0,
                      help='Which anion to add in the system')

    def _insertAllSteps(self):
      self._insertFunctionStep('solvateStep')
      self._insertFunctionStep('createOutputStep')


    def solvateStep(self):
      recFile = self.getReceptorFilename()
      molFile = self.getSpecifiedMolFile() if self.inputFrom.get() == LIGAND else None

      with open(self.getParamsFile(), 'w') as f:
        f.write(f'receptorFile :: {recFile}\n')
        if molFile:
          f.write(f'ligandFile :: {molFile}\n')
          f.write(f'ligandFF :: {self.getLigandFFVersion()}\n')

        mFF, wFF = self.getFFFiles()
        f.write(f'mFF :: {mFF}\nwFF :: {wFF}\n')
        f.write(f'nonbondedMethod :: {self.getEnumText("nonbondedMethod")}\n')
        f.write(f'nonbondedCutoff :: {self.nonbondedCutoff.get()}\n')
        f.write(f'constraints :: {self.getEnumText("constraints")}\n')

        wModel = self.getWaterModel(wFF)
        f.write(f'wModel :: {wModel}\n')

        f.write(f'addH :: {self.addH.get()}\n')
        if self.addH.get():
          f.write(f'hPH :: {self.hPH.get()}\n')

        if self.sizeType.get() == 0:
          f.write(f'boxSize :: {self.distA.get()}, {self.distB.get()}, {self.distC.get()}\n')
        else:
          f.write(f'padDist :: {self.padDist.get()}\n')

        f.write(f'saltConc :: {self.saltConc.get()}\n')
        f.write(f'neutralize :: {self.neutralize.get()}\n')
        f.write(f'cationType :: {self.getEnumText("cationType")}\n')
        f.write(f'anionType :: {self.getEnumText("anionType")}\n')

      Plugin.runScript(self, 'openmmPrepareSystem.py', args=self.getParamsFile(), env=OPENMM_DIC,
                             cwd=self._getPath())


    def createOutputStep(self):
      systemBasename = self.getSystemName()
      outStructFile, outSystemFile = self._getPath(f'{systemBasename}_system.pdb'), \
                                     self._getPath(f'{systemBasename}_system.xml')

      ligName = self.inputLigand.get() if self.inputFrom.get() == LIGAND else None
      mFF, wFF = self.getFFFiles()
      outSystem = OpenMMSystem(filename=outStructFile, oriStructFile=outStructFile, serieFile=outSystemFile,
                               ff=mFF, wff=wFF, ligName=ligName)

      self._defineOutputs(outputSystem=outSystem)
      # self._defineSourceRelation(self.inputStructure, outSystem)


    def getWaterModel(self, wFF):
      model = 'tip3p'
      if 'spce' in wFF:
        model = 'spce'
      elif 'tip4p' in wFF:
        model = 'tip4pew'
      elif 'tip5p' in wFF:
        model = 'tip5p'
      return model

    def getFFFiles(self):
      if self.ffType.get() == 0:
        mFF = 'amber14-all.xml' if self.ffAmberType.get() == 0 \
          else 'amber14/{}.xml'.format(self.getEnumText('ffAmberType'))
        wFF = 'amber14/{}.xml'.format(self.getEnumText('ffAmberWaterType').lower())

      elif self.ffType.get() == 1:
        mFF = 'charmm36.xml'
        wFF = 'charmm36/{}.xml'.format(self.getEnumText('ffCHARMMWaterType').lower())

      # elif self.ffType.get() == 2:
      #   mFF = 'amoeba{}.xml'.format(self.getEnumText('ffAMOEBAType'))
      #   # wFF = '{}.xml'.format(self.getEnumText('ffWaterType'))
      #   ffs = [mFF]

      elif self.ffType.get() == 2:
        mFF = '{}.xml'.format(self.getEnumText('ffOldType'))
        wFF = '{}.xml'.format(self.getEnumText('ffWaterType'))

      return mFF, wFF

    def getLigandFFVersion(self):
      ffOption = self.ffSmallType.get()
      if ffOption == 0:
        return self.getEnumText('gaffVersion')
      elif ffOption == 1:
        return self.getEnumText('smirnoffVersion')
      elif ffOption == 2:
        return self.getEnumText('espalomaVersion')

    def getParamsFile(self):
      return os.path.abspath(self._getExtraPath('solvationParams.txt'))

    def getReceptorFilename(self):
      if self.inputFrom.get() == STRUCTURE:
          proteinFile = self.inputStructure.get().getFileName()
      elif self.inputFrom.get() == LIGAND:
          proteinFile = self.inputSetOfMols.get().getProteinFile()
      return os.path.abspath(proteinFile)

    def getSystemName(self):
      return getBaseName(self.getReceptorFilename())

    def getSpecifiedMolFile(self):
        myMol = None
        for mol in self.inputSetOfMols.get():
          if mol.__str__() == self.inputLigand.get():
            myMol = mol.clone()
            break
        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            molFile = myMol.getPoseFile()
            return convertToSdf(self, molFile)
