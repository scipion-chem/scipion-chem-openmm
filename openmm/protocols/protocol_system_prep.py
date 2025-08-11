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

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC
from pwchem.utils import getBaseName, convertToSdf

from .. import Plugin
from ..constants import OPENMM_DIC
from ..objects import OpenMMSystem

scriptLigPrepName = 'rdkit_addHydrogens.py'

STRUCTURE, LIGAND = 0, 1
GAFF_Vs = ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11']
SMIRNOFF_Vs = ['openff-1.0.1', 'openff-1.1.1', 'openff-1.0.0-RC1', 'openff-1.2.0', 'openff-1.1.0', 'openff-1.0.0', 'openff-1.0.0-RC2', 'smirnoff99Frosst-1.0.2', 'smirnoff99Frosst-1.0.0', 'smirnoff99Frosst-1.1.0', 'smirnoff99Frosst-1.0.4', 'smirnoff99Frosst-1.0.8', 'smirnoff99Frosst-1.0.6', 'smirnoff99Frosst-1.0.3', 'smirnoff99Frosst-1.0.1', 'smirnoff99Frosst-1.0.5', 'smirnoff99Frosst-1.0.9', 'smirnoff99Frosst-1.0.7']
SMIRNOFF_Vs.sort()
ESPALOMA_Vs = ['espaloma-0.3.2']

CATION_NAMES = ['Cs+', 'K+', 'Li+', 'Na+', 'Rb+']
ANION_NAMES = ['Cl-', 'Br-', 'F-', 'I-']

LIG_INPUT = f'inputFrom == {LIGAND}'

class ProtOpenMMSystemPrep(EMProtocol):
    """
    This protocol will start a Molecular Dynamics preparation. It will create the system
    and the topology, structure, and position restriction files

    It is necessary to insert a cleaned PDB structure from Protocol Import Atomic Structure
    or other similar protocols.
    """
    _label = 'system preparation'

    def _defineFFParams(self, form, ligandCondition='True'):
        form.addParam('ffType', params.EnumParam, default=0, choices=['Amber14', 'CHARMM36', 'Old'],
                      label="Main atomic force field: ", help='Main force field to use')
        form.addParam('ffAmberType', params.EnumParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                      condition='ffType==0', label="Amber atomic force field: ",
                      choices=['All', 'protein.ff14SB', 'protein.ff15ipq', 'DNA.OL15', 'DNA.bsc1', 'RNA.OL3', 'lipid17'],
                      help='Amber main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')
        form.addParam('ffAmberWaterType', params.EnumParam, default=3, condition='ffType==0',
                      label="Amber water force field: ",
                      choices=['SPCE', 'OPC', 'OPC3', 'tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb'],
                      help='Water amber force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amber14')

        form.addParam('ffCHARMMWaterType', params.EnumParam, default=0, condition='ffType==1',
                      label="CHARMM water force field: ", expertLevel=params.LEVEL_ADVANCED,
                      choices=['Water', 'SPCE', 'tip3p-pme-b', 'tip3p-pme-f', 'tip4pew', 'tip4p2005', 'tip5p', 'tip5pew'],
                      help='Water CHARMM force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#charmm36')

        # ffGroup.addParam('ffAMOEBAType', params.EnumParam, default=0, expertLevel=params.LEVEL_ADVANCED,
        #                  choices=['2018', '2013', '2009'], condition='ffType==2', label="AMOEBA atomic force field: ",
        #                  help='AMOEBA main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#amoeba')
        # ffGroup.addParam('useAMOEBAImplicit', params.BooleanParam, default=False,
        #                  label='Use AMOEBA implicit solvent: ', condition='ffType==2',
        #                  help='Whether to use the implicit or explicit AMOEBA solvent model')

        form.addParam('ffOldType', params.EnumParam, default=0,
                      choices=['amber96', 'amber99sb', 'amber99sbildn', 'amber99sbnmr', 'amber03', 'amber10', 'charmm_polar_2013'],
                      condition='ffType==2', label="Older force field: ",
                      help='Select an older main force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#older-force-fields')
        # ffGroup.addParam('useOldImplicit', params.BooleanParam, default=False,
        #                  label='Use implicit solvent: ', condition='ffType==3',
        #                  help='Whether to use the implicit or explicit solvent model for Amber old force fields')

        form.addParam('ffWaterType', params.EnumParam, default=0,
                      choices=['tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb', 'tip5p', 'spce', 'swm4ndp', 'opc', 'opc3'],
                      condition='ffType==2', label="Water force field: ",
                      help='Select an water force field to use. http://docs.openmm.org/latest/userguide/application/02_running_sims.html#water-models')

        form.addParam('ffSmallType', params.EnumParam, default=2, choices=['GAFF', 'SMIRNOFF', 'ESPALOMA'],
                      condition=ligandCondition, label="Small molecules force field: ",
                      help='Small molecules force field to use')
        form.addParam('gaffVersion', params.EnumParam, default=4, choices=GAFF_Vs, expertLevel=params.LEVEL_ADVANCED,
                      condition=f'{ligandCondition} and ffSmallType==0', label="GAFF force field: ",
                      help='GAFF force field to use')
        form.addParam('smirnoffVersion', params.EnumParam, default=6, choices=SMIRNOFF_Vs,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition=f'{ligandCondition} and ffSmallType==1', label="SMIRNOFF force field: ",
                      help='SMIRNOFF force field to use')
        form.addParam('espalomaVersion', params.EnumParam, default=0, choices=ESPALOMA_Vs,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition=f'{ligandCondition} and ffSmallType==2', label="ESPALOMA force field: ",
                      help='ESPALOMA force field to use')

        form.addParam('constraints', params.EnumParam, default=1, label="Forcefield constraints: ",
                      choices=['None', 'HBonds', 'AllBonds', 'HAngles'],
                      help='You can optionally tell OpenMM to constrain certain bond lengths and angles.'
                           'http://docs.openmm.org/latest/userguide/application/02_running_sims.html#constraints')
        return form

    def _defineNonBondedParams(self, form):
        form.addParam('nonbondedMethod', params.EnumParam, default=0,
                      choices=['NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'Ewald', 'PME', 'LJPME'],
                      label="Non bonded method: ",
                      help='Non bonded method to simulate the non bonded atom interactions')
        form.addParam('nonbondedCutoff', params.FloatParam, default=1.0, expertLevel=params.LEVEL_ADVANCED,
                      label='Distance cutoff for non bonded interactions (nm): ', condition='nonbondedMethod!=0',
                      help='TThe cutoff distance to use for nonbonded interactions')
        return form

    def _defineHydrogenParams(self, form):
        form.addParam('addH', params.BooleanParam, default=False,
                      label='Add hydrogens to the system: ', help='Add hydrogens to the system')
        form.addParam('hPH', params.FloatParam, default=7.0, expertLevel=params.LEVEL_ADVANCED,
                      label='PH for hydrogen addition: ', help='The pH based on which to select variants')
        return form

    def _defineBoxParams(self, form):
        form.addParam('sizeType', params.EnumParam, label="System size type: ", default=1,
                      choices=['Absolute', 'Padding'], display=params.EnumParam.DISPLAY_HLIST,
                      help='Absolute: absolute size of the box (diameter)\n'
                             'Buffer: distance from the solute to the edge of the box\n')
        line = form.addLine('Box size (nm):', condition='sizeType == 0',
                            help='Distances of the bounding box (nm).\nIf BSS, then it will be the value of the '
                                 'image distance')
        line.addParam('distA', params.FloatParam, default=5.0, label='a: ')
        line.addParam('distB', params.FloatParam, default=5.0, label='b: ')
        line.addParam('distC', params.FloatParam, default=5.0, label='c: ')
        form.addParam('padDist', params.FloatParam, condition='sizeType == 1',
                      default=1.0, label='Padding distance: ',
                      help='Distance (nm) from the solute to the edge of the box.')
        return form

    def _defineSaltParams(self, form):
        form.addParam('saltConc', params.FloatParam, default=0, label='Salt concentration (M): ',
                      help='Ionic strength to prepare the system')

        form.addParam('neutralize', params.BooleanParam, default=True, label='Neutralize system: ',
                      help='Whether to add ions to the system until neutralize.')

        form.addParam('cationType', params.EnumParam,
                      label='Cation to add: ', choices=CATION_NAMES, default=3,
                      help='Which cation to add in the system')

        form.addParam('anionType', params.EnumParam,
                      label='Anions to add: ', choices=ANION_NAMES, default=0,
                      help='Which anion to add in the system')
        return form

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        iGroup = form.addGroup('Input')
        iGroup.addParam('inputFrom', params.EnumParam, default=STRUCTURE,
                        label='Input from: ', choices=['AtomStruct', 'SetOfSmallMolecules'],
                        help='Type of input you want to use')
        iGroup.addParam('inputStructure', params.PointerParam, pointerClass='SchrodingerAtomStruct, AtomStruct',
                        label='Input structure to be prepared for MD:', condition='inputFrom==0',
                        help='Atomic structure to be prepared for MD by solvation, ions addition etc')
        iGroup.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input set of molecules:', condition=LIG_INPUT,
                        help='Input set of docked molecules. One of them will be prepared together with its target')
        iGroup.addParam('inputLigand', params.StringParam, condition=LIG_INPUT,
                        label='Ligand to prepare: ',
                        help='Specific ligand to prepare in the system')

        ffGroup = form.addGroup('System force fields')
        self._defineFFParams(ffGroup, ligandCondition=LIG_INPUT)

        ffGroup = form.addGroup('Non bonded interactions')
        self._defineNonBondedParams(ffGroup)

        ffGroup = form.addGroup('Hydrogens')
        self._defineHydrogenParams(ffGroup)

        # todo: allow the use of variants

        form.addSection(label='Solvent box')
        sGroup = form.addGroup('Boundary box')
        self._defineBoxParams(sGroup)

        iGroup = form.addGroup('Ions')
        self._defineSaltParams(iGroup)


    def _insertAllSteps(self):
      self._insertFunctionStep(self.solvateStep)
      self._insertFunctionStep(self.createOutputStep)


    def solvateStep(self):
      recFile = self.getReceptorPDB()
      molFile = self.getSpecifiedMolFile() if self.inputFrom.get() == LIGAND else None

      with open(self.getParamsFile(), 'w') as f:
        f.write(f'receptorFile :: {recFile}\n')
        if molFile:
          f.write(f'ligandFile :: {molFile}\n')
          f.write(f'ligandFF :: {self.getLigandFFVersion()}\n')

        f.write(self.getFFParams())

      Plugin.runScript(self, 'openmmPrepareSystem.py', args=self.getParamsFile(), env=OPENMM_DIC, cwd=self._getPath())

    def createOutputStep(self):
      systemBasename = self.getSystemName()
      outStructFile = self._getPath(f'{systemBasename}_system.pdb')
      outCifFile = self._getPath(f'{systemBasename}_system.cif')
      outSystemFile = self._getPath(f'{systemBasename}_system.xml')

      ligName = self.inputLigand.get() if self.inputFrom.get() == LIGAND else None
      mFF, wFF = self.getFFFiles()
      outSystem = OpenMMSystem(filename=outStructFile, oriStructFile=outStructFile,
                               cifFile=outCifFile, serieFile=outSystemFile,
                               ff=mFF, wff=wFF, ligName=ligName)

      self._defineOutputs(outputSystem=outSystem)
      # self._defineSourceRelation(self.inputStructure, outSystem)

    def getFFParams(self):
        ffStr = ''
        mFF, wFF = self.getFFFiles()
        ffStr += f'mFF :: {mFF}\nwFF :: {wFF}\n'
        ffStr += f'nonbondedMethod :: {self.getEnumText("nonbondedMethod")}\n'
        ffStr += f'nonbondedCutoff :: {self.nonbondedCutoff.get()}\n'
        ffStr += f'constraints :: {self.getEnumText("constraints")}\n'

        wModel = self.getWaterModel(wFF)
        ffStr += f'wModel :: {wModel}\n'

        ffStr += f'addH :: {self.addH.get()}\n'
        if self.addH.get():
          ffStr += f'hPH :: {self.hPH.get()}\n'

        if self.sizeType.get() == 0:
          ffStr += f'boxSize :: {self.distA.get()}, {self.distB.get()}, {self.distC.get()}\n'
        else:
          ffStr += f'padDist :: {self.padDist.get()}\n'

        ffStr += f'saltConc :: {self.saltConc.get()}\n'
        ffStr += f'neutralize :: {self.neutralize.get()}\n'
        ffStr += f'cationType :: {self.getEnumText("cationType")}\n'
        ffStr += f'anionType :: {self.getEnumText("anionType")}\n'
        return ffStr

    def getLigandFileDir(self):
      lDir = os.path.abspath(self._getExtraPath('ligand'))
      if not os.path.exists(lDir):
        os.mkdir(lDir)
      return lDir

    def getLigParamFile(self):
      return os.path.abspath(self._getExtraPath('addHydrogens.txt'))

    def writePrepParamsFile(self, molFiles):
        paramsFile = self.getLigParamFile()
        with open(paramsFile, 'w') as f:
            molFilesStr = ' '.join(molFiles)
            f.write(f"ligandFiles: {molFilesStr}\n")

            f.write(f'outputDir: {self.getLigandFileDir()}\n')
            f.write('doHydrogens: True\n')
            f.write('doGasteiger: False\n')
            f.write('sanitize: False\n')
        return paramsFile

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

    def getReceptorPDB(self):
      recPDB = os.path.abspath(self._getExtraPath(f'{self.getSystemName()}.pdb'))
      if not os.path.exists(recPDB):
        recFile = self.getReceptorFilename()
        args = f'{recFile} --output {recPDB}'
        pwchemPlugin.runOPENBABEL(self, 'pdbfixer', args=args, cwd=self._getExtraPath())
      return recPDB

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
            sdfFile = convertToSdf(self, molFile)
            paramFile = self.writePrepParamsFile([sdfFile])
            pwchemPlugin.runScript(self, scriptLigPrepName, paramFile, env=RDKIT_DIC, cwd=self._getPath())
            return os.path.join(self.getLigandFileDir(), os.listdir(self.getLigandFileDir())[0])

    def _warnings(self):
      ws = []
      if self.constraints.get() == 0:
        ws.append('Running the simulation without restraints might lead to errors in the simulation.\n')
      return ws
