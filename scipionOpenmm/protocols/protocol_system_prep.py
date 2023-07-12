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
from pwem.convert import AtomicStructHandler

from pwchem.utils import runOpenBabel, getBaseName
from pwchem import Plugin as pwchemPlugin

from scipionOpenmm import Plugin as openmmPlugin
from scipionOpenmm.constants import OPENMM_DIC
from scipionOpenmm.objects import OpenMMSystem


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
        form.addParam('inputStructure', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='AtomStruct', help='Atom structure to convert to OpenMM system')

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
      outDir = os.path.abspath(self._getExtraPath())
      inFile = self.getSystemFilename()
      systemBasename = self.getSystemName()

      with open(self.getParamsFile(), 'w') as f:
        f.write('inputFile :: {}\n'.format(inFile))
        mFF, wFF = self.getFFFiles()
        f.write('mFF :: {}\nwFF :: {}\n'.format(mFF, wFF))

        wModel = self.getWaterModel(wFF)
        f.write('wModel :: {}\n'.format(wModel))

        f.write('addH :: {}\n'.format(self.addH.get()))
        if self.addH.get():
          f.write('hPH :: {}\n'.format(self.hPH.get()))

        if self.sizeType.get() == 0:
          f.write('boxSize :: {}, {}, {}\n'.format(self.distA.get(), self.distB.get(), self.distC.get()))
        else:
          f.write('padDist :: {}\n'.format(self.padDist.get()))

        f.write('saltConc :: {}\n'.format(self.saltConc.get()))
        f.write('neutralize :: {}\n'.format(self.neutralize.get()))
        f.write('cationType :: {}\n'.format(self.getEnumText('cationType')))
        f.write('anionType :: {}\n'.format(self.getEnumText('anionType')))

      openmmPlugin.runScript(self, 'openmmPrepareSystem.py', args=self.getParamsFile(), env=OPENMM_DIC,
                             cwd=self._getPath())


    def createOutputStep(self):
      systemBasename = self.getSystemName()
      outSystemFile = self._getPath('{}_system.pdb'.format(systemBasename))

      mFF, wFF = self.getFFFiles()
      outSystem = OpenMMSystem(filename=outSystemFile, ff=mFF, wff=wFF,
                               nonbondedMethod=self.getEnumText('nonbondedMethod'),
                               nonbondedCutoff=self.nonbondedCutoff.get())

      self._defineOutputs(outputSystem=outSystem)
      self._defineSourceRelation(self.inputStructure, outSystem)


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

    def getParamsFile(self):
      return os.path.abspath(self._getExtraPath('solvationParams.txt'))

    def getSystemFilename(self):
      return os.path.abspath(self.inputStructure.get().getFileName())

    def getSystemName(self):
      return getBaseName(self.getSystemFilename())
