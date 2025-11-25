# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
This module will run the simulation for a constant pH
"""
import os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem import Plugin as pwchemPlugin

from .. import Plugin
from ..constants import OPENMM_DIC
from ..objects import OpenMMSystem

from openmm import *
from openmm.app import *
from openmm.unit import *


class ProtOpenMMSystemSimulationConstantPH(EMProtocol):
    """
    This protocol will start a Molecular Dynamics simulation with constant pH specified by user.
    """
    _label = 'constant pH system simulation'
    stepsExecutionMode = params.STEPS_PARALLEL
    ASP_VAR = {1: ['ASP', 'ASH']}
    GLU_VAR = {1: ['GLU', 'GLH']}
    CYS_VAR = {1: ['CYS', 'CYX']}
    HID_VAR = {1: ['HIP', 'HID']}
    HIE_VAR = {1: ['HIP', 'HIE']}
    LYS_VAR = {1: ['LYS', 'LYN']}

    EXPLICIT_PARAMS = dict(nonbondedMethod=PME,
                           nonbondedCutoff=0.9 * nanometers,
                           constraints=HBonds,
                           hydrogenMass=1.5 * amu)

    IMPLICIT_PARAMS = dict(nonbondedMethod=CutoffNonPeriodic,
                           nonbondedCutoff=2.0 * nanometers,
                           constraints=HBonds)

    # -------------------------- DEFINE param functions ----------------------
    def _defineMinimization(self, form):
        form.addParam('addMinimization', params.BooleanParam, default=True, label="Add minimization: ",
                      help='Add energy minimization')
        form.addParam('minimTol', params.FloatParam, default=10, label="Minimization tolerance (kJ/mol): ",
                      condition='addMinimization',
                      help='This specifies how precisely the energy minimum must be located.  Minimization is halted '
                           'once the root-mean-square value of all force components reaches this tolerance.')
        form.addParam('maxIter', params.IntParam, default=10000, label="Maximum iterations: ",
                      condition='addMinimization',
                      help='The maximum number of iterations to perform.  If this is 0, minimization is continued until'
                           ' the results converge without regard to how many iterations it takes.')
        return form

    def _defineIntegrator(self, form):
        form.addParam('integrator', params.EnumParam, default=1, label="Simulation integrator: ",
                      choices=['Verlet', 'Langevin', 'LangevinMiddle', 'NoseHoover', 'Brownian', 'VariableVerlet',
                               'VariableLangevin'],
                      help='http://docs.openmm.org/latest/userguide/theory/04_integrators.html')

        form.addParam('stepSize', params.FloatParam, default=0.004, label="Step size for integration (ps): ",
                      condition='not integrator in [5, 6]',
                      help='The step size with which to integrate the system (in picoseconds)')
        form.addParam('fricCoef', params.FloatParam, default=1, label="Friction coefficient (1/ps): ",
                      condition='integrator in [1, 2, 4, 6]',
                      help='The friction coefficient which couples the system to the heat bath (in inverse picoseconds)')
        form.addParam('temperature', params.FloatParam, default=300, label="Simulation temperature (K): ",
                      condition='integrator in [1, 2, 3, 4, 6]', help='Temperature for the simulation')
        form.addParam('colFreq', params.FloatParam, default=1, label="Collision frequency (1/ps): ",
                      condition='integrator in [3]',
                      help='The friction coefficient which couples the system to the heat bath (in inverse picoseconds)')

        form.addParam('errTol', params.FloatParam, default=0.001, label="Error tolerance: ",
                      condition='integrator in [5, 6]', help='The error tolerance')
        return form

    def _defineBarostat(self, form):
        form.addParam('addBarostat', params.BooleanParam, default=False, label="Add barostat: ",
                      help='Add MonteCarlo Barostat to run a NPT simulation')
        form.addParam('pressure', params.FloatParam, default=1, label="Pressure (bar): ", condition='addBarostat',
                      help='The default pressure acting on the system (in bar)')
        form.addParam('barFreq', params.IntParam, default=25, label="Barostat frequency: ",
                      condition='addBarostat',
                      help='The frequency at which Monte Carlo pressure changes should be attempted (in time steps)')
        return form

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution: ",
                       help="This protocol has both CPU and GPU implementation.\
                                                 Select the one you want to use.")
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSystem', params.PointerParam, label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='OpenMMSystem', help='OpenMMSystem to execute the simulation over')
        form.addParam('nSteps', params.IntParam, default=10000, label="Number of simulation steps: ",
                      help='Number of steps for simulation')

        phGroup = form.addGroup('pH ')
        phGroup.addParam('singlePH', params.BooleanParam, default=True, label="Use one single pH value: ",
                      help='Choose whether to use one single pH value or an array of them.')
        phGroup.addParam('onePH', params.FloatParam, default=7.5, label="pH value: ", condition='singlePH',
                      help='The pH value to use.')
        phGroup.addParam('manyPH', params.StringParam, default='6.5, 7.0, 7.5, 8.0, 8.5', label="pH values: ", condition='not singlePH',
                      help='The pH values to use, separated with commas.')

        mGroup = form.addGroup('Minimization')
        self._defineMinimization(mGroup)

        iGroup = form.addGroup('Integrator')
        self._defineIntegrator(iGroup)

        bGroup = form.addGroup('Barostat')
        self._defineBarostat(bGroup)


        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
      self._insertFunctionStep(self.prepareTitrattionStep)
      #if self.useOpenmmdl.get() and self.inputSystem.get().getLigTopologyFile():
      #  self._insertFunctionStep(self.analyzeStep)
      #self._insertFunctionStep(self.createOutputStep)
      pass

    def prepareTitrationStep(self):
        home = Plugin.getVar(OPENMM_DIC['home'])
        repoPath = os.path.join(home, 'openmm-cph')
        if repoPath not in sys.path:
            sys.path.insert(0, repoPath)

        from constantph import ConstantPH
        from reference_energy import ReferenceEnergyFinder

        cifFile = self.inputSystem.get().getCifFile()
        structure = MMCIFFile(cif_file)
        topology = structure.topology

        variants, referenceEnergies = self.getVarsAndRefEnergies(topology)
        
        #ph values
        if self.singlePH.get() :
            ph = self.onePH.get()
        else:
            ph = self.getListOfPH()
        #force fields
        explicitFFfiles, implicitFFfiles = self.getFFFiles()
        explicitFF = ForceField(*explicitFFfiles)
        implicitFF = ForceField(*implicitFFfiles)
        explicitParams = self.EXPLICIT_PARAMS
        implicitParams = self.IMPLICIT_PARAMS
        #integrators
        integrator, relaxationIntegrator = self.getIntegrators()
        
        cph = ConstantPH(topology, structure.positions, ph, explicitFF, implicitFF,
                         variants, referenceEnergies, 100, explicitParams, implicitParams, integrator, relaxationIntegrator)

        #barostat
        if self.addBarostat.get():
            cph.simulation.system.addForce(MonteCarloBarostat(self.pressure.get() * bar,
                                                              temperature))
            cph.simulation.context.reinitialize(preserveState=True)

        #minimize
        if self.addMinimization.get():
            print("Minimizing energy...")
            cph.simulation.minimizeEnergy(tolerance=self.minimTol.get() * kilojoules_per_mole,
                                          maxIterations=self.maxIter.get())
            print("Energy minimization complete.")


    def createOutputStep(self):
      systemName = self.getSystemName()
      systemFile = os.path.relpath(self.getSystemFile())
      outTopFile, outDcdFile = self._getPath(f'{systemName}.pdb'), self._getPath(f'{systemName}.dcd')
      outCifFile = self._getPath(f'{systemName}.cif')

      mFF, wFF = self.getFFFiles()
      nFrames = self.getNFrames()
      nTime = nFrames * self.stepSize.get()
      outSystem = OpenMMSystem(filename=outTopFile, serieFile=systemFile, cifFile=outCifFile,
                               repFile=self._getPath('md_log.txt'),
                               ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)
      outSystem.setTrajectoryFile(outDcdFile)

      ligFile = self.inputSystem.get().getLigTopologyFile()
      if ligFile:
        outSystem.setLigTopologyFile(ligFile)

      anaDir = self._getExtraPath('OpenMMDL')
      if os.path.exists(anaDir):
        outSystem.setOpenmmdlDir(anaDir)

      self._defineOutputs(outputSystem=outSystem)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
      summary = []
      return summary

    def _methods(self):
      methods = []
      return methods

    def _warnings(self):
      ws = []
      if not self.addMinimization.get():
        ws.append('Running the simulation without a prior minimization might lead to errors in the simulation.\n')
      return ws

    # --------------------------- UTILS functions -----------------------------------
    def getListOfPH(self):
      userInput = self.manyPH.get()
      listPH = userInput.split(',')
      return listPH

    def getIntegrators(self):
        temperature = self.temperature.get() * kelvin
        stepSize = self.stepSize.get() * picoseconds
        fricCoef = self.fricCoef.get() / picosecond
        integrator = LangevinIntegrator(temperature, fricCoef, stepSize)
        relaxationIntegrator = LangevinIntegrator(temperature, 10.0 / picosecond, 0.002 * picoseconds)
        return integrator, relaxationIntegrator

    def computeReferenceEnergies(self, pdbFile, variantsDict, targetpKa):
        print(f"Computing reference energies for {pdbFile} (target pKa={targetpKa})")
        pdb = PDBFile(pdbFile)
        referenceEnergies = {index: [0.0] * len(states) for index, states in variantsDict.items()}

        integrator, relaxationIntegrator = self.getIntegrators()

        cph = ConstantPH(pdb.topology, pdb.positions, 7.0,
                         explicitFF, implicitFF,
                         variantsDict, referenceEnergies, 250,
                         self.EXPLICIT_PARAMS, self.IMPLICIT_PARAMS,
                         integrator, relaxationIntegrator)

        finder = ReferenceEnergyFinder(cph, targetpKa, TEMPERATURE)

        total_iterations = 20000
        chunk = 200
        for start in range(0, total_iterations, chunk):
            print(f"  Iterations {start}?{start + chunk}...")
            finder.findReferenceEnergies(iterations=chunk, substeps=10)

        ref_energies = {index: cph.titrations[index].referenceEnergies for index in variantsDict}
        print(f"Reference energies for {pdbFile} computed.")
        return ref_energies

    def getVarsAndRefEnergies(self, structure):
        referenceEnergies = {}
        variants = {}
        #get the variants energies
        home = Plugin.getVar(OPENMM_DIC['home'])
        aspPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/ASP.pdb')
        aspRefEnergies = computeReferenceEnergies(aspPDB, self.ASP_VAR, 3.9)
        gluPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/GLU.pdb')
        gluRefEnergies = computeReferenceEnergies(glupPDB, self.GLU_VAR, 4.2)
        cysPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/CYS.pdb')
        cysRefEnergies = computeReferenceEnergies(cysPDB, self.CYS_VAR, 8.3)
        hisPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/HIS.pdb')
        hidRefEnergies = computeReferenceEnergies(hisPDB, self.HID_VAR, 7.1)
        hieRefEnergies = computeReferenceEnergies(hiePDB, self.HIE_VAR, 6.5)
        lysPDB = os.path.abspath(f'{home}/openmm-cph/model-compounds/LYS.pdb')
        lysRefEnergies = computeReferenceEnergies(lysPDB, self.LYS_VAR, 10.5)
        
        #prepare titration
        for residue in structure.residues:
            if residue.name == 'ASP':
                variants[residue.index] = ['ASP', 'ASH']
                referenceEnergies[residue.index] = aspRefEnergies[1]
            elif residue.name == 'GLU':
                variants[residue.index] = ['GLU', 'GLH']
                referenceEnergies[residue.index] = gluRefEnergies[1]
            elif residue.name == 'CYS':
                variants[residue.index] = ['CYS', 'CYX']
                referenceEnergies[residue.index] = cysRefEnergies[1]
            elif residue.name == 'HIS':
                variants[residue.index] = ['HIP', 'HID', 'HIE']
                referenceEnergies[residue.index] = [0.0*kilojoules_per_mole, 
                                                    hidRefEnergies[1][1], 
                                                    hieRefEnergies[1][1]]
            elif residue.name == 'LYS':
                variants[residue.index] = ['LYS', 'LYN']
                referenceEnergies[residue.index] = lysRefEnergies[1]
            
        return variants, referenceEnergies

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
      system = self.inputSystem.get()
      return system.getForceField(), system.getWaterForceField()

    def getNBParams(self):
      system = self.inputSystem.get()
      return system._nbMethod.get(), system._nbCutoff.get()

    def getParamsFile(self):
      return os.path.abspath(self._getExtraPath('simulationParams.txt'))

    def getStructureFile(self):
      return os.path.abspath(self.inputSystem.get().getCifFile())

    def getSystemFile(self):
      return os.path.abspath(self.inputSystem.get().getSerieFile())

    def getSystemName(self):
      return self.inputSystem.get().getSystemName()
