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
This protocol will calculate the Coulomb and LJ energies of interaction between a protein receptor
and a docked ligand. It can be used to get the average interaction energy in a trajectory.
"""
import os
import numpy as np

from pyworkflow.protocol import params
from pyworkflow.utils import Message
import pyworkflow.object as pwobj

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import getBaseName, convertToSdf
from pwchem.constants import RDKIT_DIC

from openmm import Plugin
from openmm.constants import OPENMM_DIC
from openmm.protocols import ProtOpenMMSystemPrep, ProtOpenMMSystemSimulation

SYSTEM, MOLSET = 0, 1
scriptLigPrepName = 'rdkit_addHydrogens.py'

class ProtOpenMMInteractionEnergy(ProtOpenMMSystemPrep, ProtOpenMMSystemSimulation):
    """
    This protocol will calculate the interaction energy of protein and ligand in a system

    AI Generated:

        ProtOpenMMInteractionEnergy - User Manual

        Overview
        --------
        This protocol calculates the interaction energy between a protein receptor and
        docked ligand(s) using OpenMM simulations. It evaluates both Coulomb (electrostatic)
        and Lennard-Jones (van der Waals) interaction energies over the course of a trajectory,
        allowing users to obtain average energies for further analysis.

        Inputs
        ------
        - **inputFrom**: Enum specifying input type:
            - `OpenMMSystem`: Use an existing OpenMMSystem for energy evaluation.
            - `SetOfSmallMolecules`: Use a set of docked small molecules that will be
              prepared, solvated, and simulated.
        - **inputSystem**: OpenMMSystem object (required if `inputFrom==OpenMMSystem`),
          containing structure, trajectory, and topology.
        - **inputSetOfMols**: SetOfSmallMolecules object (required if `inputFrom==SetOfSmallMolecules`),
          representing docked ligands.

        Options
        -------
        - **Minimization**: Add energy minimization for systems without trajectories.
        - **Integrator parameters**: Define integration scheme, step size, temperature,
          and friction coefficient for MD simulations.
        - **Forcefield parameters**: Specify force fields for receptor and ligands,
          non-bonded interaction methods, and hydrogen treatment.
        - **Solvent box and ions**: Configure solvation box, ion concentrations, and
          boundary conditions when preparing systems from small molecules.

        Workflow
        --------
        1. **Input Conversion** (if using SetOfSmallMolecules):
           - Converts ligand poses to SDF format.
           - Prepares ligands with hydrogen atoms using RDKit.

        2. **System Preparation**:
           - Generates solvated systems with force field parameters.
           - Optionally adds ions and solvent box.

        3. **Simulation**:
           - Runs MD simulation on prepared systems.
           - Computes Coulomb and Lennard-Jones interaction energies.
           - For existing OpenMMSystems, it can operate directly on trajectories.

        4. **Output Creation**:
           - For OpenMMSystem inputs: updates the report file with interaction energies.
           - For small molecule inputs: creates a copy of each molecule and stores
             Coulomb and LJ energies as Float objects.

        Outputs
        -------
        - **OpenMMSystem** (if inputFrom==OpenMMSystem): system with updated report file
          including interaction energies.
        - **SetOfSmallMolecules** (if inputFrom==SetOfSmallMolecules): copy of input molecules
          with `coulomb_Interaction` and `lj_Interaction` properties populated.

        Practical Recommendations
        -------------------------
        - Use this protocol after docking to evaluate binding energies of ligands.
        - When using large sets of molecules, ensure sufficient computational resources
          as solvation and simulation are computationally intensive.
        - Verify that molecules are docked before running interaction energy calculations.

        Summary & Interpretation
        ------------------------
        - Average Coulomb and LJ energies are computed from the trajectory or simulation output.
        - Standard deviations are reported if multiple frames are available.
        - Outputs can be used for ranking ligands or further post-processing analyses.

        Warnings
        --------
        - Simulating more than 50 molecules from SetOfSmallMolecules can be computationally
          expensive.
        - Ensure molecules are docked; otherwise, interaction energies cannot be calculated.
    """
    _label = 'system interaction energy'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputFrom', params.EnumParam, default=SYSTEM,
                      label='Input from: ', choices=['OpenMMSystem', 'SetOfSmallMolecules'],
                      help='Type of input you want to use')
        form.addParam('inputSystem', params.PointerParam, label="Input structure: ", condition=f'inputFrom=={SYSTEM}',
                      important=True, pointerClass='OpenMMSystem', help='OpenMMSystem to execute the calculation over')
        form.addParam('inputSetOfMols', params.PointerParam, label="Input docked molecules: ",
                      condition=f'inputFrom=={MOLSET}', important=True, pointerClass='SetOfSmallMolecules',
                      help='Input docked molecules to execute the interaction analysis over')

        mGroup = form.addGroup('Minimization',
                               condition=f'(inputSystem and not inputSystem.hasTrajectory()) or inputFrom=={MOLSET}',
                               help='Add energy minimization to the original system if there is no trajectory')
        self._defineMinimization(mGroup)

        iGroup = form.addGroup('Integrator')
        self._defineIntegrator(iGroup)

        form.addSection("System preparation forcefield")
        ffGroup = form.addGroup('System force fields', condition=f'inputFrom=={MOLSET}')
        self._defineFFParams(ffGroup)

        ffGroup = form.addGroup('Non bonded interactions', condition=f'inputFrom=={MOLSET}')
        self._defineNonBondedParams(ffGroup)

        ffGroup = form.addGroup('Hydrogens', condition=f'inputFrom=={MOLSET}')
        self._defineHydrogenParams(ffGroup)

        form.addSection("System preparation solvent box")
        sGroup = form.addGroup('Boundary box', condition=f'inputFrom=={MOLSET}')
        self._defineBoxParams(sGroup)

        iGroup = form.addGroup('Ions', condition=f'inputFrom=={MOLSET}')
        self._defineSaltParams(iGroup)

        form.addParallelSection(threads=4, mpi=1)


    def _insertAllSteps(self):
        simSteps = []
        if self.inputFrom.get() == MOLSET:
            cStep = self._insertFunctionStep(self.convertInputStep)

            for mol in self.inputSetOfMols.get():
                molFile = mol.getPoseFile()
                molFile = os.path.join(self.getLigandFileDir(), getBaseName(molFile) + '.sdf')
                sStep = self._insertFunctionStep(self.solvateStep, molFile, prerequisites=[cStep])
                simSteps.append(self._insertFunctionStep(self.simulateStep, molFile, prerequisites=[sStep]))
        else:
            simSteps.append(self._insertFunctionStep(self.simulateStep))
        self._insertFunctionStep(self.createOutputStep, prerequisites=simSteps)

    def convertInputStep(self):
        sdfFiles = []
        for mol in self.inputSetOfMols.get():
            molFile = mol.getPoseFile()
            sdfFiles.append(os.path.abspath(convertToSdf(self, molFile)))

        paramFile = self.writePrepParamsFile(sdfFiles)
        pwchemPlugin.runScript(self, scriptLigPrepName, paramFile, env=RDKIT_DIC, cwd=self._getPath())

    def solvateStep(self, molFile):
        if not os.path.exists(molFile):
            convFile = self._getTmpPath(getBaseName(molFile) + '.sdf')
            if os.path.exists(convFile):
                os.rename(convFile, molFile)
            else:
                print(f'Molecule {getBaseName(molFile)} could not be managed')

        if os.path.exists(molFile):
            recFile = self.getReceptorPDB()
            molBase = getBaseName(molFile)
            oDir = self._getExtraPath(molBase)
            os.mkdir(oDir)

            paramsFile = self.getSolvateParamsFile(oDir)
            with open(paramsFile, 'w') as f:
                f.write(f'receptorFile :: {recFile}\n')
                f.write(f'ligandFile :: {molFile}\n')
                f.write(f'ligandFF :: {self.getLigandFFVersion()}\n')

                f.write(self.getFFParams())

            Plugin.runScript(self, 'openmmPrepareSystem.py', args=paramsFile, env=OPENMM_DIC, cwd=oDir)


    def simulateStep(self, molFile=None):
        if not molFile or os.path.exists(molFile):
            if molFile:
                molBase = getBaseName(molFile)
                oDir = self._getExtraPath(molBase)
            else:
                oDir = self._getPath()

            sysFile, structFile = self.getSerieFile(molFile), self.getStructureFile(molFile)
            trajFile = self.getSystemTrajFile()

            paramsFile = self.getInteractionParamsFile(oDir)
            with open(paramsFile, 'w') as f:
              f.write(f'systemFile :: {sysFile}\n')
              f.write(f'structureFile :: {structFile}\n')
              if trajFile:
                  f.write(f'trajFile :: {os.path.abspath(trajFile)}\n')
              else:
                  f.write(f'addMin :: {self.addMinimization.get()}\n')
                  if self.addMinimization.get():
                    f.write(f'minimTol :: {self.minimTol.get()}\n')
                    f.write(f'maxIter :: {self.maxIter.get()}\n')

              f.write(f'mFF :: {self.getSystemFF()}\n')

              integrator = self.getEnumText('integrator')
              f.write(f'integrator :: {integrator}\n')
              if self.integrator.get() not in [0, 5]:
                  f.write(f'temperature :: {self.temperature.get()}\n')

              if self.integrator.get() not in [5, 6]:
                  f.write(f'stepSize :: {self.stepSize.get()}\n')

              if self.integrator.get() not in [0, 3, 5]:
                  f.write(f'fricCoef :: {self.fricCoef.get()}\n')

            Plugin.runScript(self, 'openmmInteractionEnergy.py', args=paramsFile, env=OPENMM_DIC, cwd=oDir)


    def createOutputStep(self):
        if self.inputFrom.get() == SYSTEM:
            outSystem = self.inputSystem.get().clone()
            if not outSystem.hasTrajectory() and self.addMinimization.get():
                outSystem.setFileName(self._getPath(f'{self.getSystemName()}.pdb'))
                outSystem.setCifFile(self._getPath(f'{self.getSystemName()}.cif'))
            elif outSystem.hasTrajectory():
                repFile = outSystem.getReportFile()

                data = np.loadtxt(repFile, delimiter=',', ndmin=2)
                cEs, ljEs = self.parseEnergies(self.getOutputFile())
                cEs, ljEs = np.array(cEs).reshape(-1, 1), np.array(ljEs).reshape(-1, 1)
                data = np.hstack((data, cEs, ljEs))

                newRepFile = self._getPath('md_log.txt')
                headerStr = self.getHeaderStr(repFile) + \
                            ',"Interaction Coulomb Energy (KJ/mol)","Interaction LJ Energy (KJ/mol)"'
                np.savetxt(newRepFile, data, delimiter=",", comments="", fmt="%f", header=headerStr)
                outSystem.setReportFile(newRepFile)

            self._defineOutputs(outputSystem=outSystem)

        else:
            inMols = self.inputSetOfMols.get()
            outputSet = inMols.createCopy(self._getPath(), copyInfo=True)
            for mol in inMols:
                molBase = getBaseName(mol.getPoseFile())
                outFile = self.getOutputFile(molBase)
                if os.path.exists(outFile):
                    cEs, ljEs = self.parseEnergies(outFile)
                    mol.coulomb_Interaction = pwobj.Float(cEs[0])
                    mol.lj_Interaction = pwobj.Float(ljEs[0])
                else:
                    mol.coulomb_Interaction = pwobj.Float(None)
                    mol.lj_Interaction = pwobj.Float(None)

                outputSet.append(mol)
            self._defineOutputs(outputSmallMolecules=outputSet)


####################### UTILS FUNCTIONS ############################

    def parseEnergies(self, resFile):
      with open(resFile) as f:
        coulombEnergies = [float(energy) for energy in f.readline().split(':')[1].split()]
        ljEnergies = [float(energy) for energy in f.readline().split(':')[1].split()]
      return coulombEnergies, ljEnergies

    def getSolvateParamsFile(self, oDir):
      paramsFile = os.path.abspath(os.path.join(oDir, 'solvationParams.txt'))
      return paramsFile

    def getInteractionParamsFile(self, oDir):
      paramsFile = os.path.abspath(os.path.join(oDir, 'interactionParams.txt'))
      return paramsFile

    def getStructureFile(self, molFile=None):
      if not molFile:
          sysFile = os.path.abspath(self.inputSystem.get().getCifFile())
      else:
          molBase = getBaseName(molFile)
          sysFile = os.path.abspath(self._getExtraPath(f'{molBase}/{self.getSystemName()}_system.cif'))
      return sysFile

    def getSerieFile(self, molFile=None):
      if not molFile:
          sysFile = os.path.abspath(self.inputSystem.get().getSerieFile())
      else:
          molBase = getBaseName(molFile)
          sysFile = os.path.abspath(self._getExtraPath(f'{molBase}/{self.getSystemName()}_system.xml'))
      return sysFile

    def getSystemName(self):
      if self.inputFrom.get() == SYSTEM:
        sysName = self.inputSystem.get().getSystemName()
      else:
        sysName = getBaseName(self.getReceptorFilename())
      return sysName

    def getSystemTrajFile(self):
      if self.inputFrom.get() == SYSTEM:
          trajFile = self.inputSystem.get().getTrajectoryFile()
      else:
          trajFile = None
      return trajFile

    def getSystemFF(self):
      if self.inputFrom.get() == SYSTEM:
          mFF = self.inputSystem.get().getForceField()
      else:
          mFF, _ = self.getFFFiles()
      return mFF

    def getHeaderStr(self, repFile):
      with open(repFile) as f:
        headerLine = f.readline().strip()
      return headerLine

####################### SUMMARY FUNCTIONS ############################

    def getSummaryStr(self, coulombEnergies, ljEnergies):
      avgCo, avgLj = np.mean(coulombEnergies), np.mean(ljEnergies)
      if len(ljEnergies) > 1:
        stdCo, stdLj = np.std(coulombEnergies), np.std(ljEnergies)
      else:
        stdCo, stdLj = 0, 0

      energyStr = f'Average Coulomb energy:\t\t{avgCo:.4f} ± {stdCo:.4f} kJ/mol\n' \
                  f'Average LJ energy:\t\t{avgLj:.4f} ± {stdLj:.4f} kJ/mol"\n'
      return energyStr

    def getOutputFile(self, molBase=None):
      if not molBase:
          outFile = self._getPath('energy_results.tsv')
      else:
          outFile = self._getExtraPath(os.path.join(molBase, 'energy_results.tsv'))
      return outFile

    def _summary(self):
      s = []
      resFile = self.getOutputFile()
      if os.path.exists(resFile):
          s = self.getSummaryStr(*self.parseEnergies(resFile))
      return s

    def _warnings(self):
      ws = []
      if self.inputFrom.get() == MOLSET and len(self.inputSetOfMols.get()) > 50:
          ws.append(f'This evaluation needs to solvate and simulate the molecules, which is computationally expensive.'
                    f'Do you really want to run the protocol over your {len(self.inputSetOfMols.get())} molecules?')
      return ws


    def _validate(self):
      vs = []
      if self.inputFrom.get() == MOLSET and not self.inputSetOfMols.get().isDocked():
          vs.append('Molecules must be docked to a receptor to calculate their interaction energy')
      return vs