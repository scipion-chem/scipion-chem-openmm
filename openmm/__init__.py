# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os

from scipion.install.funcs import InstallHelper

import pwchem

from .constants import *

_version_ = '0.1'
_logo = "openmm_logo.png"
_references = ['']

class Plugin(pwchem.Plugin):
    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(OPENMM_DIC['home'], '{}-{}'.format(OPENMM_DIC['name'], OPENMM_DIC['version']))
        cls._defineVar("OPENMM_ENV_ACTIVATION", cls.getEnvActivationCommand(OPENMM_DIC))

    @classmethod
    def defineBinaries(cls, env):
        cls.addOPENMMPackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def addOPENMMPackage(cls, env, default=True):
        installer = InstallHelper(OPENMM_DIC['name'], packageHome=cls.getVar(OPENMM_DIC['home']),
                                  packageVersion=OPENMM_DIC['version'])

        condaPackages = ['openmm={}'.format(OPENMM_DIC['version']), 'pdbfixer']

        installer.getCondaEnvCommand(requirementsFile=False). \
            addCondaPackages(condaPackages, channel='conda-forge'). \
            addPackage(env, dependencies=['conda'], default=default)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def getPluginHome(cls, path=""):
        import openmm
        fnDir = os.path.split(openmm.__file__)[0]
        return os.path.join(fnDir, path)

    @classmethod
    def runOpenMM(cls, protocol, program, args, cwd=None):
        """ Run Ambertools command from a given protocol. """
        fullProgram = ' %s && %s' % (cls.getEnvActivationCommand(OPENMM_DIC), program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runOpenMMScript(cls, protocol, program, args, cwd=None):
        """ Run Ambertools command from a given protocol. """
        fullProgram = ' %s && %s' % (cls.getEnvActivationCommand(OPENMM_DIC), program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

