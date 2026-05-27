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
        cls._defineEmVar(OPENMM_DIC['home'], cls.getEnvName(OPENMM_DIC))
        cls._defineVar("OPENMM_ENV_ACTIVATION", cls.getEnvActivationCommand(OPENMM_DIC))

        cls._defineEmVar(ODUCK_DIC['home'], cls.getEnvName(ODUCK_DIC))

    @classmethod
    def defineBinaries(cls, env):
        cls.addOPENMMPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addODUCKPackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def addOPENMMPackage(cls, env, default=True):
        """ This function installs Espaloma package. """
        installer = InstallHelper(OPENMM_DIC['name'], packageHome=cls.getVar(OPENMM_DIC['home']),
                                  packageVersion=OPENMM_DIC['version'])
        home = cls.getEnvName(OPENMM_DIC)
        # Installing package
        installer.addCommand(
            f'conda create -n {home} -c conda-forge espaloma=0.4.0 openmm=8.3 cuda-version=12.8 openmmdl=1.2.0 '
            f'-y ', 'OPENMM_ENV_CREATED'
        ).addCommand(
            f'wget {cls.getEspalomaModelUrl()} -O {cls.getEspalomaModelFile()} ',
            'ESPALOMA_MODEL_DOWNLOADED'
        ).addCommand(
            f"cd {cls.getVar(OPENMM_DIC['home'])} && git clone https://github.com/openmm/openmm-cph.git",
            'CPH_REPO_CLONED'
        ).addPackage(env, dependencies=['conda'], default=default)


    @classmethod
    def addODUCKPackage(cls, env, default=True):
        """ This function installs Espaloma package. """
        installer = InstallHelper(ODUCK_DIC['name'], packageHome=cls.getVar(ODUCK_DIC['home']),
                                  packageVersion=ODUCK_DIC['version'])

        # Installing package
        installer.getCloneCommand(cls.getOpenDuckGithub(), targeName='ODUCK_CLONED'). \
            addCommand(f'{cls.getEnvActivationCommand(OPENMM_DIC)} && cd openduck && python setup.py install',
                       'ODUCK_INSTALLED'). \
            addCommand(f'python {cls.getScriptsDir("_updateOpenMMImports.py")} {ODUCK_DIC["name"]}',
                       'ODUCK_OPENMM_UPDATED'). \
            addPackage(env, dependencies=['conda'], default=default)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def getPluginHome(cls, path=""):
        import openmm
        fnDir = os.path.split(openmm.__file__)[0]
        return os.path.join(fnDir, path)

    @classmethod
    def runOpenMM(cls, protocol, program, args, cwd=None):
        """ Run OpenMM command from a given protocol. """
        fullProgram = f' {cls.getEnvActivationCommand(OPENMM_DIC)} && {program}'
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def getEspalomaModelUrl(cls):
        v = ESPALOMA_DIC["version"]
        return f'https://github.com/choderalab/espaloma/releases/download/{v}/espaloma-latest.pt'

    @classmethod
    def getEspalomaModelFile(cls):
        return cls.getPluginHome(f"models/espaloma-{ESPALOMA_DIC['version']}.pt")

    @classmethod
    def getOpenDuckGithub(cls):
        return "https://github.com/CBDD/openduck.git"

