# **************************************************************************
# *
# * Authors: Daniel Del Hoyo Gomez
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


OPENMM_DIC = {'name': 'openmm',    'version': '8.4', 'home': 'OPENMM_HOME'}
ESPALOMA_DIC = {'name': 'espaloma', 'version': '0.4.0', 'home': 'ESPALOMA_HOME'}
ODUCK_DIC = {'name': 'openduck',    'version': '0.2.0', 'home': 'ODUCK_HOME'}

# Open Free Energy (openfe) defaults for the RBFE/ABFE protocols
# Default values reproduce the benchmark protocol of Baumann et al., J. Chem. Inf. Model. 2026, 66
DEFAULT_TEMPERATURE = 298.15      # K
DEFAULT_SOLVENT_PADDING = 1.5     # nm
DEFAULT_PROTOCOL_REPEATS = 3
DEFAULT_EQUIL_LENGTH = 1.0        # ns of NPT equilibration per window
DEFAULT_RBFE_PRODUCTION = 5.0     # ns per lambda window (neutral transformations)
DEFAULT_RBFE_N_REPLICAS = 11      # lambda windows (neutral transformations)
DEFAULT_SMALL_MOL_FF = 'openff-2.2.0'   # Open Force Field Sage 2.2.0
DEFAULT_CHARGE_METHOD = 'am1bcc'        # AM1-BCC via AmberTools/Antechamber
DEFAULT_MINIMIZATION_STEPS = 5000       # openfe's own per-window default

# Floor for the solvent padding
MIN_SOLVENT_PADDING = 1.3               # nm

# `openfe gather --report dg` has TWO hard requirements
MIN_MLE_EDGES = 3
MIN_MLE_REPEATS = 2
