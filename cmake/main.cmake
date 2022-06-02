#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================
# Author: Maksym Brilenkov
#================================================================================
# Description: This script contains general instructions on how to fetch and build 
# Commander3 and all its dependencies. It is split into three parts, each containing 
# its set of instaructions/variables. It is done for easier maintenance. 
#================================================================================

# TODO:
# [x] Split this file into several different files, containing corresponding instructions (versions, variables, toolchains etc.);
# [ ] Change URL_MD5 to URL_HASH of every project;
# [x] Change compiler variables from list to string (but leave APPEND); <= doesn't work this way
# [ ] Remove include_directory() and use target_include_directory() instead (for commander3 target);
# [x] Add one variable which will force all libraries to be recompiled;
# [x] Change CFitsIO to CMake installation
# [ ] Write your own CFitsIO or modify an existing one, because it doesn't seem to work with version 4.0.0
# [x] Change FFTW to CMake installation
# [x] Finish Custom CMake module for CAMB as it is now works only with Intel compilers

#------------------------------------------------------------------------------
# including links where to download library sources
include(sources)
# including compiler definitions
include(toolchains)
# including project defined variables
include(variables)
# custom functions
include(functions)
# including required and dependent libraries
include(libraries)
# including configuration summary of installation
#include(summary)
#------------------------------------------------------------------------------

#message(STATUS "Finished looking for packages.")
#message(STATUS "---------------------------------------------------------------")
