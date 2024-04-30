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
# Description: This script installs CAMB as Commander3 subproject. 
#================================================================================

#------------------------------------------------------------------------------
# Note: the explicit splitting for download and install step is done on purpose
# to avoid errors when you want to recompile libraries for different owls etc.
# In addition, this will allow us to download sources only once and then just 
# reuse it whenever possible.
#------------------------------------------------------------------------------
# Getting CAMB from source
#------------------------------------------------------------------------------
# Checking whether we have source directory and this directory is not empty.
if(NOT EXISTS "${CAMB_SOURCE_DIR}/CMakeLists.txt")
	message(STATUS "No CAMB sources were found; thus, will download it from source:\n${camb_git_url}")
	ExternalProject_Add(
		camb_src
		GIT_REPOSITORY		"${camb_git_url}"
		GIT_TAG						"${camb_git_tag}"
		PREFIX						"${LIBS_BUILD_DIR}"
		DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
		SOURCE_DIR				"${CAMB_SOURCE_DIR}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD			ON
		# commands how to build the project
		CONFIGURE_COMMAND ""
		BUILD_COMMAND			""
		INSTALL_COMMAND		""
		)
else()
	message(STATUS "Found an existing CAMB sources inside:\n${CAMB_SOURCE_DIR}")
	add_custom_target(camb_src
		ALL ""
		)
endif()
#------------------------------------------------------------------------------
# Compiling and installing CAMB
#------------------------------------------------------------------------------
ExternalProject_Add(
	camb
	DEPENDS						required_libraries 
										curl
										cfitsio
										healpix	
										camb_src
	PREFIX						"${LIBS_BUILD_DIR}"
	SOURCE_DIR				"${CAMB_SOURCE_DIR}"
	INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}" 
	LOG_DIR						"${CMAKE_LOG_DIR}"
	LOG_CONFIGURE			ON
	LOG_BUILD					ON
	LOG_INSTALL				ON 
	# Commadns to build the project
	DOWNLOAD_COMMAND	""
	CMAKE_ARGS
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		# Specifying installations paths for binaries and libraries
		-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		-DBUILD_SHARED_LIBS:BOOL=OFF
		# Specifying compilers
		-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER}
		-DMPI_C_COMPILER=${MPI_C_COMPILER}
		-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
		# Check submodules during build
		-DGIT_SUBMODULE:BOOL=ON
		# CFitsIO paths
		#-DCFITSIO_FOUND:BOOL=TRUE
		-DCFITSIO_LIBRARIES:FILEPATH=${CFITSIO_LIBRARIES}
		-DCFITSIO_INCLUDE_DIRS:PATH=${CFITSIO_INCLUDE_DIRS}
		# HEALPix Paths
		-DHEALPIX_LIBRARIES:FILEPATH=${HEALPIX_LIBRARIES}
		-DHEALPIX_INCLUDE_DIRS:PATH=${HEALPIX_INCLUDE_DIRS}
	)
#------------------------------------------------------------------------------
# Including CAMB directories and Libraries
set(CAMB_INCLUDE_DIR
	"${CMAKE_INSTALL_PREFIX}/mod"
	)
set(CAMB_LIBRARIES
	${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}camb${CMAKE_STATIC_LIBRARY_SUFFIX}
	${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}forutils${CMAKE_STATIC_LIBRARY_SUFFIX}
	)
include_directories(${CAMB_INCLUDE_DIR})
#------------------------------------------------------------------------------
