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
# Module to find HEALPix on the system
# It looks for HEALPix components, such as:
# sharp, f90, cxx, c 
#
# Usage:
#   find_package(HEALPIX [version number] [REQUIRED] [QUIET] [COMPONENTS] [list of all components])
#
# Supported components are:
# SHARP, Fortran, CXX, C
#
# The following variables are defined:
#   HEALPIX_FOUND               set to TRUE if HEALPix is found on the system;
#   HEALPIX_[COMPONENT]_FOUND   set to TRUE if component is found on the system;
#   HEALPIX_LIBRARIES           complete list (with all found components) of HEALPIX libs' directory paths;
#   HEALPIX_[COMPONENT]_LIB     full path to one of the components;
#   HEALPIX_INCLUDE_DIRS        HEALPIX include directory paths; 
#================================================================================
# TODO: This file doesn't really work because we need to add *.mod and *.o files
# into includes (aka Fortran includes).
# Reference to CMake docs:
# https://cmake.org/cmake/help/v3.17/module/FindPackageHandleStandardArgs.html
include(FindPackageHandleStandardArgs)
include(SelectLibraryConfigurations)
#==============================================================================
# First check whether the HEALPIX environment variable was set by the user:
# To run Commander3 (and probably other healpix dependent codes), you need to
# export HEALPIX environment variable, which will point to the root folder of
# HEALPIX with all source and compiled codes. Thus, we first looking if such
# variable is defined.
if(NOT HEALPIX_ROOT)
	set(HEALPIX_ROOT "$ENV{HEALPIX}")
endif()
# TODO: 
# - Wrap it around if(Linux) statement as PkgConfig is Linux thing
# HEALPIX comes with PkgConfig file, so, in principle, we can also check
# for existing package on the system installed with PkgConfig. 
find_package(PkgConfig QUIET)
# pkg_search_module, which checks for the package and uses the first available.
# pkg_check_modules, which check for all the corresponding packages.
# Link to pkg.org: https://pkgs.org/search/?q=libhealpix
if(PKG_CONFIG_FOUND)
	pkg_check_modules(PKG_HEALPIX QUIET 
		# Fortran
		healpix
		libhealpix0
		libhealpix-dev
		# C
		chealpix
		# C++
		healpix_cxx
		libhealpix_cxx2 
		libhealpix_cxx-devel
		libhealpix-cxx-dev
		libhealpix-cxx2
		# Java
		#libhealpix-java 
		)
endif()
# Most probably Healpix won't be found with the way above (or it may be, 
# but you will need more libraries etc), so we switch to manual search.
#==============================================================================
# Setting up the list of valid components
list(APPEND _HEALPIX_VALID_COMPONENTS_
	SHARP
	Fortran	
	C
	CXX
	)
# Setting the list of library names to look for
list(APPEND _HEALPIX_LIB_NAMES_
	sharp
	healpix
	chealpix
	healpix_cxx
	)
# Healpix is a collection of different codes
# and not simply one mixture of c, c++, fortran etc. Not all of them (e.g.
# Fortran) have header files => look for other files as well..
list(APPEND _HEALPIX_SHARP_HEADERS_
	sharp.h	
	)
list(APPEND _HEALPIX_Fortran_HEADERS_
	healpix_modules.mod 
	healpix_types.mod 
	healpix_fft.mod
	healpix_sharp_f90.mod	
	)
list(APPEND _HEALPIX_C_HEADERS_
	chealpix.h	
	)	
list(APPEND _HEALPIX_CXX_HEADERS_
	healpix_tables.h
	healpix_base.h
	healpix_data_io.h
	healpix_map.h
	healpix_map_fitsio.h
	)
# HEALPix source contains the data folder as well, which is used
# by various projects and contains pixel window functions etc.
# So, we look for it as well.
list(APPEND _HEALPIX_DATA_FILES_
	pixel_window_n0002.fits
	pixel_window_n0004.fits
	pixel_window_n0008.fits
	pixel_window_n0016.fits
	pixel_window_n0032.fits
	pixel_window_n0064.fits
	pixel_window_n0128.fits
	pixel_window_n0256.fits
	pixel_window_n0512.fits
	pixel_window_n1024.fits
	pixel_window_n2048.fits
	pixel_window_n4096.fits
	pixel_window_n8192.fits
	weight_pixel_n00016.fits
	weight_pixel_n00032.fits
	weight_pixel_n00064.fits
	weight_pixel_n00128.fits
	weight_pixel_n00256.fits
	weight_pixel_n00512.fits
	weight_pixel_n01024.fits
	weight_pixel_n02048.fits
	weight_ring_n00002.fits
	weight_ring_n00004.fits
	weight_ring_n00008.fits
	weight_ring_n00016.fits
	weight_ring_n00032.fits
	weight_ring_n00064.fits
	weight_ring_n00128.fits
	weight_ring_n00256.fits
	weight_ring_n00512.fits
	weight_ring_n01024.fits
	weight_ring_n02048.fits
	weight_ring_n04096.fits
	weight_ring_n08192.fits
	)
# To get HEALPix version
list(APPEND _HEALPIX_VERSION_FILE_
	Version
	)
# Paths were to search for HEALPix components
# TODO: Figure out why this doesn't work, whereas simple listys works
# It finds it only with root and install prefix variables....
list(APPEND _HEALPIX_LIB_PATHS_
	"/usr"
	"/usr/local"
	"~/local"
	"~/.local"
	"${HEALPIX_ROOT}"
	${HEALPIX_INSTALL_PREFIX}
	)
#message(${HEALPIX_ROOT})
#==============================================================================
# Creating two lists - variables and corresponding libraries
# Searching for libraries - looping through the lists above in parallel 
foreach(component lib_name IN ZIP_LISTS _HEALPIX_VALID_COMPONENTS_ _HEALPIX_LIB_NAMES_)
	#message(STATUS "component=HEALPIX_${component}_LIB, name=${lib_name}, path=${_HEALPIX_LIBS_PATHS_}")
	# Include Dirs
	find_path(HEALPIX_${component}_INCLUDE_DIR
		NAMES "${_HEALPIX_${component}_HEADERS_}" 
		PATHS ${HEALPIX_ROOT} ${HEALPIX_INSTALL_PREFIX} #"${_HEALPIX_LIB_PATHS_}"
		PATH_SUFFIXES include include/healpix_cxx include/libsharp
		)
	# Libraries
	# In version 3.70 (and probably earlier), libs for sharp, f90, c, cxx
	# are stored in "lib" folder
	find_library(HEALPIX_${component}_LIB
		NAMES ${lib_name}
		PATHS ${HEALPIX_ROOT} ${HEALPIX_INSTALL_PREFIX} #"${_HEALPIX_LIB_PATHS_}" 
		PATH_SUFFIXES lib lib64 libs .lib .libs .lib64 
		)
	#message("${HEALPIX_ROOT}")
	#message("${_HEALPIX_PATHS_}")
	#message("${_HEALPIX_${component}_HEADERS_}")
endforeach()
# Looking for data files as well
find_path(HEALPIX_DATA_DIR
	NAMES "${_HEALPIX_DATA_FILES_}"
	PATHS ${HEALPIX_ROOT} ${HEALPIX_INSTALL_PREFIX} #"${_HEALPIX_LIB_PATHS_}" 
	PATH_SUFFIXES data
	)
#==============================================================================
# We have a bunch of components which we need to loop through
# Defining boolean variables based on components found
# i.e. HEALPIX_[COMPONENT]_FOUND
#
# Setting HEALPIX_[COMPONENT]_FOUND to FALSE by default
# only after user specifies the component he/she wants to use
# these variables will be set to TRUE.
foreach(component IN LISTS _HEALPIX_VALID_COMPONENTS_)
	set(HEALPIX_${component}_FOUND FALSE)
	#message(STATUS "HEALPIX_${component}_FOUND is ${HEALPIX_${component}_FOUND}")
endforeach()

# We need to loop through each component user specified inside find_package().
# If no component was specified, we set default value to Sharp+Fortran.
# If wrong component was specified (e.g. typo) then will give an error.
if(NOT HEALPIX_FIND_COMPONENTS)
	#message("Debug message")
	if(HEALPIX_Fortran_LIB AND HEALPIX_SHARP_FOUND)
		set(HEALPIX_Fortran_FOUND TRUE)
		set(HEALPIX_LIBRARIES ${HEALPIX_LIBRARIES} ${HEALPIX_SHARP_LIB} ${HEALPIX_Fortran_LIB})
	endif()
else()
	#message(${HEALPIX_FIND_COMPONENTS})
	foreach(component IN LISTS HEALPIX_FIND_COMPONENTS)
		# Returns the index of the element specified in the list or -1 if it wasn’t found.
		list(FIND _HEALPIX_VALID_COMPONENTS_ ${component} component_index)
		if(NOT component_index EQUAL -1)
			#message("Debug message, for component ${component}")
			#message("${HEALPIX_${component}_LIB}")
			# checking whether the library was found
			if(HEALPIX_${component}_LIB)
				#message("My component found")
				set(HEALPIX_${component}_FOUND TRUE)
				list(APPEND HEALPIX_LIBRARIES "${HEALPIX_${component}_LIB}")
			endif()
			# checking whether the headers were found
			if(HEALPIX_${component}_INCLUDE_DIR)
				list(APPEND HEALPIX_INCLUDE_DIRS "${HEALPIX_${component}_INCLUDE_DIR}")
			endif()
		else()
			message(FATAL_ERROR "${component} is not a valid HEALPIX component! Valid are: ${_HEALPIX_VALID_COMPONENTS_}")
		endif()
	endforeach()
endif()
#==============================================================================
# HEALPix comes with the "Version" file, where it places its version. This file 
# is located inside root folder. So, we just need to retrieve it via regular 
# expressions:
if(HEALPIX_FOUND)
	if(NOT HEALPIX_ROOT)
		find_path(_HEALPIX_ROOT_
			NAMES Version
			PATHS ${_HEALPIX_PATHS_}
			)
	else()
		set(_HEALPIX_ROOT_ ${HEALPIX_ROOT})
	endif()

	file(STRINGS "${_HEALPIX_ROOT_}/Version" _VERSION_CONTENT_) 
	foreach(version IN LISTS _VERSION_CONTENT_)
		set(HEALPIX_VERSION ${version})
	endforeach()
	message("${HEALPIX_VERSION}")

	# Checking if version was specified correctly:
	if(HEALPIX_VERSION MATCHES "^([0-9]+)\\.([0-9]+)$")
		message(STATUS "HEALPix version is: ${HEALPIX_VERSION}")
	endif()
endif()

	#message(${HEALPIX_Fortran_FOUND})
#==============================================================================
find_package_handle_standard_args(HEALPIX
	FOUND_VAR			HEALPIX_FOUND 
	REQUIRED_VARS HEALPIX_LIBRARIES 
								HEALPIX_INCLUDE_DIRS
								#HEALPIX_DATA_DIR
	VERSION_VAR		HEALPIX_VERSION 
	HANDLE_COMPONENTS
	)


# Note: component names comes from variables HEALPIX_Fortran_FOUND, e.g.:
# HEALPIX_[COMPONENT]_FOUND
# This option will just hide those variables 
mark_as_advanced(
	# general things
	HEALPIX_LIBRARIES
	HEALPIX_INCLUDE_DIRS
	# Components
	HEALPIX_SHARP_LIB
	HEALPIX_Fortran_LIB
	HEALPIX_C_LIB
	HEALPIX_CXX_LIB
	)

