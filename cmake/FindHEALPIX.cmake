#==============================================================================
# Module to find HEALPix on the system
# It looks for DOUBLE, FLOAT and LONG components, such as:
# serial, openmp, threads and mpi
# Author: Maksym Brilenkov
#
# This file is part of Commander. For  details, see <LINK>.  The
# full copyright notice is contained in the file COPYRIGHT located at the root
# of the distribution or at <LICENSE>.
#
# Usage:
#   find_package(HEALPIX [version number] [REQUIRED] [QUIET] [COMPONENTS] [list of all components])
#
# Supported components are:
# - Float:       FLOAT, FLOAT_THREADS, FLOAT_OPENMP, FLOAT_MPI
# - Double:      DOUBLE, DOUBLE_THREADS, DOUBLE_OPENMP, DOUBLE_MPI
# - Long Double: LONG, LONG_THREADS, LONG_OPENMP, LONG_MPI
#
# The following variables are defined:
#   FFTW_FOUND                  set to TRUE if FFTW is found on the system;
#   FFTW_[COMPONENT]_FOUND      set to TRUE if component is found on the system;
#   FFTW_LIBRARIES              complete list (with all found components) of FFTW libs' directory paths;
#   FFTW_[COMPONENT]_LIB        full path to one of the components;
#   FFTW_INCLUDE_DIRS           FFTW include directory paths; 
#   FFTW_MPI_INCLUDE_DIRS       FFTW MPI include directory paths; 
#==============================================================================
# Reference to CMake docs:
# https://cmake.org/cmake/help/v3.17/module/FindPackageHandleStandardArgs.html
include(FindPackageHandleStandardArgs)
include(SelectLibraryConfigurations)
# setting the list of valid components
#==============================================================================
# First check whether the HEALPIX environment variable was set by the user:
# To run commander (and probably other healpix dependent codes), you need to
# export HEALPIX environment variable, which will point to the root folder of
# HEALPIX with all source and compiled codes. Thus, we first looking if such
# variable is defined.
if(NOT HEALPIX_ROOT)
	set(HEALPIX_ROOT "$ENV{HEALPIX}")
endif()
# TODO: 
# - Figure out how make it work;
# - Wrap it around if(Linux) statement as PkgConfig is Linux thing
# HEALPIX comes with PkgConfig file, so, in principle, we can also check
# for existing  package on the system installed with PkgConfig. 
find_package(PkgConfig QUIET)
# pkg_search_module, which checks for the package and uses the first available.
# pkg_check_modules, which check for all the corresponding packages.
# Link to pkg.org: https://pkgs.org/search/?q=libhealpix
if(PKG_CONFIG_FOUND)
	pkg_check_modules(PKG_HEALPIX QUIET 
		# Fortran
		libhealpix0
		libhealpix-dev
		# C++
		libhealpix_cxx2 
		libhealpix_cxx-devel
		libhealpix-cxx-dev
		libhealpix-cxx2
		# Java
		libhealpix-java 
		)
endif()
# Most probably Healpix won't be found with the way above (or it may be, 
# but you will need more libraries etc), so we switch to manual search.
#==============================================================================
# HEALPix comes with the "Version" file, where it places its version. We just
# need to retrieve it via regular expressions:

# Checking if version was specified correctly:
if(HEALPIX_VERSION MATCHES "^([0-9]+)\\.([0-9]+)$")
	message(STATUS "HEALPix version is: ${HEALPIX_VERSION}")
endif()
#==============================================================================
# TODO: Add support for different systems, as they have different paths
# Define variables for paths where to find "include" and "lib" folders
list(APPEND HEALPIX_INCLUDE_DIRS_HINTS 
	"/usr" 
	"/usr/healpix" 
	"/usr/local" 
	"~/local/include"
	"${FFTW_ROOT}/include"
	"${CMAKE_INSTALL_PREFIX}/include"
	#"${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include"
	)
# FFTW MPI has a separate header, so we need to look for that one as well
# Usually, it should be located in the same folder as "fftw3.h", but who knows
#list(APPEND FFTW_MPI_INCLUDE_DIRS_HINTS
#	"/usr/include" 
#	"/usr/local/include" 
#	"~/local/include"
#	"${FFTW_ROOT}/include"
#	"${CMAKE_INSTALL_PREFIX}/include"
#	#"${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include"
#	)
# Specifying the location for FFTW libraries
# Usually libraries compiled as double with other custom specifications
# like --enable-threads etc. But user can recompile it like long or float
# with the same specification. In addition, these can be placed in 3 
# separate locations (although this is highly unlikely, but ...). Thus,
# we need 3 separate paths.
# TODO: Figure out whether we really need 3 separate paths, as right now
# you simply use the same paths 3 times.
# DOUBLE
#list(APPEND FFTW_DOUBLE_LIBS_PATHS
#	"/usr/lib"
#	#"/usr/lib64"
#	"/usr/local"
#	#"/usr/local/lib64"
#	"~/local"
#	#"~/local/lib64"	
#	"${FFTW_ROOT}"
#	#"${FFTW_ROOT}/.libs"
#	"${CMAKE_INSTALL_PREFIX}"
#	#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
#	"${PKG_FFTW_LIBRARY_DIRS}"
#	)
# FLOAT
#list(APPEND FFTW_FLOAT_LIBS_PATHS
#	"/usr"
#	#"/usr/lib64"
#	"/usr/local"
#	#"/usr/local/lib64"
#	"~/local"
#	#"~/local/lib64"	
#	"${FFTW_ROOT}"
#	#"${FFTW_ROOT}/.libs"
#	"${CMAKE_INSTALL_PREFIX}"
#	#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
#	"${PKG_FFTW_LIBRARY_DIRS}"
#	)
# LONG
list(APPEND FFTW_LONG_LIBS_PATHS
	"/usr"
	#"/usr/lib64"
	"/usr/local"
	#"/usr/local/lib64"
	"~/local"
	#"~/local/lib64"	
	"${FFTW_ROOT}"
	#"${FFTW_ROOT}/.libs"
	"${CMAKE_INSTALL_PREFIX}"
	#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
	"${PKG_FFTW_LIBRARY_DIRS}"
	)
#==============================================================================
# Search for fftw3.h, fftw3-mpi.h files on the system. If header files were 
# successfully found, FFTW_INCLUDE_DIRS and FFTW_MPI_INCLUDE_DIRS are set to
# appropriate locations. 
if(NOT FFTW_ROOT)
	find_path(FFTW_INCLUDE_DIRS 
		NAMES fftw3.h fftw3-mpi.h
		PATHS ${FFTW_ROOT}
		#PATH_SUFFIXES include
		HINTS ${FFTW_INCLUDE_DIRS_HINTS} 
		)
	find_path(FFTW_MPI_INCLUDE_DIRS
		NAMES fftw3-mpi.h
		PATHS ${FFTW_ROOT}
		HINTS ${FFTW_MPI_INCLUDE_DIRS_HINTS}
		)
else()
	set(FFTW_INCLUDE_DIRS			"${FFTW_ROOT}")
	set(FFTW_MPI_INCLUDE_DIRS "${FFTW_ROOT}")
endif()

# TODO: try to figure out why it doesn't work, i.e. how to find version 3.3.8 
# Version information should be provided inside Header files, but it is not for FFTW
# nothing we can do about it :( 

#==============================================================================
# Creating two lists - variables and corresponding libraries
# Setting up the list of valid components
list(APPEND _FFTW_VALID_COMPONENTS_
	# Double
	DOUBLE
	DOUBLE_THREADS
	DOUBLE_OPENMP
	DOUBLE_MPI
	# Float
	FLOAT
	FLOAT_THREADS
	FLOAT_OPENMP
	FLOAT_MPI
	# Long Double
	LONG
	LONG_THREADS
	LONG_OPENMP
	LONG_MPI
	)
list(APPEND _FFTW_LIB_NAMES_
	# Double
	""
	"_threads"	
	"_omp"
	"_mpi"
	# Float
	"f"
	"f_threads"	
	"f_omp"
	"f_mpi"
	# Long Double
	"l"
	"l_threads"	
	"l_omp"
	"l_mpi"
	)
list(APPEND _FFTW_LIB_PATHS_
	# Double
	_DOUBLE_LIBS_
	_DOUBLE_LIBS_
	_DOUBLE_LIBS_
	_DOUBLE_LIBS_
	# Float
	_FLOAT_LIBS_
	_FLOAT_LIBS_
	_FLOAT_LIBS_
	_FLOAT_LIBS_
	# long Double
	_LONG_LIBS_
	_LONG_LIBS_
	_LONG_LIBS_
	_LONG_LIBS_
	)
# Searching for libraries - looping through the lists above in parallel 
foreach(_component_ _name_ _path_ IN ZIP_LISTS _FFTW_VALID_COMPONENTS_ _FFTW_LIB_NAMES_ _FFTW_LIB_PATHS_)
	#message(STATUS "component=FFTW_${_component_}_LIB, name=fftw3${_name_}, path=FFTW${_path_}PATHS")
	#message(STATUS "${FFTW${_path_}PATHS}")
	#message(STATUS "${PKG_FFTW_LIBRARY_DIRS}")
	find_library(FFTW_${_component_}_LIB
		NAMES "fftw3${_name_}"
		PATHS "${FFTW${_path_}PATHS}" 
		#HINTS ${FFTW${_path_}PATHS} 
		PATH_SUFFIXES lib lib64 libs .lib .libs .lib64 
		#HINTS ${FFTW_DOUBLE_LIBS_HINTS} 
		)
endforeach()
#==============================================================================
# We have a bunch of components which we need to loop through
# Defining boolean variables based on components found
# i.e. FFTW_[COMPONENT]_FOUND
#==============================================================================
# Setting FFTW_[COMPONENT]_FOUND to FALSE by default
# only after user specifies the component he/she wants to use
# these variables will be set to TRUE.
foreach(component IN LISTS _FFTW_VALID_COMPONENTS_)
	set(FFTW_${component}_FOUND FALSE)
	#message(STATUS "FFTW_${component}_FOUND is ${FFTW_${component}_FOUND}")
endforeach()

# We need to loop thrpugh each component user specified inside find_package().
# If no component was specified, we set default value to DOUBLE.
# If wrong component was specified (e.g. typo) then will give an error.
if(NOT FFTW_FIND_COMPONENTS)
	if(FFTW_DOUBLE_LIB)
		set(FFTW_DOUBLE_FOUND TRUE)
		set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_LIB})
	endif()
else()
	foreach(component IN LISTS FFTW_FIND_COMPONENTS)
		# Returns the index of the element specified in the list or -1 if it wasn’t found.
		list(FIND _FFTW_VALID_COMPONENTS_ ${component} component_index)
		if(NOT component_index EQUAL -1)
			# checking whether the library was found
			if(FFTW_${component}_LIB)
				set(FFTW_${component}_FOUND TRUE)
				list(APPEND FFTW_LIBRARIES "${FFTW_${component}_LIB}")
			endif()
		else()
			message(FATAL_ERROR "${component} is not a valid FFTW component! Valid are: ${FFTW_VALID_COMPONENTS}")
		endif()
	endforeach()
endif()

find_package_handle_standard_args(HEALPIX
	FOUND_VAR			HEALPIX_FOUND 
	REQUIRED_VARS HEALPIX_LIBRARIES #FFTW_INCLUDE_DIRS
	VERSION_VAR		HEALPIX_VERSION 
	HANDLE_COMPONENTS
	)


# Note: component names comes from variables HEALPIX_DOUBLE_THREADS_FOUND, e.g.:
# HEALPIX_[COMPONENT]_FOUND
# This option will just hide those variables 
#mark_as_advanced(
#	# general things
#	FFTW_LIBRARIES
#	FFTW_INCLUDE_DIRS
#	FFTW_MPI_INCLUDE_DIRS
#	# Components
#	# DOUBLE
#	FFTW_DOUBLE_LIB
#	FFTW_DOUBLE_THREADS_LIB
#	FFTW_DOUBLE_OPENMP_LIB
#	FFTW_DOUBLE_MPI_LIB
#	# FLOAT
#	FFTW_FLOAT_LIB
#	FFTW_FLOAT_THREADS_LIB
#	FFTW_FLOAT_OPENMP_LIB
#	FFTW_FLOAT_MPI_LIB
#	# LONG
#	FFTW_LONG_LIB
#	FFTW_LONG_THREADS_LIB
#	FFTW_LONG_OPENMP_LIB
#	FFTW_LONG_MPI_LIB
#	)

