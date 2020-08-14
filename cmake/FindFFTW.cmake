#==============================================================================
# Module to find FFTW3 on the system
# It looks for DOUBLE, FLOAT and LONG components, such as:
# serial, openmp, threads and mpi
# Author: Maksym Brilenkov
#
# This file is  part of Commander. For  details, see <LINK>.  The
# full copyright notice is contained in the file COPYRIGHT located at the root
# of the distribution or at <LICENSE>.
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] [COMPONENTS] [list of all components])
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
# First check whether the FFTW_ROOT variable was set by the user,
# as this variable can be used to guide detection of the FFTW lib
# to a non-standard installation directory
if(NOT FFTW_ROOT)
	set(FFTW_ROOT "$ENV{FFTW_ROOT}")
endif()
# TODO: Wrap it around if(Linux) statement as PkgConfig is Linux thing
# Next, we can also check for a PkgConfig
find_package(PkgConfig QUIET)
# pkg_search_module, which checks for the package and uses the first available.
# pkg_check_modules, which check for all the corresponding packages.
if(PKG_CONFIG_FOUND)
	pkg_check_modules(PKG_FFTW QUIET "fftw3")
endif()
#==============================================================================
# TODO: Add support for different systems, as they have different paths
# Define variables for paths where to find "include" and "lib" folders
list(APPEND FFTW_INCLUDE_DIRS_HINTS 
	"/usr/include" 
	"/usr/local/include" 
	"~/local/include"
	"${FFTW_ROOT}/include"
	"${CMAKE_INSTALL_PREFIX}/include"
	"${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include"
	)
# FFTW MPI has a separate header, so we need to look for that one is well
# Usually, it should be located in the same folder as "fftw3.h", but who knows
list(APPEND FFTW_MPI_INCLUDE_DIRS_HINTS
	"/usr/include" 
	"/usr/local/include" 
	"~/local/include"
	"${FFTW_ROOT}/include"
	"${CMAKE_INSTALL_PREFIX}/include"
	"${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include"
	)
# Specifying the location for FFTW libraries
# Usually libraries compiled as double with other custom specifications
# like --enable-threads etc. But user can recompile it like long or float
# with the same specification, so we need 3 separate paths
# DOUBLE
list(APPEND FFTW_DOUBLE_LIBS_PATHS
	"/usr/lib"
	"/usr/lib64"
	"/usr/local/lib"
	"/usr/local/lib64"
	"~/local/lib"
	"~/local/lib64"	
	"${FFTW_ROOT}/lib"
	"${FFTW_ROOT}/.libs"
	"${CMAKE_INSTALL_PREFIX}/lib"
	"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
	"${PKG_FFTW_LIBRARY_DIRS}"
	)
# FLOAT
list(APPEND FFTW_FLOAT_LIBS_PATHS
	"/usr/lib"
	"/usr/lib64"
	"/usr/local/lib"
	"/usr/local/lib64"
	"~/local/lib"
	"~/local/lib64"	
	"${FFTW_ROOT}/lib"
	"${FFTW_ROOT}/.libs"
	"${CMAKE_INSTALL_PREFIX}/lib"
	"${CMAKE_INSTALL_OUTPUT_DIRECTORY}/lib"
	"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
	"${PKG_FFTW_LIBRARY_DIRS}"
	)
# LONG
list(APPEND FFTW_LONG_LIBS_PATHS
	"/usr/lib"
	"/usr/lib64"
	"/usr/local/lib"
	"/usr/local/lib64"
	"~/local/lib"
	"~/local/lib64"	
	"${FFTW_ROOT}/lib"
	"${FFTW_ROOT}/.libs"
	"${CMAKE_INSTALL_PREFIX}/lib"
	"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
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
# Now, we define FFTW version via string manipulations and regular expressions
#set(_FFTW_H ${FFTW_INCLUDE_DIRS}/fftw3.h)

#function(_fftwver_EXTRACT _FFTW_VER_COMPONENT _FFTW_VER_OUTPUT)
#  set(CMAKE_MATCH_1 "1")
#	set(_FFTW_expr "^[ \\t]*#define[ \\t]+${_FFTW_VER_COMPONENT}[ \\t]+([0-9]+)$")
#	message(${_FFTW_expr})
#	file(STRINGS "${_FFTW_H}" _FFTW_ver REGEX "${_FFTW_expr}")
#	string(REGEX MATCH "${_FFTW_expr}" FFTW_ver "${_FFTW_ver}")
#	set(${_FFTW_VER_OUTPUT} "${CMAKE_MATCH_1}" PARENT_SCOPE)
#	message(${_FFTW_H})
#	message(${_FFTW_VER_OUTPUT})
#endfunction()

#_fftwver_EXTRACT("FFTW_VERSION_MAJOR" FFTW_VERSION_MAJOR)
#_fftwver_EXTRACT("FFTW_VERSION_MINOR" FFTW_VERSION_MINOR)
#_fftwver_EXTRACT("FFTW_VERSION_PATCH" FFTW_VERSION_PATCH)

#message(STATUS "FFTW version: ${FFTW_VERSION_MAJOR}.${FFTW_VERSION_MINOR}.${FFTW_VERSION_PATCH}")

#if(FFTW_FIND_VERSION_COUNT GREATER 2)
#	set(FFTW_VERSION "${FFTW_VERSION_MAJOR}.${FFTW_VERSION_MINOR}.${FFTW_VERSION_PATCH}")
#else()
#	set(FFTW_VERSION "${FFTW_VERSION_MAJOR}.${FFTW_VERSION_MINOR}")
#endif()

#find_library(ZeroMQ_LIBRARIES
#      NAMES
#        libzmq
#        "libzmq-mt-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}"
#        "libzmq-${CMAKE_VS_PLATFORM_TOOLSET}-mt-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}"
#        libzmq_d
#        "libzmq-mt-gd-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}
#“find_library(ZeroMQ_LIBRARIES"
#      NAMES
#        libzmq
#        "libzmq-mt-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}"
#        "libzmq-${CMAKE_VS_PLATFORM_TOOLSET}-mt-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}"
#        libzmq_d
#        "libzmq-mt-gd-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}"
#"libzmq-${CMAKE_VS_PLATFORM_TOOLSET}-mt-gd-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}"
#      HINTS
#       ${_ZeroMQ_ROOT}/lib
#      )
#endif()


#==============================================================================
# Searching for libraries
# DOUBLE
find_library(FFTW_DOUBLE_LIB
	NAMES "fftw3"
	PATHS ${FFTW_DOUBLE_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS ${FFTW_DOUBLE_LIBS_HINTS} 
	)

find_library(FFTW_DOUBLE_THREADS_LIB
	NAMES "fftw3_threads"
	PATHS ${FFTW_DOUBLE_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	HINTS "${FFTW_DOUBLE_LIBS_HINTS}"
	)

find_library(FFTW_DOUBLE_OPENMP_LIB
	NAMES "fftw3_omp"
	PATHS ${FFTW_DOUBLE_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_DOUBLE_LIBS_HINTS}"
	)

find_library(FFTW_DOUBLE_MPI_LIB
	NAMES "fftw3_mpi"
	PATHS ${FFTW_DOUBLE_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_DOUBLE_LIBS_HINTS}"
	)

# FLOAT
find_library(FFTW_FLOAT_LIB
	NAMES "fftw3f" libfftw3f-3
	PATHS ${FFTW_FLOAT_LIBS_PATHS} #${FFTW_ROOT} ${FFTW_FLOAT_LIBS_HINTS}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_FLOAT_LIBS_HINTS}"
	)
find_library(FFTW_FLOAT_THREADS_LIB
	NAMES "fftw3f_threads"
	PATHS ${FFTW_FLOAT_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_FLOAT_LIBS_HINTS}"
	)
find_library(FFTW_FLOAT_OPENMP_LIB
	NAMES "fftw3f_omp"
	PATHS ${FFTW_FLOAT_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_FLOAT_LIBS_HINTS}"
	)
find_library(FFTW_FLOAT_MPI_LIB
	NAMES "fftw3f_mpi"
	PATHS ${FFTW_FLOAT_LIBS_PATHS} #${FFTW_ROOT} ${FFTW_FLOAT_LIBS_HINTS}
	#PATH_SUFFIXES lib lib64 
	#HINTS ${FFTW_FLOAT_LIBS_PATHS} #"${FFTW_FLOAT_LIBS_HINTS}"
	)

# LONG DOUBLE
find_library(FFTW_LONG_LIB
	NAMES "fftw3l"
	PATHS ${FFTW_LONG_LIBS_PATHS} # ${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_LONG_LIBS_HINTS}"
	)

find_library(FFTW_LONG_THREADS_LIB
	NAMES "fftw3l_threads"
	PATHS ${FFTW_LONG_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS ${FFTW_LONG_LIBS_PATHS} #"${FFTW_LONG_LIBS_HINTS}"
	)

find_library(FFTW_LONG_OPENMP_LIB
	NAMES "fftw3l_omp"
	PATHS ${FFTW_LONG_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_LONG_LIBS_HINTS}"
	)

find_library(FFTW_LONG_MPI_LIB
	NAMES "fftw3l_mpi"
	PATHS ${FFTW_LONG_LIBS_PATHS} #${FFTW_ROOT}
	#PATH_SUFFIXES lib lib64 
	#HINTS "${FFTW_LONG_LIBS_HINTS}"
	)
#==============================================================================
# We have a bunch of components which we need to loop through
# Defining boolean variables based on components found
# i.e. FFTW_[COMPONENT]_FOUND
#==============================================================================
# Setting up the list of valid components
set(FFTW_VALID_COMPONENTS 
	DOUBLE
	DOUBLE_THREADS
	DOUBLE_OPENMP
	DOUBLE_MPI
	FLOAT
	FLOAT_THREADS
	FLOAT_OPENMP
	FLOAT_MPI
	LONG
	LONG_THREADS
	LONG_OPENMP
	LONG_MPI
	)
# Setting FFTW_[COMPONENT]_FOUND to FALSE by default
# only after user specifies the component he/she wants to use
# these variables will be set to TRUE.
foreach(component IN LISTS FFTW_VALID_COMPONENTS)
	set(FFTW_${component}_FOUND FALSE)
	#message("FFTW_${component}_FOUND is ${FFTW_${component}_FOUND}")
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
		list(FIND FFTW_VALID_COMPONENTS ${component} component_index)
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

find_package_handle_standard_args(FFTW
	FOUND_VAR			FFTW_FOUND 
	REQUIRED_VARS FFTW_INCLUDE_DIRS
	VERSION_VAR		FFTW_VERSION 
	HANDLE_COMPONENTS
	)


# Note: component names comes from variables FFTW_DOUBLE_THREADS_FOUND, e.g.:
# FFTW_[COMPONENT]_FOUND
# This option will just hide those variables 
mark_as_advanced(
	# general things
	FFTW_LIBRARIES
	FFTW_INCLUDE_DIRS
	FFTW_MPI_INCLUDE_DIRS
	# Components
	# DOUBLE
	FFTW_DOUBLE_LIB
	FFTW_DOUBLE_THREADS_LIB
	FFTW_DOUBLE_OPENMP_LIB
	FFTW_DOUBLE_MPI_LIB
	# FLOAT
	FFTW_FLOAT_LIB
	FFTW_FLOAT_THREADS_LIB
	FFTW_FLOAT_OPENMP_LIB
	FFTW_FLOAT_MPI_LIB
	# LONG
	FFTW_LONG_LIB
	FFTW_LONG_THREADS_LIB
	FFTW_LONG_OPENMP_LIB
	FFTW_LONG_MPI_LIB
	)

