#==============================================================================
# This file contains general instructions how to
# fetch and build the commander dependencies
# Author: Maksym Brilenkov
#==============================================================================


#==============================================================================
# COMPILE AND LINK OPTIONS 
#==============================================================================
# Setting project configuration -- variable is empty by default
# - "Debug" builds library/executable w/o optimization and w/ debug symbols; 
# - "Release" builds library/executable w/ optimization and w/o debug symbols;
# - "RelWithDebInfo" builds library/executable w/ less aggressive optimizations and w/ debug symbols;
# - "MinSizeRel" builds library/executable w/ optimizations that do not increase object code size. 
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release
		CACHE STRING
		"Specifies the Build type. Available options are: Release, Debug, RelWithDebInfo, MinSizeRel. Default: Release." FORCE)
endif()
# Currently supported Compilers (CMake v3.17).
# More details in: 
# - AppleClang: Apple Clang for Xcode versions 4.4+.
# - Clang: Clang compiler versions 2.9+.
# - GNU: GNU compiler versions 4.4+.
# - MSVC: Microsoft Visual Studio versions 2010+.
# - SunPro: Oracle SolarisStudio versions 12.4+.
# - Intel: Intel compiler versions 12.1+.
# - Cray: Cray Compiler Environment version 8.1+.
# - PGI: PGI version 12.10+.
# - XL: IBM XL version 10.1+.
# setting custom compile/link flags for each release type 
list(APPEND COMMANDER3_COMPILER_FLAGS "${CMAKE_Fortran_FLAGS}")
# for some reason it doesn't want to work with variables, 
# so need to specify flags manually 
#list(APPEND COMMANDER3_COMPILER_FLAGS_RELEASE "")#"-O3" "-DNDEBUG")#"${CMAKE_Fortran_FLAGS_RELEASE}")
#list(APPEND COMMANDER3_COMPILER_FLAGS_DEBUG "-g")#"${CMAKE_Fortran_FLAGS_DEBUG}")
#list(APPEND COMMANDER3_COMPILER_FLAGS_RELWITHDEBINFO "-O2" "-g" "-DNDEBUG")#"${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
#list(APPEND COMMANDER3_COMPILER_FLAGS_MINSIZEREL "-Os" "-DNDEBUG")#"${CMAKE_Fortran_FLAGS_MINSIZEREL}")

# Specific flags for each subproject
# According to installation guidelines, we need to 
# specify slightly different flags for different compilers
#list(APPEND SHARP2_C_FLAGS "-DUSE_MPI" "-O3" "-std=c99" "-ffast-math")
# For easier debugging, I have created compiler variables
set(COMMANDER3_C_COMPILER "${MPI_C_COMPILER}")
set(COMMANDER3_CXX_COMPILER "${MPI_CXX_COMPILER}")
set(COMMANDER3_Fortran_COMPILER "${MPI_Fortran_COMPILER}")
# To get rid of "C preprocessor fails sanity check" error
# we need to link the CPP as follows for each subproject
set(COMMANDER3_CPP_COMPILER "${MPI_CXX_COMPILER} -E")

#message(${COMMANDER3_Fortran_COMPILER})

# Specifying flags per Fortran compiler
# Intel
# Note: from https://www.osc.edu/documentation/knowledge_base/compilation_guide
# With the Intel compilers, use -xHost and -O2 or higher. 
# With the gnu compilers, use -march=native and -O3. 
# The PGI compilers by default use the highest available instruction set, so no additional flags are necessary.
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	# Compiler flags
	list(APPEND COMMANDER3_COMPILER_FLAGS "")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELEASE "-Ofast" "-ipo" "-xHost" "-parallel" "-qopenmp" "-qopt-matmul")#"-fast" "-parallel")#"-qopt-matmul" "-heap-arrays 16384 -fpe0 -CB")
	list(APPEND COMMANDER3_COMPILER_FLAGS_DEBUG "-O0" "-g" "-traceback")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELWITHDEBINFO "-O2" "-g" "-DNDEBUG")
	list(APPEND COMMANDER3_COMPILER_FLAGS_MINSIZEREL "-Os" "-DNDEBUG")
	# Linker flags
	list(APPEND COMMANDER3_LINKER_FLAGS "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELEASE "-qopt-matmul")
	list(APPEND COMMANDER3_LINKER_FLAGS_DEBUG "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELWITHDEBINFO "")
	list(APPEND COMMANDER3_LINKER_FLAGS_MINSIZEREL "")
# GNU - 9.3 - 10.x needs different flags
# setting different flags for different version
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	#message(STATUS "${CMAKE_Fortran_COMPILER_VERSION}")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELEASE "-O3" "-march=native" "-flto" "-fopenmp")
	list(APPEND COMMANDER3_COMPILER_FLAGS_DEBUG "-O0" "-g3" "-Wall" "-Wextra" "-Wconversion" "-pedantic" "-fbacktrace" "-fcheck=bounds" "-ffpe-trap=zero,overflow,underflow" "-ffunction-sections" "-pipe")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELWITHDEBINFO "-O2" "-g" "-DNDEBUG" "-fopenmp")
	list(APPEND COMMANDER3_COMPILER_FLAGS_MINSIZEREL "-Os" "-DNDEBUG" "-fopenmp")
	# adding different flags depending on the compiler version
	if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL "10")
		list(APPEND COMMANDER3_COMPILER_FLAGS "-ffree-line-length-none" "-fallow-argument-mismatch")
	else()
		list(APPEND COMMANDER3_COMPILER_FLAGS "-ffree-line-length-none" "-Wno-argument-mismatch")
	endif()
	# Linker flags
	list(APPEND COMMANDER3_LINKER_FLAGS "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELEASE "-flto")
	list(APPEND COMMANDER3_LINKER_FLAGS_DEBUG "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELWITHDEBINFO "")
	list(APPEND COMMANDER3_LINKER_FLAGS_MINSIZEREL "")
# PGI	
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
	# Compiler flags
	list(APPEND COMMANDER3_COMPILER_FLAGS "")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELEASE "-O4" "-fast" "-mp=all")
	list(APPEND COMMANDER3_COMPILER_FLAGS_DEBUG "-O0" "-g" "-traceback" "-Minfo")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELWITHDEBINFO "-fast")
	list(APPEND COMMANDER3_COMPILER_FLAGS_MINSIZEREL "-fast")
	# Linker flags
	list(APPEND COMMANDER3_LINKER_FLAGS "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELEASE "")
	list(APPEND COMMANDER3_LINKER_FLAGS_DEBUG "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELWITHDEBINFO "")
	list(APPEND COMMANDER3_LINKER_FLAGS_MINSIZEREL "")
# Flang
# TODO: need to figure out why healpix doesn't compile with flang
# and then add support for flang
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Flang)
	# Compiler flags
	list(APPEND COMMANDER3_COMPILER_FLAGS "")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELEASE "")
	list(APPEND COMMANDER3_COMPILER_FLAGS_DEBUG "")
	list(APPEND COMMANDER3_COMPILER_FLAGS_RELWITHDEBINFO "")
	list(APPEND COMMANDER3_COMPILER_FLAGS_MINSIZEREL "")
	# Linker flags
	list(APPEND COMMANDER3_LINKER_FLAGS "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELEASE "")
	list(APPEND COMMANDER3_LINKER_FLAGS_DEBUG "")
	list(APPEND COMMANDER3_LINKER_FLAGS_RELWITHDEBINFO "")
	list(APPEND COMMANDER3_LINKER_FLAGS_MINSIZEREL "")
endif()
# Making a summary of compiler location and compile flags
message(STATUS "---------------------------------------------------------------")
message(STATUS "SUMMARY ON COMPILERS:")
message(STATUS "---------------------------------------------------------------")
message(STATUS "Fortran Compiler is: ${CMAKE_Fortran_COMPILER}")
message(STATUS "C Compiler is: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler is: ${CMAKE_CXX_COMPILER}")
message(STATUS "Commander3 configuration is: ${CMAKE_BUILD_TYPE}. Compiler flags to be applied:")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	message(STATUS "${COMMANDER3_COMPILER_FLAGS_DEBUG};${COMMANDER3_COMPILER_FLAGS}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
	message(STATUS "${COMMANDER3_COMPILER_FLAGS_RELEASE};${COMMANDER3_COMPILER_FLAGS}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
	message(STATUS "${COMMANDER3_COMPILER_FLAGS_RELWITHDEBINFO};${COMMANDER3_COMPILER_FLAGS}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel")
	message(STATUS "${COMMANDER3_COMPILER_FLAGS_MINSIZEREL};${COMMANDER3_COMPILER_FLAGS}")
endif()

# defining the compillation procedure depending on the system
# TODO: do I really need this?
# link to the wiki: https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Write-Platform-Checks
if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	message(STATUS "You seem to be running Linux!!")
endif()
message(STATUS "${CMAKE_SYSTEM}")
message(STATUS "${CMAKE_SYSTEM_NAME}")

#if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
#  list(APPEND CXX_FLAGS "-fno-rtti" "-fno-exceptions")
#  list(APPEND CXX_FLAGS_DEBUG "-Wsuggest-final-types" "-Wsuggest-final-methods" "-Wsuggest-override")
#  list(APPEND CXX_FLAGS_RELEASE "-O3" "-Wno-unused")
#endif()
#if(CMAKE_Fortran_COMPILER_ID MATCHES Clang)
#  list(APPEND CXX_FLAGS "-fno-rtti" "-fno-exceptions" "-Qunused-arguments" "-fcolor-diagnostics")
#  list(APPEND CXX_FLAGS_DEBUG "-Wdocumentation")
#  list(APPEND CXX_FLAGS_RELEASE "-O3" "-Wno-unused")
#endif()
#==============================================================================
# PROJECT'S DOWNLOAD/INSTALL/OUTPUT DIRECTORIES 
#==============================================================================
# Download dir
set(CMAKE_DOWNLOAD_DIRECTORY "${CMAKE_SOURCE_DIR}/build/downloads"
	CACHE STRING
	"Directory where to download commander dependencies' source files")
# Output dir - originally it is CMAKE_INSTALL_PREFIX, but we need to set a
# different default value and it isn't advisable to temper w/ this variable
set(CMAKE_INSTALL_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/install"
	CACHE STRING
	"Directory where to install commander dependencies")
# where to output libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_OUTPUT_DIRECTORY}/lib")
# where to output executable(s)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_OUTPUT_DIRECTORY}/bin")
set(DOXYGEN_BUILD_DOCS OFF 
	CACHE BOOL
	"Determine whether to use doxygen or not")
# Commander source dir
set(COMMANDER3_SOURCE_DIR "${CMAKE_SOURCE_DIR}/commander3/src")
# tempita source dir
set(TEMPITA_DIR ${CMAKE_SOURCE_DIR}/commander3/python)
# adding custom cmake modules directory, e.g. for FindSomething.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
# output of the summary into the screen
message(STATUS "---------------------------------------------------------------")
message(STATUS "SUMMARY ON INSTALLATION:")
message(STATUS "---------------------------------------------------------------")
message(STATUS "Projects will be downloaded into: ${CMAKE_DOWNLOAD_DIRECTORY}")
message(STATUS "Projects will be installed into: ${CMAKE_INSTALL_OUTPUT_DIRECTORY}")
#==============================================================================
# MAIN PROJECT DEPENDENCIES
#==============================================================================
message(STATUS "---------------------------------------------------------------")
message(STATUS "Looking for packages...")
# use this to write your own find_package
find_package(PkgConfig)
# We will be using Git to download some dependencies, so we need to check if git available
#include(FindGit)
find_package(Git REQUIRED)
# finding math library
message(STATUS "math (m) libraries are: ${LIBM_LIBRARIES}")
# printing out the dl libs, which are also required on some unix systems
message(STATUS "dl libs are: ${CMAKE_DL_LIBS}")

unset(projects)
# project names
list(APPEND projects 
	tempita
	blas # blas-lapack module 
	mpi
	openmp
	curl
	zlib
	sharp2
	fftw
	cfitsio
	hdf5
	doxygen
	healpix
	)
#==============================================================================
# PROJECTS' URL SOURCES, MD5 HASHES AND CONFIGURE COMMANDS
#==============================================================================
# cURL - needed by CFitsio and HEALPix
# need to specify command separately, othewise it won't work
set(curl_url "https://github.com/curl/curl/releases/download/curl-7_69_0/curl-7.69.0.zip")#"https://github.com/curl/curl/releases/download/curl-7_69_1/curl-7.69.1.tar.gz")
set(curl_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>")
#------------------------------------------------------------------------------
# FFTW
set(fftw_url "http://fftw.org/fftw-3.3.8.tar.gz")
set(fftw_md5 "8aac833c943d8e90d51b697b27d4384d")
# we need to compile fftw twice
# float
set(fftw_f_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>" "--enable-float" "--enable-threads" "--enable-openmp")
# default
set(fftw_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>" "--enable-threads" "--enable-openmp")
#------------------------------------------------------------------------------
# HDF5
set(hdf5_url "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz")
set(hdf5_md5 "e115eeb66e944fa7814482415dd21cc4")
set(hdf5_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>" "--enable-fortran")
#------------------------------------------------------------------------------
# LibSharp2
set(sharp2_url "https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz")#"https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz") #"https://github.com/Libsharp/libsharp/archive/v1.0.0.tar.gz")
set(SHARP2_C_FLAGS "-DUSE_MPI -std=c99 -O3 -ffast-math")
set(sharp2_configure_command "autoreconf" "-i" "&&" "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "CFLAGS=${SHARP2_C_FLAGS}" "./configure")
#------------------------------------------------------------------------------
# CFitsio
set(cfitsio_url "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz")
set(cfitsio_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>" "--disable-curl")
#------------------------------------------------------------------------------
# HEALPix
set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.tar.gz/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.zip/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")#"https://sourceforge.net/projects/healpix/files/latest/download")
set(healpix_md5 "ed7c9a3d7593577628ed1286fa7a9250")
set(healpix_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure")
#------------------------------------------------------------------------------
# Doxygen
# there is some weird errors appearing for doxygen v1.8.17 and above, so will stick with this one
set(doxygen_url "https://github.com/doxygen/doxygen/archive/Release_1_8_16.tar.gz")#"https://github.com/doxygen/doxygen/archive/Release_1_8_18.tar.gz")#"https://github.com/doxygen/doxygen.git")
# compile doxygen with cmake itself, so no need to specify configure options here
set(flex_url "http://sourceforge.net/projects/flex/files/flex-2.5.39.tar.gz/download")#"https://sourceforge.net/projects/flex/files/flex-2.6.0.tar.xz/download")
set(flex_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>")
set(bison_url "http://ftp.gnu.org/gnu/bison/bison-3.6.tar.gz")#"http://ftp.gnu.org/gnu/bison/bison-3.6.2.tar.gz")
set(bison_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>")
#------------------------------------------------------------------------------

# projects configure commands (by default it is assumes to be cmake)
#set(tempita_configure_command "")

# projects install commands
#set(fftw_install_command "make install")
##set(hdf5_install_command )
##set(sharp2_install_command "")
#set(cfitsio_install_command )

# include all project configuration files
foreach(project ${projects})
	include("${project}")
	#target_link_libraries(${${project}_lib} ${MPI_Fortran_LIBRARIES})
endforeach()
# 
include_directories(${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include)
message(STATUS "---------------------------------------------------------------")
#include_directories(${out_install_dir}/lib)
#include_directories(${out_install_dir}/mod)
# this one here is not advised to use, but desperate times need desperate measures
#link_directories(${out_install_dir}/mod)
