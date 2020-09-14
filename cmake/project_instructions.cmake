#==============================================================================
# This file contains general instructions how to
# fetch and build the Commander dependencies
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
# setting custom compile/link flags for each release type. These are
# additional flags, which user can define if he/she is not satisfied
# with the existing ones.
set(COMMANDER3_Fortran_COMPILER_FLAGS "" #"${CMAKE_Fortran_FLAGS}"
	CACHE STRING
	"List of all additional flags user wants to add to configuration."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS ""
	CACHE STRING
	"List of additional linker flags user wants to add to configuration."
	)
# setting default compiler flags, but user will be able to overwrite them
# (although it is not really recommended to do so).
set(COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE ""
	CACHE STRING
	"List of compiler flags for Release version."
	)
set(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG ""
	CACHE STRING
	"List of compiler flags for Debug version."
	)
set(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO ""
	CACHE STRING
	"List of compiler flags for RelWithDebInfo version."
	)
set(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL ""
	CACHE STRING
	"List of compiler flags for MinSizeRel version."
	)

# the same as the above but linker flags
set(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE ""
	CACHE STRING
	"List of linker flags for Release version."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS_DEBUG ""
	CACHE STRING
	"List of linker flags for Debug version."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO ""
	CACHE STRING
	"List of linker flags for RelWithDebInfo version."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL ""
	CACHE STRING
	"List of linker flags for MinSizeRel version."
	)
# Commander also has one cpp file, so I add one flag
# similarly to the ones specified in the config files
set(COMMANDER3_CXX_COMPILER_FLAGS "-O3")




#set(COMMANDER3_COMPILER_FLAGS_ADDITIONAL "" #"${CMAKE_Fortran_FLAGS}"
#	CACHE STRING
#	"List of all additional flags user wants to add to configuration."
#	)
#set(COMMANDER3_LINKER_FLAGS_ADDITIONAL ""
#	CACHE STRING
#	"List of additional linker flags user wants to add to configuration."
#	)
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
# Different OS requires different linking
if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	set(COMMANDER3_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
	set(COMMANDER3_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
endif()


# Specifying flags per Fortran compiler
# Intel
# Note: from https://www.osc.edu/documentation/knowledge_base/compilation_guide
# With the Intel compilers, use -xHost and -O2 or higher. 
# With the gnu compilers, use -march=native and -O3. 
# The PGI compilers by default use the highest available instruction set, so no additional flags are necessary.
#------------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	# Compiler flags
	# If user has not specified compilation flag, we use default configuration
	if (COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
			"-Ofast" 
			"-ipo" 
			"-xHost" 
			"-parallel" 
			"-qopenmp" 
			"-qopt-matmul" 
			"-assume" "byterecl" 
			"-heap-arrays" "16384"
			"-fpe0"
			)#"-fast" "-parallel")#"-qopt-matmul" "-heap-arrays 16384 -fpe0 -CB")
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
			"-O0" 
			"-g" 
			"-traceback" 
			"-parallel" 
			"-qopenmp" 
			"-C" 
			"-assume" "byterecl" 
			"-heap-arrays" "16384"
			"-fpe0"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO 
			"-O2" 
			"-g" 
			"-DNDEBUG" 
			"-parallel" 
			"-qopenmp" 
			"-C"
			"-assume" "byterecl" 
			"-heap-arrays" "16384"
			"-fpe0"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL 
			"-Os" 
			"-DNDEBUG" 
			"-parallel" 
			"-qopenmp" 
			"-C"
			"-assume" "byterecl" 
			"-heap-arrays" "16384"
			"-fpe0"
			)
	endif()

	# Linker flags
	# the same logic as with compiler flags
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "-qopt-matmul")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG "")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO "")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL "")
	endif()
#------------------------------------------------------------------------------
# GNU - 9.3 - 10.x needs different flags
# setting different flags for different version
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	if (COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
			"-O3" 
			"-march=native" 
			"-flto" 
			"-fopenmp"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
			"-O0" 
			"-g3" 
			"-Wall" 
			"-Wextra" 
			"-Wconversion" 
			"-C" 
			"-pedantic" 
			"-fbacktrace" 
			"-fcheck=bounds" 
			"-ffpe-trap=zero,overflow,underflow" 
			"-ffunction-sections" 
			"-pipe"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO 
			"-O2" 
			"-g" 
			"-DNDEBUG" 
			"-fopenmp" 
			"-C"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL 
			"-Os" 
			"-DNDEBUG" 
			"-fopenmp" 
			"-C"
			)
	endif()
	# adding different flags depending on the compiler version
	if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL "10")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS "-ffree-line-length-none" "-fallow-argument-mismatch")
	else()
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS "-ffree-line-length-none" "-Wno-argument-mismatch")
	endif()

	# Linker flags
	# the same logic as with compiler flags
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "-flto")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG "")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO "")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL "")
	endif()
#------------------------------------------------------------------------------
# PGI	
# TODO: For some reason commander dependencies crashes on this one
# so figure out how to make it work
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
	# Compiler flags
	if (COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
			"-O4" 
			"-fast" 
			"-mp=all"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
			"-O0" 
			"-g" 
			"-traceback" 
			"-Minfo" 
			"-C"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO 
			"-fast" 
			"-C"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL 
			"-fast" 
			"-C"
			)
	endif()

	# Linker flags
	# the same logic as with compiler flags
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG "")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO "")
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL "")
	endif()
#------------------------------------------------------------------------------
# NVIDIA bought PGI compilers and now they are NVIDIA CUDA compilers
#elseif(CMAKE_Fortran_COMPILER_ID MATCHES NVIDIA)
#	message(STATUS "This is a DEBUG MESSAGE FO NVIDIA COMPILERS")
#------------------------------------------------------------------------------
# Flang
# TODO: need to figure out why healpix doesn't compile with flang
# and then add support for flang
#elseif(CMAKE_Fortran_COMPILER_ID MATCHES Flang)
#	# Compiler flags
#	#list(APPEND COMMANDER3_COMPILER_FLAGS "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL "")
#	# Linker flags
#	#list(APPEND COMMANDER3_LINKER_FLAGS "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL "")
endif()
#------------------------------------------------------------------------------
# Making a summary of compiler location and compile flags
message(STATUS "---------------------------------------------------------------")
message(STATUS "SUMMARY ON COMPILERS:")
message(STATUS "---------------------------------------------------------------")
message(STATUS "Your system is: ${CMAKE_SYSTEM_NAME}, ${CMAKE_SYSTEM}")
message(STATUS "Fortran Compiler is: ${CMAKE_Fortran_COMPILER}")
message(STATUS "C Compiler is: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler is: ${CMAKE_CXX_COMPILER}")
message(STATUS "Commander3 configuration is: ${CMAKE_BUILD_TYPE}. Compiler flags to be applied:")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG};${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_Fortran_COMPILER_FLAGS_ADDITIONAL}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE};${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_COMPILER_FLAGS_ADDITIONAL}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO};${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_COMPILER_FLAGS_ADDITIONAL}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL};${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_COMPILER_FLAGS_ADDITIONAL}")
endif()

# defining the compillation procedure depending on the system
# TODO: do I really need this?
# link to the wiki: https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Write-Platform-Checks
#if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
#	message(STATUS "You seem to be running Linux!!")
#endif()
#message(STATUS "${CMAKE_SYSTEM}")
#message(STATUS "${CMAKE_SYSTEM_NAME}")

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
# setting default values if this variable wasn't defined
set(CMAKE_DOWNLOAD_DIRECTORY "${CMAKE_SOURCE_DIR}/build/downloads"
	CACHE STRING
	"Directory where to download commander dependencies' source files"
	)
# Output dir - originally it is CMAKE_INSTALL_PREFIX, but we need to set a
# different default value and it isn't advisable to temper w/ this variable
#set(CMAKE_INSTALL_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/install"
#	CACHE STRING
#	"Directory where to install commander dependencies"
#	)
# where to output libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib"
	CACHE STRING
	"Directory where to install all the libraries."
	)
# where to output executable(s)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/bin"
	CACHE STRING
	"Directory where to install all the executables."
	)
# HEALPix install (root) dir - by default we will install it in "healpix" dir
set(HEALPIX_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/healpix"
	CACHE STRING
	"Directory where to install (copy compiled and data files of) HEALPix."
	)
set(DOXYGEN_BUILD_DOCS OFF 
	CACHE BOOL
	"Determine whether to use doxygen or not."
	)
#------------------------------------------------------------------------------
# If any problems with installation will occur, which cannot be fixed quickly,
# these variables will force a fresh installation for every specified library.
# forces fresh installation of HDF5 to avoid some errors with old versions
set(HDF5_FORCE_COMPILE FALSE
  CACHE BOOL
  "Forces fresh installation of HDF5."
  )
set(FFTW_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of FFTW."
  )
set(CFITSIO_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of CFITSIO."
  )
set(HEALPIX_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of HEALPIX."
  )
set(CURL_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of CURL."
  )
set(DOXYGEN_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of DOXYGEN."
  )
set(FLEX_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of FLEX."
  )
set(BISON_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of BISON."
  )

#------------------------------------------------------------------------------
# Commander source dir
set(COMMANDER3_SOURCE_DIR "${CMAKE_SOURCE_DIR}/commander3/src")
# tempita source dir
set(TEMPITA_DIR ${CMAKE_SOURCE_DIR}/commander3/python)
# adding custom cmake modules directory, e.g. for FindSomething.cmake
# Note: It should be already inside root CmakeLists.txt, so 
# don't need to include in here
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
#------------------------------------------------------------------------------
# output of the summary into the screen
message(STATUS "---------------------------------------------------------------")
message(STATUS "SUMMARY ON INSTALLATION:")
message(STATUS "---------------------------------------------------------------")
message(STATUS "Projects will be downloaded into: ${CMAKE_DOWNLOAD_DIRECTORY}")
#message(STATUS "Projects will be installed into: ${CMAKE_INSTALL_OUTPUT_DIRECTORY}")
message(STATUS "Projects will be installed into: ${CMAKE_INSTALL_PREFIX}")
#==============================================================================
# MAIN PROJECT DEPENDENCIES
#==============================================================================
message(STATUS "---------------------------------------------------------------")
message(STATUS "Looking for packages...")
# use this to write your own find_package
find_package(PkgConfig)
# We will be using Git to download some dependencies, so we need to check if git available
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
	##sharp2
	fftw
	cfitsio
	hdf5
	doxygen
	healpix
	commander3
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
#set(fftw_f_configure_command 
#	"${CMAKE_COMMAND}" "-E" "env" 
#	"FC=${COMMANDER3_Fortran_COMPILER}" 
#	"CXX=${COMMANDER3_CXX_COMPILER}" 
#	"CPP=${COMMANDER3_CPP_COMPILER}" 
#	"CC=${COMMANDER3_C_COMPILER}" 
#	"MPICC=${COMMANDER3_C_COMPILER}" 
#	"./configure" 
#	"--prefix=<INSTALL_DIR>" 
#	"--enable-float" 
#	"--enable-threads" 
#	"--enable-openmp" 
#	"--enable-mpi"
#	)
# double
#set(fftw_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "MPICC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>" "--enable-threads" "--enable-openmp" "--enable-mpi")
#------------------------------------------------------------------------------
# HDF5
set(hdf5_url "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz")
set(hdf5_md5 "e115eeb66e944fa7814482415dd21cc4")
set(hdf5_configure_command 
	"${CMAKE_COMMAND}" "-E" "env" 
	"FC=${COMMANDER3_Fortran_COMPILER}" 
	"CXX=${COMMANDER3_CXX_COMPILER}" 
	"CPP=${COMMANDER3_CPP_COMPILER}" 
	"CC=${COMMANDER3_C_COMPILER}" 
	"./configure" 
	"--prefix=<INSTALL_DIR>" 
	"--enable-fortran"
	)
#------------------------------------------------------------------------------
# LibSharp2
# Note: libsharp2 comes as a native part of Healpix 3.70 and thus, we do not
# need to compile it independently. But, I will leave this just in case we need
# to switch to older versions (for whatever reason).
set(sharp2_url "https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz")#"https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz") #"https://github.com/Libsharp/libsharp/archive/v1.0.0.tar.gz")
set(SHARP2_C_FLAGS "-DUSE_MPI -std=c99 -O3 -ffast-math")
set(sharp2_configure_command 
	"autoreconf" "-i" "&&" 
	"${CMAKE_COMMAND}" "-E" "env" 
	"FC=${COMMANDER3_Fortran_COMPILER}" 
	"CXX=${COMMANDER3_CXX_COMPILER}" 
	"CPP=${COMMANDER3_CPP_COMPILER}" 
	"CC=${COMMANDER3_C_COMPILER}" 
	"CFLAGS=${SHARP2_C_FLAGS}" 
	"./configure"
	)
#------------------------------------------------------------------------------
# CFitsio
set(cfitsio_url "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz")
set(cfitsio_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure" "--prefix=<INSTALL_DIR>" "--disable-curl")
#------------------------------------------------------------------------------
# HEALPix
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	set(HEALPIX_SHARP2_C_FLAGS "-O3 -ffast-math -march=native -std=c99 -DUSE_MPI -qopenmp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	set(HEALPIX_SHARP2_C_FLAGS "-O3 -ffast-math -march=native -std=c99 -DUSE_MPI -fopenmp")
endif()
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.tar.gz/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.zip/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")#"https://sourceforge.net/projects/healpix/files/latest/download")
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.tar.gz/download")
set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.70/Healpix_3.70_2020Jul23.tar.gz/download")
#set(healpix_md5 "ed7c9a3d7593577628ed1286fa7a9250")
#set(healpix_md5 "540b243406596205a7a82434d99af41e")
#set(healpix_md5 "9b51b2fc919f4e70076d296826eebee0")
set(healpix_md5 "bdcc2a4b1ede3ed5a07be57e4aec01d2")
# this command is for healpix 3.50 and below
#set(healpix_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure")
#set(healpix_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure")
#set(healpix_components
	#profile
	#sharp
	#profile
	#healpy
	#profile
	#f90
	#profile
	#c
	#profile
	#cxx
	#	)
set(healpix_configure_command 
	"${CMAKE_COMMAND}" "-E" "env" 
	"FITSDIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
	"FITSINC=${CMAKE_INSTALL_PREFIX}/include"
	"FC=${COMMANDER3_Fortran_COMPILER}" 
	"CXX=${COMMANDER3_CXX_COMPILER}" 
	"CPP=${COMMANDER3_CPP_COMPILER}" 
	"CC=${COMMANDER3_C_COMPILER}" 
	"SHARP_COPT=${HEALPIX_SHARP2_C_FLAGS}"
	"./configure" 
	"--auto=f90" #${healpix_components}" #profile,f90,c,cxx;" 
	)
#set(healpix_configure_command 
#	"${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER} CXX=${COMMANDER3_CXX_COMPILER} CPP=${COMMANDER3_CPP_COMPILER} CC=${COMMANDER3_C_COMPILER} FITSDIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY} FITSINC=${CMAKE_INSTALL_PREFIX}/include SHARP_COPT=${HEALPIX_SHARP2_C_FLAGS} ./configure --auto=all" 
#	)
#message(STATUS "Healpix configure command is: ${healpix_configure_command}")
#------------------------------------------------------------------------------
# Doxygen
# there is some weird errors appearing for doxygen v1.8.17 and above, so will stick with this one
set(doxygen_url "https://github.com/doxygen/doxygen/archive/Release_1_8_16.tar.gz")#"https://github.com/doxygen/doxygen/archive/Release_1_8_18.tar.gz")#"https://github.com/doxygen/doxygen.git")
# compile doxygen with cmake itself, so no need to specify configure options here
set(flex_url "http://sourceforge.net/projects/flex/files/flex-2.5.39.tar.gz/download")#"https://sourceforge.net/projects/flex/files/flex-2.6.0.tar.xz/download")
set(flex_configure_command 
	"${CMAKE_COMMAND}" "-E" "env" 
	"FC=${COMMANDER3_Fortran_COMPILER}" 
	"CXX=${COMMANDER3_CXX_COMPILER}" 
	"CPP=${COMMANDER3_CPP_COMPILER}" 
	"CC=${COMMANDER3_C_COMPILER}" 
	"./configure" 
	"--prefix=<INSTALL_DIR>"
	)
set(bison_url "http://ftp.gnu.org/gnu/bison/bison-3.6.tar.gz")#"http://ftp.gnu.org/gnu/bison/bison-3.6.2.tar.gz")
set(bison_configure_command 
	"${CMAKE_COMMAND}" "-E" "env" 
	"FC=${COMMANDER3_Fortran_COMPILER}" 
	"CXX=${COMMANDER3_CXX_COMPILER}" 
	"CPP=${COMMANDER3_CPP_COMPILER}" 
	"CC=${COMMANDER3_C_COMPILER}" 
	"./configure" 
	"--prefix=<INSTALL_DIR>"
	)
#------------------------------------------------------------------------------

# include all project configuration files
foreach(project ${projects})
	include("${project}")
	#target_link_libraries(${${project}_lib} ${MPI_Fortran_LIBRARIES})
endforeach()
# 
#include_directories(${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/include)
message(STATUS "---------------------------------------------------------------")
#include_directories(${out_install_dir}/lib)
#include_directories(${out_install_dir}/mod)
# this one here is not advised to use, but desperate times need desperate measures
#link_directories(${out_install_dir}/mod)
