# Project: FFTW
# File which contains setup for current project 
# Author: Maksym Brilenkov

# To get rid of "C preprocessor fails sanity check" error
# we need to link the CPP as follows
#set(FFTW_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
message(STATUS "---------------------------------------------------------------")
# TODO: make it so components will matter because now it install everything because 
# I gave the command to add appropriate configure suboptions to configure command
if(NOT FFTW_FORCE_COMPILE)
	find_package(FFTW 
		COMPONENTS 
		DOUBLE 
		FLOAT 
		#FLOAT_MPI 
		FLOAT_OPENMP
		FLOAT_THREADS
		)
endif()
# Is TRUE if one of the components were missing
if(NOT FFTW_FOUND)
	# Configure command to compile FFTW from source
	# DOUBLE (default one)
	list(APPEND fftw_double_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${COMMANDER3_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${COMMANDER3_C_COMPILER}" 
		"MPICC=${COMMANDER3_C_COMPILER}" 
		"./configure" 
		"--prefix=<INSTALL_DIR>")
	# FLOAT
	list(APPEND fftw_float_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${COMMANDER3_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${COMMANDER3_C_COMPILER}" 
		"MPICC=${COMMANDER3_C_COMPILER}" 
		"./configure" 
		"--prefix=<INSTALL_DIR>")
	# First, we determine which component were found, so we can link them
	# others will be compiled from source
	# double component is default one, so no configuration option required
	if(NOT FFTW_DOUBLE_FOUND)
		message(STATUS "Missing component - DOUBLE - will be compiled from source")	
	else()
		message(STATUS "Found FFTW_DOUBLE_LIB: ${FFTW_DOUBLE_LIB}")
	endif()
	if(NOT FFTW_FLOAT_FOUND)
		message(STATUS "Missing component - FLOAT - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-float")
	else()
		message(STATUS "Found FFTW_FLOAT_LIB: ${FFTW_FLOAT_LIB}")
	endif()
	if(NOT FFTW_FLOAT_THREADS_FOUND)
		message(STATUS "Missing component - FLOAT_THREADS - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-threads")
	else()
		message(STATUS "Found FFTW_FLOAT_THREADS_LIB: ${FFTW_FLOAT_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_OPENMP_FOUND)
		message(STATUS "Missing component - FLOAT_OPENMP - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-openmp")
	else()
		message(STATUS "Found FFTW_FLOAT_OPENMP_LIB: ${FFTW_FLOAT_OPENMP_LIB}")
	endif()
	if(NOT FFTW_FLOAT_MPI_FOUND)
		message(STATUS "Missing component - FLOAT_MPI - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-mpi")
	else()
		message(STATUS "Found FFTW_FLOAT_MPI_LIB: ${FFTW_FLOAT_MPI_LIB}")
	endif()
	# Splitting external project add into 3 steps:
	# 1. To download the project
	# 2. To compile with single and double precision - requiores by GNU compilers
	ExternalProject_Add(${project}
		URL "${${project}_url}"
		URL_MD5 "${${project}_md5}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		# Ommiting Configuration, build and install steps
		CONFIGURE_COMMAND ""
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		# need to copy fftw to anotehr dir to run it in paralle
		COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_float" 
		)
	# FLOAT
	ExternalProject_Add(${project}_float
		DEPENDS ${project}	
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_float"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_float"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${${project}_float_configure_command}"
		)
	# DOUBLE
	ExternalProject_Add(${project}_double
		# specifying that this project depends on the previous one
		DEPENDS ${project}
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${${project}_double_configure_command}"
		)

	# adding fftw3, fftw3_threads, fftw3_mpi and fftws3_omp into a library variable
	# Defining this variable just to not to overwrite FFTW_LIBRARIES created by FindFFTW
	if(NOT FFTW_DOUBLE_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_DOUBLE_LIB}")
	endif()
	if(NOT FFTW_FLOAT_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_LIB}")
	endif()
	if(NOT FFTW_FLOAT_THREADS_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_OPENMP_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_OPENMP_LIB}")
	endif()
	if(NOT FFTW_FLOAT_MPI_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_mpi${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_MPI_LIB}")
	endif()
	#set(FFTW3_LIBRARIES 
	#	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX} 
	#	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f${CMAKE_STATIC_LIBRARY_SUFFIX} 
	#	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}
	#	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}
	#	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_mpi${CMAKE_STATIC_LIBRARY_SUFFIX}
	#	)
else()
	# adding empty targets in case FFTW was found on the system
	add_custom_target(${project} ALL "")
	add_custom_target(fftw_double ALL "")
	add_custom_target(fftw_float ALL "")
	set(FFTW3_LIBRARIES
		${FFTW_DOUBLE_LIB}
		${FFTW_FLOAT_LIB}
		${FFTW_FLOAT_OPENMP_LIB}
		${FFTW_FLOAT_THREADS_LIB}
		#${FFTW_FLOAT_MPI_LIB}
		)
endif()
