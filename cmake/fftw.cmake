# Project: FFTW
# File which contains setup for current project 
# Author: Maksym Brilenkov

# looking for FFTW3 library on the system
#find_package(FFTW3 3.3.8)
#if(NOT FFTW3_FOUND)
#	ExternalProject_Add(${project}
#		URL "${${project}_url}"
#		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
#		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
#		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
#		INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
#		# commands how to build the project
#		CONFIGURE_COMMAND "${${project}_configure_command}"
#		COMMAND ./configure --prefix=<INSTALL_DIR> --enable-float --enable-threads --enable-openmp --enable-mpi 
#		#BUILD_IN_SOURCE 1	
#		)
#	set(FFTW3_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX})
#else()
#	add_custom_target(${project} ALL "")
#	message(STATUS "FFTW3 LFLAGS are: ${FFTW3_LFLAGS}")
#	# this is an include dirs provided by custom FindFFTW3.cmake
#	include_directories(${FFTW3_INCLUDES})
#endif()

# To get rid of "C preprocessor fails sanity check" error
# we need to link the CPP as follows
#set(FFTW_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
# splitting external project add into 3 steps:
# 1. To download the project
# 2. To compile with single and double precision - requiores by GNU compilers
ExternalProject_Add(${project}
	URL "${${project}_url}"
	URL_MD5 "${${project}_md5}"
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
	DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
	BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
	INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
	# Ommiting Configuration, build and install steps
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
	# need to copy fftw to anotehr dir to run it in paralle
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_f" 
	#COMMAND ./configure --prefix=<INSTALL_DIR> #--enable-float --enable-threads --enable-openmp --enable-mpi 
	# setting the compiler environment variable (gives weird errors for intel compiler otherwise)
	#COMMAND ${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} MPCC=${CMAKE_C_COMPILER} MPFC=${CMAKE_Fortran_COMPILER} MPCXX=${CMAKE_CXX_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-float --enable-threads --enable-openmp --enable-mpi 

	#COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-float --enable-threads --enable-openmp --enable-mpi #--with-sgimp #--enable-mpi 
	#COMMAND ${CMAKE_COMMAND} --build . -j #make
	#COMMAND ${CMAKE_COMMAND} --target install #make install
	#COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-threads --enable-openmp --enable-mpi #--with-sgimp #--enable-mpi 
	#COMMAND ${CMAKE_COMMAND} --build . -j #make
	#COMMAND ${CMAKE_COMMAND} --target install #make install
	#BUILD_COMMAND ""
	#INSTALL_COMMAND ""
	)
ExternalProject_Add(${project}_f
	DEPENDS ${project}	
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
	SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_f"
	BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_f"
	INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
	# Disabling download
	DOWNLOAD_COMMAND ""
	BUILD_ALWAYS FALSE
	# Commands to configure, build and install the project
	#CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-float --enable-threads --enable-openmp #--enable-mpi
	# we need to copy source directory into another one, as configuration will create a Makefile
	# for some reason --enable-mpi options starts lookin g for gcc mpi and not ifort 
	# unless you export those variables in .bashrc. I have no idea why this happens,
	# so I decided to just let it be without mpi. We are not using it anyway.
	#CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CPP=${COMMANDER3_CPP_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-float --enable-threads --enable-openmp #--enable-mpi
	CONFIGURE_COMMAND "${${project}_f_configure_command}"
	)
ExternalProject_Add(${project}_
	# specifying that this project depends on the previous one
	DEPENDS ${project}
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
	SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
	BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
	INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
	# Disabling download
	DOWNLOAD_COMMAND ""
	# Commands to configure, build and install the project
	#CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-threads --enable-openmp #--enable-mpi
	#CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CPP=${COMMANDER3_CPP_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-threads --enable-openmp #--enable-mpi
	CONFIGURE_COMMAND "${${project}_configure_command}"
	BUILD_ALWAYS FALSE
	)


# adding fftw3, fftw3_threads, fftw3_mpi and fftws3_omp into a library variable
#set(FFTW3_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3_omp${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3_threads${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3_mpi${CMAKE_STATIC_LIBRARY_SUFFIX})
set(FFTW3_LIBRARIES 
	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX} 
	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f${CMAKE_STATIC_LIBRARY_SUFFIX} 
	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}
   	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}
		#${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_mpi${CMAKE_STATIC_LIBRARY_SUFFIX}
	)

#message(STATUS
#"
#FFTW VERSION: ${FFTW3_VERSION}
#FFTW INCLUDE DIRS: ${FFTW3_INCLUDE_DIRS}
#FFTW LIBRARIES: ${FFTW3_LIBRARIES}
#")
#message("${${project}_configure_command}")




# adding the library to the project
#add_library(${project}_lib STATIC IMPORTED)
#set(${${project}_lib}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(${${project}_lib} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${${${project}_lib}_name}")
#message("The ${${${project}_lib}_name} Path is " ${out_install_dir}/lib/${${${project}_lib}_name})
#include_directories(${out_install_dir}/include)

#add_library(${project} STATIC IMPORTED)
#set(lib_${project}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(${project} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${lib_${project}_name}") 
