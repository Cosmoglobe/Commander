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
ExternalProject_Add(${project}
	URL "${${project}_url}"
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
	DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
	BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
	INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
	# commands how to build the project
	CONFIGURE_COMMAND "${${project}_configure_command}"
	#COMMAND ./configure --prefix=<INSTALL_DIR> #--enable-float --enable-threads --enable-openmp --enable-mpi 
	# setting the compiler environment variable (gives weird errors for intel compiler otherwise)
	COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} MPCC=${CMAKE_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-float --enable-threads --enable-openmp --with-sgimp #--enable-mpi 
	COMMAND make
	COMMAND make install
	COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} MPCC=${CMAKE_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-threads --enable-openmp --with-sgimp #--enable-mpi 
	COMMAND make
	COMMAND make install
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
	#BUILD_IN_SOURCE 1	
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
