# Project: CFitsio 
# File which contains setup for current project 
# Author: Maksym Brilenkov


message(STATUS "---------------------------------------------------------------")
#find_package(FITSIO)
#set($ENV{PATH} ${out_install_dir}/bin/)
#execute_process(COMMAND export PATH=$PATH:/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild2/build/install/bin) #"${${project}_configure_command}"
#add_custom_command(
#	OUTPUT ""
#	COMMAND export "PATH=$PATH:/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild2/build/install/bin"
#	) #"${${project}_configure_command}"
# seeting an environment variable for cfitsio to find curl
#add_custom_target(environment_command
#	COMMAND ${CMAKE_COMMAND} -E env "PATH=$ENV{PATH}" COMMAND export PATH=/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild2/build/install/bin #$ENV{PATH} #<real_command> args...
#)
ExternalProject_Add(${project}
	# specifying that cfitsio depends on the curl project and should be built after it
	DEPENDS curl #curl_lib curl_bin #"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/downloads/curl/src/curl/lib/curl_config.h"#curl_lib
	URL "${${project}_url}"
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
	DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}" #"${download_dir}"
	BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
	INSTALL_DIR "${CMAKE_INSTALL_PREFIX}" #"${out_install_dir}"
	# commands how to build the project
	CONFIGURE_COMMAND "${${project}_configure_command}"
	#COMMAND ${CMAKE_COMMAND} -E env --unset=PATH PATH=$ENV{PATH} ./configure --prefix=<INSTALL_DIR> 
	#COMMAND ./configure --prefix=<INSTALL_DIR> #"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild2/build/install/bin"
	#COMMAND export PATH=${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl-config #"${${project}_configure_command}"
	#COMMAND export PATH=$PATH:/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild2/build/install/bin #"${${project}_configure_command}"
	#COMMAND ${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} MPCC=${CMAKE_C_COMPILER} MPFC=${CMAKE_Fortran_COMPILER} MPCXX=${CMAKE_CXX_COMPILER} ./configure --prefix=<INSTALL_DIR> --disable-curl # <= if specified manually, the cmake args will not work
	#COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CPP=${COMMANDER3_CPP_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --disable-curl # <= if specified manually, the cmake args will not work
	#CMAKE_ARGS
	# specifying where to find curl library
	#-DCURL_INCLUDE_DIR:PATH=${CURL_INCLUDE_DIR}
	#-DCURL_BINARY_DIR:PATH=${CURL_BINARY_DIR}
	#-DCURL_SOURCE_DIR:PATH=${CURL_SOURCE_DIR}
	#-DCURL_LIBRARIES:PATH=${CURL_LIBRARIES}
	# specifying where to install CFitsio
	# (and how to install it)
	#-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
	#-DBUILD_SHARED_LIBS:BOOL=OFF
	# specifying build command
	#BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config Debug #--target INSTALL
	#BUILD_IN_SOURCE 1	
	)


set(CFITSIO_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
message(STATUS "CFITSIO LIBRARIES will be: ${CFITSIO_LIBRARIES}")
