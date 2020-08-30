# Project: HEALPix 
# File which contains setup for current project 
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")
#add_dependencies(${project} cfitsio sharp)
#target_link_libraries(${project} cfitsio)

#set(CMAKE_INSTALL_RPATH "${out_install_dir}")
#message(${CMAKE_INSTALL_RPATH})
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

ExternalProject_Add(${project}
	URL "${${project}_url}"
	URL_MD5 "${${project}_md5}"
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
	DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
	#SOURCE_DIR "${download_dir}/${project}/src/${project}"
	BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" 
	INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
	# commands how to build the project
	# setting mpi wrappers as compilers
	#CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} MPCC=${CMAKE_C_COMPILER} MPFC=${CMAKE_Fortran_COMPILER} MPCXX=${CMAKE_CXX_COMPILER} ./configure #COMMAND cd <SOURCE_DIR> 
	#CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CPP=${COMMANDER3_CPP_COMPILER} CC=${MPI_C_COMPILER} ./configure 
	
	CONFIGURE_COMMAND "${${project}_configure_command}"
	#COMMAND bash "echo 3"
	#COMMENT "3"	
	#COMMAND ${CMAKE_COMMAND} cd "${download_dir}/${project}/src/${project}" #"${${project}_configure_command}"
	#COMMAND ${download_dir}/${project}/src/${project}/configure --auto=f90 --prefix=<INSTALL_DIR>
	#COMMAND ./configure #--prefix=<INSTALL_DIR>
	# making healpix to be installed the last before commander3
	DEPENDS cfitsio 
					hdf5 
					#sharp2 
					fftw 
					fftw_double 
					fftw_float 
					doxygen 
					tempita 
					blas 
					openmp 
					curl 
					mpi 
					zlib
	#
	INSTALL_COMMAND ""
	# copying Healpix and all its files (src and compiled) into CMAKE_INSTALL_PREFIX directory
	#COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}/build" "${CMAKE_INSTALL_PREFIX}/healpix_build"
	#COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}/lib" "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"	
	#COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}/data" "${CMAKE_INSTALL_PREFIX}/healpix_build/data"
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${CMAKE_INSTALL_PREFIX}/healpix"
	#DEPENDS cfitsio # <= it still relies on cfitsio which is inside /lib/64 or something and not currently installed one
	)

#ExternalProject_Add_Step(${project} 3
#		COMMAND 3
#		DEPENDEES configure
#		)

#execute_process(COMMAND "${${project}_configure_command}")

#add_custom_command(OUTPUT var1
#	MAIN_DEPENDENCY ${download_dir}/${project}/src/${project}/hpxconfig_functions.sh
#	DEPENDS
#	COMMAND yes | ./hpxconfig_functions.sh #echo ${download_dir}/${project}/src/${project}
	#COMMAND cd "${download_dir}/${project}/src/${project}"
	#COMMAND ${download_dir}/${project}/src/${project}/configure
	#	)

	#add_custom_target(varxxx1 ALL
	#	COMMAND echo "Some random stuff"
	#DEPENDS var1
	#)

#execute_process(COMMAND cd ${download_dir}/${project}/src/${project})
#execute_process(COMMAND ./configure)
#add_custom_target(${project} ALL)

#set(HEALPIX_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(HEALPIX_LIBRARIES 
	${CMAKE_INSTALL_PREFIX}/healpix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sharp${CMAKE_STATIC_LIBRARY_SUFFIX}
	${CMAKE_INSTALL_PREFIX}/healpix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX}
	)
include_directories("${CMAKE_INSTALL_PREFIX}/healpix/include")
message(STATUS "HEALPIX LIBRARIES will be: ${HEALPIX_LIBRARIES}")

#set(CMAKE_Fortran_MODULE_DIRECTORY ${out_install_dir}/mod)
#add_library(${project}_lib STATIC IMPORTED)
#set(${${project}_lib}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(${${project}_lib} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${${${project}_lib}_name}")
#message("The ${${${project}_lib}_name} Path is " ${out_install_dir}/lib/${${${project}_lib}_name})
