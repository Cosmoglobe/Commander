# Project: HDF5 
# File which contains setup for current project 
# Author: Maksym Brilenkov

# asking for an exact hdf5 version
#find_package(HDF5 1.10.5 EXACT COMPONENTS Fortran) #Fortran_HL)
find_package(HDF5 1.10.0 COMPONENTS Fortran) #Fortran_HL)
if(NOT HDF5_FOUND)
	message(STATUS "Will download HDF5v1.10.5 from source.")
	ExternalProject_Add(${project}
		DEPENDS zlib
		URL "${${project}_url}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
		# commands how to build the project
		CONFIGURE_COMMAND "${${project}_configure_command}"
		COMMAND ./configure --prefix=<INSTALL_DIR> --enable-fortran #--enable-cxx #--enable-parallel
		#BUILD_IN_SOURCE 1	
		)
	# adding hdf5_fortran and hdf5 into a variable (to be consistent with cmake docs)
	set(HDF5_Fortran_LIBRARIES 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhdf5_fortran.a" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhdf5.a" 
		"${ZLIB_LIBRARIES}")
	set(HDF5_Fortran_INCLUDE_DIRS "${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include/")
	include_directories(${HDF5_Fortran_INCLUDE_DIRS})
	message(STATUS "HDF5 Fortran LIBRARIES will be: ${HDF5_Fortran_LIBRARIES}")
	message(STATUS "HDF5 Fortran INCLUDE DIRS will be: ${HDF5_Fortran_INCLUDE_DIRS}")
	#set($ENV{PATH} "${out_lib_dir}/")
else()
	add_custom_target(hdf5 ""
		ALL
		DEPENDS zlib
		)
	message(STATUS "FOUND HDF5 Fortran!")
	include_directories(${HDF5_Fortran_INCLUDE_DIRS})
	message(STATUS "HDF5 Fortran INCLUDE DIRS are: ${HDF5_Fortran_INCLUDE_DIRS}")
	message(STATUS "HDF5 Fortran DEFINITIONS are: ${HDF5_Fortran_DEFINITIONS}")
	message(STATUS "HDF5 Fortran LIBRARIES are: ${HDF5_Fortran_LIBRARIES}")
endif()

# foor some reason it doesn't work by simply linking to ${HDF5_Fortran_LIBRARIES}
#add_library(hdf5_fortran STATIC IMPORTED GLOBAL)
#set_target_properties(hdf5_fortran PROPERTIES IMPORTED_LOCATION ${HDF5_Fortran_LIBRARIES})

#if(HDF5_FOUND)
#	message(STATUS
#	"
#	HDF5 VERSION: ${HDF5_VERSION}
#	HDF5 INCLUDE DIRS: ${HDF5_INCLUDE_DIR}
#	HDF5 LIBRARIES: ${HDF5_LIBRARIES}
#	"
	#	)
	#	include_directories(${HDF5_INCLUDE_DIRS})
	#else()
	#ExternalProject_Add(${project}
	#		URL "${${project}_url}"
	#	PREFIX "${download_dir}/${project}"
	#	DOWNLOAD_DIR "${download_dir}"
	#	BINARY_DIR "${download_dir}/${project}/src/${project}"
	#	INSTALL_DIR "${out_install_dir}"
		# commands how to build the project
		#	CONFIGURE_COMMAND "${${project}_configure_command}"
		#COMMAND ./configure --prefix=<INSTALL_DIR> --enable-fortran #--enable-parallel
		#)#BUILD_IN_SOURCE 1	
		#endif()

# define this variable here for easier reference in the future


#set(CMAKE_INSTALL_RPATH "${out_install_dir}")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#ExternalProject_Add(${project}
#	DEPENDS zlib
#	URL "${${project}_url}"
#	PREFIX "${download_dir}/${project}"
#	DOWNLOAD_DIR "${download_dir}"
#	BINARY_DIR "${download_dir}/${project}/src/${project}"
#	INSTALL_DIR "${out_install_dir}"
	# commands how to build the project
	#	CONFIGURE_COMMAND "${${project}_configure_command}"
	#	COMMAND ./configure --prefix=<INSTALL_DIR> --enable-fortran #--enable-cxx #--enable-parallel
	#BUILD_IN_SOURCE 1	
	#	)

# adding this library to the project
#add_library(${project}_lib STATIC IMPORTED GLOBAL)
#set(${${project}_lib}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(${${project}_lib} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${${${project}_lib}_name}")
#message("The ${${${project}_lib}_name} Path is " ${out_install_dir}/lib/${${${project}_lib}_name})

#add_library(hdf5_lib STATIC IMPORTED GLOBAL)
#set(hdf5_lib_name ${CMAKE_STATIC_LIBRARY_PREFIX}hdf5${CMAKE_STATIC_LIBRARY_SUFFIX})
#set(hdf5_fortran_lib_name ${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_fortran${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(hdf5_lib 
#	PROPERTIES IMPORTED_LOCATION 
#	"${out_install_dir}/lib/${hdf5_lib_name}" 
	#"${out_install_dir}/lib/${hdf5_fortran_lib_name}"
	#	)
	#message("The hdf5_lib_name Path is " ${out_install_dir}/lib/hdf5_lib_name)
