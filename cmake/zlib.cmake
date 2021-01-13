# Project: ZLIB 
# Project link: https://zlib.net/ 
# Description: A Massively Spiffy Yet Delicately Unobtrusive Compression Library
# Required by: HDF5 
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")
find_package(ZLIB REQUIRED)

message(STATUS "ZLIB LIBRARIES ARE: ${ZLIB_LIBRARIES}")
message(STATUS "ZLIB INCLUDE DIRS ARE: ${ZLIB_INCLUDE_DIRS}")

add_custom_target(${project} ALL "")
include_directories(${ZLIB_INCLUDE_DIRS})
add_library(zlib_lib SHARED IMPORTED GLOBAL) #${ZLIB_LIBRARIES})#ZLIB::ZLIB)
set_target_properties(zlib_lib PROPERTIES IMPORTED_LOCATION ${ZLIB_LIBRARIES})


# looking for curl in the system and download it if it is not present
#find_package(CURL)
#if(NOT CURL_FOUND)
#ExternalProject_Add(${project}
#		URL "${${project}_url}"
#		PREFIX "${download_dir}/${project}"
#		DOWNLOAD_DIR "${download_dir}"
#		BINARY_DIR "${download_dir}/${project}/src/${project}"
#		INSTALL_DIR "${out_install_dir}"
		#PATCH_COMMAND ./buildconf
		#		CONFIGURE_COMMAND "${${${project}_configure}_command}"
		#COMMAND ./configure --prefix=<INSTALL_DIR>
		#LOG_DOWNLOAD ON
		#LOG_UPDATE ON
		#LOG_CONFIGURE ON
		#LOG_BUILD ON
		#LOG_TEST ON
		#LOG_INSTALL ON
		#)
	# getting curl directories
	#ExternalProject_Get_Property(${project} source_dir)
	#ExternalProject_Get_Property(${project} install_dir)
	# specifying curl libraries and binaries
	#set(CURL_SOURCE_DIR ${source_dir})
	#set(CURL_BINARY_DIR ${install_dir}/bin)
	#set(CURL_INCLUDE_DIR ${install_dir}/include)#/${project})
	#set(CURL_LIBRARIES ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
	# including curl into a project
	#include_directories(${CURL_SOURCE_DIR})
	#include_directories(${CURL_BINARY_DIR})
	#include_directories(${CURL_INCLUDE_DIR})
	#message(STATUS "Curl INCLUDE DIR is ${CURL_INCLUDE_DIR}")
	#message(STATUS "Curl BINARY DIR is ${CURL_BINARY_DIR}")
	#message(STATUS "Curl SOURCE DIR is ${CURL_SOURCE_DIR}")
	#message(STATUS "Curl LIBRARIES are ${CURL_LIBRARIES}")

	#LIST(APPEND CMAKE_PROGRAM_PATH  "${CURL_BINARY_DIR}")
	#set(${project}_bin ${install_dir}/bin)
	# adding curl as an external library
	#add_library(${project}_lib STATIC IMPORTED GLOBAL)
	# asking cmake to identify its name for us
	#set(${${project}_lib}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
	# Specifying the directory
	#set_target_properties(${${project}_lib} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${${${project}_lib}_name}")
	#add_custom_target(${project}_bin WORKING_DIRECTORY ${out_install_dir}/bin)
	# exporting curl so, cfitsio will be able to find it
	#execute_process(COMMAND export PATH=$PATH:/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib)
	# setting an environment variable for cfitsio to find curl library
	#set($ENV{PATH} "${out_install_dir}/bin/curl")
	#else()
	#set(${project} "")
	#add_custom_target(curl "")
	#include_directories(${CURL_INCLUDE_DIR})
	#include_directories(${CURL_BINARY_DIR})
	#message(STATUS "Curl INCLUDE DIR is ${CURL_INCLUDE_DIR}")
	#message(STATUS "Curl BINARY DIR is ${CURL_BINARY_DIR}")
	#message(STATUS "Curl SOURCE DIR is ${CURL_SOURCE_DIR}")
	#message(STATUS "Curl LIBRARIES are ${CURL_LIBRARIES}")
	#endif()

