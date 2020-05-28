# Project: Doxygen 
# File determines the location of Doxygen on the system.
# Doxygen is a documentation file generator
# Author: Maksym Brilenkov

if(DOXYGEN_BUILD_DOCS)
	set(ENV{PATH}
		 #$ENV{PATH}:${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl
		 $ENV{PATH}:${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
		 )
	find_package(FLEX)
	if(NOT FLEX_FOUND)
		message(STATUS "Will download flex from source.")
		ExternalProject_Add(flex
			URL ${flex_url}
			PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/flex"
			DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
			BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/flex/src/flex"
			INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
			CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR>
			)
		set(FLEX_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/flex)
		message(STATUS "${FLEX_EXECUTABLE}")
	else()
		message(STATUS "Found Flex")
		add_custom_target(flex ALL "")
	endif()

	find_package(BISON)
	if(NOT BISON_FOUND)
		message(STATUS "Will download Bison from source")
		ExternalProject_Add(bison
			URL ${bison_url}
			PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/bison"
			DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
			BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/bison/src/bison"
			INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
			CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR>
			)
		set(BISON_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bison)	
	else()
		message(STATUS "Found Bison")
		add_custom_target(bison ALL "")
	endif()

	find_package(Doxygen)
	if(NOT DOXYGEN_FOUND)
		message(STATUS "Will download Doxygen from source.")
		ExternalProject_Add(${project}
			DEPENDS flex bison
			URL "${${project}_url}"
			#GIT_REPOSITORY "${${project}_url}"
			DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}" #"${download_dir}"
			BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
			INSTALL_DIR ${CMAKE_INSTALL_OUTPUT_DIRECTORY} #"${out_install_dir}"
			CMAKE_ARGS
			-DFLEX_EXECUTABLE:PATH=${FLEX_EXECUTABLE}
			-DBISON_EXECUTABLE:PATH=${BISON_EXECUTABLE}
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_BUILD_TYPE=Release
			)
		# refer to cmake docs about doxygen variables and options
		#message(STATUS "Building Commander3 documentation from source code with Doxygen.")
		#set(DOXYGEN_PROJECT_NAME "Commander3")
		#set(DOXYGEN_PROJECT_NUMBER "${CMAKE_PROJECT_VERSION}")
		#set(DOXYGEN_PROJECT_BRIEF "${CMAKE_PROJECT_DESCRIPTION}")
		#set(DOXYGEN_GENERATE_HTML YES)
		#set(DOXYGEN_GENERATE_LATEX YES)
		#set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build/docs)
		
		#doxygen_add_docs(${project}_docs
		#	ALL
		#	${COMMANDER3_SOURCE_DIR}
		#	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
		#	)
	else()
		add_custom_target(${project} ALL "")
		# refer to cmake docs about doxygen variables and options
		message(STATUS "Building Commander3 documentation from source code with Doxygen.")
		set(DOXYGEN_PROJECT_NAME "Commander3")
		set(DOXYGEN_PROJECT_NUMBER "${CMAKE_PROJECT_VERSION}")
		set(DOXYGEN_PROJECT_BRIEF "${CMAKE_PROJECT_DESCRIPTION}")
		set(DOXYGEN_GENERATE_HTML YES)
		set(DOXYGEN_GENERATE_LATEX YES)
		set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build/docs)
		
		doxygen_add_docs(${project}_docs
			ALL
			${COMMANDER3_SOURCE_DIR}
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
			)
	endif()

else()
	# making an empty target so the project will compile regardless of doxygen
	add_custom_target(${project} ALL "")
endif()

#message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
