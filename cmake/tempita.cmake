# Project: Tempita 
# File which contains setup for current project 
# Tempita used to deal with comm_hdf_mod.f90.in
# Author: Maksym Brilenkov

add_custom_target(${project} ALL "")
set(comm_hdf_mod "${COMMANDER3_SOURCE_DIR}/comm_hdf_mod.f90")
# running python command at configure time
execute_process(
	COMMAND ${TEMPITA_DIR}/tempita_proc.py < $< > $@
	INPUT_FILE ${comm_hdf_mod}.in
	OUTPUT_FILE ${comm_hdf_mod}
	)

#message(${comm_hdf_mod})
# the Makefile equivalent of this command is:
# OUTPUT: MAIN_DEPENDENCY DEPENDS
#        COMMAND
#add_custom_command(OUTPUT ${comm_hdf_mod} #"${CMAKE_SOURCE_DIR}/src/commander/comm_hdf_mod.f90"
#	MAIN_DEPENDENCY ${comm_hdf_mod}.in #"${CMAKE_SOURCE_DIR}/src/commander/comm_hdf_mod.f90.in"
#	DEPENDS 
#	COMMAND echo ${comm_hdf_mod}.in #${tempita_dir}/tempita_proc.py < $< > $@ 
#	COMMAND ${tempita_dir}/tempita_proc.py < $< > $@ 
#	)

# target zoo is always built
#add_custom_target(${project} ALL
#    COMMAND echo "This is ALL target '${project}', and it depends on ${comm_hdf_mod}"
#    # If the file exists, then commands related to that file won't be executed
    # DONOT let other target depends on the same OUTPUT as current target,
    #   or it may be bad when doing parallel make
	#    DEPENDS ${comm_hdf_mod}
    # to make quotes printable,for example
	#VERBATIM
	#)
