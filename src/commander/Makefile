
COBJS  := sharp.o sort_utils.o comm_system_backend.o comm_system_mod.o comm_map_mod.o \
          comm_defs.o comm_utils.o \
          comm_shared_output_mod.o comm_status_mod.o \
	  comm_bp_utils.o comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o \
	  comm_ptsrc_comp_mod.o comm_cmb_comp_mod.o comm_powlaw_comp_mod.o \
	  comm_signal_mod.o comm_N_mod.o comm_N_rms_mod.o \
          comm_B_mod.o comm_B_bl_mod.o comm_beam_mod.o \
          comm_data_mod.o comm_bp_mod.o

all : commander

comm_map_mod.o :           sharp.o
comm_system_mod.o :        comm_system_backend.o
comm_shared_output_mod.o : comm_system_mod.o
comm_status_mod.o :        comm_shared_output_mod.o
comm_utils.o :             sort_utils.o comm_system_mod.o sharp.o
comm_bp_utils.o :          comm_utils.o
comm_param_mod.o :         comm_utils.o
comm_comp_mod.o :          comm_bp_mod.o comm_param_mod.o
comm_diffuse_comp_mod.o :  comm_comp_mod.o
comm_ptsrc_comp_mod.o :    comm_comp_mod.o
comm_cmb_comp_mod.o :      comm_diffuse_comp_mod.o
comm_powlaw_comp_mod.o :   comm_diffuse_comp_mod.o	
comm_signal_mod.o :        comm_powlaw_comp_mod.o comm_cmb_comp_mod.o \
                           comm_ptsrc_comp_mod.o
comm_N_mod.o :             comm_param_mod.o
comm_N_rms_mod.o :         comm_N_mod.o
comm_noise_mod.o :         comm_N_rms_mod.o
comm_bp_mod.o :            comm_bp_utils.o comm_param_mod.o
comm_fg_component_mod.o :  comm_bp_mod.o
comm_data_mod.o :          comm_map_mod.o comm_noise_mod.o comm_bp_mod.o comm_beam_mod.o
commander.o :              comm_signal_mod.o comm_data_mod.o


commander : libcommander.a commander.o
	$(MPF90) -o commander commander.o $(LINK) $(MPFCLIBS)

libcommander.a : $(COBJS)
	$(AR) $(ARFLAGS) libcommander.a $(COBJS) 
	$(RANLIB) libcommander.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.cpp
	$(MPCXX) $(CXXCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ commander
