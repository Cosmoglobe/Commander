
COBJS  := comm_system_backend.o comm_system_mod.o comm_utils.o \
	  comm_bp_mod.o comm_task_mod.o comm_Cl_util_mod.o comm_genvec_mod.o \
          comm_S_mult_mod.o comm_beam_mod.o comm_N_mult_mod.o comm_S_mult_mod.o \
          comm_fg_component_mod.o comm_data_mod.o comm_chisq_mod.o comm_fg_mod.o  \
          comm_cgd_matmul_mod.o \
          comm_cgd_precond_mod.o comm_mp_mod.o \
          comm_cgd_mod.o comm_direct_sol_mod.o comm_Cl_sampling_mod.o comm_signal_mod.o \
          comm_mcmc_mod.o comm_noisesamp_mod.o comm_global_par_mod.o #comm_cl_spline_mod2.o

all : commander

comm_system_mod.o : comm_system_backend.o
comm_utils.o : comm_system_mod.o ../include/libcommander_healpix.a
comm_task_mod.o : comm_utils.o
comm_bp_mod.o : comm_utils.o
comm_Cl_util_mod.o : comm_utils.o
comm_genvec_mod.o : comm_utils.o
comm_S_mult_mod.o : comm_utils.o
comm_beam_mod.o : comm_utils.o
comm_N_mult_mod.o : comm_utils.o
comm_S_mult_mod.o : comm_utils.o
comm_fg_component_mod.o :: comm_bp_mod.o
comm_data_mod.o : comm_fg_component_mod.o comm_genvec_mod.o comm_N_mult_mod.o comm_beam_mod.o
comm_chisq_mod.o : comm_data_mod.o
comm_fg_mod.o    : comm_data_mod.o comm_task_mod.o
comm_cgd_matmul_mod.o : comm_data_mod.o comm_S_mult_mod.o
comm_cgd_precond_mod.o :comm_cgd_matmul_mod.o 
comm_mp_mod.o : comm_fg_mod.o comm_chisq_mod.o comm_cgd_precond_mod.o comm_cgd_matmul_mod.o
comm_cgd_mod.o : comm_mp_mod.o
comm_signal_mod.o : comm_cgd_mod.o
comm_direct_sol_mod.o : comm_mp_mod.o
comm_Cl_sampling_mod.o : comm_mp_mod.o 
comm_mcmc_mod.o : comm_mp_mod.o
comm_noisesamp_mod.o : comm_chisq_mod.o
comm_global_par_mod.o : comm_mp_mod.o
commander.o : comm_Cl_sampling_mod.o comm_signal_mod.o comm_mcmc_mod.o comm_noisesamp_mod.o 


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

