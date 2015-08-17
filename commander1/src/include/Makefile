
HOBJS  := d1mach.o drc3jj.o \
          math_tools.o \
          mpi_alm_tools.o sort_utils.o stat_gof_mod.o ARS_mod.o \
	  locate_mod.o spline_1D_mod.o spline_2D_mod.o powell_mod.o InvSamp_mod.o \
	  quasi_newton_mod.o


libcommander_healpix.a : $(HOBJS)
	$(AR) $(ARFLAGS) libcommander_healpix.a $(HOBJS)
	$(RANLIB) libcommander_healpix.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~

