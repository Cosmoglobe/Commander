
COBJS  := comm_camb_mod.o cambtest.o

all : cambtest


comm_camb_mod.o :
cambtest.o      : comm_camb_mod.o


cambtest : $(COBJS)
	$(MPF90) -o cambtest $(COBJS) $(HEALPIX_LINK)

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.cpp
	$(MPCXX) $(CXXCOMP) -c $< 

%.f90 : %.f90.in
	$(TEMPITA) < "$<" > "$@"

clean :
	@rm -f *.o *.mod *.MOD *.a *~ cambtest
