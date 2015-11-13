
COBJS  := messcomm.o

all : messcomm

messcomm : messcomm.o
	$(MPF90) -o messcomm messcomm.o -L../commander $(LINK) $(MPFCLIBS)

%.o : %.F90
	$(MPF90) $(F90COMP) -I../commander -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -I../commander -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.cpp
	$(MPCXX) $(CXXCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ messcomm
