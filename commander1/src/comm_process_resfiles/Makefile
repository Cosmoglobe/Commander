
#WMAP_DIR = /mn/stornext/u2/hke/owl/local/src/wmap_likelihood_v5
#WMAP = -DWMAP -L$(WMAP_DIR) -lwmap9 -I$(WMAP_DIR)
WMAP = 

COBJS := comm_like_utils.o comm_proc_utils.o comm_lowl_mod.o comm_br_old_mod.o comm_br_mod.o comm_gauss_br_mod.o comm_like_model_mod.o comm_like_optimization_mod.o comm_like_timing_mod.o

all: comm_process_resfiles comm_like_tools comm_like_sampler comm_like_tools_mpi

comm_like_tools.o : comm_gauss_br_mod.o 

comm_process_resfiles : $(COBJS) comm_process_resfiles.o
	$(MPF90) -qopenmp -o comm_process_resfiles comm_process_resfiles.o $(COBJS) $(LINK) 

comm_like_tools : $(COBJS) comm_like_tools.o
	$(MPF90) -o comm_like_tools comm_like_tools.o $(COBJS) $(LINK) $(WMAP)

comm_like_tools_mpi : $(COBJS) comm_like_tools_mpi.o
	$(MPF90) -o comm_like_tools_mpi comm_like_tools_mpi.o $(COBJS) $(LINK) $(WMAP)

comm_like_sampler : $(COBJS) comm_like_sampler.o
	$(MPF90) -o comm_like_sampler comm_like_sampler.o $(COBJS) $(LINK)

%.o : %.F90
	$(MPF90) $(F90COMP) -qopenmp -c $< $(WMAP) 

%.o : %.f90
	$(MPF90) $(F90COMP) -qopenmp -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ comm_process_resfiles comm_like_tools comm_like_sampler

