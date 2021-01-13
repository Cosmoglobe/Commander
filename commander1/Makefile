# Commander Makefile for use with EITHER the planck
# module build system OR a set of stand-alone
# config files.  Default make target prints help.
#
# This makefile requires the use of GNU make or
# compatible tools that allow the substitution
# operator ":=".
#

NOW := $(shell date +%Y%m%d-%H%M)
DIR := commander_$(NOW)

TOPDIR := $(shell pwd)
export TOPDIR


ifdef COMMANDER
	include $(TOPDIR)/config/config.$(COMMANDER)
	ifndef INSTALL
		INSTALL := $(TOPDIR)/install_$(COMMANDER)
	endif
else
	$(error COMMANDER undefined)UNDEFINED
endif

ifndef MAKE
	export MAKE := make
endif
ifndef AR
	export AR := ar
endif
ifndef ARFLAGS
	export ARFLAGS := crv
endif
ifndef RANLIB
	export RANLIB := ranlib
endif
ifndef CTAR
	export CTAR := tar czvf
endif
ifndef F90
	export F90 := f90
endif
ifndef MPF90
	export MPF90 := mpif90
endif
ifndef F90FLAGS
	export F90FLAGS := -g -O2
endif
ifndef MPF77
	export MPF77 := mpif77
endif
ifndef FFLAGS
	export FFLAGS := -g -O2
endif
ifndef MPCC
	export MPCC := cc
endif
ifndef CFLAGS
	export CFLAGS := -g -O2
endif
ifndef MPCXX
	export MPCXX := c++
endif
ifndef CXXFLAGS
	export CXXFLAGS := -g -O2
endif
ifndef LDFLAGS
	export LDFLAGS := -lm
endif
ifndef FORTRAN_UPPER
	export FORTRAN_UPPER := 0
endif
ifndef CFITSIO_LINK
	export CFITSIO_LINK := -L/usr/local/lib -lcfitsio
endif
ifndef LAPACK_LINK
	export LAPACK_LINK := -L/usr/local/lib -llapack -lblas
endif


export F90COMP := $(F90FLAGS) -I../include $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export FCOMP := $(FFLAGS)  $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export CCOMP := $(CFLAGS)  $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export LINK := -L../commander -L../include -lcommander -lcommander_healpix $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK) $(LDFLAGS) $(F90OMPFLAGS)


all : commander comm_process_resfiles

help :
	@echo ' '
	@echo '  This Makefile is used to build Commander in a way that is'
	@echo '  customized to your system.  You must export the COMMANDER'
	@echo '  environment variable and set it to the name of your platform.'
	@echo '  Then you must create a config file named config/config.<$COMMANDER>.'
	@echo '  I suggest copying the example config file and modifying the'
	@echo '  parameters for your platform.'
	@echo ' '
	@echo '  The following make targets are supported:'
	@echo ' '
	@echo '    make         : build everything'
	@echo '    make help    : print this help screen'
	@echo '    make install : install everything'
	@echo '    make clean   : remove build files'
	@echo '    make dist    : construct a date-stamped source tarball'
	@echo ' '

install : all
	@echo 'installing into:' $(INSTALL)
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	@mkdir -p $(INSTALL)/bin
	@cp src/include/libcommander_healpix.a $(INSTALL)/lib
	@cp src/commander/libcommander.a $(INSTALL)/lib
	@cp src/commander/commander $(INSTALL)/bin
	@cp src/comm_process_resfiles/comm_process_resfiles $(INSTALL)/bin
	@cp src/comm_process_resfiles/comm_like_tools $(INSTALL)/bin
	@cp src/comm_process_resfiles/comm_like_sampler $(INSTALL)/bin
	@if test $(FORTRAN_UPPER) = 1; then \
	cp src/commander/*.MOD $(INSTALL)/include; \
	cp src/include/*.MOD $(INSTALL)/include; \
	else \
	cp src/commander/*.mod $(INSTALL)/include; \
	cp src/include/*.mod $(INSTALL)/include; \
	fi

commander : libcommander_healpix
	@cd src/commander; $(MAKE)

libcommander_healpix :
	@cd src/include; $(MAKE)

comm_process_resfiles :
	@cd src/comm_process_resfiles; $(MAKE)

clean :
	@cd src/commander; $(MAKE) clean
	@cd src/comm_process_resfiles; $(MAKE) clean
	@cd src/include; $(MAKE) clean

dist : clean
	@mkdir $(DIR)
	@mkdir -p $(DIR)/src/commander
	@mkdir -p $(DIR)/src/comm_process_resfiles
	@mkdir -p $(DIR)/src/include
	@cp -r commander_test.tar.gz config Makefile $(DIR)
	@cp src/commander/*.f90 src/commander/Makefile $(DIR)/src/commander
	@cp src/comm_process_resfiles/*.f90 src/comm_process_resfiles/Makefile $(DIR)/src/comm_process_resfiles
	@cp src/include/*.f90 src/include/*.f src/include/*.c src/include/Makefile $(DIR)/src/include
	@rm -rf $(DIR)/config/.svn
	@$(CTAR) $(DIR).tar.gz $(DIR)
	@rm -rf $(DIR)

