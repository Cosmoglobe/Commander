! Builds facility FGA_GET_ARCHIVE.
!
! Changes:
!	Conversion to facility FGA.  R. Kummerer, STX, November 14, 1987.
!	Addition of error handler FGA_ERROR.  F. Shuman, STX, 1989 Jan 4.
!	SPR 6011, VMS 5.2 options file changes.  R. Kummerer, STX, 1990 Feb 12.
!
!       SPR 8372, Add new and modefied rdl's. N. Gonzales/STX, 1991 Sept 10.  
!                 Subroutines: fga_list_facr, fga_list_faker, fga_list_fakel,
!                              fga_list_gainr, fga_list_gainl, fga_list_gtran,
!                              fga_list_grawt, fga_list_mtmsweep, 
!                              fga_list_vabsaa
!       SPR 9099, Add new subroutine fga_list_mincoadd, 
!                 N. Gonzales/HughesSTX, 1991 Dec 12.
!
!       SPR 9795, Remove obsolete datasets from FCI, FES, FCS, FFC, FPR, 
!                                               FPS and FSF.
!                 N. Gonzales/HughesSTX, 1992 July 1.
!
!	SPR 10487, Remove references to FGA_LIST_CTH and FGA_LIST_NOS.
!		   S. Alexander/Hughes STX, 1993, January 21.
!
!	SPR 11741, Remove references to FGA_LIST_GRTSWT.
!   		   S. Alexander/Hughes STX, 1994, May 2.

FGA : FGABLD.TLB, FGA.EXE
  !Facility FGA is up-to-date.

! Define Macros:

TEXT_LIBS = FGABLD.TLB/LIB + CSDR$LIBRARY:FUTLIB.TLB/LIB + -
	    CSDR$LIBRARY:CSDRLIB.TLB/LIB

! Define rules:

.FOR.OBJ
  $(FORT)/EXTEND_SOURCE/CONTINUATIONS=99 $(FFLAGS) $(MMS$SOURCE) + $(TEXT_LIBS)

inc_files =	fga_error.txt

!Build FGA.

FGA_OBJ = FGA.OBJ, -
	  FGA_LIST_ENG.OBJ, -
	  FGA_LIST_ENGANLG.OBJ, -
	  FGA_LIST_ENGSTAT.OBJ, -
	  FGA_LIST_HKP.OBJ, -
	  FGA_LIST_IDX.OBJ, -
	  FGA_LIST_SCI.OBJ, -
	  FGA_LIST_EMF.OBJ, -
	  FGA_LIST_SSC.OBJ, -
	  FGA_LIST_EXT.OBJ, -
	  FGA_LIST_ANC.OBJ, -
	  FGA_LIST_ATT.OBJ, -
	  FGA_LIST_FACR.OBJ, -
	  FGA_LIST_FAKER.OBJ, -
	  FGA_LIST_FAKEL.OBJ, -
	  FGA_LIST_GAINR.OBJ, -
	  FGA_LIST_GAINL.OBJ, -
	  FGA_LIST_GTRAN.OBJ, -
	  FGA_LIST_GRAWT.OBJ, -
	  FGA_LIST_MTMSWEEP.OBJ, -
	  FGA_LIST_VABSAA.OBJ, -
	  FGA_LIST_MINCOADD.OBJ, -
	  FGA_LIST_SCILIM.OBJ, -
	  FGA_LIST_ENGLIM.OBJ, -
	  FGA_LIST_IDXFLAGS.OBJ, -
	  FGA_LIST_IDXTOLS.OBJ, -
	  FGA_LIST_LIMFLAGS.OBJ, -
	  FGA_ERROR.OBJ, -
	  FGA_MSG.OBJ

FGA_MSG.OBJ : FGA_MSG.MSG
 MESSAGE FGA_MSG

! Build the libraries.

fgabld.tlb : $(inc_files)
  LIBRARY/CREATE/TEXT fgabld $(inc_files)
  ! FGABLD.TLB has been built.


FGA.EXE :  FGABLD.OLB($(FGA_OBJ)), -
	   CSDR$LIBRARY:FUTLIB.OLB, -
	   CSDR$LIBRARY:V5.OPT
 LINK $(linkflags) FGABLD.OLB/LIB/INC=(FGA), -
			CSDR$LIBRARY:FUTLIB.OLB/LIB,-
			CSDR$LIBRARY:CSDRLIB/LIB,-
			CSDR$LIBRARY:V5.OPT/OPTION
 !FGA has been built.
 !

!
! Further compilation dependencies
!
FGA.OBJ		 : FDQ_ENG^, NFS_HKP^, -
		   NFS_SDF^, NFS_EMF^, FDQ_IDX^, -
		   FNT_NOISE^
FGA_LIST_ENG.OBJ : FDQ_ENG^
FGA_LIST_ENGSTAT.OBJ : FUT_ENGSTAT^
FGA_LIST_ENGANLG.OBJ : FUT_ENGANLG^
FGA_LIST_HKP.OBJ : NFS_HKP^
FGA_LIST_SCI.OBJ : NFS_SDF^
FGA_LIST_EMF.OBJ : NFS_EMF^
FGA_LIST_IDX.OBJ : FDQ_IDX^
FGA_LIST_SSC.OBJ : FEC_SSCAL^
FGA_LIST_EXT.OBJ : FXT_ENG_XTRM^
FGA_LIST_ANC.OBJ : NFS_ANC^
FGA_LIST_ATT.OBJ : FUT_ATTIT^
FGA_LIST_FACR.OBJ : FEX_AV_CALRS^
FGA_LIST_FAKER.OBJ : FEX_FAKEIT^
FGA_LIST_FAKEL.OBJ : FEX_FAKEIT^
FGA_LIST_GAINR.OBJ : FEX_GAIN^
FGA_LIST_GAINL.OBJ : FEX_GAIN^
FGA_LIST_GTRAN.OBJ : FEX_GRTTRANS^
FGA_LIST_GRAWT.OBJ : FEX_GRTRAWWT^
FGA_LIST_MTMSWEEP.OBJ : FEX_MTMSWEEP^
FGA_LIST_VABSAA.OBJ : FEX_VABSAA^
FGA_LIST_MINCOADD.OBJ : FEX_MINCOADD^
FGA_LIST_SCILIM.OBJ : FEX_SCILIM^
FGA_LIST_ENGLIM.OBJ : FEX_ENGLIM^
FGA_LIST_IDXFLAGS.OBJ : FEX_IDX_FLAG^
FGA_LIST_IDX_TOLS.OBJ : FEX_IDX_TOLS^
FGA_LIST_LIMFLAGS.OBJ : FEX_LIMFLAGS^
