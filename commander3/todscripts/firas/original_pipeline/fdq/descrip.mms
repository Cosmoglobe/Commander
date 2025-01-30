! Builds facility FDQ.
!
! Author: Jeff Durachta
!          February 18, 1987
!	   ST Systems Corporation
!
!         Adapted from work done by Rob Kummerer and Wes Young
!		and Doug Ward for FDQ.WORK, October 1987
! Modified: 
!	    Shirley M. Read
!	    January 1988
!	    ST Systems Corporation
!
!	    Added new routines for Build 3.1, modified routines with 
!	    include files and standardized names to begin with FDQ.
!
!	    Shirley M. Read
!	    March 1988
!	    ST Systems Corporation
!
!	    Added new routines and deleted others for Build 3.2, modified 
!	    routines using the new COBETRIEVE FORTRAN user open and query
!	    catalog routines. 
!
!	    Shirley M. Read
!	    June 1988
!	    ST Systems Corporation
!
!	    Added new routines and replaced others for Build 4.0.
!
!	    Shirley M. Read
!	    October 19, 1988
!	    ST Systems Corporation
!
!	    Added dependency for FDQ_Set_Limflags.
!
!	    R. Kummerer
!	    February 12, 1990
!	    ST Systems Corporation
!
!	    SPR 6006, VMS 5.2 options file changes.
!
!	    R. Kummerer
!	    March 2, 1990
!	    ST Systems Corporation
!
!	    SPR 6353, Prevent MMS abort on warning "Multiply defined symbol
!		      WAIT".
!    Modified by: Harte Wang, STX
!            FDQ new requirement 1/29/91
!
!           S. Alexander, April 29, 1991: Provisionally remove H_CONV.
!           S. Alexander, April 27, 1992: Permanently remove H_CONV; now
!                                         moved to facility FHC. SPR 9642.
!
.suffixes
.suffixes .exe .olb .obj .for .cld .msg .tlb .txt

FDQ : FDQBLD.TLB, FDQ.EXE
  !Facility FDQ is up-to-date.

INC_FILES =     FDQ_NEW_IDX_PARS.TXT,-
		FDQ_OPTIONS.TXT,-
		FDQ_STATBUFF.TXT

TEXT_LIBS = FDQBLD.TLB/LIB + CSDR$LIBRARY:CSDRLIB.TLB/LIB + -
	    CSDR$LIBRARY:FUTLIB.TLB/LIB

! Define Macros:
!
! *** Text Libraries ( .TLB files) were not created because the source
! *** files had hardcoded references to complete file names, including
! *** .CMN, .INC and other extensions.
!
! Define rules:

.FOR.OBJ
  $(FORT) /EXTEND_SOURCE/CONTINUATIONS=99 $(FFLAGS) $(MMS$SOURCE) + $(TEXT_LIBS)

!Build the DATA_QUALIFY program.

FDQ_OBJ =  	FDQ, -
                FDQ_CLOSE_ARCV, -
                FDQ_CONVERT, -
                FDQ_FORM_STAT, -
                FDQ_GEN_STAT, -
                FDQ_GET_BRK_HSK, -
                FDQ_GET_ENG_FLDS, -
                FDQ_GET_HSK_FLDS, -
		FDQ_GET_TIME_RANGE, -
                FDQ_GET_NXT_SEGMENT, -
		FDQ_GET_NXT_SET, -
		FDQ_GET_OPTIONS, -
                FDQ_GET_SCI_FLDS, -
		FDQ_GET_SC_CONFIG, -
		FDQ_GET_VALID_SEGMENTS, -
                FDQ_INTERPOLATE, -
                FDQ_LOAD_CONV, -
                FDQ_MAKE_IDX, -
                FDQ_NEW_IDX, -
		FDQ_OL_GRT_CVT, -
		FDQ_OL_POLY_CVT, -
                FDQ_OPEN_ARCV, -
		FDQ_PROC_CUR_SET, -
		FDQ_SAVE_CAT_INFO, -
		FDQ_SET_LIMFLAGS, -
		FDQ_STAT_INIT, -
		FDQ_WITHIN, -
		FDQ_MSG

FDQBLD.TLB : $(INC_FILES)
  LIBRARY/CREATE/TEXT FDQBLD $(INC_FILES)
  ! FDQBLD.TLB has been built.

FDQ_MSG.OBJ : FDQ_MSG.MSG
  MESSAGE FDQ_MSG.MSG

FDQ.EXE :  FDQBLD.OLB($(FDQ_OBJ)), -             
	   CSDR$LIBRARY:FUTLIB.OLB, -
	   CSDR$LIBRARY:V5.OPT
 - $(LINK) $(LINKFLAGS) FDQBLD/LIB/INC=(FDQ), -
		        CSDR$LIBRARY:FUTLIB.OLB/LIB, -
		        CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		        CSDR$LIBRARY:CSDRMSG.OLB/INC=(CCT_MSG,CUT_MSG), -
			CSDR$LIBRARY:CSDRIMSL.OLB/LIB, -
			IMSL$DIR:IMSL.OLB/LIB, -
		        CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FDQ_DATA_QUALIFY has been built.
 !
!
! FOR FDQ
!
FDQ.OBJ              :  FDQBLD.TLB(FDQ_OPTIONS), -
		        CT$LIBRARY:CTUSER.INC, -
		        CSDR$LIBRARY:FUTLIB.TLB
FDQ_CLOSE_ARCV.OBJ   :  CT$LIBRARY:CTUSER.INC
FDQ_CONVERT.OBJ      :  CSDR$LIBRARY:FUTLIB.TLB
FDQ_FORM_STAT.OBJ    :  CT$LIBRARY:CTUSER.INC, -
                        FDQBLD.TLB(FDQ_STATBUFF)
FDQ_GEN_STAT.OBJ     :  CT$LIBRARY:CTUSER.INC, -
                        FDQBLD.TLB(FDQ_STATBUFF)
FDQ_GET_BRK_HSK.OBJ  :  CT$LIBRARY:CTUSER.INC
FDQ_GET_NXT_SEGMENT.OBJ  :  FDQBLD.TLB(FDQ_OPTIONS), -
			CT$LIBRARY:CTUSER.INC
FDQ_GET_NXT_SET.OBJ  :	CT$LIBRARY:CTUSER.INC
FDQ_GET_OPTIONS.OBJ  :	FDQBLD.TLB(FDQ_OPTIONS), -
			CSDR$LIBRARY:FUTLIB.TLB
FDQ_INTERPOLATE.OBJ  :  CT$LIBRARY:CTUSER.INC
FDQ_LOAD_CONV        :  CSDR$LIBRARY:FUTLIB.TLB
FDQ_MAKE_IDX.OBJ     :  CT$LIBRARY:CTUSER.INC
FDQ_NEW_IDX.OBJ      :  FDQBLD.TLB(FDQ_NEW_IDX_PARS), -
			CSDR$LIBRARY:FUTLIB.TLB
FDQ_OL_GRT_CVT.OBJ   :  CSDR$LIBRARY:FUTLIB.TLB
FDQ_OL_POLY_CVT.OBJ  :  CSDR$LIBRARY:FUTLIB.TLB
FDQ_OPEN_ARCV.OBJ    :  CT$LIBRARY:CTUSER.INC
FDQ_PROC_CUR_SET.OBJ :	CT$LIBRARY:CTUSER.INC, -
		        CSDR$LIBRARY:FUTLIB.TLB
FDQ_SAVE_CAT_INFO.OBJ :	
FDQ_SET_LIMFLAGS.OBJ :  CSDR$LIBRARY:FUTLIB.TLB
FDQ_STAT_INIT.OBJ    :  CT$LIBRARY:CTUSER.INC, -
			FDQBLD.TLB(FDQ_STATBUFF)
FDQ_GET_TIME_RANGE.OBJ    :  
