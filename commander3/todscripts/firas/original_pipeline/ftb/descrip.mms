! Builds facility FTB.
!
! Author: R. Kummerer
!         STX
!	  November 25, 1986
!
!         Adapted from work done by W. Young.
!
! Changes:
!	Include FTB_MAKE_RSE. Object library name should be FTBBLD, not
!		FTBLIB. R. Kummerer, February 9, 1987.
!	Include FTB_SNAPSHOT. R. Kummerer, February 18, 1987.
!	OPT file for SNAPSHOT,SEE,CTVIEWER. R. Kummerer, March 13, 1987.
!	Remove GETRSE references. R. Kummerer, March 24, 1987.
!	Include FTB_TESTLOG. R. Kummerer, April 1, 1987.
!       Include JDATE.  L. Olson, May 21, 1987
!       Add files created from the one module per file rule.
!		Reid Wilson, November 22, 1987
!       Include FTB_HISTO. R. Kummerer, Sept 1, 1988
!	Link all programs with FTB_Msg. Shirley M. Read, September 2, 1988.
!       Link Programs which use PLT to Xanadu Libs and remove FUTILIB.
!		Shirley M. Read, September 8, 1988 .
!       Added FTB_ACCUM_HISTO to the FTB_Histo Program. Deleted FTB_READ_HISTO.
!		Rick Shafer, October 7, 1988.
!       Version 4.1.1 10/25/88, SPR 2537, Shirley M. Read, STX
!	        All FIRAS facilities must be able to run in batch mode.
!	        FTB_Ctviewer needed this capability. Thus FTB_Get_Rse was added.
!       Include FTB_TEMPCONTROL dependencies. R. Kummerer, May 18, 1989.
!       Repartition FTB_WORD31/FTB_TEMPCONTROL dependencies. F. Shuman,
!               May 23, 1989.
!	Remove FTB_MAKE_RSE. SPR 3161. R.I.P.  R. Kummerer, July 18, 1989.
!	Remove FTB_SNAPSHOT. SER 3493. R.I.P.  R. Kummerer, August 18, 1989.
!	Add FTB_GAPFLAG. SER 4023. R. Kummerer, September 6, 1989.
!	Add FTB_OVERLAP. SER 4210. R. Kummerer, October 13, 1989.
!	VMS 5.2 options file changes. SPR 6021. R. Kummerer, February 12, 1990.
!	Prevent MMS abort on warning "Multiply defined symbol WAIT".  SPR 6357,
!		R. Kummerer, March 2, 1990. 
!	Add FTB_BINFLD. SER 6852. N. GONZALES, June 1, 1990.
!	Add FTB_CHECKSUM. SER 6852. N. GONZALES, June 1, 1990.
!	Add FTB_VALID_TIME. SER 6852. N. GONZALES, June 1, 1990.
!       Add FTB_CHECK_DQ. SER 6852. N. GONZALES, June 14, 1990
!       Remove FTB_GET_RSE. SPR 9430. N. GONZALES, January 23, 1992.
!       Remove FTB_GAPFLAG. SPR 9726. S. Alexander, June 29, 1992.
!       Remove obsolete datasets from FTB_CTVIEWER and FTB_Gettim; 
!              Ref. SPR 9798. N. Gonzales/Hughes STX, June 29, 1992.
!	Add FTB_TRANS. SER 10764, 10772. S. Alexander, April 5, 1993.
!	Add FTB_LTRANS. SPR 12285. S. Brodd, December 20, 1995.
!
FTB : FTB_CONVTIME.EXE, FTB_CTVIEWER.EXE, FTB_WORD31.EXE, FTB_TEMPCONTROL.EXE, -
      FTB_GETTIM.EXE, FTB_MTM_JITTER.EXE, FTB_TESTLOG.EXE, FTB_JDATE.EXE, -
      FTB_HISTO.EXE, FTB_OVERLAP.EXE, FTB_BINFLD.EXE, -
      FTB_CHECKSUM.EXE, FTB_VALID_TIME.EXE, FTB_CHECK_DQ.EXE, FTB_TRANS.EXE, -
      FTB_LTRANS.EXE
  !Facility FTB is up-to-date.

TEXT_LIBS = CSDR$LIBRARY:FUTLIB.TLB/LIB + -
	    CSDR$LIBRARY:CSDRLIB.TLB/LIB

.FOR.OBJ
    $(FORT) /EXTEND_SOURCE $(FFLAGS) $(MMS$SOURCE) + $(TEXT_LIBS)

!BUILD THE GETTIM.EXE PROGRAM

FTB_GETTIM_OBJ = FTB_GETTIM, FTB_MSG

FTB_GETTIM.EXE : FTBBLD.OLB($(FTB_GETTIM_OBJ)), -
	         CSDR$LIBRARY:FUTLIB.OLB, -
	         CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_GETTIM), -
	               CSDR$LIBRARY:FUTLIB.OLB/LIB, -
                       CSDR$LIBRARY:CSDRLIB/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_GETTIM has been built.
 !

!BUILD THE FTB_WORD31.EXE PROGRAM

FTB_WORD31_OBJ = FTB_WORD31, FTB_FETCH_DATA, FTB_READ_TEMPC, -
		 FTB_READ_WORD31, FTB_GRAPH_WORD31, FTB_FILL_WORD31, -
		 FTB_MSG

FTB_WORD31.EXE : FTBBLD.OLB($(FTB_WORD31_OBJ)), -
		 CSDR$LIBRARY:FUTLIB.OLB, -
		 CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_WORD31), -
	               CSDR$LIBRARY:FUTLIB.OLB/LIB, -
                       CSDR$LIBRARY:CSDRLIB/LIB, -
	 	       CSDR$LIBRARY:CSDRMSG/LIB, -
	 	       XANADU:[LIB]VILIB/LIB, -
		       XANADU:[LIB]XLIB/LIB, -
		       GRAPHICS:GRPSHR/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_WORD31 has been built.
 !

!BUILD THE FTB_TEMPCONTROL PROGRAM.

FTB_TEMPCONTROL_OBJ =	FTB_TEMPCONTROL, FTB_FETCH_DATA, FTB_READ_WORD31, -
                        FTB_READ_TEMPC, FTB_GRAPH_TEMPC, FTB_FILL_TEMPC, -
                        FTB_MSG

FTB_TEMPCONTROL.EXE :	FTBBLD.OLB( $(FTB_TEMPCONTROL_OBJ) ), -
			CSDR$LIBRARY:FUTLIB.OLB, -
	                CSDR$LIBRARY:V5.OPT
 - $(link) $(linkflags)  FTBBLD.OLB/LIB/INC=(FTB_TEMPCONTROL), -
		         CSDR$LIBRARY:FUTLIB.OLB/LIB, -
			 CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		         XANADU:[LIB]VILIB/LIB, -
		         XANADU:[LIB]XLIB/LIB, -
		         GRAPHICS:GRPSHR/LIB, -
	                 CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_TEMPCONTROL has been built.
 !

!BUILD THE FTB_JDATE.EXE PROGRAM

FTB_JDATE_OBJ = FTB_JDATE, FTB_J_TO_G, FTB_G_TO_J, FTB_MSG

FTB_JDATE.EXE : FTBBLD.OLB($(FTB_JDATE_OBJ)), -
		CSDR$LIBRARY:FUTLIB.OLB, -
	        CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_JDATE), -
	               CSDR$LIBRARY:FUTLIB.OLB/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_JDATE has been built.
 !

!BUILD THE CTVIEWER.EXE PROGRAM.

FTB_CTVIEWER_OBJ = FTB_CTVIEWER, FTB_MSG

FTB_CTVIEWER.EXE : FTBBLD.OLB($(FTB_CTVIEWER_OBJ)), -
		   CSDR$LIBRARY:FUTLIB.OLB, -
	           CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_CTVIEWER), -
                       CSDR$LIBRARY:FUTLIB.OLB/LIB, -
                       CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		       XANADU:[LIB]VILIB/LIB, -
		       XANADU:[LIB]XLIB/LIB, -
		       GRAPHICS:GRPSHR/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_CTVIEWER has been built.
 !

!BUILD THE CONVTIME.EXE PROGRAM.

FTB_CONVTIME_OBJ = FTB_CONVTIME, FTB_MSG

FTB_CONVTIME.EXE : FTBBLD.OLB($(FTB_CONVTIME_OBJ)), -
		   CSDR$LIBRARY:FUTLIB.OLB, -
	           CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_CONVTIME), -
	               CSDR$LIBRARY:FUTLIB.OLB/LIB, -
                       CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_CONVTIME has been built.
 !

!BUILD THE MTM_JITTER.EXE PROGRAM.

FTB_MTM_JITTER_OBJ = FTB_MTM_JITTER, FTB_MTM_SUMPLOT, -
		     FTB_PLOTSPEC, FTB_MSG

FTB_MTM_JITTER.EXE : FTBBLD.OLB($(FTB_MTM_JITTER_OBJ)), -
		     CSDR$LIBRARY:FUTLIB.OLB, -
	             CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_MTM_JITTER), -
		       IMSL$DIR:IMSL.OLB/LIB, -
                       CSDR$LIBRARY:FUTLIB.OLB/LIB, -
                       CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
	               XANADU:[LIB]VILIB.OLB/LIB, -
	               XANADU:[LIB]XLIB.OLB/LIB, -
		       GRAPHICS:GRPSHR/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_MTM_JITTER has been built.
 !

!BUILD THE FTB_TESTLOG PROGRAM.

FTB_TESTLOG_OBJ =	FTB_TESTLOG, FTB_MSG

FTB_TESTLOG.EXE :	FTBBLD.OLB($(FTB_TESTLOG_OBJ)), -
			FTB_TESTLOG.OPT, -
			DTR$LIBRARY:TERMSERVE.OLB, -
			CSDR$LIBRARY:FUTLIB.OLB, -
			SYS$SHARE:DTRSHR.EXE, -
			CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_TESTLOG), -
		       FTB_TESTLOG.OPT/OPTION, -
                       CSDR$LIBRARY:FUTLIB.OLB/LIB, -
	 	       CSDR$LIBRARY:CSDRLIB/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_TESTLOG has been built.
 !

!BUILD THE FTB_HISTO PROGRAM.

FTB_HISTO_OBJ =		FTB_HISTO, FTB_ACCUM_HISTO, -
                        FTB_MSG

FTB_HISTO.EXE :		FTBBLD.OLB( $(FTB_HISTO_OBJ) ), -
			CSDR$LIBRARY:FUTLIB.OLB, -
	                CSDR$LIBRARY:V5.OPT
 - $(link) $(linkflags)  FTBBLD.OLB/LIB/INC=(FTB_HISTO), -
		         CSDR$LIBRARY:FUTLIB.OLB/LIB, -
			 CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		         XANADU:[LIB]VILIB/LIB, -
		         XANADU:[LIB]XLIB/LIB, -
		         GRAPHICS:GRPSHR/LIB, -
	                 CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_HISTO has been built.
 !

!Build the FTB_OVERLAP program.

FTB_OVERLAP_OBJ = FTB_OVERLAP.OBJ, -
		  FTB_OVERLAP_INIT.OBJ, -
		  FTB_OVERLAP_SEARCH.OBJ, -
		  FTB_REMOVE_OVERLAP.OBJ, -
		  FTB_OVERLAP_SET_START.OBJ,-
		  FTB_OVERLAP_SET_STOP.OBJ,-
		  FTB_REFINE_OVERLAP.OBJ, -
		  FTB_SORT_CATALOG.OBJ, -
		  FTB_MSG.OBJ

FTB_OVERLAP.EXE :  FTBBLD.OLB($(FTB_OVERLAP_OBJ)), -
		   CSDR$LIBRARY:V5.OPT
 - $(LINK) $(LINKFLAGS) /EXE=FTB_OVERLAP.EXE -
		FTBBLD.OLB/LIB/INC=(FTB_OVERLAP), -
                CSDR$LIBRARY:FUTLIB.OLB/LIB, -
	   	CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_OVERLAP has been built.
 !

!BUILD MESSAGES.

FTB_MSG.OBJ : FTB_MSG.MSG
  MESSAGE FTB_MSG

!DEFINE THE COMPILATION DEPENDENCIES.

FTB_CONVTIME.OBJ	: CSDR$LIBRARY:CTUSER.INC
FTB_CTVIEWER.OBJ	: CSDR$LIBRARY:CTUSER.INC, CSDR$LIBRARY:FUTLIB.TLB,-
			  NFS_SDF^, FNT_NOISE^
FTB_TESTLOG.OBJ		:
FTB_GETTIM.OBJ		: CSDR$LIBRARY:CTUSER.INC, CSDR$LIBRARY:FUTLIB.TLB, -
			  FDQ_ENG^, NFS_HKP^, NFS_SDF^, NFS_EMF^, FDQ_IDX^, -
			  FNT_NOISE^
FTB_JDATE.OBJ           :
FTB_FETCH_DATA.OBJ	: CSDR$LIBRARY:CTUSER.INC
FTB_READ_DATA.OBJ	: CSDR$LIBRARY:CTUSER.INC
FTB_HISTO.OBJ		: CT$LIBRARY:CTUSER.INC, NFS_SDF^
FTB_ACCUM_HISTO.OBJ	: CT$LIBRARY:CTUSER.INC, CSDR$LIBRARY:FUTLIB.TLB,-
                          NFS_SDF^
FTB_MTM_JITTER.OBJ      : CSDR$LIBRARY:FUTLIB.TLB,-
			  NFS_EMF^
FTB_PLOTSPEC.OBJ        :
FTB_MTM_SUMPLOT.OBJ     :
FTB_FILL_TEMPC.OBJ	: NFS_HKP^
FTB_FILL_WORD31.OBJ	: NFS_ANC^
FTB_GRAPH_WORD31.OBJ	:
FTB_READ_TEMPC.OBJ	: CSDR$LIBRARY:CTUSER.INC, NFS_HKP^
FTB_READ_WORD31.OBJ	: CSDR$LIBRARY:CTUSER.INC, NFS_ANC^

!Build the FTB_BINFLD program.

FTB_BINFLD_OBJ = FTB_BINFLD.OBJ, -
		  FTB_PARSE_BINFLD.OBJ, -
		  FTB_FLD_IDX.OBJ, -
		  FTB_GET_FLDVAL.OBJ,-
		  FTB_GET_MTM_MODE.OBJ,-
		  FTB_PRINT_ENG_INFO.OBJ,-
		  FTB_ACCUM_ENG_TABLE.OBJ,-
		  FTB_PRINT_ENG_TABLE.OBJ,-
		  FTB_MSG.OBJ

FTB_BINFLD.EXE :  FTBBLD.OLB($(FTB_BINFLD_OBJ)), -
	           CSDR$LIBRARY:V5.OPT
 - $(LINK) $(LINKFLAGS) /EXE=FTB_BINFLD.EXE -
		FTBBLD.OLB/LIB/INC=(FTB_BINFLD), -
                CSDR$LIBRARY:FUTLIB.OLB/LIB, -
	   	CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		CSDR$LIBRARY:V5.OPT/OPTION
!
! The "Multiply defined symbol WAIT" warning may be ignored.
!
!FTB_BINFLD has been built.

!Build the FTB_CHECKSUM program.

FTB_CHECKSUM_OBJ = FTB_CHECKSUM.OBJ, -
		  FTB_PARSE_CHECKSUM.OBJ, -
		  FTB_INIT_CHECK_REPORT.OBJ, -
		  FTB_PRINT_SCI_INFO.OBJ,-
		  FTB_PRINT_BADSCI_INFO.OBJ,-
		  FTB_MSG.OBJ

FTB_CHECKSUM.EXE :  FTBBLD.OLB($(FTB_CHECKSUM_OBJ)), -
	           CSDR$LIBRARY:V5.OPT
 - $(LINK) $(LINKFLAGS) /EXE=FTB_CHECKSUM.EXE -
		FTBBLD.OLB/LIB/INC=(FTB_CHECKSUM), -
                CSDR$LIBRARY:FUTLIB.OLB/LIB, -
	   	CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		CSDR$LIBRARY:V5.OPT/OPTION
!
! The "Multiply defined symbol WAIT" warning may be ignored.
!
!FTB_CHECKSUM has been built.

!Build the FTB_VALID_TIME program.

FTB_VALID_TIME_OBJ = FTB_VALID_TIME.OBJ, -
		  FTB_PARSE_VALTIME.OBJ, -
		  FTB_INIT_REPORT.OBJ, -
		  FTB_CHECK_SCI_TIME.OBJ,-
		  FTB_CHECK_PRE_IT_SCI_TIME.OBJ,-
		  FTB_LIST_SCI_TIME.OBJ,-
		  FTB_CHECK_ENG_TIME.OBJ,-
		  FTB_CHECK_HKP_TIME.OBJ,-
		  FTB_DUMP_ENG_TIME.OBJ,-
     	          FTB_DUMP_HKP_TIME.OBJ,-
                  FTB_DUMP_SCI_DATA.OBJ,-
                  FTB_DUMP_SCI_TIME.OBJ,-
                  FTB_IF_MIDPOINT_TIME.OBJ,-
                  FTB_SUM_IFG.OBJ,-
		  FTB_MSG.OBJ

FTB_VALID_TIME.EXE :  FTBBLD.OLB($(FTB_VALID_TIME_OBJ)), -
	           CSDR$LIBRARY:V5.OPT
 - $(LINK) $(LINKFLAGS) /EXE=FTB_VALID_TIME.EXE -
		FTBBLD.OLB/LIB/INC=(FTB_VALID_TIME), -
                CSDR$LIBRARY:FUTLIB.OLB/LIB, -
	   	CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		CSDR$LIBRARY:V5.OPT/OPTION
!
! The "Multiply defined symbol WAIT" warning may be ignored.
!
!FTB_VALID_TIME has been built.

!Build the FTB_CHECK_DQ program.

FTB_CHECK_DQ_OBJ = FTB_CHECK_DQ.OBJ, -
		  FTB_LIST_SCI.OBJ, -
		  FTB_LIST_ENGSTAT.OBJ, -
		  FTB_MSG.OBJ

FTB_CHECK_DQ.EXE :  FTBBLD.OLB($(FTB_CHECK_DQ_OBJ)), -
	           CSDR$LIBRARY:V5.OPT
 - $(LINK) $(LINKFLAGS) /EXE=FTB_CHECK_DQ.EXE -
		FTBBLD.OLB/LIB/INC=(FTB_CHECK_DQ), -
                CSDR$LIBRARY:FUTLIB.OLB/LIB, -
	   	CSDR$LIBRARY:CSDRLIB.OLB/LIB, -
		CSDR$LIBRARY:V5.OPT/OPTION
!
! The "Multiply defined symbol WAIT" warning may be ignored.
!
!FTB_CHECK_DQ has been built.

!BUILD THE TRANS.EXE PROGRAM

FTB_TRANS_OBJ = FTB_TRANS, FTB_MSG

FTB_TRANS.EXE : FTBBLD.OLB($(FTB_TRANS_OBJ)), -
	         CSDR$LIBRARY:FUTLIB.OLB, -
	         CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_TRANS), -
	               CSDR$LIBRARY:FUTLIB.OLB/LIB, -
                       CSDR$LIBRARY:CSDRLIB/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_TRANS has been built.
 !

!BUILD THE LTRANS.EXE PROGRAM

FTB_LTRANS_OBJ = FTB_LTRANS, FTB_MSG

FTB_LTRANS.EXE : FTBBLD.OLB($(FTB_LTRANS_OBJ)), -
	         CSDR$LIBRARY:FUTLIB.OLB, -
	         CSDR$LIBRARY:V5.OPT
 - $(link)$(linkflags) FTBBLD.OLB/LIB/INC=(FTB_LTRANS), -
	               CSDR$LIBRARY:FUTLIB.OLB/LIB, -
                       CSDR$LIBRARY:CSDRLIB/LIB, -
	               CSDR$LIBRARY:V5.OPT/OPTION
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 !FTB_LTRANS has been built.
 !
