!============================================================================
!
!	[FREF]DESCRIP.MMS
!
! Purpose:  Builds the FIRAS reference data files.
!
!	There are three sources for FIRAS reference dataset handled by
!	this one MMS file:
!
!	    (1) Create the reference dataset from a FIRAS FRD facility
!		program.  This sometimes requires an input text file.
!		If so, the text file must be in the directory from which
!		the MMS file is run.
!
!	    (2) Receive a user supplied record selection expression for
!		the RSE reference datasets (FEX_FCIRSE.DAT, FEX_FECRSE.TXT,
!		FEX_FESRSE.DAT, FEX_FFCRSE.DAT, FEX_FPRRSE.DAT, FEX_FTBRSE.DAT).
!
!	    (3) Take a snapshot of the pertinent STOL database files from
!		DB$.  These are the FDB reference datasets (FDB_GRTCAL.DAT,
!		FDB_GRTRICH.DAT, and FDB_LIMCUVKY.DAT).
!
!	All input text files required for this MMS have the extension TXT.
!	The text file must reside in the directory from which the MMS is
!	performed.
!
!	All reference files created by this MMS have the extension DAT.  The
!	INSTALL places the DAT files in CSDR$FIRAS_REFERENCE.  It is set
!	CSDR$FIRAS_REFERENCE prior to running this MMS.
!
!	For case three, logical DB$ defines the STOL database from which the
!	snapshot is taken.  Therefore it is important to have DB$ properly
!	defined prior to running this MMS.
!
!
! Author:   R. Kummerer
!           STX
!           September 30, 1989
!
! Modification Log:
! 
!	  Author	  Date		Modification
!	  ------	  ----		------------
!         H. Wang         3/7/90        Add FEX_CALRES.DAT
!         S. Read         9/13/91       Modify structure to allow all target
!                                       reference files to be built as the 
!                                       default or specify one or more target
!                                       files on MMS invocation. Add new 
!                                       MTM_Sweep and VAB_SAA reference data 
!                                       sets. Add the time ordered reference 
!                                       files: gain, fakeit and average 
!                                       calibrator resistors, created from the 
!                                       mission raw housekeeping and ancillary 
!                                       archive to this Descrip.MMS file. 
!                                       Also, add to the ".FIRST" block the 
!                                       additional archive definitions needed 
!                                       by the FRD programs which create these 
!                                       time ordered files.
!
!                                       Note: all the reference files will NOT
!                                             be delivered for Phase 1. Targets 
!                                             for these files will have to be 
!                                             delivered later.
!
!                                     Note 2: FRD_GFG and FRD_AV_CAL_RES should
!                                             NOT be run unless requested. They
!                                             are separate targets and must be
!                                             specified as targets along with 
!                                             FREF if they are to be built.
!	Modified: 1991 4 October - Rick Shafer
!		Changed the number of runs for the CALRES to 4 time ranges
!		to include the post cryogen runout period.
!       Modified: 1991 9 October - Steve Alexander
!               FRD programs that create a .DAT file from a .TXT file will fail
!               unless the .TXT file is in the physical directory from which
!               the FRD program is run.  This is not true of search lists.
!               Add a COPY of the .TXT file before running the program.
!	Certified: 1991 9 October - Rick Shafer Pass 2 Phase FPPFDQ
!
!       Modified: 1991 22 November - Steve Alexander
!               Add new target FEX_ENGLIM.DAT_90165_90194 to accomodate fix
!               to DB$:FIRLIMS.INP file and reprocess this time range of raw
!               data.
!
!	Certified: 1991 22 November - Rick Shafer Pass 2 Phase FPPFDQ
!                                     Reprocessing of Time Range 90165-90194
!
!       Modified: 1991 12 December - Steve Alexander SER 9339
!               Add new target FEX_MINCOADD.DAT as reference dataset input to
!               FEC and FSS processing.
!	Certified:  1991 Dec 12 - Rick Shafer Pass 2 Phase FSS
!					MINCOADD reference dataset
!
!       Modified: 1992 5 August - Steve Alexander SER 9865
!               Add new target FEX_SAMPRATE.DAT as reference dataset input to
!               FDQ processing; remove temporary target
!               FEX_ENGLIM.DAT_90165_90194.
!	Certified: 1992 Aug 5 - Rick Shafer Pass 2A Phase FPPFDQ,FSS,FEC
!
!	Modified: 1993 11 February - Steve Alexander SER 10565
!		Add new targets FEX_BASIS.DAT, FEX_CMDGAIN.DAT, FEX_CTH.DAT,
!		FEX_DTRF.DAT, FEX_ETF.DAT, FEX_GLTCHPRO.DAT,
!		FEX_GRTCOAWT.DAT, FEX_NYQUIST.DAT, FEX_REFTEMPS.DAT as
!		reference dataset inputs to FIC processing.
!	Certified: 1992 feb 22 - Rick Shafer Pass 2a Phase FIC sky
!
!	Modified: 1993 5 April - Steve Alexander SER 8836 
!		Add new target FEX_APOD.DAT as reference dataset input 
!		to FFI processing.
!	Certified: 1993 Apr 6 - John Mather Pass 2a Phase FFI
!
!	Modified: 1993 18 May - Steve Alexander SER 10954
!		Add new targets FEX_MOD_XXYY.F16_93HYBRID as reference 
!               datasets for calibration models for input to FCF processing.
!	Certified: 1993 May 18 - John Mather Pass 2a Phase FCF
!
!	Modified: 1993 24 May - Steve Alexander SER 8292
!		Add new target FEX_VIBCORR.DAT as reference dataset for input 
!               to FCF processing.
!	Certified: 1993 May 24 - John Mather Pass 2a Phase FCF
!
!	Modified: 1994 4 April - Steve Alexander SER 11527
!		Add new targets FEX_MOD_RLSF.F16_93HYBRID, 
!               FEX_MOD_LLSF.F16_93HYBRID, FEX_MOD_RLFL.F16_93HYBRID,
!               FEX_MOD_LLFL.F16_93HYBRID  as reference datasets for the final
!               calibration models.
!	Certified: 1994 April 4 - John Mather
!
!	Modified: 1994 18 April - Steve Alexander SER 11702
!		Add new target FEX_GLTCHCOR.DAT as reference dataset for input
!               to FDS processing.
!	Certified: 1994 April 18 - John Mather Pass 3 Phase FDS
!
!	Modified: 1995 28 August - Steve Brodd SER 12244
!		Add new targets FEX_APODL.DAT, FEX_ETFL.DAT, and 
!               FEX_NYQUISTL.DAT as reference datasets for input
!               to FIL processing.
!	Certified: 1995 August 28 - John Mather Pass 4 Phase FIL
!
!	Modified: 1996 2 Jan - Steve Brodd SPR 12287
!		Add new target FEX_VIBCORRL.DAT as reference dataset for 
!               input to FSL processing.
!
!	Modified: 1996 2 Jan - Steve Brodd SPR 12282
!		Add new targets FEX_MOD_XXYY.PASS4 as reference 
!               datasets for calibration models for input to FSL processing.
!============================================================================

.FIRST
	SET COMMAND CSDR$CLD:FRD
	@ WRITE SYS$OUTPUT "Note: FDB_*.DAT reference files taken from:"
	@ WRITE SYS$OUTPUT "      DB$     = ''F$TRNLNM(""DB$"")'"
	@ WRITE SYS$OUTPUT "      DB$ROOT = ''F$TRNLNM(""DB$ROOT"")'"
        @ Write SYS$OUTPUT " "
        DEFINE CSDR$FIRAS_RAW CSDR$FIRAS_EDIT
        DEFINE CSDR$FIRAS_IN []
        DEFINE CSDR$FIRAS_OUT []
        DEFINE CSDR$FIRAS_REF CSDR$FIRAS_REFERENCE
        @ Write SYS$OUTPUT "CSDR$FIRAS_RAW=''F$TRNLNM(""CSDR$FIRAS_RAW"")'"
        @ Write SYS$OUTPUT "CSDR$FIRAS_IN=''F$TRNLNM(""CSDR$FIRAS_IN"")'"
        @ Write SYS$OUTPUT "CSDR$FIRAS_OUT=''F$TRNLNM(""CSDR$FIRAS_OUT"")'"
        @ Write SYS$OUTPUT "CSDR$FIRAS_REF=''F$TRNLNM(""CSDR$FIRAS_REF"")'"
        @ Write SYS$OUTPUT " "

FREF : -
	FDB_GRTCAL.DAT        ,-
	FDB_GRTRICH.DAT       ,-
	FDB_LIMCUVKY.DAT      ,-
	FEX_MTMSWEEP.DAT      ,-	
	FEX_VABSAA.DAT        ,-
	FEX_IDX_FLAG.DAT      ,-
	FEX_IDX_TOLS.DAT      ,-
	FEX_ENGLIM.DAT        ,-
	FEX_SCILIM.DAT        ,-
	FEX_LIMFLAGS.DAT      ,-
	FEX_GRTTRANS.DAT      ,-
	FEX_GRTRAWWT.DAT      ,-
	FEX_CALRES.DAT        ,-
	FEX_SAMPRATE.DAT      ,-
        FEX_MINCOADD.DAT      ,-
	FEX_BASIS.DAT         ,-
	FEX_CMDGAIN.DAT       ,-
	FEX_CTH.DAT           ,-
	FEX_DTRF.DAT          ,-
	FEX_ETF.DAT           ,-
	FEX_ETFL.DAT           ,-
	FEX_GLTCHPRO.DAT      ,-
	FEX_GRTCOAWT.DAT      ,-
	FEX_NYQUIST.DAT       ,-
	FEX_NYQUISTL.DAT       ,-
	FEX_REFTEMPS.DAT      ,-
	FEX_APOD.DAT          ,-
	FEX_APODL.DAT          ,-
	FEX_VIBCORR.DAT	      ,-
	FEX_VIBCORRL.DAT       ,-
	FEX_GLTCHCOR.DAT      ,-
        FEX_MOD_RHSS.F16_93HYBRID,-
        FEX_MOD_RHSF.F16_93HYBRID,-
        FEX_MOD_RHLF.F16_93HYBRID,-
        FEX_MOD_LHSS.F16_93HYBRID,-
        FEX_MOD_LHSF.F16_93HYBRID,-
        FEX_MOD_LHLF.F16_93HYBRID,-
        FEX_MOD_RLSS.F16_93HYBRID,-
        FEX_MOD_RLSF.F16_93HYBRID,-
        FEX_MOD_RLLF.F16_93HYBRID,-
        FEX_MOD_RLFL.F16_93HYBRID,-
        FEX_MOD_LLSS.F16_93HYBRID,-
        FEX_MOD_LLSF.F16_93HYBRID,-
        FEX_MOD_LLLF.F16_93HYBRID,-
        FEX_MOD_LLFL.F16_93HYBRID,-
        FEX_MOD_RHSS.PASS4,-
        FEX_MOD_RHSF.PASS4,-
        FEX_MOD_RHLF.PASS4,-
        FEX_MOD_LHSS.PASS4,-
        FEX_MOD_LHSF.PASS4,-
        FEX_MOD_LHLF.PASS4,-
        FEX_MOD_RLSS.PASS4,-
        FEX_MOD_RLLF.PASS4,-
        FEX_MOD_RLFS.PASS4,-
        FEX_MOD_RLFL.PASS4,-
        FEX_MOD_LLSS.PASS4,-
        FEX_MOD_LLLF.PASS4,-
        FEX_MOD_LLFS.PASS4,-
        FEX_MOD_LLFL.PASS4
 ! The FIRAS Reference Data files have been built (except FGFG and FACR).

!
! Calibrator Resistor Counts
!
FDB_GRTCAL.DAT		: DB$:GRTCAL.DB
       COPY $(MMS$SOURCE) $(MMS$TARGET);0
 ! FDB_GRTCAL.DAT has been built.

!
! GRT Conversion Tables
!
FDB_GRTRICH.DAT		: DB$:GRTRICH.DB
       COPY $(MMS$SOURCE) $(MMS$TARGET);0
 ! FDB_GRTRICH.DAT has been built.

!
! Engineering Conversion Tables
!
FDB_LIMCUVKY.DAT	: DB$:LIMCUVKEY.DB
       COPY $(MMS$SOURCE) $(MMS$TARGET);0
 ! FDB_LIMCUVKY.DAT has been built.

!
! MTM Sweep Times
!
FEX_MTMSWEEP.DAT	: FEX_MTMSWEEP.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
	MTMSWEEP
 ! FEX_MTMSWEEP.DAT has been built.

!
! VAB and SAA Boundaries
!
VABSAA_TXT = FEX_VABN.TXT, FEX_VABS.TXT, FEX_SAA.TXT

FEX_VABSAA.DAT		: $(VABSAA_TXT)
        COPY $(VABSAA_TXT) []
	VABSAA
 ! FEX_VABSAA.DAT has been built.

!
! Engineering Index Flags
!
FEX_IDX_FLAG.DAT	: FEX_IDX_FLAG.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        IDX_PARAMS/INDEX=FLAGS/REPORT
 ! FEX_IDX_FLAG.DAT has been built.

!
! Engineering Index Tolerances
!
FEX_IDX_TOLS.DAT	: FEX_IDX_TOLS.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        IDX_PARAMS/INDEX=TOLS/REPORT
 ! FEX_IDX_TOLS.DAT has been built.

!
! Engineering Limits
!
FEX_ENGLIM.DAT		: FEX_GRTLIM.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LIMITS/ENG
 ! FEX_ENGLIM.DAT has been built.

!
! Science Limits
!
FEX_SCILIM.DAT		: FEX_SCILIM.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LIMITS/SCI
 ! FEX_SCILIM.DAT has been built.

!
! Engineering Limits Flags
!
FEX_LIMFLAGS.DAT	: FEX_LIMFLAGS.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LIMITS/FLAG
 ! FEX_LIMFLAGS.DAT has been built.

!
! GRT Low to High Current Transition Temperatures and Half-Widths
!
FEX_GRTTRANS.DAT	: FEX_GRTTRANS.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        GTRAN/OUTPUT=TRAN
 ! FEX_GRTTRANS.DAT has been built.

!
! GRT Weights for Combining A and B GRT Temperatures for Raw Data Input
!
FEX_GRTRAWWT.DAT	: FEX_GRTRAWWT.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        GTRAN/OUTPUT=RAW
 ! FEX_GRTRAWWT.DAT has been built.

!
! Dwellplot calibration resistance counts
!
FEX_CALRES.DAT		: FEX_CALRES.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        CALRES
 ! FEX_CALRES.DAT has been built.

!
! Sampling rate
!
FEX_SAMPRATE.DAT		: FEX_SAMPRATE.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        SAMPRATE
 ! FEX_SAMPRATE.DAT has been built.

!
! Minimum IFGs needed to coadd
!
FEX_MINCOADD.DAT      	: FEX_MINCOADD.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        MINCOADD
 ! FEX_MINCOADD.DAT has been built.

!
! Baseline subtraction basis vectors
!
FEX_BASIS.DAT      	: 
        BASIS
 ! FEX_BASIS.DAT has been built.

!
! Commanded gains
!
FEX_CMDGAIN.DAT      	: FEX_CMDGAIN.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        CMDGAIN
 ! FEX_CMDGAIN.DAT has been built.

!
! Consistency thresholds
!
FEX_CTH.DAT      	: FEX_CTH.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        CCTHRESHOLDS
 ! FEX_CTH.DAT has been built.

!
! Digital transient subtraction functions
!
FEX_DTRF.DAT      	: 
        DTRANS
 ! FEX_DTRF.DAT has been built.

!
! Electronics transfer functions
!
FEX_ETF.DAT      	: FEX_SAMPRATE.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        ELEX_TRANSFCN
 ! FEX_ETF.DAT has been built.

!
! Long electronics transfer functions
!
FEX_ETFL.DAT      	: FEX_SAMPRATE.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LELEX_TRANSFCN
 ! FEX_ETFL.DAT has been built.

!
! Glitch profiles
!
FEX_GLTCHPRO.DAT      	: 
        GPROFILE
 ! FEX_GLTCHPRO.DAT has been built.

!
! GRT Weights for Combining A and B GRT Temperatures for Coadd Data Input
!
FEX_GRTCOAWT.DAT	: FEX_GRTCOAWT.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        GTRAN/OUTPUT=COAD
 ! FEX_GRTCOAWT.DAT has been built.

!
! Nyquist frequencies
!
NYQUIST_TXT = FEX_NYQUIST.TXT, FEX_SAMPRATE.TXT

FEX_NYQUIST.DAT		: $(NYQUIST_TXT)
        COPY $(NYQUIST_TXT) []
	NYQUIST
 ! FEX_NYQUIST.DAT has been built.

!
! Long Nyquist frequencies
!
FEX_NYQUISTL.DAT		: $(NYQUIST_TXT)
        COPY $(NYQUIST_TXT) []
	LNYQUIST
 ! FEX_NYQUISTL.DAT has been built.

!
! Reference templates
!
FEX_REFTEMPS.DAT      	: FEX_REFTEMPS.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        REFTEMPS
 ! FEX_REFTEMPS.DAT has been built.

!
! Apodization functions
!
FEX_APOD.DAT      	: 
        APOD
 ! FEX_APOD.DAT has been built.

!
! Long apodization functions
!
FEX_APODL.DAT      	: 
        LAPOD
 ! FEX_APODL.DAT has been built.

!
! Vibration correction
!
FEX_VIBCORR.DAT	: FEX_VIBCORR.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        VIBCORR
 ! FEX_VIBCORR.DAT has been built.

!
! Long vibration correction
!
FEX_VIBCORRL.DAT	: FEX_VIBCORRL.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LVIBCORR
 ! FEX_VIBCORRL.DAT has been built.

!
! Glitch rate correction
!
FEX_GLTCHCOR.DAT : FEX_GLTCHCOR.TXT
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        GCORR
 ! FEX_GLTCHCOR.DAT has been built.

!
! Calibration model RHSS
!
FEX_MOD_RHSS.F16_93HYBRID : FEX_MOD_RHSS.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=RH/SCAN=SS
 ! FEX_MOD_RHSS.F16_93HYBRID has been built.

!
! Calibration model RHSF
!
FEX_MOD_RHSF.F16_93HYBRID : FEX_MOD_RHSF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=RH/SCAN=SF
 ! FEX_MOD_RHSF.F16_93HYBRID has been built.

!
! Calibration model RHLF
!
FEX_MOD_RHLF.F16_93HYBRID : FEX_MOD_RHLF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=RH/SCAN=LF
 ! FEX_MOD_RHLF.F16_93HYBRID has been built.

!
! Calibration model LHSS
!
FEX_MOD_LHSS.F16_93HYBRID : FEX_MOD_LHSS.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=LH/SCAN=SS
 ! FEX_MOD_LHSS.F16_93HYBRID has been built.

!
! Calibration model LHSF
!
FEX_MOD_LHSF.F16_93HYBRID : FEX_MOD_LHSF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=LH/SCAN=SF
 ! FEX_MOD_LHSF.F16_93HYBRID has been built.

!
! Calibration model LHLF
!
FEX_MOD_LHLF.F16_93HYBRID : FEX_MOD_LHLF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=LH/SCAN=LF
 ! FEX_MOD_LHLF.F16_93HYBRID has been built.

!
! Calibration model RLSS
!
FEX_MOD_RLSS.F16_93HYBRID : FEX_MOD_RLSS.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=RL/SCAN=SS
 ! FEX_MOD_RLSS.F16_93HYBRID has been built.

!
! Calibration model RLSF
!
FEX_MOD_RLSF.F16_93HYBRID : FEX_MOD_RLSF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=RL/SCAN=SF
 ! FEX_MOD_RLSF.F16_93HYBRID has been built.

!
! Calibration model RLLF
!
FEX_MOD_RLLF.F16_93HYBRID : FEX_MOD_RLLF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=RL/SCAN=LF
 ! FEX_MOD_RLLF.F16_93HYBRID has been built.

!
! Calibration model RLFL
!
FEX_MOD_RLFL.F16_93HYBRID : FEX_MOD_RLFL.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=RL/SCAN=FL
 ! FEX_MOD_RLFL.F16_93HYBRID has been built.

!
! Calibration model LLSS
!
FEX_MOD_LLSS.F16_93HYBRID : FEX_MOD_LLSS.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=LL/SCAN=SS
 ! FEX_MOD_LLSS.F16_93HYBRID has been built.

!
! Calibration model LLSF
!
FEX_MOD_LLSF.F16_93HYBRID : FEX_MOD_LLSF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=LL/SCAN=SF
 ! FEX_MOD_LLSF.F16_93HYBRID has been built.

!
! Calibration model LLLF
!
FEX_MOD_LLLF.F16_93HYBRID : FEX_MOD_LLLF.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=LL/SCAN=LF
 ! FEX_MOD_LLLF.F16_93HYBRID has been built.

!
! Calibration model LLFL
!
FEX_MOD_LLFL.F16_93HYBRID : FEX_MOD_LLFL.TXT_F16_93HYBRID
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        EXTRACT_MODEL/FLIGHT/FILE=F16_93HYBRID/CHAN=LL/SCAN=FL
 ! FEX_MOD_LLFL.F16_93HYBRID has been built.

!
! Long calibration model RHSS
!
FEX_MOD_RHSS.PASS4 : FEX_MOD_RHSS.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=RH/SCAN=SS
 ! FEX_MOD_RHSS.PASS4 has been built.

!
! Long calibration model RHSF
!
FEX_MOD_RHSF.PASS4 : FEX_MOD_RHSF.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=RH/SCAN=SF
 ! FEX_MOD_RHSF.PASS4 has been built.

!
! Long calibration model RHLF
!
FEX_MOD_RHLF.PASS4 : FEX_MOD_RHLF.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=RH/SCAN=LF
 ! FEX_MOD_RHLF.PASS4 has been built.

!
! Long calibration model LHSS
!
FEX_MOD_LHSS.PASS4 : FEX_MOD_LHSS.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=LH/SCAN=SS
 ! FEX_MOD_LHSS.PASS4 has been built.

!
! Long calibration model LHSF
!
FEX_MOD_LHSF.PASS4 : FEX_MOD_LHSF.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=LH/SCAN=SF
 ! FEX_MOD_LHSF.PASS4 has been built.

!
! Long calibration model LHLF
!
FEX_MOD_LHLF.PASS4 : FEX_MOD_LHLF.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=LH/SCAN=LF
 ! FEX_MOD_LHLF.PASS4 has been built.

!
! Long calibration model RLSS
!
FEX_MOD_RLSS.PASS4 : FEX_MOD_RLSS.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=RL/SCAN=SS
 ! FEX_MOD_RLSS.PASS4 has been built.

!
! Long calibration model RLLF
!
FEX_MOD_RLLF.PASS4 : FEX_MOD_RLLF.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=RL/SCAN=LF
 ! FEX_MOD_RLLF.PASS4 has been built.

!
! Long calibration model RLFS
!
FEX_MOD_RLFS.PASS4 : FEX_MOD_RLFS.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=RL/SCAN=FS
 ! FEX_MOD_RLFS.PASS4 has been built.

!
! Long calibration model RLFL
!
FEX_MOD_RLFL.PASS4 : FEX_MOD_RLFL.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=RL/SCAN=FL
 ! FEX_MOD_RLFL.PASS4 has been built.

!
! Long calibration model LLSS
!
FEX_MOD_LLSS.PASS4 : FEX_MOD_LLSS.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=LL/SCAN=SS
 ! FEX_MOD_LLSS.PASS4 has been built.

!
! Long calibration model LLLF
!
FEX_MOD_LLLF.PASS4 : FEX_MOD_LLLF.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=LL/SCAN=LF
 ! FEX_MOD_LLLF.PASS4 has been built.

!
! Long calibration model LLFS
!
FEX_MOD_LLFS.PASS4 : FEX_MOD_LLFS.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=LL/SCAN=FS
 ! FEX_MOD_LLFS.PASS4 has been built.

!
! Long calibration model LLFL
!
FEX_MOD_LLFL.PASS4 : FEX_MOD_LLFL.TXT_PASS4
        COPY $(MMS$SOURCE) $(MMS$SOURCE)
        LEXT/FLIGHT/FILE=PASS4/CHAN=LL/SCAN=FL
 ! FEX_MOD_LLFL.PASS4 has been built.

!
! Gain and Fakeit Status Changes and Telemetry Gaps in the Housekeeping Data.
! Must have all NFS_HKP files for the mission on line.
!
FRD_GFG :
 FGFG/JSTART=89324000000000/JSTOP=90264121043671
 !FRD_GFG is completed.

!
! Average Calibrator Resistor Values For Specified Time Ranges Over the Mission.
! Must have all NFS_HKP files for the mission on line.
!                     (See SER 8937)
FRD_AVE_CALRES :
 FACR/JSTART=89324000000000/JSTOP=90006222822554
 FACR/JSTART=90006222822554/JSTOP=90141094958599
 FACR/JSTART=90141094958599/JSTOP=90264130000000
 FACR/JSTART=90264160000000/JSTOP=90282150807000
 !FRD_AVE_CALRES is completed.
