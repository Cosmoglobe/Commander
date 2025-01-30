	PROGRAM FGA
 
C-------------------------------------------------------------------------
C
C	PROGRAM NAME:
C	  FGA_GET_ARCHIVE
C
C	PROGRAM DESCRIPTION:
C	  This program will inquire a user for an interval of FIRAS
C	  archive records, and then reads and formats the records
C	  read in through COBETRIEVE.
C	AUTHOR:
C	  E.FUNG
C	  GSFC
C	  October 31, 1985
C
C-------------------------------------------------------------------------
C         This program is modified to give users capability to dump
C         the FIRAS archive records either in batch mode, or in interactive
C         mode as default.
C    
C       Modified by: Nilo G. Gonzales/Hughes STX, December 17, 1991.
C         Ref. SPR 9355.
C
C       - Example of Interactive process; default mode
C         $Type FGA, then hit return key 
C
C       - Example of batch mode command format:
C         FGA/NOINT/DATASET=FDQ_SDF_RH/JSTART=YYDDDTTHHSS/JSTOP=YYDDTTHHSS/NORSE
C          
C    PDL:
C       Begin
C            Set clear screen
C            Set Status to success
C            Establish condition handler and display program banner
C 
C            Get operational mode; default is "INTERACTIVE"
C	     Interactive = FAC_Present
C	     If ((UPM_Present('INTERACTIVE') .EQ. UPM_Pres) .OR.
C	*       (UPM_Present('INTERACTIVE') .EQ. UPM_Defaulted)) Then
C	        Set Interactive to FAC_Present
C	     Elseif (UPM_Present('INTERACTIVE') .EQ. UPM_Negated) Then
C	        Set Interactive FAC_Not_Present
C	     Endif
C
C       If batch mode, get the DAFS dataset name.
C            
C            If (UPM_Present('DATASET')) Then
C               Restat = UPM_Get_Value for Dataset
C            Endif
C
C            If (UPM_Present('JSTART')) Then
C               Restat = UPM_Get_Value for Jstart
C            Endif
C
C            If (UPM_Present('JSTOP')) Then
C               Restat = UPM_Get_Value for Jstop
C            Endif
C          
C            Call CT_Init ;initialize COBETRIEVE
C 
C       Initialize for NoRse setting
C       If batch mode, get name of the RSE file.
C            If (UPM_Present('RSE')) then
C               Restat = UPM_Get_Value for RSE file; default FEX_FTBRSE
C            Endif
C           
C            If (Interactive .Eq. FAC_Present) Then
C               Call CT_Query_Bld
C            Else        ;batch mode
C               Restat = CUT_Read_Dafs_record_Igse_Key
C            Endif
C            Check Dataset_Id against FIRAS database names 
C            Get Lun (Write_Lun)
C            Open New file (Unit=Write_Lun)
C            Write translation of Logical names
C            Open archive file(s)
C            Call CT_Query_Arcv
C            
C            If Process is Normal, Then
C               Do While (process is normal for each datasets)
C                  If (Datatype .Eq. FIRAS database_name(s)) Then
C                     Restat = FUT_Qget_XX ;subroutine to read archive data
C                     If (Restat .Eq. Normal) Then
C                         Record = Record + 1
C                         Call FGA_List_XX ;subroutine to dump FIRAS dataset(s)
C                     Else
C                         Set flag "Normal" to false
C                     Endif
C                  Endif
C               Enddo
C            Endif
C            
C            Call CT_Close_Arcv ; close COBETRIEVE archive
C            Close (Write_Lun, iostat = istat)
C       End 
C    END PDL:
C             
C	CALLING SEQUENCE:
C	  This is a main program
C
C	INPUT/OUTPUT PARAMETERS:
C	  NONE
C
C	INPUT FILES:
C	  FIRAS ARCHIVE
C
C	OUTPUT FILES:
C	  List files.
C
C	INCLUDE FILES USED:
C	  CTUSER.INC
C	  FUT_PARAMS
C	  CCT_STATUS_RECORD
C         UPM_STAT_MSG
C	  $SSDEF
C
C	SUBROUTINES CALLED:
C	  FGA_LIST_EDF
C	  FGA_LIST_ENG
C	  FGA_LIST_HKP
C	  FGA_LIST_IDX
C	  FGA_LIST_SCI
C	  FGA_LIST_EMF
C	  FGA_LIST_ENGLIM
C	  FGA_LIST_LIMFLAGS
C	  FGA_LIST_SCILIM
C	  FGA_LIST_IDXTOLS
C	  FGA_LIST_IDXFLAGS
C         FGA_LIST_FACR
C         FGA_LIST_FAKER
C         FGA_LIST_FAKEL
C         FGA_LIST_GAINR
C         FGA_LIST_GAINL
C         FGA_LIST_GTRAN
C         FGA_LIST_GRAWT
C         FGA_LIST_MTMSWEEP
C         FGA_LIST_VABSAA
C         FGA_LIST_MINCOADD
C	  FUT_QGET_ENG
C	  FUT_QGET_HKP
C	  FUT_QGET_HSK
C	  FUT_QGET_IDX
C	  FUT_QGET_SCI
C	  FUT_QGET_EMF
C	  FUT_QGET_NOS
C         FUT_QGET_LIMFLAGS
C         FUT_QGET_ENGLIM
C         FUT_QGET_SCILIM
C         FUT_QGET_IDXTOLS
C         FUT_QGET_IDXFLAGS
C         FUT_QGET_CTH
C         FUT_QGET_FACR
C         FUT_QGET_FAKER
C         FUT_QGET_FAKEL
C         FUT_QGET_GAINR
C         FUT_QGET_GAINL
C         FUT_QGET_GTRAN
C         FUT_QGET_GRAWT
C	  FUT_QGET_MTMSWEEP
C	  FUT_QGET_VABSAA
C         FUT_QGET_MINCOADD
C	  LIB$GET_LUN
C	  LIB$FREE_LUN
C	  STR$UPCASE
C
C---------------------------------------------------------------------------
C
C CHANGES:
C
C	R. Kummerer, August 25, 1987.  Dump the latest IDX and HKP records.
C	R. Kummerer, August 31, 1987.  Dump the latest ENG records.
C	R. Kummerer, November 13, 1987.  Convert to standards.
C	R. Kummerer, December 10, 1987. SPR 1785, Dump model covariance and
C		phase-corrected spectra residuals.
C	R. Kummerer, December 17, 1987. SPR 1796, FUT_OPEN_ARCHIVES supersedes
C		N_OPEN_ARCHIVES.
C	R. Kummerer, December 23, 1987. Dump FPR science catalog and matrix.
C	F. Shuman,   March    10, 1988. Coadd split overhaul.
C	R. Kummerer, March    30, 1988. Dump FPR coadd catalog and matrix,
C					short science archives.
C	F. Shuman,   August   23, 1988. Dump new (Time-Ordered) Ph Cor'd Spec
C					and Uncomb'd Cal'd Spec.
C	J. Bonnell,  September14, 1988. dump new archived house keeping
C					_eng and _idx files
C	F. Shuman,   October  20, 1988. Dump LMAC temperature fields--
C					Subroutines _List_ENG and _List_HKP.
C	F. Shuman,   January   5, 1989. SPR 2552.  Added error handler FGA_ERROR
C					to trap err msgs to screen.
C	R. Kummerer, April    27, 1989. SER 2956. Dump remaining archives.
C	R. Kummerer, May       4, 1989. SPR 3622. Appropriate dump file names.
C	R. Kummerer, July     22, 1989. SER 4173. Dump new archives FPP_SDF_xx
C					and FDQ_SDF_xx where xx is the channel
C					id (RH, RL, LH, LL).
C	Q. Chung,    August   28, 1989. SER 3306. Provide version number to
C					track software update.
C       Q. Chung,    October      1989. SPR 4750. Formatted DUMP FEX_LIMFLAGS,
C                                       FEX_SCILIM,FEX_ENGLIM FILES,
C                    November  8, 1989. FEX_IDX_FLAG,FEX_IDX_TOLS,
C                    November  9, 1989. FEX_CTH,
C                    November 16, 1989. FEX_GRTSWT.
C	R. Kummerer, January  30, 1990. SPR 4530. Gracefully exit FGA if no
C					dataset is selected for dumping by
C					the query builder.
C	R. Kummerer, May      25, 1990. SPR 6158. Remove inappropriate
C					questions prompted prior to dumping
C					reference datasets.
C	R. Kummerer, May      25, 1990. SPR 6159. Version 6.1. Clean up
C					reference data set lists.
C       N. Gonzales, August    8, 1990. SPR 7245. If CT_KEYREAD_ARCV return
C                                       is not CTP_NORMAL, signal a fatal error.
C       N. Gonzales, September 4, 1991. SPR 8372, Updated the code to handle
C                                       new and revised rdl's.
C	S. Alexander, October 15, 1991. SPR 9163, Add translation of archive
C					logical name CSDR$FIRAS_ARCHIVE to
C					output dump file.
C       N. Gonzales, December 12, 1991. SPR 9281, Replaced calls to CT_Open_ARCV
C                               with FORTRAN open statements using useropen and
C                               also Add ref. file FEX_Mincoadd, SPR 9099.
C       N. Gonzales, June 29, 1992. SPR 9795. Remove obsolete FIRAS datasets
C                               from FCI, FES, FCS, FFC, FPR, FPS and FSF.
C	S. Alexander, January 21, 1993. SPR 10487. Remove references to
C				FGA_LIST_CTH and FGA_LIST_NOS.
C	S. Alexander, May 2, 1994. SPR 11741. Remove references to 
C                                             FGA_LIST_GRTSWT.
C-------------------------------------------------------------------------

	IMPLICIT	NONE

	INCLUDE		'CT$LIBRARY:CTUSER.INC'
	INCLUDE		'(CCT_STATUS_RECORD)'
	INCLUDE         '(UPM_STAT_MSG)'
	INCLUDE		'(FUT_PARAMS)'
	INCLUDE		'(FGA_ERROR)'
	INCLUDE		'($SSDEF)'

	LOGICAL		*1	NORMAL
	CHARACTER	*80	COMMENT
	CHARACTER	*16	FNAME
	CHARACTER	*80	LINE
        CHARACTER       *64     ENGLIM_FILE
        CHARACTER       *64     LIMFLAGS_FILE
        CHARACTER       *64     SCILIM_FILE
        CHARACTER       *64     IDXTOLS_FILE
        CHARACTER       *64     IDXFLAG_FILE
	CHARACTER       *64     CTH_FILE
	CHARACTER       *64     FACR_FILE
	CHARACTER       *64     FAKER_FILE
	CHARACTER       *64     FAKEL_FILE
	CHARACTER       *64     GAINR_FILE
	CHARACTER       *64     GAINL_FILE
	CHARACTER       *64     GTRAN_FILE
	CHARACTER       *64     GRAWT_FILE
	CHARACTER       *64     MTMSWEEP_FILE
	CHARACTER       *64     VABSAA_FILE
	CHARACTER       *64     MINCOADD_FILE
        CHARACTER       *21     ARCH_ID/'CSDR$FIRAS_REFERENCE:'/
	CHARACTER       *72	ARCH_IN
        CHARACTER       *10     SCILIM_DATASET
        CHARACTER       *10     ENGLIM_DATASET
        CHARACTER       *12     LIMFLAGS_DATASET
	CHARACTER       *12     IDXFLAG_DATASET
        CHARACTER       *12     IDXTOLS_DATASET
        CHARACTER       *7      CTH_DATASET
	CHARACTER       *12     FACR_DATASET
	CHARACTER       *12     FAKER_DATASET
	CHARACTER       *12     FAKEL_DATASET
	CHARACTER       *12     GAINR_DATASET
	CHARACTER       *12     GAINL_DATASET
	CHARACTER       *12     GTRAN_DATASET
	CHARACTER       *12     GRAWT_DATASET
	CHARACTER       *12     MTMSWEEP_DATASET
	CHARACTER       *12     VABSAA_DATASET
	CHARACTER       *12     MINCOADD_DATASET

	CHARACTER	*14	STARTING_TIME
	CHARACTER	*14	ENDING_TIME
	CHARACTER       *14     KEY_TIME
	CHARACTER	*30	TIME_RANGE
	INTEGER		*2	DATASET_ID
	CHARACTER       *12     DATASET_NAME
	INTEGER         *2      NAMLEN
        INTEGER         *2      DATABASE_ID
	INTEGER         *4	PSTATUS		      ! Status of FGA processing
	INTEGER         *4	SUCCESS/1/, ERROR/2/  ! Local status values
	INTEGER	        *4      INTERACTIVE           ! Oper. mode indicator
	LOGICAL         *1	PROMPT_TIMES
	INTEGER         *4	LSTART
	INTEGER         *4	LSTOP

	CHARACTER	*128	RSE(16)
	LOGICAL		*1	OLDRSE/.FALSE./
        LOGICAL         *1      EOF
	CHARACTER       *64     RSE_FILE       !FILENAME FOR RSE FILE
	INTEGER	        *4      RSE_PRESENT
	INTEGER         *4      GMT_JSTART(2)
	INTEGER         *4      GMT_JSTOP(2)
	INTEGER         *4      FUT_QUERY_RSE
	INTEGER         *4      FUT_GET_RSE
	INTEGER         *4      NEG1 / 'ffffffff'x /
	INTEGER         *4	SFCN
	INTEGER         *4	WINDOW
	INTEGER         *4	I	       ! A COUNTER

	CHARACTER	*8	ANS
	CHARACTER       *64     FILENAME
	INTEGER         *4      LEN
	INTEGER         *4      LENGTH
	INTEGER		*2	DATATYP
        INTEGER         *4      IOS

	INTEGER		*4	RSTATUS
	CHARACTER	*72	LOGN, FLOGN, TLOGN
	INTEGER		*4	FLEN, TLEN, LARC
	
	BYTE			WRITE_SCI
	BYTE			WRITE_IFG
	BYTE			WRITE_NUM
	BYTE			WRITE_WTS
	BYTE			WRITE_SPC
	BYTE			WRITE_SIG
	BYTE			WRITE_MNS
	BYTE			WRITE_MOD
	BYTE			WRITE_VAR
	BYTE			WRITE_CHI
	BYTE			WRITE_COV
	BYTE			WRITE_RES
	BYTE			WRITE_ANC
	BYTE			WRITE_APO

	INTEGER		*2	CT_LUN
	INTEGER		*2	SCI_LUN
	INTEGER		*2	SCI_ID
	INTEGER		*2	CT_STAT(20)
	INTEGER		*4	CAT_ENTRY
	INTEGER		*4	RESTAT
	INTEGER		*4	CHAN
	INTEGER		*2	NL
	INTEGER		*2	ISTAT
	INTEGER         *4      STATUS
	INTEGER         *4      IO_STAT
	INTEGER         *4      IO_NORMAL
	INTEGER		*2	NMJFR
	INTEGER		*2	NREC/0/
	INTEGER		*4	WRITE_LUN

	INTEGER		*4	LUN_OUT/6/
	INTEGER		*4	NUM_80/80/
	INTEGER		*4	NUM_132/132/
	CHARACTER	*6	VERSION
	PARAMETER   	(VERSION='11.9')

        INTEGER         *4      CT_CONNECT_READ
	INTEGER         *4      CT_CONNECT_QUERY_NONRAW
	INTEGER         *4      CT_GMT_TO_BINARY
	INTEGER         *4      UPM_GET_VALUE
	INTEGER         *4      UPM_PRESENT
	INTEGER         *4      STR$TRIM
	INTEGER		*4  	CUT_REGISTER_VERSION
	INTEGER		*4	CUT_DISPLAY_BANNER
	INTEGER		*4 	CUT_TRANSLATE_ARCHIVE_ID
	INTEGER		*4	CUT_READ_DAFS_RECORD_IGSE_KEY
	INTEGER         *4      CUT_READ_DAFS_RECORD_CSDR_KEY
	INTEGER		*4	FUT_QGET_ENG
	INTEGER		*4	FUT_QGET_HKP
	INTEGER		*4	FUT_QGET_HSK
	INTEGER		*4	FUT_QGET_IDX
	INTEGER		*4	FUT_QGET_SCI
	INTEGER		*4	FUT_QGET_EMF
	INTEGER		*4	FUT_QGET_NOS
	INTEGER		*4	FUT_QGET_SSC
	INTEGER		*4	FUT_QGET_EXT
	INTEGER		*4	FUT_QGET_ANC
        INTEGER         *4      FUT_QGET_ENGLIM
        INTEGER         *4      FUT_QGET_LIMFLAGS
        INTEGER         *4      FUT_QGET_SCILIM
        INTEGER         *4      FUT_QGET_IDXTOLS
        INTEGER         *4      FUT_QGET_IDXFLAGS
	INTEGER         *4      FUT_QGET_CTH
	INTEGER         *4      FUT_QGET_FACR
	INTEGER         *4      FUT_QGET_FAKER
	INTEGER         *4      FUT_QGET_FAKEL
	INTEGER         *4      FUT_QGET_GAINR
	INTEGER         *4      FUT_QGET_GAINL
	INTEGER         *4      FUT_QGET_GTRAN
	INTEGER         *4      FUT_QGET_GRAWT
	INTEGER         *4      FUT_QGET_MTMSWEEP
	INTEGER         *4      FUT_QGET_VABSAA
	INTEGER         *4      FUT_QGET_MINCOADD

	INTEGER         *4	FGA_LIST_ENGLIM
	INTEGER         *4	FGA_LIST_LIMFLAGS
	INTEGER         *4	FGA_LIST_SCILIM
	INTEGER         *4	FGA_LIST_IDXTOLS
	INTEGER         *4	FGA_LIST_IDXFLAGS
	INTEGER         *4      FGA_LIST_FACR
	INTEGER         *4      FGA_LIST_FAKER
	INTEGER         *4      FGA_LIST_FAKEL
	INTEGER         *4      FGA_LIST_GAINR
	INTEGER         *4      FGA_LIST_GAINL
	INTEGER         *4      FGA_LIST_GTRAN
	INTEGER         *4      FGA_LIST_GRAWT
	INTEGER         *4      FGA_LIST_MTMSWEEP
	INTEGER         *4      FGA_LIST_VABSAA
	INTEGER         *4      FGA_LIST_MINCOADD

	INTEGER		*4	LIB$GET_LUN
	INTEGER		*4	LIB$FREE_LUN
	INTEGER		*4	STR$UPCASE

	RECORD /CCT_STATUS/ STRUC

	DICTIONARY 'CDU_DAFS'
	RECORD /CDU_DAFS/ DAFS_REC

	DICTIONARY 'NFS_HKP'
	RECORD /NFS_HKP/ HKP_REC
	DICTIONARY 'NFS_SDF'
	RECORD /NFS_SDF/ SCI_REC
	DICTIONARY 'NFS_EMF'
	RECORD /NFS_EMF/ EMF_REC
	DICTIONARY 'NFS_ANC'
	RECORD /NFS_ANC/ ANC_REC

	DICTIONARY 'FDQ_IDX'
	RECORD /FDQ_IDX/ IDX_REC
	DICTIONARY 'FDQ_ENG'
	RECORD /FDQ_ENG/ ENG_REC
	DICTIONARY 'FDQ_ETR'
	RECORD /FDQ_ETR/ ETR_REC
	DICTIONARY 'FXT_ENG_XTRM'
	RECORD /FXT_ENG_XTRM/ EXT_REC

	DICTIONARY 'FEC_SSCAL'
	RECORD /FEC_SSCAL/ SSC_REC

	DICTIONARY 'FNT_NOISE'
	RECORD /FNT_NOISE/ NOS_REC

        DICTIONARY 'FEX_ENGLIM'
        RECORD /FEX_ENGLIM/ ENGLIM_REC
        DICTIONARY 'FEX_LIMFLAGS'
        RECORD /FEX_LIMFLAGS/ LIMFLAGS
        DICTIONARY 'FEX_SCILIM'
        RECORD /FEX_SCILIM/ SCILIM
        DICTIONARY 'FEX_IDX_TOLS'
        RECORD /FEX_IDX_TOLS/ IDXTOLS
        DICTIONARY 'FEX_IDX_FLAG'
        RECORD /FEX_IDX_FLAG/ IDXFLAG
        DICTIONARY 'FEX_CTH'
        RECORD /FEX_CTH/ CTH_REC

        DICTIONARY 'FEX_AV_CALRS'
        RECORD /FEX_AV_CALRS/ FACR_REC
        DICTIONARY 'FEX_FAKEIT'
        RECORD /FEX_FAKEIT/ FAKER_REC, FAKEL_REC
        DICTIONARY 'FEX_GAIN'
        RECORD /FEX_GAIN/ GAINR_REC, GAINL_REC
        DICTIONARY 'FEX_GRTTRANS'
        RECORD /FEX_GRTTRANS/ GTRAN_REC
        DICTIONARY 'FEX_GRTRAWWT'
        RECORD /FEX_GRTRAWWT/ GRAWT_REC
        DICTIONARY 'FEX_MTMSWEEP'
        RECORD /FEX_MTMSWEEP/ MTMSWEEP_REC
        DICTIONARY 'FEX_VABSAA'
        RECORD /FEX_VABSAA/ VABSAA_REC
        DICTIONARY 'FEX_MINCOADD'
        RECORD /FEX_MINCOADD/ MINCOADD_REC

	EXTERNAL        CT_CONNECT_READ
	EXTERNAL        CT_CONNECT_QUERY_NONRAW
	EXTERNAL  	CUT_TRANSLATE_ARCHIVE_ID
        EXTERNAL	CUT_NORMAL
	EXTERNAL	FUT_NORMAL
	EXTERNAL	FUT_EOF
	EXTERNAL	FGA_NORMAL
	EXTERNAL        FGA_DSNAMLEN
	EXTERNAL        FGA_MISSDSNAME
	EXTERNAL	FGA_NODATA
	EXTERNAL	FGA_OPEN_ERR
	EXTERNAL	FGA_WRITE_ERR
	EXTERNAL	FGA_CLOSE_ERR
	EXTERNAL	FGA_CTINIT_ERR
	EXTERNAL	FGA_CTOPEN_ERR
	EXTERNAL	FGA_CTQBLD_ERR
	EXTERNAL	FGA_CTCLOSE_ERR
	EXTERNAL	FGA_CTRARC_ERR
	EXTERNAL	FGA_ERROR

5	CALL LIB$ERASE_PAGE (1, 1)

c
c Set status for FGA processing to success.
c
	PSTATUS = SUCCESS
c
c Establish condition handler and display program banner.

	CALL LIB$ESTABLISH (FGA_ERROR)

	RESTAT = CUT_REGISTER_VERSION(VERSION)
	RESTAT = CUT_DISPLAY_BANNER(LUN_OUT,NUM_80,
	1			'FIRAS Facility FGA_Get_Archive')
c
c Get the operational mode. the default is 'interactive'.
c
	INTERACTIVE = FAC_PRESENT
	IF ((UPM_PRESENT('INTERACTIVE') .EQ. UPM_Pres) .OR.
	1   (UPM_PRESENT('INTERACTIVE') .EQ. UPM_Defaulted)) THEN
	  INTERACTIVE = FAC_PRESENT
	ELSEIF (UPM_PRESENT('INTERACTIVE') .EQ. UPM_Negated) THEN
	  INTERACTIVE = FAC_NOT_PRESENT
	ENDIF
c
c If batch mode, get the DAFS dataset name.
c
	IF (UPM_PRESENT('DATASET')) THEN
	  RESTAT = UPM_GET_VALUE('DATASET', DATASET_NAME, NAMLEN)
	  IF (NAMLEN .GT. 12) THEN
	    PSTATUS = ERROR
	    CALL LIB$SIGNAL(FGA_DSNAMLEN,%VAL(1),DATASET_NAME(1:12))
	  ENDIF
	ELSE IF (UPM_PRESENT('DATASET') .NE. UPM_Absent) Then
	  PSTATUS = ERROR
	  CALL LIB$SIGNAL(FGA_MISSDSNAME)
	ENDIF
c
c Get the JSTART and JSTOP times if running from an automation command file:
c
	STARTING_TIME = FAC_JSTART_DEFAULT
	ENDING_TIME = FAC_JSTOP_DEFAULT
	PROMPT_TIMES = .TRUE.

	IF (UPM_PRESENT ( 'JSTART' ) .EQ. UPM_Pres) THEN
	  RESTAT = UPM_GET_VALUE ( 'JSTART', STARTING_TIME, NAMLEN )

	  IF (RESTAT .EQ. SS$_NORMAL) THEN
	    LSTART = INDEX(STARTING_TIME,' ') - 1
	    IF (LSTART .EQ. -1) LSTART=14
	    STARTING_TIME = STARTING_TIME(1:LSTART) //
	2                   FAC_JSTART_DEFAULT(LSTART+1:)
	  END IF

	  PROMPT_TIMES = .FALSE.
	END IF

	IF (UPM_PRESENT ( 'JSTOP' ) .EQ. UPM_Pres) THEN
	  RESTAT = UPM_GET_VALUE ( 'JSTOP', ENDING_TIME, NAMLEN )

	  IF (RESTAT .EQ. SS$_NORMAL) THEN
	    LSTOP = INDEX(ENDING_TIME,' ') - 1
	    IF (LSTOP .EQ. -1)LSTOP=14
	    ENDING_TIME = ENDING_TIME(1:LSTOP) //
	2                 FAC_JSTOP_DEFAULT(LSTOP+1:)
	  END IF

	  PROMPT_TIMES = .FALSE.
	END IF

	TIME_RANGE = STARTING_TIME//';'//ENDING_TIME//';'
	NORMAL = .TRUE.
c
c Initialize:
c
	IF (PSTATUS .EQ. SUCCESS) THEN
	  CALL CT_INIT(STRUC)
c
c If batch mode, get the name of the Rse file.
c
	  RSE_Present = fac_not_present
	  DO I = 1,16
	     RSE(I) = ' '
	  ENDDO
	  DO I = 1,16
	     CALL LIB$MOVC3(4,NEG1,%REF(RSE(I)(125:128)))
	  ENDDO

	  IF (UPM_PRESENT('RSE') .EQ. UPM_Pres) THEN
	     RSE_Present = FAC_PRESENT
	     RESTAT = UPM_GET_VALUE('RSE', RSE_FILE, NAMLEN)
             IF (RSE_FILE .EQ. ' ') THEN
               CALL CT_GMT_TO_BINARY(STARTING_TIME,GMT_JSTART)
               RESTAT = FUT_QUERY_RSE('FEX_FTBRSE', GMT_JSTART, GMT_JSTART, RSE)
	       IF (.NOT. RESTAT) THEN
	          PSTATUS = ERROR
	       ENDIF
             ELSE
                RESTAT = FUT_GET_RSE (RSE_FILE, RSE)
	        IF (.NOT. RESTAT) THEN
                   PSTATUS = ERROR
                ENDIF
	     END IF             ! (RSE_FILE .EQ. ' ')
	     SFCN = 1
	     WINDOW = 1
	  END IF
	END IF

	IF (((STRUC.CTERR .EQ. CTP_NORMAL) .OR. (PSTATUS .EQ. SUCCESS))
	1  .AND. (RESTAT)) THEN 
           IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	      CALL CT_QUERY_BLD(CTU_$FIRAS, TIME_RANGE, RSE, STRUC, OLDRSE)
	      IF (STRUC.CTERR .NE. CTP_NORMAL) THEN 
	         CALL LIB$SIGNAL (FGA_CTQBLD_ERR, %VAL(1), %VAL(STRUC.CTERR))
	         NORMAL = .FALSE.
	      END IF
	      DATASET_ID=STRUC.DATASET_ID
	      RESTAT=CUT_READ_DAFS_RECORD_IGSE_KEY(CTU_$FIRAS,DATASET_ID,
	1                                       DAFS_REC)
	   ELSE
	      RESTAT = CUT_READ_DAFS_RECORD_CSDR_KEY (DATASET_NAME,DAFS_REC)
	      IF ( .NOT. RESTAT ) THEN
	         PSTATUS = ERROR
	         CALL LIB$SIGNAL(%VAL(RESTAT))
	         NORMAL = .FALSE.
	      ENDIF
	      DATASET_ID = DAFS_REC.IGSE_ID.IGSE_DATASET_ID
	   END IF
	END IF
	IF (((STRUC.CTERR .EQ. CTP_NORMAL) .OR. (PSTATUS .EQ. SUCCESS))
	1  .AND. (RESTAT)) THEN 
	    NL = INDEX(DAFS_REC.DATASET,' ')
	    IF (NL .EQ. 0) THEN
	       NL = LEN(DAFS_REC.DATASET)
	    ELSE
	       NL = NL - 1
	    END IF

	    FNAME = DAFS_REC.DATASET(1:NL) // '.LIS'
	    DATATYP = -999

	    IF (DATASET_ID .EQ. CTU_$FIR_HKP) THEN
	       DATATYP = CTU_$FIR_HKP

               IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	          TYPE 6701
6701	          FORMAT (/' List both major frames? ([Y]/N) ', $)
	          ACCEPT 5001,ANS
5001	          FORMAT (A)
	          RESTAT = STR$UPCASE(ANS,ANS)
	          IF (ANS(1:1) .EQ. 'N') THEN
	             TYPE 6702
6702	             FORMAT (/' Only the 1st major frame in each record will',
	2                 ' be listed'/)
	             NMJFR = 1
	          ELSE
	             NMJFR = 2
	          END IF
	       ELSE
                   NMJFR = 1
	       END IF

	    ELSE IF ((DATASET_ID .GE. CTU_$FIR_SS1) .AND.
	2            (DATASET_ID .LE. CTU_$FIR_SS4)) THEN
	       DATATYP = CTU_$FIR_SS1
               IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	          TYPE 8830
8830	          FORMAT (' Dump the full science record? (Y/[N]) ', $)
	          ACCEPT 5001,ANS
	          RESTAT = STR$UPCASE(ANS,ANS)
	          IF (ANS(1:1) .EQ. 'Y') THEN
	             WRITE_SCI = .TRUE.
	          ELSE
	             WRITE_SCI = .FALSE.
	          END IF
	       ELSE
	           WRITE_SCI = .FALSE.
	       END IF

	       IF (WRITE_SCI) THEN
	         IF (DATASET_ID .EQ. CTU_$FIR_SS1) THEN
		   SCI_ID = CTU_$FIR_RS1
	         ELSE IF (DATASET_ID .EQ. CTU_$FIR_SS2) THEN
		   SCI_ID = CTU_$FIR_RS2
	         ELSE IF (DATASET_ID .EQ. CTU_$FIR_SS3) THEN
		   SCI_ID = CTU_$FIR_RS3
	         ELSE IF (DATASET_ID .EQ. CTU_$FIR_SS4) THEN
		   SCI_ID = CTU_$FIR_RS4
	         END IF

                 IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	            TYPE 6601
6601	            FORMAT (' Dump IFG? (Y/[N]) ', $)
	            ACCEPT 5001,ANS
	            RESTAT = STR$UPCASE(ANS,ANS)
	            IF (ANS(1:1) .EQ. 'Y') THEN
	               WRITE_IFG = .TRUE.
	            ELSE
	               WRITE_IFG = .FALSE.
	            END IF
	         ELSE
	             WRITE_IFG = .FALSE.
	         END IF
	       END IF

	     ELSE IF ((DATASET_ID .GE. CTU_$FIR_RS1 .AND.
	2              DATASET_ID .LE. CTU_$FIR_RS4) .OR.
	3             (DATASET_ID .GE. FAC_FPP_RH .AND.
	4              DATASET_ID .LE. FAC_FPP_LL) .OR.
	5             (DATASET_ID .GE. FAC_FDQ_RH .AND.
	6              DATASET_ID .LE. FAC_FDQ_LL)) THEN
	       DATATYP = CTU_$FIR_RS1
               IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	          TYPE 6601
	          ACCEPT 5001,ANS
	          RESTAT = STR$UPCASE(ANS,ANS)
	          IF (ANS(1:1) .EQ. 'Y') THEN
	             WRITE_IFG = .TRUE.
	          ELSE
	             WRITE_IFG = .FALSE.
	          END IF
	       ELSE
	          WRITE_IFG = .FALSE.
	       END IF

	     ELSE IF (DATASET_ID .GE. CTU_$FIR_ED1 .AND.
	2             DATASET_ID .LE. CTU_$FIR_ED4) THEN
	       DATATYP = CTU_$FIR_ED1
               IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	          TYPE 6610
6610	          FORMAT (' Dump engineering dump data? (Y/[N]) ', $)
	          ACCEPT 5001,ANS
	          RESTAT = STR$UPCASE(ANS,ANS)
	          IF (ANS(1:1) .EQ. 'Y') THEN
	             WRITE_IFG = .TRUE.
	          ELSE
	             WRITE_IFG = .FALSE.
	          END IF
	        ELSE
	            WRITE_IFG = .FALSE.
	        END IF

	     ELSE IF (DATASET_ID .EQ. CTU_$FIR_EDF .or.
	1	      DATASET_ID .EQ. FAC_HKP_EDF) THEN
	       DATATYP = CTU_$FIR_EDF

	     ELSE IF (DATASET_ID .EQ. CTU_$FIR_ETR) THEN
	       DATATYP = CTU_$FIR_ETR

	     ELSE IF (DATASET_ID .EQ. FAC_ENG_XTRM) THEN
	       DATATYP = FAC_ENG_XTRM

	     ELSE IF (DATASET_ID .EQ. CTU_$FIR_IDX .OR.
	1	      DATASET_ID .EQ. FAC_HKP_IDX) THEN
	       DATATYP = CTU_$FIR_IDX

	     ELSE IF (DATASET_ID .GE. CTU_$FIR_NRH .AND.
	2             DATASET_ID .LE. CTU_$FIR_NLL) THEN
	       DATATYP = CTU_$FIR_NRH
               IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	          TYPE 6602
6602   	          FORMAT (' Dump Spectrum? (Y/[N]) ',$)	
	          ACCEPT 5001,ANS
 	          RESTAT = STR$UPCASE(ANS,ANS)
	          IF (ANS(1:1) .EQ. 'Y') THEN
	             WRITE_SPC = .TRUE.
	          ELSE
	             WRITE_SPC = .FALSE.
	          END IF

	          TYPE 8810
8810   	          FORMAT (' Dump Spectrum Sigmas? (Y/[N]) ',$)	
	          ACCEPT 5001,ANS
	          RESTAT = STR$UPCASE(ANS,ANS)
 	          IF (ANS(1:1) .EQ. 'Y') THEN
	             WRITE_SIG = .TRUE.
 	          ELSE
	             WRITE_SIG = .FALSE.
	          END IF
 	       ELSE
	           WRITE_SPC = .FALSE.
	           WRITE_SIG = .FALSE.
	       END IF

	     ELSE IF (DATASET_ID .EQ. FAC_FIR_ANC) THEN
	       DATATYP = FAC_FIR_ANC
               IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
   	          TYPE 6603
6603	          FORMAT (' Dump Ancillary Housekeeping Data? (Y/[N]) ', $)
	          ACCEPT 5001,ANS
	          RESTAT = STR$UPCASE(ANS,ANS)
	          IF (ANS(1:1) .EQ. 'Y') THEN
	             WRITE_ANC = .TRUE.
	          ELSE
 	             WRITE_ANC = .FALSE.
 	          END IF
	       ELSE
 	          WRITE_ANC = .FALSE.
 	       END IF

	     ELSE IF (DATASET_ID .EQ. FAC_FIR_ENGLIM) THEN
	       DATATYP = FAC_FIR_ENGLIM

	     ELSE IF (DATASET_ID .EQ. FAC_FIR_LIMFLAG) THEN
	       DATATYP = FAC_FIR_LIMFLAG

	     ELSE IF (DATASET_ID .EQ. FAC_FIR_SCILIM) THEN
	       DATATYP = FAC_FIR_SCILIM

	     ELSE IF (DATASET_ID .EQ. FAC_FIR_IDX_TOLS) THEN
	       DATATYP = FAC_FIR_IDX_TOLS

	     ELSE IF (DATASET_ID .EQ. FAC_FIR_IDX_FLAG) THEN
	       DATATYP = FAC_FIR_IDX_FLAG

	     ELSE IF (DATASET_ID .EQ. FAC_FIR_CTH) THEN
	       DATATYP = FAC_FIR_CTH

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_FACR) THEN
	       DATATYP = FAC_FRD_FACR

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_FAKER) THEN
	       DATATYP = FAC_FRD_FAKER

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_FAKEL) THEN
	       DATATYP = FAC_FRD_FAKEL

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_GAINR) THEN
	       DATATYP = FAC_FRD_GAINR

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_GAINL) THEN
	       DATATYP = FAC_FRD_GAINL

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_GTRAN) THEN
	       DATATYP = FAC_FRD_GTRAN

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_GRAWT) THEN
	       DATATYP = FAC_FRD_GRAWT

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_MTMSWEEP) THEN
	       DATATYP = FAC_FRD_MTMSWEEP

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_VABSAA) THEN
	       DATATYP = FAC_FRD_VABSAA

	     ELSE IF (DATASET_ID .EQ. FAC_FRD_MINCOADD) THEN
	       DATATYP = FAC_FRD_MINCOADD

	     ELSE
	       TYPE 6002
	       GOTO 5
	     END IF

6002	     FORMAT (//' Invalid response, please try again!')


C
C Dump the records in the selected timerange.
C
	     RESTAT = LIB$GET_LUN (WRITE_LUN)

	     IF (((STRUC.CTERR .EQ. CTP_NORMAL) .OR. (PSTATUS .EQ. SUCCESS))
	1      .AND. (RESTAT)) THEN 

	       OPEN (UNIT=WRITE_LUN, FILE=FNAME, TYPE='NEW', IOSTAT=ISTAT)

	       IF (ISTAT .EQ. 0) THEN

		 RESTAT = CUT_REGISTER_VERSION(VERSION)
		 RESTAT = CUT_DISPLAY_BANNER(WRITE_LUN,NUM_132,
	1				'FIRAS Facility FGA_Get_Archive')
		 WRITE (WRITE_LUN,61)

C  Write translation of logical names.

                 ARCH_IN = 'CSDR$FIRAS_ARCHIVE'
	         CALL STR$UPCASE(ARCH_IN,ARCH_IN)
	         LARC = LEN(ARCH_IN)
	         LOGN(1:LARC) = ARCH_IN
	         RSTATUS = 
     .              CUT_TRANSLATE_ARCHIVE_ID(LOGN, FLOGN, FLEN, TLOGN, TLEN)
	         WRITE (WRITE_LUN,40) ARCH_IN, TLOGN
40  	FORMAT(1X,'Logical Translation for Archive:',4X,A,' = ',/,25X,A)

		 WRITE (WRITE_LUN,61)
61		 FORMAT (//)
  
                 IF (INTERACTIVE .EQ. FAC_PRESENT) THEN
	            TYPE 6003
6003	            FORMAT (/' Enter a comment line for the list '
	2                 'file (max 80 characters): ')
	            ACCEPT 5001,LINE
	         END IF
               END IF
C
C   Open archive(s)
C
	      RESTAT = STR$TRIM(DAFS_REC.DATASET, DAFS_REC.DATASET, LENGTH)
	      FILENAME = 'CSDR$FIRAS_ARCHIVE:' // DAFS_REC.DATASET(1:LENGTH)
	2                // '/' // TIME_RANGE

	      RESTAT = LIB$Get_LUN(CT_LUN)
	      OPEN (UNIT = CT_LUN,
	2          FILE = FILENAME,
	3          STATUS = 'OLD',
	4          IOSTAT = IO_STAT,
	5          USEROPEN = CT_CONNECT_QUERY_NONRAW)

	         IF (DATATYP .EQ. CTU_$FIR_SS1 .AND. WRITE_SCI .AND.
	2	     STRUC.CTERR .EQ. CTP_NORMAL) THEN
	             FILENAME = 'CSDR$FIRAS_RAW:' // 'FDQ_SDF_' //
	2            DAFS_REC.DATASET(LENGTH-1:LENGTH) // '/' // TIME_RANGE
	             RESTAT = LIB$Get_LUN(SCI_LUN)
	             OPEN (UNIT = SCI_LUN,
	2                 FILE = FILENAME,
	3                 STATUS = 'OLD',
	4                 IOSTAT = IO_STAT,
	5                 USEROPEN = CT_CONNECT_READ)
c
c	Failed to open archive.
c	If the io_stat is not normal, signal cobetrieve open error.
c
           	     IF (IO_STAT .NE. IO_NORMAL) THEN
                        CALL LIB$SIGNAL(FGA_CTOPEN_ERR,%VAL(1),%VAL(IO_STAT))
	             END IF
	         END IF

	         IF (((STRUC.CTERR .EQ. CTP_NORMAL) .OR. (PSTATUS .EQ. SUCCESS))
	1           .AND. (RESTAT)) THEN 

	           CALL CT_QUERY_ARCV(, CT_LUN, RSE, STRUC)

	           IF (((STRUC.CTERR .EQ. CTP_NORMAL) .OR. (PSTATUS .EQ. SUCCESS))
	1              .AND. (RESTAT)) THEN 

	             DO WHILE (NORMAL)
C
C Dump housekeeping.
C
	               IF (DATATYP .EQ. CTU_$FIR_HKP) THEN
	                 RESTAT = FUT_QGET_HKP (CT_LUN, HKP_REC)
	                 IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC =  NREC + 1
	                   CALL FGA_LIST_HKP (WRITE_LUN, HKP_REC, NREC, NMJFR,
	2                                     LINE, CAT_ENTRY)
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump raw science.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_RS1) THEN
	                 RESTAT = FUT_QGET_SCI (CT_LUN, SCI_REC)
	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_SCI (WRITE_LUN, SCI_REC, NREC, LINE,
	2                                     CAT_ENTRY, WRITE_IFG)
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump short science.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_SS1) THEN
	                 RESTAT= FUT_QGET_SSC (CT_LUN, SSC_REC)
	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_SSC (WRITE_LUN, SSC_REC, NREC, LINE,
	2                                     CAT_ENTRY)
			   IF (WRITE_SCI) THEN
	                     CALL CT_KEYEDREAD_ARCV (, SCI_LUN, SCI_REC,
	2					     SSC_REC.TIME(1),
	3					     STRUC.CTERR)
	                     IF (STRUC.CTERR .EQ. CTP_NORMAL) THEN
	                       CALL FGA_LIST_SCI (WRITE_LUN, SCI_REC, NREC,
	2                                         LINE,	CAT_ENTRY, WRITE_IFG)
	                     ELSE
                                 CALL CT_BINARY_TO_GMT(SSC_REC.TIME, KEY_TIME)

	                     CALL LIB$SIGNAL (FGA_CTRARC_ERR, %VAL(1),key_time)

	                     END IF
	  		   END IF
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump engineering.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_EDF) THEN
	                 RESTAT= FUT_QGET_ENG (CT_LUN, ENG_REC)
	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_ENG (WRITE_LUN, ENG_REC, NREC, LINE,
	2                                     CAT_ENTRY)
	                 ELSE
	                   NORMAL = .FALSE.
    	                 END IF
C
C Dump engineering trends.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_ETR) THEN
	                 RESTAT= FUT_QGET_ENG (CT_LUN, ETR_REC)
	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_ENG (WRITE_LUN, ETR_REC, NREC, LINE,
	2                                     CAT_ENTRY)
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump engineeering mode.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_ED1) THEN
	                 RESTAT= FUT_QGET_EMF (CT_LUN, EMF_REC)
	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_EMF (WRITE_LUN, EMF_REC, NREC, LINE,
	2                                     CAT_ENTRY, WRITE_IFG)
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump engineering extrema.
C
	               ELSE IF (DATATYP .EQ. FAC_ENG_XTRM) THEN
	                 RESTAT= FUT_QGET_EXT (CT_LUN, EXT_REC)
	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_EXT (WRITE_LUN, EXT_REC, NREC, LINE,
	2                                     CAT_ENTRY)
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump engineering index.
C
	               ELSE IF (DATATYP .EQ. CTU_$FIR_IDX) THEN
	                 RESTAT= FUT_QGET_IDX (CT_LUN, IDX_REC)
 	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_IDX (WRITE_LUN, IDX_REC, NREC, LINE,
	2		 	              CAT_ENTRY)
	                 ELSE
	                   NORMAL = .FALSE.
	                 END IF
C
C Dump ancillary housekeeping data.
C
	               ELSE IF (DATATYP .EQ. FAC_FIR_ANC) THEN
	                 RESTAT= FUT_QGET_ANC (CT_LUN, ANC_REC)
	                 IF (RESTAT.EQ. %LOC(FUT_NORMAL)) THEN
	                   NREC = NREC + 1
	                   CALL FGA_LIST_ANC (WRITE_LUN, ANC_REC, NREC, LINE,
	2                                     CAT_ENTRY, WRITE_ANC)
	                 ELSE
	                   NORMAL = .FALSE.
 	                 END IF
C
C Dump engineering limits reference data.
C
                	ELSEIF (DATATYP .EQ. FAC_FIR_ENGLIM) THEN
                  	  ENGLIM_DATASET = DAFS_REC.DATASET(1:NL)
                  	  ENGLIM_FILE = ARCH_ID // ENGLIM_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_ENGLIM(CT_LUN,ENGLIM_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_ENGLIM (WRITE_LUN, ENGLIM_REC,
	1                                             ENGLIM_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump engineering limits flags reference data.
C
          	        ELSEIF (DATATYP .EQ. FAC_FIR_LIMFLAG) THEN
                  	  LIMFLAGS_DATASET = DAFS_REC.DATASET(1:NL)
                          LIMFLAGS_FILE = ARCH_ID // LIMFLAGS_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_LIMFLAGS(CT_LUN,LIMFLAGS)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_LIMFLAGS (WRITE_LUN, LIMFLAGS,
	1                                               LIMFLAGS_FILE, LINE,
	2						NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump science limits reference data.
C
                        ELSEIF (DATATYP .EQ. FAC_FIR_SCILIM) THEN
                          SCILIM_DATASET = DAFS_REC.DATASET(1:NL)
                          SCILIM_FILE = ARCH_ID // SCILIM_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_SCILIM(CT_LUN,SCILIM)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_SCILIM (WRITE_LUN, SCILIM,
	1                                             SCILIM_FILE, LINE, NREC)
                          ELSE
                    	    NORMAL = .FALSE.
                          END IF
C
C Dump engineering index tolerances reference data.
C
	                ELSEIF (DATATYP .EQ. FAC_FIR_IDX_TOLS) THEN
                          IDXTOLS_DATASET = DAFS_REC.DATASET(1:NL)
                          IDXTOLS_FILE = ARCH_ID // IDXTOLS_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_IDXTOLS(CT_LUN,IDXTOLS)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_IDXTOLS (WRITE_LUN, IDXTOLS,
	1                                              IDXTOLS_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump engineering index flags reference data.
C
               		  ELSEIF (DATATYP .EQ. FAC_FIR_IDX_FLAG) THEN
                  	    IDXFLAG_DATASET = DAFS_REC.DATASET(1:NL)
                  	    IDXFLAG_FILE = ARCH_ID // IDXFLAG_DATASET // '/' 
	1                      // TIME_RANGE
                            RESTAT = FUT_QGET_IDXFLAGS(CT_LUN,IDXFLAG)
		            IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                              NREC = NREC + 1
                              RESTAT = FGA_LIST_IDXFLAGS (WRITE_LUN, IDXFLAG,
	1                                                 IDXFLAG_FILE, LINE,
	2						  NREC)
                            ELSE
                              NORMAL = .FALSE.
                            END IF
C
C Dump average calibration resistances count.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_FACR) THEN
                          FACR_DATASET = DAFS_REC.DATASET(1:NL)
                          FACR_FILE = ARCH_ID // FACR_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_FACR(CT_LUN,FACR_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_FACR (WRITE_LUN, FACR_REC,
	1					   FACR_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump fex_fakeit_r reference dataset for fakeit statuses and gaps.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_FAKER) THEN
                          FAKER_DATASET = DAFS_REC.DATASET(1:NL)
                          FAKER_FILE = ARCH_ID // FAKER_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_FAKER(CT_LUN,FAKER_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_FAKER (WRITE_LUN, FAKER_REC,
	1					   FAKER_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump fex_fakeit_l reference dataset for fakeit statuses and gaps.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_FAKEL) THEN
                          FAKEL_DATASET = DAFS_REC.DATASET(1:NL)
                          FAKEL_FILE = ARCH_ID // FAKEL_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_FAKEL(CT_LUN,FAKEL_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_FAKEL (WRITE_LUN, FAKEL_REC,
	1					   FAKEL_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump fex_gain_r reference dataset for gain change and gaps.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_GAINR) THEN
                          GAINR_DATASET = DAFS_REC.DATASET(1:NL)
                          GAINR_FILE = ARCH_ID // GAINR_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_GAINR(CT_LUN,GAINR_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_GAINR (WRITE_LUN, GAINR_REC,
	1					   GAINR_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump fex_gain_l reference dataset for gain change and gaps.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_GAINL) THEN
                          GAINL_DATASET = DAFS_REC.DATASET(1:NL)
                          GAINL_FILE = ARCH_ID // GAINL_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_GAINL(CT_LUN,GAINL_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_GAINL (WRITE_LUN, GAINL_REC,
	1					   GAINL_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump GRT transition temperature and half-width.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_GTRAN) THEN
                          GTRAN_DATASET = DAFS_REC.DATASET(1:NL)
                          GTRAN_FILE = ARCH_ID // GTRAN_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_GTRAN(CT_LUN,GTRAN_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_GTRAN (WRITE_LUN, GTRAN_REC,
	1					   GTRAN_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump GRT weights for FDQ raw data.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_GRAWT) THEN
                          GRAWT_DATASET = DAFS_REC.DATASET(1:NL)
                          GRAWT_FILE = ARCH_ID // GRAWT_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_GRAWT(CT_LUN,GRAWT_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_GRAWT (WRITE_LUN, GRAWT_REC,
	1					   GRAWT_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump reference file FEX_MtmSweep.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_MTMSWEEP) THEN
                          MTMSWEEP_DATASET = DAFS_REC.DATASET(1:NL)
                          MTMSWEEP_FILE = ARCH_ID // MTMSWEEP_DATASET // '/' 
	1                      // TIME_RANGE
                          RESTAT = FUT_QGET_MTMSWEEP(CT_LUN,MTMSWEEP_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_MTMSWEEP (WRITE_LUN, 
	1		             MTMSWEEP_REC, MTMSWEEP_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump reference file FEX_VABSAA.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_VABSAA) THEN
                          VABSAA_DATASET = DAFS_REC.DATASET(1:NL)
                          VABSAA_FILE = ARCH_ID // VABSAA_DATASET // '/' 
	1                               // TIME_RANGE
                          RESTAT = FUT_QGET_VABSAA(CT_LUN, VABSAA_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_VABSAA (WRITE_LUN, 
	1		             VABSAA_REC, VABSAA_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF
C
C Dump reference file FEX_MINCOADD.
C
                        ELSEIF (DATATYP .EQ. FAC_FRD_MINCOADD) THEN
                          MINCOADD_DATASET = DAFS_REC.DATASET(1:NL)
                          MINCOADD_FILE = ARCH_ID // MINCOADD_DATASET // '/' 
	1                               // TIME_RANGE
                          RESTAT = FUT_QGET_MINCOADD(CT_LUN, MINCOADD_REC)
		          IF (RESTAT .EQ. %LOC(FUT_NORMAL)) THEN
                            NREC = NREC + 1
                            RESTAT = FGA_LIST_MINCOADD (WRITE_LUN, 
	1		             MINCOADD_REC, MINCOADD_FILE, LINE, NREC)
                          ELSE
                            NORMAL = .FALSE.
                          END IF

	                END IF

	             END DO

		   END IF

	            IF (RESTAT.EQ. %LOC(FUT_EOF) .AND. NREC .EQ. 0) THEN
	              CALL LIB$SIGNAL (FGA_NODATA)
	            END IF

	            IF (RESTAT.NE. %LOC(FUT_EOF)) THEN
	              CALL LIB$SIGNAL (RESTAT)
	            END IF

	            CALL CT_CLOSE_ARCV (, CT_LUN, STRUC)

	            IF (STRUC.CTERR .NE. CTP_NORMAL) THEN
	              CALL LIB$SIGNAL (FGA_CTCLOSE_ERR, %VAL(1),
	2				%VAL(STRUC.CTERR))
 	              NORMAL = .FALSE.
	            END IF

	          ELSE

	            CALL LIB$SIGNAL (FGA_CTOPEN_ERR, %VAL(1), %VAL(STRUC.CTERR))
	            NORMAL = .FALSE.

	          END IF

	          CLOSE (WRITE_LUN, IOSTAT=ISTAT)

	          IF (ISTAT .NE. 0) THEN
	            CALL LIB$SIGNAL (FGA_CLOSE_ERR, %VAL(1), %VAL(ISTAT))
	          END IF

	        ELSE

	          CALL LIB$SIGNAL (FGA_OPEN_ERR, %VAL(1), %VAL(ISTAT))

	        END IF

	        RESTAT = LIB$FREE_LUN(WRITE_LUN)

	        IF (RESTAT .NE. SS$_NORMAL) THEN
	          CALL LIB$SIGNAL (%VAL(RESTAT))
	        END IF

	   ELSE

	      IF (STRUC.DATASET_ID .NE. -1) THEN
	        CALL LIB$SIGNAL (%VAL(RESTAT))
	        NORMAL = .FALSE.
	      END IF

	   END IF

	IF (NREC .GT. 0) THEN
	 WRITE(6, 111) FNAME
111	 FORMAT (/' List file is named ', A  //)
	END IF

	IF (NORMAL) THEN
	 CALL LIB$SIGNAL(FGA_NORMAL)
	END IF

	END
