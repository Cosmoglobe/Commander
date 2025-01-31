	INTEGER*4 FUNCTION FDQ_OPEN_ARCV ( FILENAMES, SCI_OPEN, TIME_RANGE,
	1				   CT_LUN, ENGLIM,outnames )
C/
C/	PROGRAM NAME:
C/	  FDQ_OPEN_ARCV
C/
C/	PROGRAM DESCRIPTION:
C/	  The routine will open the necessary archive segments for processing
C/	  one set of science file segments.
C/	  This new routine replaces the old routine OPEN_ARCHIVES, which
C/	  was under the filename FDQ_OPEN_ARCV as well.
C/
C/	AUTHOR:
C/	  Edwin H. Fung
C/	  GSFC
C/	  May 6, 1987
C/
C/	MODIFIED BY:
C/	  Edwin H. Fung
C/	  GSFC
C/	  June 25, 1987
C/	  REASON:	Not to call CT_MODIFY_ARCV if the catalog entry of
C/			a Science file is -1 (no legal entry for that channel).
C/
C/      MODIFIED BY
C/         J. W. Durachta
C/         ARC
C/         July 24,1987
C/         REASON:       Reference to HSK replaced by HKP in accordance with
C/                       notch filter alterations.
C/
C/      MODIFIED BY
C/         J. W. Durachta
C/         ARC
C/         August 25,1987
C/         REASON:       Reference to RSI replaced by IDX in accordance with
C/                       changes to the index record structure.
C/
C/      MODIFIED BY
C/         J. W. Durachta
C/         ARC
C/         Sept. 3,1987
C/         REASON:       Open for engineering statistics (HES and DES) commented
C/                       out. These should be removed when convenient.
C/
C/      MODIFIED BY
C/         D. WARD
C/         GSFC/STX
C/         October 20, 1987
C/         REASON:      Replace HES and DES with ETR.
C/
C/      MODIFIED BY
C/	  Shirley M. Read
C/	  STX
C/	  January 6, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal. Removed return status from
C/			calling sequence since function value serves the same
C/			purpose.
C/
C/      MODIFIED BY
C/	  Shirley M. Read
C/	  STX
C/	  March, 1988
C/	  REASON: 	Converted to new COBETRIEVE FORTRAN user open using
C/	                filenames instead of catalog numbers. 
C/
C/      MODIFIED BY
C/	  Shirley M. Read
C/	  STX
C/	  May, 1988
C/	  REASON: 	Added comment Ct Unit 10 is reserved for the report
C/			file (not currently archived).
C/
C/      MODIFIED BY
C/	  R. Kummerer
C/	  STX
C/	  June, 1988
C/	  REASON: 	Add segment tracking file opens.
C/
ch
ch	version 4.1.1 12/01/88, ser 2379, J.T.Bonnell, STX@GSFC
ch		This program was modified to refer to the
ch		new firas archive logical names in response
ch		to SER 2379 (csdr$firas_in, _out, _raw,
ch		_ref, _uref, and _cal).
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH		There have been problems during test operations for FIRAS
CH		processing due to the required clean-up of the archives after
CH		an FPP or FDQ abort. The FPR tracking system compounds the
CH		problems. Files with non-matching version numbers seem often
CH		to result from improper clean-up. Bad record times cause
CH		SEGCTL to abort and mess up the tracking system. It was
CH		decided to change the modify of the science records in FPP
CH		and FDQ to a simple COBETRIEVE read of the existing records
CH		from a dataset and write a modifed dataset with the same
CH		information which was entered on the modify. Two new science
CH		data sets will be required: a science dataset of raw science
CH		data plus FPP input and a science dataset with FPP science
CH		data plus FDQ input. These datasets will be FPP_SDF_xx, where
CH		xx is the channel id (RH, RL, LH or LL) and FDQ_SDF_xx, where
CH		xx is the channel id. The new datasets must be opened and
CH		processed in FPP and FDQ. 
CH	Version 4.4.1 08/20/89, SER 4210, R. Kummerer, STX
CH		Prevent overlaps in raw science segments related changes.
CH	Version 4.4.3 11/14/89, SPR 5032, R. Kummerer STX
CH		Skip processing raw science segments from missing channels.
CH	Version 4.4.04 11/29/89, SPR 5207, R. Kummerer STX
CH		Skip closes on missing channel segments.
CH	Version 5.0 12/11/89, SPR 5313, R. Kummerer STX
CH		Perform appropriate raw science and IFG tracker archive closes.
CH      Version ??? 1/17/91, NFS_HKP only be need to open once
CH              H. Wang, STX.
CH  
CH      Modified by: H. Wang, STX, 1/29/91
CH      Reason: New requirements for FDQ
CH              
C/	CALLING SEQUENCE:
C/	  STATUS = FDQ_OPEN_ARCV (FILENAMES, SCI_OPEN, TIME_RANGE,
C/				  CT_LUN, ENGLIM )
C/
C/	INPUT PARAMETERS:
C/	  FILENAMES(5)  C*39    Filenames for 4 science and 1 HKP datasets.
C/	  TIME_RANGE	C*30	Timerange for housekeeping data.
C/	  ENGLIM	I*4	Engineering limits violation counts.
C/
C/	OUTPUT PARAMETERS:
C/	  SCI_OPEN(4)   L*1	Flags indicating which raw science archives
C/				have been opened.
C/	  CT_LUN(22)	I*2	Array of COBETRIEVE logical units for archival
C/				access:
C/					CT_LUN(1)  --	RS1 Input
C/					CT_LUN(2)  --	RS2
C/					CT_LUN(3)  --	RS3
C/					CT_LUN(4)  --	RS4
C/					CT_LUN(5)  --   HKP
C/					CT_LUN(6)  --   IDX
C/					CT_LUN(7)  --   ENG
C/					CT_LUN(8)  --   ETR
C/					CT_LUN(9)  --   CAL
C/					CT_LUN(10) --   RPT -- Report File
C/					CT_LUN(11)  --	SM1
C/					CT_LUN(12)  --	SM2
C/					CT_LUN(13)  --	SM3
C/					CT_LUN(14)  --	SM4
C/					CT_LUN(15)  --  SC1
C/					CT_LUN(16)  --  SC2
C/					CT_LUN(17)  --  SC3
C/					CT_LUN(18)  --  SC4
C/					CT_LUN(19)  --  ORS1 Output
C/					CT_LUN(20)  --  ORS2
C/					CT_LUN(21)  --  ORS3
C/					CT_LUN(22)  --  ORS4
C/
C/	INPUT/OUTPUT FILES:
C/	  FIRAS INDEX ARCHIVE
C/	  FIRAS HOUSEKEEPING ARCHIVE
C/	  FIRAS SCIENCE ARCHIVES
C/	  FIRAS ENGINEERING ARCHIVE
C/	  FIRAS ENGINEERING STATISTICS ARCHIVE (hourly and daily)
C/
C/	INCLUDE FILES USED:
C/	  CT$LIBRARY:CTUSER.INC
C/
C/	SUBROUTINES CALLED:
C/ 	  FUT_GET_LUN
C/ 	  LIB$GET_LUN (from system library)
C/	  LIB$SIGNAL (from system library)
C/	  LIB$LOCC (from system library)
C/
C/
C/	ERROR HANDLING:
C/	  Message to Lib$Signal and bad return in RETSTAT
C/
C/	METHOD USED:
C/	  The following is the PDL --
C/
C/	    If (first time around) then
C/	      Set FIRST to false;
C/	      Get logical unit numbers for files;
C/	      Open daily and hourly engineering statistic archives for write
C/		   using the file extension of the housekeeping file;
C/	    Else
C/	      Close all 4 science archives;
C/	      Close HKP archive;
C/	      Close ENG archive;
C/	      Close IDX (index) archive;
C/	    Endif;
C/
C/	    Do for all 4 channels
C/	      Use FILENAMES to open raw science for CT modify access;
C/	    Enddo;
C/
C/	    Use FILENAMES to open HKP archive for sequential read access;
C/	    Use FILENAMES to open FEX_AV_CALRS ref. archive for sequential read access;
C/
C/	    Open ENG archive using the same file extension as the housekeeping;
C/	    Open IDX (index) archive using the same file extension as the 
C/		 housekeeping;
C/
C/	    Return;
C/	    End.
C/

	IMPLICIT	NONE

!	Passed Parameters


	CHARACTER*39    FILENAMES(5)
	Character*39    OUTNAMES(7)
	LOGICAL*1	SCI_OPEN(4)
	LOGICAL*1	LAST_SCI_OPEN(4)
	CHARACTER*30	TIME_RANGE
	INTEGER*2	CT_LUN(22)
	INTEGER*4	LUN
	INTEGER*4	ENGLIM

!	Include Files and External Parameters

	INCLUDE 'CT$LIBRARY:CTUSER.INC'
	INCLUDE '(FUT_PARAMS)'
	INCLUDE '($SSDEF)'

	integer*4 	CT_Connect_Write
	external        CT_Connect_Write
	integer*4 	CT_Connect_Read
	external        CT_Connect_Read

	EXTERNAL   	FUT_NORMAL
	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR
	EXTERNAL	FDQ_LUNGETERR
	EXTERNAL	FDQ_LUNFREER
	EXTERNAL	FDQ_CTOPENERR
	EXTERNAL	FDQ_RMSOPENCTL
	EXTERNAL	FDQ_RMSREADCTL
	EXTERNAL	FDQ_RMSWRITECTL
	EXTERNAL	FDQ_RMSCLOSECTL
	EXTERNAL	FDQ_CTCLOSERR
	EXTERNAL	FDQ_NOFILEX
	EXTERNAL	FDQ_BADFILNAM

!	Functions

	integer*4 fut_get_lun    !Get unit number
	integer*4 lib$locc       !String position (from system library)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			  !
!     Local variables     !
!			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

	logical*1	first /.true./, found
	integer*2	chan, ct_stat(20)
	integer*4       ix		! index
	integer*4       ipos1           ! string position 
	integer*4 	retstat		! Return status
	integer*4 	success / 1 /, error / 2 /  ! Values for status
	integer*4       status	        ! System status
	integer*4       iostatus        ! FORTRAN I/O status
	integer*4       zero  / 0 /, ISDF

	character*15     arc_id_r /'CSDR$FIRAS_RAW:'/
	character*14     arc_id_in /'CSDR$FIRAS_IN:'/
	character*15     arc_id_out /'CSDR$FIRAS_OUT:'/
	character*15     arc_ref /'CSDR$FIRAS_REF:'/
	character*60     arc_file	! Storage for complete filename
	character*7      hkp_rdl /'NFS_HKP'/  ! NFS_HKP RDL name
	character*7	 etr_rdl /'FDQ_ETR'/  ! FDQ_ETR RDL name
	character*7	 eng_rdl /'FDQ_ENG'/  ! FDQ_ENG RDL name 
	character*7      idx_rdl /'FDQ_IDX'/  ! FDQ_IDX RDL name
	character*12      CAL_rdl /'FEX_AV_CALRS'/  ! FEX_AV_CALRS RDL name
	character*8      out_sci_rdl /'FDQ_SDF_'/  ! FDQ_SDF RDL name
	character*32     file_ext        ! Storage for file extension
	character*39     blank		 ! Blank string
	character*1      blanks(39) / 39 * ' ' / 
	equivalence      (blanks(1), blank)

!!!!!!!!!!!!!!!!!!!!!!
!		     !
!     parameters     !
!		     !
!!!!!!!!!!!!!!!!!!!!!!

	INTEGER*2  	GOOD, BAD, RS1, RS2, RS3, RS4, HKP, 
	1		IDX, ENG, ETR, CAL,RPT,
	2		SM1, SM2, SM3, SM4,
	3		SC1, SC2, SC3, SC4,
	4		ORS1, ORS2, ORS3, ORS4
	PARAMETER	(GOOD = 1)
	PARAMETER	(BAD = 2)

	PARAMETER	(RS1 = 1)
	PARAMETER	(RS2 = 2)
	PARAMETER	(RS3 = 3)
	PARAMETER	(RS4 = 4)
	PARAMETER	(HKP = 5)
	PARAMETER	(IDX = 6)
	PARAMETER	(ENG = 7)
	PARAMETER	(ETR = 8)
	PARAMETER	(CAL = 9)
	PARAMETER	(RPT = 10)
	PARAMETER	(SM1 = 11)
	PARAMETER	(SM2 = 12)
	PARAMETER	(SM3 = 13)
	PARAMETER	(SM4 = 14)
	PARAMETER	(SC1 = 15)
	PARAMETER	(SC2 = 16)
	PARAMETER	(SC3 = 17)
	PARAMETER	(SC4 = 18)
	PARAMETER	(ORS1 = 19)
	PARAMETER	(ORS2 = 20)
	PARAMETER	(ORS3 = 21)
	PARAMETER	(ORS4 = 22)

	Character*12	Datasets(22)
	Data Datasets / 'FPP_SDF_RH',
	1		'FPP_SDF_RL',
	1		'FPP_SDF_LH',
	1		'FPP_SDF_LL',
	1		'NFS_HKP',
	1		'FDQ_IDX',
	1		'FDQ_ENG',
	1		'FDQ_ETR',
	1		'FEX_AV_CALRS',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		' ',
	1		'FDQ_SDF_RH',
	1		'FDQ_SDF_RL',
	1		'FDQ_SDF_LH',
	1		'FDQ_SDF_LL' /

	SAVE LAST_SCI_OPEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Set return status to success.

	RETSTAT = SUCCESS

!	Extract the file extension name for the segment to build all
!	corresponding file names.

	  ipos1 = lib$locc('.',filenames(5))
	  if (ipos1 .eq. 0 ) then
	      retstat = error
	      call lib$signal(FDQ_NOFILEX, %val(1), filenames(5))
	  elseif (ipos1 .ne. 8 ) then
	      retstat = error
	      call lib$signal(FDQ_BADFILNAM, %val(1), filenames(5))
          endif
	  ipos1 = lib$locc('.',filenames(1))
           ISDF = 1
          do while ((ipos1 .ne. 11) .and. (isdf .lt. 4))
            isdf = isdf + 1  
            ipos1 = lib$locc('.',filenames(isdf))
          enddo
	  if (ipos1 .ne. 11 ) then
	      retstat = error
	      call lib$signal(FDQ_BADFILNAM, %val(1), filenames(5))
          endif
	  file_ext = filenames(isdf)(ipos1:39)

	IF (FIRST) THEN
	  FIRST = .FALSE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                        !
!	Get logical unit numbers for all datasets to be accessed by FDQ. !
!                                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  DO IX = 1, 9

	    if ( retstat .eq. success ) then 
	      status = fut_get_lun(lun)
	      ct_lun(ix)=lun
	      if ( status .ne. %Loc(FUT_Normal) ) then
	         retstat = error
	         call lib$signal (FDQ_LUNGETERR, %val(1), %val(status))
	      endif
	    endif		! retstat .eq.success

	  ENDDO

	  DO IX = 11, 22

	    if ( retstat .eq. success ) then 
	      status = fut_get_lun(lun)
	      ct_lun(ix)=lun
	      if ( status .ne. %Loc(FUT_Normal) ) then
	         retstat = error
	         call lib$signal (FDQ_LUNGETERR, %val(1), %val(status))
	      endif
	    endif		! retstat .eq.success

	  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					                        !
!     Open Engineering Statistics Archives for write-access     !
!					                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  If ( retstat .eq. success ) then

	    arc_file = arc_id_out // etr_rdl // file_ext
            outnames(6)=arc_file
	   
	    open ( unit=ct_lun(etr),
	1	   file=arc_file, shared,
	2	   status='NEW', iostat=iostatus,
	3	   useropen=CT_Connect_Write)
 
	    if (iostatus .ne. zero) then
	      retstat = error
	      call lib$signal(FDQ_CTOPENERR,%val(2),%val(iostatus),
	1	datasets(etr))
	    ENDIF
	  ENDIF  		! retstat is success
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								 !
!     Now open FIRAS HKP file for read.                          !
!								 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  If ( retstat .eq. success ) Then

	    arc_file = arc_id_r // hkp_rdl // '/' // TIME_RANGE
	    OPEN ( unit=ct_lun(hkp),
	1	   file=arc_file, 
	2	   status='OLD', iostat=iostatus,
	3	   useropen=CT_Connect_Read)
 
	    if (iostatus .ne. zero) then
	      retstat = error
	      call lib$signal(FDQ_CTOPENERR,%val(2),%val(iostatus),
	1	DATASETS(hkp))
	    endif
	  Endif		! retstat is success
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								 !
!     Now open FIRAS FEX_AV_CALRS file for read.                          !
!								 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  If ( retstat .eq. success ) Then

	    arc_file = arc_ref // cal_rdl // '/' // TIME_RANGE
	    OPEN ( unit=ct_lun(cal),
	1	   file=arc_file, 
	2	   status='OLD', iostat=iostatus,
	3	   useropen=CT_Connect_Read)
 
	    if (iostatus .ne. zero) then
	      retstat = error
	      call lib$signal(FDQ_CTOPENERR,%val(2),%val(iostatus),
	1	DATASETS(cal))
	    endif
	  Endif		! retstat is success

	ELSE				! Not the first time around


	   CALL CT_CLOSE_ARCV (, CT_LUN(IDX), CT_STAT)
	   IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	      RETSTAT = ERROR
	      STATUS = CT_STAT(1)
	      CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(IDX))
	   ENDIF

	   CALL CT_CLOSE_ARCV (, CT_LUN(ENG), CT_STAT)
	   IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	      RETSTAT = ERROR
	      STATUS = CT_STAT(1)
	      CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1	DATASETS(ENG))
	   ENDIF

	   DO CHAN=1,4

	      IF (LAST_SCI_OPEN(CHAN)) THEN

	        CALL CT_CLOSE_ARCV (, CT_LUN(CHAN), CT_STAT)
	        IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	          RETSTAT = ERROR
	          STATUS = CT_STAT(1)
	          CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1			DATASETS(CHAN))
	        ENDIF

	        CALL CT_CLOSE_ARCV (, CT_LUN(CHAN+18), CT_STAT)
	        IF (CT_STAT(1) .NE. CTP_NORMAL) THEN
	          RETSTAT = ERROR
	          STATUS = CT_STAT(1)
	          CALL LIB$SIGNAL(FDQ_CTCLOSERR,%VAL(2),%VAL(STATUS),
	1			DATASETS(CHAN+18))
	        ENDIF

	        
	      ENDIF

	   END DO

	ENDIF		! FIRST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                    !
!	Open all science files for modify.           !
!                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	DO	CHAN = 1, 4
	  SCI_OPEN(CHAN) = .FALSE.
	ENDDO

	DO	CHAN = 1, 4
	  If (filenames(chan) .ne. blank .and. retstat.eq.success) Then

	    SCI_OPEN(CHAN) = .TRUE.

	    arc_file = arc_id_in // filenames(chan)	    
	    OPEN ( unit=ct_lun(chan),
	1	   file=arc_file,
	2	   status='OLD', iostat=iostatus,
	3	   useropen=CT_Connect_Read)
 
	    if (iostatus .ne. zero) then
	      retstat = error
	      call lib$signal(FDQ_CTOPENERR,%val(2),%val(iostatus),
	1	DATASETS(chan))
	    endif

	    arc_file = arc_id_out // out_sci_rdl // fac_channel_ids(chan) //
	1		file_ext
            outnames(CHAN) = arc_file
	    OPEN ( unit=ct_lun(chan+18),
	1	   file=arc_file,
	2	   status='NEW', iostat=iostatus,
	3	   useropen=CT_Connect_Write)
 
	    if (iostatus .ne. zero) then
	      retstat = error
	      call lib$signal(FDQ_CTOPENERR,%val(2),%val(iostatus),
	1	DATASETS(chan+18))
	    endif


	  Endif			! filename is not blank & retstat is success

	ENDDO

	DO	CHAN = 1, 4
	  LAST_SCI_OPEN(CHAN) = SCI_OPEN(CHAN)
	ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!     Open IDX (Raw Science Index) for write     !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  If ( retstat .eq. success ) Then

	    arc_file = arc_id_out // idx_rdl // file_ext
            outnames(7)=arc_file
	    OPEN ( unit=ct_lun(idx),
	1	   file=arc_file, shared,
	2	   status='NEW', iostat=iostatus,
	3	   useropen=CT_Connect_Write)
 
	    if (iostatus .ne. zero) then
	      retstat = error
	      call lib$signal(FDQ_CTOPENERR,%val(2),%val(iostatus),
	1	DATASETS(idx))
	    endif
	  Endif			! retstat is success

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!							!
!     Now open Engineering Archive for write-access     !
!							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  If ( retstat .eq. success ) Then

	    arc_file = arc_id_out // eng_rdl // file_ext
            outnames(5)=arc_file
	    OPEN ( unit=ct_lun(eng),
	1	   file=arc_file, shared,
	2	   status='NEW', iostat=iostatus,
	3	   useropen=CT_Connect_Write)
 
	    if (iostatus .ne. zero) then
	      retstat = error
	      call lib$signal(FDQ_CTOPENERR,%val(2),%val(iostatus),
	1	DATASETS(eng))
	    endif
	  Endif			! retstat is success

!	Set function to return status

	if (retstat.eq.success) then
	  FDQ_OPEN_ARCV = %loc(FDQ_NORMAL)
	else
	  FDQ_OPEN_ARCV = %loc(FDQ_ABERR)
	endif

	return
	end
