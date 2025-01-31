	Integer*4 Function FSD_Astroplots_Read(scitype,chan,file_seg,
	2				       eof,sci_rec)

C---------------------------------------------------------------------
C
C	This subroutine reads science data from archive
C
C---------------------------------------------------------------------
C       Configured in CSDR environment by Reid Wilson,  STX Inc., 31-JAN-1988
C
C       version 4.2.1 11/24/88, SPR 2310, QUOC CHUNG, STX
C       BRING FSD UP TO STANDARD ERROR MESSAGE TRAPPING,
C       AND EXIT STATUS.
C
C	version 4.2.1  Enable plotting vs sky.  (SER 2373)
C	Fred Shuman, STX, 1989 Jan 25
C
C	SPR 5197, Archive open failure due to using old interface, CT_OPEN_ARCV.
C	    CT_OPEN_ARCV does not know about archives FPP_SDF and FDQ_SDF.
C	    R. Kummerer, STX / 1990 January 30.
C---------------------------------------------------------------------

	Implicit None

	Dictionary 'NFS_SDF'
	Record /NFS_SDF/ sci_rec

	Include   'ct$library:ctuser.inc'
	Include   '(CCT_Status_Record)'
	Include   '($SSDEF)'
	Include   '(FSD_Astroplots)'

	Record / CCT_Status / struc

	Integer   *  4      chan
	Character *  39     file_seg
	Logical   *  4      eof
	Character *  3      scitype

	Character *  39     blank / '                                       ' /
	Character *  64     filespec
	Integer   *  2      ios           !io status
	Integer   *  2      ctstat(20)	  !CT return status
	Integer   *  4      ct_nr
	Logical   *  4      opensuccess
	Integer   *  4      status
	Integer   *  4      LIB$Get_Lun
	Integer   *  4      CT_Connect_Read

	External            FSD_Normal
	External            FSD_CTOpen
	External            FSD_CTRead
	External            FSD_CTClos
	External            CT_Connect_Read


	FSD_Astroplots_Read = %loc(FSD_Normal)

	opensuccess = .True.

	If (begin) Then
c
c   Open the proper archive.
c
	   status = LIB$Get_Lun(ct_nr)
	   If (status .Ne. SS$_Normal) Then
	      Call Exit(status)
	   End If

	   If (file_seg .Eq. blank) Then
c
c   User supplied a timerange, not a filename.  Open proper archive by CT_open.
c
	      filespec = 'CSDR$FIRAS_ARCHIVE:' // scitype // '_SDF_' //
	2		 fac_channel_ids(chan) // '/' // time_range

	      Open (UNIT = ct_nr,
	2           FILE = filespec,
	3           STATUS = 'old',
	4           USEROPEN = CT_Connect_Read,
	5           IOSTAT = ios)

	      status = ios

	      If (ios .Ne. 0) Then
	         opensuccess = .False.
	      End If

	   Else
c
c   User supplied a filename, not a timerange.  Open with USEROPEN.
c
	      filespec = 'CSDR$FIRAS_ARCHIVE:' // file_seg

	      Open (UNIT = ct_nr,
	2           FILE = filespec,
	3           STATUS = 'old',
	4           USEROPEN = CT_Connect_Read,
	5           IOSTAT = ios)

	      status = ios

	      If (ios .Ne. 0) Then
	         opensuccess = .False.
	      End If

	   End If  !(file_seg .Eq. blank)

	End If     !begin

	If (more) Then
c
c   Get next data record.
c
	   If (opensuccess) Then
	      begin = .False.
	      Call CT_Read_Arcv(,ct_nr, sci_rec, ctstat)
c
c   If CT_Read wasn't normal, we will assume it encountered an End-Of-File.
c
	      If (ctstat(1) .Ne. CTP_Normal) Then
	         more = .False.
	         eof = .True.
	      End If
	   Else
	      more = .False.
	      Call Lib$Signal(FSD_CTOpen,%val(2),%val(chan),%val(status))
	      status = SS$_abort
	      Call Exit(status)
	   End If

	End If     ! (more)
c
c   End of timerange or section.
c
	If (.Not. more) Then

	   Call CT_Close_Arcv(, ct_nr, struc)
	   If (struc.cterr .Ne. CTP_Normal) Then
	      status = struc.cterr
	      Call Lib$signal(FSD_CTClos,%val(1),%val(status))
	   Endif

	End If     ! (.Not. more)

	End
