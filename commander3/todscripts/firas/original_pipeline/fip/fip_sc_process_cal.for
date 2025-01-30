	Integer*4  Function  FIP_SC_Process_Cal ( input, filexts, fnum, archin,
	1                                         archout, nrec )

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Function  FIP_SC_PROCESS_CAL.FOR
C
C  Function to convert the data in the case of calibration data.  Input is
C  either FCF_CAL_xxxx, FCF_DCL_xxxx, or FIC_CAL_xxxx.  Uses Cobetrieve to
C  read and write data.  If input is FIC_CAL, then call use FIP_REFORMAT_IFG
C  else use FIP_REFORMAT_SPECTRUM.
C
C  Passed parameters:
C     Input:
C	Character*12	input			! Input filename base
C	Character*20	filexts (fac_max_num)	! Input file extensions
C	Integer*2	fnum			! Number of input files
C	Character*13	archin			! Input archive name
C	Character*14	archout			! Output archive name
C     Output:
C	Integer*4	nrec			! Number of records processed
C
C  Author:  Larry P. Rosen, Hughes STX, 30 November 1994.
C
C PDL:
C  Open the input and output data files with FIP_SC_Open_Cal.
C  Loop over all files.
C     Loop over all records.
C        Set the frequency limits the first time through.
C        Read the record.
C        Reformat the record.  If coadds use FIC structures and
C        FIP_REFORMAT_IFG, else use fcf structures and FIP_REFORMAT_SPECTRUM.
C        Write out the reformated records to the ADB skymap.
C  Close the files with FIP_SC_Close_Cal.
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Implicit None

C  Include files:
	Include		'(fut_error)'
	Include		'(fut_params)'
	Include		'(fip_invoc_sky)'
	Include		'ct$library:ctuser.inc'

C  Passed parameters:

	Character*12	input			! Input filename base
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			! Number of input files
	Character*13	archin			! Input archive name
	Character*14	archout			! Output archive name
	Integer*4	nrec			! Number of records processed

C  Functions:
	Integer*4	FIP_SC_Open_Cal
	Integer*4	FIP_Reformat_Spectrum
	Integer*4	FIP_Reformat_IFG
	Integer*4	FIP_Frequency_Cut
	Integer*4	FIP_SC_Close_Cal

C  Externals:
	External	FIP_Normal
	External	FIP_Maxrec
	External	FIP_RMSRead
	External	FIP_Writerr

C  Local:
	Integer*4	status		! Processing status
	Integer*4	inlun		! Input lun
	Integer*4	outlun		! Output lun
	Logical*1	first_time	! first time flag for frequency cut
	Real*4		fnyq_hz		! Nyquist frequency in hz
	Real*4		fnyq_icm	! Nyquist frequency in icm
	Integer*4	j		! Counter
	Integer*4	fil		! File counter
	Logical*1	eof		! End of file flag
	Integer*2	ct_stat (20)	! CT return status
	Character*33	infile			!  Input file name
	Character*1	blanks(33) /33*' '/
	Character*33	blankname		!  Blank input file name
	Equivalence	(blanks, blankname)
	Integer*4	io_status

C  Record structures for input and output data.  Note that fcf_cal and fcf_dcl
C  use the same record structure (fcf_cal), and that their corresponding output
C  files, fip_cal and fip_dcl, use the same data structure as well.

	Dictionary 'fcf_cal'
	Record /fcf_cal/	in_fcf
	Dictionary 'fic_cal'
	Record /fic_cal/	in_fic
	Dictionary 'fip_ccl'
	Record /fip_ccl/	out_fcf
	Dictionary 'fip_icl'
	Record /fip_icl/	out_fic

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Begin

	FIP_SC_Process_Cal = %loc (FIP_Normal)
	nrec = 0

C  Open the input and output data files.

	status = FIP_SC_Open_Cal ( input, filexts, fnum, archin, archout,
	1                          inlun, outlun )

C  Main loop over files.

	first_time = .TRUE.
	fil = 1
	Do While ( fil .LE. fnum .AND. status .EQ. %loc(FIP_Normal) )

C  Loop over records in each file.

	   infile = blankname
	   infile = input // '.' // filexts (fil)
	   eof = .FALSE.
	   Do While ( status .EQ. %loc(FIP_Normal) .AND. .NOT. eof )
	      If (fcc_coadd) Then			! FIC SKY data
	         Call CT_Read_Arcv (, inlun, in_fic, ct_stat )
	      Else
	         Call CT_Read_Arcv (, inlun, in_fcf, ct_stat )
	      EndIf
	      If (ct_stat(1) .NE. ctp_normal) Then
	         If (ct_stat(1) .EQ. ctp_endoffile ) Then
	            eof = .TRUE.
	         Else		! read error
	            status = ct_stat(1)
	            Call Lib$Signal (FIP_rmsread,%val(2),infile,%loc(status))
	         Endif		! read error
	      Else		! normal return from ct_read_arcv
	         nrec = nrec + 1

C  Set the frequency limits the first time through.

	         If (first_time) Then
	            first_time = .FALSE.
	            If (fcc_coadd) Then		! FIC CAL data
	               fcc_destriped = 0
	               fnyq_icm = in_fic.coad_spec_data.nyquist_icm
	               fnyq_hz = in_fic.coad_spec_data.nyquist_hertz
	            Else				! FCF CAL or DCL data
	               fcc_destriped = in_fcf.spec_data.destriped
	               fnyq_icm = in_fcf.coad_spec_data.nyquist_icm
	               fnyq_hz = in_fcf.coad_spec_data.nyquist_hertz
	            EndIf
	            status = FIP_Frequency_Cut (fnyq_icm, fnyq_hz)
	         EndIf
	         If (fcc_coadd) Then			! FIC CAL data
	            status = FIP_Reformat_IFG (in_fic, out_fic)
	         Else				! FCF CAL or DCL data
	            status = FIP_Reformat_Spectrum (in_fcf, out_fcf)
	         EndIf
	         If (status .EQ. %loc(FIP_Normal)) Then

C Write record to output file.

	            If (fcc_coadd) Then			! FIC CAL data
	               Write ( outlun, iostat=io_status ) out_fic
	            Else				! FCF CAL or DCL data
	               Write ( outlun, iostat=io_status ) out_fcf
	            EndIf
	            If (io_status .NE. 0) Then
	               status = %loc (io_status)
	               Call Lib$Signal (fip_writerr, %Val(2),
	1                               fcc_outfile(1:fcc_outlen),
	2                               %loc (io_status))
	            EndIf
	         EndIf
	      EndIf
	   EndDo			! Do for all records in file
	   fil = fil + 1
	EndDo				! Do for all input files

C  Close the files.

	If (status .EQ. %loc(FIP_Normal)) Then
	   status = FIP_SC_Close_Cal ( input, filexts, fnum, archin, inlun,
	1                              outlun)
	EndIf

	FIP_SC_Process_Cal = status

	Return
	End
