	Integer*4  Function  FFP_SC_Process_Cal ( input, filexts, fnum, archin,
	1                                         archout, nrec )
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Function  FFP_SC_PROCESS_CAL.FOR
C
C  Function to convert the data in the case of calibration data.  Input is
C  either FSL_CAL_xxxx, FSL_DCL_xxxx, or FIL_CAL_xxxx.  Uses Cobetrieve to
C  read and write data.  If input is FIL_CAL, then use FFP_REFORMAT_IFG
C  else use FFP_REFORMAT_SPECTRUM.
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
C  Author:  S. Brodd, HSTX, 3/21/96
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Implicit None

C  Include files:
	Include		'(fut_error)'
	Include		'(fut_params)'
	Include		'(ffp_invoc_sky)'
	Include		'ct$library:ctuser.inc'

C  Passed parameters:

	Character*12	input			! Input filename base
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			! Number of input files
	Character*13	archin			! Input archive name
	Character*14	archout			! Output archive name
	Integer*4	nrec			! Number of records processed

C  Functions:
	Integer*4	FFP_SC_Open_Cal
	Integer*4	FFP_Reformat_Spectrum
	Integer*4	FFP_Reformat_IFG
	Integer*4	FFP_Frequency_Cut
	Integer*4	FFP_SC_Close_Cal

C  Externals:
	External	FFP_Normal
	External	FFP_Maxrec
	External	FFP_RMSRead
	External	FFP_Writerr

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

C  Record structures for input and output data.  Note that fsl_cal and fsl_dcl
C  use the same record structure (fsl_cal), and that their corresponding output
C  files, ffp_cal and ffp_dcl, use the same data structure as well.

	Dictionary 'fsl_sky'
	Record /fsl_sky/	in_fsl
	Dictionary 'fil_sky'
	Record /fil_sky/	in_fil
	Dictionary 'ffp_spc'
	Record /ffp_spc/	out_fsl
	Dictionary 'ffp_ifg'
	Record /ffp_ifg/	out_fil

C  Begin

	FFP_SC_Process_Cal = %loc (FFP_Normal)
	nrec = 0

C  Open the input and output data files.

	status = FFP_SC_Open_Cal ( input, filexts, fnum, archin, archout,
	1                          inlun, outlun )

C  Main loop over files.

	first_time = .TRUE.
	fil = 1
	Do While ( fil .LE. fnum .AND. status .EQ. %loc(FFP_Normal) )

C  Loop over records in each file.

	   infile = blankname
	   infile = input // '.' // filexts (fil)
	   eof = .FALSE.
	   Do While ( status .EQ. %loc(FFP_Normal) .AND. .NOT. eof )
	      If (fcc_coadd) Then			! FIL CAL data
	         Call CT_Read_Arcv (, inlun, in_fil, ct_stat )
	      Else
	         Call CT_Read_Arcv (, inlun, in_fsl, ct_stat )
	      EndIf
	      If (ct_stat(1) .NE. ctp_normal) Then
	         If (ct_stat(1) .EQ. ctp_endoffile ) Then
	            eof = .TRUE.
	         Else		! read error
	            status = ct_stat(1)
	            Call Lib$Signal (FFP_rmsread,%val(2),infile,%loc(status))
	         Endif		! read error
	      Else		! normal return from ct_read_arcv
	         nrec = nrec + 1

C  Set the frequency limits the first time through.

	         If (first_time) Then
	            first_time = .FALSE.
	            If (fcc_coadd) Then		! FIL CAL data
	               fcc_destriped = 0
	               fnyq_icm = in_fil.coad_spec_data.nyquist_icm
	               fnyq_hz = in_fil.coad_spec_data.nyquist_hertz
	            Else				! FSL CAL or DCL data
	               fcc_destriped = in_fsl.spec_data.destriped
	               fnyq_icm = in_fsl.coad_spec_data.nyquist_icm
	               fnyq_hz = in_fsl.coad_spec_data.nyquist_hertz
	            EndIf
	            status = FFP_Frequency_Cut (fnyq_icm, fnyq_hz)
	         EndIf
	         If (fcc_coadd) Then			! FIL CAL data
	            status = FFP_Reformat_IFG (in_fil, out_fil)
	         Else				! FSL CAL or DCL data
	            status = FFP_Reformat_Spectrum (in_fsl, out_fsl)
	         EndIf
	         If (status .EQ. %loc(FFP_Normal)) Then

C Write record to output file.

	            If (fcc_coadd) Then			! FIL CAL data
	               Write ( outlun, iostat=io_status ) out_fil
	            Else				! FSL CAL or DCL data
	               Write ( outlun, iostat=io_status ) out_fsl
	            EndIf
	            If (io_status .NE. 0) Then
	               status = %loc (io_status)
	               Call Lib$Signal (ffp_writerr, %Val(2),
	1                               fcc_outfile(1:fcc_outlen),
	2                               %loc (io_status))
	            EndIf
	         EndIf
	      EndIf
	   EndDo			! Do for all records in file
	   fil = fil + 1
	EndDo				! Do for all input files

C  Close the files.

	If (status .EQ. %loc(FFP_Normal)) Then
	   status = FFP_SC_Close_Cal ( input, filexts, fnum, archin, inlun,
	1                              outlun)
	EndIf

	FFP_SC_Process_Cal = status

	Return
	End
