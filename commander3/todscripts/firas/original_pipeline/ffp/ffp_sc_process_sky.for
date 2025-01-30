	Integer*4  Function  FFP_SC_Process_Sky ( input, filexts, fnum, archin,
	1                                         archout, nrec, npix )
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Function  FFP_SC_PROCESS_SKY.FOR
C
C  Function to convert the data in the case of skymap data.  Input is either
C  FSL_SKY_xxxx, FSL_DSK_xxxx, or FIL_SKY_xxxx.  Uses CSA routines to read and
C  write skymap data.  If input is FIL_SKY, then use FFP_REFORMAT_IFG
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
C	Integer*4	npix			! Number of pixels processed
C
C  Author:  S. Brodd, HSTX, 3/21/96
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Implicit None

C  Include files:
	Include 	'(csa_pixel_input_rec)'
	Include 	'(csa_pixel_output_rec)'
	Include 	'ct$library:ctuser.inc'
	Include 	'(fut_error)'
	Include 	'(fut_params)'
	Include		'(ffp_invoc_sky)'

C  Passed parameters:

	Character*12	input			! Input filename base
	Character*20	filexts (fac_max_num)	! Input file extensions
	Integer*2	fnum			! Number of input files
	Character*13	archin			! Input archive name
	Character*14	archout			! Output archive name
	Integer*4	nrec			! Number of records processed
	Integer*4	npix			! Number of pixels processed

C  Functions:
	Integer*4	FFP_SC_Open_Sky
	Integer*4	FFP_Galactic_Cut
	Integer*4	CSA_Read_Pixels
	Integer*4	FFP_Reformat_Spectrum
	Integer*4	FFP_Reformat_IFG
	Integer*4	FFP_Frequency_Cut
	Integer*4	CSA_Write_Pixels
	Integer*4	FFP_SC_Close_Sky

C  Externals:
	External	FFP_Normal
	External	CSA_Normal
	External	FFP_Maxrec
	External	FFP_CSAWrite
	External	FFP_CSARead

C  Local:
	Integer*4	status		! Processing status
	Integer*4	inlun (fac_max_num)	! Input lun's
	Integer*4	outlun		! Output lun
	Integer*4	cstatus		! CSA return status
	Integer*4	pstatus		! pixel status
	Integer*4	blocks		! block count for CSA
	Logical*1	first_time	! first time flag for frequency cut
	record /pixel_input_list/  inlist
	record /pixel_output_list/ outlist
	Integer*4	max_in		! maximun number of pixels to read
					!    in one call to CSA
	Integer*4	num_in		! number of pixels to be read
					!    in one call to CSA
	Integer*4	num_out		! number of records to be written
					!    in one call to CSA
	Integer*4	pixel_no	! pixel number to be read by CSA
	Integer*4	num_read	! number of pixels read by CSA
	Real*4		fnyq_hz		! Nyquist frequency in hz
	Real*4		fnyq_icm	! Nyquist frequency in icm
	Integer*4	fil, j          ! Counters
	Character*33	infile		! Input file name

C  Record structures for input and output data.  Note that fsl_sky and fsl_dsk
C  use the same record structure (fsl_sky), and that their corresponding output
C  files, ffp_csk and ffp_dsk, use the same data structure as well.

	Dictionary 'fsl_sky'
	Record /fsl_sky/	in_fsl (fac_max_skymap_recs)
	Dictionary 'fil_sky'
	Record /fil_sky/	in_fil (fac_max_skymap_recs)
	Dictionary 'ffp_spc'
	Record /ffp_spc/	out_fsl (fac_max_skymap_recs)
	Dictionary 'ffp_ifg'
	Record /ffp_ifg/	out_fil (fac_max_skymap_recs)

C  Begin

	FFP_SC_Process_Sky = %loc (FFP_Normal)
	npix = 0
	nrec = 0

C  Open the input and output data files.

	status = FFP_SC_Open_Sky ( input, filexts, fnum, archin, archout,
	1                          inlun, outlun )
	If (status .NE. %loc (FFP_Normal)) Then
	   FFP_SC_Process_Sky = status
	   Return
	EndIf

C  Loop over the pixels, reading and processing all records in each pixel.

	cstatus = %loc(CSA_Normal)
	pstatus = %loc(FFP_Normal)
	blocks = 0
	first_time = .TRUE.
	inlist.level_no = fac_skymap_level
	max_in  = fac_max_skymap_recs
	num_in  = 1
	num_out = 1
	pixel_no = -1

c Main loop over pixels

	Do While ( (cstatus .EQ. %loc(CSA_Normal)) .AND.
	1          (status .EQ. %loc(FFP_Normal))  .AND.
	2          (pstatus .EQ. %loc(FFP_Normal)) .AND.
	3          (pixel_no .LT. 6143) )

	   pixel_no = pixel_no + 1
	   If (fcc_galexc .EQ. fac_present) Then
	      pstatus = FFP_Galactic_Cut (pixel_no)
	   EndIf
	   If (pstatus .EQ. %loc(FFP_Normal)) Then
	      inlist.pixel_no = pixel_no

c File Loop: do for each input skymap file.

	      fil = 1
	      Do While ( fil .LE. fnum .AND. status .EQ. %loc(FFP_Normal) )

c Read any input records in this pixel

	         If (fcc_coadd) Then			! FIL SKY data
	            cstatus = CSA_Read_Pixels ( inlun (fil), inlist, num_in,
	1              in_fil, max_in, outlist, num_read, blocks )
	         Else					! FSL SKY or DSK data
	            cstatus = CSA_Read_Pixels ( inlun (fil), inlist, num_in,
	1              in_fsl, max_in, outlist, num_read, blocks )
	         EndIf
	         If (outlist.no_records .GT. fac_max_skymap_recs) Then
	            status = %loc(FFP_Maxrec)
	            Call Lib$Signal ( FFP_Maxrec, %val(2),
	1                             %val(outlist.no_records),
	2                             %val(pixel_no) )
	         ElseIf (outlist.no_records .NE. 0) Then
	            If (cstatus .EQ. %loc(CSA_Normal)) Then
	               npix = npix + 1
	               nrec = nrec + outlist.no_records

C  Set the frequency limits the first time through.

	               If (first_time) Then
	                  first_time = .FALSE.
	                  If (fcc_coadd) Then		! FIL SKY data
	                     fcc_destriped = 0
	                     fnyq_icm = in_fil(1).coad_spec_data.nyquist_icm
	                     fnyq_hz = in_fil(1).coad_spec_data.nyquist_hertz
	                  Else				! FSL SKY or DSK data
	                     fcc_destriped = in_fsl(1).spec_data.destriped
	                     fnyq_icm = in_fsl(1).coad_spec_data.nyquist_icm
	                     fnyq_hz = in_fsl(1).coad_spec_data.nyquist_hertz
	                  EndIf
	                  status = FFP_Frequency_Cut (fnyq_icm, fnyq_hz)
	               EndIf

C  Reformat the records.  If coadds use FIL structures and FFP_REFORMAT_IFG,
C  else use fsl structures and FFP_REFORMAT_SPECTRUM.

	               If (fcc_coadd) Then		! FIL SKY data
	                  j = 1
	                  Do While ( j .LE. outlist.no_records .AND.
	1                            status .EQ. %loc(FFP_Normal) )
	                     status = FFP_Reformat_IFG (in_fil(j), out_fil(j))
	                     j = j + 1
	                  EndDo

	                  If (status .EQ. %loc(FFP_Normal)) Then

C  Write out the reformated records to the ADB skymap.

	                     j = 1
	                     Do While ( (cstatus .EQ. %loc(CSA_Normal)) .AND.
	1                               (j .LE. outlist.no_records) )
	                        cstatus = CSA_Write_Pixels ( outlun,
	1                                                    out_fil(j),
	2                                                    num_out, blocks )
	                        j = j + 1
	                     EndDo
	                     If (cstatus .NE. %loc(CSA_Normal)) Then
	                        status = %loc(FFP_CSAWrite)
	                        call Lib$Signal ( FFP_CSAWrite, %val(2),
	1                                         fcc_outfile(1:fcc_outlen),
	2                                         %val(cstatus) )
	                     EndIf
	                  EndIf
	               Else				! FSL SKY or DSK data
	                  j = 1
	                  Do While ( j .LE. outlist.no_records .AND.
	1                            status .EQ. %loc(FFP_Normal) )
	                     status = FFP_Reformat_Spectrum (in_fsl(j),
	1                                                    out_fsl(j))
	                     j = j + 1
	                  EndDo

	                  If (status .EQ. %loc(FFP_Normal)) Then

C  Write out the reformated records to the ADB skymap.

	                     j = 1
	                     Do While ( (cstatus .EQ. %loc(CSA_Normal)) .AND.
	1                               (j .LE. outlist.no_records) )
	                        cstatus = CSA_Write_Pixels ( outlun,
	1                                                    out_fsl(j),
	2                                                    num_out, blocks )
	                        j = j + 1
	                     EndDo
	                     If (cstatus .NE. %loc(CSA_Normal)) Then
	                        status = %loc(FFP_CSAWrite)
	                        call Lib$Signal ( FFP_CSAWrite, %val(2),
	1                                         fcc_outfile(1:fcc_outlen),
	2                                         %val(cstatus) )
	                     EndIf
	                  EndIf
	               EndIf		!  fcc_coadd
	            Else
	               infile = input // '.' // filexts (fil)
	               status = %loc(FFP_CSARead)
	               call Lib$Signal ( FFP_CSARead, %val(2), infile,
	1                                %val(cstatus) )
	            EndIf		!  return status from read pixels
	         EndIf			!  outlist ne 0
	         fil = fil + 1
	      EndDo			!  Do for all input files
	   EndIf			!  If pixel exluded by galactic cut
	EndDo				!  Do for all pixels

C  Close the skymaps.

	If (status .EQ. %loc(FFP_Normal)) Then
	   status = FFP_SC_Close_Sky (inlun, outlun, fnum, input, filexts)
	EndIf

	FFP_SC_Process_Sky = status

	Return
	End
