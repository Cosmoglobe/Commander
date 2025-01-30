	Program FFL

c-------------------------------------------------------------------------------
c
c	FFL_FISHINPUT_LONG
c
c	This program drives the FISH spectrum generation routines.  The program
c	reads coadded calibration IFGs from the Cobetrieve archives and writes
c	FISH-format voltage spectra to an RMS file.  The program runs for a
c	single channel and scan mode.  The program can read coadds from two
c	input archives and write a hybrid output to the RMS file.
c
c	Author:   Gene Eplee
c		  General Sciences Corporation
c		  513-7768
c		  9 March 1993
c		  SER 10763
c
c-------------------------------------------------------------------------------
c	Subroutines called:
c		FFL_Parse
c		FFL_Initialize_Report
c		FFL_Get_Reference
c		FFL_Open_Coadd
c		FFL_Read_Hybrid_Coadd
c		FFL_Read_Coadd
c		FFL_Produce_Spectra
c		FFL_Update_Report
c		cct_close_config
c		ct_init
c		cut_display_banner
c		cut_register_version
c		fut_free_lun
c		LIB$Establish
c		LIB$Signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fut_error.txt
c		fut_params.txt
c		ffl_config.txt
c		ffl_invoc.txt
c		ffl_spec.txt
c-------------------------------------------------------------------------------
c   Changes:
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11690
c
c	Name of facility changed from FFI to FFL.
c	Fred Shuman, HSTX, 1995 May 19. SPR 12272.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c
c	Eliminate routine FFL_Write_Spec by incorporating its code into here;
c	   phase_shift and var_norm deleted from include file ffl_spec.txt.
c	Fred Shuman, HSTX, 1995 June 30.
c
c	Subsumed the FFL_SPEC.TXT include file into FFL_CONFIG.TXT .
c	Fred Shuman, HSTX, 1995 July 21.
c
c	Removed array rwksp, Common /worksp/ rwksp, and call to iwkin, allowing
c	   IMSL to allocate workspace automatically.
c	Fred Shuman, HSTX, 1995 Aug 11.
c-------------------------------------------------------------------------------

	Implicit None

	Include 'ct$library:ctuser.inc'
	Include '(fut_error)'
	Include '(fut_params)'		! defines FIRAS parameters, fac_*
	Include '(ffl_config)'		! defines record fflc (structure config)
		!  and certain other vbles; Includes CCT_Get_Config; Dictionary
		!  structures FEX_GRTCOAWT, FEX_GRTTRANS.
		!defines structure fishspec and parameters vrecsize and vspecsiz
	Include '(ffl_invoc)'		! defines record ffli (structure invoc)
					!   and parameter "version"
	Record /fishspec/ ffls(fac_max_coad)

	Character * 79	cmd_line(3)	!  command line invocation
	Character * 14	current_gmt	!  GMT time of invocation

	Integer *  2	ct_stat(20)	!  Cobetrieve return status

	Integer *  4	cmdlen(3)	!  length of command lines in invocation
	Integer *  4	cstatus		!  Close_Config return status
	Integer *  4	ct_lun1		!  coadd file lun
	Integer *  4	ct_lun2		!  coadd file lun
	Integer *  4	io_stat		!  I/O return status
	Integer *  4	j		!  a counter
	Integer *  4	lun_out /6/	!  terminal lun
	Integer *  4	ncmd		!  number of command lines in invocation
	Integer *  4	num		!  number of current coadds read
	Integer *  4	parse_status	!  cmd line parse return status
	Integer *  4	read_status	!  coadd read return status
	Integer *  4	rstatus, status	!  return statuses
	Integer *  4	vs_lun		!  RMS file lun

	Logical *  1	ref_open_dir	!  reference dataset open flag
	Logical *  1	ref_open_seq	!  reference dataset open flag

	Integer * 4	FFL_Parse
	Integer * 4	FFL_Initialize_Report
	Integer * 4	FFL_Get_Reference
	Integer * 4	FFL_Open_Coadd
	Integer * 4	FFL_Read_Hybrid_Coadd
	Integer * 4	FFL_Read_Coadd
	Integer * 4	FFL_Produce_Spectra
	Integer * 4	FFL_Update_Report

	Integer * 4	cct_close_config
	Integer * 4	cut_display_banner, cut_register_version
	Integer * 4	fut_free_lun

	Dictionary 'fil_sky'
	Record /fil_sky/ coadd_recs(fac_max_coad)

	External	cct_normal
	External	ffl_chanerr
	External	ffl_cfgdirclose
	External	ffl_cfgseqclose
	External	ffl_ctinit
	External	ffl_eof
	External	ffl_failure
	External	ffl_normal
	External	ffl_repclose
	External	ffl_rmsclose
	External	ffl_rmswrite
	External	fut_normal, fut_error

C
C  Initialize the program.
C

c
c  Parse the command line.
c
	parse_status = FFL_Parse (current_gmt, ncmd, cmd_line, cmdlen, ffli)

c
c  Print the banner.
c
	rstatus = cut_register_version (version)
	rstatus = cut_display_banner (lun_out, 80,
     &					'FIRAS Facility FFL_Fishinput_Long')
	Write(lun_out,10)
 10	Format (/)

c
c  Initialize the processing report.
c
	If (ffli.report .Eq. fac_present) Then
	   status = FFL_Initialize_Report (current_gmt, ncmd, cmd_line, cmdlen,
     &					   parse_status, ffli)
	   if (status .Eq. %Loc(ffl_normal)) Call LIB$Establish (fut_error)
	Else
	   status = %Loc(ffl_normal)
	EndIf

	If (parse_status .Eq. %Loc(ffl_normal)  .And.
     &	    status .Eq. %Loc(ffl_normal)) Then
c
c  Initialize Cobetrieve.
c
	   Call ct_init(ct_stat)

	   If (ct_stat(1) .Eq. ctp_normal) Then
c
c  Get the reference datasets.
c
	      status = FFL_Get_Reference(ref_open_dir, ref_open_seq, fflc, ffli)


	      If (status .Eq. %Loc(ffl_normal)) Then
C
C  Process the spectra for the specified channel and scan mode.
C

c
c  Open the coadd file.
c
	         status = FFL_Open_Coadd (ct_lun1, ct_lun2, vs_lun, ffli)

	         If (status .Eq. %Loc(ffl_normal)) Then
c
c  Initialize the coadd read.
c
	            ffli.nspec  = 0
	            ffli.nspec1 = 0
	            ffli.nspec2 = 0
	            num         = 0
	            read_status  = %Loc(ffl_normal)

	            type 20, fac_channel_ids(ffli.chan),
     &			     fac_scan_mode_idsl(ffli.smode)
 20	            Format (x, 'Processing spectra for channel ', a,
     &			       ', scan mode ', a, /)

c
c  Read in a set of data for a particular scan mode.
c
	            Do While (read_status .Eq. %Loc(ffl_normal))

	               If (ffli.hybrid .Eq. fac_present) Then
	                  read_status = FFL_Read_Hybrid_Coadd (ct_lun1, ct_lun2,
     &						    num, coadd_recs, ffli)
	               Else
	                  read_status = FFL_Read_Coadd (ct_lun1, num,
     &							coadd_recs, ffli)
	               EndIf

	               If ( read_status .Eq. %Loc(ffl_normal)  .Or.
     &                     (read_status .Eq. %Loc(ffl_eof)  .And.  num .Gt. 0) )
     &			    Then

	                  status = %Loc(ffl_normal)
	                  j = 0
	                  ffli.nspec = ffli.nspec + num

	                  Do While (status .Eq. %Loc(ffl_normal)  .And.
     &                              j .Lt. num)
	                     j = j + 1
c
c  Produce the voltage spectra.
c
	                     status = FFL_Produce_Spectra (coadd_recs(j),
     &                                                     fflc, ffli, ffls(j))

	                  EndDo		!  (While j .Lt. num

c
c  Write the voltage spectra to the FISH-format RMS file.
c
	                  If (status .Eq. %Loc(ffl_normal)) Then

	                     j = 1
	                     Do While (status .Eq. %Loc(ffl_normal)  .And.
     &                                 j .Le. num)

	                        Write (unit=vs_lun, iostat=io_stat) ffls(j)

	                        If (io_stat .Ne. 0) Then
	                           status = %Loc(ffl_rmswrite)
	                           Call LIB$Signal (ffl_rmswrite, %Val(2),
     &						    ffli.outfile(1:ffli.outlen),
     &						    %Val(io_stat))
	                        EndIf

	                        j = j + 1
	                     End Do

	                     If (read_status .Eq. %Loc(ffl_eof)) Then
	                        Close (unit=vs_lun, iostat=io_stat)
	                        If (io_stat .Ne. 0) Then
	                           status = %Loc(ffl_rmsclose)
	                           Call LIB$Signal (ffl_rmsclose, %Val(2),
     &						    ffli.outfile(1:ffli.outlen),
     &						    %Val(io_stat))
	                        EndIf
	                        rstatus = fut_free_lun (vs_lun)
	                        If (rstatus .Ne. %Loc(fut_normal)) Then
	                           Call LIB$Signal (%Val(rstatus))
	                        EndIf
	                     EndIf

	                  EndIf		!  (status from processing spectra

c
c  Close the FISH-format files if no records were read the last time.
c
	               ElseIf (read_status .Eq. %Loc(ffl_eof)  .And.
     &			       num .Eq. 0) Then
	                  Close (unit=vs_lun, iostat=io_stat)
	                  If (io_stat .Ne. 0) Then
	                     status = %Loc(ffl_rmsclose)
	                     Call LIB$Signal (ffl_rmsclose, %Val(2),
     &					      ffli.outfile(1:ffli.outlen),
     &					      %Val(io_stat))
	                  EndIf
	                  rstatus = fut_free_lun (vs_lun)
	                  If (rstatus .Ne. %Loc(fut_normal)) Then
	                     Call LIB$Signal (%Val(rstatus))
	                  EndIf
	               Else
	                  status = %Loc(ffl_chanerr)
	                  Call LIB$Signal (ffl_chanerr)
	               EndIf	!  (read_status from FFL_Read_Coadd

	            EndDo	!  (Do While read_status


C
C  Close out the program.
C

c
c  Update the processing report with the number of spectra processed.
c
	            If (status .Eq. %Loc(ffl_normal)  .And.
     &	                ffli.report .Eq. fac_present) Then
	               status = FFL_Update_Report (ffli)
	            EndIf

c
c  Write out number of spectra processed.
c
	            If (status .Eq. %Loc(ffl_normal)) Then
	               If (ffli.hybrid .Eq. fac_present) Then
	                  Write (lun_out,30,iostat=io_stat) ffli.nspec,
     &							ffli.nspec1, ffli.nspec2
	               Else
	                  Write (lun_out,40,iostat=io_stat) ffli.nspec
	               EndIf
	            EndIf
  30	Format (/, 4x, 'Number of calibration spectra processed:  ', I5,
     &		/, 4x, '    Number from input 1:                  ', I5,
     &		/, 4x, '    Number from input 2:                  ', I5, //)
  40	Format (/, 4x, 'Number of calibration spectra processed:  ', I5, //)

	         EndIf	!	(status from FFL_Open_Coadd

	      EndIf	!	(status from FFL_Get_Reference

c
c  Close the configuration files.
c
	      If (ref_open_seq) Then
	         cstatus = cct_close_config (ndset_tod, fflc.tod_lun,
     &                                       fflc.tod_index)
	         If (cstatus .Ne. %Loc(cct_normal)) Then
	            status = %Loc(ffl_cfgseqclose)
	            Call LIB$Signal (ffl_cfgseqclose, %Val(1), %Val(cstatus))
	         EndIf
	      EndIf

	      If (ref_open_dir) Then
	         cstatus = cct_close_config (ndset_dir, fflc.dir_lun,
     &                                       fflc.dir_index)
	         If (cstatus .Ne. %Loc(cct_normal)) Then
	            status = %Loc(ffl_cfgdirclose)
	            Call LIB$Signal (ffl_cfgdirclose, %Val(1), %Val(cstatus))
	         EndIf
	      EndIf

	   Else
	      status = %Loc(ffl_ctinit)
	      Call LIB$Signal (ffl_ctinit, %Val(1), %Val(ct_stat(1)))
	   EndIf		!  (ct_stat from ct_init

	EndIf			!  (status from FFL_Parse and FFL_Init*_Report

c
c  Signal the program completion status.
c
	If (status .Eq. %Loc(ffl_normal)  .Or.  status .Eq. %Loc(ffl_chanerr))
     &	      Then
	   Call LIB$Signal (ffl_normal)
	Else
	   Call LIB$Signal (ffl_failure)
	EndIf

c
c  Close the processing report file.
c
	If (ffli.report .Eq. fac_present) Then
	   Close (unit=fut_report_lun, iostat=io_stat)
	   If (io_stat .Ne. 0) Then
	      Call LIB$Signal (ffl_repclose, %Val(2),
     &			       ffli.report_file(1:ffli.replen), %Val(io_stat))
	   EndIf
	   rstatus = fut_free_lun (fut_report_lun)
	   If (rstatus .Ne. %Loc(fut_normal)) Then
	      Call LIB$Signal (%Val(rstatus))
	   EndIf
	EndIf


	End
