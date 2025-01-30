	Integer * 4 Function  FFL_Read_Hybrid_Coadd (ct_lun1, ct_lun2, num,
     &						     coadd_recs, ffli)

c-------------------------------------------------------------------------------
c
c	Function FFL_READ_HYBRID_COADD
c
c	This function reads in coadded IFGs from the Cobetrieve calibration
c	archives.  Each time through, the function reads one to FAC_MAX_COAD
c	coadd records from the cal archive.  The routine decides which input
c	stream coadd to place in the input buffer.  The last time through, the
c	function closes the Cobetrieve archive.
c
c	Author:	 Gene Eplee
c		 General Sciences Corp.
c		 513-7768
c		 9 March 1993
c		 SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c		ct_lun1		Integer * 4	coadd file lun
c		ct_lun2		Integer * 4	coadd file lun
c		ffli		Record Structure defined in FFL_Invoc.txt
c
c	Output:
c		num		Integer * 4	number of coadds read
c		coadd_recs	coadd records	input coadd records
c		ffli		Record Structure defined in FFL_Invoc.txt
c
c	Subroutines called:
c		FFL_Bracket_Coadd
c		fut_free_lun
c		ct_close_arcv
c		LIB$Signal
c
c	Include files:
c		ffl_invoc.txt
c		fut_params.txt
c		ct$library:ctuser.inc
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Name of facility changed from FFI to FFL.  Fred Shuman, HSTX
c	1995 May 19.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c
c	Read peakpos from the 1st coadd record rather than calling
c	   fut_default_peak to compute it.  Fred Shuman, HSTX, 1995 June 29.
c
c-------------------------------------------------------------------------------

	Implicit None

	Include 'ct$library:ctuser.inc'
	Include '(fut_params)'
	Include '(ffl_invoc)'		! defines call argument (structure) ffli
c
c  Other call arguments:
c
	Integer * 4	ct_lun1		!  CT logical unit number
	Integer * 4	ct_lun2		!  CT logical unit number
	Integer * 4	num		!  number of coadd records
	Dictionary 'fil_sky'
	Record /fil_sky/ coadd_rec1
	Record /fil_sky/ coadd_rec2
	Record /fil_sky/ coadd_recs(fac_max_coad)
c
c  All other variables and functions:
c
	Integer * 2	ct_stat1(20)	!  CT return status
	Integer * 2	ct_stat2(20)	!  CT return status

	Integer * 4	rstatus		!  return status
	Integer * 4	status		!  return status
	Integer * 4	peakpos		!  peak position read from FIL record

	Logical * 1	bracket /.True./ !  bracketing coadd flag

	Real	* 4	peakht		!  ifg peak height

	Integer * 4	FFL_Bracket_Coadd
	Integer * 4	fut_free_lun

	External	ffl_ctifgclose
	External	ffl_eof
	External	ffl_normal
	External	fut_normal

C
C  Read and filter the coadded IFGs.
C
	num = 0
	status = %Loc(ffl_normal)

	Do While ((status .Eq. %Loc(ffl_normal))  .And.
     &		  (num .Lt. fac_max_coad))

c
c  Read the coadd records.
c
	   status = FFL_Bracket_Coadd (ct_lun1, ct_lun2, coadd_rec1,
     &				       coadd_rec2, bracket, ffli)

	   If (status .Eq. %Loc(ffl_normal)  .And.
     &	       coadd_rec1.coad_spec_data.dq_summary_flag .Le.
     &						    ffli.instr_qual  .And.
     &	       coadd_rec1.coad_spec_data.att_summary_flag .Le.
     &						    ffli.attit_qual) Then

c
c  Put the coadds (from the correct stream) into the coadd buffer.
c
	      num = num + 1
	      If ((bracket .Eq. .True.)  .And.
     &	          (coadd_rec2.coad_spec_data.dq_summary_flag .Le.
     &					     ffli.instr_qual)  .And.
     &	          (coadd_rec2.coad_spec_data.att_summary_flag .Le.
     &					     ffli.attit_qual)) Then
	         peakpos = coadd_rec1.coad_spec_data.peak_pos
	         peakht = Abs(coadd_rec1.coad_data.ifg(peakpos))
	         If (peakht .Gt. ffli.threshold) Then
	            coadd_recs(num) = coadd_rec2
	            ffli.nspec2 = ffli.nspec2 + 1
	         Else
	            coadd_recs(num) = coadd_rec1
	            ffli.nspec1 = ffli.nspec1 + 1
	         EndIf
	      Else
	         coadd_recs(num) = coadd_rec1
	         ffli.nspec1 = ffli.nspec1 + 1
	      EndIf

	   EndIf

	EndDo	!  While (status from read


	If (status .Eq. %Loc(ffl_eof)) Then
C
C  Close the archive after the final read.
C
	   Call ct_close_arcv (, ct_lun1, ct_stat1)
	   If (ct_stat1(1) .Ne. ctp_normal) Then
	      status = %Loc(ffl_ctifgclose)
	      Call LIB$Signal (ffl_ctifgclose, %Val(2),
     &			       ffli.infile1(1:ffli.inlen1), %Val(ct_stat1(1)))
	   EndIf

	   Call ct_close_arcv (, ct_lun2, ct_stat2)
	   If (ct_stat2(1) .Ne. ctp_normal) Then
	      status = %Loc(ffl_ctifgclose)
	      Call LIB$Signal (ffl_ctifgclose, %Val(2),
     &			       ffli.infile2(1:ffli.inlen2), %Val(ct_stat2(1)))
	   EndIf

c
c  Free the Cobetrieve logical unit number.
c
	   rstatus = fut_free_lun (ct_lun1)
	   rstatus = fut_free_lun (ct_lun2)
	   If (rstatus .Ne. %Loc(fut_normal)) Then
	      Call LIB$Signal (%Val(rstatus))
	   EndIf

	EndIf		!  (end of file


	FFL_Read_Hybrid_Coadd = status

	Return
	End
