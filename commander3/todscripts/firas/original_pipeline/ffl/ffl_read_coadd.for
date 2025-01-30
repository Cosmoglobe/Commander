	Integer * 4 Function  FFL_Read_Coadd (ct_lun1, num, coadd_recs, ffli)

c-------------------------------------------------------------------------------
c
c	Function FFL_READ_COADD
c
c	This function reads in coadded IFGs from the Cobetrieve calibration
c	archive.  Each time through, the function reads one up to FAC_MAX_COAD
c	coadd records from the cal archive.  The last time through, the function
c	closes the Cobetrieve archive.
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
c		ct_lun1		Integer * 4		coadd file lun
c		ffli		Record Structure defined in FFL_Invoc.txt
c
c	Output:
c		num		Integer * 4		number of coadds read
c		coadd_recs	coadd records		input coadd records
c		ffli		Record Structure defined in FFL_Invoc.txt
c
c	Subroutines called:
c		ct_close_arcv
c		ct_read_arcv
c		fut_free_lun
c		LIB$Signal
c
c	Include files:
c		ct$library:ctuser.inc
c		fut_params.txt
c		ffl_invoc.txt
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
c-------------------------------------------------------------------------------

	Implicit None

	Include 'ct$library:ctuser.inc'
	Include '(fut_params)'
	Include '(ffl_invoc)'
c
c  Call arguments:
c
	Integer * 4	ct_lun1		!  CT logical unit number
	Integer * 4	num		!  number of coadd records
	Dictionary 'fil_sky'
	Record /fil_sky/ coadd_recs(fac_max_coad)
c
c  All other variables and functions:
c
	Record /fil_sky/ coadd_rec
	Integer * 4	rstatus		!  return status
	Integer * 4	status		!  return status
	Integer * 4	cstatus		!  return status
	Integer * 2	ct_stat(20)	!  CT return status

	Integer * 4	fut_free_lun


	External	ffl_ctifgclose
	External	ffl_ctifgread
	External	ffl_eof
	External	ffl_normal
	External	fut_normal

C
C  Read and filter the coadded IFGs.
C
	num = 0
	ct_stat(1) = ctp_normal

	Do While (ct_stat(1) .Eq. ctp_normal  .And.  num .Lt. fac_max_coad)

c
c  Read the coadd records.
c
	   Call ct_read_arcv (, ct_lun1, coadd_rec, ct_stat)

	   If (ct_stat(1) .Eq. ctp_normal  .And.
     &	       coadd_rec.coad_spec_data.dq_summary_flag .Le.
     &						    ffli.instr_qual  .And.
     &	       coadd_rec.coad_spec_data.att_summary_flag .Le.
     &						     ffli.attit_qual) Then

c
c  Put the coadds into the coadd buffer.
c
	      num = num + 1
	      coadd_recs(num) = coadd_rec

	   EndIf

	EndDo	!  While (ct_stat


C
C  Check the status of the read.
C

	If (ct_stat(1) .Eq. ctp_normal) Then
	   status = %Loc(ffl_normal)
	ElseIf (ct_stat(1) .Eq. ctp_endoffile) Then
	   status = %Loc(ffl_eof)
	Else
	   status = %Loc(ffl_ctifgread)
	   Call LIB$Signal (ffl_ctifgread, %Val(2), ffli.infile1(1:ffli.inlen1),
     &					   %Val(ct_stat(1)))
	EndIf


	If (status .Eq. %Loc(ffl_eof)) Then
C
C  Close the archive after the final read.
C
	   Call ct_close_arcv (, ct_lun1, ct_stat)
	   If (ct_stat(1) .Ne. ctp_normal) Then
	      status = %Loc(ffl_ctifgclose)
	      Call LIB$Signal (ffl_ctifgclose, %Val(2),
     &			       ffli.infile1(1:ffli.inlen1), %Val(cstatus))
	   EndIf

c
c  Free the Cobetrieve logical unit number.
c
	   rstatus = fut_free_lun (ct_lun1)
	   If (rstatus .Ne. %Loc(fut_normal)) Then
	      Call LIB$Signal (%Val(rstatus))
	   EndIf

	EndIf		!  (end of file


	FFL_Read_Coadd = status

	Return
	End
