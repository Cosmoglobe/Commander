	Integer * 4 Function  ffl_bracket_coadd (ct_lun1, ct_lun2, coadd_rec1,
     &						 coadd_rec2, bracket, ffli)

c-------------------------------------------------------------------------------
c
c	Function FFL_BRACKET_COADD
c
c	This function reads in coadded IFGs from the Cobetrieve calibration
c	archives.  The function matches coadds from the two input streams, if
c	possible.
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
c		ct_lun2		Integer * 4		coadd file lun
c		ffli		Record Structure defined in FFL_Invoc.txt
c
c	Output:
c		coadd_rec1	coadd record		stream 1 input coadd
c		coadd_rec2	coadd record		stream 2 input coadd
c		bracket		Logical * 1		bracketing coadd flag
c		ffli		Record Structure defined in FFL_Invoc.txt
c
c	Subroutines called:
c		ct_read_arcv
c		time_gt
c		time_lt
c		LIB$Signal
c
c	Include files:
c		ct$library:ctuser.inc
c		ffl_invoc.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Name of facility changed from FFI to FFL.  Fred Shuman, HSTX
c	1995 May 19.
c
c-------------------------------------------------------------------------------

	Implicit None

	Include 'ct$library:ctuser.inc'
	Include '(ffl_invoc)'
c
c  Call arguments:
c
	Integer * 4	ct_lun1		!  CT logical unit number
	Integer * 4	ct_lun2		!  CT logical unit number
	Dictionary 'fil_sky'
	Record /fil_sky/ coadd_rec1
	Record /fil_sky/ coadd_rec2
	Logical * 1	bracket
c
c  All other variables and functions:
c
	Record /fil_sky/ buffer
	Integer * 4	j		!  a counter
	Integer * 4	status		!  return status
	Integer * 2	ct_stat1(20)	!  CT return status
	Integer * 2	ct_stat2(20)	!  CT return status
	Integer * 4	time1a(2)	!  first ifg time of coadd 1
	Integer * 4	time1b(2)	!  last ifg time of coadd 1
	Integer * 4	time2(2)	!  coadd time of coadd 2

	Logical * 1	stream /.True./

	Logical * 1	time_gt
	Logical * 1	time_lt

	External	ffl_ctifgread
	External	ffl_eof
	External	ffl_normal

C
C  Read the coadd records.
C

	ct_stat1(1) = ctp_normal
	ct_stat2(1) = ctp_normal
	status = %Loc(ffl_normal)
c
c  Read the records from the two input streams
c
	Call ct_read_arcv (, ct_lun1, coadd_rec1, ct_stat1)
	If (stream .Eq. .True.) Then
	   If (bracket .Eq. .False.) Then
	      bracket = .True.
	      coadd_rec2 = buffer
	   Else
	      Call ct_read_arcv (, ct_lun2, coadd_rec2, ct_stat2)
	   EndIf

	   If ((ct_stat1(1) .Eq. ctp_normal)  .And.
     &	       (ct_stat2(1) .Eq. ctp_normal)) Then
c
c  Check timetags for bracketing records.
c
	      Do j=1,2
	         time1a(j) = coadd_rec1.coad_spec_head.first_time(j)
	         time1b(j) = coadd_rec1.coad_spec_head.last_time(j)
	         time2(j)  = coadd_rec2.ct_head.time(j)
	      EndDo

	      If (time_lt(time2,time1a)) Then
	         Do While ((time_lt(time2,time1a))  .And.
     &			   (ct_stat2(1) .Eq. ctp_normal))
	            Call ct_read_arcv (, ct_lun2, coadd_rec2, ct_stat2)
	            Do j=1,2
	               time2(j)  = coadd_rec2.ct_head.time(j)
	            EndDo
	         EndDo
	      EndIf

	      If (ct_stat2(1) .Eq. ctp_normal) Then
	         If (time_gt(time2,time1b)) Then
	            bracket = .False.
	            buffer = coadd_rec2
	         Else
	            bracket = .True.
	         EndIf
	      EndIf
	   EndIf	!  (ct_stat from initial read

	EndIf		!  (stream


C
C  Check the status of the read.
C

	If (ct_stat1(1) .Eq. ctp_normal) Then
	   status = %Loc(ffl_normal)
	ElseIf (ct_stat1(1) .Eq. ctp_endoffile) Then
	   status = %Loc(ffl_eof)
	Else
	   status = %Loc(ffl_ctifgread)
	   Call LIB$Signal (ffl_ctifgread, %Val(2), ffli.infile1(1:ffli.inlen1),
     &					   %Val(ct_stat1(1)))
	EndIf

	If (stream .Eq. .True.) Then
	   If (ct_stat2(1) .Eq. ctp_normal) Then
	      stream = .True.
	   ElseIf (ct_stat2(1) .Eq. ctp_endoffile) Then
	      stream = .False.
	      bracket = .False.
	   Else
	      status = %Loc(ffl_ctifgread)
	      Call LIB$Signal (ffl_ctifgread, %Val(2),
     &			       ffli.infile2(1:ffli.inlen2), %Val(ct_stat2(1)))
	   EndIf
	EndIf


	ffl_bracket_coadd = status

	Return
	End
