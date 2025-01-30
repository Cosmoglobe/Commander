	Integer*4  Function  FAD_Close_Archives  ( lun_in, lun_out, lun_rpt,
	1                                          report, filein, fileout )

c------------------------------------------------------------------------------
c   Purpose: Close the input, reference, and output files.
c
c   Input Parameters:
c      integer*4     lun_in        -- logical unit number for input file
c      integer*4     lun_out       -- logical unit number for output file
c      integer*4     lun_rpt       -- logical unit number for report file
c      logical*1     report        -- true if report is to be written (default)
c      character*50  filein, fileout -- file names
c
c   Output Parameters:
c      none
c
c   Include Files:
c      FAD_msg       -- External message params needed for FAD
c      FUT_params    -- FIRAS parameters
c
c   Functions:
c      Lib$Signal
c      CSA_Close_Skymap
c
c   Author:  Larry Paris Rosen, Hughes STX, 20 April 1993
c   Modified: L. Rosen, February 1994.  Add FEX_GN file.
c   Modified: L. Rosen, May 1994.  No gain file anymore.
c------------------------------------------------------------------------------
	Implicit None

c Passed Parameters

	Integer*4	lun_in, lun_out
	Integer*4	lun_rpt
	Logical*1	report
	Character*50	filein, fileout 			! file names

c Include

	Include		'(fad_msg)'
	Include		'(fut_params)'

c Functions

	Integer*4	CSA_Close_Skymap

c External

	External	csa_normal

c Local

	Integer*4	rstat

c------------------------------------------------------------------------------
c Begin

	FAD_Close_Archives = %LOC (FAD_Normal)

	rstat = CSA_Close_Skymap (lun_in, fac_skymap_no_levels)
	If (rstat .NE. %Loc (CSA_Normal)) Then
	   FAD_Close_Archives = %Loc (FAD_Abort)
	   Call Lib$Signal (fad_csaclose, %Val(2), filein, %Val(rstat))
	Else
	   If (report) Then
	      Write (lun_rpt,*)
	      Write (lun_rpt, 10) filein
  10	      Format (1X, 'Successfully closed: ',A)
	   Endif
	   rstat = CSA_Close_Skymap (lun_out, fac_skymap_no_levels)
	   If (rstat .NE. %Loc (CSA_Normal)) Then
	      FAD_Close_Archives = %Loc (FAD_Abort)
	      Call Lib$Signal (fad_csaclose, %Val(2), fileout, %Val(rstat))
	   Else
	      If (report) Then
	         Write (lun_rpt, 10) fileout
	      Endif
	   Endif
	Endif
	Return
	End
