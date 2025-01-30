	Program FSD_Poserr

C-----------------------------------------------------------------------------
C
C ***   PROGRAM FOR STUDYING THE POSITION ERROR IN DATA *******
C ***   ACCEPTS  DATA FROM RMS FILE(MINIPIPE OUTPUT) OR ARCHIVE FILE
C
C ***   J.H.DAVE'  ARC, NOV 86
C
C	COMPUTE MEAN,RMSD OF IFGS, DERIVATIVE OF IFGS
C		 "   "    OF SIN FFT, DERIVATIVE OF SIN FFT
C		 "   "    OF COS FFT, DERIVATIVE OF COS FFT
C-----------------------------------------------------------------------------
C   Changes:
C
C	9/4/87   INCLUDED PLOT OF  ZERO PATH DIFFERENCE POSITION
C                AND INTENSITY AT ZPD OF INTERFEROGRAM.
C	9/10/87  USE TEMPLATE ROUTINES INSTEAD OF NEWVECPLT FOR PLOTS
C
C ***   Shirley M. Read  STX, July 1988
C
C	Added capability of smoothing IFGs and making a template IFG for
C	comparisons. The generated plots contain the template IFG, the
C	difference between the individual smoothed IFG and the template IFG,
C	and the derivative of the template IFG.
C
C       version 4.2.1 11/30/88, SPR 2310, QUOC CHUNG, STX
C       BRING FSD UP TO STANDARD ERROR MESSAGE TRAPPING,
C       AND EXIT STATUS.
C
C       version 4.2.1 2/10/89 SPR 3135, R. KUMMERER, STX
C	DITHER NOW REMOVED IN FSD_POSERR_SELECT PRIOR TO FUT_FIND_IFG_CTR.
C	REMOVE DUPLICATE DITHER REMOVAL.
C
C	VERSION 4.4.1 SER 3306 Q. CHUNG STX 08/10/89
C                    PROVIDE VERSION NUMBER TO TRACK SOFTWARE UPDATE.
C
C	SPR 4178, Allow PosErr to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_PosErr, _Select, _Spectrum,
C	    _ZPDDisp, _PosDisp, and FSD.CLD.  Fred Shuman, STX / 1989 Sep 8.
C
C       Version 4.4.1 9/25/89 SER 4568, R. Kummerer, STX
C		Use PLT instead of TEMPLATE for graphics.
C
C       Version 4.4.1 9/27/89 SPR 4635, R. Kummerer, STX
C		Flag process run as aborted only when an error actually
C		occurs.
C
C	SER 5728, Steven Alexander, STX, March 7, 1990.  Add capability
C		to pass in a PLT command file.
C
C-----------------------------------------------------------------------------

	Implicit None

	Integer  * 4  CUT_Register_Version
	Integer  * 4  CUT_Display_Banner
	Integer  * 4  cli$present
	Integer  * 4  cli$get_value
	Integer  * 4  num_vol/80/
	Integer  * 4  lun_out/6/
	Character* 6  version
	Parameter    (version = '5.8')

C   Include Files.

	Include '(FUT_Params)'
	Include '(FSD_Poserr)'
	Include '($SSDef)'

C   Functions Invoked.

	Integer  * 4  FSD_Poserr_Smooth_IFGs  !Smooths IFGs with boxcar
	                                      !function as low pass filter
	Integer  * 4  FUT_Swg_Make_Template   !Creates the template IFG from
	                                      !the individual IFGs.
	Integer  * 4  FSD_Poserr_Compare_Disp !Computes the differences of the
	                                      !IFGs and the template. Computes
	                                      !the derivatives at each point
	                                      !of the template. Displays the
	                                      !template, the differences and
	                                      !the derivatives.
C   Declarations.

	Real     * 4  fgs(512,100), stn(513,100), ctn(513,100)
	Real     * 4  smooth_ifgs(512,100)    !Smoothed IFGs
	Real     * 4  dummy_ifgs(512,100)     !Dummy IFGs for function arg
	Real     * 4  template_ifg(512)       !Template IFG from smoothed IFGs
	Real     * 4  sigmas(100)             !RMS noise from IFGs - template
	Real     * 4  ctr(100), ctr_int(100)  !IFG center pk positions & values
	Real     * 4  buffer(512)             !Dummy IFG buffer for sort
	Real     * 4  avg_sigma           !Median of elements in array "Sigmas"
	Character* 1  reply, term, spec
	Integer  * 4  np
	Logical  * 4  i_spec, valnum
	Real     * 4  mean(513), sd(513), dsd(513), dmean(513)
	Integer  * 4  status
	Integer  * 4  i, j                    !Indices

	External   FSD_Normal
	External   FSD_AbErr2
	External   FSD_ErrSIFG
	External   FSD_ErrCIFG
	External   cli$_present

1	Format(a)

	status = CUT_Register_Version(version)
	status = CUT_Display_Banner(lun_out, num_vol,
	1                           'FIRAS Facility FSD_Poserr')

	init = .True.
	status = %Loc(FSD_Normal)

	interactive = fac_present
	zp = fac_not_present
	plot_device = 'PLT_DEVICE'

	plt_com = fac_not_present
	if (cli$present('PLTFILE') .eq. %loc(cli$_present)) then
	   plt_com = fac_present
	   status = cli$get_value('PLTFILE',plt_com_file)
	end if

	Do i=1,513
	   x(i) = floatj(i)
	End Do

	reply = 'Y'
	Do While (reply .Eq. 'Y')
C
C ***   READ (M*N) IFGS
C
	   Call FSD_Poserr_Select(fgs, ctr, ctr_int)
	   If (num .Eq. 1) Then
	      Print *, 'Number of IFGs = 1, No position error test.'
	   Else If (num .Ge. 2) Then
C
C **    PROCESS IFG OR SPECTRUM
C
	      i_spec = .False.
	      term = 'I'

	      Do While (term .Ne. 'Q')
	         Write (6,5)
5		 Format (/, ' Position error test:')
	         Type 10
10	         Format(' Enter I to test IFG data ', /,
	2               '       S to test SPECTRUM data ', /,
	3               '       C to smooth IFGs and compare to template', /,
	4               '       Q to Quit current channel ', /,
	5                       10x, ' > ', $)
	         Accept 1,term
	         Call STR$UpCase(term,term)

	         If (term .Eq. 'I') Then
	            reply(1:1) = 'Y'
	            Do While (reply(1:1) .Ne. 'Q')
	               Write (6,15)
15		       Format (/, ' Testing for position error in IFG:')
	               Type 20
20	               Format(' Enter Z to check ZPD position ', /,
	2                     '       C to check coadded IFG data ', /,
	3                     '       Q to quit processing IFG data ', /,
	4                             10x, ' > ', $)
	               Accept 1,reply
	               Call STR$UpCase(reply,reply)
	               If (reply(1:1) .Eq. 'Z') Then
	                  Call FSD_Poserr_Zpddisp(ctr, ctr_int)
	               Else If (reply(1:1) .Eq. 'C') Then
C
C ***   MEAN AND VARIANCE OF IFGS
C
	                  Call FSD_Poserr_Pos(fgs, 512, mean, sd, dmean, dsd)

	                  Call FSD_Poserr_Posdisp(mean, sd, dmean, dsd, 512)
	               End If
	            End Do

	         Else If (term .Eq. 'S') Then

	            If (.Not. i_spec) Then
	               Call FSD_Poserr_Spectrum(fgs, ctr, ctr_int, stn, ctn)
	               i_spec = .True.
	            End If

	            Write (6,25)
25		    Format (/, ' Testing for position error in spectra:')

	            spec = 'R'
	            Do While (spec .Ne. 'Q')
	               Type 30
30	               Format(' Enter R for REAL part of spectrum ', /,
	2                     '       I for IMAGINARY part of spectrum ', /,
	3                     '       Q to quit ', /,
	4                             10x, ' > ', $)

	               Accept 1,spec
	               Call STR$UpCase(spec,spec)

	               If (spec .Eq. 'I') Then

	                  Call FSD_Poserr_Pos(stn, 513, mean, sd, dmean, dsd)
	                  Call FSD_Poserr_Posdisp(mean, sd, dmean, dsd, 513)

	               Else If (spec .Eq. 'R') Then

	                  Call FSD_Poserr_Pos(ctn, 513, mean, sd, dmean, dsd)
	                  Call FSD_Poserr_Posdisp(mean, sd, dmean, dsd, 513)

	               End If

	            End Do

	         Else If (term .Eq. 'C') Then

	            valnum = .False.
	            Do While (.Not. valnum )
	               Type 40
40                     Format(/, 1x, 'Enter N, the number of points to average on each side'/
	2                        1x, '         in the boxcar low pass filter'/
	3                        1x, '         i - N <-- i --> i + N   '/
	4                        10x, ' > ', $)
	               Accept *,np
	               If ((np .Lt. 0) .Or. (np .Gt. 255)) Then
	                  Print *, 'Invalid choice. Please enter 0 to 255.'
	               Else
	                  valnum = .True.
	               End If
	            End Do

	            Write (6,45)
45		    Format (/, ' Smoothing IFGS and comparing to template.')

	            status = FSD_Poserr_Smooth_IFGs ( fgs, np, smooth_ifgs )
	            If (status .Ne. %Loc(FSD_Normal)) Then
	               term = 'Q'
	               Call LIB$Signal(FSD_Errsifg, %Val(1), %Val(num))
	            Else

C   Use dummy argument for smoothed IFGs. Function changes input values.

	               Do i=1,num
	                  Do j=1,512
	                     dummy_ifgs(j,i) = smooth_ifgs(j,i)
	                  End Do
	               End Do

	               status = FUT_Swg_Make_Template ( dummy_ifgs, num,
	2                              template_ifg, sigmas, avg_sigma )
	               Write (6,60) avg_sigma
60		       Format (' Median of sigmas for IFGs =', f9.3, '.', /)

C  Invoke function to compare and display the IFGs.

	               status = FSD_Poserr_Compare_Disp (smooth_ifgs,
	2                                                template_ifg, sigmas )
	               If (status .Ne. %Loc(FSD_Normal)) Then
	                  term = 'Q'
	                  Call LIB$Signal(FSD_Errcifg, %Val(1), %Val(status))
	               End If
	            End If

	         Else

	            If (term .Ne. 'Q') Then
	               Print *, 'Improper Choice.  Try again.'
	            End If

	         End If
	      End Do

	   Else    !r_data

	      Print *, 'No data, quit this channel.'

	   End If  !r_data

	   Write (6,55)
55	   Format (/, ' Do you want to examine another channel? (Y/[N]) >', $)
	   Accept 1,reply
	   Call STR$UpCase(reply,reply)

	End Do   !reply

	Call LIB$Erase_Page ( 1, 1 )

C
C **** exit with proper message
C
	If (status) Then
	   Call LIB$Signal(FSD_Normal)
	   Call Exit (SS$_Normal)
	Else
	   Call LIB$Signal(FSD_Aberr2)
	   Call Exit (SS$_Abort)
	End If

	Stop
	End
