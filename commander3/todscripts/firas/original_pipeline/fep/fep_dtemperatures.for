	subroutine fep_dtemperatures (nrec, timetags,values, dwell_loc, jside,
	2                             config_lun,ngrt, grtnums, icurrent, ier)

C------------------------------------------------------------------------------
C
C    Converts the array 'values' (a time-series of resistances) to Kelvin and
C	returns results in the same array.  Any portion of the time series that
C	corresponds to a calibration resistor enters this routine in counts and
C	is unaltered, exiting in counts.
C
C------------------------------------------------------------------------------
C	Change Log:
C	  Shirley M. Read, STX, September 1988.
C		Added the pause after the type statement so that the user
C	        could read the message and then type a continue character.
C
C	  R. Kummerer, January 12, 1989.
C		SPR 3132, Dwellplot fails to plot XCAL S5 and S6.
C
C	  Fred Shuman, STX,  1989 Feb 17.
C		SPR 2910, Plot all fields dwelled in the timerange rather than
C		just the one with the plurality of hits on each side (A/B).
C
C	  Fred Shuman, STX,  1989 Feb 20.
C		SPR 2738, Add the ability to plot cal resistors.  Plot them in
C		counts.
C
C	  Harte Wang, STX,  1989 Dec. 8.
C		SPR 5130, Spurious temperature jump in plots.
C
C	  Harte Wang, STX,  1990 Feb. 22.
C		rename fep_grt_correction to Grt_correction
C
C          Spr 3253, Dwellplot has to fetch the reference data from the 
C                    reference data archive.
C                    H. Wang, STX, Mar. 2, 1990
C	SPR 8956, removal of grt_correction routine.
C		L. Rosen, STX, 6 Sept. 1991
C------------------------------------------------------------------------------
	Implicit None

	integer*2 nrec, dwell_loc(1)
        Integer*4 config_lun, index, Btime(2)
	integer*4 jside, ier, dt, grtold, ngrt, grtnums(1), icurrent
	integer*4 status, irec, minor, npts, k, iwhichgrt
	real*4 values(1), ohms_tbl(400), temps_tbl(400), grttemp
        Character*14 timetags(1000)
	character*40 grttable
	character*72 dbstring
	logical*1 foundit, newgrt, calres, Signaled,first,new_grt  
	character*1 ans
	Integer*4 Lib$Signal, FEP_GRT_Lookup
        Integer*4 Cct_get_config_idx_tod
        Integer*4 Fep_get_grt_conv
        Include '(Fep_dwell_dbword)'
        Include '(cct_get_config)'
        Record /config_status/stat

	External FEP_TBLEDGE
	External FEP_Getconfigerr

	Signaled = .False.
	ngrt = 0
	grtold = -64
        first  = .true. 
        Ier = 0
	Do irec=1,nrec
c
c  First figure out which GRT we are dwelling on
c
	   icurrent = dwell_loc(irec)/16		! 0=lo current, 1=hi
	   iwhichgrt = dwell_loc(irec) + 1 - 16*icurrent
	   iwhichgrt = iwhichgrt + 16*jside		! 1-16=A, 17-32=B

	   If (iwhichgrt .Eq. grtold) Then
	      newgrt = .False.
	   Else
	      newgrt = .True.
	      grtold = iwhichgrt
	      ngrt = ngrt + 1
	      grtnums(ngrt) = iwhichgrt
	      If ((iwhichgrt .Ge. 11  .And.  iwhichgrt .Le. 14) .Or.
	2         (iwhichgrt .Ge. 27  .And.  iwhichgrt .Le. 30)
	3	        .or. (iwhichgrt .eq. 32)) Then
	         calres = .True.
	         type *, ' This Dwell Field is a CAL RESISTOR and cannot be' //
	2                ' temperature-converted.'
	         type *, '    It will be plotted in counts.'
	      Else
	         calres = .False.
	      End If
	   End If
c
c  Now open the appropriate GRT lookup table and read it
c
	   If (.Not. calres) Then
              Call ct_gmt_to_binary(timetags(nrec),Btime)
              status = cct_get_config_idx_tod(btime,1,config_lun,index,
	1	        new_grt,stat)
              If (.not. status) then
                call lib$signal(Fep_GetConfigErr,%val(1),%val(status))
                Ier = 1
                Return
              Endif
C                
 	      If (newgrt .or. new_grt) Then
                 Status = fep_get_grt_conv(config_lun,dbword(iwhichgrt),
	1		          Ohms_tbl,temps_tbl,npts)
              If (.not. status) then
                Ier = 1
                Return
              Endif
  	    End If
c
c  Do the conversion, point by point. 
c
	      Do minor=1,64
	         k = (irec-1)*64 + minor
	         status = FEP_GRT_Lookup(values(k), 0, ohms_tbl, temps_tbl,
	2                                first,npts, grttemp, dt)
                 first = .false.
	         If (status .Eq. %Loc(FEP_TblEdge)) Then
	            If (.Not. Signaled) Then
	               Call Lib$Signal(FEP_TblEdge)
	               Signaled = .True.
	            End If
	            values(k) = 0.
	         Else
	            values(k) = grttemp
	         End If
	      End Do       ! minor=1,64
	   End If        ! not a cal resistor
	End Do         ! irec=1,nrec
	Return
	End
