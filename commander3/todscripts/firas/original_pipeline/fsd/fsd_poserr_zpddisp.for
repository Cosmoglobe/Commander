	Subroutine FSD_Poserr_ZPDDisp(ctr, ctr_int)

C-----------------------------------------------------------------------------
C	This subroutine displays zero path difference(ZPD) position of ifg
C	and intensity at zpd position
C
C-----------------------------------------------------------------------------
C   Changes:
C
C	SPR 4178, Allow PosErr to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_PosErr, _Select, _Spectrum,
C	    _ZPDDisp, _PosDisp, and FSD.CLD.  Fred Shuman, STX / 1989 Sep 8.
C
C       Version 4.4.1 9/25/89 SER 4568, R. Kummerer, STX
C		Use PLT instead of TEMPLATE for graphics.
C
C	SER 5728, Steven Alexander, STX, March 7, 1990.  Add capability
C		to pass in a PLT command file.
C
C-----------------------------------------------------------------------------

	Implicit None

	Include '(FUT_Params)'
	Include '(FSD_PosErr)'

	Real        * 4       ctr(100), ctr_int(100)

	Real        * 4       y(100)
	Logical     * 4       z_flag
	Character   *16       xlabel
	Character   *16       ylabel
	Integer     *4        i

	Write (label,50) gmt(1), gmt(num)
50	Format (1x, a14, ' - ', a14)

	z_flag = .True.
	Do While (z_flag)

	   Write (6,5)
5	   Format (/, ' Display ZPD data:')
	   Type 10
10	   Format(' Enter P to display ZPD position', /,
	2         '       I to display ZPD intensity', /,
	3         '       Q to quit plotting ZPD values', /,
	4                  10x, '= > ', $)
	   Accept 1, ans
1	   Format(a)
	   Call STR$UpCase(ans, ans)
	   xlabel = 'Number'

	   If ((ans(1:1) .Eq. 'P') .Or. (ans(1:1) .Eq. 'I')) Then
	      If (ans(1:1) .Eq. 'P') Then
	         ylabel = 'Position'
		 plot_label = 'ZPD Position'
	         Call lib$movc3(num*4, ctr(1), y(1))
	      Else If (ans(1:1) .Eq. 'I') Then
	         ylabel = 'Intensity'
		 plot_label = 'ZPD Intensity'
	         Call lib$movc3(num*4, ctr_int(1), y(1))
	      End If

	      startx = -1
	      spacex = 1

	      Do i=1,num
		 difg(i) = y(i)
	      Enddo

	      Call FUT_Plot_Title(label, ngroup, nsweep, speed,
	2			  length, chan_id, upmode, xcal_pos,
	3			  fakeit, gain, plot_label, ptitle)
	      Call FUT_Plot(difg, num, startx, spacex, ptitle,
	2		    xlabel, ylabel, zp,
	3		    interactive, plot_device, plt_com, plt_com_file)
	   Else
	      If (ans .Eq. 'Q') Then
	         z_flag = .False.
	      End If
	   End If      ! ((ans(1:1) .Eq. 'P') .Or. (ans(1:1) .Eq. 'I'))
	End Do         ! (z_flag)

	Return
	End
