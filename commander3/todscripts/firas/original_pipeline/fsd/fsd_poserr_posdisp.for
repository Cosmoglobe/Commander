	Subroutine FSD_PosErr_PosDisp(mean, sd, dmean, dsd, m)

C-----------------------------------------------------------------------------
C
C ***	DISPLAY ROUTINE
C
C	INPUT:
C		MEAN    MEAN VALUE
C		SD      RMSD
C		DMEAN   MEAN OF DERIVATIVE
C	        DSD     RMSD OF DERIVATIVE
C		M       # OF POINTS
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

	Integer   * 4    m
	Real      * 4    mean(m), sd(m), dmean(m), dsd(m)

	Integer   * 4    i, ii
	Integer   * 4    incx
	Integer   * 4    imx
	Integer   * 4    np
	Integer   * 4    npt
	Integer   * 4    isamax
	Integer   * 4    numpts
	Integer   * 4    status
	Real      * 4    y(513), dummy(513)
	Real      * 4    hist(101), max_val, con
	Character *60    xlabel
	Character *60    ylabel
	Integer   * 2    iplot
	Logical   * 4    p_flag

1	Format(a)

	Write (label,50) gmt(1), gmt(num)
50	Format (1x, a14, ' - ', a14)

	p_flag = .True.
	Do While (p_flag)
	   Write (6,5)
5	   Format (/, ' Check coadded data:')
	   Type 10
10	   Format(' Enter 1 Mean vs Sample Number', /,
	2         '       2 RMSD vs Sample Number', /,
	3         '       3 Histogram of RMSD', /,
	4         '       4 RMSD vs Intensity', /,
	5         '       5 RMSD vs Derivative of Intensity', /,
	6         '       6 Quit', /,
	7                 10x, '> ', $)

	   Accept *,iplot
C
C Display mean vs position.
C
	   If (iplot .Eq. 1) Then
	      xlabel = 'Sample Number'
	      ylabel = 'Intensity'
	      plot_label = 'Mean vs Sample Number'
	      startx = -1
	      spacex = 1
	      npt = m
	      Call LIB$MovC3(m*4, mean(1), y(1))
C
C Display RMSD vs position.
C
	   Else If (iplot .Eq. 2) Then
	      xlabel = 'Sample Number'
	      ylabel = 'RMSD'
	      plot_label = 'RMSD vs Sample Number'
	      startx = -1
	      spacex = 1
	      npt = m
	      Call LIB$MovC3(m*4, sd(1), y(1))
C
C Display histogram of RMSD.
C
	   Else If (iplot .Eq. 3) Then
	      xlabel = 'RMSD'
	      ylabel = 'Number'
	      plot_label = 'Histogram of RMSD'
	      Call LIB$MovC5(0, ,0, 101*4, hist)
C
C Compute the histogram; duplicate consecutive points to get
C appearance of histogram.
C
	      Call LIB$MovC3(m*4, sd(1), y(1))

	      np = m - 1
	      incx = 1
	      imx = isamax(np, y(1), incx)
	      max_val = abs(y(imx))

	      If (max_val .Ne. 0.0) Then
	         con = 100./max_val
	      Else
	         con = 1
	      End If

	      Do i=1,m
	         ii = jint(y(i)*con)+1
	         hist(ii) = hist(ii)+1.
	      End Do

	      startx = -1/(2*con)
	      spacex =  1/(2*con)
	      npt = 202

	      Do i=1,m
		ii = (i-1)*2 + 1
		y(ii) = hist(i)
		y(ii+1) = hist(i)
	      End Do
C
C Display mean vs RMSD.
C
	   Else If  (iplot .Eq. 4) Then
	      xlabel = 'Intensity'
	      ylabel = 'RMSD'
	      plot_label = 'RMSD vs Intensity'
	      startx = -1
	      spacex =  1
	      npt = m
	      Call LIB$MovC3(m*4, mean(1), dummy(1))
	      Call LIB$MovC3(m*4, sd(1), y(1))
C
C Display mean vs RMSD of the derivative.
C
	   Else If (iplot .Eq. 5) Then
	      xlabel = 'Derivative of Intensity'
	      ylabel = 'RMSD'
	      plot_label = 'RMSD vs Derivative of Intensity'
	      startx = -1
	      spacex =  1
	      npt = m
	      Call LIB$MovC3(m*4, dmean(1), dummy(1))
	      Call LIB$MovC3(m*4, dsd(1), y(1))
	   Else
	      p_flag = .False.
	   End If

C
C Make the plot.
C
	   If (p_flag) Then
	      Call FUT_Plot_Title(label, ngroup, nsweep, speed,
	2			  length, chan_id, upmode, xcal_pos,
	3			  fakeit, gain, plot_label, ptitle)
	      If (iplot .Ne. 4 .And. iplot .Ne. 5) Then
	         Do i=1,npt
	           difg(i) = y(i)
	         Enddo
	         Call FUT_Plot(difg, npt, startx, spacex, ptitle,
	2		       xlabel, ylabel, zp,
	3		       interactive, plot_device,
	4		       plt_com, plt_com_file)
	      Else
	         Call FUT_Scatter_Plot(dummy, y, npt, ptitle,
	2			       xlabel, ylabel, zp,
	3			       interactive, plot_device,
	4			       plt_com, plt_com_file)
	      End If
	   End If

	End Do

	Return
	End
