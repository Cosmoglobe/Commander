	subroutine fsd_glitchmap_lin_disp(ans)

C--------------------------------------------------------------------
C
C       VERSION 4.2.1 SPR 3057, QUOC CHUNG STX, 12/29/88
C       Fixed the problem of the variable ANS did not pass back
C       to the main program to set up proper termination of the
C       processing channel loop.
C
C	Version 4.4.1 SER 4567 R. Kummerer STX 09/15/89
C	Use PLT instead of TEMPLATE for graphics.
C
C	Version 4.4.2 SER 5728 S. Alexander STX 03/07/90
C	Add capability to pass in a PLT command file.
C
C---------------------------------------------------------------------

        implicit none

	include '(fsd_glitchmap)'

	integer   * 4      status
	integer   * 4      k
        integer   * 4      incx
        integer   * 4      imx
        integer   * 4      np
	integer   * 4      mi
	integer   * 4      ii
	integer   * 4      idup
	integer   * 4      ik
	integer   * 4      iopt
	integer   * 4      nsum
	integer   * 4      numpts
	
	integer   * 2      ibin

	integer   * 4      ppts
	real      * 4      pstart
	real      * 4      pbin

	complex   * 8      difg(1024)
	character * 60     plot_label
	character *100     title(3)
	integer   * 4      zp
	integer   * 4      interactive

	real      * 4      rifg(512)
	real      * 4      hist(512)
	real      * 4      rout(512)
	real      * 4      dummy1(512)
	real      * 4      dummy2(512)

        real      * 4      avg,sum
        real      * 4      con,cons
        real      * 4      smax,smin
        real      * 4      gmin,gmax

	logical   * 4      i_flag

        character * 1      answer
        character *60      label

	integer   * 4      STR$UpCase
        integer	  * 4      ISMAX
        integer	  * 4      ISMIN

c
c Initialize.
c
	zp = fac_not_present
	if (batch) then
	   interactive = fac_not_present
	else
	   interactive = fac_present
	end if

c
c Prompt for an instruction for display.
c
	ans='P'

	Do While (ans(1:1) .ne. 'Q')

	   Type *
	   Type 40
40	   Format ( ' P to see individual glitch map ', /,
	1	    ' G to see final group glitch map ', /,
	2	    ' Q to quit plotting', /,
	3	    10x,' => ',$)

	   Accept 1, ans
1	   Format(a)

	   status = STR$UpCase(ans,ans)

c
c Plot the data...
c
	   If (ans(1:1) .Eq. 'P')then
	      Type 50, num
50	      Format (' Choose one of ', i2, ' individual glitch maps => ', $)
	      Accept *,mi
	      If (mi .le. num) Then
	         answer = 'Y'
	         avg = avg_g(mi)
	         Call LIB$MovC3(512*4,his(1,mi),dummy1(1))
	         If (ifg_flag) Call LIB$MovC3(512*4,sci_ifg(1,mi),dummy2(1))
	      Else
	         Type *, ' Wrong choice. Try again.'
	         answer = 'N'
	      Endif
	   Elseif (ans(1:1) .Eq. 'G') Then
	      answer = 'Y'
	      avg = 0.0
	      Do i=1,num
	         avg = avg + avg_g(i)
	      Enddo
	      avg = avg/num
	      Call LIB$MovC3(512*4,g_hist,dummy1)
	      If (ifg_flag) Call LIB$MovC3(512*4,g_ifg,dummy2)
	   Else
	      answer = 'N'
	   Endif

	   Do While (answer .Eq. 'Y')

	      Call LIB$MovC5(0,,0,2048,rifg)
	      Call LIB$MovC5(0,,0,2048,rout)
	      Call LIB$MovC5(0,,0,2048,hist)

	      i_flag = .False.

	      Type 60
60	      Format (' Enter bin width => ', $)
	      Accept *, ibin

  	      If (ibin .le.0) ibin = 1
	      If ((ibin .Eq. 1) .and. (ifg_flag)) i_flag = .True.
	      ii = 0
              Do i=1,512,ibin
                 ii = ii + 1
                 sum=0
                 Do j=1,ibin	
	 	    ik = i + j - 1
	            If (ik .le. 512) Then
                       sum = sum + dummy1(ik)
		    Endif
                 Enddo
	         hist(ii) = hist(ii) + sum
	         If (i_flag) rifg(ii) = rifg(ii) + dummy2(ii)
	      Enddo	

	      imx = ISMAX(512,hist,1)
	      gmax = hist(imx)
	      imx = ISMIN(512,hist,1)
	      gmin = hist(imx)

              If (gmax .Eq. 0.) then
            	 iopt = 2
                 Print *, 'No glitches in this region.'
              Else
	         If (i_flag) then
	            iopt = 1
		    imx = ISMAX(512,rifg,1)
		    smax = rifg(imx)
		    imx = ISMIN(512,rifg,1)
		    smin = rifg(imx)
	            con = abs(smin)
	            cons = gmax / (smax + abs(smin))
	            Do i=1,512
	               rout(i) = (rifg(i) + con) * cons
	            Enddo
                 Else
	            iopt = 0
                 Endif
              Endif

	      If (iopt .ne. 2) Then

	         nsum = 0
	         Do i=1,ii
	            nsum = nsum + hist(i)
                 Enddo
c
c Duplicate consecutive data points to make plot look more like a histogram.
c
	         Do i=1,ii
		    idup = (i-1)*2 + 1
	            difg(idup) = Cmplx(hist(i),rout(i))
	            difg(idup+1) = Cmplx(hist(i+1),rout(i+1))
                 Enddo
		 difg(idup+1) = Cmplx(hist(ii),rout(ii))

	 	 ppts = 2*ii
		 pbin = FloatI(ibin)/2.
		 pstart = -pbin

		 If (ans(1:1) .Eq. 'P') Then

	            Encode (60,70,label) gmt(mi)(1:11), ibin, nsum,
	1				 avg, n_over(mi), degl_thre(mi)
70      	    Format(a11,' Bn=',i2,' Gl=',i5,' AvGl=',f5.1,
	1		   ' Ovfl=',z4,' DgTh=',i2)

		    plot_label = 'Raw Science Interferogram Glitch Map'

            	    Call FUT_Plot_Title(label,ngroup,nsweep,mtm_speed,
	1			mtm_length,chan_id,upmode,xcal_pos,fakeit,
	2			gain,plot_label,title)
		    Call FUT_Plot(difg,ppts,pstart,pbin,title,
	1			'Sample Number','Number of Glitches',zp,
	2			interactive,plot_device,plt_com,plt_com_file)

		 Elseif (ans(1:1) .Eq. 'G') Then

          	    Encode (60,80,label) gmt(1)(1:9), gmt(num)(1:9),
	1				 num, ibin, nsum, avg
80		    Format (a9,'-',a9,' IFG',i3,' Bn=',i2,
	1		    ' Gl=',i5,' AvGl=',F5.1)

		    plot_label = 'Coadded Raw Science Interferogram Glitch Map'

            	    Call FUT_Plot_Title(label,ngroup,nsweep,mtm_speed,
	1			mtm_length,chan_id,upmode,xcal_pos,fakeit,
	2			gain,plot_label,title)
		    Call FUT_Plot(difg,ppts,pstart,pbin,title,
	1			'Sample Number','Number of Glitches',zp,
	2			interactive,plot_device,plt_com,plt_com_file)

		 Endif

	      Endif

	      Type 90
90	      Format (' Do want to consider another bin width (Y/[N])? ', $)

	      Accept 1, answer

	      status = STR$UpCase(answer,answer)

	   Enddo

	Enddo

	Return
	End
