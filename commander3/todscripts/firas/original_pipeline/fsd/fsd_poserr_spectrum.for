	SUBROUTINE FSD_POSERR_SPECTRUM(FGS,CTR,CTR_INT,STN,CTN)

C-----------------------------------------------------------------------------
C
C ***    TEST POSITION ERROR IN SIN,COS SPECTRUM
C
C	INPUT:
C		FGS	IFGS(512,NUM)
C
C      OUTPUT:
C	       STN     IMAGINARY  SPECTRUM(513,NUM)
C	       CTN     REAL SPECTRUM (513,NUM)
C-----------------------------------------------------------------------------
C
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

	include '(fut_params)'
	include '(fsd_poserr)'

	real*4 f(512),ff(1024),fgs(512,100),y(513),xout(513),yout(513)
	real*4 ctr(100),ctr_int(100)
	real*4 stn(513,100),ctn(513,100),st(513),ct(513),wk(10)
	integer iwk(10)
	logical d_flag
	character reply,xlabel*16,ylabel*16

1	format(a)
	fact=sqrt(512.)

	Write (6,5)
5	Format (/, ' Computing sine and cosine components of spectra.')

	do j=1,num

	   call lib$movc3(512*4,fgs(1,j),f)
c
c ***	zeropad data for fft
c
	   ictr=(ctr(j)+0.5)

	   do i=ictr,1024
    	      if (i .gt. 512)then
	         ff(i-ictr+1)=0.0
	      else	
	         ff(i-ictr+1)=f(i)
	      endif
 	   enddo	

	   do i=1,ictr-1
	      ff(1024-ictr+1+i)=f(i)
	   enddo
c
c ***   Compute fourier transform
c
	   call fftrf (1024,ff,ff)
c
c ***	Normalize the spectrum
c
	   st(1)=ff(1)/fact
	   ct(1)=0.0

	   do i=2,1022,2
	      k=i/2+1
	      st(k)=ff(i)/fact
	      ct(k)=-ff(i+1)/fact
	   enddo

	   st(513)=ff(1024)/fact
	   ct(513)=0.0

	   call lib$movc3(513*4,st,stn(1,j))
	   call lib$movc3(513*4,ct,ctn(1,j))
	enddo

	Write (6,25)
25	Format (/, ' Display spectrum:')

	reply='I'
	do while (reply .ne. 'Q')
	   d_flag=.false.
	   type 30
30	   format(' Enter R to display REAL part of spectrum',/,
     .            '       I to display IMAGINARY part of spectrum',/,
     .            '       P to display PHASE',/,
     .            '       Q to quit',/,
     .                    10x,' > ',$)

	   accept 1,reply
	   call str$upcase(reply,reply)

	   if (reply .ne. 'Q')then
	      xlabel = 'Frequency'
	      ylabel = 'Counts'

	      if ((reply .eq. 'I') .or.(reply .eq. 'R') .or.
     .            (reply .eq. 'P')) d_flag=.true.

	      do while (d_flag)
	         ans='Y'
	         Write (6,35)
35		 Format (' Enter serial number of IFG (0 to quit) > ', $)
	         accept *,j

	         if ((j.gt.0).and. (j.le.num))then
		    encode (60,50,label)gmt(j),ctr(j),ctr_int(j),j
50                  format(1x, a14, ' ', 'IFG Ctr=', f7.3, ' ',
     .                     'Peak Ht=', f8.1, ' ', 'IFG Num=', i3)
	            if (reply .eq. 'I')then
	               call lib$movc3(513*4,stn(1,j),y)
		       plot_label = 'Sine Spectrum'
	            elseif (reply .eq. 'R')then
	               call lib$movc3(513*4,ctn(1,j),y)
		       plot_label = 'Cosine Spectrum'
	            elseif (reply .eq. 'P')then
                       call lib$movc5(0,,0,2048,y)
	               call lib$movc3(513*4,stn(1,j),st)
	               call lib$movc3(513*4,ctn(1,j),ct)
	               do i=1,513
                          if ((ct(i) .eq. 0.0).and.(st(i).eq. 0.0))then
	                     y(i)=0.0
    	                  else
  	                     y(i)=atan2(st(i),ct(i))
	                  endif 
	               enddo
		       plot_label = 'Phase'
	               ylabel = 'Degrees'
	            endif  

		    startx = -1
		    spacex = 1

		    Do i=1,513
		      difg(i) = y(i)
		    Enddo

		    Call FUT_Plot_Title(label, ngroup, nsweep, speed,
	2				length, chan_id, upmode, xcal_pos,
	3				fakeit, gain, plot_label, ptitle)
		    Call FUT_Plot(difg, 513, startx, spacex, ptitle,
	2			  xlabel, ylabel, zp,
	3			  interactive, plot_device,
	4			  plt_com, plt_com_file)

            	 elseif (j .eq. 0)then
	              d_flag=.false.
	         else
	              print *,'Improper choise of spectum number. Try again.'

	         endif !J

 	      enddo  !d_flag

            endif !reply
	enddo !reply

	return
	end
