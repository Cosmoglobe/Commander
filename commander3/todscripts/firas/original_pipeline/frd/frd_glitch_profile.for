	program frd_glitch_profile
c
c  Use electronics xfer fcn and detector time constant to make predicted 
c  glitch profile.  The output is a direct access file in which there
c  are N records associated with every important value of adds_per_group=N.
c  The organization of the direct-access file is as follows: 
c
c	 MTM 	 channel
c
c              |  RH   (26 records)
c	 Slow  |  RL      "
c              |  LH      "
c              |  LL      "
c
c              |  RH   (26 records)
c	 Fast  |  RL      "
c              |  LH      "
c              |  LL      "
c
c              |  RH   (26 records)
c	 Fakeit|  RL      "
c              |  LH      "
c              |  LL      "
c
c  ...where each group of 26 records consists of a single record for 
c  adds_per_group=1, followed by 2 records for adds=2, and so on for 3, 8, 
c  and 12.  Each record is 512 real*4 values long.  The first 510 points 
c  are the glitch profile, the 511th is the index of the glitch peak 
c  position, and the 512th is the offset from that integer position to
c  the "true" peak position obtained by parabolic interpolation.
c
c  Written by:		Rich Isaacman (286-5758)
c			General Sciences Corp.
c			2 January 1991
c  Modified by John Sims (STX) 11 June 1991, SER 7985
c                       brought up to csdr coding standards
c------------------------------------------------------------------------- 
	implicit none

	include '($ssdef)'
	
	integer*4	status
	integer*4	mtmspd 
	integer*4 	ngroup
	integer*4	npts
	integer*4	ipk_posn
	integer*4	ichan 
	integer*4	iadds
	integer*4	iadds_per_group(5)	/1, 2, 3, 8, 12/
	integer*4	irecno
	integer*4	ishift
	integer*4	k
        integer*4	lun
	integer*4	ios
	integer*4    	nrec
    	character *16	out_file /'fex_gltchpro.dat'/

	complex*8	elex(12288)
	complex*8	twopij /(0.,6.2831857)/

	real*4		df
	real*4		ampl
	real*4		fnyq(3)	/340.71, 340.71, 256./
	real*4		shape(12288)
	real*4		tlmout(512)
	real*4 		tau_det(4) /0.0055, 0.0438, 0.0055, 0.0438/

	real*4		rwksp(78467)
	common	/worksp/ rwksp

	integer *4	lib$get_lun
	integer *4	lib$free_lun

	external	frd_normal
	external	frd_rmsopen
	external	frd_rmswrite
	external	frd_rmsclose

	call iwkin (78467)

	status = lib$get_lun(lun)

	if (status .eq. ss$_normal) then

	open (unit=lun, name=out_file, status='new',
     .			form='unformatted', access='direct',recl=512)

	if (ios .eq. 0) then

	nrec = 0                               
c
c  Generate elex xfer fcn up to 681.4 Hz (raw sample rate) for MTM 
c  scans, or 512 Hz for fakeit data.  Loop over the 3 sample rates 
c  (slow, fast, fakeit), then the channels, then the adds_per_group.
c
	do mtmspd = 1,3
	   do ichan = 1,4
	      do ngroup = 1,5
		 iadds = iadds_per_group(ngroup)
	   	 df = fnyq(mtmspd)/256./iadds
	   	 npts = 512 * iadds
	         call make_xfer_fcn (npts, df, ichan, 4, 
     .						2.*fnyq(mtmspd), elex)
c
c.....  convolve with exponential glitch shape at detector, i.e. 
c.....  multiply by 1/(1+ j*omega*tau) using detector time constant
c
		 do k=1, npts
	   	    elex(k) = elex(k)/(1. + twopij * df * (k-1) * 
     .							tau_det(ichan))
		 enddo
c
c.....  make the transfer function Hermitian, then inverse FFT into time 
c.....  domain
c
		 do k = npts+2, 2*npts
	   	    elex(k) = conjg (elex(2*npts+2-k))
		 enddo
		 elex(npts+1) = (0.,0.)
		 call fftcb (2*npts, elex, elex)
c
c.....  copy it into a real array in preparation for data compression, 
c.....  keeping track of the position and height of the peak.
c
		 ampl = -1.
		 do k=1,2*npts
		    shape(k) = real(elex(k))
		    if (abs (shape(k)) .gt. ampl) then
			ampl = abs (shape(k))
			ipk_posn = k
		    endif
		 enddo
c
c.....  Compress data down to 512 points.  For iadds = N, do it N times, 
c.....  sliding the window over by one sample each time to simulate 
c.....  glitch hits at all possible positions within the averaging 
c.....  window.  The compressed profile is returned in tlmout.
c
		 do ishift = 1, iadds
		    call squeeze512 (ishift, iadds, npts, ipk_posn,
     .	 						shape, tlmout)
		    nrec = nrec + 1
		    write (lun, rec=nrec, iostat=ios) tlmout
			if(ios .ne. 0)then
     			  call lib$signal(frd_rmswrite,%val(2),
     .				out_file,%val(ios))
	                  status = %loc(frd_rmswrite)
			endif
		 enddo
	      enddo
	   enddo
	enddo
	if (status .ne. ss$_normal) call lib$signal(%val(status))

	close (lun, iostat=ios)
	if(ios .ne. 0) then
		status = %loc(frd_rmsclose)
		call lib$signal(frd_rmsclose,%val(2),out_file,%val(ios))

	endif
	
	else		!ios .eq. 0
		call lib$signal(frd_rmsopen,%val(2),out_file,%val(ios))
	endif		!ios .eq. 0

	status = lib$free_lun(lun)

	else		!ss$_normal
		call lib$signal(%val(status))
	endif		!ss$_normal

	if(status .eq. ss$_normal)then
		call lib$signal(frd_normal)
	else
		call lib$signal(%val(ss$_abort))
	endif                                    
	end


	subroutine make_xfer_fcn (npts, dfhz, ichan, micro, 
     .						fsamp, ztrans)
	complex *8	ztrans(10000)
	complex *8	zbesl
	complex *8	ztboost
	complex *8	zdcblock
	complex *8	zdigfil
	complex *8	zsmooth
	complex *8	zanalog
	complex *8	zdigital
	complex *8	zxfer

	integer *4	ichan
	integer *4	micro
	integer *4	k,npts


	real *4		tau(4)
	data tau /0.00647,				!RH treble boost
     .		  0.03619,				!RL   "     "
     .		  0.00722,				!LH   "     "
     .		  0.04022/				!LL   "     "

	real *4		bes3db
	parameter (bes3db = 100.)			!Bessel corner freq
	real *4		fixed_gain
	parameter (fixed_gain = 31. * 1.3823)		!Preamp fixed gain
	real *4		bessel_gain
	parameter (bessel_gain = 1.2255 * 1.9099)	!Bessel DC gain

	real *4		dcgain
	real *4		tauboost
	real *4		dfhz
	real *4		ampl
	real *4		freqhz
c
c  generate electronics transfer function
c
	dcgain = bessel_gain * fixed_gain
	tauboost = tau(ichan)
	do k=1,npts
	   freqhz = (k-1) * dfhz
	   call bessel (freqhz, bes3db, zbesl)
	   call tboost (freqhz, tauboost, ztboost)
	   call dcblock (freqhz, zdcblock)
	   call digfltr (freqhz, micro, ichan, fsamp, zdigfil)
	   zsmooth = (1.,0.)

	   zanalog = zbesl * ztboost * zdcblock
	   zdigital = zdigfil * zsmooth
	   ztrans(k) = dcgain * zanalog * zdigital
	enddo
	return
	end


	subroutine bessel (fhz, bes3db, zbesl)
	implicit complex*8 (z)
	zfb = cmplx (0.,fhz/bes3db)
	zfb2 = zfb * zfb
	zfb3 = zfb * zfb2
	zfb4 = zfb * zfb3
	zfb5 = zfb * zfb4
	zbesl = 1./(1. + 2.4275*zfb  + 2.6189*zfb2 + 1.5894*zfb3 + 
     .			 0.5511*zfb4 + 0.0892*zfb5)
	return
	end


	subroutine tboost (fhz, tau, ztboost)
	implicit complex*8 (z)
	ztboost = 1. + (0.,6.2831853)*fhz*tau
	return
	end


	subroutine dcblock (fhz, zdcblock)
	implicit complex*8 (z)
	data twopi/6.2831853/
	zs = (0.,3.2) * twopi * fhz
	zdcblock = (zs/(1. + zs))**5 * zs/(2. + zs)
	return
	end

	
	subroutine digfltr (freqhz, micromode, ichan, samplrate, zdigfil)
	implicit complex*8 (z)
	zi = cmplx(0.,-6.2831853)
	z = cexp (zi*freqhz/samplrate)
	zdigfil = 1.
	if (micromode.eq.3 .or. micromode.eq.1 .or. 
     .					micromode.eq.0) return
	if (ichan.eq.2 .or. ichan.eq.4) then
c	   zdigfil = (1. + z*z)**2/(8. - 12.5*z*z + 5.*z**4)/8.     !Low freq
	   zdigfil = ((1. + 2.*z*(1.+ z + z*z) + z**4)/
     .			(8. - 8.*z*z + z**4))/8.
	else
	   zdigfil = (1. + z)**2/(8. - 11.*z + 5.*z*z)/2.           !High freq
	endif
	return
	end


	subroutine squeeze512 (ishift, iadds, npts, ipk_posn, shape, 
     .								tlmdat)
c
c  Do the "adds-per-grouping" down to 512 points, with peak shifted by 
c  amount ishift.  (The box should be centered on the peak when ishift = 
c  iadds/2.)  Note that 
c  we have to double the averaging window size since we calculated the 
c  spectrum up to 2 x the Nyquist frequency, thereby halving the effective 
c  time interval between raw samples.
c
	integer*4	ishift
	integer*4	iadds
	integer*4	iadd2
	integer*4	leftover
	integer*4	ipk_posn
	integer*4	npts
	integer*4	istartpoint
	integer*4	ifirst
	integer*4	j, k, n

	real*4		shape(npts)
	real*4		tlmdat(512)
	real*4		ampl
	real*4		dpk
	real*4		sum
c
c  Peak may not be an integral number of boxes from the origin; find out 
c  the remainder.
c
	iadd2 = 2 * iadds
	leftover = mod (ipk_posn, iadd2)
c
c  Leftover=1 means peak is at left edge of averaging box (ishift=1); to 
c  put peak at position M within box, start averaging at array position 
c  leftover + 1 - M.  Because of this, the first and last bins may not 
c  have all iadd2 points and must be handled separately.
c
	istartpoint = leftover + 2 - 2*ishift
	if (istartpoint .le. 0) then
	   ifirst = iadd2 + istartpoint - 1
	   istartpoint = istartpoint + iadd2
	else
	   ifirst = istartpoint - 1
	endif
c
c  Sum the first bin, move through the rest of the array, then sum the 
c  last bin.
c
	sum = 0
	do j=1, istartpoint-1
	   sum = sum + shape(j)
	enddo
	tlmdat(1) = sum

	do k=2,511
           sum = 0.
	   do n=1,iadd2
	      sum = sum + shape(iadd2*(k-2)+n+istartpoint-1)	!boxcar avg
	   enddo
	   tlmdat(k) = sum
	enddo

	sum = 0.
	do j=npts-(iadd2-istartpoint), npts
	   sum = sum + shape(j)
	enddo
	tlmdat(512) = sum
c
c  Normalize everything to peak value=1 and pass it back, along with the 
c  new peak position info in indices 511 and 512.
c
	ampl = -1
	do j=1,512
	   if (abs(tlmdat(j)) .gt. ampl) then
	      ampl = tlmdat(j)
	      ipk = j
	   endif
	enddo

	do k=1,510
	   tlmdat(k) = tlmdat(k)/ampl
	enddo
c
c  Get the "true" (i.e. noninteger) peak position by interpolation.  
c  Stuff the integer index of the peak in posn 511 and put the offset from 
c  there to the "true" peak in 512.
c
	dpk = (tlmdat(ipk-1) - tlmdat(ipk+1))/
     .			(tlmdat(ipk-1) - 2.*tlmdat(ipk) + tlmdat(ipk+1))
	tlmdat(511) = ipk 
	tlmdat(512) = dpk/2.
	return
	end
