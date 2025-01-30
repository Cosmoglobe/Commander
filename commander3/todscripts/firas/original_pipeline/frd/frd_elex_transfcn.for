	Program frd_elex_transfcn
c
c  This routine creates the data file ELEX_TRANSFCN.DAT containing the FIRAS
c  post-detector electronics transfer function. The function is derived from the
c  analytical expressions of the various analog and digital filters; these are
c  (1) a treble boost, (2) a 5-pole Bessel filter, (3) a 6-pole DC-blocking 
c  filter, (4) a digital lowpass filter, and (5) smoothing due to data 
c  compression.
c
c  Each record of the output file is an array of 257 complex*8 numbers
c  representing the complex transfer function for a particular scan mode at
c  frequencies ranging from DC up to the Nyquist frequency associated with
c  the data compression.
c
c  There are 288 records in an unformatted direct-access output file. These
c  are:
c		3 sample rates (fakeit, fast scan, and slow scan)
c	    x	4 channels
c	    x   2 microprocessor science modes (digital filters on or off)
c	    x  12 possible amounts of data compression ("adds per group" =1to12)
c
c  Written by Rich Isaacman (Applied Research Corp.)
c  Modified 7 July 1986 to make DC gain 240 instead of 100 per measurement
c  Modified 28 October 1987 to treat all four channels separately and
c	update digital filters and fixed gain.
c  Modified by J.T.Bonnell (SSTX@GSFC) nov.1988
c		brought up to standards under spr 1814
c  Modified by R.Kummerer (STX@GSFC) may 1989 Remove useless TYPE statement.
c  SPR 5291, Generate electronics transfer function for new sampling rate.
c		R. Kummerer, Dec 15, 1989, STX@GSFC.
c  SPR 8919, Add qualifier /SAMPLE=INT or, MISSION ; If qualifier = INT, then
c            sampling rate is 684.00, else sampling rate is 681.43.
c               Nilo G. Gonzales/STX, August 26, 1991.
c  Converted to double precision.  Gene Eplee, GSC, 12 September 1991.
c  Removed complex conjugate of ETF.  Gene Eplee, GSC, 23 January 1992.
c  Updated the code to read new FEX_SAMPRATE text file. SPR 9846.
c               Nilo G. Gonzales/Hughes STX, July 31, 1992.
c  Fixed reading of FEX_SAMPRATE text file and parsing of SAMPLE qualifier.
c     SPR 10575, Steve Alexander, Hughes STX, 2/12/93.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	implicit none

	include '($ssdef)'
	include '(upm_stat_msg)'

	complex *16	ztrans(257)
	complex *16	zbesl
	complex *16	ztboost
	complex *16	zdcblock
	complex *16	zdigfil
	complex *16	zsmooth
	complex *16	zanalog
	complex *16	zdigital
	complex *16	zxfer

	integer *4	status
	integer *4      upm_stat
	integer *4	ios
	integer *4	lun,lun_samp
	integer *4	nrec
	integer *4	isampl
	integer *4	ichan
	integer *4	micro
	integer *4	ncompress
	integer *4	k
	integer *2      len
	logical *1      mission /.false./, int /.false./
	character *6    sample
	character *7    qual
	character *16   in_file /'fex_samprate.txt'/
	character *11	out_file /'fex_etf.dat'/
	real *8		int_samp_hz(3)
	real *8		mis_samp_hz(3)
	real *8		int_samp_read, mis_samp_read

	real *8		tau(4)
	data tau /0.00647D0,				!RH treble boost
     .		  0.03619D0,				!RL   "     "
     .		  0.00722D0,				!LH   "     "
     .		  0.04022D0/				!LL   "     "

	real *8		bes3db
	parameter (bes3db = 100.0D0)			!Bessel corner freq
	real *8		fixed_gain
	parameter (fixed_gain = 31.0D0 * 1.3823D0)	!Preamp fixed gain
	real *8		bessel_gain
	parameter (bessel_gain = 1.2255D0 * 1.9099D0)	!Bessel DC gain

	real *8		dcgain
	real *8		fsamp
	real *8		tauboost
	real *8		dfhz
	real *8		ampl
	real *8		ampmax
	real *8		freqhz

	integer *4	lib$get_lun
	integer *4	lib$free_lun
	integer *4      upm_get_value
	integer *4      upm_present

	external	frd_normal
	external	frd_rmsopen
	external	frd_rmswrite
	external	frd_rmsclose

c Parse the command line.

        mission = .true.
	if (upm_present('sample') .eq. upm_pres) then
           upm_stat = upm_get_value('sample',qual,len)
	   if (upm_stat .ne. ss$_normal) then
	      call lib$signal (%val(ss$_abort))
	   end if
	   if (qual .eq. 'INT') then
	      int = .true.
	      mission = .false.
	   end if
	end if	
c
c open fex_samprate text file for mtm sampling rate (INT, or MISSION) 
c
	status = lib$get_lun(lun_samp)
	if (status .eq. ss$_normal) then
	   open (unit=lun_samp, name=in_file, status='old',
     .	         iostat=ios, readonly, shared)
	      if (ios .ne. 0) then
	         status = %loc(FRD_RMSOpen)
	         call lib$signal(FRD_RMSOpen,%val(2),in_file,%val(ios))
	      end if
	else
	    call lib$signal(%val(status))
	end if
c
c read rms file
c
        read (lun_samp,*) int_samp_read
        read (lun_samp,*) mis_samp_read

	if (int) then
           int_samp_hz(1) = int_samp_read
           int_samp_hz(2) = int_samp_read
           int_samp_hz(3) = 512.00D0
	else
           mis_samp_hz(1) = mis_samp_read
	   mis_samp_hz(2) = mis_samp_read
           mis_samp_hz(3) = 512.00D0 
	end if
c
c  open output file
c
	status = lib$get_lun(lun)

	if (status .eq. ss$_normal) then

	open (unit=lun,name=out_file,status='new',form='unformatted',
     .			iostat=ios,access='direct',recl=1028)

	if (ios .eq. 0) then

c
c  generate electronics transfer function
c

	dcgain = bessel_gain * fixed_gain
	nrec = 0
	do isampl=1,3				!  1,2=MTM modes, 3=Fakeit 
	   if (int) then
	      fsamp = int_samp_hz(isampl)       !  INT sampling rate
	   else
	      fsamp = mis_samp_hz(isampl)       !  Mission sampling rate
	   end if

	   do ichan=1,4				!  RH,RL,LH,LL channels
	   tauboost = tau(ichan)

	      do micro=1,2			!  1=dig fltr on, 2=off

		 do ncompress=1,12		!  adds per group
c
c  stop processing on error
c
		 if(status .eq. ss$_normal)then

		    nrec = nrec + 1
		    dfhz = fsamp/ncompress/512.0D0
		    ampmax = -1.0D0

		    do k=1,257
		       freqhz = dble(k-1) * dfhz
		       call bessel (freqhz, bes3db, zbesl)
		       call tboost (freqhz, tauboost, ztboost)
		       call dcblock (freqhz, zdcblock)
		       call digfltr (freqhz, micro, ichan, fsamp, zdigfil)
		       call compress (freqhz, ncompress, fsamp, zsmooth)

		       zanalog = zbesl * ztboost * zdcblock
		       zdigital = zdigfil * zsmooth
		       zxfer = dcgain * zanalog * zdigital

	   	       ampl = dsqrt (dreal (zxfer * dconjg (zxfer)))
	   	       ampmax = dmax1 (ampl, ampmax)
	   	       if (ampl .lt. ampmax/1000.0D0) then
	       		  ztrans(k) = 100000.0D0
		       else
			  ztrans(k) = -zxfer
		       endif
		    enddo

!	   	    type *, nrec
	   	    write (lun, rec=nrec, iostat=ios) ztrans
			if(ios .ne. 0)then
			  call lib$signal(frd_rmswrite,%val(2),
	1			out_file,%val(ios))
			  status = %loc(frd_rmswrite)
			endif

	      endif	!status .ne. ss$_normal

		 enddo	!k=1,257
	      enddo	!ncompress
	   enddo	!micro
	enddo		!ichan

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

	stop
	end

	subroutine bessel (fhz, bes3db, zbesl)
	implicit complex*16 (z)
	implicit real*8 (a-h,o-y)
	zfb = dcmplx (0.0D0,fhz/bes3db)
	zfb2 = zfb * zfb
	zfb3 = zfb * zfb2
	zfb4 = zfb * zfb3
	zfb5 = zfb * zfb4
	zbesl = 1.0D0/(1.0D0 + 2.4275D0*zfb  + 2.6189D0*zfb2 + 1.5894D0*zfb3 + 
     .			 0.5511D0*zfb4 + 0.0892D0*zfb5)
	return
	end


	subroutine tboost (fhz, tau, ztboost)
	implicit complex*16 (z)
	implicit real*8 (a-h, o-y)
	ztboost = 1.0D0 + dcmplx(0.0D0,6.2831853D0)*fhz*tau
	return
	end


	subroutine dcblock (fhz, zdcblock)
	implicit complex*16 (z)
	implicit real*8 (a-h, o-y)
	data twopi/6.2831853D0/
	zs = dcmplx(0.0D0,3.2D0) * twopi * fhz
	zdcblock = (zs/(1.0D0 + zs))**5 * zs/(2.0D0 + zs)
	return
	end

	
	subroutine digfltr (freqhz, micromode, ichan, samplrate, zdigfil)
	implicit complex*16 (z)
	implicit real*8 (a-h, o-y)
	zi = dcmplx(0.0D0,-6.2831853D0)
	z = cdexp (zi*freqhz/samplrate)
	zdigfil = dcmplx(1.0D0, 0.0D0)
	if (micromode .eq. 2) return
	if (ichan.eq.2 .or. ichan.eq.4) then
c	   zdigfil = (1.0D0 + z*z)**2/(8.0D0 - 12.5D0*z*z + 5.0D0*z**4)/8.0D0
	   zdigfil = ((1.0D0 + 2.0D0*z*(1.0D0+ z + z*z) + z**4)/   !Low freq
     .			(8.0D0 - 8.0D0*z*z + z**4))/8.0D0
	else
	   zdigfil = (1.0D0 + z)**2/(8.0D0 - 11.0D0*z + 5.0D0*z*z)/2.0D0
								    !High freq
	endif
	return
	end


	subroutine compress  (fhz, ncompress, samplrate, zsmooth)
	implicit complex*16 (z)
	implicit real*8 (a-h, o-y)
	pif = 3.141592654D0  * fhz/samplrate
	zsmooth = dcmplx(1.0D0,0.0D0)
	if (pif .lt. 0.001D0) return
	ampl = dsin(pif*dble(ncompress))/dsin(pif)
	zphase = dcmplx (0.0D0, pif*dble(1-ncompress))
	zsmooth = ampl/dble(ncompress) * cdexp(zphase)
	return
	end
