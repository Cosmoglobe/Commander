	integer * 4 function fnt_convert_to_volts(chan_id,spec_rec)

c--------------------------------------------------------------------------
c	Author: W. K. Young
c		STI Inc.
c		July 1986
c
c	Input:
c		Channel ID
c		Coadded spectrum
c	Output:
c		Converted spectrum
c
c	Include files:
c		fnt_invoc.txt
c		fut_params.txt
c	Calling Sequence
c		status = fnt_convert_to_volts ( chan_id,
c						spec_rec )
c	Gist of routine
c		IF gain, adds per group, and sweeps non zero THEN
c		  DO for all frequencies
c		      scale spectrum to come out with V/Hz**0.5
c		  END DO
c	        ELSE
c		  SET return status to ERROR
c		END IF
c--------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fnt_invoc)'

	integer * 4 i			! a counter
	integer * 4 sweeps		! sweeps per IFG
	integer * 4 chan_id		! channel ID
	real * 4 actual_gain(8)		!actual gain values
	real * 4 fudge_fac		!factor that produces right results
	integer * 4 gain
	integer * 4 ngroup

	dictionary 'fnt_noise'
	record /fnt_noise/ spec_rec

	data actual_gain/1.,3.,10.,30.,100.,300.,1000.,3000./

	fnt_convert_to_volts = fac_normal
        sweeps = spec_rec.chan.sweeps
	if(sweeps .eq. 0) then
	   sweeps = 1
	end if
	gain = spec_rec.chan.gain
	ngroup = spec_rec.chan.ngroup
	if( gain .ge. 0 .and. gain .lt. 8 .and. sweeps .gt. 0 .and.
     .		   ngroup .gt. 0)then
           if (ftc_preamp .eq. fac_present) then
	     fudge_fac = sqrt(float(sweeps))/(204.8*actual_gain(gain+1))
	   else
	     fudge_fac = sqrt(float(sweeps))/204.8
	   end if
	   do i = 1,257
	      spec_rec.chan.spec(i) =
     .			spec_rec.chan.spec(i)*fudge_fac
	   end do
	else
	   type *,'Didn''t convert to volts'
	   type *,'gain was ',gain
	   type *,'sweeps were ',sweeps
	   type *,'adds per group were ',ngroup
	   fnt_convert_to_volts = 13
	end if
	return
	end

