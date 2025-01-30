	integer * 4 function fnt_write_noise(chan,archive_type,
	1				     first_IFG, last_IFG,
	2				     spec_rec)

c--------------------------------------------------------------------------
c
c	Function fnt_write_noise
c	This function writes spectral data to the archive in the standard
c		CT spectrum format.
c
c	Author: W. K. Young
c		STI Inc.
c		July 1986
c
c	Calling sequence:
c		status = fnt_write_noise(channel id,
c					 spectrum)
c	Input:	channel id
c		spectrum
c
c	Output:	return status
c
c	Include Files:
c		fut_params.txt
c		fnt_invoc.txt
c		ct$library:ctuser.inc
c
c	Gist of routine:
c		SET return status to normal
c		ASSIGN data from temp variables to spectral record format
c			variables
c		OPEN up output file
c		IF no error occurs THEN
c		   WRITE record to output file
c		   IF an error occurs THEN
c		      SET return status to error
c		      CALL error reporting
c		   END IF
c		ELSE
c                  SET return status to error
c	           CALL error reporting
c               END IF
c               RETURN
c--------------------------------------------------------------------------
c
c Changes:
c
c	Write to archive capability. R. Kummerer, May 5, 1987.
c
c	Fill dataset ID. R. Kummerer, Sep 23, 1987.
c
c	Archive split into separate files, one for each channel.
c		F. Shuman, 1987 Dec 11.
c
c	Remove Write-to-RMS.  F. Shuman, 1988 Mar 29.
c
c	Version 4.4.02, SPR 5063, Nov 16, 1989, R. Kummerer, STX
c		Direct output with CSDR$FIRAS_OUT.
c
c--------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fnt_invoc)'
	include 'ct$library:ctuser.inc'

	integer * 4	chan		!current channel # = 1 to 4
	integer * 2	fir_chan	!chan # from ctuser params file
	integer * 4	ct_unit		!logical write unit number
	integer * 2	ct_stat(20)	!ct return status
	character * 2	archive_type	!RL, TD, ED etc. archive
	integer * 4	first_IFG(2)	!timetag of first IFG in spectrum
	integer * 4	last_IFG(2)	!timetag of last IFG in spectrum
	character * 14	start_IFG
	character * 14	stop_IFG

	integer * 4	ios
	character * 64	noise_file
	integer * 4	CT_Connect_Write
	external	CT_Connect_Write

	dictionary 'FNT_NOISE'
	record /FNT_NOISE/ spec_rec	!power density spectra

C Create an output file name.

	fnt_write_noise = fac_normal

	if (ftc_write .eq. fac_present) then

	   Call CT_Binary_To_GMT ( first_IFG, start_IFG )
	   Call CT_Binary_To_GMT ( last_IFG, stop_IFG )

	   noise_file = 'CSDR$FIRAS_OUT:FNT_NOISE_' // fac_channel_ids(chan) //
	1		'.' // archive_type // '_' // start_IFG(1:7) // '_' //
	2		stop_IFG(1:7)

	   if (chan .eq. 1) then
	      fir_chan = ctu_$fir_nrh
	   else if (chan .eq. 2) then
	      fir_chan = ctu_$fir_nrl
	   else if (chan .eq. 3) then
	      fir_chan = ctu_$fir_nlh
	   else if (chan .eq. 4) then
	      fir_chan = ctu_$fir_nll
	   end if

	   spec_rec.ct_head.dataset_id = fir_chan

C Write to archives.

	   Call LIB$Get_LUN(ct_unit)

	   Open ( Unit=ct_unit, File=noise_file, Status='New', IOStat=ios,
	1	  UserOpen=CT_Connect_Write )

	   if (ios .eq. 0) then
	      call ct_write_arcv(,ct_unit,spec_rec,ct_stat)
	      if(ct_stat(1) .ne. ctp_normal)then
	         type *,'Noise spectrum archive write error, Status= ',ct_stat(1)
	      end if
	      call ct_close_arcv(,ct_unit,ct_stat)
	      if(ct_stat(1) .ne. ctp_normal)then
	         type *,'Noise spectrum archive close error, Status= ',ct_stat(1)
	      end if
	   else
	      type *,'Noise spectrum archive open error, Status= ', ios
	   end if

	end if
	
	return
	end
