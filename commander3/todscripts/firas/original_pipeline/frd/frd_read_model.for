	integer * 4 function  frd_read_model ()

c-------------------------------------------------------------------------------
c
c	Function FRD_READ_MODEL
c
c	This function reads Dale Fixsen's Firas calibration model solution from
c	the RMS file FEX_MOD_CCSS.TXT_VVV_XXXXXXXXXX, where CC is the channel,
c	SS is the scan mode, VVV is the version of the calibration program that
c	generated the solution, and XXXXXXXXXX is the model solution label. 
c	The file extension is specified by the user via the command line
c	qualifier /FILE_EXT=VVV_XXXXXXXXXX.
c
c	The model solution bolometer parameters are renormalized from the those
c	produced by the minimization routines.  The emissivities are shifted by
c	one sample so that the first sample is for a frequency of zero icm.
c	The emissivities are normalized by the Nyquist frequency and by the
c	throughput of the instrument.  Kirchoff's law is applied to calculate
c	the emissivity of the bolometer.  Finally, the model solution is
c	inserted into the RDL FEX_MOD and written into the binary file
c	FEX_MOD_CCSS.VVV_XXXXXXXXXX.
c
c	The input text file is read from the directory CSDR$FIRAS_IN, while the
c	output binary file is written to the directory CSDR$FIRAS_OUT.
c
c	Author:
c		Gene Eplee
c		General Sciences Corp.
c		513-7768
c		2 October 1992
c		SER 8178
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		none
c
c	Subroutines called:
c		fut_free_lun
c		fut_get_lun
c		str$trim
c		str$upcase
c
c	Include files:
c		frd_model_config.txt
c		frd_model_invoc.txt
c		frd_model_soln.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 22 October 1993
c	SER 11394
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(frd_model_config)'
	include '(frd_model_invoc)'
	include '(frd_model_soln)'

	character * 23	datestr			!  model solution time tag
	character * 60	in_file			!  input parameter file name
	character * 80	label			!  model file label

	complex * 16	sum_emiss		!  sum of emissivities

	integer *  2	inlen			!  input file name length
	integer *  2	outlen			!  output file extension length

	integer *  4	colon			!  colon label index
	integer *  4	fstart			!  final file extension
						!    start label index
	integer *  4	fstop			!  final file extension
						!    stop label index
	integer *  4	istart			!  initial file extension
						!    start label index
	integer *  4	istop			!  initial file extension
						!    stop label index
	integer *  4	exclam			!  exclamation point label index
	integer *  4	frec			!  Nyquist frequency record
						!    number
	integer *  4	i			!  a counter
	integer *  4	in_lun			!  input file lun
	integer *  4	io_stat			!  I/O status
	integer *  4	is			!  upper limit index for soln
	integer *  4	j			!  a counter
	integer *  4	k			!  a counter
	integer *  4	lstart			!  starting label index
	integer *  4	lstop			!  stopping label index
	integer *  4	nlines			!  number of lines in model file
	integer *  4	npar			!  number of bolometer params
	integer *  4	rstatus			!  return status
	integer *  4	smode			!  mtm scan mode
	integer *  4	speed			!  mtm scan speed
	integer *  4	status			!  return status

	real	*  8	buff(4)			!  model solution buffer
	real	*  8	etendu			!  instrument throughput
	parameter	(etendu = 1.5D0)
	real	*  8	fnyq			!  Nyquist frequency
	real	*  8	normalize		!  emiss normalization factor
	real	*  8	param(mpar)		!  input bolometer parameters
	real	*  8	xferi(6,257)		!  imag part of emissivities
	real	*  8	xferr(6,257)		!  real part of emissivities

	integer *  4	fut_free_lun
	integer *  4	fut_get_lun

	external	frd_normal
	external	frd_rmsclose
	external	frd_rmsopen
	external	frd_rmsread
	external	frd_rmsrew
	external	frd_rmssize
	external	fut_normal

	status = %loc(frd_normal)

C
C  Set the instrument state parameters.
C
	smode = fcc_smode
	if (smode .eq. 5) smode = 2
	is = 3*(fcc_chan-1) + (fcc_speed+1) + fcc_length


C
C  Read in the calibration model.
C

c
c  Determine the model parameter file name.
c
	in_file = 'csdr$firas_in:fex_mod_' // fac_channel_ids(fcc_chan) //
     .				 fac_scan_mode_ids(fcc_smode) // '.TXT_' //
     .				 fcc_infile_ext
	call str$trim (in_file, in_file, inlen)

c
c  Open the model parameter file.
c
	rstatus = fut_get_lun(in_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal(rstatus)
	endif

	open (unit=in_lun, file=in_file, status='old', form='formatted',
     .        access='sequential', readonly, iostat=io_stat)

	if (io_stat .ne. 0) then
	   status = %loc(frd_rmsopen)
	   call lib$signal (frd_rmsopen, %val(2), in_file(1:inlen),
     .					 %val(io_stat))
	endif

	if (status .eq. %loc(frd_normal)) then
c
c  Read the model parameter file to determine the number of bolometer
c  	parameters.
c
	   nlines = 0
	   read (in_lun, *, iostat=io_stat)
	   read (in_lun, *, iostat=io_stat)
	   do while (io_stat .eq. 0)
	      read (in_lun, 10, iostat=io_stat) buff
	      if (io_stat .eq. 0) nlines = nlines + 1
	   enddo
  10	   format (3(e19.13,x),e19.13)
	   if (io_stat .lt. 0) then
	      npar = 4 * (nlines - (3 * uplim(is)))
	      rewind (in_lun, iostat=io_stat)
	      if (io_stat .ne. 0) then
	         status = %loc(frd_rmsrew)
	         call lib$signal (frd_rmsrew, %val(2),
     .			          in_file(1:inlen), %val(io_stat))
	      endif
	   else
	      status = %loc(frd_rmssize)
	      call lib$signal (frd_rmssize, %val(2),
     .			       in_file(1:inlen), %val(io_stat))
	   endif

	   if (status .eq. %loc(frd_normal)) then
c
c  Read in the model parameters from the variable-length parameter files.
c
	      read (in_lun, 20, iostat=io_stat) label
	      read (in_lun, *, iostat=io_stat)
	      read (in_lun, 30, iostat=io_stat) (param(i),i=1,npar),
     .		 ((xferr(j,k),xferi(j,k),j=1,6),k=1,uplim(is))
  20	      format (a)
  30	      format (3(e19.13,x),e19.13)

	      if (io_stat .ne. 0) then
	         status = %loc(frd_rmsread)
	         call lib$signal (frd_rmsread, %val(2), in_file(1:inlen),
     .					       %val(io_stat))
	      endif


	      if (status .eq. %loc(frd_normal)) then
C
C  Fill in the model solution header information
C

c
c  Get the initial part of the output file name extension
c
                 colon = index(label, ':')
	         istop = colon - 3
	         do while ((label(istop:istop) .eq. ' ')  .and.
     .			   (istop .gt. 2))
	            istop = istop - 1
	         enddo
	         istart = istop - 1
	         do while ((label(istart:istart) .ne. ' ')  .and.
     .			   (istart .gt. 2))
	            istart = istart - 1
	         enddo
	         istart = istart + 1

c
c  Get the final part of the output file name extension.
c
	         fstop = istart - 1
	         do while ((label(fstop:fstop) .eq. ' ')  .and.
     .			   (fstop .gt. 2))
	            fstop = fstop - 1
	         enddo
	         fstart = fstop - 1
	         do while ((label(fstart:fstart) .ne. ' ')  .and.
     .			   (fstart .gt. 2))
	            fstart = fstart - 1
	         enddo
	         fstart = fstart + 1
	         fcc_outfile_ext = label(istart:istop) // '_' //
     .				   label(fstart:fstop)
	         call str$upcase (fcc_outfile_ext, fcc_outfile_ext)
	         call str$trim (fcc_outfile_ext, fcc_outfile_ext, outlen)

c
c  Get the model header label.
c
	         exclam = index(label, '!')
	         lstart = exclam + 1
	         do while ((label(lstart:lstart) .eq. ' ')  .and.
     .			   (lstart .lt. 40))
	            lstart = lstart + 1
	         enddo
	         fex_model.mod_head.label = label(lstart:fstop) // '    ' //
     .					    label(istart:istop)
	         call str$upcase (fex_model.mod_head.label,
     .				  fex_model.mod_head.label)

c
c  Fill in the remaining model header information.
c
	         fex_model.mod_head.chan_id        = fcc_chan
	         fex_model.mod_head.mtm_length     = fcc_length
	         fex_model.mod_head.mtm_speed      = fcc_speed
	         fex_model.mod_head.adds_per_group = fcc_ngroup

c
c  Get the model solution time tag.
c
	         lstart = colon + 2
	         lstop = colon + 24
	         datestr = label(lstart:lstop)
	         call sys$bintim (datestr, fex_model.ct_head.time)
	         call ct_binary_to_gmt (fex_model.ct_head.time,
     .					fex_model.ct_head.gmt)
	         fex_model.mod_head.gmt = fex_model.ct_head.gmt


C
C  Renormalize the calibration model.
C

c
c  Renormalize the bolometer parameters.
c
	         do j = 1, npar
	            fex_model.bolparm(j) = param(j) * renormalize(j)
	         enddo

c
c  Normalize the emissivities by the Nyquist frequency and the etendu.
c	Shift the emissivities so that the first sample is for zero frequency.
c	Switch the sign of the emissivities for the right side model solutions.
c
	         frec = 4*jmod((fcc_chan-1),2) + smode
	         fnyq = dble(nyquist.icm(frec))
	         if (fcc_chan .le. 2) then
                    normalize = - fnyq * etendu
	         elseif (fcc_chan .ge. 3) then
                    normalize = fnyq * etendu
	         endif
	         do j = 1, 6
	            do k = 1, 256
	               fex_model.emissivity(k+1,j) =
     .					 dcmplx(xferr(j,k),xferi(j,k))/normalize
	            enddo
	         enddo

c
c  Calculate the bolometer emissivity.
c
	         do j = 1, 257
	            sum_emiss = dcmplx(0.0D0,0.0D0)
	            do k = 1, 6
	               sum_emiss = sum_emiss + fex_model.emissivity(j,k)
	            enddo
	            fex_model.emissivity(j,7) = - sum_emiss
	         enddo

	      endif	!  (status from reading of model

	   endif	!  (status from sizing of model

c
c  Close the model parameter file.
c
	   close (in_lun, iostat=io_stat)
	   if (io_stat .ne. 0) then
	      status = %loc(frd_rmsclose)
	      call lib$signal (frd_rmsclose, %val(2), in_file, %val(io_stat))
	   endif
	   rstatus = fut_free_lun(in_lun)
	   if (rstatus .ne. %loc(fut_normal)) then
	      call lib$signal(rstatus)
	   endif

	endif		!  (status from opening of model


	frd_read_model = status

	return
	end
