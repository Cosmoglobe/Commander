	integer * 4 function	frd_read_variances (label, nifgs, nspec,
     .						    sum_var, sum_nk)

c-------------------------------------------------------------------------------
c
c	Function FRD_READ_VARIANCES
c
c	This function reads the real-to-real spectrum variances for pixels
c	at latitudes greater than the galactic cutoff in the FCF skymaps.  The
c	function sums the variances * N * (N-K) and sums (N-K), where N is the
c	number of ifgs in the spectra and K is the number of templates that
c	were subtracted from the antecedent coadded ifg (1 or 2).
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  3 September 1993
c
c-------------------------------------------------------------------------------
c
c	Input:
c		none
c
c	Output:
c		label		character * 40  calibration model solution label
c		nifgs		integer * 4	number of ifgs in spectra
c		nspec		integer * 4	number of spectra read
c		sum_var		real	* 8	sum of variances
c		sum_nk		real	* 8	sum of (N-K)
c
c	Subroutines called:
c		csa_close_skymap
c		csa_open_skymap
c		csa_read_pixels
c		frd_galactic_cut
c		fut_free_lun
c		fut_get_lun
c		lib$movc5
c		lib$signal
c
c	Include files:
c		csa_pixel_input_rec.txt
c		csa_pixel_output_rec.txt
c		csdr$library:ctuser.inc
c		frd_invoc_variances.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------

	implicit none

	include '(csa_pixel_input_rec)'
	include '(csa_pixel_output_rec)'
	include 'csdr$library:ctuser.inc'
	include '(fut_params)'
	include '(frd_invoc_variances)'

	character * 40  label			!  cal model solution label

	integer * 4	cstatus			!  CSA return status
	integer * 4	ct_lun			!  Cobetrieve lun
	integer * 4	io_stat			!  I/O return status
	integer * 4	j			!  a counter
	integer * 4	k			!  a counter
	integer * 4	map_ptr			!  skymap pointer
	integer * 4	nifgs			!  number of ifgs in spectra
	integer * 4	nspec			!  number of spectra read
	integer * 4	num_out			!  number of records read
	integer * 4	pixel			!  pixel number
	integer * 4	pstatus			!  return status
	integer * 4	rstatus			!  return status
	integer * 4	status			!  return status

	logical * 1	first_time /.true./	!  logical flag

	real	* 8	n_k			!  (n-k) for a spectrum
	real	* 8	sum_nk			!  sum of (n-k)
	real	* 8	sum_var(257)		!  sum of variances

	integer * 4	csa_close_skymap
	integer * 4	csa_open_skymap
	external	csa_open_skymap
	integer * 4	csa_read_pixels
	integer * 4	frd_galactic_cut
	integer * 4	fut_free_lun
	integer * 4	fut_get_lun

	dictionary 'fcf_sky'
	record /fcf_sky/ in_recs(fac_max_coad)

	record /pixel_input_list/  inlist
	record /pixel_output_list/ outlist

	external	csa_normal
	external	frd_csaclose
	external	frd_csaopen
	external	frd_csaread
	external	frd_normal
	external	fut_normal

c
c  Get the logical unit number.
c
	rstatus = fut_get_lun(ct_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif

c
c  Read in the records for each skymap.
c
	call lib$movc5 (0,,0,2056,sum_var)
	sum_nk = 0.0D0
	nspec = 0
	nifgs = 0
	map_ptr = 0
	status = %loc(frd_normal)

	do while ((status .eq. %loc(frd_normal))  .and.
     .		  (map_ptr .lt. fcc_mapnum))
	   map_ptr = map_ptr + 1

	   type 10, fcc_infile(map_ptr)
  10	   format (x, 'Processing spectra for skymap ', a)
	   open (unit=ct_lun, file=fcc_infile(map_ptr), status='old',
     .		 iostat=io_stat, form='unformatted', recordtype='fixed',
     .		 readonly, useropen=csa_open_skymap)

	   if (io_stat .eq. 0) then
c
c  Read the pixels.
c
	      cstatus = %loc(csa_normal)
	      pstatus = %loc(frd_normal)
	      inlist.level_no = fac_skymap_level
	      pixel = -1

	      do while ((cstatus .eq. %loc(csa_normal))  .and.
     .			(pstatus .eq. %loc(frd_normal))  .and.
     .			(pixel .lt. 6143))
	         pixel = pixel + 1
	         if (fcc_galexc .eq. fac_present) then
	            pstatus = frd_galactic_cut (pixel)
	         endif
	         if (pstatus .eq. %loc(frd_normal)) then
	            inlist.pixel_no = pixel
	            cstatus = csa_read_pixels (ct_lun, inlist, 1, in_recs,
     .					      fac_max_coad, outlist, num_out, 0)

	            if (outlist.no_records .gt. 0) then
	               if (cstatus .eq. %loc(csa_normal)) then
c
c  Get the calibrtion model solution label the first time through.
c
	                  if (first_time .eq. .true.) then
	                     first_time = .false.
	                     label = in_recs(1).spec_data.model_label
	                  endif

c
c  Sum the real-to-real variances.
c
	                  do j = 1,outlist.no_records
	                     nifgs = nifgs + in_recs(j).coad_spec_head.num_ifgs
	                     n_k = dble(in_recs(j).coad_spec_head.num_ifgs - 1 -
     .			      in_recs(j).coad_spec_data.sec_template.subtracted)
	                     do k = 1,257
	                        sum_var(k) = sum_var(k) +
     .				  (dble(in_recs(j).spec_data.real_var(k)) *
     .				 dble(in_recs(j).coad_spec_head.num_ifgs) * n_k)
	                     enddo
	                     sum_nk = sum_nk + n_k
	                  enddo
	                  nspec = nspec + outlist.no_records
	               else
	                  status = %loc(frd_csaread)
	                  call lib$signal (frd_csaread, %val(2),
     .				      fcc_infile(map_ptr)(1:fcc_inlen(map_ptr)),
     .							%val(cstatus))
	               endif	!  return status from read
	            endif	!  no_records ne 0
	         endif		!  return status from galactic cut
	      enddo		!  loop over pixels

c
c  Close the input skymap file.
c
	      cstatus = csa_close_skymap (ct_lun, fac_skymap_no_levels)
	      if (cstatus .ne. %loc(csa_normal)) then
	         status = %loc(frd_csaclose)
	         call lib$signal (frd_csaclose, %val(2),
     .				  fcc_infile(map_ptr)(1:fcc_inlen(map_ptr)),
     .						%val(cstatus))
	      endif

	   else
	      status = %loc(frd_csaopen)
	      call lib$signal (frd_csaopen, %val(2),
     .			       fcc_infile(map_ptr)(1:fcc_inlen(map_ptr)),
     .					    %val(io_stat))
	   endif		!  (io_stat from input skymap open

	enddo			!  loop over skymaps.

c
c  Free the logical unit number for the input skymaps.
c
	rstatus = fut_free_lun(ct_lun)
	if (rstatus .ne. %loc(fut_normal)) then
	   call lib$signal (%val(rstatus))
	endif


	frd_read_variances = status

	return
	end
