module comm_zodi_samp_mod
   use comm_zodi_mod
   use comm_chisq_mod
   use powell_mod
   implicit none

   private
   public initialize_zodi_samp_mod, downsamp_invariant_structs, project_and_downsamp_sky, compute_downsamp_zodi
   public minimize_zodi_with_powell, get_chisq_priors, precompute_lowres_zodi_lookups, create_zodi_glitch_mask
   public apply_zodi_glitch_mask, sample_static_zodi_map!, sample_static_zodi_amps
   
   real(dp), allocatable :: chisq_previous, step_size, prior_vec(:, :), prior_vec_powell(:, :), step_sizes_emissivity(:, :), step_sizes_albedo(:, :), step_sizes_ipd(:), step_sizes_n0(:), theta_0(:)
   real(dp), allocatable :: powell_emissivity(:, :), powell_albedo(:, :)
   integer(i4b) :: n_samp_bands, reference_emissivity_band
   real(dp) :: EPS = TINY(1.0_dp)
   real(dp), dimension(2) :: emissivity_prior, albedo_prior
   logical(lgt), allocatable :: powell_included_params(:), ref_band(:)
   character(len=128), allocatable, dimension(:) :: implemented_sampling_algorithms
   
contains
  subroutine initialize_zodi_samp_mod(cpar)
    implicit none
      ! Initialize the zodi sampling module.
      !
      ! Parameters
      ! ----------
      ! cpar: comm_params
      !    Parameter file variables.

      type(comm_params), intent(inout) :: cpar
      integer(i4b) :: i, j, idx_start, idx_stop, ref_band_count, ierr, group_idx
      real(dp), allocatable :: param_vec(:)
      character(len=128), allocatable :: labels(:)
      integer(i4b), allocatable :: indices(:)
      ! Figure out how many sampling bands there are and initialize the tod step sizes

      implemented_sampling_algorithms = ["powell", "mh"]
      if (.not. any(implemented_sampling_algorithms == cpar%zs_sample_method)) then
         if (cpar%myid == 0) then 
            print *, "Error: invalid sampling method for zodi, must be one of: ", [(trim(adjustl(implemented_sampling_algorithms(i)))//", ", i=1, size(implemented_sampling_algorithms))]
            stop
         end if
      end if
      n_samp_bands = 0
      do i = 1, numband
         if (data(i)%tod_type == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle
         n_samp_bands = n_samp_bands + 1
      end do
      ref_band = cpar%ds_zodi_reference_band
      ref_band_count = count(cpar%ds_zodi_reference_band == .true.)
      if (trim(adjustl(cpar%zs_sample_method)) == "mh") then
         if (ref_band_count > 1) then
            stop "Error: cannot have more than one reference band for zodi emissivity."
         else if (ref_band_count == 0) then
            stop "Error: cannot sample zodi without the reference band being active."
         end if
      end if
      
    end subroutine initialize_zodi_samp_mod

   function get_boxwidth(samprate_lowres, samprate) result(box_width)
      ! Returns the boxcar width for downsampling the zodi tods
      real(dp), intent(in) :: samprate_lowres, samprate
      real(dp) :: box_width
      box_width = max(real(nint(samprate/samprate_lowres),dp),1.d0)
   end function



   function get_chisq_priors(params, samp_group) result(chisq_prior)
     implicit none 
     real(dp), intent(in) :: params(:)
     integer(i4b), intent(in) :: samp_group
      
      real(dp) :: chisq_prior
      integer(i4b) :: i, j
      logical(lgt) :: prior_is_violated

      chisq_prior = 0.d0
      j = 1
      do i = 1, zodi_model%npar_tot
         if (zodi_model%theta_stat(i,samp_group) /= 0) cycle

         if (params(j) < zodi_model%theta_prior(1,i) .or. &
              & params(j) > zodi_model%theta_prior(2,i)) then
            write(*,fmt='(a,a, 3e16.8)') 'Out of bounds -- ', trim(zodi_model%par_labels(i))//':', zodi_model%theta_prior(1,i), params(j), zodi_model%theta_prior(2,i)
            !write(*,*) 'Parameter out of bounds', i, zodi_model%theta_prior(1,i), params(j), zodi_model%theta_prior(2,i)
            chisq_prior = 1.d30
            return
         end if
         
         if (zodi_model%theta_prior(4,i) > 0.d0) then
            chisq_prior = chisq_prior + (params(j) - zodi_model%theta_prior(3,i))**2 / zodi_model%theta_prior(4,i)**2
         end if

         j = j+1
      end do

      if (trim(zodi_model%phasefunc_type) == 'Hong') then
         ! Check that sum of weights is smaller than 1
         if (params(4)+params(5) > 1.d0) then
            chisq_prior = 1.d30
            return
         end if
      end if
    end function get_chisq_priors

   subroutine downsamp_invariant_structs(cpar)
      ! Downsamples pointing, tod and caches the zodi mask (sky mask + flags)
      type(comm_params), intent(in) :: cpar

      integer(i4b) :: i, j, k, l, g, scan, npix, nmaps, ndelta, ext(2), padding, ierr, ntod, ndet, nhorn, ndownsamp, box_halfwidth
      real(dp) :: box_width, dt_tod, theta, phi, vec0(3), vec1(3), M_ecl2gal(3,3), day, lambda_solar, lat, lon
      real(sp), allocatable :: tod(:), mask(:), downsamp_vec(:, :), s_static(:)
      integer(i4b), allocatable :: pix(:, :), psi(:, :), flag(:)
      real(dp), allocatable, dimension(:, :) :: m_buf, vec
      integer(i4b), allocatable, dimension(:) :: downsamp_pix
      logical(lgt), allocatable, dimension(:) :: downsamp_mask_idx
      real(sp), allocatable, dimension(:) :: downsamp_mask, downsamp_tod, downsamp_obs_time, obs_time
      real(sp), allocatable, dimension(:) :: procmask_zodi
      type(hdf_file) :: tod_file
      character(len=4) :: scan_str

      padding = 5
      if (cpar%myid == cpar%root) print *, "downsampling tod and pointing"
      ! For each zodi band, create downsampled residual time stream

      call ecl_to_gal_rot_mat(M_ecl2gal)
      
      do i = 1, numband
         ! Only generate downsampled arrays for tod bands and bands where we have zodi
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         ndelta = 1
         npix = 12*data(i)%map%info%nside**2
         nmaps = data(i)%map%info%nmaps
         ndet = data(i)%tod%ndet
         nhorn = data(i)%tod%nhorn

         ! Build zodi mask
         allocate (m_buf(0:npix - 1, nmaps))
         allocate (procmask_zodi(0:npix - 1))
         call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:, 1)
         box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
         box_halfwidth = anint(box_width / 2.)
         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            ntod = data(i)%tod%scans(scan)%ntod
            ! Allocate temporary variables to be downsampled
            allocate (pix(ntod, nhorn))
            allocate (psi(ntod, nhorn))
            allocate (flag(ntod))
            allocate (tod(ntod))
            allocate (s_static(ntod))
            allocate (mask(ntod))
            allocate (vec(3, ntod))

            ! Get shape of downsampled structs
            call data(i)%tod%downsample_tod(tod, ext, width=box_width)

            ! Allocate downsampled strucs
            allocate (downsamp_tod(ext(1):ext(2)))
            allocate (downsamp_mask(ext(1):ext(2)))
            allocate (downsamp_pix(ext(1):ext(2)))
            allocate (downsamp_mask_idx(ext(1):ext(2)))
            allocate (downsamp_vec(3, ext(1):ext(2)))

            ! downsample obs_time
            allocate (obs_time(ntod))
            dt_tod = (1./data(i)%tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
            do k = 1, ntod
               obs_time(k) = data(i)%tod%scans(scan)%t0(1) + (k - 1)*dt_tod
            end do
            allocate (downsamp_obs_time(ext(1):ext(2)))
            call data(i)%tod%downsample_tod(obs_time, ext, downsamp_obs_time, width=box_width)

            do j = 1, data(i)%tod%ndet
               pix = 0.
               psi = 0.
               flag = 0
               tod = 0.
               mask = 0.

               downsamp_tod = 0.
               downsamp_mask = 0.
               downsamp_pix = 0
               downsamp_mask_idx = .true.
               downsamp_vec = 0. 

               call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)

               ! Get TOD after subtracting static zodi
               if (data(i)%tod%compressed_tod) then
                  call data(i)%tod%decompress_tod(scan, j, tod)
               else
                  tod = data(i)%tod%scans(scan)%d(j)%tod
               end if
               call get_s_tot_zodi(zodi_model, data(i)%tod, j, scan, s_static)
               tod = tod - data(i)%tod%scans(scan)%d(j)%gain * s_static

               do k = 1, data(i)%tod%scans(scan)%ntod
                  mask(k) = procmask_zodi(pix(k, 1))
                  if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
                  vec(:, k) = data(i)%tod%ind2vec(:, data(i)%tod%pix2ind(pix(k, 1)))
               end do

               where (mask > 0.5) 
                  mask = 1.
               elsewhere
                  mask = 0.
               end where

               do k = 1, 3
                  call data(i)%tod%downsample_tod(real(vec(k, :), sp), ext, downsamp_vec(k, :), width=box_width)
               end do

               do k = 0, ext(2)-padding
                  call vec2pix_ring(data(i)%tod%nside, real(downsamp_vec(:, k), dp), downsamp_pix(k))
               end do

               ! Downsample the mask
               call data(i)%tod%downsample_tod(mask, ext, downsamp_mask, width=box_width)

               ! Downsample the tods
               call data(i)%tod%downsample_tod(tod, ext, downsamp_tod, width=box_width)

               ! Construct the downsampled pix array (pick central pixel in each box)

               ! do k = 1, int(size(tod)/box_width) + 1
               !    l = floor(min(((k-1)*box_width + box_halfwidth), real(ntod, dp)))
               !    downsamp_pix(k) = pix(l, 1)
               ! end do

               where (downsamp_mask < 1.) downsamp_mask_idx = .false.

               ! Find upper bound and truncate the downsampled arrays by where they were padded and
               ! filter out the masked values. Store the result in the tod objects for use in sampling.
               data(i)%tod%scans(scan)%d(j)%downsamp_pix = pack(downsamp_pix(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               data(i)%tod%scans(scan)%d(j)%downsamp_tod = pack(downsamp_tod(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               if (j == 1) data(i)%tod%scans(scan)%downsamp_obs_time = pack(downsamp_obs_time(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               ! write timestreams to files
               ! call int2string(data(i)%tod%scanid(scan), scan_str)
               ! call open_hdf_file(trim(adjustl("/mn/stornext/u3/metins/dirbe/chains/chains_downsamp/dtod_"//scan_str//".h5")), tod_file, 'w')
               ! call write_hdf(tod_file, '/dtod', data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               ! call write_hdf(tod_file, '/vec', vec)
               ! call write_hdf(tod_file, '/dvec', downsamp_vec)
               ! call write_hdf(tod_file, '/pix', pix(:, 1))
               ! call write_hdf(tod_file, '/dpix', downsamp_pix)
               ! call close_hdf_file(tod_file)

               ! Allocate other downsampled quantities with same shape
               ndownsamp = size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
               if (.not. allocated(data(i)%tod%scans(scan)%d(j)%downsamp_scat)) then
                  allocate (data(i)%tod%scans(scan)%d(j)%downsamp_zodi(ndownsamp))
                  allocate (data(i)%tod%scans(scan)%d(j)%downsamp_scat(ndownsamp, zodi_model%n_comps))
                  allocate (data(i)%tod%scans(scan)%d(j)%downsamp_therm(ndownsamp, zodi_model%n_comps))

               ! Compute downsampled pointing in various coordinates
                  allocate(data(i)%tod%scans(scan)%d(j)%downsamp_point(ndownsamp,5))
               end if
               do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
!!$                  call pix2ang_ring(data(i)%tod%nside, data(i)%tod%scans(scan)%d(j)%downsamp_pix(k), lat, lon)
!!$                  phi = phi * 180.d0/pi
!!$                  lon = 90.d0 - lon * 180.d0/pi
!!$                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,1) = lon    ! Ecliptic longitude
!!$                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,2) = lat    ! Ecliptic latitude
!!$                  day          = data(i)%tod%scanid(scan)
!!$                  lambda_solar = (-80.598349 + 0.98564736 * day + 1.912 * cos(2.d0*pi/365.25 * (day-94.8))) * pi/180.d0
!!$                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,5) = acos(cos(lat*pi/180.d0) * cos(lon*pi/180.d0 - lambda_solar)) * 180.d0/pi ! Solar elongation
!!$                  call pix2vec_ring(data(i)%tod%nside, data(i)%tod%scans(scan)%d(j)%downsamp_pix(k), vec1)
!!$                  !vec0 = -0.5d0 * (data(i)%tod%scans(scan)%x0_obs + data(i)%tod%scans(scan)%x1_obs)
!!$                  !vec0 = vec0/sqrt(sum(vec0**2))
!!$                  !theta = dot_product(vec0, vec1)
!!$                  !data(i)%tod%scans(scan)%d(j)%downsamp_point(k,5) = acos(theta)*180.d0/pi                   ! Solar elongation
!!$                  vec1 = matmul(M_ecl2gal, vec1)
!!$                  call vec2ang(vec1, theta, phi)
!!$                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,3) = phi*180.d0/pi                 ! Galactic longitude
!!$                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,4) = (0.5d0*pi-theta)*180.d0/pi    ! Galactic latitude

                  call pix2ang_ring(data(i)%tod%nside, data(i)%tod%scans(scan)%d(j)%downsamp_pix(k), lat, lon)
                  lon = lon * 180.d0/pi
                  if (lon > 180.d0) lon = lon - 360.d0
                  lat = 90.d0 - lat * 180.d0/pi
                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,3) = lon    ! Galactic longitude
                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,4) = lat    ! Galactic latitude
                  call pix2vec_ring(data(i)%tod%nside, data(i)%tod%scans(scan)%d(j)%downsamp_pix(k), vec1)
                  vec1 = matmul(transpose(M_ecl2gal), vec1)
                  call vec2ang(vec1, lat, lon)
                  lon = lon * 180.d0/pi
                  if (lon > 180.d0) lon = lon - 360.d0
                  lat = 90.d0 - lat * 180.d0/pi
                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,1) = lon    ! Galactic longitude
                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,2) = lat    ! Galactic latitude
                  day          = data(i)%tod%scanid(scan)
                  lambda_solar = (-80.598349 + 0.98564736 * day + 1.912 * cos(2.d0*pi/365.25 * (day-94.8))) * pi/180.d0
                  data(i)%tod%scans(scan)%d(j)%downsamp_point(k,5) = acos(cos(lat*pi/180.d0) * cos(lon*pi/180.d0 - lambda_solar)) * 180.d0/pi ! Solar elongation
               end do
               
            end do
            deallocate (obs_time, downsamp_obs_time)
            deallocate (pix, psi, vec, flag, tod, s_static,  mask, downsamp_tod, downsamp_mask, downsamp_pix, downsamp_mask_idx, downsamp_vec)
         end do

         deallocate (m_buf, procmask_zodi)
      end do
   end subroutine
   

   subroutine precompute_lowres_zodi_lookups(cpar)
      ! Loop over each band with zodi and precompute lookup tables for lowres zodi caching
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, l, scan, n_lowres_obs, nobs_downsamp, nobs_lowres, pix_high, pix_low, ierr
      integer(i4b), allocatable :: pix2ind_highres(:), ind2pix_highres(:)
      real(dp), allocatable :: ind2vec_zodi_temp(:, :)
      real(dp) :: rotation_matrix(3, 3)

      call ecl_to_gal_rot_mat(rotation_matrix)

      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         allocate(pix2ind_highres(0:12*data(i)%tod%nside**2-1))         
         pix2ind_highres = 0
         
         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            do j = 1, data(i)%tod%ndet
               do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
                  pix2ind_highres(data(i)%tod%scans(scan)%d(j)%downsamp_pix(k)) = 1
               end do
            end do
         end do
         
         nobs_downsamp = count(pix2ind_highres == 1)
         
         allocate(ind2vec_zodi_temp(3, nobs_downsamp))
         allocate(ind2pix_highres(nobs_downsamp))

         j = 1
         do k = 0, 12*data(i)%tod%nside**2-1
            if (pix2ind_highres(k) == 1) then
               ind2pix_highres(j) = k
               j = j+1
            end if
         end do

         allocate(data(i)%tod%pix2ind_lowres(0:12*zodi_nside**2-1))
         data(i)%tod%pix2ind_lowres = 0
         ind2vec_zodi_temp = 0.

         j = 1
         do k = 1, nobs_downsamp
            pix_high = ind2pix_highres(k)
            pix_low = data(i)%tod%udgrade_pix_zodi(pix_high)
            if (data(i)%tod%pix2ind_lowres(pix_low) == 0) then
               data(i)%tod%pix2ind_lowres(pix_low) = j
               call pix2vec_ring(zodi_nside, pix_low, ind2vec_zodi_temp(:, j))
            end if
            j =  j + 1
         end do

         nobs_lowres = j - 1
         allocate(data(i)%tod%ind2vec_ecl_lowres(3, nobs_lowres))
         data(i)%tod%ind2vec_ecl_lowres = ind2vec_zodi_temp(:, 1:nobs_lowres)
         
         allocate(data(i)%tod%zodi_scat_cache_lowres(nobs_lowres, data(i)%tod%zodi_n_comps, data(i)%tod%ndet))
         allocate(data(i)%tod%zodi_therm_cache_lowres(nobs_lowres, data(i)%tod%zodi_n_comps, data(i)%tod%ndet))
         data(i)%tod%zodi_scat_cache_lowres = -1.d0
         data(i)%tod%zodi_therm_cache_lowres = -1.d0
         
         do k = 1, nobs_lowres
            data(i)%tod%ind2vec_ecl_lowres(:, k) = matmul(data(i)%tod%ind2vec_ecl_lowres(:, k), rotation_matrix)
         end do   
         deallocate(ind2vec_zodi_temp, pix2ind_highres, ind2pix_highres)
      end do
   end subroutine

   subroutine project_and_downsamp_sky(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, scan, ext(2), upper_bound, padding, ierr, ntod, nhorn, npix, ndet, nmaps, ndelta
      real(dp) :: box_width
      type(map_ptr), allocatable, dimension(:, :) :: sky_signal
      real(sp), allocatable, dimension(:) :: downsamp_sky, downsamp_mask
      logical(lgt), allocatable, dimension(:) :: downsamp_mask_idx
      real(sp), allocatable, dimension(:, :, :, :) :: map_sky
      real(dp), allocatable, dimension(:, :) :: m_buf
      type(hdf_file) :: tod_file
      integer(i4b), allocatable :: pix(:, :), psi(:, :), flag(:)
      real(sp), allocatable :: mask(:), sky(:)
      real(sp), allocatable, dimension(:) :: procmask_zodi

      if (cpar%myid == cpar%root) print *, "projecting downsampled sky"

      padding = 5
      ! For each zodi band, create downsampled residual time stream
      do i = 1, numband
         ! Only generate downsampled arrays for tod bands and bands where we have zodi
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         npix = 12*data(i)%map%info%nside**2
         ndelta = 1
         nmaps = data(i)%map%info%nmaps
         ndet = data(i)%tod%ndet
         nhorn = data(i)%tod%nhorn
         box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)

         allocate (sky_signal(data(i)%tod%ndet, ndelta))
         do j = 1, data(i)%tod%ndet
            call get_sky_signal(i, j, sky_signal(j, ndelta)%p, mono=.false.)
         end do

         allocate (map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
         call distribute_sky_maps(data(i)%tod, sky_signal, 1.e0, map_sky)

         allocate (m_buf(0:npix - 1, nmaps))
         allocate (procmask_zodi(0:npix - 1))
         call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf)
         procmask_zodi = m_buf(:, 1)

         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            ntod = data(i)%tod%scans(scan)%ntod

            do j = 1, ndet
               allocate (pix(ntod, nhorn))
               allocate (psi(ntod, nhorn))
               allocate (flag(ntod))
               allocate (mask(ntod))
               allocate (sky(ntod))
               pix = 0
               psi = 0
               flag = 0
               mask = 0.
               sky = 0.

               ! Project sky map down to downsampled timestream and cache
               call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)

               do k = 1, ntod
                  mask(k) = procmask_zodi(pix(k, 1))
                  if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
               end do
               where (mask > 0.) mask = 1. ! make sure mask is binary before downsampling
               
               do k = 1, ntod
                  sky(k) = map_sky(1, data(i)%tod%pix2ind(pix(k, 1)), j, 1)  !zodi is only temperature (for now)
               end do
               
               ! Get downsample shape
               call data(i)%tod%downsample_tod(sky, ext, width=box_width)

               allocate (downsamp_sky(ext(1):ext(2)))
               allocate (downsamp_mask(ext(1):ext(2)))
               allocate (downsamp_mask_idx(ext(1):ext(2)))
               downsamp_mask_idx = .true.

               call data(i)%tod%downsample_tod(sky, ext, downsamp_sky, width=box_width)

               call data(i)%tod%downsample_tod(mask, ext, downsamp_mask, width=box_width)
               where (downsamp_mask < 1.) downsamp_mask_idx = .false.

               data(i)%tod%scans(scan)%d(j)%downsamp_sky = pack(downsamp_sky(0:ext(2)-padding), downsamp_mask_idx(0:ext(2)-padding))
               deallocate(pix, psi, flag, mask, sky)
               deallocate(downsamp_sky, downsamp_mask_idx, downsamp_mask)
            end do

         end do
         deallocate (map_sky, sky_signal, m_buf, procmask_zodi)
      end do
   end subroutine

   subroutine compute_downsamp_zodi(cpar, model)
      type(comm_params), intent(in) :: cpar
      type(ZodiModel), intent(inout) :: model
      integer(i4b) :: i, j, scan, ierr, ndet

      if (cpar%myid == cpar%root) print *, "computing downsampled zodi"

      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle

            ndet = data(i)%tod%ndet
            do j = 1, ndet
               call get_zodi_emission(&
                   & tod=data(i)%tod, &
                   & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                   & scan=scan, &
                   & det=j, &
                   & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                   & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                   & model=model, &
                   & use_lowres_pointing=.true. &
               &)
               call get_s_zodi(i, &
                   & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                   & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                   & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                   &)
            end do
         end do
      end do
   end subroutine

   subroutine create_zodi_glitch_mask(cpar, handle)
     type(comm_params), intent(in) :: cpar
     type(planck_rng), intent(inout) :: handle
      integer(i4b) :: i, j, k, scan, ierr, non_glitch_size, thinstep, offset
      real(dp) :: box_width, rms, frac
      real(sp), allocatable :: res(:)
      real(sp), allocatable :: downsamp_scat_comp(:, :), downsamp_therm_comp(:, :)

      thinstep = nint(1.d0/cpar%zs_tod_thin_factor)
      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         do scan = 1, data(i)%tod%nscan
            offset    = rand_uni(handle) * thinstep
            box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
            do j = 1, data(i)%tod%ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

               ! Search for strong outliers
               res = data(i)%tod%scans(scan)%d(j)%downsamp_tod - data(i)%tod%scans(scan)%d(j)%gain * (data(i)%tod%scans(scan)%d(j)%downsamp_sky + data(i)%tod%scans(scan)%d(j)%downsamp_zodi)
               rms = sqrt(mean(real(res**2, dp)))
               data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask = abs(res) > 10. * real(rms, sp)
               frac = count(data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)/real(size(data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask),dp)
               if (frac > 0.01d0) then
                  write(*,*) 'Warning: Removing high fraction of glitches = ', frac
               end if

               ! Apply TOD thinning
               !allocate(data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask(data(i)%tod%scans(scan)%ntod))
               do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  if (mod(k+offset,thinstep) /= 0) then
                     data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask(k) = .true.
                  end if
               end do
            end do
         end do
      end do
   end subroutine

   subroutine apply_zodi_glitch_mask(cpar)
      type(comm_params), intent(in) :: cpar
      integer(i4b) :: i, j, k, scan, ierr, non_glitch_size
      real(sp), allocatable :: res(:)
      real(sp), allocatable :: downsamp_scat_comp(:, :), downsamp_therm_comp(:, :)
      integer(i4b), allocatable :: downsamp_pix(:)

      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         do scan = 1, data(i)%tod%nscan
            do j = 1, data(i)%tod%ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
               !if (size(data(i)%tod%scans(scan)%d(j)%downsamp_tod) < 10) cycle
               
               ! Remove flagged samples
               non_glitch_size = count(.not. data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)
               if (size(data(i)%tod%scans(scan)%d(j)%downsamp_pix) /= non_glitch_size) data(i)%tod%scans(scan)%d(j)%downsamp_pix = pack(data(i)%tod%scans(scan)%d(j)%downsamp_pix, .not. data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)
               if (size(data(i)%tod%scans(scan)%d(j)%downsamp_tod) /= non_glitch_size) data(i)%tod%scans(scan)%d(j)%downsamp_tod = pack(data(i)%tod%scans(scan)%d(j)%downsamp_tod, .not. data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)
               if (size(data(i)%tod%scans(scan)%d(j)%downsamp_sky) /= non_glitch_size) data(i)%tod%scans(scan)%d(j)%downsamp_sky = pack(data(i)%tod%scans(scan)%d(j)%downsamp_sky, .not. data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)
               if (size(data(i)%tod%scans(scan)%d(j)%downsamp_zodi) /= non_glitch_size) data(i)%tod%scans(scan)%d(j)%downsamp_zodi = pack(data(i)%tod%scans(scan)%d(j)%downsamp_zodi, .not. data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)

               ! pack doesnt work on multidimensional arrays so here we manually reallocate the zodi caches to the new sizes
               if (size(data(i)%tod%scans(scan)%d(j)%downsamp_therm,1) /= non_glitch_size) then
                  allocate(downsamp_scat_comp(non_glitch_size, zodi_model%n_comps))
                  allocate(downsamp_therm_comp(non_glitch_size, zodi_model%n_comps))
                  do k = 1, zodi_model%n_comps
                     downsamp_scat_comp(:, k) = pack(data(i)%tod%scans(scan)%d(j)%downsamp_scat(:, k), .not. data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)
                     downsamp_therm_comp(:, k) = pack(data(i)%tod%scans(scan)%d(j)%downsamp_therm(:, k), .not. data(i)%tod%scans(scan)%d(j)%zodi_glitch_mask)
                  end do
                  deallocate(data(i)%tod%scans(scan)%d(j)%downsamp_scat)
                  deallocate(data(i)%tod%scans(scan)%d(j)%downsamp_therm)
                  allocate(data(i)%tod%scans(scan)%d(j)%downsamp_scat(non_glitch_size, zodi_model%n_comps))
                  allocate(data(i)%tod%scans(scan)%d(j)%downsamp_therm(non_glitch_size, zodi_model%n_comps))
                  data(i)%tod%scans(scan)%d(j)%downsamp_therm = downsamp_therm_comp
                  data(i)%tod%scans(scan)%d(j)%downsamp_scat = downsamp_scat_comp
                  deallocate(downsamp_scat_comp, downsamp_therm_comp)
               end if


            end do
         end do
      end do
   end subroutine

   subroutine minimize_zodi_with_powell(cpar, iter, handle, samp_group)
      implicit none 
      type(comm_params), intent(in)      :: cpar
      integer(i4b),      intent(in)      :: iter
      type(planck_rng),  intent(inout)   :: handle
      integer(i4b),      intent(in)      :: samp_group

      logical(lgt) :: accept
      real(dp), allocatable :: theta(:), theta_new(:), theta_old(:), scale(:)
      real(dp), allocatable, dimension(:) :: theta_prev, chisq_prev
      integer(i4b) :: i, j, k, ind, ierr, flag, ntot, npar, unit
      real(dp) :: chisq_old, chisq_new
      character(len=6) :: iter_string
      character(len=6) :: sgroup
      character(len=512) :: filename

      if (cpar%myid == 0) print *, "minimizing zodi parameters with powell, samp_group =", samp_group
      
      npar = count(zodi_model%theta_stat(:,samp_group)==0)
      ntot = zodi_model%npar_tot
      allocate(theta_old(npar), theta_new(npar), theta(npar))
      allocate(theta_prev(npar), chisq_prev(numband), scale(npar))

      ! Initialize active monopoles
      do i = 1, numband
         ind = zodi_model%get_par_ind(mono_band=i)
         band_monopole(i) = get_monopole_amp(data(i)%label)
      end do

      if (cpar%myid_chain == 0) then
         call int2string(samp_group, sgroup)
         call int2string(iter, iter_string)
         filename = trim(cpar%outdir)//'/zodi_powell_sg'//sgroup//'_k'//iter_string//'.dat'
         unit     = getlun()
         open(unit, file=trim(filename), recl=10000)
         write(unit, '(a)', advance="no") "# chisq_red "
         do i = 1, zodi_model%npar_tot
            if (zodi_model%theta_stat(i,samp_group)==0) then
               write(unit, "(a,a)", advance="no") trim(adjustl(zodi_model%par_labels_full(i))), " "
            end if
         end do
         write(unit,*)
      end if

      ! Get chisq of old point
      scale = pack(zodi_model%theta_scale(:,1), zodi_model%theta_stat(:,samp_group)==0)
      call zodi_model%model_to_params(theta_old, samp_group)
!!$      if (cpar%myid == cpar%root) then
!!$         do i = 1, npar
!!$            write(*,*) i, theta_old(i)
!!$         end do
!!$      end if
!!$      call mpi_finalize(ierr)
!!$      stop

      ! Enforce priors; rms = 0
      call randomize_zodi_init(theta_old, samp_group, cpar, handle, rms=0.d0)
      
      theta_prev = 0.d0
      chisq_prev = 0.d0
      if (cpar%myid == cpar%root) then
         theta = theta_old/scale
         chisq_old  = lnL_zodi(theta)
      else
         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
         chisq_old  = lnL_zodi()
      end if
      if (cpar%zs_output_tod_res) then
         call mpi_finalize(ierr)
         stop
      end if

      ! Initialize new point
      theta_new = theta_old
      call randomize_zodi_init(theta_new, samp_group, cpar, handle)

      ! Rescale 
      theta = theta_new/scale
         
      if (cpar%myid == cpar%root) then
         ! Perform search
         call powell(theta, lnL_zodi, ierr, tolerance=1d-5)
         chisq_new = lnL_zodi(theta)
         flag = 0
         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)

         ! Apply approximate Metropolis rule, using reduced chisq instead of chisq
         if (chisq_new < chisq_old) then
            accept = .true. 
         else
            accept = rand_uni(handle) < exp(-0.5d0*(chisq_new-chisq_old)/0.02d0)
         end if
         if (accept) then
            ! Accept new point; update
            if (cpar%myid == cpar%root) write(*,fmt='(a,f8.2,a,f8.2)') 'Zodi sample accepted, chisq_new =', chisq_new, ', chisq_old = ', chisq_old
            theta_new = theta*scale ! Convert to physical units
         else
            ! Reject new solution, reset to previous solution
            if (cpar%myid == cpar%root) write(*,fmt='(a,f8.2,a,f8.2)') 'Zodi sample rejected, chisq_new =', chisq_new, ', chisq_old = ', chisq_old
            theta_new = theta_old
         end if
      else
         do while (.true.)
             call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
             if (flag == 1) then 
                 chisq_new = lnL_zodi()
             else
                 exit
              end if
           end do
          !call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)          
      end if
      if (cpar%myid_chain == 0) close(unit)
      
      ! Distribute final solution
      call mpi_bcast(theta_new, size(theta_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
      
      ! update model with final parameters
      call zodi_model%params_to_model(theta_new, samp_group)

      ! Update monopole for requested bands
      do i = 1, numband
         if (band_update_monopole(i,samp_group)) call set_monopole_amp(data(i)%label, band_monopole(i))
      end do

      deallocate(theta_old, theta_new, theta, theta_prev, chisq_prev, scale)
      
    contains
      
      function lnL_zodi(p)
      use healpix_types
      implicit none
      real(dp), dimension(:), intent(in), optional :: p
      real(dp)                                     :: lnL_zodi

      real(dp), allocatable :: theta(:)
      real(dp) :: chisq, chisq_tot, box_width, t1, t2, t3, t4, mono
      integer(i4b) :: i, j, k, scan, ntod, ndet, nscan, flag, ierr, ndof, ndof_tot
      logical(lgt) :: accept
      logical(lgt), dimension(numband) :: update_band
      character(len=4) :: scan_str
      type(hdf_file) :: tod_file

      call wall_time(t1)
      
      allocate(theta(npar))
      if (cpar%myid_chain == 0) then
         flag = 1
         call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
         theta = p*scale
      end if
      call mpi_bcast(theta, size(theta), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)
      
      ! Check which parameters have changed
      update_band = .true.
!!$      update_band = .false.
!!$      j = 0
!!$      do i = 1, zodi_model%npar_tot
!!$         if (zodi_model%theta_stat(i,samp_group) /= 0) cycle
!!$         j = j+1
!!$         if (theta(j) /= theta_prev(j)) then
!!$            !if (data(1)%tod%myid == 0) write(*,*) j, theta(j), theta_prev(j), zodi_model%theta2band(i)
!!$            if (zodi_model%theta2band(i) == 0) then
!!$               update_band = .true.
!!$               exit
!!$            else
!!$               update_band(zodi_model%theta2band(i)) = .true.
!!$            end if
!!$         end if
!!$      end do
      !if (data(1)%tod%myid == 0) write(*,*) data(1)%tod%myid, ' -- update =', update_band, all(theta == theta_prev)
      
      ! Check priors
      if (cpar%myid_chain == 0) then
         chisq_tot = get_chisq_priors(theta, samp_group)
         accept    = chisq_tot < 1.d30
      else
         chisq_tot = 0.d0
      end if
      call mpi_bcast(accept, 1, MPI_LOGICAL, 0, data(1)%tod%comm, ierr)

      if (.not. accept) then
         deallocate(theta)
         lnL_zodi = 1.d30
         return
      end if

      call zodi_model%params_to_model(theta, samp_group)
      if (cpar%myid_chain == 0) call print_zodi_model(theta, samp_group)
      
      ndof = 0
      do i = 1, numband
         if (data(i)%tod_type == "none") cycle
         if (.not. data(i)%tod%sample_zodi) cycle
         if (.not. zodi_model%sampgroup_active_band(i,samp_group)) cycle
         ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
         if (chisq_tot >= 1.d30) exit

         ndet = data(i)%tod%ndet
         nscan = data(i)%tod%nscan

         if (.not. update_band(i)) then
            !write(*,*) 'skipping', i, chisq_prev(i)
            if (cpar%myid_chain == 0) write(*,*) 'skipping band ', i, chisq_prev(i)
            chisq_tot = chisq_tot + chisq_prev(i)
            do scan = 1, nscan
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                  ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               end do
            end do
            cycle
         end if

         ! Get monopole
         mono = band_monopole(i) !get_monopole_amp(data(i)%label)
!!$         if (data(1)%tod%myid == 0) then
!!$            write(*,*) 'theta =', theta
!!$            write(*,*) 'mono =', mono
!!$         end if
         
         box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)

         ! Make sure that the zodi cache is cleared before each new band
         call data(i)%tod%clear_zodi_cache()

         ! Evaluate zodi model with newly proposed values for each band and calculate chisq
         chisq_prev(i) = 0.d0
         do scan = 1, nscan
            ! Skip scan if no accepted data
            do j = 1, ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) then
                  !write(*,*) '   Scan rejected:', data(i)%tod%scanid(scan), data(1)%tod%myid
                  cycle
               end if

               call wall_time(t3)
               call get_zodi_emission(&
                   & tod=data(i)%tod, &
                   & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                   & scan=scan, &
                   & det=j, &
                   & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                   & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                   & model=zodi_model, &
                   & use_lowres_pointing=.true. &
                   &)
               call wall_time(t4)
               !if (data(1)%tod%myid == 10) write(*,*) ' CPU1 = ', t4-t3
               call get_s_zodi(i, &
                   & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                   & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                   & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                   &)
               call wall_time(t3)
               !if (data(1)%tod%myid == 10) write(*,*) ' CPU2 = ', t3-t4

               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'tod  ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_tod))
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'sky  ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_sky))
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'zodi ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_zodi))
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'mono ', mono
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'noise', data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
               chisq = sum( &
                  & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                  &   - data(i)%tod%scans(scan)%d(j)%gain * &
                  &     (data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                  &   + data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                  &   + mono) &
                  & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
                  &)
               chisq_tot     = chisq_tot     + chisq
               chisq_prev(i) = chisq_prev(i) + chisq
               call wall_time(t4)
               !if (data(1)%tod%myid == 10) write(*,*) ' CPU3 = ', t4-t3
               
!!$               write(*,*) 'a', i, scan!, allocated(data(i)%tod%scans(scan)%d(j)%downsa1mp_tod)
!!$               write(*,*) 'b', j!, allocated(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               !write(*,*) 'do not remove -- memory corruption bug "fix"', ndof!, allocated(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               if (chisq_tot >= 1.d30) exit
               ! call int2string(data(i)%tod%scanid(scan), scan_str)
               ! call open_hdf_file(trim(adjustl("/mn/stornext/u3/metins/dirbe/chains/chains_downsamp/dtodlnl_"//scan_str//".h5")), tod_file, 'w')
               ! call write_hdf(tod_file, '/dtod', data(i)%tod%scans(scan)%d(j)%downsamp_tod)
               ! call write_hdf(tod_file, '/dzodi', data(i)%tod%scans(scan)%d(j)%downsamp_zodi)
               ! call write_hdf(tod_file, '/dsky', data(i)%tod%scans(scan)%d(j)%downsamp_sky)
               ! call write_hdf(tod_file, '/dpix', data(i)%tod%scans(scan)%d(j)%downsamp_pix)
               ! call close_hdf_file(tod_file)

               !if (.false. .and. data(1)%tod%myid == 0 .and. scan == 1) then
               if (cpar%zs_output_tod_res) then
                  call int2string(data(i)%tod%scanid(scan), scan_str)
                  !write(*,*) "scan = ", data(i)%tod%scanid(scan), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_tod)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_sky)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_zodi)), data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
                  open(58,file=trim(cpar%outdir)//'/todres_'//trim(data(i)%tod%freq)//'_'//scan_str//'.dat', recl=2048)
                  do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                     write(58,*) data(i)%tod%scans(scan)%d(j)%downsamp_point(k,:), data(i)%tod%scans(scan)%d(j)%downsamp_tod(k) &
                          &   - data(i)%tod%scans(scan)%d(j)%gain * &
                          &     (data(i)%tod%scans(scan)%d(j)%downsamp_sky(k) &
                          &   + data(i)%tod%scans(scan)%d(j)%downsamp_zodi(k) &
                          &   + mono)
                  end do
                  close(58)
               end if
               call wall_time(t3)
               !if (data(1)%tod%myid == 10) write(*,*) ' CPU4 = ', t3-t4

            end do
         end do
      end do

      !write(*,*) data(1)%tod%myid, ' -- recomputed =', chisq_prev
      
      ! Reduce chisq to root process
      call mpi_reduce(chisq_tot, chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)
      call mpi_reduce(ndof,      ndof_tot, 1, MPI_INTEGER, MPI_SUM, 0, data(1)%tod%comm, ierr)

      call wall_time(t4)
      !if (data(1)%tod%myid == 0) write(*,*) ' CPU5 = ', t4-t3
      
      if (cpar%myid_chain == 0) then
         lnL_zodi = chisq/ndof_tot
         call wall_time(t2)
         if (ndof_tot > 0) write(*,fmt='(a,e16.8,a,f10.4,a,f8.3)') "chisq_zodi = ", chisq, ", chisq_red = ", chisq/ndof_tot, ", time = ", t2-t1
         write(*,*)
         write(unit,*) chisq/ndof_tot, real(theta,sp)
      end if

      theta_prev = theta
   end function

 end subroutine minimize_zodi_with_powell


!!$   subroutine minimize_zodi_with_powell2(cpar, handle)
!!$      type(comm_params), intent(in) :: cpar
!!$      type(planck_rng),  intent(inout)   :: handle
!!$      logical(lgt), save :: first_call = .true.
!!$      logical(lgt) :: accept
!!$      real(dp), allocatable :: theta(:), theta_phys(:), theta_new(:), theta_final(:)
!!$      character(len=128), allocatable :: labels(:)
!!$      integer(i4b) :: i, j, k, ierr, flag, group_idx, end_idx, n_bands
!!$      real(dp) :: chisq_lnL, chisq_red
!!$      real(dp), allocatable, dimension(:) :: prior_vec_powell_min, prior_vec_powell_max, prior_vec_powell_type 
!!$      integer(i4b), allocatable :: indices(:)
!!$
!!$      if (allocated(param_vec_current)) then
!!$         ! Start at a new random point 
!!$         call randomize_zodi_init(cpar, handle)
!!$      else
!!$         allocate(param_vec_current(zodi_model%n_params))
!!$         chisq_red_current = 1.d30
!!$      end if
!!$
!!$      allocate(theta(zodi_model%n_params))
!!$      call zodi_model%model_to_params2(theta, labels)
!!$      if (cpar%myid == 0) print *, "minimizing zodi parameters with powell"
!!$      
!!$      ! filter out parameters for powell search
!!$      allocate(powell_included_params(zodi_model%n_params))
!!$      powell_included_params = .false.
!!$
!!$      ! Get parameters to sample from sampling groups
!!$      do group_idx = 1, cpar%zs_num_samp_groups
!!$         call parse_samp_group_strings(cpar%zs_samp_groups(group_idx), labels, indices)
!!$         do i = 1, size(indices)
!!$            powell_included_params(indices(i)) = .true.
!!$         end do
!!$      end do
!!$      prior_vec_powell_min = prior_vec(:, 1)
!!$      prior_vec_powell_max = prior_vec(:, 2)
!!$      prior_vec_powell_type = prior_vec(:, 3)
!!$
!!$      ! ! append emissivities and albedoes to the vector given to powell.
!!$      do i = 1, numband
!!$         if (trim(data(i)%tod_type) == 'none') cycle
!!$         if (.not. data(i)%tod%subtract_zodi) cycle
!!$         if (cpar%ds_zodi_reference_band(data(i)%id_abs)) cycle
!!$         theta = [theta, data(i)%tod%zodi_emissivity]
!!$         powell_included_params = [powell_included_params, [(.true. , j=1, size(data(i)%tod%zodi_emissivity))]]
!!$         prior_vec_powell_min = [prior_vec_powell_min, [(emissivity_prior(1) , j=1, size(data(i)%tod%zodi_emissivity))]]
!!$         prior_vec_powell_max = [prior_vec_powell_max, [(emissivity_prior(2) , j=1, size(data(i)%tod%zodi_emissivity))]]
!!$         prior_vec_powell_type = [prior_vec_powell_type, [(real(0, dp), j=1, size(data(i)%tod%zodi_emissivity))]]
!!$      end do
!!$      
!!$      do i = 1, numband
!!$         if (trim(data(i)%tod_type) == 'none') cycle
!!$         if (.not. data(i)%tod%subtract_zodi) cycle
!!$         if (.not. any(data(i)%tod%zodi_albedo > EPS)) cycle
!!$         theta = [theta, data(i)%tod%zodi_albedo]
!!$         powell_included_params = [powell_included_params, [(.true. , j=1, size(data(i)%tod%zodi_albedo))]]
!!$         prior_vec_powell_min = [prior_vec_powell_min, [(albedo_prior(1) , j=1, size(data(i)%tod%zodi_albedo))]]
!!$         prior_vec_powell_max = [prior_vec_powell_max, [(albedo_prior(2) , j=1, size(data(i)%tod%zodi_albedo))]]
!!$         prior_vec_powell_type = [prior_vec_powell_type, [(real(0, dp), j=1, size(data(i)%tod%zodi_albedo))]]
!!$      end do
!!$      
!!$      allocate(prior_vec_powell(size(prior_vec_powell_min), 3))
!!$      prior_vec_powell(:, 1) = prior_vec_powell_min
!!$      prior_vec_powell(:, 2) = prior_vec_powell_max
!!$      prior_vec_powell(:, 3) = prior_vec_powell_type
!!$
!!$      if (any(pack(theta, powell_included_params) == 0.)) then
!!$         do i = 1, size(theta)
!!$            write(*,*) i, theta(i)
!!$         end do
!!$         stop "theta_0 contains zeros in zodi powell sampling. Cant compute physical units due to theta/theta_0"
!!$      end if
!!$
!!$      allocate(theta_0(size(theta)))
!!$      theta_0 = theta
!!$         
!!$      allocate(theta_phys(count(powell_included_params)))
!!$      
!!$      if (cpar%myid == cpar%root) then
!!$         !filter out N_0 parameters and scale to physical units
!!$         theta_phys = pack(theta / theta_0, powell_included_params)
!!$         !call powell(theta_phys, lnL_zodi, ierr)
!!$         do k = 1, 2
!!$            call powell(theta_phys, lnL_zodi2, ierr, tolerance=1d-3)
!!$         end do
!!$         flag = 0
!!$         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
!!$         chisq_red = lnL_zodi2(theta_phys)
!!$      else
!!$         do while (.true.)
!!$             call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
!!$             if (flag == 1) then 
!!$                 chisq_lnL = lnL_zodi2()
!!$             else
!!$                 exit
!!$             end if
!!$          end do
!!$          call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
!!$          chisq_red = lnL_zodi2() 
!!$      end if
!!$
!!$      call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$      call mpi_bcast(theta_phys, size(theta_phys), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
!!$      call mpi_bcast(chisq_red, 1, MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
!!$      
!!$      allocate(theta_final(zodi_model%n_params))
!!$      theta_final = theta(1:zodi_model%n_params)
!!$
!!$      j = 1
!!$      do i = 1, size(theta)
!!$         if (powell_included_params(i)) then
!!$            theta(i) = theta_phys(j) * theta_0(i)
!!$            j = j + 1
!!$         end if
!!$      end do
!!$
!!$      !call min_max_param_vec_with_priors(theta, prior_vec_powell)
!!$      
!!$
!!$      j = 1
!!$      do i = 1, size(theta)
!!$         if (powell_included_params(i)) then
!!$            if (i <= zodi_model%n_params) then 
!!$               theta_final(i) = theta_phys(j) * theta_0(i)
!!$            end if
!!$            j = j + 1
!!$         end if
!!$      end do
!!$
!!$      ! Apply approximate Metropolis rule, using reduced chisq instead of chisq
!!$      if (cpar%myid == cpar%root) then
!!$         if (chisq_red_current == 1.d30) then
!!$            accept = .true.
!!$         else
!!$            accept = rand_uni(handle) < exp(-0.5d0*(chisq_red-chisq_red_current)/0.1d0)
!!$         end if
!!$      end if
!!$      call mpi_bcast(accept, 1, MPI_LOGICAL, cpar%root, cpar%comm_chain, ierr)
!!$      
!!$      if (accept) then
!!$         ! Accept new point; update
!!$         if (cpar%myid == cpar%root) write(*,fmt='(a,f8.2,a,f8.2)') 'Zodi sample accepted, chisq_new =', chisq_red, ', chisq_old = ', chisq_red_current
!!$         param_vec_current = theta_final
!!$         chisq_red_current = chisq_red
!!$      else
!!$         ! Reject new solution, reset to previous solution
!!$         if (cpar%myid == cpar%root) write(*,fmt='(a,f8.2,a,f8.2)') 'Zodi sample rejected, chisq_new =', chisq_red, ', chisq_old = ', chisq_red_current
!!$         call zodi_model%params_to_model2(param_vec_current)         
!!$         deallocate(prior_vec_powell, theta_0, powell_included_params)
!!$         return
!!$      end if
!!$      
!!$      ! update model with final parameters
!!$      call zodi_model%params_to_model2(theta_final)
!!$
!!$      ! update emissivities and albedos with final parameters
!!$      j = 0
!!$      do i = numband, 1, -1
!!$         if (data(i)%tod_type == "none") cycle
!!$         if (.not. data(i)%tod%sample_zodi) cycle
!!$         if (any(data(i)%tod%zodi_albedo > EPS)) then
!!$            j = j + 1
!!$            data(i)%tod%zodi_albedo = theta(size(theta)-(j*zodi_model%n_comps)+1:size(theta)-(j-1)*zodi_model%n_comps)
!!$         else 
!!$            data(i)%tod%zodi_albedo = data(i)%tod%zodi_albedo
!!$         end if 
!!$      end do
!!$      do i = numband, 1, -1
!!$         if (data(i)%tod_type == "none") cycle
!!$         if (.not. data(i)%tod%sample_zodi) cycle
!!$         if (.not. ref_band(data(i)%id_abs)) then
!!$            j = j + 1
!!$            data(i)%tod%zodi_emissivity = theta(size(theta)-(j*zodi_model%n_comps)+1:size(theta)-(j-1)*zodi_model%n_comps)
!!$         else
!!$            data(i)%tod%zodi_emissivity = 1.
!!$         end if
!!$      end do
!!$      deallocate(prior_vec_powell, theta_0, powell_included_params)
!!$   end subroutine
!!$
!!$   function lnL_zodi2(p)
!!$      use healpix_types
!!$      implicit none
!!$      real(dp), dimension(:), intent(in), optional :: p
!!$      real(dp) :: lnL_zodi2
!!$      type(ZodiModel) :: model
!!$
!!$      real(dp), allocatable :: theta(:), theta_phys(:), theta_full(:)
!!$      character(len=128), allocatable :: labels(:)
!!$      real(dp) :: chisq, chisq_tod, chisq_prior, box_width, t1, t2, t3, t4, prior
!!$      real(dp), allocatable, dimension(:), save :: theta_prev, chisq_prev
!!$      integer(i4b) :: i, j, k, scan, ntod, ndet, nscan, flag, ierr, n_bands, ndof, ndof_tot
!!$      logical(lgt), dimension(numband) :: update_band
!!$      character(len=4) :: scan_str
!!$      type(hdf_file) :: tod_file
!!$
!!$      call wall_time(t1)
!!$      model = zodi_model
!!$      allocate(theta_phys(count(powell_included_params)))
!!$      
!!$      if (data(1)%tod%myid == 0) then
!!$         flag = 1
!!$         call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
!!$         theta_phys = p
!!$      end if
!!$      call mpi_bcast(theta_phys, size(theta_phys), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)
!!$      
!!$
!!$      allocate(theta(model%n_params))
!!$      call model%model_to_params2(theta, labels)
!!$
!!$      allocate(theta_full(size(powell_included_params)))
!!$      theta_full(1:model%n_params) = theta
!!$      
!!$      j = 1
!!$      do i = 1, size(theta_full)
!!$         if (powell_included_params(i)) then
!!$            if (i <= model%n_params) then 
!!$               theta(i) = theta_phys(j) * theta_0(i)
!!$            end if
!!$            theta_full(i) = theta_phys(j) * theta_0(i)
!!$            j = j + 1
!!$         end if
!!$      end do
!!$
!!$      call model%params_to_model2(theta)
!!$
!!$      ! if (data(1)%tod%myid == 0) call print_zodi_model(theta, labels)
!!$
!!$
!!$      j = 0
!!$      do i = numband, 1, -1
!!$         if (data(i)%tod_type == "none") cycle
!!$         if (.not. data(i)%tod%sample_zodi) cycle
!!$         if (any(data(i)%tod%zodi_albedo > EPS)) then
!!$            j = j + 1
!!$            powell_albedo(i, :) = theta_full(size(theta_full)-(j*zodi_model%n_comps)+1:size(theta_full)-(j-1)*zodi_model%n_comps)
!!$         else 
!!$            powell_albedo(i, :) = data(i)%tod%zodi_albedo
!!$         end if 
!!$      end do
!!$      do i = numband, 1, -1
!!$         if (data(i)%tod_type == "none") cycle
!!$         if (.not. data(i)%tod%sample_zodi) cycle
!!$         if (.not. ref_band(data(i)%id_abs)) then
!!$            j = j + 1
!!$            powell_emissivity(i, :) = theta_full(size(theta_full)-(j*zodi_model%n_comps)+1:size(theta_full)-(j-1)*zodi_model%n_comps)
!!$         else
!!$            powell_emissivity(i, :) = 1.
!!$         end if
!!$      end do
!!$
!!$      ! rescale parameters 
!!$      call model%model_to_params2(theta)
!!$      theta_full(1:model%n_params) = theta
!!$
!!$      ! Check which parameters have changed
!!$      update_band = .false.
!!$      if (.not. allocated(theta_prev)) then
!!$         allocate(theta_prev(model%n_params), chisq_prev(numband))
!!$         theta_prev  = theta_full
!!$         chisq_prev  = 0.d0
!!$         update_band = .true.
!!$      else
!!$         do i = 1, model%n_params
!!$            if (theta_full(i) /= theta_prev(i)) then
!!$               if (param2band(i) == 0) then
!!$                  update_band = .true.
!!$                  exit
!!$               else
!!$                  update_band(param2band(i)) = .true.
!!$               end if
!!$            end if
!!$         end do
!!$      end if
!!$      
!!$      ! Check priors
!!$      if (data(1)%tod%myid == 0) then
!!$         call print_zodi_model(theta, labels, chisq)
!!$         chisq_tod = get_chisq_priors(theta_full, prior_vec_powell)
!!$         prior     = chisq_tod
!!$      else
!!$         chisq_tod = 0.d0
!!$      end if
!!$      call mpi_bcast(prior, 1, MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)
!!$
!!$      if (prior >= 1.d30) then
!!$         lnL_zodi2 = 1.d30
!!$         return
!!$      end if
!!$      
!!$      ndof = 0
!!$      do i = 1, numband
!!$         if (data(i)%tod_type == "none") cycle
!!$         if (.not. data(i)%tod%sample_zodi) cycle
!!$         ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
!!$         if (chisq_tod >= 1.d30) exit
!!$
!!$         ndet = data(i)%tod%ndet
!!$         nscan = data(i)%tod%nscan
!!$
!!$         if (.not. update_band(i)) then
!!$            chisq_tod = chisq_tod + chisq_prev(i)
!!$            do scan = 1, nscan
!!$               do j = 1, ndet
!!$                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
!!$                  ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
!!$               end do
!!$            end do
!!$            cycle
!!$         end if
!!$         
!!$         box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)
!!$
!!$         ! Make sure that the zodi cache is cleared before each new band
!!$         call data(i)%tod%clear_zodi_cache()
!!$
!!$         ! Evaluate zodi model with newly proposed values for each band and calculate chisq
!!$         do scan = 1, nscan
!!$            ! Skip scan if no accepted data
!!$            do j = 1, ndet
!!$               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
!!$
!!$               call wall_time(t3)
!!$               call get_zodi_emission(&
!!$                   & tod=data(i)%tod, &
!!$                   & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
!!$                   & scan=scan, &
!!$                   & det=j, &
!!$                   & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
!!$                   & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
!!$                   & model=model, &
!!$                   & use_lowres_pointing=.true. &
!!$                   &)
!!$               call wall_time(t4)
!!$               !if (data(1)%tod%myid == 0) write(*,*) ' CPU1 = ', t4-t3
!!$               call get_s_zodi(&
!!$                   & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
!!$                   & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
!!$                   & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi, &
!!$                   & emissivity=powell_emissivity(i, :), &
!!$                   & albedo=powell_albedo(i, :) &
!!$                   &)
!!$               call wall_time(t3)
!!$               !if (data(1)%tod%myid == 0) write(*,*) ' CPU2 = ', t3-t4
!!$               
!!$               chisq = sum( &
!!$                  & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
!!$                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky &
!!$                  &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
!!$                  & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &
!!$                  &)
!!$               chisq_tod = chisq_tod + chisq
!!$               chisq_prev(i) = chisq
!!$               call wall_time(t4)
!!$               !if (data(1)%tod%myid == 0) write(*,*) ' CPU3 = ', t4-t3
!!$
!!$
!!$               
!!$               ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
!!$               if (chisq_tod >= 1.d30) exit
!!$               ! call int2string(data(i)%tod%scanid(scan), scan_str)
!!$               ! call open_hdf_file(trim(adjustl("/mn/stornext/u3/metins/dirbe/chains/chains_downsamp/dtodlnl_"//scan_str//".h5")), tod_file, 'w')
!!$               ! call write_hdf(tod_file, '/dtod', data(i)%tod%scans(scan)%d(j)%downsamp_tod)
!!$               ! call write_hdf(tod_file, '/dzodi', data(i)%tod%scans(scan)%d(j)%downsamp_zodi)
!!$               ! call write_hdf(tod_file, '/dsky', data(i)%tod%scans(scan)%d(j)%downsamp_sky)
!!$               ! call write_hdf(tod_file, '/dpix', data(i)%tod%scans(scan)%d(j)%downsamp_pix)
!!$               ! call close_hdf_file(tod_file)
!!$
!!$               if (data(1)%tod%myid == 0 .and. scan == 1) then
!!$                  !write(*,*) "scan = ", data(i)%tod%scanid(scan), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_tod)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_sky)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_zodi)), data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
!!$                  open(58,file='res'//trim(data(i)%tod%freq)//'.dat')
!!$                  do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
!!$                     write(58,*) data(i)%tod%scans(scan)%d(j)%downsamp_tod(k) &
!!$                          &   - data(i)%tod%scans(scan)%d(j)%downsamp_sky(k) &
!!$                          &   - data(i)%tod%scans(scan)%d(j)%downsamp_zodi(k)
!!$                  end do
!!$                  close(58)
!!$               end if
!!$               call wall_time(t3)
!!$               !if (data(1)%tod%myid == 0) write(*,*) ' CPU4 = ', t3-t4
!!$
!!$            end do
!!$         end do
!!$      end do
!!$      ! call mpi_barrier(MPI_COMM_WORLD, ierr)
!!$      ! stop
!!$
!!$      ! Reduce chisq to root process
!!$      call mpi_reduce(chisq_tod, chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)
!!$      call mpi_reduce(ndof,      ndof_tot, 1, MPI_INTEGER, MPI_SUM, 0, data(1)%tod%comm, ierr)
!!$
!!$      call wall_time(t4)
!!$      !if (data(1)%tod%myid == 0) write(*,*) ' CPU5 = ', t4-t3
!!$
!!$      if (data(1)%tod%myid == 0) then
!!$         lnL_zodi2 = chisq/ndof_tot
!!$         call wall_time(t2)
!!$         if (ndof_tot > 0) write(*,fmt='(a,e16.8,a,f10.4,a,f8.3)') "chisq_zodi = ", chisq, ", chisq_red = ", chisq/ndof_tot, ", time = ", t2-t1
!!$         write(*,*)
!!$         ! call print_zodi_model(theta, labels, chisq)
!!$      end if
!!$   end function

   subroutine get_s_zodi_with_n0(s_therm, s_scat, s_zodi, emissivity, albedo, n_0_comp_ratio, comp)
      real(sp), dimension(:, :), intent(in) :: s_scat, s_therm
      real(sp), dimension(:), intent(inout) :: s_zodi
      real(dp), dimension(:), intent(in) :: emissivity, albedo, n_0_comp_ratio
      integer(i4b) :: comp
      integer(i4b) :: i

      s_zodi = 0.
      do i = 1, size(s_therm, dim= 2)
         s_zodi = s_zodi + (s_scat(:, i) * albedo(i) + (1. - albedo(i)) * emissivity(i) * s_therm(:, i)) * n_0_comp_ratio(i)
      end do
   end subroutine get_s_zodi_with_n0

   subroutine parse_samp_group_strings(samp_group_str, param_labels, param_indices)
      character(len=*), intent(in) :: samp_group_str
      character(len=*), intent(in) :: param_labels(:)
      integer(i4b), allocatable, intent(inout) :: param_indices(:)
      character(len=128) :: tokens(100), comp_param(2), label, param_label_tokens(10)
      character(len=128), allocatable :: tokens_trunc(:)
      integer(i4b) :: i, j, n_params
      logical(lgt) :: found


      if (allocated(param_indices)) deallocate(param_indices)

      call get_tokens(samp_group_str, ',', tokens, n_params) 
      tokens_trunc = tokens(1:n_params)
      allocate(param_indices(n_params))
      param_indices = -1
      do i = 1, size(tokens_trunc)
         call get_tokens(tokens_trunc(i), ':', comp_param) 
         call toupper(comp_param(1))
         call toupper(comp_param(2))
         if (trim(adjustl(comp_param(2))) == "ALL") then
            do j = 1, size(param_labels)
               call get_tokens(param_labels(j), "_", param_label_tokens) 
               call toupper(param_label_tokens(1))
               if (trim(adjustl(param_label_tokens(2))) == "N_0") then ! dont add N_0 to sampling groups with comp:all
                  cycle
               end if
               if (trim(adjustl(comp_param(1))) == param_label_tokens(1)) then
                  if (.not. any(param_indices == j)) param_indices = [param_indices , j]
               end if
            end do
         else
            found = .false.
            label = trim(adjustl(comp_param(1)))//"_"//trim(adjustl(comp_param(2)))
            do j = 1, size(param_labels)
               if (label == param_labels(j)) then
                  if (.not. any(param_indices == j)) param_indices(i) = j
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) then
               print *, "Error: invalid zodi sampling group parameter label :" // trim(adjustl(label)) 
               stop
            end if 
         end if
      end do
      param_indices = pack(param_indices, param_indices > 0)
   end subroutine


!!$   subroutine zodi_model_to_ascii(cpar, model, filename, overwrite)
!!$      ! Dumps the zodi model to an ascii file on the format {COMP}_{PARAM} = {VALUE}.
!!$      class(ZodiModel), target, intent(in) :: model
!!$      type(comm_params), intent(in) :: cpar
!!$      character(len=*), intent(in) :: filename
!!$      logical(lgt), intent(in), optional :: overwrite
!!$
!!$      integer(i4b) :: io, i, j, running_idx
!!$      logical(lgt) :: exists, overwrite_
!!$      real(dp), allocatable :: params(:)
!!$      integer(i4b), allocatable :: comp_switch_indices(:)
!!$      character(len=128), allocatable :: labels(:)
!!$      character(len=512) :: concatenated_string, val
!!$
!!$      if (present(overwrite)) then
!!$         overwrite_ = overwrite
!!$      else
!!$         overwrite_ = .false.
!!$      end if
!!$
!!$      if (cpar%myid_chain /= cpar%root) return
!!$      inquire(file=trim(adjustl(filename)), exist=exists)
!!$      if (exists .and. (.not. overwrite_)) then
!!$         print *, "zodi asciifile: " // trim(adjustl(filename)) // " exists and overwrite = .false."
!!$         stop
!!$      end if
!!$
!!$      open(newunit=io, file=trim(adjustl(filename)), action="write")
!!$      allocate(params(model%n_params))
!!$      call model%model_to_params2(params, labels=labels)
!!$
!!$      allocate(comp_switch_indices(model%n_comps))
!!$
!!$      running_idx = 0
!!$      do i = 1, model%n_comps
!!$         running_idx = running_idx + size(model%comps(i)%labels)
!!$         comp_switch_indices(i) = running_idx
!!$      end do
!!$
!!$      do i = 1, model%n_params
!!$         if (any(comp_switch_indices == i)) then
!!$               write(io, fmt='(a, T25, a, ES12.5, a)') trim(adjustl(labels(i))), "= ", params(i), new_line('a')
!!$            else
!!$               write(io, fmt='(a, T25, a, ES12.5)') trim(adjustl(labels(i))), "= ", params(i)
!!$         end if
!!$      end do
!!$
!!$      write(io, fmt='(a)') ''
!!$      do i = 1, numband
!!$         if (trim(data(i)%tod_type) == 'none') cycle
!!$         if (.not. data(i)%tod%subtract_zodi) cycle
!!$         concatenated_string = ""
!!$         do j = 1, model%n_comps
!!$            write(val, fmt='(ES12.5)') model%comps(j)%c%emissivity(i) 
!!$            concatenated_string = trim(adjustl(concatenated_string)) // "," // trim(adjustl(val))
!!$         end do
!!$         write(io, fmt='(a, T25, a, a)') trim(adjustl("EMISSIVITY_"//trim(adjustl(data(i)%tod%freq)))), "= ", trim(adjustl(concatenated_string(2:)))
!!$      end do
!!$
!!$      write(io, fmt='(a)') ''
!!$      do i = 1, numband
!!$         if (trim(data(i)%tod_type) == 'none') cycle
!!$         if (.not. data(i)%tod%subtract_zodi) cycle
!!$         concatenated_string = ""
!!$         do j = 1, model%n_comps
!!$            write(val, fmt='(ES12.5)') model%comps(j)%c%albedo(i)
!!$            concatenated_string = trim(adjustl(concatenated_string)) // "," // trim(adjustl(val))
!!$         end do
!!$         write(io, fmt='(a, T25, a, a)') trim(adjustl("ALBEDO_"//trim(adjustl(data(i)%tod%freq)))), "= ", trim(adjustl(concatenated_string(2:)))
!!$      end do
!!$
!!$      close(io)
!!$   end subroutine
!!$
!!$   subroutine ascii_to_zodi_model(cpar, model, filename)
!!$      ! Reads in and updates the zodi model from an ascii file on the format {COMP}_{PARAM} = {VALUE}.
!!$      class(ZodiModel), target, intent(inout) :: model
!!$      type(comm_params), intent(in) :: cpar
!!$      character(len=*), intent(in) :: filename
!!$      type(hash_tbl_sll) :: htbl
!!$
!!$      integer(i4b) :: i, j, io, io_status, ierr, n_comps
!!$      logical(lgt) :: exists
!!$      character(len=512) :: key, val, line
!!$      character(len=128), allocatable :: labels(:)
!!$      characteR(len=128) :: toks(100)
!!$      characteR(len=512) :: concatenated_string
!!$      real(dp), allocatable :: params(:)
!!$
!!$      allocate(params(model%n_params))
!!$      !if (cpar%myid_chain == cpar%root) then
!!$         inquire(file=trim(adjustl(filename)), exist=exists)
!!$         if (.not. exists) then
!!$            print *, "zodi asciifile: " // trim(adjustl(filename)) // " does not exist"
!!$            stop
!!$         end if
!!$         
!!$         call init_hash_tbl_sll(htbl, tbl_len=500)
!!$         
!!$         open(newunit=io, file=trim(adjustl(filename)), action="read")
!!$         io_status = 0
!!$         do while (io_status == 0)
!!$            read(io, "(a)", iostat=io_status) line
!!$            if (io_status == 0 .and. line /= "") then
!!$               j = index(line, "=")
!!$               if (j == 0) then
!!$                  print *, "Error: invalid line in ascii file: ", trim(adjustl(line))
!!$                  close(io)
!!$                  stop
!!$               end if
!!$
!!$               key = trim(adjustl(line(:j-1)))
!!$               val = trim(adjustl(line(j+1:)))
!!$               call tolower(key)
!!$               call put_hash_tbl_sll(htbl, trim(adjustl(key)), trim(adjustl(val))) 
!!$            end if
!!$         end do
!!$         close(io)
!!$
!!$         call model%model_to_params2(params, labels)
!!$         params = 0.
!!$         if (size(labels) /= size(params)) stop "Error: size of labels and params do not match"
!!$         do i = 1, size(labels)
!!$            call get_parameter_hashtable(htbl, labels(i), par_dp=params(i))
!!$         end do
!!$      !end if
!!$
!!$      !call mpi_bcast(params, size(params), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
!!$      call model%params_to_model2(params)
!!$
!!$      do i = 1, numband
!!$         if (trim(data(i)%tod_type) == 'none') cycle
!!$         if (.not. data(i)%tod%subtract_zodi) cycle
!!$         !if (cpar%myid == 0) then
!!$            call get_parameter_hashtable(htbl, trim(adjustl("EMISSIVITY_"//trim(adjustl(data(i)%tod%freq)))), par_string=concatenated_string)
!!$            call get_tokens(trim(adjustl(concatenated_string)), ',', toks, n_comps)
!!$            if (n_comps /= model%n_comps) stop "Error: number of components in ascii file does not match model emissivity"
!!$            do j = 1, n_comps
!!$               read(toks(j), *) model%comps(j)%c%emissivity(i)
!!$            end do
!!$
!!$            call get_parameter_hashtable(htbl, trim(adjustl("ALBEDO_"//trim(adjustl(data(i)%tod%freq)))), par_string=concatenated_string)
!!$            call get_tokens(trim(adjustl(concatenated_string)), ',', toks, n_comps)
!!$            if (n_comps /= model%n_comps) stop "Error: number of components in ascii file does not match model albedo"
!!$            do j = 1, n_comps
!!$               read(toks(j), *) model%comps(j)%c%albedo(i)
!!$            end do
!!$         !end if
!!$         !call mpi_bcast(data(i)%tod%zodi_emissivity, size(data(i)%tod%zodi_emissivity), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
!!$         !call mpi_bcast(data(i)%tod%zodi_albedo, size(data(i)%tod%zodi_albedo), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
!!$      end do
!!$   end subroutine
!!$
!!$
!!$   subroutine print_zodi_model(theta, samp_group)
!!$     implicit none
!!$     real(dp),     allocatable, intent(in) :: theta(:)
!!$     integer(i4b),              intent(in) :: samp_group
!!$      
!!$      integer(i4b) :: i, j, k, l, n, idx, n_cols, col, comp, n_params
!!$      logical(lgt) :: newline
!!$      
!!$      n_cols = 5
!!$      n_params = size(theta)
!!$
!!$      ! General parameters
!!$      k = 1
!!$      if (any(zodi_model%theta_stat(1:zodi_model%n_general_params,samp_group)==0)) then
!!$         write(*, "(a)", advance="no") 'General: '
!!$         col = 1
!!$         do i = 1, zodi_model%n_general_params
!!$            if (zodi_model%theta_stat(i,samp_group)==0) then
!!$               if (col > n_cols .or. i == zodi_model%n_general_params) then
!!$                  write(*, "(a,a,g0.4,a)") trim(adjustl(zodi_model%par_labels(i))), "=", theta(k), ",  "
!!$                  col = 1
!!$               else
!!$                  write(*, "(a,a,g0.4,a)", advance="no") trim(adjustl(zodi_model%par_labels(i))), "=", theta(k), ",  "
!!$                  col = col+1
!!$               end if
!!$               k = k+1
!!$            end if
!!$         end do
!!$      end if
!!$
!!$      ! Component parameters
!!$      do j = 1, zodi_model%n_comps
!!$         idx = zodi_model%comps(j)%start_ind
!!$         n   = zodi_model%comps(j)%npar + 2*numband
!!$         if (all(zodi_model%theta_stat(idx:idx+n-1,samp_group)/=0)) cycle
!!$         write(*, "(a,a)", advance="no") trim(adjustl(zodi_model%comp_labels(j))),': '
!!$         col = 1
!!$         do i = idx, idx+n-1
!!$            newline = (i==idx+n-1) .or. col == n_cols-1
!!$            if (zodi_model%theta_stat(i,samp_group)==0) then
!!$               write(*, "(a,a,a,g0.4,a)", advance="no") "  ", trim(adjustl(zodi_model%par_labels(i))), "=", theta(k), ",  "
!!$               k   = k+1
!!$               col = col+1
!!$            end if
!!$            if (newline .and. col>1) then
!!$               write(*,*)
!!$               col = 1
!!$            end if
!!$         end do
!!$      end do
!!$
!!$      ! Monopoles
!!$      col = 1
!!$      write(*, "(a)", advance="no") 'Mono : '
!!$      do j = 1, numband
!!$         idx = zodi_model%npar_tot - numband + j
!!$         if (zodi_model%theta_stat(idx,samp_group)==0) then
!!$            write(*, "(a,a,a,g0.4,a)", advance="no") "  ", trim(adjustl(data(j)%label)), "=", theta(k), ",  "
!!$            k   = k+1
!!$            col = col+1
!!$         end if
!!$         if ((col == n_cols-1 .and. col>1) .or. j == numband) then
!!$            write(*,*)
!!$            col = 1
!!$         end if
!!$      end do
!!$      
!!$   end subroutine
 
   subroutine min_max_param_vec_with_priors(params, priors)
      real(dp), intent(inout) :: params(:)
      real(dp), intent(in) :: priors(:, :)
      integer(i4b) :: i

      if (size(params) /= size(priors,1)) then
         write(*,*) "Error: params and priors must have the same size", size(params), size(priors,1)
         stop
      end if
      do i = 1, size(params)
         if (priors(i, 3) == 1) cycle ! (0 is uniform, 1 is gaussian)
         if (params(i) < priors(i, 1)) then
            params(i) = priors(i, 1)
         else if (params(i) > priors(i, 2)) then
            params(i) = priors(i, 2)
         end if
      end do

   end subroutine

   subroutine randomize_zodi_init(x, samp_group, cpar, handle, rms)
     implicit none
     real(dp),          dimension(:), intent(inout) :: x
     integer(i4b),                    intent(in)    :: samp_group
     type(comm_params),               intent(in)    :: cpar
     type(planck_rng),                intent(inout) :: handle
     real(dp),                        intent(in), optional :: rms
     
     integer(i4b) :: i, j, ierr
     real(dp)     :: eps

     eps = cpar%zs_randomize_rms; if (present(rms)) eps = rms
     if (eps < 0.d0) return     

     if (cpar%myid == 0) then
        j = 0
        do i = 1, zodi_model%npar_tot
           if (zodi_model%theta_stat(i,samp_group) == 0) then
              j = j+1
              x(j) = x(j) + eps*zodi_model%theta_scale(i,2)*rand_gauss(handle)
              !x(j) = x(j) * (1.d0 + eps*rand_gauss(handle))
              x(j) = max(x(j), zodi_model%theta_prior(1,i))
              x(j) = min(x(j), zodi_model%theta_prior(2,i))
           end if
        end do
     end if
     call mpi_bcast(x, size(x), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)

   end subroutine randomize_zodi_init
   
   subroutine sample_static_zodi_map(cpar, handle, map_id)
     implicit none
      type(comm_params), intent(inout) :: cpar
      type(planck_rng),  intent(inout) :: handle
      character(len=*),  intent(in)    :: map_id

      integer(i4b) :: band, i, j, k, ndet, scan, nscan, npix, nmaps, p, ierr, ntod, nhorn, npix_band, ncomp, nactive
      real(dp)     :: res, w, vec(3), elon, amp
      character(len=512) :: model
      !type(comm_scandata) :: sd
      real(dp),      allocatable, dimension(:)       :: A, b
      real(sp),      allocatable, dimension(:)       :: s_sky
      real(sp),      allocatable, dimension(:,:)     :: s_scat, s_therm
      real(sp),      allocatable, dimension(:)       :: s_zodi
      real(sp),      allocatable, dimension(:,:,:,:) :: map_sky
      type(map_ptr), allocatable, dimension(:,:)     :: sky_signal
      real(sp),      allocatable, dimension(:)       :: tod, mask, procmask
      real(dp),      allocatable, dimension(:,:)     :: m_buf
      integer(i4b),  allocatable, dimension(:,:)     :: pix, psi
      integer(i4b),  allocatable, dimension(:)       :: flag
      character(len=128), dimension(100)             :: active_bands
      logical(lgt), allocatable, dimension(:)        :: active

      if (cpar%myid == 0) then
         if (trim(map_id) == 'solar') then
            write(*,*) '   Sampling solar centric model maps'
         else if (trim(map_id) == 'moon') then
            write(*,*) '   Sampling Moon centric model maps'
         else if (trim(map_id) == 'earth') then
            write(*,*) '   Sampling Earth centric model maps'
         else
            write(*,*) '   Unknown static map type = ', trim(map_id)
            stop
         end if
      end if

      ncomp = zodi_model%n_comps
      nmaps = 1
      
      do band = 1, numband
         if (trim(data(band)%tod_type) == 'none') cycle
         if (trim(map_id) == 'solar') then
            model = cpar%ds_tod_solar_model(data(band)%tod%band)
         else if (trim(map_id) == 'moon') then
            model = cpar%ds_tod_moon_model(data(band)%tod%band)
         else if (trim(map_id) == 'earth') then
            model = cpar%ds_tod_earth_model(data(band)%tod%band)
         end if
         if (trim(model) == 'none') cycle
         if (model(1:1) == '>') then
            do i = 1, numband
               if (trim(data(i)%label) == trim(model(2:))) then
                  if (trim(map_id) == 'solar') then
                     data(band)%tod%map_solar => data(i)%tod%map_solar
                  else if (trim(map_id) == 'moon') then
                     data(band)%tod%map_moon => data(i)%tod%map_moon
                  else if (trim(map_id) == 'earth') then
                     data(band)%tod%map_earth => data(i)%tod%map_earth
                  end if
                  exit
               end if
            end do
            cycle
         end if

         if (trim(map_id) == 'earth') then
            npix  = NBIN_EARTH_ELON
         else
            npix  = 12*data(band)%info%nside**2
         end if

         ! Allocate temporary map structures
         allocate(A(0:npix-1), b(0:npix-1))
         A = 0.d0; b = 0.d0

         ! Find active bands
         allocate(active(numband))
         call get_tokens(model, ',', active_bands, nactive)
         active = .false.
         do j = 1, nactive
            do i = 1, numband
               if (trim(active_bands(j)) == trim(data(i)%label)) then
                  active(i) = .true.
                  exit
               end if
            end do
         end do
      
         ! Add up contributions from all active bands
         do i = 1, numband
            if (data(i)%tod_type == "none") cycle
            if (.not. data(i)%tod%sample_zodi) cycle
            if (.not. active(i)) cycle
            ndet      = data(i)%tod%ndet
            nscan     = data(i)%tod%nscan
            nhorn     = data(i)%tod%nhorn
            npix_band = 12*data(i)%map%info%nside**2

            ! Get and distribute sky signal
            allocate(sky_signal(data(i)%tod%ndet,1))
            do j = 1, data(i)%tod%ndet
               call get_sky_signal(i, j, sky_signal(j,1)%p, mono=.true.)
            end do
            allocate (map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
            call distribute_sky_maps(data(i)%tod, sky_signal, 1.e0, map_sky)
            
            ! Initialize frequency-specific mask
            allocate(m_buf(0:npix_band-1, nmaps), procmask(0:npix_band-1))
            !call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:, 1)
            call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask = m_buf(:, 1)
            deallocate(m_buf)
         
            do scan = 1, nscan
               ntod = data(i)%tod%scans(scan)%ntod
               allocate(s_scat(ntod,ncomp), s_therm(ntod,ncomp), s_zodi(ntod), s_sky(ntod))
            
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                  
                  ! Get data and pointing
                  allocate(pix(ntod, nhorn), psi(ntod, nhorn), flag(ntod), tod(ntod), mask(ntod))
                  if (data(i)%tod%compressed_tod) then
                     call data(i)%tod%decompress_tod(scan, j, tod)
                  else
                     tod = data(i)%tod%scans(scan)%d(j)%tod
                  end if
                  
                  ! Set up mask; remove flagged samples and foreground contaminated regions
                  call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)
                  do k = 1, data(i)%tod%scans(scan)%ntod
                     mask(k) = procmask(pix(k, 1))
                     if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
                  end do
                  where (mask > 0.5) 
                     mask = 1.
                  elsewhere
                     mask = 0.
                  end where

                  ! Compute non-stationary zodi TOD; exclude current static component
                  call get_s_tot_zodi(zodi_model, data(i)%tod, j, scan, s_zodi, pix_dynamic=pix, exclude_static=map_id)
                  
                  ! Add residual to mapmaking equation in solar centric coordinates
                  w  = 1.d0/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0**2
                  amp = 1.d0 !zodi_model%amp_static(i)
                  do k = 1, ntod
                     if (mask(k) == 0) cycle
                     if (trim(map_id) == 'solar') then
                        p = data(i)%tod%scans(scan)%d(j)%pix_sol(k,1)
                     else if (trim(map_id) == 'moon') then
                        p = data(i)%tod%scans(scan)%d(j)%pix_moon(k,1)
                     else if (trim(map_id) == 'earth') then
                        p = max(min(int(data(i)%tod%scans(scan)%d(j)%earth_elon(k,1) / (pi/NBIN_EARTH_ELON)), NBIN_EARTH_ELON),1)
                     end if

                     !call pix2vec_ring(data(i)%tod%nside, p, vec)
                     !elon = acos(min(max(vec(1),-1.d0),1.d0)) * 180.d0/pi
                     
                     s_sky(k) = map_sky(1, data(i)%tod%pix2ind(pix(k, 1)), j, 1)  ! zodi is only temperature (for now)
                     res      = tod(k) / data(i)%tod%scans(scan)%d(j)%gain - (s_zodi(k)+s_sky(k))
                     A(p)     = A(p) + w * amp * amp
                     b(p)     = b(p) + w * amp * res
                  end do
                  deallocate(pix, psi, flag, tod, mask)
               end do
               deallocate(s_scat, s_therm, s_zodi, s_sky)
            end do
            deallocate(procmask,map_sky, sky_signal)
         end do
         
         ! Gather information across cores
         call mpi_allreduce(MPI_IN_PLACE, A, size(A), MPI_DOUBLE_PRECISION, MPI_SUM, cpar%comm_chain, ierr)
         call mpi_allreduce(MPI_IN_PLACE, b, size(b), MPI_DOUBLE_PRECISION, MPI_SUM, cpar%comm_chain, ierr)
         
         ! Solve for best-fit map
         if (trim(map_id) == 'solar' .and. .not. associated(data(band)%tod%map_solar)) allocate(data(band)%tod%map_solar(0:12*data(band)%info%nside**2-1,1))
         if (trim(map_id) == 'moon'  .and. .not. associated(data(band)%tod%map_moon))  allocate(data(band)%tod%map_moon(0:12*data(band)%info%nside**2-1,1))
         if (trim(map_id) == 'earth' .and. .not. associated(data(band)%tod%map_earth)) allocate(data(band)%tod%map_earth(NBIN_EARTH_ELON))

         if (trim(map_id) == 'solar') then
            where (A > 0.d0)
               data(band)%tod%map_solar(:,1) = b/A
            elsewhere
               data(band)%tod%map_solar(:,1) = -1.6375d30
            end where
         end if

         if (trim(map_id) == 'moon') then
            where (A > 0.d0)
               data(band)%tod%map_moon(:,1) = b/A
            elsewhere
               data(band)%tod%map_moon(:,1) = -1.6375d30
            end where
         end if
         
         if (trim(map_id) == 'earth') then
            where (A > 0.d0)
               data(band)%tod%map_earth = b/A
            elsewhere
               data(band)%tod%map_earth = -1.6375d30
            end where
         end if

!!$         if (cpar%myid_chain == 0) then
!!$            call write_map2('static_'//trim(data(band)%label)//'.fits', real(data(band)%tod%map_solar,dp))
!!$         end if
         
         ! Clean up
         deallocate(A, b, active)
      end do

!!$      call mpi_finalize(ierr)
!!$      stop
      
    end subroutine sample_static_zodi_map


!!$    subroutine sample_static_zodi_amps(cpar, handle)
!!$     implicit none
!!$      type(comm_params), intent(inout) :: cpar
!!$      type(planck_rng), intent(inout) :: handle
!!$
!!$      integer(i4b) :: band, i, j, k, ndet, scan, nscan, npix, nmaps, p, ierr, ntod, nhorn, npix_band, ncomp, refband
!!$      real(dp)     :: res, w, vec(3), elon
!!$      character(len=128) :: reflabel
!!$      real(dp),      allocatable, dimension(:)       :: A, b
!!$      real(sp),      allocatable, dimension(:)       :: s_sky
!!$      real(sp),      allocatable, dimension(:,:)     :: s_scat, s_therm
!!$      real(sp),      allocatable, dimension(:)       :: s_zodi, s_static
!!$      real(sp),      allocatable, dimension(:,:,:,:) :: map_sky
!!$      type(map_ptr), allocatable, dimension(:,:)     :: sky_signal
!!$      real(sp),      allocatable, dimension(:)       :: tod, mask, procmask_zodi
!!$      real(dp),      allocatable, dimension(:,:)     :: m_buf
!!$      integer(i4b),  allocatable, dimension(:,:)     :: pix, psi
!!$      integer(i4b),  allocatable, dimension(:)       :: flag
!!$
!!$      if (cpar%myid == 0) then
!!$         write(*,*) '   Sampling static zodi amplitudes'
!!$      end if
!!$      
!!$      npix  = 12*cpar%zodi_solar_nside**2
!!$      nmaps = 1
!!$      ncomp = zodi_model%n_comps
!!$
!!$      ! Find reference band
!!$      refband  = 0
!!$      reflabel = 'none'
!!$      do i = 1, numband
!!$         if (trim(data(i)%label) == trim(cpar%zs_refband)) then
!!$            refband  = i
!!$            reflabel = data(refband)%instlabel
!!$            exit
!!$         end if
!!$      end do
!!$      
!!$      ! Allocate temporary map structures
!!$      allocate(A(numband), b(numband), amp_old(numband))
!!$      A = 0.d0; b = 0.d0
!!$      amp_old = zodi_model%amp_static
!!$
!!$      ! Add up contributions from all bands
!!$      do i = 1, numband
!!$         if (data(i)%tod_type == "none") cycle
!!$         if (.not. data(i)%tod%sample_zodi) cycle
!!$         if (trim(data(i)%label) == trim(cpar%zs_refband)) cycle
!!$         if (trim(data(i)%instlabel) == trim(reflabel)) cycle
!!$         
!!$         ndet      = data(i)%tod%ndet
!!$         nscan     = data(i)%tod%nscan
!!$         nhorn     = data(i)%tod%nhorn
!!$         npix_band = 12*data(i)%map%info%nside**2
!!$
!!$         ! Set amplitude to 1 during these calculations
!!$         zodi_model%amp_static(i) = 1.d0
!!$
!!$         ! Get and distribute sky signal
!!$         allocate(sky_signal(data(i)%tod%ndet,1))
!!$         do j = 1, data(i)%tod%ndet
!!$            call get_sky_signal(i, j, sky_signal(j,1)%p, mono=.true.)
!!$         end do
!!$         allocate (map_sky(nmaps, data(i)%tod%nobs, 0:data(i)%tod%ndet, 1))
!!$         call distribute_sky_maps(data(i)%tod, sky_signal, 1.e0, map_sky)
!!$
!!$         ! Initialize frequency-specific mask
!!$         allocate(m_buf(0:npix_band-1, nmaps), procmask_zodi(0:npix_band-1))
!!$         call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:, 1)
!!$         deallocate(m_buf)
!!$
!!$         ! Sum up contributions
!!$         do scan = 1, nscan
!!$            ntod = data(i)%tod%scans(scan)%ntod
!!$            allocate(s_scat(ntod,ncomp), s_therm(ntod,ncomp), s_zodi(ntod), s_sky(ntod), s_static(ntod))
!!$            
!!$            do j = 1, ndet
!!$               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
!!$               
!!$               ! Get data and pointing
!!$               allocate(pix(ntod, nhorn), psi(ntod, nhorn), flag(ntod), tod(ntod), mask(ntod))
!!$               if (data(i)%tod%compressed_tod) then
!!$                  call data(i)%tod%decompress_tod(scan, j, tod)
!!$               else
!!$                  tod = data(i)%tod%scans(scan)%d(j)%tod
!!$               end if
!!$
!!$               ! Set up mask; remove flagged samples and foreground contaminated regions
!!$               call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)
!!$               do k = 1, data(i)%tod%scans(scan)%ntod
!!$                  mask(k) = procmask_zodi(pix(k, 1))
!!$                  if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
!!$               end do
!!$               where (mask > 0.5) 
!!$                  mask = 1.
!!$               elsewhere
!!$                  mask = 0.
!!$               end where
!!$
!!$               ! Compute zodi TODs
!!$               call get_s_tot_zodi(zodi_model, data(i)%tod, j, scan, s_zodi,   pix_dynamic=pix)
!!$               call get_s_tot_zodi(zodi_model, data(i)%tod, j, scan, s_static, pix_static=data(i)%tod%scans(scan)%d(j)%pix_sol)
!!$               
!!$               ! Add residual to normal equations
!!$               w = 1.d0/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0**2
!!$               do k = 1, ntod
!!$                  if (mask(k) == 0) cycle
!!$                  s_sky(k) = map_sky(1, data(i)%tod%pix2ind(pix(k, 1)), j, 1)  ! zodi is only temperature (for now)
!!$                  A(i)     = A(i) + w * s_static(k) * s_static(k) 
!!$                  b(i)     = b(i) + w * s_static(k) * (tod(k) - (s_zodi(k)+s_sky(k)))
!!$               end do
!!$               deallocate(pix, psi, flag, tod, mask)
!!$            end do
!!$            deallocate(s_scat, s_therm, s_zodi, s_sky, s_static)
!!$         end do
!!$         deallocate(procmask_zodi,map_sky, sky_signal)
!!$      end do
!!$
!!$      ! Gather information across cores
!!$      call mpi_allreduce(MPI_IN_PLACE, A, size(A), MPI_DOUBLE_PRECISION, MPI_SUM, cpar%comm_chain, ierr)
!!$      call mpi_allreduce(MPI_IN_PLACE, b, size(b), MPI_DOUBLE_PRECISION, MPI_SUM, cpar%comm_chain, ierr)
!!$
!!$      ! Solve for best-fit amplitudes
!!$      if (cpar%myid_chain == 0) then
!!$         ! Collect data from bands with same instrument
!!$
!!$         do i = 1, numband
!!$            do j = i+1, numband
!!$               if (trim(data(j)%instlabel) == trim(data(i)%instlabel)) then
!!$                  A(i) = A(i) + A(j)
!!$                  b(i) = b(i) + b(j)
!!$               end if
!!$            end do
!!$         end do
!!$
!!$         ! Solve for new amplitudes
!!$         do i = 1, numband
!!$            if (trim(data(i)%instlabel) == trim(reflabel)) cycle
!!$            if (A(i) > 0.d0) then
!!$               if (trim(cpar%operation) == 'optimize') then
!!$                  zodi_model%amp_static(i) = b(i)/A(i)
!!$               else if (trim(cpar%operation) == 'sample') then
!!$                  zodi_model%amp_static(i) = b(i)/A(i) + rand_gauss(handle) / sqrt(A(i))
!!$               else
!!$                  write(*,*) 'Unknown operation in sample_static_zodi_amps = ', trim(cpar%operation)
!!$               end if
!!$            else
!!$               zodi_model%amp_static(i) = 0.d0
!!$            end if
!!$            zodi_model%amp_static(i) = max(zodi_model%amp_static(i), 0.d0)
!!$         end do
!!$
!!$         ! Synchronize amplitudes for bands with same instrument
!!$         do i = 1, numband
!!$            do j = i+1, numband
!!$               if (trim(data(j)%instlabel) == trim(data(i)%instlabel)) then
!!$                  zodi_model%amp_static(j) = zodi_model%amp_static(i)
!!$               end if
!!$            end do
!!$         end do
!!$
!!$         
!!$      end if
!!$
!!$      call mpi_bcast(zodi_model%amp_static, numband, MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)
!!$      
!!$      ! Output to screen
!!$      if (cpar%myid == 0) then
!!$         do i = 1, numband
!!$            if (zodi_model%amp_static(i) /= amp_old(i) .and. zodi_model%amp_static(i) /= 0.d0) then
!!$               write(*,fmt='(a,a,a,f8.3,a,f8.3)') '  Static amp: Band = ', trim(data(i)%label), ', old = ', amp_old(i), ', new = ', zodi_model%amp_static(i)
!!$            end if
!!$         end do
!!$      end if
!!$
!!$      ! Clean up
!!$      deallocate(A, b, amp_old)
!!$      
!!$      call mpi_finalize(ierr)
!!$      stop
!!$      
!!$    end subroutine sample_static_zodi_amps
    
end module
