module comm_zodi_samp_mod
   use comm_zodi_mod
   use comm_chisq_mod
   use powell_mod
   implicit none

   private
   public initialize_zodi_samp_mod, downsamp_invariant_structs, project_and_downsamp_sky
   public minimize_zodi_with_powell, get_chisq_priors, create_zodi_sampgroup_mask
   public apply_zodi_sampgroup_mask, sample_static_zodi_map!, sample_static_zodi_amps
   
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


   function get_chisq_priors(params, samp_group) result(chisq_prior)
     implicit none 
     real(dp), intent(in) :: params(:)
     integer(i4b), intent(in) :: samp_group
      
      real(dp) :: chisq_prior
      integer(i4b) :: i, j, ind1, ind2
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

      ! Check that sum of phase function weights is smaller than 1
      if (trim(zodi_model%phasefunc_type) == 'Hong') then
         ind1 = zodi_model%get_par_ind(param="w2")
         if (zodi_model%theta_stat(ind1,samp_group) == 0) then
            ! Find index of w2 in current sampgroup; assume both w2 and w3 are fitted
            j = 1
            do i = 1, ind1-1
               if (zodi_model%theta_stat(i,samp_group) == 0) j = j+1
            end do
            if (params(j)+params(j+1) > 1.d0) then
               !write(*,*) 'Error in phase function prior', params(j), params(j+1)
               chisq_prior = 1.d30
               return
            end if
         end if
      end if
    end function get_chisq_priors

    ! Downsamples pointing, tod and caches the zodi mask (sky mask + flags)
    subroutine downsamp_invariant_structs(cpar)!, map_sky, procmask)
      implicit none
      type(comm_params), intent(in) :: cpar

      integer(i4b) :: i, j, h, k, kp, l, g, scan, nside, npix, nmaps, ext(2), padding, ierr, ntod, ndet, nhorn, ndownsamp, box_halfwidth, thin, ntod_lowres, oper
      real(dp) :: dt_tod, theta, phi, vec0(3), vec1(3), M_ecl2gal(3,3), day, lambda_solar, lat, lon
      real(sp), allocatable :: tod(:), downsamp_vec(:, :), s_static(:)
      real(dp), allocatable, dimension(:, :) :: vec
      integer(i4b), allocatable, dimension(:,:) :: downsamp_pix
      logical(lgt), allocatable, dimension(:) :: downsamp_mask_idx
      real(sp), allocatable, dimension(:) :: downsamp_mask, downsamp_tod, downsamp_obs_time, obs_time
      !real(sp), allocatable, dimension(:) :: procmask_zodi
      logical(lgt), allocatable, dimension(:) :: mask
      type(hdf_file) :: tod_file
      character(len=4) :: scan_str
      type(comm_detscan), pointer :: d
      class(comm_scandata), allocatable :: sd

      if (cpar%myid == cpar%root) print *, "downsampling tod and pointing"
      ! For each zodi band, create downsampled residual time stream

      call ecl_to_gal_rot_mat(M_ecl2gal)

      oper = get_sd_operation_code([SD_BASE,SD_IND,SD_TOD])
      
      do i = 1, numband
         ! Only generate downsampled arrays for tod bands and bands where we have zodi
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         nside = data(i)%map%info%nside
         npix  = 12*nside**2
         nmaps = data(i)%map%info%nmaps
         ndet  = data(i)%tod%ndet
         nhorn = data(i)%tod%nhorn

         ! Distribute zodi mask
!!$         allocate(m_buf(0:npix-1, nmaps))
!!$         allocate(procmask_zodi(0:npix-1))
!!$         call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:, 1)
         
         thin   = max(nint(data(i)%tod%samprate/data(i)%tod%samprate_lowres), 1)
         dt_tod = (thin/data(i)%tod%samprate)*SECOND_TO_DAY ! dt between two samples in units of days (assumes equispaced tods)
         do scan = 1, data(i)%tod%nscan
            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
            ntod        = data(i)%tod%scans(scan)%ntod
            ntod_lowres = ntod/thin
            
            ! Initialize detector specific downsampled quantities
            do j = 1, data(i)%tod%ndet
               d => data(i)%tod%scans(scan)%d(j)               
               allocate(d%downsamp_pix_full(ntod_lowres,nhorn))
               allocate(d%downsamp_tod_full(ntod_lowres))
               allocate(d%downsamp_obs_time_full(ntod_lowres))
               allocate(mask(ntod_lowres))

               ! Get downsampled pixel, mask and obstime
               call init_scan_data(data(i)%tod, scan, oper, TODMASK_ZODI, sd)
               do k = 1, ntod_lowres
                  kp                          = (k-1)*thin + 1
                  do h = 1, nhorn
                     d%downsamp_pix_full(k,h)  = sd%pix(kp,j,h)
                  end do
                  d%downsamp_obs_time_full(k) = data(i)%tod%scans(scan)%t0(1) + (kp-1)*dt_tod
                  mask(k)                     = sd%mask(kp,j)
               end do

               ! Get TOD after subtracting static zodi
               do k = 1, ntod_lowres
                  kp = (k-1)*thin + 1
                  if (nhorn == 1) then
                     d%downsamp_tod_full(k) = sd%tod(kp,j) - data(i)%tod%scans(scan)%d(j)%gain * sd%s_objctr(kp,j,1)
                  else
                     write(*,*) 'Multi-horn zodi fitting not yet implemented'
                     stop
                  end if
               end do
               call dealloc_scan_data(sd)

               ! Remove masked samples
               ntod_lowres = count(mask)
               allocate(downsamp_pix(ntod_lowres,nhorn))
               do h = 1, nhorn
                  downsamp_pix(:,h) = pack(d%downsamp_pix_full(:,h), mask)
               end do
               d%downsamp_pix_full      = downsamp_pix
               d%downsamp_tod_full      = pack(d%downsamp_tod_full, mask)
               d%downsamp_obs_time_full = pack(d%downsamp_obs_time_full, mask)
               deallocate(mask, downsamp_pix)

               ! Allocate other downsampled quantities with same shape
               allocate(d%downsamp_point_full(ntod_lowres,nhorn,5))
               do k = 1, ntod_lowres
                  if (any(d%downsamp_pix_full(k,:)<0) .or. any(d%downsamp_pix_full(k,:)> 12*nside**2-1)) then
                     write(*,*) 'a', data(i)%tod%scanid(scan), j, nside, k, ntod_lowres, d%downsamp_pix_full(k,:), d%downsamp_tod_full(k), d%downsamp_obs_time_full(k)
                  end if
                  do h = 1, nhorn
                     call pix2ang_ring(nside, d%downsamp_pix_full(k,h), lat, lon)
                     !write(*,*) 'b'
                     lon = lon * 180.d0/pi
                     if (lon > 180.d0) lon = lon - 360.d0
                     lat = 90.d0 - lat * 180.d0/pi
                     d%downsamp_point_full(k,h,3) = lon    ! Galactic longitude
                     d%downsamp_point_full(k,h,4) = lat    ! Galactic latitude
                     call pix2vec_ring(nside, d%downsamp_pix_full(k,h), vec1)
                     vec1 = matmul(transpose(M_ecl2gal), vec1)
                     call vec2ang(vec1, lat, lon)
                     lon = lon * 180.d0/pi
                     if (lon > 180.d0) lon = lon - 360.d0
                     lat = 90.d0 - lat * 180.d0/pi
                     d%downsamp_point_full(k,h,1) = lon    ! Galactic longitude
                     d%downsamp_point_full(k,h,2) = lat    ! Galactic latitude
                     day          = data(i)%tod%scanid(scan)
                     lambda_solar = (-80.598349 + 0.98564736 * day + 1.912 * cos(2.d0*pi/365.25 * (day-94.8))) * pi/180.d0
                     d%downsamp_point_full(k,h,5) = acos(cos(lat*pi/180.d0) * cos(lon*pi/180.d0 - lambda_solar)) * 180.d0/pi ! Solar elongation
                  end do
               end do      
            end do
         end do
      end do

    end subroutine downsamp_invariant_structs
   

!!$    ! Loop over each band with zodi and precompute lookup tables for lowres zodi caching
!!$    subroutine precompute_lowres_zodi_lookups(cpar)
!!$      implicit none 
!!$      type(comm_params), intent(in) :: cpar
!!$      
!!$      integer(i4b) :: i, j, k, l, scan, n_lowres_obs, nobs_downsamp, nobs_lowres, pix_high, pix_low, ierr
!!$      integer(i4b), allocatable :: pix2ind_highres(:), ind2pix_highres(:)
!!$      real(dp), allocatable :: ind2vec_zodi_temp(:, :)
!!$      real(dp) :: rotation_matrix(3, 3)
!!$
!!$      call ecl_to_gal_rot_mat(rotation_matrix)
!!$write(*,*) 'pre1'
!!$      
!!$      do i = 1, numband
!!$         if (trim(data(i)%tod_type) == 'none') cycle
!!$         if (.not. data(i)%tod%subtract_zodi) cycle
!!$
!!$         allocate(pix2ind_highres(0:12*data(i)%tod%nside**2-1))         
!!$         pix2ind_highres = 0
!!$         
!!$         do scan = 1, data(i)%tod%nscan
!!$            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
!!$            do j = 1, data(i)%tod%ndet
!!$               do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_pix_full)
!!$                  pix2ind_highres(data(i)%tod%scans(scan)%d(j)%downsamp_pix_full(k)) = 1
!!$               end do
!!$            end do
!!$         end do
!!$         
!!$         nobs_downsamp = count(pix2ind_highres == 1)
!!$         data(i)%tod%zodi_cache_nobs_lowres = nobs_downsamp
!!$         
!!$         allocate(ind2vec_zodi_temp(3, nobs_downsamp))
!!$         allocate(ind2pix_highres(nobs_downsamp))
!!$
!!$         j = 1
!!$         do k = 0, 12*data(i)%tod%nside**2-1
!!$            if (pix2ind_highres(k) == 1) then
!!$               ind2pix_highres(j) = k
!!$               j = j+1
!!$            end if
!!$         end do
!!$
!!$         allocate(data(i)%tod%pix2ind_lowres(0:12*zodi_nside**2-1))
!!$         data(i)%tod%pix2ind_lowres = 0
!!$         ind2vec_zodi_temp = 0.
!!$
!!$         j = 1
!!$         do k = 1, nobs_downsamp
!!$            pix_high = ind2pix_highres(k)
!!$            pix_low = data(i)%tod%udgrade_pix_zodi(pix_high)
!!$            if (data(i)%tod%pix2ind_lowres(pix_low) == 0) then
!!$               data(i)%tod%pix2ind_lowres(pix_low) = j
!!$               write(*,*) zodi_nside, pix_low
!!$               call pix2vec_ring(zodi_nside, pix_low, ind2vec_zodi_temp(:, j))
!!$            end if
!!$            j =  j + 1
!!$         end do
!!$
!!$         nobs_lowres = j-1
!!$         allocate(data(i)%tod%ind2vec_ecl_lowres(3, nobs_lowres))
!!$         data(i)%tod%ind2vec_ecl_lowres = ind2vec_zodi_temp(:, 1:nobs_lowres)
!!$         
!!$         do k = 1, nobs_lowres
!!$            data(i)%tod%ind2vec_ecl_lowres(:, k) = matmul(data(i)%tod%ind2vec_ecl_lowres(:, k), rotation_matrix)
!!$         end do   
!!$         deallocate(ind2vec_zodi_temp, pix2ind_highres, ind2pix_highres)
!!$      end do
!!$
!!$write(*,*) 'pre2'
!!$    end subroutine precompute_lowres_zodi_lookups

   subroutine project_and_downsamp_sky(cpar)
     implicit none
     type(comm_params), intent(in) :: cpar
     
     integer(i4b) :: i, j, k, h, scan, nhorn, ndet, ntod_lowres, pix

     if (cpar%myid == cpar%root) print *, "projecting downsampled sky"

     ! For each zodi band, create downsampled residual time stream
     do i = 1, numband
        if (trim(data(i)%tod_type) == 'none') cycle
        if (.not. data(i)%tod%subtract_zodi) cycle

        nhorn = data(i)%tod%nhorn
        ndet  = data(i)%tod%ndet

        ! Project sky signal into already downsampled data structures
        do scan = 1, data(i)%tod%nscan
           if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
           do j = 1, ndet
              ntod_lowres = size(data(i)%tod%scans(scan)%d(j)%downsamp_pix_full)
              if (.not. allocated(data(i)%tod%scans(scan)%d(j)%downsamp_sky_full)) then
                 allocate(data(i)%tod%scans(scan)%d(j)%downsamp_sky_full(ntod_lowres,nhorn))
              end if
              do h = 1, nhorn
                 do k = 1, ntod_lowres
                    pix = data(i)%tod%scans(scan)%d(j)%downsamp_pix_full(k,h)
                    data(i)%tod%scans(scan)%d(j)%downsamp_sky_full(k,h) = &
                         & data(i)%tod%pixcache%map_sky(1, &
                         & data(i)%tod%pixcache%pix2ind(pix), j, 1)  ! T only
                 end do
              end do
           end do
        end do
     end do
   end subroutine project_and_downsamp_sky

!!$   subroutine compute_downsamp_zodi(cpar, model)
!!$      type(comm_params), intent(in) :: cpar
!!$      type(ZodiModel), intent(inout) :: model
!!$      integer(i4b) :: i, j, scan, ierr, ndet
!!$
!!$      if (cpar%myid == cpar%root) print *, "computing downsampled zodi"
!!$
!!$      do i = 1, numband
!!$         if (trim(data(i)%tod_type) == 'none') cycle
!!$         if (.not. data(i)%tod%subtract_zodi) cycle
!!$
!!$         do scan = 1, data(i)%tod%nscan
!!$            if (.not. any(data(i)%tod%scans(scan)%d%accept)) cycle
!!$
!!$            ndet = data(i)%tod%ndet
!!$            do j = 1, ndet
!!$               call get_zodi_emission(&
!!$                   & tod=data(i)%tod, &
!!$                   & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
!!$                   & scan=scan, &
!!$                   & det=j, &
!!$                   & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
!!$                   & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
!!$                   & model=model, &
!!$                   & use_lowres_pointing=.true. &
!!$               &)
!!$               call get_s_zodi(i, &
!!$                   & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
!!$                   & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
!!$                   & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
!!$                   &)
!!$            end do
!!$         end do
!!$      end do
!!$   end subroutine

   subroutine create_zodi_sampgroup_mask(cpar, handle)
     type(comm_params), intent(in) :: cpar
     type(planck_rng), intent(inout) :: handle
     
      integer(i4b) :: i, j, k, l, h, scan, ierr, non_glitch_size, thinstep, offset, ntod, ngood
      real(dp) :: rms, frac, thin_frac
      real(sp), allocatable :: res(:)
      real(sp), allocatable :: downsamp_scat_comp(:, :), downsamp_therm_comp(:, :)

      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         do scan = 1, data(i)%tod%nscan
            do j = 1, data(i)%tod%ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

               ! Initialize masks
               ntod = size(data(i)%tod%scans(scan)%d(j)%downsamp_tod_full)
               allocate(data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(ntod,cpar%zs_num_samp_groups))
               data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask = .true.
                                 
               !allocate(basemask(ntod))
               
               ! Search for strong outliers
!!$               res = data(i)%tod%scans(scan)%d(j)%downsamp_tod - data(i)%tod%scans(scan)%d(j)%gain * (data(i)%tod%scans(scan)%d(j)%downsamp_sky + data(i)%tod%scans(scan)%d(j)%downsamp_zodi)
!!$               rms = sqrt(mean(real(res**2, dp)))
!!$               basemask = abs(res) > 10. * real(rms, sp)
!!$               frac = count(basemask)/real(ntod,dp)
!!$               if (frac > 0.01d0) write(*,*) 'Warning: Removing high fraction of glitches = ', frac

               ! Apply TOD thinning; may be different for each sampling group
               do l = 1, cpar%zs_num_samp_groups

                  ! Remove samples at high latitudes, if requested by user
                  if (cpar%zs_samp_group_max_b_ecl(l) < 90.d0) then
                     do k = 1, ntod
                        if (all(abs(data(i)%tod%scans(scan)%d(j)%downsamp_point_full(k,:,3)) > cpar%zs_samp_group_max_b_ecl(l))) then
                           data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(k,l) = .false.
                        end if
                     end do
                  end if

                  ! Apply thinning factor, adjusted for already masked samples
                  ngood    = count(data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(:,l))
                  thinstep = max(nint(real(ngood,sp)/real(ntod,sp) / cpar%zs_tod_thin_factor),1)
                  k        = 1
                  ngood    = 0
                  do while (k <= ntod)
                     if (data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(k,l)) then
                        ngood = ngood + 1
                        if (mod(ngood,thinstep) /= 0) data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(k,l) = .false.
                     end if
                     k = k+1
                  end do
               end do
            end do
         end do
      end do
    end subroutine create_zodi_sampgroup_mask

   subroutine apply_zodi_sampgroup_mask(cpar, samp_group)
     implicit none
     type(comm_params), intent(in) :: cpar
     integer(i4b),      intent(in) :: samp_group
     
      integer(i4b) :: i, j, k, h, scan, ngood, nhorn
      
      do i = 1, numband
         if (trim(data(i)%tod_type) == 'none') cycle
         if (.not. data(i)%tod%subtract_zodi) cycle

         nhorn = data(i)%tod%nhorn
         do scan = 1, data(i)%tod%nscan
            do j = 1, data(i)%tod%ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
               
               ! Remove flagged samples
               ngood = count(data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(:,samp_group))
               if (allocated(data(i)%tod%scans(scan)%d(j)%downsamp_pix)) then
                  deallocate(data(i)%tod%scans(scan)%d(j)%downsamp_pix)
                  deallocate(data(i)%tod%scans(scan)%d(j)%downsamp_sky)
               end if
               if (.not. allocated(data(i)%tod%scans(scan)%d(j)%downsamp_pix)) then
                  allocate(data(i)%tod%scans(scan)%d(j)%downsamp_pix(ngood,nhorn))
                  allocate(data(i)%tod%scans(scan)%d(j)%downsamp_sky(ngood,nhorn))
               end if
               do h = 1, nhorn
                  data(i)%tod%scans(scan)%d(j)%downsamp_pix(:,h) = pack(data(i)%tod%scans(scan)%d(j)%downsamp_pix_full(:,h), data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(:,samp_group))
                  data(i)%tod%scans(scan)%d(j)%downsamp_sky(:,h) = pack(data(i)%tod%scans(scan)%d(j)%downsamp_sky_full(:,h),  data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(:,samp_group))
               end do
               data(i)%tod%scans(scan)%d(j)%downsamp_tod  = pack(data(i)%tod%scans(scan)%d(j)%downsamp_tod_full,  data(i)%tod%scans(scan)%d(j)%zodi_sampgroup_mask(:,samp_group))
               
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
         if (ierr /= 0) write(*,*) 'powell failed, ierr =', ierr
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
      real(dp) :: chisq, chisq_tot, t1, t2, t3, t4, mono
      integer(i4b) :: i, j, k, h, scan, ntod, ndet, nscan, flag, ierr, ndof, ndof_tot, nhorn
      logical(lgt) :: accept
      character(len=4) :: scan_str
      type(hdf_file) :: tod_file
      real(sp), allocatable, dimension(:,:) :: s_zodi

      call wall_time(t1)
      
      allocate(theta(npar))
      if (cpar%myid_chain == 0) then
         flag = 1
         call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
         theta = p*scale
      end if
      call mpi_bcast(theta, size(theta), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)
            
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
         nhorn = data(i)%tod%nhorn

         ! Get monopole
         mono = band_monopole(i)
         
         ! Evaluate zodi model with newly proposed values for each band and calculate chisq
         do scan = 1, nscan
            ! Skip scan if no accepted data
            do j = 1, ndet
               if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

               allocate(s_zodi(size(data(i)%tod%scans(scan)%d(j)%downsamp_tod),nhorn))
               
               call wall_time(t3)
               do h = 1, nhorn
                  call get_zodi_emission(data(i)%tod, data(i)%tod%scans(scan)%d(j)%downsamp_pix(:,h), &
                       & scan, j, zodi_model, s_zodi(:,h), use_lowres_pointing=.true.)
               end do
               call wall_time(t4)

               if (nhorn > 1) then
                  ! Coadd horns; store result in first column of s_zodi
                  write(*,*) 'Multi-horn zodi fitting not yet implemented'
                  stop
               end if
               !if (data(1)%tod%myid == 10) write(*,*) ' CPU1 = ', t4-t3

               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'tod  ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_tod))
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'sky  ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_sky))
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'zodi ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_zodi))
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'mono ', mono
               !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'noise', data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
               chisq = sum( &
                    & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod - data(i)%tod%scans(scan)%d(j)%gain * &
                    & (data(i)%tod%scans(scan)%d(j)%downsamp_sky(:,1) + s_zodi(:,1) + mono))/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
               chisq_tot     = chisq_tot     + chisq
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
                     write(58,*) data(i)%tod%scans(scan)%d(j)%downsamp_point(k,1,:), data(i)%tod%scans(scan)%d(j)%downsamp_tod(k) &
                          &   - data(i)%tod%scans(scan)%d(j)%gain * &
                          &     (data(i)%tod%scans(scan)%d(j)%downsamp_sky(k,1) &
                          &   + s_zodi(k,1) &
                          &   + mono)
                  end do
                  close(58)
               end if
               call wall_time(t3)
               !if (data(1)%tod%myid == 10) write(*,*) ' CPU4 = ', t3-t4

               deallocate(s_zodi)
            end do
         end do
      end do

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

   end function

 end subroutine minimize_zodi_with_powell

 subroutine parse_samp_group_strings(samp_group_str, param_labels, param_indices)
   implicit none
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

      integer(i4b) :: band, i, j, k, ndet, scan, nscan, npix, nmaps, p, ierr, ntod, nhorn, npix_band, ncomp, nactive, oper
      real(dp)     :: res, w, vec(3), elon, amp
      character(len=512) :: model
      class(comm_scandata), allocatable :: sd
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

      oper = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
           & SD_SKY,SD_BP,SD_SL,SD_ORB,SD_INST,SD_ZODI])
      
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

!!$            ! Get and distribute sky signal
!!$            allocate(sky_signal(data(i)%tod%ndet,1))
!!$            do j = 1, data(i)%tod%ndet
!!$               call get_sky_signal(i, j, sky_signal(j,1)%p, mono=.true.)
!!$            end do
!!$            allocate (map_sky(nmaps, data(i)%tod%pixcache%nobs, 0:data(i)%tod%ndet, 1))
!!$            call distribute_sky_maps(data(i)%tod, sky_signal, 1.e0, map_sky)
!!$            
!!$            ! Initialize frequency-specific mask
!!$            allocate(m_buf(0:npix_band-1, nmaps), procmask(0:npix_band-1))
!!$            !call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask_zodi = m_buf(:, 1)
!!$            call data(i)%tod%procmask_zodi%bcast_fullsky_map(m_buf); procmask = m_buf(:, 1)
!!$            deallocate(m_buf)
         
            do scan = 1, nscan
               ntod = data(i)%tod%scans(scan)%ntod
               !allocate(s_scat(ntod,ncomp), s_therm(ntod,ncomp), s_zodi(ntod), s_sky(ntod))
               call init_scan_data(data(i)%tod, scan, oper, TODMASK_ZODI, sd)
               
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                  
                  ! Get data and pointing
!!$                  allocate(pix(ntod, nhorn), psi(ntod, nhorn), flag(ntod), tod(ntod), mask(ntod))
!!$                  if (data(i)%tod%compressed_tod) then
!!$                     call data(i)%tod%decompress_tod(scan, j, tod)
!!$                  else
!!$                     tod = data(i)%tod%scans(scan)%d(j)%tod
!!$                  end if
!!$                  
!!$                  ! Set up mask; remove flagged samples and foreground contaminated regions
!!$                  call data(i)%tod%decompress_pointing_and_flags(scan, j, pix, psi, flag)
!!$                  do k = 1, data(i)%tod%scans(scan)%ntod
!!$                     mask(k) = procmask(pix(k, 1))
!!$                     if (iand(flag(k), data(i)%tod%flag0) .ne. 0) mask(k) = 0.
!!$                  end do
!!$                  where (mask > 0.5) 
!!$                     mask = 1.
!!$                  elsewhere
!!$                     mask = 0.
!!$                  end where
!!$
!!$                  ! Compute non-stationary zodi TOD; exclude current static component
!!$                  call get_s_tot_zodi(zodi_model, data(i)%tod, j, scan, s_zodi, pix_dynamic=pix, exclude_static=map_id)
                  
                  ! Add residual to mapmaking equation in solar centric coordinates
                  w  = 1.d0/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0**2
                  amp = 1.d0 !zodi_model%amp_static(i)
                  do k = 1, ntod
                     if (sd%mask(k,j) == 0) cycle
                     if (trim(map_id) == 'solar') then
                        p = data(i)%tod%scans(scan)%d(j)%pix_sol(k,1)
                     else if (trim(map_id) == 'moon') then
                        p = data(i)%tod%scans(scan)%d(j)%pix_moon(k,1)
                     else if (trim(map_id) == 'earth') then
                        p = max(min(int(data(i)%tod%scans(scan)%d(j)%earth_elon(k,1) / (pi/NBIN_EARTH_ELON)), NBIN_EARTH_ELON),1)
                     end if

                     !call pix2vec_ring(data(i)%tod%nside, p, vec)
                     !elon = acos(min(max(vec(1),-1.d0),1.d0)) * 180.d0/pi
                     res      = sd%tod(k,j) / data(i)%tod%scans(scan)%d(j)%gain - sd%s_tot(k,j,0,1)
                     A(p)     = A(p) + w * amp * amp
                     b(p)     = b(p) + w * amp * res
                  end do
               end do
               call dealloc_scan_data(sd)
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

  end module comm_zodi_samp_mod
