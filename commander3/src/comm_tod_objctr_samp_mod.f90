module comm_tod_objctr_samp_mod
  use comm_data_mod
  implicit none
  
  private
  public sample_objctr_map

contains

  subroutine sample_objctr_map(cpar, handle, map_id)
    implicit none
    type(comm_params), intent(inout) :: cpar
    type(planck_rng),  intent(inout) :: handle
    character(len=*),  intent(in)    :: map_id
    
    integer(i4b) :: band, i, j, k, ndet, scan, nscan, npix, nmaps, p, ierr, ntod, nhorn, npix_band, nactive, oper
    real(dp)     :: res, w, vec(3), elon, amp
    character(len=512) :: model
    type(comm_scandata) :: sd
    real(dp),      allocatable, dimension(:)       :: A, b
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
    
    nmaps = 1
    
    oper = get_sd_operation_code([SD_TOT,SD_BASE,SD_IND,SD_MASK,SD_TOD,&
         & SD_SKY,SD_INST,SD_ZODI,SD_SPUR])
    
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
          
          do scan = 1, nscan
             ntod = data(i)%tod%scans(scan)%ntod
             call init_scan_data(data(i)%tod, scan, oper, TODMASK_ZODI, sd)
             
             do j = 1, ndet
                if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                
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
                   res      = (sd%tod(k,j)-sd%s_spur(k,j)) / data(i)%tod%scans(scan)%d(j)%gain - sd%s_tot(k,j,0,1)
                   A(p)     = A(p) + w * amp * amp
                   b(p)     = b(p) + w * amp * res
                end do
             end do
             call dealloc_scan_data(sd)
          end do
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
       
       ! Clean up
       deallocate(A, b, active)
    end do
    
  end subroutine sample_objctr_map
  
end module comm_tod_objctr_samp_mod
