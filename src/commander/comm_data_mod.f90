module comm_data_mod
  use comm_param_mod
  use comm_bp_mod
  use comm_noise_mod
  use comm_beam_mod
  use comm_map_mod
  use comm_tod_mod
  use comm_tod_LFI_mod
  use locate_mod
  implicit none

  type comm_data_set
     character(len=512)           :: label, unit, comp_sens
     integer(i4b)                 :: period, id_abs
     logical(lgt)                 :: sample_gain
     real(dp)                     :: gain
     character(len=128)           :: gain_comp
     integer(i4b)                 :: gain_lmin, gain_lmax

     class(comm_mapinfo), pointer :: info
     class(comm_map),     pointer :: map
     class(comm_map),     pointer :: res
     class(comm_map),     pointer :: mask
     class(comm_map),     pointer :: procmask
     class(comm_map),     pointer :: gainmask
     class(comm_tod),     pointer :: tod
     class(comm_N),       pointer :: N
     class(comm_bp),      pointer :: bp
     class(comm_B),       pointer :: B
     type(comm_B_bl_ptr),  allocatable, dimension(:) :: B_smooth
     type(comm_B_bl_ptr),  allocatable, dimension(:) :: B_postproc
     type(comm_N_rms_ptr), allocatable, dimension(:) :: N_smooth
   contains
     procedure :: RJ2data
     procedure :: apply_proc_mask
  end type comm_data_set

  integer(i4b) :: numband
  type(comm_data_set), allocatable, dimension(:) :: data
  integer(i4b),        allocatable, dimension(:) :: ind_ds

  
contains

  subroutine initialize_data_mod(cpar, handle)
    implicit none
    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b)       :: i, j, n, m, nmaps, ierr, numband_tot
    real(dp)           :: t1, t2
    character(len=512) :: dir, mapfile
    class(comm_N), pointer  :: tmp
    class(comm_mapinfo), pointer :: info_smooth, info_postproc
    real(dp), allocatable, dimension(:)   :: nu
    real(dp), allocatable, dimension(:,:) :: map, regnoise, mask_misspix

    ! Read all data sets
    numband_tot = cpar%numband
    dir = trim(cpar%datadir) // '/'
    allocate(data(numband_tot))
    n = 0
    do i = 1, numband_tot
       if (.not. cpar%ds_active(i)) cycle
       n                   = n+1
       data(n)%id_abs      = i
       data(n)%label       = cpar%ds_label(i)
       data(n)%period      = cpar%ds_period(i)
       data(n)%unit        = cpar%ds_unit(i)
       data(n)%sample_gain = cpar%ds_sample_gain(i)
       data(n)%gain_comp   = cpar%ds_gain_calib_comp(i)
       data(n)%gain_lmin   = cpar%ds_gain_lmin(i)
       data(n)%gain_lmax   = cpar%ds_gain_lmax(i)
       data(n)%comp_sens    = cpar%ds_component_sensitivity(i)
       if (cpar%myid == 0 .and. cpar%verbosity > 0) &
            & write(*,fmt='(a,i5,a,a)') '  Reading data set ', i, ' : ', trim(data(n)%label)
       call update_status(status, "data_"//trim(data(n)%label))

       ! Initialize map structures
       nmaps = 1; if (cpar%ds_polarization(i)) nmaps = 3
       data(n)%info => comm_mapinfo(cpar%comm_chain, cpar%ds_nside(i), cpar%ds_lmax(i), &
            & nmaps, cpar%ds_polarization(i))
       call get_mapfile(cpar, i, mapfile)
       !data(n)%map  => comm_map(data(n)%info, trim(dir)//trim(cpar%ds_mapfile(i)), mask_misspix=mask_misspix)
       data(n)%map  => comm_map(data(n)%info, trim(dir)//trim(mapfile), mask_misspix=mask_misspix)
       if (cpar%only_pol) data(n)%map%map(:,1) = 0.d0

!       call data(n)%map%writeFITS('data.fits')

!!$       data(n)%res => comm_map(data(n)%map)
!!$       call data(n)%res%writeFITS('res1.fits')
!!$       call data(n)%res%YtW()
!!$       call data(n)%res%Y()
!!$       call data(n)%res%writeFITS('res2.fits')
!!$       data(n)%res%map = data(n)%map%map - data(n)%res%map
!!$       call data(n)%res%writeFITS('res3.fits')
!       call mpi_finalize(ierr)
!       stop

       ! Read processing mask
       if (trim(cpar%ds_procmask) /= 'none') then
          data(n)%procmask => comm_map(data(n)%info, trim(cpar%datadir)//'/'//trim(cpar%ds_procmask), &
               & udgrade=.true.)
          !data(n)%map%map = data(n)%map%map * data(n)%procmask%map
          !call smooth_inside_procmask(data(n), cpar%ds_fwhm_proc)
       end if
       data(n)%res  => comm_map(data(n)%map)
       call update_status(status, "data_map")

       if (data(n)%sample_gain) then
          ! Read calibration mask
          if (trim(cpar%ds_maskfile_calib(i)) /= 'fullsky') then
             data(n)%gainmask => comm_map(data(n)%info, trim(cpar%datadir)//'/'//trim(cpar%ds_maskfile_calib(i)), &
                  & udgrade=.true.)
          end if
       end if

       ! Initialize beam structures
       select case (trim(cpar%ds_beamtype(i)))
       case ('b_l')
          data(n)%B => comm_B_bl(cpar, data(n)%info, n, i)
       case default
          call report_error("Unknown beam format: " // trim(cpar%ds_noise_format(i)))
       end select
       call update_status(status, "data_beam")
       
       ! Read default gain from instrument parameter file
       call read_instrument_file(trim(cpar%datadir)//'/'//trim(cpar%cs_inst_parfile), &
            & 'gain', cpar%ds_label(i), 1.d0, data(n)%gain)

       ! Initialize mask structures
       if (trim(cpar%ds_maskfile(i)) == 'fullsky') then
          data(n)%mask     => comm_map(data(n)%info)
          data(n)%mask%map =  1.d0
       else
          data(n)%mask  => comm_map(data(n)%info, trim(dir)//trim(cpar%ds_maskfile(i)))
          where (data(n)%mask%map > 0.5d0)
             data(n)%mask%map = 1.d0
          elsewhere
             data(n)%mask%map = 0.d0
          end where
       end if
       data(n)%mask%map = data(n)%mask%map * mask_misspix
       if (trim(cpar%ds_sourcemask) /= 'none') then
          call apply_source_mask(data(n)%mask, trim(cpar%datadir)//'/'//trim(cpar%ds_sourcemask), &
               & data(n)%B%r_max)
       end if
       if (cpar%only_pol) data(n)%mask%map(:,1) = 0.d0
       data(n)%map%map  = data(n)%map%map  * data(n)%mask%map ! Apply mask to map to avoid surprises
       call update_status(status, "data_mask")
       deallocate(mask_misspix)

       ! Initialize noise structures
       select case (trim(cpar%ds_noise_format(i)))
       case ('rms') 
          allocate(regnoise(0:data(n)%info%np-1,data(n)%info%nmaps))
          if (associated(data(n)%procmask)) then
             data(n)%N       => comm_N_rms(cpar, data(n)%info, n, i, 0, data(n)%mask, handle, regnoise, &
                  & data(n)%procmask)
          else
             data(n)%N       => comm_N_rms(cpar, data(n)%info, n, i, 0, data(n)%mask, handle, regnoise)
          end if
          data(n)%map%map = data(n)%map%map + regnoise  ! Add regularization noise
          deallocate(regnoise)
       case default
          call report_error("Unknown noise format: " // trim(cpar%ds_noise_format(i)))
       end select
       call update_status(status, "data_N")

       ! Initialize bandpass structures
       data(n)%bp => comm_bp(cpar, n, i)
       call update_status(status, "data_bp")

       ! Initialize smoothed data structures
       allocate(data(n)%B_smooth(cpar%num_smooth_scales))
       allocate(data(n)%B_postproc(cpar%num_smooth_scales))
       allocate(data(n)%N_smooth(cpar%num_smooth_scales))
       do j = 1, cpar%num_smooth_scales
          if (cpar%fwhm_smooth(j) > 0.d0) then
             info_smooth => comm_mapinfo(data(n)%info%comm, data(n)%info%nside, cpar%lmax_smooth(j), &
                  & data(n)%info%nmaps, data(n)%info%pol)
             data(n)%B_smooth(j)%p => &
               & comm_B_bl(cpar, info_smooth, n, i, fwhm=cpar%fwhm_smooth(j), pixwin=cpar%pixwin_smooth(j), &
               & init_realspace=.false.)
          else
             nullify(data(n)%B_smooth(j)%p)
          end if
          if (cpar%fwhm_postproc_smooth(j) > 0.d0) then
             info_postproc => comm_mapinfo(data(n)%info%comm, cpar%nside_smooth(j), cpar%lmax_smooth(j), &
                  & data(n)%info%nmaps, data(n)%info%pol)
             data(n)%B_postproc(j)%p => &
                  & comm_B_bl(cpar, info_postproc, n, i, fwhm=cpar%fwhm_postproc_smooth(j),&
                  & init_realspace=.false.)
          else
             nullify(data(n)%B_postproc(j)%p)
          end if
          if (trim(cpar%ds_noise_rms_smooth(i,j)) == 'native') then
             tmp => data(n)%N
             select type (tmp)
             class is (comm_N_rms)
                data(n)%N_smooth(j)%p => tmp
             end select
          else if (trim(cpar%ds_noise_rms_smooth(i,j)) /= 'none') then
             data(n)%N_smooth(j)%p => comm_N_rms(cpar, data(n)%info, n, i, j, data(n)%mask, handle, regnoise)
          else
             nullify(data(n)%N_smooth(j)%p)
          end if
       end do

       ! Initialize TOD structures
       if (cpar%enable_TOD_analysis) then
          if (trim(cpar%ds_tod_type(n)) == 'LFI') then
             data(n)%tod => comm_LFI_tod(cpar)             
          else
             write(*,*) 'Unrecognized TOD experiment type = ', trim(cpar%ds_tod_type(n))
             stop
          end if
       end if
    end do
    numband = n
    if (cpar%myid == 0 .and. cpar%verbosity > 0) &
         & write(*,fmt='(a,i5)') '  Number of active data sets = ', numband
    
    ! Sort bands according to nominal frequency
    allocate(ind_ds(numband), nu(numband))
    do i = 1, numband
       ind_ds(i)  = i
       nu(i)        = data(i)%bp%nu_c
    end do
    call QuickSort(ind_ds, nu)
    deallocate(nu)

    ! Dump unit conversion factors to file
    if (cpar%myid == 0) call dump_unit_conversion(cpar%outdir)

  end subroutine initialize_data_mod


  function RJ2data(self)
    implicit none

    class(comm_data_set), intent(in) :: self
    real(dp)                         :: RJ2data

    select case (trim(self%unit))
    case ('uK_cmb') 
       RJ2data = self%bp%a2t
    case ('MJy/sr') 
       RJ2data = self%bp%a2t / self%bp%f2t
    case ('y_SZ') 
       RJ2data = self%bp%a2sz
    case ('uK_RJ') 
       RJ2data = 1.d0
    case ('K km/s') ! NEW
       RJ2data = 1.d0
    end select
    
  end function RJ2data

  subroutine dump_unit_conversion(dir)
    implicit none
    character(len=*), intent(in) :: dir
    integer(i4b) :: i, q, unit

    unit = getlun()
    open(unit, file=trim(dir)//'/unit_conversions.dat', recl=1024)
    write(unit,*) '# Band   BP type   Nu_c (GHz)  a2t [K_cmb/K_RJ]' // &
         & '  t2f [MJy/K_cmb] a2sz [y_sz/K_RJ]'
    do i = 1, numband
       q = ind_ds(i)
       write(unit,fmt='(a7,a10,f10.1,3e16.5)') trim(data(q)%label), trim(data(q)%bp%type), &
            & data(q)%bp%nu_c/1.d9, data(q)%bp%a2t, 1.d0/data(q)%bp%f2t*1e6, data(q)%bp%a2sz * 1.d6
    end do
    close(unit)
  end subroutine dump_unit_conversion

  subroutine apply_source_mask(mask, sourcefile, r_max)
    implicit none
    class(comm_map),  intent(inout) :: mask
    character(len=*), intent(in)    :: sourcefile
    real(dp),         intent(in)    :: r_max

    integer(i4b)       :: i, j, l, p, unit, nlist, nmax=10000, itmp
    real(dp)           :: lon, lat, rad, vec(3)
    character(len=512) :: line
    integer(i4b), allocatable, dimension(:) :: listpix

    unit = getlun()
    open(unit,file=trim(sourcefile),recl=1024,status='old')
    do while (.true.)
       read(unit,fmt='(a)',end=23) line
       line = trim(adjustl(line))
       if (line(1:1) == '#' .or. line(1:1) == ' ') cycle
       read(line,*) lon, lat, rad
       call ang2vec(0.5d0*pi-lat*DEG2RAD, lon*DEG2RAD, vec)
       allocate(listpix(0:nmax-1))
       call query_disc(mask%info%nside, vec, r_max*rad, listpix, nlist)

       ! Sort pixel list according to increasing pixel number
       do j = 1, nlist-1
          itmp          = listpix(j)
          l             = j-1
          do while (l >= 0)
             if (listpix(l) <= itmp) exit
             listpix(l+1) = listpix(l)
             l            = l-1
          end do
          listpix(l+1)    = itmp
       end do

       ! Mask pixels belonging to current processor
       i    = 0
       j    = locate(mask%info%pix, listpix(i))
       if (j > 0) then
          do while (.true.)
             if (listpix(i) == mask%info%pix(j)) then
                mask%map(j,:) = 0.d0
                i             = i+1
                j             = j+1
             else if (listpix(i) < mask%info%pix(j)) then
                i               = i+1
             else
                j               = j+1
             end if
             if (i > nlist-1) exit
             if (j > mask%info%np) exit
          end do
       end if

       deallocate(listpix)
    end do
23  close(unit)

  end subroutine apply_source_mask

  subroutine smooth_inside_procmask(data, fwhm)
    implicit none
    class(comm_data_set), intent(in) :: data
    real(dp),             intent(in) :: fwhm

    integer(i4b) :: i, j
    real(dp)     :: w
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer :: map

    map => comm_map(data%map)
    call map%smooth(fwhm)
    do j = 1, data%map%info%nmaps
       do i = 0, data%map%info%np-1
          w = data%procmask%map(i,j)
          data%map%map(i,j) = w * data%map%map(i,j) + (1.d0-w) * map%map(i,j)
       end do
    end do
    deallocate(map)

  end subroutine smooth_inside_procmask

  subroutine apply_proc_mask(self, map)
    implicit none
    class(comm_data_set), intent(in)    :: self
    class(comm_map),      intent(inout) :: map

    if (.not. associated(self%procmask)) return
    where (self%procmask%map < 1.d0)  ! Apply processing mask
       map%map = -1.6375d30
    end where

  end subroutine apply_proc_mask

  subroutine smooth_map(info, alms_in, bl_in, map_in, bl_out, map_out)
    implicit none
    class(comm_mapinfo),                      intent(in),   target :: info
    logical(lgt),                             intent(in)           :: alms_in
    real(dp),            dimension(0:,1:),    intent(in)           :: bl_in, bl_out
    class(comm_map),                          intent(inout)        :: map_in
    class(comm_map),                          intent(out), pointer :: map_out

    integer(i4b) :: i, j, l, lmax

    map_out => comm_map(info)

    if (.not. alms_in) then
       !map_out%map = map_in%map
       call map_in%udgrade(map_out)
       call map_out%YtW
    else
       call map_in%alm_equal(map_out)
    end if


    ! Deconvolve old beam, and convolve with new beam
    lmax  = min(size(bl_in,1)-1, size(bl_out,1)-1)
    do i = 0, info%nalm-1
       l = info%lm(1,i)
       if (l > lmax) then
          map_out%alm(i,:) = 0.d0
          cycle
       end if
       do j = 1, map_out%info%nmaps
          if (bl_in(l,j) > 1.d-12) then
             map_out%alm(i,j) = map_out%alm(i,j) * bl_out(l,j) / bl_in(l,j)
          else
             map_out%alm(i,j) = 0.d0
          end if
       end do
    end do    

    ! Recompose map
    call map_out%Y

  end subroutine smooth_map

  subroutine get_mapfile(cpar, band, mapfile)
    implicit none
    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: band
    character(len=*),  intent(out)   :: mapfile

    integer(i4b)       :: i, n, unit
    character(len=512) :: filename

    filename = trim(adjustl(cpar%ds_mapfile(band)))
    n        = len(trim(adjustl(filename)))
    
    if (filename(n-3:n) == 'fits') then
       mapfile = filename
    else if (filename(n-2:n) == 'txt') then
       ! Assume filelist; pick the number given by mychain. 
       unit = getlun()
       open(unit, file=trim(cpar%datadir)//'/'//trim(filename), recl=1024)
       do i = 1, cpar%nskip_filelist + cpar%mychain
          read(unit,'(a)') mapfile
       end do
       close(unit)
    end if

  end subroutine get_mapfile

end module comm_data_mod
