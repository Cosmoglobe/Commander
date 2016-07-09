module comm_data_mod
  use comm_param_mod
  use comm_bp_mod
  use comm_noise_mod
  use comm_beam_mod
  use comm_map_mod
  use locate_mod
  implicit none

  type comm_data_set
     character(len=512)           :: label, unit
     integer(i4b)                 :: period
     real(dp)                     :: gain

     class(comm_mapinfo), pointer :: info
     class(comm_map),     pointer :: map
     class(comm_map),     pointer :: res
     class(comm_map),     pointer :: mask
     class(comm_map),     pointer :: procmask
     class(comm_N),       pointer :: N
     class(comm_bp),      pointer :: bp
     class(comm_B),       pointer :: B
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
    character(len=512) :: dir
    real(dp), allocatable, dimension(:)   :: nu
    real(dp), allocatable, dimension(:,:) :: map, regnoise, mask_misspix

    ! Read all data sets
    numband_tot = cpar%numband
    dir = trim(cpar%datadir) // '/'
    allocate(data(numband_tot))
    n = 0
    do i = 1, numband_tot
       if (.not. cpar%ds_active(i)) cycle
       n              = n+1
       data(n)%label  = cpar%ds_label(i)
       data(n)%period = cpar%ds_period(i)
       data(n)%unit   = cpar%ds_unit(i)
       if (cpar%myid == 0 .and. cpar%verbosity > 0) &
            & write(*,fmt='(a,i5,a,a)') '  Reading data set ', i, ' : ', trim(data(n)%label)
       call update_status(status, "data_"//trim(data(n)%label))

       ! Initialize map structures
       nmaps = 1; if (cpar%ds_polarization(i)) nmaps = 3
       data(n)%info => comm_mapinfo(cpar%comm_chain, cpar%ds_nside(i), cpar%ds_lmax(i), &
            & nmaps, cpar%ds_polarization(i))
       data(n)%map  => comm_map(data(n)%info, trim(dir)//trim(cpar%ds_mapfile(i)), mask_misspix=mask_misspix)

       ! Read processing mask
       if (trim(cpar%ds_procmask) /= 'none') then
          data(n)%procmask => comm_map(data(n)%info, trim(cpar%datadir)//'/'//trim(cpar%ds_procmask), &
               & udgrade=.true.)
          !data(n)%map%map = data(n)%map%map * data(n)%procmask%map
          call smooth_inside_procmask(data(n), cpar%ds_fwhm_proc)
       end if
       data(n)%res  => comm_map(data(n)%map)
       call update_status(status, "data_map")

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
       call update_status(status, "data_mask")
       deallocate(mask_misspix)

       ! Initialize noise structures
       select case (trim(cpar%ds_noise_format(i)))
       case ('rms') 
          allocate(regnoise(0:data(n)%info%np-1,data(n)%info%nmaps))
          data(n)%N       => comm_N_rms(cpar, data(n)%info, n, i, data(n)%mask, handle, regnoise)
          data(n)%map%map = data(n)%map%map + regnoise  ! Add regularization noise
          deallocate(regnoise)
       case default
          call report_error("Unknown noise format: " // trim(cpar%ds_noise_format(i)))
       end select
       call update_status(status, "data_N")

       ! Initialize bandpass structures
       data(n)%bp => comm_bp(cpar, n, i)
       call update_status(status, "data_bp")
       
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

end module comm_data_mod
