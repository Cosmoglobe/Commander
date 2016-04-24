module comm_data_mod
  use comm_param_mod
  use comm_bp_mod
  use comm_noise_mod
  use comm_beam_mod
  use comm_map_mod
  implicit none

  type comm_data_set
     logical(lgt)                 :: active
     character(len=512)           :: label, unit
     integer(i4b)                 :: period
     real(dp)                     :: gain

     class(comm_mapinfo), pointer :: info
     class(comm_map),     pointer :: map
     class(comm_map),     pointer :: res
     class(comm_map),     pointer :: mask
     class(comm_N),       pointer :: N
     class(comm_bp),      pointer :: bp
     class(comm_B),       pointer :: B
   contains
     procedure :: RJ2data
  end type comm_data_set

  integer(i4b) :: numband
  type(comm_data_set), allocatable, dimension(:) :: data
  integer(i4b),        allocatable, dimension(:) :: ind_ds

  
contains

  subroutine initialize_data_mod(cpar, handle)
    implicit none
    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b)       :: i, j, m, nmaps
    character(len=512) :: dir
    real(dp), allocatable, dimension(:)   :: nu
    real(dp), allocatable, dimension(:,:) :: map, regnoise

    ! Read all data sets
    numband = cpar%numband
    dir = trim(cpar%datadir) // '/'
    allocate(data(numband))
    do i = 1, numband
       data(i)%active       = cpar%ds_active(i)
       data(i)%label        = cpar%ds_label(i)
       data(i)%period       = cpar%ds_period(i)
       data(i)%unit         = cpar%ds_unit(i)
       if (cpar%myid == 0 .and. cpar%verbosity > 0) &
            & write(*,fmt='(a,i5,a,a)') '  Reading data set ', i, ' : ', trim(data(i)%label)
       call update_status(status, "data_"//trim(data(i)%label))

       ! Initialize map structures
       nmaps = 1; if (cpar%ds_polarization(i)) nmaps = 3
       data(i)%info => comm_mapinfo(cpar%comm_chain, cpar%ds_nside(i), cpar%ds_lmax(i), &
            & nmaps, cpar%ds_polarization(i))
       data(i)%map  => comm_map(data(i)%info, trim(dir)//trim(cpar%ds_mapfile(i)))
       data(i)%res  => comm_map(data(i)%map)
       call update_status(status, "data_map")

       ! Read default gain from instrument parameter file
       call read_instrument_file(trim(cpar%datadir)//'/'//trim(cpar%cs_inst_parfile), &
            & 'gain', cpar%ds_label(i), 1.d0, data(i)%gain)

       ! Initialize mask structures
       if (trim(cpar%ds_maskfile(i)) == 'fullsky') then
          data(i)%mask     => comm_map(data(i)%info)
          data(i)%mask%map =  1.d0
       else
          data(i)%mask  => comm_map(data(i)%info, trim(dir)//trim(cpar%ds_maskfile(i)))
          where (data(i)%mask%map > 0.5d0)
             data(i)%mask%map = 1.d0
          elsewhere
             data(i)%mask%map = 0.d0
          end where
       end if
       call update_status(status, "data_mask")

       ! Initialize noise structures
       select case (trim(cpar%ds_noise_format(i)))
       case ('rms') 
          allocate(regnoise(0:data(i)%info%np-1,data(i)%info%nmaps))
          data(i)%N       => comm_N_rms(cpar, data(i)%info, i, data(i)%mask, handle, regnoise)
          data(i)%map%map = data(i)%map%map + regnoise  ! Add regularization noise
          deallocate(regnoise)
       case default
          call report_error("Unknown noise format: " // trim(cpar%ds_noise_format(i)))
       end select
       call update_status(status, "data_N")

       ! Initialize bandpass structures
       data(i)%bp => comm_bp(cpar, i)
       call update_status(status, "data_bp")
       
       ! Initialize beam structures
       select case (trim(cpar%ds_beamtype(i)))
       case ('b_l')
          data(i)%B => comm_B_bl(cpar, data(i)%info, i)
       case default
          call report_error("Unknown beam format: " // trim(cpar%ds_noise_format(i)))
       end select
       call update_status(status, "data_beam")

    end do
    
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

end module comm_data_mod
