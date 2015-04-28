module comm_data_mod
  use comm_param_mod
  use comm_bp_mod
  use comm_noise_mod
  use comm_beam_mod
  use comm_map_mod
  use comm_F_mod
  implicit none

  type comm_data_set
     logical(lgt)                 :: active
     character(len=512)           :: label, unit
     integer(i4b)                 :: period

     class(comm_mapinfo), pointer :: info
     class(comm_map),     pointer :: map
     class(comm_map),     pointer :: mask
     class(comm_N),       pointer :: N
     class(comm_bp),      pointer :: bp
     class(comm_B),       pointer :: B
     class(comm_F),       pointer :: F     ! Linked mixmat list
   contains
     procedure :: RJ2data
  end type comm_data_set

  integer(i4b) :: numband
  type(comm_data_set), allocatable, dimension(:) :: data
  integer(i4b),        allocatable, dimension(:) :: ind_ds

  
contains

  subroutine initialize_data_mod(cpar)
    implicit none

    type(comm_params), intent(in) :: cpar

    integer(i4b)       :: i, j, m, nmaps
    character(len=512) :: dir
    real(dp), allocatable, dimension(:)   :: nu
    real(dp), allocatable, dimension(:,:) :: map

    ! Read all data sets
    numband = cpar%numband
    dir = trim(cpar%datadir) // '/'
    allocate(data(numband))
    do i = 1, numband
       data(i)%active       = cpar%ds_active(i)
       data(i)%label        = cpar%ds_label(i)
       data(i)%period       = cpar%ds_period(i)
       data(i)%unit         = cpar%ds_unit(i)

       ! Initialize map structures
       nmaps = 1; if (cpar%ds_polarization(i)) nmaps = 3
       data(i)%info => comm_mapinfo(cpar%comm_chain, cpar%ds_nside(i), cpar%ds_lmax(i), &
            & nmaps, cpar%ds_polarization(i))
       data(i)%map  => comm_map(data(i)%info, trim(dir)//trim(cpar%ds_mapfile(i)))
       call update_status(status, "data_map")

       ! Initialize mask structures
       if (trim(cpar%ds_maskfile(i)) == 'fullsky') then
          data(i)%mask     => comm_map(data(i)%info)
          data(i)%mask%map =  1.d0
       else
          data(i)%mask  => comm_map(data(i)%info, trim(dir)//trim(cpar%ds_maskfile(i)))
       end if
       call update_status(status, "data_mask")

       ! Initialize noise structures
       select case (trim(cpar%ds_noise_format(i)))
       case ('rms') 
          data(i)%N => comm_N_rms(cpar, data(i)%info, i, data(i)%mask)
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
    case ('K km/s') 
       RJ2data = self%bp%a2t / self%bp%co2t
    case ('y_SZ') 
       RJ2data = self%bp%a2sz
    case ('uK_RJ') 
       RJ2data = 1.d0
    end select
    
  end function RJ2data

end module comm_data_mod
