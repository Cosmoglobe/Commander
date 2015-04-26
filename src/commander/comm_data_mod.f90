module comm_data_mod
  use comm_param_mod
  use comm_bp_mod
  use comm_noise_mod
  use comm_map_mod
  implicit none

  type comm_data_set
     logical(lgt)                 :: active, pol
     character(len=512)           :: label, unit, beamtype
     integer(i4b)                 :: period

     class(comm_mapinfo), pointer :: info
     class(comm_map),     pointer :: map
     class(comm_map),     pointer :: mask
     class(comm_N),       pointer :: N
     class(comm_bp),      pointer :: bp

     real(dp),     allocatable, dimension(:,:)   :: f
     real(dp),     allocatable, dimension(:,:)   :: b_l
   contains
     procedure :: RJ2data
  end type comm_data_set

  type(comm_data_set), allocatable, dimension(:) :: data
  integer(i4b),        allocatable, dimension(:) :: i2f

  
contains

  subroutine initialize_data_mod(cpar)
    implicit none

    type(comm_params), intent(in) :: cpar

    integer(i4b)       :: i, j, n, m
    character(len=512) :: dir
    real(dp), allocatable, dimension(:,:) :: map

    ! Read all data sets
    n = cpar%numband
    dir = trim(cpar%datadir) // '/'
    allocate(data(n))
    do i = 1, cpar%numband
       data(i)%active       = cpar%ds_active(i)
       data(i)%label        = cpar%ds_label(i)
       data(i)%period       = cpar%ds_period(i)
       data(i)%unit         = cpar%ds_unit(i)

       ! Initialize map structures
       data(i)%info => comm_mapinfo(cpar%comm_chain, cpar%ds_nside(i), cpar%ds_lmax(i), &
            & cpar%ds_polarization(i))
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
          call report_error("Unknown file format: " // trim(cpar%ds_noise_format(i)))
       end select
       call update_status(status, "data_N")

       ! Initialize bandpass structures
       data(i)%bp => comm_bp(cpar, i)
       call update_status(status, "data_bp")
       
       ! Initialize beam structures
!       data(i)%beamtype = cpar%ds_beamtype(i)
!       call read_beam(data(i)%lmax, data(i)%nmaps, &
!            & data(i)%b_l, beamfile=trim(dir)//trim(cpar%ds_blfile(i)), &
!            & pixwin=trim(dir)//trim(cpar%ds_pixwin(i)))
!       call update_status(status, "data_beam")

    end do
    
!!$    ! Sort bands according to nominal frequency
!!$    allocate(i2f(numband), freq(numband))
!!$    do i = 1, numband
!!$       i2f(i)  = i
!!$       freq(i) = bp(i)%nu_c
!!$    end do
!!$    call QuickSort(i2f,freq)
!!$    deallocate(freq)

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
