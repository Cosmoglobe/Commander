module comm_data_mod
  use comm_param_mod
  use comm_bp_mod
  use comm_noise_mod
  implicit none

  type comm_data_set
     class(comm_N),  pointer :: N
     class(comm_bp), pointer :: bp

     logical(lgt)           :: active, pol
     character(len=512)     :: label, unit, beamtype
     integer(i4b)           :: nside, nmaps, npix, np, lmax, ntemp, period
     real(dp),     allocatable, dimension(:,:)   :: map
     real(dp),     allocatable, dimension(:,:)   :: mask
     real(dp),     allocatable, dimension(:,:)   :: mask_calib
     real(dp),     allocatable, dimension(:,:)   :: f
     real(dp),     allocatable, dimension(:,:)   :: b_l
     integer(i4b), allocatable, dimension(:)     :: rings
     integer(i4b), allocatable, dimension(:)     :: pix
     logical(lgt), allocatable, dimension(:)     :: samptemp
     real(dp),     allocatable, dimension(:,:,:) :: T
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
    character(len=512) :: dir, noise_format
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
       data(i)%nside        = cpar%ds_nside(i)
       data(i)%lmax         = cpar%ds_lmax(i)
       data(i)%pol          = cpar%ds_polarization(i)
       data(i)%nmaps        = 1; if (data(i)%pol) data(i)%nmaps = 3
       noise_format         = cpar%ds_noise_format(i)

       ! Initialize map structures
       call allocate_map(cpar%comm_chain, data(i)%nside, data(i)%nmaps, 0, &
            & data(i)%map, data(i)%np, data(i)%rings, data(i)%pix, trim(dir)//cpar%ds_mapfile(i))
       if (trim(cpar%ds_maskfile(i)) == 'fullsky') then
          call allocate_map(cpar%comm_chain, data(i)%nside, data(i)%nmaps, 0, &
               & data(i)%mask)
          data(i)%mask = 1.d0
       else
          call allocate_map(cpar%comm_chain, data(i)%nside, data(i)%nmaps, 0, &
               & data(i)%mask, filename=trim(dir)//cpar%ds_maskfile(i))
       end if
       if (trim(cpar%ds_maskfile_calib(i)) == 'fullsky') then
          call allocate_map(cpar%comm_chain, data(i)%nside, data(i)%nmaps, 0, &
               & data(i)%mask_calib)
          data(i)%mask_calib = 1.d0
       else
          call allocate_map(cpar%comm_chain, data(i)%nside, data(i)%nmaps, 0, &
               & data(i)%mask_calib, filename=trim(dir)//cpar%ds_maskfile_calib(i))
       end if

       ! Initialize noise object
       if (trim(noise_format) == 'rms') then
          data(i)%N => comm_N_rms(cpar, i, data(i)%mask)
       else
          call report_error(cpar%myid, "Unknown file format: " // trim(noise_format))
       end if

       ! Initialize bandpass
       data(i)%bp => comm_bp(cpar, i)
       
       ! Initialize beam structures
       data(i)%beamtype = cpar%ds_beamtype(i)
       call read_beam(data(i)%lmax, data(i)%nmaps, &
            & data(i)%b_l, beamfile=trim(dir)//trim(cpar%ds_blfile(i)), &
            & pixwin=trim(dir)//trim(cpar%ds_pixwin(i)))

       ! Initialize template structures
       data(i)%ntemp = 0
       if (cpar%ds_samp_monopole(i)) data(i)%ntemp = data(i)%ntemp+1
       if (cpar%ds_samp_dipole(i))   data(i)%ntemp = data(i)%ntemp+3
                                     data(i)%ntemp = data(i)%ntemp+cpar%ds_numtemp(i)
       allocate(data(i)%T(data(i)%np,data(i)%nmaps,data(i)%ntemp), data(i)%samptemp(data(i)%ntemp))
       data(i)%T = 0.d0
       m         = 1
       if (cpar%ds_samp_monopole(i)) then
          call initialize_mono_dipole(data(i)%nside, data(i)%pix, monopole=data(i)%T(:,1,m))
          data(i)%samptemp(m) = cpar%ds_samp_monopole(i)
          m = m+1
       end if
       if (cpar%ds_samp_dipole(i)) then
          call initialize_mono_dipole(data(i)%nside, data(i)%pix, dipole=data(i)%T(:,1,m:m+2))
          data(i)%samptemp(m:m+2) = cpar%ds_samp_dipole(i)
          m = m+3
       end if
       do j = 1, cpar%ds_numtemp(i)
          call allocate_map(cpar%comm_chain, data(i)%nside, data(i)%nmaps, 0, &
               & map, filename=trim(dir)//cpar%ds_tempname(i,j))
          data(i)%T(:,:,m)    = map
          data(i)%samptemp(m) = cpar%ds_samptemp(i,j)
          m                   = m+1
          deallocate(map)
       end do

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
