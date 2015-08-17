module comm_genvec_mod
  use healpix_types
  use comm_utils
  implicit none

  integer(i4b), private :: nmaps, numcomp
  integer(i4b), private :: num_signal, npix, numband, lmax
  integer(i4b), private :: num_fg_temp
  logical(lgt), private :: sample_fg_pix, sample_fg_temp, reduce_pix

  type genvec
     real(dp), allocatable, dimension(:,:)     :: cmb_amp
     real(dp), allocatable, dimension(:,:,:)   :: fg_amp
     real(dp), allocatable, dimension(:,:)     :: temp_amp
  end type genvec

  logical(lgt), allocatable, dimension(:,:)    :: fix_temp

contains

  subroutine initialize_genvec_mod(paramfile)
    implicit none

    character(len=128), intent(in) :: paramfile
    
    integer(i4b) :: i, j, nside
    logical(lgt) :: polarization, sample_template, sample_T_modes, ok, enforce_zero_cl
    character(len=2)   :: i_text, j_text
    character(len=128) :: param_text, comp_type

    call get_parameter(paramfile, 'LMAX',         par_int=lmax)    
    numcomp = lmax2numcomp(lmax)

    call get_parameter(paramfile, 'NSIDE',           par_int=nside)    
    call get_parameter(paramfile, 'NUMBAND',         par_int=numband)    
    call get_parameter(paramfile, 'SAMPLE_TT_MODES', par_lgt=sample_T_modes)
    call get_parameter(paramfile, 'ENFORCE_ZERO_CL', par_lgt=enforce_zero_cl)
    npix = 12*nside**2

    call get_parameter(paramfile, 'POLARIZATION',     par_lgt=polarization)
    if ((.not. polarization) .and. (.not. sample_T_modes)) then
       write(*,*) 'ERROR: Neither POLARIZATION nor SAMPLE_TT_MODES are selected. Aborting.'
       stop
    end if

    if (polarization) then
       nmaps = 3
    else
       nmaps = 1
    end if

    call get_parameter(paramfile, 'NUM_FG_SIGNAL_COMPONENTS', par_int=num_signal)
    sample_fg_pix = num_signal > 0

    call get_parameter(paramfile, 'NUM_FG_TEMP', par_int=num_fg_temp)
    sample_fg_temp = num_fg_temp > 0
    if (sample_T_modes) num_fg_temp = num_fg_temp + 4

    ! Find the number of constraints on the free template amplitudes
    allocate(fix_temp(num_fg_temp,numband))
    do i = 1, numband
       call int2string(i, i_text)

       if (sample_T_modes) then
          param_text = 'FIX_MONOPOLE' // i_text
          call get_parameter(paramfile, trim(param_text), par_lgt=fix_temp(1,i))
          param_text = 'FIX_DIPOLE' // i_text
          call get_parameter(paramfile, trim(param_text), par_lgt=fix_temp(2,i))
          fix_temp(3:4,i) = fix_temp(2,i)

          if (sample_fg_temp) then
             do j = 5, num_fg_temp
                call int2string(j-4, j_text)
                param_text = 'FIX_TEMP' // i_text // '_' // j_text
                call get_parameter(paramfile, trim(param_text), par_lgt=fix_temp(j,i))
             end do
          end if

       else
          do j = 1, num_fg_temp
             call int2string(j, j_text)
             param_text = 'FIX_TEMP' // i_text // '_' // j_text
             call get_parameter(paramfile, trim(param_text), par_lgt=fix_temp(j,i))
          end do
       end if
    end do

    reduce_pix = .false.
    do i = 1, num_signal
       call int2string(i, i_text)
       param_text= 'ENFORCE_POSITIVE_AMPLITUDE' // i_text
       call get_parameter(paramfile, param_text, par_lgt=ok)
       if (ok) cycle
       param_text= 'COMP_TYPE' // i_text
       call get_parameter(paramfile, param_text, par_string=comp_type)
       if (trim(comp_type)=='cmb' .and. .not. enforce_zero_cl) cycle
       reduce_pix = .true.
    end do

  end subroutine initialize_genvec_mod


  ! ****************************
  ! *          Utilities       *
  ! ****************************

  function genvec_dot_product(u,v)
    implicit none

    type(genvec), intent(in) :: u, v
    real(dp)                 :: genvec_dot_product

    integer(i4b) :: i, j, k, l, m

    genvec_dot_product = 0.d0
    do j = 1, nmaps
       do i = 1, numcomp
          genvec_dot_product = genvec_dot_product + u%cmb_amp(i,j) * v%cmb_amp(i,j)
       end do
    end do

    if (num_fg_temp > 0) then
       genvec_dot_product = genvec_dot_product + sum(u%temp_amp * v%temp_amp)
    end if

    if (sample_fg_pix .and. reduce_pix) then
       genvec_dot_product = genvec_dot_product + sum(u%fg_amp * v%fg_amp)
    end if

  end function genvec_dot_product


  subroutine allocate_genvec(vector)
    implicit none

    type(genvec), intent(inout) :: vector

    allocate(vector%cmb_amp(numcomp, nmaps))
    vector%cmb_amp          = 0.d0

    if (num_fg_temp > 0) then
       allocate(vector%temp_amp(num_fg_temp, numband))
       vector%temp_amp    = 0.d0
    end if

    if (sample_fg_pix) then
       allocate(vector%fg_amp(0:npix-1, nmaps, num_signal))
       vector%fg_amp        = 0.d0
    end if

  end subroutine allocate_genvec


  subroutine deallocate_genvec(vector)
    implicit none

    type(genvec), intent(inout) :: vector

    if (allocated(vector%cmb_amp))             deallocate(vector%cmb_amp)
    if (allocated(vector%fg_amp))              deallocate(vector%fg_amp)
    if (allocated(vector%temp_amp))            deallocate(vector%temp_amp)

  end subroutine deallocate_genvec



  ! ****************************
  ! *          Operator        *
  ! ****************************

  subroutine genvec_set_equal(v1, v2)
    implicit none

    type(genvec), intent(in)    :: v1
    type(genvec), intent(inout) :: v2

    v2%cmb_amp                             = v1%cmb_amp
    if (num_fg_temp > 0)      v2%temp_amp  = v1%temp_amp
    if (sample_fg_pix)        v2%fg_amp    = v1%fg_amp

  end subroutine genvec_set_equal

  subroutine genvec_plus(v1, v2, v3) 
    implicit none

    type(genvec), intent(in)    :: v1, v2
    type(genvec), intent(inout) :: v3

    v3%cmb_amp                            = v1%cmb_amp             + v2%cmb_amp
    if (num_fg_temp > 0)      v3%temp_amp = v1%temp_amp            + v2%temp_amp
    if (sample_fg_pix .and. reduce_pix)        v3%fg_amp   = v1%fg_amp              + v2%fg_amp

  end subroutine genvec_plus


  subroutine genvec_minus(v1, v2, v3) 
    implicit none

    type(genvec), intent(in)    :: v1, v2
    type(genvec), intent(inout) :: v3

    v3%cmb_amp                            = v1%cmb_amp  - v2%cmb_amp
    if (num_fg_temp > 0)      v3%temp_amp = v1%temp_amp - v2%temp_amp
    if (sample_fg_pix .and. reduce_pix)        v3%fg_amp   = v1%fg_amp   - v2%fg_amp

  end subroutine genvec_minus


  subroutine genvec_times(a, v1, v2) 
    implicit none

    real(dp),     intent(in)    :: a
    type(genvec), intent(in)    :: v1
    type(genvec), intent(inout) :: v2

    v2%cmb_amp                            = a * v1%cmb_amp      
    if (num_fg_temp > 0)      v2%temp_amp = a * v1%temp_amp
    if (sample_fg_pix .and. reduce_pix)        v2%fg_amp   = a * v1%fg_amp       

  end subroutine genvec_times


  subroutine genvec_plus_times(v1, alpha, v2, v3) 
    implicit none

    real(dp),     intent(in)    :: alpha
    type(genvec), intent(in)    :: v1, v2
    type(genvec), intent(inout) :: v3

    v3%cmb_amp                            = v1%cmb_amp  + alpha * v2%cmb_amp
    if (num_fg_temp > 0)      v3%temp_amp = v1%temp_amp + alpha * v2%temp_amp
    if (sample_fg_pix .and. reduce_pix)        v3%fg_amp   = v1%fg_amp   + alpha * v2%fg_amp

  end subroutine genvec_plus_times


  subroutine nullify_genvec(v)
    implicit none

    type(genvec), intent(inout)    :: v

    v%cmb_amp                            = 0.d0
    if (num_fg_temp > 0)      v%temp_amp = 0.d0
    if (sample_fg_pix)        v%fg_amp   = 0.d0

  end subroutine nullify_genvec

  subroutine genvec2linvec(a, b)
    implicit none

    type(genvec),                        intent(in) :: a
    real(dp),     allocatable, dimension(:)         :: b

    integer(i4b) :: n, m

    if (.not. allocated(b)) then
       n = size(a%cmb_amp)
       if (num_fg_temp > 0)  n = n + num_fg_temp
       if (sample_fg_pix)    n = n + size(a%fg_amp)
       allocate(b(n))
    end if
    b = 0.d0

    n = 1
    m = size(a%cmb_amp)
    b(n:n+m-1) = reshape(a%cmb_amp, shape(b))
    n = n + m

    if (num_fg_temp > 0) then
       m = num_fg_temp
       b(n:n+m-1) = reshape(a%temp_amp, shape(b))
       n = n + m
    end if

    if (sample_fg_pix) then
       m = size(a%fg_amp)
       b(n:n+m-1) = reshape(a%fg_amp, shape(b))
       n = n + m
    end if

  end subroutine genvec2linvec

  subroutine linvec2genvec(b, a)
    implicit none

    type(genvec)                             :: a
    real(dp),     allocatable, dimension(:)  :: b

    integer(i4b) :: n, m

    call nullify_genvec(a)

    n = 1
    m = size(a%cmb_amp)
    a%cmb_amp = reshape(b(n:n+m-1), shape(a%cmb_amp))
    n = n + m

    if (num_fg_temp > 0) then
       m = num_fg_temp
       a%temp_amp = reshape(b(n:n+m-1), shape(a%temp_amp))
       n = n + m
    end if

    if (sample_fg_pix) then
       m = size(a%fg_amp)
       a%fg_amp = reshape(b(n:n+m-1), shape(a%fg_amp))
       n = n + m
    end if

  end subroutine linvec2genvec

  subroutine bcast_genvec(comm, v)
    implicit none

    integer(i4b)  :: comm, root=0, ierr
    type(genvec)  :: v

    if (allocated(v%cmb_amp))              call mpi_bcast(v%cmb_amp,              size(v%cmb_amp), &
             & MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (allocated(v%fg_amp))               call mpi_bcast(v%fg_amp,               size(v%fg_amp), &
             & MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (allocated(v%temp_amp))             call mpi_bcast(v%temp_amp,             size(v%temp_amp), &
             & MPI_DOUBLE_PRECISION, root, comm, ierr)

  end subroutine bcast_genvec

  subroutine reduce_genvec(comm, v1, v2)
    implicit none

    integer(i4b)  :: comm, root=0, ierr
    type(genvec)             :: v1
    type(genvec), optional   :: v2

    if (present(v2)) then
       if (allocated(v1%cmb_amp))              call mpi_reduce(v1%cmb_amp,             v2%cmb_amp, &
            & size(v1%cmb_amp),             MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr)
       if (allocated(v1%fg_amp) .and. reduce_pix)               call mpi_reduce(v1%fg_amp,              v2%fg_amp, &
            & size(v1%fg_amp),              MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr)
       if (num_fg_temp > 0)                 call mpi_reduce(v1%temp_amp, v2%temp_amp, &
            & size(v1%temp_amp), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr)
    else
       if (allocated(v1%cmb_amp))              call mpi_reduce(MPI_IN_PLACE, v1%cmb_amp, &
            & size(v1%cmb_amp),             MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr)
       if (allocated(v1%fg_amp) .and. reduce_pix)               call mpi_reduce(MPI_IN_PLACE, v1%fg_amp, &
            & size(v1%fg_amp),              MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr)
       if (allocated(v1%temp_amp))             call mpi_reduce(MPI_IN_PLACE, v1%temp_amp, &
            & size(v1%temp_amp), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr)
    end if

  end subroutine reduce_genvec
  
end module comm_genvec_mod
