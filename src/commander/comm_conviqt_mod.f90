module comm_conviqt_mod
  use sharp
  use healpix_types
  use pix_tools
  use iso_c_binding, only : c_ptr, c_double, c_int
  use comm_map_mod
  use comm_utils
  implicit none

  private
  public comm_conviqt, conviqt_ptr

  type :: comm_conviqt
    integer(i4b) :: lmax, mmax, bmax, nside, comm, optim, psisteps, win
    real(dp), allocatable, dimension(:,:)      :: c
    real(dp) :: psires
  contains
    procedure     :: interp
    procedure     :: precompute_sky
    procedure     :: dealloc
  end type comm_conviqt

  interface comm_conviqt
    procedure constructor
  end interface comm_conviqt

  type conviqt_ptr
     class(comm_conviqt), pointer :: p
  end type conviqt_ptr

contains

  !Constructor

  function constructor(beam, map, optim)
    implicit none
    class(comm_map),              intent(inout) :: beam, map
    integer(i4b),                 intent(in)    :: optim ! desired optimization flags
    class(comm_conviqt), pointer                :: constructor
    integer(i4b)                                :: ierr, winsize, disp_unit
    type(c_ptr)                                 :: baseptr
    integer(i4b), allocatable, dimension(:)     :: arrayshape

    !inquire(file='precomp.unf', exist=exist)
    !if (exist) then
    !   open(59,file='precomp.unf', form='unformatted')
    !   read(59) array
    !   close(59)
    !else
       

    !   open(59,file='precomp.unf', form='unformatted')
    !   write(59) array
    !   close(59
    !end if

    allocate(constructor)

    constructor%lmax  = min(beam%info%lmax, map%info%lmax)
    constructor%mmax  = beam%info%mmax
    constructor%bmax  = constructor%lmax
    constructor%nside = map%info%nside
    constructor%comm  = map%info%comm
    constructor%optim = optim

    !current optimization levels:
    !0 - no optimization
    !1 - use lower resolution psi grid
    !2 - don't interpolate between samples, just use closest

    !make this 2^n
    constructor%psires = 2.d0*pi/constructor%bmax
    constructor%psisteps = int(2.d0*pi/constructor%psires)

    !TODO: move this to a more useful spot so it can be used elsewhere
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, constructor%comm, ierr)

    winsize = 0
    if(map%info%myid == 0) then
      winsize = (12*map%info%nside**2 *constructor%psisteps)*8
    end if
    disp_unit = 1
    call MPI_Win_allocate_shared(winsize, disp_unit, MPI_INFO_NULL, constructor%comm, baseptr, constructor%win, ierr)

    if (map%info%myid /= 0) then
      call MPI_Win_shared_query(constructor%win, 0, winsize, disp_unit, baseptr, ierr)
    end if

    allocate(arrayshape(2))
    arrayshape(1) = 12*map%info%nside**2
    arrayshape(2) = constructor%psisteps

    call C_F_POINTER(baseptr, constructor%c, arrayshape)

    call MPI_WIN_FENCE(0, constructor%win, ierr)

    deallocate(arrayshape)

    !allocate(constructor%c(12*map%info%nside**2, constructor%psisteps))

    call constructor%precompute_sky(beam, map)

  end function constructor


  !interpolates the precomputed map to the psi angle
  function interp(self, theta, phi, psi)
    implicit none
    class(comm_conviqt),       intent(in) :: self
    real(dp),                  intent(in) :: theta, phi, psi
    real(sp)                              :: interp

    real(sp)       :: x0, x1, unwrap
    integer(i4b)   :: pixnum, psii, psiu

    ! Get pixel number
    call ang2pix_ring(self%nside, theta, modulo(phi + pi, 2.d0 *pi), pixnum)

    !unwrap psi
    unwrap = modulo(psi, 2.d0*pi)
    !if (psi < 0) then
    !  unwrap = psi + 2*pi
    !else
    !  unwrap = psi
    !end if

    if (self%optim == 2) then
      interp = self%c(pixnum+1, nint(unwrap / self%psires) + 1)
      return
    end if

    !index of the lower bound of psi in the array
    psii = int(unwrap / self%psires)
    !index of the upper bound of psi
    psiu = psii + 1
    if (psiu >= self%psisteps) then !catch rollover
      psiu = 0
    end if
    
    !linear interpolation between evaluated samples in psi
    ! y = (y_0 (x_1 - x) + y_1 (x - x_0))/(x_1 - x_0)
    x0     = psii * self%psires
    x1     = psiu * self%psires
    interp = ((self%c(pixnum+1, psii+1)) * (x1 - unwrap) + self%c(pixnum+1, psiu+1) * (unwrap - x0))/(x1 - x0)

    !interp = 0.d0

    !if( isNaN(interp)) then
    !  write(*,*) self%c(pixnum+1, psii), self%c(pixnum+1, psiu), interp, theta, phi,psi
    !end if

  end function interp

  !precomputes the maps and stores to the allocated array
  !TODO: handle polarized sidelobes?
  subroutine precompute_sky(self, beam, map)
    implicit none
    class(comm_conviqt),                          intent(inout) :: self
    class(comm_map),                              intent(in)    :: beam, map ! Must contain alms

    integer(i4b) :: i,j,k,np,l,m_s, m_b, pixNum, q,u, ierr
    real(dp)     :: theta, phi, psi, one_over_root_two
    real(dp)     :: alm_b_r, alm_s_r, alm_b_c, alm_s_c
    integer*8    :: fft_plan
    real(dp),       allocatable, dimension(:,:)   :: marray
    real(c_double), allocatable, dimension(:,:) :: alm    
    real(c_double), allocatable, dimension(:,:)   :: mout
    real(dp),       allocatable, dimension(:)     :: dt
    double complex,   allocatable, dimension(:)     :: dv

    if(map%info%myid == 0) then
      self%c = 0
    end if

    np = map%info%np

    allocate(marray(np, -self%bmax:self%bmax))
    marray = 0

    allocate(alm(0:map%info%nalm-1, 2))
    alm = 0
    allocate(mout(0:map%info%np-1, 2))

    !check lmaxes
    if (beam%info%lmax > map%info%lmax) then
      write(*,*) "possible problem with not enought ms in precompute_sky"
      stop
    endif

    !copy to alm array 
    one_over_root_two = 1.d0/sqrt(2.d0)

    !m_b = 0 case
    do i=0, map%info%nalm-1
      l   = map%info%lm(1,i)
      m_s = map%info%lm(2,i)
      !write(*,*) i, l, m_s, q
      alm_s_r = map%alm(i, 1)
      if(m_s /= 0) then
        alm_s_r = alm_s_r * one_over_root_two
      end if
      call beam%info%lm2i(l, 0, q)
      if(q == -1) then
        alm_b_r = 0.d0
      else
        alm_b_r = beam%alm(q, 1)
      end if
      !if(m_s == 0) then
      !  write(*, *) l, 0, alm_b_r
      !end if
      alm(i, 1) = -1 ** m_s * sqrt(4*pi/(2*l + 1)) * (alm_b_r * alm_s_r) !real part
    end do
    mout = 0
    !if(map%info%myid == 0) then
    !  write(*,*) 'About to call sharp execute m=0', map%info%nalm
    !end if

     !Compute m = 0 case
    call sharp_execute(SHARP_Y, 0, 1, alm(:, 1:1), map%info%alm_info, mout(:,1:1), &
         & map%info%geom_info_T, comm=self%comm)
    marray(:,0) = mout(:,1)

    if(map%info%myid == 0 .and. .false.) then
      write(*,*) 'Sharp  execute, m=', 0, sum(mout(:,1:1)), sum(alm(:,1:1)), alm_b_r, alm_s_r
    end if

    ! Compute m /= 0 cases
    do j=1, beam%info%mmax
      alm = 0
 
      do i=0, map%info%nalm-1
        l   = map%info%lm(1,i)
        m_s = map%info%lm(2,i)
        if(j > l) cycle !catch case where m_b > l_sky
        alm_s_r = map%alm(i, 1)
        call map%info%lm2i(l, -m_s, q)
        alm_s_c = map%alm(q, 1)
        if(m_s /= 0) then
          alm_s_r = alm_s_r * one_over_root_two
          alm_s_c = alm_s_c * one_over_root_two
        end if

        call beam%info%lm2i(l, j, q)
        if(q < 0) then
          alm_b_r = 0.d0
        else
          alm_b_r = beam%alm(q, 1) * one_over_root_two
        end if

        call beam%info%lm2i(l, -j, u)
        if(u < 0) then
          alm_b_c = 0.d0
        else
          alm_b_c = beam%alm(u, 1) * one_over_root_two
        end if

        !if(j == 10 .and. m_s == 0) then
        !  write(*,*) l, alm_b_r, alm_b_c
        !end if

        alm(i, 1) = -1 ** m_s * sqrt(4*pi/(2*l + 1)) * ((alm_b_r * alm_s_r) + (alm_b_c * alm_s_c)) !real part
        alm(i, 2) = -1 ** m_s * sqrt(4*pi/(2*l + 1)) * ((alm_b_r * alm_s_c) - (alm_b_c * alm_s_r)) !imag part
    
      end do

      mout = 0.d0
      !if(map%info%myid == 0) then
      !  write(*,*) 'Before execute, m=', j, sum(mout(:,1:2)), sum(alm(:,1:2)) 
      !end if
      call sharp_execute(SHARP_Y, j, 2, alm(:,1:2), map%info%alm_info, mout(:,1:2), &
           & map%info%geom_info_T, comm=self%comm)
      if(map%info%myid == 0 .and. .false.) then
        write(*,*) 'Sharp  execute, m=', j, sum(mout(:,1:2)), sum(alm(:,1:2))
      end if
      marray(:,-j) = mout(:,2)
      marray(:, j) = mout(:,1)
    end do

    !if(map%info%myid == 0) then
    !  write(*,*) 'Done sharp execute', sum(marray)
    !end if

    ! Fourier transform in psi direction
    allocate(dt(self%psisteps), dv(0:self%psisteps/2))
    call dfftw_plan_dft_c2r_1d(fft_plan, self%psisteps, dv, dt, fftw_estimate + fftw_unaligned)
    if(fft_plan == 0) then
      write(*,*) 'Failed to create fftw plan, thread ', map%info%myid
    end if
    self%c = 0.d0
    do i=1, np
      !do fft of data, store to dt
      pixNum = map%info%pix(i)
      dv = 0
      dt = 0
      dv(0)  = marray(i, 0)
      do j=1, self%mmax
        dv(j) = cmplx(marray(i, j), marray(i, -j))
      end do

      !if(map%info%myid == 6 .and. i == 0) then
      !  write(*,*) size(dv), self%mmax, size(marray(0,:))
      !end if


      call dfftw_execute_dft_c2r(fft_plan, dv, dt)
      !if(any(isNaN(dt))) then
      !  write(*,*) dv, dt, map%info%myid
      !end if

      dt = 2.d0 * pi / (real(self%psisteps)**2) * dt

      !copy to c
      self%c(pixNum+1,:) = dt(1:self%psisteps)
      if(map%info%myid == 2 .and. .false.) then
        write(*,*) 'Executed fft, i = ', i, sum(dv), sum(dt), sum(self%c), size(self%c(pixNum+1,:)), size(dt), dv(0), sum(marray(:,0))
      end if

    end do
 
    deallocate(marray, alm, mout) 
 
    call dfftw_destroy_plan(fft_plan)
    deallocate(dt)
    deallocate(dv)
 
    !bcast to all cores

    !call mpi_allreduce(MPI_IN_PLACE, self%c, size(self%c), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)

    call MPI_WIN_FENCE(0, self%win, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if(map%info%myid == 0) then
      write(*,*) 'Allreduced, sum = ', sum(self%c)
    end if

  end subroutine precompute_sky

  !deallocate memory used by the class
  subroutine dealloc(self)
    implicit none
    class(comm_conviqt),                          intent(inout) :: self

    deallocate(self%c)
  end subroutine dealloc

end module comm_conviqt_mod
