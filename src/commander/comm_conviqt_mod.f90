module comm_conviqt_mod
  use sharp
  use healpix_types
  use pix_tools
  use iso_c_binding, only : c_ptr, c_double, c_int
  use comm_map_mod
  use comm_utils
  use comm_shared_arr_mod
  implicit none

  private
  public comm_conviqt, conviqt_ptr

  type :: comm_conviqt
    integer(i4b) :: lmax, mmax, nmaps, bmax, nside, npix, comm, optim, psisteps, win
    real(dp), allocatable, dimension(:)        :: lnorm
    type(shared_2d_sp) :: c
    type(shared_1d_int) :: pixLookup
    real(sp) :: psires
    class(comm_mapinfo), pointer :: info => null()
    type(shared_2d_spc) :: alm_beam
  contains
    procedure     :: interp
    procedure     :: precompute_sky
    procedure     :: get_alms
    procedure     :: dealloc
  end type comm_conviqt

  interface comm_conviqt
    procedure constructor
  end interface comm_conviqt

  type conviqt_ptr
     class(comm_conviqt), pointer :: p => null()
  end type conviqt_ptr

contains

  !Constructor

  function constructor(myid_shared, comm_shared, myid_inter, comm_inter, nside, lmax, nmaps, bmax, beam, map, optim)
    implicit none
    integer(i4b),                 intent(in) :: nside, lmax, bmax, nmaps
    integer(i4b),                 intent(in) :: myid_shared, comm_shared, myid_inter, comm_inter
    class(comm_map),              intent(in) :: beam, map
    integer(i4b),                 intent(in) :: optim ! desired optimization flags
    class(comm_conviqt), pointer             :: constructor

    integer(i4b) :: i, j, k, l, m, ierr, nalm_tot, nalm
    integer(i4b), allocatable, dimension(:)   :: ind
    complex(spc), allocatable, dimension(:,:) :: alm
    real(dp) :: theta, phi


    allocate(constructor)
    constructor%lmax  = lmax
    constructor%mmax  = min(beam%info%mmax,lmax)
    constructor%bmax  = bmax
    constructor%nside = nside
    constructor%nmaps = nmaps
    constructor%npix  = 12*nside**2
    constructor%comm  = map%info%comm
    constructor%optim = optim
    nalm_tot          = (lmax+1)*(lmax+2)/2

    !current optimization levels:
    !0 - no optimization
    !1 - use lower resolution psi grid
    !2 - don't interpolate between samples, just use closest

    !make this 2^n
    constructor%psisteps = 2*constructor%bmax
    constructor%psires   = 2*pi/constructor%psisteps

    allocate(constructor%lnorm(0:lmax))
    do l = 0, lmax
       constructor%lnorm(l) = 0.5d0*sqrt(4.d0*pi/(2.d0*l+1.d0))
    end do

    constructor%info => comm_mapinfo(map%info%comm, nside, lmax, nmaps, nmaps==3)
 
    !Set up shared beam
    nalm = count(constructor%info%lm(2,:) >= 0)  

    allocate(ind(nalm), alm(nmaps,nalm))
    alm = 0
    j = 1
    do i = 0, constructor%info%nalm-1
       l = constructor%info%lm(1,i)
       m = constructor%info%lm(2,i)
       if (m < 0) cycle
       call beam%info%lm2i(l,m,k)
       if (k == -1) cycle
       ind(j) = (l-1)*l/2 + l + m + 1
       if (m == 0) then
          alm(:,j) = cmplx(real(beam%alm(k,:),sp),0.0) 
       else
          alm(:,j) = cmplx(real(beam%alm(k,:),sp),real(beam%alm(k+1,:),sp))/sqrt(2.0) 
       end if
       j = j+1
    end do
    call init_shared_2d_spc(myid_shared, comm_shared, myid_inter, comm_inter, &
         & [nmaps,nalm_tot], constructor%alm_beam)
    call sync_shared_2d_spc_alm(constructor%alm_beam, ind, alm)

    deallocate(ind, alm)

    ! Precompute convolution cube
    call init_shared_2d_sp(myid_shared, comm_shared, myid_inter, comm_inter, &
         & [constructor%npix,constructor%psisteps], constructor%c)

    call constructor%precompute_sky(map)

    !if (map%info%myid == 0) then
    !   do i = 0, constructor%psisteps-1
    !      call int2string(i,id)
    !      call write_map2('c'//id//'.fits', constructor%c%a(:,i:i+1))
    !   end do    
    !end if

    !precompute pixel->pixel lookup table
    call init_shared_1d_int(myid_shared, comm_shared, myid_inter, comm_inter, &
        & [12*(map%info%nside**2)], constructor%pixLookup)

    if(myid_shared == 0) then
      do i=0, 12*(map%info%nside**2) -1
        call pix2ang_ring(map%info%nside, i, theta, phi)
        call ang2pix_ring(constructor%nside, theta, phi, j)
        constructor%pixLookup%a(i+1) = j
      end do
    end if

    call mpi_win_fence(0, constructor%pixLookup%win, ierr)    

  end function constructor


  !interpolates the precomputed map to the psi angle
  function interp(self, pix, psi)
    implicit none
    class(comm_conviqt),       intent(in) :: self
    integer(i4b),              intent(in) :: pix
    real(dp),                  intent(in) :: psi
    real(sp)                              :: interp

    real(sp)       :: x0, x1, unwrap
    integer(i4b)   :: pixnum, psii, psiu, bpsi

    ! Get pixel number
    pixnum = pix  !self%pixLookup%a(pix+1)

    !unwrap psi
    unwrap = modulo(-real(psi,sp), real(2.d0*pi,sp))
    !unwrap = modulo(psi, 2.d0*pi)

    if (self%optim == 2) then
       bpsi = max(nint(unwrap / self%psires),0)
       if (bpsi == self%psisteps) bpsi = 0
      interp = self%c%a(pixnum+1, bpsi+1)
      if (isNaN(interp)) then
         write(*,*) 'nan in sl interpolation', pix, psi, pixnum, unwrap, interp
         stop
      end if
      !write(*,*) "Interp: ", interp
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
    interp = (self%c%a(pixnum+1, psii+1) * (x1 - unwrap) + &
         & self%c%a(pixnum+1, psiu+1) * (unwrap - x0))/(x1 - x0)

  end function interp

  !precomputes the maps and stores to the allocated array
  !TODO: handle polarized sidelobes?
  subroutine precompute_sky(self, map)
    implicit none
    class(comm_conviqt),                          intent(inout) :: self
    class(comm_map),                              intent(in)    :: map ! Must contain alms

    integer(i4b) :: i, j, np, ierr
    integer*8    :: fft_plan
    real(dp),       allocatable, dimension(:,:)   :: marray
    real(c_double), allocatable, dimension(:,:) :: alm    
    real(c_double), allocatable, dimension(:,:)   :: mout
    real(dp),       allocatable, dimension(:)     :: dt
    complex(dpc),   allocatable, dimension(:)     :: dv

    np = self%info%np

    allocate(marray(np, -self%bmax:self%bmax))
    marray = 0

    allocate(alm(0:self%info%nalm-1, 2))
    allocate(mout(0:self%info%np-1, 2))

    ! Compute maps for each m
    do j=0, self%bmax

       call self%get_alms(j, map, alm)

       if (j == 0) then
          call sharp_execute(SHARP_Y, 0, 1, alm(:,1:1), self%info%alm_info, &
               & mout(:,1:1), self%info%geom_info_T, comm=self%comm)
          marray(:, 0) = mout(:,1)
       else
          call sharp_execute(SHARP_Y, j, 2, alm(:,1:2), self%info%alm_info, &
               & mout(:,1:2), self%info%geom_info_T, comm=self%comm)

          ! Rotate by 90 degrees in psi
!!$          sign1 = 1; sign2 = -1; quadrant = modulo(j,4)
!!$          if (iand(quadrant,1)) then
!!$             marray(:,j) = mout(:,1)
!!$             mout(:,1)   = mout(:,2)
!!$             mout(:,2)   = marray(:,j)
!!$             sign2       = sign1
!!$             sign1       = -sign1
!!$          end if
!!$          if (quadrant==1 .or. quadrant == 2) sign1 = -sign1
!!$          if (quadrant==2 .or. quadrant == 3) sign2 = -sign2
!!$          if (sign1 /= 1) mout(:,1) = -mout(:,1)
!!$          if (sign2 /= 1) mout(:,2) = -mout(:,2)

          marray(:, j) = mout(:,1)
          marray(:,-j) = mout(:,2)

       end if
      if(map%info%myid == 0 .and. .false.) then
        write(*,*) 'Sharp  execute, m=', j, sum(mout(:,1:2)), sum(alm(:,1:2))
      end if
      
    end do

    ! Fourier transform in psi direction
    allocate(dt(self%psisteps), dv(0:self%psisteps/2))
    call dfftw_plan_dft_c2r_1d(fft_plan, self%psisteps, dv, dt, fftw_estimate + fftw_unaligned)
    if(fft_plan == 0) then
      write(*,*) 'Failed to create fftw plan, thread ', map%info%myid
    end if
    do i=1, np

      !do fft of data, store to dt
      dv(0)  = marray(i, 0)
      do j=1, self%bmax
        dv(j) = dcmplx(marray(i, j), marray(i, -j))
      end do

      call dfftw_execute_dft_c2r(fft_plan, dv, dt)

      self%c%a(self%info%pix(i)+1,:) = real(dt(1:self%psisteps),sp)

    end do
    call mpi_win_fence(0, self%c%win, ierr)
 
    deallocate(marray, alm, mout) 
 
    call dfftw_destroy_plan(fft_plan)
    deallocate(dt)
    deallocate(dv)
 
  end subroutine precompute_sky

  ! Note: lnorm is set to 1 here. Todo: Figure out why 
  subroutine get_alms(self, m_b, map, alm)
    implicit none
    class(comm_conviqt),        intent(in)  :: self
    integer(i4b),               intent(in)  :: m_b
    class(comm_map),            intent(in)  :: map
    real(dp), dimension(0:,1:), intent(out) :: alm

    integer(i4b) :: i, ineg, j, l, m
    real(dp)     :: spinsign, mfac, sqrt_two
    complex(dpc) :: v1, v2, alm_b(3), alm_s(3), almc

    spinsign = 1; if (m_b /= 0) spinsign = -1
    sqrt_two = sqrt(2.d0)

    alm = 0.d0
    do i = 0, self%info%nalm-1
        l  = self%info%lm(1,i)
        m  = self%info%lm(2,i)
        if (m < 0) cycle ! Do negative m at the same time as positive m
        call self%info%lm2i(l,-m,ineg)
        if (l < m_b) cycle
        mfac = 1; if (iand(m,1) /= 0) mfac = -1 

        j = (l-1)*l/2 + l + m_b + 1
        alm_b = self%alm_beam%a(:,j)
        call map%get_alm_TEB( l, m,   .true., alm_s)
        !if (l == 0 .and. m == 0) write(*,*) alm_b, alm_s
        v1 = sum(      alm_s *alm_b)
        !THIS WAS THE LINE THAT WAS WRONG
        !I FOUND IT
        !HAHAHAHAHAHA
        !v2 = 0; if (m_b /= 0) 
        v2 = sum(conjg(alm_s)*alm_b)*mfac
        
        ! Positive spin
        almc = spinsign*self%lnorm(l) * (v1+conjg(v2)*mfac)
        if (m == 0) then
           alm(i,1) = real(almc)
        else
           alm(i,1) = real(almc) * sqrt_two
           alm(ineg,1) = imag(almc) * sqrt_two
        end if

        if (m_b > 0) then
           ! Negative spin
           almc = -dcmplx(0,1)*spinsign*self%lnorm(l) * (v1-conjg(v2)*mfac)
           if (m == 0) then
              alm(i,2) = real(almc)
           else
              alm(i,2) = real(almc) * sqrt_two
              alm(ineg,2) = imag(almc) * sqrt_two
           end if
        end if
     end do

  end subroutine get_alms


  !deallocate memory used by the class
  subroutine dealloc(self)
    implicit none
    class(comm_conviqt), intent(inout) :: self

    deallocate(self%lnorm)
    call dealloc_shared_2d_spc(self%alm_beam)
    call dealloc_shared_2d_sp(self%c)
    call dealloc_shared_1d_int(self%pixLookup)

  end subroutine dealloc

end module comm_conviqt_mod
