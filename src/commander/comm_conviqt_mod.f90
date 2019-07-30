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
    integer(i4b) :: lmax, mmax, nmaps, bmax, nside, npix, comm, optim, psisteps
    real(dp), allocatable, dimension(:)        :: lnorm
    real(dp), allocatable, dimension(:,:)      :: c
    real(dp) :: psires
    class(comm_mapinfo), pointer :: info
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
     class(comm_conviqt), pointer :: p
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

    logical(lgt) :: exist
    integer(i4b) :: i, j, k, l, m, ierr, nalm_tot, nalm
    integer(i4b), allocatable, dimension(:)   :: ind
    complex(spc), allocatable, dimension(:,:) :: alm
    character(len=4) :: id

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

    constructor%lmax  = lmax
    constructor%mmax  = min(beam%info%mmax,lmax)
    constructor%bmax  = bmax
    constructor%nside = nside
    constructor%nmaps = nmaps
    constructor%npix  = 12*nside**2
    constructor%comm  = map%info%comm
    constructor%optim = optim
    nalm_tot          = (lmax+1)*(lmax+2)/2
!    if(map%info%myid == 0) write(*,*) constructor%lmax, constructor%bmax, constructor%nside, constructor%mmax, constructor%optim

 !   call mpi_finalize(ierr)
  !  stop

    !current optimization levels:
    !0 - no optimization
    !1 - use lower resolution psi grid
    !2 - don't interpolate between samples, just use closest

    !make this 2^n
    constructor%psisteps = 2*constructor%bmax
    constructor%psires   = 2.d0*pi/constructor%psisteps

    allocate(constructor%lnorm(0:lmax))
    do l = 0, lmax
       !constructor%lnorm(l) = 0.5d0*sqrt(4.d0*pi/(2.d0*l+1.d0))
       constructor%lnorm(l) = 1.d0 
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
          alm(:,j) = beam%alm(k,:) 
       else
          alm(:,j) = cmplx(beam%alm(k,:),beam%alm(k+1,:))/sqrt(2.d0) 
       end if
       j = j+1
    end do
    call init_shared_2d_spc(myid_shared, comm_shared, myid_inter, comm_inter, &
         & [nmaps,nalm_tot], constructor%alm_beam)
    call sync_shared_2d_spc_alm(constructor%alm_beam, ind, alm)

!    write(*,*) constructor%info%myid, minval(abs(constructor%alm_beam%a)), maxval(abs(constructor%alm_beam%a))
    deallocate(ind, alm)

    ! Precompute convolution cube
    allocate(constructor%c(0:constructor%npix-1, 0:constructor%psisteps-1))

    call constructor%precompute_sky(beam, map)
!!$    if (map%info%myid == 0) then
!!$       do i = 0, constructor%psisteps-1
!!$          call int2string(i,id)
!!$          call write_map2('c'//id//'.fits', constructor%c(:,i:i))
!!$       end do
!!$    end if
!!$    call mpi_finalize(ierr)
!!$    stop

  end function constructor


  !interpolates the precomputed map to the psi angle
  function interp(self, theta, phi, psi)
    implicit none
    class(comm_conviqt),       intent(in) :: self
    real(dp),                  intent(in) :: theta, phi, psi
    real(sp)                              :: interp

    real(sp)       :: x0, x1, unwrap, bpsi
    integer(i4b)   :: pixnum, psii, psiu

    ! Get pixel number
    call ang2pix_ring(self%nside, theta, modulo(phi, 2.d0 *pi), pixnum)
!    call ang2pix_ring(self%nside, theta, modulo(phi + pi, 2.d0 *pi), pixnum)

    !unwrap psi
    unwrap = modulo(psi, 2.d0*pi)
    !if (psi < 0) then
    !  unwrap = psi + 2*pi
    !else
    !  unwrap = psi
    !end if

    if (self%optim == 2) then
       bpsi = max(nint(unwrap / self%psires),0)
       if (bpsi == self%psisteps) bpsi = 0
      interp = self%c(pixnum, bpsi)
      if (isNaN(interp)) then
         write(*,*) 'nan', theta, phi, psi, pixnum, unwrap, interp
      end if
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
    interp = ((self%c(pixnum, psii)) * (x1 - unwrap) + self%c(pixnum, psiu) * (unwrap - x0))/(x1 - x0)

    !interp = 0.d0

    !if( isNaN(interp)) then
    !  write(*,*) self%c(pixnum, psii), self%c(pixnum, psiu), interp, theta, phi,psi
    !end if

  end function interp

  !precomputes the maps and stores to the allocated array
  !TODO: handle polarized sidelobes?
  subroutine precompute_sky(self, beam, map)
    implicit none
    class(comm_conviqt),                          intent(inout) :: self
    class(comm_map),                              intent(in)    :: beam, map ! Must contain alms

    integer(i4b) :: i,j,k,np,l,m,m_s, m_b, pixNum, q,u, ierr
    integer(i4b) :: sign1, sign2, quadrant
    real(dp)     :: theta, phi, psi, one_over_root_two
    real(dp)     :: alm_b_r, alm_s_r, alm_b_c, alm_s_c
    integer*8    :: fft_plan
    real(dp),       allocatable, dimension(:,:)   :: marray
    real(c_double), allocatable, dimension(:,:) :: alm    
    real(c_double), allocatable, dimension(:,:)   :: mout
    real(dp),       allocatable, dimension(:)     :: dt
    complex(dpc),   allocatable, dimension(:)     :: dv

    self%c = 0

    np = self%info%np


    allocate(marray(np, -self%bmax:self%bmax))
    marray = 0

    allocate(alm(0:self%info%nalm-1, 2))
    allocate(mout(0:self%info%np-1, 2))

    !check lmaxes
!!$    if (beam%info%lmax > map%info%lmax) then
!!$      write(*,*) "possible problem with not enought ms in precompute_sky"
!!$      stop
!!$    endif

    ! Compute maps for each m
    do j=0, self%bmax

       call self%get_alms(j, map, alm)

       do k = 0, self%info%nalm-1
!          write(*,*) j, self%info%lm(:,k), alm(k,1)
       end do

       if (j == 0) then
          call sharp_execute(SHARP_Y, j, 1, alm(:,1:1), self%info%alm_info, &
               & mout(:,1:1), self%info%geom_info_T, comm=self%comm)
          marray(:, j) = mout(:,1)
       else
          call sharp_execute(SHARP_Y, j, 2, alm(:,1:2), self%info%alm_info, &
               & mout(:,1:2), self%info%geom_info_P, comm=self%comm)

!!$          ! Rotate by 90 degrees in psi
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

    if(map%info%myid == 0) then
!      write(*,*) 'Done sharp execute', sum(marray)
    end if

    ! Fourier transform in psi direction
    allocate(dt(self%psisteps), dv(0:self%psisteps/2))
    call dfftw_plan_dft_c2r_1d(fft_plan, self%psisteps, dv, dt, fftw_estimate + fftw_unaligned)
    if(fft_plan == 0) then
      write(*,*) 'Failed to create fftw plan, thread ', map%info%myid
    end if
    
    self%c = 0.d0
    do i=1, np
        if(map%info%myid == 0) then
!          write(*,*) 'i=', i, np
        end if


      !do fft of data, store to dt
      dv(0)  = marray(i, 0)
      do j=1, self%bmax
        dv(j) = cmplx(marray(i, j), marray(i, -j))
      end do

      !if(map%info%myid == 6 .and. i == 0) then
      !  write(*,*) size(dv), self%mmax, size(marray(0,:))
      !end if


      call dfftw_execute_dft_c2r(fft_plan, dv, dt)
      !if(any(isNaN(dt))) then
      !  write(*,*) dv, dt, map%info%myid
      !end if

!      dt = dt / self%psisteps

      !copy to c
      self%c(self%info%pix(i),0:self%psisteps -1) = dt(1:self%psisteps)
      !if(map%info%myid == 2 .and. sum(self%c) /= 0.d0) then
      !  write(*,*) 'Executed fft, i = ', i, sum(dv), sum(dt), sum(self%c), size(self%c(pixNum,:)), size(dt)
      !end if

    end do
 
    deallocate(marray, alm, mout) 
 
    call dfftw_destroy_plan(fft_plan)
    deallocate(dt)
    deallocate(dv)
 
    !bcast to all cores

    if(.false.) then
      write(*,*) ' Pre allreduce, sum = ', sum(self%c), size(self%c), map%info%myid
    end if


!    if(map%info%myid == 0) write(*,*) 'b'
    call mpi_allreduce(MPI_IN_PLACE, self%c, size(self%c), MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
!    if(map%info%myid == 0) write(*,*) 'c'

    if(map%info%myid == 0) then
!      write(*,*) 'Allreduced, sum = ', sum(self%c)
    end if

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

    do i = 0, self%info%nalm-1
        l  = self%info%lm(1,i)
        m  = self%info%lm(2,i)
        if (m < 0) cycle ! Do negative m at the same time as positive m
        call self%info%lm2i(l,-m,ineg)
        alm(i,   :) = 0.d0
        alm(ineg,:) = 0.d0
        if (l < m_b) cycle
        mfac = 1; if (iand(m,1) /= 0) mfac = -1 

        j = (l-1)*l/2 + l + m_b + 1
        alm_b = self%alm_beam%a(:,j)
        call map%get_alm_TEB( l, m,   .true., alm_s)
        v1 = sum(      alm_s *alm_b)
        v2 = 0; if (m_b /= 0) v2 = sum(conjg(alm_s)*alm_b)*mfac
        
!        if (self%info%myid==23) then
!!$        if (isNaN(abs(v1))) then
!           write(*,*) 'v1', l, m, v1, alm_s, alm_b, mfac
!!$        end if
!!$        if (isNaN(abs(v2))) then
!!$           write(*,*) 'v2', i, l, m, v2
!!$        end if
!     end if


        ! Positive spin
        almc = spinsign*self%lnorm(l) * (v1+conjg(v2)*mfac)
!!$        if (m == 64) then
!!$           write(*,*) i, l, m, m_b, almc
!!$           write(*,*) spinsign, self%lnorm(l)
!!$           write(*,*) alm_b
!!$           write(*,*) alm_s
!!$        end if
        if (m == 0) then
           alm(i,1) = real(almc)
        else
           alm(i,1)    = real(almc) * sqrt_two
           alm(ineg,1) = imag(almc) * sqrt_two
        end if

        if (m_b > 0) then
           ! Negative spin
           almc = -spinsign*self%lnorm(l) * (v1-conjg(v2)*mfac)
           !almc = -cmplx(0,1)*spinsign*self%lnorm(l) * (v1-conjg(v2)*mfac)
           if (m == 0) then
              alm(i,2) = real(almc)
           else
              alm(i,2)    = real(almc) * sqrt_two
              alm(ineg,2) = imag(almc) * sqrt_two
           end if
        end if

!!$        call map%info%lm2i(l, m, j)        
!!$        write(*,*) l, m, alm(i,1)/(map%alm(j,1)*exp(-0.5d0*l*(l+1)*(420.d0*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2))
!!$        !alm(i,1) = map%alm(j,1)*exp(-0.5d0*l*(l+1)*(420.d0*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
!!$
!!$        call map%info%lm2i(l, -m, j)        
!!$        write(*,*) l, -m, alm(ineg,1)/(map%alm(j,1)*exp(-0.5d0*l*(l+1)*(420.d0*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)) 
!!$        !alm(ineg,1) = map%alm(j,1)*exp(-0.5d0*l*(l+1)*(420.d0*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)


     end do

  end subroutine get_alms


  !deallocate memory used by the class
  subroutine dealloc(self)
    implicit none
    class(comm_conviqt), intent(inout) :: self

    deallocate(self%c)
    deallocate(self%lnorm)
    call dealloc_shared_2d_spc(self%alm_beam)

  end subroutine dealloc

end module comm_conviqt_mod
