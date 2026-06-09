!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_tod_cgmap_mod
  use comm_map_mod
  use comm_tod_mod
  use comm_tod_Tbol_mod
  use comm_tod_crosstalk_mod
  use comm_fft_mod
  implicit none

  private
  public comm_cgmap, dealloc_cgmap
  
  type comm_cgmap
     type(MPI_Comm) :: comm
     integer(i4b)   :: myid, nprocs
     integer(i4b)   :: maptype, ndet, ncol, nside, npix, nscan, ntod, nobs
     type(C_PTR)    :: plan_fwd, plan_back
     class(comm_tod), pointer :: tod0
     logical(lgt),            allocatable, dimension(:,:)   :: accept      ! Accept status, (ndet,nscan)
     integer(i4b),            allocatable, dimension(:)     :: ind_scan    ! Index range for each scan (0:nscan)
     real(dp),                allocatable, dimension(:,:)   :: x           ! Final full-sky map; only stored by root processor
     real(dp),                allocatable, dimension(:,:,:) :: invM        ! Diagonal preconditioner, (ncol,ncol,npix)
     real(dp),                allocatable, dimension(:,:)   :: xi          ! Pixel index based mini-map; separate for each core
     real(sp),                allocatable, dimension(:,:)   :: tod         ! Stacked, calibrated and in-painted TOD, (nscan*ntod,ntod)
     real(dp),                allocatable, dimension(:,:)   :: invN        ! 1/sigma**2 per scan, (ndet,nscan)
     integer(i4b),            allocatable, dimension(:,:)   :: ind         ! Stacked pixel index numbers, (nscan*ntod,ndet)
     integer(i4b),            allocatable, dimension(:,:)   :: psi         ! Stacked psi index numbers, (nscan*ntod,ndet)
     integer(i4b),            allocatable, dimension(:,:)   :: mask        ! Stacked TOD processing mask, (nscan*ntod,ndet)
     type(comm_crosstalk),                 pointer          :: W_S         ! Signal crosstalk matrix
     type(comm_crosstalk),                 pointer          :: W_N         ! Noise crosstalk matrix
     type(Tbol_ptr),          allocatable, dimension(:)     :: T           ! Bolometer transfer function, (ndet)
     integer(i4b),            allocatable, dimension(:)     :: ind2pix     ! Conversion table from reduced to full pixel index; individual for each processor
   contains
     procedure :: load_data => load_calibrated_data
     procedure :: solve     => solve_cgmap
     procedure :: A         => multiply_Amat
     procedure :: compute_rhs
     procedure :: apply_precond
     procedure :: init_precond
     procedure :: x_bcast
     procedure :: x_reduce
     procedure :: P         => apply_P
     procedure :: Pt        => apply_Pt
     procedure :: multiply_invN
     procedure :: convolve_T
  end type comm_cgmap

  interface comm_cgmap
     procedure constructor_cgmap
  end interface comm_cgmap
    
contains

  function constructor_cgmap(tod0, maptype, W_S, W_N) result(c)
    implicit none
    class(comm_tod),      target,       intent(in) :: tod0
    integer(i4b),                       intent(in) :: maptype
    type(comm_crosstalk), target,       intent(in), optional :: W_S
    type(comm_crosstalk), target,       intent(in), optional :: W_N
    class(comm_cgmap), pointer                  :: c

    integer(i4b) :: i, ierr, nfft, n
    real(sp),     allocatable, dimension(:)   :: dt
    complex(spc), allocatable, dimension(:)   :: dv    
    
    allocate(c)
    c%comm      = tod0%comm
    call mpi_comm_rank(c%comm, c%myid, ierr)
    call mpi_comm_size(c%comm, c%nprocs, ierr) 

    c%tod0      => tod0
    c%maptype   = maptype
    c%nside     = tod0%nside
    c%npix      = 12*c%nside**2
    c%ndet      = tod0%ndet
    c%nscan     = tod0%nscan
    c%ntod      = sum(tod0%scans%ntod)
    c%nobs      = size(tod0%pixcache%ind2pix)

    if (c%maptype == 1) then
       ! Temperature only
       c%ncol = 1
    else if (c%maptype == 2) then
       ! TQU
       c%ncol = 3
    else
       write(*,*) 'CGmap type not supported = ', maptype
       stop
    end if
    
    allocate(c%ind2pix(c%nobs))
    c%ind2pix = tod0%pixcache%ind2pix
    
    ! Find start index for each scan; each scan is stored in (ind_scan(i-1)+1):ind_scan(i)
    allocate(c%ind_scan(0:c%nscan))
    c%ind_scan(0) = 0
    do i = 1, c%nscan
       c%ind_scan(i) = c%ind_scan(i-1) + tod0%scans(i)%ntod
    end do
    
    allocate(c%tod(c%ntod,c%ndet))
    allocate(c%ind(c%ntod,c%ndet))
    allocate(c%mask(c%ntod,c%ndet))
    allocate(c%invN(c%ndet,c%nscan))
    allocate(c%xi(c%ncol,c%nobs))
    allocate(c%accept(c%ndet,c%nscan))
    if (c%myid == 0) allocate(c%x(c%ncol,0:c%npix-1))
    if (c%myid == 0) allocate(c%invM(c%ncol,c%ncol,0:c%npix-1))
    c%accept = .false.
    c%tod    = 0.

    if (c%maptype > 1) then ! Polarization
       allocate(c%psi(c%ntod,c%ndet))
    end if
    
    ! Check for optional matrix operators
    if (present(W_S))  c%W_S  => W_S
    if (present(W_N))  c%W_N  => W_N
    if (tod0%correct_Tbol) then
       allocate(c%T(c%ndet))
       do i = 1, c%ndet
          c%T(i)%p => tod0%Tbol(i)%p
       end do

       nfft     = get_closest_fft_pow2(c%ntod)
       n        = nfft / 2 + 1
       allocate(dt(nfft), dv(0:n-1))
       c%plan_fwd  = fftwf_plan_dft_r2c_1d(nfft, dt, dv, fftw_estimate + fftw_unaligned)
       c%plan_back = fftwf_plan_dft_c2r_1d(nfft, dv, dt, fftw_estimate + fftw_unaligned)
       deallocate(dt, dv)
    end if
    
  end function constructor_cgmap

  subroutine dealloc_cgmap(c)
    implicit none
    class(comm_cgmap), pointer, intent(inout) :: c
    if (allocated(c%ind_scan))  deallocate(c%ind_scan)
    if (allocated(c%tod))       deallocate(c%tod)
    if (allocated(c%ind))       deallocate(c%ind)
    if (allocated(c%psi))       deallocate(c%psi)
    if (allocated(c%mask))      deallocate(c%mask)
    if (allocated(c%invM))      deallocate(c%invM)
    if (allocated(c%x))         deallocate(c%x)
    if (allocated(c%xi))        deallocate(c%xi)
    if (allocated(c%ind2pix))   deallocate(c%ind2pix)
    if (allocated(c%accept))    deallocate(c%accept)
    call fftw_destroy_plan(c%plan_fwd)
    call fftw_destroy_plan(c%plan_back)
    deallocate(c)
  end subroutine dealloc_cgmap

  subroutine load_calibrated_data(self, scan, tod, sd, d_calib)
    implicit none
    class(comm_cgmap),                   intent(inout)  :: self
    integer(i4b),                        intent(in)     :: scan
    class(comm_tod),                     intent(in)     :: tod
    class(comm_scandata),                intent(in)     :: sd
    real(sp),          dimension(1:,1:), intent(in)     :: d_calib !(ntod, ndet)

    integer(i4b) :: i, j, det
    
    i = self%ind_scan(scan-1)+1
    j = self%ind_scan(scan)
    self%tod(i:j,:)     = d_calib
    self%ind(i:j,:)     = sd%ind(:,:,1)
    self%mask(i:j,:)    = iand(sd%flag,tod%flag0)
    do det = 1, self%ndet
       self%accept(det,scan) = tod%scans(scan)%d(det)%accept
       self%invN(det,scan)   = 1./(tod%scans(scan)%d(det)%N_psd%sigma0/&
            & tod%scans(scan)%d(det)%gain)**2
    end do

    if (self%maptype > 1) then
       self%psi(i:j,:) = sd%psi(:,:,1)
    end if

    ! Inpaint flagged samples
    do det = 1, self%ndet
       if (self%accept(det,scan)) then
          where (self%mask(i:j,det) /= 0.)
             self%tod(i:j,det) = tod%scans(scan)%d(det)%gain * sd%s_sky(:,det,0,1)
          end where
       else
          self%tod(i:j,det) = tod%scans(scan)%d(det)%gain * sd%s_sky(:,det,0,1)
       end if
    end do
    
  end subroutine load_calibrated_data

  subroutine solve_cgmap(self, map_out)
    implicit none
    class(comm_cgmap), intent(inout) :: self
    class(comm_map),   intent(inout) :: map_out

    integer(i4b) :: i, j, oper, MAXITER=3000, ierr
    real(dp)     :: delta_old, delta_new, delta0, eps = 1d-6, lim_convergence, dq, alpha, beta, t1, t2
    real(dp), allocatable, dimension(:,:) :: b, invMb, r, d, q, s

    ! Initialize preconditioner
    call self%init_precond
    
    if (self%myid == 0) then


       if (sum(abs(self%invM)) == 0.d0) then
          write(*,*) "CGmap error: Preconditioner is zero; no accepted data; . Exiting."
          stop
       end if
       
       ! Allocate temporary arrays
       allocate(b(self%ncol,self%npix))
       allocate(invMb(self%ncol,self%npix))
       allocate(r(self%ncol,self%npix))
       allocate(d(self%ncol,self%npix))
       allocate(q(self%ncol,self%npix))
       allocate(s(self%ncol,self%npix))       
       
       ! Initialize RHS
       call self%compute_RHS(b)
       
       ! Initialize map
       self%x  = 0.d0
       
       ! Initialize CG variables
       r  = b 
       call self%apply_precond(r, d)
       
       delta_new = sum(r*d)
       call self%apply_precond(b, invMb)
       delta0    = sum(b*invMb)
       lim_convergence = eps*delta0
       
       ! Run CG search
       do i = 1, MAXITER
          call wall_time(t1)
          
          ! Check convergence
          if (delta_new < lim_convergence) exit
          
          call self%A(d, q)
          dq        = sum(d*q)
          alpha     = delta_new / dq
          self%x    = self%x + alpha * d
          r         = r      - alpha * q
          call self%apply_precond(r, s)
          delta_old = delta_new 
          delta_new = sum(r*s)
          beta      = delta_new / delta_old
          d         = s + beta * d

          call wall_time(t2)
          write(*,fmt='(a,i5,a,e13.5,a,e13.5,a,f8.2)') ' |  CG iter. ', i, ' -- res = ', &
               & min(delta_new,1d30), ', tol = ', real(lim_convergence,sp), &
               & ', time = ', real(t2-t1,sp)
       end do
       oper = 0
       call mpi_bcast(oper, 1, MPI_INTEGER, 0, self%comm, ierr)  ! Release slaves
       
       call wall_time(t2)
       write(*,fmt='(a,i5,a,e13.5,a,e13.5,a,f8.2)') ' |  Final CG iter ', i, ' -- res = ', &
            & real(delta_new,sp), ', tol = ', real(lim_convergence,sp)

       deallocate(b, invMb, r, d, q, s)
    else

       ! Slave processor
       oper = 1
       do while (oper > 0)
          call mpi_bcast(oper, 1, MPI_INTEGER, 0, self%comm, ierr)
          if (oper == 1) then
             call self%A
          else if (oper == 2) then
             call self%compute_rhs
          end if
       end do

    end if

    ! Distribute full-sky map
    call map_out%bcast_fullsky_from_root(transpose(self%x))
    
  end subroutine solve_cgmap

  subroutine multiply_Amat(self, x, Ax)
    implicit none
    class(comm_cgmap),                 intent(inout)         :: self
    real(dp),          dimension(:,:), intent(in),  optional :: x
    real(dp),          dimension(:,:), intent(out), optional :: Ax

    integer(i4b) :: ierr, oper

    ! Activate slaves
    oper = 1
    if (self%myid /= 0) oper = 0
    call mpi_bcast(oper, 1, MPI_INTEGER, 0, self%comm, ierr)

    ! Distrbute x; stored in self%xi on each core
    call self%x_bcast(x)

    ! Compute y = P * x, overwrite self%tod
    call self%P

    ! Compute y = T P * x
    if (allocated(self%T)) call self%convolve_T(.false.)

    ! Compute y = W_S T P * x
    if (associated(self%W_S)) call self%W_S%multiply(.false., self%tod)
    
    ! Compute y = W_N^t * invN * W_N * d
    if (associated(self%W_N)) call self%W_N%multiply(.false., self%tod)
    call self%multiply_invN
    if (associated(self%W_N)) call self%W_N%multiply(.true., self%tod)

    ! Compute y = W_S^t invN d
    if (associated(self%W_S)) call self%W_S%multiply(.true., self%tod)

        ! Compute y = T^t W_S^t * invN d
    if (allocated(self%T)) call self%convolve_T(.true.)

    ! Compute y = P^t invN d
    call self%Pt

    ! Collect results from all cores; stored in self%x on root
    call self%x_reduce(Ax)
    
  end subroutine multiply_Amat

  ! Compute rhs = P^t T^t W_S^t invN d
  subroutine compute_RHS(self, b)
    implicit none
    class(comm_cgmap),                 intent(inout)           :: self
    real(dp),          dimension(:,:), intent(inout), optional :: b
    
    integer(i4b) ::  ierr, oper

    ! Activate slaves
    oper = 2
    if (self%myid /= 0) oper = 0
    !call mpi_bcast(oper, 1, MPI_INTEGER, 0, self%comm, ierr)
    call mpi_bcast(oper, 1, MPI_INTEGER, 0, self%comm, ierr)

    ! Compute y = W_N^t * invN * W_N * d
    if (associated(self%W_N)) call self%W_N%multiply(.false., self%tod)
    call self%multiply_invN
    if (associated(self%W_N)) call self%W_N%multiply(.true., self%tod)

    ! Compute y = W_S^t invN d
    if (associated(self%W_S)) call self%W_S%multiply(.true., self%tod)

    ! Compute y = T^t W_S^t * invN d
    if (allocated(self%T)) call self%convolve_T(.true.)

    ! Compute y = P^t invN d
    call self%Pt

    ! Collect results from all cores
    call self%x_reduce(b)

  end subroutine compute_RHS
 
  subroutine apply_precond(self, x, invMx)
    implicit none
    class(comm_cgmap),                   intent(in)  :: self
    real(dp),          dimension(1:,0:), intent(in)  :: x
    real(dp),          dimension(1:,0:), intent(out) :: invMx

    integer(i4b) :: i
    do i = 0, self%npix-1
       invMx(:,i) = matmul(self%invM(:,:,i), x(:,i))
    end do
  end subroutine apply_precond

  subroutine init_precond(self)
    implicit none
    class(comm_cgmap),                 intent(inout)    :: self

    integer(i4b) :: i, j, col, scan, det, t
    real(dp)     :: w, eig(3)

    do col = 1, self%ncol
       ! Build the preconditioner column-by-column, in order to reuse the MPI infrastructure    
       self%xi = 0.d0
       do scan = 1, self%nscan
          do det = 1, self%ndet
             if (.not. self%accept(det,scan)) cycle
             do t = self%ind_scan(scan-1)+1, self%ind_scan(scan)
                if (self%mask(t,det) /= 0) cycle
                if (self%maptype == 1) then
                   self%xi(1,self%ind(t,det)) = self%xi(1,self%ind(t,det)) + self%invN(det,scan)
                else if (self%maptype == 2) then
                   if (col == 1) then
                      w = 1.
                   else if (col == 2) then
                      w = self%tod0%pixcache%cos2psi(self%psi(t,det))
                   else if (col == 3) then
                      w = self%tod0%pixcache%sin2psi(self%psi(t,det))
                   end if
                   self%xi(1,self%ind(t,det)) = self%xi(1,self%ind(t,det)) + self%invN(det,scan) * w 
                   self%xi(2,self%ind(t,det)) = self%xi(2,self%ind(t,det)) + self%invN(det,scan) * w * self%tod0%pixcache%cos2psi(self%psi(t,det))
                   self%xi(3,self%ind(t,det)) = self%xi(3,self%ind(t,det)) + self%invN(det,scan) * w * self%tod0%pixcache%sin2psi(self%psi(t,det))
                end if
             end do
          end do
       end do

       call self%x_reduce(self%invM(:,col,:))
    end do
    
    if (self%myid == 0) then
       do i = 0, self%npix-1
          if (self%invM(1,1,i) > 0.) then
             call invert_singular_matrix(self%invM(:,:,i), 1.d-8)
          else
             self%invM(:,:,i) = 0.d0
          end if
       end do
    end if
    
  end subroutine init_precond
  
  subroutine multiply_invN(self)
    implicit none
    class(comm_cgmap),                 intent(inout)    :: self

    integer(i4b) :: i, j, scan, det
    
    do scan = 1, self%nscan
       i = self%ind_scan(scan-1)+1
       j = self%ind_scan(scan)
       do det = 1, self%ndet
          if (.not. self%accept(det,scan)) cycle
          where (self%mask(i:j,det) == 0)
             self%tod(i:j,det) = self%tod(i:j,det) * self%invN(det,scan)
          elsewhere
             self%tod(i:j,det) = 0.  ! Remove flagged samples
          end where
       end do
    end do
    
  end subroutine multiply_invN

  subroutine apply_P(self)
    implicit none
    class(comm_cgmap),                 intent(inout)    :: self

    integer(i4b) :: t, det, scan
    
    do scan = 1, self%nscan
       do det = 1, self%ndet
          if (.not. self%accept(det,scan)) cycle
          do t = self%ind_scan(scan-1)+1, self%ind_scan(scan)
             if (self%mask(t,det) /= 0) cycle
             if (self%maptype == 1) then
                self%tod(t,det) = self%xi(1,self%ind(t,det))
             else if (self%maptype == 2) then
                self%tod(t,det) = self%xi(1,self%ind(t,det)) +                                                 & ! T
                               & self%xi(2,self%ind(t,det)) * self%tod0%pixcache%cos2psi(self%psi(t,det)) + & ! Q
                               & self%xi(3,self%ind(t,det)) * self%tod0%pixcache%sin2psi(self%psi(t,det))     ! U
             end if
          end do
       end do
    end do
    
  end subroutine apply_P

  subroutine apply_Pt(self)
    implicit none
    class(comm_cgmap),                 intent(inout) :: self

    integer(i4b) :: t, det, scan

    self%xi = 0.d0
    do scan = 1, self%nscan
       do det = 1, self%ndet
          if (.not. self%accept(det,scan)) cycle
          do t = self%ind_scan(scan-1)+1, self%ind_scan(scan)
             if (self%mask(t,det) /= 0) cycle
             if (self%maptype == 1) then
                self%xi(1,self%ind(t,det)) = self%xi(1,self%ind(t,det)) + self%tod(t,det)
             else if (self%maptype == 2) then
                self%xi(1,self%ind(t,det)) = self%xi(1,self%ind(t,det)) + self%tod(t,det)                                                 ! T
                self%xi(2,self%ind(t,det)) = self%xi(2,self%ind(t,det)) + self%tod(t,det) * self%tod0%pixcache%cos2psi(self%psi(t,det)) ! Q
                self%xi(3,self%ind(t,det)) = self%xi(3,self%ind(t,det)) + self%tod(t,det) * self%tod0%pixcache%sin2psi(self%psi(t,det)) ! U
             end if
          end do
       end do
    end do
    
  end subroutine apply_Pt

  subroutine convolve_T(self, trans)
    implicit none
    class(comm_cgmap),                 intent(inout) :: self
    logical(lgt),                      intent(in)    :: trans

    integer(i4b) :: det

    do det = 1, self%ndet
       if (trans) then
          call self%T(det)%p%convolve("cgmap", "transpose", self%tod(:,det), self%plan_fwd, self%plan_back)
       else
          call self%T(det)%p%convolve("cgmap", "regular",   self%tod(:,det), self%plan_fwd, self%plan_back)
       end if
    end do
    
  end subroutine convolve_T

  subroutine x_bcast(self, x)
    implicit none
    class(comm_cgmap),                   intent(inout)          :: self
    real(dp),          dimension(1:,0:), intent(in),   optional :: x

    integer(i4b) :: i, j, scan, det, nobs, ierr
    integer(i4b), allocatable, dimension(:)   :: p
    real(dp),     allocatable, dimension(:,:) :: buffer
    type(MPI_Status)  :: mpistat
    
    if (self%myid == 0) then
       self%xi = x(:,self%ind2pix)
       do i = 1, self%nprocs-1
          call mpi_recv(nobs,       1, MPI_INTEGER, i, 98, self%comm, mpistat, ierr)
          allocate(p(nobs), buffer(self%ncol,nobs))
          call mpi_recv(p,      nobs, MPI_INTEGER, i, 98, self%comm, mpistat, ierr)
          buffer = x(:,p)
          call mpi_send(buffer,      size(buffer), MPI_DOUBLE_PRECISION, i, 98, &
               & self%comm, ierr)
          deallocate(p, buffer)
       end do
    else
       call mpi_send(self%nobs,    1,             MPI_INTEGER, 0, 98, self%comm, ierr)
       call mpi_send(self%ind2pix, self%nobs,     MPI_INTEGER, 0, 98, self%comm, ierr)
       call mpi_recv(self%xi,      size(self%xi), &
            & MPI_DOUBLE_PRECISION, 0, 98, self%comm, mpistat, ierr)
    end if

  end subroutine x_bcast

  subroutine x_reduce(self, x)
    implicit none
    class(comm_cgmap),                   intent(in)            :: self
    real(dp),          dimension(1:,0:), intent(out), optional :: x

    integer(i4b) :: i, j, scan, det, nobs, ierr
    integer(i4b), allocatable, dimension(:)   :: p
    real(dp),     allocatable, dimension(:,:) :: buffer
    type(MPI_Status)  :: mpistat
    
    if (self%myid == 0) then
       x                 = 0.d0
       x(:,self%ind2pix) = self%xi 
       do i = 1, self%nprocs-1
          call mpi_recv(nobs,       1, MPI_INTEGER, i, 98, self%comm, mpistat, ierr)
          allocate(p(nobs), buffer(self%ncol,nobs))
          call mpi_recv(p,      nobs, MPI_INTEGER, i, 98, self%comm, mpistat, ierr)
          call mpi_recv(buffer, size(buffer), &
               & MPI_DOUBLE_PRECISION, i, 98, self%comm, mpistat, ierr)
          x(:,p) = x(:,p) + buffer
          deallocate(p, buffer)
       end do
    else
       call mpi_send(self%nobs,    1,             MPI_INTEGER, 0, 98, self%comm, ierr)
       call mpi_send(self%ind2pix, self%nobs,     MPI_INTEGER, 0, 98, self%comm, ierr)
       call mpi_send(self%xi,      size(self%xi), MPI_DOUBLE_PRECISION, 0, 98, &
            & self%comm, ierr)
    end if
    
  end subroutine x_reduce

  
end module comm_tod_cgmap_mod
