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
  implicit none

  private
  public comm_cgmap, dealloc_cgmap
  
  type comm_cgmap
     integer(i4b)  :: comm, myid, nprocs
     integer(i4b)  :: ndet, ncol, nside, npix, nscan, ntod, nobs
     logical(lgt),            allocatable, dimension(:,:) :: accept      ! Accept status, (ndet,nscan)
     integer(i4b),            allocatable, dimension(:)   :: col_def     ! Column Stokes definition; T_full=1, Q_full=2, U_full=3, S_det=det+3, T_det=-det
     integer(i4b),            allocatable, dimension(:)   :: ind_scan    ! Index range for each scan (0:nscan)
     real(dp),                allocatable, dimension(:,:) :: x           ! Final full-sky map; only stored by root processor
     real(dp),                allocatable, dimension(:,:) :: invM        ! Diagonal preconditioner
     real(dp),                allocatable, dimension(:,:) :: xi          ! Pixel index based mini-map; separate for each core
     real(sp),                allocatable, dimension(:,:) :: tod         ! Stacked, calibrated and in-painted TOD, (ndet,nscan*ntod)
     real(sp),                allocatable, dimension(:,:) :: invN        ! 1/sigma**2 per scan, (ndet,nscan)
     integer(i4b),            allocatable, dimension(:,:) :: ind         ! Stacked pixel index numbers, (ndet,nscan*ntod)
     integer(i4b),            allocatable, dimension(:,:) :: mask        ! Stacked TOD processing mask, (ndet,nscan*ntod)
     type(comm_crosstalk),                 pointer        :: W_S         ! Signal crosstalk matrix
     type(comm_crosstalk),                 pointer        :: W_N         ! Noise crosstalk matrix
     type(Tbol_ptr),          allocatable, dimension(:)   :: T           ! Bolometer transfer function, (ndet)
     integer(i4b),            allocatable, dimension(:)   :: ind2pix     ! Conversion table from reduced to full pixel index; individual for each processor
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

  function constructor_cgmap(comm, nside, ndet, ntod_scan, col_def, ind2pix, W_S, W_N, Tbol) result(c)
    implicit none
    integer(i4b),                       intent(in) :: comm, nside, ndet
    integer(i4b),         dimension(:), intent(in) :: ntod_scan
    integer(i4b),         dimension(:), intent(in) :: col_def
    integer(i4b),         dimension(:), intent(in) :: ind2pix
    type(comm_crosstalk), target,       intent(in), optional :: W_S
    type(comm_crosstalk), target,       intent(in), optional :: W_N
    type(Tbol_ptr),       dimension(:), intent(in), optional :: Tbol
    class(comm_cgmap), pointer                  :: c

    integer(i4b) :: i, ierr
    
    allocate(c)
    c%comm      = comm
    call mpi_comm_rank(comm, c%myid, ierr)
    call mpi_comm_size(comm, c%nprocs, ierr) 
    
    c%nside     = nside
    c%npix      = 12*c%nside**2
    c%ndet      = ndet
    c%nscan     = size(ntod_scan)
    c%ncol      = size(col_def)
    c%ntod      = sum(ntod_scan)
    c%nobs      = size(ind2pix)
    
    allocate(c%col_def(c%ncol))
    allocate(c%ind2pix(c%nobs))
    c%col_def   = col_def
    c%ind2pix   = ind2pix
    
    ! Find start index for each scan; each scan is stored in (ind_scan(i-1)+1):ind_scan(i)
    allocate(c%ind_scan(0:c%nscan))
    c%ind_scan(0) = 0
    do i = 1, c%nscan
       c%ind_scan(i) = c%ind_scan(i-1) + ntod_scan(i)
    end do
    
    allocate(c%tod(c%ndet,c%ntod))
    allocate(c%ind(c%ndet,c%ntod))
    allocate(c%mask(c%ndet,c%ntod))
    allocate(c%invN(c%ndet,c%nscan))
    allocate(c%xi(c%ncol,c%nobs))
    allocate(c%accept(c%ndet,c%nscan))
    if (c%myid == 0) allocate(c%x(c%ncol,c%npix))
    if (c%myid == 0) allocate(c%invM(c%ncol,c%npix))
    c%accept = .false.
    
    ! Check for optional matrix operators
    if (present(W_S))  c%W_S  => W_S
    if (present(W_N))  c%W_N  => W_N
    if (present(Tbol)) then
       allocate(c%T(c%ndet))
       do i = 1, c%ndet
          c%T(i)%p => Tbol(i)%p
       end do
    end if
    
  end function constructor_cgmap

  subroutine dealloc_cgmap(c)
    implicit none
    class(comm_cgmap), pointer, intent(inout) :: c
    if (allocated(c%col_def))   deallocate(c%col_def)
    if (allocated(c%ind_scan))  deallocate(c%ind_scan)
    if (allocated(c%tod))       deallocate(c%tod)
    if (allocated(c%ind))       deallocate(c%ind)
    if (allocated(c%mask))      deallocate(c%mask)
    if (allocated(c%invM))      deallocate(c%invM)
    if (allocated(c%x))         deallocate(c%x)
    if (allocated(c%xi))        deallocate(c%xi)
    if (allocated(c%ind2pix))   deallocate(c%ind2pix)
    if (allocated(c%accept))   deallocate(c%accept)
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
    self%tod(:,i:j)     = transpose(d_calib)
    self%ind(:,i:j)     = transpose(sd%ind(:,:,1))
    self%mask(:,i:j)    = transpose(iand(sd%flag,tod%flag0))
    do det = 1, self%ndet
       self%accept(det,scan) = tod%scans(scan)%d(det)%accept
       self%invN(det,scan)   = 1./(tod%scans(scan)%d(det)%N_psd%sigma0/&
            & tod%scans(scan)%d(det)%gain)**2
    end do
    
  end subroutine load_calibrated_data

  subroutine solve_cgmap(self, map_out)
    implicit none
    class(comm_cgmap), intent(inout) :: self
    class(comm_map),   intent(inout) :: map_out

    integer(i4b) :: i, j, oper, MAXITER=30, ierr
    real(dp)     :: delta_old, delta_new, delta0, eps = 1d-6, lim_convergence, dq, alpha, beta, t1, t2
    real(dp), allocatable, dimension(:,:) :: b, invMb, r, d, q, s

    ! Initialize preconditioner
    call self%init_precond
    
    if (self%myid == 0) then

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
       call mpi_bcast(0, 1, MPI_INTEGER, 0, self%comm, ierr)  ! Release slaves
       
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

    integer(i4b) :: ierr

    ! Activate slaves
    if (self%myid == 0) call mpi_bcast(1, 1, MPI_INTEGER, 0, self%comm, ierr)

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
    
    integer(i4b) ::  ierr

    ! Activate slaves
    if (self%myid == 0) call mpi_bcast(2, 1, MPI_INTEGER, 0, self%comm, ierr)

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
    class(comm_cgmap),                 intent(in)  :: self
    real(dp),          dimension(:,:), intent(in)  :: x
    real(dp),          dimension(:,:), intent(out) :: invMx
    invMx = self%invM * x
  end subroutine apply_precond

  subroutine init_precond(self)
    implicit none
    class(comm_cgmap),                 intent(inout)    :: self

    integer(i4b) :: i, j, scan, det, t

    self%xi = 0.d0
    do scan = 1, self%nscan
       do det = 1, self%ndet
          if (.not. self%accept(det,scan)) cycle
          do t = self%ind_scan(scan-1)+1, self%ind_scan(scan)
             if (self%mask(det,t) == 0) then
                self%xi(1,self%ind(det,t)) = self%xi(1,self%ind(det,t)) + self%invN(det,scan)
             end if
          end do
       end do
    end do

    call self%x_reduce(self%invM)
    if (self%myid == 0) then
       where (self%invM > 0.d0) 
          self%invM = 1.d0 / self%invM
       end where
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
          where (self%mask(det,i:j) == 0)
             self%tod(det,i:j) = self%tod(det,i:j) * self%invN(det,scan)
          elsewhere
             self%tod(det,i:j) = 0.  ! Remove flagged samples
          end where
       end do
    end do
    
  end subroutine multiply_invN

  subroutine apply_P(self)
    implicit none
    class(comm_cgmap),                 intent(inout)    :: self

    integer(i4b) :: i, j, det, scan
    
    do scan = 1, self%nscan
       i = self%ind_scan(scan-1)+1
       j = self%ind_scan(scan)
       do det = 1, self%ndet
          if (.not. self%accept(det,scan)) cycle
          self%tod(det,i:j) = self%xi(1,self%ind(det,i:j))
       end do
    end do
    
  end subroutine apply_P

  subroutine apply_Pt(self)
    implicit none
    class(comm_cgmap),                 intent(inout) :: self

    integer(i4b) :: i, j, det, scan

    self%xi = 0.d0
    do scan = 1, self%nscan
       i = self%ind_scan(scan-1)+1
       j = self%ind_scan(scan)
       do det = 1, self%ndet
          if (.not. self%accept(det,scan)) cycle
          self%xi(1,self%ind(det,i:j)) = self%xi(1,self%ind(det,i:j)) + self%tod(det,i:j)
       end do
    end do
    
  end subroutine apply_Pt

  subroutine convolve_T(self, trans)
    implicit none
    class(comm_cgmap),                 intent(inout) :: self
    logical(lgt),                      intent(in)    :: trans

    integer(i4b) :: i, j, scan, det   

    do scan = 1, self%nscan
       i = self%ind_scan(scan-1)+1
       j = self%ind_scan(scan)
       do det = 1, self%ndet
          if (.not. self%accept(det,scan)) cycle
          call self%T(det)%p%convolve(self%tod(det,i:j))
       end do
    end do

  end subroutine convolve_T

  subroutine x_bcast(self, x)
    implicit none
    class(comm_cgmap),                   intent(inout)          :: self
    real(dp),          dimension(1:,0:), intent(in),   optional :: x

    integer(i4b) :: i, j, scan, det, nobs, ierr
    integer(i4b), allocatable, dimension(:)   :: p
    real(dp),     allocatable, dimension(:,:) :: buffer
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat
    
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
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: mpistat
    
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
