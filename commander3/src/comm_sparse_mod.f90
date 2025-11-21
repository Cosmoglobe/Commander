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
module comm_sparse_mod
  use healpix_types
  use math_tools
  implicit none

  private
  public sparse_system, sparse_ptr, sparse_system_tester

  type sparse_system
     integer(i8b) :: pt(64)
     integer(i4b) :: n, ni, nj, n_a_max, icurr, iparm(64)
     integer(i4b) :: maxfct, mnum, mtype, nrhs, msglvl, idum
     real(dp)     :: ddum
     integer(i4b), allocatable, dimension(:) :: ia, ja
     real(dp),     allocatable, dimension(:) :: a, b, x, y
   contains
     procedure :: decomp    => sparse_system_decomp
     procedure :: solve     => sparse_system_solve
     procedure :: dealloc   => sparse_system_dealloc
     procedure :: set_A     => sparse_system_add_matrix_element
     procedure :: set_rhs   => sparse_system_set_rhs
     procedure :: add_diag  => sparse_system_add_diag
     procedure :: scale_row => sparse_system_scale_row
     procedure :: scale_col => sparse_system_scale_col
     procedure :: print     => sparse_system_print
  end type sparse_system

  interface sparse_system
     procedure constructor_sparse_system
  end interface sparse_system

  type sparse_ptr
     class(sparse_system), pointer :: p => null()
  end type sparse_ptr
  
contains

  function constructor_sparse_system(n, n_A_max) result(c)
    implicit none
    integer(i4b),         intent(in) :: n, n_A_max
    class(sparse_system), pointer    :: c

    ! Assign memory
    allocate(c)

    ! General parameters
    c%n       = n
    c%n_a_max = n_a_max
    c%nrhs    = 1
    c%mtype   = 2        ! symmetric matrix, positive definite
    c%msglvl  = 1         ! output statistical information
    c%maxfct  = 1
    c%mnum    = 1

    ! Initialize arrays
    c%ni    = 0
    c%nj    = 0
    c%icurr = 0
    allocate(c%ia(c%n_a_max), c%ja(c%n_a_max), c%a(c%n_a_max), c%b(n), c%x(n), c%y(n))
    c%ia    = 0
    c%ja    = 0
    c%a     = 0.d0
    c%b     = 0.d0

  end function constructor_sparse_system

  subroutine sparse_system_add_matrix_element(self, i, j, A_ij)
    implicit none
    class(sparse_system),     intent(inout) :: self
    integer(i4b),             intent(in) :: i, j
    real(dp),                 intent(in) :: A_ij     

    integer(i4b), allocatable, dimension(:) :: ibuff
    real(dp),     allocatable, dimension(:) :: fbuff

    ! Check if array is already full; if so, double its size
    if (self%nj == self%n_a_max) then
       self%n_a_max = 2*self%n_a_max
       allocate(ibuff(self%nj), fbuff(self%nj))
       ! Extend row array
       ibuff(1:self%ni) = self%ia(1:self%ni)
       deallocate(self%ia); allocate(self%ia(self%n_a_max))
       self%ia(1:self%ni) = ibuff(1:self%ni)
       ! Extend column array
       ibuff(1:self%nj) = self%ja(1:self%nj)
       deallocate(self%ja); allocate(self%ja(self%n_a_max))
       self%ja(1:self%nj) = ibuff(1:self%nj)
       ! Extend coefficient array
       fbuff(1:self%nj)   = self%a(1:self%nj)
       deallocate(self%a); allocate(self%a(self%n_a_max))
       self%a(1:self%nj)  = fbuff(1:self%nj)
       self%a(self%nj+1:) = 0.d0
       deallocate(ibuff, fbuff)
    end if

    ! Insert new element
    self%nj          = self%nj+1
    self%ja(self%nj) = j
    self%a(self%nj)  = A_ij
    if (i /= self%icurr) then
       self%icurr         = i
       self%ni            = self%ni+1
       self%ia(self%ni)   = self%nj
       self%ia(self%ni+1) = self%nj+1 ! Sentinel value indicating end of array
    end if
    
  end subroutine sparse_system_add_matrix_element

  subroutine sparse_system_set_rhs(self, b)
    implicit none
    class(sparse_system),               intent(inout) :: self
    real(dp),             dimension(:), intent(in)    :: b
    if (size(b) /= self%n) then
       write(*,*) 'sparse_matrix_set_rhs: Incorrect size of input array = ', size(b), self%n
       stop
    else
       self%b = b
    end if
  end subroutine sparse_system_set_rhs
  
  subroutine sparse_system_decomp(self, verbosity)
    implicit none
    class(sparse_system),     intent(inout) :: self
    integer(i4b),             intent(in), optional :: verbosity

    integer(i4b) :: solver, error, phase, msglvl

    ! Initialize PARDISO solver
    solver     =  10  ! use sparse direct method
    call pardisoinit(self%pt, self%mtype, self%iparm)
    self%iparm(3) = 1    ! OMP_NUM_THREADS

    ! Symbolic factorization
    phase          = 11     ! only reordering and symbolic factorization
    self%iparm(33) = 1      ! compute determinant
    msglvl         = self%msglvl; if (present(verbosity)) msglvl = verbosity
    call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, phase, self%n, &
         & self%a(1:self%nj), self%ia(1:self%ni+1), self%ja(1:self%nj), &
         & self%idum, self%nrhs, self%iparm, msglvl, self%ddum, self%ddum, error)
    if (error /= 0) THEN
       write(*,*) 'pardiso symbolic factorization error: ', error
       stop
    else
       !write(*,*) 'pardiso -- number of nonzeros   = ', self%iparm(18)
       !write(*,*) 'pardiso -- MFLOPS               = ', self%iparm(19)
    end if

    ! Perform factorization
    phase     = 22
    call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, phase, self%n, &
         & self%a(1:self%nj), self%ia(1:self%ni+1), self%ja(1:self%nj), &
         & self%idum, self%nrhs, self%iparm, msglvl, self%ddum, self%ddum, error) 
    if (error /= 0) THEN
       write(*,*) 'pardiso factorization error: ', error
       stop
    end if
    
  end subroutine sparse_system_decomp
  
  subroutine sparse_system_solve(self, x)
    implicit none
    class(sparse_system),     intent(inout) :: self
    real(dp), dimension(:),   intent(out)   :: x

    integer(i4b) :: phase, error
    real(dp)     :: normb, normr
    
    ! Solve system
    self%iparm(8)  = 1   ! direct solution
    phase       = 33  ! only solve
    CALL pardiso (self%pt, self%maxfct, self%mnum, self%mtype, phase, self%n, self%a(1:self%nj), self%ia(1:self%ni+1), self%ja(1:self%nj), &
         & self%idum, self%nrhs, self%iparm, 0, self%b, self%x, error) 
    if (error /= 0) THEN
       write(*,*) 'pardiso solver error: ', error
       stop
    else
       x = self%x
    end if

  end subroutine sparse_system_solve

  subroutine sparse_system_dealloc(self)
    implicit none
    class(sparse_system),     intent(inout) :: self

    integer(i4b) :: phase, error

    
    ! Clean up pardiso
    phase     = -1   ! deallocate pardiso memory
    call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, phase, self%n, self%ddum, self%idum, self%idum, &
         & self%idum, self%nrhs, self%iparm, self%msglvl, self%ddum, self%ddum, error) 

    ! Deallocate object arrays
    deallocate(self%ia, self%ja, self%a, self%b, self%x, self%y)

    ! Nullify array sizes
    self%n  = 0
    self%ni = 0
    self%nj = 0
    
  end subroutine sparse_system_dealloc

  subroutine sparse_system_add_diag(self, val, row)
    implicit none
    class(sparse_system),               intent(inout) :: self
    real(dp),                           intent(in)    :: val
    integer(i4b),                       intent(in), optional :: row
    integer(i4b) :: i
    if (present(row)) then
       self%a(self%ia(row)) = self%a(self%ia(row)) + val
    else
       do i = 1, self%n
          self%a(self%ia(i)) = self%a(self%ia(i)) + val
       end do
    end if
  end subroutine sparse_system_add_diag

  subroutine sparse_system_scale_row(self, row, val)
    implicit none
    class(sparse_system),               intent(inout) :: self
    integer(i4b),                       intent(in)    :: row
    real(dp),                           intent(in)    :: val
    integer(i4b) :: c1, c2
    c1 = self%ia(row)
    c2 = self%ia(row+1)-1
    self%a(c1:c2) = self%a(c1:c2) * val
  end subroutine sparse_system_scale_row

  subroutine sparse_system_scale_col(self, col, val)
    implicit none
    class(sparse_system),               intent(inout) :: self
    integer(i4b),                       intent(in)    :: col
    real(dp),                           intent(in)    :: val
    where (self%ja == col)
       self%a = self%a * val
    end where
  end subroutine sparse_system_scale_col
  
  subroutine sparse_system_print(self)
    implicit none
    class(sparse_system),               intent(in) :: self

    write(*,*) 'n = ', self%n, ', ni = ', self%ni, ', nj = ', self%nj
    write(*,*) 'ia = ', self%ia(1:self%ni+1)
    write(*,*) 'ja = ', self%ja(1:self%nj)
    write(*,*) 'A  = ', self%a(1:self%nj)
    write(*,*) 'b  = ', self%b(1:self%n)
        
  end subroutine sparse_system_print

  subroutine sparse_system_tester
    implicit none

    integer(i4b) :: i, j, n
    real(dp), allocatable, dimension(:)   :: b, x
    real(dp), allocatable, dimension(:,:) :: A
    class(sparse_system), pointer :: ss

    ! Allocate system
    n = 8
    ss => sparse_system(n, 4)

    ! Set up coefficient matrix
    allocate(A(n,n)); A = 0.d0
    call ss%set_A(1,1,  7.d0); A(1,1) = 7.d0
    call ss%set_A(1,3,  1.d0); A(1,3) = 1.d0; A(3,1) = 1.d0
    call ss%set_A(1,6,  2.d0); A(1,6) = 2.d0; A(6,1) = 2.d0
    call ss%set_A(1,7,  7.d0); A(1,7) = 7.d0; A(7,1) = 7.d0
    call ss%set_A(2,2, -4.d0); A(2,2) = -4.d0
    call ss%set_A(2,3,  8.d0); A(2,3) = 8.d0; A(3,2) = 8.d0
    call ss%set_A(2,5,  2.d0); A(2,5) = 2.d0; A(5,2) = 2.d0
    call ss%set_A(3,3,  1.d0); A(3,3) = 1.d0
    call ss%set_A(3,8,  5.d0); A(3,8) = 5.d0; A(8,3) = 5.d0
    call ss%set_A(4,4,  7.d0); A(4,4) = 7.d0
    call ss%set_A(4,7,  9.d0); A(4,7) = 9.d0; A(7,4) = 9.d0
    call ss%set_A(5,5,  5.d0); A(5,5) = 5.d0
    call ss%set_A(5,6,  1.d0); A(5,6) = 1.d0; A(6,5) = 1.d0
    call ss%set_A(5,7,  5.d0); A(5,7) = 5.d0; A(7,5) = 5.d0
    call ss%set_A(6,6,  0.d0); A(6,6) = 0.d0
    call ss%set_A(6,8,  5.d0); A(6,8) = 5.d0; A(8,6) = 5.d0
    call ss%set_A(7,7, 11.d0); A(7,7) = 11.d0
    call ss%set_A(8,8,  5.d0); A(8,8) = 5.d0

!!$    write(*,*) 'n = ', n, ', ni = ', ss%ni, ', nj = ', ss%nj
!!$    write(*,*) 'ia = ', ss%ia(1:ss%ni+1)
!!$    write(*,*) 'ja = ', ss%ja(1:ss%nj)
!!$    write(*,*) 'A  = ', ss%a(1:ss%nj)
!!$    write(*,*) 'b  = ', ss%b(1:ss%n)
!!$    stop
    
    ! Decompose matrix
    call ss%decomp()

    ! Set up RHS
    allocate(b(n), x(n))
    b = 1.d0
    call ss%set_rhs(b)
    
    ! Solve system
    call ss%solve(x)

    ! Output solution
    write(*,*) 'x(sparse) = ', x

    ! Solve brute-force
    call solve_system_real(A, x, b)
    write(*,*) 'x(dense) = ', x
    
    ! Clean up
    call ss%dealloc()
        
  end subroutine sparse_system_tester

end module comm_sparse_mod
