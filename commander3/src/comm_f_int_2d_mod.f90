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
module comm_F_int_2D_mod
  use comm_F_int_mod
  use comm_comp_mod
  use comm_bp_mod
  use spline_2D_mod
  use comm_utils
  implicit none

  private
  public comm_F_int_2D

  type, extends (comm_F_int) :: comm_F_int_2D
     class(comm_comp), pointer :: comp => null()
     real(dp), allocatable, dimension(:)       :: x, y
     real(dp), allocatable, dimension(:,:,:,:) :: coeff
   contains
     ! Data procedures
     procedure :: eval       => evalIntF
     procedure :: eval_deriv => evalIntdF
     procedure :: update     => updateIntF
  end type comm_F_int_2D

  interface comm_F_int_2D
     procedure constructor
  end interface comm_F_int_2D

  integer(i4b) :: n = 100
  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(comp, bp, pol)
    implicit none
    class(comm_comp),     intent(in), target :: comp
    class(comm_bp),       intent(in), target :: bp
    integer(i4b),         intent(in)         :: pol
    class(comm_F_int_2D), pointer            :: constructor

    integer(i4b) :: i

    allocate(constructor)
    constructor%comp => comp
    constructor%bp   => bp

    allocate(constructor%x(n), constructor%y(n), constructor%coeff(4,4,n,n))
    do i = 1, n
       constructor%x(i) = constructor%comp%p_uni(1,1) + (constructor%comp%p_uni(2,1)-constructor%comp%p_uni(1,1))/(n-1) * (i-1)
       constructor%y(i) = constructor%comp%p_uni(1,2) + (constructor%comp%p_uni(2,2)-constructor%comp%p_uni(1,2))/(n-1) * (i-1)
    end do

    call constructor%update(pol=pol)

    
  end function constructor

  
  ! Evaluate SED in brightness temperature normalized to reference frequency
  function evalIntF(self, theta)
    class(comm_F_int_2D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    real(dp)                                        :: evalIntF
    evalIntF = splin2_full_precomp(self%x, self%y, self%coeff, theta(1), theta(2)) 
  end function evalIntF

  ! Evaluate partial derivative of SED in brightness temperature normalized to reference frequency
  function evalIntdF(self, theta, par)
    class(comm_F_int_2D),                intent(in) :: self
    real(dp),             dimension(1:), intent(in) :: theta
    integer(i4b),                        intent(in) :: par
    real(dp)                                        :: evalIntdF

    real(dp) :: p(2), delta = 1.d-10, f1, f2

    p      = theta
    p(par) = p(par) + delta
    f1     = self%eval(theta)
    f2     = self%eval(p)
    evalIntdF = (f2-f1)/delta
  end function evalIntdF

  ! Compute/update integration look-up tables
  subroutine updateIntF(self, f_precomp, pol)
    class(comm_F_int_2D), intent(inout)           :: self
    real(dp),             intent(in),    optional :: f_precomp
    integer(i4b),         intent(in),    optional :: pol

    integer(i4b) :: i, j, k, m, ierr
    real(dp)     :: t1, t2, t3, t4
    real(dp), allocatable, dimension(:)   :: s
    real(dp), allocatable, dimension(:,:) :: F, Fp

    call wall_time(t1)

    m = self%bp%n
    allocate(F(n,n), Fp(n,n), s(m))

    ! Evaluate the bandpass integrated SED over all relevant parameters
    F = 0.d0
    call wall_time(t3)
    do i = 1+self%comp%myid, n, self%comp%numprocs
       do j = 1, n
          do k = 1, m
             s(k) = self%comp%S(nu=self%bp%nu(k), pol=pol, theta=[self%x(i),self%y(j)])
          end do
          F(i,j) = self%bp%SED2F(s)
       end do
    end do
    call wall_time(t4)
    !if (self%comp%myid == 0) write(*,*) 'a = ', real(t4-t3,sp)
    call wall_time(t3)
    call mpi_allreduce(F, Fp,            n*n, MPI_DOUBLE_PRECISION, MPI_SUM, self%comp%comm, ierr)
    call wall_time(t4)
    !if (self%comp%myid == 0) write(*,*) 'b = ', real(t4-t3,sp)

    ! Precompute spline object
    call wall_time(t3)
    call splie2_full_precomp_mpi(self%comp%comm, self%x, self%y, Fp, self%coeff)
    call wall_time(t4)
    !if (self%comp%myid == 0) write(*,*) 'c = ', real(t4-t3,sp)
    
    call wall_time(t2)
    !if (comp%myid == 0) write(*,*) 'sum = ', sum(abs(constructor%coeff)), real(t2-t1,sp)
    !call mpi_finalize(i)
    !stop

    deallocate(F, Fp, s)

  end subroutine updateIntF


end module comm_F_int_2D_mod
