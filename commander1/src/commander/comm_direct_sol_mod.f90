module comm_direct_sol_mod
  use comm_mp_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b),     private :: n

contains


  ! Routine for solving the big, difficult linear system,   Au = v,  Ax = rhs
  subroutine solve_fundamental_system_direct(rhs, x_out, stat)
    implicit none
    
    type(genvec), intent(in)    :: rhs
    type(genvec), intent(inout) :: x_out
    integer(i4b), intent(inout) :: stat

    integer(i4b)     :: i, j, k, l, m
    type(genvec)     :: u, Au
    real(dp)         :: epsilon

    real(dp), allocatable, dimension(:,:)   :: A_mat, A_mat2
    real(dp), allocatable, dimension(:)     :: b, x, y

    logical(lgt), save :: first_call = .true.
    integer(i4b), allocatable, dimension(:), save :: ind

    if (first_call) then
       m          = n
       allocate(ind(n))
       do i = 1, n
          ind(i) = i
       end do
    else
       m          = size(ind)
    end if

    ! Set up the coefficient matrix
    allocate(A_mat(m,m), x(m), y(n), b(m))
    call allocate_genvec(u)
    call allocate_genvec(Au)
       
    A_mat = 0.d0
    do i = 1, m
       if (mod(i,1000) == 0) write(*,fmt='(a,i8,a,i8)') 'Computing column ', i, ' of ', m
       x      = 0.d0; y = 0.d0
       x(i)   = 1.d0
       y(ind) = x
       call linvec2genvec(y, u)
       call compute_Au(u, Au)
       call genvec2linvec(Au,y)
       A_mat(i,:) = y(ind)
    end do

    if (first_call) then
       m = 0
       do i = 1, n
          if (A_mat(i,i) > 0.d0) m = m+1
       end do
       deallocate(ind)
       allocate(ind(m))
       m = 0
       do i = 1, n
          if (A_mat(i,i) > 0.d0) then
             m = m+1
             ind(m) = i
          end if
       end do
       allocate(A_mat2(m,m))
       A_mat2 = A_mat(ind,ind)
       deallocate(A_mat, b, x)
       allocate(A_mat(m,m), b(m), x(m))
       A_mat = A_mat2
       deallocate(A_mat2)
       first_call = .false.
    end if

!!$    ! Check for symmetry
!!$    do i = 1, m
!!$       do j = 1, i
!!$          epsilon = (A_mat(i,j)-A_mat(j,i)) / abs(A_mat(i,j) + A_mat(j,i))
!!$          if (epsilon > 1.d-6) then
!!$             write(*,*) i, j, real(A_mat(i,j),sp), real(abs(A_mat(i,j) + A_mat(j,i)),sp), &
!!$                  & real(epsilon,sp)
!!$          end if
!!$       end do
!!$    end do

    ! Set up right-hand side
    call genvec2linvec(rhs, y)
    b = y(ind)

    ! Compute Cholesky solution
    write(*,*) '   Solving system by Cholesky decomposition'
    call cholesky_decompose_single(A_mat, stat)
    if (stat /= 0) write(*,*) 'Error in Cholesky -- stat = ', stat
    call cholesky_solve(A_mat, b, x)
    y(ind) = x
    call linvec2genvec(y, x_out)

    if (.true.) then
       open(58,file='sol_brute.dat')
       do j = 1, n
          write(58,*) j, y(j)
       end do
       close(58)
       !stop
    end if

    ! Output results and change to physical variables
    call multiply_by_sqrt_S(.false., x_out%cmb_amp)
    !call convert_x2amp(x_out%fg_amp)
    !call apply_P_trans_constraint(x_out)

    ! Clean up
    call deallocate_genvec(u)
    call deallocate_genvec(Au)
    deallocate(x, y, b, A_mat)
    
  end subroutine solve_fundamental_system_direct


  subroutine initialize_direct_solution_mod(paramfile)
    implicit none
    
    character(len=128), intent(in) :: paramfile

    real(dp), allocatable, dimension(:) :: x
    type(genvec) :: s

    ! Find number of elements in unknown vector
    call allocate_genvec(s)
    call genvec2linvec(s, x)
    n = size(x)
    deallocate(x)
    call deallocate_genvec(s)

  end subroutine initialize_direct_solution_mod

end module comm_direct_sol_mod
