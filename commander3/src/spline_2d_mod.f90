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
module spline_2D_mod
  use healpix_types
  use spline_1D_mod
  use math_tools
  use comm_utils
  implicit none


contains

  ! -----------------------------------------------------------    
  ! Routines for full precomputation 2D spline interpolation
  !-------------------------------------------------------------

  subroutine splie2_full_precomp_mpi(comm, x, y, f, coeff)
    implicit none
     
    integer(i4b),                     intent(in)  :: comm
    real(dp), dimension(1:),          intent(in)  :: x, y
    real(dp), dimension(1:,1:),       intent(in)  :: f
    real(dp), dimension(1:,1:,1:,1:), intent(out) :: coeff

    real(dp)     :: h_x, h_y, t1, t2
    integer(i4b) :: i, j, k, m, n, numprocs, myid, ierr
    real(dp), allocatable, dimension(:)   :: y2
    real(dp), allocatable, dimension(:,:) :: df_dx, df_dy, ddf_dxdy
    real(dp), allocatable, dimension(:,:,:) :: buffer

    real(dp),              dimension(16,16) :: a
    real(dp),              dimension(16)    :: b, z

    call mpi_comm_rank(comm, myid,     ierr)
    call mpi_comm_size(comm, numprocs, ierr)

    m = size(x)
    n = size(y)

    ! Compute df/dx at all nodes
    call wall_time(t1)
    allocate(y2(m))
    allocate(df_dx(m,n))
    df_dx = 0.d0
    do i = 1+myid, n, numprocs
       call spline(x, f(:,i), 1.d30, 1.d30, y2)
       call splint_deriv_all_nodes(x, f(:,i), y2, df_dx(:,i))
    end do
    deallocate(y2)
    call mpi_allreduce(MPI_IN_PLACE, df_dx, m*n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    call wall_time(t2)
    !if (myid == 0) write(*,*) 'q1 = ', real(t2-t1,sp)

    ! Compute df/dy at all nodes
    call wall_time(t1)
    allocate(y2(n))
    allocate(df_dy(m,n))
    df_dy = 0.d0
    do i = 1+myid, m, numprocs
       call spline(y, f(i,:), 1.d30, 1.d30, y2)
       call splint_deriv_all_nodes(y, f(i,:), y2, df_dy(i,:))
    end do
    deallocate(y2)
    call mpi_allreduce(MPI_IN_PLACE, df_dy, m*n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    call wall_time(t2)
    !if (myid == 0) write(*,*) 'q2 = ', real(t2-t1,sp)


    ! Compute ddf/dxdy at all nodes
    call wall_time(t2)
    allocate(y2(m))
    allocate(ddf_dxdy(m,n))
    ddf_dxdy = 0.d0
    do i = 1+myid, n, numprocs
       call spline(x, df_dy(:,i), 1.d30, 1.d30, y2)
       call splint_deriv_all_nodes(x, df_dy(:,i), y2, ddf_dxdy(:,i))
    end do
    deallocate(y2)
    call mpi_allreduce(MPI_IN_PLACE, ddf_dxdy, m*n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    call wall_time(t2)
    !if (myid == 0) write(*,*) 'q3 = ', real(t2-t1,sp)

    ! Solve for spline coefficients in each grid cell
    call wall_time(t1)
    !coeff = 0.d0
    do i = 1+myid, m-1, numprocs
       do j = 1, n-1

          ! See http://en.wikipedia.org/wiki/Bicubic_interpolation for equations
          A        = 0.d0
          h_x      = x(i+1)-x(i)
          h_y      = y(j+1)-y(j)

          ! Matching function values at nodes
          A(1,1)   = 1.d0
          b(1)     = f(i,j)

          A(2,1)   = 1.d0
          A(2,2)   = 1.d0
          A(2,3)   = 1.d0
          A(2,4)   = 1.d0
          b(2)     = f(i+1,j)

          A(3,1)   = 1.d0
          A(3,5)   = 1.d0
          A(3,9)   = 1.d0
          A(3,13)  = 1.d0
          b(3)     = f(i,j+1)
          
          A(4,:)   = 1.d0
          b(4)     = f(i+1,j+1)
          
          ! Matching x-derivatives at each node
          A(5,2)   = 1.d0
          b(5)     = df_dx(i,j)      * h_x

          A(6,2)   = 1.d0
          A(6,3)   = 2.d0
          A(6,4)   = 3.d0
          b(6)     = df_dx(i+1,j)    * h_x

          A(7,2)   = 1.d0
          A(7,6)   = 1.d0
          A(7,10)  = 1.d0
          A(7,14)  = 1.d0
          b(7)     = df_dx(i,j+1)    * h_x
          
          do k = 1, 16
             A(8,k)   = mod(k-1,4)
          end do
          b(8)     = df_dx(i+1,j+1)  * h_x

          ! Matching y-derivatives at each node
          A(9,5)   = 1.d0
          b(9)     = df_dy(i,j)      * h_y

          A(10,5)   = 1.d0
          A(10,6)   = 1.d0
          A(10,7)   = 1.d0
          A(10,8)   = 1.d0
          b(10)     = df_dy(i+1,j)   * h_y

          A(11,5)   = 1.d0
          A(11,9)   = 2.d0
          A(11,13)  = 3.d0
          b(11)     = df_dy(i,j+1)   * h_y
          
          do k = 1, 16
             A(12,k)   = (k-1)/4
          end do
          b(12)     = df_dy(i+1,j+1) * h_y
          
          ! Matching dxdy derivatives at each node
          A(13,6)   = 1.d0
          b(13)     = ddf_dxdy(i,j)     * (h_x*h_y)

          A(14,6)   = 1.d0
          A(14,7)   = 2.d0
          A(14,8)   = 3.d0
          b(14)     = ddf_dxdy(i+1,j)   * (h_x*h_y)

          A(15,6)    = 1.d0
          A(15,10)   = 2.d0
          A(15,14)   = 3.d0
          b(15)      = ddf_dxdy(i,j+1)  * (h_x*h_y)

          do k = 1, 16
             A(16,k)   = mod(k-1,4) * ((k-1)/4)
          end do
          b(16)     = ddf_dxdy(i+1,j+1) * (h_x*h_y)
          

          ! Solve linear system
          call solve_system_real(A, z, b)

          ! Re-organize coefficients into 4x4 ordering
          do k = 1, 4
             coeff(:,k,i,j) = z(4*(k-1)+1:4*k)
          end do

       end do
    end do
    call wall_time(t2)
    !if (myid == 0) write(*,*) 'q4 = ', real(t2-t1,sp)

    call wall_time(t1)    
    allocate(buffer(4,4,n))
    !call mpi_allreduce(coeff, buffer, size(coeff), MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    do j = 0, numprocs-1
       do i = 1+j, m-1, numprocs
          if (myid == j) buffer = coeff(:,:,i,:)
          call mpi_bcast(buffer, 16*n, MPI_DOUBLE_PRECISION, j, comm, ierr)
          coeff(:,:,i,:) = buffer
       end do
    end do
    !coeff = buffer
    call wall_time(t2)
    deallocate(buffer)
    !if (myid == 0) write(*,*) 'q5 = ', real(t2-t1,sp)



    deallocate(df_dx)
    deallocate(df_dy)
    deallocate(ddf_dxdy)
    
  end subroutine splie2_full_precomp_mpi


  subroutine splie2_full_precomp(x, y, f, coeff)
    implicit none
     
    real(dp), dimension(1:),          intent(in)  :: x, y
    real(dp), dimension(1:,1:),       intent(in)  :: f
    real(dp), dimension(1:,1:,1:,1:), intent(out) :: coeff

    real(dp)     :: h_x, h_y
    integer(i4b) :: i, j, k, m, n
    real(dp), allocatable, dimension(:)   :: y2
    real(dp), allocatable, dimension(:,:) :: df_dx, df_dy, ddf_dxdy

    real(dp),              dimension(16,16) :: a
    real(dp),              dimension(16)    :: b, z

    m = size(x)
    n = size(y)

    ! Compute df/dx at all nodes
    allocate(y2(m))
    allocate(df_dx(m,n))
    do i = 1, n
       call spline(x, f(:,i), 1.d30, 1.d30, y2)
       call splint_deriv_all_nodes(x, f(:,i), y2, df_dx(:,i))
    end do
    deallocate(y2)

    ! Compute df/dy at all nodes
    allocate(y2(n))
    allocate(df_dy(m,n))
    do i = 1, m
       call spline(y, f(i,:), 1.d30, 1.d30, y2)
       call splint_deriv_all_nodes(y, f(i,:), y2, df_dy(i,:))
    end do
    deallocate(y2)


    ! Compute ddf/dxdy at all nodes
    allocate(y2(m))
    allocate(ddf_dxdy(m,n))
    do i = 1, n
       call spline(x, df_dy(:,i), 1.d30, 1.d30, y2)
       call splint_deriv_all_nodes(x, df_dy(:,i), y2, ddf_dxdy(:,i))
    end do
    deallocate(y2)


    ! Solve for spline coefficients in each grid cell
    do i = 1, m-1
       do j = 1, n-1

          ! See http://en.wikipedia.org/wiki/Bicubic_interpolation for equations
          A        = 0.d0
          h_x      = x(i+1)-x(i)
          h_y      = y(j+1)-y(j)

          ! Matching function values at nodes
          A(1,1)   = 1.d0
          b(1)     = f(i,j)

          A(2,1)   = 1.d0
          A(2,2)   = 1.d0
          A(2,3)   = 1.d0
          A(2,4)   = 1.d0
          b(2)     = f(i+1,j)

          A(3,1)   = 1.d0
          A(3,5)   = 1.d0
          A(3,9)   = 1.d0
          A(3,13)  = 1.d0
          b(3)     = f(i,j+1)
          
          A(4,:)   = 1.d0
          b(4)     = f(i+1,j+1)
          
          ! Matching x-derivatives at each node
          A(5,2)   = 1.d0
          b(5)     = df_dx(i,j)      * h_x

          A(6,2)   = 1.d0
          A(6,3)   = 2.d0
          A(6,4)   = 3.d0
          b(6)     = df_dx(i+1,j)    * h_x

          A(7,2)   = 1.d0
          A(7,6)   = 1.d0
          A(7,10)  = 1.d0
          A(7,14)  = 1.d0
          b(7)     = df_dx(i,j+1)    * h_x
          
          do k = 1, 16
             A(8,k)   = mod(k-1,4)
          end do
          b(8)     = df_dx(i+1,j+1)  * h_x

          ! Matching y-derivatives at each node
          A(9,5)   = 1.d0
          b(9)     = df_dy(i,j)      * h_y

          A(10,5)   = 1.d0
          A(10,6)   = 1.d0
          A(10,7)   = 1.d0
          A(10,8)   = 1.d0
          b(10)     = df_dy(i+1,j)   * h_y

          A(11,5)   = 1.d0
          A(11,9)   = 2.d0
          A(11,13)  = 3.d0
          b(11)     = df_dy(i,j+1)   * h_y
          
          do k = 1, 16
             A(12,k)   = (k-1)/4
          end do
          b(12)     = df_dy(i+1,j+1) * h_y
          
          ! Matching dxdy derivatives at each node
          A(13,6)   = 1.d0
          b(13)     = ddf_dxdy(i,j)     * (h_x*h_y)

          A(14,6)   = 1.d0
          A(14,7)   = 2.d0
          A(14,8)   = 3.d0
          b(14)     = ddf_dxdy(i+1,j)   * (h_x*h_y)

          A(15,6)    = 1.d0
          A(15,10)   = 2.d0
          A(15,14)   = 3.d0
          b(15)      = ddf_dxdy(i,j+1)  * (h_x*h_y)

          do k = 1, 16
             A(16,k)   = mod(k-1,4) * ((k-1)/4)
          end do
          b(16)     = ddf_dxdy(i+1,j+1) * (h_x*h_y)
          

          ! Solve linear system
          call solve_system_real(A, z, b)

          ! Re-organize coefficients into 4x4 ordering
          do k = 1, 4
             coeff(:,k,i,j) = z(4*(k-1)+1:4*k)
          end do

       end do
    end do


    deallocate(df_dx)
    deallocate(df_dy)
    deallocate(ddf_dxdy)
    
  end subroutine splie2_full_precomp


  function splin2_full_precomp(x, y, coeff, x0, y0)
    implicit none
     
    real(dp), dimension(1:),          intent(in)  :: x, y
    real(dp), dimension(1:,1:,1:,1:), intent(in)  :: coeff
    real(dp),                         intent(in)  :: x0, y0
    real(dp)                                      :: splin2_full_precomp

    integer(i4b) :: i, j, b_x, b_y, m, n
    real(dp)     :: u, v, s, inv_h_x, inv_h_y
    real(dp), dimension(4) :: u_pow, v_pow

    m = size(x)
    n = size(y)

    inv_h_x = 1.d0 / (x(2)-x(1))
    inv_h_y = 1.d0 / (y(2)-y(1))

!    b_x = locate(x, x0)
!    b_y = locate(y, y0)
    b_x = max(min(int((x0-x(1))*inv_h_x)+1,m-1),1)
    b_y = max(min(int((y0-y(1))*inv_h_y)+1,n-1),1)


    if (x0 < x(1) .or. x0-x(m) > 1.d-12) then
       !write(*,fmt='(a,3f8.3)') 'splin2_full_precomp -- Warning: x0 out of bounds = ', x(1), x0, x(m),
       write(*,*) 'splin2_full_precomp -- Warning: x0 out of bounds = ', x(1), x0, x(m), x0-x(m)
       splin2_full_precomp = 0.d0
       return
    end if

    if (y0 < y(1) .or. y0-y(n) > 1d-12) then
       write(*,fmt='(a,3f8.3)') 'splin2_full_precomp -- Warning: y0 out of bounds = ', y(1), y0, y(n)
       splin2_full_precomp = 0.d0
       return
    end if

    ! Compute spline value
    u = (x0-x(b_x)) * inv_h_x
    v = (y0-y(b_y)) * inv_h_y

    u_pow(1) = 1.d0
    u_pow(2) = u
    u_pow(3) = u_pow(2) * u
    u_pow(4) = u_pow(3) * u

    v_pow(1) = 1.d0
    v_pow(2) = v
    v_pow(3) = v_pow(2) * v
    v_pow(4) = v_pow(3) * v

    s = 0.d0
    do j = 1, 4
       do i = 1, 4
          s = s + coeff(i,j,b_x,b_y) * u_pow(i) * v_pow(j)
       end do
    end do

    splin2_full_precomp = s

  end function splin2_full_precomp


  function splin2_full_precomp_irreg(x, y, coeff, x0, y0)
    implicit none
     
    real(dp), dimension(1:),          intent(in)  :: x, y
    real(dp), dimension(1:,1:,1:,1:), intent(in)  :: coeff
    real(dp),                         intent(in)  :: x0, y0
    real(dp)                                      :: splin2_full_precomp_irreg

    integer(i4b) :: i, j, b_x, b_y, m, n
    real(dp)     :: u, v, s, inv_h_x, inv_h_y
    real(dp), dimension(4) :: u_pow, v_pow

    m = size(x)
    n = size(y)

    b_x = locate(x, x0)
    b_y = locate(y, y0)
    !b_x = max(min(int((x0-x(1))*inv_h_x)+1,m-1),1)
    !b_y = max(min(int((y0-y(1))*inv_h_y)+1,n-1),1)


    if (x0 < x(1) .or. x0 > x(m)) then
       write(*,fmt='(a,3f8.3)') 'splin2_full_precomp_irreg -- Warning: x0 out of bounds = ', x(1), x0, x(m)
       splin2_full_precomp_irreg = 0.d0
       return
    end if

    if (y0 < y(1) .or. y0 > y(n)) then
       write(*,fmt='(a,3f8.3)') 'splin2_full_precomp_irreg -- Warning: y0 out of bounds = ', y(1), y0, y(n)
       splin2_full_precomp_irreg = 0.d0
       return
    end if

    !inv_h_x = 1.d0 / (x(2)-x(1))
    !inv_h_y = 1.d0 / (y(2)-y(1))
    if (b_x == m) then
       inv_h_x = 1.d0 / (x(b_x)-x(b_x-1))
    else
       inv_h_x = 1.d0 / (x(b_x+1)-x(b_x))
    end if

    if (b_y == n) then
       inv_h_y = 1.d0 / (y(b_y)-y(b_y-1))
    else
       inv_h_y = 1.d0 / (y(b_y+1)-y(b_y))
    end if

    ! Compute spline value
    u = (x0-x(b_x)) * inv_h_x
    v = (y0-y(b_y)) * inv_h_y

    u_pow(1) = 1.d0
    u_pow(2) = u
    u_pow(3) = u_pow(2) * u
    u_pow(4) = u_pow(3) * u

    v_pow(1) = 1.d0
    v_pow(2) = v
    v_pow(3) = v_pow(2) * v
    v_pow(4) = v_pow(3) * v

    s = 0.d0
    do j = 1, 4
       do i = 1, 4
          s = s + coeff(i,j,b_x,b_y) * u_pow(i) * v_pow(j)
       end do
    end do

    splin2_full_precomp_irreg = s

  end function splin2_full_precomp_irreg

end module spline_2D_mod
