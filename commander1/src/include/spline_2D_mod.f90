module spline_2D_mod
  use healpix_types
  use spline_1D_mod
  use math_tools
  implicit none


contains

  ! -----------------------------------------------------------    
  ! Routines for full precomputation 2D spline interpolation
  !-------------------------------------------------------------

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
    real(dp), dimension(4) :: u_pow, v_pow, Av

    m = size(x)
    n = size(y)

!    inv_h_x = 1.d0 / (x(2)-x(1))
!    inv_h_y = 1.d0 / (y(2)-y(1))

    b_x = locate(x, x0)
    b_y = locate(y, y0)
!    b_x = max(min(int((x0-x(1))*inv_h_x)+1,m-1),1)
!    b_y = max(min(int((y0-y(1))*inv_h_y)+1,n-1),1)

    inv_h_x = 1.d0 / (x(b_x+1)-x(b_x))
    inv_h_y = 1.d0 / (y(b_y+1)-y(b_y))

!   if (x0 < x(1) .or. x0 > x(n) then
    if (b_x <= 0 .or. b_x >= m) then
       write(*,*) 'splin2_full_precomp -- Warning: x0 out of bounds = ', x0
       splin2_full_precomp = 0.d0
       return
    end if

!   if (y0 < y(1) .or. y0 > y(n) then
    if (b_y < 0 .or. b_y > n) then
       write(*,*) 'splin2_full_precomp -- Warning: y0 out of bounds = ', y0
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


  ! -----------------------------------------------------------    
  ! Routines for direct low-precomputation 2D spline interpolation
  !-------------------------------------------------------------

  function splin2(x1a, x2a, ya, y2a, x1, x2)
    implicit none
     
    real(dp), dimension(1:),    intent(in)  :: x1a, x2a
    real(dp), dimension(1:,1:), intent(in)  :: ya, y2a
    real(dp),                   intent(in)  :: x1, x2
    real(dp)                                :: splin2

    integer(i4b) :: j, k, n, m
    real(dp), allocatable, dimension(:) :: y2tmp, yytmp

    m = size(ya(:,1))
    n = size(ya(1,:))
    
    allocate(yytmp(m))
    allocate(y2tmp(m))

    do j = 1, m
       yytmp(j) = splint(x2a, ya(j,:), y2a(j,:), x2)
    end do

    call spline(x1a, yytmp, 1.d30, 1.d30, y2tmp)
    splin2 = splint(x1a, yytmp, y2tmp, x1)

    deallocate(yytmp)
    deallocate(y2tmp)
    
  end function splin2


  subroutine splie2(x1a, x2a, ya, y2a)
    implicit none
     
    real(dp), dimension(1:),    intent(in)  :: x1a, x2a
    real(dp), dimension(1:,1:), intent(in)  :: ya
    real(dp), dimension(1:,1:), intent(out) :: y2a

    integer(i4b) :: j

    do j = 1, size(ya(:,1))
       call spline(x2a, ya(j,:), 1.d30, 1.d30, y2a(j,:))
    end do
    
  end subroutine splie2


end module spline_2D_mod
