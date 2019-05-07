module comm_S_mult_mod
  use healpix_types
  use math_tools
  use comm_utils
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b),                                private :: lmax, nmaps
  real(dp),     allocatable, dimension(:,:,:), private :: S, sqrt_S, inv_sqrt_S

contains

  subroutine initialize_S_mult_mod(lmax_in, nmaps_in)
    implicit none

    integer(i4b), intent(in) :: lmax_in, nmaps_in

    lmax  = lmax_in
    nmaps = nmaps_in

    allocate(S(nmaps, nmaps, 0:lmax))
    allocate(sqrt_S(nmaps, nmaps, 0:lmax))
    allocate(inv_sqrt_S(nmaps, nmaps, 0:lmax))

  end subroutine initialize_S_mult_mod

  subroutine clean_up_S_mult_mod
    implicit none

    if (allocated(S))          deallocate(S)
    if (allocated(sqrt_S))     deallocate(sqrt_S)
    if (allocated(inv_sqrt_S)) deallocate(inv_sqrt_S)

  end subroutine clean_up_S_mult_mod


  subroutine update_S_mult_mod(cls)
    implicit none

    real(dp), dimension(0:,1:), intent(in) :: cls

    integer(i4b) :: i, j, l, info
    real(dp), allocatable, dimension(:)     :: w
    real(dp), allocatable, dimension(:,:)   :: A
    real(dp),              dimension(1:300) :: work

    logical(lgt), dimension(3) :: acceptable_mode

    allocate(w(nmaps))
    allocate(A(nmaps, nmaps))

    S          = 0.d0
    sqrt_S     = 0.d0
    inv_sqrt_S = 0.d0
    do l = 0, lmax
       if (l <= 1) then
          if (any(cls(:,1) /= 0.d0)) then
             S(1,1,l) = 1.d12
             !S(1,1,l) = 0.d0
          else
             S(1,1,l) = 0.d0
          end if
       else
          call cl2s(cls(l,:), S(:,:,l))
          S(:,:,l) = S(:,:,l) / (l*(l+1)/(2.d0*pi))
       end if

       if (nmaps == 1) then
          sqrt_S(1,1,l)     = sqrt(S(1,1,l))
          if (sqrt_S(1,1,l) > 0.d0) inv_sqrt_S(1,1,l) = 1.d0 / sqrt_S(1,1,l)
       else

          acceptable_mode = .true.
          do i = 1, nmaps
             if (S(i,i,l) <= 0.d0) then
                acceptable_mode(i) = .false.
                S(i,:,l) = 0.d0
                S(:,i,l) = 0.d0
                S(i,i,l) = 1.d0
             end if
          end do

          sqrt_S(:,:,l) = S(:,:,l)
          call compute_hermitian_root(sqrt_S(:,:,l), 0.5d0)
          inv_sqrt_S(:,:,l) = sqrt_S(:,:,l)
          call invert_matrix(inv_sqrt_S(:,:,l))

          do i = 1, nmaps
             if (.not. acceptable_mode(i)) then
                S(i,:,l)          = 0.d0
                S(:,i,l)          = 0.d0
                sqrt_S(i,:,l)     = 0.d0
                sqrt_S(:,i,l)     = 0.d0
                inv_sqrt_S(i,:,l) = 0.d0
                inv_sqrt_S(:,i,l) = 0.d0
             end if
          end do

       end if
    end do

    deallocate(w)
    deallocate(A)

  end subroutine update_S_mult_mod


  subroutine multiply_by_diag_S(alms_in, diag, alms_out)
    implicit none

    integer(i4b),                intent(in)  :: diag
    real(dp),     dimension(1:), intent(in)  :: alms_in
    real(dp),     dimension(1:), intent(out) :: alms_out

    integer(i4b) :: l, m, ind

    do l = 0, lmax
       do m = -l, l
          ind = lm2ind(l,m)
          alms_out(ind) = S(diag,diag,l) * alms_in(ind)
       end do
    end do

  end subroutine multiply_by_diag_S


  subroutine multiply_by_inv_sqrt_S_single_l(l, alms)
    implicit none

    integer(i4b),                intent(in)    :: l
    real(dp),     dimension(1:), intent(inout) :: alms

    alms = matmul(inv_sqrt_S(:,:,l), alms)

  end subroutine multiply_by_inv_sqrt_S_single_l

  subroutine multiply_by_inv_sqrt_S(trans, alms_in, alms_out, lmax_int)
    implicit none

    logical(lgt),                   intent(in)              :: trans
    real(dp),     dimension(1:,1:), intent(inout)           :: alms_in
    real(dp),     dimension(1:,1:), intent(inout), optional :: alms_out
    integer(i4b),                   intent(in),    optional :: lmax_int

    integer(i4b) :: l, m, n, ind, ind1, ind2, lmax_
    real(dp),     allocatable, dimension(:,:) :: y
    integer(i4b), allocatable, dimension(:)   :: mask

    if (present(lmax_int)) then
       lmax_ = lmax_int
    else
       lmax_ = lmax
    end if

    allocate(y(2*lmax_+1,nmaps))

    n = 0
    do m = 1, nmaps
       if (any(inv_sqrt_S(m,m,:) /= 0.d0)) n = n+1
    end do

    if (n == 0) then
       if (present(alms_out)) then
          alms_out = 0.d0
       else
          alms_in  = 0.d0
       end if
       return
    end if

    allocate(mask(n))
    n = 0
    do m = 1, nmaps
       if (any(inv_sqrt_S(m,m,:) /= 0.d0)) then
          n = n+1
          mask(n) = m
       end if
    end do

    ! Note: Matrices are assumed symmetric!
    ind = 0
    do l = 0, lmax_
       ind1 = lm2ind(l,-l)
       ind2 = lm2ind(l, l)
       y(1:2*l+1,:) = 0.d0
       y(1:2*l+1,mask) = matmul(alms_in(ind1:ind2,mask), inv_sqrt_S(mask,mask,l))
       if (present(alms_out)) then
          alms_out(ind1:ind2,:) = y(1:2*l+1,:)
       else
          alms_in(ind1:ind2,:) = y(1:2*l+1,:)
       end if
    end do

    deallocate(y)

  end subroutine multiply_by_inv_sqrt_S




  subroutine multiply_by_sqrt_S(trans, alms_in, alms_out, lmax_int)
    implicit none

    logical(lgt),                   intent(in)              :: trans
    real(dp),     dimension(1:,1:), intent(inout)           :: alms_in
    real(dp),     dimension(1:,1:), intent(inout), optional :: alms_out
    integer(i4b),                   intent(in),    optional :: lmax_int

    integer(i4b) :: l, m, n, ind, ind1, ind2, lmax_
    real(dp),     allocatable, dimension(:,:) :: y
    integer(i4b), allocatable, dimension(:)   :: mask

    if (present(lmax_int)) then
       lmax_ = lmax_int
    else
       lmax_ = lmax
    end if

    allocate(y(2*lmax_+1,nmaps))

    n = 0
    do m = 1, nmaps
       if (any(sqrt_S(m,m,:) /= 0.d0)) n = n+1
    end do

    if (n == 0) then
       if (present(alms_out)) then
          alms_out = 0.d0
       else
          alms_in  = 0.d0
       end if
!       write(*,*) 'a should not be here...'
       return
    end if

    allocate(mask(n))
    n = 0
    do m = 1, nmaps
       if (any(sqrt_S(m,m,:) /= 0.d0)) then
          n = n+1
          mask(n) = m
       end if
    end do

    ! Note: Matrices are assumed symmetric!
    ind = 0
    do l = 0, lmax_
       ind1 = lm2ind(l,-l)
       ind2 = lm2ind(l, l)
       y(1:2*l+1,:) = 0.d0
       y(1:2*l+1,mask) = matmul(alms_in(ind1:ind2,mask), sqrt_S(mask,mask,l))
       if (present(alms_out)) then
          alms_out(ind1:ind2,:) = y(1:2*l+1,:)
       else
          alms_in(ind1:ind2,:) = y(1:2*l+1,:)
       end if
    end do

    deallocate(y)

  end subroutine multiply_by_sqrt_S


  subroutine compute_Lt_A_L(l, A)
    implicit none

    integer(i4b),                   intent(in)     :: l
    real(dp),     dimension(1:,1:), intent(inout)  :: A

    integer(i4b) :: ind

    A = matmul(A, sqrt_S(:,:,l))
    A = matmul(transpose(sqrt_S(:,:,l)), A)

  end subroutine compute_Lt_A_L


  function get_Cl(l,i,j)
    implicit none
    
    integer(i4b), intent(in) :: l, i, j
    real(dp)                 :: get_Cl

    get_Cl = S(i,j,l)

  end function get_Cl


end module comm_S_mult_mod
