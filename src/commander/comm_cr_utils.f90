module comm_cr_utils
  use comm_utils
  implicit none

  integer(i4b)  :: ncr
  integer(i4b), allocatable, dimension(:,:) :: ind_comp    ! (start, length, nmaps)

  interface cr_insert_comp
     module procedure cr_insert_comp_1d, cr_insert_comp_2d
  end interface cr_insert_comp

  interface cr_extract_comp
     module procedure cr_extract_comp_1d, cr_extract_comp_2d
  end interface cr_extract_comp
  
contains

  subroutine cr_insert_comp_1d(id, add, x_in, x_out)
    implicit none

    integer(i4b),                  intent(in)    :: id
    logical(lgt),                  intent(in)    :: add
    real(dp),     dimension(1:),   intent(in)    :: x_in
    real(dp),     dimension(1:),   intent(inout) :: x_out

    integer(i4b) :: i, n, pos

    n   = size(x_in,1)
    pos = ind_comp(id,1)
    if (add) then
       x_out(pos:pos+n-1) = x_out(pos:pos+n-1) + x_in
    else
       x_out(pos:pos+n-1) = x_in
    end if
    
  end subroutine cr_insert_comp_1d

  subroutine cr_insert_comp_2d(id, add, x_in, x_out)
    implicit none

    integer(i4b),                   intent(in)    :: id
    logical(lgt),                   intent(in)    :: add
    real(dp),     dimension(0:,1:), intent(in)    :: x_in
    real(dp),     dimension(1:),    intent(inout) :: x_out

    integer(i4b) :: i, n, pos

    n   = size(x_in,1)
    pos = ind_comp(id,1)
    do i = 1, size(x_in,2)
       if (add) then
          x_out(pos:pos+n-1) = x_out(pos:pos+n-1) + x_in(:,i)
       else
          x_out(pos:pos+n-1) = x_in(:,i)
       end if
       pos                = pos+n
    end do
    
  end subroutine cr_insert_comp_2d


  subroutine cr_extract_comp_1d(id, x_in, x_out)
    implicit none

    integer(i4b),                            intent(in)  :: id
    real(dp),     dimension(:),              intent(in)  :: x_in
    real(dp),     dimension(:), allocatable, intent(out) :: x_out

    integer(i4b) :: i, n, pos, nmaps

    pos   = ind_comp(id,1)
    n     = ind_comp(id,2)
    allocate(x_out(n))
    x_out = x_in(pos:pos+n-1)
    
  end subroutine cr_extract_comp_1d

  subroutine cr_extract_comp_2d(id, x_in, x_out)
    implicit none

    integer(i4b),                              intent(in)  :: id
    real(dp),     dimension(:),                intent(in)  :: x_in
    real(dp),     dimension(:,:), allocatable, intent(out) :: x_out

    integer(i4b) :: i, n, pos, nmaps

    pos   = ind_comp(id,1)
    nmaps = ind_comp(id,3)
    n     = ind_comp(id,2)/nmaps
    allocate(x_out(0:n-1,nmaps))
    do i = 1, nmaps
       x_out(:,i) = x_in(pos:pos+n-1)
       pos        = pos + n
    end do
    
  end subroutine cr_extract_comp_2d
  
end module comm_cr_utils
