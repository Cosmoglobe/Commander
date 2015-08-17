module sort_utils
  use healpix_types
  implicit none

contains

!********************
!* Sorting routines *
!********************


! =====================================================================
! Subroutines for sorting. It is a combined quicksort algorithm,
! with cutoff at 10, after which an insertion-sort completes the
! sorting
! =====================================================================

  ! This routine sorts an integer array, according to distances given
  ! by the additional array dist
  subroutine QuickSort(numbers, dist)
    implicit none

    integer(i4b), dimension(:)  :: numbers
    real(dp),     dimension(:)  :: dist

    call quick_sort(numbers, dist, 1, size(numbers))
    call insertion_sort(numbers, dist)
  end subroutine QuickSort

  recursive subroutine quick_sort(numbers, dist, left, right)
    implicit none
    
    integer(i4b), dimension(:)   :: numbers
    real(dp),     dimension(:)   :: dist
    integer(i4b)                 :: left, right

    integer(i4b)  :: i,j, itemp
    real(dp)      :: pivot, rtemp

    if (left+10 < right) then

       call median3(numbers, dist, left, right, pivot)

       i = left; j = right-1

       i = i+1
       do while (dist(i)<pivot)
          i = i+1
       end do
       
       j = j-1
       do while (dist(j)>pivot)
          j = j-1
       end do
       
       itemp = numbers(i)
       rtemp = dist(i)
          
       numbers(i) = numbers(j)
       dist(i)    = dist(j)

       numbers(j) = itemp
       dist(j)    = rtemp

       do while (j>i)
          i = i+1
          do while (dist(i)<pivot)
             i = i+1
          end do
    
          j = j-1
          do while (dist(j)>pivot)
             j = j-1
          end do

          itemp = numbers(i)
          rtemp = dist(i)

          numbers(i) = numbers(j)
          dist(i)    = dist(j)

          numbers(j) = itemp
          dist(j)    = rtemp
       end do

       ! Undo last swap
       itemp = numbers(i)
       rtemp = dist(i)

       numbers(i) = numbers(j)
       dist(i)    = dist(j)

       numbers(j) = itemp
       dist(j)    = rtemp

       ! Restore pivot  
       itemp = numbers(i)
       rtemp = dist(i)

       numbers(i) = numbers(right-1)
       dist(i)    = dist(right-1)

       numbers(right-1) = itemp
       dist(right-1)    = rtemp
          
       call quick_sort(numbers, dist, left, i-1)
       call quick_sort(numbers, dist, i+1, right)
    end if
  end subroutine quick_sort

  subroutine median3(numbers, dist, left, right, pivot)
    implicit none

    integer(i4b)                  :: left, right
    real(dp)                      :: pivot
    integer(i4b), dimension(:)    :: numbers
    real(dp),     dimension(:)    :: dist

    integer(i4b)                  :: center, itemp
    real(dp)                      :: rtemp

    center = (left+right)/2

    if (dist(left)>dist(center)) then
       itemp = numbers(left)
       rtemp = dist(left)

       numbers(left) = numbers(center)
       dist(left)    = dist(center)

       numbers(center) = itemp
       dist(center)    = rtemp
    end if
    
    if (dist(left) > dist(right)) then
       itemp = numbers(left)
       rtemp = dist(left)

       numbers(left) = numbers(right)
       dist(left)    = dist(right)

       numbers(right) = itemp
       dist(right)    = rtemp
    end if

    if (dist(center) > dist(right)) then
       itemp = numbers(center)
       rtemp = dist(center)

       numbers(center) = numbers(right)
       dist(center)    = dist(right)

       numbers(right) = itemp
       dist(right)    = rtemp
    end if

    pivot = dist(center)

    ! Swap the pivot away
    itemp = numbers(center)
    rtemp = dist(center)

    numbers(center) = numbers(right-1)
    dist(center)    = dist(right-1)

    numbers(right-1) = itemp
    dist(right-1)    = rtemp
  end subroutine median3

  subroutine insertion_sort(numbers, dist)
    implicit none

    integer(i4b), dimension(:)   :: numbers
    real(dp),     dimension(:)   :: dist

    integer(i4b)  :: length, i, j
    integer(i4b)  :: itemp
    real(dp)      :: rtemp

    length = size(numbers)

    do i = 2, length
       j = i

       itemp = numbers(i)
       rtemp = dist(i)

       do while (rtemp < dist(j-1))
          dist(j) = dist(j-1)
          numbers(j) = numbers(j-1)
          j = j-1

          if (j == 1) then
             exit
          end if
       end do

       dist(j) = rtemp
       numbers(j) = itemp
    end do
  end subroutine insertion_sort


  !*****************************************************************************


  ! This routine sorts a double-precision array, according to distances given
  ! by the additional array dist
  subroutine QuickSort_dp_dist(numbers, dist)
    implicit none

    real(dp),     dimension(:)  :: numbers, dist

    call quick_sort_dp_dist(numbers, dist, 1, size(numbers))
    call insertion_sort_dp_dist(numbers, dist)
  end subroutine QuickSort_dp_dist

  recursive subroutine quick_sort_dp_dist(numbers, dist, left, right)
    implicit none
    
    real(dp),     dimension(:)   :: numbers, dist
    integer(i4b)                 :: left, right

    integer(i4b)  :: i,j
    real(dp)      :: pivot, rtemp1, rtemp2

    if (left+10 < right) then

       call median3_dp_dist(numbers, dist, left, right, pivot)

       i = left; j = right-1

       i = i+1
       do while (dist(i)<pivot)
          i = i+1
       end do
       
       j = j-1
       do while (dist(j)>pivot)
          j = j-1
       end do
       
       rtemp1 = numbers(i)
       rtemp2 = dist(i)
          
       numbers(i) = numbers(j)
       dist(i)    = dist(j)

       numbers(j) = rtemp1
       dist(j)    = rtemp2

       do while (j>i)
          i = i+1
          do while (dist(i)<pivot)
             i = i+1
          end do
    
          j = j-1
          do while (dist(j)>pivot)
             j = j-1
          end do

          rtemp1 = numbers(i)
          rtemp2 = dist(i)

          numbers(i) = numbers(j)
          dist(i)    = dist(j)

          numbers(j) = rtemp1
          dist(j)    = rtemp2
       end do

       ! Undo last swap
       rtemp1 = numbers(i)
       rtemp2 = dist(i)

       numbers(i) = numbers(j)
       dist(i)    = dist(j)

       numbers(j) = rtemp1
       dist(j)    = rtemp2

       ! Restore pivot  
       rtemp1 = numbers(i)
       rtemp2 = dist(i)

       numbers(i) = numbers(right-1)
       dist(i)    = dist(right-1)

       numbers(right-1) = rtemp1
       dist(right-1)    = rtemp2
          
       call quick_sort_dp_dist(numbers, dist, left, i-1)
       call quick_sort_dp_dist(numbers, dist, i+1, right)
    end if
  end subroutine quick_sort_dp_dist

  subroutine median3_dp_dist(numbers, dist, left, right, pivot)
    implicit none

    integer(i4b)                  :: left, right
    real(dp)                      :: pivot
    real(dp),     dimension(:)    :: numbers, dist

    integer(i4b)                  :: center
    real(dp)                      :: rtemp1, rtemp2

    center = (left+right)/2

    if (dist(left)>dist(center)) then
       rtemp1 = numbers(left)
       rtemp2 = dist(left)

       numbers(left) = numbers(center)
       dist(left)    = dist(center)

       numbers(center) = rtemp1
       dist(center)    = rtemp2
    end if
    
    if (dist(left) > dist(right)) then
       rtemp1 = numbers(left)
       rtemp2 = dist(left)

       numbers(left) = numbers(right)
       dist(left)    = dist(right)

       numbers(right) = rtemp1
       dist(right)    = rtemp2
    end if

    if (dist(center) > dist(right)) then
       rtemp1 = numbers(center)
       rtemp2 = dist(center)

       numbers(center) = numbers(right)
       dist(center)    = dist(right)

       numbers(right) = rtemp1
       dist(right)    = rtemp2
    end if

    pivot = dist(center)

    ! Swap the pivot away
    rtemp1 = numbers(center)
    rtemp2 = dist(center)

    numbers(center) = numbers(right-1)
    dist(center)    = dist(right-1)

    numbers(right-1) = rtemp1
    dist(right-1)    = rtemp2

  end subroutine median3_dp_dist

  subroutine insertion_sort_dp_dist(numbers, dist)
    implicit none

    real(dp),     dimension(:)   :: numbers, dist

    integer(i4b)  :: length, i, j
    real(dp)      :: rtemp1, rtemp2

    length = size(numbers)

    do i = 2, length
       j = i

       rtemp1 = numbers(i)
       rtemp2 = dist(i)

       do while (rtemp2 < dist(j-1))
          dist(j) = dist(j-1)
          numbers(j) = numbers(j-1)
          j = j-1

          if (j == 1) then
             exit
          end if
       end do

       dist(j)    = rtemp2
       numbers(j) = rtemp1
    end do
  end subroutine insertion_sort_dp_dist



  !******************************************************************************

  ! This routine sorts an integer array
  subroutine QuickSort_int(numbers)
    implicit none

    integer(i4b), dimension(:)  :: numbers

    call quick_sort_int(numbers, 1, size(numbers))
    call insertion_sort_int(numbers)
  end subroutine QuickSort_int

  recursive subroutine quick_sort_int(numbers, left, right)
    implicit none
    
    integer(i4b), dimension(:)   :: numbers
    integer(i4b)                 :: left, right

    integer(i4b)  :: i,j, itemp, pivot

    if (left+10 < right) then

       call median3_int(numbers, left, right, pivot)

       i = left; j = right-1

       i = i+1
       do while (numbers(i)<pivot)
          i = i+1
       end do
       
       j = j-1
       do while (numbers(j)>pivot)
          j = j-1
       end do
       
       itemp = numbers(i)
       numbers(i) = numbers(j)
       numbers(j) = itemp

       do while (j>i)
          i = i+1
          do while (numbers(i)<pivot)
             i = i+1
          end do
    
          j = j-1
          do while (numbers(j)>pivot)
             j = j-1
          end do

          itemp = numbers(i)
          numbers(i) = numbers(j)
          numbers(j) = itemp
       end do

       ! Undo last swap
       itemp = numbers(i)
       numbers(i) = numbers(j)
       numbers(j) = itemp

       ! Restore pivot  
       itemp = numbers(i)
       numbers(i) = numbers(right-1)
       numbers(right-1) = itemp
          
       call quick_sort_int(numbers, left, i-1)
       call quick_sort_int(numbers, i+1, right)
    end if
  end subroutine quick_sort_int

  subroutine median3_int(numbers, left, right, pivot)
    implicit none

    integer(i4b), dimension(:)    :: numbers
    integer(i4b)                  :: left, right, pivot    

    integer(i4b)                  :: center, itemp

    center = (left+right)/2

    if (numbers(left) > numbers(center)) then
       itemp = numbers(left)
       numbers(left) = numbers(center)
       numbers(center) = itemp
    end if
    
    if (numbers(left) > numbers(right)) then
       itemp = numbers(left)
       numbers(left) = numbers(right)
       numbers(right) = itemp
    end if

    if (numbers(center) > numbers(right)) then
       itemp = numbers(center)
       numbers(center) = numbers(right)
       numbers(right) = itemp
    end if

    pivot = numbers(center)

    ! Swap the pivot away
    itemp = numbers(center)
    numbers(center) = numbers(right-1)
    numbers(right-1) = itemp
  end subroutine median3_int

  subroutine insertion_sort_int(numbers)
    implicit none

    integer(i4b), dimension(:)   :: numbers

    integer(i4b)  :: length, i, j
    integer(i4b)  :: itemp

    length = size(numbers)

    do i = 2, length
       j = i

       itemp = numbers(i)

       do while (itemp < numbers(j-1)) 
          numbers(j) = numbers(j-1)
          j = j-1

          if (j == 1) then
             exit
          end if
       end do

       numbers(j) = itemp
    end do
  end subroutine insertion_sort_int

! ****************************************************************

  ! This routine sorts a real array
  subroutine QuickSort_real(numbers)
    implicit none

    real(dp), dimension(:)  :: numbers

    call quick_sort_real(numbers, 1, size(numbers))
    call insertion_sort_real(numbers)
  end subroutine QuickSort_real

  recursive subroutine quick_sort_real(numbers, left, right)
    implicit none
    
    real(dp), dimension(:)   :: numbers
    integer(i4b)             :: left, right

    integer(i4b)  :: i,j
    real(dp)      :: rtemp, pivot

    if (left+10 < right) then

       call median3_real(numbers, left, right, pivot)

       i = left; j = right-1

       i = i+1
       do while (numbers(i)<pivot)
          i = i+1
       end do
       
       j = j-1
       do while (numbers(j)>pivot)
          j = j-1
       end do
       
       rtemp = numbers(i)
       numbers(i) = numbers(j)
       numbers(j) = rtemp

       do while (j>i)
          i = i+1
          do while (numbers(i)<pivot)
             i = i+1
          end do
    
          j = j-1
          do while (numbers(j)>pivot)
             j = j-1
          end do

          rtemp = numbers(i)
          numbers(i) = numbers(j)
          numbers(j) = rtemp
       end do

       ! Undo last swap
       rtemp = numbers(i)
       numbers(i) = numbers(j)
       numbers(j) = rtemp

       ! Restore pivot  
       rtemp = numbers(i)
       numbers(i) = numbers(right-1)
       numbers(right-1) = rtemp
          
       call quick_sort_real(numbers, left, i-1)
       call quick_sort_real(numbers, i+1, right)
    end if
  end subroutine quick_sort_real

  subroutine median3_real(numbers, left, right, pivot)
    implicit none

    real(dp), dimension(:)    :: numbers
    integer(i4b)              :: left, right
    real(dp)                  :: pivot    

    integer(i4b)              :: center
    real(dp)                  :: rtemp

    center = (left+right)/2

    if (numbers(left) > numbers(center)) then
       rtemp = numbers(left)
       numbers(left) = numbers(center)
       numbers(center) = rtemp
    end if
    
    if (numbers(left) > numbers(right)) then
       rtemp = numbers(left)
       numbers(left) = numbers(right)
       numbers(right) = rtemp
    end if

    if (numbers(center) > numbers(right)) then
       rtemp = numbers(center)
       numbers(center) = numbers(right)
       numbers(right) = rtemp
    end if

    pivot = numbers(center)

    ! Swap the pivot away
    rtemp = numbers(center)
    numbers(center) = numbers(right-1)
    numbers(right-1) = rtemp
  end subroutine median3_real

  subroutine insertion_sort_real(numbers)
    implicit none

    real(dp), dimension(:)   :: numbers

    integer(i4b)  :: length, i, j
    real(dp)      :: rtemp

    length = size(numbers)

    do i = 2, length
       j = i

       rtemp = numbers(i)

       do while (rtemp < numbers(j-1)) 
          numbers(j) = numbers(j-1)
          j = j-1

          if (j == 1) then
             exit
          end if
       end do

       numbers(j) = rtemp
    end do
  end subroutine insertion_sort_real

! ********************************************************

  subroutine bucket_sort(numbers, bindist, pos, sort_temp)
    implicit none

    integer(i4b), dimension(:),   intent(in)     :: bindist
    integer(i4b), dimension(:),   intent(inout)  :: numbers
    integer(i4b), dimension(:,:), intent(out)    :: pos, sort_temp

    integer(i4b)   :: i, j, current, temp_i
    integer(i4b), allocatable, dimension(:)   :: num

    allocate(num(size(pos(:,1))))
    num = 0

    ! Go through all points
    do i = 1, size(numbers)
       if (bindist(i) > -1) then
          temp_i = bindist(i)+1
          num(temp_i) = num(temp_i) + 1
          sort_temp(num(temp_i), temp_i) = numbers(i)
       end if
    end do

    ! Fill in the output arrays
    current = 1
    do i = 1, size(pos(:,1))
       if (num(i) > 0) then
          numbers(current:current+num(i)-1) = sort_temp(1:num(i),i)
          pos(i,1) = current-1
          pos(i,2) = current+num(i)-2

          current = current + num(i)
       else
          pos(i,1) = -1
          pos(i,2) = -1
       end if
    end do

    deallocate(num)

  end subroutine bucket_sort

! **********************************************************************

!!$  subroutine get_string_from_int(text, number, length)
!!$    implicit none
!!$
!!$    integer(i4b)  :: number, length
!!$    character(len=length) :: text, temp_text
!!$
!!$    if (length == 4) then
!!$       write(temp_text,'(I4)') number
!!$       temp_text = adjustl(temp_text); temp_text = trim(temp_text)
!!$       text = repeat('0',4-len_trim(temp_text)) // temp_text
!!$    else if (length == 5) then
!!$       write(temp_text,'(I5)') number
!!$       temp_text = adjustl(temp_text); temp_text = trim(temp_text)
!!$       text = repeat('0',5-len_trim(temp_text)) // temp_text
!!$    else if (length == 9) then
!!$       write(temp_text,'(I9)') number
!!$       temp_text = adjustl(temp_text); temp_text = trim(temp_text)
!!$       text = repeat('0',9-len_trim(temp_text)) // temp_text
!!$    end if
!!$  end subroutine get_string_from_int


end module sort_utils
