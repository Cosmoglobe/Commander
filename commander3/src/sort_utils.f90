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
module sort_utils
  use healpix_types
  implicit none

INTERFACE swap
  MODULE PROCEDURE swap_i, swap_r, swap_rv, swap_c,&
    swap_cv, swap_cm, swap_z, swap_zv, swap_zm,&
    masked_swap_is,masked_swap_iv,masked_swap_im,&
    masked_swap_rs,masked_swap_rv,masked_swap_rm
END INTERFACE


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



  ! This routine sorts an integer array, according to distances given
  ! by the additional array dist
  subroutine QuickSort_int_dist(numbers, dist)
    implicit none

    integer(i4b), dimension(:)  :: numbers
    integer(i4b), dimension(:)  :: dist

    call quick_sort_int_dist(numbers, dist, 1, size(numbers))
    call insertion_sort_int_dist(numbers, dist)
  end subroutine QuickSort_int_dist

  recursive subroutine quick_sort_int_dist(numbers, dist, left, right)
    implicit none
    
    integer(i4b), dimension(:)   :: numbers
    integer(i4b), dimension(:)   :: dist
    integer(i4b)                 :: left, right

    integer(i4b)  :: i,j, itemp
    integer(i4b)  :: pivot, rtemp

    if (left+10 < right) then

       call median3_int_dist(numbers, dist, left, right, pivot)

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
          
       call quick_sort_int_dist(numbers, dist, left, i-1)
       call quick_sort_int_dist(numbers, dist, i+1, right)
    end if
  end subroutine quick_sort_int_dist

  subroutine median3_int_dist(numbers, dist, left, right, pivot)
    implicit none

    integer(i4b)                  :: left, right
    integer(i4b)                  :: pivot
    integer(i4b), dimension(:)    :: numbers
    integer(i4b), dimension(:)    :: dist

    integer(i4b)                  :: center, itemp
    integer(i4b)                  :: rtemp

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
  end subroutine median3_int_dist

  subroutine insertion_sort_int_dist(numbers, dist)
    implicit none

    integer(i4b), dimension(:)   :: numbers
    integer(i4b), dimension(:)   :: dist

    integer(i4b)  :: length, i, j
    integer(i4b)  :: itemp
    integer(i4b)  :: rtemp

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
  end subroutine insertion_sort_int_dist


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

    call quick_sort_int(numbers, 1, size(numbers), 0)
    call insertion_sort_int(numbers)
  end subroutine QuickSort_int

  recursive subroutine quick_sort_int(numbers, left, right, depth)
    implicit none
    
    integer(i4b), dimension(:)   :: numbers
    integer(i4b)                 :: left, right, depth

    integer(i4b), parameter      :: max_depth = 1000
    integer(i4b)  :: i,j, itemp, pivot

    if (depth > max_depth) then
      ! Returns after max_depth recursive calls
      return
    end if


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
          
       call quick_sort_int(numbers, left, i-1, depth+1)
       call quick_sort_int(numbers, i+1, right, depth+1)
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

    integer(i4b)   :: i, current, temp_i
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


  ! Quicksort implementation from Numerical Recipes Fortran 90 implementation. 
  ! Gets around the very large stack for large integer arrays.
  SUBROUTINE quicksort_nr_sp(arr)
    IMPLICIT NONE
  
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
  
    INTEGER(I4B), PARAMETER :: NN = 15        ! Cutoff for insertion sort
    INTEGER(I4B), PARAMETER :: NSTACK = 50    ! Size of manual stack for subarrays
  
    REAL(SP) :: a                            ! Temporary storage for insertion sort and pivot
    INTEGER(I4B) :: n, k                    ! n = size of array, k = middle index
    INTEGER(I4B) :: i, j                    ! Loop variables and partition pointers
    INTEGER(I4B) :: jstack                  ! Stack pointer
    INTEGER(I4B) :: l, r                    ! Left and right bounds of current subarray
    INTEGER(I4B), DIMENSION(NSTACK) :: istack  ! Manual stack to hold subarray bounds
  
    n = size(arr)                           ! Array size
    jstack = 0                             ! Initialize stack pointer
    l = 1                                 ! Start index
    r = n                                 ! End index
  
    DO
      ! If subarray is small enough, use insertion sort
      IF (r - l < NN) THEN
        DO j = l + 1, r
          a = arr(j)
          DO i = j - 1, l, -1
            IF (arr(i) <= a) EXIT
            arr(i + 1) = arr(i)         ! Shift element right
          END DO
          arr(i + 1) = a                ! Insert element
        END DO
  
        ! If stack is empty, sorting is done
        IF (jstack == 0) RETURN
  
        ! Pop the next subarray bounds from the stack
        r = istack(jstack)
        l = istack(jstack - 1)
        jstack = jstack - 2
  
      ELSE
        ! Median-of-three pivot selection and reordering
        k = (l + r) / 2
  
        CALL swap(arr(k), arr(l + 1))
        CALL swap(arr(l), arr(r), arr(l) > arr(r))
        CALL swap(arr(l + 1), arr(r), arr(l + 1) > arr(r))
        CALL swap(arr(l), arr(l + 1), arr(l) > arr(l + 1))
  
        ! Initialize partition pointers
        i = l + 1
        j = r
        a = arr(l + 1)   ! Pivot element
  
        ! Partition loop
        DO
          ! Move i right until arr(i) >= pivot
          DO
            i = i + 1
            IF (arr(i) >= a) EXIT
          END DO
  
          ! Move j left until arr(j) <= pivot
          DO
            j = j - 1
            IF (arr(j) <= a) EXIT
          END DO
  
          IF (j < i) EXIT  ! Pointers crossed, partitioning done
  
          CALL swap(arr(i), arr(j))   ! Swap elements
        END DO
  
        ! Place pivot in correct position
        arr(l + 1) = arr(j)
        arr(j) = a
  
        ! Push larger subarray bounds onto stack; process smaller subarray immediately
        jstack = jstack + 2
        IF (jstack > NSTACK) write(*,*) 'quicksort_nr_sp: NSTACK too small'
  
        IF (r - i + 1 >= j - l) THEN
          istack(jstack) = r
          istack(jstack - 1) = i
          r = j - 1
        ELSE
          istack(jstack) = j - 1
          istack(jstack - 1) = l
          l = i
        END IF
      END IF
    END DO
  
  END SUBROUTINE quicksort_nr_sp
  
  SUBROUTINE quicksort_nr_int(arr)
      IMPLICIT NONE
  
      INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER(I4B), PARAMETER :: NN = 15, NSTACK = 50
  
      INTEGER(I4B) :: a
      INTEGER(I4B) :: n, k, i, j, jstack, l, r
      INTEGER(I4B), DIMENSION(NSTACK) :: istack
  
      n = SIZE(arr)
      jstack = 0
      l = 1
      r = n
  
      DO
          ! Insertion sort when subarray small enough
          IF (r - l < NN) THEN
              DO j = l + 1, r
                  a = arr(j)
                  DO i = j - 1, l, -1
                      IF (arr(i) <= a) EXIT
                      arr(i + 1) = arr(i)
                  END DO
                  arr(i + 1) = a
              END DO
  
              IF (jstack == 0) RETURN
              r = istack(jstack)
              l = istack(jstack - 1)
              jstack = jstack - 2
  
          ELSE
              ! Median-of-three partitioning
              k = (l + r) / 2
              CALL swap(arr(k), arr(l + 1))
              CALL swap(arr(l), arr(r), arr(l) > arr(r))
              CALL swap(arr(l + 1), arr(r), arr(l + 1) > arr(r))
              CALL swap(arr(l), arr(l + 1), arr(l) > arr(l + 1))
  
              i = l + 1
              j = r
              a = arr(l + 1)
  
              ! Partitioning loop
              DO
                  DO
                      i = i + 1
                      IF (arr(i) >= a) EXIT
                  END DO
                  DO
                      j = j - 1
                      IF (arr(j) <= a) EXIT
                  END DO
                  IF (j < i) EXIT
                  CALL swap(arr(i), arr(j))
              END DO
  
              arr(l + 1) = arr(j)
              arr(j) = a
  
              ! Push larger subarray on stack
              jstack = jstack + 2
              IF (jstack > NSTACK) write(*,*) "quicksort_nr_int: NSTACK too small"
  
              IF (r - i + 1 >= j - l) THEN
                  istack(jstack)     = r
                  istack(jstack - 1) = i
                  r = j - 1
              ELSE
                  istack(jstack)     = j - 1
                  istack(jstack - 1) = l
                  l = i
              END IF
          END IF
      END DO
  END SUBROUTINE quicksort_nr_int
  
  SUBROUTINE swap_i(a, b)
      INTEGER(I4B), INTENT(INOUT) :: a, b
      INTEGER(I4B) :: dum
      dum = a
      a   = b
      b   = dum
  END SUBROUTINE swap_i
  
  SUBROUTINE masked_swap_is(a, b, mask)
      INTEGER(I4B), INTENT(INOUT) :: a, b
      LOGICAL(LGT), INTENT(IN)    :: mask
      INTEGER(I4B) :: swp
      IF (mask) THEN
          swp = a
          a   = b
          b   = swp
      END IF
  END SUBROUTINE masked_swap_is
  
  SUBROUTINE masked_swap_iv(a, b, mask)
      INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: a, b
      LOGICAL(LGT),  DIMENSION(:), INTENT(IN)   :: mask
      INTEGER(I4B), DIMENSION(SIZE(a)) :: swp
      WHERE (mask)
          swp = a
          a   = b
          b   = swp
      END WHERE
  END SUBROUTINE masked_swap_iv
  
  SUBROUTINE masked_swap_im(a, b, mask)
      INTEGER(I4B), DIMENSION(:,:), INTENT(INOUT) :: a, b
      LOGICAL(LGT),  DIMENSION(:,:), INTENT(IN)   :: mask
      INTEGER(I4B), DIMENSION(SIZE(a,1), SIZE(a,2)) :: swp
      WHERE (mask)
          swp = a
          a   = b
          b   = swp
      END WHERE
  END SUBROUTINE masked_swap_im
  
  SUBROUTINE masked_swap_rs(a, b, mask)
    REAL(SP), INTENT(INOUT) :: a, b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp
  
    IF (mask) THEN
      swp = a
      a = b
      b = swp
    END IF
  END SUBROUTINE masked_swap_rs
  
  SUBROUTINE masked_swap_rv(a, b, mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a, b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(SIZE(a)) :: swp
  
    WHERE (mask)
      swp = a
      a = b
      b = swp
    END WHERE
  END SUBROUTINE masked_swap_rv
  
  SUBROUTINE masked_swap_rm(a, b, mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a, b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(SIZE(a,1), SIZE(a,2)) :: swp
  
    WHERE (mask)
      swp = a
      a = b
      b = swp
    END WHERE
  END SUBROUTINE masked_swap_rm
  
  SUBROUTINE swap_r(a, b)
    REAL(SP), INTENT(INOUT) :: a, b
    REAL(SP) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_r
  
  SUBROUTINE swap_rv(a, b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a, b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_rv
  
  SUBROUTINE swap_c(a, b)
    COMPLEX(SPC), INTENT(INOUT) :: a, b
    COMPLEX(SPC) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_c
  
  SUBROUTINE swap_cv(a, b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a, b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_cv
  
  SUBROUTINE swap_cm(a, b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a, b
    COMPLEX(SPC), DIMENSION(SIZE(a,1), SIZE(a,2)) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_cm
  
  SUBROUTINE swap_z(a, b)
    COMPLEX(DPC), INTENT(INOUT) :: a, b
    COMPLEX(DPC) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_z
  
  SUBROUTINE swap_zv(a, b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a, b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_zv

  SUBROUTINE swap_zm(a, b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a, b
    COMPLEX(DPC), DIMENSION(SIZE(a,1), SIZE(a,2)) :: dum
  
    dum = a
    a = b
    b = dum
  END SUBROUTINE swap_zm




end module sort_utils
