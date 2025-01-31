	integer*4 function FEC_SORT(
     M                              array,
     R                              length)
c
c      Written by Connie Lau, SASC Technologies Inc., May 1986
c
c      Commented by Reid Wilson, SASC Technologies Inc.,  23-AUG-1986
c
c      Description:
c            This function sorts a single dimension array with 'length'
c            elements, ordered smallest to largest.
c
c
c      Format:
c            Ret_Code = FEC_Sort (
c     M                           Array,
c     R                           Length)
c
c      Input Parameters:
c            Length                  = size of the array (number of elements)
c                                      integer*4
c      
c      Output Parameters:
c            Ret_Code                = status of operation
c                                      integer*4
c
c      Input/Output Parameters:
c            Array                   = array of values to be sorted
c                                      real*4 Array(*)
c
c      Modules Called:
c            None
c
c      Include Files:
c            FEC_INC                 = FEC parameters 
c            FEC_MSG                 = FEC message definition
c
c
c      *** BEGIN ***
c
        implicit none
c
        include '(FEC_INC)'
        include '(FEC_MSG)'
c
	real*4		array(*)            ! array containing values to sort
	real*4		temp                ! used for switching
c
	integer*4	length              ! number of elements in array
        integer*2       i                   ! loop counter
        integer*2       j                   ! temp variable
c
	logical*1	alldone/.FALSE./    ! flag for when done sorting
c
c
        do i = 1, length
           do j = i + 1, length
              if( array(j) .lt. array(i)) then
                 temp = array(i)
                 array(i) = array(j)
                 array(j) = temp
              end if
           end do
        end do
        FEC_Sort = %loc(FEC_NORMAL)
	RETURN
	END
