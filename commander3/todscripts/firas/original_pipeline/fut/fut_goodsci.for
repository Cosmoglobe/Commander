	logical*4 function fut_goodsci(sci_data, alltests, tests, flags)

c--------------------------------------------------------------------------
c	Author: RAShafer
c		685
c		1988 Oct 7
c	Function that assesses the quality of a science record using ONLY
c	information contained in the raw science.  If TRUE the science
c	record passed all the tests requested.
c
c_________________________________________________________________________

	implicit none
	Dictionary 'NFS_SDF'
	integer*4 numtests
	parameter (numtests=2)		! The number of tests that can
					!	be performed
c
c	Arguments:
c		Returned Value:
c	Logical*4 FUT_GOODSCI	: TRUE - The record passed the tests
c				: FALSE - The record failed at least one
c
	Record /NFS_SDF/ sci_data	! Input: science record
	logical*4 alltests		! Input: If TRUE all tests are
					!	performed
	logical*1 tests(numtests)	! Input: If tests(i) is TRUE, then
					!	the Ith test is performed,
					!	even if ALLTESTS is FALSE
	logical*1 flags(numtests)	! Returned: FLAGS(I) is set to TRUE
					!	if the Ith test is performed
					!	AND the science record passed
					!	the test, otherwise FALSE
c
c	Currently there is only two tests made.
c
c	Test 1:
c		That all 114 data ready bits for the IFG are set.
c	Test 2:
c		That the data checksum comes out correctly.
c
c	Future tests:
c		The Status bits are all zero (with exception of synch?)
c
c
	integer*4	i,j		! counters
	integer*2	ready(8)	! Correct data ready status
	data ready/7*'FFFF'x,'0003'x/

	integer*4 checksum	!   Function Call to calculate checksum
	integer*2 checkvalue	!   Returned checksum
	integer*4 status	!   Returned Code

	fut_goodsci = .true.

	if (alltests .or. tests(1)) then
c **	    ** data ready tests
	    flags(1) = .true.
	    i = 1
	    do while ( ready(i) .eq. sci_data.sci_head.data_ready(i)
     &	     .and. i.le.8)
		i = i+1
	      end do
	    flags(1) = i.ge. 9	! All of the data ready words match
	    fut_goodsci = flags(1)
	else
	    flags(1) = .false.
	  end if

	if (alltests .or. tests(2)) then
c **	    ** data checksum test
	    status = checksum(512, sci_data.ifg_data.ifg, checkvalue)
	    flags(2) = checkvalue .eq. sci_data.sci_head.sc_head6
	    fut_goodsci = fut_goodsci .and. flags(2)
	else
	    flags(2) = .false.
	  end if
	return
	end
