	integer * 4 function fut_get_rse ( rse_file, rse )

c------------------------------------------------------------------------------
c
c	Function fut_get_rse
c
c	This function reads an RSE file; the file name is either specified 
c       by the user on the command line or defaulted to Rse.Dat.
c	
c	Author: Shirley M. Read
c		ST Systems Corporation
c		October 25, 1988
c
c	Calling sequence: status = fut_get_rse ( rse_file, rse )
c
c	Input:
c		Rse_File -- file containing the Rse
c
c	Output:
c		Rse -- record selection expression
c
c	Subroutines called:
c		none
c
c	Include files:
c		$ssdef
c
c       Modifications:
c   
c       Moved to FUT; identical to old FTB_GET_RSE subroutine.  
c       SPR 9430; N. Gonzales/Hughes STX, January 23, 1992.
c
c------------------------------------------------------------------------------

	implicit none

c	Include Files.

	include '($ssdef)'

c	Passed Parameters.

	character * 64  rse_file
	character * 128 rse(16)

c	Local Declarations.

	integer * 4	nr
	integer * 4	ix
	integer * 4	r_status

	integer	* 4	lib$get_lun
	integer * 4	lib$free_lun

	external	fut_rseopen
	external	fut_rseread
	external	fut_rseclose
	external	fut_normal

	fut_get_rse = %loc(fut_normal)

	do ix = 1,16
	   rse(ix) = ' '
	end do

        r_status = lib$get_lun(nr)

	if (r_status .eq. ss$_normal) then

	   open(unit=nr,file=rse_file,status='old',
     .			iostat=r_status,shared,readonly)

	   if (r_status .eq. 0) then

	      read (nr,10,iostat=r_status) (rse(ix),ix=1,16)
10	      format (a128)
	      if (r_status .ne. 0) then
		 fut_get_rse = %loc(fut_rseread)
	         call lib$signal(fut_rseread,%val(1),%val(r_status))
	      end if

	      close(nr,iostat=r_status)
	      if (r_status .ne. 0) then
		 fut_get_rse = %loc(fut_rseclose)
	         call lib$signal(fut_rseclose,%val(1),%val(r_status))
	      end if

	      r_status = lib$free_lun(nr)

	      if (r_status .ne. ss$_normal) then
                 call lib$signal(r_status)
	      end if

	   else
	      fut_get_rse = %loc(fut_rseopen)
              call lib$signal(fut_rseopen,%val(1),%val(r_status))
           end if

	else
	   fut_get_rse = r_status
           call lib$signal(r_status)
	end if

	return
	end
