	integer * 4 function ftb_get_rse ( rse_file, rse )

c------------------------------------------------------------------------------
c
c	Function ftb_get_rse
c
c	This function reads an RSE file created by utility FTB_MAKE_RSE.
c	The file name is either specified by the user on the command line
c	or defaulted to Rse.Dat.
c	
c	This function is written in response to SPR 2537. It is also a derived
c	requirement from the specifications in the Brown Book as restated in the
c	following SPR description:
c
c	Version 4.1.1 10/25/88, SPR 2537, Shirley M. Read, STX
c		Teach FTB_Ctviewer to run in batch. All FIRAS facilities are
c	 	required to run in batch mode. The interactive/batch option,
c		a dataset qualifier and an Rse file qualifier must be added 
c		to the command line. Code switches nust be implemented for
c		the two operational modes.
c
c	Author: Shirley M. Read
c		ST Systems Corporation
c		October 25, 1988
c
c	Calling sequence: status = ftb_get_rse ( rse_file, rse )
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

	external	ftb_rseopen
	external	ftb_rseread
	external	ftb_rseclose
	external	ftb_normal

	ftb_get_rse = %loc(ftb_normal)

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
		 ftb_get_rse = %loc(ftb_rseread)
	         call lib$signal(ftb_rseread,%val(1),%val(r_status))
	      end if

	      close(nr,iostat=r_status)
	      if (r_status .ne. 0) then
		 ftb_get_rse = %loc(ftb_rseclose)
	         call lib$signal(ftb_rseclose,%val(1),%val(r_status))
	      end if

	      r_status = lib$free_lun(nr)

	      if (r_status .ne. ss$_normal) then
                 call lib$signal(r_status)
	      end if

	   else
	      ftb_get_rse = %loc(ftb_rseopen)
              call lib$signal(ftb_rseopen,%val(1),%val(r_status))
           end if

	else
	   ftb_get_rse = r_status
           call lib$signal(r_status)
	end if

	return
	end
