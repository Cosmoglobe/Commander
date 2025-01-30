	Program FRD_MINCOADD

C------------------------------------------------------------------------
C    PURPOSE:
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: H. Wang
C            Hughes STX Corporation
C            Dec. 4, 1991
C            SER 9338
C
C    INVOCATION: 
C
C    INPUT PARAMETERS:
C
C    OUTPUT PARAMETERS: 
C
C    SUBROUTINES CALLED: 
C
C    COMMON VARIABLES USED:
C
C    INCLUDE FILES: 
C
C    PROCESSING METHOD:
C  
C----------------------------------------------------------------------
C
C Changes:
C
C
C----------------------------------------------------------------------

	implicit none

	include '($ssdef)'

        integer * 4     ios
	integer * 4	status
	integer * 4     lun
	integer * 4     k
	character * 16  in_file /'FEX_MINCOADD.TXT'/
	character * 16  out_file /'FEX_MINCOADD.DAT'/

	integer * 4	current_time(2)
	character * 14	gmt

        integer * 4	min_number(4)

        integer * 4	lib$get_lun
        integer * 4	lib$free_lun

	dictionary 'FEX_MINCOADD'
	record /FEX_MINCOADD/ MINCOADD

	external	FRD_Normal
	external	FRD_RMSOpen
	external	FRD_RMSRead
	external	FRD_RMSWrite
	external	FRD_RMSClose

c
c Retrieve the consistency check mincoadd.
c
	status = lib$get_lun(lun)

	if (status .eq. ss$_normal) then

	   open (unit=lun, name=in_file, status='old',
     .			iostat=ios, readonly, shared)

	   if (ios .eq. 0) then

c
c Interpret the RMS mincoadd file.
c
	      read (lun,*) (min_number(k),k=1,4)
c
c Build the output record.
c
	      do k=1,4
		 mincoadd.min_ifg_coadd(k) = min_number(k)
	      end do
	      call sys$gettim(current_time)
	      call ct_binary_to_gmt(current_time,gmt)

	      mincoadd.ct_head.gmt = gmt

	      do k=1,2
	         mincoadd.ct_head.time(k) = current_time(k)
	      end do

	      close (lun, iostat=ios)

	      if (ios .ne. 0) then
		 status = %loc(FRD_RMSClose)
	         call lib$signal(FRD_RMSClose,%val(2),in_file,%val(ios))
	      end if

c
c Write the output record.
c
	      open (unit=lun, name=out_file, status='new',
     .			form='unformatted', access='sequential',
     .			organization='sequential', recordsize=32,
     .			recordtype='fixed', iostat=ios, shared)

	      if (ios .eq. 0) then

	         write (lun, iostat=ios) mincoadd

	         if (ios .ne. 0) then
		    status = %loc(FRD_RMSWrite)
	            call lib$signal(FRD_RMSWrite,%val(2),out_file,%val(ios))
	  	 end if

	      else
		 status = %loc(FRD_RMSOpen)
	         call lib$signal(FRD_RMSOpen,%val(2),out_file,%val(ios))
	      end if

	      close (lun, iostat=ios)

	      if (ios .ne. 0) then
		 status = %loc(FRD_RMSClose)
	         call lib$signal(FRD_RMSClose,%val(2),out_file,%val(ios))
	      end if

	   else
	     status = %loc(FRD_RMSOpen)
	     call lib$signal(FRD_RMSOpen,%val(2),in_file,%val(ios))
	   end if

	   status = lib$free_lun(lun)

	else
	   call lib$signal(%val(status))
	end if

	if (status) then
	   call lib$signal(FRD_Normal)
	else
	   call lib$signal(%val(ss$_abort))
	end if

	stop
	end
