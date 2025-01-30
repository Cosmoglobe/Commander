	Program FRD_CALRES

C------------------------------------------------------------------------
C    PURPOSE:
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: H. WANG
C            ST Systems Corporation
C            Feb. 28, 1990
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
	integer * 4     i,j, k
	character * 14  in_file /'FEX_CALRES.TXT'/
	character * 14  out_file /'FEX_CALRES.DAT'/

	integer * 4	current_time(2)
	character * 14	gmt
        integer * 4	lib$get_lun
        integer * 4	lib$free_lun

	dictionary 'FEX_CALRES'
	record /FEX_CALRES/ calres

	external	FRD_Normal
	external	FRD_RMSOpen
	external	FRD_RMSRead
	external	FRD_RMSWrite
	external	FRD_RMSClose

c
c Retrieve the consistency check calres.
c
	status = lib$get_lun(lun)

	if (status .eq. ss$_normal) then

	   open (unit=lun, name=in_file, status='old',
     .			iostat=ios, readonly, shared)

	   if (ios .eq. 0) then

c
c Interpret the RMS calres file.
c
            Do j = 1,4
               read(lun,*) (Calres.kal_counts(j,i),i=1,4)
            Enddo

c
c Build the output record.
c

	      call sys$gettim(current_time)
	      call ct_binary_to_gmt(current_time,gmt)

	      calres.ct_head.gmt = gmt

	      do k=1,2
	         calres.ct_head.time(k) = current_time(k)
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
     .			organization='sequential', recordsize=64,
     .			recordtype='fixed', iostat=ios, shared)

	      if (ios .eq. 0) then

	         write (lun, iostat=ios) calres

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
