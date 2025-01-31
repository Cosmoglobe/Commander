	Program FRD_GTRANS

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    PURPOSE: This program creates three reference datasets (FEX_GRTTRANS.DAT,
c             FEX_GRTRAWWT.DAT, FEX_GRTCOAWT.DAT) which will be used by 
c             the FDQ facility.
c
c    Written by: Nilo G. Gonzales/STX, August 23, 1991
c
c    PDL:
c       Begin
c            If output qualifier is present, then
c               If qualifier is equal 'TRAN', then
c                  set tran flag to true
c               Else If qualifier is equal 'RAW', then
c                       set raw flag to true
c               Else If qualifier is equal 'COAD', then
c                       set coadd qualifier to true
c               End If
c            Else 
c                set default flag 'ALL' to true
c            End If
c
c            For each type of file
c               If flag is true, then
c                  get lun for input text file                   
c                  open input file
c                  read file
c                  build output record
c                  open output file
c                  write output record
c                  close lun
c               End if
c               stop
c	      End for
c       End 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	implicit none

	include '($ssdef)'
	include '(upm_stat_msg)'

c Functions:

	integer * 4     upm_present
	integer * 4     upm_get_value
	integer * 4     str$upcase
        integer * 4	lib$get_lun
        integer * 4	lib$free_lun

c Externals:

	external	FRD_Normal
	external        FRD_Abort
	external	FRD_RMSOpen
	external	FRD_RMSRead
	external	FRD_RMSWrite
	external	FRD_RMSClose

c Local Declarations:
 
	dictionary 'FEX_GRTTRANS'
	record /FEX_GRTTRANS/ gtran_rec

	dictionary 'FEX_GRTRAWWT'
	record /FEX_GRTRAWWT/ grawt_rec

	dictionary 'FEX_GRTCOAWT'
	record /FEX_GRTCOAWT/ gcoat_rec

	character * 16  in_file1 /'FEX_GRTTRANS.TXT'/
	character * 16  in_file2 /'FEX_GRTRAWWT.TXT'/
	character * 16  in_file3 /'FEX_GRTCOAWT.TXT'/
	character * 16  out_file1 /'FEX_GRTTRANS.DAT'/
	character * 16  out_file2 /'FEX_GRTRAWWT.DAT'/
	character * 16  out_file3 /'FEX_GRTCOAWT.DAT'/
	character * 6   output
	character * 4   ovalue
	integer   * 2   len
	integer   * 4   upm_stat
        integer   * 4   ios
	integer   * 4	status
	integer   * 4   lun_gtran, lun_rawt, lun_coat
	integer   * 4   k
	logical   * 1   trans /.false./, rawt /.false./
	logical   * 1   coad /.false./, all /.false./
	integer   * 4	current_time(2)
	character * 14	gmt

	real      * 4	grt_a_trans_temp(16)
	real      * 4	grt_a_trans_hwid(16)
        real      * 4   grt_a_weight(16)

	real      * 4	grt_b_trans_temp(16)
	real      * 4	grt_b_trans_hwid(16)
        real      * 4   grt_b_weight(16)


c Parse the command line.

	if (upm_present('output') .eq. upm_pres) then
           upm_stat = upm_get_value('output',ovalue,len)
	   if (upm_stat .ne. ss$_normal) then
	      call lib$signal (%val(ss$_abort))
	   end if
	   if (ovalue .eq. 'TRAN') then
	      trans = .true.
	   else if (ovalue .eq. 'RAW') then
	      rawt = .true.
	   else if (ovalue .eq. 'COAD') then
	      coad = .true.
	   end if
	else
	   all = .true.
	end if	
c
c Retrieve the GRT temperature transition.
c
	if (trans .or. all) then
	   status = lib$get_lun(lun_gtran)
	   if (status .eq. ss$_normal) then
	      open (unit=lun_gtran, name=in_file1, status='old',
     .	           iostat=ios, readonly, shared)
	      if (ios .eq. 0) then
c
c Interpret the RMS  file.
c
	         do k=1,16
	            read (lun_gtran,*) grt_a_trans_temp(k),
	1                 grt_a_trans_hwid(k)
	         end do

	         do k=1,16
	            read (lun_gtran,*), grt_b_trans_temp(k),
	1                 grt_b_trans_hwid(k)
	         end do
c
c Build the output record.
c
	         gtran_rec.num_grt = 16

	         do k=1,16
		    gtran_rec.grt_a_trans_temp(k) = grt_a_trans_temp(k)
		    gtran_rec.grt_a_trans_hwid(k) = grt_a_trans_hwid(k)
	         end do

	         do k=1,16
		    gtran_rec.grt_b_trans_temp(k) = grt_b_trans_temp(k)
		    gtran_rec.grt_b_trans_hwid(k) = grt_b_trans_hwid(k)
	         end do

	         call sys$gettim(current_time)
	         call ct_binary_to_gmt(current_time,gmt)

	         gtran_rec.ct_head.gmt = gmt

	         do k=1,2
	            gtran_rec.ct_head.time(k) = current_time(k)
	         end do

	         close (lun_gtran, iostat=ios)

	         if (ios .ne. 0) then
		    status = %loc(FRD_RMSClose)
	            call lib$signal(FRD_RMSClose,%val(2),in_file1,%val(ios))
	         end if
c
c Write the output record.
c
	         open (unit=lun_gtran, name=out_file1, status='new',
     .			form='unformatted', access='sequential',
     .			organization='sequential', recordsize=96,
     .			recordtype='fixed', iostat=ios, shared)

	         if (ios .eq. 0) then

	            write (lun_gtran, iostat=ios) gtran_rec

	            if (ios .ne. 0) then
		       status = %loc(FRD_RMSWrite)
	               call lib$signal(FRD_RMSWrite,%val(2),out_file1,%val(ios))
	  	    end if

	         else
		    status = %loc(FRD_RMSOpen)
	            call lib$signal(FRD_RMSOpen,%val(2),out_file1,%val(ios))
	         end if

	         close (lun_gtran, iostat=ios)

	         if (ios .ne. 0) then
		    status = %loc(FRD_RMSClose)
	            call lib$signal(FRD_RMSClose,%val(2),out_file1,%val(ios))
	         end if

	      else
	         status = %loc(FRD_RMSOpen)
	         call lib$signal(FRD_RMSOpen,%val(2),in_file1,%val(ios))
	      end if

	      status = lib$free_lun(lun_gtran)

	   else
	      call lib$signal(%val(status))
	   end if

	   if (status) then
	      call lib$signal(FRD_Normal)
	   else
	      call lib$signal(%val(ss$_abort))
	   end if
	end if     ! output qualifier = tran
c
c Retrieve weights for Raw data.
c
	if (rawt .or. all) then
	   status = lib$get_lun(lun_rawt)
	   if (status .eq. ss$_normal) then
	       open (unit=lun_rawt, name=in_file2, status='old',
     .	  	   iostat=ios, readonly, shared)
	       if (ios .eq. 0) then
c
c Interpret the RMS file.
c
	          do k=1,16
	           read (lun_rawt,*) grt_a_weight(k)
	          end do
	          do k=1,16
	             read (lun_rawt,*), grt_b_weight(k)
	          end do
c
c Build the output record.    
c
	          grawt_rec.num_grt = 16

	          do k=1,16
		     grawt_rec.grt_a_weight(k) = grt_a_weight(k)
	          end do

	          do k=1,16
		     grawt_rec.grt_b_weight(k) = grt_b_weight(k)
	          end do

	          call sys$gettim(current_time)
	          call ct_binary_to_gmt(current_time,gmt)

	          grawt_rec.ct_head.gmt = gmt

	          do k=1,2
	             grawt_rec.ct_head.time(k) = current_time(k)
	          end do

	          close (lun_rawt, iostat=ios)

	          if (ios .ne. 0) then
		     status = %loc(FRD_RMSClose)
	             call lib$signal(FRD_RMSClose,%val(2),in_file2,%val(ios))
	          end if
c
c Write the output record.
c
	          open (unit=lun_rawt, name=out_file2, status='new',
     .			form='unformatted', access='sequential',
     .			organization='sequential', recordsize=64,
     .			recordtype='fixed', iostat=ios, shared)

	          if (ios .eq. 0) then
	             write (lun_rawt, iostat=ios) grawt_rec
	             if (ios .ne. 0) then
		        status = %loc(FRD_RMSWrite)
	                call lib$signal(FRD_RMSWrite,%val(2),
	1                    out_file2,%val(ios))
	  	     end if
	          else
		     status = %loc(FRD_RMSOpen)
	             call lib$signal(FRD_RMSOpen,%val(2),out_file2,%val(ios))
	          end if

	          close (lun_rawt, iostat=ios)

	          if (ios .ne. 0) then
		     status = %loc(FRD_RMSClose)
	             call lib$signal(FRD_RMSClose,%val(2),out_file2,%val(ios))
	          end if
	       else
	          status = %loc(FRD_RMSOpen)
	          call lib$signal(FRD_RMSOpen,%val(2),in_file2,%val(ios))
	       end if

	       status = lib$free_lun(lun_rawt)
	    else
	       call lib$signal(%val(status))
	    end if

	    if (status) then
	       call lib$signal(FRD_Normal)
	    else
	       call lib$signal(%val(ss$_abort))
	    end if
	end if       ! output qualifier = raw
c
c Retrieve weights for Coadd data
c
	if (coad .or. all) then
	   status = lib$get_lun(lun_coat)
	   if (status .eq. ss$_normal) then
	      open (unit=lun_coat, name=in_file3, status='old',
     .	           iostat=ios, readonly, shared)

	      if (ios .eq. 0) then
c
c Interpret the RMS weights for coad file.
c
	         do k=1,16
	            read (lun_coat,*) grt_a_weight(k)
	         end do

	         do k=1,16
	            read (lun_coat,*), grt_b_weight(k)
	         end do
c
c Build the output record.    
c
	         gcoat_rec.num_grt = 16

	         do k=1,16
		    gcoat_rec.grt_a_weight(k) = grt_a_weight(k)
	         end do

	         do k=1,16
		    gcoat_rec.grt_b_weight(k) = grt_b_weight(k)
	         end do

	         call sys$gettim(current_time)
	         call ct_binary_to_gmt(current_time,gmt)

	         gcoat_rec.ct_head.gmt = gmt

	         do k=1,2
	            gcoat_rec.ct_head.time(k) = current_time(k)
	         end do

	         close (lun_coat, iostat=ios)

	         if (ios .ne. 0) then
		    status = %loc(FRD_RMSClose)
	            call lib$signal(FRD_RMSClose,%val(2),in_file3,%val(ios))
	         end if
c
c Write the output record.
c
	         open (unit=lun_coat, name=out_file3, status='new',
     .			form='unformatted', access='sequential',
     .			organization='sequential', recordsize=64,
     .			recordtype='fixed', iostat=ios, shared)

	         if (ios .eq. 0) then
	            write (lun_coat, iostat=ios) gcoat_rec
	            if (ios .ne. 0) then
		       status = %loc(FRD_RMSWrite)
	               call lib$signal(FRD_RMSWrite,%val(2),out_file3,%val(ios))
	  	    end if

	         else
		    status = %loc(FRD_RMSOpen)
	            call lib$signal(FRD_RMSOpen,%val(2),out_file3,%val(ios))
	         end if

	         close (lun_coat, iostat=ios)

	         if (ios .ne. 0) then
		    status = %loc(FRD_RMSClose)
	            call lib$signal(FRD_RMSClose,%val(2),out_file3,%val(ios))
	         end if
	      else
	         status = %loc(FRD_RMSOpen)
	         call lib$signal(FRD_RMSOpen,%val(2),in_file3,%val(ios))
	      end if

	      status = lib$free_lun(lun_coat)

	   else
	      call lib$signal(%val(status))
	   end if

	   if (status) then
	      call lib$signal(FRD_Normal)
	   else
	      call lib$signal(%val(ss$_abort))
	   end if
	end if     ! output qualifier = coad
	stop
	end
