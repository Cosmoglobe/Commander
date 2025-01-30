	Program FRD_CC_Thresholds

C------------------------------------------------------------------------
C    PURPOSE: Put data from FEX_CTH.TXT reference text file into FEX_CTH
C             record.
C
C       INPUT DATA: FEX_CTH.TXT
C
C       OUTPUT DATA: FEX_CTH.DAT
C
C    AUTHOR: R. Kummerer
C            ST Systems Corporation
C            October 14, 1988
C
C    INCLUDE FILES: $ssdef
C
C----------------------------------------------------------------------
C
C Changes:
C
C	SPR 3736, Use input text file FEX_CTH.TXT. R. Kummerer,
C		May 1, 1989.
C
C       Updated for new requirement for FIC. SER 7985.
C              Hughes STX, H. Wang, Dec 4, 1991
C       
C       Add galactic latitude cutoffs for neighbor selection and
C              dihedral temperature bin. S. Brodd, 8/28/95, SER 12244.
C----------------------------------------------------------------------

	implicit none

	include '($ssdef)'

        integer * 4	lib$get_lun
        integer * 4	lib$free_lun

        integer * 4     ios, i
	integer * 4	status
	integer * 4     lun
	integer * 4     k
	character * 11  in_file /'FEX_CTH.TXT'/
	character * 11  out_file /'FEX_CTH.DAT'/

	integer * 4	current_time(2)
	character * 14	gmt

        Real *4         Bol_VOL_TOL(4)
        Real *4         GRT_Tol(10)
        Real *4         Spa_Grad(5,3)
        Integer *2      Temp_Grad(20)
        Real *4         Gal_Lats(4)
        Real *4         Up_Bounds(4)
        Real *4         Lo_Bounds(4)
        Integer *2      Peak_Pts(4)
        real * 4        prim_temp_amp(4,4)
        real * 4        prim_temp_snr(4,4)
        real * 4        sec_temp_amp(4,4)
        real * 4        sec_temp_snr(4,4)
        real * 4	noise_min(4)
	real * 4	noise_max(4)
	Integer * 4	mask_ifg(2,4)
        real * 4	Max_pt_deviation(4)
	integer * 2	Max_bad_pts(4)
        integer * 2	Min_IFG_shp_number(4)

        real * 4        temp_Xcal
        real * 4        temp_Ical
        real * 4        avg_SkyREF_temp
        real * 4        temp_Dihed
        real * 4	Bol_temp(4)
        real * 4        Gli_Rat(4)
        integer *2      Abs_Gal_Lat 
        Character*14    GMT_SAA
        Integer*4       SAA_Time(2)

	dictionary 'FEX_CTH'
	record /FEX_CTH/ thresholds

	external	FRD_Normal
	external	FRD_RMSOpen
	external	FRD_RMSRead
	external	FRD_RMSWrite
	external	FRD_RMSClose
c
c Retrieve the consistency check thresholds.
c
	status = lib$get_lun(lun)

	if (status .eq. ss$_normal) then

	   open (unit=lun, name=in_file, status='old',
     .			iostat=ios, readonly, shared)

	   if (ios .eq. 0) then
c
c Interpret the RMS thresholds file.
c
	      read (lun,*) (Bol_Vol_Tol(k),k=1,4)
	      read (lun,*) (Grt_Tol(k),k=1,10)
              Do I =1, 5
	        read (lun,*) (Spa_Grad(I,k),k=1,3)
              Enddo
	      read (lun,*) (Temp_Grad(k),k=1,10)
	      read (lun,*) (Temp_Grad(k),k=11,20)
	      read (lun,*) (Gal_Lats(k),k=1,4)
	      read (lun,*) (Up_Bounds(k),k=1,4)
	      read (lun,*) (Lo_Bounds(k),k=1,4)
	      read (lun,*) (Peak_pts(k),k=1,4)
              Do I =1, 4
	        read (lun,*) (prim_temp_amp(I,k),k=1,4)
              Enddo
              Do I =1, 4
	        read (lun,*) (prim_temp_snr(I,k),k=1,4)
              Enddo
              Do I =1, 4
	        read (lun,*) (sec_temp_amp(I,k),k=1,4)
              Enddo
              Do I =1, 4
	        read (lun,*) (sec_temp_snr(I,k),k=1,4)
              Enddo
	      read (lun,*) (noise_min(k),k=1,4)
	      read (lun,*) (noise_max(k),k=1,4)
              Do i =1, 2
	        read (lun,*) (mask_ifg(i,k),k=1,4)
              Enddo
	      read (lun,*) (Max_pt_deviation(k),k=1,4)
	      read (lun,*) (Max_bad_pts(k),k=1,4)
	      read (lun,*) (Min_IFG_shp_number(k),k=1,4)
 	      read (lun,*) (Temp_Xcal)
 	      read (lun,*) (Temp_Ical)
	      read (lun,*) (Avg_SkyREF_temp)
 	      read (lun,*) (Temp_Dihed)
	      read (lun,*) (Bol_temp(k),k=1,4)
 	      read (lun,*) (Gli_rat(k),k=1,4)
	      read (lun,*) (Abs_Gal_Lat)
	      read (lun,14) (Gmt_SAA)
14	      format(a14)
	      call ct_gmt_to_binary(Gmt_SAA, SAA_TIME)
c
c Build the output record.
c
	      do k=1,4
		 thresholds.Bolometer_Voltage_tolerances(k) = Bol_Vol_tol(k)
                 thresholds.gal_lat_cutoff(k) = gal_lats(k)
		 thresholds.upper_bounds(k) = up_bounds(k)
		 thresholds.Lower_bounds(k) = Lo_bounds(k)
		 thresholds.peak_points(k) = peak_pts(k)
		 thresholds.min_ifg_noise(k) = noise_min(k)
		 thresholds.max_ifg_noise(k) = noise_max(k)
		 thresholds.max_point_deviation(k) = max_pt_deviation(k)
		 thresholds.max_bad_points(k) = Max_bad_pts(k)
		 thresholds.min_IFG_SHAPE_COADD(k) = min_IFG_SHP_NUMBER(k)
		 thresholds.Bolometer_Temp(k) = Bol_temp(k)
		 thresholds.Glitch_Rate(k) = Gli_Rat(k)
	      end do

	      do k=1,10
		 thresholds.grt_tolerances(k) = Grt_Tol(k)
	      end do
              Do i=1, 5
                Do  k = 1, 3
		 thresholds.Spatial_Gradients(i,k) = Spa_Grad(i,k)
                Enddo
              enddo
	      do k=1,20
		 thresholds.Temporal_gradients(k) = Temp_grad(k)
	      end do
              Do i=1, 4
                Do  k = 1, 4
		 thresholds.prim_temp_amp(i,k) = prim_temp_amp(i,k)
                Enddo
              enddo
              Do i=1, 4
                Do  k = 1, 4
		 thresholds.prim_temp_snr(i,k) = prim_temp_snr(i,k)
                Enddo
              enddo
              Do i=1, 4
                Do  k = 1, 4
		 thresholds.sec_temp_amp(i,k) = sec_temp_amp(i,k)
                Enddo
              enddo
              Do i=1, 4
                Do  k = 1, 4
		 thresholds.sec_temp_snr(i,k) = sec_temp_snr(i,k)
                Enddo
              enddo
              Do i=1, 2
                Do  k = 1, 4
		 thresholds.Mask_IFG(i,k) = mask_ifg(i,k)
                Enddo
              enddo
              thresholds.xcal_temp = temp_xcal
              thresholds.ical_temp = temp_ical
              thresholds.Skyhorn_Refhorn_avg = avg_skyref_temp 
              thresholds.dihedral_temp = temp_dihed
              thresholds.abs_galactic_lat = Abs_Gal_Lat
	      thresholds.GMT_TIME_SAA = Gmt_SAA

	      call sys$gettim(current_time)
	      call ct_binary_to_gmt(current_time,gmt)

	      thresholds.ct_head.gmt = gmt

	      do k=1,2
		 thresholds.TIME_SAA(k) = SAA_TIME(k)
	         thresholds.ct_head.time(k) = current_time(k)
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
     .			organization='sequential', recordsize=256,
     .			recordtype='fixed', iostat=ios, shared)

	      if (ios .eq. 0) then

	         write (lun, iostat=ios) thresholds

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
