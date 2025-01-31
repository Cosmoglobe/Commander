
      Integer * 4 function fnt_power_density(chan_id,num,sci_data,spec_rec,
	1                        fake_it,dif_flag,number,size,lun,cindex,
	2                        config, new_segment, stat)

c--------------------------------------------------------------------------
c	Author: W. K. Young
c		STI Inc.
c		July 1986
c
c	Input:	channel #
c		number of ifgs
c		IFGS
c		fake-it bit
c
c	Output:	Coadded spectrum
c
c	Include files:
c		fnt_invoc.txt
c		fut_params.txt
c		ct$library:ctuser.inc
c
c	Calling Sequence:
c		status = fnt_power_density(channel #
c					   number of ifgs,
c					   ifgs,
c					   coadded spectrum,
c					   fake-it bit, dif_flag,
c                                          number, size, lun, cindex,
c                                          config, new_segment, stat)
c				
c	Gist of routine
c		IF number of ifgs greater than ZERO THEN
c		   DO for all ifgs
c		      FIND dither value
c		      SUBTRACT dither from data
c		      APODIZE data
c                     RESET 25 bins to zero
c		      FFT data
c		      SUM the square of the FFT
c		   END DO
c	 	   CALCULATE the power density
c		   CALL fnt_get_trnsfcn to get the appropriate electronix
c			 transfer function
c		   IF return status is an ERROR THEN
c		      SET status to ERROR
c		   ELSE
c	  	      REMOVE the effects of the electronics transfer function
c		   END IF
c		ELSE
c		   SET status to bad
c		END IF
c		RETURN
c
c--------------------------------------------------------------------------
c
c  Modifications:
c	19 Nov 1987 Isaacman added apodization function and put lower
c		    bound on usable amplitude of electronics xfer fcn
c	24 Dec 1987 Isaacman corrected IFG magnitude for number of MTM sweeps
c	 1 Feb 1988 F. Shuman made code consistent with coadd split
c	13 Mar 1988 R. Kummerer did IMSL convertion, SPR 1355.
c--------------------------------------------------------------------------
c *************************************************
c * Modification: 1988 08 15                      *
c * Purpose: Related to SPR 2158 .Mods made affect*
c * 2 primary items. First, the average time com- *****
c * putation is moved to just before FNT_get_transfcn *
c * CALL. Second, an argument with the average times  *
c * is added to the FNT_get_transfcn CALL.            *
c * Programmer: Nick Iascone, with lots of help from  *
c *             Bob Kummerer.                         *
c *****************************************************
c *****************************************************
c * Modification: 1988 10 08                          *
c * Purpose: To address SPR 2532. Noisetest processes *
c *          one channel at a time OK, but when       *
c *          /CHAN=ALL is used, it crashes.           *
c * Cause of problem: GAIN is not reset to 0 at the   ******
c *                   beginning of each FNT_Power_Density  *
c *                   call.                           ******
c * Programmer: Nick Iascone, STX/COBE                *
c *****************************************************
c *****************************************************
c * Modification: 1988 10 21                          *
c * Purpose: To address SPR 2658. The field           *
c *          .CHAN.FAKEIT is not being filled in the  *
c *          noise spectrum output record.            *
c * Solution: Since the FAKEIT value is already passed*
c *           to this function, the addition of the   *
c *           necessary assignment statements and as- *
c *           sociated Equivalence statement are the  *
c *           only requirement.                       *
c * Programmer: Nick Iascone, STX/COBE                *
c *****************************************************
c *****************************************************
c * Modification: 1988 11 01                          *
c * Purpose: To address SPR 2709. the Field           *
c *          .CHAN.CHAN_ID is not being filled in the *
c *          noise spectrum record.                   *
c * Solution: Since the CHAN_ID value is already      *
c *           passed to this function, the addition of*
c *           the necessary assignment statement(s) is*
c *           the only requirement.                   *
c * Programmer: Nick Iascone, STX/COBE                *
c *****************************************************
c ****************************************************************************
c                                                                            *
ch Change Log;                                                               *
ch                                                                           *
ch 	Version 4.2.1 12/05/88, Spr 2896, Gene Eplee, ARC.                   *
ch		Documentation in prolog for configuration: S. Read, STX.     *
ch    		FNT_Power_Density does not properly calculate the average    *
ch		number of sweeps to be written to the noise spectrum record. *
ch		The number of sweeps written to the record is always zero.   *
ch	        This happens because the variable 'sweeps' is reset before   *
ch              it is written to the record. This code is now corrected by   *
ch              adding 'avg_sweeps' and 'sum_sweeps' to calculate the average*
ch		number of sweeps and write it to the record.                 *
ch                                                                           *
c ****************************************************************************
c ********************************************************
c * Modification: 1989 01 13 (Friday the Thirteenth)     *
c * Purpose: To address SPR 2900. FNT needs to write     *
c *          several parameters to the noise spectrum    *
c *          records which are not currently written     *
c *          there.                                      *
c * Solution: Take the following steps:                  *
c *        ----- In FNT_power_density -----              *
c *           1. Assign Start_time & Stop_time, which are*
c *              already determined, to the proper       *
c *              Spec_rec variable.                      *
c *           2. Assign Num, which is already determined,*
c *              to the proper Spec_rec variable.        *
c *           3. Assign Bin_time_ave, which is already   *
c *              determined, to the proper Spec_rec vari-*
c *              able.                                   *
c *           4. Handle MTM length in a manner similar to*
c *              MTM speed. Then assign MTM_length to the*
c *              proper Spec_rec variable.               *
c *    ----- In FNT_power_density & FNT_read_noise ----- *
c *           5. Pass the commanded bolometer biases from*
c *              FNT_read_noise to FNT_power_density and *
c *              then assign them to the proper Spec_rec *
c *              variable.                               *
c * Programmer: Nick Iascone, STX/COBE                   *
c ********************************************************
c
c SER 3493, Use PLT graphics for IFG/spectra plots.
c August 18, 1989. R. Kummerer.
c *****************************************************
c * Modification: 1990 02 05                          *
c * Purpose: To take succesive differences of IFG's   *
c *	     to get noise spectra directly.  This     *
c *	     option is implemented with the "/diff"   *
c *	     switch in the command line.              *
c * Programmer: Joel Gales, ARC                       *
c *****************************************************
C   SER# 6625, 6626, FNT should call firas lib routine for nyquist freq.
C                    Fnt difference in spectra needed for noise monitoring
C              H. Wang, STX, 4/12/90
c--------------------------------------------------------------------------
c   Add statement to reset 25 bins to zero after ifg has been de-dithered 
c   and apodized.  SPR# 7590, Nov. 2, 1990, Nilo Gonzales/STX
c
c   SER 8919, S. Alexander, HSTX, August 3, 1992.  Electronics transfer
c             function changed to double precision.
c   SPR 9583,9790, Update FNT to use new FEX_Nyquist
c             Nilo G. Gonzales, Hughes STX, 1992 August 3.
c--------------------------------------------------------------------------
      IMPLICIT NONE

	include '(fut_params)'
	include '(fnt_invoc)'
	include 'ct$library:ctuser.inc'
	include '(cct_get_config)'

	integer * 4 	i		! a counter
	integer * 4	j		!another counter
	integer * 4	u_mode		!micro processor mode
	integer * 4	sweeps		!sweeps per IFG
	integer * 4	gain		!gain code
	integer * 4	ngroup		!adds per group
	integer * 4	mtm_speed	!mtm speed
	integer * 4	chan_id		!channel id
	integer * 4	num		!number of ifgs to be processed
	integer * 4	avg_group	!average of adds per group
c *  avg_sweeps added by Gene Eplee per SPR 2896.
        integer * 4	avg_sweeps	!average sweeps
	integer * 4	fake_it		!fake it stuff
        logical * 1     Dif_flag
c *************************************
c * Items related to mod 1988 10 21   *
c * SPR 2658                          *
c *                                   *
      Integer * 4 I4_fakeit
      Byte        Byte_fakeit(4)
      Equivalence(I4_fakeit, Byte_fakeit(1))
c *                                   *
c *************************************
c ***********************************************************
c * Items related to Mod 1989 01 13/SPR 2900.               *
c *                                                         *
c * For the purpose of receiving the value of Bol_cmd_biases*
c * from FNT_read_noise.                                    *
c *                                                         *
      Byte Bol_cmd_biases
      COMMON /BOL_BLOCK/ Bol_cmd_biases
      Integer * 2 I2_dummy
      Byte Byte_dummy(2)
      Equivalence (Byte_dummy(1), I2_dummy)
c *                                                         *
c * For the purpose of computing & assigning MTM_length.    *
c *                                                         *
      Integer * 2 MTM_length
      Integer * 4 Mtm_length4
c *                                                         *
c * For the purpose of Time conversion with ADT record ele- *
c * ment.                                                   *
c *                                                         *
      Integer * 4 Time_dummy(2)
c *                                                         *
c ***********************************************************
	integer * 4	r_status	!return status from get xferfcn
	
	real	* 4	apod		!apodization function
	real    * 4	dith		!dither level
	complex * 8	tifg(512)	!dedithered ifg
	complex *16	xferfcn(257)	!transfer function
	real	* 4	xampl		!amplitude of transfer fcn
	real	* 4	nyq_freq	!nyquist frequency
	real    * 4	sum_group	!summation of adds per group
c *  sum_sweeps added by Gene Eplee per SPR 2896
        real	* 4	sum_sweeps	!summation of sweeps

	complex	* 8	srcs	
	character * 14	start_time	!starting time of interval
	character * 14	stop_time	!ending time of interval
	character *  3	ifgno		!label for number of ifgs

	integer * 4     scan_mode
	integer * 4     index
	integer * 4     jtime(2)
	Integer * 4     number   	! number of data sets
	Integer * 4     size            ! size of data sets in bytes
	Integer * 4     lun(1)		! logical unit numbers
	Integer * 4     cindex(1)       ! initial cache pointers
	Logical * 1     new_segment(1)	! flag for new segments

	integer * 4	bin_time_array(2,fac_max_num)
	integer * 4	bin_ave_array(2,4)/0,0,0,0,0,0,0,0/
	integer * 4	bin_time_ave(2)
        Integer * 4 ave_num  ! Added per mode 1988 08 15 to
        Data ave_num /2/     ! Indicate number of ave times
	integer * 4	weights(fac_max_num)
	integer * 4	num_spec
	integer * 4	Iunits

	complex * 8	tmp_spec(512)	!temporary spectrum
	complex * 8	spec(257)	!a spectrum

	character * 8	fake_label	!label for fake it mode
	character * 10  preamp_label	!label for preamp
	character * 60	label		!noisetest generic label
	character * 14  string

	integer * 4	fnt_get_transfcn
	integer * 4	fut_ave_bin_times
c
c	items related to mod 1990 02 05 "succesive diff"
c
	logical * 2	diff_status
	complex * 8	diff_tifg(512)
	complex * 8	prev_tifg(512)

	integer*4	cli$present
	external	cli$_present
	external        fnt_getconfigerr
	integer * 4     cct_get_config_tod

c
c Records
c
	Dictionary 'fex_nyquist'
	Structure /config_data/
	   Record /fex_nyquist/fex_nyquist
	Endstructure

	Record /config_data/ config
	Record /config_status/ stat(1)

	dictionary 'NFS_SDF'
	record /NFS_SDF/ sci_data(num)	!science records
	dictionary 'FNT_NOISE'
	record /FNT_NOISE/ spec_rec
c------------------------------------------------------------------------
        	save bin_ave_array
c------------------------------------------------------------------------
c
        dif_flag = .false.
        MTM_length4 = 0 !Per SPR 2900 1989/01/13 DPI
        GAIN = 0       !Per SPR 2532 1988/10/08 DPI
c
	fnt_power_density = fac_normal
	do j=1,257
	   spec(j) = (0.0,0.0)		!initialize summing spectra
	end do
        num_spec = 0
	if(num .gt. 0)then
	   sum_group = 0.0		!initialize sum of adds per group to 0
c *  Initialization of sum_sweeps added by Gene Eplee per SPR 2896
           sum_sweeps = 0.0		! initialize sum of sweeps to 0
c
c	check for "\diff" switch in command line
c
	   diff_status = (cli$present('diff') .eq. %loc(cli$_present))	
c
	   do i=1,num
	      dith = 0.0
	      do j=1,512
	         dith = dith + sci_data(i).ifg_data.ifg(j)
	      end do
	      dith = dith/512.			!find dither value
	      do j=1,512
	         tifg(j) = sci_data(i).ifg_data.ifg(j)-dith	!remove dither
	      end do
c
c	Apodize and rescale the dedithered IFG. Rescaling is necessary to
c	conserve energy after apodization. Also, normalize for number
c	of sweeps in ifg.
c
	      sweeps = sci_data(i).sci_head.sc_head11
	      do j=1,512
		 apod = (1. - ((j - 256.5)/255.5)**2)**2
		 tifg(j) = tifg(j)/sweeps * apod * 1.57
	      enddo
c       
c             Reset 25 bins to zero after ifg has been de-dithered and
c             apodized.

              do j=1,25
	         tifg(j) = 0.0
	      enddo
c
c	if "diff" switch then take successive differences of ifg's
c
c
c		if 1st ifg then just store in prev_tifg else
c		take difference
c
c
	      if (diff_status) then
                Dif_flag = .true. 
		if (i.eq.1) then
		 do j=1,512
		 prev_tifg(j) = tifg(j)
	         end do
		else
		 do j=1,512
		 diff_tifg(j) = (tifg(j) - prev_tifg(j)) / sqrt(2.)
		 prev_tifg(j) = tifg(j)
		 tifg(j) = diff_tifg(j)
		 end do
		end if	
	      end if	
c
c	skip accumulation of ft if diff and i=1
c 
	if ((diff_status.eq..false.).or.(i.gt.1)) then	
c
c	Find the fft
c
	      call fftcf(512,tifg,tmp_spec)
c
c	Add things quadratically
c
	      do j=1,257
	         spec(j) = (tmp_spec(j) * conjg(tmp_spec(j))) +
     .						spec(j)
	      end do

	      ngroup = sci_data(i).sci_head.sc_head9
	      if (ngroup .eq. 0) then
	  	ngroup = 1
	      end if
	      sum_group = sum_group + ngroup
	      mtm_speed = mtm_speed + sci_data(i).sci_head.mtm_speed
       MTM_length4 =
     & MTM_length4 + sci_data(i).sci_head.MTM_length !Per SPR 2900
	      u_mode = u_mode + sci_data(i).sci_head.sc_head1a
	      gain = gain + sci_data(i).sci_head.gain
	      start_time = sci_data(i).ct_head.gmt
c *  Summation over sweeps modified by Gene Eplee per SPR 2896
              sum_sweeps = sum_sweeps + sci_data(i).sci_head.sc_head11
              bin_time_array(1,i) = sci_data(i).ct_head.time(1)
              bin_time_array(2,i) = sci_data(i).ct_head.time(2)
              weights(i) = 1
	      num_spec = num_spec + 1
c
	end if
c
	   end do
c
c	if diff then decrement # of spectra
c
	   if (diff_status) num = num - 1

	   avg_group = sum_group/num
	   mtm_speed = mtm_speed/num
      MTM_length = MTM_length4/num  !Per SPR 2900
	   u_mode = u_mode/num
c *  Calculation of avg_sweeps modified by Gene Eplee per SPR 2896
           avg_sweeps = sum_sweeps/num
	   gain = gain/num

	   Call ct_gmt_to_binary (start_time, jtime)
c
c Access FEX_Nquist reference dataset using CCT_Get_Config. 
c
	   r_status = cct_get_config_tod ( jtime, number, size,
	1                         lun, cindex, config, new_segment, stat )

	   if (.not. r_status) then
	      call lib$signal(fnt_getconfigerr,%val(1),%val(r_status))
	   endif

	   if (fake_it .eq. 1) then
	       if (avg_group .eq. 1) then
	          index = 9
	       else if (avg_group .eq. 2) then
	              index = 10
	       else if (avg_group .eq. 3) then
	              index = 11
	       else if (avg_group .eq. 8) then
	              index = 12
	       else 
	              index = 13
	       end if
	   else
	       if ((chan_id .eq. 1) .or. (chan_id .eq. 3)) then
                   index = 0
	       else
                   index = 4
	       end if
	   end if 
	   scan_mode = mtm_length * 2 + mtm_speed + 1
           index = index + scan_mode
	   nyq_freq = config.fex_nyquist.hz(index)

      r_status = fut_ave_bin_times(bin_time_array,
     .           num_spec,
     .           weights,bin_time_ave)

              spec_rec.ct_head.time(1) = bin_time_ave(1)
              spec_rec.ct_head.time(2) = bin_time_ave(2)
              call ct_binary_to_gmt(bin_time_ave(1), string)
	      spec_rec.ct_head.gmt = string
c -                                                          -
c ------------------------------------------------------------
c -   Arguments bin_time_ave & ave_num added per mod         -
c -   1988 08 15 .                                           -
c -                                                          -
	   r_status = fnt_get_transfcn(Bin_time_ave, ave_num,
     &                chan_id,mtm_speed,fake_it,
     &                u_mode,avg_group,xferfcn,ftc_preamp)
c -                                                          -
c ------------------------------------------------------------
	   if(r_status .eq. fac_normal)then

	      start_time = sci_data(1).ct_head.gmt
	      stop_time = sci_data(num).ct_head.gmt

	      spec_rec.chan.ngroup = avg_group
	      spec_rec.chan.mtm_speed = mtm_speed
	      spec_rec.chan.sci_mode = u_mode
	      spec_rec.chan.gain = gain
       I2_dummy = MTM_length                    !Per SPR 2900
       spec_rec.chan.mtm_length = Byte_dummy(1) !Per SPR 2900
       I2_dummy = 0                             !Per SPR 2900
       Byte_dummy(1) = Bol_cmd_biases           !Per SPR 2900
       spec_rec.chan.bol_cmd_bias = I2_dummy    !Per SPR 2900
       spec_rec.coa_head.first_gmt = start_time !Per SPR 2900
       spec_rec.coa_head.last_gmt =  stop_time  !Per SPR 2900
       CALL CT_GMT_to_binary(start_time, time_dummy)      !Per SPR 2900
       spec_rec.coa_head.first_bin_gmt(1) = time_dummy(1) !Per SPR 2900
       spec_rec.coa_head.first_bin_gmt(2) = time_dummy(2) !Per SPR 2900
       CALL CT_GMT_to_binary(stop_time, time_dummy)       !Per SPR 2900
       spec_rec.coa_head.last_bin_gmt(1)  = time_dummy(1) !Per SPR 2900
       spec_rec.coa_head.last_bin_gmt(2)  = time_dummy(2) !Per SPR 2900
       spec_rec.coa_head.num_ifgs = num                   !Per SPR 2900
       spec_rec.chan.av_collect_time(1) = bin_time_ave(1) !Per SPR 2900
       spec_rec.chan.av_collect_time(2) = bin_time_ave(2) !Per SPR 2900
c *  Write avg_sweeps to the record modified by Gene Eplee per SPR 2896.
              spec_rec.chan.sweeps = avg_sweeps
              I4_fakeit = fake_it                   !Per mod 1988 10 21/SPR 2658
              spec_rec.chan.fakeit = Byte_fakeit(1) !Per mod 1988 10 21/SPR 2658
              I4_fakeit = chan_id                   !Per mod 1988 11 01/SPR 2709
              spec_rec.chan.chan_id = Byte_fakeit(1)!Per mod 1988 11 01/SPR 2709

              if(fake_it .eq. 1)then
	            fake_label = 'Fk Mode'
	      else 
	            fake_label = 'MTM Mode'
	      end if
	      if(ftc_preamp .eq. fac_present)then
	            preamp_label = 'Preamp In'
	      else
	            preamp_label = 'Preamp Out'
  	      end if

	      write(ifgno,40)num
40	      format(i3)
	      label = start_time(1:11)//'_'//stop_time(1:11)//'_'//ifgno//
     .			' '//fake_label//' '//preamp_label

	      spec_rec.coa_head.label = label


c
c	start from two to avoid divide by zero. Also, don't use
c	values of the transfer function that are smaller than
c	~1% of the DC gain
c

	      do j=2,257
		 xampl = abs (xferfcn(j))
		 if (xampl.lt.3. .or. xampl.gt.1.e04) then
     		    xferfcn(j) = xferfcn(j-1)
		 endif

	         srcs =
     .			sqrt(spec(j)/(num*(conjg(xferfcn(j))*
     .			     xferfcn(j))*nyq_freq*512))

		 spec_rec.chan.spec(j) = srcs
	      end do
	   else
	      fnt_power_density = 13
	   end if
	else
	   fnt_power_density = 13
        endif 
	return
	end
