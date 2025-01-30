      program fnt

c--------------------------------------------------------------------------
c
c	This program does the noise analysis for FIRAS
c
c	Author: W. K. Young
c		STI Inc.
c		July 1986
c
c	Input:
c		None
c	Output:
c		None
c
c	Include files:
c		fnt_noise_invoc.inc
c	Calling Sequence
c		Main program
c	Gist of routine
c		CALL fnt_parse_noise to parse the command line
c		IF status equal success THEN
c                  Call Ct_Init
c               End If
c               If First_time Then
c                  Call CCT_Open_Config for reference file FEX_Nyquist
c               End If
c		DO until no more noise tests are desired and process flag
c		        is set to OK
c		   DO for all desired channels
c		      CALL fnt_read_noise to collect science data
c		      IF no error occurs THEN
c		         CALL fnt_power_density to dedither series and 
c		                   calculate power density spectrum
c		         IF no error occurs THEN
c		            CALL fnt_convert_to_volts to normalize spectrum
c		            IF no error occurs THEN
c		               CALL fnt_display_noise to make pretty plots
c		*           ELSE
c		*              SET processing flag to FALSE
c		            END IF
c		*         ELSE
c		*           SET processing flag to FALSE
c		          END IF
c		*     ELSE
c		*         SET processing flag to FALSE
c		      END IF
c		   END DO
c		END DO
c		END routine
c
c    (* These error-checking sections not yet in place.  F. Shuman 88Mar29)
c--------------------------------------------------------------------------
c
c Changes:
c
c	Add FUT_SETDEV, FUT_DEADEV to allow multi-device plots. R. Kummerer,
c	May 5, 1987.
c
c	Coadd split.   F. Shuman, 1987 Dec 8.
c
c	Remove RMS write feature.   F. Shuman, 1988 Mar 29.
c
c       SPR 5129,6538, Correction of INPUT default of FNT command line
C           H. WANG, STX, 3/29/90
c
c       SPR 6625,6626, FNT Difference in spectra needed for noise monitoring
c        , FNT should call the firas lib. routine for nyquist frequency
C           H. WANG, STX, 4/12/90
c
c       SPR 9583,9790, Update FNT to use new FEX_Nyquist and FUT_Setxax.
c           Nilo G. Gonzales, Hughes STX, 1992 August 5.
c--------------------------------------------------------------------------
c
c ******************************************************
c * Modification: 1988 08 29 by Dominick P. Iascone    *
c * Purpose: To modify per SPR 1729 for Release 4.2    *
c * SPR Comments: "Lacks error trapping"               *
c ******************************************************
c ******************************************************
c * Modification: 1988 08 30 by Dominick P. Iascone    *
c * Purpose: To modify per SPR 1573 for Release 4.2    *
c * SPR Comments: All FIRAS facilities should tell of  *
c *               successful completion.               *
c ******************************************************
c ******************************************************
c * Modification: 1988 08 30 by Dominick P. Iascone    *
c * Purpose: To modify per SPR 1760 for Release 4.2    *
c * SPR Comments: FIRAS facilities signal errors via   *
c *               $Status for batch processing.        *
c ******************************************************
c ******************************************************
c * Modification: 1988 10 08 by Dominick P. Iascone    *
c * Purpose: To modify per SPR 2569                    *
c * SPR Comments: While looking at Noisetest data for  *
c *               very narrow time range, an adjustable*
c *               array error may occur.               *
c *************************************************************************
c June 23, 1988   STX   D. Bouler  SPR 3992 : Access violation
c      FNT_READ_NOISE was returning 10 for condition 
c      FNT_RNOISE_NODATA and 15 for condition FNT_RNOISE_INCONSIST.
c      Lib$signal was being called with the status value
c      rather than the address of the message;
c      this caused the access violation. Since the 
c      message was already being signaled from within FNT_READ_NOISE,
c      the call to lib$signal was deleted.
c*************************************************************************
c
c	Version 4.4.1 SER 3306 Q. CHUNG STX 08/10/89
c                    PROVIDE VERSION NUMBER TO TRACK SOFTWARE UPDATE.
c	Version 4.4.1 SER 3493, R. Kummerer, August 18, 1989.  Use new PLT
c		graphics to display IFGs and spectra.
c	Version 4.4.2 SPR 5041, 5057, R. Kummerer, Nov 15, 1989.  Inappropriate
c		FIRAS logical for fetching raw data; change CSDR$FIRAS_ARCHIVE 
c		to CSDR$FIRAS_RAW.  Correct interpretation of /INPUT.
c	Version 4.4.02, SPR 5063, Nov 16, 1989, R. Kummerer, STX
c		Direct output with CSDR$FIRAS_OUT.
c
c =========================================================================

	IMPLICIT NONE

	include '(fut_params)'
	include '(fnt_invoc)'
	include '(cct_get_config)'
	include 'ct$library:ctuser.inc'
	Include '($ssdef)'
c
c Routines called
c
	integer * 4 fnt_parse_noise
	integer * 4 fnt_display_noise
	integer * 4 fnt_convert_to_volts
	integer * 4 fnt_power_density
	integer * 4 fnt_read_noise
	integer * 4 fnt_write_noise
c
c Functions
c
	Integer*4  Cut_Register_Version
	Integer*4  Cut_Display_Banner
	Integer*4  CCT_Open_Config
c
c Records
c
	dictionary 'nfs_sdf'
	record /nfs_sdf/ sci_data(fac_max_num)
	dictionary 'fnt_noise'
	record /fnt_noise/ spec_rec
c
c GET_CONFIG variables.
c
	Dictionary 'fex_nyquist'

	Structure /config_data/
	   Record /fex_nyquist/ fex_nyquist
	Endstructure

	Record /config_data/ config
	Record /config_status/ stat(1)

	Character *  1	access_mode/' '/! data set access mode
	Integer   *  4  number/1/	! number of data sets
	Integer   *  4	size(1)/128/    ! size of data sets in bytes
	                                !  (record length).
	Character * 32  name(1)		! names of data sets
	Character * 14  ref_jstart /'86001000000000'/
	Character * 14  ref_jstop /'99365235958990'/
	Integer   *  4  jstart(2), jstop(2)

	Data name(1)/'csdr$firas_ref:fex_nyquist'/

	Integer * 4 lun(1)		! logical unit numbers
	Integer * 4 cindex(1)           ! initial cache pointers
	Logical * 1 new_segment(1)	! flag for new segments
	Logical * 1 first_time /.true./
	Integer * 4 ncache/1/
	Integer * 4 ref_count
c
c Local variables
c
	integer * 2 ct_stat(20)         !cobetrieve status return array
	integer * 4 chan_id		!channel counter
	integer * 4 start_chan		!start channel number
	integer * 4 stop_chan		!stop channel number
	integer * 4 chan_skip		!channel increment
	integer * 4 num			!number of ifgs in ensemble
	integer * 4 fake_it		!fake it bit
	integer * 4 mtm_speed		!mtm_speed
	integer * 4 u_mode		!micro processor mode
	integer * 4 r_status		!return status from function calls
        INTEGER * 4 READNOISE_status         ! SPR 2569
	integer * 4 cpy_lun		!SNAPSHOT mode lun
	logical * 1 process		!processing flag
        Logical * 1 dif_flag
	character * 2 archive_type
	integer * 4 first_IFG(2)
	integer * 4 last_IFG(2)
	integer * 4 Status              !Status of FNT processing
      	integer * 4 Success, Error      !Local status values
      	Data Success /1/, Error /2/
	Integer*4  num_vol/80/
	Integer*4  rstatus
	Integer*4  lun_out/6/
	Character*6 Version
	Parameter   (version='9.8')
c
c ******************************************
c * Added per mod 1988 08 15: SPR 2158     *
c *                                        *
      Logical * 4 transfcn_first_time
      Common /Logical_first/ transfcn_first_time
c *                                        *
c ******************************************
c 
c ******************************************
c * Added per mod 1988 08 30:              *
c * SPR 1573, 1729, 1760                   *
c *                                        *
c * FNT messages                           *
c *                                        *
	EXTERNAL FNT_normal
        EXTERNAL FNT_aberr
        EXTERNAL FUT_error
        EXTERNAL FNT_CON_VOLTS
        EXTERNAL FNT_POW_DENSITY
        EXTERNAL FNT_R_NOISE
        EXTERNAL FNT_P_NOISE
        EXTERNAL FNT_RNOISE_NODATA
        EXTERNAL FNT_RNOISE_INCONSIST !SPR 2649
	EXTERNAL FNT_CTInit
	EXTERNAL FNT_OPNCONFIGERR
		
c **********************************************
c * Items related to mod 1988 10 08/SPR 2569   *
c *                                            *
      INTEGER * 4 NO_read_chan_total, NO_chans
      INTEGER * 4 No_read_chans(4)
c *                                            *
c **********************************************
c
c-------------------------------------------------

c
c Initialize.
c
        Dif_flag = .false.
	Rstatus = Cut_Register_Version(version)
	Rstatus = Cut_Display_Banner(lun_out,Num_vol,
	1			'FIRAS Facility FNT_NoiseTest')
	write(lun_out,61)
  61	format(//)
c
      READNOISE_status = 0      !Mod 1988 10 08/SPR 2569
      NO_read_chan_total = 0 !Mod 1988 10 08/SPR 2569
c
c --- Added per mod 1988 08 15: SPR 2158 ---
      transfcn_first_time = .TRUE.
c ---                                    ---
c
c --- Added per mod 1988 08 30: SPR 1573, 1729, 1760  ---
      Status = Success
      CALL Lib$establish( FUT_error )  !Establish condition handler
c ---                                    ---
c
c	Parse command line to get channel stuff
c

	ftc_display = .true.
	ftc_jump = 0
	r_status = fnt_parse_noise(start_chan,stop_chan,chan_skip)
	If (status .eq. success) Then
	   Call CT_Init(CT_Stat)
	   If (CT_Stat(1) .Ne. CTP_Normal) Then
	      Call LIB$Signal(FNT_CTInit, %val(1), %val(ct_stat(1)))
	   End If
	End If

	call ct_gmt_to_binary (ref_jstart,jstart)
	call ct_gmt_to_binary (ref_jstop,jstop)

	If (First_time) Then
	    r_status = cct_open_config (jstart, jstop, number, name,
	1   		                size, access_mode, ncache, lun,
	2			        cindex, stat, ref_count)
	    first_time = .false.
	    If (.not. r_status) Then
	       call lib$signal(fnt_opnconfigerr,%val(1),%val(r_status))
	    End if
	End If

	do chan_id=start_chan,stop_chan,chan_skip
	   r_status = fnt_read_noise(chan_id,archive_type,first_IFG,
	1				last_IFG,num,sci_data,fake_it)
c vvv--- SPR 2569 items ---vvv
            IF (r_status .EQ. 10 .OR. r_status .EQ. 15) THEN
                READNOISE_status = r_status
                NO_read_chan_total = NO_read_chan_total + 1
                NO_read_chans(NO_read_chan_total) = chan_id
            END IF
c ^^^--- SPR 2569 items ---^^^
	    if(r_status .eq. fac_normal)then
	       r_status = fnt_power_density(chan_id,num,sci_data,spec_rec,
	1				fake_it,dif_flag,number,size,lun,
	2                               cindex,config,new_segment,stat)
	       if(r_status .eq. fac_normal)then
	          r_status = fnt_convert_to_volts(chan_id,spec_rec)
	          if(r_status .eq. fac_normal)then
	             ftc_jump = ftc_jump - 1
		     if(ftc_plots .ne. fac_not_present .and.
	1		ftc_jump .lt. 0 .and.  ftc_display)then
	                r_status = fnt_display_noise(chan_id,
	1                                       spec_rec,fake_it,dif_flag,
	2	                                number,size,lun,cindex,config,
	3	                                new_segment,stat)
		     end if
	             if(ftc_write .eq. fac_present)then
	                  r_status = fnt_write_noise(chan_id,archive_type,
	1					     first_IFG,last_IFG,
	2					     spec_rec)
	             end if
                  ELSE
                     CALL LIB$signal(%val(r_status)) !SPR 1729, 1760
                     CALL Lib$signal(FNT_con_volts)  !SPR 1729, 1760
	          end if
                 ELSE
                   CALL LIB$signal(%val(r_status))     !SPR 1729, 1760
                   CALL Lib$signal(FNT_pow_density)    !SPR 1729, 1760
	         end if
            ELSE
               CALL Lib$signal(FNT_R_noise)           !SPR 1729, 1760
	    end if
	   end do
c ***********************************************
c * Items related to mod 1988 08 30: SPR 1573   *
c *                                             *
      IF (Status .EQ. Success) THEN
        IF (READNOISE_status .EQ. 10) THEN        !SPR 2569
           CALL LIB$signal(FNT_RNOISE_NODATA)     !SPR 2569
           DO NO_chans = 1, NO_read_chan_total    !SPR 2569
              TYPE *,'No READ for channel ',NO_read_chans(NO_chans) !SPR 2569
           END DO                                 !SPR 2569
        END IF                                    !SPR 2569
c *                                             *
        IF (READNOISE_status .EQ. 15) THEN        !SPR 2649
           CALL LIB$signal(FNT_RNOISE_INCONSIST)  !SPR 2649
           DO NO_chans = 1, NO_read_chan_total    !SPR 2649
              TYPE *,'Inconsistent for channel ', fac_channel_ids(NO_chans) !SPR 2649
           END DO                                 !SPR 2649
        END IF                                    !SPR 2649
c *                                             *
        CALL Lib$signal(FNT_normal)
        CALL Exit(SS$_normal)
c *                                             *
      ELSE
        CALL Lib$signal(FNT_aberr)
        CALL Exit(SS$_abort)
      ENDIF
c *                                             *
c ***********************************************
c
        END
