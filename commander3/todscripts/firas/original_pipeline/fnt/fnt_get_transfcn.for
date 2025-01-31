      INTEGER * 4 function fnt_get_transfcn(Bin_time_ave, ave_num,
     &                     chan_id, mtm_speed,
     &                     fake_it,u_mode,ngroup,transfcn,preamp)

c--------------------------------------------------------------------------
c
c	This function gets the correct electronix transfer function
c		from the ELEX_TRANSFCN.DAT file
c
c	Stolen by: W. K. Young
c		   STI Inc.
c		   July 1986
c	From: 	   R. Isaacman's routine get_recnum
c
c--------------------------------------------------------------------------
c
c Changes:
c
c	R. Kummerer, Apr 20, 1987 to point to FUT_GET_RECNUM.
c
c--------------------------------------------------------------------------
c ************************************************
c * Modification: 1988 08 15 for Build 4.1 .     *
c * Purpose: To implement GET_Config per SPR 2158*
c * Programmer: Nick Iascone (with lots of       *
c *             Shirley Read help).              *
c ************************************************
c ************************************************
c * Modification: 1988 09 28                     *
c * Programmer: Nick Iascone, STX, COBE, Bldg 15 *
c * Purpose: To implement change phoned in by    *
c *          Rich Isaacman at approx 1500 (local)*
c *          on 28 September 1988.               *
c * Mod: To reverse the order of file names in   *
c *      diset. This was accomplished by an ex-  *
c *      change of index values in the DATA      *
c *      STATEMENTS for diset.                   *
c *      Apparently, while making mods for       *
c *      SPR 2158 (Mod 1988 08 15), I made the   *
c *      file name assignments backwards.        *
c ************************************************
c ************************************************
c * Modification: 1989 09 13                     *
c * Programmer: Harte Wang, STX, COBE, Bldg 15   *
c * Purpose: To fix the SPR # 4519               *   
c *                                              *
c * Mod:  change the pointer reference from      *
c *       csdr$firas_archive to csdr$firas_ref,  *
c *       also add new message when read error   *
c *       occurred.                              *
c ************************************************
c   SER 8919, S. Alexander, HSTX, August 3, 1992.  Electronics transfer
c             function changed to double precision.
c

      IMPLICIT NONE

      include '(fut_params)'
      include '(fnt_invoc)'  ! To acquire FTC_start, FTC_stop for
                             ! Get_config use.

	integer * 4	chan_id			!channel id
	integer * 4	status			!return status
	integer * 4	er			!logical unit number
	integer * 4	ios			!return status from open&read
	integer * 4	mtm_speed		!mtm speed
	integer * 4	fake_it			!fake it mode flag
	integer * 4	u_mode			!micro mode
	integer * 4	ngroup			!adds per group
	integer * 4	nrec_elx		!electonics transfer function id
	complex *16	transfcn(257)		!electronics transfer function

	integer * 4	fut_get_recnum
	logical * 4	get_xfcnfile
	integer * 4	preamp

	external	fut_normal
c
c ***************************************************
c * Items related to SPR 2158 mods (use GET_Config) *
c *                                                 *
        External FNT_GETconfigerr
        External FNT_OPNconfigerr
        External FNT_READERR
        Include '(cct_get_config)'
      Include '($ssdef)'
c *                                                 *
      RECORD /Config_status/ stat2(2)
c *                                                 *
      Integer * 4 Size2(2)    ! Size of datasets in bytes
      Data Size2 /4112, 2056/
c *                                                 *
cc      Character * 6 Access_mode2(2)  ! Access mode
cc      Data Access_mode2 /2 * 'DIRECT'/
      Character * 6 Access_mode2  ! Access mode
      Data Access_mode2 /'DIRECT'/
c *                                                 *
      Character * 32 diset(2) ! Name of datasets
      Data diset(2) /'CSDR$FIRAS_REF:fex_DTF'/
      Data diset(1) /'CSDR$FIRAS_REF:fex_ETF'/
c *                                                 *
      Integer * 4 nidset /2/        ! Number of datasets
      Integer * 4 config_lun(2)     ! Logical unit numbers
      Integer * 4 Logical_unit      ! Assumes one of config_lun's vals
      Integer * 4 current_lun       ! Current logical unit number
      Integer * 4 Index2(2)         ! Initial cache pointers
      Logical * 1 new_segment2(2)   ! Flag for new segments accesses
      Integer * 4 ncache /1/        ! Number of caches - not used here
      Integer * 4 ref_count         ! Reference counter for cache
      Integer * 4 retstat           ! Status returned by function
      Integer * 4 success, error_st ! Local status values
      Data success /1/, error_st /2/
c *                                                 *
      Integer * 4 CCT_OPEN_CONFIG
      Integer * 4 CCT_GET_CONFIG_IDX_TOD
c *                                                 *
      Integer * 4 Start_time(2), End_time(2)  ! Complete data range for
                                               ! the run
      Common /Times/ Start_time, End_time
      Integer * 4 AVG_time(2)                  ! Time for a particular
                                               ! record being processed
      Integer * 4 ave_num  ! The number of average times
      Integer * 4 bin_time_ave(ave_num)  ! The average times
      Logical * 4 Transfcn_first_time
      Common /Logical_first/ Transfcn_first_time
c *                                                 *
c *                                                 *
c ***************************************************

	
	fnt_get_transfcn = fac_normal
c

c --- Set status for FNT processing to success. ---
      Status = Success

c
c --------------------------------------------------------------
c - Replace with CCT_OPEN_CONFIG call                          -
c -                                                            -
cc	er = 99			!logical unit assignment
cc	inquire(file='csdr$firas:elex_transfcn.dat',
cc     .			opened=get_xfcnfile)
cc	if(.not.get_xfcnfile)then
cc	   if(preamp .eq. fac_present)then
cc
cc   	      open(unit=er,file='csdr$firas:elex_transfcn.dat',
cc     .		status ='OLD',
cc     .		FORM='UNFORMATTED', ACCESS='DIRECT', RECL=514,
cc     .		READONLY,iostat=ios)
c -                                                            -
c --- OPEN Configuration files. Time range is passed from    ---
c --- FNT_parse_noise using Common /FTC/ .                   ---
c ---                                                        ---
      IF (status .EQ. success) THEN
        IF (Transfcn_first_time) THEN
          Transfcn_first_time = .FALSE.
          RETSTAT = CCT_OPEN_CONFIG(start_time, end_time,
     &              nidset, diset, size2, access_mode2, ncache,
     &              config_lun, index2, stat2, ref_count)
          IF (.NOT. RETSTAT) THEN
            status = Error_st
            CALL Lib$signal(FNT_OPNCONFIGERR, %val(1), %val(RETSTAT))
          END IF
        END IF
      END IF
      
c --------------------------------------------------------------
cc	      if(ios .ne. 0)then
cc	         type *,'Failed to open electronics transfer function file'
cc	         fnt_get_transfcn = 13
cc	      end if
cc	   else
c --------------------------------------------------------------
c - Replaced by CCT_OPEN_CONFIG call above, which OPENed both  -
c - datasets.                                                  -
c -                                                            -
cc           open(unit=er,file='CSDR$FIRAS:DGTL_TRANSFCN.DAT',
cc     .          status ='OLD',
cc     .          FORM='UNFORMATTED', ACCESS='DIRECT', RECL=514,
cc     .          READONLY,iostat=ios)
c -                                                            -
c --------------------------------------------------------------
cc	      if(ios .ne. 0)then
cc	         type *,'Failed to open digital filter file'
cc	         fnt_get_transfcn = 13
cc	      end if
cc	   end if
cc	end if
c
c *************************************************************
c * Added per mod 1988 08 15                                  *
c * CALL COBETRIEVE function to get all conversion datasets.  *
c *                                                           *
      IF (status .EQ. success) THEN
        RETSTAT = CCT_GET_CONFIG_IDX_TOD(Bin_time_ave, nidset,
     &            config_lun, Index2, new_segment2, stat2)
        IF (.NOT. RETSTAT) THEN
          Status = error_st
          CALL Lib$signal(FNT_getconfigerr, %val(1),
     &     %val(RETSTAT))
        END IF
      END IF
c *                                                           *
c *************************************************************
      
c  Retrieve appropriate transfer fcn based on instrument
c  configuration, scan mode, etc.
c
	IF (ngroup .EQ. 0) ngroup = 1
	status = FUT_GET_RECNUM (FAKE_IT, MTM_SPeeD, CHAN_id, u_mode,
     .					ngroup, NREC_ELX)
	if(status .eq. %loc(fut_normal))then
c ----------------------------------------------------------
c - This replaced by decision on which dataset to be READ. -
c - This caused because new GET_config approach causes both-
c - possible datasets to be OPENed at once. So, instead of -
c - a decision being made at OPEN, it is made at READ.     -
c -                                                        -
          Logical_unit = config_lun(2)
          IF (preamp .eq. fac_present) Logical_unit = config_lun(1)
c
          READ (Logical_unit,rec=NREC_ELX,iostat=ios) TRANSFCN
          if (ios .ne. 0) then
           CALL Lib$signal(FNT_readerr)        
cc	type *, 'Failed to read electronics transfer function'
           fnt_get_transfcn = 13
          end if
      else
	  fnt_get_transfcn = 13
      end if
	return
	end
