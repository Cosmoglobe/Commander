  	PROGRAM FEP_DWELLPLOT

C/----------------------------------------------------------------------------
C/	PROGRAM DESCRIPTION:
C/	  This program reads DWELL-mode housekeeping records in a specified
C/	  time range and converts temperatures and bias voltages into eng'g
C/	  values for purposes of making quick plots. This obviates the
C/	  need to run Data Qualify whenever a temperature or bias plot
C/	  is desired. The program was adapted from GETARCV and ENGPLOTS.
C/
C/	AUTHOR:
C/	  R. Isaacman
C/	  Applied Research Corp.
C/	  24 March 1986
C/
C/	INPUT FILES:
C/	  FIRAS HOUSEKEEPING ARCHIVE
C/
C/	OUTPUT FILES:
C/	  None
C/
C/	INCLUDE FILES USED:
C/	  $SSDEF
C/
C/	SUBROUTINES CALLED:
C/	  FUT_TIMERANGE
C/	  FUT_OPEN_ARCHIVES
C/	  FEP_Read_FIRAS_HKP
C/	  FEP_Dwell_Address
C/	  FEP_Convert_to_Ohms-<FEP_GRTCoeffs-<FEP_Get_CalRes-<FEP_Decode_DBName
C/	  FEP_DMakePlot- - - _/FEP_DTemperatures- _/FEP_GRT_Lookup- -<FEP_Invert
C/	                      \FEP_DPlotter        \GRT_Correction
C/	  LIB$ESTABLISH			            \FEP_GET_GRT_CONV		
C/	  CT_INIT
C/	  LIB$ERASE_PAGE
C/	  CT_CLOSE_ARCV
C/        STR$UPCASE
C/	  LIB$SIGNAL
C/	  EXIT
C/
C/	ERROR HANDLING:
C/	  NONE
C/
C/	METHOD USED:
C/	  Trivial
C/
c----------------------------------------------------------------------
c	Modified by W. Young
c		    STI corp formerly SASC TECH
c		    May 14, 1986
c	Reason:
c		Updated CT return status and changed calling sequence for
c		READ_FIRAS_HSK to conform with the new version
c
C	November 4, 1987, R. Kummerer.
C	   Treat FEP_POLYSLNERR as an informational error instead of fatal.
C
C	December 11, 1987, R. Kummerer.
C	   SPR 1507, Signal TBLEDGE message once per plot.
C
c	L. Olson  6/7/88:  Changed reference to N_OPEN_ARCHIVES (library FUTI)
c	   to FUT_OPEN_ARCHIVES in FUT
c
c	F. Shuman  1988 08 01 : Convert module DPLOTTER to PLT (SER 2065)
c
c	F. Shuman  1988 08 10 : Replace non-functioning LIB$MovC3 with
c	   UNION / MAP...ENDMAP / to equivalence the BUFF array with the
c	   housekeeping record structure.  Also, loop on user's choice of
c	   a new timerange, as is done in ENGPLOTS.
c
c	Shirley M. Read, STX, 09/09/88: Changed Oper_Entry call to
c          FUT_Timerange. The FUTILIB can now be removed from the link.
c          SPRs 1566 and 1763 requested that all FIRAS facilities
c          should tell of successful completion and signal errors
c	   via $status for batch processing.
C
C	R. Kummerer, STX, January 12, 1989.
C	   SPR 3132: Dwellplot fails to plot XCAL S5 and S6.
C
C	Fred Shuman, STX, 1989 Feb 10.
C	   SPR 2922: Spurious and misplaced values.
C
C	Fred Shuman, STX,  1989 Feb 17.
C	   SPR 2910, Plot all fields dwelled in the timerange rather than just
C	   the one with the plurality of hits on each side (A/B).
C
C	Fred Shuman, STX,  1989 Feb 20.
C	   SPR 2738, Add the ability to plot cal resistors (in counts).
C
C	Fred Shuman, STX,  1989 Aug 30.
C	   SER 3306, Call CUT_Register_Version and CUT_Display_Banner to
C	   generate a 'banner'.
C
C       H. Wang, STX, 1990 Feb. 22
C          Change version to 5.7
C
C       H. Wang, STX, 1990 Mar. 2
C          Change version to 5.8
C          Spr 3253, Dwellplot has to fetch calres params from reference
C          data archive. 
C
C       H. Wang, STX, 1990 Mar. 5
C          Change version to 5.8
C          Spr 6408, Dwellplot Currently plots fake dwell data
C
C       H. Wang, STX, 1991 Jul. 10
C          Change version to 9.0
C          Spr 8731, Dwellplot users obsolete routine fut_open_archives
C           Make the small changes to FEP_DWELLPLOT to remove the call to
C           FUT_OPEN_ARCHIVES   
C
C------------------------------------------------------------------------

	IMPLICIT	NONE

	DICTIONARY	'NFS_HKP'

	STRUCTURE /HOUSEKP/
	   UNION
	      MAP
	         RECORD	/NFS_HKP/ HKP_REC
	      END MAP
	      MAP
	         BYTE BUFF(576)
	      END MAP
	   END UNION
	END STRUCTURE

	RECORD /HOUSEKP/ HK

	Integer		*4	CUT_Register_Version
	Integer		*4	CUT_Display_Banner
	Integer		*4	num_vol/80/
	Integer		*4	lun_out/6/
	Character	*6	version
	Parameter		(version='9.0')

	CHARACTER	COMMENT*80, END_TIME*14,
	1		FNAME*7, LINE*80, START_TIME*14,
	1		another_time*1

	integer*4       Fut_Timerange		!FUT library request time
	integer*4       STR$UPCASE
	integer*4       binstart(2), binstop(2)
        Character*32    diset(2), dset
        Data Diset(1)   /'Csdr$Firas_Ref:FDB_GRTCAL'/
        Data Diset(2)   /'Csdr$Firas_Ref:FDB_GRTRICH'/
        Data Dset       /'Csdr$Firas_Ref:FEX_CALRES'/
        Integer*4       Nidset/2/
        Integer*4       Ndset/1/
        Integer*4       Config_lun(2)
        Integer*4       Con_Lun
        Integer*4       Index(2)
        Integer*4       Index1
        Integer*4       Ncache/1/
        Integer*4       Ref_Count
        Integer*4       Size(2)/20,7120/
        Integer*4       Size1/256/
        Character*5     Access_mode/'KEYED'/
        Character*5     Access_mode1/' '/
        Integer*4       CCT_OPEN_CONFIG
        Integer*4       CCT_CLOSE_CONFIG   

	LOGICAL *1	MORE
	LOGICAL *1	MOREA
	LOGICAL *1	MOREB

	INTEGER*2	COA_LUN_EQV, CT_LUN(40),
     .			CT_STAT(20), ENG_LUN_EQV,
     .			HKP_LUN_EQV, IDX_LUN_EQV, ISTAT,
     .			adwellst_mf1, bdwellst_mf1,
     .			adwellst_mf2, bdwellst_mf2,
     .			dwella(1000), dwellb(1000),
     .			dwadressa1, dwadressb1, dwadressa2, dwadressb2
	integer*2	dwellflaga, dwellflagb, nrec, Telm_mode,
     .			nreca, nrecb, SCI_LUN_EQV, Switch, WRITE_LUN
        integer * 4     cat_entry, status, istatus, TR$UpCase, rcode
	integer*4	success/1/, error/2/  ! Local status values

	real*4		ohmsa(64000), ohmsb(64000)

	character*14	timetaga(1000), timetagb(1000)
	PARAMETER	ASK_OPER = 0
	PARAMETER	GET_HSK = 1
	PARAMETER	GET_IDX = 2
	PARAMETER	GET_SCI = 3
	PARAMETER	GET_ENG = 4
	PARAMETER	GET_COA = 5

	EQUIVALENCE	(CT_LUN(1),	HKP_LUN_EQV),
     .			(CT_LUN(2),	IDX_LUN_EQV),
     .			(CT_LUN(3),	SCI_LUN_EQV),
     .			(CT_LUN(4),	ENG_LUN_EQV),
     .			(CT_LUN(5),	COA_LUN_EQV)
	parameter (adwellst_mf1 = 89)
	parameter (bdwellst_mf1 = 90)
	parameter (adwellst_mf2 = 305)
	parameter (bdwellst_mf2 = 306)

c	FEP messages

	external	fep_normal
	external	fep_aberr
	external	fep_opnconfigerr
	external	fep_clsconfigerr

c	Condition Handler

	external        fut_error
	external        fut_normal

	include         '($ssdef)'
        include         '(CCT_GET_CONFIG)'
C
        Record/Config_Status/ Stat(2)
        Record/Config_Status/ Stat1

        Character*64 Infile
        Character*30 time_range
        integer*4   ct_connect_read,ios
        external    Ct_connect_read            
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	status = CUT_Register_Version(version)
	status = CUT_Display_Banner(lun_out, num_vol,
	1                           'FIRAS Facility  FEP_Dwellplot ')
	Write(lun_out,61)
61	Format(//)

c     Establish condition handler.

	call lib$establish ( fut_error )

c     Set status to success.

	istatus = success

	call ct_init (ct_stat)

	another_time = 'Y'
	Do While (another_time(1:1) .Eq. 'Y')


	   nrec  = 0
	   nreca = 0
	   nrecb = 0
	   switch = get_hsk
c
c  Get time range for extracting housekeeping data
c
10	   STATUS= FUT_TIMERANGE (START_TIME, binstart, END_TIME, binstop)
	   IF (STATUS .NE. %loc(FUT_NORMAL)) THEN
	     TYPE *,' Please try again!'
	     GO TO 10
	   ENDIF


	   time_range = start_time // ';' //
	2                             end_time // ';'
c Open the housekeeping archive.
c
	         infile = 'CSDR$FIRAS_RAW:NFS_HKP/' // time_range
	  	 ct_lun(switch) = 1

	         Open ( UNIT = ct_lun(switch),
	1               FILE = infile,
	2               STATUS = 'old', IOSTAT = ios,
	3               USEROPEN = CT_Connect_Read )

                 If (ios .ne. 0) then
   	          TYPE *,' FILE OPEN FAILED,  Please try again!'
	          GO TO 10
                 Endif       
           Status = CCT_OPEN_Config(Binstart, Binstop,Nidset,Diset, Size,
	1	    access_mode,ncache,config_lun, index, stat, ref_count)
           If (.not. status) then
             Call Lib$Signal(Fep_Opnconfigerr,%val(1),%val(status))
             Call Exit(SS$_Abort)
           Endif
           Status = CCT_OPEN_Config(Binstart, Binstop,Ndset,Dset, Size1,
	1	    access_mode1,ncache,con_lun, index1, stat1, ref_count)
           If (.not. status) then
             Call Lib$Signal(Fep_Opnconfigerr,%val(1),%val(status))
             Call Exit(SS$_Abort)
           Endif
	   MORE = .TRUE.
	   MOREA = .TRUE.
	   MOREB = .TRUE.
	   DO WHILE (MORE)

	      CALL FEP_READ_FIRAS_HKP (HKP_LUN_EQV, HK.HKP_REC, CAT_ENTRY,
	2                                                           ISTAT)
C
C  Recall that HKP_REC and HK.BUFF are co-mapped via the structure HOUSEKP
C    at the top of this program.
C
              Telm_mode = 0
	      IF (ISTAT .EQ. 1) THEN
	         nrec = nrec + 1
                 If ((HK.BUFF(523) .ne. 0) .or. (Hk.buff(555) .ne. 0))
	1	          Telm_mode = 1
                 If (Telm_mode .EQ. 0) then
	         call fep_dwell_address (HK.BUFF(adwellst_mf1), dwadressa1)
	         call fep_dwell_address (HK.BUFF(bdwellst_mf1), dwadressb1)
	         call fep_dwell_address (HK.BUFF(adwellst_mf2), dwadressa2)
	         call fep_dwell_address (HK.BUFF(bdwellst_mf2), dwadressb2)
	         dwellflaga = (dwadressa1 + dwadressa2)/2
	         dwellflagb = (dwadressb1 + dwadressb2)/2
	         If (dwellflaga .ge. 0 .and. morea) Then
	            nreca = nreca + 1
	         if (nreca.gt.1000 .and. istat.eq.1) then
	         type 20
20            format (///' WARNING: ** Truncating data set after 1000 records',
	1               'for A side')
                 Morea = .false.
                 else
	            dwella(nreca) = dwellflaga
	            Call FEP_Convert_to_Ohms (HK.buff, nreca, ohmsa,
	1	                             config_lun(1),
	2                                  con_lun,dwella(nreca), timetaga, 0,
	3                                  rcode )
                 If (rcode .ne. 0) then
                    Call Exit(SS$_Abort)
                 Endif
	         End If
                 Endif
	         If (dwellflagb .ge. 0 .and. moreb) Then
	            nrecb = nrecb + 1
	         if (nrecb.gt.1000 .and. istat.eq.1) then
	         type 40
40            format (///' WARNING: ** Truncating data set after 1000 records',
	1               'for B side')
                 Moreb = .false.
                 else
	            dwellb(nrecb) = dwellflagb
	            Call FEP_Convert_to_Ohms (HK.buff, nrecb, ohmsb,
	1  	                           config_lun(1),
	2                                  con_lun,dwellb(nrecb), timetagb, 1,
	3	                           rcode )
                 If (rcode .ne. 0) then
                    Call Exit(SS$_Abort)
                 Endif
	         End If
                endif 
                Endif
	      ELSEIF (ISTAT .EQ. -1) THEN	! End of file
	         MORE = .FALSE.
              Elseif ((.not. morea) .and. (.not. moreb)) then
                 More = .false.
  	      ELSE
	         type *, 'Archive read status =', istat,'. DWELLPLOT aborts.'
	         istatus = error
		 more = .false.
		 another_time(1:1) = 'N'
	      ENDIF

	   ENDDO

	   Call CT_Close_Arcv( , CT_LUN, CT_Stat)

	   if ( istatus .eq. success ) then
	     if (nreca .gt. 0) then
	       call fep_dmakeplot (config_lun(2),nreca, ohmsa, dwella, 
	1	       timetaga, 0)
	     else
	       type *, '**** There were no DWELL-mode records on side A. ****'
	     endif
	   endif

	   if ( istatus .eq. success ) then
	     if (nrecb .gt. 0) then
	       call fep_dmakeplot (config_lun(2),nrecb, ohmsb, dwellb, 
	1	      timetagb, 1)
	     else
	       type *, '**** There were no DWELL-mode records on side B. ****'
	     endif
	   endif
           Status = CCT_Close_Config(Nidset,Config_lun,Index)
           If (.not. Status) then
             Call Lib$Signal(Fep_Clsconfigerr,%val(1),%val(status))
             Call Exit(SS$_Abort)
           Endif   
           Status = CCT_Close_Config(Ndset,Con_lun,Index1)
           If (.not. Status) then
             Call Lib$Signal(Fep_Clsconfigerr,%val(1),%val(status))
             Call Exit(SS$_Abort)
           Endif   
C           CALL LIB$ERASE_PAGE (1,1)

	   if ( istatus .eq. success ) then
	     Type *, ' '
	     Type *, 'Would you like another timerange? (Y/[N])'
	     Accept 50, another_time
50	     format (a)
	     status = STR$UpCase(another_time, another_time)
	   endif
	End Do

c       Exit the program.

	if ( istatus .eq. success ) then
	  call lib$signal(fep_normal)
	  call exit(ss$_normal)
	else
	  call lib$signal(fep_aberr)
	  call exit(ss$_abort)
	endif
	END
