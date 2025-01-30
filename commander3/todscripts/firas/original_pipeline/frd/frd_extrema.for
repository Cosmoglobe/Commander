	Program FRD_EXTREMA

C-----------------------------------------------------------------------------
C
C 	Program Name:
C 	  FRD_EXTREMA
C 
C 	Program Description:
C         This program is used for updating extrema values as recorded
C         in the extrema file. The current date is used when the record
C         is updated, and a new FXT file is written with a time tag 
C         from the update time to Cobetrieve infinity (day 365 of year 1999).
C         Only the most recent extrema file will be read in and updated.
C
C	Programmers:
C	  David B. Bouler, STX, Sept, 1990.
C         (Screen management adapted from FEP by Bob Kummerer)
C 
C 	Calling Sequence:
C 	  This is a main program.
C 
C 	Input/Output Parameters:
C 	  None
C 
C 	Input Files:
C         FIRAS log extrema file, FXT_ENG_XTRM
C       
C 	Output Files:
C         FIRAS log extrema file, FXT_ENG_XTRM
C 	  Report file for run :   FRD_EXTREMA_REPORT.LIS
C 
C 	Include Files Used:
C 	  CT$Library:CTUser.Inc (COBETRIEVE return status definitions)
C         $SMGDEF               (SCREEN MANAGER DEFINITIONS)
C         $SSDEF                
C         cct_filespec_fields_record (get info on old FXT file)
C 
C 	Subroutines and Functions Called:
C         CCT_Get_Filespec_fields
C	  CT_Read_Arcv
C	  CT_Write_Arcv
C	  Lib$Establish
C	  Lib$Signal
C         OTS$Cvt_TI_L
C         OTS$Cvt_T_F
C	  SMG$Create_Virtual_Display
C	  SMG$Create_Virtual_Keyboard
C	  SMG$Delete_Virtual_Keyboard
C	  SMG$Paste_Virtual_Display
C	  SMG$Draw_Rectangle
C	  SMG$Draw_Line
C	  SMG$Put_Chars
C	  SMG$Ring_Bell
C	  SMG$Read_String
C	  SMG$Set_Cursor_Abs
C	  SMG$Create_PasteBoard
C	  SMG$Erase_PasteBoard
C	  SMG$erase_chars
C	  SMG$Delete_PasteBoard
C	  SMG$Delete_Virtual_Display
C
C 	Error Handling:
C	  Establish the condition handler. Signal all errors via Lib$Signal.
C 	  Error messages to SYS$OUTPUT (terminal if run interactively, log
C 	  file if run batch) 
C 
C       Change Log:
C
C
C
C----------------------------------------------------------------------------- 

       implicit none

       include  '($smgdef)'
       include  '($ssdef)'
       include  'ct$library:ctuser.inc'
       include  '(cct_filespec_fields_record)'

       record    /filespec_fields/  specs
c
c Variables for menus
c
       integer*4    minmax                       !Whether min or max is being
                                                 !updated (Min = 1, Max = 2)

       integer*4    max_lines                    !Max number of lines on a menu
       parameter    (max_lines = 20)

       integer*4    num_menus                    !Number of menus 
       parameter    (num_menus = 10)

       integer*4    num_fields(num_menus)        !Number of fields in menu

       real*4       update_values(2,num_menus,max_lines)  !Updated values
       real*4       old_values(2,num_menus,max_lines)     !Values before update

       logical*1    update_flags(2,num_menus,max_lines)   !(value updated?)

       character*14 old_gmts (2,num_menus,max_lines)      !Times before update

       character*65 field_names(2,num_menus,max_lines)    !Names of menu fields
       character*65 menu_titles(2,num_menus)              !Titles of menus
       character*65 main_menu(num_menus+3)                !Main menu text
       character*65 ret_line                              !Text for main return

       data ret_line  / '99   Return to main menu' /  

       data main_menu / '1    Update minima (default)         ',
     .                  '2    Update maxima                   ',
     .                  '3    GRT temps side A low  current   ' ,
     .                  '4    GRT temps side A high current   ' ,
     .                  '5    GRT temps side B low  current   ' ,
     .                  '6    GRT temps side B high current   ' ,
     .                  '7    Misc temps and currents         ' ,  
     .                  '8    Bolometer voltages              ' , 
     .                  '9    IPDU digital and analog voltages' , 
     .                  '10   IPDU other voltages             ',
     .                  '11   IPDU currents                   ' , 
     .                  '12   LMAC temps                      ' ,  
     .                  '99   Quit update of extrema          '   /

       data menu_titles(1, 1) / 'GRT temps side A low  current' /
       data menu_titles(1, 2) / 'GRT temps side A high current' /
       data menu_titles(1, 3) / 'GRT temps side B low  current' /
       data menu_titles(1, 4) / 'GRT temps side B high current' /
       data menu_titles(1, 5) / 'Misc temps and currents      ' / 
       data menu_titles(1, 6) / 'Bolometer voltages           ' /
       data menu_titles(1, 7) / 'IPDU digital and analog voltages' /
       data menu_titles(1, 8) / 'IPDU other voltages          ' /
       data menu_titles(1, 9) / 'IPDU currents                ' /
       data menu_titles(1,10) / 'LMAC temps                   ' / 

       data num_fields( 1) /16/  ! Side A low current GRTs 
       data num_fields( 2) /16/  ! Side A hi  current GRTs 
       data num_fields( 3) /16/  ! Side B low current GRTs 
       data num_fields( 4) /16/  ! Side B hi  current GRTs 
       data num_fields( 5) /16/  ! Misc temps and currents
       data num_fields( 6)  /4/  ! Bolomoter voltages
       data num_fields( 7) /10/  ! IPDU digital and analog voltages
       data num_fields( 8) /10/  ! IPDU other voltages
       data num_fields( 9) /12/  ! IPDU currents
       data num_fields(10)  /2/  ! LMAC temps 

!      GRTs Side A, Low Current

       data field_names( 1, 1, 1) / '1    XCAL_TEMP_AL' /
       data field_names( 1, 1, 2) / '2    SKY_HORN_TEMP_AL' /
       data field_names( 1, 1, 3) / '3    REFERENCE_HORN_TEMP_AL' /
       data field_names( 1, 1, 4) / '4    INTERNAL_CALIBRATOR_TEMP_AL' /
       data field_names( 1, 1, 5) / '5    DIHEDRAL_TEMP_AL' /
       data field_names( 1, 1, 6) / '6    BOLOMETER_ASSEMBLY_RH_TEMP_AL' /
       data field_names( 1, 1, 7) / '7    BOLOMETER_ASSEMBLY_RL_TEMP_AL' /
       data field_names( 1, 1, 8) / '8    BOLOMETER_ASSEMBLY_LH_TEMP_AL' /
       data field_names( 1, 1, 9) / '9    BOLOMETER_ASSEMBLY_LL_TEMP_AL' /
       data field_names( 1, 1,10) / '10   MIRROR_TEMP_AL' /
       data field_names( 1, 1,11) / '11   CAL_RESISTORS_RH_TEMP_AL' /
       data field_names( 1, 1,12) / '12   CAL_RESISTORS_RL_TEMP_AL' /
       data field_names( 1, 1,13) / '13   CAL_RESISTORS_LH_TEMP_AL' /
       data field_names( 1, 1,14) / '14   CAL_RESISTORS_LL_TEMP_A' /
       data field_names( 1, 1,15) / '15   XCAL_TEMP_S5_AL' /
       data field_names( 1, 1,16) / '16   COLLIMATOR_TEMP_AL' /  

!      GRTs Side A, High current.
           
       data field_names( 1, 2, 1) / '1    XCAL_TEMP_AH' /
       data field_names( 1, 2, 2) / '2    SKY_HORN_TEMP_AH' /
       data field_names( 1, 2, 3) / '3    REFERENCE_HORN_TEMP_AH' /
       data field_names( 1, 2, 4) / '4    INTERNAL_CALIBRATOR_TEMP_AH' /
       data field_names( 1, 2, 5) / '5    DIHEDRAL_TEMP_AH' /
       data field_names( 1, 2, 6) / '6    BOLOMETER_ASSEMBLY_RH_TEMP_AH' /
       data field_names( 1, 2, 7) / '7    BOLOMETER_ASSEMBLY_RL_TEMP_AH' /
       data field_names( 1, 2, 8) / '8    BOLOMETER_ASSEMBLY_LH_TEMP_AH' /
       data field_names( 1, 2, 9) / '9    BOLOMETER_ASSEMBLY_LL_TEMP_AH' /
       data field_names( 1, 2,10) / '10   MIRROR_TEMP_AH' /
       data field_names( 1, 2,11) / '11   CAL_RESISTORS_RH_TEMP_AH' /
       data field_names( 1, 2,12) / '12   CAL_RESISTORS_RL_TEMP_AH' /
       data field_names( 1, 2,13) / '13   CAL_RESISTORS_LH_TEMP_AH' /
       data field_names( 1, 2,14) / '14   CAL_RESISTORS_LL_TEMP_AH' /
       data field_names( 1, 2,15) / '15   XCAL_TEMP_S5_AH' /
       data field_names( 1, 2,16) / '16   COLLIMATOR_TEMP_AH' /  

!      GRTs Side B, Low current.

       data field_names( 1, 3, 1) / '1    XCAL_TEMP_BL' /
       data field_names( 1, 3, 2) / '2    SKY_HORN_TEMP_BL' /
       data field_names( 1, 3, 3) / '3    REFERENCE_HORN_TEMP_BL' /
       data field_names( 1, 3, 4) / '4    INTERNAL_CALIBRATOR_TEMP_BL' /
       data field_names( 1, 3, 5) / '5    DIHEDRAL_TEMP_BL' /
       data field_names( 1, 3, 6) / '6    BOLOMETER_ASSEMBLY_RH_TEMP_BL' /
       data field_names( 1, 3, 7) / '7    BOLOMETER_ASSEMBLY_RL_TEMP_BL' /
       data field_names( 1, 3, 8) / '8    BOLOMETER_ASSEMBLY_LH_TEMP_BL' /
       data field_names( 1, 3, 9) / '9    BOLOMETER_ASSEMBLY_LL_TEMP_BL' /
       data field_names( 1, 3,10) / '10   MIRROR_TEMP_BL' /
       data field_names( 1, 3,11) / '11   CAL_RESISTORS_RH_TEMP_BL' /
       data field_names( 1, 3,12) / '12   CAL_RESISTORS_RL_TEMP_BL' /
       data field_names( 1, 3,13) / '13   CAL_RESISTORS_LH_TEMP_BL' /
       data field_names( 1, 3,14) / '14   CAL_RESISTORS_LL_TEMP_BL' /
       data field_names( 1, 3,15) / '15   XCAL_TEMP_S6_BL' /
       data field_names( 1, 3,16) / '16   COLLIMATOR_TEMP_BL' /  

!      Side B, High current.

       data field_names( 1, 4, 1) / '1    XCAL_TEMP_BH' /
       data field_names( 1, 4, 2) / '2    SKY_HORN_TEMP_BH' /
       data field_names( 1, 4, 3) / '3    REFERENCE_HORN_TEMP_BH' /
       data field_names( 1, 4, 4) / '4    INTERNAL_CALIBRATOR_TEMP_BH' /
       data field_names( 1, 4, 5) / '5    DIHEDRAL_TEMP_BH' /
       data field_names( 1, 4, 6) / '6    BOLOMETER_ASSEMBLY_RH_TEMP_BH' /
       data field_names( 1, 4, 7) / '7    BOLOMETER_ASSEMBLY_RL_TEMP_BH' /
       data field_names( 1, 4, 8) / '8    BOLOMETER_ASSEMBLY_LH_TEMP_BH' /
       data field_names( 1, 4, 9) / '9    BOLOMETER_ASSEMBLY_LL_TEMP_BH' /
       data field_names( 1, 4,10) / '10   MIRROR_TEMP_BH' /
       data field_names( 1, 4,11) / '11   CAL_RESISTORS_RH_TEMP_BH' /
       data field_names( 1, 4,12) / '12   CAL_RESISTORS_RL_TEMP_BH' /
       data field_names( 1, 4,13) / '13   CAL_RESISTORS_LH_TEMP_BH' /
       data field_names( 1, 4,14) / '14   CAL_RESISTORS_LL_TEMP_BH' /
       data field_names( 1, 4,15) / '15   XCAL_TEMP_S6_BH' /
       data field_names( 1, 4,16) / '16   COLLIMATOR_TEMP_BH' /
 
!      Temperatures and Currents

       data field_names( 1, 5, 1) / '1    IPDU_TEMP_A' /
       data field_names( 1, 5, 2) / '2    IPDU_TEMP_B' /
       data field_names( 1, 5, 3) / '3    CHANNEL_TEMP_RH' /
       data field_names( 1, 5, 4) / '4    CHANNEL_TEMP_RL' /
       data field_names( 1, 5, 5) / '5    CHANNEL_TEMP_LH' /
       data field_names( 1, 5, 6) / '6    CHANNEL_TEMP_LL' /
       data field_names( 1, 5, 7) / '7    DRIVE_BOX_TEMP_A' /
       data field_names( 1, 5, 8) / '8    DRIVE_BOX_TEMP_B' /
       data field_names( 1, 5, 9) / '9    STATUS_MONITOR_TEMP_A' /
       data field_names( 1, 5,10) / '10   STATUS_MONITOR_TEMP_B' /
       data field_names( 1, 5,11) / '11   CHAN_PREAMP_TEMP' /
       data field_names( 1, 5,12) / '12   OPTICAL_PREAMP' /
       data field_names( 1, 5,13) / '13   HOT_SPOT_CURRENT_A' /
       data field_names( 1, 5,14) / '14   HOT_SPOT_CURRENT_B' /
       data field_names( 1, 5,15) / '15   MTM_CAL_MOTOR_A' /
       data field_names( 1, 5,16) / '16   MTM_CAL_MOTOR_B' /

!      Bolometer voltages

       data field_names( 1, 6, 1) / '1    BOLOMETER_BIAS_VOLTAGE_RH' /
       data field_names( 1, 6, 2) / '2    BOLOMETER_BIAS_VOLTAGE_RL' /
       data field_names( 1, 6, 3) / '3    BOLOMETER_BIAS_VOLTAGE_LH' /
       data field_names( 1, 6, 4) / '4    BOLOMETER_BIAS_VOLTAGE_LL' /
 
!      IPDU digital and analog voltages.

       data field_names( 1, 7, 1) / '1    DIGITAL_CONV_N15V_A' /
       data field_names( 1, 7, 2) / '2    DIGITAL_CONV_N15V_B' /
       data field_names( 1, 7, 3) / '3    DIGITAL_CONV_P15V_A' /
       data field_names( 1, 7, 4) / '4    DIGITAL_CONV_P15V_B' /
       data field_names( 1, 7, 5) / '5    DIGITAL_CONV_P5V_A' /
       data field_names( 1, 7, 6) / '6    DIGITAL_CONV_P5V_B' /
       data field_names( 1, 7, 7) / '7    ANALOG_CONV_P15V_A' /
       data field_names( 1, 7, 8) / '8    ANALOG_CONV_P15V_B' /
       data field_names( 1, 7, 9) / '9    ANALOG_CONV_N15V_A' /
       data field_names( 1, 7,10) / '10   ANALOG_CONV_N15V_B' /

!      IPDU other voltages.

       data field_names( 1, 8, 1) / '1    BIAS_PRE_REG_P25V_A' /
       data field_names( 1, 8, 2) / '2    BIAS_PRE_REG_P25V_B' /
       data field_names( 1, 8, 3) / '3    INT_PS_P28V_A' /
       data field_names( 1, 8, 4) / '4    INT_PS_P28V_B' /
       data field_names( 1, 8, 5) / '5    INT_PS_P15V_A' /
       data field_names( 1, 8, 6) / '6    INT_PS_P15V_B' /
       data field_names( 1, 8, 7) / '7    INT_PS_N15V_A' /
       data field_names( 1, 8, 8) / '8    INT_PS_N15V_B' /
       data field_names( 1, 8, 9) / '9    INT_PS_P5V_A' /
       data field_names( 1, 8,10) / '10   INT_PS_P5V_B' /

!      IPDU currents.

       data field_names( 1, 9, 1) / '1    BIAS_PRE_REG_A' /
       data field_names( 1, 9, 2) / '2    BIAS_PRE_REG_B' /
       data field_names( 1, 9, 3) / '3    ANALOG_CONV_A' /
       data field_names( 1, 9, 4) / '4    ANALOG_CONV_B' /
       data field_names( 1, 9, 5) / '5    DIGITAL_CONV_A' /
       data field_names( 1, 9, 6) / '6    DIGITAL_CONV_B' /
       data field_names( 1, 9, 7) / '7    CON_CURRENT_RH' /
       data field_names( 1, 9, 8) / '8    CON_CURRENT_LH' /
       data field_names( 1, 9, 9) / '9    CON_CURRENT_RL' /
       data field_names( 1, 9,10) / '10   CON_CURRENT_LL' /
       data field_names( 1, 9,11) / '11   CON_INT_CONV_A' /
       data field_names( 1, 9,12) / '12   CON_INT_CONV_B' /

!      LMAC temperatures.

       data field_names( 1,10, 1) / '1    LMAC_ANALOG_TEMP' /
       data field_names( 1,10, 2) / '2    LMAC_DIGITAL_TEMP' /
c
c Other variables of interest
c
       integer  *2      ct_status(20)       ! Status array for COBETRIEVE

       Integer	*4	status              ! general status return
       Integer	*4	irow                ! rows for your terminal
       Integer	*4	icol                ! columns for your terminal
       Integer	*4	pbid                ! ID for pasteboard
       Integer	*4	main_vdid           ! ID for main menu
       Integer  *4      vdid(2, num_menus)  ! ID for other menus
       Integer	*4	kbid                ! ID for virtual keyboard
       integer  *4      i,j,k,l             ! Loop counters
       integer  *4      menu_choice         ! menu number user selected
       integer  *4      offset              ! which extrema is being updated
       integer  *4      exit_num            ! 99 on all menus means EXIT
       integer  *4      sub_menu            ! which menu was selected from main
       integer  *4      curr_menu           ! current menu for updating
       integer  *4      frd_status          ! program status reported at end
       integer  *4      inp_len             ! length of input string
       integer  *4      rpt_lun             ! LUN for report
       integer  *4      fxt_lun             ! LUN to read extrema file
       integer  *4      lines               ! lines written in report
       integer  *4      pages               ! Pages written in report
       integer  *4      bin_jstart(2)       ! binary JSTART time
       integer  *4      bin_jstop(2)        ! binary JSTOP  time
       integer  *4      today_adt(2)        ! Todays date in ADT format
      
       integer  *4      cct_get_filespec_fields
       Integer	*4	OTS$Cvt_TI_L
       Integer	*4	OTS$Cvt_T_F
       integer  *4      SMG$Create_PasteBoard
       Integer	*4	SMG$Create_Virtual_Display
       Integer	*4	SMG$Create_Virtual_Keyboard
       Integer	*4	SMG$Delete_Virtual_Keyboard
       Integer	*4	SMG$Paste_Virtual_Display
       Integer	*4	SMG$Unpaste_Virtual_Display
       Integer	*4	SMG$Draw_Rectangle
       Integer	*4	SMG$Draw_Line
       Integer	*4	SMG$Put_Chars
       Integer	*4	SMG$Ring_Bell
       Integer	*4	SMG$Read_String
       Integer	*4	SMG$Set_Cursor_Abs
       Integer	*4	SMG$Set_Cursor_Rel
       Integer	*4	SMG$Erase_PasteBoard
       Integer	*4	SMG$erase_chars
       Integer	*4	SMG$Delete_Virtual_Display
       Integer	*4	SMG$Delete_Pasteboard
       integer  *4      lib$get_lun
       integer  *4      lib$free_lun
       integer  *4      cct_set_ttg_time_range
       Integer  *4      Cut_Register_Version
       Integer  *4      Cut_Display_Banner

       logical  *1      exit                  ! Was exit selected from menu ?
       logical  *1      display_main          ! Should main menu be on screen ?
       logical  *1      display_submenu       ! Should sub-menu be on screen ?
       logical  *1      found_name            ! Found name of updated field
       logical  *1      got_date              ! Got a good JSTART date       
       logical  *1      valid_entry           ! Was entry a number?

       character  *1    date_ok               ! input date OK (Y/N)
       character  *2    char_choice           ! main menu choice string input 
       character  *3    ext                   ! File extension of old fxt file
       character  *6    update_type(2)        ! 'MINIMA' or 'MAXIMA'
       character  *14   jstart                ! jstart value of new extrema file
       character  *14   old_gmt               ! Old update time
       character  *14   old_st                ! jstart value of old extrema file
       character  *14   default_jstart        ! jstart default
       character  *14   jstop                 ! jstop  default
       character  *15   field_value           ! update value for extremum
       character  *30   today                 ! today's date for report
       character  *40   main_head             ! header line for main menu
       character  *65   menu_line             ! a line on the menu
       character  *128  rpt_name              ! name of report file
       character  *128  fxt_name              ! name of extrema file

       Character  *5    version               ! release number
       Parameter       (version='8.7')

       real*8           real_value            ! update value for extremum

       External         FRD_Normal
       external         frd_aberr
       external         ct_init
       external         ct_connect_read
       external         ct_connect_write
       external         ct_read_arcv
       external         ct_write_arcv
       external         cct_set_ttg_time_range
       external         ct_gmt_to_binary
       external         ct_binary_to_gmt

       dictionary 'fxt_eng_xtrm'
       record / fxt_eng_xtrm / fxt_rec

c
c*********************** Begin code ***************************************
c
c
c Initialize variables. 
c
       frd_status      = %loc(frd_normal)
       exit            = .false.
       minmax          = 1
       display_main    = .true.
       display_submenu = .false.
       exit_num        = 99       
       update_type(1)  = 'MINIMA'
       update_type(2)  = 'MAXIMA'
       call sys$gettim(today_adt)   ! Get today's date in ADT format.
c
c Display banner.
c
       status = Cut_Register_Version(version)
       status = Cut_Display_Banner(6,80,'FIRAS Facility FRD_Extrema')
       If ( .not. status ) Then
            frd_status = status
	    Call Lib$Signal(%val(frd_status)) 
            goto 9999
       Endif
c
c Initize COBETRIEVE.
c
       call ct_init(ct_status)
       if (ct_status(1) .ne. ctp_normal) then
           frd_status = ct_status(1)
           call lib$signal( %val(ct_status(1)) )
           goto 9999
       endif
c
c Prompt for start time of file
c
       default_jstart = '00000000000000'
       jstop          = '99365235959999'
       old_st         = default_jstart
       got_date       = .false.
       date_ok        = ' '
       write(6,fmt='(x,/////)')
       do while (.not. got_date)
          write(6,fmt='(1x,a42,$)') 'Enter extrema file JSTART (YYDDDHHMMSS) : '
          read(5,fmt='(q,a14)') inp_len,old_st(1:inp_len)
          old_st(inp_len+1:14) = default_jstart(inp_len+1:14)       
          write(6,fmt='(x,a8,a14,a17,$)') 'Date is ',old_st,' Correct? (Y/N): '
          read(5,fmt='(a1)') date_ok
          if ((date_ok .eq. 'y') .or. (date_ok .eq. 'Y') ) got_date = .true. 
       enddo
c
c Open archive for extrema file
c
       call lib$get_lun(fxt_lun)

       fxt_name='csdr$firas_archive:fxt_eng_xtrm/' // old_st // ';'//jstop//';' 

       open (fxt_lun,file=fxt_name,status='old',
     .               iostat=status,useropen=ct_connect_read)
	If (status .Ne. 0) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If
c
c Read extrema file
c
       call ct_read_arcv(, fxt_lun, fxt_rec, ct_status)
       If (ct_status(1) .ne. ctp_normal) Then
           frd_status = status
           Call LIB$Signal(%val(frd_status))
           goto 9999
       End If
c
c Get file extension of old FXT file to use for New one
c
       status = cct_get_filespec_fields(fxt_lun, specs)
       If (.not. status) Then
            frd_status = status
            Call LIB$Signal(%val(frd_status))
            goto 9999
       End If
 
       ext = specs.filename_extension(1:3)
c
c Close archive for old extrema file
c
       call ct_close_arcv(, fxt_lun, ct_status)
       If (ct_status(1) .ne. ctp_normal) Then
           frd_status = status
           Call LIB$Signal(%val(frd_status))
           goto 9999
       End If

       call lib$free_lun(fxt_lun)
c
c Set up text and titles of maximum sub-menus from the minimum sub-menu values.
c
       do i = 1, num_menus
          menu_titles(2,i) = menu_titles(1,i)
          do j = 1, num_fields(i)
             field_names(2,i,j) = field_names(1,i,j)
          end do
       end do
c
c Save present values of extrema file.
c
       do i = 1 , 2

c         GRTs side A low current (menu 1)
 
          offset = 0
          k      = 1
          do j = 1,16
             old_values(i,k,j) = fxt_rec.minmax(i).grt(offset+j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).grt_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         GRTs side A high current (memu 2)
 
          offset = 16
          k      = 2
          do j = 1,16
             old_values(i,k,j) = fxt_rec.minmax(i).grt(offset+j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).grt_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         GRTs side B low current  (menu 3)
 
          offset = 32
          k      = 3
          do j = 1,16
             old_values(i,k,j) = fxt_rec.minmax(i).grt(offset+j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).grt_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         GRTs side B high current  (menu 4)
 
          offset = 48
          k      = 4
          do j = 1,16
             old_values(i,k,j) = fxt_rec.minmax(i).grt(offset+j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).grt_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         Misc temps and currents  (menu 5)

          offset = 0
          k      = 5
          do j = 1,16
             old_values(i,k,j) = fxt_rec.minmax(i).t_and_i(offset+j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).t_and_i_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         Bolometer voltages  (menu 6)

          offset = 0
          k      = 6
          do j = 1,4
             old_values(i,k,j) = fxt_rec.minmax(i).v_and_i(j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).v_and_i_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         IPDU digital and analog voltages  (menu 7)

          offset = 4
          k      = 7
          do j = 1,10
             old_values(i,k,j) = fxt_rec.minmax(i).v_and_i(j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).v_and_i_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         IPDU other voltages  (menu 8)

          offset = 14
          k      = 8
          do j = 1,10
             old_values(i,k,j) = fxt_rec.minmax(i).v_and_i(j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).v_and_i_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         IPDU currents  (menu 9)

          offset = 24
          k      = 9
          do j = 1,12
             old_values(i,k,j) = fxt_rec.minmax(i).v_and_i(j)
             call ct_binary_to_gmt(fxt_rec.minmax(i).v_and_i_time(1,offset+j),
     .                             old_gmts(i,k,j) )
          end do

c         LMAC temps  (menu 10)

          old_values(i,10,1) = fxt_rec.minmax(i).lmac_analog_temp
          old_values(i,10,2) = fxt_rec.minmax(i).lmac_digital_temp
          call ct_binary_to_gmt(fxt_rec.minmax(i).lmac_time(1,1),
     .                          old_gmts(i,10,1) )
          call ct_binary_to_gmt(fxt_rec.minmax(i).lmac_time(1,2),
     .                          old_gmts(i,10,2) )

       end do
c
c Create pasteboard and virtual displays (menu screens).
c Pasteboard is used as background on which screens can be displayed.
c Also create a virtual keyboard for input.
c
	status = SMG$Create_PasteBoard ( pbid, 'SYS$OUTPUT', irow, icol, )
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If

	status = SMG$Create_Virtual_Keyboard ( kbid )
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
           goto 9999
	End If

	status = SMG$Create_Virtual_Display ( irow, icol, main_vdid )
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
           goto 9999
	End If


        do i = 1,2
           do j = 1, num_menus

              status = SMG$Create_Virtual_Display ( irow, icol, vdid(i,j) )
              If (status .Ne. SS$_Normal) Then
	          frd_status = status
	          Call LIB$Signal(%val(frd_status))
                  goto 9999
	      End If

           end do
        end do
c
c Set up main menu with text.
c
        main_head = 'Main Menu            UPDATE = ' // update_type(minmax)

	status = SMG$Put_Chars ( main_vdid, main_head, 2, 25 )
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If

	status = SMG$Draw_Line ( main_vdid, 3, 1, 3, icol )
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If

	Do i=1,num_menus+3

	   status = SMG$Put_Chars ( main_vdid, MAIN_menu(i), i+3, 10 )
	   If (status .Ne. SS$_Normal) Then
	      frd_status = status
	      Call LIB$Signal(%val(frd_status))
	      goto 9999
	   End If

	End Do

	status = SMG$Draw_rectangle ( main_vdid, 1, 1, irow, icol)
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If

c
c Set up other menus with text.
c
        do i = 1, 2
           do j = 1, num_menus

           menu_titles(i,j)(41:59) = '  UPDATE = ' // update_type(i)

	   status = SMG$Put_Chars ( vdid(i,j), menu_titles(i,j), 2, 10 )
	   If (status .Ne. SS$_Normal) Then
	       frd_status = status
	       Call LIB$Signal(%val(frd_status))
	       goto 9999
	   End If

	   status = SMG$Draw_Line ( vdid(i,j), 3, 1, 3, icol )
	   If (status .Ne. SS$_Normal) Then
	       frd_status = status
	       Call LIB$Signal(%val(frd_status))
	       goto 9999
	   End If

	   Do k = 1,num_fields(j)+1 

              if ( k .lt. (num_fields(j)+1) ) then
                   write(field_value,fmt='(e14.5)',iostat=status,err=50)
     .                   old_values(i,j,k)
50                 continue
                   field_names(i,j,k)(40:55) = field_value
              else
                   field_names(i,j,k) = ret_line
              endif

	      status = SMG$Put_Chars ( vdid(i,j), field_names(i,j,k), k+3, 10 )
	      If (status .Ne. SS$_Normal) Then
	          frd_status = status
	          Call LIB$Signal(%val(frd_status))
	          goto 9999
	      End If

   	    End Do   !  (set up menu text)

	    status = SMG$Draw_rectangle ( vdid(i,j), 1, 1, irow, icol)
	    If (status .Ne. SS$_Normal) Then
                frd_status = status
                Call LIB$Signal(%val(frd_status))
                goto 9999
            End If

	    status = SMG$Set_Cursor_Abs ( vdid(i,j), irow-2, 3 )
	    If (status .Ne. SS$_Normal) Then
	        frd_status = status
	        Call LIB$Signal(%val(frd_status))
	        goto 9999
	    End If

        end do       !  (j = 1,num_fields {loop for sub menus})
        end do       !  (i = 1,2          {loop for min and max})

c
c Display main menu and prompt for choice.
c
        status = smg$paste_virtual_display( main_vdid, pbid, 1, 1 )
	If (status .Ne. SS$_Normal) Then
	    frd_status = status
	    Call LIB$Signal(%val(frd_status))
	    goto 9999
	End If

	Do While (display_main .and. (.not. exit) )

           main_head = 'Main Menu            UPDATE = ' // update_type(minmax)

	   status = SMG$Put_Chars ( main_vdid, main_head, 2, 25 )
	   If (status .Ne. SS$_Normal) Then
	       frd_status = status
	       Call LIB$Signal(%val(frd_status))
	       goto 9999
	   End If

           char_choice = ' '
           sub_menu  = 0
           curr_menu = 0

	   status = SMG$Set_Cursor_Abs ( main_vdid, irow-2, 3 )
	   If (status .Ne. SS$_Normal) Then
	       frd_status = status
	       Call LIB$Signal(%val(frd_status))
	       goto 9999
	   End If

	   status = SMG$Read_String ( kbid, char_choice, 'Select a menu item: ',
     .				      ,,,,,, main_vdid )
	   If (status .Ne. SS$_Normal) Then
	      frd_status = status
	      Call LIB$Signal(%val(frd_status))
	      goto 9999
	   End If

	   status = SMG$erase_chars ( main_vdid, icol-25, irow-2, 23 )
	   If (status .Ne. SS$_Normal) Then
	       frd_status = status
	       Call LIB$Signal(%val(frd_status))
	       goto 9999
	   End If

           status = ots$cvt_ti_l(char_choice,sub_menu,%val(4),%val(17))
           if (.not. status) then
	       status = SMG$Ring_Bell ( main_vdid, 3 )
               valid_entry = .false.
           else
               valid_entry = .true.
           endif
c
c Figure out what action is appropriate from entry.
c Since the first two selections are to choose maxes or mins to update,
c we have to subtract two from the selection number to get correct submenu. 
c
           if (valid_entry) then
               if (sub_menu .eq. exit_num) then
                   exit = .true.    
	       else if ((sub_menu .lt. 1) .Or. (sub_menu .gt. num_menus+2))Then
	           status = SMG$Ring_Bell ( main_vdid, 3 )
	       else if (sub_menu .gt. 2) then
	           display_main    = .False.
                   display_submenu = .true.
                   sub_menu = sub_menu - 2
                   curr_menu = vdid(minmax,sub_menu)
               else if (sub_menu .le. 2) then
                   minmax = sub_menu
	       End If
           endif
c
c Got a sub-menu choice; now display that menu for updating
c
           if (display_submenu .and. .not. exit) then
	       status = SMG$paste_virtual_display ( curr_menu, pbid, 1, 1 )
	       If (status .Ne. SS$_Normal) Then
	           frd_status = status
	           Call LIB$Signal(%val(frd_status))
	           goto 9999
	       End If
           endif
c
c Prompt for an extrema field to update
c
           do while ( display_submenu .and. (.not. exit) )

              char_choice = '  '
              menu_choice = 0

	      status = SMG$Set_Cursor_Abs ( curr_menu, irow-2, 3 )
	      If (status .Ne. SS$_Normal) Then
	          frd_status = status
	          Call LIB$Signal(%val(frd_status))
	          goto 9999
	      End If

	      status = SMG$Read_String(kbid,char_choice,'Select a menu item: ',
     .	                                           ,,,,,, curr_menu )
	      If (status .Ne. SS$_Normal) Then
	          frd_status = status
	          Call LIB$Signal(%val(frd_status))
	          goto 9999
	      End If

	      status = SMG$erase_chars ( curr_menu, icol-25, irow-2, 23 )
	      If (status .Ne. SS$_Normal) Then
	          frd_status = status
	          Call LIB$Signal(%val(frd_status))
	          goto 9999
	      End If

              status = ots$cvt_ti_l(char_choice,menu_choice,%val(4),%val(17))
              if (.not. status) then
                  status = SMG$Ring_bell( curr_menu, 3 )
                  valid_entry = .false.
              else
                  valid_entry = .true.
              endif

              if (valid_entry) then
                  if ( menu_choice .eq. exit_num ) then
                       display_submenu = .false.
                       display_main    = .true.
	          else If ((menu_choice .lt. 1) .Or. 
     .                     (menu_choice .gt. num_fields(sub_menu))) Then
	               status = SMG$Ring_Bell ( curr_menu, 3 )
                  endif
              endif
c
c If choice is to exit sub-menu, go back to main menu; otherwise,
c get the new value for the extremum, then update the screen.
c
              if (menu_choice .eq. exit_num) then
                  status = SMG$unpaste_virtual_display( curr_menu,pbid )
                  If (status .Ne. SS$_Normal) Then
                      frd_status = status
                      Call LIB$Signal(%val(frd_status))
                      goto 9999
 	           End If
              endif

              if ( display_submenu .and. (menu_choice .ne. exit_num) .and. 
     .            (menu_choice .ge. 1) .and. (valid_entry) .and. 
     .            (menu_choice .le. num_fields(sub_menu) ) ) then

                      field_value = '              '
                      real_value  = 0.0
                      valid_entry = .false.
c
c Get a new value for the chosen field.
c
                      do while (.not. valid_entry)
                      
	                 status = SMG$Set_Cursor_Abs ( curr_menu, irow-1, 3 )
	                 If (status .Ne. SS$_Normal) Then
	                     frd_status = status
	                     Call LIB$Signal(%val(frd_status))
	                     goto 9999
	                  End If

	                  status = SMG$Read_String (kbid,field_value,
     .                                              'Enter  a new value: ',
     .		                                    ,,,,,, curr_menu )
	                  If (status .Ne. SS$_Normal) Then
	                      frd_status = status
	                      Call LIB$Signal(%val(frd_status))
	                      goto 9999
	                  End If
c
c Convert character input into a real number;
c Write back into ascii string in e14.5 format
c
                          status=ots$cvt_t_f(field_value,real_value,,,%val(23))
                          if (status) write(field_value,fmt='(e14.5)',
     .                        iostat=status,err=100) real_value
100                       continue

                          if (status .ne. 0) then
                              status = SMG$Ring_Bell( curr_menu, 3 )
                              valid_entry = .false.
                          else
                              valid_entry = .true.
                          endif
                 
                          if (valid_entry) then

                              field_names(minmax,sub_menu,menu_choice)(40:55) = 
     .                        field_value
          
	                      status = SMG$Put_Chars(curr_menu, 
     .                                 field_names(minmax,sub_menu,menu_choice),
     .                                 menu_choice+3, 10)
	                      If (status .Ne. SS$_Normal) Then
	                          frd_status = status
	                          Call LIB$Signal(%val(frd_status))
	                          goto 9999
	                      End If

	                      status=SMG$erase_chars(curr_menu,icol-25,
     .                                               irow-1,23)
	                      If (status .Ne. SS$_Normal) Then
	                          frd_status = status
	                          Call LIB$Signal(%val(frd_status))
	                          goto 9999
	                      End If
c
c Set update flag for value and save new value.
c
                          update_flags(minmax,sub_menu,menu_choice) = .true.
                          update_values(minmax,sub_menu,menu_choice)=real_value
                        
                          endif  ! (valid_entry)

                     enddo       ! (while .not. valid_entry)

                  endif   !  ( display_submenu .and. menu_choice.ne.exit_num)

               end do     !  ( while display_submenu .and. .not. exit)
            
	End Do            !  (while display_main .and. .not. exit)

c
c Delete pasteboard , keyboard, and virtual displays.
c
	status = SMG$Delete_Virtual_Display ( main_vdid )
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If

        do i = 1, 2
           do j = 1, num_menus
	      status = SMG$Delete_Virtual_Display ( vdid(i,j) )
   	      If (status .Ne. SS$_Normal) Then
	          frd_status = status
	          Call LIB$Signal(%val(frd_status))
	          goto 9999
	      End If
            end do
        end do

	status = SMG$Delete_Virtual_Keyboard ( kbid)
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If
       
	status = SMG$Delete_Pasteboard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If
c
c Get date for start time tag of new extrema file.
c
       default_jstart = '00000000000000'
       jstart         = default_jstart
       jstop          = '99365235959999'
       got_date       = .false.
       date_ok        = ' '
       write(6,*) ' '
       write(6,*) 'If you enter a start date for the new extrema file'
       write(6,*) 'that is before today, FXT will have to be rerun from'
       write(6,*) 'the entered date to supercede existing extrema files.'
       write(6,*) ' '
       do while (.not. got_date)
          write(6,fmt='(1x,a42,$)') 'Enter JSTART for NEW file (YYDDDHHMMSS) : '
          read(5,fmt='(q,a14)') inp_len,jstart(1:inp_len)
          jstart(inp_len+1:14) = default_jstart(inp_len+1:14)       
          write(6,fmt='(x,a8,a14,a17,$)') 'Date is ',jstart,' Correct? (Y/N): '
          read(5,fmt='(a1)') date_ok
          if ((date_ok .eq. 'y') .or. (date_ok .eq. 'Y') ) got_date = .true. 
       enddo
c
c Open archive for new extrema file.
c
       call lib$get_lun(fxt_lun)

       fxt_name = 'csdr$firas_archive:fxt_eng_xtrm.' // ext // 
     .             jstart(1:7) // '_' // jstop(1:7) 

       open (fxt_lun,file=fxt_name,status='new',
     .               iostat=status,useropen=ct_connect_write)
       If (status .Ne. 0) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
       End If
c
c Update the start time in the COBETRIEVE header. 
c
       fxt_rec.ct_head.gmt = jstart
       call ct_gmt_to_binary(jstart, bin_jstart)
       call ct_gmt_to_binary(jstop , bin_jstop )
       fxt_rec.ct_head.time(1) = bin_jstart(1)
       fxt_rec.ct_head.time(2) = bin_jstart(2)
c
c Update extrema record with new values. 
c
        do i = 1, 2
           do j = 1, num_menus
              do k = 1, num_fields(j)
                 if (update_flags(i,j,k)) then
                     if (j .le. 4) then              
                         l = ((j-1) * 16) + k 
                         fxt_rec.minmax(i).grt(l) = update_values(i,j,k)
                         fxt_rec.minmax(i).grt_time(l,1) = today_adt(1)
                         fxt_rec.minmax(i).grt_time(l,2) = today_adt(2)
                     else if ( j .eq. 5 ) then
                         l = k
                         fxt_rec.minmax(i).t_and_i(l) = update_values(i,j,k)
                         fxt_rec.minmax(i).t_and_i_time(l,1) = today_adt(1)
                         fxt_rec.minmax(i).t_and_i_time(l,2) = today_adt(2)
                     else if ( j .eq. 6 ) then
                         l = k
                         fxt_rec.minmax(i).v_and_i(l) = update_values(i,j,k)
                         fxt_rec.minmax(i).v_and_i_time(l,1) = today_adt(1)
                         fxt_rec.minmax(i).v_and_i_time(l,2) = today_adt(2)
                     else if ( j .eq. 7 ) then
                         l = 4 + k
                         fxt_rec.minmax(i).v_and_i(l) = update_values(i,j,k)
                         fxt_rec.minmax(i).v_and_i_time(l,1) = today_adt(1)
                         fxt_rec.minmax(i).v_and_i_time(l,2) = today_adt(2)
                     else if ( j .eq. 8 ) then
                         l = 14 + k
                         fxt_rec.minmax(i).v_and_i(l) = update_values(i,j,k)
                         fxt_rec.minmax(i).v_and_i_time(l,1) = today_adt(1)
                         fxt_rec.minmax(i).v_and_i_time(l,2) = today_adt(2)
                     else if ( j .eq. 9 ) then
                         l = 24 + k
                         fxt_rec.minmax(i).v_and_i(l) = update_values(i,j,k)
                         fxt_rec.minmax(i).v_and_i_time(l,1) = today_adt(1)
                         fxt_rec.minmax(i).v_and_i_time(l,2) = today_adt(2)
                      else if ( j .eq. 10) then
                         l = k
                         fxt_rec.minmax(i).lmacs(l) = update_values(i,j,k)
                         fxt_rec.minmax(i).lmac_time(l,1) = today_adt(1)
                         fxt_rec.minmax(i).lmac_time(l,2) = today_adt(2)

                      endif    ! ( j .eq. 1)  {cases for menus 1 - 10}
                  endif        ! ( update_flag(i,j,k) {field was updated}
              enddo            ! ( k = 1,num_fields )
           enddo               ! ( j = 1,num_menus )
        enddo                  ! ( i = 1,2)
c
c Write the updated extrema record to the new extrema file.
c
       call ct_write_arcv(, fxt_lun, fxt_rec, ct_status)
       If (ct_status(1) .ne. ctp_normal) Then
	   frd_status = ct_status(1)
	   Call LIB$Signal(%Val(frd_status))
	   goto 9999
	End If
c
c Set the time tag for the new extrema file.
c
       status = cct_set_ttg_time_range(fxt_lun, bin_jstart, bin_jstop)       
       If (.not. status ) Then
	   frd_status = status
	   Call LIB$Signal(%Val(frd_status))
	   goto 9999
       End If
c
c Close the archive.
c
       call ct_close_arcv(, fxt_lun, ct_status)
       If (ct_status(1) .ne. ctp_normal) Then
           frd_status = status
           Call LIB$Signal(%val(frd_status))
           goto 9999
       End If

       call lib$free_lun(fxt_lun)
c
c Report the updated values
c
        status = lib$get_lun(rpt_lun)
	If (status .Ne. SS$_Normal) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If

        rpt_name = 'frd_extrema_report.lis'

        open(rpt_lun,file=rpt_name,iostat=status,status='new')
	If (status .Ne. 0) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If
c
c Get today's date in ASCII format and write page header.
c Format 810 is for first header without page feed;
c Format 820 does a page feed.
c
        call lib$date_time(today)
        pages = 1

        write(rpt_lun,810,iostat=status) 
     .        pages,old_st,jstart,today,update_type(1)
	If (status .Ne. 0) Then
	   frd_status = status
	   Call LIB$Signal(%val(frd_status))
	   goto 9999
	End If

810     format(' ',t24,'REPORT OF UPDATED EXTREMA VALUES',t68,'Page ',i6,//,
     .     t7,'Old start time: ',t23,a14,t42,'New start time: ',t59,a14,/,
     .     t7,'Run time: ',a30,'    Type: ',a6,//, 
     .     ' Field Name',t34,'Old value',t47,'New value',t59,'Old Time',/)

820     format('1',t24,'REPORT OF UPDATED EXTREMA VALUES',t68,'Page ',i6,//,
     .     t7,'Old start time: ',t23,a14,t42,'New start time: ',t59,a14,/,
     .     t7,'Run time: ',a30,'    Type: ',a6,//, 
     .     ' Field Name',t34,'Old value',t47,'New value',t59,'Old Time',/)

        lines = 7
        pages = 2

        do i = 1,2

           if (i .eq. 2) lines = 60    ! New page for start of maxima.
c
c Write out the title for the menu.
c
           do j = 1, num_menus         

              if (lines .ge. 55) then

                  write(rpt_lun,820,iostat=status)
     .                  pages,old_st,jstart,today,update_type(i) 
	          If (status .Ne. 0) Then
	              frd_status = status
	              Call LIB$Signal(%val(frd_status))
	              goto 9999
	          End If

                  lines = 7
                  pages = pages + 1
              endif

              write(rpt_lun,fmt='(x,/,t20,a34,/)',iostat=status) 
     .              menu_titles(1,j)(1:34)
	      If (status .Ne. 0) Then
	          frd_status = status
	          Call LIB$Signal(%val(frd_status))
	          goto 9999
	      End If

              lines = lines + 3
c
c Write out fields if updated.
c
              do k = 1, num_fields(j)  

                 if (update_flags(i,j,k)) then   

                     if (lines .ge. 55) then
                         write(rpt_lun,820,iostat=status)
     .                         pages,old_st,jstart,today,update_type(i)
	                 If (status .Ne. 0) Then
	                     frd_status = status
	                     Call LIB$Signal(%val(frd_status))
	                     goto 9999
	                 End If
                         lines = 7
                         pages = pages + 1
                     endif                  

                     write(rpt_lun,830,iostat=status) field_names(i,j,k)(6:36),
     .                                                old_values(i,j,k),
     .                                                update_values(i,j,k),
     .                                                old_gmts(i,j,k)
830                  format(' ',a30,t33,e12.5,t46,e12.5,t59,a14)

                     lines = lines + 1

	             If (status .Ne. 0) Then
	                 frd_status = status
	                 Call LIB$Signal(%val(frd_status))
	                 goto 9999
	             End If

                  endif   !  ( update_flags(i,j,k) )
               enddo      !  ( k = 1,num_fields)
           enddo          !  ( j=1,fxt_fields )
        enddo             !  ( i=1,2 )

        close(rpt_lun)
        call lib$free_lun(rpt_lun)

9999    continue
c
C Signal the processing status and exit.
c
        If (frd_status) Then
            Call Lib$Signal(FRD_Normal)
        Else
            Call Lib$Signal(FRD_Aberr)
        Endif

       End
