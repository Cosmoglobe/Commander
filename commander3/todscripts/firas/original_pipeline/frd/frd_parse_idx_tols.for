
	Integer*4 Function Frd_parse_idx_tols(tol_rec,report,lun,rpt_lun,newlun)
C
C       This routine reads the Fex_idx_tols.txt (RMS file), and
C       get the idx name and values of each record, use the values  
C       to update the Frd_idx_tols.dat .
! 
!	PROGRAM NAME : FRD_Parse_idx_tols
! 
!           
! 	AUTHOR:
! 	  Quoc C Chung
! 	  STX
! 	  April 7 1989
! 
! 	MODIFIED BY:
! 	  N/A
! 
! 	INPUT PARAMETERS:
! 	  NONE
! 
! 	OUTPUT PARAMETERS:
! 	  The record structure buffer Tols
! 
!       INPUT FILES:
! 	  FEX_Idx_Tols.Txt
! 
! 	OUPUT FILES:
! 	  FEX_Idx_Tols.Dat
! 
! 	INCLUDE FILES USED:
!        None
! 
! 	SUBROUTINES USED:
!	  Lib$Skpc
!	  FRD_L functions:
!		FRD_L_Remtab
! 		FRD_Extract_Idx_Tolparms
! 
! 	ERROR HANDLING:
! 	  Passed back in output parameter Status
! 
! 	METHOD USED:
! 
!	Description :   FRD_Parse_idx_tols reads FEX_idx_tols.Txt, containing 
!		      the Field name ,Class Number and the % Value, into 
!		      an 80 character buffer. It calls FRD_Extract_idx_tols to 
!		      extract the Class Number and % values and store
!		      the % values into a record arrays. When the reading
!		      of the file is complete, the % values are moved into the
!		      Fex_idx_tol record.
!	              Finally , write Tol_rec to a binary file (RMS file),
!                     FEX_idx_tols.Dat.
!
!----------------------------------------------------------------------
!       VERSION 4.4  QUOC C CHUNG, STX, MAY 9, 1989
!       SPR 3788  FRD_DEFINE_IDX_PARAMS CAN NOT BE INTERACTIVE.
!---------------------------------------------------------------------
!       VERSION 4.4  QUOC C CHUNG, STX, MAY 10, 1989
!       SPR 3773  FRD_DEFINE_IDX_PARAMS WRITES REPORT TO FOR000.DAT.
!----------------------------------------------------------------------
	Implicit None
!
!	Include Files

	Include '($SSdef)'
!	
!	Pass parameter

	dictionary 'Fex_idx_tols'

        Record    /Fex_Idx_Tols/Tol_Rec

        Integer*4 lun,rpt_lun,newlun
        
!	Functions

	Integer*4 FRD_L_Remtab
	Integer*4 FRD_Extract_idx_tolparms
	Integer*4 Lib$Skpc	     ! Library function to skip characters
	
!	Local variables

	Character*80	buffer 	     ! Input character string
 
	Character*32 grp_names(22)/'Grt_controlled_low_current      ',
     a                             'Grt_controlled_high_current     ',
     b                             'Grt_dihedrals_low_current       ',
     c                             'Grt_dihedrals_high_current      ',
     d                             'Grt_bolometers_low_current      ',
     e                             'Grt_bolometers_high_current     ',
     f                             'Grt_mirrors_low_current         ',
     g                             'Grt_mirrors_high_current        ',
     h                             'Temperature_controllers         ',
     i                             'Other_temperatures              ',
     j                             'Bolometer_voltage_readouts      ',
     k                             'Ipdu_digital_convert_voltages   ',
     l                             'Ipdu_analog_convert_voltages    ',
     m                             'Ipdu_bias_pre_regulator_voltages',
     n                             'Ipdu_28v_internal_pwr_voltages  ',
     o                             'Ipdu_15v_internal_pwr_voltages  ',
     p                             'Ipdu_5v_internal_pwr_voltages   ',
     q                             'Ipdu_bias_pre_regulator_currents',
     r                             'Ipdu_analog_converter_currents  ',
     s                             'Ipdu_digital_converter_currents ',
     t                             'Ipdu_constant_currents          ',
     u                             'Ipdu_internal_converter_currents'/
!
        Character Ans                ! user response
!
	Integer*4 lpos               ! Position on line
	Integer*4 retstat 	     ! Program Processing status
	Integer*4 status 	     ! Returned system status
	Integer*4 systatus           ! I/O status 
	Integer*4 success / 1 /, error / 2 /
        Integer*4 ix
!
	Integer*2 Tols(256)/256*0/   ! percentage values
        Integer*2 perc_number        ! converted percentage number
        Integer*2 number             ! converted Class number
        Integer*2 class_no(22)       ! Group class number
!
	Logical*1 eof                ! End of file
        Logical*1 Found              ! matching idx tols field name
        Logical*1 report

        External Frd_itemnofnd       ! item field not found in the
                                     ! table.
!	Set status to success.

	retstat = success



!	Read the Fex_idx_tols.Txt file and extract the data.

        If ( retstat .eq. success ) then
	  eof = .false.
	  do while (( .not. eof ) .and. (retstat .eq. success ))
            Read (unit=lun,iostat=systatus,end=500,err=4000,fmt=300)
	1	  buffer
	    lpos = lib$skpc(' ',buffer)
	    if ( lpos .ne. 0 ) then		! Line is not blank
	      if ( buffer(1:1) .ne. '*'.and. buffer(1:1) .ne. '!' ) then 
	        status = FRD_l_Remtab ( buffer )
	        if ( status .ne. success ) then
	          retstat = error
	        endif
	        if ( retstat .eq. success ) then
	          status = FRD_Extract_idx_tolparms ( buffer, perc_number,
	1                            rpt_lun,report,grp_names, number , found)
	          if ( status .eq. success .and. found) then
                     If ( number .gt. 0 .or. number .le. 22) then
                       Class_no(number) = number                  
                       Tols(number) = perc_number
                     else
                      retstat = error
                      if (report) write (rpt_lun,10) number
                     endif  ! Number
		  else
	            retstat = error
                    Call lib$signal(frd_itemnofnd)
	          endif  ! Status and Found
	        endif    ! Retstat is success
	      endif	! Buffer is not a comment
	    endif	! Buffer is not blank
	  enddo         ! end do while
       endif            ! Retstat is success

4000       If ( systatus .ne. 0 ) then
	    retstat = error
	    if (report) write(rpt_lun,200) status
	  Endif
500        If ( retstat .eq. success ) then
	    if (report) write(rpt_lun,250) 
	  Endif
!
! Display the Fex_Idx_Tols.txt for verification
!
          Write(6,6098)
          if (report) write(rpt_lun,6098)
          Do ix = 1 , 22
            write(6,6100) Grp_names(ix),Class_no(ix),Tols(ix)
            if (report)write(rpt_lun,6100) Grp_names(ix),Class_no(ix),Tols(ix)
          Enddo
 6098   format(13x,'Field name',13x,'   Choice   percentage'/,80('-'))
 6100	Format(4X,a,5x,i2,5x,i3)

!       Open the output file, Fex_Idx_Tols.Dat 

	  if ( retstat .eq. success ) then

	    Open (unit=newlun, file='Fex_Idx_Tols.Dat',
	1 	  status='new', form='unformatted', 
	2	  organization='sequential', iostat=status,
	3	  recordsize=256, recordtype='fixed')

	    if (status .ne. 0) then
	      if (report) Write(Rpt_lun,100) status
	      retstat = error
	    endif


! 	Move the idx tols from the array to the idx tols Record.

	If ( retstat .eq. success ) then
	  do ix = 1 , 256
	     tol_rec.idx_tols.tols(ix) = tols(ix)
          enddo
	Endif

!	Close Fex_idx_tols.Txt. If error, keep going.

	  close(unit=lun,iostat=systatus)
          if (systatus .ne. 0 ) then
            if (report) write (rpt_lun,275) status
	  endif

	Endif  	! Retstat is success	  

!	If success,  write Fex_idx_tols.Dat.

	If ( retstat .eq. success ) then
             write(unit= newlun,iostat=systatus) tol_rec
	    if ( systatus .ne. 0 ) then
	      if (report) write(rpt_lun,350) 
	      retstat = error
	    else
              if (report) Write(rpt_lun,400)
	    endif
	Endif


!	Set function return status.

	FRD_Parse_idx_tols = retstat

	Return

!	Formats
!
  10	format(1x,' Invalid class number: ',i4)
 100	format(1x,'Error: Opening Fex_Idx_Tols.Dat. Status= ',z8.8)
 200    format(1x,'Error : Reading Fex_idx_tols.Txt. Status= ',z8.8)
 250    format(1x,'Success : Reading of Fex_idx_tols.Txt completed.')
 275    format(1x,'Error : Closing Fex_idx_tols.Txt. Status= ',z8.8)
 300    format(a)
 350    format(1x,'fail: Fex_idx_tols.Dat file written.')
 400    format(1x,'success: Fex_idx_tols.Dat file written.')
	end
