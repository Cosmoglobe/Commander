
	PROGRAM FRD_L_DEFINE_lIMITS
! 
! 	PROGRAM NAME : FRD_L_Define_Limits
! 
! 	PROGRAM DESCRIPTION:
! 	    This subroutine reads the Limits text files (RMS disk files)  
! 	  and loads the red and yellow limits into structured buffers.
! 	  These limits will be used for checking both science and converted 
! 	  housekeeping data and setting quality flags in the science records. 
! 	  The red and yellow (when applicable) limiting values for each 
! 	  entity corresponding to one of 110 science record quality flags will
!         be compared to the values stored in the Raw Science Record, the 
! 	  Housekeeping Record, the Engineering Record or the Index Record.
!         The quality flags will be set in the Science record according to 
! 	  the results of the comparisons. 
! 	    In order to do this, new database record structures were defined.
! 	  The FEX_SCILIM is organized like the NFS_SDF, but omits the 
! 	  structures with fields that are not involved in the quality flag
!         checks. It contains all limits needed for comparison with data in  
! 	  the Science Record. There will be one record, with the science field
!         structure dimensioned by two, in the file. One dimension of the array 
! 	  contains red limits and one contains yellow limits. This file
! 	  will be maintained offline by the Firas Subsystem. FRD_L_DEFINE_LIMITS 
! 	  will read an edit file, FEX_Scilim.Txt, and possibly some
! 	  interactive user input and produce an RMS binary data file, 
! 	  FEX_Scilim.Dat, to be used by FDQ. Eventually the file will be
! 	  archived under Cobetrieve by a utility program. The program will be
!         rerun whenever there is a change in any of the limits.  
!           The FEX_ENGLIM contains all RDL substructures in FDQ_ENG.
! 	  The red and yellow limits for the GRT low and high currents will 
!         be read from an edit file, FEX_Grtlim.Txt, by this program.
! 	  The limits for the remainder of the engineering analogs will be
! 	  obtained from the IGSE database file, DB$:Firlims.DB, which already
! 	  contains the red and yellow limits for these analog database points.
! 	  The binary engineering file produced by the program is FEX_Englim.Dat.
! 	  There will be one record containing a field structure dimensioned by
!         four for the four sets: yellow low, yellow high, red low and red
! 	  high. There may also be interactive user inputs for limits. This 
! 	  program will also be run by the Firas Subsystem whenever there is a
! 	  change in any of the limits. This file will also be put in the 
!	  COBETRIEVE archive. 
! 	    The program also reads the Limits enable/disable file, 
!         FEX_Limflags.Txt, which contains flags to switch each limit check on
!         or off. It will also be put in the COBETRIEVE archive.
!           
! 	AUTHOR:
! 	  Shirley M. Read
! 	  STX
! 	  January 1988
! 
! 	MODIFIED BY:
! 	  Shirley M. Read
! 	  STX
! 	  May 1988
!	    All limits datasets were redesigned to be easily accessable via a
!         new COBETRIEVE routine which gets the datasets with time tags which
!         cover the time of the data. All files now contain only one record
!         with dimensioned arrays of structures. The Limits Flags was completely
!         redesigned to enable limit checking on individual FIRAS instrument 
!         components which may have different readouts depending on which side
!         of the spacecraft is powering the FIRAS components and which side
!         of the FIRAS instrument is powering the components. The GRT input
!	  file was updated to include the engineering record GRT field array
!         name on each input line as well as the position within the array.
!         The current pipeline software does not convert the calibrator 
!	  resistors to engineering units nor does it provide for limits. 
!         This program rejects any input lines with the array positions 
!         corresponding to the calibrator resistors.
! 
! 	INPUT PARAMETERS:
! 	  NONE
! 
! 	OUTPUT PARAMETERS:
! 	  The record structure buffers SCILIM_REC and ENGLIM_REC are filled.
! 	  The record structure buffer LIMFLAGS is filled.
! 
!       INPUT FILES:
! 	  FEX_Scilim.Txt
! 	  DB$:Firlims.Inp
!         FEX_Grtlims.Txt
! 	  FEX_Limflags.Txt
! 
! 	OUPUT FILES:
! 	  FEX_Scilim.Dat
! 	  FEX_Englim.Dat
! 	  FEX_Limflags.Dat
! 
! 	INCLUDE FILES USED:
!        $SSDef
! 
! 	SUBROUTINES USED:
! 	  Lib$Get_Lun (from system library)
! 	  Lib$Free_Lun (from system library)
! 	  Sys$Gettim (from system library)
! 	  Binary_To_GMT
! 	  FRD_L Functions:
!               FRD_L_Process_Sci
!               FRD_L_Process_Eng
!               FRD_L_Process_Flags
!
! 	ERROR HANDLING:
! 	  Passed back in output parameter STATUS
! 
! 	METHOD USED:
! 
! 	FRD_L_DEFINE_LIMITS reads the limits text files, consisting of red and
! 	yellow limits for science/attitude and engineering data, and loads these
! 	limits into record structured buffers. The limits will be used for
! 	setting the 110 data quality flags in the science record. The routine
! 	then reads the limits enable/disable flags corresponding to the 
!       FIRAS instrument components and selected telemetry information for 
!       which data quality flags are defined.
! 
!
!	Description :   FRD_L_Define_Limits reads text files containing 
!		      parameters and their Red and Yellow Limits,
!		      extracts the parm and limits, converts the limits
!	              to real*4 or integer, depending on the type, and
!		      creates binary FEX_SCILIM.DAT and FEX_ENGLIM.DAT
!		      limits files for use in FDQ. It also reads a text 
!		      file containing flags for enabling/disabling limit
!		      checks and creates binary file FEX_LIMFLAGS.DAT. 
!			Each input record for all of the limit text files
!		      is a string of 80 characters. The parameters on the 
!		      string are different for each type of limits text file.
!		        For the science text files, the first parameter is
!		      the RDL field name of one of the science or attitude
!		      fields. After the science or attitude field name two 
!		      integer limits, separated by blanks, must follow on 
!		      the same line : Red and Yellow. 
!		        For the engineering analog text files, the first 
!		      parameter is assumed to be the name of an IGSE 
!		      engineering database word. After the database word are 
!		      the key words 'RANGE' and either 'ENGR' or 'COUNTS', 
!		      meaning engineering units or unconverted counts,
!		      respectively. Four real limits, separated by blanks,
!		      must follow on the next line. The order is Red Low, 
!		      Yellow Low, Yellow High and Red High. The IGSE database
!		      name is used to reference the position of the analog
!		      limit in the analog fields of the FEX_ENGLIM record. 
!		      The matching of position to database name will be done
!		      by using FDQ_Firnames as an include file.
!			The GRT limits file consists of 64 GRT sets of a GRT
!		      field name followed by the array position number and 
!		      then followed by the red and yellow limits in the 
!		      same order as the analog limits. The number corresponds
!		      to the position of the GRT in the GRT record fields
!		      of the FEX_ENGLIM record. 
! 			The individual limit flags for the science type or 
!		      engineering type limits are extracted from the records 
!                     containing the flag name in the Limflags arrays followed
!                     by the value of the flag for single dimension flags and
!                     the ordered values for the arrays. Most limit flags have
!		      a value of 0 or -1 ( false or true ). Others are bit
!                     masks for individual limit checks which set only one
!                     bit each in the flag.
!
!----------------------------------------------------------------------

	Implicit None

!	Passed parameters for functions

	Dictionary	'Fex_Scilim'
	Dictionary	'Fex_Englim'
	Dictionary	'Fex_Limflags'

	Record		/Fex_Scilim/scilim_rec
	Record		/Fex_Englim/englim_rec
	Record	     	/Fex_Limflags/limflags

	Integer*4 lun, lun2     ! Logical unit numbers

!	Include files

	Include '($SSdef)'

!	Functions

	Integer*4 Lib$Get_Lun    !Get unit number (from system library)
	Integer*4 Lib$Free_Lun   !Free unit number (from system library)
	Integer*4 Sys$Gettim
	Integer*4 Upm_Present
	Integer*4 FRD_L_Process_Sci
	Integer*4 FRD_L_Process_Eng
	Integer*4 FRD_L_Process_Flags
	
!	Local variables

	Integer*4 retstat	! Program Processing status
	Integer*4 status 	! Returned system status
	Integer*4 success / 1 /, error / 2 /
	Logical*1 sci /.false./  ! Flags to make science, engineering or 
			         ! limit flag files. Default is false.
	Logical*1 eng /.false./
	Logical*1 flag /.false./
	Integer*4 curtime(2)	! Current time quadword
	Character*14 curgmt	! Ascii GMT
	Integer*4 ix            ! Index

!	Set status to success.
	
	retstat = success

!	Get a logical unit for reading and writing files.

	status = Lib$Get_Lun(lun)
	If ( status .ne. SS$_Normal ) Then
	  retstat = error
	  Write(6,40) status
	Endif
	If ( retstat .eq. success ) Then
	  status = lib$get_lun(lun2)
	  If ( status .ne. SS$_Normal ) Then
	    retstat = error
	    Write(6,40) status
	  Endif
	Endif

!	Get the current time for CT_Head in record.

	If (retstat .eq. success) Then
	    status = Sys$Gettim ( curtime )
	    If ( status .eq. SS$_Normal ) Then
	      call Binary_to_GMT ( curtime, %ref(curgmt))
	    Else
	      retstat = error
              Write(6,50) status
	    Endif
        Endif

!	Get the user input options.

	If (UPM_Present('SCI')) sci = .true.

	If (UPM_Present('ENG')) eng = .true.

	If (UPM_Present('FLAG')) flag = .true.

	Write(6,60) curgmt

!  	Put current time in the records for output.

	  If ( (retstat .eq. success) .and. (sci) ) Then
	      scilim_rec.ct_head.gmt(1:14) = curgmt(1:14)
	      Do ix = 1, 2
                 scilim_rec.ct_head.time(ix)= curtime(ix)
              Enddo
	  Endif     ! Retstat is success

	  If ( (retstat .eq. success) .and. (eng) ) Then
	      englim_rec.ct_head.gmt(1:14) = curgmt(1:14)
	      Do ix = 1, 2
                 englim_rec.ct_head.time(ix)= curtime(ix)
              Enddo
	  Endif     ! Retstat is success

	  If ( (retstat .eq. success) .and. (flag) ) Then
	      limflags.ct_head.gmt(1:14) = curgmt(1:14)
	      Do ix = 1, 2
                 limflags.ct_head.time(ix)= curtime(ix)
              Enddo
	  Endif     ! Retstat is success

!	If user requested a new science limits file, call FRD_L_Process_Sci.

	If ( (retstat .eq. success) .and. (sci)) Then
	  Write(6,100)
	  status = FRD_L_Process_Sci ( Scilim_Rec, Lun, Lun2 )
          If ( status .ne. success ) Then
	    retstat = error
	  Endif
	Endif

!	If user requested a new engineering limits file, call FRD_L_Process_Eng.

	If ( (retstat .eq. success) .and. (eng)) Then
	  Write(6,200)
	  status = FRD_L_Process_Eng ( Englim_Rec, Lun, Lun2 )
	  If ( status .ne. success ) Then
	    retstat = error
	  Endif
	Endif

!	If user requested a new limits flag file, call FRD_L_Process_Flags.

	If ( (retstat .eq. success) .and. (flag)) Then
	  Write(6,300)
	  status = FRD_L_Process_Flags ( Limflags, Lun, Lun2 )
	  If ( status .ne. success ) Then
	    retstat = error
	  Endif
	Endif

	Call Exit(retstat)

!	Formats:
 40	format(1x,'Error: Failed to get unit number. Status= ',z8.8)
 50     format(1x,'Error: Bad return from Sys$Gettim. Status= ',z8.8)  
 60     format(10x,'*****  FRD_L_Define_Limts  *****'//
	1      10x,'Current Time: ',a/ )
 100    format(10x,'Create Science Limits Database '/)     
 200    format(10x,'Create Engineering Limits Database '/)     
 300    format(10x,'Create Limit Flags Database '/)     
	End
