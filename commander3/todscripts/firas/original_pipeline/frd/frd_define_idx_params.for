	program frd_define_idx_params
c/
c/	program name:
c/	  frd_define_idx_params
c/
c/	program description:
c/	  This program allows the user to set or change the tolerances
c/	  for various class of items in the FIRAS Index File.  These
c/	  are the tolerances over which DATA QUALifY will write out a
c/	  new Index record.
c/	  It also allows the user to set or change the index flags for
c/	  the same class of items.
c/
c/	author:
c/	  E. Fung
c/	  GSFC
c/	  April 9, 1986
c/
c/	modified by:
c/	  J. Durachta
c/        STX
c/        December 31, 1986
c/        Reason: Capability to create and alter index flags file added.
c/
c/	  J. Durachta
c/        ARC
c/        Nov 18,1987
c/        Reason: Remove READONLY from open statements. FDQ_TEXT_ITEM included.
c/                Program made more user friendly.
c/
c/	  Shirley M. Read
c/        STX
c/        June 2, 1988
c/        Reason:  Moved program to the FRD Facility; prefixed name with FRD.
c/	           To create the Index Params files in the new format to 
c/	           be archived as a time tagged file under COBETRIEVE. The
c/	           COBETRIEVE type header was added and a database structure
c/	           was defined for this dataset. The new files will be 
c/		   named Fex_Idx_Tols.Dat and Fex_Idx_Flag.Dat.
c/
c/----------------------------------------------------------------------------
c/	  Quoc C Chung
c/        STX
c/        April 6, 1989
c/        Reason:  fixed FRD Facility of SPR 3048 ; FRD_Define_idx_params
c/                 still read obsolete files.
c/----------------------------------------------------------------------------
c/        version 4.4 : Quoc C Chung, STX, May 9, 1989
c/        SPR 3773  Frd_Define_Idx_Params writes report to FOR000.DAT.
c/
c/        SPR 9047  Add qualifier INDEX with a string values of TOLS or FLAGS
c/                  to generate reference files. If /INDEX=TOLS, then program
c/                  output will be FEX_IDX_TOLS, ElseIf /INDEX=FLAGS, then
c/                  output will be FEX_IDX_FLAGS. If program runs without
c/                  qualifiers, then both reference files will be created plus
c/                  a report files.
c/                  Nilo G. Gonzales/STX, September 19, 1991.
c/----------------------------------------------------------------------------
c/
c/	calling sequence:
c/	  This is a main program
c/
c/	input/output parameters:
c/	  n/a
c/
c/      Input files :
c/       Fex_Idx_Tols.txt
c/       Fex_Idx_FlagS.txt
c/
c/	Output files:
c/	  Fex_Idx_Tols.Dat
c/        Fex_Idx_Flag.Dat
c/
c/	subroutines called:
c/ 	  Lib$Get_Lun (from system library)
c/ 	  Sys$Gettim (from system library)
c/	  CT_Binary_To_GMT
c/
c/	include files used:
c/	  none
c/
c/	error handling:
c/	  tbd
c/
c/	method used:
c/	  PDL --
c/
c/        If qualifier INDEX is present, then
c/           If qualifier is equal 'TOLS', then
c/              set Tol flag to true
c/           Else If qualifier is equal 'FLAGS', then
c/              set Flags flag to true
c/           Endif
c/         Else 
c/             set flag All to true
c/             set flag Tol to true
c/             set flag Flags to true
c/         End If
c/
c/         If qualifier REPORT is present, then
c/            set report flag to true
c/         Else If flag (All) is true
c/               set report flag to true
c/         Else
c/               set report flag to false
c/         End If
c/
c/	  Inquire if tolerance file exists
c/	  if not then
c/	    Open tolerance file as a new file;
c/	  else
c/	    Open tolerance file as an old file;
c/	    Read in old tolerance values into buffer;
c/	  endif;
c/
c/	  Display "old" values;
c/	  Inquire if changes needed;
c/	  if yes then
c/	    Ask user to supply new value(s);
c/	    Write updated value(s) to tolerance file;
c/	  endif
c/	  

	Implicit	None

!	Include files

	Include '($SSdef)'
        INCLUDE '(UPM_STAT_MSG)'

	Dictionary	'Fex_Idx_Tols'
	Dictionary	'Fex_Idx_Flag'

	Record		/Fex_Idx_Tols/Tol_Rec
	Record		/Fex_Idx_Flag/Flag_Rec

!	Functions

	Integer*4 Lib$Get_Lun    !Get unit number (from system library)
	Integer*4 Lib$Free_Lun   !Free unit number (from system library)
	Integer*4 Sys$Gettim     !Get time
        Integer*4 Frd_Parse_Idx_Tols  ! routine to process Idx Tols
        Integer*4 Frd_parse_Idx_flags ! routine to process Idx flags
        Integer*4 Upm_present
	Integer*4 Upm_get_value	
	Integer*4 status 	! Returned system status
	Integer*4 retstat       ! Status of processing
	Integer*4 success / 1 /, error / 2 /
	Integer*4 curtime(2)	! Current time quadword
	Character*14 curgmt	! Ascii GMT
	Integer*4 ix            ! Index

	logical*1	file_there1,file_there2, blocks,report
        logical*1       idx_flags(282)/282*.true./, more/.true./
        Logical*1       Tol/.false./,flags/.false./,all/.false./

        character*1     flag, ans
        Character*15    file
        Character*30    rpt_file1, rpt_file2
	character*5     index
	character*5     oqual
	integer*2       len
        integer*4       upm_stat

	integer*2	choice, i, tol_val, tols(256) / 256*0 /
        integer*2       upper, lower, num

	integer*4	ist, prt_idx_flags, num_false/0/
        integer*4       rpt_lun1, rpt_lun2
	integer*4       lun1, newlun1, lun2, newlun2
        integer*4       type
        Integer*4       zero/0/

	parameter	bol_rdout_tol	= 1
	parameter	grt_tol		= 2
	parameter	tc_tol		= 3
	parameter	analog_tol	= 4	! This may be sub-divided later
	parameter	nclass		= 22	! # of classes

        External        Frd_abort
        External        Frd_normal
        External        Frd_getluerr
        External        Frd_filemiss
        External        Frd_rmsopen
        External        Frd_badtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Set status to success.
	
	retstat = success
c
c Parse the command line.
c
        if (upm_present ('Index') .eq. upm_pres) then
	   upm_stat = upm_get_value('index',oqual,len)
	   if (upm_stat .ne. ss$_normal) then
	      call lib$signal (%val(ss$_abort))
	   end if
	   if (oqual .eq. 'TOLS') then
               Tol = .true.
	   else if (oqual .eq. 'FLAGS') then
                    Flags = .true.
           end if
	else
	       All = .true.
	       Tol = .true.
	       Flags = .true.
	end if
c
c **** get logical unit for report file
c
        if (upm_present ('REPORT') .eq. upm_pres) then
           report = .true.
        elseif (All) then
	       report = .true.
        else
           report = .false.
	endif  ! UPm_present report

	if ((Tol .and. Report) .or. All) then         
 	   status = Lib$Get_Lun(rpt_lun1)
           if ( Status .ne. SS$_normal) then
             retstat = error
             call lib$signal(frd_getluerr,%val(1),%val(status))
           endif
        endif

	if ((Flags .and. Report) .or. All) then         
 	   status = Lib$Get_Lun(rpt_lun2)
           if ( Status .ne. SS$_normal) then
             retstat = error
             call lib$signal(frd_getluerr,%val(1),%val(status))
           endif
        endif
c
c	Get the current time for CT_Head in record.
c
	status = Sys$Gettim ( curtime )
	if ( status .eq. SS$_Normal ) Then
	   call CT_Binary_to_GMT ( curtime, curgmt )
	   Write(6,60) curgmt
	else
	   retstat = error
           Write(6,50) status
	endif ! Status

        if ((Tol .and. Report) .or. All) then
            Rpt_File1 = 'FEX_idx_tols_rpt.' // curgmt(1:9)
	endif

        if ((Flags .and. Report) .or. All) then
            Rpt_File2 = 'FEX_idx_flag_rpt.' // curgmt(1:9)
        else if (.not. report) then
             print *,' No report file set.'
        endif

	if ((Tol .and. Report) .or. All) then
           open (unit=rpt_lun1, file=Rpt_File1,
	1  type='new',form='formatted',iostat=ist)
           if (ist .ne. zero) then
              call Lib$signal(Frd_rmsopen,%val(2),%val(rpt_lun1),
	1                       %val(ist) )
              retstat= error
           endif

           if (report) write (rpt_lun1,10) curgmt
           if (report) write (rpt_lun1,11) rpt_file1
  10          format(22X,' ***** FEX_IDX_Tols report file at: ',a,' *****'//)
  11	      format(' report file : ',a,' has been opened successfully.'/)
 	endif 
c      
c Get a logical unit for reading and writing files.
c
	if (retstat .eq. success .and. (Tol .or. All)) Then
           status = Lib$get_lun(lun1)
	   if (status .ne. SS$_Normal) Then
	      retstat = error
              call lib$signal(frd_getluerr,%val(1),%val(status))
              if (report) Write(rpt_lun1,40) status
	   endif
        endif
c
	if (retstat .eq. success .and. (Tol .or. All)) Then
           status = lib$get_lun(newlun1)
	   if (status .ne. SS$_Normal) Then
	       retstat = error
               call lib$signal(frd_getluerr,%val(1),%val(status))
               if (report) Write(rpt_lun1,40) status
	   endif
	endif

c  Put current time in the records for output.

	if ( retstat .eq. success .and. (TOL .or. ALL) ) Then
	     tol_rec.ct_head.gmt(1:14) = curgmt(1:14)
	     Do ix = 1, 2
                tol_rec.ct_head.time(ix)= curtime(ix)
             Enddo
	endif     ! Retstat is success

	If (retstat .eq. success .and. (Tol .or. All)) then
c
c *** get which type to process TOLS/Flag
c
           Inquire (file='Fex_idx_tols.txt', 
	1	   exist=file_there1)

	   if (file_there1) then	   
	       open (unit=lun1, file='FEX_idx_tols.txt',
	1	     type='old',form='formatted',recl=80, iostat=ist,
	1	     readonly, shared)
	       if (ist .ne. zero) then
	          write (6, 6001) ist
                  Call lib$signal(FRD_rmsopen,%val(2),%val(lun1),
	1                          %val(ist) )
                  retstat = error
	       endif  ! Ist

               If (Report) write(rpt_lun1,62)
 62	       format(1x,' Fex_Idx_Tols.Txt has been opened successfully.'/)
c
c **** parse the Index Tolerances parameters for changing
c
               if (retstat .eq. success .and. (Tol .or. All)) then
                  status = Frd_parse_idx_tols (tol_rec,report,lun1,rpt_lun1,
	1                                     newlun1)
                  If ( .not. Status) then
                     if (report)  write(rpt_lun1,6002) status
                        Retstat = error
                  Endif !status
               endif  ! Retstat is success
           else
                Call lib$signal(FRD_filemiss)
                Retstat = error
           endif    ! File_there1

	endif      !retstat and Tol or All

	if ((Flags .and. Report) .or. All) then
           open (unit=rpt_lun2, file=Rpt_File2,
	1  type='new',form='formatted',iostat=ist)
           if (ist .ne. zero) then
              call Lib$signal(Frd_rmsopen,%val(2),%val(rpt_lun2),
	1                       %val(ist) )
              retstat= error
           endif

           if (report) write (rpt_lun2,13) curgmt
           if (report) write (rpt_lun2,14) rpt_file2
  13          format(22X,' ***** FEX_IDX_Flags report file at: ',a,' *****'//)
  14	      format(' report file : ',a,' has been opened successfully.'/)
 	endif 
c      
c Get a logical unit for reading and writing files.
c
	if (retstat .eq. success .and. (Flags .or. All)) Then
           status = Lib$get_lun(lun2)
	   if (status .ne. SS$_Normal) Then
	      retstat = error
              call lib$signal(frd_getluerr,%val(1),%val(status))
              if (report) Write(rpt_lun2,40) status
	   endif
        endif
c
	if (retstat .eq. success .and. (Flags .or. All)) Then
           status = lib$get_lun(newlun2)
	   if (status .ne. SS$_Normal) Then
	      retstat = error
              call lib$signal(frd_getluerr,%val(1),%val(status))
              if (report) Write(rpt_lun2,40) status
	   endif
	endif

c  Put current time in the records for output.

	if ( retstat .eq. success .and. (Flags .or. ALL) ) Then
	     flag_rec.ct_head.gmt(1:14) = curgmt(1:14)
	     Do ix = 1, 2
                flag_rec.ct_head.time(ix)= curtime(ix)
             Enddo
	endif     ! Retstat is success
c
c  process idx flags data
c
	if ( retstat .eq. success .and. ( Flags .or. All) ) then 

	   Inquire (file='FEx_idx_flag.txt', 
	1           exist=file_there2)

	   If (file_there2) then
	       open (unit=lun2, file='FEx_idx_flag.txt',
	1      form='formatted', recl=80, iostat=ist,
	1      type='old',readonly, shared)
	       if (ist .ne. 0) then
	          write (6, 6103) ist
                  retstat = error
	       endif
               if( report) write(rpt_lun2,61)
 61	       format(1x,' Fex_Idx_Flag.Txt has been opened successfully.'/)
c
c start read and process the idx flags
c
               If (retstat .eq. success .and. ( Flags .or. All)) then
                  Status = Frd_parse_idx_flags(flag_rec ,lun2 ,report,
 	1                                      rpt_lun2, newlun2)
                  if ( .not. status) then
                     if ( report) write(rpt_lun2,6104) status
                          retstat = error
                  endif
	        endif
	   Else
               Call lib$signal(FRD_filemiss)
               Retstat = error
           Endif    ! File_there2

	Endif  ! retstat and flags .or. all

 40	format(1x,'Error: Failed to get unit number. Status= ',z8.8)
 50     format(1x,'Error: Bad return from Sys$Gettim. Status= ',z8.8)  
 60     format(1x,'Current GMT : ',a)
 400	format(1x,'Error: Closing Fex_Idx_Tols.Dat. Status= ',
	1	z8.8)
 401	format(1x,'Error: Closing Fex_Idx_Tols.txt. Status= ',
	1	z8.8)
 700    format(1x,'Success: Fex_Idx_Tols.Dat Created')
6001	format (/' Cannot open tolerance text file. IOSTAT: ', I4)
6002	format (/' Failed on process/read tolerance text file. STATUS: ',
	1          Z8.8)
 600	format(1x,'Error: Closing Fex_Idx_Flag.Dat. Status= ',
	1	z8.8)
 601	format(1x,'Error: Closing Fex_Idx_Flag.txt. Status= ',
	1	z8.8)
 800    format(1x,'Success: Fex_Idx_Flag.Dat Created')
6103	format (/' Cannot open index flag file. IOSTAT: ', I10)
6104	format (/' Failed on process/read index flag file. STATUS: ', Z8.8)

c
c   exit the program with a status
c
	if ( retstat .eq. error ) then
           Close (unit=newlun1, iostat=status)
           Close (unit=newlun2, iostat=status)
           close(lun1, iostat=ist)
           close(lun2, iostat=ist)
           call lib$signal(frd_abort,%val(1),%val(retstat))
           status = SS$_abort
           call exit(status)
        Elseif (retstat .eq. success) then
c
c close all files before exit
c
        if (Tol .or. All ) then
	   Close (unit=newlun1, iostat=status)
 	   if ( status .ne. 0 ) then
	      if (report) Write(rpt_lun1,400) status
	         retstat = error
           endif
	   close (lun1, iostat=ist)
           if (report) Write(rpt_lun1,700)
	endif  !Tol
c
	if (Flags .or. All ) then
	   Close (unit=newlun2, iostat=status)
 	       if ( status.ne.0 ) then
	          if (report) Write(6,600) status
	             retstat = error
	          else
	             Write(6,800)
	          endif
                  close(lun2, iostat=ist)
 	          if ( ist .ne. 0 ) then
	             if ( report) Write(rpt_lun2,601) status
	          endif
	       endif  ! Flags
               call lib$signal(frd_normal,%val(1),%val(retstat))
               status = SS$_normal
               call exit(status)
        endif
        stop
	end

