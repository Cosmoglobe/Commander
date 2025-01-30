      integer*4 function FEC_Analyze (Report_Flag,
     .                                Report_All_Flag,
     .                                LU_Report,
     .                                CT_Unit,
     .                                LU_Short_Sci,
     .                                Data_Array,
     .                                Data_Array_Size,
     .                                Plat_Array,
     .                                Plat_Array_Size,
     .                                Window_Size,
     .                                Tol_type,
     .                                Temp_Size,
     .                                Min_DiHed,
     .                                Max_DiHed,
     .				      Start_Chan,
     .				      Stop_Chan,
     .				      Skip_Chan,
     .				      Write_Flag,
     .				      Instr_Qual,
     .				      Attit_Qual)
c
c     By Ken Jensen, STX, 25-FEB-1991
c
c     Revision of FEC_ANALYZE_STABILITY
c     By Reid Wilson, SASC Technologies Inc,  12-AUG-1986
c
c     Revision made to distinguish final requirements version of FEC
c     from the older version.
c
c     Description:
c         This function examines the plateaus previously found and checks
c         each for stability.  Stable plateaus are output to LU_Short_Sci,
c         the stable plateau records file.  A line to summarize why a plateau
c         is not stable or how long a plateau is stable is output to LU_Report,
c         the report file.
c
c     Format:
c         Ret_Code = FEC_Analyze (Report_Flag,
c                                 Report_All_Flag,
c                                 LU_Report,
c                                 CT_Unit,
c                                 Data_Array,
c                                 Data_Array_Size,
c                                 Plat_Array,
c                                 Plat_Array_Size,
c                                 Window_Size,
c                                 Tol_type,
c                                 Temp_Size,
c                                 Min_DiHed,
c                                 Max_DiHed,
c                                 Start_Chan,
c                                 Stop_Chan,
c                                 Skip_Chan,
c                                 Write_flag,
c                                 Instr_Qual,
c                                 Attit_Qual)
c
c     Input Parameters:
c
c         Report_Flag                = flag to report presence of REPORT
c                                      integer*4
c         Report_All_Flag            = flag to report null and hot spot plateaus
c                                      integer*4
c         LU_Report                  = logical unit connected to report file
c                                      integer*2
c         CT_Unit                    = logical unit connected to eng file
c                                      integer*2
c         LU_Short_Sci               = logical unit connected to short science
c                                      files  integer*4 LU_Short_Sci(4)
c         Data_Array                 = array of engineering records
c                                      RECORD /FDQ_ENG_DATA/ Data_Array(*)
c         Data_Array_Size            = number of elements actually in Data_Array
c                                      integer*4
c         Window_Size                = Number of points in lookahead window
c                                      integer*4 Window_Size
c         Tol_Type                   = Specifies Absolute or Relative tolerance
c                                      integer * 4
c         Temp_Size                  = Amount (Absolute or Relative) that
c                                      temperature may deviate from the mean
c                                      real*4 Temp_Size(4)
c         Min_DiHed                  = Minimum allowable Dihedral Temperature
c                                      real*4
c         Max_DiHed                  = Maximum allowable Dihedral Temperature
c                                      real*4
c	  Start_Chan		     = Start channel
c                                      integer*4
c	  Stop_Chan		     = Stop channel
c                                      integer*4
c	  Skip_Chan		     = Increment channel
c                                      integer*4
c	  Write_flag		     = To write or not to write to the short
c				       science archive; Integer*4
c	  Instr_Qual		     = Threshold of acceptable data quality
c                                      integer*4
c	  Attit_Qual		     = Threshold of acceptable attitude quality
c                                      integer*4
c
c     Input/Output Parameters:
c         Plat_Array                 = array containing rightmost points of
c                                      plateaus found earlier.  all non-stable
c                                      plateaus are replaced by FEC_UNSTABLE,
c                                      while good plateaus contain the number
c                                      of good records in the plateau. these
c                                      numbers are actually placed one element
c                                      before where the plateau is found.
c                                      integer*4 Plat_Array(0:MAX_PLATEAUS)
c         Plat_Array_Size            = number of elements actually in Plat_Array
c                                      integer*4
c 
c     Output Parameters:
c         Ret_Code                   = status of operation
c                                      integer*4
c
c     Modules Called:
c         FEC_SORT                   = routine to sort an array
c         FEC_Report_Plateau         = writes plateau report
c         FEC_Write_Records          = write index file: all stable points
c         FUT_Temp_List              = calculate temperatures from eng. record
c
c     Include Files:
c	  FUT_PARAMS		     = FIRAS parameters
c         FEC_INC                    = parameters used by FEC
c         FEC_MSG                    = message definition for FEC
c         FUT_ERROR
c
c     Changes:
c
c	SER 2379, This program was modified to refer to the new
c		FIRAS archive logical names (csdr$firas_in, _out,
c		_raw, _ref, _uref, and _cal). J.T.Bonnell,
c		December 1, 1988, STX@GSFC
c
c	SPR 2716, Build output file extensions for FEC_SSCAL datasets
c		that reflect the time of the data being analyzed.
c		R. Kummerer, January 5, 1989.
c
c	SPR 2980, Fetch reference data using GET_CONFIG.  R. Kummerer,
c		January 5, 1989.
c
c	SPR 3098, Always returned a successful status even on a failure.
c	R. Kummerer, January 6, 1989.
c
c	SPR 2890, FUT GET_CONFIG related changes to FUT_TEMPERATURE_LIST.
c	R. Kummerer, January 10, 1989.
c
c       SPR 2944, add report_flag to pass presence/absence of REPORT qualifier.
c       D. Bouler, Jan 20, 1989.
c
c	SER 3299, Select data by commandable data quality.
c	R. Kummerer, March 8, 1989.
c
c       SPR 3972, Replace LIB$GET_LUN / LIB$FREE_LUN with analogous FUT
c	routines.  R. Kummerer, June 6, 1989.
c
c       SPR 4051, Error in calculating mean temperatures of plateau if plateau
c	has an odd number of points.  R. Kummerer, June 20, 1989.
c
c       SER 6851, SPR 6721
c         FEC incorrectly interpreting Hot Spot heater commands
c         H. Wang, June 13, 1990, STX
c
c       SPR 9428, MAX_DIHED qualifier added to argument list, and passed
c         to subroutine FEC_WRITE_RECORDS.
c         K. Jensen, STX, 01/23/92
c
c       SPR 9429, The declared sizes of the reference data sets FEX_GRTRAWWT
c         and FEX_GRTTRANS needed correction.
c         K. Jensen, STX, 01/23/92
c
c       SPR 12270, MIN_DIHED qualifier added to argument list, and passed
c         to subroutine FEC_WRITE_RECORDS.
c         K. Jensen, Hughes STX, 10/20/95
c
      implicit none
c
c     Include files:
c
      include '(FUT_PARAMS)'
      include '(FEC_INC)'
      include '(FEC_MSG)'
      include '(FUT_ERROR)'
      include '($SSDEF)'
      include 'ct$library:ctuser.inc'
      include '(cct_filespec_fields_record)'
      include '(cct_get_config)'

      record /filespec_fields/ file_spec
c
c     Formal Parameters:
c
      dictionary 'FDQ_ENG'
      record /FDQ_ENG/           Data_Array(*)
      integer*4                  LU_Report,
     .                           Report_Flag,
     .                           Report_All_Flag,
     .                           LU_Short_Sci(4),
     .				 CT_Unit
      integer*4                  Data_Array_Size,
     .                           Plat_Array_Size,
     .                           Plat_Array(0:MAX_PLATEAUS),
     .                           Window_Size,
     .				 Tol_type,
     .				 Write_flag,
     .				 Start_Chan,
     .				 Stop_Chan,
     .				 Skip_Chan,
     .				 Instr_Qual,
     .				 Attit_Qual,
     .				 First_Good_Point
      real*4                     Temp_Size(Num_Temps_Used)
      real*4                     Min_DiHed
      real*4                     Max_DiHed
      character*64		 out_file
      integer*4                  status
      integer*4                  offgrts
c
c     External Routines:
c
      integer*4 FEC_Write_Records,
     .          FEC_Report_Plateau,
     .          FUT_Temp_List,
     .          FEC_Sort,
     .		CCT_Open_Config,
     .		CCT_Get_Config_TOD,
     .		CCT_Close_Config,
     .		FUT_Get_LUN,
     .		SYS$GetTim,
     .		CCT_Get_FileSpec_Fields,
     .		CT_Connect_Write
      external  CT_Connect_Write

c
c     Local Variables:
c
      real*4    Temp_Array(FEC_MAX_AVERAGE_SIZE, FAC_Struct),
     .          Temp(10),
     .          SigTemp(10),
     .          Sort_Array(FEC_MAX_AVERAGE_SIZE),
     .          Mean(Num_Temps_Used),
     .          Mean_A(Num_Temps_Used),
     .          Mean_B(Num_Temps_Used)
      integer*4 loop1,
     .          loop2,
     .          loop3,
     .          Ret_Code,
     .          Stable_counter,
     .          Plat_Counter,
     .          Report_Counter,
     .          Counter,
     .          Num_in_Window,
     .          Num_Good,
     .          Num_Bad_Tol,
     .          Num_Bad_Qual(4),
     .          Min_No_of_IFGS,
     .          LU_CC_Thresholds,
     .          CC_Min,
     .          Min_Numbers(4), !as defined in include file in FCI
     .          Min,
     .		Temp_Count,
     .          Temp_Flag
      logical*1 		 Stable
      logical*1 		 First_Time
      logical*4                  Singlifg/.TRUE./
      character*32		 file_ext
      character*14 		 eff_start
      character*14 		 eff_stop

c
c     GET_CONFIG variables.
c
      dictionary 'fex_mincoadd'
      dictionary 'fex_grtrawwt'
      dictionary 'fex_grttrans'
      structure /CONFIG_DATA/
         record /fex_mincoadd/ fex_mincoadd
         record /fex_grtrawwt/ fex_grtrawwt
         record /fex_grttrans/ fex_grttrans
      endstructure

      record /CONFIG_DATA/ CONFIG
      record /config_status/ stat(3)

      character*1		 access_mode/' '/	! data set access mode
      integer*4			 number/3/		! number of data sets
      integer*4			 size(3)/128,256,384/	! size of data sets

      character*32		 name(3)		! names of data sets
      data name(1)/'csdr$firas_ref:fex_mincoadd'/
      data name(2)/'csdr$firas_ref:fex_grtrawwt'/
      data name(3)/'csdr$firas_ref:fex_grttrans'/

      integer*4			 lun(3)			! logical unit numbers
      integer*4			 index(3)		! initial cache pointers
      logical*1			 new_segment(3)		! flag for new segments
      integer*4			 ncache/1/
      integer*4			 ref_count

      dictionary 'FUT_EngAnlg'
      record /FUT_EngAnlg/       Dummy_Sigs
c
c     *** BEGIN ***
c
      Plat_Array(0) = 0
      Plat_Counter  = 0

c
c Get consistency thresholds using FEX_MINCOADD configuration file.
c
      Ret_Code = CCT_Open_Config ( Data_Array(1).CT_Head.Time,
     .				   Data_Array(Data_Array_Size).CT_Head.Time,
     .				   number, name, size, access_mode,
     .				   ncache, lun, index, stat, ref_count )
      if (.not. Ret_Code) then
	 call lib$signal(fec_opnconfigerr,%val(1),%val(Ret_Code))
         FEC_Analyze = Ret_Code
         return
      endif

C
C Open output files
C
      If (Write_flag .Eq. 1) Then

           Ret_Code = cct_get_filespec_fields(CT_Unit,file_spec)

           if (.not. Ret_Code) then
	         call lib$signal(%val(Ret_Code))
                 FEC_Analyze = Ret_Code
                 return
           endif
       
	   call ct_binary_to_gmt ( data_array(1).ct_head.time(1), eff_start )
	   call ct_binary_to_gmt ( data_array(data_array_size).ct_head.time(1),
     .					eff_stop )

 	   file_ext = '.' // file_spec.filename_extension(1:3) //
     .				eff_start(1:7) // '_' // eff_stop(1:7)

	   Do Counter = Start_Chan, Stop_Chan, Skip_Chan

	      Ret_Code = fut_get_lun(LU_Short_Sci(Counter))
	      if (.Not. Ret_Code)then
	         call lib$signal(%val(Ret_Code))
                 FEC_Analyze = Ret_Code
                 return
	      end if

	      out_file = 'CSDR$FIRAS_OUT:fec_sscal_' //
     .				fac_channel_ids(Counter) //
     .				file_ext
	      open(unit=LU_Short_Sci(Counter), file=out_file, status='new',
     .			iostat=status, useropen=ct_connect_write)

	      if (status .ne. 0) then
	         Ret_Code = %loc(fec_ctopenss)
		 call lib$signal(fec_ctopenss,%val(1),%val(counter))
                 FEC_Analyze = Ret_Code
                 return
	      end if

	   End Do

      End If

c
c Find the calibration plateaus.
c
      report_counter = 0
      loop1 = 1 
      do while (loop1 .le. Plat_Array_Size)
C        
        if( Ret_Code) then
         if ( plat_array(loop1) .eq. 0 .and. report_flag .eq. 1) then
c
c               If /REPORTALL, Then Report Hot Spot Plateau
c
            if( report_all_flag .eq. 1) then

                   Ret_Code = FEC_Report_Plateau (
     .                            LU_Report,
     .                            -1,
     .			          0,
     .                            Data_Array(Plat_Array(loop1-1)+1),
     .                            Data_Array(Plat_Array(loop1+1)),
     .                            Data_Array(Plat_Array(loop1-1)+1),
     .                            Data_Array(Plat_Array(loop1+1)),
     .                            Plat_Array(loop1+1)-Plat_Array(loop1-1),
     .                            0,
     .                            Num_Bad_Qual,
     .                            Plat_Array(loop1+1)-Plat_Array(loop1-1),
     .                            Tol_Type,
     .                            Temp_Size,
     .                            Mean_A,
     .                            Mean_B)
                   if (.Not. Ret_Code)then
	              call lib$signal(%val(Ret_Code))
                      FEC_Analyze = Ret_Code
                      return
	           end if
            end if
         end if
         if ( plat_array(loop1) .ne. 0) then
          Min_No_Of_IFGS = Plat_Array(loop1) - Plat_Array(loop1-1)

c
c	Get consistency thresholds reference data using GET_CONFIG.
c
	  Ret_Code = CCT_Get_Config_TOD ( Data_Array(Plat_Array(loop1-1)+1).CT_Head.Time,
     .					  number, size, lun, index,
     .					  CONFIG, new_segment, stat )

	  if (.not. Ret_Code) then
             call lib$signal(fec_getconfigerr,%val(1),%val(Ret_Code))
             FEC_Analyze = Ret_Code
             return
	  end if

	  if (new_segment(1)) then
	     do loop2 = 1, 4
	        min_numbers(loop2) = CONFIG.FEX_MINCOADD.Min_IFG_Coadd(loop2)
	     end do
             CC_Min = Min ( Min_Numbers(1), Min_Numbers(2),
     .				Min_Numbers(3), Min_Numbers(4) )
	  end if

          if(Min_No_of_IFGS .lt. CC_Min) then
                do loop2 = 1, Num_Temps_Used
                   Mean_A(loop2) = 0.0
                   Mean_B(loop2) = 0.0
                end do              
                do loop2 = 1, 4
                   Num_Bad_Qual(loop2) = 0
                end do              
         	if( Ret_Code .and. (report_flag .eq. 1) )then
c
c               If /REPORTALL, Then Report Null Plateau
c
                   If (report_all_flag .eq. 1) Then
                      Ret_Code = FEC_Report_Plateau (
     .                               LU_Report,
     .                               0,
     .	         		     0,
     .                               Data_Array(Plat_Array(loop1-1)+1),
     .                               Data_Array(Plat_Array(loop1)),
     .                               Data_Array(Plat_Array(loop1-1)+1),
     .                               Data_Array(Plat_Array(loop1)),
     .                               Plat_Array(loop1)-Plat_Array(loop1-1),
     .                               0,
     .                               Num_Bad_Qual,
     .                               Plat_Array(loop1)-Plat_Array(loop1-1),
     .                               Tol_Type,
     .                               Temp_Size,
     .                               Mean_A,
     .                               Mean_B)
                      if (.Not. Ret_Code)then
	                 call lib$signal(%val(Ret_Code))
                         FEC_Analyze = Ret_Code
                         return
                      end if
	           end if
                end if
                Plat_Array(loop1 - 1) = FEC_UNSTABLE  ! Indicates unstable plateau  
          else
c
c         else if Min_No_of_IFGs . ge. CC_Min
c
	      First_Time = .TRUE.
              if(Min_No_of_IFGS .gt. FEC_MAX_AVERAGE_SIZE) 
     -            Min_No_of_IFGS = FEC_MAX_AVERAGE_SIZE

c             Calculate the mean temperatures ( A side )
c
              do loop2 = 0, Min_No_of_IFGS-1
		  If(Ret_Code)
     .            Offgrts = FUT_Temp_List(
     .                  Data_Array(Plat_Array(loop1)-loop2).EN_ANALOG,
     .                  Dummy_Sigs,
     .			CONFIG.FEX_GRTRAWWT,
     .			CONFIG.FEX_GRTTRANS,
     .                  1,Singlifg,Temp,SigTemp)
                  do loop3 = 1, Num_Temps_Used
                      Temp_Array(loop2+1, loop3) = Temp(loop3)
                  end do
              end do
	      Temp_Count = Min_No_of_IFGS - 2 * (Min_No_of_IFGS / 4)
              do loop2 = 1, Num_Temps_Used
                  Mean(loop2) = 0.0
                  do loop3 = 1, Min_No_of_IFGS
                      Sort_Array(loop3) = Temp_Array(loop3, loop2)
                  end do
		  If(Ret_Code)then
                     Ret_Code = FEC_Sort (Sort_Array, Min_No_of_IFGS)
                     if (.Not. Ret_Code)then
	                call lib$signal(%val(Ret_Code))
                        FEC_Analyze = Ret_Code
                        return
	             end if
                  End if
                  Temp_Flag = 0
                  do loop3 = (Min_no_of_IFGS/4)+1,
     .                       Min_no_of_IFGS-(Min_No_of_IFGS/4)
                      Mean(loop2) = Mean(loop2) + Sort_Array(loop3)
                      if(Sort_Array(loop3) .eq. -9999.)Temp_Flag = 1
                  end do
                  if(Temp_Flag .eq. 0)then
                     Mean_A(loop2) = Mean(loop2) / Temp_Count
                  else
                     Mean_A(loop2) = 0.0
                  endif
              end do
c
c             Calculate the mean temperatures ( B side )
c
              do loop2 = 0, Min_No_of_IFGS-1
		  If(Ret_Code)
     .            Offgrts = FUT_Temp_List(
     .                  Data_Array(Plat_Array(loop1)-loop2).EN_ANALOG,
     .                  Dummy_Sigs,
     .			CONFIG.FEX_GRTRAWWT,
     .			CONFIG.FEX_GRTTRANS,
     .                  2,Singlifg,Temp,SigTemp)
                  do loop3 = 1, Num_Temps_Used
                      Temp_Array(loop2+1, loop3) = Temp(loop3)
                  end do
              end do
	      Temp_Count = Min_No_of_IFGS - 2 * (Min_No_of_IFGS / 4)
              do loop2 = 1, Num_Temps_Used
                  Mean(loop2) = 0.0
                  do loop3 = 1, Min_No_of_IFGS
                      Sort_Array(loop3) = Temp_Array(loop3, loop2)
                  end do
		  If(Ret_Code)Then
                     Ret_Code = FEC_Sort (Sort_Array, Min_No_of_IFGS)
                     if (.Not. Ret_Code)then
	                call lib$signal(%val(Ret_Code))
                        FEC_Analyze = Ret_Code
                        return
	             end if
                  Endif
                  Temp_Flag = 0
                  do loop3 = (Min_no_of_IFGS/4)+1,
     .                       Min_no_of_IFGS-(Min_No_of_IFGS/4)
                      Mean(loop2) = Mean(loop2) + Sort_Array(loop3)
                      if(Sort_Array(loop3) .eq. -9999.)Temp_Flag = 1
                  end do
                  if(Temp_Flag .eq. 0)then
                     Mean_B(loop2) = Mean(loop2) / Temp_Count
                  else
                     Mean_B(loop2) = 0.0
                  endif
              end do
c
c         Check for stability: Min_No_of_IFGS points in a row, all
c         GRTs ( A side and B side ) within Tolerance of the Mean.   
c
              Stable = .TRUE.              
              Stable_Counter = 0
              Num_Bad_Tol = 0
              loop2 = Plat_Array(loop1)
              do while (Stable .and. (loop2.gt.Plat_Array(loop1-1)))

c     A side

		  If(Ret_Code)
     .            Offgrts = FUT_Temp_List(
     .                                    Data_Array(loop2).EN_ANALOG,
     .                                    Dummy_Sigs,
     .					  CONFIG.FEX_GRTRAWWT,
     .                                    CONFIG.FEX_GRTTRANS,
     .                                    1,Singlifg,Temp,SigTemp)
                  do loop3 = 1, Num_Temps_Used
                   if (mean_A(loop3) .gt. 0.0 .and. Stable) then
                    if (tol_type .eq. 1) then 
                       if(abs(Mean_A(loop3)-Temp(loop3)) .gt. Temp_Size(loop3)) then
                          Stable = .FALSE.
                       end if
                    elseif(tol_type .eq. 2) then
                       if((abs(Mean_A(loop3)-Temp(loop3)))/mean_A(loop3)
	1	           .gt. Temp_Size(loop3)) then
                           Stable = .FALSE.
                       end if
                    endif
                   endif 
                  enddo

c     B side   (If A side is stable)

                  If( Stable )Then
                   If(Ret_Code)
     .             Offgrts = FUT_Temp_List(
     .                                     Data_Array(loop2).EN_ANALOG,
     .                                     Dummy_Sigs,
     .	                                   CONFIG.FEX_GRTRAWWT,
     .                                     CONFIG.FEX_GRTTRANS,
     .                                     2,Singlifg,Temp,SigTemp)
                   do loop3 = 1, Num_Temps_Used
                    if (Mean_B(loop3) .gt. 0.0) then
                     if (tol_type .eq. 1) then 
                        if(abs(Mean_B(loop3)-Temp(loop3)) .gt. Temp_Size(loop3)) then
                           Stable = .FALSE.
                        end if
                     elseif(tol_type .eq. 2) then
                        if((abs(Mean_B(loop3)-Temp(loop3)))/mean_B(loop3)
	1	            .gt. Temp_Size(loop3)) then
                            Stable = .FALSE.
                        end if
                     endif
                    endif
                   end do
                  Endif

                  if(Stable) then
                      Stable_Counter = Stable_Counter + 1
		      if(First_Time) then
			  First_Good_Point = loop2
			  First_Time = .FALSE.
		      end if
c
                  else
c
c                     Use window to look ahead for good data and determine whether
c                     to continue or not.
c
                      Num_in_Window = 1
                      Num_Good = 0
                      Num_Bad_Tol = Num_Bad_Tol + 1
                      do while ((Num_in_Window .le. Window_Size) .and.
     .                   ((loop2 - Num_in_Window) .gt. Plat_Array(loop1 - 1)) .and.
     .                   Ret_Code)
                         Offgrts = FUT_Temp_List(
     .                            Data_Array(loop2-Num_in_Window).EN_ANALOG,
     .                            Dummy_Sigs,
     .                            CONFIG.FEX_GRTRAWWT,
     .                            CONFIG.FEX_GRTTRANS,
     .                            1,Singlifg,Temp,SigTemp)
                         if( Ret_Code ) then
                             Stable = .TRUE.
                             do loop3 = 1, Num_Temps_Used
                               if (Mean_A(loop3) .gt. 0.0) then
                                if (tol_type .eq. 1) then
                                   if(abs(Mean_A(loop3)-Temp(loop3)) .gt.
     .                                Temp_Size(loop3))  then
                                      Stable = .FALSE.
                                   end if
                                elseif(tol_type .eq. 2) then
                                   if((abs(Mean_A(loop3)-Temp(loop3)))
	1	                   /mean_A(loop3) .gt. Temp_Size(loop3)) then
                                      Stable = .FALSE.
                                   end if
                                end if
                               end if
                             end do
                         endif 
  
                         Offgrts = FUT_Temp_List(
     .                            Data_Array(loop2-Num_in_Window).EN_ANALOG,
     .                            Dummy_Sigs,
     .                            CONFIG.FEX_GRTRAWWT,
     .                            CONFIG.FEX_GRTTRANS,
     .                            2,Singlifg,Temp,SigTemp)
                         if( Ret_Code ) then
                             do loop3 = 1, Num_Temps_Used
                               if (Mean_B(loop3) .gt. 0.0) then
                                if (tol_type .eq. 1) then
                                   if(abs(Mean_B(loop3)-Temp(loop3)) .gt.
     .                                Temp_Size(loop3))  then
                                      Stable = .FALSE.
                                   end if
                                elseif(tol_type .eq. 2) then
                                   if((abs(Mean_B(loop3)-Temp(loop3)))
	1	                   /mean_B(loop3) .gt. Temp_Size(loop3)) then
                                      Stable = .FALSE.
                                   end if
                                end if
                               end if
                             end do
                             if ( Stable ) then
                                Num_Good = Num_Good + 1
                             endif
                         endif 
                         Num_in_Window = Num_in_Window + 1
                      end do

c
                      if( Ret_Code .and. (Num_Good .ge. int((Window_Size+1)/2))) then
                          Stable = .TRUE.
c
c                 Set flag in record to mark a bad data point
c
                          Data_Array(loop2).CT_HEAD.TIME(1) = 0
                          Data_Array(loop2).CT_HEAD.TIME(2) = 0
                      else
                          Stable = .FALSE.
		          Num_Bad_Tol = Num_Bad_Tol - 1
                      end if
                  end if
                  if( Stable ) loop2 = loop2 - 1
              end do !do while (Stable .and. (loop2.gt.Plat_Array(loop1-1)))
              if(Stable_Counter.ge.Min_No_of_IFGS) then
                  Plat_Counter = Plat_Counter + 1
c
c     ** Plat_Array(loop1)-Plat_Array(loop1-1)-Stable_Counter-Num_Bad_Tol is the
c        number of points that were stranded because a fatally bad point
c        occurred before the end of the plateau.
c
		  If( Ret_Code)then
                      Report_Counter = Report_Counter + 1
                      Ret_Code = FEC_Write_Records(LU_Short_Sci,
     .                                             Data_Array,
     .                                             Config,
     .                                             loop2+1,
     .                                             First_Good_Point,
     .		                                   Start_Chan,
     .					           Stop_Chan,
     .					           Skip_Chan,
     .					           Write_flag,
     .					           Instr_Qual,
     .					           Attit_Qual,
     .					           Min_DiHed,
     .					           Max_DiHed,
     .					           Num_Bad_Qual)
                      if (.Not. Ret_Code)then
                         call lib$signal(%val(Ret_Code))
                         FEC_Analyze = Ret_Code
                         return
                      end if
	          End if

		  If( Ret_Code .and. (report_flag .eq. 1))then
c
c                     Report Short Science Plateau
c
                      Ret_Code = FEC_Report_Plateau (
     .                                   LU_Report,
     .                                   Report_Counter,
     .					 Stable_Counter,
     .                                   Data_Array(Plat_Array(loop1-1)+1),
     .                                   Data_Array(Plat_Array(loop1)),
     .                                   Data_Array(loop2+1),
     .                                   Data_Array(First_Good_Point),
     .                                   Plat_Array(loop1)-Plat_Array(loop1-1),
     .                                   Num_Bad_Tol,
     .                                   Num_Bad_Qual,
     .                                   Plat_Array(loop1)-Plat_Array(loop1-1)
     _                                       -Stable_Counter-Num_Bad_Tol,
     .                                   tol_type,
     .                                   Temp_Size,
     .                                   Mean_A,
     .                                   Mean_B)
                      if (.Not. Ret_Code)then
	                 call lib$signal(%val(Ret_Code))
                         FEC_Analyze = Ret_Code
                         return
                      end if
	          End if


                  Plat_Array(loop1 - 1) = Stable_Counter !Number of stable
                                                         !points in plateau
              else

                  do loop3 = 1, 4
                     Num_Bad_Qual(loop3) = 0
                  end do              

		  If( Ret_Code .and. (report_flag .eq. 1))then
c
c                     If /REPORTALL, Then Report Null Plateau
c
                    If ( report_all_flag .eq. 1) then
                      Ret_Code = FEC_Report_Plateau (
     .                                   LU_Report,
     .                                   0,
     .					 Stable_Counter,
     .                                   Data_Array(Plat_Array(loop1-1)+1),
     .                                   Data_Array(Plat_Array(loop1)),
     .                                   Data_Array(Plat_Array(loop1-1)+1),
     .                                   Data_Array(Plat_Array(loop1)),
     .                                   Plat_Array(loop1)-Plat_Array(loop1-1),
     .                                   Num_Bad_Tol,
     .                                   Num_Bad_Qual,
     .                                   Plat_Array(loop1)-Plat_Array(loop1-1)
     _                                       -Stable_Counter-Num_Bad_Tol,
     .                                   tol_type,
     .                                   Temp_Size,
     .                                   Mean_A,
     .                                   Mean_B)

                      if (.Not. Ret_Code)then
	                  call lib$signal(%val(Ret_Code))
                          FEC_Analyze = Ret_Code
                          return
                      end if
	            End if
                  End if

                  Plat_Array(loop1 - 1) = FEC_UNSTABLE  ! Indicates unstable plateau  

              end if !if(Stable_Counter.ge.Min_No_of_IFGs) then 

          end if !if(Min_No_of_IFGs .lt. CC_Min) then
          loop1 = loop1 + 1
         else
           loop1 = loop1 + 2
         endif !if( Plat_array(loop1) .ne. 0) then
        end if !if( Ret_Code ) then
      end do  !do while (loop1 .le. Plat_Array_Size)

c
c Close configuration file.
c
      status = CCT_Close_Config ( number, lun, index )
      if (.not. status) then
         Ret_Code=%loc(fec_clsconfigerr)
	 call lib$signal(fec_clsconfigerr,%val(1),%val(status))
         FEC_Analyze = Ret_Code
         return
      endif

      FEC_Analyze = Ret_Code

      Return
      End
