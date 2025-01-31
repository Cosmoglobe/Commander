      integer*4 function FEC_Write_Records(
     .                                     LU_Short_Sci,
     .                                     Array,
     .                                     Config,
     .                                     Leftmost,
     .                                     Rightmost,
     .					   Start_Chan,
     .					   Stop_Chan,
     .					   Skip_Chan,
     .					   Write_Flag,
     .					   Instr_Qual,
     .					   Attit_Qual,
     .					   Min_DiHed,
     .					   Max_DiHed,
     .					   Num_Bad_Qual)
c
c      By Reid Wilson, SASC Technologies Inc.,  22-OCT-1986
c
c      Description:
c            This routine writes all of the engineering records of a stable plateau.
c
c      Format:
c            Ret_Code = FEC_Write_Records(
c                                           LU_Short_Sci,
c                                           Array,
c                                           Config,
c                                           Leftmost,
c                                           Rightmost,
c                                           Start_Chan,
c                                           Stop_Chan,
c                                           Skip_Chan,
c					    Write_Flag,
c     					    Instr_Qual,
c     					    Attit_Qual,
c     					    Min_DiHed,
c     					    Max_DiHed,
c					    Num_Bad_Qual)
c
c      Input Parameters:
c            LU_Short_Sci                = logical units connected to short
c                                          science output files
c                                          integer*4 LU_Short_Sci(4)
c
c            Array                       = array containing engineering records
c                                          record /FDQ_ENG/ Array(*)
c
c            Config                      = array containing configuration data
c                                          record /CONFIG_DATA/ Config
c
c            Leftmost                    = leftmost point of plateau
c                                          integer*4
c
c            Rightmost                   = rightmost point of plateau
c                                          integer*4
c
c            Start_Chan			 = start channel; integer*4
c
c            Stop_Chan			 = stop channel; integer*4
c
c            Skip_Chan			 = increment channel; integer*4
c
c	     Write_Flag			 = write to archive control flag
c                                          integer*4
c
c	     Instr_Qual			 = Instrument quality threshold
c                                          integer*4
c
c	     Attit_Qual			 = Attitude quality threshold
c                                          integer*4
c
c	     Min_DiHed			 = Minimum allowable Dihedral
c                                          Temperature
c                                          real*4
c
c	     Max_DiHed			 = Maximum allowable Dihedral
c                                          Temperature
c                                          real*4
c
c      Output Parameters:
c	     Num_Bad_Qual		 = number of data points rejected
c					   due to bad quality, etc
c                                          integer*4
c            Ret_Code                    = return status of operation
c                                          integer*4
c      
c      Modules Called:
c            none
c
c      Include Files:
c            FEC_INC                    = facility definitions for FEC
c            FEC_MSG                    = contains message definitions 
c
c      Modification History:
c      
c      12/07/87 Reid Wilson
c               Modified to write short science records instead of engineering.
c       2/09/88 Reid Wilson
c               Replace I/O to RMS files with I/O to short science archives
c       2/07/89 R. Kummerer
c		SPR 3249, Keep FEC from destroying time-order of data. Remove
c		call to FEC_SORT_BY_SM.
c       3/08/89 R. Kummerer
c		SPR 3243, Query on instrument modes. Science mode must be in
c		short science record so it can be queried.
c       3/08/89 R. Kummerer
c		SER 3299, Select data by commandable quality.
c       1/23/92 K. Jensen, STX
c               SPR 9428, Add quality check on Dihedral Temperature.
c       8/19/92 K. Jensen, Hughes STX
c               SPR 9903, Break up large plateaus into equal sized plateaus
c       8/19/92 K. Jensen, Hughes STX
c               SPR 9921, Transition flags to mark gain changes
c       10/20/95 K. Jensen, Hughes STX
c                SPR 12270, Add Minimum Dihedral Temperature.
c
c
c     ******* BEGIN *******
c
      implicit none
c
c     Include Files:
c
      dictionary 'FDQ_ENG' !logical name for file containing
                           !structure of engineering record.
      dictionary 'FEC_SSCAL'  !structure of short science record

      dictionary 'FUT_EngAnlg'

      dictionary 'Fex_MinCoadd'

      dictionary 'Fex_GrtRawWt'

      dictionary 'Fex_GrtTrans'

      include '(FEC_INC)'                !contains normal return value,
                                         !values of Tolerance
      include '(FEC_MSG)'                !contains FEC message definitions
      include '(FUT_PARAMS)'             !FIRAS global parameters
      include 'CT$Library:CTUser.Inc'
c
c     Formal Parameters:
c
      RECORD /FDQ_ENG/             Array(*)
      RECORD /FUT_EngAnlg/         Dummy_Sigs
 
      structure /CONFIG_DATA/
         record /fex_mincoadd/ fex_mincoadd
         record /fex_grtrawwt/ fex_grtrawwt
         record /fex_grttrans/ fex_grttrans
      endstructure

      record /CONFIG_DATA/         CONFIG

      integer*4                    LU_Short_Sci(4)
      integer*4                    Leftmost,
     .                             Rightmost

      
c
c     Called Routines:
c
      integer*4 fut_temp_list
c
c     Local Variables:
c
      integer * 2     ct_status(20)     !I/O status for COBETRIEVE read
      logical * 4     singlifg/.TRUE./
      integer*4 
     .      Bad_Qflags,
     .      Eng_Counter,
     .      ios,
     .      Num_Recs,
     .      Num_Group,
     .      Num_Div,
     .      Group_Size,
     .      ix,
     .      ix0,
     .      Ret_Code,
     .      SS_Index,
     .      SS_Counter,
     .      Start_Chan,
     .      Stop_Chan,
     .      Skip_Chan,
     .      Write_Flag,
     .      Instr_Qual,
     .      Attit_Qual,
     .      Offgrts,
     .      j,
     .      Num_Bad_Qual(4),
     .      TGain(1000)
      real*4      Temps(10),
     .            SigTemps(10)
      real*4      Min_DiHed
      real*4      Max_DiHed
                
      record /FEC_SSCAL/ Short_Sci(MAX_DATA)
c
c     *** BEGIN ***
c
      do SS_Counter = Start_Chan, Stop_Chan, Skip_Chan
          Bad_Qflags = 0
	  Num_Bad_Qual(SS_Counter) = 0
          do Eng_Counter = Leftmost, Rightmost
c
c           Acquire Dihedral Temperature for Quality Check
c
               Offgrts = FUT_Temp_List(Array(Eng_Counter).En_Analog,
     .                                 Dummy_Sigs,
     .                                 Config.Fex_GrtRawWt,
     .                                 Config.Fex_GrtTrans,
     .                                 0,Singlifg,Temps,SigTemps)
c      
               if( Temps(5) .gt. Min_DiHed .and. Temps(5) .le. Max_DiHed .and.
     .             Array(Eng_Counter).En_Head.DQ_Summary_Flag(SS_Counter) .Le. Instr_Qual .And.
     .             Array(Eng_Counter).En_Head.Att_Summary_Flag(SS_Counter) .Le. Attit_Qual .And.
     .	           (Array(Eng_Counter).CT_Head.Time(1) .ne. 0 .or.
     .              Array(Eng_Counter).CT_Head.Time(2) .ne. 0) .and.
     .             (Array(Eng_Counter).En_Head.Sci_Time(SS_Counter).Bin_Time(1) .ne. 0 .or.
     .              Array(Eng_Counter).En_Head.Sci_Time(SS_Counter).Bin_Time(2) .ne. 0)) then
                  SS_Index = Eng_Counter - Leftmost - Bad_Qflags + 1
                  Short_Sci( SS_Index ).Time(1) = 
     .			Array(Eng_Counter).En_Head.Sci_Time(SS_Counter).Bin_Time(1)
                  Short_Sci( SS_Index ).Time(2) =
     .			Array(Eng_Counter).En_Head.Sci_Time(SS_Counter).Bin_Time(2)
                  Short_Sci( SS_Index ).Channel_Id = SS_Counter
                  Short_Sci( SS_Index ).Pixel_No = -1
                  Short_Sci( SS_Index ).MTM_Scan_Speed = 
     .                  Array(Eng_Counter).Chan(SS_Counter).XMIT_MTM_Speed
                  Short_Sci( SS_Index ).MTM_Scan_Length = 
     .                  Array(Eng_Counter).Chan(SS_Counter).XMIT_MTM_Len
                  Short_Sci( SS_Index ).Sci_Mode =
     .                  Array(Eng_Counter).Chan(SS_Counter).uP_Sci_Mode
                  Short_Sci( SS_Index ).Adds_Per_Group =
     .                  Array(Eng_Counter).Chan(SS_Counter).uP_Adds_Per_Group
                  Short_Sci( SS_Index ).Transition_Flag = .FALSE.
                  Short_Sci( SS_Index ).Pixel_Definition = 'q'
                  Short_Sci( SS_Index ).Skymap_Index = 0
                  Short_Sci( SS_Index ).Data_Quality =
     .                  Array(Eng_Counter).En_Head.DQ_Summary_Flag(SS_Counter)
                  Short_Sci( SS_Index ).Attitude_Quality =
     .                  Array(Eng_Counter).En_Head.Att_Summary_Flag(SS_Counter)
                  Offgrts = FUT_Temp_List(Array(Eng_Counter).En_Analog,
     .                                    Dummy_Sigs,
     .                                    Config.Fex_GrtRawWt,
     .                                    Config.Fex_GrtTrans,
     .                                    0,Singlifg,Temps,SigTemps)
                  Short_Sci( SS_Index ).Xcal_Temp = Temps(1)
                  Short_Sci( SS_Index ).Skyhorn_Temp = Temps(3)
                  Short_Sci( SS_Index ).Refhorn_Temp = Temps(4)
                  Short_Sci( SS_Index ).Ical_Temp = Temps(2)
                  Short_Sci( SS_Index ).Dihedral_Temp = Temps(5)
                  Short_Sci( SS_Index ).Bolometer_Temp = Temps(6+SS_Counter)
                  Short_Sci( SS_Index ).Bol_Cmd_Bias =
     .                  Array(Eng_Counter).En_Stat.Bol_Cmd_Bias(SS_Counter)
                  j = 1
                  do while (j .le. 16)
                     Short_Sci( SS_Index ).Spares(j) = 0
                     j = j + 1
                  enddo

                  TGain(SS_Index) = Array(Eng_Counter).Chan(SS_Counter).Sci_Gain

              else
                  Bad_Qflags = Bad_Qflags + 1
     	          if( (Array(Eng_Counter).CT_Head.Time(1) .ne. 0 .or.
     .                 Array(Eng_Counter).CT_Head.Time(2) .ne. 0) ) then
		     Num_Bad_Qual(SS_Counter) = Num_Bad_Qual(SS_Counter) + 1
		  end if
              end if
          end do

          Num_Recs = Rightmost - Leftmost - Bad_Qflags + 1

	  if (Num_Recs .gt. 0) then
              Short_Sci(Num_Recs).Transition_Flag = .TRUE.

              do Eng_Counter = 1, Num_Recs - 1
                                    
                  if( (Eng_Counter .ge. 1) .and.
     .                 ((Short_Sci(Eng_Counter).MTM_Scan_Speed .ne.
     .  		Short_Sci(Eng_Counter + 1).MTM_Scan_Speed) .or.
     .                  (Short_Sci(Eng_Counter).MTM_Scan_Length .ne.
     .			Short_Sci(Eng_Counter + 1).MTM_Scan_Length))) then
                      Short_Sci(Eng_Counter).Transition_Flag = .TRUE.
                  end if
  
                  if( (Eng_Counter .ge. 1) .and.
     .                 (TGain(Eng_Counter) .ne. TGain(Eng_Counter + 1)) ) then
                      Short_Sci(Eng_Counter).Transition_Flag = .TRUE.
                  end if
  
              end do

              Num_Group = 0
              ix0 = 0

              do Eng_Counter = 1, Num_Recs
                                    
                  Num_Group = Num_Group + 1

                  if( Short_Sci(Eng_Counter).Transition_Flag .eq. .TRUE.) then
                      Num_Div = (Num_Group - 1)/fac_max_num
                      if( Num_Div .gt. 0 ) then
                        Group_Size = ( (Eng_Counter - ix0) / (Num_Div + 1) ) + 1
                        if(Group_Size .gt. fac_max_num)Group_Size = fac_max_num
                        do j = 1, Num_Div
                           ix = ix0 + j * Group_Size
                           Short_Sci(ix).Transition_Flag = .TRUE.
                        end do
                      endif
                      Num_Group = 0
                      ix0 = Eng_Counter
                  endif

              end do

              do Eng_Counter = 1, Num_Recs - 1
                                    
		  if (Write_Flag .eq. 1) then

                     call ct_write_arcv (, 
     .                                LU_Short_Sci(SS_Counter), 
     .                                Short_Sci(Eng_Counter),
     .                                ct_status)

		     if (ct_status(1) .ne. ctp_normal) then
		        call lib$signal(fec_ctwriteerror,%val(1),%val(ct_status(1)))
		        fec_write_records = %loc(fec_ctwriteerror)
	 	        return
	 	     end if

		  end if

              end do

	      if (Write_Flag .eq. 1) then

                 call ct_write_arcv (, 
     .                                LU_Short_Sci(SS_Counter), 
     .                                Short_Sci(Num_Recs),
     .                                ct_status)

		 if (ct_status(1) .ne. ctp_normal) then
		   call lib$signal(fec_ctwriteerror,%val(1),%val(ct_status(1)))
		   fec_write_records = %loc(fec_ctwriteerror)
	 	   return
	 	 end if

	      end if

          end if

      end do

      FEC_Write_Records = %loc(FEC_NORMAL)

      Return
      End
