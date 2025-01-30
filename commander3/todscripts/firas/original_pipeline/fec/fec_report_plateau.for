      integer*4 function FEC_Report_Plateau (   LU_Report,
     .                                          Plat_Number,
     .                                          Num_Good,
     .                                          First_Com_Point,
     .                                          Last_Com_Point,
     .                                          First_Point,
     .                                          Last_Point,
     .                                          Num_Total,
     .                                          Num_Bad_Tol,
     .                                          Num_Bad_Qual,
     .                                          Num_Stranded,
     .                                          Tol_type,
     .                                          Tolerance,
     .                                          Mean_A,
     .                                          Mean_B)
c
c      By Reid Wilson, SASC Technologies Inc.,  22-OCT-1986
c
c      Description:
c            This routine writes the message that a plateau is stable to the
c            report file, including information on the parameters.
c
c      Format:
c            Ret_Code = FEC_Report_Plateau (   LU_Report,
c                                              Plat_Number,
c                                              Num_Good,
c                                              First_Com_Point,
c                                              Last_Com_Point,
c                                              First_Point,
c                                              Last_Point,
c                                              Num_Total,
c                                              Num_Bad_Tol,
c                                              Num_Bad_Qual,
c                                              Num_Stranded,
c                                              Tol_type,
c                                              Tolerance,
c                                              Mean_A,
c                                              Mean_B)
c
c      Input Parameters:
c            LU_Report             = logical unit of report file
c                                    integer*2
c
c            Plat_Number           = index of plateau 
c                                    integer*4
c
c            Num_Good              = number of good data points in plateau
c                                    integer*4
c
c            First_Com_Point       = First Commanded Plateau point
c                                    record /FDQ_ENG/ First_Com_Point
c
c            Last_Com_Point        = Last Commanded Plateau point
c                                    record /FDQ_ENG/ Last_Com_Point
c
c            First_Point           = the first data point in the plateau
c                                    record /FDQ_ENG/ First_Point
c
c            Last_Point            = the last data point in the plateau
c                                    record /FDQ_ENG/ Last_Point
c
c            Mean_A                = A side mean temperatures used for plateau
c                                    real*4 Mean(Num_Temps_Used)
c
c            Mean_B                = B side mean temperatures used for plateau
c                                    real*4 Mean(Num_Temps_Used)
c
c            Num_Total             = total number of data points in plateau
c                                    integer*4
c
c            Num_Bad_Tol           = number of bad data points in plateau
c                                    integer*4
c
c            Num_Bad_Qual          = number of bad quality data points in
c				     plateau; integer*4
c
c            Num_Stranded          = number of data points not examined
c                                    integer*4
c
c            Tolerance             = tolerance allowed between actual and
c                                    average temperatures
c                                    real*4 Tolerance(Num_Temps_Used)
c
c      Output Parameters:
c            Ret_Code              = return status of operation
c                                    integer*4
c
c      Modules Called:
c            None
c
c      Include Files:
c            FEC_INC               = contains definitions used in facility
c            FEC_MSG               = contains message definitions
c
c Changes:
c	SPR 3823, Ensure command lines longer than 132 can be written to the
c	report file.  R. Kummerer, May 16, 1989.
c
c       SER 6851, Modifier report file to include the absolute and relative
c                 difference tolerance method
c                 H. Wang, June 13, 1990, STX
c
c     ******* BEGIN *******
c
      implicit none
c
c     Include Files:
c
      include '(FEC_INC)'          !contains normal return value,
                                   !values of Tolerance
      include '(FEC_MSG)'          ! contains FEC message definitions
      include '($JPIDef)'
      include '(FUT_Params)'
c
c     Formal Parameters:
c
      integer*4                    LU_Report
      integer*4                    Num_Good
      DICTIONARY 'FDQ_ENG'         
      record /FDQ_ENG/             First_Com_Point
      record /FDQ_ENG/             Last_Com_Point
      record /FDQ_ENG/             First_Point
      record /FDQ_ENG/             Last_Point
      integer*4                    Plat_Number
      integer*4                    Num_Total
      integer*4                    Num_Bad_Tol
      integer*4                    Num_Bad_Qual(4)
      integer*4                    Num_Stranded
      integer*4                    Tol_type 
      real*4                       Mean_A(Num_Temps_Used)
      real*4                       Mean_B(Num_Temps_Used)
      real*4                       Tolerance(Num_Temps_Used)
      integer*4                    itol(Num_Temps_Used)
      real*4                       rtol(Num_Temps_Used)
      integer*4                    iloop
      integer*4                    i
      integer*4			   Ret_Code
      logical*1			   first_time/.true./  !Need this ???

      do iloop = 1, Num_Temps_Used
          if (tol_type .eq. 1) then
            itol(iloop) = tolerance(iloop) * 1000
          else
            rtol(iloop) = tolerance(iloop) * 100.0
          endif
      end do
      If (Plat_Number .gt. 0) Write (LU_Report, 100) Plat_Number
      If (Plat_Number .eq. 0) Write (LU_Report, 101) Plat_Number
      If (Plat_Number .lt. 0) Write (LU_Report, 102) Plat_Number
      Write (LU_Report, 110) First_Com_Point.CT_Head.GMT(1:11),
     .                       Last_Com_Point.CT_Head.GMT(1:11)
      Write (LU_Report, 115) First_Point.CT_Head.GMT(1:11),
     .                       Last_Point.CT_Head.GMT(1:11)
      If (Plat_Number .lt. 0) Write (LU_Report, 116)
      If (Plat_Number .ge. 0) Then
        Write (LU_Report, 120) Num_Total,
     .                         Num_Good,
     .                         Num_Bad_Tol,
     .                         Num_Stranded
        Write (LU_Report, 130) ((Num_Bad_Qual(i),fac_channel_ids(i)),i=1,4)
        if (tol_type .eq. 1) then
          Write (LU_Report, 140) Mean_A,Mean_B,Itol
        endif
        if (tol_type .eq. 2) then
          Write (LU_Report, 141) Mean_A,Mean_B,Rtol
        endif
      Endif

100   Format (///,' Plateau Number ',I4,/)
101   Format (///,' Plateau Number ',I4,'   *** REJECTED PLATEAU ***',/)
102   Format (///,' Plateau Number ',I4,'   *** HOT SPOT PLATEAU ***',/)
110   Format ('   Commanded Plateau Timerange:  ',A11,' to ',A11,/)
115   Format ('   Plateau Timerange:  ',A11,' to ',A11,/)
116   Format ('   Hot Spot Plateau : No Analysis Performed',/)
120   Format ('   Number of Points:   ',I5,' Total ',I5,' Good ',I5,
     .	      ' Bad Tolerance ',I5,' Stranded')
130   Format ('                       ',4(I5,1X,A2,1X),' Bad Quality',/)
140   Format (1X,32X,'XCAL',5X,'ICAL',5X,'SKYH',5X,'REFH',/,
     .        1X,'  Mean Temperature (A side)',4(1X,F8.3),/,
     .        1X,'  Mean Temperature (B side)',4(1X,F8.3),/,
     .        1X,'  Absolute Tolerance(mK)',3x,4(5X,I4),/)
141   Format (1X,32X,'XCAL',5X,'ICAL',5X,'SKYH',5X,'REFH',/,
     .        1X,'  Mean Temperature (A side)',4(1X,F8.3),/,
     .        1X,'  Mean Temperature (B side)',4(1X,F8.3),/,
     .        1X,'  Relative Tolerance(%)',4x,4(2X,F7.3),/)

      FEC_Report_Plateau = %loc(FEC_NORMAL)
      Return
      End
