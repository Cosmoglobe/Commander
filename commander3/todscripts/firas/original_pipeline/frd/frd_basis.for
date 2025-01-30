      PROGRAM  FRD_BASIS
C*******************************************************************************
C     Purpose: The following program generates Legendre polynomials for  
C     512 index points. It then places the resulting values in a 512x5  
C     array that is written out to as a single record to the reference
C     file, FEX_BASIS.DAT
C      
C     Author: John Sims, STX, 15 July 1991, SER 7985
C
C*******************************************************************************
      Implicit None

      Dictionary          'FEX_BASIS'
      Record /FEX_BASIS/  FEX_Rec
      Integer*4           Lun
      Character*14        Out_file / 'FEX_BASIS.DAT;' /
      Integer*4           Istat
      External            FRD_RMSOpen
      integer*4           i
      real*8              A

      do i=1,512
         A = (i/256.D0) - 1.D0
         FEX_Rec.LEG_POLY(i,1) = 1.D0
         FEX_Rec.LEG_POLY(i,2) = 1.D0*A
         FEX_Rec.LEG_POLY(i,3) = -0.5D0 + 1.5D0*(A**2)
         FEX_Rec.LEG_POLY(i,4) = -1.5D0*A + 2.5D0*(A**3)
         FEX_Rec.LEG_POLY(i,5) = 0.375D0 - 3.75D0*(A**2) + 4.375D0*(A**4)  
      end do        
 
      Call Lib$Get_Lun (Lun)     
      Open (UNIT=Lun, FILE=Out_File, STATUS='NEW', FORM='UNFORMATTED',
     .    ACCESS='SEQUENTIAL', RECORDSIZE=5120, ORGANIZATION='SEQUENTIAL',
     .    RECORDTYPE='FIXED', IOSTAT=Istat, SHARED)  
      If (Istat) Then
         Call Lib$Signal (FRD_RMSOpen, %Val(2), Out_File, %Val(Istat))
      Else
         Write (Lun) FEX_Rec
         Close (Lun)
      Endif
      End
