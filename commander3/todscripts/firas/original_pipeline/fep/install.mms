! Installs facility FEP.
!
! Author: Rob Kummerer
!	  April 20, 1987
!	  ST Systems Corporation

.suffixes
.suffixes .olb .obj .msg .tlb .cld .exe .hlb .hlp 

fep : csdr$cld:fep.cld, -
      csdr$library:csdrmsg.olb(fep_msg), -
      csdr$help:csdrhelp.hlb(firas_eng_plots=fep.hlp), -
      csdr$system:fep_engplots.exe, -
      csdr$system:fep_dwellplot.exe
  ! Facility FEP has been installed.

csdr$cld:fep.cld : fep.cld
  COPY fep.cld csdr$cld:fep.cld;0
  ! FEP.CLD has been installed.

csdr$system:fep_engplots.exe : fep_engplots.exe
  COPY fep_engplots.exe csdr$system:fep_engplots.exe;0
  ! FEP_ENGPLOTS.EXE has been installed.

csdr$system:fep_dwellplot.exe : fep_dwellplot.exe
  COPY fep_dwellplot.exe csdr$system:fep_dwellplot.exe;0
  ! FEP_DWELLPLOT.EXE has been installed.
