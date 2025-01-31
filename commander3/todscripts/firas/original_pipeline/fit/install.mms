! Installs facility FIT.
!
! Author: Fred Shuman
!	  November 20, 1987
!	  ST Systems Corporation

.suffixes
.suffixes .olb .obj .msg .tlb .cld .exe .hlb .hlp

fit : csdr$cld:fit.cld, -
      csdr$library:csdrmsg.olb(fit_msg), -
      csdr$help:csdrhelp.hlb(fit=fit.hlp), -
      csdr$system:fit.exe
  ! Facility FIT has been installed.

csdr$cld:fit.cld : fit.cld
  COPY fit.cld csdr$cld:fit.cld;0
  ! FIT.CLD has been installed.

csdr$system:fit.exe : fit.exe
  COPY fit.exe csdr$system:fit.exe;0
  ! FIT.EXE has been installed.
