! Installs facility fpp.
!
! Author: Quoc C. chung
!	  March 6, 1989
!	  STX

.suffixes
.suffixes .olb .obj .msg .tlb .cld .exe .hlb .hlp

fpp : csdr$cld:fpp.cld, -
      csdr$library:csdrmsg.olb(fpp_msg), -
      csdr$help:csdrhelp.hlb(fpp.hlp), -
      csdr$system:fpp.exe
  ! Facility fpp has been installed.

csdr$cld:fpp.cld : fpp.cld
  COPY fpp.cld csdr$cld:fpp.cld;0
  ! fpp.CLD has been installed.

csdr$system:fpp.exe : fpp.exe
  COPY fpp.exe csdr$system:fpp.exe;0
  ! fpp.EXE has been installed.
