!***************************************************************
! Installs facility FSS
!
! Author: D. Bouler, STX, CDAC, Jan, 1991
!         Modified from FES install file by Reid Wilson
!****************************************************************
! Changes:
!
!****************************************************************

.suffixes
.suffixes .olb .obj .msg .msg~ .tlb .cld .cld~ .hlb .hlp .hlp~ .exe

fss : csdr$cld:fss.cld, -
      csdr$library:csdrmsg.olb(fss_msg), -
      csdr$help:csdrhelp.hlb(fss.hlp), -
      csdr$system:fss.exe
  ! Facility fss has been installed.

csdr$cld:fss.cld : fss.cld
  COPY fss.cld csdr$cld:fss.cld;0
  ! fss.CLD has been installed.

csdr$system:fss.exe : fss.exe
  COPY fss.exe csdr$system:fss.exe;0
  ! fss.EXE has been installed.
