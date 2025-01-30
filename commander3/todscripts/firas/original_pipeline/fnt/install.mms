! Installs facility FNT.
!
! Author: Rob Kummerer
!	  December 19, 1986
!	  STX

.suffixes
.suffixes .olb .obj .msg .tlb .cld .exe .hlb .hlp

fnt : csdr$cld:fnt.cld, -
      csdr$library:csdrmsg.olb(fnt_msg), -
      csdr$help:csdrhelp.hlb(fnt.hlp), -
      csdr$system:fnt.exe
  ! Facility FNT has been installed.

csdr$cld:fnt.cld : fnt.cld
  COPY fnt.cld csdr$cld:fnt.cld;0
  ! FNT.CLD has been installed.

csdr$system:fnt.exe : fnt.exe
  COPY fnt.exe csdr$system:fnt.exe;0
  ! FNT.EXE has been installed.
