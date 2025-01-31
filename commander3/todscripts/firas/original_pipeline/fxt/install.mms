! Installs facility FXT.
!
! Author: Shirley M. Read
!	  January 10, 1989
!	  STX Incorporated

.suffixes
.suffixes .olb .obj .msg .tlb .cld .exe .hlb .hlp

fxt : csdr$cld:fxt.cld, -
      csdr$library:csdrmsg.olb(fxt_msg), -
      csdr$help:csdrhelp.hlb(fxt=fxt.hlp), -
      csdr$system:fxt.exe
  ! Facility FXT has been installed.

csdr$cld:fxt.cld : fxt.cld
  COPY fxt.cld csdr$cld:fxt.cld;0
  ! FXT.CLD has been installed.

csdr$system:fxt.exe : fxt.exe
  COPY fxt.exe csdr$system:fxt.exe;0
  ! FXT.EXE has been installed.
