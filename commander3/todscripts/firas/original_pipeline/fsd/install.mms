! Installs facility FSD
!
! Author: Reid Wilson
!	  May 20, 1987
!	  STX Incorporated

.suffixes
.suffixes .olb .obj .msg .tlb .cld .hlb .hlp .exe

.cld.tlb
  $(libr) $(librflags) $(mms$target) $(mms$source)

fsd : CSDR$SYSTEM:FSD_GLITCHMAP.EXE, -
      CSDR$SYSTEM:FSD_POSERR.EXE, -
      CSDR$SYSTEM:FSD_ASTROPLOTS.EXE, -
      CSDR$LIBRARY:CSDRMSG.OLB(FSD_MSG), -
      CSDR$HELP:CSDRHELP.HLB(FSD.HLP), -
      CSDR$CLD:FSD.CLD
  ! Facility FSD has been installed.

csdr$system:fsd_glitchmap.exe : fsd_glitchmap.exe
  COPY fsd_glitchmap.exe csdr$system:fsd_glitchmap.exe;0
  ! FSD_GLITCHMAP.EXE has been installed.

csdr$system:fsd_astroplots.exe : fsd_astroplots.exe
  COPY fsd_astroplots.exe csdr$system:fsd_astroplots.exe;0
  ! FSD_ASTROPLOTS.EXE has been installed.

csdr$system:fsd_poserr.exe : fsd_poserr.exe
  COPY fsd_poserr.exe csdr$system:fsd_poserr.exe;0
  ! FSD_POSERR.EXE has been installed.

csdr$cld:fsd.cld : fsd.cld
  COPY fsd.cld csdr$cld:fsd.cld
  ! FSD.CLD has been installed.
