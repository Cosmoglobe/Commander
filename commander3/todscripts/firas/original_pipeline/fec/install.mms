! Installs facility FEC_EXTRACT_CALIBRATION      
!
! Author: Reid Wilson
!	  October 29, 1986
!	  STI Incorporated

.suffixes
.suffixes .olb .obj .msg .tlb .cld .hlb .hlp .exe

.cld.tlb
  $(libr) $(librflags) $(mms$target) $(mms$source)

fec : csdr$cld:fec.cld, -
      csdr$library:csdrmsg.olb(fec_msg), -
      csdr$help:csdrhelp.hlb(fec.hlp), -
      csdr$system:fec.exe
  ! Facility FEC has been installed.

csdr$cld:fec.cld : fec.cld
  COPY fec.cld csdr$cld:fec.cld;0
  ! FEC.CLD has been installed.

csdr$system:fec.exe : fec.exe
  COPY fec.exe csdr$system:fec.exe;0
  ! FEC.EXE has been installed.
