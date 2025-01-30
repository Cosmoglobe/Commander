!   FIRAS_Calibrate_Covariances (FCC) Install MMS File
!
!   Purpose: Install facility fcc_calibrate_covariances
!
!   Author: S. Alexander, HSTX, 7/93, SER 11189
!
!   Modifications:
!   

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fcc : csdr$system:fcc.exe,-
      csdr$cld:fcc.cld,-
      csdr$library:csdrmsg.olb(fcc_msg),-
      csdr$help:csdrhelp.hlb(fcc)
 ! Facility fcc has been installed.

csdr$system:fcc.exe : fcc.exe
 copy fcc.exe csdr$system:fcc.exe;0

csdr$cld:fcc.cld : fcc.cld
 copy fcc.cld csdr$cld:fcc.cld;0
