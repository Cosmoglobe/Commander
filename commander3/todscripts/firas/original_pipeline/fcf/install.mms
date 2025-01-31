!
!   Purpose: Install facility fcf_calibrate_firas.
!
!   Author: S. Alexander, STX, 3/19/92, SER 8292
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fcf : csdr$system:fcf.exe,-
      csdr$cld:fcf.cld,-
      csdr$library:csdrmsg.olb(fcf_msg),-
      csdr$help:csdrhelp.hlb(fcf)
 ! Facility fcf has been installed.

csdr$system:fcf.exe : fcf.exe
 copy fcf.exe csdr$system:fcf.exe;0

csdr$cld:fcf.cld : fcf.cld
 copy fcf.cld csdr$cld:fcf.cld;0
