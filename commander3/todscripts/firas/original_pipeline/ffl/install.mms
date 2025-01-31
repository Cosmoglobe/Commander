!
!   Purpose: Install facility ffl_fishinput_long
!
!   Author: Gene Eplee, GSC, 9/23/92, SER 10763
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

ffl : csdr$system:ffl.exe,-
      csdr$cld:ffl.cld,-
      csdr$library:csdrmsg.olb(ffl_msg)
 ! Facility ffl has been installed.

csdr$system:ffl.exe : ffl.exe
 copy ffl.exe csdr$system:ffl.exe;0

csdr$cld:ffl.cld : ffl.cld
 copy ffl.cld csdr$cld:ffl.cld;0
