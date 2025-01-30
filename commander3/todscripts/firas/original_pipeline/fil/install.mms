!   FIRAS_Interferogram_Long (FIL) Install MMS File
!
!   Purpose: Install facility fil_interferogram_long.
!
!   Author: S. Brodd, HSTX, 4/95, SER 12244
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fil : csdr$system:fil.exe,-
      csdr$cld:fil.cld,-
      csdr$library:csdrmsg.olb(fil_msg),-
      csdr$help:csdrhelp.hlb(fil)
 ! Facility fil has been installed.

csdr$system:fil.exe : fil.exe
 copy fil.exe csdr$system:fil.exe;0

csdr$cld:fil.cld : fil.cld
 copy fil.cld csdr$cld:fil.cld;0
