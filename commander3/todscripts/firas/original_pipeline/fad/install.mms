!
!   Purpose: Install facility fad_apply_destriper
!
!   Author: S. Alexander, STX, 5/5/93, SER 10798
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fad : csdr$system:fad.exe,-
      csdr$cld:fad.cld,-
      csdr$help:csdrhelp.hlb(fad.hlp),-
      csdr$library:csdrmsg.olb(fad_msg)
 ! Facility fad has been installed.

csdr$system:fad.exe : fad.exe
 copy fad.exe csdr$system:fad.exe;0

csdr$cld:fad.cld : fad.cld
 copy fad.cld csdr$cld:fad.cld;0
