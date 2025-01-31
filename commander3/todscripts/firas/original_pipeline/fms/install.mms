!
!   Purpose: Install facility fms_merge_skymaps.
!
!   Author: S. Alexander, HSTX, 4/20/94, SER 11418
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fms : csdr$system:fms.exe,-
      csdr$cld:fms.cld,-
      csdr$library:csdrmsg.olb(fms_msg),-
      csdr$help:csdrhelp.hlb(fms)
 ! Facility fms has been installed.

csdr$system:fms.exe : fms.exe
 copy fms.exe csdr$system:fms.exe;0

csdr$cld:fms.cld : fms.cld
 copy fms.cld csdr$cld:fms.cld;0
