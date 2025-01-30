!
!   Purpose: Install facility fcs_combine_spectra.
!
!   Author: S. Alexander, STX, 3/19/92, SER 10960
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fcs : csdr$system:fcs.exe,-
      csdr$cld:fcs.cld,-
      csdr$library:csdrmsg.olb(fcs_msg),-
      csdr$help:csdrhelp.hlb(fcs)
 ! Facility fcs has been installed.

csdr$system:fcs.exe : fcs.exe
 copy fcs.exe csdr$system:fcs.exe;0

csdr$cld:fcs.cld : fcs.cld
 copy fcs.cld csdr$cld:fcs.cld;0
