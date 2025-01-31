!
!   Purpose: Install facility ffi_fish_input
!
!   Author: Gene Eplee, GSC, 9/23/92, SER 10763
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

ffi : csdr$system:ffi.exe,-
      csdr$cld:ffi.cld,-
      csdr$library:csdrmsg.olb(ffi_msg),-
      csdr$help:csdrhelp.hlb(ffi)
 ! Facility ffi has been installed.

csdr$system:ffi.exe : ffi.exe
 copy ffi.exe csdr$system:ffi.exe;0

csdr$cld:ffi.cld : ffi.cld
 copy ffi.cld csdr$cld:ffi.cld;0
