!
!   Purpose: Install facility FSL_Spectra_Long
!
!   Author: Shirley M. Read, Hughes STX Corporation, 8/15/95, SPR 12288.
!
!   Modifications:

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fsl : csdr$system:fsl.exe,-
      csdr$cld:fsl.cld,-
      csdr$library:csdrmsg.olb(fsl_msg),-
      csdr$help:csdrhelp.hlb(fsl)
 ! Facility fsl has been installed.

csdr$system:fsl.exe : fsl.exe
 copy fsl.exe csdr$system:fsl.exe;0

csdr$cld:fsl.cld : fsl.cld
 copy fsl.cld csdr$cld:fsl.cld;0
