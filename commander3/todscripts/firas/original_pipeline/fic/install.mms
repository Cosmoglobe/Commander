!   FIRAS_Interferogram_Coaddition (FIC) Install MMS File
!
!   Purpose: Install facility fic_interferogram_coadd.
!
!   Author: J. Durachta, STI, 10/86
!
!   Modifications:
!   
!   SER 7985, S. Alexander, 9/1/91: Implement new fic requirements.

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fic : csdr$system:fic.exe,-
      csdr$cld:fic.cld,-
      csdr$library:csdrmsg.olb(fic_msg),-
      csdr$help:csdrhelp.hlb(fic)
 ! Facility fic has been installed.

csdr$system:fic.exe : fic.exe
 copy fic.exe csdr$system:fic.exe;0

csdr$cld:fic.cld : fic.cld
 copy fic.cld csdr$cld:fic.cld;0
