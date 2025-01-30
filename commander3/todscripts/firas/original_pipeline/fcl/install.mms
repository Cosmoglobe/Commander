!   FIRAS_CalibrateCovariances_Long (FCL) Install MMS File
!
!   Purpose: Install facility fcl_calibratecovariances_long
!
!   Author: S. Brodd, HSTX, 12/95, SPR 12291
!
!   Modifications:
!   

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg

fcl : csdr$system:fcl.exe,-
      csdr$cld:fcl.cld,-
      csdr$library:csdrmsg.olb(fcl_msg),-
      csdr$help:csdrhelp.hlb(fcl)
 ! Facility fcl has been installed.

csdr$system:fcl.exe : fcl.exe
 copy fcl.exe csdr$system:fcl.exe;0

csdr$cld:fcl.cld : fcl.cld
 copy fcl.cld csdr$cld:fcl.cld;0
