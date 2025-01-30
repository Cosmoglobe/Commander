!*******************************************************************!
!***                                                             ***!
!***                      INSTALL.MMS for FDP                    ***!
!***                      Written by Reid Wilson                 ***!
!***                      12-30-86                               ***!
!***                      STX Inc.                               ***!
!***                                                             ***!
!*******************************************************************!

.suffixes
.suffixes .olb .obj .msg .tlb .cld .exe .hlb .hlp

FGA : csdr$cld:FGA.CLD,-
      csdr$help:csdrhelp.hlb(FGA.HLP), -
      csdr$library:csdrmsg.olb(FGA_MSG), -
      csdr$system:FGA.EXE
 ! The FGA facility has been installed

csdr$cld:FGA.CLD : FGA.CLD
 copy FGA.CLD csdr$cld:FGA.CLD;0
 ! FGA.CLD has been installed

csdr$system:FGA.EXE : FGA.EXE
 copy FGA.EXE csdr$system:FGA.EXE;0
 ! FGA.EXE has been installed
