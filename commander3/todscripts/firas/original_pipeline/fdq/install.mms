!*******************************************************************!
!***                                                             ***!
!***                      INSTALL.MMS for FDQ                    ***!
!***                      Written by Reid Wilson                 ***!
!***                      12-31-86                               ***!
!***                      STX Inc.                               ***!
!***                                                             ***!
!***			  Modified by Shirley M. Read            ***!
!***			  12-09-87				 ***!
!***			  STX					 ***!
!***                                                             ***!
!***                      Steven Alexander, April 27, 1992:      ***!
!***                             remove H_CONV to new facility   ***!
!***                             FHC. SPR 9642.                  ***!
!***                                                             ***!
!*******************************************************************!

.suffixes
.suffixes .olb .obj .msg .msg~ .tlb .cld .cld~  .hlb .hlp .hlp~ .exe

.cld.tlb
 $(libr) $(librflags) $(mms$target) $(mms$source)

fdq : -
      csdr$system:fdq.exe, -
      csdr$cld:fdq.cld,-
      csdr$library:csdrmsg.olb(fdq_msg), -
      csdr$help:csdrhelp.hlb(fdq.hlp)
 ! The fdq facility has been installed

csdr$system:fdq.exe : fdq.exe
 copy fdq.exe csdr$system:fdq.exe
 ! fdq.exe has been installed

csdr$cld:fdq.cld : fdq.cld
 copy $(mms$source) $(mms$target);0
 ! fdq.cld has been installed
