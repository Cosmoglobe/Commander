!
!   Purpose: Install facility ffp_final_product
!
!   Author: S. Brodd, HSTX, 3/21/96, SPR 12316
!
!   Modifications:
!

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg .tlb .exe

ffp : csdr$cld:ffp.cld,-
      csdr$library:csdrmsg.olb(ffp_msg),-
      csdr$help:csdrhelp.hlb(ffp),-
      csdr$system:ffp_spectra_coadd.exe

 ! Facility ffp has been installed.

csdr$cld:ffp.cld : ffp.cld
 copy ffp.cld csdr$cld:ffp.cld;0
 ! ffp.cld has been installed.

csdr$system:ffp_spectra_coadd.exe : ffp_spectra_coadd.exe
 copy ffp_spectra_coadd.exe csdr$system:ffp_spectra_coadd.exe;0
 ! ffp_spectra_coadd.exe has been installed
