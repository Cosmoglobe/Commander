!
!   Purpose: Install facility fip_initial_product
!
!   Author: Gene Eplee, GSC, 7/2/93, SER 11058
!
!   Modifications: Gene Eplee, GSC, 11/1/93, SER 11414
!                  Added new executable FIP_LINES.
!                  Larry Rosen, HSTX, 12/6/93, SER 11259
!                  Added new executable FIP_COVAR.
!		Gene Eplee, GSC, 03/02/94
!		Added new executable FIP_EJGN.
!		Gene Eplee, GSC, 05/09/94
!		Added new executable FIP_ERR.  SER 11704
!               Replace FIP_EJGN executable with FIP_EJV executable.
!                  Gene Eplee, GSC, 07/08/94  SER 11618
!               Added new executable FIP_FEF.
!                  Gene Eplee, GSC, 09/08/94
!               Added new executable FIP_DUST.
!                  Gene Eplee, GSC, 10/05/94, SER 11936
!               Larry Rosen, HSTX, 10/94, SER 11704
!                  Added new rdls and FIP_FCF for fcf file conversion.
!               Added new executable FIP_CVS.
!                  Gene Eplee, GSC, 11/17/94, SER 11936
!               FIP_FCF has been replaced by FIP_SPECTRA_COADD (FIPA verb)
!               that processes FCF_SKY, FCF_DSK, FCF_CAL, FCF_DCL, FIC_SKY,
!               and FIC_CAL.
!                  Larry P. Rosen, HSTX, 12/15/94, SER 11936

.suffixes
.suffixes .cld .olb .obj .hlb .hlp .msg .tlb .exe

fip : csdr$cld:fip.cld,-
      csdr$library:csdrmsg.olb(fip_msg),-
      csdr$help:csdrhelp.hlb(fip),-
      csdr$system:fip_sky.exe,-
      csdr$system:fip_model.exe,-
      csdr$system:fip_lines.exe,-
      csdr$system:fip_covar.exe,-
      csdr$system:fip_ejv.exe,-
      csdr$system:fip_err.exe,-
      csdr$system:fip_fef.exe,-
      csdr$system:fip_dust.exe,-
      csdr$system:fip_cvs.exe,-
      csdr$system:fip_spectra_coadd.exe

 ! Facility fip has been installed.

csdr$cld:fip.cld : fip.cld
 copy fip.cld csdr$cld:fip.cld;0
 ! fip.cld has been installed.

csdr$system:fip_sky.exe : fip_sky.exe
 copy fip_sky.exe csdr$system:fip_sky.exe;0
 ! fip_sky.exe has been installed

csdr$system:fip_model.exe : fip_model.exe
 copy fip_model.exe csdr$system:fip_model.exe;0
 ! fip_model.exe has been installed

csdr$system:fip_lines.exe : fip_lines.exe
 copy fip_lines.exe csdr$system:fip_lines.exe;0
 ! fip_lines.exe has been installed

csdr$system:fip_covar.exe : fip_covar.exe
 copy fip_covar.exe csdr$system:fip_covar.exe;0
 ! fip_covar.exe has been installed

csdr$system:fip_ejv.exe : fip_ejv.exe
 copy fip_ejv.exe csdr$system:fip_ejv.exe;0
 ! fip_ejv.exe has been installed

csdr$system:fip_err.exe : fip_err.exe
 copy fip_err.exe csdr$system:fip_err.exe;0
 ! fip_err.exe has been installed

csdr$system:fip_fef.exe : fip_fef.exe
 copy fip_fef.exe csdr$system:fip_fef.exe;0
 ! fip_fef.exe has been installed

csdr$system:fip_dust.exe : fip_dust.exe
 copy fip_dust.exe csdr$system:fip_dust.exe;0
 ! fip_dust.exe has been installed

csdr$system:fip_cvs.exe : fip_cvs.exe
 copy fip_cvs.exe csdr$system:fip_cvs.exe;0
 ! fip_cvs.exe has been installed

csdr$system:fip_spectra_coadd.exe : fip_spectra_coadd.exe
 copy fip_spectra_coadd.exe csdr$system:fip_spectra_coadd.exe;0
 ! fip_spectra_coadd.exe has been installed
