!
!   Purpose: Install facility fef_extra_factors.
!
!   Author: S. Brodd, HSTX, 9/16/94, SER 11895
!
!   Modifications: S. Brodd, HSTX, 12/19/94, SER 11979
!                  Enhancements for correct FMS operation.
!

.suffixes
.suffixes .pro

.default
      copy/log $(mms$source) $(mms$target);0

fef : csdr$idl:fef_destriper.pro,-
      csdr$idl:fef_1.pro,-
      csdr$idl:fef_2.pro,-
      csdr$idl:fef_3.pro,-
      csdr$idl:fef_4.pro,-
      csdr$idl:fef_5.pro,-
      csdr$idl:fef_6.pro,-
      csdr$idl:fef_7.pro,-
      csdr$idl:fef_8.pro,-
      csdr$idl:fef_9.pro,-
      csdr$idl:fef_dstripe_sky.pro,-
      csdr$idl:fef_spectra.pro,-
      csdr$idl:fef_spectra_x.pro,-
      csdr$idl:fef_spectra_xx.pro,-
      csdr$idl:fef_write.pro,-
      csdr$idl:fef_write_spc.pro,-
      csdr$idl:fef_write_spc_x.pro
 ! Facility fef has been installed.

csdr$idl:fef_destriper.pro   : fef_destriper.pro
csdr$idl:fef_1.pro           : fef_1.pro
csdr$idl:fef_2.pro           : fef_2.pro
csdr$idl:fef_3.pro           : fef_3.pro
csdr$idl:fef_4.pro           : fef_4.pro
csdr$idl:fef_5.pro           : fef_5.pro
csdr$idl:fef_6.pro           : fef_6.pro
csdr$idl:fef_7.pro           : fef_7.pro
csdr$idl:fef_8.pro           : fef_8.pro
csdr$idl:fef_9.pro           : fef_9.pro
csdr$idl:fef_dstripe_sky.pro : fef_dstripe_sky.pro
csdr$idl:fef_spectra.pro     : fef_spectra.pro
csdr$idl:fef_spectra_x.pro     : fef_spectra_x.pro
csdr$idl:fef_spectra_xx.pro     : fef_spectra_xx.pro
csdr$idl:fef_write.pro       : fef_write.pro
csdr$idl:fef_write_spc.pro   : fef_write_spc.pro
csdr$idl:fef_write_spc_x.pro   : fef_write_spc_x.pro
