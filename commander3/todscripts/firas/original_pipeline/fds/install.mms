!
!   Purpose: Install facility fds_destriper.
!
!   Author: S. Alexander, STX, 3/19/92, SER 10797
!
!   Modifications: S. Alexander, HSTX, 6/9/93, SER 11039. Add routine
!                  fds_cvector.pro to produce input for fip.
!                  S. Alexander, HSTX, 2/25/94, SER 11402. Change ASCII
!                  C vector to a binary reference dataset.
!                  S. Alexander, HSTX, 4/94, SER 11702. New gain and
!                  offset destriper.
!                  S. Brodd, HSTX, 8/11/94, SER 11871. New procedures
!                  FDS_ERRORS.PRO, FDS_ERRMAT.PRO, FDS_CHOLESKEY.PRO, 
!                  FDS_FIP_DSE.PRO, FDS_FIP_DSQ.PRO

.suffixes
.suffixes .pro

.default
      copy/log $(mms$source) $(mms$target);0

fds : csdr$idl:fds_destriper.pro,-
      csdr$idl:fds_dstripe_sky.pro,-
      csdr$idl:fds_process_sky.pro,-
      csdr$idl:fds_process_var.pro,-
      csdr$idl:fds_read_data.pro,-
      csdr$idl:fds_variance.pro,-
      csdr$idl:fds_cvector.pro,-
      csdr$idl:fds_errors.pro,-
      csdr$idl:fds_errmat.pro,-
      csdr$idl:fds_choleskey.pro,-
      csdr$idl:fds_fex_cvs.pro,-
      csdr$idl:fds_fip_dse.pro,-
      csdr$idl:fds_fip_dsq.pro
 ! Facility fds has been installed.

csdr$idl:fds_destriper.pro     : fds_destriper.pro
csdr$idl:fds_dstripe_sky.pro   : fds_dstripe_sky.pro
csdr$idl:fds_process_sky.pro   : fds_process_sky.pro
csdr$idl:fds_process_var.pro   : fds_process_var.pro
csdr$idl:fds_read_data.pro     : fds_read_data.pro
csdr$idl:fds_variance.pro      : fds_variance.pro
csdr$idl:fds_cvector.pro       : fds_cvector.pro
csdr$idl:fds_errors.pro        : fds_errors.pro
csdr$idl:fds_errmat.pro	       : fds_errmat.pro
csdr$idl:fds_choleskey.pro     : fds_choleskey.pro
csdr$idl:fds_fex_cvs.pro       : fds_fex_cvs.pro
csdr$idl:fds_fip_dse.pro       : fds_fip_dse.pro
csdr$idl:fds_fip_dsq.pro       : fds_fip_dsq.pro
