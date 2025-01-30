!
!   Purpose: Install facility fel_fish_errors_longspectra.
!
!   Author: S. Brodd, HSTX, 6/6/96, SER 12332
!
!   Modifications:
!
.suffixes
.suffixes .pro

.default
      copy/log $(mms$source) $(mms$target);0

fel : csdr$idl:fel_avg_gain.pro,-
      csdr$idl:fel_avg_off.pro,-
      csdr$idl:fel_errors.pro,-
      csdr$idl:fel_err_str.pro,-
      csdr$idl:fel_fsl_dsk_st.pro,-
      csdr$idl:fel_pep.pro,-
      csdr$idl:fel_process_jcj.pro,-
      csdr$idl:fel_ptp.pro,-
      csdr$idl:fel_pup.pro,-
      csdr$idl:fel_read_jcj.pro
 ! Facility fel has been installed.

csdr$idl:fel_avg_gain.pro    : fel_avg_gain.pro
csdr$idl:fel_avg_off.pro     : fel_avg_off.pro
csdr$idl:fel_errors.pro      : fel_errors.pro
csdr$idl:fel_err_str.pro     : fel_err_str.pro
csdr$idl:fel_fsl_dsk_st.pro  : fel_fsl_dsk_st.pro
csdr$idl:fel_pep.pro         : fel_pep.pro
csdr$idl:fel_process_jcj.pro : fel_process_jcj.pro
csdr$idl:fel_ptp.pro         : fel_ptp.pro
csdr$idl:fel_pup.pro         : fel_pup.pro
csdr$idl:fel_read_jcj.pro    : fel_read_jcj.pro
