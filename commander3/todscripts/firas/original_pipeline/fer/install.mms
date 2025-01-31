!
!   Purpose: Install facility fer_fish_errors.
!
!   Author: S. Brodd, HSTX, 9/20/94, SER 11894
!
!   Modifications:
!

.suffixes
.suffixes .pro .awk

.default
      copy/log $(mms$source) $(mms$target);0

fer : csdr$idl:fer_avg_gain.pro,-
      csdr$idl:fer_avg_off.pro,-
      csdr$idl:fer_errors.pro,-
      csdr$idl:fer_err_str.pro,-
      csdr$idl:fer_fex_cvs.pro,-
      csdr$idl:fer_fex_var.pro,-
      csdr$idl:fer_get_cvs.pro,-
      csdr$idl:fer_get_var.pro,-
      csdr$idl:fer_pep.pro,-
      csdr$idl:fer_process_jcj.pro,-
      csdr$idl:fer_ptp.pro,-
      csdr$idl:fer_pup.pro,-
      csdr$idl:fer_read_jcj.pro,-
      csdr$system:fer_rhss.awk,-      
      csdr$system:fer_rhss0.awk,-     
      csdr$system:fer_rhlf.awk,-     
      csdr$system:fer_rhlf0.awk,-     
      csdr$system:fer_rlss.awk,-     
      csdr$system:fer_rlss0.awk,-     
      csdr$system:fer_rlsf.awk,-      
      csdr$system:fer_rlsf0.awk,-     
      csdr$system:fer_rllf.awk,-     
      csdr$system:fer_rllf0.awk,-     
      csdr$system:fer_lhss.awk,-      
      csdr$system:fer_lhss0.awk,-    
      csdr$system:fer_lhlf.awk,-      
      csdr$system:fer_lhlf0.awk,-     
      csdr$system:fer_llss.awk,-      
      csdr$system:fer_llss0.awk,-     
      csdr$system:fer_llsf.awk,-      
      csdr$system:fer_llsf0.awk,-    
      csdr$system:fer_lllf.awk,-      
      csdr$system:fer_lllf0.awk
 ! Facility fer has been installed.

csdr$idl:fer_avg_gain.pro    : fer_avg_gain.pro
csdr$idl:fer_avg_off.pro     : fer_avg_off.pro
csdr$idl:fer_errors.pro      : fer_errors.pro
csdr$idl:fer_err_str.pro     : fer_err_str.pro
csdr$idl:fer_fex_cvs.pro     : fer_fex_cvs.pro
csdr$idl:fer_fex_var.pro     : fer_fex_var.pro
csdr$idl:fer_get_cvs.pro     : fer_get_cvs.pro
csdr$idl:fer_get_var.pro     : fer_get_var.pro
csdr$idl:fer_pep.pro         : fer_pep.pro
csdr$idl:fer_process_jcj.pro : fer_process_jcj.pro
csdr$idl:fer_ptp.pro         : fer_ptp.pro
csdr$idl:fer_pup.pro         : fer_pup.pro
csdr$idl:fer_read_jcj.pro    : fer_read_jcj.pro
csdr$system:fer_rhss.awk     : fer_rhss.awk
csdr$system:fer_rhss0.awk    : fer_rhss0.awk
csdr$system:fer_rhlf.awk     : fer_rhlf.awk
csdr$system:fer_rhlf0.awk    : fer_rhlf0.awk
csdr$system:fer_rlss.awk     : fer_rlss.awk
csdr$system:fer_rlss0.awk    : fer_rlss0.awk
csdr$system:fer_rlsf.awk     : fer_rlsf.awk
csdr$system:fer_rlsf0.awk    : fer_rlsf0.awk
csdr$system:fer_rllf.awk     : fer_rllf.awk
csdr$system:fer_rllf0.awk    : fer_rllf0.awk
csdr$system:fer_lhss.awk     : fer_lhss.awk
csdr$system:fer_lhss0.awk    : fer_lhss0.awk
csdr$system:fer_lhlf.awk     : fer_lhlf.awk
csdr$system:fer_lhlf0.awk    : fer_lhlf0.awk
csdr$system:fer_llss.awk     : fer_llss.awk
csdr$system:fer_llss0.awk    : fer_llss0.awk
csdr$system:fer_llsf.awk     : fer_llsf.awk
csdr$system:fer_llsf0.awk    : fer_llsf0.awk
csdr$system:fer_lllf.awk     : fer_lllf.awk
csdr$system:fer_lllf0.awk    : fer_lllf0.awk
