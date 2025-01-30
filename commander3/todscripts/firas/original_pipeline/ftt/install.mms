!
!   Purpose: Install facility ftt_tennis_tree.
!
!   Author: S. Brodd, HSTX, 11/22/94, SER 11979
!
!   Modifications: 
!

.suffixes
.suffixes .pro

.default
      copy/log $(mms$source) $(mms$target);0

ftt : csdr$idl:ftt_hifa.pro,-
      csdr$idl:ftt_high.pro,-
      csdr$idl:ftt_hisl.pro,-
      csdr$idl:ftt_hres.pro,-
      csdr$idl:ftt_lofa.pro,-
      csdr$idl:ftt_losl.pro,-
      csdr$idl:ftt_lowf.pro,-
      csdr$idl:ftt_lres.pro
 ! Facility ftt has been installed.

csdr$idl:ftt_hifa.pro : ftt_hifa.pro
csdr$idl:ftt_high.pro : ftt_high.pro
csdr$idl:ftt_hisl.pro : ftt_hisl.pro
csdr$idl:ftt_hres.pro : ftt_hres.pro
csdr$idl:ftt_lofa.pro : ftt_lofa.pro
csdr$idl:ftt_losl.pro : ftt_losl.pro
csdr$idl:ftt_lowf.pro : ftt_lowf.pro
csdr$idl:ftt_lres.pro : ftt_lres.pro
