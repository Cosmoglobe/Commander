Pro FMD_DEFLT_QUALS
;

;
;  FMD_DEFLT_QUALS creates IDL save set containing FMD default qualifiers.
;  The default qualifiers are HARD-CODED in this program
;
;  ARGUMENTS  :  None
;
;  PROGRAMS Called  :  None
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory where FMD_QUALS_DEFAULT.ISS will be sent.
;
;
;  HISTORY  :  Written by Ken Jensen, Hughes STX, 23-Apr-1997.
;              Modified by K.Jensen, 27-May-97,  revised BADCOADD indices,
;              required because the coadd sort algorithm in FMD_READ
;              was changed.
;              Modified by K.Jensen, 30-May-97,  revised BADCOADD indices.
;      
;
;


; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,' '
;

del_temp = 2.25E-3
freq_corr = DOUBLE(1.00159)
dihed_cut_min = 1.9
;

; Bad Coadd Defaults
; ------------------
badcoadd_lhs = [1471,4988,5039,20306,22193,29834,33011]  &  badcal_lhs = [-1]
badcoadd_rhs= [1960,10290,13529,20469,20878,21903,21918,29944,30481,32168,32170]
badcal_rhs = [-1]
badcoadd_lhf = [-1]  &  badcal_lhf = [-1]
badcoadd_rhf = [15091,17424]  &  badcal_rhf = [-1]
badcoadd_lls = [1435,1777,21842,26706,27898,32076,34692]  &  badcal_lls = [-1]
badcoadd_rls = [1460,19994,21870,33954]  &  badcal_rls = [-1]
badcoadd_lsf = [-1]  &  badcal_lsf = [-1]
badcoadd_rsf = [-1]  &  badcal_rsf = [-1]
badcoadd_llf = [101,102,4800,7587,13381]  &  badcal_llf = [-1]
badcoadd_rlf = [14480,23751,23788]  &  badcal_rlf = [-1]
;

; Reference File of FMD Default Qualifiers
; ----------------------------------------
sname = 'csdr$firas_ref:fmd_quals_default.iss'
SAVE,filename=sname,del_temp,freq_corr,dihed_cut_min,$
                    badcoadd_rhs,badcoadd_rhf,badcoadd_lhs,badcoadd_lhf,$
                    badcal_rhs,badcal_rhf,badcal_lhs,badcal_lhf,$
                    badcoadd_rls,badcoadd_rsf,badcoadd_lls,badcoadd_lsf,$
                    badcal_rls,badcal_rsf,badcal_lls,badcal_lsf,$
                    badcoadd_llf,badcoadd_rlf,badcal_llf,badcal_rlf
PRINT,' '
PRINT,'IDL Save Set "'+reftrans(0)+'FMD_QUALS_DEFAULT.ISS" Created.'
PRINT,' '
;

RETURN
END
