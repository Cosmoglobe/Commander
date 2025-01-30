Pro FMD_CVECTORS_HI
;

;
;  FMD_CVECTORS_HI Makes the CVECTORS_HI Reference Data Set
;
;
;  Required Logicals :
;
;     CSDR$FIRAS_OUT  : Directory containing the xxx_n_CVECTOR.ISS .
;
;     CSDR$FIRAS_REF  : Directory containing CVECTORS_HI_0.ISS, and
;                       where CVECTORS_HI.ISS will be sent.
;
;
;  Written by : Ken Jensen, Hughes STX, 26-Mar-1997
;
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
outtrans = STRUPCASE(outtrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,'CSDR$FIRAS_OUT    == ' + outtrans
PRINT,' '
;

; Restore the Reference C-Vectors
; -------------------------------
RESTORE,'csdr$firas_ref:cvectors_hi_0.iss'
;

; Restore the Destriped LHS C-Vectors
; -----------------------------------
RESTORE,'csdr$firas_out:lhs_2_cvector.iss'
RESTORE,'csdr$firas_out:lhs_3_cvector.iss'
RESTORE,'csdr$firas_out:lhs_4_cvector.iss'
;

; Restore the Destriped RHS C-Vectors
; -----------------------------------
RESTORE,'csdr$firas_out:rhs_2_cvector.iss'
RESTORE,'csdr$firas_out:rhs_3_cvector.iss'
RESTORE,'csdr$firas_out:rhs_4_cvector.iss'
;

; Restore the Destriped LHF C-Vectors
; -----------------------------------
RESTORE,'csdr$firas_out:lhf_2_cvector.iss'
RESTORE,'csdr$firas_out:lhf_3_cvector.iss'
RESTORE,'csdr$firas_out:lhf_4_cvector.iss'
;

; Restore the Destriped RHF C-Vectors
; -----------------------------------
RESTORE,'csdr$firas_out:rhf_2_cvector.iss'
RESTORE,'csdr$firas_out:rhf_3_cvector.iss'
RESTORE,'csdr$firas_out:rhf_4_cvector.iss'
;

; Concatenate the C-Vectors
; -------------------------
cvec_lhs = [cvec_lhs(0:39),cvec_lhs_2,cvec_lhs_3,cvec_lhs_4]
cvec_rhs = [cvec_rhs(0:39),cvec_rhs_2,cvec_rhs_3,cvec_rhs_4]
cvec_lhf = [cvec_lhf(0:39),cvec_lhf_2,cvec_lhf_3,cvec_lhf_4]
cvec_rhf = [cvec_rhf(0:39),cvec_rhf_2,cvec_rhf_3,cvec_rhf_4]
;

; Make the CVECTORS_HI Save Set
; -----------------------------
sname = 'csdr$firas_ref:cvectors_hi.iss'
SAVE,filename=sname,f_hi,cvec_lhs,cvec_rhs,cvec_lhf,cvec_rhf
;
PRINT,' '
PRINT,'IDL Save Set "' + reftrans(0) + 'CVECTORS_HI.ISS" Created.'
PRINT,' '
;

chanscan=0 & stripes_descrip=0 & freq_dependence=0 & cvec_cut=0
freq_lhs_2=0 & freq_rhs_2=0 & freq_rhf_2=0 & freq_lhf_2=0
ndf_cvec_lhs_2=0 & ndf_cvec_rhs_2=0 & ndf_cvec_rhf_2=0 & ndf_cvec_lhf_2=0
freq_lhs_3=0 & freq_rhs_3=0 & freq_rhf_3=0 & freq_lhf_3=0
ndf_cvec_lhs_3=0 & ndf_cvec_rhs_3=0 & ndf_cvec_rhf_3=0 & ndf_cvec_lhf_3=0
freq_lhs_4=0 & freq_rhs_4=0 & freq_rhf_4=0 & freq_lhf_4=0
ndf_cvec_lhs_4=0 & ndf_cvec_rhs_4=0 & ndf_cvec_rhf_4=0 & ndf_cvec_lhf_4=0
;

RETURN
END
