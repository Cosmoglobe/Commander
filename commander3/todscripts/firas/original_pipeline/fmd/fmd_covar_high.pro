Pro FMD_COVAR_HIGH,error
;

;
;  FMD_COVAR_HIGH drives the FMD_COVAR procedure to create IDL
;  save sets of combined high-channel covariance.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_2_COVAR.ISS,
;                        HIGH_23_COVAR.ISS, HIGH_24_COVAR.ISS,
;                        HIGH_3_COVAR.ISS, and HIGH_34_COVAR.ISS, and
;                        HIGH_4_COVAR.ISS, and where HIGH_COVAR.ISS
;                        will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_COVAR_HIGH,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 14-Apr-97.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;; BEGIN
;
; Initialize Return Error
; -----------------------
error = 1
;

; Procedure Invoked Correctly ?  If not, signal and RETURN
; --------------------------------------------------------
IF N_Params() ne 1 THEN BEGIN 
 PRINT,'FMD_COVAR_HIGH : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_COVAR_HIGH,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 RETURN
ENDIF
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
outtrans = STRUPCASE(outtrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_OUT    == ' + outtrans
PRINT,' '
;

; Restore Covariance Matrices
; ---------------------------
;
PRINT,' '
PRINT,'Restoring HIGH_2_COVAR.ISS'
RESTORE,'csdr$firas_out:high_2_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
cut0 = cvec_cut
max0 = max_frac
;
PRINT,'Restoring HIGH_23_COVAR.ISS'
RESTORE,'csdr$firas_out:high_23_covar.iss'
;
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HIGH_24_COVAR.ISS'
RESTORE,'csdr$firas_out:high_24_covar.iss'
;
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HIGH_3_COVAR.ISS'
RESTORE,'csdr$firas_out:high_3_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HIGH_34_COVAR.ISS'
RESTORE,'csdr$firas_out:high_34_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HIGH_4_COVAR.ISS'
RESTORE,'csdr$firas_out:high_4_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HIGH : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;

; Initialize Output Matrices
; --------------------------
covar_high = FLTARR(170,170)
ndf_covar_high = 0. * covar_high
;

; Fill in Output Matrices
; -----------------------
covar_high(0:54,0:54) = covar_high_2(0:54,0:54)
covar_high(0:54,55:109) = covar_high_23(0:54,0:54)
covar_high(0:54,110:169) = covar_high_24(0:54,0:59)
covar_high(55:109,0:54) = TRANSPOSE(covar_high_23(0:54,0:54))
covar_high(55:109,55:109) = covar_high_3(0:54,0:54)
covar_high(55:109,110:169) = covar_high_34(0:54,0:59)
covar_high(110:169,0:54) = TRANSPOSE(covar_high_24(0:54,0:59))
covar_high(110:169,55:109) = TRANSPOSE(covar_high_34(0:54,0:59))
covar_high(110:169,110:169) = covar_high_4(0:59,0:59)
;
ndf_covar_high(0:54,0:54) = FLOAT(ndf_covar_high2(0))
ndf_covar_high(0:54,55:109) = ndf_covar_hi_23
ndf_covar_high(0:54,110:169) = ndf_covar_hi_24
ndf_covar_high(55:109,0:54) = ndf_covar_hi_23
ndf_covar_high(55:109,55:109) = FLOAT(ndf_covar_high3(0))
ndf_covar_high(55:109,110:169) = ndf_covar_hi_34
ndf_covar_high(110:169,0:54) = ndf_covar_hi_24
ndf_covar_high(110:169,55:109) = ndf_covar_hi_34
ndf_covar_high(110:169,110:169) = FLOAT(ndf_covar_high4(0))
;

; Create HIGH_COVAR Save Set
; --------------------------
chanscan = 'HIGH'
sname = 'csdr$firas_out:high_covar.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,max_frac,freq_high_2,$
                    freq_high_3,freq_high_4,covar_high,ndf_covar_high
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_COVAR.ISS" Created.'
PRINT,' '
;

n_stripes=0
error = 0
;

RETURN
END
