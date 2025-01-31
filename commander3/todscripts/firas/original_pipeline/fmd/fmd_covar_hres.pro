Pro FMD_COVAR_HRES,error
;

;
;  FMD_COVAR_HRES patches segments of the HRES C-Matrix to produce
;  the HRES covariance.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HRES_1_COVAR.ISS,
;                        HRES_12_COVAR.ISS, HRES_13_COVAR.ISS,
;                        HRES_2_COVAR.ISS, and HRES_23_COVAR.ISS, and
;                        HRES_3_COVAR.ISS, and where HRES_COVAR.ISS
;                        will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_COVAR_HRES,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 13-May-97.
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
 PRINT,'FMD_COVAR_HRES : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_COVAR_HRES,error'
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
PRINT,'Restoring HRES_1_COVAR.ISS'
RESTORE,'csdr$firas_out:hres_1_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
cut0 = cvec_cut
max0 = max_frac
;
PRINT,'Restoring HRES_12_COVAR.ISS'
RESTORE,'csdr$firas_out:hres_12_covar.iss'
;
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HRES_13_COVAR.ISS'
RESTORE,'csdr$firas_out:hres_13_covar.iss'
;
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HRES_2_COVAR.ISS'
RESTORE,'csdr$firas_out:hres_2_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HRES_23_COVAR.ISS'
RESTORE,'csdr$firas_out:hres_23_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;
PRINT,'Restoring HRES_3_COVAR.ISS'
RESTORE,'csdr$firas_out:hres_3_covar.iss'
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Frequency-Dependent Input C-Matrix !'
 RETURN
ENDIF
IF (MAX(ABS(cut0-cvec_cut)) NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched CVEC_CUT !'
 RETURN
ENDIF
IF (max0 NE max_frac) THEN BEGIN
 PRINT,'FMD_COVAR_HRES : Mismatched MAX_FRAC !'
 RETURN
ENDIF
;

; Initialize Output Matrices
; --------------------------
covar_hres = FLTARR(182,182)
ndf_covar_hres = 0. * covar_hres
;

; Fill in Output Matrices
; -----------------------
covar_hres(0:59,0:59) = covar_hres_1(0:59,0:59)
covar_hres(0:59,60:119) = covar_hres_12(0:59,0:59)
covar_hres(0:59,120:181) = covar_hres_13(0:59,0:61)
covar_hres(60:119,0:59) = TRANSPOSE(covar_hres_12(0:59,0:59))
covar_hres(60:119,60:119) = covar_hres_2(0:59,0:59)
covar_hres(60:119,120:181) = covar_hres_23(0:59,0:61)
covar_hres(120:181,0:59) = TRANSPOSE(covar_hres_13(0:59,0:61))
covar_hres(120:181,60:119) = TRANSPOSE(covar_hres_23(0:59,0:61))
covar_hres(120:181,120:181) = covar_hres_3(0:61,0:61)
;
ndf_covar_hres(0:59,0:59) = FLOAT(ndf_covar_hres1)
ndf_covar_hres(0:59,60:119) = ndf_covar_hr_12
ndf_covar_hres(0:59,120:181) = ndf_covar_hr_13
ndf_covar_hres(60:119,0:59) = ndf_covar_hr_12
ndf_covar_hres(60:119,50:119) = FLOAT(ndf_covar_hres2)
ndf_covar_hres(60:119,120:181) = ndf_covar_hr_23
ndf_covar_hres(120:181,0:59) = ndf_covar_hr_13
ndf_covar_hres(120:181,60:119) = ndf_covar_hr_23
ndf_covar_hres(120:181,120:181) = FLOAT(ndf_covar_hres3)
;

; Create HRES_COVAR Save Set
; --------------------------
chanscan = 'HRES'
freq_hres = [freq_hr_1,freq_hr_2,freq_hr_3]
;
sname = 'csdr$firas_out:hres_covar.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,max_frac,freq_hres,$
                    covar_hres,ndf_covar_hres
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_COVAR.ISS" Created.'
PRINT,' '
;

n_stripes=0
error = 0
;

RETURN
END
