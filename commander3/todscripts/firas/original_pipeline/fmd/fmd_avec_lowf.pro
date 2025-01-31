Pro FMD_AVEC_LOWF,error
;

;
;  FMD_AVEC_LOWF drives the FMD_AVECTOR procedure to create an IDL
;  save set of combined low channel A-Vector.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing LOWF_COVAR.ISS,
;                        and where LOWF_AVECTOR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_AVEC_LOWF,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 22-Apr-97.
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
 PRINT,'FMD_AVEC_LOWF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_AVEC_LOWF,error'
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

; Restore Covariance Matrix
; -------------------------
PRINT,' '
PRINT,'Restoring LOWF_COVAR.ISS'
RESTORE,'csdr$firas_out:lowf_covar.iss'
covar = DOUBLE(covar_lowf)
;

; Compute A-Vector
; ----------------
avec_lowf = DBLARR(43)
;
FOR i=0,42 DO BEGIN
 ;
 asum = 0.
 ;
 FOR j=0,42-i DO asum = asum + $
     ( covar(j,j+i) / SQRT(covar(j,j)) / SQRT(covar(j+i,j+i)) )
 ;
 avec_lowf(i) = asum / (43.-i)
 ;
ENDFOR
;

; 
; Create LOWF_COVAR Save Set
; --------------------------
chanscan = 'LOWF'
delf = (MAX(freq_lowf) - MIN(freq_lowf)) / 42.
delta_freq_lowf = FINDGEN(43) * delf
;
sname = 'csdr$firas_out:lowf_avector.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,max_frac,freq_lowf,$
                    delta_freq_lowf,avec_lowf
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LOWF_AVECTOR.ISS" Created.'
PRINT,' '
;

n_stripes=0 & ndf_covar_lowf=0
error = 0
;

RETURN
END
