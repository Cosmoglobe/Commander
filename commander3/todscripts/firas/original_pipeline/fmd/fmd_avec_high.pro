Pro FMD_AVEC_HIGH,error
;

;
;  FMD_AVEC_HIGH drives the FMD_AVECTOR procedure to create an IDL
;  save set of combined high-channel A-Vector.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_COVAR.ISS,
;                        and where HIGH_AVECTOR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_AVEC_HIGH,error
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
 PRINT,'FMD_AVEC_HIGH : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_AVEC_HIGH,error'
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
PRINT,'Restoring HIGH_COVAR.ISS'
RESTORE,'csdr$firas_out:high_covar.iss'
covar = DOUBLE(covar_high)
;

; Compute A-Vector
; ----------------
avec_high = DBLARR(170)
;
FOR i=0,169 DO BEGIN
 ;
 asum = 0.
 ;
 FOR j=0,169-i DO asum = asum + $
     ( covar(j,j+i) / SQRT(covar(j,j)) / SQRT(covar(j+i,j+i)) )
 ;
 avec_high(i) = asum / (170.-i)
 ;
ENDFOR
;

; 
; Create HIGH_COVAR Save Set
; --------------------------
chanscan = 'HIGH'
freq_high = [freq_high_2,freq_high_3,freq_high_4]
delf = (MAX(freq_high) - MIN(freq_high)) / 169.
delta_freq_high = FINDGEN(170) * delf
;
sname = 'csdr$firas_out:high_avector.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,max_frac,freq_high,$
                    delta_freq_high,avec_high
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_AVECTOR.ISS" Created.'
PRINT,' '
;

n_stripes=0 & ndf_covar_high=0
error = 0
;

RETURN
END
