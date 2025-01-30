Pro FMD_AVEC_HIGH_4,error
;

;
;  FMD_AVEC_HIGH_4 drives the FMD_AVECTOR procedure to create an IDL
;  save set of combined high-channel Band_4 A-Vector.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_COVAR_4.ISS,
;                        and where HIGH_4_AVECTOR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_AVEC_HIGH_4,error
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
 PRINT,'FMD_AVEC_HIGH_4 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_AVEC_HIGH_4,error'
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
PRINT,'Restoring HIGH_4_COVAR.ISS'
RESTORE,'csdr$firas_out:high_4_covar.iss'
covar = DOUBLE(covar_high_4)
;

; Compute A-Vector
; ----------------
avec_high_4 = DBLARR(60)
;
FOR i=0,59 DO BEGIN
 ;
 asum = 0.
 ;
 FOR j=0,59-i DO asum = asum + $
     ( covar(j,j+i) / SQRT(covar(j,j)) / SQRT(covar(j+i,j+i)) )
 ;
 avec_high_4(i) = asum / (60.-i)
 ;
ENDFOR
;

; 
; Create HIGH_COVAR Save Set
; --------------------------
chanscan = 'HIGH_4'
delf = (MAX(freq_high_4) - MIN(freq_high_4)) / 59.
delta_freq_hi_4 = FINDGEN(60) * delf
;
sname = 'csdr$firas_out:high_4_avector.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,max_frac,freq_high_4,$
                    delta_freq_hi_4,avec_high_4
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_4_AVECTOR.ISS" Created.'
PRINT,' '
;

n_stripes=0 & ndf_covar_high4=0
error = 0
;

RETURN
END
