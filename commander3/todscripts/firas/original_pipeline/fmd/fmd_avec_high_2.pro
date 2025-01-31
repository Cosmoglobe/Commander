Pro FMD_AVEC_HIGH_2,error
;

;
;  FMD_AVEC_HIGH_2 drives the FMD_AVECTOR procedure to create an IDL
;  save set of combined high-channel Band_2 A-Vector.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_COVAR_2.ISS,
;                        and where HIGH_2_AVECTOR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_AVEC_HIGH_2,error
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
 PRINT,'FMD_AVEC_HIGH_2 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_AVEC_HIGH_2,error'
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
PRINT,'Restoring HIGH_2_COVAR.ISS'
RESTORE,'csdr$firas_out:high_2_covar.iss'
covar = DOUBLE(covar_high_2)
;

; Compute A-Vector
; ----------------
avec_high_2 = DBLARR(55)
;
FOR i=0,54 DO BEGIN
 ;
 asum = 0.
 ;
 FOR j=0,54-i DO asum = asum + $
     ( covar(j,j+i) / SQRT(covar(j,j)) / SQRT(covar(j+i,j+i)) )
 ;
 avec_high_2(i) = asum / (55.-i)
 ;
ENDFOR
;

; 
; Create HIGH_COVAR Save Set
; --------------------------
chanscan = 'HIGH_2'
delf = (MAX(freq_high_2) - MIN(freq_high_2)) / 54.
delta_freq_hi_2 = FINDGEN(55) * delf
;
sname = 'csdr$firas_out:high_2_avector.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,max_frac,freq_high_2,$
                    delta_freq_hi_2,avec_high_2
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_2_AVECTOR.ISS" Created.'
PRINT,' '
;

n_stripes=0 & ndf_covar_high2=0
error = 0
;

RETURN
END
