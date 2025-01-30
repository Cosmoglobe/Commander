Pro FMD_AVEC_HRES,error
;

;
;  FMD_AVEC_HRES drives the FMD_AVECTOR procedure to create an IDL
;  save set of combined HRES A-Vector.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HRES_COVAR.ISS,
;                        and where HRES_AVECTOR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_AVEC_HRES,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 09-May-97.
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
 PRINT,'FMD_AVEC_HRES : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_AVEC_HRES,error'
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
PRINT,'Restoring HRES_COVAR.ISS'
RESTORE,'csdr$firas_out:hres_covar.iss'
covar = DOUBLE(covar_hres)
;

; Compute A-Vector
; ----------------
avec_hres = DBLARR(182)
;
FOR i=0,181 DO BEGIN
 ;
 asum = 0.
 ;
 FOR j=0,181-i DO asum = asum + $
     ( covar(j,j+i) / SQRT(covar(j,j)) / SQRT(covar(j+i,j+i)) )
 ;
 avec_hres(i) = asum / (182.-i)
 ;
ENDFOR
;

; 
; Create HRES_COVAR Save Set
; --------------------------
chanscan = 'HRES'
delf = (MAX(freq_hres) - MIN(freq_hres)) / 181.
delta_freq_hres = FINDGEN(182) * delf
;
sname = 'csdr$firas_out:hres_avector.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,max_frac,freq_hres,$
                    delta_freq_hres,avec_hres
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_AVECTOR.ISS" Created.'
PRINT,' '
;

n_stripes=0 & ndf_covar_hres=0
error = 0
;

RETURN
END
