Pro FMD_BAD_HR,error
;

;
;  FMD_BAD_HR creates an IDL save set of bad coadd indices for
;             the LLLF and RLLF .
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_BAD_LLF
;                        FMD_BAD_RLF
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_DEFAULT.ISS,
;                        and where FMD_BAD_COADD_HR.ISS will be sent.
;
;    CSDR$FIRAS_OUT   =  Directory containing xxxX_SKY_CHI2.ISS and
;                        and xxxX_CAL_CHI2.ISS .
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_BAD_HR,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 09-May-1997.
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
 PRINT,'FMD_BAD_HR : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_BAD_HR,error'
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
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_OUT    == ' + outtrans
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,' '
;

; Bad LLF Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD_LLF'
FMD_BAD_LLF,badcoadd_llf,badcal_llf,error
IF (error NE 0) then begin
 print,'FMD_BAD_HR : Error Returned from FMD_BAD_LLF !'
 RETURN
ENDIF
;

; Bad RLF Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD_RLF'
FMD_BAD_RLF,badcoadd_rlf,badcal_rlf,error
IF (error NE 0) then begin
 print,'FMD_BAD_HR : Error Returned from FMD_BAD_RLF !'
 RETURN
ENDIF
;

; Make FMD_BAD_COADD_HR Save Set
; ------------------------------
sname='csdr$firas_ref:fmd_bad_coadd_hr.iss'
SAVE,filename=sname,badcoadd_llf,badcoadd_rlf,badcal_llf,badcal_rlf
PRINT,' '
PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_BAD_COADD_HR.ISS" Created.'
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
