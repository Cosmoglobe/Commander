Pro FMD_BAD_LO,error
;

;
;  FMD_BAD_LO creates an IDL save set of bad coadd indices for
;             the four low CHANSCANs.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_BAD_LLS
;                        FMD_BAD_RLS
;                        FMD_BAD_LSF
;                        FMD_BAD_RSF
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_DEFAULT.ISS,
;                        and where FMD_BAD_COADD_LO.ISS will be sent.
;
;    CSDR$FIRAS_OUT   =  Directory containing xxxX_SKY_CHI2.ISS and
;                        and xxxX_CAL_CHI2.ISS .
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_BAD_LO,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 23-Apr-1997.
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
 PRINT,'FMD_BAD_LO : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_BAD_LO,error'
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

; Bad LLS Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD_LLS'
FMD_BAD_LLS,badcoadd_lls,badcal_lls,error
IF (error NE 0) then begin
 print,'FMD_BAD_LO : Error Returned from FMD_BAD_LLS !'
 RETURN
ENDIF
;

; Bad RLS Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD_RLS'
FMD_BAD_RLS,badcoadd_rls,badcal_rls,error
IF (error NE 0) then begin
 print,'FMD_BAD_LO : Error Returned from FMD_BAD_RLS !'
 RETURN
ENDIF
;

; Bad LSF Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD_LSF'
FMD_BAD_LSF,badcoadd_lsf,badcal_lsf,error
IF (error NE 0) then begin
 print,'FMD_BAD_LO : Error Returned from FMD_BAD_LSF !'
 RETURN
ENDIF
;

; Bad RSF Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD_RSF'
FMD_BAD_RSF,badcoadd_rsf,badcal_rsf,error
IF (error NE 0) then begin
 print,'FMD_BAD_LO : Error Returned from FMD_BAD_RSF !'
 RETURN
ENDIF
;

; Make FMD_BAD_COADD_LO Save Set
; ------------------------------
sname='csdr$firas_ref:fmd_bad_coadd_lo.iss'
SAVE,filename=sname,badcoadd_lls,badcoadd_rls,badcoadd_lsf,badcoadd_rsf, $
                    badcal_lls,badcal_rls,badcal_lsf,badcal_rsf
PRINT,' '
PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_BAD_COADD_LO.ISS" Created.'
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
