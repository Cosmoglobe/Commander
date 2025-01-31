Pro FMD_BAD2_HI,error
;

;
;  FMD_BAD2_HI creates an IDL save set of bad coadd indices for
;             the four high CHANSCANs.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_BAD2_LHS
;                        FMD_BAD2_RHS
;                        FMD_BAD2_LHF
;                        FMD_BAD2_RHF
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_BAD_COADD_HI.ISS,
;                        and where new FMD_BAD_COADD_HI.ISS will be sent.
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_n_SKY_CHI2.ISS and
;                        and HIGH_n_CAL_CHI2.ISS .
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_BAD2_HI,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 08-Apr-1997.
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
 PRINT,'FMD_BAD2_HI : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_BAD2_HI,error'
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

; Bad LHS Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD2_LHS'
FMD_BAD2_LHS,badcoadd_lhs,badcal_lhs,error
IF (error NE 0) then begin
 print,'FMD_BAD2_HI : Error Returned from FMD_BAD2_LHS !'
 RETURN
ENDIF
;

; Bad RHS Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD2_RHS'
FMD_BAD2_RHS,badcoadd_rhs,badcal_rhs,error
IF (error NE 0) then begin
 print,'FMD_BAD2_HI : Error Returned from FMD_BAD2_RHS !'
 RETURN
ENDIF
;

; Bad LHF Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD2_LHF'
FMD_BAD2_LHF,badcoadd_lhf,badcal_lhf,error
IF (error NE 0) then begin
 print,'FMD_BAD2_HI : Error Returned from FMD_BAD2_LHF !'
 RETURN
ENDIF
;

; Bad RHF Coadds
; --------------
PRINT,' '
PRINT,'Calling FMD_BAD2_RHF'
FMD_BAD2_RHF,badcoadd_rhf,badcal_rhf,error
IF (error NE 0) then begin
 print,'FMD_BAD2_HI : Error Returned from FMD_BAD2_RHF !'
 RETURN
ENDIF
;

; Make FMD_BAD_COADD_HI Save Set
; ------------------------------
sname='csdr$firas_ref:fmd_bad_coadd_hi.iss'
SAVE,filename=sname,badcoadd_lhs,badcoadd_rhs,badcoadd_lhf,badcoadd_rhf, $
                    badcal_lhs,badcal_rhs,badcal_lhf,badcal_rhf
PRINT,' '
PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_BAD_COADD_HI.ISS" Created.'
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
