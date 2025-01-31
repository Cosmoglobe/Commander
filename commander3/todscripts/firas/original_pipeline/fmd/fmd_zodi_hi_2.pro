Pro FMD_ZODI_HI_2,chanscan_array,error
;

;
;  FMD_ZODI_HI_2 creates an IDL save set of High Channel Band_2 zodi spectra.
;
;
;  ARGUMENTS (I/O)      :
;
;   CHANSCAN_ARRAY (I)  :  "Y(es)" to process LHS_RHS_LHF_RHF
;
;   ERROR (O)           :  Return Error Status
;
;
;  PROGRAMS Called    :  None
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_2.ISS.
;
;    CSDR$FIRAS_IN    =  Directory containing HIGH_ZODISPEC.ISS and
;                        xxx_2_CALSPEC.ISS, and whereoutput IDL save set
;                        xxx_2_ZODISPEC.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_ZODI_HI_2,chanscan_array,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 26-Mar-1997.
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
IF N_Params() ne 2 THEN BEGIN 
 PRINT,'FMD_ZODI_HI_2 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_ZODI_HI_2,chanscan_array,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 RETURN
ENDIF
;

chanscan_array = STRUPCASE(chanscan_array)
good = WHERE(chanscan_array EQ 'Y',cgood)
IF (cgood LE 0) THEN BEGIN
 PRINT,'FMD_ZODI_HI_2 : No CHANSCANs Enabled !'
 RETURN
ENDIF
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_in',intrans,/full,/issue_error)
intrans = STRUPCASE(intrans)
;
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_IN     == ' + intrans
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,' '
;

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_hi_2.iss'
;

; Restore ZODI Model
; ------------------
RESTORE,'csdr$firas_in:high_zodispec.iss'
;

idd = ['LHS','RHS','LHF','RHF']
;

found = 0
;

FOR ch=0,3 DO BEGIN
 ;
 IF (chanscan_array(ch) EQ 'Y') THEN BEGIN
  ;

  IF (found EQ 0) THEN BEGIN
   ;

   ; Restore IDL Save Set of Undestriped Spectra
   ; -------------------------------------------
   PRINT,' '
   PRINT,'Restoring ' + intrans(0) + idd(ch) + '_CALSPEC.ISS'
   RESTORE,'csdr$firas_in:' + idd(ch) + '_calspec.iss'
   ;

   ; Frequency Correction
   ; --------------------
   f_hi2 = f * freq_corr
   ;

   ; Frequency Cut
   ; -------------
   nf = WHERE((f_hi2 ge freq_band(0))and(f_hi2 le freq_band(1)),cf)
   IF (cf LE 0) THEN BEGIN
    PRINT,''
    PRINT,'FMD_ZODI_HI_2 : Error in Frequency Array !'
    PRINT,''
    RETURN
   ENDIF
   ;
   f_hi2 = f_hi2(nf)
   zspec = TEMPORARY(zodi_spec(nf,*))
   ;
 
  ENDIF   ; IF(found EQ 0)
  ;

  ; CHANSCAN Band_2 Zodi Spectra
  ; ----------------------------
  zodi_spec = zspec(*,WHERE(sky_idx EQ ch))
  ;

  ; Make xxx_2_ZODISPEC Save Set
  ; ----------------------------
  chan_label = idd(ch) + '_2'
  sname='csdr$firas_in:' + chan_label + '_zodispec.iss'
  SAVE,filename=sname,f_hi2,zodi_spec
  PRINT,' '
  PRINT,'IDL Save Set "' + intrans(0) + idd(ch)+'_2_ZODISPEC.ISS" Created.'
  PRINT,' '
  ;
  
  found = 1
  ;

 ENDIF   ; IF(chanscan_array...
 ;

ENDFOR
;

; Re-Define Unused Restored Parameters
; ------------------------------------
del_temp=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0 & cmn_tm=0
step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0
ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & zodi_array=0
good_sky_dihd=0 & tmin=0 & tmax=0 & good_cal_dihd=0 & cal_tm=0 & cal_sp=0
hot_cal=0 & gain_convg=0 & gain_iter=0 & cvec_cut=0 & latcut=0 & loncut=0
n_stripes=0 & max_frac=0 & alpha=0 & chanscan=0 & delbeta=0 & zodi_freq=0
wdelbeta=0 & zeta=0 & zodi_model=0 & zodimodpars=0 & px=0 & tm=0 & sp=0
badcoadd_rhs=0 & badcoadd_lhs=0 & badcoadd_rhf=0 & badcoadd_lhf=0
badcal_rhs=0 & badcal_lhs=0 & badcal_rhf=0 & badcal_lhf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
