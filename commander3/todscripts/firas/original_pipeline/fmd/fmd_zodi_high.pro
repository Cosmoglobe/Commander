Pro FMD_Zodi_High
;

;
;  FMD_Zodi_High Makes the HIGH_ZODISPEC Data Set
;  for later processing by FMD .
;
;
;  Required Logicals :
;
;     CSDR$FIRAS_IN  : Directory containing LHS.ISS, RHS.ISS, LHF.ISS,
;                      RHF.ISS, and FMD_CZM_HIGH.ISS, and where
;                      HIGH_ZODISPEC.ISS will be sent.
;
;
;  Written by : Ken Jensen, Hughes STX, 17-Jun-1997
;
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_in',intrans,/full,/issue_error)
intrans = STRUPCASE(intrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_IN     == ' + intrans
PRINT,' '
;

RESTORE,'csdr$firas_in:fmd_czm_high.iss'
zodi_freq = freq
zodi_spec = TEMPORARY(zodimodel)
zodi_model = STRUPCASE(zodimodid)
;

RESTORE,'csdr$firas_in:lhs.iss'
sky_idx = 0*px+0B
;
RESTORE,'csdr$firas_in:rhs.iss'
sky_idx = [sky_idx,0*px+1B]
;
RESTORE,'csdr$firas_in:lhf.iss'
sky_idx = [sky_idx,0*px+2B]
;
RESTORE,'csdr$firas_in:rhf.iss'
sky_idx = [sky_idx,0*px+3B]
;

sname = 'csdr$firas_in:high_zodispec.iss'
SAVE,filename=sname,chanscan,chan_label,sky_idx,zodi_freq,$
                    zodi_spec,zodi_model,zodimodpars
PRINT,' '
PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_ZODISPEC.ISS" Created.'
PRINT,' '
;

inputdatass=0 & zoditime=0 & tm=0 & nifgs=0 & glon=0 & glat=0 & scan=0 & time=0
st_sub=0 & solution=0 & f=0 & galcut=0 & cal_nifgs=0 & cal_tm=0 & xcal=0
ical=0 & refh=0 & skyh=0 & dihd=0 & sky_glitch=0 & cal_glitch=0 & sky_wgts=0
cal_wgts=0 & sky_dihd=0 & sky_s0=0 & cal_s0=0 & cal_lbl=0 & sky_lbl=0
sl_weight_rat=0 & fsl_idx=0
n_lhs=0 & b_lhs=0 & l_lhs=0 & d_lhs=0 & n_rhs=0 & b_rhs=0 & l_rhs=0 & d_rhs=0
n_lhf=0 & b_lhf=0 & l_lhf=0 & d_lhf=0 & n_rhf=0 & b_rhf=0 & l_rhf=0 & d_rhf=0
; 

RETURN
END
