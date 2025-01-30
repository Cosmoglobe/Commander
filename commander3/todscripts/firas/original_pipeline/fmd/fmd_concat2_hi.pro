Pro FMD_CONCAT2_HI,error
;
;  FMD_CONCAT2_HI concatenates C-Vector, skymap spectra, zodi spectra,
;  and coadd chi-squared for the three bands of combined HIGH destriped
;  data, and stores them in IDL save sets .
;
;
;  ARGUMENTS (I/O)    :  ERROR (O)   :  Return Error Status
;
;
;  PROGRAMS Called    :  None
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_n_SKYMAP.ISS,
;                        HIGH_n_SKY_CHI2.ISS, and HIGH_n_CAL_CHI2.ISS,
;                        and where IDL save sets HIGH_SKYMAP.ISS,
;                        HIGH_SKY_CHI2.ISS, and HIGH_CAL_CHI2.ISS
;                        will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CONCAT2_HI,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 15-May-1997
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
 PRINT,'FMD_CONCAT2_HI : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CONCAT2_HI,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 RETURN
ENDIF
;

; Logical translations
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

; Skymap Spectra, Zodi Spectra, and C-Vector
; ------------------------------------------
RESTORE,'csdr$firas_out:high_2_skymap.iss'
freq = TEMPORARY(freq_high_2)
cvec = TEMPORARY(c_high_2)
;
RESTORE,'csdr$firas_out:high_3_skymap.iss'
freq = [freq,TEMPORARY(freq_high_3)]
cvec = [cvec,TEMPORARY(c_high_3)]
;
RESTORE,'csdr$firas_out:high_4_skymap.iss'
freq_high = [freq,TEMPORARY(freq_high_4)]
c_high = [cvec,TEMPORARY(c_high_4)]
;

s_high = DBLARR(96,64,170)
FOR i=0,54 DO s_high(*,*,i) = s_high_2(*,*,i)
FOR i=0,54 DO s_high(*,*,i+55) = s_high_3(*,*,i)
FOR i=0,59 DO s_high(*,*,i+110) = s_high_4(*,*,i)
;

z_high = FLTARR(96,64,170)
FOR i=0,54 DO z_high(*,*,i) = z_high_2(*,*,i)
FOR i=0,54 DO z_high(*,*,i+55) = z_high_3(*,*,i)
FOR i=0,59 DO z_high(*,*,i+110) = z_high_4(*,*,i)
;

; SKY_CHI2
; --------
RESTORE,'csdr$firas_out:high_2_sky_chi2.iss'
sky_chi2 = TRANSPOSE(TEMPORARY(sky_chi2_high_2))
cvec = TEMPORARY(cvec_high_2)
;
RESTORE,'csdr$firas_out:high_3_sky_chi2.iss'
sky_chi2 = [[sky_chi2],[TRANSPOSE(TEMPORARY(sky_chi2_high_3))]]
cvec = [cvec,TEMPORARY(cvec_high_3)]
;
RESTORE,'csdr$firas_out:high_4_sky_chi2.iss'
sky_chi2_high = TRANSPOSE([[sky_chi2],[TRANSPOSE(TEMPORARY(sky_chi2_high_4))]])
cvec_high = [cvec,TEMPORARY(cvec_high_4)]
;

; CAL_CHI2
; --------
RESTORE,'csdr$firas_out:high_2_cal_chi2.iss'
cal_chi2 = TRANSPOSE(TEMPORARY(cal_chi2_high_2))
;
RESTORE,'csdr$firas_out:high_3_cal_chi2.iss'
cal_chi2 = [[cal_chi2],[TRANSPOSE(TEMPORARY(cal_chi2_high_3))]]
;
RESTORE,'csdr$firas_out:high_4_cal_chi2.iss'
cal_chi2_high = TRANSPOSE([[cal_chi2],[TRANSPOSE(TEMPORARY(cal_chi2_high_4))]])
;

; Make IDL Save Sets
; ------------------
;
chanscan = 'HIGH'
chan_label = 'LHS_RHS_LHF_RHF'
;

; Make HIGH_SKYMAP.ISS Save Set
; -----------------------------
sname = 'csdr$firas_out:high_skymap.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,$
                    freq_high,s_high,z_high,c_high,n_high,b_high,l_high
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_SKYMAP.ISS" Created.'
PRINT,' '
;

; Make HIGH_SKY_CHI2.ISS Save Set
; -------------------------------
sname = 'csdr$firas_out:high_sky_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,$
                    sky_idx,sky_wgts_ds,frac_wgt,sky_mask,freq_high,$
                    cvec_high,sky_chi2_high
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_SKY_CHI2.ISS" Created.'
PRINT,' '
;

; Make HIGH_CAL_CHI2.ISS Save Set
; -------------------------------
sname = 'csdr$firas_out:high_cal_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,cal_idx,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_high,cvec_high,$
                    cal_chi2_high
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_CAL_CHI2.ISS" Created.'
PRINT,' '
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
