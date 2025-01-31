Pro FMD_CONCAT2_HR,error
;
;  FMD_CONCAT2_HR concatenates C-Vector, skymap spectra, cal coadd residuals, 
;  and coadd chi-squared for the three "bands" of combined HRES destriped data.
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
;    CSDR$FIRAS_OUT   =  Directory containing HRES_n_CVECTOR.ISS,
;                        HRES_n_CAL_CHI2.ISS, HRES_n_CAL_RESID.ISS,
;                        HRES_n_SKY_CHI2.ISS, and HRES_n_SKYMAP_F.ISS
;                        or HRES_n_SKYMAP_2.ISS, and where IDL save sets 
;                        HRES_CVECTOR.ISS, HRES_SKYMAP_F.ISS
;                        or HRES_SKYMAP_2.ISS, HRES_SKY_CHI2.ISS,
;                        HRES_CAL_RESID.ISS, and HRES_CAL_CHI2.ISS
;                        will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CONCAT2_HR,band,error
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
 PRINT,'FMD_CONCAT2_HR : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CONCAT2_HR,error'
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

; CVECTOR
; -------
RESTORE,'csdr$firas_out:hres_1_cvector.iss'
freq = TEMPORARY(freq_hres)
cvec = TEMPORARY(cvec_hres)
;
RESTORE,'csdr$firas_out:hres_2_cvector.iss'
freq = [freq,TEMPORARY(freq_hres)]
cvec = [cvec,TEMPORARY(cvec_hres)]
;
RESTORE,'csdr$firas_out:hres_3_cvector.iss'
freq = [freq,TEMPORARY(freq_hres)]
cvec = [cvec,TEMPORARY(cvec_hres)]
;

; SKY_CHI2
; --------
RESTORE,'csdr$firas_out:hres_1_sky_chi2.iss'
sky_chi2 = TRANSPOSE(TEMPORARY(sky_chi2_hres))
;
RESTORE,'csdr$firas_out:hres_2_sky_chi2.iss'
sky_chi2 = [[sky_chi2],[TRANSPOSE(TEMPORARY(sky_chi2_hres))]]
;
RESTORE,'csdr$firas_out:hres_3_sky_chi2.iss'
sky_chi2 = [[sky_chi2],[TRANSPOSE(TEMPORARY(sky_chi2_hres))]]
;

; CAL_CHI2
; --------
RESTORE,'csdr$firas_out:hres_1_cal_chi2.iss'
cal_chi2 = TRANSPOSE(TEMPORARY(cal_chi2_hres))
;
RESTORE,'csdr$firas_out:hres_2_cal_chi2.iss'
cal_chi2 = [[cal_chi2],[TRANSPOSE(TEMPORARY(cal_chi2_hres))]]
;
RESTORE,'csdr$firas_out:hres_3_cal_chi2.iss'
cal_chi2 = [[cal_chi2],[TRANSPOSE(TEMPORARY(cal_chi2_hres))]]
;

; CAL_RESID
; ---------
RESTORE,'csdr$firas_out:hres_1_cal_resid.iss'
cal_resid = TRANSPOSE(TEMPORARY(cal_resid_hres))
;
RESTORE,'csdr$firas_out:hres_2_cal_resid.iss'
cal_resid = [[cal_resid],[TRANSPOSE(TEMPORARY(cal_resid_hres))]]
;
RESTORE,'csdr$firas_out:hres_3_cal_resid.iss'
cal_resid = [[cal_resid],[TRANSPOSE(TEMPORARY(cal_resid_hres))]]
;

; SKYMAP
; ------
sp = FLTARR(96,64,182)
;
IF (freq_dependence EQ 'Y') THEN BEGIN
 ;
 RESTORE,'csdr$firas_out:hres_1_skymap_f.iss'
 FOR i=0,59 DO sp(*,*,i) = TEMPORARY(sp_hres(*,*,i))
 ;
 RESTORE,'csdr$firas_out:hres_2_skymap_f.iss'
 FOR i=0,59 DO sp(*,*,i+60) = TEMPORARY(sp_hres(*,*,i))
 ;
 RESTORE,'csdr$firas_out:hres_3_skymap_f.iss'
 FOR i=0,61 DO sp(*,*,i+120) = TEMPORARY(sp_hres(*,*,i))
 ;
ENDIF
;
IF (freq_dependence NE 'Y') THEN BEGIN
 ;
 RESTORE,'csdr$firas_out:hres_1_skymap_2.iss'
 FOR i=0,59 DO sp(*,*,i) = TEMPORARY(sp_hres(*,*,i))
 ;
 RESTORE,'csdr$firas_out:hres_2_skymap_2.iss'
 FOR i=0,59 DO sp(*,*,i+60) = TEMPORARY(sp_hres(*,*,i))
 ;
 RESTORE,'csdr$firas_out:hres_3_skymap_2.iss'
 FOR i=0,61 DO sp(*,*,i+120) = TEMPORARY(sp_hres(*,*,i))
 ;
ENDIF
;

freq_hres = TEMPORARY(freq)
cvec_hres = TEMPORARY(cvec)
sky_chi2_hres = TRANSPOSE(TEMPORARY(sky_chi2))
cal_chi2_hres = TRANSPOSE(TEMPORARY(cal_chi2))
cal_resid_hres = TRANSPOSE(TEMPORARY(cal_resid))
sp_hres = TEMPORARY(sp)
;

; Make IDL Save Sets
; ------------------
;
chanscan = 'HRES_'
chan_label = 'LLF_RLF'
;

IF (freq_dependence EQ 'Y') THEN BEGIN
 ;
 ; Make HRES_SKYMAP_F.ISS Save Set
 ; -------------------------------
 sname = 'csdr$firas_out:hres_skymap_f.iss'
 SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,$
                     freq_hres,sp_hres,cvec_hres
 PRINT,' '
 PRINT,'IDL Save Set "'+outtrans(0)+'HRES_SKYMAP_F.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (freq_dependence NE 'Y') THEN BEGIN
 ;
 ; Make HRES_SKYMAP_2.ISS Save Set
 ; -------------------------------
 sname = 'csdr$firas_out:hres_skymap_2.iss'
 SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,$
                     freq_hres,sp_hres,cvec_hres
 PRINT,' '
 PRINT,'IDL Save Set "'+outtrans(0)+'HRES_SKYMAP_2.ISS" Created.'
 PRINT,' '
;
ENDIF
;

; Make HRES_CVECTOR.ISS Save Set
; ------------------------------
sname = 'csdr$firas_out:hres_cvector.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,cvec_cut,$
                    freq_hres,cvec_hres,ndf_cvec_hres
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_CVECTOR.ISS" Created.'
PRINT,' '
;

; Make HRES_SKY_CHI2.ISS Save Set
; -------------------------------
sname = 'csdr$firas_out:hres_sky_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,$
                    sky_idx,sky_wgts_ds,frac_wgt,sky_mask,freq_hres,$
                    cvec_hres,sky_chi2_hres
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_SKY_CHI2.ISS" Created.'
PRINT,' '
;

; Make HRES_CAL_RESID.ISS Save Set
; --------------------------------
sname = 'csdr$firas_out:hres_cal_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,cal_idx,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_hres,cvec_hres,$
                    cal_resid_hres
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_CAL_RESID.ISS" Created.'
PRINT,' '
;

; Make HRES_CAL_CHI2.ISS Save Set
; -------------------------------
sname = 'csdr$firas_out:hres_cal_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,cal_idx,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_hres,cvec_hres,$
                    cal_chi2_hres
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_CAL_CHI2.ISS" Created.'
PRINT,' '
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
