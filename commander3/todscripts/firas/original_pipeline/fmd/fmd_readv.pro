Pro FMD_ReadV,idd,error
;

;  FMD_READ creates IDL save sets of coadd calibrated variances, by
;  driving IDL procedures which read data from the FSL_SKY and FSL_CAL files .
;
;
;  ARGUMENTS :  IDD       (I) = Channel/Scan-Mode Identifier (e.g. 'RHS')
;               ERROR     (O) = Return Status
;
;
;  PROGRAMS CALLED   :   FMD_VARIANCE
;                        FMD_VARIANCE2
;
;
;  EXAMPLE   :  $ upridl
;               UPRIDL> FMD_READV,'RHS',error
;
;
;  Required Logicals :
;
;   CSDR$FIRAS_REF    == Directory containing FMD_ARCV_TIMES_xxx.TXT files
;
;   CSDR$FIRAS_INSKY  == Directory containing FSL_SKY_xxxx files.
;
;   CSDR$FIRAS_INCAL  == Directory containing FSL_CAL_xxxx files.
;
;   CSDR$FIRAS_IN     == Directory where xxx.ISS, xxx_CALSPEC.ISS,
;                        and xxx_VAR.ISS will be written.
;
;
;  Written by :  Ken Jensen,  Hughes STX,  05-May-97
;
;

; Initialize Return Status
; ------------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
IF N_Params() ne 2 THEN BEGIN
 PRINT,'FMD_ReadV,idd,error'
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
ret = TRNLOG('csdr$firas_insky',inskytrans,/full,/issue_error)
inskytrans = STRUPCASE(inskytrans)
;
ret = TRNLOG('csdr$firas_incal',incaltrans,/full,/issue_error)
incaltrans = STRUPCASE(incaltrans)
;

PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,'CSDR$FIRAS_INSKY   == ' + inskytrans
PRINT,'CSDR$FIRAS_INCAL   == ' + incaltrans
PRINT,'CSDR$FIRAS_IN     == ' + intrans
PRINT,' '
;

; Default Options
; ---------------
splen = 320
galcut_hi = 20
galcut_lo = 15
arc_file = 'fmd_arcv_times_ssl'
;

IF (idd eq 'RHF') THEN BEGIN
 ;
 arc_file = ['fmd_arcv_times_lfl','fmd_arcv_times_lfs']
 chsc = ['rhlf','rhsf']
 galcut = galcut_hi
 ;
 galcut_rhlf = galcut
 FMD_VARIANCE,arc_file(0),chsc(0),splen,d_rhlf,n_rhlf,f_rhlf,var_rhlf,var_tm,$
              cal_var_rhlf,cal_var_tm,st_sub,galcut=galcut_rhlf
 ;
 galcut_rhsf = galcut
 FMD_VARIANCE,arc_file(1),chsc(1),splen,d_rhsf,n_rhsf,f_rhsf,var_rhsf,var_tm,$
              cal_var_rhsf,cal_var_tm,st_sub,galcut=galcut_rhsf
 ;
 PRINT,' '
 ;
 sl_weight_rat = TOTAL(1/(d_rhsf^2))/TOTAL(1/(d_rhlf^2))
 ;
 nl_adj = n_rhlf / SQRT(sl_weight_rat)
 ns_adj = n_rhsf * SQRT(sl_weight_rat)
 d_rhf = SQRT((nl_adj*(d_rhlf^2) + ns_adj*(d_rhsf^2))/(nl_adj+ns_adj))
 ;
 FMD_VARIANCE2,arc_file,chsc,splen,f_rhf,var_rhf,var_tm,var_lbl,$
               cal_var_rhf,cal_var_tm,cal_var_lbl,st_sub,galcut=galcut
 ;
 sname='csdr$firas_in:rhf_var.iss'
 SAVE,filename=sname,galcut_rhlf,galcut_rhsf,f_rhlf,f_rhsf,d_rhlf,d_rhsf,$
                     d_rhf,sl_weight_rat,var_rhlf,var_rhsf,cal_var_rhlf,$
                     cal_var_rhsf,f_rhf,var_rhf,cal_var_rhf,var_tm,$
                     cal_var_tm,var_lbl,cal_var_lbl,galcut
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RHF_VAR.ISS" Created.'
 ;
 PRINT,' '
 ;
ENDIF
;

IF (idd eq 'LHF') THEN BEGIN
 ;
 arc_file = ['fmd_arcv_times_lfl','fmd_arcv_times_lfs']
 chsc = ['lhlf','lhsf']
 galcut = galcut_hi
 ;
 galcut_lhlf = galcut
 FMD_VARIANCE,arc_file(0),chsc(0),splen,d_lhlf,n_lhlf,f_lhlf,var_lhlf,var_tm,$
              cal_var_lhlf,cal_var_tm,st_sub,galcut=galcut_lhlf
 ;
 galcut_lhsf = galcut
 FMD_VARIANCE,arc_file(1),chsc(1),splen,d_lhsf,n_lhsf,f_lhsf,var_lhsf,var_tm,$
              cal_var_lhsf,cal_var_tm,st_sub,galcut=galcut_lhsf
 ;
 PRINT,' '
 ;
 sl_weight_rat = TOTAL(1/(d_lhsf^2))/TOTAL(1/(d_lhlf^2))
 ;
 nl_adj = n_lhlf / SQRT(sl_weight_rat)
 ns_adj = n_lhsf * SQRT(sl_weight_rat)
 d_lhf = SQRT((nl_adj*(d_lhlf^2) + ns_adj*(d_lhsf^2))/(nl_adj+ns_adj))
 ;
 FMD_VARIANCE2,arc_file,chsc,splen,f_lhf,var_lhf,var_tm,var_lbl,$
               cal_var_lhf,cal_var_tm,cal_var_lbl,st_sub,galcut=galcut
 ;
 sname='csdr$firas_in:lhf_var.iss'
 SAVE,filename=sname,galcut_lhlf,galcut_lhsf,f_lhlf,f_lhsf,d_lhlf,d_lhsf,$
                     d_lhf,sl_weight_rat,var_lhlf,var_lhsf,cal_var_lhlf,$
                     cal_var_lhsf,f_lhf,var_lhf,cal_var_lhf,var_tm,$
                     cal_var_tm,var_lbl,cal_var_lbl,galcut
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LHF_VAR.ISS" Created.'
 PRINT,' '
 ;
ENDIF
;

; Return Status
; -------------
error = 0
PRINT,' '
PRINT,'FMD_READV : Successful Completion'
PRINT,' '
;

RETURN
END
