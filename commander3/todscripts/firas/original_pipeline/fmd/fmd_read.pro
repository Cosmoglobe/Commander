Pro FMD_Read,idd,error
;

;  FMD_READ creates IDL save sets of coadd data, coadd calibrated variances,
;  and coadd calibrated spectra, by driving IDL procedures which read data
;  from the FSL_SKY and FSL_CAL files .
;
;
;  ARGUMENTS :  IDD       (I) = Channel/Scan-Mode Identifier (e.g. 'RHS')
;               ERROR     (O) = Return Status
;
;
;  PROGRAMS CALLED   :   FMD_VARIANCE
;                        FMD_VARIANCE2
;                        FMD_READ_DATA
;
;
;  EXAMPLE   :  $ upridl
;               UPRIDL> FMD_READ,'RHS',error
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
;  Written by :  Ken Jensen,  Hughes STX,  23-Apr-97
;
;  Modifications : K.Jensen, HSTX, 05-May-97 ;  if LHF or RHF, do not
;                  call FMD_VARIANCE2 and do not save xxx_VAR.ISS .
;

; Initialize Return Status
; ------------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
IF N_Params() ne 2 THEN BEGIN
 PRINT,'FMD_Read,idd,error'
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

IF (idd eq 'RHS') THEN BEGIN
 ;
 galcut = galcut_hi
 ;
 galcut_rhss = galcut
 FMD_VARIANCE,arc_file,'rhss',splen,d_rhss,n_rhss,f_rhss,var_rhs,var_tm,$
                        cal_var_rhs,cal_var_tm,st_sub,galcut=galcut_rhss
;
 f_rhs = f_rhss  &  d_rhs = d_rhss  &  galcut_rhs = galcut_rhss
 sname = 'csdr$firas_in:rhs_var.iss'
 SAVE,filename=sname,f_rhs,d_rhs,var_tm,var_rhs,cal_var_tm,$
                     cal_var_rhs,galcut_rhs
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RHS_VAR.iss" Created.'
 PRINT,' '
 ;
 FMD_READ_DATA,arc_file,'rhss',px,tm,sp, $
 	       nifgs,glon,glat,solution,f,n_rhs,b_rhs,l_rhs, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:rhs.iss',px,tm,nifgs,glon,glat,scan,time, $
                st_sub,solution,f,n_rhs,b_rhs,l_rhs,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,fsl_idx, $
     	        sky_wgts,cal_wgts,sky_dihd,sky_s0,cal_s0,cal_lbl,sky_lbl
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RHS.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:rhs_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RHS_CALSPEC.ISS" Created.'
 PRINT,' '
 ;
ENDIF
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
 PRINT,'SL_WEIGHT_RAT =',sl_weight_rat
 ;
 FMD_READ_DATA,arc_file,chsc,px,tm,sp, $
 	       nifgs,glon,glat,solution,f,n_rhf,b_rhf,l_rhf, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
               weight_cor=sl_weight_rat, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:rhf.iss',px,tm,nifgs,glon,glat,scan,time, $
                st_sub,solution,f,n_rhf,b_rhf,l_rhf,d_rhf,galcut,cal_nifgs, $
                cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch, $
                sky_wgts,fsl_idx,cal_wgts,sky_dihd,sky_s0,cal_s0, $
 	        cal_lbl,sky_lbl,sl_weight_rat
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RHF.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:rhf_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RHF_CALSPEC.ISS" Created.'
 PRINT,' '
 ;
ENDIF
;


IF (idd eq 'LHS') THEN BEGIN
 ;
 galcut = galcut_hi
 ;
 galcut_lhss = galcut
 FMD_VARIANCE,arc_file,'lhss',splen,d_lhss,n_lhss,f_lhss,var_lhs,var_tm,$
                        cal_var_lhs,cal_var_tm,st_sub,galcut=galcut_lhss
 ;
 f_lhs = f_lhss  &  d_lhs = d_lhss  &  galcut_lhs = galcut_lhss
 sname = 'csdr$firas_in:lhs_var.iss'
 SAVE,filename=sname,f_lhs,d_lhs,var_tm,var_lhs,cal_var_tm,$
                     cal_var_lhs,galcut_lhs
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LHS_VAR.ISS" Created.'
 ;
 FMD_READ_DATA,arc_file,'lhss',px,tm,sp, $
  	       nifgs,glon,glat,solution,f,n_lhs,b_lhs,l_lhs, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:lhs.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_lhs,b_lhs,l_lhs,d_lhs,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
                fsl_idx,cal_wgts,sky_dihd,sky_s0,cal_s0,cal_lbl,sky_lbl
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LHS.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:lhs_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LHS_CALSPEC.ISS" Created.'
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
 PRINT,'SL_WEIGHT_RAT =',sl_weight_rat
 ;
 FMD_READ_DATA,arc_file,chsc,px,tm,sp, $
 	       nifgs,glon,glat,solution,f,n_lhf,b_lhf,l_lhf, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
               weight_cor=sl_weight_rat, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:lhf.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_lhf,b_lhf,l_lhf,d_lhf,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
                fsl_idx,cal_wgts,sky_dihd,sky_s0,cal_s0, $
   	        cal_lbl,sky_lbl,sl_weight_rat
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LHF.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:lhf_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LHF_CALSPEC.ISS" Created.'
 ;
ENDIF
;


IF (idd eq 'RLS') THEN BEGIN
 ;
 galcut = galcut_lo
 ;
 galcut_rlss = galcut
 FMD_VARIANCE,arc_file,'rlss',splen,d_rlss,n_rlss,f_rlss,var_rls,var_tm,$
                        cal_var_rls,cal_var_tm,st_sub,galcut=galcut_rlss
 ;
 f_rls = f_rlss  &  d_rls = d_rlss  &  galcut_rls = galcut_rlss
 sname = 'csdr$firas_in:rls_var.iss'
 SAVE,filename=sname,f_rls,d_rls,var_tm,var_rls,cal_var_tm,$
                     cal_var_rls,galcut_rls
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RLS_VAR.iss" Created.'
 PRINT,' '
 ;
 FMD_READ_DATA,arc_file,'rlss',px,tm,sp, $
  	       nifgs,glon,glat,solution,f,n_rls,b_rls,l_rls, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:rls.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_rls,b_rls,l_rls,d_rls,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,fsl_idx, $
   	        sky_wgts,cal_wgts,sky_dihd,sky_s0,cal_s0,cal_lbl,sky_lbl
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RLS.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:rls_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RLS_CALSPEC.ISS" Created.'
 PRINT,' '
 ;
ENDIF
;

IF (idd eq 'RLF') THEN BEGIN
 ;
 arc_file = 'fmd_arcv_times_lfl'
 galcut = galcut_lo
 ;
 galcut_rllf = galcut
 FMD_VARIANCE,arc_file,'rllf',splen,d_rllf,n_rllf,f_rllf,var_rlf,var_tm,$
                        cal_var_rlf,cal_var_tm,st_sub,galcut=galcut_rllf
 ;
 f_rlf = f_rllf  &  d_rlf = d_rllf  &  galcut_rlf = galcut_rllf
 sname = 'csdr$firas_in:rlf_var.iss'
 SAVE,filename=sname,f_rlf,d_rlf,var_tm,var_rlf,cal_var_tm,$
                     cal_var_rlf,galcut_rlf
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RLF_VAR.ISS" Created.'
 PRINT,' '
 ;
 FMD_READ_DATA,arc_file,'rllf',px,tm,sp, $
  	       nifgs,glon,glat,solution,f,n_rlf,b_rlf,l_rlf, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:rlf.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_rlf,b_rlf,l_rlf,d_rlf,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,fsl_idx, $
   	        sky_wgts,cal_wgts,sky_dihd,sky_s0,cal_s0,cal_lbl,sky_lbl
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RLF.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:rlf_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RLF_CALSPEC.ISS" Created.'
 ;
ENDIF
;


IF (idd eq 'RSF') THEN BEGIN
 ;
 splen = 80
 arc_file = ['fmd_arcv_times_lfl','fmd_arcv_times_lfs']
 chsc = ['rlfl','rlfs']
 galcut = galcut_lo
 ;
 galcut_rlfl = galcut
 FMD_VARIANCE,arc_file(0),chsc(0),splen,d_rlfl,n_rlfl,f_rlfl,var_rlfl,var_tm,$
              cal_var_rlfl,cal_var_tm,st_sub,galcut=galcut_rlfl
 ;
 galcut_rlfs = galcut
 FMD_VARIANCE,arc_file(1),chsc(1),splen,d_rlfs,n_rlfs,f_rlfs,var_rlfs,var_tm,$
              cal_var_rlfs,cal_var_tm,st_sub,galcut=galcut_rlfs
 ;
 PRINT,' '
 ;
 sl_weight_rat = TOTAL(1/(d_rlfs^2))/TOTAL(1/(d_rlfl^2))
 ;
 nl_adj = n_rlfl / SQRT(sl_weight_rat)
 ns_adj = n_rlfs * SQRT(sl_weight_rat)
 d_rsf = SQRT((nl_adj*(d_rlfl^2) + ns_adj*(d_rlfs^2))/(nl_adj+ns_adj))
 ;
 FMD_VARIANCE2,arc_file,chsc,splen,f_rsf,var_rsf,var_tm,var_lbl,$
               cal_var_rsf,cal_var_tm,cal_var_lbl,st_sub,galcut=galcut
 ;
 sname='csdr$firas_in:rsf_var.iss'
 SAVE,filename=sname,galcut_rlfl,galcut_rlfs,f_rlfl,f_rlfs,d_rlfl,d_rlfs,$
                     d_rsf,sl_weight_rat,var_rlfl,var_rlfs,cal_var_rlfl,$
                     cal_var_rlfs,f_rsf,var_rsf,cal_var_rsf,var_tm,$
                     cal_var_tm,var_lbl,cal_var_lbl,galcut
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RSF_VAR.ISS" Created.'
 PRINT,' '
 ;
 PRINT,'SL_WEIGHT_RAT =',sl_weight_rat
 ;
 FMD_READ_DATA,arc_file,chsc,px,tm,sp, $
  	       nifgs,glon,glat,solution,f,n_rsf,b_rsf,l_rsf, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
               weight_cor=sl_weight_rat, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:rsf.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_rsf,b_rsf,l_rsf,d_rsf,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
                fsl_idx,cal_wgts,sky_dihd,sky_s0,cal_s0,$
     	        cal_lbl,sky_lbl,sl_weight_rat
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RSF.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:rsf_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'RSF_CALSPEC.ISS" Created.'
 PRINT,' '
 ;
ENDIF
;


IF (idd eq 'LLS') THEN BEGIN
 ;
 galcut = galcut_lo
 ;
 galcut_llss = galcut
 FMD_VARIANCE,arc_file,'llss',splen,d_llss,n_llss,f_llss,var_lls,var_tm,$
                        cal_var_lls,cal_var_tm,st_sub,galcut=galcut_llss
 ;
 f_lls = f_llss  &  d_lls = d_llss  &  galcut_lls = galcut_llss
 sname = 'csdr$firas_in:lls_var.iss'
 SAVE,filename=sname,f_lls,d_lls,var_tm,var_lls,cal_var_tm,$
                     cal_var_lls,galcut_lls
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LLS_VAR.ISS" Created.'
 PRINT,' '
 ;
 FMD_READ_DATA,arc_file,'llss',px,tm,sp, $
  	       nifgs,glon,glat,solution,f,n_lls,b_lls,l_lls, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:lls.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_lls,b_lls,l_lls,d_lls,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,fsl_idx, $
    	        sky_wgts,cal_wgts,sky_dihd,sky_s0,cal_s0,cal_lbl,sky_lbl
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LLS.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:lls_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LLS_CALSPEC.ISS" Created.'
 PRINT,' '
 ;
ENDIF
;


IF (idd eq 'LLF') THEN BEGIN
 ;
 arc_file = 'fmd_arcv_times_lfl'
 galcut = galcut_lo
 ;
 galcut_lllf = galcut
 FMD_VARIANCE,arc_file,'lllf',splen,d_lllf,n_lllf,f_lllf,var_llf,var_tm,$
                        cal_var_llf,cal_var_tm,st_sub,galcut=galcut_lllf
 ;
 f_llf = f_lllf  &  d_llf = d_lllf  &  galcut_llf = galcut_lllf
 sname = 'csdr$firas_in:llf_var.iss'
 SAVE,filename=sname,f_llf,d_llf,var_tm,var_llf,cal_var_tm,$
                     cal_var_llf,galcut_llf
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LLF_VAR.ISS" Created.'
 ;
 FMD_READ_DATA,arc_file,'lllf',px,tm,sp, $
  	       nifgs,glon,glat,solution,f,n_llf,b_llf,l_llf, $
 	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:llf.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_llf,b_llf,l_llf,d_llf,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
                fsl_idx,cal_wgts,sky_dihd,sky_s0,cal_s0,cal_lbl,sky_lbl
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LLF.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:llf_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LLF_CALSPEC.ISS" Created.'
 PRINT,' '
 ;
ENDIF
;


IF (idd eq 'LSF') THEN BEGIN
 ;
 splen = 80
 arc_file = ['fmd_arcv_times_lfl','fmd_arcv_times_lfs']
 chsc = ['llfl','llfs']
 galcut = galcut_lo
 ;
 galcut_llfl = galcut
 FMD_VARIANCE,arc_file(0),chsc(0),splen,d_llfl,n_llfl,f_llfl,var_llfl,var_tm,$
              cal_var_llfl,cal_var_tm,st_sub,galcut=galcut_llfl
 ;
 galcut_llfs = galcut
 FMD_VARIANCE,arc_file(1),chsc(1),splen,d_llfs,n_llfs,f_llfs,var_llfs,var_tm,$
              cal_var_llfs,cal_var_tm,st_sub,galcut=galcut_llfs
 ;
 sl_weight_rat = TOTAL(1/(d_llfs^2))/TOTAL(1/(d_llfl^2))
 ;
 nl_adj = n_llfl / SQRT(sl_weight_rat)
 ns_adj = n_llfs * SQRT(sl_weight_rat)
 d_lsf = SQRT((nl_adj*(d_llfl^2) + ns_adj*(d_llfs^2))/(nl_adj+ns_adj))
 ;
 FMD_VARIANCE2,arc_file,chsc,splen,f_lsf,var_lsf,var_tm,var_lbl,$
               cal_var_lsf,cal_var_tm,cal_var_lbl,st_sub,galcut=galcut
 ;
 sname='csdr$firas_in:lsf_var.iss'
 SAVE,filename=sname,galcut_llfl,galcut_llfs,f_llfl,f_llfs,d_llfl,d_llfs,$
                     d_lsf,sl_weight_rat,var_llfl,var_llfs,cal_var_llfl,$
                     cal_var_llfs,f_lsf,var_lsf,cal_var_lsf,var_tm,$
                     cal_var_tm,var_lbl,cal_var_lbl,galcut
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LSF_VAR.iss" Created.'
 PRINT,' '
 ;
 PRINT,'SL_WEIGHT_RAT =',sl_weight_rat
 ;
 FMD_READ_DATA,arc_file,chsc,px,tm,sp, $
 	       nifgs,glon,glat,solution,f,n_lsf,b_lsf,l_lsf, $
	       cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
               dihd,galcut=galcut,sky_glitch=sky_glitch, $
	       cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	       sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
               weight_cor=sl_weight_rat, $
	       sky_s0=sky_s0,cal_s0=cal_s0,time=time, $
               sky_dihd=sky_dihd,scan=scan,fsl_idx=fsl_idx
 ;
 SAVE,filename='csdr$firas_in:lsf.iss',px,tm,nifgs,glon,glat,scan,time,st_sub, $
                solution,f,n_lsf,b_lsf,l_lsf,d_lsf,galcut,cal_nifgs,cal_tm, $
                xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
                fsl_idx,cal_wgts,sky_dihd,sky_s0,cal_s0, $
   	        cal_lbl,sky_lbl,sl_weight_rat
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LSF.ISS" Created.'
 ;
 SAVE,filename='csdr$firas_in:lsf_calspec.iss',f,px,tm,sp,cal_tm,cal_sp
 ;
 PRINT,'IDL Save Set "' + intrans(0) + 'LSF_CALSPEC.ISS" Created.'
 PRINT,' '
 ;
ENDIF
;

; Return Status
; -------------
error = 0
PRINT,' '
PRINT,'FMD_READ : Successful Completion'
PRINT,' '
;

RETURN
END
