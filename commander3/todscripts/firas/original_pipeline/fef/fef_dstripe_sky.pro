;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FEF_DSTRIPE_SKY
;
;DESCRIPTION:  
;     This command file drives the 'FEF_Destriper' procedure which computes
;     the correction spectra for the calibrated input spectra.
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;   CHANSCAN  (I) =  Channel/Scan-Mode ID  (e.g. 'RHSS')  = String
;   TAU_ARRAY (I) =  UFO Time Constants (days)            = Float(nufo)
;   STEP_UP   (I) =  OffSet Period Start Times (Zulu)     = String(noff)
;   STEP_DN   (I) =  OffSet Period End Times (Zulu)       = String(noff)
;   S_xxx     (O) =  DeStriped Skymap (Eplees)            = Complex(96,64,nfreq)
;   N_xxx     (O) =  Skymap Weights                       = Float(96,64)
;   Cxxx      (O) =  Real Part of C-Vector (Eplees)       = Float(nfreq)
;   PCVR      (O) =  Q-Inverse                            = Float(ns,ns)
;   DIAG      (O) =  Diagonal Elements                    = Float(npix)
;   RECT      (O) =  Rectangular Elements                 = Float(npix,ns)
;   ERROR     (O) =  Error Status                         = 0 or 1
;
;WARNING:
;   The Logicals CSDR$FIRAS_IN and CSDR$FIRAS_SAVE must be previously defined.
;
;   CSDR$FIRAS_IN must contain the input xxx.ISS IDL Save Set.
;   CSDR$FIRAS_SAVE will contain the output FEF_xxx_y.ISS IDL Save Set,
;   ('xxx' identifies the channel/scan-mode, e.g. 'RHS' )
;   ('y' identifies the destriping model, e.g. '55_4MP' )
;
;EXAMPLE:
;   $ UPRIDL
;   chanscan='LLSS' & tau_array=[55.]
;   step_up = ['89328','901391535','901931850','902081120']
;   step_dn = ['901391535','901931850','902081120','9026409']
;   FEF_DSTRIPE_SKY,chanscan,tau_array,step_up,step_dn,,s_xxx,n_xxx,cxxx,error
;
;#
;COMMON BLOCKS:
;   None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;   Calls FEF_DESTRIPER.PRO
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;   None
;  
;MODIFICATION HISTORY:
;   Written by Ken Jensen, Hughes STX, Sept 14, 1994
;   Modified by Ken Jensen, Hughes STX, Sept 26, 1994 (SPR 11918)
; 
;-
;___________________________________________________________________
;
Pro FEF_DStripe_Sky,chanscan,tau_array,step_up,step_dn,s_xxx,n_xxx,cxxx,$
                    pcvr,diag,rect,error
;

; Set Error Status
; ----------------
error=1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 11 then begin
 print,'FEF_DStripe_Sky,chanscan,tau_array,step_up,step_dn,s_xxx,n_xxx,cxxx,pcvr,diag,rect,error'
 return
endif
;

; Channel/Scan-Mode Identifiers
; -----------------------------
chanscan=strupcase(chanscan)
if(chanscan eq 'RHSS')then idd='RHS'
if((chanscan eq 'RHSF')or(chanscan eq 'RHLF'))then idd='RHF'
if(chanscan eq 'RLSS')then idd='RLS'
if((chanscan eq 'RLSF')or(chanscan eq 'RLFL'))then idd='RSF'
if(chanscan eq 'RLLF')then idd='RLF'
if(chanscan eq 'LHSS')then idd='LHS'
if((chanscan eq 'LHSF')or(chanscan eq 'LHLF'))then idd='LHF'
if(chanscan eq 'LLSS')then idd='LLS'
if((chanscan eq 'LLSF')or(chanscan eq 'LLFL'))then idd='LSF'
if(chanscan eq 'LLLF')then idd='LLF'
;
id_array=chanscan
if((chanscan eq 'RHSF')or(chanscan eq 'RHLF'))then id_array=['rhsf','rhlf']
if((chanscan eq 'RLSF')or(chanscan eq 'RLFL'))then id_array=['rlsf','rlfl']
if((chanscan eq 'LHSF')or(chanscan eq 'LHLF'))then id_array=['lhsf','lhlf']
if((chanscan eq 'LLSF')or(chanscan eq 'LLLF'))then id_array=['llsf','llfl']
;

; Number of "Stripes"
; -------------------
nufo = n_elements(tau_array)
if(nufo lt 1)then begin
 print,'FEF_DSTRIPE_SKY : AT Least One TAU Required !!!'
 return
endif
;
noff = n_elements(step_up)
if(n_elements(step_dn) ne noff)then begin
 print,'FEF_DSTRIPE_SKY : STEP_UP and STEP_DN Have Unequal Size !!!'
 return
endif
;
nstripes = nufo + noff
;

; Correction Index (1=UFO, 2=Vib, 3=Offset)
; -----------------------------------------
corr_index = intarr(nstripes) + 3
corr_index(0:nufo-1) = 1
;

; Restore sky and cal spectra data
; ---------------------------------
PRINT,'Restoring Uncorrected SKY and CAL Spectra'
restore,'csdr$firas_in:'+idd+'.iss'
;

; Verify Galactic Plane Exclusion
; -------------------------------
if(galcut ne 10)then begin
 print,'!!! FEF_DSTRIPE_SKY : GALCUT NE 10 !!!'
 return
endif
;

; Screen out Cal Coadds, except for sky-like nulls
; ------------------------------------------------
cal_wgts_ds = cal_wgts
;
tl = 2.6
tu = 2.8
bad = WHERE((xcal LE  tl) OR (xcal GE tu) OR (ical LE 2.7) OR (ical GE tu)$
 OR (skyh LE tl) OR (skyh GE tu) OR (refh LE tl) OR (refh GE tu) OR (dihd GE 3.5))
;
cal_wgts_ds(bad) = 0.
;

; Good Sky Coadds
; ---------------
sky_wgts_ds = sky_wgts
;
ll = coorconv(px,infmt='p',inco='f',outfmt='l',outco='g')
bad = WHERE(ABS(ll(*,1)) LE galcut,cbad)
if(cbad gt 0)then sky_wgts_ds(bad) = 0.
;
list=[391,401,473,496,501,2139,2151,2160,2161,2163,2168,2169,2171,2208,2234]
pixcut=[list,2615,2616,2799,4096,4108,4244,4278,4284,4286,4628,4735,4821,5715]
npx=n_elements(pixcut)
for i=0,npx-1 do begin
 bad = WHERE(px eq pixcut(i),cx)
 if(cx gt 0)then sky_wgts_ds(bad) = 0.
endfor
;

; Define Inputs for Destriper
; ---------------------------
del_temp = 0.004               ; XCAL adjustment
index = FIX(0*f)               ; Offset Destriping
tau_xxx=tau_array
;

; Run Destriper for the appropriate ID
; ------------------------------------
print,'Calling FEF_DESTRIPER'
cxxx = FEF_DESTRIPER(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id=id_array, $
                      tm=tm,cal_tm=cal_tm,tau=tau_xxx,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_xxx,afp=s_xxx,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
;

; Save saveset
; ------------
fext=strtrim(strcompress(string(fix(tau_xxx(0)))),2)
if(nufo gt 1)then begin
 for i=1,nufo-1 do fext=fext+'_'+strtrim(strcompress(string(fix(tau_xxx(i)))),2)
endif
soff=strtrim(strcompress(string(noff)),2)
fext='_'+fext+'_'+soff+'mp'
;
sname='csdr$firas_save:fef_'+idd+fext+'.iss'
;
if(idd eq 'RHS')then begin
 tau_rhs=tau_xxx & ejv_rhs=ejv_xxx & crhs=cxxx & s_rhs=s_xxx & n_xxx=n_rhs
 save,filename=sname,idd,tau_rhs,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_rhs,ndf_ejv,crhs,ejv_rhs,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'RHF')then begin
 tau_rhf=tau_xxx & ejv_rhf=ejv_xxx & crhf=cxxx & s_rhf=s_xxx & n_xxx=n_rhf
 save,filename=sname,idd,tau_rhf,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_rhf,ndf_ejv,crhf,ejv_rhf,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'RLS')then begin
 tau_rls=tau_xxx & ejv_rls=ejv_xxx & crls=cxxx & s_rls=s_xxx & n_xxx=n_rls
 save,filename=sname,idd,tau_rls,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_rls,ndf_ejv,crls,ejv_rls,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'RSF')then begin
 tau_rsf=tau_xxx & ejv_rsf=ejv_xxx & crsf=cxxx & s_rsf=s_xxx & n_xxx=n_rsf
 save,filename=sname,idd,tau_rsf,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_rsf,ndf_ejv,crsf,ejv_rsf,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'RLF')then begin
 tau_rlf=tau_xxx & ejv_rlf=ejv_xxx & crlf=cxxx & s_rlf=s_xxx & n_xxx=n_rlf
 save,filename=sname,idd,tau_rlf,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_rlf,ndf_ejv,crlf,ejv_rlf,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'LHS')then begin
 tau_lhs=tau_xxx & ejv_lhs=ejv_xxx & clhs=cxxx & s_lhs=s_xxx & n_xxx=n_lhs
 save,filename=sname,idd,tau_lhs,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_lhs,ndf_ejv,clhs,ejv_lhs,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'LHF')then begin
 tau_lhf=tau_xxx & ejv_lhf=ejv_xxx & clhf=cxxx & s_lhf=s_xxx & n_xxx=n_lhf
 save,filename=sname,idd,tau_lhf,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_lhf,ndf_ejv,clhf,ejv_lhf,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'LLS')then begin
 tau_lls=tau_xxx & ejv_lls=ejv_xxx & clls=cxxx & s_lls=s_xxx & n_xxx=n_lls
 save,filename=sname,idd,tau_lls,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_lls,ndf_ejv,clls,ejv_lls,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'LSF')then begin
 tau_lsf=tau_xxx & ejv_lsf=ejv_xxx & clsf=cxxx & s_lsf=s_xxx & n_xxx=n_lsf
 save,filename=sname,idd,tau_lsf,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_lsf,ndf_ejv,clsf,ejv_lsf,sqr_mat,rect,diag,pcvr,corr_index,index
endif
if(idd eq 'LLF')then begin
 tau_llf=tau_xxx & ejv_llf=ejv_xxx & cllf=cxxx & s_llf=s_xxx & n_xxx=n_llf
 save,filename=sname,idd,tau_llf,step_up,step_dn,galcut,solution,sky_wgts_ds, $
      f,s_llf,ndf_ejv,cllf,ejv_llf,sqr_mat,rect,diag,pcvr,corr_index,index
endif
;

print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '

; Real Part of C-Vector
; ---------------------
cxxx = cxxx(*,0)
;

; Unused Parameters Redefined
; ---------------------------
nifgs=0 & glon=0 & glat=0 & cal_nifgs=0 & ical=0 & refh=0 & skyh=0 & dihd=0
sky_glitch=0 & cal_glitch=0 & cal_lbl=0 & sky_lbl=0 & sl_weight_rat=0
;
n_rhs=0 & b_rhs=0 & l_rhs=0 & urhs=0 & d_rhs=0
n_rhf=0 & b_rhf=0 & l_rhf=0 & urhf=0 & d_rhf=0 & d_rhsf=0 & d_rhlf=0
n_rls=0 & b_rls=0 & l_rls=0 & urls=0 & d_rls=0
n_rsf=0 & b_rsf=0 & l_rsf=0 & ursf=0 & d_rsf=0 & d_rlsf=0 & d_rlfl=0
n_rlf=0 & b_rlf=0 & l_rlf=0 & urlf=0 & d_rlf=0
n_lhs=0 & b_lhs=0 & l_lhs=0 & ulhs=0 & d_lhs=0
n_lhf=0 & b_lhf=0 & l_lhf=0 & ulhf=0 & d_lhf=0 & d_lhsf=0 & d_lhlf=0
n_lls=0 & b_lls=0 & l_lls=0 & ulls=0 & d_lls=0
n_lsf=0 & b_lsf=0 & l_lsf=0 & ulsf=0 & d_lsf=0 & d_llsf=0 & d_llfl=0
n_llf=0 & b_llf=0 & l_llf=0 & ullf=0 & d_llf=0
;

; Set Return Status to NO ERROR
; -----------------------------
error=0
;

return
end
