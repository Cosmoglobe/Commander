;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FDS_DSTRIPE_SKY
;
;DESCRIPTION:  
;     This command file drives the 'fds_destriper' procedure which computes
;     the correction spectra for the calibrated input spectra.
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;   CHANSCAN  (I) =  Channel/MTM Scan Mode (e.g. 'LLSS')         =  String
;   ERROR     (O) =  Error Status
;
;WARNING:
;     The Logicals CSDR$FIRAS_REF, CSDR$FIRAS_IN, CSDR$FIRAS_OUT, 
;     and CSDR$FIRAS_SAVE must be previously defined.
;
;     CSDR$FIRAS_REF should contain the reference data sets FEX_MOD_XXXX
;                    and FDS_PARAMS_XXXX.
;     CSDR$FIRAS_IN must contain the input XXX.ISS IDL Save Set.
;     CSDR$FIRAS_OUT will contain the output reference data set FEX_EJV_XXXX.
;     CSDR$FIRAS_SAVE will contain the output XXX_DESTRIPED.ISS IDL Save Set.
;
;EXAMPLE:
;      $ UIDL
;      chanscan = 'LLSS'
;      fds_dstripe_sky,chanscan,error
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;     Calls FDS_DESTRIPER.PRO
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;    Written by Ken Jensen, Hughes STX, June 17, 1994
;
;    Modification by Ken Jensen, Hughes STX, June 21, 1994,  SKY_WGTS_DS
;      set to zero for RHSS coadds with sky pixel numbers contained in
;      the list PIXCUT. Keyword SKY_WGTS_DS added to call to FDS_DESTRIPER.
;
;    Modification by Ken Jensen, Hughes STX, June 23, 1994, TAU_ARRAY,
;      VIB_FLAG, STEP_UP, and STEP_DN  hard-coded rather than passed
;      as input arguments.  SPR 11813.
;
;    Modification by Ken Jensen, Hughes STX, June 30, 1994, STEP_UP and
;      STEP_DN Arrays modified to make separate offset periods before
;      and after the hot horn periods.
;
;    Modification by Ken Jensen, Hughes STX, July 6, 1994, Parameters
;      TAU_ARRAY, VIB_FLAG, STEP_UP, STEP_DN read from a reference data set.
;
;    Modification by Ken Jensen, Hughes STX, Aug 9, 1994,  INDEX array
;      elements all set to zero to enable destriping at all frequencies,
;      (SER 11865)
;
;-
;___________________________________________________________________
;
Pro FDS_Dstripe_Sky,chanscan,error
;

; Set Error Status
; ----------------
error=1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 2 then begin
 print,'FDS_Dstripe_Sky,chanscan,error'
 return
endif
;

;  Define the FDS 3-letter ID (RHSS -> RHS etc.)
;  ---------------------------------------------
id = strupcase(strmid(chanscan,0,2))+strupcase(strmid(chanscan,2,2))
;
if((id ne 'RLSF')and(id ne 'LLSF'))then begin
 id = strupcase(strmid(chanscan,0,2))+strupcase(strmid(chanscan,3,1))
endif
;
if(id eq 'RLSF')then id = 'RSF'
if(id eq 'LLSF')then id = 'LSF'
;


; Destriping Parameters
; ---------------------
OPENR,1,'csdr$firas_ref:fds_params_'+chanscan+'.txt'
;
nufo = 0
READF,1,nufo
if(nufo lt 1)then begin
 print,'FDS_DSTRIPE_SKY : NUFO LT 1 !!!'
 return
endif
tau_array = fltarr(nufo)
for i=0,nufo-1 do begin
 aa = 0.
 READF,1,aa
 tau_array(i)=aa
endfor
;
vib_flag = ''
READF,1,vib_flag
;
noff = 0
READF,1,noff
if(noff lt 1)then begin
 print,'FDS_DSTRIPE_SKY : NOFF LT 1 !!!'
 return
endif
aa = ''
READF,1,aa
step_up = aa
if(noff gt 1)then begin
 for i=1,noff-1 do begin
  READF,1,aa
  step_up = [step_up,aa]
 endfor
endif
;
READF,1,aa
step_dn = aa
if(noff gt 1)then begin
 for i=1,noff-1 do begin
  READF,1,aa
  step_dn = [step_dn,aa]
 endfor
endif
;
CLOSE,1
;


; Number of "Stripes"
; -------------------
nufo = n_elements(tau_array)
if(nufo lt 1)then begin
 print,'FDS_DSTRIPE_SKY : AT Least One TAU Required !!!'
 return
endif
;
nvib = 1
if(vib_flag ne 'Y')then nvib=0
;
noff = n_elements(step_up)
if(n_elements(step_dn) ne noff)then begin
 print,'FDS_DSTRIPE_SKY : STEP_UP and STEP_DN Have Unequal Size !!!'
 return
endif
;
nstripes = nufo + nvib + noff


; Correction Index (1=UFO, 2=Vib, 3=Offset)
; -----------------------------------------
corr_index = intarr(nstripes) + 3
corr_index(0:nufo-1) = 1
if(vib_flag eq 'Y')then corr_index(nufo) = 2
;


; IF(VIB_FLAG) THEN Read vibration correction parameters
; ------------------------------------------------------
vib = dblarr(5)
if(nvib eq 1)then begin
 OPENR,1,'csdr$firas_ref:fex_mod_'+chanscan+'.f16_93hybrid'
 rec = DBLARR(32+20)
 READU,1,rec
 CLOSE,1
 vib = rec(32+15:32+19)
 vib = vib / max(abs(vib))  ;  Normalize Vibration Coefficients
endif


; Restore sky and cal spectra data
; ---------------------------------
PRINT,'Restoring Uncorrected SKY and CAL Spectra'
restore,'csdr$firas_in:'+id+'.iss'
;


; Verify Galactic Plane Exclusion
; -------------------------------
if(galcut ne 10)then begin
 print,'!!! FDS_DSTRIPE_SKY : GALCUT NE 10 !!!'
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

;
; Run Destriper and save data for the appropriate ID
; --------------------------------------------------
;
; RHS Data
; --------
if ( id eq 'RHS' ) then begin

 tau_rhs=tau_array
;

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 crhs = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rhss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rhs,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
    		      ndf_off=ndf_ejv,ejv=ejv_rhs,afp=s_rhs,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 crhs = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rhss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rhs,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
    		      ndf_off=ndf_ejv,ejv=ejv_rhs,afp=s_rhs,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:rhs_destriped.iss', $
              tau_rhs,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_rhs,csp,step_up,step_dn,ndf_ejv,crhs,ejv_rhs, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_rhs=0 & b_rhs=0 & l_rhs=0 & urhs=0 & d_rhs=0

endif
;
; RHF Data
; --------
if ( id eq 'RHF' ) then begin

 tau_rhf=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 crhf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rhlf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rhf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
   		      ndf_off=ndf_ejv,ejv=ejv_rhf,afp=s_rhf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
  		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 crhf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rhlf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rhf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
   		      ndf_off=ndf_ejv,ejv=ejv_rhf,afp=s_rhf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
  		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:rhf_destriped.iss', $
              tau_rhf,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_rhf,csp,step_up,step_dn,ndf_ejv,crhf,ejv_rhf, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_rhf=0 & b_rhf=0 & l_rhf=0 & urhf=0 & d_rhf=0
cal_lbl=0 & sky_lbl=0 & sl_weight_rat=0 & d_rhsf=0 & d_rhlf=0

endif
;
; LHS Data
; --------
if ( id eq 'LHS' ) then begin

 tau_lhs=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 clhs = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='lhss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lhs,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
   		      ndf_off=ndf_ejv,ejv=ejv_lhs,afp=s_lhs,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 clhs = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='lhss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lhs,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
   		      ndf_off=ndf_ejv,ejv=ejv_lhs,afp=s_lhs,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:lhs_destriped.iss', $
              tau_lhs,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_lhs,csp,step_up,step_dn,ndf_ejv,clhs,ejv_lhs, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_lhs=0 & b_lhs=0 & l_lhs=0 & ulhs=0 & d_lhs=0

endif
;
; LHF Data
; --------
if ( id eq 'LHF' ) then begin

 tau_lhf=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 clhf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='lhlf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lhf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_lhf,afp=s_lhf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 clhf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='lhlf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lhf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_lhf,afp=s_lhf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:lhf_destriped.iss', $
              tau_lhf,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_lhf,csp,step_up,step_dn,ndf_ejv,clhf,ejv_lhf, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_lhf=0 & b_lhf=0 & l_lhf=0 & ulhf=0 & d_lhf=0
cal_lbl=0 & sky_lbl=0 & sl_weight_rat=0 & d_lhsf=0 & d_lhlf=0

endif
;
; RLS Data
; --------
if ( id eq 'RLS' ) then begin

 tau_rls=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 crls = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rlss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rls,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_rls,afp=s_rls,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 crls = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rlss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rls,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_rls,afp=s_rls,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:rls_destriped.iss', $
              tau_rls,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_rls,csp,step_up,step_dn,ndf_ejv,crls,ejv_rls, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_rls=0 & b_rls=0 & l_rls=0 & urls=0 & d_rls=0

endif
;
; RLF Data
; --------
if ( id eq 'RLF' ) then begin

 tau_rlf=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 crlf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rllf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rlf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_rlf,afp=s_rlf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 crlf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rllf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rlf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_rlf,afp=s_rlf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:rlf_destriped.iss', $
              tau_rlf,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_rlf,csp,step_up,step_dn,ndf_ejv,crlf,ejv_rlf, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_rlf=0 & b_rlf=0 & l_rlf=0 & urlf=0 & d_rlf=0

endif
;
; RSF Data
; --------
if ( id eq 'RSF' ) then begin

 tau_rsf=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 crsf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rlsf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rsf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_rsf,afp=s_rsf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

;
if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 crsf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='rlsf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_rsf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_rsf,afp=s_rsf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:rsf_destriped.iss', $
              tau_rsf,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_rsf,csp,step_up,step_dn,ndf_ejv,crsf,ejv_rsf, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_rsf=0 & b_rsf=0 & l_rsf=0 & ursf=0 & d_rsf=0
cal_lbl=0 & sky_lbl=0 & sl_weight_rat=0 & d_rlsf=0 & d_rlfl=0

endif
;
; LLS Data
; --------
if ( id eq 'LLS' ) then begin

 tau_lls=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 clls = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='llss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lls,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_lls,afp=s_lls,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 clls = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='llss', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lls,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_lls,afp=s_lls,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:lls_destriped.iss', $
              tau_lls,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_lls,csp,step_up,step_dn,ndf_ejv,clls,ejv_lls, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_lls=0 & b_lls=0 & l_lls=0 & ulls=0 & d_lls=0

endif
;
; LLF Data
; --------
if ( id eq 'LLF' ) then begin

 tau_llf=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 cllf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='lllf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_llf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_llf,afp=s_llf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 cllf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='lllf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_llf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_llf,afp=s_llf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:llf_destriped.iss', $
              tau_llf,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_llf,csp,step_up,step_dn,ndf_ejv,cllf,ejv_llf, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_llf=0 & b_llf=0 & l_llf=0 & ullf=0 & d_llf=0

endif
;
; LSF Data
; --------
if ( id eq 'LSF' ) then begin

 tau_lsf=tau_array

; Run DeStriper
; -------------
;
if (vib_flag eq 'Y') then begin
 print,'Calling FDS_DESTRIPER With VIB Keyword'
 clsf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='llsf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lsf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_lsf,afp=s_lsf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,vib=vib,sky_wgts_ds=sky_wgts_ds)
endif

if (vib_flag ne 'Y') then begin
 print,'Calling FDS_DESTRIPER Without VIB Keyword'
 clsf = fds_destriper(f,px,sp,cal_sp,sky_wgts,cal_wgts_ds,xcal, $
                      del_temp=del_temp,index=index,id='llsf', $
                      tm=tm,cal_tm=cal_tm,tau=tau_lsf,solution=solution, $
                      galcut=galcut,step_up=step_up,step_dn=step_dn, $
  		      ndf_off=ndf_ejv,ejv=ejv_lsf,afp=s_lsf,pcvr=pcvr, $
                      sqr_mat=sqr_mat,rect=rect,diag=diag,c_calsp=c_calsp, $
 		      ndf_none=ndf_none,csp=csp,sky_wgts_ds=sky_wgts_ds)
endif

; Save saveset
; ------------
save,filename='csdr$firas_save:lsf_destriped.iss', $
              tau_lsf,galcut,solution,f,del_temp,xcal,sky_wgts_ds, $
              s_lsf,csp,step_up,step_dn,ndf_ejv,clsf,ejv_lsf, $
              sqr_mat,rect,diag,c_calsp,ndf_none, $
              cal_wgts_ds,pcvr,vib,corr_index,index

; Unused Parameters Redefined
; ---------------------------
n_lsf=0 & b_lsf=0 & l_lsf=0 & ulsf=0 & d_lsf=0
cal_lbl=0 & sky_lbl=0 & sl_weight_rat=0 & d_llsf=0 & d_llfl=0

endif
;

; Unused Parameters Redefined
; ---------------------------
nifgs=0 & glon=0 & glat=0 & cal_nifgs=0 & ical=0 & refh=0 & skyh=0
sky_glitch=0 & cal_glitch=0
;
; Set Return Status to NO ERROR .
error=0
;
return
end
