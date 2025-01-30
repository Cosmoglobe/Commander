;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FEF_6
;
;DESCRIPTION:  
;     A driver to process FEF Slow-Fast Scan Mode Comparison
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;   CHANSCAN  (I) =  Channel/Scan-Mode ID  (e.g. 'RHSS')   = String
;   ERROR     (O) =  Error Status                          = 0 or 1
;
;WARNING:
;   The Logicals CSDR$FIRAS_ERRORS, CSDR$FIRAS_IN, and
;                CSDR$FIRAS_SAVE must be previously defined.
;
;   CSDR$FIRAS_IN must contain the input PEP_ERRORS.xxxx files. To use
;                 production data, directory is FIRCOADD:[ERRORS].
;   CSDR$FIRAS_IN must contain the input xxx.ISS IDL Save Sets. To use
;                 production data, directory is FIRCOADD:[CALMOD].
;   CSDR$FIRAS_SAVE will contain the output xxx_FEF_6.ISS IDL Save Set.
;
;EXAMPLE:
;   $ UPRIDL
;   chanscan='LLSS'
;   FEF_6,chanscan,error
;
;#
;COMMON BLOCKS:
;   None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;   Calls  FEF_SPECTRA.PRO
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;   None
;  
;MODIFICATION HISTORY:
;   Written by Ken Jensen, Hughes STX, Sept 15, 1994
;   Modified by Ken Jensen, Hughes STX, Sept 26, 1994 (SPR 11918)
; 
;-
;___________________________________________________________________
;
Pro FEF_6,chanscan,error
;

; Set Error Status
; ----------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 2 then begin
 print,'FEF_6,chanscan,error'
 return
endif
;

; Restore Save Sets
; -----------------
channel=strupcase(strmid(chanscan,0,2))
;
if (channel eq 'RH') then begin
   id1 = 'RHS'  &  id2 = 'RHF'
;
   restore,'csdr$firas_in:rhs.iss'
   restore,'csdr$firas_in:rhs_destriped.iss'
   n1 = n_rhs  &  s1 = s_rhs  &  c1 = crhs(*,0)
   pcvr1 = pcvr  &  diag1 = diag  &  rect1 = rect
   restore,'csdr$firas_errors:pep_errors.rhss'
   pep1 = gain
;
   restore,'csdr$firas_in:rhf.iss'
   restore,'csdr$firas_in:rhf_destriped.iss'
   n2 = n_rhf  &  s2 = s_rhf  &  c2 = crhf(*,0)
   pcvr2 = pcvr  &  diag2 = diag  &  rect2 = rect
   restore,'csdr$firas_errors:pep_errors.rhlf'
   pep2 = gain
;
endif
;
;
if (channel eq 'RL') then begin
   id1 = 'RLS'  &  id2 = 'RSF'
;
   restore,'csdr$firas_in:rls.iss'
   restore,'csdr$firas_in:rls_destriped.iss'
   n1 = n_rls  &  s1 = s_rls  &  c1 = crls(*,0)
   pcvr1 = pcvr  &  diag1 = diag  &  rect1 = rect
   restore,'csdr$firas_errors:pep_errors.rlss'
   pep1 = gain
;
   restore,'csdr$firas_in:rsf.iss'
   restore,'csdr$firas_in:rsf_destriped.iss'
   n2 = n_rsf  &  s2 = s_rsf  &  c2 = crsf(*,0)
   pcvr2 = pcvr  &  diag2 = diag  &  rect2 = rect
   restore,'csdr$firas_errors:pep_errors.rlsf'
   pep2 = gain
;
endif
;
;
if (channel eq 'LH') then begin
   id1 = 'LHS'  &  id2 = 'LHF'
;
   restore,'csdr$firas_in:lhs.iss'
   restore,'csdr$firas_in:lhs_destriped.iss'
   n1 = n_lhs  &  s1 = s_lhs  &  c1 = clhs(*,0)
   pcvr1 = pcvr  &  diag1 = diag  &  rect1 = rect
   restore,'csdr$firas_errors:pep_errors.lhss'
   pep1 = gain
;
   restore,'csdr$firas_in:lhf.iss'
   restore,'csdr$firas_in:lhf_destriped.iss'
   n2 = n_lhf  &  s2 = s_lhf  &  c2 = clhf(*,0)
   pcvr2 = pcvr  &  diag2 = diag  &  rect2 = rect
   restore,'csdr$firas_errors:pep_errors.lhlf'
   pep2 = gain
;
endif
;
;
if (channel eq 'LL') then begin
   id1 = 'LLS'  &  id2 = 'LSF'
;
   restore,'csdr$firas_in:lls.iss'
   restore,'csdr$firas_in:lls_destriped.iss'
   n1 = n_lls  &  s1 = s_lls  &  c1 = clls(*,0)
   pcvr1 = pcvr  &  diag1 = diag  &  rect1 = rect
   restore,'csdr$firas_errors:pep_errors.llss'
   pep1 = gain
;
   restore,'csdr$firas_in:lsf.iss'
   restore,'csdr$firas_in:lsf_destriped.iss'
   n2 = n_lsf  &  s2 = s_lsf  &  c2 = clsf(*,0)
   pcvr2 = pcvr  &  diag2 = diag  &  rect2 = rect
   restore,'csdr$firas_errors:pep_errors.llsf'
   pep2 = gain
;
endif
;
;
;

; Call FEF_SPECTRA
; ----------------
FEF_SPECTRA,f,s1,s2,n1,n2,c1,c2,diff=diff,d_unc=d_unc,ratio=ratio,r_unc=r_unc,$
      pcvr1=pcvr1,pcvr2=pcvr2,diag1=diag1,diag2=diag2,rect1=rect1,rect2=rect2,$
      pep1=pep1,pep2=pep2
;

; FEF IDL Save Set
; ----------------
tau = [61.,153.]
sname='csdr$firas_save:'+id1+'_fef_6.iss'
save,filename=sname,id1,id2,f,tau,step_up,step_dn,diff,d_unc,ratio,r_unc
;

; Re-Define Unused Restored Parameters
; ------------------------------------
galcut=0 & solution=0 & del_temp=0 & xcal=0 & sky_wgts_ds=0 & csp=0 & c_calsp=0
step_up=0 & step_dn=0 & ndf_ejv=0 & ndf_none=0 & sqr_mat=0 & rect=0 & diag=0
cal_wgts_ds=0 & pcvr=0 & vib=0 & corr_index=0 & index=0 & sl_weight_rat=0
px=0 & tm=0 & sp=0 & nifgs=0 & glon=0 & glat=0 & cal_sp=0 & cal_nifgs=0
cal_tm=0 & ical=0 & refh=0 & skyh=0 & dihd=0 & sky_glitch=0 & cal_glitch=0
sky_wgts=0 & cal_wgts=0
;
ejv_rhs=0 & ejv_rhf=0 & ejv_rls=0 & ejv_rsf=0 & ejv_rlf=0
ejv_lhs=0 & ejv_lhf=0 & ejv_lls=0 & ejv_lsf=0 & ejv_llf=0
tau_rhs=0 & tau_rhf=0 & tau_rls=0 & tau_rsf=0 & tau_rlf=0
tau_lhs=0 & tau_lhf=0 & tau_lls=0 & tau_lsf=0 & tau_llf=0
b_rhs=0 & b_rhf=0 & b_rls=0 & b_rsf=0 & b_rlf=0
l_rhs=0 & l_rhf=0 & l_rls=0 & l_rsf=0 & l_rlf=0
urhs=0 & urhf=0 & urls=0 & ursf=0 & urlf=0
d_rhs=0 & d_rhf=0 & d_rls=0 & d_rsf=0 & d_rlf=0
d_rhsf=0 & d_rhlf=0 & d_rlsf=0 & d_rlfl=0
b_lhs=0 & b_lhf=0 & b_lls=0 & b_lsf=0 & b_llf=0
l_lhs=0 & l_lhf=0 & l_lls=0 & l_lsf=0 & l_llf=0
ulhs=0 & ulhf=0 & ulls=0 & ulsf=0 & ullf=0
d_lhs=0 & d_lhf=0 & d_lls=0 & d_lsf=0 & d_llf=0
d_lhsf=0 & d_lhlf=0 & d_llsf=0 & d_llfl=0
;

; Return Status
; -------------
error = 0
;

;
return
end
