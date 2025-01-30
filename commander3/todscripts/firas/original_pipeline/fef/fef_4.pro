;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FEF_4
;
;DESCRIPTION:  
;     A driver to process FEF Alternative_4 ( 3 Mission Periods, 1 UFO )
;     There are three sets of input TAU, determined by maximizing
;     Delta_Chi_Squared over elected frequency ranges.
;     The first TAU = 1790 days, determined by maximizing LLSS
;     Delta_Chi_Squared over band [2.3-4.5]icm, is applied to
;     Low-Channel combinations.
;     The second TAU = 40 days, determined by maximizing LLSS
;     Delta_Chi_Squared over band [7.9-10.2]icm, is applied to
;     Low-Channel combinations.
;     The third TAU = 20 days, determined by maximizing RHSS
;     Delta_Chi_Squared over band [28.9-33.4]icm, is applied to
;     High-Channel combinations.
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;   CHANSCAN  (I) =  Channel/Scan-Mode ID  (e.g. 'RHSS')   = String
;   ERROR     (O) =  Error Status                          = 0 or 1
;
;WARNING:
;   The Logicals CSDR$FIRAS_ERRORS , CSDR$FIRAS_IN, and
;                CSDR$FIRAS_SAVE must be previously defined.
;
;   CSDR$FIRAS_ERRORS must contain the input PEP_ERRORS.xxxx files. To use
;                 production data, directory is FIRCOADD:[ERRORS].
;   CSDR$FIRAS_IN must contain the input xxx.ISS IDL Save Set. To use
;                 production data, directory is FIRCOADD:[CALMOD].
;   CSDR$FIRAS_SAVE will contain the output xxx_FEF_y.ISS IDL Save Sets,
;   ('y' identifies the destriping model, e.g. '1790_3MP' )
;
;EXAMPLE:
;   $ UPRIDL
;   chanscan='LLSS'
;   FEF_4,chanscan,error
;
;#
;COMMON BLOCKS:
;   None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;   Calls  FEF_DSTRIPE_SKY.PRO
;          FEF_SPECTRA.PRO
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
Pro FEF_4,chanscan,error
;

; Set Error Status
; ----------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 2 then begin
 print,'FEF_4,chanscan,error'
 return
endif
;

chanscan=strupcase(chanscan)
if(chanscan eq 'RHSF')then chanscan='RHLF'
if(chanscan eq 'RLFL')then chanscan='RLSF'
if(chanscan eq 'LHSF')then chanscan='LHLF'
if(chanscan eq 'LLFL')then chanscan='LLSF'
;

; Acquire PEP Gain Errors
; -----------------------
restore,'csdr$firas_errors:pep_errors.'+chanscan
pep1 = gain  &  pep2 = gain
;

; Mission Periods
; ---------------
step_up1 = ['89328','901391535','901931850']
step_dn1 = ['9026409','901931850','902081120']
;
if((chanscan eq 'RLLF')or(chanscan eq 'LLLF'))then begin
 step_up1 = ['89328','901391535']
 step_dn1 = ['9026409','901931850']
endif
;

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

; Restore Comparison FDS Skymap
; -----------------------------
restore,'csdr$firas_in:'+idd+'_destriped.iss'
;

; Comparison Skymap Spectra (Eplees)
; ----------------------------------
if(idd eq 'RHS')then s2=s_rhs
if(idd eq 'RHF')then s2=s_rhf
if(idd eq 'RLS')then s2=s_rls
if(idd eq 'RSF')then s2=s_rsf
if(idd eq 'RLF')then s2=s_rlf
if(idd eq 'LHS')then s2=s_lhs
if(idd eq 'LHF')then s2=s_lhf
if(idd eq 'LLS')then s2=s_lls
if(idd eq 'LSF')then s2=s_lsf
if(idd eq 'LLF')then s2=s_llf
;

; Comparison C-Vector (Eplees)
; ----------------------------
if(idd eq 'RHS')then c2=crhs(*,0)
if(idd eq 'RHF')then c2=crhf(*,0)
if(idd eq 'RLS')then c2=crls(*,0)
if(idd eq 'RSF')then c2=crsf(*,0)
if(idd eq 'RLF')then c2=crlf(*,0)
if(idd eq 'LHS')then c2=clhs(*,0)
if(idd eq 'LHF')then c2=clhf(*,0)
if(idd eq 'LLS')then c2=clls(*,0)
if(idd eq 'LSF')then c2=clsf(*,0)
if(idd eq 'LLF')then c2=cllf(*,0)
;

; Comparison Error Elements
; -------------------------
pcvr2=pcvr & diag2=diag & rect2=rect
;

; First TAU
; ---------
tau1 = [1790.]
;

; Low Channels Only
; -----------------
if(strmid(chanscan,1,1) eq 'L')then begin  ; Low Channels
;

; Call FEF_DSTRIPE_SKY
; --------------------
 FEF_DSTRIPE_SKY,chanscan,tau1,step_up1,step_dn1,s1,n1,c1,pcvr1,diag1,rect1,error
;
; S1 = Destriped Skymap Spectra (Eplees)
; N1 = Skymap Weights
; C1 = C-Vector (Eplees)
;

; Return IF(Error)
; ----------------
 if (error ne 0) then return

; Call FEF_SPECTRA
; ----------------
 FEF_SPECTRA,f,s1,s2,n1,n1,c1,c2,diff=diff,d_unc=d_unc,ratio=ratio,r_unc=r_unc,$
       pcvr1=pcvr1,pcvr2=pcvr2,diag1=diag1,diag2=diag2,rect1=rect1,rect2=rect2,$
       pep1=pep1,pep2=pep2
;

; FEF IDL Save Set
; ----------------
 sname='csdr$firas_save:'+idd+'_fef_4a.iss'
 save,filename=sname,idd,tau1,step_up1,step_dn1,f,diff,d_unc,ratio,r_unc
;

endif
;


; Second TAU
; ----------
tau1 = [40.]
;

; Low Channels Only
; -----------------
if(strmid(chanscan,1,1) eq 'L')then begin  ; Low Channels
;

; Call FEF_DSTRIPE_SKY
; --------------------
 FEF_DSTRIPE_SKY,chanscan,tau1,step_up1,step_dn1,s1,n1,c1,pcvr1,diag1,rect1,error
;
; S1 = Destriped Skymap Spectra (Eplees)
; N1 = Skymap Weights
; C1 = C-Vector (Eplees)
;

; Return IF(Error)
; ----------------
 if (error ne 0) then return

; Call FEF_SPECTRA
; ----------------
 FEF_SPECTRA,f,s1,s2,n1,n1,c1,c2,diff=diff,d_unc=d_unc,ratio=ratio,r_unc=r_unc,$
       pcvr1=pcvr1,pcvr2=pcvr2,diag1=diag1,diag2=diag2,rect1=rect1,rect2=rect2,$
       pep1=pep1,pep2=pep2
;

; FEF IDL Save Set
; ----------------
 sname='csdr$firas_save:'+idd+'_fef_4b.iss'
 save,filename=sname,idd,tau1,step_up1,step_dn1,f,diff,d_unc,ratio,r_unc
;

endif
;



; Third TAU
; ----------
tau1 = [20.]
;

; High Channels Only
; ------------------
if(strmid(chanscan,1,1) eq 'L')then return
;

; Call FEF_DSTRIPE_SKY
; --------------------
FEF_DSTRIPE_SKY,chanscan,tau1,step_up1,step_dn1,s1,n1,c1,pcvr1,diag1,rect1,error
;
; S1 = Destriped Skymap Spectra (Eplees)
; N1 = Skymap Weights
; C1 = C-Vector (Eplees)
;

; Return IF(Error)
; ----------------
if (error ne 0) then return

; Call FEF_SPECTRA
; ----------------
FEF_SPECTRA,f,s1,s2,n1,n1,c1,c2,diff=diff,d_unc=d_unc,ratio=ratio,r_unc=r_unc,$
      pcvr1=pcvr1,pcvr2=pcvr2,diag1=diag1,diag2=diag2,rect1=rect1,rect2=rect2,$
      pep1=pep1,pep2=pep2
;

; FEF IDL Save Set
; ----------------
sname='csdr$firas_save:'+idd+'_fef_4c.iss'
save,filename=sname,idd,tau1,step_up1,step_dn1,f,diff,d_unc,ratio,r_unc
;

; Re-Define Unused Restored Parameters
; ------------------------------------
galcut=0 & solution=0 & del_temp=0 & xcal=0 & sky_wgts_ds=0 & csp=0 & c_calsp=0
step_up=0 & step_dn=0 & ndf_ejv=0 & ndf_none=0 & sqr_mat=0 & rect=0 & diag=0
cal_wgts_ds=0 & pcvr=0 & vib=0 & corr_index=0 & index=0
;
ejv_rhs=0 & ejv_rhf=0 & ejv_rls=0 & ejv_rsf=0 & ejv_rlf=0
ejv_lhs=0 & ejv_lhf=0 & ejv_lls=0 & ejv_lsf=0 & ejv_llf=0
tau_rhs=0 & tau_rhf=0 & tau_rls=0 & tau_rsf=0 & tau_rlf=0
tau_lhs=0 & tau_lhf=0 & tau_lls=0 & tau_lsf=0 & tau_llf=0
;

;
return
end
