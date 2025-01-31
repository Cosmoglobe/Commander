;_______________________________________________________________________
;
;+NAME/BRIEFDESCRIPTION OF ROUTINE:
;     FDS_PROCESS_VAR computes the sky sigmas and makes an IDL save set.
;
;DESCRIPTION:  
;    This command file drives the 'fds_variance' procedure.  It reads 
;    in the coadd variances from the FCF_SKY files, computes the sigmas, 
;    and stores then in the 'xxy_var.iss' IDL saveset, where xx = channel
;    and y = scan speed.  It also generates the 'd_xxx' sigma vector.
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     CHANSCAN (I)  Channel/MTM Scan Mode String (e.g. 'LLSS')
;     ERROR    (O)  Error Status
;
;WARNINGS:
;     The Logicals CSDR$FIRAS_INSKY, CSDR$FIRAS_REF, and CSDR$FIRAS_SAVE
;     must be previously defined.
;     CSDR$FIRAS_INSKY should contain the desired FCF_SKY records.
;     CSDR$FIRAS_REF should contain the desired FDS_ARCV_TIMES file.
;     CSDR$FIRAS_SAVE will contain the output XXX_VAR.ISS IDL Save Set.
;
;EXAMPLE:
;     $ UIDL
;     chanscan = 'RHSS'
;     fds_process_var,chanscan,error
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;     Calls FDS_VARIANCE.PRO
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;    Written by Joel Gales, ARC, April 1993
;    Revised by Shirley M. Read, Hughes STX, April 26, 1993 to run
;        under the same setup as the main FDS Destriper procedures
;        which will be run in L1 production.
;    Modified by Shirley M. Read, Hughes STX, May 12, 1993 to process
;        combinations of channel and MTM scan modes. The short fast and
;        and long fast data will be combined for the high channels. The
;        short slow data for the high channels will be processed alone
;        since there is no long slow data taken during the mission. The
;        current updates are for the high channels only. Thus FDS_Variance
;        will be invoked twice for fast data and the resulting variances
;        and sigmas combined for short and long stroke data. In addition
;        the output variables for the IDL save sets must be given unique
;        names so that further analysis can be done after destriping.
;        The spectra processed in this updated version are now complex.
;    Modified by Ken Jensen, Hughes STX, May 17, 1993 to correct the
;        statements which extract substrings from SMODE.
;    Modified by Ken Jensen, Hughes STX, July 21, 1993 to handle the
;        low frequency channels. Eight separate data sets are now
;        defined ( RHS, RHF, LHS, LHF, RLS, RLF, LLS, LLF ) and processed
;        accordingly.
;    Modified by Ken Jensen, Hughes STX, March 1, 1994 to handle the
;        RSF and LSF combinations.  Ten channel/scan-mode combinations
;        ( RHS, RHF, LHS, LHF, RLS, RLF, LLS, LLF, LSF, RSF ) can be
;        processed separately.
;    Modified by Ken Jensen, Hughes STX, March 30, 1994. Additional
;        documentation.
;    Modified by Ken Jensen, Hughes STX, July 26, 1994. GALCUT added to
;        arguments in call to FDS_VARIANCE (SER 11859)
;    Modified by Ken Jensen, Hughes STX, Aug 8, 1994.  SF/LF and SF/FL
;        combined D-Vectors computed with new relative weights (SER 11860)
;
;-
;______________________________________________________________________
;
Pro FDS_Process_Var,chanscan,error
;
; Initialize Return Error Status
; ------------------------------
error=1

; Procedure Invoked Correctly ?
; -----------------------------
if N_params() ne 2 then begin
 print,'FDS_Process_Var,chanscan,error'
 return
endif
;
;  Define the FDS 3-letter ID (RHSS -> RHS etc.)
chanscan=strupcase(chanscan)
if ( (chanscan ne 'RLSF') and (chanscan ne 'LLSF') ) then begin
 id = strupcase(strmid(chanscan,0,2)+strmid(chanscan,3,1))
endif
if (chanscan eq 'RLSF') then id = 'RSF'
if (chanscan eq 'LLSF') then id = 'LSF'
idd=strmid(id,1,2)

; Define GALCUT
; -------------
galcut = 10

;
; Define FDS_ARCV_TIMES filename
; ------------------------------
;
arcv_name='fds_arcv_times_'+strmid(chanscan,2,2)
;
;  Call FDS_VARIANCE.PRO . This procedure reads in the
;  coadd variances from the FCF_SKY files and computes the D-Vector. 
;  The variances and D-Vector are stored in the 'XXX_VAR.ISS' IDL saveset.
;
; For RHS, LHS, RLF, LLF, RLS, LLS, FDS_VARIANCE is called once.
;
if ( ( idd ne 'HF' ) and ( idd ne 'SF' ) ) then begin
;

; Call FDS_VARIANCE
; -----------------
 fds_variance,arcv_name,'csdr$firas_insky:',chanscan,vr,d_xxx,n,galcut
 ntot = n
;
; RHS
; ---
 if(id eq 'RHS')then begin
  d_rhs=d_xxx
; Save data in saveset
  save,filename='csdr$firas_save:rhs_var.iss',vr,d_rhs,ntot,galcut
 endif
;
; LHS
; ---
 if(id eq 'LHS')then begin
  d_lhs=d_xxx
; Save data in saveset
  save,filename='csdr$firas_save:lhs_var.iss',vr,d_lhs,ntot,galcut
 endif
;
; RLS
; ---
 if(id eq 'RLS')then begin
  d_rls=d_xxx
; Save data in saveset
  save,filename='csdr$firas_save:rls_var.iss',vr,d_rls,ntot,galcut
 endif
;
; RLF
; ---
 if(id eq 'RLF')then begin
  d_rlf=d_xxx
; Save data in saveset
  save,filename='csdr$firas_save:rlf_var.iss',vr,d_rlf,ntot,galcut
 endif
;
; LLS
; ---
 if(id eq 'LLS')then begin
  d_lls=d_xxx
; Save data in saveset
  save,filename='csdr$firas_save:lls_var.iss',vr,d_lls,ntot,galcut
 endif
;
; LLF
; ---
 if(id eq 'LLF')then begin
  d_llf=d_xxx
; Save data in saveset
  save,filename='csdr$firas_save:llf_var.iss',vr,d_llf,ntot,galcut
 endif
;
endif

;
; For RHF, LHF, RSF, LSF, FDS_VARIANCE is called twice.
; The weighting ratio between the two scan modes is computed, as is
; a weighted combined D-Vector.
;
; RHF
; ---
if (id eq 'RHF') then begin
;
;  Call FDS_VARIANCE separately for LF and SF
 fds_variance,'fds_arcv_times_lf','csdr$firas_insky:','rhlf',vr_l,d_rhlf,nl,galcut
 fds_variance,'fds_arcv_times_sf','csdr$firas_insky:','rhsf',vr_s,d_rhsf,ns,galcut
 ntot = nl + ns
;
;  Compute SF/LF Weight Ratio
 sl_weight_rat = TOTAL(1/(d_rhsf^2))/TOTAL(1/(d_rhlf^2))
;
;  Compute combined coadd variance
 nl_adj = nl / sl_weight_rat
 ns_adj = ns * sl_weight_rat
 d_rhf = SQRT((nl_adj*(d_rhlf^2.) + ns_adj*(d_rhsf^2.))/(nl_adj+ns_adj))

;  Save Data in Saveset.
 save,filename='csdr$firas_save:rhf_var.iss',vr_l,vr_s, $
 	d_rhsf,d_rhlf,d_rhf,sl_weight_rat,nl,ns,ntot,galcut
;
endif
;
;
; LHF
; ---
if ( id eq 'LHF' ) then begin
;
;  Call FDS_VARIANCE separately for LF and SF
 fds_variance,'fds_arcv_times_lf','csdr$firas_insky:','lhlf',vr_l,d_lhlf,nl,galcut
 fds_variance,'fds_arcv_times_sf','csdr$firas_insky:','lhsf',vr_s,d_lhsf,ns,galcut
 ntot = nl + ns
;
;  Compute SF/LF Weight Ratio
 sl_weight_rat = TOTAL(1/(d_lhsf^2))/TOTAL(1/(d_lhlf^2))
;
;  Compute combined coadd variance
 nl_adj = nl / sl_weight_rat
 ns_adj = ns * sl_weight_rat
 d_lhf = SQRT((nl_adj*(d_lhlf^2.) + ns_adj*(d_lhsf^2.))/(nl_adj+ns_adj))

;  Save Data in Saveset.
 save,filename='csdr$firas_save:lhf_var.iss',vr_l,vr_s, $
 	d_lhsf,d_lhlf,d_lhf,sl_weight_rat,nl,ns,ntot,galcut
;
endif
;
; RSF
; ---
if ( id eq 'RSF' ) then begin
;
;  Call FDS_VARIANCE separately for LF and SF
 fds_variance,'fds_arcv_times_lf','csdr$firas_insky:','rlfl',vr_l,d_rlfl,nl,galcut
 fds_variance,'fds_arcv_times_sf','csdr$firas_insky:','rlsf',vr_s,d_rlsf,ns,galcut
 ntot = nl + ns

;  Compute SF/LF Weight Ratio
 sl_weight_rat = TOTAL(1/(d_rlsf^2))/TOTAL(1/(d_rlfl^2))

;  Compute combined coadd variance
 nl_adj = nl / sl_weight_rat
 ns_adj = ns * sl_weight_rat
 d_rsf = SQRT((nl_adj*(d_rlfl^2.) + ns_adj*(d_rlsf^2.))/(nl_adj+ns_adj))

;  Save data in savesets
 save,filename='csdr$firas_save:rsf_var.iss',vr_l,vr_s, $
       	        d_rlsf,d_rlfl,d_rsf,sl_weight_rat,nl,ns,ntot,galcut

endif
;
; LSF
; ---
if ( id eq 'LSF' ) then begin
;
;  Call FDS_VARIANCE separately for LF and SF
 fds_variance,'fds_arcv_times_lf','csdr$firas_insky:','llfl',vr_l,d_llfl,nl,galcut
 fds_variance,'fds_arcv_times_sf','csdr$firas_insky:','llsf',vr_s,d_llsf,ns,galcut
 ntot = nl + ns

;  Compute SF/LF Weight Ratio
 sl_weight_rat = TOTAL(1/(d_llsf^2))/TOTAL(1/(d_llfl^2))

;  Compute combined coadd variance
 nl_adj = nl / sl_weight_rat
 ns_adj = ns * sl_weight_rat
 d_lsf = SQRT((nl_adj*(d_llfl^2.) + ns_adj*(d_llsf^2.))/(nl_adj+ns_adj))

;  Save data in savesets
 save,filename='csdr$firas_save:lsf_var.iss',vr_l,vr_s, $
       	        d_llsf,d_llfl,d_lsf,sl_weight_rat,nl,ns,ntot,galcut

endif
;
;
; Set Return Status to NO ERROR .
error=0
;
return
end
