;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FDS_PROCESS_SKY 
;
;DESCRIPTION:  
;     This command file drives the 'fds_read_data' procedure which 
;     reads in data from the FCF_SKY files and stores it in an IDL saveset 
;     for later processing by 'fds_destriper'.  
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     CHANSCAN (I)   Channel/MTM Scan Mode  String  (e.g. 'LLSS')
;     ERROR    (O)   Error Status
;
;WARNINGS:
;     The Logicals CSDR$FIRAS_INSKY, CSDR$FIRAS_INCAL, CSDR$FIRAS_REF,
;     and CSDR$FIRAS_SAVE must be previously defined.
;     CSDR$FIRAS_INSKY should contain the desired FCF_SKY records.
;     CSDR$FIRAS_INCAL should contain the desired FCF_CAL records.
;     CSDR$FIRAS_REF should contain the desired FDS_ARCV_TIMES and
;                    FEX_GLTCHCOR reference data sets.
;     CSDR$FIRAS_SAVE must contain the input XXX_VAR.ISS IDL Save Set.
;     CSDR$FIRAS_SAVE will contain the output XXX.ISS IDL Save Set.
;
;EXAMPLE:
;      $ UIDL
;      chanscan = 'RHSS'
;      fds_process_sky,chanscan,error
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;     Calls FDS_READ_DATA.PRO
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;    Written by Joel Gales, ARC, April 1993
;    Revised by Shirley M. Read, Hughes STX, April 26,1993 to run in
;         L1 production.
;    Revised by Ken Jensen, Hughes STX, 30-APR-1993, to add ARC_FILE,
;         WEIGHT_COR, and SC_LEN to the FDS_READ_DATA calling statement.
;    Revised by Ken Jensen, Hughes STX, 30-APR-1993, to add code defining
;         ARC_FILE
;    Modified by Shirley M. Read, Hughes STX, May 12, 1993 to process
;        combinations of channel and MTM scan modes. The short fast and
;        and long fast data will be combined for the high channels. The
;        short slow data for the high channels will be processed alone
;        since there is no long slow data taken during the mission. The
;        current updates are for the high channels only.  The output 
;        variables for the IDL save sets must be given unique names so 
;        that further analysis can be done after destriping. The spectra
;        processed in this updated version are now complex.
;    Modified by Ken Jensen, Hughes STX, May 17, 1993 to correct the
;        statements which extract substrings from SMODE.
;    Modified by Ken Jensen, Hughes STX, May 20, 1993 to restore save sets
;        from csdr$firas_save rather than local directory.
;    Modified by Ken Jensen, Hughes STX, July 21, 1993 to handle the
;        low frequency channels. Eight separate data sets are now
;        defined ( RHS, RHF, LHS, LHF, RLS, RLF, LLS, LLF ) and processed
;        accordingly.
;    Modified by Ken Jensen, Hughes STX, February 16, 1994 to pass the
;        channel/scan-mode string into FDS_READ_DATA explicitly rather than
;        as CHANSCAN. This prevents CHANSCAN from being passed back to the
;        main program as a string array.
;    Modified by Ken Jensen, Hughes STX, March 30, 1994. Additional
;        documentation added.
;    Modified by Ken Jensen, Hughes STX, May 12, 1994. GALCUT changed
;        from 7 to 10.
;    Modified by Ken Jensen, Hughes STX, June 2, 1994. Dihedral temperature
;        DIHD included in call to FDS_READ_DATA, and saved in save set.
;    Modified by Ken Jensen, Hughes STX, June 20, 1994. VR_L and VR_S
;        re-defined.
;    Modified by Ken Jensen, Hughes STX, Aug 11, 1994. NL NS, and NTOT
;        re-defined.
;
;-
;___________________________________________________________________________
;
Pro FDS_Process_Sky,chanscan,error
;
;  Set Return Error Status
;  -----------------------
error=1

;  Procedure Invoked Correctly ?
;  -----------------------------
if N_Params() ne 2 then begin
 print,'FDS_Process_Sky,chanscan,error'
 return
endif
;
;  Set GALCUT ( = Galactic Plane Exclusion )
;  -----------------------------------------
galcut = 10
;
; Process the Input CHANSCAN
; --------------------------
chanscan=strupcase(chanscan)
;

if ( chanscan eq 'RHSS' ) then begin

; RHS Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:rhs_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,'fds_arcv_times_ss','rhss',px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_rhs,b_rhs,l_rhs, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=urhs,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:rhs.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_rhs,b_rhs,l_rhs,urhs,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,d_rhs

endif
;

if ( chanscan eq 'RHLF' ) then begin

; RHF Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:rhf_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,['fds_arcv_times_lf','fds_arcv_times_sf'], $ 
              ['rhlf','rhsf'],px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_rhf,b_rhf,l_rhf, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=urhf,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl,weight_cor=sl_weight_rat

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:rhf.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_rhf,b_rhf,l_rhf,urhf,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,cal_lbl,sky_lbl,d_rhsf,d_rhlf,d_rhf,sl_weight_rat

endif
;

if ( chanscan eq 'LHSS' ) then begin

; LHS Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:lhs_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,'fds_arcv_times_ss','lhss',px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_lhs,b_lhs,l_lhs, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=ulhs,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:lhs.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_lhs,b_lhs,l_lhs,ulhs,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,d_lhs

endif
;

if ( chanscan eq 'LHLF' ) then begin

; LHF Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:lhf_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,['fds_arcv_times_lf','fds_arcv_times_sf'], $ 
              ['lhlf','lhsf'],px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_lhf,b_lhf,l_lhf, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=ulhf,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl,weight_cor=sl_weight_rat

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:lhf.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_lhf,b_lhf,l_lhf,ulhf,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,cal_lbl,sky_lbl,d_lhsf,d_lhlf,d_lhf,sl_weight_rat

endif
;

if ( chanscan eq 'RLSS' ) then begin

; RLS Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:rls_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,'fds_arcv_times_ss','rlss',px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_rls,b_rls,l_rls, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=urls,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:rls.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_rls,b_rls,l_rls,urls,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,d_rls

endif
;

if ( chanscan eq 'RLLF' ) then begin

; RLF Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:rlf_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,'fds_arcv_times_lf','rllf',px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_rlf,b_rlf,l_rlf, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=urlf,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:rlf.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_rlf,b_rlf,l_rlf,urlf,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,d_rlf

endif

;
if ( chanscan eq 'LLSS' ) then begin

; LLS Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:lls_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,'fds_arcv_times_ss','llss',px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_lls,b_lls,l_lls, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=ulls,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:lls.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_lls,b_lls,l_lls,ulls,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,d_lls

endif
;

if ( chanscan eq 'LLLF' ) then begin

; LLF Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:llf_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,'fds_arcv_times_lf','lllf',px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_llf,b_llf,l_llf, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=ullf,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:llf.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_llf,b_llf,l_llf,ullf,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,d_llf

endif

if ( chanscan eq 'RLSF' ) then begin

; RSF Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:rsf_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,['fds_arcv_times_lf','fds_arcv_times_sf'], $ 
              ['rlfl','rlsf'],px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_rsf,b_rsf,l_rsf, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=ursf,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl,weight_cor=sl_weight_rat

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:rsf.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_rsf,b_rsf,l_rsf,ursf,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,cal_lbl,sky_lbl,d_rlsf,d_rlfl,d_rsf,sl_weight_rat

endif

if ( chanscan eq 'LLSF' ) then begin

; LSF Processing
; --------------

; Restore variance data
; ---------------------
 restore,'csdr$firas_save:lsf_var.iss'

; Read SKY and CAL spectral records
; ---------------------------------
fds_read_data,['fds_arcv_times_lf','fds_arcv_times_sf'], $ 
              ['llfl','llsf'],px,tm,sp, $
	      nifgs,glon,glat,solution,f,n_lsf,b_lsf,l_lsf, $
	      cal_sp,cal_nifgs,cal_tm,xcal,ical,refh,skyh, $
              dihd,uxxx=ulsf,galcut=galcut,sky_glitch=sky_glitch, $
	      cal_glitch=cal_glitch,sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
	      sky_lbl=sky_lbl,cal_lbl=cal_lbl,weight_cor=sl_weight_rat

; Save data to saveset
; --------------------
save,filename='csdr$firas_save:lsf.iss',px,tm,sp,nifgs,glon,glat, $
           solution,f,n_lsf,b_lsf,l_lsf,ulsf,galcut,cal_sp,cal_nifgs, $
           cal_tm,xcal,ical,refh,skyh,dihd,sky_glitch,cal_glitch,sky_wgts, $
	   cal_wgts,cal_lbl,sky_lbl,d_llsf,d_llfl,d_lsf,sl_weight_rat

endif

; Unused Arrays Redefined
; -------------------------
vr=0 & vr_l=0 & vr_s=0
nl=0 & ns=0 & ntot=0
;
; Set Return Status to NO ERROR
; -----------------------------
error=0

;
return
end
