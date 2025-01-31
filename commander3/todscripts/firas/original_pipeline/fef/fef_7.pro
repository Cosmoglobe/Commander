;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FEF_7
;
;DESCRIPTION:  
;     A driver to process FEF HIGH-LOWF Comparison
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;   ERROR     (O) =  Error Status                          = 0 or 1
;
;WARNING:
;   The Logicals CSDR$FIRAS_IN and CSDR$FIRAS_SAVE must be previously defined.
;
;   CSDR$FIRAS_IN must contain the input HIGH.ISS and LOWF.ISS IDL Save Sets.
;                 These must be created by FTT prior to executing this code.
;   CSDR$FIRAS_SAVE will contain the output FEF_7.ISS IDL Save Set.
;
;EXAMPLE:
;   $ UPRIDL
;   FEF_7,error
;
;#
;COMMON BLOCKS:
;   None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;   Calls  FEF_SPECTRA_X.PRO
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;   None
;  
;MODIFICATION HISTORY:
;   Written by Ken Jensen, Hughes STX, Nov 15, 1994
; 
;-
;___________________________________________________________________
;
Pro FEF_7,error
;

; Set Error Status
; ----------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FEF_7,error'
 return
endif
;

; Restore Save Sets
; -----------------
idx=indgen(34)
restore,'csdr$firas_in:high.iss'
n1 = n_high  &  s1 = s_high(*,*,idx(*))  &  c1 = c_high(idx(*))
;
restore,'csdr$firas_in:lowf.iss'
n2 = n_lowf  &  s2 = s_lowf  &  c2 = c_lowf
;

; Call FEF_SPECTRA_X
; ------------------
FEF_SPECTRA_X,s1,s2,n1,n2,c1,c2,diff=diff
;

; FEF IDL Save Set
; ----------------
id1='HIGH' & id2='LOWF'
sname='csdr$firas_save:fef_7.iss'
save,filename=sname,id1,id2,diff
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
b_high=0 & d_high=0 & b_lowf=0 & d_lowf=0
;

; Return Status
; -------------
error = 0
;

;
return
end
