;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FEF_9
;
;DESCRIPTION:  
;     A driver to process FEF LOSL vs LOFA Comparison
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
;   CSDR$FIRAS_IN must contain the input LOSL.ISS and LOFA.ISS IDL Save Sets.
;                 These must be created by FTT prior to executing this code.
;   CSDR$FIRAS_SAVE will contain the output FEF_9.ISS IDL Save Set.
;
;EXAMPLE:
;   $ UPRIDL
;   FEF_9,error
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
;   Written by Ken Jensen, Hughes STX, Nov 28, 1994
; 
;-
;___________________________________________________________________
;
Pro FEF_9,error
;

; Set Error Status
; ----------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FEF_9,error'
 return
endif
;

; Restore Save Sets
; -----------------
restore,'csdr$firas_in:losl.iss'
n1 = n_losl  &  s1 = s_losl  &  c1 = c_losl
;
restore,'csdr$firas_in:lofa.iss'
n2 = n_lofa  &  s2 = s_lofa  &  c2 = c_lofa
;

; Call FEF_SPECTRA_X
; ------------------
FEF_SPECTRA_X,s1,s2,n1,n2,c1,c2,diff=diff
;

; FEF IDL Save Set
; ----------------
id1='LOSL' & id2='LOFA'
sname='csdr$firas_save:fef_9.iss'
save,filename=sname,id1,id2,diff
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
b_losl=0 & d_losl=0 & b_lofa=0 & d_lofa=0
;

; Return Status
; -------------
error = 0
;

;
return
end
