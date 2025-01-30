Pro FMD_HIGH_SKYMAP,error
;
;  FMD_HIGH_SKYMAP concatenates C-Vector, skymap spectra, and zodi spectra
;  for the three bands of combined HIGH destriped data, and stores them in
;  an IDL save set .
;
;
;  ARGUMENTS (I/O)    :  ERROR (O)   :  Return Error Status
;
;
;  PROGRAMS Called    :  None
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_n_SKYMAP.ISS
;                        and where HIGH_SKYMAP.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_HIGH_SKYMAP,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 23-May-1997
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;; BEGIN
;
; Initialize Return Error
; -----------------------
error = 1
;

; Procedure Invoked Correctly ?  If not, signal and RETURN
; --------------------------------------------------------
IF N_Params() ne 1 THEN BEGIN 
 PRINT,'FMD_HIGH_SKYMAP : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_HIGH_SKYMAP,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 RETURN
ENDIF
;

; Logical translations
; --------------------
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
outtrans = STRUPCASE(outtrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_OUT    == ' + outtrans
PRINT,' '
;

; Skymap Spectra, Zodi Spectra, and C-Vector
; ------------------------------------------
RESTORE,'csdr$firas_out:high_2_skymap.iss'
freq = TEMPORARY(freq_high_2)
cvec = TEMPORARY(c_high_2)
rect_2 = TEMPORARY(rect)
pcvr_2 = TEMPORARY(pcvr)
str_descrip_2 = stripes_descrip
;
RESTORE,'csdr$firas_out:high_3_skymap.iss'
freq = [freq,TEMPORARY(freq_high_3)]
cvec = [cvec,TEMPORARY(c_high_3)]
rect_3 = TEMPORARY(rect)
pcvr_3 = TEMPORARY(pcvr)
str_descrip_3 = stripes_descrip
;
RESTORE,'csdr$firas_out:high_4_skymap.iss'
freq_high = [freq,TEMPORARY(freq_high_4)]
c_high = [cvec,TEMPORARY(c_high_4)]
rect_4 = TEMPORARY(rect)
pcvr_4 = TEMPORARY(pcvr)
str_descrip_4 = stripes_descrip
;

chi2_high = chi2_high_2 + chi2_high_3 + chi2_high_4
;

s_high = DBLARR(96,64,170)
FOR i=0,54 DO s_high(*,*,i) = s_high_2(*,*,i)
FOR i=0,54 DO s_high(*,*,i+55) = s_high_3(*,*,i)
FOR i=0,59 DO s_high(*,*,i+110) = s_high_4(*,*,i)
;

z_high = FLTARR(96,64,170)
FOR i=0,54 DO z_high(*,*,i) = z_high_2(*,*,i)
FOR i=0,54 DO z_high(*,*,i+55) = z_high_3(*,*,i)
FOR i=0,59 DO z_high(*,*,i+110) = z_high_4(*,*,i)
;

; Make HIGH_SKYMAP.ISS Save Set
; -----------------------------
chanscan = 'HIGH'
sname = 'csdr$firas_out:high_skymap.iss'
SAVE,filename=sname,chanscan,str_descrip_2,str_descrip_3,str_descrip_4,$
                    cmp_px,rect_2,rect_3,rect_4,pcvr_2,pcvr_3,pcvr_4,$
                    freq_high,s_high,z_high,c_high,n_high,b_high,l_high,$
                    nc_high,chi2_high,stripe_contrib
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_SKYMAP.ISS" Created.'
PRINT,' '
;

beta=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
