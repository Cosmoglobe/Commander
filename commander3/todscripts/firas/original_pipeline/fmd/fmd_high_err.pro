Pro FMD_HIGH_ERR,error
;
;
;  FMD_HIGH_ERR computes and saves BETA for the high data combined over
;  all frequency bands.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  None
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_n_ERRORS.ISS and
;                        HIGH_n_EJG.ISS, and where HIGH_ERRORS.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_HIGH_ERR,error
;
;
;    HISTORY : Written by Dale Fixsen, Hughes STX, 06-Jun-1997, as ER4.PRO.
;
;              Modified by Ken Jensen, HSTX, 19-Jun-1997 ; renamed FMD_HIGH_ERR
;              and modified for compatibility with the FMD facility.
;
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
 PRINT,'FMD_HIGH_ERR : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_HIGH_ERR,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 RETURN
ENDIF
;

; Logical Translations
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


; Restore HIGH C-Vector
; ---------------------
PRINT,'Restoring HIGH_2_EJG.ISS'
RESTORE,'csdr$firas_out:high_2_ejg.iss'
;
PRINT,'Restoring HIGH_3_EJG.ISS'
RESTORE,'csdr$firas_out:high_3_ejg.iss'
;
PRINT,'Restoring HIGH_4_EJG.ISS'
RESTORE,'csdr$firas_out:high_4_ejg.iss'
;
ch = [c_high_2,c_high_3,c_high_4]
;

; Relative Band Weights
; ---------------------
w = ([TOTAL(1/ch(0:54)^2),TOTAL(1/ch(55:109)^2),TOTAL(1/ch(110:*)^2)])
w = SQRT(w/TOTAL(w)) & w=[0,0,w]
;

; Restore HIGH Band_2 Data
; ------------------------
PRINT,'Restoring HIGH_2_ERRORS.ISS'
RESTORE,'csdr$firas_out:high_2_errors.iss'
;

; Band_2 Dimensionless BETA
; -------------------------
b2 = beta
FOR i=0,16 DO b2(*,i) = w(2) * b2(*,i) * SQRT(diag)
;

; Band_2 Weighted Gamma
; ---------------------
g2 = gamma_high_2 / w(2)
;
nb = REFORM(diag#beta)
e2 = gamma_high_2
FOR i=0,16 DO e2(*,i) = e2(*,i) * nb(i)
;

; Restore HIGH Band_3 Data
; ------------------------
PRINT,'Restoring HIGH_3_ERRORS.ISS'
RESTORE,'csdr$firas_out:high_3_errors.iss'
;

; Band_3 Dimensionless BETA
; -------------------------
b3 = beta
FOR i=0,18 DO b3(*,i) = w(3) * b3(*,i) * SQRT(diag)
;

; Band_3 Weighted Gamma
; ---------------------
g3 = gamma_high_3 / w(3)
nb = REFORM(diag#beta)
e3 = gamma_high_3
FOR i=0,18 DO e3(*,i) = e3(*,i) * nb(i)
;

; Restore HIGH Band_4 Data
; ------------------------
PRINT,'Restoring HIGH_4_ERRORS.ISS'
RESTORE,'csdr$firas_out:high_4_errors.iss'
;

; Band_4 Dimensionless BETA
; -------------------------
b4 = beta
FOR i=0,5 DO b4(*,i) = w(4) * b4(*,i) * SQRT(diag)
;

; Band_4 Weighted Gamma
; ---------------------
g4 = gamma_high_4 / w(4)
nb = REFORM(diag#beta)
e4 = gamma_high_4
FOR i=0,5 DO e4(*,i) = e4(*,i) * nb(i)
;

; Build Weighted Dimensionless BETAs
; ----------------------------------
bf = [[b2],[b3],[b4]] ;The weighted combination=betas*sqrt(N)
;

; SVD for the 10 Most Significant Stripes
; ---------------------------------------
NR_SVD,bf,e,r,v
;

; Combined BETA and GAMMA
; -----------------------
r = r(*,0:9)
v = v(*,0:9)
FOR i=0,9 DO r(*,i) = r(*,i) * e(i)
;
beta = r
FOR i=0,9 DO beta(*,i) = r(*,i) / SQRT(diag)
;
gamma_high = [g2#v(0:16,*),g3#v(17:35,*),g4#v(36:*,*)]
;

f_hi = [f_hi2,f_hi3,f_hi4]
c_high = ch
chanscan = 'HIGH'
;
sname='csdr$firas_out:high_errors.iss'
SAVE,filename=sname,chanscan,cmp_px,beta,diag,gamma_high,f_hi,c_high,$
                    stripe_contrib
                              
PRINT,' '
PRINT,'IDL Save Set "' + outtrans(0) + 'HIGH_ERRORS.ISS" Created.'
PRINT,' '
;

pcvr=0 & rect=0 & stripes_descrip=0 & omega=0 & square=0
d_inv=0 & stripe_conv=0 & destriper_wgt=0 & lmat=0
f_hi2=0 & f_hi3=0 & f_hi4=0
;

; Status = NO Error
; -----------------
error = 0
;

chan_label=0 & stripe_order=0 & stripe_id=0 & dirbe_array=0 & lp_order=0
xtm=0 & cmn_tm=0 & dihed_pow=0 & ref_dihd=0 & dihed_cut=0 & cmn_dihd=0
ref_s0=0 & rms_s0=0 & cmn_bol=0 & step_up=0 & step_dn=0 & n_krnl=0 & n_cmn=0
dirbe_cut=0 & latcut=0 & loncut=0 & del_temp=0 & freq_corr=0 & good_sky_dihd=0
sky_wgts_ds=0 & sky_mask=0 & cvec_mask=0 & cvec_cut=0 & dirbe_mask=0
tmin=0 & tmax=0 & good_cal_dihd=0 & cal_wgts_ds=0 & dihed_cut_min=0
max_frac=0 & n_stripes=0 & freq_band=0
c_high_2=0 & ejg_high_2=0 & ndf_ejg_high_2=0 & chi2_hi2ch=0
c_high_3=0 & ejg_high_3=0 & ndf_ejg_high_3=0 & chi2_hi3ch=0
c_high_4=0 & ejg_high_4=0 & ndf_ejg_high_4=0 & chi2_hi4ch=0
;

RETURN
END
