Pro FEF_Write_SPC_X

;  Written by Ken Jensen, Hughes STX, 5 Dec 1994


; Comparison Taus and Mission Periods (34 and 167)
; ------------------------------------------------
tau = [61.,153.]
step_up = ['89328','901391535','901931850']
step_dn = ['9026409','901931850','902081120']


rec_len = 3572
	; fixed length record size in bytes


; HIGH Correction Spectra
; -----------------------

OPENW, 10, 'csdr$firas_out:fef_spc_high.f16_93hybrid', rec_len, /FIXED

descrip1 = 'HISL HIFA Comparison  Std. Destriper'

restore,'csdr$firas_save:fef_8.iss'
chsn1='HISL' & chsn2='HIFA'
d_unc = 0. * diff
r_unc = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn

CLOSE,10
;


; LOWF Correction Spectra
; -----------------------

OPENW, 10, 'csdr$firas_out:fef_spc_lowf.f16_93hybrid', rec_len, /FIXED

descrip1 = 'LOSL LOFA Comparison  Std. Destriper'

restore,'csdr$firas_save:fef_9.iss'
chsn1='LOSL' & chsn2='LOFA'
d_unc = 0. * diff
ratio = 0. * float(diff)
r_unc = ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn

CLOSE,10
;


; LRES Correction Spectra
; -----------------------

OPENW, 10, 'csdr$firas_out:fef_spc_lres.f16_93hybrid', rec_len, /FIXED

descrip1 = 'HIGH LOWF Comparison  Std. Destriper'

restore,'csdr$firas_save:fef_7.iss'
chsn1='HIGH' & chsn2='LOWF'
d_unc = 0. * diff
ratio = 0. * float(diff)
r_unc = ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn

CLOSE,10
;


; Re-Define Restored Parameters
; -----------------------------
idd=0 & f=0

;
RETURN
END
