Pro FEF_Write_SPC

;  Written by Ken Jensen, Hughes STX, 19 Sept 1994
;  Modifications by Ken Jensen, Hughes STX, 23 Sept 1994,
;    - Change Filename Extension from DAT to F16_93HYBRID
;    - Reverse order of 148 FEF_1 and FEF_2 so LLLF is first.
;    - Change the arguments to FEF_WRITE from TAU2 to TAU, from
;        STEP_UP2 to STEP_UP, and from STEP_DN2 to STEP_DN, for
;        FEF_5 and FEF_6. 
;  Modification by Ken Jensen, Hughes STX, 29 Sept 1994, RATIO
;        array set to zero for Low Channels. (SPR 11925)


; Comparison Taus and Mission Periods (34 and 167)
; ------------------------------------------------
tau2 = [61.,153.]
step_up2 = ['89328','901391535','901931850']
step_dn2 = ['9026409','901931850','902081120']


; 34 Frequency Spectra
; --------------------

rec_len = 3572
	; fixed length record size in bytes

OPENW, 10, 'csdr$firas_out:fef_spc_llss.f16_93hybrid', rec_len, /FIXED

; FEF_1 Descriptions
; ------------------
desc0 = '4MP 1 tau vs 3MP 2 tau Comparison'
descrip1 = desc0 + '  LLSS'
descrip2 = desc0 + '  RLSS'
descrip3 = desc0 + '  LLFA'
descrip4 = desc0 + '  RLFA'

restore,'csdr$firas_save:lls_fef_1.iss'
chsn1='LLSS'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rls_fef_1.iss'
chsn1='RLSS'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lsf_fef_1.iss'
chsn1='LLFA'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rsf_fef_1.iss'
chsn1='RLFA'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_2 Descriptions
; ------------------
desc0 = '4MP 2 tau vs 3MP 2 tau Comparison'
descrip1 = desc0 + '  LLSS'
descrip2 = desc0 + '  RLSS'
descrip3 = desc0 + '  LLFA'
descrip4 = desc0 + '  RLFA'

restore,'csdr$firas_save:lls_fef_2.iss'
chsn1='LLSS'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rls_fef_2.iss'
chsn1='RLSS'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lsf_fef_2.iss'
chsn1='LLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rsf_fef_2.iss'
chsn1='RLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_3A Descriptions
; -------------------
desc0 = 'New Taus ( 50.2 , 140. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLSS'
descrip2 = desc0 + '  RLSS'
descrip3 = desc0 + '  LLFA'
descrip4 = desc0 + '  RLFA'

restore,'csdr$firas_save:lls_fef_3a.iss'
chsn1='LLSS'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rls_fef_3a.iss'
chsn1='RLSS'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lsf_fef_3a.iss'
chsn1='LLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rsf_fef_3a.iss'
chsn1='RLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_3B Descriptions
; -------------------
desc0 = 'New Taus ( 70.7 , 168.7 ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLSS'
descrip2 = desc0 + '  RLSS'
descrip3 = desc0 + '  LLFA'
descrip4 = desc0 + '  RLFA'

restore,'csdr$firas_save:lls_fef_3b.iss'
chsn1='LLSS'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rls_fef_3b.iss'
chsn1='RLSS'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lsf_fef_3b.iss'
chsn1='LLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rsf_fef_3b.iss'
chsn1='RLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_4A Descriptions
; -------------------
desc0 = '(2.3-4.5) icm Optimized Tau ( 1790. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLSS'
descrip2 = desc0 + '  RLSS'
descrip3 = desc0 + '  LLFA'
descrip4 = desc0 + '  RLFA'

restore,'csdr$firas_save:lls_fef_4a.iss'
chsn1='LLSS'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rls_fef_4a.iss'
chsn1='RLSS'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lsf_fef_4a.iss'
chsn1='LLFA'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rsf_fef_4a.iss'
chsn1='RLFA'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_4B Descriptions
; -------------------
desc0 = '(7.9-10.2) icm Optimized Tau ( 40. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLSS'
descrip2 = desc0 + '  RLSS'
descrip3 = desc0 + '  LLFA'
descrip4 = desc0 + '  RLFA'

restore,'csdr$firas_save:lls_fef_4b.iss'
chsn1='LLSS'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rls_fef_4b.iss'
chsn1='RLSS'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lsf_fef_4b.iss'
chsn1='LLFA'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rsf_fef_4b.iss'
chsn1='RLFA'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_5 Descriptions
; ------------------
descrip1 = 'RLSS LLSS Comparison  Std. Destriper'
descrip2 = 'RLFA LLFA Comparison  Std. Destriper'

restore,'csdr$firas_save:rls_fef_5.iss'
chsn1='RLSS' & chsn2='LLSS'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn

restore,'csdr$firas_save:rsf_fef_5.iss'
chsn1='RLFA' & chsn2='LLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn


; FEF_6 Descriptions
; ------------------
descrip1 = 'LLSS LLFA Comparison  Std. Destriper'
descrip2 = 'RLSS RLFA Comparison  Std. Destriper'

restore,'csdr$firas_save:lls_fef_6.iss'
chsn1='LLSS' & chsn2='LLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn

restore,'csdr$firas_save:rls_fef_6.iss'
chsn1='RLSS' & chsn2='RLFA'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn


CLOSE,10


; 167 Frequency Spectra
; ---------------------

OPENW, 10, 'csdr$firas_out:fef_spc_rhss.f16_93hybrid', rec_len, /FIXED


; FEF_1 Descriptions
; ------------------
desc0 = '4MP 1 tau vs 3MP 2 tau Comparison'
descrip1 = desc0 + '  LHSS'
descrip2 = desc0 + '  RHSS'
descrip3 = desc0 + '  LHFA'
descrip4 = desc0 + '  RHFA'

restore,'csdr$firas_save:lhs_fef_1.iss'
chsn1='LHSS'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhs_fef_1.iss'
chsn1='RHSS'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lhf_fef_1.iss'
chsn1='LHFA'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhf_fef_1.iss'
chsn1='RHFA'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_2 Descriptions
; -------------------
desc0 = '4MP 2 tau vs 3MP 2 tau Comparison'
descrip1 = desc0 + '  LHSS'
descrip2 = desc0 + '  RHSS'
descrip3 = desc0 + '  LHFA'
descrip4 = desc0 + '  RHFA'

restore,'csdr$firas_save:lhs_fef_2.iss'
chsn1='LHSS'
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhs_fef_2.iss'
chsn1='RHSS'
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lhf_fef_2.iss'
chsn1='LHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhf_fef_2.iss'
chsn1='RHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_3A Descriptions
; -------------------
desc0 = 'New Taus ( 50.2 , 140. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LHSS'
descrip2 = desc0 + '  RHSS'
descrip3 = desc0 + '  LHFA'
descrip4 = desc0 + '  RHFA'

restore,'csdr$firas_save:lhs_fef_3a.iss'
chsn1='LHSS'
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhs_fef_3a.iss'
chsn1='RHSS'
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lhf_fef_3a.iss'
chsn1='LHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhf_fef_3a.iss'
chsn1='RHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_3B Descriptions
; -------------------
desc0 = 'New Taus ( 70.7 , 168.7 ) vs Std Taus Comparison'
descrip1 = desc0 + '  LHSS'
descrip2 = desc0 + '  RHSS'
descrip3 = desc0 + '  LHFA'
descrip4 = desc0 + '  RHFA'

restore,'csdr$firas_save:lhs_fef_3b.iss'
chsn1='LHSS'
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhs_fef_3b.iss'
chsn1='RHSS'
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lhf_fef_3b.iss'
chsn1='LHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhf_fef_3b.iss'
chsn1='RHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_4C Descriptions
; -------------------
desc0 = '(28.9-33.4) icm Optimized Tau ( 20. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LHSS'
descrip2 = desc0 + '  RHSS'
descrip3 = desc0 + '  LHFA'
descrip4 = desc0 + '  RHFA'

restore,'csdr$firas_save:lhs_fef_4c.iss'
chsn1='LHSS'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhs_fef_4c.iss'
chsn1='RHSS'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:lhf_fef_4c.iss'
chsn1='LHFA'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip3,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rhf_fef_4c.iss'
chsn1='RHFA'  &  tau1 = [tau1,0.]
fef_write,diff,d_unc,ratio,r_unc,descrip4,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_5 Descriptions
; ------------------
descrip1 = 'RHSS LHSS Comparison  Std. Destriper'
descrip2 = 'RHFA LHFA Comparison  Std. Destriper'

restore,'csdr$firas_save:rhs_fef_5.iss'
chsn1='RHSS' & chsn2='LHSS'
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn

restore,'csdr$firas_save:rhf_fef_5.iss'
chsn1='RHFA' & chsn2='LHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn


; FEF_6 Descriptions
; ------------------
descrip1 = 'LHSS LHFA Comparison  Std. Destriper'
descrip2 = 'RHSS RHFA Comparison  Std. Destriper'

restore,'csdr$firas_save:lhs_fef_6.iss'
chsn1='LHSS' & chsn2='LHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn

restore,'csdr$firas_save:rhs_fef_6.iss'
chsn1='RHSS' & chsn2='RHFA'
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn


CLOSE,10


; Comparison Taus and Mission Periods (148)
; -----------------------------------------
tau2 = [61.,153.]
step_up2 = ['89328','901391535']
step_dn2 = ['9026409','901931850']


; 148 Frequency Spectra
; ---------------------

OPENW, 10, 'csdr$firas_out:fef_spc_lllf.f16_93hybrid', rec_len, /FIXED


; FEF_1 Descriptions
; ------------------
desc0 = '3MP 1 tau vs 2MP 2 tau Comparison'
descrip1 = desc0 + '  LLLF'
descrip2 = desc0 + '  RLLF'

restore,'csdr$firas_save:llf_fef_1.iss'
chsn1='LLLF'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rlf_fef_1.iss'
chsn1='RLLF'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_2 Descriptions
; ------------------
desc0 = '3MP 2 tau vs 2MP 2 tau Comparison'
descrip1 = desc0 + '  LLLF'
descrip2 = desc0 + '  RLLF'

restore,'csdr$firas_save:llf_fef_2.iss'
chsn1='LLLF'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rlf_fef_2.iss'
chsn1='RLLF'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_3A Descriptions
; -------------------
desc0 = 'New Taus ( 50.2 , 140. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLLF'
descrip2 = desc0 + '  RLLF'

restore,'csdr$firas_save:llf_fef_3a.iss'
chsn1='LLLF'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rlf_fef_3a.iss'
chsn1='RLLF'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_3B Descriptions
; -------------------
desc0 = 'New Taus ( 70.7 , 168.7 ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLLF'
descrip2 = desc0 + '  RLLF'

restore,'csdr$firas_save:llf_fef_3b.iss'
chsn1='LLLF'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rlf_fef_3b.iss'
chsn1='RLLF'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_4A Descriptions
; -------------------
desc0 = '(2.3-4.5) icm Optimized Tau ( 1790. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLLF'
descrip2 = desc0 + '  RLLF'

restore,'csdr$firas_save:llf_fef_4a.iss'
chsn1='LLLF'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rlf_fef_4a.iss'
chsn1='RLLF'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_4B Descriptions
; -------------------
desc0 = '(7.9-10.2) icm Optimized Tau ( 40. ) vs Std Taus Comparison'
descrip1 = desc0 + '  LLLF'
descrip2 = desc0 + '  RLLF'

restore,'csdr$firas_save:llf_fef_4b.iss'
chsn1='LLLF'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip1,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

restore,'csdr$firas_save:rlf_fef_4b.iss'
chsn1='RLLF'  &  tau1 = [tau1,0.]
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip2,chsn1,chsn1, $
          tau1,tau2,step_up1,step_dn1,step_up2,step_dn2


; FEF_5 Descriptions
; ------------------
descrip = 'RLLF LLLF Comparison  Std. Destriper'

restore,'csdr$firas_save:rlf_fef_5.iss'
chsn1='RLLF' & chsn2='LLLF'
ratio = 0. * ratio
fef_write,diff,d_unc,ratio,r_unc,descrip,chsn1,chsn2, $
          tau,tau,step_up,step_dn,step_up,step_dn


CLOSE,10


; Re-Define Restored Parameters
; -----------------------------
idd=0 & f=0

;
RETURN
END
