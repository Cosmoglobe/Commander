PRO fef_write,diff,d_unc,ratio,r_unc,descrip,chsn1,chsn2, $
	      tau1,tau2,step_up1,step_dn1,step_up2,step_dn2

;
;  Modification History:
;
;   Written by Joel Gales, ARC, Sept 1994
;
;   Ken Jensen, Hughes STX, 19 Sept 1994,  FEF_STRUCT.TAU2 = TAU2
;       replaces FEF_STRUCT.TAU1 = TAU2 statement.
;
;   Ken Jensen, Hughes STX, 29 Sept 1994,  R_UNC field is zeroed.
;       RATIO field is zeroed for Low Channels. (SPR 11925)



; Define, Load and Store FEF dataset structure
; --------------------------------------------

fef_struct = {DESCRIP: STRING(' ',FORMAT='(A80)'), $
	      CHANSCAN1: STRING(' ',FORMAT='(A4)'), $
              START_TIMES1: LONARR(2,4), $
              STOP_TIMES1: LONARR(2,4), $
	      TAU1: FLTARR(2), $
	      CHANSCAN2: STRING(' ',FORMAT='(A4)'), $
              START_TIMES2: LONARR(2,4), $
              STOP_TIMES2: LONARR(2,4), $
	      TAU2: FLTARR(2), $
	      DIFF_SPEC: COMPLEXARR(167), $
	      DIFF_SPEC_SIG: FLTARR(167), $
	      RATIO_SPEC: FLTARR(167), $
	      RATIO_SPEC_SIG: FLTARR(167)}



fef_struct.descrip = STRING(descrip,FORMAT='(A80)')

fef_struct.chanscan1 = STRING(chsn1,FORMAT='(A4)')
fef_struct.start_times1 = timeconv(step_up1,infmt='z',outfmt='a')
fef_struct.stop_times1 = timeconv(step_dn1,infmt='z',outfmt='a')
fef_struct.tau1 = tau1

fef_struct.chanscan2 = STRING(chsn2,FORMAT='(A4)')
fef_struct.start_times2 = timeconv(step_up2,infmt='z',outfmt='a')
fef_struct.stop_times2 = timeconv(step_dn2,infmt='z',outfmt='a')
fef_struct.tau2 = tau2

fef_struct.diff_spec = diff
fef_struct.diff_spec_sig = d_unc
fef_struct.ratio_spec = ratio
fef_struct.ratio_spec_sig = r_unc * 0.


WRITEU, 10, fef_struct


RETURN
END
