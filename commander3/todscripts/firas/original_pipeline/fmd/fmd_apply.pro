FUNCTION FMD_APPLY, freq=freq,func=func,n_krnl=n_krnl,n_cmn=n_cmn, $
                    spec=spec,cal_spec=cal_spec, $
                    diff_spec=diff_spec,diff_cal_spec=diff_cal_spec, $
 		    ejg=ejg,sky_idx=sky_idx,cal_idx=cal_idx, $
                    dsp=dsp,cal_dsp=cal_dsp,pure_off=pure_off

;
; FMD_APPLY applies the destriper corrections to the coadd spectra.
;
; Written by Ken Jensen, Hughes STX, 3-Oct-96
;
;

; Set up various utility variables
; --------------------------------
n_freq = N_ELEMENTS(freq)	        ; # of frequency points
n_sky = N_ELEMENTS(spec) / n_freq       ; # of sky spectra
n_cal = N_ELEMENTS(cal_spec) / n_freq	; # of cal spectra
n_obs = n_sky + n_cal    	        ; Total # of spectra
n_stripes = N_ELEMENTS(ejg) / n_freq    ; # of stripes (incl. gain stripes)
;

; Determine total # of functions
; ------------------------------
;
IF (KEYWORD_SET(sky_idx)) THEN $
    n_chan = MAX(sky_idx) + 1	        ; # of channels
IF (NOT KEYWORD_SET(sky_idx)) THEN $
    n_chan = 1                    	; # of channels
;
n_uncmn = n_krnl - n_cmn
n_off = n_uncmn * n_chan + n_cmn
;
n_gain = 0
IF (NOT KEYWORD_SET(pure_off)) THEN n_gain = n_chan
;
n_func = n_off + n_gain
;
IF (NOT KEYWORD_SET(sky_idx)) THEN sky_idx = BYTARR(n_sky)
IF (NOT KEYWORD_SET(cal_idx)) THEN cal_idx = BYTARR(n_cal)
;

; Offset EJG
; ----------
ejg_off = FLTARR(n_freq,n_off)
;
FOR i=0,n_off-1 DO ejg_off(*,i) = ejg(*,i)
;

; Gain EJG
; --------
ejg_gain = FLTARR(n_freq,n_chan)
FOR i=0,n_chan-1 DO ejg_gain(*,i) = ejg(*,n_off+i)
;

; DSP = Destriped Sky Spectra
; ---------------------------
sp_temp = TRANSPOSE(spec)
dsp = 0. * sp_temp
;

; CAL_DSP = Destriped Cal Spectra
; -------------------------------
csp_temp = TRANSPOSE(cal_spec)
cal_dsp = 0. * csp_temp
;

; EJG1 = Offset EJG for one CHANSCAN
; ----------------------------------
ejg1 = FLTARR(n_freq,n_krnl)
;

; EJG2 = Gain EJG for one CHANSCAN
; --------------------------------
ejg2 = FLTARR(n_freq)
;

; Common Stripes
; --------------
IF(n_cmn ne 0) THEN FOR j=0,n_cmn-1 DO ejg1(*,j) = ejg_off(*,j)
;

; Loop over CHANSCANs
; -------------------
;
FOR ch=0,n_chan-1 DO BEGIN
;
; Channel-Specific Offset Stripes
; -------------------------------
 FOR j=n_cmn,n_krnl-1 DO ejg1(*,j) = ejg_off(*,j+(n_krnl-n_cmn)*ch)
;

; Gain Stripe
; -----------
 ejg2(*) = ejg(*,n_off+ch)
;

; Apply Stripes
; -------------
 nq = WHERE(cal_idx eq ch,cq)
 FOR i=long(0),long(cq)-1 DO BEGIN
  cal_dsp(*,nq(i)) = csp_temp(*,nq(i)) - ejg1 # func(*,n_sky+nq(i))
  IF (NOT KEYWORD_SET(pure_off)) THEN $
     cal_dsp(*,nq(i)) = cal_dsp(*,nq(i)) - (1.-ejg2(*)) * diff_cal_spec(nq(i),*)
 ENDFOR
;
 nq = WHERE(sky_idx eq ch,cq)
 FOR i=long(0),long(cq)-1 DO BEGIN
  dsp(*,nq(i)) = sp_temp(*,nq(i)) - ejg1 # func(*,nq(i))
  IF (NOT KEYWORD_SET(pure_off)) THEN $
     dsp(*,nq(i)) = dsp(*,nq(i)) - (1.-ejg2(*)) * diff_spec(nq(i),*)
 ENDFOR
;

ENDFOR
;

RETURN,n_freq
END
