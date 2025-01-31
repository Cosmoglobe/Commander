FUNCTION FMD_FUNC, g8=g8,g9=g9,g10=g10,dirbe_array=dirbe_array,$
                   tm=tm,cal_tm=cal_tm,lp_order=lp_order,xtm=xtm,$
                   step_up=step_up,step_dn=step_dn,$
  		   sky_dihd=sky_dihd,cal_dihd=cal_dihd,ref_dihd=ref_dihd,$
	           dihed_pow=dihed_pow,dihed_cut=dihed_cut,$
                   sky_s0=sky_s0,cal_s0=cal_s0,ref_s0=ref_s0,rms_s0=rms_s0,$
                   sky_idx=sky_idx,cal_idx=cal_idx,$
                   n_dirbe=n_dirbe,n_time=n_time,n_tophat=n_tophat,$
                   n_dihd=n_dihd,n_bol=n_bol,n_chan=n_chan,$
                   stripe_order=stripe_order

;
; FMD_FUNC builds the destriper kernels
;
; Calls      :  TIMECONV
;
; Written by :  Ken Jensen,  Hughes STX,  03-Oct-96
;
;

; Sanity Checks and defaults
; --------------------------
;
IF (NOT KEYWORD_SET(step_up)) THEN BEGIN
    PRINT,'FMD_FUNC : Keyword STEP_UP Must Be Present'
    func = -1
    RETURN,func
ENDIF
;
IF (NOT KEYWORD_SET(step_dn)) THEN BEGIN
    PRINT,'FMD_FUNC : Keyword STEP_DN Must Be Present'
    func = -1
    RETURN,func
ENDIF
;
IF (N_ELEMENTS(step_up) NE N_ELEMENTS(step_dn)) THEN BEGIN
    PRINT,'FMD_FUNC : STEP_UP and STEP_DN Must Have Equal Dimension'
    func = -1
    RETURN,func
ENDIF
;

IF (NOT KEYWORD_SET(dirbe_array)) THEN dirbe_array = [0,0,0]
;
dummy = WHERE(lp_order GT 0,n_tm)
;
n_sky = N_ELEMENTS(tm)
	; # of sky spectra

n_cal = N_ELEMENTS(cal_tm)
	; # of cal spectra

n_obs = n_sky + n_cal
	; Total # of spectra

n_type = 0  ; Number of Kernel "types" (DIRBE,TIME,DIHEDRAL,BOLOMETER,TOPHAT)
;

; DIRBE Gradient Kernels
; ----------------------
n_dirbe = 0   ;  Initialize Number of DIRBE kernels = 0
;
IF (MAX(dirbe_array) gt 0) THEN BEGIN  ;  Make DIRBE kernels
;
 n_type = n_type + 1
;
 IF (dirbe_array(0) eq 1) THEN BEGIN         ; Band 10
  func = [g10,0*cal_tm]                ; Band 10 kernel
  n_dirbe = n_dirbe + 1      ; Increment Number of DIRBE kernels
 ENDIF
;
 IF (dirbe_array(1) eq 1) THEN BEGIN   ; Band 9
  IF (n_dirbe eq 0) THEN func = [g9,0*cal_tm]   ; Band 9 kernel
  IF (n_dirbe eq 1) THEN $
   func = [[func],[g9,0*cal_tm]]          ; Add Band 9 kernel
  n_dirbe = n_dirbe + 1      ; Increment Number of DIRBE kernels
 ENDIF
;
 IF (dirbe_array(2) eq 1) THEN BEGIN   ; Band 8
  IF (n_dirbe eq 0) THEN func = [g8,0*cal_tm]   ; Band 8 kernel
  IF (n_dirbe ge 1) THEN $
   func = [[func],[g8,0*cal_tm]]          ; Add Band 8 kernel
  n_dirbe = n_dirbe + 1      ; Increment Number of DIRBE kernels
 ENDIF
;
 IF (stripe_order(0) EQ 1) THEN func1 = TEMPORARY(func)
 IF (stripe_order(0) EQ 2) THEN func2 = TEMPORARY(func)
 IF (stripe_order(0) EQ 3) THEN func3 = TEMPORARY(func)
 IF (stripe_order(0) EQ 4) THEN func4 = TEMPORARY(func)
 IF (stripe_order(0) EQ 5) THEN func5 = TEMPORARY(func)
;
ENDIF
;


; Time Kernels
; ------------
n_time = 0
;
; Legendre Polynomial Kernels
; ---------------------------
IF (n_tm GT 0) THEN BEGIN
;
 n_type = n_type + 1
;
 tme = [tm,cal_tm]            ; Observation Times (Days)
;
 xtm = 1.*xtm
 tmx = xtm(0) + tme * (xtm(1)-xtm(0)) / MAX(tme)
   ; TMX varies linearly  from XTM(0) at TME = 0  to XTM(1) at TME = MAX(TME)
;

; Maximum Order of legendre Polynomials
; -------------------------------------
 lpmax = MAX(lp_order)
;

; First-Order and Second-Order Legendre Polynomials
; -------------------------------------------------
 p1 = tmx
 p2 = 0.5*(3.*tmx^2.-1.)
;
 found = 0     ; Sets to 1 when the first kernal is built
 lp_count = 0  ; Kernal counter
;
 np = WHERE(lp_order EQ 1,cnp)
 IF (cnp EQ 1) THEN BEGIN
  func = p1
  found = 1
  lp_count = 1
 ENDIF
;

 np = WHERE(lp_order EQ 2,cnp)
 IF (cnp EQ 1) THEN BEGIN
  IF (found EQ 1) THEN func = [[func],[p2]]
  IF (found EQ 0) THEN BEGIN
   func = p2
   found = 1
  ENDIF
  lp_count = lp_count + 1
 ENDIF
;

; IF (lpmax GT 2) Append additional Legendre Polynomials
; ------------------------------------------------------
 IF (lpmax GT 2) THEN BEGIN
;
  FOR n=3,lpmax DO BEGIN
;
   p3 = 2.*tmx*p2 - p1 - (tmx*p2 - p1)/n  ; Nth-Order Legendre Polynomial
;
;  IF Nth-Order is in LP_ORDER, Append the Nth-Order Polynomial
;  ------------------------------------------------------------
   np = WHERE(lp_order EQ n,cnp)
   IF (cnp EQ 1) THEN BEGIN
    IF (found EQ 1) THEN func = [[func],[p3]]
    IF (found EQ 0) THEN BEGIN
     func = p3
     found = 1
    ENDIF
    lp_count = lp_count + 1
   ENDIF
;
   p1=p2  ; Redefine P1 for next recursion
   p2=p3  ; Redefine P2 for next recursion
;
  ENDFOR
;
 ENDIF
;

; Check on Time Kernals
; ---------------------
 IF (lp_count NE n_tm) THEN BEGIN
  PRINT,'FMD_FUNC : LP_COUNT NE N_TM !!!'
  RETURN,func
 ENDIF
;

; Number of Time stripes
; ----------------------
 n_time = n_time + n_tm 
;

; Place Time Kernals in Correct-Order
; -----------------------------------
 IF (stripe_order(1) EQ 1) THEN func1 = TEMPORARY(func)
 IF (stripe_order(1) EQ 2) THEN func2 = TEMPORARY(func)
 IF (stripe_order(1) EQ 3) THEN func3 = TEMPORARY(func)
 IF (stripe_order(1) EQ 4) THEN func4 = TEMPORARY(func)
 IF (stripe_order(1) EQ 5) THEN func5 = TEMPORARY(func)
;

ENDIF
;


; DIHEDRAL Kernels
; ----------------
;
; Dihedral Cut
; ------------
IF (KEYWORD_SET(dihed_cut)) THEN $
	split_diheds = 1 $
ELSE $
	split_diheds = 0
;
n_dihd = 0
;
IF (TOTAL(dihed_pow) GT 0.) THEN BEGIN
;
 n_type = n_type + 1
;
 n_dihd = N_ELEMENTS(dihed_pow)
 dihd0 = [sky_dihd,cal_dihd]
 dihd = dihd0^dihed_pow(0) - ref_dihd^dihed_pow(0)
 IF (n_dihd GT 1) THEN BEGIN
  FOR i=1,n_dihd-1 DO BEGIN
   dihdx = dihd0^dihed_pow(i)
   dihd = [[dihd],[dihdx - ref_dihd^dihed_pow(i)]]
  ENDFOR
 ENDIF
;
; If splitting the dihedrals is desired, do it.
;
 IF (split_diheds) THEN BEGIN
  cool_diheds = where (dihd0 lt dihed_cut)
  warm_diheds = where (dihd0 ge dihed_cut)
  n_dihd0 = n_dihd + n_dihd
;
;  Put in a little insurance just in case the cut was nonsense.
  IF ((cool_diheds(0) eq -1) OR (warm_diheds(0) eq -1)) THEN BEGIN
   split_diheds = 0
   n_dihd0 = n_dihd
   print, ' The dihedral cut you requested did not split the data.'
   print, '   Consequently only one dihedral kernel will be used.' 
  ENDIF
;
  n_dihd = n_dihd0
;
 ENDIF
;
 IF (split_diheds) THEN $
        func = [[dihd*(dihd0 LT dihed_cut)],$
                [dihd*(dihd0 GE dihed_cut)]] $
 ELSE $
        func = TEMPORARY(dihd)
;
 IF (stripe_order(2) EQ 1) THEN func1 = TEMPORARY(func)
 IF (stripe_order(2) EQ 2) THEN func2 = TEMPORARY(func)
 IF (stripe_order(2) EQ 3) THEN func3 = TEMPORARY(func)
 IF (stripe_order(2) EQ 4) THEN func4 = TEMPORARY(func)
 IF (stripe_order(2) EQ 5) THEN func5 = TEMPORARY(func)
;
ENDIF
;


; BOLOMETER Kernels
; -----------------
IF (n_bol GT 0) THEN BEGIN
;
 n_type = n_type + 1
;
 bol = [sky_s0,cal_s0]
 refb = 0.*bol
;
 IF (n_chan GT 1) THEN BEGIN
  idx = [sky_idx,cal_idx]
  refb(*) = ref_s0(idx(*))
 ENDIF
;
 IF (n_chan EQ 1) THEN refb(*) = ref_s0(0)
 func0 = (bol - refb) / rms_s0
;
 IF (n_bol EQ 2) THEN $
         func = [[func0*((idx eq 0)or(idx EQ 2))],$
                 [func0*((idx eq 1)or(idx EQ 3))]] $
 ELSE $
         func = TEMPORARY(func0)
;
 IF (stripe_order(3) EQ 1) THEN func1 = TEMPORARY(func)
 IF (stripe_order(3) EQ 2) THEN func2 = TEMPORARY(func)
 IF (stripe_order(3) EQ 3) THEN func3 = TEMPORARY(func)
 IF (stripe_order(3) EQ 4) THEN func4 = TEMPORARY(func)
 IF (stripe_order(3) EQ 5) THEN func5 = TEMPORARY(func)
;
ENDIF
;


; TOPHAT Kernels
; --------------
;
tophat_tme = [tm,cal_tm]     ; TM Array for Tophat Selection
;
eject = TIMECONV('893251118',infmt='z',outfmt='s')
sup = (TIMECONV(step_up,infmt='z',outfmt='s') - eject) / 86400.
sdn = (TIMECONV(step_dn,infmt='z',outfmt='s') - eject) / 86400.
;
n_tophat = N_ELEMENTS(step_up)
func = DBLARR(n_obs,n_tophat)
FOR i=0,n_tophat-1 DO BEGIN
 nq = WHERE( (tophat_tme GE sup(i)) AND (tophat_tme LT sdn(i)), cq )
 IF (cq GT 0) THEN func(nq,i) = 1
ENDFOR
;
IF (stripe_order(4) EQ 1) THEN func1 = TEMPORARY(func)
IF (stripe_order(4) EQ 2) THEN func2 = TEMPORARY(func)
IF (stripe_order(4) EQ 3) THEN func3 = TEMPORARY(func)
IF (stripe_order(4) EQ 4) THEN func4 = TEMPORARY(func)
IF (stripe_order(4) EQ 5) THEN func5 = TEMPORARY(func)
;
n_type = n_type + 1
;

; Build FUNC array from DIRBE, Time, Dihedral, Bolometer, and Tophat Kernels
; --------------------------------------------------------------------------
;
IF (n_type EQ 1) THEN func = TEMPORARY(func1)
;
IF (n_type EQ 2) THEN func = [[TEMPORARY(func1)],[TEMPORARY(func2)]]
;
IF (n_type EQ 3) THEN $
    func = [[TEMPORARY(func1)],[TEMPORARY(func2)],[TEMPORARY(func3)]]
;
IF (n_type EQ 4) THEN $
    func = [[TEMPORARY(func1)],[TEMPORARY(func2)],$
            [TEMPORARY(func3)],[TEMPORARY(func4)]]
;
IF (n_type EQ 5) THEN $
    func = [[TEMPORARY(func1)],[TEMPORARY(func2)],[TEMPORARY(func3)],$
            [TEMPORARY(func4)],[TEMPORARY(func5)]]
;

	; Note: FUNC is a "compressed" version of the "true" kernel
	; array of the form:
	;
        ;    | FUNC(cmn)     0         0         0    .      0    |
	;    |	   0	 FUNC(ch1)     0         0    .      0    |
	;    |	   .	     0     FUNC(ch2)     0    .      0    |
	;    |	   .	     .         .         .    .      .    |
	;    |	   0	     0         0         0    .  FUNC(chN)|
	;
	; ie:
	;
	;	FUNC = 	| FUNC(cmn) FUNC(ch1) FUNC(ch2) FUNC(ch3) . FUNC(chN)|


RETURN,func
END
