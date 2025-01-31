FUNCTION FMD_LEG_FUNC,x=x,lp_order=lp_order
;

; Legendre Polynomial Kernels
; ---------------------------
;
; Zeroth-Order
; ------------
func = FLTARR(N_ELEMENTS(x)) + 1.
;

IF (NOT KEYWORD_SET(lp_order)) THEN RETURN,func
;

dummy = WHERE(lp_order GT 0,n_leg)
;
IF (n_leg LE 0) THEN RETURN,func
;

lp_order = lp_order(dummy)
;
lpmax = MAX(lp_order)       ; Maximum Order of Legendre Polynomials
;

; First-Order and Second-Order Legendre Polynomials
; -------------------------------------------------
p1 = x
p2 = 0.5*(3.*x^2.-1.)
;

lp_count = 0  ; Kernal counter
;

np = WHERE(lp_order EQ 1,cnp)
IF (cnp EQ 1) THEN BEGIN
 func = [[func],[p1]]
 lp_count = 1
ENDIF
;
np = WHERE(lp_order EQ 2,cnp)
IF (cnp EQ 1) THEN BEGIN
 func = [[func],[p2]]
 lp_count = lp_count + 1
ENDIF
;

; IF (lpmax GT 2) Append additional Legendre Polynomials
; ------------------------------------------------------
IF (lpmax GT 2) THEN BEGIN
 ;
 FOR n=3,lpmax DO BEGIN
  ;
  p3 = 2.*x*p2 - p1 - (x*p2 - p1)/n  ; Nth-Order Legendre Polynomial
  ;
  ; IF Nth-Order is in LP_ORDER, Append the Nth-Order Polynomial
  ; ------------------------------------------------------------
  np = WHERE(lp_order EQ n,cnp)
  IF (cnp EQ 1) THEN BEGIN
   func = [[func],[p3]]
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

; Check on Legendre Kernals
; -------------------------
IF (lp_count NE n_leg) THEN BEGIN
 PRINT,'FMD_LEG_FUNC  :  LP_COUNT NE N_LEG !!!'
 RETURN,func
ENDIF
;


RETURN,func
END
