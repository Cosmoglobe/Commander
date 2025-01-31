FUNCTION FMD_MODEL_WGT,var2n=var2n,good=good,nifgs=nifgs,tm=tm,sky_s0=sky_s0,$
                       sky_glitch=sky_glitch,scan=scan,cal_nifgs=cal_nifgs,$
                       cal_tm=cal_tm,cal_glitch=cal_glitch,cal_s0=cal_s0,$
                       lp_order=lp_order,tim_up=tim_up,tim_dn=tim_dn,tcut=tcut,$
                       gpow=gpow,bpow=bpow,nscan=nscan,delt=delt,func=func,$
                       cal_func=cal_func,condition=condition,inv_stat=inv_stat,$
                       dtmod=dtmod,sky_wgts=sky_wgts,cal_wgts=cal_wgts,$
                       weightit=weightit,glcut=glcut,skut=skut,$
                       dscmod=dscmod,vfunc_descrip=vfunc_descrip
;


weightf = 'N'
IF (KEYWORD_SET(weightit)) THEN weightf = STRUPCASE(STRMID(weightit,0,1))
;

found=0
;

vfunc_descrip = ' All_Mission_Step'  ;  Initialize Variance Function descriptor
;

; IF (LP_ORDER) Get Time Functions (Legendre Polynomials)
; -------------------------------------------------------
;
IF (KEYWORD_SET(lp_order)) THEN BEGIN
 ;
 IF (MAX(lp_order) GT 0) THEN BEGIN
  ;
  smax = STRCOMPRESS(STRING(N_ELEMENTS(lp_order)+1),/remove_all)
  vfunc_descrip = ' ' + smax + '_Legendre'
  ;
  ; Coadd Times on a range of [-1,+1]
  ; ---------------------------------
  xtm = -1. + 2. * (tm - MIN(tm)) / (MAX(tm) - MIN(tm))
  func = FMD_LEG_FUNC(x=xtm,lp_order=lp_order)
  ;
  cxtm = -1. + 2. * (cal_tm - MIN(tm)) / (MAX(tm) - MIN(tm))
  cal_func = FMD_LEG_FUNC(x=cxtm,lp_order=lp_order)
  ;
  found = 1
  ;
 ENDIF
 ;
ENDIF
;

; IF (TIM_UP) Get Time Functions (Tophats)
; ----------------------------------------
;
IF (KEYWORD_SET(tim_up)) THEN BEGIN
 ;
 nstep = 0
 FOR i=0,N_ELEMENTS(tim_up)-1 DO BEGIN
  ;
  j = WHERE((tm(good) GE tim_up(i))and(tm(good) LE tim_dn(i)),cj1)
  print,'Number of Good Coadds in Special Time Range  = ' + STRCOMPRESS(STRING(cj1))
  ;
  j = WHERE((tm(good) LT tim_up(i))or(tm(good) GT tim_dn(i)),cj2)
  print,'Number of Good Coadds Outside Special Time Range  = ' + STRCOMPRESS(STRING(cj2))
  ;
  IF ((cj1 GT 100)and(cj2 GT 100)) THEN BEGIN
   ;
   IF (found EQ 0) THEN vfunc_descrip = ' '
   ;
   nstep = nstep + 1
   ;
   funcx = 0.*tm
   j = WHERE((tm GT tim_up(i))and(tm LE tim_dn(i)),cj)
   print,'Number of Coadds in Time Range  = ' + STRCOMPRESS(STRING(cj)) 
   IF (cj GT 0) THEN funcx(j) = 1.
   ;
   funcy = 0.*cal_tm
   j = WHERE((cal_tm GT tim_up(i))and(cal_tm LE tim_dn(i)),cj)
   IF (cj GT 0) THEN funcy(j) = 1.
   ;
   IF (found EQ 1) THEN func = [[func],[funcx]]
   IF (found EQ 1) THEN cal_func = [[cal_func],[funcy]]
   IF (found EQ 0) THEN BEGIN
    func = funcx
    cal_func = funcy
    found = 1
   ENDIF
   ;
  ENDIF
  ;
 ENDFOR
 ;
 sstep = STRCOMPRESS(STRING(nstep),/remove_all)
 IF (nstep GT 1) THEN vfunc_descrip = vfunc_descrip + '_+_' + sstep + '_Steps'
 IF (nstep EQ 1) THEN vfunc_descrip = vfunc_descrip + '_+_1_Step'
 ;
ENDIF
;

IF (KEYWORD_SET(nscan)) THEN BEGIN
 ;
 IF (nscan GT 0) THEN BEGIN
  ;
  vfunc_descrip = vfunc_descrip + ',' + $
    STRCOMPRESS(STRING(2*nscan+1),/remove_all) + '_Scan_Terms_per_Time_Func'
  ;
  xsc = !pi * scan / 180.
  ;
  IF (found EQ 1) THEN BEGIN
   ;
   funcz = func
   ;
   n_tfunc = N_ELEMENTS(funcz) / N_ELEMENTS(tm)
   funcx = 0. * funcz
   funcy = 0. * funcz
   FOR k=1,nscan DO BEGIN
    FOR kk=0,n_tfunc-1 DO funcx(*,kk) = funcz(*,kk) * COS(FLOAT(k)*xsc)
    FOR kk=0,n_tfunc-1 DO funcy(*,kk) = funcz(*,kk) * SIN(FLOAT(k)*xsc)
    func = [[func],[funcx],[funcy]]
   ENDFOR
   ;
   funcz = 0. * cal_func
   FOR k = 1,nscan DO cal_func = [[cal_func],[funcz],[funcz]]
   ;
  ENDIF
  ;
  IF (found EQ 0) THEN BEGIN
   ;
   func = 0. * tm + 1.
   FOR k = 1,nscan DO func = [[func],[COS(FLOAT(k)*xsc)],[SIN(FLOAT(k)*xsc)]]
   ;
   cal_func = 0. * cal_tm + 1.
   funcx = 0. * cal_func
   FOR k = 1,nscan DO cal_func = [[cal_func],[funcx],[funcx]]
   ;
   found = 1
   ;
  ENDIF
  ;
 ENDIF
 ;
ENDIF
;

; IF (GPOW) THEN Append Glitch Functions
; --------------------------------------
IF (KEYWORD_SET(gpow)) THEN BEGIN
 ;
 nq = WHERE(gpow GT 0.,cq)
 IF (cq GT 0) THEN BEGIN  ;  Glitch_Rate Modelling Enabled.
  ;
  IF (cq GT 1) THEN vfunc_descrip = vfunc_descrip + ',' + $
    STRCOMPRESS(STRING(cq),/remove_all) + '_Glitch_Terms'
  ;
  IF (cq EQ 1) THEN vfunc_descrip = vfunc_descrip + ',1_Glitch_Term'
  ;
  gpow = gpow(nq)
  ;
  IF (found EQ 0) THEN BEGIN
   func = 0. * tm + 1.
   cal_func = 0. * cal_tm + 1.
   found = 1
  ENDIF
  ;
  FOR i=0,N_ELEMENTS(gpow)-1 DO func = [[func],[sky_glitch^gpow(i)]]
  FOR i=0,N_ELEMENTS(gpow)-1 DO cal_func = [[cal_func],[cal_glitch^gpow(i)]]
  ;

  nq = WHERE((tm(good) GE tcut(0))and(tm(good) LE tcut(1)),cq)
  ;
  IF (cq GE 3) THEN BEGIN      ; Special High Glitch Functions
   ;
   funcx = 0.*sky_glitch
   nq = WHERE((tm GE tcut(0))and(tm LE tcut(1)),cq)
   funcx(nq) = sky_glitch(nq)
   func = [[func],[funcx],[funcx^2.]]
   ;
   funcx = 0. * cal_glitch
   nq = WHERE((cal_tm GE tcut(0))and(cal_tm LE tcut(1)),ccq)
   IF (ccq GE 1) THEN funcx(nq) = cal_glitch(nq)
   cal_func = [[cal_func],[funcx],[funcx^2.]]
   ;
  ENDIF
  ;

  IF ((cq LT 3)and(KEYWORD_SET(glcut))) THEN BEGIN  ; High Glitch Functions
   ;
   nq = WHERE(sky_glitch(good) GT glcut,cq)
   IF (cq GE 3) THEN BEGIN
    ;
    vfunc_descrip = vfunc_descrip + ',2_Special_Glitch_Terms'
    ;
    funcx = 0.*sky_glitch
    nq = WHERE(sky_glitch GT glcut,cq)
    funcx(nq) = sky_glitch(nq) - glcut
    func = [[func],[funcx],[funcx^2.]]
    ;
    funcx = 0. * cal_glitch
    nq = WHERE(cal_glitch GT glcut,cq)
    IF (cq GE 1) THEN funcx(nq) = cal_glitch(nq)
    cal_func = [[cal_func],[funcx],[funcx^2.]]
    ;
   ENDIF
   ;
  ENDIF
  ;
 ENDIF
 ;
ENDIF
;

; IF (BPOW) THEN Append Bolometer Functions
; -----------------------------------------
IF (KEYWORD_SET(bpow)) THEN BEGIN
 ;
 nq = WHERE(bpow GT 0.,cq)
 IF (cq GT 0) THEN BEGIN  ;  Bolometer Modelling Enabled.
  ;
  IF (cq GT 1) THEN vfunc_descrip = vfunc_descrip + ',' + $
    STRCOMPRESS(STRING(cq),/remove_all) + '_Bolometer_Terms'
  ;
  IF (cq EQ 1) THEN vfunc_descrip = vfunc_descrip + ',1_Bolometer_Term'
  ;
  bpow = bpow(nq)
  ;
  IF (found EQ 0) THEN BEGIN
   func = 0. * tm + 1.
   cal_func = 0. * cal_tm + 1.
   found = 1
  ENDIF
  ;
  funcx = ABS(sky_s0)
  avx = TOTAL(funcx) / N_ELEMENTS(funcx)
  ;
  funcy = (funcx - MIN(funcx)) / (avx - MIN(funcx))
  funcz = (ABS(cal_s0) - MIN(funcx)) / (avx - MIN(funcx))
  ;
  FOR i=0,N_ELEMENTS(bpow)-1 DO func = [[func],[funcy^bpow(i)]]
  FOR i=0,N_ELEMENTS(bpow)-1 DO cal_func = [[cal_func],[funcz^bpow(i)]]
  ;
 ENDIF
 ;
ENDIF
;

n_func = N_ELEMENTS(func) / N_ELEMENTS(tm)
;

; Matrix Elements
; ---------------
mm = DBLARR(n_func,n_func)
;
IF (weightf EQ 'Y') THEN $
    FOR i=0,n_func-1 DO FOR j=i,n_func-1 DO $
        mm(i,j) = TOTAL(func(good,i)*func(good,j)*nifgs(good))
;
IF (weightf NE 'Y') THEN $
    FOR i=0,n_func-1 DO FOR j=i,n_func-1 DO $
        mm(i,j) = TOTAL(func(good,i)*func(good,j))
;
FOR i=0,n_func-1 DO mm(i,i) = 0.5 * mm(i,i)
mm = mm + TRANSPOSE(mm)
;

; Variance Vector
; ---------------
vec = DBLARR(n_func)
;
IF (weightf EQ 'Y') THEN $
    FOR i=0,n_func-1 DO vec(i) = TOTAL(var2n*func(good,i)*nifgs(good))
;
IF (weightf NE 'Y') THEN $
    FOR i=0,n_func-1 DO vec(i) = TOTAL(var2n*func(good,i))
;

; Matrix Condition Number
; -----------------------
NR_SVD,mm,w,u,v,/double
condition = MAX(w)/MIN(w)
;

; Matrix Inversion
; ----------------
mi = INVERT(mm,inv_stat)
;

; Solution Vector
; ---------------
avec = mi # vec
;

; Variance Model
; --------------
vmod = func # avec
;

; Variance - Model Functions (Good Coadds)
; ---------------------------------------
vres = var2n - vmod(good)
;

; Small Timescale Correction
; --------------------------
dtmod = DOUBLE(tm) * 0.
;
IF (KEYWORD_SET(delt)) THEN BEGIN
 ;
 IF (delt GT 0.) THEN BEGIN
  ;
  vfunc_descrip = vfunc_descrip + ',' + $
    STRCOMPRESS(STRING(FIX(delt)),/remove_all) + '_Day_Time_Correction'
  ;
  tmg = tm(good)
  FOR i=long(0),long(N_ELEMENTS(tm))-1 DO BEGIN
   nx = WHERE((ABS(tmg-tm(i))) LE delt,cx)
   IF (cx GT 0) THEN dtmod(i) = TOTAL(vres(nx))/cx
  ENDFOR
  ;
 ENDIF
 ;
ENDIF
;

; Special Scan_Angle Correction
; -----------------------------
dscmod = DOUBLE(scan) * 0.
;
IF (KEYWORD_SET(skut)) THEN BEGIN
 ;
 IF (skut(4) GT 0) THEN BEGIN
  ;
  tmg = tm(good)
  scg = scan(good)
  nx = $
   WHERE((tmg GT skut(0))and(tmg LT skut(1))and(scg GT skut(2))and(scg LT skut(3)),cx)
  IF (cx GE skut(4)) THEN BEGIN
   ;
   PRINT,' '
   PRINT,'Special Scan_Angle Correction'
   PRINT,' '
   PRINT,STRCOMPRESS(STRING(cx)) + '  Good Variance Coadds in Range'
   PRINT,' '
   ;
   vfunc_descrip = vfunc_descrip + ',Special_Scan_Angle_Correction'
   ;
   scg = scan(good(nx))
   vresg = vres(nx)
   ;
   nx = $
    WHERE((tm GT skut(0))and(tm LT skut(1))and(scan GT skut(2))and(scan LT skut(3)),cx)
   PRINT,STRCOMPRESS(STRING(cx),/remove_all) + '  Good Coadds in Range'
   PRINT,' '
   scx = scan(nx)
   ;
   PRINT,STRCOMPRESS(STRING(skut(4))) + '  Good Variance Coadds Averaged'
   PRINT,'     for each Good Coadd Correction'
   PRINT,' '
   ;
   FOR i=0,cx-1 DO BEGIN
    dels = ABS(scg - scx(i))
    sx = SORT(dels)
    dscmod(nx(i)) = TOTAL(vresg(sx(0:skut(4)-1))) / skut(4)
   ENDFOR
   ;
  ENDIF
  ;
 ENDIF
 ;
ENDIF
;

; Adjusted Coadd Weights
; ----------------------
ng = WHERE(sky_wgts GT 0.)
sky_wgts(ng) = FLOAT(nifgs(ng)) / (vmod(ng) + dtmod(ng) + dscmod(ng))
;

; Adjusted Cal Weights
; --------------------
vmod = cal_func # avec
ng = WHERE(cal_wgts GT 0.)
cal_wgts(ng) = FLOAT(cal_nifgs(ng)) / vmod(ng)
;

RETURN,avec
END
