FUNCTION fla_fitdust,freq,spec,weight,f_char=f_char,parm0=parm0, $
         active=active,mask=mask,out_parm=out_parm,parm_sig=parm_sig, $
	 err_mat=err_mat
;
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Feb 1995
;    Initial Delivery  SPR XXXXX
;-


sz = SIZE(spec)
sz_weight = SIZE(weight)
n_freq = N_ELEMENTS(freq)


tolerance = 1.e-4

out_parm = FLTARR(7)
out_fit = FLTARR(n_freq)
parm = parm0

res_units = 'E'
funits = 'I'

m2ep = 2.99792458d-7	; megajanskies to eplees
act_parm = WHERE(active GT 0)
		; get active (variable) paramater indicies

n_parm = N_ELEMENTS(act_parm)
n_freq = N_ELEMENTS(freq)
ndf = n_freq - n_parm
		; get # active parms, # frequencies, # dof

s_fq = freq / f_char

sqrt_wt = SQRT(weight)
IF (N_ELEMENTS(mask) GT 0) THEN sqrt_wt(mask) = 0
		; scale frequencies
		; get sqrt(weights)
		; deweight desired points


plnck = DBLARR(n_freq,3)
dbdt = DBLARR(n_freq,3)
wt_derv = DBLARR(n_freq,n_parm)
parm_sig = DBLARR(7)
		; allocate stored planck, dbdt, weighted derv 
		; and parameter uncertainties arrays

FOR i=0,2 DO BEGIN
	IF (parm(i*i) NE 0.0) THEN BEGIN
		CASE STRUPCASE(STRMID(res_units,0,1)) OF

		'E' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy) * m2ep
			dbdt(*,i) = db * m2ep
			END

		'W' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/wcm2)
			dbdt(*,i) = db
			END

		'M' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy)
			dbdt(*,i) = db
			END

		ENDCASE

	ENDIF
ENDFOR
		; for temp parms if temp ne 0 then calculate and
		; store planck and dbdt


ft = plnck(*,0) + plnck(*,1) * s_fq^parm(3) * parm(2) + $
                  plnck(*,2) * s_fq^parm(6) * parm(5)
		; calculate initial fit

chi2 = TOTAL(((spec-ft)*sqrt_wt)^2)/ndf
		; calculate initial chi_squared


n_iter = 0
zero_temp = 0
REPEAT BEGIN
		; interation fit loop

n_iter = n_iter + 1

FOR i=0,n_parm-1 DO BEGIN
	CASE act_parm(i) OF
		0: derv = dbdt(*,0)
		1: derv = parm(2) * s_fq^parm(3) * dbdt(*,1)
		2: derv = s_fq^parm(3) * plnck(*,1)
		3: derv = parm(2) * s_fq^parm(3) * plnck(*,1) * alog(s_fq)
		4: derv = parm(5) * s_fq^parm(6) * dbdt(*,2)
		5: derv = s_fq^parm(6) * plnck(*,2)
		6: derv = parm(5) * s_fq^parm(6) * plnck(*,2) * alog(s_fq)
	ENDCASE
	wt_derv(*,i) = derv * sqrt_wt
ENDFOR
		; calculate d(ft)/d(parm) for all active parms

wt_spec_diff = REFORM((spec-ft)*sqrt_wt,1,n_freq)
		; calculate weighted difference between spec and fit

beta = REFORM(wt_spec_diff # wt_derv)
alpha = TRANSPOSE(wt_derv) # wt_derv
		; build fit vector (beta) and array (alpha)

IF (n_parm EQ 1) THEN BEGIN
	err_mat = 1 / alpha(0,0)
	del_parm = beta(0) * err_mat
ENDIF ELSE BEGIN
	err_mat = INVERT(alpha)
	del_parm = err_mat # beta
ENDELSE
		; calculate error matrix (=inverse of alpha)
		; calculate change in parameters

FOR j=0,n_parm-1 DO parm(act_parm(j)) = parm(act_parm(j)) + del_parm(j)
		; calculate new parameters

IF (parm(1) LT 0) THEN zero_temp = -1

FOR i=0,2 DO BEGIN
	IF (active(i*i) NE 0) THEN BEGIN
		CASE STRUPCASE(STRMID(res_units,0,1)) OF

		'E' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy) * m2ep
			dbdt(*,i) = db * m2ep
			END

		'W' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/wcm2)
			dbdt(*,i) = db
			END

		'M' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy)
			dbdt(*,i) = db
			END

		ENDCASE
	ENDIF
ENDFOR
		; for temp parms if active parm then calculate and
		; store planck and dbdt


ft = plnck(*,0) + plnck(*,1) * s_fq^parm(3) * parm(2) + $
                  plnck(*,2) * s_fq^parm(6) * parm(5)
		; calculate new fit

last_chi2 = chi2
chi2 = TOTAL(((spec-ft)*sqrt_wt)^2)/ndf
		; calculate new chi_squared

ENDREP UNTIL (ABS((chi2-last_chi2)/last_chi2) LT tolerance OR $
	      n_iter GE 20 OR zero_temp EQ -1)
		; repeat if relative change in chi2 gt tolerance or
		; dust temperature becomes negative


FOR j=0,n_parm-1 DO parm_sig(act_parm(j)) = SQRT(err_mat(j,j))
		; calculate parameter uncertainties


y_fit = ft

out_fit = y_fit
out_parm = parm


RETURN,out_fit
END
