FUNCTION FLA_FIT,freq,spec,line_wgt,fit_func=fit_func, $
		 nifgs=nifgs,dst_sig=dst_sig,n_lines=n_lines, $
		 mask=mask,dst_wgt=dst_wgt,active=active,parm0=parm0, $
		 lne_parms=lne_parms,dst_parms=dst_parms,dst_covar=dst_covar, $
		 f_char=f_char,pcvr=pcvr,dust_spec=dust_spec

;
; This program generates fluxes at a given set of line frequencies.
; A CMB and dust continuum model is initially fit the line_subtracted
; spectra and this fit the subtracted.  The residual spectra are then
; fit to a set of line profiles plus a "smooth" baseline.  Both the 
; real and imaginary parts of the spectra are fit simultaneously and
; the fit is weighted via a weights matrix that can include off-diagonal
; correlations.
;
; The user must supply the frequency vector, the spectra, the weights
; matrix and the line frequencies.  The user must also supply the
; order of the baseline.
;
; In order to fit the dust continuum the user must supply a mask vector
; giving those frequency bins that should be completely deweighted for
; this fit along with a weights vector, initial dust fitting parameters,
; and and vector specifying which parameters are fixed and which vary.
;
; The output consists of the CMB and dust continuum parmeters (stored
; in dst_PARMS), the line fluxes and baseline parameters, (stored in
; LNE_PARMS), and the total chi_squared (summed over frequency) for each
; spectra, passes as the output of this function.
;
; 
;  Modification History:
;           Written by Joel Gales, ARC, October 1994
;
;           Modified to call FLA_FITDUST instead of FITSPEC in order
;           to return covariance as well as variance.
;           J.M. Gales, ARC, February 1995
;
;           K.Jensen, Hughes STX, 2-Mar-1995, eliminate most of the debug
;           messages.
;-



; Get input array size info
; -------------------------
sz = SIZE(spec)
IF (sz(0) EQ 2) THEN BEGIN
	n_freq = sz(1)
	n_spec = sz(2)
ENDIF ELSE BEGIN
	n_freq = sz(1)
	n_spec = 1
ENDELSE

type = sz(sz(0)+1)



; Allocate output arrays
; ----------------------
sz = SIZE(fit_func)
n_parms = sz(2)

ndf = 2*n_freq - n_parms

dst_parms = FLTARR(4,n_spec)
dst_sig = FLTARR(4,n_spec)
dst_covar = FLTARR(n_spec)

lne_parms = FLTARR(n_parms,n_spec)
chi_sqr = FLTARR(n_spec)
chi_pep = FLTARR(n_spec)



; Generate curvature matrix (fwf) and its inverse (pcvr)
; ------------------------------------------------------
hc = dpc_math(fit_func,code='h')

fw = dpc_math(hc,[[line_wgt],[0*line_wgt]],code='g')
fwf = dpc_math(fw,fit_func,code='m')

pcvr = dpc_math(fwf,code='i')


; Main fitting loop
; -----------------
FOR i=0,n_spec-1 DO BEGIN

	dst_model = FLA_FITDUST(freq,FLOAT(spec(*,i)),dst_wgt,mask=mask,$
			        parm0=parm0,active=active,out=out,$
			        f_char=f_char,parm_sig=p_sig, $
				err_mat=err_mat)


	dst_parms(*,i) = out(0:3)
	dst_sig(*,i) = p_sig(0:3) / SQRT(nifgs(i))
		; Store best-fit CBR and dust parameters and uncertainties

	IF (N_ELEMENTS(err_mat) GT 1) THEN $
		dst_covar(i) = err_mat(1) / nifgs(i)
		; Get tau-temp covariance



	; Generate residual spectra
	; -------------------------
	spc_r = DOUBLE(spec(*,i)) - dst_model

	IF (type NE 6) THEN spc_i = 0 * spc_r ELSE $
			    spc_i = DOUBLE(IMAGINARY(spec(*,i)))
	spc =[[spc_r],[spc_i]]



	; Fit for spectral line fluxes
	; ----------------------------
	fwy = dpc_math(fw,spc,code='m')
	prm = dpc_math(pcvr,fwy,code='m')

	lne_parms(*,i) = prm
		; Store line fluxes




	; Calculate chi-squared
	; ---------------------
	ywy = TOTAL(line_wgt # (spc*spc))

	awa = dpc_math(prm,fwf,prm,code='n')

	chi_sqr(i) = nifgs(i) * (ywy - awa(0))
		; Compute chi-squared


	IF (i/100 EQ i/100.) THEN print,i,chi_sqr(i)/ndf,ndf

ENDFOR




; Generate line-subtracted dust spectra
; -------------------------------------
dust_spec = spec - fit_func(*,0:n_lines-1) # lne_parms(0:n_lines-1,*)



RETURN,chi_sqr
END
