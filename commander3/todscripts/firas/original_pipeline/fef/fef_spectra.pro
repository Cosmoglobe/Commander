PRO fef_spectra,f,map1,map2,n_map1,n_map2,cvec1,cvec2, $
		pcvr1=pcvr1,diag1=diag1,rect1=rect1, $
		pcvr2=pcvr2,diag2=diag2,rect2=rect2, $
		pep1=pep1,pep2=pep2, $
                diff=diff,d_unc=d_unc,ratio=ratio,r_unc=r_unc

;
;
;   Written by Joel Gales, ARC, Sept 1994
;   Modification by Ken Jensen, Hughes STX, 18 Sept 1994, Correction in
;    calculation of Ratio Weight for second skymap (N2).
;   Modification by Joel Gales, ARC, 21 Sept 1994, Use common weight for
;    spectra 1 and 2
;   Modification by Joel Gales, ARC, 23 Sept 1994, Compute ratio spectra
;    using least-squared fit gain model.
;   Modification by Ken Jensen, Hughes STX, 23 Sept 1994, The variable
;    DIFF is replaced by DIFFX in the code which calculates Ratio Spectra.
;    Necessary to preserve the Difference Spectra.
;   Modification by Joel Gales, ARC, 26 Sept 1994, (SPR 11918)
;	1) Modify convergence criterion to use delta(chi^2) rather than
;	   delta(gain).
;	2) Use destriper errors (beta) to compute difference uncertainty.
;	3) Use PEP gain errors in ratio uncertainty.
;   Modification by Ken Jensen, Hughes STX, 27 Sept 1994, N_STRIPES is
;       defined for EACH skymap, since they can have a different number
;       of stripes.
;   Modification by Joel Gales, ARC, 27 Sept 1994, (SPR 11921)
;        Switch datasets, if necessary, for ratio spectrum so that
;        "noisy" is 1 and "clean" is 2.

n_freq = N_ELEMENTS(f)


; Difference Spectrum
; -------------------
pole = BYTARR(6144)
ll = coorconv(indgen(6144),infmt='p',inco='f',outfmt='l',outco='g')
pole(WHERE(ABS(ll(*,1)) GE 60)) = 1
	; Mask out pixels |b| < 60

diff = COMPLEXARR(n_freq)
d_unc  = FLTARR(n_freq)
	; Allocate space for difference and difference uncertainty

w1 = pix2dat(pix=INDGEN(6144),ras=n_map1)
px1 = WHERE(w1 NE 0)
w1 = w1 * pole + 1d-12

w2 = pix2dat(pix=INDGEN(6144),ras=n_map2)
px2 = WHERE(w2 NE 0)
w2 = w2 * pole + 1d-12
	; Individual weights (Add small number to avoid divide by 0)

spec1 = pix2dat(pix=INDGEN(6144),ras=map1)
spec2 = pix2dat(pix=INDGEN(6144),ras=map2)
	; Spectra (1 & 2)



; Beta
; ----
sz = SIZE(pcvr1)
n_stripes = sz(1)

one = INTARR(n_stripes) + 1

d_1 = diag1 # one
beta1 = (rect1 / d_1) # fds_choleskey(pcvr1)

sz = SIZE(pcvr2)
n_stripes = sz(1)

one = INTARR(n_stripes) + 1

d_2 = diag2 # one
beta2 = (rect2 / d_2) # fds_choleskey(pcvr2)



FOR i=0,n_freq-1 DO BEGIN

	c1 = cvec1(i) / 1d-6
	c2 = cvec2(i) / 1d-6
		; C - vector in units of 1e-6 eplees (1 & 2)

	b1 = beta1 * c1
	b2 = beta2 * c2
		; Destriper errors (beta)

	wgt = 1 / ((c1*c1/w1) + (c2*c2/w2))
	tot_wgt = TOTAL(wgt)
		; Combined weight

	diff(i) = TOTAL((spec1(*,i) - spec2(*,i)) * wgt) / tot_wgt
		; Difference

	d_unc(i) = 1d-6 * SQRT((1 / tot_wgt) + $
			       (TOTAL((wgt(px1) # b1)^2) + $
			        TOTAL((wgt(px2) # b2)^2)) / tot_wgt^2)
		; Uncertainty

ENDFOR



; Ratio Spectrum
; --------------
conv = 0.01
	; Convergence criterion
max_iter = 50
	; Max # of iterations

gain = FLTARR(n_freq)
r_unc = FLTARR(n_freq)
	; Allocate space for output

n_eff1 = TOTAL(1/cvec1^2)
n_eff2 = TOTAL(1/cvec2^2)

IF (n_eff1 GT n_eff2) THEN BEGIN
	PRINT,'Datasets switched for ratio spectrum'
	switch = 1
ENDIF ELSE switch = 0


; Effective N's
; -------------
n1 = pix2dat(pix=INDGEN(6144),ras=n_map1)
n2 = pix2dat(pix=INDGEN(6144),ras=n_map2)

nz = WHERE(n1 NE 0 OR n2 NE 0)
	; Non-Empty pixels

n1 = n1(nz)
n2 = n2(nz)

spec1 = spec1(nz,*)
spec2 = spec2(nz,*)
	; Get spectra from these pixels


FOR i=0,n_freq-1 DO BEGIN

	gn = 1
	diff_chi2 = conv + 1
	n_iter = 0
	chi2 = 0

	w1 = n1 / cvec1(i)^2
	w2 = n2 / cvec2(i)^2

	s1 = DOUBLE(spec1(*,i))
	s2 = DOUBLE(spec2(*,i))

	IF (switch EQ 1) THEN BEGIN

		tmp = w2
		w2 = w1
		w1 = tmp

		tmp = s2
		s2 = s1
		s1 = tmp

	ENDIF

	sw1 = s1 * w1
	sw2 = s2 * w2

	WHILE(ABS(diff_chi2) GT conv AND n_iter LE max_iter) DO BEGIN

		old_chi2 = chi2

		m_pp = gn ^2 * w1 + w2
		s_f = gn * sw1 + sw2
		a_p = s_f/ m_pp

		m_gg = TOTAL(a_p^2 * w1)
		gn = TOTAL(sw1 * a_p) / m_gg

		chi2 = gn^2 * m_gg + TOTAL(a_p^2 * w2) - 2 * TOTAL(s_f * a_p)
		diff_chi2 = chi2 - old_chi2
		n_iter = n_iter + 1

	ENDWHILE

	gain(i) = gn

	m_gp = gn * a_p * w1
	r_var0 = 1 / (m_gg - TOTAL(m_gp^2/m_pp))

	IF (switch EQ 1) THEN r_var0 = r_var0 / gn^4

	r_unc(i) = SQRT(r_var0 + pep1(i)^2 + pep2(i)^2)
		; Uncertainty in Ratio		   

ENDFOR

ratio = FLOAT(gain)

IF (switch EQ 1) THEN ratio = 1 / ratio

	
RETURN
END
