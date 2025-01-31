PRO fef_spectra_xx,map1,map2,n_map1,n_map2,cvec1,cvec2,diff=diff,ratio=ratio
;
;
;   Written by Joel Gales, ARC, Nov 1994
;
;   Modified from FEF_SPECTRA to compute difference spectra
;   without computing uncertainties.
;

n_freq = N_ELEMENTS(cvec1)

; Difference Spectrum
; -------------------
pole = BYTARR(6144)
ll = coorconv(indgen(6144),infmt='p',inco='f',outfmt='l',outco='g')
pole(WHERE(ABS(ll(*,1)) GE 60)) = 1
	; Mask out pixels |b| < 60

diff = COMPLEXARR(n_freq)
	; Allocate space for difference

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


FOR i=0,n_freq-1 DO BEGIN

	c1 = cvec1(i) / 1d-6
	c2 = cvec2(i) / 1d-6
		; C - vector in units of 1e-6 eplees (1 & 2)

	wgt = 1 / ((c1*c1/w1) + (c2*c2/w2))
	tot_wgt = TOTAL(wgt)
		; Combined weight

	diff(i) = TOTAL((spec1(*,i) - spec2(*,i)) * wgt) / tot_wgt
		; Difference

ENDFOR


; Ratio Spectrum
; --------------
conv = 0.01
	; Convergence criterion
max_iter = 50
	; Max # of iterations

gain = FLTARR(n_freq)
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


ENDFOR

ratio = FLOAT(gain)

IF (switch EQ 1) THEN ratio = 1 / ratio

	
RETURN
END
