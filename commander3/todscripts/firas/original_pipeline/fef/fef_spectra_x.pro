PRO fef_spectra_x,map1,map2,n_map1,n_map2,cvec1,cvec2,diff=diff
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

	
RETURN
END
