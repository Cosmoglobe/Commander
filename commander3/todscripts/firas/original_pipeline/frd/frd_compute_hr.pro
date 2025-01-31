;_______________________________________________________________________________
;
;+NAME/ONE-LINE DESCRIPTION:
;    FRD_COMPUTE_HR computes the merge correction spectra for HRES.
;
; DESCRIPTION:
;    FRD_COMPUTE_HR computes the FIRAS merge offset and gain correction
;    spectra for the low frequency long fast skymap merge from the FEF
;    difference and ratio spectra and from the relative weights of the
;    calibration model solutions.
;
; CALLING SEQUENCE:
;    FRD_COMPUTE_HR, s01_d, s01_r, $
;                    cw0, cw1, $
;                    off_cor, gain_cor
;
; ARGUMENTS (I= input, O=output)
;    s01_d     I  complex arr  difference spectrum
;    s01_r     I  complex arr  ratio spectrum
;    cw0       I  float arr    first cal weight
;    cw1       I  float arr    second cal weight
;    off_cor   O  complex arr  offset correction spectra
;    gain_cor  O  float arr    gain correction spectra
;
; WARNINGS:  none
;
;#
; COMMON BLOCKS:  none
;
; PROCEDURE:
;       (1)  Compute the offset correction spectra from the difference spectra
;       (2)  Iteratively compute the gain correction spectra from the ratio
;              spectra
;
; LIBRARY CALLS:  none
;
; MODIFICATION HISTORY:
;    Written by Gene Eplee, General Sciences Corp., 27 October 1994
;    Modified by Gene Eplee, General Sciences Corp., 4 November 1994
;        to work with the HRES combination only.  Renamed from
;        FRD_COMPUTE_RES.PRO.
;    Renamed FRD_COMPUTE_HR.PRO by Ken Jensen, Hughes STX, 7 Nov 1994.
;-
;______________________________________________________________________________
;
PRO frd_compute_hr,s01_d,s01_r, $
                   cw0,cw1, $
                   off_cor,gain_cor

;
; Initialize the arrays and matrices.
; -----------------------------------
w0 = DOUBLE([1,cw1/cw0])

q0 = [[1,-1], $
      [-1,1]]

n_freq = N_ELEMENTS(s01_d)


;
; Compute the Offset Correction
; -----------------------------
q = [[q0],[w0]]
	; Set up function matrix

mat = q # TRANSPOSE(q)
	; Compute curvature matrix
	
diff = [[s01_d],[-s01_d],[COMPLEXARR(n_freq)]]
b = q # TRANSPOSE(diff)

off_cor = TRANSPOSE(INVERT(mat) # b)
	; Compute offset solution


;
; Comput the Gain Correction
; --------------------------
w = w0
a = [0,0]
	; Initalize solution

gain_cor = DBLARR(2,n_freq)

rat = [[ALOG10(s01_r)],[ALOG10(1.0/s01_r)], $
       [DBLARR(n_freq)]]

FOR i=0,n_freq-1 DO BEGIN

	REPEAT BEGIN

		last_a = a
		q = [[q0],[w]]
		mat = q # TRANSPOSE(q)
		b = q # TRANSPOSE(rat(i,*))
		a = INVERT(mat) # b

		gn = 10 ^ a
		gn0 = gn / gn(0)
		w = w0 / gn0^2

	ENDREP UNTIL (TOTAL(ABS(a - last_a)) LE 1E-5)

	gain_cor(*,i) = gn

ENDFOR

gain_cor = FLOAT(TRANSPOSE(gain_cor))


RETURN
END
