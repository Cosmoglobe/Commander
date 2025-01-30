;_______________________________________________________________________________
;
;+NAME/ONE-LINE DESCRIPTION:
;    FRD_COMPUTE_CSP computes the merge correction spectra.
;
; DESCRIPTION:
;    FRD_COMPUTE_CSP computes the FIRAS merge offset and gain correction spectra
;    for the low or high frequency channels from the FEF difference and ratio
;    spectra and from the relative weights of the calibration model solutions.
;
; CALLING SEQUENCE:
;    FRD_COMPUTE_CSP, s01_d, s12_d, s23_d, s30_d, $
;                     s01_r, s12_r, s23_r, s30_r, $
;                     cw0, cw1, cw2, cw3, $
;                     off_cor, gain_cor
;
; ARGUMENTS (I= input, O=output)
;    s01_d     I  complex arr  first difference spectrum
;    s12_d     I  complex arr  second difference spectrum
;    s23_d     I  complex arr  third difference spectrum
;    s30_d     I  complex arr  fourth difference spectrum
;    s01_r     I  float arr    first ratio spectrum
;    s12_r     I  float arr    second ratio spectrum
;    s23_r     I  float arr    third ratio spectrum
;    s30_r     I  float arr    fourth ratio spectrum
;    cw0       I  float arr    first cal weight
;    cw1       I  float arr    second cal weight
;    cw2       I  float arr    third cal weight
;    cw3       I  float arr    fourth cal weight
;    off_cor   O  complex arr  offset correction spectra
;    gain_cor  O  float arr    gain correction spectra
;
; WARNINGS:
;    This routine does not compute the corrections for the low frequency long
;    fast scan modes.
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
;    Written by Joel Gales, Applied Research Corp., 25 October 1994
;    Redrafted by Gene Eplee, General Sciences Corp., 27 October 1994
;    Modified by Ken Jensen, Hughes STX, 2 November 1994,  arrays S31_D and
;      S31_R  renamed to  S30_D nad S30_R  for clarity.
;-
;______________________________________________________________________________
;
PRO frd_compute_csp,s01_d,s12_d,s23_d,s30_d, $
                    s01_r,s12_r,s23_r,s30_r, $
                    cw0,cw1,cw2,cw3, $
                    off_cor,gain_cor

;
; Initialize the arrays and matrices.
; -----------------------------------
w0 = DOUBLE([1,cw1/cw0,cw2/cw0,cw3/cw0])

q0 = [[1,-1,0,0], $
      [0,1,-1,0], $
      [0,0,1,-1], $
      [-1,0,0,1]]

n_freq = N_ELEMENTS(s01_d)


;
; Compute the Offset Correction
; -----------------------------
q = [[q0],[w0]]
	; Set up function matrix

mat = q # TRANSPOSE(q)
	; Compute curvature matrix
	
diff = [[s01_d],[s12_d],[s23_d],[s30_d],[COMPLEXARR(n_freq)]]
b = q # TRANSPOSE(diff)

off_cor = TRANSPOSE(INVERT(mat) # b)
	; Compute offset solution


;
; Comput the Gain Correction
; --------------------------
w = w0
a = [0,0,0,0]
	; Initalize solution

gain_cor = DBLARR(4,n_freq)

rat = [[ALOG10(s01_r)],[ALOG10(s12_r)], $
       [ALOG10(s23_r)],[ALOG10(s30_r)], $
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
