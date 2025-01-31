;_______________________________________________________________________________
;
;+NAME/ONE-LINE DESCRIPTION:
;    FRD_COMPUTE_LR computes the merge correction spectra for LRES.
;
; DESCRIPTION:
;    FRD_COMPUTE_LR computes the FIRAS merge offset correction spectra
;    for the high and lowf skymap merge from the FEF difference spectra.
;    The gain correction spectra are set to zero.
;
; CALLING SEQUENCE:
;    FRD_COMPUTE_LR, s01_d, s01_r, $
;                    off_cor, gain_cor
;
; ARGUMENTS (I= input, O=output)
;    s01_d     I  complex arr  difference spectrum
;    s01_r     I  complex arr  ratio spectrum
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
;       (2)  Set the gain correction spectra to zero.
;
; LIBRARY CALLS:  none
;
; MODIFICATION HISTORY:
;    Written by Gene Eplee, General Sciences Corp., 4 November 1994
;    Renamed FRD_COMPUTE_LR.PRO by Ken Jensen, Hughes STX, 7 Nov 1994.
;-
;______________________________________________________________________________
;
PRO frd_compute_lr,s01_d,s01_r, $
                   off_cor,gain_cor

;
; Initialize the arrays.
; -----------------------------------
nsize = size(s01_d)
off_cor = complexarr(nsize(1),2)
gain_cor = fltarr(nsize(1),2)

;
; Assume that the low frequency spectrum is correct and
;    convert the difference spectrum into a high frequency correction spectrum.
; -----------------------------
for j=0,nsize(1)-1 do off_cor(j,0) = s01_d(j)


RETURN
END
