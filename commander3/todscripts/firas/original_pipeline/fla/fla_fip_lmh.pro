PRO FLA_FIP_LMH,b_xxx=b_xxx,l_xxx=l_xxx,n_xxx=n_xxx,covar=covar,galcut=galcut, $
 line_parms=line_parms,dst_parms=dst_parms,chanscan=chanscan,dst_sig=dst_sig, $
 dst_covar=dst_covar,n_lines=n_lines,line_idx=line_idx

;
; This program writes the FIP_LMP skymap records.
;
; The 'b_xxx', 'l_xxx', and 'n_xxx' arrays contain the galactic latitude,
; galactic longitude, and number of IFGs within each pixel and are
; stored in the xxx.ISS raw data saveset.
;
; The 'covar' array contains the line fit parameter covariance matrix.
; The 'line_parms' and 'dst_parms' arrays contain the line fluxes and
; baseline paramaters and the continuum dust parameters.  These
; variables are stored in the 'lines_xxx_xxxx.ISS' saveset.
; The 'dst_covar' array contains the dust temperature-tau covariance.
;
;
;  Modification History:
;	Written by Joel Gales, ARC, October 1994 as FLA_FIP_LMP
;       Modified by Ken Jensen, HSTX, 22-Feb-1995, and renamed FLA_FIP_LMH,
;       customized to write the FIP_LMH_xxxx binary file.
;

IF ( N_ELEMENTS(line_idx) gt n_lines) THEN BEGIN
 print,'LINE_IDX Size Error'
 RETURN
ENDIF
;

m2ep = 2.99792458d-7	; megajanskies to eplees
ep2m = 1 / m2ep		; eplees to megajanskies


; Determine Non-zero Galactic Pixels
; ----------------------------------
n = pix2dat(pix=INDGEN(6144),raster=n_xxx)
pixel = WHERE(n GT 0)
nifg_in_pix = n(pixel)
	; get non-zero pixels

IF (KEYWORD_SET(galcut) EQ 0) THEN galcut = 90
ll = coorconv(pixel,infmt='p',inco='f',outfmt='l',outco='g')
gal = WHERE(abs(ll(*,1)) lt galcut)
pix = pixel(gal)
nifg_in_pix = nifg_in_pix(gal)
n_pix = N_ELEMENTS(pix)
	; get non-zero pixels within 'gal_cut' degrees of gal plane




; Extract Average Longitude and Latitude
; --------------------------------------
glon = pix2dat(pix=pix,raster=l_xxx)
glat = pix2dat(pix=pix,raster=b_xxx)
	; get pixel average galactic lon/lat

ll_e = coorconv([[glon],[glat]],infmt='l',inco='g',outfmt='l',outco='e')
	; convert to ecliptic coordinates

ll_q = coorconv([[glon],[glat]],infmt='l',inco='g',outfmt='l',outco='q')
	; convert to equatorial coordinates


; Extract variances from covariance matrix
; ----------------------------------------
sz = SIZE(covar)
n_base = (sz(1) - n_lines) - 1
var = FLTARR(sz(1))
FOR i=0,sz(1)-1 DO var(i) = covar(i,i)
;

; Extract line fluxes
; -------------------
lflux = line_parms(8:15,gal)
;

; Extract dust continuum parameters
; ---------------------------------
c_parms = dst_parms(0:3,gal)
c_sig = dst_sig(0:3,gal)
c_covar = dst_covar(gal)
;

; Define FIP_LMP structure
; ------------------------
lmp_struct = {PIXEL: LONARR(1), $
	      ECLON: FLTARR(1),$
	      ECLAT: FLTARR(1),$
	      LINE_FLUX: FLTARR(8),$
	      LINE_FLUX_SIGMA: FLTARR(8),$
	      DUST_TEMP: FLTARR(1),$
	      DUST_TEMP_SIGMA: FLTARR(1),$
	      DUST_TAU: FLTARR(1),$
	      DUST_TAU_SIGMA: FLTARR(1),$
              DUST_COVAR: FLTARR(1),$
              DUST_INDEX: FLTARR(1),$
              DUST_INDEX_SIG: FLTARR(1),$
	      CHANSCAN: STRING(' ',FORMAT='(A4)'), $
	      NUM_IFGS: FLTARR(1),$
	      GALON: FLTARR(1),$
	      GALAT: FLTARR(1),$
	      RA: FLTARR(1),$
	      DEC: FLTARR(1)}
;
rec_len = ( 15 + 2*8 ) * 4 + 4
;


; Build FIP_LMP filename
; ----------------------
filename = 'csdr$firas_out:fip_lmh_' + chanscan + '.f16_93hybrid'
;

; Write records to FIP_LMy_xxxx.F16_93HYBRID
; ------------------------------------------
OPENW, 10, filename, rec_len, /FIXED

nz=WHERE(line_idx eq 0,cz)

FOR i=0,n_pix-1 DO BEGIN

	IF (i/500 EQ i/500.) THEN PRINT,'Writing Record',i

	lmp_struct.pixel = pix(i)
	lmp_struct.eclon = ll_e(i,0)
	lmp_struct.eclat = ll_e(i,1)

	lf = FLOAT(lflux(*,i))
        IF (cz gt 0) THEN lf(nz) = 0.
	lmp_struct.line_flux = lf * 1e-3

	lfsig = SQRT(var(8:15)/nifg_in_pix(i))
        IF (cz gt 0) THEN lfsig(nz) = 0.
	lmp_struct.line_flux_sigma = lfsig * 1e-3

        IF ((c_parms(1,i) le 10.)or(c_parms(2,i) le 0.)) THEN BEGIN
 	 lmp_struct.dust_temp = 0.
	 lmp_struct.dust_temp_sigma = 0.
 	 lmp_struct.dust_tau = 0.
	 lmp_struct.dust_tau_sigma = 0.
	 lmp_struct.dust_covar = 0.
        ENDIF

        IF ((c_parms(1,i) gt 10.)and(c_parms(2,i) gt 0.)) THEN BEGIN
 	 lmp_struct.dust_temp = c_parms(1,i)
	 lmp_struct.dust_temp_sigma = c_sig(1,i)
 	 lmp_struct.dust_tau = c_parms(2,i)
	 lmp_struct.dust_tau_sigma = c_sig(2,i)
	 lmp_struct.dust_covar = c_covar(i)
        ENDIF

        lmp_struct.dust_index = c_parms(3,i)
        lmp_struct.dust_index_sig = c_sig(3,i)

	lmp_struct.chanscan = STRING(STRUPCASE(chanscan),FORMAT='(A4)')
	lmp_struct.num_ifgs = nifg_in_pix(i)

	lmp_struct.galon = glon(i)
	lmp_struct.galat = glat(i)

	lmp_struct.ra =  ll_q(i,0)
	lmp_struct.dec = ll_q(i,1)


	WRITEU, 10, lmp_struct

ENDFOR


CLOSE, 10

PRINT,' '
PRINT,n_pix,' records written of',rec_len,' bytes in length'

RETURN
END
