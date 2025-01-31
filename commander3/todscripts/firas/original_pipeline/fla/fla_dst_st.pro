PRO fla_dst_st,b_xxx=b_xxx,l_xxx=l_xxx,n_xxx=n_xxx,spec=spec, $
	    chanscan=chanscan

;
; This program writes the FLA_DST skymap records.
;
; The 'b_xxx', 'l_xxx', and 'n_xxx' arrays contain the galactic latitude,
; galactic longitude, and number of IFGs within each pixel and are
; stored in the xxx.ISS raw data saveset.
;
; The 'covar' array contains the line fit parameter covariance matrix.
; The 'line_parms' and 'dst_parms' arrays contain the line fluxes and
; baseline paramaters and the continuum dust parameters.  These
; variables are stored in the 'lines_xxx_xxxx.ISS' saveset.
;
;
;  Modification History:
;	Written by Joel Gales, ARC, October 1994 as FLA_DST.PRO
;	Modified by Ken Jensen, Hughes STX, Oct 24, 1994, and renamed
;	    FLA_DST_ST.PRO.
;	Modified by Fred Shuman, Hughes STX, 1994 Dec 16.  Corrected array
;	    index typos in final storage:
;From:	dst_struct.eclat = ll_e(i,0)	to:	dst_struct.eclat = ll_e(i,1)
;From:	dst_struct.dec = ll_q(i,0)	to:	dst_struct.dec = ll_q(i,1)



; Determine Non-zero Pixels
; -------------------------
n = pix2dat(pix=INDGEN(6144),raster=n_xxx)
pixel = WHERE(n GT 0)
nifg_in_pix = n(pixel)
	; get non-zero pixels

ll = coorconv(pixel,infmt='p',inco='f',outfmt='l',outco='g')
n_pix = N_ELEMENTS(pixel)
	; get non-zero pixels within 'gal_cut' degrees of gal plane




; Extract Average Longitude and Latitude
; --------------------------------------
glon = pix2dat(pix=pixel,raster=l_xxx)
glat = pix2dat(pix=pixel,raster=b_xxx)
	; get pixel average galactic lon/lat


ll_e = coorconv([[glon],[glat]],infmt='l',inco='g',outfmt='l',outco='e')
	; convert to ecliptic coordinates

ll_q = coorconv([[glon],[glat]],infmt='l',inco='g',outfmt='l',outco='q')
	; convert to equatorial coordinates





; Define FLA_DST structure
; ------------------------
dst_struct = {PIXEL: LONARR(1), $
	      ECLON: FLTARR(1),$
	      ECLAT: FLTARR(1),$
	      REAL_SPECTRUM: FLTARR(167),$
	      IMAG_SPECTRUM: FLTARR(167),$
	      CHANSCAN: STRING(' ',FORMAT='(A4)'), $
	      NUM_IFGS: FLTARR(1),$
	      GALON: FLTARR(1),$
	      GALAT: FLTARR(1),$
	      RA: FLTARR(1),$
	      DEC: FLTARR(1)}

rec_len = (8 + 2*167) * 4 + 4






; Build FLA_DST filename
; ----------------------
filename = 'csdr$firas_out:fla_dst_' + chanscan + '.f16_93hybrid'




; Write records to FLA_DST_xxxx.F16_93HYBRID
; ------------------------------------------
OPENW, 10, filename, rec_len, /FIXED


FOR i=0,n_pix-1 DO BEGIN

	IF (i/100 EQ i/100.) THEN PRINT,'Writing Record',i

	dst_struct.pixel = pixel(i)
	dst_struct.eclon = ll_e(i,0)
	dst_struct.eclat = ll_e(i,1)

	dst_struct.real_spectrum = FLOAT(spec(*,i))
	dst_struct.imag_spectrum = IMAGINARY(spec(*,i))

	dst_struct.chanscan = STRING(STRUPCASE(chanscan),FORMAT='(A4)')
	dst_struct.num_ifgs = nifg_in_pix(i)

	dst_struct.galon = glon(i)
	dst_struct.galat = glat(i)

	dst_struct.ra =  ll_q(i,0)
	dst_struct.dec = ll_q(i,1)

	WRITEU, 10, dst_struct

ENDFOR


CLOSE, 10

PRINT,' '
PRINT,n_pix,' records written of',rec_len,' bytes in length'

RETURN
END
