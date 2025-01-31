FUNCTION FMD_pix_sum,pixel=pixel,data=data,cmp_px=cmp_px,wgt=wgt

;
; Returns the weighted pixel average of input array DATA.
;
; Calls  :  PIXAVG
;
; Modified from FDS_PIX_SUM (J.Gales), K.Jensen, HSTX, 23-Apr-97
;
;

h = HISTOGRAM(pixel,min=0)
h = h(WHERE(h GT 0))

sent = -1e30
srt = SORT(pixel)



IF (N_ELEMENTS(wgt) NE 0) THEN BEGIN

	; Weighted Pixel Average
	; ----------------------
	px = LONG(pixel(srt))
	w = wgt(srt)
	PIXAVG,px,w,sent
	w = w * h

	px = LONG(pixel(srt))
	d = data(srt) * wgt(srt)
	PIXAVG,px,d,sent

	d = d * h / w

ENDIF ELSE BEGIN

	; Pixel Sum
	; ---------
	px = LONG(pixel(srt))
	d = data(srt)
	PIXAVG,px,d,sent
	d = d * h

ENDELSE

cmp_px = px

RETURN,d
END
