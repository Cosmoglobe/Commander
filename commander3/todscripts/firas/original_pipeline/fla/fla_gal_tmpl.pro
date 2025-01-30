PRO FLA_GAL_TMPL


; Build Galactic Template
; -----------------------
restore,'csdr$firas_ftt:high.iss'

w = 1 / (c_high(*,0)^2)
w(0:40) = 0
	; Deweight frequencies below 25 icm

spec = pix2dat(pix=INDGEN(6144),ras=DOUBLE(s_high))

pow = spec # w
pow = 6144 * pow / TOTAL(pow)

pix2xy,INDGEN(6144),data=pow,res=6,/six,ras=gal_tmpl

save,file='csdr$firas_save:gal_tmpl.iss',gal_tmpl

RETURN
END
