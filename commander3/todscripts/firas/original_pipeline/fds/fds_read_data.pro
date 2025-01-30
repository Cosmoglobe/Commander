;
;________________________________________________________________________
;
;+NAME/BRIEF DESCRIPTION OF ROUTINE:
;     FDS_READ_DATA reads and concatenates the sky spectra for the destriper.
;
;MODIFICATION HISTORY:
;     Written by Joel Gales, ARC, April 1993
;     Modified by Shirley M. Read, Hughes STX, May 18, 1993 to correct the
;         declaration of pix_spec from DBLARR to COMPLEXARR.
;     Modified by Ken Jensen, Hughes STX, March 1, 1994 , Logical pointer
;         to FCF_SKY records changed to CSDR$FIRAS_INSKY .
;     Modified by Ken Jensen, Hughes STX, April 7, 1994 , Corrected
;     Modified by Ken Jensen, Hughes STX, April 11, 1994 , Row and
;         Column for Cal Glitch Corrections obtained from Model_Label
;         instead of input ID.
;     Modified by Joel Gales, ARC, June 2, 1994, Read dihedral temp
;     Modified by Joel Gales, ARC, Aug  9, 1994, Renormalize cal
;         weights using sky weight renormalization factor (SPR 11866)
;     Modified by Ken Jensen, Hughes STX, Aug 11, 1994 , adjustment to
;         cal weights normalization (SPR 11866)
;-
;____________________________________________________________________________

PRO fds_read_data,arc_file,id,pixel,obs_time,spec,nifgs, $
	          glon,glat,model,f,n_xxx,b_xxx,l_xxx, $
		  cal_spec,cal_nifgs,cal_obs_time,xcal,ical,refh,skyh, $
		  dihd,galcut=galcut,uxxx=uxxx,weight_cor=weight_cor, $
	          sky_glitch=sky_glitch,cal_glitch=cal_glitch, $
		  sky_wgts=sky_wgts,cal_wgts=cal_wgts,v_xxx=v_xxx, $
		  sky_lbl=sky_lbl,cal_lbl=cal_lbl


eject = timeconv('893251118',infmt='z',outfmt='s')
		; get cover eject time in seconds


r2d = 180 / !dpi

IF (N_ELEMENTS(weight_cor) EQ 0) THEN  weight_cor = 1
print, 'Setting IFG weight to', strcompress( weight_cor )




; Open FCF file list
; ------------------
print,'Reading FCF archive times listing'
file = STRARR(40)
n_files = 0
a = ''

arc_file = [arc_file]
id = [id]
	; make array if scalar

FOR k=0,N_ELEMENTS(arc_file)-1 DO BEGIN

	OPENR,1,'csdr$firas_ref:'+arc_file(k)+'.txt'

	WHILE (NOT EOF(1)) DO BEGIN
		READF,1,a
		file(n_files) = 'fcf_sky_' + id(k) + '.' + a
		n_files = n_files + 1
	ENDWHILE
	CLOSE,1

ENDFOR

file = file(0:n_files-1)
file = STRUPCASE(file)





; Determine calibration model and nyquist frequency
; -------------------------------------------------
print,'Reading calibration model solution and nyquist frequency'
fld = 'spec_data.model_label,coad_spec_data.nyquist_icm'
stat = read_skymap('csdr$firas_insky:'+file(0),fld,pix,label,nyq)

model = label(0)
nyquist = nyq(0)



; Build frequency vector
; ---------------------
IF (STRUPCASE(STRMID(model,1,1)) EQ 'H') THEN BEGIN
	bins = '5:171'
	n_freq = 167
	f = (4 + FINDGEN(167))*nyquist/256
ENDIF ELSE IF (STRUPCASE(STRMID(model,1,2)) EQ 'LS') THEN BEGIN
	bins = '5:38'
	n_freq = 34
	f = (4 + FINDGEN(34))*nyquist/256
ENDIF ELSE IF (STRUPCASE(STRMID(model,2,2)) EQ 'FL') THEN BEGIN
	bins = '5:38'
	n_freq = 34
	f = (4 + FINDGEN(34))*nyquist/256
ENDIF ELSE IF (STRUPCASE(STRMID(model,1,2)) EQ 'LL') THEN BEGIN
	bins = '9:156'
	n_freq = 148
	f = (8 + FINDGEN(148))*nyquist/256
ENDIF

model = STRMID(label(0),5,100)
		; trim off scan mode info 




; Read glitch correction parameters
; ---------------------------------
sl_int = DBLARR(24)
OPENR,1,'csdr$firas_ref:fex_gltchcor.dat'
READU,1,sl_int
CLOSE,1
slpe = REFORM(sl_int(0:11),3,4)
intc = REFORM(sl_int(12:*),3,4)
	; Read slope and intercept of glitch tweak for all scan modes





; Read Data from FCF skymap files
; -------------------------------
flds = 'ct_head.gmt'
flds = flds + ',spec_data.spec(' + bins + ')'
flds = flds + ',coad_spec_head.num_ifgs,attitude.galactic_longitude'
flds = flds + ',attitude.galactic_latitude'
flds = flds + ',spec_data.model_label'
flds = flds + ',coad_spec_data.glitch_rate'
FOR j=0,n_files-1 DO BEGIN

	PRINT,'Reading ',file(j)
	stat = read_skymap('csdr$firas_insky:'+file(j),flds,pix,gmt,$
			    s,n,lon,lat,label,gl)

	tme = timeconv(gmt,infmt='z',outfmt='s')

	IF (j EQ 0) THEN BEGIN

		pixel = pix
		spec = s

		obs_time = (tme - eject) / 86400.

		sky_lbl = STRMID(label,0,4)
		nifgs = n

		row = STRPOS('RH_RL_LH_LL',STRMID(label(0),0,2))/3
		col = STRPOS('SS_SF_LF_FL',STRMID(label(0),2,2))/3
		col = col < 2

		sl = slpe(col,row)
		intcpt = intc(col,row)
		sky_wgts = n / (sl*gl+intcpt)

		glon = lon*r2d*1.e-4
		glat = lat*r2d*1.e-4
		sky_glitch = gl

	ENDIF ELSE BEGIN

		pixel = [pixel,pix]
		spec = [[spec],[s]]

		obs_time = [obs_time,(tme - eject) / 86400.]

		nifgs = [nifgs,n]
		sky_lbl = [sky_lbl,STRMID(label,0,4)]

		row = STRPOS('RH_RL_LH_LL',STRMID(label(0),0,2))/3
		col = STRPOS('SS_SF_LF_FL',STRMID(label(0),2,2))/3
		col = col < 2

		sl = slpe(col,row)
		intcpt = intc(col,row)
		sky_wgts = [sky_wgts,n/(sl*gl+intcpt)]

		glon = [glon,lon*r2d*1.e-4]
		glat = [glat,lat*r2d*1.e-4]
		sky_glitch = [sky_glitch,gl]
	ENDELSE

ENDFOR

obs_time = FLOAT(obs_time)


i = WHERE(STRMID(sky_lbl,2,2) EQ 'SS')
IF (i(0) NE -1) THEN BEGIN
	fac_ss = TOTAL(nifgs(i)) / TOTAL(sky_wgts(i))
	sky_wgts(i) = sky_wgts(i) * fac_ss
ENDIF

i = WHERE(STRMID(sky_lbl,2,2) EQ 'SF')
IF (i(0) NE -1) THEN BEGIN
	fac_sf = TOTAL(nifgs(i)) / TOTAL(sky_wgts(i))
	sky_wgts(i) = sky_wgts(i) * fac_sf
ENDIF

i = WHERE(STRMID(sky_lbl,2,2) EQ 'LF')
IF (i(0) NE -1) THEN BEGIN
	fac_lf = TOTAL(nifgs(i)) / TOTAL(sky_wgts(i))
	sky_wgts(i) = sky_wgts(i) * fac_lf
ENDIF

i = WHERE(STRMID(sky_lbl,2,2) EQ 'FL')
IF (i(0) NE -1) THEN BEGIN
	fac_fl = TOTAL(nifgs(i)) / TOTAL(sky_wgts(i))
	sky_wgts(i) = sky_wgts(i) * fac_fl
ENDIF
	; Renormalize weights




; Reweight weights for scan mode
; ------------------------------
FOR i=0,N_ELEMENTS(sky_wgts)-1 DO BEGIN
	sc_mode = STRPOS('SF_LF_FL',STRMID(sky_lbl(i),2,2))/3
	sc_mode = sc_mode < 1
	sky_wgts(i) = sky_wgts(i) * (weight_cor ^ (0.5 - sc_mode))
ENDFOR




; Sort by pixel number
; --------------------
print,'Sorting by pixel number'
i = SORT(pixel)

pixel = pixel(i)
spec = spec(*,i)
obs_time = obs_time(i)
nifgs = nifgs(i)
sky_wgts = sky_wgts(i)
sky_glitch = sky_glitch(i)
sky_lbl = sky_lbl(i)
glon = glon(i)
glat = glat(i)




; Generate sixpack of weights
; ---------------------------

pix_dist = WHERE(HISTOGRAM(MIN=0,pixel) GT 0)
n_dist = N_ELEMENTS(pix_dist)
		; get # of distinct pixels

num_pix = FLTARR(n_dist)
j = 0
num_pix(0) = sky_wgts(0)

FOR k=1,N_ELEMENTS(pixel)-1 DO BEGIN	; loop over all obs
	IF (pixel(k) NE pixel(k-1)) THEN j = j + 1	; new pixel
	num_pix(j) = num_pix(j) + sky_wgts(k)
ENDFOR

print, 'Total # of IFGs' + Strcompress( total(sky_wgts) )
print, 'Total from num_pix' + strcompress( total(num_pix) )

pix2xy,/six,pix_dist,data=num_pix,res=6,raster=n_xxx




; Generate sixpack of weighted galactic longitude and latitude
; ------------------------------------------------------------
print,'b_xxx'
uv = coorconv([[glon],[glat]],infmt='l',outfmt='u')

acm_x = FLTARR(n_dist)
acm_y = acm_x
acm_z = acm_x
j = 0
acm_x(0) = uv(0,0)*sky_wgts(0)
acm_y(0) = uv(0,1)*sky_wgts(0)
acm_z(0) = uv(0,2)*sky_wgts(0)

FOR k=1,N_ELEMENTS(pixel)-1 DO BEGIN	; loop over all obs
	IF (pixel(k) NE pixel(k-1)) THEN j = j + 1	; new pixel
	acm_x(j) = acm_x(j) + uv(k,0)*sky_wgts(k)
	acm_y(j) = acm_y(j) + uv(k,1)*sky_wgts(k)
	acm_z(j) = acm_z(j) + uv(k,2)*sky_wgts(k)
ENDFOR

acm_x = acm_x / num_pix
acm_y = acm_y / num_pix
acm_z = acm_z / num_pix

temp = coorconv([[acm_x],[acm_y],[acm_z]],infmt='u',outfmt='l')

pix2xy,/six,pix_dist,data=temp(*,0),res=6,raster=l_xxx
pix2xy,/six,pix_dist,data=temp(*,1),res=6,raster=b_xxx




; Calculating uncorrected skymap
; ------------------------------
max_pix = MAX(pixel)
n_obs = N_ELEMENTS(pixel)


print,'Determining spatial (pixel) parameters'
diag = DBLARR(max_pix+1)
pix_spec = COMPLEXARR(n_freq,max_pix+1)
hist = HISTOGRAM(MIN=0,pixel)

last_pix = 0
FOR i=0,max_pix DO BEGIN
	IF (hist(i) NE 0) THEN BEGIN
		j = WHERE(pixel(last_pix:last_pix+hist(i)-1) EQ i)
		j = j + last_pix
		spt = spec(*,j)
		nm = sky_wgts(j)
		diag(i) = TOTAL(nm)
		pix_spec(*,i) = spt # nm
	ENDIF
last_pix = last_pix + hist(i)
ENDFOR

nz = WHERE(diag NE 0)

afp = pix_spec 
FOR i=0,n_freq-1 DO pix_spec(i,nz) = pix_spec(i,nz) / diag(nz)

pix2xy,/six,nz,data=pix_spec(*,nz),res=6,raster=v_xxx




; Calculating uncorrected sigma
; -----------------------------
print,'Calculating sigma(f)'

rr = FLTARR(n_freq)
ri = FLTARR(n_freq)

mask = 0 * sky_wgts
mask(WHERE(ABS(glat) GE galcut)) = 1

FOR j=0,n_obs-1 DO begin
	rr = rr + (FLOAT(spec(*,j)-pix_spec(*,pixel(j))))^2 * sky_wgts(j) * mask(j)
	ri = ri + (IMAGINARY(spec(*,j)-pix_spec(*,pixel(j))))^2 * sky_wgts(j) * mask(j)
ENDFOR

pix_nogal = pixel(WHERE(mask EQ 1))
n_dist = N_ELEMENTS(WHERE(HISTOGRAM(pix_nogal) GT 0))
ndf = (TOTAL(mask) - n_dist)
uxxx = COMPLEX(SQRT(rr/ndf),SQRT(ri/ndf))





; Read Data from FCF cal files
; ----------------------------
PRINT,'Reading Cal Data'

flds = 'spec_data.spec(' + bins + ')'
flds = flds + ',ct_head.gmt,coad_spec_head.num_ifgs'
flds = flds + ',coad_spec_data.xcal,coad_spec_data.ical'
flds = flds + ',coad_spec_data.skyhorn,coad_spec_data.refhorn'
flds = flds + ',coad_spec_data.glitch_rate'

flds2 = 'spec_data.model_label,coad_spec_data.dihedral'

FOR j=0,N_ELEMENTS(id)-1 DO BEGIN

	cfile = 'csdr$firas_incal:fcf_cal_' + id(j)

	status = read_tod(cfile,flds,'89001','90365',s,gmt,n,$
			  xc,ic,sh,rh,gl,maxrec=500)

	status = read_tod(cfile,flds2,'89001','90365',label,dh,maxrec=500)


	IF (j EQ 0) THEN BEGIN

		tme = timeconv(gmt,infmt='z',outfmt='s')
		cal_obs_time = (tme - eject) / 86400.

		cal_spec = s
		cal_nifgs = n

		xcal = xc
		ical = ic
		skyh = sh
		refh = rh
		dihd = dh

		row = STRPOS('RH_RL_LH_LL',STRMID(label(0),0,2))/3
		col = STRPOS('SS_SF_LF_FL',STRMID(label(0),2,2))/3
		col = col < 2

		sl = slpe(col,row)
		intcpt = intc(col,row)
		cw = n / (sl*gl+intcpt)

		cal_wgts = cw
		cal_glitch = gl

		cal_lbl = STRMID(label,0,4)

	ENDIF ELSE BEGIN

		tme = timeconv(gmt,infmt='z',outfmt='s')
		cal_obs_time = [cal_obs_time,(tme - eject) / 86400.]

		cal_spec = [[cal_spec],[s]]
		cal_nifgs = [cal_nifgs,n]

		xcal = [xcal,xc]
		ical = [ical,ic]
		skyh = [skyh,sh]
		refh = [refh,rh]
		dihd = [dihd,dh]

		row = STRPOS('RH_RL_LH_LL',STRMID(label(0),0,2))/3
		col = STRPOS('SS_SF_LF_FL',STRMID(label(0),2,2))/3
		col = col < 2

		sl = slpe(col,row)
		intcpt = intc(col,row)
		cw = n / (sl*gl+intcpt)

		cal_wgts = [cal_wgts,cw]

		cal_glitch = [cal_glitch,gl]

		cal_lbl = [cal_lbl,STRMID(label,0,4)]

	ENDELSE


ENDFOR


; Disabling the Temperature-dependent Re-Normalization (SPR 11866)
;cal_wgts = cal_wgts * SQRT(22)*2.73*2.73 / $
;		SQRT(10*ical^4 + 10*xcal^4 + skyh^4 + refh^4)
;



; Renormalize cal weights by sky weight factor
; --------------------------------------------
i = WHERE(STRMID(cal_lbl,2,2) EQ 'SS')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_ss

i = WHERE(STRMID(cal_lbl,2,2) EQ 'SF')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_sf

i = WHERE(STRMID(cal_lbl,2,2) EQ 'LF')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_lf

i = WHERE(STRMID(cal_lbl,2,2) EQ 'FL')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_fl





FOR i=0,N_ELEMENTS(cal_wgts)-1 DO BEGIN
	sc_mode = STRPOS('SF_LF_FL',STRMID(cal_lbl(i),2,2))/3
	sc_mode = sc_mode < 1
	cal_wgts(i) = cal_wgts(i) * (weight_cor ^ (0.5 - sc_mode))
ENDFOR


RETURN
END
