;________________________________________________________________________
;
;+NAME/BRIEF DESCRIPTION OF ROUTINE:
;     FMD_READ_DATA reads and concatenates the sky spectra and other
;                   coadd fields for the destriper.
;
;CALLS  :  FMD_ADT_SEC
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
;     Modified by Joel Gales, ARC, Dec 1995
;         Read FSL records instead of FCF records
;	  Replace calls to READ_SKY_MAP by RSM (C linkimage)
;         Read detector responsivity field
;	  Deweight "bad" responsivity coadds for xLSF & xLFL
;         Read Dihedral-binned cal records
;     Modified by Alice Trenholme, FSC, Jan 1996
;         Spectrum length is now passed in
;         Read bolometer voltage
;         Deweight bad coadds for high channels
;     Renamed FMD_READ_DATA,  K.Jensen, Hughes STX, 27-Feb-96
;     Replace CAL_BOL_VOLT, SKY_BOL_VOLT with SCAN, TIME, K.Jensen, 07-Jan-97
;     Calls FMD_ADT_SEC, K.Jensen, 23-Apr-97
;     New Sorting for compatibility with FFP_CSK files,  K.Jensen, 27-May-97.
;     Reads and sorts long and short data separately prior to
;         concatenation, K.Jensen, 27-May-97 .
;     Keyword FSL_IDX added,  K.Jensen, 27-May-97.
;
;-
;____________________________________________________________________________

PRO FMD_READ_DATA,arc_file,id,pixel,obs_time,spec,nifgs, $
	          glon,glat,model,f,n_xxx,b_xxx,l_xxx, $
		  cal_spec,cal_nifgs,cal_obs_time,xcal,ical,refh,skyh, $
		  dihd,galcut=galcut,weight_cor=weight_cor, $
	          sky_glitch=sky_glitch,cal_glitch=cal_glitch, $
		  sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
		  sky_lbl=sky_lbl,cal_lbl=cal_lbl,splen=splen, $
		  sky_dihd=sky_dihd,sky_s0=sky_s0,cal_s0=cal_s0, $
                  time=time,scan=scan,fsl_idx=fsl_idx


eject = timeconv('893251118',infmt='z',outfmt='s')
		; get cover eject time in seconds

r2d = 180 / !dpi

IF (N_ELEMENTS(weight_cor) EQ 0) THEN  weight_cor = 1
print, 'Setting IFG weight to', strcompress( weight_cor )


; Open fsl file list
; ------------------
print,'Reading fsl archive times listing'

arc_file = [arc_file]
id = [STRUPCASE(id)]
	; make array if scalar

FOR k=0,N_ELEMENTS(arc_file)-1 DO BEGIN

  file = STRARR(40)
  n_files = 0
  a = ''

  OPENR,1,'csdr$firas_ref:'+arc_file(k)+'.txt'

  WHILE (NOT EOF(1)) DO BEGIN
    READF,1,a
    file(n_files) = 'csdr$firas_insky:fsl_sky_' + id(k) + '.' + a
    n_files = n_files + 1
  ENDWHILE
  CLOSE,1

  file = file(0:n_files-1)
  ;

  OPENR,2,file(0)
  finfo = FSTAT(2)
  nrec0 = finfo.size / finfo.rec_len
  start_idx = INDGEN(nrec0)
  CLOSE,2
  ;
  FOR i=1,n_files-1 DO BEGIN
   OPENR,2,file(i)
   finfo = FSTAT(2)
   nrec0 = finfo.size / finfo.rec_len
   start_idx = [start_idx,INDGEN(nrec0)]
   CLOSE,2
  ENDFOR
  ;

  file = STRUPCASE(file)

  ; Determine calibration model and nyquist frequency
  ; -------------------------------------------------
  IF (k EQ 0) THEN BEGIN
   PRINT,'Reading calibration model solution and nyquist frequency'
   bufsz = 11264L
   ndim = [0,0,1]
   dimlist = [40,361L]
   off = [2512,462,2660L]
   type = [7,4,6]

   nrec = RSM(file(0),bufsz,ndim,dimlist,off,type,label,nyq,s)

   model = label(0)
   nyquist = nyq(0)
   nzero = WHERE(FLOAT(s(*,0)) NE 0)

   ; Build frequency vector
   ; ----------------------
   n_freq = N_ELEMENTS(nzero)
   bins = STRCOMPRESS(STRING(MIN(nzero)+1) + ':' + STRING(MAX(nzero)+1),/REM)
   f = (MIN(nzero) + FINDGEN(n_freq)) * nyquist / splen

  ENDIF
  ;

  ; Read Data from fsl skymap files
  ; -------------------------------
  ndim = [0,0,1,1,0,0,0,0,0,0,0,0]
  dimlist = [14L,2,n_freq]
  off = [11184L,0,14,2660+nzero(0)*8,132,134,11238,11240,11250,450,2310,2632]
  type = [3,7,3,6,2,4,2,2,2,4,4,4]

  nrec = RSM(file,bufsz,ndim,dimlist,off,type,pixel0,gmt,adt,spec0,nifgs0, $
  	          sky_wgts0,lon,lat,scan0,sky_glitch0,sky_dihd0,sky_s00,fnum)

  tme = timeconv(gmt,infmt='z',outfmt='s')

  ; Special Sort by pixel number, file number, record number
  ; --------------------------------------------------------
  dummy = 220000L * pixel0 + 20000L * fnum + start_idx
  pix_sort = SORT(dummy)
  
  IF (k EQ 0) THEN BEGIN

    fsl_idx = LINDGEN(N_ELEMENTS(pix_sort))
    obs_time = (tme(pix_sort) - eject) / 86400.
    time = FMD_ADT_SEC(adt)
    time = time(pix_sort)
    glon = lon(pix_sort) * r2d * 1.e-4
    glat = lat(pix_sort) * r2d * 1.e-4
    scan = scan0(pix_sort) * r2d * 1.e-4
    sky_lbl = STRMID(file(fnum(pix_sort)),STRPOS(file(0),'.ED')-4,4)
    pixel = pixel0(pix_sort)
    nifgs = nifgs0(pix_sort)
    sky_wgts = sky_wgts0(pix_sort)
    sky_glitch = sky_glitch0(pix_sort)
    sky_dihd = sky_dihd0(pix_sort)
    sky_s0 = sky_s00(pix_sort)
 
    ; Sort spectra frequency by frequency
    spec0 = FLOAT(spec0)
    FOR j=0,n_freq-1 DO BEGIN
	slice = spec0(j,*)
  	spec0(j,*) = slice(pix_sort)
    ENDFOR
    ;
    spec = TEMPORARY(spec0)
    ;
 
  ENDIF

  IF (k GT 0) THEN BEGIN

    fsl_idx = [fsl_idx,LINDGEN(N_ELEMENTS(pix_sort))]
    obs_time = [obs_time,(tme(pix_sort) - eject) / 86400.]
    time0 = FMD_ADT_SEC(adt)
    time = [time,time0(pix_sort)]
    glon = [glon,lon(pix_sort) * r2d * 1.e-4]
    glat = [glat,lat(pix_sort) * r2d * 1.e-4]
    scan = [scan,scan0(pix_sort) * r2d * 1.e-4]
    sky_lbl = [sky_lbl,STRMID(file(fnum(pix_sort)),STRPOS(file(0),'.ED')-4,4)]
    pixel = [pixel,pixel0(pix_sort)]
    nifgs = [nifgs,nifgs0(pix_sort)]
    sky_wgts = [sky_wgts,sky_wgts0(pix_sort)]
    sky_glitch = [sky_glitch,sky_glitch0(pix_sort)]
    sky_dihd = [sky_dihd,sky_dihd0(pix_sort)]
    sky_s0 = [sky_s0,sky_s00(pix_sort)]
 
    ; Sort spectra frequency by frequency
    spec0 = FLOAT(spec0)
    FOR j=0,n_freq-1 DO BEGIN
	slice = spec0(j,*)
  	spec0(j,*) = slice(pix_sort)
    ENDFOR
    ;
    spec = [[spec],[TEMPORARY(spec0)]]
    ;
 
  ENDIF

ENDFOR
;


; Deweight low S0 coadds
; ----------------------
i = WHERE((ABS(sky_s0) LT 1 AND (STRUPCASE(STRMID(sky_lbl,1,3)) EQ 'LFS' OR $
                                STRUPCASE(STRMID(sky_lbl,1,3)) EQ 'LFL')) OR $
          (ABS(sky_s0) LT .23 AND (STRUPCASE(STRMID(sky_lbl,1,3)) EQ 'HSF' OR $
                                STRUPCASE(STRMID(sky_lbl,1,3)) EQ 'HLF')))
IF (i(0) NE -1) THEN sky_wgts(i) = 0


; Adjust Sky Weights
; ------------------
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

i = WHERE(STRMID(sky_lbl,2,2) EQ 'FS')
IF (i(0) NE -1) THEN BEGIN
	fac_fs = TOTAL(nifgs(i)) / TOTAL(sky_wgts(i))
	sky_wgts(i) = sky_wgts(i) * fac_fs
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
FOR i=0L,N_ELEMENTS(sky_wgts)-1 DO BEGIN
	sc_mode = STRPOS('SF_FS_LF_FL',STRMID(sky_lbl(i),2,2))/6
	sc_mode = sc_mode < 1
	sky_wgts(i) = sky_wgts(i) * (weight_cor ^ (0.5 - sc_mode))
ENDFOR
;

; Generate sixpack of weights
; ---------------------------
pix_dist = WHERE(HISTOGRAM(MIN=0,pixel) GT 0)
n_dist = N_ELEMENTS(pix_dist)
		; get # of distinct pixels

num_pix = FLTARR(n_dist)

sq = SORT(pixel)
pxq = pixel(sq)
wq = sky_wgts(sq)

j = 0
num_pix(0) = wq(0)


FOR k=1L,N_ELEMENTS(pixel)-1 DO BEGIN	; loop over all obs
	IF (pxq(k) NE pxq(k-1)) THEN j = j + 1	; new pixel
	num_pix(j) = num_pix(j) + wq(k)
ENDFOR

print, 'Total # of IFGs' + Strcompress( total(sky_wgts) )
print, 'Total from num_pix' + strcompress( total(num_pix) )

pix2xy,/six,pix_dist,data=num_pix,res=6,raster=n_xxx
;

; Generate sixpack of weighted galactic longitude and latitude
; ------------------------------------------------------------
print,'b_xxx'
uv = coorconv([[glon(sq)],[glat(sq)]],infmt='l',outfmt='u')

acm_x = FLTARR(n_dist)
acm_y = acm_x
acm_z = acm_x
j = 0
acm_x(0) = uv(0,0)*wq(0)
acm_y(0) = uv(0,1)*wq(0)
acm_z(0) = uv(0,2)*wq(0)

FOR k=1L,N_ELEMENTS(pixel)-1 DO BEGIN	; loop over all obs
	IF (pxq(k) NE pxq(k-1)) THEN j = j + 1	; new pixel
	acm_x(j) = acm_x(j) + uv(k,0)*wq(k)
	acm_y(j) = acm_y(j) + uv(k,1)*wq(k)
	acm_z(j) = acm_z(j) + uv(k,2)*wq(k)
ENDFOR

acm_x = acm_x / num_pix
acm_y = acm_y / num_pix
acm_z = acm_z / num_pix

temp = coorconv([[acm_x],[acm_y],[acm_z]],infmt='u',outfmt='l')

pix2xy,/six,pix_dist,data=temp(*,0),res=6,raster=l_xxx
pix2xy,/six,pix_dist,data=temp(*,1),res=6,raster=b_xxx
;

; Read Data from fsl cal files
; ----------------------------
PRINT,'Reading Cal Data'

ndim = [1,0,0,0,0,0,0,0,0,0,0]
dimlist = [n_freq,14L]
off = [2660L+nzero(0)*8,0,132,134,2294,2298,2302,2306,2310,450,2632]
type = [6,7,2,4,4,4,4,4,4,4,4]


FOR j=0,N_ELEMENTS(id)-1 DO BEGIN

    FOR k=1,7 DO BEGIN

	cfile = 'csdr$firas_incal:fsl_cal_' + id(j) + '.dihed_' $
		+ STRING(k)
	cfile = STRCOMPRESS(cfile,/REMOVE_ALL)

	nrec = RSM(cfile,bufsz,ndim,dimlist,off,type,s,gmt,n,n_adj, $
		   xc,ic,sh,rh,dh,gl,s0)

	IF (j+k EQ 1) THEN BEGIN

		tme = timeconv(gmt,infmt='z',outfmt='s')
		cal_obs_time = (tme - eject) / 86400.

		cal_spec = FLOAT(s)
		cal_nifgs = n

		xcal = xc
		ical = ic
		skyh = sh
		refh = rh
		dihd = dh

		cal_lbl = REPLICATE(id(j),nrec)
		cal_glitch = gl

		cal_wgts = n_adj

		cal_s0 = s0

	ENDIF ELSE BEGIN

		tme = timeconv(gmt,infmt='z',outfmt='s')
		cal_obs_time = [cal_obs_time,(tme - eject) / 86400.]

		cal_spec = [[cal_spec],[FLOAT(s)]]
		cal_nifgs = [cal_nifgs,n]

		xcal = [xcal,xc]
		ical = [ical,ic]
		skyh = [skyh,sh]
		refh = [refh,rh]
		dihd = [dihd,dh]

		cal_lbl = [cal_lbl,REPLICATE(id(j),nrec)]
		cal_glitch = [cal_glitch,gl]

		cal_wgts = [cal_wgts,n_adj]

		cal_s0 = [cal_s0,s0]

	ENDELSE

    ENDFOR
ENDFOR


; Deweight low S0 coadds
; ----------------------
i = WHERE((ABS(cal_s0) LT 1 AND (STRUPCASE(STRMID(cal_lbl,1,3)) EQ 'LFS' OR $
                                STRUPCASE(STRMID(cal_lbl,1,3)) EQ 'LFL')) OR $
          (ABS(cal_s0) LT .23 AND (STRUPCASE(STRMID(cal_lbl,1,3)) EQ 'HSF' OR $
                                STRUPCASE(STRMID(cal_lbl,1,3)) EQ 'HLF')))
IF (i(0) NE -1) THEN cal_wgts(i) = 0



; Renormalize cal weights by sky weight factor
; --------------------------------------------
i = WHERE(STRMID(cal_lbl,2,2) EQ 'SS')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_ss

i = WHERE(STRMID(cal_lbl,2,2) EQ 'SF')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_sf

i = WHERE(STRMID(cal_lbl,2,2) EQ 'FS')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_fs

i = WHERE(STRMID(cal_lbl,2,2) EQ 'LF')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_lf

i = WHERE(STRMID(cal_lbl,2,2) EQ 'FL')
IF (i(0) NE -1) THEN cal_wgts(i) = cal_wgts(i) * fac_fl





FOR i=0L,N_ELEMENTS(cal_wgts)-1 DO BEGIN
	sc_mode = STRPOS('SF_FS_LF_FL',STRMID(cal_lbl(i),2,2))/6
	sc_mode = sc_mode < 1
	cal_wgts(i) = cal_wgts(i) * (weight_cor ^ (0.5 - sc_mode))
ENDFOR


RETURN
END
 
