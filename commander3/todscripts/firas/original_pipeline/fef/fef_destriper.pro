FUNCTION FEF_Destriper, freq,pixel,spec,cal_spec,nifgs,cal_nifgs,xcal, $
	                 tm=tm,cal_tm=cal_tm,tau=tau,galcut=galcut,vib=vib, $
	 		 step_up=step_up,step_dn=step_dn,ejv=ejv,afp=afp, $
			 index=index,id=id,solution=solution,ndf_off=ndf_off, $
                         ndf_none=ndf_none,csp=csp,del_temp=del_temp, $
                         pcvr=pcvr,sqr_mat=sqr_mat,rect=rect,diag=diag, $
                         c_calsp=c_calsp,sky_wgts_ds=sky_wgts_ds

;
;  Modification History:
;   Written by Ken Jensen, Hughes STX, 08 Sept 1994, as a modification
;     of FDS_DESTRIPER.PRO,  modification is that the code for writing
;     the FEX_EJV data set is deleted.
;

; Temporary Arrays for Uncorrected Spectra
; ----------------------------------------
sp_temp = spec
cal_sp_temp = cal_spec


m2ep = 2.99792458d-7
eject = timeconv('893251118',infmt='z',outfmt='s')


; Transpose input spectral arrays to get frequencies along rows
; -------------------------------------------------------------
spec = TRANSPOSE(spec)
cal_spec = TRANSPOSE(cal_spec)


; Get # of frequencies, # of observations, # of pixels
; ----------------------------------------------------
sz = SIZE(spec)
n_sky = sz(1)
n_freq = n_ELEMENTS(freq)
n_pixels = MAX(pixel)+1

sz = SIZE(cal_spec)
n_cal = sz(1)

n_obs = n_sky + n_cal



; Number of offset Heaviside (tophat) input functions
; ---------------------------------------------------
n_tophat = N_ELEMENTS(step_up)


; Number of UFO (exponential) input functions
; -------------------------------------------
n_ufo = N_ELEMENTS(tau)


; Number of vibration corrections
; -------------------------------
n_vib = KEYWORD_SET(vib)
if(n_vib eq 0)then vib = dblarr(5)

n_ufo_v = n_ufo + n_vib


; Total number of offset fit functions (tophat + UFO + Vibration)
; ---------------------------------------------------------------
n_f  = n_tophat + n_ufo + n_vib


; Spectral weights
; ----------------
nifg_ = [nifgs,cal_nifgs]


; Allocate space for Planck function
; ----------------------------------
plnk = DBLARR(n_obs)


; Build galactic mask, Get NDF
; ----------------------------
ll = coorconv(pixel,infmt='p',inco='f',outfmt='l',outco='g')
mask = 0 * nifgs
mask(WHERE((ABS(ll(*,1)) GT galcut)AND(sky_wgts_ds gt 0.))) = 1
	; Mask out observations within galcut degrees of galactic plane
	; Mask out observations with zero weight

nmsk = nifgs * mask

pix_nogal = pixel(WHERE(mask EQ 1))
hist_nogal = HISTOGRAM(min=0,max=MAX(pixel),pix_nogal)
n_dist = n_ELEMENTS(WHERE(hist_nogal GT 0))
	; Number of observed pixels


; Number of degrees of freedom for pure offset destriper
; ------------------------------------------------------
ndf_off = TOTAL(mask) - n_dist - n_f


; Number of degrees of freedom for no destriper
; ---------------------------------------------
ndf_none = TOTAL(mask) - n_dist


; Allocate space for sigma arrays
; -------------------------------
sig = FLTARR(n_freq,2)


; Modify XCAL temp by "official" temperature offset
; -------------------------------------------------
temp = xcal - del_temp


; Include all cal data
; --------------------
mask = [mask,INTARR(n_cal)+1]


; Generate Orthogonalized Exponentials & Vibration
; ------------------------------------------------
tme = DOUBLE([tm,cal_tm])

IF (n_ufo_v NE 0) THEN BEGIN

	orth = DBLARR(N_ELEMENTS(tme),n_ufo_v)
	ufo_vib = DBLARR(N_ELEMENTS(tme),n_ufo_v)

	FOR i=0,n_ufo-1 DO ufo_vib(*,i) = exp(-[tm,cal_tm]/tau(i))

	IF (n_vib EQ 1) THEN BEGIN
		vib_fnc = vib(4)
		FOR i=3,0,-1 DO vib_fnc = vib_fnc * (tme/365) + vib(i)
		ufo_vib(*,n_ufo) = vib_fnc
	ENDIF

	FOR i=0,n_ufo_v-1 DO orth(*,i) = ufo_vib(*,i) / $
						SQRT(TOTAL(ufo_vib(*,i)^2))

	FOR i=1,n_ufo_v-1 DO BEGIN
		FOR j=0,i-1 DO orth(*,i) = orth(*,i) - $
			TOTAL(orth(*,i)*orth(*,j))*orth(*,j)
		orth(*,i) = orth(*,i) / SQRT(TOTAL(orth(*,i)^2))
	ENDFOR

ENDIF


; Start of each mission period
; ----------------------------
sup = (timeconv(step_up,infmt='z',outfmt='s') - eject) / 86400.

; Stop of each mission period
; ---------------------------
sdn = (timeconv(step_dn,infmt='z',outfmt='s') - eject) / 86400.


; Generate Heaviside for Offset
; -----------------------------
tophat = DBLARR(n_obs,n_tophat)
FOR i=0,n_tophat-1 DO $
	tophat(WHERE( (tme GE sup(i)) AND (tme LT sdn(i)) ),i) = 1
	; Set to 1 if observation within ith mission period


; Generate Functions for Offset
; -----------------------------
IF (n_ufo_v NE 0) THEN BEGIN
	func = orth
	FOR i=0,n_tophat-1 DO func = [[func],[tophat(*,i)]]
		; Append tophat

	IF (n_ufo_v EQ 1) THEN f_conv = 1 / (TRANSPOSE(orth) # ufo_vib) $
			  ELSE f_conv = INVERT(TRANSPOSE(orth) # ufo_vib)
		; Generate orthog to non-orthog conversion matricies

ENDIF ELSE func = tophat



; Schematic of Least-Square Fit Matrix
; ------------------------------------


; Evaluating _f, f_, f_f
; ----------------------
_n = SQRT(nifg_)

_f = func
FOR k=0,n_f-1 DO _f(*,k) = _f(*,k) * _n
f_ = TRANSPOSE(_f)

f_f = f_ # _f

f_f_off = f_(*,WHERE(mask ne 0)) # _f(WHERE(mask ne 0),*)


; Evaluating a_a, f_a, a_f
; ------------------------
hist = HISTOGRAM(min=0,pixel)
sent = FLOAT(-1e7)
num_pix = hist(WHERE(hist GT 0))
one = INTARR(2*n_f) + 1

pix = LONG(pixel)
dat = _n^2
pixavg,pix,dat,sent
diag = dat * num_pix

nz_pix = pix
n_nz = n_ELEMENTS(nz_pix)
nog = WHERE(hist_nogal(nz_pix) NE 0)

a_f = DBLARR(n_nz,n_f)
FOR i=0,n_f-1 DO BEGIN
	pix = LONG(pixel)
	dat = _f(*,i) * _n
	pixavg,pix,dat,sent
	a_f(*,i) = dat * num_pix
ENDFOR

f_a = TRANSPOSE(a_f)


; Get obs -> compressed pixel index code
; -------------------------------------- 
pix2xy,nz_pix,data=INDGEN(n_nz),res=6,/six,ras=code
code = pix2dat(pix=INDGEN(n_pixels),ras=code)
code = code(pixel)


; Compute Pure Offset Inverse
; ---------------------------
d_inv = 1 / (diag # one)
q_off = f_f_off - (f_a(*,nog) # (d_inv(nog,*) * a_f(nog,*)))
q_inv = INVERT(q_off)
q_inv_off = [[q_inv],[0*q_inv]]
pcvr = q_inv


; Spectrum Dependent Terms
; ------------------------
off  = COMPLEXARR(n_f,n_freq)
	; Allocate destriper correction spectra arrays

pix_spec = COMPLEXARR(n_nz,n_freq)
	; Allocate pixel spectrum arrays


FOR l=0,n_freq-1 DO BEGIN
print,'freq ',l


	; Extract spectra at given frequency, calculate Planck, b_f_off
	; -------------------------------------------------------------
	sp = spec(*,l)
	cal_sp = cal_spec(*,l)

	; Real part of spectra * sqrt(wgts)
        ; ---------------------------------
	snr = DOUBLE(FLOAT([sp,cal_sp])) * _n

	; Imag part of spectra * sqrt(wgts)
        ; ---------------------------------
	sni = DOUBLE(IMAGINARY([sp,cal_sp])) * _n


	; Planck at (adjusted) XCAL temp for each obs
        ; -------------------------------------------
	FOR i=0,n_cal-1 DO plnk(n_sky+i) = $
		planck(temp(i),freq(l),units='i',/mjy) * m2ep

	_p  = plnk * _n

	b_f_off = [[f_ # ((snr -_p)*mask)],[f_ # (sni*mask)]]


	; Compute b_a_off
	; ---------------
	pix = LONG(pixel)
	dat = snr * _n
	pixavg,pix,dat,sent
	b_a_off = dat * num_pix

	pix = LONG(pixel)
	dat = sni * _n
	pixavg,pix,dat,sent
	b_a_off = [[b_a_off],[dat * num_pix]]


	; Compute pure offset parameters
	; ------------------------------
	off0 = dpc_math(b_f_off - $
			f_a(*,nog) # (d_inv(nog,*) * b_a_off(nog,*)), $
		        q_inv_off,code='m')

	off0 = REFORM(off0,n_f,2)

	IF (index(l) EQ -1) THEN off0 = 0 * off0
		; If no destriping then set offset to 0

	ap_off = dpc_math( $
			[[(b_a_off(*,0) - (off0(*,0) # f_a)) * d_inv], $
			 [(b_a_off(*,1) - (off0(*,1) # f_a)) * d_inv]],$
		 code='d2c')

	off0   = dpc_math(off0,code='d2c')
			; Convert to complex


	; Pure Offset
        ; -----------
	off_set = func # off0
	c_spc = spec(*,l) - off_set
	res_r = TOTAL(FLOAT(c_spc-ap_off(code))^2 * nmsk)
	res_i = TOTAL(IMAGINARY(c_spc-ap_off(code))^2 * nmsk)

	IF (index(l) EQ 0) THEN $
		sig(l,*) = [SQRT(res_r/ndf_off),SQRT(res_i/ndf_off)] $
	ELSE $
		sig(l,*) = [SQRT(res_r/ndf_none),SQRT(res_i/ndf_none)]



	IF (n_ufo_v NE 0) THEN BEGIN
		off(0:n_ufo_v-1,l) = f_conv # off0(0:n_ufo_v-1)
		off(n_ufo_v:*,l)   = off0(n_ufo_v:*)
	ENDIF ELSE off(*,l) = off0



	; Extract pixel spectra from desired destriper
	; --------------------------------------------
	pix_spec(*,l) = ap_off

ENDFOR
	; Frequency Loop


; Extract corrections from desired destriper
; ------------------------------------------
ejv = off

FOR i=0,n_freq-1 DO BEGIN
	IF (index(i) EQ -1) THEN ejv(*,i) = 0
ENDFOR

ejv = TRANSPOSE(ejv)


; Rasterize Pixel Spectra
; -----------------------
pix2xy,/six,nz_pix,data=pix_spec,res=6,raster=afp



; Compute corrected spectra
; -------------------------

IF (n_ufo_v NE 0) THEN func = ufo_vib
func = [[func],[tophat]]
func = TRANSPOSE(func)


csp = 0 * sp_temp
FOR i=0,N_ELEMENTS(tm)-1 DO csp(*,i) = sp_temp(*,i) - ejv # func(*,i)


c_calsp = 0 * TRANSPOSE(cal_spec)
FOR i=0,N_ELEMENTS(cal_tm)-1 DO c_calsp(*,i) = cal_sp_temp(*,i) - $
						  ejv # func(*,n_sky+i)



; Correct Parameter Covariance Matrix for Orthogonalization
; ---------------------------------------------------------
IF (n_ufo_v NE 0) THEN BEGIN
	f_conv0 = DBLARR(n_f,n_f)
	f_conv0(INDGEN(n_f)*(n_f+1)) = 1
	f_conv0(0:n_ufo_v-1,0:n_ufo_v-1) = f_conv

	pcvr = f_conv0 # pcvr # TRANSPOSE(f_conv0)
	rect = a_f # INVERT(f_conv0)
	sqr_mat = TRANSPOSE(INVERT(f_conv0)) # f_f_off # INVERT(f_conv0)
ENDIF


RETURN,sig
END
