FUNCTION FMD_DESTRIPER, freq=freq,pixel=pixel,tm=tm,cal_tm=cal_tm, $
		        spec=spec,cal_spec=cal_spec, $
                        diff_spec=diff_spec,diff_cal_spec=diff_cal_spec, $
		        sky_wgts=sky_wgts,cal_wgts=cal_wgts, $
                        frac_wgt=frac_wgt,max_frac=max_frac, $
			sky_mask=sky_mask,cvec_mask=cvec_mask, $
                        sky_idx=sky_idx,cal_idx=cal_idx, $
			xcal=xcal,del_temp=del_temp, $
                        convg=convg,iter=iter,n_iter=n_iter, $
			ejg=ejg,afp=afp,chi2_chan=chi2_chan, $
                        func=func,sky_func=sky_func,cal_func=cal_func, $
			pcvr=pcvr,rect=rect,diag=diag,chi2_map=chi2_map, $
                        n_krnl=n_krnl,n_cmn=n_cmn,square=square,d_inv=d_inv, $
			num_coadd=num_coadd,ndf=ndf,pure_off=pure_off,mjy=mjy

;
;  Modification History:
;        FDS_DESTRIPER Delivered by Joel Gales, ARC, Dec 1995
;        Renamed FMD_DESTRIPER by K.Jensen, Hughes STX, 12/11/95
;        Added SKY_MASK keyword.  K.Jensen, 4/3/96.
;        Added FUNC keyword, all kernels now passed in. K.Jensen, 6/7/96.
;        Added N_KRNL and N_CMN keywords. K.Jensen, 7/31/96.
;        Removed orthogonalization of kernels. Destriping will be performed
;                in the original basis. K.Jensen, 9/17/96.
;        Added CVEC_MASK keyword.  K.Jensen, 10/3/96.
;        IF Non-Zero FRAC_WGT, NDF = NOBS - N_STRIPES,  K.Jensen, 4/7/97.
;        Add NUM_COADD Keyword and calculate NUM_COADD,  K.Jensen, 5/23/97.
;        CHI2_MAP is now total chi-squared per pixel, previously was
;             chi-squared per degree of freedom,  K.Jensen, 5/23/97.
; 

; Define Planck Function Conversion Factor
; ----------------------------------------
IF (KEYWORD_SET(mjy)) THEN $
	unit_conv = 1 $
ELSE $
	unit_conv = 2.99792458d-7


IF (KEYWORD_SET(pure_off)) THEN convg = 0.0

IF (NOT KEYWORD_SET(frac_wgt)) THEN frac_wgt = 0. * sky_wgts

sz = SIZE(spec)
n_sky = sz(2)
	; # of sky spectra

sz = SIZE(cal_spec)
n_cal = sz(2)
	; # of cal spectra

; Build channel index arrays if not defined (Single channel mode)
; ---------------------------------------------------------------
IF (N_ELEMENTS(sky_idx) EQ 0) THEN sky_idx = BYTARR(n_sky)
IF (N_ELEMENTS(cal_idx) EQ 0) THEN cal_idx = BYTARR(n_cal)

idx = [sky_idx,cal_idx]


; Set up various utility variables
; --------------------------------
n_freq = N_ELEMENTS(freq)
	; # of frequency points

n_chan = MAX(sky_idx) + 1
	; # of channels

n_obs = n_sky + n_cal
	; Total # of spectra

n_uncmn = n_krnl - n_cmn            ; # of Uncommon Kernels
;
n_func = n_uncmn * n_chan + n_cmn   ; # of Stripes
;

max_iter = 20
IF (NOT KEYWORD_SET(pure_off)) THEN $
 max_iter = iter
	; Maximum number of gain iterations


; Transpose spectrum arrays
; -------------------------
spec = TRANSPOSE(spec)
cal_spec = TRANSPOSE(cal_spec)
	; Columns - coadds, Rows - freqs

no_diff = 1
IF (NOT KEYWORD_SET(pure_off)) THEN BEGIN
 IF (N_ELEMENTS(diff_spec) NE 0) THEN BEGIN
	no_diff = 0
	diff_spec = TRANSPOSE(diff_spec)
	diff_cal_spec = TRANSPOSE(diff_cal_spec)
 ENDIF ELSE no_diff = 1
	; If no differential spectra then set no_diff flag
ENDIF


; Find Pixels with Non-zero Weight
; --------------------------------
good = WHERE(sky_wgts GT 0.)
nz_pix = WHERE(HISTOGRAM(min=0,pixel(good)) NE 0)
n_nz = N_ELEMENTS(nz_pix)
;
IF (n_nz LE 0) THEN BEGIN
 PRINT,'FMD_DESTRIPER : Error : No Pixels with non-zero Weight !'
 c_vec = FLTARR(n_freq) - 1.
 RETURN,c_vec
ENDIF
;

; Find Pixels with Non-zero Weight*Mask
; -------------------------------------
goodm = WHERE(sky_wgts*sky_mask GT 0.,ncoadd)
nzm_pix = WHERE(HISTOGRAM(min=0,pixel(goodm)) NE 0)
n_nzm = N_ELEMENTS(nzm_pix)
;
IF (n_nzm LE 0) THEN BEGIN
 PRINT,'FMD_DESTRIPER : Error : No Pixels with non-zero Weight*Mask !'
 c_vec = FLTARR(n_freq) - 1.
 RETURN,c_vec
ENDIF
;

; Destriping Degrees of Freedom
; -----------------------------
ndf = 1. * ncoadd - n_nzm - n_func
IF (NOT KEYWORD_SET(pure_off)) THEN ndf = ndf - n_chan
;

; Destriping Degrees of Freedom per CHANSCAN
; ------------------------------------------
ndf_chan = INTARR(n_chan)
FOR i=0,n_chan-1 DO BEGIN
 dummy = WHERE((sky_idx EQ i)and(sky_wgts*sky_mask GT 0.),ncoadd)
 dummy2 = WHERE(HISTOGRAM(min=0,pixel(dummy)) NE 0)
 ndf_chan(i) = ncoadd - N_ELEMENTS(dummy2) - n_krnl
ENDFOR
;

; Coadds with Non-zero Weight*Sky_Mask and Good Frac_Wgt
; ------------------------------------------------------
goodf = WHERE((sky_wgts*sky_mask GT 0.)and(frac_wgt LE max_frac),n_frac)
IF (n_frac LE 0) THEN BEGIN
 PRINT,'FMD_DESTRIPER : No Coadds with non-zero Weight*CVec_Mask and Good Frac_Wgt !'
 c_vec = FLTARR(n_freq) - 1.
 RETURN,c_vec
ENDIF
nzf_pix = WHERE(HISTOGRAM(min=0,pixel(goodf)) NE 0)
;

; Pixels with Non-zero Weight*CVec_Mask and Good Frac_Wgt
; -------------------------------------------------------
goodc = WHERE((sky_wgts*cvec_mask GT 0.)and(frac_wgt LE max_frac),n_cvec)
nzc_pix = WHERE(HISTOGRAM(min=0,pixel(goodc)) NE 0)
n_nzc = N_ELEMENTS(nzc_pix)
;
IF (n_nzc LE 0) THEN BEGIN
 PRINT,'FMD_DESTRIPER : Error : No Pixels with non-zero Weight*CVec_Mask !'
 c_vec = FLTARR(n_freq) - 1.
 RETURN,c_vec
ENDIF
;

; C-Vector Degrees of Freedom
; ---------------------------
ndf_cvec = n_cvec - n_nzc - n_func
IF (TOTAL(frac_wgt) GT 0.) THEN ndf_cvec = n_cvec - n_func
IF (NOT KEYWORD_SET(pure_off)) THEN ndf_cvec = ndf_cvec - n_chan
;

; Number of Coadds per Pixel (for Chi-Squared)
; --------------------------------------------
np = FIX(FMD_PIX_SUM(pixel=pixel(goodf),data=pixel(goodf)*0+1))
;

; Extend calibration observations to mask
; ---------------------------------------
mask = [sky_mask,BYTARR(n_cal)+1B]


; Compute _f, f_
; --------------
_f = func
_wgt = SQRT([sky_wgts,cal_wgts])
FOR k=0,n_krnl-1 DO _f(*,k) = _f(*,k) * _wgt
f_ = TRANSPOSE(_f)
	; _f = func * SQRT(wgt)
	; f_ = TRANSPOSE(_f)
	;
	; Note: Variables with leading "_" have been multiplied by _wgt



; Compute SQUARE (= F # wgt # F)
; ------------------------------
square = DBLARR(n_func,n_func)
	; Square part of curvature matrix

one = INTARR(n_func) + 1
	; Replication utility vector: a # one = [[a],...[a]]


; Case of No Common Kernels
; -------------------------
IF (n_cmn LT 1) THEN BEGIN
;
	FOR i=0,n_chan-1 DO BEGIN
		j = WHERE(idx EQ i)
		k = i * n_krnl
		l = k + n_krnl - 1
		square(k:l,k:l) = f_(*,j) # (_f(j,*) * (mask(j)#one))
	ENDFOR
;
ENDIF


; Case of Common Kernels
; ----------------------
IF (n_cmn GE 1) THEN BEGIN
;
	square(0:n_cmn-1,0:n_cmn-1) = f_(0:n_cmn-1,*) # $
					(_f(*,0:n_cmn-1) * (mask#one(0:n_cmn-1)))
			;  Common # wgt # Common terms
;
	FOR i=0,n_chan-1 DO BEGIN
		j = WHERE(idx EQ i)
		k = i * n_uncmn + n_cmn
		l = k + n_uncmn - 1
		square(0:n_cmn-1,k:l) = f_(0:n_cmn-1,j) # $
					(_f(j,n_cmn:*) * (mask(j)#one(n_cmn:*)))
			; Common # wgt # Uncommon terms
;
		square(k:l,0:n_cmn-1) = TRANSPOSE(square(0:n_cmn-1,k:l))
			; Uncommon # wgt # Common terms
;
		square(k:l,k:l) = f_(n_cmn:*,j) # $
					(_f(j,n_cmn:*) * (mask(j)#one(n_cmn:*)))
			; Uncommon # wgt # Uncommon terms
;
	ENDFOR
;
ENDIF
;

; Allocate output and utility arrays
; ----------------------------------
ejg = DBLARR(n_func + n_chan,n_freq)
pix_spec = DBLARR(n_nz,n_freq)
n_iter = INTARR(n_freq)
c_vec = FLTARR(n_freq)
chi2_pix = FLTARR(n_nz)
chi2_chan = FLTARR(n_freq,n_chan)
_pl_ref = DBLARR(n_obs)
;

wmask = sky_wgts * sky_mask
;

; Loop over frequencies
; ---------------------
FOR l=0,n_freq-1 DO BEGIN
print,'freq ',l,freq(l)

;  Subtract Mean Sky Spectrum from Sky Coadds
;  ------------------------------------------
    mean_sp = TOTAL(spec(*,l) * wmask) / TOTAL(wmask)
    sp = DOUBLE(spec(*,l)) - mean_sp
    cal_sp = DOUBLE(cal_spec(*,l))

    IF (no_diff) THEN BEGIN
       _spc = [sp,cal_sp] * _wgt
       _ref = 0 * _wgt
    ENDIF ELSE BEGIN
       _spc = DOUBLE([diff_spec(*,l),diff_cal_spec(*,l)]) * _wgt
       _ref = [sp,cal_sp] * _wgt - _spc
    ENDELSE
    ; If differential spectra not used (pure offset)
    ; Then combine full sky and cal spectra and set ref to 0
    ; Else combine diff sky and cal spectra and compute ref


    ; Compute Planck - Reference ([bP - R] in design document)
    ; --------------------------------------------------------
    _pl_ref(0:n_sky-1) = -_ref(0:n_sky-1)
	; No planck for sky spectra

    FOR i=0,n_cal-1 DO _pl_ref(n_sky+i) = $
        PLANCK(xcal(i)-del_temp,freq(l),units='i',/mjy) * $
 	 unit_conv * _wgt(n_sky+i) - _ref(n_sky+i)


    ; Initialize gain, chi_squared
    ; ----------------------------
    chan_gain = FLTARR(n_chan) + 1
    gain = chan_gain(idx)
    chi2 = 1e9
       ; CHAN_GAIN - gain for each channel
       ; GAIN - gain for each observation


    ; Gain Iteration
    ; --------------
    REPEAT BEGIN

        chi2_last = chi2   ; Store previous chi-squared

        gain_sky = gain(0:n_sky-1)  ; Get gains for sky observations


        ; If first frequency or gain destriping 
        ; then compute RECT, DIAG, Inverses
        ; -------------------------------------
        IF (l EQ 0 OR KEYWORD_SET(pure_off) EQ 0) THEN BEGIN

	    ; Compute RECT
	    ; ------------
	    rect = DBLARR(n_nz,n_func)

	    IF (n_cmn GE 1) THEN BEGIN
	        FOR i=0,n_cmn-1 DO BEGIN
	            data = gain * _wgt * _f(*,i)
	            dat = FMD_PIX_SUM(pixel=pixel(good),data=data(good),cmp_px=cmp_px)
     		    PIX2XY,cmp_px,data=dat,res=6,/six,ras=ras
		    rect(*,i) = PIX2DAT(pix=nz_pix,ras=ras)
	        ENDFOR
	    ENDIF

	    FOR i=0,n_chan-1 DO BEGIN
	        j = WHERE((sky_wgts GT 0) and(sky_idx EQ i))
	        k0 = i * (n_krnl - n_cmn) + n_cmn
	        FOR k=0,n_krnl-n_cmn-1 DO BEGIN
	            data = gain(j) * _wgt(j) * _f(j,k+n_cmn)
		    dat = FMD_PIX_SUM(pixel=pixel(j),data=data,cmp_px=cmp_px)
                    PIX2XY,cmp_px,data=dat,res=6,/six,ras=ras
                    rect(*,k0+k) = PIX2DAT(pix=nz_pix,ras=ras)
	        ENDFOR
	    ENDFOR
	    ; For all channels
	    ;     Find sky observations for that channel
	    ;	  Compute RECT row index offset
 	    ;	  For all non-common kernel functions
	    ;	      Compute: G(c) * wgt * F
	    ;	      Sum over multiple obs within pixel
	    ;	      Build six-pack raster
	    ;	      Store in corresponding rows of RECT
	    ;	  Endfor
	    ; Endfor
	    ;
	    ; Note: CMP_PX is a compressed pixel list.
	    ;	    If PIXEL = [0,1,1,3,4,4,4,6] then 
	    ;	    CMP_PX =   [0,1,3,4,6]
	    ;	    If DATA =  [2,5,1,3,2,5,4,1] then 
	    ;	    DAT =      [2,6,3,11,1]
	    ;
	    ;
	    ; Not all pixels are observed by all channels and some
	    ; pixels are not observed by any channel. The columns
	    ; of RECT correspond to a sorted pixel list including
	    ; all pixels observed by at least one channel.  The
	    ; PIX2XY and PIX2DAT procedures place the pixel-summed
	    ; data array DAT into the correct columns of RECT.
	    ;
	    ; If CMP_PX1 = [0,1,3,4,6] and CMP_PX2 = [0,2,4,7]
	    ; and DAT1 = [2,6,3,11,1] and DAT2 = [1,9,2,6]
	    ; then RECT = [[2,6,0,3,11,1,0],
	    ;	          [1,0,9,0, 2,0,6]]
	    ;
	    ; Note that there is no column in RECT for unobserved
	    ; pixel number 5.


	    ; Compute DIAG (= G * G * wgt)
	    ; ---------------------------------
	    data = (gain_sky * _wgt)^2
	    dat = FMD_PIX_SUM(pixel=pixel(good),data=data(good),cmp_px=cmp_px)
	    PIX2XY,cmp_px,data=dat,res=6,/six,ras=ras
	    diag = PIX2DAT(pix=nz_pix,ras=ras)


	    ; Compute Q INVERSE
	    ; -----------------
            
            ; D INVERSE
	    d_inv = 1 / (diag # one) 

            ; Zero out masked pixels
            dat = FMD_PIX_SUM(pixel=pixel(good),data=sky_mask(good))
            nog_pix = WHERE(dat LE 0)
            IF (nog_pix(0) NE -1) THEN d_inv(nog_pix,*) = 0

	    q0 = square - (TRANSPOSE(rect) # (d_inv*rect))   ;  Q

            q0 = 0.5 * (q0 + TRANSPOSE(q0))                  ; Symmetrize Q

            q_inv = INVERT(q0,q_stat)                        ; Q INVERSE

            IF (q_stat NE 0) THEN $
                PRINT,'FMD_DESTRIPER : Q_STAT =' + STRCOMPRESS(STRING(q_stat))

        ENDIF


        ; Compute B_A (= [S + G * R] * G * wgt)
        ; -------------------------------------
        data = (_spc + (gain_sky * _ref)) * (gain_sky * _wgt)
        dat = FMD_PIX_SUM(pixel=pixel(good),data=data(good),cmp_px=cmp_px)
        PIX2XY,cmp_px,data=dat,res=6,/six,ras=ras
        b_a = PIX2DAT(pix=nz_pix,ras=ras)


        ; Compute B_F
        ; -----------
        b_f = 0d0    ; Build "stub" for B_F

        ; Common Kernels
        FOR i=0,n_cmn-1 DO BEGIN
            b_f0 = _spc - (gain * _pl_ref)
	    b_f0 = TOTAL((b_f0 * mask) * _f(*,i))
	    b_f = [b_f,b_f0]
        ENDFOR

        ; UnCommon Kernels
        IF (n_uncmn GE 1) THEN BEGIN
 	    FOR i=0,n_chan-1 DO BEGIN
	        j = WHERE(idx EQ i)
	        b_f0 = (_spc(j) - (gain(j) * _pl_ref(j)))*mask(j)
;	        b_f0 = (b_f0 * mask(j)) # _f(j,n_cmn:*)
	        ff = _f(j,n_cmn:*)
                IF (n_uncmn GT 1) THEN b_f0 = REFORM(b_f0 # ff)
                IF (n_uncmn EQ 1) THEN b_f0 = TOTAL(b_f0 * ff)
                b_f = [b_f,b_f0]
	    ENDFOR

	    ; For all channels
	    ; 	  Find observations for that channel
	    ;	  Compute (S - G * [bP - R]) * wgt * func
	    ; Endfor

        ENDIF
;

        b_f = b_f(1:*)	; Cut off "stub" 


        ; Compute correction and pixel spectra
        ; ------------------------------------
        ej = REFORM((b_f - b_a # (d_inv*rect)) # q_inv)
        big_a = (b_a - (rect # ej)) / diag


        ; Compute _APR [= (aA - bP - R) * SQRT(wgt)]
        ; ------------------------------------------
        PIX2XY,cmp_px,data=big_a,res=6,ras=ras
        _a = [PIX2DAT(pixel=pixel,ras=ras) * _wgt,DBLARR(n_cal)]

        _apr = _a + _pl_ref


        ; Compute _S_CORR [= (S - SUM[E*u] - SUM[J*T]) * SQRT(wgt)]
        ; ---------------------------------------------------------
        _s_corr = _spc

        IF (n_cmn GE 1) THEN $
            _s_corr = _s_corr - _f(*,0:n_cmn-1) # ej(0:n_cmn-1)

        FOR i=0,n_chan-1 DO BEGIN
            j = WHERE(idx EQ i)
	    k0 = i * (n_krnl - n_cmn) + n_cmn
	    _s_corr(j) = _s_corr(j) - _f(j,n_cmn:*) # ej(k0:k0+(n_krnl-n_cmn)-1)
        ENDFOR

        ; Compute G (= [S_CORR * APR] / APR^2)
        ; ------------------------------------
        IF (NOT KEYWORD_SET(pure_off)) THEN BEGIN
            FOR i=0,n_chan-1 DO BEGIN
	        j = WHERE(idx EQ i)
	        _apr_m = _apr(j) * mask(j)
	        numer = TOTAL(_apr_m * _s_corr(j))
	        denom = TOTAL(_apr_m^2)
	        chan_gain(i) = numer / denom
	    ENDFOR
	    gain = chan_gain(idx)
        ENDIF

        ; Compute New Chi^2
        ; -----------------		
        resid = _s_corr - gain * _apr

        chi2 = TOTAL(resid^2)

        n_iter(l) = n_iter(l) + 1   ; Increment number of iterations

        IF (NOT KEYWORD_SET(pure_off)) THEN $
     	    PRINT,n_iter(l),ABS(chi2_last-chi2)/chi2_last
	 	  ; Print relative change in chi-squared

    ENDREP UNTIL (ABS(chi2_last-chi2)/chi2_last LT convg/(2*ndf) OR $
		      n_iter(l) EQ max_iter OR KEYWORD_SET(pure_off))

    ; Exit gain iteration if
    ; (1) relative change in chi^2 less than 
    ;     convergence criterion, or
    ; (2) # of iterations EQ max # of iterations
    ; (3) Pure Offset Destriping
    ;
    ; Note: If del(chi^2) < convg and chi^2 ~ ndf then
    ;       del(chi^2) / chi^2 ~ convg / ndf.
    ;       Since chi^2 is usually > ndf we use 2 * ndf.

    ejg(*,l) = [ej,chan_gain]      ; Convert and store correction spectra

    pix_spec(*,l) = big_a + mean_sp
          ; Store pixel spectra (with Mean Sky spectrum added back in)



    ; Compute C-Vector from sky observations outside CVEC_MASK
    ; --------------------------------------------------------
    rsd = resid(goodc)^2. * cvec_mask(goodc) / (1. - frac_wgt(goodc))
    c_vec(l) = SQRT(TOTAL(rsd) / ndf_cvec)
    ;

    ; Compute Reduced Chi-Squared per Frequency for each CHANSCAN
    ; -----------------------------------------------------------
    rsd = resid(goodf)^2. * sky_mask(goodf) / (1. - frac_wgt(goodf))
    FOR i=0,n_chan-1 DO BEGIN
	  j = WHERE(sky_idx EQ i)
	  chi2_chan(l,i) = TOTAL((rsd(j)/c_vec(l))^2) / ndf_chan(i)
    ENDFOR
    ;

    ; Compute Total Chi-Squared for each Pixel
    ; ----------------------------------------
    rsd = resid(goodf)^2. / (1. - frac_wgt(goodf))
    chi2_pix = chi2_pix + $
	   FMD_PIX_SUM(pixel=pixel(goodf),data=rsd) / c_vec(l)^2

ENDFOR
	; Frequency Loop



; Rasterize Pixel Spectra
; -----------------------
PIX2XY,/six,nz_pix,data=pix_spec,res=6,raster=afp
;

; Rasterize Number of Good Coadds per Pixel
; -----------------------------------------
PIX2XY,/six,nzf_pix,data=np,res=6,raster=num_coadd
;

; Rasterize Pixel Chi-Squared
; ---------------------------
PIX2XY,/six,nzf_pix,data=chi2_pix,res=6,raster=chi2_map
;

; PCVR in Orthogonal Basis
; ------------------------
pcvr = q_inv
;

; Transpose EJG -> (n_freq,n_func)
; -----------------------------------
ejg = TRANSPOSE(ejg)
;

; Transpose FUNC
; --------------
func = TRANSPOSE(func)
;

; Make SKY_FUNC and CAL_FUNC Arrays
; ---------------------------------
n_sky = N_ELEMENTS(sky_idx)
n_cal = N_ELEMENTS(cal_idx)
n_obs = n_sky + n_cal
sky_func = func(*,0:long(n_sky)-1)
cal_func = func(*,n_sky:n_obs-1)
;

RETURN,c_vec
END
