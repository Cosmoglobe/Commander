FUNCTION fla_get_tcbr,freq,pix,num,l_xxx,b_xxx,s_xxx,cxxx


RESTORE,'csdr$firas_restore:gal_tmpl.iss'
	; Restore galactic template


d2r = !dpi / 180
glon = pix2dat(pix=pix,ras=l_xxx)
glat = pix2dat(pix=pix,ras=b_xxx)
pow = pix2dat(pix=pix,ras=gal_tmpl)
	; Extract galactic longitude/latitude and gal template for
	; non-empty pixels


sqrt_w = (ABS(glat) GT 15) * SQRT(num)
	; Exclude pixels within 15 deg of galactic plane


spec = pix2dat(pix=pix,ras=DOUBLE(s_xxx))
	; Extract spectra 


spc_w = spec
FOR i=0,N_ELEMENTS(freq)-1 DO spc_w(*,i) = spc_w(*,i) * sqrt_w
	; Compute weighted spectra

spec = TRANSPOSE(spec)
spc_w = TRANSPOSE(spc_w)
	; Transpose to make frequency first index

unf = sqrt_w
dx = COS(glat * d2r) * COS(glon * d2r) * sqrt_w
dy = COS(glat * d2r) * SIN(glon * d2r) * sqrt_w
dz = SIN(glat * d2r) * sqrt_w
template = pow * sqrt_w
	; Build (weighted) uniform, dipole, and galactic templates




; Compute Galactic Spectrum from Least-Squared Spectral Fit
; ---------------------------------------------------------
func = [[unf],[dx],[dy],[dz],[template]]

alpha = TRANSPOSE(func) # func
covar = INVERT(alpha)
beta = spc_w # func

spc = beta # covar
gal_spc = spc(*,4)



; Compute CBR Temperature
; -----------------------
pl = planck(2.726d0,freq,dpdt,unit='i',/mjy) * 2.99792458d-7
dpdt = dpdt * 2.99792458d-7
	; Get Planck and dP/dT

rspec = 0 * spec
sz = SIZE(rspec)
FOR i=0,sz(2)-1 DO rspec(*,i) = (spec(*,i) - pl) / cxxx(*,0)
	; Compute (weighted) CBR-subtracted spectra



; Fit to Planck, dP/dT, and Galactic spectrum
; -------------------------------------------
fnc = [[dpdt/cxxx(*,0)],[gal_spc/cxxx(*,0)]]

alpha = TRANSPOSE(fnc) # fnc
covar = INVERT(alpha)
beta = TRANSPOSE(fnc) # rspec
out = covar # beta

tcbr = 2.726d0 + REFORM(out(0,*))
unc_tcbr = SQRT(covar(0,0) / num)
	; Compute CBR temp and uncertainty

RETURN,[[tcbr],[unc_tcbr]]
END
