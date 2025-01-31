Pro FLA_Fit_Hi,chanscan,error
;
; 1) Calls FLA_GET_TCBR
; 2) Computes CBR-subtracted spectra to be fitted for lines and dust.
;
error = 1
;
if N_Params() ne 2 then begin
 print,'FLA_Fit_Hi,chanscan,error'
 return
endif
;

; CHANSCAN and IDD
; ----------------
chanscan = strupcase(chanscan)
idd = strmid(chanscan,0,2) + strmid(chanscan,3,1)
if((idd ne 'RHS')and(idd ne 'RHF')and(idd ne 'LHS')and(idd ne 'LHF'))then begin
 print,'Error in CHANSCAN ID !'
 return
endif
;

; Restore DeStriper Save Sets
; ---------------------------
RESTORE,'csdr$firas_raw:'+idd+'.iss'
sp = 0
freq = f
n_freq=n_elements(freq)
if(n_freq ne 167)then begin
 print,'Error in Frequency Array !'
 return
endif
;
RESTORE,'csdr$firas_ds:'+idd+'_destriped.iss'
csp = 0
;

if(idd eq 'RHS')then begin
 cxxx=crhs & s_xxx=s_rhs & b_xxx=b_rhs & l_xxx=l_rhs & n_xxx=n_rhs
 tau_rhs=0 & ejv_rhs=0 & d_rhs=0 & urhs=0
endif
if(idd eq 'RHF')then begin
 cxxx=crhf & s_xxx=s_rhf & b_xxx=b_rhf & l_xxx=l_rhf & n_xxx=n_rhf
 tau_rhf=0 & ejv_rhf=0 & d_rhf=0 & d_rhsf=0 & d_rhlf=0 & urhf=0
endif
if(idd eq 'LHS')then begin
 cxxx=clhs & s_xxx=s_lhs & b_xxx=b_lhs & l_xxx=l_lhs & n_xxx=n_lhs
 tau_lhs=0 & ejv_lhs=0 & d_lhs=0 & ulhs=0
endif
if(idd eq 'LHF')then begin
 cxxx=clhf & s_xxx=s_lhf & b_xxx=b_lhf & l_xxx=l_lhf & n_xxx=n_lhf
 tau_lhf=0 & ejv_lhf=0 & d_lhf=0 & d_lhsf=0 & d_lhlf=0 & ulhf=0
endif
;

; Determining Non-zero Pixels
; ---------------------------
PRINT,'Determining Non-zero Pixels'
n = pix2dat(pix=INDGEN(6144),raster=n_xxx)
pix = WHERE(n GT 0,npix)
nifg_in_pix = n(pix)
;

; Generate Line & Dust Weights
; ----------------------------
dst_wgt = 1 / DOUBLE(cxxx(*,0))^2
line_wgt = dst_wgt
;

; Extracting Spectra from Destriped Dataset
; -----------------------------------------
PRINT,'Extracting spectra from map'
spec = pix2dat(pixel=pix,raster=s_xxx)
spec = TRANSPOSE(spec)
;

; Set line profile parameters
; ---------------------------
mode = chanscan
fex_gain_ext = 'f16_93hybrid'
;

; Set baseline order for line fits
; --------------------------------
n_base = 10
;

; Line frequencies (in icm)
; -------------------------
co_04  =  3.845
co_08  =  7.690
co_11  = 11.535
o2_14  = 14.168
co_15  = 15.379
ci_16  = 16.419
h2o_18 = 18.576
co_19  = 19.222
;
line_freq = [co_04,co_08,co_11,o2_14,co_15,ci_16,h2o_18,co_19]
;
co_freq = 23.06
c1_freq = 27.00
h2o_freq = 37.136
n2a_freq = 48.738
c2_freq = 63.395
o1_freq = 68.716
si1_freq = 77.11
n2b_freq = 82.036
;
line_freq = [line_freq,co_freq,c1_freq,h2o_freq,n2a_freq, $
             c2_freq,o1_freq,si1_freq,n2b_freq]
;

; Set line profile frequencies & generate line profiles
; -----------------------------------------------------
n_lines = N_ELEMENTS(line_freq)
;
line_pro = COMPLEXARR(n_freq,n_lines)
;
FOR i=0,n_lines-1 DO $
	line_pro(*,i) = fla_fir_line(line_freq(i),mode,fex_gain_ext, $
				     gauss=1./156,sinc=0)
;

; Generate Legendre Polynomial Fitting Vectors
; --------------------------------------------
PRINT,'Generating Legendre Polynomial Fitting Vectors'

base_line = DBLARR(n_freq,n_base+1)

avg_freq = (MAX(freq) + MIN(freq)) / 2
f_scl = MAX(freq) - avg_freq

x = (freq - avg_freq) / f_scl
base_line(*,0) = 1
base_line(*,1) = x

FOR i=1,n_base-1 DO base_line(*,i+1)= $
	((2*i+1) * x * base_line(*,i) - i * base_line(*,i-1))/(i+1)

fit_func = [[DOUBLE(line_pro)],[base_line], $
	    [DOUBLE(IMAGINARY(line_pro))],[0*base_line]]
;


; Compute CBR Temperature and Uncertainty and CBR-subtracted spectra
; ------------------------------------------------------------------
m2ep = 2.99792458d-7
tu = FLA_GET_TCBR(freq,pix,nifg_in_pix,l_xxx,b_xxx,s_xxx,cxxx)
tcbr = tu(*,0)
tcbr_unc = tu(*,1)
;

; Cosmic background subtraction
; -----------------------------
FOR i=0,npix-1 DO $
	spec(*,i) = spec(*,i) - planck(tcbr(i),freq,un='i',/mjy) * m2ep
;

; ZODI subtraction
; ----------------
RESTORE,'csdr$firas_restore:modzodi.iss'
FLA_MATCH,indgen(6144),pix,nz,nx
nnpix = N_ELEMENTS(nz)
FOR i=0,nnpix-1 DO $
        spec(*,nx(i)) = spec(*,nx(i)) - ( modzodi(*,nz(i)) * m2ep )
;


; Initial Guesses for Dust Fit
; ----------------------------
alpha = 1.55                ; Dust Emissivity Index
f_char = 1800./29.9792458   ; Characteristic Frequency (icm)
;

; Initial Guesses for Dust Model
; ------------------------------
parm0 = [0,22.,0.1,alpha,0,0,0]
;

; Active Parameter Flags
; ----------------------
act = [0,1,1,0,0,0,0]   ; Fit Temperature and Optical Depth
;

; Frequency Mask
; --------------
dew = [36,37,43,44,61,62,81,82,83,107,108,109,118,132,141]
;

; Pixel-by-Pixel Linear least-squares Fit
; ---------------------------------------
csqr = FLA_FIT(freq,spec,line_wgt,nifg=nifg_in_pix, $
	       mask=dew,dst_wgt=dst_wgt,dst_sig=dst_sig, $
	       lne_parms=lfx,dst_parms=dst,act=act,parm0=parm0, $
	       fit_func=fit_func,f_char=f_char,n_lines=n_lines, $
	       pcvr=pcvr,dst_covar=dst_covar,dust_spec=dust_spec)
;

; Store CMBR temps and uncertainties
; ---------------------------------
dst(0,*) = tcbr
dst_sig(0,*) = tcbr_unc
;

; IDL Save Set
; ------------
sname='csdr$firas_save:lines_'+idd+'.iss'
SAVE,filename=sname,chanscan,mode,f,lfx,nifg_in_pix,pix,dst,dst_sig,dst_covar, $
                    alpha,f_char,n_base,csqr,fit_func,pcvr,freq,b_xxx,l_xxx
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Restore the background spectrum to the continuum spectrum
; ---------------------------------------------------------
FOR i=0,npix-1 DO $
    dust_spec(*,i) = dust_spec(*,i) + planck(tcbr(i),freq,un='i',/mjy) * m2ep
;

; Restore the ZODI spectrum to the continuum spectrum
; ---------------------------------------------------
FOR i=0,nnpix-1 DO $
    dust_spec(*,nx(i)) = dust_spec(*,nx(i)) + ( modzodi(*,nz(i)) * m2ep )
;

; Generate sixpack dust skympap
; -----------------------------
PIX2XY,pix,data=TRANSPOSE(dust_spec),res=6,/six,ras=dust_map
;

; DUST_xxx IDL Save Set
; ---------------------
sname='csdr$firas_save:dust_'+idd+'.iss'
SAVE,filename=sname,chanscan,f,dust_map
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Line Mask
; ---------
line_idx = [1,1,1,1,1,1,1,1]
;

; FIP_LMH_xxxx binary file
; ------------------------
print,'Writing FIP_LMH_'+chanscan+'.F16_93HYBRID'
FLA_FIP_LMH,b_xxx=b_xxx,l_xxx=l_xxx,n_xxx=n_xxx,covar=pcvr, $
  line_parms=lfx,dst_parms=dst,chanscan=chanscan,dst_sig=dst_sig, $
  dst_covar=dst_covar,n_lines=n_lines,line_idx=line_idx
;

; FLA_DST_xxxx binary file
; ------------------------
print,'Writing FLA_DST_'+chanscan+'.F16_93HYBRID'
FLA_DST_ST,b_xxx=b_xxx,l_xxx=l_xxx,n_xxx=n_xxx,spec=dust_spec,chanscan=chanscan
;

; Re-Define Unused Restored Parameters
; ------------------------------------
px=0 & tm=0 & sp=0 & nifgs=0 & glon=0 & glat=0 & solution=0 & galcut=0
cal_sp=0 & cal_nifgs=0 & cal_tm=0 & xcal=0 & ical=0 & skyh=0 & refh=0 & dihd=0
sky_glitch=0 & cal_glitch=0 & sky_wgts=0 & cal_wgts=0 & del_temp=0 & csp=0
sky_wgts_ds=0 & step_up=0 & step_dn=0 & ndf_ejv=0 & sqr_mat=0 & rect=0 & diag=0
c_calsp=0 & ndf_none=0 & cal_wgts_ds=0 & vib=0 & corr_index=0 & index=0
sl_weight_rat=0 & tcovar=0 & goodzodipix=0 & evec=0 & zodimodel=0
basezeta=0 & basealpha=0 & allzodi=0 & deltaalpha=0 & deltazeta=0
;

print,' '
error = 0
;

RETURN
END
