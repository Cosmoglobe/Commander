Pro FLA_Fit_LHR,chanscan,error
;
; 1) Calls FLA_GET_TCBR
; 2) Computes CBR-subtracted spectra to be fitted for lines and dust.
;
error = 1
;
if N_Params() ne 2 then begin
 print,'FLA_Fit_LHR,chanscan,error'
 return
endif
;

; CHANSCAN and IDD
; ----------------
chanscan = strupcase(chanscan)
idd=strmid(chanscan,0,1)+'LF'
if((chanscan ne 'RLLF')and(chanscan ne 'LLLF'))then begin
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
if(n_freq ne 148)then begin
 print,'Error in Frequency Array !'
 return
endif
;
RESTORE,'csdr$firas_ds:'+idd+'_destriped.iss'
csp = 0
;

if(idd eq 'RLF')then begin
 cxxx=crlf & s_xxx=s_rlf & b_xxx=b_rlf & l_xxx=l_rlf & n_xxx=n_rlf
 tau_rlf=0 & ejv_rlf=0 & d_rlf=0 & urlf=0
endif
if(idd eq 'LLF')then begin
 cxxx=cllf & s_xxx=s_llf & b_xxx=b_llf & l_xxx=l_llf & n_xxx=n_llf
 tau_llf=0 & ejv_llf=0 & d_llf=0 & ullf=0
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
n_base = 4
;

; Line Frequencies (in icm)
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

; Set line profile frequencies & generate line profiles
; -----------------------------------------------------
n_lines = N_ELEMENTS(line_freq)
;
line_pro = COMPLEXARR(n_freq,n_lines)
;
FOR i=0,n_lines-1 DO $
	line_pro(*,i) = FLA_FIR_LINE(line_freq(i),mode,fex_gain_ext, $
				     gauss=1./156,sinc=0)
;

; Generate Legendre Polynomial Fitting Vectors
; --------------------------------------------
PRINT,'Generating Legendre Polynomial Fitting Vectors'
;
base_line = DBLARR(n_freq,n_base+1)
;
avg_freq = (MAX(freq) + MIN(freq)) / 2
f_scl = MAX(freq) - avg_freq
;
x = (freq - avg_freq) / f_scl
base_line(*,0) = 1
base_line(*,1) = x
;
FOR i=1,n_base-1 DO base_line(*,i+1)= $
	((2*i+1) * x * base_line(*,i) - i * base_line(*,i-1))/(i+1)
;
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

; Initial Guesses for Dust Fit
; ----------------------------
alpha = 1.55                ; Dust Emissivity Index
f_char = 1800./29.9792458   ; Characteristic Frequency (icm)
;

; Initial Guesses for Dust Model
; ------------------------------
parm0 = [0,18.,0.1,alpha,0,0,0]
;

; Active Parameter Flags
; ----------------------
act = [0,0,1,0,0,0,0]   ; Fix Temperature, Fit Optical Depth
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

; FIP_LML_xxxx binary file
; ------------------------
print,'Writing FIP_LML_'+chanscan+'.F16_93HYBRID'
FLA_FIP_LML,b_xxx=b_xxx,l_xxx=l_xxx,n_xxx=n_xxx,covar=pcvr, $
   	    line_parms=lfx,chanscan=chanscan,n_lines=n_lines, $
            line_idx=line_idx
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
sl_weight_rat=0 & tcovar=0
;

print,' '
error = 0
;

RETURN
END
