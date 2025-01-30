PRO FMD_variance2,arc_file,id,splen,f,var,var_tm,var_lbl,cal_var,$
        cal_var_tm,cal_var_lbl,st_sub,galcut=galcut
;
;+NAME/BRIEF DESCRIPTION OF ROUTINE:
;     FMD_VARIANCE2 reads in the coadd variances for Long and Short
;
;MODIFICATION HISTORY:
;     Written by Joel Gales, ARC, May 1993, as FDS_VARIANCE.
;     Modified by Joel Gales, ARC, July 1994
;	Exclude galactic pixels from d-vector calculation (SPR 11859)
;     Modified by Joel Gales, ARC, Dec 1995
;       Replace calls to READ_SKY_MAP by RSM (C linkimage)
;     Modified by Ken Jensen, Hughes STX, Dec 1995, renamed FMD_VARIANCE.
;     Modified by Ken Jensen, Hughes STX, 28 Feb 1996, frequency vector
;       computed and returned.
;     Modified by Ken Jensen, Hughes STX, 23 Jan 1997, eliminate flagged
;       variances prior to D-vector computation
;     Modified by Ken Jensen, Hughes STX, 24 Jan 1997, read Calibration
;       variances
;     Modified by Ken Jensen, Hughes STX, 20 Feb 1997, and renamed
;       FMD_VARIANCE2,  concatenates Long and Short BEFORE pixel sort,
;       for consistency with FMD_READ_DATA.
;     Modified by Ken Jensen, Hughes STX, 21 Feb 1997, read GMT and
;       save VAR_TM , CAL_VAR_TM .
;     Modified by Ken Jensen, Hughes STX, 14 Mar 1997, include
;       secondary template subtraction in calculation of D-Vector
;     Modified by Ken Jensen, Hughes STX, 11 Apr 1997, return
;       secondary template subtraction flag
;     Modified by Ken Jensen, Hughes STX, 27 Apr 1997, new algorithm
;       for sorting the arrays, for compatibility with FFP_CSK files.
;


eject = TIMECONV('893251118',infmt='z',outfmt='s')
		; get cover eject time in seconds

; Open FSL file list
; ------------------

file = STRARR(20)
;

OPENR,1,'csdr$firas_ref:' + arc_file(0) + '.txt'

n_files = 0
a = ''
WHILE (NOT EOF(1)) DO BEGIN
	READF,1,a
	file(n_files) = a
	n_files = n_files + 1
ENDWHILE
CLOSE,1
;
file = file(0:n_files-1)
file = 'csdr$firas_insky:fsl_sky_' + id(0) + '.' + file
;

file2 = STRARR(20)
;

OPENR,1,'csdr$firas_ref:' + arc_file(1) + '.txt'

n_files = 0
a = ''
WHILE (NOT EOF(1)) DO BEGIN
	READF,1,a
	file2(n_files) = a
	n_files = n_files + 1
ENDWHILE
CLOSE,1
;
file2 = file2(0:n_files-1)
file2 = 'csdr$firas_insky:fsl_sky_' + id(1) + '.' + file2
;

;file = STRUPCASE([file,file2])
;

; Determine calibration model and nyquist frequency
; -------------------------------------------------
print,'Reading calibration model solution and nyquist frequency'
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


; Read Data from FSL skymap files
; -------------------------------
bufsz = 11264L
ndim = [0,0,1,0,0,0,0]
dimlist = [14L,n_freq]
off = [11184L,0,5548L+nzero(0)*4,132,134,1681,2632]
type = [3,7,4,2,4,1,4]

nrec = RSM(file,bufsz,ndim,dimlist,off,type,pixel0,gmt,var0, $
		nifg0,n_adj0,st_sub0,sky_s00,fnum)

n_files = N_ELEMENTS(file)
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

dummy = 220000L * pixel0 + 20000L * fnum + start_idx
i = SORT(dummy)

var_lbl = STRMID(file(fnum(i)),STRPOS(file(0),'.ED')-4,4)

tme = timeconv(gmt(i),infmt='z',outfmt='s')
var_tm = (tme - eject) / 86400.

pixel = pixel0(i)
nifg = nifg0(i)
n_adj = n_adj0(i)
var = var0(*,i)
st_sub = st_sub0(i)
sky_s0 = sky_s00(i)

nrec = RSM(file2,bufsz,ndim,dimlist,off,type,pixel0,gmt,var0, $
		nifg0,n_adj0,st_sub0,sky_s00,fnum)


n_files = N_ELEMENTS(file2)
;
OPENR,2,file2(0)
finfo = FSTAT(2)
nrec0 = finfo.size / finfo.rec_len
start_idx = INDGEN(nrec0)
CLOSE,2
;
FOR i=1,n_files-1 DO BEGIN
 OPENR,2,file2(i)
 finfo = FSTAT(2)
 nrec0 = finfo.size / finfo.rec_len
 start_idx = [start_idx,INDGEN(nrec0)]
 CLOSE,2
ENDFOR
;

dummy = 220000L * pixel0 + 20000L * fnum + start_idx
i = SORT(dummy)

var_lbl = [var_lbl,STRMID(file2(fnum(i)),STRPOS(file2(0),'.ED')-4,4)]

tme = timeconv(gmt(i),infmt='z',outfmt='s')
var_tm = [var_tm,(tme - eject) / 86400.]

pixel = [pixel,pixel0(i)]
nifg = [nifg,nifg0(i)]
n_adj = [n_adj,n_adj0(i)]
var = [[var],[var0(*,i)]]
st_sub = [st_sub,st_sub0(i)]
sky_s0 = [sky_s0,sky_s00(i)]


; Compute D-Vector
; ----------------

; Index of Good Coadds for D-Vector
; ---------------------------------
;
vt = TOTAL(var,1) / n_freq
;
ll = COORCONV(pixel,infmt='p',inco='f',outfmt='l',outco='g')
;
IF (STRMID(STRUPCASE(id(0)),1,3) EQ 'LFL' OR $
    STRMID(STRUPCASE(id(0)),1,3) EQ 'LFS') THEN begin
     i = WHERE(vt GT -9998. AND ABS(ll(*,1)) GE galcut AND ABS(sky_s0) GE 1) 
ENDIF ELSE BEGIN
  IF (STRMID(STRUPCASE(id(0)),1,3) EQ 'HLF' OR $
    STRMID(STRUPCASE(id(0)),1,3) EQ 'HSF') THEN $
     i = WHERE(vt GT -9998. AND ABS(ll(*,1)) GE galcut AND ABS(sky_s0) GE .23) $
  ELSE $
	i = WHERE(vt GT -9998. AND ABS(ll(*,1)) GE galcut)
ENDELSE


; D-Vector, excluding galactic pixels and default coadds
; ------------------------------------------------------
num = var(*,i) # ((nifg(i)-1-st_sub(i))*n_adj(i))
den = TOTAL(nifg(i)-1-st_sub(i))
d = FLOAT(SQRT(num/den))

; Total Weight in D-Vector Calculation
; ------------------------------------
ntot = TOTAL(n_adj(i))


; Read Variances from fsl cal files
; ---------------------------------
PRINT,' '
PRINT,'Reading Cal Data'

bufsz = 11264L
ndim = [0,1]
dimlist = [14L,n_freq]
off = [0,5548L+nzero(0)*4]
type = [7,4]

FOR j=0,N_ELEMENTS(id)-1 DO BEGIN

    FOR k=1,7 DO BEGIN

	cfile = 'csdr$firas_incal:fsl_cal_' + id(j) + '.dihed_' $
		+ STRING(k)
	cfile = STRCOMPRESS(cfile,/REMOVE_ALL)

        nrec = RSM(cfile,bufsz,ndim,dimlist,off,type,gmt,cvar)

        tme = timeconv(gmt,infmt='z',outfmt='s')
        obs_time = (tme - eject) / 86400.

	IF (j+k EQ 1) THEN BEGIN
            cal_var = cvar
            cal_var_lbl = REPLICATE(id(j),nrec)
            cal_var_tm = obs_time
        ENDIF

	IF (j+k GT 1) THEN BEGIN
            cal_var = [[cal_var],[cvar]]
            cal_var_lbl = [cal_var_lbl,REPLICATE(id(j),nrec)]
            cal_var_tm = [cal_var_tm,obs_time]
        ENDIF

    ENDFOR

ENDFOR

cal_var_lbl = STRUPCASE(cal_var_lbl)


RETURN
END
