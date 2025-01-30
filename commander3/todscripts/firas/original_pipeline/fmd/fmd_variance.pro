PRO FMD_VARIANCE,arc_file,id,splen,d,ntot,f,var,var_tm,cal_var,$
                 cal_var_tm,st_sub,galcut=galcut
;
;+NAME/BRIEF DESCRIPTION OF ROUTINE:
;     FMD_VARIANCE reads in the coadd variances and computes the d-vector
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
;     Modified by Ken Jensen, Hughes STX, 21 Feb 1997, read GMT and
;       save VAR_TM, CAL_VAR_TM .
;     Modified by Ken Jensen, Hughes STX, 14 Mar 1997, Include Secondary
;       Template subtraction in calculation of D-vector .
;     Modified by Ken Jensen, Hughes STX, 11 Apr 1997, Return Secondary
;       Template Subtraction Flag .
;     Modified by Ken Jensen, Hughes STX, 23 May 1997, new algorithm
;       for sorting arrays for compatibility with FFP_CSK files.
;


eject = TIMECONV('893251118',infmt='z',outfmt='s')
		; get cover eject time in seconds


; Open FSL file list
; ------------------

file = STRARR(20)
OPENR,1,'csdr$firas_ref:' + arc_file + '.txt'

n_files = 0
a = ''
WHILE (NOT EOF(1)) DO BEGIN
	READF,1,a
	file(n_files) = a
	n_files = n_files + 1
ENDWHILE
CLOSE,1
file = file(0:n_files-1)
file = 'csdr$firas_insky:fsl_sky_' + id + '.' + file

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

nrec = RSM(file,bufsz,ndim,dimlist,off,type,pixel,gmt,var, $
		nifg,n_adj,st_sub,sky_s0,fnum)

tme = timeconv(gmt,infmt='z',outfmt='s')
var_tm = (tme - eject) / 86400.

; Special Sort by pixel number, file number, record number
; --------------------------------------------------------
PRINT,'Sorting SKY Records'
dummy = 220000L * pixel + 20000L * fnum + start_idx
i = SORT(dummy)
pixel = pixel(i)
var_tm = var_tm(i)
nifg = nifg(i)
n_adj = n_adj(i)
var = var(*,i)
st_sub = st_sub(i)
sky_s0 = sky_s0(i)


; Compute D-Vector
; ----------------

; Index of Good Coadds for D-Vector
; ---------------------------------
;
vt = TOTAL(var,1) / n_freq
;
ll = COORCONV(pixel,infmt='p',inco='f',outfmt='l',outco='g')
;
IF (STRMID(STRUPCASE(id),1,3) EQ 'LFL' OR $
    STRMID(STRUPCASE(id),1,3) EQ 'LFS') THEN begin
     i = WHERE(vt GT -9998. AND ABS(ll(*,1)) GE galcut AND ABS(sky_s0) GE 1) 
ENDIF ELSE BEGIN
  IF (STRMID(STRUPCASE(id),1,3) EQ 'HLF' OR $
    STRMID(STRUPCASE(id),1,3) EQ 'HSF') THEN $
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

        tme = TIMECONV(gmt,infmt='z',outfmt='s')

        obs_time = (tme - eject) / 86400.

	IF (j+k EQ 1) THEN cal_var = cvar
        IF (j+k EQ 1) THEN cal_var_tm = obs_time

	IF (j+k GT 1) THEN cal_var = [[cal_var],[cvar]]
        IF (j+k GT 1) THEN cal_var_tm = [cal_var_tm,obs_time]

    ENDFOR

ENDFOR


RETURN
END
