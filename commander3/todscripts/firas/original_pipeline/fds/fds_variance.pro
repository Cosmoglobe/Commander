PRO fds_variance,arc_file,input_dir,id,var,d,ntot,galcut
;+NAME/BRIEF DESCRIPTION OF ROUTINE:
;     FDS_VARIANCE reads in the coad variances and computes the d-vector
;
;MODIFICATION HISTORY:
;     Written by Joel Gales, ARC, May 1993
;     Modified by Joel Gales, ARC, July 1994
;	Exclude galactic pixels from d-vector calculation (SPR 11859)


; Open FCF file list
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
file = 'fcf_sky_' + id + '.' + file
file = STRUPCASE(file)




IF (STRUPCASE(STRMID(id,1,1)) EQ 'H') THEN BEGIN
	bins = '5:171'
	n_freq = 167
ENDIF ELSE IF (STRUPCASE(STRMID(id,1,2)) EQ 'LS') THEN BEGIN
	bins = '5:38'
	n_freq = 34
ENDIF ELSE IF (STRUPCASE(STRMID(id,2,2)) EQ 'FL') THEN BEGIN
	bins = '5:38'
	n_freq = 34
ENDIF ELSE IF (STRUPCASE(STRMID(id,1,2)) EQ 'LL') THEN BEGIN
	bins = '9:156'
	n_freq = 148
ENDIF



; Load Output Arrays
; ------------------

flds = 'ct_head.gmt'
flds = flds + ',spec_data.real_var(' + bins + ')'
flds = flds + ',coad_spec_head.num_ifgs'
flds = flds + ',coad_spec_data.sec_template.subtracted'

FOR j=0,n_files-1 DO BEGIN

	PRINT,'Reading ',file(j)
	stat = read_skymap(input_dir+file(j),flds,pix,gmt,v,n,sub)

	IF (j EQ 0) THEN BEGIN

		pixel = pix
		nifgs = n
		var = v
		st_sub = sub

	ENDIF ELSE BEGIN

		pixel = [pixel,pix]
		nifgs = [nifgs,n]
		var = [[var],[v]]
		st_sub = [st_sub,sub]

	ENDELSE

ENDFOR



i = SORT(pixel)
pixel = pixel(i)
nifgs = nifgs(i)
var = var(*,i)
st_sub = st_sub(i)

ll = coorconv(pixel,infmt='p',inco='f',outfmt='l',outco='g')
i = WHERE(ABS(ll(*,1)) GE galcut)

d = FLOAT(SQRT(var(*,i)#(nifgs(i)*(nifgs(i)-1-st_sub(i))) / $
		TOTAL(nifgs(i)-1-st_sub(i))))
	; Compute d-vector excluding galactic pixels

ntot = TOTAL(nifgs(i))

RETURN
END
