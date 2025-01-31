PRO fds_errmat,pcvr=pcvr,rect=rect,diag=diag,beta=beta,omega=omega, $
                strp_conv=strp_conv

; This procedure builds the auxillary error matrices, BETA and OMEGA.
; From these, the full pixel covariance and weights matrices can be
; generated.  It also generates the conversion matrix from the old
; stripe basis to the new (beta) stripe basis.
;
; Cov = (1/DIAG) + BETA # TRANSPOSE(BETA)
;
; Wgt = DIAG - OMEGA # TRANSPOSE(OMEGA)
;
; BETA = (1/DIAG) # RECT # STRP_CONV
;
; Written by Joel Gales, ARC, Sep 1994
;
; This routine calls FDS_CHOLESKY to compute the "square root" of
; certain symmetric matrices.  It is written by R. Shafer (GSFC)
;

sz = SIZE(pcvr)
n_stripes = sz(1)

iden = INTARR(n_stripes,n_stripes)
d_elem = INDGEN(n_stripes) * (n_stripes + 1)
iden(d_elem) = 1
	; Build Identity Matrix

one = INTARR(n_stripes) + 1
d_1 = diag # one
d_inv = 1 / d_1



; Omega
; -----
mat = SQRT(d_inv) * rect
qrdcmp,mat,umat,tmat

bmat = iden - INVERT(iden + (tmat # pcvr # TRANSPOSE(tmat)))

lmat = fds_choleskey(bmat)
	; Decompose BMAT

omega = (SQRT(d_1)*umat) # lmat



; Beta
; ----
strp_conv = fds_choleskey(pcvr)
	; Decompose PCVR

beta = (d_inv*rect) # strp_conv


RETURN
END
