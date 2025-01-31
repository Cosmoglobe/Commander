PRO FMD_errmat,pcvr=pcvr,rect=rect,diag=diag,beta=beta,omega=omega, $
                stripe_conv=stripe_conv,ejg=ejg,gamma=gamma,lmat=lmat

; This procedure builds the auxillary error matrices, BETA and OMEGA.
; From these, the full pixel covariance and weights matrices can be
; generated.  It also generates the conversion matrix from the old
; stripe basis to the new (beta) stripe basis, and the orthogonalized
; stripes.
;
; Cov = (1/DIAG) + BETA # TRANSPOSE(BETA)
;
; Wgt = DIAG - OMEGA # TRANSPOSE(OMEGA)
;
; BETA = (1/DIAG) # RECT # STRIPE_CONV
;
; GAMMA = INVERSE(STRIPE_CONV) # EJG
;
;
; Written by Joel Gales, ARC, Sep 1994
;
; Renamed FMD_ERRMAT, K.Jensen, Hughes STX, 21-Jun-96
; Calculation of GAMMA,  K.Jensen,  21-Jun-96
; Return LMAT, K.Jensen, 3-Sep-96
;
; This routine calls FMD_CHOLESKY to compute the "square root" of
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
QRDCMP,mat,umat,tmat

bmat = iden - INVERT(iden + (tmat # pcvr # TRANSPOSE(tmat)))

lmat = FMD_CHOLESKEY(bmat)
	; Decompose BMAT

omega = (SQRT(d_1)*umat) # lmat



; Beta
; ----
stripe_conv = FMD_CHOLESKEY(pcvr)
	; Decompose PCVR

beta = (d_inv*rect) # stripe_conv

;
kk = INVERT(stripe_conv)
;
; kk # ejg = gamma
;

sz = SIZE(ejg)
n_freq = sz(1)

gammax = FLTARR(n_freq,n_stripes)
idx = INDGEN(n_stripes)
ejg0 = dblarr(n_stripes)
;
FOR i=0,n_freq-1 DO BEGIN
 ejg0(*) = ejg(i,idx(*))
 ejgg = kk # ejg0
 gammax(i,*) = ejgg(*)
ENDFOR
;

gamma = gammax
;

RETURN
END
