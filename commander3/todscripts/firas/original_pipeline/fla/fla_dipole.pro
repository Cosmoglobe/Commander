;______________________________________________________________________________
;
;+NAME/ONE-LINE DESCRIPTION:
;    FLA_DIPOLE subtracts a dipole term from a scalar sky map.
;
; DESCRIPTION:
;    An input intensity map (with optional weights) is the basis of
;    a least-squares fit for the position, monopole value, and dipole
;    amplitude of the CMBR.  The routine returns these values as well 
;    as a residual pixel list with the dipole subtracted.  The user can also
;    specify a range of galactic latitude to exclude from the fit.
;
; CALLING SEQUENCE:
;    FLA_DIPOLE, intemps, glon, glat, $
;                Tcmbr, Tcmbr_sig, diamp, diamp_sig, diglon, diglat, dig_sig, $
;                [weights=weights], [galexc=...], [badval=badval], [resid=resid]
;
; ARGUMENTS (I=input, O=output, []=optional):
;    intemps      I    flt arr   pixel list of input temperatures
;    glon         I    flt arr   galactic longitudes of input temperatures
;    glat         I    flt arr   galactic latitudes of input temperatures
;    Tcmbr        O    dble      Mean CMBR temperature 
;    Tcmbr_sig    O    dble      Sigma of CMBR temperature
;    diamp        O    dble      Dipole amplitude
;    diamp_sig    O    dble      Sigma of dipole amplitude
;    diglon       O    dble      Galactic longitude of dipole hot pole (deg)
;    diglat       O    dble      Galactic latitude of dipole hot pole (deg)
;    digl_sig     O    dble      Sigma of dipole direction (deg)
;    [weights]   [I]   flt arr   pixel list of weights of input temperatures,
;                                      default=uniform wts
;    [galexc]    [I]   flt       Absolute value of minimum galactic latitude
;                                      to include in dipole fit (deg),
;                                      default=0.0
;    [badval]    [I]   flt       Flag value for bad pixels, default=0.
;    [resid]     [O]   flt arr   pixel list of residual temperatures
;
; WARNINGS: 
;    1.  It is generally wise to specify a galactic latitude exclusion
;           of 20 degrees or more for believable results.
;    2.  Flagged pixels will show up as zeros in the residual pixel list.
;
; EXAMPLE:
;    If "intemps" is a map of temperatures in which empty pixels have
;    a value of zero then a simple call for an unweighted fit would be
;
;    UIDL> fla_dipole,intemps,glon,glat, $
;                     avtemp,tsig,diamp,ampsig,dlon,dlat,dir_sig,galexc=20
;    UIDL> print,avtemp,tsig,diamp,ampsig,dlon,dlat,dir_sig
;             2.7256  0.000005  0.003345  0.000016  265.60  48.34  0.24
;
;   If, e.g., "resid=residual" had also been included in the call 
;   then the variable "residual" would contain a pixel list of residuals
;   after subtracting a 2.7256 K average and 3.345 mK dipole from intemps.
;
;   If in addition to intemps we have a second pixel list, e.g. "sigmas"
;   that contains the weight at each pixel, then we could also have
;   included ", weights=sigmas" in the call.
;#
; COMMON BLOCKS:  none
;
; PROCEDURE:
;       (1)  Cull pixels with flag values from the list
;       (2)  Cull low galactic latitude pixels from the list
;       (3)  Convert galactic coordinates of remaining pixels to unit vectors
;       (4)  Using unit vectors, do a linear fit to the dipole in
;               Cartesian coordinates, i.e. Model = T0 + Td.u, where
;               the second term is the dot product of the dipole (of
;               amplitude Td) and the unit vector u.
;       (5)  Compute the errors for the fit.  For uniform weights, use the
;               assumption that the chi-square per degree of freedom is 1 to
;               compute sigmas for the fit.
;       (6)  Convert the best-fit unit vector into galactic coords
;       (7)  Compute the sigma for the direction by computing the angle
;               between the dipole and the dipole + sigmas.
;       (8)  Create a dipole pixel list and subtract from the input temperatures
;               if temperature residuals are desired.
;
; LIBRARY CALLS:
;     coorconv
;
; MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  May 1992
;    Modified by Celine Groden, USRA,  December 1992.  SPR # 10327
;    Modified by Gene Eplee, General Sciences Corp, Dec 1992
;        to compute sigmas for fit.
;    Modified by use in FLA by Gene Eplee, GSC, Oct 1994
;    Modified to correct dipole uncertainties by Joel Gales, ARC, Dec 1994
;-
;______________________________________________________________________________
;
PRO fla_dipole, intemps, glon, glat, $
                tcb, sigt, di_amp, sigdp, di_glon, di_glat, sigdpd, $
                weights=weights, galexc=galexc, badval=badval, resid=resid
;
; Set Error Status
;
ON_ERROR, 2
;
;  Set defaults for optional parameters (e.g. uniform weighting) 
;
IF (N_ELEMENTS(badval) EQ 0) THEN badval = 0.0
IF (N_ELEMENTS(galexc) EQ 0) THEN galexc = 0.0                       
IF (N_ELEMENTS(weights) EQ 0) THEN begin
    weights = 0.0*intemps + 1.0
    x2=1
endif else begin
    x2=0
endelse
;        
;  Get the 1-D index of all non-flag data points, and create 
;  corresponding 1-D pixel-ordered lists of good data, weights,
;  and the galactic coordinates.
;
tsz=size(intemps)
galcoin = fltarr(tsz(1),2)
galcoin(*,0) = glon
galcoin(*,1) = glat
gooddata = WHERE (intemps NE badval)
cbtemps = intemps(gooddata)                        ; good data
wtvec = weights(gooddata)                          ; weights
galco = galcoin(gooddata,*)
;
;  Filter out the pixels whose latitude is too low.
;
highlat = WHERE (ABS(glat) GE galexc)
bsz=size(highlat)
IF (bsz(0) EQ 0) THEN BEGIN
   PRINT,'There is no good data above the specified galactic latitude'
   RETURN
ENDIF
galco = galco(highlat,*)
cbtemps = cbtemps(highlat)                        ; culled data
wtvec = wtvec(highlat)                            ; culled weights
avtemp = total(cbtemps)/n_elements(cbtemps)
cbtemps = cbtemps - avtemp                        ; create a zero mean
;
;  Create the unit vector for each pixel, then fill the normal matrix
;  QMAT for the fit. The unit vectors are the derivatives of the model, and
;  the fit equation is  QMAT # DPARMS = RHS
;
uvecs = coorconv (galco, infmt='l', inco='g', outfmt='u', outco='g')
uvecs = DOUBLE(uvecs)
;
qmat = DBLARR(4,4)
qmat(0,0) = TOTAL(wtvec)
qmat(1,0) = TOTAL(uvecs(*,0) * wtvec)
qmat(2,0) = TOTAL(uvecs(*,1) * wtvec)
qmat(3,0) = TOTAL(uvecs(*,2) * wtvec)
qmat(0,1) = qmat(1,0)
qmat(1,1) = TOTAL(uvecs(*,0)^2 * wtvec)
qmat(2,1) = TOTAL(uvecs(*,0) * uvecs(*,1) * wtvec)
qmat(3,1) = TOTAL(uvecs(*,0) * uvecs(*,2) * wtvec)  
qmat(0,2) = qmat(2,0)
qmat(1,2) = qmat(2,1)
qmat(2,2) = TOTAL(uvecs(*,1)^2 * wtvec)
qmat(3,2) = TOTAL(uvecs(*,1) * uvecs(*,2) * wtvec)
qmat(0,3) = qmat(3,0)
qmat(1,3) = qmat(3,1)
qmat(2,3) = qmat(3,2)
qmat(3,3) = TOTAL(uvecs(*,2)^2 * wtvec)
;
rhs = DBLARR(4)
rhs(0) = TOTAL(cbtemps * wtvec)
rhs(1) = TOTAL(cbtemps * uvecs(*,0) * wtvec)
rhs(2) = TOTAL(cbtemps * uvecs(*,1) * wtvec)
rhs(3) = TOTAL(cbtemps * uvecs(*,2) * wtvec)
;
;  DPARMS holds the 4 fit parameters.  Add the average temperature back
;  in to restore the monopole value.
;
invmat = INVERT(qmat)
dparms = invmat # rhs
tcb = dparms(0) + avtemp
;
;  Normalize the cartesian representation of the dipole to get the unit
;  vector of its direction, and convert to galactic coordinates.  The 
;  squared modulus of the three componenets is the dipole amplitude.
;
di_amp = SQRT (TOTAL (dparms(1:3)^2))
di_uvec = dparms(1:3)/di_amp
galcoords = coorconv(di_uvec, infmt='u', outfmt='l')
di_glon = galcoords(0)
di_glat = galcoords(1)
;
;  Compute the errors of the fit parameters.
;    For uniform weights, force the chi-square per degree of freedom of the fit
;    to be 1.
;
chi2 = total(wtvec*((dparms(0) + (uvecs # dparms(1:3)) - cbtemps)^2))
dof = double(n_elements(cbtemps)) - 4.0d
x2dof = chi2/dof
variances = dblarr(4)
sigmas = dblarr(4)
if (x2 eq 1) then begin
   for j=0,3 do begin
      variances(j) = invmat(j,j)*x2dof
      sigmas(j) = sqrt(variances(j))
   endfor
endif else begin
   for j=0,3 do begin
      variances(j) = invmat(j,j)
      sigmas(j) = sqrt(variances(j))
   endfor
endelse
;
;  Error in monopole amplitude.
;
sigt = sigmas(0)
;
;  Error in dipole amplitude and direction
;
d2r = !dpi / 180

abs_dpl = ABS(dparms(1:3))
	; Get absolute values of dipole parameters

di_adj = abs_dpl + sigmas(1:3)
	; Add dipole uncertainties

amp_adj = SQRT(TOTAL(di_adj^2))
	; Get magnitude of adjusted dipole
 
sigdp = TOTAL(abs_dpl * sigmas(1:3)) / di_amp
	; Get uncertainty in dipole magnitude

sigdpd = ACOS(TOTAL(di_adj * abs_dpl) / (di_amp * amp_adj)) / d2r
	; Get uncertainty in dipole direction
;
;  We now have all fit parameters, so create an all-sky map holding only
;  the dipole.  First create the unit vector for each input pixel.
;
alluvecs = coorconv (galcoin, infmt='l', inco='g', outfmt='u', outco='g')
;
;  Now create a model sky by taking the dot product of each pixel's
;  unit vector with the dipole, and adding in the monopole.
;
outtemps = intemps * 0.0
for j=0,tsz(1)-1 do outtemps(j) = tcb + total(alluvecs(j,*) * dparms(1:3))
;
;  Flag as bad any pixels in the model sky that are bad in the original data.
;  Then difference it with the original data to create a residual pixel list.
;  (Residuals will be zero at flagged points.)
;
if ((where(intemps eq badval)) ne -1) then $
     outtemps(where(intemps eq badval)) = badval
resid = intemps - outtemps
;
return
end
