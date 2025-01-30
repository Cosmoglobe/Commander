      integer*4 function fil_baseline_sub(basis,coa_rec)
c----------------------------------------------------------------------------
c
c     Purpose: Calculate and subtract the baseline from the coadded
c              interferogram. A fourth-order Legendre polynomial is fit to the 
c              last 492 points of the interferogram using a least squares 
c              criterion implemented with the imsl routine "dlsqrr". The 
c              baseline is subtracted from the coadded interferogram and the 
c              result is returned in the coadd record. The coefficients of the 
c              fit are also returned. The first twenty points of the 
c              interferogram are excluded from the fit in order to reduce 
c              possible transient effects.
c
c     Authors: R. Pisarski, ARC, 3/85
c              A. Trenholme, GSC, 9/91
c              S. Brodd, HSTX, 4/95
c
c     Input: basis              rec  Legendre polynomial basis vectors.
c            coa_rec            rec  Coadd record.
c
c     Output: coa_rec            rec  Coadd record.
c
c     Modifications:
c
c----------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fil_normal
c
c     Input parameters.
c
      dictionary 'fex_basis'
      record /fex_basis/ basis
c
c     Input/output parameters.
c
      dictionary 'fil_sky'
      record /fil_sky/ coa_rec
c
c     Local variables.
c
      integer*4 ifg_pos,trunc_ifg_pos,vec,num_trunc_ifg_pts,num_basis_cols
      real*8 difg(492),trunc_leg_poly(492,5),tol,coeffs(5),resid(492),baseline

      fil_baseline_sub = %loc(fil_normal)
c
c     Fill intermediate arrays with the interferogram and the polynomials.
c
      do ifg_pos = 21,512
         trunc_ifg_pos = ifg_pos - 20
         difg(trunc_ifg_pos) = coa_rec.coad_data.ifg(ifg_pos)
         do vec = 1,5
            trunc_leg_poly(trunc_ifg_pos,vec) = basis.leg_poly(ifg_pos,vec)
         end do
      end do
c
c     Use the imsl routine "dlsqrr" to least-squares fit the baseline of the 
c     truncated interferogram with the fourth order Legendre polynomial.
c
      num_trunc_ifg_pts = 492
      tol = 1.0d-7

      call dlsqrr(num_trunc_ifg_pts,5,trunc_leg_poly,num_trunc_ifg_pts,difg,tol,
     .            coeffs,resid,num_basis_cols)
c
c     Calculate the baseline for all 512 points using the coefficients 
c     from "dlsqrr" and subtract it pointwise from the interferogram.
c
      do ifg_pos = 1,512
         baseline = 0.0
         do vec = 1,5
            baseline = baseline + coeffs(vec) * basis.leg_poly(ifg_pos,vec)
         end do
         coa_rec.coad_data.ifg(ifg_pos) = coa_rec.coad_data.ifg(ifg_pos) - 
     .                                    baseline
      end do
c
c     Record the baseline coefficients in the coadd record.
c
      do vec = 1,5
         coa_rec.coad_spec_data.bl_coeffs(vec) = coeffs(vec)
      end do

      return
      end
