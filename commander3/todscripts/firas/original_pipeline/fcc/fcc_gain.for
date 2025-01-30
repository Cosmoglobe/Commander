      integer*4 function fcc_gain(bol_avg,cbias_avg,volt_avg,nyquist_hz,model,
     .                            s0,tau,tbol,gain)
c-------------------------------------------------------------------------------
c
c     Purpose: Calculate detector responsivity, time constant, derived
c              bolometer temperature, and gain function.
c
c     Authors: R. Eplee, GSC, 3/93
c              S. Alexander, HSTX, 7/93, SER 11189
c
c     Input: bol_avg     r*4  Average bolometer temperature in degrees K.
c            cbias_avg   r*4  Average commanded bias in counts.
c            volt_avg    r*4  Average readout voltage in volts.
c            nyquist_hz  r*4  Nyquist frequency in Hertz.
c            model       rec  Calibration model.
c
c     Output: s0          r*8  Detector responsivity.
c             tau         r*8  Detector time constant.
c             tbol        r*8  Derived bolometer temperature.
c             gain        c*16(257)  Gain function.
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return statuses.
c
      external fcc_normal
c
c     Input parameters.
c
      real*4 bol_avg,cbias_avg,volt_avg,nyquist_hz

      dictionary 'fex_mod'
      record /fex_mod/ model
c
c     Output parameters.
c
      real*8 s0,tau,tbol
      complex*16 gain(257)
c
c     Local variables.
c
      real*8 r0,t0,g1,beta,rho,c3,c1,jo,jg
      real*8 bias_avg,v,r,x,z,g,c,dt,w,rtrans,itrans
      integer*4 pos

      fcc_gain = %loc(fcc_normal)
c
c     Extract the bolometer parameters from the calibration model record:
c     r0 is the detector resistance at infinite temperature; t0 is the 
c     characteristic temperature for the detector resistance function; g1 is the
c     coefficient of the detector thermal conductance; beta is the index of the 
c     temperature dependence of the detector thermal conductance; rho is the
c     electric field dependence of the detector resistance; c3 is the 
c     coefficient of the cubic heat capacity term; c1 is the coefficient of 
c     the linear heat capacity term; jo is the jfet offset; and jg is the
c     jfet gain.
c
      r0   = model.bolparm(1)
      t0   = model.bolparm(2)
      g1   = model.bolparm(3)
      beta = model.bolparm(4)
      rho  = model.bolparm(5)
      c3   = model.bolparm(6)
      c1   = model.bolparm(7)
      jo   = model.bolparm(8)
      jg   = model.bolparm(9)
c
c     Convert the average commanded bias from counts to volts.
c
      bias_avg = cbias_avg/fac_count_to_volt
c
c     Calculate the bolometer state:  v is the total detector bias; r is the
c     total detector resistance; x is the non-ideal electric field term; and
c     z is the detector impedance.
c
      v    = (volt_avg - jo) / jg
      r    = fac_load_resist*v / (bias_avg - v)
      x    = v*rho
      z    = r/r0/x
      tbol = bol_avg
c
c     Iterate to find the actual detector temperature.
c
      do pos = 1, 8
         tbol = t0 / dlog(z*tbol*dsinh(x/tbol))**2
      enddo
c
c     Find the detector responsivity and time constant:  g is the total
c     thermal conductance; c is the total heat capacity; and dt is the detector
c     temperature derivative.
c
      x   = tbol/x * dtanh(x/tbol)
      g   = g1*tbol**beta
      c   = c3*tbol**3 + c1*tbol
      dt  = 1.0/x - 1.0 - 0.5*dsqrt(t0/tbol)
      z   = (g*tbol*r + dt*v**2) / (g*tbol*r/x - dt*v**2)
      s0  = fac_erg_to_watt*r*(z-x) / (v*(z*r/fac_load_resist + 1.0)*(x+1.0))
      tau = c/g * (z+1.0)*(r*x+fac_load_resist)/((z*r+fac_load_resist)*(x+1.0))
c
c     Calculate the gain function.
c
      do pos = 1,257
         w = 2.0 * fac_pi * (pos-1) * (nyquist_hz/256.0)
         if (abs(model.transfer(pos)) .lt. 1.0D-10) then
            gain(pos) = 0.0
         else
            rtrans = dreal(model.transfer(pos))
            itrans = dimag(model.transfer(pos))
            gain(pos) = dcmplx((rtrans+w*tau*itrans)/(s0*(rtrans**2+itrans**2)),
     .                         (w*tau*rtrans-itrans)/(s0*(rtrans**2+itrans**2)))
         end if
      end do

      return
      end
