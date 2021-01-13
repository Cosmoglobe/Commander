!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
module healpix_fft
  use healpix_types
  use extension, only: exit_with_status
  implicit none
  private

  ! module for FFT operations
  ! edited for compatibility with Absoft Pro, G95 and GFORTRAN
  !
  public real_fft2, complex_fft2, make_fft2_plan, destroy_fft2_plan, &
       & planck_fft2_plan, init_fft2_plan, complex_fft, real_fft

  ! routines of Healpix 2.0, compatible with FFTW 2.x

  logical, parameter, public :: fft2_forward=.false., fft2_backward=.true.

  type planck_fft2_plan
     logical      :: direction
     integer(i4b) :: length
  end type planck_fft2_plan

  interface complex_fft2
     module procedure d_c_complex_fft2, d_r_complex_fft2
  end interface

  interface real_fft2
     module procedure s_real_fft2, d_real_fft2
  end interface

  ! routines of Healpix 1.2

  interface complex_fft
     module procedure complex_fft_orig, complex_fft_alt
  end interface

  interface real_fft
     module procedure s_real_fft, d_real_fft
  end interface

contains

  subroutine fft_gpd (data,nn,backward,onlyreal)
    real(dp), intent(inout) :: data(*)
    integer(i4b), intent(in) :: nn(:)
    logical, intent(in) :: backward, onlyreal

!  real(dp) :: work(2*maxval(nn)) ! not valid with Absoft
   real(dp) :: work(2*nn(1)) ! assumes 1 dimension
  real(dp) :: difr, difi, oldsr, oldsi, sumr, sumi, sinth
  real(dp) :: twowr, t2r,t2i, t3r, t3i, t4r, t4i
  real(dp) :: theta
  real(dp) :: tempr, tempi
  real(dp) :: u1r, u1i, u2r, u2i, u3r, u3i, u4r, u4i
  real(dp) :: wr, wi, w2r, w2i, w3r, w3i, wstpr, wstpi
  integer(i4b) :: i, imin, i1, i2, i3, imax, i1max, i2max, idiv, idim

  integer(i4b) :: ifact(32), iff, ifp1, ifp2, ipar, i1rng, icase, irem, iquot
  integer(i4b) :: j, jmin, j1, j2, j3, j1min, j2min, j2stp
  integer(i4b) :: jmax, j1max, j2max, j3max, j1cnj, j1rng, j1rg2
  integer(i4b) :: k, kconj, krang, kmin, k1, k2, k3, k4, kdif,kstep
  integer(i4b) :: l, lmax
  integer(i4b) :: m, mmax
  integer(i4b) :: n, nhalf, np0, np1, np2, np1hf, np2hf, non2, non2t, ntwo
  integer(i4b) :: ntot, nprev
  integer(i4b) :: ndim

!=======================================================================
  ndim = size(nn)

  nprev = 0
  np0 = 0
  wr = 0d0
  wi = 0d0
  w2r = 0d0
  w2i = 0d0
  w3r = 0d0
  w3i = 0d0
  wstpr = 0d0
  wstpi = 0d0

  if (ndim<1) return
  ntot=2
  do idim=1,ndim
    if (nn(idim)<=0) return
    if (2*nn(idim) > size(work)) call exit_with_status(1,"FFT_GDP: work array too small")
    ntot=ntot*nn(idim)
  end do

!---
!--- MAIN LOOP FOR EACH DIMENSION
!---
  np1 = 2
  dimloop: do idim=1, ndim
    n = nn(idim)
    np2 = np1*n
    if (n==1) goto 900
!---
!--- FACTOR N
!---
    m = n
    ntwo = np1
    iff = 1
    idiv = 2
    do
      iquot = m/idiv
      irem = m-idiv*iquot
      if (iquot<idiv) then
        if (irem/=0) then
          ifact(iff) = m
        else
          ntwo = ntwo+ntwo
        endif
        goto 70
      endif
      if (irem/=0) exit
      ntwo = ntwo+ntwo
      m = iquot
    end do

    idiv = 3
    do
      iquot = m/idiv
      irem = m-idiv*iquot
      if (iquot<idiv) then
        ifact(iff) = m
        goto 70
      endif
      if (irem==0) then
        ifact (iff) = idiv
        iff = iff+1
        m = iquot
      else
        idiv = idiv+2
      endif
    end do

!---
!--- SEPARATE FOUR CASES--
!1. COMPLEX TRANSFORM OR REAL TRANSFORM FOR THE 4TH, 5TH,ETC.
!   DIMENSIONS.
!2. REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION.  METHOD--
!   TRANSFORM HALF THE DATA, SUPPLYING THE OTHER HALF BY CON-
!   JUGATE SYMMETRY.
!3. REAL TRANSFORM FOR THE 1ST DIMENSION, N ODD.  METHOD--
!   TRANSFORM HALF THE DATA AT EACH STAGE, SUPPLYING THE OTHER
!   HALF BY CONJUGATE SYMMETRY.
!4. REAL TRANSFORM FOR THE 1ST DIMENSION, N EVEN.  METHOD--
!   TRANSFORM A COMPLEX ARRAY OF LENGTH N/2 WHOSE REAL PARTS
!   ARE THE EVEN NUMBERED REAL VALUES AND WHOSE IMAGINARY PARTS
!   ARE THE ODD NUMBERED REAL VALUES.  SEPARATE AND SUPPLY
!   THE SECOND HALF BY CONJUGATE SYMMETRY.
!---
70  non2 = np1*(np2/ntwo)
    if ((idim>=4) .or. (.not. onlyreal)) then
      icase = 1
    elseif (idim>1) then
      icase = 2
    elseif (ntwo<=np1) then
      icase = 3
    else
      icase = 4
      ntwo = ntwo/2
      n = n/2
      np2 = np2/2
      ntot = ntot/2
      i = 3
      do j=2,ntot
        data(j) = data(i)
        i = i+2
      end do
    endif
    i1rng = np1
    if (icase==2) i1rng = np0*(1+nprev/2)
!---
!--- SHUFFLE ON THE FACTORS OF TWO IN N.  AS THE SHUFFLING
!--- CAN BE DONE BY SIMPLE INTERCHANGE, NO WORKING ARRAY IS NEEDED
!---
    if (ntwo<=np1) goto 600
    np2hf = np2/2
    j = 1
    do i2=1,np2,non2
      if (j<i2) then
        i1max = i2+non2-2
        do i1=i2,i1max,2
          do i3=i1,ntot,np2
            j3 = j+i3-i2
            tempr = data(i3)
            tempi = data(i3+1)
            data(i3) = data(j3)
            data(i3+1) = data(j3+1)
            data(j3) = tempr
            data(j3+1) = tempi
          end do
        end do
      endif
      m = np2hf
      do
        if (j<=m) exit
        j = j-m
        m = m/2
        if (m<non2) exit
      end do
      j = j+m
    end do
!---
!--- MAIN LOOP FOR FACTORS OF TWO.  PERFORM FOURIER TRANSFORMS OF
!--- LENGTH FOUR, WITH ONE OF LENGTH TWO IF NEEDED.  THE TWIDDLE FACTOR
!--- W=EXP(ISIGN*2*PI*SQRT(-1)*M/(4*MMAX)).  CHECK FOR W=ISIGN*SQRT(-1)
!--- AND REPEAT FOR W=ISIGN*SQRT(-1)*CONJUGATE(W).
!---
    non2t = non2+non2
    ipar = ntwo/np1
    do while (ipar>2)
      ipar = ipar/4
    end do
    if (ipar==2) then
      do i1=1,i1rng,2
        do j3=i1,non2,np1
          do k1=j3,ntot,non2t
            k2 = k1+non2
            tempr = data(k2)
            tempi = data(k2+1)
            data(k2) = data(k1)-tempr
            data(k2+1) = data(k1+1)-tempi
            data(k1) = data(k1)+tempr
            data(k1+1) = data(k1+1)+tempi
          end do
        end do
      end do
    endif
    mmax = non2

    do while (mmax<np2hf)
      lmax = max(non2t,mmax/2)
      if (mmax>non2) then
        theta = -twopi*real(non2,kind=dp)/real(4*mmax,kind=dp)
        if (backward) theta = -theta
        wr = cos(theta)
        wi = sin(theta)
        wstpr = -2d0*wi*wi
        wstpi = 2d0*wr*wi
      endif
      do l=non2,lmax,non2t
        m = l
        if (mmax<=non2) goto 420
410     w2r = wr*wr-wi*wi
        w2i = 2d0*wr*wi
        w3r = w2r*wr-w2i*wi
        w3i = w2r*wi+w2i*wr
420     do i1=1,i1rng,2
          do j3=i1,non2,np1
            kmin = j3+ipar*m
            if (mmax<=non2) kmin = j3
            kdif = ipar*mmax
            do
              kstep = 4*kdif
              do k1=kmin,ntot,kstep
                k2 = k1+kdif
                k3 = k2+kdif
                k4 = k3+kdif
                if (mmax<=non2) then
                  u1r = data(k1)+data(k2)
                  u1i = data(k1+1)+data(k2+1)
                  u2r = data(k3)+data(k4)
                  u2i = data(k3+1)+data(k4+1)
                  u3r = data(k1)-data(k2)
                  u3i = data(k1+1)-data(k2+1)
                  if (.not. backward) then
                    u4r = data(k3+1)-data(k4+1)
                    u4i = data(k4)-data(k3)
                  else
                    u4r = data(k4+1)-data(k3+1)
                    u4i = data(k3)-data(k4)
                  endif
                else
                  t2r = w2r*data(k2)-w2i*data(k2+1)
                  t2i = w2r*data(k2+1)+w2i*data(k2)
                  t3r = wr*data(k3)-wi*data(k3+1)
                  t3i = wr*data(k3+1)+wi*data(k3)
                  t4r = w3r*data(k4)-w3i*data(k4+1)
                  t4i = w3r*data(k4+1)+w3i*data(k4)
                  u1r = data(k1)+t2r
                  u1i = data(k1+1)+t2i
                  u2r = t3r+t4r
                  u2i = t3i+t4i
                  u3r = data(k1)-t2r
                  u3i = data(k1+1)-t2i
                  if (.not. backward) then
                    u4r = t3i-t4i
                    u4i = t4r-t3r
                  else
                    u4r = t4i-t3i
                    u4i = t3r-t4r
                  endif
                endif
                data(k1) = u1r+u2r
                data(k1+1) = u1i+u2i
                data(k2) = u3r+u4r
                data(k2+1) = u3i+u4i
                data(k3) = u1r-u2r
                data(k3+1) = u1i-u2i
                data(k4) = u3r-u4r
                data(k4+1) = u3i-u4i
              end do
              kmin = 4*(kmin-j3)+j3
              kdif = kstep
              if (kdif>=np2) exit
            end do
          end do
        end do
        m = mmax-m
        if (.not. backward) then
          tempr = wr
          wr = -wi
          wi = -tempr
        else
          tempr = wr
          wr = wi
          wi = tempr
        endif
        if (m>lmax) goto 410
        tempr = wr
        wr = wr*wstpr-wi*wstpi+wr
        wi = wi*wstpr+tempr*wstpi+wi
      end do
      ipar = 3-ipar
      mmax = mmax+mmax
    end do
!---
!--- MAIN LOOP FOR FACTORS NOT EQUAL TO TWO.  APPLY THE TWIDDLE FACTOR
!--- W=EXP(ISIGN*2*PI*SQRT(-1)*(J2-1)*(J1-J2)/(NP2*IFP1)), THEN
!--- PERFORM A FOURIER TRANSFORM OF LENGTH IFACT(IFf), MAKING USE OF
!--- CONJUGATE SYMMETRIES.
!---
600 if (ntwo>=np2) goto 700
    ifp1 = non2
    iff = 1
    np1hf = np1/2
610 ifp2 = ifp1/ifact(iff)
    j1rng = np2
    if (icase==3) then
      j1rng = (np2+ifp1)/2
      j2stp = np2/ifact(iff)
      j1rg2 = (j2stp+ifp2)/2
    endif
    j2min = 1+ifp2
    if (ifp1<np2) then
      do j2=j2min,ifp1,ifp2
        theta = -twopi*real(j2-1,kind=dp)/real(np2,kind=dp)
        if (backward) theta = -theta
        sinth = sin(theta/2d0)
        wstpr = -2d0*sinth*sinth
        wstpi = sin(theta)
        wr = wstpr+1d0
        wi = wstpi
        j1min = j2+ifp1
        do j1=j1min,j1rng,ifp1
          i1max = j1+i1rng-2
          do i1=j1,i1max,2
            do i3=i1,ntot,np2
              j3max = i3+ifp2-np1
              do j3=i3,j3max,np1
                tempr = data(j3)
                data(j3) = data(j3)*wr-data(j3+1)*wi
                data(j3+1) = tempr*wi+data(j3+1)*wr
              end do
            end do
          end do
          tempr = wr
          wr = wr*wstpr-wi*wstpi+wr
          wi = tempr*wstpi+wi*wstpr+wi
        end do
      end do
    endif
    theta = -twopi/real(ifact(iff),kind=dp)
    if (backward) theta = -theta
    sinth = sin(theta/2d0)
    wstpr = -2d0*sinth*sinth
    wstpi = sin(theta)
    kstep = 2*n/ifact(iff)
    krang = kstep*(ifact(iff)/2)+1
    do i1=1,i1rng,2
      do i3=i1,ntot,np2
        do kmin=1,krang,kstep
          j1max = i3+j1rng-ifp1
          do j1=i3,j1max,ifp1
            j3max = j1+ifp2-np1
            do j3=j1,j3max,np1
              j2max = j3+ifp1-ifp2
              k = kmin+(j3-j1+(j1-i3)/ifact(iff))/np1hf
              if (kmin<=1) then
                sumr = 0d0
                sumi = 0d0
                do j2=j3,j2max,ifp2
                  sumr = sumr+data(j2)
                  sumi = sumi+data(j2+1)
                end do
                work(k) = sumr
                work(k+1) = sumi
              else
                kconj = k+2*(n-kmin+1)
                j2 = j2max
                sumr = data(j2)
                sumi = data(j2+1)
                oldsr = 0d0
                oldsi = 0d0
                j2 = j2-ifp2
                do
                  tempr = sumr
                  tempi = sumi
                  sumr = twowr*sumr-oldsr+data(j2)
                  sumi = twowr*sumi-oldsi+data(j2+1)
                  oldsr = tempr
                  oldsi = tempi
                  j2 = j2-ifp2
                  if (j2<=j3) exit
                end do
                tempr = wr*sumr-oldsr+data(j2)
                tempi = wi*sumi
                work(k) = tempr-tempi
                work(kconj) = tempr+tempi
                tempr = wr*sumi-oldsi+data(j2+1)
                tempi = wi*sumr
                work(k+1) = tempr+tempi
                work(kconj+1) = tempr-tempi
              endif
            end do
          end do
          if (kmin<=1) then
            wr = wstpr+1d0
            wi = wstpi
          else
            tempr = wr
            wr = wr*wstpr-wi*wstpi+wr
            wi = tempr*wstpi+wi*wstpr+wi
          endif
          twowr=wr+wr
        end do
        if ((icase/=3) .or. (ifp1>=np2)) then
          k = 1
          i2max = i3+np2-np1
          do i2=i3,i2max,np1
            data(i2) = work(k)
            data(i2+1) = work(k+1)
            k = k+2
          end do
        else
!---
!--- COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N ODD, BY CON-
!--- JUGATE SYMMETRIES AT EACH STAGE.
!---
          j3max = i3+ifp2-np1
          do j3=i3,j3max,np1
            j2max = j3+np2-j2stp
            do j2=j3,j2max,j2stp
              j1max = j2+j1rg2-ifp2
              j1cnj = j3+j2max+j2stp-j2
              do j1=j2,j1max,ifp2
                k = 1+j1-i3
                data(j1) = work(k)
                data(j1+1) = work(k+1)
                if (j1>j2) then
                  data(j1cnj) = work(k)
                  data(j1cnj+1)= -work(k+1)
                endif
                j1cnj = j1cnj-ifp2
              end do
            end do
          end do
        endif
      end do
    end do
    iff = iff+1
    ifp1 = ifp2
    if (ifp1>np1) goto 610
!---
!--- COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N EVEN, BY CON-
!--- JUGATE SYMMETRIES.
!---

700 if ((icase==1) .or. (icase==3)) goto 900
    if (icase==2) goto 800
    nhalf = n
    n = n+n
    theta = -twopi/real(n,kind=dp)
    if (backward) theta = -theta
    sinth = sin (theta/2d0)
    wstpr = -2d0*sinth*sinth
    wstpi = sin(theta)
    wr = wstpr+1d0
    wi = wstpi
    imin = 3
    jmin = 2*nhalf-1
    do while (imin<jmin)
      j = jmin
      do i=imin,ntot,np2
        sumr = (data(i)+data(j))/2d0
        sumi = (data(i+1)+data(j+1))/2d0
        difr = (data(i)-data(j))/2d0
        difi = (data(i+1)-data(j+1))/2d0
        tempr = wr*sumi+wi*difr
        tempi = wi*sumi-wr*difr
        data(i) = sumr+tempr
        data(i+1) = difi+tempi
        data(j) = sumr-tempr
        data(j+1) = -difi+tempi
        j = j+np2
      end do
      imin = imin+2
      jmin = jmin-2
      tempr = wr
      wr = wr*wstpr-wi*wstpi+wr
      wi = tempr*wstpi+wi*wstpr+wi
    end do
    if ((imin<=jmin) .and. (.not. backward)) then
      do i=imin,ntot,np2
        data(i+1) = -data(i+1)
      end do
    endif
    np2 = np2+np2
    ntot = ntot+ntot
    j = ntot+1
    imax = ntot/2+1
    do
      imin = imax-2*nhalf
      i = imin

      do
        i = i+2
        j = j-2
        if (i>=imax) exit
        data(j) = data(i)
        data(j+1) = -data(i+1)
      end do

      data(j) = data(imin)-data(imin+1)
      data(j+1) = 0d0

      if (i>=j) exit

      do
        i = i-2
        j = j-2
        if (i<=imin) exit
        data(j) = data(i)
        data(j+1) = data(i+1)
      end do
      data(j) = data(imin)+data(imin+1)
      data(j+1) = 0d0
      imax = imin
    end do
    data(1) = data(1)+data(2)
    data(2) = 0d0
    goto 900
!---
!--- COMPLETE A REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION BY
!--- CONJUGATE SYMMETRIES.
!---
800 if (i1rng>=np1) goto 900
    do i3=1,ntot,np2
      i2max = i3+np2-np1
      do i2=i3,i2max,np1
        imin = i2+i1rng
        imax = i2+np1-2
        jmax = 2*i3+np1-imin
        if (i2>i3) jmax = jmax+np2
        if (idim>2) then
          j = jmax+np0
          do i=imin,imax,2
            data(i) = data(j)
            data(i+1) = -data(j+1)
            j = j-2
          end do
        endif
        j = jmax
        do i=imin,imax,np0
          data(i) = data(j)
          data(i+1) = -data(j+1)
          j = j-np0
        end do
      end do
    end do
!---
!--- END OF LOOP ON EACH DIMENSION
!---
900 np0 = np1
    np1 = np2
    nprev = n
  end do dimloop
end subroutine fft_gpd

!================================================
subroutine init_fft2_plan(plan)
  ! sets initial values of the plan structure
  type(planck_fft2_plan), intent(inout) :: plan

  plan%direction = fft2_forward
  plan%length    = -1
end subroutine init_fft2_plan


subroutine sanity_check (plan, len)
  type(planck_fft2_plan), intent(in) :: plan
  integer, intent(in) :: len

  if (len/=plan%length) &
    call exit_with_status(1,"planck_fft: invalid plan for this transform")
end subroutine sanity_check

subroutine make_fft2_plan (plan,length,direction)
  type(planck_fft2_plan), intent(out) :: plan
  integer, intent(in) :: length
  logical, intent(in) :: direction

  plan%length=length
  plan%direction=direction
end subroutine make_fft2_plan

subroutine destroy_fft2_plan (plan)
  type(planck_fft2_plan), intent(out) :: plan

  plan%direction=fft2_forward
  plan%length=-1
end subroutine destroy_fft2_plan

subroutine d_c_complex_fft2 (plan, data)
  ! replaced TRANSFER functions by loops for compatibility with GFORTRAN
  type(planck_fft2_plan), intent(in) :: plan
  complex(dp), intent(inout) :: data(:)

  real(dp) data2(2*size(data))
  integer :: i, lb

  lb = lbound(data,   dim = 1)
  call sanity_check (plan, size(data))
  do i=1,size(data)
    data2(2*i-1)=real (data(i+lb-1),kind=dp)
    data2(2*i)  =aimag(data(i+lb-1))
  end do
  call fft_gpd (data2,(/size(data)/),plan%direction,.false.)
  do i=1,size(data)
    data(i+lb-1) = cmplx(data2(2*i-1), data2(2*i), kind=dp)
  end do
end subroutine d_c_complex_fft2

subroutine d_r_complex_fft2 (plan, data)
  type(planck_fft2_plan), intent(in) :: plan
  real(dp), intent(inout) :: data(:)

  call sanity_check (plan, size(data)/2)
  call fft_gpd (data,(/size(data)/2/),plan%direction,.false.)
end subroutine d_r_complex_fft2

subroutine d_real_fft2 (plan, data)
  type(planck_fft2_plan), intent(in) :: plan
  real(dp), intent(inout) :: data(:)

  real(dp) data2(2*size(data))
  integer n, i

  call sanity_check (plan, size(data))
  n = size(data)

  if (plan%direction .eqv. fft2_forward) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),plan%direction,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),plan%direction,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine d_real_fft2

subroutine s_real_fft2 (plan, data)
  type(planck_fft2_plan), intent(in) :: plan
  real(sp), intent(inout) :: data(:)

  real(dp) data2(2*size(data))
  integer n, i

  call sanity_check (plan, size(data))
  n = size(data)

  if (plan%direction .eqv. fft2_forward) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),plan%direction,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),plan%direction,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine s_real_fft2
!======================================================================

subroutine complex_fft_orig (data, backward, onlyreal)
  complex(dp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward, onlyreal

  real(dp) data2(2*size(data))
  logical or, bw
  or = .false.
  if (present(onlyreal)) or=onlyreal
  bw = .false.
  if (present(backward)) bw=backward

  data2 = transfer (data, data2, size=size(data2))
  call fft_gpd (data2,(/size(data)/),bw,or)
  data = transfer (data2, data, size=size(data))
end subroutine

subroutine complex_fft_alt (data, backward, onlyreal)
  real(dp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward, onlyreal

  logical or, bw

  or = .false.
  if (present(onlyreal)) or=onlyreal
  bw = .false.
  if (present(backward)) bw=backward

  call fft_gpd (data,(/size(data)/2/),bw,or)
end subroutine

subroutine d_real_fft (data, backward)
  real(dp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward

  real(dp) data2(2*size(data))
  logical bw
  integer n, i

  n = size(data)
  bw = .false.
  if (present(backward)) bw=backward

  if (.not. bw) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),bw,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),bw,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine

subroutine s_real_fft (data, backward)
  real(sp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward

  real(dp) data2(2*size(data))
  logical bw
  integer n, i

  n = size(data)
  bw = .false.
  if (present(backward)) bw=backward

  if (.not. bw) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),bw,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),bw,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine
!==================================================================
end module healpix_fft
