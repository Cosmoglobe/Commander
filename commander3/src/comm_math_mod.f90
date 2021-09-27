!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_math_mod
  implicit none

contains

  function calc_corr_coeff(spec_arr, n_spec, n_lim)
    implicit none
    integer(i4b),               intent(in) :: n_spec, n_lim
    real(dp),     dimension(:), intent(in) :: spec_arr
    real(dp)                               :: calc_corr_coeff

    integer(i4b) :: i, j, maxcorr, ns
    real(dp)     :: x_mean, y_mean, sig_x, sig_y, covarr
    real(dp), dimension(:), allocatable :: xarr, yarr
    !
    !  Function to calculate the Pearson correlation coefficient between samples in a set of a sampled spectral parameter.
    !
    !  Arguments:
    !  ------------------
    !  spec_arr: double precision array, unknown length
    !     An array containing the sampled spectral parameters
    !  n_spec: integer(i4b)
    !     The last sample number/index
    !  n_lim: integer(i4b)
    !     The maximum index to consider when computing the correlation coefficient
    !  Returns:
    !  -----------------
    !  calc_corr_coeff: double precision real number
    !     The computed correlation coefficient of the sample set.
    !

    !default if no corr coeff is found, return a value that is non-physical.
    !correlation coefficient is -1 <= coeff <= 1, by definition
    calc_corr_coeff = -2.d0 

    ns=n_spec
    if (ns > n_lim) ns = n_lim

    if (ns > 2) then
       allocate(xarr(ns),yarr(ns))
       do i = 1,ns
          xarr(i)=i*1.d0
          yarr(i)=spec_arr(i)
       end do
       x_mean=sum(xarr)/ns
       y_mean=sum(yarr)/ns
       sig_x = sqrt(sum((xarr-x_mean)**2)/ns)
       sig_y = sqrt(sum((yarr-y_mean)**2)/ns)
       if (sig_x > 0.d0 .and. sig_y > 0.d0) then
          covarr = sum((xarr-x_mean)*(yarr-y_mean))/ns
          calc_corr_coeff = covarr/(sig_x*sig_y)
       end if
       deallocate(xarr,yarr)
    end if

  end function calc_corr_coeff


  function calc_corr_len(spec_arr,n_spec) 
    implicit none
    integer(i4b),               intent(in) :: n_spec
    real(dp),     dimension(:), intent(in) :: spec_arr
    integer(i4b)                           :: calc_corr_len

    integer(i4b) :: i, j, maxcorr, ns
    real(dp)     :: x_mean, y_mean, sig_x, sig_y
    real(dp), dimension(:), allocatable :: xarr, yarr,covarr
    !
    !  Function to calculate the sampling correlation length of a spectral parameter sampling.
    !
    !  Arguments:
    !  ------------------
    !  spec_arr: double precision array, unknown length
    !     An array containing the sampled spectralparameters
    !  n_spec: integer(i4b)
    !     the total length of spec_arr of which to compute the correlation length. must be 0 < n_spec <= len(spec_arr).
    !
    !  Returns:
    !  -----------------
    !  calc_corr_len: integer(i4b)
    !     The computed sample correlation length of the sampled spectral parameters
    !

    calc_corr_len = -1 !default if no corr length is found

    maxcorr = n_spec/2
    allocate(xarr(maxcorr),yarr(maxcorr),covarr(maxcorr))

    do i = 1,maxcorr
       ns=maxcorr-i
       xarr=spec_arr(1:ns)
       yarr=spec_arr(1+i:ns+i)
       x_mean=sum(xarr)/ns
       y_mean=sum(yarr)/ns
       sig_x = sqrt(sum((xarr-x_mean)**2)/ns)
       sig_y = sqrt(sum((yarr-y_mean)**2)/ns)
       covarr(i) = sum((xarr-x_mean)*(yarr-y_mean))/(ns*sig_x*sig_y)
       if (covarr(i) < 0.d0) then
          calc_corr_len = i
          exit
       end if
    end do

    deallocate(xarr,yarr,covarr)

  end function calc_corr_len

