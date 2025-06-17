!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3
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
module comm_dust_extinction_mod
  use comm_utils
  use comm_map_mod
  implicit none

  private
  public initialize_dust_extinction_mod, get_dust_attenuation_map, get_dust_attenuation_pos, active_dust_ext_model

  integer(i4b), parameter :: NPAR_EXT = 47
  
  integer(i4b), parameter :: UV_C1_A  = 1
  integer(i4b), parameter :: UV_C1_B  = 2
  integer(i4b), parameter :: UV_C2_A  = 3
  integer(i4b), parameter :: UV_C2_B  = 4
  integer(i4b), parameter :: UV_C3_A  = 5
  integer(i4b), parameter :: UV_C3_B  = 6
  integer(i4b), parameter :: UV_C4_A  = 7
  integer(i4b), parameter :: UV_C4_B  = 8
  integer(i4b), parameter :: UV_X0    = 9
  integer(i4b), parameter :: UV_GAMMA = 10

  integer(i4b), parameter :: OPT_E0_A  = 11
  integer(i4b), parameter :: OPT_E0_B  = 12
  integer(i4b), parameter :: OPT_E1_A  = 13
  integer(i4b), parameter :: OPT_E1_B  = 14
  integer(i4b), parameter :: OPT_E2_A  = 15
  integer(i4b), parameter :: OPT_E2_B  = 16
  integer(i4b), parameter :: OPT_E3_A  = 17
  integer(i4b), parameter :: OPT_E3_B  = 18
  integer(i4b), parameter :: OPT_E4_A  = 19
  integer(i4b), parameter :: OPT_E4_B  = 20
  integer(i4b), parameter :: OPT_F1_A  = 21
  integer(i4b), parameter :: OPT_F1_B  = 22
  integer(i4b), parameter :: OPT_F2_A  = 23
  integer(i4b), parameter :: OPT_F2_B  = 24
  integer(i4b), parameter :: OPT_F3_A  = 25
  integer(i4b), parameter :: OPT_F3_B  = 26
  integer(i4b), parameter :: OPT_X1    = 27
  integer(i4b), parameter :: OPT_X2    = 28
  integer(i4b), parameter :: OPT_X3    = 29
  integer(i4b), parameter :: OPT_GAMMA1 = 30
  integer(i4b), parameter :: OPT_GAMMA2 = 31
  integer(i4b), parameter :: OPT_GAMMA3 = 32

  integer(i4b), parameter :: IR_G1_A     = 33
  integer(i4b), parameter :: IR_G1_B     = 34
  integer(i4b), parameter :: IR_ALPHA1_A = 35
  integer(i4b), parameter :: IR_ALPHA1_B = 36
  integer(i4b), parameter :: IR_ALPHA2   = 37
  integer(i4b), parameter :: IR_LAMBDAB  = 38
  integer(i4b), parameter :: IR_DELTA    = 39
  integer(i4b), parameter :: IR_S1       = 40
  integer(i4b), parameter :: IR_LAMBDAO1 = 41
  integer(i4b), parameter :: IR_GAMMAO1  = 42
  integer(i4b), parameter :: IR_A1       = 43
  integer(i4b), parameter :: IR_S2       = 44
  integer(i4b), parameter :: IR_LAMBDAO2 = 45
  integer(i4b), parameter :: IR_GAMMAO2  = 46
  integer(i4b), parameter :: IR_A2       = 47
  
  type comm_dustext_model
     integer(i4b) :: nside = 512
     real(dp),              dimension(NPAR_EXT) :: p
     real(sp), allocatable, dimension(:,:)      :: EBV
  end type comm_dustext_model

  real(dp) :: nu_V, nu_B, R
  type(comm_dustext_model) :: ext_model
  
contains

  subroutine initialize_dust_extinction_mod(cpar)
    implicit none
    type(comm_params), intent(in)    :: cpar

    real(dp) :: aB, bB, aV, bV

    if (trim(cpar%EBVmap) ==  'none') return
    
    ! Initialize sky model
    allocate(ext_model%EBV(0:12*ext_model%nside**2-1,1))
    call read_map(trim(cpar%EBVmap), ext_model%EBV)
    !info    => comm_mapinfo(cpar%comm_chain, nside, lmax, 1, .false.)
    !EBV_map => comm_map(info, trim(cpar%EBVmap))

    ! Initialize parametric model
    ext_model%p(UV_C1_A)     =  0.81297d0
    ext_model%p(UV_C1_B)     = -2.97868d0
    ext_model%p(UV_C2_A)     =  0.2775d0
    ext_model%p(UV_C2_B)     =  1.89808d0
    ext_model%p(UV_C3_A)     =  1.06295d0
    ext_model%p(UV_C3_B)     =  3.10334d0
    ext_model%p(UV_C4_A)     =  0.11303d0
    ext_model%p(UV_C4_B)     =  0.65484d0
    ext_model%p(UV_X0)       =  4.60d0
    ext_model%p(UV_GAMMA)    =  0.99d0

    ext_model%p(OPT_E0_A)    = -0.35848d0
    ext_model%p(OPT_E0_B)    =  0.12354d0
    ext_model%p(OPT_E1_A)    =  0.7122d0
    ext_model%p(OPT_E1_B)    = -2.68335d0
    ext_model%p(OPT_E2_A)    =  0.08746d0
    ext_model%p(OPT_E2_B)    =  2.01901d0
    ext_model%p(OPT_E3_A)    = -0.05403d0
    ext_model%p(OPT_E3_B)    = -0.39299d0
    ext_model%p(OPT_E4_A)    =  0.00674d0
    ext_model%p(OPT_E4_B)    =  0.03355d0
    ext_model%p(OPT_F1_A)    =  0.03893d0
    ext_model%p(OPT_F1_B)    =  0.18453d0
    ext_model%p(OPT_F2_A)    =  0.02965d0
    ext_model%p(OPT_F2_B)    =  0.19728d0
    ext_model%p(OPT_F3_A)    =  0.01747d0
    ext_model%p(OPT_F3_B)    =  0.1713d0
    ext_model%p(OPT_X1)      =  2.288d0
    ext_model%p(OPT_X2)      =  2.054d0
    ext_model%p(OPT_X3)      =  1.587d0
    ext_model%p(OPT_GAMMA1)  =  0.243d0
    ext_model%p(OPT_GAMMA2)  =  0.179d0
    ext_model%p(OPT_GAMMA3)  =  0.243d0

    ! Typo in paper --  extra minus in exponent in alpha1?
    ext_model%p(IR_G1_A)     =  0.38526d0
    ext_model%p(IR_G1_B)     = -1.01251d0
    ext_model%p(IR_ALPHA1_A) =  1.68467d0 
    ext_model%p(IR_ALPHA1_B) =  1.06099d0 ! Sign change
    ext_model%p(IR_ALPHA2)   =  0.78791d0
    ext_model%p(IR_LAMBDAB)  =  4.30578d0
    ext_model%p(IR_DELTA)    =  4.78338d0
    ext_model%p(IR_S1)       =  0.06652d0
    ext_model%p(IR_LAMBDAO1) =  9.8434d0
    ext_model%p(IR_GAMMAO1)  =  2.21205d0
    ext_model%p(IR_A1)       = -0.24703d0
    ext_model%p(IR_S2)       =  0.0267d0
    ext_model%p(IR_LAMBDAO2) =  19.258294d0
    ext_model%p(IR_GAMMAO2)  =  17.0d0
    ext_model%p(IR_A2)       = -0.27d0

    nu_V = c/551d-9
    nu_B = c/445d-9
    call compute_dust_ext_linfit(nu_V, aV, bV)
    call compute_dust_ext_linfit(nu_B, aB, bB)
    !R = (1.d0-(bB-bV))/((aB-aV)-(bB-bV)/3.1d0)
    R = 3.1d0  !(1.d0-bB)/(aB-bB/3.1d0)
    !R = 6.0d0  !(1.d0-bB)/(aB-bB/3.1d0)
!!$    write(*,*) 'a, b, V = ', aV, bV
!!$    write(*,*) 'a, b, B = ', aB, bB
!!$    write(*,*) 'R = ', R
    
  end subroutine initialize_dust_extinction_mod

  subroutine get_dust_attenuation_map(nu, A_ext)
    implicit none
    real(dp),                 intent(in)    :: nu
    class(comm_map),          intent(inout) :: A_ext

    integer(i4b) :: i
    real(dp)     :: f, a, b
    real(sp), allocatable, dimension(:) :: EBV
    
    ! Get A_V = R*E(B-V) at correct resolution
    if (A_ext%info%nside /= ext_model%nside) then
       allocate(EBV(0:A_ext%info%npix-1))
       call udgrade_ring(ext_model%EBV(:,1), ext_model%nside, EBV, A_ext%info%nside)
       A_ext%map(:,1) = R*EBV(A_ext%info%pix)
    else
       A_ext%map(:,1) = R*ext_model%EBV(A_ext%info%pix,1)
    end if

    ! Convert to attenuation
    call compute_dust_ext_linfit(nu, a, b)
    f = a + b * (1d0/R - 1d0/3.1d0)
    A_ext%map = exp(-A_ext%map*f/2.5d0)

!!$    call A_ext%writeFITS('A.fits')
!!$    call mpi_finalize(i)
!!$    stop
  end subroutine get_dust_attenuation_map

  subroutine get_dust_attenuation_pos(vec, nu, A_ext)
    implicit none
    real(dp), intent(in)    :: vec(3)
    real(dp), intent(in)    :: nu
    real(dp), intent(out)   :: A_ext

    integer(i4b) :: i, pix
    real(dp) :: A_V, a, b

    if (.not. allocated(ext_model%EBV)) then
       A_ext = 1.d0
       return
    end if
    
    call vec2pix_ring(ext_model%nside, vec, pix)
    A_V = R * ext_model%EBV(pix,1)
    call compute_dust_ext_linfit(nu, a, b)
    A_ext = exp(-0.4d0 * A_V * (a + b * (1d0/R - 1d0/3.1d0)))    
  end subroutine get_dust_attenuation_pos

  
  subroutine compute_dust_ext_linfit(nu, a, b)
    implicit none
    real(dp), intent(in)    :: nu
    real(dp), intent(out)   :: a, b
    
    real(dp) :: lambda, lam0, lamb, a_uv, b_uv, a_opt, b_opt, a_ir, b_ir, x, D, W, Dm1, Dm2, z, g1, F, delta, gamma

    lambda = c/nu * 1d6 ! wavelength in micron
    x      = 1d0/lambda
    D      = x*x / ((x-ext_model%p(UV_X0))**2 + (x*ext_model%p(UV_GAMMA))**2)
    
    ! Compute basic linear coefficients for each range
    if (lambda >= 0.091d0 .and. lambda < 0.33d0) then
       ! Ultra-violet
       F    = 0.5392d0*(x-5.9d0)**2 + 0.05644d0*(x-5.9d0)**3
       a_uv = ext_model%p(UV_C1_A)   + ext_model%p(UV_C2_A)*x + &
            & ext_model%p(UV_C3_A)*D + ext_model%p(UV_C4_A)*F
       b_uv = ext_model%p(UV_C1_B)   + ext_model%p(UV_C2_B)*x + &
            & ext_model%p(UV_C3_B)*D + ext_model%p(UV_C4_B)*F
    end if
    if (lambda >= 0.3d0 .and. lambda < 1.1d0) then
       ! Optical
       a_opt = ext_model%p(OPT_E0_A)      + ext_model%p(OPT_E1_A)*x    + &
            &  ext_model%p(OPT_E2_A)*x**2 + ext_model%p(OPT_E3_A)*x**3 + &
            &  ext_model%p(OPT_E4_A)*x**4 + &
            &  ext_model%p(OPT_F1_A)*D*ext_model%p(OPT_GAMMA1) + &
            &  ext_model%p(OPT_F2_A)*D*ext_model%p(OPT_GAMMA2) + &
            &  ext_model%p(OPT_F3_A)*D*ext_model%p(OPT_GAMMA3) 
       b_opt = ext_model%p(OPT_E0_B)      + ext_model%p(OPT_E1_B)*x    + &
            &  ext_model%p(OPT_E2_B)*x**2 + ext_model%p(OPT_E3_B)*x**3 + &
            &  ext_model%p(OPT_E4_B)*x**4 + &
            &  ext_model%p(OPT_F1_B)*D*ext_model%p(OPT_GAMMA1)**2 + &
            &  ext_model%p(OPT_F2_B)*D*ext_model%p(OPT_GAMMA2)**2 + &
            &  ext_model%p(OPT_F3_B)*D*ext_model%p(OPT_GAMMA3)**2
    end if
    if (lambda >= 0.9d0 .and. lambda < 32d0) then
       ! Infrared
       lamb  = ext_model%p(IR_LAMBDAB)
       lam0  = ext_model%p(IR_LAMBDAO1)
       delta = ext_model%p(IR_DELTA)
       gamma = 2.d0*ext_model%p(IR_GAMMAO1)/(1.d0+exp(ext_model%p(IR_A1)*(lambda-lam0)))
       Dm1   = (gamma/lam0)**2 / ((lambda/lam0-lam0/lambda)**2 + (gamma/lam0)**2)
       lam0  = ext_model%p(IR_LAMBDAO2)
       gamma = 2.d0*ext_model%p(IR_GAMMAO2)/(1.d0+exp(ext_model%p(IR_A2)*(lambda-lam0)))
       Dm2   = (gamma/lam0)**2 / ((lambda/lam0-lam0/lambda)**2 + (gamma/lam0)**2)
       W     = get_W_smooth_step(lambda, lamb, delta)
       g1    = ext_model%p(IR_G1_A)
       a_ir  = g1*lambda**(-ext_model%p(IR_ALPHA1_A))*(1.d0-W) + &
            & g1*lamb**(ext_model%p(IR_ALPHA2)-ext_model%p(IR_ALPHA1_A)) * &
            & lambda**(-ext_model%p(IR_ALPHA2)) * W + &
            & ext_model%p(IR_S1)*Dm1 + ext_model%p(IR_S2)*Dm2
       b_ir  = ext_model%p(IR_G1_B)*lambda**(-ext_model%p(IR_ALPHA1_B)) 
    end if

    ! Generate smooth function
    if (lambda >= 0.091d0 .and. lambda < 0.30d0) then
       a = a_uv
       b = b_uv
    else if (lambda >= 0.30d0 .and. lambda < 0.33d0) then
       W = get_W_smooth_step(lambda, 0.315d0, 0.03d0)
       a = (1.d0-W)*a_uv + W*a_opt
       b = (1.d0-W)*b_uv + W*b_opt
    else if (lambda >= 0.33d0 .and. lambda < 0.90d0) then
       a = a_opt
       b = b_opt
    else if (lambda >= 0.90d0 .and. lambda < 1.1d0) then
       W = get_W_smooth_step(lambda, 1.d0, 0.2d0)
       a = (1.d0-W)*a_opt + W*a_ir
       b = (1.d0-W)*b_opt + W*b_ir
    else if (lambda >= 1.10d0 .and. lambda < 32d0) then
       a = a_ir
       b = b_ir
    else
       a = 0.d0
       b = 0.d0
    end if
  end subroutine compute_dust_ext_linfit

  function get_W_smooth_step(lambda, lambda_b, delta)
    implicit none
    real(dp), intent(in) :: lambda, lambda_b, delta
    real(dp)             :: get_W_smooth_step

    real(dp) :: z

    z = (lambda - lambda_b + 0.5d0*delta)/delta
    if (z >= 0.d0 .and. z <= 1.d0) then
       get_W_smooth_step = 3.d0*z**2 -2.d0*z**3
    else
       get_W_smooth_step = 0.d0
    end if
  end function get_W_smooth_step

  function active_dust_ext_model(nu)
    implicit none
    real(dp), intent(in) :: nu
    logical(lgt)         :: active_dust_ext_model

    real(dp) :: lambda

    lambda = c/nu * 1d6 ! wavelength in microns    
    if (allocated(ext_model%EBV) .and. lambda >= 0.091d0 .and. lambda < 32.d0) then
       active_dust_ext_model = .true.
    else
       active_dust_ext_model = .false.
    end if
    
  end function active_dust_ext_model

end module comm_dust_extinction_mod
