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
module comm_tod_noise_psd_mod
  ! 
  !   Module for defining noise PSD models
  !
  !   Main Classes
  !   ------------
  !   comm_noise_psd
  !       Class for basic 1/f noise model
  !   comm_noise_psd_2oof
  !       Class for two-component 1/f noise model
  use comm_utils
  implicit none

  private
  public comm_noise_psd, comm_noise_psd_2oof, comm_noise_psd_oof_gauss, comm_noise_psd_oof_f

  integer(i4b), parameter :: SIGMA0 = 1
  integer(i4b), parameter :: FKNEE  = 2
  integer(i4b), parameter :: ALPHA  = 3
  integer(i4b), parameter :: FKNEE2 = 4
  integer(i4b), parameter :: ALPHA2 = 5
  integer(i4b), parameter :: G_AMP  = 4
  integer(i4b), parameter :: G_LOC  = 5
  integer(i4b), parameter :: G_SIG  = 6

  type :: comm_noise_psd
     ! 
     ! Class definition for basic 1/f noise PSD model
     !
     integer(i4b) :: npar                                            ! Number of free parameters
     real(sp),     allocatable, dimension(:,:)  :: nu_fit            ! Frequency range used to fit non-linear parameters
     real(sp),     pointer :: sigma0                                 ! Pointer to xi_n(1)
     real(sp),     allocatable, dimension(:)    :: xi_n              ! Active sampling parameters, xi_n(1) = sigma0
     real(sp),     allocatable, dimension(:,:)  :: P_uni             ! Uniform prior on xi_n (n_xi,lower/upper)
     real(sp),     allocatable, dimension(:,:)  :: P_active          ! Informative prior on xi_n (n_xi, mean/rms)
     logical(lgt), allocatable, dimension(:)    :: P_lognorm         ! true = lognorm prior; false = Gaussian prior
   contains
     procedure :: eval_full   => eval_noise_psd_full
     procedure :: eval_corr   => eval_noise_psd_corr
  end type comm_noise_psd

  interface comm_noise_psd
     procedure constructor_oof
  end interface comm_noise_psd


  type, extends(comm_noise_psd) :: comm_noise_psd_2oof
     ! 
     ! Class definition for 2-component 1/f noise PSD model
     !
   contains
     procedure :: eval_full   => eval_noise_psd_2oof_full
     procedure :: eval_corr   => eval_noise_psd_2oof_corr
  end type comm_noise_psd_2oof

  type, extends(comm_noise_psd) :: comm_noise_psd_oof_gauss
     ! 
     ! Class definition for 2-component 1/f + Gauss noise PSD model
     !
   contains
     procedure :: eval_full   => eval_noise_psd_oof_gauss_full
     procedure :: eval_corr   => eval_noise_psd_oof_gauss_corr
  end type comm_noise_psd_oof_gauss

  type, extends(comm_noise_psd) :: comm_noise_psd_oof_f
     ! 
     ! Class definition for 2-component 1/f + Gauss noise PSD model
     !
   contains
     procedure :: eval_full   => eval_noise_psd_oof_f_full
     procedure :: eval_corr   => eval_noise_psd_oof_f_corr
  end type comm_noise_psd_oof_f

  interface comm_noise_psd_2oof
     procedure constructor_2oof
  end interface comm_noise_psd_2oof

  interface comm_noise_psd_oof_gauss
     procedure constructor_oof_gauss
  end interface comm_noise_psd_oof_gauss

  interface comm_noise_psd_oof_f
     procedure constructor_oof_f
  end interface comm_noise_psd_oof_f

contains

  function constructor_oof(P_active_mean, P_active_rms, P_uni, nu_fit)
    ! 
    ! Constructor for basic 1/f noise PSD object, where
    !     
    !     P(nu) = sigma0^2 * (1 + (nu/fknee)^alpha)
    ! 
    ! Arguments
    ! --------- 
    ! xi_n_def: sp (array)
    !          3-element array containing default {sigma0, alpha, fknee}, where
    !          [sigma0] = du/volts/tod unit, [alpha] = 1, and [fknee] = Hz
    ! P_uni: sp (2D array)
    !          Array containing absolute upper and lower limits for each parameter (npar,upper/lower)
    ! P_active: sp (2D array)
    !          Array containing informative priors for each parameter (npar,mean/rms)
    ! nu_fit: sp (2-element array)
    !          Array with [nu_min,nu_max] in Hz, defining ranged used for fittig non-linear parameters
    ! 
    implicit none
    real(sp),              dimension(:),   intent(in)      :: P_active_mean
    real(sp),              dimension(:),   intent(in)      :: P_active_rms
    real(sp),              dimension(:,:), intent(in)      :: P_uni
    real(sp),              dimension(:,:), intent(in)      :: nu_fit
    class(comm_noise_psd), pointer                         :: constructor_oof

    allocate(constructor_oof)

    if (P_active_mean(FKNEE) <= 0.0)     write(*,*) 'comm_noise_psd error: Default fknee less than zero'
    if (P_uni(FKNEE,1) <= 0.0)           write(*,*) 'comm_noise_psd error: Lower fknee prior less than zero'
    if (P_uni(FKNEE,1) > P_uni(FKNEE,2)) write(*,*) 'comm_noise_psd error: Lower fknee prior higher than upper prior'
    if (P_uni(ALPHA,1) > P_uni(ALPHA,2)) write(*,*) 'comm_noise_psd error: Lower alpha prior higher than upper prior'

    constructor_oof%npar = 3
    allocate(constructor_oof%xi_n(constructor_oof%npar))
    allocate(constructor_oof%nu_fit(constructor_oof%npar, 2))
    allocate(constructor_oof%P_uni(constructor_oof%npar,2))
    allocate(constructor_oof%P_active(constructor_oof%npar,2))
    allocate(constructor_oof%P_lognorm(constructor_oof%npar))

    constructor_oof%xi_n          = P_active_mean
    constructor_oof%P_uni         = P_uni
    constructor_oof%P_active(:,1) = P_active_mean
    constructor_oof%P_active(:,2) = P_active_rms
    constructor_oof%nu_fit        = nu_fit
    constructor_oof%P_lognorm     = [.false., .true., .false.] ! [sigma0, fknee, alpha]

    constructor_oof%sigma0 => constructor_oof%xi_n(1)

  end function constructor_oof

  function eval_noise_psd_full(self, nu)
    ! 
    ! Evaluation routine for basic 1/f noise PSD object
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd),               intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_full

    eval_noise_psd_full = self%xi_n(SIGMA0)**2 * (1. + (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA))

  end function eval_noise_psd_full

  function eval_noise_psd_corr(self, nu)
    ! 
    ! Evaluation routine for basic 1/f noise PSD object; correlated noise only
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd),               intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_corr

    eval_noise_psd_corr = self%xi_n(SIGMA0)**2 * (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA)

  end function eval_noise_psd_corr

  function constructor_2oof(P_active_mean, P_active_rms, P_uni, nu_fit)
    ! 
    ! Constructor for two-component 1/f noise PSD object, where
    !     
    !     P(nu) = sigma0^2 * (1 + (nu/fknee)^alpha + (nu/fknee2)^alpha2)
    ! 
    ! Arguments
    ! --------- 
    ! xi_n_def: sp (array)
    !          5-element array containing default {sigma0, fknee, alpha, fknee2, alpha2}, where
    !          [sigma0] = du/volts/tod unit, [alpha] = 1, and [fknee] = Hz
    ! P_uni: sp (2D array)
    !          Array containing absolute upper and lower limits for each parameter (npar,upper/lower)
    ! P_active: sp (2D array)
    !          Array containing informative priors for each parameter (npar,mean/rms)
    ! nu_fit: sp (2-element array)
    !          Array with [nu_min,nu_max] in Hz, defining ranged used for fittig non-linear parameters
    ! 
    implicit none
    real(sp),                   dimension(:),   intent(in)      :: P_active_mean
    real(sp),                   dimension(:),   intent(in)      :: P_active_rms
    real(sp),                   dimension(:,:), intent(in)      :: P_uni
    real(sp),                   dimension(:,:), intent(in)      :: nu_fit
    class(comm_noise_psd_2oof), pointer                         :: constructor_2oof

    allocate(constructor_2oof)

    if (P_active_mean(FKNEE) <= 0.0)       write(*,*) 'comm_noise_psd error: fknee prior mean less than zero'
    if (P_active_mean(FKNEE2) <= 0.0)      write(*,*) 'comm_noise_psd error: fknee2 prior mean less than zero'
    if (P_uni(FKNEE,1) <= 0.0)             write(*,*) 'comm_noise_psd error: Lower fknee prior less than zero'
    if (P_uni(FKNEE2,1) <= 0.0)            write(*,*) 'comm_noise_psd error: Lower fknee2 prior less than zero'
    if (P_uni(FKNEE,1) > P_uni(FKNEE,2))   write(*,*) 'comm_noise_psd error: Lower fknee prior higher than upper prior'
    if (P_uni(ALPHA,1) > P_uni(ALPHA,2))   write(*,*) 'comm_noise_psd error: Lower alpha prior higher than upper prior'
    if (P_uni(FKNEE2,1) > P_uni(FKNEE2,2)) write(*,*) 'comm_noise_psd error: Lower fknee2 prior higher than upper prior'
    if (P_uni(ALPHA2,1) > P_uni(ALPHA2,2)) write(*,*) 'comm_noise_psd error: Lower alpha2 prior higher than upper prior'

    constructor_2oof%npar = 5
    allocate(constructor_2oof%xi_n(constructor_2oof%npar))
    allocate(constructor_2oof%nu_fit(constructor_2oof%npar, 2))
    allocate(constructor_2oof%P_uni(constructor_2oof%npar,2))
    allocate(constructor_2oof%P_active(constructor_2oof%npar,2))
    allocate(constructor_2oof%P_lognorm(constructor_2oof%npar))

    constructor_2oof%xi_n          = P_active_mean
    constructor_2oof%P_uni         = P_uni
    constructor_2oof%P_active(:,1) = P_active_mean
    constructor_2oof%P_active(:,2) = P_active_rms
    constructor_2oof%nu_fit        = nu_fit
    constructor_2oof%P_lognorm     = [.false., .true., .false., .true., .false.] !  [sigma0, fknee, alpha, fknee2, alpha2]

    constructor_2oof%sigma0 => constructor_2oof%xi_n(1)

  end function constructor_2oof
  
  function eval_noise_psd_2oof_full(self, nu)
    ! 
    ! Evaluation routine for 2-component 1/f noise PSD object
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_2oof),          intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_2oof_full

    eval_noise_psd_2oof_full = self%xi_n(SIGMA0)**2 * (1. + (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA) + (nu/self%xi_n(FKNEE2))**self%xi_n(ALPHA2))

  end function eval_noise_psd_2oof_full

  function eval_noise_psd_2oof_corr(self, nu)
    ! 
    ! Evaluation routine for 2-component 1/f noise PSD object; correlated noise only
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_2oof),          intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_2oof_corr

    eval_noise_psd_2oof_corr = self%xi_n(SIGMA0)**2 * ((nu/self%xi_n(FKNEE))**self%xi_n(ALPHA) + (nu/self%xi_n(FKNEE2))**self%xi_n(ALPHA2))

  end function eval_noise_psd_2oof_corr


  function constructor_oof_gauss(P_active_mean, P_active_rms, P_uni, nu_fit)
    ! 
    ! Constructor for two-component 1/f noise PSD object, where
    !     
    !     P(nu) = sigma0^2 * (1 + (nu/fknee)^alpha + (nu/fknee2)^alpha2)
    ! 
    ! Arguments
    ! --------- 
    ! xi_n_def: sp (array)
    !          5-element array containing default {sigma0, fknee, alpha, fknee2, alpha2}, where
    !          [sigma0] = du/volts/tod unit, [alpha] = 1, and [fknee] = Hz
    ! P_uni: sp (2D array)
    !          Array containing absolute upper and lower limits for each parameter (npar,upper/lower)
    ! P_active: sp (2D array)
    !          Array containing informative priors for each parameter (npar,mean/rms)
    ! nu_fit: sp (2-element array)
    !          Array with [nu_min,nu_max] in Hz, defining ranged used for fittig non-linear parameters
    ! 
    implicit none
    real(sp),                        dimension(:),   intent(in)      :: P_active_mean
    real(sp),                        dimension(:),   intent(in)      :: P_active_rms
    real(sp),                        dimension(:,:), intent(in)      :: P_uni
    real(sp),                        dimension(:,:), intent(in)      :: nu_fit
    class(comm_noise_psd_oof_gauss), pointer                         :: constructor_oof_gauss

    allocate(constructor_oof_gauss)

    if (P_active_mean(FKNEE) <= 0.0)       write(*,*) 'comm_noise_psd error: fknee prior mean less than zero'
    if (P_active_mean(G_SIG) <= 0.0)     write(*,*) 'comm_noise_psd error: g_sig prior mean less than zero'
    if (P_uni(FKNEE,1) <= 0.0)             write(*,*) 'comm_noise_psd error: Lower fknee prior less than zero'
    if (P_uni(FKNEE,1) > P_uni(FKNEE,2))   write(*,*) 'comm_noise_psd error: Lower fknee prior higher than upper prior'
    if (P_uni(ALPHA,1) > P_uni(ALPHA,2))   write(*,*) 'comm_noise_psd error: Lower alpha prior higher than upper prior'

    constructor_oof_gauss%npar = 6
    allocate(constructor_oof_gauss%xi_n(constructor_oof_gauss%npar))
    allocate(constructor_oof_gauss%nu_fit(constructor_oof_gauss%npar, 2))
    allocate(constructor_oof_gauss%P_uni(constructor_oof_gauss%npar,2))
    allocate(constructor_oof_gauss%P_active(constructor_oof_gauss%npar,2))
    allocate(constructor_oof_gauss%P_lognorm(constructor_oof_gauss%npar))

    constructor_oof_gauss%xi_n          = P_active_mean
    !write(*,*) size(constructor_oof_gauss%P_uni, 1), size(constructor_oof_gauss%P_uni, 2), size(P_uni, 1), size(P_uni,2)
    !write(*,*) P_uni
    constructor_oof_gauss%P_uni         = P_uni
    constructor_oof_gauss%P_active(:,1) = P_active_mean
    constructor_oof_gauss%P_active(:,2) = P_active_rms
    constructor_oof_gauss%nu_fit        = nu_fit
    constructor_oof_gauss%P_lognorm     = [.false., .true., .false., .false., .false., .false.] !  [sigma0, fknee, alpha, fknee2, alpha2]

    constructor_oof_gauss%sigma0 => constructor_oof_gauss%xi_n(1)

  end function constructor_oof_gauss
  
  function eval_noise_psd_oof_gauss_full(self, nu)
    ! 
    ! Evaluation routine for 2-component 1/f noise PSD object
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_oof_gauss),     intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_gauss_full
 
!!$    real(sp) :: S1, S2
!!$
!!$    S1 = self%xi_n(SIGMA0)**2 * (1. + (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA))
!!$    S2 = self%xi_n(SIGMA0)**2 * self%xi_n(G_AMP) / nu * exp(-0.5 * ((log10(nu) - log10(self%xi_n(G_LOC))/self%xi_n(G_SIG)))**2 ) 
!!$    eval_noise_psd_oof_gauss_full = S1 + S2

    eval_noise_psd_oof_gauss_full = self%xi_n(SIGMA0)**2 * (1. + (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA)) + self%xi_n(SIGMA0)**2 * self%xi_n(G_AMP) / nu * exp(-0.5 * ((log10(nu) - log10(self%xi_n(G_LOC))/self%xi_n(G_SIG)))**2 ) 

  end function eval_noise_psd_oof_gauss_full

  function eval_noise_psd_oof_gauss_corr(self, nu)
    ! 
    ! Evaluation routine for 2-component 1/f noise PSD object; correlated noise only
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_oof_gauss),          intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_gauss_corr

    real(sp) :: S1, S2

    S1 = self%xi_n(SIGMA0)**2 * (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA)
    S2 = self%xi_n(SIGMA0)**2 * self%xi_n(G_AMP) / nu * exp(-0.5 * (log10(nu) - log10(self%xi_n(G_LOC))/self%xi_n(G_SIG)) ** 2) 

    eval_noise_psd_oof_gauss_corr = S1 + S2

  end function eval_noise_psd_oof_gauss_corr

  function constructor_oof_f(P_active_mean, P_active_rms, P_uni, nu_fit)
    ! 
    ! Constructor for two-component 1/f noise PSD object, where
    !     
    !     P(nu) = sigma0^2 * (1 + (nu/fknee)^alpha + g*(nu/fknee))
    ! 
    ! Arguments
    ! --------- 
    ! xi_n_def: sp (array)
    !          5-element array containing default {sigma0, fknee, alpha, fknee2, alpha2}, where
    !          [sigma0] = du/volts/tod unit, [alpha] = 1, and [fknee] = Hz
    ! P_uni: sp (2D array)
    !          Array containing absolute upper and lower limits for each parameter (npar,upper/lower)
    ! P_active: sp (2D array)
    !          Array containing informative priors for each parameter (npar,mean/rms)
    ! nu_fit: sp (2-element array)
    !          Array with [nu_min,nu_max] in Hz, defining ranged used for fittig non-linear parameters
    ! 
    implicit none
    real(sp),                        dimension(:),   intent(in)      :: P_active_mean
    real(sp),                        dimension(:),   intent(in)      :: P_active_rms
    real(sp),                        dimension(:,:), intent(in)      :: P_uni
    real(sp),                        dimension(:,:), intent(in)      :: nu_fit
    class(comm_noise_psd_oof_f),     pointer                         :: constructor_oof_f

    allocate(constructor_oof_f)

    if (P_active_mean(FKNEE) <= 0.0)       write(*,*) 'comm_noise_psd error: fknee prior mean less than zero'
    if (P_uni(FKNEE,1) <= 0.0)             write(*,*) 'comm_noise_psd error: Lower fknee prior less than zero'
    if (P_uni(FKNEE,1) > P_uni(FKNEE,2))   write(*,*) 'comm_noise_psd error: Lower fknee prior higher than upper prior'
    if (P_uni(ALPHA,1) > P_uni(ALPHA,2))   write(*,*) 'comm_noise_psd error: Lower alpha prior higher than upper prior'

    constructor_oof_f%npar = 4
    allocate(constructor_oof_f%xi_n(constructor_oof_f%npar))
    allocate(constructor_oof_f%nu_fit(constructor_oof_f%npar, 2))
    allocate(constructor_oof_f%P_uni(constructor_oof_f%npar,2))
    allocate(constructor_oof_f%P_active(constructor_oof_f%npar,2))
    allocate(constructor_oof_f%P_lognorm(constructor_oof_f%npar))

    constructor_oof_f%xi_n          = P_active_mean
    !write(*,*) size(constructor_oof_f%P_uni, 1), size(constructor_oof_f%P_uni, 2), size(P_uni, 1), size(P_uni,2)
    !write(*,*) P_uni
    constructor_oof_f%P_uni         = P_uni
    constructor_oof_f%P_active(:,1) = P_active_mean
    constructor_oof_f%P_active(:,2) = P_active_rms
    constructor_oof_f%nu_fit        = nu_fit
    constructor_oof_f%P_lognorm     = [.false., .true., .false., .false.] !  [sigma0, fknee, alpha, gamma]

    constructor_oof_f%sigma0 => constructor_oof_f%xi_n(1)

  end function constructor_oof_f
  
  function eval_noise_psd_oof_f_full(self, nu)
    ! 
    ! Evaluation routine for 2-component 1/f  + gamma*f noise PSD object
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_oof_f),         intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_f_full
 
    eval_noise_psd_oof_f_full = self%xi_n(SIGMA0)**2 + self%eval_corr(nu)

  end function eval_noise_psd_oof_f_full

  function eval_noise_psd_oof_f_corr(self, nu)
    ! 
    ! Evaluation routine for 1/f noise + gamma*f PSD object; correlated noise only
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_oof_f),          intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_f_corr

    real(sp) :: S1, S2

    S1 = self%xi_n(SIGMA0)**2 * (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA)
    S2 = self%xi_n(SIGMA0)**2 * self%xi_n(G_AMP) * (nu/self%xi_n(FKNEE))

    eval_noise_psd_oof_f_corr = S1 + S2

  end function eval_noise_psd_oof_f_corr
  
end module comm_tod_noise_psd_mod
