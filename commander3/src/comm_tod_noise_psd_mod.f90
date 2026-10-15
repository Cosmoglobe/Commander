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
  public comm_noise_psd, comm_noise_psd_white, comm_noise_psd_oof, comm_noise_psd_2oof, comm_noise_psd_oof_gauss, comm_noise_psd_oof_quad, comm_noise_psd_spline

  integer(i4b), parameter :: SIGMA0 = 1
  integer(i4b), parameter :: FKNEE  = 2
  integer(i4b), parameter :: ALPHA  = 3
  integer(i4b), parameter :: FKNEE2 = 4
  integer(i4b), parameter :: ALPHA2 = 5
  integer(i4b), parameter :: G_AMP  = 4
  integer(i4b), parameter :: G_LOC  = 5
  integer(i4b), parameter :: G_SIG  = 6
  integer(i4b), parameter :: SLOPE  = 4
  integer(i4b), parameter :: QUADRATIC = 5

  type, abstract :: comm_noise_psd
     ! 
     ! Class definition for basic 1/f noise PSD model
     !
     integer(i4b) :: npar                                            ! Number of free parameters
     real(sp),     allocatable, dimension(:,:)  :: nu_fit            ! Frequency range used to fit non-linear parameters

     real(sp),     pointer :: sigma0                                 ! Pointer to xi_n(1)
     real(sp)              :: sigma0_preproc                                 ! Sigma0 for preprocessing
     real(sp),     allocatable, dimension(:)    :: xi_n              ! Active sampling parameters, xi_n(1) = sigma0
     real(sp),     allocatable, dimension(:,:)  :: P_uni             ! Uniform prior on xi_n (n_xi,lower/upper)
     real(sp),     allocatable, dimension(:,:)  :: P_active          ! Informative prior on xi_n (n_xi, mean/rms)
     logical(lgt), allocatable, dimension(:)    :: P_lognorm         ! true = lognorm prior; false = Gaussian prior
     logical(lgt)                               :: apply_filter      ! If we should apply the spline filter to the noise output
     type(spline_type)                          :: modulation_filter ! The filter multiplied to the output noise
   contains
     procedure(eval_noise_psd_full), deferred :: eval_full
     procedure(eval_noise_psd_corr), deferred :: eval_corr
     procedure :: init_common
  end type comm_noise_psd

  abstract interface
     function eval_noise_psd_full(self, nu)
       import comm_noise_psd, sp
       ! 
       ! Evaluation routine for general noise PSD object; both correlated and uncorrelated components
       ! 
       ! Arguments
       ! ---------
       ! self:    derived type (comm_noise_psd)
       !          Basic noise PSD object
       ! nu:      sp (scalar)
       !          Frequency (in Hz) at which to evaluate PSD
       ! 
       implicit none
       class(comm_noise_psd),         intent(in)      :: self
       real(sp),                      intent(in)      :: nu
       real(sp)                                       :: eval_noise_psd_full
     end function eval_noise_psd_full
  end interface

  abstract interface
     function eval_noise_psd_corr(self, nu)
       import comm_noise_psd, sp
       ! 
       ! Evaluation routine for general noise PSD object; correlated noise only
       ! 
       ! Arguments
       ! ---------
       ! self:    derived type (comm_noise_psd)
       !          Basic noise PSD object
       ! nu:      sp (scalar)
       !          Frequency (in Hz) at which to evaluate PSD
       ! 
       implicit none
       class(comm_noise_psd),         intent(in)      :: self
       real(sp),                      intent(in)      :: nu
       real(sp)                                       :: eval_noise_psd_corr
     end function eval_noise_psd_corr
  end interface

  ! #######################################
  !     Specific noise class definition
  ! #######################################
  
  type, extends(comm_noise_psd) :: comm_noise_psd_white
     ! 
     ! Class definition for white noise PSD model
     !
   contains
     procedure :: eval_full   => eval_noise_psd_white_full
     procedure :: eval_corr   => eval_noise_psd_white_corr
  end type comm_noise_psd_white
  
  type, extends(comm_noise_psd) :: comm_noise_psd_oof
     ! 
     ! Class definition for basic 1/f noise PSD model
     !
   contains
     procedure :: eval_full   => eval_noise_psd_oof_full
     procedure :: eval_corr   => eval_noise_psd_oof_corr
  end type comm_noise_psd_oof
  
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

  type, extends(comm_noise_psd) :: comm_noise_psd_oof_quad
     ! 
     ! Class definition for 2-component 1/f + linear noise PSD model
     !
   contains
     procedure :: eval_full   => eval_noise_psd_oof_quad_full
     procedure :: eval_corr   => eval_noise_psd_oof_quad_corr
  end type comm_noise_psd_oof_quad

  type, extends(comm_noise_psd) :: comm_noise_psd_spline
     !
     ! Class definition for splined noise PSD model
     !
     real(sp) :: pcorr_offset = 0.d0
     type(spline_type) :: spline_profile    ! Spline profile for 'spline' type PSD
     integer(i4b) :: n_common_points ! Number of common points
     real(sp), allocatable, dimension(:,:) :: common_points ! coordinates of common points
   contains
     procedure :: eval_full   => eval_noise_psd_spline_full
     procedure :: eval_corr   => eval_noise_psd_spline_corr
     procedure :: update_spline
  end type comm_noise_psd_spline

  interface comm_noise_psd_white
     procedure constructor_white
  end interface comm_noise_psd_white

  interface comm_noise_psd_oof
     procedure constructor_oof
  end interface comm_noise_psd_oof

  interface comm_noise_psd_2oof
     procedure constructor_2oof
  end interface comm_noise_psd_2oof

  interface comm_noise_psd_oof_gauss
     procedure constructor_oof_gauss
  end interface comm_noise_psd_oof_gauss

  interface comm_noise_psd_oof_quad
     procedure constructor_oof_quad
  end interface comm_noise_psd_oof_quad
  
  interface comm_noise_psd_spline
     procedure constructor_spline
  end interface comm_noise_psd_spline

contains

  subroutine init_common(self, P_active_mean, P_active_rms, P_uni, nu_fit, filter)
    ! Contains common initialization for all classes of noise model
    !
    ! Arguments:
    !
    ! self : comm_tod_noise_psd_mod
    !        The noise model object to initialize
    implicit none
    class(comm_noise_psd), target,         intent(inout)   :: self
    real(sp),              dimension(:),   intent(in)      :: P_active_mean
    real(sp),              dimension(:),   intent(in)      :: P_active_rms
    real(sp),              dimension(:,:), intent(in)      :: P_uni
    real(sp),              dimension(:,:),  intent(in)      :: nu_fit
    real(dp),     optional, dimension(:,:), intent(in)   :: filter

    allocate(self%xi_n(self%npar))
    allocate(self%P_uni(self%npar,2))
    allocate(self%P_active(self%npar,2))
    allocate(self%P_lognorm(self%npar))

    self%xi_n          = P_active_mean
    self%P_uni         = P_uni
    self%P_active(:,1) = P_active_mean
    self%P_active(:,2) = P_active_rms
    self%nu_fit        = nu_fit

    self%sigma0 => self%xi_n(1)

    if(.not. present(filter)) then
      self%apply_filter = .false.
    else
      self%apply_filter = .true.
    end if

    if(self%apply_filter) then
      call spline(self%modulation_filter, filter(1,:), filter(2,:))
    end if
  end subroutine
  
  function constructor_white(P_active_mean, P_active_rms, P_uni, nu_fit, filter)
    ! 
    ! Constructor for basic white noise PSD object, where
    !     
    !     P(nu) = sigma0^2 
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
    ! filter : dp of size (2, :), optional
    !          The noise filter to apply, in the form (x(:), y(:)). Regions 
    !          outside x(1) -> x(n) will not be filtered. Omit this if you
    !          don't want filtering

    implicit none
    real(sp),               dimension(:),   intent(in)      :: P_active_mean
    real(sp),               dimension(:),   intent(in)      :: P_active_rms
    real(sp),               dimension(:,:), intent(in)      :: P_uni
    real(sp),               dimension(:,:), intent(in)      :: nu_fit
    real(dp),     optional, dimension(:,:), intent(in)      :: filter
    class(comm_noise_psd_white), pointer                    :: constructor_white

    allocate(constructor_white)

    constructor_white%npar = 1

    call constructor_white%init_common(P_active_mean, P_active_rms, P_uni, nu_fit, filter)

    constructor_white%P_lognorm     = [.false.] ! [sigma0]

  end function constructor_white

  function eval_noise_psd_white_full(self, nu)
    ! 
    ! Evaluation routine for basic white noise PSD object
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_white),         intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_white_full

    eval_noise_psd_white_full = self%xi_n(SIGMA0)**2 

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_white_full = eval_noise_psd_white_full * splint(self%modulation_filter, dble(nu))
      end if
    end if

  end function eval_noise_psd_white_full

  function eval_noise_psd_white_corr(self, nu)
    ! 
    ! Evaluation routine for basic white noise PSD object; correlated noise only, so this returns zero
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    ! 
    implicit none
    class(comm_noise_psd_white),         intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_white_corr

    eval_noise_psd_white_corr = 0.

  end function eval_noise_psd_white_corr
  
  function constructor_oof(P_active_mean, P_active_rms, P_uni, nu_fit, filter)
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
    ! filter : dp of size (2, :), optional
    !          The noise filter to apply, in the form (x(:), y(:)). Regions 
    !          outside x(1) -> x(n) will not be filtered. Omit this if you
    !          don't want filtering

    implicit none
    real(sp),               dimension(:),   intent(in)      :: P_active_mean
    real(sp),               dimension(:),   intent(in)      :: P_active_rms
    real(sp),               dimension(:,:), intent(in)      :: P_uni
    real(sp),               dimension(:,:), intent(in)      :: nu_fit
    real(dp),     optional, dimension(:,:), intent(in)      :: filter
    class(comm_noise_psd_oof), pointer                      :: constructor_oof

    allocate(constructor_oof)

    if (P_active_mean(FKNEE) <= 0.0)     write(*,*) 'comm_noise_psd error: Default fknee less than zero'
    if (P_uni(FKNEE,1) <= 0.0)           write(*,*) 'comm_noise_psd error: Lower fknee prior less than zero'
    if (P_uni(FKNEE,1) > P_uni(FKNEE,2)) write(*,*) 'comm_noise_psd error: Lower fknee prior higher than upper prior'
    if (P_uni(ALPHA,1) > P_uni(ALPHA,2)) write(*,*) 'comm_noise_psd error: Lower alpha prior higher than upper prior'

    constructor_oof%npar = 3

    call constructor_oof%init_common(P_active_mean, P_active_rms, P_uni, nu_fit, filter)

    constructor_oof%P_lognorm     = [.false., .true., .false.] ! [sigma0, fknee, alpha]

  end function constructor_oof  

  function eval_noise_psd_oof_full(self, nu)
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
    class(comm_noise_psd_oof),           intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_full

    eval_noise_psd_oof_full = self%xi_n(SIGMA0)**2 * (1. + (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA))

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_oof_full = eval_noise_psd_oof_full * splint(self%modulation_filter, dble(nu))
      end if
    end if

  end function eval_noise_psd_oof_full

  function eval_noise_psd_oof_corr(self, nu)
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
    class(comm_noise_psd_oof),           intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_corr

    eval_noise_psd_oof_corr = self%xi_n(SIGMA0)**2 * (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA)

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_oof_corr = eval_noise_psd_oof_corr * splint(self%modulation_filter, dble(nu))
      end if
    end if

  end function eval_noise_psd_oof_corr

  function constructor_2oof(P_active_mean, P_active_rms, P_uni, nu_fit, filter)
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
    ! filter : dp of size (2, :), optional
    !          The noise filter to apply, in the form (x(:), y(:)). Regions 
    !          outside x(1) -> x(n) will not be filtered. Omit this if you
    !          don't want filtering
 
    implicit none
    real(sp),                   dimension(:),   intent(in)      :: P_active_mean
    real(sp),                   dimension(:),   intent(in)      :: P_active_rms
    real(sp),                   dimension(:,:), intent(in)      :: P_uni
    real(sp),                   dimension(:,:),   intent(in)      :: nu_fit
    real(dp),     optional,     dimension(:,:), intent(in)      :: filter
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

    call constructor_2oof%init_common(P_active_mean, P_active_rms, P_uni, nu_fit, filter)

    constructor_2oof%P_lognorm     = [.false., .true., .false., .true., .false.] !  [sigma0, fknee, alpha, fknee2, alpha2]

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

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_2oof_full = eval_noise_psd_2oof_full * splint(self%modulation_filter, dble(nu))
      end if
    end if

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

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_2oof_corr = eval_noise_psd_2oof_corr * splint(self%modulation_filter, dble(nu))
      end if
    end if


  end function eval_noise_psd_2oof_corr


  function constructor_oof_gauss(P_active_mean, P_active_rms, P_uni, nu_fit, filter)
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
    ! filter : dp of size (2, :), optional
    !          The noise filter to apply, in the form (x(:), y(:)). Regions 
    !          outside x(1) -> x(n) will not be filtered. Omit this if you
    !          don't want filtering

    implicit none
    real(sp),                        dimension(:),   intent(in)      :: P_active_mean
    real(sp),                        dimension(:),   intent(in)      :: P_active_rms
    real(sp),                        dimension(:,:), intent(in)      :: P_uni
    real(sp),                        dimension(:,:),   intent(in)      :: nu_fit
    real(dp),     optional,          dimension(:,:), intent(in)      :: filter
    class(comm_noise_psd_oof_gauss), pointer                         :: constructor_oof_gauss

    allocate(constructor_oof_gauss)

    if (P_active_mean(FKNEE) <= 0.0)       write(*,*) 'comm_noise_psd error: fknee prior mean less than zero'
    if (P_active_mean(G_SIG) <= 0.0)     write(*,*) 'comm_noise_psd error: g_sig prior mean less than zero'
    if (P_uni(FKNEE,1) <= 0.0)             write(*,*) 'comm_noise_psd error: Lower fknee prior less than zero'
    if (P_uni(FKNEE,1) > P_uni(FKNEE,2))   write(*,*) 'comm_noise_psd error: Lower fknee prior higher than upper prior'
    if (P_uni(ALPHA,1) > P_uni(ALPHA,2))   write(*,*) 'comm_noise_psd error: Lower alpha prior higher than upper prior'

    constructor_oof_gauss%npar = 6

    call constructor_oof_gauss%init_common(P_active_mean, P_active_rms, P_uni, nu_fit, filter)

    !write(*,*) size(constructor_oof_gauss%P_uni, 1), size(constructor_oof_gauss%P_uni, 2), size(P_uni, 1), size(P_uni,2)
    !write(*,*) P_uni
    constructor_oof_gauss%P_lognorm     = [.false., .true., .false., .false., .false., .false.] !  [sigma0, fknee, alpha, fknee2, alpha2]

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

    eval_noise_psd_oof_gauss_full = self%xi_n(SIGMA0)**2 * (1. + (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA)) + self%xi_n(SIGMA0)**2 *self%xi_n(G_AMP) / nu * exp(-0.5 * ((log10(nu) - log10(self%xi_n(G_LOC)))/self%xi_n(G_SIG))**2 ) 

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_oof_gauss_full = eval_noise_psd_oof_gauss_full * splint(self%modulation_filter, dble(nu))
      end if
    end if


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
    S2 = self%xi_n(SIGMA0)**2 * self%xi_n(G_AMP) / nu * exp(-0.5 * ((log10(nu) - log10(self%xi_n(G_LOC)))/self%xi_n(G_SIG)) ** 2) 

    eval_noise_psd_oof_gauss_corr = S1 + S2

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_oof_gauss_corr = eval_noise_psd_oof_gauss_corr * splint(self%modulation_filter, dble(nu))
      end if
    end if

  end function eval_noise_psd_oof_gauss_corr

  function constructor_oof_quad(P_active_mean, P_active_rms, P_uni, nu_fit, filter)
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
    ! filter : dp of size (2, :), optional
    !          The noise filter to apply, in the form (x(:), y(:)). Regions 
    !          outside x(1) -> x(n) will not be filtered. Omit this if you
    !          don't want filtering
 
    implicit none
    real(sp),                        dimension(:),   intent(in)      :: P_active_mean
    real(sp),                        dimension(:),   intent(in)      :: P_active_rms
    real(sp),                        dimension(:,:), intent(in)      :: P_uni
    real(sp),                        dimension(:,:),   intent(in)      :: nu_fit
    real(dp),     optional,          dimension(:,:), intent(in)      :: filter
    class(comm_noise_psd_oof_quad),  pointer                         :: constructor_oof_quad

    allocate(constructor_oof_quad)

    if (P_active_mean(FKNEE) <= 0.0)       write(*,*) 'comm_noise_psd error: fknee prior mean less than zero'
    if (P_uni(FKNEE,1) <= 0.0)             write(*,*) 'comm_noise_psd error: Lower fknee prior less than zero'
    if (P_uni(FKNEE,1) > P_uni(FKNEE,2))   write(*,*) 'comm_noise_psd error: Lower fknee prior higher than upper prior'
    if (P_uni(ALPHA,1) > P_uni(ALPHA,2))   write(*,*) 'comm_noise_psd error: Lower alpha prior higher than upper prior'

    constructor_oof_quad%npar = 5

    call constructor_oof_quad%init_common(P_active_mean, P_active_rms, P_uni, nu_fit, filter)

    !write(*,*) size(constructor_oof_f%P_uni, 1), size(constructor_oof_f%P_uni, 2), size(P_uni, 1), size(P_uni,2)
    !write(*,*) P_uni
    constructor_oof_quad%P_lognorm  = [.false., .true., .false., .false., .false.] !  [sigma0, fknee, alpha, slope, interecept]


  end function constructor_oof_quad
  
  function eval_noise_psd_oof_quad_full(self, nu)
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
    class(comm_noise_psd_oof_quad),         intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_quad_full
 
    eval_noise_psd_oof_quad_full = self%xi_n(SIGMA0)**2 + self%eval_corr(nu)

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_oof_quad_full = eval_noise_psd_oof_quad_full * splint(self%modulation_filter, dble(nu))
      end if
    end if

  end function eval_noise_psd_oof_quad_full

  function eval_noise_psd_oof_quad_corr(self, nu)
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
    class(comm_noise_psd_oof_quad),         intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_oof_quad_corr

    real(sp) :: S1, S2, S3


    S1 = self%xi_n(SIGMA0)**2 * (nu/self%xi_n(FKNEE))**self%xi_n(ALPHA)
    S2 = self%xi_n(SLOPE) * nu + self%xi_n(QUADRATIC) * nu**2
    eval_noise_psd_oof_quad_corr = S1 + S2


    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_oof_quad_corr = eval_noise_psd_oof_quad_corr * splint(self%modulation_filter, dble(nu))
      end if
    end if


  end function eval_noise_psd_oof_quad_corr

  function constructor_spline(P_active_mean, P_active_rms, P_uni, nu_fit, filter)
    !
    ! Constructor for splined noise PSD object
    !
    ! Arguments
    ! ---------
    ! xi_n_def: sp (array)
    !          n-elements array containing white-noise level and the list of spline nodes
    !          {sigma0, n1,n2,n3,...}, where [sigma0] = du/volts/tod unit, [nodes] = Hz
    ! P_uni: sp (2D array)
    !          Array containing absolute upper and lower limits for each parameter (npar,upper/lower)
    ! P_active: sp (2D array)
    !          Array containing informative priors for each parameter (npar,mean/rms)
    ! nu_fit: sp (2-element array)
    !          Array with [nu_min,nu_max] in Hz, defining ranged used for fittig non-linear parameters
    !
    ! filter : dp of size (2, :), optional
    !          The noise filter to apply, in the form (x(:), y(:)). Regions
    !          outside x(1) -> x(n) will not be filtered. Omit this if you
    !          don't want filtering

    implicit none
    real(sp),               dimension(:),   intent(in)      :: P_active_mean
    real(sp),               dimension(:),   intent(in)      :: P_active_rms
    real(sp),               dimension(:,:), intent(in)      :: P_uni
    real(sp),               dimension(:,:), intent(in)      :: nu_fit
    real(dp),     optional, dimension(:,:), intent(in)      :: filter

    integer(i4b),                 parameter :: n_points = 19
    integer(i4b) :: i
    class(comm_noise_psd_spline), pointer                   :: constructor_spline

    allocate(constructor_spline)

    constructor_spline%npar = size(P_active_mean)
    call constructor_spline%init_common(P_active_mean, P_active_rms, P_uni, nu_fit, filter)
    constructor_spline%P_lognorm = .false.

    constructor_spline%n_common_points = n_points
    allocate(constructor_spline%common_points(n_points,2))
    constructor_spline%common_points = 0.d0
    constructor_spline%common_points(1:5,1) = (/1.d-5, 1.d-4, 1.d-3, 1.d-2, 1.d-1/) ! log points
    constructor_spline%common_points(6:,1)  = (/1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0/) ! lin points

  end function constructor_spline


  subroutine update_spline(self,x,y,handle)
    !
    ! Routine to update spline profile
    !
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! y:       sp (array)
    !          PSD values
    ! x:       sp (array)
    !          PSD Frequencies (Hz)
    implicit none
    class(comm_noise_psd_spline), target, intent(inout) :: self
    real(sp),           dimension(:),     intent(in)    :: x, y
    type(planck_rng),           optional, intent(inout) :: handle ! pass it to Gibb sample the nodes

    integer(i4b) :: i, j, n, nbin, n_dec
    real(sp) :: thr, dx_lin, wn_level, slope_lowf, f_nu
    real(dp) :: nu
    real(sp), allocatable, dimension(:)   :: sigmas_full_range
    real(sp), allocatable, dimension(:,:) :: bin_spec_full_range
    real(sp),     allocatable, dimension(:) :: xbin, pbin
    integer(i4b), allocatable, dimension(:) :: nbins

    n = size(x)
    if (size(y) /= n) stop 'update_spline: x and y have different sizes'

    n_dec  = 4 ! log bins per decade
    thr    = 1.d0 ! Hz  ! Frequency dividing logarithmic and linear regimes
    dx_lin = 8.d0 ! Hz
    call bin_psd_spline(x, y, thr, n_dec, dx_lin, 20, xbin, pbin, nbins) 

    nbin = size(xbin)
    if (any(pbin /= pbin)) stop 'update_spline: nan in binned PSD'
    if (nbin < 2) stop 'update_spline: less than 2 bins'


    ! Gibbs ssample the bins from the conditional posterior
    if (present(handle)) then
       do i = 1, nbin
          pbin(i) = draw_bin_power(pbin(i), nbins(i), handle)
       end do
    end if


    wn_level = minval(pbin, mask = xbin >= thr) ! setting wn level on lowest point after threshold

    ! Check if lower frequencies have less power. In that case inject regularization noise
    do i = 1, nbin
       if (xbin(i) < thr .and. pbin(i) < wn_level) pbin(i) = pbin(i) + wn_level
    end do

    allocate(bin_spec_full_range(nbin+2,2), sigmas_full_range(nbin+2)) ! extend spline over first and last bins

    ! Extend to low frequencies using the slope between first two bins. 
    ! Instead, at high frequencies it's fine to just go straight
    slope_lowf = (log10(pbin(2)) - log10(pbin(1))) / (log10(xbin(2)) - log10(xbin(1)))
    slope_lowf = min(slope_lowf, 0.0) ! can only increase to low frequencies

    ! rescale for logarithmic spline
    bin_spec_full_range(1,1)        = log10(x(1))
    bin_spec_full_range(2:nbin+1,1) = log10(xbin)
    bin_spec_full_range(nbin+2,1)   = log10(x(n) + 0.5)

    bin_spec_full_range(1,2)        = log10(pbin(1)) + slope_lowf * (log10(x(1)) - log10(xbin(1)))
    bin_spec_full_range(2:nbin+1,2) = log10(pbin)
    bin_spec_full_range(nbin+2,2)   = log10(pbin(nbin))

!    sigmas_full_range(1)        = sigmas(1)
!    sigmas_full_range(2:nbin+1) = sigmas
!    sigmas_full_range(nbin+2)   = sigmas(nbin)
    deallocate(xbin, pbin,nbins)

    call free_spline(self%spline_profile)
    call spline(self%spline_profile, real(bin_spec_full_range(:,1),dp), real(bin_spec_full_range(:,2),dp))
    deallocate(bin_spec_full_range)

    ! Resetting the white noise level on the lowest point
    do i = 1, n
       nu = real(log10(x(i)),dp)
       f_nu = real(10**splint(self%spline_profile, nu), sp)
       wn_level = min(f_nu, wn_level)
    end do

    ! Rebinning in fixed bins which stay common through different scans
    ! for remove_outlier routine
    do i = 1, self%n_common_points
       nu = real(log10(self%common_points(i,1)), dp)
       f_nu = real(splint(self%spline_profile, nu), sp)
       self%common_points(i,2) = f_nu
    end do

    if (self%npar /= nbin + 1) then ! spline nodes (but first and last) + sigma0
       deallocate(self%xi_n,self%P_active,self%P_uni,self%nu_fit,self%P_lognorm)
       self%npar = nbin + 1
       allocate(self%xi_n(self%npar))
       allocate(self%P_uni(self%npar,2))
       allocate(self%P_active(self%npar,2))
       allocate(self%P_lognorm(self%npar))
       allocate(self%nu_fit(self%npar,2))
    end if

    ! Sigma0
    self%xi_n(1)       = sqrt(0.99 * wn_level) ! 0.99 to avoid singularity in correlated noise
    self%sigma0 => self%xi_n(1)
    self%pcorr_offset = self%xi_n(1)**2 ! P_corr(f) = P_spline(f) - pcorr_offset, doesn't break if sigma0 is sampled

    self%P_uni(1,:)    = [0.5, 2.0] * self%xi_n(1)
    self%P_active(1,1) = self%xi_n(1)
    self%P_active(1,2) = 50.d0
    self%nu_fit(1,1)   = 5.d0 ! only fit sigma0 from 5 Hz on
    self%nu_fit(1,2)   = x(n)
    self%P_lognorm = .false.

    ! Correlated noise
    do i = 2, self%npar
       self%P_uni(i,1)    = -6.0
       self%P_uni(i,2)    = 12.0
       self%P_active(i,2) = 0 

       if (i==2) then
          self%nu_fit(i,1) = 10**self%spline_profile%x(i)
          self%nu_fit(i,2) = (10**self%spline_profile%x(i+1) + 10**self%spline_profile%x(i)) / 2
       else if (i==self%npar) then
          self%nu_fit(i,1) = (10**self%spline_profile%x(i) + 10**self%spline_profile%x(i-1)) / 2
          self%nu_fit(i,2) = 10**self%spline_profile%x(i)
       else
          self%nu_fit(i,1) = (10**self%spline_profile%x(i) + 10**self%spline_profile%x(i-1)) / 2
          self%nu_fit(i,2) = (10**self%spline_profile%x(i+1) + 10**self%spline_profile%x(i)) / 2
       end if
    end do

    !self%spline_profile%y = 10**self%spline_profile%y - self%xi_n(1)**2 ! Correlated noise only!!
    !self%xi_n(1) = max(self%xi_n(1),self%sigma0_preproc) ! avoid too low wn_level estimation from binning
    !self%sigma0 => self%xi_n(1)
    !self%spline_profile%y = log10(self%spline_profile%y)
    self%xi_n(2:) = self%spline_profile%y(2:self%npar)
    self%P_active(2:,1) = self%spline_profile%y(2:self%npar)

    deallocate(sigmas_full_range)

  end subroutine update_spline

  function eval_noise_psd_spline_full(self, nu)
    !
    ! Evaluation routine for splined noise PSD object
    !
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    !
    implicit none
    class(comm_noise_psd_spline),        intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_spline_full

    eval_noise_psd_spline_full = self%xi_n(SIGMA0)**2 + self%eval_corr(nu)

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_spline_full = eval_noise_psd_spline_full * splint(self%modulation_filter, dble(nu))
      end if
    end if

  end function eval_noise_psd_spline_full

  function eval_noise_psd_spline_corr(self, nu)
    !
    ! Evaluation routine for splined noise PSD object; correlated noise only
    !
    ! Arguments
    ! ---------
    ! self:    derived type (comm_noise_psd)
    !          Basic noise PSD object
    ! nu:      sp (scalar)
    !          Frequency (in Hz) at which to evaluate PSD
    !
    implicit none
    class(comm_noise_psd_spline),        intent(in)      :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_spline_corr

    real(sp) :: nu_, p

    nu_ = log10(max(nu, real(10.d0**self%spline_profile%x(1), sp))) ! Avoid evaluation below the spline domain
    p = 10**(splint(self%spline_profile, dble(nu_))) - self%pcorr_offset
    eval_noise_psd_spline_corr = max(p, 1.e-3 * self%pcorr_offset)  ! P_corr is strictly positive

    if(self%apply_filter) then
      if(nu >= self%modulation_filter%x(1) .and. nu <= self%modulation_filter%x(size(self%modulation_filter%x))) then
        eval_noise_psd_spline_corr = eval_noise_psd_spline_corr * splint(self%modulation_filter, dble(nu))
      end if
    end if

  end function eval_noise_psd_spline_corr

  subroutine bin_psd_spline(x, ps, thr, n_dec, dx_lin, n_min, xbin, pbin, nbin)
    !
    ! Bin a PSD for the spline noise model:
    ! - logarithmic bins (n_dec per decade) from x(1) up to thr
    ! - linear bins of width dx_lin above thr. The last one reaches x(n)
    !
    implicit none
    real(sp),     dimension(:),              intent(in)  :: x, ps
    real(sp),                                intent(in)  :: thr, dx_lin
    integer(i4b),                            intent(in)  :: n_dec, n_min
    real(sp),     allocatable, dimension(:), intent(out) :: xbin, pbin
    integer(i4b), allocatable, dimension(:), intent(out) :: nbin

    integer(i4b) :: i, k, n, nlog, nlin, nedge, nout, cnt
    real(dp)     :: logmin, logthr, s_ps, s_logx
    real(dp),     allocatable, dimension(:) :: edges, sum_ps, sum_logx
    integer(i4b), allocatable, dimension(:) :: n_out

    n = size(x)
    if (size(ps) /= n) stop 'bin_psd_spline: x and ps have different sizes'

    ! Bin edges
    logmin  = log10(real(x(1),dp))
    logthr  = log10(real(thr,dp))
    nlog  = max(1, ceiling((logthr - logmin) * n_dec))
    nlin  = max(1, floor((x(n) - thr) / dx_lin))  ! last bin absorbs the remaining frequencies
    nedge = nlog + nlin + 1
    allocate(edges(nedge))
    do k = 0, nlog
       edges(k+1) = 10.d0**(logmin + (logthr - logmin) * real(k,dp) / real(nlog,dp))
    end do
    do k = 1, nlin
       edges(nlog+1+k) = real(thr,dp) + real(k,dp) * real(dx_lin,dp)
    end do
    edges(nedge) = real(x(n),dp) * (1.d0 + 1.d-6)        ! include the last point

    allocate(sum_ps(nedge-1), sum_logx(nedge-1), n_out(nedge-1))
    nout = 0; cnt = 0; s_ps = 0.d0; s_logx = 0.d0; k = 1

    do i = 1, n
       do while (k < nedge-1)
          if (x(i) < edges(k+1)) exit
          call close_bin()
          k = k + 1
       end do
       s_ps = s_ps + real(ps(i),dp)
       s_logx = s_logx + log10(real(x(i),dp))
       cnt  = cnt + 1
    end do
    call close_bin()

    ! Leftover points (fewer than n_min) are merged into the last bin
    if (cnt > 0) then
       if (nout == 0) then
          nout = 1; sum_ps(1) = 0.d0; sum_logx(1) = 0.d0; n_out(1) = 0
       end if
       sum_ps(nout) = sum_ps(nout) + s_ps
       sum_logx(nout) = sum_logx(nout) + s_logx
       n_out(nout)  = n_out(nout)  + cnt
    end if

    allocate(xbin(nout), pbin(nout), nbin(nout))
    nbin = n_out(1:nout)
    pbin = real(sum_ps(1:nout) / nbin, sp)
    xbin = real(10.d0**(sum_logx(1:nout) / nbin), sp)
    deallocate(edges, sum_ps, sum_logx, n_out)

  contains

    subroutine close_bin()
      ! Close the bin if it has enough points, otherwise keep
      ! them for the next bin
      if (cnt >= n_min) then
         nout = nout + 1
         sum_ps(nout) = s_ps
         sum_logx(nout) = s_logx
         n_out(nout)  = cnt
         cnt = 0; s_ps = 0.d0; s_logx = 0.d0
      end if
    end subroutine close_bin

  end subroutine bin_psd_spline
                  
  function draw_bin_power(pmean, k, handle) result(p)
    !
    ! Draw the PSD node from its conditional posterior.
    ! Chi2 distribution comes from choosing 1/p as prior
    !
    implicit none
    real(sp),         intent(in)    :: pmean
    integer(i4b),     intent(in)    :: k
    type(planck_rng), intent(inout) :: handle
    real(sp)                        :: p
    integer(i4b) :: i
    real(dp)     :: chi2, z, y

    if (k <= 200) then ! sum of k squared Gaussians
       chi2 = 0.d0
       do i = 1, k
          chi2 = chi2 + rand_gauss(handle)**2
       end do
    else  ! Wilson-Hilferty approximation of chi2_k
       z = rand_gauss(handle)
       y = 1.d0 - 2.d0/(9.d0*k) + z*sqrt(2.d0/(9.d0*k))
       chi2 = k * y**3
    end if
    p = real(k * real(pmean,dp) / chi2, sp)
  end function draw_bin_power

  
end module comm_tod_noise_psd_mod
