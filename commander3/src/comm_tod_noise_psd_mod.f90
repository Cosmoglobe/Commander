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
  public comm_noise_psd, comm_noise_psd_2oof

  type :: comm_noise_psd
     ! 
     ! Class definition for basic 1/f noise PSD model
     !
     real(sp) :: sigma0,     alpha,     fknee                 ! Main dynamic/sampling parameters
     real(sp) :: sigma0_def, alpha_def, fknee_def             ! Default parameters/priors
     real(dp), allocatable, dimension(:)    :: log_n_psd      ! Noise power spectrum density; in uncalibrated units
     real(dp), allocatable, dimension(:)    :: log_n_psd2     ! Second derivative (for spline)
     real(dp), allocatable, dimension(:)    :: log_nu         ! Noise power spectrum bins; in Hz
   contains
     procedure :: eval   => eval_noise_psd
     procedure :: update => update_noise_psd
     procedure :: xi_n   => get_xi_n_noise_psd
  end type comm_noise_psd

  type, extends(comm_noise_psd) :: comm_noise_psd_2oof
     ! 
     ! Class definition for 2-component 1/f noise PSD model
     !
     real(sp) :: alpha2,     fknee2
     real(sp) :: alpha2_def, fknee2_def
   contains
     procedure :: eval   => eval_noise_psd_2oof
     procedure :: update => update_noise_psd_2oof
     procedure :: xi_n   => get_xi_n_noise_psd_2oof
  end type comm_noise_psd_2oof

  interface comm_noise_psd
     procedure constructor_oof
  end interface comm_noise_psd

  interface comm_noise_psd_2oof
     procedure constructor_2oof
  end interface comm_noise_psd_2oof

contains

  function constructor_oof(xi_n, xi_n_def)
    ! 
    ! Constructor for basic 1/f noise PSD object, where
    !     
    !     P(nu) = sigma0^2 * (1 + (nu/fknee)^alpha)
    ! 
    ! Arguments
    ! --------- 
    ! xi_n:    sp (array)
    !          3-element array containing {sigma0, alpha, fknee}, where
    !          [sigma0] = du/volts/tod unit, [alpha] = 1, and [fknee] = Hz
    ! xi_n_def: sp (array)
    !          3-element array containing default parameters/prior
    ! 
    implicit none
    real(sp),              dimension(:), intent(in)      :: xi_n
    real(sp),              dimension(:), intent(in)      :: xi_n_def
    class(comm_noise_psd), pointer                       :: constructor_oof

    allocate(constructor_oof)
    constructor_oof%sigma0     = xi_n(1)
    constructor_oof%fknee      = xi_n(2)
    constructor_oof%alpha      = xi_n(3)
    constructor_oof%sigma0_def = xi_n_def(1)
    constructor_oof%fknee_def  = xi_n_def(2)
    constructor_oof%alpha_def  = xi_n_def(3)

  end function constructor_oof

  function eval_noise_psd(self, nu)
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
    class(comm_noise_psd),               intent(inout)   :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd

    eval_noise_psd = self%sigma0**2 * (1. + (nu/self%fknee)**self%alpha)

  end function eval_noise_psd

  subroutine update_noise_psd(self, xi_n)
    ! 
    ! Routine to update parameters in basic 1/f noise PSD object
    ! 
    ! Arguments
    ! --------- 
    ! self:    derived class (comm_noise_psd)
    !          Object to be updated
    ! xi_n:    sp (array)
    !          3-element array containing {sigma0, fknee, alpha}, where
    !          [sigma0] = du/volts/tod unit, [alpha] = 1, and [fknee] = Hz
    !
    implicit none
    class(comm_noise_psd),               intent(inout)   :: self
    real(sp),              dimension(:), intent(in)      :: xi_n

    self%sigma0     = xi_n(1)
    self%fknee      = xi_n(2)
    self%alpha      = xi_n(3)

  end subroutine update_noise_psd

  subroutine get_xi_n_noise_psd(self, xi_n)
    ! 
    ! Routine to return parameters in basic 1/f noise PSD object
    ! 
    ! Arguments
    ! --------- 
    ! self:    derived class (comm_noise_psd)
    !          Object to be updated
    ! 
    ! Returns
    ! -------
    ! xi_n:    sp (array)
    !          3-element array containing {sigma0, fknee, alpha}
    !
    implicit none
    class(comm_noise_psd),               intent(in)      :: self
    real(sp),              dimension(:), intent(out)     :: xi_n

    xi_n(1) = self%sigma0
    xi_n(2) = self%fknee
    xi_n(3) = self%alpha

  end subroutine get_xi_n_noise_psd



  function constructor_2oof(xi_n, xi_n_def)
    ! 
    ! Constructor for two-component 1/f noise PSD object, where
    !     
    !     P(nu) = sigma0^2 * (1 + (nu/fknee)^alpha + (nu/fknee2)^alpha2)
    ! 
    ! Arguments
    ! ---------
    ! xi_n:    sp (array)
    !          5-element array containing {sigma0, fknee, alpha, fknee2, alpha2}, where
    !          [sigma0] = du/volts/tod unit, [alpha,alpha2] = 1, and [fknee,fknee2] = Hz
    ! xi_n_def: sp (array)
    !          5-element array containing default parameters/prior values
    ! 
    implicit none
    real(sp),                   dimension(:), intent(in)      :: xi_n
    real(sp),                   dimension(:), intent(in)      :: xi_n_def
    class(comm_noise_psd_2oof), pointer                       :: constructor_2oof

    allocate(constructor_2oof)
    constructor_2oof%sigma0     = xi_n(1)
    constructor_2oof%fknee      = xi_n(2)
    constructor_2oof%alpha      = xi_n(3)
    constructor_2oof%fknee2     = xi_n(4)
    constructor_2oof%alpha2     = xi_n(5)
    constructor_2oof%sigma0_def = xi_n_def(1)
    constructor_2oof%fknee_def  = xi_n_def(2)
    constructor_2oof%alpha_def  = xi_n_def(3)
    constructor_2oof%fknee2_def = xi_n_def(4)
    constructor_2oof%alpha2_def = xi_n_def(5)

  end function constructor_2oof
  
  function eval_noise_psd_2oof(self, nu)
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
    class(comm_noise_psd_2oof),          intent(inout)   :: self
    real(sp),                            intent(in)      :: nu
    real(sp)                                             :: eval_noise_psd_2oof

    eval_noise_psd_2oof = self%sigma0**2 * (1. + (nu/self%fknee)**self%alpha + (nu/self%fknee2)**self%alpha2)

  end function eval_noise_psd_2oof

  subroutine update_noise_psd_2oof(self, xi_n)
    ! 
    ! Routine to update parameters in two-component 1/f noise PSD object
    ! 
    ! Arguments
    ! --------- 
    ! self:    derived class (comm_noise_psd_2oof)
    !          Object to be updated
    ! xi_n:    sp (array)
    !          5-element array containing {sigma0, fknee, alpha, fknee2,alpha2}, where
    !          [sigma0] = du/volts/tod unit, [alpha] = 1, and [fknee] = Hz
    !
    implicit none
    class(comm_noise_psd_2oof),               intent(inout)   :: self
    real(sp),                   dimension(:), intent(in)      :: xi_n

    self%sigma0     = xi_n(1)
    self%fknee      = xi_n(2)
    self%alpha      = xi_n(3)
    if (size(xi_n) == 5) then
       self%fknee2     = xi_n(4)
       self%alpha2     = xi_n(5)
    end if

  end subroutine update_noise_psd_2oof

  subroutine get_xi_n_noise_psd_2oof(self, xi_n)
    ! 
    ! Routine to return parameters in 2-component 1/f noise PSD object
    ! 
    ! Arguments
    ! --------- 
    ! self:    derived class (comm_noise_psd)
    !          Object to be updated
    ! 
    ! Returns
    ! -------
    ! xi_n:    sp (array)
    !          5-element array containing {sigma0, fknee, alpha, fknee2, alpha2}
    !
    implicit none
    class(comm_noise_psd_2oof),               intent(in)      :: self
    real(sp),                   dimension(:), intent(out)     :: xi_n

    xi_n(1) = self%sigma0
    xi_n(2) = self%fknee
    xi_n(3) = self%alpha
    xi_n(4) = self%fknee2
    xi_n(5) = self%alpha2

  end subroutine get_xi_n_noise_psd_2oof


end module comm_tod_noise_psd_mod
