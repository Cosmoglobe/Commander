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
module comm_camb_eval_mod
  use CAMB
  implicit none

  private
  public get_c_l_from_camb

contains


  subroutine get_c_l_from_camb(params, Cl)
    ! 
    ! Gets TT, EE, and TE power spectra from camb using the cosmological
    ! parameters in theta.
    !
    ! Arguments
    ! ---------
    ! params: List of CAMB parameters
    !
    ! Returns
    ! -------
    ! Cl: Power spectrum
    ! 
    implicit none
    real(8), dimension(1:),    intent(in) :: params
    real(8), dimension(1:,0:), intent(out) :: Cl
    
    !class(comm_camb),                           intent(inout) :: self
    !type(comm_camb_sample) :: new_sample 
    real(8), dimension(6) :: cosmo_param
    integer(4) :: l, k, lmax
    
    type(CAMBparams) P
    type(CAMBdata) camb_data
    real(8) :: CMB_outputscale, pi

    !cosmo_param = new_sample%theta
    call CAMB_SetDefParams(P)

    pi      = 3.141592653589793238462643383279502884197_8
    lmax    = size(Cl,2)-1
    
    P%ombh2 = cosmo_param(1)
    P%omch2 = cosmo_param(2)
    P%omk = 0.d0
    P%H0= cosmo_param(3)
    !call P%set_H0_for_theta(0.0104d0)!cosmo_param(4)
    select type(InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
       InitPower%As = exp(cosmo_param(5))*1e-10
       InitPower%ns = cosmo_param(6)
       InitPower%r = 0.0
    end select
    
    select type(Reion=>P%Reion)
    class is (TTanhReionization)
       Reion%use_optical_depth = .true.
       Reion%optical_depth = cosmo_param(4)
    end select
    
    P%WantScalars           = .true.
    P%WantTensors           = .true.
    P%Accuracy%AccurateBB   = .true.
    P%WantCls               = .true.
    P%DoLensing             = .false.

    P%Max_l=lmax+2
    P%Max_eta_k=6000
    P%Max_l_tensor=lmax+2
    P%Max_eta_k_tensor=6000

    ! From the CAMB documentation you need this to get micro K^2.
    ! This is (2.726 K * 10^6)^2
    CMB_outputscale = 7.4311e12
    
    call CAMB_GetResults(camb_data, P)
    
    ! Set TT, EE, and TE
    !c_l = 0.d0
    !write(*,*) 'camb', shape(camb_data%CLData%Cl_scalar)
    !write(*,*) 'camb k=1', camb_data%CLData%Cl_scalar(2:4, 1)*CMB_outputscale
    !write(*,*) 'camb k=2', camb_data%CLData%Cl_scalar(2:4, 2)*CMB_outputscale
    !write(*,*) 'camb k=3', camb_data%CLData%Cl_scalar(2:4, 3)*CMB_outputscale
    !write(*, *) 'my cl', shape(c_l)
    DO k = 1, 3
       Cl(k, 0) = 0.d0
       Cl(k, 1) = 0.d0
       DO l = 2, lmax
          Cl(k, l) = 2.d0 * pi / (l * (l + 1)) * camb_data%CLData%Cl_scalar(l, k) * CMB_outputscale
       END DO
    END DO
    !new_sample%c_l = c_l
    
  end subroutine get_c_l_from_camb


end module comm_camb_eval_mod
