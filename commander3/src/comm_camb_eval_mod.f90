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
  use NonLinear
  implicit none

  private
  public get_c_l_from_camb

contains


  subroutine get_c_l_from_camb(cosmo_param, Cl)
    ! 
    ! Gets TT, EE, and TE power spectra from camb using the cosmological
    ! parameters in theta.
    !
    ! Arguments
    ! ---------
    ! cosmo_param: List of CAMB parameters
    !
    ! Returns
    ! -------
    ! Cl: Power spectrum
    ! 
    implicit none
    real(8), dimension(1:,0:), intent(out) :: Cl
    real(8), dimension(6), intent(in) :: cosmo_param
    
    !class(comm_camb),                           intent(inout) :: self
    !type(comm_camb_sample) :: new_sample 

    integer(4) :: l, k, lmax
    logical :: param_file_read
    type(CAMBparams) P
    type(CAMBdata) camb_data
    real(8) :: CMB_outputscale, pi

    !cosmo_param = new_sample%theta
    call CAMB_SetDefParams(P)
    param_file_read = CAMB_ReadParamFile(P, 'params.ini', 3000)
    write(*,*) 'reading param file:', param_file_read

    pi      = 3.141592653589793238462643383279502884197_8
    lmax    = size(Cl,2)-1
    
    P%ombh2 = cosmo_param(1)
    P%omch2 = cosmo_param(2)
    P%omk = 0.d0
    P%H0 = cosmo_param(3)

    select type(InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
       InitPower%As = cosmo_param(5)*1e-9
       InitPower%ns = cosmo_param(6)
       InitPower%r = 0.0
    end select
    
    select type(Reion=>P%Reion)
    class is (TTanhReionization)
       Reion%use_optical_depth = .true.
       Reion%optical_depth = cosmo_param(4)
    end select

    !P%WantScalars           = .true.
    !P%WantTensors           = .true.
    !P%WantCls               = .true.
    !P%WantTransfer          = .true.
    !P%DoLensing             = .true.
    !P%Want_CMB_lensing      = .true.
    !P%NonLinear = NonLinear_both
    !P%NonLinear             = 3
    !P%halofit_version       = 3

    !P%Max_l = 1000
    !P%Max_eta_k = 18000
    !P%Max_l_tensor = 1000
    !P%Max_eta_k_tensor = 6000


    !P%Accuracy%AccuracyBoost = 2
    !P%Accuracy%lAccuracyBoost = 2
    !P%Accuracy%LensingBoost = 2
    !P%Accuracy%lSampleBoost = 2
    !P%Accuracy%AccurateReionization = .true.
    !P%Accuracy%AccurateBB   = .true.
    !P%Accuracy%AccuratePolarization = .true.


    ! From the CAMB documentation you need this to get micro K^2.
    ! This is (2.726 K * 10^6)^2
    CMB_outputscale = 7.4311e12
    
    call CAMB_GetResults(camb_data, P)
    
    ! Set TT, EE, TE and BB
    DO k = 1, 4
       Cl(k, 0) = 0.d0
       Cl(k, 1) = 0.d0
       DO l = 2, lmax
          !if (k == 3) then ! We are doing scalar, so no BB
          !  Cl(k, l) = 1.0e-10
          !else if (k < 3) then!Cl_Lensing
          Cl(k, l) = 2.d0 * pi / (l * (l + 1)) * camb_data%CLData%Cl_Lensed(l, k) * CMB_outputscale!+camb_data%CLData%Cl_tensor(l, k)) * CMB_outputscale
          !else
          !  Cl(k, l) = 2.d0 * pi / (l * (l + 1)) * camb_data%CLData%Cl_Scalar(l, k-1) * CMB_outputscale
          !end if    
        END DO
    END DO
    !write(*,*) 'C^TT ell=100', Cl(1, 100), Cl(1, 100)*100*101/(2.d0*pi)
    open(unit=1, file='cosmo_cl_out.dat', position="append", action='write')
      write(1, '( 6(2X, ES14.6) )') Cl(3, :)
    close(1)
    write(*,*) 'cl written out!!!'
    
  end subroutine get_c_l_from_camb


end module comm_camb_eval_mod
