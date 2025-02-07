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
module comm_bp_mod
  use comm_param_mod
  use comm_bp_utils
  implicit none 

  type :: comm_bp
     ! Data variables
     character(len=512) :: type, model
     logical(lgt)       :: sample_bandpass
     integer(i4b)       :: n, npar
     real(dp)           :: threshold
     real(dp)           :: nu_c, a2t, f2t, a2sz, unit_scale, nu_eff, a2f
     real(dp), allocatable, dimension(:) :: nu0, nu, tau0, tau, delta, a2f_arr
   contains
     ! Data procedures
     procedure     :: update_tau
     procedure     :: SED2F
     procedure     :: lineAmp_RJ
  end type comm_bp

  type comm_bp_ptr
     class(comm_bp), pointer :: p => null()
  end type comm_bp_ptr


  interface comm_bp
     procedure constructor_bp
  end interface comm_bp

  real(dp) :: ind_iras = 1.d0
  
contains

  subroutine initialize_bp_mod(cpar)
    implicit none

    type(comm_params), intent(in) :: cpar

    
  end subroutine initialize_bp_mod
  
  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_bp(cpar, id, id_abs, detlabel, subdets) result(c)
    !
    ! Initialization routine (constructor) for bandpass objects. Reads in bandpass 
    ! data, and precomputes default unit conversions etc.
    !
    ! Arguments:
    ! ----------
    ! cpar: type(comm_params)
    !    Commander parameter structure
    ! id: int (scalar)
    !    Frequency channel ID/counter, only counting active bands
    ! id_abs: int (scalar)
    !    Frequency channel ID/counter, counting all bands (as defined in parameter file)
    ! detlabel: string (scalar; optional)
    !    Detector label; typically used for full-frequency bands 
    ! subdets: string (scalar; optional)
    !    Comma-separated string with sub-detector labels. Used for TOD-type bands
    !
    ! Returns:
    ! --------
    ! constructor: class(comm_bp)
    !    Pointer to new comm_bp object
    ! 
    implicit none
    type(comm_params),                intent(in)           :: cpar
    integer(i4b),                     intent(in)           :: id, id_abs
    character(len=*),                 intent(in), optional :: detlabel, subdets
    class(comm_bp),     pointer                            :: c

    integer(i4b)       :: i, j, ndet
    character(len=512) :: label
    character(len=25)  :: dets(1500)
    real(dp), allocatable, dimension(:) :: nu0, tau0
    
    label = cpar%ds_label(id_abs)
    
    ! General parameters
    allocate(c)
    
    
  end function constructor_bp
  

  
  subroutine update_tau(self, delta)
    implicit none

    class(comm_bp),                       intent(inout) :: self
    real(dp),       dimension(self%npar), intent(in)    :: delta
    real(dp), allocatable, dimension(:)  :: a, bnu_prime, bnu_prime_RJ, sz

    integer(i4b) :: i, n


  end subroutine update_tau

  function SED2F(self, f)
    implicit none

    ! Implementation of the mixing matrix (SED2F = M, f = frequency scaling of a component).
    ! Depending on the units of the bandpass a different function which converts from K_RJ sed to DATA units bandpass integrated.
    ! See BP1 footnote 7 for K_RJ -> MJy/sr (additional factor of 2.d0*k_B*self%nu**2/c**2).

    class(comm_bp),               intent(in) :: self
    real(dp),       dimension(:), intent(in) :: f
    real(dp)                                 :: SED2F

    integer(i4b) :: i, j
    real(dp)     :: a2f1, a2f2, a2fc, K, Inu0

    SED2F = 0d0

  end function SED2F

  function lineAmp_RJ(self, nu)
    implicit none

    class(comm_bp), intent(in) :: self
    real(dp),       intent(in) :: nu
    real(dp)                   :: lineAmp_RJ

    integer(i4b) :: i
    real(dp)     :: x, tau
    real(dp), allocatable, dimension(:) :: bnu_prime_RJ, a2t


    lineAmp_RJ = 0d0

  end function lineAmp_RJ


end module comm_bp_mod
