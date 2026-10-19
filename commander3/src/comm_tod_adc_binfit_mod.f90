!===============================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University
! of Oslo.
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
!===============================================================================

! This module handles corrections to time ordered data due to adc issues

module comm_tod_adc_binfit_mod
  use math_tools
  use comm_param_mod
  use spline_1d_mod
  use locate_mod
  implicit none

  private
  public comm_adc_binfit, adc_binfit_pointer

  type :: comm_adc_binfit
     integer(i4b)       :: comm, myid, npar_adc
     character(len=128) :: label
     character(len=2048) :: datadir, outdir
     integer(i4b)       :: min_adu, max_adu, min_coadd, max_coadd, ncoadd, nbit
     integer(i4b), allocatable, dimension(:,:) :: param_adc ! (code, width, global mod/local)
     real(dp),   allocatable, dimension(:) :: p ! (npar_adc)
     real(dp),   allocatable, dimension(:) :: Q     ! ADU bin widths (fast samp)
     real(dp),   allocatable, dimension(:) :: v     ! Voltage boundaries
     real(dp),   allocatable, dimension(:) :: DNL   ! ADU bin widths (fast samp)
     real(dp),   allocatable, dimension(:) :: INL   ! Integral non-linearity
     type(spline_type)                     :: A_s   ! Noise weighted transfer functio
     type(spline_type)                     :: F     ! ADC transfer function (coadded)
     real(dp),   allocatable, dimension(:) :: invF  ! ADC correction function
     real(dp),   allocatable, dimension(:,:) :: invF_dpc ! DPC ADC correction function; two parities
  contains
    procedure :: adc_correct
    procedure :: param2Q
    procedure :: Q2As
    procedure :: As2F
    procedure :: mcmc_sample_adc
    procedure :: powell_adc
    procedure :: outputASCII
  end type comm_adc_binfit

  interface comm_adc_binfit
    procedure constructor_adc_binfit
  end interface comm_adc_binfit

  type adc_binfit_pointer
    class(comm_adc_binfit), pointer :: p => null()
  end type adc_binfit_pointer

  real(dp),     parameter :: A_s_sigma = 5.d0     ! Noise integration limit in sigma
  integer(i4b), parameter :: A_s_oversampling = 3 ! Number of evaluation points per ADU
  
interface

  module function constructor_adc_binfit(comm, datadir, outdir, label, nbit, min_adu, max_adu, ncoadd) result (c)
    ! ====================================================================
    ! Sets up an adc correction object that maps 
    !
    ! Inputs:
    ! 
    ! comm   : integer
    !          mpi communicator
    !
    ! myid   : integer
    !          mpi identifier
    !
    ! nbins  : integer
    !          number of bins used for building the adc correction tables
    !
    ! Returns:
    ! --------
    !
    ! constructor_internal: pointer
    !    contains all of the bins needed for computing adc corrections
    !    and the actual correction tables
    ! ====================================================================
    implicit none
    character(len=*),       intent(in) :: label, datadir, outdir
    integer(i4b),           intent(in) :: comm, nbit, min_adu, max_adu, ncoadd
    class(comm_adc_binfit), pointer    :: c
  end function constructor_adc_binfit

  module subroutine adc_correct(self, tod, mask, mod_phase)
    !=========================================================================
    ! Adc corrects a timestream 
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! tod : double precision array
    !    The tod that is to be corrected
    ! mask : real
    !    Flag mask
    ! 
    ! Outputs : 
    !
    ! tod_out : float array
    !    The adc corrected version of tod_in    
    ! ====================================================================
    implicit none
    class(comm_adc_binfit),          intent(in)    :: self
    real(dp), dimension(:),          intent(inout) :: tod
    logical(lgt), dimension(:),      intent(in)    :: mask
    real(sp),                        intent(in)    :: mod_phase
  end subroutine adc_correct

  module subroutine Q2As(self, sigma0)
    !=========================================================================
    ! Compute coadded transfer function from bin widths
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! sigma0 : real
    !    White noise level (coadded)
    ! 
    ! Outputs : 
    ! ====================================================================
    implicit none
    class(comm_adc_binfit),          intent(inout) :: self
    real(sp),                        intent(in)    :: sigma0
  end subroutine Q2AS

  module subroutine As2F(self, gain, offset)
    !=========================================================================
    ! Compute coadded transfer function given (optional) fast gain and offset
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! gain : real
    !    40-sample vector with gain; should have unity mean
    ! offset: real
    !    40-sample vector with offset in ADU; for HFI, this is typically 4K lines at >180 Hz
    ! 
    ! Outputs : 
    ! ====================================================================
    implicit none
    class(comm_adc_binfit),  intent(inout) :: self
    real(sp), dimension(40), intent(in), optional :: gain
    real(sp), dimension(40), intent(in), optional :: offset
  end subroutine As2F

  module subroutine param2Q(self, p)
    !=========================================================================
    ! 
    ! 
    ! Inputs:
    !
    ! self:     comm_adc object
    !           Base ADC object
    ! p:        real(sp)
    !           Parameter vector; must be of length npar_adc
    implicit none
    class(comm_adc_binfit),             intent(inout) :: self
    real(dp), dimension(self%npar_adc), intent(in), optional    :: p
  end subroutine param2Q

  module subroutine powell_adc(self, powell_chisq_adc)
    !=========================================================================
    ! Sample bin widths for single scan; coadd into posterior
    ! 
    ! Inputs:
    !
    ! self:     comm_adc object
    !           Base ADC object
    implicit none
    class(comm_adc_binfit),          intent(inout) :: self
    interface
       function powell_chisq_adc(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: powell_chisq_adc
       end function powell_chisq_adc
    end interface
  end subroutine powell_adc
  
  module subroutine mcmc_sample_adc(self, handle, chisq_adc)
    !=========================================================================
    ! Sample bin widths for single scan; coadd into posterior
    ! 
    ! Inputs:
    !
    ! self:     comm_adc object
    !           Base ADC object
    implicit none
    class(comm_adc_binfit),          intent(inout) :: self
    type(planck_rng),                intent(inout) :: handle
    interface
       function chisq_adc(p, ndof)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in),  optional :: p
         integer(i8b),           intent(out), optional :: ndof
         real(dp)                                      :: chisq_adc
       end function chisq_adc
    end interface
  end subroutine mcmc_sample_adc

  module subroutine outputASCII(self)
    !=========================================================================
    ! Sample bin widths for single scan; coadd into posterior
    ! 
    ! Inputs:
    !
    ! self:     comm_adc object
    !           Base ADC object
    implicit none
    class(comm_adc_binfit),          intent(inout) :: self
  end subroutine outputASCII

  
end interface

end module comm_tod_adc_binfit_mod
