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

module comm_tod_adc_mod

  use math_tools
  use spline_1D_mod
  use comm_tod_mod
  use comm_map_mod
  use comm_param_mod
  use InvSamp_mod
  use comm_hdf_mod
  implicit none

  private
  public comm_adc, adc_pointer

  type :: comm_adc
    real(sp),          allocatable, dimension(:) :: adc_in, adc_out, rms_bins, rms2_bins, v_bins, vbin_edges
    real(sp),          allocatable, dimension(:) :: dpc_in, dpc_out, err_bins
    integer(i4b),      allocatable, dimension(:) :: nval
    integer(i4b)                                 :: comm, myid, nbins, window
    real(sp)                                     :: v_min, v_max ! Global variable for the experiment determined in the parameter file 
    class(comm_mapinfo), pointer                 :: info => null()    ! Map definition
    character(len=512)                           :: outdir
    type(spline_type)                            :: sadc

    real(sp),          allocatable, dimension(:) :: rms_bins2
    integer(i4b),      allocatable, dimension(:) :: nval2

  contains
    procedure :: adc_correct
    procedure :: build_table
    procedure :: bin_scan_rms
    procedure :: find_horn_min_max
    procedure :: construct_voltage_bins
    
    procedure :: corr_rms_out
  end type comm_adc

  interface comm_adc
    procedure constructor_internal, constructor_precomp
  end interface comm_adc

  type adc_pointer
    class(comm_adc), pointer :: p => null()
  end type adc_pointer

interface

  !=========================================================================
  ! This code is based off of the DPC ADC corrections (thanks Bob).
  ! 
  ! Those corrections procede as such:
  !
  ! The detector has some characteristic response function V' = R(V)V
  ! which is assumed to be near unity everywhere.
  !
  ! As the sky signal changes, the change in the measured detector
  ! voltage is given by 
  ! 
  !             \deltaV' = (V * dR(V)/dV + R(V) )* \deltaV
  !
  ! We aim to find the derivative dR(V)/dV of the response function to 
  ! correct for the non-linearties described in 
  ! Planck 2013: III https://arxiv.org/abs/1303.5064
  !
  ! DCH
  !=========================================================================

  module function constructor_internal(cpar, info, nbins) result(res)
    ! ====================================================================
    ! Sets up an adc correction object that maps input and output voltages
    ! Also initializes the bins used for the actual correction model
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
    integer(i4b),           intent(in) :: nbins
    class(comm_mapinfo),    target     :: info
    class(comm_adc),        pointer    :: res
    type(comm_params),      intent(in) :: cpar

  end function constructor_internal

  module function constructor_precomp(instfile, path, load) result(res)
    ! ====================================================================
    ! Sets up an adc correction object that maps input and output voltages
    ! Also initializes the bins used for the actual correction model
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
    ! constructor: pointer
    !    contains all of the bins needed for computing adc corrections
    !    and the actual correction tables
    ! ====================================================================
    
    implicit none
    type(hdf_file),     intent(in) :: instfile
    character(len=512), intent(in) :: path
    logical(lgt),       intent(in) :: load
    class(comm_adc),    pointer    :: res
  end function constructor_precomp

  module subroutine adc_correct(self, tod_in, tod_out, scan, det, di)
    !=========================================================================
    ! Adc corrects a timestream 
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! tod_in : float array
    !    The tod that is to be corrected
    ! scan   : integer
    !    tod scan number
    ! det    : integer
    !    instrument detector identifier
    ! di     : integer
    !    detector diode identifier
    ! 
    ! Outputs : 
    !
    ! tod_out : float array
    !    The adc corrected version of tod_in    
    ! ====================================================================
    implicit none
    class(comm_adc),                 intent(inout) :: self
    real(sp), dimension(:),          intent(in)    :: tod_in
    real(sp), dimension(:),          intent(out)   :: tod_out
    integer(i4b), intent(in), optional :: scan, det, di
  end subroutine adc_correct

  module subroutine build_table(self,handle,name)
    !=========================================================================
    ! Adc corrects a timestream by fitting regular Gaussian dips to the binned
    ! RMS values of the 
    ! 
    ! Inputs:
    !
    ! self:     comm_adc object
    !           Defines the adc correction that should be applied
    ! handle:   type(planck_rng)
    !           Healpix random number type
    ! name:     string
    !           diode name for output file names
    ! Outputs: 
    !
    ! volt_in:  float array
    !           Array of the input voltages
    !
    ! volt_out: float
    !           array of the corrected voltages
    implicit none
    class(comm_adc),                 intent(inout) :: self
    type(planck_rng),                intent(inout) :: handle
    character(len=50),               intent(in)    :: name
  end subroutine build_table

  module subroutine construct_voltage_bins(self)

    implicit none
    class(comm_adc), intent(inout) :: self
  end subroutine construct_voltage_bins
  
  module subroutine find_horn_min_max(self,tod_in,flag,flag0)
    ! ==================================================================
    ! This subroutine loops through the diode data to find the global
    ! minimum and maximum voltages
    !
    ! Inputs:
    !
    ! self     : comm_adc object
    !    Defines the adc correction that should be applied
    !
    ! tod_in   : float array (sp)
    !    The chunk of data which we will be looping over
    ! 
    ! flag     : integer array
    !    data flagging corresponding to the chunk of data 
    !
    ! flag0    : integer
    !    something I don't really know what it does
    ! 
    ! "Outputs":
    !
    ! self%v_min : float (sp)
    !    updated global minimum voltage
    !
    ! self%v_max : float (sp)
    !    updated global maximum voltage
    !
    ! ====================================================================
    implicit none
    class(comm_adc),                   intent(inout) :: self
    real(sp),     dimension(:),        intent(in)    :: tod_in
    integer(i4b), dimension(:),        intent(in)    :: flag
    integer(i4b),                      intent(in)    :: flag0
  end subroutine find_horn_min_max

  module subroutine bin_scan_rms(self,tod_in,flag,flag0,corr)
    ! ====================================================================
    ! This subroutine takes in a chunk of data (from wherever it lives) and puts
    ! it in the appropriate bins
    !
    ! Inputs:
    !
    ! self     : comm_adc object
    !    Defines the adc correction that should be applied
    !
    ! tod_in   : float array (sp)
    !    The chunk of data which we will be binning here
    !
    !
    ! "Outputs":
    !
    !   self%nval: integer array
    !     Counts the number of entries in each bin
    ! 
    !   self%binval: float array (sp)
    !     Sum of all rms values for each bin
    ! ====================================================================
    
    implicit none
    class(comm_adc),               intent(inout) :: self
    real(sp),     dimension(:),    intent(in)    :: tod_in
    integer(i4b), dimension(:),    intent(in)    :: flag
    integer(i4b),                  intent(in)    :: flag0

    logical(lgt),      optional,   intent(in)    :: corr ! follow slightly different procedure if data is corrected
  end subroutine bin_scan_rms

  module subroutine corr_rms_out(self,name)
    ! ====================================================================
    !
    ! Test function to output the rms of the corrected TODs
    !
    ! Input: 
    !
    ! self : class
    !        comm_adc
    ! name : character
    !        diode name
    ! ====================================================================
    implicit none
    class(comm_adc),   intent(inout) :: self
    character(len=50), intent(in)    :: name
  end subroutine corr_rms_out

  module subroutine mask_bins(vbins,rms,nval,mask)
    ! ====================================================================
    ! This subroutine iterates masks out bins with spuriously large deviations
    ! in the white noise level (spikes) and bins with 0 entries
    !
    ! Inputs:
    !
    ! vbins:     float array
    !            the 'x-axis': the voltage bins
    !
    ! rms:       float array
    !            the 'y-axis': the rms bins
    !
    ! nval:      integer array
    !            the count of entries for each rms bin
    !
    ! mask:      integer array
    !            binary array which determines which bins contribute to the fitting procedure
    !
    ! Outputs:
    !
    ! mask:      integer array
    !            binary array which determines which bins contribute to the fitting procedure
    ! ====================================================================
    implicit none
    real(sp),     dimension(:), intent(in)    :: vbins
    integer(i4b), dimension(:), intent(in)    :: nval
    real(sp),     dimension(:), intent(inout) :: rms
    integer(i4b), dimension(:), intent(inout) :: mask
  end subroutine mask_bins
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! DP FUNCTIONS !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  module subroutine return_linreg_dp(x,y,mask,slope,offset,trim)
    !=========================================================================
    ! Very simple function that fits the slope and offset for an x-y pair
    !
    ! Inputs:
    !
    ! x      : float array (sp)
    !          The independent variable
    ! y      : float array (sp)
    !          The dependent variable
    ! mask   : integer array
    !          which array values do not contribute to the fit? 
    ! trim   : logical (optional)
    !          are large deviations disregarded in the linear fit?
    !
    ! Output:
    !
    ! slope  : float (sp)
    !          The slope of the linear fit
    ! offset : float (sp)
    !          The offset of the linear fit
    !=========================================================================    
    implicit none
    
    real(dp),     dimension(:),  intent(in)    :: x, y
    integer(i4b), dimension(:),  intent(in)    :: mask
    logical(lgt), optional,      intent(in)    :: trim
    real(dp),                    intent(inout) :: slope, offset
  end subroutine return_linreg_dp

  module subroutine return_dips_dp(x,y,mask,diprange,res,name)
    ! ====================================================================
    ! Dip identifying subroutine that returns the first dip location and
    ! the distance between dips in the RMS
    !
    ! Inputs:
    !
    ! x       : float array
    !            input voltage axis
    ! 
    ! y       : float array
    !           binned rms array with the linear part removed
    !
    ! mask    : integer array
    !           mask for binned rms array - 0 if y(i) = NaN or has been trimmed
    !
    ! diprange: integer
    !           how many x axis indices we want to search for dips within
    !
    ! Outputs:
    !
    ! dips    : integer array
    !           array of the dip locations in the voltage bin array
    ! ====================================================================
    implicit none
    
    real(dp),     dimension(:), intent(in)    :: x, y
    integer(i4b), dimension(:), intent(inout) :: mask
    integer(i4b),               intent(in)    :: diprange
    character(len=50),          intent(in)    :: name
    integer(i4b), dimension(:), allocatable   :: res
  end subroutine return_dips_dp

  module function return_gaussian_dp(x, pars) result(y)
    ! ====================================================================
    ! Super simple function which returns a guassian function given the parameters
    !
    ! Inputs: 
    !
    ! x     : float array
    !         the x-array over which the gaussian is evaluated
    ! pars  : float array (size = 3)
    !         mean, stddev, and amplitude of the gaussian
    !
    ! Outputs:
    !
    ! y     : float array
    !         returned gaussian function
    ! ====================================================================
    
    implicit none
    
    real(dp), dimension(:), intent(in)         :: x
    real(dp), dimension(3), intent(in)         :: pars
    real(dp), allocatable, dimension(:)        :: y
  end function return_gaussian_dp

  module function return_gauss_lin_model_dp(x, y, mask, bincount, a, b, dips, handle, name) result(model)
    !=========================================================================
    ! Fits a gaussian function to each recognized dip in the 
    ! white noise level
    !
    ! Inputs:
    ! 
    ! x:       float array
    !          voltage bins from the tod_in binning
    ! y:       float array
    !          rms bins from the rms estimates as a function of voltage
    ! mask:    integer array
    !          mask for binned rms array - 0 if y(i) = NaN or has been trimmed
    ! dips:    integer array
    !          index value for the location of the dips to be fit
    ! handle:  type(planck_rng)
    !          Healpix random number type
    !
    ! Outputs: 
    !
    ! idrf:    float array 
    !          array of length(x) - the inverse differential response function
    !          which will be integrated
    !=========================================================================
    
    implicit none
    
    real(dp),     dimension(:), intent(in)    :: x, y
    real(dp),                   intent(inout) :: a, b ! slope and offset of the model fit
    integer(i4b), dimension(:), intent(in)    :: mask
    integer(i4b), dimension(:), intent(in)    :: bincount
    integer(i4b), dimension(:), intent(in)    :: dips
    type(planck_rng),           intent(inout) :: handle
    real(dp),     dimension(:), allocatable   :: model
    character(len=50),          intent(in)    :: name

  ! contains

  !   ! Grid out and solve for maximum likelihood parameter value
  !   function find_maxlike_gauss_par(par_i) result(gpar)
  !     use healpix_types
  !     implicit none
  !     ! real(sp),   intent(inout) :: gpar
  !     integer(i4b), intent(in)  :: par_i
  !     real(dp), dimension(1000) :: lnL, grid
  !     real(dp)                  :: gpar, tmp
  !     integer(i4b)              :: l, i, ind
  !     real(dp), dimension(last-first) :: gauss

  !     lnL(:) = 0.0

  !     do l = 1, 1000
  !        grid(l) = (P_uni(2)-P_uni(1))*(l-1)/1000 + P_uni(1)
         
  !        tmp = pars(par_i)
  !        pars(par_i) = grid(l)

  !        gauss = return_gaussian_dp(x_tmp,pars)
  !        do i = 1, last-first
  !           if (mask_tmp(i) == 0) cycle
  !           lnL(l) = lnL(l) - 0.50*(y_tmp(i) - gauss(i))**2/count_tmp(i)
  !        end do
         
  !     end do
      
  !     ind = maxloc(lnL,dim=1)

  !     gpar = grid(ind)

  !   end function find_maxlike_gauss_par
    
  end function return_gauss_lin_model_dp
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! SP FUNCTIONS !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine return_linreg_sp(x,y,mask,slope,offset,trim)
    !=========================================================================
    ! Very simple function that fits the slope and offset for an x-y pair
    !
    ! Inputs:
    !
    ! x      : float array (sp)
    !          The independent variable
    ! y      : float array (sp)
    !          The dependent variable
    ! mask   : integer array
    !          which array values do not contribute to the fit? 
    ! trim   : logical (optional)
    !          are large deviations disregarded in the linear fit?
    !
    ! Output:
    !
    ! slope  : float (sp)
    !          The slope of the linear fit
    ! offset : float (sp)
    !          The offset of the linear fit
    !=========================================================================    
    implicit none
    
    real(sp),     dimension(:),  intent(in)    :: x, y
    integer(i4b), dimension(:),  intent(in)    :: mask
    logical(lgt), optional,      intent(in)    :: trim
    real(sp),                    intent(inout) :: slope, offset
  end subroutine return_linreg_sp

end interface

end module comm_tod_adc_mod
