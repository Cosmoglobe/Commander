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

  implicit none

  private
  public comm_adc, adc_pointer

  type :: comm_adc
    real(sp),          allocatable, dimension(:) :: adc_in, adc_out, rms_bins, v_bins, vbin_edges
    real(sp),          allocatable, dimension(:) :: dpc_in, dpc_out
    integer(i4b),      allocatable, dimension(:) :: nval
    integer(i4b)                                 :: comm, myid, nbins, window
    real(sp)                                     :: v_min, v_max ! Global variable for the experiment determined in the parameter file 
    class(comm_mapinfo), pointer                 :: info => null()    ! Map definition
    character(len=512)                           :: outdir
    type(spline_type)                            :: sadc
  contains
    procedure :: adc_correct
    procedure :: build_table
    procedure :: bin_scan_rms
    procedure :: load_dpc
    procedure :: find_horn_min_max
    procedure :: construct_voltage_bins
  end type comm_adc

  interface comm_adc
    procedure constructor
  end interface comm_adc

  type adc_pointer
    class(comm_adc), pointer :: p => null()
  end type adc_pointer

contains

  function constructor(cpar, info, nbins)
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
    integer(i4b),           intent(in) :: nbins
    class(comm_mapinfo),    target     :: info
    class(comm_adc),        pointer    :: constructor
    type(comm_params),      intent(in) :: cpar
    
    real(sp)     :: diff

    integer(i4b) :: i, ierr
    
    allocate(constructor)
    
    constructor%window  =  10
    constructor%info    => info
    constructor%myid    =  cpar%myid_chain
    constructor%comm    =  cpar%comm_chain
    constructor%outdir  =  cpar%outdir
    constructor%nbins   =  nbins
    
    allocate(constructor%adc_in(constructor%nbins), constructor%adc_out(constructor%nbins))
    allocate(constructor%rms_bins(constructor%nbins), constructor%v_bins(constructor%nbins))
    allocate(constructor%nval(constructor%nbins), constructor%vbin_edges(constructor%nbins+1))

    constructor%adc_in(:)     = 0.0
    constructor%adc_out(:)    = 0.0
    constructor%vbin_edges(:) = 0.0
    constructor%v_bins(:)     = 0.0
    constructor%rms_bins(:)   = 0.0
    constructor%nval(:)       = 0

    ! Initialize v_min and v_max on obscenely wrong numbers
    constructor%v_max = 0.0
    constructor%v_min = 100000.0

  end function constructor

  subroutine load_dpc(self, table_in, table_out)
    !=========================================================================
    ! Load the dpc correction tables into the diode adc object
    ! 
    ! Inputs:
    !
    ! self    : comm_adc object
    !    
    ! dpc_in  : float array
    !           dpc voltage in grid
    ! dpc_out : float array
    !           dpc voltage out grid

    implicit none
    class(comm_adc),        intent(inout) :: self
    real(dp), dimension(:), intent(in)    :: table_in
    real(dp), dimension(:), intent(in)    :: table_out

    integer(i4b)                          :: leng, ierr

    leng = size(table_in)

    allocate(self%dpc_in(leng))
    allocate(self%dpc_out(leng))

    self%dpc_in(:)  = real(table_in(:))
    self%dpc_out(:) = real(table_out(:))

    ! call mpi_bcast(self%dpc_in,  leng, MPI_REAL, 0, self%comm, ierr)
    ! call mpi_bcast(self%dpc_out, leng, MPI_REAL, 0, self%comm, ierr)
    
  end subroutine load_dpc

  subroutine adc_correct(self, tod_in, tod_out, mode)
    !=========================================================================
    ! Adc corrects a timestream 
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! tod_in : float array
    !    The tod that is to be corrected
    ! 
    ! Outputs : 
    !
    ! tod_out : float array
    !    The adc corrected version of tod_in
    
    implicit none
    class(comm_adc),                 intent(inout) :: self
    real(sp), dimension(:),          intent(in)    :: tod_in
    real(sp), dimension(:),          intent(out)   :: tod_out
    character(len=*),                intent(in)    :: mode
    integer(i4b)                                   :: i, leng

    !real(dp), dimension(:), allocatable            :: dbl_in, dbl_out
    !real(dp), dimension(:), allocatable            :: in_buff, out_buff
    
    leng = size(tod_in)

    ! allocate(self%sadc(leng))
    allocate(dbl_in(self%nbins),dbl_out(self%nbins))

    allocate(in_buff(leng),out_buff(leng))

    in_buff  = dble(tod_in)
    out_buff = 0.d0

    if (trim(mode) == 'gauss') then
       dbl_in   = dble(self%adc_in)
       dbl_out  = dble(self%adc_out)
    else if (trim(mode) == 'dpc') then
       dbl_in   = dble(self%dpc_in)
       dbl_out  = dble(self%dpc_out)
    else if (trim(mode) == 'none') then
       tod_out = tod_in
       return
    end if
    call spline(self%sadc, dbl_in, dbl_out, regular=.true.)
    do i = 1, leng
       out_buff(i) = splint(self%sadc,in_buff(i))
    end do

    !tod_out = real(out_buff)

    ! If adc_correct_type == 'dpc' then
    ! tod_out = tod_in

    deallocate(dbl_in,dbl_out)
    deallocate(in_buff,out_buff)
    
  end subroutine adc_correct

  subroutine construct_voltage_bins(self)

    implicit none
    class(comm_adc), intent(inout) :: self
    integer(i4b)                   :: i, ierr

    if (self%myid == 0) then

       ! Declare bin edges
       do i = 1, self%nbins+1
          self%vbin_edges(i) = (self%v_max-self%v_min)*(i-1)/self%nbins + self%v_min
       end do
       ! Declare bins
       do i = 1, self%nbins
          self%v_bins(i) = (self%vbin_edges(i) + self%vbin_edges(i+1))/2.0
       end do
    end if
    
    call mpi_bcast(self%vbin_edges,self%nbins, MPI_REAL, 0, self%comm, ierr)
    call mpi_bcast(self%v_bins,    self%nbins, MPI_REAL, 0, self%comm, ierr)

  end subroutine construct_voltage_bins
  
  subroutine find_horn_min_max(self,tod_in,flag,flag0)

    implicit none
    class(comm_adc),                   intent(inout) :: self
    real(sp),     dimension(:),        intent(in)    :: tod_in
    integer(i4b), dimension(:),        intent(in)    :: flag
    integer(i4b),                      intent(in)    :: flag0

    integer(i4b)                                     :: i, ierr, leng

    leng = size(tod_in)

    do i = 1, leng
       if (iand(flag(i),flag0) .ne. 0) cycle 
       if (tod_in(i) < self%v_min) self%v_min = tod_in(i)
       if (tod_in(i) > self%v_max) self%v_max = tod_in(i)
    end do

  end subroutine find_horn_min_max

  subroutine bin_scan_rms(self,tod_in,flag,flag0)
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
    class(comm_adc),                   intent(inout) :: self
    real(sp),     dimension(:),        intent(in)    :: tod_in
    integer(i4b), dimension(:),        intent(in)    :: flag
    integer(i4b),                      intent(in)    :: flag0

    real(sp)                                     :: sum
    real(sp), allocatable,          dimension(:) :: rt, rms, tod_trim
    integer(i4b)                                 :: leng, i, j, j_min, j_max
    
    leng = size(tod_in)
    
    allocate(rt(leng-1))
    allocate(rms(leng-1-self%window))
    allocate(tod_trim(leng-1-self%window))
    
    ! Trim the data according to account for windows - this is used for pointing to the correct bins 
    ! This ensures that tod_trim(i) corresponds to rms(i)
    tod_trim = tod_in(int(self%window/2):leng-1-int(self%window/2))
    
    ! Compute pair-wise difference
    do i = 1, leng-1
       rt(i) = tod_in(i+1)-tod_in(i)
    end do
       
    ! Compute the rms within a window around each rt sample (excluding the ends)
    do i = int(self%window/2)+1, leng-1-int(self%window/2)
       sum = 0.d0
       j_min = max(i-int(self%window/2),1)
       j_max = min(i+int(self%window/2), leng-1)
       do j = j_min, j_max
          if (iand(flag(j),flag0) .ne. 0 .or. iand(flag(j+1),flag0) .ne. 0) cycle 
          sum = sum + rt(j)**2
       end do
       rms(i-int(self%window/2)) = sqrt(sum/(j_max-j_min+1))
    end do
    
    ! Bin the rms values as a function of input voltage, and take the mean
    do i = 1, leng-1-self%window
       do j = 1, self%nbins
          if (tod_trim(i) .ge. self%vbin_edges(j) .and. tod_trim(i) .lt. self%vbin_edges(j+1)) then
             self%nval(j)     = self%nval(j) + 1
             self%rms_bins(j) = self%rms_bins(j) + rms(i)
          end if
       end do
    end do

    deallocate(rt)
    deallocate(rms)
    deallocate(tod_trim)
    
  end subroutine bin_scan_rms

  subroutine build_table(self,handle,name)
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
    !=========================================================================
    
    
    implicit none
    class(comm_adc),                 intent(inout) :: self
    type(planck_rng),                intent(inout) :: handle
    character(len=50),               intent(in)    :: name
    real(sp),     dimension(:),      allocatable   :: tod_trim, idrf, rirf, flatrirf
    real(sp),     dimension(:),      allocatable   :: linrms, flatrms, model
    real(sp),     dimension(:),      allocatable   :: x1, x2, y1, y2
    integer(i4b), dimension(:),      allocatable   :: binmask
    integer(i4b), dimension(:),      allocatable   :: dummymask
    integer(i4b), dimension(:),      allocatable   :: v_dips
    integer(i4b)                                   :: i, j, leng
    integer(i4b)                                   :: ierr, trims
    integer(i4b)                                   :: dip1, v_off, diprange
    real(sp)                                       :: sum, slope, offset
    real(sp)                                       :: middle_mean, middle_std
    
    ! Combine together all of the bins determined from chunk adding

    call mpi_allreduce(mpi_in_place,self%rms_bins,self%nbins,MPI_REAL, MPI_SUM, self%comm, ierr)
    call mpi_allreduce(mpi_in_place,self%nval,self%nbins,MPI_INTEGER, MPI_SUM, self%comm, ierr)
    
    ! The rest should be light enough to do on a single core
    if (self%myid == 0) then
       allocate(binmask(self%nbins))
       allocate(dummymask(self%nbins))

       dummymask(:) = 1
       
       ! Initialize the middle mean and std
       middle_mean = 0.0
       middle_std  = 0.0

       ! Mask out bins that have no entries - otherwise return mean rms for each bin
       i = 0
       do j = 1, self%nbins
          if(self%nval(j) == 0) then
             binmask(j) = 0
             cycle
          end if
          if(self%rms_bins(j) == 0) then
             binmask(j) = 0
             cycle
          end if
          binmask(j)    = 1
          self%rms_bins(j) = self%rms_bins(j)/self%nval(j)
          
          ! Calculate the mean and std of the middle 200 bins
          if (j > 150 .and. j < 350) then
             if (binmask(j) == 0) cycle
             i = i + 1
             middle_mean = middle_mean + self%rms_bins(j)
          end if
       end do

       !! This is some seg-faulty territory
       if (i == 0) then
          middle_mean = 0.d0
       else
          middle_mean = middle_mean/i
       end if
       
       i = 0
       do j = 150, 350
          if (binmask(j) == 0) cycle
          i = i + 1
          middle_std = middle_std + (middle_mean-self%rms_bins(j))**2
       end do
       if (i == 0) then
          middle_std = 0.d0
       else
          middle_std = sqrt(middle_std/i)
       end if

       ! Mask out large deviations
       ! let large deviations below the mean to stay so we can find the dips
       do j = 1, self%nbins
          if (self%rms_bins(j) < 0.5*middle_mean) binmask(j) = 0
          if (self%rms_bins(j) > middle_mean+2.0*middle_std) binmask(j) = 0
       end do

       ! How many edge bins do we remove for our fit?
       trims = 10

       ! Trim 10 off the "bottom"
       i = 0
       do j = 1, self%nbins
          if (binmask(j) /= 0) then
             binmask(j) = 0
             i = i + 1
          end if
          if (i == trims) exit
       end do

       ! Trim 10 off the "top"
       i = 0
       do j = self%nbins, 1, -1
          if (binmask(j) /= 0) then
             binmask(j) = 0
             i = i + 1
          end if
          if (i == trims) exit
       end do

       ! Remove the linear term from V vs RMS before fitting for the dips
       call return_linreg(self%v_bins, self%rms_bins, binmask, slope, offset,trim=.true.)

       ! Allocate and intialize everything
       allocate(linrms(self%nbins))
       allocate(flatrms(self%nbins))
       allocate(idrf(self%nbins),rirf(self%nbins),model(self%nbins))
       !allocate(rirf(self%nbins),model(self%nbins))
       allocate(flatrirf(self%nbins))

       rirf(:)     = 0.0
       idrf(:)     = 0.0

       ! Remove toe linear time 
       linrms  = slope*self%v_bins + offset
       flatrms = self%rms_bins - linrms

       ! After we remove the linear term, let's fit Gaussians to return the inverse differential resposne function

       ! First step is to identify dips
       dip1  = 0
       v_off = 0
       diprange = 10

       v_dips = return_dips(self%v_bins, flatrms, binmask, diprange)

       open(44, file=trim(self%outdir)//'/adc_binned_rms_'//trim(name)//'_flat.dat')
       open(45, file=trim(self%outdir)//'/adc_linear_term_'//trim(name)//'.dat')
       ! open(46, file=trim(self%outdir)//'/adc_response_function_'//trim(name)//'.dat')
       ! open(47, file=trim(self%outdir)//'/adc_rirf_'//trim(name)//'.dat')
       ! open(48, file=trim(self%outdir)//'/adc_idrf_'//trim(name)//'.dat')
       ! open(49, file=trim(self%outdir)//'/adc_model_'//trim(name)//'.dat') 
       open(50, file=trim(self%outdir)//'/adc_binned_rms_'//trim(name)//'.dat')
       open(51, file=trim(self%outdir)//'/adc_voltage_bins_'//trim(name)//'.dat')
       ! open(53, file=trim(self%outdir)//'/adc_out_'//trim(name)//'.dat') 
       open(54, file=trim(self%outdir)//'/adc_binmask_'//trim(name)//'.dat') 
       do i = 1, self%nbins
          write(44, fmt='(f30.8)') flatrms(i)
          write(45, fmt='(f16.8)') linrms(i)
          ! write(46, fmt='(f16.8)') flatrirf(i)
          ! write(47, fmt='(f16.8)') rirf(i)
          ! write(48,fmt='(f16.8)')  idrf(i)
          ! write(49,fmt='(f16.8)')  model(i)
          write(50, fmt='(f16.8)') self%rms_bins(i)
          write(51, fmt='(f16.8)') self%v_bins(i)
          ! write(53, fmt='(f16.8)') self%adc_out(i)
          write(54, fmt='(i1)')    binmask(i)
       end do
       close(44)
       close(45)
       ! close(46)
       ! close(47)
       ! close(48)
       ! close(49)
       close(50)
       close(51)
       ! close(53)
       close(54)

       ! Return the Inverse Differential Response Function (idrf)
       idrf = return_gaussian_idrf_dips(self%v_bins, flatrms, binmask, v_dips, handle,name)

       ! Also create a composition of our model to compare with the binned RMS
       model = linrms - idrf

       ! open(44, file=trim(self%outdir)//'/adc_binned_rms_'//trim(name)//'_flat.dat')
       ! open(45, file=trim(self%outdir)//'/adc_linear_term_'//trim(name)//'.dat')
       ! open(46, file=trim(self%outdir)//'/adc_response_function_'//trim(name)//'.dat')
       ! open(47, file=trim(self%outdir)//'/adc_rirf_'//trim(name)//'.dat')
       open(48, file=trim(self%outdir)//'/adc_idrf_'//trim(name)//'.dat')
       open(49, file=trim(self%outdir)//'/adc_model_'//trim(name)//'.dat') 
       ! open(50, file=trim(self%outdir)//'/adc_binned_rms_'//trim(name)//'.dat')
       ! open(51, file=trim(self%outdir)//'/adc_voltage_bins_'//trim(name)//'.dat')
       ! open(53, file=trim(self%outdir)//'/adc_out_'//trim(name)//'.dat') 
       ! open(54, file=trim(self%outdir)//'/adc_binmask_'//trim(name)//'.dat') 
       do i = 1, self%nbins
          ! write(44, fmt='(f30.8)') flatrms(i)
          ! write(45, fmt='(f16.8)') linrms(i)
          ! write(46, fmt='(f16.8)') flatrirf(i)
          ! write(47, fmt='(f16.8)') rirf(i)
          write(48,fmt='(f16.8)')  idrf(i)
          write(49,fmt='(f16.8)')  model(i)
          ! write(50, fmt='(f16.8)') self%rms_bins(i)
          ! write(51, fmt='(f16.8)') self%v_bins(i)
          ! write(53, fmt='(f16.8)') self%adc_out(i)
          ! write(54, fmt='(i1)')    binmask(i)
       end do
       ! close(44)
       ! close(45)
       ! close(46)
       ! close(47)
       close(48)
       close(49)
       ! close(50)
       ! close(51)
       ! close(53)
       ! close(54)


       ! Integrate up that idrf to receive the Reconstructed Inverse Response Function (rirf)
       ! This function is what we actually use as our correction function
       do i = 1, self%nbins
          rirf(i) = sum(idrf(:i))
       end do

       ! Can't forget to remove the linear part of rirf so as to not be degenerate in gain       
       if (dip1 /= 0) then 
          call return_linreg(self%v_bins, rirf, dummymask, slope, offset)
          do i = 1, self%nbins
             flatrirf(i) = rirf(i) - slope*self%v_bins(i) - offset
          end do
       end if

       ! Now finally declare the actual adc_in and adc_out tables
       self%adc_in  = self%v_bins
       self%adc_out = self%v_bins
       do i = 1, self%nbins
          self%adc_out(i) = self%adc_out(i) + flatrirf(i) - slope*self%v_bins(i) - offset
       end do
       
       ! Write to file binned rms, voltages, and response function to files

       ! open(44, file=trim(self%outdir)//'/adc_binned_rms_'//trim(name)//'_flat.dat')
       ! open(45, file=trim(self%outdir)//'/adc_linear_term_'//trim(name)//'.dat')
       open(46, file=trim(self%outdir)//'/adc_response_function_'//trim(name)//'.dat')
       open(47, file=trim(self%outdir)//'/adc_rirf_'//trim(name)//'.dat')
       ! open(48, file=trim(self%outdir)//'/adc_idrf_'//trim(name)//'.dat')
       ! open(49, file=trim(self%outdir)//'/adc_model_'//trim(name)//'.dat') 
       ! open(50, file=trim(self%outdir)//'/adc_binned_rms_'//trim(name)//'.dat')
       ! open(51, file=trim(self%outdir)//'/adc_voltage_bins_'//trim(name)//'.dat')
       open(53, file=trim(self%outdir)//'/adc_out_'//trim(name)//'.dat') 
       ! open(54, file=trim(self%outdir)//'/adc_binmask_'//trim(name)//'.dat') 
       do i = 1, self%nbins
          ! write(44, fmt='(f30.8)') flatrms(i)
          ! write(45, fmt='(f16.8)') linrms(i)
          write(46, fmt='(f16.8)') flatrirf(i)
          write(47, fmt='(f16.8)') rirf(i)
          ! write(48,fmt='(f16.8)')  idrf(i)
          ! write(49,fmt='(f16.8)')  model(i)
          ! write(50, fmt='(f16.8)') self%rms_bins(i)
          ! write(51, fmt='(f16.8)') self%v_bins(i)
          write(53, fmt='(f16.8)') self%adc_out(i)
          ! write(54, fmt='(i1)')    binmask(i)
       end do
       ! close(44)
       ! close(45)
       close(46)
       close(47)
       ! close(48)
       ! close(49)
       ! close(50)
       ! close(51)
       close(53)
       ! close(54)

       deallocate(linrms)
       deallocate(flatrms)
       deallocate(idrf,rirf,model)
       deallocate(flatrirf)
    end if

    ! mpi_bcast the tables to all other cores
    call mpi_bcast(self%adc_in,  self%nbins, MPI_REAL, 0, self%comm, ierr) 
    call mpi_bcast(self%adc_out, self%nbins, MPI_REAL, 0, self%comm, ierr) 
    
  end subroutine build_table
    
  subroutine return_linreg(x,y,mask,slope,offset,trim)
  ! subroutine return_linreg(x,y,slope,offset)
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
    real(sp), dimension(:,:),    allocatable   :: adata, bdata
    real(sp), dimension(:),      allocatable   :: work
    real(sp), dimension(:),      allocatable   :: x2, y2
    integer(i4b)                               :: i, leng, count
    real(sp)                                   :: y_mean, y_var, y_std
    integer(i4b)                               :: info
    
    leng = size(x)

    y_mean = 0.0
    y_var  = 0.0
    y_std  = 0.0

    count = 0
    do i = 1, leng
       if (mask(i) == 0) cycle
       count = count + 1
       y_mean = y_mean + y(i)
    end do
    if (count == 0) then
       slope = 0.d0
       offset = 0.d0
       return
    end if
    y_mean = y_mean/count
    
    do i = 1, leng
       if (mask(i) == 0) cycle
       y_var  = y_var + (y(i)-y_mean)**2
    end do
    
    y_var = y_var/(count-1)
    
    y_std  = sqrt(y_var)

    allocate(x2(leng),y2(leng))

    x2(:) = 0.0
    y2(:) = 0.0

    ! Mask out bins that have no entries
    count = 0
    do i = 1, leng
       ! Mask out outliers (if (y < mean-std .or. y > mean + std))
       if (present(trim)) then
          if (y(i) < y_mean-y_std) cycle
          if (y(i) > y_mean+y_std) cycle
       end if
       if (mask(i) == 0) cycle
       count     = count + 1
       x2(count) = x(i)
       y2(count) = y(i) 
    end do

    allocate(adata(count,2), bdata(count,1), work(2*count))
    adata(:,:) = 0.0
    bdata(:,:) = 0.0
    work(:)    = 0.0
    
    adata(:,1) = 1.0
    adata(:,2) = x2(1:count)
    bdata(:,1) = y2(1:count)
    
    call sgels('N', count, 2, 1, adata, count, bdata, count, work, 2*count, info)
    
    slope  = bdata(2,1)
    offset = bdata(1,1)
    
    deallocate(adata)
    deallocate(bdata)
    deallocate(work)
    
  end subroutine return_linreg

  function return_dips(x,y,mask,diprange)
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
    ! vrange  : integer
    !           how many x axis indices correspond to 1mV
    !
    ! Outputs:
    !
    ! dip1    : float
    !           voltage value - location of the "first" dip
    !
    ! v_off   : float
    !           voltage value - estimated distance between dips
    ! ====================================================================
    
    implicit none
    
    real(sp),     dimension(:), intent(in)    :: x, y
    integer(i4b), dimension(:), intent(in)    :: mask
    integer(i4b),               intent(in)    :: diprange
    logical(lgt), dimension(:), allocatable   :: truths
    integer(i4b), dimension(:), allocatable   :: dips
    integer(i4b), dimension(:), allocatable   :: return_dips
    integer(i4b)                              :: leng, i, j, count
    integer(i4b)                              :: ndips
    real(sp)                                  :: y_mean, y_var, y_std
    
    leng = size(x)
    
    allocate(truths(diprange))
    allocate(dips(leng))

    dips(:)   = 0
    ndips     = 0
    truths(:) = .false.
    y_mean    = 0.0
    y_var     = 0.0

    ! Determine mean and standard deviation of the input y-array
    count = 0
    do i = 1, leng
       if (mask(i) == 0) cycle
       count  = count + 1
       y_mean = y_mean + y(i)
    end do
    if (count == 0) then
       dip1 = 0.d0
       v_off = 0.d0
       deallocate(truths, dips)
       return
    end if
    y_mean = y_mean/count
    
    count = 0
    do i = 1, leng
       if (mask(i) == 0) cycle
       count = count + 1
       y_var  = y_var + (y(i)-y_mean)**2
    end do
    
    y_var = y_var/(count-1)
    
    y_std  = sqrt(y_var)

    ! Since the linear portion has been removed, mean should be near zero,
    ! so dips are identified first by finding y-values where y < -1.0*y_std
    
    do i = 1, leng
       if (mask(i) == 0) cycle
       ! Only consider variations
       if (y(i) < y_mean-1.0*y_std .and. y(i-1) < y_mean-1.0*y_std .and. y(i+1) < y_mean-1.0*y_std) then
          truths(:) = .false.
          ! search local range
          do j = 1, diprange
             if (i+j == leng) then
                exit
             end if
             ! if lower than your neighbors, share the good news!
             if (y(i) < y(i-j) .and. y(i) < y(i+j)) then
                truths(j) = .true.
             else
                truths(j) = .false.
             end if
          end do
          ! If lower than all your neighbors
          if (all(truths)) then
             ! append a dip location
             ndips       = ndips + 1
             dips(ndips) = i
          end if
       end if
    end do
    
    return_dips = dips(1:ndips)

    ! dip1 = dips(1)

    ! if (ndips .ne. 1) then
    !    do i = 1, ndips-1
    !       v_off = v_off + dips(i+1) - dips(i)
    !    end do
    !    v_off = v_off
    !    write(*,fmt='(a,i2,a,i4)') 'ndips = ', ndips, ', v_off = ', v_off
    ! else
    !    write(*,*) 'ndips = 1, v_off = 0'
    !    v_off = 0
    ! end if
    deallocate(truths)
    deallocate(dips)

  end function return_dips


  ! Going to change this up to be a bit more modular for the inversion sampler
  ! function return_gaussian(x, mu, sigma, amp) result(y)
  function return_gaussian(x, pars) result(y)
    ! ====================================================================
    ! Super simple function which returns a guassian function given the parameters
    !
    ! Inputs: 
    !
    ! x     : float array
    !         the x-array over which the gaussian is evaluated
    ! mu    : float
    !         mean (center) of the gaussian
    ! sigma : float
    !         width of the guassian
    ! amp   : float
    !         amplitude of the gaussian peak
    !
    ! Outputs:
    !
    ! y     : float array
    !         returned gaussian function
    ! ====================================================================
    
    implicit none
    
    real(sp), dimension(:), intent(in)         :: x
    real(sp), dimension(3), intent(in)         :: pars
    real(sp), allocatable, dimension(:)        :: y
    real(sp)                                   :: mu, sigma, amp
    integer(i4b)                               :: leng, i

    leng = size(x)
    
    allocate(y(leng))

    mu    = pars(1)
    sigma = pars(2)
    amp   = pars(3)

    do i = 1, leng
       y(i) = amp * exp(-((x(i)-mu)/sigma)**2)
    end do
    
  end function return_gaussian

  function return_gaussian_idrf_dips(x, y, mask, dips, handle,name) result(idrf)
    !=========================================================================
    ! Fits a single gaussian function to all recognized dips in the 
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
    
    real(sp), dimension(:),     intent(in)    :: x, y
    integer(i4b), dimension(:), intent(in)    :: mask
    integer(i4b), dimension(:), intent(in)    :: dips
    type(planck_rng),           intent(inout) :: handle
    real(sp), allocatable, dimension(:)       :: idrf, newy, model
    real(sp)                                  :: sigma, amp, mean, fwhm
    real(sp)                                  :: sigma_est, amp_est, mean_est
    integer(i4b)                              :: leng, i, j, ndips
    integer(i4b)                              :: fit_range

    character(len=2)                          :: dip_str
    character(len=50),          intent(in)    :: name

    real(sp), dimension(:), allocatable       :: x_tmp, y_tmp, ret_mod

    integer(i4b)                              :: currpar, first, last, ngibbs, k

    real(dp), dimension(2)                    :: P_par
    real(dp), dimension(3)                    :: x_par
    real(sp), dimension(3)                    :: pars
    real(sp), dimension(3)                    :: par_est
    
    ! declare the goodies
    fit_range = 30
    ngibbs    = 10
    leng      = size(x)
    
    ! allocate all relevant arrays
    allocate(idrf(leng))
    allocate(newy(leng))
    allocate(model(leng))

    ! init cumulative arrays
    idrf(:)  = 0.0
    newy(:)  = 0.0

    ndips = size(dips)

    do j = 1, ndips
       sigma        = 0.0
       mean         = 0.0
       amp          = 0.0

       ! Define first and last for indices - range to fit Gaussian to dip
       first = dips(j)!dip1+(j-1)*v_off
       last  = dips(j)!dip1+(j-1)*v_off
       if (dips(j) - fit_range < 1) then
          first = 1
       else if (dips(j) + fit_range > leng) then
          last = leng
       else 
          do i = dips(j) - fit_range, dips(j) + fit_range
             if (mask(i) == 0) cycle
             first = min(first,i)
             last  = max(last,i)
          end do
       end if

       ! allocate temporary array for voltage bins and rms bins
       allocate(x_tmp(last-first),y_tmp(last-first))

       allocate(ret_mod(last-first))
       ! Flip the dip!

       do i = dips(j)-fit_range, dips(j)+fit_range
          ! skip over incides that are masked or outside the index range
          if (i < 1) cycle
          if (i > leng) cycle
          if (mask(i) == 0) cycle
          newy(i)  = -1.0*y(i)
       end do

       x_tmp(:) = x(first:last)
       y_tmp(:) = newy(first:last)

       write(dip_str,'(i0.2)') j
       
       ! Define estimates to the Gaussian parameters
       do i = dips(j), dips(j)+fit_range
          if (mask(i) == 0) cycle
          if (newy(i) < newy(dips(j))/2.0) then
             fwhm = x(i)-x(dips(j))
             exit
          end if
       end do

       sigma   = max(2.0*(fwhm/2.355), 0.00001)
       mean    = x(dips(j))
       amp     = maxval(y_tmp)

       pars(1) = mean
       pars(2) = sigma
       pars(3) = amp

       ! par_est = pars

       ! write(*,*) 'initial estimates:'
       ! write(*,*) 'mean  = ', pars(1) 
       ! write(*,*) 'sigma = ', pars(2) 
       ! write(*,*) 'amp   = ', pars(3) 

       ! ! With our estimates, let's sample for the parameters
       ! do k = 1, 1!ngibbs
       !    do i = 1, 1!3
       !       currpar  = i
       !       x_par(:) = 0.0
             
       !       ! define parameter uniform prior ranges
       !       if (i == 1) then
       !          ! Ensure the mean value is within the dip range
       !          P_par(1) = minval(x_tmp)
       !          P_par(2) = maxval(x_tmp)
       !       else if (i == 2) then
       !          P_par(1) = 0.75*par_est(i)
       !          P_par(2) = 1.25*par_est(i)
       !       else if (i == 3) then
       !          ! Ensure amplitude is always greater than 0
       !          P_par(1) = 0.0
       !          P_par(2) = maxval(y_tmp)*1.5
       !       end if
             
       !       ! Grid out parameter estimates
       !       x_par(1) = max(par_est(i) - 0.5 * abs(par_est(i)), P_par(1))
       !       x_par(3) = min(par_est(i) + 0.5 * abs(par_est(i)), P_par(2))
       !       x_par(3) = max(x_par(3), x_par(1)+1.d-3*(P_par(2)-P_par(1)))
       !       x_par(2) = 0.5 * (x_par(1) + x_par(3))
             
       !       pars(i) = real(sample_InvSamp(handle, x_par, lnL_dip_n, P_par))

       !    end do
       ! end do

       ! write(*,*) 'from the inversion sampler:'
       ! write(*,*) 'mean  = ', pars(1) 
       ! write(*,*) 'sigma = ', pars(2) 
       ! write(*,*) 'amp   = ', pars(3) 

       ! For overwriting the inversion sampling pars
       ! pars(1) = x(dips(j))
       ! pars(2) = max(2.0*(fwhm/2.355), 0.00001)
       ! pars(3) = maxval(y_tmp)

       ! write(*,*) 'Dip ', j
       ! write(*,*) 'mean  = ', pars(1) 
       ! write(*,*) 'sigma = ', pars(2) 
       ! write(*,*) 'amp   = ', pars(3) 
       ! write(*,*) '-----------------'
       ! write(*,*) ''

       ! ret_mod = return_gaussian(x_tmp,pars)
       ! open(60,file='dip_xs_'//dip_str//'_'//trim(name)//'.dat')
       ! open(61,file='dip_ys_'//dip_str//'_'//trim(name)//'.dat')
       ! open(62,file='dip_model_'//dip_str//'_'//trim(name)//'.dat')
       ! do i = 1, last-first
       !    write(60,*) x_tmp(i)
       !    write(61,*) y_tmp(i)
       !    write(62,*) ret_mod(i)
       ! end do
       ! close(60)
       ! close(61)
       ! close(62)

       deallocate(x_tmp,y_tmp)
       ! deallocate(ret_mod)

       if (amp < 0.0) cycle
       if (sigma < 0.0) cycle

       model = return_gaussian(x,pars)
       idrf  = idrf + model
       
       deallocate(model)
    end do    

  contains
    
    ! A function to evaluate the log-likelihood
    function lnL_dip_n(samp)
      use healpix_types
      implicit none
      real(dp), intent(in) :: samp
      real(dp)             :: lnL_dip_n
      real(sp)             :: tmp

      real(sp), dimension(last-first) :: gauss

      integer(i4b)         :: i

      ! Save old parameter value
      tmp = pars(currpar)
      pars(currpar) = real(samp)

      lnL_dip_n = 0.0

      ! Return the gaussian given the parameters
      gauss = return_gaussian(x_tmp,pars)

      lnl_dip_n = -0.5d0 * dble(sum((y_tmp-gauss)**2))

      ! do i = 1, last-first
      !    lnL_dip_n = lnL_dip_n -0.5d0*dble((y_tmp(i)-gauss(i))**2.0)
      ! end do

      write(*,*) 'samp, lnl_dip_n'
      write(*,*) samp, lnl_dip_n

      ! Put that old value back
      pars(currpar) = tmp

    end function lnL_dip_n
  end function return_gaussian_idrf_dips
  
end module comm_tod_adc_mod
