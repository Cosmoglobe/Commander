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

  use spline_1D_mod

  implicit none

  private
  public comm_adc, adc_pointer

  type :: comm_adc
    real(dp), dimension(:), allocatable :: adc_in, adc_out
  contains
    procedure :: adc_correct
  end type comm_adc

  interface comm_adc
    procedure constructor
  end interface comm_adc

  type adc_pointer
    class(comm_adc), pointer :: p => null()
  end type adc_pointer

contains

   function constructor(adc_i, adc_o)
     ! ====================================================================
     ! Sets up an adc correction object that maps input and output voltages
     !
     ! Inputs:
     ! 
     ! adc_i : float array
     !   The array of input voltages
     !
     ! adc_o : float array
     !   The array of output voltages
     ! ====================================================================

     implicit none
     real(dp), dimension(:), intent(in)          :: adc_i, adc_o
     class(comm_adc), pointer                    :: constructor    

     allocate(constructor)

     allocate(constructor%adc_in(size(adc_i)), constructor%adc_out(size(adc_o)))

     constructor%adc_in  = adc_i
     constructor%adc_out = adc_o
 
   end function constructor

   subroutine return_linreg(x,y,mask,slope,offset)
     implicit none
     
     real(sp),     dimension(:),  intent(in)    :: x, y
     integer(i4b), dimension(:),  intent(in)    :: mask
     real(sp),                    intent(inout) :: slope, offset
     integer(i4b)                               :: i, len
     real(sp)                                   :: sumx, sumy, sumxy, sumx2

     len = size(x)

     sumx  = 0.d0
     sumy  = 0.d0
     sumxy = 0.d0
     sumx2 = 0.d0

     ! Mask out bins that have no entries
     do i = 1, len
        sumx  = sumx  + x(i)
        sumy  = sumy  + y(i)
        sumxy = sumxy + x(i)*y(i)
        sumx2 = sumx2 + x(i)*x(i)
     end do

     slope  = (sumy*sumx2 - sumx*sumxy)/(len*sumx2-sumx**2)
     offset = (len*(sumxy) - sumx*sumy)/(len*sumx2-sumx**2)

   end subroutine return_linreg

   subroutine solve_lin_eq_via_inv(A,x,b)
     ! ====================================================================
     ! Solve a simple linear set of equations through matrix inversion
     ! For use to find the Gaussian parameters for the dips
     !
     ! Inputs:
     !
     ! A      : float array - 3x3
     !          the input matrix 
     !
     ! b      : float array - size 3
     !          input RHS
     !
     ! Outputs:
     !
     ! x      : float array - size 3
     !          output Gaussian parameters
     !
     ! ====================================================================

     implicit none
     real(sp), dimension(3,3), intent(inout) :: A
     real(sp), dimension(3),   intent(in)    :: b
     real(sp), dimension(3),   intent(inout) :: x
     real(sp), dimension(3,3)                :: a1
     integer(i4b)                            :: i,j,k,l,n
     real(sp)                                :: z

     n = 3
 
     a1(1,1) = 1.0
     a1(2,2) = 1.0
     a1(3,3) = 1.0

     do i = 1, n
        z = a(i,i)
        do j = 1, n
           a(i,j)  = a(i,j)/z
           a1(i,j) = a1(i,j)/z
        end do
        do j = i+1,n
           z = a(j,i)
           do k = 1,n
              a(j,k)  = a(j,k)-z*a(i,k)
              a1(j,k) = a1(j,k)-z*a1(i,k)
           end do
        end do
     end do

     do i = 1, n-1
        do j = i+1,n
           z = a(i,j)
           do l=1, n
              a(i,l)  = a(i,l)-z*a(j,l)
              a1(i,l) = a1(i,l) - z*a1(j,l)
           end do
        end do
     end do

     ! Now multiply A1 and b to get x
     do i=1,n
        do j=1,n
           x(i) = 0
           do k =1,n
              x(i) = x(i) + a1(i,k)*b(k)
           end do
        end do
     end do

   end subroutine solve_lin_eq_via_inv

   subroutine return_v_off(x,y,binmask,vrange,dip1,v_off)
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
     ! binmask : integer array
     !           mask for binned rms array - 0 if y(i) = NaN
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
     integer(i4b), dimension(:), intent(in)    :: binmask
     integer(i4b),               intent(in)    :: vrange
     real(sp),                   intent(inout) :: dip1, v_off
     logical(lgt), dimension(:), allocatable   :: truths
     integer(i4b), dimension(:), allocatable   :: dips
     integer(i4b)                              :: len, i, j
     integer(i4b)                              :: diprange, ndips
     real(dp)                                  :: y_mean, y_var, y_std

     diprange = vrange/2
     len = size(x)

     allocate(truths(diprange))
     allocate(dips(len))

     dip1      = 0.0
     v_off     = 0.0
     dips(:)   = 0
     truths(:) = .false.

     ! Determine mean and standard deviation of the input y-array
     y_mean = sum(y)/len

     do i = 1, len
        y_var  = y_var + (y(i)-y_mean)**2
     end do

     y_var = y_var/(len-1)
     
     y_std  = sqrt(y_var)

     ! Since the linear portion has been removed, mean should be near zero,
     ! so dips are identified first by finding y-values where y < -1.0*y_std

     do i = 1, len
        ! Only consider variations
        if (y(i) < -1.0*y_std) then
           truths(:) = .false.
           ! search local range
           do j = 1, diprange
              if (i+j == len) then
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
              dips(i) = 1
              ndips   = ndips + 1
           end if
        end if
     end do

     dip1 = x(dips(1))
     
     do i = 1, ndips-1
        v_off = v_off + x(dips(i+1)) - x(dips(i))
     end do
     v_off = v_off/ndips

   end subroutine return_v_off

   function return_gaussian(x, mu, sigma, amp) result(y)
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
     real(sp), allocatable, dimension(:)        :: y
     real(sp)                                   :: mu, sigma, amp
     integer(i4b)                               :: len, i

     allocate(y(len))

     do i = 1, len
        y(i) = amp * exp(((x(i)-mu)/sigma)**2)
     end do

   end function return_gaussian

   function return_gaussian_idrf(x, y, binmask, mu_est, v_off) result(idrf)
     !=========================================================================
     ! Fits a single gaussian function to all recognized dips in the 
     ! white noise level
     !
     ! Inputs:
     ! 
     ! self : comm_adc object (probably)
     !
     ! x         : float array
     !             voltage bins from the tod_in binning
     ! y         : float array
     !             rms bins from the rms estimates as a function of voltage
     ! mu_est    : float
     !             estimated first gaussian dip locations - a voltage value 
     !
     ! Outputs: 
     !
     ! idrf      : float array
     !             array of length(x) - the inverse differential response function
     !             which will be integrated
     !=========================================================================

     implicit none

     real(sp), dimension(:),     intent(in)    :: x, y
     real(sp),                   intent(inout) :: mu_est
     real(sp),                   intent(in)    :: v_off
     integer(i4b), dimension(:), intent(in)    :: binmask
     real(sp), allocatable, dimension(:)       :: y_gauss, idrf, sol, rhs 
     real(sp), allocatable, dimension(:)       :: newx, newy, model
     real(sp), allocatable, dimension(:,:)     :: mat
     real(sp)                                  :: sigma, amp, mean
     integer(i4b)                              :: len, ndips, i, vr2 
     integer(i4b)                              :: minus_dips
 
     len = size(x)

     allocate(y_gauss(len),idrf(len))

     ! Set up linear system with a 3x3 matrix to solve for the Gaussian fit parameters
     
     allocate(mat(3,3),sol(3),rhs(3))
     allocate(newx(len),newy(len))
     allocate(model(len))

     ! Define a new x-array by taking the modulus of x, v_off. This way we fit a single
     ! gaussian to all of the data points simultaneously

     do i = 1, len
        newx(i) = mod(x(i),v_off)
        newy(i) = -1.0*newy(i)
     end do

     ! Using Caruana's algorithm for gaussian parameter determination:
     ! https://arxiv.org/pdf/1907.07241.pdf
     !
     ! In other words, solving the 3 dimensional linear system:
     !
     !  [   N      sum(x)   sum(x^2) ] [ a ]   [   sum(ln(y))   ]
     !  [ sum(x)   sum(x^2) sum(x^3) ] [ b ] = [  sum(x*ln(y))  ]
     !  [ sum(x^2) sum(x^3) sum(x^4) ] [ c ]   [ sum(x^2*ln(y)) ]
     !

     do i = 1, len
        if(binmask(i) == 0) cycle
        mat(1,1) = mat(1,1) + 1.0
        mat(1,2) = mat(1,2) + newx(i)
        mat(2,2) = mat(2,2) + newx(i)*newx(i)
        mat(2,3) = mat(2,3) + newx(i)**3.0
        mat(3,3) = mat(3,3) + newx(i)**4.0
        
        rhs(1)   = rhs(1) + log(newy(i))
        rhs(2)   = rhs(2) + newx(i)*log(newy(i))
        rhs(3)   = rhs(3) + newx(i)*newx(i)*log(newy(i))
     end do
     
     mat(2,1) = mat(1,2)
     mat(1,3) = mat(2,2)
     mat(3,1) = mat(2,2)
     mat(3,2) = mat(2,3)

     call solve_lin_eq_via_inv(mat,sol,rhs)

     ! Take the output array a, b, c to determine the parameters as follows:
     amp    = exp(sol(1) - (sol(2)**2 / 4*sol(3)))
     sigma  = sqrt(-1/(2*sol(3))) 
     mu_est = -sol(2)/(2*sol(3))

     ! I think what would be nice is to just make a new grid from v_min to v_max
     ! and solve for an idrf that spans the whole region

     minus_dips = -1.0*int((mu_est-minval(x))/v_off)

     ndips      = 1 + int((maxval(x)-mu_est)/v_off) + int((mu_est-minval(x))/v_off)

     idrf(:)   = 0.d0

     do i = minus_dips, ndips+minus_dips
        mean   = mu_est + (i-1)*v_off
        model  = return_gaussian(x, mean, sigma, amp)
        idrf   = idrf + model
     end do

   end function return_gaussian_idrf

   subroutine return_diode_response_function(self,tod_in,window,tod_out)
     !=========================================================================
     ! Adc corrects a timestream 
     ! 
     ! Inputs:
     !
     ! self : comm_adc object
     !    Defines the adc correction that should be applied
     ! tod_in : float array
     !    The tod that is to be corrected
     ! window : integer
     !    Window size for the sliding mean used to determine the data RMS
     ! 
     ! Outputs : 
     !
     ! response_function : float array
     !    Reconstructed inverse response function correcting the non-linearity 
     !    identified by the dips in the white-noise level
     !=========================================================================
     

     implicit none
     class(comm_adc),                 intent(inout) :: self
     real(sp),     dimension(:),      intent(in)    :: tod_in
     integer(i4b),                    intent(in)    :: window
     real(sp),     dimension(:),      intent(out)   :: tod_out
     type(spline_type)                              :: sresponse
     real(sp),     dimension(:),      allocatable   :: tod_trim, idrf, rirf, bins
     real(sp),     dimension(:),      allocatable   :: rt, rms, binval, bin_edges
     real(sp),     dimension(:),      allocatable   :: binned_rms, linrms, flatrms
     real(dp),     dimension(:),      allocatable   :: tod_buf, spline_buf, bins_buf, tod_out_buf
     integer(i4b), dimension(:),      allocatable   :: nval, binmask
     integer(i4b)                                   :: i, j, len, vrange
     integer(i4b)                                   :: nbins
     real(sp)                                       :: dip1, sum, v_off, slope, offset

     reaL(sp) :: dummy

     len = size(tod_in)

     allocate(rt(len-1))
     allocate(rms(len-1-window))
     allocate(bin_edges(nbins+1))
     allocate(bins(nbins))
     allocate(binmask(nbins))
     allocate(binval(nbins))
     allocate(nval(nbins))
     allocate(binned_rms(nbins),linrms(nbins),flatrms(nbins))

     ! Declare bin edges
     do i = 1, nbins+1
        bin_edges(i) = (maxval(tod_in)-minval(tod_in))*(i-1)/nbins + minval(tod_in)
     end do

     ! Declare bins
     do i = 1, nbins
        bins(i) = (bin_edges(i) + bin_edges(i+1))/2.0
     end do

     ! Trim the data according to account for windows
     tod_trim = tod_in(int(window/2):len-1-int(window/2))

     ! Compute pair-wise difference
     rt = tod_in(1:)-tod_in(:len-1)

     ! Compute the rms within a window around each rt sample (excluding the ends)
     do i = int(window/2), len-1-int(window/2)
        sum = 0.d0
        do j = i-int(window/2), i+int(window/2)
           sum = sum + rt(j)**2
        end do
        rms(i) = sqrt(sum/window)
     end do

     ! Bin the rms values as a function of input voltage, and take the mean
     do i = 1, len-1-window
        do j = 1, nbins
           if (tod_trim(i) > bin_edges(j) .and. tod_trim(i) < bin_edges(j+1)) then
              nval(j)   = nval(j) + 1
              binval(j) = binval(j) + rms(i)
           end if
        end do
     end do

     ! Mask out bins that have no entries
     do j = 1, nbins
        if(nval(j) == 0) then
           ! write(*,*) "YELL VERY LOUDLY"
           binmask(i) = 0
           cycle
        end if
        binmask(j)    = 1
        binned_rms(j) = binval(j)/nval(j)
     end do

     call return_linreg(bins, binned_rms, binmask, slope, offset)
     
     linrms  = slope*bins + offset
     flatrms = binned_rms - linrms

     ! After we remove the linear term, let's fit Gaussians to return the
     ! inverse differential resposne function

     ! First step here is to estimate the voltage space between the dips
     
     ! Let's look over a window of ~1mV
     ! Define this range in terms of indices
     do i = 1, nbins
        if (bins(i) - bins(1) > 0.001) then
           vrange = i
           exit
        end if
     end do

     call return_v_off(bins, flatrms, binmask, vrange, dip1, v_off)

     ! Return the Inverse Differential Response Function (idrf)
     idrf = return_gaussian_idrf(bins, flatrms, binmask, dip1, v_off)

     ! Integrate up that idrf to receive the Reconstructed Inverse Response Function (rirf)
     ! This function is what we actually use as our correction function
     do i = 1, nbins
        do j = 1, i
           rirf(i) = rirf(i) + idrf(j)
        end do
     end do

     ! Can't forget to remove the linear part of rirf so as to not be degenerate in gain

     call return_linreg(bins, rirf, binmask, slope, offset)

     rirf = rirf - slope*bins

     ! Write to file binned rms, voltages, and response function to files
     ! if (self%myid == 0) then
     !    open(50, file='adc_binned_rms.dat')
     !    open(51, file='adc_voltage_bins.dat')
     !    open(52, file='adc_response_function.dat')
     !    do i = 1, nbins
     !       write(50, fmt='(f16.8)') binned_rms(i)
     !       write(51, fmt='(f16.8)') bins(i)
     !       write(52, fmt='(f16.8)') rirf(i)
     !    end do
     !    close(50)
     !    close(51)
     !    close(52)
     ! end if

     allocate(tod_buf(size(tod_in)), tod_out_buf(size(tod_out)), spline_buf(size(bins)), bins_buf(size(bins)))     

     spline_buf = bins + rirf
     bins_buf = bins

     call spline(sresponse, bins_buf, spline_buf, regular=.true.)

     tod_buf = tod_in

     call splint_simple_multi(sresponse, tod_buf, tod_out_buf)

     tod_out = tod_out_buf

     call free_spline(sresponse)
     deallocate(tod_buf, spline_buf, bins_buf, tod_out_buf)
     
   end subroutine return_diode_response_function

   subroutine adc_correct(self, tod_in, correct_tod)
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
     ! correct_tod : float array
     !    The adc corrected version of tod_in

     implicit none
     class(comm_adc),                 intent(inout) :: self
     real(sp), dimension(:),          intent(in)    :: tod_in
     real(sp), dimension(:),          intent(out)   :: correct_tod
     type(spline_type)                              :: sadc
     integer(i4b)                                   :: i, len

     len = size(self%adc_in)

     !call spline(sadc, self%adc_in, self%adc_out, regular=.true.)

     ! write(*,*) tod_in

     ! stop

     !TODO: figure out the correct algorithm and implement it
     !-------------------------------------------------------

     ! To start we'll the just spline the DPC adc correction table
     ! Must check units!!

     !correct_tod = splint(sadc,tod_in)
     correct_tod = tod_in

   end subroutine adc_correct

end module comm_tod_adc_mod
