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


!================================================================================
! This is a module which doesn't actually do anything, but is just lined
! out as a template for adding a new diffuse sky component to Commander
!================================================================================
module comm_brand_new_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
  implicit none

  private
  public comm_brand_new_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_brand_new_com
   contains
     procedure :: S    => evalSED
  end type comm_brand_new_com

  interface comm_brand_new_com
     procedure constructor
  end interface comm_brand_new_com

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),           intent(in) :: cpar
    integer(i4b),                intent(in) :: id, id_abs
    class(comm_brand_new_com), pointer   :: constructor

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),     allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop

    real(dp)                     :: par_dp
    character(len=512)           :: temptxt, partxt
    integer(i4b)                 :: smooth_scale, p_min, p_max
    type(comm_mapinfo),  pointer :: info => null()
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 

    ! General parameters
    allocate(constructor)

    ! How many spectral parameters are there for this component?
    constructor%npar         = 1
    
    ! Allocate polarization types for each parameter
    allocate(constructor%poltype(constructor%npar))
    
    ! Read the poltypes from the Commander parameters
    do i = 1, constructor%npar
       constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do

    ! Initialize the LMAX for each spectral index
    call constructor%initLmaxSpecind(cpar, id, id_abs)

    ! Initialize this diffuse component
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(constructor%theta_def(2), constructor%p_gauss(2,2), constructor%p_uni(2,2))
    allocate(constructor%indlabel(2))
    ! What are the frequency ranges which we care about for each parameter
    allocate(constructor%nu_min_ind(2), constructor%nu_max_ind(2))
    ! Point all of the parameters from cpar to the component object
    do i = 1, 2
       constructor%theta_def(i) = cpar%cs_theta_def(i,id_abs)
       constructor%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
       constructor%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
       constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
       constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
    end do
    
    ! What are the index labels for each index?
    constructor%indlabel = ['beta']

    ! Initialize spectral index map by creating a map object
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    ! Allocate your spectral parameters
    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar
       ! Do we read in an initialization map for the spectral parameters?
       ! If not, create a map object for each 
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(1)
       else
          ! Read map from FITS file
          constructor%theta(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(i,id_abs)))
       end if

       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (constructor%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call constructor%theta(i)%p%YtW_scalar
       end if
    end do

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    ! For each map
    do k = 1, 3
       ! For each band
       do i = 1, numband
          ! and each detector
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (constructor%nu_ref(k) == constructor%nu_ref(k-1)) then
                   constructor%F_int(k,i,j)%p => constructor%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             constructor%F_int(k,i,j)%p => comm_F_int_1D(constructor, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    ! Call the pixel region sampler initialization
    call constructor%initPixregSampling(cpar, id, id_abs)

    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_brand_new_com), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    ! Here is where you define the funcitonal form for the evaluation of the
    ! spectral energy distribution given a frequency, and your spectral parameters
    !
    ! Here is an example for the modified blackbody, who's funcitonal form
    ! in Rayleigh-Jeans units is:
    ! 
    ! I_nu = (nu/nu_ref)^beta_d * B_nu(T_d)
    ! 
    ! In Commander form:
    !
    ! beta    = theta(1)
    ! T       = theta(2)
    ! x       = h*nu               / (k_b*T)
    ! x_ref   = h*self%nu_ref(pol) / (k_b*T)
    ! evalSED = (nu/self%nu_ref(pol))**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)
    !

    ! Code your own SED here:
    evalSED = 1.d0

  end function evalSED
  
end module comm_brand_new_com_mod
