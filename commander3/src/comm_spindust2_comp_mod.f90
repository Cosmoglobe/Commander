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
module comm_spindust2_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_spindust2_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_spindust2_comp
     real(dp)          :: nu_p0, nu_min_SED, nu_max_SED
     type(spline_type) :: SED_spline
   contains
     procedure :: S    => evalSED_spindust2
  end type comm_spindust2_comp

  interface comm_spindust2_comp
     procedure constructor_spindust2
  end interface comm_spindust2_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_spindust2(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_spindust2_comp), pointer   :: c

    integer(i4b) :: ind(1)
    real(dp), allocatable, dimension(:,:) :: SED
    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 

    ! General parameters
    allocate(c)

    c%npar         = 2
    allocate(c%poltype(c%npar))
    do i = 1, c%npar
       c%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call c%initLmaxSpecind(cpar, id, id_abs)

    call c%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters for 2 parameter model
    allocate(c%theta_def(2), c%p_gauss(2,2), c%p_uni(2,2))
    allocate(c%theta_steplen(2,cpar%mcmc_num_samp_groups))
    allocate(c%indlabel(2))
    allocate(c%nu_min_ind(2), c%nu_max_ind(2))
    do i = 1, 2
       c%theta_def(i) = cpar%cs_theta_def(i,id_abs)
       c%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
       c%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
       c%nu_min_ind(i) = cpar%cs_nu_min_beta(id_abs,i)
       c%nu_max_ind(i) = cpar%cs_nu_max_beta(id_abs,i)
    end do
    c%theta_steplen = 0d0
    c%indlabel  = ['nu_p ','alpha']

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, c%nside, c%lmax_ind, &
         & c%nmaps, c%pol)

    allocate(c%theta(c%npar))
    do i = 1, c%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
          c%theta(i)%p => comm_map(info)
          c%theta(i)%p%map = c%theta_def(i)
       else
          ! Read map from FITS file, and convert to alms
          c%theta(i)%p => comm_map(info, trim(cpar%cs_input_ind(i,id_abs)))
       end if

       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (c%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call c%theta(i)%p%YtW_scalar
       end if
    end do

    ! Initialize spectral template !CHANGE??
    call read_spectrum(trim(cpar%cs_SED_template(1,id_abs)), SED)
    ind                    = maxloc(SED(:,2))
    c%nu_p0      = SED(ind(1),1)
    c%nu_min_SED = minval(SED(:,1))
    c%nu_max_SED = maxval(SED(:,1))
    SED                    = log(SED)
    call spline(c%SED_spline, SED(:,1), SED(:,2))
    deallocate(SED)

    ! Precompute mixmat integrator for each band
    allocate(c%F_int(3,numband,0:c%ndet))
    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (c%nu_ref(k) == c%nu_ref(k-1)) then
                   c%F_int(k,i,j)%p => c%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             c%F_int(k,i,j)%p => comm_F_int_2D(c, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    call c%initPixregSampling(cpar, id, id_abs)
    ! Init alm 
    if (c%lmax_ind >= 0) call c%initSpecindProp(cpar, id, id_abs)

    ! Initialize mixing matrix
    call c%updateMixmat

  end function constructor_spindust2

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED_spindust2(self, nu, band, pol, theta)
    implicit none
    class(comm_spindust2_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_spindust2

    real(dp) :: scale, nu_p, alpha

    nu_p    = theta(1)
    alpha   = theta(2)
    scale   = self%nu_p0 / (nu_p*1.d9) ! nu_p is in GHz
    
    if (scale*nu < self%nu_min_SED .or. scale*nu > self%nu_max_SED) then
       evalSED_spindust2 = 0.d0
    else
       evalSED_spindust2 = exp(splint(self%SED_spline, log(scale*nu))) / &
               & exp(splint(self%SED_spline, log(scale*self%nu_ref(pol)))) * (self%nu_ref(pol)/nu)**(2.d0-alpha)
    end if

  end function evalSED_spindust2
  
end module comm_spindust2_comp_mod
