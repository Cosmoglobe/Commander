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
module comm_powlaw_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_powlaw_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_powlaw_comp
   contains
     procedure :: S    => evalSED_powlaw
  end type comm_powlaw_comp

  interface comm_powlaw_comp
     procedure constructor_powlaw
  end interface comm_powlaw_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_powlaw(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_powlaw_comp), pointer   :: c

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

    c%npar         = 1
    allocate(c%poltype(c%npar))
    do i = 1, c%npar
       c%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call c%initLmaxSpecind(cpar, id, id_abs)

    call c%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(c%theta_def(1), c%p_gauss(2,1), c%p_uni(2,1))
    allocate(c%theta_steplen(1,cpar%mcmc_num_samp_groups))
    allocate(c%indlabel(1))
    allocate(c%nu_min_ind(1), c%nu_max_ind(1))
    c%theta_def(1) = cpar%cs_theta_def(1,id_abs)
    c%p_uni(:,1)   = cpar%cs_p_uni(id_abs,:,1)
    c%p_gauss(:,1) = cpar%cs_p_gauss(id_abs,:,1)
    c%nu_min_ind(1) = cpar%cs_nu_min_beta(id_abs,1)
    c%nu_max_ind(1) = cpar%cs_nu_max_beta(id_abs,1)
    c%apply_jeffreys = cpar%cs_apply_jeffreys(id_abs)
    c%theta_steplen = 0d0
    c%indlabel(1)  = 'beta'

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, c%nside, c%lmax_ind, &
         & c%nmaps, c%pol)

    allocate(c%theta(c%npar))
    do i = 1, c%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
          c%theta(i)%p => comm_map(info)
          c%theta(i)%p%map = c%theta_def(1)
       else
          ! Read map from FITS file
          c%theta(i)%p => comm_map(info, trim(cpar%cs_input_ind(i,id_abs)))
       end if

       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (c%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call c%theta(i)%p%YtW_scalar
       end if
    end do

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
             c%F_int(k,i,j)%p => comm_F_int_1D(c, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    call c%initPixregSampling(cpar, id, id_abs)
    ! Init alm 
    if (c%lmax_ind >= 0) call c%initSpecindProp(cpar, id, id_abs)

    ! Initialize mixing matrix
    call c%updateMixmat

  end function constructor_powlaw

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED_powlaw(self, nu, band, pol, theta)
    implicit none
    class(comm_powlaw_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_powlaw

    evalSED_powlaw = (nu/self%nu_ref(pol))**theta(1)

  end function evalSED_powlaw
  
end module comm_powlaw_comp_mod
