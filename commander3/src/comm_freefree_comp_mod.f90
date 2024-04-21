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
module comm_freefree_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_freefree_comp

  !**************************************************
  !      Free-free component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_freefree_comp
   contains
     procedure :: S    => evalSED_freefree
  end type comm_freefree_comp

  interface comm_freefree_comp
     procedure constructor_freefree
  end interface comm_freefree_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_freefree(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_freefree_comp), pointer   :: c

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
    do i = 1, 1
       c%theta_def(i) = cpar%cs_theta_def(i,id_abs)
       c%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
       c%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
       c%nu_min_ind(i) = cpar%cs_nu_min_beta(id_abs,i)
       c%nu_max_ind(i) = cpar%cs_nu_max_beta(id_abs,i)
    end do
    c%theta_steplen = 0d0
    c%indlabel  = ['Te']

    !c%npar         = 1
    !allocate(c%theta_def(1), c%p_gauss(1,1), c%p_uni(1,1))
    !allocate(c%poltype(1), c%indlabel(1))
    !allocate(c%nu_min_ind(1), c%nu_max_ind(1))
    !i = 1
    !c%poltype(i)   = cpar%cs_poltype(i,id_abs)
    !c%theta_def(i) = cpar%cs_theta_def(i,id_abs)
    !c%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
    !c%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
    !c%nu_min_ind(i) = cpar%cs_nu_min_beta(id_abs,i)
    !c%nu_max_ind(i) = cpar%cs_nu_max_beta(id_abs,i)
    
    !c%indlabel  = ['Te']

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, c%nside, c%lmax_ind, &
         & c%nmaps, c%pol)

    allocate(c%theta(c%npar))
    do i = 1, c%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
          c%theta(i)%p => comm_map(info)
          c%theta(i)%p%map = c%theta_def(1)
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

  end function constructor_freefree

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)
  ! where 
  !    beta = theta(1)
  function evalSED_freefree(self, nu, band, pol, theta)
    implicit none
    class(comm_freefree_comp),    intent(in)      :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_freefree
    real(dp)     :: S, S_ref, EM, T_e
    real(dp)     :: g, g_ref, Z_i, tau, tau_ref, EM1, Te

!!$    EM      = theta(1)
!!$    !EM1 = 1.d0 
!!$    Te      = theta(2)
!!$    Z_i     = 1.d0
!!$    g       = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * nu/1.d9          * (Te/1.d4)**(-1.5d0))) + 2.71828d0)
!!$    !g_ref   = log(exp(5.960d0 - sqrt(3.d0)/pi * log(Z_i * self%nu_ref/1.d9 * (Te/1.d4)**(-1.5d0))) + 2.71828d0)
!!$    tau     = 5.468d-2 * Te**(-1.5d0) * (nu/1.d9)**(-2)          * EM * g
!!$    !tau_ref = 5.468d-2 * Te**(-1.5d0) * (self%nu_ref/1.d9)**(-2) * EM * g_ref
!!$
!!$    evalSED_freefree = 1.d6 * Te * (1.d0 - exp(-tau)) !/ (1.d0 - exp(-tau_ref)) 
!!$    
!!$    return
!!$    !write(*,*) "1:", evalSED_freefree


    !EM    = theta(1) ! Not used
    T_e   = theta(1)
    S     = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * nu    /1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    S_ref = log(exp(5.960d0 - sqrt(3.d0)/pi * log(1.d0 * self%nu_ref(pol)/1.d9 * (T_e/1.d4)**(-1.5d0))) + 2.71828d0)
    !evalSED_freefree = S/S_ref * exp(-h*(nu-self%nu_ref(pol))/k_b/T_e) * (nu/self%nu_ref(pol))**(-2)
    !evalSED_freefree = S/S_ref * (nu/self%nu_ref(pol))**-2
    evalSED_freefree = S/S_ref * (nu/self%nu_ref(pol))**(-2)
    !write(*,*) "2",evalSED_freefree
    
  end function evalSED_freefree
  
end module comm_freefree_comp_mod
