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
module comm_powlaw_break_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_powlaw_break_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_powlaw_break_comp
     real(dp), allocatable, dimension(:) :: nu_break
   contains
     procedure :: S    => evalSED_powlaw_break
  end type comm_powlaw_break_comp

  interface comm_powlaw_break_comp
     procedure constructor_powlaw_break
  end interface comm_powlaw_break_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_powlaw_break(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),           intent(in) :: cpar
    integer(i4b),                intent(in) :: id, id_abs
    class(comm_powlaw_break_comp), pointer   :: c

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
    allocate(c)

    ! Specific nu_break parameter
    allocate(c%nu_break(3))

    c%nu_break = cpar%cs_nu_break(id_abs,:)

    c%npar         = 2
    allocate(c%poltype(c%npar))
    do i = 1, c%npar
       c%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call c%initLmaxSpecind(cpar, id, id_abs)

    call c%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(c%theta_def(c%npar))
    allocate(c%p_gauss(2,c%npar))
    allocate(c%p_uni(2,c%npar))
    allocate(c%indlabel(c%npar))
    allocate(c%theta_steplen(c%npar,cpar%mcmc_num_samp_groups))
    allocate(c%nu_min_ind(c%npar), c%nu_max_ind(c%npar))
    do i = 1, c%npar ! We don't want to initialize priors and stuff for nu_break
       c%theta_def(i) = cpar%cs_theta_def(i,id_abs)
       if (i == 3) then
          c%p_uni(:,i)   = cpar%cs_theta_def(i,id_abs)
          c%p_gauss(1,i) = cpar%cs_theta_def(i,id_abs)
          c%p_gauss(2,i) = 0.d0 ! Hard coded to 0 so that we never sample
       else
          c%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
          c%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
          c%nu_min_ind(i) = cpar%cs_nu_min_beta(id_abs,i)
          c%nu_max_ind(i) = cpar%cs_nu_max_beta(id_abs,i)
       end if
    end do
    c%theta_steplen = 0d0
    c%indlabel = ['beta','dbeta']

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, c%nside, c%lmax_ind, &
         & c%nmaps, c%pol)

    allocate(c%theta(c%npar))
    do i = 1, c%npar !We don't want to initialize maps and stuff for nu_break
       if (i == 3) then
          c%theta(i)%p => comm_map(info)
          c%theta(i)%p%map = c%theta_def(1)
       else
          ! Do we read in an initialization map for the 
          if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
             c%theta(i)%p => comm_map(info)
             c%theta(i)%p%map = c%theta_def(i)
          else
             ! Read map from FITS file
             c%theta(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(i,id_abs)))
          end if
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
             c%F_int(k,i,j)%p => comm_F_int_2D(c, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    call c%initPixregSampling(cpar, id, id_abs)
    ! Init alm 
    if (c%lmax_ind >= 0) call c%initSpecindProp(cpar, id, id_abs)

    ! Initialize mixing matrix
    call c%updateMixmat

  end function constructor_powlaw_break

  ! Definition:
  !    SED  = (nu/nu_ref)**(beta+dbeta) if nu < nu_break
  !    SED  = (nu/nu_ref)**(beta) if nu => nu_break
  ! where 
  !    beta = theta(1), dbeta = theta(2)
  function evalSED_powlaw_break(self, nu, band, pol, theta)
    implicit none
    class(comm_powlaw_break_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: beta,dbeta,nu_break
    real(dp)                                      :: evalSED_powlaw_break

    beta     = theta(1)
    dbeta    = theta(2)
    nu_break = self%nu_break(pol)

    if (nu < nu_break) then
       evalSED_powlaw_break = (nu/self%nu_ref(pol))**(beta)
    else
       evalSED_powlaw_break = (nu/self%nu_ref(pol))**(beta+dbeta)
    end if

  end function evalSED_powlaw_break
  
end module comm_powlaw_break_comp_mod
