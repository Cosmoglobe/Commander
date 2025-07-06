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
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_2D_mod
  use comm_data_mod
  implicit none

  private
  public comm_powlaw_break_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_powlaw_break_comp
     real(dp), allocatable, dimension(:) :: nu_break
   contains
     procedure :: S    => evalSED
  end type comm_powlaw_break_comp

  interface comm_powlaw_break_comp
     procedure constructor
  end interface comm_powlaw_break_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),           intent(in) :: cpar
    integer(i4b),                intent(in) :: id, id_abs
    class(comm_powlaw_break_comp), pointer   :: constructor

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

    ! Specific nu_break parameter
    allocate(constructor%nu_break(3))

    constructor%nu_break = cpar%cs_nu_break(id_abs,:)

    constructor%npar         = 2
    allocate(constructor%poltype(constructor%npar))
    do i = 1, constructor%npar
       constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call constructor%initLmaxSpecind(cpar, id, id_abs)

    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(constructor%theta_def(constructor%npar))
    allocate(constructor%p_gauss(2,constructor%npar))
    allocate(constructor%p_uni(2,constructor%npar))
    allocate(constructor%indlabel(constructor%npar))
    allocate(constructor%nu_min_ind(constructor%npar), constructor%nu_max_ind(constructor%npar))
    do i = 1, constructor%npar ! We don't want to initialize priors and stuff for nu_break
       constructor%theta_def(i) = cpar%cs_theta_def(i,id_abs)
       if (i == 3) then
          constructor%p_uni(:,i)   = cpar%cs_theta_def(i,id_abs)
          constructor%p_gauss(1,i) = cpar%cs_theta_def(i,id_abs)
          constructor%p_gauss(2,i) = 0.d0 ! Hard coded to 0 so that we never sample
       else
          constructor%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
          constructor%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
          constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
          constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
       end if
    end do
    constructor%indlabel = ['beta ','dbeta']

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar !We don't want to initialize maps and stuff for nu_break
       if (i == 3) then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(1)
       else
          ! Do we read in an initialization map for the 
          if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
             constructor%theta(i)%p => comm_map(info)
             constructor%theta(i)%p%map = constructor%theta_def(i)
          else
             ! Read map from FITS file
             constructor%theta(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(i,id_abs)))
          end if
       end if
       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (constructor%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call constructor%theta(i)%p%YtW_scalar
       end if
    end do

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (constructor%nu_ref(k) == constructor%nu_ref(k-1)) then
                   constructor%F_int(k,i,j)%p => constructor%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             constructor%F_int(k,i,j)%p => comm_F_int_2D(constructor, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    call constructor%initPixregSampling(cpar, id, id_abs)
    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !    SED  = (nu/nu_ref)**(beta+dbeta) if nu < nu_break
  !    SED  = (nu/nu_ref)**(beta) if nu => nu_break
  ! where 
  !    beta = theta(1), dbeta = theta(2)
  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_powlaw_break_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: beta,dbeta,nu_break
    real(dp)                                      :: evalSED

    beta     = theta(1)
    dbeta    = theta(2)
    nu_break = self%nu_break(pol)

    if (nu < nu_break) then
       evalSED = (nu/self%nu_ref(pol))**(beta)
    else
       evalSED = (nu/self%nu_ref(pol))**(beta+dbeta)
    end if

  end function evalSED
  
end module comm_powlaw_break_comp_mod
