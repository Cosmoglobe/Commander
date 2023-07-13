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
module comm_ame_lognormal_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_2D_mod
  use comm_data_mod
  implicit none

  private
  public comm_ame_lognormal_comp

  !**************************************************
  !        AME - lognormal width component
  !  https://arxiv.org/pdf/2001.07159.pdf (Eqn. 6)
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_ame_lognormal_comp
     contains
       procedure :: s => evalSED
  end type comm_ame_lognormal_comp

  interface comm_ame_lognormal_comp
     procedure constructor
  end interface comm_ame_lognormal_comp

contains

  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),  intent(in) :: cpar
    integer(i4b),       intent(in) :: id, id_abs
    class(comm_ame_lognormal_comp), pointer :: constructor

    integer(i4b) :: i, j, k
    type(comm_mapinfo), pointer :: info => null()

    allocate(constructor)

    constructor%npar     = 2
    allocate(constructor%poltype(constructor%npar))
    do i = 1, constructor%npar
       constructor%poltype(i)  = cpar%cs_poltype(i,id_abs)
    end do
    call constructor%initLmaxSpecInd(cpar, id, id_abs)

    call constructor%initDiffuse(cpar,id,id_abs)

    ! Component specific parameters

    allocate(constructor%theta_def(2), constructor%p_gauss(2,2), constructor%p_uni(2,2))
    allocate(constructor%indlabel(2))
    allocate(constructor%nu_min_ind(2), constructor%nu_max_ind(2))
    do i = 1, 2
       constructor%theta_def(i)  = cpar%cs_theta_def(i,id_abs)
       constructor%p_uni(:,i)    = cpar%cs_p_uni(id_abs,:,i)
       constructor%p_gauss(:,i)  = cpar%cs_p_gauss(id_abs,:,i)
       constructor%nu_min_ind(i) = cpar%cs_nu_min_beta(id_abs,i)
       constructor%nu_max_ind(i) = cpar%cs_nu_max_beta(id_abs,i)
    end do
    constructor%indlabel = ['nu_p', 'W_AME']

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' ) then
          constructor%theta(i)%p     => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(i)
       else
          ! Read map from FITS file, and convert to alms
          constructor%theta(i)%p => comm_map(info, trim(cpar%cs_input_ind(i,id_abs)))
       end if

       ! convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (constructor%lmax_ind >= 0) then
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

    ! Init mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !
  ! S_AME(nu) = A_AME * exp{ -1/2 * [ln(nu/nu_AME)/W_AME]**2}
  !
  ! where nu_AME is the peak frequency and W_AME is the width 
  ! of the spectrum

  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_ame_lognormal_comp),    intent(in) :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED
    real(dp)                                      :: peak
    real(dp)                                      :: nu_AME, W_AME

    nu_AME  = theta(1)
    W_AME   = theta(2)

    peak    = nu_AME*1d9

    ! Moving from flux density to T_RJ requires an multiplicative (1/nu^2) 
    
    evalSED = exp(-0.5*(log(nu/peak)/W_AME)**2)*(self%nu_ref(pol)/nu)**2
  end function evalSED
end module comm_ame_lognormal_mod
