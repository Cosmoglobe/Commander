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
module comm_pah_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_0D_mod
  use comm_data_mod
  use spline_2D_mod
  implicit none

  private
  public comm_pah_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_pah_comp
     real(dp)          :: nu_p0, nu_min, nu_max
     type(spline_type) :: SED_spline
   contains
     procedure :: S    => evalSED
  end type comm_pah_comp

  interface comm_pah_comp
     procedure constructor
  end interface comm_pah_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_pah_comp), pointer   :: constructor

    integer(i4b) :: ind(1), i, j, k
    real(dp), allocatable, dimension(:,:) :: SED


    ! General parameters
    allocate(constructor)
    constructor%npar = 0
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Initialize spectral template
    call read_spectrum(trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(1,id_abs)), SED)
    call spline(constructor%SED_spline, SED(:,1), SED(:,2))
    deallocate(SED)

    ! Precompute mixmat integrator for each band
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
             constructor%F_int(k,i,j)%p => comm_F_int_0D(constructor, data(i)%bp(j)%p, k)
          end do
       end do
    end do
    
    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_pah_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    evalSED = (splint(self%SED_spline, nu)) / &
            & (splint(self%SED_spline, self%nu_ref(pol)))! UNITS?? * (self%nu_ref(pol)/nu)**(2.d0)

  end function evalSED
  
end module comm_pah_comp_mod
