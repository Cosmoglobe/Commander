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
  use comm_comp_interface_mod
  implicit none

  private
  public comm_pah_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_pah_comp
     real(dp)          :: nu_p0
     type(spline_type) :: SED_spline
   contains
     procedure :: S    => evalSED_pah
  end type comm_pah_comp

  interface comm_pah_comp
     procedure constructor_pah
  end interface comm_pah_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_pah(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_pah_comp), pointer   :: c

    integer(i4b) :: ind(1), i, j, k
    real(dp), allocatable, dimension(:,:) :: SED


    ! General parameters
    allocate(c)
    c%npar = 0
    call c%initDiffuse(cpar, id, id_abs)

    ! Initialize spectral template
    call read_spectrum(trim(cpar%cs_SED_template(1,id_abs)), SED)
    call spline(c%SED_spline, SED(:,1), SED(:,2))
    deallocate(SED)

    ! Precompute mixmat integrator for each band
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
             c%F_int(k,i,j)%p => comm_F_int_0D(c, data(i)%bp(j)%p, k)
          end do
       end do
    end do
    
    ! Initialize mixing matrix
    call c%updateMixmat

  end function constructor_pah

  function evalSED_pah(self, nu, band, pol, theta)
    implicit none
    class(comm_pah_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_pah

    if (nu < minval(self%SED_spline%x) .or.nu > maxval(self%SED_spline%x)) then
       evalSED_pah = 0.d0
    else
       evalSED_pah = (splint(self%SED_spline, nu)) / (splint(self%SED_spline, self%nu_ref(1))) !* (self%nu_ref(1)/nu)**(2.d0)
    end if
  end function evalSED_pah
  
end module comm_pah_comp_mod
