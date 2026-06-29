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
module comm_lcut_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_lcut_comp

  !**************************************************
  !           Lcut component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_lcut_comp
     integer(i4b)                            :: ref_band
   contains
     procedure :: S    => evalSED_lcut
     procedure :: update_F_int => updateIntF_lcut
  end type comm_lcut_comp

  interface comm_lcut_comp
     procedure constructor_lcut
  end interface comm_lcut_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_lcut(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_lcut_comp), pointer   :: c

    integer(i4b) :: i, ierr

    ! General parameters
    allocate(c)
    c%npar = 0 
    allocate(c%poltype(1))
    c%poltype  = 1
    call c%initLmaxSpecind(cpar, id, id_abs)
    call c%initDiffuse(cpar, id, id_abs)

    ! Find active band
    c%ref_band = -1
    do i = 1, numband
       if (trim(cpar%cs_band_ref(id_abs)) == trim(data(i)%label)) then
          c%ref_band = i
          exit
       end if
    end do
    if (c%ref_band == -1) then
       write(*,*) 'Warning: Reference band for lcut component ', id_abs, ' not found: ', trim(cpar%cs_band_ref(id_abs))
    end if

    ! Nullify other bands
    do i = 1, numband
       if (i /= c%ref_band) then
          c%F_null(i,:) = .true.          
       end if
    end do
    
    ! Precompute mixmat integrator for each band
    allocate(c%F_int(3,numband,0:c%ndet))
    call c%update_F_int
    
    ! Initialize mixing matrix
    call c%updateMixmat

  end function constructor_lcut

  ! Definition:
  !    SED  = delta_{band}
  function evalSED_lcut(self, nu, band, pol, theta)
    class(comm_lcut_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_lcut

    integer(i4b) :: i, ind

    if (band == self%ref_band) then
       evalSED_lcut = 1.d0
    else 
       evalSED_lcut = 0.d0
    end if

  end function evalSED_lcut

  ! Update band integration lookup tables
  subroutine updateIntF_lcut(self, band)
    implicit none
    class(comm_lcut_comp),                    intent(inout)        :: self
    integer(i4b),                            intent(in), optional :: band

    integer(i4b) :: i, j, k
    real(dp)     :: f

    do k = 1, 3
       i = self%ref_band
       do j = 0, data(i)%ndet
          if (k > 1) then
             if (self%nu_ref(k) == self%nu_ref(k-1)) then
                self%F_int(k,i,j)%p => self%F_int(k-1,i,j)%p
                cycle
             end if
          end if
          f = data(i)%bp(j)%p%RJ2data
             !if (.not. associated(self%F_int(k,i,j)%p)) then
          self%F_int(k,i,j)%p => comm_F_int_0D(self, data(i)%bp(j)%p, k, f_precomp=f)
             !else
             !   call self%F_int(k,i,j)%p%update(f_precomp=f)
             !end if
       end do
    end do

  end subroutine updateIntF_lcut
  
 end module comm_lcut_comp_mod
 
