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
module comm_tsz_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_tsz_comp

  !**************************************************
  !      thermal Sunyaev-Zeldovich (MBB) component - subclass under diffuse class
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_tsz_comp
   contains
     procedure :: S    => evalSED_tsz
     procedure :: update_F_int => updateIntF_tsz
  end type comm_tsz_comp

  interface comm_tsz_comp
     procedure constructor_tsz
  end interface comm_tsz_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_tsz(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_MBB_comp), pointer   :: c

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    ! character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()

    ! General parameters
    allocate(c)

  function constructor_cmb(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_cmb_comp), pointer   :: c

    integer(i4b) :: i, j, k
    real(dp)     :: f
    

    ! General parameters
    allocate(c)
    c%npar         = 0
    allocate(c%poltype(1))
    c%poltype  = 1
    call c%initLmaxSpecind(cpar, id, id_abs)
    call c%initDiffuse(cpar, id, id_abs)


    ! Precompute mixmat integrator for each band
    allocate(c%F_int(3,numband,0:c%ndet))
    call c%update_F_int
    
    ! Initialize mixing matrix
    call c%updateMixmat

  end function constructor_tsz

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  =  x[(exp(x)+1)/(exp(x)-1)-4]*y (y is the sz map)

  function evalSED_tsz(self, nu, band, pol, theta)
    implicit none
    class(comm_MBB_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_tsz

    real(dp) :: x 
    
    
    if (nu>self%nu_max .or. nu<self%nu_min) then
      evalSED_mbbtab = 0.d0
      return
    end if 


    x   = h*nu / (k_b*T_CMB)
!! NEED to normalize to the reference frequency??!!!
    if (x < EXP_OVERFLOW) then
         evalSED_tsz = T_CMB * (x * (exp(x)+1.d0)/(exp(x)-1.d0)-4.d0) * (x * x * exp(x))/((exp(x)-1)*(exp(x)-1))
    else
         ! evalSED_tsz = 0.d0
         ! try to avoid overflow here but should probably have the comp_nu_max such that you don't 
         ! need this
         evalSED_tsz= T_CMB * (x-4.d0)
    endif

  end function evalSED_tsz
  

  !!! RAELYN WHAT IS THIS? IS THIS DOING THE RIGHT UNIT CONVERSION?
  ! Update band integration lookup tables
  subroutine updateIntF_tsz(self, band)
    implicit none
    class(comm_tsz_comp),                    intent(inout)        :: self
    integer(i4b),                            intent(in), optional :: band

    integer(i4b) :: i, j, k
    real(dp)     :: f

    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (self%nu_ref(k) == self%nu_ref(k-1)) then
                   self%F_int(k,i,j)%p => self%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             f = comp_a2t(self%nu_ref(k)) / data(i)%bp(j)%p%a2t * data(i)%bp(j)%p%RJ2data

             self%F_int(k,i,j)%p => comm_F_int_0D(self, data(i)%bp(j)%p, k, f_precomp=f)
          end do
       end do
    end do

  end subroutine updateIntF_tsz


end module comm_tsz_comp_mod
