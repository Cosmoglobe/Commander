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
module comm_B_bl_mod
  use comm_param_mod
  use comm_map_mod
  use comm_B_mod
  use comm_utils
  use spline_1D_mod
  implicit none

  private
  public comm_B_bl, comm_B_bl_ptr

  type, extends (comm_B) :: comm_B_bl
   contains
     ! Data procedures
     procedure :: conv           => matmulB
     procedure :: deconv         => matmulInvB
     procedure :: update         => updateBeam
  end type comm_B_bl

  interface comm_B_bl
     procedure constructor
  end interface comm_B_bl

  type comm_B_bl_ptr
     type(comm_B_bl), pointer :: p => null()
  end type comm_B_bl_ptr

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, fwhm, mb_eff, nside, pixwin, init_realspace)
    implicit none
    type(comm_params),                  intent(in)           :: cpar
    type(comm_mapinfo), target,         intent(in)           :: info
    integer(i4b),                       intent(in)           :: id, id_abs
    real(dp),                           intent(in), optional :: fwhm
    real(dp),                           intent(in), optional :: mb_eff
    character(len=*),                   intent(in), optional :: pixwin
    integer(i4b),                       intent(in), optional :: nside
    logical(lgt),                       intent(in), optional :: init_realspace
    class(comm_B_bl),   pointer                              :: constructor

    character(len=4)   :: nside_text
    character(len=512) :: dir
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    ! Component specific parameters
    constructor%type        =  'b_l'    
    constructor%info        => info
    constructor%almFromConv = .true.
    if (present(fwhm)) then
       allocate(constructor%b_l(0:constructor%info%lmax,constructor%info%nmaps))
!!$       do l = 0, constructor%info%lmax
!!$          constructor%b_l(l,:) = exp(-0.5d0*l*(l+1)*(fwhm*pi/180.d0/60/sqrt(8.d0*log(2.d0)))**2)
!!$       end do
       if (present(nside)) then
          call int2string(nside, nside_text)
          call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, fwhm=fwhm, &
               & pixwin=trim(dir)//'/pixel_window_n'//nside_text//'.fits')
       else if (present(pixwin)) then
          call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, fwhm=fwhm, &
               & pixwin=trim(dir)//'/'//trim(pixwin))
       else
          call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, fwhm=fwhm)
       end if
    else
       call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, &
            & beamfile=trim(cpar%ds_blfile(id_abs)), &
            & pixwin=trim(cpar%ds_pixwin(id_abs)))
    end if

    ! Multiply with main beam filling factor
    constructor%mb_eff = 1.d0; if (present(mb_eff)) constructor%mb_eff = mb_eff

    ! Initialize real-space profile
!!$    init_real = .true.; if (present(init_realspace)) init_real = init_realspace
!!$    if (init_real) then
!!$       call compute_radial_beam(constructor%info%nmaps, constructor%b_l, constructor%b_theta)
!!$       constructor%r_max = maxval(constructor%b_theta(1)%x)
!!$    end if
    
  end function constructor
  
  subroutine matmulB(self, trans, map)
    implicit none
    class(comm_B_bl), intent(in)    :: self
    logical(lgt),     intent(in)    :: trans
    class(comm_map),  intent(inout) :: map

    integer(i4b) :: i, j, l

    do i = 0, map%info%nalm-1
       l = map%info%lm(1,i)
       if (l <= self%info%lmax) then
          do j = 1, min(self%info%nmaps, map%info%nmaps)
             map%alm(i,j) = map%alm(i,j) * self%b_l(l,j) * self%mb_eff
          end do
       else
          map%alm(i,:) = 0.d0
       end if
    end do

  end subroutine matmulB

  subroutine matmulInvB(self, trans, map)
    implicit none
    class(comm_B_bl), intent(in)    :: self
    logical(lgt),     intent(in)    :: trans
    class(comm_map),  intent(inout) :: map

    integer(i4b) :: i, l, j

    do i = 0, map%info%nalm-1
       l = map%info%lm(1,i)
       if (l <= self%info%lmax) then
          do j = 1, min(self%info%nmaps, map%info%nmaps)
             if (self%b_l(l,j) > 1.d-12) then
                map%alm(i,j) = map%alm(i,j) / self%b_l(l,j) / self%mb_eff
             else
                map%alm(i,j) = 0.d0
             end if
          end do
       else
          map%alm(i,:) = 0.d0
       end if
    end do

  end subroutine matmulInvB

  subroutine updateBeam(self, b_l_norm, mb_eff) 
    implicit none
    class(comm_B_bl),                   intent(inout)           :: self
    real(dp),         dimension(0:,1:), intent(in),    optional :: b_l_norm
    real(dp),                           intent(in),    optional :: mb_eff

    if (present(mb_eff)) self%mb_eff = mb_eff
    
  end subroutine updateBeam
  
  
end module comm_B_bl_mod
