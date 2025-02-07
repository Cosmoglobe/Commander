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
module comm_line_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_line_comp

  !**************************************************
  !           Line emission component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_line_comp
     integer(i4b)                            :: ref_band
     real(dp)                                :: line2RJ_ref
     integer(i4b), allocatable, dimension(:) :: ind2band
     real(dp),     allocatable, dimension(:) :: line2RJ
   contains
     procedure :: S    => evalSED_line
     procedure :: sampleSpecInd => sampleLineRatios
  end type comm_line_comp

  interface comm_line_comp
     procedure constructor_line
  end interface comm_line_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_line(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_line_comp), pointer   :: c

    integer(i4b) :: i, j, k, l, m, nline, b, n, ierr
    real(dp)     :: f
    logical(lgt) :: ref_exist
    character(len=512), allocatable, dimension(:) :: label
    real(dp),           allocatable, dimension(:) :: mu, sigma, line2RJ
    integer(i4b),       allocatable, dimension(:) :: poltype
    type(comm_mapinfo), pointer :: info
    
    ! General parameters
    allocate(c)

  end function constructor_line

  ! Definition:
  !    SED  = delta_{band,
  function evalSED_line(self, nu, band, pol, theta)
    class(comm_line_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_line

    integer(i4b) :: i, ind

    if (band == self%ref_band) then
       evalSED_line = 1.d0
    else 
       do i = 1, self%npar
          if (band == self%ind2band(i)) exit
       end do
       if (i > self%npar) then
          evalSED_line = 0.d0
       else
          evalSED_line = theta(i) * self%line2RJ(i) / self%line2RJ_ref
       end if
    end if

  end function evalSED_line

  subroutine read_line_template(filename, nline, label, mu, sigma, line2RJ, poltype)
    implicit none
    character(len=*),                              intent(in)  :: filename
    integer(i4b),                                  intent(out) :: nline
    character(len=512), allocatable, dimension(:), intent(out) :: label
    real(dp),           allocatable, dimension(:), intent(out) :: mu, sigma, line2RJ
    integer(i4b),       allocatable, dimension(:), intent(out) :: poltype

    integer(i4b)        :: i, unit
    character(len=1024) :: line

    unit  = getlun()

    ! Find number of lines
    nline = 0
    open(unit, file=trim(filename), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
    end do
1   close(unit)
    
    allocate(label(nline), mu(nline), sigma(nline), line2RJ(nline), poltype(nline))
    open(unit, file=trim(filename), recl=1024)
    nline = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
       read(line,*) label(nline), mu(nline), sigma(nline), line2RJ(nline), poltype(nline)
    end do
2   close(unit)

  end subroutine read_line_template

  ! Sample line ratios
  subroutine sampleLineRatios(self, cpar, handle, id, iter)
    implicit none
    class(comm_line_comp),                   intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b)    :: i, j, l, n, m, band, ierr
    real(dp)        :: A, b, mu, sigma, par, sigma_p, scale, w
    class(comm_map), pointer :: invN_amp, amp
    character(len=2) :: id_text
    

  end subroutine sampleLineRatios
  
end module comm_line_comp_mod
