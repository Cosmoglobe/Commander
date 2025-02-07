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
module comm_MBBtab_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_MBBtab_comp

  !**************************************************
  !      Modified Black Body (MBBtab) component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_MBBtab_comp
     character(len=128) :: mbbtab_type
     integer(i4b) :: npar_tab
     !real(dp), allocatable, dimension(:,:) :: SEDtab
     !real(dp), allocatable, dimension(:,:) :: SEDtab_buff
   contains
     procedure :: S    => evalSED_mbbtab
     procedure :: read_SED_table
  end type comm_MBBtab_comp

  interface comm_MBBtab_comp
     procedure constructor_mbbtab
  end interface comm_MBBtab_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_mbbtab(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_MBBtab_comp), pointer   :: c

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()

    ! General parameters
    allocate(c)

  end function constructor_mbbtab

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)
  ! where 
  !    beta = theta(1)
  function evalSED_mbbtab(self, nu, band, pol, theta)
    implicit none
    class(comm_MBBtab_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_mbbtab

    integer(i4b) :: i
    real(dp) :: x, x_ref, beta, T


    
    evalSED_mbbtab = 0

  end function evalSED_mbbtab


  ! Read precomputed SED table
  !   Each line in the file should contain {nu_min, nu_max, SED}
  subroutine read_SED_table(self, filename)
    implicit none
    class(comm_MBBtab_comp),    intent(inout)   :: self
    character(len=*),           intent(in)      :: filename
    
    integer(i4b) :: i, unit
    character(len=1024) :: line

    unit = getlun()
    self%ntab = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       self%ntab = self%ntab+1
    end do
1   close(unit)

    allocate(self%SEDtab(2+self%npar_tab,self%ntab))
    allocate(self%SEDtab_buff(2+self%npar_tab,self%ntab))
    open(unit, file=trim(filename))
    i = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       i = i+1
       read(line,*) self%SEDtab(:,i)
    end do
2   close(unit)

    self%SEDtab(1:2,:) = self%SEDtab(1:2,:) * 1d9
    self%SEDtab_buff = self%SEDtab

  end subroutine read_SED_table

  
end module comm_MBBtab_comp_mod
