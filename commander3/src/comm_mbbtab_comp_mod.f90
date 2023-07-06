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
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_2D_mod
  use comm_data_mod
  implicit none

  private
  public comm_MBBtab_comp

  !**************************************************
  !      Modified Black Body (MBBtab) component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_MBBtab_comp
     integer(i4b) :: ntab
     real(dp), allocatable, dimension(:,:) :: SEDtab
   contains
     procedure :: S    => evalSED
     procedure :: read_SED_table
  end type comm_MBBtab_comp

  interface comm_MBBtab_comp
     procedure constructor
  end interface comm_MBBtab_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_MBBtab_comp), pointer   :: constructor

    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()

    ! General parameters
    allocate(constructor)


    constructor%npar         = 2
    allocate(constructor%poltype(constructor%npar))
    do i = 1, constructor%npar
       constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call constructor%initLmaxSpecind(cpar, id, id_abs)

    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(constructor%theta_def(2), constructor%p_gauss(2,2), constructor%p_uni(2,2))
    allocate(constructor%indlabel(2))
    allocate(constructor%nu_min_ind(2), constructor%nu_max_ind(2))
    do i = 1, 2
       constructor%theta_def(i)  = cpar%cs_theta_def(i,id_abs)
       constructor%p_uni(:,i)    = cpar%cs_p_uni(id_abs,:,i)
       constructor%p_gauss(:,i)  = cpar%cs_p_gauss(id_abs,:,i)
       constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
       constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
    end do
    constructor%indlabel  = ['beta', 'T   ']

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

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(i)
       else
          ! Read map from FITS file, and convert to alms
          constructor%theta(i)%p => comm_map(info, trim(cpar%cs_input_ind(i,id_abs)))
       end if

       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (constructor%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call constructor%theta(i)%p%YtW_scalar
       end if
    end do

    call constructor%initPixregSampling(cpar, id, id_abs)
    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)

    ! Read SED table
    call constructor%read_SED_table(cpar%cs_SED_template(1,id_abs))

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_MBBtab_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i
    real(dp) :: x, x_ref, beta, T

    ! First check if requested frequency is in tabulated range
    do i = 1, self%ntab
       if (nu > self%SEDtab(1,i) .and. nu <= self%SEDtab(2,i)) then
          ! Table defined in flux density
          evalSED = self%SEDtab(3,i) * (self%nu_ref(pol)/nu)**2
          return
       end if
    end do
    
    ! If not, use 
    beta    = theta(1)
    T       = theta(2)
    x       = h*nu               / (k_b*T)
    if (x > EXP_OVERFLOW) then
      evalSED = 0.d0
      return
    end if
    x_ref   = h*self%nu_ref(pol) / (k_b*T)
    evalSED = (nu/self%nu_ref(pol))**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)

  end function evalSED


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

    allocate(self%SEDtab(3,self%ntab))
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

  end subroutine read_SED_table

  
end module comm_MBBtab_comp_mod
