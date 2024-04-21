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
     !integer(i4b) :: ntab
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


    c%npar         = 2
    allocate(c%poltype(c%npar))
    !if (cpar%myid == 0) write(*,*) cpar%cs_poltype(:,id_abs)
    do i = 1, c%npar
       c%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call c%initLmaxSpecind(cpar, id, id_abs)

    call c%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(c%theta_def(2), c%p_gauss(2,2), c%p_uni(2,2))
    allocate(c%indlabel(2))
    allocate(c%nu_min_ind(2), c%nu_max_ind(2))
    do i = 1, 2
       c%theta_def(i)  = cpar%cs_theta_def(i,id_abs)
       c%p_uni(:,i)    = cpar%cs_p_uni(id_abs,:,i)
       c%p_gauss(:,i)  = cpar%cs_p_gauss(id_abs,:,i)
       c%nu_min_ind(i) = cpar%cs_nu_min_beta(id_abs,i)
       c%nu_max_ind(i) = cpar%cs_nu_max_beta(id_abs,i)
    end do
    c%indlabel  = ['beta', 'T   ']

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
             c%F_int(k,i,j)%p => comm_F_int_2D(c, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, c%nside, c%lmax_ind, &
         & c%nmaps, c%pol)

    allocate(c%theta(c%npar))
    do i = 1, c%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
          c%theta(i)%p => comm_map(info)
          c%theta(i)%p%map = c%theta_def(i)
       else
          ! Read map from FITS file, and convert to alms
          c%theta(i)%p => comm_map(info, trim(cpar%cs_input_ind(i,id_abs)))
       end if

       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (c%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call c%theta(i)%p%YtW_scalar
       end if
    end do

    call c%initPixregSampling(cpar, id, id_abs)
    ! Init alm 
    if (c%lmax_ind >= 0) call c%initSpecindProp(cpar, id, id_abs)

    ! Read SED table
    call c%read_SED_table(cpar%cs_SED_template(1,id_abs))

    allocate(c%theta_steplen(2+c%ntab, cpar%mcmc_num_samp_groups))
    c%theta_steplen = 0d0

    ! Initialize SED priors
    c%SEDtab_prior = cpar%cs_SED_prior(id_abs)

    ! Initialize mixing matrix
    call c%updateMixmat

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


    ! First check if requested frequency is in tabulated range
    do i = 1, self%ntab
       if (nu > self%SEDtab(1,i) .and. nu <= self%SEDtab(2,i)) then
          ! Table defined in flux density
          evalSED_mbbtab = self%SEDtab(3,i) * (self%nu_ref(pol)/nu)**2
          return
       end if
    end do
    
    ! If not, use 
    beta    = theta(1)
    T       = theta(2)
    x       = h*nu               / (k_b*T)
    if (x > EXP_OVERFLOW) then
      evalSED_mbbtab = 0.d0
      return
    end if
    x_ref   = h*self%nu_ref(pol) / (k_b*T)
    evalSED_mbbtab = (nu/self%nu_ref(pol))**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)

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

    allocate(self%SEDtab(3,self%ntab))
    allocate(self%SEDtab_buff(3,self%ntab))
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
