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
  use spline_1d_mod
  implicit none

  private
  public comm_MBBtab_comp

  !**************************************************
  !      Modified Black Body (MBBtab) component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_MBBtab_comp
     character(len=128) :: mbbtab_type
     integer(i4b) :: npar_tab,posneg
     !real(dp), allocatable, dimension(:,:) :: SEDtab
     !real(dp), allocatable, dimension(:,:) :: SEDtab_buff
     type(spline_type) :: spl
     type(spline_type) :: spl_buff
   contains
     procedure :: S    => evalSED_mbbtab
     procedure :: read_SED_table
     procedure :: update_spline  !!procedure to update the spline 
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


    c%npar         = 2 !!first 2 npar for the T and beta?
    allocate(c%poltype(c%npar))
    !if (cpar%myid == 0) write(*,*) cpar%cs_poltype(:,id_abs)
    do i = 1, c%npar
       c%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call c%initLmaxSpecind(cpar, id, id_abs)

    call c%initDiffuse(cpar, id, id_abs)
    
    ! Set up MBBtab type
    c%mbbtab_type  = cpar%cs_mbbtab_type(id_abs)
    if (trim(c%mbbtab_type) == 'binned') then
       c%npar_tab = 1
    else if (trim(c%mbbtab_type) == 'linear') then
       c%npar_tab = 2
    else if (trim(c%mbbtab_type) == 'spline_log') then
       c%npar_tab = 3        
    else
       write(*,*) 'Error: Unknown MBBtab type =', trim(c%mbbtab_type)
       stop
    end if
    
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

    ! Make the initial spline if the right type
    if ((trim(c%mbbtab_type) == 'spline_log')) then
      call c%update_spline()
      c%spl_buff=c%spl
    end if
    
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
    real(dp) :: x, x_ref, beta, T,maxnu,minnu


    ! First check if requested frequency is in tabulated range
    ! SED table has nu_min nu_max SED
    maxnu=-1.0d0
    minnu=-1.0d0
    do i = 1, self%ntab
       if (maxnu < self%SEDtab(2,i)) then
         maxnu=self%SEDtab(2,i)
       end if

       if (minnu > self%SEDtab(1,i)) then
         minnu=self%SEDtab(1,i)
       end if 

       if (nu > self%SEDtab(1,i) .and. nu <= self%SEDtab(2,i)) then
          ! Table defined in flux density
          if (trim(self%mbbtab_type) == 'binned') then
             evalSED_mbbtab = self%SEDtab(3,i)
          else if (trim(self%mbbtab_type) == 'linear') then
             evalSED_mbbtab = self%SEDtab(3,i) + self%SEDtab(4,i) * (nu-self%SEDtab(1,i))/(self%SEDtab(2,i)-self%SEDtab(1,i))
          else if (trim(self%mbbtab_type) == 'spline_log') then
             evalSED_mbbtab = self%posneg*exp(splint(self%spl, log(nu))) 
          else 
            write(*,*) 'Warning: MBBtab type =', trim(self%mbbtab_type)
            write(*,*) 'Evaluating as if a MBB'
          end if
          !!should the units of the table be checked first?
          evalSED_mbbtab = evalSED_mbbtab * (self%nu_ref(pol)/nu)**2
          return
       end if
    end do
    
    ! RAELYN ADD QUIT FOR NU BIGGER? Or follow the Spline? Maybe checks and balances
    ! if we reach here then nu is not in a tabulated range
    if (nu <= minnu) then
       write(*,*) 'Warning: using MBB, nu less than minimum tabulated value = ', minnu
    else if (nu > maxnu) then
       write(*,*) 'Error: nu greater than the tabulated values = ', maxnu
       write(*,*) 'Should have tabulated values up to the desired nu in table'
       stop
    else
       write(*,*) 'Error: nu greater or less than the tabulated values (???), using MBB'
       write(*,*) 'Should not reach this point, quitting'
       stop
    end if


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
  !   The units should be Mj/sr/map_units where map_units are usually uK_rj@545
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


  subroutine update_spline(self)
    implicit none
    class(comm_MBBtab_comp),    intent(inout)   :: self
    
    integer :: i
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: xnu, xnu_ref, beta, T,nu,nu_ref,pol,MBBbound,y_linear

    allocate(x(self%ntab+1), y(self%ntab+1))


    if (self%mbbtab_type == 'spline_log') then
      ! first match the low frequency edge of the first bin to the associated modified blackbody 
      ! (this assumes the bins are sorted from the lowest to highest frequencies)
      ! computes the spline in log space
      
      nu=self%SEDtab(1,1)
      x(1)=log(nu)
      
      beta    = self%theta_def(1)
      T       = self%theta_def(2)
      write(*,*) 'beta, T = ', beta ,T
      pol     = self%pol
      nu_ref  = self%nu_ref(pol)
      xnu       = h*self%SEDtab(1,1) / (k_b*T)
      if (xnu > EXP_OVERFLOW) then
        y(1) = 0.d0
        write(*,*) 'Error: MBB overflow in exponent'
        return
      end if
      xnu_ref   = h*nu_ref / (k_b*T)
      !normalize the MBB in Mj/sr by the value at the reference frequency 

      y_linear=(nu)**(beta+3.d0)/(exp(xnu)-1.d0)/((nu_ref)**(beta+3.d0)/(exp(xnu_ref)-1.d0))
      y(1)=log(y_linear)
      !!Checks if the table is negative dust or not, only checks the first column, if the dust
      !changes signs this will not work (would need a new spline function to handle this zero crossing)
      self%posneg=INT(SIGN(1.0,self%SEDtab(3,1)))
      ! choose the logarithmic midpoint of the bins for the spline nodes x
      ! and the logarithm of the SED tabulated values
      do i = 1, self%ntab
         if (abs(self%SEDtab(3,i))>1e-16) then !check for non-zero elements, don't want zeros in the spline? 
            x(i+1) = log(0.5d0*(self%SEDtab(1,i) + self%SEDtab(2,i)))
            y(i+1) = log(abs(self%SEDtab(3,i)))
         end if 
      end do

      ! the left boundary should match the first derivative between the MBB and the tabulated values
      ! the right boundary can be left to natural, set by 1e30
      MBBbound=(beta + 3.d0) - xnu * exp(xnu) / (exp(xnu) - 1.0)

      call spline(self%spl, x, y, &
                  boundary=[MBBbound,1d30], &
                  regular=.false., &
                  linear=.false.)
    else 
      write(*,*) 'Warning: MBBtab spline type unknown = ', trim(self%mbbtab_type)
    end if
    
    deallocate(x,y)



    !First derivative
    ! (exp(x_ref)-1.d0)/(self%nu_ref(pol)**(beta+1.d0))*(nu**beta)/(exp(x)-1.d0)
    ! *(
    ! (beta+1.d0)-x*exp(x)/(exp(x)-1)
    ! )

    !Second derivative
    !(exp(x_ref)-1.d0)/(self%nu_ref(pol)**(beta+1.d0))*(nu**(beta-1.d0))/(exp(x)-1.d0)
    ! *(
    !((beta+1.d0) * beta)
    ! -2.d0*(beta+1.d0)*x*exp(x)/(exp(x)-1.d0)
    ! +x**2.d0 * exp(x) * (exp(x)+1.d0) / (exp(x)-1.d0)**2.d0
    ! )
  end subroutine  update_spline

  
end module comm_MBBtab_comp_mod

