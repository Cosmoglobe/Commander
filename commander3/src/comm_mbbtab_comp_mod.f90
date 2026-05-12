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
  use locate_mod 
  implicit none

  private
  public comm_MBBtab_comp

  !**************************************************
  !      Modified Black Body (MBBtab) component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_MBBtab_comp
     character(len=128) :: mbbtab_type
     integer(i4b) :: npar_tab, posneg  !npar_tab - how many columns in the table minus 2 
     real(dp)          :: nu_join ! frequency where to join the MBB and the tabulated values
     type(spline_type) :: spl
     type(spline_type) :: spl_buff
   contains
     procedure :: S    => evalSED_mbbtab
     procedure :: read_SED_table
     procedure :: read_astrodust_table
     procedure :: update_spline  !!procedure to update the spline 
     procedure :: update_spline_astrodust
     procedure :: write_spline !!procedure to print the spline 
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
    
    ! Set up MBBtab type
    c%mbbtab_type  = cpar%cs_mbbtab_type(id_abs)
    if (trim(c%mbbtab_type) == 'binned') then
       c%npar_tab = 1 ! 3 column table, nu_min, nu_max sed
       c%npar         = 2 !['beta', 'T   ']
    else if (trim(c%mbbtab_type) == 'linear') then
       c%npar_tab = 2  !4 column table
       c%npar         = 2 !['beta', 'T   ']
    else if (trim(c%mbbtab_type) == 'spline_log') then
       c%npar_tab = 1  ! 3 column table, nu_min, nu_max sed, same setup as binned
       c%npar = 2 !['beta', 'T   ']
    else if  (trim(c%mbbtab_type) == 'spline_astrodust') then 
       c%npar_tab = 0 ! 2 column table, nu_central, sed
       c%npar = 3 ! ['beta', 'T   ', 'adScale']
       c%nu_join = cpar%cs_nu_join(id_abs) * 1d9
    else
       write(*,*) 'Error: Unknown MBBtab type =', trim(c%mbbtab_type)
       stop
    end if


    allocate(c%poltype(c%npar))
    do i = 1, c%npar 
       c%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call c%initLmaxSpecind(cpar, id, id_abs)

    call c%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(c%theta_def(c%npar), c%p_gauss(2,c%npar), c%p_uni(2,c%npar))
    allocate(c%indlabel(c%npar))
    allocate(c%nu_min_ind(c%npar), c%nu_max_ind(c%npar))
    do i = 1, c%npar
       c%theta_def(i)  = cpar%cs_theta_def(i,id_abs)
       c%p_uni(:,i)    = cpar%cs_p_uni(id_abs,:,i)
       c%p_gauss(:,i)  = cpar%cs_p_gauss(id_abs,:,i)
       c%nu_min_ind(i) = cpar%cs_nu_min_beta(id_abs,i)
       c%nu_max_ind(i) = cpar%cs_nu_max_beta(id_abs,i)
    end do

    if  (trim(c%mbbtab_type) == 'spline_astrodust') then 
      c%indlabel  = ['beta', 'T   ', 'adScale']
    else
      c%indlabel  = ['beta', 'T   ']
    end if 

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


    !!!!!! RAELYN TO DO !!!!!!!!!
   if  (trim(c%mbbtab_type) == 'spline_astrodust') then
      call c%read_astrodust_table(cpar%cs_SED_template(2,id_abs))
   else
      c%nastrotab = 0
   end if 

    ! Make the initial spline if the right type
    if (trim(c%mbbtab_type) == 'spline_log') then
      call c%update_spline(c%theta_def(1), c%theta_def(2), 1) !1 is for pol, not set up for polarization currently
      c%spl_buff=c%spl
    else if (trim(c%mbbtab_type) == 'spline_astrodust')  then
      call c%update_spline_astrodust(c%theta_def(1), c%theta_def(2), c%theta_def(3), 1)
      c%spl_buff=c%spl
    end if

    allocate(c%theta_steplen(c%npar+c%ntab, cpar%mcmc_num_samp_groups))
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
    real(dp) :: x, x_ref, beta, T,maxnu,minnu, val,maxnu_ast
    
   ! nu, pol and theta are not in fact optional for this code so we should check we actually are getting those
    if (.not. present(nu)) then
         write(*,*) "nu missing in mbbtab"
         stop
      end if

      if (.not. present(pol)) then
         write(*,*) "pol missing in mbbtab"
         stop
      end if

      if (.not. present(theta)) then
         write(*,*) "theta missing in mbbtab"
         stop
      end if

    ! First check if requested frequency is in tabulated range
    ! SED table has nu_min nu_max SED or nu_central, SED for spline_astrodust
    ! SED table should be ordered from lowest to highest 
    if (self%ntab .ne. 0) then
       !SEDtab(1,1) !! minimum frequency in the spline, less than this uses MBB
       !! maximum frequency in the spline, no astrodust
      if (trim(self%mbbtab_type) == 'spline_astrodust') then 
         minnu=self%nu_join
         maxnu=self%astrotab(self%npar_tab+1,1)!self%SEDtab(self%npar_tab+1,self%ntab)
         maxnu_ast=self%astrotab(self%npar_tab+1,self%nastrotab) !! maximum frequency in the spline, with astrodust
      else 
         minnu=self%SEDtab(1,1)
         maxnu=self%SEDtab(self%npar_tab+1,self%ntab)
         maxnu_ast=maxnu-1 !set like this so that the evaluation into astrodust never occurs for non-astrodust cases
      end if 

      if (trim(self%mbbtab_type) == 'spline_log' .or. trim(self%mbbtab_type) == 'spline_astrodust') then
         if (nu<=minnu) then
            ! evaluates the low frequency range as a MBB
            beta    = theta(1)
            T       = theta(2)
            x       = h*nu               / (k_b*T)
            if (x > EXP_OVERFLOW) then
               evalSED_mbbtab = 0.d0
               return
            end if
            x_ref   = h*self%nu_ref(pol) / (k_b*T)
            evalSED_mbbtab = self%posneg*(nu/self%nu_ref(pol))**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)
            return
         else if (nu > minnu .and. nu <= maxnu) then 
            ! evaluates the spline for the tabulated values
            val = splint(self%spl, log(nu))
            if (.not. ieee_is_finite(val)) then
               if (self%x%info%myid == 0) then
                  write(*,*) "SPLINE RETURNED NAN, that is probably a bad thing"
                  write(*,*) "val = ", val
                  write(*,*) "nu =", nu
                  write(*,*) "log(nu) =", log(nu)
                  write(*,*) 'minnu = ', minnu
                  write(*,*) 'maxnu = ', maxnu
                  write(*,*) 'nuref = ', self%nu_ref(pol)
                  write(*,*) 'pol = ', pol
               end if 
               evalSED_mbbtab = 0.d0
            end if
            if (val > EXP_OVERFLOW) then
               evalSED_mbbtab = HUGE(0.d0)
               if (self%x%info%myid == 0) write(*,*) 'Warning, dust spline value is huge, possible unstable spline behaviour.'
            else if (val < -16) then
               evalSED_mbbtab = 0.d0
               if (self%x%info%myid == 0) write(*,*) 'Warning, dust spline value is very small, possible unstable spline behaviour.'
            else
               evalSED_mbbtab = self%posneg*exp(val) 
            end if 
            evalSED_mbbtab = evalSED_mbbtab * (self%nu_ref(pol)/nu)**2
            return
         else if (nu > maxnu .and. nu <= maxnu_ast) then 
            ! evaluates the astrotab values with linear between the bins, scaled by the astrodust scaling and posneg
            ! this should not evaluate if masnu_ast is set to maxnu-1, in which case there is no astrotab
            ! and higher than the maximum tabulated value will return zero. 
            i = locate_dp(self%astrotab(1,1:),nu)
            evalSED_mbbtab = self%posneg*theta(3)*(self%nu_ref(pol)/nu)**2 * (self%astrotab(2,i) + (nu - self%astrotab(1,i))*(self%astrotab(2,i) - self%astrotab(2,i+1))/(self%astrotab(1,i) - self%astrotab(1,i+1)))
            return
         else
            evalSED_mbbtab = 0.d0
            return
         end if
      else
         do i = 1, self%ntab
            if (nu > self%SEDtab(1,i) .and. nu <= self%SEDtab(2,i)) then
               ! Table defined in flux density
               if (trim(self%mbbtab_type) == 'binned') then
                  evalSED_mbbtab = self%SEDtab(3,i)
               else if (trim(self%mbbtab_type) == 'linear') then
                  evalSED_mbbtab = self%SEDtab(3,i) + self%SEDtab(4,i) * (nu-self%SEDtab(1,i))/(self%SEDtab(2,i)-self%SEDtab(1,i))
               else 
                  if (self%x%info%myid == 0) write(*,*) 'Error: MBBtab type unknown, =', trim(self%mbbtab_type)
                  stop
               end if
               evalSED_mbbtab = evalSED_mbbtab * (self%nu_ref(pol)/nu)**2
               return
            end if
         end do
      end if 
   end if

    beta    = theta(1)
    T       = theta(2)
    x       = h*nu               / (k_b*T)
    if (x > EXP_OVERFLOW) then
      evalSED_mbbtab = 0.d0
      return
    end if
    x_ref   = h*self%nu_ref(pol) / (k_b*T)
    evalSED_mbbtab = self%posneg*(nu/self%nu_ref(pol))**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)

  end function evalSED_mbbtab


  ! Read precomputed SED table
  !   Each line in the file should contain {nu_min, nu_max, SED}
  !   The units should be  Mj/sr/map_units where map_units are usually uK_rj@545
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
       if (len_trim(line) == 0) cycle ! fixes crash with empty lines in the table
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
       if (len_trim(line) == 0) cycle ! fixes crash with empty lines in the table
       i = i+1
       read(line,*) self%SEDtab(:,i)
    end do
2   close(unit)

    !!! astrodust type SED table only have the central frequency node and the amplitude, not two frequencies and the amplitude
    if (self%npar_tab==0) then 
      self%SEDtab(1,:) = self%SEDtab(1,:) * 1d9
    else 
      self%SEDtab(1:2,:) = self%SEDtab(1:2,:) * 1d9
    end if 
    self%SEDtab_buff = self%SEDtab
    self%posneg=1
  end subroutine read_SED_table

  subroutine read_astrodust_table(self, filename)
    implicit none
    class(comm_MBBtab_comp),    intent(inout)   :: self
    character(len=*),           intent(in)      :: filename
    
    integer(i4b) :: i, unit
    character(len=1024) :: line

    unit = getlun()
    self%nastrotab = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       if (len_trim(line) == 0) cycle 
       self%nastrotab = self%nastrotab+1
    end do
1   close(unit)

    allocate(self%astrotab(2,self%nastrotab)) 
    ! npar_tab is zero for astrodust and we want 2 columns, frequency and SED amplitude
    open(unit, file=trim(filename))
    i = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       if (len_trim(line) == 0) cycle 
       i = i+1
       read(line,*) self%astrotab(:,i)
    end do
2   close(unit)

    !!! astrodust type SED table only have the central frequency node and the amplitude, not two frequencies and the amplitude
    self%astrotab(1,:) = self%astrotab(1,:) * 1d9
  end subroutine read_astrodust_table



  subroutine update_spline(self,beta,T,pol)
    implicit none
    class(comm_MBBtab_comp),    intent(inout)   :: self
    real(dp), intent(in)                        :: T, beta
    integer(i4b),            intent(in)         :: pol
 
    
    integer :: i, n_pts
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: xnu, xnu_ref,nu,nu_ref,MBBbound,y_linear,Astbound
    character(512) :: filename


    ! If ntab < 3, we add 1 extra point to ensure spline stability
    if (self%ntab == 1) then
       n_pts = 3  ! One bin -> [Left, Mid, Right]
    else if (self%ntab == 2) then
       n_pts = 4  ! Two bins -> [Left1, Mid1, Mid2, Right2]
    else if (self%mbbtab_type == 'spline_log') then
       n_pts = self%ntab + 1 ! one extra point to link the MBB
    end if
    allocate(x(n_pts), y(n_pts))

   ! first match the low frequency edge of the first bin to the associated modified blackbody 
   ! the bins should be sorted from the lowest to highest frequencies
   ! you should not have zeros in the table, but if they are there you will get a warning and it will be set to 
   ! a very small value instead for the spline calculation - spline will be irratic though very likely

   ! computes the spline in log space
   
   nu=self%SEDtab(1,1)
   x(1)=log(nu)
   
   ! beta    = theta(1)
   ! T       = theta(2)
   nu_ref  = self%nu_ref(pol)
   xnu       = h*self%SEDtab(1,1) / (k_b*T)
   if (xnu > EXP_OVERFLOW) then
      y(1) = 0.d0
      if (self%x%info%myid == 0) then
      write(*,*) 'Error: MBB overflow in exponent, mbbTab'
      end if
      return
   end if
   xnu_ref   = h*nu_ref / (k_b*T)
   ! normalize the MBB in Mj/sr by the value at the reference frequency 

   y_linear=((nu)**(beta+3.d0)/(exp(xnu)-1.d0))/((nu_ref)**(beta+3.d0)/(exp(xnu_ref)-1.d0))
   y(1)=log(y_linear)
   
   ! Checks if the table is negative dust or not, only checks the first column, spline cannot cross the zero line
   ! so if the first line is negative they all should be
   self%posneg=INT(SIGN(1.0,self%SEDtab(3,1)))

   ! check if there are enough rows in the table, of there are only 1 or 2 we add extra bins of equal height 
   ! to help stabilize the spline, this is probably not optimal
   select case (self%ntab)
   case (1)
      ! Sample the only bin: Midpoint and Right Edge
      x(2) = log(0.5d0*(self%SEDtab(1,1) + self%SEDtab(2,1)))
      y(2) = log(max(1d-16, abs(self%SEDtab(3,1))))
      x(3) = log(self%SEDtab(2,1))
      y(3) = y(2)
   case (2)
      ! Sample both bins: Midpoint 1, Midpoint 2, and Right Edge 2
      x(2) = log(0.5d0*(self%SEDtab(1,1) + self%SEDtab(2,1)))
      y(2) = log(max(1d-16, abs(self%SEDtab(3,1))))
      x(3) = log(0.5d0*(self%SEDtab(1,2) + self%SEDtab(2,2)))
      y(3) = log(max(1d-16, abs(self%SEDtab(3,2))))
      x(4) = log(self%SEDtab(2,2))
      y(4) = y(3)

   case default
      ! choose the log of the midpoint of the bins for the spline nodes x
      ! and the logarithm of the SED tabulated values
      do i = 1, self%ntab-1
         if (abs(self%SEDtab(3,i))>1e-16) then !check for non-zero elements, don't want zeros in the spline? 
            x(i+1) = log(0.5d0*(self%SEDtab(1,i) + self%SEDtab(2,i)))
            y(i+1) = log(abs(self%SEDtab(3,i)))
         else
            x(i+1) = log(0.5d0*(self%SEDtab(1,i) + self%SEDtab(2,i)))
            y(i+1) = log(1d-16)
            if (self%x%info%myid == 0) then
               write(*,*) 'Warning, dust spline value is very small, did you forget a zero in your table? Possible unstable spline behaviour.'
            end if
         end if 
         !uncomment to debug if necessary
         ! if (self%x%info%myid == 0) then
         !    write(*,*) "spline nodes ", i+1, x(i+1), y(i+1)
         ! end if
      end do
      ! the last original point
      x(self%ntab+1)=log((self%SEDtab(2,self%ntab)))
      y(self%ntab+1) = log(abs(self%SEDtab(3,self%ntab)))
   end select


   ! the left boundary should match the first derivative between the MBB and the tabulated values
   ! the right boundary can be left to natural, set by 1e30
   MBBbound=(beta + 3.d0) - xnu * exp(xnu) / (exp(xnu) - 1.0)
   ! write(*,*) "MBBbound", MBBbound

   do i = 2, size(x)
      if (x(i) <= x(i-1)) then
         if (self%x%info%myid == 0) then
            write(*,*) "ERROR: mbbtab grid not strictly increasing, this will break the spline."
            write(*,*) "i=", i, "x(i-1)=", x(i-1), "x(i)=", x(i)
         end if
         stop
      end if
   end do

   call spline(self%spl, x, y, boundary=[MBBbound,1d30], regular=.false., linear=.false.)    
   deallocate(x,y)

    ! to debug set to true, warning! will output to the run folder
   !  if (.false.) then 
   !    filename = 'SEDdebug_.dat'
   !    call self%write_spline(filename,1,beta,T)
   !  end if 


  end subroutine  update_spline


  subroutine update_spline_astrodust(self,beta,T,adScale,pol)
    implicit none
    class(comm_MBBtab_comp),    intent(inout)   :: self
    real(dp), intent(in)                        :: T, beta, adScale
    integer(i4b),            intent(in)         :: pol
 
    
    integer :: i, n_pts,ind
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: xnu, xnu_ref,nu,nu_ref,MBBbound,y_linear,Astbound
    character(512) :: filename

    ind=1
    n_pts = self%ntab + 3 ! one extra point to link the MBB and 2 points to link the astrotab

    allocate(x(n_pts), y(n_pts))

    nu=self%nu_join !SEDtab(1,1) - 1e9 !shift the joining MBB frequency backwards by 1GHz to add to the spline. 
    x(ind)=log(nu)
    nu_ref  = self%nu_ref(pol)
    xnu       = h*nu / (k_b*T)
    if (xnu > EXP_OVERFLOW) then
         y(1) = 0.d0
         if (self%x%info%myid == 0) write(*,*) 'Error: MBB overflow in exponent, mbbTab'
         return
    end if
    xnu_ref   = h*nu_ref / (k_b*T)
    ! normalize the MBB in Mj/sr by the value at the reference frequency 

    y_linear=((nu)**(beta+3.d0)/(exp(xnu)-1.d0))/((nu_ref)**(beta+3.d0)/(exp(xnu_ref)-1.d0))
    y(ind)=log(y_linear)
      
    ! Checks if the table is negative dust or not, only checks the first column, spline cannot cross the zero line
    ! so if the first line is negative they all should be
    self%posneg=INT(SIGN(1.0,self%SEDtab(2,1)))

    do i = 1, self%ntab
        ind=ind+1
        if (abs(self%SEDtab(2,i))>1e-16) then !check for non-zero elements, don't want zeros in the spline? 
            x(ind) = log(self%SEDtab(1,i))
            y(ind) = log(abs(self%SEDtab(2,i)))
        else
            x(ind) = log(self%SEDtab(1,i))
            y(ind) = log(1d-16)
            if (self%x%info%myid == 0) then
               write(*,*) 'Warning, dust spline value is very small, did you forget a zero in your table? Possible unstable spline behaviour.'
            end if
        end if 
        if (x(ind) <= x(ind-1)) then 
            if (self%x%info%myid == 0) then
               write(*,*) "ERROR: mbbtab grid not strictly increasing, this will break the spline."
               write(*,*) "Also make sure nu_join is less than the start of the tabulated values."
               write(*,*) "i=", i, "x(i+1)=", x(ind), "x(i)=", x(ind-1), "nu_join=", self%nu_join
            end if
            stop
        end if 
    end do

      !add 2 points of the astrodust to the spline fit, and use those for the derivative of the RHS
    do i = 1, 2
        ind=ind+1
        if (abs(self%astrotab(2,i)*adScale)>1e-16) then !check for non-zero elements, don't want zeros in the spline? 
            x(ind) = log(self%astrotab(1,i))
            y(ind) = log(abs(self%astrotab(2,i)*adScale))
            ! if (self%x%info%myid == 0) write(*,*) 'x', self%astrotab(1,i), 'y' , self%astrotab(2,i)*adScale
        else
            x(ind) = log(self%astrotab(1,i))
            y(ind) = log(1d-16)
            if (self%x%info%myid == 0) then
               write(*,*) 'Warning, dust spline value is very small, did you forget a zero in your astrodust table? Possible unstable spline behaviour.'
            end if
        end if 
        if (x(ind) <= x(ind-1)) then 
            if (self%x%info%myid == 0) then
               write(*,*) "ERROR: Astrotab from mbbTab grid not strictly increasing, this will break the spline."
               write(*,*) "Probably there is illegal overlap between Astrotab and mbbTab."
               write(*,*) "i=", i, "x(i+1)=", x(ind), "x(i)=", x(ind-1)
            end if
            stop
        end if 
    end do

   !  if (self%x%info%myid == 0) then
   !    write(*,*) "Spline nodes"
   !    write(*,*) 'x', exp(x)
   !    write(*,*) 'y', exp(y)
   !    write(*,*) 'T,beta,adScale', T, beta, adScale
   !  end if

    ! the left boundary should match the first derivative between the MBB and the tabulated values
    ! the right boundary should match the slope of the line of the first two points of the astrodust table
    ! there should not be any frequency overlap between the two tables
    MBBbound = (beta + 3.d0) - xnu * exp(xnu) / (exp(xnu) - 1.0)
    Astbound = (y(ind)-y(ind-1))/(x(ind)-x(ind-1))



    call spline(self%spl, x, y, boundary=[MBBbound,Astbound], regular=.false., linear=.false.)
    

    
    deallocate(x,y)

    ! to debug set to true, warning! will output to the run folder
   !  if (.false.) then 
   !    filename = 'SEDdebug_.dat'
   !    call self%write_spline(filename,1,beta,T)
   !  end if 


  end subroutine  update_spline_astrodust




  !subroutine to write the spline, not currently used
  subroutine write_spline(self, filename, pol, beta, T, npts, nu_min, nu_max)
      implicit none
      class(comm_MBBtab_comp), intent(in) :: self
      character(len=*),        intent(in) :: filename
      integer(i4b),            intent(in) :: pol
      real(dp),                intent(in) :: beta, T
      integer(i4b), optional,  intent(in) :: npts
      real(dp),    optional,   intent(in) :: nu_min, nu_max

      integer(i4b) :: unit, i, N
      real(dp) :: nu1, nu2, dlognu, nu, sed
      real(dp), dimension(2)              :: theta

      theta(1) = beta
      theta(2) = T
      ! Defaults
      if (present(npts)) then
         N = npts
      else
         N = 500
      end if

      if (present(nu_min)) then
         nu1 = nu_min
      else
         nu1 = 30d0*1d9 !self%SEDtab(1,1)
      end if

      if (present(nu_max)) then
         nu2 = nu_max
      else
         nu2 = self%SEDtab(2,self%ntab)
      end if

      unit = getlun()
      open(unit, file=trim(filename), status='replace')

      write(unit,'(a)') '# nu[Hz]    SED'

      dlognu = (log(nu2) - log(nu1)) / real(N-1,dp)

      do i = 0, N-1
         nu  = exp(log(nu1) + dlognu*i)
         sed = self%S(nu=nu, pol=pol, theta=theta)
         write(unit,'(2E20.10)') nu, sed
      end do

      close(unit)

  end subroutine write_spline


end module comm_MBBtab_comp_mod
