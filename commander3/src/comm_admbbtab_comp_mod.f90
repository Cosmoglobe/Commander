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
module comm_adMBBtab_comp_mod
  use comm_comp_interface_mod
  use spline_1d_mod
  use locate_mod 
  implicit none

  private
  public comm_adMBBtab_comp

  !**************************************************
  !   Modified Black Body + nodes + astrodust (adMBBtab) component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_adMBBtab_comp
   !   character(len=128) :: mbbtab_type
     integer(i4b) :: npar_tab, posneg  !npar_tab - how many columns in the table minus 2 
     real(dp)          :: nu_join,adScale,adscale_buff
     type(spline_type) :: spl
     type(spline_type) :: spl_buff

   contains
     procedure :: S    => evalSED_admbbtab
     procedure :: read_SED_table
     procedure :: read_astrodust_table
     procedure :: update_spline_astrodust
  end type comm_adMBBtab_comp

  interface comm_adMBBtab_comp
     procedure constructor_admbbtab
  end interface comm_adMBBtab_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_admbbtab(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_adMBBtab_comp), pointer   :: c

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
    c%npar_tab = 0 ! 2 column table, nu_central, sed
    c%npar = 2 ! ['beta', 'T   ']
    c%nu_join = cpar%cs_nu_join(id_abs) * 1d9   ! GHz -> Hz. <= 0 means "no join": pure MBB
    c%adscale = cpar%cs_adscale(id_abs)
    c%adscale_buff = cpar%cs_adscale(id_abs)

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

    c%indlabel  = ['beta', 'T   ']


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
    

    ! Read SED tables, the nodes and the astrodust model. Either can be 'none'
    call c%read_SED_table(cpar%cs_SED_template(1,id_abs))
    call c%read_astrodust_table(cpar%cs_SED_template(2,id_abs))

    ! Check that the three pieces fit together, and that the tables are consistent with nu_join
    if (c%nu_join <= 0.d0) then
       if (c%ntab > 0 .or. c%nastrotab > 0) then
          if (c%x%info%myid == 0) write(*,*) 'Warning: nu_join is not set, so NO spline or astrodust ', &
               & 'model is applied to the dust model (pure MBB). The SED tables are ignored.'
       end if
    else
       if (c%ntab == 0 .and. c%nastrotab == 0) then
          if (c%x%info%myid == 0) write(*,*) 'Error: nu_join is set but neither a nodes table nor ', &
               & 'an astrodust table is given, so there is nothing above nu_join.'
          stop
       end if
       if (c%ntab > 0) then
          if (c%nu_join >= c%SEDtab(1,1)) then
             if (c%x%info%myid == 0) write(*,*) 'Error: nu_join must be less than the smallest ', &
                  & 'frequency in the nodes table.'
             stop
          end if
       end if
       if (c%nastrotab > 0) then
          if (c%ntab > 0) then
             if (c%SEDtab(1,c%ntab) >= c%astrotab(1,1)) then
                if (c%x%info%myid == 0) write(*,*) 'Error: the highest frequency in the nodes table ', &
                     & 'must be smaller than the lowest frequency in the astrodust table.'
                stop
             end if
          else
             if (c%nu_join >= c%astrotab(1,1)) then
                if (c%x%info%myid == 0) write(*,*) 'Error: nu_join must be less than the lowest ', &
                     & 'frequency in the astrodust table.'
                stop
             end if
          end if
       end if
       if (c%nastrotab == 0 .and. c%ntab == 1) then
          if (c%x%info%myid == 0) write(*,*) 'Warning: only one node and no astrodust table, so no ', &
               & 'spline is computed; the tabulated value is used as a constant above nu_join.'
       end if
       ! The spline is normalized with nu_ref(pol), and is only built for one pol,
       ! so nu_ref has to be the same for all the Stokes parameters of this component
       do k = 2, 3
          if (c%nmaps >= k) then
             if (c%nu_ref(k) /= c%nu_ref(1)) then
                if (c%x%info%myid == 0) write(*,*) 'Error: adMBBtab needs the same nu_ref for all polarizations.'
                stop
             end if
          end if
       end do
    end if

    call c%update_spline_astrodust(c%theta_def(1), c%theta_def(2), c%adscale, 1)
    if (c%nu_join > 0.d0 .and. (c%nastrotab > 0 .or. c%ntab > 1)) c%spl_buff=c%spl
    allocate(c%theta_steplen(c%npar+c%ntab+1, cpar%mcmc_num_samp_groups))
 

    
    c%theta_steplen = 0d0

    ! Initialize SED priors
    c%SEDtab_prior = cpar%cs_SED_prior(id_abs)

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
    
    ! Initialize mixing matrix
    call c%updateMixmat

  end function constructor_admbbtab

  ! Definition:
  !      x  = h*nu/(k_b*T)
  !    SED  = (nu/nu_ref)**(beta+1) * (exp(x_ref)-1)/(exp(x)-1)       (nu <= nu_join)
  ! where 
  !    beta = theta(1)
  !    T    = theta(2)
  ! Above nu_join the spline and the astrodust table are used
  ! The spline is NOT rebuilt here: it uses the beta, T and adScale of the last call to
  ! update_spline_astrodust, while the astrodust part uses self%adscale.
  function evalSED_admbbtab(self, nu, band, pol, theta)
    implicit none
    class(comm_adMBBtab_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_admbbtab

    integer(i4b) :: i
    real(dp) :: x, x_ref, beta, T,maxnu,minnu, val,maxnu_ast, xmax, frac
    logical, save :: spline_warning_printed = .false.
    
   ! nu, pol and theta are not in fact optional for this code so we should check we actually are getting those
    if (.not. present(nu)) then
         write(*,*) "nu missing in admbbtab"
         stop
      end if

      if (.not. present(pol)) then
         write(*,*) "pol missing in admbbtab"
         stop
      end if

      if (.not. present(theta)) then
         write(*,*) "theta missing in admbbtab"
         stop
      end if

    if (nu>self%nu_max .or. nu<self%nu_min) then
      evalSED_admbbtab = 0.d0
      return
    end if 

    minnu=self%nu_join !! up to nu_join is MBB

    ! Modified blackbody: everywhere if there is no join frequency, otherwise up to nu_join
    if (minnu <= 0.d0 .or. nu <= minnu) then
       beta    = theta(1)
       T       = theta(2)
       x       = h*nu               / (k_b*T)
       if (x > EXP_OVERFLOW) then
          evalSED_admbbtab = 0.d0
          return
       end if
       x_ref   = h*self%nu_ref(pol) / (k_b*T)
       evalSED_admbbtab = (nu/self%nu_ref(pol))**(beta+1.d0) * (exp(x_ref)-1.d0)/(exp(x)-1.d0)
       return
    end if

    ! Above nu_join. The spline extends up to the first astrodust point, or, without an astrodust
    ! table, up to the last node. Beyond the last node the value is held constant.
    if (self%nastrotab > 0) then
       maxnu     = self%astrotab(1,1)                    !! up to the first point in astrodust is the spline
       maxnu_ast = self%astrotab(1,self%nastrotab)       ! maximum frequency in astrodust
       xmax      = log(maxnu)
    else
       maxnu     = HUGE(0.d0)                            ! no astrodust: spline region goes on to nu_max
       maxnu_ast = HUGE(0.d0)
       xmax      = log(self%SEDtab(1,self%ntab))
    end if

    if (nu <= maxnu) then
       if (self%nastrotab == 0 .and. self%ntab == 1) then
          ! one node and no astrodust: no spline, constant tabulated value
          evalSED_admbbtab = self%SEDtab(2,1)
       else
          ! evaluates the spline for the tabulated values, which holds asinh(amplitude)
          val = splint(self%spl, min(log(nu), xmax))
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
             evalSED_admbbtab = 0.d0
             return
          end if
          ! sinh is odd, so overflow can happen on either side. There is no "too small" case:
          ! val near 0 just means amplitude near 0, which is fine. The flag only limits printing.
          if (abs(val) > EXP_OVERFLOW) then
             evalSED_admbbtab = SIGN(HUGE(0.d0), val)
             if (self%x%info%myid == 0 .and. .not. spline_warning_printed) &
                  & write(*,*) 'Warning, dust spline value is huge, possible unstable spline behaviour.'
             spline_warning_printed = .true.
          else
             evalSED_admbbtab = sinh(val)
          end if
       end if
       evalSED_admbbtab = evalSED_admbbtab * (self%nu_ref(pol)/nu)**2
       return
    else if (nu <= maxnu_ast) then 
       ! evaluates the astrotab values with linear interpolation between the points, scaled by the
       ! astrodust scaling. Higher than the maximum tabulated frequency returns zero.
       i = locate_dp(self%astrotab(1,1:),nu)
       i = max(1, min(i, self%nastrotab-1))    ! also safe at the very ends of the table
       frac = (nu - self%astrotab(1,i))/(self%astrotab(1,i+1) - self%astrotab(1,i))
       evalSED_admbbtab = self%adscale*(self%nu_ref(pol)/nu)**2 * &
            & (self%astrotab(2,i) + frac*(self%astrotab(2,i+1) - self%astrotab(2,i)))
       return
    else
       evalSED_admbbtab = 0.d0
       return
    end if

  end function evalSED_admbbtab


  ! Read precomputed SED table of spline nodes
  !   Each line in the file should contain {nu, SED}, nu in GHz
  !   The units should be  Mj/sr/map_units where map_units are usually uK_rj@545
  !   The file name 'none' (or empty) means no nodes
  !   The table is sorted by frequency (with a warning) if it isn't already
  subroutine read_SED_table(self, filename)
    implicit none
    class(comm_adMBBtab_comp),    intent(inout)   :: self
    character(len=*),           intent(in)      :: filename
    
    self%posneg=1
    if (trim(filename) == 'none' .or. len_trim(filename) == 0) then
       self%ntab = 0
       allocate(self%SEDtab(2,0))
       allocate(self%SEDtab_buff(2,0))
       return
    end if

    call count_table_rows(filename, self%ntab)

    allocate(self%SEDtab(2,self%ntab)) ! freq, SED
    allocate(self%SEDtab_buff(2,self%ntab))
    call read_table_rows(filename, 'nodes', self%SEDtab)

    !!! astrodust type SED table only have the central frequency node and the amplitude, not two frequencies and the amplitude
    self%SEDtab(1,:) = self%SEDtab(1,:) * 1d9
    call check_and_sort_table(self%SEDtab, 'nodes', self%x%info%myid)
    self%SEDtab_buff = self%SEDtab
  end subroutine read_SED_table

  ! Read the astrodust table. Each line in the file should contain {nu, SED}, nu in GHz
  !   The file name 'none' (or empty) means no astrodust
  !   Needs at least two rows, and is sorted by frequency (with a warning) if it isn't already
  subroutine read_astrodust_table(self, filename)
    implicit none
    class(comm_adMBBtab_comp),    intent(inout)   :: self
    character(len=*),           intent(in)      :: filename
    
    if (trim(filename) == 'none' .or. len_trim(filename) == 0) then
       self%nastrotab = 0
       allocate(self%astrotab(2,0))
       return
    end if

    call count_table_rows(filename, self%nastrotab)

    allocate(self%astrotab(2,self%nastrotab)) 
    ! 2 columns, frequency and SED amplitude
    call read_table_rows(filename, 'astrodust', self%astrotab)

    !!! astrodust type SED table only have the central frequency node and the amplitude, not two frequencies and the amplitude
    self%astrotab(1,:) = self%astrotab(1,:) * 1d9

    ! the slope of the first two points is used to join the spline to the astrodust model
    if (self%nastrotab < 2) then
       if (self%x%info%myid == 0) write(*,*) 'Error: the astrodust table needs at least two rows.'
       stop
    end if
  end subroutine read_astrodust_table


  ! Build the spline that connects the MBB to the astrodust model, or to the end of the nodes.
  ! Has to be called again whenever beta, T or adScale change (proposals, and restoring after a
  ! rejection: copy spl_buff back and reset adscale from adscale_buff).
  ! The spline lives in (ln nu, asinh(amplitude)) space. The asinh keeps it well behaved over many
  ! orders of magnitude and, unlike the log, it preserves the sign of each node, so tables that
  ! change sign are fine. S undoes it with sinh.
  ! Knots: the MBB at nu_join, all the nodes, and the first astrodust point (if there is astrodust).
  ! Boundary conditions (first derivative in ln nu):
  !    left  = slope of the MBB at nu_join
  !    right = slope of the scaled astrodust table at its first point (the slope of the linear
  !            interpolation used in S), or 0 if there is no astrodust table
  ! Nothing is built without nu_join, or for a single node without astrodust (constant in S).
  subroutine update_spline_astrodust(self,beta,T,adScale,pol)
    implicit none
    class(comm_adMBBtab_comp),    intent(inout)   :: self
    real(dp), intent(in)                        :: T, beta, adScale
    integer(i4b),            intent(in)         :: pol
 
    
    integer(i4b) :: i, n_pts
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: xnu, xnu_ref, nu, nu_ref, f, dlnf_dlnnu, left_slope, right_slope
    real(dp) :: nu0, nu1, I0, I1, dI_dnu

    if (self%nu_join <= 0.d0) return
    if (self%nastrotab == 0 .and. self%ntab <= 1) return

    nu      = self%nu_join
    nu_ref  = self%nu_ref(pol)
    xnu     = h*nu     / (k_b*T)
    xnu_ref = h*nu_ref / (k_b*T)
    if (max(xnu, xnu_ref) > EXP_OVERFLOW) then
         if (self%x%info%myid == 0) write(*,*) 'Error: MBB overflow in exponent, adMBBtab'
         stop
    end if

    ! MBB at nu_join relative to nu_ref, and the derivative of asinh(f) with respect to ln(nu):
    !   d asinh(f)/d ln(nu) = f' / sqrt(1+f**2),  f' = f * d ln(f)/d ln(nu)
    f          = (nu/nu_ref)**(beta+3.d0) * (exp(xnu_ref)-1.d0)/(exp(xnu)-1.d0)
    dlnf_dlnnu = (beta + 3.d0) - xnu * exp(xnu) / (exp(xnu) - 1.d0)
    left_slope = f*dlnf_dlnnu / sqrt(1.d0 + f**2)

    n_pts = 1 + self%ntab
    if (self%nastrotab > 0) n_pts = n_pts + 1  ! the first astrodust point closes the spline
    allocate(x(n_pts), y(n_pts))

    ! first knot: the MBB at the join frequency
    x(1) = log(nu)
    y(1) = asinh(f)

    ! middle knots: the nodes (sorted and checked to lie above nu_join when they were read)
    do i = 1, self%ntab
       x(1+i) = log(self%SEDtab(1,i))
       y(1+i) = asinh(self%SEDtab(2,i))
    end do

    if (self%nastrotab > 0) then
       ! last knot: the first astrodust point, with the slope of the linear interpolation just above it
       nu0 = self%astrotab(1,1)
       nu1 = self%astrotab(1,2)
       I0  = self%astrotab(2,1)
       I1  = self%astrotab(2,2)
       dI_dnu = (I1 - I0)/(nu1 - nu0)
       right_slope = adScale*nu0*dI_dnu / sqrt(1.d0 + (adScale*I0)**2)
       x(n_pts) = log(nu0)
       y(n_pts) = asinh(adScale*I0)
    else
       right_slope = 0.d0   ! no astrodust: the spline ends flat
    end if

    call spline(self%spl, x, y, boundary=[left_slope,right_slope], regular=.false., linear=.false.)
    
    deallocate(x,y)

  end subroutine  update_spline_astrodust


  !**************************************************
  !     Table helpers (module private)
  !**************************************************

  ! Number of data rows in a text table (comment lines starting with # and empty lines are skipped)
  subroutine count_table_rows(filename, n)
    implicit none
    character(len=*), intent(in)  :: filename
    integer(i4b),     intent(out) :: n

    integer(i4b) :: unit
    character(len=1024) :: line

    unit = getlun()
    n = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       if (len_trim(line) == 0) cycle ! fixes crash with empty lines in the table
       n = n+1
    end do
1   close(unit)
  end subroutine count_table_rows


  ! Read a two-column text table into tab(2,n), n being the count from count_table_rows.
  ! Stops if a line doesn't have exactly two numbers (e.g. an old {nu_min, nu_max, SED} table
  ! would otherwise be read with nu_max as the amplitude)
  subroutine read_table_rows(filename, label, tab)
    implicit none
    character(len=*), intent(in)  :: filename, label
    real(dp),         intent(out) :: tab(:,:)

    integer(i4b) :: i, unit, ios
    character(len=1024) :: line
    real(dp) :: extra(3)

    unit = getlun()
    open(unit, file=trim(filename))
    i = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       if (len_trim(line) == 0) cycle ! fixes crash with empty lines in the table
       i = i+1
       read(line,*,iostat=ios) tab(:,i)
       if (ios /= 0) then
          write(*,*) 'Error: could not read two numbers (frequency, amplitude) from the ', trim(label), &
               & ' table ', trim(filename), ', line: ', trim(line)
          stop
       end if
       read(line,*,iostat=ios) extra
       if (ios == 0) then
          write(*,*) 'Error: the ', trim(label), ' table ', trim(filename), &
               & ' must have exactly two columns (frequency, amplitude), line: ', trim(line)
          stop
       end if
    end do
2   close(unit)
  end subroutine read_table_rows


  ! Frequencies (row 1) must be positive and unique. Sorts the table by frequency, with a warning,
  ! if it is not ordered already
  subroutine check_and_sort_table(tab, label, myid)
    implicit none
    real(dp),         intent(inout) :: tab(:,:)
    character(len=*), intent(in)    :: label
    integer(i4b),     intent(in)    :: myid

    integer(i4b) :: i, j, n
    real(dp)     :: tmp(2)

    n = size(tab,2)
    if (any(tab(1,:) <= 0.d0)) then
       if (myid == 0) write(*,*) 'Error: the frequencies in the ', trim(label), ' table must be positive.'
       stop
    end if

    if (.not. all(tab(1,2:n) >= tab(1,1:n-1))) then
       if (myid == 0) write(*,*) 'Warning: the ', trim(label), ' table is not ordered by frequency, sorting it.'
       do i = 2, n   ! insertion sort
          tmp = tab(:,i)
          j = i-1
          do while (j >= 1)
             if (tab(1,j) <= tmp(1)) exit
             tab(:,j+1) = tab(:,j)
             j = j-1
          end do
          tab(:,j+1) = tmp
       end do
    end if

    if (any(tab(1,2:n) <= tab(1,1:n-1))) then
       if (myid == 0) write(*,*) 'Error: the ', trim(label), ' table has duplicate frequencies.'
       stop
    end if
  end subroutine check_and_sort_table


end module comm_adMBBtab_comp_mod
