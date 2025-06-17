!===============================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University
! of Oslo.
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
!===============================================================================

! This module handles corrections to time ordered data due to adc issues

submodule (comm_tod_adc_binfit_mod) comm_tod_adc_binfit_smod
  contains
  
  module function constructor_adc_binfit(comm, label, nbit, min_adu, max_adu, ncoadd) result(c)
    ! ====================================================================
    ! Sets up an adc correction object that maps 
    !
    ! Inputs:
    ! 
    ! comm   : integer
    !          mpi communicator
    !
    ! myid   : integer
    !          mpi identifier
    !
    ! nbins  : integer
    !          number of bins used for building the adc correction tables
    !
    ! Returns:
    ! --------
    !
    ! constructor_internal: pointer
    !    contains all of the bins needed for computing adc corrections
    !    and the actual correction tables
    ! ====================================================================
    implicit none
    character(len=*),       intent(in) :: label
    integer(i4b),           intent(in) :: comm, nbit, min_adu, max_adu, ncoadd
    class(comm_adc_binfit), pointer    :: c

    integer(i4b) :: i, j, c0, ierr
    real(dp)     :: delta
    
    allocate(c)
    c%comm      = comm
    call mpi_comm_rank(comm, c%myid, ierr)
    c%label     = label
    c%nbit      = nbit
    c%min_adu   = min_adu
    c%max_adu   = max_adu
    c%ncoadd    = ncoadd
    c%min_coadd = min_adu*ncoadd
    c%max_coadd = max_adu*ncoadd

    allocate(c%Q(min_adu:max_adu))
    allocate(c%v(min_adu:max_adu+1)) ! Lower limit of each bin
    allocate(c%INL(min_adu:max_adu))
    allocate(c%DNL(min_adu:max_adu))
    allocate(c%invF(c%min_coadd:c%max_coadd))

    ! Count number of free parameters
    c%npar_adc = 63              ! Global repeated 64-code pattern
    c%npar_adc = c%npar_adc + 64 ! 64 individual codes after midpoint
    if (mod(c%min_adu,64) == 0) then ! Start counting at the first 64-code boundary
       c0 = c%min_adu
    else
       c0 = c%min_adu - mod(c%min_adu,64) + 64
    end if
    do i = c0, c%max_adu, 64
       if (i == 32768) cycle
       c%npar_adc = c%npar_adc + 1 ! Jump between 64-code patterns
    end do
    if (c%myid == 0) write(*,*) '  Number of ADC parameters = ', c%npar_adc, trim(label)
    
    ! Initialize parameter structure
    allocate(c%param_adc(c%npar_adc,3), c%p(c%npar_adc))
    do i = 1, 63
       c%param_adc(i,:) = [i,1,-64]  ! Global repeated 64-code pattern
    end do
    do i = 0, 63
       c%param_adc(64+i,:) = [32768+i,1,1]  ! 64 individual codes after midpoint
    end do
    j = 1
     if (mod(c%min_adu,64) == 0) then ! Start counting at the first 64-code boundary
       c0 = c%min_adu
    else
       c0 = c%min_adu - mod(c%min_adu,64) + 64
    end if
    do i = c0, c%max_adu, 64
       if (i == 32768) cycle
       c%param_adc(127+j,:) = [i,1,1]  ! Jump between 64-code patterns 
       j            = j+1
    end do
    
    ! Initialize bin widths
    c%p     = 1.
    !c%p(64) = 0.2
    call c%param2Q(c%p)

    ! Testing
    if (.false.) then
       if (c%myid == 0) then
          !c%Q(32768) = 0.2
          call c%Q2As(92.)

          open(58, file='adc_fast_'//trim(label)//'.dat', recl=1024)
          write(58,*) '# ADU        Q          v_min         DNL         INL'
          do i = min_adu, max_adu
             write(58,*) i, c%Q(i), c%v(i), c%DNL(i), c%INL(i)
          end do
          close(58)

          open(58, file='adc_As_'//trim(label)//'.dat', recl=1024)
          write(58,*) '# v_in       v_out'
          do i = 1, size(c%A_s%x)
             write(58,*) c%A_s%x(i), c%A_s%y(i)
          end do
          close(58)

          call c%As2F
          
          open(58, file='adc_F_'//trim(label)//'.dat', recl=1024)
          do i = c%min_coadd, c%max_coadd
             delta  = c%invF(int(splint(c%F, real(i,dp))))-i
             write(58,*) i, c%F%y(i-c%min_coadd+1), c%invF(i), delta
          end do
          close(58)          
          
       end if
       call mpi_finalize(ierr)
       stop
    end if
    
  end function constructor_adc_binfit

  module subroutine adc_correct(self, scan, det, tod, sigma0, mask)
    !=========================================================================
    ! Adc corrects a timestream 
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! tod : float array
    !    The tod that is to be corrected
    ! sigma0 : real
    !    White noise rms, coadded
    ! 
    ! Outputs : 
    !
    ! tod_out : float array
    !    The adc corrected version of tod_in    
    ! ====================================================================
    implicit none
    class(comm_adc_binfit),          intent(inout) :: self
    integer(i4b),                    intent(in)    :: scan, det
    real(sp),     dimension(:),      intent(inout) :: tod
    real(sp),                        intent(in)    :: sigma0
    logical(lgt), dimension(:),      intent(in), optional :: mask

    integer(i4b) :: i
    real(dp)     :: delta
    character(len=6) :: scan_text
    character(len=1) :: det_text
    
    ! Initialize transfer function
    !call self%Q2As(sigma0)
    call self%Q2As(92.)
    call self%As2F

!!$    call int2string(scan, scan_text)
!!$    call int2string(det, det_text)
!!$    open(58, file='adc_F_'//trim(self%label)//'_'//scan_text//'_'//det_text//'.dat', recl=1024)
!!$    do i = self%min_coadd, self%max_coadd
!!$       delta  = self%invF(int(splint(self%F, real(i,dp))))-i
!!$       write(58,*) i, self%F%y(i-self%min_coadd+1), self%invF(i), delta
!!$    end do
!!$    close(58)          

    ! Apply correction
    if (present(mask)) then
       do i = 1, size(tod)
          if (mask(i)) then
             tod(i) = self%invF(nint(tod(i)))
          end if
       end do
    else
       do i = 1, size(tod)
          tod(i) = self%invF(nint(tod(i)))
       end do
    end if
  end subroutine adc_correct

  module subroutine param2Q(self, p)
    !=========================================================================
    ! 
    ! 
    ! Inputs:
    !
    ! self:     comm_adc object
    !           Base ADC object
    ! p:        real(sp)
    !           Parameter vector; must be of length npar_adc
    implicit none
    class(comm_adc_binfit),             intent(inout) :: self
    real(sp), dimension(self%npar_adc), intent(in), optional    :: p

    integer(i4b) :: i, m, w, c, c0, c1, b
    
    do i = 1, self%npar_adc
       if (self%param_adc(i,3) < 0) then
          ! Global parameter; affects all indices with modulo m
          c0 =  self%param_adc(i,1)
          w  =  self%param_adc(i,2)
          m  = -self%param_adc(i,3)
          do b = 0, w-1
             if (mod(self%min_adu,m) == 0) then
                c1 = self%min_adu + c0 + b
                !if (self%myid == 0) write(*,*) 'a', i, c0, w, m, b, c1, self%min_adu
             else
                c1 = self%min_adu - mod(self%min_adu,m) + c0 + b
                !if (self%myid == 0) write(*,*) 'b', i, c0, w, m, b, c1, self%min_adu
             end if
             do c = c1, self%max_adu, m
                if (c < self%min_adu) cycle
                if (present(p)) then
                   self%Q(c) = p(i)
                else
                   self%Q(c) = self%p(i)
                end if
             end do
          end do
       else
          ! Local parameter; affects only specified codes
          c0                = self%param_adc(i,1)
          w                 = self%param_adc(i,2)
          if (present(p)) then
             self%Q(c0:c0+w-1) = p(i)
          else
             self%Q(c0:c0+w-1) = self%p(i)
          end if
       end if
    end do
!!$    do i = self%min_adu, self%max_adu
!!$       if (self%myid == 0) write(*,*) i, self%Q(i)
!!$    end do

    ! Normalize widths to unit average
    self%Q = self%Q/sum(self%Q) * size(self%Q) 
    
  end subroutine param2Q
  
  module subroutine Q2As(self, sigma0)
    !=========================================================================
    ! Compute noise weighted transfer function from bin widths
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! sigma0 : real
    !    White noise level (coadded)
    ! 
    ! Outputs : 
    ! ====================================================================
    implicit none
    class(comm_adc_binfit),          intent(inout) :: self
    real(sp),                        intent(in)    :: sigma0

    integer(i4b) :: i, j, j_min, j_max, n
    real(dp)     :: Phi1, Phi2, sigma_sum, s, sigma_limit, s_min, s_max, ds, w, w_tot, sig0
    real(dp), allocatable, dimension(:) :: x, y

    sig0 = sigma0/sqrt(40.)
    
    ! Compute bin boundaries, v
    self%v(self%min_adu) = self%min_adu ! Assume uniform bins at lower values
    do i = self%min_adu+1, self%max_adu+1
       self%v(i) = self%v(i-1) + self%Q(i-1)
    end do

    ! Compute differential nonlinearity, DNL = Q(k)/Q-1
    do i = self%min_adu, self%max_adu
       self%DNL(i) = self%Q(i) - 1.
    end do

    ! Compute integrated nonlinearity (INL = Q(k)/Q -1)
    self%INL(self%min_adu) = self%Q(self%min_adu) - 1.
    do i = self%min_adu+1, self%max_adu
       self%INL(i) = self%INL(i-1) + self%Q(i) - 1.
    end do

    ! Compute noise-weighted transfer function
    ! A_s = int N(x) A(s+x) dx, u = x+s
    !     = int N(u-s) A(u) du
    !     = sum_i int_v(i)^v(i+1) N(u-s) A(u) du
    !     = sum_i i * int_v(i)^v(i+1) N(u-s) du, x = u-s
    !     = sum_i i * int_v(i)-s^v(i+1)-s N(x) dx
    !     = sum_i i * (erf(v(i+1)-s)-erf(i)-s)
    n  = (self%max_adu-self%min_adu) * A_s_oversampling + 1
    ds = (self%max_adu-self%min_adu) / real(n-1,dp)
    allocate(x(n), y(n))
    do i = 1, n
       s     = self%min_adu + ds*(i-1)
       s_min = s - A_s_sigma*sig0
       s_max = s + A_s_sigma*sig0
       x(i)  = s
       
       ! Determine integration range
       j_min = max(locate(self%v, real(s_min,sp))+self%min_adu-1,self%min_adu)
       j_max = j_min
       do while (s_max > self%v(j_max) .and. j_max <= self%max_adu)
          j_max = j_max+1
       end do

       ! Perform integral
       y(i)  = 0.d0
       w_tot = 0.d0
       Phi2  = 0.5d0*(1.d0+corr_erf(real((self%v(j_min)-s)/sig0,dp)/sqrt(2.d0)))
!!$       if (self%myid == 0 .and. abs(s-32768) < 1e-2) then
!!$          write(*,*) 's = ', s, s_min, s_max
!!$          write(*,*) 'range = ', j_min, j_max
!!$          write(*,*) 'v = ', self%v(j_min), self%v(j_max)
!!$       end if
       do j = j_min, j_max-1
          Phi1  = Phi2
          Phi2  = 0.5d0*(1.d0+corr_erf(real((self%v(j+1)-s)/sig0,dp)/sqrt(2.d0)))
          w     = Phi2-Phi1
          y(i)  = y(i) + j * w
          w_tot = w_tot + w
!!$          if (self%myid == 0 .and. abs(s-32768) < 1e-2) then
!!$             write(*,*) j, w, Phi1, Phi2, y(i)
!!$          end if
       end do
       y(i) = y(i) / w_tot ! Guard against truncation and round-off errors
    end do
    call spline_simple(self%A_s, x, y, linear=.true.)
    deallocate(x,y)

  end subroutine Q2As

  module subroutine As2F(self, gain, offset)
    !=========================================================================
    ! Compute coadded transfer function given (optional) fast gain and offset
    ! 
    ! Inputs:
    !
    ! self : comm_adc object
    !    Defines the adc correction that should be applied
    ! gain : real
    !    40-sample vector with gain; should have unity mean
    ! offset: real
    !    40-sample vector with offset in ADU; for HFI, this is typically 4K lines at >180 Hz
    ! 
    ! Outputs : 
    ! ====================================================================
    implicit none
    class(comm_adc_binfit),  intent(inout) :: self
    real(sp), dimension(40), intent(in), optional :: gain
    real(sp), dimension(40), intent(in), optional :: offset

    integer(i4b) :: i, j, j_min, j_max, n
    real(dp)     :: s, s0
    real(dp), allocatable, dimension(:) :: x, y    

    n = self%max_coadd-self%min_coadd+1
    allocate(x(n), y(n))
    do i = 1, n
       x(i) = self%min_coadd + i-1
    end do
    y = 0.d0
    
    if (present(gain) .and. present(offset)) then
       ! Perform full sum with both gain and offset
       do i = 1, n  
          s0 = (self%min_coadd+i-1)/40.d0
          do j = 1, 40
             s    = gain(j)*s0 + offset(j)
             y(i) = y(i) + splint_simple(self%A_s, s)
          end do
       end do       
    else if (present(gain)) then
       ! Perform sum with gain
       do i = 1, n  
          s0 = (self%min_coadd+i-1)/40.d0
          do j = 1, 40
             s    = gain(j)*s0 
             y(i) = y(i) + splint_simple(self%A_s, s)
          end do
       end do       
    else if (present(offset)) then
       ! Perform sum with offset
       do i = 1, n  
          s0 = (self%min_coadd+i-1)/40.d0
          do j = 1, 40
             s    = s0 + offset(j)
             y(i) = y(i) + splint_simple(self%A_s, s)
          end do
       end do       
    else
       ! Assume identical response for all fast samples
       do i = 1, n  
          s0   = (self%min_coadd+i-1)/40.d0
          y(i) = 40.*splint_simple(self%A_s, s0)
       end do
    end if
    call spline_simple(self%F, x, y, linear=.true.)
    
    ! Compute inverse transfer function
    do i = self%min_coadd, self%max_coadd
       j = locate(y, real(i,dp))
       if (j < 1) then
          self%invF(i) = i+0.5d0
       else if (j >= n) then
          self%invF(i) = i + x(n)-y(n) + 0.5d0 ! Uniform bins after available range
       else
          self%invF(i) = zriddr(self%F, x(j), x(j+1), real(i,dp), 1d-5)+0.5d0
       end if
    end do
    
  end subroutine As2F

  module subroutine mcmc_sample_adc(self, handle, chisq_adc)
    !=========================================================================
    ! Sample bin widths for single scan; coadd into posterior
    ! 
    ! Inputs:
    !
    ! self:     comm_adc object
    !           Base ADC object
    implicit none
    class(comm_adc_binfit),          intent(inout) :: self
    type(planck_rng),                intent(inout) :: handle
    interface
       function chisq_adc(p, ndof)
         use healpix_types
         implicit none
         real(sp), dimension(:), intent(in),  optional :: p
         integer(i8b),           intent(out), optional :: ndof
         real(dp)                                      :: chisq_adc
       end function chisq_adc
    end interface
    
    integer(i4b) :: i, j, k, n_gibbs, npar, blocksize, b0, b1
    integer(i8b) :: ndof
    logical(lgt) :: accept
    real(sp)     :: eta, chisq0, chisq_old, chisq_prop, a, mu, sigma
    integer(i4b), allocatable, dimension(:)   :: n_accept
    real(sp),     allocatable, dimension(:)   :: x_prop, rms_prop
    real(sp),     allocatable, dimension(:,:) :: x

    npar      = self%npar_adc
    n_gibbs   = 3
    blocksize = 5

    allocate(x(0:n_gibbs,npar), x_prop(npar), rms_prop(npar), n_accept(npar))
    rms_prop = 0.01
    n_accept = 0

    x(0,:)    = self%p
    chisq_old = chisq_adc(x(0,:), ndof)
    chisq0    = chisq_old
    open(58,file='adc_mcmc.dat', recl=8192)
    do i = 1, n_gibbs
       x_prop  = x(i-1,:)
       do j = 1, npar, blocksize
          b0 = j
          b1 = min(j+blocksize-1,npar)
          
          ! Draw proposal
          do k = b0, b1
             x_prop(k) = x_prop(k) + rand_gauss(handle) * rms_prop(k)
          end do

          if (any(x_prop(b0:b1) < 0.)) then
             x_prop(b0:b1) = x(i-1,b0:b1) ! Reset to previous sample
             cycle
          end if
          
          ! Compute chisq, and apply Metropolis rule
          chisq_prop = chisq_adc(x_prop)
          if (chisq_prop < chisq_old) then
             accept = .true.
          else
             a      = exp(-0.5d0*(chisq_prop-chisq_old))
             accept = rand_uni(handle) < a
          end if

          write(58,*) i, j, accept, chisq_old/ndof, chisq_prop/ndof, x_prop
          write(*,fmt='(i4,i4,l,4f8.3)')  i, j, accept, chisq_old/ndof, chisq_prop/ndof, x(i-1,b0), x_prop(b0)
          
          if (accept) then
             chisq_old       = chisq_prop
             n_accept(b0:b1) = n_accept(b0:b1) + 1
          else
             x_prop(b0:b1)   = x(i-1,b0:b1) ! Reset to previous sample
          end if
       end do
       x(i,:) = x_prop
    end do
    close(58)

    ! Update main parameter
    self%p = x(n_gibbs,:)
    
    write(*,fmt='(a,f8.3,a,f8.3,a,f8.3)') ' adc mcmc -- X^2 old = ', chisq0/ndof, ', new =', chisq_old/ndof, ', mean accept = ', mean(n_accept/real(n_gibbs,dp))

    ! Clean up
    deallocate(x, x_prop, rms_prop, n_accept)
  end subroutine mcmc_sample_adc

end submodule comm_tod_adc_binfit_smod
