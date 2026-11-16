!================================================================================
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
!================================================================================
module comm_tod_cray_mod
  use comm_param_mod
  use spline_1D_mod
  implicit none
  
  private
  public constructor_cray, comm_tod_cray, cray_ptr
  
  integer(i4b), parameter :: CR_CHUNK_SIZE = 100000  ! Number of CRs to add in each array reallocstiom
  
  ! Cosmic ray types
  integer(i4b), parameter :: CR_SHORT   =  1
  integer(i4b), parameter :: CR_LONG    =  2
  integer(i4b), parameter :: CR_SLOW    =  3
  integer(i4b), parameter :: CR_BRIGHT  = -1
  
  ! Class for single cosmic ray event
  type :: cray_event
     integer(i4b) :: type    ! CR type, defined above
     integer(i4b) :: nsamp   ! Maximum TOD segment length
     real(dp)     :: fsamp   ! Sampling rate in Hz
     integer(i4b) :: nspline ! Number of baseline spline nodes 
     integer(i4b) :: mask(2) ! Start and end sample of mask
     real(dp), allocatable, dimension(:) :: p_cr   ! Free CR shape parameters
     real(dp), allocatable, dimension(:) :: p_base ! Free baseline parameters
     real(sp), allocatable, dimension(:) :: T_cr   ! CR shape template; normalized to 1
     real(sp), allocatable, dimension(:) :: T_base ! Baseline template; absolute normalization
   contains
     procedure :: get_baseline_spline_nodes
     procedure :: build_cr_template
     procedure :: build_baseline_template
     procedure :: fit_bright_event_with_baseline
     procedure :: dealloc => deallocate_event
  end type cray_event
  
  type :: comm_tod_cray
     integer(i4b)      :: n        ! Number of events in current scan
     integer(i4b)      :: nmax     ! Maximum number of events in current arrays
     integer(i4b)      :: ntypes   ! Number of active base types
     character(len=64) :: freq     ! Frequency label
     integer(i4b)      :: det      ! Detector ID
     integer(i4b)      :: scanid   ! Absolute scan ID
     integer(i4b)      :: mod_phase ! 1 => odd = +1; -1 => even = +1
     real(dp)          :: fsamp   ! Sampling rate in Hz
     logical(lgt)      :: first_call  ! Set to false after detecting faint CRs
     type(cray_event), pointer,     dimension(:) :: T_base       ! Default templates
     type(cray_ptr),   allocatable, dimension(:) :: event_list   ! List of individual event objects; 
     real(sp),         allocatable, dimension(:) :: amp          ! CR amplitudes
     integer(i4b),     allocatable, dimension(:) :: event_start  ! Starting sample number
     integer(i4b),     allocatable, dimension(:) :: event_length ! TOD truncation length of each event
   contains
     procedure :: detect_bright_events
     procedure :: detect_events
     procedure :: add_event
     procedure :: select_cr_type
     procedure :: fit_local_shape_params
     procedure :: fit_amps
     procedure :: generate => generate_cray_correction
     procedure :: params2cache                            ! IO communication routine
     procedure :: set_mod_phase
  end type comm_tod_cray
  
  interface comm_tod_cray 
     procedure constructor_cray
  end interface comm_tod_cray
  
  interface cray_event 
     procedure constructor_cray_event
  end interface cray_event
  
  type cray_ptr
     class(cray_event), pointer :: p => null()
  end type cray_ptr
  
contains

  ! Constructor
  function constructor_cray(freq, det, scanid, active_cr_types, fsamp) result(c)
    implicit none
    character(len=*),           intent(in) :: freq             ! Frequency label
    integer(i4b),               intent(in) :: det              ! Detector ID
    integer(i4b),               intent(in) :: scanid             ! Scan ID (absolute)
    integer(i4b), dimension(:), intent(in) :: active_cr_types  ! List of active types
    real(dp),                   intent(in) :: fsamp            ! Sampling rate in Hz
    class(comm_tod_cray), pointer          :: c
    
    allocate(c)
    c%freq       = freq
    c%det        = det
    c%scanid     = scanid
    c%ntypes     = size(active_cr_types)
    c%fsamp      = fsamp
    c%mod_phase   = -1000000
    c%first_call = .true.
    
    ! Initialize arrays; events should br sorted according to event_start
    c%n    = 0
    c%nmax = CR_CHUNK_SIZE
    allocate(c%event_list(c%nmax))
    allocate(c%amp(c%nmax))
    allocate(c%event_start(c%nmax))
    allocate(c%event_length(c%nmax))
    
  end function constructor_cray
  
  ! Constructor for individual events
  function constructor_cray_event(type, nsamp, fsamp, p_cr, p_base) result(c)
    implicit none
    integer(i4b),                    intent(in)           :: type   ! CR type ID 
    integer(i4b),                    intent(in)           :: nsamp  ! Number of samples
    real(dp),                        intent(in)           :: fsamp  ! Sampling rate in Hz
    real(dp),          dimension(:), intent(in)           :: p_cr   ! CR parameters
    real(dp),          dimension(:), intent(in), optional :: p_base ! Baseline amplitudes
    class(cray_event), pointer                            :: c
    
    allocate(c)
    c%type       = type
    c%nsamp      = nsamp
    c%fsamp      = fsamp
    c%mask       = -1     ! Initalize without mask
    
    ! Initialize CR template
    allocate(c%p_cr(size(p_cr)), c%T_cr(nsamp))
    c%p_cr = p_cr
    call c%build_cr_template()
    
    ! Initialize baseline template
    if (c%type < 0) then
       c%nspline = size(p_base)
       allocate(c%p_base(c%nspline), c%T_base(nsamp))
       c%p_base = p_base
       call c%build_baseline_template(p_cr(5))
    end if
    
  end function constructor_cray_event

  subroutine deallocate_event(self)
    implicit none
    class(cray_event), intent(inout) :: self

    if (allocated(self%p_cr))   deallocate(self%p_cr)
    if (allocated(self%p_base)) deallocate(self%p_base)
    if (allocated(self%T_cr))   deallocate(self%T_cr)
    if (allocated(self%T_base)) deallocate(self%T_base)
        
  end subroutine deallocate_event
  
  ! Routine for constructing TOD shape template for current event/base model
  subroutine build_cr_template(self)
    implicit none
    class(cray_event), intent(inout) :: self
    
    integer(i4b) :: i
    real(dp)     :: t, dt, t_dep
    real(dp)     :: tau0 = 0.002d0
    
    t     = 0.d0
    dt    = 1.d0/self%fsamp
    t_dep = self%p_cr(5)
    do i = 1, self%nsamp
       if (t < t_dep) then
          self%T_cr(i) = 0.d0
       else
          self%T_cr(i) =  self%p_cr(1) * (exp(-(t-t_dep)/self%p_cr(3)) - exp(-(t-t_dep)/tau0)) + &
                        & self%p_cr(2) * (exp(-(t-t_dep)/self%p_cr(4)) - exp(-(t-t_dep)/tau0))
       end if
       t = t + dt
    end do
    
  end subroutine build_cr_template
  
  ! Routine for constructing TOD shape template for current event/base model
  subroutine build_baseline_template(self, t_dep)
    implicit none
    class(cray_event),              intent(inout) :: self
    real(dp),                       intent(in)    :: t_dep
    
    integer(i4b) :: i
    real(dp)     :: t, dt
    real(dp), allocatable, dimension(:) :: x, y
    type(spline_type) :: b
    
    call self%get_baseline_spline_nodes(x) 
    
    ! Subtract glitch baseline
    self%nspline = size(x)-2
    allocate(y(self%nspline+2))
    y(1:self%nspline)  = self%p_base
    y(self%nspline+1:) = 0.d0
    call spline(b, x, y)
    t = 0.d0; dt = 1.d0/self%fsamp
    do i = 1, self%nsamp
       if (t < t_dep .or. t > t_dep+x(self%nspline)) then
          self%T_base(i) = 0.d0
       else
          self%T_base(i) = splint(b, t-t_dep)
       end if
       t = t + dt
    end do
    
    deallocate(x, y)
    call free_spline(b)
    
  end subroutine build_baseline_template

  ! builds a set of templates for cosmic ray response from the TOD
  function nsamp2nspline(nsamp) result(nspline)
    implicit none
    integer(i4b), intent(in)  :: nsamp
    integer(i4b)              :: nspline
    nspline = 3                                                 ! Step size = 3
    if (nsamp > 10) nspline = nspline + min(nsamp-10, 20)/5 + 1 ! Step size = 5
    if (nsamp > 30) nspline = nspline + (nsamp-30)/10 + 1       ! Step size = 10
  end function nsamp2nspline

  
  ! builds a set of templates for cosmic ray response from the TOD
  subroutine get_baseline_spline_nodes(self, x)
    implicit none
    class(cray_event),                            intent(inout) :: self
    real(dp),          allocatable, dimension(:), intent(out)   :: x
    
    integer(i4b) :: j, k
    real(dp)     :: dt
    
    dt = 1.d0/self%fsamp
        
    ! Initialize spline nodes
    allocate(x(self%nspline+2)) ! Allow two dead nodes to go smoothly to zero
    x(1) = 0.d0
    x(2) = 3.d0*dt
    x(3) = 6.d0*dt
    k = 10; j = 4
    do while (k < 30)
       x(j) = k*dt
       k    = k+5
       j    = j+1
    end do
    do while (j <= self%nspline+2)
       x(j) = k*dt
       k    = k+10
       j    = j+1
    end do
    
  end subroutine get_baseline_spline_nodes
  
  ! builds a set of templates for cosmic ray response from the TOD
  subroutine set_mod_phase(self, mod_phase)
    implicit none
    class(comm_tod_cray), intent(inout) :: self
    integer(i4b),         intent(in)    :: mod_phase
    self%mod_phase = mod_phase
  end subroutine set_mod_phase
  
  ! Search for bright events that changes the baseline
  ! This routine should only be called once per main run
  subroutine detect_bright_events(self, tod, s_tot, sigma)
    implicit none
    class(comm_tod_cray),   intent(inout) :: self
    real(sp), dimension(:), intent(in)    :: tod   ! Raw modulated TOD
    real(sp), dimension(:), intent(in)    :: s_tot ! Total signal in ADU
    real(sp),               intent(in)    :: sigma ! White noise rms in ADU
    
    integer(i4b) :: i, j, len, min_len, n_base, ntod, npar, nspline
    real(dp)     :: threshold, b1, b2, mu1, mu2, rms, rms2, chisq_corr, chisq_uncorr
    real(dp), allocatable, dimension(:) :: res, p
    class(cray_event), pointer :: event
    
    if (self%mod_phase == -1000000) then
       write(*,*) 'comm_tod_cray_mod: mod_phase not set'
       stop
    end if
    
    threshold = 30.*sigma  ! Detection trigger 
    min_len   = 10   ! Minimum number of samples for a bright event
    n_base    = 10000 ! Length for median-based baseline
    ntod      = size(tod)
    
    ! Subtract baselines
    allocate(res(ntod))
    b1 = median(tod(1:n_base:2))
    b2 = median(tod(2:n_base:2))
    res(1:ntod:2) = tod(1:ntod:2) - b1
    res(2:ntod:2) = tod(2:ntod:2) - b2
    
    ! Subtract modulated signal model
    do i = 1, ntod
       if (mod(i,2) == 1) then
          res(i) = res(i) - self%mod_phase*s_tot(i)
       else
          res(i) = res(i) + self%mod_phase*s_tot(i)
       end if
    end do
    
    ! Search for bright CR candidates
    i = 0
    do while (i < ntod-min_len)
       i = i+1
       if (abs(res(i)) < threshold) cycle
       ! Sample i is a candidate starting sample
       
       ! 1) Find distance to local rms has returned to normal levels
       j   = i+min_len-1
       rms = sqrt(variance(res(j-min_len+1:j)))
       do while (rms > 2.d0*sigma .and. j < ntod-min_len+1)
          j   = j + min_len
          rms = sqrt(variance(res(j-min_len+1:j)))
       end do
       
       ! 2) Check that odd and even samples move together after instability period
!!$       rms  = sqrt(variance(res(i+n_mask:j:2)-res(i+n_mask+1:j:2))) ! Odd - even
!!$       rms2 = sqrt(variance(res(i+n_mask:j)))                       ! Full
!!$
!!$       write(*,*) 'a', i, j, rms, rms2
!!$       
!!$       ! 3) Reject if too small difference; glitch should be handled elsewhere
!!$       if (rms > 0.2d0 * rms2) then
!!$          i = j ! Skip forward to stable range
!!$          cycle
!!$       end if

       ! 4) Create a bright cray event object
       nspline = nsamp2nspline(j-i+4)
       npar    = nspline+5
       allocate(p(npar))
       p(1:nspline) = 0.d0         ! Spline amps
       p(nspline+1) = 1.d0         ! a1 CR
       p(nspline+2) = 1.d0         ! a2 CR
       p(nspline+3) = 0.040        ! tau2 CR
       p(nspline+4) = 0.020        ! tau4 CR
       p(nspline+5) = 3/self%fsamp ! t_dep CR
       event => cray_event(-1, j-i+4, self%fsamp, p(nspline+1:), p(1:nspline))

       ! 5) Fit CR + baseline model; measure chisq improvement
       call event%fit_bright_event_with_baseline(res(i-3:j), sigma, i-3, &
            & self%mod_phase, chisq_corr, chisq_uncorr)

       ! 6) Accept if the chisq improvement is large (10 times Akaike IC measure)
       !    Start of current event is at i-3; end is at j
       if (chisq_corr < chisq_uncorr) then
          write(*,fmt='(a,2i10,2f16.5)') 'Bright CR = ', i, j-i+4, chisq_corr, chisq_uncorr
          call self%add_event(i-3, j-i+4, 1.0, event)
          i = j ! Skip forward
       else
          call event%dealloc()
       end if
       !stop
       
       deallocate(p)
    end do

    write(*,*) 'Bright cray -- ', trim(self%freq), self%det, self%scanid, ', n = ', self%n
    
  end subroutine detect_bright_events
  
  ! Fit bright CRs + baseline model
  subroutine fit_bright_event_with_baseline(self, res, sigma, i0, mod_phase, chisq_corr, chisq_uncorr)
    implicit none
    class(cray_event),                   intent(inout) :: self
    real(dp),              dimension(:), intent(in)    :: res
    real(sp),                            intent(in)    :: sigma
    integer(i4b),                        intent(in)    :: i0        ! First sample
    integer(i4b),                        intent(in)    :: mod_phase
    real(dp),                            intent(out)   :: chisq_corr, chisq_uncorr
    
    integer(i4b) :: i, j, k, m, err, npar, ind(1)
    real(dp)     :: fsamp, dt, threshold
    logical(lgt) :: incr_low, incr_high
    real(dp), allocatable, dimension(:) :: p, p0, corr

    threshold = 5.d0  ! mask threshold in sigma 
    m = size(res) ! Number of slow samples in template
        ! Initialize parameters
    npar    = 5 + self%nspline
    allocate(p(npar), p0(npar), corr(m))
    p(self%nspline+1:) = self%p_cr
    p(1:self%nspline)  = self%p_base
    
    ! Perform fit
    call powell(p, chisq_cr, err)

    ! Check residual; increase mask if necessary
    call correct_cr(res, p, corr)
    ind        = maxloc(corr**2)
    chisq_corr = 0.d0
    !write(*,*) 'a', corr(ind(1)), sigma, corr(ind(1))**2 / sigma**2
    if (corr(ind(1))**2 / sigma**2 > threshold) then
       j = 3; k = ind(1)
       do while (corr(j)**2 / sigma**2 > threshold .or. &
               & corr(k)**2 / sigma**2 > threshold .or. &
               & chisq_corr > 2)
          incr_low  = corr(j)**2 / sigma**2 > threshold .or. chisq_corr > 2
          incr_high = corr(k)**2 / sigma**2 > threshold .or. chisq_corr > 2
          if (incr_low)       j = j-2
          if (incr_high)      k = k+2
          if (j <= 0)         j = 1
          if (k > self%nsamp) k = self%nsamp
          if ((.not. incr_low .or. j == 1) .and. (.not. incr_high .or. k == self%nsamp)) then
             self%mask = i0 + [j-1,k-1]
             exit
          else
             self%mask = i0 + [j+1,k-1] - 1
             call powell(p, chisq_cr, err)
             chisq_corr = chisq_cr(p)
             call correct_cr(res, p, corr)
             !write(*,*) 'b0', corr
             !write(*,*) 'b', corr(j), corr(k), chisq_corr, self%mask
          end if
       end do
       ! Increase mask by two samples in each direction for good measure
       self%mask(1) = max(self%mask(1)-2,1)
       self%mask(2) = min(self%mask(2)+2,i0+self%nsamp)
    end if
    
    ! Chisq without any correction
    p0 = p
    p0(1:self%nspline+2)  = 0.d0
    chisq_uncorr = chisq_cr(p0)

    ! Chisq for best-fit model
    chisq_corr = chisq_cr(p)
        
    deallocate(p, p0)
    
  contains
    
    subroutine correct_cr(tod, p, corr)
      implicit none
      real(dp), dimension(:), intent(in)  :: tod
      real(dp), dimension(:), intent(in)  :: p
      real(dp), dimension(:), intent(out) :: corr
      
      real(dp)          :: s, t, t_dep
      integer(i4b)      :: i
      type(spline_type) :: b
      
      corr  = tod
      t_dep = p(self%nspline+5)
      
      ! Subtract glitch baseline
      self%p_cr   = p(self%nspline+1:)
      self%p_base = p(1:self%nspline)
      call self%build_cr_template()
      call self%build_baseline_template(t_dep)
    
      ! Subtract baseline template
      corr = corr - self%T_base
      
      ! Subtract modulated CR model
      do i = 1, m
         if (i0+i-1 >= self%mask(1) .and. i0+i-1 <= self%mask(2)) then
            ! Masked sample
            corr(i) = 0.d0
         else
            if (mod(i0+i-1,2) == 1) then
               corr(i) = corr(i) - mod_phase*self%T_cr(i)
            else
               corr(i) = corr(i) + mod_phase*self%T_cr(i)
            end if
         end if
      end do
      
    end subroutine correct_cr
    
    function chisq_cr(p) result(chisq)
      implicit none
      real(dp), dimension(:), intent(in), optional :: p
      real(dp)     :: chisq
      
      integer(i4b) :: i, ndof
      integer(i4b), save :: counter = 0
      real(dp), allocatable, dimension(:) :: r
      
      if (any(p(self%nspline+3:self%nspline+5) < 0.d0) .or. &
           & p(self%nspline+3) < 0.01d0 .or. p(self%nspline+4) < 0.005d0) then
         chisq = 1d30
         return
      end if
      
      allocate(r(m))
      call correct_cr(res, p, r)

      ! Compute chisq
      chisq = 0.d0
      do i = 1, self%nsamp
         if (self%mask(1) == -1 .or. i0+i-1 < self%mask(1) .or. i0+i-1 > self%mask(2)) then
            chisq = chisq + r(i)**2 / sigma**2
         end if
         !write(*,*) i, i0, self%mask(1), i0+i-1, self%mask(2), chisq, r(i), sigma
      end do
      !stop

      ! Compute reduced chisq
      ndof = self%nsamp
      if (self%mask(1) /= -1) ndof = ndof - (self%mask(2)-self%mask(1)+1)
      chisq = chisq / ndof
      
      !counter = counter+1
      !if (mod(counter,10) == 0) write(*,*) counter, chisq, p
      deallocate(r)
      
    end function chisq_cr
    
  end subroutine fit_bright_event_with_baseline

  ! Add new event to main data structure
  subroutine add_event(self, start, len, amp, event)
    implicit none
    class(comm_tod_cray),     intent(inout) :: self
    integer(i4b),             intent(in)    :: start, len
    real(sp),                 intent(in)    :: amp
    type(cray_event), target, intent(in)    :: event
    
    integer(i4b) :: i, j
    
    ! Expand arrays if full
    if (self%n == self%nmax) then
       
    end if
    
    ! Search for new position
    if (self%n == 0) then
       j = 1
    else if (start > self%event_start(self%n)) then
       j = self%n+1
    else
       j = locate(self%event_start(1:self%n), start)+1
    end if
    
    ! Create space for new point
    self%n = self%n + 1
    self%event_start(j+1:self%n)  =  self%event_start(j:self%n-1)
    self%event_length(j+1:self%n) =  self%event_length(j:self%n-1)
    do i = self%n, j+1, -1
       self%event_list(i)%p => self%event_list(i-1)%p
    end do
    self%amp(j+1:self%n)          =  self%amp(j:self%n-1)
    
    ! Add new point
    self%event_start(j)  =  start
    self%event_length(j) =  len
    self%amp(j)          =  amp
    self%event_list(j)%p => event
    
  end subroutine add_event

  
  ! Search for events; store list of impact samples; allocate event arrays;
  ! initialize defaults; define local length
  ! This routine should only be called once per main run
  subroutine detect_events(self, tod)
    implicit none
    class(comm_tod_cray),                      intent(inout) :: self
    real(sp), dimension(:),                    intent(in)    :: tod
    
    integer(i4b) :: i
    
    self%first_call = .false.
  end subroutine detect_events
  
  
  ! Search for events; store list of impact samples; allocate event arrays;
  ! initialize defaults
  subroutine select_cr_type(self, tod)
    implicit none
    class(comm_tod_cray),                      intent(inout) :: self
    real(sp), dimension(:),                    intent(in)    :: tod
    
    integer(i4b) :: i
    
  end subroutine select_cr_type
  
  
    ! fits the constructed templates to the cosmic rays in the timestreams
  subroutine fit_local_shape_params(self, tod)
    implicit none
    class(comm_tod_cray),                      intent(inout) :: self
    real(sp), dimension(:),                    intent(in)    :: tod
    
    integer(i4b) :: i
    
    ! Fit parameters for individual events
    do i = 1, self%n
       if (self%event_list(i)%p%type > 0) cycle  ! Skip elements modelled by global templates
    end do
    
  end subroutine fit_local_shape_params
  
  ! fits the constructed templates to the cosmic rays in the timestreams
  subroutine fit_amps(self, tod)
    implicit none
    class(comm_tod_cray),                      intent(inout) :: self
    real(sp), dimension(:),                    intent(in)    :: tod
    self%amp = 0.
  end subroutine fit_amps
  
  ! Routine for constructing total CR model for a given scan
  subroutine generate_cray_correction(self, s_cray)
    implicit none
    class(comm_tod_cray),               intent(in)  :: self
    real(sp),             dimension(:), intent(out) :: s_cray
    
    integer(i4b) :: i
    
    s_cray = 0.
    do i = 1, self%n
       ! Add contribution for current event to s_cray
    end do
    
  end subroutine generate_cray_correction
  
  ! IO communication routine
  subroutine params2cache(self, cache)
    implicit none
    class(comm_tod_cray),                         intent(in)    :: self
    real(sp),          allocatable, dimension(:), intent(inout) :: cache
    
    if (allocated(cache)) then
       ! Initialize model parameters from cache
    else
       ! Store model parameters in cache
    end if
    
  end subroutine params2cache
  
end module comm_tod_cray_mod
