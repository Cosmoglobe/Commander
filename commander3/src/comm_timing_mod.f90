module comm_timing_mod
  use comm_utils
  implicit none

  ! Global parameters
  integer(i4b), parameter, public :: NUM_GLOBAL    =  10
  integer(i4b), parameter, public :: TOT_RUNTIME   =  1
  integer(i4b), parameter, public :: TOT_INIT      =  2
  integer(i4b), parameter, public :: TOT_FFT       =  3
  integer(i4b), parameter, public :: TOT_SHT       =  4
  integer(i4b), parameter, public :: TOT_GIBBSSAMP =  5
  integer(i4b), parameter, public :: TOT_AMPSAMP   =  6
  integer(i4b), parameter, public :: TOT_TODPROC   =  7
  integer(i4b), parameter, public :: TOT_SPECIND   =  8
  integer(i4b), parameter, public :: TOT_CLS       =  9
  integer(i4b), parameter, public :: TOT_OUTPUT    =  10

  ! Channel specific parameters
  integer(i4b), parameter, public :: NUM_TOD       = 26
  integer(i4b), parameter, public :: TOD_TOT       =  1
  integer(i4b), parameter, public :: TOD_INIT      =  2
  integer(i4b), parameter, public :: TOD_SL_PRE    =  3
  integer(i4b), parameter, public :: TOD_SL_INT    =  4
  integer(i4b), parameter, public :: TOD_PROJECT   =  5
  integer(i4b), parameter, public :: TOD_ORBITAL   =  6
  integer(i4b), parameter, public :: TOD_DECOMP    =  7
  integer(i4b), parameter, public :: TOD_ABSCAL    =  8
  integer(i4b), parameter, public :: TOD_RELCAL    =  9
  integer(i4b), parameter, public :: TOD_DELTAG    = 10
  integer(i4b), parameter, public :: TOD_NCORR     = 11
  integer(i4b), parameter, public :: TOD_XI_N      = 12
  integer(i4b), parameter, public :: TOD_MAPBIN    = 13
  integer(i4b), parameter, public :: TOD_MAPSOLVE  = 14
  integer(i4b), parameter, public :: TOD_ZODI      = 15
  integer(i4b), parameter, public :: TOD_IMBAL     = 16
  integer(i4b), parameter, public :: TOD_1HZ       = 17
  integer(i4b), parameter, public :: TOD_4D        = 18
  integer(i4b), parameter, public :: TOD_CHISQ     = 19
  integer(i4b), parameter, public :: TOD_BP        = 20
  integer(i4b), parameter, public :: TOD_WAIT      = 21
  integer(i4b), parameter, public :: TOD_MPI       = 22
  integer(i4b), parameter, public :: TOD_BASELINE  = 23
  integer(i4b), parameter, public :: TOD_ALLOC     = 24
  integer(i4b), parameter, public :: TOD_INSTCORR  = 25
  integer(i4b), parameter, public :: TOD_WRITE     = 26

  private
  public comm_timing

  type comm_timing
     integer(i4b) :: numband, comm, myid, n_tot
     integer(i4b), allocatable, dimension(:)   :: numsamp
     real(dp),     allocatable, dimension(:)   :: t     ! Accumulated time for global timers (NUM_GLOBAL+NUM_TOD*numband)
     real(dp),     allocatable, dimension(:)   :: t1    ! Start time for currently active timers for global timers (NUM_GLOBAL+NUM_TOD*numband)
   contains
     procedure :: start        => comm_timer_start
     procedure :: stop         => comm_timer_stop
     procedure :: incr_numsamp => comm_timer_incr_numsamp
     procedure :: dumpASCII    => comm_timer_dumpASCII
  end type comm_timing

  interface comm_timing
     procedure constructor
  end interface comm_timing

contains

  function constructor(numband, comm) result(res)
    ! 
    ! Constructor routine for timer object
    ! 
    ! Input variables:
    !    numband  = number of frequency channels in current run
    !    comm     = MPI communicator
    !
    ! Output variable:
    !    res      = pointer to comm_timing object
    ! 
    implicit none
    integer(i4b),               intent(in) :: numband, comm
    type(comm_timing), pointer             :: res
    
    integer(i4b) :: ierr

    allocate(res)
    res%numband = numband
    res%n_tot   = NUM_GLOBAL + NUM_TOD * numband
    res%comm    = comm
    call mpi_comm_rank(comm, res%myid, ierr) 

    allocate(res%t(res%n_tot), res%t1(res%n_tot))
    allocate(res%numsamp(0:numband))
    res%numsamp = 0
    res%t       = 0.d0
    res%t1      = 0.d0

  end function constructor


  subroutine comm_timer_start(self, timer_id, band)
    ! 
    ! Routine for starting timers
    ! 
    ! Input variables:
    !    self     = comm_timing object
    !    timer_id = timer ID
    !    bands    = band ID (optional)
    ! 
    implicit none
    class(comm_timing), intent(inout)          :: self
    integer(i4b),       intent(in)             :: timer_id
    integer(i4b),       intent(in),   optional :: band
    
    integer(i4b) :: i, timer
    real(dp)     :: t1

    call wall_time(t1)
    timer = timer_id; if (present(band)) timer = NUM_GLOBAL + NUM_TOD*(band-1) + timer_id
    self%t1(timer) = t1

  end subroutine comm_timer_start

  subroutine comm_timer_stop(self, timer_id, band)
    ! 
    ! Routine for stopping timers; difference between t2 and t1 will be 
    ! accumulated into self%t_tot. Only currently active timers are accumulated
    ! 
    ! Input variables:
    !    self     = comm_timing object
    !    timer_id = timer ID
    !    bands    = band ID (optional)
    ! 
    implicit none
    class(comm_timing), intent(inout)           :: self
    integer(i4b),       intent(in)              :: timer_id
    integer(i4b),       intent(in),   optional  :: band
    
    integer(i4b) :: i, timer
    real(dp)     :: t2

    timer = timer_id; if (present(band)) timer = NUM_GLOBAL + NUM_TOD*(band-1) + timer_id
    if (self%t1(timer) > 0) then
       call wall_time(t2)
       self%t(timer)  = self%t(timer) + t2 - self%t1(timer)
       self%t1(timer) = 0.d0
    end if

  end subroutine comm_timer_stop

  subroutine comm_timer_incr_numsamp(self, band)
    ! 
    ! Routine for incrementing sample counter; used to output time per sample
    ! 
    ! Input variables:
    !    self     = comm_timing object
    ! 
    implicit none
    class(comm_timing),               intent(inout) :: self
    integer(i4b),       intent(in)   :: band

    self%numsamp(band) = self%numsamp(band) + 1

  end subroutine comm_timer_incr_numsamp


  subroutine comm_timer_dumpASCII(self, labels, filename)
    ! 
    ! Routine for outputting timing information
    ! 
    ! Input variables:
    !    self     = comm_timing object
    !    labels   = frequency channel labels
    !    filename = output filename
    ! 
    implicit none
    class(comm_timing),               intent(inout) :: self
    character(len=*),   dimension(:), intent(in)    :: labels
    character(len=*),                 intent(in)    :: filename

    integer(i4b) :: unit, ierr, band, b
    real(dp), dimension(self%n_tot) :: t

    if (self%numsamp(0) == 0) return

    call self%stop(TOT_RUNTIME)
    call self%start(TOT_RUNTIME)

    call mpi_reduce(self%t, t, self%n_tot, MPI_DOUBLE_PRECISION, &
         & MPI_SUM, 0, self%comm, ierr)
    t = t/3600 ! CPU-hours

    if (self%myid == 0) then
       unit = getlun()
       open(unit,file=trim(filename), recl=1024)
       write(unit,*) 'Timing summary'
       write(unit,*) ''
       write(unit,*) '   Numbers are given in (CPU-hours,%)'
       write(unit,*) ''
       write(unit,*) '   Global total timers:'
       write(unit,*) '      Number of samples             = ', self%numsamp(0)
       write(unit,fmt='(a,f12.3,"h")') '       Total runtime                 = ', t(TOT_RUNTIME)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Initialization                = ', t(TOT_INIT), 100*t(TOT_INIT)/t(TOT_RUNTIME)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Total FFT                     = ', t(TOT_FFT), 100*t(TOT_FFT)/t(TOT_RUNTIME)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Total SHT                     = ', t(TOT_SHT), 100*t(TOT_SHT)/t(TOT_RUNTIME)
       write(unit,*) ''
       write(unit,*) '   Global per-sample timers:'
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Chain output                  = ', t(TOT_OUTPUT)  / self%numsamp(0), 100*t(TOT_OUTPUT)/t(TOT_GIBBSSAMP)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Amplitude sampling            = ', t(TOT_AMPSAMP) / self%numsamp(0), 100*t(TOT_AMPSAMP)/t(TOT_GIBBSSAMP)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Spectral index sampling       = ', t(TOT_SPECIND) / self%numsamp(0), 100*t(TOT_SPECIND)/t(TOT_GIBBSSAMP)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Cls sampling                  = ', t(TOT_CLS)     / self%numsamp(0), 100*t(TOT_CLS)/t(TOT_GIBBSSAMP)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       TOD processing                = ', t(TOT_TODPROC) / self%numsamp(0), 100*t(TOT_TODPROC)/T(TOT_GIBBSSAMP)
       write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '       Other                         = ', (t(TOT_GIBBSSAMP)-sum(t(6:10))) / self%numsamp(0), 100*(t(TOT_GIBBSSAMP)-sum(t(6:10)))/t(TOT_GIBBSSAMP)
          write(unit,fmt='(a,f12.3,"h")')        '       Total cost per Gibbs sample   = ', t(TOT_GIBBSSAMP)     / self%numsamp(0)
       write(unit,*) ''
       write(unit,*) '   Channel-specific global timers:'

       write(*, *) "numsamp:", self%numsamp

       do band = 1, self%numband
          b = NUM_GLOBAL + (band-1)*NUM_TOD
          if (all(t(b+1:b+NUM_TOD) == 0.d0) .or. T(b+TOD_TOT) == 0.d0) cycle
          write(unit,*) 
          write(unit,*) '     Channel                      = ', trim(labels(band))
          write(unit,fmt='(a,f12.3,"h")') '      TOD initialization           = ', t(b+TOD_INIT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD sidelobe precomputation  = ', t(b+TOD_SL_PRE)   / self%numsamp(band), 100*t(b+TOD_SL_PRE)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD sidelobe interpolation   = ', t(b+TOD_SL_INT)   / self%numsamp(band), 100*t(b+TOD_SL_INT)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD sky-to-tod projection    = ', t(b+TOD_PROJECT)  / self%numsamp(band), 100*t(b+TOD_PROJECT)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD orbital dipole           = ', t(b+TOD_ORBITAL)  / self%numsamp(band), 100*t(b+TOD_ORBITAL)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD decompression            = ', t(b+TOD_DECOMP)   / self%numsamp(band), 100*t(b+TOD_DECOMP)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD absolute calibration     = ', t(b+TOD_ABSCAL)   / self%numsamp(band), 100*t(b+TOD_ABSCAL)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD relative calibration     = ', t(b+TOD_RELCAL)   / self%numsamp(band), 100*t(b+TOD_RELCAL)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD delta G calibration      = ', t(b+TOD_DELTAG)   / self%numsamp(band), 100*t(b+TOD_DELTAG)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD transmission imbalance   = ', t(b+TOD_IMBAL)    / self%numsamp(band), 100*t(b+TOD_IMBAL)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD baseline sampling        = ', t(b+TOD_BASELINE)  /  self%numsamp(band),   100*t(b+TOD_BASELINE)/t(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD 1 Hz spikes              = ', t(b+TOD_1HZ)    / self%numsamp(band), 100*t(b+TOD_1HZ)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD correlated noise         = ', t(b+TOD_NCORR)    / self%numsamp(band), 100*t(b+TOD_NCORR)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD corr noise PSD           = ', t(b+TOD_XI_N)     / self%numsamp(band), 100*t(b+TOD_XI_N)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD binning                  = ', t(b+TOD_MAPBIN)   / self%numsamp(band), 100*t(b+TOD_MAPBIN)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD chisq                    = ', t(b+TOD_CHISQ)   / self%numsamp(band), 100*t(b+TOD_CHISQ)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD bandpass                 = ', t(b+TOD_BP)   / self%numsamp(band), 100*t(b+TOD_BP)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD map solution             = ', t(b+TOD_MAPSOLVE) / self%numsamp(band), 100*t(b+TOD_MAPSOLVE)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD load-balancing           = ', t(b+TOD_WAIT)  /  self%numsamp(band),   100*t(b+TOD_WAIT)/t(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD MPI operations           = ', t(b+TOD_MPI)  /  self%numsamp(band),   100*t(b+TOD_MPI)/t(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD scan (de)allocation      = ', t(b+TOD_ALLOC)  /  self%numsamp(band),   100*t(b+TOD_ALLOC)/t(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      TOD instrument correction    = ', t(b+TOD_INSTCORR)  /  self%numsamp(band),   100*t(b+TOD_INSTCORR)/t(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      Zodiacal Light model         = ', t(b+TOD_ZODI)     / self%numsamp(band), 100*t(b+TOD_ZODI)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      4D map output                = ', t(b+TOD_4D)     / self%numsamp(band), 100*t(b+TOD_4D)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      Miscellaneous file writing   = ', t(b+TOD_WRITE)     / self%numsamp(band), 100*t(b+TOD_WRITE)/T(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h",f10.2,"%")') '      Other                        = ', (t(b+TOD_TOT)-sum(t(b+3:b+NUM_TOD))) / self%numsamp(band), 100*(t(b+TOD_TOT)-sum(t(b+3:b+NUM_TOD)))/t(b+TOD_TOT)
          write(unit,fmt='(a,f12.3,"h")') '      Total TOD                    = ', t(b+TOD_TOT)     / self%numsamp(band)
       end do
       close(unit)
    end if

  end subroutine comm_timer_dumpASCII

end module comm_timing_mod
