module comm_timing_mod
  use comm_utils
  implicit none

  ! Global parameters
  integer(i4b), parameter, public :: NUM_GLOBAL    =  10
  integer(i4b), parameter, public :: TOT_RUNTIME   =  1
  integer(i4b), parameter, public :: TOT_AMPSAMP   =  2
  integer(i4b), parameter, public :: TOT_TODPROC   =  3
  integer(i4b), parameter, public :: TOT_SPECIND   =  4
  integer(i4b), parameter, public :: TOT_INIT      =  5
  integer(i4b), parameter, public :: TOT_FFT       =  6
  integer(i4b), parameter, public :: TOT_SHT       =  7
  integer(i4b), parameter, public :: TOT_OUTPUT    =  8
  integer(i4b), parameter, public :: TOT_GIBBSSAMP =  9
  integer(i4b), parameter, public :: TOT_CLS       =  10

  ! Channel specific parameters
  integer(i4b), parameter, public :: NUM_TOD       = 16
  integer(i4b), parameter, public :: TOD_INIT      =  1
  integer(i4b), parameter, public :: TOD_SL_PRE    =  2
  integer(i4b), parameter, public :: TOD_SL_INT    =  3
  integer(i4b), parameter, public :: TOD_PROJECT   =  4
  integer(i4b), parameter, public :: TOD_ORBITAL   =  5
  integer(i4b), parameter, public :: TOD_DECOMP    =  6
  integer(i4b), parameter, public :: TOD_ABSCAL    =  7
  integer(i4b), parameter, public :: TOD_RELCAL    =  8
  integer(i4b), parameter, public :: TOD_DELTAG    =  9
  integer(i4b), parameter, public :: TOD_NCORR     = 10
  integer(i4b), parameter, public :: TOD_XI_N      = 11
  integer(i4b), parameter, public :: TOD_MAPBIN    = 12
  integer(i4b), parameter, public :: TOD_MAPSOLVE  = 13
  integer(i4b), parameter, public :: TOD_ZODI      = 14
  integer(i4b), parameter, public :: TOD_IMBAL     = 15
  integer(i4b), parameter, public :: TOD_TOT       = 16

  private
  public comm_timing

  type comm_timing
     integer(i4b) :: numband, numsamp, comm, myid, n_tot
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
       self%t(timer)  = t2 - self%t1(timer)
       self%t1(timer) = 0.d0
    end if

  end subroutine comm_timer_stop

  subroutine comm_timer_incr_numsamp(self)
    ! 
    ! Routine for incrementing sample counter; used to output time per sample
    ! 
    ! Input variables:
    !    self     = comm_timing object
    ! 
    implicit none
    class(comm_timing),               intent(inout) :: self

    self%numsamp = self%numsamp + 1

  end subroutine comm_timer_incr_numsamp


  subroutine comm_timer_dumpASCII(self, filename)
    ! 
    ! Routine for outputting timing information
    ! 
    ! Input variables:
    !    self     = comm_timing object
    !    filename = output filename
    ! 
    implicit none
    class(comm_timing),              intent(inout) :: self
    character(len=*),                intent(in)    :: filename

    integer(i4b) :: unit, ierr, band, b
    real(dp), dimension(self%n_tot) :: t

    if (self%numsamp == 0) return

    call self%stop(TOT_RUNTIME)
    call self%start(TOT_RUNTIME)

    call mpi_reduce(self%t, t, self%n_tot, MPI_DOUBLE_PRECISION, &
         & MPI_SUM, 0, self%comm, ierr)
     
    if (self%myid == 0) then
       unit = getlun()
       open(unit,file=trim(filename), recl=1024)
       write(unit,*) 'Timing summary'
       write(unit,*) ''
       write(unit,*) '   Global total timers:'
       write(unit,*) '       Number of samples             = ', self%numsamp
       write(unit,*) '       Total runtime                 = ', t(TOT_RUNTIME)
       write(unit,*) '       Initialization                = ', t(TOT_INIT)
       write(unit,*) ''
       write(unit,*) '   Global per-sample timers:'
       write(unit,*) '       Chain output                  = ', t(TOT_OUTPUT)  / self%numsamp
       write(unit,*) '       Amplitude sampling            = ', t(TOT_AMPSAMP) / self%numsamp
       write(unit,*) '       Spectral index sampling       = ', t(TOT_SPECIND) / self%numsamp
       write(unit,*) '       Cls sampling                  = ', t(TOT_CLS)     / self%numsamp
       write(unit,*) '       TOD processing                = ', t(TOT_TODPROC) / self%numsamp
       write(unit,*) '       Total FFT                     = ', t(TOT_FFT)     / self%numsamp
       write(unit,*) '       Total SHT                     = ', t(TOT_SHT)     / self%numsamp
       write(unit,*) ''
       write(unit,*) '   Channel-specific global timers:'

       do band = 1, self%numband
          b = NUM_GLOBAL + (band-1)*NUM_TOD
          if (all(t(b+1:b+NUM_TOD) == 0.d0)) cycle
          write(unit,*) 
          write(unit,*) '     Channel ID                   = ', band
          write(unit,*) '     TOD initialization           = ', t(b+TOD_INIT)
          write(unit,*) '     TOD sidelobe precomputation  = ', t(b+TOD_SL_PRE)   / self%numsamp
          write(unit,*) '     TOD sidelobe interpolation   = ', t(b+TOD_SL_INT)   / self%numsamp
          write(unit,*) '     TOD sky-to-tod projection    = ', t(b+TOD_PROJECT)  / self%numsamp
          write(unit,*) '     TOD orbital dipole           = ', t(b+TOD_ORBITAL)  / self%numsamp
          write(unit,*) '     TOD decompression            = ', t(b+TOD_DECOMP)   / self%numsamp
          write(unit,*) '     TOD absolute calibration     = ', t(b+TOD_ABSCAL)   / self%numsamp
          write(unit,*) '     TOD relative calibration     = ', t(b+TOD_RELCAL)   / self%numsamp
          write(unit,*) '     TOD delta G calibration      = ', t(b+TOD_DELTAG)   / self%numsamp
          write(unit,*) '     TOD transmission imbalance   = ', t(b+TOD_IMBAL)    / self%numsamp
          write(unit,*) '     TOD correlated noise         = ', t(b+TOD_NCORR)    / self%numsamp
          write(unit,*) '     TOD corr noise PSD           = ', t(b+TOD_XI_N)     / self%numsamp
          write(unit,*) '     TOD binning                  = ', t(b+TOD_MAPBIN)   / self%numsamp
          write(unit,*) '     TOD map solution             = ', t(b+TOD_MAPSOLVE) / self%numsamp
          write(unit,*) '     Zodiacal Light model         = ', t(b+TOD_ZODI)     / self%numsamp
          write(unit,*) '     Total TOD                    = ', t(b+TOD_TOT)     / self%numsamp
       end do
       close(unit)
    end if

  end subroutine comm_timer_dumpASCII

end module comm_timing_mod
