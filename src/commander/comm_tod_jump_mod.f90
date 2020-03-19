module comm_tod_jump_mod
  use healpix_types
  use comm_utils
  use rngmod
  implicit none

  interface tod2file
    module procedure tod2file_sp
    module procedure tod2file_dp
    module procedure tod2file_int
  end interface tod2file
  
contains

  subroutine tod2file_sp(filename,d)
    implicit none
    character(len=*),                 intent(in)            :: filename
    real(sp),           dimension(:), intent(in)            :: d

    integer(i4b)                                            :: unit, io_error, length, i

    write(*,*) "Writing TOD to file - ", trim(filename)
    unit = 22

    length = size(d)

    open(unit,file=trim(filename),status='replace',action='write',iostat=io_error)
    do i = 1, length
      write(unit,*) d(i)
    end do
    
    close(unit)

  end subroutine tod2file_sp


  subroutine tod2file_dp(filename,d)
    implicit none
    character(len=*),                 intent(in)            :: filename
    real(dp),           dimension(:), intent(in)            :: d

    integer(i4b)                                            :: unit, io_error, length, i

    write(*,*) "Writing TOD to file - ", trim(filename)
    unit = 22

    length = size(d)

    open(unit,file=trim(filename),status='replace',action='write',iostat=io_error)
    do i = 1, length
      write(unit,*) d(i)
    end do
    
    close(unit)
  end subroutine tod2file_dp


  subroutine tod2file_int(filename,d)
    implicit none
    character(len=*),                 intent(in)            :: filename
    integer(i4b),       dimension(:), intent(in)            :: d

    integer(i4b)                                            :: unit, io_error, length, i

    write(*,*) "Writing TOD to file - ", trim(filename)
    unit = 22

    length = size(d)

    open(unit,file=trim(filename),status='replace',action='write',iostat=io_error)
    do i = 1, length
      write(unit,*) d(i)
    end do
    
    close(unit)
  end subroutine tod2file_int


  real(dp) function welford_std(x_new, x_old, std_old, mean_new, mean_old, N)
    implicit none

    real(dp), intent(in)               :: x_new, x_old, std_old, mean_old
    integer(i4b), intent(in)           :: N
    real(dp), intent(inout)            :: mean_new


    mean_new = mean_old + (x_new-x_old)/N

    welford_std = sqrt(std_old**2 + (x_new**2-x_old**2-N*(mean_new**2-mean_old**2))/(N-1))

  end function welford_std


  real(dp) function std(x)
    implicit none

    real(sp), intent(in), dimension(:)   :: x
    real(dp)                             :: mean
    integer(i4b)                         :: n 

    n = size(x)
    mean = sum(x(1:n))/n
    std = sqrt(sum((x(1:n)-mean)**2)/(n-1))

  end function std


  real(dp) function std_nan(x)
    implicit none

    real(sp), intent(in), dimension(:)   :: x
    real(dp)                             :: mean
    integer(i4b)                         :: n 

    real(dp)                             :: st_dev
    integer(i4b)                         :: counter, i
    real(sp), allocatable, dimension(:)  :: x_no_nan

    n = size(x)
    counter = 0

    st_dev = std(x)
    if (isnan(st_dev)) then
      do i=1, n
         if (.not. isnan(x(i))) then
            counter = counter + 1
         end if
      end do
      allocate(x_no_nan(counter))
      counter = 1
      do i=1, n
         if (.not. isnan(x(i))) then
            x_no_nan(counter) = x(i)
            counter = counter + 1
         end if
      end do
      st_dev = std(x_no_nan)
    end if

    std_nan = st_dev

  end function std_nan


  real(dp) function mean_nan(x)
   implicit none

   real(sp), dimension(:), intent(in)    :: x

   integer(i4b)                          :: n, i, counter, counter2
   real(dp)                              :: avg
   real(sp), allocatable, dimension(:)   :: x_no_nan

   counter = 0
   n = size(x)
   avg = sum(x)/n
   
   if (isnan(avg)) then
      do i=1, n
         if (.not. isnan(x(i))) then
            counter = counter + 1
         end if
      end do
      allocate(x_no_nan(counter))
      counter2 = 1
      do i=1, n
         if (.not. isnan(x(i))) then
            x_no_nan(counter2) = x(i)
            counter2 = counter2 + 1
         end if
      end do
      avg = sum(x_no_nan)/counter
   end if

   mean_nan = avg

  end function mean_nan


  subroutine gap_fill_linear(tod,flag,tod_gapfill,handle,noise)
    implicit none
    real(sp),     dimension(:), intent(in)     :: tod
    integer(i4b), dimension(:), intent(in)     :: flag
    real(sp),     dimension(:), intent(inout)  :: tod_gapfill
    logical,                    intent(in)     :: noise
    type(planck_rng),            intent(inout) :: handle
 
    character(len=50)                          :: filename
    integer(i4b)                               :: tod_len, i, j, counter, marker, N
    logical                                    :: switch
    real(dp)                                   :: mean_low, mean_high
    real(sp)                                   :: std1, std2, std_mean
 
    ! noise: If "true", then white noise is added ontop of the interpolation
 
    N = 100
    write(*,*) "Routine: Gap fill"
    filename = 'gapfill/gapfill_avg.txt'
 
    tod_gapfill = tod
    where (flag==1) tod_gapfill = nan   
 
    tod_len = size(tod)
    marker = 0
 
    do i=1, tod_len
       if (i<=marker) cycle
       if (isnan(tod_gapfill(i))) then
          counter = 0
          switch = .true.
          do while (switch .and. ((i+counter)<=tod_len))
 
             if ((i+counter==tod_len) .and. isnan(tod_gapfill(i+counter))) then
                tod_gapfill(i:counter) = 0
             end if
              
             if (.not. isnan(tod_gapfill(i+counter))) then
                if (i==1) then
                   tod_gapfill(i:counter) = 0
                else
 
                   mean_low  = mean_nan(tod_gapfill(i-1-N:i-1))
                   mean_high = mean_nan(tod_gapfill(i+counter:i+counter+N))
                   
                   if (noise) then
                      std1 = std_nan(tod_gapfill(i-1-N:i-1))
                      std2 = std_nan(tod_gapfill(i+counter:i+counter+N))
                      std_mean = (std1+std2)/2.0
                      if (isnan(std_mean)) std_mean = 0
                   end if
                   
                   do j=1, counter
                      tod_gapfill(i+j-1) = mean_low + (j)*(mean_high-mean_low)/(counter+1)
                      if (noise) tod_gapfill(i+j-1) = tod_gapfill(i+j-1) + std_mean*rand_gauss(handle)
                   end do   
                end if
                switch = .false.
                marker = i+counter-1
             else
                counter = counter + 1
             end if
              
          end do
       end if
    end do
 
  end subroutine gap_fill_linear


  subroutine jump_scan(tod,flag,jumps,offset_range,offset_level,handle)
    implicit none
    real(sp),     dimension(:),   intent(in)    :: tod
    integer(i4b), dimension(:),   intent(inout) :: flag
    integer(i4b), dimension(:),   intent(inout) :: jumps
    integer(i4b), allocatable,    dimension(:,:), intent(inout) :: offset_range
    real(sp),    allocatable,     dimension(:),   intent(inout) :: offset_level
    type(planck_rng),             intent(inout) :: handle
 
    real(sp), allocatable, dimension(:)        :: tod_gapfill
    real(dp), allocatable, dimension(:)        :: rolling_std, rolling_std_flagged
    integer(i4b)                               :: tod_len, N, i, threshold, num_offsets, counter, low, high, N_delta
    real(dp)                                   :: std_old, mean_old, mean_new, x_new, x_old, st_dev, std_test, med, delta, delta_l, delta_r
    character(len=100)                         :: filename
    logical                                    :: switch, first_call   
 
    write(*,*) 'Routine: Jump scan'
 
    N = 100
    threshold = 2
    tod_len = size(tod)
    jumps(:) = 0
 
    ! Interpolate cosmic ray gaps
    allocate(tod_gapfill(tod_len))
    call gap_fill_linear(tod,flag,tod_gapfill,handle,.false.)
 
    ! Compute rolling standard deviation
    allocate(rolling_std(tod_len))
    allocate(rolling_std_flagged(tod_len))
    rolling_std = nan
    
    do i=N+1, tod_len-N
       if ((i==N+1) .or. (modulo(i,100)==0)) then
          std_old = std(tod_gapfill(i-N:i+N))
          mean_old = sum(tod_gapfill(i-N:i+N))/(2*N+1)
       else
          x_new = tod_gapfill(i+N)
          x_old = tod_gapfill(i-N-1)
 
          st_dev = welford_std(x_new, x_old, std_old, mean_new, mean_old, 2*N+1)
         
          mean_old = mean_new
          std_old = st_dev
       end if
       rolling_std(i) = std_old
    end do
 
    ! Compute median
    rolling_std_flagged = rolling_std
    where (flag==1) rolling_std_flagged = nan
    med = median(rolling_std_flagged)
    
    ! Do the flagging
    where (rolling_std > (threshold*med)) jumps = 1
 
    ! Find number of offsets so that offset list can be allocated
    num_offsets = 0
    switch = .false.
    do i=1, tod_len
       if (jumps(i)==1 .and. switch) then
          num_offsets = num_offsets + 1
          switch = .false.
       elseif (jumps(i)==0) then
          switch = .true.
       end if
    end do
    if (jumps(tod_len)==0) num_offsets = num_offsets + 1
    
 
 
    allocate(offset_range(num_offsets,2))
    allocate(offset_level(num_offsets))
 
    ! Define offset regions
    switch = .true.
    first_call = .true.
    counter = 1
    do i=1, tod_len
       if (first_call .and. jumps(i)==0) then
          offset_range(counter,1) = i
          first_call = .false.
       elseif (jumps(i)==1 .and. switch .and. first_call==.false.) then
          offset_range(counter,2) = i-1
          switch = .false.
          counter = counter + 1
       elseif (jumps(i)==0 .and. switch==.false.) then
          offset_range(counter,1) = i
          switch = .true.
       end if
    end do
    if (jumps(tod_len)==0) offset_range(counter,2) = tod_len
 
    
    ! Compute average and monopole template
    N_delta = 2000
    do i=1, num_offsets
       low  = offset_range(i,1)
       high = offset_range(i,2)
       offset_level(i) = sum(tod_gapfill(low:high))/(high-low+1)
 
       ! Compute correction for drift
       if (i==1) then
          delta = 0
       else
          delta_l = sum(tod_gapfill(low:low+N_delta))/(N_delta+1) - offset_level(i)
          delta = delta_r - delta_l
       end if
 
       offset_level(i) = offset_level(i) - delta
 
       delta_r = sum(tod_gapfill(high-N_delta:high))/(N_delta+1) - offset_level(i)
    end do
    
    ! Add jump regions to full flags
    where (jumps==1) flag=1
  
  end subroutine jump_scan


  subroutine expand_offset_list(offset_range,offset_level,s_jump,scan,det)
    implicit none
    integer(i4b), dimension(:,:), intent(in)     :: offset_range
    real(sp),     dimension(:),   intent(in)     :: offset_level
    real(sp),     dimension(:),   intent(inout)  :: s_jump
    integer(i4b),                 intent(in)     :: scan, det
 
    integer(i4b)                                 :: n_offsets, i, i_scan, i_det, i_low, i_high
 
    write(*,*) "Routine: Expand Offset List"
 
    n_offsets = size(offset_level)
 
    s_jump = 0
 
    do i=1, n_offsets
       i_scan = offset_range(i,1)
       i_det  = offset_range(i,2)
 
       if ((i_scan==scan) .and. (i_det==det)) then
          i_low  = offset_range(i,3)
          i_high = offset_range(i,4)
 
          s_jump(i_low:i_high) = offset_level(i)      
       end if
 
    end do
 
  end subroutine expand_offset_list





end module comm_tod_jump_mod
