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
    if (n==1) then
      std = 0
      return
    end if
    mean = sum(x(1:n))/n
    std = sqrt(sum((x(1:n)-mean)**2)/(n-1))

  end function std


  real(dp) function std_nan(x)
    implicit none
    real(sp), intent(in), dimension(:)   :: x
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


  real(dp) function std_flagged(x,flag)
    implicit none
    real(sp),     dimension(:), intent(in)   :: x
    integer(i4b), dimension(:), intent(in)   :: flag

    integer(i4b)                         :: n 
    integer(i4b)                         :: counter, i
    real(sp), allocatable, dimension(:)  :: x_not_flagged

    n = size(x)
    if (sum(flag)==0) then
      std_flagged = std(x)
    else
      counter = 0
      do i=1, n
         if (flag(i)==0) then
            counter = counter + 1
         end if
      end do
      allocate(x_not_flagged(counter))
      counter = 1
      do i=1, n
         if (flag(i)==0) then
            x_not_flagged(counter) = x(i)
            counter = counter + 1
         end if
      end do
      std_flagged = std(x_not_flagged)
    end if

  end function std_flagged


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


  real(dp) function mean_flagged(x,flag)
   implicit none
   real(sp),     dimension(:), intent(in)    :: x
   integer(i4b), dimension(:), intent(in)    :: flag

   integer(i4b)                          :: n, i, counter, counter2
   real(sp), allocatable, dimension(:)   :: x_not_flagged

   n = size(x)
   if (sum(flag)==0) then
    mean_flagged = sum(x)/n
   else
    counter = 0
    do i=1, n
      if (flag(i)==0) then
          counter = counter + 1
      end if
    end do
    allocate(x_not_flagged(counter))
    counter2 = 0
    do i=1, n
      if (flag(i)==0) then
          counter2 = counter2 + 1
          x_not_flagged(counter2) = x(i)
      end if
    end do
    mean_flagged = sum(x_not_flagged)/counter
   end if

  end function mean_flagged


  real(dp) function median_flagged(x,flag)
   implicit none
   real(dp),     dimension(:), intent(in)    :: x
   integer(i4b), dimension(:), intent(in)    :: flag

   real(dp), allocatable, dimension(:)   :: tmp
   integer(i4b)                          :: n, i, counter, counter2
   real(dp), allocatable, dimension(:)   :: x_not_flagged

   n = size(x)
   if (sum(flag)==0) then
    allocate(tmp(n))
    tmp = x
    call QuickSort_real(tmp)
    median_flagged = tmp(n/2+1)
    deallocate(tmp)
   else
    counter = 0
    do i=1, n
      if (flag(i)==0) then
          counter = counter + 1
      end if
    end do
    allocate(x_not_flagged(counter))
    counter2 = 0
    do i=1, n
      if (flag(i)==0) then
          counter2 = counter2 + 1
          x_not_flagged(counter2) = x(i)
      end if
    end do
    call QuickSort_real(x_not_flagged)
    median_flagged = x_not_flagged(counter2/2+1)
   end if

  end function median_flagged


  subroutine gap_fill_linear_deprecated(tod,flag,tod_gapfill,handle,noise)
    implicit none
    real(sp),     dimension(:), intent(in)     :: tod
    integer(i4b), dimension(:), intent(in)     :: flag
    real(sp),     dimension(:), intent(inout)  :: tod_gapfill
    type(planck_rng),           intent(inout)  :: handle
    logical,                    intent(in)     :: noise

    character(len=50)                          :: filename
    integer(i4b)                               :: tod_len, i, j, counter, marker, N
    logical                                    :: switch
    real(dp)                                   :: mean_low, mean_high
    real(sp)                                   :: std1, std2, std_mean
 
    ! noise: If "true", then white noise is added ontop of the interpolation
 
    N = 100
    write(*,*) "Routine: Gap fill deprecated"
 
    tod_gapfill = tod
    tod_len = size(tod)
    marker = 0
 
    do i=1, tod_len
       if (i<=marker) cycle
       if (flag(i)==1) then
          counter = 0
          switch = .true.
          do while (switch .and. ((i+counter)<=tod_len))
 
             if ((i+counter==tod_len) .and. (flag(i+counter)==1)) then
                tod_gapfill(i:counter) = 0 
             end if
              
             if (flag(i+counter)==0) then
                if (i==1) then
                   tod_gapfill(i:counter) = 0 
                else
 
                   mean_low  = mean_flagged(tod_gapfill(i-1-N:i-1),flag(i-1-N:i-1))
                   mean_high = mean_flagged(tod_gapfill(i+counter:i+counter+N),flag(i+counter:i+counter+N))
                   
                   if (noise) then
                      std1 = std_flagged(tod_gapfill(i-1-N:i-1),flag(i-1-N:i-1))
                      std2 = std_flagged(tod_gapfill(i+counter:i+counter+N),flag(i+counter:i+counter+N))
                      std_mean = (std1+std2)/2.0
                      if (isnan(std_mean)) std_mean = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
 
  end subroutine gap_fill_linear_deprecated


  subroutine gap_fill_linear(tod,flag,tod_gapfill,handle,noise)
   implicit none
   real(sp),     dimension(:), intent(in)    :: tod
   integer(i4b), dimension(:), intent(in)    :: flag
   real(sp),     dimension(:), intent(inout) :: tod_gapfill
   type(planck_rng),           intent(inout) :: handle
   logical,                    intent(in)    :: noise

   character(len=50)                      :: filename
   integer(i4b)                           :: N, tod_len, i, counter, marker, j, low, high
   logical                                :: counting
   real(dp)                               :: mean_low, mean_high
   real(sp)                               :: std_low, std_high, std_combined

   ! noise: If "true", then white noise is added ontop of the interpolation
   !write(*,*) "Routine: Gap fill"
   N = 100

   tod_gapfill = tod
   tod_len = size(tod)
   counting = .false.
   do i=1, tod_len

      ! If you reach a flagged stretch, start counting and mark position
      if ((flag(i)==1) .and. (.not. counting)) then
         marker = i
         counting = .true.
         counter = 1 
         
         
      ! When leaving flagged stretch or flagged stretch reaches the end of tod -> do the interpolation
      else if (((flag(i)==0) .or. ((flag(i)==1) .and. (i==tod_len))) .and. counting) then
         counting = .false.
         
         ! Compute mean before and after the gap
         low  = max(1,marker-N)
         high = min(i+N-1,tod_len)
         if (marker==1) then
            mean_high = mean_flagged(tod_gapfill(i:high),flag(i:high))
            mean_low = mean_high
         else if (i==tod_len) then
            mean_low  = mean_flagged(tod_gapfill(low:marker-1),flag(low:marker-1))
            mean_high = mean_low
         else
            mean_low  = mean_flagged(tod_gapfill(low:marker-1),flag(low:marker-1))
            mean_high = mean_flagged(tod_gapfill(i:high),flag(i:high))
         end if
         
         ! Compute std before and after the gap
         if (noise) then
            if (marker==1) then
               std_high = std_flagged(tod_gapfill(i:high),flag(i:high))
               std_low = std_high
            else if (i==tod_len) then
               std_low = std_flagged(tod_gapfill(low:marker-1),flag(low:marker-1))
               std_high = std_low
            else            
               std_low  = std_flagged(tod_gapfill(low:marker-1),flag(low:marker-1))
               std_high = std_flagged(tod_gapfill(i:high),flag(i:high))
            end if
            std_combined = (std_low + std_high)/2
            if ((std_low==0) .or. (std_high==0)) std_combined = max(std_low,std_high)
         end if
         
         ! Do the interpolation
         do j=marker, i-1
            tod_gapfill(j) = mean_low + (j-marker+1)*(mean_high-mean_low)/(counter+1)
            if (noise) tod_gapfill(j) = tod_gapfill(j) + std_combined*rand_gauss(handle)
         end do
         
      ! When counting, keep counting
      else if ((flag(i)==1) .and. counting) then
         counter = counter + 1
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
    real(dp), allocatable, dimension(:)        :: rolling_std
    integer(i4b)                               :: tod_len, N, i, threshold, num_offsets, counter, low, high, N_delta, marker, len_min
    real(dp)                                   :: std_old, mean_old, mean_new, x_new, x_old, st_dev, std_test, med, delta, delta_l, delta_r
    character(len=100)                         :: filename
    logical                                    :: switch, first_call, counting   
 
    !write(*,*) 'Routine: Jump scan'
 
    N = 100
    threshold = 2
    tod_len = size(tod)
    jumps(:) = 0
    len_min = 2000
 
    ! Interpolate cosmic ray gaps
    allocate(tod_gapfill(tod_len))
    call gap_fill_linear(tod,flag,tod_gapfill,handle,.false.)
 
    ! Compute rolling standard deviation
    allocate(rolling_std(tod_len))
    rolling_std = 0
   
    do i=N+1, tod_len-N
       if ((i==N+1) .or. (modulo(i,1)==0)) then
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
    med = median_flagged(rolling_std,flag)
    
    ! Do the flagging
    where (rolling_std > (threshold*med)) jumps = 1

    ! Remove regions that are shorter than the minimum allowed tod length
    counter = 0
    counting = .false.
    do i=1, tod_len
      if ((jumps(i)==0) .and. (.not. counting)) then
         counter = counter + 1
         counting = .true.
         marker = i
      elseif ((jumps(i)==0) .and. counting) then
         counter = counter + 1
         if ((i==tod_len) .and. (counter<len_min)) jumps(marker:i) = 1
      elseif ((jumps(i)==1) .and. counting) then
         counting = .false.
         if (counter<len_min) jumps(marker:i-1) = 1
         counter = 0
      end if
    end do
         

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
    do i=1, num_offsets
       low  = offset_range(i,1)
       high = offset_range(i,2)
       offset_level(i) = sum(tod_gapfill(low:high))/(high-low+1)
 
       ! Compute correction for drift
       if (i==1) then
          delta = 0
       else
          delta_l = sum(tod_gapfill(low:low+len_min))/(len_min+1) - offset_level(i)
          delta = delta_r - delta_l
       end if
 
       offset_level(i) = offset_level(i) - delta
 
       if (i==num_offsets) exit
       delta_r = sum(tod_gapfill(high-len_min:high))/(len_min+1) - offset_level(i)
    end do
    
    ! Add jump regions to full flags
    where (jumps==1) flag=1
  
  end subroutine jump_scan


  subroutine expand_offset_list(offset_range,offset_level,s_jump)
    implicit none
    integer(i4b), dimension(:,:), intent(in)     :: offset_range
    real(sp),     dimension(:),   intent(in)     :: offset_level
    real(sp),     dimension(:),   intent(inout)  :: s_jump
 
    integer(i4b)                                 :: n_offsets, i, i_low, i_high
 
    !write(*,*) "Routine: Expand Offset List"
 
    n_offsets = size(offset_level)
 
    s_jump = 0
 
    do i=1, n_offsets
       i_low  = offset_range(i,1)
       i_high = offset_range(i,2)

       s_jump(i_low:i_high) = s_jump(i_low:i_high) + offset_level(i)      
    end do
 
  end subroutine expand_offset_list


  subroutine update_offset_list(offset_range_new,offset_level_new,offset_range_global,offset_level_global)
   implicit none
   integer(i4b),              dimension(:,:), intent(in)    :: offset_range_new
   real(sp),                  dimension(:),   intent(in)    :: offset_level_new
   integer(i4b), allocatable, dimension(:,:), intent(inout) :: offset_range_global
   real(sp),     allocatable, dimension(:),   intent(inout) :: offset_level_global

   integer(i4b)                                             :: n_row_new, n_row_global
   integer(i4b), allocatable, dimension(:,:)                :: offset_range_temp
   real(sp),     allocatable, dimension(:)                  :: offset_level_temp

   n_row_new    = size(offset_level_new)
   n_row_global = size(offset_level_global)

   if (n_row_new==1) return

   ! Update offset range
   allocate(offset_range_temp(n_row_global,2))
   offset_range_temp = offset_range_global
   deallocate(offset_range_global)
   allocate(offset_range_global(n_row_global+n_row_new,2))
   offset_range_global(1:n_row_global,:)                        = offset_range_temp
   offset_range_global(n_row_global+1:n_row_global+n_row_new,:) = offset_range_new
   deallocate(offset_range_temp)

   ! Update offset level
   allocate(offset_level_temp(n_row_global))
   offset_level_temp = offset_level_global
   deallocate(offset_level_global)
   allocate(offset_level_global(n_row_global+n_row_new))
   offset_level_global(1:n_row_global)                        = offset_level_temp
   offset_level_global(n_row_global+1:n_row_global+n_row_new) = offset_level_new
   deallocate(offset_level_temp)

  end subroutine update_offset_list


end module comm_tod_jump_mod
