module comm_nonlin_mod
  use comm_param_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_chisq_mod
  implicit none

contains

  subroutine sample_nonlin_params(cpar, handle)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    type(planck_rng),   intent(inout) :: handle    

    integer(i4b) :: i
    real(dp)     :: t1, t2
    class(comm_map), pointer :: res
    class(comm_comp),   pointer                    :: c
    real(dp),          allocatable, dimension(:,:) :: m

    call wall_time(t1)
    
    ! Initialize residual maps
    do i = 1, numband
       res             => compute_residual(i)
       data(i)%res%map =  res%map
       call res%dealloc()
    end do
    
    ! Sample spectral parameters for each signal component
    c => compList
    do while (associated(c))
       if (c%npar == 0) then
          c => c%next()
          cycle
       end if
       if (all(c%p_gauss(2,:) == 0.d0)) then
          c => c%next()
          cycle
       end if
       
       ! Add current component back into residual
       if (trim(c%class) /= 'ptsrc') then
          do i = 1, numband
             allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
             m               = c%getBand(i)
             data(i)%res%map = data(i)%res%map + m
             deallocate(m)
          end do
       end if

       ! Sample spectral parameters
       call c%sampleSpecInd(handle)

       ! Subtract updated component from residual
       if (trim(c%class) /= 'ptsrc') then
          do i = 1, numband
             allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
             m               = c%getBand(i)
             data(i)%res%map = data(i)%res%map - m
             deallocate(m)
          end do
       end if

       c => c%next()
    end do

    call wall_time(t2)
    if (cpar%myid == 0) write(*,*) 'CPU time specind = ', real(t2-t1,sp)
    
  end subroutine sample_nonlin_params

end module comm_nonlin_mod
