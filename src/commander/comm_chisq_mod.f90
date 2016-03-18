module comm_chisq_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  implicit none


contains

  subroutine compute_chisq(chisq_map, chisq_fullsky)
    implicit none
    class(comm_map), intent(inout), optional :: chisq_map
    real(dp),        intent(out),   optional :: chisq_fullsky

    integer(i4b) :: i, j, k, p
    real(dp)     :: t1, t2
    class(comm_map), pointer :: res, chisq_sub

    if (present(chisq_fullsky) .or. present(chisq_map)) then
       if (present(chisq_fullsky)) chisq_fullsky = 0.d0
       if (present(chisq_map))     chisq_map%map = 0.d0
       do i = 1, numband
          res => compute_residual(i)
          call data(i)%N%sqrtInvN(res)
          res%map = res%map**2
          
          if (present(chisq_map)) then
             chisq_sub => comm_map(chisq_map%info)
             call res%udgrade(chisq_sub)
             chisq_map%map = chisq_map%map + chisq_sub%map * (res%info%npix/chisq_sub%info%npix)
          end if
          if (present(chisq_fullsky)) chisq_fullsky = chisq_fullsky + sum(res%map)
          call res%dealloc()
       end do
    end if

  end subroutine compute_chisq

  function compute_residual(band, exclude_comps) result (res)
    implicit none

    integer(i4b),                     intent(in)           :: band
    character(len=512), dimension(:), intent(in), optional :: exclude_comps
    class(comm_map),    pointer                            :: res

    real(dp) :: t1, t2, t3, t4
    class(comm_comp),    pointer :: c
    real(dp),     allocatable, dimension(:,:) :: map, alm
    integer(i4b), allocatable, dimension(:)   :: pix
    integer(i4b) :: ierr
    
    ! Initialize to full data set
    res => comm_map(data(band)%info)

    ! Compute predicted signal for this band
    c => compList
    !call wall_time(t1)
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          !allocate(map(0:res%info%np-1,res%info%nmaps))          
          !map     = c%getBand(band)
          !res%map = res%map - map
          !deallocate(map)
          allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))          
          !call wall_time(t3)
          alm     = c%getBand(band, alm_out=.true.)
          !call wall_time(t4)
          !if (data(band)%info%myid == 0) write(*,*) 'CPU getband = ', t4-t3
          !call wall_time(t3)
          call res%add_alm(alm, c%x%info)
          !call wall_time(t4)
          !if (data(band)%info%myid == 0) write(*,*) 'CPU add = ', t4-t3
          deallocate(alm)
       end select
       c => c%next()
    end do
    !call wall_time(t2)
    !if (data(band)%info%myid == 0) write(*,*) 'CPU build = ', t2-t1

    !call wall_time(t1)
    call res%Y()
    !call wall_time(t2)
    !if (data(band)%info%myid == 0) write(*,*) 'CPU Y = ', t2-t1
    
    ! Compute residual map
    !call wall_time(t1)
    res%map = data(band)%map%map - res%map
    !call wall_time(t2)
    !if (data(band)%info%myid == 0) write(*,*) 'CPU subtract = ', t2-t1

    ! Clean up
    nullify(c)

  end function compute_residual
    

end module comm_chisq_mod
