module comm_chisq_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  implicit none


contains

  subroutine compute_chisq(chisq_fullsky)
    implicit none

    real(dp),                              intent(out), optional :: chisq_fullsky

    integer(i4b) :: i
    real(dp)     :: t1, t2
    class(comm_map), pointer :: res

    if (present(chisq_fullsky)) then
       chisq_fullsky = 0.d0
       do i = 1, numband
          res => compute_residual(i)
          call data(i)%N%sqrtInvN(res)
          chisq_fullsky = chisq_fullsky + sum(res%map**2)
       end do
    end if

  end subroutine compute_chisq

  function compute_residual(band, exclude_comps) result (res)
    implicit none

    integer(i4b),                     intent(in)           :: band
    character(len=512), dimension(:), intent(in), optional :: exclude_comps
    class(comm_map),    pointer                            :: res

    class(comm_comp),    pointer :: c
    real(dp),     allocatable, dimension(:,:) :: map
    integer(i4b), allocatable, dimension(:)   :: pix
    integer(i4b) :: ierr
    
    ! Initialize to full data set
    res => comm_map(data(band)%map)

    ! Subtract all components
    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          map     = c%getBand(band)
          res%map = res%map - map
       end select
       c => c%next()
    end do

    nullify(c)

  end function compute_residual
    

end module comm_chisq_mod
