module comm_chisq_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  implicit none


contains

  subroutine compute_chisq(chisq_fullsky)
    implicit none

    real(dp), intent(out) :: chisq_fullsky

    


  end subroutine compute_chisq

  function compute_residual(band) result (res)
    implicit none

    integer(i4b),             intent(in) :: band
    class(comm_map), pointer             :: res

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

    call data(band)%map%writeFITS('map.fits')
    call res%writeFITS('res.fits')
    call mpi_finalize(ierr)
    stop

  end function compute_residual
    

end module comm_chisq_mod
