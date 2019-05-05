module comm_chisq_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_ptsrc_comp_mod
  use comm_template_comp_mod
  implicit none


contains

  subroutine compute_chisq(comm, chisq_map, chisq_fullsky)
    implicit none
    integer(i4b),    intent(in)              :: comm
    class(comm_map), intent(inout), optional :: chisq_map
    real(dp),        intent(out),   optional :: chisq_fullsky

    integer(i4b) :: i, j, k, p, ierr
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

    if (present(chisq_fullsky)) then
       call mpi_allreduce(MPI_IN_PLACE, chisq_fullsky, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    end if

  end subroutine compute_chisq

  function compute_residual(band, exclude_comps, cg_samp_group) result (res)
    implicit none

    integer(i4b),                     intent(in)           :: band
    character(len=512), dimension(:), intent(in), optional :: exclude_comps
    integer(i4b),                     intent(in), optional :: cg_samp_group
    class(comm_map),    pointer                            :: res

    integer(i4b) :: i
    logical(lgt) :: skip
    real(dp)     :: t1, t2, t3, t4
    class(comm_comp),    pointer :: c
    real(dp),     allocatable, dimension(:,:) :: map, alm
    integer(i4b), allocatable, dimension(:)   :: pix
    integer(i4b) :: ierr
    logical(lgt) :: nonzero
    class(comm_map), pointer :: ptsrc
    
    ! Initialize to full data set
    res   => comm_map(data(band)%info)  ! Diffuse
    ptsrc => comm_map(data(band)%info)  ! Compact

    ! Compute predicted signal for this band
    c => compList
    nonzero = .false.
    do while (associated(c))
       skip = .false.
       if (present(exclude_comps)) then
          ! Skip if the component is requested to be excluded
          do i = 1, size(exclude_comps)
             if (trim(c%label) == trim(exclude_comps(i))) skip = .true.
          end do
       end if
       if (present(cg_samp_group)) then
          if (c%cg_samp_group == cg_samp_group) skip = .true.
       end if
       if (skip) then
          c => c%next()
          cycle
       end if

       select type (c)
       class is (comm_diffuse_comp)
          allocate(alm(0:data(band)%info%nalm-1,data(band)%info%nmaps))          
          alm     = c%getBand(band, alm_out=.true.)
          res%alm = res%alm + alm
          !call res%add_alm(alm, c%x%info)
          deallocate(alm)
          nonzero = .true.
       class is (comm_ptsrc_comp)
          allocate(map(0:data(band)%info%np-1,data(band)%info%nmaps))
          map       = c%getBand(band)
          ptsrc%map = ptsrc%map + map
          deallocate(map)
       class is (comm_template_comp)
          allocate(map(0:data(band)%info%np-1,data(band)%info%nmaps))
          map       = c%getBand(band)
          ptsrc%map = ptsrc%map + map
          deallocate(map)
       end select
       c => c%next()
    end do
    if (nonzero) call res%Y()

    ! Compute residual map
    res%map = data(band)%map%map - res%map - ptsrc%map

    ! Clean up
    nullify(c)
    call ptsrc%dealloc()

  end function compute_residual

  subroutine output_signals_per_band(outdir, postfix)
    implicit none
    character(len=*), intent(in) :: outdir, postfix
    
    integer(i4b) :: i
    logical(lgt) :: skip
    character(len=1024) :: filename
    class(comm_comp), pointer :: c
    class(comm_map),  pointer :: out
    real(dp),     allocatable, dimension(:,:) :: map, alm
    integer(i4b), allocatable, dimension(:)   :: pix
    
    do i = 1, numband
       out => comm_map(data(i)%info)  

       ! Compute predicted signal for this band
       c => compList
       do while (associated(c))
          if (trim(c%type) == 'md') then
             c => c%next()
             cycle
          end if

          skip    = .false.
          out%alm = 0.d0
          out%map = 0.d0
          select type (c)
          class is (comm_diffuse_comp)
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))          
             alm     = c%getBand(i, alm_out=.true.)
             call out%add_alm(alm, c%x%info)
             call out%Y()
             deallocate(alm)
          class is (comm_ptsrc_comp)
             allocate(map(0:data(i)%info%np-1,data(i)%info%nmaps))
             map     = c%getBand(i)
             out%map = out%map + map
             deallocate(map)
          class is (comm_template_comp)
             if (c%band /= i) skip = .true.
             if (.not. skip) then
                allocate(map(0:data(i)%info%np-1,data(i)%info%nmaps))
                map     = c%getBand(i)
                out%map = out%map + map
                deallocate(map)
             end if
          end select
          filename = trim(outdir)//'/'//trim(c%label)//'_'//trim(data(i)%label)//'_'//trim(postfix)//'.fits'
          !call data(i)%apply_proc_mask(out)
          if (.not. skip) call out%writeFITS(filename)
          c => c%next()
       end do
    end do

    ! Clean up
    nullify(c)
    call out%dealloc

  end subroutine output_signals_per_band

  subroutine get_sky_signal(band, det, map_out)
    implicit none
    integer(i4b),    intent(in)     :: band, det
    class(comm_map), pointer        :: map_out

    logical(lgt) :: skip
    class(comm_comp), pointer :: c
    real(dp),     allocatable, dimension(:,:) :: map, alm

    ! Allocate map
    map_out => comm_map(data(band)%info)  

    ! Compute predicted signal for this band
    c => compList
    do while (associated(c))
       if (trim(c%type) == 'md') then
          c => c%next()
          cycle
       end if

       skip    = .false.
       map_out%alm = 0.d0
       map_out%map = 0.d0
       select type (c)
       class is (comm_diffuse_comp)
          allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))          
          alm     = c%getBand(band, alm_out=.true.)
          call map_out%add_alm(alm, c%x%info)
          call map_out%Y()
          deallocate(alm)
       class is (comm_ptsrc_comp)
          allocate(map(0:data(band)%info%np-1,data(band)%info%nmaps))
          map         = c%getBand(band)
          map_out%map = map_out%map + map
          deallocate(map)
       class is (comm_template_comp)
          if (c%band /= band) skip = .true.
          if (.not. skip) then
             allocate(map(0:data(band)%info%np-1,data(band)%info%nmaps))
             map         = c%getBand(band)
             map_out%map = map_out%map + map
             deallocate(map)
          end if
       end select
       c => c%next()
    end do

  end subroutine get_sky_signal


end module comm_chisq_mod
