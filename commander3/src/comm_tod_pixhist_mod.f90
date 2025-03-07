module comm_tod_pixhist_mod
  use comm_tod_mod
  use comm_tod_driver_mod
  implicit none



contains

  subroutine compute_min_max_nhit_per_pix(tod, map_sky, map_gain, procmask, procmask2, sub_zodi, &
      & map_min, map_max, map_nhit, polang)
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_gain
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2
    logical(lgt),                                 intent(in)    :: sub_zodi
    real(sp),            dimension(0:),           intent(out)   :: map_min, map_max, map_nhit
    real(dp),                                     intent(in),   optional :: polang

    integer(i4b) :: i, j, k, ierr, pix
    real(sp)     :: val
    type(comm_scandata) :: sd

    if (tod%nhorn /= 1) then
       write(*,*) 'compute_min_max_nhit_per_pix does not yet support multi-horn data'
       stop
    end if
    
    map_min  =  1e30
    map_max  = -1e30
    map_nhit =  0.
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle

       ! Prepare data
       call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2)

       ! Find min, max and nhits
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          do k = 1, sd%ntod
             if (iand(sd%flag(i,j), tod%flag0) .ne. 0) cycle
             pix = sd%pix(k,j,1)
             val = sd%tod(k,j)
             if (sub_zodi) val = val - tod%scans(i)%d(j)%gain * sd%s_zodi(k,j)
             map_min(pix)  = min(map_min(pix), val)
             map_max(pix)  = max(map_min(pix), val)
             map_nhit(pix) = map_nhit(pix) + 1.0
          end do
       end do

       ! Clean up
       call dealloc_scan_data(sd)
    end do
    
  end subroutine compute_min_max_nhit_per_pix


end module comm_tod_pixhist_mod
