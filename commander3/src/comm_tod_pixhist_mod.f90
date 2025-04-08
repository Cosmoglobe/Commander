module comm_tod_pixhist_mod
  use comm_tod_mod
  use comm_tod_driver_mod
  implicit none

  private
  public compute_tod_pixhist

  integer(i4b), parameter :: NBIN_HIST      = 100
  integer(i4b), parameter :: MAX_ITER_HIST  = 20
  integer(i4b), parameter :: CUTOFF_N_EMPTY = 20

contains

  subroutine compute_tod_pixhist(tod, map_sky, map_gain, procmask, procmask2)
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_gain
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j, k, ierr, pix, npix_hist, q, bin(1), iter, nhit, hit, n_empty, j_cut
    real(sp)     :: val, center, mu, sigma
    logical(lgt) :: refine
    type(comm_scandata) :: sd
    character(len=6) :: pix_text
    integer(i4b), allocatable, dimension(:,:) :: hist
    real(sp), allocatable, dimension(:) :: delta
    real(sp),              dimension(NBIN_HIST) :: x, P
    
    if (tod%nhorn /= 1) then
       write(*,*) 'compute_tod_pixhist does not support multi-horn data'
       stop
    end if

    npix_hist = 12*tod%nside_pixhist**2
    q         = (tod%nside / tod%nside_pixhist)**2

    ! Find absolute min and max per low-res pixel
    allocate(tod%pixhist(5,0:npix_hist-1))  ! (mu, rms, nhit, min, max)
    allocate(delta(0:npix_hist-1))  ! (mu, rms, nhit, min, max)
    allocate(hist(0:NBIN_HIST,0:npix_hist-1))
    tod%pixhist       = 0.
    tod%pixhist(4,:)  =  1e30
    tod%pixhist(5,:)  = -1e30
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2, skip_zodi=.true.)
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          do k = 1, sd%ntod
             if (iand(sd%flag(k,j), tod%flag0) .ne. 0) cycle
             call ring2nest(tod%nside, sd%pix(k,j,1), pix)
             pix = pix / q
             tod%pixhist(4,pix)  = min(tod%pixhist(4,pix), sd%tod(k,j))
             tod%pixhist(5,pix)  = max(tod%pixhist(5,pix), sd%tod(k,j))
          end do
       end do

       ! Clean up
       call dealloc_scan_data(sd)
    end do
    call mpi_allreduce(MPI_IN_PLACE, tod%pixhist(4,:),  size(tod%pixhist(4,:)),  MPI_REAL, MPI_MIN, tod%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, tod%pixhist(5,:),  size(tod%pixhist(5,:)),  MPI_REAL, MPI_MAX, tod%comm, ierr)
    delta = (tod%pixhist(5,:)-tod%pixhist(4,:))/NBIN_HIST

    refine = .true.
    iter   = 0
    do while (refine)
       iter = iter + 1
       if (tod%myid == 0) write(*,*) 'Pixhist, band = ', trim(tod%freq), ', iter = ', iter
    
       ! Build histogram
       hist = 0.
       do i = 1, tod%nscan
          if (.not. any(tod%scans(i)%d%accept)) cycle
          call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2, skip_zodi=.true.)
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             do k = 1, sd%ntod
                if (iand(sd%flag(k,j), tod%flag0) .ne. 0) cycle
                call ring2nest(tod%nside, sd%pix(k,j,1), pix)
                pix = pix / q
                if (delta(pix) == 0.) cycle ! 0 or 1 samples in current pixel
                bin = int((sd%tod(k,j)-tod%pixhist(4,pix))/delta(pix))+1
                if (bin(1) < 1 .or. bin(1) > NBIN_HIST) bin(1) = 0 ! Sample out of range; discarded
                hist(bin(1),pix) = hist(bin(1),pix) + 1
             end do
          end do
          
          ! Clean up
          call dealloc_scan_data(sd)
       end do
       call mpi_allreduce(MPI_IN_PLACE, hist,  size(hist),  MPI_INTEGER, MPI_SUM, tod%comm, ierr)

       if (iter >= MAX_ITER_HIST) exit
       
       ! Check histogram quality, and refine limits
       refine = .false.
       do i = 0, npix_hist-1
          if (delta(i) == 0.) cycle ! 0 or 1 samples in current pixel
          
          nhit = sum(hist(1:NBIN_HIST,i))
          if (nhit == 0) cycle
          ! Cut bins from low side
          j       = 1
          j_cut   = -1
          hit     = 0
          n_empty = 0
          do while (hit < nhit/4)
             if (hist(j,i) == 0) then
                n_empty = n_empty + 1
             else
                if (n_empty >= CUTOFF_N_EMPTY) then
                   j_cut   = j
                   refine  = .true.
                end if
                n_empty = 0
             end if
             hit = hit + hist(j,i)
             j   = j+1
          end do
          if (j_cut /= -1) tod%pixhist(4,i) = tod%pixhist(4,i) + (j_cut-2)*delta(i)

          ! Cut bins from high side
          j       = NBIN_HIST
          j_cut   = -1
          hit     = 0
          n_empty = 0
          do while (hit < nhit/4)
             if (hist(j,i) == 0) then
                n_empty = n_empty + 1
             else
                if (n_empty >= CUTOFF_N_EMPTY) then
                   j_cut   = j
                   refine  = .true.
                end if
                n_empty = 0
             end if
             hit = hit + hist(j,i)
             j   = j-1
          end do
          if (j_cut /= -1) tod%pixhist(5,i) = tod%pixhist(5,i) - (NBIN_HIST-j_cut-1)*delta(i)
          
          ! Update bin width
          delta(i) = (tod%pixhist(5,i)-tod%pixhist(4,i))/NBIN_HIST
          !if (delta(i) == 0.) write(*,*) "zero", i, tod%pixhist(4,i), tod%pixhist(5,i), hist(:,i)
       end do
    end do

    ! Define cuts and output some histograms to disk
    do i = 0, npix_hist-1
       if (delta(i) == 0.) then
          ! 0 or 1 samples in current pixel; accept anything
          tod%pixhist(1,i) = 0.
          tod%pixhist(2,i) = 0.
          tod%pixhist(3,i) = 0.
          tod%pixhist(4,i) = -1e30
          tod%pixhist(5,i) =  1e30
          cycle
       end if
       
       nhit = sum(hist(1:NBIN_HIST,i))
       if (nhit == 0) cycle
       do j = 1, NBIN_HIST
          x(j) = tod%pixhist(4,i) + (j-0.5) * delta(i)
       end do
       P = real(hist(1:NBIN_HIST,i),sp) / real(nhit,sp) / delta(i)
       
       mu    =      sum(P *  x)        * delta(i)
       sigma = sqrt(sum(P * (x-mu)**2) * delta(i))
       tod%pixhist(1,i) = mu
       tod%pixhist(2,i) = sigma
       tod%pixhist(3,i) = nhit
       tod%pixhist(4,i) = mu - 4.0*sigma
       tod%pixhist(5,i) = mu + 4.0*sigma

       if (tod%myid == 0 .and. mod(i,1000)==0) then
          call int2string(i,pix_text)
          open(58,file='pixhist'//pix_text//'.dat', recl=1024)
          write(58,*) '# pixhist =', tod%pixhist(:,i)
          do j = 1, NBIN_HIST
             write(58,*) x(j), P(j)
          end do
          write(58,*)
          write(58,*) tod%pixhist(4,i), 0.
          write(58,*) tod%pixhist(4,i), 2.*maxval(P)
          write(58,*)
          write(58,*) tod%pixhist(5,i), 0.
          write(58,*) tod%pixhist(5,i), 2.*maxval(P)
       end if
    end do

    deallocate(delta, hist)
    
  end subroutine compute_tod_pixhist



end module comm_tod_pixhist_mod
