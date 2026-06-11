module comm_tod_pixhist_mod
  use comm_tod_mod
  use comm_tod_driver_mod
  implicit none

  private
  public compute_tod_pixhist

  integer(i4b), parameter :: NBIN_HIST      = 100
  integer(i4b), parameter :: MAX_ITER_HIST  = 20
  integer(i4b), parameter :: CUTOFF_N_EMPTY = 20
  real(sp),     parameter :: FRAC_THRESHOLD = 0.0001 ! Defines empty bin; fraction of most populated bin

contains

  subroutine compute_tod_pixhist(tod)
    implicit none
    class(comm_tod),                              intent(inout) :: tod

    integer(i4b) :: i, j, k, ierr, pix, npix_hist, bin, iter, nhit, hit, n_empty, j_cut1, j_cut2, ndet, det, oper, ind
    real(sp)     :: val, center, mu, sigma, x0, x1, delta0, f_threshold
    real(dp)     :: bin_tmp
    logical(lgt) :: refine
    type(comm_detdata) :: dd
    character(len=6) :: pix_text
    integer(i4b), allocatable, dimension(:,:) :: hist
    real(sp),     allocatable, dimension(:) :: delta, tmp
    real(sp),                 dimension(NBIN_HIST) :: x, P
    
    if (tod%nhorn /= 1) then
       write(*,*) 'compute_tod_pixhist does not support multi-horn data'
       stop
    end if

    ndet      = tod%ndet
    npix_hist = 12*tod%nside_pixhist**2
    oper      = get_sd_operation_code([SD_BASE,SD_TOD,SD_MASK,SD_IND])

    ! Find absolute min and max per low-res pixel
    allocate(tod%pixhist(5,0:npix_hist-1,ndet))  ! (mu, rms, nhit, min, max)
    allocate(delta(0:npix_hist-1), tmp(0:npix_hist-1))  ! (mu, rms, nhit, min, max)
    allocate(hist(0:NBIN_HIST,0:npix_hist-1))
    tod%pixhist         = 0.
    tod%pixhist(4,:,:)  =  1e30
    tod%pixhist(5,:,:)  = -1e30

    ! Loop over detectors
    do det = 1, ndet
    
       ! Set up full-mission scan data
       call init_det_data(tod, det, oper, -1, -tod%nside_pixhist, .false., dd, nonlin_level=0, spur_level=0)
       call timer%start(TOD_PIXHIST, tod%band)
       
!!$       if (tod%myid == 0 .and. det == 2) then
!!$          open(58,file='detdata.dat', recl=1024)
!!$          do j = 1, dd%ntod
!!$             write(58,*) j, dd%tod(j), dd%pix(j)
!!$          end do
!!$          close(58)
!!$       end if
       
       ! Find absolute min and max for each pixel
       do k = 1, dd%ntod
          pix = dd%pix(k)
          tod%pixhist(4,pix,det)  = min(tod%pixhist(4,pix,det), dd%tod(k))
          tod%pixhist(5,pix,det)  = max(tod%pixhist(5,pix,det), dd%tod(k))
       end do
       tmp = tod%pixhist(4,:,det)
       call mpi_allreduce(MPI_IN_PLACE, tmp,  size(tod%pixhist(4,:,det)),  MPI_REAL, MPI_MIN, tod%comm, ierr)
       tod%pixhist(4,:,det) = tmp

       tmp = tod%pixhist(5,:,det)
       call mpi_allreduce(MPI_IN_PLACE, tmp,  size(tod%pixhist(5,:,det)),  MPI_REAL, MPI_MAX, tod%comm, ierr)
       tod%pixhist(5,:,det) = tmp
       delta = (tod%pixhist(5,:,det)-tod%pixhist(4,:,det))/NBIN_HIST
       !if (tod%myid == 0) write(*,*) 'a', tod%pixhist(4:5,42527,1), delta(42527,1)
    
!!$    write(*,*) 'a', minval(delta, delta > -1e27)
!!$    call mpi_finalize(ierr)
!!$    stop

       refine = .true.
       iter   = 0
       do while (refine)
          iter = iter + 1
          !if (tod%myid == 0) write(*,*) 'Pixhist, band = ', trim(tod%freq), ', iter = ', iter
    
          ! Build histogram
          hist = 0.
          do k = 1, dd%ntod
             pix = dd%pix(k)
             if (delta(pix) == 0.) cycle ! 0 or 1 samples in current pixel
             bin_tmp = (dd%tod(k)-tod%pixhist(4,pix,det))/delta(pix)
             ! BUGFIX: this guard previously only caught large positive values; after the window
             ! limits have been refined, samples far below the new minimum give hugely negative
             ! ratios, which overflow the int() below just the same. Use abs().
             if (abs(bin_tmp) > 1d6) cycle ! Too large values causes issues with int() below.
             bin = int((dd%tod(k)-tod%pixhist(4,pix,det))/delta(pix))+1
             if (bin < 1 .or. bin > NBIN_HIST) bin = 0 ! Sample out of range; discarded
             hist(bin,pix) = hist(bin,pix) + 1
          end do
          call mpi_allreduce(MPI_IN_PLACE, hist,  size(hist),  MPI_INTEGER, MPI_SUM, tod%comm, ierr)
          
          if (iter >= MAX_ITER_HIST) exit
       
          ! Check histogram quality, and refine limits
          refine = .false.
          do i = 0, npix_hist-1
             if (delta(i) == 0.) cycle ! 0 or 1 samples in current pixel

             x0     = tod%pixhist(4,i,det)
             x1     = tod%pixhist(5,i,det)
             delta0 = delta(i)
             nhit   = sum(hist(1:NBIN_HIST,i))
             if (nhit == 0) cycle
             ! Cut bins from low side
             f_threshold = FRAC_THRESHOLD * maxval(hist(1:NBIN_HIST,i))
             j       = 1
             j_cut1  = -1
             hit     = 0
             n_empty = 0
             do while (hit < nhit/4)
                if (hist(j,i) < f_threshold) then
                   n_empty = n_empty + 1
                else
                   if (n_empty >= CUTOFF_N_EMPTY) then
                      j_cut1  = j
                      refine  = .true.
                   end if
                   n_empty = 0
                end if
                hit = hit + hist(j,i)
                j   = j+1
             end do
             if (j_cut1 /= -1) tod%pixhist(4,i,det) = x0 + (j_cut1-2)*delta(i)

             ! Cut bins from high side
             j       = NBIN_HIST
             j_cut2   = -1
             hit     = 0
             n_empty = 0
             do while (hit < nhit/4)
                if (hist(j,i) < f_threshold) then
                   n_empty = n_empty + 1
                else
                   if (n_empty >= CUTOFF_N_EMPTY) then
                      j_cut2  = j
                      refine  = .true.
                   end if
                   n_empty = 0
                end if
                hit = hit + hist(j,i)
                j   = j-1
             end do
             if (j_cut2 /= -1) tod%pixhist(5,i,det) = x0 + (j_cut2+2)*delta(i)
             
             ! Update bin width
             delta(i) = (tod%pixhist(5,i,det)-tod%pixhist(4,i,det))/NBIN_HIST
             if (delta(i) < 0.) write(*,*) "neg", i, det, j_cut1, j_cut2, x0, x1, delta0, tod%pixhist(4,i,det), tod%pixhist(5,i,det), delta(i), hist(:,i)
          end do
       end do
       call dealloc_det_data(dd)
       
       ! Define cuts and output some histograms to disk
       do i = 0, npix_hist-1
          if (delta(i) == 0.) then
             ! 0 or 1 samples in current pixel; accept anything
             tod%pixhist(1,i,det) = 0.
             tod%pixhist(2,i,det) = 0.
             tod%pixhist(3,i,det) = 0.
             tod%pixhist(4,i,det) = -1e30
             tod%pixhist(5,i,det) =  1e30
             cycle
          end if
          
          nhit = sum(hist(1:NBIN_HIST,i))
          if (nhit == 0) cycle
          do j = 1, NBIN_HIST
             x(j) = tod%pixhist(4,i,det) + (j-0.5) * delta(i)
          end do
          P = real(hist(1:NBIN_HIST,i),sp) / real(nhit,sp) / delta(i)
          
          mu    =      sum(P *  x)        * delta(i)
          sigma = sqrt(sum(P * (x-mu)**2) * delta(i))
          tod%pixhist(1,i,det) = mu
          tod%pixhist(2,i,det) = sigma
          tod%pixhist(3,i,det) = nhit
          tod%pixhist(4,i,det) = mu - 4.0*sigma
          tod%pixhist(5,i,det) = mu + 4.0*sigma
          
          if (.false. .and. tod%myid == 0 .and. mod(i,1000)==0 .and. det == 10) then
             call int2string(i,pix_text)
             open(58,file='pixhist'//pix_text//'.dat', recl=1024)
             write(58,*) '# pixhist =', tod%pixhist(:,i,det)
             do j = 1, NBIN_HIST
                write(58,*) x(j), P(j)
             end do
             write(58,*)
             write(58,*) tod%pixhist(4,i,det), 0.
             write(58,*) tod%pixhist(4,i,det), 2.*maxval(P)
             write(58,*)
             write(58,*) tod%pixhist(5,i,det), 0.
             write(58,*) tod%pixhist(5,i,det), 2.*maxval(P)
             close(58)
          end if
       end do

       call timer%stop(TOD_PIXHIST, tod%band)
       if (tod%myid == 0) write(*,fmt='(a,a,a,i4,a,i6)') '    --> Pixhist, band = ', trim(tod%freq), ', det = ', det, ', numiter = ', iter
    end do
    
    deallocate(delta, hist, tmp)
    
  end subroutine compute_tod_pixhist



  subroutine compute_tod_pixhist2(tod)
    implicit none
    class(comm_tod),                              intent(inout) :: tod

    integer(i4b) :: i, j, k, ierr, pix, npix_hist, q, bin, iter, nhit, hit, n_empty, j_cut1, j_cut2, ndet, det, oper, ind
    real(sp)     :: val, center, mu, sigma, x0, x1, delta0, f_threshold
    logical(lgt) :: refine
    type(comm_scandata) :: sd
    character(len=6) :: pix_text
    integer(i4b), allocatable, dimension(:,:,:) :: hist
    real(sp), allocatable, dimension(:,:) :: delta, buff
    real(sp),              dimension(NBIN_HIST) :: x, P
    
    if (tod%nhorn /= 1) then
       write(*,*) 'compute_tod_pixhist does not support multi-horn data'
       stop
    end if

    if (tod%myid == 0) write(*,fmt='(a,a)') '    --> Pixhist, band = ', trim(tod%freq)

    ndet      = tod%ndet
    npix_hist = 12*tod%nside_pixhist**2
    q         = (tod%nside / tod%nside_pixhist)**2
    oper      = get_sd_operation_code([SD_BASE,SD_TOD])

!!$    if (tod%myid == 0) then
!!$       open(58,file='tod.dat')
!!$       do i = 1, tod%nscan
!!$          if (.not. any(tod%scans(i)%d%accept)) cycle
!!$          call init_scan_data_singlehorn(sd, tod, i, map_sky, map_gain, procmask, procmask2, skip_zodi=.true.,handle_=handle)
!!$          do j = 1, tod%ndet
!!$             if (.not. tod%scans(i)%d(j)%accept) cycle
!!$             do k = 1, sd%ntod
!!$                if (iand(sd%flag(k,j), tod%flag0) .ne. 0) cycle
!!$                write(58,*) i, j, k, sd%tod(k,j)
!!$             end do
!!$          end do
!!$          call dealloc_scan_data(sd)
!!$       end do
!!$       close(58)
!!$    end if
!!$    call mpi_finalize(ierr)
!!$    stop
    
    ! Find absolute min and max per low-res pixel
    allocate(tod%pixhist(5,0:npix_hist-1,ndet))  ! (mu, rms, nhit, min, max)
    allocate(delta(0:npix_hist-1,ndet))  ! (mu, rms, nhit, min, max)
    allocate(hist(0:NBIN_HIST,0:npix_hist-1,ndet))
    tod%pixhist         = 0.
    tod%pixhist(4,:,:)  =  1e30
    tod%pixhist(5,:,:)  = -1e30
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call init_scan_data(tod, i, oper, -1, sd)
       call timer%start(TOD_PIXHIST, tod%band)
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          do k = 1, sd%ntod
             if (iand(sd%flag(k,j), tod%flag0) .ne. 0) cycle
             call ring2nest(tod%nside, sd%pix(k,j,1), pix)
             pix = pix / q
             tod%pixhist(4,pix,j)  = min(tod%pixhist(4,pix,j), sd%tod(k,j))
             tod%pixhist(5,pix,j)  = max(tod%pixhist(5,pix,j), sd%tod(k,j))
          end do
       end do
       call dealloc_scan_data(sd)
       call timer%stop(TOD_PIXHIST, tod%band)
    end do
    allocate(buff(0:npix_hist-1, ndet))
    buff = tod%pixhist(4,:,:)
    call mpi_allreduce(MPI_IN_PLACE, buff,  size(buff),  MPI_REAL, MPI_MIN, tod%comm, ierr)
    tod%pixhist(4,:,:) = buff
    buff = tod%pixhist(5,:,:)
    call mpi_allreduce(MPI_IN_PLACE, buff,  size(buff),  MPI_REAL, MPI_MAX, tod%comm, ierr)
    tod%pixhist(5,:,:) = buff
    deallocate(buff)
    delta = (tod%pixhist(5,:,:)-tod%pixhist(4,:,:))/NBIN_HIST
    call timer%stop(TOD_PIXHIST, tod%band)

    !if (tod%myid == 0) write(*,*) 'a', tod%pixhist(4:5,42527,1), delta(42527,1)
    
!!$    write(*,*) 'a', minval(delta, delta > -1e27)
!!$    call mpi_finalize(ierr)
!!$    stop


    refine = .true.
    iter   = 0
    do while (refine)
       iter = iter + 1
       !if (tod%myid == 0) write(*,*) 'Pixhist, band = ', trim(tod%freq), ', iter = ', iter
    
       ! Build histogram
       hist = 0.
       do i = 1, tod%nscan
          if (.not. any(tod%scans(i)%d%accept)) cycle
          call init_scan_data(tod, i, oper, -1, sd)
          call timer%start(TOD_PIXHIST, tod%band)
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             do k = 1, sd%ntod
                if (iand(sd%flag(k,j), tod%flag0) .ne. 0) cycle
                call ring2nest(tod%nside, sd%pix(k,j,1), pix)
                pix = pix / q
                if (delta(pix,j) == 0.) cycle ! 0 or 1 samples in current pixel
!!$                if (sd%tod(k,j) /= sd%tod(k,j) .or.tod%pixhist(4,pix,j) /= tod%pixhist(4,pix,j) .or. delta(pix,j) /= delta(pix,j)) then
                !write(*,*) 'nan', tod%scanid(i), j, k, pix, sd%tod(k,j), tod%pixhist(4,pix,j), tod%pixhist(5,pix,j), delta(pix,j) 
!!$                end if
!!$                if (abs(sd%tod(k,j)) > 1e30 .or.abs(tod%pixhist(4,pix,j)) > 1e30  .or. abs(delta(pix,j)) < 1e-10) then
!!$                   write(*,*) 'ovf', tod%scanid(i), j, k, sd%tod(k,j), tod%pixhist(4,pix,j), delta(pix,j) 
!!$                end if

                bin = int((sd%tod(k,j)-tod%pixhist(4,pix,j))/delta(pix,j))+1
                if (bin < 1 .or. bin > NBIN_HIST) bin = 0 ! Sample out of range; discarded
                hist(bin,pix,j) = hist(bin,pix,j) + 1

             end do
          end do
          call dealloc_scan_data(sd)
          call timer%stop(TOD_PIXHIST, tod%band)
       end do
       call timer%start(TOD_PIXHIST, tod%band)
       call mpi_allreduce(MPI_IN_PLACE, hist,  size(hist),  MPI_INTEGER, MPI_SUM, tod%comm, ierr)

       if (iter >= MAX_ITER_HIST) exit
       
       ! Check histogram quality, and refine limits
       refine = .false.
       do det = 1, ndet
          do i = 0, npix_hist-1
             if (delta(i,det) == 0.) cycle ! 0 or 1 samples in current pixel


             x0     = tod%pixhist(4,i,det)
             x1     = tod%pixhist(5,i,det)
             delta0 = delta(i,det)
             nhit   = sum(hist(1:NBIN_HIST,i,det))
             if (nhit == 0) cycle
             ! Cut bins from low side
             f_threshold = FRAC_THRESHOLD * maxval(hist(1:NBIN_HIST,i,det))
             j       = 1
             j_cut1  = -1
             hit     = 0
             n_empty = 0
             do while (hit < nhit/4)
                if (hist(j,i,det) < f_threshold) then
                   n_empty = n_empty + 1
                else
                   if (n_empty >= CUTOFF_N_EMPTY) then
                      j_cut1  = j
                      refine  = .true.
                   end if
                   n_empty = 0
                end if
                hit = hit + hist(j,i,det)
                j   = j+1
             end do
             if (j_cut1 /= -1) tod%pixhist(4,i,det) = x0 + (j_cut1-2)*delta(i,det)

             ! Cut bins from high side
             j       = NBIN_HIST
             j_cut2   = -1
             hit     = 0
             n_empty = 0
             do while (hit < nhit/4)
                if (hist(j,i,det) < f_threshold) then
                   n_empty = n_empty + 1
                else
                   if (n_empty >= CUTOFF_N_EMPTY) then
                      j_cut2  = j
                      refine  = .true.
                   end if
                   n_empty = 0
                end if
                hit = hit + hist(j,i,det)
                j   = j-1
             end do
             if (j_cut2 /= -1) tod%pixhist(5,i,det) = x0 + (j_cut2+2)*delta(i,det)
             
             ! Update bin width
             delta(i,det) = (tod%pixhist(5,i,det)-tod%pixhist(4,i,det))/NBIN_HIST
             if (delta(i,det) < 0.) write(*,*) "neg", i, det, j_cut1, j_cut2, x0, x1, delta0, tod%pixhist(4,i,det), tod%pixhist(5,i,det), delta(i,det), hist(:,i,det)
             !if (tod%myid == 0) write(*,*) 'b', iter, i, tod%pixhist(4:5,42527,1), delta(42527,1)
          end do
       end do
       call timer%stop(TOD_PIXHIST, tod%band)
    end do

    
    ! Define cuts and output some histograms to disk
    call timer%start(TOD_PIXHIST, tod%band)
    do det = 1, ndet
       do i = 0, npix_hist-1
          if (delta(i,det) == 0.) then
             ! 0 or 1 samples in current pixel; accept anything
             tod%pixhist(1,i,det) = 0.
             tod%pixhist(2,i,det) = 0.
             tod%pixhist(3,i,det) = 0.
             tod%pixhist(4,i,det) = -1e30
             tod%pixhist(5,i,det) =  1e30
             cycle
          end if
          
          nhit = sum(hist(1:NBIN_HIST,i,det))
          if (nhit == 0) cycle
          do j = 1, NBIN_HIST
             x(j) = tod%pixhist(4,i,det) + (j-0.5) * delta(i,det)
          end do
          P = real(hist(1:NBIN_HIST,i,det),sp) / real(nhit,sp) / delta(i,det)
          
          mu    =      sum(P *  x)        * delta(i,det)
          sigma = sqrt(sum(P * (x-mu)**2) * delta(i,det))
          tod%pixhist(1,i,det) = mu
          tod%pixhist(2,i,det) = sigma
          tod%pixhist(3,i,det) = nhit
          tod%pixhist(4,i,det) = mu - 4.0*sigma
          tod%pixhist(5,i,det) = mu + 4.0*sigma
          
!!$          if (tod%myid == 0 .and. mod(i,1000)==0 .and. det == 1) then
!!$             call int2string(i,pix_text)
!!$             open(58,file='pixhist'//pix_text//'.dat', recl=1024)
!!$             write(58,*) '# pixhist =', tod%pixhist(:,i,det)
!!$             do j = 1, NBIN_HIST
!!$                write(58,*) x(j), P(j)
!!$             end do
!!$             write(58,*)
!!$             write(58,*) tod%pixhist(4,i,det), 0.
!!$             write(58,*) tod%pixhist(4,i,det), 2.*maxval(P)
!!$             write(58,*)
!!$             write(58,*) tod%pixhist(5,i,det), 0.
!!$             write(58,*) tod%pixhist(5,i,det), 2.*maxval(P)
!!$          end if
       end do
    end do
    
    if (tod%myid == 0) write(*,fmt='(a,a,a,i6)') '    --> Pixhist, band = ', trim(tod%freq), ', numiter = ', iter

    deallocate(delta, hist)
    call timer%stop(TOD_PIXHIST, tod%band)
    
  end subroutine compute_tod_pixhist2



end module comm_tod_pixhist_mod
