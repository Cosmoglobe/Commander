module comm_tod_pointing_mod
  use comm_tod_mod
  use comm_map_mod
  use comm_utils
  implicit none


contains

  ! Sky signal template
  subroutine project_sky(tod, map, pix_in, psi_in, flag, pmask, scan_id, &
       & s_sky, tmask, s_bp)
    implicit none
    class(comm_tod),                          intent(in)  :: tod
    integer(i4b),        dimension(0:),       intent(in)  :: pmask
    real(sp),            dimension(1:,1:,0:), intent(in)  :: map
    !type(shared_2d_sp),  dimension(0:),     intent(in)  :: map
    integer(i4b),        dimension(:,:,:),      intent(in)  :: pix_in, psi_in
    integer(i4b),        dimension(:,:),      intent(in)  :: flag
    integer(i4b),                             intent(in)  :: scan_id
    real(sp),            dimension(:,:),      intent(out) :: s_sky, tmask
    real(sp),            dimension(:,:),      intent(out), optional :: s_bp

    integer(i4b) :: i, p, det
    real(sp)     :: s
    integer(i4b), allocatable, dimension(:,:)   :: pix, psi

    if(size(pix, 2) /= 1 .or. size(psi,2) /= 1) then
      write(*,*) "Call to project sky with nhorn /= 1. You probably want project_sky_differential."
      return
    end if

    pix = pix_in(:,:,1)
    psi = psi_in(:,:,1)
    ! s = T + Q * cos(2 * psi) + U * sin(2 * psi)
    ! T - temperature; Q, U - Stoke's parameters
    do det = 1, tod%ndet
       if (.not. tod%scans(scan_id)%d(det)%accept) then
          s_sky(:,det) = 0.d0
          tmask(:,det) = 0.d0
          cycle
       end if
       do i = 1, tod%scans(scan_id)%ntod
          p = tod%pix2ind(pix(i,det))
          s_sky(i,det) = map(1,p,det) + &
                       & map(2,p,det) * tod%cos2psi(psi(i,det)) + &
                       & map(3,p,det) * tod%sin2psi(psi(i,det))
!!$          s_sky(i,det) = map(det)%a(1,pix(i,det)+1) + &
!!$                       & map(det)%a(2,pix(i,det)+1) * tod%cos2psi(psi(i,det)) + &
!!$                       & map(det)%a(3,pix(i,det)+1) * tod%sin2psi(psi(i,det))
!          if (s_sky(i,det) /= s_sky(i,det)) then
!             write(*,*) det, i, map(det)%a(:,pix(i,det)+1), tod%cos2psi(psi(i,det)), tod%sin2psi(psi(i,det))
!             stop
!          end if

          tmask(i,det) = pmask(pix(i,det)) 
          if (iand(flag(i,det),tod%flag0) .ne. 0) tmask(i,det) = 0.
       end do
    end do

    if (present(s_bp)) then
       do det = 1, tod%ndet
          if (.not. tod%scans(scan_id)%d(det)%accept) then
             s_bp(:,det) = 0.d0
             cycle
          end if
          do i = 1, tod%scans(scan_id)%ntod
             p = tod%pix2ind(pix(i,det))
             s =    map(1,p,0) + &
                  & map(2,p,0) * tod%cos2psi(psi(i,det)) + &
                  & map(3,p,0) * tod%sin2psi(psi(i,det))
!!$             s =    map(0)%a(1,pix(i,det)+1) + &
!!$                  & map(0)%a(2,pix(i,det)+1) * tod%cos2psi(psi(i,det)) + &
!!$                  & map(0)%a(3,pix(i,det)+1) * tod%sin2psi(psi(i,det))
             s_bp(i,det)  = s_sky(i,det) - s
          end do
       end do
    end if

  end subroutine project_sky


    ! Sky signal template
  subroutine project_sky_differential(tod, map, pix, psi, x_im, dx_im, flag, pmask, scan_id,&
       & s_sky, tmask, s_bp)
    implicit none
    class(comm_tod),                          intent(in)  :: tod
    integer(i4b),        dimension(0:),       intent(in)  :: pmask
    real(sp),            dimension(1:,1:,0:), intent(in)  :: map
    real(sp),                                 intent(in)  :: x_im, dx_im
    !type(shared_2d_sp),  dimension(0:),     intent(in)  :: map
    integer(i4b),        dimension(:,:,:),      intent(in)  :: pix, psi
    integer(i4b),        dimension(:,:),      intent(in)  :: flag
    integer(i4b),                             intent(in)  :: scan_id
    real(sp),            dimension(:,:),      intent(out) :: s_sky, tmask
    real(sp),            dimension(:,:),      intent(out), optional :: s_bp

    integer(i4b) :: i, j, lpoint, rpoint


    do i = 1, tod%ndet
      if (.not. tod%scans(scan_id)%d(i)%accept) then
        s_sky(:,i) = 0.d0
        tmask(:,i) = 0.d0
        cycle
      end if

      do j = 1, tod%scans(scan_id)%ntod
        lpoint = tod%pix2ind(pix(j,i,1))
        rpoint = tod%pix2ind(pix(j,i,2))
        !todo: add correct coefficients here when they are loaded in

        s_sky(j,i) = (map(1,lpoint,i) + &
                    & map(2,lpoint,i) * tod%cos2psi(psi(j,i,1)) + &
                    & map(3,lpoint,i) * tod%sin2psi(psi(j,i,1)) - &
                    & map(1,rpoint,i) - &
                    & map(2,rpoint,i) * tod%cos2psi(psi(j,i,2)) - &
                    & map(3,rpoint,i) * tod%sin2psi(psi(j,i,2)))
        ! The gain imbalance parameters x are different for each radiometer.
        ! Could be another parameter that would be fit. How would that be
        ! included? Conversely, could just use the fixed values from Bennett et
        ! al. (2013). I believe they're supposed to be constant throughout the
        ! experiment.
        ! d1 = (1+x1)*[T(pA) + P(pA,gA) + S(pA)] 
        !     -(1-x1)*[T(pB) + P(pB,gB) + S(pB)]

        ! Actually, is it possible and / or necessary to do this for WMAP? You
        ! can create sky templates for each data stream in the differencing
        ! assembly, would this make sense?
        s_sky(j,i) = ((1+x_im) * (map(1,lpoint,i) + &
                   &              map(2,lpoint,i) * tod%cos2psi(psi(j,i,1)) + &
                   &              map(3,lpoint,i) * tod%sin2psi(psi(j,i,1)) + &
                   &              map(4,lpoint,i))                          - &
                   &  (1-x_im) * (map(1,rpoint,i) + &
                   &              map(2,rpoint,i) * tod%cos2psi(psi(j,i,2)) + &
                   &              map(3,rpoint,i) * tod%sin2psi(psi(j,i,2)) + &
                   &              map(4,rpoint,i)))                          

        if (iand(flag(j,i),tod%flag0) .ne. 0) tmask(j,i) = 0.
        if(pmask(lpoint) .or. pmask(rpoint)) then
          tmask(j,i) = 1
        end if
      end do
    end do

  end subroutine project_sky_differential

end module comm_tod_pointing_mod
