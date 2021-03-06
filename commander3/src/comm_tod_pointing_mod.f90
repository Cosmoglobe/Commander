!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_tod_pointing_mod
  use comm_tod_mod
  use comm_map_mod
  use comm_utils
  implicit none


contains

  ! Sky signal template
  subroutine project_sky(tod, map, pix, psi, flag, pmask, scan_id, &
       & s_sky, tmask, s_bp)
    implicit none
    class(comm_tod),                          intent(in)  :: tod
    integer(i4b),        dimension(0:),       intent(in)  :: pmask
    real(sp),            dimension(1:,1:,0:), intent(in)  :: map
    !type(shared_2d_sp),  dimension(0:),     intent(in)  :: map
    integer(i4b),        dimension(:,:),      intent(in)  :: pix, psi
    integer(i4b),        dimension(:,:),      intent(in)  :: flag
    integer(i4b),                             intent(in)  :: scan_id
    real(sp),            dimension(:,:),      intent(out) :: s_sky, tmask
    real(sp),            dimension(:,:),      intent(out), optional :: s_bp

    integer(i4b) :: i, p, det
    real(sp)     :: s

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


end module comm_tod_pointing_mod
