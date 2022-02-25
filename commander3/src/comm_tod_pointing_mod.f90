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
      class(comm_tod),                   intent(in)             :: tod
      real(sp),     dimension(0:),       intent(in)             :: pmask
      real(sp),     dimension(1:,1:,0:), intent(in)             :: map
      integer(i4b), dimension(:,:),      intent(in)             :: pix, psi
      integer(i4b), dimension(:,:),      intent(in)             :: flag
      integer(i4b),                      intent(in)             :: scan_id
      real(sp),     dimension(:,:),      intent(out)            :: s_sky, tmask
      real(sp),     dimension(:,:),      intent(out), optional  :: s_bp

      integer(i4b)                                      :: i, p, det, nmap
      real(sp)                                          :: s

      ! s = T + Q * cos(2 * psi) + U * sin(2 * psi)
      ! T - temperature; Q, U - Stoke's parameters
!      if (tod%myid == 78) write(*,*) 'c611', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

      nmap = SIZE(map, 1)
      do det = 1, tod%ndet
!         if (tod%myid == 78) write(*,*) 'c6111', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
         if (.not. tod%scans(scan_id)%d(det)%accept) then
            s_sky(:, det) = 0.d0
            tmask(:, det) = 0.d0
            cycle
         end if
         !if (tod%myid == 78) write(*,*) 'c6112', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
         do i = 1, tod%scans(scan_id)%ntod
            p = tod%pix2ind(pix(i,det))
            !if (tod%myid == 78 .and. p == 7863) write(*,*) 'c61121', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires, i, p
            
            if (nmap == 3) then
                s_sky(i,det) = map(1,p,det) + &
                         & map(2,p,det) * tod%cos2psi(psi(i,det)) + &
                         & map(3,p,det) * tod%sin2psi(psi(i,det))
            else if (nmap == 1) then
                s_sky(i,det) = map(1,p,det)  ! HFI 545 thing
            end if  


            !if (tod%myid == 78 .and. p == 7863) write(*,*) 'c61122', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires, i, p
            tmask(i,det) = pmask(pix(i,det))
            !if (tod%myid == 78 .and. p == 7863) write(*,*) 'c61123', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires, i, p
            if (iand(flag(i,det), tod%flag0) .ne. 0) tmask(i,det) = 0.
            !if (tod%myid == 78 .and. p == 7863) write(*,*) 'c61124', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires, i, p
         end do
         !if (tod%myid == 78) write(*,*) 'c6113', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
      end do
      !if (tod%myid == 78) write(*,*) 'c612', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

      if (present(s_bp)) then
         do det = 1, tod%ndet
            if (.not. tod%scans(scan_id)%d(det)%accept) then
               s_bp(:,det) = 0.d0
               cycle
            end if
            do i = 1, tod%scans(scan_id)%ntod
               p = tod%pix2ind(pix(i,det))
               if (nmap == 3) then
                   s = map(1,p,0) + &
                      & map(2,p,0) * tod%cos2psi(psi(i,det)) + &
                      & map(3,p,0) * tod%sin2psi(psi(i,det))
               else if (nmap == 1) then
                   s = map(1,p,0)
               end if
               s_bp(i,det) = s_sky(i,det) - s
            end do
         end do
      end if
      !if (tod%myid == 78) write(*,*) 'c613', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

   end subroutine project_sky

   ! Sky signal template
   subroutine project_sky_differential(tod, map, pix, psi, flag, pmask, scan_id,&
        & s_skyA, s_skyB, tmask, s_bpA, s_bpB)
      implicit none
      !class(comm_tod), intent(in)  :: tod
      ! It is only inout for simulating data
      class(comm_tod),                      intent(inout)  :: tod
      real(sp),        dimension(0:),       intent(in)     :: pmask
      real(sp),        dimension(1:,1:,0:), intent(in)     :: map
      integer(i4b),    dimension(:,:),      intent(in)     :: pix, psi
      integer(i4b),    dimension(:),        intent(in)     :: flag
      integer(i4b),                         intent(in)     :: scan_id
      real(sp),        dimension(:,:),      intent(out)    :: s_skyA, s_skyB, tmask
      real(sp),        dimension(:,:),      intent(out), optional :: s_bpA, s_bpB

      integer(i4b) :: i, j, lpoint, rpoint, det
      real(sp)     :: sA, sB, tr, pr, tl, pl
      real(sp), dimension(4) :: sgn=[1., 1., -1., -1.]


      if (any(.not. tod%scans(scan_id)%d(:)%accept)) then
         s_skyA = 0.
         s_skyB = 0.
         tmask  = 0.
         if (present(s_bpA) .and. present(s_bpB)) then
            s_bpA = 0.
            s_bpB = 0.
         end if
         return
      end if

      do i = 1, tod%ndet
         do j = 1, tod%scans(scan_id)%ntod
            lpoint = tod%pix2ind(pix(j, 1))
            rpoint = tod%pix2ind(pix(j, 2))
            ! The gain imbalance parameters x are different for each radiometer.
            ! d13 = (1+x1)*[T(pA) + P(pA,gA) + S(pA)]
            !      -(1-x1)*[T(pB) + P(pB,gB) + S(pB)]
            ! We need to make sure that the imbalance parameters are redundant,
            ! i.e., d13 and d14 have the same model,
            ! d14 = (1+x1)*[T(pA) + P(pA,gA) + S(pA)]
            !      -(1-x1)*[T(pB) + P(pB,gB) + S(pB)]
            ! but d23 and d24 have different models,
            ! i.e., d13 and d14 have the same model,
            ! d23 = (1+x2)*[T(pA) - P(pA,gA) - S(pA)]
            !      -(1-x2)*[T(pB) - P(pB,gB) - S(pB)]

            s_skyA(j, i) = map(1, lpoint, i) + &
                       &  sgn(i)*( &
                       &  map(2, lpoint, i)*tod%cos2psi(psi(j, 1)) + &
                       &  map(3, lpoint, i)*tod%sin2psi(psi(j, 1))) 
            s_skyB(j, i) = map(1, rpoint, i) + &
                       &  sgn(i) *( &
                       &  map(2, rpoint, i)*tod%cos2psi(psi(j, 2)) + &
                       &  map(3, rpoint, i)*tod%sin2psi(psi(j, 2)))
                    
            if (i == 1) then
               if (iand(flag(j), tod%flag0) .ne. 0) then
                  tmask(j, :) = 0.
               else
                  tmask(j, :) = pmask(pix(j, 1))*pmask(pix(j,2))
               end if
            end if
         end do
      end do

      if (present(s_bpA) .and. present(s_bpB)) then
         do i = 1, tod%scans(scan_id)%ntod
            lpoint = tod%pix2ind(pix(i, 1))
            rpoint = tod%pix2ind(pix(i, 2))
            tl     = map(1, lpoint, 0) 
            tr     = map(1, rpoint, 0) 
            pl     = map(2, lpoint, 0)*tod%cos2psi(psi(i, 1)) + &
                  &  map(3, lpoint, 0)*tod%sin2psi(psi(i, 1))
            pr     = map(2, rpoint, 0)*tod%cos2psi(psi(i, 2)) + &
                  &  map(3, rpoint, 0)*tod%sin2psi(psi(i, 2))

            do det = 1, tod%ndet      
               s_bpA(i, det) = s_skyA(i, det) - (tl + sgn(det)*pl) 
               s_bpB(i, det) = s_skyB(i, det) - (tr + sgn(det)*pr) 
            end do
         end do
      end if

   end subroutine project_sky_differential

end module comm_tod_pointing_mod
