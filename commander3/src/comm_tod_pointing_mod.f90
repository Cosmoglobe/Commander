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
  use comm_utils
  implicit none
   
contains

  ! Sky signal template for single detector
  subroutine compute_scan_pix2ind(tod, sd)
     implicit none
     class(comm_tod),      intent(in)    :: tod
     class(comm_scandata), intent(inout) :: sd
 
     integer(i4b) :: i, j, k, d, h, p, nmaps, scan, ndet

     scan  = sd%scan
     nmaps = tod%nmaps
     ndet  = sd%ndet

     do j = 1, ndet ! Loop over detectors
        d = sd%det(j)
        
        if (.not. tod%scans(scan)%d(d)%accept) then
           sd%ind(:,j,:) = -999
           cycle
        end if

        do i = 1, tod%scans(scan)%ntod ! Loop over time samples
           do h = 1, sd%nhorn          ! Loop over horns
              sd%ind(i,j,h) = tod%pixcache%pix2ind(sd%pix(i,j,h))
              if (sd%pix(i,j,h) == -1) then
                 write(*,*) "missing pixel", sd%pix(i,j,h), tod%pixcache%ind2pix(tod%pixcache%nobs-2:tod%pixcache%nobs)
              end if
           end do
        end do
     end do
   end subroutine compute_scan_pix2ind

  
  ! Sky signal template for single detector
  subroutine project_sky(tod, sd)
     implicit none
     class(comm_tod),      intent(in)             :: tod
     class(comm_scandata), intent(inout)          :: sd
 
     integer(i4b) :: i, j, d_cache, k, d, h, hp, p, nmaps, scan, ndet
     logical(lgt) :: do_gain, do_bp
     real(sp)     :: s, eff

     do_gain = btest(sd%oper,SD_GAIN) .and. allocated(tod%pixcache%map_gain)
     do_bp   = btest(sd%oper,SD_BP)
     
     ! if (tod%myid==0) write(*,*) 'd', tod%pixcache%map_sky(:,109952,1,1)
     
     ! s = T + eff*(Q * cos(2*psi) + U * sin(2*psi))
     ! T - temperature; Q, U - Stoke's parameters
     scan  = sd%scan
     nmaps = tod%nmaps
     ndet  = sd%ndet

     do j = 1, ndet ! Loop over detectors
        d = sd%det(j)
        d_cache = d; if (tod%equal_det_bp_beam) d_cache = 0

        if (.not. tod%scans(scan)%d(d)%accept) then
           sd%s_sky(:,j,:,:) = 0.
           if (allocated(tod%pixcache%map_gain)) sd%s_gain(:,j,:) = 0.
           cycle
        end if

        eff = 1.0 !tod%pol_sign(d) * tod%pol_eff(d)
        do h = 1, sd%nhorn          ! Loop over horns
           hp = h; if (sd%nhorn == 1) hp = 0
           do i = 1, tod%scans(scan)%ntod ! Loop over time samples
              p = sd%ind(i,j,h)
              if ((sd%psi(i,j,h) > tod%npsi)) then
                 write(*,*) 'Polarization angle is wrong', d, tod%scanid(scan), sd%psi(i,j,h)
                 cycle
              end if
              do k = 1, sd%nbp ! Loop over bandpass models
                 if (nmaps == 3) then
                    sd%s_sky(i,j,hp,k) = tod%pixcache%map_sky(1,p,d_cache,k) + &
                         & eff*(tod%pixcache%map_sky(2,p,d_cache,k) * tod%pixcache%cos2psi(sd%psi(i,j,h)) + &
                         &      tod%pixcache%map_sky(3,p,d_cache,k) * tod%pixcache%sin2psi(sd%psi(i,j,h)))
                   ! write(*,*) j, i, p, sd%s_sky(i,j,hp,k), eff, tod%pixcache%map_sky(1:3,p,d,k), tod%pixcache%cos2psi(sd%psi(i,d,h)), tod%pixcache%sin2psi(sd%psi(i,d,h))
                    if (do_gain .and. k == 1) then
                       sd%s_gain(i,j,hp) = tod%pixcache%map_gain(1,p,d_cache) + &
                            & eff*(tod%pixcache%map_gain(2,p,d_cache) * tod%pixcache%cos2psi(sd%psi(i,j,h)) + &
                            &      tod%pixcache%map_gain(3,p,d_cache) * tod%pixcache%sin2psi(sd%psi(i,j,h)))
                    end if
                 else 
                    sd%s_sky(i,j,hp,k) = tod%pixcache%map_sky(1,p,d_cache,k)
                    if (do_gain .and. k == 1) sd%s_gain(i,j,hp) = tod%pixcache%map_gain(1,p,d_cache)
                 end if

                 if (do_bp) then
                    if (nmaps == 3) then
                       s = tod%pixcache%map_sky(1,p,0,k) + eff * &
                            & (tod%pixcache%map_sky(2,p,0,k) * tod%pixcache%cos2psi(sd%psi(i,j,h)) + &
                            &  tod%pixcache%map_sky(3,p,0,k) * tod%pixcache%sin2psi(sd%psi(i,j,h)))
                    else if (nmaps == 1) then
                       s = tod%pixcache%map_sky(1,p,0,k)
                    end if
                    if (allocated(sd%s_bp)) sd%s_bp(i,j,hp,k) = sd%s_sky(i,j,hp,k) - s
                 end if
              end do
           end do
        end do
     end do
   end subroutine project_sky

   ! Sky signal template for single detector
   subroutine project_mask(tod, bitmask0, sd)
     implicit none
     class(comm_tod),      intent(in)    :: tod
     integer(i4b),         intent(in)    :: bitmask0
     class(comm_scandata), intent(inout) :: sd
 
     integer(i4b) :: i, j, k, d, h, p, nmaps, scan, ndet

     scan        = sd%scan
     ndet        = sd%ndet
     sd%bitmask0 = bitmask0
     
     do j = 1, ndet ! Loop over detectors
        d = sd%det(j)
         
        if (.not. tod%scans(scan)%d(d)%accept) then
           sd%mask(:,j)  = 0.
           cycle
        end if

        do h = 1, sd%nhorn          ! Loop over horns
           do i = 1, tod%scans(scan)%ntod ! Loop over time samples
              p = sd%ind(i,j,h)
              if (h == 1) then
                 if (iand(sd%flag(i,j), tod%flag0) .ne. 0) then
                    sd%mask(i,j) = 0.
                 else
                    sd%mask(i,j) = 1.
                 end if
              end if
              if (sd%bitmask0 >= 0 .and. sd%mask(i,j) == 1.) then
                 if (btest(tod%pixcache%bitmask(p), sd%bitmask0)) sd%mask(i,j) = 0.
              end if
           end do
        end do

     end do
   end subroutine project_mask
    
 end module comm_tod_pointing_mod
