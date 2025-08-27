!===============================================================================
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
submodule (comm_tod_mod) comm_tod_pixcache_smod
contains

  module function constructor_tod_pixcache(nside, nside_sl, nmaps, fullsky) result(c)
    implicit none
    integer(i4b),             intent(in) :: nside, nside_sl, nmaps
    logical(lgt),             intent(in) :: fullsky
    class(comm_tod_pixcache), pointer    :: c

    allocate(c)
    c%nside        = nside
    c%nside_sl     = nside_sl
    c%nmaps        = nmaps
    c%fullsky      = fullsky
    c%nmax         = 1
    c%nobs         = 0
    allocate(c%ind2pix(c%nmax))
    
  end function constructor_tod_pixcache

  module function pix2ind(self, pix, flag_missing) result(ind)
    implicit none
    class(comm_tod_pixcache), intent(in)           :: self
    integer(i4b),             intent(in)           :: pix
    logical(lgt),             intent(in), optional :: flag_missing
    integer(i4b)                                   :: ind
    logical(lgt) :: flag_missing_
    if (self%fullsky) then
       ind = pix+1
    else
       if (self%nobs == 0) then
          ind = 0
          flag_missing_ = .true.; if (present(flag_missing)) flag_missing_ = flag_missing
          if (flag_missing_) ind = -1
       else
          ind = locate(self%ind2pix(1:self%nobs), pix)
          flag_missing_ = .true.; if (present(flag_missing)) flag_missing_ = flag_missing
          if (flag_missing_) then
             if (.not. self%ind2pix(ind) == pix) then
                write(*,*) "pix2ind", pix, ind, self%nobs, self%ind2pix(ind-1:ind+1)
                ind = -1
             end if
          end if
       end if
    end if
  end function pix2ind
  
  module subroutine add_pixels(self, pix)
    implicit none
    class(comm_tod_pixcache),               intent(inout) :: self
    integer(i4b),             dimension(:), intent(in)    :: pix

    integer(i4b) :: i, j, k, p, ind, n
    real(dp) :: t1, t2
    integer(i4b), allocatable, dimension(:) :: newpix, oldpix

    if (self%fullsky) return
!!$    call wall_time(t1)

    !if (any(pix == 46051400)) write(*,*) 'hit'
    
    allocate(newpix(size(pix)), oldpix(self%nobs))
    newpix = pix
    oldpix = self%ind2pix(1:self%nobs)
    call QuickSort_int(newpix)

    ! Count unique pixels
    n = 1
    do p = 2, size(newpix)
       if (newpix(p) /= newpix(p-1)) n = n+1
    end do
    call self%expand_storage(n)
!    write(*,*) "new", n, self%nmax
    
    ! Merge new sorted list into existing
    i = 1; j = 1; k = 1
    do while (i <= self%nobs .and. j <= size(newpix))
!!$       if (any(pix == 46051400) .and. self%nobs==22036) then
!!$          write(*,*) k, i, j, oldpix(i), newpix(j)
!!$       end if
       if (oldpix(i) == newpix(j)) then
          self%ind2pix(k) = oldpix(i)
          i = i+1; j = j+1; k = k+1
       else if (oldpix(i) < newpix(j)) then
          self%ind2pix(k) = oldpix(i)
          i = i+1; k = k+1
       else if (oldpix(i) > newpix(j)) then
          self%ind2pix(k) = newpix(j)
          j = j+1; k = k+1
       end if
       ! Skip duplicates
       if (j > size(newpix)) exit
       if (j == 1) cycle
       do while (newpix(j) == newpix(j-1))
          j = j+1
          if (j > size(newpix)) exit
       end do
    end do

!!$    if(any(pix == 46051400)) then
!!$       write(*,*) 'check leftover', i, self%nobs
!!$       write(*,*) 'check leftover', j, size(newpix)
!!$
!!$    end if
    
    if (i <= self%nobs) then
       ! Append left-over pixels in old array
       !if(any(pix == 46051400)) write(*,*) 'old leftover'
       self%ind2pix(k:k+self%nobs-i) = oldpix(i:self%nobs)
       self%nobs = k+self%nobs-i
    else if (j <= size(newpix)) then
       ! Append left-over pixels in new array
       !if(any(pix == 46051400)) write(*,*) 'new leftover'
       do while (j <= size(newpix))
          self%ind2pix(k) = newpix(j)
          j = j+1; k = k+1
          ! Skip duplicates
          if (j > size(newpix)) exit
          do while (newpix(j) == newpix(j-1))
             j = j+1
             if (j > size(newpix)) exit
          end do
       end do
       self%nobs = k-1
    else
       ! No leftovers; just update nobs
       self%nobs = k-1
    end if
!!$    call wall_time(t2)
!!$    write(*,*) "pix", t2-t1, self%nobs, size(newpix)
!!$    if(any(pix == 46051400)) then
!!$       write(*,*) 'merge', any(pix == 46051400), any(self%ind2pix(1:self%nobs) == 46051400)
!!$       write(*,*) 'merge2', self%nobs
!!$    end if
    
  end subroutine add_pixels

  module subroutine get_ind_range(self, pix1, pix2, ind1, ind2)
    implicit none
    class(comm_tod_pixcache), intent(in)    :: self
    integer(i4b),             intent(in)    :: pix1, pix2
    integer(i4b),             intent(out)   :: ind1, ind2

    if (self%nobs == 0) then
       ind1 = -1
       ind2 = -1
    else
       if (self%ind2pix(self%nobs) < pix1 .or. self%ind2pix(1) > pix2) then
          ind1 = -2
          ind2 = -2
       else
          ind1 = self%pix2ind(pix1, flag_missing=.false.)
          ind2 = self%pix2ind(pix2, flag_missing=.false.)
       end if
    end if
  end subroutine get_ind_range
  
  module subroutine expand_storage(self, n, trim_unused)
    implicit none
    class(comm_tod_pixcache), intent(inout)          :: self
    integer(i4b),             intent(in),   optional :: n
    logical(lgt),             intent(in),   optional :: trim_unused

    integer(i4b) :: m
    integer(i4b), allocatable, dimension(:) :: buffer
    allocate(buffer(self%nobs))
    buffer    = self%ind2pix(1:self%nobs)
    if (present(n)) then
       self%nmax = self%nobs + n
    else if (present(trim_unused)) then
       self%nmax = self%nobs
    else
       self%nmax = 2*self%nobs
    end if
    deallocate(self%ind2pix)
    allocate(self%ind2pix(self%nmax))
    self%ind2pix(1:self%nobs) = buffer
    deallocate(buffer)
  end subroutine expand_storage

  module subroutine precomp_aux(self, npsi)
    implicit none
    class(comm_tod_pixcache), intent(inout) :: self
    integer(i4b),             intent(in)    :: npsi

    integer(i4b) :: i, pix
    real(sp)     :: psi
    real(dp)     :: theta, phi, rotation_matrix(3,3)
    
    allocate(self%ind2pix_nest(self%nobs))
    allocate(self%ind2ang(2,self%nobs))
    allocate(self%ind2vec(3,self%nobs))
    allocate(self%ind2sl(self%nobs))
    allocate(self%ind2vec_ecl(3, self%nobs))
    
    ! Precompute pixel-based lookup tables
    call ecl_to_gal_rot_mat(rotation_matrix)
    do i = 1, self%nobs
       pix = self%ind2pix(i)
       call ring2nest(self%nside, pix, self%ind2pix_nest(i))
       call pix2ang_ring(self%nside, pix, theta, phi)
       self%ind2ang(:,i) = [theta, phi]
       call pix2vec_ring(self%nside, pix, self%ind2vec(:,i))
       call ang2pix_ring(self%nside_sl, theta, phi, self%ind2sl(i))
       self%ind2vec_ecl(:,i) = matmul(self%ind2vec(:,i), rotation_matrix)
    end do
    
    ! Precompute trigonometric functions
    self%npsi = npsi
    allocate(self%sin2psi(self%npsi), self%cos2psi(self%npsi))
    allocate(self%psi(self%npsi))
    do i = 1, self%npsi
       psi             = (i-0.5)*2.0*pi/real(self%npsi,sp)
       self%psi(i)     = psi
       self%sin2psi(i) = sin(2.0*psi)
       self%cos2psi(i) = cos(2.0*psi)
    end do

  end subroutine precomp_aux

  module subroutine init_map_mask(self, map_sky, bitmask, map_gain, scale)
    implicit none
    class(comm_tod_pixcache),                 intent(inout) :: self
    type(map_ptr),       dimension(1:,1:),    intent(in)    :: map_sky ! (ndet,ndelta)
    class(comm_map),     pointer,             intent(in)    :: bitmask 
    type(map_ptr),       dimension(1:),       intent(in), optional :: map_gain
    real(sp),                                 intent(in), optional :: scale

    integer(i4b) :: i, j, k, l, ndet, ndelta
    real(sp), allocatable, dimension(:,:) :: buffer

    ndet   = size(map_sky,1)
    ndelta = size(map_sky,2)

    ! Allocate storage in first call
    if (.not. allocated(self%map_sky)) then
       allocate(self%map_sky(self%nmaps,self%nobs,0:ndet,ndelta))
       allocate(self%bitmask(self%nobs))
       if (present(map_gain)) then
          allocate(self%map_gain(self%nmaps,self%nobs,0:ndet))
       end if
    end if

    ! Distribute sky and (optionally) gain maps
    do j = 1, ndelta
       do i = 1, ndet
          call map_sky(i,j)%p%map2pix(self%ind2pix, self%map_sky(:,:,i,j))
          if (present(map_gain) .and. j == 1) then
             call map_gain(i)%p%map2pix(self%ind2pix, self%map_gain(:,:,i))
          end if
!!$          if (self%nobs> 0) write(*,*) 'sly', j, i, sum(abs(self%map_sky(:,:,i,j)))
!!$          if (self%nobs>0) write(*,*) 'a1', self%map_sky(:,109952,1,1)
       end do

       ! Compute detector average
       do k = 1, self%nobs
          do l = 1, self%nmaps
             self%map_sky(l,k,0,j) = sum(self%map_sky(l,k,1:ndet,j))/ndet
             if (present(map_gain) .and. j == 1) then
                self%map_gain(l,k,0) = sum(self%map_gain(l,k,1:ndet))/ndet
             end if
          end do
       end do
!!$       write(*,*) 'nobs', j, self%nobs, ndet, ndelta


    end do
    
    if (present(scale)) then
       self%map_sky = scale * self%map_sky
       if (present(map_gain)) self%map_gain = scale * self%map_gain
    end if
    
    ! Distribute bitmask
    allocate(buffer(bitmask%info%nmaps,self%nobs))
    call bitmask%map2pix(self%ind2pix, buffer)
    self%bitmask = 0
    do i = 0, 6
       self%bitmask = self%bitmask + nint(1.0-buffer(1,:)) * 2**i
    end do
    deallocate(buffer)
    
  end subroutine init_map_mask
  
end submodule comm_tod_pixcache_smod
