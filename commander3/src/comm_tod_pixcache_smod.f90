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

  module function constructor_tod_pixcache(nside, nside_sl, nmaps, fullsky, equal_dets, nside_zodi) result(c)
    implicit none
    integer(i4b),             intent(in)           :: nside, nside_sl, nmaps
    integer(i4b),             intent(in), optional :: nside_zodi
    logical(lgt),             intent(in)           :: fullsky, equal_dets
    class(comm_tod_pixcache), pointer              :: c

    allocate(c)
    c%nside        = nside
    c%nside_sl     = nside_sl
    c%nmaps        = nmaps
    c%fullsky      = fullsky
    c%nmax         = 1
    c%nobs         = 0
    c%equal_dets   = equal_dets
    allocate(c%ind2pix(c%nmax))

    if (present(nside_zodi)) then
      c%nside_zodi   = nside_zodi
      c%nobs_zodi    = 0
      c%nmax_zodi    = 1
      allocate(c%ind2pix_zodi(c%nmax_zodi))
    end if
    
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
          if (ind1 == -1) ind1 = 1
          if (self%ind2pix(ind1) < pix1) ind1 = ind1+1
          if (ind2 == -1) ind2 = self%nobs
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

    if (self%equal_dets) then
       ndet = 0
    else
       ndet   = size(map_sky,1)
    end if
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
       if (self%equal_dets) then
          ! Fill detector average only
          call map_sky(1,j)%p%map2pix(self%ind2pix, self%map_sky(:,:,0,j))
          if (present(map_gain) .and. j == 1) then
             call map_gain(1)%p%map2pix(self%ind2pix, self%map_gain(:,:,0))
          end if
       else
          do i = 1, ndet
             call map_sky(i,j)%p%map2pix(self%ind2pix, self%map_sky(:,:,i,j))
             if (present(map_gain) .and. j == 1) then
                call map_gain(i)%p%map2pix(self%ind2pix, self%map_gain(:,:,i))
             end if
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
       end if
    end do
    
    if (present(scale)) then
       self%map_sky = scale * self%map_sky
       if (present(map_gain)) self%map_gain = scale * self%map_gain
    end if
    
    ! Distribute bitmask
    allocate(buffer(bitmask%info%nmaps,self%nobs))
    call bitmask%map2pix(self%ind2pix, buffer)
    self%bitmask = buffer(1,:)
    deallocate(buffer)
    
!    allocate(buffer(bitmask%info%nmaps,self%nobs))
!    call bitmask%map2pix(self%ind2pix, buffer)
!    self%bitmask = 0
!    do i = 0, 6
!       self%bitmask = self%bitmask + nint(1.0-buffer(1,:)) * 2**i
!    end do
!    deallocate(buffer)
    
  end subroutine init_map_mask

  module function pix2ind_zodi(self, pix_solar, flag_missing) result(ind)
    implicit none
    class(comm_tod_pixcache), intent(in)           :: self
    integer(i4b),             intent(in)           :: pix_solar
    logical(lgt),             intent(in), optional :: flag_missing
    integer(i4b)                                   :: ind
    logical(lgt) :: flag_missing_
    if (self%fullsky) then
       ind = pix_solar+1
    else
       if (self%nobs_zodi == 0) then
          ind = 0
          flag_missing_ = .true.; if (present(flag_missing)) flag_missing_ = flag_missing
          if (flag_missing_) ind = -1
       else
          ind = locate(self%ind2pix_zodi(1:self%nobs_zodi), pix_solar)
          flag_missing_ = .true.; if (present(flag_missing)) flag_missing_ = flag_missing
          if (flag_missing_) then
             if (.not. self%ind2pix_zodi(ind) == pix_solar) then
                write(*,*) "pix2ind_zodi", pix_solar, ind, self%nobs_zodi, self%ind2pix_zodi(ind-1:ind+1)
                ind = -1
             end if
          end if
       end if
    end if
  end function pix2ind_zodi

  module subroutine add_pixels_zodi(self, pix_solar)  !uses solar_centric pointings
    implicit none
    class(comm_tod_pixcache),               intent(inout) :: self
    integer(i4b),             dimension(:), intent(in)    :: pix_solar

    integer(i4b) :: i, j, k, p, n, a, pix_nest, downg_divide 
    integer(i4b), allocatable, dimension(:) :: newpix, oldpix, pix_downg

    if (self%fullsky) return

    !pix_solar are pointings at the native nside of the hdf files, 
    !I want to store them with nside=nside_zodi so I have to downgrade them
    allocate(pix_downg(size(pix_solar)))
    if (self%nside /= self%nside_zodi) then !Angela:check that self%nside is the nside of pix_solar
      if (self%nside > self%nside_zodi) then 
        downg_divide = (self%nside/self%nside_zodi)**2
        do a = 1, size(pix_solar)
           call ring2nest(self%nside, pix_solar(a), pix_nest)
           pix_nest = pix_nest/downg_divide 
           call nest2ring(self%pixcache%nside_zodi, pix_nest, pix_downg(a)) 
           !Angela: check whether this is ok or better to use vector space
        end do
      else 
        write(*,*) "Error in pixcache%add_pixels_zodi: nside_zodi larger than nside."
      end if
    end if 

    allocate(newpix(size(pix_downg)), oldpix(self%nobs_zodi))
    newpix = pix_downg
    oldpix = self%ind2pix_zodi(1:self%nobs_zodi)
    call QuickSort_int(newpix)

    ! Count unique pixels
    n = 1
    do p = 2, size(newpix)
       if (newpix(p) /= newpix(p-1)) n = n+1
    end do
    call self%expand_storage_zodi(n)
    
    ! Merge new sorted list into existing
    i = 1; j = 1; k = 1
    do while (i <= self%nobs_zodi .and. j <= size(newpix))
       if (oldpix(i) == newpix(j)) then
          self%ind2pix_zodi(k) = oldpix(i)
          i = i+1; j = j+1; k = k+1
       else if (oldpix(i) < newpix(j)) then
          self%ind2pix_zodi(k) = oldpix(i)
          i = i+1; k = k+1
       else if (oldpix(i) > newpix(j)) then
          self%ind2pix_zodi(k) = newpix(j)
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
    
    if (i <= self%nobs_zodi) then
       ! Append left-over pixels in old array
       self%ind2pix_zodi(k:k+self%nobs_zodi-i) = oldpix(i:self%nobs_zodi)
       self%nobs_zodi = k+self%nobs_zodi-i
    else if (j <= size(newpix)) then
       ! Append left-over pixels in new array
       do while (j <= size(newpix))
          self%ind2pix_zodi(k) = newpix(j)
          j = j+1; k = k+1
          ! Skip duplicates
          if (j > size(newpix)) exit
          do while (newpix(j) == newpix(j-1))
             j = j+1
             if (j > size(newpix)) exit
          end do
       end do
       self%nobs_zodi = k-1
    else
       ! No leftovers; just update nobs_zodi
       self%nobs_zodi = k-1
    end if    

    deallocate(newpix, oldpix, pix_downg)

  end subroutine add_pixels_zodi

  module subroutine expand_storage_zodi(self, n, trim_unused)
    implicit none
    class(comm_tod_pixcache), intent(inout)          :: self
    integer(i4b),             intent(in),   optional :: n
    logical(lgt),             intent(in),   optional :: trim_unused

    integer(i4b), allocatable, dimension(:) :: buffer

    allocate(buffer(self%nobs_zodi))
    buffer    = self%ind2pix_zodi(1:self%nobs_zodi)
    if (present(n)) then
       self%nmax_zodi = self%nobs_zodi + n
    else if (present(trim_unused)) then
       self%nmax_zodi = self%nobs_zodi
    else
       self%nmax_zodi = 2*self%nobs_zodi
    end if
    deallocate(self%ind2pix_zodi)
    allocate(self%ind2pix_zodi(self%nmax_zodi))
    self%ind2pix_zodi(1:self%nobs_zodi) = buffer
    deallocate(buffer)
  end subroutine expand_storage_zodi

  module subroutine init_map_zodi(self, start_time, end_time, timeslice, n_detectors)
    implicit none
    class(comm_tod_pixcache),         intent(inout) :: self
    integer(i4b)                                    :: n_detectors
    real(sp),                         intent(in)    :: start_time, end_time, timeslice 
    !begin of the mission, end of the mission, time between one zodi lookup maps and the other, all in MJD

    integer(i4b)                                    :: n_slices

    ! Allocate storage in first call
    if (.not. allocated(self%map_zodi)) then
      n_slices = ((end_time - start_time)/timeslice) + 2 !start and end time + all slices in the middle, the last slice might be smaller !Angela: check whether this rounding is okay
      allocate(self%map_zodi(self%nobs_zodi,n_slices,n_detectors))
    end if

    ! Distribute sky and (optionally) gain maps !Angela: this is done in init_map_mask, relevant also here?

  end subroutine init_map_zodi
  
end submodule comm_tod_pixcache_smod
