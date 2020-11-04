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
module comm_chisq_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_ptsrc_comp_mod
  use comm_template_comp_mod
  implicit none


contains

  subroutine compute_chisq(comm, chisq_map, chisq_fullsky, mask, lowres_eval, evalpol)
    implicit none
    integer(i4b),                   intent(in)              :: comm
    logical(lgt),                   intent(in),    optional :: lowres_eval
    logical(lgt),                   intent(in),    optional :: evalpol
    !character(len=512),             intent(in),    optional :: evalsig
    !logical(lgt),                   intent(in),    optional :: udgrade_chisq
    class(comm_map),                intent(inout), optional :: chisq_map
    real(dp),                       intent(out),   optional :: chisq_fullsky
    type(map_ptr),   dimension(1:), intent(in),    optional :: mask

    integer(i4b) :: i, j, k, p, ierr, nmaps
    real(dp)     :: t1, t2
    logical(lgt) :: apply_mask, lowres
    class(comm_map), pointer :: res, res_lowres => null(), res_lowres_temp, chisq_sub
    class(comm_mapinfo), pointer :: info, info_lowres

    if (present(chisq_fullsky) .or. present(chisq_map)) then
       if (present(chisq_fullsky)) chisq_fullsky = 0.d0
       if (present(chisq_map))     chisq_map%map = 0.d0
       do i = 1, numband
          
          ! Skip non-essential chisq evaluation
          if (present(evalpol)) then
             if (evalpol) then
                if (.not. data(i)%info%pol) cycle
             else
                if (data(i)%info%pol) cycle
             end if
          end if

          res => compute_residual(i)

          apply_mask = present(mask)
          if (apply_mask) apply_mask = associated(mask(i)%p)
          if (apply_mask) then
             res%map = res%map * mask(i)%p%map
             !call res%writeFITS("chisq.fits")
             !call mask(i)%p%writeFITS("mask.fits")
             !call mpi_finalize(j)
             !stop
          end if
          
          if (trim(data(i)%N%type) == "rms" .and. data(i)%N%nside_chisq_lowres < res%info%nside .and. present(chisq_fullsky) .and. present(lowres_eval)) then
             if (lowres_eval) then
                lowres = .true.
                info_lowres  => comm_mapinfo(data(i)%info%comm, data(i)%N%nside_chisq_lowres, 0, data(i)%info%nmaps, data(i)%info%nmaps==3)

                res_lowres => comm_map(info_lowres)
                res_lowres_temp => comm_map(info_lowres)

                call res%udgrade(res_lowres)
                res_lowres_temp%map = res_lowres%map ! Save temporarily

                call data(i)%N%invN_lowres(res_lowres) ! invN*res
                res_lowres%map = res_lowres_temp%map*res_lowres%map ! res*(invN*res)

                call res_lowres_temp%dealloc(); deallocate(res_lowres_temp)
             end if
          else
             lowres=.false.
             call data(i)%N%sqrtInvN(res) 
             res%map = res%map**2 !(sqrtInvN*res)**2 = res*invN*res
          end if

          
          if (present(chisq_map)) then
             info  => comm_mapinfo(data(i)%info%comm, chisq_map%info%nside, 0, data(i)%info%nmaps, data(i)%info%nmaps==3)
             chisq_sub => comm_map(info)
             call res%udgrade(chisq_sub)
             do j = 1, data(i)%info%nmaps
                chisq_map%map(:,j) = chisq_map%map(:,j) + chisq_sub%map(:,j) * (res%info%npix/chisq_sub%info%npix)
             end do
             call chisq_sub%dealloc(); deallocate(chisq_sub)
          end if
          if (present(chisq_fullsky)) then
             if (lowres) then
                chisq_fullsky = chisq_fullsky + sum(res_lowres%map)
             else
                chisq_fullsky = chisq_fullsky + sum(res%map)
             end if
          end if

          if (associated(res_lowres)) then
             call res_lowres%dealloc(); deallocate(res_lowres)
             nullify(res_lowres)
          end if
          call res%dealloc(); deallocate(res)
       end do
    end if

    if (present(chisq_fullsky)) then
       call mpi_allreduce(MPI_IN_PLACE, chisq_fullsky, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    end if

  end subroutine compute_chisq


  subroutine compute_jeffreys_prior(c, df, pol, par, chisq_jeffreys)
    implicit none
    class(comm_diffuse_comp),       intent(in)              :: c
    type(map_ptr),   dimension(1:), intent(in)              :: df
    integer(i4b),                   intent(in)              :: pol, par
    real(dp),                       intent(out)             :: chisq_jeffreys

    integer(i4b) :: i, j, k, p, ierr, nmaps
    class(comm_map),     pointer :: map
    class(comm_mapinfo), pointer :: info 


    chisq_jeffreys = 0.d0
    do i = 1, numband   

       if (c%F_null(i,0)) cycle

       ! Initialize amplitude map
       nmaps =  min(data(i)%info%nmaps, c%nmaps)
       info  => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, data(i)%info%lmax, nmaps, nmaps==3)
       map   => comm_map(info)
       call c%x%alm_equal(map)
       call map%Y()

       ! Multiply with derivative of mixing matrix
       map%map = map%map * df(i)%p%map

       ! Nullify inactive polarization components, depending on polarization type
       if (c%poltype(par) == 2) then
          if (pol == 1) then
             if (map%info%nmaps > 1) map%map(:,2:3) = 0.d0
          else
             map%map(:,1) = 0.d0
          end if
       else if (c%poltype(par) == 3) then
          do k = 1, map%info%nmaps
             if (k /= pol) map%map(:,k) = 0.d0
          end do
       end if

       ! Convolve with band-specific beam
       call map%YtW()
       call data(i)%B(0)%p%conv(trans=.false., map=map)
       call map%Y()

       ! Apply mask if requested
       if (allocated(c%indmask)) map%map = map%map * c%indmask(i)%p%map
          
       ! Compute Jeffreys prior, df*invN*df
       call data(i)%N%sqrtInvN(map)
       chisq_jeffreys = chisq_jeffreys + sum(map%map**2)

       call map%dealloc(); deallocate(map)
    end do

    call mpi_allreduce(MPI_IN_PLACE, chisq_jeffreys, 1, MPI_DOUBLE_PRECISION, MPI_SUM, c%comm, ierr)    
    if (chisq_jeffreys > 0.d0) then
       chisq_jeffreys = -2.d0*log(sqrt(chisq_jeffreys))
    else
       chisq_jeffreys = 0.d0
    end if

  end subroutine compute_jeffreys_prior



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
          if (c%active_samp_group(cg_samp_group)) skip = .true.
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
    call ptsrc%dealloc(); deallocate(ptsrc)

  end function compute_residual

  subroutine subtract_fiducial_CMB_dipole(band, map)
    implicit none
    integer(i4b),    intent(in)    :: band
    class(comm_map), intent(inout) :: map

    integer(i4b)        :: i, j, l, m
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer     :: dipole
    real(dp),            allocatable, dimension(:,:) :: alm
    class(comm_comp),    pointer :: c

    ! Compute predicted signal for this band
    c => compList
    do while (associated(c))
       if (trim(c%type) /= 'cmb') then
          c => c%next()
          cycle
       end if
       
       select type (c)
       class is (comm_diffuse_comp)
          dipole => comm_map(data(band)%info)
          allocate(alm(0:data(band)%info%nalm-1,data(band)%info%nmaps))
          alm = 0.d0
          do j = 0, data(band)%info%nalm-1
             l = data(band)%info%lm(1,j)
             m = data(band)%info%lm(2,j)
             if (l == 1 .and. m == -1) then
                alm(j,1) = -4.54107d3 / c%RJ2unit_(1)
             else if (l==1 .and. m == 0) then
                alm(j,1) = 5.119744d3 / c%RJ2unit_(1)
             else if (l == 1 .and. m == 1) then
                alm(j,1) = 4.848587d2 / c%RJ2unit_(1)
             end if
          end do
          dipole%alm = c%getBand(band, amp_in=alm, alm_out=.true.)
          call dipole%Y()
          map%map = map%map - dipole%map
          deallocate(alm)
          call dipole%dealloc(); deallocate(dipole)
       end select
       c => c%next()
    end do

    ! Clean up
    nullify(c)

  end subroutine subtract_fiducial_CMB_dipole

  subroutine add_fiducial_CMB_dipole(info, RJ2unit, alm)
    implicit none
    class(comm_mapinfo),                   intent(in)    :: info
    real(dp),                              intent(in)    :: RJ2unit
    real(dp),            dimension(0:,1:), intent(inout) :: alm

    integer(i4b)        :: i, j, l, m

    do j = 0, info%nalm-1
       l = info%lm(1,j)
       m = info%lm(2,j)
       if (l == 1 .and. m == -1) then
          alm(j,1) = alm(j,1) - 4.54107d3 / RJ2unit
       else if (l==1 .and. m == 0) then
          alm(j,1) = alm(j,1) + 5.11974d3 / RJ2unit
       else if (l == 1 .and. m == 1) then
          alm(j,1) = alm(j,1) + 4.84858d2 / RJ2unit
       end if
    end do

  end subroutine add_fiducial_CMB_dipole

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
             !allocate(alm(0:data(i)%info%nalm-1,data(i)%info%nmaps))
             !allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))          
             out%alm = c%getBand(i, alm_out=.true.)
             !call out%add_alm(alm, c%x%info)
             call out%Y()
             !deallocate(alm)
          class is (comm_ptsrc_comp)
             !allocate(map(0:data(i)%info%np-1,data(i)%info%nmaps))
             out%map     = c%getBand(i)
             !out%map = out%map + map
             !deallocate(map)
          class is (comm_template_comp)
              if (c%band /= i) skip = .true.
             if (.not. skip) then
                !allocate(map(0:data(i)%info%np-1,data(i)%info%nmaps))
                out%map     = c%getBand(i)
                !out%map = out%map + map
                !deallocate(map)
             end if
          end select
          filename = trim(outdir)//'/'//trim(c%label)//'_'//trim(data(i)%label)//'_'//trim(postfix)//'.fits'
          !call data(i)%apply_proc_mask(out)
          if (.not. skip) call out%writeFITS(filename)
          c => c%next()
       end do
       call out%dealloc; deallocate(out)
    end do

    ! Clean up
    nullify(c)

  end subroutine output_signals_per_band

  subroutine get_sky_signal(band, det, map_out, mono, cmb_pol)
    implicit none
    integer(i4b),    intent(in)     :: band, det
    class(comm_map), pointer        :: map_out
    logical(lgt), optional          :: mono 
    logical(lgt), optional          :: cmb_pol

    integer(i4b) :: i, j, k
    logical(lgt) :: skip, mono_
    class(comm_map),  pointer :: map_diff
    class(comm_comp), pointer :: c
    real(dp),     allocatable, dimension(:,:) :: map, alm
    real(dp),                  dimension(5)   :: P_quad
    
    mono_ = .true.; if (present(mono)) mono_=mono 

    P_quad = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
    !P_quad = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0]  ! NPOPE
    !P_quad = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0]  ! SROLL2

    ! Allocate map
    map_out  => comm_map(data(band)%info)  
    map_diff => comm_map(data(band)%info)

    ! Compute predicted signal for this band
    c => compList
    map_out%alm  = 0.d0
    map_out%map  = 0.d0
    map_diff%alm = 0.d0
    do while (associated(c))
       if (.not. mono_ .and. trim(c%type)=="md") then
          c => c%next()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
          !allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))          
          alm     = c%getBand(band, alm_out=.true., det=det)
!!$          if (c%x%info%myid == 0) then
!!$             write(*,*) c%label
!!$             write(*,*) shape(alm)
!!$             write(*,*) shape(map_out%alm)
!!$          end if
          !write(*,*) c%x%info%nalm, map_diff%info%nalm, c%x%info%nmaps, map_diff%info%nmaps
!          call map_diff%add_alm(alm, c%x%info)
          if (present(cmb_pol) .and. trim(c%label) == 'cmb') then
             !map_diff%alm(:,1) = map_diff%alm(:,1) + alm(:,1)
             do j = 1, data(band)%info%nmaps
                do i = 0, data(band)%info%nalm-1
                   if (j == 1 .or. data(band)%info%lm(1,i) > 2) then
                      map_diff%alm(i,j) = map_diff%alm(i,j) + alm(i+1,j)
                   else if (j == 2 .and. data(band)%info%lm(1,i) == 2) then
                      ! Apply external E quadrupole prior
                      k = 3+data(band)%info%lm(2,i)
                      map_diff%alm(i,j) = map_diff%alm(i,j) + P_quad(k)
                   end if
                end do
             end do
          else
             map_diff%alm = map_diff%alm + alm
          end if
          deallocate(alm)
       class is (comm_ptsrc_comp)
          allocate(map(0:data(band)%info%np-1,data(band)%info%nmaps))
          map         = c%getBand(band, det=det)
          map_out%map = map_out%map + map
          deallocate(map)
       class is (comm_template_comp)
          allocate(map(0:data(band)%info%np-1,data(band)%info%nmaps))
          map         = c%getBand(band, det=det)
          map_out%map = map_out%map + map
          deallocate(map)
       end select
       c => c%next()
    end do
    
    call map_diff%Y()

    ! Compute residual map
    map_out%map = map_out%map + map_diff%map

    ! Clean up
    nullify(c)
    call map_diff%dealloc; deallocate(map_diff)

  end subroutine get_sky_signal


  subroutine compute_marginal(mixing, data, invN, marg_map, marg_fullsky)
    implicit none
    
    real(c_double),  intent(in),    dimension(:,:,:) :: mixing   !(nbands,ncomp,npix) mixing matrix
    real(c_double),  intent(in),    dimension(:,:)   :: invN     !(nbands,npix) inverse noise matrix
    real(c_double),  intent(in),    dimension(:,:)   :: data     !(nbands,npix) data matrix
    class(comm_map), intent(inout), optional         :: marg_map
    real(dp),        intent(out),   optional         :: marg_fullsky
    integer :: i, j, k, l, p, ierr, nb, npix, nc
    logical :: temp_bool
    double precision     :: temp_marg
    double precision, allocatable, dimension(:)    :: MNd    ! (M.T*invN*M)
    double precision, allocatable, dimension(:)    :: M_d    ! (M.T*invN*M)^-1 * (M.T*invN*d)
    double precision, allocatable, dimension(:,:)  :: MN     ! M.T*ivnN
    double precision, allocatable, dimension(:,:)  :: MNM    ! M.T*ivnN*M (and its inverse)
    double precision, allocatable, dimension(:,:)  :: invmat ! matrix to invert MNM
    double precision, allocatable, dimension(:)    :: temp_arr ! array to flip rows/columns in matrices

    if (present(marg_fullsky) .or. present(marg_map)) then
       if (present(marg_fullsky)) marg_fullsky = 0.d0
       if (present(marg_map))     marg_map%map = 0.d0

       ! pixel last to speed up lookup time (this can be easily changed if needed)
       nb   = size(mixing(:,1,1)) !we assume 1st dimension of mixing matrix to be nbands
       nc   = size(mixing(1,:,1)) !we assume 2nd dimension of mixing matrix to be ncomp
       npix = size(mixing(1,1,:)) !we assume 3rd dimension of mixing matrix to be npix

       ! allocate temporary arrays and matrices 
       allocate(MN(nc,nb),MND(nc),MNM(nc,nc),M_d(nc),invmat(nc,2*nc),temp_arr(2*nc))

       ! for each pixel
       do p = 0,npix-1
          ! calc M.T*invN
          do i = 1,nb
             MN(:,i) = mixing(i,:,p)*invN(i,p)
          end do

          ! calc M.T*invN*d
          do i = 1,nc
             MNd(i) = sum(MN(i,:)*data(:,p))
          end do

          ! calc M.T*invN*M
          do i = 1,nc
             do j = 1,nc
                MNM(i,j) = sum(MN(i,:)*mixing(:,j,p))
             end do
          end do

          ! invert MNM
          if (nc==1) then
             MNM = 1.d0/MNM
          else
             !!! some function to compute the invese of a matrix 
             !!! Need to consider potential zeroes! We are scaling many orders of magnitude
             !!! (need to aviod division by zero among other concerns)
             invmat(:,:)=0.d0
             invmat(:,:nc)=MNM
             do j=1,nc
                invmat(j,j+nc) = 1.d0
             end do

             do j = 1,nc
                if (invmat(j,j)==0.d0) then
                   temp_bool = .true.
                   k = j+1
                   do while ((k <= nc) .and. temp_bool)
                      if (invmat(k,j) /= 0.d0) then !flip with row j
                         temp_arr(:) = invmat(j,:)
                         invmat(j,:) = invmat(k,:)
                         invmat(k,:) = temp_arr(:)
                         temp_bool   = .false.
                      end if
                      k = k + 1
                   end do
                   !if (temp_bool==.true.) then !not possible to invert matrix
                   if (temp_bool) then !not possible to invert matrix
                      temp_marg = -1.d30
                      goto 1
                   end if
                end if

                invmat(j,j+1:) = invmat(j,j+1:)/invmat(j,j) !normalize row with the first non-zero digit of the row (i.e. j)
                invmat(j,j)    = 1.d0 !escape problem with precision

                do k = j+1,nc
                   ! for each row after row j, subtract row j * the digit in column j of that row
                   invmat(k,:)=invmat(k,:)-invmat(k,j)*invmat(j,:)
                   invmat(k,j)=0.d0 !to escape later problems with precision
                end do
             end do
             ! now the matrix should look like this
             ! |1 x x x y y y y|
             ! |0 1 x x y y y y|
             ! |0 0 1 x y y y y| 
             ! |0 0 0 1 y y y y|
             ! 
             ! go the other way back up, from the bottom
             do j = nc,1,-1
                do k = 1,j-1
                   invmat(k,:)=invmat(k,:)-invmat(k,j)*invmat(j,:)
                   invmat(k,j)=0.d0 !to escape later problems with precision                   
                end do
             end do
             !set MNM equal to its inverse
             MNM(:,:)=invmat(:,nc+1:) !the y part of the matrix above

          end if

          ! calc (M.T*invN*M)^-1 (M.T*invN*d)
          do i = 1,nc
             M_d(i) = sum(MNM(i,:)*MNd(:))
          end do

          ! calc final value
          temp_marg = -sum(MNd(:)*M_d(:))

1         if (present(marg_map))     marg_map%map(p,1) = temp_marg
          if (present(marg_fullsky)) marg_fullsky  = marg_fullsky + temp_marg
       end do

       deallocate(MN,MND,MNM,M_d,invmat,temp_arr)
    end if
  end subroutine compute_marginal

end module comm_chisq_mod
