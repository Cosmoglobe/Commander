module rem_tools
  use healpix_types
  use pix_tools
  use math_tools
  use mask_tools
!  use corr_fileutils

  ! *********************************************************************
  ! *                                                                   *
  ! *          Module for removing low-l components from a map          *
  ! *                                                                   *
  ! *               Written by Hans Kristian Eriksen, 2003              *
  ! *                                                                   *
  ! *                Copyright 2003, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************


  integer(i4b),                              private :: map_type, lmax, ordering
  integer(i4b),                              private :: nside, npix, npix_map
  complex(dpc), allocatable, dimension(:,:), private :: coupling_mat
  complex(dpc), allocatable, dimension(:),   private :: exp_iphi
  logical(lgt), allocatable, dimension(:),   private :: mask
  integer(i4b), allocatable, dimension(:),   private :: mask2map, ringnum, numpix_in_ring
  integer(i4b), allocatable, dimension(:,:), private :: ring_pixels
  real(sp),     allocatable, dimension(:),   private :: weights
  real(dp),     allocatable, dimension(:),   private :: phis, thetas
  real(dp),     allocatable, dimension(:,:), private :: P_lms

  interface rem_multi
     module procedure rem_multi_maps
  end interface
  
contains

  subroutine clean_rem_multi
    implicit none

    if (allocated(coupling_mat))   deallocate(coupling_mat)
    if (allocated(exp_iphi))       deallocate(exp_iphi)
    if (allocated(mask))           deallocate(mask)
    if (allocated(mask2map))       deallocate(mask2map)
    if (allocated(ringnum))        deallocate(ringnum)
    if (allocated(numpix_in_ring)) deallocate(numpix_in_ring)
    if (allocated(ring_pixels))    deallocate(ring_pixels)
    if (allocated(weights))        deallocate(weights)
    if (allocated(phis))           deallocate(phis)
    if (allocated(thetas))         deallocate(thetas)
    if (allocated(P_lms))          deallocate(P_lms)

  end subroutine clean_rem_multi

  subroutine init_rem_multi(maskfile, map_type_in, lmax_in)
    implicit none

    character(len=*), intent(in) :: maskfile
    integer(i4b),     intent(in) :: map_type_in, lmax_in

    integer(i4b)   :: numcomp, totnumcomp, numpix, temp_i, unit
    integer(i4b)   :: i, j, l, m, ind, ind1, ind2, l1, m1, l2, m2, iz
    integer(i4b)   :: counter
    complex(dpc)   :: im_i, temp_cmplx
    real(dp)       :: theta, phi, dt, t1, t2
    real(dp),     allocatable, dimension(:) :: plm
    integer(i4b), allocatable, dimension(:) :: map2mask, pixlist

    complex(dpc), allocatable, dimension(:,:) :: sum_exp_iphi

    unit = 59

    im_i = cmplx(0., 1.)
    lmax = lmax_in
    map_type = map_type_in

    ! Read the mask
    call get_maskparam(unit, maskfile, nside, ordering, npix)
    allocate(mask(0:12*nside**2-1))
    call read_mask(unit, maskfile, mask)
    
    if (map_type == 1) then
       npix_map = npix
    else
       npix_map = 12*nside**2
    end if
    
    allocate(weights(0:npix_map-1))
    weights = 1.
    
    ! Create lookup-tables for conversion between map based pixel numbering
    ! and mask based numbering
    allocate(mask2map(0:npix-1))
    allocate(map2mask(0:12*nside**2-1))
    j = 0
    do i = 0, 12*nside**2-1
       if (mask(i)) then
          mask2map(j) = i
          map2mask(i) = j
          j = j+1
       end if
    end do
    
    ! Precalculate the P_lms and phis
    allocate(ringnum(0:npix-1))
    ringnum = -1
    
    allocate(pixlist(0:4*nside-1))
    allocate(numpix_in_ring(4*nside-1))
    allocate(ring_pixels(0:4*nside-1, 4*nside-1))

    do i = 1, 4*nside-1
       
       call in_ring(nside, i, 0.d0, pi, pixlist, numpix)
       call in_ring(nside, i, 0.d0, pi, ring_pixels(:,i), numpix_in_ring(i))
       
       if (ordering == 2) then
          do j = 0, numpix-1
             call ring2nest(nside, pixlist(j), temp_i)
             pixlist(j) = temp_i
             ring_pixels(j,i) = temp_i
          end do
       end if
       
       do j = 0, numpix-1
          if (mask(pixlist(j))) then
             ringnum(map2mask(pixlist(j))) = i
          end if
       end do
       
       j = 0
       do while (j < numpix_in_ring(i))
          if (mask(ring_pixels(j,i))) then
             ring_pixels(j,i) = map2mask(ring_pixels(j,i))
             j = j+1
          else
             ring_pixels(j,i) = ring_pixels(numpix_in_ring(i)-1,i)
             numpix_in_ring(i) = numpix_in_ring(i) - 1
          end if
       end do
       
    end do

    allocate(phis(0:npix-1))
    allocate(exp_iphi(0:npix-1))
    do i = 0, npix-1
       if (ordering == 1) then
          call pix2ang_ring(nside, mask2map(i), theta, phis(i))
       else
          call pix2ang_nest(nside, mask2map(i), theta, phis(i))
       end if
       
       exp_iphi(i) = exp(im_i * phis(i))
    end do
    
    numcomp = (lmax+1)**2
    
    allocate(P_lms(4*nside-1, numcomp))
    allocate(plm(0:lmax))
    
    dt = pi / real(4*nside-1,dp)
    do i = 1, 4*nside-1
       
       call in_ring(nside, i, 0.d0, pi, pixlist, numpix)
       if (ordering == 1) then
          call pix2ang_ring(nside, pixlist(0), theta, phi)
       else
          call ring2nest(nside, pixlist(0), temp_i)
          call pix2ang_nest(nside, temp_i, theta, phi)
       end if
       
       do m = 0, lmax
          
          call comp_normalised_Plm(lmax, m, theta, plm)
          
          do l = m, lmax
             
             ! First the positive m
             ind = l**2 + l + m + 1 
             P_lms(i, ind) = plm(l)
             
             if (m > 0) then
                
                ! Then the negative m
                ind = l**2 + l - m + 1 
                
                if (mod(m,2) == 0) then
                   P_lms(i, ind) = plm(l)
                else
                   P_lms(i, ind) = -plm(l)
                end if
                
             end if
             
          end do
          
       end do
    end do
       
    deallocate(pixlist)
       
    ! Precompute sum over exp_iphi
    allocate(sum_exp_iphi(-2*lmax:2*lmax, 4*nside-1))
    sum_exp_iphi = cmplx(0., 0.)
    do i = 1, 4*nside-1
       do m = -2*lmax, 2*lmax
          do j = 0, numpix_in_ring(i)-1          
             sum_exp_iphi(m,i) = sum_exp_iphi(m,i) + exp_iphi(ring_pixels(j,i))**m
          end do
       end do
    end do
    
    ! Compute the geometric coupling matrix
    allocate(coupling_mat(numcomp, numcomp))
    coupling_mat = cmplx(0.d0, 0.d0)
    
    do l1 = 0, lmax
       do m1 = -l1, l1
          
          ind1 = l1**2 + l1 + m1 + 1 
          
          do l2 = 0, lmax
             do m2 = -l2, l2
                
                ind2 = l2**2 + l2 + m2 + 1 
                
                if (ind1 >= ind2) then
                   
                   do i = 1, 4*nside-1
                      coupling_mat(ind1,ind2) = coupling_mat(ind1,ind2) + &
                           & sum_exp_iphi(m1-m2,i) * P_lms(i,ind1) * P_lms(i,ind2)
                   end do
                   
                   if (ind1 > ind2) then
                      coupling_mat(ind2,ind1) = conjg(coupling_mat(ind1,ind2))
                   end if
                   
                end if
                
             end do
          end do
          
       end do
    end do
    
    ! Multiply with the pixel area
    coupling_mat = coupling_mat * 4.d0*pi / cmplx(real(12*nside**2,dp), 0.d0)

    ! Invert the coupling matrix
!    call invert_matrix2(coupling_mat)

    call invert_matrix(coupling_mat, numcomp)

    deallocate(sum_exp_iphi)
    
  end subroutine init_rem_multi
  


  ! NB: Make sure the input map has the same ordering as the mask before
  !     using this routine
  subroutine rem_multi_maps(map, outmap)
    Implicit none

    real(sp), dimension(1:,0:), intent(in)  :: map
    real(sp), dimension(1:,0:), intent(out) :: outmap

    integer(i4b) :: i, j, k, l, m, nside_mask, ind, temp_i
    integer(i4b) :: numcomp, ordering, nmaps, ind_pos, ind_neg
    real(dp)     :: theta, phi, t1, t2
    complex(dpc) :: im_i

    complex(dpc), allocatable, dimension(:)   :: multipoles, b_lms
    complex(dpc), allocatable, dimension(:)   :: rem_map, subsum1
    complex(dpc), allocatable, dimension(:,:) :: sum_exp_iphi, subsum2

    numcomp  = (lmax+1)**2
    im_i     = cmplx(0., 1.)
    nmaps    = size(map(:,1))
    npix_map = size(map(1,:))

    allocate(rem_map(0:npix_map-1))

    ! Precompute sum over T * exp_i*m*phi
    allocate(sum_exp_iphi(0:lmax, 4*nside-1))
    sum_exp_iphi = cmplx(0., 0.)
    if (map_type == 1) then
       do i = 1, 4*nside-1
          do m = 0, lmax
             do j = 0, numpix_in_ring(i)-1
                sum_exp_iphi(m,i) = sum_exp_iphi(m,i) + &
                     & map(1,ring_pixels(j,i)) * exp_iphi(ring_pixels(j,i))**m
             end do
          end do
       end do
    else
       do i = 1, 4*nside-1
          do m = 0, lmax
             do j = 0, numpix_in_ring(i)-1
                sum_exp_iphi(m,i) = sum_exp_iphi(m,i) + &
                     & map(1,mask2map(ring_pixels(j,i))) * &
                     & exp_iphi(ring_pixels(j,i))**m
             end do
          end do
       end do
    end if

    ! Do the integrations
    allocate(multipoles(numcomp))
    multipoles = cmplx(0.d0,0.d0)

    do l = 0, lmax
       do m = 0, l
          ! First the positive m
          ind_pos = l**2 + l + m + 1 

          do i = 1, 4*nside-1
             ! conjg instead of -m in sum above
             multipoles(ind_pos) = multipoles(ind_pos) + &
                  & conjg(sum_exp_iphi(m,i)) * P_lms(i,ind_pos)
          end do

          if (m > 0) then
             ! Then the negative m
             ind_neg = l**2 + l - m + 1 
 
             if (mod(m,2) == 0) then
                multipoles(ind_neg) = conjg(multipoles(ind_pos))
             else
                multipoles(ind_neg) = -conjg(multipoles(ind_pos))
             end if
          end if

       end do
    end do

    ! Multiply with the pixel area
    multipoles = multipoles * 4.d0*pi / real(12*nside**2,dp)

    allocate(b_lms(numcomp))
    b_lms = cmplx(0.,0.)
    do i = 1, numcomp
       do k = 1, numcomp
          b_lms(k) = b_lms(k) + conjg(coupling_mat(k,i)) * multipoles(i)
       end do
    end do

    write(*,*) b_lms
!    stop

    ! Precalculate sum(b_lm P_lm)
    allocate(subsum1(4*nside-1))
    subsum1 = cmplx(0., 0.)
    do l = 1, lmax
       ind = l**2 + l + 1 
       do i = 1, 4*nside-1
          subsum1(i) = subsum1(i) + b_lms(ind) * P_lms(i, ind)
       end do
    end do

    allocate(subsum2(0:lmax, 4*nside-1))
    subsum2 = cmplx(0.,0.)
    do m = 1, lmax
       do i = 1, 4*nside-1
          do l = m, lmax
             ind = l**2 + l + m + 1 
             subsum2(m,i) = subsum2(m,i) + b_lms(ind) * P_lms(i,ind)
          end do
       end do
    end do


    ! Create the low-order map
    rem_map = cmplx(0.d0, 0.d0)
    if (map_type == 1) then

       do i = 0, npix-1
          
          do m = 1, lmax
             rem_map(i) = rem_map(i) + &
                  & exp_iphi(i)**m * subsum2(m,ringnum(i))
          end do

          rem_map(i) = rem_map(i) + conjg(rem_map(i))

          rem_map(i) = rem_map(i) + b_lms(1) * P_lms(ringnum(i),1) + subsum1(ringnum(i))

       end do

    else

       do i = 0, npix-1
          
          temp_i = mask2map(i)

          do m = 1, lmax
             rem_map(temp_i) = rem_map(temp_i) + &
                  & exp_iphi(i)**m * subsum2(m,ringnum(i))
          end do

          rem_map(temp_i) = rem_map(temp_i) + conjg(rem_map(temp_i))

          rem_map(temp_i) = rem_map(temp_i) + &
               & b_lms(1) * P_lms(ringnum(i),1) + subsum1(ringnum(i))

       end do

    end if

    if (nmaps > 1) then
       outmap(2:nmaps,:) = map(2:nmaps,:)
    end if

    if (map_type == 1) then
       outmap(1,:) = outmap(1,:) - real(rem_map,sp)
    else
       do i = 0, 12*nside**2-1
          if (mask(i)) then
             outmap(1,i) = map(1,i) - real(rem_map(i),sp)
          else
             outmap(1,i) = -1.6375e30
          end if
       end do
    end if

    deallocate(subsum1)
    deallocate(subsum2)
    deallocate(b_lms)
    deallocate(sum_exp_iphi)
    deallocate(rem_map)
    
  end subroutine rem_multi_maps



end module rem_tools
