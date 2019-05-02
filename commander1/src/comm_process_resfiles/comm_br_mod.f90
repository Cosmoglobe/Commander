module comm_br_mod
  use comm_like_utils
  implicit none

  ! *********************************************************************
  ! *  comm_br_mod -- An F90 module for computing the Blackwell-Rao     *
  ! *                 estimator given signal samples from the posterior *
  ! *                                                                   *
  ! *           H. K. Eriksen, E. Gjerløw (University of Oslo)          *
  ! *                                                                   *   
  ! *                               and                                 *
  ! *                                                                   *   
  ! *                         I. K. Wehus (JPL)                         *   
  ! *                                                                   *
  ! *                                                                   *
  ! *  Please cite the following papers when using this code:           *
  ! *                                                                   *
  ! *  - H. K. Eriksen et al. 2008, ApJ, 676, 10  (Commander)           *
  ! *  - M. Chu et al. 2005, Phys. Rev. D, 71, 103002 (Blackwell-Rao)   *
  ! *                                                                   *
  ! *  History:                                                         *
  ! *      March 20th, 2013 -- Planck 2013 release                      *
  ! *      July  12th, 2013 -- Generalized to polarization              *
  ! *      Nov   22nd, 2013 -- Support for multiple BR elements with    *
  ! *                          overlap                                  *
  ! *                                                                   *
  ! *                                                                   *
  ! *********************************************************************

  ! ================================================================================================
  ! User routines:
  !
  !  ** subroutine comm_br_initialize_object(initfile, handle)
  !
  !    Input parameters:
  !      initfile   (character=*) :: BR initialization file
  ! 
  !    Output parameters:
  !      handle      (i4b)         :: this code supports multiple BR objects (optional; default=1), and
  !                                   the handle identifies specific BR data sets
  !
  ! 
  !  ** function comm_br_compute_lnL(cls, handle)
  !
  !    Input parameters:
  !      cls         (dp)          :: array containing at least cls(2:lmax) in units of l(l+1)/2pi
  !      handle      (i4b)         :: handle selecting which BR estimator to use (optional; default = 1)
  !    
  !  
  !  ** subroutine comm_br_deallocate_object(handle)
  ! 
  !    Input parameters:
  !      handle      (i4b)         :: handle of object to deallocate (optional; default=remove all)
  !
  ! ================================================================================================

  integer(i4b),       parameter       :: MAX_N_BR   = 5
  integer(i4b),       parameter       :: MAX_N_FACT = 1000
  integer(i4b),       parameter       :: MIN_NUM_ACTIVE_SAMPLES = 1

  type comm_br_factor
     integer(i4b) :: lmin, lmax, p, n, d
     character(len=1), allocatable, dimension(:,:)   :: status
     integer(i4b),     allocatable, dimension(:)     :: ind
     real(dp),         allocatable, dimension(:)     :: log_det_sigma
     real(dp),         allocatable, dimension(:,:)   :: cl
     real(dp),         allocatable, dimension(:,:,:) :: sigma
  end type comm_br_factor

  type comm_P_struct
     integer(i4b) :: nbin, nfactor, nsamp, nspec, nmaps
     real(dp)     :: offset, loglike_weight
     type(comm_br_factor), allocatable, dimension(:) :: factor
  end type comm_P_struct

  type comm_br_data
     logical(lgt) :: initialized=.false.
     integer(i4b) :: n_p
     type(comm_P_struct), allocatable, dimension(:) :: P
  end type comm_br_data

  type(comm_br_data), dimension(MAX_N_BR) :: comm_br
  

contains

  ! Initialization routines
  subroutine comm_br_initialize_object(BRfile, handle, ell_min, ell_max)
    implicit none

    character(len=*),  intent(in)  :: BRfile
    integer(i4b),      intent(out), optional :: handle
    integer(i4b),      intent(in),  optional :: ell_min, ell_max

    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep
    integer(i4b)       :: b, i, j, k, l, m, n, q, s, c1, c2, f, ind, d, n_h
    integer(i4b)       :: unit, numsamples, numchains, lmax_chain, id, bin(2), n_p, n_g, nmode
    integer(i4b)       :: col, p, nmaps, nspec, nsamp, lmax_cl, lmin_bin, lmax_bin
    logical(lgt)       :: pattern(3,3), polarization, exist
    real(dp)           :: threshold, lnL
    character(len=1)   :: status(6)
    character(len=512) :: line, s1, s2
    character(len=128) :: sigmafile, clfile, offsetfile
    real(dp),     allocatable, dimension(:,:)     :: cls, cls_offset
    integer(i4b), allocatable, dimension(:)       :: i2p
    real(dp),     allocatable, dimension(:,:,:,:) :: sigma
    real(dp),     allocatable, dimension(:,:,:)   :: sigma_1D

    id = 1
    unit = comm_getlun()
    if (present(handle)) then
       do while (comm_br(id)%initialized .and. id < MAX_N_BR)
          id = id+1
       end do
       handle = id
    end if
    if (id < 1 .or. id > MAX_N_BR) stop 'Error -- planck_br_mod: BR id out of range'

    ! Find number of unique posterior distributions to include
    inquire(file=trim(BRfile), exist=exist)
    if (.not. exist) then
       write(*,*) 'comm_br_mod -- Error: ', trim(BRfile), ' does not exist'
       stop
    end if
    open(unit,file=trim(BRfile))
    read(unit,fmt='(a)') line
    do while (line(1:1) == '#' .or. trim(line) == '') 
       read(unit,fmt='(a)') line
       cycle ! Skip header comments
    end do
    read(line,*) s1, s2, n_p
    read(unit,*) s1, s2, offsetfile
    call read_fiducial_spectrum(offsetfile, cls_offset)

    ! Do BR posteriors
    if (n_p > 0) then
       comm_br(id)%n_p = n_p
       allocate(comm_br(id)%P(n_p))
       
       ! Initialize all distributions
       comm_br(id)%initialized = .true.
       do i = 1, n_p
          
          ! Allocate factor structures
          allocate(comm_br(id)%P(i)%factor(MAX_N_FACT))
          
          ! Skip header lines
          read(unit,fmt='(a)') line
          do while (line(1:1) == '#' .or. trim(line) == '') 
             read(unit,fmt='(a)') line
             cycle        
          end do
          
          ! Read info from BR file
          read(line,*) sigmafile, comm_br(id)%P(i)%loglike_weight, comm_br(id)%P(i)%nbin, &
               & firstchain, lastchain, firstsample, lastsample, thinstep

          ! Read sigma samples
          inquire(file=trim(sigmafile), exist=exist)
          if (.not. exist) then
             write(*,*) 'Error -- sigma file = ', trim(sigmafile), ' does not exist'
             stop
          end if
          call comm_read_chain(sigmafile, unit, lmax_chain, numchains, numsamples, sigma)
          nspec = size(sigma,2)
          nmaps = 1; if (nspec == 6) nmaps = 3
          comm_br(id)%P(i)%nspec = nspec
          comm_br(id)%P(i)%nmaps = nmaps

          ! Consistency check
          if (numchains < lastchain) then
             write(*,*) 'Error -- comm_br_mod: Requested parameters inconsistent with sigma file'
             write(*,*) ''
             write(*,*) '    sigmafile  = ', trim(sigmafile)
             write(*,*) '    numchains  = ', numchains,  ' vs requested lastchain  = ', lastchain
             stop
          end if

          ! Find number of valid samples
          nsamp = 0
          do j = firstchain, lastchain
             do k = firstsample, min(lastsample,numsamples), thinstep
                if (any(sigma(:,:,j,k) /= 0.d0)) nsamp = nsamp+1
             end do
          end do
          allocate(sigma_1D(0:lmax_chain,nspec,nsamp))
          n = 0
          do j = firstchain, lastchain
             do k = firstsample, min(lastsample,numsamples), thinstep
                if (any(sigma(:,:,j,k) /= 0.d0)) then
                   n     = n+1
                   sigma_1D(:,:,n) = sigma(:,:,j,k)
                end if
             end do
          end do

          ! Use sigma posterior mean for normalization
          allocate(cls(0:lmax_chain,nspec))
          cls(0:lmax_chain,:) = cls_offset(0:lmax_chain,1:nspec)
!!$          do j = 1, nspec
!!$             do l = 2, lmax_chain
!!$                cls(l,j) = sum(sigma_1D(l,j,:)) / nsamp
!!$             end do
!!$          end do
          
          
          ! Read bins and copy samples over to posterior structure
          comm_br(id)%P(i)%nsamp   = nsamp
          comm_br(id)%P(i)%nfactor = 0
          do b = 1, comm_br(id)%P(i)%nbin
             read(unit,*) bin, status, d
             if (present(ell_min) .and. present(ell_max)) then
                if (bin(2) < ell_max .or. bin(1) > ell_min) then
                   cycle
                else
                   bin(1) = max(bin(1), ell_min)
                   bin(2) = min(bin(2), ell_max)
                end if
             end if

             ! Initialize covariance matrix structure
             n = 1
             pattern = .false.
             do j = 1, nmaps
                do k = j, nmaps
                   pattern(j,k) = status(n) /= '0' 
                   pattern(k,j) = pattern(j,k)
                   n            = n+1
                end do
             end do

             do while (any(pattern))
                comm_br(id)%P(i)%nfactor        = comm_br(id)%P(i)%nfactor + 1
                f                               = comm_br(id)%P(i)%nfactor 
                comm_br(id)%P(i)%factor(f)%lmin = bin(1)
                comm_br(id)%P(i)%factor(f)%lmax = bin(2)
                comm_br(id)%P(i)%factor(f)%d    = d
                
                ! Find which elements to include this time
                do col = 1, nmaps
                   if (any(pattern(1:nmaps,col))) exit
                end do
                
                ! Extract the appropriate segment
                p = count(pattern(1:nmaps,col))
                comm_br(id)%P(i)%factor(f)%p = p
                allocate(comm_br(id)%P(i)%factor(f)%ind(p), comm_br(id)%P(i)%factor(f)%status(p,p))
                n = 1
                do j = 1, nmaps
                   if (pattern(j,col)) then
                      comm_br(id)%P(i)%factor(f)%ind(n) = j
                      n                            = n+1
                   end if
                end do
                
                do c1 = 1, p
                   do c2 = c1, p
                      ind = comm_br(id)%P(i)%factor(f)%ind(c1)*(1-comm_br(id)%P(i)%factor(f)%ind(c1))/2 + &
                           & (comm_br(id)%P(i)%factor(f)%ind(c1)-1)*nmaps + &
                           & comm_br(id)%P(i)%factor(f)%ind(c2)
                      comm_br(id)%P(i)%factor(f)%status(c1,c2) = status(ind)
                   end do
                end do
                
                ! Mark current sub-matrix as completed
                pattern(comm_br(id)%P(i)%factor(f)%ind,comm_br(id)%P(i)%factor(f)%ind) = .false.
             end do
             
          end do

          ! Set up fiducial C_l's, and compute offsets
          do f = 1, comm_br(id)%P(i)%nfactor

             ! Allocate data structures
             p = comm_br(id)%P(i)%factor(f)%p
             allocate(comm_br(id)%P(i)%factor(f)%sigma(p,p,nsamp))
             allocate(comm_br(id)%P(i)%factor(f)%cl(p,p))
             allocate(comm_br(id)%P(i)%factor(f)%log_det_sigma(nsamp))
             
             ! Copy sigma samples into data structure
             lmin_bin                                = comm_br(id)%P(i)%factor(f)%lmin
             lmax_bin                                = comm_br(id)%P(i)%factor(f)%lmax
             comm_br(id)%P(i)%factor(f)%sigma(:,:,:) = 0.d0
             comm_br(id)%P(i)%factor(f)%n            = (lmax_bin+1)**2 - lmin_bin**2
             do k = 1, nsamp
                do c1 = 1, p
                   do c2 = c1, p
                      ind = comm_br(id)%P(i)%factor(f)%ind(c1)*(1-comm_br(id)%P(i)%factor(f)%ind(c1))/2 + &
                           & (comm_br(id)%P(i)%factor(f)%ind(c1)-1)*nmaps + &
                           & comm_br(id)%P(i)%factor(f)%ind(c2)
                      do l = lmin_bin, lmax_bin
                         comm_br(id)%P(i)%factor(f)%sigma(c1,c2,k) = &
                              & comm_br(id)%P(i)%factor(f)%sigma(c1,c2,k) + (2*l+1)*sigma_1D(l,ind,k)
                      end do
                      comm_br(id)%P(i)%factor(f)%sigma(c2,c1,k) = comm_br(id)%P(i)%factor(f)%sigma(c1,c2,k)
                   end do
                end do
                comm_br(id)%P(i)%factor(f)%log_det_sigma(k) = &
                     & log(comm_det(comm_br(id)%P(i)%factor(f)%sigma(:,:,k)))
             end do
             
             ! Set up fiducial C_l's
             call comm_bin_cls(cls, comm_br(id)%P(i)%factor(f)%lmin, &
                  & comm_br(id)%P(i)%factor(f)%lmax, &
                  & comm_br(id)%P(i)%factor(f)%ind, comm_br(id)%P(i)%factor(f)%cl)
          end do
          
          ! Compute offset
          comm_br(id)%P(i)%offset = comm_br_compute_offset(comm_br(id)%P(i))
          
          deallocate(sigma, sigma_1D, cls)
          
       end do
       close(unit)
       
    end if

  end subroutine comm_br_initialize_object


  subroutine comm_br_deallocate_object(handle)
    implicit none
    
    integer(i4b), optional :: handle

    integer(i4b) :: i, j, k, id, id_min, id_max

    id_min = 1; id_max = MAX_N_BR; 
    if (present(handle)) then 
       id_min = handle
       id_max = handle
    end if

    do id = id_min, id_max
       if (comm_br(id)%initialized) then
          do i = 1, comm_br(id)%n_p
             do j = 1, comm_br(id)%P(i)%nfactor
                if (allocated(comm_br(id)%P(i)%factor(j)%status))        deallocate(comm_br(id)%P(i)%factor(j)%status)
                if (allocated(comm_br(id)%P(i)%factor(j)%ind))           deallocate(comm_br(id)%P(i)%factor(j)%ind)
                if (allocated(comm_br(id)%P(i)%factor(j)%log_det_sigma)) deallocate(comm_br(id)%P(i)%factor(j)%log_det_sigma)
                if (allocated(comm_br(id)%P(i)%factor(j)%cl))            deallocate(comm_br(id)%P(i)%factor(j)%cl)
                if (allocated(comm_br(id)%P(i)%factor(j)%sigma))         deallocate(comm_br(id)%P(i)%factor(j)%sigma)
                comm_br(id)%P(i)%factor(j)%lmin = -1
                comm_br(id)%P(i)%factor(j)%lmax = -1
                comm_br(id)%P(i)%factor(j)%p    = -1
                comm_br(id)%P(i)%factor(j)%n    = -1
             end do

             comm_br(id)%P(i)%nsamp       = 0
             comm_br(id)%P(i)%nfactor     = 0
             comm_br(id)%P(i)%offset      = -1.6375d30
          end do
          if (allocated(comm_br(id)%P)) deallocate(comm_br(id)%P)
          comm_br(id)%n_p         = 0
          comm_br(id)%initialized = .false.
       end if
    end do

  end subroutine comm_br_deallocate_object


  function comm_br_compute_lnL(cls, ierr, handle, P_id, check_posdef)
    implicit none

    real(dp),     dimension(0:,1:), intent(in)           :: cls
    integer(i4b),                   intent(out)          :: ierr
    real(dp)                                             :: comm_br_compute_lnL
    integer(i4b),                   intent(in), optional :: handle, P_id
    logical(lgt),                   intent(in), optional :: check_posdef

    integer(i4b) :: i, l, id
    logical(lgt) :: posdef, check_posdef_
    real(dp)     :: lnL

    id = 1;                 if (present(handle))       id            = handle
    check_posdef_ = .true.; if (present(check_posdef)) check_posdef_ = check_posdef

    ! Check that likelihood structure is initialized
    if (.not. comm_br(id)%initialized) then
       write(*,*) 'Error -- comm_br_mod: Requested BR handle ', id, ' is not initialized'
       stop
    end if

    ! Check that spectrum is positive definite
    posdef = .true.
    do l = 2, size(cls,1)-1
       if (size(cls,2) == 1) then
          posdef = cls(l,1) >= 0.d0
       else
          posdef = cls(l,1)*cls(l,4)-cls(l,2)**2 >= 0.d0 .and. cls(l,6) >= 0.d0
       end if
       if (.not. posdef .and. check_posdef_) then
          comm_br_compute_lnL = -1.d30
          return
       end if
    end do

    ! Add contributions from all BR posteriors
    lnL  = 0.d0
    ierr = 0
    if (present(P_id)) then
       lnL = comm_br(id)%P(P_id)%loglike_weight * &
            & comm_br_compute_single_lnL(cls, comm_br(id)%P(P_id), ierr)       
    else
       do i = 1, comm_br(id)%n_p
          lnL = lnL + comm_br(id)%P(i)%loglike_weight * &
               & comm_br_compute_single_lnL(cls, comm_br(id)%P(i), ierr)       
          if (ierr /= 0) exit
       end do
    end if

    if (ierr == 0) then
       comm_br_compute_lnL = lnL
    else
       comm_br_compute_lnL = -1.d30
    end if

  end function comm_br_compute_lnL

  ! Base computation routine
  function comm_br_compute_single_lnL(cls, P, ierr)
    implicit none

    real(dp),            dimension(0:,1:), intent(in)   :: cls
    type(comm_P_struct),                   intent(in)   :: P
    integer(i4b),                          intent(out)  :: ierr
    real(dp)                                            :: comm_br_compute_single_lnL

    logical(lgt) :: posdef
    integer(i4b) :: i, j, l, m, n, d, f, k, id, num_active
    real(dp)     :: subtotal, x, lnL, log_det_Cl, Cl(3,3), lnL_hist
    real(dp), allocatable, dimension(:) :: subsum

    ! Compute the Blackwell-Rao estimator
    allocate(subsum(P%nsamp))
    subsum = 0.d0
    ierr = 0
    do f = 1, P%nfactor
       m = P%factor(f)%p
       n = P%factor(f)%n
       call comm_bin_cls(cls, P%factor(f)%lmin, P%factor(f)%lmax, P%factor(f)%ind, Cl(1:m,1:m))
       log_det_Cl = comm_det(Cl(1:m,1:m))
       if (log_det_Cl > 0.d0) then
          log_det_Cl = log(log_det_Cl)
          do k = 1, P%nsamp
             n        = P%factor(f)%n
             m        = P%factor(f)%p
             d        = P%factor(f)%d
             if (m == 2) then
                if (P%factor(f)%status(1,1)=='S' .and. P%factor(f)%status(1,2)=='S' .and. &
                     & .not. P%factor(f)%status(2,2)=='S') then
                   ! TT-TE marginalized over EE
                   subsum(k) = subsum(k) + 0.5d0 * ((n-1)*P%factor(f)%log_det_sigma(k) + &
                        & (n-2) * log(Cl(1,1)) - &
                        & (n-1) * log(P%factor(f)%sigma(1,1,k)*Cl(1,2)**2 - &
                        & 2.d0*P%factor(f)%sigma(1,2,k)*Cl(1,2)*Cl(1,1)+ &
                        & P%factor(f)%sigma(2,2,k)*Cl(1,1)**2) - &
                        & P%factor(f)%sigma(1,1,k)/Cl(1,1))
                else
                   subsum(k) = subsum(k) + 0.5d0 * ((n-1-2*d+m)*P%factor(f)%log_det_sigma(k) &
                        & - (n-2*(d-m)) * log_det_Cl &
                        & -  comm_tr_A_invB(P%factor(f)%sigma(:,:,k),Cl(1:m,1:m)))
                end if
             else
                subsum(k) = subsum(k) + 0.5d0 * ((n-1-2*d+m)*P%factor(f)%log_det_sigma(k) &
                     & - (n-2*(d-m)) * log_det_Cl &
                     & -  comm_tr_A_invB(P%factor(f)%sigma(:,:,k),Cl(1:m,1:m)))
             end if
             !subsum(k) = subsum(k) + 0.5d0 * ((n-1-m)*P%factor(f)%log_det_sigma(k) &
             !     & -  n * log_det_Cl &
             !     & -  comm_tr_A_invB(P%factor(f)%sigma(:,:,k),Cl(1:m,1:m)))
          end do
       else
          ierr = 1
          exit
       end if
    end do
    
    if (P%nfactor == 0) then
       lnL = 0.d0
    else if (ierr == 0) then
       subsum = subsum - P%offset
       num_active = count(subsum-maxval(subsum) > -4.605d0)
       if (num_active >= MIN_NUM_ACTIVE_SAMPLES) then
          lnL = sum(exp(subsum))
          if (lnL > 1d-20) then
!          write(*,*) 'Fraction of contributing samples = ', &
!              & count(subsum-maxval(subsum) > -4.605d0)/real(size(subsum),sp), &
!              & count(subsum-maxval(subsum) > -4.605d0)
             !write(*,*) lnL
             lnL = log(lnL)
          else
             lnL    = log(1d-30)
             ierr = 2
          end if
       else
          lnL  = log(1d-30)
          ierr = 3
       end if
    else
       lnL    = log(1d-30)
    end if

    comm_br_compute_single_lnL = lnL

    deallocate(subsum)

  end function comm_br_compute_single_lnL


  ! Utility routine for initializing the offset to be subtracted from each term
  ! to avoid overflow errors. Only called with the first power spectrum
  function comm_br_compute_offset(P)
    implicit none

    type(comm_P_struct), intent(in) :: P
    real(dp)                        :: comm_br_compute_offset

    integer(i4b) :: i, j, l, m, n, f, k, d
    real(dp)     :: subtotal, x, offset, inv_Cl(3,3)

    ! Compute the largest contribution to the Blackwell-Rao estimator from a given sample
    offset = -1.6375d30
    do k = 1, P%nsamp
       subtotal = 0.d0
       do f = 1, P%nfactor
          n        = P%factor(f)%n
          m        = P%factor(f)%p
          d        = P%factor(f)%d
          if (m == 2) then
             if (P%factor(f)%status(1,1)=='S' .and. P%factor(f)%status(1,2)=='S' .and. &
                  & .not. P%factor(f)%status(2,2)=='S') then
                subtotal = subtotal + 0.5d0 * ((n-1)*P%factor(f)%log_det_sigma(k) + &
                     & (n-2) * log(P%factor(f)%cl(1,1)) - &
                     & (n-1) * log(P%factor(f)%sigma(1,1,k)*P%factor(f)%cl(1,2)**2 - &
                     & 2.d0*P%factor(f)%sigma(1,2,k)*P%factor(f)%cl(1,2)*P%factor(f)%cl(1,1)+ &
                     & P%factor(f)%sigma(2,2,k)*P%factor(f)%cl(1,1)**2) - &
                     & P%factor(f)%sigma(1,1,k)/P%factor(f)%cl(1,1))
             else
                subtotal = subtotal + 0.5d0 * ((n-1-2*d+m)*P%factor(f)%log_det_sigma(k) &
                     & - (n-2*(d-m)) * log(comm_det(P%factor(f)%cl)) &
                     & -  comm_tr_A_invB(P%factor(f)%sigma(:,:,k),P%factor(f)%cl))
             end if
          else
             subtotal = subtotal + 0.5d0 * ((n-1-2*d+m)*P%factor(f)%log_det_sigma(k) &
                  & - (n-2*(d-m)) * log(comm_det(P%factor(f)%cl)) &
                  & -  comm_tr_A_invB(P%factor(f)%sigma(:,:,k),P%factor(f)%cl))
          end if
       end do
       offset = max(offset,subtotal)
    end do

    if (offset < -1.637e30) then
       write(*,*) 'Error -- comm_mod: Offset evaluation failed', offset
       stop
    endif

    comm_br_compute_offset = offset

  end function comm_br_compute_offset

end module comm_br_mod
