module comm_mcmc_mod
  use comm_mp_mod
  use ARS_mod
  use comm_Cl_util_mod
  implicit none

  real(dp),                                  private :: joint_threshold
  integer(i4b),                              private :: comm_master, myid_master, root, numprocs
  integer(i4b),                              private :: nstep
  integer(i4b),                              private :: verbosity
  integer(i4b),                              private :: lmin_bin, lmax_bin, spec_bin
  character(len=128),                        private :: operation, cdir

  real(dp),           allocatable, dimension(:,:),   private :: cls_prop

  type(genvec), private :: s_prop


contains

  subroutine initialize_mcmc_mod(paramfile, comm_master_in)
    implicit none

    integer(i4b),                            intent(in) :: comm_master_in
    character(len=128),                      intent(in) :: paramfile

    integer(i4b)       :: i, j, k, l, l1, l2, m, n, ierr
    real(dp)           :: mcmc_scale_fac, noise, cmb, a, b, c
    character(len=128) :: filename
    logical(lgt)       :: polarization, cl_binning, enforce_zero_cl
    real(dp),           allocatable, dimension(:)       :: input_vals
    real(dp),           allocatable, dimension(:,:,:,:) :: stddev
    character(len=128), allocatable, dimension(:,:,:)   :: prop_type

    comm_master = comm_master_in
    call mpi_comm_rank(MPI_COMM_WORLD, myid_master, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

    call get_parameter(paramfile, 'OPERATION',                     par_string=operation)
    call get_parameter(paramfile, 'VERBOSITY',                     par_int=verbosity)
    call get_parameter(paramfile, 'NUM_JOINT_ALMS_CL_STEPS',       par_int=nstep)
    call get_parameter(paramfile, 'JOINT_ALMS_CL_CHISQ_THRESHOLD', par_dp=joint_threshold)
    call get_parameter(paramfile, 'CHAIN_DIRECTORY',               par_string=cdir)

  end subroutine initialize_mcmc_mod


  subroutine cleanup_mcmc_mod
    implicit none

  end subroutine cleanup_mcmc_mod



  subroutine sample_cls_and_alms_by_mcmc(rng_handle, cl, s)  
    implicit none

    type(planck_rng),                   intent(inout) :: rng_handle
    real(dp),         dimension(0:,1:), intent(inout) :: cl
    type(genvec),                       intent(inout) :: s

    real(dp)     :: cl_bounds(2), my_cl, cl_init(3), cl_min, chisq0, chisq1, chisq_test
    integer(i4b) :: i, j, k, l, l1, l2, m, ind, bin, iter, niter, neval, status, unit
    integer(i4b), save :: numcall = 0
    logical(lgt), allocatable, dimension(:,:), save :: sample_joint

    if (nstep == 0) return

    numcall = numcall+1
    if (numcall == 1) then
       allocate(sample_joint(num_Cl_bin, nspec))
       do i = 1, num_cl_bin
          do j = 1, nspec
             sample_joint(i,j) = Cl_bin_stat(i,j) == 'S' .or. Cl_bin_stat(i,j) == 'M'
          end do
       end do
    end if

    call allocate_genvec(s_prop)
    allocate(cls_prop(0:size(cl,1)-1, size(cl,2)))
    cls_prop = cl
    s_prop   = s

    do iter = 1, nstep
       do i = 1, num_Cl_bin
          do j = 1, nspec
             l = Cl_bins(i,1)
             if (sample_joint(i,j)) then
                if (j == 1) then
                   cl_init   = cls_prop(l,j) * [0.5d0, 1.d0, 2.d0]
                   cl_min    = 0.d0
                   if (cls_prop(l,2) /= 0.d0) cl_min = max(cl_min, cls_prop(l,2)**2/cls_prop(l,4))
                   if (cls_prop(l,3) /= 0.d0) cl_min = max(cl_min, cls_prop(l,3)**2/cls_prop(l,6))
                   cl_bounds = [cl_min, 2.d4]
                else if (j == 2) then
                   cl_init   = cls_prop(l,j) + [-1.d0, 0.d0, 1.d0]
                   cl_bounds = [-sqrt(cls_prop(l,1)*cls_prop(l,4)), sqrt(cls_prop(l,1)*cls_prop(l,4))]
                else if (j == 3) then
                   cl_init   = cls_prop(l,j) + 0.1d0*[-1.d0, 0.d0, 1.d0]
                   cl_bounds = [-sqrt(cls_prop(l,1)*cls_prop(l,6)), sqrt(cls_prop(l,1)*cls_prop(l,6))]
                else if (j == 4) then
                   cl_init   = cls_prop(l,j) * [0.5d0, 1.d0, 2.d0]
                   cl_min = 1.d-6 !0.d0
                   if (cls_prop(l,2) /= 0.d0) cl_min = max(cl_min, cls_prop(l,2)**2/cls_prop(l,1))
                   if (cls_prop(l,5) /= 0.d0) cl_min = max(cl_min, cls_prop(l,3)**2/cls_prop(l,6))
                   cl_bounds = [cl_min, 2.d0]
                else if (j == 5) then
                   cl_init   = cls_prop(l,j) + 0.1d0*[-1.d0, 0.d0, 1.d0]
                   cl_bounds = [-sqrt(cls_prop(l,4)*cls_prop(l,6)), sqrt(cls_prop(l,4)*cls_prop(l,6))]
                else if (j == 6) then
                   cl_init   = cls_prop(l,j) * [0.5d0, 1.d0, 2d0]
                   cl_min = 1.d-6 !0.d0
                   if (cls_prop(l,3) /= 0.d0) cl_min = max(cl_min, cls_prop(l,2)**2/cls_prop(l,1))
                   if (cls_prop(l,5) /= 0.d0) cl_min = max(cl_min, cls_prop(l,3)**2/cls_prop(l,4))
                   cl_bounds = [cl_min, 2.d0]
                end if
                if (cl_init(1) <= cl_bounds(1)) cl_init(1) = 0.5d0*(cls_prop(l,j) + cl_bounds(1))
                if (cl_init(3) >= cl_bounds(2)) cl_init(3) = 0.5d0*(cls_prop(l,j) + cl_bounds(2))
                !cl_init(1) = max(cl_init(1), cl_bounds(1)+1.d-7)
                !cl_init(3) = min(cl_init(3), cl_bounds(2)-1.d-7)

                lmin_bin      = Cl_bins(i,1)
                lmax_bin      = Cl_bins(i,2)
                spec_bin      = j

                if (lmin_bin == 2 .and. spec_bin == 4 .and. numcall == 1) then
                   open(58,file='p.dat')
                   do m = 1, 1000
                      my_cl = cl_bounds(1) + (cl_bounds(2)-cl_bounds(1))/(1000.-1.d0) * (m-1)
                      !my_cl = cl_bounds(1) + (1d-2-cl_bounds(1))/(1000.-1.d0) * (m-1)
                      write(58,*) my_cl, lnL_joint_alms_Cl(my_cl)
                   end do
                   close(58)
                   call mpi_finalize(m)
                   stop

                   write(*,*) 'hei'
                   open(58,file='samp2.dat')
                   do m = 1, 1000
                      write(*,*) m
                      my_cl = sample_ARS(rng_handle, lnL_joint_alms_Cl, cl_init, cl_bounds, &
                           & neval=neval, status=status)
                      write(58,*) my_cl, -2*lnL_joint_alms_Cl(my_cl)
                   end do
                   close(58)
                   call mpi_finalize(m)
                   stop
                end if

                if (numcall == 1) then
                   chisq0 = lnL_joint_alms_Cl(cls_prop(l,j))
                   chisq1 = lnL_joint_alms_Cl(max(cl_bounds(1), 0.d0))
                   if (chisq1-chisq0 > joint_threshold) then
                      sample_joint(i,j) = .false.
                      cycle
                   end if
                end if

                !write(*,*) lmin_bin, spec_bin, real(chisq0-chisq1,sp)

!                if (lmin_bin == 2 .and. spec_bin == 6) then
!                   write(*,*) real(cls_prop(lmin_bin,:),sp)
!                   write(*,*) cl_bounds
!                   write(*,*) 'init', real(cl_init,sp)
!!$                   open(58,file='p.dat')
!!$                   do m = 1, 1000
!!$                      write(58,*) 0.1d0*m, lnL_joint_alms_Cl(0.1d0*m)
!!$                   end do
!!$                   close(58)
                   !stop
!                end if

!!$                 write(*,*) '********************************************************'
!!$                write(*,*) 'init  ', real(cl_init,sp)
!!$                write(*,*) 'bounds', real(cl_bounds,sp)
                my_cl         = sample_ARS(rng_handle, lnL_joint_alms_Cl, cl_init, cl_bounds, &
                     & neval=neval, status=status)
                write(*,*) lmin_bin, lmax_bin, spec_bin,my_cl, status
!                if (spec_bin == 6 .and. lmin_bin == 2) then
                   call mpi_finalize(i)
                   stop
!                end if
                if (status == 0) then
!                   write(*,*) 'success'
!                   write(*,*) lmin_bin, spec_bin, real(cls_prop(i,j),sp), real(my_cl,sp), neval
!!$                if (lmin_bin == 4 .and. spec_bin == 2) then
!!$                   open(58,file='samp.dat')
!!$                   do m = 1, 1000
!!$                      write(58,*) sample_ARS(rng_handle, lnL_joint_alms_Cl, cl_init, cl_bounds)
!!$                   end do
!!$                   close(58)
!!$                   stop
!!$                end if
                   call rescale_alms(lmin_bin, lmax_bin, spec_bin, my_cl, s_prop%cmb_amp)
                   cls_prop(lmin_bin:lmax_bin,j) = my_cl
                   call compute_chisq(map_id, s_prop, chisq_highlat=chisq_test)
!                   write(*,*) 'chisq etter = ', chisq_test
                else
!                   write(*,*) 'failed', status
                end if
             end if
          end do
       end do
    end do

    if (numcall == 1 .and. myid_master == root) then
       unit = getlun()
       open(unit, file=trim(cdir)//'/joint_sampling_scheme.dat')
       do i = 1, num_Cl_bin
          write(unit,*) Cl_bins(i,:), sample_joint(i,:)
       end do
       close(unit)
    end if

    cl = cls_prop
    call genvec_set_equal(s_prop, s)

    call deallocate_genvec(s_prop)
    deallocate(cls_prop)

  end subroutine sample_cls_and_alms_by_mcmc

  function lnL_joint_alms_Cl(Cl)
    use healpix_types
    implicit none
    real(dp), intent(in) :: Cl
    real(dp)             :: lnL_joint_alms_Cl

    integer(i4b) :: l, m, ind
    real(dp)     :: chisq
    type(genvec) :: s_int

!    if (Cl == 0.d0) then
!       lnL_joint_alms_Cl = -1.d30
!       return
!    end if

    call allocate_genvec(s_int)
    call genvec_set_equal(s_prop, s_int)
    !call compute_chisq(map_id, s_int, chisq_highlat=chisq)
    call rescale_alms(lmin_bin, lmax_bin, spec_bin, Cl, s_int%cmb_amp)
    call compute_chisq(map_id, s_int, chisq_highlat=chisq)

    !write(*,*) chisq, ((lmax_bin+1)**2-lmin_bin**2)*log(Cl/cls_prop(lmin_bin,spec_bin)), ((lmax_bin+1)**2-lmin_bin**2)
    lnL_joint_alms_Cl = -0.5d0 * chisq !- 0.5d0 * ((lmax_bin+1)**2-lmin_bin**2)*log(Cl)
!!$    do l = lmin_bin, lmax_bin
!!$       do m = -l, l
!!$          ind = lm2ind(l,m)
!!$          lnL_joint_alms_Cl = lnL_joint_alms_Cl - 0.5d0 * s_prop%cmb_amp(ind,2)**2 / Cl - 0.5d0 * log(Cl)
!!$       end do
!!$    end do

    call deallocate_genvec(s_int)

  end function lnL_joint_alms_Cl

  subroutine rescale_alms(lmin, lmax, spec, Cl, alms)
    implicit none

    integer(i4b),                   intent(in)    :: lmin, lmax, spec
    real(dp),                       intent(in)    :: Cl
    real(dp),     dimension(1:,1:), intent(inout) :: alms

    integer(i4b) :: i, l, m, ind(3), n
    real(dp),     dimension(nspec)       :: this_cl
    real(dp),     allocatable, dimension(:,:) :: mat
    real(dp),     dimension(nmaps,nmaps) :: inv_sqrt_S, sqrt_S_prop

    do l = lmin, lmax
       this_cl       = cls_prop(lmin,:)
       call cl2s(this_cl, inv_sqrt_S) ! Don't need the l(l+1)/2pi factor, since it cancels

       n = 0
       do i = 1, nmaps
          if (inv_sqrt_S(i,i) > 0.d0) then
             n = n+1
             ind(n) = i
          end if
       end do
       allocate(mat(n,n))
       mat = inv_sqrt_S(ind(1:n),ind(1:n))
       call compute_hermitian_root(mat, -0.5d0)
       inv_sqrt_S(ind(1:n),ind(1:n)) = mat
       this_cl(spec) = Cl
       call cl2s(this_cl, sqrt_S_prop) 
       mat = sqrt_S_prop(ind(1:n),ind(1:n))
       call compute_hermitian_root(mat, 0.5d0)
       sqrt_S_prop(ind(1:n),ind(1:n)) = mat
       do m = -l, l
          i = lm2ind(l,m)
          alms(i,:) = matmul(sqrt_S_prop, matmul(inv_sqrt_S, alms(i,:)))
       end do
       deallocate(mat)
    end do

  end subroutine rescale_alms

end module comm_mcmc_mod
