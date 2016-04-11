module comm_cr_mod
  use comm_comp_mod
  use comm_data_mod
  use comm_param_mod
  use comm_diffuse_comp_mod
  use comm_ptsrc_comp_mod
  use comm_template_comp_mod
  use rngmod
  use comm_cr_utils
  use comm_cr_precond_mod
  use math_tools
  use comm_output_mod
  implicit none

  private
  public solve_cr_eqn_by_CG, cr_amp2x, cr_x2amp, cr_computeRHS, cr_matmulA, cr_invM

  interface cr_amp2x
     module procedure cr_amp2x_full
  end interface cr_amp2x

  interface cr_x2amp
     module procedure cr_x2amp_full
  end interface cr_x2amp

contains

  subroutine solve_cr_eqn_by_CG(cpar, x, b, stat)
    implicit none

    type(comm_params),                intent(in)    :: cpar
    real(dp),          dimension(1:), intent(out)   :: x
    real(dp),          dimension(1:), intent(in)    :: b
    integer(i4b),                     intent(out)   :: stat

    integer(i4b) :: i, j, k, l, m, n, maxiter, root, ierr
    real(dp)     :: eps, tol, delta0, delta_new, delta_old, alpha, beta, t1, t2, t3, t4
    real(dp), allocatable, dimension(:)   :: Ax, r, d, q, temp_vec, s
    real(dp), allocatable, dimension(:,:) :: alm
    class(comm_comp),   pointer :: c

    root    = 0
    maxiter = cpar%cg_maxiter
    eps     = cpar%cg_tol
    n       = size(x)

    ! Allocate temporary data vectors
    allocate(Ax(n), r(n), d(n), q(n), s(n))

    ! Update preconditioner
    call wall_time(t1)
    call update_precond
    call wall_time(t2)
    if (cpar%myid == root .and. cpar%verbosity > 2) then
       write(*,fmt='(a,f8.2)') 'CG initialize preconditioner, time = ', real(t2-t1,sp)
    end if
    
!!$    if (cpar%myid == root) write(*,*) P_cr%invM_src(1,1)%M(1,:)
!!$    if (cpar%myid == root) write(*,*) P_cr%invM_src(1,1)%M(2,:)
!!$    if (cpar%myid == root) write(*,*)
!!$    if (cpar%myid == root) x    = 0.d0
!!$    if (cpar%myid == root) x(1) = 1.d0
!!$    q     = cr_matmulA(x)
!!$    if (cpar%myid == root) write(*,*) q
!!$    if (cpar%myid == root) x    = 0.d0
!!$    if (cpar%myid == root) x(2) = 1.d0
!!$    q     = cr_matmulA(x)
!!$    if (cpar%myid == root) write(*,*) q
!!$
!!$    call mpi_finalize(ierr)
!!$    stop
!!$    
    ! Initialize the CG search
    x  = 0.d0
    r  = b ! - A(x)   ! x is zero
    d  = cr_invM(r)

    delta_new = mpi_dot_product(cpar%comm_chain,r,d)
    delta0    = delta_new
    do i = 1, maxiter
       call wall_time(t1)
       
       if (delta_new < eps * delta0) exit

       call wall_time(t3)
       q     = cr_matmulA(d)
       call wall_time(t4)
       !if (cpar%myid == root .and. cpar%verbosity > 2) write(*,fmt='(a,f8.2)') 'A time = ', real(t4-t3,sp)
       alpha = delta_new / mpi_dot_product(cpar%comm_chain, d, q)
       x     = x + alpha * d

       ! Restart every 50th iteration to suppress numerical errors
       if (.false. .and. mod(i,50) == 0) then
          r = b - cr_matmulA(x)
       else
          r = r - alpha*q
       end if

       call wall_time(t3)
       s         = cr_invM(r)
       call wall_time(t4)
       !if (cpar%myid == root .and. cpar%verbosity > 2) write(*,fmt='(a,f8.2)') 'invM time = ', real(t4-t3,sp)
       delta_old = delta_new 
       delta_new = mpi_dot_product(cpar%comm_chain, r, s)
       beta      = delta_new / delta_old
       d         = s + beta * d

       if (cpar%output_cg_freq > 0) then
          if (mod(i,cpar%output_cg_freq) == 0) then
             ! Multiply with sqrt(S)     
             c       => compList
             do while (associated(c))
                select type (c)
                class is (comm_diffuse_comp)
                   if (trim(c%cltype) /= 'none') then
                      allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
                      call cr_extract_comp(c%id, x, alm)
                      call c%Cl%sqrtS(alm=alm, info=c%x%info) ! Multiply with sqrt(Cl)  
                      call cr_insert_comp(c%id, .false., alm, x)
                      deallocate(alm)
                   end if
                end select
                c => c%next()
             end do
             call cr_x2amp(x)
             call output_FITS_sample(cpar, i)
          end if
       end if

       if (cpar%myid == root) write(*,*) x(size(x)-1:size(x))

       call wall_time(t2)
       if (cpar%myid == root .and. cpar%verbosity > 2) then
          write(*,fmt='(a,i5,a,e13.5,a,e13.5,a,f8.2)') 'CG iter. ', i, ' -- res = ', &
               & real(delta_new,sp), ', tol = ', real(eps * delta0,sp), &
               & ', time = ', real(t2-t1,sp)
       end if

    end do

    ! Multiply with sqrt(S), and insert into right object
    c       => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, x, alm)
             call c%Cl%sqrtS(alm=alm, info=c%x%info) ! Multiply with sqrt(Cl)
             call cr_insert_comp(c%id, .false., alm, x)
             deallocate(alm)
          end if
       end select
       c => c%next()
    end do

    if (i >= maxiter) then
       write(*,*) 'ERROR: Convergence in CG search not reached within maximum'
       write(*,*) '       number of iterations = ', maxiter
       stat = stat + 1
    else
       if (cpar%myid == root .and. cpar%verbosity > 1) then
          write(*,fmt='(a,i5,a,e13.5,a,e13.5,a,f8.2)') 'Final CG iter ', i-1, ' -- res = ', &
               & real(delta_new,sp), ', tol = ', real(eps * delta0,sp)
       end if
    end if

    deallocate(Ax, r, d, q, s)
    
  end subroutine solve_cr_eqn_by_CG

  function cr_amp2x_full()
    implicit none

    real(dp), allocatable, dimension(:) :: cr_amp2x_full

    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c

    ! Stack parameters linearly
    allocate(cr_amp2x_full(ncr))
    ind = 1
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp) 
          do i = 1, c%x%info%nmaps
             cr_amp2x_full(ind:ind+c%x%info%nalm-1) = c%x%alm(:,i)
             ind = ind + c%x%info%nalm
          end do
       class is (comm_ptsrc_comp) 
          do i = 1, c%nmaps
             cr_amp2x_full(ind:ind+c%nsrc-1) = c%x(:,i)
             ind = ind + c%nsrc
          end do
       end select
       c => c%next()
    end do

  end function cr_amp2x_full

  subroutine cr_x2amp_full(x)
    implicit none

    real(dp), dimension(:), intent(in) :: x

    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c

    ind = 1
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          do i = 1, c%x%info%nmaps
             c%x%alm(:,i) = x(ind:ind+c%x%info%nalm-1)
             ind = ind + c%x%info%nalm
          end do
       class is (comm_ptsrc_comp)
          do i = 1, c%nmaps
             if (c%myid == 0) then
                c%x(:,i) = x(ind:ind+c%nsrc-1)
                ind = ind + c%nsrc
             end if
          end do
       end select
       c => c%next()
    end do
    
  end subroutine cr_x2amp_full

  ! ---------------------------
  ! Definition of linear system
  ! ---------------------------

  subroutine cr_computeRHS(handle, rhs)
    implicit none

    type(planck_rng),                            intent(inout)          :: handle
    real(dp),         allocatable, dimension(:), intent(out)            :: rhs

    integer(i4b) :: i, j, l, m, k, n, ierr
    class(comm_map),     pointer                 :: map, Tm, mu
    class(comm_comp),    pointer                 :: c
    class(comm_mapinfo), pointer                 :: info
    real(dp),        allocatable, dimension(:,:) :: eta, Tp

    call rand_init(handle, 10)
    
    ! Initialize output vector
    allocate(rhs(ncr))
    rhs = 0.d0

    ! Add channel dependent terms
    do i = 1, numband

       ! Set up Wiener filter term
       map => comm_map(data(i)%map)
       call data(i)%N%sqrtInvN(map)

       ! Add channel-dependent white noise fluctuation
       do k = 1, map%info%nmaps
          do j = 0, map%info%np-1
!             map%map(j,k) = map%map(j,k) + rand_gauss(handle)
          end do
       end do

       ! Multiply with sqrt(invN)
       call data(i)%N%sqrtInvN(map)

       ! Convolve with transpose beam
       call map%Yt()
       call data(i)%B%conv(alm_in=.false., alm_out=.false., trans=.true., map=map)

!!$       call map%Y()
!!$       call map%writeFITS('test.fits')
!!$       call mpi_finalize(ierr)
!!$       stop
       
       ! Multiply with (transpose and component specific) mixing matrix, and
       ! insert into correct segment
       c => compList
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             info  => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, c%lmax_amp, &
                  & c%nmaps, data(i)%info%pol)
             Tm     => comm_map(info)
             if (c%F_null(i)) then
                Tm%alm = 0.d0
             else
                call map%alm_equal(Tm)
                if (c%lmax_ind == 0) then
                   do j = 1, min(c%x%info%nmaps, Tm%info%nmaps)
                      Tm%alm(:,j) = Tm%alm(:,j) * c%F_mean(i,j)
                   end do
                else
                   call Tm%Y()
                   Tm%map = c%F(i)%p%map * Tm%map
                   call Tm%YtW()
                end if
             end if
             call c%Cl%sqrtS(map=Tm) ! Multiply with sqrt(Cl)
             call cr_insert_comp(c%id, .true., Tm%alm, rhs)
             deallocate(Tm)
          class is (comm_ptsrc_comp)
             allocate(Tp(c%nsrc,c%nmaps))
             Tp = c%projectBand(i,map)
             if (c%myid == 0) call cr_insert_comp(c%id, .true., Tp, rhs)
             deallocate(Tp)
          end select
          c => c%next()
       end do

       deallocate(map)
    end do

!!$    write(*,*) real(rhs(1:100),sp)
!!$    call mpi_finalize(ierr)
!!$    stop

    ! Add prior terms
    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             n = ind_comp(c%id,2)
             allocate(eta(0:c%x%info%nalm-1,c%x%info%nmaps))
             eta = 0.d0
             ! Variance term
             do j = 1, c%x%info%nmaps
                do i = 0, c%x%info%nalm-1
!                   eta(i,j) = rand_gauss(handle)
                end do
             end do
             ! Mean term
             if (associated(c%mu)) then
                mu => comm_map(c%mu)
                call c%Cl%sqrtInvS(map=mu)
                eta = eta + mu%alm
                call mu%dealloc()
             end if
             call cr_insert_comp(c%id, .true., eta, rhs)
             deallocate(eta)
          end if
       end select
       c => c%next()
    end do
    nullify(c)

!    write(*,*) real(rhs(1:100),sp)

  end subroutine cr_computeRHS

  recursive function cr_matmulA(x)
    implicit none

    real(dp), dimension(1:),     intent(in)  :: x
    real(dp), dimension(size(x))             :: cr_matmulA

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, myid
    class(comm_map),  pointer :: map, pmap
    class(comm_comp), pointer :: c
    real(dp),        allocatable, dimension(:)   :: y, sqrtS_x
    real(dp),        allocatable, dimension(:,:) :: alm, m, pamp

    ! Initialize output array
    allocate(y(ncr), sqrtS_x(ncr))
    y = 0.d0
    myid = data(1)%map%info%myid

    ! Multiply with sqrt(S)
    call wall_time(t1)
    sqrtS_x = x
    c       => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, sqrtS_x, alm)
             call c%Cl%sqrtS(alm=alm, info=c%x%info) ! Multiply with sqrt(Cl)
             call cr_insert_comp(c%id, .false., alm, sqrtS_x)
             deallocate(alm)
          end if
       end select
       c => c%next()
    end do
    call wall_time(t2)
    !if (myid == 0) write(*,fmt='(a,f8.2)') 'sqrtS time = ', real(t2-t1,sp)

    
    ! Add frequency dependent terms
    do i = 1, numband

       ! Compute component-summed map, ie., column-wise matrix elements
       call wall_time(t1)
       map  => comm_map(data(i)%info)   ! For diffuse components
       pmap => comm_map(data(i)%info)   ! For point-source components
       c   => compList
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             call cr_extract_comp(c%id, sqrtS_x, alm)
             allocate(m(0:c%x%info%nalm-1,c%x%info%nmaps))
             m = c%getBand(i, amp_in=alm, alm_out=.true.)
             call map%add_alm(m, c%x%info)
             deallocate(alm, m)
          class is (comm_ptsrc_comp)
             call cr_extract_comp(c%id, sqrtS_x, pamp)
             allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
             m = c%getBand(i, amp_in=pamp)
             pmap%map = pmap%map + m
             deallocate(pamp, m)
          end select
          c => c%next()
       end do
       call map%Y()                    ! Diffuse components
       map%map = map%map + pmap%map    ! Add compact objects
       call wall_time(t2)
       !if (myid == 0) write(*,fmt='(a,f8.2)') 'getBand time = ', real(t2-t1,sp)

       ! Multiply with invN
       call wall_time(t1)
       call data(i)%N%InvN(map)
       call wall_time(t2)
       !if (myid == 0) write(*,fmt='(a,f8.2)') 'invN time = ', real(t2-t1,sp)

       ! Project summed map into components, ie., row-wise matrix elements
       call wall_time(t1)
       c   => compList
       call map%Yt()             ! Prepare for diffuse components
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             alm = c%projectBand(i, map, alm_in=.true.)
             call cr_insert_comp(c%id, .true., alm, y)
             deallocate(alm)
          class is (comm_ptsrc_comp)
             allocate(pamp(0:c%nsrc-1,c%nmaps))
             pamp = c%projectBand(i, map)
             if (c%myid == 0) call cr_insert_comp(c%id, .true., pamp, y)
             deallocate(pamp)
          end select
          c => c%next()
       end do
       call wall_time(t2)
       !if (myid == 0) write(*,fmt='(a,f8.2)') 'projBand time = ', real(t2-t1,sp)

       deallocate(map,pmap)
       
    end do

    ! Add prior term and multiply with sqrt(S) for relevant components
    call wall_time(t1)
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             ! Multiply with sqrt(Cl)
             call cr_extract_comp(c%id, y, alm)
             call c%Cl%sqrtS(alm=alm, info=c%x%info)
             call cr_insert_comp(c%id, .false., alm, y)
             ! Add (unity) prior term
             call cr_extract_comp(c%id, x, alm)
             call cr_insert_comp(c%id, .true., alm, y)
             deallocate(alm)
          end if
       end select
       c => c%next()
    end do
    nullify(c)
    call wall_time(t2)
    !if (myid == 0) write(*,fmt='(a,f8.2)') 'prior time = ', real(t2-t1,sp)
    !call mpi_finalize(i)
    !stop

    ! Return result and clean up
    call wall_time(t1)
    cr_matmulA = y
    deallocate(y, sqrtS_x)
    call wall_time(t2)
    !    write(*,*) 'f', t2-t1
    
  end function cr_matmulA

  function cr_invM(x)
    implicit none
    real(dp),              dimension(:), intent(in) :: x
    real(dp), allocatable, dimension(:)             :: cr_invM

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:,:) :: alm
    class(comm_comp), pointer :: c

    if (.not. allocated(cr_invM)) allocate(cr_invM(size(x)))
    cr_invM = x
    call applyDiffPrecond(cr_invM)
    call applyPtsrcPrecond(cr_invM)
    
  end function cr_invM

  subroutine update_precond
    implicit none

    call updateDiffPrecond
    call updatePtsrcPrecond

  end subroutine update_precond
  
  
end module comm_cr_mod
