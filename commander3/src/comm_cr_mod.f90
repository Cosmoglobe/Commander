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
module comm_cr_mod
  use comm_output_mod
  implicit none

!  private
!  public solve_cr_eqn_by_CG, cr_amp2x, cr_x2amp, cr_computeRHS, cr_matmulA, cr_invM

  interface cr_amp2x
     module procedure cr_amp2x_full
  end interface cr_amp2x

  interface cr_x2amp
     module procedure cr_x2amp_full
  end interface cr_x2amp

contains

  subroutine solve_cr_eqn_by_CG(cpar, samp_group, x, b, stat)
    implicit none
    type(comm_params),                intent(in)    :: cpar
    integer(i4b),                     intent(in)    :: samp_group
    real(dp),          dimension(1:), intent(out)   :: x
    real(dp),          dimension(1:), intent(in)    :: b
    integer(i4b),                     intent(out)   :: stat

    integer(i4b) :: i, j, k, l, m, n, maxiter, root, ierr
    integer(i4b), save :: samp_group_prev
    real(dp)     :: eps, tol, delta0, delta_new, delta_old, alpha, beta, t1, t2, t3, t4
    real(dp)     :: lim_convergence, val_convergence, chisq, chisq_prev, buff
    real(dp), allocatable, dimension(:)   :: Ax, r, d, q, temp_vec, s, x_out
    real(dp), allocatable, dimension(:,:) :: alm, pamp
    class(comm_comp),   pointer :: c => null()

    root    = 0
    maxiter = cpar%cg_samp_group_maxiter(samp_group)
    eps     = cpar%cg_tol
    n       = size(x)

    ! Allocate temporary data vectors
    call update_status(status, "cr1")
    allocate(Ax(n), r(n), d(n), q(n), s(n))

    ! Update preconditioner
    call wall_time(t1)
    !call update_status(status, "cr2")
    call update_precond(samp_group, samp_group /= samp_group_prev)
    samp_group_prev = samp_group
    !call update_status(status, "cr3")
    call wall_time(t2)
    if (cpar%myid_chain == root .and. cpar%verbosity > 2) then
       write(*,fmt='(a,f8.2)') ' |  CG initialize preconditioner, time = ', real(t2-t1,sp)
    end if


    ! Initialize the CG search
    if (cpar%cg_init_zero) then
       x  = 0.d0
    else 
       call cr_amp2x(samp_group, x)
       ! Multiply with sqrt(invS)
       c       => compList
       do while (associated(c))
          if (.not. c%active_samp_group(samp_group)) then
             c => c%nextComp()
             cycle
          end if
          select type (c)
          class is (comm_diffuse_comp)
             if (trim(c%cltype) /= 'none') then
                allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
                call cr_extract_comp(c%id, x, alm)
                call c%Cl%sqrtInvS(alm=alm, info=c%x%info) ! Multiply with sqrt(inv(Cl))
                call cr_insert_comp(c%id, .false., alm, x)
                deallocate(alm)
             end if
          class is (comm_ptsrc_comp)
             if (c%myid == 0 .and. .not. c%precomputed_amps) then
                call cr_extract_comp(c%id, x, pamp)
                do j = 1, c%nmaps
                   do i = 1, c%nsrc
                      pamp(i-1,j) = pamp(i-1,j) / c%src(i)%P_x(j,2) ! Multiply with sqrtInvS
                   end do
                end do
                call cr_insert_comp(c%id, .false., pamp, x)
                deallocate(pamp)
             end if
          class is (comm_template_comp)
             if (c%myid == 0) then
                call cr_extract_comp(c%id, x, pamp)
                pamp = pamp / c%P_cg(2) ! Multiply with sqrtInvS
                call cr_insert_comp(c%id, .false., pamp, x)
                deallocate(pamp)
             end if
          end select
          c => c%nextComp()
       end do
    end if
!!$    eps = 0
!!$    do i = 1, 100000
!!$       call update_status(status, "cr4")
!!$       r  = cr_matmulA(x, samp_group) 
!!$       c   => compList
!!$       do while (associated(c))
!!$          if (.not. c%active_samp_group(samp_group)) then
!!$             c => c%nextComp()
!!$             cycle
!!$          end if
!!$          select type (c)
!!$          class is (comm_diffuse_comp)
!!$             allocate(mp(0:data(1)%info%nalm-1,data(1)%info%nmaps))
!!$             mp = c%getBand(1, amp_in=data(1)%map%alm, alm_out=.true.)
!!$             eps = eps + sum(abs(mp))
!!$             deallocate(mp)
!!$          end select
!!$          c => c%nextComp()
!!$       end do
!!$    end do
!!$    write(*,*) eps
!!$    call mpi_finalize(ierr)
!!$    stop

    call update_status(status, "cr4")
    r  = b - cr_matmulA(x, samp_group)   ! x is zero
    call update_status(status, "cr5")
    d  = cr_invM(cpar%comm_chain, r, samp_group)
    call update_status(status, "cr6")

    delta_new = mpi_dot_product(cpar%comm_chain,r,d)
    call update_status(status, "cr7")
    delta0    = mpi_dot_product(cpar%comm_chain,b,cr_invM(cpar%comm_chain, b, samp_group))
    call update_status(status, "cr8")
    if (delta0 > 1d30) then
       if(cpar%myid == root) then
          write(*,*) 'CR warning: Large initial residual = ', delta0
       end if
!!$       call mpi_finalize(ierr)
!!$       stop
    end if

    ! Set up convergence criterion
    if (trim(cpar%cg_conv_crit) == 'residual' .or. trim(cpar%cg_conv_crit) == 'fixed_iter') then
       lim_convergence = eps*delta0
       val_convergence = 1.d2*lim_convergence
    else if (trim(cpar%cg_conv_crit) == 'chisq') then
       lim_convergence = eps
       val_convergence = 1.d0
       call cr_compute_chisq(cpar%comm_chain, samp_group, x, chisq)
    else
       write(*,*) 'Error: Unsupported convergence criterion = ', trim(cpar%cg_conv_crit)
    end if
    do i = 1, maxiter
       call wall_time(t1)

       call update_status(status, "cg1")
       
       ! Check convergence
       if (mod(i,cpar%cg_check_conv_freq) == 0) then
          if (trim(cpar%cg_conv_crit) == 'residual' .or. trim(cpar%cg_conv_crit) == 'fixed_iter') then
             val_convergence = delta_new
          else if (trim(cpar%cg_conv_crit) == 'chisq') then
             chisq_prev      = chisq
             call cr_compute_chisq(cpar%comm_chain, samp_group, x, chisq)
             val_convergence = abs((chisq_prev-chisq)/chisq)
          end if
          if (val_convergence < lim_convergence .and. &
               & (i >= cpar%cg_miniter .or. delta_new <= 1d-40 * delta0) .and. &
               & trim(cpar%cg_conv_crit) /= 'fixed_iter') exit
          if (delta_new <= 1d-40 * delta0 .and. &
               & trim(cpar%cg_conv_crit) == 'fixed_iter') exit
       end if
       
       call update_status(status, "cg2")
   
       !if (delta_new < eps * delta0 .and. (i >= cpar%cg_miniter .or. delta_new <= 1d-30 * delta0)) exit

       q     = cr_matmulA(d, samp_group)
       alpha = delta_new / mpi_dot_product(cpar%comm_chain, d, q)
       x     = x + alpha * d

       ! Restart every 50th iteration to suppress numerical errors
       if (.false. .and. mod(i,50) == 0) then
          r = b - cr_matmulA(x, samp_group)
       else
          r = r - alpha*q
       end if

      call update_status(status, "cg3")
       call wall_time(t3)
       s         = cr_invM(cpar%comm_chain, r, samp_group)
       call wall_time(t4)
       !if (cpar%myid == root .and. cpar%verbosity > 2) write(*,fmt='(a,f8.2)') 'invM time = ', real(t4-t3,sp)
       delta_old = delta_new 
       delta_new = mpi_dot_product(cpar%comm_chain, r, s)
       beta      = delta_new / delta_old
       d         = s + beta * d
       call update_status(status, "cg4")

!call mpi_finalize(ierr)
!stop

       if (cpar%output_cg_freq > 0) then
          if (mod(i,cpar%output_cg_freq) == 0) then
             ! Multiply with sqrt(S)     
             allocate(x_out(n))
             x_out = x
             c       => compList
             do while (associated(c))
                if (.not. c%active_samp_group(samp_group)) then
                   c => c%nextComp()
                   cycle
                end if
                select type (c)
                class is (comm_diffuse_comp)
                   if (trim(c%cltype) /= 'none') then
                      allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
                      call cr_extract_comp(c%id, x_out, alm)
                      call c%Cl%sqrtS(alm=alm, info=c%x%info) ! Multiply with sqrt(Cl)  
                      call cr_insert_comp(c%id, .false., alm, x_out)
                      deallocate(alm)
                   end if
                class is (comm_ptsrc_comp)
                   if (c%myid == 0 .and. .not. c%precomputed_amps) then
                      call cr_extract_comp(c%id, x_out, pamp)
                      do j = 1, c%nmaps
                         do k = 1, c%nsrc
                            pamp(k-1,j) = pamp(k-1,j) * c%src(k)%P_x(j,2) ! Multiply with sqrtInvS
                         end do
                      end do
                      call cr_insert_comp(c%id, .false., pamp, x_out)
                      deallocate(pamp)
                   end if
                class is (comm_template_comp)
                   if (c%myid == 0) then
                      call cr_extract_comp(c%id, x_out, pamp)
                      pamp = pamp * c%P_cg(2) ! Multiply with sqrtInvS
                      call cr_insert_comp(c%id, .false., pamp, x_out)
                      deallocate(pamp)
                   end if
                end select
                c => c%nextComp()
             end do
             call cr_x2amp(samp_group, x_out)
             call output_FITS_sample(cpar, i, .false.)
             deallocate(x_out)
             call cr_x2amp(samp_group, x)
          end if
       end if
       call update_status(status, "cg5")

       !if (cpar%myid == root) write(*,*) x(size(x)-1:size(x))

       call wall_time(t2)
       if (cpar%myid_chain == root .and. cpar%verbosity > 2) then
          if (trim(cpar%cg_conv_crit) == 'residual' .or. trim(cpar%cg_conv_crit) == 'fixed_iter') then
!!$             write(*,*) '  CG iter. ', i, ' -- res = ', &
!!$                  & val_convergence, ', tol = ', lim_convergence, &
!!$                  & ', time = ', real(t2-t1,sp)
             buff = min(val_convergence,1d30)
             write(*,fmt='(a,i5,a,e13.5,a,e13.5,a,f8.2)') ' |  CG iter. ', i, ' -- res = ', &
                  & buff, ', tol = ', real(lim_convergence,sp), &
                  & ', time = ', real(t2-t1,sp)
          else if (trim(cpar%cg_conv_crit) == 'chisq') then
!             write(*,fmt='(a,i5,a,e13.5,a,f7.4,a,f8.2)') '  CG iter. ', i, ' -- chisq = ', &
!                  & real(chisq,sp), ', delta = ', real(val_convergence,sp), &
!                  & ', time = ', real(t2-t1,sp)
             write(*,*) '|  CG iter. ', i, ' -- chisq = ', &
                  & chisq, ', delta = ', val_convergence, &
                  & ', time = ', real(t2-t1,sp)
          end if
       end if

       call update_status(status, "cg6")

    end do

    ! Multiply with sqrt(S), and insert into right object
    c       => compList
    do while (associated(c))
       if (.not. c%active_samp_group(samp_group)) then
          c => c%nextComp()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, x, alm)
             call c%Cl%sqrtS(alm=alm, info=c%x%info) ! Multiply with sqrt(Cl)
             ! Add CMB dipole back again
             if (cpar%resamp_CMB .and. trim(c%type) == 'cmb' .and. .not. cpar%only_pol) &
                  & call add_fiducial_CMB_dipole(c%x%info, c%RJ2unit_(1), alm)
             call cr_insert_comp(c%id, .false., alm, x)
             deallocate(alm)
          end if
       class is (comm_ptsrc_comp)
          if (c%myid == 0 .and. .not. c%precomputed_amps) then
             call cr_extract_comp(c%id, x, pamp)
             do j = 1, c%nmaps
                do k = 1, c%nsrc
                   pamp(k-1,j) = pamp(k-1,j) * c%src(k)%P_x(j,2) ! Multiply with sqrtS
                end do
             end do
             call cr_insert_comp(c%id, .false., pamp, x)
             deallocate(pamp)
          end if
       class is (comm_template_comp)
          if (c%myid == 0) then
             call cr_extract_comp(c%id, x, pamp)
             pamp = pamp * c%P_cg(2) ! Multiply with sqrtS
             call cr_insert_comp(c%id, .false., pamp, x)
             deallocate(pamp)
          end if
       end select
       c => c%nextComp()
    end do
    !call update_status(status, "cr8")

    if (i >= maxiter .and. trim(cpar%cg_conv_crit) /= 'fixed_iter') then
       write(*,*) 'ERROR: Convergence in CG search not reached within maximum'
       write(*,*) '       number of iterations = ', maxiter
       stat = stat + 1
    else
       if (cpar%myid_chain == root .and. cpar%verbosity > 1) then
          write(*,fmt='(a,i5,a,e13.5,a,e13.5,a,f8.2)') ' |  Final CG iter ', i, ' -- res = ', &
               & real(val_convergence,sp), ', tol = ', real(lim_convergence,sp)
       end if
    end if

    deallocate(Ax, r, d, q, s)
    call update_status(status, "cr9")
    
  end subroutine solve_cr_eqn_by_CG

  subroutine cr_compute_chisq(comm, samp_group, x, chisq)
    implicit none
    integer(i4b),            intent(in)  :: comm, samp_group
    real(dp), dimension(1:), intent(in)  :: x
    real(dp),                intent(out) :: chisq

    integer(i4b) :: j, k, n
    real(dp), allocatable, dimension(:)   :: x_out
    real(dp), allocatable, dimension(:,:) :: alm, pamp
    class(comm_comp),   pointer :: c => null()

    n = size(x)

    ! Multiply with sqrt(S)     
    allocate(x_out(n))
    x_out = x
    c       => compList
    do while (associated(c))
       if (.not. c%active_samp_group(samp_group)) then
          c => c%nextComp()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, x_out, alm)
             call c%Cl%sqrtS(alm=alm, info=c%x%info) ! Multiply with sqrt(Cl)  
             call cr_insert_comp(c%id, .false., alm, x_out)
             deallocate(alm)
          end if
       class is (comm_ptsrc_comp)
          if (c%myid == 0 .and. .not. c%precomputed_amps) then
             call cr_extract_comp(c%id, x_out, pamp)
             do j = 1, c%nmaps
                do k = 1, c%nsrc
                   pamp(k-1,j) = pamp(k-1,j) * c%src(k)%P_x(j,2) ! Multiply with sqrtInvS
                end do
             end do
             call cr_insert_comp(c%id, .false., pamp, x_out)
             deallocate(pamp)
          end if
       class is (comm_template_comp)
          if (c%myid == 0) then
             call cr_extract_comp(c%id, x_out, pamp)
             pamp = pamp * c%P_cg(2) ! Multiply with sqrtInvS
             call cr_insert_comp(c%id, .false., pamp, x_out)
             deallocate(pamp)
          end if
       end select
       c => c%nextComp()
    end do
    call cr_x2amp(samp_group, x_out)
    call compute_chisq(comm, chisq_fullsky=chisq)
    deallocate(x_out)
    call cr_x2amp(samp_group, x)

  end subroutine cr_compute_chisq

  subroutine cr_amp2x_full(samp_group, x) 
    implicit none
    integer(i4b),           intent(in)  :: samp_group
    real(dp), dimension(:), intent(out) :: x

    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c => null()

    ! Stack parameters linearly
    ind = 1
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp) 
          do i = 1, c%x%info%nmaps
             if (c%active_samp_group(samp_group)) x(ind:ind+c%x%info%nalm-1) = c%x%alm(:,i)
             ind = ind + c%x%info%nalm
          end do
       class is (comm_ptsrc_comp)
          if (c%myid == 0 .and. .not. c%precomputed_amps) then
             do i = 1, c%nmaps
                if (c%active_samp_group(samp_group)) x(ind:ind+c%nsrc-1) = c%x(:,i)
                ind = ind + c%nsrc
             end do
          end if
       class is (comm_template_comp)
          if (c%myid == 0) then
             if (c%active_samp_group(samp_group)) x(ind) = c%x(1,1)
             ind    = ind + 1
          end if
       end select
       c => c%nextComp()
    end do

  end subroutine cr_amp2x_full

  subroutine cr_x2amp_full(samp_group, x)
    implicit none
    integer(i4b),           intent(in) :: samp_group
    real(dp), dimension(:), intent(in) :: x

    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c => null()

    ind = 1
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          do i = 1, c%x%info%nmaps
             if (c%active_samp_group(samp_group)) c%x%alm(:,i) = x(ind:ind+c%x%info%nalm-1)
             ind = ind + c%x%info%nalm
          end do
       class is (comm_ptsrc_comp)
          do i = 1, c%nmaps
            if (c%myid == 0) then
              if (c%active_samp_group(samp_group)) c%x(:,i) = x(ind:ind+c%nsrc-1)
              ind = ind + c%nsrc
            end if
          end do
       class is (comm_template_comp)
          if (c%myid == 0) then
             if (c%active_samp_group(samp_group)) c%x(1,1) = x(ind)
             ind      = ind + 1
          end if
       end select
       c => c%nextComp()
    end do
    
  end subroutine cr_x2amp_full

  ! ---------------------------
  ! Definition of linear system
  ! ---------------------------

  subroutine cr_computeRHS(operation, resamp_cmb, only_pol, handle, handle_noise, mask, samp_group, rhs)
    implicit none
    character(len=*),                            intent(in)             :: operation
    logical(lgt),                                intent(in)             :: resamp_cmb, only_pol
    type(planck_rng),                            intent(inout)          :: handle, handle_noise
    integer(i4b),                                intent(in)             :: samp_group
    real(dp),         allocatable, dimension(:), intent(in)             :: mask
    real(dp),         allocatable, dimension(:), intent(out)            :: rhs

    integer(i4b) :: i, j, l, m, k, n, ierr
    real(dp)     :: tmp
    class(comm_map),     pointer                 :: map  => null()
    class(comm_map),     pointer                 :: Tm   => null()
    class(comm_map),     pointer                 :: mu   => null()
    class(comm_comp),    pointer                 :: c    => null()
    class(comm_mapinfo), pointer                 :: info => null()
    real(dp),        allocatable, dimension(:,:) :: eta, Tp

    ! Initialize output vector
    allocate(rhs(ncr))
    rhs = 0.d0

    ! Add channel dependent terms
    do i = 1, numband

       if (.not. data(i)%cr_active) cycle
       
       ! Set up Wiener filter term
       map => compute_residual(i, cg_samp_group=samp_group) 

!!$       if (map%info%myid == 0) write(*,*) sum(abs(map%map)), 'pre sqrtinvN'
!!$       call data(i)%N%sqrtInvN(map, samp_group=samp_group)
!!$       if (map%info%myid == 0) write(*,*) sum(abs(map%map)), 'post sqrtinvN'
!!$       call data(i)%N%sqrtInvN(map, samp_group=samp_group)
!!$       if (map%info%myid == 0) write(*,*) sum(abs(map%map)), 'post sqrtinvN2'
!!$       call map%dealloc()
!!$
!!$       map => compute_residual(i, cg_samp_group=samp_group) 
!!$
!!$       if (map%info%myid == 0) write(*,*) sum(abs(map%map)), 'pre sqrtinvN'
!!$       call data(i)%N%invN(map, samp_group=samp_group)
!!$       if (map%info%myid == 0) write(*,*) sum(abs(map%map)), 'post invN'
!!$
!!$       call data(i)%N%N(map, samp_group=samp_group)
!!$       if (map%info%myid == 0) write(*,*) sum(abs(map%map)), 'post N (should b back to above)'
!!$
!!$
!!$
!!$       call mpi_finalize(ierr)
!!$       stop

       ! Subtract CMB dipole if resamp mode, to avoid large condition numbers; add back later
       if (resamp_cmb .and. .not. only_pol) call subtract_fiducial_CMB_dipole(i, map)

       ! Apply projection matrix, ie., mask in pixel space and multipoles above lmax in harmonic space
!!$       map%map = map%map * data(i)%mask%map
!!$       call map%YtW
!!$       call map%Y

       ! Add channel-dependent white noise fluctuation
       if (trim(operation) == 'sample') then
          call data(i)%N%sqrtInvN(map, samp_group=samp_group)           ! Multiply with sqrt(invN)
          do k = 1, map%info%nmaps
             do j = 0, map%info%np-1
                map%map(j,k) = map%map(j,k) + rand_gauss(handle)
                !tmp          = rand_gauss(handle)
                !map%map(j,k) = map%map(j,k) + rand_gauss(handle_noise)
             end do
          end do
          call data(i)%N%sqrtInvN(map, samp_group=samp_group)          ! Multiply with sqrt(invN)
       else
          call data(i)%N%invN(map, samp_group=samp_group)          ! Multiply with sqrt(invN)
       end if

       ! Convolve with transpose beam
       call map%Yt()
       call data(i)%B(0)%p%conv(trans=.true., map=map)

       ! Multiply with (transpose and component specific) mixing matrix, and
       ! insert into correct segment
       c => compList
       do while (associated(c))
          if (.not. c%active_samp_group(samp_group)) then
             c => c%nextComp()
             cycle
          end if
          select type (c)
          class is (comm_diffuse_comp)
             info  => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, c%lmax_amp, &
                  & c%nmaps, c%nmaps==3)
             Tm     => comm_map(info)
             if (c%F_null(i,0)) then
                Tm%alm = 0.d0
             else
                call map%alm_equal(Tm)
                !need to check all relevant polarizations for lmax_ind == 0
                if (all(c%lmax_ind_mix(1:min(c%nmaps,data(i)%info%nmaps),:) == 0)) then 
                   do j = 1, c%nmaps
                      Tm%alm(:,j) = Tm%alm(:,j) * c%F_mean(i,0,j)
                   end do
                else
                   call Tm%Y()
                   do j = 1, c%nmaps
                      if (j <= data(i)%info%nmaps) then
                         Tm%map(:,j) = c%F(i,0)%p%map(:,j) * Tm%map(:,j)
                      else
                         Tm%map(:,j) = 0.d0
                      end if
                   end do
                   call Tm%YtW()
                end if
             end if
             call c%Cl%sqrtS(map=Tm) ! Multiply with sqrt(Cl)

             ! Apply projection operator
             do l = 0, Tm%info%nalm-1
                if (Tm%info%lm(1,l) > data(i)%info%lmax) Tm%alm(l,:) = 0.d0
             end do

             call cr_insert_comp(c%id, .true., Tm%alm, rhs)
             call Tm%dealloc(); deallocate(Tm)
             nullify(info)
          class is (comm_ptsrc_comp)
             if(.not. c%precomputed_amps) then 
                 
               allocate(Tp(c%nsrc,c%nmaps))
               Tp = c%projectBand(i,map)
               if (c%myid == 0) then
                  do j = 1, c%nmaps
                     do k = 1, c%nsrc
                        Tp(k,j) = Tp(k,j) * c%src(k)%P_x(j,2)
                     end do
                  end do
                  call cr_insert_comp(c%id, .true., Tp, rhs)
               end if
               deallocate(Tp)
             end if
          class is (comm_template_comp)
             allocate(Tp(1,1))
             Tp = c%projectBand(i,map)
             if (c%myid == 0) then
                Tp = Tp * c%P_cg(2)
                call cr_insert_comp(c%id, .true., Tp, rhs)
             end if
             deallocate(Tp)
          end select
          c => c%nextComp()
       end do

       call map%dealloc(); deallocate(map)
    end do

    ! Add prior terms
    c => compList
    do while (associated(c))
       if (.not. c%active_samp_group(samp_group)) then
          c => c%nextComp()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             n = ind_comp(c%id,2)
             allocate(eta(0:c%x%info%nalm-1,c%x%info%nmaps))
             eta = 0.d0
             ! Variance term
             if (trim(operation) == 'sample') then
                do j = 1, c%x%info%nmaps
                   if (j == 1 .and. only_pol) cycle
                   do i = 0, c%x%info%nalm-1
                      eta(i,j) = rand_gauss(handle)
                   end do
                end do
             end if
             ! Mean term
             if (associated(c%mu)) then
                mu => comm_map(c%mu)
                call c%Cl%sqrtInvS(map=mu)
                do j = 1, c%x%info%nmaps
                   do i = 0, c%x%info%nalm-1
                      !if (mu%info%lm(1,i) <= c%lmax_prior) then
                         !write(*,*) j, i, mu%info%lm(i,1), c%lmax_prior
                         eta(i,j) = eta(i,j) + mu%alm(i,j)
                      !end if
                   end do
                end do
                !eta = eta + mu%alm
                call mu%dealloc(); deallocate(mu)
             end if
             call cr_insert_comp(c%id, .true., eta, rhs)
             deallocate(eta)
          end if
       class is (comm_ptsrc_comp)
          if (.not. c%precomputed_amps) then 
            if (c%myid == 0) then
               allocate(eta(1:c%nsrc,c%nmaps))
               eta = 0.d0
               ! Variance term
               if (trim(operation) == 'sample') then
                  do j = 1, c%nmaps
                     do i = 1, c%nsrc
                        eta(i,j) = rand_gauss(handle)
                     end do
                  end do
               end if
               ! Mean term
               do j = 1, c%nmaps
                  do i = 1, c%nsrc
                     eta(i,j) = eta(i,j) + c%src(i)%P_x(j,1)/c%src(i)%P_x(j,2)
                  end do
               end do
               call cr_insert_comp(c%id, .true., eta, rhs)
               deallocate(eta)
            end if
          end if
       class is (comm_template_comp)
          if (c%myid == 0) then
             allocate(eta(1,1))
             eta = 0.d0
             ! Variance term
             if (trim(operation) == 'sample') then
                eta(1,1) = rand_gauss(handle)
             end if
             ! Mean term
             eta(1,1) = eta(1,1) + c%P_cg(1)/c%P_cg(2)
             call cr_insert_comp(c%id, .true., eta, rhs)
             deallocate(eta)
          end if
       end select
       c => c%nextComp()

    end do
    nullify(c)

  end subroutine cr_computeRHS

  function cr_matmulA(x, samp_group)
    implicit none

    real(dp),     dimension(1:),     intent(in)  :: x
    integer(i4b),                    intent(in)  :: samp_group
    real(dp),     dimension(size(x))             :: cr_matmulA

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, l, myid, lmax
    class(comm_map),  pointer :: map => null(), pmap => null(), map_buff => null()
    class(comm_comp), pointer :: c => null()
    real(dp),        allocatable, dimension(:)   :: y, sqrtS_x
    real(dp),        allocatable, dimension(:,:) :: alm, m, pamp
    class(comm_mapinfo), pointer :: info => null()

!    write(*,*) 'df1'
    ! Initialize output array
    !call update_status(status, "A1")
    allocate(y(ncr), sqrtS_x(ncr))
    y = 0.d0
    myid = data(1)%map%info%myid

    ! Multiply with sqrt(S)
    call wall_time(t1)
    !call update_status(status, "A2")
    sqrtS_x = x
    !call update_status(status, "A3")
    c       => compList
    lmax    = -1
!    !write(*,*) 'df2'
    do while (associated(c))
       if (.not. c%active_samp_group(samp_group)) then
          c => c%nextComp()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%cltype) /= 'none') then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, sqrtS_x, alm)
             !call update_status(status, "A4")
             call c%Cl%sqrtS(alm=alm, info=c%x%info) ! Multiply with sqrt(Cl)
             !call update_status(status, "A5")
             call cr_insert_comp(c%id, .false., alm, sqrtS_x)
             deallocate(alm)
          end if
          lmax = max(max(lmax, c%lmax_amp),2)
       class is (comm_ptsrc_comp)
          if (c%myid == 0 .and. .not. c%precomputed_amps) then
             call cr_extract_comp(c%id, sqrtS_x, pamp)
             do j = 1, c%nmaps
                do i = 1, c%nsrc
                   pamp(i-1,j) = pamp(i-1,j) * c%src(i)%P_x(j,2) ! Multiply with sqrtS
                end do
             end do
             call cr_insert_comp(c%id, .false., pamp, sqrtS_x)
             deallocate(pamp)
          end if
       class is (comm_template_comp)
          if (c%myid == 0) then
             call cr_extract_comp(c%id, sqrtS_x, pamp)
             pamp = pamp * c%P_cg(2) ! Multiply with sqrtS
             call cr_insert_comp(c%id, .false., pamp, sqrtS_x)
             deallocate(pamp)
          end if
       end select
       c => c%nextComp()
    end do
    call wall_time(t2)
    !if (myid == 0) write(*,fmt='(a,f8.2)') 'sqrtS time = ', real(t2-t1,sp)
    !write(*,*) 'df3' 

    
    ! Add frequency dependent terms
    do i = 1, numband

       if (.not. data(i)%cr_active) cycle
       
       ! Compute component-summed map, ie., column-wise matrix elements
       call wall_time(t1)
       map  => comm_map(data(i)%info)   ! For diffuse components
       pmap => comm_map(data(i)%info)   ! For point-source components and alm-buffer for diffuse components
       c   => compList
    !write(*,*) 'df41'
       do while (associated(c))
          if (.not. c%active_samp_group(samp_group)) then
             c => c%nextComp()
             cycle
          end if
          select type (c)
          class is (comm_diffuse_comp)
    !write(*,*) 'df42'
             call cr_extract_comp(c%id, sqrtS_x, alm)
             do l = 0, c%x%info%nalm-1
                if (c%x%info%lm(1,l) > data(i)%info%lmax) alm(l,:) = 0.d0
             end do
    !write(*,*) 'df43'
             call pmap%set_alm(alm,c%x%info)
    !write(*,*) 'df44'
             allocate(m(0:data(i)%info%nalm-1,data(i)%info%nmaps))
             !allocate(m(0:c%x%info%nalm-1,c%x%info%nmaps))
             !call update_status(status, "A6")
             m = c%getBand(i, amp_in=pmap%alm, alm_out=.true.)
    !write(*,*) 'df45'
             !call update_status(status, "A7")
             map%alm = map%alm + m
             !call map%add_alm(m, c%x%info)
             !call update_status(status, "A8")
             deallocate(alm, m)
    !write(*,*) 'df6'
          class is (comm_ptsrc_comp)
             if(.not. c%precomputed_amps) then
               call cr_extract_comp(c%id, sqrtS_x, pamp)
               allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
               m = c%getBand(i, amp_in=pamp)
               pmap%map = pmap%map + m
               deallocate(pamp, m)
             end if
          class is (comm_template_comp)
             call cr_extract_comp(c%id, sqrtS_x, pamp)
             allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
             m = c%getBand(i, amp_in=pamp)
             pmap%map = pmap%map + m
             deallocate(pamp, m)
          end select
          c => c%nextComp()
       end do
       !call update_status(status, "A9")
    !write(*,*) 'df49'
       if (lmax > -1) then
          info     => comm_mapinfo(map%info%comm, map%info%nside, lmax, map%info%nmaps, map%info%nmaps==3)
          map_buff => comm_map(info)
    !write(*,*) 'df410'
          call map%alm_equal(map_buff)
          call map_buff%Y()                    ! Diffuse components
          map%map = map_buff%map
    !write(*,*) 'df411'
       else
    !write(*,*) 'df412'
          map%map = 0.d0
       end if
    !write(*,*) 'df413'
       !call update_status(status, "A10")
       map%map = map%map + pmap%map    ! Add compact objects
       !call update_status(status, "A11")
       !write(*,*) 'c', sum(abs(pmap%map))
       call wall_time(t2)
       !if (myid == 0) !write(*,fmt='(a,f8.2)') 'getBand time = ', real(t2-t1,sp)

    !write(*,*) 'df5'
       ! Multiply with invN
       call wall_time(t1)
       call data(i)%N%InvN(map, samp_group=samp_group)
       call wall_time(t2)
       !call update_status(status, "A12")
       !if (myid == 0) write(*,fmt='(a,f8.2)') 'invN time = ', real(t2-t1,sp)

       ! Project summed map into components, ie., row-wise matrix elements
       call wall_time(t1)
       c   => compList
       if (lmax > -1) then
          map_buff%map = map%map
          call map_buff%Yt()             ! Prepare for diffuse components
          call map_buff%alm_equal(map)
          call map_buff%dealloc(); deallocate(map_buff)
       end if
       !call update_status(status, "A13")
       do while (associated(c))
          if (.not. c%active_samp_group(samp_group)) then
             c => c%nextComp()
             cycle
          end if
          select type (c)
          class is (comm_diffuse_comp)
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             !call update_status(status, "A14")
             alm = c%projectBand(i, map, alm_in=.true.)
             !call update_status(status, "A15")
             do l = 0, c%x%info%nalm-1
                if (c%x%info%lm(1,l) > data(i)%info%lmax) alm(l,:) = 0.d0
             end do
             call cr_insert_comp(c%id, .true., alm, y)
             deallocate(alm)
          class is (comm_ptsrc_comp)
             if(.not. c%precomputed_amps) then
               allocate(pamp(0:c%nsrc-1,c%nmaps))
               pamp = c%projectBand(i, map)
               if (c%myid == 0) call cr_insert_comp(c%id, .true., pamp, y)
               deallocate(pamp)
             end if
          class is (comm_template_comp)
             allocate(pamp(1,1))
             pamp = c%projectBand(i, map)
             if (c%myid == 0) call cr_insert_comp(c%id, .true., pamp, y)
             deallocate(pamp)
          end select
          c => c%nextComp()
       end do
       call wall_time(t2)
       !if (myid == 0) write(*,fmt='(a,f8.2)') 'projBand time = ', real(t2-t1,sp)

       call map%dealloc(); deallocate(map)
       call pmap%dealloc(); deallocate(pmap)
    end do
    !call update_status(status, "A16")
    !write(*,*) 'df6'

    ! Add prior term and multiply with sqrt(S) for relevant components
    call wall_time(t1)
    c   => compList
    do while (associated(c))
       if (.not. c%active_samp_group(samp_group)) then
          c => c%nextComp()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
    !write(*,*) 'df8'
          if (trim(c%cltype) /= 'none') then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             ! Multiply with sqrt(Cl)
             call cr_extract_comp(c%id, y, alm)
             !call update_status(status, "A17")
             call c%Cl%sqrtS(alm=alm, info=c%x%info)
             !call update_status(status, "A18")
             call cr_insert_comp(c%id, .false., alm, y)
             ! Add (unity) prior term
             call cr_extract_comp(c%id, x, alm)
             call cr_insert_comp(c%id, .true., alm, y)
             deallocate(alm)
          end if
    !write(*,*) 'df9'
       class is (comm_ptsrc_comp)
          if (c%myid == 0 .and. .not. c%precomputed_amps) then
             ! Multiply with sqrt(Cl)
             call cr_extract_comp(c%id, y, pamp)
             do j = 1, c%nmaps
                do i = 1, c%nsrc
                   pamp(i-1,j) = pamp(i-1,j) * c%src(i)%P_x(j,2) ! Multiply with sqrtS
                end do
             end do
             call cr_insert_comp(c%id, .false., pamp, y)
             ! Add (unity) prior term
             call cr_extract_comp(c%id, x, pamp)
             call cr_insert_comp(c%id, .true., pamp, y)
             deallocate(pamp)
          end if
       class is (comm_template_comp)
    !write(*,*) 'df7'
          if (c%myid == 0) then
             ! Multiply with sqrt(Cl)
             call cr_extract_comp(c%id, y, pamp)
             pamp = pamp * c%P_cg(2) ! Multiply with sqrtS
             call cr_insert_comp(c%id, .false., pamp, y)
             ! Add (unity) prior term
             call cr_extract_comp(c%id, x, pamp)
             call cr_insert_comp(c%id, .true., pamp, y)
             deallocate(pamp)
          end if
       end select
       c => c%nextComp()
    end do
    nullify(c)
    call wall_time(t2)
    !write(*,*) 'df10'
    !if (myid == 0) write(*,fmt='(a,f8.2)') 'prior time = ', real(t2-t1,sp)
    !call mpi_finalize(i)
    !stop

    ! Return result and clean up
    call wall_time(t1)
    !call update_status(status, "A19")
    cr_matmulA = y
    !call update_status(status, "A20")
    deallocate(y, sqrtS_x)
    call wall_time(t2)
    !    write(*,*) 'f', t2-t1
        !write(*,*) 'df11'
  end function cr_matmulA

  function cr_invM(comm, x, samp_group)
    implicit none
    integer(i4b),                        intent(in) :: comm, samp_group
    real(dp),              dimension(:), intent(in) :: x
    real(dp), allocatable, dimension(:)             :: cr_invM

    integer(i4b) :: ierr
    logical(lgt) :: Q_is_active
    real(dp), allocatable, dimension(:,:) :: alm, alm0
    real(dp), allocatable, dimension(:)   :: Qx
    class(comm_comp), pointer :: c => null()

    if (.not. allocated(cr_invM)) allocate(cr_invM(size(x)))
    allocate(Qx(size(x)))
    cr_invM = x

    call applyDeflatePrecond(cr_invM, Qx)
    Q_is_active = .false.  !any(Qx /= 0.d0)
    !call mpi_allreduce(MPI_IN_PLACE, Q_is_active, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
    if (Q_is_active) cr_invM = x - cr_matmulA(Qx, samp_group)   

!!$    if (size(x) > 0) write(*,*) sum(abs(x)), sum(abs(cr_invM)), sum(abs(Qx))
!!$    call mpi_finalize(ierr)
!!$    stop

    call applyDiffPrecond(cr_invM)
    call applyPtsrcPrecond(cr_invM)
    call applyTemplatePrecond(cr_invM)

    if (Q_is_active) cr_invM = cr_invM + Qx
    
    ! Apply low-l preconditioner
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (c%lmax_pre_lowl > -1) then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             allocate(alm0(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, x, alm)
             call cr_extract_comp(c%id, cr_invM, alm0)
             call c%applyLowlPrecond(alm, alm0)
             call cr_insert_comp(c%id, .false., alm, cr_invM)
             deallocate(alm, alm0)
          end if
       end select
       c => c%nextComp()
    end do

    deallocate(Qx)

  end function cr_invM


  subroutine applyDeflatePrecond(x, Qx)
    real(dp), dimension(:), intent(in)  :: x
    real(dp), dimension(:), intent(out) :: Qx

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:,:) :: alm, Qalm
    class(comm_comp), pointer :: c => null()

    Qx = 0.d0

    return

    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (c%lmax_def > -1) then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             allocate(Qalm(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, x, alm)
             call c%applyDeflatePrecond(alm, Qalm)
             call cr_insert_comp(c%id, .false., Qalm, Qx)
             deallocate(alm, Qalm)
          end if
       end select
       c => c%nextComp()
    end do

  end subroutine applyDeflatePrecond


  subroutine update_precond(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update
    class(comm_comp), pointer :: c => null()
    logical(lgt), save :: first_call = .true.

    ! Set up deflation preconditioner for CMB+diagonal only
    if (.false.) then
       c   => compList
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             call c%updateDeflatePrecond()
          end select
          c => c%nextComp()
       end do
    end if

    call updateDiffPrecond(samp_group, force_update)
    call updatePtsrcPrecond(samp_group)
    call updateTemplatePrecond(samp_group)

    !if (.not. first_call) return

    ! Set up low-l preconditioner for CMB+diagonal only
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (c%lmax_pre_lowl > -1) then
             call c%updateLowlPrecond()
          end if
       end select
       c => c%nextComp()
    end do

    first_call = .false.

  end subroutine update_precond
  
end module comm_cr_mod
