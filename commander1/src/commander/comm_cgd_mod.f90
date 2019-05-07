module comm_cgd_mod
  use comm_mp_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  real(dp),         private :: rel_tol
  integer(i4b),     private :: verbosity, maxiter
  integer(i4b),     private :: chain, root = 0, ierr
  logical(lgt),     private :: enforce_zero_cl

contains


  ! Routine for solving the big, difficult linear system,   Au = v,  Ax = b
  subroutine draw_constrained_realization_by_CG(precond_type, b, x_out, stat, init)
    implicit none
    
    integer(i4b), intent(in)    :: precond_type
    type(genvec), intent(in)    :: b
    type(genvec), intent(inout) :: x_out
    integer(i4b), intent(inout) :: stat
    type(genvec), intent(in), optional :: init

    real(dp)     :: t1, t2, t3, t4, tol
    integer(i4b) :: i, j, k, l, m
    real(dp)     :: alpha, beta, r_invM_r_i, r_invM_r_ip1, d_A_d, cl, delta_new, delta_0, delta_old, r_max

    type(genvec) :: Ax, x, r, d, invM_r, temp_vec, q, s, x_old

    character(len=4) :: iter_text
    real(dp), allocatable, dimension(:,:,:) :: cls0, cls1
    real(dp), allocatable, dimension(:)     :: x_lin


    ! Allocate data vectors
    call allocate_genvec(Ax)
    call allocate_genvec(x)
    call allocate_genvec(r)
    call allocate_genvec(d)
    call allocate_genvec(q)
    call allocate_genvec(invM_r)
    call allocate_genvec(temp_vec)
    call allocate_genvec(s)
    if (enforce_zero_cl) call allocate_genvec(x_old)
    allocate(cls0(0:lmax,nmaps,nmaps), cls1(0:lmax,nmaps,nmaps))

!!$       call genvec_set_equal(b, x)
!!$       x%cmb_amp(1:4,1) = 0.d0
!!$       do l = 2, 3
!!$          do m = -l, l
!!$             i = lm2ind(l,m)
!!$             write(*,*) x%cmb_amp(i,1)
!!$          end do
!!$       end do
!!$       call compute_Au(x, r)
!!$       write(*,*)
!!$       do l = 2, 3
!!$          do m = -l, l
!!$             i = lm2ind(l,m)
!!$             write(*,*) r%cmb_amp(i,1)
!!$          end do
!!$       end do
!!$       call mpi_finalize(ierr)
!!$       stop

    
    if (present(init)) then
       call genvec_set_equal(init, x)
       call multiply_by_inv_sqrt_S(.false., x%cmb_amp)
    end if

    ! Set up convergence criterion
!    tol    = sqrt(abs(genvec_dot_product(v,v))) * rel_tol

    ! Initialize preconditioner
    if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- initializing preconditioner'
    call wall_time(t1)
    call init_precond(precond_type)
    call wall_time(t2)
    if (verbosity > 2) then
       write(*,*) 'Chain no. ', chain, ' -- init precond time = ', real(t2-t1,sp)
    end if


    ! Do the Conjugate Gradient search -- see Shewchuck 1994 ("Conjugate
    ! Gradients without the Agonizing Pain"), page 46
    if (verbosity > 1) write(*,*) 'Chain no. ', chain, ' -- Starting CG search'
!    open(98,file='cls.dat')

!!$    write(*,*) 'precond_type = ', precond_type
!!$    write(*,*) 'mask_state   = ', mask_state == OUTSIDE_MASK
!!$    do i = 1, numcomp
!!$       !call nullify_genvec(x)
!!$       !call nullify_genvec(Ax)
!!$       x%cmb_amp      = 0.d0
!!$       x%cmb_amp(i,3) = 1.d0
!!$!    write(*,*) 'a', x%cmb_amp(6,:)
!!$!    write(*,*) 'a', sum(abs(x%cmb_amp(6,:)))
!!$       call compute_Au(x, Ax)
!!$!    write(*,*) 'b', Ax%cmb_amp(6,:)
!!$!    write(*,*) 'b', sum(abs(Ax%cmb_amp(6,:)))
!!$       call compute_invM_u(Ax, x)
!!$       write(*,*) i, x%cmb_amp(i,3), sum(abs(x%cmb_amp)), sum(abs(x%fg_amp)), sum(abs(x%temp_amp))
!!$    end do
!!$    call mpi_finalize(ierr)
!!$    stop


!    call genvec_set_equal(x, x_out)
!    call multiply_by_sqrt_S(.false., x_out%cmb_amp)
!    call output_maps_from_iteration(256, 750, 0, 0, 'chains_cls', x_out%cmb_amp)
!    call int2string(i, iter_text)
    !call write_map('chains_cls/dust'//iter_text//'.fits', x%fg_amp(:,:,1))
!    k = 0
!    do l = 0, 750
!       cl = 0.d0
!       do m = -l, l
!          cl = cl + x_out%cmb_amp(k,1)**2
!          k = k+1
!       end do
!       write(98,*) l, cl/real(2*l+1,dp) * l*(l+1)/2./pi
!    end do
!    write(98,*)

    ! Initialize CG
    call compute_Au(x, Ax)
    call genvec_minus(b, Ax, r)
    call compute_invM_u(r, d)

    delta_new = genvec_dot_product(r,d)
    delta_0   = delta_new
    i         = 0
    r_max     = 1.d30    
    do while (i < 2 .or. (i < maxiter .and. (r_max > rel_tol)))

       call wall_time(t1)
       call compute_Au(d, q)
       call wall_time(t2)
!       write(*,*) 'wall Au = ', t2-t1

       alpha = delta_new / genvec_dot_product(d, q)

       call genvec_plus_times(x, alpha, d, temp_vec)
       call genvec_set_equal(temp_vec, x)

       if (mod(i,50) == 0) then
          call compute_Au(x, Ax)
          call genvec_minus(b, Ax, r)
       else
          call genvec_plus_times(r, -alpha, q, temp_vec)
          call genvec_set_equal(temp_vec, r)
       end if

       call compute_invM_u(r, s)
       delta_old = delta_new
       delta_new = genvec_dot_product(r, s)
       beta      = delta_new / delta_old
       call genvec_plus_times(s, beta, d, d)       
       i = i+1

       ! Output results and change to physical variables
       if (.false.) then
          call genvec_set_equal(x, x_out)
          call multiply_by_sqrt_S(.false., x_out%cmb_amp)
          call output_maps_from_iteration(32, 64, 0, i, 'chains', x_out%cmb_amp, 420.d0)
          !call int2string(i, iter_text)
          !call write_map('chains/dust'//iter_text//'.fits', x%fg_amp(:,:,1))
!!$          k = 0
!!$          do l = 0, 750
!!$             cl = 0.d0
!!$             do m = -l, l
!!$                cl = cl + x_out%cmb_amp(k,1)**2
!!$                k = k+1
!!$             end do
!!$             write(98,*) l, cl/real(2*l+1,dp) * l*(l+1)/2./pi
!!$          end do
!!$          write(98,*)
       end if

       ! Check convergence 
       if (i == 1) then
          if (enforce_zero_cl) then
             call genvec_set_equal(x, x_old)
          else
             call compute_sigma_l(x%cmb_amp, cls0)
          end if          
       else
          if (enforce_zero_cl) then
             r_max = 0.d0
             do k = 1, size(x%fg_amp,3)
                do j = 1, size(x%fg_amp,2)
                   do l = 0, size(x%fg_amp,1)-1
                      if (x_old%fg_amp(l,j,k) /= 0.d0) then
                         r_max = max(r_max, abs((x_old%fg_amp(l,j,k) - x%fg_amp(l,j,k)) / x_old%fg_amp(l,j,k)))
                      end if
                   end do
                end do
             end do
             call genvec_set_equal(x, x_old)
          else
             call compute_sigma_l(x%cmb_amp, cls1)
             r_max = 0.d0
             do k = 1, nmaps
                do j = 1, nmaps
                   do l = 0, lmax
                      if (cls1(l,j,k) /= 0) then
                         r_max = max(r_max, abs((cls1(l,j,k)-cls0(l,j,k)) / cls1(l,j,k)))
!                         write(*,*) real(cls1(l,j,k),sp), real(cls0(l,j,k),sp), real(abs((cls1(l,j,k)-cls0(l,j,k)) / cls1(l,j,k)),sp)
                      end if
                   end do
                end do
             end do
             cls0  = cls1
          end if
       end if

       call wall_time(t2)
       if (verbosity > 2) then
          write(*,*) 'Chain no. ', chain, ' -- res = ', &
               & real(r_max,sp), ', tol = ',  real(rel_tol,sp)
          write(*,*) 'Chain no. ', chain, ' -- iter. ', i, ', wall time = ', real(t2-t1,sp)
       end if

    end do
!    close(98)

    if (i == maxiter) then
       write(*,*) 'ERROR: Convergence in CG search not reached within maximum'
       write(*,*) '       number of iterations = ', maxiter
       stat = stat + 1
    else

       if (verbosity > 1) then
          write(*,*) 'Chain no. ', chain, ' -- final res = ', &
               & real(r_max,sp), ', tol = ',  real(rel_tol,sp)
          write(*,*) 'Chain no. ', chain,' -- number of CGD iterations  = ', i
       end if
       
       ! Output results and change to physical variables
       call genvec_set_equal(x, x_out)
       call multiply_by_sqrt_S(.false., x_out%cmb_amp)
       
       ! Apply constraints on free template amplitudes
       !call apply_P_trans_constraint(x_out)

    end if

    call deallocate_genvec(Ax)
    call deallocate_genvec(x)
    call deallocate_genvec(r)
    call deallocate_genvec(d)
    call deallocate_genvec(invM_r)
    call deallocate_genvec(temp_vec)
    call deallocate_genvec(q)
    call deallocate_genvec(s)
    if (enforce_zero_cl) call deallocate_genvec(x_old)
    deallocate(cls0, cls1)
    
  end subroutine draw_constrained_realization_by_CG


  subroutine initialize_cgd_mod(chain_in, comm_chain, comm_data, comm_alms, &
       & paramfile)
    implicit none
    
    integer(i4b),                      intent(in)           :: chain_in
    integer(i4b),                      intent(in)           :: comm_chain, comm_data, comm_alms
    character(len=128),                intent(in)           :: paramfile

    logical(lgt) :: polarization

    chain = chain_in
    call get_parameter(paramfile, 'VERBOSITY',     par_int=verbosity)
    call get_parameter(paramfile, 'CGD_MAXITER',   par_int=maxiter)
    call get_parameter(paramfile, 'CGD_TOLERANCE', par_dp=rel_tol)
    call get_parameter(paramfile, 'ENFORCE_ZERO_CL', par_lgt=enforce_zero_cl)
    call get_parameter(paramfile, 'LMAX',            par_int=lmax)
    call get_parameter(paramfile, 'POLARIZATION',    par_lgt=polarization)
    if (polarization) then
       nmaps = 3
    else
       nmaps = 1
    end if
    
  end subroutine initialize_cgd_mod


end module comm_cgd_mod
 
