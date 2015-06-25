module comm_cr_mod
  use comm_comp_mod
  use comm_data_mod
  use comm_param_mod
  implicit none

contains


  recursive subroutine solve_cr_eqn_by_CG(cpar, A, invM, x, b, stat, Nscale)
    implicit none

    type(comm_params),                intent(in)             :: cpar
    real(dp),          dimension(1:), intent(out)            :: x
    real(dp),          dimension(1:), intent(in)             :: b
    integer(i4b),                     intent(out)            :: stat
    real(dp),                         intent(in),   optional :: Nscale

    interface
       recursive function A(x, Nscale)
         use healpix_types
         implicit none
         real(dp), dimension(:),       intent(in)           :: x
         real(dp), dimension(size(x))                       :: A
         real(dp),                     intent(in), optional :: Nscale
       end function A

       recursive function invM(x, Nscale)
         use healpix_types
         implicit none
         real(dp), dimension(:),      intent(in)           :: x
         real(dp), dimension(size(x))                      :: invM
         real(dp),                    intent(in), optional :: Nscale
       end function invM
    end interface

    integer(i4b) :: i, j, k, l, m, n, maxiter, root
    real(dp)     :: eps, tol, delta0, delta_new, delta_old, alpha, beta, t1, t2
    real(dp), allocatable, dimension(:) :: Ax, r, d, q, invM_r, temp_vec, s

    root    = 0
    maxiter = cpar%cg_maxiter
    eps     = cpar%cg_tol
    n       = size(x)

    ! Allocate temporary data vectors
    allocate(Ax(n), r(n), d(n), q(n), invM_r(n), s(n))

    ! Initialize the CG search
    x  = 0.d0
    r  = b-A(x)
    d  = invM(r)

    delta_new = dot_product(r,d)
    delta0    = delta_new
    do i = 1, maxiter
       call wall_time(t1)
       
       if (delta_new < eps**2 * delta0) exit

       q     = A(d)
       alpha = delta_new / dot_product(d, q)
       x     = x + alpha * d

       ! Restart every 50th iteration to suppress numerical errors
       if (mod(i,50) == 0) then
          r = b - A(x)
       else
          r = r - alpha*q
       end if

       s         = invM(r)
       delta_old = delta_new 
       delta_new = dot_product(r, s)
       beta      = delta_new / delta_old
       d         = s + beta * d

       call wall_time(t2)
       if (cpar%myid == root .and. cpar%verbosity > 2) then
          write(*,fmt='(a,i5,a,e8.3,a,e8.3,a,f8.3)') 'CG iter. ', i, ' -- res = ', &
               & real(delta_new,sp), ', tol = ', real(eps**2 * delta0,sp), &
               & ', wall time = ', real(t2-t1,sp)
       end if

    end do

    if (i >= maxiter) then
       write(*,*) 'ERROR: Convergence in CG search not reached within maximum'
       write(*,*) '       number of iterations = ', maxiter
       stat = stat + 1
    else
       if (cpar%myid == root .and. cpar%verbosity > 1) then
          write(*,fmt='(a,i5,a,e8.3,a,e8.3,a,f8.3)') 'Final CG iter ', i-1, ' -- res = ', &
               & real(delta_new,sp), ', tol = ', real(eps**2 * delta0,sp)
       end if
       
!!$       ! Output results and change to physical variables
!!$       call genvec_set_equal(x, x_out)
!!$       call multiply_by_sqrt_S(.false., x_out%cmb_amp)
    end if

    deallocate(Ax, r, d, q, invM_r, s)
    
  end subroutine solve_cr_eqn_by_CG

  function cr_amp2x()
    implicit none

    real(dp), allocatable, dimension(:) :: cr_amp2x

    integer(i4b) :: n
    class(comm_comp), pointer :: c

    ! Find number of free CR parameters for current processor
    n = 0
    c => compList
    do while (associated(c))
       n = n + c%ncr
       c => c%next()
    end do

    write(*,*) 'n = ', n
    allocate(cr_amp2x(n))
    cr_amp2x = 1.d0

  end function cr_amp2x
  
end module comm_cr_mod
