module powell_mod
  use healpix_types
  implicit none

  real(dp), parameter, private :: upper_bound = 1d30

contains

  function f1lin(x, p, xi, func)
    implicit none

    real(dp)   :: x
    real(dp), dimension(:)  :: p, xi
    real(dp)                :: f1lin
    interface
       function func(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                                     :: func
       end function func
    end interface

    real(dp), allocatable, dimension(:)  :: xt

    allocate(xt(0:size(p)-1))
    
    xt = p + x * xi

    f1lin = func(xt)
    
    deallocate(xt)

  end function f1lin


  subroutine powell(p, func, ierr, niter, tolerance)
    implicit none
    real(dp),     dimension(0:),   intent(inout)  :: p
    integer(i4b),                  intent(out)    :: ierr
    integer(i4b), optional,        intent(in)     :: niter
    real(dp),     optional,        intent(in)     :: tolerance
    interface
       function func(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: func
       end function func
    end interface

    real(dp)        :: epsilon, fret, fpt, fptt, delta, t, fp
    integer(i4b)    :: i, j, k, l, n, ibig, iter
    integer(i4b)    :: numbin_cl, numbin, maxiter
    real(dp), allocatable, dimension(:)    :: pt, ptt, xit
    real(dp), allocatable, dimension(:,:)  :: xi

    ierr    = 0
    epsilon = 1.d-10; if(present(tolerance)) epsilon = tolerance
    maxiter = 100;    if(present(niter))     maxiter = niter
    n       = size(p)


    ! Initilize the directions
    allocate(pt(0:n-1))
    allocate(ptt(0:n-1))
    allocate(xit(0:n-1))
    allocate(xi(0:n-1, 0:n-1))
    xi = 0.d0
    do i = 0, n-1
       xi(i,i) = 1.d0
    end do

    pt = p
    fret = func(p)

    iter = 0
    do while (iter <= maxiter)

       fp = fret
       ibig = 0
       delta = 0.d0
       
       do j = 0, n-1
          xit = xi(:,j)
          fptt = fret
          call linmin(p, xit, fret, func, ierr)
          if (ierr /= 0) return
          if (fptt-fret > delta) then
             delta = fptt-fret
             ibig = j+1
          end if

       end do

       ! The second test here, fret < 1, is *only* for Anderson-Darling fitting
       if (2.d0*(fp-fret) <= epsilon*(abs(fp)+abs(fret))) then
          ierr = 0
          return
       end if

       if (iter == maxiter) then
!          write(*,*) 'Maximum number of iterations reached. Exiting.'
          ierr = 2
          return
       end if

       do j = 0, n-1
          ptt(j) = 2.d0*p(j) - pt(j)
          xit(j) = p(j) - pt(j)
          pt(j)  = p(j)
       end do

       fptt = func(ptt)

       if (fptt < fp) then
          t = 2.d0*(fp-2.d0*fret+fptt) * sqrt(fp-fret-delta) - &
               & delta*sqrt(fp-fptt)

          if (t < 0.d0) then
             call linmin(p, xit, fret, func, ierr)
             if (ierr /= 0) return
             do j = 0, n-1
                xi(j,ibig-1)      = xi(j,n-1)
                xi(j,n-1) = xit(j)
             end do
          end if
       end if

       iter = iter+1

    !write(*,*) 'PARAMS for iter ', iter, ' is ', p
    end do

    deallocate(pt)
    deallocate(ptt)
    deallocate(xi)
    deallocate(xit)

  end subroutine powell


  subroutine linmin(p, xi, fret, func, ierr, xmin_out, eps)
    implicit none

    real(dp), intent(in),  optional :: eps
    real(dp), intent(out), optional :: xmin_out

    real(dp)   :: fret
    real(dp), dimension(:)  :: p, xi
    integer(i4b), intent(out), optional :: ierr
    interface
       function func(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: func
       end function func
    end interface

    integer(i4b)   :: j
    real(dp)       :: epsilon, xx, xmin, fx, fb, fa, bx, ax

    if (present(ierr)) ierr = 0
    epsilon = 1d-4; if (present(eps)) epsilon = eps

    ax = 0.d0
    xx = 1.d0

    call mnbrak(ax, xx, bx, fa, fx, fb, p, xi, func, ierr)
    if (present(ierr)) then
       if (ierr /= 0) return
    end if
    call brent(ax, xx, bx, epsilon, xmin, fret, p, xi, func, ierr)
    if (present(ierr)) then
       if (ierr /= 0) return
    end if

    xi = xi * xmin
    p  = p + xi

    if (present(xmin_out)) xmin_out = xmin

  end subroutine linmin

  subroutine brent(ax, bx, cx, epsilon, xmin, fret, pp, xi, func, ierr)
    implicit none

    real(dp)   :: ax, bx, cx, epsilon, xmin, fret
    real(dp), dimension(:)  :: xi, pp
    integer(i4b), intent(out), optional :: ierr
    interface
       function func(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: func
       end function func
    end interface

    integer(i4b)  :: itmax, iter
    real(dp)      :: cgold, a, b, d, etemp, fu, fv, fw, fx, zeps
    real(dp)      :: p, q, r, tol, tol1, tol2, u, v, w, x, xm, e

    if (present(ierr)) ierr = 0

    d = 0.d0
    e = 0.d0
    cgold = 0.3819660d0
    itmax = 100
    zeps = 1d-6
    tol = 1d-4

    if (ax < cx) then
       a = ax
       b = cx
    else
       a = cx
       b = ax
    end if

    x = bx; w = bx; v = bx
    fw = f1lin(x, pp, xi, func); fv = fw; fx = fw
    do iter = 0, itmax

       xm = 0.5d0*(a+b)
       tol1 = tol*abs(x)+zeps
       tol2 = 2.d0*tol1
       
       if (abs(x-xm) <= (tol2-0.5d0*(b-a))) then
          xmin = x
          fret = fx
          return       ! Quit minimizing
       end if
       
       if (abs(e) > tol1) then
          r = (x-w)*(fx-fv)
          q = (x-v)*(fx-fw)
          p = (x-v)*q-(x-w)*r
          q = 2.d0*(q-r)
          if (q > 0.d0) p = -p
          q = abs(q)
          etemp = e
          e = d
          if ((abs(p) >= abs(0.5d0*q*etemp)) .or. p <= q*(a-x) .or. &
               & p >= q*(b-x)) then
             if (x >= xm) then
                e = a-x
             else
                e = b-x
             end if
             d = cgold * e
          else
             d = p/q
             u = x+d
             if ((u-a < tol2) .or. (b-u < tol2)) then
                d = sign(tol1, xm-x)
             end if
          end if
       else 
          if (x > xm) then
             e = a-x
          else
             e = b-x
          end if
          d = cgold *e
       end if

       if (abs(d) >= tol1) then
          u = x+d
       else
          u = x+sign(tol1,d)
       end if

       fu = f1lin(u,pp,xi,func)

       if (fu <= fx) then
          if (u >= x) then
             a = x
          else
             b = x
          end if

          v = w; w = x; x = u
          fv = fw; fw = fx; fx = fu
       else
          if (u < x) then
             a = u
          else
             b = u
          end if

          if ((fu <= fw) .or. (w == x)) then
             v = w; w = u
             fv = fw; fw = fu
          else if ((fu <= fv) .or. (v == x) .or. (v == w)) then
             v = u
             fv = fu
          end if
       end if
    end do

    if (present(ierr)) then
       ierr = 3
    else
       write(*,*) 'Too many iterations in brent!'
    end if
    xmin = x
    fret = fx

  end subroutine brent

  subroutine mnbrak(ax, bx, cx, fa, fb, fc, p, xi, func, ierr)
    implicit none
    
    real(dp)  :: ax, bx, cx, fa, fb, fc
    real(dp), dimension(:)  :: p, xi
    integer(i4b), optional, intent(out) :: ierr
    interface
       function func(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: func
       end function func
    end interface

    real(dp)  :: gold, glimit, tiny
    real(dp)  :: ulim, u, r, q, fu, rt

    if (present(ierr)) ierr = 0

    gold = 1.618034d0
    tiny = 1.0d-20
    glimit = 100.d0

    fa = f1lin(ax, p, xi, func)
    fb = f1lin(bx, p, xi, func)

    if (fb > fa) then
       rt = ax ! Swap ax and bx
       ax = bx
       bx = rt

       rt = fa
       fa = fb
       fb = rt
    end if
    
    cx = bx + gold*(bx-ax)
    fc = f1lin(cx, p, xi, func)

    do while (fb > fc)

       if ((abs(fb) > upper_bound) .or. (abs(fc) > upper_bound)) then
          ! There is no minimum here..
          if (present(ierr)) then
             ierr = 4
             return
          else
             write(*,*) 'Unable to bracket a minimum. Exiting.'
             stop
          end if
       end if

       r = (bx-ax)*(fb-fc)
       q = (bx-cx)*(fb-fa)
       u = bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),tiny),q-r))
       ulim = bx+glimit*(cx-bx)
       
       if ((bx-u)*(u-cx) > 0.d0) then
          fu = f1lin(u, p, xi, func)
          if (fu < fc) then
             ax = bx
             bx = u
             fa = fb
             fb = fu
             return
          else if (fu > fb) then
             cx = u
             fc = fu
             return
          end if
          u = cx+gold*(cx-bx)
          fu = f1lin(u, p, xi, func)

       else if ((bx-u)*(u-cx) > 0.d0) then
          fu = f1lin(u, p, xi, func)
          if (fu < fc) then
             bx = cx
             cx = u
             u = cx+gold*(cx-bx)

             fb = fc
             fc = fu
             fu = f1lin(u, p, xi, func)
          end if

       else if ((u-ulim)*(ulim-cx) >= 0.d0) then
          u = ulim
          fu = f1lin(u, p, xi, func)
       else
          u = cx+gold*(cx-bx)
          fu = f1lin(u, p, xi, func)
       end if

       ax = bx
       bx = cx
       cx = u

       fa = fb
       fb = fc
       fc = fu
    end do

  end subroutine mnbrak
  
end module powell_mod


