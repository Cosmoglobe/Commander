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
! but WITHout ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
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

  function df1lin(x, p, xi, dfunc)
    implicit none

    real(dp)   :: x
    real(dp), dimension(:)  :: p, xi
    real(dp)                :: df1lin
    interface
       function dfunc(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                                     :: dfunc
       end function dfunc
    end interface

    real(dp), allocatable, dimension(:)  :: xt, df

    allocate(xt(0:size(p)-1), df(0:size(p)-1))
    
    xt = p + x * xi
    df = dfunc(xt)
    df1lin = dot_product(df, xi)
    
    deallocate(xt, df)

  end function df1lin


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

    real(dp)        :: epsilon, fret, fptt, delta, t, fp
    integer(i4b)    :: i, j, n, ibig, iter
    integer(i4b)    :: maxiter
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

  subroutine dbrent(ax, bx, cx, eps, xmin, fret, pp, xi, func, dfunc, ierr)
    implicit none
    real(dp) :: ax, bx, cx, eps, xmin, fret
    real(dp), dimension(:) :: xi, pp
    integer(i4b), intent(out), optional :: ierr
  
    interface
       function func(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: func
       end function func

       function dfunc(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: dfunc
       end function dfunc
    end interface
  
    integer(i4b), parameter :: itmax = 100
    real(dp), parameter :: ZEPS = 1.0e-3_dp * epsilon(ax)
  
    integer(i4b) :: iter
    real(dp) :: a, b, d, d1, d2, du, dv, dw, dx, e, fu, fv, fw, fx, olde, tol, tol1, tol2
    real(dp) :: u, u1, u2, v, w, x, xm
  
    logical :: ok1, ok2  ! Will be used as flags for whether proposed steps are acceptable or not.
  
    a = min(ax, cx)
    b = max(ax, cx)
    v = bx
    w = v
    x = v
    e = 0.0
    tol = 1d-4
    fx = f1lin(x, pp, xi, func)
    fv = fx
    fw = fx
    dx = df1lin(x, pp, xi, dfunc)
  
    ! All our housekeeping chores are doubled by the necessity of moving derivative values around as well as function values.
    dv = dx
    dw = dx
  
    do iter = 1, itmax
      xm = 0.5_dp * (a + b)
      tol1 = tol * abs(x) + ZEPS
      tol2 = 2.0_dp * tol1
  
      if (abs(x - xm) <= (tol2 - 0.5_dp * (b - a))) exit
  
      if (abs(e) > tol1) then
        d1 = 2.0_dp * (b - a)  ! Initialize these d’s to an out-of-bracket value.
        d2 = d1
        if (dw /= dx) d1 = (w - x) * dx / (dx - dw)  ! Secant method with each point.
        if (dv /= dx) d2 = (v - x) * dx / (dx - dv)  ! Which of these two estimates of d shall we take?
        ! We will insist that they be within the bracket, and on the side pointed to by the derivative at x:
        u1 = x + d1
        u2 = x + d2
        ok1 = ((a - u1) * (u1 - b) > 0.0) .AND. (dx * d1 <= 0.0)
        ok2 = ((a - u2) * (u2 - b) > 0.0) .AND. (dx * d2 <= 0.0)
        olde = e  ! Movement on the step before last.
        if (ok1 .or. ok2) then  ! Take only an acceptable d, and if both are acceptable, then take the smallest one.
          if (ok1 .AND. ok2) then
            d = merge(d1, d2, abs(d1) < abs(d2))
          else
            d = merge(d1, d2, ok1)
          end if
          if (abs(d) <= abs(0.5_dp * olde)) then
            u = x + d
            if (u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)
          else
            e = merge(a, b, dx >= 0.0) - x  ! Decide which segment by the sign of the derivative.
            d = 0.5_dp * e  ! Bisect, not golden section.
          end if
        else
          e = merge(a, b, dx >= 0.0) - x  ! Decide which segment by the sign of the derivative.
          d = 0.5_dp * e  ! Bisect, not golden section.
        end if
      end if
  
      if (abs(d) >= tol1) then
        u = x + d
        fu = f1lin(u, pp, xi, func)
      else
        u = x + sign(tol1, d)
        fu = f1lin(u, pp, xi, func)
      end if
  
      ! If the minimum step in the downhill direction takes us uphill, then we are done.
      if (fu > fx) exit
  
      du = df1lin(u, pp, xi, dfunc)
  
      ! Now all the housekeeping, sigh.
      if (fu <= fx) then
        if (u >= x) then
          a = x
        else
          b = x
        end if
        call mov3(v, fv, dv, w, fw, dw)
        call mov3(w, fw, dw, x, fx, dx)
        call mov3(x, fx, dx, u, fu, du)
      else
        if (u < x) then
          a = u
        else
          b = u
        end if
        if (fu <= fw .or. w == x) then
          call mov3(v, fv, dv, w, fw, dw)
          call mov3(w, fw, dw, u, fu, du)
        else if (fu <= fv .or. v == x .or. v == w) then
          call mov3(v, fv, dv, u, fu, du)
        end if
      end if
    end do
  
    if (iter > itmax) then
      if (present(ierr)) then
        ierr = 3
      else
        write(*,*) 'dbrent: exceeded maximum iterations'
      end if
    end if
    xmin = x
    fret = fx
  
  contains
  
    subroutine mov3(a, b, c, d, e, f)
      real(dp), intent(in) :: d, e, f
      real(dp), intent(out) :: a, b, c
      a = d
      b = e
      c = f
    end subroutine mov3
  
  end subroutine dbrent
  
  
  
  
  
  subroutine dlinmin(p, xi, fret, func, dfunc, ierr, xmin_out, eps)
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
  
       function dfunc(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                           :: dfunc
       end function dfunc
    end interface
  
    real(dp)       :: epsilon, xx, xmin, fx, fb, fa, bx, ax
  
    if (present(ierr)) ierr = 0
    epsilon = 1d-4; if (present(eps)) epsilon = eps
  
    ax = 0.d0
    xx = 1.d0
  
    call mnbrak(ax, xx, bx, fa, fx, fb, p, xi, func, ierr)
    if (present(ierr)) then
       if (ierr /= 0) return
    end if
    call dbrent(ax, xx, bx, epsilon, xmin, fret, p, xi, func, dfunc, ierr)
    if (present(ierr)) then
       if (ierr /= 0) return
    end if
  
    xi = xi * xmin
    p  = p + xi
  
    if (present(xmin_out)) xmin_out = xmin
  
  end subroutine dlinmin
  
end module powell_mod
