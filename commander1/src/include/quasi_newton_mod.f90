module quasi_newton_mod
  use healpix_types
  implicit none


contains

  function outerprod(x, y) result(z)
    implicit none

    real(dp), dimension(:)   :: x, y
    real(dp), dimension(size(x),size(y)) :: z

    integer(i4b) :: i, j

    do i = 1, size(x)
       do j = 1, size(y)
          z(i,j) = x(i)*y(j)
       end do
    end do

  end function outerprod


  SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func,ierr)
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN) :: xold,g
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
    REAL(dp), INTENT(IN) :: fold,stpmax
    REAL(dp), DIMENSION(:), INTENT(OUT) :: x
    integer(i4b),           intent(out), optional :: ierr
    REAL(dp), INTENT(OUT) :: f
    LOGICAL(LGT), INTENT(OUT) :: check
    INTERFACE
       FUNCTION func(x)
         USE healpix_types
         IMPLICIT NONE
         REAL(dp) :: func
         REAL(dp), DIMENSION(:), INTENT(IN), optional :: x
       END FUNCTION func
    END INTERFACE
    REAL(dp), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
    INTEGER(I4B) :: ndum
    REAL(dp) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
!    ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
    check=.false.
    pabs=sqrt(sum(p(:)**2))
    if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
    slope=sum(g*p)
    !if (slope >= 0.d0) write(*,*) 'roundoff problem in lnsrch', slope
    alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))
    alam=1.d0
    do
       x(:)=xold(:)+alam*p(:)
       f=func(x)
       if (f /= f) then
          if (present(ierr)) ierr = 1
          return
       end if
       if (alam < alamin) then
          x(:)=xold(:)
          check=.true.
          RETURN
       else if (f <= fold+ALF*alam*slope) then
          RETURN
       else
          if (alam == 1.d0) then
             tmplam=-slope/(2.0_dp*(f-fold-slope))
          else
             rhs1=f-fold-alam*slope
             rhs2=f2-fold-alam2*slope
             a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
             b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                  (alam-alam2)
             if (a == 0.d0) then
                tmplam=-slope/(2.0_dp*b)
             else
                disc=b*b-3.0_dp*a*slope
                if (disc < 0.d0) then
                   tmplam=0.5_dp*alam
                else if (b <= 0.d0) then
                   tmplam=(-b+sqrt(disc))/(3.0_dp*a)
                else
                   tmplam=-slope/(b+sqrt(disc))
                end if
             end if
             if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
          end if
       end if
       alam2=alam
       f2=f
       alam=max(tmplam,0.1_dp*alam)
    end do
  END SUBROUTINE lnsrch






  SUBROUTINE dfpmin(p,gtol,iter,fret,func,dfunc,ierr)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(OUT) :: iter
    REAL(dp), INTENT(IN) :: gtol
    REAL(dp), INTENT(OUT) :: fret
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
    integer(i4b), optional, intent(out) :: ierr
    INTERFACE
       FUNCTION func(p)
         USE healpix_types
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN), optional :: p
         REAL(dp) :: func
       END FUNCTION func
       FUNCTION dfunc(p)
         USE healpix_types
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: p
         REAL(dp), DIMENSION(size(p)) :: dfunc
       END FUNCTION dfunc
    END INTERFACE
    INTEGER(I4B), PARAMETER :: ITMAX=1000
    REAL(dp), PARAMETER :: STPMX=100.0_dp,EPS=epsilon(p),TOLX=4.0_dp*EPS
    INTEGER(I4B) :: its, i, j, ierr_int
    LOGICAL :: check
    REAL(dp) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
    REAL(dp), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi
    REAL(dp), DIMENSION(size(p),size(p)) :: hessin
    if (present(ierr)) ierr = 0
    fp=func(p)
    if (fp /= fp) then
       if (present(ierr)) ierr = 1
       return
    end if
    g=dfunc(p)
    hessin = 0.d0
    do i = 1, size(p)
       hessin(i,i) = 1.d0
    end do
    xi=-g
    stpmax=STPMX*max(sqrt(sum(p*p)),real(size(p),dp))
    do its=1,ITMAX
       iter=its
       call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func,ierr_int)
       if (fret /= fret) then
          if (present(ierr)) ierr = ierr_int
          return
       end if
       fp=fret
       xi=pnew-p
       p=pnew
       if (maxval(abs(xi)/max(abs(p),1.0_dp)) < TOLX) RETURN
       dg=g
       g=dfunc(p)
       den=max(abs(fret),1.0_dp)
       if (maxval(abs(g)*max(abs(p),1.0_dp)/den) < gtol) RETURN
       dg=g-dg
       hdg=matmul(hessin,dg)
       fac=sum(dg*xi)
       fae=sum(dg*hdg)
       sumdg=sum(dg*dg)
       sumxi=sum(xi*xi)
       if (fac > sqrt(EPS*sumdg*sumxi)) then
          fac=1.0_dp/fac
          fad=1.0_dp/fae
          dg=fac*xi-fad*hdg
          hessin=hessin+fac*outerprod(xi,xi)-&
               fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
       end if
       xi=-matmul(hessin,g)
    end do
    if (present(ierr)) then
       ierr = 2
    else
       write(*,*) 'dfpmin: too many iterations'
    end if
    
  END SUBROUTINE dfpmin





!!$
!!$
!!$
!!$  subroutine lnsrch(n, xold, fold, g, p, x, f, stpmax, check, func)
!!$    implicit none
!!$
!!$    integer(i4b) :: n
!!$    logical(lgt) :: check
!!$    real(dp)     :: f, fold, stpmax, g(n), p(n), x(n), xold(n), func, ALF, TOLX
!!$    parameter (ALF=1.d-4, TOLX=1.d-7)
!!$    external func
!!$
!!$    integer(i4b) :: i
!!$    real(dp)     :: a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp, test, tmplam
!!$    
!!$    check = .false.
!!$    sum   = 0.d0
!!$    do i = 1, n
!!$       sum = sum + p(i)*p(i)
!!$    end do
!!$    sum = sqrt(sum)
!!$    if (sum > stpmax) then
!!$       do i = 1, n
!!$          p(i) = p(i) * stpmax/sum
!!$       end do
!!$    end if
!!$    slope = 0.d0
!!$    do i = 1, n
!!$       slope = slope + g(i)*p(i)
!!$    end do
!!$    if (slope > 0.d0) write(*,*) 'Roundoff problem in lnsrch'
!!$    test = 0.d0
!!$    do i = 1, n
!!$       temp = abs(p(i)) / max(abs(xold(i)),1.d0)
!!$       if (temp > test) test = temp
!!$    end do
!!$    alamin = TOLX / test
!!$    alam = 1.d0
!!$1   continue
!!$    do i = 1, n
!!$       x(i) = xold(i) + alam*p(i)
!!$    end do
!!$    f = func(x)
!!$    if (alam < alamin) then
!!$       do i = 1, n
!!$          x(i) = xold(i)
!!$       end do
!!$       check = .true.
!!$       return
!!$    else if ( f < fold + ALF*alam*slope) then
!!$       return
!!$    else 
!!$       if (alam == 1.d0) then
!!$          tmplam = -slope/(2.d0*(f-fold-slope))
!!$       else
!!$          rhs1 = f-fold-alam*slope
!!$          rhs2 = f2-fold-alam2*slope
!!$          a = (rhs1/alam**2-rhs2/alam2**2) / (alam-alam2)
!!$          b = (-alam2*rhs1/alam**2 + alam*rhs2/alam2**2) / (alam-alam2)
!!$          if (a == 0.d0) then
!!$             tmplam = -slope / (2.d0*b)
!!$          else
!!$             disc = b*b - 3.d0*a*slope
!!$             if (disc < 0.d0) then
!!$                tmplam = 0.5d0*alam
!!$             else if (b < 0.d0) then
!!$                tmplam = (-b+sqrt(disc)) / (3.d0*a)
!!$             else
!!$                tmplam = -slope / (b+sqrt(disc))
!!$             end if
!!$          end if
!!$          if (tmplam > 0.5d0*alam) tmplam = 0.5d0*alam
!!$       end if
!!$    end if
!!$    alam2 = alam
!!$    f2 = f
!!$    alam = max(tmplam, 0.1d0*alam)
!!$    goto 1
!!$  end subroutine lnsrch
!!$
!!$
!!$
!!$  subroutine dfpmin(p, n, gtol, iter, fret, func, dfunc)
!!$    implicit none
!!$
!!$    integer(i4b) :: iter, n, NMAX, ITMAX
!!$    real(dp)     :: fret, gtol, p(n), EPS, STPMX, TOLX, func
!!$    parameter (NMAX=50, ITMAX=200, STPMX=100.d0, EPS = 3.d-8, TOLX=4.d0*EPS)
!!$    external func, dfunc
!!$    integer(i4b) :: i, its, j
!!$    logical(lgt) :: check
!!$    real(dp)     :: den, fac, fad, fae, fp, stpmax, sum, sumdg, sumxi, temp, test
!!$    real(dp)     :: dg(NMAX), g(NMAX), hdg(NMAX), hessin(NMAX,NMAX), pnew(NMAX), xi(NMAX)
!!$
!!$    write(*,*) p
!!$    write(*,*) size(p)
!!$    fp = func(p)
!!$    write(*,*) fp
!!$    call dfunc(p, g)
!!$    write(*,*) g
!!$    sum = 0.d0
!!$    do i = 1, n
!!$       do j = 1, n
!!$          hessin(i,i) = 1.d0
!!$       end do
!!$       xi(i)       = -g(i)
!!$       sum         = sum+p(i)**2
!!$    end do
!!$    stpmax = STPMX*max(sqrt(sum), float(n))
!!$    do its = 1, ITMAX
!!$       iter = its
!!$       call lnsrch(n, p, fp, g, xi, pnew, fret, stpmax, check, func)
!!$       fp = fret
!!$       do i = 1, n
!!$          xi(i) = pnew(i) - p(i)
!!$          p(i)  = pnew(i)
!!$       end do
!!$       test = 0.d0
!!$       do i = 1, n
!!$          temp = abs(xi(i)) / max(abs(p(i)), 1.d0)
!!$          if (temp > test) test = temp
!!$       end do
!!$       if (test < TOLX) return
!!$       do i = 1, n
!!$          dg(i) = g(i)
!!$       end do
!!$       call dfunc(p, g)
!!$       test = 0.d0
!!$       den = max(fret, 1.d0)
!!$       do i = 1, n
!!$          temp = abs(g(i)) * max(abs(p(i)), 1.d0) / den
!!$          if (temp > test) test = temp
!!$       end do
!!$       if (test < gtol) return
!!$       do i = 1, n
!!$          dg(i) = g(i) - dg(i)
!!$       end do
!!$       do i = 1, n
!!$          hdg(i) = 0.d0
!!$          do j = 1, n
!!$             hdg(i) = hdg(i) + hessin(i,j) * dg(j)
!!$          end do
!!$       end do
!!$       fac = 0.d0
!!$       fae = 0.d0
!!$       sumdg = 0.d0
!!$       sumxi = 0.d0
!!$       do i = 1, n
!!$          fac = fac + dg(i) * xi(i)
!!$          fae = fae + dg(i) * hdg(i)
!!$          sumdg = sumdg + dg(i)**2
!!$          sumxi = sumxi + xi(i)**2
!!$       end do
!!$       if (fac > sqrt(EPS*sumdg*sumxi)) then
!!$          fac = 1.d0 / fac
!!$          fad = 1.d0 / fae
!!$          do i = 1, n
!!$             dg(i) = fac*xi(i) - fad*hdg(i)
!!$          end do
!!$          do i = 1, n
!!$             do j = i, n
!!$                hessin(i,j) = hessin(i,j) + fac*xi(i)*xi(j) &
!!$                     & -fad*hdg(i)*hdg(j) + fae*dg(i)*dg(j)
!!$                hessin(j,i) = hessin(i,j)
!!$             end do
!!$          end do
!!$       end if
!!$       do i = 1, n
!!$          xi(i) = 0.d0
!!$          do j = 1, n
!!$             xi(i) = xi(i) - hessin(i,j)*g(j)
!!$          end do
!!$       end do
!!$    end do
!!$    write(*,*) 'Too many iterations in dfpmin'
!!$    return
!!$
!!$  end subroutine dfpmin


end module quasi_newton_mod
