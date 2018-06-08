module comm_like_utils
  use healpix_types
  use math_tools
  implicit none

  ! *********************************************************************
  ! *  comm_like_mod -- Support module for Commander likelihood modules *
  ! *                                                                   *
  ! *          H. K. Eriksen, E. Gjerløw (University of Oslo)           *
  ! *                                                                   *   
  ! *                               and                                 *
  ! *                                                                   *   
  ! *                       I. K. Wehus  (JPL)                          *   
  ! *                                                                   *
  ! *  History:                                                         *
  ! *      January 10th, 2014 -- First fully functional version         *
  ! *                                                                   *
  ! *********************************************************************

!!$  integer, parameter, private :: i4b = selected_int_kind(9)
!!$  integer, parameter, private :: sp  = selected_real_kind(5,30)
!!$  integer, parameter, private :: dp  = selected_real_kind(12,200)
!!$  integer, parameter, private :: lgt = kind(.true.)
!!$
!!$  real(dp), parameter, private :: pi = 3.141592653589793238462643383279502d0

  ! Spectrum order definition
  integer(i4b), parameter :: TT = 1
  integer(i4b), parameter :: TE = 2
  integer(i4b), parameter :: TB = 3
  integer(i4b), parameter :: EE = 4
  integer(i4b), parameter :: EB = 5
  integer(i4b), parameter :: BB = 6
  character(len=2), dimension(6), parameter :: spec_label = ['TT','TE','TB','EE','EB','BB']

contains

  subroutine get_invfisher(cls, w_l, PY, N_cov, bins, invF, inv_error)
    implicit none

    real(dp),     dimension(0:,1:),    intent(in)  :: cls
    real(dp),     dimension(0:,1:,1:), intent(in)  :: w_l
    real(dp),     dimension(1:,1:),    intent(in)  :: PY, N_cov
    integer(i4b), dimension(1:,1:),    intent(in)  :: bins
    real(dp),     dimension(1:,1:),    intent(out) :: invF
    logical(lgt),                      intent(out) :: inv_error
    
    integer(i4b) :: i, j, k, lmax, ierr, n, n_h, l, m, ind, col, row, ncomp, s, npar, nmode, nspec
    integer(i4b) :: nmaps
    real(dp)     :: t1, t2, t3, t4, var
    real(dp), allocatable, dimension(:,:) :: cls_binned, C, invC, d, Pt_invC_D, A, invC_P, Cl
    real(dp), allocatable, dimension(:,:) :: Pt_invC_P, Pt_invC_P_dC, map, V, test, P_sub, PYS
    real(dp), allocatable, dimension(:)   :: W
    integer(i4b), allocatable, dimension(:) :: indmap
    integer(i4b), dimension(6) :: ind1 = [1, 1, 1, 2, 2, 3]
    integer(i4b), dimension(6) :: ind2 = [1, 2, 3, 2, 3, 3]
    real(dp),     dimension(6) :: mult = [1.d0, 2.d0, 2.d0, 1.d0, 2.d0, 1.d0]

    ! Compute inverse covariance matrix
    n       = size(PY,1)
    n_h     = size(PY,2)
    lmax    = size(w_l,1)-1
    nmaps   = size(w_l,2)
    nspec   = size(cls,2)
    ncomp   = n_h / nmaps
    npar    = size(bins,1)

    allocate(cls_binned(0:lmax,nspec))
    cls_binned = cls
    do i = 1, npar
       cls_binned(bins(i,1):bins(i,2),bins(i,3)) = sum(cls(bins(i,1):bins(i,2),bins(i,3)))/&
            & (bins(i,2)-bins(i,1)+1)
    end do

    ! Compute C = S + N
    allocate(C(n,n), PYs(n,n_h), Cl(nmaps,nmaps))
    PYS            = 0.d0
    ind            = 5
    do l = 2, lmax
       call cl2s(cls_binned(l,:), Cl)      
       Cl = Cl * w_l(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
       do m = -l, l
          call dgemm('N', 'N', n, nmaps, nmaps, 1.d0, PY(:,ind:n_h:ncomp), n, &
               & Cl, nmaps, 0.d0, PYS(:,ind:n_h:ncomp), n)
          ind  = ind  + 1
       end do
    end do
    call dgemm('N', 'T', n, n, n_h, 1.d0, PYS, n, PY, n, 0.d0, C, n)
    C = C + N_cov
    deallocate(PYS, Cl)

    ! Precompute P^t * invC * P
    allocate(invC_P(n,n_h), Pt_invC_P(n_h,n_h), A(n_h,n_h), Pt_invC_P_dC(n_h,n_h), indmap(n_h))

!    call invert_matrix(C)
!    Pt_invC_P = matmul(transpose(PY), matmul(C, PY))
    call wall_time(t1)
    call dpotrf('L', n, C, n, ierr)
    invC_P = PY
    call dpotrs('L', n, n_h, C, n, invC_P, n, ierr)    
    call dgemm('T', 'N', n_h, n_h, n, 1.d0, PY, n, &
         & invC_P, n, 0.d0, Pt_invC_P, n_h)
    call wall_time(t2)
    write(*,*) '     Computing P*invC*P^t -- wall time = ', real(t2-t1,sp), ' sec'

    ! Compute Fisher matrix
    invF = 0.d0
    call wall_time(t1)
    do i = 1, npar

       call wall_time(t3)
       ! Precompute Pt*invC*P*dC
       s   = bins(i,3)
       row  = (ind1(s)-1)*ncomp + bins(i,1)**2 + 1
       col  = (ind2(s)-1)*ncomp + bins(i,1)**2 + 1
       Pt_invC_P_dC = 0.d0
       do l = bins(i,1), bins(i,2)
          do m = -l, l
             do j = 1, n_h
                Pt_invC_P_dC(j,col) = Pt_invC_P(j,row) * 2.d0*pi/real(l*(l+1),dp) * &
                     & w_l(l,ind1(s),ind2(s))
                if (row /= col) then
                   Pt_invC_P_dC(j,row) = Pt_invC_P(j,col) * 2.d0*pi/real(l*(l+1),dp) * &
                        & w_l(l,ind1(s),ind2(s))
                end if
             end do
             row = row + 1
             col = col + 1
          end do
       end do

       nmode = 0
       do j = 1, n_h
          if (any(Pt_invC_P_dC(:,j) /= 0.d0)) then
             nmode         = nmode+1
             indmap(nmode) = j
          end if
       end do
       if (nmode == 0) then
          inv_error = .True.
          write(*, *) 'No modes'
          return
       end if
       call wall_time(t4)
!       write(*,*) 'a', t4-t3

       call wall_time(t3)
       call dgemm('N', 'N', n_h, n_h, nmode, 1.d0, Pt_invC_P_dC(:,indmap(1:nmode)), n_h, &
            & Pt_invC_P(indmap(1:nmode),:), nmode, 0.d0, A, n_h)
       call wall_time(t4)
!       write(*,*) 'b', t4-t3

       call wall_time(t3)
       do j = i, npar
          s   = bins(j,3)
          row = (ind1(s)-1)*ncomp + bins(j,1)**2 + 1
          col = (ind2(s)-1)*ncomp + bins(j,1)**2 + 1
          do l = bins(j,1), bins(j,2)
             do m = -l, l
                invF(i,j) = invF(i,j) + 0.5d0 * mult(s) * A(row,col) * w_l(l,ind1(s),ind2(s)) * &
                     & 2.d0*pi/real(l*(l+1),dp)
                row = row + 1
                col = col + 1
             end do
          end do
          invF(j,i) = invF(i,j)
       end do
       call wall_time(t4)
!       write(*,*) 'c', t4-t3

    end do
    call wall_time(t2)
    write(*,fmt='(a,f8.3,a)') &
         & '      Computing Fisher matrix -- wall time = ', real(t2-t1,sp), ' sec'

!!$    allocate(W(npar))
!!$    call get_eigenvalues(invF, W)
!!$    write(*,fmt='(a,3e16.8)') 'eigen d = ', maxval(W)/minval(W), minval(W), maxval(W)
!!$    deallocate(W)

    allocate(W(npar))
    call get_eigenvalues(invF, W)
    if (any(W <= 0.d0) .or. W(npar)/W(1) > 1.d12) then
       write(*,*) 'Error -- Fisher matrix not invertible'
       write(*,*) '         Condition number = ', W(npar)/W(1)
       write(*,*) '         Check cls and S+N condition number.'
       open(48,file='fisher_W.dat')
       do i = 1, npar
          write(48,*) i, W(i)
       end do
       close(48)

       open(48,file='fisher_cls.dat', recl=1024)
       do l = 2, lmax
          write(48,fmt='(i8,6f10.3)') l, cls(l,:)
       end do
       close(48)

       ! Output error bars
       open(48,file='fisher_full.dat', recl=10000)
       do j = 1, npar
          write(48,fmt='(i8)',advance='no') j
          do i = 1, npar
             write(48,fmt='(e16.8)',advance='no') invF(j,i)
          end do
          write(48,*)
       end do
       close(48)

       inv_error = .True.
    else
      inv_error = .False.
      call invert_matrix(invF)
    end if

    deallocate(cls_binned, invC_P, Pt_invC_P, A, indmap, W, Pt_invC_P_dC)

  end subroutine get_invfisher

  subroutine get_invfisher_safe(cls, w_l, PY, N_cov, bins, invF)
    implicit none

    real(dp),     dimension(0:,1:),    intent(in)  :: cls
    real(dp),     dimension(0:,1:,1:), intent(in)  :: w_l
    real(dp),     dimension(1:,1:),    intent(in)  :: PY, N_cov
    integer(i4b), dimension(1:,1:),    intent(in)  :: bins
    real(dp),     dimension(1:,1:),    intent(out) :: invF
    
    integer(i4b) :: i, j, k, lmax, ierr, n, n_h, l, m, ind, col, row, ncomp, s, npar, nmode, nspec
    integer(i4b) :: nmaps
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:,:)   :: cls_binned, C, d, A, Cl
    real(dp), allocatable, dimension(:,:)   :: map, V, test, PYS, dC, buffer, dC_n
    real(dp), allocatable, dimension(:,:,:) :: invC_dC
    real(dp), allocatable, dimension(:)   :: W
    integer(i4b), allocatable, dimension(:) :: indmap
    integer(i4b), dimension(6) :: ind1 = [1, 1, 1, 2, 2, 3]
    integer(i4b), dimension(6) :: ind2 = [1, 2, 3, 2, 3, 3]
    real(dp),     dimension(6) :: mult = [1.d0, 2.d0, 2.d0, 1.d0, 2.d0, 1.d0]

    ! Compute inverse covariance matrix
    n       = size(PY,1)
    n_h     = size(PY,2)
    lmax    = size(w_l,1)-1
    nmaps   = size(w_l,2)
    nspec   = size(cls,2)
    ncomp   = n_h / nmaps
    npar    = size(bins,1)

    allocate(cls_binned(0:lmax,nspec))
    cls_binned = cls
    do i = 1, npar
       cls_binned(bins(i,1):bins(i,2),bins(i,3)) = sum(cls(bins(i,1):bins(i,2),bins(i,3)))/&
            & (bins(i,2)-bins(i,1)+1)
    end do

    ! Compute C = S + N
    allocate(C(n,n), PYS(n,n_h), Cl(nmaps,nmaps))
    PYS            = 0.d0
    ind            = 5
    do l = 2, lmax
       call cl2s(cls_binned(l,:), Cl)      
       Cl = Cl * w_l(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
       do m = -l, l
          call dgemm('N', 'N', n, nmaps, nmaps, 1.d0, PY(:,ind:n_h:ncomp), n, &
               & Cl, nmaps, 0.d0, PYS(:,ind:n_h:ncomp), n)
          ind  = ind  + 1
       end do
    end do
    call dgemm('N', 'T', n, n, n_h, 1.d0, PYS, n, PY, n, 0.d0, C, n)
    C = C + N_cov
    deallocate(PYS, Cl)

    ! Invert covariance
    call invert_matrix(C)

    ! Precompute invC * dC
    allocate(invC_dC(n,n,npar), dC(n_h,n_h), buffer(n,n_h), dC_n(n,n))
    do i = 1, npar
       write(*,*) i, npar
       s    = bins(i,3)
       row  = (ind1(s)-1)*ncomp + bins(i,1)**2 + 1
       col  = (ind2(s)-1)*ncomp + bins(i,1)**2 + 1
       dC   = 0.d0
       do l = bins(i,1), bins(i,2)
          do m = -l, l
             dC(row,col) = 2.d0*pi/real(l*(l+1),dp) * w_l(l,ind1(s),ind2(s))
             dC(col,row) = dC(row,col)
             row = row + 1
             col = col + 1
          end do
       end do
       !dC             = matmul(PY, matmul(dC, transpose(PY)))
       call dgemm('N', 'N', n, n_h, n_h, 1.d0, PY, n, dC, n_h, 0.d0, buffer, n)
       call dgemm('N', 'T', n, n, n_h, 1.d0, buffer, n, PY, n, 0.d0, dC_n, n)
       !invC_dC(:,:,i) = matmul(C, dC)
       call dgemm('N', 'N', n, n, n, 1.d0, C, n, dC_n, n, 0.d0, invC_dC(:,:,i), n)
    end do
    deallocate(buffer, dC, dC_n)


    ! Compute Fisher matrix
    invF = 0.d0
    call wall_time(t1)
    do i = 1, npar
       write(*,*) i, npar
       do j = 1, npar
          do k = 1, n
             invF(i,j) = invF(i,j) + 0.5d0 * sum(invC_dC(k,:,i)*invC_dC(:,k,j))
          end do
       end do
    end do
    call wall_time(t2)
    write(*,fmt='(a,f8.3,a)') &
         & '      Computing Fisher matrix -- wall time = ', real(t2-t1,sp), ' sec'

    allocate(W(npar))
    call get_eigenvalues(invF, W)
    if (any(W <= 0.d0) .or. W(npar)/W(1) > 1.d12) then
       write(*,*) 'Error -- Fisher matrix not invertible'
       write(*,*) '         Condition number = ', W(npar)/W(1)
       write(*,*) '         Either remove spectrum bins, or add more basis vectors'
       open(48,file='fisher_W.dat')
       do i = 1, npar
          write(48,*) i, W(i)
       end do
       close(48)

       open(48,file='fisher_cls.dat', recl=1024)
       do l = 2, lmax
          write(48,fmt='(i8,6f10.3)') l, cls(l,:)
       end do
       close(48)

       ! Output error bars
       open(48,file='fisher_full.dat', recl=10000)
       do j = 1, npar
          write(48,fmt='(i8)',advance='no') j
          do i = 1, npar
             write(48,fmt='(e16.8)',advance='no') invF(j,i)
          end do
          write(48,*)
       end do
       close(48)

       stop
    end if
    call invert_matrix(invF)

    deallocate(cls_binned, invC_dC, W)

  end subroutine get_invfisher_safe

  subroutine get_invfisher_analytic(cls, w_l, N_l, bins, invF)
    implicit none

    real(dp),     dimension(0:,1:),    intent(in)  :: cls
    real(dp),     dimension(0:,1:,1:), intent(in)  :: w_l
    real(dp),                          intent(in)  :: N_l
    integer(i4b), dimension(1:,1:),    intent(in)  :: bins
    real(dp),     dimension(1:,1:),    intent(out) :: invF
    
    integer(i4b) :: i, j, k, lmax, ierr, n, n_h, l, m, ind, col, row, ncomp, s, npar, nmode, nspec
    integer(i4b) :: nmaps
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:,:)   :: cls_binned, C, d, A, Cl
    real(dp), allocatable, dimension(:,:)   :: map, V, test, PYS, dC, buffer, dC_n
    real(dp), allocatable, dimension(:,:,:) :: invC_dC
    real(dp), allocatable, dimension(:)   :: W
    integer(i4b), allocatable, dimension(:) :: indmap
    integer(i4b), dimension(6) :: ind1 = [1, 1, 1, 2, 2, 3]
    integer(i4b), dimension(6) :: ind2 = [1, 2, 3, 2, 3, 3]
    real(dp),     dimension(6) :: mult = [1.d0, 2.d0, 2.d0, 1.d0, 2.d0, 1.d0]

    ! Compute inverse covariance matrix
    lmax    = size(w_l,1)-1
    nmaps   = 1 !size(w_l,2)
    nspec   = 1 !size(cls,2)
    npar    = size(bins,1)
    n_h     = nmaps*(lmax+1)**2

    allocate(cls_binned(0:lmax,nspec))
    cls_binned = cls(:,1:1)
!    do i = 1, npar
!       cls_binned(bins(i,1):bins(i,2),bins(i,3)) = sum(cls(bins(i,1):bins(i,2),bins(i,3)))/&
!            & (bins(i,2)-bins(i,1)+1)
!    end do

    ! Compute C = S + N
    allocate(C(n_h,n_h), dC(n_h,n_h), Cl(nmaps,nmaps))
    C  = 0.d0
    ind = 1
    do l = 0, lmax
       Cl = 0.d0
       if (l >= 2) then
          call cl2s(cls_binned(l,:), Cl)      
          Cl = Cl * w_l(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
       end if
       write(*,*) real(Cl,sp)
       do m = -l, l
          do i = 1, nmaps
             do j = 1, nmaps
                if (i == j) then
                   C((i-1)*nmaps+ind,(j-1)*nmaps+ind) = Cl(i,j) + N_l
                else
                   C((i-1)*nmaps+ind,(j-1)*nmaps+ind) = Cl(i,j)
                end if
             end do
          end do
          ind = ind+1
       end do
    end do

    ind = 1
    do l = 0, lmax
       do m = -l, l
          write(*,fmt='(i6,i6,a,3f16.8)') l, m, ', C = ', C(ind,ind), sum(abs(C(ind,:))), sum(abs(C(:,ind)))
          ind = ind+1
       end do
    end do
    
    ! Invert covariance
    call invert_matrix(C)

    ind = 1
    do l = 0, lmax
       do m = -l, l
          write(*,*) l, m, ', invC = ', C(ind,ind)
          ind = ind+1
       end do
    end do

    ! Precompute invC * dC
    allocate(invC_dC(n_h,n_h,npar))
    do i = 1, npar
       write(*,*) i, npar
       s    = bins(i,3)
       row  = (ind1(s)-1)*ncomp + bins(i,1)**2 + 1
       col  = (ind2(s)-1)*ncomp + bins(i,1)**2 + 1
       dC   = 0.d0
       do l = bins(i,1), bins(i,2)
          do m = -l, l
             dC(row,col) = 2.d0*pi/real(l*(l+1),dp) * w_l(l,ind1(s),ind2(s))
             dC(col,row) = dC(row,col)
             row         = row + 1
             col         = col + 1
          end do
       end do
       invC_dC(:,:,i) = matmul(C, dC)
    end do
    write(*,*) 'invC_dC = ', invC_dC

    ! Compute Fisher matrix
    invF = 0.d0
    call wall_time(t1)
    do i = 1, npar
       write(*,*) i, npar
       do j = 1, npar
          do k = 1, n_h
             invF(i,j) = invF(i,j) + 0.5d0 * sum(invC_dC(k,:,i)*invC_dC(:,k,j))
          end do
       end do
    end do
    call wall_time(t2)
    write(*,fmt='(a,f8.3,a)') &
         & '      Computing Fisher matrix -- wall time = ', real(t2-t1,sp), ' sec'

    write(*,*) 'F = ', invF
    call invert_matrix(invF)
    write(*,*) 'invF = ', invF

    deallocate(cls_binned, invC_dC, dC)

  end subroutine get_invfisher_analytic

  ! **********************************************************************
  !                       File utilities
  ! **********************************************************************
  ! *****************************************************************************************
  !
  ! Routine for reading one parameter from a parameter file
  !     Example usage:   call comm_get_parameter(21, "parfile.txt", "NSIDE", par_int=nside)
  !
  ! *****************************************************************************************

  subroutine comm_get_parameter(parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt)
    implicit none

    character(len=*),  intent(in)  :: parfile, parname
    integer(i4b),      intent(out), optional :: par_int
    logical(lgt),      intent(out), optional :: par_lgt
    character(len=1),  intent(out), optional :: par_char
    character(len=*),  intent(out), optional :: par_string
    real(sp),          intent(out), optional :: par_sp
    real(dp),          intent(out), optional :: par_dp


    integer(i4b)        :: i, unit
    character(len=128)  :: string, variable, value
    character(len=1)    :: equals

    unit = comm_getlun()
    open(unit, file=trim(parfile))

    do while (.true.)
       read(unit,*,end=1) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       read(unit,*) variable, equals, value

       if (trim(variable) == trim(parname)) then

          if (present(par_int)) then
             read(value,*) par_int
          else if (present(par_char)) then
             read(value,*) par_char
          else if (present(par_string)) then
             read(value,'(a)') par_string
          else if (present(par_sp)) then
             read(value,*) par_sp
          else if (present(par_dp)) then
             read(value,*) par_dp
          else if (present(par_lgt)) then
             read(value,*) par_lgt
          end if

          close(unit)
          return

       end if

    end do

1   write(*,*) 'GET_PARAMETER:    Critical error -- parameter not found'
    write(*,*) 'GET_PARAMETER:       Parameter file = ', trim(parfile)
    write(*,*) 'GET_PARAMETER:       Parameter name = ', trim(parname)

    close(unit)
    stop

  end subroutine comm_get_parameter

  function comm_getlun()
    implicit none
    integer(i4b) :: comm_getlun
    logical(lgt) :: exists, isopen
    comm_getlun = 9
    do
       comm_getlun = comm_getlun+1
       inquire(unit=comm_getlun,exist=exists)
       if(exists) then
          inquire(unit=comm_getlun,opened=isopen)
          if(.not. isopen) return
       end if
    end do
  end function comm_getlun
  
  ! Routine for reading the Gibbs sigma samples 
  subroutine comm_read_chain(filename, unit, lmax, numchains, numsamples, data)
    implicit none

    character(len=*),                              intent(in)  :: filename
    integer(i4b),                                  intent(in)  :: unit
    integer(i4b),                                  intent(out) :: lmax, numchains, numsamples
    real(dp),         allocatable, dimension(:,:,:,:)          :: data

    integer(i4b)         :: l, status, blocksize, readwrite, numspec, i, j, k
    integer(i4b)         :: fpixel, group, numargs
    logical(lgt)         :: anyf
    real(dp)             :: nullval
    character(len=80)    :: comment
    integer(i4b),          dimension(4)     :: naxes
    real(dp),     pointer, dimension(:,:,:,:) :: indata

    status = 0
    readwrite = 0
    nullval = 0.

    numargs = 1
    !numargs = 0

    ! Open the result file
    call ftopen(unit,trim(filename),readwrite,blocksize,status)

    ! Read keywords
    call ftgkyj(unit,'LMAX',     lmax,       comment,status)
    call ftgkyj(unit,'NUMSAMP',  numsamples, comment,status)
    call ftgkyj(unit,'NUMCHAIN', numchains,  comment,status)
    call ftgkyj(unit,'NUMSPEC',  numspec,    comment,status)

    allocate(data(0:lmax,numspec,numchains,numsamples))
    allocate(indata(0:lmax,numspec,numchains,numargs+numsamples))

    ! Read the binned power spectrum array
    group  = 1
    fpixel = 1
    call ftgpvd(unit,group,fpixel,size(indata),nullval,indata,anyf,status)
    call ftclos(unit,status)
    data = indata(:,:,:,numargs+1:numargs+numsamples)
    deallocate(indata)

  end subroutine comm_read_chain

  ! **********************************************************************
  !                        Utilities
  ! **********************************************************************
  subroutine comm_bin_cls(cls_in, lmin, lmax, ind, cl_out)
    implicit none

    integer(i4b),                   intent(in)  :: lmin, lmax
    integer(i4b), dimension(1:),    intent(in)  :: ind
    real(dp),     dimension(0:,1:), intent(in)  :: cls_in
    real(dp),     dimension(1:,1:), intent(out) :: cl_out

    integer(i4b) :: i, j, l, n, p, q, nmaps

    p = size(ind)
    if (size(cls_in,2) == 1) then
       nmaps = 1
    else
       nmaps = 3
    end if

    cl_out = 0.d0
    do i = 1, p
       do j = i, p
          q = ind(i)*(1-ind(i))/2 + (ind(i)-1)*nmaps + ind(j)
          n = 0
          do l = lmin, lmax
             cl_out(i,j) = cl_out(i,j) + (2*l+1) * cls_in(l,q)
             n           = n           + (2*l+1)
          end do
          cl_out(i,j) = cl_out(i,j) / n
          cl_out(j,i) = cl_out(i,j)
       end do
    end do

  end subroutine comm_bin_cls

  function comm_det(M)
    implicit none

    real(dp), dimension(:,:), intent(in) :: M
    real(dp)                             :: comm_det

    integer(i4b) :: p

    p = size(M,1)

    if (p == 1) then
       comm_det = M(1,1)
    else if (p == 2) then
       comm_det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
    else
       write(*,*) 'Higher than 2D determinants not yet implemented'
       stop
    end if

  end function comm_det

  function comm_tr_A_invB(A, B)
    implicit none

    real(dp), dimension(:,:), intent(in) :: A, B
    real(dp)                             :: comm_tr_A_invB

    integer(i4b) :: p

    p = size(A,1)
    if (p == 1) then
       comm_tr_A_invB = A(1,1) / B(1,1)
    else if (p == 2) then
       comm_tr_A_invB = (A(1,1)*B(2,2) - A(1,2)*B(1,2) - A(2,1)*B(2,1) + A(2,2)*B(1,1)) / &
            & (B(1,1)*B(2,2) - B(1,2)*B(2,1))
    else 
       write(*,*) 'Higher than 2D traces not yet implemented'
       stop
    end if

  end function comm_tr_A_invB


  ! Small utility for converting an integer to a string
  subroutine comm_int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine comm_int2string

  subroutine cl2s(cl, S)
    implicit none

    real(dp), dimension(1:),    intent(in)  :: cl
    real(dp), dimension(1:,1:), intent(out) :: S

    integer(i4b) :: i, j, k, nmaps

    nmaps = size(S,1)

    S = 0.d0
    k = 1
    do i = 1, nmaps
       do j = i, nmaps
          S(i,j) = cl(k)
          S(j,i) = cl(k)
          k      = k+1
       end do
    end do

  end subroutine cl2s

  subroutine read_fiducial_spectrum(clfile, cls)
    implicit none

    character(len=*),                              intent(in)  :: clfile
    real(dp),         allocatable, dimension(:,:), intent(out) :: cls

    integer(i4b)        :: unit, l, lmax
    real(dp)            :: cl_in(4)
    character(len=2048) :: line

    ! Read fiducial spectrum file
    unit = comm_getlun()
    open(unit,file=trim(clfile))
    lmax = -1
    do while (.true.)
       read(unit,'(a)',end=89) line
       line = trim(line)
       if (line(1:1) == '#') cycle
       read(line,*) l
       lmax = max(lmax, l)
    end do
89  close(unit)
    if (lmax > -1) then
       allocate(cls(0:lmax,6))
       cls = 0.d0
       open(unit,file=trim(clfile))
       do while (.true.)
          read(unit,'(a)',end=90) line
          line = trim(line)
          if (line(1:1) == '#') cycle
          read(line,*) l, cl_in
          cls(l,1) = cl_in(1) ! Assume (TT, EE, BB, TE) ordering for now
          cls(l,2) = cl_in(4)          
          cls(l,4) = cl_in(2)
          cls(l,6) = cl_in(3)
       end do
90     close(unit)
    else
       write(*,*) 'Error -- comm_mod: No valid entries in clfile = ', trim(clfile)
       stop
    end if

  end subroutine read_fiducial_spectrum

  subroutine cholesky_decompose_with_mask_dp(M, ierr)
    implicit none

    real(dp),     dimension(1:,1:), intent(inout) :: M
    integer(i4b),                   intent(out), optional :: ierr

    integer(i4b)     :: i, j, q, n, info
    integer(i4b), allocatable, dimension(:)   :: indmap
    real(dp),     allocatable, dimension(:,:) :: L

    if (present(ierr)) ierr = 0

    q     = size(M,1)
    allocate(indmap(q))
    n = 0
    do i = 1, q
       if (M(i,i) > 0.d0) then
          n = n+1
          indmap(n) = i
       else if (M(i,i) < 0.d0) then
          M    = 0.d0
          ierr = 1
          deallocate(indmap)
          return
       end if
    end do
    if (n == 0) then
       M = 0.d0
       deallocate(indmap)
       return
    end if

    allocate(L(n,n))
    L = M(indmap(1:n),indmap(1:n))

    call dpotrf('U', n, L, n, info)
    if (info /= 0) then
       if (present(ierr)) then
          ierr = 1
       else
          write(*,*) 'DPOTRF: Cholesky factorization failed. Info = ', info
          stop
       end if
    end if

    ! Nullify lower triangle and missing rows/columns
    do i = 1, n
       L(i,1:i-1) = 0.d0
    end do

    ! Copy back to main matrix
    M                          = 0.d0
    M(indmap(1:n),indmap(1:n)) = L

    deallocate(indmap, L)

  end subroutine cholesky_decompose_with_mask_dp

  subroutine convert_real_to_complex_alms(alms_real, alms_cmplx)
    implicit none

    real(sp),     dimension(1:,1:),    intent(in)  :: alms_real
    complex(spc), dimension(1:,0:,0:), intent(out) :: alms_cmplx

    integer(i4b) :: ind, l, m, lmax
    real(dp)     :: one_over_sqrt_two 

    lmax              = size(alms_cmplx(1,:,0))-1
    one_over_sqrt_two = 1.d0 / sqrt(2.d0)

    alms_cmplx = cmplx(0.d0, 0.d0)
    do l = 0, lmax
       ind = l**2 + l + 1

       alms_cmplx(:,l,0) = cmplx(alms_real(ind,:),0.d0)
       do m = 1, l
          alms_cmplx(:,l,m) = one_over_sqrt_two * cmplx(alms_real(ind+m,:), alms_real(ind-m,:))
       end do
    end do

  end subroutine convert_real_to_complex_alms

  subroutine convert_real_to_complex_alms_dp(alms_real, alms_cmplx)
    implicit none

    real(dp),     dimension(1:,1:),    intent(in)  :: alms_real
    complex(dpc), dimension(1:,0:,0:), intent(out) :: alms_cmplx

    integer(i4b) :: ind, l, m, lmax
    real(dp)     :: one_over_sqrt_two 

    lmax              = size(alms_cmplx(1,:,0))-1
    one_over_sqrt_two = 1.d0 / sqrt(2.d0)

    alms_cmplx = cmplx(0.d0, 0.d0)
    do l = 0, lmax
       ind = l**2 + l + 1

       alms_cmplx(:,l,0) = cmplx(alms_real(ind,:),0.d0)
       do m = 1, l
          alms_cmplx(:,l,m) = one_over_sqrt_two * cmplx(alms_real(ind+m,:), alms_real(ind-m,:))
       end do
    end do

  end subroutine convert_real_to_complex_alms_dp

  subroutine convert_complex_to_real_alms(alms_cmplx, alms_real)
    implicit none

    complex(dpc), dimension(1:,0:,0:), intent(in)  :: alms_cmplx
    real(dp),     dimension(1:,1:),    intent(out) :: alms_real

    integer(i4b) :: ind, l, m, lmax
    real(dp)     :: sqrt_two

    lmax     = size(alms_cmplx(1,:,0))-1
    sqrt_two = sqrt(2.d0)

    alms_real = cmplx(0.d0, 0.d0)
    do l = 0, lmax
       ind = l**2 + l + 1

       alms_real(ind,:) = real(alms_cmplx(:,l,0),dp)
       do m = 1, l
          alms_real(ind-m,:) = sqrt_two * imag(alms_cmplx(:,l,m))
          alms_real(ind+m,:) = sqrt_two * real(alms_cmplx(:,l,m))
       end do
    end do

  end subroutine convert_complex_to_real_alms

  function comm_cls_posdef(cls)
    implicit none

    real(dp),     dimension(:,:), intent(in) :: cls
    logical(lgt)                             :: comm_cls_posdef

    integer(i4b) :: i, nmaps
    logical(lgt) :: ok
    real(dp), allocatable, dimension(:)   :: W
    real(dp), allocatable, dimension(:,:) :: S

    nmaps = size(cls,2)
    comm_cls_posdef = .true.
    do i = 1, size(cls,1)
       if (nmaps == 1) then
          ! Only temperature
          ok = cls(i,1) >= 0.d0
       else if (cls(i,3) == 0.d0 .and. cls(i,5) == 0.d0) then
          ! TB = EB = 0 => 2x2 + 1
          ok = cls(i,1)*cls(i,4)-cls(i,2)**2 >= 0.d0 .and. cls(i,6) >= 0.d0
       else
          ! All six spectra are included
          allocate(S(nmaps,nmaps), W(nmaps))
          call cl2S(cls(i,:), S)
          call comm_get_eigenvalues(S, W)
          ok = all(W > 0.d0)
          deallocate(S, W)
       end if
       if (.not. ok) then
          comm_cls_posdef = .false.
          return
       end if
    end do

  end function comm_cls_posdef

  subroutine comm_get_eigenvalues(A, eigenvals)
    
    real(dp), dimension(1:,1:), intent(in)            :: A
    real(dp), dimension(1:),    intent(out)           :: eigenvals

    integer(i4b)     :: n, liwork, lwork, lda, info
    character(len=1) :: job, uplo
    real(dp),     allocatable, dimension(:,:) :: A_copy
    real(dp),     allocatable, dimension(:)   :: W, work
    integer(i4b), allocatable, dimension(:)   :: iwork    

    n      = size(eigenvals)

    if (n == 1) then

       eigenvals(1) = A(1,1)

    else if (n == 2) then

       eigenvals(1) = 0.5d0*(A(1,1)+A(2,2) + sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))
       eigenvals(2) = 0.5d0*(A(1,1)+A(2,2) - sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))

    else

       job    = 'n'
       uplo   = 'l'
       lda    = n
       liwork = 1
       lwork  = 2*n+1

       ! Perform eigenvalue decomposition
       allocate(work(lwork))
       allocate(iwork(liwork))
       allocate(A_copy(n,n))
       A_copy = A
       call dsyevd(job, uplo, n, A_copy, lda, eigenvals, work, lwork, iwork, liwork, info)
       if (info /= 0) write(*,*) 'get_eigenvalues -- dsyevd info = ', info
       
       deallocate(work)
       deallocate(iwork)
       deallocate(A_copy)

    end if

  end subroutine comm_get_eigenvalues


  ! Convention: First psi around z, then theta around y, then phi around z
  subroutine compute_euler_matrix_zyz(phi, theta, psi, euler_matrix)
    implicit none

    real(dp),                 intent(in)  :: phi, theta, psi
    real(dp), dimension(3,3), intent(out) :: euler_matrix

    real(dp) :: sphi, cphi, sth, cth, spsi, cpsi

    sphi = sin(phi)
    cphi = cos(phi)

    sth  = sin(theta)
    cth  = cos(theta)

    spsi = sin(psi)
    cpsi = cos(psi)

    euler_matrix(1,1) = -sphi * spsi + cth * cphi * cpsi
    euler_matrix(1,2) = -sphi * cpsi - cth * cphi * spsi
    euler_matrix(1,3) =                sth * cphi
    euler_matrix(2,1) =  cphi * spsi + cth * sphi * cpsi
    euler_matrix(2,2) =  cphi * cpsi - cth * sphi * spsi
    euler_matrix(2,3) =                sth * sphi
    euler_matrix(3,1) =              - sth * cpsi
    euler_matrix(3,2) =                sth * spsi
    euler_matrix(3,3) =                cth

  end subroutine compute_euler_matrix_zyz

   subroutine get_basis_S(BBamp, power, cls, beam, Y, sqrtS)
    implicit none

    real(dp),                   intent(in)  :: BBamp, power
    real(dp), dimension(0:,1:), intent(in)  :: cls, beam
    real(dp), dimension(1:,1:), intent(in)  :: Y
    real(dp), dimension(1:,1:), intent(out) :: sqrtS

    real(dp) :: var, sqrt_sqrt_Cl
    integer(i4b) :: i, j, k, l, m, lmax, n_p, nmaps, n_h, numcomp
    real(dp), allocatable, dimension(:,:) :: buffer

    lmax    = size(beam,1)-1
    nmaps   = size(beam,2)
    n_p     = size(Y,1)
    n_h     = size(Y,2)
    numcomp = n_h / nmaps

    allocate(buffer(n_p, n_h))
    
    k = 1
    do i = 1, nmaps
       do l = 0, lmax
          if (l <= 1) then
             sqrt_sqrt_Cl = 0.d0
          else
             if (i == 1) then
                sqrt_sqrt_Cl = sqrt((cls(l,1) * 2.d0*pi/real(l*(l+1),dp)*beam(l,1)**2)**power)
             else if (i == 2) then
                sqrt_sqrt_Cl = sqrt((cls(l,4) * 2.d0*pi/real(l*(l+1),dp)*beam(l,2)**2)**power)
             else if (i == 3) then
                sqrt_sqrt_Cl = sqrt((BBamp*cls(l,4) * 2.d0*pi/real(l*(l+1),dp)*beam(l,2)**2)**power)
             end if
          end if
          do m = -l, l
             buffer(:,k) = Y(:,k) * sqrt_sqrt_Cl 
             k           = k+1
          end do
       end do
    end do
    call dsyrk('L','N',n_p,n_h,1.d0,buffer,n_p,0.d0,sqrtS,n_p)

    ! Symmetrize matrix
    do i = 1, n_p
       do j = i+1, n_p
          sqrtS(i,j) = sqrtS(j,i)
       end do
    end do

    deallocate(buffer)

  end subroutine get_basis_S

  subroutine compute_asymmetric_errors(x, y, peak, upper, lower)
    implicit none

    real(dp), dimension(:), intent(in)  :: x, y
    real(dp),               intent(out) :: peak, upper, lower

    integer(i4b) :: i, n, posmax(1), left, right
    real(dp)     :: integral, totint

    n        = size(x)
    posmax   = maxloc(y)
    left     = posmax(1)
    right    = posmax(1)
    integral = 0.d0
    totint   = comm_lowl_tsum(x,y)

    do while (integral < 0.68d0*totint)
       if (left == 1) then
          right = right+1
       else if (right == n) then
          left  = left-1
       else
          if (y(left) > y(right)) then
             left = left-1
          else
             right = right+1
          end if
       end if
       integral = comm_lowl_tsum(x(left:right), y(left:right))
    end do

    peak  = x(posmax(1))
    upper = x(right)-peak
    lower = peak-x(left)
!    write(*,*) x(left), x(posmax(1)), x(right)
!    write(*,*) peak, upper, lower

  end subroutine compute_asymmetric_errors

  subroutine remove_duplicates(x, y, n, epsilon)
    implicit none

    real(dp),     dimension(:),   intent(inout) :: x
    real(dp),     dimension(:,:), intent(inout) :: y
    integer(i4b),                 intent(out)   :: n
    real(dp),                     intent(in)    :: epsilon

    integer(i4b) :: i, m

    n = size(x)
    i = 2
    do while (i <= n)
       if (abs(x(i)-x(i-1)) < epsilon*abs(x(i))) then
          write(*,*) i, x(i-1), x(i)
          x(i:n-1)   = x(i+1:n)
          y(i:n-1,:) = y(i+1:n,:)
          n          = n-1
       else
          i          = i+1
       end if
    end do

  end subroutine remove_duplicates

  subroutine insert_entry(x_new, y_new, x, y, n, epsilon)
    implicit none

    real(dp),                 intent(in)    :: x_new, epsilon
    real(dp), dimension(:),   intent(in)    :: y_new
    real(dp), dimension(:),   intent(inout) :: x
    real(dp), dimension(:,:), intent(inout) :: y
    integer(i4b),             intent(inout) :: n

    integer(i4b) :: i
    logical(lgt) :: accept

    if (n == size(x)) then
       write(*,*) 'Error in insert_entry -- array already full'
    end if

    if (n == 0) then
       n      = 1
       x(n)   = x_new
       y(n,:) = y_new
    else
       i = 1
       do while (x_new > x(i) .and. i <= n)
          i = i+1
       end do
       
       if (i == 1) then
          accept = abs(x_new-x(1)) > epsilon*abs(x(1))
       else if (i == n) then
          accept = abs(x_new-x(n)) > epsilon*abs(x(n))
       else
          accept = abs(x_new-x(i)) > epsilon*abs(x(i)) .and. abs(x_new-x(i-1)) > epsilon*abs(x(i-1))
       end if

       if (accept) then
          x(i+1:n+1)   = x(i:n)
          y(i+1:n+1,:) = y(i:n,:)
          x(i)         = x_new
          y(i,:)       = y_new
          n            = n+1
       end if
    end if

  end subroutine insert_entry

  function comm_lowl_tsum(x, y)
    implicit none

    real(dp), dimension(:), intent(in) :: x, y
    real(dp)                           :: comm_lowl_tsum

    integer(i4b) :: i

    comm_lowl_tsum = 0.d0
    do i = 1, size(x)-1
       comm_lowl_tsum = comm_lowl_tsum + 0.5d0 * (y(i)+y(i+1)) * (x(i+1)-x(i))
    end do

  end function comm_lowl_tsum

  subroutine write_lowl_datafile(filename, d, C, beam, P_harm, mapfile, maskfile, &
       & covfile, beamfile, clfile, basis, lcut_T, lcut_P, epsilon_T, epsilon_P)
    implicit none

    character(len=*),              intent(in) :: filename
    real(dp), dimension(1:,1:),    intent(in) :: d, C, P_harm
    real(dp), dimension(0:,1:),    intent(in) :: beam
    character(len=*),              intent(in) :: mapfile, maskfile, covfile, beamfile, clfile
    character(len=*),              intent(in) :: basis
    integer(i4b),                  intent(in) :: lcut_T, lcut_P
    real(dp),                      intent(in) :: epsilon_T, epsilon_P
    
    integer(i4b)         :: lmax, nmaps, n_d, n_h, n, m, unit
    integer(i4b)         :: status, blocksize, readwrite
    integer(i4b)         :: nelements, fpixel, group, bitpix, naxis, hdutype
    logical(lgt)         :: simple, extend, anyf, exist
    real(dp)             :: nullval
    character(len=80)    :: comment, errorline
    integer(i4b), dimension(2) :: naxes

    unit  = comm_getlun()
    lmax  = size(beam,1)-1
    nmaps = size(beam,2)
    n     = size(d,1)       ! Number of modes
    n_d   = size(d,2)       ! Number of independent data sets
    n_h   = size(P_harm,1)  ! Number of harmonics
    m     = len(trim(filename))

    write(*,*) 'Writing data objects to ', trim(filename)
    if (filename(m-2:m) == 'unf') then

       open(unit,file=trim(filename),form='unformatted')
       write(unit) lmax, nmaps
       write(unit) n_d, n_h, n
       write(unit) d
       write(unit) C
       write(unit) beam
       write(unit) P_harm
       close(unit)

    else if (filename(m-3:m) == 'fits') then

       ! Open a new fitsfile
       status    = 0
       bitpix    = -64
       simple    = .true.
       extend    = .true.
       blocksize = 1
       group     = 1
       fpixel    = 1
       call ftinit(unit,trim(filename),blocksize,status)

       ! -----------------------------------------
       ! Write data vector
       ! -----------------------------------------
       naxis     = 2
       naxes(1)  = n
       naxes(2)  = n_d
       nelements = naxes(1)*naxes(2)
       call ftcrhd(unit, status)
       call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
       call ftpkys(unit,'BASIS',    basis, 'Basis type', status)
       call ftpkys(unit,'MAPFILE',  mapfile, 'Map filename', status)
       call ftpkys(unit,'COVFILE',  covfile, 'Covariance matrix filename', status)
       call ftpkys(unit,'BEAMFILE', beamfile, 'Beam profile filename', status)
       call ftpkys(unit,'CLFILE',   clfile, 'Fiducial C_l filename', status)
       call ftpkyj(unit,'LMAX',  lmax, 'Maximum multipole moment', status)
       call ftpkyj(unit,'NMAPS', nmaps, 'Maximum multipole moment', status)
       call ftpkyj(unit,'N_D',   n_d,  'Number of data vectors',   status)
       call ftpkyj(unit,'N',     n,    'Number of basis vectors',  status)
       call ftpkyj(unit,'N_H',   n_h,  'Number of harmonics',      status)
       call ftpkyj(unit,'L_T',   lcut_T,  'Truncation multipole (temperature)',  status)
       call ftpkyj(unit,'L_P',   lcut_P,  'Truncation multipole (polarization)', status)
       call ftpkyd(unit,'EPS_T', epsilon_T,  0, 'Eigenvalue threshold (temperature)', status)
       call ftpkyd(unit,'EPS_P', epsilon_P,  0, 'Eigenvalue threshold (polarization)', status)
       call ftpcom(unit, "------------------------------------", status)
       call ftpcom(unit, "This file is generated with comm_like_tools,", status)
       call ftpcom(unit, "available in the Commander SVN repository at", status)
       call ftpcom(unit, "https://codeforge.lbl.gov/svn/planck/trunk/commander", status)
       call ftpcom(unit, "Contact Hans Kristian Eriksen (h.k.k.eriksen@astro.uio.no)", status)
       call ftpcom(unit, "if you want access to this repository", status)
       call ftpcom(unit, "------------------------------------", status)
       call ftpprd(unit, group, fpixel, nelements, d, status)

       ! -----------------------------------------
       ! Write covariance matrix
       ! -----------------------------------------
       naxis     = 2
       naxes(1)  = n
       naxes(2)  = n
       nelements = naxes(1)*naxes(2)
       call ftcrhd(unit, status)
       call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
       call ftpprd(unit, group, fpixel, nelements, C, status)

       ! -----------------------------------------
       ! Write beam
       ! -----------------------------------------
       naxis     = 2
       naxes(1)  = lmax+1
       naxes(2)  = nmaps
       nelements = naxes(1)*naxes(2)
       call ftcrhd(unit, status)
       call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
       call ftpprd(unit, group, fpixel, nelements, beam, status)

       ! -----------------------------------------
       ! Write P_harm
       ! -----------------------------------------
       naxis     = 2
       naxes(1)  = n_h
       naxes(2)  = n
       nelements = naxes(1)*naxes(2)
       call ftcrhd(unit, status)
       call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
       call ftpprd(unit, group, fpixel, nelements, P_harm, status)
       
       ! Write the current date
       call ftpdat(unit,status)
       
       ! close the file and free the unit number
       call ftclos(unit, status)
       call ftfiou(unit, status)
       
       if (status /= 0) then
          write(*,*) 'Error in FITS operations -- status = ', status
          call ftgmsg(errorline)
          do while (trim(errorline) /= '')
             write(*,*) errorline
             call ftgmsg(errorline)
          end do
       end if

    else
       write(*,*) 'Error -- unknown filetype = ', trim(filename)
       stop
    end if

  end subroutine write_lowl_datafile

  subroutine read_lowl_datafile(filename, n, n_d, n_h, lmax, nmaps, d, C, beam, P_harm, window)
    implicit none

    character(len=*),                        intent(in)  :: filename
    integer(i4b),                            intent(out) :: n, n_d, n_h, lmax, nmaps
    real(dp), allocatable, dimension(:,:),   intent(out) :: d, C, P_harm
    real(dp), allocatable, dimension(:,:),   intent(out) :: beam
    real(dp), allocatable, dimension(:,:,:), intent(out) :: window
    
    integer(i4b)         :: i, j, l, m, unit
    integer(i4b)         :: status, blocksize, readwrite, hdutype
    integer(i4b)         :: nelements, fpixel, group, bitpix
    logical(lgt)         :: simple, extend, anyf, exist
    real(dp)             :: nullval
    character(len=80)    :: comment, errorline

    unit  = comm_getlun()
    m     = len(trim(filename))

    if (filename(m-2:m) == 'unf') then

       open(unit,file=trim(filename),form='unformatted')
       read(unit) lmax, nmaps
       read(unit) n_d, n_h, n
       allocate(window(0:lmax,nmaps,nmaps), beam(0:lmax,nmaps))
       allocate(d(n,n_d), C(n,n), P_harm(n_h,n))
       read(unit) d
       read(unit) C
       read(unit) beam
       read(unit) P_harm
       close(unit)

    else if (filename(m-3:m) == 'fits') then

       status    = 0
       readwrite = 0
       nullval   = 0.d0
       group     = 1
       fpixel    = 1

       ! Open the data file
       call ftopen(unit, trim(filename), readwrite, blocksize, status)

       ! Read data vector and keywords
       call ftmahd(unit, 1, hdutype, status)
       call ftgkyj(unit, 'LMAX',  lmax, comment,status)
       call ftgkyj(unit, 'NMAPS', nmaps, comment,status)
       call ftgkyj(unit, 'N_D',   n_d,  comment,status)
       call ftgkyj(unit, 'N',     n,    comment,status)
       call ftgkyj(unit, 'N_H',   n_h,  comment,status)
       allocate(window(0:lmax,nmaps,nmaps), beam(0:lmax,nmaps))
       allocate(d(n,n_d), C(n,n), P_harm(n_h,n))
       call ftgpvd(unit, group, fpixel, size(d), nullval, d, anyf, status)

       ! Read covariance matrix
       call ftmahd(unit, 2, hdutype, status)
       call ftgpvd(unit, group, fpixel, size(C), nullval, C, anyf, status)

       ! Read beam
       call ftmahd(unit, 3, hdutype, status)
       call ftgpvd(unit, group, fpixel, size(beam), nullval, beam, anyf, status)

       ! Read projection operator
       call ftmahd(unit, 4, hdutype, status)
       call ftgpvd(unit, group, fpixel, size(P_harm), nullval, P_harm, anyf, status)

       ! Close file
       call ftclos(unit, status)
       
       if (status /= 0) then
          write(*,*) 'Error in FITS operations -- status = ', status
          call ftgmsg(errorline)
          do while (trim(errorline) /= '')
             write(*,*) errorline
             call ftgmsg(errorline)
          end do
       end if

    else
       write(*,*) 'Error -- unknown filetype = ', trim(filename)
       stop
    end if

    do l = 0, lmax
       do i = 1, nmaps
          do j = 1, nmaps
             window(l,i,j) = beam(l,i)*beam(l,j)
          end do
       end do
    end do
    
  end subroutine read_lowl_datafile


  subroutine write_gauss_BR_datafile(filename, lmin, lmax, cl2x, mu, sigma, mu_sigma)
    implicit none

    character(len=*),                 intent(in) :: filename
    integer(i4b),                     intent(in) :: lmin, lmax
    real(dp), dimension(lmin:,1:,1:), intent(in) :: cl2x
    real(dp), dimension(lmin:),       intent(in) :: mu, mu_sigma
    real(dp), dimension(lmin:,lmin:), intent(in) :: sigma
    
    integer(i4b)         :: nmaps, n_d, n_h, n, m, unit, nbin
    integer(i4b)         :: status, blocksize, readwrite
    integer(i4b)         :: nelements, fpixel, group, bitpix, naxis, hdutype
    logical(lgt)         :: simple, extend, anyf, exist
    real(dp)             :: nullval
    character(len=80)    :: comment, errorline
    integer(i4b), dimension(3) :: naxes

    unit  = comm_getlun()

    write(*,*) 'Writing data objects to ', trim(filename)

    ! Open a new fitsfile
    status    = 0
    bitpix    = -64
    simple    = .true.
    extend    = .true.
    blocksize = 1
    group     = 1
    fpixel    = 1
    call ftinit(unit,trim(filename),blocksize,status)

    ! -----------------------------------------
    ! Write conversion table
    ! -----------------------------------------
    nbin      = size(cl2x,1)
    naxis     = 3
    naxes(1)  = nbin
    naxes(2)  = lmax-lmin+1
    naxes(3)  = 3
    nelements = naxes(1)*naxes(2)*naxes(3)
    call ftcrhd(unit, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
    call ftpkyj(unit,'LMIN',    lmin, 'Lmin', status)
    call ftpkyj(unit,'LMAX',    lmax, 'Lmax', status)
    call ftpkyj(unit,'NBIN',    nbin, 'Nbin', status)
    call ftpcom(unit, "------------------------------------", status)
    call ftpcom(unit, "This file is generated with comm_like_tools,", status)
    call ftpcom(unit, "available in the Commander SVN repository at", status)
    call ftpcom(unit, "https://codeforge.lbl.gov/svn/planck/trunk/commander", status)
    call ftpcom(unit, "Contact Hans Kristian Eriksen (h.k.k.eriksen@astro.uio.no)", status)
    call ftpcom(unit, "if you want access to this repository", status)
    call ftpcom(unit, "------------------------------------", status)
    call ftpprd(unit, group, fpixel, nelements, cl2x, status)

    ! -----------------------------------------
    ! Write mean vector
    ! -----------------------------------------
    naxis     = 1
    naxes(1)  = lmax-lmin+1
    nelements = naxes(1)
    call ftcrhd(unit, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes(1:1), 0, 1, extend, status)
    call ftpprd(unit, group, fpixel, nelements, mu, status)

    ! -----------------------------------------
    ! Write covariance matrix
    ! -----------------------------------------
    naxis     = 2
    naxes(1)  = lmax-lmin+1
    naxes(2)  = lmax-lmin+1
    nelements = naxes(1)*naxes(2)
    call ftcrhd(unit, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes(1:2), 0, 1, extend, status)
    call ftpprd(unit, group, fpixel, nelements, sigma, status)

    ! -----------------------------------------
    ! Write mean sigma vector
    ! -----------------------------------------
    naxis     = 1
    naxes(1)  = lmax-lmin+1
    nelements = naxes(1)
    call ftcrhd(unit, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes(1:1), 0, 1, extend, status)
    call ftpprd(unit, group, fpixel, nelements, mu_sigma, status)

    ! Write the current date
    call ftpdat(unit,status)
    
    ! close the file and free the unit number
    call ftclos(unit, status)
    call ftfiou(unit, status)
    
    if (status /= 0) then
       write(*,*) 'Error in FITS operations -- status = ', status
       call ftgmsg(errorline)
       do while (trim(errorline) /= '')
          write(*,*) errorline
          call ftgmsg(errorline)
       end do 
    end if
    
  end subroutine write_gauss_BR_datafile

  subroutine read_gauss_BR_datafile(filename, lmin, lmax, mu, cov, cl2x, mu_sigma)
    implicit none

    character(len=*),                        intent(in)  :: filename
    integer(i4b),                            intent(in)  :: lmin, lmax
    real(dp), allocatable, dimension(:),     intent(out) :: mu, mu_sigma
    real(dp), allocatable, dimension(:,:),   intent(out) :: cov
    real(dp), allocatable, dimension(:,:,:), intent(out) :: cl2x
    
    integer(i4b)         :: i, j, l, m, unit
    integer(i4b)         :: status, blocksize, readwrite, hdutype, nbin
    integer(i4b)         :: nelements, fpixel, group, bitpix, lmax_in, lmin_in
    logical(lgt)         :: simple, extend, anyf, exist
    real(dp)             :: nullval
    character(len=80)    :: comment, errorline
    real(dp), allocatable, dimension(:)     :: mu_in, sigma_in
    real(dp), allocatable, dimension(:,:)   :: cov_in
    real(dp), allocatable, dimension(:,:,:) :: cl2x_in

    unit  = comm_getlun()
    m     = len(trim(filename))

    status    = 0
    readwrite = 0
    nullval   = 0.d0
    group     = 1
    fpixel    = 1

    ! Open the data file
    call ftopen(unit, trim(filename), readwrite, blocksize, status)

    ! Read data vector and keywords
    call ftmahd(unit, 1, hdutype, status)
    call ftgkyj(unit, 'LMIN',  lmin_in, comment,status)
    call ftgkyj(unit, 'LMAX',  lmax_in, comment,status)
    call ftgkyj(unit, 'NBIN',  nbin,    comment,status)
    allocate(cl2x_in(nbin,lmin_in:lmax_in,3), mu_in(lmin_in:lmax_in), sigma_in(lmin_in:lmax_in))
    call ftgpvd(unit, group, fpixel, size(cl2x_in), nullval, cl2x_in, anyf, status)

    ! Read mean vector
    call ftmahd(unit, 2, hdutype, status)
    call ftgpvd(unit, group, fpixel, size(mu_in), nullval, mu_in, anyf, status)

    ! Read covariance matrix
    call ftmahd(unit, 3, hdutype, status)
    call ftgpvd(unit, group, fpixel, size(cov_in), nullval, cov_in, anyf, status)

    ! Read mean vector
    call ftmahd(unit, 4, hdutype, status)
    call ftgpvd(unit, group, fpixel, size(sigma_in), nullval, sigma_in, anyf, status)

    ! Close file
    call ftclos(unit, status)
       
    if (status /= 0) then
       write(*,*) 'Error in FITS operations -- status = ', status
       call ftgmsg(errorline)
       do while (trim(errorline) /= '')
          write(*,*) errorline
          call ftgmsg(errorline)
       end do
    end if

    if (lmax > lmax_in) then
       write(*,*) 'Error: Requested lmax greater than that supported by ', trim(filename)
       stop
    end if

    if (lmin > lmin_in) then
       write(*,*) 'Error: Requested lmin smaller than that supported by ', trim(filename)
       stop
    end if

    ! Copy relevant data into output data structures
    allocate(cl2x(lmin:lmax,nbin,3), mu(lmin:lmax), cov(lmin:lmax,lmin:lmax), mu_sigma(lmin:lmax))
    cl2x     = cl2x_in(lmin:lmax,:,:)
    mu       = mu_in(lmin:lmax)
    cov      = cov_in(lmin:lmax,lmin:lmax)
    mu_sigma = sigma_in(lmin:lmax)

    deallocate(cl2x_in, mu_in, cov_in)
    
  end subroutine read_gauss_BR_datafile

end module comm_like_utils
