module comm_ptsrc_comp_mod
  use math_tools
  use comm_param_mod
  use comm_comp_mod
  use comm_F_int_mod
  use comm_F_int_0D_mod
  use comm_F_int_2D_mod
  use comm_data_mod
  use pix_tools
  use comm_hdf_mod
  use comm_cr_utils
  use comm_cr_precond_mod
  use locate_mod
  implicit none

  private
  public comm_ptsrc_comp, initPtsrcPrecond, updatePtsrcPrecond, applyPtsrcPrecond

  !**************************************************
  !            Compact object class
  !**************************************************
  type Tnu
     integer(i4b) :: nside, np, nmaps
     integer(i4b), allocatable, dimension(:,:) :: pix     ! Pixel list, both absolute and relative
     real(dp),     allocatable, dimension(:,:) :: map     ! (0:np-1,nmaps)
     real(dp),     allocatable, dimension(:)   :: F       ! Mixing matrix (nmaps)
     real(dp),     allocatable, dimension(:)   :: Omega_b ! Solid angle (nmaps
  end type Tnu

  type ptsrc
     character(len=512) :: outprefix
     character(len=512) :: id_src
     real(dp)           :: glon, glat, f_beam
     type(Tnu), allocatable, dimension(:)   :: T      ! Spatial template (nband)
     real(dp),  allocatable, dimension(:,:) :: theta  ! Spectral parameters (npar,nmaps)
  end type ptsrc
  
  type, extends (comm_comp) :: comm_ptsrc_comp
     character(len=512) :: outprefix
     real(dp)           :: cg_scale
     integer(i4b)       :: nside, nsrc
     real(dp),        allocatable, dimension(:,:) :: x      ! Amplitudes (sum(nsrc),nmaps)
     type(F_int_ptr), allocatable, dimension(:)   :: F_int  ! SED integrator (numband)
     type(ptsrc),     allocatable, dimension(:)   :: src    ! Source template (nsrc)
   contains
     procedure :: dumpFITS => dumpPtsrcToFITS
     procedure :: getBand  => evalPtsrcBand
     procedure :: projectBand  => projectPtsrcBand
     procedure :: updateF
     procedure :: S => evalSED
     procedure :: getScale
  end type comm_ptsrc_comp

  interface comm_ptsrc_comp
     procedure constructor
  end interface comm_ptsrc_comp

  type ptsrc_ptr
     class(comm_ptsrc_comp), pointer :: p
  end type ptsrc_ptr
  
  integer(i4b) :: ncomp_pre =   0
  integer(i4b) :: npre      =   0
  integer(i4b) :: nmaps_pre =  -1
  class(ptsrc_ptr), allocatable, dimension(:) :: ptsrcComps
  
contains

  function constructor(cpar, id, id_abs)
    implicit none
    class(comm_params),       intent(in) :: cpar
    integer(i4b),             intent(in) :: id, id_abs
    class(comm_ptsrc_comp),   pointer    :: constructor

    integer(i4b) :: i, j, k, nlist, npix, listpix(0:10000-1), hits(10000)
    real(dp)     :: vec0(3), vec(3), r
    
    ! General parameters
    allocate(constructor)

    ! Initialize general parameters
    constructor%class     = cpar%cs_class(id_abs)
    constructor%type      = cpar%cs_type(id_abs)
    constructor%label     = cpar%cs_label(id_abs)
    constructor%id        = id
    constructor%nmaps     = 1; if (cpar%cs_polarization(id_abs)) constructor%nmaps = 3
    constructor%nu_ref    = cpar%cs_nu_ref(id_abs)
    constructor%nside     = cpar%cs_nside(id_abs)
    constructor%outprefix = trim(cpar%cs_label(id_abs))
    constructor%cg_scale  = cpar%cs_cg_scale(id_abs)
    allocate(constructor%poltype(1))
    constructor%poltype   = cpar%cs_poltype(1,id_abs)
    ncomp_pre                 = ncomp_pre + 1

    ! Initialize frequency scaling parameters
    allocate(constructor%F_int(numband))
    select case (trim(constructor%type))
    case ("radio")
       constructor%npar = 2   ! (alpha, beta)
       allocate(constructor%p_uni(2,constructor%npar), constructor%p_gauss(2,constructor%npar))
       allocate(constructor%theta_def(constructor%npar))
       constructor%p_uni     = cpar%cs_p_uni(id_abs,:,:)
       constructor%p_gauss   = cpar%cs_p_gauss(id_abs,:,:)
       constructor%theta_def = cpar%cs_theta_def(1:2,id_abs)
       do i = 1, numband
          constructor%F_int(i)%p => comm_F_int_2D(constructor, data(i)%bp)
       end do
    case ("fir")
       constructor%npar = 2   ! (beta, T_d)
       allocate(constructor%p_uni(2,constructor%npar), constructor%p_gauss(2,constructor%npar))
       allocate(constructor%theta_def(constructor%npar))
       constructor%p_uni     = cpar%cs_p_uni(id_abs,:,:)
       constructor%p_gauss   = cpar%cs_p_gauss(id_abs,:,:)
       constructor%theta_def = cpar%cs_theta_def(1:2,id_abs)
       do i = 1, numband
          constructor%F_int(i)%p => comm_F_int_2D(constructor, data(i)%bp)
       end do
       allocate(constructor%p_uni(2,constructor%npar), constructor%p_gauss(2,constructor%npar))
    case ("sz")
       constructor%npar = 0   ! (none)
       do i = 1, numband
          constructor%F_int(i)%p => comm_F_int_0D(constructor, data(i)%bp)
       end do
    case default
       call report_error("Unknown point source model: " // trim(constructor%type))
    end select

    ! Read and allocate source structures
    call read_sources(constructor, cpar, id, id_abs)

    ! Update mixing matrix
    call constructor%updateF

  end function constructor



  subroutine updateF(self, beta)
    implicit none
    class(comm_ptsrc_comp),                   intent(inout)        :: self
    real(dp),               dimension(:,:,:), intent(in), optional :: beta  ! (npar,nmaps,nsrc)

    integer(i4b) :: i, j
    
    do j = 1, self%nsrc
       if (present(beta)) then
          self%src(j)%theta = beta(:,:,j)
       else
          do i = 1, self%nmaps
             self%src(j)%theta(:,i) = self%theta_def
          end do
       end if
       do i = 1, numband
          ! Temperature
          self%src(j)%T(i)%F(1) = &
               & self%F_int(i)%p%eval(self%src(j)%theta(:,1)) * data(i)%gain * self%cg_scale

          ! Polarization
          if (self%nmaps == 3) then
             ! Stokes Q
             if (self%poltype(1) < 2) then
                self%src(j)%T(i)%F(2) = self%src(j)%T(i)%F(1)
             else
                self%src(j)%T(i)%F(2) = &
                     & self%F_int(i)%p%eval(self%src(j)%theta(:,2)) * data(i)%gain * self%cg_scale
             end if
          
             ! Stokes U
             if (self%poltype(1) < 3) then
                self%src(j)%T(i)%F(3) = self%src(j)%T(i)%F(2)
             else
                self%src(j)%T(i)%F(3) = &
                     & self%F_int(i)%p%eval(self%src(j)%theta(:,3)) * data(i)%gain * self%cg_scale
             end if
          end if
       end do
    end do
    
  end subroutine updateF

  function evalSED(self, nu, band, theta)
    class(comm_ptsrc_comp),    intent(in)           :: self
    real(dp),                  intent(in), optional :: nu
    integer(i4b),              intent(in), optional :: band
    real(dp), dimension(1:),   intent(in), optional :: theta
    real(dp)                                        :: evalSED

    real(dp) :: x
    
    select case (trim(self%type))
    case ("radio")
       !evalSED = exp(theta(1) * (nu/self%nu_ref) + theta(2) * (log(nu/self%nu_ref))**2) * &
       !     & (self%nu_ref/nu)**2
       evalSED = (self%nu_ref/nu)**(2.d0+theta(1)) 
    case ("fir")
       x = h/(k_B*theta(2))
       evalSED = (exp(x*self%nu_ref)-1.d0)/(exp(x*nu)-1.d0) * (nu/self%nu_ref)**(theta(1)+1.d0)
    case ("sz")
       evalSED = 0.d0
       call report_error('SZ not implemented yet')
    end select
    
  end function evalSED

  function evalPtsrcBand(self, band, amp_in, pix, alm_out)
    implicit none
    class(comm_ptsrc_comp),                       intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    real(dp),        dimension(:,:), allocatable                        :: evalPtsrcBand

    integer(i4b) :: i, j, p, q
    real(dp)     :: amp

    if (.not. allocated(evalPtsrcBand)) &
         & allocate(evalPtsrcBand(0:data(band)%info%np-1,data(band)%info%nmaps))

    ! Loop over sources
    evalPtsrcBand = 0.d0
    do i = 1, self%nsrc
       do j = 1, self%src(i)%T(band)%nmaps
          if (present(amp_in)) then
             amp = amp_in(i,j)
          else
             amp = self%x(i,j)
          end if
          
          ! Scale to correct frequency through multiplication with mixing matrix
          amp = self%getScale(band,i,j) * self%src(i)%T(band)%F(j) *  amp
          
          ! Project with beam
          do q = 1, self%src(i)%T(band)%np
             p = self%src(i)%T(band)%pix(q,1)
             evalPtsrcBand(p,j) = evalPtsrcBand(p,j) + amp * self%src(i)%T(band)%map(q,j)
          end do
       end do
    end do
    
  end function evalPtsrcBand
  
  ! Return component projected from map
  function projectPtsrcBand(self, band, map, alm_in)
    implicit none
    class(comm_ptsrc_comp),                       intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    real(dp),        dimension(:,:), allocatable                        :: projectPtsrcBand

    integer(i4b) :: i, j, q, p
    real(dp)     :: val
    
    if (.not. allocated(projectPtsrcBand)) &
         & allocate(projectPtsrcBand(self%nsrc*self%nmaps,self%nmaps))

    ! Loop over sources
    projectPtsrcBand = 0.d0
    do i = 1, self%nsrc
       do j = 1, self%src(i)%T(band)%nmaps
          val = 0.d0
          do q = 1, self%src(i)%T(band)%np
             p   = self%src(i)%T(band)%pix(q,1)
             val = val + self%src(i)%T(band)%map(q,j) * map%map(p,j)
          end do

          ! Scale to correct frequency through multiplication with mixing matrix
          val = self%getScale(band,i,j) * self%src(i)%T(band)%F(j) * val

          ! Return value
          projectPtsrcBand(i,j) = val
       end do
    end do
    
  end function projectPtsrcBand
  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpPtsrcToFITS(self, postfix, dir)
    class(comm_ptsrc_comp),                  intent(in)           :: self
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i, l, m, ierr, unit
    real(dp)           :: vals(10)
    logical(lgt)       :: exist, first_call = .true.
    character(len=512) :: filename
    class(comm_map), pointer :: map

    ! Output point source maps for each frequency
    do i = 1, numband
       map => comm_map(data(i)%info)
       map%map = self%getBand(i) * self%cg_scale
       filename = trim(self%label) // '_' // trim(data(i)%label) // '_' // trim(postfix) // '.fits'
       call map%writeFITS(trim(dir)//'/'//trim(filename))
       deallocate(map)
    end do

  end subroutine dumpPtsrcToFITS


  subroutine read_sources(self, cpar, id, id_abs)
    implicit none
    class(comm_ptsrc_comp), intent(inout) :: self
    class(comm_params),     intent(in)    :: cpar
    integer(i4b),           intent(in)    :: id, id_abs

    integer(i4b)        :: unit, i, j, npar, nmaps, pix, nside, n
    real(dp)            :: glon, glat, nu_ref
    logical(lgt)        :: pol
    character(len=1024) :: line, filename
    character(len=128)  :: id_ptsrc, flabel
    real(dp), allocatable, dimension(:)   :: amp
    real(dp), allocatable, dimension(:,:) :: beta

    unit = getlun()

    nmaps = 1; if (cpar%cs_polarization(id_abs)) nmaps = 3
    select case (trim(cpar%cs_type(id_abs)))
    case ("radio")
       npar = 2
    case ("fir")
       npar = 2
    case ("sz")
       npar = 0
    end select
    allocate(amp(nmaps), beta(npar,nmaps))

    ! Count number of valid sources
    open(unit,file=trim(cpar%datadir) // '/' // trim(cpar%cs_catalog(id_abs)),recl=1024)
    self%nsrc = 0
    self%ncr  = 0
    do while (.true.)
       read(unit,'(a)',end=1) line
       line = trim(line)
       if (line(1:1) == '#' .or. trim(line) == '') then
          cycle
       else
          self%nsrc = self%nsrc + 1
          self%ncr  = self%ncr  + nmaps
          npre      = npre + 1
          nmaps_pre = max(nmaps_pre, nmaps)
       end if
    end do 
1   close(unit)
    
    ! Initialize point sources based on catalog information
    allocate(self%x(self%nsrc,self%nmaps), self%src(self%nsrc))
    open(unit,file=trim(cpar%datadir) // '/' // trim(cpar%cs_catalog(id_abs)),recl=1024)
    i = 0
    do while (.true.)
       read(unit,'(a)',end=2) line
       line = trim(line)
       if (line(1:1) == '#' .or. trim(line) == '') cycle
       read(line,*) flabel, pix, nside, glon, glat, amp, beta, id_ptsrc
       i                             = i+1
       allocate(self%src(i)%theta(self%npar,self%nmaps), self%src(i)%T(numband))
       self%src(i)%id_src   = id_ptsrc
       self%src(i)%glon     = glon * DEG2RAD
       self%src(i)%glat     = glat * DEG2RAD
       self%src(i)%theta    = beta
       self%x(i,:)          = amp / self%cg_scale
    end do 
2   close(unit)


    ! Initialize beam templates
    do j = 1, self%nsrc
       do i = 1, numband
          self%src(j)%T(i)%nside   = data(i)%info%nside
          self%src(j)%T(i)%nmaps   = min(data(i)%info%nmaps, self%nmaps)
          allocate(self%src(j)%T(i)%F(self%src(j)%T(i)%nmaps))
          self%src(j)%T(i)%F       = 0.d0

          ! Get pixel space template
          filename = trim(cpar%datadir)//'/'//trim(cpar%ds_btheta_file(i))
          n        = len(trim(adjustl(filename)))
          if (n == 0 .or. trim(filename) == 'none') then
             ! Build template internally from b_l
             call report_error('Bl ptsrc beam not yet implemented')
             call compute_symmetric_beam(self%src(j)%glon, self%src(j)%glat, &
                  & self%src(j)%T(i), bl=data(i)%B%b_l)
          else if (filename(n-2:n) == '.dat') then
             ! Build template internally from b_l
             call report_error('Radial ptsrc beam profile not yet implemented')
          else if (filename(n-2:n) == '.h5' .or. filename(n-3:n) == '.hdf') then
             ! Read precomputed Febecop beam from HDF file
             call read_febecop_beam(filename, self%src(j)%glon, self%src(j)%glat, i, &
                  & self%src(j)%T(i))
          else
             call report_error('Unsupported point source template = '//trim(filename))
          end if
       end do
    end do
    
  end subroutine read_sources


  subroutine read_febecop_beam(filename, glon, glat, band, T)
    implicit none
    character(len=*), intent(in)    :: filename
    real(dp),         intent(in)    :: glon, glat
    integer(i4b),     intent(in)    :: band
    type(Tnu),        intent(inout) :: T

    integer(i4b)      :: i, j, n, pix, ext(1)
    character(len=16) :: itext
    type(hdf_file)    :: file
    integer(i4b), allocatable, dimension(:)   :: ind
    integer(i4b), allocatable, dimension(:,:) :: mypix
    real(dp),     allocatable, dimension(:,:) :: b, mybeam

    ! Find center pixel number for current source
    call ang2pix_ring(T%nside, 0.5d0*pi-glat, glon, pix)

    ! Find number of pixels in beam
    write(itext,*) pix
    call open_hdf_file(filename, file, 'r')
    call get_size_hdf(file, trim(adjustl(itext))//'/indices', ext)
    n = ext(1)

    ! Read full beam from file
    allocate(ind(n), b(n,T%nmaps), mypix(n,2), mybeam(n,T%nmaps))
    call read_hdf(file, trim(adjustl(itext))//'/indices', ind)
    call read_hdf(file, trim(adjustl(itext))//'/values',  b(:,1))
    call close_hdf_file(file)

    ! Find number of pixels belonging to current processor
    T%np = 0
    i    = 1
    j    = locate(data(band)%info%pix, ind(i))
    if (j > -1) then
       do while (.true.)
          if (ind(i) == data(band)%info%pix(j)) then
             T%np            = T%np + 1
             mypix(T%np,1)   = j
             mypix(T%np,2)   = data(band)%info%pix(j)
             mybeam(T%np,:)  = b(i,:)
             i               = i+1
             j               = j+1
          else if (ind(i) < data(band)%info%pix(j)) then
             i               = i+1
          else
             j               = j+1
          end if
          if (i > n) exit
          if (j > data(band)%info%np) exit
       end do
    end if

    ! Store pixels that belong to current processor
    allocate(T%pix(T%np,2), T%map(T%np,T%nmaps), T%Omega_b(T%nmaps))
    T%pix = mypix(1:T%np,:)
    do i = 1, T%nmaps
       T%map(:,i)   = mybeam(1:T%np,i) / maxval(b(:,i))
       T%Omega_b(i) = sum(b(:,i))/maxval(b(:,i)) * 4.d0*pi/(12.d0*T%nside**2)
    end do

    deallocate(ind, b, mypix, mybeam)
    
  end subroutine read_febecop_beam


  subroutine compute_symmetric_beam(glon, glat, T, bl, br)
    implicit none
    real(dp),  intent(in)    :: glon, glat
    type(Tnu), intent(inout) :: T
    real(dp),  dimension(0:,1:), intent(in), optional :: bl
    real(dp),  dimension(1:),    intent(in), optional :: br
          
!!$    ! Search for pixels controlled by current processor
!!$    npix = 12*constructor%T(i)%nside**2
!!$    call ang2vec(0.5d0*pi-constructor%glat, constructor%glon, vec0)
!!$    call query_disc(constructor%T(i)%nside, vec0, data(i)%B%r_max, listpix, nlist)
!!$    constructor%T(i)%np = 0
!!$    k                   = 1   ! map counter
!!$    j                   = 0   ! source counter
!!$    do while (k <= data(i)%info%np .and. j <= nlist-1)
!!$       if (listpix(j) == data(i)%info%pix(k)) then
!!$          constructor%T(i)%np       = constructor%T(i)%np + 1
!!$          hits(constructor%T(i)%np) = k
!!$          j                         = j+1
!!$          k                         = k+1
!!$       else if (listpix(j) < data(i)%info%pix(k)) then
!!$          j = j+1
!!$       else
!!$          k = k+1
!!$       end if
!!$    end do
!!$    
!!$    allocate(constructor%T(i)%pix(constructor%T(i)%np,2))
!!$    allocate(constructor%T(i)%map(constructor%T(i)%np,constructor%T(i)%nmaps))
!!$    allocate(constructor%T(i)%F(constructor%T(i)%nmaps))
!!$    constructor%T(i)%pix(:,2) = hits(1:constructor%T(i)%np)
!!$    do j = 1, constructor%T(i)%np
!!$       constructor%T(i)%pix(j,1) = data(i)%info%pix(constructor%T(i)%pix(j,2))
!!$       call pix2vec_ring(constructor%T(i)%nside, constructor%T(i)%pix(j,1), vec)
!!$       call angdist(vec0, vec, r)
!!$       constructor%T(i)%map(j,:) = r
!!$    end do

  end subroutine compute_symmetric_beam

  subroutine initPtsrcPrecond(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    integer(i4b) :: i, i1, i2, j, j1, j2, k1, k2, q, l, m, n, p, p1, p2, n1, n2, myid, ierr
    real(dp)     :: t1, t2
    class(comm_comp),         pointer :: c, c1, c2
    class(comm_ptsrc_comp),   pointer :: pt1, pt2
    real(dp),     allocatable, dimension(:,:) :: mat, mat2

    if (ncomp_pre == 0) return

    call mpi_comm_rank(comm, myid, ierr)
    
!!$    if (.not. allocated(ptsrcComps)) then
!!$       ! Set up an array of all the diffuse components
!!$       write(*,*) ncomp_pre
!!$       allocate(ptsrcComps(ncomp_pre))
!!$       c => compList
!!$       i =  1
!!$       do while (associated(c))
!!$          select type (c)
!!$          class is (comm_ptsrc_comp)
!!$             ptsrcComps(i)%p => c
!!$             write(*,*) i
!!$             write(*,*) 'a', c%src(1)%T(1)%np
!!$             write(*,*) 'b', c%src(1)%T(2)%np
!!$             write(*,*) 'c', c%src(2)%T(1)%np
!!$             write(*,*) 'd', c%src(2)%T(2)%np
!!$             write(*,*) 'a', ptsrcComps(i)%p%src(1)%T(1)%np
!!$             write(*,*) 'b', ptsrcComps(i)%p%src(1)%T(2)%np
!!$             write(*,*) 'c', ptsrcComps(i)%p%src(2)%T(1)%np
!!$             write(*,*) 'd', ptsrcComps(i)%p%src(2)%T(2)%np
!!$             i               =  i+1
!!$          end select
!!$          c => c%next()
!!$       end do
!!$    end if
    
    ! Build frequency-dependent part of preconditioner
    call wall_time(t1)
    allocate(P_cr%invM_src(1,nmaps_pre))
    allocate(mat(npre,npre), mat2(npre,npre))
    do j = 1, nmaps_pre

       mat = 0.d0
       i1  = 0
       c1 => compList
       do while (associated(c1))
          select type (c1)
          class is (comm_ptsrc_comp)
             pt1 => c1
          class default 
             c1 => c1%next()
             cycle
          end select
          if (j > pt1%nmaps) then
             c1 => c1%next()
             cycle
          end if
          do k1 = 1, pt1%nsrc
             !write(*,*) k1, pt1%nsrc             
             i1 = i1+1

             i2 = 0
             c2 => compList
             do while (associated(c2))
                !do j2 = 1, ncomp_pre
                select type (c2)
                class is (comm_ptsrc_comp)
                   pt2 => c2
                class default 
                   c2 => c2%next()
                   cycle
                end select
                
                if (j > pt2%nmaps) then
                   c2 => c2%next()
                   cycle
                end if
                do k2 = 1, pt2%nsrc
                   !write(*,*) k2, pt2%nsrc
                   i2 = i2+1
                   if (i2 < i1) cycle

                   do l = 1, numband
                      n1 = pt1%src(k1)%T(l)%np
                      n2 = pt2%src(k2)%T(l)%np

                      ! Search for common pixels; skip if no pixel overlap
                      if (pt1%src(k1)%T(l)%pix(1,1)  > pt2%src(k2)%T(l)%pix(n2,1)) cycle
                      if (pt1%src(k1)%T(l)%pix(n1,1) < pt2%src(k2)%T(l)%pix(1,1))  cycle

                      p1 = 1
                      p2 = 1
                      do while (.true.)
                         if (pt1%src(k1)%T(l)%pix(p1,1) == pt2%src(k2)%T(l)%pix(p2,1)) then
                            p  = pt1%src(k1)%T(l)%pix(p1,1)
                            mat(i1,i2) = mat(i1,i2) + &
                                 & data(l)%N%invN_diag%map(p,j) * &          ! invN_{p,p}
                                 & pt1%src(k1)%T(l)%map(p1,j) * & ! B_1
                                 & pt2%src(k2)%T(l)%map(p2,j) * & ! B_2
                                 & pt1%src(k1)%T(l)%F(j)      * & ! F_1
                                 & pt2%src(k2)%T(l)%F(j)      * & ! F_2
                                 & pt1%getScale(l,k1,j)       * & ! Unit 1
                                 & pt2%getScale(l,k2,j)           ! Unit 2
                            p1 = p1+1
                            p2 = p2+1                            
                         else if (pt1%src(k1)%T(l)%pix(p1,1) < pt2%src(k2)%T(l)%pix(p2,1)) then
                            p1 = p1+1
                         else
                            p2 = p2+1
                         end if
                         if (p1 >= n1 .or. p2 >= n2) exit
                         if (pt1%src(k1)%T(l)%pix(p1,1) > pt2%src(k2)%T(l)%pix(n2,1)) exit
                         if (pt1%src(k1)%T(l)%pix(n1,1) < pt2%src(k2)%T(l)%pix(p2,1)) exit
                      end do
                   end do
                   mat(i2,i1) = mat(i1,i2)
                end do
                c2 => c2%next()
             end do
          end do
          c1 => c1%next()
       end do

       ! Collect contributions from all cores
       call mpi_reduce(mat, mat2, size(mat2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
       if (myid == 0) then
          call invert_matrix_with_mask(mat2)
          allocate(P_cr%invM_src(1,j)%M(npre,npre))
          P_cr%invM_src(1,j)%M = mat2
       end if
    end do
    call wall_time(t2)
    write(*,*) 'ptsrc precond init = ', real(t2-t1,sp)

    deallocate(mat,mat2)
    
  end subroutine initPtsrcPrecond

  subroutine updatePtsrcPrecond
    implicit none

    ! Placeholder for now; already fully initialized
    if (npre == 0) return
       
  end subroutine updatePtsrcPrecond


  subroutine applyPtsrcPrecond(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, j, k, l, m, nmaps
    real(dp), allocatable, dimension(:,:) :: amp
    real(dp), allocatable, dimension(:,:) :: y
    class(comm_comp),       pointer :: c
    class(comm_ptsrc_comp), pointer :: pt

    if (npre == 0) return
    
    ! Reformat linear array into y(npre,nalm,nmaps) structure
    allocate(y(npre,nmaps_pre))
    y = 0.d0
    l = 1
    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_ptsrc_comp)
          pt => c
       class default 
          c => c%next()
          cycle
       end select
       call cr_extract_comp(pt%id, x, amp)
       do k = 1, pt%nmaps
          y(l:l+pt%nsrc-1,k) = amp(:,k)
       end do
       l  = l + pt%nsrc
       c => c%next()
       deallocate(amp)
    end do

    ! Multiply with preconditioner
    do j = 1, nmaps_pre
       y(:,j) = matmul(P_cr%invM_src(1,j)%M, y(:,j))
    end do

    ! Reformat y(npre,nmaps) structure back into linear array
    l = 1
    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_ptsrc_comp)
          pt => c
       class default 
          c => c%next()
          cycle
       end select
       allocate(amp(pt%nsrc,pt%nmaps))
       do k = 1, pt%nmaps
          amp(:,k) = y(l:l+pt%nsrc-1,k)
       end do
       call cr_insert_comp(pt%id, .false., amp, x)
       l = l + pt%nsrc
       c => c%next()
       deallocate(amp)
    end do
        
    deallocate(y)

  end subroutine applyPtsrcPrecond

  function getScale(self, band, id, pol)
    implicit none
    class(comm_ptsrc_comp), intent(in) :: self
    integer(i4b),           intent(in) :: band, id, pol
    real(dp)                           :: getScale

    if (trim(self%type) == 'radio') then
       getScale = 1.d-23 * (c/self%nu_ref)**2 / (2.d0*k_b*self%src(id)%T(band)%Omega_b(pol))
    end if

  end function getScale
  
end module comm_ptsrc_comp_mod
