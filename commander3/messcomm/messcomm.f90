program messcomm
  use comm_utils
  use comm_param_mod
  use alm_tools
  use math_tools
  implicit none

  type data_type
     integer(i4b) :: nside, npix, nmaps, lmax
     logical(lgt) :: pol
     character(len=16) :: label
     real(dp),     allocatable, dimension(:)     :: rms_uni, N_ell_uni
     real(dp),     allocatable, dimension(:,:)   :: map, sig, rms_tot, rms_res, T
     complex(dpc), allocatable, dimension(:,:,:) :: T_lm
     real(dp),     allocatable, dimension(:,:)   :: bl
     real(dp),     allocatable, dimension(:,:)   :: F
  end type data_type

  type comp_type
     character(len=16) :: label
     integer(i4b)      :: nside, npix, nmaps, lmax
     real(dp)          :: fwhm
     logical(lgt)      :: pol
     real(dp),     allocatable, dimension(:,:)   :: map, bl
     complex(dpc), allocatable, dimension(:,:,:) :: alm
  end type comp_type

  type(data_type), allocatable, dimension(:) :: data
  type(comp_type), allocatable, dimension(:) :: comp

  integer(i4b)       :: i, j, k, l, m
  integer(i4b)       :: nband, ncomp, niter, seed, lmax, nmaps
  real(dp)           :: noisefrac
  character(len=2)   :: itext2
  character(len=3)   :: itext3
  character(len=512) :: parfile, filename, pixwin
  type(planck_rng)   :: handle
  real(dp), allocatable, dimension(:,:) :: tmp


  call getarg(1, parfile)

  call get_parameter(parfile, 'NUMBAND', par_int=nband)
  call get_parameter(parfile, 'NUMCOMP', par_int=ncomp)
  call get_parameter(parfile, 'NUMITER', par_int=niter)
  call get_parameter(parfile, 'SEED',    par_int=seed)

  call get_parameter(parfile, 'MESSENGER_NOISE_FRACTION', par_dp=noisefrac)

  call rand_init(handle, seed)

  ! Initialize data structure
  allocate(data(nband))
  do i = 1, ncomp
     call int2string(i, itext3)

     call get_parameter(parfile, 'LABEL'//itext3, par_string=data(i)%label)
     call get_parameter(parfile, 'NSIDE'//itext3, par_int=data(i)%nside)
     call get_parameter(parfile, 'LMAX'//itext3,  par_int=data(i)%lmax)
     call get_parameter(parfile, 'POLAR'//itext3, par_lgt=data(i)%pol)
     data(i)%nmaps = 1; if (data(i)%pol) data(i)%nmaps = 3
     data(i)%npix = 12*data(i)%nside**2

     allocate(data(i)%map(0:data(i)%npix-1,data(i)%nmaps))
     allocate(data(i)%sig(0:data(i)%npix-1,data(i)%nmaps))
     allocate(data(i)%rms_tot(0:data(i)%npix-1,data(i)%nmaps))
     allocate(data(i)%rms_uni(data(i)%nmaps))
     allocate(data(i)%rms_res(0:data(i)%npix-1,data(i)%nmaps))
     allocate(data(i)%bl(0:data(i)%lmax,data(i)%nmaps))
     allocate(data(i)%F(data(i)%nmaps,ncomp))
     allocate(data(i)%T(0:data(i)%npix-1,data(i)%nmaps))
     allocate(data(i)%T_lm(data(i)%nmaps,0:data(i)%lmax,0:data(i)%lmax))
     data(i)%sig = 0.d0

     ! Read map file
     call get_parameter(parfile, 'MAPFILE'//itext3, par_string=filename)
     call read_map(filename, data(i)%map)

     ! Read noise
     call get_parameter(parfile, 'RMSFILE'//itext3, par_string=filename)
     call read_map(filename, data(i)%rms_tot)
     do j = 1, data(i)%nmaps
        data(i)%rms_uni(j)   = minval(data(i)%rms_tot(:,j)) * noisefrac
        data(i)%rms_res(:,j) = sqrt(data(i)%rms_tot(:,j)**2 - data(i)%rms_uni(j)**2)
     end do
     data(i)%N_ell_uni = data(i)%rms_uni**2 * 4.d0*pi/data(i)%npix

     ! Read beam
     call get_parameter(parfile, 'BEAMFILE'//itext3, par_string=filename)
     call get_parameter(parfile, 'PIXWIN'//itext3,   par_string=pixwin)
     call read_beam(data(i)%lmax, data(i)%nmaps, data(i)%bl, filename, pixwin=pixwin)

     ! Read mixing matrices
     allocate(tmp(0:data(i)%npix-1,data(i)%nmaps))
     do j = 1, ncomp
        call int2string(j, itext2)
        call get_parameter(parfile, 'MIXMAT'//itext3//'_'//itext2, par_string=filename)
        call read_map(filename, tmp)
        do k = 1, data(i)%nmaps
           data(i)%F(k,j) = sum(tmp(:,k))/size(tmp(:,k))
        end do
     end do
     deallocate(tmp)
  end do

  ! Initialize component structure
  allocate(comp(ncomp))
  lmax  = 0
  nmaps = 0
  do i = 1, ncomp
     call int2string(i, itext2)
     call get_parameter(parfile, 'COMP_LABEL'//itext2, par_string=comp(i)%label)
     call get_parameter(parfile, 'COMP_NSIDE'//itext2, par_int=comp(i)%nside)
     call get_parameter(parfile, 'COMP_LMAX'//itext2,  par_int=comp(i)%lmax)
     call get_parameter(parfile, 'COMP_FWHM'//itext2,  par_dp=comp(i)%fwhm)
     call get_parameter(parfile, 'COMP_POLAR'//itext2, par_lgt=comp(i)%pol)
     comp(i)%nmaps = 1; if (comp(i)%pol) comp(i)%nmaps = 3
     comp(i)%npix = 12*comp(i)%nside**2
     lmax         = max(lmax,  comp(i)%lmax)
     nmaps        = max(nmaps, comp(i)%nmaps)

     allocate(comp(i)%map(0:comp(i)%npix-1,comp(i)%nmaps))
     allocate(comp(i)%bl(0:comp(i)%lmax,comp(i)%nmaps))
     allocate(comp(i)%alm(comp(i)%nmaps,0:comp(i)%lmax,0:comp(i)%lmax))
     comp(i)%map = 0.d0
     comp(i)%alm = cmplx(0.d0, 0.d0)
     call read_beam(comp(i)%lmax, comp(i)%nmaps, comp(i)%bl, fwhm=comp(i)%fwhm)

  end do


  ! Run Gibbs chain
  do i = 1, niter
     write(*,*) 'Iter ', i, ' of ', niter
     call sample_T_given_s
     call sample_s_given_T
     call output_sample(i)
  end do

contains

  subroutine sample_T_given_s
    implicit none

    integer(i4b) :: i, j, k
    real(dp)     :: mu, var
    real(dp), allocatable, dimension(:,:) :: weights

    do i = 1, nband
       do k = 1, data(i)%nmaps
          do j = 0, data(i)%npix-1
             var = 1.d0 / (1.d0/data(i)%rms_uni(k)**2 + 1.d0/data(i)%rms_res(j,k)**2)
             mu  = var * (data(i)%map(j,k)/data(i)%rms_res(j,k)**2 + &
                        & data(i)%sig(j,k)/data(i)%rms_uni(k)**2)
             data(i)%T(j,k) = mu + sqrt(var) * rand_gauss(handle)
          end do
       end do

       allocate(weights(2*data(i)%nside,data(i)%nmaps))
       weights = 1.d0
       if (data(i)%nmaps == 1) then
          call map2alm(data(i)%nside, data(i)%lmax, data(i)%lmax, &
               & data(i)%T(:,1), data(i)%T_lm, [0.d0, 0.d0], weights)
       else
          call map2alm(data(i)%nside, data(i)%lmax, data(i)%lmax, &
               & data(i)%T, data(i)%T_lm, [0.d0, 0.d0], weights)
       end if
    end do

  end subroutine sample_T_given_s


  subroutine sample_s_given_T
    implicit none

    integer(i4b) :: i, j, k, l, m, p
    real(dp)     :: val
    real(dp),     allocatable, dimension(:  )   :: mu, eta
    real(dp),     allocatable, dimension(:,:)   :: cov, L_cov, map
    complex(dpc), allocatable, dimension(:,:,:) :: alm

    do i = 1, ncomp
       comp(i)%alm = 0.d0
    end do

    
    allocate(mu(ncomp), cov(ncomp,ncomp), eta(ncomp), L_cov(ncomp, ncomp))
    do p = 1, nmaps
       do l = 0, lmax
          cov = 0.d0
          do k = 1, nband
             do i = 1, ncomp
                do j = 1, ncomp
                   cov(i,j) = cov(i,j) + data(k)%F(p,i) * data(k)%bl(l,p) * &
                        & 1.d0/data(k)%N_ell_uni(p)**2 * data(k)%bl(l,p) * data(k)%F(p,j)
                end do
             end do
          end do
          call invert_matrix(cov)             
          call cholesky_decompose(cov, L_cov)

          do m = -l, l
             mu  = 0.d0
             do k = 1, nband
                do i = 1, ncomp
                   if (m < 0) then
                      val = imag(data(k)%T_lm(p,l,-m))
                   else
                      val = real(data(k)%T_lm(p,l,m))
                   end if
                   mu(i) = mu(i) + data(k)%F(p,i) * data(k)%bl(l,p) * &
                        & 1.d0/data(k)%N_ell_uni(p)**2 * val
                end do
             end do
             do k = 1, ncomp
                eta(k) = rand_gauss(handle)
             end do
             mu = matmul(cov, mu) + matmul(L_cov, eta)
             do k = 1, ncomp
                if (m < 0) then
                   comp(k)%alm(p,l,-m) = comp(k)%alm(p,l,-m) + cmplx(0.d0,mu(k))
                else
                   comp(k)%alm(p,l,m) = comp(k)%alm(p,l,m) + real(mu(k))
                end if
             end do
          end do
       end do
    end do
    deallocate(mu, cov, eta, L_cov)

    ! Sum up signal
    do i = 1, nband
       data(i)%sig = 0.d0
       allocate(map(0:data(i)%npix-1,data(i)%nmaps), &
            & alm(data(i)%nmaps,0:data(i)%lmax,0:data(i)%lmax))
       do j = 1, ncomp
          alm = comp(j)%alm(1:data(i)%nmaps,0:data(i)%lmax,0:data(i)%lmax)
          do l = 0, data(i)%lmax
             do m = 0, l
                alm(:,l,m) = alm(:,l,m) * data(i)%bl(l,:)
             end do
          end do
          if (data(i)%nmaps == 1) then
             call alm2map(data(i)%nside, data(i)%lmax, data(i)%lmax, alm, map(:,1))
          else
             call alm2map(data(i)%nside, data(i)%lmax, data(i)%lmax, alm, map)
          end if
          do p = 1, data(i)%nmaps
             data(i)%sig(:,p) = data(i)%sig(:,p) + data(i)%F(p,j) * map(:,p)
          end do
       end do
       deallocate(map, alm)
    end do

    ! Compute component maps
    do i = 1, ncomp
       allocate(alm(comp(i)%nmaps,0:comp(i)%lmax,0:comp(i)%lmax))
       alm = comp(i)%alm
       do l = 0, data(i)%lmax
          do m = 0, l
             alm(:,l,m) = alm(:,l,m) * comp(i)%bl(l,:)
          end do
       end do
       if (comp(i)%nmaps == 1) then
          call alm2map(comp(i)%nside, comp(i)%lmax, comp(i)%lmax, alm, comp(i)%map(:,1))
       else
          call alm2map(comp(i)%nside, comp(i)%lmax, comp(i)%lmax, alm, comp(i)%map)
       end if
    end do
    
  end subroutine sample_s_given_T


  subroutine output_sample(iter)
    implicit none
    integer(i4b), intent(in) :: iter

    integer(i4b) :: i
    character(len=512) :: filename
    character(len=4)   :: iter_text

    call int2string(iter, iter_text)
    do i = 1, ncomp
       filename = trim(comp(i)%label)//'_c'//iter_text//'.fits'
       call write_map(filename, comp(i)%map)
    end do

    do i = 1, nband
       filename = 'res_'//trim(data(i)%label)//'_c'//iter_text//'.fits'
       call write_map(filename, data(i)%map-data(i)%sig)

       filename = 't_'//trim(data(i)%label)//'_c'//iter_text//'.fits'
       call write_map(filename, data(i)%map-data(i)%T)
    end do

  end subroutine output_sample

  subroutine write_map(filename, map, comptype, nu_ref, unit, ttype, spectrumfile)
    implicit none

    character(len=*),                   intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(in)  :: map
    character(len=*),                   intent(in), optional :: comptype, unit, spectrumfile, ttype
    real(dp),                           intent(in), optional :: nu_ref

    integer(i4b)   :: npix, nlheader, nmaps, i, nside
    logical(lgt)   :: exist, polarization

    character(len=80), dimension(1:120)    :: header
    character(len=16) :: unit_, ttype_

    npix         = size(map(:,1))
    nside        = nint(sqrt(real(npix,sp)/12.))
    nmaps        = size(map(0,:))
    polarization = (nmaps == 3)
    unit_        = '';       if (present(unit)) unit_  = unit
    ttype_       = 'Stokes'; if (present(unit)) ttype_ = ttype


    !-----------------------------------------------------------------------
    !                      write the map to FITS file
    !  This is copied from the synfast.f90 file in the Healpix package
    !-----------------------------------------------------------------------
    
    nlheader = SIZE(header)
    do i=1,nlheader
       header(i) = ""
    enddo

    ! start putting information relative to this code and run
    call add_card(header)
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Sky Map Pixelisation Specific Keywords    ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")
    call add_card(header,"ORDERING","RING",  "Pixel ordering scheme, either RING or NESTED")
    call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
    call add_card(header,"FIRSTPIX",0,"First pixel # (0 based)")
    call add_card(header,"LASTPIX",npix-1,"Last pixel # (0 based)")
    call add_card(header,"BAD_DATA",  HPX_DBADVAL ,"Sentinel value given to bad pixels")
    call add_card(header) ! blank line
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Data Description Specific Keywords       ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"POLCCONV","COSMO"," Coord. convention for polarisation (COSMO/IAU)")
    call add_card(header,"INDXSCHM","IMPLICIT"," Indexing : IMPLICIT or EXPLICIT")
    call add_card(header,"GRAIN", 0, " Grain of pixel indexing")
    call add_card(header,"COMMENT","GRAIN=0 : no indexing of pixel data                         (IMPLICIT)")
    call add_card(header,"COMMENT","GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)")
    call add_card(header,"COMMENT","GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)")
    call add_card(header) ! blank line
    call add_card(header,"POLAR",polarization," Polarisation included (True/False)")

    call add_card(header) ! blank line
    call add_card(header,"TTYPE1", "I_"//ttype_,"Stokes I")
    call add_card(header,"TUNIT1", unit_,"Map unit")
    call add_card(header)

    if (polarization) then
       call add_card(header,"TTYPE2", "Q_"//ttype_,"Stokes Q")
       call add_card(header,"TUNIT2", unit_,"Map unit")
       call add_card(header)
       
       call add_card(header,"TTYPE3", "U_"//ttype_,"Stokes U")
       call add_card(header,"TUNIT3", unit_,"Map unit")
       call add_card(header)
    endif
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Commander Keywords                        ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","Commander is a code for global CMB analysis    ")
    call add_card(header,"COMMENT","developed in collaboration between the University")
    call add_card(header,"COMMENT","of Oslo and Jet Propulsion Laboratory (NASA).  ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    if (present(comptype)) call add_card(header,"COMPTYPE",trim(comptype), "Component type")
    if (present(nu_ref))   call add_card(header,"NU_REF",  nu_ref,         "Reference frequency")
    if (present(spectrumfile)) call add_card(header,"SPECFILE",  trim(spectrumfile), &
         & "Reference spectrum")
    call add_card(header,"COMMENT","-----------------------------------------------")


    !call write_bintab(map, npix, nmaps, header, nlheader, filename)
    call output_map(map, header, "!"//trim(filename))

  end subroutine write_map
  
end program messcomm
