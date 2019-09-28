module comm_4D_map_mod
  use comm_utils
  use hashtbl_4Dmap
  implicit none

  private
  public comm_4D_map, output_4D_maps

  type ::  comm_4D_map
     integer(i4b) :: n, nside, npsi, ndet
     character(len=16), allocatable, dimension(:)   :: detlabel
     real(dp),          allocatable, dimension(:)   :: psi0
     integer(i4b),      allocatable, dimension(:)   :: pixel, ipsi
     real(sp),          allocatable, dimension(:)   :: weight
     real(sp),          allocatable, dimension(:,:) :: map
  end type comm_4D_map

  interface comm_4D_map
     procedure constructor
  end interface comm_4D_map

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(nside, npsi, detlabel, psi0, pixel, psi, tod, mask)
    implicit none
    integer(i4b),                          intent(in) :: nside, npsi
    character(len=*),      dimension(:),   intent(in) :: detlabel
    real(sp),              dimension(:),   intent(in) :: psi0
    integer(i4b),          dimension(:),   intent(in) :: pixel, psi
    real(sp),              dimension(:,:), intent(in) :: tod
    integer(i4b),          dimension(:,:), intent(in) :: mask
    class(comm_4D_map),    pointer                    :: constructor

    integer(i4b) :: i, j, n, ndet
    type(hash_tbl_4Dmap_sll) :: hashtbl

    ! Initialize object
    allocate(constructor)
    ndet                        = size(detlabel)
    constructor%ndet     = ndet
    allocate(constructor%detlabel(ndet), constructor%psi0(ndet))
    constructor%nside    = nside
    constructor%npsi     = npsi
    constructor%detlabel = detlabel
    constructor%psi0     = psi0

    ! Construct hashed 4D map
    call init_hash_tbl_4Dmap_sll(hashtbl)
    do i = 1, size(pixel)
       if (any(mask(i,:) /= 0)) cycle
       call hashtbl%put([pixel(i),psi(i)], tod(i,:))
    end do
    n = hashtbl%get_n_elements()
    constructor%n = n

    ! Populate 4D map with hashed values
    allocate(constructor%pixel(n), constructor%ipsi(n))
    allocate(constructor%weight(n), constructor%map(n,ndet))
    call hashtbl%linearize(constructor%pixel, constructor%ipsi, &
         & constructor%weight, constructor%map)

    ! Clean up
    call free_hash_tbl_4Dmap_sll(hashtbl)

  end function constructor


  !**************************************************
  !             Output routines
  !**************************************************
  subroutine output_4D_maps(prefix, postfix, nside, npsi, detlabel, horn_id, psi0, sigma0, &
       & pixel, psi, tod, mask)
    implicit none
    character(len=*),                      intent(in) :: prefix, postfix
    integer(i4b),                          intent(in) :: nside, npsi
    character(len=*),      dimension(:),   intent(in) :: detlabel
    integer(i4b),          dimension(:),   intent(in) :: horn_id
    real(sp),              dimension(:),   intent(in) :: psi0, sigma0
    integer(i4b),          dimension(:,:), intent(in) :: pixel, psi
    real(sp),              dimension(:,:), intent(in) :: tod
    integer(i4b),          dimension(:,:), intent(in) :: mask

    integer(i4b) :: i, j, h, horn, nhorn, ndet
    integer(i4b), allocatable, dimension(:) :: d
    character(len=1)   :: itext
    character(len=16)  :: dlabel
    character(len=512) :: filename
    character(len=80), dimension(180) :: header
    class(comm_4D_map), pointer       :: map4D
    
    integer   :: status, unit, readwrite, blocksize, hdutype, tfields, nrows, bitpix, naxis, naxes(0)
    integer   :: varidat, colnum,frow,felem
    logical   :: simple, extend
    character :: extname*16
    character(len=16), allocatable, dimension(:) :: ttype, tform, tunit

    ! Find number of unique horns
    ndet  = size(detlabel)
    nhorn = 0
    do horn = minval(horn_id), maxval(horn_id)
       if (count(horn_id == horn) == 0) cycle
       nhorn = nhorn + 1       
    end do

    ! Construct 4D maps for each horn, and output each in a separate HDU
    h = 0
    do horn = minval(horn_id), maxval(horn_id)
       ndet = count(horn_id == horn)
       if (ndet == 0) cycle
       h = h+1

       ! Find which detectors belong to current horn
       allocate(d(ndet))
       ndet = 0
       do i = 1, size(horn_id)
          if (horn_id(i) == horn) then
             ndet       = ndet+1
             d(ndet) = i
          end if
       end do
       
       if (psi0(d(2)) /= 0.d0) then
          write(*,*) 'Error: psi0 for second detector is non-zero'
          stop
       end if
       
       ! Construct 4D map for current horn
       map4D => comm_4D_map(nside, npsi, detlabel(d), psi0(d), pixel(:,d(2)), psi(:,d(2)), &
            & tod(:,d), mask(:,d))
       
       ! Output file
       dlabel   = detlabel(d(1))
       filename = trim(prefix) // "_" // dlabel(1:2) // trim(postfix)

       allocate(ttype(ndet+3), tform(ndet+3), tunit(ndet+3))
       ttype(1:3) = ['pixel', 'ipsi', 'weight']
       tform(1:3) = ['1J',    '1I',   '1E'    ]
       tunit(1:3) = ['none',  'none', 'none'  ]
       do j = 1, ndet
          call int2string(j, itext)
          ttype(3+j) = 'signal'//itext
          tform(3+j) = '1E'
          tunit(3+j) = 'none'
       end do

       status    = 0
       readwrite = 1
       blocksize = 1
       bitpix    = 8
       naxis     = 0
       unit      = getlun()
       call ftinit(unit,filename,blocksize,status)
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       ! Output general information in HDU1
       call ftprec(unit,'----------------------------------------------------', status)
       call ftprec(unit,'                   CODE                             ', status)
       call ftprec(unit,'----------------------------------------------------', status)
       call ftpkys(unit,"CREATOR", "Commander", "Software creating FITS file",  status)
       call ftprec(unit,'----------------------------------------------------', status)
       call ftprec(unit,'             GENERAL PARAMETERS                     ', status)
       call ftprec(unit,'----------------------------------------------------', status)
       call ftpkys(unit,"PIXTYPE",  "HEALPIX",  "",                             status)
       call ftpkys(unit,"COORDSYS", "GALACTIC", "Coordinate system",            status)
       call ftpkys(unit,"POLCCONV", "COSMO",    "Polarization convention",      status)
       call ftpkys(unit,"ORDERING", "RING",     "Healpix ordering",             status)
       call ftpkyj(unit,"NSIDE",    nside,      "Healpix Nside",                status)
       call ftprec(unit,'----------------------------------------------------', status)
       call ftprec(unit,'             DETECTOR INFORMATION                   ', status)
       call ftprec(unit,'----------------------------------------------------', status)
       call ftpkyj(unit,"NODET",    ndet,       "Number of detectors",          status)
       
       do i = 1, ndet
          call int2string(i, itext)
          call ftpkys(unit,"DETNAM"//itext, trim(detlabel(i)), "Detector name",           status)
          call ftpkye(unit,"PSIPOL"//itext, psi0(i),    6,     "Psi_pol/deg",             status)
          !call ftpkyj(unit,"IPOINT"//itext, horn_id(i),        "Pointing ID",             status)
          call ftpkye(unit,"SIG"//itext,    sigma0(i),  6,     "White noise RMS (K_cmb)", status)
       end do


       call ftmahd(unit,1,hdutype,status)
       call ftcrhd(unit,status)
       
       tfields = 3+ndet
       nrows   = map4D%n
       extname = 'xtension'
       varidat = 0
       call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

       frow   = 1
       felem  = 1
       colnum = 1
       call ftpclj(unit,colnum,frow,felem,nrows,map4D%pixel,status)  
       colnum = 2
       call ftpcli(unit,colnum,frow,felem,nrows,map4D%ipsi,status)  
       colnum = 3
       call ftpcle(unit,colnum,frow,felem,nrows,map4D%weight,status)  
       do j = 1, ndet
          colnum = 3+j
          call ftpcle(unit,colnum,frow,felem,nrows,map4D%map(:,j),status)  
       end do
       call ftpkys(unit,"Ordering", "RING",       "Pixel ordering scheme",                     status)
       call ftpkyj(unit,"Nside",    nside,        "Resolution parameter for HEALPIX",          status)
       call ftpkyj(unit,"Npsi",     npsi,         "Resolution parameter for beam orientation", status)
       call ftpkye(unit,"psi0",     0.,        6, "Lower boundary for first psi bin",          status)
       call ftpkye(unit,"sigma1",   sigma0(1), 6, "Sigma for detector 1",                      status)
       call ftpkye(unit,"sigma2",   sigma0(2), 6, "Sigma for detector 2",                      status)
       call ftpkye(unit,"psipol1",  psi0(1),   6, "Sigma for detector 1",                      status)
       call ftpkye(unit,"psipol2",  psi0(2),   6, "Sigma for detector 2",                      status)
       call ftpkys(unit,"Coordsys", 'G',          "Coordinate system",                         status)
       call ftpkyl(unit,"Calibrated", .true.,     "",                                          status)
       call ftclos(unit, status)
       call ftfiou(unit, status)

       deallocate(d, ttype, tform, tunit)
    end do

  end subroutine output_4D_maps

end module comm_4D_map_mod
