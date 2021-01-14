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
module comm_4D_map_mod
  use comm_utils
  use hashtbl_4Dmap
  use comm_hdf_mod
  implicit none

  private
  public comm_4D_map, output_4D_maps, output_4D_maps_hdf

  type ::  comm_4D_map
     integer(i4b) :: n, nside, npsi, ndet
     character(len=16), allocatable, dimension(:)   :: detlabel
     real(dp),          allocatable, dimension(:)   :: psi0
     integer(i4b),      allocatable, dimension(:)   :: pixel, ipsi
     real(sp),          allocatable, dimension(:)   :: weight
     real(sp),          allocatable, dimension(:,:) :: map
   contains
     procedure :: dealloc => deallocate_4D_map
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

    integer(i4b) :: i, n, ndet
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
  subroutine output_4D_maps_hdf(filename, path, scanid, nside, npsi, detlabel, horn_id, psi0, sigma0, &
       & pixel, psi, tod, mask, accept)
    implicit none
    character(len=*),                      intent(in) :: filename, path
    integer(i4b),                          intent(in) :: scanid, nside, npsi
    character(len=*),      dimension(:),   intent(in) :: detlabel
    integer(i4b),          dimension(:),   intent(in) :: horn_id
    real(sp),              dimension(:),   intent(in) :: psi0, sigma0
    integer(i4b),          dimension(:,:), intent(in) :: pixel, psi
    real(sp),              dimension(:,:), intent(in) :: tod
    integer(i4b),          dimension(:,:), intent(in) :: mask
    logical(lgt),          dimension(:),   intent(in) :: accept

    integer(i4b) :: i, j, h, horn, nhorn, ndet, pid(1), nsamp(1), hdferr
    integer(i4b), allocatable, dimension(:) :: d
    logical(lgt) :: exist
    character(len=1)   :: itext
    character(len=6)   :: scantext
    character(len=16)  :: dlabel
    character(len=512) :: path_scanid, fullpath
    class(comm_4D_map), pointer       :: map4D => null()
    
    type(hdf_file)     :: file
    TYPE(h5o_info_t) :: object_info

    inquire(file=trim(filename), exist=exist)
    if (exist) then
       call open_hdf_file(filename, file, 'b')
    else
       call open_hdf_file(filename, file, 'w')
    end if

    call int2string(scanid, scantext)
    path_scanid = trim(adjustl(path))//'/'//scantext

    ! Create group if it doesnt already exists
    call h5eset_auto_f(0, hdferr)
    call h5oget_info_by_name_f(file%filehandle, trim(adjustl(path)), object_info, hdferr)
    if (hdferr /= 0) call create_hdf_group(file, trim(adjustl(path)))
    call h5oget_info_by_name_f(file%filehandle, trim(adjustl(path_scanid)), object_info, hdferr)
    if (hdferr /= 0) call create_hdf_group(file, trim(adjustl(path_scanid)))

!!$    if (hdferr == 0) call h5gunlink_f(file%filehandle, trim(adjustl(path)), hdferr)
!!$
!!$    call create_hdf_group(file, trim(adjustl(path)))
!!$    call create_hdf_group(file, trim(adjustl(path_scanid)))
    
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
       
       ! Check that current detectors are accepted
       if (any(.not. accept(d))) then
          deallocate(d)
          cycle
       end if
          
       ! Construct 4D map for current horn
       map4D => comm_4D_map(nside, npsi, detlabel(d), psi0(d), pixel(:,d(2)), psi(:,d(2)), &
            & tod(:,d), mask(:,d))
       
       ! Output file
       dlabel   = detlabel(d(1))
       fullpath = trim(adjustl(path_scanid))//'/'//dlabel(1:2)

       call h5eset_auto_f(0, hdferr)
       call h5oget_info_by_name_f(file%filehandle, trim(adjustl(fullpath)), object_info, hdferr)
       if (hdferr == 0) call h5gunlink_f(file%filehandle, trim(adjustl(fullpath)), hdferr)
       call create_hdf_group(file, trim(fullpath))

       ! Output general information in HDU1
       call write_hdf(file, trim(adjustl(fullpath))//'/nside', nside)
       call write_hdf(file, trim(adjustl(fullpath))//'/nodet', ndet)
       call write_hdf(file, trim(adjustl(fullpath))//"/npsi",  npsi)
       call write_hdf(file, trim(adjustl(fullpath))//"/psi0",  0.d0)
       call write_hdf(file, trim(adjustl(fullpath))//"/horn",  horn)
       call write_hdf(file, trim(adjustl(fullpath))//"/firstPointingID",  scanid)
       call write_hdf(file, trim(adjustl(fullpath))//"/lastPointingID",  scanid)

       do i = 1, ndet
          call int2string(i, itext)
          call write_hdf(file, trim(adjustl(fullpath))//"/detnam"//itext, trim(detlabel(d(i))))
          call write_hdf(file, trim(adjustl(fullpath))//"/psipol"//itext, psi0(d(i)))
          call write_hdf(file, trim(adjustl(fullpath))//"/sig"//itext,    sigma0(d(i)))
       end do

       call write_hdf(file, trim(adjustl(fullpath))//"/pixel",    map4D%pixel)
       call write_hdf(file, trim(adjustl(fullpath))//"/ipsi",    map4D%ipsi)
       call write_hdf(file, trim(adjustl(fullpath))//"/weight",    map4D%weight)
       call write_hdf(file, trim(adjustl(fullpath))//"/map",    map4D%map)

       call write_hdf(file, trim(adjustl(fullpath))//"/pid",    scanid)
       call write_hdf(file, trim(adjustl(fullpath))//"/sample_offset",    0)
       call write_hdf(file, trim(adjustl(fullpath))//"/nsamp",    map4d%n)
       call write_hdf(file, trim(adjustl(fullpath))//"/scet",    0.d0)
       call write_hdf(file, trim(adjustl(fullpath))//"/ecet",    0.d0)

       call map4D%dealloc; deallocate(map4D)
       deallocate(d)
    end do

    call close_hdf_file(file)

  end subroutine output_4D_maps_hdf


  subroutine output_4D_maps(prefix, postfix, scanid, nside, npsi, detlabel, horn_id, psi0, sigma0, &
       & pixel, psi, tod, mask, accept)
    implicit none
    character(len=*),                      intent(in) :: prefix, postfix
    integer(i4b),                          intent(in) :: scanid, nside, npsi
    character(len=*),      dimension(:),   intent(in) :: detlabel
    integer(i4b),          dimension(:),   intent(in) :: horn_id
    real(sp),              dimension(:),   intent(in) :: psi0, sigma0
    integer(i4b),          dimension(:,:), intent(in) :: pixel, psi
    real(sp),              dimension(:,:), intent(in) :: tod
    integer(i4b),          dimension(:,:), intent(in) :: mask
    logical(lgt),          dimension(:),   intent(in) :: accept

    integer(i4b) :: i, j, h, horn, nhorn, ndet, pid(1), nsamp(1)
    integer(i8b) :: sample_offset(1)
    real(dp)     :: scet(1), ecet(1)
    integer(i4b), allocatable, dimension(:) :: d
    character(len=1)   :: itext
    character(len=16)  :: dlabel
    character(len=512) :: filename
    !character(len=80), dimension(180) :: header
    class(comm_4D_map), pointer       :: map4D => null()
    
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
       
       ! Check that current detectors are accepted
       if (any(.not. accept(d))) then
          deallocate(d)
          cycle
       end if
          
!!$       if (psi0(d(2)) /= 0.d0) then
!!$          write(*,*) 'Error: psi0 for second detector is non-zero'
!!$          stop
!!$       end if
       
!!$       if (any(pixel(:,d(2)) < 0)) then
!!$          write(*,*) scanid
!!$          write(*,*) pixel(:,d(2))
!!$       end if

       ! Construct 4D map for current horn
       map4D => comm_4D_map(nside, npsi, detlabel(d), psi0(d), pixel(:,d(2)), psi(:,d(2)), &
            & tod(:,d), mask(:,d))
       
       ! Output file
       dlabel   = detlabel(d(1))
       filename = trim(prefix) // "_" // dlabel(1:2) // trim(postfix)

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

       allocate(ttype(ndet+3), tform(ndet+3), tunit(ndet+3))
       ttype(1:3) = ['pixel ', 'ipsi  ', 'weight']
       tform(1:3) = ['1J',    '1I',   '1E'    ]
       tunit(1:3) = ['none',  'none', 'none'  ]
       do j = 1, ndet
          call int2string(j, itext)
          ttype(3+j) = 'signal'//itext
          tform(3+j) = '1E'
          tunit(3+j) = 'none'
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
       call ftpcli(unit,colnum,frow,felem,nrows,int(map4D%ipsi,i2b),status)  
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
       call ftpkye(unit,"psipol1",  psi0(1),   6, "Polarization angle for detector 1",                      status)
       call ftpkye(unit,"psipol2",  psi0(2),   6, "Polarization angle for detector 2",                      status)
       call ftpkys(unit,"Coordsys", 'G',          "Coordinate system",                         status)
       call ftpkyl(unit,"Calibrated", .true.,     "",                                          status)
       call ftpkyj(unit,"horn",     horn,         "Internal horn ID number",                   status)
       call ftpkyj(unit,"firstPointingID", scanid, "",                                         status)
       call ftpkyj(unit,"lastPointingID",  scanid, "",                                         status)


       !call ftmahd(unit,2,hdutype,status)
       call ftcrhd(unit,status)

       ttype(1:5) = ['pointingID   ', 'sample_offset', 'nsamples_PID ', 'start_SCET   ', 'end_SCET     ']
       tform(1:5) = ['1J',         '1K',            '1J',           '1D',         '1D'      ]
       tunit(1:5) = ['none',       'none',          'none',         'mus ',        'mus '     ]
             
       tfields = 5
       nrows   = 1
       extname = 'xtension'
       varidat = 0
       call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

       frow   = 1
       felem  = 1
       colnum = 1
       pid = scanid
       call ftpclj(unit,colnum,frow,felem,nrows,pid,status)  
       colnum = 2
       sample_offset = 0
       call ftpclk(unit,colnum,frow,felem,nrows,sample_offset,status)  
       colnum = 3
       nsamp  = map4D%n
       call ftpclj(unit,colnum,frow,felem,nrows,nsamp,status)  
       colnum = 4
       scet(1)  = 0
       call ftpcld(unit,colnum,frow,felem,nrows,scet,status)  
       colnum = 5
       ecet(1)  = 0
       call ftpcld(unit,colnum,frow,felem,nrows,ecet,status)  

       call ftclos(unit, status)
       call ftfiou(unit, status)

       call map4D%dealloc; deallocate(map4D)
       deallocate(d, ttype, tform, tunit)
    end do

  end subroutine output_4D_maps

  subroutine deallocate_4D_map(self)
    implicit none
    class(comm_4D_map), intent(inout)  :: self

    if (allocated(self%detlabel)) deallocate(self%detlabel)
    if (allocated(self%psi0))     deallocate(self%psi0)
    if (allocated(self%pixel))    deallocate(self%pixel)
    if (allocated(self%ipsi))     deallocate(self%ipsi)
    if (allocated(self%weight))   deallocate(self%weight)
    if (allocated(self%map))      deallocate(self%map)

  end subroutine deallocate_4D_map



end module comm_4D_map_mod
