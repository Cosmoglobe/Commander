
! pixelization/scan definitions

subroutine f2c_set_pixelization( pixchoice, npars, parint, pixelization)

  USE healpix_types
  USE s2hat_types
  USE s2hat_pixelization

  implicit none
   
  !input 
  integer(i4b) :: pixchoice, npars
  integer(i4b), dimension(1:npars) :: parint
  ! output
  type( pixeltype), pointer :: pixelization
  ! internal
  type( pixparameters) :: pixpars

  pixpars%par1 = parint(1)
  pixpars%par2 = parint(2)

  pixelization => globalPixelization

  ! write(*,*)pixelization%npixsall, pixelization%nringsall

  call set_pixelization( pixchoice, pixpars, globalPixelization)

  ! write(*,*)pixelization%fpix(3), pixelization%nph(3), pixelization%cth(1), pixelization%sth(1)

  return  

end subroutine f2c_set_pixelization

subroutine f2c_destroy_pixelization( pixelization)

   USE s2hat_types
   USE s2hat_pixelization

   implicit none

   !input
   type(pixeltype), pointer :: pixelization

   pixelization => globalPixelization

   call destroy_pixelization( globalPixelization)

   return

end subroutine f2c_destroy_pixelization

subroutine f2c_zbounds2scan( zbounds, pixelization, scan)

  USE healpix_types
  USE s2hat_types
  USE s2hat_pixelization

  implicit none

  !input
  real(dp), dimension(2) :: zbounds
  type( pixeltype), pointer :: pixelization

  ! output
  type( scandef), pointer :: scan

  ! write(*,*)'pre',pixelization%npixsall, pixelization%nringsall

  pixelization => globalPixelization
  scan => globalScan

  ! write(*,*)'post',pixelization%npixsall, pixelization%nringsall

  call zbounds2scan( zbounds, globalPixelization, globalScan)
  
  ! write(*,*)scan%npixsobs, scan%nringsobs

  return

end subroutine f2c_zbounds2scan

subroutine f2c_zbounds2mask( zbounds, pixelization, mask)

  USE healpix_types
  USE s2hat_types
  USE s2hat_pixelization

  implicit none

  !input
  real(dp), dimension(2) :: zbounds
  type( pixeltype), pointer :: pixelization

  ! output
  integer(i4b), dimension(0:pixelization%npixsall-1), target :: mask

  ! internal
  integer(i4b), dimension(:), pointer :: mask_ptr  

  ! write(*,*)'pre',pixelization%npixsall, pixelization%nringsall

  pixelization => globalPixelization
  mask_ptr => mask

  ! write(*,*)'post',pixelization%npixsall, pixelization%nringsall

  call zbounds2mask( zbounds, globalPixelization, mask_ptr)
  
  ! write(*,*)scan%npixsobs, scan%nringsobs

  return

end subroutine f2c_zbounds2mask

subroutine f2c_mask2scan( mask, pixelization, scan)

  USE healpix_types
  USE s2hat_types
  USE s2hat_pixelization

  implicit none

  ! input
  type(pixeltype), pointer :: pixelization
  integer(i4b), dimension(0:pixelization%npixsall-1), target :: mask
  ! output
  type(scandef), pointer :: scan

  ! internal
  integer(i4b), dimension(:), pointer :: mask_ptr  

  pixelization => globalPixelization
  scan => globalScan

  mask_ptr => mask

  call mask2scan( mask_ptr, globalPixelization, globalScan)

  return
  
end subroutine f2c_mask2scan

subroutine f2c_destroy_scan( scan)

   USE s2hat_types
   USE s2hat_pixelization

   implicit none

   ! input
   type(scandef), pointer :: scan

   scan => globalScan

   call destroy_scan( globalScan)

   return

end subroutine f2c_destroy_scan

subroutine f2c_mpi_scanbcast( pixelization, scan, root, my_rank, current_comm)

   USE s2hat_types
   USE s2hat_pixelization

  implicit none

  ! input
  type( pixeltype), pointer :: pixelization
  integer(i4b) :: root, my_rank, current_comm

  ! output (and input on the root)
  type( scandef), pointer :: scan

  pixelization => globalPixelization
  scan => globalScan

  call mpi_scanbcast( globalPixelization, globalScan, root, my_rank, current_comm)

end subroutine f2c_mpi_scanbcast

subroutine f2c_mpi_pixelizationbcast( pixelization, root, my_rank, current_comm)

   USE s2hat_types
   USE s2hat_pixelization

  implicit none

  ! input
  integer(i4b) :: root, my_rank, current_comm

  type( pixeltype), pointer :: pixelization

  pixelization => globalPixelization

  call mpi_pixelizationbcast( globalPixelization, root, my_rank, current_comm)

end subroutine f2c_mpi_pixelizationbcast

function totalringsnumber( pixelization)
  USE healpix_types
  USE s2hat_types

  implicit none

  ! input
  type( pixeltype), pointer :: pixelization

  integer(i4b) :: totalringsnumber

  totalringsnumber = pixelization%nringsall

  return

end function totalringsnumber

subroutine pixelization2vectors( pixelization, ivec, ringInts, ringDPs)

  USE healpix_types
  USE s2hat_types

  implicit none

  ! input
  type( pixeltype), pointer :: pixelization
  ! output
  integer(i8b), dimension( 4) :: ivec
  integer(i8b), dimension(1:pixelization%nringsall, 1:2) ::  ringInts
  real(dp), dimension(1:pixelization%nringsall, 1:6) :: ringDPs

  integer(i4b) :: iring

  ivec(1) = pixelization%type
  ivec(2) = pixelization%npixsall
  ivec(3) = pixelization%nringsall
  ivec(4) = pixelization%nphmx

  do iring = 1, pixelization%nringsall

     ringInts( iring, 1) = pixelization%fpix(iring)
     ringInts( iring, 2) = pixelization%nph(iring)

     ringDPs( iring, 5) = pixelization%kphi(iring)
     ringDPs( iring, 6) = pixelization%qwght(iring)
     ringDPs( iring, 1) = pixelization%pixphi(iring)
     ringDPs( iring, 2) = pixelization%parea(iring)
     ringDPs( iring, 3) = pixelization%cth(iring)
     ringDPs( iring, 4) = pixelization%sth(iring)

  enddo

end subroutine pixelization2vectors

subroutine scan2vectors( pixelization, scan, ivec, ringFlags)

  USE healpix_types
  USE s2hat_types

  implicit none

  ! input
  type( pixeltype), pointer :: pixelization
  type( scandef), pointer :: scan

  ! output
  integer(i8b) :: ivec(2)
  integer(i8b), dimension( 1:pixelization%nringsall, 1:3) :: ringFlags

  integer(i4b) :: iring

  ivec(1) = scan%npixsobs
  ivec(2) = scan%nringsobs

  do iring = 1, pixelization%nringsall

     ringFlags( iring, 1) = scan%nfl( iring)
     ringFlags( iring, 2) = scan%sfl( iring)
     ringFlags( iring, 3) = scan%fl( iring)

  end do

end subroutine scan2vectors

subroutine vectors2pixelization( nringsall, ivec, ringInts, ringDPs, pixelization)

  USE healpix_types
  USE s2hat_types

  implicit none

  ! input
  integer(i8b) :: nringsall
  integer(i8b), dimension( 4) :: ivec
  integer(i8b), dimension(1:nringsall, 1:2) ::  ringInts
  real(dp), dimension(1:nringsall, 1:6) :: ringDPs

  ! output
  type( pixeltype), pointer :: pixelization

  integer(i4b) :: iring

  globalPixelization%type = ivec(1)
  globalPixelization%npixsall = ivec(2)
  globalPixelization%nringsall = ivec(3)
  globalPixelization%nphmx = ivec(4)

  allocate( globalPixelization%fpix(1:nringsall))
  allocate( globalPixelization%nph(1:nringsall))
  allocate( globalPixelization%kphi(1:nringsall))
  allocate( globalPixelization%qwght(1:nringsall))
  allocate( globalPixelization%pixphi(1:nringsall))
  allocate( globalPixelization%parea(1:nringsall))
  allocate( globalPixelization%cth(1:nringsall))
  allocate( globalPixelization%sth(1:nringsall))
  
  do iring = 1, globalPixelization%nringsall

     globalPixelization%fpix(iring) = ringInts( iring, 1)
     globalPixelization%nph(iring) = ringInts( iring, 2)

     globalPixelization%kphi(iring) = ringDPs( iring, 5)
     globalPixelization%qwght(iring) = ringDPs( iring, 6)

     globalPixelization%pixphi(iring) = ringDPs( iring, 1)
     globalPixelization%parea(iring) = ringDPs( iring, 2)
     globalPixelization%cth(iring) = ringDPs( iring, 3)
     globalPixelization%sth(iring) = ringDPs( iring, 4)

  enddo

  pixelization => globalPixelization

end subroutine vectors2pixelization

subroutine vectors2scan( nringsall, ivec, ringFlags, scan)

  USE healpix_types
  USE s2hat_types

  implicit none

  ! input
  integer(i8b) :: nringsall
  integer(i8b), dimension(2) :: ivec
  integer(i8b), dimension( 1:nringsall, 1:3) :: ringFlags

  ! output
  type( scandef), pointer :: scan

  !internal

  integer(i4b) :: iring

  globalScan%npixsobs = ivec(1)
  globalScan%nringsobs = ivec(2)

  allocate( globalScan%nfl(1:nringsall))
  allocate( globalScan%sfl(1:nringsall))
  allocate( globalScan%fl(1:nringsall))

  do iring = 1, nringsall

     globalScan%nfl( iring) = ringFlags( iring, 1)
     globalScan%sfl( iring) = ringFlags( iring, 2)
     globalScan%fl( iring)  = ringFlags( iring, 3)

  end do

  scan => globalScan

end subroutine vectors2scan

! data distribution/collection subroutine

subroutine f2c_distribute_local_data_objects_map2alm( precompute_plms, pixelization, scan, nlmax, nmmax, nmaps, nstokes, map, map_size, local_map, &
                       & nmvals, mvals, nplm, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, local_plm, first_ring, last_ring, local_w8ring, w8ring, myid, numprocs, root, comm)

    ! distributes and/or precomputes the data on all procs - a simple global "do-it-all" driver - (modified HEALPix routine)

    USE healpix_types
    USE s2hat_types
    USE s2hat_pixelization
    USE s2hat_toolbox_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b) :: precompute_plms, nmaps, nstokes, map_size, nlmax, nmmax, nmvals, first_ring, last_ring, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up
    integer(i4b) :: myid, numprocs, root, comm
    integer(i8b) :: nplm
    real(dp), dimension(0:pixelization%npixsall-1,1:nstokes,1:nmaps), target :: map
    real(dp), dimension(1:pixelization%nringsall,1:nstokes), target :: w8ring

    ! output
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(dp), dimension(0:map_size-1,1:nstokes,1:nmaps), target :: local_map
    real(dp), dimension(1:last_ring-first_ring+1,1:nstokes), target :: local_w8ring
    real(dp), dimension(ldim1_dwn:ldim1_up,ldim2_dwn:ldim2_up), target :: local_plm

    ! internal
    real(dp), dimension(:,:,:), pointer :: map_ptr
    real(dp), dimension(:,:), pointer :: w8ring_ptr
    real(dp), dimension(:,:,:), pointer :: local_map_ptr
    real(dp), dimension(:,:), pointer :: local_w8ring_ptr
    real(dp), dimension(:,:), pointer :: local_plm_ptr
	
    map_ptr => map
    w8ring_ptr => w8ring
    local_map_ptr => local_map
    local_w8ring_ptr => local_w8ring
    local_plm_ptr => local_plm

    pixelization => globalPixelization
    scan => globalScan
	
    call distribute_local_data_objects_map2alm( precompute_plms, globalPixelization, globalScan, nlmax, nmmax, nmaps, nstokes, map_ptr, map_size, local_map_ptr, &
                           & nmvals, mvals, ldim1_up, nplm, local_plm_ptr, first_ring, last_ring, local_w8ring_ptr, w8ring_ptr, myid, numprocs, root, comm) 

    return
	
end subroutine f2c_distribute_local_data_objects_map2alm

subroutine f2c_distribute_local_data_objects_alm2map( precompute_plms, pixelization, scan, nlmax, nmmax, nmaps, nstokes, ldim1_dwn,ldim1_up,ldim2_dwn,ldim2_up,ldim3_dwn,ldim3_up, &
                                                    & alms, nmvals, mvals, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up, local_alm, nplm, &
                                                    & ldim7_dwn, ldim7_up, ldim8_dwn, ldim8_up, local_plm, myid, numprocs, root, comm)

    ! distributes and/or precomputes the data on all procs - a simple global "do-it-all" driver - (modified HEALPix routine)

    USE healpix_types
    USE s2hat_types
    USE s2hat_pixelization
    USE s2hat_toolbox_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b) :: precompute_plms, nmaps, nstokes, nlmax, nmmax, nmvals
    integer(i4b) :: myid, numprocs, root, comm
    integer(i4b) :: ldim1_dwn,ldim1_up,ldim2_dwn,ldim2_up,ldim3_dwn,ldim3_up,ldim4_dwn,ldim4_up,ldim5_dwn,ldim5_up,ldim6_dwn,ldim6_up,ldim7_dwn,ldim7_up,ldim8_dwn,ldim8_up
    integer(i8b) :: nplm
    complex(DP),  dimension(ldim1_dwn:ldim1_up,ldim2_dwn:ldim2_up,ldim3_dwn:ldim3_up,1:nmaps), target :: alms

    !output
    complex(DP),  dimension(ldim4_dwn:ldim4_up,ldim5_dwn:ldim5_up,ldim6_dwn:ldim6_up,1:nmaps), target :: local_alm
    integer(i4b), dimension(0:nmvals-1), intent( out) :: mvals
    real(dp), dimension(ldim7_dwn:ldim7_up,ldim8_dwn:ldim8_up), target :: local_plm
	
    ! internal
    complex(dp), dimension(:,:,:,:), pointer :: alms_ptr
    complex(dp), dimension(:,:,:,:), pointer :: local_alm_ptr
    real(dp), dimension(:,:), pointer :: local_plm_ptr

    alms_ptr => alms
    local_alm_ptr => local_alm
    local_plm_ptr => local_plm

    pixelization => globalPixelization
    scan => globalScan
	
    call distribute_local_data_objects_alm2map( precompute_plms, globalPixelization, globalScan, nlmax, nmmax, nmaps, nstokes, ldim1_up, alms_ptr, nmvals, mvals, local_alm_ptr, &
                                             &  nplm, local_plm_ptr, myid, numprocs, root, comm)

    return
																																				
end subroutine f2c_distribute_local_data_objects_alm2map

subroutine f2c_get_local_data_sizes( precompute_plms, pixelization, scan, nlmax, nmmax, myid, numprocs, nmvals, &
                                 & first_ring, last_ring, map_size, nplm, root, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none
    
    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b) :: precompute_plms, nlmax, nmmax, myid, numprocs, root, comm
    ! output
    integer(i4b) :: first_ring, last_ring, map_size, nmvals
    integer(i8b) :: nplm

    pixelization => globalPixelization
    scan => globalScan

    call get_local_data_sizes( precompute_plms, globalPixelization, globalScan, nlmax, nmmax, myid, numprocs, nmvals, &
			     & first_ring, last_ring, map_size, nplm, root, comm)

    return

end subroutine f2c_get_local_data_sizes

subroutine f2c_find_mvalues( myid, numprocs, nmmax, nmvals, mvals)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none
  
    ! input
    integer(i4b) :: myid, numprocs, nmmax, nmvals
    ! output
    integer(i4b), dimension(0:nmvals-1) :: mvals
	
    call find_mvalues(myid, numprocs, nmmax, nmvals, mvals)

    return

end subroutine f2c_find_mvalues

function f2c_nummvalues( myid, numprocs, nmmax)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input 
    integer(i4b), intent(in) :: myid, numprocs, nmmax

    ! output 
    integer(i4b) :: f2c_nummvalues


    f2c_nummvalues = nummvalues( myid, numprocs, nmmax) 

end function f2c_nummvalues

function f2c_nummmodes( nlmax, nmvals, mvals)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input
    integer(i4b) :: nlmax, nmvals
    integer(i4b), dimension(0:nmvals-1) :: mvals
    ! output
    integer(i4b) :: f2c_nummmodes

    f2c_nummmodes = nummmodes( nlmax, nmvals, mvals)
	
    return

end function f2c_nummmodes

subroutine f2c_find_ring_range(pixelization, scan, nmmax, myid, numprocs, first_ring, last_ring, outerror)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

   implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b), intent(in) :: nmmax, myid, numprocs
    ! output
    integer(i4b), intent(out) :: first_ring, last_ring, outerror

    pixelization => globalPixelization
    scan => globalScan

    call find_ring_range( globalPixelization, globalScan, nmmax, myid, numprocs, first_ring, last_ring, outerror)

    return

end subroutine f2c_find_ring_range

subroutine f2c_distribute_map( pixelization, nmaps, mapnum, nstokes, first_ring, last_ring, map_size, local_map, map, myid, numprocs, root, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    integer(i4b) :: nmaps, mapnum, nstokes, numprocs, myid, root, comm
    integer(i4b) :: first_ring, last_ring, map_size
    real(DP), dimension(0:pixelization%npixsall-1,1:nstokes), target :: map
    ! output
    real(DP), dimension(0:map_size-1,1:nstokes,1:nmaps), target :: local_map

    ! internal
    integer(i4b) :: my_mapnum
    real(DP), dimension(:,:,:), pointer :: local_map_ptr
    real(DP), dimension(:,:), pointer :: map_ptr

    map_ptr => map
    local_map_ptr => local_map
    pixelization => globalPixelization

    my_mapnum = mapnum + 1

    call distribute_map( globalPixelization, nmaps, my_mapnum, nstokes, first_ring, last_ring, map_size, local_map_ptr, map_ptr, myid, numprocs, root, comm)

    return

end subroutine f2c_distribute_map

subroutine f2c_distribute_w8ring(npol, first_ring, last_ring, local_w8ring, nringsall, w8ring, myid, numprocs, root, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input
    integer(i4b) :: npol, numprocs, myid, root, comm
    integer(i4b) :: first_ring, last_ring, nringsall
    real(DP), dimension(1:nringsall,1:npol), target :: w8ring
    ! output
    real(DP), dimension(1:last_ring-first_ring+1,1:npol), target :: local_w8ring

    ! internal

    real(DP), dimension(:,:), pointer :: local_w8ring_ptr
    real(DP), dimension(:,:), pointer :: w8ring_ptr

    w8ring_ptr => w8ring
    local_w8ring_ptr => local_w8ring

    call distribute_w8ring(npol, first_ring, last_ring, local_w8ring_ptr, w8ring_ptr, myid, numprocs, root, comm)

    return

end subroutine f2c_distribute_w8ring

subroutine f2c_distribute_alms( nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, mvals, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, &
                              & local_alm, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up, alms, myid, numprocs, root, comm)

    ! distributes alms from proc root to all the others. It is m value-based distribution with the
    ! m values for each proc defined by mvals

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    !input
    integer(i4b) :: nmaps, mapnum, nlmax, nmmax, nstokes, nmvals, myid, numprocs, root, comm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    integer(i4b) :: ldim1_dwn,ldim1_up,ldim2_dwn,ldim2_up,ldim3_dwn,ldim3_up,ldim4_dwn,ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up
    complex(DP),  dimension(ldim4_dwn:ldim4_up,ldim5_dwn:ldim5_up,ldim6_dwn:ldim6_up),  target :: alms

    !output
    complex(DP),  dimension(ldim1_dwn:ldim1_up,ldim2_dwn:ldim2_up,ldim3_dwn:ldim3_up,1:nmaps), target :: local_alm
	
    ! internal
    integer(i4b) :: my_mapnum
    complex(DP),  dimension(:,:,:), pointer :: alms_ptr
    complex(DP),  dimension(:,:,:,:), pointer :: local_alm_ptr
	
    local_alm_ptr => local_alm
    alms_ptr => alms	
	
    my_mapnum = mapnum + 1

    call distribute_alms( nlmax, nmmax, nmaps, my_mapnum, nstokes, nmvals, mvals, ldim1_up, local_alm_ptr, alms_ptr, myid, numprocs, root, comm)

    return 
	
end subroutine f2c_distribute_alms

subroutine f2c_collect_map( pixelization, nmaps, mapnum, nstokes, map, first_ring, last_ring, map_size, local_map, myid, numprocs, root, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    integer(i4b) :: nmaps, mapnum, nstokes, numprocs, myid, root, comm
    integer(i4b) :: first_ring, last_ring, map_size
    real(DP), dimension(0:map_size-1,1:nstokes,1:nmaps), target :: local_map

    ! output
    real(DP), dimension(0:pixelization%npixsall-1,1:nstokes), target :: map

    ! internal

    integer(i4b) :: my_mapnum
    real(DP), dimension(:,:,:), pointer :: local_map_ptr
    real(DP), dimension(:,:), pointer :: map_ptr

    map_ptr => map
    local_map_ptr => local_map
    pixelization => globalPixelization

    my_mapnum = mapnum + 1

    call collect_map( globalPixelization, nmaps, my_mapnum, nstokes, map_ptr, first_ring, last_ring, map_size, local_map_ptr, myid, numprocs, root, comm)

    return

end subroutine f2c_collect_map

subroutine f2c_collect_alms( nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, mvals, ldim1_dwn,ldim1_up,ldim2_dwn,ldim2_up,ldim3_dwn,ldim3_up, &
                           & local_alm, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up, alms, myid, numprocs, root, comm)

    ! collects all the alms distributed over proc on a proc root

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    !input
    integer(i4b) :: nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, myid, numprocs, root, comm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    integer(i4b) :: ldim1_dwn,ldim1_up,ldim2_dwn,ldim2_up,ldim3_dwn,ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up
    complex(DP),  dimension(ldim1_dwn:ldim1_up,ldim2_dwn:ldim2_up,ldim3_dwn:ldim3_up,1:nmaps), target :: local_alm

    !output
    complex(DP),  dimension(ldim4_dwn:ldim4_up,ldim5_dwn:ldim5_up,ldim6_dwn:ldim6_up),  target :: alms

    ! internal
    integer(i4b) :: my_mapnum
    complex(DP),  dimension(:,:,:), pointer :: alms_ptr
    complex(DP),  dimension(:,:,:,:), pointer :: local_alm_ptr
	
    alms_ptr => alms
    local_alm_ptr => local_alm
	
    my_mapnum = mapnum + 1

    call collect_alms( nlmax, nmmax, nmaps, my_mapnum, nstokes, nmvals, mvals, ldim1_up, local_alm_ptr, alms_ptr, myid, numprocs, root, comm)

    return

end subroutine f2c_collect_alms

subroutine f2c_collect_cls( nmaps, mapnum, nstokes, nlmax, nmvals, mvals, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, local_alm, nspec, cl, myid, numprocs, root, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input
    integer(i4b) :: nmaps, mapnum, nstokes, numprocs, myid, root, comm
    integer(i4b) :: nspec, nlmax, nmvals, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up;
    integer(i4b), dimension(0:nmvals-1) :: mvals
    complex(DP), dimension( ldim1_dwn:ldim1_up, ldim2_dwn:ldim2_up, ldim3_dwn:ldim3_up, 1:nmaps), target :: local_alm

    ! output
    real(DP), dimension(0:nlmax,1:nspec) :: cl

    ! internal

    integer(i4b) :: my_mapnum
    complex(DP), dimension(:,:,:,:), pointer :: local_alm_ptr

    local_alm_ptr => local_alm

    my_mapnum = mapnum + 1

    call collect_cls( nmaps, my_mapnum, nstokes, nlmax, nmvals, mvals, ldim1_up, local_alm_ptr, nspec, cl, myid, numprocs, root, comm)

    return

end subroutine f2c_collect_cls

subroutine f2c_collect_xls( nmaps1, mapnum1, nmaps2, mapnum2, nstokes, nlmax, nmvals, mvals, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, &
                          & local_alm1, local_alm2, nspec, xcl, myid, numprocs, root, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input
    integer(i4b) :: nmaps1, nmaps2, mapnum1, mapnum2, nstokes, numprocs, myid, root, comm
    integer(i4b) :: nspec, nlmax, nmvals, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up;
    integer(i4b), dimension(0:nmvals-1) :: mvals
    complex(DP), dimension( ldim1_dwn:ldim1_up, ldim2_dwn:ldim2_up, ldim3_dwn:ldim3_up, 1:nmaps1), target :: local_alm1
    complex(DP), dimension( ldim1_dwn:ldim1_up, ldim2_dwn:ldim2_up, ldim3_dwn:ldim3_up, 1:nmaps2), target :: local_alm2

    ! output
    real(DP), dimension(0:nlmax,1:nspec) :: xcl

    ! internal

    integer(i4b) :: my_mapnum1, my_mapnum2
    complex(DP), dimension(:,:,:,:), pointer :: local_alm1_ptr, local_alm2_ptr

    local_alm1_ptr => local_alm1
    local_alm2_ptr => local_alm2	

    my_mapnum1 = mapnum1 + 1
    my_mapnum2 = mapnum2 + 1	

    call collect_xls( nmaps1, my_mapnum1, nmaps2, my_mapnum2, nstokes, nlmax, nmvals, mvals, ldim1_up, local_alm1_ptr, local_alm2_ptr, nspec, xcl, myid, numprocs, root, comm)

    return

end subroutine f2c_collect_xls

! transform routines

subroutine f2c_s2hat_alm2map( precompute_plms, pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring, &
                         & map_size, local_map, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, local_alm, &
                         & nplm, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, local_plm, numprocs, myid, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_alm2map_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b)  :: nlmax, nmmax, nmvals, nmaps, nstokes, first_ring, last_ring, map_size, numprocs, myid, comm, precompute_plms
    integer(i4b) :: ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up
    integer(i8b)  :: nplm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(kind=dp), dimension(ldim4_dwn:ldim4_up, ldim5_dwn:ldim5_up), target :: local_plm
    complex(kind=dp), dimension(ldim1_dwn:ldim1_up, ldim2_dwn:ldim2_up, ldim3_dwn:ldim3_up, 1:nmaps), target :: local_alm


    ! output
    real(kind=dp), dimension(0:map_size-1,1:nstokes,1:nmaps), target :: local_map

    ! internal

    real(kind=dp), dimension(:,:,:), pointer :: local_map_ptr
    real(kind=dp), dimension(:,:), pointer :: local_plm_ptr
    complex(kind=dp), dimension(:,:,:,:), pointer :: local_alm_ptr

    local_map_ptr => local_map
    local_plm_ptr => local_plm
    local_alm_ptr => local_alm
    pixelization => globalPixelization
    scan => globalScan

    call s2hat_alm2map( precompute_plms, globalPixelization, globalScan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring, &
                       & map_size, local_map_ptr, ldim1_up, local_alm_ptr, nplm, local_plm_ptr, numprocs, myid, comm)

    return

end subroutine f2c_s2hat_alm2map

subroutine f2c_s2hat_alm2map_spin( pixelization, scan, spin, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, &
                                 & map_size, local_map, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, local_alm, numprocs, myid, comm)
 
    USE healpix_types
    USE s2hat_types
    USE s2hat_alm2map_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b)  :: spin, nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, map_size, numprocs, myid, comm
    integer(i4b) :: ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up
    integer(i4b), dimension(0:nmvals-1) :: mvals
    complex(kind=dp), dimension(ldim1_dwn:ldim1_up, ldim2_dwn:ldim2_up, ldim3_dwn:ldim3_up, 1:nmaps), target :: local_alm

    ! output
    real(kind=dp), dimension(0:map_size-1,1:2,1:nmaps), target :: local_map

    ! internal

    real(kind=dp), dimension(:,:,:), pointer :: local_map_ptr
    complex(kind=dp), dimension(:,:,:,:), pointer :: local_alm_ptr

    local_map_ptr => local_map
    local_alm_ptr => local_alm
    pixelization => globalPixelization
    scan => globalScan

    call s2hat_alm2map_spin( globalPixelization, globalScan, spin, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, &
                           & map_size, local_map_ptr, ldim1_up, local_alm_ptr, numprocs, myid, comm)

    return

end subroutine f2c_s2hat_alm2map_spin

subroutine f2c_s2hat_map2alm( precompute_plms, pixelization, scan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring, local_w8ring, &
                           &  map_size, local_map, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, local_alm, nplm, ldim4_dwn, ldim4_up, &
                              ldim5_dwn, ldim5_up, local_plm, numprocs, myid, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_map2alm_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b) :: nlmax, nmmax, nmvals, nmaps, nstokes, first_ring, last_ring, map_size, numprocs, myid, comm, precompute_plms
    integer(i4b) :: ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up
    integer(i8b) :: nplm
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(kind=dp), dimension(0:map_size-1, 1:nstokes, 1:nmaps), target :: local_map
    real(kind=dp), dimension(ldim4_dwn:ldim4_up, ldim5_dwn:ldim5_up), target :: local_plm
    real(kind=dp), dimension(1:last_ring-first_ring+1, 1:nstokes), target :: local_w8ring

    ! output
    complex(kind=dp), dimension(ldim1_dwn:ldim1_up, ldim2_dwn:ldim2_up, ldim3_dwn:ldim3_up, 1:nmaps), target :: local_alm

    ! internal

    real(kind=dp), dimension(:,:), pointer :: local_w8ring_ptr   
    real(kind=dp), dimension(:,:,:), pointer :: local_map_ptr
    real(kind=dp), dimension(:,:), pointer :: local_plm_ptr
    complex(kind=dp), dimension(:,:,:,:), pointer :: local_alm_ptr

    local_map_ptr => local_map
    local_w8ring_ptr => local_w8ring
    local_plm_ptr => local_plm
    local_alm_ptr => local_alm

    pixelization => globalPixelization
    scan => globalScan

    call s2hat_map2alm( precompute_plms, globalPixelization, globalScan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring, local_w8ring_ptr, &
                       & map_size, local_map_ptr, ldim1_up, local_alm_ptr, nplm, local_plm_ptr, numprocs, myid, comm)

    return
  
end subroutine f2c_s2hat_map2alm

subroutine f2c_s2hat_map2alm_spin( pixelization, scan, spin, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, local_w8ring, &
                                &  map_size, local_map, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, local_alm, numprocs, myid, comm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_map2alm_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization
    type( scandef), pointer :: scan
    integer(i4b) :: spin, nlmax, nmmax, nmvals, nmaps, first_ring, last_ring, map_size, lda, numprocs, myid, comm
    integer(i4b) :: ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up
    integer(i4b), dimension(0:nmvals-1) :: mvals
    real(kind=dp), dimension(0:map_size-1, 1:2, 1:nmaps), target :: local_map
    real(kind=dp), dimension(1:last_ring-first_ring+1, 1:2), target :: local_w8ring

    ! output
    complex(kind=dp), dimension(ldim1_dwn:ldim1_up, ldim2_dwn:ldim2_up, ldim3_dwn:ldim3_up, 1:nmaps), target :: local_alm

    ! internal

    real(kind=dp), dimension(:,:), pointer :: local_w8ring_ptr   
    real(kind=dp), dimension(:,:,:), pointer :: local_map_ptr
    complex(kind=dp), dimension(:,:,:,:), pointer :: local_alm_ptr

    local_map_ptr => local_map
    local_w8ring_ptr => local_w8ring
    local_alm_ptr => local_alm

    pixelization => globalPixelization
    scan => globalScan

    call s2hat_map2alm_spin( globalPixelization, globalScan, spin, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring, local_w8ring_ptr, &
                           & map_size, local_map_ptr, ldim1_up, local_alm_ptr, numprocs, myid, comm)

    return
  
end subroutine f2c_s2hat_map2alm_spin

subroutine f2c_plm_mvalues_gen( pixelization, scan, npols, nlmax, nmmax, nmvals, mvals, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, nplm, local_plm)

    USE healpix_types
    USE s2hat_types
    USE s2hat_toolbox_mod

    implicit none

    ! input
    type( pixeltype), pointer :: pixelization  
    type( scandef), pointer :: scan 
    integer(i4b)   :: nlmax, nmmax, nmvals, npols, ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up
    integer(i4b), dimension( 0:nmvals-1) :: mvals
    integer(i8b) :: nplm
    !output
    real(dp), dimension(ldim1_dwn:ldim1_up,ldim2_dwn:ldim2_up), target :: local_plm
	
    !internal
    real(dp), dimension(:,:), pointer :: local_plm_ptr
	
    local_plm_ptr => local_plm

    pixelization => globalPixelization
    scan => globalScan
	
    call plm_mvalues_gen( globalPixelization,  globalScan,  npols, nlmax, nmmax, nmvals, mvals, ldim1_up, nplm, local_plm_ptr)

    return

end subroutine f2c_plm_mvalues_gen
