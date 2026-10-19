module comm_tod_objctr_mod
  use comm_tod_mod
  use comm_utils
  use comm_param_mod
  implicit none
  
  private
  public compute_obj_centric_signal, initialize_objctr_mod

  ! Define a comm_solsys_obj object here (or come up with a better name)
  ! Should contain things like name, reference flux density, phase (where relevant), albedo (?), splined {x,y,z} components in heliocentric coordinates as a function of time
  ! Perhaps a local map object in rectalinear coordinates?
  
  ! The Sun, Moon and Earth should be treated as special cases, as their properties are very different from distant planets, comets, asteroids etc.
  
  ! Put global variables for the ephemeris here


contains

  subroutine initialize_objctr_mod(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar

    ! Define a list of objects in comm_defs
    
    ! Define a new ephemeris ASCII filename in comm_param_mod
    
    ! Initialize ephemeris; read ephemeris file

    ! Spline each objects

  end subroutine initialize_objctr_mod

  ! Task: Compute predicted TOD for all objects listed in tod%active_objctr, except (possibly) the one in exclude_obj
  ! Put the result in sd%s_objctr
  ! By default, do it for all detectors. However, if det is present, there is only one, with the detector ID listed in det
  subroutine compute_obj_centric_signal(tod, sd, det, exclude_obj)
    implicit none
    class(comm_tod),      intent(in)             :: tod
    class(comm_scandata), intent(inout)          :: sd
    integer(i4b),         intent(in),   optional :: det
    character(len=*),     intent(in),   optional :: exclude_obj
    
    integer(i4b) :: comp, i, j, k, d, h, hp, ntod, nhorn, ndet, scan
    real(sp)     :: w, s
    real(dp)     :: t1, t2
    logical(lgt) :: incl_sun, incl_moon, incl_earth

    sd%s_objctr = 0.
    if (.not. allocated(tod%incl_objctr)) return

    scan  = sd%scan
    ntod  = sd%ntod
    ndet  = tod%ndet; if (present(det)) ndet = 1
    nhorn = tod%nhorn
        
    do comp = 1, size(tod%incl_objctr)
       if (present(exclude_obj)) then
          if (trim(tod%incl_objctr(comp)) == trim(exclude_obj)) cycle
       end if

       do j = 1, sd%ndet
          d = j; if (present(det)) d = det
          do h = 1, tod%nhorn
             hp = h; if (nhorn == 1) hp = 0
             do i = 1, ntod
                if (trim(tod%incl_objctr(comp)) == "sun") then
                   k = tod%scans(scan)%d(d)%pix_sol(i,h)
                   s = tod%map_solar(k,1)
                else if (trim(tod%incl_objctr(comp)) == "moon") then
                   k = tod%scans(scan)%d(d)%pix_moon(i,h)
                   s = tod%map_moon(k,1)
                else if (trim(tod%incl_objctr(comp)) == "earth") then
                   k = max(min(int(tod%scans(scan)%d(d)%earth_elon(i,h)/(pi/NBIN_EARTH_ELON)),NBIN_EARTH_ELON),1)
                   s = tod%map_earth(k)
                else
                   ! Solar system object; computing pointing on the fly relative to ephemeris; fill in predicted signal given current model
                   s = 0.
                end if
                
                if (s > -1.d30) sd%s_objctr(i,j,hp) = sd%s_objctr(i,j,hp) + s
             end do
          end do
       end do
    end do
    
  end subroutine compute_obj_centric_signal  
  
end module comm_tod_objctr_mod
