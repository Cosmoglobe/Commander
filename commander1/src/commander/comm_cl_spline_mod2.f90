!************************************************************************
!                                                                       *
!                                                                       *       
!       spline_mod.f90: fortran module to sample spline                 *
!                        coefficients for power spectra                 *
!                          Written by Eirik Gjerlow                     *       
!                                                                       *
!               Typically one calls initialize_cl_spline mod, then      *
!                sample_splined_spectrum, then cleanup_cl_spline_mod    *
!                                                                       *
!                                                                       *
!************************************************************************


module comm_cl_spline_mod
    use comm_utils
    !use corr_fileutils
    use math_tools
    use rngmod
    !use interpolation_mod
    use healpix_types
    use comm_mp_mod
    implicit none

    private

    public :: initialize_cl_spline_mod
    public :: sample_splined_spectrum
    public :: cleanup_cl_spline_mod

    type spline_data
        integer(i4b)                                        :: num_nodes
        integer(i4b), allocatable, dimension(:)             :: num_iter
        real(dp), allocatable, dimension(:)                 :: nodes
        real(dp), allocatable, dimension(:)                 :: propose_sd
    end type spline_data

    type spline_amps
        real(dp), allocatable, dimension(:)                 :: amps
    end type spline_amps

    real(dp),           allocatable,    dimension(:, :) :: bb_templ
    real(dp),           allocatable,    dimension(:)    :: bb_propose_sd
    integer(i4b),       allocatable,    dimension(:, :) :: no_acc_strikes
    integer(i4b),       allocatable,    dimension(:)    :: l_max, l_min
    integer(i4b),       allocatable,    dimension(:)    :: bb_num_iter
    integer(i4b)                                        :: l_max_max
    integer(i4b)                                        :: l_min_min
    integer(i4b)                                        :: l_max_min
    integer(i4b)                                        :: l_min_max
    integer(i4b)                                        :: num_iter_main
    integer(i4b)                                        :: num_modes
    integer(i4b)                                        :: num_modeel
    integer(i4b)                                        :: num_splinemodes
    integer(i4b)                                        :: num_bbtempl
    integer(i4b)                                        :: npix, nmaps
    real(dp)                                            :: tol
    logical(lgt)                                        :: firstcall
    logical(lgt)                                        :: autonodes
    logical(lgt)                                        :: bb
    type(spline_data), allocatable, dimension(:)        :: cl_spline_data

    contains

    subroutine initialize_cl_spline_mod(unit, paramfile)
        implicit none
        character(len=128)                              :: paramfile, nodefile
        character(len=128), allocatable, dimension(:)   :: bb_templfile
        integer(i4b)                                    :: unit, i, j, l, nside
        character(len=1)                                :: num_bbtempl_string
        real(dp)                                        :: dummy
        real(dp)                                        :: tempel
        logical(lgt)                                    :: provide_nodes, polarization

        call get_parameter(unit, paramfile, 'NSIDE',        par_int=nside)
        call get_parameter(unit, paramfile, 'POLARIZATION', par_lgt=polarization)
        nmaps = 1; if (polarization) nmaps = 3
        npix = 12*nside**2

        bb = .false.
        call get_parameter(unit, paramfile, 'NUM_MODES', par_int=num_modes)
        if (num_modes == 1) then 
            num_modeel = 1
            num_splinemodes = 1
        else if (num_modes == 2) then
            num_modeel = 3
            num_splinemodes = 3
        else if (num_modes == 3) then
            num_modeel = 4
            num_splinemodes = 3
            bb = .true.
        end if

        allocate(cl_spline_data(num_splinemodes))
        allocate(l_min(num_modeel))
        allocate(l_max(num_modeel))

        call get_parameter(unit, paramfile, 'LMAX_TT',       par_int=l_max(1))
        call get_parameter(unit, paramfile, 'LMIN_TT',       par_int=l_min(1))
        if (num_modes > 1) then
            call get_parameter(unit, paramfile, 'LMAX_TE',  par_int=l_max(2))
            call get_parameter(unit, paramfile, 'LMIN_TE',  par_int=l_min(2))
            call get_parameter(unit, paramfile, 'LMAX_EE',  par_int=l_max(3))
            call get_parameter(unit, paramfile, 'LMIN_EE',  par_int=l_min(3))

            if (l_max(2) > l_max(1) .or. l_max(2) > l_max(3)) then
                stop 'LMAX_TE > LMAX_TT or LMAX_EE. Stopping.'
            else if (l_min(2) < l_min(1) .or. l_min(2) < l_min(3)) then
                stop 'LMIN_TE < LMIN_TT or LMIN_EE. Stopping.'
            end if

            if (num_modes > 2) then
                call get_parameter(unit, paramfile, 'LMAX_BB', par_int=l_max(4))
                call get_parameter(unit, paramfile, 'LMIN_BB', par_int=l_min(4))
                call get_parameter(unit, paramfile, 'NUM_BBTEMPL', &
                   &  par_int=num_bbtempl)
                allocate(bb_templfile(num_bbtempl))
                allocate(bb_templ(l_min(4):l_max(4), num_bbtempl))
                allocate(bb_propose_sd(num_bbtempl))
                allocate(bb_num_iter(num_bbtempl))

                do i = 1, num_bbtempl
                    call int2string(i, num_bbtempl_string)
                    call get_parameter(unit, paramfile, 'BB_TEMPL_' & 
                        & // num_bbtempl_string, par_string=bb_templfile(i))
                end do
            end if
        end if

        call get_parameter(unit, paramfile, 'NUMITERMAIN', &
           &  par_int=num_iter_main)
        call get_parameter(unit, paramfile, 'AUTONODES', par_lgt=autonodes)
        if (autonodes) then
            call get_parameter(unit, paramfile, 'PROVIDE_NODES', &
               &  par_lgt=provide_nodes)
        else
            provide_nodes = .false.
        end if

        if (provide_nodes) then
            call get_parameter(unit, paramfile, 'NODEFILE', par_string=nodefile)
        end if

        call get_parameter(unit, paramfile, 'TOLERANCE',   par_dp=tol)

        l_max_max = maxval(l_max)
        l_min_min = minval(l_min)
        l_max_min = minval(l_max)
        l_min_max = maxval(l_min)

        if (provide_nodes) then
            open(unit, file=trim(nodefile), action='read')
            if (num_modeel == 1) then
                read(unit, *) cl_spline_data(1)%num_nodes
            else
                read(unit, *) cl_spline_data(1)%num_nodes, &
                   & cl_spline_data(2)%num_nodes, cl_spline_data(3)%num_nodes
            end if

            do i = 1, num_splinemodes
                allocate(cl_spline_data(i)%nodes(cl_spline_data(i)%num_nodes))
                allocate(cl_spline_data(i)%num_iter(cl_spline_data(i)%num_nodes))
                allocate(cl_spline_data(i)%propose_sd(cl_spline_data(i)%num_nodes))
            end do
            do i = 1, num_splinemodes
                do j = 1, cl_spline_data(i)%num_nodes
                    read(unit, *) cl_spline_data(i)%nodes(j)
                end do
            end do
            close(unit)
        else
            do i = 1, num_splinemodes
                cl_spline_data(i)%num_nodes = 2
                allocate(cl_spline_data(i)%nodes(2))
                allocate(cl_spline_data(i)%num_iter(2))
                allocate(cl_spline_data(i)%propose_sd(2))
                cl_spline_data(i)%nodes(1) = l_min(i)
                cl_spline_data(i)%nodes(2) = l_max(i)
            end do
        end if

        if (bb) then
            do i = 1, num_bbtempl
                open(unit, file = bb_templfile(i), action = 'read')
                do l = l_min(4), l_max(4)
                    read(unit, *) dummy, bb_templ(l, i)
                end do
                close(unit)
            end do
            bb_propose_sd = 1.d0
        end if

        firstcall = .true.
        do i = 1, num_splinemodes
            cl_spline_data(i)%propose_sd = 1.d0
        end do

        if (allocated(bb_templfile)) deallocate(bb_templfile)

    end subroutine initialize_cl_spline_mod

    !Nodenum can either be the node number (for TT, TE, EE) or template
    !number (for BB)
    subroutine draw_nodeamp(oldamp, amp, nodenum, modenum, rng_handle)
        implicit none
        real(dp), intent(in)            :: oldamp
        real(dp), intent(out)           :: amp
        integer(i4b), intent(in)        :: nodenum, modenum
        type(planck_rng)                :: rng_handle

        if (modenum <= 3) then
            amp  = oldamp + rand_gauss(rng_handle) &
                & * cl_spline_data(modenum)%propose_sd(nodenum)
        else
            amp = oldamp + rand_gauss(rng_handle)*bb_propose_sd(nodenum)
        end if

    end subroutine draw_nodeamp

    subroutine get_splined_cls(mode, amps, cl2, cl)
        implicit none
        real(dp), dimension(:), intent(in)   :: amps, cl2
        integer(i4b)                         :: mode
        real(dp), dimension(l_min(mode):), intent(out)  :: cl
        integer(i4b)                         :: l

        do l=l_min(mode), l_max(mode)
            cl(l) =  splint(cl_spline_data(mode)%nodes, amps, cl2, real(l, dp))
        end do

    end subroutine get_splined_cls

    subroutine adjust_nodes(cls, amps)
        implicit none
        real(dp),       dimension(l_min_min:, :), intent(in)    :: cls
        type(spline_amps), dimension(num_splinemodes),  intent(inout) :: amps
        real(dp),       allocatable, dimension(:, :)    :: splined_cls, rel_diff
        type(spline_amps), dimension(num_splinemodes)   :: cl2
        real(dp), allocatable, dimension(:)             :: tempnodes
        integer(i4b)                                    :: i, j, k, newnode
        real(dp)                                        :: lp1, lpn
        logical(lgt), dimension(num_splinemodes)        :: ok
        logical(lgt)                                    :: okprior
        character(len=128)                              :: filename

        lp1 = 1.d30
        lpn = 1.d30

        allocate(splined_cls(l_min_min:l_max_max, num_splinemodes))
        allocate(rel_diff(l_min_min:l_max_max, num_splinemodes))

        splined_cls = 0.d0
        rel_diff = 0.d0

        do i = 1, num_splinemodes
            allocate(cl2(i)%amps(cl_spline_data(i)%num_nodes))
            call spline(cl_spline_data(i)%nodes, amps(i)%amps, lp1, lpn, &
               & cl2(i)%amps)
            call get_splined_cls(i, amps(i)%amps, cl2(i)%amps, &
               & splined_cls(l_min(i):l_max(i), i))
        end do

        call check_spline(splined_cls, cls, ok, rel_diff)

        do while(any(.not. ok))
            do i = 1, num_splinemodes
                if (.not. ok(i)) then
                    newnode = maxloc(abs(rel_diff(l_min(i):l_max(i), i)), 1) &
                        &  + l_min(i) - 1
                    call enforce_deltal_prior(newnode, i, okprior)
                    do while (.not. okprior)
                        rel_diff(newnode, i) = 0.d0
                        newnode = maxloc(abs(rel_diff(l_min(i):l_max(i), i)), &
                            & 1)  + l_min(i) - 1
                        call enforce_deltal_prior(newnode, i, okprior)
                    end do

                    allocate(tempnodes(cl_spline_data(i)%num_nodes))
                    cl_spline_data(i)%num_nodes = cl_spline_data(i)%num_nodes &
                       &  + 1
                    tempnodes = cl_spline_data(i)%nodes
                    deallocate(cl_spline_data(i)%nodes)
                    deallocate(amps(i)%amps)
                    deallocate(cl2(i)%amps)
                    allocate(cl_spline_data(i)%nodes(cl_spline_data(i)%num_nodes))
                    allocate(amps(i)%amps(cl_spline_data(i)%num_nodes))
                    allocate(cl2(i)%amps(cl_spline_data(i)%num_nodes))
                    k = 1
                    cl_spline_data(i)%nodes = 0.d0
                    do j = 1, cl_spline_data(i)%num_nodes
                        if (k == cl_spline_data(i)%num_nodes) then
                            cl_spline_data(i)%nodes(j) = real(newnode, dp)
                        else if (k == j .and. newnode < tempnodes(k)) then
                            cl_spline_data(i)%nodes(j) = real(newnode, dp)
                        else
                            cl_spline_data(i)%nodes(j) = tempnodes(k)
                            k = k + 1
                        end if
                    end do

                    amps(i)%amps = cls(int(cl_spline_data(i)%nodes), i)

                    call spline(cl_spline_data(i)%nodes, amps(i)%amps, &
                        & lp1, lpn, cl2(i)%amps)
                    call get_splined_cls(i, amps(i)%amps, cl2(i)%amps, &
                        & splined_cls(l_min(i):l_max(i), i))
                    if (isnan_cls(splined_cls, num_splinemodes)) then
                        filename = 'nan_splinecls.dat'
                        call print_cls(splined_cls, filename, num_splinemodes)
                        filename = 'nan_splineamps.dat'
                        call print_amps(amps, filename)
                        stop 'NaN in splined cls. Stopping'
                    end if
                    deallocate(tempnodes)
                end if
            end do
            call check_spline(splined_cls, cls, ok, rel_diff)
        end do

        deallocate(rel_diff)
        deallocate(splined_cls)
        do i = 1, num_splinemodes
            deallocate(cl2(i)%amps)
        end do
        
    end subroutine adjust_nodes

    !Now enforces the exclusion of l=3. Uncomment below to include 
    subroutine enforce_deltal_prior(node, mode, ok)
        integer(i4b), intent(inout)     :: node           
        integer(i4b), intent(in)        :: mode
        integer(i4b)                    :: savednode
        logical(lgt)                    :: ok
        real(dp)                        :: deltalogl, loglength
        logical(lgt)                    :: rising, falling
        integer(i4b)                    :: i, j

        ok = .true.
        if (node == 3) then
            do j = 1, cl_spline_data(mode)%num_nodes
                if (cl_spline_data(mode)%nodes(j) == 4) then
                    ok = .false.
                    return
                end if
            end do
            node = 4
        end if

        do j = 1, cl_spline_data(mode)%num_nodes
            if (cl_spline_data(mode)%nodes(j) == node) then
                stop 'Tolerance for relative difference is too low. Stopping.'
            end if
        end do

        !loglength = 0.d0
        !savednode = node
        !ok = .false.
        !rising = .false.
        !falling = .false.

        !do j = 1, cl_spline_data(mode)%num_nodes
        !    if (node == cl_spline_data(mode)%nodes(j)) then
        !        stop 'Proposed node is equal to an already existing node. Either your tolerance is too low or the spacing between nodes is set too high.'
        !    end if
        !end do

        !do while (.not. ok)
        !    ok = .true.
        !    do j = 1, cl_spline_data(mode)%num_nodes
        !        deltalogl = log(real(node, dp)) - &
        !           & log(cl_spline_data(mode)%nodes(j))
        !        if (deltalogl > 0) then
        !            if (.not. (deltalogl >= loglength)) then
        !                ok = .false.
        !                rising = .true.
        !                node = node + 1   
        !            end if
        !        else
        !            if (.not. (deltalogl <= -loglength)) then
        !                ok = .false.
        !                falling = .true.
        !                node = node - 1
        !            end if
        !        end if

        !        if (node == l_min(mode) .or. node == l_max(mode)) then
        !            ok = .false.
        !            node = savednode
        !            return
        !        end if

        !        if (rising .and. falling) then
        !            ok = .false.
        !            node = savednode
        !            return
        !        end if
        !    end do
        !end do
    end subroutine enforce_deltal_prior

    function isnan_cls(cls, nmodes)
        implicit none
        real(dp), dimension(l_min_min:, :)      :: cls
        integer(i4b), optional                  :: nmodes
        integer(i4b)                            :: whichmodes
        logical(lgt)                            :: isnan_cls
        integer(i4b)                            :: i, j, l

        if (present(nmodes)) then
            whichmodes = nmodes
        else
            whichmodes = num_modeel
        end if

        isnan_cls = .false.

        do i = 1, whichmodes
            do l = l_min(i), l_max(i)
               ! HKE: Removed isnan for now, since pg compiler doesn't like it
                if (.false.) then
!                if (isnan(cls(l, i))) then
                    isnan_cls = .true.
                end if
            end do
        end do

    end function isnan_cls

    !An initial chain to adjust proposal density
    subroutine adjust_proposal_sd(amps, s, bb_amps)
        implicit none
        type(genvec),                                   intent(inout) :: s
        type(spline_amps),  dimension(num_splinemodes), intent(inout) :: amps
        real(dp),       dimension(:), optional, intent(inout)         :: bb_amps
        type(spline_amps), dimension(num_splinemodes)              :: cl2
        real(dp),       allocatable, dimension(:, :)    :: currcl, newcl, cl
        real(dp)                                        :: loglike
        real(dp)                                        :: currloglike
        real(dp)                                        :: newloglike
        real(dp)                                        :: savedval
        real(dp)                                        :: currval
        real(dp)                                        :: propval
        real(dp)                                        :: acc_rat
        real(dp)                                        :: acc, prop
        real(dp)                                        :: lp1, lpn
        integer(i4b), parameter                         :: numsamp = 100
        real(dp), dimension(numsamp)                    :: chain
        real(dp)                                        :: sd
        character(len=128)                              :: filename
        type(planck_rng)                                :: rng_handle

        integer(i4b)                                    :: i, j, k, reps

        call rand_init(rng_handle, 12837)

        do i = 1, num_splinemodes
            allocate(cl2(i)%amps(cl_spline_data(i)%num_nodes))
        end do
        if (bb) then
            if (.not. present(bb_amps)) then
                stop 'Error (adjust_proposal_sd): num_modes = 3 but  bb_amps not present. '
            end if
        end if
        allocate(currcl(l_min_min:l_max_max, num_modeel))
        allocate(newcl(l_min_min:l_max_max, num_modeel))
        allocate(cl(l_min_min:l_max_max, num_modeel))

        currcl = 0.d0
        newcl = 0.d0
        cl = 0.d0

        lp1 = 1.d30
        lpn = 1.d30
        
        do i = 1, num_splinemodes
            call spline(cl_spline_data(i)%nodes, amps(i)%amps, lp1, lpn, &
               & cl2(i)%amps)
            call get_splined_cls(i, amps(i)%amps, cl2(i)%amps, & 
               & cl(l_min(i):l_max(i), i))
        end do

        if (bb) then
            call get_bbcl(bb_amps, cl(l_min(4):l_max(4), 4))
        end if

        loglike = eval_loglike(s, cl, cl)

        if (loglike == -1.d30) then
            filename = 'zerocls.dat'
            call print_cls(cl, filename)
            stop 'Initial likelihood (in adjust_proposal_sd) is zero'
        end if

        currcl = cl
        currloglike =  loglike

        do i = 1, num_modeel
            print *, 'Mode: ', i
            if (i <= 3) then
                do j = 1, cl_spline_data(i)%num_nodes
                    print *, 'Node: ', j
                    acc_rat = 1.d0
                    reps = 0
                    savedval = amps(i)%amps(j)
                    do while (acc_rat > 0.8 .or. acc_rat < 0.2)
                        reps = reps + 1
                        currcl = cl
                        newcl = cl
                        acc = 0.d0
                        prop = 0.d0
                        currloglike = loglike
                        currval = savedval
                        do k = 1, numsamp
                            prop = prop + 1.d0
                            call draw_nodeamp(amps(i)%amps(j), propval, j, i, rng_handle)
                            amps(i)%amps(j) = propval
                            call spline(cl_spline_data(i)%nodes, &
                               & amps(i)%amps, lp1, lpn, cl2(i)%amps)
                            call get_splined_cls(i, amps(i)%amps, &
                               & cl2(i)%amps, newcl(l_min(i):l_max(i), i))
                            newloglike = eval_loglike(s, currcl, newcl)
                            if  (exp(newloglike-currloglike) > &
                               &  rand_uni(rng_handle) .and. &
                               & newloglike .ne. -1.d30) then
                                currloglike = newloglike
                                call scale_s(s, currcl, newcl)
                                currcl = newcl
                                acc = acc + 1.d0
                                currval = propval
                            end if
                            amps(i)%amps(j) = currval
                            chain(k) = amps(i)%amps(j)
                        end do
                        acc_rat = acc/prop
                        if ((reps == 3 .or. reps == 6) .and. &
                        !if (mod(reps, 3) == 0 .and. &
                            &(acc_rat < 0.2 .or. acc_rat > 0.8)) then
                            sd = get_sd(chain)
                            if (sd .ne. 0.d0) then
                                cl_spline_data(i)%propose_sd(j) = sd
                            end if
                        else
                            if (acc_rat < 0.2) then
                                cl_spline_data(i)%propose_sd(j) = cl_spline_data(i)%propose_sd(j)*0.5d0
                            else if (acc_rat > 0.8) then
                                cl_spline_data(i)%propose_sd(j) = cl_spline_data(i)%propose_sd(j)*2.d0
                            end if
                        end if
                        amps(i)%amps(j) = savedval
                    end do
                    write(*,*) 'Accept rate = ', acc_rat
                end do
            else
                do j = 1, num_bbtempl
                    acc_rat = 1.d0
                    reps = 0
                    savedval = bb_amps(j)
                    do while (acc_rat > 0.8 .or. acc_rat < 0.2)
                        reps = reps + 1
                        currloglike = loglike
                        currcl = cl
                        newcl = cl
                        acc = 0.d0
                        prop = 0.d0
                        currval = savedval
                        do k = 1, numsamp
                            prop = prop + 1.d0
                            call draw_nodeamp(bb_amps(j), propval, &
                               & j,i, rng_handle)
                            bb_amps(j) = propval
                            call get_bbcl(bb_amps, & 
                                & newcl(l_min(i):l_max(i), i))
                            newloglike = eval_loglike(s, currcl, newcl)
                            if  (exp(newloglike-currloglike) > &
                               &  rand_uni(rng_handle) .and. &
                               & newloglike .ne. -1.d30) then
                                currloglike = newloglike
                                currval = propval
                                call scale_s(s, currcl, newcl)
                                currcl = newcl
                                acc = acc + 1.d0
                            end if
                            bb_amps(j) = currval
                            chain(k) = bb_amps(j)
                        end do
                        acc_rat = acc/prop
                        !if (mod(reps, 3) == 0 .and. &
                        if ((reps == 3 .or. reps == 6) .and. &
                            &(acc_rat < 0.2 .or. acc_rat > 0.8)) then
                            sd = get_sd(chain)
                            if (sd .ne. 0.d0) then
                                bb_propose_sd(j) = sd
                            end if
                        else
                            if (acc_rat < 0.2) then
                                bb_propose_sd(j) = bb_propose_sd(j)*0.5d0
                            else if (acc_rat > 0.8) then
                                bb_propose_sd(j) = bb_propose_sd(j)*2.d0
                            end if
                        end if
                        bb_amps(j) = savedval
                    end do
                end do
            end if
        end do

        do i = 1, num_splinemodes
            deallocate(cl2(i)%amps)
        end do
        deallocate(currcl)
        deallocate(newcl)

    end subroutine adjust_proposal_sd

    subroutine get_bbcl(amps, cl)
        implicit none
        real(dp), dimension(num_bbtempl), intent(in)    :: amps
        real(dp), dimension(l_min(4):), intent(out)     :: cl
        integer(i4b)                                    :: i

        cl(l_min(4):l_max(4)) = 0.d0
        do i = 1, num_bbtempl
            cl = cl + amps(i)*bb_templ(:, i)
        end do

    end subroutine get_bbcl
    
    function get_sd(chain)
        implicit none
        real(dp)                        :: get_sd, var, mean
        real(dp),       dimension(:)    :: chain
        integer(i4b)                    :: i, length

        length = size(chain)
        var = 0.d0

        mean = sum(chain)/max(1.d0, real(length, dp))
        do i = 1, length
            var = var + chain(i)*chain(i)
        end do
        var = var/max(1.d0, real(length, dp)) - mean*mean
        if (var < 0.d0) then
            get_sd = 0.d0
            return
        end if

        get_sd = sqrt(var)

    end function get_sd

    !auxilliary routine to dump cls to file
    subroutine print_cls(cls, filename, modes)
        implicit none
        real(dp),       dimension(l_min_min:, :)        :: cls
        integer(i4b), optional                          :: modes
        integer(i4b)                                    :: unit, l, howmany
        character(len=128)                              :: filename

        unit = 87

        if (present(modes)) then
            howmany = modes
        else
            howmany = num_modeel
        end if

        open(unit, file=trim(filename), recl=512)
        if (howmany == 1) then
            do l = l_min_min, l_max_max
                write(unit, *) l, cls(l, 1)
            end do
        else if (howmany == 3) then
            do l = l_min_min, l_max_max
                write(unit, *) l, cls(l, 1), cls(l, 2), cls(l, 3), cls(l,1)*cls(l,3)-cls(l,2)**2
            end do
        else if (howmany == 4) then
            do l = l_min_min, l_max_max
                write(unit, *) l, cls(l, 1), cls(l, 2), cls(l, 3), cls(l, 4)
            end do
        end if

        close(unit)

    end subroutine print_cls

    subroutine print_amps(amps, filename, bbamps)
        implicit none
        type(spline_amps),       dimension(num_splinemodes)     :: amps
        real(dp),       dimension(num_bbtempl), optional        :: bbamps
        integer(i4b)                                            :: unit, i, l, j
        character(len=128), optional                    :: filename
        character(len=128)                              :: filename1, filename2

        if (present(filename)) then
            filename1 = filename
        else 
            filename1 = 'dumpamp.log'
        end if
        filename2 = 'dumpnodes.log'
        unit = 88

        open(unit, file=trim(filename1))
        do i = 1, num_splinemodes
            write(unit, *) cl_spline_data(i)%num_nodes
            do j = 1, cl_spline_data(i)%num_nodes
                write(unit, *) amps(i)%amps(j)
            end do
        end do

        if (bb .and. present(bbamps)) then
            do j = 1, num_bbtempl
                write(unit, *) bbamps(j)
            end do
        end if

        close(unit)

        open(unit, file=trim(filename2))
        do i = 1, num_splinemodes
            write(unit, *) cl_spline_data(i)%num_nodes
            do j = 1, cl_spline_data(i)%num_nodes
                write(unit, *), cl_spline_data(i)%nodes(j)
            end do
        end do
        close(unit)
    end subroutine print_amps

    !Routine to adjust number of iterations to reduce correlation length
    subroutine adjust_iternum(amps, s, bb_amps)
        implicit none
        type(spline_amps),  dimension(num_splinemodes), intent(inout)    :: amps
        type(genvec),                                   intent(inout)    :: s
        real(dp),       dimension(num_bbtempl), optional        :: bb_amps
        real(dp),       allocatable, dimension(:, :)    :: currcl, newcl, cl
        type(spline_amps),  dimension(num_splinemodes)  :: cl2
        real(dp),       allocatable, dimension(:)       :: chain
        real(dp)                                        :: lp1, lpn
        real(dp)                                        :: currloglike, newloglike
        real(dp)                                        :: loglike
        real(dp)                                        :: savedval
        real(dp)                                        :: currval
        real(dp)                                        :: propval
        logical(lgt)                                    :: ok

        integer(i4b)                                    :: i, j, k
        integer(i4b)                                    :: chainlength
        type(planck_rng)                                :: rng_handle

        call rand_init(rng_handle, 12837)

        allocate(currcl(l_min_min:l_max_max, num_modeel))
        allocate(newcl(l_min_min:l_max_max, num_modeel))
        allocate(cl(l_min_min:l_max_max, num_modeel))
        do i = 1, num_splinemodes
            allocate(cl2(i)%amps(cl_spline_data(i)%num_nodes))
        end do

        currcl = 0.d0
        newcl = 0.d0
        cl = 0.d0

        lp1 = 1.d30
        lpn = 1.d30
        
        do i = 1, num_splinemodes
            call spline(cl_spline_data(i)%nodes, amps(i)%amps, lp1, lpn, cl2(i)%amps)
            call get_splined_cls(i, amps(i)%amps, cl2(i)%amps, cl(l_min(i):l_max(i), i))
        end do

        if (bb) then
            call get_bbcl(bb_amps, cl(l_min(4):l_max(4), 4))
        end if

        loglike =  eval_loglike(s, cl, cl)

        if (loglike == -1.d30) then
            stop 'Initial likelihood is 0 (adjust_iternum)'
        end if

        currloglike = loglike
        currcl = cl

        do i = 1, num_modeel
            print *, 'Mode: ', i
            if (i <= 3) then
                do j = 1, cl_spline_data(i)%num_nodes
                    print *, 'Node: ', j
                    ok = .false.
                    savedval = amps(i)%amps(j)
                    chainlength = 100
                    allocate(chain(chainlength))
                    do while (.not. ok)
                        currval = savedval
                        currcl = cl
                        newcl = cl
                        currloglike = loglike
                        do k = 1, chainlength
                            call draw_nodeamp(amps(i)%amps(j), & 
                               & propval, j, i, rng_handle)
                            amps(i)%amps(j) = propval
                            call spline(cl_spline_data(i)%nodes, &
                               & amps(i)%amps, lp1, lpn, cl2(i)%amps)
                            call get_splined_cls(i, amps(i)%amps, &
                               & cl2(i)%amps, newcl(l_min(i):l_max(i), i))
                            newloglike = eval_loglike(s, currcl, newcl)
                            if  (exp(newloglike-currloglike) > & 
                               & rand_uni(rng_handle) .and. newloglike .ne. -1.d30) &
                               & then
                                currval = propval
                                currloglike = newloglike
                                call scale_s(s, currcl, newcl)
                                currcl = newcl
                            end if
                            amps(i)%amps(j) = currval
                            chain(k) = amps(i)%amps(j)
                        end do
                        call set_iternum(chain, j, i, ok)
                        if (.not. ok) then
                            chainlength = chainlength*10
                            deallocate(chain)
                            allocate(chain(chainlength))
                        end if
                        amps(i)%amps(j) = savedval
                    end do
                    deallocate(chain)
                end do
            else
                do j = 1, num_bbtempl
                    ok = .false.
                    chainlength = 100
                    allocate(chain(chainlength))
                    savedval = bb_amps(j)
                    do while (.not. ok)
                        currcl = cl
                        newcl = cl
                        currloglike = loglike
                        currval = savedval
                        do k = 1, chainlength
                            call draw_nodeamp(bb_amps(j), & 
                                & propval, j, i, rng_handle)
                            bb_amps(j) = propval
                            call get_bbcl(bb_amps, & 
                                & newcl(l_min(i):l_max(i), i))
                            newloglike = eval_loglike(s, currcl, newcl)
                            if  (exp(newloglike-currloglike) > & 
                               & rand_uni(rng_handle) .and. newloglike .ne. -1.d30) &
                               & then
                                currloglike = newloglike
                                currval = propval
                                call scale_s(s, currcl, newcl)
                                currcl = newcl
                            end if
                            bb_amps(j) = currval
                            chain(k) = bb_amps(j)
                        end do
                        call set_iternum(chain, j, i, ok)
                        if (.not. ok) then
                            chainlength = chainlength*10
                            deallocate(chain)
                            allocate(chain(chainlength))
                        end if
                        bb_amps(j) = savedval
                    end do
                    deallocate(chain)
                end do
            end if
        end do

        deallocate(currcl)
        deallocate(newcl)
        do i = 1, num_splinemodes
            deallocate(cl2(i)%amps)
        end do

    end subroutine adjust_iternum

    !node_num can either be node number (for TT, TE and EE) or template number
    !(for BB)
    subroutine set_iternum(chain, node_num, mode_num, ok)
        implicit none
        real(dp),       dimension(:), intent(in)            :: chain
        real(dp), allocatable,  dimension(:)    :: chainminmean
        integer(i4b), intent(in)                :: node_num, mode_num
        integer(i4b)                            :: chainlength, corrsize
        integer(i4b)                            :: i, j, k
        real(dp), allocatable, dimension(:)     :: corr
        real(dp)                                :: auto, mean
        logical(lgt), intent(out)               :: ok

        ok = .true.

        chainlength = size(chain)
        mean = sum(chain)/max(1, chainlength)
        corrsize = chainlength/2
        allocate(corr(corrsize))
        allocate(chainminmean(chainlength))

        chainminmean = chain - mean

        auto = 0.d0
        corr = 0.d0

        do i = 1, chainlength
            auto = auto + (chainminmean(i))*(chainminmean(i))
            do j = 1, corrsize
                if ((i + j) <= chainlength) then
                    k = i + j
                else
                    k = i + j - chainlength
                end if
                corr(j) = corr(j) + (chainminmean(i))*(chainminmean(k))
            end do
        end do

        corr = corr/auto

        i = 1
        do while (corr(i) > 0.1)
            i = i + 1
            if (i > corrsize) then
                print *, 'Warning: correlation length is greater than ', corrsize
                ok = .false.
                return
            end if
        end do

        if (mode_num <= 3) then
            cl_spline_data(mode_num)%num_iter(node_num) = i
        else
            bb_num_iter(node_num) = i
        end if

        deallocate(corr)
        deallocate(chainminmean)
        
    end subroutine set_iternum

    function get_detcl(l, cl, submat)
        implicit none
        real(dp)                                        :: get_detcl
        real(dp), dimension(:), intent(in)              :: cl
        integer(i4b), optional, intent(in)              :: submat
        integer(i4b)                                    :: mat
        integer(i4b), intent(in)                        :: l

        if (.not. present(submat)) then
            mat = num_modes
        else
            mat = submat
        end if


        if (num_modes == 1) then
            get_detcl = cl(1)
        else if (num_modes == 2 .and. mat == num_modes) then
            if (l >= l_min_max .and. l <= l_max_min) then
                get_detcl = cl(1)*cl(3)-cl(2)*cl(2)
            else
                if (in_range(l, 1) .and. (.not. in_range(l, 3))) then
                    get_detcl = cl(1)
                else
                    get_detcl = cl(3)
                end if
            end if
        else if (num_modes == 2 .and. mat == 1) then
            if (in_range(l, 1)) then
                get_detcl = cl(1)
            else
                get_detcl = cl(3)
            end if
        else if (num_modes == 3 .and. mat == num_modes) then
            if (l >= l_min_max .and. l <= l_max_min) then
                get_detcl = cl(4)*(cl(1)*cl(3) - cl(2)*cl(2))
            else
                if (in_range(l, 1) .and. (.not. in_range(l, 3)) .and. &
                    & (.not. in_range(l, 4))) then
                    get_detcl = cl(1)
                else if (in_range(l, 1) .and. in_range(l, 3) .and. &
                    & (.not. in_range(l, 4))) then
                    get_detcl = cl(1)*cl(3)-cl(2)*cl(2)
                else if (in_range(l, 1) .and. (.not. in_range(l, 3)) .and. &
                    & in_range(l, 4)) then
                    get_detcl = cl(1)*cl(4)
                else if ((.not. in_range(l, 1)) .and. in_range(l, 3) .and. &
                    & in_range(l, 4)) then
                    get_detcl = cl(3)*cl(4)
                else if ((.not. in_range(l, 1)) .and. (.not. in_range(l, 3)) &
                    & .and. in_range(l, 4)) then
                    get_detcl = cl(4)
                else
                    get_detcl = cl(3)
                end if
            end if
        else if (num_modes == 3 .and. mat == 2) then
            if ((l >= l_min_max .and. l <= l_max_min) .or. (in_range(l, 1) &
                & .and. in_range(l, 2))) then
                get_detcl = cl(1)*cl(3)-cl(2)*cl(2)
            else if (in_range(l, 1) .and. (.not. in_range(l, 2)) &
                & .and. in_range(l, 3)) then
                get_detcl = cl(1)*cl(4)
            else if (in_range(l, 1) .and. (.not. in_range(l, 2)) .and. &
                & (.not. in_range(l, 3))) then
                get_detcl = cl(1)
            else if ((.not. in_range(l, 1)) .and. in_range(l, 2) .and. &
                & in_range(l, 3)) then
                get_detcl = cl(3)*cl(4)
            else if ((.not. in_range(l, 1)) .and. in_range(l, 2) .and. &
                & (.not. in_range(l, 3))) then
                get_detcl = cl(3)
            else if ((.not. in_range(l, 1)) .and. (.not. in_range(l, 2)) .and. &
                & in_range(l, 3)) then
                get_detcl = cl(4)
            end if
        else if (num_modes == 3 .and. mat == 1) then
            if (in_range(l, 1)) then
                get_detcl = cl(1)
            else if (in_range(l, 2)) then
                get_detcl = cl(2)
            else if (in_range(l, 3)) then
                get_detcl = cl(3)
            end if
        end if

    end function get_detcl

    function is_pos_def(l, cl)
        implicit none
        integer(i4b)                    :: l
        real(dp), dimension(:)          :: cl
        logical(lgt)                    :: is_pos_def

        is_pos_def = .false.

        if (num_modes == 1) then
            if (get_detcl(l, cl) > 0) then
                is_pos_def = .true.
            end if
        else if (num_modes == 2) then
            if (get_detcl(l, cl) > 0 .and. get_detcl(l, cl, 1) > 0) then
                is_pos_def = .true.
            end if
        else if (num_modes == 3) then
            if (get_detcl(l, cl) > 0 .and. get_detcl(l, cl, 2) > 0 .and. &
                & get_detcl(l, cl, 1) > 0 ) then
                is_pos_def = .true.
            end if
        end if

    end function is_pos_def

    function in_range(l, mode)
        implicit none
        logical(lgt)            :: in_range
        integer(i4b)            :: l, mode

        if (l >= l_min(mode) .and. l <= l_max(mode)) then
            in_range = .true.
        else
            in_range = .false.
        end if
    end function in_range

    function eval_loglike(s, cl_old, cl_new)
        implicit none
        real(dp), dimension(l_min_min:, :), intent(in)     :: cl_old, cl_new
        type(genvec),                       intent(inout)  :: s
        real(dp)                                        :: loglike
        real(dp)                                        :: trinv
        real(dp)                                        :: detcl, const
        real(dp)                                        :: eval_loglike
        integer(i4b)                                    :: l, i

        type(genvec) :: s_scaled
        integer(i4b) :: ierr
        real(dp),     allocatable, dimension(:,:)     :: chisq_map
        
        allocate(chisq_map(0:npix-1,1))
        call allocate_genvec(s_scaled)
        call genvec_set_equal(s, s_scaled)
        call scale_s(s_scaled, cl_old, cl_new, ierr)
        if (ierr == 0) then
           call compute_chisq(map_id, s_scaled, chisq_map=chisq_map)
           eval_loglike = -0.5d0*sum(chisq_map)
        else
           eval_loglike = -1.d30
        end if
        call deallocate_genvec(s_scaled)
        deallocate(chisq_map)

    end function eval_loglike

    subroutine scale_s(s, cl_old, cl_new, ierr)
      implicit none

        real(dp), dimension(l_min_min:, :), intent(in)     :: cl_old, cl_new
        type(genvec),                       intent(inout)  :: s
        integer(i4b),                       intent(out), optional :: ierr

        integer(i4b) :: i, l, m, lmax, err
        real(dp), allocatable, dimension(:,:,:)   :: inv_sqrt_S_old, sqrt_S_new

        if (present(ierr)) ierr = 0

        allocate(inv_sqrt_S_old(nmaps, nmaps, 0:l_max_max))
        allocate(sqrt_S_new(nmaps, nmaps, 0:l_max_max))

        if (any(cl_old /= cl_new)) then
           call cls2sqrt_covar(2, cl_old, inv_sqrt_S_old, err)
           if (err /= 0) then
              deallocate(inv_sqrt_S_old, sqrt_S_new)
              ierr = err
              return
           end if
           call cls2sqrt_covar(1, cl_new, sqrt_S_new, err)
           do l = l_min_min, l_max_max
              do m = -l, l
                 i = lm2ind(l, m)
                 s%cmb_amp(i,:) = matmul(sqrt_S_new(:,:,l), matmul(inv_sqrt_S_old(:,:,l), s%cmb_amp(i,:)))
              end do
           end do
        end if
        
        deallocate(inv_sqrt_S_old)
        deallocate(sqrt_S_new)

    end subroutine scale_s

    ! Mode = 1   =>   compute sqrt(S)
    ! Mode = 2   =>   compute inv(sqrt(S))
    subroutine cls2sqrt_covar(mode, cls, covar, ierr)
      implicit none
      
      integer(i4b),                              intent(in)  :: mode
      real(dp),     dimension(l_min_min:,1:),    intent(in)  :: cls
      real(dp),     dimension(1:,1:,0:),         intent(out) :: covar
      integer(i4b),                              intent(out), optional :: ierr
      
      integer(i4b) :: i, j, k, l, info
      logical(lgt),              dimension(3)     :: acceptable_mode
      real(dp),     allocatable, dimension(:)     :: w
      real(dp),     allocatable, dimension(:,:)   :: A
      real(dp),                  dimension(1:300) :: work

      if (present(ierr)) ierr = 0
      
      allocate(w(nmaps), A(nmaps, nmaps))
      
      covar(:,:,0:1) = 0.d0
      do l = l_min_min, l_max_max

         if (nmaps == 1) then
            if (mode == 1) then
               covar(1,1,l) = sqrt(cls(l,1))
            else if (mode == 2) then
               covar(1,1,l) = 1.d0 / sqrt(cls(l,1))             
            end if
         else
         
            A      = 0.d0
            A(1,1) = cls(l,1)
            A(1,2) = cls(l,2)
            A(2,1) = cls(l,2)
            A(2,2) = cls(l,3)
            if (size(cls(l,:)) == 4) A(3,3) = cls(l,4)
            
            acceptable_mode = .true.
            do i = 1, nmaps
               if (A(i,i) <= 0.d0) then
                  acceptable_mode(i) = .false.
                  A(i,:) = 0.d0
                  A(:,i) = 0.d0
                  A(i,i) = 1.d0
               end if
            end do
            
            call dsyev('V', 'L', nmaps, A, nmaps, w, work, size(work), info)
            if (minval(w) >= 0.d0) then
               if (mode == 1) then
                  w = sqrt(w)
               else if (mode == 2) then
                  w = 1.d0 / sqrt(w)
               end if
            else
               if (present(ierr)) then
                  !write(*,*) l, w
                  ierr = 1
                  deallocate(w, A)
                  return
               else
                  write(*,*) 'Error: Cl has negative eigenvalues.'
                  write(*,*) '       Multipole l = ', l
                  stop
               end if
            end if
            covar(:,:,l) = 0.d0
            do i = 1, nmaps
               do j = 1, nmaps
                  covar(i,j,l) = covar(i,j,l) + w(1) * A(i,1) * A(j,1)
                  covar(i,j,l) = covar(i,j,l) + w(2) * A(i,2) * A(j,2)
                  covar(i,j,l) = covar(i,j,l) + w(3) * A(i,3) * A(j,3)
               end do
            end do
            
            do i = 1, nmaps
               if (.not. acceptable_mode(i)) then
                  covar(i,:,l) = 0.d0
                  covar(:,i,l) = 0.d0
               end if
            end do
            
         end if
      end do
      
      deallocate(w, A)

    end subroutine cls2sqrt_covar




    subroutine check_spline(spline, cls, ok, rel_diff)
        implicit none
        real(dp), dimension(l_min_min:, :), intent(in)          :: spline, cls
        real(dp), dimension(l_min_min:, :), intent(out)         :: rel_diff
        real(dp)                                                :: norm
        integer(i4b)                                            :: i, l, j
        logical(lgt), dimension(num_modeel)                     :: ok
        character(len=128)                                      :: filename

        ok = .true.

        do i = 1, num_splinemodes
            norm = sqrt(sum(cls(l_min(i):l_max(i), i)*cls(l_min(i):l_max(i),i)))
            norm = norm/(sqrt(real(l_max(i), dp)-real(l_min(i), dp)))
            rel_diff(:, i) = (spline(:, i) - cls(:, i))/norm
        end do

        !rel_diff = (spline-cls)/cls
        !rel_diff = (spline-cls)

        do i = 1, num_splinemodes
            do l = l_min(i), l_max(i)
                if (abs(rel_diff(l, i)) > tol) then
                    ok(i) = .false.
                    exit
                end if
            end do
        end do

    end subroutine check_spline

    subroutine get_bb_inamps(bb_cl, bb_amps)
        implicit none
        real(dp), dimension(l_min(4):l_max(4)), intent(in)      :: bb_cl
        real(dp), dimension(num_bbtempl),       intent(out)     :: bb_amps
        real(dp), dimension(num_bbtempl, num_bbtempl)           :: templmat
        real(dp), dimension(num_bbtempl)                        :: b
        integer(i4b)                                            :: i, j, l

        templmat = 0.d0
        b = 0.d0

        do i = 1, num_bbtempl
            do j = 1, num_bbtempl
                templmat(i, j) = sum(bb_templ(l_min(4):l_max(4), i)*bb_templ(l_min(4):l_max(4), j))
            end do
            b(i) = sum(bb_templ(l_min(4):l_max(4), i)*bb_cl(l_min(4):l_max(4)))
        end do

        call solve_system_real(templmat, bb_amps, b)

    end subroutine get_bb_inamps
        

    subroutine sample_splined_spectrum(s, rng_handle, cl)
        implicit none

        type(genvec),                  intent(inout)    :: s
        real(dp), dimension(0:, :, :), intent(inout)    :: cl
        real(dp), allocatable, dimension(:, :)          :: currcl, newcl
        real(dp), allocatable, dimension(:, :)          :: rel_diff
        type(spline_amps), dimension(num_splinemodes)   :: amps
        type(spline_amps), dimension(num_splinemodes)   :: cl2
        real(dp),  dimension(num_bbtempl)               :: bb_amps
        real(dp)                                        :: currloglike
        real(dp)                                        :: newloglike
        real(dp)                                        :: currval
        real(dp)                                        :: propval
        real(dp)                                        :: lp1, lpn
        integer(i4b), allocatable, dimension(:, :)      :: acc
        integer(i4b)                                    :: maxnodenum

        integer(i4b)                                    :: i, j, k, m, l
        logical(lgt), dimension(num_modeel)              :: ok

        character(len=128)                              :: filename
        type(planck_rng)                                :: rng_handle

        if (.not. allocated(l_max)) then
            stop 'Module not initialized. Call initialize_cl_spline_mod prior to using sample_splined_spectrum.'
        end if

        do i = 1, num_splinemodes
            allocate(amps(i)%amps(cl_spline_data(i)%num_nodes))
            allocate(cl2(i)%amps(cl_spline_data(i)%num_nodes))
        end do
        allocate(currcl(l_min_min:l_max_max, num_modeel))
        allocate(newcl(l_min_min:l_max_max, num_modeel))
        allocate(rel_diff(l_min_min:l_max_max, num_modeel))

        currcl = 0.d0
        newcl = 0.d0
        rel_diff = 0.d0
        bb_amps = 0.d0

        !currcl(l_min(1):l_max(1), 1) = cl(l_min(1):l_max(1), 1, 1)     
        do l = l_min(1), l_max(1)
            currcl(l, 1) = cl(l, 1, 1)*real(l*(l+1), dp)/(2.d0*pi)
        end do
        if (num_modes > 1) then
            !currcl(l_min(2):l_max(2), 2) = cl(l_min(2):l_max(2), 2, 1)     
            !currcl(l_min(3):l_max(3), 3) = cl(l_min(3):l_max(3), 2, 2)     
            do l = l_min(2), l_max(2)
                currcl(l, 2) = cl(l, 2, 1)*real(l*(l+1), dp)/(2.d0*pi)
            end do
            do l = l_min(3), l_max(3)
                currcl(l, 3) = cl(l, 2, 2)*real(l*(l+1), dp)/(2.d0*pi)
            end do
            if (num_modes > 2) then 
                do l = l_min(4), l_max(4)
                    currcl(l, 4) = cl(l, 3, 3)*real(l*(l+1), dp)/(2.d0*pi)
                end do
                !currcl(l_min(4):l_max(4), 4) = cl(l_min(4):l_max(4), 3, 3)
            end if
        end if

        !filename = 'incls.dat'
        !call print_cls(currcl, filename)
        do i = 1, num_modeel
            if (i <= 3) then
                do j = 1, cl_spline_data(i)%num_nodes
                    amps(i)%amps(j) = currcl(int(cl_spline_data(i)%nodes(j)), i)
                end do
            else
                call get_bb_inamps(currcl(l_min(4):l_max(4), 4), bb_amps)
            end if
        end do

        if (autonodes .and. firstcall) then
            call adjust_nodes(currcl(l_min_min:, 1:num_splinemodes), amps)

            do i = 1, num_splinemodes
                if (size(cl2(i)%amps) /= cl_spline_data(i)%num_nodes) then
                    deallocate(cl2(i)%amps)
                    allocate(cl2(i)%amps(cl_spline_data(i)%num_nodes))
                    deallocate(cl_spline_data(i)%propose_sd)
                    allocate(cl_spline_data(i)%propose_sd(cl_spline_data(i)%num_nodes))
                    cl_spline_data(i)%propose_sd = 1.d0
                    deallocate(cl_spline_data(i)%num_iter)
                    allocate(cl_spline_data(i)%num_iter(cl_spline_data(i)%num_nodes))
                end if
                !print *, cl_spline_data(i)%nodes
            end do
        end if

        lp1 = 1.d30
        lpn = 1.d30

        do i = 1, num_modeel
            if (i <= 3) then
                call spline(cl_spline_data(i)%nodes, &
                   & amps(i)%amps, lp1, lpn, cl2(i)%amps)
                call get_splined_cls(i, amps(i)%amps, &
                   & cl2(i)%amps, currcl(l_min(i):l_max(i), i))
            else
                call get_bbcl(bb_amps, currcl(l_min(i):l_max(i), i))
            end if
        end do

        currloglike = eval_loglike(s, currcl, currcl)

        if (currloglike == 1) then
            filename = 'zerocls.dat'
            call print_cls(currcl, filename)
            do l = l_min_min, l_max_max
               if (.not. is_pos_def(l, currcl(l, :))) then
                  write(*,*) 'posdef', l, currcl(l,:)
               end if
            end do
            stop 'Initial likelihood is zero'
        end if

        if (firstcall) then
            !print *, 'Adjusting proposal density'
            call adjust_proposal_sd(amps, s, bb_amps)
            !print *, 'Adjusting iteration number'
            call adjust_iternum(amps, s, bb_amps)
            !DEBUG
            !filename = 'incls.dat'
            !call print_cls(currcl, filename)
        end if

        maxnodenum = 0
        do i = 1, num_modeel
            if (i <= 3) then
                if (cl_spline_data(i)%num_nodes > maxnodenum) then
                    maxnodenum = cl_spline_data(i)%num_nodes
                end if
            else
                if (num_bbtempl > maxnodenum) then
                    maxnodenum = num_bbtempl
                end if
            end if
        end do 

        allocate(acc(maxnodenum, num_modeel))
        if (firstcall) then
            firstcall = .false.
            allocate(no_acc_strikes(maxnodenum, num_modeel))
            no_acc_strikes = 0
        end if
        acc = 0

        do m = 1, num_iter_main
            do i = 1, num_modeel
                if (i <= 3) then
                    do j = 1, cl_spline_data(i)%num_nodes
                        newcl = currcl
                        currval = amps(i)%amps(j)
                        do k = 1, cl_spline_data(i)%num_iter(j)
                            call draw_nodeamp(amps(i)%amps(j), &
                                & propval, j, i, rng_handle)
                            amps(i)%amps(j) = propval
                            call spline(cl_spline_data(i)%nodes, &
                                & amps(i)%amps, lp1, lpn, cl2(i)%amps)
                            call get_splined_cls(i, amps(i)%amps, &
                                & cl2(i)%amps, newcl(l_min(i):l_max(i), i))
                            newloglike = eval_loglike(s, currcl, newcl)
                            if  (exp(newloglike-currloglike) > &
                                & rand_uni(rng_handle) .and. newloglike &
                                & .ne. -1.d30) then
                                currloglike = newloglike
                                currval = propval
                                call scale_s(s, currcl, newcl)
                                currcl = newcl
                                acc(j, i) = acc(j, i) + 1
                            end if
                            amps(i)%amps(j) = currval
                        end do
                    end do
                else
                    do j = 1, num_bbtempl
                        newcl = currcl
                        currval = bb_amps(j)
                        do k = 1, bb_num_iter(j)
                            call draw_nodeamp(bb_amps(j), &
                                & propval, j, i, rng_handle)
                            bb_amps(j) = propval
                            call get_bbcl(bb_amps, &
                                & newcl(l_min(i):l_max(i), i))
                            newloglike = eval_loglike(s, currcl, newcl)
                            if  (exp(newloglike-currloglike) > &
                                & rand_uni(rng_handle) .and. newloglike &
                                & .ne. -1.d30) then
                                currloglike = newloglike
                                currval = propval
                                call scale_s(s, currcl, newcl)
                                currcl = newcl
                                acc(j, i) = acc(j, i) + 1
                            end if
                            bb_amps(j) = currval
                        end do
                    end do
                end if
            end do
        end do

        do i = 1, num_modeel
            if (i <= 3) then
                do j = 1, cl_spline_data(i)%num_nodes
                    if (acc(j, i) == 0) then
                        no_acc_strikes(j, i) = no_acc_strikes(j, i) + 1
                        if (no_acc_strikes(j, i) == 300) then
                            print *, 'No samples accepted for mode ', i, ', node ', j, ' for the third time in a row. Stopping.' 
                            stop 
                        end if
                    else
                        if (no_acc_strikes(j, i) /= 0) then
                            no_acc_strikes(j, i) = 0
                        end if
                    end if
                end do
            end if
        end do

        !cl(l_min(1):l_max(1), 1, 1) = currcl(l_min(1):l_max(1), 1)
        !if (num_modes > 1) then
        !    cl(l_min(2):l_max(2), 1, 2) = currcl(l_min(2):l_max(2), 2)
        !    cl(l_min(2):l_max(2), 2, 1) = currcl(l_min(2):l_max(2), 2)
        !    cl(l_min(3):l_max(3), 2, 2) = currcl(l_min(3):l_max(3), 3)
        !    if (num_modes > 2) then
        !        cl(l_min(4):l_max(4), 3, 3) = currcl(l_min(4):l_max(4), 4)
        !    end if
        !end if
        do l = l_min(1), l_max(1)
            cl(l, 1, 1) = currcl(l, 1)/real(l*(l+1), dp)*2.d0*pi
        end do
        if (num_modes > 1) then
            do l = l_min(2), l_max(2)
                cl(l, 1, 2) = currcl(l, 2)/real(l*(l+1), dp)*2.d0*pi
                cl(l, 2, 1) = currcl(l, 2)/real(l*(l+1), dp)*2.d0*pi
                cl(l, 2, 2) = currcl(l, 3)/real(l*(l+1), dp)*2.d0*pi
            end do
            if (num_modes > 2) then
                do l = l_min(4), l_max(4)
                    cl(l, 3, 3) = currcl(l, 4)/real(l*(l+1), dp)*2.d0*pi
                end do
            end if
        end if

!        print *, 'Sample complete'
        do i = 1, num_splinemodes
            deallocate(amps(i)%amps)
            deallocate(cl2(i)%amps)
        end do
        deallocate(currcl)
        deallocate(newcl)

    end subroutine sample_splined_spectrum

    subroutine cleanup_cl_spline_mod()
        implicit none
        integer(i4b)            :: i

        do i = 1, num_splinemodes
            deallocate(cl_spline_data(i)%nodes)
            deallocate(cl_spline_data(i)%num_iter)
            deallocate(cl_spline_data(i)%propose_sd)
        end do
        
        if (bb) then
            deallocate(bb_templ)
            deallocate(bb_propose_sd)
            deallocate(bb_num_iter)
        end if

    end subroutine cleanup_cl_spline_mod

end module comm_cl_spline_mod
