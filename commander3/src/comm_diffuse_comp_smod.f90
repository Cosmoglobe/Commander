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
submodule (comm_diffuse_comp_mod) comm_diffuse_comp_smod
  
contains

  module subroutine initDiffuse(self, cpar, id, id_abs)
    !
    ! Routine that initializes a diffuse type component. 
    !
    ! Arguments:
    ! self: comm_diffuse_comp 
    !       Diffuse type component
    !
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! id: integer
    !       Integer ID of the diffuse component with respect to the active components
    !
    ! id_abs: integer
    !       Integer ID of the diffuse component with respect to all components defined in the parameter file
    !       (and also the id in the 'cpar' parameter)
    !
    ! Returns:
    !       The diffuse component parameter is returned (self).
    !       Any other changes are done internally
    !
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    
  end subroutine initDiffuse

  module subroutine initLmaxSpecind(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    integer(i4b) :: i, j, k, l, m, nmaps, p, ierr



  end subroutine initLmaxSpecind

  module subroutine initPixregSampling(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    integer(i4b) :: i, j, k, l, m, n, p, pr, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    character(len=16),         dimension(1000) :: pixreg_label
    character(len=16),         dimension(1000) :: pixreg_prior
    character(len=512),        dimension(1000) :: tokens
    integer(i4b), allocatable, dimension(:)    :: sum_pix
    real(dp),     allocatable, dimension(:)    :: sum_theta, sum_proplen, sum_nprop
    real(dp),     allocatable, dimension(:,:)  :: m_in, m_out, buffer
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    logical(lgt) :: all_fixed
    class(comm_mapinfo), pointer :: info2 => null(), info_ud => null(), info3 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 


  end subroutine initPixregSampling

  module subroutine initSpecindProp(self,cpar, id, id_abs)
    !
    !  Subroutine that initializes the spectral index proposal matrix
    !
    !  Arguments:
    !  ----------
    !  self: comm_diffuse_component
    !       Object that contains all of the diffuse component properties
    !  cpar: comm_params
    !       Object that contains parameters input to Commander from parameter
    !       file.
    !  id: int
    !       index of the components, out of only the ones being used
    !  id_abs: int
    !       index of the component, out of all, whether they are used or not. 
    !  Returns:
    !  --------
    !  None, but modifies self by changing the proposal matrices.
    !
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    
  end subroutine initSpecindProp


  module subroutine initDiffPrecond(comm, samp_group)
    implicit none
    integer(i4b),                intent(in) :: comm, samp_group


  end subroutine initDiffPrecond

  module subroutine initDiffPrecond_diagonal(comm, samp_group)
    implicit none
    integer(i4b),                intent(in) :: comm
    integer(i4b),                intent(in) :: samp_group


  end subroutine initDiffPrecond_diagonal


  module subroutine initDiffPrecond_pseudoinv(comm, samp_group)
    implicit none
    integer(i4b),                intent(in) :: comm, samp_group

    integer(i4b) :: i, i1, i2, j, k1, k2, q, l, m, n
    real(dp)     :: t1, t2
    integer(i4b), allocatable, dimension(:) :: ind
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: p1 => null(), p2 => null()
    real(dp),     allocatable, dimension(:,:) :: mat

    if (npre == 0) return
    if (allocated(P_cr(samp_group)%invM_diff)) return
    
    if (.not. allocated(diffComps)) then
       ! Set up an array of all the diffuse components
       allocate(diffComps(npre))
       c => compList
       i =  1
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             diffComps(i)%p => c
             i              =  i+1
          end select
          c => c%nextComp()
       end do
       info_pre => comm_mapinfo(comm, nside_pre, lmax_pre, nmaps_pre, nmaps_pre==3)
    end if
    
    ! Allocate space for pseudo-inverse of U
    call wall_time(t1)
    allocate(P_cr(samp_group)%invM_diff(0:lmax_pre,info_pre%nmaps))
    do j = 1, info_pre%nmaps
       do l = 0, lmax_pre
          allocate(P_cr(samp_group)%invM_diff(l,j)%M(npre,numband+npre))
       end do
    end do

  end subroutine initDiffPrecond_pseudoinv


  module subroutine updateDiffPrecond(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update


  end subroutine updateDiffPrecond

  module subroutine updateDiffPrecond_diagonal(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    integer(i4b) :: i, j, k, k1, k2, l, m, n, ierr, unit, p, q
    real(dp)     :: W_ref, t1, t2
    logical(lgt), save :: first_call = .true.
    integer(i4b), allocatable, dimension(:) :: ind
    real(dp),     allocatable, dimension(:) :: W
    real(dp),     allocatable, dimension(:) :: cond, W_min
    real(dp),     allocatable, dimension(:,:) :: alm

       
  end subroutine updateDiffPrecond_diagonal


  module subroutine updateDiffPrecond_pseudoinv(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    integer(i4b) :: i, j, k, l, n, ierr, p, q
    real(dp)     :: t1, t2, Cl
    real(dp),     allocatable, dimension(:,:) :: mat


  end subroutine updateDiffPrecond_pseudoinv
    
  
  ! Evaluate amplitude map in brightness temperature at reference frequency
  module subroutine updateDiffuseMixmat(self, theta, beta, band, df, par)
    implicit none
    class(comm_diffuse_comp),                  intent(inout)           :: self
    class(comm_map),           dimension(:),   intent(in),    optional :: theta
    real(dp),  dimension(:,:,:),               intent(in),    optional :: beta  ! Not used here
    integer(i4b),                              intent(in),    optional :: band
    class(map_ptr), dimension(:),              intent(inout), optional :: df    ! Derivative of mixmat with respect to parameter par; for Jeffreys prior
    integer(i4b),                              intent(in),    optional :: par   ! Parameter ID for derivative

    integer(i4b) :: i, j, k, l, n, p, p_min, p_max, nmaps, ierr
    real(dp)     :: lat, lon, t1, t2
    logical(lgt) :: precomp, mixmatnull, bad ! NEW
    character(len=2) :: ctext
    real(dp),        allocatable, dimension(:,:,:) :: theta_p
    real(dp),        allocatable, dimension(:)     :: nu, s, buffer, buff2
    class(comm_mapinfo),          pointer          :: info => null(), info_tp => null()
    class(comm_map),              pointer          :: t => null(), tp => null(), td => null()
    class(map_ptr),  allocatable, dimension(:)     :: theta_prev


    if (trim(self%type) == 'md') return

    
    ! Copy over alms from input structure, and compute pixel-space parameter maps
    if (present(theta)) then
       do i = 1, self%npar
          self%theta(i)%p%alm = theta(i)%alm          
       end do
    end if
    
    ! Compute mixing matrix
    do i = 1, numband
       
       ! Only update requested band if present
       if (present(band)) then
          if (i /= band) cycle
       end if

       ! Compute spectral parameters at the correct resolution for this channel
       if (self%npar > 0) then
          nmaps = min(data(i)%info%nmaps, self%theta(1)%p%info%nmaps)
          allocate(theta_p(0:data(i)%info%np-1,nmaps,self%npar))
          !if (self%myid==0) write(*,*) 'ee1 npix nmaps npar', data(i)%info%np-1,nmaps,self%npar
          
          do j = 1,self%npar
             info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, &
                  & self%theta(j)%p%info%lmax, nmaps, data(i)%info%pol)
             td => comm_map(info)
             
             ! if any polarization is alm sampled. Only use alms to set polarizations with alm sampling
             if (any(self%lmax_ind_pol(1:self%poltype(j),j) >= 0)) then
                t => comm_map(info)
                t%alm(:,1:nmaps) = self%theta(j)%p%alm(:,1:nmaps)
                call t%Y_scalar()
                do p = 1,self%poltype(j)
                   if (self%lmax_ind_pol(p,j) < 0) cycle
                   if (self%poltype(j) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(j)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(j)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(j) 
                      stop
                   end if

                   do k = p_min,p_max
                      td%map(:,k) = t%map(:,k)
                   end do
                end do
                call t%dealloc(); deallocate(t)
             end if

             ! if any polarization is local sampled. Only set theta using polarizations with local sampling
             if (any(self%lmax_ind_pol(1:self%poltype(j),j) < 0)) then
                info_tp => comm_mapinfo(self%theta(j)%p%info%comm, self%theta(j)%p%info%nside, &
                     & 2*self%theta(j)%p%info%nside, nmaps, .false.)
                tp    => comm_map(info_tp)
                tp%map(:,1:nmaps) = self%theta(j)%p%map(:,1:nmaps)
                call tp%YtW_scalar()
                info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, &
                     & data(i)%info%lmax, nmaps, .false.)
                t    => comm_map(info)
                call tp%alm_equal(t)
                call t%Y_scalar()
                where (t%map < self%p_uni(1,j))
                   t%map = self%p_uni(1,j)
                end where
                where (t%map > self%p_uni(2,j))
                   t%map = self%p_uni(2,j)
                end where
                call tp%dealloc(); deallocate(tp)

                
                call wall_time(t2)
                do p = 1,self%poltype(j)
                   if (self%lmax_ind_pol(p,j) >= 0) cycle
                   if (self%poltype(j) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(j)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(j)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(j) 
                      stop
                   end if

                   do k = p_min,p_max
                      td%map(:,k) = t%map(:,k)
                   end do
                end do

                call t%dealloc(); deallocate(t)

                !if (info%myid == 0) write(*,*) 'udgrade = ', t2-t1
             end if
             theta_p(:,:,j) = td%map
             call td%dealloc(); deallocate(td)
          end do
       end if

       do l = 0, data(i)%ndet
          
          ! Don't update null mixing matrices
          if (self%F_null(i,l)) then
             if (present(df)) df(i)%p%map = 0.d0
             cycle
          end if
          
          
          ! Loop over all pixels, computing mixing matrix for each
          !allocate(theta_p(self%npar,self%nmaps))
          call wall_time(t1)
          do j = 0, self%F(i,l)%p%info%np-1
             if (self%latmask >= 0.d0) then
                call pix2ang_ring(data(i)%info%nside, data(i)%info%pix(j+1), lat, lon)
                lat = 0.5d0*pi-lat
                if (abs(lat) < self%latmask) then
                   self%F(i,l)%p%map(j,:) = 0.d0
                   cycle
                end if
             end if

!          if (all(self%lmax_ind_mix(1:min(self%nmaps,data(i)%info%nmaps)) == 0)) then  !if (self%lmax_ind == 0) then
!             cycle
!          end if

             ! NEW ! Check band sensitivity before mixing matrix update
             ! Possible labels are "broadband", "cmb", "synch", "dust", "co10", "co21", "co32", "ff", "ame"
             if (data(i)%comp_sens == "broadband") then
                ! If broadband, calculate mixing matrix
                mixmatnull = .false.
             else
                ! If component sensitivity, only calculate mixmat on that component.
                mixmatnull = .true.
                If (data(i)%comp_sens == self%label) then
                   mixmatnull = .false.
                end if
             end if
             
             ! Temperature
             if (.not. only_pol) then
                if (self%npar > 0) then
                   if (mixmatnull) then
                      self%F(i,l)%p%map(j,1) = 0.0
                   else
                      if (trim(self%label) == 'dust' .and. any(theta_p(j,1,:)==0.d0)) then
                         write(*,*) 'dust beta and T can not be null, crashing'
                         write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                         stop !debug, replace by proper stop and error message
                      end if
                      self%F(i,l)%p%map(j,1) = self%F_int(1,i,l)%p%eval(theta_p(j,1,:)) * data(i)%gain * self%cg_scale(1)
                      !write(*,*) i, j, theta_p(j,1,:), self%F_int(1,i,l)%p%eval(theta_p(j,1,:)), self%F(i,l)%p%map(j,1)
                   end if
                else
                   if (mixmatnull) then 
                      self%F(i,l)%p%map(j,1) = 0.0
                   else
                      self%F(i,l)%p%map(j,1) = self%F_int(1,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale(1)
                   end if
                end if
             end if
             
             ! Polarization
             if (self%nmaps == 3 .and. data(i)%info%nmaps == 3) then
                ! Stokes Q
                if (self%npar == 0) then
                   self%F(i,l)%p%map(j,2) = self%F(i,l)%p%map(j,1) 
                else if (all(self%poltype < 2)) then
                   self%F(i,l)%p%map(j,2) = self%F(i,l)%p%map(j,1) 
                else
                   if (trim(self%label) == 'dust' .and. any(theta_p(j,2,:)==0.d0)) then
                      write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                   end if
                   if (self%npar > 0) then
                      self%F(i,l)%p%map(j,2) = self%F_int(2,i,l)%p%eval(theta_p(j,2,:)) * data(i)%gain * self%cg_scale(2)
                   else
                      self%F(i,l)%p%map(j,2) = self%F_int(2,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale(2)
                   end if
                end if
                
                ! Stokes U
                if (self%npar == 0) then
                   self%F(i,l)%p%map(j,3) = self%F(i,l)%p%map(j,2) 
                else if (all(self%poltype < 3)) then
                   self%F(i,l)%p%map(j,3) = self%F(i,l)%p%map(j,2) 
                else
                   if (trim(self%label) == 'dust' .and. any(theta_p(j,1,:)==0.d0)) then
                      write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                   end if
                   if (self%npar > 0) then
                      self%F(i,l)%p%map(j,3) = self%F_int(3,i,l)%p%eval(theta_p(j,3,:)) * data(i)%gain * self%cg_scale(3)
                   else
                      self%F(i,l)%p%map(j,3) = self%F_int(3,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale(3)
                   end if
                end if
             end if
          
             ! Compute derivative if requested
             if (present(df) .and. l == 0) then
                if (self%npar > 0) then
                   do k = 1, nmaps
                      if (k <= self%poltype(par)) then
                         df(i)%p%map(j,k) = self%F_int(k,i,l)%p%eval_deriv(theta_p(j,k,:),par) * data(i)%gain * self%cg_scale(k)
                      end if
                   end do
                else
                   df(i)%p%map = 0.d0
                end if
             end if

          end do

          call wall_time(t2)
          !if (self%x%info%myid == 0) write(*,*) 'eval = ', t2-t1
                
          ! Compute mixing matrix average; for preconditioning only
          allocate(buffer(self%nmaps), buff2(self%nmaps))
          do j = 1, min(self%nmaps, data(i)%info%nmaps)
             self%F_mean(i,l,j) = sum(self%F(i,l)%p%map(:,j))
          end do
          buff2 = self%F_mean(i,l,:)
          !call mpi_barrier(mpi_comm_world, ierr)
          call mpi_allreduce(buff2, buffer, self%nmaps, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
          self%F_mean(i,l,:) = buffer / self%F(i,l)%p%info%npix
          deallocate(buffer,buff2)
    
       end do
       if (allocated(theta_p)) deallocate(theta_p)
    end do


    ! Request preconditioner update
    recompute_diffuse_precond = .true.

  end subroutine updateDiffuseMixmat



  module function evalDiffuseBand(self, band, amp_in, pix, alm_out, det) result(res)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    integer(i4b),                                 intent(in),  optional :: det
    real(dp),        dimension(:,:), allocatable                        :: res

    res = 0d0


  end function evalDiffuseBand

  ! Return component projected from map
  module function projectDiffuseBand(self, band, map, alm_in, det) result(res)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    integer(i4b),                                 intent(in), optional  :: det
    real(dp),        dimension(:,:), allocatable                        :: res

    integer(i4b) :: i, nmaps, d
    logical(lgt) :: alm_in_
    class(comm_mapinfo), pointer :: info_in => null(), info_out => null()
    class(comm_map),     pointer :: m       => null(), m_out    => null()

    if (self%F_null(band,0)) then
       if (.not. allocated(res)) allocate(res(0:self%x%info%nalm-1,self%x%info%nmaps))
       res = 0.d0
       return
    end if

    res = 0.d0


  end function projectDiffuseBand


  module subroutine applyDiffPrecond(x, samp_group)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x
    integer(i4b),                     intent(in)    :: samp_group


  end subroutine applyDiffPrecond


  module subroutine applyDiffPrecond_diagonal(x, samp_group)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x
    integer(i4b),                     intent(in)    :: samp_group

    integer(i4b)              :: i, j, k, l, m, nmaps
    real(dp), allocatable, dimension(:,:)   :: alm
    real(dp), allocatable, dimension(:,:,:) :: y

    if (npre == 0) return
    
    ! Reformat linear array into y(npre,nalm,nmaps) structure
    allocate(y(npre,0:info_pre%nalm-1,info_pre%nmaps))
    y = 0.d0
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       call cr_extract_comp(diffComps(i)%p%id, x, alm)
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          y(i,k,1:nmaps) = alm(j,1:nmaps)
       end do
       deallocate(alm)
    end do

    ! Multiply with preconditioner
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
          if (P_cr(samp_group)%invM_diff(i,j)%n == 0) cycle
          y(P_cr(samp_group)%invM_diff(i,j)%ind,i,j) = &
               & matmul(P_cr(samp_group)%invM_diff(i,j)%M, y(P_cr(samp_group)%invM_diff(i,j)%ind,i,j))
       end do
    end do

    ! Reformat y(npre,nalm,nmaps) structure into linear array
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       allocate(alm(0:diffComps(i)%p%x%info%nalm-1,nmaps))
       alm = 0.d0
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          alm(j,1:nmaps) = y(i,k,1:nmaps)
       end do
       call cr_insert_comp(diffComps(i)%p%id, .false., alm, x)
       deallocate(alm)
    end do
    
    deallocate(y)

  end subroutine applyDiffPrecond_diagonal


  module subroutine applyDiffPrecond_pseudoinv(x, samp_group)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x
    integer(i4b),                     intent(in)    :: samp_group

    integer(i4b)              :: i, ii, j, k, l, m, p, q, qq, nmaps, npre_int
    real(dp)                  :: t1, t2
    real(dp), allocatable, dimension(:)     :: w, w2
    real(dp), allocatable, dimension(:,:)   :: alm
    real(dp), allocatable, dimension(:,:,:) :: y, z
    class(comm_map), pointer                :: invN_x => null()

    if (.not. allocated(ind_pre)) return

    npre_int = size(ind_pre)
    if (npre_int <= 0) return
    
    ! Reformat linear array into y(npre,nalm,nmaps) structure
    allocate(y(npre_int,0:info_pre%nalm-1,info_pre%nmaps))
    allocate(z(npre_int,0:info_pre%nalm-1,info_pre%nmaps))
    y = 0.d0
    do i = 1, npre_int
       ii = ind_pre(i)
       nmaps = diffComps(ii)%p%x%info%nmaps
       call cr_extract_comp(diffComps(ii)%p%id, x, alm)
       do j = 0, diffComps(ii)%p%x%info%nalm-1
          call diffComps(ii)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          y(i,k,1:nmaps) = alm(j,1:nmaps)
       end do
       deallocate(alm)
    end do

    ! Frequency-dependent terms
    z = 0.d0
    do k = 1, numband
       invN_x => comm_map(data(k)%info)
       nmaps  =  data(k)%info%nmaps
       
       ! Sum over (U^plus)^t
       !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,m,j,q,qq,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, data(k)%info%nalm-1
          call data(k)%info%i2lm(i, l, m)
          if (l > info_pre%lmax) cycle
          call info_pre%lm2i(l,m,j)
          do q = 1, npre_int
             qq = ind_pre(q)
             do p = 1, nmaps
                invN_x%alm(i,p) = invN_x%alm(i,p) + P_cr(samp_group)%invM_diff(l,p)%M(qq,k) * y(q,j,p)
             end do
          end do
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL

       ! Multiply by T
       call invN_x%WY
       call data(k)%N%N(invN_x)
       call invN_x%YtW
       do i = 1, nmaps
          invN_x%alm(:,i) = invN_x%alm(:,i) * data(k)%N%alpha_nu(i)**2
       end do

       ! Sum over U^plus
       !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,m,j,q,qq,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, data(k)%info%nalm-1
          call data(k)%info%i2lm(i, l, m)
          if (l > info_pre%lmax) cycle
          call info_pre%lm2i(l,m,j)
          do q = 1, npre_int
             qq = ind_pre(q)
             do p = 1, nmaps
                z(q,j,p) = z(q,j,p) + P_cr(samp_group)%invM_diff(l,p)%M(qq,k) * invN_x%alm(i,p)
             end do
          end do
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL

       call invN_x%dealloc(); deallocate(invN_x)
    end do

    ! Prior terms
    call wall_time(t1)
    !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,k,j,m,w,w2,p)
    allocate(w(npre_int), w2(npre_int))
    !!$OMP DO SCHEDULE(guided)
    do i = 0, info_pre%nalm-1
       do p = 1, info_pre%nmaps
          !call info_pre%i2lm(i, l, m)
          l = info_pre%lm(1,i)
          m = info_pre%lm(2,i)
          w        = y(:,i,p)
          w2       = 0.d0
          do j = 1, npre_int
             do k = 1, npre_int
                w2(j) = w2(j) + P_cr(samp_group)%invM_diff(l,p)%M(ind_pre(k),numband+ind_pre(j))*w(k)
             end do
          end do
          w       = 0.d0
          do j = 1, npre_int
             do k = 1, npre_int
                w(j) = w(j) + P_cr(samp_group)%invM_diff(l,p)%M(ind_pre(j),numband+ind_pre(k))*w2(k)
             end do
          end do
          z(:,i,p) = z(:,i,p) + w
       end do
    end do
    deallocate(w, w2)
    !!$OMP END DO
    !!$OMP END PARALLEL
    call wall_time(t2)
    !if (info_pre%myid == 0 .or. info_pre%myid == 25) write(*,*) info_pre%myid, ', nalm = ', info_pre%nalm, real(t2-t1,sp)

    ! Reformat z(npre,nalm,nmaps) structure into linear array
    do i = 1, npre_int
       ii = ind_pre(i)
       nmaps = diffComps(ii)%p%x%info%nmaps
       allocate(alm(0:diffComps(ii)%p%x%info%nalm-1,nmaps))
       alm = 0.d0
       do j = 0, diffComps(ii)%p%x%info%nalm-1
          call diffComps(ii)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          alm(j,1:nmaps) = z(i,k,1:nmaps)
       end do
       call cr_insert_comp(diffComps(ii)%p%id, .false., alm, x)
       deallocate(alm)
    end do
    
    deallocate(y, z)

  end subroutine applyDiffPrecond_pseudoinv


  
  ! Dump current sample to HEALPix FITS file
  module subroutine dumpDiffuseToFITS(self, iter, chainfile, output_hdf, postfix, dir)
    !
    ! Routine that writes a diffuce component to FITS (and HDF) files. 
    !
    ! Arguments:
    ! self: comm_diffuse_comp 
    !       Diffuse type component
    !
    ! iter: integer
    !       Sample number in the Gibb's chain.
    !
    ! chainfile: hdf_file
    !       HDF file to write the component to
    !
    ! output_hdf: logical
    !       Logical parameter to tell whether or not to write the component to the specified HDF file
    !
    ! postfix: string
    !       A string label to be added to the end of FITS-files.
    !       (default format: cXXXX_kYYYYYY; XXXX = chain number, YYYYYY = sample number)
    !
    ! dir: string
    !       Output directory to which output is written
    !
    ! Returns:
    !       The diffuse component parameter is returned (self).
    !       Any other changes are done internally
    !
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i, l, j, k, m, ierr, unit
    integer(i4b)       :: p, p_min, p_max, npr, npol
    real(dp)           :: vals(10)
    logical(lgt)       :: exist, first_call = .true.
    character(len=6)   :: itext
    character(len=512) :: filename, path
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map), pointer :: map => null(), tp => null()
    real(dp), allocatable, dimension(:,:) :: sigma_l
    real(dp),     allocatable, dimension(:,:) :: dp_pixreg
    integer(i4b), allocatable, dimension(:,:) :: int_pixreg


  end subroutine dumpDiffuseToFITS

  ! Dump current sample to HEALPix FITS file
  module subroutine initDiffuseHDF(self, cpar, hdffile, hdfpath)
    implicit none
    class(comm_diffuse_comp),  intent(inout) :: self
    type(comm_params),         intent(in)    :: cpar    
    type(hdf_file),            intent(in)    :: hdffile
    character(len=*),          intent(in)    :: hdfpath

    integer(i4b)       :: i, j, l, m, ierr
    integer(i4b)       :: p, p_min, p_max, npr, npol
    real(dp)           :: md(4)
    character(len=512) :: path
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map), pointer     :: tp => null()
    real(dp),     allocatable, dimension(:,:) :: dp_pixreg
    integer(i4b), allocatable, dimension(:,:) :: int_pixreg
    

  end subroutine initDiffuseHDF
  
  module subroutine add_to_npre(n, nside, lmax, nmaps)
    implicit none
    integer(i4b), intent(in) :: n, nside, lmax, nmaps
    npre      = npre + n
    nside_pre = min(nside_pre, nside)
    lmax_pre  = max(lmax_pre, lmax)
    nmaps_pre = max(nmaps_pre, nmaps)
  end subroutine add_to_npre

  ! Sample spectral parameters. This subroutine is obsolete, it is defined and used in comm_nonlin_mod instead
  module subroutine sampleDiffuseSpecInd(self, cpar, handle, id, iter)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter

    integer(i4b) :: i, j, k, l, n, p, q, pix, ierr, ind(1), counter, n_ok
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp), pointer :: c => null()
    real(dp),     allocatable, dimension(:)   :: lnL, P_tot, F, theta, a_curr
    real(dp),     allocatable, dimension(:,:) :: amp, buffer, alm_old
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, accept_rate
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len
    logical(lgt) :: first_sample
    class(comm_mapinfo),            pointer :: info => null()

    return 

  end subroutine sampleDiffuseSpecInd
  
  module subroutine print_precond_mat(samp_group)
    implicit none
    integer(i4b), intent(in) :: samp_group

    integer(i4b) :: l, m, i, j, p
    real(dp), allocatable, dimension(:)   :: W
    real(dp), allocatable, dimension(:,:) :: mat

    if (info_pre%myid /= 0) return
    p = 2

    open(58,file=trim(outdir)//'/precond_W.dat',recl=1024)
    if (trim(precond_type) == 'diagonal') then
       do l = 0, info_pre%lmax
          call info_pre%lm2i(l, 0, i)
          write(*,*) 
          write(*,*) l 
          do j = 1, size(P_cr(samp_group)%invM_diff(i,p)%M(j,:),1)
             write(*,*) real(P_cr(samp_group)%invM_diff(i,p)%M(j,:),sp)
          end do
          allocate(W(P_cr(samp_group)%invM_diff(i,p)%n))
          call get_eigenvalues(P_cr(samp_group)%invM_diff(i,p)%M, W)
          write(58,*) l, real(W,sp)
          deallocate(W)
       end do
    else
       allocate(mat(npre,npre))
       do l = 0, info_pre%lmax
          mat = matmul(P_cr(samp_group)%invM_diff(l,p)%M, transpose(P_cr(samp_group)%invM_diff(l,p)%M))
          write(*,*) 
          write(*,*) l 
          do j = 1, npre
             write(*,*) real(mat(j,:),sp)
          end do
          allocate(W(npre))
          call get_eigenvalues(mat, W)
          write(58,*) l, real(W,sp)
          deallocate(W)
       end do
       deallocate(mat)
    end if
    close(58)


  end subroutine print_precond_mat

  module subroutine updateDiffuseFInt(self, band)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
    integer(i4b),             intent(in),   optional :: band

    integer(i4b) :: i, j, k

    if (present(band)) then
       do i = 1, data(band)%info%nmaps
          do j = 0, data(band)%ndet
             call self%F_int(i,band,j)%p%update(pol=i)
          end do
       end do
    else
       do k = 1, numband
          do i = 1, data(k)%info%nmaps
             do j = 0, data(k)%ndet
                call self%F_int(i,k,j)%p%update(pol=i)
             end do
          end do
       end do
    end if

  end subroutine updateDiffuseFInt

  module subroutine updateLowlPrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, q, lp, mp, myid, nalm, ierr, nmaps
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()
    real(dp),        allocatable, dimension(:,:) :: invM, buffer



  end subroutine updateLowlPrecond

  module subroutine applyLowlPrecond(self, alm, alm0)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(inout)  :: alm
    real(dp),                 dimension(0:,1:), intent(in)     :: alm0

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, q, myid, nalm, ntot, ierr
    real(dp), allocatable, dimension(:) :: y, yloc, buffer


  end subroutine applyLowlPrecond


  module subroutine updateDeflatePrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self

    real(dp)                  :: t1, t2, norm
    integer(i4b)              :: i, j, k, l, m, p, q, lp, mp, myid, nalm, ierr, nmaps, nmax, ndef, nside
    logical(lgt)              :: add
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()
    real(dp),        allocatable, dimension(:,:) :: invM, buffer, Z
    real(dp),        allocatable, dimension(:) :: W


  end subroutine updateDeflatePrecond



  module subroutine applyDeflatePrecond(self, alm, Qalm)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(in)     :: alm
    real(dp),                 dimension(0:,1:), intent(out)    :: Qalm

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, myid, nalm, ntot, ierr, ndef
    real(dp), allocatable, dimension(:) :: y, ytot
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()

      Qalm = 0d0


  end subroutine applyDeflatePrecond

  module subroutine setup_needlets(info, nside_def, defmask, Z, ndef)
    implicit none
    class(comm_mapinfo), pointer,          intent(in)  :: info
    integer(i4b),                          intent(in)  :: nside_def
    class(comm_map),                       intent(in)  :: defmask
    real(dp),            dimension(0:,1:), intent(out) :: Z
    integer(i4b),                          intent(out) :: ndef


    z = 0d0
    ndef = 0

  end subroutine setup_needlets

  module subroutine applyMonoDipolePrior(self, handle)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
    type(planck_rng),         intent(inout)          :: handle

    integer(i4b) :: i, j, k, l, m, ierr, pix
    real(dp)     :: mu(0:3), a, b, Amat(0:3,0:3), bmat(0:3), v(0:3), corr_res(3)
    class(comm_map), pointer :: map, lr_map 
    class(comm_mapinfo), pointer :: info => null()
    real(dp), dimension(:), allocatable :: mask_list, corr_list, amp_list, intersect, all_thetas
    real(dp)     :: mean_intersect, std_intersect 
    real(dp)     :: diff_mono, diff_comp, mono_mix, comp_mix, prior_vals(2)
    class(comm_comp), pointer :: c => null()

    
  end subroutine applyMonoDipolePrior

  module subroutine nullify_monopole_amp(band)
    implicit none
    character(len=*), intent(in) :: band

    integer(i4b) :: i
    class(comm_comp), pointer :: c => null()


  end subroutine nullify_monopole_amp


  module function get_monopole_amp(band)
    implicit none
    character(len=*), intent(in) :: band
    real(dp)                     :: get_monopole_amp

    integer(i4b) :: l, m, ierr
    real(dp)     :: mono
    class(comm_comp), pointer :: c => null()

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%label) == trim(band)) then
             mono = -1.d100
             if (c%x%info%nalm > 0) then
                call c%x%info%i2lm(0,l,m)
                if (l == 0) then
                   mono = 1.d0/sqrt(4.d0*pi) * c%x%alm(0,1) * c%RJ2unit_(1)
                end if
             end if
             call mpi_allreduce(MPI_IN_PLACE, mono, 1, MPI_DOUBLE_PRECISION, MPI_MAX, c%x%info%comm, ierr)
             get_monopole_amp = mono
             return
          end if
       end select
       c => c%nextComp()
    end do

  end function get_monopole_amp

  module subroutine set_monopole_amp(band, mono)
    implicit none
    character(len=*), intent(in) :: band
    real(dp),         intent(in) :: mono

    integer(i4b) :: l, m
    class(comm_comp), pointer :: c => null()

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%label) == trim(band)) then
             if (c%x%info%nalm > 0) then
                call c%x%info%i2lm(0,l,m)
                if (l == 0) then
                   c%x%alm(0,1) = sqrt(4.d0*pi) * mono /c%RJ2unit_(1)
                end if
             end if
             return
          end if
       end select
       c => c%nextComp()
    end do

  end subroutine set_monopole_amp

end submodule comm_diffuse_comp_smod
