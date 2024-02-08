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
module comm_mh_specind_mod
  use comm_param_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_chisq_mod
  use comm_gain_mod
  use comm_line_comp_mod
  use comm_diffuse_comp_mod
  use comm_signal_mod
  use comm_utils
  use InvSamp_mod
  use powell_mod
  use comm_output_mod
  implicit none

contains


  subroutine sample_specind_mh_sample(outdir, cpar, handle, handle_noise)
    implicit none
    character(len=*),               intent(in)    :: outdir
    type(planck_rng),               intent(inout) :: handle, handle_noise
    type(comm_params) :: cpar

    real(dp)     :: chisq, my_chisq, chisq_old, chisq_new, chisq_prop
    integer(i4b) :: band, samp_group, ierr, i, j, k, pol, pix
    logical(lgt)  :: include_comp, reject, todo
    character(len=512) :: tokens(10), str_buff, operation
    class(comm_comp),   pointer           :: c => null()
    class(comm_map), pointer              :: invN_res => null(), map => null(), sig => null(), res => null()



    ! Calculate initial chisq
    todo = .true.
    c => compList
    chisq_old = 0d0
    do while (associated(c))
       select type (c)
       type is (comm_md_comp)
         continue
       class is (comm_diffuse_comp)
         if (todo) then
           if (allocated(c%indmask)) then
             call compute_chisq(c%comm, chisq_fullsky=chisq_old, mask=c%indmask)
           end if
           todo = .false.
         end if
       end select
       c => c%next()
    end do

    if (cpar%myid_chain .eq. 0) then
      write(*,*) '| Old chisq is ', chisq_old
    end if

    ! Propose spectral index MH step - for a given component or all components?
    ! Let's assume we are doing this akin to the sample_specind_local arguments.


    ! From the sample_nonlin_params section
    c => compList
    do while (associated(c))
       if (c%npar == 0) then
          c => c%next()
          cycle
       end if
                    
       do j = 1, c%npar
          if (c%p_gauss(2,j) == 0.d0) then
            cycle
          end if
          select type (c)
          class is (comm_diffuse_comp)


            if (cpar%myid_chain .eq. 0) then
              write(*,*) trim(c%label), j
              do pol = 1, 3
                 c%theta_pixreg_buff(1,pol,j) = c%theta_pixreg(1,pol,j)
                 c%theta_pixreg(1,pol,j) = c%theta_pixreg(1,pol,j) + rand_gauss(handle) * c%p_gauss(2,j)
              end do
              write(*,*) 'theta_pixreg original', c%theta_pixreg_buff(1,:,j)
              write(*,*) 'theta_pixreg proposal', c%theta_pixreg(1,:,j)
            end if


            call mpi_bcast(c%theta_pixreg_buff, size(c%theta_pixreg_buff), MPI_DOUBLE_PRECISION, &
              & 0, data(1)%info%comm, ierr)
            call mpi_bcast(c%theta_pixreg,      size(c%theta_pixreg),      MPI_DOUBLE_PRECISION, &
              & 0, data(1)%info%comm, ierr)

            !do pol = 1,3
               do pix = 0, c%theta(j)%p%info%np-1
                  !if (cpar%myid_chain==0) write(*,*) c%ind_pixreg_arr(pix,pol,j), c%theta(j)%p%map(pix,1), cpar%myid_chain, pol, pix
                  c%theta(j)%p%map(pix,1) = c%theta_pixreg(1,c%ind_pixreg_arr(pix,1,j),j)
               end do
            !end do
            !call mpi_allreduce(theta_smooth%info%nalm, nalm_tot_reg, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)

            c%theta(j)%p%map(:,1) = c%theta_pixreg(1,1,j)


            call c%theta(j)%p%YtW()


          end select

       end do
       
       !go to next component
       c => c%next()
           
    end do

    ! Update mixing matrices
    c => compList
    do while (associated(c))
       call c%updateMixmat
       c => c%next()
    end do

    ! Perform component separation
    call timer%start(TOT_AMPSAMP)
    do samp_group = 1, cpar%cg_num_user_samp_groups
       if (cpar%myid_chain == 0) then
          write(*,fmt='(a,i4,a,i4,a,i4)') ' |  Chain = ', cpar%mychain, ' -- CG sample group = ', &
               & samp_group, ' of ', cpar%cg_num_user_samp_groups
       end if
       call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

       if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

    end do
    call timer%stop(TOT_AMPSAMP)


    todo = .true.
    chisq_prop = 0d0
    c => compList
    do while (associated(c))
       select type (c)
       type is (comm_md_comp)
         continue
       class is (comm_diffuse_comp)
         if (todo) then
           if (allocated(c%indmask)) then
             call compute_chisq(c%comm, chisq_fullsky=chisq_prop, mask=c%indmask)
           end if
           todo = .false.
         end if
       end select
       c => c%next()
    end do

    if (cpar%myid_chain .eq. 0) then
      write(*,*) "|    Proposal chisq is ", chisq_prop
      write(*,*) "|    Delta chi^2 is    ", chisq_prop - chisq_old
    end if

    ! Check MH statistic
    reject = log(rand_uni(handle)) > (chisq_old - chisq_prop)/2
    call mpi_bcast(reject, 1, MPI_LOGICAL, 0, data(1)%info%comm, ierr)


    if (reject) then
      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| MH step rejected, returning to original spectral indices.'
        write(*,*) '| Rejected chisq is ', chisq_old
        write(*,*) '| '
      end if

      c => compList
      do while (associated(c))
         if (c%npar == 0) then
            c => c%next()
            cycle
         end if
                      
         do j = 1, c%npar
            if (c%p_gauss(2,j) == 0.d0) then
              cycle
            end if
            select type (c)
            class is (comm_diffuse_comp)

              if (cpar%myid_chain .eq. 0) then
                do k = 1, 3
                   c%theta_pixreg(1,k,j) = c%theta_pixreg_buff(1,k,j)
                end do
              end if

              call mpi_bcast(c%theta_pixreg,      size(c%theta_pixreg),      MPI_DOUBLE_PRECISION, &
                & 0, data(1)%info%comm, ierr)

              c%theta(j)%p%map(:,1) = c%theta_pixreg_buff(1,1,j)


              call c%theta(j)%p%YtW()
            end select

         end do
         
         !go to next component
         c => c%next()
             
      end do

      ! Update mixing matrices
      c => compList
      do while (associated(c))
         call c%updateMixmat
         c => c%next()
      end do

      ! Perform component separation
      call timer%start(TOT_AMPSAMP)
      do samp_group = 1, cpar%cg_num_user_samp_groups
         if (cpar%myid_chain == 0) then
            write(*,fmt='(a,i4,a,i4,a,i4)') ' |  Chain = ', cpar%mychain, ' -- CG sample group = ', &
                 & samp_group, ' of ', cpar%cg_num_user_samp_groups
         end if
         call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

         if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

      end do
      call timer%stop(TOT_AMPSAMP)

      todo = .true.
      c => compList
      chisq_old = 0d0
      do while (associated(c))
         select type (c)
         type is (comm_md_comp)
           continue
         class is (comm_diffuse_comp)
           if (todo) then
             if (allocated(c%indmask)) then
               call compute_chisq(c%comm, chisq_fullsky=chisq_old, mask=c%indmask)
             end if
             todo = .false.
           end if
         end select
         c => c%next()
      end do

      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| Current MH chisq is ', chisq_old
        write(*,*) '| '
      end if

    else
      if (cpar%myid_chain == 0) then
        write(*,*) '| '
        write(*,*) '| MH step accepted'
        write(*,*) '| Current MH chisq is ', chisq_prop
      end if
    end if

  end subroutine sample_specind_mh_sample

end module comm_mh_specind_mod
