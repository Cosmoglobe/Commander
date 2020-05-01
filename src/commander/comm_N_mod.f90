module comm_N_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_N, compute_invN_lm
  
  type, abstract :: comm_N
     ! Data variables
     character(len=512)       :: type
     integer(i4b)             :: nside, nside_chisq_lowres, nmaps, np, npix, myid, comm, nprocs
     logical(lgt)             :: pol, pol_only
     logical(lgt)             :: set_noise_to_mean
     character(len=512)       :: cg_precond
     real(dp)                 :: uni_fsky
     real(dp), allocatable, dimension(:) :: alpha_nu ! (T,Q,U)
     class(comm_map),     pointer :: invN_diag => null()
     class(comm_map),     pointer :: rms_reg   => null()
     class(comm_mapinfo), pointer :: info      => null()
   contains
     ! Data procedures
     procedure(matmulInvN),       deferred :: invN
     procedure(matmulInvNlowres), deferred :: invN_lowres
     procedure(matmulN),          deferred :: N
     procedure(matmulSqrtInvN),   deferred :: sqrtInvN
     procedure(returnRMS),        deferred :: rms
     procedure(returnRMSpix),     deferred :: rms_pix
     procedure(update_N),         deferred :: update_N
  end type comm_N

  abstract interface
     ! Return map_out = invN * map
     subroutine matmulInvN(self, map)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
     end subroutine matmulInvN

     ! Return map_out = invN * map
     subroutine matmulInvNlowres(self, map)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
     end subroutine matmulInvNlowres

     ! Return map_out = N * map
     subroutine matmulN(self, map)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
     end subroutine matmulN

     ! Return map_out = sqrtInvN * map
     subroutine matmulSqrtInvN(self, map)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
     end subroutine matmulSqrtInvN

     ! Return rms map
     subroutine returnRMS(self, res)
       import comm_map, comm_N, dp
       implicit none
       class(comm_N),   intent(in)    :: self
       class(comm_map), intent(inout) :: res
     end subroutine returnRMS

     ! Return rms map
     function returnRMSpix(self, pix, pol)
       import i4b, comm_N, dp
       implicit none
       class(comm_N),   intent(in)    :: self
       integer(i4b),    intent(in)    :: pix, pol
       real(dp)                       :: returnRMSpix
     end function returnRMSpix

     ! Update noise model
     subroutine update_N(self, info, handle, mask, regnoise, procmask, noisefile, map)
       import comm_N, comm_mapinfo, comm_map, dp, planck_rng
       implicit none
       class(comm_N),                      intent(inout)          :: self
       type(planck_rng),                   intent(inout)          :: handle
       class(comm_mapinfo),                intent(in)             :: info
       class(comm_map),                    intent(in),   optional :: mask
       real(dp),         dimension(0:,1:), intent(out),  optional :: regnoise
       class(comm_map),                    intent(in),   optional :: procmask
       character(len=*),                   intent(in),   optional :: noisefile
       class(comm_map),                    intent(in),   optional :: map
     end subroutine update_N

  end interface

contains

  subroutine compute_invN_lm(invN_diag)
    implicit none

    class(comm_map),  intent(inout) :: invN_diag

    real(dp)     :: l0min, l0max, l1min, l1max, npix, t1, t2
    integer(i4b) :: j, l, m, lp, l2, ier, twolmaxp2, pos, lmax
    complex(dpc) :: val(invN_diag%info%nmaps)
    real(dp), allocatable, dimension(:)   :: threej_symbols, threej_symbols_m0
    real(dp), allocatable, dimension(:,:) :: N_lm, a_l0

    lmax      = invN_diag%info%lmax
    twolmaxp2 = 2*lmax+2
    npix      = real(invN_diag%info%npix,dp)

    allocate(N_lm(0:invN_diag%info%nalm-1,invN_diag%info%nmaps))
    allocate(a_l0(0:lmax,invN_diag%info%nmaps))
    call invN_diag%YtW_scalar

    if (invN_diag%info%myid == 0) then
       ! Seriously ugly hack. Should figure out how to compute N_{lm,l'm'} properly
       ! for polarization
       a_l0(:,:) = invN_diag%alm(0:lmax,:)
    end if
    call mpi_bcast(a_l0, size(a_l0), MPI_DOUBLE_PRECISION, 0, invN_diag%info%comm, ier)

    call wall_time(t1)
    !$OMP PARALLEL PRIVATE(pos,j,m,l,threej_symbols_m0,threej_symbols,ier,val,l2,lp,l0min,l0max,l1min,l1max)
    allocate(threej_symbols(twolmaxp2))
    allocate(threej_symbols_m0(twolmaxp2))
    !$OMP DO SCHEDULE(guided)
    do j = 1, invN_diag%info%nm
       m = invN_diag%info%ms(j)
       do l = m, lmax
       
          call DRC3JJ(real(l,dp), real(l,dp), 0.d0, 0.d0, l0min, l0max, &
               & threej_symbols_m0, twolmaxp2, ier)
          call DRC3JJ(real(l,dp), real(l,dp), real(-m,dp), real(m,dp), l1min, l1max, &
               & threej_symbols, twolmaxp2, ier)
             
          val = dcmplx(0.d0,0.d0)
          do l2 = 1, nint(l1max-l1min)+1
             lp = nint(l1min) + l2 - 1
             if (lp > lmax) exit
             ! Only m3 = 0 contributes, because m2 = -m1 => m3 = 0
             val = val + a_l0(lp,:) * sqrt(2.d0*lp+1.d0) * &
                  & threej_symbols(l2) * threej_symbols_m0(lp-nint(l0min)+1)
          end do
          val = val * (2*l+1) / sqrt(4.d0*pi) * npix / (4.d0*pi)
          if (mod(m,2)/=0) val = -val

          call invN_diag%info%lm2i(l,m,pos)
          if (m == 0) then
             N_lm(pos,:)   = real(val,dp)
          else
             N_lm(pos,:)   = real(val,dp)
             N_lm(pos+1,:) = real(val,dp)
          end if
       end do
    end do
    !$OMP END DO
    deallocate(threej_symbols, threej_symbols_m0)
    !$OMP END PARALLEL
    call wall_time(t2)

    invN_diag%alm = N_lm
    !write(*,*) sum(abs(invN_diag%alm))

    deallocate(N_lm, a_l0)
    
  end subroutine compute_invN_lm
  
end module comm_N_mod
