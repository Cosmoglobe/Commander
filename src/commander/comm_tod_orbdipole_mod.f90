module comm_tod_orbdipole_mod
  use comm_tod_mod
  use comm_utils
  implicit none


contains

  ! Compute map with white noise assumption from correlated noise 
  ! corrected and calibrated data, d' = (d-n_corr-n_temp)/gain 
  subroutine compute_orbital_dipole(tod, ind, pix, psi, s_orb)
    implicit none
    class(comm_tod),                     intent(in)  :: tod
    integer(i4b),                        intent(in)  :: ind !scan nr/index
    integer(i4b),        dimension(:,:), intent(in)  :: pix, psi
    real(sp),            dimension(:,:), intent(out) :: s_orb
    integer(i4b)         :: i, j, p
    real(dp)             :: x, T_0, q, pix_dir(3), b, b_dot, summation
    real(dp)             :: theta, phi, psi_d
    real(dp), dimension(3,3) :: rot_mat
    real(dp), dimension(3)   :: vnorm
    real(dp), parameter  :: h = 6.62607015d-34   ! Planck's constant [Js]

    !T_0 = T_CMB*k_b/h                           ! T_0 = T_CMB frequency
    b = sqrt(sum(tod%scans(ind)%v_sun**2))/c   ! beta for the given scan

    !these are the npipe paper definitions
    !TODO: maybe also use the bandpass shift to modify the central frequency? 
    !will that matter?
    x = h * tod%central_freq/(k_B * T_CMB)
    q = (x/2.d0)*(exp(x)+1)/(exp(x) -1)

    !if (trim(tod%label(1)) == '27M') open(58,file='orb_27M.dat')
    do i = 1,tod%ndet 
       if (.not. tod%scans(ind)%d(i)%accept) cycle
       do j=1,tod%scans(ind)%ntod !length of the tod
          b_dot = dot_product(tod%scans(ind)%v_sun, tod%pix2vec(:,pix(j,i)))/c
          !s_orb(j,i) = real(T_CMB  * b_dot,sp) !* tod%mb_eff(i) !only dipole, 1.d6 to make it uK, as [T_CMB] = K
          !s_orb(j,i) = T_CMB  * 1.d6 * b_dot !only dipole, 1.d6 to make it uK, as [T_CMB] = K
          !s_orb(j,i) = T_CMB * 1.d6 * (b_dot + q*b_dot**2) ! with quadrupole
          s_orb(j,i) = T_CMB * (b_dot + q*(b_dot**2 - b**2/3.)) ! net zero monopole

          !TODO: add sl contribution to orbital dipole here
          call pix2ang_ring(tod%info%nside, pix(j,i), theta, phi)
          !rotate v_sun into frame where pointing is along z axis
          !write(*,*) -phi, -theta, -tod%psi(psi(j,i)), psi(j,i)
          p     = tod%pix2ind(pix(j,i))
          theta = tod%ind2ang(1,p)
          phi   = tod%ind2ang(2,p)
          psi_d = tod%psi(psi(j,i))
          !write(*,*), j, phi, theta, psi_d, rot_mat 
          !call compute_euler_matrix_zyz(-phi, -theta, -psi_d, rot_mat)
          call compute_euler_matrix_zyz(-psi_d, -theta, -phi, rot_mat)
          vnorm = matmul(rot_mat, tod%scans(ind)%v_sun)
          vnorm = vnorm / sum(vnorm**2)
          summation = vnorm(1)*tod%orb_dp_s(i,1)+vnorm(2)*tod%orb_dp_s(i,2)+& 
            & vnorm(3)*tod%orb_dp_s(i,3)+vnorm(1)*vnorm(1)*tod%orb_dp_s(i,4)+&
            & vnorm(1)*vnorm(2)*tod%orb_dp_s(i,5) + vnorm(1)*vnorm(3)* &
            & tod%orb_dp_s(i,6) + vnorm(2)*vnorm(2)*tod%orb_dp_s(i,7) + &
            & vnorm(2)*vnorm(2)*tod%orb_dp_s(i,8) + vnorm(3)*vnorm(3)*&
            & tod%orb_dp_s(i,9) 
          !if (trim(tod%label(i)) == '27M') write(*,*) j, T_CMB *summation, vnorm(1), vnorm(2), vnorm(3), tod%orb_dp_s(i,1), tod%scans(ind)%v_sun
          s_orb(j,i) = s_orb(j,i) + T_CMB *summation

       end do
   end do
   !if (trim(tod%label(1)) == '27M') close(58)

  end subroutine compute_orbital_dipole


end module comm_tod_orbdipole_mod
