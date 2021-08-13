REAL FUNCTION fun_harmonic(K1, K2)
! return the Harmonic average of K1, K2 but adjusted to take into account
!   negative K1, K2 that are used to set a zero diffusivity.
   REAL K1, K2
   REAL, PARAMETER :: toler = 1.e-15
   fun_harmonic = 2. * ABS(K1 * K2) / MAX(toler, K1 + K2)
END FUNCTION fun_harmonic


SUBROUTINE viral_4_eqn(T_new, I1_new, I2_new, V_new, T, I1, I2, V, beta, k, delta_c, p, c, dt, ntime)
! python call looks like:
!    T_new, I1_new, I2_new, V_new = viral_4_eqn(T, I1, I2, V, beta, k, delta_c, p, c, dt, ntime)

! *******************************************************************************************
! This subroutine solve the classical 4 eqn VIRAL model (The model with an eclipse phase) using an iterative explicit
! method and as such the time step size can not be too small e.g. 10 mins or dt = 10. * 60. = 600.
! don't forget there needs to be a seed of an infected person in the E group.
! have a look at eqn 2 of Modeling-the-viral-dynamics-of-SARS-CoV-2-infect.pdf
! *******************************************************************************************

! T_new, I1_new, I2_new, V_new the calculated solution variables after ntime time steps output
! T, I1, I2, V initial solution variables sent down to this subroutine
! beta, k, delta_c, p, c are the coefficients in the VIRAL equations sent down into this subroutine
! dt = time step size
! ntime = no of time steps

   IMPLICIT NONE
! dummy variables
   REAL, INTENT(OUT) :: T_new, I1_new, I2_new, V_new
   REAL, INTENT(IN) :: T, I1, I2, V, beta, k, delta_c, p, c, dt
   INTEGER, INTENT(IN) :: ntime

! local variables
   INTEGER its, nits, itime ! no of non-linear iterations
   REAL T_old, I1_old, I2_old, V_old ! previous time step values of T, I1, I2, V

! execution
   T_new = T; I1_new = I1; I2_new = I2; V_new = V
   nits = 3

   DO itime = 1, ntime
      T_old = T_new; I1_old = I1_new; I2_old = I2_new; V_old = V_new
      
      DO its = 1, nits
         T_new = T_old + dt * (-beta * V_new * T_new)
         I1_new = I1_old + dt * (beta * V_new * T_new - k * I1_new)
         I2_new = I2_old + dt * (k * I1_new - delta_c * I2_new) 
         V_new = V_old + dt * (p * I2_new - c * V_new)
      END DO ! DO its = 1, nits
   END DO
   RETURN
END SUBROUTINE viral_4_eqn


SUBROUTINE viral_4_eqn_secondary_target(T1_new, T2_new, I_new, V_new, T1, T2, I, V, beta, lambda, delta_c, w, p, c, dt, ntime)
! python call looks like:
!    T_new, I1_new, I2_new, V_new = viral_4_eqn_secondary_target(T1, T2, I, V, beta, lambda, delta_c, w, p, c, dt, ntime)

! *******************************************************************************************
! This subroutine solve the classical 4 eqn VIRAL model based on a secondary target using an iterative explicit
! method and as such the time step size can not be too small e.g. 10 mins or dt = 10. * 60. = 600.
! don't forget there needs to be a seed of an infected person in the E group.
! have a look at eqn 3 of Modeling-the-viral-dynamics-of-SARS-CoV-2-infect.pdf
! *******************************************************************************************

! T1_new, T2_new, I_new, V_new the calculated solution variables after ntime time steps output
! T1, T2, I, V initial solution variables sent down to this subroutine
! beta, lambda, delta_c, w, p, c are the coefficients in the VIRAL equation sent down into this subroutine
! dt = time step size
! ntime = no of time steps.

   IMPLICIT NONE
! dummy variables
   REAL, INTENT(OUT) :: T1_new, T2_new, I_new, V_new
   REAL, INTENT(IN) :: T1, T2, I, V, beta, lambda, delta_c, w, p, c, dt
   INTEGER, INTENT(in) :: ntime

! local variables
   INTEGER its, nits, itime ! no of non-linear iterations
   REAL T1_old, T2_old, I_old, V_old ! previous time step values of T1, T2, I, V

! execution
   T1_new = T1; T2_new = T2; I_new = I; V_new = V
   nits = 3

   DO itime = 1, ntime
      T1_old = T1_new; T2_old = T2_new; I_old = I_new; V_old = V_new
      
      DO its = 1, nits
         T1_new = T1_old + dt * (-beta * V_new * T1_new)
         T2_new = T2_old + dt * (lambda - beta * V_new * T2_new)
         I_new = I_old + dt * (beta * V_new * (T1_new + T2_new) - (delta_c + w * T2_new) * I_new)
         V_new = V_old + dt * (p * I_new - c * V_new)
      END DO ! DO its = 1, nits
   END DO
   RETURN
END SUBROUTINE viral_4_eqn_secondary_target


SUBROUTINE sim_time_steping_diffusion_calc(T_new, T_initial, ntime, nits, nits_solv_ng, relax, &
   error_solv, error_solv_ng, dx, dy, dz, dt, V_n, sigmat, sigma_s_off, kdiff, S, U, I_upwind, &
   I_harmonic, ndim_vel, nx, ny, nz, ng, ng2) 
! python call looks like:
!    T_new = sim_time_steping_diffusion_calc(T_initial, ntime, nits, nits_solv_ng, relax, error_solv,
!                                            error_solv_ng, dx, dy, dz, dt, V_n, sigmat, sigma_s_off,
!                                            kdiff, S, U, I_upwind, I_harmonic, ndim_vel, nx, ny, nz, ng, ng2)

! *******************************************************************************************
! A simplified version of the `time_steping_diffusion_calc` for not so many groups. 
! *******************************************************************************************

! ntime = no of time steps
! nits = no of non-linear iterations
! dt = time step size (can use dt = 1.e15 for steady state cals) 
! NITS,NITS_NG are the max no of single group iterations and overall iterations. 
! ERROR_SOLV,ERROR_SOLV_NG are the error tolerances of the solver of single group iterations and overall iterations.
! ng2 = 0 if ncol_sigma_s .ne. 0

   IMPLICIT NONE
! dummy variables
   INTEGER, INTENT(IN) :: I_upwind, I_harmonic, ndim_vel, nx, ny, nz, ng, ng2
   INTEGER, INTENT(IN) :: ntime, nits, nits_solv_ng
   REAL, INTENT(IN) :: relax, error_solv, error_solv_ng, dx, dy, dz, dt, V_n(ng)
   REAL, INTENT(OUT) :: T_new(nx, ny, nz, ng)
   REAL, INTENT(IN) :: T_initial(nx, ny, nz, ng)
   REAL, INTENT(IN) :: sigmat(nx, ny, nz, ng), sigma_s_off(nx, ny, nz, ng, ng2), kdiff(nx, ny, nz, ng)
   REAL, INTENT(IN) :: S(nx, ny, nz, ng), U(ndim_vel, nx + 1, ny + 1, nz + 1, ng)

! local variables
   INTEGER :: ncol_sigma_s
   REAL, ALLOCATABLE :: sigma_s_off_pter(:, :, :, :) 
   INTEGER, ALLOCATABLE :: fin_sigma_s(:), col_sigma_s(:)

! execution
   ncol_sigma_s = 0
   ALLOCATE(fin_sigma_s(ng + 1))
   ALLOCATE(col_sigma_s(ncol_sigma_s))
   ALLOCATE(sigma_s_off_pter(nx, ny, nz, ncol_sigma_s))
   fin_sigma_s = 0

   CALL time_steping_diffusion_calc(T_new, T_initial, ntime, nits, nits_solv_ng, relax, error_solv, &
      error_solv_ng, dx, dy, dz, dt, V_n, sigmat, sigma_s_off, sigma_s_off_pter, fin_sigma_s,       &
      col_sigma_s, kdiff, S, U, I_upwind, I_harmonic, ndim_vel, ncol_sigma_s, nx, ny, nz, ng, ng2)
   RETURN 
END SUBROUTINE sim_time_steping_diffusion_calc


SUBROUTINE time_steping_diffusion_calc(T_new, T_initial, ntime, nits, nits_solv_ng, relax, error_solv,  &
   error_solv_ng, dx, dy, dz, dt, V_n, sigmat, sigma_s_off, sigma_s_off_pter, fin_sigma_s, col_sigma_s, &
   kdiff, S, U, I_upwind, I_harmonic, ndim_vel, ncol_sigma_s, nx, ny, nz, ng, ng2)
! python call looks like:
!    T_new = time_steping_diffusion_calc(T_initial, ntime, nits, nits_solv_ng, relax, error_solv,
!                                        error_solv_ng, dx, dy, dz, dt, V_n, sigmat, sigma_s_off,
!                                        sigma_s_off_pter, fin_sigma_s, col_sigma_s, kdiff, S, U,
!                                        I_upwind, I_harmonic, ndim_vel, ncol_sigma_s, nx, ny, nz, ng, ng2)

! *******************************************************************************************
! *******************************************************************************************

   IMPLICIT NONE
! dummy variables
   INTEGER, INTENT(IN) :: I_upwind, I_harmonic, ndim_vel, ncol_sigma_s, nx, ny, nz, ng, ng2
   INTEGER, INTENT(IN) :: ntime, nits, nits_solv_ng
   REAL, INTENT(IN) :: relax, error_solv, error_solv_ng, dx, dy, dz, dt, V_n(ng)
   REAL, INTENT(OUT) :: T_new(nx, ny, nz, ng)
   REAL, INTENT(IN) :: T_initial(nx, ny, nz, ng)
   REAL, INTENT(IN) :: sigmat(nx, ny, nz, ng), sigma_s_off(nx, ny, nz, ng, ng2)
   REAL, INTENT(IN) :: sigma_s_off_pter(nx, ny, nz, ncol_sigma_s)
   INTEGER, INTENT(IN) :: fin_sigma_s(ng + 1), col_sigma_s(ncol_sigma_s)
   REAL, INTENT(IN) :: kdiff(nx, ny, nz, ng), S(nx, ny, nz, ng), U(ndim_vel, nx + 1, ny + 1, nz + 1, ng)

! Local variables...
   REAL, PARAMETER :: relax_ng = 1. ! Multi-group relaxation term.
   REAL, PARAMETER :: toler = 1.e-10, xi_limit = 2. ! defines TVD curve on the NVD diagram.
   LOGICAL upwind1st ! apply 1st order upwinding to advection.
   LOGICAL harmonic_ave ! perform a Harmonic average of the diffusion along the CV boundary
   INTEGER, PARAMETER :: iSwitch_FBGS = 1 ! Switch on FBGS (ISWITCH_FBGS=1 switched on FBGS and ISWITCH_FBGS=0 used forward block Gauss Siedel.)
   REAL, ALLOCATABLE :: T_old(:, :, :, :), A_mat_i_j_k(:, :, :, :)
   REAL, ALLOCATABLE :: A_mat_ip1_j_k(:, :, :, :), A_mat_im1_j_k(:, :, :, :)
   REAL, ALLOCATABLE :: A_mat_i_jp1_k(:, :, :, :), A_mat_i_jm1_k(:, :, :, :)
   REAL, ALLOCATABLE :: A_mat_i_j_kp1(:, :, :, :), A_mat_i_j_km1(:, :, :, :)
   REAL, ALLOCATABLE :: RHS_scat(:, :, :), T_new_temp_ng(:, :, :), S_inc_lim(:, :, :, :)
   INTEGER itime, its, its_ng, ii, i, jj, j, kk, k, ig, jg
   REAL diff_i, R_T_new, R_T_max_dif, R_T_max_dif_ng
   REAL A_i_j_k, A_ip1_j_k, A_im1_j_k, A_i_jp1_k, A_i_jm1_k, A_i_j_kp1, A_i_j_km1
   INTEGER iigg, ig_start, ig_finish, ig_step
   INTEGER switch_up1(3), switch_up2(3), im, jm, km, iDim, iSwitch, jSwitch, kSwitch
   INTEGER i1, i2, i3, j1, j2, j3, k1, k2, k3
   INTEGER ii1, ii2, ii3, jj1, jj2, jj3, kk1, kk2, kk3, count, iLast
   REAL denoIn1, denoIn2, ctilIn1, ctilIn2, ftilIn1, ftilIn2, T_limit1, T_limit2, U_dim(3)
   REAL RH_sig_u, T_upwind(3), T_limit(3), H_sig_u(3)

! function
   REAL fun_harmonic

   upwind1st = (I_upwind == 1)
   harmonic_ave = (I_harmonic == 1)
   ALLOCATE(T_old(nx, ny, nz, ng), A_mat_i_j_k(nx, ny, nz, ng),     &
      A_mat_ip1_j_k(nx, ny, nz, ng), A_mat_im1_j_k(nx, ny, nz, ng), &
      A_mat_i_jp1_k(nx, ny, nz, ng), A_mat_i_jm1_k(nx, ny, nz, ng), &
      A_mat_i_j_kp1(nx, ny, nz, ng), A_mat_i_j_km1(nx, ny, nz, ng), &
      RHS_scat(nx, ny, nz), T_new_temp_ng(nx, ny, nz), S_inc_lim(nx, ny, nz, ng))
   ! print *,'here1'

! Define matricies
   DO ig = 1, ng
      DO k = 2, nz - 1
         DO j = 2, ny - 1
            DO i = 2, nx - 1
               IF (harmonic_ave) THEN
                  ! max{(D_{i-1,j}*D_{I,j})/(D_{i-1,j}+D_{I,j}),0}
                  diff_i = .5 * (fun_harmonic(kdiff(i, j, k, ig), kdiff(i - 1, j, k, ig)) + &
                                 fun_harmonic(kdiff(i, j, k, ig), kdiff(i + 1, j, k, ig))) / (dx ** 2) &
                         + .5 * (fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j - 1, k, ig)) + &
                                 fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j + 1, k, ig))) / (dy ** 2) & 
                         + .5 * (fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j, k - 1, ig)) + &
                                 fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j, k + 1, ig))) / (dz ** 2)

                  A_mat_i_j_k(i, j, k, ig) = diff_i + sigmat(i, j, k, ig) + (1. / (dt * V_n(ig)))
                  A_mat_ip1_j_k(i, j, k, ig) = -.5 * fun_harmonic(kdiff(i, j, k, ig), kdiff(i + 1, j, k, ig)) / (dx ** 2)
                  A_mat_im1_j_k(i, j, k, ig) = -.5 * fun_harmonic(kdiff(i, j, k, ig), kdiff(i - 1, j, k, ig)) / (dx ** 2)
                  A_mat_i_jp1_k(i, j, k, ig) = -.5 * fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j + 1, k, ig)) / (dy ** 2)
                  A_mat_i_jm1_k(i, j, k, ig) = -.5 * fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j - 1, k, ig)) / (dy ** 2)
                  A_mat_i_j_kp1(i, j, k, ig) = -.5 * fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j, k + 1, ig)) / (dz ** 2)
                  A_mat_i_j_km1(i, j, k, ig) = -.5 * fun_harmonic(kdiff(i, j, k, ig), kdiff(i, j, k - 1, ig)) / (dz ** 2)
               ELSE ! simple mean...
                  diff_i = .5 * (MAX(kdiff(i, j, k, ig) + kdiff(i - 1, j, k, ig), 0.0) + &
                                 MAX(kdiff(i, j, k, ig) + kdiff(i + 1, j, k, ig), 0.0)) / (dx ** 2) &
                         + .5 * (MAX(kdiff(i, j, k, ig) + kdiff(i, j - 1, k, ig), 0.0) + &
                                 MAX(kdiff(i, j, k, ig) + kdiff(i, j + 1, k, ig), 0.0)) / (dy ** 2) &
                         + .5 * (MAX(kdiff(i, j, k, ig) + kdiff(i, j, k - 1, ig), 0.0) + &
                                 MAX(kdiff(i, j, k, ig) + kdiff(i, j, k + 1, ig), 0.0)) / (dz ** 2)

                  A_mat_i_j_k(i, j, k, ig) = diff_i + sigmat(i, j, k, ig) + (1. / (dt * V_n(ig))) 
                  A_mat_ip1_j_k(i, j, k, ig) = -.5 * MAX(kdiff(i, j, k, ig) + kdiff(i + 1, j, k, ig), 0.0) / (dx ** 2)
                  A_mat_im1_j_k(i, j, k, ig) = -.5 * MAX(kdiff(i, j, k, ig) + kdiff(i - 1, j, k, ig), 0.0) / (dx ** 2)
                  A_mat_i_jp1_k(i, j, k, ig) = -.5 * MAX(kdiff(i, j, k, ig) + kdiff(i, j + 1, k, ig), 0.0) / (dy ** 2)
                  A_mat_i_jm1_k(i, j, k, ig) = -.5 * MAX(kdiff(i, j, k, ig) + kdiff(i, j - 1, k, ig), 0.0) / (dy ** 2)
                  A_mat_i_j_kp1(i, j, k, ig) = -.5 * MAX(kdiff(i, j, k, ig) + kdiff(i, j, k + 1, ig), 0.0) / (dz ** 2)
                  A_mat_i_j_km1(i, j, k, ig) = -.5 * MAX(kdiff(i, j, k, ig) + kdiff(i, j, k - 1, ig), 0.0) / (dz ** 2)
               ENDIF
            END DO
         END DO
      END DO
   END DO

! Discretize advection...
   IF (ndim_vel .NE. 0) THEN
      IF ((nz == 3) .AND. (ndim_vel == 3)) THEN
         PRINT *, 'incompatable nz and ndim_vel'
         STOP 292
      ENDIF

      U_dim = 0.0; T_upwind = 0.0; T_limit = 0.0; H_sig_u = 0.0; S_inc_lim = 0.0

      DO ig = 1, ng
         DO k = 2, nz
            switch_up1(3) = 1 - MIN(1, ABS(k - 2))
            switch_up2(3) = 1 - MIN(1, ABS(k - nz))
            km = k + switch_up1(3) - switch_up2(3)
            DO j = 2, ny
               switch_up1(2) = 1 - MIN(1, ABS(j - 2))
               switch_up2(2) = 1 - MIN(1, ABS(j - ny))
               jm = j + switch_up1(2) - switch_up2(2)
               DO i = 2, nx
                  ! advection upwards
                  ! Now adjust to use upwind neer the boundaries...
                  switch_up1(1) = 1 - MIN(1, ABS(i - 2))
                  switch_up2(1) = 1 - MIN(1, ABS(i - nx))
                  im = i + switch_up1(1) - switch_up2(1)
                  DO iDim = 1, ndim_vel
                     ISWITCH = 1 - MIN(1,ABS(IDIM-1)); JSWITCH = 1 - MIN(1,ABS(IDIM-2)); KSWITCH = 1 - MIN(1,ABS(IDIM-3))

                     I1 =IM - 2*ISWITCH; I2 =IM - ISWITCH; I3 =IM
                     II1=IM + ISWITCH;   II2=IM;           II3=IM - ISWITCH

                     J1 =JM - 2*JSWITCH; J2 =JM - JSWITCH; J3 =JM
                     JJ1=JM + JSWITCH;   JJ2=JM;           JJ3=JM - JSWITCH

                     K1 =KM - 2*KSWITCH; K2 =KM - KSWITCH; K3 =KM
                     KK1=KM + KSWITCH;   KK2=KM;           KK3=KM - KSWITCH

                     DENOIN1 = SIGN( MAX(ABS( T_NEW(I3,J3,K3,IG)  - T_NEW(I1,J1,K1,IG) ), TOLER),  &
                        T_NEW(I3,J3,K3,IG)  - T_NEW(I1,J1,K1,IG) )
                     DENOIN2 = SIGN( MAX(ABS( T_NEW(II3,JJ3,KK3,IG) - T_NEW(II1,JJ1,KK1,IG) ), TOLER), &
                        T_NEW(II3,JJ3,KK3,IG) - T_NEW(II1,JJ1,KK1,IG) )
                     CTILIN1 = ( T_NEW(I2,J2,K2,IG) - T_NEW(I1,J1,K1,IG) ) / DENOIN1
                     CTILIN2 = ( T_NEW(II2,JJ2,KK2,IG) - T_NEW(II1,JJ1,KK1,IG) ) / DENOIN2

                     FTILIN1 = ( 0.5*(T_NEW(I2,J2,K2,IG)+T_NEW(I3,J3,K3,IG))   - T_NEW(I1,J1,K1,IG) ) / DENOIN1 ! the mean with 0.5 coeff is the high order value of T.
                     FTILIN2 = ( 0.5*(T_NEW(II2,JJ2,KK2,IG)+T_NEW(II3,JJ3,KK3,IG)) - T_NEW(II1,JJ1,KK1,IG) ) / DENOIN2 ! the mean with 0.5 coeff is the high order value of T.

                     ! T_LIMIT(I+1) =          T(I-1) + MAX(  MIN(FTILIN, XI_LIMIT*CTILIN, 1.0), CTILIN) * DENOIN
                     ! T_LIMIT1(I+1) =          T(I-1) + MAX(  MIN(FTILIN, XI_LIMIT*CTILIN, 1.0), CTILIN) * ( T(I+1) - T(I-1) )
                     T_LIMIT1 =  T_NEW(I1,J1,K1,IG)    + MAX(  MIN(FTILIN1, XI_LIMIT*CTILIN1, 1.0), CTILIN1) &
                                 * ( T_NEW(I3,J3,K3,IG)    - T_NEW(I1,J1,K1,IG) )
                     T_LIMIT2 =  T_NEW(II1,JJ1,KK1,IG) + MAX(  MIN(FTILIN2, XI_LIMIT*CTILIN2, 1.0), CTILIN2) &
                                 * ( T_NEW(II3,JJ3,KK3,IG) - T_NEW(II1,JJ1,KK1,IG) )

                     U_DIM(IDIM)=U(IDIM,I,J,K,IG)
                     RH_SIG_U = 0.5 + 0.5*SIGN( 1.0, U_DIM(IDIM) )
                     T_UPWIND(IDIM) = RH_SIG_U * T_NEW(I-ISWITCH,J-JSWITCH,K-KSWITCH,IG)  +  (1.0-RH_SIG_U) * T_NEW(I,J,K,IG) 
                     T_LIMIT(IDIM)  = RH_SIG_U * ( T_UPWIND(IDIM)*REAL(switch_up1(IDIM)) +T_LIMIT1*(1.0-REAL(switch_up1(IDIM))) ) &
                        +  (1.0-RH_SIG_U) * ( T_UPWIND(IDIM)*REAL(switch_up2(IDIM)) +T_LIMIT2*(1.0-REAL(switch_up2(IDIM))))

                     ! T_COEF(IDIM) = U_DIM(IDIM)*T_LIMIT/ ( SIGN(1.0,T_UPWIND)* MAX(TOLER, ABS(T_UPWIND) ) )
                     H_SIG_U(IDIM)=RH_SIG_U
                  END DO ! DO IDIM = 1, NDIM_VEL

                  ! A_MAT(I,J,K,IG,0,0,0)  = A_MAT(I,J,K,IG,0,0,0)   &
                  !          - (1.-H_SIG_U(1))*U_DIM(1)/DX - (1.-H_SIG_U(2))*U_DIM(2)/DY - (1.-H_SIG_U(3))*U_DIM(3)/DZ
                  ! A_MAT(I,J,K,IG,-1,0,0) = A_MAT(I,J,K,IG,-1,0,0) - H_SIG_U(1)*U_DIM(1)/DX
                  ! A_MAT(I,J,K,IG,0,-1,0) = A_MAT(I,J,K,IG,0,-1,0) - H_SIG_U(2)*U_DIM(2)/DY
                  ! A_MAT(I,J,K,IG,0,0,-1) = A_MAT(I,J,K,IG,0,0,-1) - H_SIG_U(3)*U_DIM(3)/DZ
         
                  ! A_MAT(I-1,J,K,IG,0,0,0)= A_MAT(I-1,J,K,IG,0,0,0)+ H_SIG_U(1)*U_DIM(1)/DX
                  ! A_MAT(I-1,J,K,IG,1,0,0)= A_MAT(I-1,J,K,IG,1,0,0)+ (1.-H_SIG_U(1))*U_DIM(1)/DX

                  ! A_MAT(I,J-1,K,IG,0,0,0)= A_MAT(I,J-1,K,IG,0,0,0)+ H_SIG_U(2)*U_DIM(2)/DY
                  ! A_MAT(I,J-1,K,IG,0,1,0)= A_MAT(I,J-1,K,IG,0,1,0)+ (1.-H_SIG_U(2))*U_DIM(2)/DY

                  ! A_MAT(I,J,K-1,IG,0,0,0)= A_MAT(I,J,K-1,IG,0,0,0)+ H_SIG_U(3)*U_DIM(3)/DZ
                  ! A_MAT(I,J,K-1,IG,0,0,1)= A_MAT(I,J,K-1,IG,0,0,1)+ (1.-H_SIG_U(3))*U_DIM(3)/DZ

                  A_MAT_I_J_K(I,J,K,IG)   = A_MAT_I_J_K(I,J,K,IG)    &
                        - (1.-H_SIG_U(1))*U_DIM(1)/DX - (1.-H_SIG_U(2))*U_DIM(2)/DY - (1.-H_SIG_U(3))*U_DIM(3)/DZ
                  A_MAT_IM1_J_K(I,J,K,IG) =  A_MAT_IM1_J_K(I,J,K,IG) - H_SIG_U(1)*U_DIM(1)/DX
                  A_MAT_I_JM1_K(I,J,K,IG) =  A_MAT_I_JM1_K(I,J,K,IG) - H_SIG_U(2)*U_DIM(2)/DY
                  A_MAT_I_J_KM1(I,J,K,IG) =  A_MAT_I_J_KM1(I,J,K,IG) - H_SIG_U(3)*U_DIM(3)/DZ
      
                  A_MAT_I_J_K(I-1,J,K,IG)   = A_MAT_I_J_K(I-1,J,K,IG)    + H_SIG_U(1)*U_DIM(1)/DX
                  A_MAT_IP1_J_K(I-1,J,K,IG) =  A_MAT_IP1_J_K(I-1,J,K,IG) + (1.-H_SIG_U(1))*U_DIM(1)/DX

                  A_MAT_I_J_K(I,J-1,K,IG)   = A_MAT_I_J_K(I,J-1,K,IG)    + H_SIG_U(2)*U_DIM(2)/DY
                  A_MAT_I_JP1_K(I,J-1,K,IG) =  A_MAT_I_JP1_K(I,J-1,K,IG) + (1.-H_SIG_U(2))*U_DIM(2)/DY
      
                  A_MAT_I_J_K(I,J,K-1,IG)   = A_MAT_I_J_K(I,J,K-1,IG)    + H_SIG_U(3)*U_DIM(3)/DZ
                  A_MAT_I_J_KP1(I,J,K-1,IG) =  A_MAT_I_J_KP1(I,J,K-1,IG) + (1.-H_SIG_U(3))*U_DIM(3)/DZ

                  ! For the vectors...
                  IF(.NOT.UPWIND1ST) THEN
                     S_INC_LIM(I,J,K,IG) = S_INC_LIM(I,J,K,IG) &
                        - U_DIM(1)*(T_UPWIND(1)-T_LIMIT(1))/DX &
                        - U_DIM(2)*(T_UPWIND(2)-T_LIMIT(2))/DY &
                        - U_DIM(3)*(T_UPWIND(3)-T_LIMIT(3))/DZ
 
                     S_INC_LIM(I-1,J,K,IG)   = S_INC_LIM(I-1,J,K,IG)  + U_DIM(1)*(T_UPWIND(1)-T_LIMIT(1))/DX

                     S_INC_LIM(I,J-1,K,IG)   = S_INC_LIM(I,J-1,K,IG)  + U_DIM(2)*(T_UPWIND(2)-T_LIMIT(2))/DY
         
                     S_INC_LIM(I,J,K-1,IG)   = S_INC_LIM(I,J,K-1,IG)  + U_DIM(3)*(T_UPWIND(3)-T_LIMIT(3))/DZ
                  ENDIF
               END DO
            END DO 
         END DO 
      END DO 
      S_INC_LIM=S_INC_LIM + S
   ELSE ! ENDOF IF(NDIM_VEL.NE.0) THEN
      ! print *,'here1.1'
      ! print *,'NCOL_SIGMA_S,NX,NY,NZ,NG,NG2,NDIM_VEL,NITS,NITS_SOLV_NG:', &
      !          NCOL_SIGMA_S,NX,NY,NZ,NG,NG2,NDIM_VEL,NITS,NITS_SOLV_NG
      ! print *,'s=',s
      ! print *,'here 1.2'
      S_INC_LIM=S
   ENDIF ! IF(NDIM_VEL.NE.0) THEN ... ELSE...
   ! print *,'here2'

! Forward backward Gauss - Seidel...
   T_OLD = T_INITIAL
   DO ITIME=1,NTIME
      ! T_NEW = T_OLD
      DO ITS_NG=1,NITS_SOLV_NG
         R_T_MAX_DIF_NG = 0.0

         DO IIGG=0,ISWITCH_FBGS ! Switch on FBGS
            IG_START =1*(1-IIGG)  + (NG-1)*IIGG
            IG_FINISH=NG*(1-IIGG) + 2*IIGG
            IG_STEP  =1*(1-IIGG)  -  1*IIGG
            ILAST = ( 1-MIN(1, ABS(ITS_NG-NITS_SOLV_NG) ) )*IIGG ! The last backwards sweep of the last group iteration.
            IG_FINISH = IG_FINISH*(1-ILAST) + 1*ILAST
            DO IG=IG_START, IG_FINISH, IG_STEP
               RHS_SCAT=0.0 
               IF(NCOL_SIGMA_S.NE.0) THEN ! Unstructured spares storage
                  DO COUNT = FIN_SIGMA_S(IG), FIN_SIGMA_S(IG+1)-1
                     JG=COL_SIGMA_S(COUNT)
                     RHS_SCAT(:,:,:) =  RHS_SCAT(:,:,:) + SIGMA_S_OFF_PTER(:,:,:,COUNT) * T_NEW(:,:,:,JG)
                  END DO
               ELSE ! Structured dense storage.
                  DO JG=1,NG
                     RHS_SCAT(:,:,:) =  RHS_SCAT(:,:,:) + SIGMA_S_OFF(:,:,:,IG,JG) * T_NEW(:,:,:,JG)
                  END DO
               ENDIF
               T_NEW_TEMP_NG(:,:,:)=T_NEW(:,:,:,IG)

               DO ITS=1,NITS
                  R_T_MAX_DIF = 0.0

                  DO KK=1,0,-1
                     DO K = 2*KK + (NZ-1)*(1-KK), 2*(KK-1) + (NZ-1)*KK , 1*KK - 1*(KK-1)
                        DO JJ=1,0,-1
                           DO J = 2*JJ + (NY-1)*(1-JJ), 2*(JJ-1) + (NY-1)*JJ , 1*JJ - 1*(JJ-1)
                              DO II=1,0,-1
                                 DO I = 2*II + (NX-1)*(1-II), 2*(II-1) + (NX-1)*II , 1*II - 1*(II-1)
                                    A_I_J_K = A_MAT_I_J_K(I,J,K,IG)
                                    A_IP1_J_K =  A_MAT_IP1_J_K(I,J,K,IG)
                                    A_IM1_J_K =  A_MAT_IM1_J_K(I,J,K,IG)
                                    A_I_JP1_K =  A_MAT_I_JP1_K(I,J,K,IG)
                                    A_I_JM1_K =  A_MAT_I_JM1_K(I,J,K,IG)

                                    A_I_J_KP1 =  A_MAT_I_J_KP1(I,J,K,IG)
                                    A_I_J_KM1=A_MAT_I_J_KM1(I,J,K,IG)

                                    R_T_NEW =                                                                                    &
                                       (-A_IP1_J_K*T_NEW(I+1,J,K,IG) - A_IM1_J_K*T_NEW(I-1,J,K,IG) - A_I_JP1_K*T_NEW(I,J+1,K,IG) &
                                       - A_I_JM1_K*T_NEW(I,J-1,K,IG) - A_I_J_KP1*T_NEW(I,J,K+1,IG) - A_I_J_KM1*T_NEW(I,J,K-1,IG) &
                                       + S_INC_LIM(I,J,K,IG) + RHS_SCAT(I,J,K) + (1./(DT*V_N(IG)))*T_OLD(I,J,K,IG) ) / A_I_J_K
                                    R_T_MAX_DIF = MAX(ABS(R_T_NEW-T_NEW(I,J,K,IG) ), R_T_MAX_DIF)
                                    ! Relax the soln.
                                    T_NEW(I,J,K,IG) = RELAX*R_T_NEW + (1.-RELAX) * T_NEW(I,J,K,IG)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
                  ! print *,'its,nits,R_T_MAX_DIF:',its,nits,R_T_MAX_DIF 

                  IF(R_T_MAX_DIF.LT.ERROR_SOLV) THEN
                     ! print *,'its,nits,R_T_MAX_DIF,IG:',its,nits,R_T_MAX_DIF,IG
                     EXIT
                  ENDIF
                  ! T_TOTAL(:,:,:,:,ITS)=T_NEW(:,:,:,:)
               END DO ! DO ITS=1,NITS
               ! print *,'IG,ITS,R_T_MAX_DIF,ERROR_SOLV:',IG,ITS,R_T_MAX_DIF,ERROR_SOLV

               T_NEW(:,:,:,IG) = RELAX_NG*T_NEW(:,:,:,IG) + (1.-RELAX_NG) * T_NEW_TEMP_NG(:,:,:)
               R_T_MAX_DIF_NG = MAX( R_T_MAX_DIF_NG,  MAXVAL(ABS(T_NEW_TEMP_NG(:,:,:)-T_NEW(:,:,:,IG) ))  )
               ! print *,'ig,r_t_max_dif_ng:',ig,r_t_max_dif_ng
            END DO ! DO IG=IG_START, IG_FINISH, IG_STEP
         END DO ! DO IIGG=0,ISWITCH_FBGS ! Switch on FBGS
         ! T_TOTAL(:,:,:,:,ITS_NG)=T_NEW(:,:,:,:)
         ! PRINT *,'ITS_NG,R_T_MAX_DIF_NG,ERROR_SOLV_NG:',ITS_NG,R_T_MAX_DIF_NG,ERROR_SOLV_NG

         IF(R_T_MAX_DIF_NG.LT.ERROR_SOLV_NG) EXIT
      END DO ! ITS_G=1,NITS_SOLV_NG

      T_OLD=T_NEW ! Prepare for next time step.
      ! print *,'2.1 sum(t_new):',sum(t_new)
      ! PRINT *,' '
      ! print *,'inside solver S:'
      ! call PRINT_NORMAL(S,NX,NY,NZ,NG)
      ! PRINT *,' '
      ! print *,'inside solver T_INITIAL:'
      ! call PRINT_NORMAL(T_INITIAL,NX,NY,NZ,NG)
      ! PRINT *,' '
      ! print *,'inside solver T_NEW:'
      ! call PRINT_NORMAL(T_NEW,NX,NY,NZ,NG)
      ! PRINT *,'FINISHEING SOLVER ***'
      ! stop 2821

   END DO ! DO ITIME=1,NTIME
   ! print *,'here3'
END SUBROUTINE time_steping_diffusion_calc
