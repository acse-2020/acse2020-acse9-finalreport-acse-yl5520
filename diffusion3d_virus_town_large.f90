

          PROGRAM MAIN
	  IMPLICIT NONE
          REAL PIE, TOLER
          PARAMETER(PIE=3.141592654, TOLER=1.E-10) 
          integer, parameter :: sp = kind(1e0), dp = kind(1d0)
	  INTEGER NX,NY,NZ,NG,NONODS_NG,NG2, NLAY_IN,NLAY_OUT, NDIM_VEL
!          parameter(NX=12,NY=12,NZ=3,NG=8,NONODS_NG=(NX-2)*(NY-2)*(NZ-2)*NG,NLAY_IN=NONODS_NG,NLAY_OUT=NONODS_NG) 
!          parameter(NX=42,NY=42,NZ=3,NG=8,NONODS_NG=(NX-2)*(NY-2)*(NZ-2)*NG,NLAY_IN=NONODS_NG,NLAY_OUT=NONODS_NG) 
          parameter(NX=130,NY=130,NZ=3,NG=260,NG2=0,NONODS_NG=(NX-2)*(NY-2)*(NZ-2)*NG,NLAY_IN=NONODS_NG,NLAY_OUT=NONODS_NG) 
          PARAMETER(NDIM_VEL=2) ! if we have velocity then NDIM=1,2, or 3. 
! NITS_SOLV=no of non-linear iterations used for anything other than THETA=0.0 NITS_SOLV>1.
          REAL T_NEW(NX,NY,NZ,NG), &
	       T_OLD(NX,NY,NZ,NG), RHS(NX,NY,NZ,NG) !, A_MAT_T(NX,NY,NZ,NG)   
          REAL T_NON_LIN(NX,NY,NZ,NG)
          REAL T_NEW_LONG(NONODS_NG), T_DECOD_LONG(NONODS_NG),T_EXACT_LONG(NONODS_NG)
          REAL N1(NX,NY,NZ),N2(NX,NY,NZ)
          REAL DT,V_N(NG),DX,DY,DZ,RELAX_SOLV
          REAL SIGMAT(NX,NY,NZ,NG),KK(NX,NY,NZ,NG),S(NX,NY,NZ,NG), SIGMA_F(NX,NY,NZ,NG,NG2)
          REAL SIGMA_S_OFF(NX,NY,NZ,NG,NG2)
          REAL U(NDIM_VEL,NX+1,NY+1,NZ+1,NG) 
          LOGICAL SPAR_RECORD_S(NG,NG), SPAR_RECORD_F(NG,NG)
! SIGMAT is diagonal sigma_t, sigma_f is total fission matrix, sigma_s_off is just offdiagonal scattering. 
!          REAL A_MAT(NX,NY,NZ,NG,-1:1, -1:1, -1:1), B_MAT(NX,NY,NZ,NG,-1:1, -1:1,-1:1),&
!               FIS_MAT(NX,NY,NZ,NG,-1:1, -1:1, -1:1)
          REAL V(NX,NY,NZ,NG) ! Temp vector used for normalizing fission as well as calculating fission source
          INTEGER MATERIAL(NX,NY,NZ) ! Material number
! ANN:
       REAL, ALLOCATABLE :: WEIGHT(:),NEURAL(:,:),IDEAL_INPUT(:,:),IDEAL_OUTPUT(:,:),IDEAL_MIDLE(:,:)
       REAL, ALLOCATABLE :: CT(:,:),ALPHA_NEW(:),ALPHA(:),NEURAL_ONE(:)
       REAL, ALLOCATABLE :: INTEGRAL_T_NEW(:,:), SAMPLE_T_NEW(:,:) 
       REAL, ALLOCATABLE :: SIGMA_S_OFF_PTER(:,:,:,:), SIGMA_F_PTER(:,:,:,:)
       INTEGER, ALLOCATABLE :: FIN_SIGMA_S(:), COL_SIGMA_S(:), FIN_SIGMA_F(:), COL_SIGMA_F(:)
       INTEGER, ALLOCATABLE :: rd_mark_x(:),rd_mark_y(:)
       INTEGER :: ISAM_PT, JSAM_PT, KSAM_PT
       REAL :: W_EXP(1000)
       INTEGER :: NLAYER, NLAY(1000)
       INTEGER :: NCOLM,NEXAMPLE
! Local...
       INTEGER :: ILAY,IEXAMPLE,I,j,k,NOD,NITS_SOLV,NITS_SOLV_NG,&
		  NTIME,NTIME_ONE, NOREXP,NBATCH,NALPHA,IALPHA,ITEST,NTEST,ITS_ROM,NITS_ROM,NO_NUR
       INTEGER :: NCOLM2, NLAYER2, NLAY_IN2,NLAY_OUT2,NO_NUR2, NDEAL_MIDLE
       INTEGER :: II,JJ, ISTART,IFINI, JSTART,JFINI, IWIDTH,JWIDTH
       LOGICAL :: GOT_WEIGHT, FOUR_REGS
       REAL :: ERROR_TOLER,ERROR_SOLV_NG,NEU_ALPHA
       REAL :: SIG_RAN_NO,K_RAN_NO, DIAG_REG_TOLER(7)
       INTEGER :: SCALE_IN,SCALE_OUT, IRAN_NO, ISTART_NEU, NITS_EIG,ITS_EIG, NITS, NITS_GLOBAL, IG,JG
       INTEGER :: ITIME, ITS, OPTION
       REAL :: RAN_NO(100000)
       logical :: RECALL, EIGEN_START, pickup_WEIGHT
          LOGICAL RECORD_S, RECORD_F
       REAL :: MAX_SCALE_IN, MIN_SCALE_IN, MAX_SCALE_OUT, MIN_SCALE_OUT, KEFF, ERROR_SOLV,&
               ERROR_EIG,KEFF_EXACT
       real :: keff2, lambda, mean_keff_error, mean_t_error
       REAL :: ONE_DAY, ACCTIM, TIME_MAT
       REAL :: FORCE, H2M
       REAL :: DURATION_INFECT, INCUBATION_DURATION, VIRUS_SIGMA, VIRUS_GAMMA
       REAL :: XI, BETA, BETA2, MU, NU
       real :: VIRUS_N1_AIM, VIRUS_N2_AIM, length, LAMBDA_H_H, LAMBDA_M_M, R_EIGEN, RDAY
       real :: VIRUS_N1, VIRUS_N2
       real :: r_ratio, rswitch, THETA_NON_LIN
       LOGICAL :: EIGEN_VALUE, EIGEN_METHOD1
       INTEGER :: NSAVE_TO_FILE,ISAVE_TO_FILE, idivid, jdivid, il,jl, ilkeep, jlkeep, idisp, jdisp
       INTEGER ::  HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED
       INTEGER :: PARK, HOSPITAL, SCHOOL, OFFICES, SHOPS
       INTEGER :: BUILDINGS
       INTEGER :: NCOL_SIGMA_S,NCOL_SIGMA_F, COUNT_S, COUNT_F, I_UPWIND,I_HARMONIC, nits_solv_ng2
       INTEGER :: IFILE_COUNT
       REAL :: R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK
       REAL :: SIGN_SMALL ! function
       LOGICAL :: FOUND_FILE
       CHARACTER CHARA*240
       CHARACTER str_trim*240
! Local variables...
!
        I_UPWIND=1 ! Appy upwind differencing for the advection
        I_HARMONIC=0 ! Apply harmonic average for the diffusion calculation.
! roads...
         allocate(rd_mark_x(nx),rd_mark_y(ny) )
! Set up problem ********************************
         EIGEN_VALUE=.false.
!         EIGEN_VALUE=.true.
         EIGEN_METHOD1=.true. ! (new method=9.1437710515488124) (old method1=9.1435851016370382)
!         NSAVE_TO_FILE=1 ! Save the 1st file and then every N'the file unless =0 and then dont save.
         NSAVE_TO_FILE=4 ! Save the 1st file and then every N'the file unless =0 and then dont save.
         ISAVE_TO_FILE=1 ! for counting the number of time steps before outputing to disc
         IFILE_COUNT=1 ! Starting file number to generate. 
! original method1=0.99937582280063098; new method=0.99937582280063209
! ratio of indoor/outdoor = 1462/57 = 25.65
!         DT=1.e+3 ! time in seconds.
!         DT=0.1e+3 ! time in seconds.
         DT=1.0e+3 ! time in seconds.
         IF(EIGEN_VALUE) DT=1.e+10
         V_N=1.0 ! velocity
         RELAX_SOLV=1.0 ! RELAX_SOLV=1 (NOD RELAXATION), =0.5 TAKE 0.5 OF PREVIOUS VALUE TO CONVERGE BETTER IF PROBLEMS CONVERGING. 
         ONE_DAY=24.*3600.0 ! A DAY IN seconds. 

         DX = 4.0e+3/REAL(NX-2) ! domain 4km
         DY = 4.0e+3/REAL(NY-2)
         DZ = 1.0/REAL(NZ)
         NTIME_ONE=1

         KEFF=1.0

         if(EIGEN_VALUE) THEN
            THETA_NON_LIN=1.0 ! How to centre the non-linear terms (=1 fullly implicit, =0 fully explicity)
            T_NEW=1.0
            T_OLD=1.0
            EIGEN_START=.FALSE.
            NTIME=1
            NITS_SOLV=4
            NITS_SOLV_NG=100
            NITS_EIG=100

            NITS_SOLV=2 ! 5 ! 1
            NITS_SOLV_NG=20
            NITS_EIG=5

!            NITS_GLOBAL=21
            NITS_GLOBAL=42
!            ERROR_SOLV=1.0E-10 ! original
!            ERROR_SOLV_NG=1.0E-10 ! original 
            ERROR_SOLV=1.0E-6
            ERROR_SOLV_NG=1.0E-6
!            ERROR_SOLV_NG=1.0E-3
         ELSE ! Time dep...
             THETA_NON_LIN=0.75 ! How to centre the non-linear terms (=1 fullly implicit, =0 fully explicity)
             T_NEW=0.0
             T_OLD=0.0
!            NTIME=1000
!            NTIME=2000
!            NTIME=23000
!            NTIME= 2.25 * one_day/dt
!            NTIME= 45.75 * one_day/dt
!            NTIME= 45.5 * one_day/dt
!            NTIME= 2.0 * one_day/dt
!            NTIME= 1.75 * one_day/dt ! 15 hours to run
!            NTIME= 0.25 * one_day/dt
!            NTIME= 0.5 * one_day/dt
!            NTIME= 0.15 * one_day/dt
            NTIME= 40.5 * one_day/dt
!            NTIME= 1.0 * one_day/dt
!            NTIME= 0.1 * one_day/dt
!            NTIME= 1
!            NTIME= 45.0 * one_day/dt
!            NTIME= 9.0 * one_day/dt
!            NTIME= 10.25 * one_day/dt
!            ntime=1
!            NTIME=20
            NITS_SOLV=20!20!20
!            NITS_SOLV_NG=10
!            NITS_SOLV=400
!            NITS_SOLV=10000
!            NITS_SOLV_NG=4
!            NITS_SOLV_NG=4
            NITS_SOLV_NG=200!120
!            NITS_SOLV_NG=100
!            ERROR_SOLV=1.0E-6!4
            ERROR_SOLV=1.0E-7!7!4
!            ERROR_SOLV_NG=1.0E-4!4
            ERROR_SOLV_NG=1.0E-5!4
!            ERROR_SOLV_NG=1.0E-7!4
            NITS_EIG=0
            EIGEN_START=.FALSE.
!            NITS_GLOBAL=10 ! no of non-linear iterations per time step.
            NITS_GLOBAL=2!1!2 ! no of non-linear iterations per time step.
         ENDIF
         ALLOCATE( INTEGRAL_T_NEW(NG,NTIME),SAMPLE_T_NEW(NG,NTIME)  ) 
         S(:,:,:,:) = 0.0 
         ERROR_EIG=1.0E-3


! **************************R0
         INQUIRE(FILE='r0-10values.csv', EXIST=FOUND_FILE) 
         IF(FOUND_FILE) THEN
            open(27, file='r0-10values.csv')
!            WRITE(27, *) '# R0_HOME,     R0_OFFICE, R0_SCHOOL'
!            WRITE(27, *) '# R0_HOSPITAL, R0_PEDEST, R0_DIFF1'
!            WRITE(27, *) '# R0_DIFF2,    R0_SHOPS,  R0_VEHC_1TO5'
!            WRITE(27, *) '# R0_PARK'
            READ(27, *) 
!            READ(27, *) 
!            READ(27, *) 
!            READ(27, *) 
            READ(27, * ) R0_HOME,     R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2,    &
                         R0_SHOPS,  R0_VEHC_1TO5,  R0_PARK
         ELSE
            R0_HOME=0.2 !5.0 !0.2 ! No of people an infected person infects.
            R0_OFFICE=5.0 !5.0 !0.2 ! No of people an infected person infects.
            R0_SCHOOL=2.0 !5.0 !0.2 ! No of people an infected person infects.
            R0_HOSPITAL=5.0
            R0_PEDEST=1.5 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
            R0_DIFF1 =2.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
            R0_DIFF2 =2.0 !Portal R0 to hospital
            R0_PARK=1.0 
            R0_SHOPS =5.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
            R0_VEHC_1TO5 =2.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
            stop 2277
         ENDIF
         PRINT *,'R0_HOME,     R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2,'  &
                , '    R0_SHOPS,  R0_VEHC_1TO5,  R0_PARK:', & 
                  R0_HOME,     R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2,    &
                  R0_SHOPS,  R0_VEHC_1TO5,  R0_PARK
! **************************R0

! Define the materials...
! define homes-roads etc*********************************************************************
         HOMES=1     
         ROADS=2
         CROSS_ROADS=3   
         HOMES_OCCUPIED=4 
         PARK=5
         HOSPITAL=6
         SCHOOL=7
         OFFICES=8
         SHOPS=9
! define what MATERIAL is assign to HOMES, ROADS etc. 
         CALL DEFINE_GROUPS_MATERIALS(HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
               HOSPITAL, SCHOOL, OFFICES, SHOPS, &
               NX,NY,NZ, NG, MATERIAL, rd_mark_x, rd_mark_y) 
! FIND THE SPARCITY OF SIGMA_S, SIGMA_F that is find SPAR_RECORD_S, SPAR_RECORD_F

          NCOL_SIGMA_S=0 
       PRINT *,'ALLOCATEING SIGMA_S_OFF_PTER:'
          ALLOCATE(SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S) ) 
          ALLOCATE(FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S) )
          NCOL_SIGMA_F=0 
          ALLOCATE(SIGMA_F_PTER(NX,NY,NZ,NCOL_SIGMA_F) )
          ALLOCATE(FIN_SIGMA_F(NG+1), COL_SIGMA_F(NCOL_SIGMA_F) )

         IF(NG2==0) THEN ! Unstructured mesh sparcity...
            RECORD_S=.true.
            RECORD_F=.true.
            RDAY=0.0
!            ACCTIM=0.0
!            TIME_MAT=0.0 
!         ACCTIM=0.0
            ACCTIM=1.*ONE_DAY/2.0
            TIME_MAT=1.*ONE_DAY/2.0

            CALL DETERMINE_X_SECTIONS( KK,  SIGMAT,  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, &
                    SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, DX,DY,DZ,DT, &
                    OPTION,ITIME,ITS,ACCTIM,TIME_MAT, RDAY,EIGEN_VALUE,EIGEN_METHOD1, R_EIGEN, &
                    HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NEW, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
                    RECORD_S, SPAR_RECORD_S, RECORD_F, SPAR_RECORD_F, &
                    R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK) 
            DEALLOCATE(COL_SIGMA_S, COL_SIGMA_F ) 
       PRINT *,'DEALLOCATEING SIGMA_S_OFF_PTER:'
            DEALLOCATE(SIGMA_S_OFF_PTER, SIGMA_F_PTER )
            COUNT_S=0
            COUNT_F=0
            DO IG=1,NG 
               DO JG=1,NG 
                  IF(SPAR_RECORD_S(IG,JG)) COUNT_S=COUNT_S+1
                  IF(SPAR_RECORD_F(IG,JG)) COUNT_F=COUNT_F+1
               END DO
            END DO
            NCOL_SIGMA_S=COUNT_S
            NCOL_SIGMA_F=COUNT_F
            ALLOCATE(COL_SIGMA_S(NCOL_SIGMA_S), COL_SIGMA_F(NCOL_SIGMA_F) ) 
            ALLOCATE(SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S), SIGMA_F_PTER(NX,NY,NZ,NCOL_SIGMA_F)) 

            COUNT_S=0
            COUNT_F=0
            DO IG=1,NG 
               FIN_SIGMA_S(IG)=COUNT_S+1
               FIN_SIGMA_F(IG)=COUNT_F+1
               DO JG=1,NG 
                  IF(SPAR_RECORD_S(IG,JG)) THEN
                     COUNT_S=COUNT_S+1
                     COL_SIGMA_S(COUNT_S)=JG
                  ENDIF
                  IF(SPAR_RECORD_F(IG,JG)) THEN
                     COUNT_F=COUNT_F+1
                     COL_SIGMA_F(COUNT_F)=JG
                  ENDIF
               END DO
            END DO
            FIN_SIGMA_S(NG+1)=COUNT_S+1
            FIN_SIGMA_F(NG+1)=COUNT_F+1
         ENDIF ! ENDOF iF(NG2==0) THEN
         RECORD_S=.false.
         RECORD_F=.false.
! define homes-roads etc*********************************************************************


!         ISAM_PT=ISTART; JSAM_PT=JSTART; KSAM_PT=2
         ISAM_PT=6*NX/14; JSAM_PT=8*NY/14; KSAM_PT=2

!         T_OLD=0.0
!         T_NEW=0.0

! time stepping...
         ACCTIM=0.0
!         ACCTIM=1.0*ONE_DAY/2.0
!         print *,'acctim:',acctim
         DO ITIME=1,NTIME
         DO ITS=1,NITS_GLOBAL
! Material at which the material is evaluated. 
            IF(EIGEN_VALUE) THEN
               TIME_MAT=1.0*ONE_DAY/2.0
               ACCTIM=1.0*ONE_DAY/2.0
               R_EIGEN=1.0
               T_OLD=T_NEW
            ELSE
               TIME_MAT=ACCTIM
               R_EIGEN=0.0
            ENDIF


          print *,'ITS,NITS_GLOBAL, itime,ntime,TIME_MAT/max(toler,one_day):', &
                   ITS,NITS_GLOBAL, itime,ntime,TIME_MAT/max(toler,one_day)
          WRITE(*, '(200E12.3)' ) (maxval(t_new(:,:,:,ig)),ig=1,ng)
!          print *,'ig,(maxval(t_new(:,:,:,ig)),ig=1,ng):',ig,(maxval(t_new(:,:,:,ig)),ig=1,ng)
          do ig=1,-ng
             print *,'ig,minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig)):',ig,minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig))
          end do
      if(.false.) then
            print *,'TIME_MAT,ONE_DAY, TIME_MAT/ONE_DAY, (TIME_MAT/ONE_DAY) * PIE * 2.0:', &
                     TIME_MAT,ONE_DAY, TIME_MAT/ONE_DAY, (TIME_MAT/ONE_DAY) * PIE * 2.0
            print *,'SIN( (TIME_MAT/ONE_DAY) * PIE * 2.0 ):',SIN( (TIME_MAT/ONE_DAY) * PIE * 2.0 )
      endif
!            LAMBDA_H_H= 0.0
!            beta=0.0
!            beta2=0.0
!            LAMBDA_M_M= 1.0/one_day
            RDAY= 0.5*SIN( (TIME_MAT/ONE_DAY) * PIE * 2.0  - 0.5*PIE) + 0.5
            OPTION=1

! Determine KK,  SIGMAT,  SIGMA_S_OFF,  SIGMA_F...
!            T_NON_LIN=0.5*(T_NEW+t_old)
!            T_NON_LIN=T_OLD
            T_NON_LIN = THETA_NON_LIN * T_NEW  +  (1.0-THETA_NON_LIN) * T_OLD  
!         print *,'1 sum(t_new):',sum(t_new)

            CALL DETERMINE_X_SECTIONS( KK,  SIGMAT,  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, &
                    SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, DX,DY,DZ,DT, &
                    OPTION,ITIME,ITS,ACCTIM,TIME_MAT, RDAY,EIGEN_VALUE,EIGEN_METHOD1, R_EIGEN, &
                    HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NON_LIN, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
!                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_OLD, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
!                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NEW, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
                    RECORD_S, SPAR_RECORD_S, RECORD_F, SPAR_RECORD_F, &
                    R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK) 
!         print *,'2 sum(t_new):',sum(t_new) 

! change the sign as we have +ve on the rhs: 
            SIGMA_S_OFF      = - SIGMA_S_OFF
            SIGMA_S_OFF_PTER = - SIGMA_S_OFF_PTER
        if(.false.) then
            print *,'minval(SIGMAT),maxval(SIGMAT):',minval(SIGMAT),maxval(SIGMAT)
            print *,'minval(SIGMA_S_OFF),maxval(SIGMA_S_OFF):',minval(SIGMA_S_OFF),maxval(SIGMA_S_OFF)
            print *,'minval(KK),maxval(KK):',minval(KK),maxval(KK)
          do ig=1,-ng
          do jg=1,ng
             print *,'ig,jg:',ig,jg
            print *,'minval(SIGMA_S_OFF(:,:,:,ig,jg)),maxval(SIGMA_S_OFF(:,:,:,ig,jg)):', &
                     minval(SIGMA_S_OFF(:,:,:,ig,jg)),maxval(SIGMA_S_OFF(:,:,:,ig,jg))
          end do
          end do
            print *,'minval(SIGMA_F),maxval(SIGMA_F):',minval(SIGMA_F),maxval(SIGMA_F)
        endif ! if(.false.) then

!            CALL GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( A_MAT, B_MAT,  &
!                         DX,DY,DZ,DT,V_N,SIGMAT,KK,U,NDIM_VEL, NX,NY,NZ,NG) 

!      print *,'itime,its=',itime,its
!      print *,'minval(t_old),maxval(t_old):',minval(t_old),maxval(t_old)
!      print *,'minval(t_new),maxval(t_new):',minval(t_new),maxval(t_new) 
!         print *,'3 sum(t_new):',sum(t_new) 
            nits_solv_ng2=NITS_SOLV_NG/5
            if(ITS==NITS_GLOBAL) nits_solv_ng2=NITS_SOLV_NG

            CALL EIG_OR_TIME_STEPING_DIFFUSION_CALC(T_NEW,KEFF, T_OLD,&
!                 NTIME_ONE,NITS_SOLV,NITS_SOLV_NG,NITS_EIG, ERROR_SOLV,ERROR_SOLV_NG,&
                 NTIME_ONE,NITS_SOLV,NITS_SOLV_NG2,NITS_EIG, ERROR_SOLV,ERROR_SOLV_NG,&
                 ERROR_EIG,EIGEN_START, RELAX_SOLV, DX,DY,DZ,DT,V_N,SIGMAT, &
                 SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S,COL_SIGMA_S, KK, SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F,COL_SIGMA_F, S,U, &
                 I_UPWIND,I_HARMONIC, NDIM_VEL, NCOL_SIGMA_S,NCOL_SIGMA_F,NX,NY,NZ,NG,NG2 )
 
!         print *,'4 sum(t_new):',sum(t_new) 
!      print *,' '
!      print *,'minval(t_old),maxval(t_old):',minval(t_old),maxval(t_old)
!      print *,' '
!      print *,'minval(t_new),maxval(t_new):',minval(t_new),maxval(t_new)
!      print *,' '

            EIGEN_START=.FALSE.

         END DO ! DO ITS=1,NITS_GLOBAL
            ACCTIM=ACCTIM+DT
            DO IG=1,NG
               INTEGRAL_T_NEW(IG,ITIME) = SUM( T_NEW(2:NX-1, 2:NY-1, 2:NZ-1, IG) ) 
               SAMPLE_T_NEW(IG,ITIME) = T_NEW( ISAM_PT, JSAM_PT, KSAM_PT, IG)
            END DO
            T_OLD=T_NEW ! Prepare for next time step.

! ************* save to disc***************
         ISAVE_TO_FILE=ISAVE_TO_FILE+1

      IF((ISAVE_TO_FILE==NSAVE_TO_FILE) .or. ((NSAVE_TO_FILE>0).and. (ITIME==1))) THEN
         ISAVE_TO_FILE=0
!         IF(ITIME==1) THEN
!            open(27, file='group-output-time.csv', STATUS='new')
!            open(27, file='group-output-time-260g.csv', status='replace')
            str_trim=chara(IFILE_COUNT); k=INDEX(str_trim,' ')-1
            open(27, file='group-output-time-260g'//str_trim(1:k)//'.csv', status='replace')
            WRITE(27, *) '# DX,DY,DZ, NX,NY,NZ:', DX,DY,DZ, NX-2,NY-2,NZ-2
            WRITE(27, *) '# there are 8 groups: HOME-S, HOME-E, HOME-I, HOME-R, MOBILE-S, MOBILE-E, MOBILE-I, MOBILE-R'
            
            WRITE(27, *) '# DO K=1,NZ'
            WRITE(27, *) '# do j=1,NY,1,-1'
            WRITE(27, *) '#    WRITE(27, "(200E14.4)" ) (T_NEW( I, J, K, IG ) ,I=1,NX) ! in single procession for printing' 
            WRITE(27, *) '# end do'
            WRITE(27, *) '# end do'
!         ELSE
!            open(27, file='group-output-time.csv', Access = 'append')
!         ENDIF
         WRITE(27, *) '# TIME=', ACCTIM
         DO IG=1,NG
            WRITE(27, *) 'GROUP=',IG
            DO K=2,NZ-1
            do j=2,ny-1
               WRITE(27, '(200E14.4)' ) (SIGN_SMALL(T_NEW( I, J, K, IG )) ,I=2,NX-1) ! in single procession for printing
            end do
            end do
!           stop 382
         END DO
         close(27) 
         str_trim=chara(IFILE_COUNT); k=INDEX(str_trim,' ')-1
!         CALL SYSTEM('gzip run/group-output-time'//str_trim(1:k)//'.csv' ) ! factor of 12 reduction in disc space.
         CALL SYSTEM('rm group-output-time-260g'//str_trim(1:k)//'.csv.gz' ) ! factor of 12 reduction in disc space.
         CALL SYSTEM('gzip group-output-time-260g'//str_trim(1:k)//'.csv' ) ! factor of 12 reduction in disc space.
         IFILE_COUNT=IFILE_COUNT+1 ! next file.
       endif
! ************* save to disc***************

         END DO ! DO ITIME=1,NTIME

         CALL MAKE_LONG(T_NEW_LONG, T_NEW, NX,NY,NZ,NG, NONODS_NG) 

! SOLVE A PROBLEM USING ROM ********************

         IF(EIGEN_VALUE) PRINT *,'KEFF=',KEFF

! ******************************TEST THE INTERPOLATION*****************
         print *,'beta,beta2:',beta,beta2
         print *,'rday,1000.0 * (1.0-RDAY) + 1000.0:',rday,1000.0 * (1.0-RDAY) + 1000.0
         print *,'ntime, acctim, acctim/one_day, ntime*(one_day/max(toler,acctim)):', &
                  ntime, acctim, acctim/one_day, ntime*(one_day/max(toler,acctim))
      DO I=2,nx-1
      DO J=2,ny-1
      DO K=2,nz-1
         n1(i,j,k)=sum(t_new(i,j,k,1:4))
         n2(i,j,k)=sum(t_new(i,j,k,5:8))
      end do
      end do
      end do
      print *,'minval(n1),maxval(n1):',minval(n1),maxval(n1)
      print *,'minval(n2),maxval(n2):',minval(n2),maxval(n2)
      print *,'sum(n1),sum(n2),sum(n1)+sum(n2):',sum(n1),sum(n2),sum(n1)+sum(n2)
         DO IG=1,NG
!       T_NEW(:,:,:, IG)=  MATERIAL(:,:,:)
            PRINT *,'GROUP=',IG
      print *,'minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig)):',minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig))
      print *,'sum(t_new(:,:,:,ig)):',sum(t_new(:,:,:,ig))
!       print *,'t_new(:,:,:,ig):',t_new(:,:,:,ig)
       IF(NX<55) THEN ! Only print small problems to screen'
            DO K=2,NZ-1
            do j=ny-1,2,-1
!               PRINT *,(T_NEW( I, J, K, IG ),I=2,NX-1)
!               PRINT *,(real(T_NEW( I, J, K, IG ), sp),I=2,NX-1) ! in single procession for printing
!               do i=2,nx-1
               WRITE(*, '(200E10.2)' ) (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
!               print *, (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
!               end do
!               WRITE(*, '(F10.6)' ) ( T_NEW( I, J, K, IG ),I=2,nx-1) ! in single procession for printing
!               i=2
!               WRITE(*, '(E12.6)' ) T_NEW( I, J, K, IG ) ! in single procession for printing
            end do
            end do
!           stop 382
      ENDIF ! IF(NX<55) THEN
         END DO
! 
         IF(EIGEN_VALUE) PRINT *,'KEFF=',KEFF

! output to disc the group data
    if(.true.) then
         open(27, file='group-output-260g-end.csv')
            WRITE(27, *) '# DX,DY,DZ, NX,NY,NZ:', DX,DY,DZ, NX-2,NY-2,NZ-2
            WRITE(27, *) '# there are 8 groups: HOME-S, HOME-E, HOME-I, HOME-R, MOBILE-S, MOBILE-E, MOBILE-I, MOBILE-R'
            
            WRITE(27, *) '# DO K=1,NZ'
            WRITE(27, *) '# do j=1,NY,1,-1'
            WRITE(27, *) '#    WRITE(27, "(200E14.4)" ) (T_NEW( I, J, K, IG ) ,I=1,NX) ! in single procession for printing' 
            WRITE(27, *) '# end do'
            WRITE(27, *) '# end do'
         WRITE(27, *) '# TIME=', ACCTIM
         DO IG=1,NG
            WRITE(27, *) 'GROUP=',IG
            DO K=2,NZ-1
            do j=2,ny-1
!               PRINT *,(T_NEW( I, J, K, IG ),I=2,NX-1)
!               PRINT *,(real(T_NEW( I, J, K, IG ), sp),I=2,NX-1) ! in single procession for printing
!               do i=2,nx-1
!               WRITE(27, '(200E12.4)' ) (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
               WRITE(27, '(200E14.4)' ) (SIGN_SMALL(T_NEW( I, J, K, IG )) ,I=2,nx-1) ! in single procession for printing
!               print *, (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
!               end do
!               WRITE(*, '(F10.6)' ) ( T_NEW( I, J, K, IG ),I=2,nx-1) ! in single procession for printing
!               i=2
!               WRITE(*, '(E12.6)' ) T_NEW( I, J, K, IG ) ! in single procession for printing
            end do
            end do
!           stop 382
         END DO
         close(27) 
         CALL SYSTEM('rm group-output-260g-end.csv.gz' )
         CALL SYSTEM('gzip group-output-260g-end.csv' )

         open(27, file='group-time-integral')
         DO IG=1,NG
            ACCTIM=0.0
            DO ITIME=1,NTIME
               ACCTIM=ACCTIM+DT
               WRITE(27, * ) ACCTIM, INTEGRAL_T_NEW(IG,ITIME) 
            END DO
            WRITE(27, * ) 
         END DO
         close(27) 

         open(27, file='group-time-sample-pt')
         DO IG=1,NG
            ACCTIM=0.0
            DO ITIME=1,NTIME
               ACCTIM=ACCTIM+DT
               WRITE(27, * ) ACCTIM, SAMPLE_T_NEW(IG,ITIME) 
            END DO
            WRITE(27, * ) 
         END DO
         close(27) 

    endif
 

       STOP
       END PROGRAM MAIN
! 
! 
       REAL FUNCTION SIGN_SMALL(T_NEW)
       REAL T_NEW
       SIGN_SMALL=T_NEW

       if(abs(T_NEW)<1.e-25) SIGN_SMALL=0.0
       END 
! 
! 
! 
! python call looks like:
!      S_NEW,E_NEW,I_NEW,R_NEW = seir_4_eqn(S,E,I,R, beta, xi, sigma, gamma_c, dt, ntime) 
      subroutine seir_4_eqn(S_NEW,E_NEW,I_NEW,R_NEW, S,E,I,R, beta, xi, sigma, gamma_c, dt, ntime) 
! *******************************************************************************************
! This subroutine solve the classical 4 eqn SEIR model using an iterative explicit 
! method and as such the time step size can not be too small e.g. 10mins or DT=10.*60.=600.0
! dont forget there needs to be a seed of an infected person in the E group. 
! have a look at eqn 4 of virus_and_neutrons.pdf
! *******************************************************************************************
! S_NEW,E_NEW,I_NEW,R_NEW the calculated soln variables after ntime time steps output
! S,E,I,R initial soln variables sent down to this sub. 
! xi, sigma, gamma_c are the coefficients in the SEIR eqns sent down into this sub.
! dt=time step size
! ntime=no of time steps. 
     IMPLICIT NONE
     REAL, intent( out ) :: S_NEW,E_NEW,I_NEW,R_NEW
     REAL, intent( in ) :: S,E,I,R, beta, xi, sigma, gamma_c, dt
     INTEGER, intent( in ) :: ntime
! Local variables
      INTEGER ITS,NITS,itime ! no of non-linear iterations.
      REAL S_OLD,E_OLD,I_OLD,R_OLD, N ! previous time step values of S,E,I,R
!      REAL S_NEW,E_NEW,I_NEW,R_NEW
! 
      S_NEW=S; E_NEW=E; I_NEW=I; R_NEW=R
      nits=3
      do itime=1,ntime
         S_OLD=S_NEW; E_OLD=E_NEW; I_OLD=I_NEW; R_OLD=R_NEW
         DO ITS=1,NITS 
            N=S_NEW+E_NEW+I_NEW+R_NEW
            S_NEW=S_OLD + DT*( - BETA*S_NEW*I_NEW/N  +XI*R_NEW)
            E_NEW=E_OLD + DT*( + BETA*S_NEW*I_NEW/N  -SIGMA*E_NEW)
            I_NEW=I_OLD + DT*( SIGMA*E_NEW - GAMMA_C*I_NEW) 
            R_NEW=R_OLD + DT*( GAMMA_C*I_NEW  -XI*R_NEW)
         END DO ! DO ITS=1,NITS
      end do
      RETURN
      END subroutine seir_4_eqn


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


      CHARACTER*240 FUNCTION CHARA(CURPRO)
!     CURPRO must be the same in calling subroutine.
      INTEGER MXNLEN
      PARAMETER(MXNLEN=5)
      INTEGER DIGIT(MXNLEN),II,J,TENPOW,CURPRO,K
      CHARACTER STRI*240
      CHARACTER*1 CHADIG(MXNLEN)
      LOGICAL SWITCH
!     MXNLEN=max no of digits in integer I.
!     This function converts integer I to a string.
      II=CURPRO
      SWITCH=.FALSE.
      K=0
      STRI=' '
      do  J=MXNLEN,1,-1! Was loop 10
         TENPOW=10**(J-1)
         DIGIT(J)=II/TENPOW
         II=II-DIGIT(J)* TENPOW
!     ewrite(3,*)'DIGIT*J)=',DIGIT(J),' TENOW=',TENPOW
!     STRI=STRI//CHAR(DIGIT(J)+48)
         IF(DIGIT(J).EQ.0) CHADIG(J)='0'
         IF(DIGIT(J).EQ.1) CHADIG(J)='1'
         IF(DIGIT(J).EQ.2) CHADIG(J)='2'
         IF(DIGIT(J).EQ.3) CHADIG(J)='3'
         IF(DIGIT(J).EQ.4) CHADIG(J)='4'
         IF(DIGIT(J).EQ.5) CHADIG(J)='5'
         IF(DIGIT(J).EQ.6) CHADIG(J)='6'
         IF(DIGIT(J).EQ.7) CHADIG(J)='7'
         IF(DIGIT(J).EQ.8) CHADIG(J)='8'
         IF(DIGIT(J).EQ.9) CHADIG(J)='9'
!         Q=INDEX(STRI,' ')-1
!         IF(Q.LT.1) STRI=CHADIG
!         IF(Q.GE.1) STRI=STRI(1:Q)//CHADIG
!         IF(DIGIT(J).NE.0) SWITCH=.TRUE.
!         IF(SWITCH) K=K+1
      end do ! Was loop 10
       k=0
      do  J=MXNLEN,1,-1! Was loop 20
         IF(SWITCH) THEN 
           STRI=STRI(1:K)//CHADIG(J)
           K=K+1
         ENDIF
         IF((CHADIG(J).NE.'0').AND.(.NOT.SWITCH)) THEN 
           SWITCH=.TRUE.
           STRI=CHADIG(J)
           K=1
         ENDIF
      end do ! Was loop 20
       IF(K.EQ.0) THEN 
         STRI='0'
          K=1
        ENDIF
!       ewrite(3,*)'STRI=',STRI,' K=',K
!       ewrite(3,*)'STRI(MXNLEN-K+1:MXNLEN)=',STRI(MXNLEN-K+2:MXNLEN+1)
!       CHARA=STRI(MXNLEN-K+2:MXNLEN+1)
        CHARA=STRI(1:K)
       END
! 
! In pythin this is:
! ALPHA_INTERP = SIMPLE_INTERPOLATION(MATERIAL_INTERP,   MATERIAL,ALPHA_COEFFS, NEXAMPLE,NALPHA,NMATERIAL) 
       SUBROUTINE SIMPLE_INTERPOLATION(ALPHA_INTERP, MATERIAL_INTERP,   MATERIAL,ALPHA_COEFFS, NEXAMPLE,NALPHA,NMATERIAL) 
! This sub calculates the interpolated ALPHA_INTERP coeffs at a point MATERIAL_INTERP in material space. 
! It does this by interpolating NEXAMPLE points in material space at positions defined by MATERIAL and 
! at these positions we have the ALPHA_COEFFS which we are interpolating. 
! NMATERIAL are the number of different materials. 
! This subroutine uses an inverse distance weighting to perform the interpolation. 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NEXAMPLE, NALPHA,NMATERIAL
       REAL, intent( out ) :: ALPHA_INTERP(NALPHA)
       REAL, intent( in ) :: MATERIAL(NMATERIAL,NEXAMPLE),ALPHA_COEFFS(NALPHA,NEXAMPLE), MATERIAL_INTERP(NMATERIAL)
! local variables...
       REAL TOLER,INFINY
       PARAMETER(TOLER=1.0E-10,INFINY=1.0E+15)
       REAL, ALLOCATABLE :: D2(:),D2_SMALL(:),D_SMALL(:),W(:)
       INTEGER, ALLOCATABLE :: EXAMP_INDEX(:)
       INTEGER IEXAMPLE,IMAT, IEXAMP_SMALL, IALPHA
       REAL D2_SMALLEST

       ALLOCATE(D2(NEXAMPLE),D2_SMALL(NMATERIAL+1),D_SMALL(NMATERIAL+1),W(NMATERIAL+1))
       ALLOCATE(EXAMP_INDEX(NMATERIAL+1))
       
       DO IEXAMPLE=1,NEXAMPLE
          D2(IEXAMPLE)=SUM((MATERIAL_INTERP(:)-MATERIAL(:,IEXAMPLE))**2)
       END DO 
       
! find NMATERIAL+1 smallest values...       
       DO IMAT=1,NMATERIAL+1
          D2_SMALLEST=INFINY 
          DO IEXAMPLE=1,NEXAMPLE
             IF(D2(IEXAMPLE).LT.D2_SMALLEST) THEN ! insert into list
                D2_SMALLEST=D2(IEXAMPLE) 
                IEXAMP_SMALL = IEXAMPLE
             ENDIF
          END DO 
          D2_SMALL(IMAT)=D2_SMALLEST
          EXAMP_INDEX(IMAT)=IEXAMP_SMALL
          D2(IEXAMP_SMALL)=INFINY ! Set to large no so not to include this in list again. 
       END DO 

       DO IMAT=1,NMATERIAL+1
          D_SMALL(IMAT)=MAX(SQRT(D2_SMALL(IMAT)),TOLER)
       END DO
       W(:)=(1./D_SMALL(:))/ SUM(1./D_SMALL(:))

       DO IALPHA=1,NALPHA
          ALPHA_INTERP(IALPHA)= SUM( W(:)*ALPHA_COEFFS(IALPHA,EXAMP_INDEX(:)) )
       END DO
       
       END SUBROUTINE SIMPLE_INTERPOLATION
! 
! 
! 
! 
       SUBROUTINE KEEP_RANDOM_NO(SIG_RAN_NO,IRAN_NO,RAN_NO,RECALL)
! GENERATE AND STORE RANDOM NUMBERS. 
       INTEGER IRAN_NO
       REAL SIG_RAN_NO, RAN_NO(*)
       LOGICAL RECALL
! 
       IRAN_NO=IRAN_NO+1
       IF(RECALL) THEN
          SIG_RAN_NO = RAN_NO(IRAN_NO) 
       ELSE
          CALL RANDOM_NUMBER(SIG_RAN_NO) 
          RAN_NO(IRAN_NO) = SIG_RAN_NO
       ENDIF
       END SUBROUTINE KEEP_RANDOM_NO
! 
! 
! 
! 
          SUBROUTINE GET_CT_T_NEW(WEIGHT,NO_NUR,NCOLM,W_EXP,NLAY,NLAYER, &
              ALPHA_NEW, NALPHA, T_NEW_LONG, CT, NONODS_NG, SCALE_OUT,MAX_SCALE_OUT,MIN_SCALE_OUT ) 
! Form CT and T_NEW
! for ann:
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! This sub calculates the neuron values. 
! If W_EXP=0.0 find exponent of NEURONS otherwise dont. 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NO_NUR,NCOLM, NLAYER
       REAL, intent( in ) :: WEIGHT(NCOLM),W_EXP(NLAYER)
       INTEGER , intent( in ):: NLAY(NLAYER)
! for this subroutine:
          INTEGER, intent( in ) :: NALPHA,NONODS_NG, SCALE_OUT
          REAL, intent( in ) :: ALPHA_NEW(NALPHA),MAX_SCALE_OUT,MIN_SCALE_OUT
          REAL, intent( out ) :: T_NEW_LONG(NONODS_NG), CT(NALPHA,NONODS_NG)
             
! LOCAL VARIABLES...
          real EPSILON_PERT
          parameter(EPSILON_PERT=1.e-4)
          INTEGER NLAY_RECALL, ISTART_NEU,ISTART_WEIGHT, NO_NUR_RECALL, NCOLM_RECALL,IALPHA
          REAL, ALLOCATABLE :: NEURAL(:), NEURAL_TEMP(:)

          ALLOCATE(NEURAL(NO_NUR),NEURAL_TEMP(NO_NUR))

!          PRINT *,'NEXAMPLE,NBATCH:',NEXAMPLE,NBATCH
          NLAY_RECALL=(NLAYER-1)/2+1
          ISTART_NEU= SUM(NLAY( 1:(NLAYER-1)/2 ))+1
          NO_NUR_RECALL = NO_NUR - ISTART_NEU +1
          NCOLM_RECALL = NCOLM/2
          ISTART_WEIGHT=NCOLM/2 +1

          NEURAL(ISTART_NEU:ISTART_NEU+NALPHA-1)=ALPHA_NEW(:)
          CALL GETNEUVALS_FAST(NEURAL(ISTART_NEU),WEIGHT(ISTART_WEIGHT),&
               NO_NUR_RECALL,NCOLM_RECALL,W_EXP((NLAYER-1)/2 +1),NLAY((NLAYER-1)/2 +1),NLAY_RECALL)
          T_NEW_LONG(1:NONODS_NG) = NEURAL(NO_NUR-NONODS_NG+1:NO_NUR)
          IF(SCALE_OUT.NE.0) T_NEW_LONG = MIN_SCALE_OUT  + T_NEW_LONG*(MAX_SCALE_OUT - MIN_SCALE_OUT)
! unscale... 

          DO IALPHA=1,NALPHA  
             NEURAL_TEMP=NEURAL
             NEURAL_TEMP(ISTART_NEU-1+IALPHA) = NEURAL_TEMP(ISTART_NEU-1+IALPHA) + EPSILON_PERT
             CALL GETNEUVALS_FAST(NEURAL_TEMP(ISTART_NEU),WEIGHT(ISTART_WEIGHT),&
                  NO_NUR_RECALL,NCOLM_RECALL,W_EXP((NLAYER-1)/2 +1),&
                  NLAY((NLAYER-1)/2 +1),NLAY_RECALL)
             CT(IALPHA,:) =  (NEURAL_TEMP(NO_NUR-NONODS_NG+1:NO_NUR) - NEURAL(NO_NUR-NONODS_NG+1:NO_NUR))/EPSILON_PERT
             IF(SCALE_OUT.NE.0) CT(IALPHA,:) =CT(IALPHA,:) * (MAX_SCALE_OUT - MIN_SCALE_OUT)
          END DO
          
          END SUBROUTINE GET_CT_T_NEW
! 
! 
! 
! In python:  
! ALPHA_NEW = ONE_AE_ITERATION( ALPHA, A_MAT,SIGMA_S_OFF,B_MAT,C, U_NEW, U_OLD, S, DIAG_REG_TOLER, NX,NY,NZ,NG, NALPHA,NONODS_NG) 
         SUBROUTINE ONE_AE_ITERATION(ALPHA_NEW,  ALPHA, A_MAT,SIGMA_S_OFF,SIGMA_S_OFF_PTER, FIN_SIGMA_S, &
                                     COL_SIGMA_S, B_MAT,CT, U_NEW, U_OLD, S,&
                                      DIAG_REG_TOLER, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, NALPHA,NONODS_NG) 
	 IMPLICIT NONE
	 INTEGER, INTENT(IN) :: NX,NY,NZ,NG,NG2, NCOL_SIGMA_S, NALPHA,NONODS_NG
! Perform matrix vector multiplication.
! B=A*X
! DIAG_REG_TOLER contains regularization terms.
         REAL, INTENT(OUT) :: ALPHA_NEW(NALPHA)
         REAL, INTENT(IN) :: U_NEW(NX,NY,NZ,NG), U_OLD(NX,NY,NZ,NG), S(NX,NY,NZ,NG), ALPHA(NALPHA), DIAG_REG_TOLER(7) 
         REAL, INTENT(IN) :: CT(NALPHA,NONODS_NG), A_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1), B_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1)
         REAL, INTENT(IN) :: SIGMA_S_OFF(NX,NY,NZ,NG,NG2) 

          REAL, INTENT(IN) :: SIGMA_S_OFF_PTER(NCOL_SIGMA_S) 
          INTEGER, INTENT(IN) :: FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S)
!          REAL, INTENT(IN) :: SIGMA_F_PTER(NCOL_SIGMA_F) 
!          INTEGER, INTENT(IN) :: FIN_SIGMA_F(NG+1), COL_SIGMA_F(NCOL_SIGMA_F)
! Local variables...
         REAL, ALLOCATABLE :: CTAC(:,:), D(:,:,:,:), A_U(:,:,:,:),P(:),C_P(:,:,:,:),A_C_P(:,:,:,:),B_U_OLD(:,:,:,:), C(:,:) 
         REAL, ALLOCATABLE :: RHS(:), S_LONG(:),A_U_LONG(:),A_C_P_LONG(:), D_LONG(:), B_U_OLD_LONG(:)
         REAL, ALLOCATABLE :: CT_A_U(:), CT_B_U_OLD(:), CT_A_C_P(:), C_P_LONG(:), DALPHA(:)
         REAL, ALLOCATABLE :: SIGMA_S_ZERO(:,:,:,:,:), SIGMA_S_ZERO_OFF_PTER(:,:,:,:)
         INTEGER I,NOD,IALPHA

         allocate( CTAC(NALPHA,NALPHA),D(NX,NY,NZ,NG), A_U(NX,NY,NZ,NG),P(NALPHA),C_P(NX,NY,NZ,NG) ) 
         allocate( A_C_P(NX,NY,NZ,NG), B_U_OLD(NX,NY,NZ,NG),  RHS(NALPHA), C(NONODS_NG,NALPHA) ) 
         allocate( S_LONG(NONODS_NG),A_U_LONG(NONODS_NG),A_C_P_LONG(NONODS_NG), D_LONG(NONODS_NG), B_U_OLD_LONG(NONODS_NG)  ) 
         allocate( CT_A_U(NALPHA),CT_B_U_OLD(NALPHA),CT_A_C_P(NALPHA), C_P_LONG(NONODS_NG), DALPHA(NALPHA) )
         ALLOCATE( SIGMA_S_ZERO(NX,NY,NZ,NG,NG2), SIGMA_S_ZERO_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S) )
         SIGMA_S_ZERO=0.0
         SIGMA_S_ZERO_OFF_PTER=0.0

         DO NOD=1,NONODS_NG
            DO IALPHA=1,NALPHA
               C(NOD,IALPHA) = CT(IALPHA,NOD) 
            END DO
         END DO
         ! print *,'here1'
       
         CALL MATRIX_VEC_MULT_TIME_STEPING( A_U, A_MAT,SIGMA_S_OFF,SIGMA_S_OFF_PTER, &
                    FIN_SIGMA_S, COL_SIGMA_S,U_NEW, NCOL_SIGMA_S,NX,NY,NZ,NG,NG2) 
         CALL MAKE_LONG(A_U_LONG, A_U, NX,NY,NZ,NG,NONODS_NG) 
         CT_A_U=matmul(CT,A_U_LONG)

         CALL MATRIX_VEC_MULT_TIME_STEPING( B_U_OLD, B_MAT,SIGMA_S_ZERO,SIGMA_S_ZERO_OFF_PTER, &
                   FIN_SIGMA_S, COL_SIGMA_S,U_OLD, NCOL_SIGMA_S,NX,NY,NZ,NG,NG2) 
         CALL MAKE_LONG(B_U_OLD_LONG, B_U_OLD, NX,NY,NZ,NG,NONODS_NG) 
         CT_B_U_OLD=matmul(CT,B_U_OLD_LONG)

         CALL MAKE_LONG(S_LONG, S,NX,NY,NZ,NG,NONODS_NG) 

         RHS=-CT_A_U + CT_B_U_OLD + matmul(CT,S_LONG)

! The matrix FORMATION:
         do i=1,nalpha
            p=0.0
            p(i)=1.0

            C_P_LONG=matmul(C,P)
            CALL MAKE_GRID( C_P, C_P_LONG, NX,NY,NZ,NG,NONODS_NG)

            CALL MATRIX_VEC_MULT_TIME_STEPING( A_C_P, A_MAT,SIGMA_S_OFF,SIGMA_S_OFF_PTER, &
                             FIN_SIGMA_S, COL_SIGMA_S,C_P, NCOL_SIGMA_S,NX,NY,NZ,NG,NG2) 

            CALL MAKE_LONG(A_C_P_LONG, A_C_P, NX,NY,NZ,NG,NONODS_NG)
            CT_A_C_P=matmul(CT,A_C_P_LONG)     
            CTAC(:,I)=CT_A_C_P(:)                  
         end do

! Regularize CTAC and RHS vector before solution:
         CALL REGULARIZE_CTAC_RHS(DIAG_REG_TOLER, CTAC, ALPHA, NALPHA) 

! Solve CTAC alpha = rhs:      
!         CALL SMLINN(CTAC,ALPHA_NEW,RHS,NALPHA,NALPHA) 
         CALL SMLINN(CTAC,DALPHA,RHS,NALPHA,NALPHA) 
         ALPHA_NEW =  DALPHA + ALPHA 

         END SUBROUTINE ONE_AE_ITERATION
! 
! 
! 
! 
         subroutine REGULARIZE_CTAC_RHS(DIAG_REG_TOLER, CTAC, ALPHA, NALPHA) 
! Regularize CTAC and RHS vector before solution:
! ALPHA is the lattest soln in compressed space. 
         implicit none
         INTEGER NALPHA
         REAL CTAC(NALPHA,NALPHA), ALPHA(NALPHA), DIAG_REG_TOLER(7) 
! local variables...
         REAL, ALLOCATABLE :: DIAG(:) 
         INTEGER IALPHA

         allocate( DIAG(NALPHA) )
! 
         DO IALPHA=1,NALPHA
            DIAG(IALPHA)= &
         + DIAG_REG_TOLER(1)                             + MAXVAL(ABS(CTAC(:,:)))*DIAG_REG_TOLER(2) &
         + MAXVAL(ABS(CTAC(IALPHA,:)))*DIAG_REG_TOLER(3) + ABS(CTAC(IALPHA,IALPHA))*DIAG_REG_TOLER(4) &
         + SUM(ABS(CTAC(:,:)))*DIAG_REG_TOLER(5)         + SUM(ABS(CTAC(IALPHA,:)))*DIAG_REG_TOLER(6) 
         END DO

         DO IALPHA=1,NALPHA
            CTAC(IALPHA,IALPHA) = CTAC(IALPHA,IALPHA) + DIAG(IALPHA)
!            RHS(IALPHA) = RHS(IALPHA) + DIAG(IALPHA)*ALPHA(IALPHA) 
         END DO     
         end subroutine REGULARIZE_CTAC_RHS
! 
! 
! 
! 
         subroutine MAKE_LONG(b_long, B_GRID, NX,NY,NZ,NG,NONODS_NG) 
         implicit none
         integer NX,NY,NZ,NG,NONODS_NG
         real B_GRID(NX,NY,NZ,NG)
         real b_long(NONODS_NG)
! Local variables...
         integer i,j,k,ig, nx2,ny2,nz2
         nx2=nx-2; ny2=ny-2; nz2=nz-2
         do ig=1,ng
            do k=1,nz2
               do j=1,ny2
                  do i=1,nx2
                     b_long(i+(j-1)*nx2 + (k-1)*nx2*ny2 + (ig-1)*nx2*ny2*nz2) = b_grid(i+1,j+1,k+1,ig) 
                  end do
               end do
            end do
         end do
         
         end subroutine MAKE_LONG
! 
! 
! 
! 
         subroutine MAKE_GRID(b_grid, B_LONG, NX,NY,NZ,NG,NONODS_NG) 
         implicit none
         integer NX,NY,NZ,NG,NONODS_NG
         real B_LONG(NONODS_NG)
         real b_grid(NX,NY,NZ,NG) 
! Local variables...
         integer i,j,k,ig, nx2,ny2,nz2
         nx2=nx-2; ny2=ny-2; nz2=nz-2
         b_grid=0.0
         do ig=1,ng
            do k=1,nz2
               do j=1,ny2
                  do i=1,nx2
                     b_grid(i+1,j+1,k+1,ig) = b_long(i+(j-1)*nx2+ (k-1)*nx2*ny2 + (ig-1)*nx2*ny2*nz2) 
                  end do
               end do
            end do
         end do
         
         end subroutine MAKE_GRID
! 
! 
!     
!	  
        SUBROUTINE SMLINN(A,X,B,NMX,N)
        IMPLICIT NONE
        INTEGER NMX,N
        REAL A(NMX,NMX),X(NMX),B(NMX)
        REAL R
        INTEGER K,I,J
!     Form X = A^{-1} B
!     Useful subroutine for inverse
!     This sub overwrites the matrix A. 
        DO K=1,N-1
           DO I=K+1,N
              A(I,K)=A(I,K)/A(K,K)
           END DO
           DO J=K+1,N
              DO I=K+1,N
                 A(I,J)=A(I,J) - A(I,K)*A(K,J)
              END DO
           END DO
        END DO
!     
!     Solve L_1 x=b
        DO I=1,N
           R=0.
           DO J=1,I-1
              R=R+A(I,J)*X(J)
           END DO
           X(I)=B(I)-R
        END DO
!     
!     Solve U x=y
        DO I=N,1,-1
           R=0.
           DO J=I+1,N
              R=R+A(I,J)*X(J)
           END DO
           X(I)=(X(I)-R)/A(I,I)
        END DO
        RETURN
        END
!     
! 
! In phython use:
! b = MATRIX_VEC_MULT_TIME_STEPING( A,SIGMA,X, NX,NY,NZ,NG) 
         SUBROUTINE MATRIX_VEC_MULT_TIME_STEPING( B, A,SIGMA_S_OFF,SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, X, &
                                                  NCOL_SIGMA_S,NX,NY,NZ,NG,NG2) 
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: NCOL_SIGMA_S, NX,NY,NZ,NG,NG2
! Perform matrix vector multiplication.
! B=A*X
! NG2=0 if NCOL_SIGMA_S .ne.0 
          REAL, INTENT(OUT) :: B(NX,NY,NZ,NG)
          REAL, INTENT(IN) :: X(NX,NY,NZ,NG), A(NX,NY,NZ,NG, -1:1, -1:1, -1:1),SIGMA_S_OFF(NX,NY,NZ,NG,NG2)
          REAL, INTENT(IN) :: SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S)
          INTEGER, INTENT(IN) :: FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S)
! Local variables...
       INTEGER I,J,K,IG, II,JJ,KK, JG, COUNT

       B=0.0
 
       IF(NCOL_SIGMA_S.NE.0) THEN ! Unstructured spares storage
          DO IG=1,NG
             DO COUNT = FIN_SIGMA_S(IG), FIN_SIGMA_S(IG+1)-1
                JG=COL_SIGMA_S(COUNT)
                B(:,:,:,IG)=B(:,:,:,IG) - SIGMA_S_OFF_PTER(:,:,:,COUNT) * X(:,:,:,JG) 
             END DO
          END DO
       ELSE ! Structured dense storage.
          DO IG=1,NG
          DO JG=1,NG
             B(:,:,:,IG)=B(:,:,:,IG) - SIGMA_S_OFF(:,:,:,IG,JG)*X(:,:,:,JG) 
          END DO
          END DO
       ENDIF


       DO KK=-1,1
       DO JJ=-1,1
       DO II=-1,1

          DO IG = 1,NG
          DO K = 2,NZ-1
          DO J = 2,NY-1  
          DO I = 2,NX-1 

            B(I,J,K,IG) = B(I,J,K,IG) + A(I,J,K,IG,II,JJ,KK)*X(I+II,J+JJ,K+KK,IG) 
          END DO
          END DO 
          END DO 
          END DO 

       END DO
       END DO 
       END DO 
       END SUBROUTINE MATRIX_VEC_MULT_TIME_STEPING
! 
! 
! 
! 
       SUBROUTINE EIG_OR_TIME_STEPING_DIFFUSION_CALC(T_NEW,KEFF, T_INITIAL,&
                  NTIME,NITS_SOLV,NITS_SOLV_NG,NITS_EIG,ERROR_SOLV,ERROR_SOLV_NG,&
                  ERROR_EIG,EIGEN_START, RELAX_SOLV, DX,DY,DZ,DT,V_N,SIGMAT,&
                  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S,COL_SIGMA_S, K, SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F,COL_SIGMA_F, S,U, &
                  I_UPWIND,I_HARMONIC, NDIM_VEL, NCOL_SIGMA_S,NCOL_SIGMA_F,NX,NY,NZ,NG,NG2 )
! This subroutine solves the advection diffusion equation - could be used for neutronics diffusion or SN calculations or Virus modelling. 
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: I_UPWIND,I_HARMONIC, NDIM_VEL, NCOL_SIGMA_S, NCOL_SIGMA_F, NX,NY,NZ,NG,NG2
! NITS=no of non-linear iterations 
! NTIME=NO OF TIME STEPS
! DT=TIME STEP SIZE (CAN USE DT=1.E+15 for steady state cacls) 
! I_UPWIND = 1 Apply 1st order upwinding to advection, =0 apply high order flux limiting.
! I_HARMONIC = 1 perform a Harmonic average of the diffusion along the CV boundary, =0 apply mean diffusion coefficient at CV boundary. 
          REAL, INTENT(OUT) :: T_NEW(NX,NY,NZ,NG),KEFF
          REAL, INTENT(IN) :: T_INITIAL(NX,NY,NZ,NG) 
          REAL, INTENT(IN) :: RELAX_SOLV,DX,DY,DZ,DT,V_N(NG), ERROR_SOLV,ERROR_SOLV_NG,ERROR_EIG
          REAL, INTENT(IN) :: SIGMAT(NX,NY,NZ,NG),K(NX,NY,NZ,NG),SIGMA_F(NX,NY,NZ,NG,NG2),&
                              SIGMA_S_OFF(NX,NY,NZ,NG,NG2),S(NX,NY,NZ,NG)
          REAL, INTENT(IN) :: U(NDIM_VEL,NX+1,NY+1,NZ+1,NG)  
          REAL, INTENT(IN) :: SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S) 
          INTEGER, INTENT(IN) :: FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S)
          REAL, INTENT(IN) :: SIGMA_F_PTER(NX,NY,NZ,NCOL_SIGMA_F) 
          INTEGER, INTENT(IN) :: FIN_SIGMA_F(NG+1), COL_SIGMA_F(NCOL_SIGMA_F)
	  INTEGER, INTENT(IN) :: NTIME,NITS_SOLV,NITS_SOLV_NG,NITS_EIG
          LOGICAL, INTENT(IN) :: EIGEN_START
! Local variables...
       REAL, ALLOCATABLE :: V(:,:,:,:),A_MAT_T(:,:,:,:),A_MAT(:,:,:,:,:,:,:), S_INC_LIM(:,:,:,:)
       INTEGER ITIME,ITS,II,I,JJ,J,IG,JG,ITS_EIG, COUNT
       REAL LAMBDA,KEFF2

       IF(NITS_EIG==0) THEN
!         print *,'here 1 s=',s
          CALL TIME_STEPING_DIFFUSION_CALC(T_NEW, T_INITIAL, NTIME,NITS_SOLV,NITS_SOLV_NG,RELAX_SOLV,ERROR_SOLV,ERROR_SOLV_NG, &
                         DX,DY,DZ,DT,V_N,SIGMAT,SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S,COL_SIGMA_S,K,S,U,  &
                         I_UPWIND,I_HARMONIC, NDIM_VEL, NCOL_SIGMA_S,NX,NY,NZ,NG,NG2)
       ELSE
          ALLOCATE(V(NX,NY,NZ,NG),A_MAT_T(NX,NY,NZ,NG), A_MAT(NX,NY,NZ,NG,-1:1,-1:1, -1:1) ) 
          ALLOCATE(S_INC_LIM(NX,NY,NZ,NG))
          S_INC_LIM=0.0 ! advection limiting rhs may be not used.
          
!          T_NEW=T_INITIAL
          IF(EIGEN_START) T_NEW=1.0
          IF(EIGEN_START) KEFF=1.E+10
          CALL GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( A_MAT, S_INC_LIM, &
                         DX,DY,DZ,DT,V_N,SIGMAT,K,T_NEW,U, I_UPWIND,I_HARMONIC, NDIM_VEL, NX,NY,NZ,NG) 
          DO ITS_EIG=1,NITS_EIG
!             V=SIGMA_F*T_NEW ! Matrix vector
             V = S_INC_LIM ! include limiting
             IF(NCOL_SIGMA_S.NE.0) THEN ! Unstructured spares storage
                DO IG=1,NG
                   DO COUNT = FIN_SIGMA_F(IG), FIN_SIGMA_F(IG+1)-1
                      JG=COL_SIGMA_F(COUNT)
                      V(:,:,:,IG) =  V(:,:,:,IG) + SIGMA_F_PTER(:,:,:,COUNT) * T_NEW(:,:,:,JG) 
                   END DO
                END DO
             ELSE ! Structured dense storage.
                DO IG=1,NG
                   DO JG=1,NG
                      V(:,:,:,IG)=V(:,:,:,IG) + SIGMA_F(:,:,:,IG,JG)*T_NEW(:,:,:,JG) 
                   END DO
                END DO
             ENDIF
                  
             CALL TIME_STEPING_DIFFUSION_CALC(T_NEW, T_INITIAL, NTIME,NITS_SOLV,NITS_SOLV_NG,RELAX_SOLV,ERROR_SOLV,ERROR_SOLV_NG, &
                         DX,DY,DZ,DT,V_N,SIGMAT,SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S,COL_SIGMA_S,K,V,U,  &
                         I_UPWIND,I_HARMONIC, NDIM_VEL, NCOL_SIGMA_S,NX,NY,NZ,NG,NG2)
! Normalize...
             V=0.0
             IF(NCOL_SIGMA_S.NE.0) THEN ! Unstructured spares storage
                DO IG=1,NG
                   DO COUNT = FIN_SIGMA_F(IG), FIN_SIGMA_F(IG+1)-1
                      JG=COL_SIGMA_F(COUNT)
                      V(:,:,:,IG) =  V(:,:,:,IG) + SIGMA_F_PTER(:,:,:,COUNT) * T_NEW(:,:,:,JG) 
                   END DO
                END DO
             ELSE ! Structured dense storage.
                DO IG=1,NG
                   DO JG=1,NG
                      V(:,:,:,IG)=V(:,:,:,IG) + SIGMA_F(:,:,:,IG,JG)*T_NEW(:,:,:,JG) 
                   END DO
                END DO
             ENDIF
             T_NEW = T_NEW / SUM(V) 
!             T_NEW = T_NEW / SUM(SIGMA_F*T_NEW) 
! Eigen-value...
             CALL MATRIX_VEC_MULT_TIME_STEPING( A_MAT_T, A_MAT,SIGMA_S_OFF,SIGMA_S_OFF_PTER, &
                                 FIN_SIGMA_S, COL_SIGMA_S,T_NEW, NCOL_SIGMA_S,NX,NY,NZ,NG,NG2) 
             A_MAT_T = A_MAT_T - S_INC_LIM ! include limiting
             V=0.0
             IF(NCOL_SIGMA_S.NE.0) THEN ! Unstructured spares storage
                DO IG=1,NG
                   DO COUNT = FIN_SIGMA_F(IG), FIN_SIGMA_F(IG+1)-1
                      JG=COL_SIGMA_F(COUNT)
                      V(:,:,:,IG) =  V(:,:,:,IG) + SIGMA_F_PTER(:,:,:,COUNT) * T_NEW(:,:,:,JG) 
                   END DO
                END DO
             ELSE ! Structured dense storage.
                DO IG=1,NG
                   DO JG=1,NG
                      V(:,:,:,IG)=V(:,:,:,IG) + SIGMA_F(:,:,:,IG,JG)*T_NEW(:,:,:,JG) 
                   END DO
                END DO
             ENDIF
	     !PRINT *, SUM(A_MAT_T), SUM(V)
             LAMBDA = SUM(A_MAT_T) / SUM(V)
!             LAMBDA = SUM(A_MAT_T) / SUM(SIGMA_F*T_NEW)
             KEFF2=1.0/MAX(1.E-10, LAMBDA) 
!              print *,'its_eig,keff2=',its_eig,keff
              print *,'SUM(A_MAT_T),  SUM(V), sum(t_new):',SUM(A_MAT_T),  SUM(V), sum(t_new)
              print *,'its_eig,NITS_EIG,keff2,keff,ABS(KEFF2-KEFF),ERROR_EIG=', &
                       its_eig,NITS_EIG,keff2,keff,ABS(KEFF2-KEFF),ERROR_EIG
		!KEFF_ALL(ITS_EIG)=KEFF2
             IF(ABS(KEFF2-KEFF).LT.ERROR_EIG) THEN
                KEFF=KEFF2
                EXIT
             ELSE
                KEFF=KEFF2
             ENDIF
          END DO ! DO ITS_EIG=1,NITS_EIG
       ENDIF
       END SUBROUTINE 

! 
! 
! 
! In phython use:
! A_MAT, B_MAT, S_INC_LIM = GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( T_INITIAL, DX,DY,DZ,DT,V_N,SIGMAT,K,S, NX,NY,NZ,NG)
         SUBROUTINE GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( A_MAT, S_INC_LIM,  &
                         DX,DY,DZ,DT,V_N,SIGMAT,KDIFF,T_NEW,U, I_UPWIND,I_HARMONIC, NDIM_VEL, NX,NY,NZ,NG) 
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: I_UPWIND,I_HARMONIC, NDIM_VEL, NX,NY,NZ,NG
! NITS=no of non-linear iterations 
! NTIME=NO OF TIME STEPS
! DT=TIME STEP SIZE (CAN USE DT=1.E+15 for steady state cacls) 
          REAL, INTENT(OUT) :: A_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1)
!          REAL, INTENT(OUT) :: B_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1)
          REAL, INTENT(OUT) :: S_INC_LIM(NX,NY,NZ,NG)
          REAL, INTENT(IN) :: DX,DY,DZ,DT,V_N(NG)
          REAL, INTENT(IN) :: SIGMAT(NX,NY,NZ,NG),KDIFF(NX,NY,NZ,NG)
          REAL, INTENT(IN) :: T_NEW(NX,NY,NZ,NG), U(NDIM_VEL,NX+1,NY+1,NZ+1,NG)
! Local variables...
       REAL, PARAMETER :: TOLER=1.0E-10, XI_LIMIT = 2.0 ! defines TVD curve on the NVD diagram. 
       LOGICAL UPWIND1ST ! Apply 1st order upwinding to advection. 
       LOGICAL HARMONIC_AVE ! perform a Harmonic average of the diffusion along the CV boundary
       REAL :: DENOIN,CTILIN,DENOOU,CTILOU,FTILIN,FTILOU
       INTEGER I,J,K, IG
       REAL DIFF_I 
       REAL A_I_J_K, A_IP1_J_K, A_IM1_J_K, A_I_JP1_K, A_I_JM1_K, A_I_J_KP1, A_I_J_KM1

       INTEGER SWITCH_UP1(3), SWITCH_UP2(3), KM, JM, IM, IDIM, ISWITCH, JSWITCH, KSWITCH
       INTEGER I1,I2,I3, J1,J2,J3, K1,K2,K3
       INTEGER II1,II2,II3, JJ1,JJ2,JJ3, KK1,KK2,KK3
       REAL DENOIN1, DENOIN2, CTILIN1, CTILIN2, FTILIN1, FTILIN2, T_LIMIT1, T_LIMIT2, U_DIM(3)
       REAL RH_SIG_U, T_UPWIND(3), T_LIMIT(3), H_SIG_U(3) 
! function
       REAL FUN_HARMONIC

       UPWIND1ST = (I_UPWIND==1) 
       HARMONIC_AVE = (I_HARMONIC==1) 

       A_MAT=0.0
!       B_MAT=0.0

	DO IG = 1,NG
	DO K = 2,NZ-1
	DO J = 2,NY-1  
	DO I = 2,NX-1 

         IF(HARMONIC_AVE) THEN
! max{(D_{i-1,j}*D_{I,j})/(D_{i-1,j}+D_{I,j}),0}
            DIFF_I = 0.5*(FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I-1,J,K,IG)) &
                                + FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I+1,J,K,IG)))/(DX**2) &
                   + 0.5*(FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J-1,K,IG)) &
                                + FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J+1,K,IG)))/(DY**2) &
                   + 0.5*(FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J,K-1,IG)) &
                                + FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J,K+1,IG)))/(DZ**2) 

            A_I_J_K = DIFF_I + SIGMAT(I,J,K,IG) + (1./(DT*V_N(IG))) 
            A_IP1_J_K =  - 0.5*FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I+1,J,K,IG)) /(DX**2) 
            A_IM1_J_K =  - 0.5*FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I-1,J,K,IG)) /(DX**2)
            A_I_JP1_K =  - 0.5*FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J+1,K,IG)) /(DY**2) 
            A_I_JM1_K =  - 0.5*FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J-1,K,IG)) /(DY**2) 
            A_I_J_KP1 =  - 0.5*FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J,K+1,IG)) /(DZ**2) 
            A_I_J_KM1 =  - 0.5*FUN_HARMONIC(KDIFF(I,J,K,IG),KDIFF(I,J,K-1,IG)) /(DZ**2) 
         ELSE ! simple mean...
            DIFF_I = 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I-1,J,K,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I+1,J,K,IG),0.0))/(DX**2) &
                   + 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I,J-1,K,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I,J+1,K,IG),0.0))/(DY**2) &
                   + 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I,J,K-1,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I,J,K+1,IG),0.0))/(DZ**2) 

            A_I_J_K = DIFF_I + SIGMAT(I,J,K,IG) + (1./(DT*V_N(IG))) 
            A_IP1_J_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I+1,J,K,IG),0.0) /(DX**2) 
            A_IM1_J_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I-1,J,K,IG),0.0) /(DX**2)
            A_I_JP1_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J+1,K,IG),0.0) /(DY**2) 
            A_I_JM1_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J-1,K,IG),0.0) /(DY**2) 
            A_I_J_KP1 =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J,K+1,IG),0.0) /(DZ**2) 
            A_I_J_KM1 =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J,K-1,IG),0.0) /(DZ**2) 
         ENDIF

            A_MAT(I,J,K,IG,0,0,0)   = A_I_J_K  
            A_MAT(I,J,K,IG,1,0,0) = A_IP1_J_K
            A_MAT(I,J,K,IG,-1,0,0) = A_IM1_J_K
            A_MAT(I,J,K,IG,0,1,0) = A_I_JP1_K
            A_MAT(I,J,K,IG,0,-1,0) = A_I_JM1_K

            A_MAT(I,J,K,IG,0,0,1)  = A_I_J_KP1
            A_MAT(I,J,K,IG,0,0,-1) = A_I_J_KM1

!            B_MAT(I,J,K,IG,0,0,0) = 1./(V_N(IG)*DT) 

	END DO
	END DO 
	END DO 
	END DO 

! Discretize advection...

! Discretize advection...
      IF(NDIM_VEL.NE.0) THEN

        IF((NZ==3).AND.(NDIM_VEL==3)) THEN
           PRINT *,'INCOMPATABLE NZ AND NEW_VEL' 
           STOP 292
        ENDIF

        U_DIM=0.0; T_UPWIND=0.0; T_LIMIT=0.0; H_SIG_U=0.0

        S_INC_LIM=0.0

	DO IG = 1,NG
	DO K = 2,NZ
!        KSWITCH_UP = MAX( 1 - MAX(1,ABS(K-(NZ-1))) , 1 - MAX(1,ABS(K-2))  )
        SWITCH_UP1(3) = 1 - MIN(1,ABS(K-2))  
        SWITCH_UP2(3) = 1 - MIN(1,ABS(K-NZ)) 
        KM = K+SWITCH_UP1(3)  - SWITCH_UP2(3)
	DO J = 2,NY
!        JSWITCH_UP = MAX( 1 - MAX(1,ABS(J-(NY-1))) , 1 - MAX(1,ABS(J-2))  )
        SWITCH_UP1(2) = 1 - MIN(1,ABS(J-2))  
        SWITCH_UP2(2) = 1 - MIN(1,ABS(J-NY)) 
        JM = J+SWITCH_UP1(2)  - SWITCH_UP2(2)
	DO I = 2,NX
! advection upwards
!          DENOIN1 = SIGN( MAX(ABS( T(I+1) - T(I-1) ), TOLER), T(I+1) - T(I-1) )
! Now adjust to use upwind neer the boundaries...
!          ISWITCH_UP = MAX( 1 - MAX(1,ABS(I-(NX-1))) , 1 - MAX(1,ABS(I-2))  )
          SWITCH_UP1(1) = 1 - MIN(1,ABS(I-2))  
          SWITCH_UP2(1) = 1 - MIN(1,ABS(I-NX)) 
          IM = I+SWITCH_UP1(1)  - SWITCH_UP2(1)

          DO IDIM = 1, NDIM_VEL

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

!          T_LIMIT(I+1) =          T(I-1) + MAX(  MIN(FTILIN, XI_LIMIT*CTILIN, 1.0), CTILIN) * DENOIN 
!          T_LIMIT1(I+1) =          T(I-1) + MAX(  MIN(FTILIN, XI_LIMIT*CTILIN, 1.0), CTILIN) * ( T(I+1) - T(I-1) )
             T_LIMIT1 =          T_NEW(I1,J1,K1,IG)    + MAX(  MIN(FTILIN1, XI_LIMIT*CTILIN1, 1.0), CTILIN1) &
                                                           * ( T_NEW(I3,J3,K3,IG)    - T_NEW(I1,J1,K1,IG) )
             T_LIMIT2 =          T_NEW(II1,JJ1,KK1,IG) + MAX(  MIN(FTILIN2, XI_LIMIT*CTILIN2, 1.0), CTILIN2) &
                                                           * ( T_NEW(II3,JJ3,KK3,IG) - T_NEW(II1,JJ1,KK1,IG) )

             U_DIM(IDIM)=U(IDIM,I,J,K,IG)
             RH_SIG_U = 0.5 + 0.5*SIGN( 1.0, U_DIM(IDIM) ) 
             T_UPWIND(IDIM) = RH_SIG_U * T_NEW(I-ISWITCH,J-JSWITCH,K-KSWITCH,IG)  +  (1.0-RH_SIG_U) * T_NEW(I,J,K,IG) 
             T_LIMIT(IDIM)  = RH_SIG_U * ( T_UPWIND(IDIM)*REAL(switch_up1(IDIM)) +T_LIMIT1*(1.0-REAL(switch_up1(IDIM))) )  &
                     +  (1.0-RH_SIG_U) * ( T_UPWIND(IDIM)*REAL(switch_up2(IDIM)) +T_LIMIT2*(1.0-REAL(switch_up2(IDIM))) ) 
! 
!             T_COEF(IDIM) = U_DIM(IDIM)*T_LIMIT/ ( SIGN(1.0,T_UPWIND)* MAX(TOLER, ABS(T_UPWIND) ) )
             H_SIG_U(IDIM)=RH_SIG_U

           END DO ! DO IDIM = 1, NDIM_VEL

            A_MAT(I,J,K,IG,0,0,0)  = A_MAT(I,J,K,IG,0,0,0)   &
                    - (1.-H_SIG_U(1))*U_DIM(1)/DX - (1.-H_SIG_U(2))*U_DIM(2)/DY - (1.-H_SIG_U(3))*U_DIM(3)/DZ
            A_MAT(I,J,K,IG,-1,0,0) = A_MAT(I,J,K,IG,-1,0,0) - H_SIG_U(1)*U_DIM(1)/DX
            A_MAT(I,J,K,IG,0,-1,0) = A_MAT(I,J,K,IG,0,-1,0) - H_SIG_U(2)*U_DIM(2)/DY 
            A_MAT(I,J,K,IG,0,0,-1) = A_MAT(I,J,K,IG,0,0,-1) - H_SIG_U(3)*U_DIM(3)/DZ
 
            A_MAT(I-1,J,K,IG,0,0,0)= A_MAT(I-1,J,K,IG,0,0,0)+ H_SIG_U(1)*U_DIM(1)/DX
            A_MAT(I-1,J,K,IG,1,0,0)= A_MAT(I-1,J,K,IG,1,0,0)+ (1.-H_SIG_U(1))*U_DIM(1)/DX
 
            A_MAT(I,J-1,K,IG,0,0,0)= A_MAT(I,J-1,K,IG,0,0,0)+ H_SIG_U(2)*U_DIM(2)/DY
            A_MAT(I,J-1,K,IG,0,1,0)= A_MAT(I,J-1,K,IG,0,1,0)+ (1.-H_SIG_U(2))*U_DIM(2)/DY
 
            A_MAT(I,J,K-1,IG,0,0,0)= A_MAT(I,J,K-1,IG,0,0,0)+ H_SIG_U(3)*U_DIM(3)/DZ
            A_MAT(I,J,K-1,IG,0,0,1)= A_MAT(I,J,K-1,IG,0,0,1)+ (1.-H_SIG_U(3))*U_DIM(3)/DZ


!            A_MAT_I_J_K(I,J,K,IG)   = A_MAT_I_J_K(I,J,K,IG)    &
!                    - (1.-H_SIG_U(1))*U_DIM(1)/DX - (1.-H_SIG_U(2))*U_DIM(2)/DY - (1.-H_SIG_U(3))*U_DIM(3)/DZ
!            A_MAT_IM1_J_K(I,J,K,IG) =  A_MAT_IM1_J_K(I,J,K,IG) - H_SIG_U(1)*U_DIM(1)/DX
!            A_MAT_I_JM1_K(I,J,K,IG) =  A_MAT_I_JM1_K(I,J,K,IG) - H_SIG_U(2)*U_DIM(2)/DY 
!            A_MAT_I_J_KM1(I,J,K,IG) =  A_MAT_I_J_KM1(I,J,K,IG) - H_SIG_U(3)*U_DIM(3)/DZ
! 
!            A_MAT_I_J_K(I-1,J,K,IG)   = A_MAT_I_J_K(I-1,J,K,IG)    + H_SIG_U(1)*U_DIM(1)/DX
!            A_MAT_IP1_J_K(I-1,J,K,IG) =  A_MAT_IP1_J_K(I-1,J,K,IG) + (1.-H_SIG_U(1))*U_DIM(1)/DX

!            A_MAT_I_J_K(I,J-1,K,IG)   = A_MAT_I_J_K(I,J-1,K,IG)    + H_SIG_U(2)*U_DIM(2)/DY
!            A_MAT_I_JP1_K(I,J-1,K,IG) =  A_MAT_I_JP1_K(I,J-1,K,IG) + (1.-H_SIG_U(2))*U_DIM(2)/DY
 
!            A_MAT_I_J_K(I,J,K-1,IG)   = A_MAT_I_J_K(I,J,K-1,IG)    + H_SIG_U(3)*U_DIM(3)/DZ
!            A_MAT_I_J_KP1(I,J,K-1,IG) =  A_MAT_I_J_KP1(I,J,K-1,IG) + (1.-H_SIG_U(3))*U_DIM(3)/DZ

! For the vectors...
         IF(.NOT.UPWIND1ST) THEN
            S_INC_LIM(I,J,K,IG)   = S_INC_LIM(I,J,K,IG)    &
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

!        S_INC_LIM=S_INC_LIM + S
!      ELSE ! ENDOF IF(NDIM_VEL.NE.0) THEN
!        S_INC_LIM=S
      ENDIF ! IF(NDIM_VEL.NE.0) THEN ... ELSE...


	  END SUBROUTINE GET_MATRIX_TIME_STEPING_DIFFUSION_CALC
! 
! 

        SUBROUTINE ONE_D_LIMIT(T_LIMIT, T, NONODS) 
    ! This sub calculates the limited face values T_LIMIT of variable T. 
    ! VEL is an input and contains the direction of velocity (+ve or -ve).  
    ! Also returned is the ratio T_RATIO which is: the limited value/ the upwind value (=1 for upwinding).  
    ! Only return the ratio for the time being. 
    !--------------------------------------------------- 
    implicit none
    INTEGER, intent( in ) :: NONODS
    LOGICAL UPWIND
    REAL, DIMENSION( NONODS+1), intent( inout ) :: T_LIMIT
    REAL, DIMENSION( NONODS ), intent( in ) :: T
    ! Local variables
    REAL, PARAMETER :: TOLER=1.0E-10, XI_LIMIT = 2.0 ! defines TVD curve on the NVD diagram. 
    REAL :: DENOIN,CTILIN,DENOOU,CTILOU,FTILIN,FTILOU
    INTEGER :: K
! cell numbering for +ve vel:
! k-1 ! k !  k+1
!   U ! C !f D
! \hat I=(I -I_U)/(I_D - I_U)  

      
! Calculate away from the boundaries (top and bottom) 
      T_limit = 0.0 ! This assumes an upwind approximation. 
      DO K=2,NONODS-1

! advection upwards
          DENOIN = SIGN( MAX(ABS( T(K+1) - T(K-1) ), TOLER), T(K+1) - T(K-1) )
          CTILIN = ( T(K) - T(K-1) ) / DENOIN

          FTILIN = ( 0.5*(T(K)+T(K+1)) - T(K-1) ) / DENOIN ! the mean with 0.5 coeff is the high order value of T. 

! Velocity is going out of element
!          T_LIMIT(K+1) =          T(K-1) + MAX(  MIN(FTILIN, XI_LIMIT*CTILIN, 1.0), CTILIN) * DENOIN 
          T_LIMIT(K+1) =          T(K-1) + MAX(  MIN(FTILIN, XI_LIMIT*CTILIN, 1.0), CTILIN) * ( T(K+1) - T(K-1) )
!          print *,'k,t(k-1),t(k),T_LIMIT(K+1),t(k+1):',k,t(k-1),t(k),T_LIMIT(K+1),t(k+1)

      END DO
!       stop 3

    RETURN

  END SUBROUTINE ONE_D_LIMIT

! 
       REAL FUNCTION FUN_HARMONIC(K1,K2)
! return the Harmonic average of K1,K2 but adjusted to take into account negative K1,K2 that are used to 
! set a zero diffusivity.
       REAL K1,K2
       REAL, PARAMETER :: TOLER=1.E-15
       FUN_HARMONIC = 2.*ABS(K1*K2)/MAX(TOLER, K1+K2)
       END FUNCTION FUN_HARMONIC


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
         PRINT *,'ITS_NG,R_T_MAX_DIF_NG,ERROR_SOLV_NG:',ITS_NG,R_T_MAX_DIF_NG,ERROR_SOLV_NG

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






! *********************************************************************************************
! *********************************************************************************************
! **************************************** ANN SUBROUTINES ************************************
! *********************************************************************************************
! *********************************************************************************************


!
!
       SUBROUTINE NEURAL_NET_BATCHES(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,NONODS,&
                  NEXAMPLE,NBATCH,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY,NITS,&
                  ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER,&
                  IDEAL_MIDLE,NDEAL_MIDLE )
! Train (IF nits>0) OR SOLVE FOR THE NEURON VALUES OF AN ANN. 
! NITS= no of iterations of the ANN
! ERROR_TOLER = tolerance of the ANN. 
! This sub forms the neural network for a set of NEXAMPLE problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! If NOREXP=0 find exponent of OUTPUT NEURONS otherwise dont.

! THE ITERATION PARAMETERS: ALPHA,NTEST,TEST_BELOW,TEST_ABOVE,ANNEAL_FRAC
! Default values: NITS=10000,ERROR_TOLER=1.e-5, ALPHA=0.001,ANNEAL_FRAC=0.999
! ALPHA= initial value of the max value to add to the weights.
! ANNEAL_FRAC=what fraction to reduce ALPHA by if we are not converging well
! WEIGHT_DECAY is the weight magnitude penalty term
! 
       IMPLICIT NONE
       INTEGER, PARAMETER :: W_TOLER=0.5
       INTEGER, intent( inout ) :: NBATCH
       INTEGER, intent( in ) :: NCOLM,NONODS,NLAYER,NLAY_IN,NLAY_OUT,NEXAMPLE,nits,NDEAL_MIDLE
       INTEGER, intent( in ) :: NLAY(NLAYER),NOREXP
       REAL, intent( in ) :: IDEAL_INPUT(NLAY_IN,NEXAMPLE),IDEAL_OUTPUT(NLAY_OUT,NEXAMPLE),&
                             IDEAL_MIDLE(NDEAL_MIDLE,NEXAMPLE),WEIGHT_DECAY
       REAL, intent( inout ) :: WEIGHT(NCOLM)
       REAL, intent( inout )  :: NEURAL(NONODS,NEXAMPLE),FF
       REAL, intent( out )  :: ERROR_TOLER
       REAL, intent( inout ) :: ALPHA
       REAL, intent( in ) :: ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER
! local variables
       INTEGER :: ITS,NPRINT,IGOT_BETTER,I,J,IEXAMPLE,ITS2,NITS2
       REAL :: FF_KEEP, ACC_TOLER, RAN_NO, ff_best_q
       LOGICAL :: SWITCH, only_best
       REAL, ALLOCATABLE :: WEIGHT_KEEP(:),RAN_WEIGHT(:)
       REAL, ALLOCATABLE :: NEURAL2(:,:)
       REAL, ALLOCATABLE :: IDEAL_INPUT2(:,:),IDEAL_OUTPUT2(:,:),IDEAL_MIDLE2(:,:)
       LOGICAL, ALLOCATABLE :: IN_EXAMPLE(:)
       INTEGER, ALLOCATABLE :: BATCH_EXAMPLE(:)

       ALLOCATE(WEIGHT_KEEP(NCOLM),RAN_WEIGHT(NCOLM))
       ALLOCATE(IN_EXAMPLE(NEXAMPLE))


       only_best=.true. ! only keep a batch if it improves the iteration

          IF((NBATCH.NE.NEXAMPLE).AND.(NBATCH.NE.0)) THEN


       NITS2=INT(SQRT(REAL(NITS)))
       DO ITS2=1,NITS2
            PRINT *,'ITS2,INT(SQRT(REAL(NITS)))=',ITS2,NITS2
       ALLOCATE(NEURAL2(NONODS,NBATCH))
       ALLOCATE(IDEAL_INPUT2(NLAY_IN,NBATCH),IDEAL_OUTPUT2(NLAY_OUT,NBATCH),IDEAL_MIDLE2(NDEAL_MIDLE,NBATCH))
       ALLOCATE(BATCH_EXAMPLE(NBATCH))
!          IF(.true.) THEN
             if(.false.) then
             DO I=1,NBATCH
                IEXAMPLE=i
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             else
             DO I=1,NBATCH
                DO J=1,10000
!                    print *,'j=',j
                   CALL RANDOM_NUMBER(RAN_NO)
                   IEXAMPLE = MIN(NEXAMPLE,MAX(1, INT(RAN_NO*REAL(NEXAMPLE)+1.0)))
                   IF(.NOT.IN_EXAMPLE(IEXAMPLE)) THEN
                      IN_EXAMPLE(IEXAMPLE) = .TRUE.
                      EXIT
                   ENDIF
                END DO
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             endif
       if(only_best) then
       CALL NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,&
            NONODS,NEXAMPLE,NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,&
            WEIGHT_DECAY,1,ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,&
            ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
       ff_best_q=ff
       weight_keep = weight
       endif

       CALL NEURAL_NET(NEURAL2,FF,IDEAL_INPUT2,IDEAL_OUTPUT2,WEIGHT,NCOLM,NONODS,&
            NBATCH,NBATCH,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY,&
            INT(SQRT(REAL(NITS))),ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,&
            ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
             NEURAL(:,BATCH_EXAMPLE(:))=NEURAL2(:,:) ! Store the results
             IN_EXAMPLE(BATCH_EXAMPLE(:))=.FALSE. ! set everything back to .false.

       if(only_best) then
       CALL NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,NONODS,NEXAMPLE,&
            NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY,&
            1,ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
           print *,'ff,ff_best_q,nbatch:',ff,ff_best_q,nbatch
       if(ff.gt.ff_best_q) weight = weight_keep
       if(ff.gt.ff_best_q) print *,'****reset the weights'
          if(ff.gt.ff_best_q) then
             nbatch=nbatch+1
             if(nbatch.gt.nexample/2) nbatch=nexample
          else
             nbatch=nbatch-1
             if(nbatch.gt.nexample/2) nbatch=nexample/2
          endif
          print *,'new nbatch=',nbatch
       endif

       DEALLOCATE(NEURAL2, IDEAL_INPUT2, IDEAL_OUTPUT2,IDEAL_MIDLE2, BATCH_EXAMPLE )

       END DO ! ITS2

          ELSE

       CALL NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,&
            NONODS,NEXAMPLE,NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY, &
            NITS,ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
          ENDIF

       RETURN
       END SUBROUTINE NEURAL_NET_BATCHES
!
!
!
!
!
!
       SUBROUTINE NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,NONODS,&
                             &NEXAMPLE,NBATCH,NLAY,NLAYER,NLAY_IN,NLAY_OUT,&
                             &NOREXP,WEIGHT_DECAY,&
                             &NITS,ERROR_TOLER,ALPHA, ANNEAL_FRAC_SMALLER,&
                             &ANNEAL_FRAC_BIGGER,IDEAL_MIDLE,NDEAL_MIDLE)
! Train (IF nits>0) OR SOLVE FOR THE NEURON VALUES OF AN ANN. 
! NITS= no of iterations of the ANN
! ERROR_TOLER = tolerance of the ANN. 
! This sub forms the neural network for a set of NEXAMPLE problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! If NOREXP=0 find exponent of OUTPUT NEURONS otherwise dont.

! THE ITERATION PARAMETERS: ALPHA,NTEST,TEST_BELOW,TEST_ABOVE,ANNEAL_FRAC
! Default values: NITS=10000,ERROR_TOLER=1.e-5, ALPHA=0.001,ANNEAL_FRAC=0.999
! ALPHA= initial value of the max value to add to the weights.
! ANNEAL_FRAC=what fraction to reduce ALPHA by if we are not converging well
! WEIGHT_DECAY is the weight magnitude penalty term
! 
       IMPLICIT NONE
       INTEGER, PARAMETER :: W_TOLER=0.5
       INTEGER, intent( in ) :: NCOLM,NONODS,NLAYER,NLAY_IN,NLAY_OUT,NEXAMPLE,NBATCH,nits,NDEAL_MIDLE
       INTEGER, intent( in ) :: NLAY(NLAYER),NOREXP
       REAL, intent( in ) :: IDEAL_INPUT(NLAY_IN,NEXAMPLE),IDEAL_OUTPUT(NLAY_OUT,NEXAMPLE),IDEAL_MIDLE(NDEAL_MIDLE,NEXAMPLE)
       REAL, intent( in ) :: WEIGHT_DECAY
       REAL, intent( inout ) :: WEIGHT(NCOLM)
       REAL, intent( inout )  :: NEURAL(NONODS,NEXAMPLE),FF
       REAL, intent( out )  :: ERROR_TOLER
       REAL, intent( inout ) :: ALPHA
       REAL, intent( in ) :: ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER
! local variables
       INTEGER :: ITS,NPRINT,IGOT_BETTER,I,J,IEXAMPLE
       REAL :: FF_KEEP, ACC_TOLER, RAN_NO
       LOGICAL :: SWITCH
       REAL, ALLOCATABLE :: WEIGHT_KEEP(:),RAN_WEIGHT(:)
       REAL, ALLOCATABLE :: NEURAL2(:,:)
       REAL, ALLOCATABLE :: IDEAL_INPUT2(:,:),IDEAL_OUTPUT2(:,:),IDEAL_MIDLE2(:,:)
       LOGICAL, ALLOCATABLE :: IN_EXAMPLE(:)
       INTEGER, ALLOCATABLE :: BATCH_EXAMPLE(:)

       ALLOCATE(WEIGHT_KEEP(NCOLM),RAN_WEIGHT(NCOLM))
       ALLOCATE(NEURAL2(NONODS,NBATCH))
       ALLOCATE(IN_EXAMPLE(NEXAMPLE))
       ALLOCATE(IDEAL_INPUT2(NLAY_IN,NBATCH),IDEAL_OUTPUT2(NLAY_OUT,NBATCH),IDEAL_MIDLE2(NDEAL_MIDLE,NBATCH))
       ALLOCATE(BATCH_EXAMPLE(NBATCH))

       IN_EXAMPLE=.FALSE.

       CALL NEURAL_EXAMPLE(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,&
            WEIGHT,NCOLM,NONODS,NEXAMPLE,NLAY,NLAYER,NLAY_IN,&
            NLAY_OUT,NOREXP,WEIGHT_DECAY, IDEAL_MIDLE,NDEAL_MIDLE)
!            print *,'original ff=',ff

       ACC_TOLER=FF
       
       IGOT_BETTER=0
       SWITCH=.FALSE.
       DO ITS=1,NITS ! Train the ANN if NITS>0
          WEIGHT_KEEP=WEIGHT
          IF(SWITCH) THEN
             WEIGHT = WEIGHT - (RAN_WEIGHT-0.5)*ALPHA
          ELSE
             CALL RANDOM_NUMBER(RAN_WEIGHT)
!            alpha=1.e-5
             WEIGHT = WEIGHT + (RAN_WEIGHT-0.5)*ALPHA
          ENDIF
!          print *,'RAN_WEIGHT-0.5:',RAN_WEIGHT-0.5
!          stop 11

          FF_KEEP = FF
          IF((NBATCH.NE.NEXAMPLE).AND.(NBATCH.NE.0)) THEN
!          IF(.true.) THEN
             if(.false.) then
             DO I=1,NBATCH
                IEXAMPLE=i
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             else
             DO I=1,NBATCH
                DO J=1,10000
!                    print *,'j=',j
                   CALL RANDOM_NUMBER(RAN_NO)
                   IEXAMPLE = MIN(NEXAMPLE,MAX(1, INT(RAN_NO*REAL(NEXAMPLE)+1.0)))
                   IF(.NOT.IN_EXAMPLE(IEXAMPLE)) THEN
                      IN_EXAMPLE(IEXAMPLE) = .TRUE.
                      EXIT
                   ENDIF
                END DO
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             endif
!             print *,'Batch_example:',batch_example
!             print *,'NEURAL2:',NEURAL2
!             print *,'IDEAL_INPUT2:',IDEAL_INPUT2
!             print *,'IDEAL_OUTPUT2:',IDEAL_OUTPUT2
             CALL NEURAL_EXAMPLE(NEURAL2,FF,IDEAL_INPUT2,IDEAL_OUTPUT2,&
                  WEIGHT,NCOLM,NONODS,NBATCH,NLAY,NLAYER,NLAY_IN,&
                  NLAY_OUT,NOREXP,WEIGHT_DECAY, IDEAL_MIDLE2,NDEAL_MIDLE)
             NEURAL(:,BATCH_EXAMPLE(:))=NEURAL2(:,:) ! Store the results
             IN_EXAMPLE(BATCH_EXAMPLE(:))=.FALSE. ! set everything back to .false.
!             ff_keep = 0.9*ff_keep + 0.1*ff
!             print *,'ff=',ff
          ELSE
             CALL NEURAL_EXAMPLE(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,&
                  WEIGHT,NCOLM,NONODS,NEXAMPLE,NLAY,NLAYER,NLAY_IN,&
                  NLAY_OUT,NOREXP,WEIGHT_DECAY, IDEAL_MIDLE,NDEAL_MIDLE)
          ENDIF

          ACC_TOLER = (1.-W_TOLER) * ACC_TOLER + W_TOLER * ABS(FF-FF_KEEP) ! Make the tolerance change slowly.

          IF(ACC_TOLER<ERROR_TOLER) CYCLE
!          print *,'its,ff,ff_keep:',its,ff,ff_keep

          IF(FF<FF_KEEP) THEN
             IF(SWITCH) ALPHA=ALPHA/ANNEAL_FRAC_BIGGER  ! make bigger
             SWITCH=.FALSE.
             IGOT_BETTER=IGOT_BETTER+1
!              print *,'****'
          ELSE
             IF(SWITCH) ALPHA=ALPHA*ANNEAL_FRAC_SMALLER  ! make smaller
             SWITCH=.NOT.SWITCH
             FF=FF_KEEP
             WEIGHT=WEIGHT_KEEP
          ENDIF
!          print *,'alpha=',alpha

!          NPRINT=100000
!          NPRINT=100000
!          NPRINT=1000
          NPRINT=100
          IF(MOD(ITS,NPRINT)==0) THEN ! REAL(IGOT_BETTER)/REAL(NPRINT)=0.66 WHEN ALPHA=SMALL NO
             print *,'its,ff,ALPHA:',its,ff,ALPHA,REAL(IGOT_BETTER)/REAL(NPRINT)
             IGOT_BETTER=0
          ENDIF

       END DO ! DO ITS=1,NITS
!
       RETURN
       END SUBROUTINE NEURAL_NET
!
!
!
!
       SUBROUTINE NEURAL_EXAMPLE(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,&
                  WEIGHT,NCOLM,NONODS,NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,&
                  NOREXP,WEIGHT_DECAY, IDEAL_MIDLE,NDEAL_MIDLE)
! This sub forms the neural network for a set of NEXAMPL problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! If NOREXP=0 find exponent of OUTPUT NEURONS otherwise dont.
! WEIGHT_DECAY is the weight magnitude penalty term
! 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NCOLM,NONODS,NLAYER,NLAY_IN,NLAY_OUT,NEXAMPLE,NDEAL_MIDLE
       INTEGER, intent( in ) :: NLAY(NLAYER),NOREXP
       REAL, intent( in ) :: IDEAL_INPUT(NLAY_IN,NEXAMPLE),IDEAL_OUTPUT(NLAY_OUT,NEXAMPLE),WEIGHT(NCOLM),&
                             WEIGHT_DECAY, IDEAL_MIDLE(NDEAL_MIDLE,NEXAMPLE)
       REAL, intent( inout )  :: NEURAL(NONODS,NEXAMPLE),FF
! Local varibales...
       INTEGER :: IEXAMPLE,ISTART
       REAL, ALLOCATABLE :: W_EXP(:)

       ALLOCATE(W_EXP(NLAYER)) 

       W_EXP(:)=0.0
       IF(NOREXP==0) W_EXP(NLAYER)=1.0 ! no exponential in final layer...
       IF(NOREXP==0) W_EXP(int(NLAYER/2) +1)=1.0 ! no exponential in final layer...

       ISTART = SUM(NLAY(1:NLAYER/2)) +1

       FF=0.0
       DO IEXAMPLE=1,NEXAMPLE
          NEURAL(1:NLAY_IN,IEXAMPLE)=IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
          CALL GETNEUVALS_FAST(NEURAL(:,IEXAMPLE),WEIGHT,NONODS,NCOLM,W_EXP,NLAY,NLAYER)

          FF=FF + SUM( (NEURAL(NONODS-NLAY_OUT+1:NONODS,IEXAMPLE)-IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE))**2 )
          FF=FF + 0.*SUM( (NEURAL(ISTART:ISTART+NDEAL_MIDLE-1,IEXAMPLE)-IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE))**2 )
       END DO ! DO IEXAMPLE=1,NEXAMPLE
       FF=FF/REAL(NEXAMPLE)
       FF=FF+(WEIGHT_DECAY/REAL(NCOLM)) * SUM(WEIGHT(:)**2) 
!
       RETURN
       END SUBROUTINE NEURAL_EXAMPLE
!
!
!
!
       SUBROUTINE GETNEUVALS_FAST(NEUVAL,WEIGHT,NONODS,NCOLM,W_EXP,NLAY,NLAYER)
! This sub forms the neural network for a set of NEXAMPL problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! This sub calculates the neuron values. 
! If W_EXP=0.0 find exponent of NEURONS otherwise dont. 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NONODS,NCOLM, NLAYER
       REAL, intent( in ) :: WEIGHT(NCOLM),W_EXP(NLAYER)
       REAL, intent( inout ) :: NEUVAL(NONODS)
       INTEGER , intent( in ):: NLAY(NLAYER)
! LOCAL VARIABLES...
       REAL :: SUMWEI_VAL
       INTEGER :: NOD,ILAY,ILAY1,ILAY2,NLAYACC_WEIT,NLAYACC_NOD,NLAYACC_NOD_PREV, IWEI_ADD
       
       NLAYACC_WEIT=0
       NLAYACC_NOD_PREV=0
       NLAYACC_NOD=NLAY(1)
       DO ILAY=2,NLAYER

          DO ILAY2=1,NLAY(ILAY) 
             SUMWEI_VAL=0.0
             IWEI_ADD = NLAYACC_WEIT + (ILAY2-1)*NLAY(ILAY-1)
             DO ILAY1=1,NLAY(ILAY-1)
                SUMWEI_VAL=SUMWEI_VAL+WEIGHT(IWEI_ADD + ILAY1)*NEUVAL(ILAY1 + NLAYACC_NOD_PREV)
             END DO
             NEUVAL(NLAYACC_NOD+ILAY2)= W_EXP(ILAY)*SUMWEI_VAL   +   (1.-W_EXP(ILAY))/(1.+EXP(-SUMWEI_VAL))
          END DO
          NLAYACC_NOD_PREV=NLAYACC_NOD_PREV+NLAY(ILAY-1)
          NLAYACC_NOD=NLAYACC_NOD+NLAY(ILAY)
          NLAYACC_WEIT=NLAYACC_WEIT + NLAY(ILAY) * NLAY(ILAY-1)

       END DO

       RETURN
       END SUBROUTINE GETNEUVALS_FAST
!
!


      subroutine road_markings(material,nx)
      integer material(nx) 

      integer count_one, count_two, count_three
      logical create_one, create_two, create_three
      integer i

      material = 0
!      print *, 'material', material

      count_one = 0
      count_two = 0
      count_three = 0

      create_one = .True.
      create_two = .False.
      create_three = .False.

      i = 3

      do 
          !i = i + 1
!          print *, "i=", i

          ! i  road, i+1 space, i+2 next location
          if (create_one) then
              material(i) = 1
              i = i + 2
              count_one = count_one + 1
          endif 


          if ((count_one==3).and.(count_two.ne.3)) then

              ! create_two
              ! i  road, i+1 road, i+2 space, i+3 next location
              count_one = 0
              material(i) = 1
              material(i+1) = 1
              i = i + 3
              count_two = count_two + 1

          else if ((count_one==3).and.(count_two==3)) then

              ! create_three
              ! i,i+1,i+2  road,, i+3 space, i+4 next location
              count_one = 0
              count_two = 0
              material(i) = 1
              material(i+1) = 1
              material(i+2) = 1
              i = i + 4
              count_three = count_three + 1

          endif



          if (i>nx-3) exit
      enddo

!       write(*,'(A, 50I2)')  'material', material
      end subroutine road_markings



       SUBROUTINE DEFINE_GROUPS_MATERIALS(HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
               HOSPITAL, SCHOOL, OFFICES, SHOPS, &
               NX,NY,NZ, NG, MATERIAL, rd_mark_x, rd_mark_y)  
! This sub defines what the regions are. 
       INTEGER, intent( in ) :: HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
               HOSPITAL, SCHOOL, OFFICES, SHOPS
       INTEGER, intent( in ) :: NX,NY,NZ, NG
       INTEGER, intent( inout ) :: MATERIAL(NX,NY,NZ),  rd_mark_x(NX), rd_mark_y(NY) 
! Local variables...
       LOGICAL SIMPLE
       PARAMETER(SIMPLE=.FALSE.) ! Turn into a simple problem...
       INTEGER I,J,K, IDISP,JDISP
       INTEGER, ALLOCATABLE :: rd_mark(:)

         allocate(rd_mark(max(nx,ny)) )
!         allocate(rd_mark_x(nx),rd_mark_y(ny) )

         call road_markings(rd_mark,max(nx,ny)) 
! modify so that there is a main road around the perimeter of the domain... 
         rd_mark_x(1:nx)=rd_mark(1:nx) 
         rd_mark_y(1:ny)=rd_mark(1:ny) 
         rd_mark_x(2:4)=1; rd_mark_x(nx-3:nx-1)=1
         rd_mark_y(2:4)=1; rd_mark_y(ny-3:ny-1)=1
         

         MATERIAL=HOMES_OCCUPIED
         DO I=2, nx-1 
! Draw roads in a cross shape. 
            IF(RD_MARK_x(I)==1) MATERIAL(I, :, :) = ROADS 
         END DO

         DO J=2,ny-1
! Draw roads in a cross shape. 
            IF(RD_MARK_y(J)==1) MATERIAL(:, J, :) = ROADS 
         END DO

         DO i=1,nx-1
         DO j=1,ny-1
! Draw roads in a cross shape. 
            IF((RD_MARK_x(I)==1).AND.(RD_MARK_y(J)==1)) MATERIAL(I,J,:)=CROSS_ROADS
         END DO
         END DO

! Defined occupied homes...
         DO I=2,NX-1
         DO J=2,NY-1
         DO K=2,NZ-1
            IF((I<=1*NX/2).AND.(J<=1*NY/2)) THEN
               IF(MATERIAL(I,J,K)==HOMES) MATERIAL(I,J,K)=HOMES_OCCUPIED
            ENDIF
         END DO
         END DO
         END DO

! define place:
         IDISP=4
         JDISP=-3
! PARK
         DO I=2,NX-1
         DO J=2,NY-1
         DO K=2,NZ-1
            IF((I>=2*NX/4-1).AND.(I<=3*NX/4)) THEN
            IF((J>=1*NY/4+0).AND.(J<=3*NY/4-4)) THEN
               MATERIAL(I+IDISP+8,J+JDISP+10,K)=PARK
            ENDIF
            ENDIF
         END DO
         END DO
         END DO

! HOSPITAL
         DO I=2,NX-1
         DO J=2,NY-1
         DO K=2,NZ-1
            IF((I>=3*NX/4).AND.(I<=7*NX/8)) THEN
            IF((J>=2*NY/4+4).AND.(J<=3*NY/4)) THEN
            IF(I.LE.NX-20) THEN
               MATERIAL(I+IDISP+12,J+JDISP+6,K)=HOSPITAL
            ENDIF
            ENDIF
            ENDIF
         END DO
         END DO
         END DO

! school
         DO I=2,NX-1
         DO J=2,NY-1
         DO K=2,NZ-1
            IF((I>=1*NY/4-1).AND.(I<2*NY/4)) THEN
            IF((J>=3*NX/4).AND.(J<=7*NX/8+1)) THEN
               MATERIAL(I+IDISP+4,J+JDISP+4,K)=SCHOOL
            ENDIF
            ENDIF
         END DO
         END DO
         END DO
! OFFICES
         DO I=2,NX-1
         DO J=2,NY-1
         DO K=2,NZ-1
            IF((I>=1*NY/4-1).AND.(I<2*NY/4)) THEN
            IF((J>=2*NX/4).AND.(J<=5*NX/8+1)) THEN
               MATERIAL(I+IDISP+4,J+JDISP+4,K)=OFFICES
            ENDIF
            ENDIF
         END DO
         END DO
         END DO
! OFFICES
         DO I=2,NX-1
         DO J=2,NY-1
         DO K=2,NZ-1
            IF((I>=1*NY/4-1).AND.(I<2*NY/4)) THEN
            IF((J>=0*NX/4).AND.(J<=1*NX/8)) THEN
            IF(J.GE.9) THEN
               MATERIAL(I+IDISP+4,J+JDISP,K)=OFFICES
            ENDIF
            ENDIF
            ENDIF
         END DO
         END DO
         END DO
! SHOPS
         DO I=2,NX-1
         DO J=2,NY-1
         DO K=2,NZ-1
            IF((I>=0*NY/4).AND.(I<1*NY/4)) THEN
            IF((J>=2*NX/4).AND.(J<=6*NX/8)) THEN
            IF(I.GE.4) THEN
               MATERIAL(I+IDISP+0,J+JDISP+4,K)=SHOPS
            ENDIF
            ENDIF
            ENDIF
         END DO
         END DO
         END DO

! Turn into a simple problem...
         IF(SIMPLE) THEN
         IF(NG==8) THEN
            DO I=2,NX-1
            DO J=2,NY-1
            DO K=2,NZ-1
              IF( (MATERIAL(I,J,K)==SHOPS) &
              .OR.(MATERIAL(I,J,K)==OFFICES) &
              .OR.(MATERIAL(I,J,K)==SCHOOL) &
              .OR.(MATERIAL(I,J,K)==HOSPITAL) ) MATERIAL(I,J,K)=HOMES_OCCUPIED
            END DO
            END DO
            END DO
         ENDIF
         ENDIF
! tweak as there is a bug near boundary so take the road away here...
       if(.false.) then
!         rd_mark_x(nx-1:nx)=0
!         rd_mark_y(1:2)=0
!         rd_mark_x(nx-4:nx)=0
!         rd_mark_y(1:4)=0
         rd_mark_x(nx-10:nx)=0
         rd_mark_y(1:10)=0
       endif
!           print *,'rd_mark_x:',rd_mark_x
!           print *,'rd_mark_y:',rd_mark_y

         END SUBROUTINE DEFINE_GROUPS_MATERIALS





         SUBROUTINE DETERMINE_X_SECTIONS( KK,  SIGMAT,  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, &
                    SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, DX,DY,DZ,DT, &
                    OPTION,ITIME,ITS,ACCTIM,TIME_MAT, RDAY,EIGEN_VALUE,EIGEN_METHOD1, R_EIGEN, &
                    HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NEW, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
                    RECORD_S, SPAR_RECORD_S, RECORD_F, SPAR_RECORD_F,  &
                    R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK )
! Determine KK,  SIGMAT,  SIGMA_S_OFF,  SIGMA_F...
! We do this for a town. 
         IMPLICIT NONE
         INTEGER, intent( in ) :: NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, NDIM_VEL
         INTEGER, intent( in ) :: FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S),FIN_SIGMA_F(NG+1), COL_SIGMA_F(NCOL_SIGMA_F)
         REAL, intent( in ) :: DX,DY,DZ,DT
         REAL, intent( inout ) :: SIGMAT(NX,NY,NZ,NG),KK(NX,NY,NZ,NG),SIGMA_F(NX,NY,NZ,NG,NG2),SIGMA_F_PTER(NX,NY,NZ,NCOL_SIGMA_F)
         REAL, intent( inout ) :: SIGMA_S_OFF(NX,NY,NZ,NG,NG2),  SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S)
         REAL, intent( inout ) :: T_OLD(NX,NY,NZ,NG),T_NEW(NX,NY,NZ,NG)
         REAL, intent( in ) :: ACCTIM,TIME_MAT, RDAY, R_EIGEN 
         REAL, intent( inout ) :: S(NX,NY,NZ,NG),U(NDIM_VEL,NX+1,NY+1,NZ+1,NG)
         LOGICAL, intent( in ) :: EIGEN_VALUE, EIGEN_METHOD1
         INTEGER, intent( in ) :: HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
               HOSPITAL, SCHOOL, OFFICES, SHOPS
         INTEGER, intent( in ) :: ITIME,ITS, OPTION
         INTEGER, intent( in ) :: MATERIAL(NX,NY,NZ)
         INTEGER, intent( in ) :: rd_mark_x(NX), rd_mark_y(NY)
         LOGICAL, intent( in ) :: RECORD_S, RECORD_F
         LOGICAL, intent( inout ) :: SPAR_RECORD_S(NG,NG), SPAR_RECORD_F(NG,NG)
         REAL, intent( in ) :: R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, &
                               R0_SHOPS, R0_VEHC_1TO5, R0_PARK
! local varibales...

         IF(NG==8) THEN
           CALL DETERMINE_X_SECTIONS_NG8( KK,  SIGMAT,  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, &
                    SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, DX,DY,DZ,DT, &
                    OPTION,ITIME,ITS,ACCTIM,TIME_MAT, RDAY,EIGEN_VALUE,EIGEN_METHOD1, R_EIGEN, &
                    HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NEW, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
                    RECORD_S, SPAR_RECORD_S, RECORD_F, SPAR_RECORD_F,  &
                    R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK )
         ELSE IF(NG==260) THEN
           CALL DETERMINE_X_SECTIONS_NG260( KK,  SIGMAT,  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, &
                    SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, DX,DY,DZ,DT, &
                    OPTION,ITIME,ITS,ACCTIM,TIME_MAT, RDAY,EIGEN_VALUE,EIGEN_METHOD1, R_EIGEN, &
                    HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NEW, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
                    RECORD_S, SPAR_RECORD_S, RECORD_F, SPAR_RECORD_F,  &
                    R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK )

         ENDIF

!               stop 292
         END SUBROUTINE DETERMINE_X_SECTIONS





         SUBROUTINE DETERMINE_X_SECTIONS_NG8( KK,  SIGMAT,  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, &
                    SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, DX,DY,DZ,DT, &
                    OPTION,ITIME,ITS,ACCTIM,TIME_MAT, RDAY,EIGEN_VALUE,EIGEN_METHOD1, R_EIGEN, &
                    HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NEW, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
                    RECORD_S, SPAR_RECORD_S, RECORD_F, SPAR_RECORD_F,  &
                    R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK )
! Determine KK,  SIGMAT,  SIGMA_S_OFF,  SIGMA_F...
! We do this for a town. 
         IMPLICIT NONE
         INTEGER, intent( in ) :: NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, NDIM_VEL
         INTEGER, intent( in ) :: FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S),FIN_SIGMA_F(NG+1), COL_SIGMA_F(NCOL_SIGMA_F)
         REAL, intent( in ) :: DX,DY,DZ,DT
         REAL, intent( inout ) :: SIGMAT(NX,NY,NZ,NG),KK(NX,NY,NZ,NG),SIGMA_F(NX,NY,NZ,NG,NG2),SIGMA_F_PTER(NX,NY,NZ,NCOL_SIGMA_F)
         REAL, intent( inout ) :: SIGMA_S_OFF(NX,NY,NZ,NG,NG2),  SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S)
         REAL, intent( inout ) :: T_OLD(NX,NY,NZ,NG),T_NEW(NX,NY,NZ,NG)
         REAL, intent( in ) :: ACCTIM,TIME_MAT, RDAY, R_EIGEN 
         REAL, intent( inout ) :: S(NX,NY,NZ,NG),U(NDIM_VEL,NX+1,NY+1,NZ+1,NG)
         LOGICAL, intent( in ) :: EIGEN_VALUE, EIGEN_METHOD1
         INTEGER, intent( in ) :: HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
               HOSPITAL, SCHOOL, OFFICES, SHOPS
         INTEGER, intent( in ) :: ITIME,ITS, OPTION
         INTEGER, intent( in ) :: MATERIAL(NX,NY,NZ)
         INTEGER, intent( in ) :: rd_mark_x(NX), rd_mark_y(NY)
         LOGICAL, intent( in ) :: RECORD_S, RECORD_F
         LOGICAL, intent( inout ) :: SPAR_RECORD_S(NG,NG), SPAR_RECORD_F(NG,NG)
         REAL, intent( in ) :: R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, &
                               R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK
!local varibales...
         LOGICAL SIMPLE
         PARAMETER(SIMPLE=.FALSE.) ! Turn into a simple problem...
         REAL TOLER
         PARAMETER(TOLER=1.E-10) 
         REAL DURATION_INFECT, INCUBATION_DURATION, VIRUS_SIGMA, VIRUS_GAMMA
         REAL XI, R0_MOBILE, MU,NU 
         REAL BETA, BETA2, LENGTH
         REAL ONE_DAY, VIRUS_N1_AIM, VIRUS_N2_AIM, LAMBDA_H_H, FORCE, H2M, LAMBDA_M_M
         REAL R_RATIO, VIRUS_N1, VIRUS_N2, RSWITCH
         real rhome_max,rbuilding_max
         INTEGER I,J,K
         LOGICAL IS_BUILDING

! Key variables...
         one_day=24.*3600.0 
         DURATION_INFECT=7.0*one_day ! 7 day infection duration. 
         INCUBATION_DURATION=4.5*one_day! 4.5 day inCUBATION duration.
         VIRUS_SIGMA=1.0/INCUBATION_DURATION
         VIRUS_GAMMA=1.0/DURATION_INFECT
! ξ (xi) is the rate which recovered individuals return to the susceptible statue due to loss of immunity.
         XI= 1./(365.0*one_day) 
!         R0_HOME=5.0 !5.0 !0.2 ! No of people an infected person infects.
         R0_MOBILE=R0_DIFF1 !20.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
         MU=1.0*1./(60.*365.0*one_day) ! BIRTH RATES - HAVE THE SAME RATE as death - life time 60 years. 
         NU=1./(60.*365.0*one_day) ! DEATH RATES HAVE THE SAME RATE

            KK=-1.0E+10
            SIGMAT=0.0
            SIGMA_S_OFF=0.0
            length=real(nx-2)*dx
            SIGMA_F=0.0 ! This is always true for time depdent problems
!        print *,'HOMES_OCCUPIED,OFFICES,HOSPITAL,SHOPS,SCHOOL:',HOMES_OCCUPIED,OFFICES,HOSPITAL,SHOPS,SCHOOL
!        print *,'rday=',rday
        rhome_max=0.
        rbuilding_max=0.0
            DO K=2,NZ-1
            DO J=2,NY-1
            DO I=2,NX-1
               IS_BUILDING=(MATERIAL(I,J,K)==HOMES_OCCUPIED) &
                     .OR.(MATERIAL(I,J,K)==OFFICES) &
                     .OR.(MATERIAL(I,J,K)==HOSPITAL) &
                     .OR.(MATERIAL(I,J,K)==SHOPS) &
                     .OR.(MATERIAL(I,J,K)==SCHOOL) 
 !            print *,'i,j,MATERIAL(I,J,K):',i,j,MATERIAL(I,J,K)

               IF(MATERIAL(I,J,K)==HOMES_OCCUPIED) THEN ! bottom middle region
                  VIRUS_N1_AIM=1000.0 * (1.0-RDAY) + 1000.0
                  rhome_max=max(rhome_max,VIRUS_N1_AIM)
                  VIRUS_N2_AIM=0.0
                  LAMBDA_H_H= 1000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.1day. 
               ELSE IF(IS_BUILDING) THEN
                  VIRUS_N1_AIM=1000.0 * RDAY + 0.0
                  rbuilding_max=max(rbuilding_max,VIRUS_N1_AIM)
                  VIRUS_N2_AIM=0.0
                  LAMBDA_H_H= 1000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.1day. 
               ELSE 
                  VIRUS_N1_AIM=0.0
                  VIRUS_N2_AIM=0.0
                  LAMBDA_H_H= 0.0/one_day ! Push some people out of their homes on time scale of 0.1day. 
               ENDIF
               if(EIGEN_VALUE) LAMBDA_H_H= 1000.0*1.0/one_day


               IF((ITIME==1).AND.(.NOT.EIGEN_VALUE)) THEN
 !                 IF(MATERIAL(I,J,K)==HOMES_OCCUPIED) THEN
                  IF(IS_BUILDING) THEN
!                     T_OLD(I,J,K,1)=(1000.0 * (1.0-RDAY) + 1000.0)*0.999 * 2./1.5
!                     T_OLD(I,J,K,2)=(1000.0 * (1.0-RDAY) + 1000.0)*0.001 * 2./1.5! start off with 0.1% of people exposed and at home.
                     T_OLD(I,J,K,1)=VIRUS_N1_AIM*0.999 * 2./1.5
                     T_OLD(I,J,K,2)=VIRUS_N1_AIM*0.001 * 2./1.5! start off with 0.1% of people exposed and at home.
                  ENDIF
                  IF(ITS==1) T_NEW(I,J,K,:)=T_OLD(I,J,K,:) 
               ENDIF

               VIRUS_N1 = MAX(TOLER,SUM( T_NEW(I,J,K,1:4) ) )
               VIRUS_N2 = MAX(TOLER,SUM( T_NEW(I,J,K,5:8) ) )

!               print *,'i,j,k,MATERIAL(I,J,K):',i,j,k,MATERIAL(I,J,K),T_NEW(I,J,K,:)

               IF( IS_BUILDING &
               .OR.(MATERIAL(I,J,K)==ROADS) &
               .OR.(MATERIAL(I,J,K)==CROSS_ROADS) ) THEN ! CROSS REGION
!                  KK(I,J,K,:)=0.001*2.5*length**2/one_day
                  KK(I,J,K,:)=2.5*length**2/one_day
!                  KK(I,J,K,:)= max(KK(I,J,K,:) ,  (1.-rday**2)*10.0*length**2/one_day )
!                  IF(MATERIAL(I,J,K)==5) KK(I,J,K,:)=0.001*length**2/one_day
               ELSE IF(MATERIAL(I,J,K)==PARK) THEN
!                  KK(I,J,K,:)=0.1*2.5*length**2/one_day
                  KK(I,J,K,:)=2.5*length**2/one_day
               ELSE
                  KK(I,J,K,:)=-1.0E+10
               ENDIF
               KK(I,J,K,1:4)=0.0 ! dont have diffusion in houses, shops, schools,hospitals, offcies. 
               BETA=VIRUS_GAMMA*R0_HOME
               BETA2=VIRUS_GAMMA*R0_MOBILE
!               if(MATERIAL(I,J,K)==-5) BETA2=VIRUS_GAMMA*R0_MOBILE*1.0 ! area on right has good social distancing. 
!               if(MATERIAL(I,J,K)==-HOMES) then
!                  BETA=0.0 ! area on right has good social distancing. 
!                  BETA2=0.0 ! area on right has good social distancing.
!                  if(EIGEN_VALUE) LAMBDA_H_H= 100.0*1.0/one_day 
!               endif
               LAMBDA_M_M=0.0
               if(EIGEN_VALUE) then
                  if(IS_BUILDING) LAMBDA_M_M= 100.0*1.0/one_day 
               endif
!               if(MATERIAL(I,J,K)==5) BETA2=VIRUS_GAMMA*R0_HOME ! area on right has good social distancing.
!               if(MATERIAL(I,J,K)==9) BETA2=VIRUS_GAMMA*R0_HOME ! area on right has good social distancing. 
!               K(I,J,K,:)=0.0
!               K(I,J,K,:)=0.1*RDAY*1.0/(one_day*length*
!               IF(MATERIAL(I,J,K)==2) THEN
                  FORCE=(VIRUS_N1-VIRUS_N1_AIM)/MAX( TOLER, VIRUS_N1, VIRUS_N1_AIM ) 
!                  FORCE2=(VIRUS_N2-VIRUS_N2_AIM(I,J,K))/MAX(TOLER,VIRUS_N2)
                  H2M=0.5+0.5*SIGN(1.0, FORCE) 

                  SIGMAT(I,J,K,1)=-MU          + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M *0.01
                  SIGMAT(I,J,K,2)=+VIRUS_SIGMA + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMAT(I,J,K,3)=+VIRUS_GAMMA + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMAT(I,J,K,4)=+XI          + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M *0.01

                  SIGMAT(I,J,K,5)=-MU          + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMAT(I,J,K,6)=+VIRUS_SIGMA + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMAT(I,J,K,7)=+VIRUS_GAMMA + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMAT(I,J,K,8)=+XI          + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
! Mu terms...
                  SIGMA_S_OFF(I,J,K,1,2)= SIGMA_S_OFF(I,J,K,1,2) - MU 
                  SIGMA_S_OFF(I,J,K,1,3)= SIGMA_S_OFF(I,J,K,1,3) - MU 
                  SIGMA_S_OFF(I,J,K,1,4)= SIGMA_S_OFF(I,J,K,1,4) - MU 

                  SIGMA_S_OFF(I,J,K,5,6)= SIGMA_S_OFF(I,J,K,5,6) - MU 
                  SIGMA_S_OFF(I,J,K,5,7)= SIGMA_S_OFF(I,J,K,5,7) - MU 
                  SIGMA_S_OFF(I,J,K,5,8)= SIGMA_S_OFF(I,J,K,5,8) - MU 


! Time dep...
                  SIGMA_S_OFF(I,J,K,1,5)= SIGMA_S_OFF(I,J,K,1,5) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMA_S_OFF(I,J,K,2,6)= SIGMA_S_OFF(I,J,K,2,6) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMA_S_OFF(I,J,K,3,7)= SIGMA_S_OFF(I,J,K,3,7) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMA_S_OFF(I,J,K,4,8)= SIGMA_S_OFF(I,J,K,4,8) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)

                  SIGMA_S_OFF(I,J,K,5,1)= SIGMA_S_OFF(I,J,K,5,1) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMA_S_OFF(I,J,K,6,2)= SIGMA_S_OFF(I,J,K,6,2) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMA_S_OFF(I,J,K,7,3)= SIGMA_S_OFF(I,J,K,7,3) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMA_S_OFF(I,J,K,8,4)= SIGMA_S_OFF(I,J,K,8,4) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01

                  SIGMA_S_OFF(I,J,K,1,3)=SIGMA_S_OFF(I,J,K,1,3) + (1.-R_EIGEN)*BETA *T_NEW(I,J,K,1)/MAX(VIRUS_N1,TOLER) 
                  SIGMA_S_OFF(I,J,K,2,3)=SIGMA_S_OFF(I,J,K,2,3) - (1.-R_EIGEN)*BETA *T_NEW(I,J,K,1)/MAX(VIRUS_N1,TOLER) 
                  SIGMA_S_OFF(I,J,K,5,7)=SIGMA_S_OFF(I,J,K,5,7) + (1.-R_EIGEN)*BETA2*T_NEW(I,J,K,5)/MAX(VIRUS_N2,TOLER)
                  SIGMA_S_OFF(I,J,K,6,7)=SIGMA_S_OFF(I,J,K,6,7) - (1.-R_EIGEN)*BETA2*T_NEW(I,J,K,5)/MAX(VIRUS_N2,TOLER) 
                  SIGMA_S_OFF(I,J,K,1,4) = SIGMA_S_OFF(I,J,K,1,4) - XI
                  SIGMA_S_OFF(I,J,K,5,8) = SIGMA_S_OFF(I,J,K,5,8) - XI

! exposed people resulting in infections...
                  SIGMA_S_OFF(I,J,K,3,2)=SIGMA_S_OFF(I,J,K,3,2) -VIRUS_SIGMA *(1.-R_EIGEN)
                  SIGMA_S_OFF(I,J,K,7,6)=SIGMA_S_OFF(I,J,K,7,6) -VIRUS_SIGMA *(1.-R_EIGEN)
!          print *,'virus_sigma:',virus_sigma
!         stop 393
! infected people that are recovered...
                  SIGMA_S_OFF(I,J,K,4,3)=SIGMA_S_OFF(I,J,K,4,3) -VIRUS_GAMMA 
                  SIGMA_S_OFF(I,J,K,8,7)=SIGMA_S_OFF(I,J,K,8,7) -VIRUS_GAMMA 
          if(EIGEN_VALUE) then
! eigen-value...
! consistent with the houses...
                  KK(I,J,K,1:8)=1.0*0.05 * KK(I,J,K,1:8)
                  r_ratio=25.65
                  rswitch=0.0
                  if(IS_BUILDING) rswitch=1.0 ! switch on where there are houses
!                  if(material(i,j,k)==3) rswitch=1.0 ! switch on where there are houses
                  LAMBDA_H_H=1.0*rswitch*1.0/one_day
                  LAMBDA_M_M=1.0*(1.0-rswitch)*10000.0/one_day

                  SIGMAT(I,J,K,1)= SIGMAT(I,J,K,1) + 10000000000000.0/one_day
                  SIGMAT(I,J,K,5)= SIGMAT(I,J,K,5) + 10000000000000.0/one_day
                  SIGMAT(I,J,K,4)= SIGMAT(I,J,K,4) + 10000000000000.0/one_day
                  SIGMAT(I,J,K,8)= SIGMAT(I,J,K,8) + 10000000000000.0/one_day

                  SIGMAT(I,J,K,1)= SIGMAT(I,J,K,1) + LAMBDA_H_H    + LAMBDA_M_M 
                  SIGMAT(I,J,K,2)= SIGMAT(I,J,K,2) + LAMBDA_H_H    + LAMBDA_M_M 
                  SIGMAT(I,J,K,3)= SIGMAT(I,J,K,3) + LAMBDA_H_H    + LAMBDA_M_M 
                  SIGMAT(I,J,K,4)= SIGMAT(I,J,K,4) + LAMBDA_H_H    + LAMBDA_M_M 

                  SIGMAT(I,J,K,5)= SIGMAT(I,J,K,5) + LAMBDA_H_H * r_ratio
                  SIGMAT(I,J,K,6)= SIGMAT(I,J,K,6) + LAMBDA_H_H * r_ratio
                  SIGMAT(I,J,K,7)= SIGMAT(I,J,K,7) + LAMBDA_H_H * r_ratio
                  SIGMAT(I,J,K,8)= SIGMAT(I,J,K,8) + LAMBDA_H_H * r_ratio

                  SIGMA_S_OFF(I,J,K,1,5)= SIGMA_S_OFF(I,J,K,1,5) - LAMBDA_H_H * r_ratio
                  SIGMA_S_OFF(I,J,K,2,6)= SIGMA_S_OFF(I,J,K,2,6) - LAMBDA_H_H * r_ratio
                  SIGMA_S_OFF(I,J,K,3,7)= SIGMA_S_OFF(I,J,K,3,7) - LAMBDA_H_H * r_ratio
                  SIGMA_S_OFF(I,J,K,4,8)= SIGMA_S_OFF(I,J,K,4,8) - LAMBDA_H_H * r_ratio

                  SIGMA_S_OFF(I,J,K,5,1)= SIGMA_S_OFF(I,J,K,5,1) - LAMBDA_H_H 
                  SIGMA_S_OFF(I,J,K,6,2)= SIGMA_S_OFF(I,J,K,6,2) - LAMBDA_H_H 
                  SIGMA_S_OFF(I,J,K,7,3)= SIGMA_S_OFF(I,J,K,7,3) - LAMBDA_H_H 
                  SIGMA_S_OFF(I,J,K,8,4)= SIGMA_S_OFF(I,J,K,8,4) - LAMBDA_H_H 
                  
! 
               if(EIGEN_METHOD1) then ! method that put the 1/keff over the I - infection eqn

                  SIGMA_F(I,J,K,3,2)=+VIRUS_SIGMA *R_EIGEN
                  SIGMA_F(I,J,K,7,6)=+VIRUS_SIGMA *R_EIGEN

                  SIGMA_S_OFF(I,J,K,1,3)=SIGMA_S_OFF(I,J,K,1,3)+R_EIGEN*BETA
                  SIGMA_S_OFF(I,J,K,5,7)=SIGMA_S_OFF(I,J,K,5,7)+R_EIGEN*BETA2
                  SIGMA_S_OFF(I,J,K,2,3)=SIGMA_S_OFF(I,J,K,2,3)-R_EIGEN*BETA
                  SIGMA_S_OFF(I,J,K,6,7)=SIGMA_S_OFF(I,J,K,6,7)-R_EIGEN*BETA2

               else
                  SIGMA_S_OFF(I,J,K,3,2)=SIGMA_S_OFF(I,J,K,3,2)-VIRUS_SIGMA *R_EIGEN
                  SIGMA_S_OFF(I,J,K,7,6)=SIGMA_S_OFF(I,J,K,7,6)-VIRUS_SIGMA *R_EIGEN

                  SIGMA_S_OFF(I,J,K,1,3)=SIGMA_S_OFF(I,J,K,1,3)+R_EIGEN*BETA
                  SIGMA_S_OFF(I,J,K,5,7)=SIGMA_S_OFF(I,J,K,5,7)+R_EIGEN*BETA2
                  SIGMA_F(I,J,K,2,3)=+R_EIGEN*BETA 
                  SIGMA_F(I,J,K,6,7)=+R_EIGEN*BETA2 
               endif
           endif ! if(EIGEN_VALUE) then
            END DO
            END DO
            END DO
!               stop 292
         END SUBROUTINE DETERMINE_X_SECTIONS_NG8



         SUBROUTINE DETERMINE_X_SECTIONS_NG260( KK,  SIGMAT,  SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, &
                    SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, DX,DY,DZ,DT, &
                    OPTION,ITIME,ITS,ACCTIM,TIME_MAT, RDAY,EIGEN_VALUE,EIGEN_METHOD1, R_EIGEN, &
                    HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
                    HOSPITAL, SCHOOL, OFFICES, SHOPS,   T_OLD,T_NEW, MATERIAL, S,U, NDIM_VEL, rd_mark_x, rd_mark_y, &
                    RECORD_S, SPAR_RECORD_S, RECORD_F, SPAR_RECORD_F, &
                    R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, R0_SHOPS, R0_VEHC_1TO5, R0_PARK ) 
! Determine KK,  SIGMAT,  SIGMA_S_OFF,  SIGMA_F...
! We do this for a town. 
         IMPLICIT NONE
         INTEGER, intent( in ) :: NCOL_SIGMA_S,NCOL_SIGMA_F, NX,NY,NZ, NG,NG2, NDIM_VEL
         INTEGER, intent( in ) :: FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S),FIN_SIGMA_F(NG+1), COL_SIGMA_F(NCOL_SIGMA_F)
         REAL, intent( in ) :: DX,DY,DZ,DT
         REAL, intent( inout ) :: SIGMAT(NX,NY,NZ,NG),KK(NX,NY,NZ,NG),SIGMA_F(NX,NY,NZ,NG,NG2),SIGMA_F_PTER(NX,NY,NZ,NCOL_SIGMA_F)
         REAL, intent( inout ) :: SIGMA_S_OFF(NX,NY,NZ,NG,NG2),  SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S)
         REAL, intent( inout ) :: T_OLD(NX,NY,NZ,NG),T_NEW(NX,NY,NZ,NG)
         REAL, intent( in ) :: ACCTIM,TIME_MAT, RDAY, R_EIGEN 
         REAL, intent( inout ) :: S(NX,NY,NZ,NG),U(NDIM_VEL,NX+1,NY+1,NZ+1,NG)
         LOGICAL, intent( in ) :: EIGEN_VALUE, EIGEN_METHOD1
         INTEGER, intent( in ) :: HOMES, ROADS, CROSS_ROADS, HOMES_OCCUPIED, PARK, &
               HOSPITAL, SCHOOL, OFFICES, SHOPS
         INTEGER, intent( in ) :: ITIME,ITS, OPTION
         INTEGER, intent( in ) :: MATERIAL(NX,NY,NZ)
         INTEGER, intent( in ) :: rd_mark_x(NX), rd_mark_y(NY)
         LOGICAL, intent( in ) :: RECORD_S, RECORD_F
         LOGICAL, intent( inout ) :: SPAR_RECORD_S(NG,NG), SPAR_RECORD_F(NG,NG)
         REAL, intent( in ) :: R0_HOME, R0_OFFICE, R0_SCHOOL, R0_HOSPITAL, R0_PEDEST, R0_DIFF1, R0_DIFF2, &
                               R0_SHOPS, R0_VEHC_1TO5, R0_PARK 
!         R0_HOME=0.2 !5.0 !0.2 ! No of people an infected person infects.
!         R0_OFFICE=5.0 !5.0 !0.2 ! No of people an infected person infects.
!         R0_SCHOOL=2.0 !5.0 !0.2 ! No of people an infected person infects.
!         R0_HOSPITAL=5.0
!         R0_PEDEST=1.5 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
!         R0_DIFF1 =2.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
!         R0_DIFF2 =2.0 !Portal R0 to hospital
!         R0_PARK=1.0 
!         R0_SHOPS =5.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
!         R0_VEHC_1TO5 =2.0 !10.0 !20.0 !5.0 ! No of people an infected person infects.
!local varibales...
         REAL TOLER, INFINY, PIE, R0_MULTIPLICATION
         LOGICAL SWITCH_ON_ROAD_CAR_WALK ! If SWITCH_ON_ROAD_CAR_WALK switch on road use by cars and walking
! beta only on works. When everything else is on it does not work. 
         LOGICAL UP_WIND_CARS_WALK
! 
         PARAMETER(TOLER=1.E-15, INFINY=1.E+15, PIE=3.141592654, SWITCH_ON_ROAD_CAR_WALK=.true.) 
!         PARAMETER(R0_MULTIPLICATION=4.0) ! multiply all values of R0 by this factor (=1 may be realistic, =4 gives big response)
         PARAMETER(R0_MULTIPLICATION=1.0) ! multiply all values of R0 by this factor (=1 may be realistic, =4 gives big response)
         PARAMETER(UP_WIND_CARS_WALK=.TRUE.) ! it seems better to biase the traffic based on upwind diffusion information. 
         LOGICAL SWITCH_ON_HOSPITAL ! Switch on hospital
         PARAMETER( SWITCH_ON_HOSPITAL=.TRUE.) 
         LOGICAL SWITCH_ON_AGEING ! Switch on aging
         PARAMETER( SWITCH_ON_AGEING=.TRUE.) 
         LOGICAL SWITCH_ON_BETA ! switch on beta term...
         PARAMETER( SWITCH_ON_BETA=.TRUE.)
         LOGICAL SWITCH_ON_DECAY_RATES ! swtich on the various decay rates...
         PARAMETER(SWITCH_ON_DECAY_RATES=.TRUE.)  
         LOGICAL SWITCH_ON_DIFFUSION ! switch on diffusion...
         PARAMETER(SWITCH_ON_DIFFUSION=.TRUE.) 
         LOGICAL SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS ! allow cars to travel through buildings. 
         PARAMETER(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS=.true.) 
         LOGICAL SWITCH_ON_ROAD_TO_BUILDINGS ! move people in cars and walking to buildings and visa versa.
         PARAMETER(SWITCH_ON_ROAD_TO_BUILDINGS=.true.) 
         REAL, ALLOCATABLE :: R0(:), VIRUS_N(:), VIRUS_N_H_AIM(:), VIRUS_PLACE_N(:)
         REAL, ALLOCATABLE :: VIRUS_SIGMA1(:),VIRUS_SIGMA2(:), VIRUS_GAMMA1(:),VIRUS_GAMMA2(:)
         REAL, ALLOCATABLE :: BETA1(:,:,:), BETA2(:,:,:), COME_AGE_SPECTRUM(:)
         REAL, ALLOCATABLE :: FIXED_N_NORMALIZE(:), N_NORMALIZE(:)
         LOGICAL, ALLOCATABLE :: IS_BIG_ROAD_KEEP(:,:,:)
         REAL DURATION_INFECT, DURATION_INFECT_CHILD, DURATION_INFECT_V_ILL
         REAL INCUBATION_DURATION, INCUBATION_DURATION_V_ILL, INCUBATION_DURATION_V_ILL_CHILD
         REAL XI, MU,NU, NU_I2 
         REAL LENGTH
         REAL ONE_DAY, LAMBDA_H_H, FORCE, H2M
         REAL LAMBDA_RD_WALK_EXCHANGE(0:1), LAMBDA_D_RD_WALK(0:1), LAMBDA_D_RD_WALK_DROP(0:1)
         REAL FRACTION_RD_WALK(0:1), LAMBDA_ILL_PORTAL
         REAL LAMBDA_PORTAL_HOSPITAL, LAMBDA_DISCHARG_HOSPITAL
         REAL DIFF_VERT1,DIFF_VERT2, DIFF_HORI1,DIFF_HORI2, SIG_VALUE1, VALUE1,VALUE2
         REAL R_EIGEN1
!         REAL T_NEW_LIMIT_RD_WALK(0:1)

         INTEGER I,J,K, IDAY_WEEK, IDAY, IWEEK, NSEIIR, NPLACE, NPEOPLE, NPLACE_NPEOPLE
         INTEGER ISIR, IPLACE, IPEOPLE, IPLACE_PEOPLE, IG, JG, IG_H, IG_DIFF1
         INTEGER IPLACE_DIFF1, IPLACE_H, IPLACE_PEOPLE_H, IPLACE_PEOPLE_DIFF1 
         INTEGER JPEOPLE, JPLACE, ICHILD, IDISP, ISIR_CHILD, JSIR, JPLACE_PEOPLE
         INTEGER IG_CHILD, IVEH, IPLEVEL, IPLACE_GO_TO, IG_GO_TO, IPLACE_PEOPLE_GO_TO
         INTEGER JPLEVEL, IPLACE_PORTAL, IG_PORTAL, IPLACE_PEOPLE_PORTAL, IIG, ILOOP
         INTEGER IPLACE_PEOPLE_RD, IPLACE_RD, IG_RD

         LOGICAL IS_BUILDING, IS_HOMES_OCCUPIED, IS_OFFICES, IS_HOSPITAL, IS_SHOPS, IS_SCHOOL
         LOGICAL IS_ROADS, IS_CROSS_ROADS, IS_PARK, IS_BETWEEN_HOMES, IS_BIG_ROAD, CROSS
         LOGICAL TOP, BOTTOM, SWAP
         LOGICAL VERTICAL_ROAD, HORIZONTAL_ROAD, SWITCH_ON_CONSTRAIN_PEOPLE_RATIO

         REAL CAR_SPEED, WALK_SPEED
         REAL H_NO_WORKING, H_NO_NOT_WORKING, H_NO_CHILDREN, H_NO_SHIELDING, AREA_H2OFFICE
         REAL OFFICE_NO_WORKING, AREA_H2SCHOOL, SCHOOL_NO_WORKING, SCHOOL_NO_CHILDREN, AREA_H2SHOP 
         REAL SHOP_NO_WORKING, SHOP_NO_WORKING_WEEKEND, SHOP_NO_NOT_WORKING, SHOP_NO_NOT_WORKING_WEEKEND, SHOP_NO_CHILDREN
         REAL ALPHA, RWEEK_END, VALUE12
         real AM7, AM9, PM5, PM9,   AM8,   PM330, PM4,   PM11,   PM530, PM7, PM6, PM8
         real RDAY_OFFICE,RDAY_SCHOOL,RDAY_SHOP_WORKING,RDAY_SHOP_NOT_WORKING,RDAY_SHOP_CHILDREN
         real one_hour, RDAY_SIN, RSIN, RDAY_CHILDREN, SUM_FIXED_N_NORMALIZE, LAMBDA_C_C
         REAL AREA_HOMES, AREA_OFFICE, AREA_HOSPITAL, AREA_SHOPS, AREA_SCHOOL, FRAC_GO_WRONG_WAY
         REAL LAMBDA_R_R
         real VALUE12_DROP, VALUE12_PICKUP
         REAL MIN_VALUE1_W, MIN_VALUE2_W, MAX_VALUE1_W, MAX_VALUE2_W
         REAL MIN_VALUE1_NW, MIN_VALUE2_NW, MAX_VALUE1_NW, MAX_VALUE2_NW, TI

         IF(RECORD_S) SPAR_RECORD_S(:,:)=.FALSE.
         IF(RECORD_F) SPAR_RECORD_F(:,:)=.FALSE.

         R_EIGEN1=0.0
         IF(EIGEN_VALUE .AND. EIGEN_METHOD1) R_EIGEN1=1.0
!         SWITCH_ON_CONSTRAIN_PEOPLE_RATIO=EIGEN_VALUE ! constrain the ratio of people
         SWITCH_ON_CONSTRAIN_PEOPLE_RATIO=.false. ! constrain the ratio of people
!
         one_day=24.*3600.0 
         one_hour=3600.0 
         AM7 = 7.*one_hour; AM9=9.*one_hour;  PM5=17.*one_hour; PM9=21.*one_hour
         AM8 = 8.*one_hour; PM330=15.5*one_hour; PM4=16.*one_hour
         PM11= 23.*one_hour
         PM530=17.5*one_hour; PM7=19.*one_hour
         PM6=18.0*one_hour; PM8=20.0*one_hour

         RSIN = SIN( (TIME_MAT/ONE_DAY) * PIE * 2.0  - 0.5*PIE) 
!         RSIN2 = SIGN( SQRT( ABS(RSIN) ), RSIN) 
         RDAY_SIN = 0.5 * RSIN + 0.5
! 
         NSEIIR=5 
         NPLACE=13
         NPEOPLE=4
         NPLACE_NPEOPLE = NPLACE*NPEOPLE
         ALLOCATE(R0(NPLACE), VIRUS_SIGMA1(NPEOPLE), VIRUS_SIGMA2(NPEOPLE) ) 
         ALLOCATE(VIRUS_GAMMA1(NPEOPLE), VIRUS_GAMMA2(NPEOPLE))  
! Key variables...
         DURATION_INFECT=7.0*one_day ! 7 day infection duration. 
         DURATION_INFECT_CHILD=2.0*one_day ! 2.0 day infection duration for children. 
         DURATION_INFECT_V_ILL=8.0*one_day ! 21 day infection duration. 
         INCUBATION_DURATION=4.5*one_day! 4.5 day inCUBATION duration.
         INCUBATION_DURATION_V_ILL=5.0*one_day! 10.0 day further inCUBATION duration of the very ill. .
         INCUBATION_DURATION_V_ILL_CHILD=100.0*one_day! 100.0 day further inCUBATION duration of the very ill children. .

         ICHILD=3

         VIRUS_SIGMA1(1:NPEOPLE)=1.0/INCUBATION_DURATION
         VIRUS_SIGMA2(1:NPEOPLE)=1.0/INCUBATION_DURATION_V_ILL
         VIRUS_SIGMA2(ICHILD)=1.0/INCUBATION_DURATION_V_ILL_CHILD

! ξ (xi) is the rate which recovered individuals return to the susceptible statue due to loss of immunity.
         XI= 1./(10.0*365.0*one_day) ! 10 year time scale

!         R0_HOME=0.2 !5.0 !0.2 ! No of people an infected person infects.
!         R0_OFFICE=5.0 !5.0 !0.2 ! No of people an infected person infects.
!         R0_SCHOOL=2.0 !5.0 !0.2 ! No of people an infected person infects.
!         R0_HOSPITAL=5.0
!         R0_PEDEST=1.5 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
!         R0_DIFF1 =2.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
!         R0_DIFF2 =2.0 !Portal R0 to hospital
!         R0_PARK=1.0 
!         R0_SHOPS =5.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
!         R0_VEHC_1TO5 =2.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 

         MU=1.0*1./(60.*365.0*one_day) ! BIRTH RATES - HAVE THE SAME RATE as death - life time 60 years. 
         NU=1./(60.*365.0*one_day) ! DEATH RATES HAVE THE SAME RATE
         NU_I2=1./(45.0*one_day) ! DEATH RATES IN the very ill group
! Time scale on which to turn off if one is walking or driving...
!         LAMBDA_RD_WALK_EXCHANGE(0)=0.01*1.0/5.0 ! Inverse time in seconds for driving to turn onto another road. 
!         LAMBDA_RD_WALK_EXCHANGE(0)=0.0/20.0 ! Inverse time in seconds for driving to turn onto another road. 
!         LAMBDA_RD_WALK_EXCHANGE(1)=0.01*1.0/20.0 ! Inverse time in seconds for walking to turn onto another road. 
! 
       if(.true.) then
! factor of 100 here means turning on a time scale of 10 seconds. 
         LAMBDA_RD_WALK_EXCHANGE(0)=1.0/10.0 ! Inverse time in seconds for driving to turn onto another road. 
!         LAMBDA_RD_WALK_EXCHANGE(0)=0.0/20.0 ! Inverse time in seconds for driving to turn onto another road. 
         LAMBDA_RD_WALK_EXCHANGE(1)=1.0/30.0  ! Inverse time in seconds for walking to turn onto another road. 
       else
! factor of 100 here means turning on a time scale of 10 seconds. 
         LAMBDA_RD_WALK_EXCHANGE(0)=100.*min(1.0/5.0, 1.0/1000.0) ! Inverse time in seconds for driving to turn onto another road. 
!         LAMBDA_RD_WALK_EXCHANGE(0)=0.0/20.0 ! Inverse time in seconds for driving to turn onto another road. 
         LAMBDA_RD_WALK_EXCHANGE(1)=100.*min(1.0/20.0, 1.0/1000.0)  ! Inverse time in seconds for walking to turn onto another road. 
       endif
         FRAC_GO_WRONG_WAY=0.2!2 ! FRAC_GO_WRONG_WAY=0.0 -no one goes wrong way on roads, =0.1 then 10% go wrong way (both values are good)
      if(.false.) then ! original...
         LAMBDA_D_RD_WALK(0) =100.*0.01*1000.0/one_day
         LAMBDA_D_RD_WALK(1) =100.*0.01*200.0/one_day
!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
         LAMBDA_D_RD_WALK_DROP(0) =1.*0.1*1000.0/one_day ! time scale for dropping people off
         LAMBDA_D_RD_WALK_DROP(1) =1.*0.1*200.0/one_day
!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK_DROP=0.0 ! Switch off pedestrians and cars
         FRACTION_RD_WALK(0)=1.0 ! 90% of people are in cars or walking
         FRACTION_RD_WALK(1)=1.0 ! 90% of people are in cars or walking
      else
! this works with a fifth in cars 8 and walking 2e+2 and diffusion 80.0 (in g260g-0-1- files):
!         LAMBDA_D_RD_WALK(0) =1.0*0.01*1000.0/one_day ! 86 second invers time scale
!         LAMBDA_D_RD_WALK(1) =100.0*0.01*200.0/one_day ! 5*86 second invers time scale
!!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
!         LAMBDA_D_RD_WALK_DROP(0) =1.0*0.01*1000.0/one_day ! time scale for dropping people off
!         LAMBDA_D_RD_WALK_DROP(1) =100.0*0.01*200.0/one_day
! this results in cars 0.9 and walking 2e+2 and diffusion 8e+1 (in g260g-0-1-- files):
!         LAMBDA_D_RD_WALK(0) =0.1*0.01*1000.0/one_day ! 86 second invers time scale
!         LAMBDA_D_RD_WALK(1) =100.0*0.01*200.0/one_day ! 5*86 second invers time scale
!!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
!         LAMBDA_D_RD_WALK_DROP(0) =0.1*0.01*1000.0/one_day ! time scale for dropping people off
!         LAMBDA_D_RD_WALK_DROP(1) =100.0*0.01*200.0/one_day
! this results in cars 2e+1 and walking 5e+0 and diffusion 2e+2  and s1 4e+2(in g260g-0-1--- files):
!         LAMBDA_D_RD_WALK(0) =1.0*0.01*1000.0/one_day ! 86 second invers time scale
!         LAMBDA_D_RD_WALK(1) =1.0*0.01*200.0/one_day ! 5*86 second invers time scale
!!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
!         LAMBDA_D_RD_WALK_DROP(0) =1.0*0.01*1000.0/one_day ! time scale for dropping people off
!         LAMBDA_D_RD_WALK_DROP(1) =1.0*0.01*200.0/one_day
! this results in cars 1e+2 and walking 2e+1 and diffusion 1e+2  and s1 2e+2(in g260g-0-1---- files):
!         LAMBDA_D_RD_WALK(0) =10.0*0.01*1000.0/one_day ! 86 second invers time scale
!         LAMBDA_D_RD_WALK(1) =10.0*0.01*200.0/one_day ! 5*86 second invers time scale
!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
!         LAMBDA_D_RD_WALK_DROP(0) =1.0*0.01*1000.0/one_day ! time scale for dropping people off
!         LAMBDA_D_RD_WALK_DROP(1) =1.0*0.01*200.0/one_day
! this results in cars 1e+2 and walking 2e+1 and diffusion 1e+2  and s1 2e+2(in g260g-0-1---- files):
!         LAMBDA_D_RD_WALK(0) =10.0*0.01*1000.0/one_day ! 86 second invers time scale
!         LAMBDA_D_RD_WALK(1) =10.0*0.01*200.0/one_day ! 5*86 second invers time scale
!!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
!         LAMBDA_D_RD_WALK_DROP(0) =100.0*0.01*1000.0/one_day ! time scale for dropping people off
!         LAMBDA_D_RD_WALK_DROP(1) =100.0*0.01*200.0/one_day
! this results in cars 1e+2 and walking 2e+1 and diffusion 1e+2  and s1 2e+2(in g260g-0-1---- files):
!         LAMBDA_D_RD_WALK(0) =10.0*0.01*1000.0/one_day ! 86 second invers time scale
!         LAMBDA_D_RD_WALK(1) =10.0*0.01*200.0/one_day ! 5*86 second invers time scale
!!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
!         LAMBDA_D_RD_WALK_DROP(0) =0.1*0.01*1000.0/one_day ! time scale for dropping people off
!         LAMBDA_D_RD_WALK_DROP(1) =0.1*0.01*200.0/one_day
! this results in cars 1e+2 and walking 2e+1 and diffusion 1e+2  and s1 2e+2(in g260g-0-1---- files)????:
         LAMBDA_D_RD_WALK(0) =10.0*0.01*1000.0/one_day ! 86 second invers time scale
         LAMBDA_D_RD_WALK(1) =10.0*0.01*200.0/one_day ! 5*86 second invers time scale
!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK=0.0 ! Switch off pedestrians and cars
         LAMBDA_D_RD_WALK_DROP(0) =1.0*10.0*0.01*1000.0/one_day ! time scale for dropping people off
         LAMBDA_D_RD_WALK_DROP(1) =1.0*10.0*0.01*200.0/one_day
!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) LAMBDA_D_RD_WALK_DROP=0.0 ! Switch off pedestrians and cars
         FRACTION_RD_WALK(0)=1.0 ! 90% of people are in cars or walking
         FRACTION_RD_WALK(1)=1.0 ! 90% of people are in cars or walking
      endif 
!         T_NEW_LIMIT_RD_WALK(0)=10.0 ! Limit of diffusion process for cars. 
!         T_NEW_LIMIT_RD_WALK(1)=10.0 ! Limit of diffusion process for walking
         LAMBDA_ILL_PORTAL=96./(one_day) ! move people to portal 
! on a time scale of an ambulance (say 15mins)... LAMBDA_ILL_PORTAL 
         LAMBDA_PORTAL_HOSPITAL=10000.0 * 96./(one_day) ! anyone in portal put in hospital.
         LAMBDA_DISCHARG_HOSPITAL= 1.0/(0.25*one_day) ! How long it takes for not very ill people to be discharged (0.25 day). 

         IDAY=1+INT(ACCTIM/ONE_DAY); IWEEK=1 + (IDAY-1) / 7; IDAY_WEEK= IDAY- (IWEEK-1)*7
         IF((IDAY_WEEK==6).OR.(IDAY_WEEK==7)) THEN
            RWEEK_END = 1.0! Start Monday every 6th and 7th day a weekend. 
         ELSE
            RWEEK_END = 0.0
         ENDIF

       if(.false.) then
         do i=1,7000
         ti=real(i-1)*1000.
         IDAY=1+INT(TI/ONE_DAY); IWEEK=1 + (IDAY-1) / 7; IDAY_WEEK= IDAY- (IWEEK-1)*7
!         IDAY=1+INT(ACCTIM/ONE_DAY); IWEEK=1 + (IDAY-1) / 7; IDAY_WEEK= IDAY- (IWEEK-1)*7
         IF((IDAY_WEEK==6).OR.(IDAY_WEEK==7)) THEN
            RWEEK_END = 1.0! Start Monday every 6th and 7th day a weekend. 
         ELSE
            RWEEK_END = 0.0
         ENDIF
           print *,ti,rweek_end
         end do
         stop 282
       endif

         KK=-1.0E+10
         SIGMAT=0.0
         if(.not.record_s) then
            if(NCOL_SIGMA_S.NE.0) THEN ! the pointers
!           print *,'NCOL_SIGMA_S, NCOL_SIGMA_F:',NCOL_SIGMA_S, NCOL_SIGMA_F
!           stop 282
              SIGMA_S_OFF_PTER=0.0
              SIGMA_F_PTER=0.0
            else
              SIGMA_S_OFF=0.0
              SIGMA_F=0.0 ! This is always true for time depdent problems
            endif
         endif
         IF(RECORD_S.NEQV.RECORD_F) STOP 3922
         length=real(nx-2)*dx
!        print *,'HOMES_OCCUPIED,OFFICES,HOSPITAL,SHOPS,SCHOOL:',HOMES_OCCUPIED,OFFICES,HOSPITAL,SHOPS,SCHOOL
!        print *,'rday=',rday
! seirr GROUPS:
!          S_G=1; E_G=2; I1_G=3; I2_G=4; R_G=5
! place GROUPS...
!          H_G=1; V1_G=2; V2_G=3; V3_G=4; V4_G=5; VV_G=6; P1_G=7; P2_G=8; P3_G=9; P4_G=10; PP_G=11; D1_G=12; D2_G=13
! people GROUPS...
!          OW_G=1; NOW=2; CH=3; SH=4

         ALLOCATE( VIRUS_N(NPLACE*NPEOPLE), VIRUS_PLACE_N(NPLACE) )
         ALLOCATE( BETA1(NPLACE,NPEOPLE,NPEOPLE), BETA2(NPLACE,NPEOPLE,NPEOPLE) )
         ALLOCATE( COME_AGE_SPECTRUM(NPEOPLE), VIRUS_N_H_AIM(NPEOPLE))

!         print *,'rd_mark_x:',rd_mark_x
!         print *,'rd_mark_y:',rd_mark_y
!         stop 21

! time scale for becomming an adult of a child (16 years)...
         ALPHA=1./(16.0*365.0*one_day) 


         CAR_SPEED= 8.94 ! m/s (20miles pr hour)
!         CAR_SPEED= 4.47 ! m/s( 10miles pr hour)
         WALK_SPEED=1.4 ! m/s
         U=0.0
         S=0.0 ! Source
! Max no of people in homes occuppied...     H place (house, office, school, shops) 
         H_NO_WORKING=1000.0
         H_NO_NOT_WORKING=500.0
         H_NO_CHILDREN=350.0
         H_NO_SHIELDING=350.0
! this is the comming of age spectrum- how we distribute children to adult groups as they come of age. 
         COME_AGE_SPECTRUM(1)=H_NO_WORKING
         COME_AGE_SPECTRUM(2)=H_NO_NOT_WORKING
         COME_AGE_SPECTRUM(3)=0.0
         COME_AGE_SPECTRUM(4)=H_NO_SHIELDING
         COME_AGE_SPECTRUM(:)=COME_AGE_SPECTRUM(:)/SUM(COME_AGE_SPECTRUM(:)) ! spectrum needs to sum to 1.
! 
         
         AREA_HOMES    = 0.0
         AREA_OFFICE   = 0.0
         AREA_HOSPITAL = 0.0
         AREA_SHOPS    = 0.0
         AREA_SCHOOL   = 0.0
         DO K=2,NZ-1
         DO J=2,NY-1
         DO I=2,NX-1
            IS_HOMES_OCCUPIED=(MATERIAL(I,J,K)==HOMES_OCCUPIED) 
            IS_OFFICES=(MATERIAL(I,J,K)==OFFICES) 
            IS_HOSPITAL=(MATERIAL(I,J,K)==HOSPITAL) 
            IS_SHOPS=(MATERIAL(I,J,K)==SHOPS) 
            IS_SCHOOL=(MATERIAL(I,J,K)==SCHOOL) 
! 
            IF(IS_HOMES_OCCUPIED) AREA_HOMES    = AREA_HOMES+1.0
            IF(IS_OFFICES)        AREA_OFFICE   = AREA_OFFICE+1.0
            IF(IS_HOSPITAL)       AREA_HOSPITAL = AREA_HOSPITAL + 1.0
            IF(IS_SHOPS)          AREA_SHOPS    = AREA_SHOPS+1.0
            IF(IS_SCHOOL)         AREA_SCHOOL   = AREA_SCHOOL+1.0
         END DO
         END DO
         END DO
! Max no of people in offices...
!         AREA_H2OFFICE=100.0 
         AREA_H2OFFICE=AREA_HOMES / AREA_OFFICE
!          print *,'AREA_H2OFFICE:',AREA_H2OFFICE
!          stop 221
         OFFICE_NO_WORKING=H_NO_WORKING * AREA_H2OFFICE
! Max no of people working at school...
!         AREA_H2SCHOOL=100.0  ! AREA RATION HOME TO SCHOOL
         AREA_H2SCHOOL=AREA_HOMES / AREA_SCHOOL
         SCHOOL_NO_WORKING=0.01*H_NO_CHILDREN * AREA_H2SCHOOL ! one adult for 20 children
! Max no of children in school...
         SCHOOL_NO_CHILDREN=1.0*H_NO_CHILDREN * AREA_H2SCHOOL
! max no of people in schops of working people during the day...
!         AREA_H2SHOP=20.0 
         AREA_H2SHOP=H_NO_WORKING / AREA_SHOPS
         SHOP_NO_WORKING=0.025*H_NO_WORKING * AREA_H2SHOP
! max no of working people during the weekend visiting shops...
         SHOP_NO_WORKING_WEEKEND=0.1*H_NO_WORKING* AREA_H2SHOP
! max no of not working people visiting shops...
!         SHOP_NO_NOT_WORKING=  0.1 * H_NO_NOT_WORKING * AREA_H2SHOP
         SHOP_NO_NOT_WORKING=  0.3 * H_NO_NOT_WORKING * AREA_H2SHOP
!         print *,'SHOP_NO_NOT_WORKING, H_NO_NOT_WORKING, AREA_H2SHOP:',SHOP_NO_NOT_WORKING, H_NO_NOT_WORKING, AREA_H2SHOP
!            stop 33
! max no of not working people during the weekend visiting shops...
         SHOP_NO_NOT_WORKING_WEEKEND = 0.1 * H_NO_NOT_WORKING * AREA_H2SHOP
! max no of CHILDREN people during the weekend visiting shops..
         SHOP_NO_CHILDREN=0.1 * H_NO_CHILDREN * AREA_H2SHOP


         ALLOCATE( IS_BIG_ROAD_KEEP(NX,NY,NZ) )
         IS_BIG_ROAD_KEEP=.FALSE.
         ALLOCATE( FIXED_N_NORMALIZE(NPEOPLE), N_NORMALIZE(NPEOPLE) )
         FIXED_N_NORMALIZE(:)=0.0
         N_NORMALIZE(:)=0.0
         DO K=2,NZ-1
         DO J=2,NY-1
         DO I=2,NX-1
               
               IS_HOMES_OCCUPIED=(MATERIAL(I,J,K)==HOMES_OCCUPIED) 
               IS_OFFICES=(MATERIAL(I,J,K)==OFFICES) 
               IS_HOSPITAL=(MATERIAL(I,J,K)==HOSPITAL) 
               IS_SHOPS=(MATERIAL(I,J,K)==SHOPS) 
               IS_SCHOOL=(MATERIAL(I,J,K)==SCHOOL) 

               IS_BUILDING=IS_HOMES_OCCUPIED &
                       .OR.IS_OFFICES &
                       .OR.IS_HOSPITAL &
                       .OR.IS_SHOPS &
                       .OR.IS_SCHOOL

               IS_PARK=(MATERIAL(I,J,K)==PARK) 

            IS_ROADS = (MATERIAL(I,J,K)==ROADS) 
!            IS_BETWEEN_HOMES = IS_ROADS .AND. &  
!                 (    (((MATERIAL(I-1,J,K)==HOMES_OCCUPIED).AND.(MATERIAL(I+1,J,K)==HOMES_OCCUPIED))  &
!                   .OR.((MATERIAL(I,J-1,K)==HOMES_OCCUPIED).AND.(MATERIAL(I,J+1,K)==HOMES_OCCUPIED)) ) &
!                .OR. ( ((MATERIAL(I-1,J-1,K)==HOMES_OCCUPIED).AND.(MATERIAL(I+1,J+1,K)==HOMES_OCCUPIED))  &
!                  .AND.((MATERIAL(I-1,J+1,K)==HOMES_OCCUPIED).AND.(MATERIAL(I+1,J-1,K)==HOMES_OCCUPIED)) )   ) 

            VERTICAL_ROAD  =.false.
            HORIZONTAL_ROAD=.false.
            if(  ((rd_mark_x(i-1)==1).and.(rd_mark_x(i)==1))  &
             .or.((rd_mark_x(i+1)==1).and.(rd_mark_x(i)==1)) ) VERTICAL_ROAD=.TRUE.
            if(  ((rd_mark_y(j-1)==1).and.(rd_mark_y(j)==1))  &
             .or.((rd_mark_y(j+1)==1).and.(rd_mark_y(j)==1)) ) HORIZONTAL_ROAD=.TRUE.

!            IS_BIG_ROAD = IS_ROADS .AND. (.NOT.IS_BETWEEN_HOMES) 
!            IS_BIG_ROAD = VERTICAL_ROAD .or. HORIZONTAL_ROAD
!            IS_BIG_ROAD_KEEP(I,J,K) = IS_BIG_ROAD
               IF(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS) THEN
                  IS_BIG_ROAD = (VERTICAL_ROAD .or. HORIZONTAL_ROAD) !.and. (.not. IS_BUILDING) 
                  IS_CROSS_ROADS = (VERTICAL_ROAD .AND. HORIZONTAL_ROAD) !.and. (.not. IS_BUILDING)
               ELSE ! IF(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS) THEN
                  IS_BIG_ROAD = (VERTICAL_ROAD .or. HORIZONTAL_ROAD) .and. (.not. IS_BUILDING) 
                  IS_CROSS_ROADS = (VERTICAL_ROAD .AND. HORIZONTAL_ROAD) .and. (.not. IS_BUILDING)
               ENDIF ! IF(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS) THEN ELSE
               IS_BIG_ROAD_KEEP(I,J,K) = IS_BIG_ROAD
!            if( IS_BIG_ROAD .neqv. (  VERTICAL_ROAD .or. HORIZONTAL_ROAD ) ) then
!               print *,'i,j,IS_BIG_ROAD, VERTICAL_ROAD.or.HORIZONTAL_ROAD,is_roads,xi,yi:',  i,j,IS_BIG_ROAD, VERTICAL_ROAD.or.HORIZONTAL_ROAD,is_roads, rd_mark_x(i), rd_mark_y(j)
!            endif
               IF((ITIME==1).AND.(.NOT.EIGEN_VALUE)) THEN
 !                 IF(MATERIAL(I,J,K)==HOMES_OCCUPIED) THEN
                  IF(IS_HOMES_OCCUPIED) THEN
! INITIAL CONDITIONS...
                     ISIR=1
                     IPLACE=1
                     DO IPEOPLE=1,NPEOPLE
                        IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                        IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                        IF(IPEOPLE==1) THEN ! workers
                           T_OLD(I,J,K, IG)  =H_NO_WORKING*0.99 
                           T_OLD(I,J,K, IG+1)=H_NO_WORKING*0.01 ! start off with 1% of people exposed and at home.
                        ELSE IF(IPEOPLE==2) THEN ! NOT workers
                           T_OLD(I,J,K, IG)  =H_NO_NOT_WORKING
                        ELSE IF(IPEOPLE==3) THEN ! CHILDREN...
                           T_OLD(I,J,K, IG)  =H_NO_CHILDREN
                        ELSE ! SHIELDING ...
                           T_OLD(I,J,K, IG)  =H_NO_SHIELDING
                        ENDIF
                     END DO
                  ENDIF
                  IF(ITS==1) T_NEW(I,J,K,:)=T_OLD(I,J,K,:) 
               ENDIF

! Normaization factors... 
                  IF(IS_HOMES_OCCUPIED) THEN
! Fixed normalization factors...
                     IPEOPLE=1 ! workers
                     FIXED_N_NORMALIZE(IPEOPLE) = FIXED_N_NORMALIZE(IPEOPLE) + H_NO_WORKING
                     IPEOPLE=2 ! NOT workers
                     FIXED_N_NORMALIZE(IPEOPLE) = FIXED_N_NORMALIZE(IPEOPLE) + H_NO_NOT_WORKING
                     IPEOPLE=3 ! CHILDREN...
                     FIXED_N_NORMALIZE(IPEOPLE) = FIXED_N_NORMALIZE(IPEOPLE) + H_NO_CHILDREN
                     IPEOPLE=4 ! SHIELDING ...
                     FIXED_N_NORMALIZE(IPEOPLE) = FIXED_N_NORMALIZE(IPEOPLE) + H_NO_SHIELDING
                  ENDIF

! normalization for this simulation...
                  DO ISIR=1,NSEIIR
                  DO IPLACE=1,NPLACE
                  DO IPEOPLE=1,NPEOPLE
                     IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                     IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                     N_NORMALIZE(IPEOPLE) = N_NORMALIZE(IPEOPLE) + T_NEW(I,J,K,IG)
                  END DO
                  END DO
                  END DO

         END DO ! DO I=2,NX-1
         END DO ! DO J=2,NY-1
         END DO ! DO K=2,NZ-1
!          stop 282

            print *,'N_NORMALIZE:',N_NORMALIZE
            print *,'FIXED_N_NORMALIZE:',FIXED_N_NORMALIZE
            print *,'N_NORMALIZE/FIXED_N_NORMALIZE:',N_NORMALIZE/FIXED_N_NORMALIZE
!         if(abs(N_NORMALIZE(2)/FIXED_N_NORMALIZE(2) -1.0).gt.0.0001)  stop 2817
!         if( (abs(N_NORMALIZE(2)).gt.0.0001) .and. (abs(N_NORMALIZE(2)/FIXED_N_NORMALIZE(2) -1.0).gt.0.0001) )  stop 2817
! 
           MIN_VALUE1_W=1.E+20
           MIN_VALUE2_W=1.E+20
           MAX_VALUE1_W=-1.E+20
           MAX_VALUE2_W=-1.E+20

           MIN_VALUE1_NW=1.E+20
           MIN_VALUE2_NW=1.E+20
           MAX_VALUE1_NW=-1.E+20
           MAX_VALUE2_NW=-1.E+20

            DO K=2,NZ-1
            DO J=2,NY-1
            DO I=2,NX-1
               
               IS_HOMES_OCCUPIED=(MATERIAL(I,J,K)==HOMES_OCCUPIED) 
               IS_OFFICES=(MATERIAL(I,J,K)==OFFICES) 
               IS_HOSPITAL=(MATERIAL(I,J,K)==HOSPITAL) 
               IS_SHOPS=(MATERIAL(I,J,K)==SHOPS) 
               IS_SCHOOL=(MATERIAL(I,J,K)==SCHOOL) 

               IS_BUILDING=IS_HOMES_OCCUPIED &
                       .OR.IS_OFFICES &
                       .OR.IS_HOSPITAL &
                       .OR.IS_SHOPS &
                       .OR.IS_SCHOOL

               IS_ROADS=(MATERIAL(I,J,K)==ROADS) 
               IS_CROSS_ROADS=(MATERIAL(I,J,K)==CROSS_ROADS) 
               IS_PARK=(MATERIAL(I,J,K)==PARK) 
!               IS_BETWEEN_HOMES=IS_ROADS .AND. &  
!                 (    (((MATERIAL(I-1,J,K)==HOMES_OCCUPIED).AND.(MATERIAL(I+1,J,K)==HOMES_OCCUPIED))  &
!                       .OR.((MATERIAL(I,J-1,K)==HOMES_OCCUPIED).AND.(MATERIAL(I,J+1,K)==HOMES_OCCUPIED)) ) &
!                .OR. ( ((MATERIAL(I-1,J-1,K)==HOMES_OCCUPIED).AND.(MATERIAL(I+1,J+1,K)==HOMES_OCCUPIED))  &
!                      .AND.((MATERIAL(I-1,J+1,K)==HOMES_OCCUPIED).AND.(MATERIAL(I+1,J-1,K)==HOMES_OCCUPIED)) )   ) 
!               IS_BIG_ROAD = IS_ROADS .AND. (.NOT.IS_BETWEEN_HOMES) 
 !            print *,'i,j,MATERIAL(I,J,K):',i,j,MATERIAL(I,J,K)

               VERTICAL_ROAD  =.false.
               HORIZONTAL_ROAD=.false.
               if(  ((rd_mark_x(i-1)==1).and.(rd_mark_x(i)==1))  &
                .or.((rd_mark_x(i+1)==1).and.(rd_mark_x(i)==1)) ) VERTICAL_ROAD=.TRUE.
               if(  ((rd_mark_y(j-1)==1).and.(rd_mark_y(j)==1))  &
                .or.((rd_mark_y(j+1)==1).and.(rd_mark_y(j)==1)) ) HORIZONTAL_ROAD=.TRUE.

               IF(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS) THEN
                  IS_BIG_ROAD = (VERTICAL_ROAD .or. HORIZONTAL_ROAD) !.and. (.not. IS_BUILDING) 
                  IS_CROSS_ROADS = (VERTICAL_ROAD .AND. HORIZONTAL_ROAD) !.and. (.not. IS_BUILDING)
                  IF(IS_BIG_ROAD) THEN
                     IS_ROADS=.TRUE.
                  ENDIF
               ELSE ! IF(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS) THEN
                  IS_BIG_ROAD = (VERTICAL_ROAD .or. HORIZONTAL_ROAD) .and. (.not. IS_BUILDING) 
                  IS_CROSS_ROADS = (VERTICAL_ROAD .AND. HORIZONTAL_ROAD) .and. (.not. IS_BUILDING)
                  IF(IS_BIG_ROAD) THEN
                     IS_ROADS=.TRUE.
                     IS_PARK=.FALSE.
! 
                     IS_HOMES_OCCUPIED=.FALSE.
                     IS_OFFICES=.FALSE. 
                     IS_HOSPITAL=.FALSE. 
                     IS_SHOPS=.FALSE.
                     IS_SCHOOL=.FALSE.

                     IS_BUILDING=.FALSE. 
                  ENDIF
               ENDIF ! IF(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS) THEN ELSE


! ***************GAMMA & BETA (BEGIN)******************************
! No of people an infected person infects. 
               VIRUS_GAMMA1(1:NPEOPLE)=1.0/DURATION_INFECT
               VIRUS_GAMMA2(1:NPEOPLE)=1.0/DURATION_INFECT
! 
               R0(1)=R0_HOME
               IF(IS_OFFICES) R0(1)=R0_OFFICE
               IF(IS_SCHOOL) R0(1)=R0_SCHOOL
               IF(IS_SHOPS)  R0(1)=R0_SHOPS
               IF(IS_HOSPITAL) R0(1)=R0_HOSPITAL
               R0(2)=R0_VEHC_1TO5
               R0(3)=R0_VEHC_1TO5
               R0(4)=R0_VEHC_1TO5
               R0(5)=R0_VEHC_1TO5
               R0(6)=R0_VEHC_1TO5
               R0(7)=R0_PEDEST
               R0(8)=R0_PEDEST
               R0(9)=R0_PEDEST
               R0(10)=R0_PEDEST
               R0(11)=R0_PEDEST
               R0(12)=R0_DIFF1
               R0(13)=R0_DIFF2
               IF(IS_PARK) R0(12)=R0_PARK
               R0(:)=R0_MULTIPLICATION * R0(:)
! BETA1(IPLACE,IPEOPLE,JPEOPLE)...
               DO IPLACE=1,NPLACE
                  DO IPEOPLE=1,NPEOPLE
                     BETA1(IPLACE,IPEOPLE,1:NPEOPLE)=VIRUS_GAMMA1(IPEOPLE)*R0(IPLACE)
                     BETA2(IPLACE,IPEOPLE,1:NPEOPLE)=VIRUS_GAMMA2(IPEOPLE)*R0(IPLACE)
                  END DO
               END DO
! redfine gamm with the correct numbers as we did not want to work out another R0:
               ICHILD=3
               VIRUS_GAMMA1(1:NPEOPLE)=1.0/DURATION_INFECT
               VIRUS_GAMMA2(1:NPEOPLE)=1.0/DURATION_INFECT_V_ILL
               VIRUS_GAMMA1(ICHILD)=1.0/DURATION_INFECT_CHILD
! ***************GAMMA & BETA (END)********************************


!               LAMBDA_H_H=0.0*1.0/one_day ! Push some people out of their homes on time scale of 0.01day. 
!               LAMBDA_H_H=100.0*1.0/one_day ! Push some people out of their homes on time scale of 0.01day. 
!               LAMBDA_H_H=1000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.01day. 
               LAMBDA_H_H=2000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.01day. 
!               LAMBDA_H_H=10000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.01day. 
!               LAMBDA_H_H=1000000.0*1.0/one_day ! BEST Push some people out of their homes on time scale of 0.01day. 
!               LAMBDA_H_H=1000000000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.01day. 
!               LAMBDA_H_H=100.*1000000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.01day. 
               IF(IS_HOMES_OCCUPIED) THEN ! bottom middle region
!                  PRINT *,'RDAY_SIN=',RDAY_SIN
!                  STOP 2921
                  CALL DAY_CYCLE(RDAY_OFFICE, TIME_MAT, AM7, AM9,   PM5, PM8) ! 0110
                  RDAY_OFFICE=1.0-RDAY_OFFICE
                  CALL DAY_CYCLE(RDAY_CHILDREN, TIME_MAT, AM8, AM9,   PM4, PM7) ! 0110
                  RDAY_CHILDREN=1.0-RDAY_CHILDREN
                  VIRUS_N_H_AIM(1) =  ( RDAY_OFFICE*(1.-RWEEK_END) + MAX(0.5,RDAY_OFFICE)* RWEEK_END )* H_NO_WORKING
                  VIRUS_N_H_AIM(2) = 0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING
!                  VIRUS_N_H_AIM(2) = 0.1*  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING
!                  VIRUS_N_H_AIM(2) = 0.0000000001*  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING ! this works ok
!                  VIRUS_N_H_AIM(2) = 0.0000001*  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING ! this works ok
!                  VIRUS_N_H_AIM(2) = 0.00001*  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING ! DIDNT work
!                  VIRUS_N_H_AIM(2) = 0.001*  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING ! DIDNT WORK
!                  VIRUS_N_H_AIM(2) = 0.1*  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING ! 
!                  VIRUS_N_H_AIM(2) = 0.000*  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING
!              print *,'VIRUS_N_H_AIM(2):',VIRUS_N_H_AIM(2)
!              stop 281
                  VIRUS_N_H_AIM(3) = ( RDAY_CHILDREN*(1.-RWEEK_END) + MAX(0.5,RDAY_CHILDREN)* RWEEK_END )* H_NO_CHILDREN
!                  VIRUS_N_H_AIM(1) =  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_WORKING
!                  VIRUS_N_H_AIM(2) =  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING
!                  VIRUS_N_H_AIM(3) = 0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_CHILDREN
                  VIRUS_N_H_AIM(4) = ( 0.9 + 0.1* 0.666*((1.0-RDAY_SIN) + 0.5) )  * H_NO_SHIELDING ! 5% of people are out during peak day
               ELSE IF(IS_OFFICES) THEN ! take hospital out of this...
                  CALL DAY_CYCLE(RDAY_OFFICE, TIME_MAT, AM7, AM9,   PM5, PM6) ! 0110
                  VIRUS_N_H_AIM(1) = RDAY_OFFICE * OFFICE_NO_WORKING * (1.0 - RWEEK_END) 
                  VIRUS_N_H_AIM(2) = 0.0
                  VIRUS_N_H_AIM(3) = 0.0
                  VIRUS_N_H_AIM(4) = 0.0
!                  LAMBDA_H_H(2)=0.0
               ELSE IF(IS_SCHOOL) THEN ! take hospital out of this...
                  CALL DAY_CYCLE(RDAY_SCHOOL, TIME_MAT, AM8, AM9,   PM330, PM4) ! 0110
                  VIRUS_N_H_AIM(1) = RDAY_SCHOOL * SCHOOL_NO_WORKING *  (1.0 - RWEEK_END)
                  VIRUS_N_H_AIM(2) = 0.0
                  VIRUS_N_H_AIM(3) = RDAY_SCHOOL * SCHOOL_NO_CHILDREN *  (1.0 - RWEEK_END)
                  VIRUS_N_H_AIM(4) = 0.0
               ELSE IF(IS_SHOPS) THEN ! take hospital out of this...
                  CALL DAY_CYCLE(RDAY_SHOP_WORKING,     TIME_MAT, AM8, AM9,   PM5,   PM8) ! 0110
                  CALL DAY_CYCLE(RDAY_SHOP_NOT_WORKING, TIME_MAT, AM8, AM9,   PM5,   PM8) ! 0110
!                  CALL DAY_CYCLE(RDAY_SHOP_CHILDREN,    TIME_MAT, PM4, PM5,   PM530, PM7) ! 0110
                  CALL DAY_CYCLE(RDAY_SHOP_CHILDREN,    TIME_MAT, PM8, PM9,   PM5, PM8) ! 0110
                  VIRUS_N_H_AIM(1) = RDAY_SHOP_WORKING     * ( SHOP_NO_WORKING *(1.-RWEEK_END)    &
                                                               + SHOP_NO_WORKING_WEEKEND    * RWEEK_END ) 
                  VIRUS_N_H_AIM(2) = RDAY_SHOP_NOT_WORKING * ( SHOP_NO_NOT_WORKING*(1.-RWEEK_END) &
                                                               + SHOP_NO_NOT_WORKING_WEEKEND* RWEEK_END ) 
                  VIRUS_N_H_AIM(3) = RDAY_SHOP_CHILDREN    * SHOP_NO_CHILDREN * RWEEK_END
                  VIRUS_N_H_AIM(4) = 0.0
!                  LAMBDA_H_H=10000.0*1.0/one_day ! Push some people into shops on time scale of 0.01day. 
               ELSE 
                  VIRUS_N_H_AIM(:)=0.0
                  LAMBDA_H_H= 0.0/one_day ! Push some people out of their homes on time scale of 0.1day. 
               ENDIF

!               IF(EIGEN_VALUE) LAMBDA_H_H=1000.0*LAMBDA_H_H ! experimental THE BEST
!               IF(EIGEN_VALUE) LAMBDA_H_H=100000.0*LAMBDA_H_H ! experimental

            if(.false.) then
                  print *,'time_mat,acctim,rday_sin,rweek_end:',time_mat,acctim,rday_sin,rweek_end
                  print *,'no of days:',time_mat/(24.*3600.)

                  print *,'H_NO_WORKING=',H_NO_WORKING
                  print *,'H_NO_NOT_WORKING=',H_NO_NOT_WORKING
                  print *,'H_NO_CHILDREN=',H_NO_CHILDREN
                  print *,'H_NO_SHIELDING=',H_NO_SHIELDING
! homes...
                  print *,'homes:'
                  CALL DAY_CYCLE(RDAY_OFFICE, TIME_MAT, AM7, AM9,   PM5, PM8) ! 0110
                  RDAY_OFFICE=1.0-RDAY_OFFICE
                  CALL DAY_CYCLE(RDAY_CHILDREN, TIME_MAT, AM8, AM9,   PM4, PM7) ! 0110
                  RDAY_CHILDREN=1.0-RDAY_CHILDREN
                  PRINT *,'VIRUS_N_H_AIM(1) =',  ( RDAY_OFFICE*(1.-RWEEK_END) + MAX(0.5,RDAY_OFFICE)* RWEEK_END )* H_NO_WORKING
                  PRINT *,'VIRUS_N_H_AIM(2) =',  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING
                  PRINT *,'VIRUS_N_H_AIM(3) =', ( RDAY_CHILDREN*(1.-RWEEK_END) + MAX(0.5,RDAY_CHILDREN)* RWEEK_END )* H_NO_CHILDREN
!                  VIRUS_N_H_AIM(1) =  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_WORKING
!                  VIRUS_N_H_AIM(2) =  0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_NOT_WORKING
!                  VIRUS_N_H_AIM(3) = 0.666*((1.0-RDAY_SIN) + 0.5) * H_NO_CHILDREN
                  PRINT *,'VIRUS_N_H_AIM(4) =', ( 0.9 + 0.1* 0.666*((1.0-RDAY_SIN) + 0.5) )  * H_NO_SHIELDING
! 
! office...
                  print *,'office:'
                  CALL DAY_CYCLE(RDAY_OFFICE, TIME_MAT, AM7, AM9,   PM5, PM6) ! 0110
                  print *,'VIRUS_N_H_AIM(1) =', RDAY_OFFICE * OFFICE_NO_WORKING * (1.0 - RWEEK_END) 
! 
! school...
                  print *,'school:'
                  CALL DAY_CYCLE(RDAY_SCHOOL, TIME_MAT, AM8, AM9,   PM330, PM4) ! 0110
                  print *,'VIRUS_N_H_AIM(1) =', RDAY_SCHOOL * SCHOOL_NO_WORKING *  (1.0 - RWEEK_END)
                  print *,'VIRUS_N_H_AIM(3) =', RDAY_SCHOOL * SCHOOL_NO_CHILDREN *  (1.0 - RWEEK_END)
! schops...
                  print *,'shops:'
                  CALL DAY_CYCLE(RDAY_SHOP_WORKING,     TIME_MAT, AM8, AM9,   PM5,   PM8) ! 0110
                  CALL DAY_CYCLE(RDAY_SHOP_NOT_WORKING, TIME_MAT, AM8, AM9,   PM5,   PM8) ! 0110
                  CALL DAY_CYCLE(RDAY_SHOP_CHILDREN,    TIME_MAT, PM4, PM5,   PM530, PM7) ! 0110
                  print *,'VIRUS_N_H_AIM(1) =', RDAY_SHOP_WORKING     * ( SHOP_NO_WORKING     &
                                                                          + SHOP_NO_WORKING_WEEKEND    * RWEEK_END ) 
                  print *,'VIRUS_N_H_AIM(2) =', RDAY_SHOP_NOT_WORKING * ( SHOP_NO_NOT_WORKING &
                                                                          + SHOP_NO_NOT_WORKING_WEEKEND* RWEEK_END ) 
                  print *,'VIRUS_N_H_AIM(3) =', RDAY_SHOP_CHILDREN    * SHOP_NO_CHILDREN * RWEEK_END
               stop 922
            endif ! if(.true.) then


!               PRINT *,'NPLACE, NPEOPLE,NSEIIR:',NPLACE, NPEOPLE,NSEIIR
!               PRINT *,'NPLACE*NPEOPLE*NSEIIR:',NPLACE*NPEOPLE*NSEIIR
               DO IPLACE_PEOPLE=1,NPLACE*NPEOPLE
                  VIRUS_N(IPLACE_PEOPLE) = MAX(TOLER, SUM( T_NEW(I,J,K,1+(IPLACE_PEOPLE-1)*NSEIIR:IPLACE_PEOPLE*NSEIIR) ) )
               END DO
               VIRUS_PLACE_N=0.0
               DO IPLACE=1,NPLACE
                  DO IPEOPLE=1,NPEOPLE
                     IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                     VIRUS_PLACE_N(IPLACE) = VIRUS_PLACE_N(IPLACE) + VIRUS_N(IPLACE_PEOPLE) 
                  END DO
               END DO

! Diffusion...
!               KK(I,J,K,:)=0.0 ! dont have diffusion in houses, shops, schools,hospitals, offcies. 
            IF(SWITCH_ON_DIFFUSION) THEN
               DO ISIR=1,NSEIIR
                  DO ILOOP=1,3
                     IF(ILOOP==1) IPLACE=1
                     IF(ILOOP==2) IPLACE=12
                     IF(ILOOP==3) IPLACE=13
                     DO IPEOPLE=1,NPEOPLE
                        IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                        IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
!                        IF(IPLACE==12) THEN
                        IF(IS_BUILDING .AND. (IPLACE==1) ) THEN
                           KK(I,J,K,IG)=0.1*2.5*length**2/one_day ! normal background transport that is slow
                        ENDIF
                        IF(IPLACE==12) THEN
! the higher the value of K the less cars are on the road.
!                           KK(I,J,K,IG)=100.*2.5*length**2/one_day ! normal background transport that is slow
!                           KK(I,J,K,IG)=2.5*length**2/one_day ! normal background transport that is slow
                           KK(I,J,K,IG)=0.1*2.5*length**2/one_day ! normal background transport that is slow
!                           KK(I,J,K,IG)=100.*2.5*length**2/one_day ! normal background transport that is slow
! good also                           KK(I,J,K,IG)=300.*2.5*length**2/one_day ! normal background transport that is slow
!                           KK(I,J,K,IG)=1000.*2.5*length**2/one_day ! normal background transport that is slow
!                           KK(I,J,K,IG)=10000.*2.5*length**2/one_day ! normal background transport that is slow
                        ENDIF
                        IF(IPLACE==13) THEN ! Portal...
!                           KK(I,J,K,IG)=100000000.*2.5*length**2/one_day ! make large for portal for people to go to hospital.
                           KK(I,J,K,IG)=10000.*2.5*length**2/one_day ! make large for portal for people to go to hospital.
                        ENDIF
                     END DO
                  END DO
               END DO
            ENDIF ! IF(SWITCH_ON_DIFFUSION) THEN


! BETA1, BETA2 coefficients...
            IF(SWITCH_ON_BETA) THEN
               ISIR=1
               DO IPLACE=1,NPLACE
                  DO IPEOPLE=1,NPEOPLE
                     DO JPEOPLE=1,NPEOPLE
                        IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                        JPLACE_PEOPLE = IPLACE + (JPEOPLE-1)*NPLACE
                        IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                        JG=ISIR + (JPLACE_PEOPLE-1)*NSEIIR 
!                        VALUE1= BETA1(IPLACE,IPEOPLE,JPEOPLE) *T_NEW(I,J,K,IG)/MAX(VIRUS_N(IPLACE_PEOPLE),TOLER) 
!                        VALUE2= BETA2(IPLACE,IPEOPLE,JPEOPLE) *T_NEW(I,J,K,IG)/MAX(VIRUS_N(IPLACE_PEOPLE),TOLER) 

                     IF(EIGEN_VALUE) THEN ! Put into fission term
                        VALUE1= BETA1(IPLACE,IPEOPLE,JPEOPLE) *VIRUS_N(IPLACE_PEOPLE)/MAX(VIRUS_PLACE_N(IPLACE),TOLER) 
                        VALUE2= BETA2(IPLACE,IPEOPLE,JPEOPLE) *VIRUS_N(IPLACE_PEOPLE)/MAX(VIRUS_PLACE_N(IPLACE),TOLER) 
                     ELSE
                        VALUE1= BETA1(IPLACE,IPEOPLE,JPEOPLE) *T_NEW(I,J,K,IG)/MAX(VIRUS_PLACE_N(IPLACE),TOLER) 
                        VALUE2= BETA2(IPLACE,IPEOPLE,JPEOPLE) *T_NEW(I,J,K,IG)/MAX(VIRUS_PLACE_N(IPLACE),TOLER) 
                     ENDIF

                     IF(EIGEN_VALUE.AND.(.NOT.EIGEN_METHOD1)) THEN ! Put into fission term -METHOD 2:
                        CALL ADD_IN_SIGMA(I,J,K, IG,JG+2, -VALUE1,   &
                           SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_F, NX,NY,NZ,NG,NG2, RECORD_F, SPAR_RECORD_F)
                        CALL ADD_IN_SIGMA(I,J,K, IG,JG+3, -VALUE2,   &
                           SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_F, NX,NY,NZ,NG,NG2, RECORD_F, SPAR_RECORD_F)
! 
                        CALL ADD_IN_SIGMA(I,J,K, IG+1,JG+2, +VALUE1,   &
                           SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_F, NX,NY,NZ,NG,NG2, RECORD_F, SPAR_RECORD_F)
                        CALL ADD_IN_SIGMA(I,J,K, IG+1,JG+3, +VALUE2,   &
                           SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_F, NX,NY,NZ,NG,NG2, RECORD_F, SPAR_RECORD_F)
                     ELSE 
!                        print *,'beta1:',beta1
!                        print *,'beta2:',beta2
!                        stop 382
                        CALL ADD_IN_SIGMA(I,J,K, IG,JG+2, VALUE1,   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                                            RECORD_S, SPAR_RECORD_S)
                        CALL ADD_IN_SIGMA(I,J,K, IG,JG+3, VALUE2,   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                                            RECORD_S, SPAR_RECORD_S)
! 
                        CALL ADD_IN_SIGMA(I,J,K, IG+1,JG+2, -VALUE1,   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                                            RECORD_S, SPAR_RECORD_S)
                        CALL ADD_IN_SIGMA(I,J,K, IG+1,JG+3, -VALUE2,   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                                            RECORD_S, SPAR_RECORD_S)
                     ENDIF
                     END DO 
                  END DO
               END DO
            ENDIF ! IF(SWITCH_ON_BETA) THEN


! Births result in more children in that group...
               ICHILD=3
               ISIR_CHILD=1
               DO IPLACE=1,NPLACE
               DO IPEOPLE=ICHILD,ICHILD 
                  IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                  IG_CHILD = ISIR_CHILD + (IPLACE_PEOPLE-1)*NSEIIR
! Mu OFF DIAGONAL terms...
                  DO JSIR=1,NSEIIR
                  DO JPEOPLE=1,NPEOPLE
                  IF(JPEOPLE.NE.ICHILD) THEN
                     JPLACE_PEOPLE = IPLACE + (JPEOPLE-1)*NPLACE
                     JG = JSIR + (JPLACE_PEOPLE-1)*NSEIIR
                     CALL ADD_IN_SIGMA(I,J,K,IG_CHILD,JG,-MU,   &
                          SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                          RECORD_S, SPAR_RECORD_S)
                  ENDIF
                  END DO
                  END DO
               END DO
               END DO

            IF(SWITCH_ON_AGEING) THEN ! switch on aging...
! Move children to adult groups with a rate alpha=1/16years and maintaining ratio of groups
               ICHILD=3
               DO ISIR=1,NSEIIR
               DO IPLACE=1,NPLACE
               DO IPEOPLE=ICHILD,ICHILD 
                  IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                  IG_CHILD = ISIR + (IPLACE_PEOPLE-1)*NSEIIR
! decay rates that go on the diagonal..
                  SIGMAT(I,J,K,IG_CHILD) = SIGMAT(I,J,K,IG_CHILD)   + ALPHA 
! Mu OFF DIAGONAL terms...
                  DO JPEOPLE=1,NPEOPLE
                  IF(JPEOPLE.NE.ICHILD) THEN
                     JPLACE_PEOPLE = IPLACE + (JPEOPLE-1)*NPLACE
                     JG = ISIR + (JPLACE_PEOPLE-1)*NSEIIR
                     CALL ADD_IN_SIGMA(I,J,K, JG,IG_CHILD, -ALPHA*COME_AGE_SPECTRUM(JPEOPLE),   &
                          SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                          RECORD_S, SPAR_RECORD_S)
                  ENDIF
                  END DO
               end do
               end do
               end do
             ENDIF ! IF(SWITCH_ON_AGEING) THEN 

             IF(SWITCH_ON_DECAY_RATES) THEN
! decay rates...
               DO IPLACE=1,NPLACE
               DO IPEOPLE=1,NPEOPLE
                  IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                  idisp = (IPLACE_PEOPLE-1)*NSEIIR

! decay rates that go on the diagonal..
                  SIGMAT(I,J,K,1+idisp) = SIGMAT(I,J,K,1+idisp)                                                  + NU 
                  SIGMAT(I,J,K,2+idisp) = SIGMAT(I,J,K,2+idisp)   +VIRUS_SIGMA1(IPEOPLE)                         + NU 
                  SIGMAT(I,J,K,3+idisp) = SIGMAT(I,J,K,3+idisp)   +VIRUS_GAMMA1(IPEOPLE) + VIRUS_SIGMA2(IPEOPLE) + NU 
                  SIGMAT(I,J,K,4+idisp) = SIGMAT(I,J,K,4+idisp)   +VIRUS_GAMMA2(IPEOPLE)                         + NU_I2 
                  SIGMAT(I,J,K,5+idisp) = SIGMAT(I,J,K,5+idisp)   +XI                                            + NU 
               end do
               end do
! The various decay rates...
                  ISIR=1
!             if(is_shops.and.(.not.RECORD_S)) then
!            print *,'VIRUS_SIGMA2(:):',VIRUS_SIGMA2(:) 
!            stop 281
!             endif
!          VIRUS_SIGMA2(:)=0.0
                  DO IPLACE=1,NPLACE
                     DO IPEOPLE=1,NPEOPLE
                        IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                        IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                        CALL ADD_IN_SIGMA(I,J,K, IG,IG+4, - XI,   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                           RECORD_S, SPAR_RECORD_S)

! exposed people resulting in infections...
                        CALL ADD_IN_SIGMA(I,J,K, IG+2,IG+1, -VIRUS_SIGMA1(IPEOPLE) *(1.-R_EIGEN1),   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                           RECORD_S, SPAR_RECORD_S)

! infected people that are very ill...
                        CALL ADD_IN_SIGMA(I,J,K, IG+3,IG+2, -VIRUS_SIGMA2(IPEOPLE),   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                           RECORD_S, SPAR_RECORD_S)

! infected people that are recovered...
                        CALL ADD_IN_SIGMA(I,J,K, IG+4,IG+2, -VIRUS_GAMMA1(IPEOPLE),   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                           RECORD_S, SPAR_RECORD_S)

! infected people that are recovered...
                        CALL ADD_IN_SIGMA(I,J,K, IG+4,IG+3, -VIRUS_GAMMA2(IPEOPLE),   &
                           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, &
                           RECORD_S, SPAR_RECORD_S)
                     END DO
                  END DO
            ENDIF ! IF(SWITCH_ON_DECAY_RATES) THEN
! 
! 
! ************from Home(office, shops, school) to diffusion field and visa versa(START)*****************************
               IPLACE_DIFF1=12
               IPLACE_H=1
               DO ISIR=1,NSEIIR
                  DO IPEOPLE=1,NPEOPLE
                     IPLACE_PEOPLE_H = IPLACE_H + (IPEOPLE-1)*NPLACE
                     IPLACE_PEOPLE_DIFF1 = IPLACE_DIFF1 + (IPEOPLE-1)*NPLACE
! Normalization is needed for eigen-value problems i.e. the magnitudes dont matter but the relative values do...
                  IF(EIGEN_VALUE) THEN
!                     FORCE = ( VIRUS_N(IPLACE_PEOPLE_H)/N_NORMALIZE(IPEOPLE) - VIRUS_N_H_AIM(IPEOPLE)/FIXED_N_NORMALIZE(IPEOPLE) ) &
!                           /MAX( TOLER, VIRUS_N(IPLACE_PEOPLE_H)/N_NORMALIZE(IPEOPLE),  VIRUS_N_H_AIM(IPEOPLE)/FIXED_N_NORMALIZE(IPEOPLE) ) 
                     FORCE = ( VIRUS_N(IPLACE_PEOPLE_H)/SUM(N_NORMALIZE(:)) - VIRUS_N_H_AIM(IPEOPLE)/SUM(FIXED_N_NORMALIZE(:)) ) &
                           /MAX( TOLER, VIRUS_N(IPLACE_PEOPLE_H)/SUM(N_NORMALIZE(:)),  &
                                                     VIRUS_N_H_AIM(IPEOPLE)/SUM(FIXED_N_NORMALIZE(:)) )  
                  ELSE
                     FORCE = ( VIRUS_N(IPLACE_PEOPLE_H) - VIRUS_N_H_AIM(IPEOPLE) ) &
                           /MAX( TOLER, VIRUS_N(IPLACE_PEOPLE_H),  VIRUS_N_H_AIM(IPEOPLE) )  
                  ENDIF
! 
                     H2M=0.5+0.5*SIGN(1.0, FORCE ) 
                     IG_H    =ISIR + (IPLACE_PEOPLE_H-1)*NSEIIR
                     IG_DIFF1=ISIR + (IPLACE_PEOPLE_DIFF1-1)*NSEIIR 

!                     VALUE1 = LAMBDA_H_H * H2M *0.01  *0.001 *1.0
!                     VALUE1 = H2M *0.001  
!                     VALUE1 = LAMBDA_H_H * H2M *0.01!*0.01
!                     VALUE1 = LAMBDA_H_H * H2M * ABS(FORCE) *0.01 ! the one
!                     VALUE1 = LAMBDA_H_H * H2M * ABS(FORCE)  ! the one
                     VALUE1 = min(LAMBDA_H_H * H2M * ABS(FORCE), 1./(1.0*dt))  ! the one
!                     VALUE1 = LAMBDA_H_H * H2M * ABS(FORCE) 
!                     VALUE1 = LAMBDA_H_H * H2M * ABS(FORCE) 
!                     VALUE2 = LAMBDA_H_H * (1.0-H2M) * ABS(FORCE) ! the one
                     VALUE2 = min(LAMBDA_H_H * (1.0-H2M) * ABS(FORCE), 1./(1.0*dt)) ! the one
!                     S(I,J,K,IG_H) = S(I,J,K,IG_H)         - VALUE1 * T_NEW(I,J,K,IG_H)
!                     S(I,J,K,IG_DIFF1) = S(I,J,K,IG_DIFF1) + VALUE1 * T_NEW(I,J,K,IG_H)
!                     VALUE2 = LAMBDA_H_H * (1.0-H2M) 
             if(.false.) then
               IF(IS_HOMES_OCCUPIED) THEN ! bottom middle region
                  IF(IPEOPLE==1) THEN
                     MIN_VALUE1_W=MIN(VALUE1, MIN_VALUE1_W)
                     MIN_VALUE2_W=MIN(VALUE2, MIN_VALUE2_W)
                     MAX_VALUE1_W=MAX(VALUE1, MAX_VALUE1_W)
                     MAX_VALUE2_W=MAX(VALUE2, MAX_VALUE2_W)
                 ENDIF
                  IF(IPEOPLE==2) THEN
                       if(MAX_VALUE2_NW + MAX_VALUE1_NW<1.e-5) then
         print *,'IPLACE_PEOPLE_H, VIRUS_N(IPLACE_H + (IPEOPLE-2)*NPLACE):',IPLACE_PEOPLE_H, VIRUS_N(IPLACE_H + (IPEOPLE-2)*NPLACE)
         print *,'force, VIRUS_N(IPLACE_PEOPLE_H), VIRUS_N_H_AIM(IPEOPLE):',force, VIRUS_N(IPLACE_PEOPLE_H), VIRUS_N_H_AIM(IPEOPLE)
                       endif 
                     MIN_VALUE1_NW=MIN(VALUE1, MIN_VALUE1_NW)
                     MIN_VALUE2_NW=MIN(VALUE2, MIN_VALUE2_NW)
                     MAX_VALUE1_NW=MAX(VALUE1, MAX_VALUE1_NW)
                     MAX_VALUE2_NW=MAX(VALUE2, MAX_VALUE2_NW)
                  ENDIF
               ENDIF 
             endif
! Move people from Home to Diffusion1 fields...
                   if(.true.) then
                     SIGMAT(I,J,K,IG_H) =  SIGMAT(I,J,K,IG_H)  +  VALUE1 ! take away from home.
!                   if(.false.) then
                     CALL ADD_IN_SIGMA(I,J,K,IG_DIFF1, IG_H, - VALUE1,   & ! Add to diffusion 1 field
             SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                   endif

                 if(.FALSE.) then
                     SIGMAT(I,J,K,IG_DIFF1) =  SIGMAT(I,J,K,IG_DIFF1)  +  VALUE1 ! take away from home.
                     CALL ADD_IN_SIGMA(I,J,K, IG_H, IG_DIFF1, - VALUE1,   & ! Add to diffusion 1 field
             SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                 endif
! 
! Move people from Diffusion1 to Home fields...
                     SIGMAT(I,J,K,IG_DIFF1) =  SIGMAT(I,J,K,IG_DIFF1)  +  VALUE2 ! take away from diffusion field.
                     CALL ADD_IN_SIGMA(I,J,K, IG_H, IG_DIFF1, - VALUE2,   & ! Add to diffusion 1 field
             SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                 if(.FALSE.) then
                     SIGMAT(I,J,K,IG_H) =  SIGMAT(I,J,K,IG_H)  +  VALUE2 ! take away from diffusion field.
                     CALL ADD_IN_SIGMA(I,J,K, IG_DIFF1, IG_H,  - VALUE2,   & ! Add to diffusion 1 field
             SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                 endif

                  end do
               end do
! ************from Home(office, shops, school) to diffusion field and visa versa(END)*******************************
! 
! 
! ************add the constraint for the eigen-value problem(START)********************
            IF(SWITCH_ON_CONSTRAIN_PEOPLE_RATIO) THEN 
               SUM_FIXED_N_NORMALIZE = SUM(FIXED_N_NORMALIZE(:))
!               LAMBDA_C_C=LAMBDA_H_H ! beetter ratio
               LAMBDA_C_C=LAMBDA_H_H*1000.0 ! beetter ratio
               IPLACE_H=1
               ISIR=1
               DO JSIR=1,NSEIIR
                  DO IPEOPLE=1,NPEOPLE
                  DO JPEOPLE=1,NPEOPLE
!                  DO JPLACE=1,NPLACE
                     JPLACE=IPLACE_H
                     IPLACE_PEOPLE_H = IPLACE_H + (IPEOPLE-1)*NPLACE
                     JPLACE_PEOPLE = JPLACE + (JPEOPLE-1)*NPLACE
! Normalization is needed for eigen-value problems i.e. the magnitudes dont matter but the relative values do...

                     IG_H    =ISIR + (IPLACE_PEOPLE_H-1)*NSEIIR
                     IIG = JSIR + (IPLACE_PEOPLE_H-1)*NSEIIR
                     JG=JSIR + (JPLACE_PEOPLE-1)*NSEIIR 

!                     VALUE1 = LAMBDA_H_H * H2M *0.01
                     VALUE1 = LAMBDA_C_C * FIXED_N_NORMALIZE(JPEOPLE)/SUM_FIXED_N_NORMALIZE
!                     VALUE1 = LAMBDA_H_H * H2M * ABS(FORCE) 
!                     VALUE1 = LAMBDA_H_H * H2M * ABS(FORCE) 
                     VALUE2 = LAMBDA_C_C * FIXED_N_NORMALIZE(IPEOPLE)/SUM_FIXED_N_NORMALIZE
! Move people from Home to Diffusion1 fields...
                  IF(IG_H==IIG) THEN
                     SIGMAT(I,J,K,IG_H) =  SIGMAT(I,J,K,IG_H)  +  VALUE1 ! take away from home.
                  ELSE
!                     CALL ADD_IN_SIGMA(I,J,K,JG, IG_H, - VALUE2,   & ! Add to diffusion 1 field
                     CALL ADD_IN_SIGMA(I,J,K, IG_H, IIG, + VALUE1,   & ! Add to diffusion 1 field
             SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                  ENDIF
                     CALL ADD_IN_SIGMA(I,J,K, IG_H, JG, - VALUE2,   & ! Add to diffusion 1 field
             SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
! 

!                  end do
                  end do
                  end do
               end do
            ENDIF ! IF(SWITCH_ON_CONSTRAIN_PEOPLE_RATIO) THEN 
! ************add the constraint for the eigen-value problem(END)**********************
! 
! 
! ******************** ROADS AND PATHWAYS START************************************
        IF(SWITCH_ON_ROAD_CAR_WALK) THEN
!        IF(SWITCH_ON_ROAD_CAR_WALK.and.(i<=nx-10).and.(j>=10)) THEN
! roads and pathways...
          DO ISIR=1,NSEIIR
          DO IPLACE=1,NPLACE
          DO IPEOPLE=1,NPEOPLE
             IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
             IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
! cars...
             if(rd_mark_x(i)==1) then
             IF(IS_BIG_ROAD_KEEP(I,J,K).AND.IS_BIG_ROAD_KEEP(I,J-1,K) ) THEN ! UP and Down
                IF(IPLACE==2) U(2, I,J,K,IG)= 1.0 * CAR_SPEED
                IF(IPLACE==3) U(2, I,J,K,IG)=-1.0 * CAR_SPEED
             ENDIF
             endif
             if(rd_mark_y(j)==1) then
             IF(IS_BIG_ROAD_KEEP(I,J,K).AND.IS_BIG_ROAD_KEEP(I-1,J,K) ) THEN ! right and left.
                IF(IPLACE==4) U(1, I,J,K,IG)= 1.0 * CAR_SPEED
                IF(IPLACE==5) U(1, I,J,K,IG)=-1.0 * CAR_SPEED
             ENDIF
             endif
! walking...
             if(rd_mark_x(i)==1) then
             IF(IS_BIG_ROAD_KEEP(I,J,K).AND.IS_BIG_ROAD_KEEP(I,J-1,K) ) THEN ! UP and Down
                IF(IPLACE==7) U(2, I,J,K,IG)= 1.0 * WALK_SPEED
                IF(IPLACE==8) U(2, I,J,K,IG)=-1.0 * WALK_SPEED
             ENDIF
             endif
             if(rd_mark_y(j)==1) then
             IF(IS_BIG_ROAD_KEEP(I,J,K).AND.IS_BIG_ROAD_KEEP(I-1,J,K) ) THEN ! right and left.
                IF(IPLACE==9) U(1, I,J,K,IG)= 1.0 * WALK_SPEED
                IF(IPLACE==10) U(1, I,J,K,IG)=-1.0 * WALK_SPEED
             ENDIF
             endif
          END DO
          END DO
          END DO
! Make sure there is no leakage out of the domain...
          U(2, I,1:2,K,IG)=0.0;  U(2, I,NY:NY+1,K,IG)=0.0
          U(1, 1:2,J,K,IG)=0.0;  U(1, NX:NX+1,J,K,IG)=0.0

       IF(SWITCH_ON_CARS_WALK_THROUGH_BUILDINGS) THEN
          cross=.false.
          if( ((rd_mark_x(i-1)==1).and.(rd_mark_x(i)==1)) &
          .or.((rd_mark_x(i)==1).and.(rd_mark_x(i+1)==1))) then !.and.IS_BIG_ROAD) then ! Travel up and down.
          if( ((rd_mark_y(j-1)==1).and.(rd_mark_y(j)==1)) &
          .or.((rd_mark_y(j)==1).and.(rd_mark_y(j+1)==1))) then !.and.IS_BIG_ROAD) then ! Travel up and down.
              cross=.true.
          endif
          endif
       else
          cross=.false.
          if( ((rd_mark_x(i-1)==1).and.(rd_mark_x(i)==1)) &
          .or.((rd_mark_x(i)==1).and.(rd_mark_x(i+1)==1)) .and.IS_BIG_ROAD) then ! Travel up and down.
          if( ((rd_mark_y(j-1)==1).and.(rd_mark_y(j)==1)) &
          .or.((rd_mark_y(j)==1).and.(rd_mark_y(j+1)==1)) .and.IS_BIG_ROAD) then ! Travel up and down.
              cross=.true.
          endif
          endif
       endif
!          road_ends=(i>nx-4).or.(i<4).or.(j>ny-4).or.(j<4) 

          if(cross) then ! swap people over between roads...
!          if(cross.and.(.not.is_building)) then ! swap people over between roads...
!          if(.false.) then ! swap people over between roads...
             IPLACE_DIFF1=12
             DO ISIR=1, NSEIIR 
! There are 3 possible directions to travel (carry on, left or right depending on size of diffusion field cells) 
                DO IVEH=0,1 ! =0 cars, =1 Pedestrains....
                DO IPLEVEL=1,4
                   IPLACE= 1 + IPLEVEL + 5*IVEH
                   DO JPLEVEL=1,4
                      JPLACE= 1 + JPLEVEL + 5*IVEH
                      TOP   =((IPLEVEL==1).OR.(IPLEVEL==2)).AND.(JPLEVEL>2)
                      BOTTOM=((IPLEVEL==3).OR.(IPLEVEL==4)).AND.(JPLEVEL<3)
                      IF(TOP.OR.BOTTOM) THEN
                         DO IPEOPLE=1,NPEOPLE
                            IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                            IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                            IPLACE_PEOPLE_DIFF1 = IPLACE_DIFF1 + (IPEOPLE-1)*NPLACE
                            JPLACE_PEOPLE = JPLACE + (IPEOPLE-1)*NPLACE
                            IG_DIFF1=ISIR + (IPLACE_PEOPLE_DIFF1-1)*NSEIIR 
                            JG=ISIR + (JPLACE_PEOPLE-1)*NSEIIR 
                         if(.not.UP_WIND_CARS_WALK) then ! original...
!                         if(.false.) then ! original...
                            DIFF_VERT1= (T_NEW(I,J+1,K,IG_DIFF1) - T_NEW(I,J,K,IG_DIFF1))/DY
                            DIFF_VERT2= (T_NEW(I,J,K,IG_DIFF1) - T_NEW(I,J-1,K,IG_DIFF1))/DY

                            DIFF_HORI1= (T_NEW(I+1,J,K,IG_DIFF1) - T_NEW(I,J,K,IG_DIFF1))/DX
                            DIFF_HORI2= (T_NEW(I,J,K,IG_DIFF1) - T_NEW(I-1,J,K,IG_DIFF1))/DX
                            IF(J==2)    DIFF_VERT2=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF(J==NY-1) DIFF_VERT1= INFINY
                            IF(I==2)    DIFF_HORI2=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF(I==NX-1) DIFF_HORI1= INFINY
                         else ! seems to be the best...
                            DIFF_VERT2= (T_NEW(I,J+1,K,IG_DIFF1) - T_NEW(I,J,K,IG_DIFF1))/DY
                            DIFF_VERT1= (T_NEW(I,J,K,IG_DIFF1) - T_NEW(I,J-1,K,IG_DIFF1))/DY

                            DIFF_HORI2= (T_NEW(I+1,J,K,IG_DIFF1) - T_NEW(I,J,K,IG_DIFF1))/DX
                            DIFF_HORI1= (T_NEW(I,J,K,IG_DIFF1) - T_NEW(I-1,J,K,IG_DIFF1))/DX
                    if(.false.) then ! original...
                            IF(J==2)    DIFF_VERT1=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF(J==NY-1) DIFF_VERT2= INFINY
                            IF(I==2)    DIFF_HORI1=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF(I==NX-1) DIFF_HORI2= INFINY
                    endif
                    if(.true.) then !
                            IF((J==2).or.(j==3).or.(j==4))    DIFF_VERT1=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF((J==NY-1).or.(j==NY-2).or.(j==NY-3)) DIFF_VERT2= INFINY
                            IF((I==2).or.(I==3).or.(I==4))    DIFF_HORI1=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF((I==NX-1).or.(I==NX-2).or.(I==NX-3)) DIFF_HORI2= INFINY
                    endif
                    if(.false.) then ! does not work so well...
                            IF(J==2)    DIFF_VERT2=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF(J==NY-1) DIFF_VERT1= INFINY
                            IF(I==2)    DIFF_HORI2=-INFINY ! Modify near boundary so cars and people dont accumulate in cells here.
                            IF(I==NX-1) DIFF_HORI1= INFINY
                    endif
                         endif

                            SWAP=.FALSE. ! SWAP OVER TRAFFIC.
                            IF(  ((IPLEVEL==1) .AND. ((DIFF_VERT1>0.).or.(J==NY-1).or.(J==NY-2).or.(J==NY-3)) ) & 
                            .or. ((IPLEVEL==2) .AND. ((DIFF_VERT2<0.).or.(J==2).or.(J==3).or.(J==4)   ) )   ) THEN ! Swap over
! Choose the most benefitial of the 2 directions...
                               IF( (DIFF_HORI1<0.0) .and. (JPLEVEL==3) ) SWAP=.TRUE.
                               IF( (DIFF_HORI2>0.0) .and. (JPLEVEL==4) ) SWAP=.TRUE.
                            ENDIF 
                            IF(  ((IPLEVEL==3) .AND. ((DIFF_HORI1>0.).or.(I==NX-1).or.(I==NX-2).or.(I==NX-3)) ) &
                            .or. ((IPLEVEL==4) .AND. ((DIFF_HORI2<0.).or.(I==2).or.(I==3).or.(I==4)   ) )   ) THEN ! Swap over
                               IF( (DIFF_VERT1<0.0) .and. (JPLEVEL==1) ) SWAP=.TRUE.
                               IF( (DIFF_VERT2>0.0) .and. (JPLEVEL==2) ) SWAP=.TRUE.
                            ENDIF 

                            VALUE1 = LAMBDA_RD_WALK_EXCHANGE(IVEH) 
!                            IF(.NOT.SWAP) VALUE1=0.1*VALUE1 
                            IF(.NOT.SWAP) VALUE1=FRAC_GO_WRONG_WAY*VALUE1 ! FRAC_GO_WRONG_WAY=0.0 -no one goes wrong way, =0.1 then 10% go wrong way
!                            IF(.NOT.SWAP) VALUE1=0.25*VALUE1
!                            IF(.NOT.SWAP) VALUE1=1.25*VALUE1
!                            IF(.NOT.SWAP) VALUE1=0.0

                            SIGMAT(I,J,K,IG) = SIGMAT(I,J,K,IG) + VALUE1 ! Move people to another direction...
! Add people into that other direction...
                            CALL ADD_IN_SIGMA(I,J,K, JG,IG, -VALUE1,   &
            SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
! make symmetric...
                          if(.false.) then
                            SIGMAT(I,J,K,JG) = SIGMAT(I,J,K,JG) + VALUE1 ! Move people to another direction...
! Add people into that other direction...
                            CALL ADD_IN_SIGMA(I,J,K, IG,JG, -VALUE1,   &
            SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                          endif
                         END DO ! DO IPEOPLE=1,NPEOPLE
                      ENDIF ! IF(TOP.OR.BOTTOM) THEN
                   END DO ! DO JPLEVEL=1,4
                END DO ! DO IPLEVEL=1,4
                END DO ! DO IVEH=0,1
            END DO ! DO ISIR=1, NSEIIR
          endif ! if(cross) then
! 
! 
! Pick people up and drop them off from diffusion field...
       IF(IS_BIG_ROAD) THEN ! only transfer when there is a road
          VERTICAL_ROAD  =.false.
          HORIZONTAL_ROAD=.false.
          if(  ((rd_mark_x(i-1)==1).and.(rd_mark_x(i)==1))  &
           .or.((rd_mark_x(i+1)==1).and.(rd_mark_x(i)==1)) ) VERTICAL_ROAD=.TRUE.
          if(  ((rd_mark_y(j-1)==1).and.(rd_mark_y(j)==1))  &
           .or.((rd_mark_y(j+1)==1).and.(rd_mark_y(j)==1)) ) HORIZONTAL_ROAD=.TRUE.

          IPLACE_DIFF1=12
          DO ISIR=1, NSEIIR 
! There are 3 possible directions to travel (carry on, left or right depending on size of diffusion field cells) 
             DO IVEH=0,1 ! =0 cars, =1 Pedestrains....
             DO IPLEVEL=1,4
             IF( ((IPLEVEL<=2).AND.VERTICAL_ROAD).OR.((IPLEVEL>=3).AND.HORIZONTAL_ROAD) ) THEN 
                IPLACE= 1 + IPLEVEL + 5*IVEH
                DO IPEOPLE=1,NPEOPLE
                   IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                   IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                   IPLACE_PEOPLE_DIFF1 = IPLACE_DIFF1 + (IPEOPLE-1)*NPLACE
                   IG_DIFF1=ISIR + (IPLACE_PEOPLE_DIFF1-1)*NSEIIR 

!                   VALUE1 = LAMBDA_D_RD &
!                          * (  T_NEW(I,J,K,IG_DIFF1) - T_NEW_LIMIT_RD_WALK(IVEH)  )/ T_NEW_LIMIT_RD_WALK(IVEH)  
!                   SIG_VALUE1 = 0.5 + 0.5*SIGN(1.0,T_NEW(I,J,K,IG_DIFF1) - T_NEW_LIMIT_RD_WALK(IVEH)) 
!                   VALUE12 = (  FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG_DIFF1) - T_NEW(I,J,K,IG)  )/ T_NEW_LIMIT_RD_WALK(IVEH)  
                   VALUE12 = (  T_NEW(I,J,K,IG_DIFF1) - FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG)  )  &
                           / MAX(TOLER, T_NEW(I,J,K,IG_DIFF1), FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG)  )
                   SIG_VALUE1 = 0.5 + 0.5*SIGN(1.0,VALUE12 ) 

                if(.false.) then
                   VALUE12_DROP   = (  2.0*T_NEW(I,J,K,IG_DIFF1) - FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG)  )  &
                                  / MAX(TOLER, T_NEW(I,J,K,IG_DIFF1), FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG)  )
                   VALUE12_PICKUP = (  0.5*T_NEW(I,J,K,IG_DIFF1) - FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG)  )  &
                                  / MAX(TOLER, T_NEW(I,J,K,IG_DIFF1), FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG)  )

                   IF(VALUE12_DROP*VALUE12_PICKUP<0.0) THEN
                       VALUE12_DROP=0.0
                       VALUE12_PICKUP=0.0
                   ENDIF

                   VALUE1 = (1.-SIG_VALUE1) * LAMBDA_D_RD_WALK_DROP(IVEH) * ABS(VALUE12_DROP)
                   VALUE2 = SIG_VALUE1 * LAMBDA_D_RD_WALK(IVEH) * ABS(VALUE12_PICKUP)
               else 

!                   SIG_VALUE1 = 0.5 + 0.5*SIGN(1.0, FRACTION_RD_WALK(IVEH)*T_NEW(I,J,K,IG_DIFF1) - T_NEW(I,J,K,IG) ) 
                   VALUE1 =min( (1.-SIG_VALUE1) * LAMBDA_D_RD_WALK_DROP(IVEH) * ABS(VALUE12), 1./dt)
                   VALUE2 =min( SIG_VALUE1 * LAMBDA_D_RD_WALK(IVEH) * ABS(VALUE12), 1./dt)
               endif 

                   SIGMAT(I,J,K,IG) = SIGMAT(I,J,K,IG) + VALUE1 ! Move from cars or walking to diffusion field (DROP OFF)...
! Add people into diffusion field...
                   CALL ADD_IN_SIGMA(I,J,K, IG_DIFF1,IG, -VALUE1,   &
           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)

                   SIGMAT(I,J,K,IG_DIFF1) = SIGMAT(I,J,K,IG_DIFF1) + VALUE2 ! Move from diffusion field to cars or walking fields...
! Add people into cars or walking fields...
                   CALL ADD_IN_SIGMA(I,J,K, IG, IG_DIFF1, -VALUE2,   &
           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)

                END DO ! DO IPEOPLE=1,NPEOPLE
             ENDIF ! IF( ((IPLEVEL<=2).AND.VERTICAL_ROAD).OR.((IPLEVEL>=3).AND.HORIZONTAL_ROAD) ) THEN
             END DO ! DO IPLEVEL=1,4
             END DO ! DO IVEH=0,1
          END DO ! DO ISIR=1, NSEIIR
       ENDIF ! IF(IS_BIG_ROAD) THEN ! only transfer when there is a road
! 
! 
! ************from Home(office, shops, school) to road (cars or walkers) and visa versa(START)*****************************
            IF(SWITCH_ON_ROAD_TO_BUILDINGS) THEN
            IF(IS_BIG_ROAD) THEN ! only transfer when there is a road
            IF(IS_BUILDING) THEN
!               LAMBDA_R_R = 10.*LAMBDA_H_H
               LAMBDA_R_R = LAMBDA_H_H
               IPLACE_H=1
               DO ISIR=1,NSEIIR
                  DO IPEOPLE=1,NPEOPLE

                     DO IVEH=0,1 ! =0 cars, =1 Pedestrains....
                     DO IPLEVEL=1,4
                     IF( ((IPLEVEL<=2).AND.VERTICAL_ROAD).OR.((IPLEVEL>=3).AND.HORIZONTAL_ROAD) ) THEN 
                        IPLACE_RD= 1 + IPLEVEL + 5*IVEH

                        IPLACE_PEOPLE_H = IPLACE_H + (IPEOPLE-1)*NPLACE
                        IPLACE_PEOPLE_RD = IPLACE_RD + (IPEOPLE-1)*NPLACE
! Normalization is needed for eigen-value problems i.e. the magnitudes dont matter but the relative values do...
                     IF(EIGEN_VALUE) THEN
!                        FORCE = ( VIRUS_N(IPLACE_PEOPLE_H)/N_NORMALIZE(IPEOPLE) - VIRUS_N_H_AIM(IPEOPLE)/FIXED_N_NORMALIZE(IPEOPLE) ) &
!                              /MAX( TOLER, VIRUS_N(IPLACE_PEOPLE_H)/N_NORMALIZE(IPEOPLE),  VIRUS_N_H_AIM(IPEOPLE)/FIXED_N_NORMALIZE(IPEOPLE) ) 
            FORCE = ( VIRUS_N(IPLACE_PEOPLE_H)/SUM(N_NORMALIZE(:)) - VIRUS_N_H_AIM(IPEOPLE)/SUM(FIXED_N_NORMALIZE(:)) ) &
             /MAX( TOLER, VIRUS_N(IPLACE_PEOPLE_H)/SUM(N_NORMALIZE(:)),  VIRUS_N_H_AIM(IPEOPLE)/SUM(FIXED_N_NORMALIZE(:)) )  
                     ELSE
                        FORCE = ( VIRUS_N(IPLACE_PEOPLE_H) - VIRUS_N_H_AIM(IPEOPLE) ) &
                              /MAX( TOLER, VIRUS_N(IPLACE_PEOPLE_H),  VIRUS_N_H_AIM(IPEOPLE) )  
                     ENDIF
! 
                        H2M=0.5+0.5*SIGN(1.0, FORCE ) 
                        IG_H    =ISIR + (IPLACE_PEOPLE_H-1)*NSEIIR
                        IG_RD=ISIR + (IPLACE_PEOPLE_RD-1)*NSEIIR 

                        VALUE1 = min(LAMBDA_R_R * H2M * ABS(FORCE), 30.*1./dt)  ! the one

!                        VALUE2 = 10.*min(LAMBDA_R_R * (1.0-H2M) * ABS(FORCE), 30.*1./dt) ! the one
!                        VALUE2 = min(10.*LAMBDA_R_R * (1.0-H2M) * ABS(FORCE),  3.*1./dt) ! the one(6)
!                        VALUE2 = min(10.*LAMBDA_R_R * (1.0-H2M) * ABS(FORCE),  1.*1./dt) ! the one(7)
!                        VALUE2 = min(10.*LAMBDA_R_R * (1.0-H2M) * ABS(FORCE), 10.*1./dt) ! the one(8)
!                        VALUE2 = min(10.*LAMBDA_R_R * (1.0-H2M) * ABS(FORCE), 30.*1./dt) ! the one (9)
!                        VALUE2 = min(20.*LAMBDA_R_R * (1.0-H2M) * ABS(FORCE), 30.*1./dt) ! the one (9)
                        VALUE2 = min(5.*LAMBDA_R_R * (1.0-H2M) * ABS(FORCE), 30.*1./dt) ! the one (9)
!                        VALUE2 = 100.*min(LAMBDA_R_R * (1.0-H2M) * ABS(FORCE), 30.*1./dt) ! the one

! Move people from Home to Road fields...
                        SIGMAT(I,J,K,IG_H) =  SIGMAT(I,J,K,IG_H)  +  VALUE1 ! take away from home.
                        CALL ADD_IN_SIGMA(I,J,K,IG_RD, IG_H, - VALUE1,   & ! Add to road (cars or walkers). 
            SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
! 
! Move people from Road to Home fields...
                        SIGMAT(I,J,K,IG_RD) =  SIGMAT(I,J,K,IG_RD)  +  VALUE2 ! take away from road (cars or walkers). 
                        CALL ADD_IN_SIGMA(I,J,K, IG_H, IG_RD, - VALUE2,   & ! Add to road (cars or walkers).
            SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                     ENDIF ! IF( ((IPLEVEL<=2).AND.VERTICAL_ROAD).OR.((IPLEVEL>=3).AND.HORIZONTAL_ROAD) ) THEN
                     end do ! DO IPLEVEL=1,4
                     end do ! DO IVEH=0,1

                  end do
               end do
            ENDIF ! IF(IS_BUILDING) THEN
            ENDIF ! IF(IS_BIG_ROAD) THEN
            ENDIF ! IF(SWITCH_ON_ROAD_TO_BUILDINGS) THEN
! ************from Home(office, shops, school) to road (cars or walkers) and visa versa(END)*******************************
! 
! 
       ENDIF ! IF(SWITCH_ON_ROAD_CAR_WALK) THEN
! ******************** ROADS AND PATHWAYS END************************************
! 
! 
! *********************HOSPITAL ATTENDANCE BEGIN*********************************
! Move anyone in virus group I2 into place group 13 which has unlimited diffusion 
! on a time scale of an ambulance (say 15mins)... LAMBDA_ILL_PORTAL 
       IF(SWITCH_ON_HOSPITAL) THEN
          ISIR=4
          IPLACE_PORTAL = 13 ! Portal...
          DO IPLACE=1,NPLACE
             IF(IPLACE.NE.IPLACE_PORTAL) THEN
                DO IPEOPLE=1,NPEOPLE
                   IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                   IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                   VALUE1 = LAMBDA_ILL_PORTAL 
                   SIGMAT(I,J,K,IG) = SIGMAT(I,J,K,IG) + VALUE1 ! Take people away into portal
! Add people into portal...
                   IPLACE_PEOPLE_PORTAL = IPLACE_PORTAL + (IPEOPLE-1)*NPLACE
                   IG_PORTAL=ISIR + (IPLACE_PEOPLE_PORTAL-1)*NSEIIR 
                   CALL ADD_IN_SIGMA(I,J,K, IG_PORTAL,IG, -VALUE1,   &
      SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                END DO ! DO IPEOPLE=1,NPEOPLE
             ENDIF ! IF(IPLACE.NE.IPLACE_GO_TO) THEN
          END DO ! DO IPLACE=1,NPLACE

! Take ALL people from the portal to the hospital (not just the very ill - THE COULD RECOVER WHILE IN PORTAL)...
          IF(IS_HOSPITAL) THEN
! take everyone in the portal to hospital...
             IPLACE_PORTAL = 13 ! Portal...
             IPLACE_H=1
             DO ISIR=1,NSEIIR
                DO IPEOPLE=1,NPEOPLE
                   IPLACE_PEOPLE = IPLACE_PORTAL + (IPEOPLE-1)*NPLACE
                   IG_PORTAL=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                   VALUE1 = LAMBDA_PORTAL_HOSPITAL
                   SIGMAT(I,J,K,IG_PORTAL) = SIGMAT(I,J,K,IG_PORTAL) + VALUE1 ! Take people away in I2 from portal
! Add people into H- GROUP, IN THIS CASE HOSPITAL...
                   IPLACE_PEOPLE_H = IPLACE_H + (IPEOPLE-1)*NPLACE
                   IG_H=ISIR + (IPLACE_PEOPLE_H-1)*NSEIIR 
                   CALL ADD_IN_SIGMA(I,J,K, IG_H, IG_PORTAL, -VALUE1,   &
       SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                END DO ! DO IPEOPLE=1,NPEOPLE
             END DO ! DO ISIR=1,NSEIIR
! 
! Release people from hospital not in group I2 and put them in place group D1. 
             IPLACE_H=1
             IPLACE_DIFF1=12
             DO ISIR=1,NSEIIR
                IF(ISIR.NE.4) THEN
                   DO IPEOPLE=1,NPEOPLE
                      IPLACE_PEOPLE_H = IPLACE_H + (IPEOPLE-1)*NPLACE
                      IG_H=ISIR + (IPLACE_PEOPLE_H-1)*NSEIIR 
                      VALUE1 = LAMBDA_DISCHARG_HOSPITAL
                      SIGMAT(I,J,K,IG_H) = SIGMAT(I,J,K,IG_H) + VALUE1 ! Take people away in I2 from portal
! Add people into H- GROUP, IN THIS CASE HOSPITAL...
                      IPLACE_PEOPLE_DIFF1 = IPLACE_DIFF1 + (IPEOPLE-1)*NPLACE
                      IG_DIFF1=ISIR + (IPLACE_PEOPLE_DIFF1-1)*NSEIIR 
                      CALL ADD_IN_SIGMA(I,J,K, IG_DIFF1, IG_H, -VALUE1,   &
           SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD_S, SPAR_RECORD_S)
                   END DO ! DO IPEOPLE=1,NPEOPLE
                ENDIF
             END DO ! DO ISIR=1,NSEIIR

          ENDIF ! IF(IS_HOSPITAL) THEN
       ENDIF ! IF(SWITCH_ON_HOSPITAL) THEN
! *********************HOSPITAL ATTENDANCE END***********************************


          if(EIGEN_VALUE) then
! eigen-value...
! consistent with the houses...
                  KK(I,J,K,:) = 0.05 * KK(I,J,K,:)

! Set 1st and 5th SIR groups to zero with large number...
                  ISIR=1
                  DO IPLACE=1,NPLACE
                     DO IPEOPLE=1,NPEOPLE
                        IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                        IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                        SIGMAT(I,J,K,IG)  = SIGMAT(I,J,K,IG) + 10000000000000.0/one_day
                        SIGMAT(I,J,K,IG+4)= SIGMAT(I,J,K,IG+4) + 10000000000000.0/one_day
                     END DO
                  END DO

! Set exchange terms between places...                  
! 
               if(EIGEN_METHOD1) then ! method that put the 1/keff over the I - infection eqn

                  ISIR=1
                  DO IPLACE=1,NPLACE
                     DO IPEOPLE=1,NPEOPLE
                        IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                        IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
! FISSION..
                        CALL ADD_IN_SIGMA(I,J,K, IG+2,IG+1, VIRUS_SIGMA1(IPEOPLE) ,   &
            SIGMA_F, SIGMA_F_PTER, FIN_SIGMA_F, COL_SIGMA_F, NCOL_SIGMA_F, NX,NY,NZ,NG,NG2, RECORD_F, SPAR_RECORD_F)
                     END DO
                  END DO

               endif ! if(EIGEN_METHOD1) then 

            endif ! if(EIGEN_VALUE) then
! 
         END DO ! DO I=2,NX-1
         END DO ! DO J=2,NY-1
         END DO ! DO K=2,NZ-1

!         IF(.NOT.SWITCH_ON_ROAD_CAR_WALK) U=0.0 ! No need to move cars people on roads. 


! CONSERVATION OF SOURCE...
!          stop 282

!            s=0.0
!            PRINT *,'AT END OF DERTMINE X_SECTIONS PRINT:'
!           CALL PRINT_NORMAL(S,NX,NY,NZ,NG)

!         PRINT *,'MIN_VALUE1_W, MIN_VALUE2_W, MAX_VALUE1_W, MAX_VALUE2_W:',MIN_VALUE1_W, MIN_VALUE2_W, MAX_VALUE1_W, MAX_VALUE2_W
!         PRINT *,'MIN_VALUE1_NW, MIN_VALUE2_NW, MAX_VALUE1_NW, MAX_VALUE2_NW:',MIN_VALUE1_NW, MIN_VALUE2_NW, MAX_VALUE1_NW, MAX_VALUE2_NW


         END SUBROUTINE DETERMINE_X_SECTIONS_NG260



         SUBROUTINE PRINT_NORMAL(S,NX,NY,NZ,NG)
         INTEGER NX,NY,NZ,NG
         REAL S(NX,NY,NZ,NG)
! Local ...
         REAL FIXED_N_NORMALIZE(4), N_NORMALIZE(4)
         INTEGER NSEIIR,NPLACE,NPEOPLE, IG, ISIR, IPLACE, IPEOPLE, IPLACE_PEOPLE

         N_NORMALIZE(:)=0.0
         NSEIIR=5
         NPLACE=13
         NPEOPLE=4
         DO K=2,NZ-1
         DO J=2,NY-1
         DO I=2,NX-1
! normalization for this simulation...
                  DO ISIR=1,NSEIIR
                  DO IPLACE=1,NPLACE
                  DO IPEOPLE=1,NPEOPLE
                     IPLACE_PEOPLE = IPLACE + (IPEOPLE-1)*NPLACE
                     IG=ISIR + (IPLACE_PEOPLE-1)*NSEIIR 
                     N_NORMALIZE(IPEOPLE) = N_NORMALIZE(IPEOPLE) + S(I,J,K,IG)
                  END DO
                  END DO
                  END DO

         END DO ! DO I=2,NX-1
         END DO ! DO J=2,NY-1
         END DO ! DO K=2,NZ-1
!          stop 282

!            s=0.0
            PRINT *,'MAXVAL(S),MINVAL(S):',MAXVAL(S),MINVAL(S)
            PRINT *,'CONSERVATION OF SOURCE:'
            print *,'N_NORMALIZE:',N_NORMALIZE
            print *,'N_NORMALIZE/SUM(N_NORMALIZE):',N_NORMALIZE/MAX(1.E-15, SUM(N_NORMALIZE) )
            print *,'N_NORMALIZE/N_NORMALIZE(1):',N_NORMALIZE/MAX(1.E-15, N_NORMALIZE(1) )
         END SUBROUTINE PRINT_NORMAL


         SUBROUTINE DAY_CYCLE(RDAY_OFFICE, time_mat, T1, T2,   T3, T4) 
! Calculate dayly cycle and output RDAY_OFFICE...
         IMPLICIT NONE
         real RDAY_OFFICE, time_mat, T1, T2,   T3, T4
! Calculate dayly cycle. 
! local variables...
         real tim, one_day

         one_day=24.*3600.0 

         tim=time_mat - one_day*int(time_mat/one_day) 

         if(tim<t1) then
            RDAY_OFFICE=0.0
         else if(tim<t2) then
            rday_office =1.0 -  (t2-tim)/(t2-t1)
         else if(tim<t3) then
            rday_office =1.0 
         else if(tim<t4) then
            rday_office =(t4-tim)/(t4-t3)
         else
            RDAY_OFFICE=0.0
         endif
             
         END SUBROUTINE DAY_CYCLE



         SUBROUTINE ADD_IN_SIGMA(I,J,K, IG,JG, VALUE1,   &
         SIGMA_S_OFF, SIGMA_S_OFF_PTER, FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NX,NY,NZ,NG,NG2, RECORD, SPAR_RECORD)
! Here we add VALUE1 into SIGMA_S_OFF at position I,J,K,  IG, JG
! We have an unstructured options for large numebrs of groups and in this case 
! if RECORD then just record the positions in SPAR_RECORD 
! iF NCOL_SIGMA_S.ne.0 then NG2=0
         IMPLICIT NONE
         INTEGER, intent( in ) :: I,J,K, IG,JG,  NCOL_SIGMA_S, NX,NY,NZ,NG,NG2
         REAL, intent( in ) :: VALUE1
         INTEGER, intent( in ) :: FIN_SIGMA_S(NG+1), COL_SIGMA_S(NCOL_SIGMA_S)
         LOGICAL, intent( in ) :: RECORD
         LOGICAL, intent( inout ) :: SPAR_RECORD(NG,NG) 
         REAL, intent( inout ) :: SIGMA_S_OFF(NX,NY,NZ,NG,NG2), SIGMA_S_OFF_PTER(NX,NY,NZ,NCOL_SIGMA_S)  
! Local variables...
         INTEGER COUNT
         IF(RECORD) THEN
            SPAR_RECORD(IG,JG)=.TRUE.
         ELSE
            IF(NCOL_SIGMA_S.NE.0) THEN
! Find where in the sparse data structure position IG, JG exists. 
               CALL SPAR_FIND(FIN_SIGMA_S, COL_SIGMA_S, NCOL_SIGMA_S, NG, IG, JG, COUNT) 
               SIGMA_S_OFF_PTER(I,J,K, COUNT) = SIGMA_S_OFF_PTER(I,J,K, COUNT) + VALUE1
            ELSE
               SIGMA_S_OFF(I,J,K, IG,JG) = SIGMA_S_OFF(I,J,K, IG,JG) + VALUE1
            ENDIF
         ENDIF 
         RETURN
         END SUBROUTINE ADD_IN_SIGMA




         SUBROUTINE SPAR_FIND( FINDRM, COLM, NCOLM, NONODS, GLOBI,GLOBJ, COUNT) 
! Find where in the sparse data structure position GLOBI,GLOBJ exists. 
         IMPLICIT NONE
         INTEGER, intent( in ) :: NCOLM, NONODS, GLOBI,GLOBJ
         INTEGER, intent( in ) :: FINDRM(nonods+1), COLM(NCOLM)
         INTEGER, intent( inout ) :: COUNT
! local variables...
         LOGICAL FAST
         PARAMETER(FAST=.TRUE.) ! If not fast then use brute force
         INTEGER LOWER, UPPER, INUM, COUNT2
         IF(FAST) THEN ! Perform rapid nested bi-section search. 
                LOWER=FINDRM(GLOBI) 
                UPPER=FINDRM(GLOBI+1)-1
7000            CONTINUE
                INUM=LOWER+(UPPER-LOWER+1)/2 
                IF(GLOBJ.GE.COLM(INUM) )  THEN 
                  LOWER=INUM
                ELSE
                  UPPER=INUM
                ENDIF
                IF(UPPER-LOWER.LE.1) THEN
                  IF(GLOBJ.EQ.COLM(LOWER)) THEN
                      COUNT=LOWER
                  ELSE
                      COUNT=UPPER
                  ENDIF
                  GOTO 9000
                ENDIF
                GOTO 7000
9000            CONTINUE
         ELSE ! Perform brute force search...
            COUNT=0
            DO COUNT2=FINDRM(GLOBI),FINDRM(GLOBI+1)-1
               IF(COLM(COUNT2)==GLOBJ) THEN
                  IF(COUNT.NE.0) STOP 2221
                  COUNT=COUNT2
               ENDIF
            END DO
            IF(COUNT==0) STOP 2982
         ENDIF
         RETURN
         END SUBROUTINE 




