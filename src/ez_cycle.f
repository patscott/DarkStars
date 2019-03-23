      MODULE ez_cycle
      USE ez_data
      USE ez_cycle_data
	  
      !Added for DarkStars by Pat Scott 23-03-08
      USE DkStrs_utils
      USE DkStrs_data
	  
      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXMODS = 10000000
      
      CONTAINS
      
      INTEGER FUNCTION Null_Check()
      USE star_controls
      Null_Check = KEEP_GOING
      END FUNCTION Null_Check
      
      SUBROUTINE Evolve ( JO, Check_Model )
      INTEGER, INTENT(OUT) :: JO
      INTERFACE
         INTEGER FUNCTION Check_Model()
         END FUNCTION Check_Model
      END INTERFACE
      CALL Evolve1 ( JO, MAXMODS, Check_Model, Null_Check )
      END SUBROUTINE Evolve
      
      SUBROUTINE Evolve1 ( JO, KP, Check_Model1, Check_Model2 )
      USE star_data
      USE ez_report
      USE ez_solve
      USE ez_solve_data
      USE star_controls
      INTEGER, INTENT(OUT) :: JO
      INTEGER, INTENT(IN) :: KP
      INTERFACE
         INTEGER FUNCTION Check_Model1()
         END FUNCTION Check_Model1
      END INTERFACE
      INTERFACE
         INTEGER FUNCTION Check_Model2()
         END FUNCTION Check_Model2
      END INTERFACE
      INTEGER :: IBKUP, CHECK1, CHECK2
      LOGICAL :: RESTARTING
      DOUBLE PRECISION :: corr_target, corr_param
      SOLV_JACOBIANS = 0; SOLV_CYCLES = 0; SOLV_CALLS = 0
      MAKE_NEW_J = .TRUE.
      IBKUP = 0 ! IBKUP is set to 1 if discover that require a backup
      REQUESTED_DTY = 0D0 ! to indicate no request
      corr_target = accuracy_target; corr_param = solver_param
      IF (JNN .GT. 0) THEN
         RESTARTING = .TRUE.
      ELSE
         IF (Check_Model1() .NE. KEEP_GOING .OR. Check_Model2() .NE. KEEP_GOING) RETURN
         CALL NEXTDT
         SOLV_ITER = solver_max_iter1
         JNN = 1
         RESTARTING = .FALSE.
      END IF
      DO WHILE ( JNN .LE. KP ) ! evolutionary loop of KP time steps
         IF (.NOT. RESTARTING) THEN
            IF ( IBKUP .EQ. 0 ) THEN
               JO = JO_AOK
               IF ( extra_Mdot_param .NE. 0D0 .OR. extra_Energy_param .NE. 0D0 ) THEN
                  SOLV_ITER = solver_max_iter3
               END IF
               CALL SOLVER ( SOLV_ITER, corr_target, corr_param, JO )
            END IF
            IF ( IBKUP .NE. 0 .OR. JO .NE. JO_AOK ) THEN ! backup and try again
               IBKUP = 0; MAKE_NEW_J = .TRUE.
               IF ( JNN .EQ. 1 ) THEN
                  CALL BACKUP1 ( JO ) ! restart from the beginning
               ELSE
                  CALL BACKUP2 ( JO ) ! restart from 2 steps back with DT decreased substantially
               END IF
               IF ( JO .GE. JO_TSTEP ) EXIT ! cannot backup because timestep already as small as can get
               CYCLE;
            END IF
            CALL Complete_Model
            CALL Check_Final ( JO )
            IF ( JO .EQ. JO_BACK ) THEN
               ! reduce timestep and try again
               IBKUP = 1; CYCLE  ! need to backup.
            END IF
            CALL UPDATE
            IF ( JO .EQ. JO_AOK ) THEN
               CHECK1 = Check_Model1()
               IF (CHECK1 .EQ. TERMINATE) RETURN
               CHECK2 = Check_Model2()
               IF (CHECK2 .EQ. TERMINATE) RETURN
               IF (CHECK1 .EQ. BACKUP .OR. CHECK2 .EQ. BACKUP) THEN
                  IBKUP = 1; CYCLE  ! need to backup.
               END IF
            END IF
            IF ( JO .GE. JO_TSTEP ) EXIT
         ELSE ! skip the stuff above first time after restarting
            RESTARTING = .FALSE.
         END IF
         CALL NEXTDT
         IF (JNN .GE. solver_iter_startup_num) THEN
            SOLV_ITER = solver_max_iter2
         END IF
         corr_target = accuracy_target; corr_param = solver_param
         JNN = JNN + 1
      END DO ! end of evolutionary loop
      END SUBROUTINE Evolve1

      SUBROUTINE Check_Final ( JO )
      USE star_extras
      USE ez_solve_data
      INTEGER, INTENT(OUT) :: JO
      CDD = CDC1 ! target value for SUM used in adjusting DT
      IF ( mass_C_Core .LE. mass_He_Core ) THEN ! sensible results for H1 & He4 abundances
         ! adjust CDD to allow larger or smaller time-increments, according to various rules-of-thumb
         IF ( center_H .LT. 1.0D-5 .AND. center_He .LT. 0.95D0 ) THEN
            ! almost no central H1 and less than 95% central He4
               CDD = CDC1*CDC2 ! extra factor for CDD during helium burning
         END IF
         IF ( center_He .LT. 0.1D0 ) THEN
            ! central He4 < 10% so has burned; He burning shell has moved out from center
                CDD = CDC1*CDC3 ! extra factor for CDD until He shell near H shell
         END IF
         IF ( mass_C_Core .GT. 0.75D0*mass_He_Core ) THEN
            ! He burning shell catching up with H burning shell
                CDD = CDC1*CDC4 ! extra factor for CDD during double shell burning
         END IF
      END IF
      IF ( center_He .LT. 0.15D0 ) accuracy_target = 1D-4 ! if central He4 abundance drops below 0.15,
         ! reset desired accuracy for SOLVER 1e-4 compared to typical accuracy_target=1e-6
      IF ( center_Degeneracy .GE. PSI_limit ) THEN
         JO = JO_PSI ! degeneracy too large for equation of state code
         WRITE(*,'(a)') ' Stopping because center degeneracy too large for equation of state code to continue.'
      ELSEIF ( center_GAM .GT. GAM_limit ) THEN
         JO = JO_GAM
         WRITE(*,'(a)') ' Stopping because center plasma interaction parameter (Gamma) exceeds limit.'
      ELSEIF ( DT .LT. dynamic_Timescale ) THEN
         JO = JO_BACK
      END IF
      END SUBROUTINE Check_Final

      SUBROUTINE NEXTDT
      USE star_data
      USE star_controls
	  DOUBLE PRECISION FAC, SUM, CT1, CT2
      INTEGER IK, IJ
      INTEGER, PARAMETER :: KN=6
      INTEGER :: KJN(KN) = (/ V_LNF, V_LNT, V_XO, V_XH, V_Q_dk, V_LNR /)
      SUM = 0D0 ! sum of modulus of increments of selected variables
      DO IK = 1, N_SHELLs ! sum for each mesh point
         DO IJ = 1, KN ! sum KN of the increments; KJN array tells which variables to include.
            SUM = SUM + DABS(DH(KJN(IJ), IK))
         END DO
      END DO
      IF ( SUM.EQ.0D0 ) THEN
         FAC = 1D0 
      ELSE
         FAC = SUM/(CDD*N_SHELLs) ! factor to make average meshpoint sum (SUM/N_SHELLs) equal to CDD
      END IF
      ! FAC > 1 means average was bigger than desired, so should reduce timestep.
      ! FAC < 1 means average was smaller than desired, so should increase timestep.
      ! Specify next DTY
      IF ( SUM .NE. 0D0 .AND. JHOLD .GT. timestep_hold ) THEN
         ! JHOLD is number of timesteps taken since last BACKUP
         ! want to adjust timestep to make SUM/N_SHELLs == CDD
         ! but don't change timestep by factor smaller than timestep_lower_limit or bigger than timestep_upper_limit 
         ! and don't change timestep at all for awhile after BACKUP
         DTY = DMAX1( timestep_lower_limit, DMIN1(timestep_upper_limit, 1D0/FAC) ) * DTY
      END IF
      IF ( REQUESTED_DTY .eq. 0 .AND. timestep_max > 0) REQUESTED_DTY = timestep_max
      IF ( REQUESTED_DTY .GT. 0 .AND. REQUESTED_DTY .LT. DTY ) THEN
         DTY = min(DTY, MAX(REQUESTED_DTY, 2D0*(dynamic_Timescale/CSY)))
      END IF
      REQUESTED_DTY = 0D0
      PREV_DT = DT
      DT = CSY*DTY ! timestep in seconds
      CALL setWIMPPop(DTY) ! Added for DarkStars by Pat Scott 23-03-08
      CT1 = timestep_lower_limit; CT2 = timestep_upper_limit; 
      IF (.NOT.(CT1 .EQ. 1D0 .AND. CT2 .EQ. 1D0) .AND. .NOT.(JHOLD .GT. timestep_hold .AND. FAC .LT. 2D0 .AND. JNN .GT. 2)) THEN
         DH = 0D0
      END IF
      END SUBROUTINE NEXTDT
      
      SUBROUTINE UPDATE
      ! Store certain current and previous values, for possible emergency backup
      ! updates AGE and JMOD (number of timesteps)
      ! increments JHOLD, the number of timesteps taken since last BACKUP
      AGE = AGE + DTY ! update model age
      JMOD = JMOD + 1 ! update number of timesteps taken
      CALL Save_State(.TRUE.)
      JHOLD = JHOLD + 1 
      END SUBROUTINE UPDATE
      
      SUBROUTINE Save_State (update)
      USE star_data
      ! save state for backup
      LOGICAL, INTENT(IN) :: update
      DOUBLE PRECISION :: DT_RATIO
      PR2 = PR1
      PR1(TV_AGE)=AGE; PR1(TV_MC)=MC; PR1(TV_SM)=SM; PR1(TV_CDD)=CDD; PR1(TV_TN)=nuc_Timescale
      JM2 = JM1; JM1 = JMOD
	  n_WIMPs_saved = n_WIMPs_prev; n_WIMPs_prev = n_WIMPs !Added for DarkStars by Pat Scott 24-03-08
      HPR2 = HPR; HPR = H
      XTRAS_PR2 = XTRAS_PR; XTRAS_PR = XTRAS
      IF (update) THEN
         DT_RATIO = DT / PREV_DT
         !DH = DH * DT_RATIO    Dont adjust DH until have at least quadratic prediction
            ! for linear prediction, actually do better without this!
         H = H + DH
      END IF
      END SUBROUTINE Save_State
      
      SUBROUTINE BACKUP2 ( JO )
      INTEGER, INTENT(OUT) :: JO
      CALL Do_Backup (JO, .FALSE.)
      END SUBROUTINE BACKUP2

      SUBROUTINE BACKUP1 ( JO )
      INTEGER, INTENT(OUT) :: JO
      CALL Do_Backup (JO, .TRUE.)
      END SUBROUTINE BACKUP1

      SUBROUTINE Do_Backup ( JO, first_step ) ! backup to previous model and reduce timestep (DT).
      USE star_data
	  INTEGER, INTENT(OUT) :: JO
      LOGICAL, INTENT(IN) :: first_step
      IF (first_step) THEN
         AGE=PR1(TV_AGE); MC=PR1(TV_MC); SM=PR1(TV_SM)
         CDD=PR1(TV_CDD); nuc_Timescale=PR1(TV_TN)
         IF (JMOD .NE. JM1) n_WIMPs_prev = n_WIMPs_saved ! Added for DarkStars by Pat Scott 23-03-08
		 JMOD = JM1 ! backup to saved JMOD model number
	  ELSE
         IF (JMOD .NE. JM2) n_WIMPs_prev = n_WIMPs_saved ! Added for DarkStars by Pat Scott 23-03-08
		 JMOD = JM2 ! backup to saved JMOD model number
         AGE=PR2(TV_AGE); MC=PR2(TV_MC); SM=PR2(TV_SM)
         CDD=PR2(TV_CDD); nuc_Timescale=PR2(TV_TN)
      END IF
      PREV_DT = DT
      DT = timestep_decrement*DT  ! reduce timestep by factor of timestep_decrement; 
	  DTY = DT/CSY
      CALL setWIMPPop(DTY) ! Added for DarkStars by Pat Scott 23-03-08
      ! Modified for DarkStars by Pat Scott 2007-07-31; allow to backup to first or second model 
      IF ( .NOT. first_step .AND. JMOD .LE. 0 ) THEN
          JO = JO_START
          write(*,*)
          write(*,*) 'Initial model would not converge; initial timestep may be too small (really!)'
      !IF ( .NOT. first_step .AND. JMOD .LE. 2 ) THEN
      !   JO = JO_START
      ELSE IF ( DT .LT. dynamic_Timescale ) THEN
         JO = JO_TSTEP ! time step too large compared to dynamic timescale
         WRITE(*,'(a)') ' Stopping because required time step smaller than dynamic timescale.'
         IF ( DT .GT. CSY ) THEN
            WRITE(*,'(2(A,E9.3),A)') '  Timestep of ', DT/CSY, ' years vs dynamic timescale of ', dynamic_Timescale/CSY, ' years.'
         ELSE
            WRITE(*,'(2(A,E9.3),A)') '  Timestep of ', DT, ' seconds vs dynamic timescale of ', dynamic_Timescale, ' seconds.'
         END IF
         RETURN
      ELSE
         JO = JO_AOK
      END IF
      ! HPR holds the previous model
      DH = 0D0
      H = HPR; XTRAS = XTRAS_PR
      IF ( .NOT. first_step ) JHOLD = -1
      END SUBROUTINE Do_Backup
      
      SUBROUTINE Reset_Age_and_Model_Number(new_age,model)
      DOUBLE PRECISION, INTENT(IN) :: new_age
      INTEGER, INTENT(IN) :: model
      ! this is needed for the Post_HE_FLASH stuff.  need to reset saved info since might backup to it.
      AGE = new_age; PR1(TV_AGE) = new_age - 1D2; PR2(TV_AGE) = new_age - 2D2
      JMOD = model; JM1 = model; JM2 = model
      END SUBROUTINE Reset_Age_and_Model_Number

      END MODULE ez_cycle
