      MODULE ez_cycle_data
						! internal save routines and variables added by Pat Scott 24-03-08
      USE star_data
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: N_XTRAS = 50
      DOUBLE PRECISION :: HPR(V_MAX,max_N_SHELLs), HPR_SAV(V_MAX,max_N_SHELLs)     ! holds previous H for backup.
      DOUBLE PRECISION :: HPR2(V_MAX,max_N_SHELLs), HPR2_SAV(V_MAX,max_N_SHELLs)   ! holds previous HPR for backup.
      DOUBLE PRECISION :: XTRAS_PR(N_XTRAS), XTRAS_PR_SAV(N_XTRAS)				   ! backup for XTRAS
      DOUBLE PRECISION :: XTRAS_PR2(N_XTRAS), XTRAS_PR2_SAV(N_XTRAS)               ! backup for XTRAS_PR

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      INTEGER :: solver_max_iter1, solver_max_iter1_SAV
      ! maximum number of equation solver iterations allowed in the first solver_iter_startup_num timesteps for a new model
      INTEGER :: solver_max_iter2, solver_max_iter2_SAV
      ! maximum number of equation solver iterations allowed after the first solver_iter_startup_num timesteps
      INTEGER :: solver_max_iter3, solver_max_iter3_SAV
      ! maximum number of equation solver iterations allowed after while making artificial changes to star mass or energy
      INTEGER :: solver_iter_startup_num, solver_iter_startup_num_SAV
      ! number of timesteps to use solver_max_iter1 before changing to solver_max_iter2
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      DOUBLE PRECISION :: DT, DT_SAV                           ! the timestep in seconds.
      DOUBLE PRECISION :: DTY, DTY_SAV                   ! timestep in years
      DOUBLE PRECISION :: PREV_DT
      DOUBLE PRECISION :: REQUESTED_DTY ! in case somebody needs a smaller timestep

      INTEGER, PARAMETER :: TV_CDD=1, TV_AGE=2, TV_MC=3, TV_SM=4, TV_SE=5, TV_TN=6, TVNUM=6
      DOUBLE PRECISION :: AGE, AGE_SAV                   ! stellar age in years.
      DOUBLE PRECISION :: CDD, CDD_SAV                   ! target value for average change in variables; used in adjusting timestep.
      DOUBLE PRECISION :: MC, MC_SAV					 ! central mass in 10^33 grams.
      DOUBLE PRECISION :: SM, SM_SAV                     ! total mass in Msolar.
      DOUBLE PRECISION :: PR1(TVNUM), PR1_SAV(TVNUM)     ! copy of the preceeding for the previous model.
      DOUBLE PRECISION :: PR2(TVNUM), PR2_SAV(TVNUM)     ! copy of PR1 for the previous model.

      INTEGER :: JHOLD, JHOLD_SAV            ! the number of timesteps taken since last BACKUP.
      INTEGER :: JMOD, JMOD_SAV				 ! current model number.
      INTEGER :: JM1, JM1_SAV                ! previous model number.
      INTEGER :: JM2, JM2_SAV                ! previous JM1.
      INTEGER :: JNN, JNN_SAV                ! number of timesteps taken.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      DOUBLE PRECISION :: CDC1, CDC1_SAV     ! value for CDD from ZAMS to He ignition.
      DOUBLE PRECISION :: CDC2, CDC2_SAV     ! extra factor for CDD during He burning.
      DOUBLE PRECISION :: CDC3, CDC3_SAV     ! extra factor for CDD until He shell near H shell.
      DOUBLE PRECISION :: CDC4, CDC4_SAV     ! extra factor for CDD during double shell burning.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! The control codes (which the program puts in the variable JO).
      INTEGER, PARAMETER :: JO_NULL = -1   ! default value while working on model.
      INTEGER, PARAMETER :: JO_AOK  = 0      ! finished okay.
      INTEGER, PARAMETER :: JO_SOLV = 1     ! equation solver ran into trouble; backup and try again.
      INTEGER, PARAMETER :: JO_TSTEP   = 2    ! time step reduced below limit; cannot backup and cannot go forward.
      INTEGER, PARAMETER :: JO_START   = 3    ! attempted backup for first or second model.
      INTEGER, PARAMETER :: JO_BACK = 4   ! request for backup.
      INTEGER, PARAMETER :: JO_PSI = 5    ! degeneracy too large for equation of state code.
      INTEGER, PARAMETER :: JO_GAM = 6    ! plasma interaction parameter too large for equation of state code.
     
      CONTAINS
      
      SUBROUTINE SAV_cycle_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) AGE, CDD, MC, SM, DTY, PR1, PR2, JHOLD, JMOD, JM1, JM2, JNN, CDC1, CDC2, CDC3, CDC4
         READ (IO_UNIT) DT, HPR, HPR2, XTRAS_PR, XTRAS_PR2
         READ (IO_UNIT) solver_max_iter1, solver_max_iter2, solver_max_iter3, solver_iter_startup_num
      ELSE
         WRITE (IO_UNIT) AGE, CDD, MC, SM, DTY, PR1, PR2, JHOLD, JMOD, JM1, JM2, JNN, CDC1, CDC2, CDC3, CDC4
         WRITE (IO_UNIT) DT, HPR, HPR2, XTRAS_PR, XTRAS_PR2
         WRITE (IO_UNIT) solver_max_iter1, solver_max_iter2, solver_max_iter3, solver_iter_startup_num
      END IF
      END SUBROUTINE SAV_cycle_data

      SUBROUTINE SAV_cycle_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         AGE=AGE_SAV; CDD=CDD_SAV; MC=MC_SAV; SM=SM_SAV; DTY=DTY_SAV; PR1=PR1_SAV; PR2=PR2_SAV; JHOLD=JHOLD_SAV
		 JMOD=JMOD_SAV; JM1=JM1_SAV; JM2=JM2_SAV; JNN=JNN_SAV; CDC1=CDC1_SAV; CDC2=CDC2_SAV; CDC3=CDC3_SAV; CDC4=CDC4_SAV
         DT=DT_SAV; HPR=HPR_SAV; HPR2=HPR2_SAV; XTRAS_PR=XTRAS_PR_SAV; XTRAS_PR2=XTRAS_PR2_SAV
         solver_max_iter1=solver_max_iter1_SAV; solver_max_iter2=solver_max_iter2_SAV 
		 solver_max_iter3=solver_max_iter3_SAV; solver_iter_startup_num=solver_iter_startup_num_SAV
      ELSE
         AGE_SAV=AGE; CDD_SAV=CDD; MC_SAV=MC; SM_SAV=SM; DTY_SAV=DTY; PR1_SAV=PR1; PR2_SAV=PR2; JHOLD_SAV=JHOLD
		 JMOD_SAV=JMOD; JM1_SAV=JM1; JM2_SAV=JM2; JNN_SAV=JNN; CDC1_SAV=CDC1; CDC2_SAV=CDC2; CDC3_SAV=CDC3; CDC4_SAV=CDC4
         DT_SAV=DT; HPR_SAV=HPR; HPR2_SAV=HPR2; XTRAS_PR_SAV=XTRAS_PR; XTRAS_PR2_SAV=XTRAS_PR2
         solver_max_iter1_SAV=solver_max_iter1; solver_max_iter2_SAV=solver_max_iter2 
		 solver_max_iter3_SAV=solver_max_iter3; solver_iter_startup_num_SAV=solver_iter_startup_num_SAV
      END IF
      END SUBROUTINE SAV_cycle_data_internal

      END MODULE ez_cycle_data

