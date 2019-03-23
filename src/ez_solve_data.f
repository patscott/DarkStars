      MODULE ez_solve_data    ! solver data
      USE star_constants
      USE star_data
      IMPLICIT NONE

      INTEGER :: SOLV_ITER, SOLV_ITER_SAV
      ! max allowed iterations for SOLVER.
      DOUBLE PRECISION :: accuracy_target, accuracy_target_SAV
      ! accuracy sought in a converged model
      ! if the root mean squared residual equation error is less than accuracy_target,
      ! then the system considers the model to be good enough.
      DOUBLE PRECISION :: solver_param, solver_param_SAV
      ! if avg_correction is greater than solver_param, then only correct by a fraction solver_param/ERR

      INTEGER, PARAMETER :: N_First_Order = 6      ! number of 1st order equations per shell
      INTEGER, PARAMETER :: N_Second_Order = 5     ! number of 2nd order equations per shell
      INTEGER, PARAMETER :: N_EQs = N_First_Order+N_Second_Order ! total number of 1st and 2nd order equations per shell
      INTEGER, PARAMETER :: N_VBs = N_EQs ! The number of variables per block = number of equations per block
      INTEGER, PARAMETER :: N_BCs = 3     ! number of boundary conditions at either surface or center
      INTEGER, PARAMETER :: N_CNTR_EQs = N_BCs ! the center equation block just has the center BCs
      INTEGER, PARAMETER :: N_SURF_EQs = N_VBs-N_CNTR_EQs
      INTEGER, PARAMETER :: max_N_TOTAL = max_N_SHELLs*N_VBs
      INTEGER, PARAMETER :: JACOB_L = N_SURF_EQs + N_EQs ! The number of (Lower) subdiagonals in the banded Jacobian
      INTEGER, PARAMETER :: JACOB_U = 2*N_EQs + N_CNTR_EQs - 1 ! The number of (Upper) superdiagonals
      INTEGER, PARAMETER :: JACOB_LD = 2*JACOB_L + JACOB_U + 1 ! The Leading Dimension of the Jacobian band array
      
      LOGICAL :: MAKE_NEW_J, MAKE_NEW_J_SAV
      
      INTEGER :: SOLV_JACOBIANS, SOLV_CYCLES, SOLV_CALLS
      INTEGER :: SOLV_JACOBIANS_SAV, SOLV_CYCLES_SAV, SOLV_CALLS_SAV
      
      CONTAINS
      
      SUBROUTINE SAV_solve_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) accuracy_target, solver_param, SOLV_ITER
         READ (IO_UNIT) MAKE_NEW_J, SOLV_JACOBIANS, SOLV_CYCLES, SOLV_CALLS
      ELSE
         WRITE (IO_UNIT) accuracy_target, solver_param, SOLV_ITER
         WRITE (IO_UNIT) MAKE_NEW_J, SOLV_JACOBIANS, SOLV_CYCLES, SOLV_CALLS
      END IF
      END SUBROUTINE SAV_solve_data

      SUBROUTINE SAV_solve_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         accuracy_target=accuracy_target_SAV; solver_param=solver_param_SAV; SOLV_ITER=SOLV_ITER_SAV
         MAKE_NEW_J=MAKE_NEW_J_SAV; SOLV_JACOBIANS=SOLV_JACOBIANS_SAV; SOLV_CYCLES=SOLV_CYCLES_SAV
		 SOLV_CALLS=SOLV_CALLS_SAV
      ELSE
         accuracy_target_SAV=accuracy_target; solver_param_SAV=solver_param; SOLV_ITER_SAV=SOLV_ITER
         MAKE_NEW_J_SAV=MAKE_NEW_J; SOLV_JACOBIANS_SAV=SOLV_JACOBIANS; SOLV_CYCLES_SAV=SOLV_CYCLES
		 SOLV_CALLS_SAV=SOLV_CALLS
      END IF
      END SUBROUTINE SAV_solve_data_internal

      END MODULE ez_solve_data
