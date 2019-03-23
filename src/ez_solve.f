      MODULE ez_solve
      USE star_data
      USE ez_data
      USE ez_solve_data
      USE ez_equations_data
      USE ez_equations
      USE ez_shell
      IMPLICIT NONE
      
      DOUBLE PRECISION :: IPIV(max_N_TOTAL) ! working storage for the matrix solver
      INTEGER :: VB_ARRAY(N_VBs) ! the permutation of the variables used by the solver
      DOUBLE PRECISION :: EQ_RHS(max_N_TOTAL) ! the RHS of the equation
      DOUBLE PRECISION :: EQ_X(max_N_TOTAL) ! the solution of the equation
      DOUBLE PRECISION :: H_DIFF(max_N_TOTAL) ! the variation for each variable in calculating difference approximations to derivatives
      DOUBLE PRECISION :: D_DH(N_VBs,max_N_SHELLs) ! the correction to the increment (which is in DH) to the state (which is in H)
      DOUBLE PRECISION :: VB_UNITs(N_VBs) ! the "natural units" for each variable based on current max(abs(values)).
      DOUBLE PRECISION :: BANDED_JACOBIAN(JACOB_LD,max_N_TOTAL) ! each row of this matrix is a band of the Jacobian matrix
      
      CONTAINS
      
      SUBROUTINE SOLVER ( ITER, corr_target, corr_param, JO )
      USE ez_cycle_data
      INTEGER, INTENT(IN) :: ITER ! max number of iterations -- actually, it is really used as the max number of Jacobian evaluations.
      DOUBLE PRECISION, INTENT(IN) :: corr_target, corr_param
      INTEGER, INTENT(OUT) :: JO ! gets   termination code if fail
      INTEGER :: i, INFO, stop_iter, K_max_corr, I_max_corr
      DOUBLE PRECISION :: avg_corr, prev_corr, max_corr
      DOUBLE PRECISION, PARAMETER :: target_corr_factor = 0.9D0
      INTEGER, PARAMETER :: corr_param_factor = 10
      LOGICAL, PARAMETER :: always_remake_Jacobian = .false.
      VB_ARRAY = (/ V_LNR, V_L, V_M, V_XH, V_XO, V_XHE, V_XC, V_XNE, V_LNF, V_LNT, V_Q_dk /)
      CALL Setup_VB_Units ! get the 'natural' units for the variables
      SOLV_CALLS = SOLV_CALLS + 1
      i = 1 ! count the number of iterations
      stop_iter = ITER ! stop when i reaches this
      
      MAKE_NEW_J = .TRUE.
      
      DO ! search for corrections to the increment that satisfy the equations
         SOLV_CYCLES = SOLV_CYCLES + 1
         CALL Build_Equation
         CALL Solve_Equation(INFO)
         !num_solves = num_solves + 1 !Commented out for DarkStars - Pat Scott Mar 27 2008
         IF (INFO .NE. 0) THEN ! singular matrix or some such disaster
            IF (MAKE_NEW_J) THEN
               EXIT ! we were using a new Jacobian and still got into trouble.  time to give up.
            END IF
            MAKE_NEW_J = .TRUE.; CYCLE ! remake the Jacobian and try again
         END IF
         avg_corr = Make_D_DH(corr_param, max_corr, K_max_corr, I_max_corr)
         ! transform from EQ_X in natural units to D_DH in normal units and return average correction
         IF ( ((avg_corr - avg_corr) .NE. 0D0) .OR. (avg_corr .GT. corr_param_factor*corr_param) ) THEN
            IF (MAKE_NEW_J) THEN
               EXIT ! we were using a new Jacobian and still got into trouble.  time to give up.
            END IF
            MAKE_NEW_J = .TRUE.; CYCLE ! remake the Jacobian and try again
         END IF
         IF (.NOT. MAKE_NEW_J) THEN ! we used an old Jacobian.  make sure it gave good results.
            IF (avg_corr .GT. target_corr_factor * prev_corr) THEN ! not good enough
               MAKE_NEW_J = .TRUE.; CYCLE ! remake the Jacobian and try again
            END IF
            stop_iter = stop_iter+1  ! it was okay to reuse the Jacobian, so increase the number of allowed iterations
         END IF
         IF ( i .GE. stop_iter ) THEN
            EXIT ! have run out of tries
         END IF
         CALL Make_New_DH ! add D_DH to DH and check for consistent composition fractions
         IF ( avg_corr .LE. corr_target ) RETURN
         ! decide if need to recalculate the Jacobian next time
         IF (always_remake_Jacobian) THEN
            MAKE_NEW_J = .TRUE.
         ELSEIF (i .EQ. 1) THEN
            MAKE_NEW_J = (avg_corr .GT. corr_param)
         ELSEIF (.NOT. MAKE_NEW_J) THEN ! currently, limit to one reuse of Jacobian
            MAKE_NEW_J = .TRUE.
         ELSE
            MAKE_NEW_J = (avg_corr .GT. target_corr_factor * prev_corr)
         END IF
         prev_corr = avg_corr
         i = i+1
      END DO
      JO = JO_SOLV ! this termination code says we failed to find a solution
      END SUBROUTINE SOLVER
      
      SUBROUTINE Build_Equation
      IF (MAKE_NEW_J) THEN ! make a new Jacobian
         SOLV_JACOBIANS = SOLV_JACOBIANS + 1
         CALL Build_Banded_Jacobian ! this creates the Jacobian and also sets up the RHS of the equation
      ELSE ! reuse the existing Jacobian
         CALL Eval_EQs ! this sets up the RHS of the equation, but leaves the Jacobian unchanged
      END IF
      END SUBROUTINE Build_Equation

      
      SUBROUTINE Solve_Equation(INFO)
      INTEGER, INTENT(OUT) :: INFO
      INTEGER :: N_TOTAL
      N_TOTAL = N_SHELLs*N_VBs
      EQ_X = -EQ_RHS
      IF (MAKE_NEW_J) THEN
         CALL DGBTRF( N_TOTAL, N_TOTAL, JACOB_L, JACOB_U, BANDED_JACOBIAN, JACOB_LD, IPIV, INFO )
         IF (INFO /= 0) RETURN
      END IF
      CALL DGBTRS( 'No transpose', N_TOTAL, JACOB_L, JACOB_U, 1, BANDED_JACOBIAN, JACOB_LD, IPIV, EQ_X, N_TOTAL, INFO )      
      END SUBROUTINE Solve_Equation
      
      
      SUBROUTINE Setup_VB_Units
      INTEGER :: IJ, I, K
      DOUBLE PRECISION :: MX, TMP
      ! corrections for variable IJ are calculated as a fraction of VB_UNITs(IJ)
      DO IJ = 1, N_VBs ! set VB_UNITs array to max(abs(values)) for the variables
         I = VB_ARRAY(IJ) ! I is the variable number
         MX = 1D0 ! this is a lower limit for the natural units so never scale-up, only scale-down.
         DO K = 1, N_SHELLs
            TMP = abs(H(I,K))
            IF (TMP  .GT.  MX) MX = TMP
         END DO
         VB_UNITs(I) = MX
      END  DO
      END SUBROUTINE Setup_VB_Units
      
      SUBROUTINE Build_Banded_Jacobian
      INTEGER :: N_TOTAL
      N_TOTAL = N_SHELLs*N_VBs
      BANDED_JACOBIAN(:,1:N_TOTAL) = 0.0D0
      !num_jacobians = num_jacobians + 1 !Commented out for DarkStars - Pat Scott Mar 27 2008
      CALL Do_Diff_EQs(.TRUE.)
      END SUBROUTINE Build_Banded_Jacobian
      
      SUBROUTINE Eval_EQs ! using the current variable values, evaluate and save the equations.
      CALL Do_Diff_EQs(.FALSE.)
      END SUBROUTINE Eval_EQs
      
      SUBROUTINE Do_Diff_EQs(JACOBIAN_FLG)
      LOGICAL, INTENT(IN) :: JACOBIAN_FLG
      INTEGER :: Q
      CALL DIFF_FNs(N_CNTR_shell, JACOBIAN_FLG)
      CALL DIFF_EQs(N_CNTR_shell+1, N_CNTR_EQs, JACOBIAN_FLG) 
      DO Q = N_CNTR_shell, N_SURF_shell+1, -1
         CALL DIFF_FNs(Q-1, JACOBIAN_FLG)
         CALL DIFF_EQs(Q, N_EQs, JACOBIAN_FLG)
      END DO
      CALL DIFF_EQs(N_SURF_shell, N_SURF_EQs, JACOBIAN_FLG) 
      END SUBROUTINE Do_Diff_EQs
      
      SUBROUTINE Make_New_DH
      INTEGER :: IJ, IK, I
      INTEGER, PARAMETER :: NV=5
      INTEGER, PARAMETER :: IX(NV)  = (/ V_XH, V_XHE, V_XC, V_XO, V_XNE /)
      DOUBLE PRECISION, PARAMETER :: X_LIM = 1D-12
      DOUBLE PRECISION :: prev_M, prev_lnR, M, lnR, diff, x
      LOGICAL :: made_change
      INTEGER :: cnt
      DH = DH + D_DH ! add correction vector (D_DH) to the increment vector (DH)
      ! make sure new composition fractions are non-negative
      DO I = 1, NV
         IJ = IX(I)
         DO IK = 1, N_SHELLs
            IF ( H(IJ,IK) + DH(IJ,IK) .LE. X_LIM ) THEN 
               H(IJ,IK) = 0D0; DH(IJ,IK) = 0D0
            END IF
         END DO
      END DO
      ! make sure mass is decreasing inward
      made_change = .TRUE.; cnt = 0
      DO WHILE (made_change .AND. cnt .LT. 10)
         made_change = .FALSE.
         prev_M = H(V_M,1) + DH(V_M,1)
         DO IK = 2, N_SHELLs
            M = H(V_M, IK) + DH(V_M,IK)
            diff = prev_M - M
            IF ( diff .LT. 0D0 ) THEN
               IF ( diff .GT. -1D-15 ) THEN
                  !WRITE(*,*) 'model', model_Number, 'merge at IK =', IK, 'prev_M', prev_M, 'M', M, 'diff', diff
                  H(V_M, IK) = prev_M
                  DH(V_M, IK) = 0D0
               ELSE
                  !WRITE(*,*) 'model', model_Number, 'swap at IK =', IK, 'prev_M', prev_M, 'M', M, 'diff', diff
                  H(V_M, IK) = prev_M; H(V_M, IK-1) = M
                  DH(V_M, IK) = 0D0; DH(V_M,IK-1) = 0D0
                  DO I = 1, NV ! swap abundances too
                     IJ = IX(I)
                     x = H(IJ, IK); H(IJ, IK) = H(IJ, IK-1); H(IJ, IK-1) = x
                  END DO
               END IF
               made_change = .TRUE.
            END IF
            prev_M = M
         END DO
         cnt = cnt + 1
      END DO
      return ! ignore lnR for now
      made_change = .TRUE.; cnt = 0
      DO WHILE (made_change .AND. cnt .LT. 10)
         made_change = .FALSE.
         prev_lnR = H(V_LNR,1) + DH(V_LNR,1)
         DO IK = 2, N_SHELLs
            lnR = H(V_LNR, IK) + DH(V_LNR,IK)
            diff = prev_lnR - lnR
            IF ( diff .LE. 0D0 ) THEN
               H(V_LNR, IK) = prev_lnR; H(V_LNR, IK-1) = lnR
               DH(V_LNR, IK) = 0D0; DH(V_LNR,IK-1) = 0D0
               made_change = .TRUE.
            END IF
            prev_lnR = lnR
         END DO
         cnt = cnt + 1
      END DO
      END SUBROUTINE Make_New_DH
      
      DOUBLE PRECISION FUNCTION Make_D_DH(corr_param, max_corr, K_max_corr, I_max_corr)
      ! D_DH is the correction vector to use in changing the increment DH
      DOUBLE PRECISION, INTENT(IN) :: corr_param
      DOUBLE PRECISION, INTENT(OUT) :: max_corr
      INTEGER, INTENT(OUT) :: K_max_corr, I_max_corr
      DOUBLE PRECISION :: avg_corr, correction_factor, TMP
      INTEGER :: IJ, K, I, K_offset
      INTEGER :: N_TOTAL
      N_TOTAL = N_SHELLs*N_VBs
      avg_corr = 0D0; max_corr = 0D0
      DO I = 1, N_VBs
         K_offset = 0
         DO K = 1, N_SHELLs
            TMP = EQ_X(K_offset+I)
            IF (TMP .LT. 0D0) TMP = -TMP
            IF (TMP .GT. max_corr) THEN
               max_corr = TMP; K_max_corr = K; I_max_corr = VB_ARRAY(I)
            END IF
            avg_corr = avg_corr + TMP
            K_offset = K_offset + N_VBs
         END DO
      END DO
      avg_corr = avg_corr/N_TOTAL ! average abs(correction)
      correction_factor = corr_param/max(avg_corr, corr_param) ! reduce corrections by this factor
      DO I = 1, N_VBs
         IJ = VB_ARRAY(I) ! IJ is the variable number
         TMP = correction_factor*VB_UNITs(IJ)
         K_offset = 0
         DO K = 1, N_SHELLs
            D_DH(IJ, K) = EQ_X(K_offset+I)*TMP ! the correction vector for this variable
            K_offset = K_offset + N_VBs
         END DO
      END DO
      Make_D_DH = avg_corr
      END FUNCTION Make_D_DH
      
      SUBROUTINE Get_FNs(K)
      INTEGER, INTENT(IN) :: K
      DVAR(1:N_VBs) = DH(1:N_VBs, K)
      VAR(1:N_VBs) = H(1:N_VBs, K) + DVAR(1:N_VBs)
      CALL FUNCS1 ( K, 0 ) ! the 0 means cache the values
      END SUBROUTINE Get_FNs
      
      SUBROUTINE DIFF_FNs( K, JACOBIAN_FLG )
      INTEGER, INTENT(IN) :: K
      LOGICAL, INTENT(IN) :: JACOBIAN_FLG
      DOUBLE PRECISION DX, DVX, VX
      DOUBLE PRECISION, PARAMETER :: DH0 = 1.0D-07
      INTEGER :: JI, I, V_offset
      CALL Get_FNs(K)
      V_offset = (K-1)*N_VBs
      IF (.NOT. JACOBIAN_FLG) RETURN
      DO I = 1, N_VBs ! modify the variables in turn and recalculate the functions for use below in calculating modified equations
         JI = VB_ARRAY(I) ! the variable index
         VX = VAR(JI) ! variable as it was last time
         DX = DH0*max(abs(VX), 1D0) ! DH0 is variation size for numerical partials (1e-7)
         DVX = DVAR(JI) ! increment for the variable for next time
         IF ( DVX .LT. 0D0 ) THEN
            DX = -DX ! 'correction' to increment has same sign as current increment
         END IF
         H_DIFF(V_offset+I) = VB_UNITs(JI)/DX
            ! this is 1 divided by the change in the scaled variable.  used for calculating derivatives wrt the variable.
            ! DX is the change in the unscaled variable.  VB_UNITs(JI) is the scale factor for this variable.
            ! the scaled variable is equal to the unscaled variable divided by the scale factor.
         VAR(JI) = VX+DX
         DVAR(JI) = DVX+DX
         CALL FUNCS1 ( K, JI ) ! the JI arg to FUNCS1 lets it know which variable we've changed
         ! restore the variable
         DVAR(JI) = DVX
         VAR(JI) = VX
      END DO
      END SUBROUTINE DIFF_FNs
      
      SUBROUTINE DIFF_EQs( Q, IQ, JACOBIAN_FLG )
      ! do the first IQ equations of equation block number Q.  assumes that the function values are already in place.
      USE ez_equations
      INTEGER, INTENT(IN) :: Q, IQ
      LOGICAL, INTENT(IN) :: JACOBIAN_FLG
      INTEGER, PARAMETER :: IC=3, IP=2, IPP=1 ! IC for current, IP for previous, IPP for past-previous.
      INTEGER :: IEE, I, V_num, Q_offset, ji, ji1, ji2
      LOGICAL :: do_IPP
      CALL EQUNS1 ( Q, 0, 0 ) ! get the unmodified equation values for equation block Q
      Q_offset = max(0, N_SURF_EQs + (Q-2)*N_EQs)
      EQ_RHS(Q_offset+1:Q_offset+IQ) = EQU(1:IQ)
      IF (.NOT. JACOBIAN_FLG) RETURN
      do_IPP = (Q .GT. N_SURF_shell) .AND. (Q .LT. N_CNTR_shell)
      DO IEE = IPP, IC ! each block of equations depends on (up to) 3 shells of variables.  IEE specifies which shell we're doing.
         IF ( (IEE .EQ. IC) .OR. (IEE .EQ. IP .AND. Q .NE. N_SHELLs+1) .OR. (IEE .EQ. IPP .AND. do_IPP) ) THEN
            DO I = 1, N_VBs ! I is the variable in the shell
               CALL EQUNS1 ( Q, VB_ARRAY(I), IEE ) ! get the modified equation values
               IF (Q .EQ. 1) THEN ! surface block is special
                  V_num = I + N_VBs*(IC-IEE)
               ELSE
                  V_num = I + N_VBs*(Q+IPP-IEE)
               END IF
               ji = JACOB_L+JACOB_U+1-V_num
               ji1 = ji+Q_offset+1
               ji2 = ji+Q_offset+IQ
               BANDED_JACOBIAN(ji1:ji2, V_num) = (EQU(1:IQ) - EQ_RHS(Q_offset+1:Q_offset+IQ))*H_DIFF(V_num)
            END DO
         END IF
      END DO
      END SUBROUTINE DIFF_EQs
      
      END MODULE ez_solve
      
      
      
      
