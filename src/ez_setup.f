      MODULE ez_setup
      USE ez_data
      USE star_data
      USE ez_vcool
      USE ez_opacity
      USE ez_nuclear_data
      USE ez_state_data
      USE ez_cycle_data
      USE ez_shell_data
      USE ez_solve_data
      USE ez_ionization_data
      USE ez_do_one_data
      USE DkStrs_data
      IMPLICIT NONE
      
      DOUBLE PRECISION :: TradeXY_target
      DOUBLE PRECISION :: TradeXY_previous
      INTEGER :: TradeXY_model
      DOUBLE PRECISION :: tradeXY_maxD=0.005D0


      INTEGER, PARAMETER :: file_N_SHELLs = 199 ! the initial models
      

      CONTAINS
      
      SUBROUTINE SETUP
      DOUBLE PRECISION CRHO, LAMC, LAMC3, CHE, EVOLT
      INTEGER JR, J, JT, I, K
      CHARACTER (LEN=strlen) :: fname
      DOUBLE PRECISION :: Qvals(20), Cvals(4)
      integer :: nz, nt, ng, n
      
      ! mesh params
      N_SHELLs = 199
      N_CNTR_shell = N_SHELLs

      ! control params
      T_CRIT = 6D6
      PSI_CRIT = 0D0 
      PRESSI_CONST = 0D0
      solver_param = 1D-2
      solver_max_iter1 = 20  
      solver_max_iter2 = 4
      solver_max_iter3 = 12
      solver_iter_startup_num = 4
      !Timestep_rescale added for DarkStars by Pat Scott somewhere in the mists of time around mid 2008.
      CDC1=2.75D-2*timestep_rescale ! value for CDD from ZAMS to He ignition.
      CDC2=.25D0 ! extra factor for CDD during He burning.
      CDC3=1D0 ! extra factor for CDD until He shell near H shell.
      CDC4=4D0 ! extra factor for CDD during double shell burning.
      ! mathematical constants
      CPI = 4D0*ATAN(1D0)  ! pi
      C4PI = 4D0*CPI       ! 4 pi
      CLN = LOG(10D0)         ! divide by CLN to convert ln to log10
      C3RD = 1D0/3D0    ! 1/3
      ! derived physical constants
      EVOLT = ECHAR/CL*1.0D8  ! (ergs/eV)
      CRAD = BOLTZM/(CL*PLANCK)
      CRAD = 8D0*CPI**5*BOLTZM*CRAD**3/15D0 ! a, the radiation constant
      CME = 1.0D6*EVOLT/AMU   ! (ergs/MeV)/(grams/amu)
      CEVB = EVOLT/BOLTZM     ! eV/kb
      CR = BOLTZM/AMU         ! kb/(grams/amu)
      CHE = AME*CL**2         ! me c^2
      CTE = BOLTZM/CHE        ! kb / (me c^2)
      LAMC = PLANCK/(AME*CL)  ! lambdaC, the Compton electron wavelength, h/(me c)
      LAMC3 = LAMC**3         ! lambdaC^3
      CRHO = 8D0*CPI/LAMC3    ! 8*pi / lambdaC^3
      CB = CRHO*CHE           ! 8 pi me c^2 / lambdaC^3
      CD = CRHO*AMU           ! (8 pi / lambdaC^3) * (grams/amu)
      CEN = PLANCK/SQRT(2D0*CPI*AMU*BOLTZM)
      CEN = CEN**3/AMU        ! h^3/(2 pi kb)^(3/2)/(grams/amu)^(5/2)
      CPL = SQRT(C4PI/(AMU*BOLTZM**3))*ECHAR**3 ! (4*pi*e^3)/(kb^3*(grams/amu))
      CG1 = 1.0D5*DSQRT(CG)   ! 10^5*G
      CG2 = CG1*0.432D0/CPI   ! 10^5*G*0.432/pi
      CGRT = 6.4D0*(CPI/4.32D4)**6*(1.0D11/CL)**5 ! 6.4*(pi/4.32e4)^6*(10^11/c)^5
      ! Read opacity, nuclear reaction and neutrino loss rate data
      fname = 'EZ_EOS.data'
      IF (.NOT. OpenToRead(fname,IO_EOS)) STOP 'FATAL ERROR: CANNOT OPEN EQUATION OF STATE DATA'
      READ (IO_EOS,*) KCSX, CZS, CH
      ! KCSX         number of subcompositions tabulated
      ! CZS          the metalicity
      ! CH           the initial hydrogen mass fraction
  990 FORMAT (1X, 10F7.3)
      READ (IO_EOS,990) CSX ! 10 numbers describing the subcompositions
      DO I = 1, KCSX ! CS(90,127,10) holds the opacity tables
         READ (IO_EOS,990) ((CS(JR, JT, I), JR = 1, 90), JT = 1, 127)
      END DO
      READ (IO_EOS,990) CNU, CRT ! neutrino loss rate and nuclear reaction rate tables
      CLOSE (IO_EOS)
      CALL init_neutrinorates ! setup tables for neutrino loss rate based on Itoh et.al., ApJS, 1996.
      ! Set up coefficients for cubic spline interpolation in opacity
      CALL OPSPLN
      ! read various nuclear data.
      fname = 'EZ_NUC.data'
      IF (.NOT. OpenToRead(fname,IO_NUC)) STOP 'FATAL ERROR: CANNOT OPEN EZ_NUC DATA'
  991 FORMAT (3(10F10.5,/),4F10.5,/,10F10.5,/,3F10.5,/,442(5f8.3,/),1P,12(10D12.4,/),32(9D12.4,/),4D12.4,/,9I12)
      READ (IO_NUC,991) xTGR,xGGR,(((xTAB(I,J,K),K=1,5),J=1,13),I=1,34),Qvals,QNT,CZA,CZB,CZC,CZD,VZ,CBN,CAN,COM,CHI,Cvals,KZN
      QPP=Qvals(1); Q33=Qvals(2); Q34=Qvals(3); QBE=Qvals(4); QBP=Qvals(5)
      QPC=Qvals(6); QPNA=Qvals(7); QPO=Qvals(8); Q3A=Qvals(9); QAC=Qvals(10)
      QAN=Qvals(11); QAO=Qvals(12); QANE=Qvals(13); QCCA=Qvals(14); QCO=Qvals(15)
      QOO=Qvals(16); QGNE=Qvals(17); QGMG=Qvals(18); QCCG=Qvals(19); QPNG=Qvals(20)
      CH2=Cvals(1); C1=Cvals(2); C2=Cvals(3); C3=Cvals(4)
      DO I = 1, NEL
         CAN1pt5(I) = CAN(I)*DSQRT(CAN(I))
         LN_CAN1pt5(I) = DLOG(CAN1pt5(I))
      END DO
      ! Read Bol. Corr, U-B, B-V table. Read nuclear reaction (QRT) and neutrino 
      ! (QNT) Q values, in MeV; constants for electron screening (CZA, CZB, CZC, 
      ! CZD, VZ); atomic parameters (CBN, KZN), with masses (CAN) consistent with 
      ! Q-values; ionization potentials (CHI) and statistical weights (COM); 
      ! molecular hydrogen parameters (CH2)
      CLOSE (IO_NUC)

      END SUBROUTINE SETUP


      LOGICAL FUNCTION Start_Read_Model ( initial_M )
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      CALL Set_Defaults
      AGE = 0D0; SM = initial_M
      Start_Read_Model = .TRUE.
      END FUNCTION Start_Read_Model
      
      LOGICAL FUNCTION Finish_Read_Model ( initial_M, Initialize_Params, adjust )
      USE star_data
      USE ez_report
      USE ez_cycle
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      INTERFACE
         SUBROUTINE Initialize_Params
         END SUBROUTINE Initialize_Params
      END INTERFACE
      LOGICAL, INTENT(IN) :: adjust
      INTEGER JO
      IF (initial_M .LE. 1D0) THEN
         DTY = 1D5
      ELSE ! use rough estimate of nuclear timescale to reduce initial time step for M > Msolar
         DTY = 1D5 / (initial_M**2.5D0)
      END IF
      DTY = timestep_decrement*DTY ! reduce a bit more as if had a forced backup.
      accuracy_target = 1D-6
      JO = JO_NULL; JNN = 0; JMOD = 0
      MC=0D0; CDD=0D0; nuc_Timescale=0D0
      CALL Initialize_Params
      DT = CSY*DTY ! DT is timestep in seconds
      PREV_DT = DT
      AGE = AGE - DTY
      CALL NEXTDT
      IF (adjust) CALL Adjust_Model ( CMSN*initial_M )
      ! store some numbers for possible restart with BACKUP
      JHOLD = 2
      CALL Save_State(.FALSE.)
      CALL Complete_Model ! calculate lots of stuff about the model just read
      CALL Check_Final ( JO )
      Finish_Read_Model = .TRUE.
      END FUNCTION Finish_Read_Model
      
      SUBROUTINE Convert_fort16
      DOUBLE PRECISION, PARAMETER :: min_Mass = 0.0999D0, max_Mass = 100D0
      INTEGER, PARAMETER :: linelen = 257
      DOUBLE PRECISION :: mass
      INTEGER :: K
      CHARACTER (LEN=strlen) :: fname, input_fname
      CHARACTER (LEN=linelen) :: line_buffer
      fname = 'new_EZ_ZAMS.data'
      input_fname = 'fort.16'
      EZ_DATA_DIR = '../metals/z001'
      IF (.NOT. OpenToRead(input_fname,IO_TMP)) RETURN
      IF (.NOT. OpenToWrite(fname,IO_ZAMS)) RETURN
      DO 
         READ (IO_TMP, *) mass
         WRITE(*,*) mass
         IF ( mass .GT. max_Mass ) EXIT
         DO K = 1, file_N_SHELLs
            READ (IO_TMP, '(A)') line_buffer
            IF ( mass .GE. min_Mass ) WRITE (IO_ZAMS, '(A)') line_buffer
         END DO
      END DO
      CLOSE (IO_ZAMS); CLOSE (IO_TMP)
      END SUBROUTINE Convert_fort16
      
      LOGICAL FUNCTION Load_ZAMS_Model ( initial_M, fname )
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      INTEGER IM
      DOUBLE PRECISION DMASS, MLO, MHI, ML
      INTEGER I, K
      CHARACTER (LEN=strlen) :: fname
      
      Load_ZAMS_Model = .FALSE.
      
      ! Added for DarkStars Z=0 protostar edition, Pat Scott 080813
      IF (metals_dir .EQ. '../metals/z0_proto') THEN
         
         fname = 'protostar.start'
         IF (.NOT. OpenToRead(fname,IO_ZAMS)) RETURN
         WRITE(*,*) 'Initial protostellar profile taken from protostar.start'
         DO K = 1, file_N_SHELLs
            READ (IO_ZAMS, *) H(1:NUMV,K)
         END DO
         DO K = 1, file_N_SHELLs
            DH(1:NUMV,K) = 0D0
            !READ (IO_ZAMS, *) DH(1:NUMV, K)
         END DO
      ELSE
      ! End DarkStars add

         MLO = -1D0; DMASS = .025D0; MHI = 2D0
         ML = log10(initial_M)
         IF ( ML .GT. MHI .OR. ML .LT. MLO ) THEN
            WRITE(*,'(a,3(f6.2,a))') 'Requested mass ', initial_M, ' out of allowed range of ', 10**MLO, ' to ', 10**MHI
            RETURN
         END IF
         IM = (ML - MLO)/DMASS + 1.501D0 ! pick which ZAMS model to use
         IF (.NOT. OpenToRead(fname,IO_ZAMS)) RETURN
         DO I = 1, IM-1
            DO K = 1, file_N_SHELLs
               READ (IO_ZAMS, *)
            END DO
         END DO
         DO K = 1, file_N_SHELLs ! this is the one we want
            DH(1:NUMV,K) = 0D0
            READ (IO_ZAMS, *) H(1:NUMV, K)
         END DO

      ENDIF

      CLOSE (IO_ZAMS)
      Load_ZAMS_Model = .TRUE.
      END FUNCTION Load_ZAMS_Model
      
      LOGICAL FUNCTION Read_ZAMS ( initial_M, Initialize_Params )
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      INTERFACE
         SUBROUTINE Initialize_Params
         END SUBROUTINE Initialize_Params
      END INTERFACE
      CHARACTER (LEN=strlen) :: fname
      fname = 'EZ_ZAMS.data'
      Read_ZAMS = .FALSE.
      IF (.NOT. Start_Read_Model ( initial_M )) RETURN
      IF (.NOT. Load_ZAMS_Model ( initial_M, fname )) RETURN
      IF (.NOT. Finish_Read_Model ( initial_M, Initialize_Params, .TRUE. )) RETURN
      Read_ZAMS = .TRUE.
      END FUNCTION Read_ZAMS
      
      SUBROUTINE Adjust_Model ( TM )  ! applies to a model just input from file, so uses file_N_SHELLs constant
      ! input TM is new total mass.  adjust model to this.  set values for added variables beyond stored ZAMS vars.
      ! if TM is zero, skip adjusting the mass.
      USE star_data
      USE ez_shell
      DOUBLE PRECISION, INTENT(IN) :: TM
      DOUBLE PRECISION :: VD, Q1, Q2, VC, HPC, VQK, QQE, QA(file_N_SHELLs)
      INTEGER IK
      QQE = 0D0
      VC = H(V_M, file_N_SHELLs)
      MC = TM
      VAR(1:NUMV) = H(1:NUMV, file_N_SHELLs)
      DVAR = 0D0
      CALL FUNCS1 ( file_N_SHELLs, -2 )
      HPC = DSQRT(P/(CG*RHO*RHO))
      MC = 3.5D-33*RHO*HPC**3 ! central core mass used in computing mass term for equations
      VD = TM/H(V_M, 1) ! TM is new total mass; VD is factor for rescaling mass
      IF (TM .GT. 0D0) H(V_M, 1:file_N_SHELLs) = VD*H(V_M, 1:file_N_SHELLs) ! rescale the mass to new TM
      H(V_Q_dk, 1:file_N_SHELLs) = 1D0 ! reset the gradient of the mesh spacing function
      DO IK = 1, file_N_SHELLs
         VAR(1:NUMV) = H(1:NUMV, IK)
         DVAR = 0D0
         CALL FUNCS1 ( IK, -2 )
         QA(IK) = MESH_Q ! value of mesh spacing function at this meshpoint
         ! Integration for mesh spacing.
         QQE = QQE + CT(2)*EX*MK*1D-8
         QA(IK) = QA(IK) + QQE
      END DO
      Q1 = QA(1)
      Q2 = (file_N_SHELLs - 1D0)/(QA(file_N_SHELLs) - Q1) ! dk/dQ
      VQK = 1D0/Q2 ! dQ/dk, the gradient of the mesh spacing function
      H(V_Q_dk, 1:file_N_SHELLs) = VQK

      END SUBROUTINE Adjust_Model
      
      INTEGER FUNCTION TradeXY_Check_Model() ! applies to the current model, so uses the variable N_SHELLs
      USE star_data
      USE star_controls
      DOUBLE PRECISION :: Y, DY, shellDY, shellDX
      INTEGER :: IK
      TradeXY_Check_Model = KEEP_GOING
      Y = star_Mass_He / star_Mass
      DY = TradeXY_target - Y
      IF ( (model_Number <= TradeXY_model .or. DABS(Y - TradeXY_previous) > 1d-6) .and. DABS(DY) .GE. 1D-3 ) THEN
         IF (DY .GT. tradeXY_maxD) DY = tradeXY_maxD
         IF (DY .LT. -tradeXY_maxD) DY = -tradeXY_maxD
         DO IK = 1, N_SHELLs
            shellDX = H(V_XH,IK)*0.5d0
            shellDY = H(V_XHE,IK)*DY
            IF (shellDX < shellDY) shellDY = shellDX
            H(V_XHE,IK) = H(V_XHE,IK) + shellDY
            H(V_XH,IK) = H(V_XH,IK) - shellDY
         END DO
      ELSE
         TradeXY_Check_Model = TERMINATE
      END IF
      TradeXY_previous = Y
      TradeXY_model = model_Number
      END FUNCTION TradeXY_Check_Model

      LOGICAL FUNCTION Trade_X_for_Y ( new_Y, Check_Model )
      USE star_data
      USE ez_report
      USE ez_cycle
      DOUBLE PRECISION, INTENT(IN) :: new_Y
      INTERFACE
         INTEGER FUNCTION Check_Model()
         END FUNCTION Check_Model
      END INTERFACE
      INTEGER :: JO
      JO = JO_AOK
      tradeXY_maxD=0.005D0
      TradeXY_previous = -1d0
      TradeXY_target = new_Y
      TradeXY_model = -1
      CALL Evolve1 ( JO, 99999, Check_Model, TradeXY_Check_Model )
      CALL Complete_Model
      CALL Check_Final ( JO )
      Trade_X_for_Y = (JO .EQ. JO_AOK)
      END FUNCTION Trade_X_for_Y
      
      FUNCTION OpenToWrite(fname,io_unit)
      LOGICAL :: OpenToWrite
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: io_unit
      OpenToWrite = OpenDataFile(fname,io_unit,.FALSE.)
      END FUNCTION OpenToWrite
      
      FUNCTION OpenToRead(fname,io_unit)
      LOGICAL :: OpenToRead
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: io_unit
      OpenToRead = OpenDataFile(fname,io_unit,.TRUE.)
      END FUNCTION OpenToRead
      
      SUBROUTINE GetDataFileName(fname,fullname)
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      CHARACTER (LEN=strlen), INTENT(OUT) :: fullname
      IF (LEN(TRIM(EZ_DATA_DIR)) .GT. 0) THEN
         fullname = TRIM(EZ_DATA_DIR) // '/' // TRIM(fname)
      ELSE
         fullname = trim(fname)
      END IF
      END SUBROUTINE GetDataFileName
      
      FUNCTION OpenDataFile(fname,io_unit,read)
      LOGICAL :: OpenDataFile
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: io_unit
      LOGICAL, INTENT(IN) :: read
      INTEGER :: ios
      CHARACTER (LEN=strlen) :: fullname
      CALL GetDataFileName(fname,fullname)
      ios = 0
      IF (read) THEN
         OPEN(UNIT=io_unit, FILE=TRIM(fullname), ACTION='READ', STATUS='OLD', IOSTAT=ios)
      ELSE
         OPEN(UNIT=io_unit, FILE=TRIM(fullname), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      END IF
      IF (ios /= 0) THEN
         WRITE(*,*) 'Failed to open file ', TRIM(fullname), read, ios
         OpenDataFile = .FALSE.
      ELSE
         OpenDataFile = .TRUE.
      END IF
      END FUNCTION OpenDataFile

      SUBROUTINE Test_logF_EoS(logT, logF, out_filename)
      USE ez_opacity
      USE ez_shell
      USE ez_state
      USE ez_state_data
      DOUBLE PRECISION, INTENT(IN) :: logT, logF
      CHARACTER (LEN=strlen), INTENT(IN) :: out_filename
      INTEGER :: iounit, ios
      DOUBLE PRECISION :: XA(NEL), kappa, temp
      CALL FUNCS1(N_SURF_shell, 0)
      XA(N_H) = XH; XA(N_HE) = XHE; XA(N_C) = XC; XA(N_N) = XN; XA(N_O) = XO; 
      XA(N_NE) = XNE; XA(N_MG) = XMG; XA(N_SI) = XSI; XA(N_FE) = XFE
      temp = 10D0**logT
      ios = 0; iounit = 44
      OPEN(UNIT=iounit, FILE=trim(out_filename), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(*,*) 'Write_logF_EoS: Failed to open file ', TRIM(out_filename)
         RETURN
      END IF
      CALL STATEF(logF*CLN, logT*CLN, temp, XA, H_HE_C_N_O_ions)
      kappa = Opacity(logT, LNRHO/CLN, XH, XHE)
      WRITE(iounit,'(999(A,1x,E18.12,1x))') "PSI", PSI, "logT", logT, "logF", logF, "RHO", RHO, "P", P, log10(P)
      CLOSE(iounit)
      END SUBROUTINE Test_logF_EoS

      SUBROUTINE Write_logF_EoS(lgTmn, lgTmx, lgTsteps, logFmin, logFmax, logFsteps, out_filename)
      USE ez_opacity
      USE ez_shell
      USE ez_state
      USE ez_state_data
      DOUBLE PRECISION, INTENT(IN) :: lgTmn, lgTmx, logFmin, logFmax
      INTEGER, INTENT(IN) :: lgTsteps, logFsteps
      CHARACTER (LEN=strlen), INTENT(IN) :: out_filename
      INTEGER :: iounit, i, j, ios
      DOUBLE PRECISION :: XA(NEL), logTdelta, logFdelta, logT, logF, kappa, temp
      ios = 0; iounit = 44
      OPEN(UNIT=iounit, FILE=trim(out_filename), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(*,*) 'Write_logF_EoS: Failed to open file ', TRIM(out_filename)
         RETURN
      END IF
      CALL FUNCS1(N_SURF_shell, 0)
      XA(N_H) = XH; XA(N_HE) = XHE; XA(N_C) = XC; XA(N_N) = XN; XA(N_O) = XO; 
      XA(N_NE) = XNE; XA(N_MG) = XMG; XA(N_SI) = XSI; XA(N_FE) = XFE
      logTdelta = (lgTmx-lgTmn)/(lgTsteps-1);
      logFdelta = (logFmax-logFmin)/(logFsteps-1);
      WRITE(iounit,'(999(E18.12,1x))') init_Z, lgTmn, lgTmx, lgTsteps, logFmin, logFmax, logFsteps
      DO j=1,lgTsteps
         logT = lgTmn+(j-1)*logTdelta
         temp = 10D0**logT
         DO i=1,logFsteps
            logF = logFmin+(i-1)*logFdelta
            CALL STATEF(logF*CLN, logT*CLN, temp, XA, H_HE_C_N_O_ions)
            kappa = Opacity(logT, LNRHO/CLN, XH, XHE)
            WRITE(iounit,'(999(E18.12,1x))') i, j, PSI, logT, logF
            WRITE(iounit,'(999(E18.12,1x))') RHO, P, kappa, NE, NE1, U, S, SCP, PCORR, log10(P)
         END DO
      END DO
      CLOSE(iounit)
      END SUBROUTINE Write_logF_EoS
      
      LOGICAL FUNCTION Get_logRHO_EoS_Info(lnF_cur, logT, temp, XA, ions, logRHO, tol, min_lnF)
      USE ez_state
      DOUBLE PRECISION, INTENT(INOUT) :: lnF_cur
      INTEGER, INTENT(IN) :: ions
      DOUBLE PRECISION, INTENT(IN) :: logT, temp, logRHO, tol, min_lnF, XA(NEL)
      DOUBLE PRECISION :: lnF_step, lnT, lnRHO_target, lnF_prev, lnF_mid
      INTEGER :: i
      INTEGER, PARAMETER :: MAXITER = 10
      Get_logRHO_EoS_Info = .FALSE.
      lnF_step = 0.5D0; lnT = logT*CLN; lnRHO_target = logRHO*CLN
      DO
         IF ( LNRHO .GE. lnRHO_target ) EXIT
         IF ( GAM .GT. 160D0 ) RETURN
         lnF_cur = lnF_cur + lnF_step
         CALL STATEF(lnF_cur, lnT, temp, XA, ions)
      END DO
      DO
         lnF_prev = lnF_cur
         lnF_cur = lnF_cur - lnF_step
         IF (lnF_cur .LT. min_lnF) RETURN
         CALL STATEF(lnF_cur, lnT, temp, XA, ions)
         IF ( LNRHO .LE. lnRHO_target ) EXIT
      END DO
      Get_logRHO_EoS_Info = .TRUE.
      lnF_mid = lnF_cur
      DO i=0,MAXITER
         IF (abs(LNRHO-lnRHO_target) .LE. tol) THEN
            lnF_cur = lnF_mid
            RETURN
         END IF
         lnF_mid = 0.5D0*(lnF_cur+lnF_prev)
         CALL STATEF(lnF_mid, lnT, temp, XA, ions)
         IF ( LNRHO .GE. lnRHO_target ) THEN
            lnF_prev = lnF_mid
         ELSE
            lnF_cur = lnF_mid
         END IF
      END DO
      WRITE(*,*) 'Get_logRHO_EoS_Info failed to converge in allowed number of steps'
      END FUNCTION Get_logRHO_EoS_Info
      
      LOGICAL FUNCTION Open_for_EoS(out_dir, outfile, iounit)
      CHARACTER (LEN=strlen), INTENT(IN) :: out_dir, outfile
      INTEGER, INTENT(IN) :: iounit
      CHARACTER (LEN=strlen) :: fullname
      INTEGER :: ios
      ios = 0
      fullname = TRIM(out_dir) // '/' // TRIM(outfile)
      OPEN(UNIT=iounit, FILE=fullname, ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(*,*) 'Write_logRHO_EoS: Failed to open file ', TRIM(outfile)
         Open_for_EoS = .FALSE.
         RETURN
      END IF
      Open_for_EoS = .TRUE.
      END FUNCTION Open_for_EoS
      
      SUBROUTINE Write_logRHO_EoS(lgTmn, lgTmx, lgTsteps, lgRHOmn, lgRHOmx, lgRHOsteps, out_dir)
      USE ez_opacity
      USE ez_shell
      USE ez_state
      USE ez_nuclear
      DOUBLE PRECISION :: lgTmn, lgTmx, lgRHOmn, lgRHOmx
      INTEGER, INTENT(IN) :: lgTsteps, lgRHOsteps
      CHARACTER (LEN=strlen), INTENT(IN) :: out_dir
      INTEGER :: iounit_P, iounit_K, iounit_I, iounit_T, iounit_RHO, iounit_S, iounit_U, iounit_PSI, iounit_GRAD_AD
      INTEGER :: i, j, k, tempSteps, rhoSteps, iounit_EX, iounit_ENX
      DOUBLE PRECISION :: NE_TOT_val, NE_free, NE_bound, ratio, F, WF
      DOUBLE PRECISION :: XA(NEL), logTdelta, logRHOdelta, logT, logRHO, lnF_cur, kappa, temp
      DOUBLE PRECISION :: Ps(999), Ks(999), Ss(999), Is(999), Us(999), PSIs(999), GRAD_ADs(999), EXs(999), ENXs(999)
      DOUBLE PRECISION, PARAMETER :: tol = 1D-3
      DOUBLE PRECISION, PARAMETER :: LNF_PSI_MAX = 9D0 ! this value of lnF gives PSI=180
      DOUBLE PRECISION, PARAMETER :: min_lnF = -1D2 ! this value of lnF gives PSI=-99.4
      LOGICAL valid
      CHARACTER (LEN=strlen) :: fname
      rhoSteps = lgRHOsteps
      IF ( rhoSteps .GT. 999 ) THEN
         rhoSteps = 999
      END IF
      tempSteps = lgTsteps
      IF ( tempSteps .GT. 999 ) THEN
         tempSteps = 999
      END IF
      iounit_P = 44; iounit_K = 45; iounit_I = 46; iounit_T = 47
      iounit_RHO = 48; iounit_S = 49; iounit_U = 50; iounit_PSI = 51; iounit_GRAD_AD = 52;
      iounit_EX = 53; iounit_ENX = 54
      fname = 'Pressure_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_P)) RETURN
      fname = 'Opacity_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_K)) RETURN
      fname = 'Ionization_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_I)) RETURN
      fname = 'logTs_for_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_T)) RETURN
      fname = 'logRHOs_for_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_RHO)) RETURN
      fname = 'Entropy_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_S)) RETURN
      fname = 'Energy_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_U)) RETURN
      fname = 'PSI_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_PSI)) RETURN
      fname = 'GRAD_AD_EoS.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_GRAD_AD)) RETURN
      fname = 'eps_nuc.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_EX)) RETURN
      fname = 'sneut.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_ENX)) RETURN
      CALL FUNCS1(N_SURF_shell, 0)
      XA(N_H) = XH; XA(N_HE) = XHE; XA(N_C) = XC; XA(N_N) = XN; XA(N_O) = XO; 
      XA(N_NE) = XNE; XA(N_MG) = XMG; XA(N_SI) = XSI; XA(N_FE) = XFE
      logTdelta = (lgTmx-lgTmn)/(tempSteps-1);
      logRHOdelta = (lgRHOmx-lgRHOmn)/(rhoSteps-1);
      k = 0

      XA(N_H) = 0.600956E+00
      XA(N_HE) = 0.383529E+00 + 0.300079D-03
      XA(N_C) = 0.220227D-04
      XA(N_N) = 0.706081D-02
      XA(N_O) = 0.670419D-02


      XA(N_H) = 0.585383E+00
      XA(N_HE) = 0.396124E+00 + 0.150393D-04
      XA(N_C) = 0.289058D-04
      XA(N_N) = 0.665889D-02
      XA(N_O) = 0.668175D-02

      XA(N_MG) = CMG*CZS
      XA(N_SI) = CSI*CZS
      XA(N_FE) = CFE*CZS
      XA(N_NE) = 1d0 - (XA(N_H) + XA(N_HE) + XA(N_C) + XA(N_N) + XA(N_O) + XA(N_MG) + XA(N_SI) + XA(N_FE))
      IF (DABS(XA(N_NE) - 0.162163D-02) > 1d-2) THEN
         WRITE(*,*) 'XA(N_NE)', XA(N_NE)
         WRITE(*,*) 'DABS(XA(N_NE) - 0.162163D-02)', DABS(XA(N_NE) - 0.162163D-02)
         STOP 'bad composition for Write_logRHO_EoS'
      END IF

      DO j=1,tempSteps
         logT = lgTmx-(j-1)*logTdelta
         WRITE(iounit_T,*) logT
         temp = 10D0**logT
         lnF_cur = LNF_PSI_MAX
         DO i=1,rhoSteps
            logRHO = lgRHOmn+(i-1)*logRHOdelta
            IF (j .EQ. 1) WRITE(iounit_RHO,*) logRHO
            CALL STATEF(lnF_cur, logT*CLN, temp, XA, H_HE_C_N_O_ions)
            valid = Get_logRHO_EoS_Info(lnF_cur, logT, temp, XA, H_HE_C_N_O_ions, logRHO, tol, min_lnF)
            IF (valid) THEN
               EX = NUCRAT ( logT*CLN, ZT, NA, NI, NE, AVM, DEXP(LNRHO), ENX )
               if (EN == 0) EN = 1D-99
               if (ENX == 0) ENX = -1D-99
               EXs(i) = LOG10(EX); ENXs(i) = LOG10(-ENX)
               logRHO = LNRHO/CLN
               kappa = Opacity(logT, logRHO, XH, XHE)
               Ps(i) = LNP/CLN
               Ks(i) = LOG10(kappa)
               Ss(i) = LOG10(S)
               Us(i) = LOG10(U)
               F = exp(lnF_cur)
               WF = sqrt(1D0 + F)
               PSIs(i) = lnF_cur + 2D0*(WF - DLOG(1D0 + WF))
               NE_TOT_val = NE
               NE_free = NE1
               NE_bound = NE_TOT_val - NE_free
               IF ( NE_free * 1D-15 .GT. NE_bound ) THEN
                  ratio = 1D-15
               ELSE
                  ratio = NE_bound/NE_free
               END IF
               Is(i) = LOG10(ratio)
               GRAD_ADs(i) = GRADA
            ELSE
               Ps(i) = -1D5
               Ks(i) = -1D5
               Is(i) = -1D5
               Ss(i) = -1D5
               Us(i) = -1D5
               PSIs(i) = -1D5
               GRAD_ADs(i) = -1D5
               EXs(i) = -1D5
               ENXs(i) = -1D5
            END IF
            k = k + 1
         END DO
         WRITE(iounit_P,*) Ps(1:rhoSteps)
         WRITE(iounit_K,*) Ks(1:rhoSteps)
         WRITE(iounit_I,*) Is(1:rhoSteps)
         WRITE(iounit_S,*) Ss(1:rhoSteps)
         WRITE(iounit_U,*) Us(1:rhoSteps)
         WRITE(iounit_PSI,*) PSIs(1:rhoSteps)
         WRITE(iounit_GRAD_AD,*) GRAD_ADs(1:rhoSteps)
         WRITE(iounit_EX,*) EXs(1:rhoSteps)
         WRITE(iounit_ENX,*) ENXs(1:rhoSteps)
      END DO
      WRITE(*,*) 'Output', k, 'EoS entries'
      CLOSE(iounit_P)
      CLOSE(iounit_K)
      CLOSE(iounit_I)
      CLOSE(iounit_S)
      CLOSE(iounit_U)
      CLOSE(iounit_PSI)
      CLOSE(iounit_RHO)
      CLOSE(iounit_T)
      CLOSE(iounit_GRAD_AD)
      CLOSE(iounit_EX)
      CLOSE(iounit_ENX)
      END SUBROUTINE Write_logRHO_EoS
      
      
      SUBROUTINE Write_EPSNUC(lgTmn, lgTmx, lgTsteps, lgRHOmn, lgRHOmx, lgRHOsteps, out_dir, XA)
      USE ez_opacity
      USE ez_shell
      USE ez_state
      USE ez_nuclear
      DOUBLE PRECISION :: XA(NEL), lgTmn, lgTmx, lgRHOmn, lgRHOmx
      INTEGER, INTENT(IN) :: lgTsteps, lgRHOsteps
      CHARACTER (LEN=strlen), INTENT(IN) :: out_dir
      INTEGER :: iounit_P, iounit_K, iounit_I, iounit_T, iounit_RHO, iounit_S, iounit_U, iounit_PSI, iounit_GRAD_AD
      INTEGER :: i, j, k, tempSteps, rhoSteps, iounit_EX, iounit_ENX
      DOUBLE PRECISION :: NE_TOT_val, NE_free, NE_bound, ratio, F, WF
      DOUBLE PRECISION :: logTdelta, logRHOdelta, logT, logRHO, lnF_cur, kappa, temp
      DOUBLE PRECISION :: Ps(999), Ks(999), Ss(999), Is(999), Us(999), PSIs(999), GRAD_ADs(999), EXs(999), ENXs(999)
      DOUBLE PRECISION, PARAMETER :: tol = 1D-3
      DOUBLE PRECISION, PARAMETER :: LNF_PSI_MAX = 9D0 ! this value of lnF gives PSI=180
      DOUBLE PRECISION, PARAMETER :: min_lnF = -1D2 ! this value of lnF gives PSI=-99.4
      LOGICAL valid
      CHARACTER (LEN=strlen) :: fname
      rhoSteps = lgRHOsteps
      IF ( rhoSteps .GT. 999 ) THEN
         rhoSteps = 999
      END IF
      tempSteps = lgTsteps
      IF ( tempSteps .GT. 999 ) THEN
         tempSteps = 999
      END IF
      iounit_T = 47; iounit_RHO = 48
      iounit_EX = 49; iounit_ENX = 50
      fname = 'logTs.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_T)) RETURN
      fname = 'logRHOs.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_RHO)) RETURN
      fname = 'eps_nuc.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_EX)) RETURN
      fname = 'sneut.data'
      IF (.NOT. Open_for_EoS(out_dir, fname, iounit_ENX)) RETURN
      CALL FUNCS1(N_SURF_shell, 0)

      logTdelta = (lgTmx-lgTmn)/(tempSteps-1);
      logRHOdelta = (lgRHOmx-lgRHOmn)/(rhoSteps-1);

      DO j=1,tempSteps
         logT = lgTmx-(j-1)*logTdelta
         WRITE(iounit_T,*) logT
         temp = 10D0**logT
         lnF_cur = LNF_PSI_MAX
         DO i=1,rhoSteps
            logRHO = lgRHOmn+(i-1)*logRHOdelta
            IF (j .EQ. 1) WRITE(iounit_RHO,*) logRHO
            CALL STATEF(lnF_cur, logT*CLN, temp, XA, H_HE_C_N_O_ions)
            valid = Get_logRHO_EoS_Info(lnF_cur, logT, temp, XA, H_HE_C_N_O_ions, logRHO, tol, min_lnF)
            IF (valid) THEN
               EX = NUCRAT ( logT*CLN, ZT, NA, NI, NE, AVM, DEXP(LNRHO), ENX )
               if (EX <= 1D-99 .or. EX-EX /= 0d0) EX = 1D-99
               if (ENX >= -1D-99 .or. ENX-ENX /= 0d0) ENX = -1D-99
               EXs(i) = LOG10(EX); ENXs(i) = LOG10(-ENX)
            ELSE
               EXs(i) = -99d0
               ENXs(i) = -99d0
            END IF
            k = k + 1
         END DO
         WRITE(iounit_EX,*) EXs(1:rhoSteps)
         WRITE(iounit_ENX,*) ENXs(1:rhoSteps)
      END DO
      WRITE(*,*) 'Output', k, 'entries'
      CLOSE(iounit_RHO)
      CLOSE(iounit_T)
      CLOSE(iounit_EX)
      CLOSE(iounit_ENX)
      END SUBROUTINE Write_EPSNUC
      
      LOGICAL FUNCTION Read_ZA_He_MS ( initial_M, Initialize_Params )
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      INTERFACE
         SUBROUTINE Initialize_Params
         END SUBROUTINE Initialize_Params
      END INTERFACE
      Read_ZA_He_MS = .FALSE.
      IF (.NOT. Start_Read_Model ( initial_M )) RETURN
      IF (.NOT. Load_ZA_He_MS_Model ( initial_M )) RETURN
      IF (.NOT. Finish_Read_Model ( initial_M, Initialize_Params, .FALSE. )) RETURN
      Read_ZA_He_MS = .TRUE.
      END FUNCTION Read_ZA_He_MS
      
      LOGICAL FUNCTION Load_ZA_He_MS_Model ( initial_M )
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      INTEGER K
      CHARACTER (LEN=strlen) :: fname
      Load_ZA_He_MS_Model = .FALSE.
      fname = 'ZA_He_MS.data'
      IF (.NOT. OpenToRead(fname,IO_ZAMS)) RETURN
      DO K = 1, N_SHELLs
         DH(1:NUMV,K) = 0D0
         READ (IO_ZAMS, *) H(1:NUMV, K)
      END DO
      CLOSE (IO_ZAMS)
      Load_ZA_He_MS_Model = .TRUE.
      END FUNCTION Load_ZA_He_MS_Model
      
       
      LOGICAL FUNCTION Read_ZAHB ( initial_M, Initialize_Params )
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      INTERFACE
         SUBROUTINE Initialize_Params
         END SUBROUTINE Initialize_Params
      END INTERFACE
      Read_ZAHB = .FALSE.
      IF (.NOT. Start_Read_Model ( initial_M )) RETURN
      IF (.NOT. Load_ZAHB_Model ( initial_M )) RETURN
      IF (.NOT. Finish_Read_Model ( initial_M, Initialize_Params, .FALSE. )) RETURN
      Read_ZAHB = .TRUE.
      END FUNCTION Read_ZAHB
     
      LOGICAL FUNCTION Load_ZAHB_Model ( initial_M )
      DOUBLE PRECISION, INTENT(IN) :: initial_M
      INTEGER K
      CHARACTER (LEN=strlen) :: fname
      Load_ZAHB_Model = .FALSE.
      fname = '../zahb.data'
      IF (.NOT. OpenToRead(fname,IO_ZAMS)) RETURN
      DO K = 1, N_SHELLs
         DH(1:NUMV,K) = 0D0
         READ (IO_ZAMS, *) H(1:NUMV, K)
      END DO
      CLOSE (IO_ZAMS)
      Load_ZAHB_Model = .TRUE.
      END FUNCTION Load_ZAHB_Model
      
      END MODULE ez_setup
