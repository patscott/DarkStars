      MODULE ez_shell
      USE ez_data
      USE ez_nuclear_data
      USE ez_state_data
      USE ez_shell_data
      USE ez_cycle_data
      USE star_constants
      USE star_controls
      USE star_data
      USE star_extras
      IMPLICIT NONE 
      
      CONTAINS

      
      SUBROUTINE FUNCS1 ( K, JI )
      ! evaluate the functions of the independent variables at meshpoint K
      ! input K   is the shell number (=1 for surface; =N_CNTR_shell for center)
      ! input JI  passed to STATEL to control the use of cached values
      ! on entry, refers to DT, VAR(*), and DVAR(*) to get args for functions
      ! on exit, has set FN1(*) with function values
      ! if JI==0, then on exit has also set SV(*) to state variable values
      ! NO OTHER EXTERNAL DEPENDENCIES on things that change during a run
      ! so if called with same K, JI args and same DT, VAR(*), and DVAR(*) globals,
      ! then can use cached results (if have them) without recalculating. 
      USE ez_convection
      USE star_controls
      INTEGER, INTENT(IN) :: K, JI 
      DOUBLE PRECISION :: APM, WC1, GRAV
      DOUBLE PRECISION :: R2, R2M, VM3
      DOUBLE PRECISION ::ATM, VRR, VPP, VMF, VMM, VTT, XMC, overshoot_param, FAC
      DOUBLE PRECISION :: CBRT, PS, VX
      CBRT(VX) = DEXP(DLOG(VX)*C3RD) ! CBRT(VX) is cube root of VX
      VX = 0 ! for the Absoft compiler.
      CALL Unload_VARs(K)
      T = exp(LNT)
      BC_Psurf = 0D0; BC_Tsurf = 0D0; BC_Msurf = 0D0
      VM3 = CBRT(M) ! M is the mass variable; VM3 is cube root of the mass
      M_dt = DM/DT ! M_dt is dM/dt at current meshpoint
      CALL Set_Composition_Vars
      ! Equation of state
      CALL STATEL ( JI )
      R = exp(LNR); R2 = R*R; R2M = 2D0/(C4PI*RHO*R) ! R2M is d_R2 / d_M
      CALL Get_GRAD(R2, GRAV, GAMMA1, APM, WC1)
      PQ = CT(4)*LNP + CT(5)*DLOG(P + CT(9)) ! PQ is pressure term for Q
      VPP = CT(4) + CT(5)*P/(P + CT(9)) ! VPP = d_PQ / d_LNP, with LNP = ln(P)
      ATM = GRAD*APM ! ATM is d_ln(T) / d_M
      TQ = CT(7)*DLOG(T/(T + CT(10))) ! TQ is temperature term for Q and for the T equation
      VTT = CT(7)*CT(10)/(T + CT(10)) ! VTT is d_TQ / d_ln(T)
      XMC = CT(6)*CBRT(MC*MC) ! proportional to center mass^(2/3)
      VMF = XMC + VM3*VM3 ! VM3*VM3 is M^(2/3)
      MQ = DLOG(XMC/VMF) ! MQ is mass term for Q
      VMM = -2D0*C3RD/(VMF*VM3) ! VMM is d_MQ / d_M
      RQ = -CT(3)*DLOG(R2/CT(8) + 1D0) ! RQ is R term for Q; as R -> 0, RQ -> 0 like R^2
      VRR = -CT(3)/(R2 + CT(8)) ! VRR is d_RQ / d_R2
      MESH_Q = PQ + TQ + MQ + RQ
      QM = VPP*APM + VTT*ATM + VMM + VRR*R2M  ! QM is d_Q / d_M
      QM = QM + CT(2)*EX*1D-8 ! EX is nuclear reaction energy release in ergs/g/sec
        ! [NB: default for CT(2) is 0D0, so this term doesn't actually contribute]
      MK = Q_dk/QM ! Q_dk is d_Q / d_k and QM is d_Q / d_M, so MK is d_M / d_k
      MQ_dk = VMM*MK     ! MQ_dk = (d_MQ / d_k) =  = (d_MQ / d_M)*(d_M / d_k)
      PQ_dk = VPP*APM*MK ! PQ_dk = (d_PQ / d_k) = (d_PQ / d_ln(P))*(d_ln(P) / d_M)*(d_M / d_k)
      TQ_dk = VTT*ATM*MK ! TQ_dk = (d_TQ / d_k) = (d_TQ / d_ln(T))*(d_ln(T) / d_M)*(d_M / d_k)
      RQ_dk = VRR*R2M*MK ! RQ_dk = (d_RQ / d_k) = (d_RQ / d_R2)*(d_R2 / d_M)*(d_M / d_k)
      WT = R2*XHI*DT/(RQ_dk*RQ_dk)
      CALL Composition_Changes
      ! convective mixing
      IF ( XH .GE. 1.0D-7 ) THEN ! still have a little H1 left in this shell
         overshoot_param = overshoot_param_H ! overshoot_param_H is convective overshoot parameter for H-burning cores
      ELSE ! hydrogen exhausted
         overshoot_param = overshoot_param_He ! overshoot_param_He is convective overshoot parameter for He-burning cores
      END IF
      FAC = 4D-22*convective_diffusion_prefactor*R2/(R2M*R2M*MK)*WC1**C3RD*XHI

      IF ( K .LT. number_of_surface_shells_to_mix .OR. K .GT. N_CNTR_shell - number_of_center_shells_to_mix ) THEN
         SG = FAC ! always do convective mixing for extremes
      ELSE
         SG = Mixing(PR,PG,APM,M,overshoot_param,GRADR,GRADA,FAC)
         IF ( JMOD .LE. 10 ) SG = SG * 1D-4
      END IF

      IF ( .false. .and. mod(K,30) == 0 .and. JI == 0 ) then
         write(*,*) 
         write(*,*) 'K = ', K
         write(*,*) 'MC23 = ', CBRT(MC*MC)*1d22
         write(*,*) 'P = ', P
         write(*,*) 'T = ', T
         write(*,*) 'R = ', R*1d11
         write(*,*) 'M = ', M*1d33
         write(*,*) 'Q = ', MESH_Q
         write(*,*) 'PQ = ', PQ
         write(*,*) 'TQ = ', TQ
         write(*,*) 'RQ = ', RQ
         write(*,*) 'MQ = ', MQ
         write(*,*) 'dPQ_dP = ', VPP/P
         write(*,*) 'dTQ_dT = ', VTT/T
         write(*,*) 'dRQ_dR = ', VRR*2*R*1d-11
         write(*,*) 'dMQ_dM = ', VMM*1d-33
         write(*,*) 
      end if

      !Altered for DarkStars; Pat Scott July 14 2007 (pat@physto.se)
      CALL DarkStars_Luminosity_Eqn(DSQRT(VMF)/VM3,APM,VMF,XMC)
      !CALL Luminosity_Eqn(DSQRT(VMF)/VM3,APM,VMF,XMC)
      IF ( K .EQ. 1 ) CALL Surface_BCs(GRAV,R2) ! do surface boundary conditions
      CALL Save_Funcs(K,JI)
      END SUBROUTINE FUNCS1
      
      SUBROUTINE Get_GRAD(R2, GRAV, G1, APM, WC1)
      USE ez_convection
      DOUBLE PRECISION, INTENT(IN) :: R2, G1 ! G1 is gamma1, for calculating the sound speed
      DOUBLE PRECISION, INTENT(OUT) :: GRAV, APM, WC1
      DOUBLE PRECISION :: HPS, HPC, HP, CSOUND
      GRAV = 1.0D11*CG*M/R2 ! GRAV is local gravitational acceleration
      APM = -1.0D11*GRAV/(C4PI*R2*P) ! APM is d_ln(P) / d_M = -(G M)/(4 Pi r^4 P)
      GRADR = 1.0D11*FK*P*L/(4D0*C4PI*CL*GRAV*R2*PR) ! radiative temp gradient
      HPS = P/(GRAV*RHO) ! HPS is pressure scale height
      CSOUND = DSQRT(G1 * P/RHO) ! sound speed
      HPC = DSQRT(P/(CG*RHO*RHO)) ! at center, grav -> 0 so HPS -> infinity, so use HPC instead
      HP = DMIN1(HPS, HPC) ! HP pressure scale height for mixing length calculation
      GRAD = Grad_Star(GRADR, GRADA, SCP, HP, T, XHI, mixing_length_alpha, CSOUND, WCV, WL, WC1)
      END SUBROUTINE Get_GRAD
      
      SUBROUTINE Luminosity_Eqn(VML,APM,VMF,XMC)
      DOUBLE PRECISION, INTENT(IN) :: VML,APM,VMF,XMC
      ! VML = sqrt(mc^(2/3)+M^(2/3))/M^(1/3)
      DOUBLE PRECISION :: WTH
      LQ = L*VML ! the luminosity equation is in terms of LQ rather than L
      WTH = 2D0*T_DS_Dt/(T_DS_Dt**4 + 1D0) ! T_DS_Dt is prefactor for certain terms in luminosity eqn
         ! T_DS_Dt usually = 1, but can set T_DS_Dt = 0 ignore the terms
         ! for T_DS_Dt = 1, WTH = 1; for T_DS_Dt = 0, WTH = 0.
      IF ( JMOD .EQ. 0 ) WTH = 0D0
      ENG = -WTH*T*(SF*DLNF + ST*DLNT)/DT ! -T*sdot, energy from gravitational contraction in ergs/g/sec
      LQ1_dk = VML*MK*(EX + EN + ENX + extra_Energy_param + ENG - C3RD*XMC/VMF*L/M)
      LQ2_dk = VML*WTH*SCP*T*APM*MK*(GRAD-GRADA)
      END SUBROUTINE Luminosity_Eqn
      
      SUBROUTINE DarkStars_Luminosity_Eqn(VML,APM,VMF,XMC)
      ! Special version of Luminosity_Eqn created for DarkStars
      ! Pat Scott July 18 2007 (pat@physto.se)
      USE DkStrs_annihilation
      USE DkStrs_transport
      DOUBLE PRECISION, INTENT(IN) :: VML,APM,VMF,XMC
      ! VML = sqrt(mc^(2/3)+M^(2/3))/M^(1/3)
      DOUBLE PRECISION :: WTH
      DOUBLE PRECISION :: E_injected
      LQ = L*VML ! the luminosity equation is in terms of LQ rather than L
      WTH = 2D0*T_DS_Dt/(T_DS_Dt**4 + 1D0) ! T_DS_Dt is prefactor for certain terms in luminosity eqn
         ! T_DS_Dt usually = 1, but can set T_DS_Dt = 0 ignore the terms
         ! for T_DS_Dt = 1, WTH = 1; for T_DS_Dt = 0, WTH = 0.
      IF ( JMOD .EQ. 0 ) WTH = 0D0
      !energy injected by WIMP annihilation 
      !plus that carried by WIMPs through heat conduction [ergs/g/sec]
      E_injected = E_annihilation(R*1.d11, RHO) + E_transport(R*1.d11, RHO)
      ENG = -WTH*T*(SF*DLNF + ST*DLNT)/DT ! -T*sdot, energy from gravitational contraction in ergs/g/sec
      LQ1_dk = VML*MK*(EX + EN + ENX + extra_Energy_param + ENG + E_injected - C3RD*XMC/VMF*L/M)
      LQ2_dk = VML*WTH*SCP*T*APM*MK*(GRAD-GRADA)
      END SUBROUTINE DarkStars_Luminosity_Eqn
      
      SUBROUTINE Unload_VARs(K)
      INTEGER, INTENT(IN) :: K
      INTEGER :: J
      DOUBLE PRECISION :: V
      Q_dk=VAR(V_Q_dk)
      LNR=VAR(V_LNR)
      L=VAR(V_L)
      M=VAR(V_M); DM=DVAR(V_M)
      LNF=VAR(V_LNF); DLNF=DVAR(V_LNF)
      LNT=VAR(V_LNT); DLNT=DVAR(V_LNT)
      XH=VAR(V_XH); DXH=DVAR(V_XH)
      XHE=VAR(V_XHE); DXHE=DVAR(V_XHE)
      XC=VAR(V_XC); DXC=DVAR(V_XC)
      XO=VAR(V_XO); DXO=DVAR(V_XO)
      XNE=VAR(V_XNE); DXNE=DVAR(V_XNE)
      END SUBROUTINE Unload_VARs
      
      SUBROUTINE Save_Funcs (K,JI)
      INTEGER, INTENT(IN) :: K,JI
      INTEGER :: J
      DOUBLE PRECISION :: V
      FN1(F_BC_Psurf)=BC_Psurf; FN1(F_BC_Tsurf)=BC_Tsurf; FN1(F_M_dt)=M_dt; FN1(F_Q_dk)=Q_dk
      FN1(F_PQ)=PQ; FN1(F_PQ_dk)=PQ_dk; FN1(F_RQ)=RQ; FN1(F_RQ_dk)=RQ_dk
      FN1(F_TQ)=TQ; FN1(F_TQ_dk)=TQ_dk; FN1(F_LQ)=LQ; FN1(F_LQ1_dk)=LQ1_dk; FN1(F_LQ2_dk)=LQ2_dk
      FN1(F_MQ)=MQ; FN1(F_MQ_dk)=MQ_dk; FN1(F_SG)=SG; FN1(F_WT)=WT; FN1(F_XH)=XH; FN1(F_XH_dt)=XHT
      FN1(F_XO)=XO; FN1(F_XO_dt)=XOT; FN1(F_XHE)=XHE; FN1(F_XHE_dt)=XHET; FN1(F_XC)=XC; FN1(F_XC_dt)=XCT
      FN1(F_XNE)=XNE; FN1(F_XNE_dt)=XNET; FN1(F_BC_Msurf)=BC_Msurf
      IF ( JI .EQ. 0 ) THEN ! save the results
         CALL Save_SV_VARs
         FN(FN_VALUE, K, 1:NUMF) = FN1(1:NUMF)
      ELSEIF ( JI .GE. 1 .AND. JI .LE. NUMV ) THEN
         FN(JI, K, 1:NUMF) = FN1(1:NUMF)
      END IF
      END SUBROUTINE Save_Funcs
      
      SUBROUTINE Surface_BCs(GRAV, R2)
      DOUBLE PRECISION, INTENT(IN) :: GRAV, R2
      DOUBLE PRECISION :: ZEP, MDTRSW
      !  pressure and temperature 
      BC_Psurf = DLOG(FK/GRAV*(PG + 0.5D0*PR)/max(2D0/3D0,tau_Photosphere))
      ! pressure at 'surface' is determined by local gravity, opacity, and photosphere optical depth
      BC_Tsurf = DLOG(1.0D11*L/(0.75D0*C4PI*CL*R2*PR)) ! effective surface temperature is determined by local P, L, and R
      MDTRSW = 4D-13*(R*L/M)/(CRSN*CLSN/CMSN)/CSY
      IF ( JMOD .EQ. 0 ) THEN
         ZEP = 0D0 ! hold off on extra perturbations until get the first model working
      ELSE
         ZEP = wind_Eta*MDTRSW ! Reimer's wind causes extra mass loss from surface
      END IF
      ! mass BC combines flow from surface shell to next-to-surface shell and flow to/from star to surroundings (e.g., winds)
      BC_Msurf = M_dt + ZEP
      ! the mass BC also includes a term for artificial gain/loss of mass from the surface
      IF ( JMOD .NE. 0 ) BC_Msurf = BC_Msurf - extra_Mdot_param*M/CSY
      END SUBROUTINE Surface_BCs
      
      SUBROUTINE Set_Composition_Vars
      DOUBLE PRECISION :: X_other
      XMG = CMG*CZS; XSI = CSI*CZS; XFE = CFE*CZS ! Mg, Si, and Fe abundances are constants
      X_other = XH + (XHE + (XC + (XO + (XNE + (XMG + XSI + XFE)))))
      XN = 1D0 - X_other
      IF (CZS == 0d0) THEN ! dump the leftovers into He
         XHE = XHE + XN
         RETURN
      END IF
      ! dump the leftovers into N
      IF ( XN .LT. 0D0 ) XN = 0D0
      IF ( LNT .GT. 19.3D0 ) THEN ! no N14 left at very high temps
         XN = 0D0
         X_other = XH + (XHE + (XC + (XO + (XNE + (XSI + XFE)))))
         XMG = 1D0 - X_other
         IF ( XMG .LT. 0D0 ) XMG = 0D0
      END IF
      END SUBROUTINE Set_Composition_Vars
      
      SUBROUTINE Composition_Changes
      ! composition equations
      ! burn_H, burn_He, and burn_Metals
      !       make it possible to inhibit composition changes while keeping energy production
      ! Xi is the fractional adundance for I
      ! XiT is the term (Ri + d_Xi/dt)*d_M/d_k for the difference equation
      ! Ri is the rate of consumption of fuel i in units of reactions/sec
      ! Multiply by number of i baryons to convert from change in
      ! relative number abundance to fractional change in abundance by mass.
      ! H1 change
      XHT = (2.0*burn_H*(RPP+RPC+RPNG+(RPN+RPO)) + DXH/DT)*MK
         ! consumers (note that each reaction consumes 2 protons; hence the factor of 2.0 in XHT)
            ! RPP, equilibrium pp chain      2 p -> 1/2 He4
            ! RPC, C12 + H1 (CN cycle)       C12 (p,positron&neutrino) C13 (p,gamma) N14
            ! RPNG, more probable N14+p      N14 (p,positron&neutrino) N15 (p,gamma) O16
            ! RPN less probable N14+p        N14 (p,positron&neutrino) N15 (p,alpha) C12
            ! RPO, O16 + H1 (ON cycle)       O16 (p,positron&neutrino) O17 (p,alpha) N14
      ! He4 change
      XHET = burn_H*(-0.5*RPP-RPN-RPO)+burn_He*(3.0*R3A+RAC+1.5*RAN+RAO+RANE)
      XHET = XHET + burn_Metals*(-RCC-RCO-2.0*ROO-RGNE-RGMG)
      XHET = (4.0*XHET + DXHE/DT)*MK
         ! producers (note that RPP produces 1/2 He4 and ROO produces 2)
            ! RPP, equilibrium pp chain      2 p -> 1/2 He4
            ! RPN, less probable N14+p       N14 (p,positron&neutrino) N15 (p,alpha) C12
            ! RPO, O16 + H1 (ON cycle)       O16 (p,positron&neutrino) O17 (p,alpha) N14
            ! RCC, C12 + C12                 C12 (C12,alpha&gamma) Ne20
            ! RCO, C12 + O16                 C12 (O16,alpha&gamma) Mg24
            ! ROO, 2 O16                     O16 (O16,alpha&gamma) Si28 (gamma,alpha) Mg24
            ! RGNE, Ne20 decay               Ne20 (gamma,alpha) O16
            ! RGMG, Mg24 decay               Mg24 (gamma,alpha) Ne20
         ! consumers (R3A consumes 3 He4 and RAN consumes 1.5)
            ! R3A, triple alpha process      He4 (alpha) Be8* (alpha,gamma) C12
            ! RAC, He4 + C12                 C12 (alpha,gamma) O16
            ! RAN, He4 + N14                 N14 (alpha,gamma) F18 (1/2 alpha,gamma) Ne20
            ! RAO, He4 + O16                 O16 (alpha,gamma) Ne20
            ! RANE, He4 + Ne20               Ne20 (alpha,gamma) Mg24
      ! C12 change 
      XCT = (12.0*(burn_H*(RPC-RPN)+burn_He*(RAC-R3A)+burn_Metals*(2.0*(RCC+RCCG)+RCO)) + DXC/DT)*MK
         ! producers
            ! RPN, less probable N14+p       N14 (p,positron&neutrino) N15 (p,alpha) C12
            ! R3A, triple alpha process      He4 (alpha) Be8* (alpha,gamma) C12
         ! consumers (RCC and RCCG each consume 2 C12's)
            ! RPC, C12 + H1 (CN cycle)       C12 (p,positron&neutrino) C13 (p,gamma) N14
            ! RAC, He4 + C12                 C12 (alpha,gamma) O16
            ! RCC, C12 + C12                 C12 (C12,alpha&gamma) Ne20
            ! RCCG, C12 + C12                C12 (C12,gamma) Mg24
            ! RCO, C12 + O16                 C12 (O16,alpha&gamma) Mg24
      ! O16 change
      XOT = (16.0*(burn_H*(RPO-RPNG)+burn_He*(RAO-RAC)+burn_Metals*(RCO+2.0*ROO-RGNE)) + DXO/DT)*MK
         ! producers
            ! RPNG, more probable N14+p      N14 (p,positron&neutrino) N15 (p,gamma) O16
            ! RAC, He4 + C12                 C12 (alpha,gamma) O16
            ! RGNE, Ne20 decay               Ne20 (gamma,alpha) O16
         ! consumers (ROO consumes 2 O16's)
            ! RPO, O16 + H1 (ON cycle)       O16 (p,positron&neutrino) O17 (p,alpha) N14
            ! RAO, He4 + O16                 O16 (alpha,gamma) Ne20
            ! RCO, C12 + O16                 O16 (C12,alpha&gamma) Mg24
            ! ROO, 2 O16                     O16 (O16,alpha&gamma) Si28 (gamma,alpha) Mg24
      ! Ne20 change
      XNET = (20.0*(burn_He*(RANE-RAN-RAO)+burn_Metals*(RGNE-RGMG-RCC)) + DXNE/DT)*MK
         ! producers
            ! RAN, He4 + N14                 N14 (alpha,gamma) F18 (1/2 alpha,gamma) Ne20
            ! RAO, He4 + O16                 O16 (alpha,gamma) Ne20
            ! RGMG, Mg24 decay               Mg24 (gamma,alpha) Ne20
            ! RCC, C12 + C12                 C12 (C12,alpha&gamma) Ne20
         ! consumers
            ! RANE, He4 + Ne20               Ne20 (alpha,gamma) Mg24
            ! RGNE, Ne20 decay               Ne20 (gamma,alpha) O16
      END SUBROUTINE Composition_Changes
      
      SUBROUTINE Save_SV_VARs ! put the current values of the state variables into the array
      SV(SV_LNP)=LNP; SV(SV_LNRHO)=LNRHO; SV(SV_U)=U; SV(SV_P)=P; SV(SV_RHO)=RHO; SV(SV_FK)=FK
      SV(SV_T)=T; SV(SV_SF)=SF; SV(SV_ST)=ST; SV(SV_GRADA)=GRADA; SV(SV_SCP)=SCP
      SV(SV_XHI)=XHI; SV(SV_S)=S; SV(SV_PR)=PR; SV(SV_PG)=PG
      SV(SV_EN)=EN; SV(SV_RPP)=RPP; SV(SV_R33)=R33; SV(SV_R34)=R34; SV(SV_RBE)=RBE
      SV(SV_RBP)=RBP; SV(SV_RPC)=RPC; SV(SV_RPN)=RPN; SV(SV_RPO)=RPO; SV(SV_R3A)=R3A
      SV(SV_RAC)=RAC; SV(SV_RAN)=RAN; SV(SV_RAO)=RAO; SV(SV_RANE)=RANE; SV(SV_RCC)=RCC
      SV(SV_RCO)=RCO; SV(SV_ROO)=ROO; SV(SV_RGNE)=RGNE; SV(SV_RGMG)=RGMG; SV(SV_RCCG)=RCCG
      SV(SV_RPNG)=RPNG; SV(SV_EX)=EX; SV(SV_ENX)=ENX; SV(SV_ENG)=ENG
      SV(SV_MK)=MK; SV(SV_GRADR)=GRADR; SV(SV_GRAD)=GRAD
      SV(SV_R)=R; SV(SV_MESH_Q)=MESH_Q; SV(SV_QM)=QM; SV(SV_WL)=WL
      SV(SV_WCV)=WCV; SV(SV_HP)=HP; SV(SV_XN)=XN; SV(SV_XMG)=XMG; SV(SV_XSI)=XSI
      SV(SV_XFE)=XFE; SV(SV_N1)=NA(N_H); SV(SV_N4)=NA(N_HE); SV(SV_N12)=NA(N_C)
      SV(SV_N14)=NA(N_N); SV(SV_O16)=NA(N_O); SV(SV_N20)=NA(N_NE); SV(SV_N24)=NA(N_MG)
      SV(SV_N28)=NA(N_SI); SV(SV_N56)=NA(N_FE); SV(SV_NE)=NE; SV(SV_NE1)=NE1
      SV(SV_NI)=NI; SV(SV_NZZ)=NZZ; SV(SV_AVM)=AVM; SV(SV_P_ION)=P_ION; SV(SV_PCORR)=PCORR
      END SUBROUTINE Save_SV_VARs
      
      SUBROUTINE Get_SV_VARs ! unload the cached state variable values
      LNP=SV(SV_LNP); LNRHO=SV(SV_LNRHO); U=SV(SV_U); P=SV(SV_P); RHO=SV(SV_RHO); FK=SV(SV_FK)
      T=SV(SV_T); SF=SV(SV_SF); ST=SV(SV_ST); GRADA=SV(SV_GRADA); SCP=SV(SV_SCP)
      XHI=SV(SV_XHI); S=SV(SV_S); PR=SV(SV_PR); PG=SV(SV_PG)
      EN=SV(SV_EN); RPP=SV(SV_RPP); R33=SV(SV_R33); R34=SV(SV_R34); RBE=SV(SV_RBE)
      RBP=SV(SV_RBP); RPC=SV(SV_RPC); RPN=SV(SV_RPN); RPO=SV(SV_RPO); R3A=SV(SV_R3A)
      RAC=SV(SV_RAC); RAN=SV(SV_RAN); RAO=SV(SV_RAO); RANE=SV(SV_RANE); RCC=SV(SV_RCC); RCO=SV(SV_RCO)
      ROO=SV(SV_ROO); RGNE=SV(SV_RGNE); RGMG=SV(SV_RGMG); RCCG=SV(SV_RCCG); RPNG=SV(SV_RPNG)
      EX=SV(SV_EX); ENX=SV(SV_ENX); ENG=SV(SV_ENG)
      END SUBROUTINE Get_SV_VARs

      SUBROUTINE Unload_SV_VARs(K,SV_VALs) ! unload the cached state variable values
      INTEGER, INTENT(IN) :: K
      DOUBLE PRECISION, INTENT(IN) :: SV_VALs(N_SHELLS,NSVARS)
      LNP=SV_VALs(K,SV_LNP); LNRHO=SV_VALs(K,SV_LNRHO); U=SV_VALs(K,SV_U); P=SV_VALs(K,SV_P)
      RHO=SV_VALs(K,SV_RHO); FK=SV_VALs(K,SV_FK)
      T=SV_VALs(K,SV_T); SF=SV_VALs(K,SV_SF); ST=SV_VALs(K,SV_ST)
      GRADA=SV_VALs(K,SV_GRADA); SCP=SV_VALs(K,SV_SCP)
      XHI=SV_VALs(K,SV_XHI); S=SV_VALs(K,SV_S); PR=SV_VALs(K,SV_PR)
      PG=SV_VALs(K,SV_PG); EN=SV_VALs(K,SV_EN); RPP=SV_VALs(K,SV_RPP); R33=SV_VALs(K,SV_R33); R34=SV_VALs(K,SV_R34)
      RBE=SV_VALs(K,SV_RBE); RBP=SV_VALs(K,SV_RBP); RPC=SV_VALs(K,SV_RPC); RPN=SV_VALs(K,SV_RPN); RPO=SV_VALs(K,SV_RPO)
      R3A=SV_VALs(K,SV_R3A); RAC=SV_VALs(K,SV_RAC); RAN=SV_VALs(K,SV_RAN); RAO=SV_VALs(K,SV_RAO); RANE=SV_VALs(K,SV_RANE)
      RCC=SV_VALs(K,SV_RCC); RCO=SV_VALs(K,SV_RCO); ROO=SV_VALs(K,SV_ROO); RGNE=SV_VALs(K,SV_RGNE); RGMG=SV_VALs(K,SV_RGMG)
      RCCG=SV_VALs(K,SV_RCCG); RPNG=SV_VALs(K,SV_RPNG); EX=SV_VALs(K,SV_EX); ENX=SV_VALs(K,SV_ENX)
      MK=SV_VALs(K,SV_MK); GRADR=SV_VALs(K,SV_GRADR); GRAD=SV_VALs(K,SV_GRAD); ENG=SV_VALs(K,SV_ENG)
      R=SV_VALs(K,SV_R); MESH_Q=SV_VALs(K,SV_MESH_Q); QM=SV_VALs(K,SV_QM); WL=SV_VALs(K,SV_WL)
      WCV=SV_VALs(K,SV_WCV); HP=SV_VALs(K,SV_HP); XN=SV_VALs(K,SV_XN); XMG=SV_VALs(K,SV_XMG); XSI=SV_VALs(K,SV_XSI)
      XFE=SV_VALs(K,SV_XFE); NA(N_H)=SV_VALs(K,SV_N1); NA(N_HE)=SV_VALs(K,SV_N4); NA(N_C)=SV_VALs(K,SV_N12)
      NA(N_N)=SV_VALs(K,SV_N14); NA(N_O)=SV_VALs(K,SV_O16); NA(N_NE)=SV_VALs(K,SV_N20); NA(N_MG)=SV_VALs(K,SV_N24)
      NA(N_SI)=SV_VALs(K,SV_N28); NA(N_FE)=SV_VALs(K,SV_N56); NE=SV_VALs(K,SV_NE); NE1=SV_VALs(K,SV_NE1)
      NI=SV_VALs(K,SV_NI); NZZ=SV_VALs(K,SV_NZZ); AVM=SV_VALs(K,SV_AVM)
      P_ION=SV_VALs(K,SV_P_ION); PCORR=SV_VALs(K,SV_PCORR)
      END SUBROUTINE Unload_SV_VARs

      SUBROUTINE STATEL ( JI )
      ! Either compute new EoS values, or recover old values from store
      ! input JI determines use of cached values
         ! JI < 0    compute EoS and don't save
         ! JI = 0    compute EoS and save
         ! JI > 0    action depends on value of JI
            ! when creating numeric partial derivatives for the Jacobian matrix in DIFRNS
            ! JI is the number of the variable whose value has been varied since we cached results
            ! if the variable doesn't effect the EoS values, then we can use the cached values
               ! use cached values when change any of the variables M, Q_dk, AR, or L
            ! if the variable does effect the EoS values, then must recompute and don't save in cache.
               ! compute and don't save for variables LOGF, LOGT, XO, XH, XHE, XC, XN, and XNE
      ! inputs LOGF and LOGT are args for STATEF and NUCRAT
      USE ez_state
      USE ez_nuclear
      USE ez_opacity
      USE ez_vcool
      USE star_controls
      INTEGER, INTENT(IN) :: JI
      INTEGER :: IFN(14) = (/2, 2, 0, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2 /)
      DOUBLE PRECISION :: XA(NEL)
      ! IFN array has codes which determine the use of cached values.  codes are 0, 1, and 2.
      ! 1 means use cached values; not = 1 means compute new value.
      ! 0 means save new value in cache; 2 means don't save.
      IF ( IFN(JI+3) .NE. 1 ) THEN ! don't use cache
         XA(N_H) = XH; XA(N_HE) = XHE; XA(N_C) = XC; XA(N_N) = XN; XA(N_O) = XO; 
         XA(N_NE) = XNE; XA(N_MG) = XMG; XA(N_SI) = XSI; XA(N_FE) = XFE
         CALL STATEF ( LNF, LNT, T, XA, ionization_level )
         EN = Neutrino_Cooling(RHO, T, LNRHO, LNT, MUE, XC, XO, XHE, eps_plasma, eps_brem, eps_pair, eps_photo)
         EN = EN * neutrino_Cooling_Factor
         FK = Opacity(LNT/CLN, LNRHO/CLN, NA(N_H)*CAN(N_H)/AVM, NA(N_HE)*CAN(N_HE)/AVM)
         XHI = 4D0*CL*PR/(FK*RHO*RHO*SCP*T)
         EX = NUCRAT ( LNT, ZT, NA, NI, NE, AVM, RHO, ENX ) ! ZT from STATEF
      ELSE
         CALL Get_SV_VARs
      END IF
      END SUBROUTINE STATEL
      
      END MODULE ez_shell

