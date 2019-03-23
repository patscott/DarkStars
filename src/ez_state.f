      MODULE ez_state
      USE ez_data
      USE ez_state_data
      USE ez_ionization_data
      IMPLICIT NONE
      
      CHARACTER (LEN=strlen) :: var_name

      CONTAINS
      
      SUBROUTINE REPORT ( val )
      DOUBLE PRECISION, INTENT(IN) :: val
      WRITE(44,'(A20,E30.14)') var_name, val
      END SUBROUTINE REPORT
      
      SUBROUTINE STATEF ( LN_F, LN_T, T, XA, ions )
      ! input LN_F is log(F).  F is the electron degeneracy variable.
      ! input LN_T is log(T).  T is temperature in K.
      ! input XA is abundance array of mass fractions.
      USE ez_ionization
      DOUBLE PRECISION, INTENT(IN) :: LN_F, LN_T, T, XA(NEL)
      INTEGER, INTENT(IN) :: ions
      DOUBLE PRECISION :: DL, RT, DPB, DB
      DOUBLE PRECISION :: Q, DSBT, DUB, DPBT, DSB, B, DSF, DST, TCR, P0, T4
      DOUBLE PRECISION :: DPA, DPLN_T, DC, DPLN_F, DUA, DSA, DSLN_T, DSLN_F, UF, WF
      DOUBLE PRECISION :: F, G, UE, SEF, SET, DSBF, DPBF, PF, PT, RF
      DOUBLE PRECISION :: DE, TI, DV, DVF, DVT, DS_CORR, PF_CORR, PT_CORR, DSF_CORR, DST_CORR, U_CORR, CORR_FAC
      DOUBLE PRECISION :: H2, H2F, H2T, SI, SIF, SIT, UI, ZH2T, NEF, NET, NI1, EN
      DOUBLE PRECISION, PARAMETER :: EPSI = 1D-30
      INTEGER :: I
      
      CORR_FAC = 1D0
      
      F = exp(LN_F)
      UF = F/(1D0 + F)
      WF = sqrt(1D0 + F)
      PSI = LN_F + 2D0*(WF - DLOG(1D0 + WF))
      G = CTE*T*WF
      CALL FDIRAC ( F, G )
      PE = G*PE
      PET = PET + 1D0
      PEF = PEF + 0.5D0*UF
      QE = QE/(RE*WF)   ! Q_star / (sqrt(1+f) * rho_star)
      SE = QE + 2D0*WF - PSI ! Q_star / (sqrt(1+f) rho_star) + 2 sqrt(1+f) - PSI
      SEF = QE*(QEF - REF - 0.5D0*UF) - 1D0/WF
      SET = QE*(QET - RET)
      UE = SE + PSI - PE/(RE*CTE*T)
      ! Evaluate some quantities that do not depend on the state of ionization
      ! Thus when set N_i = X_i / CBN(i), get n_i / (rho * N_avo).
      AVM = 0D0; NI = 0D0; NE = 0D0; NZZ = 0D0
      DO I = NEL,1,-1
         NA(I) = XA(I)/CBN(I)
         AVM = AVM + CAN(I)*NA(I)
         NI = NI + NA(I)
         NE = NE + KZN(I)*NA(I)
         NZZ = NZZ + KZN(I)*KZN(I)*NA(I)
      END DO
      ! CAN's are atomic weights
      ! CAN*X/CBN converts to relative abundances by atomic weight
      ! AVM is sum of these, so dividing CAN*X/CBN's by AVM gives new fractional abundance by mass
      ! NI is sum of N's, so dividing N's by NI gives fractional abundance by number 
      TI = CEVB/T
      DE = RE*CD

      IF (CORR_FAC .GT. 1D-6) THEN
         CALL PRESSI ( 1, TI, PSI, DE, REF, RET, F, DC, DVT, DVF, DPA, DPLN_T, DPLN_F, DSA, DSLN_T, DSLN_F, DUA )
         DC = CORR_FAC*DC
         DVT = CORR_FAC*DVT
         DVF = CORR_FAC*DVF
      ELSE
         DC = 0D0
         DVT = 0D0
         DVF = 0D0
      END IF

      DV = DC - PSI ! DV is the corrected PSI
      DVF = DVF - WF ! DVF is d_DV / d_F
      ! calculate ionization
      CALL EoS_Ionization(NE1,NA,T,RET,EN,NI,REF, DE, TI, DV,DVF,DVT, H2,H2F,H2T, SI,SIF,SIT, UI, ZH2T,NEF,NET,NI1,ions)
      DB = EN*DE
      DL = DLOG(DB)
      RHO = DB*AVM ! RHO is the mass density in g/cm^3
      LNRHO = DLOG(RHO)
      RT = RET - EN*NET
      RF = REF - EN*NEF
      DE = DB*NE
      ! pressure terms
      PE = CB*PE
      TCR = T*CR
      P0 = TCR*DB
      P_ION = NI1*P0
      T4 = T*T*T*T
      PR = CRAD*T4/3D0
      B = 4D0*PR/P0
      
      IF (CORR_FAC .GT. 1D-6) THEN
         CALL PRESSI ( 0, TI, PSI, DE, RF, RT, F, DC, DVT, DVF, DPB, DPBT, DPBF, DSB, DSBT, DSBF, DUB )
         ! second call on PRESSI: provides corrections to pressure, entropy, energy and their derivatives
         ! so that pressure ionization effects go away in full ionization limit (see PTEH, eq 28)
         ! use the following: DPB, DPBT, DPBF, DSB, DSBT, DSBF, DUB
         PCORR = CORR_FAC*TCR*(DPA - DPB)
         PT_CORR = CORR_FAC*TCR*(DPLN_T - DPBT)
         PF_CORR = CORR_FAC*TCR*(DPLN_F - DPBF)
         DS_CORR = CORR_FAC*(DSA*NE1 - DSB*NE)
         DSF_CORR = CORR_FAC*(NE1*DSLN_F - NE*DSBF)
         DST_CORR = CORR_FAC*(NE1*DSLN_T - NE*DSBT)
         U_CORR = CORR_FAC*(DUA*NE1 - DUB*NE)
      ELSE
         PCORR = 0D0
         PT_CORR = 0D0
         PF_CORR = 0D0
         DS_CORR = 0D0
         DSF_CORR = 0D0
         DST_CORR = 0D0
         U_CORR = 0D0
      END IF

      PG = PE + P_ION + PCORR 
      P = PG + PR
      LNP = DLOG(P)
      PF = (PE*PEF      + P_ION*RF - H2F*P0 + PF_CORR)/P
      PT = (PE*PET + P_ION + P_ION*RT - H2T*P0 + PT_CORR + PR*4D0)/P
      ! specific entropy, in erg/g/K
      DSF = NEF*DSA + DSF_CORR - RF*B
      DST = NET*DSA + DST_CORR - (RT - 3D0)*B
      SF = CR*(-NI1*RF           + NEF*SE + NE1*SEF + SIF + DSF)/AVM ! SF is d_S / d_ln(F)
      ST = CR*( NI1*(1.5D0 - RT) + NET*SE + NE1*SET + SIT + DST)/AVM ! ST is d_S / d_ln(T)
      S = CR*(SE*NE1 + DS_CORR + B + (1.5D0*LN_T - DL + 2.5D0 - DLOG(CEN))*NI1 + SI)/AVM
      ! specific internal energy, in erg/g
      U = TCR*(UE*NE1 + U_CORR + 1.5D0*NI1 + ZH2T*H2 + 0.75D0*B + TI*UI)/AVM
      ! other thermodynamic quantities
      Q = PT*SF - PF*ST
      SCP = -Q/PF
      GRADA = SF/Q
      GAMMA1 = Q/(RT*SF-RF*ST) ! adiabatic exponent: (d_ln(P)/d_ln(rho)) at constant entropy.
      ZT = DSQRT(DABS((NE1*REF/WF + NZZ)/NI1))
      END SUBROUTINE STATEF

      SUBROUTINE FDIRAC ( F, G )
      ! Evaluate Fermi-Dirac integrals (Eggleton, Faulkner & Flannery 1973).
      ! to get electron gas density, pressure, entropy, and internal energy.
      DOUBLE PRECISION, INTENT(IN) :: F, G
      DOUBLE PRECISION :: WW, WV, VG, VF, UG, FDF, UF
      DOUBLE PRECISION :: VX11, VX12, VX13, VX14, VX21, VX22, VX23, VX31, VX32, VX33, VX41
      DOUBLE PRECISION :: VW11, VW12, VW13, VW14, VW21, VW22, VW23, VW31, VW32, VW41
      INTEGER :: n, m
      DOUBLE PRECISION :: FF(4), GG(4), C_rho(4,4), C_P(4,4), C_Q(4,4)
         ! C_rho is 4x4 matrix of coefficients for rho_star
         ! C_P is 4x4 matrix of coefficients for P_star
         ! C_Q is 4x4 matrix of coefficients for Q_star
         
      C_rho(:,1) = (/  2.315472D0,  7.128660D0,  7.504998D0,  2.665350D0 /) 
      C_rho(:,2) = (/  7.837752D0, 23.507934D0, 23.311317D0,  7.987465D0 /) 
      C_rho(:,3) = (/  9.215560D0, 26.834068D0, 25.082745D0,  8.020509D0 /) 
      C_rho(:,4) = (/  3.693280D0, 10.333176D0,  9.168960D0,  2.668248D0 /)
      
      C_P(:,1) = (/  2.315472D0,  6.748104D0,  6.564912D0,  2.132280D0 /) 
      C_P(:,2) = (/  7.837752D0, 21.439740D0, 19.080088D0,  5.478100D0 /) 
      C_P(:,3) = (/  9.215560D0, 23.551504D0, 19.015888D0,  4.679944D0 /) 
      C_P(:,4) = (/  3.693280D0,  8.859868D0,  6.500712D0,  1.334124D0 /)
      
      C_Q(:,1) = (/  1.157736D0,  3.770676D0,  4.015224D0,  1.402284D0 /) 
      C_Q(:,2) = (/  8.283420D0, 26.184486D0, 28.211372D0, 10.310306D0 /) 
      C_Q(:,3) = (/ 14.755480D0, 45.031658D0, 46.909420D0, 16.633242D0 /) 
      C_Q(:,4) = (/  7.386560D0, 22.159680D0, 22.438048D0,  7.664928D0 /)
      
      VF = 1D0/(1D0 + F)      ! 1/(1+f)
      VG = 1D0/(1D0 + G)      ! 1/(1+g)
      UF = F*VF                  ! f/(1+f)
      UG = G*VG                  ! g/(1+g)
      FDF = G*G + G              ! g*(1+g)
      FDF = UF*FDF*DSQRT(FDF)    ! f/(1+f)*(g/(1+g))^(3/2)
      FF(1) = 1D0; GG(1) = 1D0
      FF(2) = F*FF(1); GG(2) = G*GG(1); FDF = FDF*VF*VG
      FF(3) = F*FF(2); GG(3) = G*GG(2); FDF = FDF*VF*VG
      FF(4) = F*FF(3); GG(4) = G*GG(3); FDF = FDF*VF*VG
      ! FF(I) = f^(I-1); GG(I) = g^(I-1)
      ! FDF = prefactor of SUM for rho_star: = f/(1+f)*(g/(1+g))^(3/2)/(1+f)^3/(1+g)^3
      ! VX<*,IM> gets (IM-1)st derivative wrt ln(F)
      ! VX<IL,*> gets (IL-1)st derivative wrt ln(T)
      ! first do rho_star.  we need 1st, 2nd, and 3rd order derivatives.
      VX11 = 0D0; VX12 = 0D0; VX13 = 0D0; VX14 = 0D0; VX21 = 0D0; VX22 = 0D0; VX23 = 0D0
      VX31 = 0D0; VX32 = 0D0; VX33 = 0D0; VX41 = 0D0
      DO n = 0, 3 ! sum n=0,3 
         DO m = 0, 3 ! sum m=0,3
            VW11 = C_rho(m+1, n+1)*GG(n+1)*FF(m+1) ! term for f^m * g^n
            VW21 = n*VW11
            VW12 = m*VW11
            VW13 = m*VW12
            VW14 = m*VW13
            VW31 = n*VW21
            VW22 = m*VW21
            VW23 = m*VW22
            VW41 = n*VW31
            VW32 = m*VW31
            VX11 = VX11 + VW11
            VX12 = VX12 + VW12
            VX13 = VX13 + VW13
            VX14 = VX14 + VW14
            VX21 = VX21 + VW21
            VX22 = VX22 + VW22
            VX23 = VX23 + VW23
            VX31 = VX31 + VW31
            VX32 = VX32 + VW32
            VX41 = VX41 + VW41
         END DO
      END DO
      WV = 1D0/VX11
      VX12 = VX12*WV
      VX13 = VX13*WV
      VX14 = VX14*WV
      VX21 = VX21*WV
      VX22 = VX22*WV
      VX23 = VX23*WV
      VX31 = VX31*WV
      VX32 = VX32*WV
      VX41 = VX41*WV
      RE = FDF*VX11
      RET = VX21 + 1.5D0 - 1.5D0*UG ! d_rho_star / d_ln(T)
      WW = 0.5D0*RET - 4D0 
      REF = VX12 + 1D0 + WW*UF      ! d_rho_star / d_ln(F)
      ! also get the second and third order derivatives of density, needed in PRESSI
      XTT = VX31 - VX21**2 - 1.5D0*UG*VG
      XFT = VX22 - VX21*VX12 + 0.5D0*UF*XTT
      XFF = VX13 - VX12**2 + UF*(XFT + WW*VF - 0.25D0*XTT*UF)
      XTTT = VX41 + VX21*(2D0*VX21**2 - 3D0*VX31) - 1.5D0*(1.0-G)*UG*VG**2
      XFTT = VX32 - 2D0*VX21*VX22 - VX12*(VX31 - 2D0*VX21**2) + 0.5D0*UF*XTTT
      XFFT = VX23 - 2D0*VX22*VX12-VX21*(VX13-2D0*VX12**2)+UF*(XFTT+0.5D0*VF*XTT-0.25D0*UF*XTTT)
      XFFF = VX14 + VX12*(2D0*VX12**2-3D0*VX13)
      XFFF = XFFF + UF*(1.5D0*(XFFT+VF*XFT)-UF*(0.75D0*(XFTT+VF*XTT)-0.125D0*UF*XTTT)+WW*(1D0-F)*VF**2)
      ! now P_star.  we only need 1st order derivatives.
      VX11 = 0D0; VX12 = 0D0; VX21 = 0D0; VX22 = 0D0
      DO n = 0, 3 ! sum n=0,3
         VW11 = C_P(0+1, n+1)*GG(n+1)*FF(1) ! term for f^0 * g^n
         VW21 = n * VW11
         VX11 = VX11 + VW11
         VX21 = VX21 + VW21
         VW11 = C_P(2, n+1)*GG(n+1)*FF(2) ! term for f^1 * g^n
         VW21 = n * VW11
         VW12 = VW11
         VX11 = VX11 + VW11
         VX12 = VX12 + VW12
         VX21 = VX21 + VW21
         VW11 = C_P(3, n+1)*GG(n+1)*FF(3) ! term for f^2 * g^n
         VW21 = n * VW11
         VW12 = 2 * VW11
         VX11 = VX11 + VW11
         VX12 = VX12 + VW12
         VX21 = VX21 + VW21
         VW11 = C_P(4, n+1)*GG(n+1)*FF(4) ! term for f^3 * g^n
         VW21 = n * VW11
         VW12 = 3 * VW11
         VX11 = VX11 + VW11
         VX12 = VX12 + VW12
         VX21 = VX21 + VW21
      END DO
      WV = 1D0/VX11
      VX12 = VX12*WV
      VX21 = VX21*WV
      PE = FDF*VX11
      PET = VX21 + 1.5D0 - 1.5D0*UG ! d_P_star / d_ln(T)
      WW = 0.5D0*PET - 4D0 
      PEF = VX12 + 1D0 + WW*UF      ! d_P_star / d_ln(F)
      ! now Q_star.  we only need 1st order derivatives.
      VX11 = 0D0; VX12 = 0D0; VX21 = 0D0; VX22 = 0D0
      DO n = 0, 3 ! sum n=0,3
         VW11 = C_Q(0+1, n+1)*GG(n+1)*FF(1) ! term for f^0 * g^n
         VW21 = n * VW11
         VX11 = VX11 + VW11
         VX21 = VX21 + VW21
         VW11 = C_Q(2, n+1)*GG(n+1)*FF(2) ! term for f^1 * g^n
         VW21 = n * VW11
         VW12 = VW11
         VX11 = VX11 + VW11
         VX12 = VX12 + VW12
         VX21 = VX21 + VW21
         VW11 = C_Q(3, n+1)*GG(n+1)*FF(3) ! term for f^2 * g^n
         VW21 = n * VW11
         VW12 = 2 * VW11
         VX11 = VX11 + VW11
         VX12 = VX12 + VW12
         VX21 = VX21 + VW21
         VW11 = C_Q(4, n+1)*GG(n+1)*FF(4) ! term for f^3 * g^n
         VW21 = n * VW11
         VW12 = 3 * VW11
         VX11 = VX11 + VW11
         VX12 = VX12 + VW12
         VX21 = VX21 + VW21
      END DO
      WV = 1D0/VX11
      VX12 = VX12*WV
      VX21 = VX21*WV
      QE = FDF*VX11
      QET = VX21 + 1.5D0 - 1.5D0*UG ! d_Q_star / d_ln(T)
      WW = 0.5D0*QET - 4D0 
      QEF= VX12 + 1D0 + WW*UF    ! d_Q_star / d_ln(F)
      END SUBROUTINE FDIRAC

      SUBROUTINE PRESSI ( IPR, TI, PSI, XI, XF_IN, XT, F, DC, DCT, DCF, DP, DPT, DPF, DS, DST, DSF, DU )
      ! PRESSI gives a model for pressure ionization and a model for Coulomb interactions.
      ! Returns corrections to the electron chemical potential, pressure, entropy, and internal energy.
      ! Computes effect (`pressure ionisation') on EOS of a free energy contribution 
      ! delta(F) = -R.T.Ne.G(X,Y), where X = XI = Ne/V = ne, and Y = YI = chiH/(R.T).
      ! delta(el. chem. pot.) = -DC = 1/(R.T) dF/dNe, delta(P) = R.T.DP = -dF/dV, 
      ! delta(S) = R.Ne.DS = -dF/dT. Works firstly with indep. vbles X, Y (called 
      ! XI, YI), which are effectively Ne/V and T, but then has to get derivatives
      ! of DC, DP, DS w.r.t. indep vbles F, T. Also does Coulomb interaction.
      ! input IPR    1 for first call, 0 for second call
      ! input TI     inverse of temperature in the form of (ergs per electron volt)/(kb*T)
      ! input PSI    electron degeneracy parameter
      ! input XI     number of electrons per unit volume (free e's first call; all e's second call?)
      ! input XF     d_rho* / d_ln(F)
      ! input XT     d_rho* / d_ln(T)
      ! input F      degeneracy variable F
      ! output DC    correction for electron chemical potential
      ! output DCT   d_DC / d_ln(T)
      ! output DCF   d_DC / d_ln(F)
      ! output DP    correction for pressure
      ! output DPT   d_DP / d_ln(T)
      ! output DPF   d_DP / d_ln(F)
      ! output DS    correction for entropy
      ! output DST   d_DS / d_ln(T)
      ! output DSF   d_DS / d_ln(F)
      ! output DU    correction for internal energy
      INTEGER, INTENT(IN) :: IPR
      DOUBLE PRECISION, INTENT(IN) :: TI, PSI, XI, XF_IN, XT, F
      DOUBLE PRECISION, INTENT(OUT) :: DC, DCT, DCF, DP, DPT, DPF, DS, DST, DSF, DU
      DOUBLE PRECISION :: XF, CA1, CA2, CA3, CP1, CP2, CP3, CP4
      DOUBLE PRECISION :: THXX, THXY, WXX, WXY, WW, THYY, THC, DGDYY, AFF, DGDXX, DGDXY, THTT
      DOUBLE PRECISION :: TOF, THFF, GAMX, GAMY, RTHC, WY, GAMYY, WWT, GAMXX, GAMXY, EEG, BEG, SRBE, GE
      DOUBLE PRECISION :: DGDG, DGDGG, GEG, AF, THF, FF2, THX, THY, THT, RXF, WX, CC2, YI, EE, FF1, TH
      DOUBLE PRECISION :: AA, BB, GI, BXY, BYY, BXX, DGDX, DGDY, CX, CXX, W2, PSXX, PSX, PSY, EEX, BX
      DOUBLE PRECISION :: PSXY, PSYY, THFT, BY, BEGG
      DOUBLE PRECISION :: CBRT, VX
      CBRT(VX) = DEXP(C3RD*DLOG(VX)) ! statement function for cube root
      VX = 0 ! for the Absoft compiler
      CA1 = 0.89752D0; CA2 = 0.768D0; CA3 = 0.208D0
      CP1 = 3D0; CP2 = 0.35D0; CP3 = 2D0; CP4 = 3.0D-2
      ! Pressure ionization
      XF = XF_IN
      IF (XF .LT. 1D-6) XF = 1D-6
      YI = 13.595D0*TI  ! YI is (H ionization energy)/(kb*T)
      EE = 1D0/(1D0 + XI/CP4)
      WX = (CP1/XI)**CP2
      CC2 = DEXP(-WX) ! CC2 is exponential factor for pressure ionization
      BB = YI + PSI - CP3*DLOG(EE)
      GI = CC2*BB
         ! extra free en. is -R.T.Ne.GI
         ! It OUGHT to have 3psi/5 for psi, but that works worse!
      FF1 = 1D0/(F + 1D0)
      ! AA = dlnf/dpsi; THeta = dlnne/dpsi; also needed in Coulomb corrn;
      ! f, T derivs, then X, Y.
      AA = DSQRT(FF1)
      TH = AA*XF
      FF2 = -0.5D0*F*FF1
      AF = AA*FF2
      THF = AF*XF + AA*XFF
      THT = AA*XFT
      RXF = 1D0/XF
      THX = THF*RXF
      THY = THX*XT - THT
      ! first and second derivatives dPSI/dlnX ... d2PSI/dlnY2
      PSX = 1D0/TH
      PSY = XT*PSX
      W2 = PSX*PSX
      PSXX = -THX*W2
      PSXY = -THY*W2
      PSYY = PSXY*XT + (XFT*XT*RXF - XTT)*PSX
      ! derivative -dlnEE/dlnX; -d2lnEE/dlnX2 = -EE*dlnEE/dlnX
      EEX = CP3*(XI/CP4)*EE
      ! derivatives of BB
      BX = PSX + EEX
      BY = YI + PSY
      BXX = PSXX + EE*EEX
      BXY = PSXY
      BYY = YI + PSYY
      ! derivatives of CC2
      CX = CP2*WX*CC2
      CXX = CP2*CX*(WX - 1D0)
      ! derivatives of GI
      DGDX = CC2*BX + CX*BB
      DGDY = CC2*BY
      DGDXX = CC2*BXX + 2D0*CX*BX + CXX*BB
      DGDXY = CC2*BXY + CX*BY
      DGDYY = CC2*BYY
      IF ( IPR.EQ.1 ) THEN ! do this on first call only
         ! Coulomb interaction.  further derivatives of AA, THeta
         AFF = 3D0*AF*FF2 + AF
         THFF = AFF*XF + 2D0*AF*XFF + AA*XFFF
         THFT = AF*XFT + AA*XFFT 
         THTT = AA*XFTT 
         ! d2TH/dlnX2, etc
         TOF = XT*RXF
         WXX = THFF - THX*XFF
            WXY = THX*XFT - THFT
         THXX = WXX*RXF*RXF
         THXY = WXY*RXF + THXX*XT
         THYY = TOF*(TOF*WXX + 2D0*WXY) + THTT - THX*XTT
         ! GAM is the plasma interaction parameter. Note that THC = ZT**2*NI/NE
         THC = TH + DABS(NZZ/NE)
         WW = 1D0/CA2
         GAM = CBRT(XI*(CPL*NE/NI)**2*C3RD)*TI/CEVB*THC
         ! new BB and EE, and their derivatives
         BB = (CA1*DSQRT(3D0/GAM))**WW
         EE = (GAM/(GAM + CA3))**WW
         SRBE = 1D0/(EE + BB)
         ! further addition GE to free en; adds to previous GI.
         ! GE = GE(GAMma), GAM = const. * X**0.33 * Y * (THeta(X,Y) + const)
         GE = (NI/NE)*CA1*GAM*SRBE**CA2
         EEG = CA3/(GAM + CA3)
         BEG = (EEG*EE - 0.5D0*BB)*SRBE
         BEGG = (EEG*EEG*(1D0 - GAM*CA2/CA3)*EE + 0.25D0*BB)*SRBE
         GEG = 1D0 - BEG
         DGDG = GE*GEG
         DGDGG = GE*(GEG*GEG + (BEG*BEG - BEGG)*WW)
         RTHC = 1D0/THC
         WX = THX*RTHC
         WY = THY*RTHC
         ! dlnGAM/dlnX, etc
         GAMX = C3RD + WX
         GAMY = 1D0 + WY
         GAMXX = THXX*RTHC - WX*WX
         GAMXY = THXY*RTHC - WX*WY
         GAMYY = THYY*RTHC - WY*WY
         GI = GI + GE
         ! derivatives w.r.t. X, Y; in effect w.r.t. Ne, V, and T, since X = Ne/V
         DGDX = DGDX + DGDG*GAMX
         DGDY = DGDY + DGDG*GAMY
         DGDXX = DGDXX + DGDGG*GAMX**2 + DGDG*GAMXX
         DGDXY = DGDXY + DGDGG*GAMX*GAMY + DGDG*GAMXY
         DGDYY = DGDYY + DGDGG*GAMY**2 + DGDG*GAMYY
      END IF
      ! evaluate changes to electron chemical potential (-DC), pressure (R.T.DP), and
      ! entropy (R.Ne.DS), and their derivatives w.r.t. log(F) and log(T)
      DC = DGDX + GI
      WW = DGDXX + DGDX
      WWT = WW*XT
      DCF = WW*XF
      DCT = WWT - DGDXY - DGDY
      DP = -XI*DGDX
      DPF = -XI*DCF
      DPT = XI*(DGDXY - WWT) + DP
      DS = GI - DGDY
      DST = DGDX - DGDXY
      DSF = DST*XF
      DST = DST*XT - DGDY + DGDYY
      DU = -DGDY
      END SUBROUTINE PRESSI

      END MODULE ez_state

