      MODULE ez_ionization
      USE ez_ionization_data
      IMPLICIT NONE
            
      CONTAINS
      
      SUBROUTINE EoS_Ionization(NE1,NA,T,RET,EN,NI,REF, DE, TI, DV,DVF,DVT, H2,H2F,H2T, SI,SIF,SIT, UI, ZH2T,NEF,NET,NI1,KION)
      DOUBLE PRECISION, INTENT(IN) :: NA(NEL), T, RET, NI, REF, DE, TI, DV, DVF, DVT
      INTEGER, INTENT(IN) :: KION
      DOUBLE PRECISION, INTENT(OUT) :: EN, NE1, H2, H2F, H2T, SI, SIF, SIT, UI, ZH2T, NEF, NET, NI1
      INTEGER :: I, J
      DOUBLE PRECISION :: H2LN_T, H2A, H2BT, ZH2TT, ZH2S, QF, QLN_F, QLN_T, QBF, QBT
      DOUBLE PRECISION :: HI, HIF, HIT, HG, HGT, HGF
      DOUBLE PRECISION :: NX, NXF, NXT, DH2TI, D1TI, D2TI, ZH2, QA, QB, QC, QD, D3TI, ZET, DZH2T, DZH2TT, SHA, SJHA, SCHA
      DOUBLE PRECISION, PARAMETER :: EPSI = 1D-30
      DOUBLE PRECISION :: HA(26), VA(26)
      DOUBLE PRECISION :: FXP, VX
      FXP(VX) = DEXP(DMAX1(-50D0, DMIN1(50D0, VX))) ! statement function for bounded EXP
      ! Contributions of the completely ionized species
      NE1 = 0D0; NEF = 0D0; NET = 0D0; SI = 0D0; SIF = 0D0; SIT = 0D0; UI = 0D0
      DO I = NEL,KION+1,-1 ! elements after KION are assumed completely ionized
         NE1 = NE1 + KZN(I)*NA(I)
         SI = SI - NA(I)*DLOG(NA(I)/CAN1pt5(I) + EPSI) ! add to entropy
      END DO
      ! Calculate ionization of the first KION of H, He, C, N, O, Ne, Mg, Si, Fe.
      DO I = KION,1,-1
         SHA = 1D0; SJHA = 0D0; SCHA = 0D0
         ! compute potentials VA and number ratios HA of ionization state J
         ! relative to the ground state
         DO J = 1, KZN(I) ! for each electron
            VA(J) = -CHI(J,I)*TI + J*DV
            IF ( J.EQ.1 ) THEN
               HA(J) = FXP(VA(J))*COM(KZN(I))/COM(KZN(I)+1)
            ELSE
               HA(J) = HA(J-1)*FXP(VA(J)-VA(J-1))*COM(KZN(I)+1-J)/COM(KZN(I)+2-J)
            END IF
            SHA = SHA + HA(J)
            SJHA = SJHA + J*HA(J)
            SCHA = SCHA + CHI(J,I)*TI*HA(J)
         END DO
         SI = SI + NA(I)*LN_CAN1pt5(I)
         IF (I .EQ. 1) EXIT ! hydrogen gets special treatment
         VX = NA(I)/SHA
         SI = SI - VX*DLOG(VX/COM(KZN(I)+1) + EPSI)
         DO J = 1, KZN(I)
            NX = HA(J)*VX
            NXF = NX*DVF*(J - SJHA/SHA)
            NXT = NXF*DVT/DVF + NX*(CHI(J,I)*TI - SCHA/SHA)
            ! keep record of ionization mass fractions for He, C, N, and O.
            IF (I .EQ. 2) THEN
               XHE_ions(J) = NX*CBN(I)
            END IF
            NE1 = NE1 + J*NX
            NEF = NEF + J*NXF
            NET = NET + J*NXT
            SI = SI - NX*DLOG(NX/COM(KZN(I) + 1 - J) + EPSI)
            SIF = SIF - VA(J)*NXF
            SIT = SIT - VA(J)*NXT
            UI = UI + CHI(J,I)*NX
         END DO
      END DO
      ! Ionization and molecular dissciation of Hydrogen.
      ! partition function for H2 from Vardya (1960), Webbink (1975)
      DH2TI = CH2*TI
      D1TI = C1*TI
      D2TI = (C2*TI)**2
      D3TI = (C3*TI)**3
      ZET = 1D0 - (1D0 + DH2TI)*DEXP(-DH2TI)
      DZH2T = -DH2TI**2*DEXP(-DH2TI)/ZET
      DZH2TT = (DH2TI - 2D0 - DZH2T)*DZH2T
      ZH2 = 6608.8D0*ZET*DH2TI**(-2.5D0)*DEXP(-D1TI - D2TI - D3TI)
      ZH2T = 2.5D0 + D1TI + 2D0*D2TI + 3D0*D3TI + DZH2T
      ZH2TT = -D1TI - 4D0*D2TI - 9D0*D3TI + DZH2TT
      ZH2S = ZH2*DSQRT(8D0)/CAN1pt5(1)
      H2A = CEN*(ZH2S/4D0)*DE/(T*DSQRT(T))*DEXP(DH2TI)
      H2BT = DH2TI + 1.5D0 - ZH2T
      H2LN_T = RET - H2BT
      ! solve for densities of H+, H, and H2
      QA = 2D0*H2A + HA(1)*(1D0 + HA(1))
      QB = NE1 + HA(N_H)*(NE1 - NA(N_H))
      QC = NA(N_H)*NE1
      HG = 2D0*QC/(DSQRT(QB*QB + 4D0*QA*QC) + QB)
      HI = HA(N_H)*HG
      XH_plus = HI*CBN(N_H) ! mass fraction of H+
      NE1 = NE1 + HI
      EN = 1D0/NE1
      MUE = EN  ! amu per free electron
      H2 = H2A*HG*HG*EN ! H2 = moles of molecular hydrogen per gram
      XH2 = H2*CBN(N_H)*2D0 ! mass fraction of molecular hydrogen
      NI1 = NI - H2 ! moles of molecules per gram = (moles of atoms per gram) - (moles of H2 per gram)
      ! derivatives w.r.t. F and T
      QA = NE1 + 4D0*HG*H2A
      QB = HA(N_H)*(NE1 - 2D0*H2)
      QD = 1D0/(QA + QB)
      QC = 2D0*H2*QD
      QLN_F = (NEF - NE1*REF)*QC
      QLN_T = (NET - NE1*H2LN_T)*QC
      QF = HG*QD
      QBF = DVF*QF
      QBT = (CHI(1,1)*TI + DVT)*QF
      HGF = QLN_F - QB*QBF
      HGT = QLN_T - QB*QBT
      HIF = HA(N_H)*(QLN_F + QA*QBF)
      HIT = HA(N_H)*(QLN_T + QA*QBT)
      NEF = NEF + HIF
      NET = NET + HIT
      H2F = H2*REF + EN*(2D0*H2A*HG*HGF - H2*NEF)
      H2T = H2*H2LN_T + EN*(2D0*H2A*HG*HGT - H2*NET)
      ! hydrogen contribution to entropy, internal energy
      SIF = SIF - VA(N_H)*HIF - H2BT*H2F
      SIT = SIT - VA(N_H)*HIT - H2BT*H2T + H2*(ZH2T + ZH2TT)
      SI = SI - HI*DLOG(DABS(HI/COM(1))+EPSI)-HG*DLOG(DABS(HG/COM(2))+EPSI)-H2*(DLOG(DABS(H2/ZH2S)+EPSI)+ZH2T)
      UI = UI + CHI(1,1)*HI + 0.5*CH2*(HI + HG)
      END SUBROUTINE EoS_Ionization
      
      END MODULE ez_ionization

