      MODULE ez_nuclear
      USE ez_nuclear_data
      IMPLICIT NONE
      
      CONTAINS
      
      DOUBLE PRECISION FUNCTION NUCRAT ( TL, ZT, N, NI, NE, AVM, RHO, ENX )
      USE star_controls
      ! Compute rates of nuclear reactions and the corresponding energy and neutrino release.
      DOUBLE PRECISION, INTENT(IN) :: TL, ZT, N(NEL), NI, NE, AVM, RHO
      DOUBLE PRECISION, INTENT(OUT) :: ENX
      DOUBLE PRECISION :: RRT(19)
      DOUBLE PRECISION :: CSA, CSD1, CXD, EX
      DOUBLE PRECISION :: STRN, SCRN, FPNG, DSTR, TT, RR, TU, PP3, PP2, QNPP, QPPT, FCCG, F34, WB
      DOUBLE PRECISION :: XB, WA, CSC, CSB, WC, RHB, TF, ZD, RN, ZA, VL, ZC, ZB, RPNA, VX, RCCA
      INTEGER :: IT, IJ
      DOUBLE PRECISION :: CBRT
      CBRT(VX) = DEXP(C3RD*DLOG(VX)) ! statement function for cube root
      VX = 0 ! for Absoft compiler
      CSA = 0.624D0; CSB = 0.316D0; CSC = 0.460D0; CSD1 = 0.38D0; CXD = 0.86D0
      RHB = RHO/AVM ! RHB is density corrected for non-integer atomic weights
      ! Electron screening theory from Graboske, DeWitt, Grossman & Cooper (1973),
      ! for strong (ZA, ZB, ZC) and intermediate screening (ZD).
      ! The reaction dependent charge parameters are stored in CZA ... CZD.
      WC = 0D0
      DO IJ=NEL,1,-1 ! important to sum in this order so little ones get a chance to add up first
         ! otherwise they get dropped in roundoff and the final result can be quite different for the model.
         WC = WC + VZ(IJ)*N(IJ)
      END DO
      WC = WC/NI
      WB = NE/NI
      WA = ZT*ZT/(WB*WB)
      XB = CBRT(DABS(WB))
      VL = CPL*DSQRT(DABS(NI)*RHB*DEXP(-3D0*TL))
      ZA = CSA*XB*DEXP(DLOG(VL)/1.5D0)
      ZB = CSB*XB*ZA
      ZC = CSC/(XB*XB)
      ZD = CSD1*WC*WA*DEXP(DLOG(VL/(WA*ZT))*CXD)
      ! Reaction rates interpolated in T, mostly from Caughlan & Fowler (1988)
      TF = TL/CLN ! TF is log10(T)
      DO IJ = 1, 19 ! set T dependent factors for rates RPP, R33, R34, ...
         RN = 0D0
         TT = 50D0*(TF - 6D0) + 1D0
         IF ( TT .GE. 1D0 ) THEN
            IT = MAX0(1, MIN0(199, INT(TT)))
            TT = TT - IT
            TU = 1D0 - TT
            RR = TU*CRT(IT, IJ) + TT*CRT(IT + 1, IJ)
            IF ( RR .GE. -50D0 ) THEN
               SCRN = ZD*CZD(IJ)
               STRN = ZA*CZA(IJ) + ZB*CZB(IJ)
               DSTR = ZC*CZC(IJ)
               IF (DSTR .LT. 0.29D0*STRN) SCRN = MIN(SCRN, STRN - DSTR)
               RN = DEXP(CLN*(RR + 20D0) + SCRN)*1.0D-20
            END IF
         END IF
         RRT(IJ) = RN ! RRT(1) is for RPP, RRT(2) is for R33, etc.
      END DO
      ! note that abundances of He3 and Be7 are not needed in equilibrium
      RPP = RHB*N(N_H)*N(N_H)*RRT(1)/2D0
      R33 = RHB*RRT(2)/2D0
      R34 = RHB*N(N_HE)*RRT(3)
      RBE = RHB*NE*RRT(4)
      RBP = RHB*N(N_H)*RRT(5)
      RPC = RHB*N(N_H)*N(N_C)*RRT(6)
      RPN = RHB*N(N_H)*N(N_N)*RRT(7)
      RPO = RHB*N(N_H)*N(N_O)*RRT(8)
      R3A = RHB*RHB*N(N_HE)*N(N_HE)*N(N_HE)*RRT(9)/6D0
      RAC = RHB*N(N_HE)*N(N_C)*RRT(10)
      RAN = RHB*N(N_HE)*N(N_N)*RRT(11)
      RAO = RHB*N(N_HE)*N(N_O)*RRT(12)
      RANE = RHB*N(N_HE)*N(N_NE)*RRT(13)
      RCC = RHB*N(N_C)*N(N_C)*RRT(14)/2D0
      RCO = RHB*N(N_C)*N(N_O)*RRT(15)
      ROO = RHB*N(N_O)*N(N_O)*RRT(16)/2D0
      RGNE = N(N_NE)*RRT(17)
      RGMG = N(N_MG)*RRT(18)
      
      
      IF (CNO_Factor /= 1d0) THEN
         RPC = CNO_Factor * RPC; RPO = CNO_Factor * RPO; RPN = CNO_Factor * RPN; RPNG = CNO_Factor * RPNG
      END IF
      
      
      ! Branching of pN and CC reactions
      FPNG = 8.0D-4
      RPNA = (1D0 - FPNG)*RPN
      RPNG = FPNG*RPN
      RPN = RPNA
      FCCG = RRT(19)
      RCCA = (1D0 - FCCG)*RCC
      RCCG = FCCG*RCC
      RCC = RCCA
      ! PP chain in equilibrium, RPP becomes effective rate of 2 H1 -> 0.5 He4
      F34 = 0D0
      IF (R34 .GT. 1.0D-20) F34 = 2D0/(1D0 + DSQRT(1D0 + 8D0*RPP*R33/(R34*R34)))
      RPP = RPP*(1D0 + F34)
      PP2 = 1D0
      IF (RBE + RBP .GT. 1.0D-20) PP2 = RBE/(RBE + RBP)
      PP3 = 1D0 - PP2
      QPPT = QPP + 0.5D0*Q33
      QNPP = (QNT(1) + F34*(QNT(4)*PP2 + QNT(5)*PP3))/(1D0 + F34)
      ! calculate total nuclear energy release and neutrino loss
      EX = QPPT*RPP + QPC*RPC + QPNA*RPN + QPO*RPO + Q3A*R3A + QAC*RAC
      EX = EX + QAN*RAN + QAO*RAO + QANE*RANE + QCCA*RCC + QCO*RCO
      EX = EX + QOO*ROO + QGNE*RGNE + QGMG*RGMG + QCCG*RCCG + QPNG*RPNG
      EX = CME*EX/AVM
      ENX = QNPP*RPP + QNT(6)*RPC + QNT(7)*RPN + QNT(8)*RPO + QNT(9)*R3A + QNT(10)*RAC
      ENX = ENX + QNT(11)*RAN + QNT(12)*RAO + QNT(13)*RANE + QNT(14)*RCC + QNT(15)*RCO
      ENX = ENX + QNT(16)*ROO + QNT(17)*RGNE + QNT(18)*RGMG + QNT(19)*RCCG + QNT(20)*RPNG
      ENX =  -CME*ENX/AVM
      ! EX is nuclear energy generation rate (in ergs/gram/sec)
      ! ENX is neutrino energy rate from nuclear reactions (in ergs/gram/sec)
      NUCRAT = EX
      END FUNCTION NUCRAT

      END MODULE ez_nuclear

