      MODULE ez_nuclear_data    ! nuclear reaction data
      USE star_constants
      IMPLICIT NONE

      DOUBLE PRECISION ::  RPP  ! equilibrium pp chain    p (p,positron+neutrino) d (p,gamma) He3
      DOUBLE PRECISION ::  R33  ! He3 + He3 (ppI chain)   He3 (He3,2p) He4
      DOUBLE PRECISION ::  R34  ! He3 + He4               He3 (He4,gamma) Be7
      DOUBLE PRECISION ::  RBP  ! Be7 + p  (ppII chain)   Be7 (p,gamma) B8 (positron+neutrino) Be8* (alpha) He4
      DOUBLE PRECISION ::  RBE  ! Be7 decay (ppIII chain) Be7 (electron,neutrino) Li7 (p,alpha) He4
      
      DOUBLE PRECISION ::  RPC  ! C12 + p (CN cycle)      C12 (p,positron+neutrino) C13 (p,gamma) N14
      DOUBLE PRECISION ::  RPO  ! O16 + p (ON cycle)      O16 (p,positron+neutrino) O17 (p,alpha) N14
      DOUBLE PRECISION ::  RPNG ! less probable N15+p     N14 (p,positron+neutrino) N15 (p,gamma) O16
      DOUBLE PRECISION ::  RPN  ! more probable N15+p     N14 (p,positron+neutrino) N15 (p,alpha) C12      

      DOUBLE PRECISION ::  R3A  ! 3 alpha                 He4 (alpha) Be8* (alpha,gamma) C12
      DOUBLE PRECISION ::  RAC  ! alpha + C12             C12 (alpha,gamma) O16
      DOUBLE PRECISION ::  RAN  ! 3/2 alpha + N14         N14 (alpha,gamma) F18 (1/2 alpha,gamma) Ne20
      DOUBLE PRECISION ::  RAO  ! alpha + O16             O16 (alpha,gamma) Ne20
      DOUBLE PRECISION ::  RANE ! alpha + Ne20            Ne20 (alpha,gamma) Mg24
      DOUBLE PRECISION ::  RCO  ! C12 + O16               C12 (O16,alpha+gamma) Mg24
      DOUBLE PRECISION ::  ROO  ! 2 O16                   O16 (O16,alpha+gamma) Si28 (gamma,alpha) Mg24
      DOUBLE PRECISION ::  RGNE ! Ne20 decay              Ne20 (gamma,alpha) O16
      DOUBLE PRECISION ::  RGMG ! Mg24 decay              Mg24 (gamma,alpha) Ne20
      DOUBLE PRECISION ::  RCCG ! less probable C+C       C12 (C12,gamma) Mg24
      DOUBLE PRECISION ::  RCC  ! more probable C+C       C12 (C12,alpha+gamma) Ne20

      DOUBLE PRECISION ::  CNU(60,41,2), CNU_SAV(60,41,2) ! table of nuclear reaction rates.
      DOUBLE PRECISION ::  CRT(200,20), CRT_SAV(200,20)  ! data used in computing the T dependent factors for nuclear reaction rates.

      DOUBLE PRECISION ::  VZ(NEL), VZ_SAV(NEL)              

      ! nuclear reaction data.  The Qs are MeVs per reaction.
      DOUBLE PRECISION ::  QPP, QPP_SAV  ! p (p,positron+neutrino) d (p,gamma) He3
      DOUBLE PRECISION ::  Q33, Q33_SAV  ! He3 (He3,2p) He4.
      DOUBLE PRECISION ::  Q34, Q34_SAV  ! He3 (He4,gamma) Be7
      DOUBLE PRECISION ::  QBE, QBE_SAV  ! 7be(e-,nu)7li(p,gamma)2 He4
      DOUBLE PRECISION ::  QBP, QBP_SAV  ! 7be(p,g)8b(,positron+neutrino)be8(,gamma)2 He4
      DOUBLE PRECISION ::  QPC, QPC_SAV  ! C12 + 2 p -> N14.
      DOUBLE PRECISION ::  QPNA, QPNA_SAV! N14 + 2 p -> C12 + He4.
      DOUBLE PRECISION ::  QPO, QPO_SAV  ! O16 + 2 p -> N14 + He4.
      DOUBLE PRECISION ::  Q3A, Q3A_SAV  ! triple alpha.
      DOUBLE PRECISION ::  QAC, QAC_SAV  ! C12 + He4 -> O16.
      DOUBLE PRECISION ::  QAN, QAN_SAV  ! N14 + 3/2 He4 -> Ne20.
      DOUBLE PRECISION ::  QAO, QAO_SAV  ! O16 + He4 -> Ne20.
      DOUBLE PRECISION ::  QANE, QANE_SAV! Ne20 + He4 -> Mg24.
      DOUBLE PRECISION ::  QCCA, QCCA_SAV! 2 C12 -> Ne20 + He4.
      DOUBLE PRECISION ::  QCO, QCO_SAV  ! C12 + O16 -> Mg24 + He4.
      DOUBLE PRECISION ::  QOO, QOO_SAV  ! 2 O16 -> Mg24 + 2 He4.
      DOUBLE PRECISION ::  QGNE, QGNE_SAV! Ne20 -> O16 + He4.
      DOUBLE PRECISION ::  QGMG, QGMG_SAV! Mg24 -> Ne20 + He4.
      DOUBLE PRECISION ::  QCCG, QCCG_SAV! 2 C12 -> Mg24.
      DOUBLE PRECISION ::  QPNG, QPNG_SAV! N14 + 2 p -> O16.
      
      INTEGER, PARAMETER :: NUMR=20 ! number of reactions in the tables below.
      ! nuclear reaction data.  The Qs are MeVs per reaction.
      DOUBLE PRECISION ::  QNT(NUMR), QNT_SAV(NUMR)   ! neutrino Q values in MeV per reaction.
      DOUBLE PRECISION ::  CZA(NUMR), CZB(NUMR), CZC(NUMR), CZD(NUMR) ! related to electron screening
      DOUBLE PRECISION ::  CZA_SAV(NUMR), CZB_SAV(NUMR), CZC_SAV(NUMR), CZD_SAV(NUMR)
	        
      CONTAINS
      
      SUBROUTINE SAV_nuclear_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) CNU, CRT, VZ
         READ (IO_UNIT) QPP, Q33, Q34, QBE, QBP, QPC, QPNA, QPO, Q3A, QAC
         READ (IO_UNIT) QAN, QAO, QANE, QCCA, QCO, QOO, QGNE, QGMG, QCCG, QPNG
         READ (IO_UNIT) QNT, CZA, CZB, CZC, CZD
      ELSE
         WRITE (IO_UNIT) CNU, CRT, VZ
         WRITE (IO_UNIT) QPP, Q33, Q34, QBE, QBP, QPC, QPNA, QPO, Q3A, QAC
         WRITE (IO_UNIT) QAN, QAO, QANE, QCCA, QCO, QOO, QGNE, QGMG, QCCG, QPNG
         WRITE (IO_UNIT) QNT, CZA, CZB, CZC, CZD
      END IF
      END SUBROUTINE SAV_nuclear_data

      SUBROUTINE SAV_nuclear_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         CNU=CNU_SAV; CRT=CRT_SAV; VZ=VZ_SAV
         QPP=QPP_SAV; Q33=Q33_SAV; Q34=Q34_SAV; QBE=QBE_SAV; QBP=QBP_SAV 
		 QPC=QPC_SAV; QPNA=QPNA_SAV; QPO=QPO_SAV; Q3A=Q3A_SAV; QAC=QAC_SAV
         QAN=QAN_SAV; QAO=QAO_SAV; QANE=QANE_SAV; QCCA=QCCA_SAV; QCO=QCO_SAV
		 QOO=QOO_SAV; QGNE=QGNE_SAV; QGMG=QGMG_SAV; QCCG=QCCG_SAV; QPNG=QPNG_SAV
         QNT=QNT_SAV; CZA=CZA_SAV; CZB=CZB_SAV; CZC=CZC_SAV; CZD=CZD_SAV
      ELSE
         CNU_SAV=CNU; CRT_SAV=CRT; VZ_SAV=VZ
         QPP_SAV=QPP; Q33_SAV=Q33; Q34_SAV=Q34; QBE_SAV=QBE; QBP_SAV=QBP 
		 QPC_SAV=QPC; QPNA_SAV=QPNA; QPO_SAV=QPO; Q3A_SAV=Q3A; QAC_SAV=QAC
         QAN_SAV=QAN; QAO_SAV=QAO; QANE_SAV=QANE; QCCA_SAV=QCCA; QCO_SAV=QCO
		 QOO_SAV=QOO; QGNE_SAV=QGNE; QGMG_SAV=QGMG; QCCG_SAV=QCCG; QPNG_SAV=QPNG
         QNT_SAV=QNT; CZA_SAV=CZA; CZB_SAV=CZB; CZC_SAV=CZC; CZD_SAV=CZD
      END IF
      END SUBROUTINE SAV_nuclear_data_internal

      END MODULE ez_nuclear_data
