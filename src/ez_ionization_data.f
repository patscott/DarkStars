      MODULE ez_ionization_data    ! equation of state data module
      USE star_constants
      IMPLICIT NONE

      DOUBLE PRECISION :: MUE
      
      DOUBLE PRECISION :: XH2 ! mass fraction for molecular hydrogen
      DOUBLE PRECISION :: XH_plus
      ! arrays of mass fractions for the various ionization levels.  I'th array entry is for I electrons removed.
      DOUBLE PRECISION :: XHE_ions(2)

      DOUBLE PRECISION :: CH2,   C1, C2, C3  ! CH2, C1, C2, C3 are molecular hydrogen parameters.
      DOUBLE PRECISION :: CH2_SAV,   C1_SAV, C2_SAV, C3_SAV
	  DOUBLE PRECISION :: CHI(26,NEL), CHI_SAV(26,NEL)			 ! CHI(J,I) is energy (in eV) for Jth ionization level of element I.
      DOUBLE PRECISION :: COM(27), COM_SAV(27)                   ! statistical weights (for Saha equation).
      DOUBLE PRECISION :: CAN(NEL), CAN_SAV(NEL)                 ! baryons per nucleus of element I; using average atomic weights for elements
      DOUBLE PRECISION :: CAN1pt5(NEL), CAN1pt5_SAV(NEL)         ! CAN^(3/2)
      DOUBLE PRECISION :: LN_CAN1pt5(NEL), LN_CAN1pt5_SAV(NEL)   ! LOG(CAN^(3/2))
      DOUBLE PRECISION :: CBN(NEL), CBN_SAV(NEL)                 ! baryons per nucleus of element I; using integer atomic weights for elements
      INTEGER :: KZN(NEL), KZN_SAV(NEL)           ! atomic number of element I.

      !      KZN   CBN     CAN
      !     H    1   1.000   1.008
      !     He   2   4.000   4.003
      !     C    6  12.000  12.000
      !     N    7  14.000  14.003
      !     O    8  16.000  15.995
      !     Ne  10  20.000  19.992
      !     Mg  12  24.000  23.985
      !     Si  14  28.000  27.977
      !     Fe  26  56.000  55.935

      CONTAINS
      
      SUBROUTINE SAV_ionization_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT)  CH2, C1, C2, C3, CHI, COM, CAN, CAN1pt5, LN_CAN1pt5, CBN, KZN
      ELSE
         WRITE (IO_UNIT) CH2, C1, C2, C3, CHI, COM, CAN, CAN1pt5, LN_CAN1pt5, CBN, KZN
      END IF
      END SUBROUTINE SAV_ionization_data

      SUBROUTINE SAV_ionization_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         CH2=CH2_SAV;C1=C1_SAV;C2=C2_SAV;C3=C3_SAV;CHI=CHI_SAV;COM=COM_SAV
		 CAN=CAN_SAV;CAN1pt5=CAN1pt5_SAV;LN_CAN1pt5=LN_CAN1pt5_SAV;CBN=CBN_SAV;KZN=KZN_SAV
      ELSE
         CH2_SAV=CH2;C1_SAV=C1;C2_SAV=C2;C3_SAV=C3;CHI_SAV=CHI;COM_SAV=COM
		 CAN_SAV=CAN;CAN1pt5_SAV=CAN1pt5;LN_CAN1pt5_SAV=LN_CAN1pt5;CBN_SAV=CBN;KZN_SAV=KZN
      END IF
      END SUBROUTINE SAV_ionization_data_internal

      END MODULE ez_ionization_data
