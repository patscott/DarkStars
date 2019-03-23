      MODULE ez_vcool_data
      IMPLICIT NONE
      
      DOUBLE PRECISION, PARAMETER :: nnu=2.0 ! number of neutrino species
      DOUBLE PRECISION, PARAMETER :: sin2thetaw=0.2319 ! weak interaction parameter
      DOUBLE PRECISION ::  CV, CAX, CVprime, CAXprime
      DOUBLE PRECISION :: LN_10, PI
      DOUBLE PRECISION ::  CV_SAV, CAX_SAV, CVprime_SAV, CAXprime_SAV
      DOUBLE PRECISION :: LN_10_SAV, PI_SAV
      
      CONTAINS
      
      SUBROUTINE SAV_vcool_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) CV, CAX, CVprime, CAXprime, LN_10, PI
      ELSE
         WRITE (IO_UNIT) CV, CAX, CVprime, CAXprime, LN_10, PI
      END IF
      END SUBROUTINE SAV_vcool_data

      SUBROUTINE SAV_vcool_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         CV=CV_SAV; CAX=CAX_SAV; CVprime=CVprime_SAV; CAXprime=CAXprime_SAV; LN_10=LN_10_SAV; PI=PI_SAV
      ELSE
         CV_SAV=CV; CAX_SAV=CAX; CVprime_SAV=CVprime; CAXprime_SAV=CAXprime; LN_10_SAV=LN_10; PI_SAV=PI
      END IF
      END SUBROUTINE SAV_vcool_data_internal

      END MODULE ez_vcool_data

