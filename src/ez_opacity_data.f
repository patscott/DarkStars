      MODULE ez_opacity_data
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: MT = 127, MR = 90 ! MT is size of temperature array; MR is size of density array
      DOUBLE PRECISION ::  FSPL(4,4,MT,MR,10), FSPL_SAV(4,4,MT,MR,10)   ! opacity interpolation tables.
      DOUBLE PRECISION ::  TFM(MT), TFM_SAV(MT)              ! temperatures for opacity interpolation.
      DOUBLE PRECISION ::  FRM(MR), FRM_SAV(MR)              ! densities for opacity interpolation.  
      DOUBLE PRECISION ::  FKLM(6), FKLM_SAV(6)
      DOUBLE PRECISION ::  FKHM(6), FKHM_SAV(6)
      DOUBLE PRECISION :: CSX(10), CSX_SAV(10)               ! parameters for subcompositions for opacities.
      DOUBLE PRECISION :: CS(MR,MT,10), CS_SAV(MR,MT,10)     ! the opacity tables.
      INTEGER :: KCSX, KCSX_SAV                     ! number of subcompositions tabulated for opacity tables.

      CONTAINS
      
      SUBROUTINE SAV_opacity_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) FSPL, TFM, FRM, FKLM, FKHM, CSX, CS, KCSX
      ELSE
         WRITE (IO_UNIT) FSPL, TFM, FRM, FKLM, FKHM, CSX, CS, KCSX
      END IF
      END SUBROUTINE SAV_opacity_data

      SUBROUTINE SAV_opacity_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         FSPL=FSPL_SAV; TFM=TFM_SAV; FRM=FRM_SAV; FKLM=FKLM_SAV; FKHM=FKHM_SAV
		 CSX=CSX_SAV; CS=CS_SAV; KCSX=KCSX_SAV
      ELSE
         FSPL_SAV=FSPL; TFM_SAV=TFM; FRM_SAV=FRM; FKLM_SAV=FKLM; FKHM_SAV=FKHM
		 CSX_SAV=CSX; CS_SAV=CS; KCSX_SAV=KCSX
      END IF
      END SUBROUTINE SAV_opacity_data_internal

      END MODULE ez_opacity_data

