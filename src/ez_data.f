      MODULE ez_data    ! the main data module for the evolution code
						! internal save routines and variables added by Pat Scott 24-03-08
      USE star_controls
      USE star_constants
      USE star_data
      USE ez_cycle_data
      IMPLICIT NONE
      
      DOUBLE PRECISION :: init_M, init_M_SAV ! initial total stellar mass (in Msolar units)
      DOUBLE PRECISION :: init_Z, init_Z_SAV ! initial metallicity
      
      ! for use in gathering statistics for debugging
      INTEGER :: stats_N
      DOUBLE PRECISION :: stats_X, stats_X2

	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      DOUBLE PRECISION :: H(V_MAX,MAX_N_SHELLs), H_SAV(V_MAX,MAX_N_SHELLs)   
												! H(J,K) is matrix of independent variables.
												!  J=1,V_MAX is the variable number.
                                                !  K=1,N_SHELLs is the shell number.
                                                !  K=1 for the surface; K=N_SHELLs for the center.
      DOUBLE PRECISION :: DH(V_MAX,MAX_N_SHELLs), DH_SAV(V_MAX,MAX_N_SHELLs)
												! DH(J,K) is increment for value of variable J at mesh point K.
      DOUBLE PRECISION :: XTRAS(N_XTRAS), XTRAS_SAV(N_XTRAS)
												! These will automatically be backed up and restored along with state in H
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      ! I/O unit numbers
      INTEGER, PARAMETER :: IO_EOS=20, IO_NUC=21, IO_CONFIG=22, IO_SAVE_RESTORE=23, IO_ZAMS=24, IO_TMP=25, IO_UBV=26
      
      CHARACTER (LEN=strlen) :: EZ_DATA_DIR ! where to look for data on ZAMS and metallicities
          
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER, PARAMETER :: SAVE_RESTORE_ID = 121 ! change this whenever modify form of save files

      CONTAINS 
      
      SUBROUTINE SysOut ( fname, ios )
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      INTEGER, INTENT(OUT) :: ios
      CHARACTER (LEN=strlen) :: name
      name = trim(fname); ios = -1
      OPEN (UNIT=IO_SAVE_RESTORE, FILE=name, ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios, FORM='UNFORMATTED')
      IF (ios .NE. 0) RETURN
      WRITE ( IO_SAVE_RESTORE ) SAVE_RESTORE_ID
      CALL Model_IO ( IO_SAVE_RESTORE, .FALSE. )
      END SUBROUTINE SysOut
      
      SUBROUTINE SysIn ( fname, ios )
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      INTEGER, INTENT(OUT) :: ios
      CHARACTER (LEN=strlen) :: name
      INTEGER ID
      name = trim(fname); ios = 0
      OPEN (UNIT=IO_SAVE_RESTORE, FILE=name, ACTION='READ', IOSTAT=ios, FORM='UNFORMATTED')
      IF (ios .NE. 0) RETURN
      READ ( IO_SAVE_RESTORE ) ID
      IF ( ID .NE. SAVE_RESTORE_ID ) THEN
         CLOSE ( SAVE_RESTORE_ID )
         ios = -10
         WRITE (*,*) name, ' is not a valid file for this version of the system.'
         RETURN
      END IF
      CALL Model_IO ( IO_SAVE_RESTORE, .TRUE. )
      END SUBROUTINE SysIn
      
      SUBROUTINE Model_IO  ( IO_UNIT, read_flag )
      USE star_constants
      USE star_controls
      USE ez_opacity_data
      USE ez_vcool_data
      USE ez_nuclear_data
      USE ez_state_data
      USE ez_solve_data
      USE ez_ionization_data
      USE ez_shell_data
      USE ez_cycle_data
      USE ez_do_one_data
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      CALL SAV_star_constants( IO_UNIT, read_flag )
      CALL SAV_star_controls( IO_UNIT, read_flag )
      CALL SAV_opacity_data( IO_UNIT, read_flag )
      CALL SAV_vcool_data( IO_UNIT, read_flag )
      CALL SAV_state_data( IO_UNIT, read_flag )
      CALL SAV_ionization_data( IO_UNIT, read_flag )
      CALL SAV_nuclear_data( IO_UNIT, read_flag )
      CALL SAV_solve_data( IO_UNIT, read_flag )
      CALL SAV_cycle_data( IO_UNIT, read_flag )
      CALL SAV_do_one_data( IO_UNIT, read_flag )
      CALL SAV_ez_data( IO_UNIT, read_flag )
      CLOSE ( IO_UNIT )
      END SUBROUTINE Model_IO

      SUBROUTINE Model_IO_internal (read_flag )
      USE star_constants
      USE star_controls
      USE ez_opacity_data
      USE ez_vcool_data
      USE ez_nuclear_data
      USE ez_state_data
      USE ez_solve_data
      USE ez_ionization_data
      USE ez_shell_data
      USE ez_cycle_data
      USE ez_do_one_data
      LOGICAL, INTENT(IN) :: read_flag
      CALL SAV_star_constants_internal(read_flag )
      CALL SAV_star_controls_internal(read_flag )
      CALL SAV_opacity_data_internal(read_flag )
      CALL SAV_vcool_data_internal(read_flag )
      CALL SAV_state_data_internal(read_flag )
      CALL SAV_ionization_data_internal(read_flag )
      CALL SAV_nuclear_data_internal(read_flag )
      CALL SAV_solve_data_internal(read_flag )
      CALL SAV_cycle_data_internal(read_flag )
      CALL SAV_do_one_data_internal(read_flag )
      CALL SAV_ez_data_internal(read_flag )
      END SUBROUTINE Model_IO_internal


      SUBROUTINE SAV_ez_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) H, DH, XTRAS, init_M, init_Z
      ELSE
         WRITE (IO_UNIT) H, DH, XTRAS, init_M, init_Z
      END IF
      END SUBROUTINE SAV_ez_data

      SUBROUTINE SAV_ez_data_internal(read_flag )
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         H=H_SAV; DH=DH_SAV; XTRAS=XTRAS_SAV; init_M=init_M_SAV; init_Z=init_Z_SAV
      ELSE
         H_SAV=H; DH_SAV=DH; XTRAS_SAV=XTRAS; init_M_SAV=init_M; init_Z_SAV=init_Z
      END IF
      END SUBROUTINE SAV_ez_data_internal

      END MODULE ez_data
