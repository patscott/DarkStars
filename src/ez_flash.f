      MODULE ez_flash
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
      USE ez_setup
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: io_flash = 97
      
      LOGICAL, PARAMETER :: building_He = .FALSE.
      ! to build EZ_HE_FLASH.data file, set this flag to .TRUE.
      ! change dispatch to call Post_He_Flash
      ! set metals_dir in Post_He_Flash for the case you are doing
      ! the data file is created in the current metals subdirectory.
      ! need to run this for each of the Z's
      ! pick the appropriate setting for age_He_Settle and for min_core_Mass
      ! it should go until it is ready to start the mass loss at which point it writes the file
      ! >>> make sure it has settled into steady helium burning by the time it quits
      ! then you should set the building_He flag to .FALSE. and run again to check it

      ! age_He_Settle is only used during the creation of the EZ_HE_FLASH.data files
      !DOUBLE PRECISION, PARAMETER :: age_He_Settle = 6.34D8 ! for Z=0.03
      !DOUBLE PRECISION, PARAMETER :: age_He_Settle = 6.1D8 ! for Z=0.02
      !DOUBLE PRECISION, PARAMETER :: age_He_Settle = 5.6D8 ! for Z=0.01
      !DOUBLE PRECISION, PARAMETER :: age_He_Settle = 5.1D8 ! for Z=0.004
      DOUBLE PRECISION, PARAMETER :: age_He_Settle = 5.0D8 ! for Z=0.001
      !DOUBLE PRECISION, PARAMETER :: age_He_Settle = 5.0D8 ! for Z=0.0003
      !DOUBLE PRECISION, PARAMETER :: age_He_Settle = 4.8D8 ! for Z=0.0001

      ! min_core_Mass is only used during the building of the EZ_HE_FLASH.data file
      !DOUBLE PRECISION, PARAMETER :: min_core_Mass = 0.37D0 ! for Z=0.03
      !DOUBLE PRECISION, PARAMETER :: min_core_Mass = 0.37D0 ! for Z=0.02
      !DOUBLE PRECISION, PARAMETER :: min_core_Mass = 0.37D0 ! for Z=0.01
      !DOUBLE PRECISION, PARAMETER :: min_core_Mass = 0.38D0 ! for Z=0.004
      DOUBLE PRECISION, PARAMETER :: min_core_Mass = 0.40D0 ! for Z=0.001
      !DOUBLE PRECISION, PARAMETER :: min_core_Mass = 0.408D0 ! for Z=0.0003
      !DOUBLE PRECISION, PARAMETER :: min_core_Mass = 0.4D0 ! for Z=0.0001
      
      DOUBLE PRECISION :: He_total_Mass, He_Settle_Age
      DOUBLE PRECISION :: Post_He_Flash_core_Mass, prev_He_core
      INTEGER :: Post_He_Flash_phase, Begin_mass_loss_model
      INTEGER, PARAMETER :: He_Pre_Flash=1, He_Flash_Settle=3, He_Begin_Mass_Loss=4
      INTEGER, PARAMETER :: He_Do_Mass_Loss=5, He_End_Mass_Loss=6, He_Core_Growth1=7, He_Core_Growth2=8
      
      CONTAINS
      
      SUBROUTINE Post_He_Initial_Params
      END SUBROUTINE Post_He_Initial_Params
      
      INTEGER FUNCTION Post_He_Flash_Check_Model()
      USE star_data
      USE star_extras
      USE ez_shell_data
      USE ez_cycle_data, ONLY : REQUESTED_DTY
      USE ez_report, ONLY : Get_Extra_Model_Data
      DOUBLE PRECISION, PARAMETER :: mass_loss_default = -0.5D-6, start_down = 1.1D0, mass_ramp_size = 10D0, core_watch = 1D-4
      DOUBLE PRECISION, PARAMETER :: log_He_temp = 7.8D0, power_He_Min = 5D0, T_DS_Dt_min = 1D-10
      DOUBLE PRECISION :: amount_to_do, amount_left, frac, mass_loss_param
      IF ( Post_He_Flash_phase .NE. He_Core_Growth1 .AND. Post_He_Flash_phase .NE. He_Core_Growth2 ) THEN
         IF ( mass_He_Core .GT. Post_He_Flash_core_Mass - core_watch ) THEN
            IF (building_He .AND. burn_H .NE. 0D0) WRITE(*,*) ' >>> Stop adding helium for now at core mass of', mass_He_Core
            burn_H = 0D0 ! stop building core for now
         END IF
      END IF
      burn_He = 0D0
      mass_loss_param = mass_loss_default
      if ( star_Mass > 3d0 ) mass_loss_param = mass_loss_param * 0.5d0
      SELECT CASE ( Post_He_Flash_phase )
      CASE ( He_Pre_Flash )
         Post_He_Flash_phase = He_Flash_Settle
      CASE ( He_Flash_Settle )
         IF ( star_Age .GT. He_Settle_Age ) THEN ! we've reached steady helium burning
            Post_He_Flash_phase = He_Core_Growth1
         END IF
      CASE ( He_Core_Growth1 )
         IF (building_He) THEN
            Post_He_Flash_Check_Model = TERMINATE
            RETURN
         END IF
         burn_H = 1D0
         IF ( mass_He_Core .GT. Post_He_Flash_core_Mass - core_watch ) THEN
            IF (building_He) WRITE(*,*) ' >>> Stop adding helium for now at core mass of', mass_He_Core
            Post_He_Flash_phase = He_Begin_Mass_Loss
            Begin_mass_loss_model = model_Number
         END IF
      CASE ( He_Begin_Mass_Loss )
         burn_H = 0D0
         IF ( extra_Mdot_param .GT. mass_loss_param ) THEN ! extra_Mdot_param starts at zero and decreases from there.
            frac = (model_Number - Begin_mass_loss_model) / mass_ramp_size
            extra_Mdot_param = mass_loss_param * frac
            T_DS_Dt = min(1D0, max(T_DS_Dt_min, frac))
         ELSE
            T_DS_Dt = 0D0
            Post_He_Flash_phase = He_Do_Mass_Loss
         END IF
      CASE ( He_Do_Mass_Loss )
         IF ( star_Mass .LE. He_total_Mass*start_down ) THEN
            Begin_mass_loss_model = model_Number
            Post_He_Flash_phase = He_End_Mass_Loss
         END IF
      CASE ( He_End_Mass_Loss )
         IF ( star_Mass .GT. He_total_Mass ) THEN
            amount_to_do = (start_down - 1D0)*He_total_Mass
            amount_left = star_Mass - He_total_Mass
            extra_Mdot_param = mass_loss_param*min(1D0, max(1D-4, amount_left / amount_to_do))
            IF ( extra_Mdot_param .NE. 0D0 ) then
               REQUESTED_DTY = max(1D-2, -amount_left / extra_Mdot_param)
            end if
         ELSEIF ( star_Mass < He_total_Mass * 0.99d0 ) THEN
            Post_He_Flash_Check_Model = BACKUP
            RETURN
         ELSE
            extra_Mdot_param = 0D0
            He_Settle_Age = star_Age * 1.01d0
            Post_He_Flash_phase = He_Core_Growth2
         END IF
      CASE ( He_Core_Growth2 )
         IF ( star_Mass .GT. He_total_Mass*1.001D0 ) THEN ! this can happen if the system backs up over the last mass loss
            Post_He_Flash_phase = He_End_Mass_Loss
         ELSEIF ( Post_He_Flash_core_Mass .GT. mass_He_Core*1.005D0 ) THEN
            amount_to_do = (Post_He_Flash_core_Mass - mass_He_Core) / Post_He_Flash_core_Mass
            burn_H = min(1D0, max(1D-1, amount_to_do))
            IF (mass_He_Core .GT. prev_He_core) THEN
               REQUESTED_DTY = max(10D0, (Post_He_Flash_core_Mass - mass_He_Core) * time_Step / (mass_He_Core - prev_He_core))
            END IF
         ELSEIF ( star_Age > He_Settle_Age ) THEN
            burn_H = 1D0
            burn_He = 1D0
            burn_metals = 1D0
            T_DS_Dt = 1D0
            Post_He_Flash_Check_Model = TERMINATE
            RETURN
         END IF
      CASE DEFAULT
         STOP 'Error in Post_He_Flash_Check_Model'
      END SELECT
      prev_He_core = mass_He_Core
      Post_He_Flash_Check_Model = KEEP_GOING
      END FUNCTION Post_He_Flash_Check_Model

      LOGICAL FUNCTION Get_Post_He_Flash ( total_Mass, core_Mass, age_Offset, Check_Model )
      USE star_data
      USE star_controls
      USE ez_report
      USE ez_cycle_data
      USE ez_cycle
      DOUBLE PRECISION, INTENT(IN) :: total_Mass, core_Mass, age_Offset
      INTERFACE
         INTEGER FUNCTION Check_Model()
         END FUNCTION Check_Model
      END INTERFACE
      INTEGER :: JO
      LOGICAL :: ok_to_evolve
      DOUBLE PRECISION :: min_start_model_Mass, start_model_Mass, mass, he_core, new_age, fac
      INTEGER :: model
      DOUBLE PRECISION :: save_CTRLs(num_CTRLs)
      CHARACTER (LEN=strlen) :: fname
      CALL Save_Control_Values(save_CTRLs)
      
      min_start_model_Mass = 2.5D0
      fac = 2.5d0
      if (core_mass < 12d0) fac = 4d0
      if (core_Mass > 0.48d0) min_start_model_Mass = 3.0d0
      start_model_Mass = MIN(100d0, MAX(total_Mass * 1.25d0, MAX(core_Mass * fac, min_start_model_Mass)))
      if (start_model_Mass > min_start_model_Mass .and. start_model_Mass < 10d0) start_model_Mass = 10d0
      
      mass = total_Mass; he_core = core_Mass
      fname = 'EZ_HE_FLASH.data'
      He_total_Mass = mass
      ok_to_evolve = .FALSE.
      new_age = AGE
      model = JMOD
      IF (building_He) THEN
         IF (.NOT. Read_ZAMS(start_model_Mass, Post_He_Initial_Params)) THEN
            WRITE(*,*)  'Get_Post_He_Flash failed to create ZAMS starting model'
            Get_Post_He_Flash = .FALSE.
            RETURN
         END IF
         ok_to_evolve = .TRUE.
         Post_He_Flash_phase = He_Pre_Flash
      ELSE IF (start_model_Mass > min_start_model_Mass) THEN
         IF ( Read_ZAMS( start_model_Mass, Post_He_Initial_Params ) ) THEN
            IF ( total_Mass == core_Mass ) THEN
               WRITE(*,*)
               IF ( .not. Trade_X_for_Y ( 0.75d0, Check_Model ) ) THEN
                  WRITE(*,*)  'Get_Post_He_Flash failed to create reduced hydrogen starting model'
                  Get_Post_He_Flash = .FALSE.
                  RETURN
               END IF
            END IF
            ok_to_evolve = .TRUE.
            Post_He_Flash_phase = He_Core_Growth1
         END IF
      ELSEIF ( Start_Read_Model(start_model_Mass) ) THEN
         IF ( Read_Model(fname) ) THEN
            IF ( Finish_Read_Model( start_model_Mass, Post_He_Initial_Params, .TRUE. ) ) THEN
               ok_to_evolve = .TRUE.
               Post_He_Flash_phase = He_Core_Growth1
            END IF
         END IF
      END IF
      IF (.NOT. ok_to_evolve) THEN
         Get_Post_He_Flash = .FALSE.
      ELSE
         JO = JO_AOK
         IF (building_He) THEN
            Post_He_Flash_core_Mass = min_core_Mass - 5D-3
         ELSE
            Post_He_Flash_core_Mass = he_core
         END IF
         DTY = timestep_decrement*DTY
         DT = CSY*DTY
         CALL Evolve1 ( JO, 10000, Check_Model, Post_He_Flash_Check_Model )
         CALL Complete_Model
         CALL Check_Final ( JO )
         Get_Post_He_Flash = JO .EQ. JO_AOK
         IF (building_He) THEN
            CALL Save_Model(fname)
            WRITE(*,*) 'Get_Post_He_Flash =', Get_Post_He_Flash
         END IF
      END IF
      CALL Restore_Control_Values(save_CTRLs)
      CALL Reset_Age_and_Model_Number(new_age + age_Offset, model+1)
      END FUNCTION Get_Post_He_Flash
      
      SUBROUTINE Save_Model(fname)
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      INTEGER, PARAMETER :: io_unit = io_flash
      INTEGER :: K
      IF (.NOT. OpenToWrite(fname,io_unit)) RETURN
      DO K = 1, N_SHELLs
         WRITE (io_unit, *) H(1:NUMV, K)
      END DO
      CLOSE (io_unit)
      END SUBROUTINE Save_Model
      
      LOGICAL FUNCTION Read_Model(fname)
      CHARACTER (LEN=strlen), INTENT(IN) :: fname
      INTEGER, PARAMETER :: io_unit = io_flash
      INTEGER :: K
      Read_Model = .FALSE.
      IF (.NOT. OpenToRead(fname,io_unit)) RETURN
      DO K = 1, N_SHELLs
         READ (io_unit, *) H(1:NUMV, K)
      END DO
      DH = 0D0
      CLOSE (io_unit)
      Read_Model = .TRUE.
      END FUNCTION Read_Model
      
      LOGICAL FUNCTION Build_He_Star ( init_Mass, Initialize_Params, Check_Model )
      USE ez_data, ONLY: init_M, init_Z
      USE ez_shell_data, ONLY: CZS
      DOUBLE PRECISION, INTENT(IN) :: init_Mass
      INTERFACE
         SUBROUTINE Initialize_Params
         END SUBROUTINE Initialize_Params
      END INTERFACE
      INTERFACE
         INTEGER FUNCTION Check_Model() ! return BACKUP, TERMINATE, or KEEP_GOING
         END FUNCTION Check_Model
      END INTERFACE
      Build_He_Star = .FALSE.
      IF (.NOT. Start_Read_Model ( init_Mass )) RETURN
      IF (.NOT. Get_Post_He_Flash ( init_Mass, init_Mass, 0d0, Check_Model )) RETURN
      IF (.NOT. Finish_Read_Model ( init_Mass, Initialize_Params, .TRUE. )) RETURN
      WRITE(*,*)
      IF ( .not. Trade_X_for_Y ( 1d0 - init_Z * 1.5, Check_Model ) ) THEN
         WRITE(*,*)  'Get_Post_He_Flash failed to create reduced hydrogen starting model'
         Build_He_Star = .FALSE.
         RETURN
      END IF
      Build_He_Star = .TRUE.
      END FUNCTION Build_He_Star
      
      END MODULE ez_flash
