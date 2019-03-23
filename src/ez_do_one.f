      MODULE ez_do_one
      USE ez_driver
      USE ez_log
      USE ez_do_one_utils
      USE star_data
      USE star_extras
      USE star_controls
      USE star_constants
      IMPLICIT NONE
      
      CONTAINS

      SUBROUTINE Get_Star_Info ( init_M, Build_Fname, Initial_Params, Before_Evolve, Check_Model )
      DOUBLE PRECISION, INTENT(IN) :: init_M
      INTERFACE
         SUBROUTINE Build_Fname(mass)
            DOUBLE PRECISION, INTENT(IN) :: mass
         END SUBROUTINE Build_Fname
      END INTERFACE
      INTERFACE
         SUBROUTINE Initial_Params
         END SUBROUTINE Initial_Params
      END INTERFACE
      INTERFACE
         SUBROUTINE Before_Evolve
         END SUBROUTINE Before_Evolve
      END INTERFACE
      INTERFACE
         INTEGER FUNCTION Check_Model() ! return BACKUP, TERMINATE, or KEEP_GOING
         END FUNCTION Check_Model
      END INTERFACE
      DOUBLE PRECISION :: mass, core, age_Offset
      DOUBLE PRECISION, PARAMETER :: age_Offset_Frac = 0.02D0
      WRITE (*,'(//,A,F5.1,2A,/)') '  Evolve ', init_M, ' Msolar star with data from ', metals_dir
      CALL Open_Logs(data_dir, init_M, Build_Fname, .TRUE.)
      CALL EZ_Start ( metals_dir )
      IF (.NOT. EZ_ZAMS(init_M, Initial_Params)) RETURN
      WRITE_PROFILE_TO_RUN = .TRUE.
      CALL Setup_Profilenames(init_M, Build_Fname)
      CALL Before_Evolve
      CALL EZ_Evolve(Check_Model)
      IF (do_post_He_flash .AND. helium_ignition) THEN
         WRITE(*,*)
         WRITE(*,*) '  Create a post Helium flash model'
         mass = star_Mass; core = mass_He_Core; age_Offset = age_Offset_Frac * star_Age
         phase_of_evolution = phase_GET_POST_FLASH
         IF (.NOT. EZ_Post_He_Flash ( mass, core, age_Offset, He_Flash_Do_One_Check_Model )) THEN
            WRITE (*,*) 'Failed in call to EZ_Post_He_Flash'
            RETURN
         END IF
         WRITE(*,'(4(A,F7.4))') '  Mass before ', mass, ', after', star_Mass, '  He core before ', core, ', after', mass_He_Core
         WRITE(*,*) '  Resume evolution post Helium flash'
         phase_of_evolution = phase_HE_IGNITION_OVER
         prev_AGE1 = -1D0
         WRITE(*,*)
         CALL EZ_Evolve(Check_Model)
      END IF
      
      IF (Model_Is_Okay()) THEN
        WRITE(*,'(a,i5)') 'Save profile for end of run:', model_Number
        CALL Record_Last_Model
      ELSE
        WRITE(*,*) 'Not saving profile for final model since it looks bogus.'
      END IF
      
      CALL Close_Logs ( .FALSE. )
      CALL Check_If_Do_Save
      IF ( WRITE_DIFFME ) CLOSE ( IO_DIFFME )
      END SUBROUTINE Get_Star_Info
      
      LOGICAL FUNCTION Model_Is_Okay()
      ! for now, just check for valid number in the final dynamic timescale
      Model_Is_Okay = ((dynamic_Timescale - dynamic_Timescale) .eq. 0d0) .and. ((dynamic_Timescale + 1d0) .gt. 1d0)
      END FUNCTION Model_Is_Okay
      
      SUBROUTINE Setup_Profilenames(init_M, Build_Fname)
      DOUBLE PRECISION, INTENT(IN) :: init_M
      INTERFACE
         SUBROUTINE Build_Fname(mass)
            DOUBLE PRECISION, INTENT(IN) :: mass
         END SUBROUTINE Build_Fname
      END INTERFACE
      IF ( WRITE_PROFILE_TO_RUN ) THEN
          IO_RUN_PROFILE = 38
          run_profilename = 'EZ_status.log'
      END IF
      IF ( WRITE_PROFILE_TO_DATA ) THEN
          IO_DATA_PROFILE = 39
          partial_name = 'model_'
          CALL Build_Fname(init_M)
          data_profile_prefix = full_name
      END IF
      IF ( WRITE_PROFILE_TO_MODELS ) THEN
          IO_MODEL_PROFILE = 40
          model_profile_prefix = 'models/model_'          
      END IF
      END SUBROUTINE Setup_Profilenames
      
      SUBROUTINE Dummy_Build_Filename(init_M)
      DOUBLE PRECISION, INTENT(IN) :: init_M
      full_name = trim(partial_name)
      END SUBROUTINE Dummy_Build_Filename
        
      INTEGER FUNCTION Do_One_Check_Model()
      LOGICAL :: must_do_log, stop_because_He_ignited
      DOUBLE PRECISION, PARAMETER :: log_He_Temp = 7.8D0, core_limit = 5d-3
      DOUBLE PRECISION, PARAMETER :: d_tau_min = 1D-2, d_tau_max = 1D0
      DOUBLE PRECISION, PARAMETER :: little_step_Factor = 10D0, little_step_Size = 10D0
      DOUBLE PRECISION :: d_tau
      INTEGER :: model
      INTEGER, PARAMETER :: tau_ramp = 50, num_little_steps_Limit = 100
      must_do_log = .FALSE.
      stop_because_He_ignited = .FALSE.
      model = model_Number
      IF ( WRITE_DIFFME .AND. model .GT. END_OF_DIFFME ) WRITE_DIFFME = .FALSE.
      IF ( model .EQ. 0 ) next_CNTR_RHO = min_CNTR_RHO
      IF (model .NE. stop_at_Model .AND. star_Age .LE. max_AGE) THEN 
         Do_One_Check_Model = KEEP_GOING
      ELSE 
         Do_One_Check_Model = TERMINATE
      END IF
      IF ( star_Age .LE. profile_AGE ) must_do_log = .TRUE. ! in case of backup, do it over
      IF ( save_helium_ignition .AND. (.NOT. helium_ignition) .AND. (log_center_Temp .GT. log_He_Temp) ) THEN
         CALL EZ_Extras
         IF ( power_Metal_burn + power_He_burn > power_Neutrinos ) THEN
            IF ( phase_of_evolution .NE. phase_GET_POST_FLASH ) THEN
               WRITE(*,'(a,i5)') 'Save profile for helium break-even reached:', model_Number
               must_do_log = .TRUE.
            END IF
            helium_ignition = .TRUE.
            phase_of_evolution = phase_HE_IGNITING
            ignition_Center_XHE = center_He
            he_Luminosity_Limit = log_Luminosity
            prev_Luminosity = log_Luminosity
            IF ( do_post_He_flash ) stop_because_He_ignited = .TRUE.
         END IF
      END IF
      IF ( (phase_of_evolution .EQ. phase_HE_IGNITION_OVER .AND. prev_AGE1 .EQ. -1D0) .OR. star_Age .LE. post_He_AGE ) THEN
         ! need to check the age since may backup and over the previous saved info
         prev_TCNTR1 = log_center_Temp; prev_TCNTR2 = prev_TCNTR1
         prev_AGE1 = star_Age; prev_AGE2 = prev_AGE1
         must_do_log = .TRUE.
         post_He_AGE = star_Age
         WRITE(*,'(a,i5)') 'Save profile for starting phase of helium burning:', model_Number
      ELSE IF ( Time_To_Profile() ) THEN
         must_do_log = .TRUE.
      ELSE IF ( model_Number .EQ. profile_model ) THEN
         must_do_log = .TRUE.
         WRITE(*,'(a,i5)') 'Save profile for model number:', model_Number
      END IF
      IF ( Log_State ( must_do_log ) ) THEN
         CALL Write_Logs
         CALL Write_Special_Log
         IF ( must_do_log ) CALL Save_Profiles
         IF (STOP_Sign()) THEN
            Do_One_Check_Model = TERMINATE
         END IF
      END IF
      IF ( stop_because_He_ignited ) Do_One_Check_Model = TERMINATE
      IF ( star_Mass - mass_He_Core .LE. core_limit .and. wind_Eta .ne. 0d0 ) THEN
         WRITE(*,'(a,i5)') 'Turning off wind:', model_Number
         wind_Eta = 0D0  ! turn off wind when envelope gone
      END IF
      IF ( target_tau_Photosphere .GT. 0D0 .AND. tau_Photosphere .NE. target_tau_Photosphere .AND. model .GT. tau_Model) THEN
         d_tau = model - tau_Model
         d_tau = (target_tau_Photosphere - tau_Photosphere) * min(d_tau_max, max(d_tau_min, d_tau / tau_ramp))
         tau_Photosphere = tau_Photosphere + d_tau
         IF ( tau_Photosphere .GE. target_tau_Photosphere ) THEN
            WRITE(*,'(I5,3X,A,1X,F8.2)') model, 'Have reached target tau_Photosphere =', tau_Photosphere
            tau_Photosphere = target_tau_Photosphere
         END IF
      END IF
      IF ( time_Step .GT. little_step_Size .OR. time_Step * CSY .GT. little_step_Factor * dynamic_Timescale ) THEN
         number_of_little_steps = 0
      ELSE
         number_of_little_steps = number_of_little_steps + 1
         IF ( number_of_little_steps .GT. num_little_steps_Limit ) THEN
            WRITE(*,*) 'Stopping because of too many little time steps'
            Do_One_Check_Model = TERMINATE
         END IF
      END IF
      extra_Mdot_param = -1e-14
      END FUNCTION Do_One_Check_Model
      
      LOGICAL FUNCTION Time_To_Profile(GoOn)
      ! end-of-run and helium break-even are always done.  this function decides on other models to be profiled.
      LOGICAL, INTENT(INOUT), OPTIONAL :: GoOn
      DOUBLE PRECISION, PARAMETER :: center_He_going = 5D-2, center_He_drop = 1D-2, surface_T_drop = 4D-2
      DOUBLE PRECISION, PARAMETER :: offset_L = 0.15D0, center_H_gone = 1D-6, center_H_going = 1D-1, lower_CNTR_RHO_limit = 3D0
      Time_To_Profile = .FALSE.
      IF ( do_CNTR_RHOs ) THEN
         IF (log_center_Density .GT. next_CNTR_RHO ) THEN
            Time_To_Profile = .TRUE.
            WRITE(*,'(a,i5)') 'Save profile for log10 center density reaching new level:', model_Number
            next_CNTR_RHO = next_CNTR_RHO + del_CNTR_RHO
            RETURN
         END IF
      END IF
      SELECT CASE ( phase_of_evolution )
      CASE ( phase_STARTING )
         IF ( center_H .LT. center_H_going ) THEN
            Time_To_Profile = .TRUE.
            prev_TSURF = log_surface_Temp
            phase_of_evolution = phase_END_ZAMS
            WRITE(*,'(a,i5)') 'Save profile for leaving main sequence:', model_Number
            IF (stop_at_Model .EQ. -5) GoOn = .false. !DarkStars add - Pat Scott 080426
         END IF
      CASE ( phase_END_ZAMS )
         IF ( center_H .LT. center_H_gone .AND. log_surface_Temp .LT. prev_TSURF-surface_T_drop ) THEN
            Time_To_Profile = .TRUE.
            phase_of_evolution = phase_WAIT_FOR_HE
            WRITE(*,'(a,i5)') 'Save profile for center hydrogen exhausted and surface temperature dropping:', model_Number
            IF (stop_at_Model .EQ. -6 .OR. stop_at_Model .EQ. -10) GoOn = .false. !DarkStars add - Pat Scott 080426
         END IF
      CASE ( phase_WAIT_FOR_HE )
      CASE ( phase_HE_IGNITING ) ! for non-flash ignition of helium core
         IF ( center_He .LE. ignition_Center_XHE-center_He_drop .AND. log_Luminosity .GT. prev_Luminosity ) THEN
            Time_To_Profile = .TRUE.
            phase_of_evolution = phase_HE_IGNITION_OVER
            WRITE(*,'(a,i5)') 'Save profile for center helium down and luminosity rising:', model_Number
            prev_TCNTR2 = prev_TCNTR1; prev_AGE2 = prev_AGE1
            prev_TCNTR1 = log_center_Temp; prev_AGE1 = star_Age
            IF ( log_Luminosity .GT. he_Luminosity_Limit ) he_Luminosity_Limit = log_Luminosity
            IF (stop_at_Model .EQ. -7) GoOn = .false. !DarkStars add - Pat Scott 080426
         END IF
         prev_Luminosity = log_Luminosity
      CASE ( phase_HE_IGNITION_OVER )
         IF ( .NOT. do_CNTR_T_drops ) THEN
            IF ( center_He .LT. center_He_going .AND. log_Luminosity .GT. he_Luminosity_Limit+offset_L ) THEN
               Time_To_Profile = .TRUE.
               phase_of_evolution = phase_FINISHING
               WRITE(*,'(a,i5)') 'Save profile for center helium low:', model_Number
               IF (stop_at_Model .EQ. -9) GoOn = .false. !DarkStars add - Pat Scott 080426
            END IF      
         ELSEIF ( log_center_Temp .LT. prev_TCNTR1 .AND. star_Age .GT. prev_AGE1 ) THEN ! center temp starting to fall
            IF ( prev_TCNTR1 .LT. prev_TCNTR2 .AND. prev_AGE1 .GT. prev_AGE2 ) THEN ! check for two in a row to avoid false alarms
               IF ( log_Luminosity .GT. he_Luminosity_Limit+offset_L ) THEN ! try to avoid text overlap on the HR diagram
                  Time_To_Profile = .TRUE.
                  phase_of_evolution = phase_FINISHING
                  WRITE(*,'(a,i5)') 'Save profile for center temperature starting to fall:', model_Number
               END IF
            END IF
         END IF
         prev_TCNTR2 = prev_TCNTR1; prev_AGE2 = prev_AGE1
         prev_TCNTR1 = log_center_Temp; prev_AGE1 = star_Age
      CASE ( phase_GET_POST_FLASH )
      CASE ( phase_FINISHING )
      CASE DEFAULT
         STOP 'Error in Time_To_Profile'
      END SELECT
      END FUNCTION Time_To_Profile
      
      SUBROUTINE Record_Last_Model
      LOGICAL :: logged
      logged = Log_State ( .TRUE. )
      CALL Write_Logs
      CALL Write_Special_Log
      CALL Save_Profiles
      END SUBROUTINE Record_Last_Model
      
      SUBROUTINE Write_Special_Log
      ! Modified for DarkStars - Pat Scott 20070905, 20080415
      USE DkStrs_data
      USE DkStrs_WIMPdens
      USE DkStrs_transport
      INTEGER, PARAMETER :: num_out1 = 7, num_out2 = 23
      DOUBLE PRECISION :: out1(num_out1), out2(num_out2), WIMPmassFrac, WIMPenergyinject, K
      IF (.TRUE.) THEN ! change to .TRUE. if you want a special log
         CALL get_condEff
         out1 = (/ star_Age, log_Luminosity, log_Radius, log_surface_Temp, star_Mass, mass_He_Core, 
     &      radius_He_Core /)
            ! galr                                                 Distance from the Galactic Centre [pc]
            ! rhowimp                                              Local ambient WIMP density [GeV/cc]
            ! v_star                                               Stellar velocity through halo [km/s]
            ! galesc                                               Local galactic escape velocity [km/s]
            ! n_WIMPs                                              Total number of WIMPs [WIMPs]
            WIMPmassFrac = n_WIMPs * mx * CGGEV / CMSN * 1.d-33!   Fraction of mass in WIMPs [dimensionless]
            ! cap                                                  Capture rate [WIMPs/yr]
            ! 2.*ann                                               Annihilation rate [WIMPs/yr]
            ! dNdt                                                 Total Net change in WIMP number [WIMPs/yr] 
            WIMPenergyinject = mx * 2. *ann * CMEV * 1.d3 / CSY!   Total luminosity from annihilations [erg/s]
            ! r_chi                                                WIMP characteristic radius [cm]
            ! Tw                                                   WIMP isothermal temperature [K]
            ! condEff                                              Dimensionless WIMP conductive effectiveness
            ! WIMPdens(0.d0)                                       Central WIMP density [WIMPs/cm^3]
            ! WIMPdens(r_chi)                                      WIMP density at characteristic radius [WIMPs/cm^3] 
            IF (DoWimpyThings) THEN
              K = mfp(1)/r_chi 
            ELSE 
              K = 0.d0
            ENDIF!                                                 Knudsen parameter [dimensionless]
            out2 = (/ out1, galr, rhowimp, v_star, galesc, n_WIMPs, WIMPmassFrac, cap, 2.*ann, dNdt, WIMPenergyinject, r_chi, Tw, 
     &       condEff, WIMPdens(0.d0), WIMPdens(r_chi), K /)
            CALL Write_Extra_Log(num_out2, out2) ! extra.log holds the output
      END IF
      END SUBROUTINE Write_Special_Log

      SUBROUTINE Save_Profiles
      IF ( WRITE_PROFILE_TO_DATA ) THEN
         CALL Setup_Profilename(data_profile_prefix, data_profilename, model_Number)
         CALL Write_Status_Info(IO_DATA_PROFILE, data_profilename, model_Number)
      END IF
      IF ( WRITE_PROFILE_TO_MODELS ) THEN
         CALL Setup_Profilename(model_profile_prefix, model_profilename, model_Number)
         CALL Write_Status_Info(IO_MODEL_PROFILE, model_profilename, model_Number)
      END IF
      CALL Write_Profiles
      profile_AGE = star_Age
      recent_profile_model = model_Number
      IF ( call_save_exp_info ) CALL Save_Experiment_Info
      END SUBROUTINE Save_Profiles
      
      SUBROUTINE Check_If_Do_Save
      INTEGER :: model
      model = model_Number
      IF ( model .LT. stop_at_Model .OR. stop_at_Model .LT. 0 ) RETURN
      IF ( .NOT. EZ_Save(save_filename) ) THEN
         WRITE(*,'(a,i5)') 'FAILED trying to SAVE at model number', model
      ELSE
         WRITE(*,'(a,i5)') 'SAVED for later debugging at model number', model
      END IF
      END SUBROUTINE Check_If_Do_Save

      SUBROUTINE DkStrs_Check_If_Do_Save
      ! DarkStars version of Check_If_Do_Save
      USE DkStrs_admin
      INTEGER :: model
      character (len=strlen) :: EZ_filename, DarkStars_filename
      model = model_Number
      IF ( model .LT. stop_at_Model .OR. stop_at_Model .LT. 0 ) RETURN
      EZ_filename = trim(main_dir)//'/'//trim(save_filename)
      DarkStars_filename = trim(EZ_filename)//trim(DarkStars_extension)
      IF (EZ_Save(EZ_filename) .AND. DarkStars_Save(DarkStars_filename) ) THEN
         WRITE(*,'(a,i5)') 'SAVED for later debugging at model number', model
      ELSE
         WRITE(*,'(a,i5)') 'FAILED trying to SAVE at model number', model
      END IF
      END SUBROUTINE DkStrs_Check_If_Do_Save
      
      INTEGER FUNCTION Bare_Bones_Check_Model()
      LOGICAL :: logged
      logged = Log_State( .FALSE. )
      Bare_Bones_Check_Model = KEEP_GOING
      END FUNCTION Bare_Bones_Check_Model
      
      INTEGER FUNCTION He_Flash_Do_One_Check_Model()
      LOGICAL :: logged
      logged = Log_State ( .FALSE. )
      He_Flash_Do_One_Check_Model = KEEP_GOING
      END FUNCTION He_Flash_Do_One_Check_Model
      
      
      LOGICAL FUNCTION STOP_Sign()
      ! Added extra run-specific STOP_sign for simultaneous runs - Pat Scott 20070905
      INTEGER :: ios, ios_specific
      STOP_Sign = .FALSE.
      OPEN(UNIT=99, FILE='stop_sign', ACTION='READ', STATUS='OLD', IOSTAT=ios)
      OPEN(UNIT=98, FILE=TRIM(data_dir)//'/../stop_sign', ACTION='READ', STATUS='OLD', IOSTAT=ios_specific)
      IF (ios == 0 .OR. ios_specific == 0) THEN
         WRITE(*,*) "Have found file named 'stop_sign', so terminate run."
         STOP_Sign = .TRUE.
      END IF
      END FUNCTION STOP_Sign
      
      SUBROUTINE Dummy_Before_Evolve
      END SUBROUTINE Dummy_Before_Evolve

      integer function Do_Delta_Check_Model()
         ! save profiles if any of the watched params has changed by more than a specified amount
         
         Do_Delta_Check_Model = KEEP_GOING
         
         if (check_one_delta(log_Luminosity, del_log_Luminosity, prv_log_Luminosity)) then
            write(*,'(a,i5)') 'Save profile because of large change in luminosity:', model_Number
            return
         end if
         if (check_one_delta(log_surface_Temp, del_log_surface_Temp, prv_log_surface_Temp)) then
            write(*,'(a,i5)') 'Save profile because of large change in surface temperature:', model_Number
            return
         end if
         if (check_one_delta(log_center_Temp, del_log_center_Temp, prv_log_center_Temp)) then
            write(*,'(a,i5)') 'Save profile because of large change in center temperature:', model_Number
            return
         end if
         if (check_one_delta(log_center_Density, del_log_center_Density, prv_log_center_Density)) then
            write(*,'(a,i5)') 'Save profile because of large change in center density:', model_Number
            return
         end if

         contains
         
         logical function check_one_delta(val, del, prv)
            double precision, intent(IN) :: val, del, prv
            
            if (abs(val-prv) < del) then
               check_one_delta = .false.
               return
            end if
            
            CALL EZ_Extras
            CALL Save_Profiles
            prv_log_Luminosity = log_Luminosity
            prv_log_surface_Temp = log_surface_Temp
            prv_log_center_Temp = log_center_Temp
            prv_log_center_Density = log_center_Density
            check_one_delta = .true.
         
         end function check_one_delta
         
      end function Do_Delta_Check_Model


      SUBROUTINE DkStrs_Get_Star_Info ( init_M, Build_Fname, Build_Extra_Fname, Initial_Params, Before_Evolve, Check_Model )
      ! DarkStars-specific version of Get_Star_Info
      DOUBLE PRECISION, INTENT(IN) :: init_M
      INTERFACE
         SUBROUTINE Build_Fname(mass)
            DOUBLE PRECISION, INTENT(IN) :: mass
         END SUBROUTINE Build_Fname
      END INTERFACE
      INTERFACE
         SUBROUTINE Build_Extra_Fname(mass, what)
            DOUBLE PRECISION, INTENT(IN) :: mass
            INTEGER, INTENT(IN) :: what
         END SUBROUTINE Build_Extra_Fname
      END INTERFACE
      INTERFACE
         SUBROUTINE Initial_Params
         END SUBROUTINE Initial_Params
      END INTERFACE
      INTERFACE
         SUBROUTINE Before_Evolve
         END SUBROUTINE Before_Evolve
      END INTERFACE
      INTERFACE
         INTEGER FUNCTION Check_Model() ! return BACKUP, TERMINATE, or KEEP_GOING
         END FUNCTION Check_Model
      END INTERFACE
      DOUBLE PRECISION :: mass, core, age_Offset
      DOUBLE PRECISION, PARAMETER :: age_Offset_Frac = 0.02D0
      WRITE (*,'(//,A,F5.1,2A,/)') '  Evolve ', init_M, ' Msolar star with data from ', metals_dir
      CALL Open_Logs(data_dir, init_M, Build_Fname, .TRUE.)
      CALL EZ_Start ( metals_dir )
      IF (.NOT. EZ_ZAMS(init_M, Initial_Params)) RETURN
      WRITE_PROFILE_TO_RUN = .TRUE.
      CALL DkStrs_Setup_Profilenames(init_M, Build_Fname, Build_Extra_Fname)
      CALL Before_Evolve
      CALL EZ_Evolve(Check_Model)
      IF (do_post_He_flash .AND. helium_ignition) THEN
         WRITE(*,*)
         WRITE(*,*) '  Create a post Helium flash model'
         mass = star_Mass; core = mass_He_Core; age_Offset = age_Offset_Frac * star_Age
         phase_of_evolution = phase_GET_POST_FLASH
         IF (.NOT. EZ_Post_He_Flash ( mass, core, age_Offset, He_Flash_Do_One_Check_Model )) THEN
            WRITE (*,*) 'Failed in call to EZ_Post_He_Flash'
            RETURN
         END IF
         WRITE(*,'(4(A,F7.4))') '  Mass before ', mass, ', after', star_Mass, '  He core before ', core, ', after', mass_He_Core
         WRITE(*,*) '  Resume evolution post Helium flash'
         phase_of_evolution = phase_HE_IGNITION_OVER
         prev_AGE1 = -1D0
         WRITE(*,*)
         CALL EZ_Evolve(Check_Model)
      END IF
      
      IF (Model_Is_Okay()) THEN
        WRITE(*,'(a,i5)') 'Save profile for end of run:', model_Number
        CALL Record_Last_Model
      ELSE
        WRITE(*,*) 'Not saving profile for final model since it looks bogus.'
      END IF
      
      CALL Close_Logs ( .TRUE. )
      CALL DkStrs_Check_If_Do_Save
      IF ( WRITE_DIFFME ) CLOSE ( IO_DIFFME )
      END SUBROUTINE DkStrs_Get_Star_Info
      
   
      SUBROUTINE DkStrs_Setup_Profilenames(init_M, Build_Fname, Build_Extra_Fname)
      ! DarkStars-specific version of Setup_Profilenames 
      DOUBLE PRECISION, INTENT(IN) :: init_M
      INTERFACE
         SUBROUTINE Build_Fname(mass)
            DOUBLE PRECISION, INTENT(IN) :: mass
         END SUBROUTINE Build_Fname
      END INTERFACE
      INTERFACE
         SUBROUTINE Build_Extra_Fname(mass, what)
            DOUBLE PRECISION, INTENT(IN) :: mass
            INTEGER, INTENT(IN) :: what
         END SUBROUTINE Build_Extra_Fname
      END INTERFACE
      IF ( WRITE_PROFILE_TO_RUN ) THEN
          IO_RUN_PROFILE = 38
          run_profilename = trim(data_dir)//'/EZ_status.log'
      END IF
      IF ( WRITE_PROFILE_TO_DATA ) THEN
          IO_DATA_PROFILE = 39
          partial_name = 'model_'
          CALL Build_Fname(init_M)
          partial_name = full_name
          CALL Build_Extra_Fname(init_M,1)
          data_profile_prefix = full_name
      END IF
      IF ( WRITE_PROFILE_TO_MODELS ) THEN
          IO_MODEL_PROFILE = 40
          partial_name = 'model_'
          CALL Build_Extra_Fname(init_M,2)
          model_profile_prefix = full_name
      END IF
      END SUBROUTINE DkStrs_Setup_Profilenames


      END MODULE ez_do_one
      
