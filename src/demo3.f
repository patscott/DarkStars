      MODULE demo3
      USE ez_driver
      USE star_data
      USE star_controls
      USE star_constants
      USE ez_log
      USE ez_do_one
      IMPLICIT NONE

      DOUBLE PRECISION :: L_SAVE ! Check_Model_Up_RGB uses this to decide when to stop and save state on way up RGB
      ! Check_Model_While_Eject uses the following to control the envelope ejection
      DOUBLE PRECISION :: L_LIMIT1 ! luminosity at which begin to turn on the mass loss
      DOUBLE PRECISION :: L_LIMIT2 ! luminosity at which finish turning on the mass loss
      DOUBLE PRECISION :: E_LIMIT1 ! envelope mass at which start turning off mass loss
      DOUBLE PRECISION :: E_LIMIT2 ! envelope mass at which stop mass loss
      DOUBLE PRECISION :: TC_LIMIT ! center temp when save the WD model
      DOUBLE PRECISION :: CMI_VALUE ! value for the artificial mass loss parameter, CMI
      DOUBLE PRECISION :: PREHEAT_TARGET, COOLING_TARGET
      DOUBLE PRECISION :: MDOT_TARGET ! for when accreting mass onto the WD
      DOUBLE PRECISION :: mdot_start_age, mdot_start_mass
      LOGICAL :: OK_TO_SAVE
      LOGICAL :: SAVING_PROFILES
      LOGICAL :: MASS_LOSS_ON
      INTEGER :: start_mass_change_model
      INTEGER, PARAMETER :: WRTNUM=9 ! max number of profiles to write to file
      DOUBLE PRECISION :: MWRITE(WRTNUM), MWRTAGE(WRTNUM)   ! control data for deciding when to write a profile

      CONTAINS
      
      INTEGER FUNCTION Raise_WD_tau_Check_Model()
      LOGICAL :: logged
      DOUBLE PRECISION, PARAMETER :: d_tau_min = 1D-2, d_tau_max = 1D0
      DOUBLE PRECISION :: d_tau
      INTEGER, PARAMETER :: tau_ramp = 50
      logged = Log_State ( .FALSE. )
      Raise_WD_tau_Check_Model = KEEP_GOING
      IF (tau_Model .EQ. 0) tau_Model = model_Number - 1
      d_tau = model_Number - tau_Model
      d_tau = (target_tau_Photosphere - tau_Photosphere) * min(d_tau_max, max(d_tau_min, d_tau / tau_ramp))
      tau_Photosphere = tau_Photosphere + d_tau
      IF ( tau_Photosphere .GE. target_tau_Photosphere ) THEN
         WRITE(*,'(I5,3X,A,1X,F8.2)') model_Number, 'Have reached target tau_Photosphere =', tau_Photosphere
         tau_Photosphere = target_tau_Photosphere
         Raise_WD_tau_Check_Model = TERMINATE
      END IF
      END FUNCTION Raise_WD_tau_Check_Model
      
      INTEGER FUNCTION Restore_WD_tau_Check_Model()
      LOGICAL :: logged
      DOUBLE PRECISION, PARAMETER :: d_tau_min = 1D-2, d_tau_max = 1D0
      DOUBLE PRECISION :: d_tau
      INTEGER, PARAMETER :: tau_ramp = 50
      logged = Log_State ( .FALSE. )
      Restore_WD_tau_Check_Model = KEEP_GOING
      IF (tau_Model .EQ. 0) tau_Model = model_Number - 1
      d_tau = model_Number - tau_Model
      d_tau = (tau_Photosphere - target_tau_Photosphere) * min(d_tau_max, max(d_tau_min, d_tau / tau_ramp))
      tau_Photosphere = tau_Photosphere - d_tau
      IF ( tau_Photosphere .LE. target_tau_Photosphere ) THEN
         WRITE(*,'(I5,3X,A,1X,F8.2)') model_Number, 'Have reached target tau_Photosphere =', tau_Photosphere
         tau_Photosphere = target_tau_Photosphere
         Restore_WD_tau_Check_Model = TERMINATE
      END IF
      END FUNCTION Restore_WD_tau_Check_Model
      
      INTEGER FUNCTION Check_Model_Up_RGB()
      LOGICAL :: logged
      logged = Log_State ( .FALSE. )
      Check_Model_Up_RGB = TERMINATE
      IF ( log_Luminosity .GT. L_SAVE ) THEN ! stop when luminosity reaches goal for saving
         OK_TO_SAVE = .TRUE.
         RETURN
      END IF
      Check_Model_Up_RGB = KEEP_GOING
      END FUNCTION Check_Model_Up_RGB
      
      INTEGER FUNCTION Do_WD_Check_Model()
      USE ez_cycle_data, ONLY : REQUESTED_DTY
      USE ez_state_data, ONLY : GAM
      USE ez_data, ONLY : H
      LOGICAL :: logged
      DOUBLE PRECISION :: mass_ramp_size, speedbump
      DOUBLE PRECISION :: fac
      Do_WD_Check_Model = TERMINATE
      logged = Log_State ( .FALSE. )
      CMI_VALUE = MDOT_TARGET / star_Mass
      speedbump = -DLOG10(MDOT_TARGET)
      mass_ramp_size = 1000D0
      !IF ( speedbump .LE. 9D0 ) mass_ramp_size = 500
      IF ( extra_Mdot_param .EQ. 0 ) THEN
         start_mass_change_model = model_Number - 1
         mdot_start_age = star_Age
         mdot_start_mass = star_Mass
         WRITE(*,*) 'mass_ramp_size', mass_ramp_size
      END IF
      fac = (model_Number-start_mass_change_model)/mass_ramp_size
      IF ( fac .LT. 1D0 ) THEN
         fac = fac*fac
         extra_Mdot_param = CMI_VALUE*fac*fac
         WRITE(*,*) 'extra_Mdot_param', extra_Mdot_param
         REQUESTED_DTY = 10D0
      ELSE
         IF (extra_Mdot_param .NE. CMI_VALUE) WRITE(*,*) 'Full speed mass gain    fac', fac
         extra_Mdot_param = CMI_VALUE
      END IF
      Do_WD_Check_Model = KEEP_GOING
      END FUNCTION Do_WD_Check_Model
      
      INTEGER FUNCTION Check_Model_To_Stop()
      LOGICAL :: logged
      logged = Log_State ( .FALSE. )
      IF ( SAVING_PROFILES ) Call Save_Profiles
      Check_Model_To_Stop = KEEP_GOING
      END FUNCTION Check_Model_To_Stop
      
      INTEGER FUNCTION Check_Model_WD()
      LOGICAL :: logged
      logged = Log_State ( .FALSE. )
      IF ( SAVING_PROFILES ) Call Save_Profiles
      IF ( log_center_Temp .LE. TC_LIMIT ) THEN
         OK_TO_SAVE = .TRUE.
         Check_Model_WD = TERMINATE
         RETURN 
      END IF
      Check_Model_WD = KEEP_GOING
      END FUNCTION Check_Model_WD
      
      INTEGER FUNCTION Preheat_WD_Check_Model()
      USE ez_cycle_data, ONLY : REQUESTED_DTY
      LOGICAL :: logged
      DOUBLE PRECISION, PARAMETER :: DX1=1D0, DX2=0.42D0
      DOUBLE PRECISION :: fac
      logged = Log_State ( .FALSE. )
      IF ( logged ) WRITE(*,*) 'power_Neutrinos', power_Neutrinos, 'neutrino_Cooling_Factor', neutrino_Cooling_Factor
      IF ( (center_Degeneracy - PREHEAT_TARGET) .LE. DX2 ) THEN
         WRITE(*,*) 'center_Degeneracy', center_Degeneracy, 'PREHEAT_TARGET', PREHEAT_TARGET
         extra_Energy_param = extra_Energy_param * 0.9
         extra_Energy_max = extra_Energy_param
         neutrino_Cooling_Factor = 1D0
         OK_TO_SAVE = .TRUE.
         IF ( extra_Energy_param .LT. 1D-2 ) Preheat_WD_Check_Model = TERMINATE
         RETURN 
      END IF
      neutrino_Cooling_Factor = -1D3
      extra_Energy_param = 1D1
      REQUESTED_DTY = 50D0
      IF ( (center_Degeneracy - PREHEAT_TARGET) .LE. DX1 ) THEN
         fac = (center_Degeneracy - PREHEAT_TARGET) / DX1
         REQUESTED_DTY = 1D2 * fac
         IF ( (center_Degeneracy - PREHEAT_TARGET) .LE. DX1*0.1D0 ) neutrino_Cooling_Factor = 1D0 - sqrt(fac)
         IF ( REQUESTED_DTY .LT. 1D1 ) REQUESTED_DTY = 1D1
         fac = fac * fac * fac
         extra_Energy_param = extra_Energy_param * fac
         IF ( extra_Energy_param .LT. 1D0 ) extra_Energy_param = 1D0
         WRITE(*,*) 'extra_Energy_param', extra_Energy_param, 'neutrino_Cooling_Factor', neutrino_Cooling_Factor
      END IF
      extra_Energy_max = extra_Energy_param
      extra_Energy_time = 1D0
      Preheat_WD_Check_Model = KEEP_GOING
      END FUNCTION Preheat_WD_Check_Model
      
      INTEGER FUNCTION After_Preheat_WD_Check_Model()
      LOGICAL :: logged
      DOUBLE PRECISION, PARAMETER :: DX=5D-1
      logged = Log_State ( .FALSE. )
      IF ( (center_Degeneracy - PREHEAT_TARGET) .GT. DX ) THEN
         WRITE(*,*) 'center_Degeneracy', center_Degeneracy, 'PREHEAT_TARGET', PREHEAT_TARGET
         OK_TO_SAVE = .TRUE.
         After_Preheat_WD_Check_Model = TERMINATE
         RETURN 
      END IF
      After_Preheat_WD_Check_Model = KEEP_GOING
      END FUNCTION After_Preheat_WD_Check_Model
      
      INTEGER FUNCTION Cool_WD_Check_Model()
      LOGICAL :: logged
      logged = Log_State ( .FALSE. )
      IF ( log_center_Temp .LE. COOLING_TARGET ) THEN
         WRITE(*,*) 'log_center_Temp', log_center_Temp, 'COOLING_TARGET', COOLING_TARGET
         OK_TO_SAVE = .TRUE.
         Cool_WD_Check_Model = TERMINATE
         RETURN 
      END IF
      CALL EZ_Extras
      Cool_WD_Check_Model = KEEP_GOING
      END FUNCTION Cool_WD_Check_Model
      
      INTEGER FUNCTION Remove_H_Check_Model()
      USE star_data
      USE star_controls
      USE ez_data
      DOUBLE PRECISION :: X, shellDX
      DOUBLE PRECISION, PARAMETER :: DX=0.5D0
      INTEGER :: IK
      LOGICAL :: logged
      Remove_H_Check_Model = KEEP_GOING
      X = star_Mass_H / star_Mass
      WRITE(*,*) 'star_Mass_H', star_Mass_H, 'X', X
      IF ( X .LE. 1D-30 ) THEN
         Remove_H_Check_Model = TERMINATE
         DO IK = 1, N_SHELLs
            H(V_XH,IK) = 0D0
            H(V_XHE,IK) = H(V_XHE,N_CNTR_shell)
            H(V_XC,IK) = H(V_XC,N_CNTR_shell)
            H(V_XO,IK) = H(V_XO,N_CNTR_shell)
            H(V_XNE,IK) = H(V_XNE,N_CNTR_shell)
         END DO
         logged = Log_State ( .TRUE. )
      ELSE
         logged = Log_State ( .FALSE. )
      END IF
      IF ( X .GE. 1D-6 ) THEN
         DO IK = 1, N_SHELLs
            shellDX = H(V_XH,IK)*DX
            H(V_XHE,IK) = H(V_XHE,IK) + shellDX
            H(V_XH,IK) = H(V_XH,IK) - shellDX
         END DO
      ELSE
         DO IK = 1, N_SHELLs
            shellDX = H(V_XH,IK)
            H(V_XHE,IK) = H(V_XHE,IK) + shellDX
            H(V_XH,IK) = 0D0
         END DO
      END IF
      DO IK = 1, N_SHELLs
         H(V_XC,IK) = H(V_XC,N_CNTR_shell)
         H(V_XO,IK) = H(V_XO,N_CNTR_shell)
         H(V_XNE,IK) = H(V_XNE,N_CNTR_shell)
      END DO
      END FUNCTION Remove_H_Check_Model
      
      INTEGER FUNCTION Check_Model_While_Eject()
      DOUBLE PRECISION :: envMass
      LOGICAL :: logged
      Check_Model_While_Eject = TERMINATE
      logged = Log_State ( .FALSE. )
      IF ( log_Luminosity .GT. L_LIMIT1 .AND. extra_Mdot_param .NE. CMI_VALUE ) THEN ! gradually turn on the mass loss
         IF ( log_Luminosity .LT. L_LIMIT2 ) THEN
            extra_Mdot_param = CMI_VALUE*(log_Luminosity-L_LIMIT1)/(L_LIMIT2-L_LIMIT1)
            IF ( .NOT. MASS_LOSS_ON ) THEN
               WRITE (*,'(A,I4,A,F6.4)') '  >>>>>>>> Turning on Mdot at model ', model_Number, ' with log(L) ', log_Luminosity
               WRITE (*,*) 'star_Age', star_Age
               MASS_LOSS_ON = .TRUE.
            END IF
         ELSE
            extra_Mdot_param = CMI_VALUE
            WRITE (*,'(A,I4,A,F6.4)') '  >>>>>>>> Mass loss fully on at model ', model_Number, ' with log(L) ', log_Luminosity
         END IF
      END IF
      envMass = star_Mass-mass_He_Core
      IF ( envMass .LT. E_LIMIT1 ) THEN ! envelope is less than limit; gradually stop the mass loss
         IF ( extra_Mdot_param .EQ. CMI_VALUE ) THEN
            L_LIMIT1 = 1D99 ! so don't try to turn mass loss on again
            IF ( MASS_LOSS_ON ) THEN
               MASS_LOSS_ON = .FALSE.
               WRITE (*,'(A,I4,A,F6.4)') '  >>>>>>>> Turning off Mdot at model ', model_Number, ' with envelope mass ', envMass
            END IF
         END IF
         IF ( extra_Mdot_param .NE. 0 ) THEN
            IF ( envMass .GT. E_LIMIT2 ) THEN ! reduce mass loss along with envelope mass
               extra_Mdot_param = CMI_VALUE*(envMass-E_LIMIT2)/(E_LIMIT1-E_LIMIT2)
            ELSE ! finished with mass loss
               WRITE (*,'(A,I4,A,F6.4)') '  >>>>>>>> Mass loss off at model ', model_Number, ' with envelope mass ', envMass
               WRITE (*,*) 'star_Age', star_Age
               extra_Mdot_param = 0D0; OK_TO_SAVE = .TRUE.; RETURN ! terminate run so can save state
            END IF
         END IF
      END IF
      !IF ( SAVING_PROFILES ) Call Save_Profiles3
      Check_Model_While_Eject = KEEP_GOING
      END FUNCTION Check_Model_While_Eject
      
      INTEGER FUNCTION Check_Model_For_Shell_Flashes()
      LOGICAL :: logged
      logged = Log_State ( .FALSE. )
      !IF ( logged ) WRITE (IO_XLOG, *) power_H_burn, size_M_Conv, star_Mdot
      Check_Model_For_Shell_Flashes = KEEP_GOING
      END FUNCTION Check_Model_For_Shell_Flashes
      
      SUBROUTINE Initial_Params
      overshoot_param_H = 0D0 ! turn off convective overshooting for these runs
      END SUBROUTINE Initial_Params
      
      SUBROUTINE Demo_3
      USE ez_data, ONLY: stats_N, stats_X, stats_X2
      USE ez_cycle_data, ONLY: PSI_limit
      USE ez_flash, ONLY: Save_Model
      CHARACTER (LEN=strlen) :: metals_dir, save_filename
      CHARACTER (LEN=strlen) :: wd_filename, wd_RGB_name, wd_eject_name, wd_adj_name, wd_tau_up_name
      CHARACTER (LEN=strlen) :: wd_cool_name, wd_tau_down_name, wd_add_h_name
      DOUBLE PRECISION, PARAMETER :: start_DL = 0.55, finish_DL = 0.50  ! luminosity parameters for turning on mass loss
      DOUBLE PRECISION, PARAMETER :: start_E=4D-3, finish_E=1.4D-3        ! envelope mass parameters for turning off mass loss
      DOUBLE PRECISION, PARAMETER :: init_M = 1D0, CMI_highwind = -0.620D-6, CMI_lowwind = -0.20D-6, CMI_shell_flashes = -0.5D-6
      DOUBLE PRECISION, PARAMETER :: eta1 = 11D0, eta2 = 13D0  ! eta1 doesn't have shell flashes.  eta2 does.
      DOUBLE PRECISION, PARAMETER :: CMI_gain = 1.25D-8   ! 1.5D-8 gives He flash at 0.89 Msun
      LOGICAL :: logged
      INTEGER :: where_saved
      DOUBLE PRECISION :: tip_L, mean
      CHARACTER (LEN=strlen) :: fname
      fname = 'EZ_HE_WD.data'
      
      CALL Init_Do_One_Utils ( init_M )
      save_filename = 'demo3.sav'
      stats_N = 0
      stats_X = 0D0
      stats_X2 = 0D0
      
      wd_RGB_name = 'wd_rgb.sav'
      wd_eject_name = 'wd_eject.sav'
      wd_adj_name = 'wd_adj.sav'
      wd_tau_up_name = 'wd_tau_up.sav'
      wd_cool_name = 'wd_cool.sav'
      wd_tau_down_name = 'wd_tau_down.sav'
      wd_add_h_name = 'wd_add_h.sav'
      
      metals_dir = '../metals/z02'
      CALL EZ_Start ( metals_dir )
      
      IF (.FALSE.) THEN ! create model near tip of RGB
         IF (.NOT. EZ_ZAMS(init_M, Initial_Params)) RETURN
         L_SAVE = 2.7D0; OK_TO_SAVE = .FALSE.
         CALL EZ_Evolve(Check_Model_Up_RGB)
         IF ( .NOT. OK_TO_SAVE ) RETURN
         IF ( .NOT. EZ_Save(wd_RGB_name) ) RETURN
         where_saved = model_Number
         WRITE (*,'(/,A,I4)') '  >>>>>>>> Save interim state along giant branch at model number', where_saved
         
         WRITE (*,'(A,/)') '  Now evolve without mass loss until get a helium flash.'
         SAVING_PROFILES = .FALSE.
         CALL EZ_Evolve(Check_Model_To_Stop)
         logged = Log_State(.TRUE.)
         tip_L = log_Luminosity ! luminosity at tip of giant branch is used to decide when to start mass loss
      ELSE
         tip_L = 3.4343
         where_saved = 514
      END IF
      wd_filename = wd_RGB_name
      
      IF (.FALSE.) THEN ! eject envelope with high wind
         IF (.NOT. EZ_Restore(wd_filename)) RETURN
         WRITE (*,'(/,A,I4,A,/)') '  Restore model ', where_saved, ' and evolve with envelope ejection to avoid the tip flash.'
         L_LIMIT1 = tip_L-start_DL; L_LIMIT2 = tip_L-finish_DL; E_LIMIT1 = start_E; E_LIMIT2 = finish_E
         SAVING_PROFILES = .FALSE.; CMI_VALUE = CMI_highwind
         CALL EZ_Evolve(Check_Model_While_Eject)
         TC_LIMIT = 7.8D0; OK_TO_SAVE = .FALSE.
         CALL EZ_Evolve(Check_Model_WD)
         IF ( .NOT. OK_TO_SAVE ) RETURN
         IF ( .NOT. EZ_Save(wd_eject_name) ) RETURN
         logged = Log_State ( .TRUE. )
      END IF
      wd_filename = wd_eject_name
      
      IF (.FALSE.) THEN ! adjust composition of the wd.  remove h and make rest uniform
         IF (.NOT. EZ_Restore(wd_filename)) RETURN
         logged = Log_State (.TRUE.)
          WRITE_PROFILE_TO_RUN = .TRUE.
          IO_RUN_PROFILE = 38
          run_profilename = 'EZ_status.log'
         CALL EZ_Evolve(Remove_H_Check_Model)
         IF ( .NOT. EZ_Save(wd_adj_name) ) RETURN
         logged = Log_State (.TRUE.)
      END IF
      wd_filename = wd_adj_name
      
      IF (.FALSE.) THEN ! increase tau_photosphere
         IF (.NOT. EZ_Restore(wd_filename)) RETURN
         logged = Log_State (.TRUE.)
          WRITE_PROFILE_TO_RUN = .TRUE.
          IO_RUN_PROFILE = 38
          run_profilename = 'EZ_status.log'
          target_tau_Photosphere = 1D2
         CALL EZ_Evolve(Raise_WD_tau_Check_Model)
         IF ( .NOT. EZ_Save(wd_tau_up_name) ) RETURN
         logged = Log_State (.TRUE.)
      END IF
      wd_filename = wd_tau_up_name
      
      ! current state of the code: trying to cool an unadjusted model
      ! instability in abundances (C/N/O) in extreme outer points
      ! look at log mass fraction vs logP plots at logP < 7
      ! attempts to avoid this by adjusting composition / tau have failed
      ! once fix this problem, duplicate Nomoto & Sugimoto (1977), case of mdot = 4e-8
      ! also, need to change to use Lars' eps3alp.f
      
      ! compare cooling rates to those in Althaus & Benvenuto Table 7
      
      IF (.FALSE.) THEN ! cool the wd
         IF (.NOT. EZ_Restore(wd_filename)) RETURN
         logged = Log_State (.TRUE.)
          WRITE_PROFILE_TO_RUN = .TRUE.
          IO_RUN_PROFILE = 38
          run_profilename = 'EZ_status.log'
         SAVING_PROFILES = .FALSE.
         COOLING_TARGET = 6.915D0
         OK_TO_SAVE = .FALSE.
         PSI_limit = 250D6
         CALL EZ_Evolve(Cool_WD_Check_Model)
         IF ( .NOT. OK_TO_SAVE ) RETURN
         OK_TO_SAVE = .FALSE.
         IF ( .NOT. EZ_Save(wd_cool_name) ) RETURN
         logged = Log_State (.TRUE.)
      END IF
      wd_filename = wd_cool_name
      
      IF (.TRUE.) THEN ! restore normal tau_photosphere
         IF (.NOT. EZ_Restore(wd_filename)) RETURN
         logged = Log_State (.TRUE.)
          WRITE_PROFILE_TO_RUN = .TRUE.
          IO_RUN_PROFILE = 38
          run_profilename = 'EZ_status.log'
          target_tau_Photosphere = 2D0/3D0
          tau_Model = 0
         CALL EZ_Evolve(Restore_WD_tau_Check_Model)
         IF ( .NOT. EZ_Save(wd_tau_down_name) ) RETURN
         logged = Log_State (.TRUE.)
      END IF
      wd_filename = wd_tau_down_name
      
      IF (.FALSE.) THEN ! add thin layer of hydrogen.  
         ! raise X to 0.7 in outer 1e-6 Msun (convection will mix)
         IF (.NOT. EZ_Restore(wd_filename)) RETURN
         logged = Log_State (.TRUE.)
          WRITE_PROFILE_TO_RUN = .TRUE.
          IO_RUN_PROFILE = 38
          run_profilename = 'EZ_status.log'
          target_tau_Photosphere = 2D0/3D0
          tau_Model = 0
         CALL EZ_Evolve(Restore_WD_tau_Check_Model)
         IF ( .NOT. EZ_Save(wd_add_h_name) ) RETURN
         logged = Log_State (.TRUE.)
      END IF
      wd_filename = wd_add_h_name
      
      IF (.FALSE.) THEN ! add mass to the wd
         IF (.NOT. EZ_Restore(wd_filename)) RETURN
         logged = Log_State (.TRUE.)
          WRITE_PROFILE_TO_RUN = .TRUE.
          IO_RUN_PROFILE = 38
          run_profilename = 'EZ_status.log'
         SAVING_PROFILES = .FALSE.
         MDOT_TARGET = 4D-8
         PSI_limit = 250D3
         CALL EZ_Evolve(Do_WD_Check_Model)
         logged = Log_State (.TRUE.)
         WRITE(*,*) 'MDOT_TARGET ', MDOT_TARGET
         WRITE(*,*) 'average mdot', (star_Mass - mdot_start_mass)/(star_Age - mdot_start_age)
         WRITE(*,*) 'star_Mass   ', star_Mass
      END IF
      
      WRITE (*,'(/,A,/)') '  End of Demo3'
      
      mean = stats_X/stats_N
      !WRITE (*,*) 'N', stats_N, 'mean', mean, 'sigma', DSQRT(stats_X2/stats_N - mean*mean)
      END SUBROUTINE Demo_3
      
      END MODULE demo3
