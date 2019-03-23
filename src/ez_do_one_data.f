      MODULE ez_do_one_data
      USE star_data
      USE star_extras
      USE star_controls
      USE star_constants
      IMPLICIT NONE
      
      LOGICAL :: WRITE_DIFFME, WRITE_DIFFME_SAV ! default is false. if true, then writes info to IO_DIFFME file
      INTEGER, PARAMETER :: END_OF_DIFFME = 75 ! stop writing after this model number
      
      LOGICAL :: WRITE_PROFILE_TO_RUN, WRITE_PROFILE_TO_RUN_SAV ! write profile info to the run directory
      CHARACTER (LEN=strlen) :: run_profilename, run_profilename_SAV

      LOGICAL :: WRITE_BRIEF_TO_RUN ! write brief version of profile info to the run directory
      CHARACTER (LEN=strlen) :: run_briefname

      LOGICAL :: WRITE_PROFILE_TO_DATA, WRITE_PROFILE_TO_DATA_SAV ! write profile info to the data directory
      CHARACTER (LEN=strlen) :: data_profilename, data_profilename_SAV
      CHARACTER (LEN=strlen) :: data_profile_prefix, data_profile_prefix_SAV
      
      LOGICAL :: WRITE_PROFILE_TO_MODELS ! write profile info to the models directory
      CHARACTER (LEN=strlen) :: model_profilename
      CHARACTER (LEN=strlen) :: model_profile_prefix
      
      INTEGER :: IO_DIFFME, IO_RUN_PROFILE, IO_DATA_PROFILE, IO_MODEL_PROFILE
      INTEGER :: IO_DIFFME_SAV, IO_RUN_PROFILE_SAV, IO_DATA_PROFILE_SAV
      
      INTEGER :: recent_profile_model ! model_Number for the most recently saved profile

      INTEGER :: profile_model, profile_model_SAV ! save profiles if model_Number equals this
      INTEGER :: stop_at_Model, stop_at_Model_SAV ! stop the run when reach this model number
      INTEGER :: tau_Model, tau_Model_SAV
      INTEGER :: number_of_little_steps, number_of_little_steps_SAV
      LOGICAL :: save_helium_ignition, save_helium_ignition_SAV
      LOGICAL :: helium_ignition, helium_ignition_SAV
      LOGICAL :: do_post_He_flash, do_post_He_flash_SAV
      LOGICAL :: do_CNTR_RHOs, do_CNTR_RHOs_SAV
      LOGICAL :: do_CNTR_T_drops, do_CNTR_T_drops_SAV
      DOUBLE PRECISION :: target_tau_Photosphere, target_tau_Photosphere_SAV
      DOUBLE PRECISION :: prev_TCNTR1, prev_TCNTR1_SAV
      DOUBLE PRECISION :: prev_AGE1, prev_AGE1_SAV
      DOUBLE PRECISION :: max_AGE, max_AGE_SAV
      DOUBLE PRECISION :: prev_TCNTR2, prev_TCNTR2_SAV
      DOUBLE PRECISION :: prev_AGE2, prev_AGE2_SAV
      DOUBLE PRECISION :: prev_TSURF, prev_TSURF_SAV
      DOUBLE PRECISION :: profile_AGE, profile_AGE_SAV
      DOUBLE PRECISION :: post_He_AGE, post_He_AGE_SAV
      DOUBLE PRECISION :: prev_Luminosity, prev_Luminosity_SAV
      DOUBLE PRECISION :: ignition_Center_XHE, ignition_Center_XHE_SAV
      DOUBLE PRECISION :: he_Luminosity_Limit, he_Luminosity_Limit_SAV
      DOUBLE PRECISION :: prev_CNTR_RHO, prev_CNTR_RHO_SAV
      DOUBLE PRECISION :: next_CNTR_RHO, next_CNTR_RHO_SAV
      DOUBLE PRECISION, PARAMETER :: del_CNTR_RHO = 1D0
      DOUBLE PRECISION, PARAMETER :: min_CNTR_RHO = 3D0
      DOUBLE PRECISION, PARAMETER :: lower_He_flash = 0.79D0
      DOUBLE PRECISION, PARAMETER :: upper_He_flash = 2.01D0
      DOUBLE PRECISION, PARAMETER :: min_Mass = 0.3D0
      DOUBLE PRECISION, PARAMETER :: no_He_ignition_limit = 0.75D0
      DOUBLE PRECISION, PARAMETER :: no_CNTR_T_drops_limit = 6.5D0
      CHARACTER (LEN=strlen) :: save_filename, save_filename_SAV
      CHARACTER (LEN=strlen) :: metals_dir, metals_dir_SAV
      CHARACTER (LEN=strlen) :: data_dir, data_dir_SAV
      LOGICAL :: call_save_exp_info
      INTEGER :: summary_cnt, summary_cnt_SAV
      INTEGER :: head_cnt, head_cnt_SAV

      INTEGER :: phase_of_evolution, phase_of_evolution_SAV
      INTEGER, PARAMETER :: phase_STARTING = 0
      INTEGER, PARAMETER :: phase_END_ZAMS = 1
      INTEGER, PARAMETER :: phase_WAIT_FOR_HE = 2
      INTEGER, PARAMETER :: phase_HE_IGNITION_OVER = 3
      INTEGER, PARAMETER :: phase_GET_POST_FLASH = 4
      INTEGER, PARAMETER :: phase_HE_IGNITING = 5
      INTEGER, PARAMETER :: phase_FINISHING = 6
        
      double precision :: del_log_Luminosity, del_log_surface_Temp, del_log_center_Temp, del_log_center_Density
      double precision :: prv_log_Luminosity, prv_log_surface_Temp, prv_log_center_Temp, prv_log_center_Density
      double precision :: del_log_Luminosity_SAV, del_log_surface_Temp_SAV, del_log_center_Temp_SAV, del_log_center_Density_SAV
      double precision :: prv_log_Luminosity_SAV, prv_log_surface_Temp_SAV, prv_log_center_Temp_SAV, prv_log_center_Density_SAV

      CONTAINS
      
      SUBROUTINE SAV_do_one_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) WRITE_DIFFME, WRITE_PROFILE_TO_RUN, WRITE_PROFILE_TO_DATA
         READ (IO_UNIT) run_profilename, data_profilename
         READ (IO_UNIT) data_profile_prefix
         READ (IO_UNIT) IO_DIFFME, IO_RUN_PROFILE, IO_DATA_PROFILE
         READ (IO_UNIT) profile_model, stop_at_Model, number_of_little_steps, tau_Model
         READ (IO_UNIT) save_helium_ignition, helium_ignition, do_post_He_flash, do_CNTR_RHOs, do_CNTR_T_drops
         READ (IO_UNIT) prev_TCNTR1, prev_AGE1, max_AGE, target_tau_Photosphere
         READ (IO_UNIT) prev_TCNTR2, prev_AGE2, prev_TSURF, profile_AGE, post_He_AGE
         READ (IO_UNIT) prev_Luminosity, ignition_Center_XHE, he_Luminosity_Limit, prev_CNTR_RHO, next_CNTR_RHO
         READ (IO_UNIT) phase_of_evolution
         READ (IO_UNIT) save_filename, metals_dir, data_dir
         READ (IO_UNIT) summary_cnt, head_cnt
         READ (IO_UNIT) del_log_Luminosity, del_log_surface_Temp, del_log_center_Temp, del_log_center_Density
         READ (IO_UNIT) prv_log_Luminosity, prv_log_surface_Temp, prv_log_center_Temp, prv_log_center_Density
      ELSE
         WRITE (IO_UNIT) WRITE_DIFFME, WRITE_PROFILE_TO_RUN, WRITE_PROFILE_TO_DATA
         WRITE (IO_UNIT) run_profilename, data_profilename
         WRITE (IO_UNIT) data_profile_prefix
         WRITE (IO_UNIT) IO_DIFFME, IO_RUN_PROFILE, IO_DATA_PROFILE
         WRITE (IO_UNIT) profile_model, stop_at_Model, number_of_little_steps, tau_Model
         WRITE (IO_UNIT) save_helium_ignition, helium_ignition, do_post_He_flash, do_CNTR_RHOs, do_CNTR_T_drops
         WRITE (IO_UNIT) prev_TCNTR1, prev_AGE1, max_AGE, target_tau_Photosphere
         WRITE (IO_UNIT) prev_TCNTR2, prev_AGE2, prev_TSURF, profile_AGE, post_He_AGE
         WRITE (IO_UNIT) prev_Luminosity, ignition_Center_XHE, he_Luminosity_Limit, prev_CNTR_RHO, next_CNTR_RHO
         WRITE (IO_UNIT) phase_of_evolution
         WRITE (IO_UNIT) save_filename, metals_dir, data_dir
         WRITE (IO_UNIT) summary_cnt, head_cnt
         WRITE (IO_UNIT) del_log_Luminosity, del_log_surface_Temp, del_log_center_Temp, del_log_center_Density
         WRITE (IO_UNIT) prv_log_Luminosity, prv_log_surface_Temp, prv_log_center_Temp, prv_log_center_Density
      END IF
      END SUBROUTINE SAV_do_one_data
      
      SUBROUTINE SAV_do_one_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         WRITE_DIFFME=WRITE_DIFFME_SAV; WRITE_PROFILE_TO_RUN=WRITE_PROFILE_TO_RUN_SAV
		 WRITE_PROFILE_TO_DATA=WRITE_PROFILE_TO_DATA_SAV; run_profilename=run_profilename_SAV
		 data_profilename=data_profilename_SAV; data_profile_prefix=data_profile_prefix_SAV
         IO_DIFFME=IO_DIFFME_SAV; IO_RUN_PROFILE=IO_RUN_PROFILE_SAV; IO_DATA_PROFILE=IO_DATA_PROFILE_SAV
         profile_model=profile_model_SAV; stop_at_Model=stop_at_Model_SAV
		 number_of_little_steps=number_of_little_steps_SAV; tau_Model=tau_Model_SAV
         save_helium_ignition=save_helium_ignition_SAV; helium_ignition=helium_ignition_SAV
		 do_post_He_flash=do_post_He_flash_SAV; do_CNTR_RHOs=do_CNTR_RHOs_SAV; do_CNTR_T_drops=do_CNTR_T_drops_SAV
         prev_TCNTR1=prev_TCNTR1_SAV; prev_AGE1=prev_AGE1_SAV; max_AGE=max_AGE_SAV
		 target_tau_Photosphere=target_tau_Photosphere_SAV; prev_TCNTR2=prev_TCNTR2_SAV; prev_AGE2=prev_AGE2_SAV
		 prev_TSURF=prev_TSURF_SAV; profile_AGE=profile_AGE_SAV; post_He_AGE=post_He_AGE_SAV
         prev_Luminosity=prev_Luminosity_SAV; ignition_Center_XHE=ignition_Center_XHE_SAV
		 he_Luminosity_Limit=he_Luminosity_Limit_SAV; prev_CNTR_RHO=prev_CNTR_RHO_SAV; next_CNTR_RHO=next_CNTR_RHO_SAV
         phase_of_evolution=phase_of_evolution_SAV; save_filename=save_filename_SAV; metals_dir=metals_dir_SAV
		 data_dir=data_dir_SAV; summary_cnt=summary_cnt_SAV; head_cnt=head_cnt_SAV
         del_log_Luminosity=del_log_Luminosity_SAV; del_log_surface_Temp=del_log_surface_Temp_SAV
		 del_log_center_Temp=del_log_center_Temp_SAV; del_log_center_Density=del_log_center_Density_SAV
         prv_log_Luminosity=prv_log_Luminosity_SAV; prv_log_surface_Temp=prv_log_surface_Temp_SAV
		 prv_log_center_Temp=prv_log_center_Temp_SAV; prv_log_center_Density=prv_log_center_Density_SAV
      ELSE
         WRITE_DIFFME_SAV=WRITE_DIFFME; WRITE_PROFILE_TO_RUN_SAV=WRITE_PROFILE_TO_RUN
		 WRITE_PROFILE_TO_DATA_SAV=WRITE_PROFILE_TO_DATA; run_profilename_SAV=run_profilename
		 data_profilename_SAV=data_profilename; data_profile_prefix_SAV=data_profile_prefix
         IO_DIFFME_SAV=IO_DIFFME; IO_RUN_PROFILE_SAV=IO_RUN_PROFILE; IO_DATA_PROFILE_SAV=IO_DATA_PROFILE
         profile_model_SAV=profile_model; stop_at_Model_SAV=stop_at_Model
		 number_of_little_steps_SAV=number_of_little_steps; tau_Model_SAV=tau_Model
         save_helium_ignition_SAV=save_helium_ignition; helium_ignition_SAV=helium_ignition
		 do_post_He_flash_SAV=do_post_He_flash; do_CNTR_RHOs_SAV=do_CNTR_RHOs; do_CNTR_T_drops_SAV=do_CNTR_T_drops
         prev_TCNTR1_SAV=prev_TCNTR1; prev_AGE1_SAV=prev_AGE1; max_AGE_SAV=max_AGE
		 target_tau_Photosphere_SAV=target_tau_Photosphere; prev_TCNTR2_SAV=prev_TCNTR2; prev_AGE2_SAV=prev_AGE2
		 prev_TSURF_SAV=prev_TSURF; profile_AGE_SAV=profile_AGE; post_He_AGE_SAV=post_He_AGE
         prev_Luminosity_SAV=prev_Luminosity; ignition_Center_XHE_SAV=ignition_Center_XHE
		 he_Luminosity_Limit_SAV=he_Luminosity_Limit; prev_CNTR_RHO_SAV=prev_CNTR_RHO; next_CNTR_RHO_SAV=next_CNTR_RHO
         phase_of_evolution_SAV=phase_of_evolution; save_filename_SAV=save_filename; metals_dir_SAV=metals_dir
		 data_dir_SAV=data_dir; summary_cnt_SAV=summary_cnt; head_cnt_SAV=head_cnt
         del_log_Luminosity_SAV=del_log_Luminosity; del_log_surface_Temp_SAV=del_log_surface_Temp
		 del_log_center_Temp_SAV=del_log_center_Temp; del_log_center_Density_SAV=del_log_center_Density
         prv_log_Luminosity_SAV=prv_log_Luminosity; prv_log_surface_Temp_SAV=prv_log_surface_Temp
		 prv_log_center_Temp_SAV=prv_log_center_Temp; prv_log_center_Density_SAV=prv_log_center_Density
	  END IF
      END SUBROUTINE SAV_do_one_data_internal

      END MODULE ez_do_one_data
      
