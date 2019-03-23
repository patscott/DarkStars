      MODULE star_controls ! This is the control parameter interface to the stellar evolution code
      IMPLICIT NONE
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      ! These are the result codes for Check_Model routines
      INTEGER, PARAMETER :: BACKUP = -1
      INTEGER, PARAMETER :: TERMINATE = 0
      INTEGER, PARAMETER :: KEEP_GOING = 1

      ! The following are parameters for the evolution code
      ! default values are set in the routine Set_Defaults at the end of this module
      ! You can change the values from your Check_Model routine and the system then will use the new settings. 
      
      ! With the REAL valued control parameters, you may find that the system chokes if you change the
      ! value abruptly.  It can be a challenge for the system to adjust to the changes, so you may
      ! need to gradually adjust the value over several successive models.
      
      DOUBLE PRECISION :: tau_Photosphere, tau_Photosphere_SAV
      ! The surface mesh point is placed where the optical depth reaches this value.
      ! A larger value moves the surface mesh point deeper into the star.
      ! The default value is 2/3, which places the "surface" near the bottom of the photosphere.
      ! I've run cases with it equal to 100 when I needed to avoid a problem with an unstable helium envelope.
      DOUBLE PRECISION :: wind_Eta, wind_Eta_SAV
      ! Reimers' wind prefactor, usually called "eta" in the literature
      ! the wind is proportional to R*L/M, so it mainly impacts giants
      DOUBLE PRECISION :: mixing_length_alpha, mixing_length_alpha_SAV
      ! mixing length ratio, alpha
      DOUBLE PRECISION :: overshoot_param_H, overshoot_param_H_SAV
      ! convective overshoot parameter for hydrogen burning shells -- zero implies no overshoot
      DOUBLE PRECISION :: overshoot_param_He, overshoot_param_He_SAV
      ! same as overshoot_param_H, but for helium burning shells -- system uses this if the shell is depleted of hydrogen
      DOUBLE PRECISION :: convective_diffusion_prefactor, convective_diffusion_prefactor_SAV
      ! control for strength of convective diffusion
      DOUBLE PRECISION :: convection_speed_limit, convection_speed_limit_SAV
      ! limit convective velocities to this fraction of the local sound speed
      DOUBLE PRECISION :: grad_star_limit, grad_star_limit_SAV
      ! limit dlnT/dlnP to this (defaults to 1 in order to prevent bogus density inversion in outer envelope)
      
      DOUBLE PRECISION :: extra_Mdot_param, extra_Mdot_param_SAV
      ! in the surface dM/dt calculation, this is multiplied by the total star mass to give an artificial mass loss or gain
      ! the meaning of extra_Mdot_param is the fraction of current star mass added per year, so extra_Mdot_param < 0 for mass loss
      ! for example, this can be used to mimic rapid envelope ejection during the common envelope phase of a binary
      ! use extra_Mdot_param > 0 to mimic mass gain from a companion,
      ! but note that the mass gained will have the composition of this star's surface
      
      DOUBLE PRECISION :: extra_Energy_param, extra_Energy_param_SAV
      ! artificial energy generation rate in ergs/g/sec -- applies uniformly to the entire star
      ! acts like another source term along with the ergs/g/sec from nuclear burning
      ! system will automatically ramp up extra_Energy_param from a small non-zero value by using the following two parameters
      DOUBLE PRECISION :: extra_Energy_max, extra_Energy_max_SAV
      ! saturation amplitude for extra_Energy_param
      ! extra_Energy_param grows exponentially to this value on a timescale given by extra_Energy_time
      DOUBLE PRECISION :: extra_Energy_time, extra_Energy_time_SAV
      ! timescale for extra_Energy_param changes
      ! in effect, dextra_Energy_param/dt = extra_Energy_param*extra_Energy_time*(1-extra_Energy_param/extra_Energy_max)
     
      INTEGER :: ionization_level, ionization_level_SAV
      ! this parameter determines which elements will have their partial ionization fractions computed.
      ! all other elements are considered completely ionized.
      ! value for IONIZATION_LEVEL must be one of the following:
      INTEGER, PARAMETER :: H_ions = 1
      INTEGER, PARAMETER :: H_HE_ions = 2
      INTEGER, PARAMETER :: H_HE_C_ions = 3
      INTEGER, PARAMETER :: H_HE_C_N_ions = 4
      INTEGER, PARAMETER :: H_HE_C_N_O_ions = 5
            
      DOUBLE PRECISION :: core_param_He, core_param_He_SAV
      ! defines core boundary in terms of minimum H1 abundance -- see mass_He_Core for mass inside this boundary
      DOUBLE PRECISION :: core_param_C, core_param_C_SAV
      ! defines core boundary in terms of minimum He4 abundance -- see mass_C_Core for mass inside this boundary
      DOUBLE PRECISION :: core_param_O, core_param_O_SAV
      ! defines core boundary in terms of minimum C12 abundance -- see mass_O_Core for mass inside this boundary
      
      DOUBLE PRECISION :: burn_H, burn_H_SAV
      ! burn_H = 0.0 turns off changes in abundance from proton consuming reactions
      ! burn_H = 1.0 is normal
      DOUBLE PRECISION :: burn_He, burn_He_SAV
      ! like burn_H, but for helium consuming reactions
      DOUBLE PRECISION :: burn_Metals, burn_Metals_SAV
      ! like burn_H, but for metals consuming reactions
      DOUBLE PRECISION :: T_DS_Dt, T_DS_Dt_SAV
      ! prefactor for T*DS/Dt as a power source -- zero means ignore this term in the luminosity equation
      DOUBLE PRECISION :: neutrino_Cooling_Factor, neutrino_Cooling_Factor_SAV
      ! neutrino_Cooling_Factor = 0.0 turns off cooling by neutrinos
      ! neutrino_Cooling_Factor = 1.0 is normal
      DOUBLE PRECISION :: CNO_Factor, CNO_Factor_SAV
      ! CNO_Factor = 1 for normal; CNO_Factor = 0 to turn off CNO reactions.

      DOUBLE PRECISION :: PSI_limit, PSI_limit_SAV
      ! upper limit for degeneracy parameter
      DOUBLE PRECISION :: GAM_limit, GAM_limit_SAV
      ! upper limit for plasma interaction parameter
      
      DOUBLE PRECISION :: timestep_max ! in years
      ! upper limit for timesteps.   value <= 0 means ignore this control.
      
      DOUBLE PRECISION :: timestep_lower_limit, timestep_lower_limit_SAV
      ! lower limit for timestep adjustments -- gives min change factor
      DOUBLE PRECISION :: timestep_upper_limit, timestep_upper_limit_SAV
      ! upper limit for timestep adjustments -- gives max change factor
      DOUBLE PRECISION :: timestep_decrement, timestep_decrement_SAV
      ! before backup, decrease timestep by this factor
      INTEGER :: timestep_hold, timestep_hold_SAV
      ! after backup, hold size of timestep constant for this many timesteps before starting to adjust it again


      INTEGER :: number_of_center_shells_to_mix, number_of_center_shells_to_mix_SAV
      ! keep constant composition for this many shells at center

      INTEGER :: number_of_surface_shells_to_mix, number_of_surface_shells_to_mix_SAV
      ! keep constant composition for this many shells at surface

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      INTEGER, PARAMETER :: CTRL_tau = 1
      INTEGER, PARAMETER :: CTRL_eta = 2
      INTEGER, PARAMETER :: CTRL_alpha = 3
      INTEGER, PARAMETER :: CTRL_over_H = 4
      INTEGER, PARAMETER :: CTRL_over_He = 5
      INTEGER, PARAMETER :: CTRL_conv_diff = 6
      INTEGER, PARAMETER :: CTRL_conv_speed = 7
      INTEGER, PARAMETER :: CTRL_grad_star_limit = 8
      INTEGER, PARAMETER :: CTRL_mdot = 9
      INTEGER, PARAMETER :: CTRL_energy = 10
      INTEGER, PARAMETER :: CTRL_etime = 11
      INTEGER, PARAMETER :: CTRL_emax = 12
      INTEGER, PARAMETER :: CTRL_ionization = 13
      INTEGER, PARAMETER :: CTRL_core_He = 14
      INTEGER, PARAMETER :: CTRL_core_C = 15
      INTEGER, PARAMETER :: CTRL_core_O = 16
      INTEGER, PARAMETER :: CTRL_burn_H = 17
      INTEGER, PARAMETER :: CTRL_burn_He = 18
      INTEGER, PARAMETER :: CTRL_burn_Metals = 19
      INTEGER, PARAMETER :: CTRL_T_DS_Dt = 20
      INTEGER, PARAMETER :: CTRL_cool = 21
      INTEGER, PARAMETER :: CTRL_CNO_Factor = 22
      INTEGER, PARAMETER :: CTRL_PSI_limit = 23
      INTEGER, PARAMETER :: CTRL_GAM_limit = 24
      INTEGER, PARAMETER :: CTRL_timestep_lower_limit = 25
      INTEGER, PARAMETER :: CTRL_timestep_upper_limit = 26
      INTEGER, PARAMETER :: CTRL_timestep_decrement = 27
      INTEGER, PARAMETER :: CTRL_timestep_hold = 28
      INTEGER, PARAMETER :: CTRL_center_mix = 29
      INTEGER, PARAMETER :: CTRL_surface_mix = 30
      
      INTEGER, PARAMETER :: num_CTRLs = 30
      
      CONTAINS
      
      SUBROUTINE Save_Control_Values(a)
      DOUBLE PRECISION, INTENT(OUT) :: a(num_CTRLs)
      a(CTRL_tau) = tau_Photosphere
      a(CTRL_eta) = wind_Eta
      a(CTRL_alpha) = mixing_length_alpha
      a(CTRL_over_H) = overshoot_param_H
      a(CTRL_over_He) = overshoot_param_He
      a(CTRL_conv_diff) = convective_diffusion_prefactor
      a(CTRL_conv_speed) = convection_speed_limit
      a(CTRL_grad_star_limit) = grad_star_limit
      a(CTRL_mdot) = extra_Mdot_param
      a(CTRL_energy) = extra_Energy_param
      a(CTRL_emax) = extra_Energy_max
      a(CTRL_etime) = extra_Energy_time
      a(CTRL_ionization) = ionization_level
      a(CTRL_core_He) = core_param_He
      a(CTRL_core_C) = core_param_C
      a(CTRL_core_O) = core_param_O
      a(CTRL_burn_H) = burn_H
      a(CTRL_burn_He) = burn_He
      a(CTRL_burn_Metals) = burn_Metals
      a(CTRL_T_DS_Dt) = T_DS_Dt
      a(CTRL_cool) = neutrino_Cooling_Factor
      a(CTRL_CNO_Factor) = CNO_Factor
      a(CTRL_PSI_limit) = PSI_limit
      a(CTRL_GAM_limit) = GAM_limit
      a(CTRL_timestep_lower_limit) = timestep_lower_limit
      a(CTRL_timestep_upper_limit) = timestep_upper_limit
      a(CTRL_timestep_decrement) = timestep_decrement
      a(CTRL_timestep_hold) = timestep_hold
      a(CTRL_center_mix) = number_of_center_shells_to_mix
      a(CTRL_surface_mix) = number_of_surface_shells_to_mix
      END SUBROUTINE Save_Control_Values
      
      SUBROUTINE Restore_Control_Values(a)
      DOUBLE PRECISION, INTENT(IN) :: a(num_CTRLs)
      tau_Photosphere = a(CTRL_tau)
      wind_Eta = a(CTRL_eta)
      mixing_length_alpha = a(CTRL_alpha)
      overshoot_param_H = a(CTRL_over_H)
      overshoot_param_He = a(CTRL_over_He)
      convective_diffusion_prefactor = a(CTRL_conv_diff)
      convection_speed_limit = a(CTRL_conv_speed)
      grad_star_limit = a(CTRL_grad_star_limit)
      extra_Mdot_param = a(CTRL_mdot) 
      extra_Energy_param = a(CTRL_energy)
      extra_Energy_max = a(CTRL_emax)
      extra_Energy_time = a(CTRL_etime)
      ionization_level = a(CTRL_ionization)
      core_param_He = a(CTRL_core_He)
      core_param_C = a(CTRL_core_C)
      core_param_O = a(CTRL_core_O)
      burn_H = a(CTRL_burn_H)
      burn_He = a(CTRL_burn_He)
      burn_Metals = a(CTRL_burn_Metals)
      T_DS_Dt = a(CTRL_T_DS_Dt)
      neutrino_Cooling_Factor = a(CTRL_cool)
      CNO_Factor = a(CTRL_CNO_Factor)
      PSI_limit = a(CTRL_PSI_limit)
      GAM_limit = a(CTRL_GAM_limit)
      timestep_lower_limit = a(CTRL_timestep_lower_limit)
      timestep_upper_limit = a(CTRL_timestep_upper_limit)
      timestep_decrement = a(CTRL_timestep_decrement)
      timestep_hold = a(CTRL_timestep_hold)
      number_of_center_shells_to_mix = a(CTRL_center_mix)
      number_of_surface_shells_to_mix = a(CTRL_surface_mix)
      END SUBROUTINE Restore_Control_Values
      
      SUBROUTINE Set_Defaults
      ! this gets called whenever you call EZ_ZAMS, so the defaults get reset each time.
      ! if you want to change some default value, you should do it in your Initialize_Params routine that gets passed to EZ_ZAMS.
      ! in addition, your Check_Model routine can change any of these parameters when it is called between evolution steps.
      tau_Photosphere = 2D0/3D0
      wind_Eta = 0D0
      mixing_length_alpha = 2D0
      overshoot_param_H = 0.12D0
      overshoot_param_He = 0.12D0
      
      convective_diffusion_prefactor = 1D-4
      convection_speed_limit = 1D6 ! turn this limit off by default.  allow supersonic speeds in convection.
      grad_star_limit = 1D6 ! turn this limit off by default and allow density inversion in superadiabatic regions
      extra_Mdot_param = 0D0
      extra_Energy_param = 0D0
      extra_Energy_max = 1D2
      extra_Energy_time = 1D-6
      ionization_level = H_HE_ions
      core_param_He = .15D0
      core_param_C = .25D0
      core_param_O = .35D0
      burn_H = 1
      burn_He = 1
      burn_Metals = 1
      T_DS_Dt = 1
      neutrino_Cooling_Factor = 1D0
      CNO_Factor = 1d0
      PSI_limit = 1000d0
      GAM_limit = 150d0
      timestep_lower_limit = 0.8D0
      timestep_upper_limit = 1.2D0
      timestep_decrement = 0.5D0
      timestep_hold = 3
      number_of_center_shells_to_mix = 15
      number_of_surface_shells_to_mix = 40
      END SUBROUTINE Set_Defaults

      SUBROUTINE SAV_star_controls( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) timestep_hold, tau_Photosphere, wind_Eta, mixing_length_alpha, overshoot_param_H, overshoot_param_He
         READ (IO_UNIT) convective_diffusion_prefactor, convection_speed_limit, grad_star_limit
         READ (IO_UNIT) extra_Mdot_param, extra_Energy_param, extra_Energy_time, extra_Energy_max, ionization_level
         READ (IO_UNIT) core_param_He, core_param_C, core_param_O
         READ (IO_UNIT) burn_H, burn_He, burn_Metals, T_DS_Dt, neutrino_Cooling_Factor
         READ (IO_UNIT) CNO_Factor, PSI_limit, GAM_limit, timestep_lower_limit, timestep_upper_limit, timestep_decrement
         READ (IO_UNIT) number_of_center_shells_to_mix, number_of_surface_shells_to_mix
      ELSE
         WRITE (IO_UNIT) timestep_hold, tau_Photosphere, wind_Eta, mixing_length_alpha, overshoot_param_H, overshoot_param_He
         WRITE (IO_UNIT) convective_diffusion_prefactor, convection_speed_limit, grad_star_limit
         WRITE (IO_UNIT) extra_Mdot_param, extra_Energy_param, extra_Energy_time, extra_Energy_max, ionization_level
         WRITE (IO_UNIT) core_param_He, core_param_C, core_param_O
         WRITE (IO_UNIT) burn_H, burn_He, burn_Metals, T_DS_Dt, neutrino_Cooling_Factor
         WRITE (IO_UNIT) CNO_Factor, PSI_limit, GAM_limit, timestep_lower_limit, timestep_upper_limit, timestep_decrement
         WRITE (IO_UNIT) number_of_center_shells_to_mix, number_of_surface_shells_to_mix
      END IF
      END SUBROUTINE SAV_star_controls

      SUBROUTINE SAV_star_controls_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         timestep_hold=timestep_hold_SAV; tau_Photosphere=tau_Photosphere_SAV; wind_Eta=wind_Eta_SAV
		 mixing_length_alpha=mixing_length_alpha_SAV; overshoot_param_H=overshoot_param_H_SAV
		 overshoot_param_He=overshoot_param_He_SAV; convective_diffusion_prefactor=convective_diffusion_prefactor_SAV
		 convection_speed_limit=convection_speed_limit_SAV; grad_star_limit=grad_star_limit_SAV
         extra_Mdot_param=extra_Mdot_param_SAV; extra_Energy_param=extra_Energy_param_SAV
		 extra_Energy_time=extra_Energy_time_SAV; extra_Energy_max=extra_Energy_max_SAV
		 ionization_level=ionization_level_SAV; core_param_He=core_param_He_SAV; core_param_C=core_param_C_SAV
		 core_param_O=core_param_O_SAV; burn_H=burn_H_SAV; burn_He=burn_He_SAV; burn_Metals=burn_Metals_SAV
		 T_DS_Dt=T_DS_Dt_SAV; neutrino_Cooling_Factor=neutrino_Cooling_Factor_SAV; CNO_Factor=CNO_Factor_SAV
		 PSI_limit=PSI_limit_SAV; GAM_limit=GAM_limit_SAV; timestep_lower_limit=timestep_lower_limit_SAV
		 timestep_upper_limit=timestep_upper_limit_SAV; timestep_decrement=timestep_decrement_SAV
         number_of_center_shells_to_mix=number_of_center_shells_to_mix_SAV
		 number_of_surface_shells_to_mix=number_of_surface_shells_to_mix_SAV
      ELSE
         timestep_hold_SAV=timestep_hold; tau_Photosphere_SAV=tau_Photosphere; wind_Eta_SAV=wind_Eta
		 mixing_length_alpha_SAV=mixing_length_alpha; overshoot_param_H_SAV=overshoot_param_H
		 overshoot_param_He_SAV=overshoot_param_He; convective_diffusion_prefactor_SAV=convective_diffusion_prefactor
		 convection_speed_limit_SAV=convection_speed_limit; grad_star_limit_SAV=grad_star_limit
         extra_Mdot_param_SAV=extra_Mdot_param; extra_Energy_param_SAV=extra_Energy_param
		 extra_Energy_time_SAV=extra_Energy_time; extra_Energy_max_SAV=extra_Energy_max
		 ionization_level_SAV=ionization_level; core_param_He_SAV=core_param_He; core_param_C_SAV=core_param_C
		 core_param_O_SAV=core_param_O; burn_H_SAV=burn_H; burn_He_SAV=burn_He; burn_Metals_SAV=burn_Metals
		 T_DS_Dt_SAV=T_DS_Dt; neutrino_Cooling_Factor_SAV=neutrino_Cooling_Factor; CNO_Factor_SAV=CNO_Factor
		 PSI_limit_SAV=PSI_limit; GAM_limit_SAV=GAM_limit; timestep_lower_limit_SAV=timestep_lower_limit
		 timestep_upper_limit_SAV=timestep_upper_limit; timestep_decrement_SAV=timestep_decrement
         number_of_center_shells_to_mix_SAV=number_of_center_shells_to_mix
		 number_of_surface_shells_to_mix_SAV=number_of_surface_shells_to_mix
      END IF
      END SUBROUTINE SAV_star_controls_internal

      END MODULE star_controls

