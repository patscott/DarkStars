      MODULE ez_do_one_utils
      USE star_data
      USE star_extras
      USE star_controls
      USE star_constants
      USE ez_do_one_data
      IMPLICIT NONE
      
      CONTAINS
      
      SUBROUTINE Write_DIFFME_File(io_diff, num_out, out)
        INTEGER, INTENT(IN) :: io_diff, num_out
        DOUBLE PRECISION, INTENT(IN) :: out(num_out)
        INTEGER :: j, diff_out
        DOUBLE PRECISION :: val, aval, vlog
        DO j=1,num_out
            val = out(j); aval = DABS(val);
            IF ( aval .LT. 1D-30 ) THEN
               WRITE(io_diff,'(I9)') 0
            ELSE
               vlog = LOG10(aval); aval = aval * 10D0**(3-floor(vlog+1D-3))
               diff_out = floor(aval + 0.5D0)
               IF ( val .LT. 0D0 ) diff_out = -diff_out
               WRITE(io_diff,'(I9)') diff_out
            END IF
        END DO
      END SUBROUTINE Write_DIFFME_File
      
      SUBROUTINE Init_Do_One_Utils( init_M )
      DOUBLE PRECISION, INTENT(IN) :: init_M

      max_AGE = 1D15
      stop_at_Model = -1
      profile_model = -1
      recent_profile_model = -1
      save_filename = 'ez.sav'
      save_helium_ignition = .TRUE.
      WRITE_DIFFME = .FALSE.
      WRITE_PROFILE_TO_RUN = .FALSE.
      WRITE_PROFILE_TO_MODELS = .FALSE.
      WRITE_BRIEF_TO_RUN = .FALSE.
      WRITE_PROFILE_TO_DATA = .FALSE.
      call_save_exp_info = .FALSE.
      tau_Model = 0
      target_tau_Photosphere = -1D0 ! not changing tau_Photosphere, so set to negative value to indicate that.
      do_post_He_flash = ( init_M .GE. lower_He_flash .AND. init_M .LE. upper_He_flash )
      do_CNTR_RHOs = init_M .LE. no_He_ignition_limit
      do_CNTR_T_drops = init_M .LE. no_CNTR_T_drops_limit
      phase_of_evolution = phase_STARTING
      post_He_AGE = -1D0
      profile_AGE = -1D0
      prev_CNTR_RHO = 1D99
      helium_ignition = .FALSE.
      number_of_little_steps = 0
      head_cnt = 200; summary_cnt = 20
      del_log_Luminosity = 1d6
      del_log_surface_Temp = 1d6
      del_log_center_Temp = 1d6
      del_log_center_Density = 1d6
      
      END SUBROUTINE Init_Do_One_Utils
      
      SUBROUTINE Write_Log_Header
      CHARACTER (LEN=300) :: log_str
      WRITE (*,*)
      log_str = '  MOD#  log RHO   log Tc   log Pc   log L    log Ts  log Age  log DT  log L_He   log R '
      log_str = trim(log_str) // '  Ttl_Mass  He_Core  Noncore     XH1     XHe4     XC12     XN14     XO16   PSIcntr  GAMcntr'
      WRITE (*,'(a)') trim(log_str)
      END SUBROUTINE Write_Log_Header
      
      LOGICAL FUNCTION Log_State ( do_write )
      ! Generates output to the terminal giving a brief summary of the current model.
      USE ez_driver
      LOGICAL, INTENT(IN) :: do_write
      INTEGER, PARAMETER :: out_size = 19
      DOUBLE PRECISION :: out(out_size), age
      INTEGER :: model
      age = star_Age
      model = model_Number
      Log_State = .FALSE.
      IF ( model .EQ. 0 ) THEN
         CALL Write_Log_Header
         RETURN
      END IF
      IF ( (.NOT. do_write) .AND. (model .NE. 1) .AND. (MOD(model,summary_cnt) .NE. 0) ) RETURN
      IF ( MOD(model,head_cnt) .EQ. 0 ) CALL Write_Log_Header
      CALL EZ_Extras
      out(1) = log_center_Density; out(2) = log_center_Temp; out(3) = log_center_Pressure; out(4) = log_Luminosity;
      out(5) = log_surface_Temp; out(6) = LOG10(1D-99+max(0D0,age))
      out(7) = LOG10(time_Step); out(8) = LOG10(1D-99+power_He_burn); out(9) = log_Radius
      out(10) = star_Mass; out(11) = mass_He_Core; out(12) = star_Mass-mass_He_Core
      out(13) = center_H; out(14) = center_He; out(15) = center_C;
      out(16) = center_N; out(17) = center_O; out(18) = SX(SX_PSI,SX_CNTR)
      out(19) = SX(SX_GAM,SX_CNTR)
      WRITE(*,'(I6,22(F9.4))') model, out
      IF ( phase_of_evolution .NE. phase_GET_POST_FLASH ) THEN
         IF ( WRITE_DIFFME ) THEN
            CALL Write_DIFFME_File(IO_DIFFME, out_size, out)
         END IF
         IF ( WRITE_PROFILE_TO_RUN ) THEN
            CALL Write_Status_Info(IO_RUN_PROFILE, run_profilename, model)
         END IF
         IF ( WRITE_PROFILE_TO_MODELS ) THEN
            CALL Setup_Profilename(model_profile_prefix, model_profilename, model)
            CALL Write_Status_Info(IO_MODEL_PROFILE, model_profilename, model)
         END IF
         IF ( WRITE_BRIEF_TO_RUN ) THEN
            CALL Write_Brief_Info(IO_RUN_PROFILE, run_briefname, model)
         END IF
      END IF
      Log_State = .TRUE.
      END FUNCTION Log_State
      
      SUBROUTINE Setup_Profilename(prefix, fullname, model)
      CHARACTER (LEN=strlen), INTENT(IN) :: prefix
      CHARACTER (LEN=strlen), INTENT(OUT) :: fullname
      INTEGER, INTENT(IN) :: model
      IF ( model_Number .LT. 10 ) THEN
         WRITE(fullname,'(A,I1,A)') trim(prefix), model_Number, '.log'
      ELSEIF ( model_Number .LT. 100 ) THEN
         WRITE(fullname,'(A,I2,A)') trim(prefix), model_Number, '.log'
      ELSEIF ( model_Number .LT. 1000 ) THEN
         WRITE(fullname,'(A,I3,A)') trim(prefix), model_Number, '.log'
      ELSEIF ( model_Number .LT. 10000 ) THEN
         WRITE(fullname,'(A,I4,A)') trim(prefix), model_Number, '.log'
      ELSE
         WRITE(fullname,'(A,I5,A)') trim(prefix), model_Number, '.log'
      END IF
      END SUBROUTINE Setup_Profilename
      
      SUBROUTINE Write_Brief_Info (io_unit, fname, model)
      INTEGER, INTENT(IN) :: io_unit, model
      CHARACTER (LEN=strlen) :: fname
      INTEGER :: I, ios
      ios = 0
      OPEN(UNIT=io_unit, FILE=trim(fname), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         RETURN
      END IF
      DO I = SX_SURF, SX_CNTR
         WRITE(io_unit,*) SX(SX_M,I), SX(SX_R,I), SX(SX_RHO,I), SX(SX_T,I), SX(SX_XH,I), SX(SX_XHE,I), SX(SX_P,I)
      END DO
      CLOSE(io_unit)
      END SUBROUTINE Write_Brief_Info
      
      SUBROUTINE Write_Status_Info (io_unit, fname, model)
      ! Added for DarkStars - Pat Scott 20070905
      USE DkStrs_data
      USE DkStrs_utils
      USE DkStrs_WIMPdens
      USE DkStrs_annihilation
      USE DkStrs_transport
      DOUBLE PRECISION :: WIMPdensSave(meshpoints-1), mfpSave(meshpoints-1)
      DOUBLE PRECISION :: E_transportSave(meshpoints-1), E_annihilationSave(meshpoints-1) 
      ! End DarkStars addition
      INTEGER, INTENT(IN) :: io_unit, model
      CHARACTER (LEN=strlen) :: fname
      INTEGER :: I, ios
      DOUBLE PRECISION :: pNe
      ios = 0
      OPEN(UNIT=io_unit, FILE=trim(fname), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         RETURN
      END IF
      WRITE(io_unit,*) model_Number, star_Age, time_Step
      WRITE(io_unit,'(99(E18.12,1x))') log_Luminosity, log_Radius, log_surface_Temp, log_center_Temp, log_center_Density
      WRITE(io_unit,'(99(E18.12,1x))') log_center_Pressure, center_Degeneracy, center_H, center_He
      WRITE(io_unit,'(99(E18.12,1x))') center_C, center_N, center_O, center_Ne
      WRITE(io_unit,'(99(E18.12,1x))') center_Mg, center_Si, center_Fe
      WRITE(io_unit,'(99(E18.12,1x))') initial_Mass, star_Mass, star_Mdot, initial_Z
      WRITE(io_unit,'(99(E18.12,1x))') star_Mass_H, star_Mass_He, star_Mass_C, star_Mass_N
      WRITE(io_unit,'(99(E18.12,1x))') star_Mass_O, star_Mass_Ne, star_Mass_Mg, star_Mass_Si, star_Mass_Fe
      WRITE(io_unit,'(99(E18.12,1x))') mass_He_Core, mass_C_Core, mass_O_Core
      WRITE(io_unit,'(99(E18.12,1x))') dynamic_Timescale, KH_Timescale,  nuc_Timescale
      WRITE(io_unit,'(99(E18.12,1x))') power_H_burn, power_He_burn, power_Metal_burn, power_Neutrinos
      pNe = power_Ne_alpha
      WRITE(io_unit,'(99(E18.12,1x))') power_PP, power_CNO, power_3_alpha, power_C_alpha, power_N_alpha, power_O_alpha, pNe
      WRITE(io_unit,'(99(E18.12,1x))') power_CC_Ne, power_CO, power_OO, power_Ne_decay, power_Mg_decay, power_CC_Mg
      WRITE(io_unit,'(99(E18.12,1x))') power_plasmon_neutrinos, power_brem_neutrinos, power_pair_neutrinos, power_photo_neutrinos
      DO I = 1, SX_NUM ! write the entire SX array
         WRITE(io_unit,'(999(E18.12,1x))') SX(I,SX_SURF:SX_CNTR)
      END DO
      !Added facility for writing WIMPy things to profile files - Pat Scott 20070905
      DO I = meshpoints, 2, -1 
         WIMPdensSave(meshpoints-I+1) = WIMPdens(starr(I)*r_star*1.d2)
         mfpSave(meshpoints-I+1) = mfp(I)
         E_transportSave(meshpoints-I+1) = E_transport(starr(I)*r_star*1.d2,SX(SX_RHO,SX_CNTR-SX_SURF-I+3))
         E_annihilationSave(meshpoints-I+1) = 2.d0 * mx * 1.d3 * CMEV * sigann / 
     &    SX(SX_RHO,SX_CNTR-SX_SURF-I+3) * WIMPdensSave(meshpoints-I+1)**2.d0 * (1.d0 - 
     &    nu_losses(SX(SX_RHO,SX_CNTR-SX_SURF-I+3)))
      ENDDO
      WRITE(io_unit,'(999(E18.12,1x))') WIMPdensSave   ! WIMP density [cm^-3]
      WRITE(io_unit,'(999(E18.12,1x))') mfpSave        ! WIMP mean free path [cm] 
      WRITE(io_unit,'(999(E18.12,1x))') L_WIMPtransport(meshpoints:2:-1)
        ! Luminosity carried conductively by WIMPs [erg/s] (negative of effective luminosity produced)
      WRITE(io_unit,'(999(E18.12,1x))') E_transportSave 
        ! Energy injection rate due to WIMP conductive energy transport [erg/g/s]
      WRITE(io_unit,'(999(E18.12,1x))') E_annihilationSave
        ! Energy injection rate due to WIMP annihilations [erg/g/s]
      CLOSE(io_unit)
      END SUBROUTINE Write_Status_Info

      


      double precision function safe_log10(x)
         double precision, intent(IN) :: x
         if (x < 1d-99) then
            safe_log10 = -99d0
         else
            safe_log10 = log10(x)
         end if
      end function safe_log10


      SUBROUTINE Save_Experiment_Info
      USE EZ_driver
      INTEGER :: ios, io_num, i
      CHARACTER (LEN=strlen) :: out_filename

      DOUBLE PRECISION, PARAMETER :: Lsun = 1d33 * CLSN
      DOUBLE PRECISION, PARAMETER :: Msun = 1d33 * CMSN
      DOUBLE PRECISION, PARAMETER :: CSDAY = 86.4d3

      ! these are the parameters of the shell
      integer, parameter :: iM=1 !  mass interior to this meshpoint (Msolar)
      integer, parameter :: ilogR=2 !  log10 radius (Rsolar)
      integer, parameter :: ilogT=3 !  log10 temperature
      integer, parameter :: ilogRHO=4 !  log10 density
      integer, parameter :: iL=5 !  luminosity (Lsolar)
      integer, parameter :: iXH=6 !   H1 mass fraction
      integer, parameter :: iXHE=7 !   He4 mass fraction
      integer, parameter :: iXC=8 !   C12 mass fraction
      integer, parameter :: iXN=9 !   N14 mass fraction
      integer, parameter :: iXO=10 !   O16 mass fraction
      integer, parameter :: iXNE=11 !   Ne20 mass fraction
      integer, parameter :: iXMG=12 !   Mg24 mass fraction
        
      integer, parameter :: num_args=12
    
      ! these are the results (cgs units)
      integer, parameter :: r_logP=1 !  pressure
      integer, parameter :: r_logKAP=2 !  opacity
      integer, parameter :: r_logU=3 !  specific energy
      integer, parameter :: r_logS=4 !  specific entropy
      integer, parameter :: r_mu=5 !  mean molecular weight per free particle
      integer, parameter :: r_free_e=6 !  mean number of free electrons per nucleon
      integer, parameter :: r_logCP=7 !  specific heat capacity at constant pressure (ergs per gram per K) 
      integer, parameter :: r_gamma1=8 ! one of the adiabatic exponents
      integer, parameter :: r_grad_ad=9 !  adiabatic temperature gradient
      integer, parameter :: r_logCS=10 !  sound speed
      integer, parameter :: r_eta=11 !  electron degeneracy
      integer, parameter :: r_grad_rad=12 !  radiative temperature gradient to carry luminosity
      integer, parameter :: r_grad_star=13 !  actual temperature gradient
      integer, parameter :: r_logDiffCoef=14 !  diffusion coefficient for convective mixing
      integer, parameter :: r_logCV=15 !  convection velocity -- zero if not in a convection zone
      integer, parameter :: r_logeps_phot=16 !  photon energy generation rate from nuclear reactions
      integer, parameter :: r_logneut_nuc=17 !  neutrino energy loss rate from nuclear reactions
      integer, parameter :: r_logEPS_PP=18 !  PP chains energy generation rate
      integer, parameter :: r_logEPS_CNO=19 !  CNO cycles energy generation rate
      integer, parameter :: r_logEPS_3alpha=20 !  triple alpha
      integer, parameter :: r_logEPS_alpha_cap=21 !  alpha captures on C, N, O, and Ne
      integer, parameter :: r_logEPS_metals=22  !   C+C, C+O, and O+O
      integer, parameter :: r_logNEU_plas=23 !  from plasmon neutrinos 
      integer, parameter :: r_logNEU_brem=24 !  from bremsstrahlung
      integer, parameter :: r_logNEU_pair=25 !  from pair annihilation
      integer, parameter :: r_logNEU_phot=26 !  from photo neutrinos
      integer, parameter :: r_t_therm=27 !  thermal timescale (yrs), := (star_Mass - M) * Cp * T / star_L
      integer, parameter :: r_t_eddy=28 !  convective eddy timescale (yrs) := mixing length / convection velocity
      integer, parameter :: r_dQ_dk=29 !  mesh spacing
      integer, parameter :: r_PQ=30 !  mesh pressure term
      integer, parameter :: r_TQ=31 !  mesh temperature term
      integer, parameter :: r_RQ=32 !  mesh radius term
      integer, parameter :: r_MQ=33 !  mesh mass term
      integer, parameter :: r_dM_dk=34 !  from mesh function
        
      integer, parameter :: num_results=34
      
      INTEGER, PARAMETER :: out_num = num_args + num_results
      DOUBLE PRECISION :: args(num_args), results(num_results), star_L

      ios = 0; io_num = 18
      
      if (model_Number < 10) then
         write(out_filename,'(a,i1,a)') 'exp_model_000', model_Number, '.data'
      else if (model_Number < 100) then
         write(out_filename,'(a,i2,a)') 'exp_model_00', model_Number, '.data'
      else if (model_Number < 1000) then
         write(out_filename,'(a,i3,a)') 'exp_model_0', model_Number, '.data'
      else if (model_Number < 10000) then
         write(out_filename,'(a,i4,a)') 'exp_model_', model_Number, '.data'
      else if (model_Number < 100000) then
         write(out_filename,'(a,i5,a)') 'exp_model_', model_Number, '.data'
      else
         write(out_filename,'(a,i6,a)') 'exp_model_', model_Number, '.data'
      end if

      OPEN(UNIT=io_num, FILE=trim(out_filename), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(*,*) 'Save_Experiment_Info: Failed to open file ', TRIM(out_filename), '     ios =', ios, '    io_num =', io_num
         RETURN
      END IF
      WRITE(*,*) '  Save experiment log info to ', TRIM(out_filename)
 101  FORMAT(99(E25.15,2x))

      WRITE(io_num,101) initial_Mass, initial_Z, star_Age, mixing_length_alpha, 
     >      overshoot_param_H, overshoot_param_He, mesh_core_mass
      
      CALL EZ_Extras
      
      star_L = SX(SX_L,N_SURF_shell)

      DO i = 1, N_SHELLs

         args(iM) = SX(SX_M,i) !  mass interior to this meshpoint (Msolar)
         args(ilogR) = safe_log10(SX(SX_R,i)) !  log10 radius (Rsolar)
         args(ilogT) = safe_log10(SX(SX_T,i)) !  log10 temperature
         args(ilogRHO) = safe_log10(SX(SX_RHO,i)) !  log10 density
         args(iL) = SX(SX_L,i) !  luminosity (Lsolar)
         args(iXH) = SX(SX_XH,i) !   H1 mass fraction
         args(iXHE) = SX(SX_XHE,i) !   He4 mass fraction
         args(iXC) = SX(SX_XC,i) !   C12 mass fraction
         args(iXN) = SX(SX_XN,i) !   N14 mass fraction
         args(iXO) = SX(SX_XO,i) !   O16 mass fraction
         args(iXNE) = SX(SX_XNE,i) !   Ne20 mass fraction
         args(iXMG) = SX(SX_XMG,i) !   Mg24 mass fraction

         results(r_logP) = safe_log10(SX(SX_P,i)) !  pressure
         results(r_logKAP) = safe_log10(SX(SX_OPACITY,i)) !  opacity
         results(r_logU) = safe_log10(SX(SX_U,i)) !  specific energy
         results(r_logS) = safe_log10(SX(SX_S,i)) !  specific entropy
         results(r_mu) = SX(SX_MU,i) !  mean molecular weight per free particle
         results(r_free_e) = SX(SX_NE,i) !  moles of free electrons per gram
         results(r_logCP) = safe_log10(SX(SX_SCP,i)) !  specific heat capacity at constant pressure (ergs per gram per K) 
         results(r_gamma1) = SX(SX_GAMMA1,i) ! one of the adiabatic exponents
         results(r_grad_ad) = SX(SX_GRAD_AD,i) !  adiabatic temperature gradient
         results(r_logCS) = LOG10(DSQRT(SX(SX_GAMMA1,i) * SX(SX_P,i) / SX(SX_RHO,i)))
         results(r_eta) = SX(SX_PSI,i)
         results(r_grad_rad) = SX(SX_GRAD_RAD,i) !  radiative temperature gradient required to carry entire luminosity without convection
         results(r_grad_star) = SX(SX_GRAD_STAR,i) !  actual temperature gradient
         results(r_logDiffCoef) = safe_log10(SX(SX_SG,i)) !  diffusion coefficient for convective mixing (if is 0, then no mixing)
         results(r_logCV) = safe_log10(SX(SX_CV,i)) !  convection velocity -- zero if not in a convection zone
         results(r_logeps_phot) = safe_log10(SX(SX_EPS_NUC,i)) !  photon energy generation rate from nuclear reactions
         results(r_logneut_nuc) = safe_log10(SX(SX_EPS_NUC_NEU,i)) !  neutrino energy loss rate from nuclear reactions
         results(r_logEPS_PP) = safe_log10(SX(SX_EPS_PP,i)) !  PP chain energy generation rate
         results(r_logEPS_CNO) = safe_log10(SX(SX_EPS_CNO,i)) !  CNO cycle energy generation rate
         results(r_logEPS_3alpha) = safe_log10(SX(SX_EPS_3A,i)) !  triple alpha
         results(r_logEPS_alpha_cap) = safe_log10(SX(SX_EPS_AC,i)+SX(SX_EPS_AN,i)+SX(SX_EPS_AO,i)+SX(SX_EPS_ANE,i))
         results(r_logEPS_metals) = safe_log10(SX(SX_EPS_CCA,i)+SX(SX_EPS_CO,i)+SX(SX_EPS_OO,i))
         results(r_logNEU_plas) = safe_log10(SX(SX_NEU_plasma,i)) !  from plasmon neutrinos 
         results(r_logNEU_brem) = safe_log10(SX(SX_NEU_brem,i)) !  from bremsstrahlung
         results(r_logNEU_pair) = safe_log10(SX(SX_NEU_pair,i)) !  from pair annihilation
         results(r_logNEU_phot) = safe_log10(SX(SX_NEU_photo,i)) !  from photo neutrinos
         results(r_t_therm) = (star_Mass - SX(SX_M,i)) * Msun * SX(SX_SCP,i) * SX(SX_T,i) / (CSDAY * Lsun * star_L)
         if (SX(SX_CV,i) == 0d0) then !  convective eddy timescale := mixing length / convection velocity
            results(r_t_eddy) = 0d0
         else
            results(r_t_eddy) = SX(SX_WL,i) / (CSDAY * SX(SX_CV,i)**2)
         end if
         results(r_dQ_dk) = SX(SX_DQ_DK,i) !  mesh spacing
         results(r_PQ) = SX(SX_PQ,i)
         results(r_TQ) = SX(SX_TQ,i)
         results(r_RQ) = SX(SX_RQ,i)
         results(r_MQ) = SX(SX_MQ,i)
         results(r_dM_dk) = SX(SX_DM_DK,i)

         WRITE(io_num,101) args, results

      END DO
      CLOSE(io_num)
      END SUBROUTINE Save_Experiment_Info
      
      END MODULE ez_do_one_utils
      
