! DarkStars
!--------------------------------------------------------------------------------
! Solves equations of WIMP capture, annihilation, stellar structure and evolution
! simultaneously for ZAMS stars of (nearly) arbitrary mass and metallicity.
!
! Includes modified parts of the EZ and DarkSUSY codes, as well as inbuilt
! copies of TSPACK, LAPACK and the zeroin function from the GO package. 
! 
! Written by Pat Scott, July 2007 & Feb/March 2008; pat@fysik.su.se
!--------------------------------------------------------------------------------

      MODULE DarkStars
      use ez_driver
      use ez_do_one
      use ez_do_one_data
      use ez_data
      use star_controls
      use DkStrs_data
      use DkStrs_admin      
      use DkStrs_init
      use DkStrs_WIMPdens
      implicit none

      logical :: DoMovie, ReconvergeModels
      character(len=strlen) :: restore_filename
      double precision :: initial_timestep, rhowimp_initial, v_star_initial
      double precision :: max_Age_in, galesc_initial, galr_initial, epsH
              
      contains


      SUBROUTINE Do_DarkStar
            
      double precision :: mass

      namelist /directories/ metals_dir, main_dir, data_dir, prof_dir, movie_dir
      namelist /switches/ capmode, interpmode, potmode, chopatvesc, altveldist,
     & WIMPdensmode, ReconvergeModels, DoTransport, DoMovie, DoWimpyThings, verbose,
     & integrator
      namelist /controls/ hold_on, initial_timestep, timestep_rescale, min_timestep, 
     & summary_cnt, stop_at_Model, max_Age, save_filename, restore_filename, save_freq
      namelist /physics/ mx, sigsi, sigsd, sigann, nuLossFactor, boost_factor, vd_3d_star, v_star,
     & rhowimp, galesc, galr, n_WIMPs, mass, K0, Knudsen_suppression_tau, tau_therm      
      namelist /abundances/ solar_C, solar_O, heavy_abuns, iso_ratios
        
!...Splash screen
      write(*,*) 
      write(*,*) '-----------------------------------------------------------------'
      write(*,*) 'Welcome to DarkStars '//trim(Darkstars_version)//' (incorportating DarkSUSY and EZ/STARS)'
      write(*,*) '-----------------------------------------------------------------'
      write(*,*) 

!...Initialize DarkStars
      write(*,*) 'Initializing DarkStars...'
      write(*,*) 
      
!...Read input parameters from input pipe
      write(*,*) '  Reading simulation input parameters...'
      read(unit=5,nml=directories)
      read(unit=5,nml=switches)
      read(unit=5,nml=controls)
      read(unit=5,nml=physics)
      read(unit=5,nml=abundances)
      metals_dir_in = metals_dir
      save_filename_in = save_filename
      stop_at_Model_in = stop_at_Model
      summary_cnt_in = summary_cnt
      if ((stop_at_Model_in.ge.0 .or. save_freq.ge.0)  .and. save_filename.eq.'') call 
     & die_quietly('  Error: incompatible save_filename and stop_at_Model/save_freq in input file.')
      if (rhowimp*v_star .le. 0.d0) call 
     & die_quietly('  Error: incompatible rhowimp and v_star in input file.')
      if (galesc*galr .le. 0.d0) call 
     & die_quietly('  Error: incompatible galesc and galr in input file.')
      if (.not.capmode .and. altveldist) call 
     & die_quietly('  Error: capmode = false and altveldist = true in input file.')
      if (nuLossFactor.gt.1.d0) call 
     & die_quietly('  Error: nuLossParameter must be <= 1.')
      if (tau_therm .gt. prec .and. WIMPdensmode) call
     & die_quietly('  Error: positive tau_therm not allowed with WIMPdensmode = true.') 
      data_dir = trim(main_dir)//'/'//trim(data_dir)
      movie_dir = trim(main_dir)//'/'//trim(movie_dir)
      prof_dir = trim(main_dir)//'/'//trim(prof_dir)
      n_WIMPs_saved = n_WIMPs
      n_WIMPs_prev = n_WIMPs
      JMOD_prev = 0
      JMOD_new_model_start = 0
      NotATest = .true.
      rhowimp_initial = rhowimp
      v_star_initial = v_star
      galesc_initial = galesc
      galr_initial = galr
      max_Age_in = max_Age
      data_dir_in = data_dir
      epsH = mass*1d-14
      write(*,*) '  done.'
      write(*,*) 

!...Initialize DarkSUSYLE
      write(*,*) '  Initializing DarkSUSYLE...'
      call dsLEinit
      write(*,*) '  done.'
      write(*,*) 
      
!...Set up other static data
      write(*,*) '  Inializing WIMP thermal diffusivities and conductivities...'
      call alpha_kappa_precalc
      write(*,*) '  done.'
      write(*,*) 
      write(*,*) '  Inializing WIMP-nucleus cross-sections...'
      call sigma_precalc
      write(*,*) '  done.'
      write(*,*) 
        
!...Set up orbits
      if (rhowimp .le. 0.d0) then
        write(*,*) '  Inializing stellar orbit...'
        call orbit_init
        write(*,*) '  done.'
        write(*,*)
      endif
      rhoonv_prev = rhowimp/v_star
        
!...Determine machine precision
      write(*,*) '  Determining machine precision...'
      prec = mach_tol_dbl()
      write(*,*) '  done.'
      write(*,*)
      
!...Finished initialization
      write(*,*) 'done.'
      write(*,*) 

      call DarkStar_init (mass)
        
      call DkStrs_Get_Star_Info(mass, Dummy_Build_Filename, DarkStar_Extra_Filename, 
     & DarkStar_Initial_Params, DarkStar_Before_Evolve, DarkStar_Check_Model)

      write(*,*)
      write(*,'(/,A,/)') 'DarkStars run completed.'

      return 
      
      END SUBROUTINE Do_DarkStar
      

      SUBROUTINE DarkStar_init(init_M)
      
      double precision, intent(IN) :: init_M

      profile_model = -1
      recent_profile_model = -1
      save_helium_ignition = .true.
      WRITE_DIFFME = .false.
      WRITE_PROFILE_TO_RUN = .true.
      WRITE_PROFILE_TO_MODELS = DoMovie
      WRITE_BRIEF_TO_RUN = .false.
      WRITE_PROFILE_TO_DATA = .true.
      call_save_exp_info = .false.
      tau_Model = 0
      target_tau_Photosphere = -1D0 ! not changing tau_Photosphere, so set to negative value to indicate that.
      do_post_He_flash = ( init_M .ge. lower_He_flash .and. init_M .le. upper_He_flash )
      do_CNTR_RHOs = init_M .le. no_He_ignition_limit
      do_CNTR_T_drops = init_M .le. no_CNTR_T_drops_limit
      phase_of_evolution = phase_STARTING
      post_He_AGE = -1D0
      profile_AGE = -1D0
      prev_CNTR_RHO = 1D99
      helium_ignition = .false.
      number_of_little_steps = 0
      del_log_Luminosity = 1d6
      del_log_surface_Temp = 1d6
      del_log_center_Temp = 1d6
      del_log_center_Density = 1d6
      cap = 0.d0
      evap = 0.d0
      effcap = 0.d0
      ann = 0.d0
      new_model = .false.
      new_model_step1 = .false.
      after_restore = .false.
        
      END SUBROUTINE DarkStar_init


      SUBROUTINE DarkStar_Extra_Filename(init_M, what)

      double precision, intent(IN) :: init_M
      integer, intent(IN) :: what

      select case (what)
        case (1)
          full_name = trim(prof_dir)//'/'//trim(partial_name)
        case (2)
          full_name = trim(movie_dir)//'/'//trim(partial_name)
        case default
          call die_quietly('Not clear what filename to build.')
      end select
 
      END SUBROUTINE DarkStar_Extra_Filename


      SUBROUTINE DarkStar_Initial_Params

      wind_Eta = 1.d0  ! set Reimers' wind eta
      timestep_lower_limit = 0.8D0
      timestep_upper_limit = 1.5D0
      timestep_decrement = 0.9D0
      timestep_hold = 3
        
      if (DoWimpyThings) then
        call init_annihilation(.true.)
        if (rhowimp.ne.0.d0) timestep_max = initial_timestep
      endif
 
      END SUBROUTINE DarkStar_Initial_Params

      
      SUBROUTINE DarkStar_Before_Evolve
        
        character (len=strlen) :: EZ_filename, DarkStars_filename
        
      if (DoWimpyThings) then

        if (WIMPdensmode) then
          write(*,*) 'WIMP radial distribution mode set to accurate.'
        else
          write(*,*) 'WIMP radial distribution mode set to approximate.'
        endif
      
        if (interpmode) then
          write(*,*) 'Interpolation mode set to advanced.'
          write(*,*) '  (Cubic or automatically-tensioned splines, depending upon application)'
        else
          write(*,*) 'Interpolation mode set to linear.'    
          write(*,*) 'This mode has yet to be fully implemented.'
          call die_quietly('DarkStars execution cancelled.')
        endif

        if (potmode) then 
          write(*,*) 'Potential mode set to explicit.'
        else
          write(*,*) 'Potential mode set to tabulate.'
        endif      

        if (capmode) then
          write(*,*) 'Halo velocity integration mode set to numerical.'
        else
          write(*,*) 'Halo velocity integration mode set to analytical.'
        endif
 
        if (chopatvesc) then
          write(*,*) 'Halo velocity distribution will be truncated at galactic escape velocity.'
        else
          write(*,*) 'Halo velocity distribution will extend to infinity.'
        endif
      
        if (altveldist) then
          write(*,*) 'Using a user-defined distribution for halo velocities.'
        else
          write(*,*) 'Using an isothermal Maxwell-Boltzmann distribution for halo velocities.'
        endif
              
      else

        write(*,*) 'All WIMP-related effects disabled; proceeding with normal evolution...'

      endif
      
      write(*,*)
      write(*,*) 'Logs to be saved to '//trim(data_dir)
      write(*,*) 'Profiles to be saved to '//trim(prof_dir)
      if (DoMovie) then
        write(*,*) 'Movie data to be saved to '//trim(movie_dir)
      else
        write(*,*) 'No movie though :('
      endif
      write(*,*)

      if(rhoWIMP_initial.lt.0.d0) then
        if(loopOrbit) then
          write(*,*) 'Stellar orbit will loop.'
        else
            write(*,*) 'Stellar orbit will not loop.'
        endif
        write(*,*)
      endif

      if(ReconvergeModels) then
        write(*,*) 'Stellar models will be reconverged at each timestep.'
      else
        write(*,*) 'Stellar models will not be reconverged.'
      endif
      write(*,*)

      head_cnt = summary_cnt * 5

      if (restore_filename.ne.'') then
        EZ_filename = trim(main_dir)//'/'//trim(restore_filename)
        DarkStars_filename = trim(EZ_filename)//trim(DarkStars_extension)      
        write(*,*) 'Restoring saved star from '//trim(EZ_filename)//' and'
        write(*,*) trim(DarkStars_filename)//'...'
        if (.not.EZ_Restore(EZ_filename)) call die_quietly('Failed to restore '//EZ_filename)
        if (.not.DarkStars_Restore(DarkStars_filename)) call die_quietly('Failed to restore '//DarkStars_filename)
        if (.not.WIMPyThingsDonePreviously) then
          new_model_step1 = .true.
          call init_annihilation(.true.)
        endif
        write(*,*) 'done.'
        write(*,*)
        if (ReconvergeModels) call Model_IO_internal(.false.)
      endif        
        
      metals_dir = metals_dir_in
      data_dir = data_dir_in
      save_filename = save_filename_in
      stop_at_Model = stop_at_Model_in
      summary_cnt = summary_cnt_in
      max_Age = max_Age_in
      if (rhowimp_initial .gt. 0.d0) then      
        rhowimp = rhowimp_initial
        v_star = v_star_initial
        if (galesc_initial .gt. 0.d0) then
          galesc = galesc_initial
          galr = galr_initial
        endif
      endif

      END SUBROUTINE DarkStar_Before_Evolve


      INTEGER FUNCTION DarkStar_Check_Model

      use DkStrs_transport

      character (len=strlen) :: save_count_str, model_str, EZ_filename, DarkStars_filename
      double precision :: equiv_t, equiv_t_advanced, equiv_t_advanced_insides, radius, rho, HPVAL
      double precision :: equil_time, ann_parameter, dscapstar, dsLEfint, H_depletion_rate
      double precision :: orbit_t_equiv, timestep_max_from_orbit
      double precision :: Hshift, rhoshift, Tshift, shift_temp_array(4), shift_temp_loc(1)
      double precision, parameter :: epsrho = 1.d-10, epsT = 1.d-10, universe_Age = 1.4d10
      integer :: i,j, exit_code

      LOGICAL :: must_do_log, stop_because_He_ignited, GoOn = .true.
      DOUBLE PRECISION, PARAMETER :: log_He_Temp = 7.8D0, core_limit = 5d-3
      DOUBLE PRECISION, PARAMETER :: d_tau_min = 1D-2, d_tau_max = 1D0
      DOUBLE PRECISION, PARAMETER :: little_step_Factor = 10D0, little_step_Size = 10D0
      DOUBLE PRECISION :: d_tau
      INTEGER :: model
      INTEGER, PARAMETER :: tau_ramp = 50, num_little_steps_Limit = 100
        
      select case (JMOD-JMOD_prev)
        case (0)
          if (DoWimpyThings .and. phase_of_evolution .eq. phase_STARTING
     &     .and. (stop_at_Model .eq. -4 .or. stop_at_Model .eq. -10)) then
            do i = 1, 4 
              centreH(i) = centreH(i+1)
              centrerho(i) = centrerho(i+1)
              centreT(i) = centreT(i+1)
              stellar_age(i) = stellar_age(i+1)
            enddo
          endif
        case (1)
          if (NotATest .and. ReconvergeModels .and. JMOD.gt.2 .and. JMOD-JMOD_new_model_start.gt.2) NotATest = .false.
          rhoonv_saved = rhoonv_prev
          rhoonv_prev = rhowimp/v_star  
          JMOD_prev = JMOD
        case default
          if (JMOD-JMOD_prev.eq.2 .and. helium_ignition .and. .not.new_model) then
            write(*,*)
            write(*,*) 'Helium flash and generation of a new He-burning model acknowledged'
            write(*,*) 'by DarkStars; existing WIMP population will be adapted to new star.'
            write(*,*)
            JMOD_prev = JMOD
            JMOD_new_model_start = JMOD
            n_WIMPs_saved = n_WIMPs
            n_WIMPs_prev = n_WIMPs
            rhoonv_prev = rhowimp/v_star
            cap = 0.d0
            evap = 0.d0
            effcap = 0.d0
            ann = 0.d0
            new_model_age = star_Age
            new_model = .true.
            new_model_step1 = .true.
            if (DoWimpyThings) call init_annihilation(.true.)
          else if (JMOD - JMOD_prev .eq. -1 .and. ReconvergeModels) then
            write(*,*) 'WARNING: extreme reconvergence backup - current model has'
            write(*,*) 'been reconverged with the WIMP properties of the next model.'
            rhoonv_prev = rhoonv_saved
            NotATest = .true.
            JMOD_prev = JMOD
            rhoonv_prev = rhowimp/v_star
            if (DoWimpyThings .and. phase_of_evolution .eq. phase_STARTING
     &       .and. (stop_at_Model .eq. -4 .or. stop_at_Model .eq. -10)) then
              do i = 1, 3 
                centreH(i) = centreH(i+2)
                centrerho(i) = centrerho(i+2)
                centreT(i) = centreT(i+2)
                stellar_age(i) = stellar_age(i+2)
              enddo
              centreH(4) = centreH(3)
              centrerho(4) = centrerho(3)
              centreT(4) = centreT(3)
              stellar_age(4) = stellar_age(3)
            endif
          else
            write(*,*) 'JMOD - JMOD_prev = ', JMOD - JMOD_prev
            call die_quietly('Strange JMOD numbers in DarkStar_Check_Model.')
          endif
      end select
      extra_Mdot_param = -1.d-14
          
      if (.not.NotATest) then
        if (DoWimpyThings) then
          call init_capture
          call init_annihilation(.false.)
          if (DoTransport) call init_transport
        endif
        call Model_IO_internal(.true.)
        n_WIMPs_prev = n_WIMPs_saved
        n_WIMPs_saved = n_WIMPs_supersaved
        NotATest = .true.
        DarkStar_Check_Model = KEEP_GOING
        return
      endif
      
      ! Exit if EZ has screwed with the initial age
      if (star_Age .lt. prec .and. JMOD .ge. 1 .and. rhowimp_initial .lt. 0.d0) 
     & call die_quietly('Results are rubbish because EZ sent star_Age negative. Reduce initial timestep.')

      ! Skip all WIMPy things if normal evolution has been selected
      if (DoWimpyThings) then
        
        ! Don't allow there to be less than zero WIMPs
        ! (in principle one should never see this message...)
        if (n_WIMPs.lt.0.) then
          write(*,*) 'WARNING: Extreme WIMP loss - setting n_WIMPs to zero...'
          n_WIMPs = 0.d0
        endif

        ! Do things associated with the stellar orbit
        if (JMOD.eq.3) prelim_Age=star_Age
        if (new_model) then
          if (JMOD.eq.JMOD_new_model_start+2) prelim_Age2 = star_Age - new_model_age 
        endif
        if (rhowimp_initial.lt.0.d0 .and. JMOD.ge.3) then
          if (new_model) then
            if (JMOD-JMOD_new_model_start.lt.2) then
              orbit_t_equiv = new_model_age-prelim_Age
            else
              orbit_t_equiv = star_Age-prelim_Age-prelim_Age2
            endif
          else
            orbit_t_equiv = star_Age-prelim_Age
          endif  
          if (loopOrbit) orbit_t_equiv = orbit_t_equiv + orbital_t(fileLen/2)
          rhowimp = ambientDens(orbit_t_equiv)
          v_star = starVelocity(orbit_t_equiv)
          if (galesc_initial.lt.0.d0) then
            galesc = galactic_vesc(orbit_t_equiv)
            galr = galactic_r(orbit_t_equiv)
          endif
          timestep_max_from_orbit = timestep_max_orbital(orbit_t_equiv)
        endif
      
        ! Compute capture, annihilation (and evaporation) rates
        call init_capture
        cap = CSY * boost_factor * dscapstar(mx,sigsi,sigsd)
        call init_annihilation(.false.)
        if (DoTransport) call init_transport
        if (n_WIMPs.ne.0.d0) then
          idtag = 'annihilation'
          ann = CSY * dsLEfint(WIMP_annihilation_integrand,0.d0,1.d2*r_star,1.d-6)
          if (ann.le.0.d0) ann = CSY * dsLEfint(WIMP_annihilation_integrand,0.d0,1.d2*r_star/32.d0,1.d-6)
          if (ann.le.0.d0) ann = CSY * dsLEfint(WIMP_annihilation_integrand,0.d0,1.d2*r_star/1024.d0,1.d-6)
          if (ann.le.0.d0) ann = CSY * real(fint_ext(W_ann_backup_integrand,0._16,1.e2_16*real(r_star,16)/1024._16,1.e-6_16),8)
          if (ann.le.0.d0) call die_quietly('Error: ann = 0 in DarkStar_Check_Model.')
        else
          ann = 0.d0
        endif
        evap = 0.d0
        effcap = cap - evap
        dNdt = effcap - 2.d0*ann
              
        ! Limit next timestep so that the WIMP population doesn't change by more than
        ! hold_on times the current population in a single step (or at least try to;
        ! this can be overidden by EZ).
        if (dNdt.ne.0.d0) then
          if (abs(ann).le.prec) then
            timestep_max = hold_on*n_WIMPs/abs(dNdt)
          else 
            ann_parameter = 2.d0*ann/n_WIMPs**2
            if (abs(effcap).le.prec) then
              timestep_max = hold_on/((hold_on+1.d0)*n_WIMPs*ann_parameter)
            else
              equil_time = 1.d0/sqrt(effcap*ann_parameter)
              equiv_t = tanhinv((n_WIMPs/effcap/equil_time)**sign(1.d0,dNdt))
              equiv_t_advanced_insides = ((1.d0+hold_on*sign(1.d0,dNdt))*n_WIMPs/effcap/equil_time)**sign(1.d0,dNdt)
              if (abs(equiv_t_advanced_insides).lt.1.d0) then
                equiv_t_advanced = tanhinv(equiv_t_advanced_insides)
                timestep_max = equil_time*(equiv_t_advanced-equiv_t)
              else
                timestep_max = -1.d0
              endif
            endif
          endif
        else
          timestep_max = -1.d0
        endif
      
        if (rhowimp_initial.lt.0.d0 .and. JMOD.ge.3 .and. ((timestep_max.eq.-1 .and. loopOrbit) .or. 
     &   (timestep_max_from_orbit.ne.-1 .and. timestep_max.gt.timestep_max_from_orbit))) timestep_max = timestep_max_from_orbit
        if (timestep_max.ne.-1.d0 .and. timestep_max.lt.initial_timestep) then
          if (JMOD.lt.2) timestep_max = initial_timestep
          if (new_model) then
           if (JMOD-JMOD_new_model_start.le.1) timestep_max = initial_timestep
          endif
        endif
      
        if (legacy) write(*,*) 'Total solar capture rate:',cap/CSY,'WIMPs per second'
        if (legacy) call die_quietly('DarkStars legacy run completed.')

        if (verbose) then
          write(*,*) 'Model number:', JMOD
          write(*,*) 'Total WIMPs in star:    ', n_WIMPs
          write(*,*) 'Capture rate:   ', cap, 'WIMPs/yr'
          write(*,*) 'Net accretion rate:  ', dNdt, 'WIMPs/yr'
          write(*,*) 'Annually injected annihilation energy: ', mx*2.*ann, 'GeV'
          write(*,*) 'WIMP thermal radius: ', r_chi/r_star*1.d-2, 'stellar radii'
          write(*,*) 'WIMP isothermal temperature: ', Tw, 'K'
          write(*,*) 'Central WIMP density: ', WIMPdens(0.d0), 'WIMPs/cm^3'
          write(*,*) 'WIMP density at thermal radius: ', WIMPdens(r_chi), 'WIMPs/cm^3'
          write(*,*) 'Current distance from centre of galaxy:   ', galr, 'pc'
          write(*,*) 'Local WIMP ambient density:   ', rhowimp, 'GeV/cm^3'
          write(*,*) 'Local galactic escape velocity:   ', galesc, 'km/s'
          write(*,*) 'Stellar velocity through halo:   ', v_star, 'km/s'
          write(*,*) 'Change in orbital rho/v (%) during last timestep: ', 1.d2*abs(rhowimp/v_star - rhoonv_prev)/rhoonv_prev 
          write(*,*) 'WIMP population change (%) during last timestep: ', 1.d2*abs(n_WIMPs - n_WIMPs_saved)/n_WIMPs_saved
          write(*,*) 'Last timestep: ', time_Step, 'yrs'
          write(*,*) 'Next requested maximum timestep: ', timestep_max, 'yrs'
          write(*,*)
        endif

      endif
      
      new_model_step1 = .false.
      after_restore = .false.
      must_do_log = .false.
      stop_because_He_ignited = .false.
      model = model_Number
      if ( model .eq. 0 ) next_CNTR_RHO = min_CNTR_RHO
  
      if ((max_AGE.gt.0.d0 .and. star_Age.gt.max_AGE) .or. model.eq.stop_at_Model) then 
        write(*,*) 'Star has reached designated maximum age - saving and quiting.'
        stop_at_Model = model
        DarkStar_Check_Model = TERMINATE
      else 
         DarkStar_Check_Model = KEEP_GOING
      endif

      if (DoWimpyThings .and. phase_of_evolution .eq. phase_STARTING
     & .and. (stop_at_Model .eq. -4 .or. stop_at_Model .eq. -10)) then
        do i = 5, 2, -1 
          centreH(i) = centreH(i-1)
          centrerho(i) = centrerho(i-1)
          centreT(i) = centreT(i-1)
          stellar_age(i) = stellar_age(i-1)
        enddo
        centreH(1) = center_H
        centrerho(1) = log_center_Density
        centreT(1) = log_center_Temp
        stellar_age(1) = star_Age
        if (JMOD .ge. 20 .and. star_Age .gt. 0.d0) then
          shift_temp_array = abs(centreH(1) - centreH(2:5))
          shift_temp_loc = maxloc(shift_temp_array)
          Hshift = maxval(shift_temp_array)
          Hshift = Hshift / (stellar_age(1) - stellar_age(1+shift_temp_loc(1)))
          shift_temp_array = abs(centrerho(1) - centrerho(2:5))
          shift_temp_loc = maxloc(shift_temp_array)
          rhoshift = maxval(shift_temp_array)
          rhoshift = rhoshift / (stellar_age(1) - stellar_age(1+shift_temp_loc(1)))
          shift_temp_array = abs(centreT(1) - centreT(2:5))
          shift_temp_loc = maxloc(shift_temp_array)
          Tshift = maxval(shift_temp_array)
          Tshift = Tshift / (stellar_age(1) - stellar_age(1+shift_temp_loc(1)))
          if (Hshift .lt. epsH .and. rhoshift .lt. epsrho .and. Tshift .lt. epsT) then
            write(*,*) 'Star has become a certified WIMP burner - saving and quitting.'
            stop_at_Model = model
            DarkStar_Check_Model = TERMINATE
          endif
        endif
      endif

      IF ( star_Age .LE. profile_AGE ) must_do_log = .TRUE. ! in case of backup, do it over
      IF ( save_helium_ignition .AND. (.NOT. helium_ignition) .AND. (log_center_Temp .GT. log_He_Temp) ) THEN
         CALL EZ_Extras
         IF (power_Metal_burn + power_He_burn > power_Neutrinos ) THEN
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
         IF (stop_at_Model .EQ. -8) THEN
            DarkStar_Check_Model = TERMINATE
            stop_at_Model = model
         ENDIF
      ELSE IF ( Time_To_Profile(GoOn) ) THEN
         must_do_log = .TRUE.
      ELSE IF ( model_Number .EQ. profile_model ) THEN
         must_do_log = .TRUE.
         WRITE(*,'(a,i5)') 'Save profile for model number:', model_Number
      END IF
      IF (.not.GoOn) THEN
        DarkStar_Check_Model = TERMINATE
        stop_at_Model = model
      ENDIF
      IF ( Log_State ( must_do_log ) ) THEN
         CALL Write_Logs
         CALL Write_Special_Log
         IF ( must_do_log ) CALL Save_Profiles
         IF (STOP_Sign()) THEN
            DarkStar_Check_Model = TERMINATE
         END IF
      END IF
      IF ( stop_because_He_ignited ) DarkStar_Check_Model = TERMINATE
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
            DarkStar_Check_Model = TERMINATE
         END IF
      END IF

      if (ReconvergeModels) call Model_IO_internal(.false.)
      n_WIMPs_supersaved = n_WIMPs_saved

      if (Time_For_Periodic_Save()) then
         write (save_count_str,'(I0)') save_count
         write (model_str,'(I0)') model
         EZ_filename = trim(main_dir)//'/'//trim(save_filename)//'_'//save_count_str
         DarkStars_filename = trim(EZ_filename)//trim(DarkStars_extension)
         if (EZ_Save(EZ_filename) .AND. DarkStars_Save(DarkStars_filename) ) then
            write(*,*) 'Periodic save '//trim(save_count_str)//' performed at model number '//trim(model_str)
         else
            write(*,*) 'FAILED trying to do periodic save at model number '//trim(model_str)
         END IF
      endif
      

      END FUNCTION DarkStar_Check_Model


      DOUBLE PRECISION FUNCTION WIMP_annihilation_integrand(r)
      ! Takes local height and returns integrand for finding total
      ! WIMP annihilation rate in a star
      ! Note that this integrand differs from Eq. 2.22 in Scott et al
      ! (2009, MNRAS 394:82) by a factor of 1/2, since the local annihilation
      ! rate should be WIMPdens^2 * sigann / 2.
      ! Input:    r   radius [cm]
      ! Output:       integrand [cm^-1]

      double precision, intent(IN) :: r      

      WIMP_annihilation_integrand = 2.d0*pi*r**2*sigann*WIMPdens(r)**2
      
      return   

      END FUNCTION WIMP_annihilation_integrand
        
        
      REAL(16) FUNCTION W_ann_backup_integrand(r)
      ! Extended precision version of WIMP_annihilation_integrand
      ! Input:    r   radius [cm]
      ! Output:       integrand [cm^-1]

      real(16), intent(IN) :: r
      real(16) :: tempWIMPdens
        
      if(LTE_invalid) then
        tempWIMPdens = WIMPdens_isothermal(real(r,8))
      else
        tempWIMPdens = real(K_suppression,16)*WIMPdens_LTE(real(r,8)) +
     &   (1._16-real(K_suppression,16))*WIMPdens_isothermal(real(r,8))
      endif
      W_ann_backup_integrand = 4._16 * real(pi,16) * r**2 * real(sigann,16) * tempWIMPdens**2
      
      return
        
      END FUNCTION W_ann_backup_integrand

      
      END MODULE DarkStars
