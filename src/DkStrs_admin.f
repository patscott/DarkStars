! Administrative module for DarkStars
!
! Pat Scott, March 2008; pat@fysik.su.se
!--------------------------------------------------------------------------------

      MODULE DkStrs_admin
      use star_controls
      use ez_driver
      use ez_report
      use DkStrs_data
      use DkStrs_utils
      use DkStrs_init

      implicit none

      contains


      DOUBLE PRECISION FUNCTION timestep_max_orbital(orbit_t_equiv)

      double precision, intent(IN) :: orbit_t_equiv
      double precision :: rhoonv, timestep_max_orbital2, t_advanced, rhoonv_advanced
      double precision, parameter :: search_scale = 0.99d0, hold_on_scale = 3.d0

      rhoonv = rhowimp/v_star

      if (loopOrbit) then

        timestep_max_orbital = orbit_len*orbit_monotonicity_scale

        do
          t_advanced = orbit_t_equiv+timestep_max_orbital
          rhoonv_advanced = ambientDens(t_advanced)/starVelocity(t_advanced)
          if (abs((rhoonv - rhoonv_advanced)/rhoonv).lt.hold_on*hold_on_scale) exit
          timestep_max_orbital = timestep_max_orbital*search_scale
        enddo

        if (sign(1.d0, drhoonvdt(orbit_t_equiv)).ne.sign(1.d0, drhoonvdt(t_advanced))) then

          timestep_max_orbital2 = timestep_max_orbital

          do
            rhoonv_advanced = ambientDens(orbit_t_equiv+timestep_max_orbital2)/starVelocity(orbit_t_equiv+timestep_max_orbital2)
            if (abs((rhoonv - rhoonv_advanced)/rhoonv).gt.hold_on*hold_on_scale) exit
            if (timestep_max_orbital2*CSY.lt.dynamic_Timescale) then
              timestep_max_orbital2 = timestep_max_orbital + 1.d0
              exit
            endif
            timestep_max_orbital2 = timestep_max_orbital2*search_scale
          enddo

          if (timestep_max_orbital2.lt.timestep_max_orbital) then
            do
              rhoonv_advanced = ambientDens(orbit_t_equiv+timestep_max_orbital2)/starVelocity(orbit_t_equiv+timestep_max_orbital2)
              if (abs((rhoonv - rhoonv_advanced)/rhoonv).lt.hold_on*hold_on_scale) exit
              timestep_max_orbital2 = timestep_max_orbital2*search_scale
            enddo
            timestep_max_orbital = timestep_max_orbital2
          endif

        endif

      else

        timestep_max_orbital = min_timestep

        do
          t_advanced = orbit_t_equiv+timestep_max_orbital
          rhoonv_advanced = ambientDens(t_advanced)/starVelocity(t_advanced)
          if (abs((rhoonv - rhoonv_advanced)/rhoonv).gt.hold_on*hold_on_scale) exit
          timestep_max_orbital = timestep_max_orbital/search_scale
          if (timestep_max_orbital.gt.orbit_len) then
            timestep_max_orbital = -1.d0
            exit
          endif
        enddo

        if(timestep_max_orbital.gt.0.d0) timestep_max_orbital = timestep_max_orbital*search_scale

      endif

      if (timestep_max_orbital.ge.0.d0 .and. timestep_max_orbital.lt.min_timestep) then
        if(verbose) then
          write(*,*) 'COMMENT: For stability, now hardsetting the requested timestep to min_timestep,'
          write(*,*) 'even though the timestep preferred here on the given orbit is less than this.'
          write(*,*) 'The influence of some short-timescale features of the orbit may be missed.'
        endif
        timestep_max_orbital = min_timestep
      endif

      return

      END FUNCTION timestep_max_orbital


      LOGICAL FUNCTION Time_For_Periodic_Save()
      !Determines whether to dump a periodic save file or not, according to the save period
      !specified in the the input file.

      logical :: needs_Save
      integer(8) :: count_rate, right_now
      character (len=8) :: currentDate
      character (len=10) :: currentTime

      if (save_freq-prec.lt.prec) then
         Time_For_Periodic_Save = .false.
         return
      endif

      if (start_time.lt.0.d0) call SYSTEM_CLOCK(start_time,count_rate)
      call SYSTEM_CLOCK(right_now,count_rate)

      if (right_now.lt.start_time) then
         write(*,*) 'System clock has ticked over - doing periodic save now and resetting counter.'
         start_time = right_now
         needs_Save = .true.
      else
         if (dble(right_now-start_time)/count_rate/3600./24. .gt. dble(save_count + 1)*save_freq) then
            save_count = save_count + 1
            needs_Save = .true.
         else
            needs_Save = .false.
         endif
      endif

      if (needs_Save) then
         call DATE_AND_TIME(currentDate, currentTime)
         write(*,*) 'Current time is ', currentTime(1:2), ':', currentTime(3:4), ':',  currentTime(5:6), ', ',
     &    currentDate(7:8), '-', currentDate(5:6), '-',  currentDate(1:4)
      endif

      Time_For_Periodic_Save = needs_Save

      END FUNCTION Time_For_Periodic_Save


      LOGICAL FUNCTION DarkStars_Save(filename)

      character (len=strlen), intent(IN) :: filename
      integer, parameter :: lun = 100
      integer :: fileStatus

      open (UNIT=lun, FILE=filename, ACTION='WRITE', STATUS='REPLACE', IOSTAT=fileStatus, FORM='UNFORMATTED')
      if (fileStatus.ne.0) call die_quietly('Opening file for saving DarkStars data failed.')

      write (lun) DoWimpyThings
      write (lun) n_WIMPs, n_WIMPs_prev, n_WIMPs_saved, n_WIMPs_supersaved
      write (lun) rhowimp, v_star, rhoonv_prev, rhoonv_saved
      write (lun) cap, ann, evap, effcap, dNdt, r_chi, Tw
      write (lun) JMOD_prev, meshpoints
      write (lun) n_species
      write (lun) new_model, LTE_invalid
      if (JMOD_prev.gt.2) write(lun) prelim_Age
      if (new_model) then
        write (lun) JMOD_new_model_start, new_model_age
        if (JMOD_prev .gt. JMOD_new_model_start+1) write (lun) prelim_Age2
      endif
      write (lun) stellar_age, centreH, centrerho, centreT
      write (lun) K_suppression, timestep_max
      write (lun) partition_func, n_WIMPs_centre_LTE
      write (lun) starr, starm, starmd2, starphi, starphid2, starrho
      write (lun) starmfr, starmfrd2
      write (lun) staralpha, staralphad2
      write (lun) startemperature, startemperatured1, temptension
      !write (lun) starcondEff, starcondEffd1, condEfftension
      if (.not.LTE_invalid) write (lun) starLTErhoWIMP, starLTErhoWIMPd1, LTErhoWIMPtension
      if (.not.LTE_invalid) write (lun) L_WIMPtransport, L_WIMPtransportd1, L_WIMPtransporttension
      if (fileStatus.ne.0) call die_quietly('Writing DarkStars data to save file failed.')

      close(lun,IOSTAT=fileStatus)
      if (fileStatus.ne.0) call die_quietly('Closing file for saving DarkStars data failed.')

      DarkStars_Save = .true.

      END FUNCTION DarkStars_Save


      LOGICAL FUNCTION DarkStars_Restore(filename)

      character (len=strlen), intent(IN) :: filename
      integer, parameter :: lun = 100
      integer :: fileStatus, meshpoints_in, n_species_in

      open (UNIT=lun, FILE=filename, ACTION='READ', IOSTAT=fileStatus, FORM='UNFORMATTED')
      if (fileStatus.ne.0) call die_quietly('Opening file for restoring DarkStars data failed.')

      read (lun) WimpyThingsDonePreviously
      read (lun) n_WIMPs, n_WIMPs_prev, n_WIMPs_saved, n_WIMPs_supersaved
      read (lun) rhowimp, v_star, rhoonv_prev, rhoonv_saved
      read (lun) cap, ann, evap, effcap, dNdt, r_chi, Tw
      read (lun) JMOD_prev, meshpoints_in
      if (meshpoints_in.ne.meshpoints) call die_quietly('Incompatible meshpoint number in restore file.')
      read (lun) n_species_in
      if (n_species_in.ne.n_species) call die_quietly('Incompatible n_species number in restore file.')
      read (lun) new_model, LTE_invalid
      if (JMOD_prev.gt.2) read (lun) prelim_Age
      if (new_model) then
        read (lun) JMOD_new_model_start, new_model_age
        if (JMOD_prev .gt. JMOD_new_model_start+1) read (lun) prelim_Age2
      endif
      read (lun) stellar_age, centreH, centrerho, centreT
      read (lun) K_suppression, timestep_max
      read (lun) partition_func, n_WIMPs_centre_LTE
      read (lun) starr, starm, starmd2, starphi, starphid2, starrho
      read (lun) starmfr, starmfrd2
      read (lun) staralpha, staralphad2
      read (lun) startemperature, startemperatured1, temptension
      !read (lun) starcondEff, starcondEffd1, condEfftension
      if (.not.LTE_invalid) read (lun) starLTErhoWIMP, starLTErhoWIMPd1, LTErhoWIMPtension
      if (.not.LTE_invalid) read (lun) L_WIMPtransport, L_WIMPtransportd1, L_WIMPtransporttension
      if (fileStatus.ne.0) call die_quietly('Reading DarkStars data from restore file failed.')

      close(lun,IOSTAT=fileStatus)
      if (fileStatus.ne.0) call die_quietly('Closing file for restoring DarkStars data failed.')

      call Complete_Model
      call EZ_Extras
      call init_radmass
      after_restore = .true.
      DarkStars_Restore = .true.

      END FUNCTION DarkStars_Restore


      END MODULE DkStrs_admin
