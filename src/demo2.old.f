      MODULE demo2
      USE ez_do_one
      USE star_extras
      IMPLICIT NONE
     
      logical :: collapsing, adjusting_extra_Energy
      double precision :: init_center_H
      integer :: save_collapse_model

      CONTAINS


      SUBROUTINE Demo_2
      USE ez_do_one_data
      USE ez_flash
      DOUBLE PRECISION :: mass
      CALL Demo_2_Info
      IF (.TRUE.) THEN ! USE THIS FOR DOING A SINGLE MASS FOR A SINGLE Z
         ! the default demo2 case is mass = 3.0 Msun and Z = 0.02
         mass = 3.0d0  ! the initial mass of the star
         
         metals_dir = '../metals/z02'  ! where to get metals data
         data_dir = '../demos/data_z02' ! where to write the output logs
         
         CALL Init_Do_One_Utils ( mass )
         CALL Open_Demo2_DIFFME()
         
         WRITE_PROFILE_TO_RUN = .TRUE.
         WRITE_PROFILE_TO_DATA = .TRUE.
         WRITE_PROFILE_TO_MODELS = .FALSE.

!        CALL Get_Star_Info( mass, Dummy_Build_Filename, Experiment_Initial_Params, Demo2_Before_Evolve, Do_One_Check_Model )

!        CALL Get_Star_Info( mass, Dummy_Build_Filename, Experiment_Initial_Params, Demo2_Before_Evolve, Bill_Check_Model )

!        CALL Get_Star_Info( mass, Dummy_Build_Filename, Demo2_Initial_Params, Demo2_Before_Evolve, Do_One_Check_Model )

!        CALL Get_Star_Info( mass, Dummy_Build_Filename, Demo2_Initial_Params, Demo2_PreMS_Before_Evolve, Do_Pre_MS_Check_Model )

         CALL Get_Star_Info( mass, Dummy_Build_Filename, stella_Initial_Params, stella_Before_Evolve, stella_Check_Model )

!        CALL Get_Star_Info( mass, Dummy_Build_Filename, Bill_Initial_Params, Bill_Before_Evolve, Bill_Check_Model )

!        CALL Get_Star_Info( mass, Dummy_Build_Filename, Phil_Initial_Params, Phil_Before_Evolve, Phil_Check_Model )

         !WRITE_PROFILE_TO_MODELS = .TRUE.
         !CALL Get_Star_Info( mass, Dummy_Build_Filename, Lars_Initial_Params, Lars_Before_Evolve, Lars_Check_Model )

         RETURN
      END IF
      IF (.FALSE.) CALL Demo_2_All ! use this for doing all masses for Z=0.02
      IF (.FALSE.) CALL Demo_2_Z ! use this for doing selected masses for all the Zs
      WRITE (*,'(/,A,/)') '  End of Demo2'
      END SUBROUTINE Demo_2

   
     
      INTEGER FUNCTION stella_Check_Model()

      character (len=256) :: filename
      integer :: i, io_unit
      integer, parameter :: svars = 13
      double precision :: shell(svars)
      io_unit = 40
      
      IF (model_Number < 882) then
         if (mod(model_Number,50)==0) write(*,*) model_Number
         stella_Check_Model = KEEP_GOING
         return
      end if
      
      call EZ_Extras
      filename = 'stella_star.data'
      open(unit=io_unit,file=trim(filename))
      do i = N_SHELLs, 1, -1
         shell(1:6) = (/ SX(SX_R,i)*CRSN*1d11, SX(SX_M,i)*CMSN*1d33, 0d0, SX(SX_L,i)*CLSN*1d33, SX(SX_RHO,i), SX(SX_T,i) /)
         shell(7:12) = ( / SX(SX_XH,i), SX(SX_XHE,i), SX(SX_XC,i), SX(SX_XN,i), SX(SX_XO,i), SX(SX_XNE,i) /) 
         shell(svars) = 1d0 - sum(shell(7:12))
         write(io_unit,'(99e28.14)') shell
      end do
      close(io_unit)
      stella_Check_Model = TERMINATE
      
      END FUNCTION stella_Check_Model


      SUBROUTINE stella_Before_Evolve
      call Demo2_Before_Evolve
      END SUBROUTINE stella_Before_Evolve

      
      SUBROUTINE stella_Initial_Params
      
      mixing_length_alpha = 1.6D0
      overshoot_param_H = 0.12D0
      overshoot_param_He = 0.12D0
      wind_Eta = 0d0
      number_of_center_shells_to_mix = 0
      
      END SUBROUTINE stella_Initial_Params

   
     
      INTEGER FUNCTION Bill_Check_Model()
      ! for project with Phil and Eric.
      integer :: k, i, flg
      double precision :: t_therm
      DOUBLE PRECISION, PARAMETER :: Lsun = 1d33 * CLSN
      DOUBLE PRECISION, PARAMETER :: Msun = 1d33 * CMSN
      DOUBLE PRECISION, PARAMETER :: CSDAY = 86.4d3
      
      flg = Do_One_Check_Model()
      IF (flg /= KEEP_GOING) then
         Bill_Check_Model = flg
         return
      end if
      Bill_Check_Model = KEEP_GOING

      if (.false.) then
      
      if (center_He > 0.57d0) timestep_max = 5d7
      if (center_He > 0.615d0) timestep_max = 1d6
      if (center_He > 0.626d0 .and. center_He < 0.627d0) then
         CALL Save_Experiment_Info
      else if (center_He >= 0.627d0) then ! this is solar
         CALL Save_Experiment_Info
         Bill_Check_Model = TERMINATE
      END IF

      end if
      
      call EZ_Extras
      if (.false.) then ! use this to force mixing at center
      do i = SX_CNTR, SX_SURF, -1
         if (SX(SX_M,i) > 0.1d0) then
            number_of_center_shells_to_mix = SX_CNTR - i
            exit
         end if
      end do
      else
      number_of_center_shells_to_mix = 0
      end if
      
      if (star_Age > 4d8) CALL Save_Profiles
      
      if (center_He > 0.9d0) Bill_Check_Model = TERMINATE
      END FUNCTION Bill_Check_Model


      SUBROUTINE Bill_Before_Evolve
      call Demo2_Before_Evolve
      head_cnt = 10; summary_cnt = 2
      END SUBROUTINE Bill_Before_Evolve

      
      SUBROUTINE Bill_Initial_Params
      
      mixing_length_alpha = 1.6D0
      overshoot_param_H = 0.12D0
      overshoot_param_He = 0.12D0
      wind_Eta = 0d0
      number_of_center_shells_to_mix = 0
      call_save_exp_info = .true.
      
      END SUBROUTINE Bill_Initial_Params



      INTEGER FUNCTION Phil_Check_Model()
      ! for project with Phil and Eric.
      integer :: k, i, flg
      double precision :: t_therm
      DOUBLE PRECISION, PARAMETER :: Lsun = 1d33 * CLSN
      DOUBLE PRECISION, PARAMETER :: Msun = 1d33 * CMSN
      DOUBLE PRECISION, PARAMETER :: CSDAY = 86.4d3
      
      flg = Do_One_Check_Model()
      IF (flg /= KEEP_GOING) then
         Phil_Check_Model = flg
         return
      end if
      Phil_Check_Model = KEEP_GOING

      if (.true.) then
      if (center_He > 0.57d0) timestep_max = 5d7
      if (center_He > 0.615d0) timestep_max = 1d6
      if (center_He > 0.626d0 .and. center_He < 0.627d0) then
         CALL Save_Experiment_Info
      else if (center_He >= 0.627d0) then ! this is solar
         CALL Save_Experiment_Info
         !Phil_Check_Model = TERMINATE
      else if (star_Age >= 0.8d9) then
         Phil_Check_Model = TERMINATE
      END IF
      end if
      
      call EZ_Extras
      if (.false.) then ! use this to force mixing at center
      do i = SX_CNTR, SX_SURF, -1
         if (SX(SX_M,i) > 0.1d0) then
            number_of_center_shells_to_mix = SX_CNTR - i
            exit
         end if
      end do
      else
      number_of_center_shells_to_mix = 0
      end if
      
      if (mod(model_Number,20) == 0 .or. log10(star_Age) > 8.5) CALL Save_Profiles
      
      if (center_He > 0.9d0) Phil_Check_Model = TERMINATE
      END FUNCTION Phil_Check_Model


      SUBROUTINE Phil_Before_Evolve
      call Demo2_Before_Evolve
      head_cnt = 10; summary_cnt = 2
      END SUBROUTINE Phil_Before_Evolve
      
      
      SUBROUTINE Experiment_Initial_Params
      timestep_decrement = 0.3D0
      timestep_hold = 3
      PSI_limit = 1d6
      GAM_limit = 155d0
      head_cnt = 100; summary_cnt = 20
      mixing_length_alpha = 1.6D0
      overshoot_param_H = 0.12D0
      overshoot_param_He = 0.12D0
      wind_Eta = 0d0
      number_of_center_shells_to_mix = 0
      do_post_He_flash = .false.
      END SUBROUTINE Experiment_Initial_Params

      
      SUBROUTINE Phil_Initial_Params
      
      mixing_length_alpha = 1.6D0
      overshoot_param_H = 0.12D0
      overshoot_param_He = 0.12D0
      wind_Eta = 0d0
      number_of_center_shells_to_mix = 0
      call_save_exp_info = .true.
      
      END SUBROUTINE Phil_Initial_Params


      
      SUBROUTINE Lars_Save_Experiment_Info
      INTEGER :: ios, io_num, i
      CHARACTER (LEN=strlen) :: out_filename
      INTEGER, PARAMETER :: out_num = 24
      DOUBLE PRECISION :: out_data(out_num)
      ios = 0; io_num = 18
      out_filename = 'exp_out.data'
      OPEN(UNIT=io_num, FILE=trim(out_filename), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(*,*) 'Init_Log: Failed to open file ', TRIM(out_filename), '     ios =', ios, '    io_num =', io_num
         RETURN
      END IF
      WRITE(*,*) '  Save experiment log info to ', TRIM(out_filename)
      CALL EZ_Extras
      DO i = 1, N_SHELLs
          out_data(1:4) = (/ SX(SX_M,i), SX(SX_R,i), SX(SX_T,i), SX(SX_RHO,i) /)
          out_data(5:11) = (/ SX(SX_XH,i), SX(SX_XHE,i), SX(SX_XC,i), SX(SX_XN,i), SX(SX_XO,i), SX(SX_XNE,i), SX(SX_XMG,i) /)
          out_data(12:14) = (/ mixing_length_alpha, overshoot_param_H, overshoot_param_He /)
          out_data(15:19) = (/ SX(SX_P,i), SX(SX_SCP,i), SX(SX_OPACITY,i), SX(SX_GRAD_AD,i), SX(SX_GRAD_RAD,i) /)
          out_data(20:out_num) = (/ SX(SX_GRAD_STAR,i), SX(SX_CV,i), SX(SX_SG,i), SX(SX_GAMMA1,i), SX(SX_L,i) /)
         WRITE(io_num,'(99(E25.15,x))') out_data
      END DO
      CLOSE(io_num)
      END SUBROUTINE Lars_Save_Experiment_Info

     
      INTEGER FUNCTION Lars_Check_Model()
      
      Lars_Check_Model = Do_One_Check_Model()
      
      !if (helium_ignition) summary_cnt = 2
      summary_cnt = 5
      
      
      END FUNCTION Lars_Check_Model

      
      SUBROUTINE Lars_Initial_Params
      
      overshoot_param_H = 0d0
      overshoot_param_He = 0d0
      wind_Eta = 0d0
      do_post_He_flash = .false.
      
      END SUBROUTINE Lars_Initial_Params


      SUBROUTINE Lars_Before_Evolve
      CHARACTER (LEN=strlen) :: filename

      timestep_decrement = 0.3D0
      timestep_hold = 3
      PSI_limit = 1d6
      GAM_limit = 155d0
      head_cnt = 100; summary_cnt = 20
      !stop_at_Model = 504
      return

      ! use the following to save and restore a selected model
      filename = 'ez.sav'
      IF (.not. EZ_Restore(filename)) stop 'failed to restore'
      head_cnt = 10; summary_cnt = 20
      stop_at_Model = -1
     END SUBROUTINE Lars_Before_Evolve

     SUBROUTINE Demo2_Before_Evolve
      CHARACTER (LEN=strlen) :: filename
      call_save_exp_info = .true.
      timestep_decrement = 0.3D0
      timestep_hold = 3
      PSI_limit = 1d6
      GAM_limit = 155d0
      head_cnt = 100; summary_cnt = 20
      !stop_at_Model = 977  use this to cause a particular model to be saved for later restore
      return
      ! use the following to save and restore a selected model
      filename = 'ez.sav'
      IF (.not. EZ_Restore(filename)) stop 'failed to restore'
      head_cnt = 10; summary_cnt = 1
      stop_at_Model = -1
     END SUBROUTINE Demo2_Before_Evolve

      ! for 0.6M premain sequence
      ! extra_Energy_max = 25d0
      ! extra_Energy_time = 1d6
      ! extra_Energy_param = 5d-2
      ! stop_at_Model = 499

      ! for 1.0M premain sequence
      ! extra_Energy_max = 32d0
      ! extra_Energy_time = 1d6
      ! extra_Energy_param = 1d-3
      ! stop_at_Model = 500

      ! for 2.0M premain sequence
      ! extra_Energy_max = 36d0
      ! extra_Energy_time = 1d6
      ! extra_Energy_param = 5d-2
      ! stop_at_Model = 499

      ! for 2.2M premain sequence
      ! extra_Energy_max = 37d0
      ! extra_Energy_time = 1d6
      ! extra_Energy_param = 1d-3
      ! stop_at_Model = 500

      ! for 3.0M premain sequence
      ! extra_Energy_max = 120d0
      ! extra_Energy_time = 1d5
      ! extra_Energy_param = 5d-3
      ! stop_at_Model = 475

      ! for 4.0M premain sequence
      ! extra_Energy_max = 218d0
      ! extra_Energy_time = 1d5
      ! extra_Energy_param = 5d-3
      ! stop_at_Model = 599

      ! for 5.0M premain sequence
      ! extra_Energy_max = 400d0
      ! extra_Energy_time = 1d5
      ! extra_Energy_param = 5d-3
      ! stop_at_Model = 599

      ! for 9.0M premain sequence
      ! extra_Energy_max = 1650d0
      ! extra_Energy_time = 1d5
      ! extra_Energy_param = 5d-3
      ! stop_at_Model = 599

      ! for 15.0M premain sequence
      ! extra_Energy_max = 4800d0
      ! extra_Energy_time = 1d5
      ! extra_Energy_param = 5d-3
      ! stop_at_Model = 599

     
     SUBROUTINE Demo2_PreMS_Before_Evolve
      use ez_cycle
      CHARACTER (LEN=strlen) :: filename
      logical :: flg
      summary_cnt = 1
      
      
      collapsing = .false.
      adjusting_extra_Energy = .false.
      timestep_decrement = 0.3D0
      timestep_hold = 3
      PSI_limit = 500d0
      GAM_limit = 150d0
      head_cnt = 100; summary_cnt = 20
      extra_Energy_max = 120d0
      extra_Energy_time = 1d5
      extra_Energy_param = 5d-3
      stop_at_Model = 475
      save_collapse_model = 508
      wind_Eta = 0
      burn_H = 0
      extra_Mdot_param = 0
      !return

      ! use the following to save and restore a selected model
      filename = 'ez.sav'
      IF (.not. EZ_Restore(filename)) stop 'failed to restore'
      head_cnt = 30; summary_cnt = 2
      extra_Energy_max = 1d-8
      extra_Energy_time = 5d2
      !extra_Energy_param = 1d-5
      burn_H = 1
      wind_Eta = 1D0
      timestep_max = 1d5
      call Reset_Age_and_Model_Number(1d5, 1)
      stop_at_Model = 2000
      save_collapse_model = 4
      profile_AGE = -1d0
      init_center_H = center_H
      collapsing = .true.
      CALL EZ_Extras
      del_log_Luminosity = 0.3d0
      del_log_surface_Temp = 0.125d0
      del_log_center_Temp = 0.3d0
      del_log_center_Density = 0.75d0
      prv_log_Luminosity = log_Luminosity
      prv_log_surface_Temp = log_surface_Temp
      prv_log_center_Temp = log_center_Temp
      prv_log_center_Density = log_center_Density
     END SUBROUTINE Demo2_PreMS_Before_Evolve


   
     
      INTEGER FUNCTION Do_Pre_MS_Check_Model()
         use ez_setup
         
         double precision, parameter :: target_log_Rho = -1.6d0
         double precision, parameter :: d_log_Rho = 0.005d0
         double precision, parameter :: XH_frac = 0.999d0
         integer, parameter :: min_model_gap = 0
         
         Do_Pre_MS_Check_Model = Do_One_Check_Model()
         
         if (.not. collapsing) then
         
            if (.not. adjusting_extra_Energy) then
               if (log_center_Density < -1.5d0) adjusting_extra_Energy = .true.
            ! Note that for the following cases, both log_center_Density and target_log_Rho are negative
            else if (log_center_Density > target_log_Rho + d_log_Rho) then ! increase extra_Energy_max
               extra_Energy_max  = extra_Energy_max * (target_log_Rho + d_log_Rho) / log_center_Density
            else if (log_center_Density < target_log_Rho - d_log_Rho) then ! decrease extra_Energy_max
               extra_Energy_max  = extra_Energy_max * (target_log_Rho - d_log_Rho) / log_center_Density
            end if
            
            if (abs(star_Mass - initial_Mass) > 1d-4) then ! adjust mass
               extra_Mdot_param = (initial_Mass - star_Mass) / timestep_max
            else ! leave it alone
               extra_Mdot_param = 0d0
            end if
            
            return
            
         end if
         
         if (Do_Pre_MS_Check_Model /= KEEP_GOING .or. recent_profile_model == model_Number) return

         if (model_Number < 3) then
            prv_log_Luminosity = log_Luminosity
            prv_log_surface_Temp = log_surface_Temp
            prv_log_center_Temp = log_center_Temp
            prv_log_center_Density = log_center_Density
         end if
         
         if (model_Number < save_collapse_model) return

         if (center_H < init_center_H * XH_frac) then
            Do_Pre_MS_Check_Model = TERMINATE
            return
         end if

         Do_Pre_MS_Check_Model = Do_Delta_Check_Model()

      END FUNCTION Do_Pre_MS_Check_Model
      
      
     SUBROUTINE Demo_2_All
     INTEGER, PARAMETER :: num_Masses = 51
     DOUBLE PRECISION :: mass, mass_Array(num_Masses)
     INTEGER :: first_to_do, last_to_do, i
     mass_Array(1:9)   = (/   0D0,   0D0, 0.3D0, 0.4D0, 0.5D0, 0.6D0, 0.7D0, 0.8D0, 0.9D0 /)
     mass_Array(10:19) = (/ 1.0D0, 1.1D0, 1.2D0, 1.3D0, 1.4D0, 1.5D0, 1.6D0, 1.7D0, 1.8D0, 1.9D0 /)
     mass_Array(20:29) = (/ 2.0D0, 2.1D0, 2.2D0, 2.3D0, 2.4D0, 2.5D0, 2.6D0, 2.7D0, 2.8D0, 2.9D0 /)
     mass_Array(30:39) = (/ 3.0D0, 3.1D0, 3.2D0, 3.3D0, 3.4D0, 3.5D0, 3.6D0, 3.7D0, 3.8D0, 3.9D0 /)
     mass_Array(40:45) = (/   4D0,   5D0,   6D0,   7D0,   8D0,   9D0 /)
     mass_Array(46:51) = (/  10D0,  12D0,  16D0,  25D0,  40D0, 100D0 /)
     metals_dir = '../metals/z02'
     data_dir = '../demos/data_z02'
     first_to_do = 1
     last_to_do = num_Masses
     DO i = first_to_do, last_to_do
        mass = mass_Array(i)
        IF ( mass .LT. min_Mass ) CYCLE
        CALL Init_Do_One_Utils( mass )
        WRITE_PROFILE_TO_DATA = .TRUE.
        CALL Get_Star_Info ( mass, Build_Demo2_Filename, Demo2_Initial_Params, Dummy_Before_Evolve, Do_One_Check_Model )
     END DO
     END SUBROUTINE Demo_2_All
     
     SUBROUTINE Demo_2_Z
     INTEGER, PARAMETER :: num_MassesZ = 7
     DOUBLE PRECISION :: mass_ArrayZ(num_MassesZ)
     mass_ArrayZ(1:num_MassesZ)   = (/ 1D0, 2D0, 3D0, 4D0, 6D0, 12D0, 25D0 /)
     CALL Demo_2_Zs(num_MassesZ, mass_ArrayZ, 1, 6, .TRUE.)
     END SUBROUTINE Demo_2_Z
     
     SUBROUTINE Do_ZRGB_Exp
     INTEGER, PARAMETER :: num_MassesZ = 5
     DOUBLE PRECISION :: mass_ArrayZ(num_MassesZ)
     
     mass_ArrayZ(1:num_MassesZ)   = (/ 0.8d0, 1.0d0, 1.2d0, 1.4d0, 1.6d0 /)
     CALL Demo_2_Zs(num_MassesZ, mass_ArrayZ, 6, 7, .FALSE.)
     END SUBROUTINE Do_ZRGB_Exp
     
     SUBROUTINE Demo_2_Zs(num_MassesZ, mass_ArrayZ, first_Z, last_Z, postHeFlash)
     INTEGER :: first_Z, last_Z, num_MassesZ
     LOGICAL :: postHeFlash
     DOUBLE PRECISION :: mass_ArrayZ(num_MassesZ), mass
     INTEGER :: first_to_doZ, last_to_doZ, i, j 
     first_to_doZ = 1; last_to_doZ = num_MassesZ
     DO j = first_Z, last_Z
        SELECT CASE (j)
        CASE (1)
           metals_dir = '../metals/z0001'
           data_dir = '../demos/data_z0001'
        CASE (2)
           metals_dir = '../metals/z0003'
           data_dir = '../demos/data_z0003'
        CASE (3)
           metals_dir = '../metals/z001'
           data_dir = '../demos/data_z001'
        CASE (4)
           metals_dir = '../metals/z004'
           data_dir = '../demos/data_z004'
        CASE (5)
           metals_dir = '../metals/z01'
           data_dir = '../demos/data_z01'
        CASE (6)
           metals_dir = '../metals/z03'
           data_dir = '../demos/data_z03'
        CASE (7)
           metals_dir = '../metals/z02'
           data_dir = '../demos/data_z02'
        CASE DEFAULT
           STOP 'error with Z choice'
        END SELECT
        DO i = first_to_doZ, last_to_doZ
          mass = mass_ArrayZ(i)
          IF ( mass .LT. min_Mass ) CYCLE
          CALL Init_Do_One_Utils( mass )
          do_post_He_flash = postHeFlash
          CALL Open_Extra( data_dir, mass, Build_Demo2_Filename)
          WRITE_PROFILE_TO_DATA = .TRUE.
          CALL Get_Star_Info ( mass, Build_Demo2_Filename, Demo2_Initial_Params, Dummy_Before_Evolve, Do_One_Check_Model )
        END DO
     END DO
     END SUBROUTINE Demo_2_Zs
     
     SUBROUTINE Demo_2_Info
     WRITE (*,*)
     WRITE (*,*)
     WRITE (*,'(a)') ' The Demo_2 tests are with convective overshooting and a Reimers'' wind (eta=1.0).'
     WRITE(*,'(a)') ' The terminal output has the following information:'
     WRITE(*,'(a)') '      MOD# is model number,'
     WRITE(*,'(a)') '      RHO is center density,'
     WRITE(*,'(a)') '      Tc is center temperature,'
     WRITE(*,'(a)') '      L is luminosity (Lsolar),'
     WRITE(*,'(a)') '      Ts is surface temperature,'
     WRITE(*,'(a)') '      Age is years from start of run,'
     WRITE(*,'(a)') '      DT is timestep in years between models,'
     WRITE(*,'(a)') '      L_He is total power from helium burning,'
     WRITE(*,'(a)') '      R is radius (Rsolar),'
     WRITE(*,'(a)') '      Ttl_Mass is total stellar mass (Msolar),'
     WRITE(*,'(a)') '      He_Core is total mass inside location where hydrogen mass abundance 1st reaches 15%,'
     WRITE(*,'(a)') '      Noncore is Ttl_Mass minus He_Core,'
     WRITE(*,'(a)') '      XH1, XHe4, XC12, XN14, and XO16 are center mass fractions,'
     WRITE(*,'(a)') '      PSIcntr is the central electron degeneracy.'
     WRITE(*,'(a)') '      GAMcntr is the central plasma interaction parameter.'
     END SUBROUTINE Demo_2_Info
     
     SUBROUTINE Demo2_Initial_Params
     wind_Eta = 1D0  ! set Reimers' wind eta
     END SUBROUTINE Demo2_Initial_Params
     
     SUBROUTINE Null_Before_Evolve
     END SUBROUTINE Null_Before_Evolve
     
     SUBROUTINE Open_Demo2_DIFFME
     INTEGER :: ios
     CHARACTER (LEN=strlen) :: diffme_filename
     diffme_filename = 'demo2.DIFFME'
     WRITE_DIFFME = .TRUE.
     IO_DIFFME = 28
     ios = 0
     OPEN(UNIT=IO_DIFFME, FILE=trim(diffme_filename), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
     IF (ios /= 0) THEN
        WRITE(*,*) 'Build_Demo2_Filename: Failed to open file ', TRIM(diffme_filename)
     END IF
     END SUBROUTINE Open_Demo2_DIFFME
     
     SUBROUTINE Build_Demo2_Filename(init_M)
     DOUBLE PRECISION, INTENT(IN) :: init_M
     INTEGER :: im, tenths
     CHARACTER (LEN = 100) :: frmt
     IF (LEN(TRIM(data_dir)) .EQ. 0) THEN
         full_name = trim(partial_name)
         RETURN
     END IF
     im = floor(init_M)
     tenths = floor(10D0*(init_M-im+.05))
     IF (im .GE. 100) THEN
        frmt = '(2A,I3,A,I1,2A)'
     ELSEIF (im .GE. 10) THEN
        frmt = '(2A,I2,A,I1,2A)'
     ELSE
        frmt = '(2A,I1,A,I1,2A)'
     END IF
     WRITE(full_name,frmt) trim(data_dir), '/demo2_', im, '_', tenths, '/', trim(partial_name)
     WRITE(*,*) full_name
     END SUBROUTINE Build_Demo2_Filename
       
     SUBROUTINE Get_logF_EoS
     USE ez_setup
     DOUBLE PRECISION :: logTmin, logTmax, logFmin, logFmax
     INTEGER :: logTpoints, logFpoints
     CHARACTER (LEN=strlen) :: metals_dir, out_filename
     metals_dir = '../metals/z02'
     out_filename = '../profile/logF_EoS.data'
     logTmin = 3D0; logTmax = 8.5D0; logTpoints = 100
     logFmin = -5D0; logFmax = 2.1D0; logFpoints = 100
     CALL EZ_Start( metals_dir )
     IF (.NOT. EZ_ZAMS(1D0, Null_Initial_Params)) RETURN
     ionization_level = 2
     IF (.FALSE.) THEN
        out_filename = 'test1_logF_EoS.data'
        WRITE(*,*) 'Writing ', trim(out_filename)
        CALL Test_logF_EoS(4.45D0, 0.464732135555556D0, out_filename)
        out_filename = 'test2_logF_EoS.data'
        WRITE(*,*) 'Writing ', trim(out_filename)
        CALL Test_logF_EoS(4.45D0, 0.464732135666667D0, out_filename)
     ELSE
        WRITE(*,*) 'Writing ', trim(out_filename)
        CALL Write_logF_EoS(logTmin, logTmax, logTpoints, logFmin, logFmax, logFpoints, out_filename)
        WRITE(*,*) 'Finished writing ', trim(out_filename)
     END IF
     
     END SUBROUTINE Get_logF_EoS
     
     SUBROUTINE Get_logRHO_EoS
     USE ez_setup
     DOUBLE PRECISION :: lgTmn, lgTmx, lgRHOmn, lgRHOmx
     INTEGER :: lgTsteps, lgRHOsteps
     CHARACTER (LEN=strlen) :: metals_dir, out_dir
     metals_dir = '../metals/z02'
     out_dir = '.'
     ! full range
     lgTmn = 6.2D0; lgTmx = 8.2D0; lgTsteps = 201
     lgRHOmn = -3.00D0; lgRHOmx = 6D0; lgRHOsteps = 201
     ! for pp_cno
     !lgTmn = 6.75D0; lgTmx = 7.5D0; lgTsteps = 201
     !lgRHOmn = -1D0; lgRHOmx = 2.50D0; lgRHOsteps = 201
     WRITE(*,*) 'Writing logRHO_EoS'
     CALL EZ_Start( metals_dir )
     IF (.NOT. EZ_ZAMS(1D0, Null_Initial_Params)) RETURN
     CALL Write_logRHO_EoS(lgTmn, lgTmx, lgTsteps, lgRHOmn, lgRHOmx, lgRHOsteps, out_dir)
     WRITE(*,*) 'Finished writing logRHO_EoS'
     END SUBROUTINE Get_logRHO_EoS
     
      SUBROUTINE Get_EPSNUC
      USE ez_setup
      DOUBLE PRECISION :: lgTmn, lgTmx, lgRHOmn, lgRHOmx, xin(NEL)
      INTEGER :: lgTsteps, lgRHOsteps
      CHARACTER (LEN=strlen) :: metals_dir, out_dir
      metals_dir = '../metals/z02'
      out_dir = 'EPSNUC_DATA'
      xin = 0
      if (.false.) then ! for pp_cno
         xin(N_H)  =   0.655186E+00    ! h1   
         xin(N_HE) =   0.312021E+00 + 0.647994D-06    ! he4 + he3 
         xin(N_C)  =   0.002725D-01   ! c12
         xin(N_N)  =   0.203101D-01 + 0.612124D-06  + 0.109305D-02  + 0.356004D-04   ! n14+n15 + o17  +c13
         xin(N_O)  =   0.094000D-01    ! o16 
         xin(N_NE) =   0.162163D-02    ! ne20 
         lgTmn = 6.2D0; lgTmx = 8.7D0; lgTsteps = 201
         lgRHOmn = -3.00D0; lgRHOmx = 7D0; lgRHOsteps = 201
      else if (.true.) then! for helium burning, early
         xin(N_H)  =   0E+00    ! h1   
         xin(N_HE) =   0.90d0    ! he4 + he3 
         xin(N_C)  =   0.01d0   ! c12
         xin(N_N)  =   0.05d0   ! n14+n15 + o17  +c13
         xin(N_O)  =   0.01d0    ! o16 
         xin(N_NE) =   0.01d0    ! ne20 
         lgTmn = 7.4D0; lgTmx = 8.7D0; lgTsteps = 201
         lgRHOmn = 0.00D0; lgRHOmx = 9D0; lgRHOsteps = 201
      else if (.false.) then! for helium burning, mid
         xin(N_H)  =   0E+00    ! h1   
         xin(N_HE) =   0.18d0    ! he4 + he3 
         xin(N_C)  =   0.40d0   ! c12
         xin(N_N)  =   0.8d-05   ! n14+n15 + o17  +c13
         xin(N_O)  =   0.40d0    ! o16 
         xin(N_NE) =   0.01D0    ! ne20 
         lgTmn = 7.4D0; lgTmx = 8.7D0; lgTsteps = 201
         lgRHOmn = 0.00D0; lgRHOmx = 9D0; lgRHOsteps = 201
      else if (.false.) then! for helium burning, late
         xin(N_H)  =   0E+00    ! h1   
         xin(N_HE) =   0.001d0    ! he4 + he3 
         xin(N_C)  =   0.49d0   ! c12
         xin(N_N)  =   0.8d-05   ! n14+n15 + o17  +c13
         xin(N_O)  =   0.49d0    ! o16 
         xin(N_NE) =   0.01D0    ! ne20 
         lgTmn = 7.4D0; lgTmx = 8.7D0; lgTsteps = 201
         lgRHOmn = 0.00D0; lgRHOmx = 9D0; lgRHOsteps = 201
      else if (.false.) then! for C/O burning, late
         xin(N_H)  =   0E+00    ! h1   
         xin(N_HE) =   0d0    ! he4 + he3 
         xin(N_C)  =   0.49d0   ! c12
         xin(N_N)  =   0.8d-05   ! n14+n15 + o17  +c13
         xin(N_O)  =   0.49d0    ! o16 
         xin(N_NE) =   0.01D0    ! ne20 
         lgTmn = 7.75D0; lgTmx = 9D0; lgTsteps = 201
         lgRHOmn = 0.00D0; lgRHOmx = 10.75D0; lgRHOsteps = 201
      end if
      WRITE(*,*) 'Writing EPSNUC_DATA'
      CALL EZ_Start( metals_dir )
      IF (.NOT. EZ_ZAMS(1D0, Null_Initial_Params)) RETURN
      CALL Write_EPSNUC(lgTmn, lgTmx, lgTsteps, lgRHOmn, lgRHOmx, lgRHOsteps, out_dir, xin)
      WRITE(*,*) 'Finished writing EPSNUC_DATA'
      END SUBROUTINE Get_EPSNUC
     
     
     END MODULE demo2
     
     
     
     
     
     
