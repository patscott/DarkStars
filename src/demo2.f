      MODULE demo2
      USE ez_do_one
      USE star_extras
      IMPLICIT NONE
     
      integer :: he_break_even_model

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
            WRITE_PROFILE_TO_MODELS = .FALSE. ! change to .TRUE. for making movies

            CALL Get_Star_Info( mass, Dummy_Build_Filename, Demo2_Initial_Params, Demo2_Before_Evolve, Do_One_Check_Model )
            !CALL Get_Star_Info( mass, Dummy_Build_Filename, Lars_Initial_Params, Lars_Before_Evolve, Lars_Check_Model )

            RETURN
         END IF
         IF (.FALSE.) CALL Demo_2_All ! use this for doing all masses for Z=0.02
         IF (.FALSE.) CALL Demo_2_Z ! use this for doing selected masses for all the Zs
         WRITE (*,'(/,A,/)') '  End of Demo2'
      END SUBROUTINE Demo_2

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

      
      SUBROUTINE Phil_Initial_Params
      
         mixing_length_alpha = 1.6D0
         overshoot_param_H = 0.12D0
         overshoot_param_He = 0.12D0
         wind_Eta = 0d0
         number_of_center_shells_to_mix = 0
         call_save_exp_info = .true.
      
      END SUBROUTINE Phil_Initial_Params


      SUBROUTINE Phil_Before_Evolve
         call Demo2_Before_Evolve
         head_cnt = 10; summary_cnt = 2
      END SUBROUTINE Phil_Before_Evolve


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
      
      SUBROUTINE Lars_Initial_Params
         call Demo2_Initial_Params  
         he_break_even_model = -1   
      END SUBROUTINE Lars_Initial_Params


      SUBROUTINE Lars_Before_Evolve
         call Demo2_Before_Evolve
      END SUBROUTINE Lars_Before_Evolve


      INTEGER FUNCTION Lars_Check_Model()
         ! for project with Lars and Elliot.
         USE ez_do_one_data
         USE star_constants
         integer :: flg
      
         flg = Do_One_Check_Model()
         IF (flg /= KEEP_GOING) then
            Lars_Check_Model = flg
            return
         end if
         Lars_Check_Model = KEEP_GOING
         if (model_Number == 10) call save_lars_model
      
         if (he_break_even_model < 0) then
            if (helium_ignition) then
               he_break_even_model = model_Number
               call save_lars_model
            end if
         else if (mod(model_Number-he_break_even_model,100) == 0) then
            call save_lars_model
         else if (center_He < 1d-5 .and. mod(model_Number-he_break_even_model,20) == 0) then
            call save_lars_model
         end if
      
         contains
      
         subroutine report_negL
      
            integer :: i
            double precision :: minL, minL_mass
         
            minL = 0
         
            do i = N_CNTR_shell, N_SURF_shell, -1
               if (SX(SX_L,i) < minL) then
                  minL = SX(SX_L,i)
                  minL_mass = SX(SX_M,i)
               end if
            end do
            if (minL < 0) write(*,*) model_number, minL, minL_mass
      
         end subroutine report_negL
      
         subroutine save_lars_model
      
            character (len=256) :: filename
            integer :: io_unit, i
            io_unit = 63
         
            if (model_Number < 1000) then
               write(filename,'(a,f3.1,a,i3,a)') 'lars_data/m_',  initial_Mass, '_mod_', model_number, '.data'
            else
               write(filename,'(a,f3.1,a,i4,a)') 'lars_data/m_',  initial_Mass, '_mod_', model_number, '.data'
            end if
            write(*,*) 'save to ', trim(filename)
            open(unit=io_unit,file=trim(filename))
            call EZ_Extras
         
            call report_negL
         
            write(io_unit,'(a8,4x,99a20)') 
     >         'mod#', 'star_Age', 'star_Mass', 'initial_Mass', 'initial_Z'
            write(io_unit,'(i8,4x,e20.6,99f20.6)') model_number, star_Age, star_Mass,
     >         initial_Mass, initial_Z
     
            write(io_unit,'(99a20)') 
     >            'mass He core', 'mass C core', 'mass O core',
     >            'radius He core', 'radius C core', 'radius O core'
            write(io_unit,'(99e20.6)') mass_He_Core, mass_C_Core, mass_O_Core,
     >         radius_He_Core, radius_C_Core, radius_O_Core
     
            write(io_unit,'(99(a11,3x))') 'm', 'r', 'L', 'rho', 'P', 'T', 'conv vel', 'eps H', 'eps He',
     >         'opacity', 'xhe', 'xc', 'xn', 'xo'
            do i = N_CNTR_shell, N_SURF_shell, -1
            write(io_unit,'(99e14.6)') SX(SX_M,i), SX(SX_R,i), SX(SX_L,i), 
     >           SX(SX_RHO,i), SX(SX_P,i), SX(SX_T,i), 
     >           SX(SX_CV,i), SX(SX_EPS_H,i), SX(SX_EPS_HE,i), SX(SX_OPACITY,i), 
     >           SX(SX_XHE,i), SX(SX_XC,i), SX(SX_XN,i), SX(SX_XO,i)
            end do
            close(io_unit)
      
         end subroutine save_lars_model
      
      END FUNCTION Lars_Check_Model


      END MODULE demo2
     
     
     
     
     
     
