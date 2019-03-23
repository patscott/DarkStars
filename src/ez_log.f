        MODULE ez_log
        USE star_data
        USE star_extras
        USE star_constants
        USE ez_log_utils
        IMPLICIT NONE
        
        ! This module takes care of producing logs to record the evolution of the star.

        ! If you need info not included in the standard logs, call Write_Extra_Log.
              
        CHARACTER (LEN=strlen) :: full_name, partial_name
        
        CONTAINS
        
        SUBROUTINE Open_Logs(log_dir, init_M, Build_Fname, open_extras)
        CHARACTER (LEN=strlen), INTENT(IN) :: log_dir ! path to directory where logs will live
        DOUBLE PRECISION, INTENT(IN) :: init_M
        INTERFACE
            SUBROUTINE Build_Fname(mass)
                DOUBLE PRECISION, INTENT(IN) :: mass
            END SUBROUTINE Build_Fname
        END INTERFACE
        LOGICAL, INTENT(IN) :: open_extras
        ! Removed superflous setting of HR_history log filename - Pat Scott 20070905 
        history_len = 0
        ! logs for various information concerning the star
        log_array(basics_log)%filename = 'basics.log' 
            CALL Init_Log(log_dir, init_M, basics_log, io_basics, Build_Fname)
        log_array(conv_log)%filename = 'convection.log'
            CALL Init_Log(log_dir, init_M, conv_log, io_conv, Build_Fname)
        log_array(masses_log)%filename = 'masses.log'
            CALL Init_Log(log_dir, init_M, masses_log, io_masses, Build_Fname)
        log_array(power_log)%filename = 'power.log'
            CALL Init_Log(log_dir, init_M, power_log, io_power, Build_Fname)
        log_array(pressure_log)%filename = 'pressure.log'
            CALL Init_Log(log_dir, init_M, pressure_log, io_pressure, Build_Fname)
        log_array(burn_log)%filename = 'burn.log' 
            CALL Init_Log(log_dir, init_M, burn_log, io_burn, Build_Fname)
        log_array(epsnuc_log)%filename = 'epsnuc.log' 
            CALL Init_Log(log_dir, init_M, epsnuc_log, io_epsnuc, Build_Fname)
        log_array(mags_log)%filename = 'mags.log' 
            CALL Init_Log(log_dir, init_M, mags_log, io_mags, Build_Fname)
        ! logs for profiles
        log_array(profiles_log)%filename = 'profiles.log'
            CALL Init_Log(log_dir, init_M, profiles_log, io_profiles, Build_Fname)
        IF ( open_extras ) CALL Open_Extra(log_dir, init_M, Build_Fname)
        END SUBROUTINE Open_Logs
        
        SUBROUTINE Open_Extra(log_dir, init_M, Build_Fname)
        CHARACTER (LEN=strlen), INTENT(IN) :: log_dir ! path to directory where logs will live
        DOUBLE PRECISION, INTENT(IN) :: init_M
        INTERFACE
            SUBROUTINE Build_Fname(mass)
                DOUBLE PRECISION, INTENT(IN) :: mass
            END SUBROUTINE Build_Fname
        END INTERFACE
         log_array(extra_log)%filename = 'extra.log' 
         CALL Init_Log(log_dir, init_M, extra_log, io_extra, Build_Fname)
        END SUBROUTINE Open_Extra
        
        SUBROUTINE Write_Logs
        CALL Write_Basics_Log
        CALL Write_Conv_Log
        CALL Write_Masses_Log
        CALL Write_Power_Log
        CALL Write_Pressure_Log
        CALL Write_Burn_Log
        CALL Write_EPS_NUC_Log
        CALL Write_Mags_Log
        END SUBROUTINE Write_Logs
        
        SUBROUTINE Write_Profiles
        CALL Write_Profiles_Log
        END SUBROUTINE Write_Profiles

        SUBROUTINE Close_Logs ( close_extras )
        LOGICAL, INTENT(IN) :: close_extras
        INTEGER :: i, io_u
        DO i=1,num_log_files
            io_u = log_array(i)%io_unit
            IF ( io_u .NE. 0 .AND. ((io_u .NE. io_extra) .OR. close_extras) ) THEN
                  CALL Flush_Log_Buffer(i)
                  CLOSE(io_u)
            END IF
        END DO
        END SUBROUTINE Close_Logs
        
        SUBROUTINE Write_Extra_Log(num_out, out)
        INTEGER, INTENT(IN) :: num_out
        DOUBLE PRECISION, INTENT(IN) :: out(num_out)
        CHARACTER (LEN=buff_line_len) :: log_string
        WRITE(log_string,*) out
        CALL Add_To_Log_Buffer(extra_log, star_Age, log_string)
        END SUBROUTINE Write_Extra_Log
      
        SUBROUTINE Init_Log(data_dir, init_M, log_num, io_num, Build_Fname)
        CHARACTER (LEN=strlen), INTENT(IN) :: data_dir
        DOUBLE PRECISION, INTENT(IN) :: init_M
        INTEGER, INTENT(IN) :: log_num, io_num
        INTERFACE
            SUBROUTINE Build_Fname(mass)
                DOUBLE PRECISION, INTENT(IN) :: mass
            END SUBROUTINE Build_Fname
        END INTERFACE
        INTEGER :: ios
        log_array(log_num)%io_unit = io_num
        IF (io_num .EQ. 0) THEN
           WRITE(*,'(a,i4)') 'Bad arg to Init_Log: io_num is 0 for log_num', log_num
           RETURN
        END IF
        partial_name = log_array(log_num)%filename
        CALL Build_Fname(init_M)
        ios = 0
        IF (LEN(TRIM(data_dir)) .GT. 0) full_name = TRIM(data_dir) // '/' // TRIM(full_name)
        OPEN(UNIT=io_num, FILE=trim(full_name), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
        IF (ios /= 0) THEN
            WRITE(*,'(3a,i4,a,i4)') 'Init_Log: Failed to open file ', TRIM(full_name), '     ios =', ios, '    io_num =', io_num
            RETURN
        END IF
        log_array(log_num)%in_use = .FALSE.
        END SUBROUTINE Init_Log
        
        END MODULE ez_log
      
      
