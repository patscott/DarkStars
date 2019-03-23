        MODULE ez_log_utils
        USE star_data
        USE star_extras
        USE star_constants
        IMPLICIT NONE
        
        INTEGER, PARAMETER :: buff_sz = 4, buff_line_len = 1024
        TYPE log_info
            CHARACTER (LEN=strlen) :: filename
            INTEGER :: io_unit
            LOGICAL :: in_use(buff_sz)
            DOUBLE PRECISION :: ages(buff_sz)
            CHARACTER(buff_line_len) :: lines(buff_sz)
        END TYPE log_info
            
        INTEGER, PARAMETER :: io_basics=50, io_conv=51, io_masses=52, io_power=53, io_burn=54, io_pressure=55
        INTEGER, PARAMETER :: io_profiles=56, io_extra=57, io_HR=58, io_epsnuc=59, io_mags=60

        INTEGER, PARAMETER :: basics_log=1, conv_log=2, masses_log=3, power_log=4, burn_log=5, pressure_log=6
        INTEGER, PARAMETER :: profiles_log=7, extra_log=8, epsnuc_log=9, mags_log=10
        INTEGER, PARAMETER :: num_log_files = 10
        TYPE (log_info) :: log_array(num_log_files)
        
        INTEGER, PARAMETER :: buff_size_HR = 10000
        DOUBLE PRECISION :: last_HR_age, h_logL(buff_size_HR), h_logTs(buff_size_HR)
        DOUBLE PRECISION :: h_logRHO(buff_size_HR), h_logTc(buff_size_HR)
        DOUBLE PRECISION :: h_age(buff_size_HR), h_mass(buff_size_HR)
        INTEGER :: h_model_number(buff_size_HR)
        INTEGER :: history_len
        CHARACTER (LEN=strlen) :: HR_history_name
        
        CONTAINS
        
        SUBROUTINE Update_HR
        USE ez_do_one_data
        INTEGER :: ios, i, n
        IF ( history_len .GE. buff_size_HR ) RETURN
        IF ( history_len .EQ. 0 .OR. last_HR_age .LT. star_Age ) history_len = history_len + 1
        h_logL(history_len) = log_Luminosity
        h_logTs(history_len) = log_surface_Temp
        h_logRHO(history_len) = log_center_Density
        h_logTc(history_len) = log_center_Temp
        h_age(history_len) = star_Age
        h_mass(history_len) = star_Mass
        h_model_number(history_len) = model_Number
        last_HR_age = star_Age
        ios = 0
        HR_history_name = TRIM(data_dir)//'/HR_history.log'
        ! Moved HR_history to logs directory for simultaneous runs - Pat Scott 20070905
        OPEN(UNIT=io_HR, FILE=TRIM(HR_history_name), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
        IF (ios /= 0) THEN
            WRITE(*,'(3a,i4,a,i4)') 'Update_HR: Failed to open ', TRIM(HR_history_name), '     ios =', ios, '    io_num =', io_HR
            RETURN
        END IF
        DO i = 1, history_len
           n = h_model_number(i)
           WRITE(io_HR,'(I6,1x,9(E18.12,1x))') n, h_age(i), h_mass(i), h_logL(i), h_logTs(i), h_logRHO(i), h_logTc(i)
        END DO
        CLOSE(io_HR)
        END SUBROUTINE Update_HR
        
        SUBROUTINE Write_Basics_Log
        INTEGER, PARAMETER :: num_out = 18
        DOUBLE PRECISION :: out(num_out)
        CHARACTER (LEN=buff_line_len) :: log_string
        out(1:7) = (/ star_Age, time_Step, log_Luminosity, log_Radius, log_surface_Temp, log_center_Temp, log_center_Density /)
        out(8:15) = (/ log_center_Pressure, center_Degeneracy, center_H, center_He, center_C, center_N, center_O, center_Ne /)
        out(16:num_out) = (/ dynamic_Timescale, KH_Timescale, nuc_Timescale /)
        WRITE(log_string,'(999(E18.12,1x))') out
        CALL Add_To_Log_Buffer(basics_log, star_Age, log_string)
        CALL Update_HR
        END SUBROUTINE Write_Basics_Log
        
        SUBROUTINE Write_Conv_log
        INTEGER, PARAMETER :: num_out = 20, null_zone = -20, num_zones = 3
        DOUBLE PRECISION :: out(num_out), CZR_s(num_zones), CZR_e(num_zones), CZM_s(num_zones), CZM_e(num_zones), CZT(num_zones)
        DOUBLE PRECISION :: inner_CZM, inner_CZR, inner_MZM, inner_MZR
        CHARACTER (LEN=buff_line_len) :: log_string
        INTEGER :: i, zone_i
        inner_CZM = SX(SX_M,SX_SURF); inner_CZR = SX(SX_R,SX_SURF) !PS DarkStars 080416
        inner_MZM = SX(SX_M,SX_SURF); inner_MZR = SX(SX_R,SX_SURF) !PS DarkStars 080416
        CZR_s = null_zone; CZR_e = null_zone; CZM_s = null_zone; CZM_e = null_zone; CZT = 0D0
        zone_i = 1
        DO i = 1, conv_boundary_count
            IF ( conv_turnover_Time(i) .GT. 0D0 ) THEN
                IF ( i .EQ. 1 ) THEN
                    CZR_s(zone_i) = -1D2
                    CZM_s(zone_i) = 0D0
                ELSE
                    CZR_s(zone_i) = log10(conv_boundaries_R(i-1) + 1D-99)
                    CZM_s(zone_i) = conv_boundaries_M(i-1)
                END IF
                CZR_e(zone_i) = log10(conv_boundaries_R(i) + 1D-99)
                CZM_e(zone_i) = conv_boundaries_M(i)
                CZT(zone_i) = conv_turnover_Time(i)
                IF ( zone_i .EQ. num_zones ) EXIT
                zone_i = zone_i + 1
            END IF
        END DO
        out(1:7) = (/ star_Age, CZR_s(1), CZR_e(1), CZR_s(2), CZR_e(2), CZR_s(3), CZR_e(3) /)
        out(8:16) = (/ CZM_s(1), CZM_e(1), CZM_s(2), CZM_e(2), CZM_s(3), CZM_e(3), CZT(1), CZT(2), CZT(3) /)
        DO i = SX_CNTR, SX_SURF, -1
            IF ( SX(SX_CV,i) .EQ. 0D0 ) THEN
                inner_CZM = SX(SX_M,i); inner_CZR = SX(SX_R,i); EXIT
            END IF
        END DO
        DO i = SX_CNTR, SX_SURF, -1
            IF ( SX(SX_SG,i) .EQ. 0D0 ) THEN
                inner_MZM = SX(SX_M,i); inner_MZR = SX(SX_R,i); EXIT
            END IF
        END DO
        out(17:num_out) = (/ inner_CZM, inner_CZR, inner_MZM, inner_MZR /)
        WRITE(log_string,'(999(E18.12,1x))') out
        CALL Add_To_Log_Buffer(conv_log, star_Age, log_string)
        END SUBROUTINE Write_Conv_log
        
        SUBROUTINE Write_Masses_Log
        INTEGER :: num_out
        DOUBLE PRECISION :: out(100)
        INTEGER :: IK
        CHARACTER (LEN=buff_line_len) :: log_string
        out(1:8) = (/ star_Age, star_Mass, star_Mdot, star_Mass_H, star_Mass_He, star_Mass_C, star_Mass_N, star_Mass_O /)
        out(9:12) = (/ star_Mass_Ne, mass_He_Core, mass_C_Core, mass_O_Core /)
        num_out = 12
        DO IK=SX_SURF+1, SX_CNTR, 10
           num_out = num_out + 1
           out(num_out) = SX(SX_M,IK)
        END DO
        WRITE(log_string,'(999(E18.12,1x))') out(1:num_out)
        CALL Add_To_Log_Buffer(masses_log, star_Age, log_string)
        END SUBROUTINE Write_Masses_Log
        
        SUBROUTINE Write_Pressure_Log
        INTEGER, PARAMETER :: num_out = 5
        DOUBLE PRECISION :: out(num_out)
        CHARACTER (LEN=buff_line_len) :: log_string
        out(1:num_out) = (/ star_Age, SX(SX_PRAD,SX_CNTR), SX(SX_PEL,SX_CNTR), SX(SX_PION,SX_CNTR), SX(SX_PCORR,SX_CNTR) /)
        WRITE(log_string,'(999(E18.12,1x))') out
        CALL Add_To_Log_Buffer(pressure_log, star_Age, log_string)
        END SUBROUTINE Write_Pressure_Log
        
        SUBROUTINE Write_Power_Log
        INTEGER, PARAMETER :: num_out = 12
        DOUBLE PRECISION :: out(num_out)
        CHARACTER (LEN=buff_line_len) :: log_string
        out(1:7) = (/ star_Age, power_H_burn, power_He_burn, power_Metal_burn, power_PP, power_CNO, power_3_alpha /)
        out(8:12)=(/ power_Neutrinos, power_plasmon_neutrinos, power_brem_neutrinos, power_pair_neutrinos, power_photo_neutrinos /)
        WRITE(log_string,'(999(E18.12,1x))') out
        CALL Add_To_Log_Buffer(power_log, star_Age, log_string)
        END SUBROUTINE Write_Power_Log
      
        SUBROUTINE Write_EPS_NUC_Log
        !        star_Age, NucBZM_1, NucBZM_2, NucBZM_3, NucBZM_4,...
        !            NucBZM_5, NucBZM_6, NucBZM_7, NucBZM_8,
        !            NucBZM_9, NucBZM_10, NucBZM_11, NucBZM_12
        INTEGER, PARAMETER :: num_out = 13
        DOUBLE PRECISION :: out(num_out), NucBZM_1, NucBZM_2, NucBZM_3, NucBZM_4
        DOUBLE PRECISION :: NucBZM_5, NucBZM_6, NucBZM_7, NucBZM_8
        DOUBLE PRECISION :: NucBZM_9, NucBZM_10, NucBZM_11, NucBZM_12
        CHARACTER (LEN=buff_line_len) :: log_string
        INTEGER :: I
        ! up to 3 zones where EPS_NUC > 1 erg/g/s
        ! for each zone have 4 numbers: start1, start2, end2, end1
        ! start1 is mass of inner edge where first goes > 1 (or -20 if none such)
        ! start2 is mass of inner edge where first zone reaches 10^3 erg/g/sec (or -20 if none such)
        ! end2 is mass of outer edge where first zone drops back below 10^3 erg/g/s
        ! end1 is mass of outer edge where first zone ends (i.e. EPS_NUC < 1)
        ! similar for second and third zones
        I = SX_CNTR
        CALL Find_EPS_NUC_Zone ( I, NucBZM_1, NucBZM_2, NucBZM_3, NucBZM_4 )
        CALL Find_EPS_NUC_Zone ( I, NucBZM_5, NucBZM_6, NucBZM_7, NucBZM_8 )
        CALL Find_EPS_NUC_Zone ( I, NucBZM_9, NucBZM_10, NucBZM_11, NucBZM_12 )
        out(1:5) = (/ star_Age, NucBZM_1, NucBZM_2, NucBZM_3, NucBZM_4 /)
        out(6:9) = (/ NucBZM_5, NucBZM_6, NucBZM_7, NucBZM_8 /)
        out(10:num_out) = (/ NucBZM_9, NucBZM_10, NucBZM_11, NucBZM_12 /)
        WRITE(log_string,'(999(E18.12,1x))') out
        CALL Add_To_Log_Buffer(epsnuc_log, star_Age, log_string)
        END SUBROUTINE Write_EPS_NUC_Log
        
         SUBROUTINE Find_EPS_NUC_Zone ( I_START, BZM_1, BZM_2, BZM_3, BZM_4 )
         USE ez_utils
         INTEGER, INTENT(INOUT) :: I_START
         DOUBLE PRECISION, INTENT(OUT) :: BZM_1, BZM_2, BZM_3, BZM_4
         DOUBLE PRECISION, PARAMETER :: null_zone = -20
         INTEGER :: I, burn_zone, I_X
         DOUBLE PRECISION :: burn_min1, burn_min2, prev_M, prev_X, cur_M, cur_X
         I_X = SX_EPS_NUC
         BZM_1 = null_zone; BZM_2 = null_zone; BZM_3 = null_zone; BZM_4 = null_zone
         burn_zone = 0 ! haven't entered the zone yet
         burn_min1 = 1D0; burn_min2 = 1D3
         IF (I_START .NE. SX_CNTR) THEN
            I = I_START+1
            prev_M = SX(SX_M,I)
            prev_X = SX(I_X,I)
         END IF
         DO I = I_START, SX_SURF, -1
            cur_M = SX(SX_M,I)
            cur_X = SX(I_X,I)
            SELECT CASE (burn_zone)
            CASE (0)
               IF ( SX(I_X,I) .GT. burn_min2 ) THEN
                  IF ( I .EQ. SX_CNTR ) THEN ! use star center as start of zone
                        BZM_2 = 0D0
                  ELSE ! interpolate to estimate where rate reached burn_min1 
                        BZM_2 = FIND0(prev_M,prev_X-burn_min2,cur_M,cur_X-burn_min2)
                  END IF
                  BZM_1 = BZM_2
                  burn_zone = 2
               ELSEIF ( SX(I_X,I) .GT. burn_min1 ) THEN
                  IF ( I .EQ. SX_CNTR ) THEN ! use star center as start of zone
                        BZM_1 = 0D0
                  ELSE ! interpolate to estimate where rate reached burn_min1 
                        BZM_1 = FIND0(prev_M,prev_X-burn_min1,cur_M,cur_X-burn_min1)
                  END IF
                  burn_zone = 1
               END IF
            CASE (1) ! in the initial eps > 1 region
               IF ( SX(I_X,I) .GT. burn_min2 ) THEN
                  BZM_2 = FIND0(prev_M,prev_X-burn_min2,cur_M,cur_X-burn_min2)
                  burn_zone = 2
               ELSEIF ( SX(I_X,I) .LT. burn_min1 ) THEN
                  BZM_4 = FIND0(prev_M,prev_X-burn_min1,cur_M,cur_X-burn_min1)
                  I_START = I
                  RETURN
               END IF
            CASE (2) ! in the initial eps > 1000 region
               IF ( SX(I_X,I) .LT. burn_min1 ) THEN
                  BZM_4 = FIND0(prev_M,prev_X-burn_min1,cur_M,cur_X-burn_min1)
                  BZM_3 = BZM_4
                  I_START = I
                  RETURN
               END IF
               IF ( SX(I_X,I) .LT. burn_min2 ) THEN
                  BZM_3 = FIND0(prev_M,prev_X-burn_min2,cur_M,cur_X-burn_min2)
                  burn_zone = 3
               END IF
            CASE (3) ! in the final eps > 1 region
               IF ( SX(I_X,I) .LT. burn_min1 ) THEN
                  BZM_4 = FIND0(prev_M,prev_X-burn_min1,cur_M,cur_X-burn_min1)
                  I_START = I
                  RETURN
               END IF
            CASE DEFAULT
               STOP 'Error in Find_EPS_NUC_Zone'
            END SELECT
            prev_M = cur_M; prev_X = cur_X
         END DO
         I_START = I
         SELECT CASE (burn_zone)
         CASE (0)
            RETURN
         CASE (1)
            BZM_4 = cur_M
         CASE (2)
            BZM_3 = cur_M
            BZM_4 = cur_M
         CASE (3)
            BZM_4 = cur_M
         CASE DEFAULT
            STOP 'Error in Find_EPS_NUC_Zone'
         END SELECT
         END SUBROUTINE Find_EPS_NUC_Zone

        SUBROUTINE Write_Burn_Log
        !        star_Age, HBZR_1, HBZR_2, HBZR_3, HBZR_4, HBZM_1, HBZM_2, HBZM_3, HBZM_4, ...
        !                   HeBZR_1, HeBZR_2, HeBZR_3, HeBZR_4, HeBZM_1, HeBZM_2, HeBZM_3, HeBZM_4
        INTEGER, PARAMETER :: num_out = 17
        DOUBLE PRECISION :: out(num_out), HBZR_1, HBZR_2, HBZR_3, HBZR_4, HBZM_1, HBZM_2, HBZM_3, HBZM_4
        DOUBLE PRECISION :: HeBZR_1, HeBZR_2, HeBZR_3, HeBZR_4, HeBZM_1, HeBZM_2, HeBZM_3, HeBZM_4
        CHARACTER (LEN=buff_line_len) :: log_string
        CALL Find_Burn_Zones ( SX_EPS_H, HBZR_1, HBZR_2, HBZR_3, HBZR_4, HBZM_1, HBZM_2, HBZM_3, HBZM_4 )
        CALL Find_Burn_Zones ( SX_EPS_HE, HeBZR_1, HeBZR_2, HeBZR_3, HeBZR_4, HeBZM_1, HeBZM_2, HeBZM_3, HeBZM_4 )
        out(1:9) = (/ star_Age, HBZR_1, HBZR_2, HBZR_3, HBZR_4, HBZM_1, HBZM_2, HBZM_3, HBZM_4 /)
        out(10:num_out) = (/ HeBZR_1, HeBZR_2, HeBZR_3, HeBZR_4, HeBZM_1, HeBZM_2, HeBZM_3, HeBZM_4 /)
        WRITE(log_string,'(999(E18.12,1x))') out
        CALL Add_To_Log_Buffer(burn_log, star_Age, log_string)
        END SUBROUTINE Write_Burn_Log
        
        SUBROUTINE Find_Burn_Zones ( I_X, BZR_1, BZR_2, BZR_3, BZR_4, BZM_1, BZM_2, BZM_3, BZM_4 )
        USE ez_utils
        INTEGER, INTENT(IN) :: I_X
        DOUBLE PRECISION, INTENT(OUT) :: BZR_1, BZR_2, BZR_3, BZR_4, BZM_1, BZM_2, BZM_3, BZM_4
        LOGICAL :: need_check 
        DOUBLE PRECISION, PARAMETER :: null_zone = -20
        INTEGER :: I, burn_zone
        DOUBLE PRECISION :: burn_min1, burn_min2, prev_M, prev_R, prev_X, cur_M, cur_R, cur_X
        BZR_1 = null_zone; BZR_2 = null_zone; BZR_3 = null_zone; BZR_4 = null_zone
        BZM_1 = null_zone; BZM_2 = null_zone; BZM_3 = null_zone; BZM_4 = null_zone
        burn_zone = 0; burn_min1 = 1D0; burn_min2 = 1D3
        DO I = SX_CNTR, SX_SURF, -1
            cur_M = SX(SX_M,I)
            cur_R = SX(SX_R,I)
            cur_X = SX(I_X,I)
            need_check = .TRUE.
            DO WHILE (need_check)
                need_check = .FALSE.
                SELECT CASE (burn_zone)
                CASE (0)
                    IF ( SX(I_X,I) .GT. burn_min2 ) THEN
                        BZM_1 = 0D0; BZR_1 = 0D0
                        IF ( I .EQ. SX_CNTR ) THEN ! use star center as start of zone
                            BZM_2 = 0D0; BZR_2 = 0D0
                        ELSE ! interpolate to estimate where rate reached burn_min1 
                            BZM_2 = FIND0(prev_M,prev_X-burn_min2,cur_M,cur_X-burn_min2)
                            BZR_2 = FIND0(prev_R,prev_X-burn_min2,cur_R,cur_X-burn_min2)
                        END IF
                        burn_zone = 2
                        need_check = .TRUE.
                    ELSEIF ( SX(I_X,I) .GT. burn_min1 ) THEN
                        IF ( I .EQ. SX_CNTR ) THEN ! use star center as start of zone
                            BZM_1 = 0D0; BZR_1 = 0D0
                        ELSE ! interpolate to estimate where rate reached burn_min1 
                            BZM_1 = FIND0(prev_M,prev_X-burn_min1,cur_M,cur_X-burn_min1)
                            BZR_1 = FIND0(prev_R,prev_X-burn_min1,cur_R,cur_X-burn_min1)
                        END IF
                        burn_zone = 1
                        need_check = .TRUE.
                    END IF
                CASE (1)
                    IF ( SX(I_X,I) .GT. burn_min2 ) THEN
                        BZM_2 = FIND0(prev_M,prev_X-burn_min2,cur_M,cur_X-burn_min2)
                        BZR_2 = FIND0(prev_R,prev_X-burn_min2,cur_R,cur_X-burn_min2)
                        burn_zone = 2
                        need_check = .TRUE.
                    ELSEIF ( SX(I_X,I) .LT. burn_min1 ) THEN
                        BZM_2 = FIND0(prev_M,prev_X-burn_min1,cur_M,cur_X-burn_min1)
                        BZR_2 = FIND0(prev_R,prev_X-burn_min1,cur_R,cur_X-burn_min1)
                        RETURN
                    END IF
                CASE (2)
                    IF ( SX(I_X,I) .LT. burn_min2 ) THEN
                        BZM_3 = FIND0(prev_M,prev_X-burn_min2,cur_M,cur_X-burn_min2)
                        BZR_3 = FIND0(prev_R,prev_X-burn_min2,cur_R,cur_X-burn_min2)
                        burn_zone = 3
                        need_check = .TRUE.
                    END IF
                CASE (3)
                    IF ( SX(I_X,I) .LT. burn_min1 ) THEN
                        BZM_4 = FIND0(prev_M,prev_X-burn_min1,cur_M,cur_X-burn_min1)
                        BZR_4 = FIND0(prev_R,prev_X-burn_min1,cur_R,cur_X-burn_min1)
                        RETURN
                    END IF
                CASE DEFAULT
                    STOP 'Error in Find_Burn_Zones'
                END SELECT
            END DO
            prev_M = cur_M; prev_R = cur_R; prev_X = cur_X
        END DO
        SELECT CASE ( burn_zone )
        CASE (0)
        CASE (1)
            BZM_2 = cur_M
            BZR_2 = cur_R
        CASE (2)
            BZM_3 = cur_M
            BZR_3 = cur_R
            BZM_4 = BZM_3
            BZR_4 = BZR_3
        CASE (3)
            BZM_4 = cur_M
            BZR_4 = cur_R
        CASE DEFAULT
            STOP 'Error in Find_Burn_Zones'
        END SELECT
        END SUBROUTINE Find_Burn_Zones
                
        SUBROUTINE Write_Mags_Log
        USE ez_magnitude_data
        INTEGER, PARAMETER :: num_out = n_mags + 1
        DOUBLE PRECISION :: out(num_out)
        CHARACTER (LEN=buff_line_len) :: log_string
        out(1) = star_Age
        out(2:num_out) = color_magnitudes(1:n_mags)
        WRITE(log_string,'(999(E18.12,1x))') out
        CALL Add_To_Log_Buffer(mags_log, star_Age, log_string)
        END SUBROUTINE Write_Mags_Log

        SUBROUTINE Write_Profiles_Log
        CHARACTER (LEN=buff_line_len) :: log_string
        WRITE(log_string,'(2(E18.12,1x),I6)') star_Age, star_Mass, model_Number
        CALL Add_To_Log_Buffer(profiles_log, star_Age, log_string)
        END SUBROUTINE Write_Profiles_Log
        
        SUBROUTINE Add_To_Log_Buffer (log_num, age, line)
        INTEGER, INTENT(IN) :: log_num
        DOUBLE PRECISION, INTENT(IN) :: age
        CHARACTER(buff_line_len), INTENT(IN) :: line
        INTEGER :: i, free_buff, min_i
        DOUBLE PRECISION :: min_age
        DO i=1,buff_sz
            IF ( log_array(log_num)%in_use(i) .AND. log_array(log_num)%ages(i) .GE. age ) log_array(log_num)%in_use(i) = .FALSE.
        END DO
        free_buff = -1
        DO i=1,buff_sz
            IF ( .NOT. log_array(log_num)%in_use(i) ) THEN
                free_buff = i; EXIT
            END IF
        END DO
        IF ( free_buff .LT. 0 ) THEN
            min_i = 1; min_age = log_array(log_num)%ages(1)
            DO i=2,buff_sz
                IF ( log_array(log_num)%ages(i) .LT. min_age ) THEN
                    min_age = log_array(log_num)%ages(i); min_i = i
                END IF
            END DO
            WRITE(log_array(log_num)%io_unit, '(a)') trim(log_array(log_num)%lines(min_i))
            log_array(log_num)%in_use(min_i) = .FALSE.
            free_buff = min_i
        END IF
        i = free_buff
        log_array(log_num)%in_use(i) = .TRUE.
        log_array(log_num)%ages(i) = age
        log_array(log_num)%lines(i) = line
        END SUBROUTINE Add_To_Log_Buffer
        
        SUBROUTINE Flush_Log_Buffer (log_num)
        INTEGER, INTENT(IN) :: log_num
        INTEGER :: i, min_i, j
        DOUBLE PRECISION :: min_age
        DO i=1,buff_sz
            min_i = -1 
            DO j=1,buff_sz
                IF ( log_array(log_num)%in_use(j) ) THEN
                    IF ( min_i .LT. 0 .OR. log_array(log_num)%ages(j) .LT. min_age ) THEN
                        min_i = j; min_age = log_array(log_num)%ages(j)
                    END IF
                END IF
            END DO
            IF ( min_i .LT. 0 ) RETURN
            WRITE(log_array(log_num)%io_unit, '(a)') trim(log_array(log_num)%lines(min_i))
            log_array(log_num)%in_use(min_i) = .FALSE.
        END DO
        END SUBROUTINE Flush_Log_Buffer
        
        END MODULE ez_log_utils
      
      
      
      
      
