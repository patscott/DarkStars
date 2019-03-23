      MODULE demo1
      USE ez_driver
      USE star_data
      USE star_extras
      USE star_controls
      USE star_constants
      IMPLICIT NONE

      INTEGER, PARAMETER :: IO_LOG = 50 ! unit number for writing log file
      INTEGER, PARAMETER :: IO_DIFF = 51 ! unit number for writing diff file
      INTEGER, PARAMETER :: IO_DATA = 52 ! unit number for writing data file
      
      CHARACTER (LEN=strlen) :: metals_path, log_file, diff_file, data_file, Z_str
      
      CONTAINS

      INTEGER FUNCTION Demo1_Check_Model()
      IF ( model_Number .LE. 1 ) THEN
         Demo1_Check_Model = KEEP_GOING
      ELSE
         Demo1_Check_Model = TERMINATE
      END IF
      END FUNCTION Demo1_Check_Model
      
      SUBROUTINE Initial_Params
      overshoot_param_H = 0D0; overshoot_param_He = 0D0 ! turn off convective overshooting
      wind_Eta = 1D0 ! turn on Reimer's wind
      END SUBROUTINE Initial_Params
      
      SUBROUTINE Sample_ZAMS(ZAMS_increment)
      USE ez_do_one
      INTEGER, INTENT(IN) :: ZAMS_increment
      INTEGER, PARAMETER :: num_out = 16, num_step1 = 30, LOG_LEN = 200
      DOUBLE PRECISION :: LM, LM_last, LM_mid, LM_step1, LM_step2, out(num_out), init_M
      DOUBLE PRECISION :: cnv_Size, size_M_Conv, top_M_Conv
      DOUBLE PRECISION, PARAMETER :: Mmin = 0.1D0, Mmax = 100D0, Mmid = 3D0
      INTEGER :: i, j, n
      CHARACTER (LEN=LOG_LEN) :: log_str
      LOGICAL :: data_flg
      data_flg = ZAMS_increment .EQ. 1
      WRITE (*,'(//,2A)') '  Sample the zero age main sequence using data from ', metals_path
      CALL EZ_Start ( metals_path )
      OPEN(UNIT=IO_LOG, FILE=trim(log_file), ACTION='WRITE', STATUS='REPLACE')
      OPEN(UNIT=IO_DIFF, FILE=trim(diff_file), ACTION='WRITE', STATUS='REPLACE')
      IF (data_flg) OPEN(UNIT=IO_DATA, FILE=trim(data_file), ACTION='WRITE', STATUS='REPLACE')
      LM = log10(Mmin); LM_last = log10(Mmax); LM_mid = log10(Mmid)
      LM_step1 = (LM_mid - LM) / num_step1;
      LM_step2 = (LM_last - LM_mid) / (num_step1 - 1);
      n = 2*num_step1
      IF (data_flg) WRITE(IO_DATA,*) n
      DO i=1,n,ZAMS_increment
         IF (MOD(i-1,10*ZAMS_increment) .EQ. 0) THEN
            log_str = '        Mass        log(L)    log(TNuc)   log(Tsurf)   log(RHO)     Tcenter   log(Pcntr)    log(R)'
            log_str = trim(log_str) // '      TopCNV      BotCNV   log(RHOsurf)    EPS PP'
            log_str = trim(log_str) // '     EPS CNO   opac cntr   opac surf   Pr/P surf'
            WRITE (*,'(/,A)') trim(log_str)
         END IF
         IF (i .EQ.1) THEN
            init_M = Mmin
         ELSE IF (i .EQ. num_step1+1) THEN
            init_M = Mmid
         ELSE IF (i .EQ. n) THEN
            init_M = Mmax
         ELSE
            init_M = 10**LM
         END IF
         IF (.NOT. EZ_ZAMS(init_M, Initial_Params)) EXIT
         CALL EZ_Evolve ( Demo1_Check_Model ) ! evolve for a step to let things settle in before output values
         CALL EZ_Extras
         size_M_Conv = 0D0
         top_M_Conv = 0D0
         DO j = 1, BSZ
            IF ( conv_turnover_Time(j) .GT. 0D0 ) THEN
               IF ( j .EQ. 1 ) THEN
                  cnv_Size = conv_boundaries_M(j)
               ELSE
                  cnv_Size = conv_boundaries_M(j) - conv_boundaries_M(j-1)
               END IF
               IF ( cnv_Size .GT. size_M_Conv ) THEN
                  size_M_Conv = cnv_Size; top_M_Conv = conv_boundaries_M(j)
               END IF
            END IF
         END DO
         out(1) = star_Mass; out(2) = log_Luminosity; out(3) = LOG10(nuc_Timescale+1D-99)
         out(4) = log_surface_Temp; out(5) = log_center_Density; out(6) = 10.0**log_center_Temp*1e-7
         out(7) = log_center_Pressure; out(8) = log_Radius
         out(9) = top_M_Conv/star_Mass; out(10) = (top_M_Conv-size_M_Conv)/star_Mass; out(11) = LOG10(SX(SX_RHO,SX_SURF))         
         out(12) = SX(SX_EPS_PP,SX_CNTR); out(13) = SX(SX_EPS_CNO,SX_CNTR)
         out(14) = SX(SX_OPACITY,SX_CNTR)
         out(15) = SX(SX_OPACITY,SX_SURF)
         out(16) = SX(SX_PRAD,SX_SURF) / SX(SX_P,SX_SURF)
         WRITE(*,'(3X,11(F12.6),4(E12.3),F12.6)') out(1:num_out)
         WRITE(IO_LOG,'(I10, 99E16.8)') i, out
         IF (data_flg) WRITE(IO_DATA,*) log_center_Density, log_center_Temp, log_Luminosity, log_surface_Temp, star_Mass
         WRITE(IO_DIFF,'(/,I9)') i
         CALL Write_DIFFME_File(IO_DIFF, num_out, out)
         IF ( i .LE. num_step1 ) THEN
            LM = LM + LM_step1*ZAMS_increment
         ELSE
            LM = LM + LM_step2*ZAMS_increment
         END IF
      END DO
      CLOSE ( IO_LOG )
      CLOSE ( IO_DIFF )
      IF (data_flg) CLOSE ( IO_DATA )
      END SUBROUTINE Sample_ZAMS
      
      SUBROUTINE Make_Filenames(num)
      INTEGER, INTENT(IN) :: num
      WRITE(metals_path,'(2A)') '../metals/', trim(Z_str)
      WRITE(log_file,'(3A,I1,1A)') 'demo1_', trim(Z_str), '_', num, '.log'
      WRITE(diff_file,'(3A,I1,1A)') 'demo1_', trim(Z_str), '_', num, '.DIFFME'
      IF (num .EQ. 1) WRITE(data_file,'(3A)') 'ZAMS_', trim(Z_str), '.data'
      END SUBROUTINE Make_Filenames

      SUBROUTINE Demo_1
      INTEGER, PARAMETER :: ZAMS_inc = 6 ! code samples (60/ZAMS_inc) models from 0.1 to 100 Msun.
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,'(/,A,/)') '  Demo1 computes zero age main sequences.'
      WRITE(*,'(a)') ' The mass, luminosity, and radius are in solar units.  The log''s are base 10.'
      WRITE(*,'(a)') ' TNuc is the nuclear burning time scale in years.  Tsurf is the surface temperature.'
      WRITE(*,'(a)') ' RHO is the central density.  Tcenter is the center temperature in 10^7 K.'
      WRITE(*,'(a)') ' TopCNV and BotCNV are the bounds of the main convection zone in units of fraction of total mass.'
      WRITE(*,'(a)') ' EPS-PP is central energy generation rate from the P-P chain (ergs/gm/sec).'
      WRITE(*,'(a)') ' EPS-CNO is central energy rate from the CNO cycle (ergs/gm/sec).'
      Z_str = 'z0001'; CALL Make_Filenames(ZAMS_inc); CALL Sample_ZAMS(ZAMS_inc)
      !Z_str = 'z0003'; CALL Make_Filenames(ZAMS_inc); CALL Sample_ZAMS(ZAMS_inc)
      !Z_str = 'z001'; CALL Make_Filenames(ZAMS_inc); CALL Sample_ZAMS(ZAMS_inc)
      Z_str = 'z004'; CALL Make_Filenames(ZAMS_inc); CALL Sample_ZAMS(ZAMS_inc)
      !Z_str = 'z01'; CALL Make_Filenames(ZAMS_inc); CALL Sample_ZAMS(ZAMS_inc)
      Z_str = 'z02'; CALL Make_Filenames(ZAMS_inc); CALL Sample_ZAMS(ZAMS_inc)
      !Z_str = 'z03'; CALL Make_Filenames(ZAMS_inc); CALL Sample_ZAMS(ZAMS_inc)
      WRITE (*,*)
      WRITE (*,*)
      END SUBROUTINE Demo_1
      
      END MODULE demo1
      
      
      
      