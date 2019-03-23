!  experiment.f


      MODULE experiment
      USE ez_driver
      USE star_data
      USE star_extras
      USE star_controls
      USE star_constants
      USE ez_do_one_utils
      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: Lsun = 1d33 * CLSN
      DOUBLE PRECISION, PARAMETER :: Rsun = 1d11 * CRSN
      DOUBLE PRECISION, PARAMETER :: Msun = 1d33 * CMSN
      DOUBLE PRECISION, PARAMETER :: CSDAY = 86.4d3

      integer, parameter :: numout = 17
      integer, parameter :: max_model_num = 2000
      double precision :: output(max_model_num, numout), delta_lnP
      integer :: num_saved, model_out(max_model_num)

      integer :: io_unit
      
      CONTAINS
      
      
      INTEGER FUNCTION Experiment_Check_Model()
      IF (center_H < 0.1d0) THEN
         WRITE(*,*)
         WRITE(*,'(a,f8.3)') '   stop with center_H  =', center_H
         WRITE(*,'(a,f8.3)') '   stop with center_He =', center_He
         WRITE(*,*)
         Experiment_Check_Model = TERMINATE
      ELSE
         Experiment_Check_Model = KEEP_GOING
      END IF
      IF (model_Number > 10 .and. log10(star_Age) >= 8) CALL Write_Info
      END FUNCTION Experiment_Check_Model

      
      subroutine Write_Info
         use interp_1D
         integer :: k, i, j, ii, ibase, iprev, info, itop, k_top
         double precision :: t_therm, t_eddy, star_L, P, lnP, alfa, beta, lnP_base, lnP1, lnP2, m, cp, T, wl, cv,
     >      diff_at_km2, diff_at_km1, diff_at_k, diff_at_kp1, diff_below, diff_above, lnP_top, k_top_root,
     >      convel, mixlen, logAge, mass_above, logT, logSCP, kappa, rad_minus_ad
         
         double precision :: dk, k_root, x(6), f(4,6)
         double precision, parameter :: tiny = 1d-6
         
         call EZ_Extras
         ibase = -1
         do k = N_CNTR_shell-2, N_SURF_shell+2, -1
            ! search from center to surface for 1st time have 2 radiative followed by 2 convective
            if (SX(SX_GRAD_AD,k)   >= SX(SX_GRAD_RAD,k) .and. 
     >          SX(SX_GRAD_AD,k-1) <  SX(SX_GRAD_RAD,k-1)) then
               ibase = k-1; exit
            end if
         end do
         
         if (ibase < 0) then
            write(*,*) 'failed to find convective zone in model', model_Number
            return
         end if
         
         if (SX(SX_GRAD_AD,k) == SX(SX_GRAD_RAD,k)) then
            k_root = k; k = k + 1
         else 
            ! find place where grad_ad == grad_rad
            do j = 1, 6
               ii = k-4+j
               x(j) = ii
               f(1,j) = SX(SX_GRAD_AD,ii) - SX(SX_GRAD_RAD,ii)
            end do
            call mkpmcub_db(x,6,f)
            ! the cubic from k to k-1 should have the zero we want
            k_root = rt_mpcub_db(x, 3, 6, f, info)
            if (info /= 0) then
               stop 'failed to find root'
            end if
         end if
         
         ! transition from radiative to convective happens between k-1 and k
         
         if (k-1 > k_root .or. k_root > k) stop 'bad k_root for base'
         dk = k_root - (k-1)
         ibase = k
         
         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = log(SX(SX_P,ii))
         end do
         call mkpmcub_db(x,6,f)
         
         if (k_root < x(3) .or. k_root > x(4)) stop 'bad xs'
         lnP_base = f(1,3) + dk*f(2,3) + dk**2*f(3,3) + dk**3*f(4,3)
         
         ! move away from the edge of the zone until ln(P) has dropped by delta_lnP
         lnP = lnP_base - delta_lnP
         lnP1 = lnP_base
         iprev = -1
         do k = ibase, N_SURF_shell+1, -1
            lnP2 = log(SX(SX_P,k-1))
            if (lnP2 <= lnP) then
               iprev = k
               i = k-1
               exit
            end if
            lnP1 = lnP2
         end do
         
         if (iprev < 0) stop 'failed to find place where logP had dropped by 1'
         
         if (lnP2 == lnP) then
            k_root = k; k = k+1
         else 
            ! find place where ln(P) == lnP
            do j = 1, 6
               ii = k-4+j
               x(j) = ii
               f(1,j) = log(SX(SX_P,ii)) - lnP
            end do
            call mkpmcub_db(x,6,f)
            ! the cubic from k to k-1 should have the zero we want
            k_root = rt_mpcub_db(x, 3, 6, f, info)
            if (info /= 0) then
               stop 'failed to find root'
            end if
         end if
         
         if (k-1 > k_root .or. k_root > k) stop 'bad k_root for lnP'
         dk = k_root - (k-1)
      
         star_L = SX(SX_L,N_SURF_shell)
         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = (star_Mass - SX(SX_M,ii))*Msun * SX(SX_SCP,ii) * SX(SX_T,ii) / (CSDAY * Lsun * star_L)
         end do
         call mkpmcub_db(x,6,f)
         t_therm = f(1,3) + dk*f(2,3) + dk**2*f(3,3) + dk**3*f(4,3)

         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = SX(SX_WL,ii)
         end do
         call mkpmcub_db(x,6,f)
         WL = f(1,3) + dk*f(2,3) + dk**2*f(3,3) + dk**3*f(4,3)

         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = SX(SX_CV,ii)
         end do
         call mkpmcub_db(x,6,f)
         convel = f(1,3) + dk * f(2,3) + dk**2 * f(3,3) + dk**3 * f(4,3)

         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = log10(SX(SX_T,ii))
         end do
         call mkpmcub_db(x,6,f)
         logT = f(1,3) + dk * f(2,3) + dk**2 * f(3,3) + dk**3 * f(4,3)

         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = log10(SX(SX_SCP,ii))
         end do
         call mkpmcub_db(x,6,f)
         logSCP = f(1,3) + dk * f(2,3) + dk**2 * f(3,3) + dk**3 * f(4,3)

         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = star_Mass - SX(SX_M,ii)
         end do
         call mkpmcub_db(x,6,f)
         mass_above = f(1,3) + dk * f(2,3) + dk**2 * f(3,3) + dk**3 * f(4,3)

         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = SX(SX_OPACITY,ii)
         end do
         call mkpmcub_db(x,6,f)
         kappa = f(1,3) + dk * f(2,3) + dk**2 * f(3,3) + dk**3 * f(4,3)

         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = SX(SX_GRAD_RAD,ii) - SX(SX_GRAD_AD,ii)
         end do
         call mkpmcub_db(x,6,f)
         rad_minus_ad = f(1,3) + dk * f(2,3) + dk**2 * f(3,3) + dk**3 * f(4,3)
         
         t_eddy = WL / (CSDAY * convel**2)
         
         logAge = log10(star_Age)
         do while (num_saved > 0 .and. logAge <= output(num_saved,1))
            num_saved = num_saved - 1
         end do
         if (num_saved == 0 .or. logAge > output(num_saved,1)) num_saved = num_saved + 1
         if (num_saved > max_model_num) stop 'make max_model_num larger and try again'

         
         ! search for top of this convection zone
         itop = -1
         do k = ibase, N_SURF_shell+2, -1
            ! search from center to surface for 1st time have 2 radiative followed by 2 convective
            if (SX(SX_GRAD_AD,k)   <  SX(SX_GRAD_RAD,k) .and. 
     >          SX(SX_GRAD_AD,k-1) >= SX(SX_GRAD_RAD,k-1)) then
               itop = k-1; exit
            end if
         end do
         
         if (itop < 0) then
            write(*,*) 'failed to find convective zone in model', model_Number
            return
         end if
         
         if (SX(SX_GRAD_AD,itop) == SX(SX_GRAD_RAD,itop)) then
            k_top = itop
         else 
            ! find place where grad_ad == grad_rad
            do j = 1, 6
               ii = k-4+j
               x(j) = ii
               f(1,j) = SX(SX_GRAD_AD,ii) - SX(SX_GRAD_RAD,ii)
            end do
            call mkpmcub_db(x,6,f)
            ! the cubic from k to k-1 should have the zero we want
            k_top_root = rt_mpcub_db(x, 3, 6, f, info)
            if (info /= 0) then
               stop 'failed to find root'
            end if
         end if
         
         ! transition from convective to radiative happens between k-1 and k
         
         if (k-1 > k_top_root .or. k_top_root > k) stop 'bad k_top for top of convection zone'
         dk = k_top_root - (k-1)
         ibase = k
         
         do j = 1, 6
            ii = k-4+j
            x(j) = ii
            f(1,j) = log(SX(SX_P,ii))
         end do
         call mkpmcub_db(x,6,f)
         
         if (k_top_root < x(3) .or. k_top_root > x(4)) stop 'bad xs'
         lnP_top = f(1,3) + dk*f(2,3) + dk**2*f(3,3) + dk**3*f(4,3)
         
               
         model_out(num_saved) = model_Number
         output(num_saved,1:numout) = (/ logAge, log10(max(1d-99,t_therm)), log10(max(1d-99,t_eddy)), 
     >         convel, WL/convel, mass_above, lnP/CLN, lnP-lnP_top, lnP_base-lnP_top, logT, logSCP, kappa, rad_minus_ad, k_root,
     >         log_Luminosity, log_Radius, log_surface_Temp /)
            
      end subroutine Write_Info
      
      
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
      END SUBROUTINE Experiment_Initial_Params
      
      DOUBLE PRECISION FUNCTION Get_Init_Mass()
      INTEGER :: ios, io_num
      CHARACTER (LEN=strlen) :: in_filename
      DOUBLE PRECISION :: init_mass
      ios = 0; io_num = 38
      in_filename = 'experiment_in.data'
      OPEN(UNIT=io_num, FILE=trim(in_filename), ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(*,*) 'Get_Init_Mass: Failed to open file ', TRIM(in_filename), '     ios =', ios, '    io_num =', io_num
         Get_Init_Mass = -1D0
         RETURN
      END IF
      READ (io_num,*) init_mass
      Get_Init_Mass = init_mass
      CLOSE(io_num)
      WRITE (*,'(A,F9.4)') '   Doing mass equal', init_mass
      END FUNCTION Get_Init_Mass
            
      SUBROUTINE Do_Experiment
      CHARACTER (LEN=strlen) :: metals_dir, fname
      DOUBLE PRECISION :: mass
      INTEGER :: ios, i
      metals_dir = '../metals/z02'  ! where to get metals data
      delta_lnP = 1
      CALL EZ_Start ( metals_dir )
      mass = Get_Init_Mass()
      io_unit = 40
      fname = 'plot.data'
      OPEN(UNIT=io_unit, FILE=TRIM(fname), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) STOP 'failed to open output file'
      num_saved = 0
      write(io_unit,'(2F6.2,A)') mass, delta_lnP, ' : mass and delta lnP from base of convection to measurements'
      write(io_unit,'(A16,A8,99A16)') 'log Age', 'model#', 'log t_therm', 'log t_eddy', 'conv vel', 'mixing len', 
     >         'mass above', 'log P', 'log T', 'log Cp', 'opacity', 'surf logL', 'surf logR', 'surf logT'      
      WRITE(*,*) '  Start of Experiment'
      IF (.NOT. EZ_ZAMS( mass, Experiment_Initial_Params )) RETURN
      CALL EZ_Evolve(Experiment_Check_Model)
      WRITE(*,*) '  End of Experiment'
      write(io_unit,*)
      do i = 1, num_saved
         write(io_unit,'(E16.6,I8,99E16.6)') output(i,1), model_out(i), output(i,2:numout)
      end do
      close(io_unit)
      fname = 'plot_names.data'
      OPEN(UNIT=io_unit, FILE=TRIM(fname), ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      write(io_unit,*) 'log star_Age'
      write(io_unit,*) 'model Number'
      write(io_unit,*) 'log t therm'
      write(io_unit,*) 'log t eddy'
      write(io_unit,*) 'conv vel'
      write(io_unit,*) 'mixing len'
      write(io_unit,*) 'mass above'
      write(io_unit,*) 'logP'
      write(io_unit,*) 'dlnP to zone top'
      write(io_unit,*) 'dlnP of zone'
      write(io_unit,*) 'logT'
      write(io_unit,*) 'log Cp'
      write(io_unit,*) 'opacity'
      write(io_unit,*) 'rad minus ad'
      write(io_unit,*) 'k_root'
      write(io_unit,*) 'surf logL'
      write(io_unit,*) 'surf logR'
      write(io_unit,*) 'surf logT'     
      close(io_unit)
      END SUBROUTINE Do_Experiment

      END MODULE experiment
      
