      MODULE ez_state_data    ! equation of state data module
      USE star_constants
      IMPLICIT NONE

      DOUBLE PRECISION ::  PRESSI_CONST, T_CRIT, PSI_CRIT ! controls 'fade out' of PRESSI corrections
      DOUBLE PRECISION ::  PRESSI_CONST_SAV, T_CRIT_SAV, PSI_CRIT_SAV

      DOUBLE PRECISION ::  PSI ! degeneracy
      
      DOUBLE PRECISION :: P             ! pressure.
      DOUBLE PRECISION :: LNP         ! ln(pressure).
      DOUBLE PRECISION :: PR           ! radiation pressure.
      DOUBLE PRECISION :: PE           ! electron pressure.
      DOUBLE PRECISION :: PG           ! gas pressure.
      DOUBLE PRECISION :: P_ION            ! ideal ion gas pressure
      DOUBLE PRECISION :: PCORR          ! correction to pressure for non-ideal gas

      DOUBLE PRECISION :: RHO       ! density.
      DOUBLE PRECISION :: LNRHO      ! ln(density).
      
      DOUBLE PRECISION :: U             ! internal energy (ergs per gram).
      DOUBLE PRECISION :: SCP         ! specific heat capacity at constant pressure (ergs per gram per K).

      DOUBLE PRECISION :: S           ! entropy per gram.
      DOUBLE PRECISION :: SF           ! dS/dln(F).
      DOUBLE PRECISION :: ST          ! dS/dln(T).

      DOUBLE PRECISION :: GRADA      ! adiabatic temperature gradient, dln(T)/dln(P).
      DOUBLE PRECISION :: GAMMA1    ! adiabatic exponent, dln(P)/dln(RHO).
      
      DOUBLE PRECISION :: GAM    ! GAM is the plasma interaction parameter.

      DOUBLE PRECISION :: NE           ! 
      DOUBLE PRECISION :: NE1         ! 
      DOUBLE PRECISION :: NI           ! 
      DOUBLE PRECISION :: NZZ         ! 
      DOUBLE PRECISION :: AVM           ! 
      
      DOUBLE PRECISION :: ZT
      DOUBLE PRECISION :: NA(NEL)

        ! The Fermi-Dirac integral approximation information.
      ! (Eggleton, Faulkner  Flannery; 1973) and (Pols, Tout, Eggleton, and Han; 1995).
      DOUBLE PRECISION ::  SE, RE, QE, RET, PET, QET, REF, PEF, QEF, XTT, XFT, XFF, XTTT, XFTT, XFFT, XFFF
      
      DOUBLE PRECISION :: xTGR(34), xGGR(13), xTAB(34,13,5) ! data for converting from luminosity and temperature to B-V and U-B.
      
      INTEGER, PARAMETER :: nzgr = 9, ntgr = 61, nggr = 11
      DOUBLE PRECISION :: zgr(nzgr), tgr(ntgr), ggr(nggr), ubv(nzgr,ntgr,nggr,5)
      DOUBLE PRECISION :: zgr_SAV(nzgr), tgr_SAV(ntgr), ggr_SAV(nggr), ubv_SAV(nzgr,ntgr,nggr,5)
      
      CONTAINS
      
      SUBROUTINE SAV_state_data( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) PRESSI_CONST, T_CRIT, PSI_CRIT, zgr, tgr, ggr, ubv
      ELSE
         WRITE (IO_UNIT) PRESSI_CONST, T_CRIT, PSI_CRIT, zgr, tgr, ggr, ubv
      END IF
      END SUBROUTINE SAV_state_data

      SUBROUTINE SAV_state_data_internal( read_flag )
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         PRESSI_CONST=PRESSI_CONST_SAV; T_CRIT=T_CRIT_SAV; PSI_CRIT=PSI_CRIT_SAV
		 zgr=zgr_SAV; tgr=tgr_SAV; ggr=ggr_SAV; ubv=ubv_SAV
      ELSE
         PRESSI_CONST_SAV=PRESSI_CONST; T_CRIT_SAV=T_CRIT; PSI_CRIT_SAV=PSI_CRIT
		 zgr_SAV=zgr; tgr_SAV=tgr; ggr_SAV=ggr; ubv_SAV=ubv
      END IF
      END SUBROUTINE SAV_state_data_internal

      END MODULE ez_state_data
