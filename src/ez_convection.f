      MODULE ez_convection ! find the diffusion coefficient for convective mixing
      USE star_constants
      USE star_controls
      IMPLICIT NONE 
      
      CONTAINS
      
      DOUBLE PRECISION FUNCTION Grad_Star(GRADR, GRADA, SCP, HP, T, XHI, alpha, CSOUND, WCV, WL, WC1)
         ! returns GRAD, the actual temperature gradient taking into account convection
      DOUBLE PRECISION, INTENT(IN) :: GRADR, GRADA, SCP, HP, T, XHI, alpha, CSOUND
         ! GRADR, radiative temperature gradient
         ! GRADA, adiabatic temperature gradient
         ! SCP, specific heat capacity at constant pressure
         ! HP, pressure scale height for mixing length calculation
         ! T, temperature
         ! alpha, parameter for mixing length
         ! XHI, thermal conductivity = 4 a c T^3/(3 kappa rho^2 c_P). (cm^2/sec)
      DOUBLE PRECISION, INTENT(OUT) :: WCV, WL, WC1
         ! WCV, convection velocity
         ! WL, convection velocity times mixing length. units (cm^2/sec) are same as thermal conductivity.
         ! WC1, ?
      DOUBLE PRECISION :: DG, S2, WC2, WC3, WC4
      DOUBLE PRECISION :: CBRT, VX
      CBRT(VX) = DEXP(DLOG(VX)*C3RD) ! CBRT(VX) is cube root of VX
      VX = 0 ! for the Absoft compiler.

      ! solve a cubic to find convection velocity
      DG = GRADR - GRADA ! DG > 0 is Schwarzschild criterion for convection
      S2 = GRADA*SCP*T
      WC1 = 0.5D0*S2*(alpha*alpha*HP/(9D0*XHI))**2
      WC2 = 546.75D0*WC1*DMAX1(0D0, DG) + 73D0
      WC3 = CBRT(WC2 + DSQRT(WC2*WC2 + 12167D0))
      WC4 = DMAX1(WC3 - 23D0/WC3 - 2D0, 0D0)

      WCV = WC4*XHI/(3D0*alpha*HP)

      IF (WCV .GT. convection_speed_limit*CSOUND) WCV = convection_speed_limit*CSOUND ! limit convection velocity
      WL = alpha*HP*WCV
      Grad_Star = GRADR - 4D0*HP*WCV**3/(alpha*S2*XHI) ! actual temperature gradient
      IF (Grad_Star .GT. grad_star_limit) Grad_Star = grad_star_limit ! to prevent density inversion caused by convection
         ! for ideal gas, create density inversion if Grad_Star > 1
      END FUNCTION Grad_Star
      
      DOUBLE PRECISION FUNCTION Mixing (PR,PG,APM,M,overshoot_param,GRADR,GRADA,FAC)
      ! returns SG, the diffusion coefficient for convective mixing
      ! APM is d_ln(P) / d_M = -(G M)/(4 Pi r^4 P)
      DOUBLE PRECISION, INTENT(IN) :: PR,PG,APM,M,overshoot_param,GRADR,GRADA,FAC
      DOUBLE PRECISION :: B, WOS, EG, SG
      DOUBLE PRECISION, PARAMETER :: CU = 0.1D0
      B = PR/PG
      WOS = (2.5D0 + B*(2D1 + 1.6D1*B))*(CU/DABS(APM*M) + 1D0)
      EG = GRADR - GRADA + overshoot_param/WOS
      IF (EG .LE. 0D0) THEN
         SG = 0D0
      ELSE
         SG = FAC*EG*EG
         IF (DABS(SG) .LT. 1D-30) SG = 0D0
      END IF
      Mixing = SG
      END FUNCTION Mixing
      
      END MODULE ez_convection

