! WIMP annihilation energy module for DarkStars
!
! Pat Scott, July 2007; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE DkStrs_annihilation
      use DkStrs_data
      use DkStrs_WIMPDens
      implicit none

      contains

      
      DOUBLE PRECISION FUNCTION E_annihilation(r,rho)
      ! Takes local height and density and returns local energy yield of
      ! WIMP annihilation 
      ! Note that this integrand differs from Eq. 2.22 in Scott et al
      ! (2009, MNRAS 394:82) by a factor of 1/2, since the local annihilation
      ! rate should be WIMPdens^2 * sigann / 2.
      ! Input:    r   radius [cm]
      !          rho  density [g/cm^3]
      ! Output:       energy generation rate from WIMP annihilation [ergs/g/s]
   
      double precision, intent(IN) :: r,rho

      E_annihilation = mx * 1.d3 * CMEV * sigann / rho * WIMPdens(r)**2
      E_annihilation = E_annihilation * (1.d0 - nu_losses(rho))

      return

      END FUNCTION E_annihilation


      DOUBLE PRECISION FUNCTION nu_losses(rho)
      ! Estimates what fraction of WIMP annihilation energy goes
      ! into neutrinos and free-streams away instead of contributing
      ! to stellar energy generation. Currently doesn't do
      ! anything except apply the factor read from the input pipe.  In
      ! principle, if this were to be made more serious many other input
      ! parameters would be required than just the local density.
      ! Input:  rho   local density [g/cm^3]
      ! Output:       percentage energy lost to neutrinos [dimensionless]

      double precision, intent(IN) :: rho

      nu_losses = nuLossFactor

      return

      END FUNCTION nu_losses


      END MODULE DkStrs_annihilation
