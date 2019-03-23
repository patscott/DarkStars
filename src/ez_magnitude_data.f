      MODULE ez_magnitude_data
      IMPLICIT none
      
      ! color indices are differences inmagnitudes in different wavelength bands
      ! as a reminder for non-experts like myself, here's how it goes
      !
      ! msun := apparent magnitude of sun is -26.81
      ! Fsun := solar flux at 1AU is 1.36e6 erg/s/cm^2
      !
      ! "apparent magnitude" m of star with flux F is m = msun - 2.5 log10(F/Fsun)
      ! "absolute magnitude" M for star of apparent magnitude m at distance d is M = m - 5 log(d/d0)
      !     where the standard distance d0 is 10pc.
      !     i.e., absolute magnitude is what the apparent magnitude would be if star were at 10 parsecs.
      !
      ! thus absolute magnitude of sun is 4.746
      ! 
      ! "bolometric magnitude" = absolute magnitude using flux integrated over all wavelengths
      !     can be derived from the current stellar luminosity using the equation
      !     log(Lstar/Lsun) = (Mbol_sun - Mbol_star)/2.5 using Mbol_sun = 4.746 (LCB)
      !
      ! "visual magnitude" = absolute magnitude only using flux in visible wavelengths
      !      more precisely, this is magnitude as measured with filter centered at 5500A, 890A width.
      !
      ! "bolometric correction" = bolometric magnitude minus visual magnitude
      !      for the sun, the bolometric correction is -0.108
      !      thus visual magnitude of sun is 4.854 = Mbol_sun - BC_sun = 4.7846 - (-0.108)
      !
      ! in order of increasing wavelength, the "color" magnitudes are as follows:
      !
      ! "U" is the ultraviolet magnitude, filter center at 365nm.
      ! "B" is the        blue magnitude, filter center at 440nm.
      ! "V" is the      visual magnitude, filter center at 550nm.
      ! "R" is the         red magnitude, filter center at 600nm.
      ! "I" is the   infra-red magnitude, filter center at 800nm.
      
      ! in addition, longer wavelength "colors" have been defined as well
      ! by order of increasing wavelength, these are J, H, K, L, and M.
      
      ! "color index" is the difference between 2 color magnitudes
      ! for example, B-V is mag_B - mag_V
      ! smaller B-V means larger brightness in blue band compared to visual band, means bluer star.


      ! color magnitude data from Lejeune, Cuisinier, Buser (1998) A&AS 130, 65 [LCB]

      ! they use [Fe/H] as a parameter; we just use log10(Z/Zsun) as an approximation for this.
            
      
      integer, parameter :: bol  = 1      ! Bolometric magnitude
      integer, parameter :: bcv  = 2      ! Bolometric correction in V (energy-weighted)
      integer, parameter :: umb  = 3      ! (U-B) color (energy-weighted)
      integer, parameter :: bmv  = 4      ! (B-V) color
      integer, parameter :: vmr  = 5      ! (V-R) color
      integer, parameter :: vmi  = 6      ! (V-I) color
      integer, parameter :: vmk  = 7      ! (V-K) color
      integer, parameter :: rmi  = 8      ! (R-I) color
      integer, parameter :: imk  = 9      ! (I-K) color
      integer, parameter :: jmh  = 10     ! (J-H) color
      integer, parameter :: hmk  = 11     ! (H-K) color
      integer, parameter :: kml  = 12     ! (K-L) color
      integer, parameter :: jmk  = 13     ! (J-K) color
      integer, parameter :: jml  = 14     ! (J-L) color
      integer, parameter :: jmlp = 15     ! (J-L') color
      integer, parameter :: kmm  = 16     ! (K-M) color
      
      integer, parameter :: n_mags  = 16

      END MODULE ez_magnitude_data

