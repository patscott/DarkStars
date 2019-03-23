      MODULE star_constants
      IMPLICIT NONE
            
      ! Here are some useful constants
      
      ! Solar
      DOUBLE PRECISION, PARAMETER :: CLSN = 3.844D0 ! Lsolar*1e-33 erg/s
      DOUBLE PRECISION, PARAMETER :: CMSN = 1.9891D0 ! Msolar*1e-33 grams
      DOUBLE PRECISION, PARAMETER :: CRSN = 0.69598D0 ! Rsolar*1e-11 cm
      DOUBLE PRECISION, PARAMETER :: CASN = 4.57D9 ! solar age in years
      
      ! A note regarding units in the code:
      ! The items in star_data and star_extras are usually given in solar units
      ! However, the calculations in the body of the code are in "working units"
      ! which are as follows
      !     mass - 10^33 gm
      !     length - 10^11 cm
      !     luminosity - 10^33 erg
      
      ! Physical
      DOUBLE PRECISION, PARAMETER :: CL = 2.99792458D10 ! speed of light in cm/sec
      DOUBLE PRECISION, PARAMETER :: CG = 6.67259D-8 ! gravitational constant
      DOUBLE PRECISION, PARAMETER :: CSY = 3.155692597D7 ! seconds per year
      DOUBLE PRECISION, PARAMETER :: CMEV = 1.602176487D-6 ! ergs per MeV
      DOUBLE PRECISION, PARAMETER :: CMP = 1.672621637D-24 ! mass of proton
      DOUBLE PRECISION, PARAMETER :: AMU = 1.6605402D-24 ! grams per atomic mass unit 
      DOUBLE PRECISION, PARAMETER :: AME = 9.1093897D-28 ! grams per electron
      DOUBLE PRECISION, PARAMETER :: PLANCK = 6.6260755D-27 ! h
      DOUBLE PRECISION, PARAMETER :: BOLTZM = 1.380658D-16 ! kb
      DOUBLE PRECISION, PARAMETER :: ECHAR = 4.8032068D-10 ! electron charge
      DOUBLE PRECISION :: CRAD, CRAD_SAV ! radiation constant, a.
      DOUBLE PRECISION :: CB, CB_SAV ! pi me c^2 / lambdaC^3, lambdaC is Compton wavelength of electron.
      DOUBLE PRECISION :: CD, CD_SAV ! $8 \pi (m_e c/h)^3$ (grams/amu).
      DOUBLE PRECISION :: CR, CR_SAV ! gas constant (erg/K/mole).
      DOUBLE PRECISION :: CTE, CTE_SAV ! $k_B/(m_e c^2)$.
      DOUBLE PRECISION :: CEVB, CEVB_SAV ! (ergs per electron volt)/kb
      DOUBLE PRECISION :: CEN, CEN_SAV ! $h^3/(2 \pi k_B)^{3/2}$/(grams/amu)$^{5/2}$.
      DOUBLE PRECISION :: CPL, CPL_SAV ! $(4 \pi e^3)/(k_B^3$(grams/amu)), where e is electron charge.
      DOUBLE PRECISION :: CME, CME_SAV ! (ergs/MeV)/(grams/amu).
      DOUBLE PRECISION :: CG1, CG1_SAV ! 1e5 $G^{1/2}$.
      DOUBLE PRECISION :: CG2, CG2_SAV ! CG1 $0.432/ \pi$.
      DOUBLE PRECISION :: CGRT, CGRT_SAV ! $6.4 (\pi/4.32e4)^6 (1e11/c)^5$.

      ! Mathematical
      DOUBLE PRECISION :: CPI, CPI_SAV ! pi
      DOUBLE PRECISION :: C4PI, C4PI_SAV ! 4 pi
      DOUBLE PRECISION :: CLN, CLN_SAV ! ln(10); divide natural log value by CLN to convert to log10
      DOUBLE PRECISION :: C3RD, C3RD_SAV ! 1/3
      
      ! Miscellaneous
      INTEGER, PARAMETER :: strlen = 128
      
      ! Composition indices
      INTEGER, PARAMETER ::  N_H=1, N_HE=2, N_C=3, N_N=4, N_O=5, N_NE=6, N_MG=7, N_SI=8, N_FE=9, NEL=9
      
      CONTAINS
      
      SUBROUTINE SAV_star_constants( IO_UNIT, read_flag )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         READ (IO_UNIT) CRAD, CB, CD, CR, CTE, CEVB, CEN, CPL, CME, CG1, CG2, CGRT, CPI, C4PI, CLN, C3RD
      ELSE
         WRITE (IO_UNIT) CRAD, CB, CD, CR, CTE, CEVB, CEN, CPL, CME, CG1, CG2, CGRT, CPI, C4PI, CLN, C3RD
      END IF
      END SUBROUTINE SAV_star_constants

      SUBROUTINE SAV_star_constants_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         CRAD=CRAD_SAV; CB=CB_SAV; CD=CD_SAV; CR=CR_SAV; CTE=CTE_SAV; CEVB=CEVB_SAV
		 CEN=CEN_SAV; CPL=CPL_SAV; CME=CME_SAV; CG1=CG1_SAV; CG2=CG2_SAV; CGRT=CGRT_SAV
		 CPI=CPI_SAV; C4PI=C4PI_SAV; CLN=CLN_SAV; C3RD=C3RD_SAV
      ELSE
         CRAD_SAV=CRAD; CB_SAV=CB; CD_SAV=CD; CR_SAV=CR; CTE_SAV=CTE; CEVB_SAV=CEVB
		 CEN_SAV=CEN; CPL_SAV=CPL; CME_SAV=CME; CG1_SAV=CG1; CG2_SAV=CG2; CGRT_SAV=CGRT
		 CPI_SAV=CPI; C4PI_SAV=C4PI; CLN_SAV=CLN; C3RD_SAV=C3RD
      END IF
      END SUBROUTINE SAV_star_constants_internal

      END MODULE star_constants

