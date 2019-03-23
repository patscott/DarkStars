! Constants module for DarkStars
!
! Pat Scott, July 2007; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE DkStrs_constants
      use star_constants
      implicit none

      ! Set true to do a solar "legacy" capture calculation, using
      ! the DarkSUSY "nt" programs (these are included
      ! in DarkSUSYLE only for the purposes of legacy capture).
      ! Uses the standard solar model of Bachall et al to get
      ! the density structure of the Sun (though additional
      ! abundances still come from the chosen input abundance file).
      ! Doing only a 16-species capture (i.e. no He3, etc) is quite
      ! different in this case, as the Bachall model includes the
      ! helium-3 shell above the core; this is left out normally since
      ! EZ doesn't follow helium-3 explicitly.  Leaving out such a zone
      ! (in a star which posesses it) raises the effective metallicity
      ! of the star, causing the capture rate to be overestimated a
      ! tad for capture dominated by the spin-independent cross-section
      ! (e.g. about 7% in the Sun).
      logical, parameter :: legacy=.false.

      ! Version tag
      character (len=strlen) :: DarkStars_version = 'v1.06'

      ! grams per GeV
      double precision, parameter :: CGGEV = 1.d3 * CMEV / CL**2
      ! cm per parsec
      double precision, parameter :: CMPC = 3.0857d18

      ! Heavy element solar abundance indexes
      integer, parameter :: n_heavy_abuns=9
      integer, parameter :: SOL_Na=1
      integer, parameter :: SOL_Al=2
      integer, parameter :: SOL_Si=3
      integer, parameter :: SOL_S =4
      integer, parameter :: SOL_Ar=5
      integer, parameter :: SOL_Ca=6
      integer, parameter :: SOL_Fe=7
      integer, parameter :: SOL_Ni=8
      integer, parameter :: SOL_Pb=9

      ! Solar/solar system isotopic ratio indexes
      integer, parameter :: n_iso_ratios=6
      integer, parameter :: ISO_4He3    =1
      integer, parameter :: ISO_12C13   =2
      integer, parameter :: ISO_16O18   =3
      integer, parameter :: ISO_58Ni60  =4
      integer, parameter :: ISO_208Pb207=5
      integer, parameter :: ISO_208Pb206=6

      character (len=strlen) :: DarkStars_extension = '.dark'

      END MODULE DkStrs_constants
