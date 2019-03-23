***********************************************************************
*** The halo velocity profile in the Maxwell-Boltzmann (Gaussian)
*** approximation.
*** input: velocity relative to earth [ km s^-1 ]
*** output: f(u) / u [ (km/ s)^(-2) ]
*** date: april 6, 1999
*** Modified: 2007-06-26
***********************************************************************

      real*8 function dscapudfgauss_star(u)
      implicit none
      include 'dscapstar.h'
      real*8 u

c...calculate f(u)/u
c...vd_3d_star is the 3D velocity dispersion (at the location of the star)
c...and v_star is the velocity of the star with respect to the non-rotating
c...halo.
      dscapudfgauss_star=
     &  sqrt(3./2.)/(vd_3d_star*v_star)/sqrt(3.1415)*(
     &  exp(-1.5*(u-v_star)**2/vd_3d_star**2)
     &  -exp(-1.5*(u+v_star)**2/vd_3d_star**2))

      return
      end
