***********************************************************************
*** dscapfoveru returns the value of f(u)/u where f(u) is the 
*** velocity distribution (at infinity) and u is the velocity
***
*** Input: u = velocity (km/s) relative to the star
*** Output: f(u) / u [ cm^-3 (cm/s)^(-2)
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2007-06-26
***********************************************************************

      real*8 function dscapfoveru(u)
      implicit none
      include 'dscapstar.h'
      real*8 u
      real*8 dscapudfgauss_star, dshmuDFnumc
      external dscapudfgauss_star 
      external dshmuDFnumc

      dscapfoveru=dscapudfgauss_star(u)
     &  *1.0d-10            ! (km/s)^(-2) -> (cm/s)^(-2)
     &  *(rhowimp/capmx)        
     
      return
      end


***********************************************************************
*** gal_foverusq returns f(u)/u^2 for an isothermal halo profile in the
*** galactic rest frame, not in the frame of the star like dscapfoveru!
***
*** Input: u = velocity (km/s) relative to the galactic rest frame
*** Output: f(u) / u^2 [cm^(-3) * (cm/s)^(-3)]
***
*** Author: Pat Scott, pat@physto.se
*** Date: 2008-03-18
***********************************************************************

      real*8 function gal_foverusq(u)
      implicit none
      include 'dscapstar.h'
      real*8 u, A

      A = 3.d0/2.d0/vd_3d_star**2
      gal_foverusq = 2.d0 / sqrt(pi) * A**(1.5d0) *
     &  exp(-A*u*u)
     &  * 1.0d-15            ! (km/s)^(-3) -> (cm/s)^(-3)
     &  * rhowimp/capmx      ! local number density of WIMPs
      
      return
      end


***********************************************************************
*** alt_foverusq is like gal_foverusq, but should point to some
*** other fancy velocity distribution which is not just a normal
*** isothermal Gaussian.
***
*** Input: u = velocity (km/s) relative to the galactic rest frame
*** Output: f(u) / u^2 [cm^(-3) * (cm/s)^(-3)]
***
*** Author: Pat Scott, pat@physto.se
*** Date: 2008-03-18
***********************************************************************

      real*8 function alt_foverusq(u)
      implicit none
      include 'dscapstar.h'
      real*8 u, malc_foverusq
      external malc_foverusq

      alt_foverusq = malc_foverusq(u) !hit up Malcolm's Tsalis-smashing vel distro
     
      return
      end

***********************************************************************
*** alt_f simply returns the alt f(u) for calculating the
*** correct renormalisation when truncating an alternative velocity
*** distribution.  Also removes the scaling for the local number 
*** density of WIMPs.
***
*** Input: u = velocity (km/s) relative to the galactic rest frame
*** Output: f(u) (km/s)^(-1)
***
*** Author: Pat Scott, pat@physto.se
*** Date: 2008-07-12
***********************************************************************

      real*8 function alt_f(u)
      implicit none
      include 'dscapstar.h'
      real*8 u, alt_foverusq
      external alt_foverusq

      alt_f = u*u*alt_foverusq(u)  !Return to f(u) rather than f(u)/u^2
     & * capmx/rhowimp             !Remove density scaling
     & * 1.0d15                    !convert to u*u*(km/s)^(-3)
     
      return
      end


***********************************************************************
*** malc_foverusq is like gal_foverusq, but returns a Fairbairn Fancypants
*** distribution instead.
***
*** Input: u = velocity (km/s) relative to the galactic rest frame
*** Output: f(u) / u^2 [cm^(-3) * (cm/s)^(-3)]
***
*** Author: Pat Scott, pat@physto.se, based on 080623 email from Malcolm Fairbairn
*** Date: 2008-07-07
***********************************************************************

      real*8 function malc_foverusq(u)
      implicit none
      include 'dscapstar.h'
      real*8 u, norm, sigma, alpha, gammconst
      parameter (alpha=0.35d0, gammconst = 37.234d0) !gammconst = gamma(1.d0+1.5d0/alpha)      
      
      sigma = 56.8d0 * sqrt(0.01d0/galr)
      norm = 2.d0**1.5d0 / 3.d0 * sigma**3
     &  * (3.d0-2.d0*alpha)**(-1.5d0/alpha) * gammconst
      malc_foverusq = exp(-(3.d0-2.d0*alpha)
     &  * (0.5d0*(u/sigma)**2)**alpha) / norm    
     &  * 1.0d-15            ! (km/s)^(-3) -> (cm/s)^(-3)
     &  * rhowimp/capmx      ! local number density of WIMPs

      return
      end
