
***********************************************************************
*** Analytical solution to velocity integral for Hydrogen given a
*** Maxwell-Boltzmann distribution of velocities (either full or
*** truncated), as seen in the frame of an observer moving with
*** speed v_star relative to the galactic rest frame. No constant of
*** integration is included.
***
*** Input:  u         velocity relative to star [km/s]
*** Output: integral                            [cm^(-3)*(km/s)*(cm/s)^-2]
***
*** Pat Scott 2008-03-11, pat@physto.se
***********************************************************************

      real*8 function MBhydrogen(u)
      implicit none
      include 'dscapstar.h'
      include 'dsparam.h'
      real*8 u,mu,mx,sigsi,sigsd,muplus,muminus,ma,costheta
      real*8 r,e0,v,A,B,C,D,E,F,part1,part2,part3,erf
      real*8 normfactor
      integer vtype
      external erf
      common/capint/mx,ma,sigsi,sigsd,r,vtype
      common/capint3/mu,muplus,muminus,v,e0
      common/capint5/costheta

      A = 1.5d0 / vd_3d_star**2.d0
      B = v_star
      C = (v/c_0)**2
      D = muminus**2/mu/c_0**2
	  
      if (chopatvesc) then
	  
        E = -2.0d0*A*C + 2.d0*A*B*B*D*costheta*costheta
        F = 2.0d-10*(rhowimp/mx)*A**1.5d0/sqrt(pi)
        part1 = E + D + D + 2.d0*A*D*(u**2 - B*u*costheta)
        part1 = part1 * exp(-2.d0 * A * B * u * costheta - 
     &   A * (B * B + u * u))
        part2 = (E + 3.d0*D) * erf(sqrt(A) * (u + B * costheta))
        part2 = part2 * sqrt(A*pi) * B * costheta * 
     &   exp(A*B*B*(costheta*costheta-1.d0))
        MBhydrogen = F * 0.25d0/A/A * (part1+part2)
			    
      else
	  
        E = sqrt(pi) * (erf(sqrt(A)*(u-B)) - erf(sqrt(A)*(u+B)))
        F = 1.0d-10*(rhowimp/mx)*sqrt(3.d0/2.d0)/
     &   (vd_3d_star*v_star)/sqrt(pi)
        part1 = C * E / sqrt(A) / 2.d0
        part2 = (1.d0 + 2.d0 * A * B**2) * E
        part3 = 2.d0 * sqrt(A) * (exp(-A * (u+B)**2) * (B - u) +
     &   (B + u) * exp(-A * (u + B)**2 + 4.d0 * A * B * u))
        MBhydrogen = F * (part1 + D / 4.d0 / A**(1.5d0)
     &   * (part3 - part2))

      endif

      return
      end
