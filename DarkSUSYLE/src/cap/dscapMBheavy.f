
***********************************************************************
*** Analytical solution to velocity integral for heavy elements given a
*** Maxwell-Boltzmann distribution of velocities (either full or
*** truncated), as seen in the frame of an observer moving with
*** speed v_star relative to the galactic rest frame.  No constant of
*** integration is included, and exponential form factors are assumed.
***
*** Input:  u         velocity relative to star [km/s]
*** Output: integral                            [cm^(-3)*(km/s)*(cm/s)^-2]
***
*** Pat Scott 2008-03-11, pat@physto.se
***********************************************************************
      
      real*8 function MBheavy(u)
      implicit none
      include 'dscapstar.h'
      include 'dsparam.h'
      real*8 u,mu,mx,sigsi,sigsd,muplus,muminus,ma,costheta
      real*8 r,e0,v,erf,A,B,E,F,G,H,J
      real*8 part1,part2,part3,part4,normfactor
      integer iel,vtype
      external erf
      common/capint/mx,ma,sigsi,sigsd,r,vtype
      common/capint3/mu,muplus,muminus,v,e0
      common/capint5/costheta

      A = 1.5d0 / vd_3d_star**2.d0
      B = v_star
      G = mx / (2.d0 * e0 * c_0**2) 
      H = mu * mx / (2.d0 * muplus**2 * e0 *c_0**2) 
      J = v
	  
      if (chopatvesc) then
	  
        E = A*A*B*B*costheta*costheta
        F = -2.0d-10*(rhowimp/mx)*A**1.5d0/sqrt(pi)
        part1 = exp(-A*B*B - (A+G)*u*u - 2.d0*A*B*u*costheta) /
     &   (2.d0*(A+G))
        part2 = A * B * sqrt(pi) * costheta * exp(-A*B*B + E/(A+G))
        part2 = part2 * erf( ((A+G)*u + A*B*costheta) / sqrt(A+G) ) 
     &   / (2.d0 * (A+G)**1.5d0)
        part3 = A * B * sqrt(pi) * costheta *
     &   erf( ((A+H)*u + A*B*costheta) / sqrt(A+H) ) *
     &   exp( E / (A+H) - A*B*B - H*J*J ) / (2.d0 * (A+H)**1.5d0)
        part4 = 0.5d0 * exp(-2.d0*A*B*u*costheta -
     &   A * (B*B + u*u) - H * (J*J + u*u)) / (A+H)
        MBheavy = F * (part1 + part2 - part3 - part4)
	
      else

        F = 1.0d-10*(rhowimp/mx)*sqrt(3.d0/2.d0)/(vd_3d_star*v_star)
     &   /sqrt(pi)
        part1 = erf((G*u+A*(u-B))/sqrt(A+G)) - erf((G*u+A*(u+B))
     &   /sqrt(A+G))
        part1 = exp(-A*B**2*G/(A+G)) * sqrt(pi) / 2.d0 / sqrt(A+G)
     &   * part1
        part2 = erf((H*u+A*(u-B))/sqrt(A+H)) - erf((H*u+A*(u+B))
     &   /sqrt(A+H))
        part2 = exp(-H*(H*J**2+A*(B**2+J**2))/(A+H)) * sqrt(pi) / 2.d0 
     &   / sqrt(A+H) * part2
        MBheavy = F * (part1 - part2)

      endif
  
      return
      end
