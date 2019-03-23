
***********************************************************************
*** Auxiliary function for r-integration
***
*** Input: radius in centimeters
***	  Depending upon switches, EITHER:
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
***       OR:
***        foveru [ cm^-3 (cm s^-1)^-3 ] external function f(u)/u^2 with
***          velocity u [ km s^-1 ] as argument.
*** Output: integrand in cm^(-4)*(km/s)*(cm/s)^(-2) cm^-1 s^-1
***
*** l.b. and j.e. 1999-04-06
*** Adapted for general star by J. Edsjo, 2007-06-26
*** Adapted to allow analytical integration over velocity distribution
***  by Pat Scott 2008-03-11
***********************************************************************

      real*8 function dscapcsint(r,foveru)
      implicit none
	  
	  include 'dscapstar.h'
	  
      real*8 mx,sigsi,sigsd,mu,muplus,muminus,ma
      real*8 v,e0,res,r,rx,umin_r,umax_r,ni
      real*8 dscapstarvesc,dscapstardenscomp,dscapcsint1,foveru
      real*8 MBhydrogen, MBheavy
      external dscapstarvesc,dscapstardenscomp,dscapcsint1,foveru
      external MBhydrogen, MBheavy
      integer iel,vtype
      common/capint/mx,ma,sigsi,sigsd,rx,vtype
      common/capint2/iel
      common/capint3/mu,muplus,muminus,v,e0
      common/capint4/umin_r,umax_r

c...determine integration limits
      rx=r
      v=dscapstarvesc(r/100.0d0)   ! escape velocity at r

c...general velocity limits
      umin_r=0.d0                  ! lower velocity limit
      umax_r=sqrt(mu/muminus**2)*v ! upper velocity limit
c...Note, we should only allow such velocities that we scatter to
c...velocities lower than the escape velocity. This means that we have to
c...make sure that (mu/muplus^2 > u^2 / (u^2 + v^2 )). This is the same as
c...the condition umax above.
        	 			    
      if (chopatvesc .or. (altveldist .and. capmode)) then
        call dshiprecint1(dscapcsint1,foveru,-1.d0,1.d0,res)
      else
        res = dscapcsint1(1.d0,foveru)
      endif

      ni = dscapstardenscomp(r/100.0d0,iel) ! target number density in star
      res = res * ni
      dscapcsint = res * r * r  ! assume spherical star

      return
      end
