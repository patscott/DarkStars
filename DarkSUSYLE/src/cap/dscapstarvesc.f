

**********************************************************************
*** dscapstarvesc gives the escape velocity in km/s as a function of
*** the radius r (in meters) from the star's core.
*** author: joakim edsjo (edsjo@physto.se)
*** input: radius in m
*** output escape velocity in km/s
*** date: 2003-11-26
*** Modified: 2007-06-26
*** Modified: 2007-07-13 Pat Scott (pat@physto.se)
**********************************************************************

      real*8 function dscapstarvesc(r)
      implicit none
      include 'dscapstar.h'

      real*8 r,dscapstarpot

      dscapstarvesc=vescsurf*sqrt(abs(dscapstarpot(r)/phisurf))
      
      return
      end
