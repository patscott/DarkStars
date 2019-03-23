**********************************************************************
*** dscapstarpotint gives the gravitational potential inside and outside
*** of a star as a function of the radius r (in meters).
*** in this routine, the actual integration is performed. for speed,
*** use dscapstarpot instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@physto.se)
*** date: april 1, 1999
*** Modified: 2007-06-26
*** Modified: 2007-07-04 Pat Scott (pat@physto.se)
**********************************************************************

      real*8 function dscapstarpotint(r)
      implicit none
      include 'dscapstar.h'

      real*8 r,dscapspfunc,dsLEfint
      external dscapspfunc

c...integrate star density

      if (r.lt.r_star) then
        dscapstarpotint=
     &    -gn * m_star/r_star  ! surface potential
     &    -dsLEfint(dscapspfunc,max(r,100.0d0),r_star,1.d-6)
      else
        dscapstarpotint=-m_star*gn/(max(r,100.0d0))
      endif

      return
      end 
