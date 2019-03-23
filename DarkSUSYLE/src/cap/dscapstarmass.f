***********************************************************************
*** dscapstarmass gives the mass of the star as a function of radius
*** the radius should be given in m and the mass is given in kg
*** up to the specified radius.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2003-11-25
*** Modified: 2007-06-26
*** Modified: 2007-07-05 Pat Scott (pat@physto.se)
***********************************************************************

      real*8 function dscapstarmass(r)
      implicit none

      include 'dscapstar.h'
      real*8 r,rpl,dscapstarmassfirst
      integer i,j
      
      if (r.ge.r_star) then
        dscapstarmass=m_star
        return
      endif

      if (r.le.0.0d0) then
        dscapstarmass=0.0d0
        return
      endif

      if (interpmode) then

        call dssplint(starr,starm,starmd2,meshpoints,r/r_star,
     &   dscapstarmass)

        if (dscapstarmass.gt.1.0d0) dscapstarmass = 1.0d0
        if (dscapstarmass.lt.0.0d0) dscapstarmass = 0.0d0
        dscapstarmass = dscapstarmass * m_star        

      else

        if (r.ge.0.0d0.and.r.lt.starr(2)*r_star) then ! do a better interpolation
          dscapstarmass=starm(2)*m_star*r**3/(starr(2)*r_star)**3
          return
        endif

        call dshunt(starr,meshpoints,r/r_star,j)
        if (j.lt.meshpoints) goto 20

        dscapstarmass=m_star
        return

 20     rpl=(r-starr(j)*r_star)/(starr(j+1)*r_star-starr(j)*r_star)

        dscapstarmass=(starm(j)*(1.0d0-rpl)+starm(j+1)*rpl)*m_star

      endif

      return

      end
