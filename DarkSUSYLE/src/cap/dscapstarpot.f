***********************************************************************
*** This routine uses a derived potential for the star
*** The data in starphi() is calculated by dscapstarpotcalc.f. 
***********************************************************************

***********************************************************************
*** dscapstarpot gives the potential in a star as a function of radius
*** the radius should be given in m and the potential is returned in
*** m^2 s^-2
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2003-11-26
*** Modified: 2007-06-26
*** Modified: 2007-07-13 Pat Scott (pat@physto.se)
***********************************************************************

      real*8 function dscapstarpot(r)
      implicit none

      include 'dscapstar.h'
      real*8 r,rpl,dscapstarpotint
      integer i,j


      if (potmode) then

        dscapstarpot = dscapstarpotint(r)
        
      else
        
        !To tabulate the potential (this is up to the user),
        !one must call dscapstarpotcalc

        if (r.ge.r_star) then
          dscapstarpot=-m_star*gn/r
          return
        endif

        if (r.le.0.0d0) then
          dscapstarpot=starphi(1)
          return
        endif

        !Interpolate in table
        if (interpmode) then

          call dssplint(starr,starphi,starphid2,meshpoints,r/r_star,
     &     dscapstarpot)
                
        else

          call dshunt(starr,meshpoints,r/r_star,j)
          if (j.lt.meshpoints) goto 20

          dscapstarpot=starphi(meshpoints)
          return

 20       rpl=(r-starr(j)*r_star)/(starr(j+1)*r_star-starr(j)*r_star)

          dscapstarpot=starphi(j)*(1.0d0-rpl)+starphi(j+1)*rpl

        endif

      endif

      return

      end
