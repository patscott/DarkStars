
***********************************************************************
*** Auxiliary function for theta integration
***
*** Input: costheta = cos(angle between incoming WIMP and stellar,
***                       where theta=0 is a head-on collision)
***       Depending upon switches, EITHER:
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
***       OR:
***        foveru [ cm^-3 (cm s^-1)^-3 ] external function f(u)/u^2 with
***          velocity u [ km s^-1 ] as argument.
*** Output: integrand [cm^(-3)*(km/s)*(cm/s)^(-2)]
***
*** Pat Scott, pat@physto.se, 2008-03-18
***********************************************************************

      real*8 function dscapcsint1(costheta,foveru)
      implicit none
	  
      include 'dscapstar.h'
	  
      real*8 mx,sigsi,sigsd,mu,muplus,muminus,ma,costheta_int
      real*8 v,e0,res,costheta,rx,umin,umax,umin_r,umax_r
      real*8 dscapcsint2,foveru
      real*8 MBhydrogen, MBheavy
      external dscapcsint2,foveru
      external MBhydrogen, MBheavy
      integer vtype,iel
      common/capint/mx,ma,sigsi,sigsd,rx,vtype
      common/capint2/iel
	  common/capint3/mu,muplus,muminus,v,e0
      common/capint4/umin_r,umax_r
      common/capint5/costheta_int
	  
      costheta_int = costheta
	  
c     theta-specific velocity limits
      umin=umin_r				 ! lower velocity limit
      if (chopatvesc) then       ! upper velocity limit
        umax=min(umax_r,sqrt(galesc**2+v_star**2+2.d0*
     &   galesc*v_star*costheta))
      else
        umax=umax_r
      endif

      if (umin.ge.umax) then
        dscapcsint1=0.0d0
      return
      endif
        	  
      if (capmode) then
c       do numerical velocity integration
        call dshiprecint2(dscapcsint2,foveru,umin,umax,res)
      else
c       do analytical velocity integration
c       Note that the value of altveldist plays no role here - 
c       if capmode is false, the analytical solution for a 
c       Maxwell-Boltzmann distribution (either full or truncated,
c       depending upon chopatvesc) is returned regardless of any attempt
c       to choose an alternative velocity distribution.
        if (iel.eq.1) then
          res = MBhydrogen(umax) - MBhydrogen(umin)
        else
          res = MBheavy(umax) - MBheavy(umin)
        endif
      endif

      dscapcsint1 = res

      if (dscapcsint1.lt.0.d0) dscapcsint1=0.d0

      return
      end
	  
