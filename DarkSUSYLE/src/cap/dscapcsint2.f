
***********************************************************************
*** Auxiliary function for velocity integration
***
*** Input: velocity relative to star [km/s]
***		  Depending upon switches, EITHER:
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
***		  OR:
***        foveru [ cm^-3 (cm s^-1)^-3 ] external function f(u)/u^2 with
***          velocity u [ km s^-1 ] as argument.
*** Output: integrand [cm^(-3)*(cm/s)^-2]
***
*** We here follow the analysis in Gould, ApJ 321 (1987) 571 and more
*** specifically the more general expressions in appendix A.
***
*** Adapted for general star by J. Edsjo, 2007-06-26
*** Minor alteration to last line by Pat Scott, 2007-07-11
*** Adapted to be compatible with higher level routines which
***  allow analytical integration over velocity distribution 
***  by Pat Scott 2008-03-11
*** Adapted some more to work with finite galactic escape velocities
***  and generic galactic-frame velocity distributions by
***  Pat Scott 2008-03-18
***********************************************************************

      real*8 function dscapcsint2(u,foveru)
      implicit none
      include 'dscapstar.h'
      include 'dsparam.h'
      real*8 u,mu,mx,sigsi,sigsd,muplus,muminus,ma
      real*8 r,e0,sla1,sla2,v,foveru,costheta
      integer iel,vtype
      external foveru
      common/capint/mx,ma,sigsi,sigsd,r,vtype
      common/capint2/iel
      common/capint3/mu,muplus,muminus,v,e0
      common/capint5/costheta

c...Note, in principle we should have a step function here to allow
c...only such scatterings that go to a velocity lower than the escape
c...velocity, Heaviside(mu/muplus^2 - u^2/(u^2+v^2)),  however, this
c...is taken care of in the limits for the u integration and is thus not
c...needed here.
c      if (u**2/(u**2+v**2).gt.mu/muplus**2) then
c        dscapcsint2=0.0d0
c        return
c      endif

      if (iel.eq.1) then  ! no form-factor suppression for Hydrogen        

        sla1=((v/c_0)**2-muminus**2/mu*(u/c_0)**2) ! gould (2.13)
        if (chopatvesc .or. altveldist) then
        !Note that foveru is actually f over u squared here
          dscapcsint2 = u * 1.d5 * foveru(sqrt(u**2 + v_star**2 + 
     &     2.d0*u*v_star*costheta))
        else
          dscapcsint2=foveru(u)
        endif
        dscapcsint2 = dscapcsint2 * sla1 ! gould (2.13) + (2.8)
		
      else

c		Below we assume an exponential form factor, which lets us
c		integrate over the energy loss analytically. For general form factors
c		we would need to generalize the expressions with an integral
c		over the energy loss (or q^2)

        sla1=exp(-mx*(u/c_0)**2/2./e0)  ! gould (A6)
        sla2=exp(-mu/muplus**2*mx*((v/c_0)**2+(u/c_0)**2)/2./e0) ! gould (A6)
		if (chopatvesc .or. altveldist) then
		!Note that foveru is actually f over u squared here
		  dscapcsint2 = u * 1.d5 * foveru(sqrt(u**2 + v_star**2 + 
     &     2.d0*u*v_star*costheta))
        else
          dscapcsint2=foveru(u) 
        endif
		dscapcsint2 = dscapcsint2 * (sla1 - sla2)  ! gould (A6) + (2.8)

      endif

      return
      end
