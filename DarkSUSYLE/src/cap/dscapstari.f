***********************************************************************
*** dscapstari gives the capture rate of neutralinos in a star
*** given a specified velocity distribution. the integrations over
*** the star's radius and over the incoming angle are
*** performed numerically, and the integral over the velocity
*** can be performed either numerically or analytically.
***
*** Input: mx [ gev ]
***        sigsi [ cm^2 ]
***        sgisd [ cm^2 ]
***		  Depending upon switches, EITHER:
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
***       OR:
***        foveru [ cm^-3 (cm s^-1)^-3 ] external function f(u)/u^2 with
***          velocity u [ km s^-1 ] as argument.
*** Output: WIMPs captured per second (s^-1)
***
*** Created:  2003-11-26 Joakim Edsjo (edsjo@physto.se)
*** Modified: 2007-06-26 Joakim Edsjo
*** Modified: 2007-07-05 Pat Scott (pat@physto.se)
*** Modified: 2008-03-11 Pat Scott
***********************************************************************

      real*8 function dscapstari(mx,sigsi,sigsd,foveru)
      implicit none

      include 'dscapstar.h'
	  include 'dsparam.h'

      real*8 mx,sigsi,sigsd,ma
      real*8 mxx,max,sigsix,sigsdx,rmin,rmax,res,rx
      real*8 mu,muplus,muminus,v,e0,sla1,siga
	  real*8 dscapcsint,foveru
      external dscapcsint,foveru
      integer i,ix,vtype,vt
      common/capint/mxx,max,sigsix,sigsdx,rx,vtype
      common/capint2/ix
      common/capint3/mu,muplus,muminus,v,e0

c...perform integration as given in Gould, ApJ 521 (1987) 571.
      mxx=mx
      sigsix=sigsi
      sigsdx=sigsd
      dscapstari=0.0d0

c...sum over elements
      do 10 i=1,n_species

        ma=starma(i) ! nucleus mass
        max=ma
        mu=mxx/max
        muplus=(mu+1.d0)/2.d0
        muminus=(mu-1.d0)/2.d0
        ix=i
	e0=1.5d0/ma*0.038938d0 ! gives units of gev; see gould (a8)
        e0=e0/(0.91d0*ma**(0.33333d0)+0.3d0)**2

c       begin radial integration
        rmin=0.d0
        rmax=r_star*100.0d0  ! star radius in centimeters
        call dshiprecint(dscapcsint,foveru,rmin,rmax,res)

        siga=sigsix*staraa(ix)**2*(mxx*max/(mxx+max))**2/
     &   (mxx*m_p/(mxx+m_p))**2
        if (ix.eq.1) then
          siga = siga + sigsdx       ! add spin-dependent for Hydrogen
        else
          sla1=muplus**2*2.*e0/mu/mx ! gould (A6)
          res = res * sla1           ! form factor constants for other nuclei
        endif

        res = res * siga
        res = res * (c_0*1.0d5)**2   ! reinstate the two long lost factors of c
        res = res * 1.0d5            ! km/s -> cm/s
        res = res * 4.d0*pi          ! put in the 4pi from the spherical integral
        if (res.lt.0.d0) res = 0.d0

        dscapstari=dscapstari+res
  
 10   continue

      return
      end
