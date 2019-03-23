      subroutine dsntset(c)
c...set parameters for neutrino telescope routines
c...  c - character string specifying choice to be made
c...author: joakim edsjo, 2000-08-16
      implicit none
      include 'dsntcom.h'
      include 'dsntdkcom.h'
      integer i
      character*(*) c

c...Make sure block datas are linked to. This is not needed anymore
c      call dsntdkcom
c      call dsntcapcom
      
c...Initialize tables
c      call dsntsunread   ! this is now done in the routines that need it

c...This is set by a block data
c...Default value for DK population
c      dklambda=1     ! amount of non-conservation of angular momentum
c                     ! 1 = angular momentum conserved to 100%.

c...approximate formulae given in jungman, kamionkowski, griest.
      if (c.eq.'jkg') then
         ntcalcmet=1           ! calculate method

c...approximate as above for the sun, but the full formulae from
c...gould, apj 321 (1987) 571 for the earth
      else if (c.eq.'gould') then
         ntcalcmet=2

c...approximate as above for the sun, but the full formulae from
c...gould, apj 321 (1987) 571 for the earth plus including the
c...damour-krauss distribution of wimps captured by the earth.
      else if (c.eq.'dknum') then
         ntcalcmet=3
         nttab=0

c...The same as above, but use tabulated results of the integration
c...instead of numerically evaluating them each time
      else if (c.eq.'dktab') then
         ntcalcmet=3
         nttab=1

c...Use full expressions in Gould, ApJ 321 (1987) 571 for the Earth, but
c...instead of using analytic expressions for a Gaussian and approximate
c...positions within the Earth for the elements, the velocity profile
c...specified in dshmset.f (i.e. by the parameter veldfearth) is
c...integrated numerically and a full integration over the Earth's
c...interior is performed. For the Sun, the same full expressions are used
c...and the numerical profile as specified in dshmset.f (i.e. by the 
c...parameter veldf) is integrated numerically over a realistic Sun
c...profile.
      else if (c.eq.'num') then
         ntcalcmet=4
         nttab=0

c...Same as above, but use tabuled versions of the numerical integrations.
c...If the tables do not exist, they are recreated with the numerical
c...routines as in option 'num' above.
      elseif (c.eq.'tab'.or.c.eq.'default') then
        ntcalcmet=4
        nttab=1

c...invalid choice
      else
         write (*,*) 'dsntset: unrecognized option ',c
         stop
      endif

      return
      end




