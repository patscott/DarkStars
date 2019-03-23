***********************************************************************
*** dsntsundenscomp gives the number density of nucleons of atomic
*** number Z per cm^3
*** input: radius - in meters
***        itype: internal element number (as defined in dsntsunread.f)
***           1 = H
***           2 = He4
***           3 = He3
***           4 = C12
***           5 = N14
***           6 = O8
***           7 = Ne
***           8 = Na
***           9 = Mg
***          10 = Al
***          11 = Si
***          12 = S
***          13 = Ar
***          14 = Ca
***          15 = Fe
***          16 = Ni
*** Author: joakim edsjo
*** Date: 2003-11-26
*** Modified: 2006-03-21 (atomic mass unit fix (was off by 6%)) JE
***********************************************************************

      real*8 function dsntsundenscomp(r,itype)
      implicit none
      include 'dsparam.h'
      real*8 r,dsntsundens,fraction,dsntsunmfrac,m_a
      integer itype

      include 'dssun.h'

      if (itype.lt.1.or.itype.gt.16) then
        write(*,*) 'WARNING in dsntsundenscomp: illegal element type: ',
     &    itype
      endif

c...m_a is the mass in atomic mass units, amu (where C12 has mass 12 amu)
c...sdma(4) is the mass of C12
      m_a = sdma(itype)/sdma(4)*12.0d0  ! JE fix 060321
      fraction=dsntsunmfrac(r,itype)

c...m_a is now in atomic units u which is what we want

      dsntsundenscomp=fraction*dsntsundens(r)/m_a*N_a ! number/cm^3

      return
      end
