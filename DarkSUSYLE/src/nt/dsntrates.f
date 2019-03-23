      subroutine dsntrates(emuth0,thmax0,rtype,rateea,
     &  ratesu,istat)
c_______________________________________________________________________
c
c            n e u t r a l i n o   b r a n c h i n g   r a t i o s
c            a n d  c a p t u r e  r a t e  i n  t h e   s u n
c            m u o n   f l u x  c a l c u l a t e d
c    november, 1995
c    uses routines by p. gondolo and j. edsjo
c    modified by l. bergstrom and j. edsjo
c    capture rate routines are written by l. bergstrom
c    input:  emuth0  - muon energy threshold in gev
c            thmax0  - muon angel cut in degrees
c            rtype   - 2 = contained events km^-3 yr^-1
c                      3 = through-going events km^-2 yr^-1
c    hidden input: ntcalcmet - 1 use jkg approximations
c                              2 use jkg for sun, full gould for earth
c                              3 use jkg for sun, full gould+dk for earth
c                              4 use full numerical calculations for Sun, Earth
c    output: rateea  - events from earth ann. per km^2(3) per yr
c            ratesu  - events from sun ann. per km^2(3) per yr
c    slightly modified by j. edsjo.
c    modified by j. edsjo 97-05-15 to match new inv. rate convention
c    modified by j. edsjo 97-12-03 to match muflux3.21 routines.
c    modified by p. gondolo 98-03-04 to detach dsntannrate from susy
c    routines.
c    modified by j. edsjo 98-09-07 to fix istat bug.
c    modified by j. edsjo 98-09-23 to use damour-krauss distributions
c      and full earth formulas.
c    modified by j. edsjo 99-03-17 to include better damour-krauss
c      velocity distributions and numerical capture rate integrations
c      for these non-gaussian distributions
c
c=======================================================================
      implicit none
      include 'dssusy.h'
      include 'dshmcom.h'
      include 'dsntcom.h'
      real*8 emuth0,thmax0,rateea,ratesu,arateea,aratesu,
     &  sig_v,mx,yield,sigsip,sigsdp,tmp1,tmp2
      real*8 dssigmav,dsntmuonyield
      integer istat,itmp,rtype

c ----------------------------------------- zero common block data

      tausu=0.0d0
      csu=0.0d0
      tauea=0.0d0
      cea=0.0d0
      ceadk=0.0d0
      gtot10=0.0d0
      ntarateea=0.0d0
      ntaratesu=0.0d0

      rateea=0.0d0
      ratesu=0.0d0

c --------------------------------------------- start calculations

      mx = mass(kn(1))
      if (rhox.eq.0.0d0) then
        rateea=0.0d0
        ratesu=0.0d0
        return
      endif
      call dsddneunuc(sigsip,tmp1,sigsdp,tmp2)

c ------------------------------------ annihilation rates at rest:

      sig_v=dssigmav(0)

c **************************************************************

      if (mx.gt.emuth0) then

        if (ntcalcmet.eq.1.or.ntcalcmet.eq.2
     &    .or.ntcalcmet.eq.4) then ! jkg and/or gould w/ Gauss or full dist.
          call dsntannrate(mx,sigsip,sigsdp,sig_v,
     &       arateea,aratesu)

        elseif (ntcalcmet.eq.3) then ! jkg for sun, gould+dk for earth
          call dsntdkannrate(mx,sigsip,sigsdp,sig_v,arateea,
     &      aratesu)
        else
          write(*,*) 'error in dsntrates: invalid option,',
     &      ' ntcalcmet = ',ntcalcmet
          return
        endif
c...arateea and aratesu in units of 10^24 yr^-1

        yield=dsntmuonyield(emuth0,thmax0,'su',rtype,istat)
c...yield is in units of 10^-30 m^-2(3).
        ratesu=yield*aratesu
c...we now have units km^-2 yr^-1 if rtype=3.
c...one more m^-1 -> km^-1 still to go for rtype=2
        if (rtype.eq.2) then
          ratesu=ratesu*1.0d3
        endif
        itmp=istat
        yield=dsntmuonyield(emuth0,thmax0,'ea',rtype,istat)
c...yield is in units of 10^-30 m^-2(3).
        rateea=yield*arateea
c...we now have units km^-2 yr^-1 if rtype=3.
c...one more m^-1 -> km^-1 still to go for rtype=2
        if (rtype.eq.2) then
          rateea=rateea*1.0d3
        endif
        istat=or(istat,itmp)
      else
        rateea=0.d0
        ratesu=0.d0
      endif

      ntarateea=arateea
      ntaratesu=aratesu

      return

      end




