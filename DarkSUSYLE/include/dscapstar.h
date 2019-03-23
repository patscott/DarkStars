*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dscapstar.h                           ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  authors: 	joakim edsjo	(edsjo@physto.se)	2007-06-26
c		pat scott	(pat@physto.se)		2007-07-11, 2008-03-xx

c...general switches etc
      logical capmode,interpmode,potmode,chopatvesc,altveldist
      integer integrator
      common /dscap/capmode,interpmode,potmode,chopatvesc,altveldist,
     &  integrator
      save /dscap/

c...parameters
c...r_star = star's radius (m)
c...m_star = star's mass (kg)
c...vd_3d_star = 3D velocity dispersion at star (km/s)
c...v_star = velocity of star with respect to halo (km/s)
c...galesc = local galactic escape velocity (km/s)
c...galr = local galactocentric distance (pc)
c...vescsurf = escape velocity at surface (km/s, calculated)
c...phisurf = potential at surface (m^2 s^-2, calculated)
      real*8 r_star,m_star,vd_3d_star,v_star,
     &  vescsurf,phisurf,capmx,rhowimp,galesc,galr
      common /dscappar/r_star,m_star,vd_3d_star,v_star,
     &  vescsurf,phisurf,capmx,rhowimp,galesc,galr
      save /dscappar/

      integer meshpoints,n_species
      parameter(meshpoints=200,n_species=22)
      real*8 starr(meshpoints),starm(meshpoints),
     &  starphi(meshpoints),starmfr(n_species,meshpoints),
     &  starmd2(meshpoints),starphid2(meshpoints),
     &  starmfrd2(n_species,meshpoints),starrho(meshpoints) 
      real*8 staraa(n_species),starza(n_species),
     &  starma(n_species),starma_atomic(n_species)
      common /dsstar/starr,starm,starphi,starmfr,
     &  starmd2,starphid2,starmfrd2,
     &  starrho,staraa,starza,starma,starma_atomic
      save /dsstar/

      real*8 pi,gn
      parameter(pi=3.141592653589793238d0)
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-2
                  
************************* end of dscapstar.h **************************
***********************************************************************





