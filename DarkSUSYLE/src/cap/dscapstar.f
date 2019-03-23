***********************************************************************
*** dscapstar calculates the capture of WIMPs on a general star
***
*** Input: mx = WIMP mass (GeV)
***        sigsi = spin-independent scattering cross section (cm^2)
***        sigsd = spin-dependent scattering cross section (cm^2)
*** Output: total capture rate (s^-1)
***
*** Input is also given as the stellar compositional density distribution
*** (starmfr), interior mass (starm) and radius (starr), as defined on
*** a mesh of points.  The array starmfr is also defined for different
*** elements.  The mesh is typically chosen as the same as the mesh over
*** which a stellar model has been calculated.  These arrays reside in
*** the common block dsstar (dscapstar.h).
***
*** There are also options in the common block dscap (dscapstar.h):
***    capmode: T=numerical integration over velocity distribution
***               and radius
***             F=analytical integration over velocity distribution
***               (Gaussian distribution assumed)
***    potmode: T=stellar interior potentials calculated explcitly every
***               time they are required
***             F=stellar interior potential tabluated at each meshpoint and
***               interpolated between
*** interpmode: T=cubic spline interpolation in tables
***				F=linear interpolation in tables
*** chopatvesc: T=truncate halo velocity distribution at the local galactic
***               escape velocity.  Works for Gaussian or user-supplied
***               distribution functions.
***             F=assume halo velocity distribution extends to infinity
*** altveldist: T=use something other than a Gaussian distribution (must be
***               supplied by the user).  Not compatible with capmode = F.
***             F=use standard approximation of an isothermal halo velocity
***               distribution (a full or truncated Gaussian, depending upon
***               chopatvesc).
***
*** Created: June 26, 2007 Joakim Edsjo (edsjo@physto.se)
*** Modified: 2007-07-14 Pat Scott (pat@physto.se)
*** Modified: 2008-02-29 Pat Scott
*** Modified: 2008-03-11 Pat Scott
***********************************************************************

      real*8 function dscapstar(mx,sigsi,sigsd)
      implicit none
      include 'dscapstar.h'

      real*8 mx,sigsi,sigsd,dscapstari,normfactor,erf,res
      real*8 dscapfoveru,gal_foverusq,alt_foverusq,dsLEfint
      real*8 starmfr_temp(meshpoints),starmfrd2_temp(meshpoints)
      integer i,j
      external dscapfoveru,gal_foverusq,alt_foverusq,erf,alt_f
      
      capmx=mx ! so that internal routines know about the mass
      phisurf=-gn * m_star/r_star  ! surface potential (m^2 s^-2)
      vescsurf=sqrt(-2.0d0*phisurf)*1.d-3 ! surface escape vel (km/s)
      if (chopatvesc) then
        if (altveldist) then
          normfactor = dsLEfint(alt_f,0.d0,galesc,1.d-6)
        else
          normfactor=erf(sqrt(3.d0/2.d0)*galesc/vd_3d_star)
     &     -sqrt(6.d0/pi)*galesc/vd_3d_star
     &     *exp(-1.5d0*galesc**2/vd_3d_star**2)
        endif
      else
          normfactor = 1.d0
      endif
      

      if (interpmode) then
        
        call dsspline(starr,starm,meshpoints,0.0d0,0.0d0,
     &   starmd2)
        do i=1,n_species
          do j=1,meshpoints
            starmfr_temp(j) = starmfr(i,j)
          enddo
        
          call dsspline(starr,starmfr_temp,meshpoints,0.0d0,
     &     0.0d0,starmfrd2_temp)
        
          do j=1,meshpoints
            starmfrd2(i,j) = starmfrd2_temp(j)
          enddo
        enddo
      
      else
      
      endif

!...If using tabulated potentials, calculate and tabulate the potential
!...inside the star
      if (.not.potmode) then 
        call dscapstarpotcalc
        if (interpmode) call dsspline(starr,starphi,
     &     meshpoints,0.0d0,0.0d0,starphid2)
      endif
	  
      if (altveldist) then
	    
        if(.not.capmode) write(*,*) 'ERROR: bad switches;
     &   ignoring altveldist and reverting to analytical integration over
     &   Maxwell-Boltzmann velocity distribution.'
        res=dscapstari(mx,sigsi,sigsd,alt_foverusq)

      else
	  
        if (chopatvesc) then
          res=dscapstari(mx,sigsi,sigsd,gal_foverusq)
        else
          res=dscapstari(mx,sigsi,sigsd,dscapfoveru)
        endif
	  
      endif
	  
      dscapstar = res / normfactor

      return
      end

