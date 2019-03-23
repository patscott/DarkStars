! Utilities module for DarkStars
!
! Pat Scott, July 2007; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE DkStrs_utils
      use star_data
      use star_extras
      use DkStrs_data
      implicit none

      double precision :: store_temp

      contains


      DOUBLE PRECISION FUNCTION alpha(r)
      ! Takes local height and returns value of WIMP thermal diffusivity (alpha)
      ! Input:    r   radius [cm]
      ! Output:       alpha [dimensionless]

      double precision, intent(IN) :: r

      if (interpmode) then
        call dssplint(starr,staralpha,staralphad2,meshpoints,r/r_star*1.d-2,alpha)
      else
      endif      

      return      

      END FUNCTION alpha


      DOUBLE PRECISION FUNCTION temperature(r)
      ! Takes local height and returns value of local temperature
      ! Input:    r   radius [cm]
      ! Output:       temperature [K]

      double precision, intent(IN) :: r
      double precision :: HVAL
      integer :: exit_code

      if (interpmode) then
        temperature = HVAL(r/r_star*1.d-2,meshpoints,starr,startemperature,
     &   startemperatured1,temptension,exit_code)
         if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in temperature...')
      else
      endif      

      return

      END FUNCTION temperature


      DOUBLE PRECISION FUNCTION tempgrad(r)
      ! Takes local height and returns value of local temperature gradient
      ! Input:    r   radius [cm]
      ! Output:       temperature gradient [K/cm]
      use ez_cycle_data

      double precision, intent(IN) :: r
      double precision :: HPVAL
      integer :: exit_code

      if (interpmode) then
      tempgrad = HPVAL(r/r_star*1.d-2,meshpoints,starr,startemperature,
     &   startemperatured1,temptension,exit_code)/r_star*1.d-2
         if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in tempgrad...')
      else
      endif

      return      

      END FUNCTION tempgrad
	  

      REAL(16) FUNCTION WIMPdens_LTE(r)
      ! Takes local height and returns value of local WIMP number density,
      ! where the WIMPs and nuclei are considered to be in LTE
      ! Input:    r   radius [cm]
      ! Output:       WIMP density [number/cm^3]

      double precision, intent(IN) :: r
      double precision :: HVAL
      integer :: exit_code

      if (interpmode) then
        WIMPdens_LTE = real(HVAL(r/r_star*1.d-2,meshpoints,starr,starLTErhoWIMP,
     &   starLTErhoWIMPd1,LTErhoWIMPtension,exit_code),16)
         if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in WIMPdens_LTE...')      
      else
      endif

      WIMPdens_LTE = WIMPdens_LTE * n_WIMPs_centre_LTE
      if (WIMPdens_LTE.lt.0.d0) WIMPdens_LTE = 0.0d0

      return      

      END FUNCTION WIMPdens_LTE

   
      REAL(16) FUNCTION WIMPdens_isothermal(r)
      ! Takes local height and returns value of local WIMP number density,
      ! where the WIMPs are assumed to exist as an isothermal sphere of
      ! temperature Tw
      ! Input:    r   radius [cm]
      ! Output:       WIMP density [number/cm^3]

      double precision, intent(IN) :: r

      WIMPdens_isothermal = real(n_WIMPs,16)*expEBeta(r)/partition_func
	  	  
      END FUNCTION WIMPdens_isothermal


      DOUBLE PRECISION FUNCTION ambientDens(t)
      ! WIMP ambient density calculator (for stars on explicit orbits)
      ! Input:    t   age of star [years]
      ! Output:       WIMP ambient density [GeV cm^-3] 
      
      double precision, intent(IN) :: t
      double precision :: t_eff, HVAL
      integer :: exit_code
	  
      if(loopOrbit .or. t.le.orbit_len) then
        t_eff = mod(t,orbit_len)
        ambientDens = HVAL(t_eff,fileLen,orbital_t(:),orbital_rho(:),
     &   orbital_rho_d1(:),orbital_rho_tension(:),exit_code)
	if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in ambientDens...')      
      else
	ambientDens = orbital_rho(fileLen)
      endif
	  
      END FUNCTION ambientDens
	  
      
      DOUBLE PRECISION FUNCTION starVelocity(t)
      ! Stellar velocity calculator (for stars on explicit orbits)
      ! Input:    t   age of star [years]
      ! Output:       Stellar velocity relative to DM halo [km s^-1] 
      
      double precision, intent(IN) :: t
      double precision :: t_eff, HVAL
      integer :: exit_code
	   
      if(loopOrbit .or. t.le.orbit_len) then
        t_eff = mod(t,orbit_len)
        starVelocity = HVAL(t_eff,fileLen,orbital_t,orbital_v,orbital_v_d1,orbital_v_tension,exit_code)
        if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in starVelocity...')
      else
        starVelocity = orbital_v(fileLen)
      endif
	  
      END FUNCTION starVelocity


      DOUBLE PRECISION FUNCTION galactic_vesc(t)
      ! Local galactic escape velocity calculator (for stars on explicit orbits)
      ! Input:    t   age of star [years]
      ! Output:       escape velocity [km s^-1] 
      
      double precision, intent(IN) :: t
      double precision :: t_eff, HVAL
      integer :: exit_code
	  
      if(loopOrbit .or. t.le.orbit_len) then
        t_eff = mod(t,orbit_len)
        galactic_vesc = HVAL(t_eff,fileLen,orbital_t(:),orbital_galesc(:),
     &   orbital_galesc_d1(:),orbital_galesc_tension(:),exit_code)
	if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in galactic_vesc...')      
      else
	galactic_vesc = orbital_galesc(fileLen)
      endif
	  
      END FUNCTION galactic_vesc

	  
      DOUBLE PRECISION FUNCTION galactic_r(t)
      ! Local galactocentric distance calculator (for stars on explicit orbits)
      ! Input:    t   age of star [years]
      ! Output:       galactocentric distance [pc] 
      
      double precision, intent(IN) :: t
      double precision :: t_eff, HVAL
      integer :: exit_code
	  
      if(loopOrbit .or. t.le.orbit_len) then
        t_eff = mod(t,orbit_len)
        galactic_r = HVAL(t_eff,fileLen,orbital_t(:),orbital_galr(:),
     &   orbital_galr_d1(:),orbital_galr_tension(:),exit_code)
	if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in galactic_r...')      
      else
	galactic_r = orbital_galr(fileLen)
      endif
	  
      END FUNCTION galactic_r


      DOUBLE PRECISION FUNCTION drhoonvdt(t)
      ! Calculates the the time derivative of rho/v, where rho is the ambient
      ! density of dark matter and v is the stellar velocity (for stars on explicit orbits).
      ! Useful for determining how long a timestep to allow the code to take along an
      ! elliptical orbit.
      !	  
      ! Input:    t   age of star [years]
      ! Output:       d(rho/v)dt [GeV cm^-3]/[km s^-1]/[yr] 

      double precision, intent(IN) :: t
      double precision :: t_eff, v, rho, drhodt, dvdt
      double precision :: HPVAL
      integer :: exit_code1, exit_code2
	  
      rho = ambientDens(t)
      v = starVelocity(t)
	  
      if(loopOrbit .or. t.le.orbit_len) then
        t_eff = mod(t,orbit_len)
        drhodt = HPVAL(t_eff,fileLen,orbital_t(:),orbital_rho(:),
     &   orbital_rho_d1(:),orbital_rho_tension(:),exit_code1)
        dvdt = HPVAL(t_eff,fileLen,orbital_t(:),orbital_v(:),
     &   orbital_v_d1(:),orbital_v_tension(:),exit_code2)
        if (exit_code1.lt.0 .or. exit_code2.lt.0) call die_quietly('Dying quietly - tensional splines snapped in drhoonvdt...')      
      else
        drhoonvdt = 0.d0
        return
      endif 
	  
      drhoonvdt = drhodt/v - rho*dvdt/v**2
	  
      return
	  
      END FUNCTION drhoonvdt


      PURE DOUBLE PRECISION FUNCTION mfp(i)
      ! Computes the value of the WIMP mean free path at a height of
      ! starr(i)*r_star*1.d2 cm
      ! Input:    i   meshpoint index [dimensionless]
      ! Output:       WIMP mean free path [cm]

      integer, intent(IN) :: i

      IF (DoWimpyThings) THEN
        mfp = 1.d0/sum(sigma_i*starmfr(:,i))
      ELSE
        mfp = 0.d0
      ENDIF

      END FUNCTION mfp      

  
      PURE DOUBLE PRECISION FUNCTION kappa(i)
      ! Computes the value of the WIMP thermal conductivity (kappa)
      ! at a height of starr(i)*r_star*1.d2 cm
      ! Input:    i   meshpoint index [dimensionless]
      ! Output:       kappa [dimensionless]

      integer, intent(IN) :: i

      kappa = 1.d0/mfp(i)/sum(sigma_i*starmfr(:,i)/kappa_i)

      END FUNCTION kappa


      DOUBLE PRECISION FUNCTION K_suppression_func(K)
      ! Knudsen-dependent suppression function for conductive WIMP luminosity
      ! Input:    K   Knudsen number [dimensionless]
      ! Output:       suppression factor [dimensionless]
   
      double precision, intent(IN) :: K
      
      K_suppression_func = 1.d0/(1.d0+(K/K0)**(1.d0/Knudsen_suppression_tau))

      return      

      END FUNCTION K_suppression_func


      DOUBLE PRECISION FUNCTION r_suppression_func(r)
      ! Radius-dependent suppression function for conductive WIMP luminosity
      ! Input:    r   radius [cm]
      ! Output:       suppression factor [dimensionless]
   
      double precision, intent(IN) :: r
      
      r_suppression_func = 1.d0 + ((r-r_chi)/r_chi)**3

      return      

      END FUNCTION r_suppression_func

    
      REAL(16) FUNCTION expEBeta(r)
      ! Auxillary function which gives exponential of gravitational
      ! potential times thermodynamic beta, ie exp(-E/kTw)
      ! Input:    r   radius [cm]
      ! Output:       exp(-E*Beta) [dimensionless] 

      double precision, intent(IN) :: r
      double precision :: dscapstarpot
      double precision :: eBeta
      real(16) :: eBeta_ext

      eBeta = mx * CGGEV * 1.d4 * (dscapstarpot(r*1.d-2) - dscapstarpot(0.d0))
     & / (BOLTZM * Tw)
      eBeta_ext = real(eBeta,16)
      expEBeta = exp(-eBeta_ext)
      if (expEBeta .lt. 2._16**minexponent(2._16)) expEBeta = 0._16

      return

      END FUNCTION expEBeta

    
      DOUBLE PRECISION FUNCTION tanhinv(x)
      ! Computes inverse hyperbolic tangent of x
      ! Input:    x   must have |x| < 1
      ! Output:       tanh^-1(x)

      double precision, intent(IN) :: x

      if (abs(x).ge.1) call die_quietly('Error: inverse hyperbolic tangent of
     & x with |x| >= 1 attempted.')
      
      tanhinv = 0.5d0*real(log((1._16+real(x,16))/(1._16-real(x,16))),8)

      return

      END FUNCTION tanhinv


      SUBROUTINE die_quietly(error_string)
      ! Kill program gracefully, citing error given by error_string

      character(*) :: error_string

      write(*,*)
      write(*,*) error_string
      write(*,*)

      stop

      END SUBROUTINE die_quietly


      DOUBLE PRECISION FUNCTION mach_tol_dbl()
	  
      double precision :: working_num, rtol
	  
      rtol = 1.d0
      working_num = 2.d0
      do while (working_num.gt.1.d0)
        store_temp = rtol + 1.d0
        working_num = store_temp
        rtol = rtol/2.d0
      enddo
      mach_tol_dbl = rtol*500.d0
	  
      return	  
	  
      END FUNCTION mach_tol_dbl


      SUBROUTINE setWIMPPop(tStep)
      ! Compute the new WIMP population on the basis of the last
      ! capture, annihilation (and evaporation) rates
      ! Input: next timestep to use (years)
	  
      double precision, intent(IN) :: tStep
      double precision :: equiv_t, equil_time, ann_parameter
	  
      if (abs(ann).le.prec) then
        n_WIMPs = n_WIMPs_prev + effcap*tStep
      else 
        ann_parameter = 2.d0*ann/n_WIMPs_prev**2
        if (abs(effcap).le.prec) then
          n_WIMPs = n_WIMPs_prev/(1.d0+n_WIMPs_prev*ann_parameter*tStep)
        else if (effcap.gt.0.d0) then
          equil_time = 1.d0/sqrt(effcap*ann_parameter)
          equiv_t = tanhinv((n_WIMPs_prev/effcap/equil_time)**sign(1.d0,dNdt))
          n_WIMPs = effcap*equil_time*(tanh(tStep/equil_time+equiv_t))**sign(1.d0,dNdt)  
        else
          !to include evaporation, insert code here for the case where evaporation outstrips capture
        endif
      endif

      return
	  
      END SUBROUTINE setWIMPPop


      END MODULE DkStrs_utils
