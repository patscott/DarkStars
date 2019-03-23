! Initialisation module for DarkStars
!
! Pat Scott, July 2007; pat@fysik.su.se
!--------------------------------------------------------------------------------      


      MODULE DkStrs_init
      use ez_driver
      use star_extras
      use DkStrs_data
      use DkStrs_utils
      use DkStrs_fint_ext

      implicit none

      contains


      SUBROUTINE init_capture
      ! Set stellar radius, mass and abundance structure for use by DarkSUSY in 
      ! determining capture rate.

      double precision :: dsntsunmfrac, dsntsundens
      double precision :: starmfr_remaining(meshpoints), heavy_numnorm, heavy_massnorm
      integer :: i, j

      call init_radmass
      
      ! Set centre values for radius, mass and elemental fractions
      starm(1)=0.0d0
      starr(1)=0.0d0
      starmfr(1,1)=center_H  ! for H1
      starmfr(2,1)=center_He ! for He4
      starmfr(3,1)=center_C  ! for C12
      starmfr(4,1)=center_N  ! for N14
      starmfr(5,1)=center_O  ! for O16
      starmfr(6,1)=center_Ne ! for Ne20
      starmfr(7,1)=center_Mg ! for Mg24
      
      ! Get stellar mass and radius in each shell
      call EZ_Extras
      starr(2:meshpoints) = SX(SX_R,SX_CNTR:SX_SURF:-1)*CRSN*1.d9/r_star
      starm(2:meshpoints) = SX(SX_M,SX_CNTR:SX_SURF:-1)/Star_mass

      ! Get EZ-derived mass fractions in each shell
      starmfr(1,2:meshpoints) = SX(SX_XH,SX_CNTR:SX_SURF:-1)
      starmfr(2,2:meshpoints) = SX(SX_XHE,SX_CNTR:SX_SURF:-1)
      starmfr(3,2:meshpoints) = SX(SX_XC,SX_CNTR:SX_SURF:-1)
      starmfr(4,2:meshpoints) = SX(SX_XN,SX_CNTR:SX_SURF:-1)
      starmfr(5,2:meshpoints) = SX(SX_XO,SX_CNTR:SX_SURF:-1)
      starmfr(6,2:meshpoints) = SX(SX_XNE,SX_CNTR:SX_SURF:-1)
      starmfr(7,2:meshpoints) = SX(SX_XMG,SX_CNTR:SX_SURF:-1)
      starmfr_remaining = 1.d0 - sum(starmfr(1:7,:), 1)

      if (legacy) then
        write(*,*)
        write(*,*) 'Legacy mode invoked: capture calculation to be'
        write(*,*) 'performed for the Sun regardless of stellar model.'
        write(*,*)
        do i=1,meshpoints 
          starmfr(1,i) = dsntsunmfrac(starr(i)*r_star,1) 
          starmfr(2,i) = dsntsunmfrac(starr(i)*r_star,2) 
          if (n_species.gt.16) starmfr(17,i) = dsntsunmfrac(starr(i)*r_star,3)
          starrho(i) = dsntsundens(starr(i)*r_star)
          do j=3,8
            starmfr(j,i) = dsntsunmfrac(starr(i)*r_star,j+1)
          end do
        end do
        starmfr(7,:) = starmfr(8,:)
        starmfr_remaining = 1.d0 - sum(starmfr(1:7,:), 1)
        if (n_species.gt.16) starmfr_remaining = starmfr_remaining - starmfr(17,:)
      endif
    
      ! Calculate remaining mass fractions in each shell
      forall(i=SOL_Na:SOL_Ni) starmfr(i+7,:) = 10.d0**heavy_abuns(i)
      heavy_numnorm = sum(starmfr(8:15,1))
      starmfr(15,:) = starmfr(15,:) * iso_ratios(ISO_58Ni60)/(1.d0+iso_ratios(ISO_58Ni60))
      starmfr(16,:) = starmfr(15,:) / iso_ratios(ISO_58Ni60)
      if (n_species.gt.16) then
        if (legacy) then
          starmfr(18,:) = 10.d0**Solar_C/iso_ratios(ISO_12C13)
          starmfr(19,:) = 10.d0**Solar_O/iso_ratios(ISO_16O18)
        else
          starmfr(17,:) = 1.d12*(star_Mass_He/starma_atomic(2)) / (star_Mass_H/starma_atomic(1))
     &     /iso_ratios(ISO_4He3)
          starmfr(18,:) = 1.d12*(star_Mass_C/starma_atomic(3)) / (star_Mass_H/starma_atomic(1))
     &     /iso_ratios(ISO_12C13)
          starmfr(19,:) = 1.d12*(star_Mass_O/starma_atomic(5)) / (star_Mass_H/starma_atomic(1))
     &     /iso_ratios(ISO_16O18)
        endif
        starmfr(20,:) = 10.d0**heavy_abuns(SOL_Pb)
        heavy_numnorm=heavy_numnorm+sum(starmfr(17:20,1))
        if (legacy) heavy_numnorm = heavy_numnorm - starmfr(17,1)
        starmfr(20,:) = starmfr(20,:)*
     &   (iso_ratios(ISO_208Pb207)+iso_ratios(ISO_208Pb206))/
     &   (1.d0+iso_ratios(ISO_208Pb207)+iso_ratios(ISO_208Pb206))
        starmfr(21,:) = starmfr(20,:) / iso_ratios(ISO_208Pb207)
        starmfr(22,:) = starmfr(20,:) / iso_ratios(ISO_208Pb206)
      endif
      forall (i=8:n_species) starmfr(i,:) = starmfr(i,:)/heavy_numnorm*starma_atomic(i)
      heavy_massnorm = sum(starmfr(8:n_species,1))
      if (n_species.gt.16 .and. legacy) heavy_massnorm = heavy_massnorm - starmfr(17,1)
      starmfr(8:n_species,:) = starmfr(8:n_species,:)/heavy_massnorm
      forall (i=8:n_species) starmfr(i,:) = starmfr(i,:) * starmfr_remaining
      
      !Convert from mass fraction to number density (per cubic cm) for DarkSUSY 
      if (legacy) then
        if (n_species.gt.16) then
          do i=1,meshpoints 
            starmfr(17,i) = dsntsunmfrac(starr(i)*r_star,3)
          end do
        endif
        forall(i = 1:n_species) starmfr(i,:) = starmfr(i,:)*
     &   starrho/starma_atomic(i)/1.d3/CMEV*CL**2 
      else
        forall(i = 1:n_species) starmfr(i,2:meshpoints) = starmfr(i,2:meshpoints)*
     &   SX(SX_RHO,SX_CNTR:SX_SURF:-1)/starma_atomic(i)/CGGEV
        starmfr(:,1) = starmfr(:,1)*10.d0**log_center_Density/starma_atomic/CGGEV
      endif    

      END SUBROUTINE init_capture


      SUBROUTINE init_radmass
      ! Get total stellar mass and radius
      m_star = Star_mass * CMSN * 1.d30         ! convert to kg for DarkSUSY
      r_star = 10.d0**log_Radius * CRSN * 1.d9  ! convert to m for DarkSUSY
      
      END SUBROUTINE init_radmass

      
      SUBROUTINE init_annihilation(before_evolve)
      ! Initialise WIMP radial density normalisation constants, alphas,
      ! temperatures, Knudsen factors and LTE WIMP densities.
      use ez_cycle_data

      logical, intent(IN) :: before_evolve
      double precision :: rho_eff, dsLEfint

      if (tau_therm .gt. prec .and. .not.before_evolve .and. JMOD.ne.0) then
        r_chi = (cap/CSY*sigann*tau_therm**2.d0)**(1.d0/3.d0)/pi**(0.5d0)
      else
        r_chi = sqrt(3.d0*BOLTZM*10.0d0**log_Center_Temp/(2.d0*pi*CG*
     &   10.d0**log_Center_Density*mx*CGGEV))
      endif      

      if (before_evolve) then
        L_WIMPtransport = 0.0d0
      else
        if (WIMPdensmode.or.DoTransport) then
          call alpha_init
          call temperature_init
          K_suppression = K_suppression_func(mfp(1)/r_chi)
        endif
        if (WIMPdensmode) then
          call WIMPdens_LTE_init
          call Find_Tw
          idtag = 'partition_f2'
          partition_func = fint_ext(partition_function_integrand,0._16,1.e2_16*r_star,1.e-6_16)
          if (partition_func.eq.0._16) then
            write(*,*) 'COMMENT: partition_func in init_annihilation being recomputed'
            partition_func = fint_ext(partition_function_integrand,0._16,1.e2_16*r_star/32._16,1.e-6_16)
            if (partition_func.eq.0._16) call die_quietly('Error: partition_func = 0 in init_annihilation.')       
          endif
        else if (tau_therm .gt. prec) then
          stardens(1) = 10.d0**log_center_Density
          stardens(2:meshpoints) = SX(SX_RHO,SX_CNTR:SX_SURF:-1)
          call dsspline(starr,stardens,meshpoints,0.0d0,0.0d0,stardensd2)
          idtag = 'rho_eff'
          rho_eff = 3.d0 * dsLEfint(rho_eff_integrand, 0.d0, 1.d0, 1.d-6)
          if (rho_eff.le.0.d0) call die_quietly('Error: rho_eff <= 0 in init_annihilation.')
          Tw = 2.d0 * pi * CG * rho_eff * mx*CGGEV *r_chi**2.d0 / (3.d0 * BOLTZM)
        endif

      endif

      END SUBROUTINE init_annihilation

   
      SUBROUTINE init_transport
      ! Populate the table of luminosities carried by WIMPs via conduction, and
      ! initialise the interpolator.
      use DkStrs_WIMPdens

      double precision :: r, transport_density, working_array(2*meshpoints-2)
      integer :: exit_code, i

      if (LTE_invalid) then
      
        L_WIMPtransport = 0.d0
      
      else

        do i=1,meshpoints
          r = starr(i) * r_star * 1.d2
          if (WIMPdensmode) then
            transport_density = real(WIMPdens_LTE(r),8)
          else
            transport_density = WIMPdens(r)
          endif
          !L_WIMPtransport(i) = 4.d0 * pi * r**2 * kappa(i) *
     &    ! transport_density * mfp(i) * abs(tempgrad(r)) *
     &    ! sqrt(BOLTZM * startemperature(i) / mx / CMEV * 1.d-3) *
     &    ! K_suppression * r_suppression_func(r)

          ! Corrected bug in transport scaling by multiplying by k*c   PS 090728
          ! Note that the c^2 in Eq. 2.31 in Scott et al (2009, MNRAS 394:82) is erroneous 
          L_WIMPtransport(i) = 4.d0 * pi * r**2 * kappa(i) *
     &     transport_density * mfp(i) * abs(tempgrad(r)) * BOLTZM *
     &     sqrt(BOLTZM * startemperature(i) / mx / CGGEV) *
     &     K_suppression * r_suppression_func(r)
        enddo
      
        L_WIMPtransportd1(1:meshpoints:meshpoints-1) = 0.d0
        if (interpmode) call TSPSI (meshpoints,starr,L_WIMPtransport,2,1,
     &   .false.,.false.,2*meshpoints-2,working_array,L_WIMPtransportd1,L_WIMPtransporttension,exit_code)
        if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in init_transport...')

      endif

      END SUBROUTINE init_transport


      SUBROUTINE alpha_kappa_precalc
      ! For each nuclear species considered, precalculate the WIMP thermal diffusivity
      ! and conductivity (alpha & kappa) for a gas consisting of just WIMPs and that nucleus 
      ! (see dscapsetup in DarkSUSY for list of nuclei).
      use Dkstrs_akraw

      integer :: i
      double precision :: alpha_rawd2(n_pts), kappa_rawd2(n_pts), current_ratio

      call dsspline(mass_ratio,alpha_raw,n_pts,0.0d0,0.0d0,alpha_rawd2)
      call dsspline(mass_ratio,kappa_raw,n_pts,1.d32,5.d0*sqrt(pi/60.d0)/64.d0,kappa_rawd2)

      do i=1,n_species
        current_ratio=mx/starma(i)
        if (current_ratio.gt.mass_ratio(n_pts)) then
          alpha_i(i) = 2.5d0
          kappa_i = sqrt(2.d0*pi*current_ratio)*5.d0/32.d0
        else
          call dssplint(mass_ratio,alpha_raw,alpha_rawd2,n_pts,current_ratio,alpha_i(i))
          call dssplint(mass_ratio,kappa_raw,kappa_rawd2,n_pts,current_ratio,kappa_i(i))
        endif
      enddo

      END SUBROUTINE alpha_kappa_precalc


      SUBROUTINE sigma_precalc
      ! Precalculate WIMP-nucleus scattering cross-section for each nucleus (see dscapsetup in
      ! DarkSUSY for list of nuclei). 

      sigma_i=sigsi*staraa**2*(mx*starma/(mx+starma))**2
     & /(mx*CMP/CGGEV/(mx+CMP/CGGEV))**2

      sigma_i(1) = sigma_i(1) + sigsd

      END SUBROUTINE sigma_precalc


      SUBROUTINE orbit_init
      ! Read in the orbital details from orbit.dat, put them in a table and initialise the interpolator
      
      integer, parameter :: lun = 100
      integer :: i, fileStatus, allocateStatus(13), exit_code1, exit_code2, exit_code3=0, exit_code4=0
      double precision, allocatable :: working_array(:)
      character (len=strlen) :: orbitFile

      orbitFile = trim(main_dir)//'/orbit.dat'
      open(UNIT=lun,FILE=orbitFile,ACCESS='SEQUENTIAL',ACTION='READ',IOSTAT=fileStatus)
      if (fileStatus.gt.0) call die_quietly('  Error: Cannot find file '//orbitFile)
      
      read(UNIT=lun,FMT=*,IOSTAT=fileStatus) fileLen
      read(UNIT=lun,FMT=*,IOSTAT=fileStatus) orbit_monotonicity_scale
      if (fileStatus.ne.0) call die_quietly('  Error: orbit.dat contains badly formatted header data.')
      allocate(orbital_t(fileLen),STAT=allocateStatus(1))
      allocate(orbital_rho(fileLen),STAT=allocateStatus(2))
      allocate(orbital_rho_d1(fileLen),STAT=allocateStatus(3))
      allocate(orbital_rho_tension(fileLen),STAT=allocateStatus(4))
      allocate(orbital_v(fileLen),STAT=allocateStatus(5))
      allocate(orbital_v_d1(fileLen),STAT=allocateStatus(6))
      allocate(orbital_v_tension(fileLen),STAT=allocateStatus(7))
      allocate(orbital_galesc(fileLen),STAT=allocateStatus(8))
      allocate(orbital_galesc_d1(fileLen),STAT=allocateStatus(9))
      allocate(orbital_galesc_tension(fileLen),STAT=allocateStatus(10))
      allocate(orbital_galr(fileLen),STAT=allocateStatus(11))
      allocate(orbital_galr_d1(fileLen),STAT=allocateStatus(12))
      allocate(orbital_galr_tension(fileLen),STAT=allocateStatus(13))
      if(any(allocateStatus(1:13).gt.0)) call die_quietly('  Error: memory allocation failed in orbit_init.')
      
      do i=1,fileLen
        if (galesc .lt. 0.d0) then
          read(UNIT=lun,FMT=*,IOSTAT=fileStatus) orbital_t(i), orbital_rho(i), orbital_v(i), 
     &     orbital_galesc(i), orbital_galr(i)
          orbital_galr(i) = orbital_galr(i)/CMPC
        else
          read(UNIT=lun,FMT=*,IOSTAT=fileStatus) orbital_t(i), orbital_rho(i), orbital_v(i)
        endif
      enddo
      if (fileStatus.ne.0) call die_quietly('  Error: orbit.dat contains badly formatted body data.')
      close(lun,IOSTAT=fileStatus)
      if (fileStatus.ne.0) call die_quietly('  Error: problem closing orbit.dat.')

      loopOrbit = (orbital_rho(1).eq.orbital_rho(fileLen)).and.(orbital_v(1).eq.orbital_v(fileLen))
      orbit_len = orbital_t(fileLen) - orbital_t(1)
      
      if (interpmode) then
      
        if(loopOrbit) then
            
          allocate(working_array(3*fileLen-3),STAT=allocateStatus(1))
          if (allocateStatus(1).gt.0) call die_quietly('  Error: memory allocation failed in orbit_init.')
          call TSPSI (fileLen,orbital_t,orbital_rho,2,2,.true.,.false.,3*fileLen-3,working_array,
     &      orbital_rho_d1,orbital_rho_tension,exit_code1)
          call TSPSI (fileLen,orbital_t,orbital_v,2,2,.true.,.false.,3*fileLen-3,working_array,
     &      orbital_v_d1,orbital_v_tension,exit_code2)
          if (galesc .lt. 0.d0) then
            orbital_galesc(fileLen) = orbital_galesc(1)
            orbital_galr(fileLen) = orbital_galr(1)
            call TSPSI (fileLen,orbital_t,orbital_galesc,2,2,.true.,.false.,3*fileLen-3,working_array,
     &       orbital_galesc_d1,orbital_galesc_tension,exit_code3)
            call TSPSI (fileLen,orbital_t,orbital_galr,2,2,.true.,.false.,3*fileLen-3,working_array,
     &       orbital_galr_d1,orbital_galr_tension,exit_code4)
          endif          

        else

          allocate(working_array(2*fileLen-2),STAT=allocateStatus(1))
          if (allocateStatus(1).gt.0) call die_quietly('  Error: memory allocation failed in orbit_init.')
          orbital_rho_d1(1:fileLen:fileLen-1) = 0.d0
          orbital_v_d1(1:fileLen:fileLen-1) = 0.d0
          call TSPSI (fileLen,orbital_t,orbital_rho,2,1,.false.,.false.,2*fileLen-2,working_array,
     &     orbital_rho_d1,orbital_rho_tension,exit_code1)
          call TSPSI (fileLen,orbital_t,orbital_v,2,1,.false.,.false.,2*fileLen-2,working_array,
     &     orbital_v_d1,orbital_v_tension,exit_code2)
          if (galesc .lt. 0.d0) then
            orbital_galesc_d1(1:fileLen:fileLen-1) = 0.d0
            orbital_galr_d1(1:fileLen:fileLen-1) = 0.d0
            call TSPSI (fileLen,orbital_t,orbital_galesc,2,1,.false.,.false.,2*fileLen-2,working_array,
     &       orbital_galesc_d1,orbital_galesc_tension,exit_code3)
            call TSPSI (fileLen,orbital_t,orbital_galr,2,1,.false.,.false.,2*fileLen-2,working_array,
     &       orbital_galr_d1,orbital_galr_tension,exit_code4)
          endif

        endif
        
        if((exit_code1.lt.0).or.(exit_code2.lt.0).or.(exit_code3.lt.0).or.(exit_code4.lt.0)) 
     &   call die_quietly('Dying quietly - tensional splines snapped in orbit_init...')

      else
      
        !Not implemented yet

      endif
      
      if (loopOrbit) then 
        rhowimp = orbital_rho(fileLen/2)
        v_star = orbital_v(fileLen/2)
        if (galesc .lt. 0.d0) then
          galesc = orbital_galesc(fileLen/2)
          galr = orbital_galr(fileLen/2)
        endif
      else
        rhowimp = orbital_rho(1)
        v_star = orbital_v(1)
        if (galesc .lt. 0.d0) then
          galesc = orbital_galesc(1)
          galr = orbital_galr(1)
        endif
      endif
      
      END SUBROUTINE orbit_init


      SUBROUTINE alpha_init
      ! Populate the table of alpha values throughout the star, and initialise the interpolator.

      integer :: i

      forall (i=1:meshpoints) staralpha(i) = sum(sigma_i*starmfr(:,i)*alpha_i)*mfp(i)
      if (interpmode) call dsspline(starr,staralpha,meshpoints,0.d0,1.d32,staralphad2)

      END SUBROUTINE alpha_init


      SUBROUTINE temperature_init
      ! Populate the table of temperature values throughout the star, and initialise the interpolator.
      
      double precision :: working_array(2*meshpoints-2)
      integer :: exit_code, i

      startemperature(1) = 10.0d0**log_Center_Temp
      startemperature(2:meshpoints) = SX(SX_T,SX_CNTR:SX_SURF:-1)
      startemperatured1(1:meshpoints:meshpoints-1) = 0.d0
      if (interpmode) call TSPSI (meshpoints,starr,startemperature,2,2,
     & .false.,.false.,2*meshpoints-2,working_array,startemperatured1,temptension,exit_code)
      if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in temperature_init...')

      END SUBROUTINE temperature_init      
      

      SUBROUTINE WIMPdens_LTE_init
      ! Populate the table of WIMP densities throughout the star assuming WIMP-nucleus
      ! LTE, and initialise the interpolator.      

      real(16) :: centre_LTE_normfactor
      double precision :: tempradius, dsLEfint
      double precision :: working_array(2*meshpoints-2), tempintegral
      integer :: exit_code, i

      starLTErhoWIMP(1) = (temperature(0.d0)/startemperature(1))**(1.5d0)  
      do i=2,meshpoints
        tempradius = starr(i)*r_star*1.d2
        starLTErhoWIMP(i) = (temperature(tempradius)/startemperature(1))**(1.5d0)
        idtag = 'WIMPdens_LTE'
        tempintegral = dsLEfint(LTEdens_integrand,0.d0,tempradius,1.d-6)
        if (tempintegral.eq.0.d0) then 
          write(*,*) 'COMMENT: tempintegral in WIMPdens_LTE_init being recomputed'
          tempintegral = dsLEfint(LTEdens_integrand,0.d0,tempradius/32.d0,1.d-6)
          if (tempintegral.eq.0.d0) call die_quietly('Error: tempintegral(r) = 0 in WIMPdens_LTE_init for r <> 0.')
        endif
        starLTErhoWIMP(i) = starLTErhoWIMP(i)*exp(-tempintegral)
      enddo

      n_WIMPs_centre_LTE = 1._16
      starLTErhoWIMPd1(1:meshpoints:meshpoints-1) = 0.d0
      if (interpmode) call TSPSI (meshpoints,starr,starLTErhoWIMP,2,1,
     & .false.,.false.,2*meshpoints-2,working_array,starLTErhoWIMPd1,LTErhoWIMPtension,exit_code)
      if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped at first try
     & in WIMPdens_LTE_init...')
       
      LTE_invalid = .false.
      idtag = 'n_WIMPs_cent'
      centre_LTE_normfactor = fint_ext(LTEdens_norm_integrand,0._16,real(r_star,16)*1.e2_16,1.e-6_16)
      if (centre_LTE_normfactor.eq.0._16) then
        centre_LTE_normfactor = fint_ext(LTEdens_norm_integrand,0._16,real(r_star,16)*1.e2_16,1.e-2_16)
        if (centre_LTE_normfactor.eq.0._16) then  
          write(*,*) 'WARNING: integral for LTE WIMP distribution normalisation factor'
          write(*,*) 'did not converge. Reverting to less sophisticated isothermal density'
          write(*,*) 'approximation for this timestep.'
          LTE_invalid = .true.
        endif
      endif
      n_WIMPs_centre_LTE = real(n_WIMPs,16)/centre_LTE_normfactor
      
      !Note that there is a subtle sleight of hand here in using extended precision for the LTE density normalisation
      !factor but not for the actual densities.  Implicit in this choice is the assumption that the density expression
      !goes to zero more quickly than the normalisation constant as R*-->infinity.  This means that as r-->R* one always
      !has rho_LTE-->0 for any choice of R* - assuming the normalisatioin factor integral converges.
      
      END SUBROUTINE WIMPdens_LTE_init


      SUBROUTINE Find_Tw
      ! Calculate the value of Tw.  Tw is the WIMP temperature at which the (incorrect) WIMP
      ! conductive energy transport scheme of Spergel & Press (1985, ApJ 294:663) produces no
      ! net energy flux into or out of the star.  This treatment draws upon the isothermal sphere
      ! approximation for the WIMP radial density, and sets Tw as the characteristic temperature
      ! of the isothermal sphere.

      double precision :: tempmax, tempmin, zeroin, SPintegral_temp
      integer :: tempmaxloc(1), i

      tempmaxloc = maxloc(startemperature)
      tempmax = startemperature(tempmaxloc(1))
      SPintegral_temp = SPintegral(tempmax)

      if (SPintegral_temp.gt.-1.d0*prec) then
        Tw = tempmax
      else            
        i=tempmaxloc(1)
        do while (SPintegral(startemperature(i)).lt.0.d0)
          i = i+1
          if (i.gt.meshpoints) call die_quietly('Error: Tw search ran too far in init_annihilation')
        enddo
        tempmin = startemperature(i)
        if (SPintegral(tempmin).eq.0.d0) then
          Tw = tempmin
        else
          if (i.ne.tempmaxloc(1)) tempmax = startemperature(i-1)
          Tw = zeroin(tempmin,tempmax,SPintegral,1.d-3)
        endif
      endif

      END SUBROUTINE Find_Tw

      
      DOUBLE PRECISION FUNCTION SPintegral(T_WIMP)
      ! Integral giving net extrastellar luminosity carried by WIMPs, according to the
      ! treatment of Spergel & Press (1985, ApJ 294:663).  Should be equal to zero if Tw
      ! is chosen correctly (and there is no WIMP evaporation).
      ! Input:    Tw  WIMP isothermal temperature [K]
      ! Output:       Luminosity carried away by WIMPs (energy stolen by WIMPs per
      !               second) from nuclear matter [erg/s (* a further (K/erg)^3/2 due to commenting out
      !               of the leading constant in SPintegrand] 

      double precision, intent(IN) :: T_WIMP
      real(16) :: SPintegral_ext
      integer :: i

      Tw = T_WIMP
      idtag = 'partition_f1'
      partition_func = fint_ext(partition_function_integrand,0._16,1.e2_16*real(r_star,16),1.e-3_16)
      if (partition_func.eq.0._16) then
        write(*,*) 'COMMENT: partition_func in SPintegral being recomputed'
        partition_func = fint_ext(partition_function_integrand,0._16,1.e2_16*real(r_star,16)/32._16,1.e-3_16)
        if (partition_func.eq.0._16) call die_quietly('Error: partition_func = 0 in SPintegral.')       
      endif
      
      idtag = 'SPintegral'
      SPintegral_ext = fint_ext(SPintegrand,0._16,1.e2_16*real(r_star,16),1.e-2_16)
      if (SPintegral_ext.eq.0._16) write(*,*) 'COMMENT: SPintegral did not converge for chosen Tw'
      !Note that due to the sharply-peaked nature of the integrand around the point of maximum temperature, 
      !the integral sometimes will only converge VERY slowly once it starts to get close to zero. This typically
      !happens only in fairly extreme cases, like ultra-low WIMP densities in giant stars burning H/He in shells.
      !To circumvent this, I've just let the integrator return zero if it finds itself running past its
      !maximum number of iterations when called on this function. PS 
        
      SPintegral = real(SPintegral_ext,8)
      if (SPintegral_ext .ne. 0._16 .and. SPintegral .eq. 0.d0) write(*,*) 'COMMENT: truncation to zero in SPintegral'

      END FUNCTION SPintegral


      REAL(16) FUNCTION SPintegrand(r_in)
      ! Integrand of above SPintegral
      ! Input:    r   radius [cm]
      ! Output:       integrand [erg/cm/s (* a further (K/erg)^3/2 due to commenting out
      !                          of the leading constant] 

      real(16), intent(IN) :: r_in
      double precision :: r, dscapstardenscomp
      integer :: i
      
      r = real(r_in,8)
      SPintegrand = 0._16
      do i=1,n_species
        SPintegrand = SPintegrand + real(dscapstardenscomp(r*1.d-2,i)*sigma_i(i)*mx*starma(i)
     &   /(mx+starma(i))*sqrt(temperature(r)/(starma(i)*CGGEV)+Tw/(mx*CGGEV)),16)
      if (abs(SPintegrand).gt.huge(SPintegrand) .or. isnan(SPintegrand)) then
        write(*,*) SPintegrand, dscapstardenscomp(r*1.d-2,i), temperature(r)/(starma(i)*CGGEV)+Tw/(mx*CGGEV)
          call die_quietly('Error: SPintegrand = NaN or Inf in middle of SPintegrand')
      endif
      enddo
      SPintegrand = SPintegrand*real(r**2*(temperature(r)-Tw),16)*WIMPdens_isothermal(r)
            
      if (abs(SPintegrand).gt.huge(SPintegrand) .or. isnan(SPintegrand)) then
        write(*,*) r*1.d-2,temperature(r),Tw,WIMPdens_isothermal(r)
        call die_quietly('Error: SPintegrand = NaN or Inf at end of SPintegrand')
      endif
      
      !Uncomment the following line to get the dimensionally correct value,
      !which is necessary if extending this treatment to incorporate evaporation.
      !SPintegrand = SPintegrand * 32._16*sqrt(2._16*real(pi,16))*real(BOLTZM,16)**(1.5_16)

      END FUNCTION SPintegrand


      DOUBLE PRECISION FUNCTION LTEdens_integrand(r)
      ! Integrand for exponential factor in WIMP density formula under
      ! the assumption of WIMP-nuclei LTE.
      ! Input:    r   radius [cm]
      ! Output:       integrand [cm^-1] 

      double precision, intent(IN) :: r
      double precision :: dscapspfunc
      
      LTEdens_integrand = (alpha(r)*tempgrad(r)*BOLTZM + mx*CGGEV*1.d2*dscapspfunc(r*1.d-2))
     & / (BOLTZM*temperature(r))

      return

      END FUNCTION LTEdens_integrand


      REAL(16) FUNCTION LTEdens_norm_integrand(r)
      ! Integrand for normalisation factor for WIMP density under
      ! the assumption of WIMP-nuclei LTE.
      ! Input:    r   radius [cm]
      ! Output:       integrand [cm^2] 

      real(16), intent(IN) :: r

      LTEdens_norm_integrand = 4._16 * real(pi,16) * r**2 * WIMPdens_LTE(real(r,8))

      return

      END FUNCTION LTEdens_norm_integrand


      REAL(16) FUNCTION partition_function_integrand(r)
      ! Integrand for gravitional partition function
      ! Input:    r   radius [cm]
      ! Output:       integrand [cm^2] 

      real(16), intent(IN) :: r

      partition_function_integrand = 4._16 * real(pi,16) * r**2 * expEBeta(real(r,8))

      return

      END FUNCTION partition_function_integrand


      REAL FUNCTION rho_eff_integrand(r)
      ! Integrand for effective stellar density in the region
      ! populated by WIMPs
      ! Input:    r   radius [dimensionless]
      ! Output:       integrand [g cm^-3] 

      double precision, intent(IN) :: r
      double precision :: local_dens

      call dssplint(starr,stardens,stardensd2,meshpoints,r,local_dens)
      rho_eff_integrand = r**2 * local_dens

      END FUNCTION rho_eff_integrand


      END MODULE DkStrs_init
