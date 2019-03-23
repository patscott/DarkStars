! WIMP radial number density module for DarkStars
!
! Pat Scott, July 2007; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE DkStrs_WIMPdens
      use ez_cycle_data
      use DkStrs_data
      use DkStrs_utils
      implicit none
   
      contains
      

      DOUBLE PRECISION FUNCTION WIMPdens(r)

      ! WIMP radial density calculator
      ! Input:    r   radius [cm]
      ! Output:       WIMP number density [cm^-3] 
      ! defaults to approximate method during EZ initialisation 
      ! and after first timestep (JMOD=0)

      double precision, intent(IN) :: r
	  
      if (r.gt.r_star*1.d2 .or. .not.DoWimpyThings) then
        WIMPdens=0.d0
        return
      endif
      	  
      if (WIMPdensmode .and. .not.new_model_step1 .and. JMOD.ge.1) then
        !work out WIMP density on the basis of density structure

        if(LTE_invalid) then
		
          WIMPdens = real(WIMPdens_isothermal(r),8)
		
        else

          WIMPdens = K_suppression*real(WIMPdens_LTE(r),8) + (1.d0-K_suppression)*real(WIMPdens_isothermal(r),8)
		
        endif
          
      else
        !use approximate expression based on assumption of constant stellar density
        !in the region populated by WIMPS (that is, use the star's central density unless tau_therm has been given).        
   
        if (r.le.0.0d0) then
          WIMPdens=pi**(-1.5d0)*n_WIMPs/r_chi**3
          return
        endif
     
        WIMPdens = pi**(-1.5d0)*n_WIMPs/r_chi**3*exp(-r**2/r_chi**2)
        
      endif

      if (WIMPdens.lt.0.d0) call die_quietly('Error: Negative WIMP density encountered (EZ probably 
     & called WIMPdens(NaN)...)')
      
      return

      END FUNCTION WIMPdens  

      
      END MODULE DkStrs_WIMPdens
