! WIMP conductive energy transport module for DarkStars
!
! Pat Scott, July 2007; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE DkStrs_transport
      use DkStrs_data
      use DkStrs_utils
      implicit none

      contains

      
      DOUBLE PRECISION FUNCTION E_transport(r,rho)
      ! Takes local height and density and returns local energy yield of
      ! WIMP conductive heat transport, using gradient of the previously 
      ! interpolated function for the WIMP-mediated luminosity.  
      ! Input:    r   radius [cm]
      !          rho  density [g/cm^3]
      ! Output:       energy generation rate from conductive heat transfer by WIMPs [ergs/g/s]
      use ez_cycle_data  

      double precision, intent(IN) :: r,rho
      double precision :: HPVAL
      integer :: exit_code

      if (sum(L_WIMPtransport).eq.0.d0 .or. .not.DoTransport
     & .or. rho.eq.0.d0 .or. r.eq.0.d0 .or. .not.DoWimpyThings) then
        E_transport = 0
        return
      endif

      if (interpmode) then

        E_transport = HPVAL(r/r_star*1.d-2,meshpoints,starr,L_WIMPtransport,
     &   L_WIMPtransportd1,L_WIMPtransporttension,exit_code)/r_star*1.d-2
        if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in E_transport...')
      
      else
      
      endif

      E_transport = -E_transport/(4.d0*pi*r**2*rho)

      return

      END FUNCTION E_transport

      
      SUBROUTINE get_condEff
      ! Calculate the dimensionless conductive transport effectiveness.
      use star_extras
      use DkStrs_fint_ext
      
      double precision :: r,working_array(2*meshpoints-2), condEff_norm
      integer :: exit_code, i

      if (sum(starcondEff).eq.0.d0 .or. .not.DoTransport .or. .not.DoWimpyThings) then
        condEff = 0.d0
        return
      endif

      starcondEff(1) = 0.d0
      starcondEff(2:meshpoints) = SX(SX_RHO,SX_CNTR:SX_SURF:-1) / SX(SX_MU,SX_CNTR:SX_SURF:-1)
      starcondEff = starcondEff * starr**2
      
      starcondEffd1(1:meshpoints:meshpoints-1) = 0.d0
      if (interpmode) call TSPSI (meshpoints,starr,starcondEff,2,1,
     & .false.,.false.,2*meshpoints-2,working_array,starcondEffd1,condEfftension,exit_code)
      if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped at get_condEff, part 1...')

      idtag = 'condEff_norm'
      condEff_norm = real(fint_ext(condEff_integrand,0._16,1._16,1.e-6_16),8)        

      do i=2,meshpoints
        r = starr(i) * r_star * 1.d2
        starcondEff(i) = starcondEff(i) * abs(E_transport(r,SX(SX_RHO,1+meshpoints-i)) / 
     &   (SX(SX_EPS_NUC,1+meshpoints-i) + SX(SX_EPS_G,1+meshpoints-i)))
      enddo

      starcondEffd1(1:meshpoints:meshpoints-1) = 0.d0
      if (interpmode) call TSPSI (meshpoints,starr,starcondEff,2,1,
     & .false.,.false.,2*meshpoints-2,working_array,starcondEffd1,condEfftension,exit_code)
      if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in get_condEff, part 2...')

      idtag = 'condEff'
      condEff = real(fint_ext(condEff_integrand,0._16,1._16,1.e-6_16),8) / condEff_norm
       
      END SUBROUTINE get_condEff


      REAL(16) FUNCTION condEff_integrand(r)
      ! Integrand for the conductive transport effectiveness
      ! number.
      ! Input:    r   radius [dimensionless]
      ! Output:       integrand [cm^-3]

      real(16), intent(IN) :: r      
      double precision :: HVAL
      integer :: exit_code

      condEff_integrand = real(HVAL(real(r,8),meshpoints,starr,starcondEff,
     & starcondEffd1,condEfftension,exit_code),16)
      if (exit_code.lt.0) call die_quietly('Dying quietly - tensional splines snapped in condEff_integrand...')


      END FUNCTION condEff_integrand 


      END MODULE DkStrs_transport
