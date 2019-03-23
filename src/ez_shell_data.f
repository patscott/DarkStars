      MODULE ez_shell_data ! data for a single shell.  also see ez_state_data for EoS info for the shell.
      USE star_data
      IMPLICIT NONE 
      
      DOUBLE PRECISION :: VAR(V_MAX) ! holds the independent variables for the current shell.
      DOUBLE PRECISION :: DVAR(V_MAX) ! holds changes in them.
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      ! These are the independent variables and their changes since the last model.
      DOUBLE PRECISION :: LNF         ! ln(F).  degeneracy parameter PSI is a function of LNF.
      DOUBLE PRECISION :: DLNF          ! change in ln(F) during the previous timestep.
      DOUBLE PRECISION :: LNT         ! ln(temperature T).
      DOUBLE PRECISION :: DLNT          ! change in ln(T) during the previous timestep.
      DOUBLE PRECISION :: M             ! mass in 10^33 g units interior to the current mesh point.
      DOUBLE PRECISION :: DM           ! change in stellar mass in 10^33 g units during the previous timestep.
      DOUBLE PRECISION :: Q_dk            ! derivative wrt meshpoint (K) of meshpoint spacing function Q.
      DOUBLE PRECISION :: LNR         ! ln(radius R in 10^11 cm units).
      DOUBLE PRECISION :: L             ! luminosity in 10^33 ergs per second.
      DOUBLE PRECISION :: XH           ! fractional abundance by mass of H1.
      DOUBLE PRECISION :: DXH         ! change in fractional abundance by mass of H1 during the previous timestep.
      DOUBLE PRECISION :: XHE         ! fractional abundance by mass of He4.
      DOUBLE PRECISION :: DXHE        ! change in fractional abundance by mass of He4 during the previous timestep.
      DOUBLE PRECISION :: XC           ! fractional abundance by mass of C12.
      DOUBLE PRECISION :: DXC         ! change in fractional abundance by mass of C12 during the previous timestep.
      DOUBLE PRECISION :: XO           ! fractional abundance by mass of O16.
      DOUBLE PRECISION :: DXO         ! change in fractional abundance by mass of O16 during the previous timestep.
      DOUBLE PRECISION :: XNE         ! fractional abundance by mass of Ne20.
      DOUBLE PRECISION :: DXNE      ! change in fractional abundance by mass of Ne20 during the previous timestep.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! These are functions of the independent variables 
      DOUBLE PRECISION :: XHT          ! dXH/dt.
      DOUBLE PRECISION :: XHET      ! dXHE/dt.
      DOUBLE PRECISION :: XCT         ! dXC/dt.
      DOUBLE PRECISION :: XOT         ! dXO/dt.
      DOUBLE PRECISION :: XNET     ! dXNE/dt.
      DOUBLE PRECISION :: BC_Psurf ! surface pressure boundary condition.
      DOUBLE PRECISION :: BC_Tsurf ! surface temperature boundary condition.
      DOUBLE PRECISION :: BC_Msurf  ! boundary condition for stellar mass.
      DOUBLE PRECISION :: PQ        ! pressure term for mesh spacing function Q.
      DOUBLE PRECISION :: PQ_dk           ! dPQ/dk. 
      DOUBLE PRECISION :: RQ        ! radius term for mesh spacing function Q.
      DOUBLE PRECISION :: RQ_dk        ! dRQ/dk.
      DOUBLE PRECISION :: TQ        ! temperature term for mesh spacing function Q.
      DOUBLE PRECISION :: TQ_dk        ! dTQ/dk.
      DOUBLE PRECISION :: LQ        ! luminosity term for mesh spacing function Q.
      DOUBLE PRECISION :: LQ1_dk          ! the part of dLQ/dk that is not multiplied by dM/dt.
      DOUBLE PRECISION :: LQ2_dk          ! the part of dLQ/dk that is multiplied by dM/dt.
      DOUBLE PRECISION :: M_dt       ! dM/dt.  change in mass coordinate of this meshpoint since the last timestep
      DOUBLE PRECISION :: MQ        ! mass term for mesh spacing function Q.
      DOUBLE PRECISION :: MQ_dk        ! dMQ/dk.
      DOUBLE PRECISION :: SG           ! diffusion coefficient for convective mixing.
      DOUBLE PRECISION :: WT           ! weighting term for evaluation of "average" values from adjacent meshpoint values.
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! These are indices in the FN array.  FN holds values that are needed to evaluate the equations.
      ! When the solver needs to create a Jacobian of partial derivatives of equation wrt variables,
      ! it does it by first collecting all the FN derivatives wrt variables, and then it gets the EQN derivatives wrt the FNs.
      INTEGER, PARAMETER :: F_XH=17
      INTEGER, PARAMETER :: F_XH_dt=18 ! dXH/dt.
      INTEGER, PARAMETER :: F_XHE=21
      INTEGER, PARAMETER :: F_XHE_dt=22   ! dXHE/dt.
      INTEGER, PARAMETER :: F_XC=23
      INTEGER, PARAMETER :: F_XC_dt=24 ! dXC/dt.
      INTEGER, PARAMETER :: F_XO=19
      INTEGER, PARAMETER :: F_XO_dt=20 ! dXO/dt.
      INTEGER, PARAMETER :: F_XNE=25
      INTEGER, PARAMETER :: F_XNE_dt=26 ! dXNE/dt.
      INTEGER, PARAMETER :: F_BC_Psurf=1 ! surface pressure boundary condition.
      INTEGER, PARAMETER :: F_BC_Tsurf=2 ! surface temperature boundary condition.
      INTEGER, PARAMETER :: F_BC_Msurf=27 ! boundary condition for stellar mass.
      INTEGER, PARAMETER :: F_PQ=3 ! pressure term for mesh spacing function Q.
      INTEGER, PARAMETER :: F_PQ_dk=4 ! dPQ/dk.
      INTEGER, PARAMETER :: F_RQ=5 ! radius term for mesh spacing function Q.
      INTEGER, PARAMETER :: F_RQ_dk=6 ! dRQ/dk.
      INTEGER, PARAMETER :: F_TQ=7 ! temperature term for mesh spacing function Q.
      INTEGER, PARAMETER :: F_TQ_dk=8 ! dTQ/dk.
      INTEGER, PARAMETER :: F_LQ=9 ! luminosity term for mesh spacing function Q.
      INTEGER, PARAMETER :: F_LQ1_dk=10 ! the part of dLQ/dk that is not multiplied by dM/dt.
      INTEGER, PARAMETER :: F_LQ2_dk=11 ! the part of dLQ/dk that is multiplied by dM/dt.
      INTEGER, PARAMETER :: F_M_dt=12 ! dM/dt.  change in mass coordinate of this meshpoint since the last timestep
      INTEGER, PARAMETER :: F_MQ=13 ! mass term for mesh spacing function Q.
      INTEGER, PARAMETER :: F_MQ_dk=14 ! dMQ/dk.
      INTEGER, PARAMETER :: F_SG=15 ! diffusion coefficient for convective mixing.
      INTEGER, PARAMETER :: F_WT=16 ! weighting term for evaluation of "average" values from adjacent meshpoint values.
      INTEGER, PARAMETER :: F_Q_dk=28 ! dQ/dk. meshpoint spacing adjusted automatically to keep this constant across the shells.
      INTEGER, PARAMETER :: NUMF=28    ! number of functions of the independent variables passed to/from the SOLVER
      INTEGER, PARAMETER :: N_FNs = NUMF ! The number of functions that the equations use directly.

      DOUBLE PRECISION :: FN1(NUMF) ! holds the function values for passing to the solver

      INTEGER, PARAMETER :: FN_VALUE = NUMV+1
      DOUBLE PRECISION :: FN(FN_VALUE, max_N_SHELLs, N_FNs)
      ! FN(J,K,F) if J=FN_VALUE, then FN(J,K,F) holds value of function F for shell K.
      ! for J between 1 and NUMV, FN(J,K,F) holds value of F for shell K with modified variable J.
      ! this is the function value used in computing the partial difference equations. 
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! These are values that are used in calculating the functions listed above, but are not directly used by the equations.
      
      DOUBLE PRECISION ::  FK        ! opacity.
      DOUBLE PRECISION ::  XHI           ! thermal conductivity = 4 a c T^3/(3 kappa rho^2 c_P). (cm^2/sec).
      DOUBLE PRECISION :: MK           ! dM/dk
      DOUBLE PRECISION :: MESH_Q          ! value of mesh spacing function.
      DOUBLE PRECISION :: QM           ! d(MESH_Q)/dM.
      DOUBLE PRECISION :: WL           ! convection velocity times mixing length.
      DOUBLE PRECISION :: WCV         ! convection velocity.
      DOUBLE PRECISION :: HP           ! pressure scale height for mixing length calculation.
      DOUBLE PRECISION :: XN           ! fractional abundance by mass of N14 -- computed from the other abundances
      DOUBLE PRECISION :: XMG     ! the abundance of Mg
      DOUBLE PRECISION :: XSI     ! the abundance of Si
      DOUBLE PRECISION :: XFE     ! the abundance of Fe
      DOUBLE PRECISION ::  EX          ! nuclear reaction energy release in ergs/g/sec.
      DOUBLE PRECISION ::  ENX           ! nuclear reaction neutrino loss in ergs/g/sec.
      DOUBLE PRECISION ::  EN          ! neutrino losses from sources other than the nuclear reaction network.
      DOUBLE PRECISION ::  eps_plasma, eps_brem, eps_pair, eps_photo
      DOUBLE PRECISION ::  ENG          ! energy from gravitational contraction in ergs/g/sec.
      DOUBLE PRECISION ::  T            ! temperature.
      DOUBLE PRECISION :: R             ! radius
      DOUBLE PRECISION :: EG           ! superadiabaticity for convection.
      DOUBLE PRECISION :: GRAD ! actual temperature gradient.
      DOUBLE PRECISION :: GRADR ! radiative temperature gradient.

      ! These are indices in the SV array of saved values which serves as a cache to accelerate the creation of the Jacobian.
     
      INTEGER, PARAMETER :: SV_R=49 ! radius
      INTEGER, PARAMETER :: SV_T=7 ! temperature.
      INTEGER, PARAMETER :: SV_MK=45 ! dM/dk
      INTEGER, PARAMETER :: SV_MESH_Q=50 ! value of mesh spacing function.
      INTEGER, PARAMETER :: SV_QM=51 ! d(MESH_Q)/dM.
      ! density
      INTEGER, PARAMETER :: SV_RHO=5 ! density (g/cm^3)
      INTEGER, PARAMETER :: SV_LNRHO=2 ! ln(density).
      ! pressure
      INTEGER, PARAMETER :: SV_P=4 ! pressure (ergs/cm^3).
      INTEGER, PARAMETER :: SV_LNP=1 ! ln(pressure).
      ! internal energy
      INTEGER, PARAMETER :: SV_U=3 ! internal energy (ergs per gram).
      INTEGER, PARAMETER :: SV_SCP=12 ! specific heat capacity at constant pressure.
      ! entropy
      INTEGER, PARAMETER :: SV_S=15 ! entropy per gram.
      INTEGER, PARAMETER :: SV_SF=8 ! d_S / d_LNF (for use in calculating dS/dt)
      INTEGER, PARAMETER :: SV_ST=9 ! d_S / d_LNT (for use in calculating dS/dt)
      ! opacity
      INTEGER, PARAMETER :: SV_FK=6 ! opacity.
      INTEGER, PARAMETER :: SV_XHI=14 ! thermal conductivity = 4 a c T^3/(3 kappa rho^2 c_P). (cm^2/sec).
      ! convection related
      INTEGER, PARAMETER :: SV_GRADA=11 ! adiabatic temperature gradient, dln(T)/dln(P).
      INTEGER, PARAMETER :: SV_GRAD=48 ! actual temperature gradient, dln(T)/dln(P).
      INTEGER, PARAMETER :: SV_GRADR=46 ! temp gradient required to carry flux by radiative diffusion alone.
      !INTEGER, PARAMETER :: SV_EG=47 ! superadiabaticity for convection.
      INTEGER, PARAMETER :: SV_WCV=53 ! convection velocity.
      INTEGER, PARAMETER :: SV_WL=52 ! convection velocity times mixing length.
      INTEGER, PARAMETER :: SV_HP=54 ! pressure scale height for mixing length calculation.
      ! breakdown of components of pressure
      INTEGER, PARAMETER :: SV_PR=16 ! radiation pressure.
      INTEGER, PARAMETER :: SV_PG=17 ! gas pressure.
      INTEGER, PARAMETER :: SV_P_ION=72 ! ideal gas ion pressure
      INTEGER, PARAMETER :: SV_PCORR=73 ! correction to pressure for non-ideal gas
      ! breakdown of components of energy
      INTEGER, PARAMETER :: SV_EN=20 ! neutrino losses from sources other than the nuclear reaction network.
      INTEGER, PARAMETER :: SV_EX=41 ! nuclear reaction energy release in ergs/g/sec.
      INTEGER, PARAMETER :: SV_ENX=42 ! nuclear reaction neutrino loss in ergs/g/sec.
      INTEGER, PARAMETER :: SV_ENG=44 ! energy from gravitational contraction in ergs/g/sec.
      ! nuclear reaction rates
      INTEGER, PARAMETER :: SV_RPP=21 ! equilibrium pp chain        2 p -> 1/2 He4
      INTEGER, PARAMETER :: SV_R33=22 ! He3 + He3 (ppI chain)      He3 (He3,2p) He4
      INTEGER, PARAMETER :: SV_R34=23 ! He3 + He4 (ppII chain)  He3 (He4,gamma) Be7
      INTEGER, PARAMETER :: SV_RBE=24 ! Be7 + H1 (ppIII chain)  Be7 (p,gamma) B8 (positron,neutrino) Be8* (alpha) He4
      INTEGER, PARAMETER :: SV_RBP=25 ! Be7 decay                 Be7 (electron,neutrino) Li7 (p,alpha) He4
      INTEGER, PARAMETER :: SV_RPC=26 ! C12 + H1 (CN cycle)    C12 (p,positron,neutrino) C13 (p,gamma) N14
      INTEGER, PARAMETER :: SV_RPO=28 ! O16 + H1 (ON cycle)    O16 (p,positron,neutrino) O17 (p,alpha) N14
      INTEGER, PARAMETER :: SV_R3A=29 ! triple alpha process      He4 (alpha) Be8* (alpha,gamma) C12
      INTEGER, PARAMETER :: SV_RAC=30 ! He4 + C12                 C12 (alpha,gamma) O16
      INTEGER, PARAMETER :: SV_RAN=31 ! He4 + N14                 N14 (alpha,gamma) F18 (1/2 alpha,gamma) Ne20
      INTEGER, PARAMETER :: SV_RAO=32 ! He4 + O16                 O16 (alpha,gamma) Ne20
      INTEGER, PARAMETER :: SV_RANE=33 ! He4 + Ne20               Ne20 (alpha,gamma) Mg24
      INTEGER, PARAMETER :: SV_RCO=35 ! C12 + O16                 C12 (O16,alpha,gamma) Mg24
      INTEGER, PARAMETER :: SV_ROO=36 ! 2 O16                        O16 (O16,alpha,gamma) Si28 (gamma,alpha) Mg24
      INTEGER, PARAMETER :: SV_RGNE=37 ! Ne20 decay               Ne20 (gamma,alpha) O16
      INTEGER, PARAMETER :: SV_RGMG=38 ! Mg24 decay               Mg24 (gamma,alpha) Ne20
      INTEGER, PARAMETER :: SV_RCCG=39 ! less probable CC         C12 (C12,gamma) Mg24
      INTEGER, PARAMETER :: SV_RCC=34 ! more probable CC       C12 (C12,alphagamma) Ne20
      INTEGER, PARAMETER :: SV_RPNG=40 ! more probable N14+p      N14 (p,positronneutrino) N15 (p,gamma) O16
      INTEGER, PARAMETER :: SV_RPN=27 ! less probable N14+p    N14 (p,positronneutrino) N15 (p,alpha) C12 
      ! information about densities of electrons and ions.
      INTEGER, PARAMETER :: SV_NE=67 ! moles of free electrons per gram of star. so 1/NE is mean molecular weight per free electron.
      INTEGER, PARAMETER :: SV_NE1=68 ! what NE would be in case of complete ionization.
      INTEGER, PARAMETER :: SV_NI=69 ! moles of molecules (including monatomic) per gram.
      INTEGER, PARAMETER :: SV_NZZ=70 ! moles of Z^2 per gram of star where Z is nuclear charge.
      INTEGER, PARAMETER :: SV_AVM=71 ! moles of baryons per gram of star.
      ! information about abundances of elements.
      INTEGER, PARAMETER :: SV_XN=74 ! fractional abundance by mass of N14 -- computed from the other abundances
      INTEGER, PARAMETER :: SV_XMG=55 ! the abundance of Mg
      INTEGER, PARAMETER :: SV_XSI=56 ! the abundance of Si
      INTEGER, PARAMETER :: SV_XFE=57 ! the abundance of Fe
      INTEGER, PARAMETER :: SV_N1=58 ! Ns hold relative abundances by number (number per gram of stellar material)
      INTEGER, PARAMETER :: SV_N4=59
      INTEGER, PARAMETER :: SV_N12=60
      INTEGER, PARAMETER :: SV_N14=61
      INTEGER, PARAMETER :: SV_O16=62
      INTEGER, PARAMETER :: SV_N20=63
      INTEGER, PARAMETER :: SV_N24=64
      INTEGER, PARAMETER :: SV_N28=65
      INTEGER, PARAMETER :: SV_N56=66
      ! max index
      INTEGER, PARAMETER :: NSVARS=74  ! number of saved variables

      DOUBLE PRECISION :: SV(NSVARS) ! saved values

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! coefficients used in the mesh spacing function
      DOUBLE PRECISION, PARAMETER :: CT(10) = (/ 0D0, 0D0, 0.05D0, 5D-2, 0.15D0, 2D-2, 0.45D0, 1D-4, 1D15, 2D4 /)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DOUBLE PRECISION, PARAMETER :: CC=.176D0, CN=.052D0, CO=.502D0, CNE=.092D0, CMG=.034D0, CSI=.072D0, CFE=.072D0
      DOUBLE PRECISION :: CZS, CZS_SAV ! metals mass fraction.
      DOUBLE PRECISION :: CH, CH_SAV ! H abundance by mass (0.76 - 3*Z)  Metals abundance Z is fixed by choice of ZAMS data directory.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      CONTAINS
      
      SUBROUTINE SAV_shell_data( IO_UNIT, read )
      INTEGER, INTENT(IN) :: IO_UNIT
      LOGICAL, INTENT(IN) :: read
      IF ( read ) THEN
         READ (IO_UNIT) CH, CZS
      ELSE
         WRITE (IO_UNIT) CH, CZS
      END IF
      END SUBROUTINE SAV_shell_data

      SUBROUTINE SAV_shell_data_internal(read_flag)
      LOGICAL, INTENT(IN) :: read_flag
      IF ( read_flag ) THEN
         CH=CH_SAV; CZS=CZS_SAV
      ELSE
		 CH_SAV=CH; CZS_SAV=CZS
	  END IF
      END SUBROUTINE SAV_shell_data_internal

      END MODULE ez_shell_data

