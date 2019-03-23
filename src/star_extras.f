      MODULE star_extras ! This is part of the data interface to the stellar evolution code 
      ! The items in here cost a little more to get than those in star_data.
      USE star_data
      USE ez_magnitude_data
      IMPLICIT NONE
      
      LOGICAL :: have_Extras ! this is set false before each call on your Check_Model
      ! it is checked by the EZ_Extras routine
      ! if this flag is false when EZ_Extras is called, the data is evaluated and the flag is set true
      ! if this flag is true when EZ_Extras is called, it simply returns without doing anything.
      ! this scheme means you don't need to worry about keeping track to avoid extra calls on EZ_Extras.
      
      ! The following variables are only computed when the EZ_Extras routine is called.
      ! NB: these are all outputs.  The system doesn't look at these values again after passing them to your Check_Model.
      
      DOUBLE PRECISION :: center_GAM ! plasma interaction parameter at center of star

      DOUBLE PRECISION :: power_H_burn ! total power from hydrogen consuming reactions (in Lsolar units)
      DOUBLE PRECISION :: power_He_burn ! total power from reactions burning helium (in Lsolar units)
      DOUBLE PRECISION :: power_Metal_burn ! total power from reactions burning metals -- but not consuming H1 or He4 (in Lsolar units)
      DOUBLE PRECISION :: power_Neutrinos ! total power loss from neutrinos (in Lsolar units)
      
      ! power_H_burn is the sum of the following
      DOUBLE PRECISION :: power_PP ! PP chain
      DOUBLE PRECISION :: power_CNO ! CNO cycle
      
      ! power_He_burn is the sum of the following
      DOUBLE PRECISION :: power_3_alpha ! triple alpha reaction
      DOUBLE PRECISION :: power_C_alpha ! C12 + He4 -> O16
      DOUBLE PRECISION :: power_N_alpha ! N14 + 3/2 He4 -> Ne20
      DOUBLE PRECISION :: power_O_alpha ! O16 + He4 -> Ne20
      DOUBLE PRECISION :: power_Ne_alpha ! Ne20 + He4 -> Mg24
      
      ! power_Metal_burn is the sum of the following
      DOUBLE PRECISION :: power_CC_Ne ! 2 C12 -> Ne20 + He4
      DOUBLE PRECISION :: power_CO ! C12 + O16 -> Mg24 + He4
      DOUBLE PRECISION :: power_OO ! 2 O16 -> Mg24 + 2 He4
      DOUBLE PRECISION :: power_Ne_decay ! Ne20 -> O16 + He4
      DOUBLE PRECISION :: power_Mg_decay ! Mg24 -> Ne20 + He4
      DOUBLE PRECISION :: power_CC_Mg ! 2 C12 -> Mg24
      
      ! power_Neutrinos is the sum of the following
      DOUBLE PRECISION :: power_plasmon_neutrinos ! plasmon neutrinos
      DOUBLE PRECISION :: power_brem_neutrinos ! bremsstrahlung neutrinos
      DOUBLE PRECISION :: power_pair_neutrinos ! pair annihilation neutrinos
      DOUBLE PRECISION :: power_photo_neutrinos ! photo neutrinos
     
      INTEGER, PARAMETER :: BSZ=8
      ! dimension of following arrays of convection boundaries
      DOUBLE PRECISION :: conv_turnover_Time(BSZ)
      ! convective turnover times for zones; 0 indicates a nonconvective region
      ! current number of valid entries is given by conv_boundary_count (below)
      DOUBLE PRECISION :: conv_boundaries_R(BSZ) ! location of top of zone in Rsun units
      ! radiative/convective boundaries (Rsolar)
      DOUBLE PRECISION :: conv_boundaries_M(BSZ) ! location of top of zone in Msun units
      ! radiative/convective boundaries (Msolar)
      INTEGER :: conv_boundary_count
      ! the actual number of radiative/convective boundaries in arrays conv_boundaries_R, conv_boundaries_M, and conv_turnover_Time

      ! estimated color magnitudes (from Lejeune, Cuisinier, Buser (1998) A&AS 130, 65)
      DOUBLE PRECISION :: color_magnitudes(n_mags)
      ! see ez_magnitudes_data.f for more information about what is available
      
      ! "core mass" as used for mesh function
      DOUBLE PRECISION :: mesh_core_mass ! in Msolar units

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
   ! the array SX holds a variety of variables for all the meshpoints of the model.
   ! it's 2 dimensional  1st is the variable index, 2nd is the meshpoint index.
      INTEGER :: SX_CNTR ! index in SX of the innermost meshpoint. values depends on how many shells there are in current model.
      INTEGER, PARAMETER :: SX_SURF = N_SURF_shell ! index in SX of the outermost meshpoint
      ! By convention, SX_SURF < SX_CNTR. i.e., shells are numbered from surface to center.
     
   ! the following defines symbolic names for the variable indices in array SX
   ! for example, SX(SX_RHO,SX_CNTR) holds the central density
   ! your code should always refer to these by symbolic name since future versions of the system will certainly change the numbers
      
   ! independent variables for stellar structure
      INTEGER, PARAMETER :: SX_R=1
      ! radius (Rsolar)
      INTEGER, PARAMETER :: SX_T=2
      ! temperature (K)
      INTEGER, PARAMETER :: SX_M=3
      ! mass interior to this meshpoint (Msolar)
      ! this is automatically adjusted during the run to best position a fixed number of meshpoints
      INTEGER, PARAMETER :: SX_L=4
      ! luminosity (Lsolar) -- net rate of outward energy flow at this meshpoint (nuclear plus thermal minus neutrino)
      ! this can be negative as a result of neutrino cooling at the center
      INTEGER, PARAMETER :: SX_PSI=5
      ! the degeneracy parameter (psi is the same as defined in Kippenhahn & Weigert)
      ! psi replaces density or pressure as an independent variable

   ! composition variables -- all are local mass fractions
      INTEGER, PARAMETER :: SX_XH=6
      ! H1
      INTEGER, PARAMETER :: SX_XHE=7
      ! He4
      INTEGER, PARAMETER :: SX_XC=8
      ! C12
      INTEGER, PARAMETER :: SX_XN=9
      ! N14
      INTEGER, PARAMETER :: SX_XO=10
      ! O16
      INTEGER, PARAMETER :: SX_XNE=11
      ! Ne20
      INTEGER, PARAMETER :: SX_XMG=12
      ! Mg24

   ! some of the major dependent variables
      INTEGER, PARAMETER :: SX_P=13
      ! pressure (dynes/cm^2)
      INTEGER, PARAMETER :: SX_RHO=14
      ! density (gm/cm^3)
      INTEGER, PARAMETER :: SX_OPACITY=15
      ! opacity
      INTEGER, PARAMETER :: SX_U=16
      ! specific energy (ergs/g)
      INTEGER, PARAMETER :: SX_S=17
      ! specific entropy (ergs/K/g)
      INTEGER, PARAMETER :: SX_MU=18
      ! grams per mole of gas particles (including electrons) -- mean particle molecular weight
      INTEGER, PARAMETER :: SX_NE=19
      ! free electron abundance (moles per gram of star) -- 1/NE is mean molecular weight per free electron
      INTEGER, PARAMETER :: SX_NE_TOT=20
      ! total electron abundance (moles per gram of star)
      INTEGER, PARAMETER :: SX_SCP=21
      ! specific heat capacity at constant pressure (ergs per gram per K)
      INTEGER, PARAMETER :: SX_WF=22
      ! weight factor for converting ergs/g/sec to total power for this shell (in Lsolar units)
      ! incorporates delta_M plus conversion to solar units

   ! the pressure is sum of the following terms
      INTEGER, PARAMETER :: SX_PRAD=23
      ! radiation pressure
      INTEGER, PARAMETER :: SX_PEL=24
      ! electron gas pressure
      INTEGER, PARAMETER :: SX_PION=25
      ! ion gas pressure
      INTEGER, PARAMETER :: SX_PCORR=26
      ! correction term for non-ideal gas effects
 
   ! variables related to temperature gradients and convection
      INTEGER, PARAMETER :: SX_GRAD_AD=27
      ! adiabatic temperature gradient
      INTEGER, PARAMETER :: SX_GRAD_RAD=28
      ! radiative temperature gradient required to carry entire luminosity without convection
      INTEGER, PARAMETER :: SX_GRAD_STAR=29
      ! actual temperature gradient
      INTEGER, PARAMETER :: SX_SG=30
      ! diffusion coefficient for convective mixing (if is 0, then no mixing)
      INTEGER, PARAMETER :: SX_CV=31
      ! convection velocity -- zero if not in a convection zone
      INTEGER, PARAMETER :: SX_WL=32
      ! convection velocity times mixing length

   ! nuclear power currently being generated everywhere inside this meshpoint (in Lsolar units)
      INTEGER, PARAMETER :: SX_L_H=33
      ! from hydrogen burning
      INTEGER, PARAMETER :: SX_L_HE=34
      ! from helium burning
      INTEGER, PARAMETER :: SX_L_Z=35
      ! from metals burning
      INTEGER, PARAMETER :: SX_L_NEU=36
      ! lost to neutrinos
     
   ! categories of local power sources at this shell (ergs/g/sec)
      INTEGER, PARAMETER :: SX_EPS_NUC=37
      ! nuclear energy generation rate (photons and neutrinos)
      INTEGER, PARAMETER :: SX_EPS_NEU=38
      ! neutrino energy loss rate from sources other than nuclear reactions (see SX_EPS_NUC_NEU from losses from reactions)

      INTEGER, PARAMETER :: SX_DQ_DK=39
      ! change in mesh spacing function Q at this shell -- i.e., approximation for dQ/dk where k is the shell index.
      ! this should be the same for all the shells when the spacing is optimal

   ! variables related to the primary hydrogen burning reactions, PP chain and CNO cycle (ergs/g/sec)
      INTEGER, PARAMETER :: SX_EPS_H=40
      ! total energy generation rate from hydrogen burning by PP and CNO at this location 
      INTEGER, PARAMETER :: SX_EPS_PP=41
      ! PP chain energy generation rate
      INTEGER, PARAMETER :: SX_EPS_CNO=42
      ! CNO cycle energy generation rate

   ! variables related to helium burning reactions (ergs/g/sec)
      INTEGER, PARAMETER :: SX_EPS_HE=43
      ! total energy generation rate from helium burning at this location
      
   ! helium burning is the sum of the following reactions (ergs/g/sec)
      INTEGER, PARAMETER :: SX_EPS_3A=44
      ! triple alpha
      INTEGER, PARAMETER :: SX_EPS_AC=45
      ! C12 + He4 -> O16
      INTEGER, PARAMETER :: SX_EPS_AN=46
      ! N14 + 3/2 He4 -> Ne20
      INTEGER, PARAMETER :: SX_EPS_AO=47
      ! O16 + He4 -> Ne20
      INTEGER, PARAMETER :: SX_EPS_ANE=48
      ! Ne20 + He4 -> Mg24
          
   ! variables related to other nuclear reactions -- metal burning that doesn't consume helium (ergs/g/sec)
      INTEGER, PARAMETER :: SX_EPS_Z=49
      ! total energy generation rate from metal burning at this location

   ! metal burning is the sum of the following reactions (ergs/g/sec)
      INTEGER, PARAMETER :: SX_EPS_CCA=50
      ! 2 C12 -> Ne20 + He4
      INTEGER, PARAMETER :: SX_EPS_CO=51
      ! C12 + O16 -> Mg24 + He4
      INTEGER, PARAMETER :: SX_EPS_OO=52
      ! 2 O16 -> Mg24 + 2 He4
      INTEGER, PARAMETER :: SX_EPS_GNE=53
      ! Ne20 -> O16 + He4
      INTEGER, PARAMETER :: SX_EPS_CCG=54
      ! 2 C12 -> Mg24

   ! photodisintegration
      INTEGER, PARAMETER :: SX_EPS_GMG=55
      ! Mg24 -> Ne20 + He4

   ! breakdown of neutrino losses (ergs/g/sec) -- only computed if ITOH_VINTAGE is at least 1996
      INTEGER, PARAMETER :: SX_NEU_plasma=56
      ! from plasmon neutrinos 
      INTEGER, PARAMETER :: SX_NEU_brem=57
      ! from bremsstrahlung
      INTEGER, PARAMETER :: SX_NEU_pair=58
      ! from pair annihilation
      INTEGER, PARAMETER :: SX_NEU_photo=59
      ! from photo neutrinos

   ! various logarithmic derivatives
      INTEGER, PARAMETER :: SX_HI1=60
      ! actual dln(RHO)/dln(P)
      INTEGER, PARAMETER :: SX_HI2=61
      ! actual -dln(R)/dln(P)
      INTEGER, PARAMETER :: SX_HI3=62
      ! actual -dln(M)/dln(P)
      INTEGER, PARAMETER :: SX_GAMMA1=63
      ! adiabatic exponent, dln(P)/dln(RHO) at constant entropy.
      
   ! mass fractions for various ionization states of H, He, C, N, and O
      INTEGER, PARAMETER :: SX_XH2=64 ! mass fraction of molecular hydrogen
      INTEGER, PARAMETER :: SX_XH_plus=65 ! neutral H mass fraction = SX_XH - SX_XH_plus - SX_XH2
      INTEGER, PARAMETER :: SX_XHE_plus1=66 ! neutral HE mass fraction = SX_XHE - SX_XHE_plus1 - SX_XHE_plus2
      INTEGER, PARAMETER :: SX_XHE_plus2=67
      
   ! misc
      INTEGER, PARAMETER :: SX_GAM=68 ! the plasma interaction parameter
      INTEGER, PARAMETER :: SX_EPS_G=69 ! thermal energy generation rate (ergs/g/sec)
      INTEGER, PARAMETER :: SX_EPS_NUC_NEU=70
      ! neutrino energy loss rate from nuclear reactions (see SX_EPS_NEU for losses from other sources)

   ! mesh function inforamtion
      INTEGER, PARAMETER :: SX_PQ=71 ! the pressure term
      INTEGER, PARAMETER :: SX_TQ=72 ! the temperature term
      INTEGER, PARAMETER :: SX_RQ=73 ! the radius term
      INTEGER, PARAMETER :: SX_MQ=74 ! the mass term
      INTEGER, PARAMETER :: SX_DM_DK=75 ! Msun per shell as derived from mesh function
     
      INTEGER, PARAMETER :: SX_NUM=75
      ! total number of entries per meshpoint
      
      DOUBLE PRECISION :: SX(SX_NUM,max_N_SHELLs)
      ! this is the array that holds all of the values
      ! for example, SX(SX_RHO,N_CNTR_shell) holds the center mesh point density.

      END MODULE star_extras

