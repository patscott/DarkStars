      MODULE star_data
      ! This is the primary data interface to the stellar evolution code
      ! (also see star_extra, star_constants, and star_controls).
      IMPLICIT NONE
      
      ! The items in this module are automatically provided for your Check_Model routine.

      ! There's lots more information of potential interest that's not listed in this module,
      !      such as densities, opacities, pressures, nuclear reaction rates, convection zones, and on and on ...
      ! All of these require evaluating the equation of state which isn't a trivial operation.
      ! Since you may not need any of that, the system doesn't automatically provide it the way it does the following.
      ! But if you want it, by all means go get it --
      !    you just need to call EZ_Extra to have it fill in the info in the star_extras module.
      
      ! In addition, the module star_constants holds some useful constants, and
      ! the module star_controls holds the control parameter interface to the stellar evolution code.
      
      ! NOTE: these are all strictly outputs from the system.
      ! The system doesn't look at them after producing them for your Check_Model.
      ! These aren't controls.  For example, you cannot change the age of the star by setting the star_Age value.
      ! The control parameters are found in the star_controls module.
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      INTEGER :: model_Number ! normally this just increments by one each timestep,
      !   but it can get smaller if the system is forced to backup.
      
      ! age
      DOUBLE PRECISION :: star_Age ! age (in years)
      DOUBLE PRECISION :: time_Step ! timestep since previous model (in years)
      
      ! surface (at the outermost mesh point)
      DOUBLE PRECISION :: log_Luminosity ! log10(total luminosity in Lsolar units)
      DOUBLE PRECISION :: log_Radius ! log10(total radius in Rsolar units)
      DOUBLE PRECISION :: log_surface_Temp ! log10(temperature at surface)
      
      ! center (actually, a little out from the center at the centermost mesh point)
      DOUBLE PRECISION :: log_center_Temp ! log10(center temperature in Kelvin)
      DOUBLE PRECISION :: log_center_Density ! log10(center density in g/cm^3)
      DOUBLE PRECISION :: log_center_Pressure ! log10(center pressure in dynes/cm^2)
      DOUBLE PRECISION :: center_Degeneracy ! PSI at the centermost meshpoint.  (PSI is the electron chemical potential in units of kT)
      DOUBLE PRECISION :: center_H  ! central fractional abundance by mass for H1
      DOUBLE PRECISION :: center_He ! for He4
      DOUBLE PRECISION :: center_C  ! for C12
      DOUBLE PRECISION :: center_N  ! for N14
      DOUBLE PRECISION :: center_O  ! for O16
      DOUBLE PRECISION :: center_Ne ! for Ne20
      DOUBLE PRECISION :: center_Mg ! for Mg24
      DOUBLE PRECISION :: center_Si ! for Si28
      DOUBLE PRECISION :: center_Fe ! for Fe56
     
      ! total mass
      DOUBLE PRECISION :: initial_Mass ! initial total stellar mass (in Msolar units)
      DOUBLE PRECISION :: star_Mass ! total stellar mass (in Msolar units)
      DOUBLE PRECISION :: star_Mdot ! d(star_Mass)/dt (in Msolar per year)
      DOUBLE PRECISION :: initial_Z ! initial metallicity
      DOUBLE PRECISION :: star_Mass_H  ! total H1
      DOUBLE PRECISION :: star_Mass_He ! total He4
      DOUBLE PRECISION :: star_Mass_C  ! total C12
      DOUBLE PRECISION :: star_Mass_N  ! total N14
      DOUBLE PRECISION :: star_Mass_O  ! total O16
      DOUBLE PRECISION :: star_Mass_Ne ! total Ne20
      DOUBLE PRECISION :: star_Mass_Mg ! total Mg24
      DOUBLE PRECISION :: star_Mass_Si ! total Si28
      DOUBLE PRECISION :: star_Mass_Fe ! total Fe56
     
      ! helium, carbon, and oxygen cores -- as defined by the abundance thresholds.
      DOUBLE PRECISION :: mass_He_Core, radius_He_Core
      ! mass_He_Core is mass location where H1 abundance XH reaches the parameter value core_param_He as go out from center (in Msolar units)
      ! gives mass of helium core (including possible inner carbon/oxygen cores and any other elements present as well)
      ! radius_He_Core is radial location (in Rsolar units) corresponding to mass_He_Core
      DOUBLE PRECISION :: mass_C_Core, radius_C_Core
      ! mass_C_Core is mass location where He4 abundance XHE reaches the parameter value core_param_C as go out from center (Msolar)
      ! gives mass of carbon core (including possible inner oxygen core and any other elements present as well)
      ! radius_C_Core is radial location (in Rsolar units) corresponding to mass_C_Core
      DOUBLE PRECISION :: mass_O_Core, radius_O_Core
      ! mass_O_Core is mass location where C12 abundance XC reaches the parameter value core_param_O as go out from center (Msolar)
      ! gives mass of oxygen core (including any other elements present as well)
      ! radius_O_Core is radial location (in Rsolar units) corresponding to mass_O_Core
     
      ! timescales
      DOUBLE PRECISION :: dynamic_Timescale ! dynamic timescale (seconds) -- estimated by SQRT(R^3/(G*M))
      DOUBLE PRECISION :: KH_Timescale ! Kelvin-Helmholtz timescale (years) -- proportional to graviational energy divided by luminosity
      DOUBLE PRECISION :: nuc_Timescale ! nuclear timescale (years) -- proportional to mass divided by luminosity

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! The star is divided into concentric "zones" with a "meshpoint" between each pair of adjacent zones.
      ! The variables such as composition, temperature, mass, and radius are defined at meshpoints.
      ! In the current configuration of the system, there are 200 zones and therefore there are 199 meshpoints.
      ! Zone 1 is at the surface and zone 200 is at the center.
      ! Meshpoint 1 is at the bottom of the surface zone and meshpoint 199 is at the top of the center zone.
      ! There is no meshpoint at the true center where R and M are zero,
      ! and there is no meshpoint at the "true surface" (whatever that might be!).
      ! The meshpoints are not assigned fixed locations in either radius or mass.
      ! Instead, the system automatically adjusts the locations so as to have a higher density of meshpoints in areas
      ! where there is lots of "activity" such as nuclear burning, and fewer zones in areas that are less active.
      !
      ! The precise locations of meshpoints are choosen as part of the solution to the stellar equations.
      ! This is accomplished by additional equations that ensure that the meshpoints are equally spaced
      ! with respect to a "mesh function", Q.  Using 'k' as the meshpoint variable, the equations state that dQ/dk
      ! is constant.  The mesh function Q varies with mass and changes sharply in areas where there is lots of
      ! activity.  A sharp change in Q forces a high density of meshpoints in order to keep dQ/dk constant.
      ! The result is the desired dynamic adjustment in mesh placement as the star evolves.
      !
      ! Note that we often use 'shell' as a synonym for 'meshpoint' in the following documentation.
      ! Some of the variables below concern overall properties of the star and some are on a shell-by-shell basis.
      ! The 'global' properties, such as age, are listed separately in the following declarations.
      ! The 'per-shell' properties are gathered together in arrays indexed by meshpoint number and variable number.
      ! 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Meshpoints are also known as 'shells' since they correspond to spherical shells between zones in the star.
      INTEGER, PARAMETER :: max_N_SHELLs = 2000 ! The max number of mesh points.      
      INTEGER, PARAMETER :: N_SURF_shell = 1 ! This is the index for the meshpoint nearest the surface.
      INTEGER :: N_SHELLs     ! The current number of mesh points.
      INTEGER :: N_CNTR_shell ! This is the current index for the meshpoint nearest the center. (same as N_SHELLs)
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      ! These are the identifiers for the independent variables at each meshpoint.  The actual values are in the MESH_Vs array declared below.
      ! composition variables
      INTEGER, PARAMETER :: V_XH=5 ! fractional abundance by mass of H1.
      INTEGER, PARAMETER :: V_XHE=9 ! fractional abundance by mass of He4.
      INTEGER, PARAMETER :: V_XC=10 ! fractional abundance by mass of C12.
      INTEGER, PARAMETER :: V_XO=3 ! fractional abundance by mass of O16.
      INTEGER, PARAMETER :: V_XNE=11 ! fractional abundance by mass of Ne20.
      ! structure variables (see star_constants for CMSN, CRSN, CLSN, and CLN for conversion of units)
      INTEGER, PARAMETER :: V_M=4 ! mass (in units of 10^33 g) interior to the mesh point.  divide M by CMSN to get mass in Msolar units.
      INTEGER, PARAMETER :: V_LNR=7 ! ln(radius R in units of 10^11 cm).  divide exp(LNR) by CRSN to get radius in Rsolar units.
      INTEGER, PARAMETER :: V_L=8 ! luminosity outward in units of 10^33 ergs per second.  divide L by CLSN to get luminosity in Lsolar units.
      INTEGER, PARAMETER :: V_LNT=2 ! ln(temperature T) at the mesh point.  divide LNT by CLN to convert to log10 of temperature.
      INTEGER, PARAMETER :: V_LNF=1 ! ln(F) at the mesh point. the degeneracy parameter PSI can be obtained by calling PSI_from_LNF
      ! mesh spacing variable
      INTEGER, PARAMETER :: V_Q_dk=6 ! derivative wrt meshpoint k of meshpoint spacing function Q.  should be same for all meshpoints.
      
      INTEGER, PARAMETER :: V_MAX=11, NUMV=11 ! number of variables per shell
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      ! These arrays hold the information about the current model and the most recent increments from the previous model.
      DOUBLE PRECISION :: MESH_Vs(V_MAX, max_N_SHELLs) ! the values of the independent variables for each meshpoint.
      ! For example, MESH_Vs(V_LNT, N_CNTR_SHELL) holds ln(T) for the centermost meshpoint.
      DOUBLE PRECISION :: MESH_DVs(V_MAX, max_N_SHELLs) ! the increments in the variables for the last timestep.
      ! For example, MESH_DVs(V_LNT, N_CNTR_SHELL) holds the change in ln(T) for the centermost meshpoint since the previous model.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      CONTAINS
      
      DOUBLE PRECISION FUNCTION PSI_from_LNF(LNF) 
      DOUBLE PRECISION, INTENT(IN) :: LNF
      DOUBLE PRECISION :: W
      W = DSQRT(1D0 + DEXP(LNF))
      PSI_from_LNF = 2D0*(W - DLOG(W + 1D0)) + LNF
      END FUNCTION PSI_from_LNF
     
      END MODULE star_data

