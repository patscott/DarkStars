      MODULE ez_report
      USE ez_shell, ONLY : FUNCS1
      USE ez_utils
      USE ez_shell_data
      USE ez_nuclear_data
      USE ez_state_data
      USE ez_solve_data
      USE ez_cycle_data
      USE ez_ionization_data
      USE star_data
      USE star_extras
      USE ez_data
      IMPLICIT NONE 
      
      CONTAINS 


      SUBROUTINE Complete_Model ! get the stuff that can be easily calculated from the meshpoint variables
      INTEGER :: IK, IKK
      DOUBLE PRECISION :: W1, Radius, DM, DR, HPC, cur_M, prev_M, cur_R, prev_R
      DOUBLE PRECISION :: cur_H, prev_H, cur_HE, prev_HE, cur_C, prev_C, luminosity

      ! clear out the variables in star_extras
      have_Extras = .FALSE.
      power_H_burn = 0D0
      power_He_burn = 0D0
      power_Metal_burn = 0D0
      power_Neutrinos = 0D0
      power_PP = 0D0
      power_CNO = 0D0
      power_3_alpha = 0D0
      power_C_alpha = 0D0
      power_N_alpha = 0D0
      power_O_alpha = 0D0
      power_Ne_alpha = 0D0
      power_CC_Ne = 0D0
      power_CO = 0D0
      power_OO = 0D0
      power_Ne_decay = 0D0
      power_Mg_decay = 0D0
      power_CC_Mg = 0D0
      power_plasmon_neutrinos = 0D0
      power_brem_neutrinos = 0D0
      power_pair_neutrinos = 0D0
      power_photo_neutrinos = 0D0
      conv_turnover_Time = 0D0
      conv_boundaries_R = 0D0
      conv_boundaries_M = 0D0
      conv_boundary_count = -1
      SX_CNTR = N_CNTR_shell
      SX(1:SX_NUM, N_SURF_shell:N_CNTR_shell) = 0D0
      ! set the variables in star_data
      IF ( JNN.EQ.0 ) THEN
         model_Number = -JMOD
      ELSEIF ( JNN.GT.0 ) THEN
         model_Number = JMOD + 1
      ELSE
         model_Number = JMOD
      END IF
      MESH_Vs = H + DH
      MESH_DVs = DH
      Radius = EXP(MESH_Vs(V_LNR,N_SURF_shell)) ! in units of 10^11 cm
      log_Radius =  DLOG10(Radius/CRSN) ! log10(stellar radius in solar units)
      luminosity = MESH_Vs(V_L,N_SURF_shell) ! in units of 10^33 ergs per second
      log_Luminosity =  DLOG10(luminosity/CLSN) ! log10(stellar luminosity in solar units)
      log_center_Temp = MESH_Vs(V_LNT,N_CNTR_shell)/CLN ! log10(temperature at center)
      center_Degeneracy = PSI_from_LNF(MESH_Vs(V_LNF,N_CNTR_shell))
      log_surface_Temp = MESH_Vs(V_LNT,N_SURF_shell)/CLN ! log10(temperature at surface)
      DTY = DT/CSY
      time_Step = DTY         ! timestep in years
      star_Age = max(0D0, AGE + time_Step)
      IF ( JNN.EQ.0 ) AGE = star_Age
      SM = MESH_Vs(V_M,N_SURF_shell)/CMSN; star_Mass = SM               ! stellar mass in solar units
      star_Mdot = MESH_DVs(V_M,N_SURF_shell)/(CMSN*DTY)        ! dM/dt in Msolar per year
      KH_Timescale = 1.0D22*CG*star_Mass**2/(Radius*MESH_Vs(V_L,N_SURF_shell)*CSY) ! Kelvin-Helmholtz timescale in years
      nuc_Timescale = 4.0D10*star_Mass/luminosity ! nuclear timescale in years (e.g., about 10^10 years for sun)
      dynamic_Timescale = DSQRT(Radius**3/(CG*star_Mass)) ! dynamic timescale
      ! extra_Energy_time sets time scale for extra_Energy_param changes
      ! W1 large means can make large change in extra_Energy_param
      ! extra_Energy_max is saturation value for extra_Energy_param
      W1 = DT*extra_Energy_time/CSY
      ! DT is timestep in seconds; DT/CSY is timestep in years
      extra_Energy_param = extra_Energy_param*(1D0 + W1)*extra_Energy_max/(extra_Energy_max + W1*extra_Energy_param)
      star_Mass_H=0D0; star_Mass_He=0D0; star_Mass_C = 0D0; star_Mass_O = 0D0; star_Mass_Ne = 0D0
      mass_He_Core=0D0; mass_C_Core=0D0; mass_O_Core=0D0
      radius_He_Core=0D0; radius_C_Core=0D0; radius_O_Core=0D0
      cur_M = 0.d0; cur_R = 0.d0; cur_H = 0.d0; cur_He = 0.d0; cur_C = 0.d0 ! Added by Pat Scott 2008-01-29
      DO IK = N_CNTR_shell, N_SURF_shell, -1 ! go through the meshpoints from center to surface
         prev_M = cur_M; prev_R = cur_R
         cur_M = MESH_Vs(V_M,IK)/CMSN ! mass in solar units
         cur_R = EXP(MESH_Vs(V_LNR,IK))/CRSN ! radius in solar units
         SX(SX_M,IK) = cur_M
         prev_H = cur_H; cur_H = MESH_Vs(V_XH,IK)
         prev_HE = cur_HE; cur_HE = MESH_Vs(V_XHE,IK)
         prev_C = cur_C; cur_C = MESH_Vs(V_XC,IK)
         IF ( IK .EQ. N_CNTR_shell ) THEN ! call FUNCS1 for the center
            DVAR(1:NUMV) = DH(1:NUMV,IK)
            VAR(1:NUMV) = H(1:NUMV,IK) + DVAR(1:NUMV)
            CALL FUNCS1 ( IK, -1 )
            log_center_Density = LNRHO/CLN
            log_center_Pressure = LNP/CLN
            HPC = DSQRT(P/(CG*RHO*RHO))
            MC = 3.5D-33*RHO*HPC**3 ! central core mass used in computing mass term for equations
            DM = cur_M; DR = cur_R
         ELSE
            IKK = IK+1 ! the previous mesh point (only valid if IK .LT. N_CNTR_SHELL)
            DM = cur_M - prev_M; DR = cur_R - prev_R
            ! locate burning shell boundaries (XH = core_param_He, XHE = core_param_C, XC = core_param_O)
            ! CXB defines core boundary for printing in terms of H1 or He4 abundance by mass
            IF ( cur_H .GT. core_param_He .AND. prev_H .LT. core_param_He ) THEN ! estimate He4 core boundary
               ! mass_He_Core is location (in Msolar) where H1 abundance XH reaches core_param_He as go out from center
               mass_He_Core = FIND0(prev_M,prev_H-core_param_He,cur_M,cur_H-core_param_He)
               radius_He_Core = prev_R + (DR)*(mass_He_Core-prev_M)/DM
            END IF
            IF ( cur_HE .GT. core_param_C .AND. prev_HE .LT. core_param_C ) THEN
               ! mass_C_Core is location (in Msolar) where He4 abundance XHE reaches core_param_C as go out from center
               mass_C_Core = FIND0(prev_M,prev_HE-core_param_C,cur_M,cur_HE-core_param_C)
               radius_C_Core = prev_R + (DR)*(mass_C_Core-prev_M)/DM
            END IF
            IF ( cur_C .GT. core_param_O .AND. prev_C .LT. core_param_O ) THEN
               ! mass_O_Core is location (in Msolar) where C12 abundance XC reaches core_param_O as go out from center
               mass_O_Core = FIND0(prev_M,prev_C-core_param_O,cur_M,cur_C-core_param_O)
               radius_O_Core = prev_R + (DR)*(mass_O_Core-prev_M)/DM
            END IF
         END IF
         ! add up the component masses
         star_Mass_H = star_Mass_H + DM*cur_H ! star_Mass_H is total mass of H1
         star_Mass_He = star_Mass_He + DM*cur_HE ! star_Mass_He is total mass of He4
         star_Mass_C = star_Mass_C + DM*cur_C ! star_Mass_C is total mass of C12
         star_Mass_O = star_Mass_O + DM*MESH_Vs(V_XO,IK) ! star_Mass_O is total mass of O16
         star_Mass_Ne = star_Mass_Ne + DM*MESH_Vs(V_XNE,IK) ! star_Mass_Ne is total mass of Ne20
      END DO
      star_Mass_Mg = star_Mass*CMG*CZS
      star_Mass_Si = star_Mass*CSI*CZS
      star_Mass_Fe = star_Mass*CFE*CZS
      star_Mass_N = max(0D0,star_Mass*CZS-(star_Mass_C+star_Mass_O+star_Mass_Ne+star_Mass_Mg+star_Mass_Si+star_Mass_Fe)) ! N is leftover of metals
      center_H = MESH_Vs(V_XH,N_CNTR_shell) ! central fractional abundance by mass for H1
      center_He = MESH_Vs(V_XHE,N_CNTR_shell) ! for He4
      center_C = MESH_Vs(V_XC,N_CNTR_shell)  ! for C12
      center_O = MESH_Vs(V_XO,N_CNTR_shell)  ! for O16
      center_Ne = MESH_Vs(V_XNE,N_CNTR_shell) ! for Ne20
      center_Mg = CMG*CZS ! for Mg24
      center_Si = CSI*CZS ! for Si28
      center_Fe = CFE*CZS ! for Fe56
      center_N =  max(0D0, 1D0 - (center_H+center_He+center_C+center_O+center_Ne+center_Mg+center_Si+center_Fe)) ! for N14
      IF ( MESH_Vs(V_XH,N_SURF_shell) .LT. core_param_He ) mass_He_Core = star_Mass ! if XH never got to core_param_He
      initial_Mass = init_M
      initial_Z = init_Z
      mesh_core_mass = MC / CMSN
      END SUBROUTINE Complete_Model
      
      SUBROUTINE Get_Extra_Shell_Data ( IK )
      INTEGER, INTENT(IN) :: IK
      VAR(1:NUMV) = H(1:NUMV,IK)
      DVAR(1:NUMV) = DH(1:NUMV,IK)
      CALL FUNCS1 ( IK, 0 )
      SX(SX_PSI,IK) = PSI_from_LNF(VAR(V_LNF))
      SX(SX_DQ_DK,IK) = Q_dk ! dQ/dk, the change in the mesh spacing function
      SX(SX_P,IK) = P ! pressure
      SX(SX_SCP,IK) = SCP ! heat capacity
      SX(SX_RHO,IK) = RHO ! density
      SX(SX_T,IK) = T ! temperature
      SX(SX_OPACITY,IK) = FK ! opacity
      SX(SX_GRAD_AD,IK) = GRADA ! adiabatic temperature gradient
      SX(SX_GRAD_RAD,IK) = GRADR ! radiative temperature gradient
      SX(SX_GRAD_STAR,IK) = GRAD ! actual temperature gradient
      SX(SX_SG,IK) = SG*MK*1d66 ! diffusion coefficient for convective mixing
      SX(SX_GAMMA1,IK) = GAMMA1 ! adiabatic exponent
      SX(SX_GAM,IK) = GAM ! the plasma interaction parameter
      IF (IK .EQ. N_CNTR_shell) center_GAM = GAM
      SX(SX_XH,IK) = CAN(N_H)*NA(N_H)/AVM ! fractional abundances corrected for non-integer atomic weights
      SX(SX_XHE,IK) = CAN(N_HE)*NA(N_HE)/AVM
      SX(SX_XC,IK) = CAN(N_C)*NA(N_C)/AVM
      SX(SX_XN,IK) = CAN(N_N)*NA(N_N)/AVM
      SX(SX_XO,IK) = CAN(N_O)*NA(N_O)/AVM
      SX(SX_XNE,IK) = CAN(N_NE)*NA(N_NE)/AVM
      SX(SX_XMG,IK) = CAN(N_MG)*NA(N_MG)/AVM
      SX(SX_R,IK) = R/CRSN ! radius in solar units
      SX(SX_L,IK) = VAR(V_L)/CLSN ! luminosity in solar units
      SX(SX_EPS_NUC,IK) = EX ! nuclear energy generation rate (ergs/g/sec)
      SX(SX_EPS_NEU,IK) = -EN ! neutrino loss rate from non-nuclear reaction sources(ergs/g/sec)
      SX(SX_EPS_NUC_NEU,IK) = -ENX ! nuclear reaction neutrino loss rate (ergs/g/sec)
      SX(SX_EPS_G,IK) = ENG ! thermal energy generation rate (ergs/g/sec)
      SX(SX_U,IK) = U ! internal energy
      SX(SX_S,IK) = S ! entropy
      SX(SX_WL,IK) = WL ! diffusion coefficient for convective mixing
      SX(SX_MU,IK) = RHO/(P_ION/(CR*T) + NE1*RHO) ! mu, mean molecular weight per gas particle
      SX(SX_NE,IK) = NE1 ! free electron abundance (moles per gram of star)
      SX(SX_NE_TOT,IK) = NE ! NE if complete ionization
      SX(SX_XH2,IK) = XH2 ! mass fraction of molecular hydrogen
      SX(SX_XH_plus,IK) = XH_plus ! mass fraction of H+
      SX(SX_XHE_plus1,IK) = XHE_ions(1)
      SX(SX_XHE_plus2,IK) = XHE_ions(2)
      SX(SX_CV,IK) = WCV ! convection velocity
      ! components of pressure
      SX(SX_PRAD,IK) = PR
      SX(SX_PEL,IK) = PE
      SX(SX_PION,IK) = P_ION
      SX(SX_PCORR,IK) = PCORR
      ! calculate energy rates
         ! CME is (ergs/MeV)/(grams/amu) = (ergs*amu)/(MeV*grams)
         ! Rxx is reactions/baryon/sec
         ! Qxx is MeV/reaction
         ! Qxx*Rxx is MeV/baryon/sec
         ! CME*Qxx*Rxx is (ergs*amu)/(gram*baryon*sec) = (amu/baryon)*(ergs/g/sec)
         ! AVM is baryon/amu (?)
         ! CME*Qxx*Rxx/AVM is (ergs/g/sec)
      SX(SX_EPS_PP,IK) = (QPP + 0.5*Q33)*RPP*CME/AVM
      SX(SX_EPS_CNO,IK) = (QPC*RPC + QPO*RPO + QPNA*RPN + QPNG*RPNG)*CME/AVM
      SX(SX_EPS_H,IK) = SX(SX_EPS_PP,IK) + SX(SX_EPS_CNO,IK)
      SX(SX_EPS_3A,IK) = Q3A*R3A*CME/AVM
      SX(SX_EPS_AC,IK) = QAC*RAC*CME/AVM
      SX(SX_EPS_AN,IK) = QAN*RAN*CME/AVM
      SX(SX_EPS_AO,IK) = QAO*RAO*CME/AVM
      SX(SX_EPS_ANE,IK) = QANE*RANE*CME/AVM
      SX(SX_EPS_HE,IK) = SX(SX_EPS_3A,IK) + SX(SX_EPS_AC,IK) + SX(SX_EPS_AN,IK) + SX(SX_EPS_AO,IK) + SX(SX_EPS_ANE,IK)
      SX(SX_EPS_CCA,IK) = QCCA*RCC*CME/AVM
      SX(SX_EPS_CO,IK) = QCO*RCO*CME/AVM
      SX(SX_EPS_OO,IK) = QOO*ROO*CME/AVM
      SX(SX_EPS_GNE,IK) = QGNE*RGNE*CME/AVM
      SX(SX_EPS_CCG,IK) = QCCG*RCCG*CME/AVM
      SX(SX_EPS_GMG,IK) = QGMG*RGMG*CME/AVM
      SX(SX_EPS_Z,IK) = SX(SX_EPS_CCA,IK) + SX(SX_EPS_CO,IK) + SX(SX_EPS_CCG,IK) + SX(SX_EPS_OO,IK) + SX(SX_EPS_GNE,IK)
      SX(SX_NEU_plasma,IK) = eps_plasma
      SX(SX_NEU_brem,IK) = eps_brem
      SX(SX_NEU_pair,IK) = eps_pair
      SX(SX_NEU_photo,IK) = eps_photo
      SX(SX_PQ,IK) = PQ
      SX(SX_TQ,IK) = TQ
      SX(SX_RQ,IK) = RQ
      SX(SX_MQ,IK) = MQ
      SX(SX_DM_DK,IK) = MK / CMSN
      END SUBROUTINE Get_Extra_Shell_Data

      SUBROUTINE Get_Extra_Model_Data ! Get the stuff that requires calling FUNCS1
      use ez_magnitude_data
      use ez_magnitude
      
      DOUBLE PRECISION :: PEG, PDG, PDG_prev, prev_M, DM, WF, conv_time, MV, BMINV, UMINB
      INTEGER IK, IKK, IC
      
      double precision, parameter :: Zsol=0.02d0
      
      IF (have_Extras) RETURN ! don't redo it
      have_Extras = .TRUE.
      IC = 1 ! IC is the radiative/convective boundary index for arrays conv_boundaries_R and conv_turnover_Time
      power_3_alpha=0; power_C_alpha=0; power_N_alpha=0; power_O_alpha=0; power_Ne_alpha=0
      power_CC_Ne=0; power_CO=0; power_OO=0; power_Ne_decay=0; power_Mg_decay=0
      power_CC_Mg=0; power_PP=0; power_CNO=0 
      power_Neutrinos=0; power_plasmon_neutrinos=0; power_brem_neutrinos=0
      power_pair_neutrinos=0; power_photo_neutrinos=0
      power_H_burn=0D0; power_He_burn=0D0; power_Metal_burn=0D0
      conv_turnover_Time = 0D0; star_Mass_N = 0D0; prev_M = 0D0
      PDG = 0.d0 ! Added by Pat Scott 2008-01-29
      DO IK = N_CNTR_shell, N_SURF_shell, -1 ! go through the meshpoints from center to surface
         CALL Get_Extra_Shell_Data ( IK )
         DM = SX(SX_M,IK) - prev_M; prev_M = SX(SX_M,IK)
         PDG_prev = PDG; PDG = GRADR - GRADA
         IF ( IK .EQ. N_SURF_shell .AND. PDG .GT. 0D0 ) PDG = 0D0
         PEG = EG; IF ( IK .EQ. N_SURF_shell .AND. PEG .GT. 0D0 ) PEG = 0D0
         IF ( IK .NE. N_CNTR_SHELL ) THEN ! not at center-most meshpoint
            IKK = IK+1 ! the previous mesh point
            ! locate convective/radiative radius boundaries (DG = 0)
            IF ( IC.NE.BSZ ) THEN ! BSZ is max number of boundaries to save
               ! PDG is grad_r - grad_a for current meshpoint.  PDG_prev is same for previous
               IF (PDG_prev*PDG .LE. 0D0 ) THEN ! at a zone boundary
                  IC = IC + 1
                  ! interpolate to estimate radius where grad_r == grad_a
                  conv_boundaries_R(IC) = FIND0(SX(SX_R,IKK),PDG_prev,SX(SX_R,IK),PDG)
                  ! interpolate to estimate mass where grad_r == grad_a
                  conv_boundaries_M(IC) = FIND0(SX(SX_M,IKK),PDG_prev,SX(SX_M,IK),PDG)
                  ! calculate the convective time for the distance to boundary from meshpoint radius
                  IF ( PDG.GT.0D0 .AND. IC.LT.BSZ ) THEN ! start a new convective zone
                     conv_turnover_Time(IC+1) = (SX(SX_R,IK) - conv_boundaries_R(IC))/WCV ! WCV is convective velocity
                  END IF
                  IF ( PDG_prev.GT.0D0 ) THEN ! finish a zone (using previous convective velocity)
                     conv_turnover_Time(IC) = conv_turnover_Time(IC) + (conv_boundaries_R(IC) - SX(SX_R, IKK))/SX(SX_CV, IKK)
                  END IF
               ! Integrate dr/(conv. vel.), for convective turnover time
               ELSE IF ( WCV.GT.0D0 .AND. IC.LT.BSZ ) THEN ! inside a convective zone
                  ! use an intermediate velocity to calculate time for convection to cross this shell
                  conv_time = (SX(SX_R,IK)-SX(SX_R,IKK))*(WCV+SX(SX_CV,IKK))/(WCV*WCV+SX(SX_CV,IKK)*(SX(SX_CV,IKK)+WCV))
                  conv_turnover_Time(IC+1) = conv_turnover_Time(IC+1) + conv_time
               END IF
            WF = 1D0/DLOG10(SX(SX_P,IK)/SX(SX_P,IKK))                !   1/d_ln(P)
            SX(SX_HI1,IK) = WF*DLOG10(SX(SX_RHO,IK)/SX(SX_RHO,IKK))     !   d_ln(rho)/d_ln(P)
            SX(SX_HI2,IK) = -WF*DLOG10(SX(SX_R,IK)/SX(SX_R,IKK))        ! - d_ln(R)/d_ln(P)
            SX(SX_HI3,IK) = -WF*DLOG10(SX(SX_M,IK)/SX(SX_M,IKK))        ! - d_ln(M)/d_ln(P)
            END IF
         END IF
         ! Some integrated quantities
         WF = DM*CMSN/CLSN
         SX(SX_WF,IK) = WF ! WF converts from (ergs/g/sec) to L/Lsolar for this shell
         ! SX(SX_DM,IK)*CMSN is mass in units of 10^33 grams
         ! multiply this by ergs/g/sec to get luminosity in 10^33 ergs/sec units
         ! divide by CLSN to convert to solar luminosity units
         star_Mass_N = star_Mass_N + DM*SX(SX_XN,IK)           ! star_Mass_N is total mass of N14
         power_PP = power_PP + WF*SX(SX_EPS_PP,IK)
         power_CNO = power_CNO + WF*SX(SX_EPS_CNO,IK)
         power_H_burn = power_H_burn + WF*SX(SX_EPS_H,IK) ! power_H_burn is total luminosity from hydrogen burning
         power_CC_Ne = power_CC_Ne + WF*SX(SX_EPS_CCA,IK) ! total luminosity from 2 C12 -> Ne20 + He4
         power_CO = power_CO + WF*SX(SX_EPS_CO,IK) ! total luminosity from C12 + O16 -> Mg24 + He4
         power_OO = power_OO + WF*SX(SX_EPS_OO,IK) ! total luminosity from 2 O16 -> Mg24 + 2 He4
         power_Ne_decay = power_Ne_decay + WF*SX(SX_EPS_GNE,IK) ! total luminosity from Ne20 -> O16 + He4
         power_Mg_decay = power_Mg_decay + WF*SX(SX_EPS_GMG,IK) ! total luminosity from Mg24 -> Ne20 + He4
         power_CC_Mg = power_CC_Mg + WF*SX(SX_EPS_CCG,IK) ! total luminosity from 2 C12 -> Mg24
         power_Metal_burn = power_Metal_burn + WF*SX(SX_EPS_Z,IK) ! power_Metal_burn is total luminosity from burning metals but not helium.
         power_plasmon_neutrinos = power_plasmon_neutrinos + WF*SX(SX_NEU_plasma,IK)
         power_brem_neutrinos = power_brem_neutrinos + WF*SX(SX_NEU_brem,IK)
         power_pair_neutrinos = power_pair_neutrinos + WF*SX(SX_NEU_pair,IK)
         power_photo_neutrinos = power_photo_neutrinos + WF*SX(SX_NEU_photo,IK)
         power_3_alpha = power_3_alpha + WF*SX(SX_EPS_3A,IK) ! power_3_alpha is total luminosity from triple alpha reaction
         power_C_alpha = power_C_alpha + WF*SX(SX_EPS_AC,IK) ! power_C_alpha is total luminosity from C12 + He4 -> O16
         power_N_alpha = power_N_alpha + WF*SX(SX_EPS_AN,IK) ! power_N_alpha is total luminosity from N14 + 3/2 He4 -> Ne20
         power_O_alpha = power_O_alpha + WF*SX(SX_EPS_AO,IK) ! power_O_alpha is total luminosity from O16 + He4 -> Ne20
         power_Ne_alpha = power_Ne_alpha + WF*SX(SX_EPS_ANE,IK) ! power_Ne_alpha is total luminosity from Ne20 + He4 -> Mg24
         power_He_burn = power_He_burn + WF*SX(SX_EPS_HE,IK) ! power_He_burn is total luminosity from helium burning
         power_Neutrinos = power_Neutrinos - WF*EN ! power_Neutrinos is total neutrino loss
         SX(SX_L_H,IK) = power_H_burn
         SX(SX_L_HE,IK) = power_He_burn
         SX(SX_L_Z,IK) = power_Metal_burn
         SX(SX_L_NEU,IK) = power_Neutrinos
      END DO
      conv_boundary_count = IC

      call Get_Magnitudes(log_surface_Temp, log_Luminosity, star_Mass, log10(initial_Z / Zsol), color_magnitudes)
      
      END SUBROUTINE Get_Extra_Model_Data


      END MODULE ez_report
