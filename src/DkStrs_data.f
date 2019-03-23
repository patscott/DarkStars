! Shared variables module for DarkStars
!
! Pat Scott, July 2007; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE DkStrs_data
      use DkStrs_constants
      implicit none

      include 'dscapstar.h'
      include 'dsidtag.h'

      character (len=strlen) :: main_dir, prof_dir, movie_dir, data_dir_in, metals_dir_in, save_filename_in
      logical :: WIMPdensmode, DoTransport, DoWimpyThings, LTE_invalid, loopOrbit, WimpyThingsDonePreviously
      logical :: new_model, new_model_step1, after_restore, NotATest, verbose
      integer :: fileLen, JMOD_prev, JMOD_new_model_start, summary_cnt_in, stop_at_Model_in
      integer(8) :: start_time=-1, save_count=0
      double precision :: save_freq
      double precision :: cap, ann, evap, effcap
      double precision :: solar_C, solar_O
      double precision :: heavy_abuns(n_heavy_abuns)
      double precision :: iso_ratios(n_iso_ratios)
      double precision :: sigsi, sigsd, sigann, mx
      double precision :: n_WIMPs, n_WIMPs_prev, n_WIMPs_saved, n_WIMPs_supersaved, dNdt 
      double precision :: r_chi, Tw, condEff, hold_on, prec
      double precision :: K_suppression, K0, Knudsen_suppression_tau
      double precision :: nuLossFactor, boost_factor, min_timestep, timestep_rescale
      double precision :: orbit_len, orbit_monotonicity_scale 
      double precision :: prelim_Age, prelim_Age2, new_model_age
      double precision :: rhoonv_prev, rhoonv_saved, tau_therm
      double precision :: stellar_age(5), centreH(5), centrerho(5), centreT(5)
      real(16)         :: partition_func, n_WIMPs_centre_LTE

      double precision, allocatable :: orbital_t(:), orbital_rho(:), orbital_v(:)
      double precision, allocatable :: orbital_galesc(:), orbital_galr(:)
      double precision, allocatable :: orbital_rho_d1(:), orbital_rho_tension(:)
      double precision, allocatable :: orbital_v_d1(:), orbital_v_tension(:)
      double precision, allocatable :: orbital_galesc_d1(:), orbital_galesc_tension(:)
      double precision, allocatable :: orbital_galr_d1(:), orbital_galr_tension(:)

      double precision :: alpha_i(n_species), kappa_i(n_species), sigma_i(n_species)

      double precision :: staralpha(meshpoints), staralphad2(meshpoints)
      double precision :: stardens(meshpoints), stardensd2(meshpoints)
      double precision :: startemperature(meshpoints), startemperatured1(meshpoints), 
     & temptension(meshpoints)
      double precision :: starLTErhoWIMP(meshpoints), starLTErhoWIMPd1(meshpoints), 
     & LTErhoWIMPtension(meshpoints)
      double precision :: L_WIMPtransport(meshpoints), L_WIMPtransportd1(meshpoints), 
     & L_WIMPtransporttension(meshpoints)
      double precision :: starcondEff(meshpoints) = 0, starcondEffd1(meshpoints) = 0, 
     & condEfftension(meshpoints) = 0

      END MODULE DkStrs_data
