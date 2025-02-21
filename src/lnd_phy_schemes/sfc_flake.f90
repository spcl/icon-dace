! The main program unit of the lake parameterization scheme FLake.
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------
!
!
! The main program unit of the lake parameterization scheme FLake.
!
! It contains all FLake procedures except for the procedures used
! to compute fluxes of momentum and of sensible and latent heat over lakes.
! (Note that the FLake flux-calculation routines can optionally be used
! to determine fluxes over lakes. By default, fluxes over lakes
! are computed in the same way as over the ocean).
! Communication between FLake and the host atmospheric model (ICON)
! occurs through the subroutines "flake_coldinit", "flake_init" and "flake_interface".
! In "flake_init" (warm start initialization procedure), FLake prognostic variables
! are initialized and some consistency checks are performed.
! In "flake_interface", a call to the FLake subroutine "flake_driver" is organized,
! where FLake variables are advanced one time step forward
! (see "flake_driver" for further details).
! The "flake_coldinit" serves to perform a cold start of FLake
! within a host model.
!
! FLake (Fresh-water Lake) is a bulk lake parameterization scheme
! capable of predicting the water temperature
! and the characteristics of ice and snow in lakes of various depth
! on the time scales from a few hours to many years.
! The scheme is based on a parametric representation of the evolving
! temperature profile, where the structure of the stratified layer between the
! upper mixed layer and the basin bottom, the lake thermocline, is described
! using the concept of self-similarity (assumed shape) of the temperature-depth curve.
! The concept was put forward by Kitaigorodskii and Miropolsky (1970) to
! describe the vertical temperature structure of the oceanic seasonal
! thermocline. It has been successfully used in geophysical applications.
! A bulk approach based on the assumed shape of the temperature-depth curve
! is also used to describe the vertical structure of the thermally active
! upper layer of bottom sediments and of the ice and snow cover.
!
! FLake incorporates the heat budget equations for the four layers
! in question, viz., snow, ice, water and bottom sediments, developed with
! due regard for the vertically distributed character of solar radiation
! heating. The entrainment equation that incorporates the spin-up term
! is used to compute the depth of a convectively-mixed layer.
! A relaxation-type equation is used to compute the wind-mixed layer depth in
! stable and neutral stratification, where a multi-limit formulation for the
! equilibrium mixed-layer depth accounts for the effects of the earth's rotation,
! of the surface buoyancy flux and of the static stability in the thermocline.
! The equations for the mixed-layer depth are developed with due regard for
! the volumetric character of the solar radiation heating.
! Simple thermodynamic arguments are invoked to develop
! the evolution equations for the ice and snow thickness.
! The heat flux through the water-bottom sediment interface is computed,
! using a parameterization proposed by Golosov et al. (1998).
! The heat flux trough the air-water interface
! (or through the air-ice or air-snow interface)
! is provided by the driving atmospheric model.
!
! Disposable constants and parameters of FLake are estimated, using
! independent empirical and numerical data. They should not be re-evaluated
! when the scheme is applied to a particular lake. The only lake-specific
! parameters are the lake depth, the optical characteristics of lake water,
! the temperature at the lower boundary of the thermally active layer of bottom
! sediments, and the depth of that layer. These parameters are not part
! the model physics, however.
! Note though that the data assimilation procedures (e.g. assimilation of observational
! data on ice fraction) utilize constants and parameters, whose estimates are essentially
! based on prior experience and common sense. Those constants and parameters
! may need to be tuned as the case requires.
!
! A detailed description of FLake is given in
!
! Mironov, D. V., 2008:
! Parameterization of lakes in numerical weather prediction. Description of a lake model.
! COSMO Technical Report, No. 11, Deutscher Wetterdienst, Offenbach am Main, Germany, 41 pp.
!
! Mironov, D., E. Heise, E. Kourzeneva, B. Ritter, N. Schneider, and A. Terzhevik, 2010:
! Implementation of the lake parameterisation scheme FLake
! into the numerical weather prediction model COSMO.
! Boreal Env. Res., 15, 218-230.
!
! A snow-ice module of FLake is discussed in
!
! Mironov, D., B. Ritter, J.-P. Schulz, M. Buchhold, M. Lange, and E. Machulskaya, 2012:
! Parameterization of sea and lake ice in numerical weather prediction models
! of the German Weather Service.
! Tellus A, 64, 17330. doi:10.3402/tellusa.v64i0.17330
!
! Further information is available at the FLake web page http://lakemodel.net.
!
! In the present configuration, the bottom sediment module of FLake is switched off
! and the heat flux at the water-bottom sediment interface is set to zero.
! Although the snowfall rate is provided by the driving atmospheric model,
! snow over lake ice is not considered explicitly.
! The effect of snow is accounted for implicitly (parametrically)
! through changes in the surface albedo with respect to solar radiation.
!
!
! Lines embraced with "!_tmp>" and "!_tmp<" contain temporary parts of the code.
! Lines embraced/marked with "!_dev>" and "!_dev<" may be replaced
! as improved formulations are developed and tested.
! Lines embraced/marked with "!_cdm>" and "!_cdm<" are DM's comments that may be helpful to a user.
! Lines embraced/marked with "!_dbg>" and "!_dbg<" are used for debugging purposes only.
! Lines starting with "!_nu" are not used.

!NEC$ options "-finline-max-function-size=1000 -finline-max-depth=3"


MODULE sfc_flake

!===================================================================================================

  USE mo_kind,           ONLY : wp     !< KIND-type parameter for real variables

  USE mo_lnd_nwp_config, ONLY : frlake_thrhld  !< fraction threshold for creating a lake grid point

!===================================================================================================

  USE sfc_flake_data, ONLY: &
    ! flake configure
    &  lflk_botsed_use   , & !< .TRUE. indicates that bottom-sediment scheme is used
    &  rflk_depth_bs_ref     !< Reference value of depth of thermally active layer 
                             !< of bottom sediments [m]

  USE sfc_flake_data, ONLY:        &
    ! derived types
    &  opticpar_medium          , & !< derived type 
    ! optical characteristics of water, ice and snow 
    &  opticpar_water_ref       , & !< Water (reference)
!_nu    &  opticpar_water_trans     , & !< Transparent water (two-band)
!_nu    &  opticpar_whiteice_ref    , & !< White ice
!_nu    &  opticpar_blueice_ref     , & !< Blue ice
!_nu    &  opticpar_drysnow_ref     , & !< Dry snow
!_nu    &  opticpar_meltingsnow_ref , & !< Melting snow
    &  opticpar_ice_opaque      , & !< Opaque ice
    &  opticpar_snow_opaque         !< Opaque snow

  USE sfc_flake_data, ONLY: &
    ! flake parameters
    &  c_cbl_1           , & !< Constant in the CBL entrainment equation [-]
    &  c_cbl_2           , & !< Constant in the CBL entrainment equation [-]
    &  c_sbl_ZM_n        , & !< Constant in the ZM1996 equation for the equilibrium SBL depth [-]
    &  c_sbl_ZM_s        , & !< Constant in the ZM1996 equation for the equilibrium SBL depth [-]
    &  c_sbl_ZM_i        , & !< Constant in the ZM1996 equation for the equilibrium SBL depth [-]
    &  c_relax_h         , & !< Constant in the relaxation equation for the SBL depth [-]
    &  c_relax_C         , & !< Constant in the relaxation equation for the shape factor 
                             !< with respect to the temperature profile in the thermocline [-]
    &  C_T_min           , & !< Minimum value of the shape factor C_T (thermocline) [-]
    &  C_T_max           , & !< Maximum value of the shape factor C_T (thermocline) [-]
    &  Phi_T_pr0_1       , & !< Constant in the expression for the T shape-function derivative [-]
    &  Phi_T_pr0_2       , & !< Constant in the expression for the T shape-function derivative [-]
    &  C_TT_1            , & !< Constant in the expression for C_TT (thermocline) [-]
    &  C_TT_2            , & !< Constant in the expression for C_TT (thermocline) [-]
    &  C_B1              , & !< Shape factor (upper layer of bottom sediments) [-]
    &  C_B2              , & !< Shape factor (lower layer of bottom sediments) [-]
    &  Phi_B1_pr0        , & !< B1 shape-function derivative [-]
    &  C_S_lin           , & !< Shape factor (linear temperature profile in the snow layer) [-]
    &  Phi_S_pr0_lin     , & !< S shape-function derivative (linear profile) [-]
    &  C_I_lin           , & !< Shape factor (linear temperature profile in the ice layer) [-]
    &  Phi_I_pr0_lin     , & !< I shape-function derivative (linear profile) [-]
    &  Phi_I_pr1_lin     , & !< I shape-function derivative (linear profile) [-]
    &  Phi_I_ast_MR      , & !< Constant in the MR2004 expression for I shape factor [-]
    &  C_I_MR            , & !< Constant in the MR2004 expression for I shape factor [-]
    &  H_Ice_max         , & !< Maximum ice thickness in the Mironov/Ritter ice model [m]
    &  h_Snow_min_flk    , & !< Minimum snow thickness [m]
    &  h_Ice_min_flk     , & !< Minimum ice thickness [m]
    &  h_ML_min_flk      , & !< Minimum mixed-layer depth [m]
    &  h_ML_max_flk      , & !< Maximum mixed-layer depth [m]
    &  H_B1_min_flk      , & !< Minimum thickness of the upper layer of bottom sediments [m]
    &  u_star_min_flk    , & !< Minimum value of the surface friction velocity [m s^{-1}]
    &  c_small_flk       , & !< A small number
    &  c_maxearg_flk         !< Maximum value of the EXP function argument [-]

!_cdm>
! Note that most physical constants are taken from the ICON module 
! "mo_physical_constants" rather than from "sfc_flake_data" 
! (the respective lines are marked with "!_nu").
! FLake-specific parameters are in fact parameters in empirical approximation formulae
! for some physical properties of different media (e.g. snow).
!_cdm<
  USE sfc_flake_data, ONLY: &
    ! flake parameters and physical constants
!_nu    &  tpl_grav          , & !< Acceleration due to gravity [m s^{-2}]
    &  tpl_T_r           , & !< Temperature of maximum density of fresh water [K]
!_nu    &  tpl_T_f           , & !< Fresh water freezing point [K]
    &  tpl_a_T           , & !< Constant in the fresh-water equation of state [K^{-2}]
!_nu    &  tpl_rho_w_r       , & !< Maximum density of fresh water [kg m^{-3}]
!_nu    &  tpl_rho_I         , & !< Density of ice [kg m^{-3}]
    &  tpl_rho_S_min     , & !< Minimum snow density [kg m^{-3}]
    &  tpl_rho_S_max     , & !< Maximum snow density [kg m^{-3}]
    &  tpl_Gamma_rho_S   , & !< Empirical parameter [kg m^{-4}] 
                             !< in the expression for the snow density
!_nu    &  tpl_L_f           , & !< Latent heat of fusion  [J kg^{-1}]
!_nu    &  tpl_c_w           , & !< Specific heat of water [J kg^{-1} K^{-1}]
!_nu    &  tpl_c_I           , & !< Specific heat of ice   [J kg^{-1} K^{-1}]
!_nu    &  tpl_c_S           , & !< Specific heat of snow  [J kg^{-1} K^{-1}]
    &  tpl_kappa_w       , & !< Molecular heat conductivity of water        [J m^{-1} s^{-1} K^{-1}]
!_nu    &  tpl_kappa_I       , & !< Molecular heat conductivity of ice          [J m^{-1} s^{-1} K^{-1}]
    &  tpl_kappa_S_min   , & !< Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
    &  tpl_kappa_S_max   , & !< Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
    &  tpl_Gamma_kappa_S     !< Empirical parameter in expression for the 
                             !< snow heat conductivity [J m^{-2} s^{-1} K^{-1}]
!_cdm>
! Most physical constants are taken from the ICON module "mo_physical_constants" 
! rather than from the FLake module "sfc_flake_data". 
! Note that values of some constants are slightly different 
! in "mo_physical_constants" and "sfc_flake_data".
! FLake: tpl_grav    = 9.81,    ICON: grav  = 9.80665.
! FLake: tpl_rho_I   = 9.1E+02, ICON: rhoi = 917.0.
! FLake: tpl_L_f     = 3.3E+05, ICON: alf  = 3.337E+05.
! FLake: tpl_c_w     = 4.2E+03, ICON: clw  = 4192.664112.
! FLake: tpl_c_I     = 2.1E+03, ICON: ci   = 2106.0.
! FLake: tpl_c_S     = 2.1E+03, ICON: ci   = 2106.0.
! FLake: tpl_kappa_I = 2.29,    ICON: ki   = 2.1656.
!_cdm<

  USE mo_physical_constants, ONLY: &
    &  tpl_grav    => grav       , & !< Acceleration due to gravity [m s^{-2}]
    &  tpl_T_f     => tmelt      , & !< Fresh-water freezing point [K]
    &  tpl_rho_w_r => rhoh2o     , & !< Maximum density of fresh water [kg m^{-3}]
    &  tpl_rho_I   => rhoi       , & !< Density of ice [kg/m^3]
    &  tpl_L_f     => alf        , & !< Latent heat of fusion  [J kg^{-1}]
    &  tpl_c_w     => clw        , & !< Specific heat of water [J kg^{-1} K^{-1}]
    &  tpl_c_I     => ci         , & !< Specific heat of ice   [J kg^{-1} K^{-1}]
    &  tpl_c_S     => ci         , & !< Specific heat of snow  [J kg^{-1} K^{-1}]
    &  tpl_kappa_I => ki             !< Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]

  USE mo_exception, ONLY:       &
                    &  finish , &  !< external procedure, finishes model run and reports the reason
                    &  message     !< external procedure, sends a message (error, warning, etc.)

  USE mo_fortran_tools, ONLY: set_acc_host_or_device

!===================================================================================================

!_cdm>
! MESSY stuff from COSMO. Currently unused within ICON.
!_cdm<
!_nu
!_nu #ifndef MESSY
!_nu   USE src_tracer,       ONLY: trcr_get, trcr_errorstr
!_nu #else
!_nu   USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
!_nu   USE messy_main_tracer_bi,     ONLY: tracer_halt
!_nu   USE messy_main_tracer,        ONLY: get_tracer
!_nu #endif

!===================================================================================================

  IMPLICIT NONE

  PRIVATE


  ! The variables declared below are accessible to all program units of the MODULE "sfc_flake".
  ! These are basically variables handled "internally" by FLake routines.
  ! All variables declared below have a suffix "flk".

  ! FLake variables of type REAL

  ! These variables are declared here so that we don't need to allocate
  ! them in OpenACC every time
  TYPE (opticpar_medium) ::  & 
  &  opticpar_water        , & !< optical characteristics of water
  &  opticpar_ice          , & !< optical characteristics of ice
  &  opticpar_snow             !< optical characteristics of snow 
  !$OMP THREADPRIVATE(opticpar_water, opticpar_ice, opticpar_snow)

  ! Entities that should be accessible from outside the present module 
  PUBLIC ::                       &
         &  flake_coldinit      , & ! procedure (cold start initialization)
         &  flake_init          , & ! procedure (initialization)
         &  flake_interface         ! procedure (time stepping)




  LOGICAL, PARAMETER :: using_cuda_graph = .FALSE.


!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

CONTAINS

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

!===================================================================================================
!  Initialize the lake parameterization scheme FLake (warm start initialization)
!---------------------------------------------------------------------------------------------------

  !> 
  !! This routine initializes prognostic variables of the lake model FLake.
  !! At present, the bottom-sediment module of FLake is switched off so that 
  !! the external parameters of the bottom-sediment module,
  !! namely, the depth of the thermally active layer of the sediment
  !! and the climatological temperature at the bottom of that layer,
  !! are not required.
  !! These parameters are set to their reference values and are not part of the input 
  !! (not read from external-parameter files).
  !! The same applies to the attenuation coefficient with respect to solar 
  !! radiation and to the wind fetch. Default values of these parameters are used.
  !! Three FLake prognostic variables, namely, 
  !! the snow thickness (always zero as snow over lake ice is not considered 
  !! explicitly at present), thickness of the upper layer of bottom sediment and 
  !! the temperature at the outer edge of that layer 
  !! (these are equal to their reference values as the bottom-sediment module is 
  !! switched off), are not included into the ICON IO list.
  !! These variable are handled internally.
  !! Observational data on ice fraction are utilized to create new ice, reduce the 
  !! ice thickness, or remove ice. To this end, an ad hoc assimilation procedure
  !! is used (that contains a good few of tuning constants). 
  !! Data on lake ice fraction are currently used for some large lakes only.
  !!

  SUBROUTINE flake_init (                                       &
                     &  nflkgb, use_iceanalysis,                &
                     &  fr_lake, depth_lk, fr_ice,              &
                     &  fetch_lk, dp_bs_lk, t_bs_lk, gamso_lk,  &
                     &  t_snow_p, h_snow_p,                     & 
                     &  t_ice_p, h_ice_p,                       & 
                     &  t_mnw_lk_p, t_wml_lk_p, t_bot_lk_p,     &
                     &  c_t_lk_p, h_ml_lk_p,                    & 
                     &  t_b1_lk_p, h_b1_lk_p,                   &           
                     &  t_scf_lk_p                             )

    IMPLICIT NONE

    ! Procedure arguments

    ! Array dimension(s)

    INTEGER, INTENT(IN) ::         &
                        &  nflkgb    !< array (vector) dimension,
                                     !< equal to the number of grid boxes within a block 
                                     !< where lakes are present 

    ! FLake external parameters

    LOGICAL, DIMENSION(:), INTENT(IN) :: use_iceanalysis !< switch whether to use ice fraction analysis
                                                         !  on a given grid point

    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  &
                        &  fr_lake          , & !< lake fraction in a grid box [-]
                        &  depth_lk         , & !< lake depth [m]
                        &  fr_ice               !< ice fraction coming from analysis

    REAL(KIND = wp)    , DIMENSION(:), INTENT(INOUT) ::  &
                        &  fetch_lk         , & !< wind fetch over lake [m]
                        &  dp_bs_lk         , & !< depth of the thermally active layer
                                                !< of bottom sediments [m]
                        &  t_bs_lk          , & !< climatological temperature at the bottom of
                                                !< the thermally active layer of sediments [K]
                        &  gamso_lk             !< attenuation coefficient of the lake water 
                                                !< with respect to solar radiation [m^{-1}]

    ! FLake prognostic variables
    ! (at the previous time step - "p", and the updated values - "n")

    REAL(KIND = wp)    , DIMENSION(:), INTENT(INOUT) ::  &
                        &  t_snow_p         , & !< temperature of snow upper surface at previous time step [K]
                        &  h_snow_p         , & !< snow thickness at previous time step [m]
                        &  t_ice_p          , & !< temperature of ice upper surface at previous time step [K]
                        &  h_ice_p          , & !< ice thickness at previous time step [m]
                        &  t_mnw_lk_p       , & !< mean temperature of the water column at previous time step [K]
                        &  t_wml_lk_p       , & !< mixed-layer temperature at previous time step [K] 
                        &  t_bot_lk_p       , & !< temperature at the water-bottom sediment interface 
                                                !< at previous time step [K] 
                        &  c_t_lk_p         , & !< shape factor with respect to the temperature profile 
                                                !< in lake thermocline at previous time step [-] 
                        &  h_ml_lk_p        , & !< thickness of the mixed-layer at previous time step [m] 
                        &  t_b1_lk_p        , & !< temperature at the bottom of the upper layer
                                                !< of the sediments at previous time step [K]  
                        &  h_b1_lk_p        , & !< thickness of the upper layer of bottom sediments 
                                                !< at previous time step [m] 
                        &  t_scf_lk_p           !< lake surface temperature at previous time step [K]
                                                !< (i.e. the temperature at the air-water, air-ice 
                                                !< or air-snow interface) 
                       
    ! Local variables

    INTEGER ::       &
            &  iflk  !< DO loop index

    CHARACTER(len=256) ::            &
                       &  nameerr  , &  !< name of procedure where an error occurs
                       &  texterr       !< error/warning message text

    LOGICAL ::             &
            &  lcallabort  !< logical switch, set .TRUE. if errors are encountered
                           !< (used to call model abort outside a DO loop)

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Logical switch, default value is .FALSE.
    lcallabort = .FALSE.

    ! Loop over grid boxes where lakes are (should be) present
    CheckFLakeExtPar: DO iflk=1, nflkgb
      ! Check lake-fraction and lake-depth fields
      IF( (fr_lake(iflk) < frlake_thrhld) .OR. (depth_lk(iflk) < 0._wp) ) THEN
        ! Lake fraction less than a minimum threshold value or negative lake depth is found
        ! Set logical switch
        lcallabort = .TRUE.

        ! Exit DO loop to call model abort
        EXIT CheckFLakeExtPar 
      END IF
    END DO CheckFLakeExtPar 

    ! Call model abort if errors are encountered
    IF( lcallabort ) THEN
      ! Send an error message
      WRITE(nameerr,*) "MODULE sfc_flake, SUBROUTINE flake_init"
      WRITE(texterr,*) "Lake fraction ", fr_lake(iflk),                          &
                    &  " is less than a minimum threshold value ", frlake_thrhld
      CALL message(TRIM(nameerr), TRIM(texterr))
      WRITE(texterr,*) " or negative lake depth ", depth_lk(iflk),               &
                    &  " Call model abort."
      CALL message(TRIM(nameerr), TRIM(texterr))
      ! Call model abort
      WRITE(nameerr,*) "sfc_flake:flake_init"
      WRITE(texterr,*) "error in lake fraction or lake depth"
      CALL finish(TRIM(nameerr), TRIM(texterr))
    END IF

    ! As different from lake fraction and lake depth 
    ! that are part of the ICON (COSMO) IO (read from GRIB, NetCDF, etc., files), 
    ! other FLake external-parameter fields are handled internally by the host atmospheric model.
    ! Set those FLake external parameters to reference values.
    SetFLakeExtPar: DO iflk=1, nflkgb
      fetch_lk(iflk) = 1.0E+04_wp        ! Use a constant long fetch [m]
      dp_bs_lk(iflk) = rflk_depth_bs_ref     ! Reference value [m]
      t_bs_lk (iflk) = tpl_T_r               ! Reference value [K]
      gamso_lk(iflk) = opticpar_water_ref%extincoef_optic(1)
                                             ! Reference value [m^{-1}]
    END DO SetFLakeExtPar 

    ! Since a loss of accuracy may occur during IO
    ! (e.g. due to GRIB encoding/decoding),
    ! FLake prognostic variables must be checked for consistency.
    ! That is, the relations between
    ! "t_wml_lk", "t_mnw_lk", "t_bot_lk", "h_ml_lk", "c_t_lk", "h_ice" and "t_ice",
    ! suggested by the assumed shape (self-similarity) of the temperature-depth curve,
    ! must hold or should be enforced.

    GridBoxesWithLakes: DO iflk=1, nflkgb

      ! FLake variables at previous time step

      ! Limit the mixed-layer depth and the temperature-profile shape factor
      h_ml_lk_p(iflk) = MIN(depth_lk(iflk), MAX(h_ml_lk_p(iflk), 0._wp))
      c_t_lk_p(iflk)  = MIN(C_T_max, MAX(c_t_lk_p(iflk), C_T_min))
      IF(h_ice_p(iflk) >= h_Ice_min_flk) THEN
        ! Ice-covered lakes
        h_ice_p(iflk)    = MIN(h_ice_p(iflk), H_Ice_max)
        t_ice_p(iflk)    = MIN(t_ice_p(iflk), tpl_T_f)
        t_wml_lk_p(iflk) = tpl_T_f
        t_mnw_lk_p(iflk) = MAX(tpl_T_f, MIN(t_mnw_lk_p(iflk), tpl_T_r))
        t_bot_lk_p(iflk) = MAX(tpl_T_f, MIN(t_bot_lk_p(iflk), tpl_T_r))
        ! Mixing down to the lake bottom or "inverse" stratification,
        ! reset mixed-layer depth, shape factor, mean temperature, and bottom temperature.
        IF( (h_ml_lk_p(iflk) >= depth_lk(iflk)-h_ML_min_flk) .OR.  &
            (t_bot_lk_p(iflk) <= t_wml_lk_p(iflk)) ) THEN
          h_ml_lk_p(iflk)  = 0._wp
          c_t_lk_p(iflk)   = C_T_min
          t_mnw_lk_p(iflk) = t_wml_lk_p(iflk)
          t_bot_lk_p(iflk) = t_wml_lk_p(iflk)
        END IF
      ELSE
        ! Ice-free lakes 
        h_ice_p(iflk)    = 0.0_wp
        t_ice_p(iflk)    = tpl_T_f
        t_wml_lk_p(iflk) = MAX(tpl_T_f, t_wml_lk_p(iflk))
        t_mnw_lk_p(iflk) = MAX(tpl_T_f, t_mnw_lk_p(iflk))
        t_bot_lk_p(iflk) = MAX(tpl_T_f, t_bot_lk_p(iflk))
        ! Avoid temperature of maximum density crossover
        IF( ((t_wml_lk_p(iflk) <= tpl_T_r) .AND. (t_bot_lk_p(iflk) > tpl_T_r)) .OR.  &
             ((t_wml_lk_p(iflk) >= tpl_T_r) .AND. (t_bot_lk_p(iflk) < tpl_T_r)) )    &
          t_bot_lk_p(iflk) = tpl_T_r
        ! Mixing down to the lake bottom or "inverse" stratification,
        ! reset mixed-layer depth, shape factor, mixed-layer temperature, and bottom temperature.
        IF( (h_ml_lk_p(iflk) >= depth_lk(iflk)-h_ML_min_flk) .OR.                               &
             ((t_wml_lk_p(iflk) <= tpl_T_r) .AND. (t_bot_lk_p(iflk) < t_wml_lk_p(iflk))) .OR.  &
             ((t_wml_lk_p(iflk) >= tpl_T_r) .AND. (t_bot_lk_p(iflk) > t_wml_lk_p(iflk))) ) THEN
          h_ml_lk_p(iflk)  = depth_lk(iflk)
          c_t_lk_p(iflk)   = C_T_min
          t_wml_lk_p(iflk) = t_mnw_lk_p(iflk)
          t_bot_lk_p(iflk) = t_mnw_lk_p(iflk)
        END IF
      END IF
      ! Snow over lake ice is not considered explicitly,
      ! set "h_snow" to zero and "t_snow" to the ice surface temperature.
      h_snow_p(iflk)   = 0.0_wp
      t_snow_p(iflk)   = t_ice_p(iflk)
      ! Bottom-sediment module is switched off 
      ! and the respective FLake variables are handled internally, 
      ! i.e. they are not part of ICON (COSMO) IO. 
      ! Set the bottom-sediment FLake variables to reference values.
      t_b1_lk_p(iflk) = tpl_T_r          
      h_b1_lk_p(iflk) = rflk_depth_bs_ref
      ! Set the lake surface temperature
      IF (h_snow_p(iflk) >= h_Snow_min_flk) THEN
        ! Snow exists, use the snow surface temperature
        t_scf_lk_p(iflk) = t_snow_p(iflk)
      ELSE IF (h_ice_p(iflk) >= h_Ice_min_flk) THEN
        ! Ice exists but there is no snow, use the ice surface temperature
        t_scf_lk_p(iflk) = t_ice_p(iflk)
      ELSE
        ! No ice-snow cover, use the water surface temperature ( = mixed-layer temperature)
        t_scf_lk_p(iflk) = t_wml_lk_p(iflk)
      END IF

      ! Adapt ice thickness and temperatures if analysis data of ice fraction are available
      AssimIceFractionData: IF (use_iceanalysis(iflk)) THEN

        IF (fr_ice(iflk) > 0.05_wp .AND. h_ice_p(iflk) < h_Ice_min_flk) THEN
          ! There was no ice in the first guess, create new ice with an amount depending on ice fraction
          h_ice_p(iflk)  = 0.025_wp*fr_ice(iflk)
          ! Set the ice surface temperature to the fresh-water freezing point
          t_ice_p(iflk)  = tpl_T_f
          ! Set the mixed-layer temperature to the fresh-water freezing point
          t_wml_lk_p(iflk)  = tpl_T_f
          ! Adjust the temperature profile 
          ! Reset h_ML and C_T as needed
          IF(h_ml_lk_p(iflk) >= (depth_lk(iflk)-h_ML_min_flk)) THEN ! h_ML=D when ice is created 
            h_ml_lk_p(iflk) = 0.0_wp           ! Set h_ML to zero 
            c_t_lk_p(iflk)  = C_T_min          ! Set C_T to its minimum value 
          END IF                               ! h_ML<D when ice is created 
          ! Limit the bottom temperature
          t_bot_lk_p(iflk) = MAX(tpl_T_f, MIN(t_bot_lk_p(iflk), tpl_T_r))
          ! Compute the mean temperature of the water column
          t_mnw_lk_p(iflk) = t_wml_lk_p(iflk) - c_t_lk_p(iflk)                             &
            &              * MAX(0.0_wp, (1.0_wp - h_ml_lk_p(iflk)/depth_lk(iflk)))        &
            &              * (t_wml_lk_p(iflk) - t_bot_lk_p(iflk))
        ELSE IF (fr_ice(iflk) >= 0.03_wp .AND. fr_ice(iflk) < 0.75_wp                      &
          &      .AND. h_ice_p(iflk) >= 0.1_wp*fr_ice(iflk)) THEN
          ! Ice exists but should be reduced in depth 
          h_ice_p(iflk)  = 0.1_wp*fr_ice(iflk)
          ! Set the ice surface temperature to the fresh-water freezing point if the ice fraction is < 0.1
          IF (fr_ice(iflk) < 0.1_wp) t_ice_p(iflk)  = tpl_T_f
        ELSE IF (fr_ice(iflk) < 0.03_wp .AND. h_ice_p(iflk) > h_Ice_min_flk) THEN
          ! Ice fraction is too small, remove ice 
          h_ice_p(iflk)  = 0.0_wp
          ! Set the ice surface temperature to the fresh-water freezing point
          t_ice_p(iflk)  = tpl_T_f
          ! Set the mixed-layer temperature to the value slightly over the fresh-water freezing point
          t_wml_lk_p(iflk)  = tpl_T_f + 0.05_wp
          ! Adjust the bottom temperature
          t_bot_lk_p(iflk) = MAX(t_wml_lk_p(iflk), MIN(t_bot_lk_p(iflk), tpl_T_r))
          ! Compute the mean temperature of the water column
          t_mnw_lk_p(iflk) = t_wml_lk_p(iflk) - c_t_lk_p(iflk)                             &
            &              * MAX(0.0_wp, (1.0_wp - h_ml_lk_p(iflk)/depth_lk(iflk)))        &
            &              * (t_wml_lk_p(iflk) - t_bot_lk_p(iflk))
        END IF
       ! Set the snow variables
        h_snow_p(iflk)   = 0.0_wp
        t_snow_p(iflk)   = t_ice_p(iflk)

      ENDIF AssimIceFractionData

    END DO GridBoxesWithLakes 


!_tmp>
! COSMO stuff, unused.
! DR
! HOWEVER! Make sure that arrays for multi-layer snow and soil models
! and the like should/should not be initialized with the respective FLake variables.
! Also check/set the values of FLake variables at grid boxes without lakes.
!_nu
!_nu ! The outermost loop is over the time levels
!_nu DO nt = 1, nztlev
!_nu DO j = jstarts, jends
!_nu !CDIR NOMOVE
!_nu DO i = istarts, iends
!_nu   Lake_points: IF(depth_lk(i,j) > 0.0_wp) THEN
!_nu     ! At lake points, set "t_s", "t_g" and "t_so" to FLake values
!_nu     t_s(i,j,nt) = t_s(i,j,nx)
!_nu     t_g(i,j,nt) = t_s(i,j,nx)
!_nu     ! Save lake surface temperature in "t_so(:,:,0,:)"
!_nu     ! if the multi-layer soil model is used
!_nu     t_so(i,j,0,nt) = t_s(i,j,nx)
!_nu     ! Set FLake variables at all time levels
!_nu     ! ("nnew, "nnow", and for three time level scheme "nold")
!_nu     ! to their values at "nx" (nx=nnew for two time level scheme, see above).
!_nu     t_wml_lk(i,j,nt) = t_wml_lk(i,j,nx)
!_nu     t_mnw_lk(i,j,nt) = t_mnw_lk(i,j,nx)
!_nu     t_bot_lk(i,j,nt) = t_bot_lk(i,j,nx)
!_nu     c_t_lk  (i,j,nt) = c_t_lk  (i,j,nx)
!_nu     h_ml_lk (i,j,nt) = h_ml_lk (i,j,nx)
!_nu     h_ice   (i,j,nt) = h_ice   (i,j,nx)
!_nu     t_ice   (i,j,nt) = t_ice   (i,j,nx)
!_nu     ! As snow over lake ice is not considered explicitly,
!_nu     ! "h_snow" is zero
!_nu     h_snow  (i,j,nt) = h_snow  (i,j,nx)
!_nu     ! and "t_snow" is equal to the ice surface temperature
!_nu     t_snow  (i,j,nt) = t_snow  (i,j,nx)
!_nu   ELSE Lake_points
!_nu     ! Set FLake variables at non-lake (sea and land) points at their reference
!_nu     ! values at all time levels. Then, all FLake variables at land and sea points
!_nu     ! remain equal to their initial values during the entire COSMO run
!_nu     ! independent of permutation of time indices.
!_nu     t_wml_lk(i,j,nt) = tpl_T_r
!_nu     t_mnw_lk(i,j,nt) = tpl_T_r
!_nu     t_bot_lk(i,j,nt) = tpl_T_r
!_nu     c_t_lk  (i,j,nt) = C_T_min
!_nu     h_ml_lk (i,j,nt) = 0.0_wp
!_nu     ! If sea ice scheme is not used,
!_nu     ! "t_ice" and "h_ice" are only used if "llake" is set true.
!_nu     ! Initial values of "t_ice" and "h_ice" are set through the surface analysis.
!_nu     IF(.NOT.lseaice) THEN
!_nu       h_ice (i,j,nt) = h_ice (i,j,nx)
!_nu       t_ice (i,j,nt) = t_ice (i,j,nx)
!_nu     ENDIF
!_nu   ENDIF Lake_points
!_nu ! Variables are not part of the input (not read from GRIB files),
!_nu ! set them to their reference values at all points.
!_nu   t_b1_lk(i,j,nt) = tpl_T_r          ! Reference value, bottom-sediment module
!_nu                                      !     is switched off
!_nu   h_b1_lk(i,j,nt) = rflk_depth_bs_ref! Reference value, bottom-sediment module
!_nu                                      !     is switched off
!_nu END DO
!_nu END DO
!_nu END DO
!_nu 
!_nu ! If multi-layer snow model is used, set "t_snow_mult" at lake points
!_nu ! to the ice surface temperature at all vertical levels.
!_nu IF (lmulti_snow) THEN
!_nu   DO nt = 1, nztlev   ! DO loop over time levels
!_nu   DO ks = 0, ke_snow  ! DO loop over snow layers
!_nu     WHERE (depth_lk(:,:) > 0.0_wp) t_snow_mult(:,:,ks,nt) = t_snow(:,:,nx)
!_nu   END DO
!_nu   END DO
!_nu ENDIF
!_nu 
!_nu ! If multi-layer soil model is used, set "t_so(:,:,:,:)" at lake points
!_nu ! to the fresh-water temperature of maximum density at all levels except 0.
!_nu ! Notice that the climatological soil temperature,
!_nu ! i.e. t_so(:,:,ke_soil+1,:), is left intact.
!_nu DO nt = 1, nztlev   ! DO loop over time levels
!_nu DO ks = 1, ke_soil  ! DO loop over soil layers
!_nu   WHERE (depth_lk(:,:) > 0.0_wp) t_so(:,:,ks,nt) = tpl_T_r
!_nu END DO
!_nu END DO
!_nu  
!_tmp<

!_dbg>
!_nu    WRITE(*,*) 'FLake: WARM START initialization completed. '
!_dbg<

    !===============================================================================================
    !  END calculations
    !----------------------------------------------------------------------------------------------

  END SUBROUTINE flake_init

!===================================================================================================
!  End of the lake parameterization scheme FLake initialization
!---------------------------------------------------------------------------------------------------




  !>
  !! Coldstart for lake parameterization scheme Flake. 
  !!
  !! Coldstart for lake parameterization scheme Flake.
  !! Note that an estimate of the surface temperature is required to initialize the 
  !! mixed-layer temperature. Also note that no lake ice is assumed at the cold start.
  !! During cold start, prognostic Flake fields are only initialized at lake-points. 
  !! Non-lake points are filled with dummy values during the warmstart procedure. 
  !! 
  SUBROUTINE flake_coldinit (                                   & 
                     &  nflkgb, idx_lst_fp,                     &
                     &  depth_lk,                               &
                     &  tskin,                                  &
                     &  t_snow_lk_p, h_snow_lk_p,               & 
                     &  t_ice_p, h_ice_p,                       & 
                     &  t_mnw_lk_p, t_wml_lk_p, t_bot_lk_p,     &
                     &  c_t_lk_p, h_ml_lk_p,                    & 
                     &  t_b1_lk_p, h_b1_lk_p,                   &
                     &  t_g_lk_p,                               &
                     &  lacc                                    )


    IMPLICIT NONE


    ! Procedure arguments

    ! Array dimension(s)

    INTEGER, INTENT(IN) ::         &
                        &  nflkgb    !< array (vector) dimension,
                                     !< equal to the number of grid boxes within a block 
                                     !< where lakes are present
    INTEGER, DIMENSION(:), INTENT(IN) ::  &     !< index list for lake covered points
                        &  idx_lst_fp           !< within a block 

    ! FLake external parameters

    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  &
                        &  depth_lk             !< lake depth [m]


    ! Additional input fields needed for initialization

    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  &
                        &  tskin                !< proper estimate of mixed-layer temperature [K]
                                                !< e.g. skin temperature


    ! FLake prognostic variables
    ! (at the previous time step - "p", and the updated values - "n")

    REAL(KIND = wp)    , DIMENSION(:), INTENT(INOUT) ::  &
                        &  t_snow_lk_p      , & !< temperature of snow upper surface at previous time step [K]
                        &  h_snow_lk_p      , & !< snow thickness at previous time step [m]
                        &  t_ice_p          , & !< temperature of ice upper surface at previous time step [K]
                        &  h_ice_p          , & !< ice thickness at previous time step [m]
                        &  t_mnw_lk_p       , & !< mean temperature of the water column at previous time step [K]
                        &  t_wml_lk_p       , & !< mixed-layer temperature at previous time step [K] 
                        &  t_bot_lk_p       , & !< temperature at the water-bottom sediment interface 
                                                !< at previous time step [K] 
                        &  c_t_lk_p         , & !< shape factor with respect to the temperature profile 
                                                !< in lake thermocline at previous time step [-] 
                        &  h_ml_lk_p        , & !< thickness of the mixed-layer at previous time step [m] 
                        &  t_b1_lk_p        , & !< temperature at the bottom of the upper layer
                                                !< of the sediments at previous time step [K]  
                        &  h_b1_lk_p        , & !< thickness of the upper layer of bottom sediments 
                                                !< at previous time step [m] 
                        &  t_g_lk_p             !< lake surface temperature at previous time [K]

    LOGICAL, OPTIONAL,   INTENT(IN)   :: lacc ! If true, use openacc



    INTEGER ::      &
            &  ic , &  !< DO loop index
            &  jc      !< icon grid cell index

    LOGICAL :: lzacc ! non-optional version of lacc


    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    ! Loop over grid boxes with lakes
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(jc)
!$NEC ivdep
    DO ic=1, nflkgb
      ! Take index from the lake index list
      jc = idx_lst_fp(ic)

      ! Mixed-layer depth is set equal to 10 m or to the lake depth, whichever is smaller 
      h_ml_lk_p(jc) = MIN(depth_lk(jc), 10._wp)


      ! Shape factor is set to its minimum value 
      c_t_lk_p(jc) = C_T_min


      ! Mixed-layer temperature is set to the surface temperature
      ! that should be provided for the cold-start initialization 
      ! (e.g. skin temperature can be taken as an estimate of t_wml_lk)
      ! Clearly, the mixed-layer temperature cannot be less than the fresh-water freezing point 
      t_wml_lk_p(jc) = MAX(tskin(jc), tpl_T_f)


      ! Bottom temperature is set equal 
      ! to the mixed-layer temperarture if the water column is mixed
      ! or to the temperarture of maximum density of fresh water  
      ! if the mixed-layer depth is less than the depth to the bottom
      t_bot_lk_p(jc) = MERGE( t_wml_lk_p(jc), tpl_T_r, h_ml_lk_p(jc)>= depth_lk(jc) )


      ! Mean temperature of the water column is computed from the formula
      ! t_mnw_lk = t_wml_lk - c_t_lk * MAX(0., 1.-h_ml_lk/depth_lk) * (t_wml_lk-t_bot_lk) 
      t_mnw_lk_p(jc) = t_wml_lk_p(jc) - c_t_lk_p(jc)                  &
        &            * MAX(0._wp, (1._wp - h_ml_lk_p(jc)/depth_lk(jc))) &
        &            * (t_wml_lk_p(jc) - t_bot_lk_p(jc))


      ! No ice is assumed at cold-start initilization
      t_ice_p(jc) = tpl_T_f      ! Fresh-water freezing point
      h_ice_p(jc) = 0._wp    ! Zero ice thickness


      ! Snow over lake ice is not considered explicitly
      t_snow_lk_p(jc) = t_ice_p(jc)                                   ! Ice surface temperature
      h_snow_lk_p(jc) = 0._wp                                     ! Zero snow thickness


      ! Bottom sediment module of FLake is switched off
      t_b1_lk_p(jc) = tpl_T_r                                         ! Reference value
      h_b1_lk_p(jc) = rflk_depth_bs_ref                               ! Reference value


      ! Set lake-surface temperature equal to the mixed-layer temperature
      t_g_lk_p(jc) = t_wml_lk_p(jc)

    END DO  ! End of loop over grid boxes with lakes
    !$ACC END PARALLEL


  END SUBROUTINE flake_coldinit




!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

!===================================================================================================
!  FLake interface. Advances FLake variables one time step forward.
!---------------------------------------------------------------------------------------------------

  !> 
  !! FLake interface is a communication routine between "flake_driver" and ICON
  !! (or any other model that uses FLake).
  !! It assigns the FLake variables at the previous time step 
  !! to their input values given by ICON,
  !! assigns several quantities to local FLake variables 
  !! (e.g. surface fluxes of heat and momentum that are computed outside FLake), 
  !! computes a few quantities required by FLake (e.g. surface friction velocity),
  !! computes various quantities related to the volumetric (vertically distributed) 
  !! solar radiation heating (see 'CALL flake_radflux' below),
  !! calls "flake_driver",
  !! and returns the updated FLake variables to ICON.
  !! Optionally, the time tendencies of FLake variables can be returned to ICON.
  !! The "flake_interface" does not contain any FLake physics. 
  !! It only serves as a convenient means to organize a call of "flake_driver".
  !! Computations, including a call of "flake_driver", are executed in a DO loop 
  !! over all ICON grid boxes that contain lakes 
  !! (lake fraction in excess of a minimum threshold value).
  !! Notice that inlining of some routines is necessary to ensure 
  !! DO loop vectorization (cf. COSMO).
  !!

  SUBROUTINE flake_interface (                                       &
                     &  dtime,                                       &
                     &  nflkgb,                                      &
                     &  coriolispar,                                 &
                     &  depth_lk,                                    &
                     &  fetch_lk, dp_bs_lk, t_bs_lk, gamso_lk,       &
                     &  qmom, qsen, qlat, qlwrnet, qsolnet,          &
                     &  t_snow_p, h_snow_p,                          & 
                     &  t_ice_p, h_ice_p,                            &
                     &  t_mnw_lk_p, t_wml_lk_p, t_bot_lk_p,          &
                     &  c_t_lk_p, h_ml_lk_p,                         & 
                     &  t_b1_lk_p, h_b1_lk_p,                        &           
                     &  t_scf_lk_p,                                  & 
                     &  t_snow_n, h_snow_n,                          & 
                     &  t_ice_n, h_ice_n,                            &
                     &  t_mnw_lk_n, t_wml_lk_n, t_bot_lk_n,          &
                     &  c_t_lk_n, h_ml_lk_n,                         &
                     &  t_b1_lk_n, h_b1_lk_n,                        &
                     &  t_scf_lk_n,                                  &
                     &  opt_dtsnowdt, opt_dhsnowdt,                  &
                     &  opt_dticedt, opt_dhicedt,                    &
                     &  opt_dtmnwlkdt, opt_dtwmllkdt, opt_dtbotlkdt, & 
                     &  opt_dctlkdt, opt_dhmllkdt,                   & 
                     &  opt_dtb1lkdt, opt_dhb1lkdt,                  & 
                     &  opt_dtsfclkdt,                               &
                     &  lacc)

    IMPLICIT NONE

    ! Procedure arguments

    ! Variables

    REAL(KIND = wp)    , INTENT(IN) ::  &
                        &  dtime          !< model time step [s]       
                     
    ! Array dimension(s)

    INTEGER, INTENT(IN) ::           &
                        &  nflkgb      !< array (vector) dimension,
                                       !< equal to the number of grid boxes within a block
                                       !< where lakes are present 

    ! Arrays

    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  & 
                        &  coriolispar          !< Coriolis parameter [s^{-1}]
    
    ! FLake external parameters            
                                           
    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  & 
                        &  depth_lk             !< lake depth [m]
    
    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  &  
                        &  fetch_lk         , & !< wind fetch over lake [m]
                        &  dp_bs_lk         , & !< depth of the thermally active layer
                                                !< of bottom sediments [m]
                        &  t_bs_lk          , & !< climatological temperature at the bottom of
                                                !< the thermally active layer of sediments [K]
                        &  gamso_lk             !< attenuation coefficient of the lake water 
                                                !< with respect to solar radiation [m^{-1}]
                                            
    ! Surface fluxes of momentum and heat (positive downward)

    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  & 
                        &  qmom             , & !< momentum flux at the surface [N m^{-2}]
                        &  qsen             , & !< sensible heat flux at the surface [W m^{-2}]
                        &  qlat             , & !< latent heat flux at the surface [W m^{-2}]
                        &  qlwrnet          , & !< net long-wave radiation flux at the surface [W m^{-2}]
                        &  qsolnet              !< net solar radiation flux at the surface [W m^{-2}]


    ! FLake prognostic variables
    ! (at the previous time step - "p", and the updated values - "n")

    REAL(KIND = wp)    , DIMENSION(:), INTENT(IN) ::  &
                        &  t_snow_p         , & !< temperature of snow upper surface at previous time step [K]
                        &  h_snow_p         , & !< snow thickness at previous time step [m]
                        &  t_ice_p          , & !< temperature of ice upper surface at previous time step [K]
                        &  h_ice_p          , & !< ice thickness at previous time step [m]
                        &  t_mnw_lk_p       , & !< mean temperature of the water column at previous time step [K]
                        &  t_wml_lk_p       , & !< mixed-layer temperature at previous time step [K]
                        &  t_bot_lk_p       , & !< temperature at the water-bottom sediment interface
                                                !< at previous time step [K]
                        &  c_t_lk_p         , & !< shape factor with respect to the temperature profile
                                                !< in lake thermocline at previous time step [-]
                        &  h_ml_lk_p        , & !< thickness of the mixed-layer at previous time step [m]
                        &  t_b1_lk_p        , & !< temperature at the bottom of the upper layer
                                                !< of the sediments at previous time step [K]
                        &  h_b1_lk_p        , & !< thickness of the upper layer of bottom sediments
                                                !< at previous time step [m]
                        &  t_scf_lk_p           !< lake surface temperature at previous time step [K]
                                                !< (i.e. the temperature at the air-water, air-ice
                                                !< or air-snow interface)

    REAL(KIND = wp)    , DIMENSION(:), INTENT(OUT) ::    &
                        &  t_snow_n         , & !< temperature of snow upper surface at new time step [K]
                        &  h_snow_n         , & !< snow thickness at new time step [m]
                        &  t_ice_n          , & !< temperature of ice upper surface at new time step [K]
                        &  h_ice_n          , & !< ice thickness at new time step [m]
                        &  t_mnw_lk_n       , & !< mean temperature of the water column at new time step [K]
                        &  t_wml_lk_n       , & !< mixed-layer temperature at new time step [K]
                        &  t_bot_lk_n       , & !< temperature at the water-bottom sediment interface
                                                !< at new time step [K]
                        &  c_t_lk_n         , & !< shape factor with respect to the temperature profile
                                                !< in lake thermocline at new time step [-]
                        &  h_ml_lk_n        , & !< thickness of the mixed-layer at new time step [m]
                        &  t_b1_lk_n        , & !< temperature at the bottom of the upper layer
                                                !< of the sediments at new time step [K]
                        &  h_b1_lk_n        , & !< thickness of the upper layer of bottom sediments
                                                !< at new time step [m]
                        &  t_scf_lk_n           !< lake surface temperature at new time step [K]
                                                !< (i.e. the temperature at the air-water, air-ice
                                                !< or air-snow interface)

    ! Tendencies of FLake variables (optional)    

    REAL(KIND = wp)    , DIMENSION(:), INTENT(OUT), OPTIONAL ::  & 
                        &  opt_dtsnowdt     , & !< time tendency of 
                                                !< snow surface temperature [K s^{-1}]
                        &  opt_dhsnowdt     , & !< time tendency of 
                                                !< snow thickness [m s^{-1}]
                        &  opt_dticedt      , & !< time tendency of 
                                                !< ice surface temperature [K s^{-1}]
                        &  opt_dhicedt      , & !< time tendency of 
                                                !< ice thickness [m s^{-1}]
                        &  opt_dtmnwlkdt    , & !< time tendency of 
                                                !< mean temperature of the water column [K s^{-1}]
                        &  opt_dtwmllkdt    , & !< time tendency of 
                                                !< mixed-layer temperature [K s^{-1}]
                        &  opt_dtbotlkdt    , & !< time tendency of 
                                                !< temperature at the water-bottom sediment 
                                                !< interface [K s^{-1}]
                        &  opt_dctlkdt      , & !< time tendency of 
                                                !< shape factor with respect to the temperature 
                                                !< profile in lake thermocline [s^{-1}]
                        &  opt_dhmllkdt     , & !< time tendency of 
                                                !< the mixed-layer thickness [m s^{-1}]
                        &  opt_dtb1lkdt     , & !< time tendency of 
                                                !< temperature at the bottom of the upper layer 
                                                !< of sediments [K s^{-1}]
                        &  opt_dhb1lkdt     , & !< time tendency of 
                                                !< thickness of the upper layer 
                                                !< of sediments [m s^{-1}]
                        &  opt_dtsfclkdt        !< time tendency of 
                                                !< lake surface temperature [K s^{-1}]

    LOGICAL, INTENT(IN), OPTIONAL :: &
                        &  lacc                 !< openACC flag

    ! Local variables and arrays

    INTEGER ::                   &
      &  iflk                      !< DO loop index

    INTEGER ::                   &
      &  izdebug                   !< debugging key 

    REAL (KIND = wp)     ::      &
      &  del_time              , & !< time step used by FLake [s]
      &  depth_w               , & !< lake depth [m]
      &  fetch                 , & !< typical wind fetch [m]
      &  depth_bs              , & !< depth of the thermally active layer 
                                   !< of the bottom sediments [m]
      &  T_bs                  , & !< temperature at the bottom 
                                   !< of the thermally active layer of sediments [K]
      &  par_Coriolis          , & !< Coriolis parameter [s^{-1}]
      &  T_sfc_p               , & !< surface temperature at the previous time step [K]  
      &  T_sfc_n               , & !< updated surface temperature [K]  
      &  r_dtime                   !< reciprocal of the time step [s^{-1}]

    REAL (KIND = wp)     ::      &
      &  albedo_water          , & !< water surface albedo with respect to the solar radiation [-]
      &  albedo_ice            , & !< ice surface albedo with respect to the solar radiation [-]
      &  albedo_snow               !< snow surface albedo with respect to the solar radiation [-]

    LOGICAL ::                   &
      &  lzacc                     !< openACC flag

!===================================================================================================

!_cdm>
! COSMO stuff. Currently unused within ICON.
!_cdm<

!_nu CHARACTER (LEN=80)   :: yzerrmsg
!_nu
!_nu !  Tracer pointers
!_nu
!_nu REAL (KIND=wp)    , POINTER :: &
!_nu   qv  (:,:,:) => NULL()      ! QV at tlev=nx
!_nu
!_nu CHARACTER(LEN=25) :: yzroutine = 'flake_interface'
!_nu
!_nu #ifdef MESSY
!_nu LOGICAL, SAVE :: l_init_pointers = .TRUE.
!_nu #endif

!===================================================================================================

    ! Tendencies of FLake variables 

    REAL(KIND = wp)    , DIMENSION(nflkgb) ::  & 
                        &  dtsnowdt          , & !< time tendency of 
                                                 !< snow surface temperature [K s^{-1}]
                        &  dhsnowdt          , & !< time tendency of 
                                                 !< snow thickness [m s^{-1}]
                        &  dticedt           , & !< time tendency of 
                                                 !< ice surface temperature [K s^{-1}]
                        &  dhicedt           , & !< time tendency of 
                                                 !< ice thickness [m s^{-1}]
                        &  dtmnwlkdt         , & !< time tendency of 
                                                 !< mean temperature of the water column [K s^{-1}]
                        &  dtwmllkdt         , & !< time tendency of 
                                                 !< mixed-layer temperature [K s^{-1}]
                        &  dtbotlkdt         , & !< time tendency of 
                                                 !< temperature at the water-bottom sediment 
                                                 !< interface [K s^{-1}]
                        &  dctlkdt           , & !< time tendency of 
                                                 !< shape factor with respect to the temperature 
                                                 !< profile in lake thermocline [s^{-1}]
                        &  dhmllkdt          , & !< time tendency of 
                                                 !< the mixed-layer thickness [m s^{-1}]
                        &  dtb1lkdt          , & !< time tendency of 
                                                 !< temperature at the bottom of the upper layer 
                                                 !< of sediments [K s^{-1}]
                        &  dhb1lkdt          , & !< time tendency of 
                                                 !< thickness of the upper layer 
                                                 !< of sediments [m s^{-1}]
                        &  dtsfclkdt             !< time tendency of 
                                                 !< lake surface temperature [K s^{-1}]


  ! ** The following local variables used to be module variables, which had to be shifted here in
  !    order to allow OpenMP parallelization ** (GZ, 2013-11-29)
  !
  ! Temperatures at the previous time step ("p"), and the updated temperatures ("n")
  REAL (KIND = wp)     ::            &
    &  T_mnw_p_flk,  T_mnw_n_flk   , & !< Mean temperature of the water column [K]
    &  T_snow_p_flk, T_snow_n_flk  , & !< Temperature at the air-snow interface [K]
    &  T_ice_p_flk,  T_ice_n_flk   , & !< Temperature at the snow-ice or air-ice interface [K]
    &  T_wML_p_flk,  T_wML_n_flk   , & !< Mixed-layer temperature [K]
    &  T_bot_p_flk,  T_bot_n_flk   , & !< Temperature at the water-bottom sediment interface [K]
    &  T_B1_p_flk,   T_B1_n_flk        !< Temperature at the bottom of the upper layer     
                                       !< of the sediments [K]

  ! Thickness of various layers at the previous time step ("p") and the updated values ("n")
  REAL (KIND = wp)     ::            &
    &  h_snow_p_flk, h_snow_n_flk  , & !< Snow thickness [m]
    &  h_ice_p_flk,  h_ice_n_flk   , & !< Ice thickness [m]
    &  h_ML_p_flk,   h_ML_n_flk    , & !< Thickness of the mixed-layer [m]
    &  H_B1_p_flk,   H_B1_n_flk        !< Thickness of the upper layer of bottom sediments [m]

  ! The shape factor(s) at the previous time step ("p") and the updated value(s) ("n")
  REAL (KIND = wp)     ::            &
    &  C_T_p_flk, C_T_n_flk        , & !< Shape factor (thermocline) [-]
    &  C_TT_flk                    , & !< Dimensionless parameter (thermocline) [-]
    &  C_Q_flk                     , & !< Shape factor with respect to heat flux (thermocline) [-]
    &  C_I_flk                     , & !< Shape factor (ice) [-]
    &  C_S_flk                         !< Shape factor (snow) [-]

  ! Derivatives of the shape functions
  REAL (KIND = wp)     ::            &
    &  Phi_T_pr0_flk               , & !< d\Phi_T(0)/d\zeta   (thermocline) [-]
    &  Phi_I_pr0_flk               , & !< d\Phi_I(0)/d\zeta_I (ice) [-]
    &  Phi_I_pr1_flk               , & !< d\Phi_I(1)/d\zeta_I (ice) [-]
    &  Phi_S_pr0_flk                   !< d\Phi_S(0)/d\zeta_S (snow) [-]

  ! Heat and radiation fluxes
  REAL (KIND = wp)     ::            &
    &  Q_snow_flk                  , & !< Heat flux through the air-snow interface [W m^{-2}]
    &  Q_ice_flk                   , & !< Heat flux through the snow-ice or air-ice 
                                       !< interface [W m^{-2}]
    &  Q_w_flk                     , & !< Heat flux through the ice-water or air-water 
                                       !< interface [W m^{-2}]
    &  Q_bot_flk                   , & !< Heat flux through the water-bottom sediment 
                                       !< interface [W m^{-2}]
    &  I_atm_flk                   , & !< Radiation flux at the lower boundary of 
                                       !< the atmosphere [W m^{-2}], i.e. the incident radiation flux 
                                       !< with no regard for the surface albedo 
    &  I_snow_flk                  , & !< Radiation flux through the air-snow interface [W m^{-2}]
    &  I_ice_flk                   , & !< Radiation flux through the snow-ice or air-ice 
                                       !< interface [W m^{-2}]
    &  I_w_flk                     , & !< Radiation flux through the ice-water or air-water 
                                       !< interface [W m^{-2}]
    &  I_h_flk                     , & !< Radiation flux through the mixed-layer-thermocline 
                                       !< interface [W m^{-2}]
    &  I_bot_flk                   , & !< Radiation flux through the water-bottom sediment 
                                       !< interface [W m^{-2}]
    &  I_intm_0_h_flk              , & !< Mean radiation flux over the mixed layer [W m^{-2}]
    &  I_intm_h_D_flk              , & !< Mean radiation flux over the thermocline [W m^{-2}]
    &  Q_star_flk                      !< A generalized heat flux scale [W m^{-2}]

  ! Velocity scales
  REAL (KIND = wp)     ::            &
    &  u_star_w_flk                , & !< Friction velocity in the surface layer 
                                       !< of lake water [m s^{-1}]
    &  w_star_sfc_flk                  !< Convective velocity scale based on
                                       !< a generalized heat flux scale [m s^{-1}]

  ! The rate of snow accumulation
  REAL (KIND = wp)     ::            &
    &  dMsnowdt_flk                    !< The rate of snow accumulation [kg m^{-2} s^{-1}]


    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)


!_cdm>
! A debugging key "izdebug" is actually not required,
! see "idbg" in subroutine "flake_driver".
! However, it is kept in the call of "flake_driver" and as a "flake_driver" argument 
! for (better) compatibility between the ICON and COSMO codes.
!_cdm<

    ! set debugging key  
    izdebug = 0

    ! set model time step
    del_time = dtime

    ! reciprocal of the model time step
    r_dtime = 1._wp/dtime

!===================================================================================================

!_cdm>
! COSMO stuff.
! Currently unused within ICON.
!_cdm<

!_nu !------------------------------------------------------------------------------
!_nu ! Retrieve pointer to required tracers
!_nu !------------------------------------------------------------------------------
!_nu 
!_nu #ifndef MESSY
!_nu CALL trcr_get(izerror, 'QV', ptr_tlev = nx, ptr = qv)
!_nu IF (izerror /= 0) THEN
!_nu   yzerrmsg = trcr_errorstr(izerror)
!_nu   CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
!_nu ENDIF
!_nu #else
!_nu   IF (.NOT.ASSOCIATED(qv)) THEN
!_nu      CALL get_tracer(izerror, GPTRSTR, 'QV',  pxtm1=qv)
!_nu      CALL tracer_halt(yzroutine,izerror)
!_nu   ENDIF
!_nu #endif

!===================================================================================================

    
    !-----------------------------------------------------------------------------------------------
    !  Set optical characteristics of the lake water, lake ice and snow
    !-----------------------------------------------------------------------------------------------

    ! Use default values
    opticpar_water = opticpar_water_ref    ! Reference value, 
                                           ! 95% of radiation is absorbed within the upper 1 m
    opticpar_ice   = opticpar_ice_opaque   ! Opaque ice
    opticpar_snow  = opticpar_snow_opaque  ! Opaque snow
    !-----------------------------------------------------------------------------------------------
    !  Set/compute albedos of the lake water, ice and snow
    !-----------------------------------------------------------------------------------------------

    ! Use default value for water  
!_nu    albedo_water = albedo_water_ref
    ! Use interpolation formula for ice
!_nu    albedo_ice   = EXP(-c_albice_MR*(tpl_T_f-T_sfc_p)/tpl_T_f)
!_nu    albedo_ice   = albedo_whiteice_ref*(1._wp-albedo_ice) + albedo_blueice_ref*albedo_ice

    ! The net solar radiation flux at the surface is computed outside FLake routines 
    ! and should already account for the surface albedo 
    ! (make sure that calculations are performed consistently in the respective ICON routines).
    ! Then, the surface albedo of the lake water, ice and snow should be set to zero here
    ! (see subroutine "flake_radflux").
    albedo_water = 0._wp
    albedo_ice   = 0._wp
    albedo_snow  = albedo_ice     ! snow is not considered explicitly

    !$ACC DATA &
    !$ACC   CREATE(dctlkdt, dhb1lkdt, dhicedt, dhmllkdt, dhsnowdt, dtb1lkdt, dtbotlkdt) &
    !$ACC   CREATE(dticedt, dtmnwlkdt, dtsfclkdt, dtsnowdt, dtwmllkdt) &
    !$ACC   COPYIN(opticpar_water, opticpar_ice, opticpar_snow) ASYNC(1) IF(lzacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR &
    !$ACC   PRIVATE(c_i_flk, c_tt_flk, c_q_flk, c_t_n_flk, c_t_p_flk, depth_w, fetch, h_b1_n_flk, h_b1_p_flk) &
    !$ACC   PRIVATE(h_ice_n_flk, h_ice_p_flk, h_ml_n_flk, h_ml_p_flk, h_snow_n_flk, h_snow_p_flk, i_atm_flk) &
    !$ACC   PRIVATE(i_bot_flk, i_h_flk, i_ice_flk, i_intm_0_h_flk, i_intm_h_d_flk, i_snow_flk, i_w_flk) &
    !$ACC   PRIVATE(par_coriolis, phi_i_pr0_flk, phi_i_pr1_flk, phi_t_pr0_flk, q_bot_flk, q_ice_flk) &
    !$ACC   PRIVATE(q_snow_flk, q_star_flk, q_w_flk, t_b1_n_flk, t_b1_p_flk, t_bot_n_flk, t_bot_p_flk, t_bs) &
    !$ACC   PRIVATE(t_ice_n_flk, t_ice_p_flk, t_mnw_n_flk, t_mnw_p_flk, t_sfc_n, t_sfc_p, t_snow_n_flk) &
    !$ACC   PRIVATE(t_snow_p_flk, t_wml_n_flk, t_wml_p_flk, u_star_w_flk, w_star_sfc_flk, depth_bs, dmsnowdt_flk)
    
    !-----------------------------------------------------------------------------------------------
    !  DO loop over grid boxes with lakes
    !-----------------------------------------------------------------------------------------------

    GridBoxesWithLakes: DO iflk=1, nflkgb

      !---------------------------------------------------------------------------------------------
      !  Set (external) parameters 
      !---------------------------------------------------------------------------------------------

      depth_w      = depth_lk   (iflk)
      fetch        = fetch_lk   (iflk)      
      depth_bs     = dp_bs_lk   (iflk)
      T_bs         = t_bs_lk    (iflk)
      par_Coriolis = coriolispar(iflk)

      !---------------------------------------------------------------------------------------------
      !  Set initial values
      !---------------------------------------------------------------------------------------------

      T_snow_p_flk = t_snow_p  (iflk)
      h_snow_p_flk = h_snow_p  (iflk)
      T_ice_p_flk  = t_ice_p   (iflk)  
      h_ice_p_flk  = h_ice_p   (iflk)
      T_mnw_p_flk  = t_mnw_lk_p(iflk) 
      T_wML_p_flk  = t_wml_lk_p(iflk)
      T_bot_p_flk  = t_bot_lk_p(iflk)
      C_T_p_flk    = c_t_lk_p  (iflk)
      h_ML_p_flk   = h_ml_lk_p (iflk)
      H_B1_p_flk   = h_b1_lk_p (iflk)
      T_B1_p_flk   = t_b1_lk_p (iflk)
      T_sfc_p      = t_scf_lk_p(iflk)    

      !---------------------------------------------------------------------------------------------
      !  Set the rate of snow accumulation
      !---------------------------------------------------------------------------------------------

      dMsnowdt_flk = 0._wp      ! snow is not considered explicitly

      !---------------------------------------------------------------------------------------------
      !  Surface solar radiation fluxes (positive downward) and related quantities 
      !  required to account for volumetric solar heating
      !---------------------------------------------------------------------------------------------

      I_atm_flk = qsolnet(iflk)  ! net surface solar radiation flux from ICON 
                                 ! with due regard fr the surface albedo

      CALL flake_radflux (  depth_w, albedo_water, albedo_ice, albedo_snow,         &
                            opticpar_water, opticpar_ice, opticpar_snow,            &
                            h_snow_p_flk, h_ice_p_flk, h_ML_p_flk, I_snow_flk,      &
                            I_bot_flk, I_ice_flk, I_w_flk, I_h_flk, I_intm_0_h_flk, &
                            I_intm_h_D_flk, I_atm_flk)

      !---------------------------------------------------------------------------------------------
      !  Compute surface heat flux (positive downward)
      !---------------------------------------------------------------------------------------------

      Q_w_flk = qlwrnet(iflk)    ! net long-wave radiation flux from ICON, 
                                 ! "qlwrnet" already includes the upward radiation from the surface
      Q_w_flk = Q_w_flk + qsen(iflk) + qlat(iflk)
                                 ! add sensible and latent heat fluxes (positive downward),
                                 ! make sure that the latent heat of fusion is used when computing "qlat"
                                 ! over ice or snow

      ! Note that the FLake flux calculation routines are not used at present.
      ! Fluxes of momentum and of sensible and latent heat are computed
      ! using transfer coefficients from the ICON turbulence routines (c/o MR). 

      !---------------------------------------------------------------------------------------------
      !  Compute surface friction velocity in the water 
      !---------------------------------------------------------------------------------------------

      u_star_w_flk = SQRT(ABS(qmom(iflk))/tpl_rho_w_r)

      !---------------------------------------------------------------------------------------------
      !  Compute heat fluxes Q_snow_flk, Q_ice_flk and Q_w_flk
      !---------------------------------------------------------------------------------------------

      IF(h_ice_p_flk .GE. h_Ice_min_flk) THEN        ! ice exists
        IF(h_snow_p_flk .GE. h_Snow_min_flk) THEN    ! there is snow above lake ice
          Q_snow_flk = Q_w_flk
          Q_ice_flk  = 0._wp
          Q_w_flk    = 0._wp
        ELSE                                         ! no snow above lake ice
          Q_snow_flk = 0._wp
          Q_ice_flk  = Q_w_flk
          Q_w_flk    = 0._wp
        END IF
      ELSE                                           ! no ice-snow cover
          Q_snow_flk = 0._wp
          Q_ice_flk  = 0._wp
      END IF

      !---------------------------------------------------------------------------------------------
      !  Advance FLake variables one time step forward 
      !---------------------------------------------------------------------------------------------
      CALL flake_driver ( depth_w, depth_bs, T_bs, par_Coriolis,         &
                          opticpar_water%extincoef_optic(1),             &
                          del_time, T_sfc_p, T_sfc_n, izdebug,           &    
                          T_snow_n_flk, T_ice_n_flk, T_mnw_n_flk,        &
                          T_bot_n_flk, T_wML_n_flk,                      &
                          T_B1_n_flk, C_T_n_flk, h_snow_n_flk,           &
                          h_ice_n_flk, h_ML_n_flk, H_B1_n_flk,           &
                          Q_bot_flk, Q_star_flk, w_star_sfc_flk,         &
                          Q_w_flk, Q_ice_flk, Q_snow_flk, T_snow_p_flk,  &
                          T_ice_p_flk, T_mnw_p_flk, T_bot_p_flk,         &
                          T_wML_p_flk, T_B1_p_flk, C_T_p_flk,            &
                          h_snow_p_flk, h_ice_p_flk, h_ML_p_flk,         &
                          H_B1_p_flk, I_bot_flk, I_snow_flk, I_ice_flk,  &
                          I_w_flk, I_h_flk, I_intm_h_D_flk,              &
                          I_intm_0_h_flk, u_star_w_flk, dMsnowdt_flk,    &
                          Phi_T_pr0_flk, Phi_I_pr1_flk, Phi_I_pr0_flk,   &
                          C_I_flk, C_TT_flk, C_Q_flk )

      !---------------------------------------------------------------------------------------------
      !  Set output values
      !---------------------------------------------------------------------------------------------

      t_snow_n  (iflk) = T_snow_n_flk
      h_snow_n  (iflk) = h_snow_n_flk   
      t_ice_n   (iflk) = T_ice_n_flk      
      h_ice_n   (iflk) = h_ice_n_flk    
      t_mnw_lk_n(iflk) = T_mnw_n_flk     
      t_wml_lk_n(iflk) = T_wML_n_flk    
      t_bot_lk_n(iflk) = T_bot_n_flk   
      c_t_lk_n  (iflk) = C_T_n_flk     
      h_ml_lk_n (iflk) = h_ML_n_flk     
      t_b1_lk_n (iflk) = T_B1_n_flk    
      h_b1_lk_n (iflk) = H_B1_n_flk    
      t_scf_lk_n(iflk) = T_sfc_n       

!===================================================================================================

!_cdm>
! Stuff from COSMO. Currently unused within ICON.
! DR
! Check if arrays for multi-layer snow and soil models
! and the like should/should not be filled with the respective FLake variables.
!_cdm<

!_nu ! Weighted surface temperature of the COSMO grid point
!_nu ! is equal to the lake surface temperature
!_nu     t_g     (i,j,nnew) = t_s (i,j,nnew)
!_nu 
!_nu ! Save lake surface temperature in t_so(:,:,0,:)
!_nu ! in case the multi-layer soil model is used.
!_nu     t_so(i,j,0,nnew) = t_s (i,j,nnew)
!_nu 
!_nu IF (lmulti_snow) THEN
!_nu   ! If multi-layer snow model is used, save updated snow temperature
!_nu   ! also in "t_snow_mult" at all vertical levels
!_nu   DO ksnow = 0, ke_snow  ! DO loop over snow layers
!_nu     WHERE (depth_lk(:,:) > 0.0_wp) t_snow_mult(:,:,ksnow,nnew) = t_snow(:,:,nnew)
!_nu   END DO
!_nu ENDIF

!===================================================================================================

      !---------------------------------------------------------------------------------------------
      !  Compute time tendencies of FLake variables (for eventual use outside FLake program units)
      !---------------------------------------------------------------------------------------------

      dtsnowdt (iflk) = ( t_snow_n  (iflk) - t_snow_p  (iflk) )*r_dtime
      dhsnowdt (iflk) = ( h_snow_n  (iflk) - h_snow_p  (iflk) )*r_dtime
      dticedt  (iflk) = ( t_ice_n   (iflk) - t_ice_p   (iflk) )*r_dtime
      dhicedt  (iflk) = ( h_ice_n   (iflk) - h_ice_p   (iflk) )*r_dtime
      dtmnwlkdt(iflk) = ( t_mnw_lk_n(iflk) - t_mnw_lk_p(iflk) )*r_dtime 
      dtwmllkdt(iflk) = ( t_wml_lk_n(iflk) - t_wml_lk_p(iflk) )*r_dtime 
      dtbotlkdt(iflk) = ( t_bot_lk_n(iflk) - t_bot_lk_p(iflk) )*r_dtime 
      dctlkdt  (iflk) = ( c_t_lk_n  (iflk) - c_t_lk_p  (iflk) )*r_dtime 
      dhmllkdt (iflk) = ( h_ml_lk_n (iflk) - h_ml_lk_p (iflk) )*r_dtime 
      dtb1lkdt (iflk) = ( t_b1_lk_n (iflk) - t_b1_lk_p (iflk) )*r_dtime 
      dhb1lkdt (iflk) = ( h_b1_lk_n (iflk) - h_b1_lk_p (iflk) )*r_dtime 
      dtsfclkdt(iflk) = ( t_scf_lk_n(iflk) - t_scf_lk_p(iflk) )*r_dtime

    !-----------------------------------------------------------------------------------------------
    !  End of DO loop over grid boxes with lakes
    !-----------------------------------------------------------------------------------------------

    END DO GridBoxesWithLakes
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------------------------------
    !  Store time tendencies of FLake variables (optional)
    !-----------------------------------------------------------------------------------------------

    IF (PRESENT(opt_dtsnowdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dtsnowdt(1:nflkgb) = dtsnowdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dtsnowdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dtsnowdt(nflkgb+1:)= 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dhsnowdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dhsnowdt(1:nflkgb) = dhsnowdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dhsnowdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dhsnowdt(nflkgb+1:)= 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dticedt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dticedt(1:nflkgb) = dticedt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dticedt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dticedt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dhicedt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dhicedt(1:nflkgb) = dhicedt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dhicedt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dhicedt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dtmnwlkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dtmnwlkdt(1:nflkgb) = dtmnwlkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dtmnwlkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dtmnwlkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dtwmllkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dtwmllkdt(1:nflkgb) = dtwmllkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dtwmllkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dtwmllkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dtbotlkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dtbotlkdt(1:nflkgb) = dtbotlkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dtbotlkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dtbotlkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dctlkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dctlkdt(1:nflkgb) = dctlkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dctlkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dctlkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dhmllkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dhmllkdt(1:nflkgb) = dhmllkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dhmllkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dhmllkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dtb1lkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dtb1lkdt(1:nflkgb) = dtb1lkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dtb1lkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dtb1lkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dhb1lkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dhb1lkdt(1:nflkgb) = dhb1lkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dhb1lkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dhb1lkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF
    IF (PRESENT(opt_dtsfclkdt)) THEN
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      opt_dtsfclkdt(1:nflkgb) = dtsfclkdt(1:nflkgb)
      !$ACC END KERNELS
      IF (nflkgb < SIZE(opt_dtsfclkdt)) THEN
        !$ACC KERNELS ASYNC(1) IF(lzacc)
        opt_dtsfclkdt(nflkgb+1:) = 0._wp
        !$ACC END KERNELS
      ENDIF
    ENDIF

  

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

    IF (.NOT. using_cuda_graph) THEN
      !$ACC WAIT(1)
    END IF
    !$ACC END DATA

END SUBROUTINE flake_interface

!===================================================================================================

!_cdm>
! All program units below (SUBROUTINE, FUNCTION) are identical within ICON and COSMO  
! (to within a few comment lines).
! Since it is these routines that contain the FLake physics,
! the actual "fresh-water lake parameterization scheme" 
! remains the same in ICON and COSMO.
! The subroutines "flake_init" and "flake_interface" are different in ICON and COSMO.
! However, these two subroutines are just communication routines 
! between FLake and a host atmospheric model 
! and are therefore host-model dependent.
!_cdm<

!===================================================================================================

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

SUBROUTINE flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow,       & 
                           opticpar_water, opticpar_ice, opticpar_snow,          &
                           h_snow_p_flk,h_ice_p_flk,h_ML_p_flk,I_snow_flk,       &
                           I_bot_flk, I_ice_flk, I_w_flk, I_h_flk,               &
                           I_intm_0_h_flk, I_intm_h_D_flk, I_atm_flk )       

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the solar radiation fluxes 
!  at the snow-ice, ice-water, air-water, 
!  mixed layer-thermocline and water column-bottom sediment interfaces,
!  the mean solar radiation flux over the mixed layer,
!  and the mean solar radiation flux over the thermocline.
!
!------------------------------------------------------------------------------

! Declarations

!  Input (procedure arguments)

REAL (KIND = wp)    , INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  albedo_water                      , & ! Albedo of the water surface 
  albedo_ice                        , & ! Albedo of the ice surface
  albedo_snow                           ! Albedo of the snow surface

TYPE (opticpar_medium), INTENT(IN) :: & 
  opticpar_water                    , & ! Optical characteristics of water
  opticpar_ice                      , & ! Optical characteristics of ice
  opticpar_snow                         ! Optical characteristics of snow 

! Thickness of various layers at the previous time step ("p") and the updated
! values ("n")
REAL (KIND = wp)    , INTENT(IN) ::   &
  &  h_snow_p_flk,                    & !< Snow thickness [m]
  &  h_ice_p_flk,                     & !< Ice thickness [m]
  &  h_ML_p_flk,                      & !< Thickness of the mixed-layer [m]
  &  I_atm_flk                          !< Radiation flux at the lower boundary of 
                                        !< the atmosphere [W m^{-2}], i.e. the
                                        !incident radiation flux 
                                        !< with no regard for the surface albedo 


REAL (KIND = wp),     INTENT(OUT)  :: &
  I_snow_flk                        , & ! Radiation flux through the air-snow
                                        ! interface [W m^{-2}]
  I_bot_flk                         , & ! Radiation flux through the water-bottom
                                        ! sediment interface [W m^{-2}]
  I_ice_flk                         , & ! Radiation flux through the snow-ice or
                                        ! air-ice interface [W m^{-2}]
  I_w_flk                           , & ! Radiation flux through the ice-water or
                                        ! air-water interface [W m^{-2}]
  I_h_flk                           , & ! Radiation flux through the mixed-layer-
                                        ! thermocline interface [W m^{-2}]
  I_intm_0_h_flk                    , & ! Mean radiation flux over the mixed layer
                                        ! [W m^{-2}]
  I_intm_h_D_flk                        ! Mean radiation flux over the thermocline
                                        ! [W m^{-2}]
!_cdm>
! In the new version of "flake_radflux",
! where the number of wave-length bands is restricted to two,
! there are no DO loops. Hence a loop index "i" is not required. 
!_cdm<

!==============================================================================

!     An n-band approximation of the exponential decay law for the solar 
!     radiation flux, where n can be any value from 1 to 10, 
!     is not used due to bad vectorization. 
!     In the new version of "flake_radflux",
!     the number of wave-length bands is restricted to two and 
!     the summation over bands is performed explicitly (no DO loops).

!  Original version (!_nu = unused) but should remain in the file for a while

!_nu !==============================================================================
!_nu !  Start calculations
!_nu !------------------------------------------------------------------------------
!_nu
!_nu   IF(h_ice_p_flk >= h_Ice_min_flk) THEN            
!_nu     ! Ice exists
!_nu     IF(h_snow_p_flk >= h_Snow_min_flk) THEN        
!_nu       ! There is snow above the ice
!_nu       I_snow_flk = I_atm_flk*(1._wp-albedo_snow) 
!_nu       I_bot_flk = 0._wp
!_nu !CDIR EXPAND=opticpar_water%nband_optic
!_nu       DO i=1, opticpar_snow%nband_optic
!_nu         I_bot_flk = I_bot_flk +                    & 
!_nu         opticpar_snow%frac_optic(i)*EXP(-opticpar_snow%extincoef_optic(i)  &
!_nu                                    *h_snow_p_flk) 
!_nu       END DO 
!_nu       I_ice_flk  = I_snow_flk*I_bot_flk
!_nu     ELSE                                           
!_nu       ! No snow above the ice 
!_nu       I_snow_flk = I_atm_flk  
!_nu       I_ice_flk  = I_atm_flk*(1._wp-albedo_ice)
!_nu     END IF 
!_nu     I_bot_flk = 0._wp
!_nu     DO i=1, opticpar_ice%nband_optic
!_nu       I_bot_flk = I_bot_flk +                      & 
!_nu       opticpar_ice%frac_optic(i)*EXP(-opticpar_ice%extincoef_optic(i)      &
!_nu                                 *h_ice_p_flk) 
!_nu     END DO 
!_nu     I_w_flk      = I_ice_flk*I_bot_flk
!_nu   ELSE                                             
!_nu     ! No ice-snow cover
!_nu     I_snow_flk   = I_atm_flk  
!_nu     I_ice_flk    = I_atm_flk
!_nu     I_w_flk      = I_atm_flk*(1._wp-albedo_water)
!_nu   END IF 
!_nu 
!_nu   IF(h_ML_p_flk >= h_ML_min_flk) THEN
!_nu     ! Radiation flux at the bottom of the mixed layer
!_nu     I_bot_flk = 0._wp
!_nu     DO i=1, opticpar_water%nband_optic
!_nu       I_bot_flk = I_bot_flk +            & 
!_nu       opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)  &
!_nu                                   *h_ML_p_flk) 
!_nu     END DO 
!_nu     I_h_flk = I_w_flk*I_bot_flk
!_nu   ELSE
!_nu     ! Mixed-layer depth is less then a minimum value
!_nu     I_h_flk = I_w_flk
!_nu   END IF
!_nu 
!_nu   ! Radiation flux at the lake bottom
!_nu   I_bot_flk = 0._wp
!_nu   DO i=1, opticpar_water%nband_optic
!_nu     I_bot_flk = I_bot_flk +              & 
!_nu     opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)    &
!_nu                                 *depth_w) 
!_nu   END DO 
!_nu   I_bot_flk = I_w_flk*I_bot_flk
!_nu 
!_nu   IF(h_ML_p_flk.GE.h_ML_min_flk) THEN
!_nu     ! Integral-mean radiation flux over the mixed layer
!_nu     I_intm_0_h_flk = 0._wp
!_nu     DO i=1, opticpar_water%nband_optic
!_nu       I_intm_0_h_flk = I_intm_0_h_flk +                                    &
!_nu       opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*      &
!_nu       (1._wp - EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk))
!_nu     END DO 
!_nu     I_intm_0_h_flk = I_w_flk*I_intm_0_h_flk/h_ML_p_flk
!_nu   ELSE
!_nu     I_intm_0_h_flk = I_h_flk
!_nu   END IF
!_nu 
!_nu   IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
!_nu     ! Integral-mean radiation flux over the thermocline
!_nu     I_intm_h_D_flk = 0._wp 
!_nu     DO i=1, opticpar_water%nband_optic
!_nu       I_intm_h_D_flk = I_intm_h_D_flk +                                    &
!_nu       opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*      &
!_nu       ( EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk)                 &
!_nu       - EXP(-opticpar_water%extincoef_optic(i)*depth_w) )
!_nu     END DO 
!_nu     I_intm_h_D_flk = I_w_flk*I_intm_h_D_flk/(depth_w-h_ML_p_flk)
!_nu   ELSE
!_nu     I_intm_h_D_flk = I_h_flk
!_nu   END IF
!_nu 
!_nu !------------------------------------------------------------------------------
!_nu !  End calculations
!_nu !------------------------------------------------------------------------------

!  New version

! We rely on Cray inlining the subroutine, otherwise compilation fails with
! derived types (that contain only scalars)

  !$ACC ROUTINE SEQ


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

  IF(h_ice_p_flk >= h_Ice_min_flk) THEN            
    ! Ice exists
    IF(h_snow_p_flk >= h_Snow_min_flk) THEN        
      ! There is snow above the ice
      I_snow_flk = I_atm_flk*(1._wp-albedo_snow) 
      I_bot_flk =                                                                 &
      opticpar_snow%frac_optic(1)                                                 &
      *EXP(-MIN(opticpar_snow%extincoef_optic(1)*h_snow_p_flk, c_maxearg_flk)) +  &
      opticpar_snow%frac_optic(2)                                                 &
      *EXP(-MIN(opticpar_snow%extincoef_optic(2)*h_snow_p_flk, c_maxearg_flk))
      I_ice_flk  = I_snow_flk*I_bot_flk
    ELSE                                           
      ! No snow above the ice 
      I_snow_flk = I_atm_flk  
      I_ice_flk  = I_atm_flk*(1._wp-albedo_ice)
    END IF 
    I_bot_flk =                                                               &
    opticpar_ice%frac_optic(1)                                                &
    *EXP(-MIN(opticpar_ice%extincoef_optic(1)*h_ice_p_flk, c_maxearg_flk)) +  &
    opticpar_ice%frac_optic(2)                                                &
    *EXP(-MIN(opticpar_ice%extincoef_optic(2)*h_ice_p_flk, c_maxearg_flk))
    I_w_flk      = I_ice_flk*I_bot_flk
  ELSE                                             
    ! No ice-snow cover
    I_snow_flk   = I_atm_flk  
    I_ice_flk    = I_atm_flk
    I_w_flk      = I_atm_flk*(1._wp-albedo_water)
  END IF 

  IF(h_ML_p_flk >= h_ML_min_flk) THEN
    ! Radiation flux at the bottom of the mixed layer
    I_bot_flk =                                                                &
    opticpar_water%frac_optic(1)                                               &
    *EXP(-MIN(opticpar_water%extincoef_optic(1)*h_ML_p_flk, c_maxearg_flk)) +  &
    opticpar_water%frac_optic(2)                                               &
    *EXP(-MIN(opticpar_water%extincoef_optic(2)*h_ML_p_flk, c_maxearg_flk))
    I_h_flk = I_w_flk*I_bot_flk
  ELSE
    ! Mixed-layer depth is less then a minimum value
    I_h_flk = I_w_flk
  END IF

  ! Radiation flux at the lake bottom
  I_bot_flk =                                                             &
  opticpar_water%frac_optic(1)                                            &
  *EXP(-MIN(opticpar_water%extincoef_optic(1)*depth_w, c_maxearg_flk)) +  &
  opticpar_water%frac_optic(2)                                            &
  *EXP(-MIN(opticpar_water%extincoef_optic(2)*depth_w, c_maxearg_flk))
  I_bot_flk = I_w_flk*I_bot_flk

  IF(h_ML_p_flk >= h_ML_min_flk) THEN
    ! Integral-mean radiation flux over the mixed layer
    I_intm_0_h_flk =                                                                        &
    opticpar_water%frac_optic(1)/opticpar_water%extincoef_optic(1)*                         &
    (1._wp - EXP(-MIN(opticpar_water%extincoef_optic(1)*h_ML_p_flk, c_maxearg_flk))) +  &
    opticpar_water%frac_optic(2)/opticpar_water%extincoef_optic(2)*                         &
    (1._wp - EXP(-MIN(opticpar_water%extincoef_optic(2)*h_ML_p_flk, c_maxearg_flk)))
    I_intm_0_h_flk = I_w_flk*I_intm_0_h_flk/h_ML_p_flk
  ELSE
    I_intm_0_h_flk = I_h_flk
  END IF

  IF(h_ML_p_flk <= depth_w-h_ML_min_flk) THEN
    ! Integral-mean radiation flux over the thermocline
    I_intm_h_D_flk =                                                           &
    opticpar_water%frac_optic(1)/opticpar_water%extincoef_optic(1)*            &
    ( EXP(-MIN(opticpar_water%extincoef_optic(1)*h_ML_p_flk, c_maxearg_flk))   &
    - EXP(-MIN(opticpar_water%extincoef_optic(1)*depth_w, c_maxearg_flk)) ) +  &
    opticpar_water%frac_optic(2)/opticpar_water%extincoef_optic(2)*            &
    ( EXP(-MIN(opticpar_water%extincoef_optic(2)*h_ML_p_flk, c_maxearg_flk))   &
    - EXP(-MIN(opticpar_water%extincoef_optic(2)*depth_w, c_maxearg_flk)) )
    I_intm_h_D_flk = I_w_flk*I_intm_h_D_flk/(depth_w-h_ML_p_flk)
  ELSE
    I_intm_h_D_flk = I_h_flk
  END IF

!------------------------------------------------------------------------------
!  End calculations
!------------------------------------------------------------------------------

END SUBROUTINE flake_radflux

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

SUBROUTINE flake_driver ( depth_w, depth_bs, T_bs, par_Coriolis,       &
                          extincoef_water_typ,                         &
                          del_time, T_sfc_p, T_sfc_n, idbg,            &
                          T_snow_n_flk, T_ice_n_flk, T_mnw_n_flk,      &
                          T_bot_n_flk, T_wML_n_flk,                    &
                          T_B1_n_flk, C_T_n_flk, h_snow_n_flk,         &
                          h_ice_n_flk, h_ML_n_flk, H_B1_n_flk,         &
                          Q_bot_flk, Q_star_flk, w_star_sfc_flk,       &
                          Q_w_flk, Q_ice_flk, Q_snow_flk, T_snow_p_flk,&
                          T_ice_p_flk, T_mnw_p_flk, T_bot_p_flk,       &
                          T_wML_p_flk, T_B1_p_flk, C_T_p_flk,          &
                          h_snow_p_flk, h_ice_p_flk, h_ML_p_flk,       &
                          H_B1_p_flk, I_bot_flk, I_snow_flk, I_ice_flk,&
                          I_w_flk, I_h_flk, I_intm_h_D_flk,            &
                          I_intm_0_h_flk, u_star_w_flk, dMsnowdt_flk,  &
                          Phi_T_pr0_flk, Phi_I_pr1_flk, Phi_I_pr0_flk,  &
                          C_I_flk, C_TT_flk, C_Q_flk ) 

!------------------------------------------------------------------------------
!
! Description:
!
!  The main driving routine of the lake parameterization scheme FLake 
!  where (most) computations are performed.
!  Advances the surface temperature
!  and other FLake variables one time step.
!  At the moment, the Euler explicit scheme is used.
!
!------------------------------------------------------------------------------

!  Input (procedure arguments)

REAL (KIND = wp)    , INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  depth_bs                          , & ! Depth of the thermally active layer 
                                        ! of bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of bottom 
                                        ! sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  extincoef_water_typ               , & ! "Typical" extinction coefficient of 
                                        ! the lake water [m^{-1}], used to 
                                        ! compute the equilibrium CBL depth
  del_time                          , & ! The model time step [s]
  T_sfc_p                               ! Surface temperature at the previous 
                                        ! time step [K]  
                                        ! (equal to either T_ice, T_snow 
                                        !  or T_wML)
REAL (KIND = wp),     INTENT(IN) ::   &
  Q_ice_flk                         , & ! Heat flux through the snow-ice or air-ice
                                        ! interface [W m^{-2}]
  Q_snow_flk                        , & ! Heat flux through the air-snow interface
                                        ! [W m^{-2}] 
  T_snow_p_flk                      , & ! Temperature at the air-snow interface [K] 
  T_ice_p_flk                       , & ! Temperature at the snow-ice or air-ice
                                        ! interface [K]
  T_mnw_p_flk                       , & ! Mean temperature of the water column [K]
  T_bot_p_flk                       , & ! Temperature at the water-bottom sediment
                                        ! interface [K]
  T_wML_p_flk                       , & ! Mixed-layer temperature [K]
  C_T_p_flk                         , & ! Shape factor (thermocline)
  h_snow_p_flk                      , & ! Snow thickness [m] 
  h_ice_p_flk                       , & ! Ice thickness [m]  
  h_ML_p_flk                        , & ! Thickness of the mixed-layer [m]
  I_bot_flk                         , & ! Radiation flux through the water-bottom
                                        ! sediment interface [W m^{-2}]
  I_snow_flk                        , & ! Radiation flux through the air-snow
                                        ! interface [W m^{-2}] 
  I_ice_flk                         , & ! Radiation flux through the snow-ice or
                                        ! air-ice interface [W m^{-2}]
  I_w_flk                           , & ! Radiation flux through the ice-water or
                                        ! air-water interface [W m^{-2}] 
  I_h_flk                           , & ! Radiation flux through the mixed-layer- 
                                        ! thermocline interface [W m^{-2}]
  I_intm_0_h_flk                    , & ! Mean radiation flux over the mixed layer
                                        ! [W m^{-2}]
  I_intm_h_D_flk                    , & ! Mean radiation flux over the thermocline
                                        ! [W m^{-2}]
  u_star_w_flk                      , & ! Friction velocity in the surface layer of
                                        ! lake water [m s^{-1}] 
  dMsnowdt_flk                          ! The rate of snow accumulation
                                        ! [kg m^{-2} s^{-1}]

!_cdm>
! Notice that "T_sfc_p" is actually not used within "flake_driver" 
! (it is required to compute the time tendency of the lake surface temperature
! but the tendency is computed in "flake_interface"). 
!_cdm<

!  Output (procedure arguments)

!  In- and Output (procedure arguments)
REAL (KIND = wp),   INTENT(INOUT) ::  &
  Q_w_flk                           , & ! Heat flux through the ice-water or
                                        ! air-water interface [W m^{-2}]
  H_B1_p_flk                        , & ! Thickness of the upper layer of bottom
                                        !  sediments [m]
  T_B1_p_flk                        , & ! Temperature at the bottom of the upper
                                        ! layer of the sediments [K] 
  Phi_T_pr0_flk                     , & !< d\Phi_T(0)/d\zeta   (thermocline) [-]
  Phi_I_pr0_flk                     , & !< d\Phi_I(0)/d\zeta_I (ice) [-]
  Phi_I_pr1_flk                     , & !< d\Phi_I(1)/d\zeta_I (ice) [-]
  C_TT_flk                          , & !< Dimensionless parameter (thermocline) [-]
  C_Q_flk                           , & !< Shape factor with respect to heat flux (thermocline) [-]
  C_I_flk                               !< Shape factor (ice) [-]



!  Output (procedure arguments)

REAL (KIND = wp),     INTENT(OUT) ::  &
  T_snow_n_flk,                       & ! Temperature at the air-snow interface [K]
  T_ice_n_flk,                        & ! Temperature at the snow-ice or air-ice
  T_mnw_n_flk,                        & ! Mean temperature of the water column [K]
  T_wML_n_flk,                        & ! Mixed-layer temperature [K]
  T_B1_n_flk,                         & ! Temperature at the bottom of the upper
                                        ! layer of the sediments [K]
  T_bot_n_flk,                        & ! Temperature at the water-bottom sediment
                                        ! interface [K] 
  C_T_n_flk,                          & ! Shape factor (thermocline)
  h_snow_n_flk,                       & ! Snow thickness [m]
  h_ice_n_flk,                        & ! Ice thickness [m]
  h_ML_n_flk,                         & ! Thickness of the mixed-layer [m]
  H_B1_n_flk,                         & ! Thickness of the upper layer of bottom
                                        ! sediments [m]
  T_sfc_n,                            & ! Updated surface temperature [K] 
                                        ! (equal to the updated value of 
                                        !  either T_ice, T_snow or T_wML)
  Q_bot_flk,                          & ! Heat flux through the water-bottom 
                                        ! sediment interface [W m^{-2}]
  Q_star_flk,                         & ! A generalized heat flux scale [W m^{-2}]
  w_star_sfc_flk                        ! Convective velocity scale, using a
                                        ! generalized heat flux scale [m s^{-1}]

!_cdm>
! A debugging key "idbg" is actually not used. 
! However, it is kept as a "flake_driver" argument 
! for (better) compatibility between the ICON and COSMO codes. 
!_cdm<
INTEGER, INTENT(IN) :: idbg   ! for debug output

!  Local variables of type LOGICAL
LOGICAL ::          &
  l_ice_create    , & ! .TRUE. = ice does not exist but should be created
  l_snow_exists   , & ! .TRUE. = there is snow above the ice
  l_ice_meltabove     ! .TRUE. = snow/ice melting from above takes place

!_cdm>
! DO loop index "i" is currently not used within "flake_driver".
!_cdm<

!  Local variables of type REAL
REAL (KIND = wp)     ::    &
  d_T_mnw_dt             , & ! Time derivative of T_mnw [K s^{-1}] 
  d_T_ice_dt             , & ! Time derivative of T_ice [K s^{-1}] 
  d_T_bot_dt             , & ! Time derivative of T_bot [K s^{-1}] 
  d_T_B1_dt              , & ! Time derivative of T_B1 [K s^{-1}] 
  d_h_snow_dt            , & ! Time derivative of h_snow [m s^{-1}]
  d_h_ice_dt             , & ! Time derivative of h_ice [m s^{-1}]
  d_h_ML_dt              , & ! Time derivative of h_ML [m s^{-1}]
  d_H_B1_dt              , & ! Time derivative of H_B1 [m s^{-1}]
  d_C_T_dt                   ! Time derivative of C_T [s^{-1}]

!  Local variables of type REAL
REAL (KIND = wp)     :: &
  N_T_mean        , & ! The mean buoyancy frequency in the thermocline [s^{-1}] 
  ZM_h_scale      , & ! The ZM96 equilibrium SBL depth scale [m] 
  conv_equil_h_scale  ! The equilibrium CBL depth scale [m]

!  Local variables of type REAL
REAL (KIND = wp)     :: &
  h_ice_threshold, & ! If h_ice<h_ice_threshold, use quasi-equilibrium ice model 
  flk_str_1      , & ! Help storage variable
  flk_str_2      , & ! Help storage variable
  R_H_icesnow    , & ! Dimensionless ratio, used to store intermediate results
  R_rho_c_icesnow, & ! Dimensionless ratio, used to store intermediate results
  R_TI_icesnow   , & ! Dimensionless ratio, used to store intermediate results
  R_Tstar_icesnow    ! Dimensionless ratio, used to store intermediate results

! We rely on Cray inlining the subroutine, otherwise compilation fails with
! derived types (that contain only scalars)

!$ACC ROUTINE SEQ


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!_cdm>
! Security. Set time-rate-of-change of prognostic variables to zero.
! Set prognostic variables to their values at the previous time step.
! (This is to avoid spurious changes of prognostic variables 
! when FLake is used within a 3D model, e.g. to avoid spurious generation of ice 
! at the neighbouring lake points as noticed by Burkhardt Rockel.)
!_cdm<

d_T_mnw_dt   = 0._wp 
d_T_ice_dt   = 0._wp 
d_T_bot_dt   = 0._wp 
d_T_B1_dt    = 0._wp 
d_h_snow_dt  = 0._wp 
d_h_ice_dt   = 0._wp 
d_h_ML_dt    = 0._wp 
d_H_B1_dt    = 0._wp 
d_C_T_dt     = 0._wp 
T_snow_n_flk = T_snow_p_flk   
T_ice_n_flk  = T_ice_p_flk    
T_wML_n_flk  = T_wML_p_flk   
T_mnw_n_flk  = T_mnw_p_flk     
T_bot_n_flk  = T_bot_p_flk  
T_B1_n_flk   = T_B1_p_flk      
h_snow_n_flk = h_snow_p_flk 
h_ice_n_flk  = h_ice_p_flk   
h_ML_n_flk   = h_ML_p_flk    
H_B1_n_flk   = H_B1_p_flk   
C_T_n_flk    = C_T_p_flk    

!------------------------------------------------------------------------------
!  Compute fluxes, using variables from the previous time step.
!------------------------------------------------------------------------------

!_cdm>
! At this point, the heat and radiation fluxes, namely,
! Q_snow_flk, Q_ice_flk, Q_w_flk, 
! I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk,     
! the mean radiation flux over the mixed layer, I_intm_0_h_flk, 
! and the mean radiation flux over the thermocline, I_intm_h_D_flk, 
! should be known.
! They are computed within "flake_interface" (or within the driving model)
! and are available to "flake_driver"
! through the above variables declared in MODULE "flake".
! If a lake is ice-covered, Q_w_flk is re-computed below.
!_cdm<

! Heat flux through the ice-water interface
IF(h_ice_p_flk >= h_Ice_min_flk) THEN    
  ! Ice exists 
  IF(h_ML_p_flk <= h_ML_min_flk) THEN    
    ! Mixed-layer depth is zero, compute flux 

    ! Flux with linear T(z) 
    Q_w_flk = -tpl_kappa_w*(T_bot_p_flk-T_wML_p_flk)/depth_w  

    ! d\Phi(0)/d\zeta (thermocline)
    Phi_T_pr0_flk = Phi_T_pr0_1*C_T_p_flk-Phi_T_pr0_2         

    ! Account for an increased d\Phi(0)/d\zeta 
    Q_w_flk = Q_w_flk*MAX(Phi_T_pr0_flk, 1._wp)           

  ELSE                    

    ! Mixed-layer depth is greater than zero, set flux to zero
    Q_w_flk = 0._wp                  

  END IF   
END IF   

! A generalized heat flux scale 
Q_star_flk = Q_w_flk + I_w_flk + I_h_flk - 2._wp*I_intm_0_h_flk

! Heat flux through the water-bottom sediment interface
IF(lflk_botsed_use) THEN
  Q_bot_flk = -tpl_kappa_w                                                   &
           *(T_B1_p_flk-T_bot_p_flk)/MAX(H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
ELSE  
  Q_bot_flk = 0._wp   ! The bottom-sediment scheme is not used
END IF


!------------------------------------------------------------------------------
!  Check if ice exists or should be created.
!  If so, compute the thickness and the temperature of ice and snow.
!------------------------------------------------------------------------------

!_cdm>
! Notice that a quasi-equilibrium ice-snow model is used 
! to avoid numerical instability when the ice is thin.
! This is always the case when new ice is created.
!_cdm<

!_dev>
! The dependence of snow density and of snow heat conductivity 
! on the snow thickness is accounted for parametrically.
! That is, the time derivatives of \rho_S and \kappa_S are neglected.
! The exception is the equation for the snow thickness 
! in case of snow accumulation and no melting, 
! where d\rho_S/dt is incorporated.
! Furthermore, some (presumably small) correction terms incorporating 
! the snow density and the snow heat conductivity are dropped out.
! Those terms may be included as better formulations 
! for \rho_S and \kappa_S are available.
!_dev<

! Default values
l_ice_create    = .FALSE.  
l_ice_meltabove = .FALSE.  

Ice_exist: IF(h_ice_p_flk < h_Ice_min_flk) THEN   
  ! Ice does not exist 

  l_ice_create = (T_wML_p_flk <= (tpl_T_f+c_small_flk)) .AND.    &
                 (Q_w_flk     <   0._wp)

  IF(l_ice_create) THEN                            
    ! Ice does not exist but should be created
    d_h_ice_dt = -Q_w_flk/tpl_rho_I/tpl_L_f                                  

    ! Advance h_ice 
    h_ice_n_flk = h_ice_p_flk + d_h_ice_dt*del_time                          

    ! Ice temperature
    T_ice_n_flk = tpl_T_f + h_ice_n_flk*Q_w_flk/tpl_kappa_I/Phi_I_pr0_lin    
    d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min 

    ! Advance h_snow
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                       

    ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin                                    & 
                  + Phi_I_ast_MR*MIN(1._wp, h_ice_n_flk/H_Ice_max)       

!CDIR NEXPAND      
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin                        &
                * tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk)    &
                * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)

    ! Snow temperature
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)           
  END IF

ELSE Ice_exist                                     
  ! Ice exists

  ! Check if there is snow above the ice
  l_snow_exists = h_snow_p_flk >= h_Snow_min_flk   

  Melting: IF(T_snow_p_flk >= (tpl_T_f-c_small_flk)) THEN  
    ! T_sfc = T_f, check for melting from above
    ! T_snow = T_ice if snow is absent 

    IF (l_snow_exists) THEN   
      ! There is snow above the ice
      flk_str_1 = Q_snow_flk + I_snow_flk - I_ice_flk    ! Atmospheric forcing
      IF(flk_str_1 >= 0._wp) THEN  ! Melting of snow and ice from above
        l_ice_meltabove = .TRUE.
        d_h_snow_dt = (-flk_str_1/tpl_L_f+dMsnowdt_flk)/             &
                                         flake_snowdensity(h_snow_p_flk)
        d_h_ice_dt  = -(I_ice_flk - I_w_flk - Q_w_flk)/tpl_L_f/tpl_rho_I 
      END IF 
    ELSE                     
      ! No snow above the ice
      ! Atmospheric forcing + heating from the water
      flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk - Q_w_flk  

      IF(flk_str_1.GE.0._wp) THEN  
        ! Melting of ice from above, snow accumulation may occur
        l_ice_meltabove = .TRUE.
        d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
        d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min
      END IF 
    END IF 
    IF(l_ice_meltabove) THEN  ! Melting from above takes place
      h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt *del_time  ! Advance h_ice
      h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time  ! Advance h_snow
      T_ice_n_flk  = tpl_T_f            ! Set T_ice to the freezing point
      T_snow_n_flk = tpl_T_f            ! Set T_snow to the freezing point
    END IF

  END IF Melting

  No_Melting: IF(.NOT.l_ice_meltabove) THEN     ! No melting from above

    d_h_snow_dt = flake_snowdensity(h_snow_p_flk)  
    IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
     flk_str_1 = h_snow_p_flk*tpl_Gamma_rho_S/tpl_rho_w_r
     flk_str_1 = flk_str_1/(1._wp-flk_str_1)
    ELSE                                     ! Snow density is equal to its 
                                             ! maximum value, d\rho_S/dt=0
     flk_str_1 = 0._wp
    END IF

    ! Snow accumulation
    d_h_snow_dt  = dMsnowdt_flk/d_h_snow_dt/(1._wp+flk_str_1)       

    ! Advance h_snow
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                         
    
    ! h_ice relative to its maximum value
    Phi_I_pr0_flk = h_ice_p_flk/H_Ice_max                              

    ! Shape factor (ice)
    C_I_flk = C_I_lin - C_I_MR*(1._wp+Phi_I_ast_MR)*Phi_I_pr0_flk  

    ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0_flk         

    ! d\Phi_I(0)/d\zeta_I (ice)
    Phi_I_pr0_flk = Phi_I_pr0_lin - Phi_I_pr0_flk                      

    h_ice_threshold = MAX(1._wp, 2._wp*C_I_flk*tpl_c_I*              &
                                               (tpl_T_f-T_ice_p_flk)/tpl_L_f)
    h_ice_threshold = Phi_I_pr0_flk/C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*   &
                                                              h_ice_threshold
    ! Threshold value of h_ice
    h_ice_threshold = SQRT(h_ice_threshold*del_time)                   

    ! h_ice(threshold) < 0.9*H_Ice_max
    h_ice_threshold = MIN(0.9_wp*H_Ice_max,                              &
                          MAX(h_ice_threshold, h_Ice_min_flk))

    IF(h_ice_p_flk < h_ice_threshold) THEN  
      ! Use a quasi-equilibrium ice model

      IF(l_snow_exists) THEN   ! Use fluxes at the air-snow interface
        flk_str_1 = Q_snow_flk + I_snow_flk - I_w_flk
      ELSE                     ! Use fluxes at the air-ice interface
        flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk
      END IF
      d_h_ice_dt = -(flk_str_1-Q_w_flk)/tpl_L_f/tpl_rho_I
      ! Advance h_ice
      h_ice_n_flk = h_ice_p_flk + d_h_ice_dt *del_time                         
      ! Ice temperature
      T_ice_n_flk = tpl_T_f + h_ice_n_flk*flk_str_1/tpl_kappa_I/Phi_I_pr0_flk

    ELSE                                     
      ! Use a complete ice model

      d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk
      d_h_ice_dt = (Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
      ! Advance h_ice
      h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt*del_time

      ! Dimensionless parameters
      R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f
      R_Tstar_icesnow = 1._wp - C_I_flk

      IF(l_snow_exists) THEN  
        ! There is snow above the ice
!CDIR NEXPAND      
        R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin                         &
                        * tpl_kappa_I/flake_snowheatconduct(h_snow_p_flk) &
                        * h_snow_p_flk/h_ice_p_flk
        R_rho_c_icesnow = flake_snowdensity(h_snow_p_flk)*                &
                          tpl_c_S/tpl_rho_I/tpl_c_I 
!_dev> 
! These terms should be included as an improved understanding of the snow 
! scheme is gained, of the effect of snow density in particular. 
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                               &
!_nu                        + (1._wp+C_S_lin*h_snow_p_flk/h_ice_p_flk)  &
!_nu                        *  R_H_icesnow*R_rho_c_icesnow
!_dev<

        ! Dimensionless parameter
        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow   

!_dev>
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                               &
!_nu                        + (1._wp-R_rho_c_icesnow)*tpl_c_I *         &
!_nu                          T_ice_p_flk/tpl_L_f
!_dev<
        ! Atmospheric fluxes
        flk_str_2 = Q_snow_flk+I_snow_flk-I_w_flk                  
        flk_str_1  = C_I_flk*h_ice_p_flk + (1._wp+C_S_lin*R_H_icesnow)  &
                            *R_rho_c_icesnow*h_snow_p_flk

        ! Effect of snow accumulation
        d_T_ice_dt = -(1._wp-2._wp*C_S_lin)*R_H_icesnow             &
                        *(tpl_T_f-T_ice_p_flk)                              &
                        * tpl_c_S*dMsnowdt_flk                          
      ELSE
        ! No snow above the ice

        ! Dimensionless parameter
        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow  

        ! Atmospheric fluxes
        flk_str_2 = Q_ice_flk+I_ice_flk-I_w_flk                    
        flk_str_1  = C_I_flk*h_ice_p_flk
        d_T_ice_dt = 0._wp

      END IF 
      ! Add flux due to heat conduction
      d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk&
                    * Phi_I_pr0_flk * (1._wp-R_Tstar_icesnow)                     
      ! Add flux from water to ice
      d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w_flk            

      ! Add atmospheric fluxes
      d_T_ice_dt = d_T_ice_dt + flk_str_2                          

      ! Total forcing
      d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    

      ! dT_ice/dt 
      d_T_ice_dt = d_T_ice_dt/flk_str_1                            

      ! Advance T_ice
      T_ice_n_flk = T_ice_p_flk + d_T_ice_dt*del_time                          
    END IF

    ! h_ice relative to its maximum value
    Phi_I_pr1_flk = MIN(1._wp, h_ice_n_flk/H_Ice_max)          

    ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1_flk     
!CDIR NEXPAND      
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin                             &
                  * tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk)       &
                  * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)

    ! Snow temperature
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)             

  END IF No_Melting

END IF Ice_exist   

! Security, limit h_ice by its maximum value
h_ice_n_flk = MIN(h_ice_n_flk, H_Ice_max)      

! Security, limit the ice and snow temperatures by the freezing point 
T_snow_n_flk = MIN(T_snow_n_flk, tpl_T_f)  
T_ice_n_flk =  MIN(T_ice_n_flk,  tpl_T_f)    

! Security, avoid too low values 
  T_snow_n_flk = MAX(T_snow_n_flk, 73.15_wp)  
  T_ice_n_flk =  MAX(T_ice_n_flk,  73.15_wp)    

! Remove too thin ice and/or snow
IF(h_ice_n_flk < h_Ice_min_flk)  THEN        ! Check ice
  h_ice_n_flk = 0._wp       ! Ice is too thin, remove it, and
  T_ice_n_flk = tpl_T_f         ! set T_ice to the freezing point.
  h_snow_n_flk = 0._wp      ! Remove snow when there is no ice, and
  T_snow_n_flk = tpl_T_f        ! set T_snow to the freezing point.
  l_ice_create = .FALSE.        ! "Exotic" case, ice has been created but 
                                ! proved to be too thin
ELSE IF(h_snow_n_flk < h_Snow_min_flk) THEN  ! Ice exists, check snow
  h_snow_n_flk = 0._wp      ! Snow is too thin, remove it, 
  T_snow_n_flk = T_ice_n_flk    ! and set the snow temperature equal to the 
                                ! ice temperature.
END IF

!------------------------------------------------------------------------------
!  Compute the mean temperature of the water column.
!------------------------------------------------------------------------------

IF (l_ice_create) THEN
  ! Ice has just been created, set Q_w to zero
  Q_w_flk = 0._wp     
ENDIF

d_T_mnw_dt = (Q_w_flk - Q_bot_flk + I_w_flk - I_bot_flk) /         &
                                              tpl_rho_w_r/tpl_c_w/depth_w
! Advance T_mnw
T_mnw_n_flk = T_mnw_p_flk + d_T_mnw_dt*del_time   

! Limit T_mnw by the freezing point 
T_mnw_n_flk = MAX(T_mnw_n_flk, tpl_T_f)           

!------------------------------------------------------------------------------
!  Compute the mixed-layer depth, the mixed-layer temperature, 
!  the bottom temperature and the shape factor
!  with respect to the temperature profile in the thermocline. 
!  Different formulations are used, depending on the regime of mixing. 
!------------------------------------------------------------------------------

HTC_Water: IF (h_ice_n_flk >= h_Ice_min_flk) THEN    ! Ice exists

  ! Limit the mean temperature under the ice by T_r 
  T_mnw_n_flk = MIN(T_mnw_n_flk, tpl_T_r) 

  ! The mixed-layer temperature is equal to the freezing point 
  T_wML_n_flk = tpl_T_f                   

  IF(l_ice_create) THEN                  ! Ice has just been created 
    IF(h_ML_p_flk.GE.depth_w-h_ML_min_flk) THEN ! h_ML=D when ice is created 
      h_ML_n_flk = 0._wp                 ! Set h_ML to zero 
      C_T_n_flk = C_T_min                ! Set C_T to its minimum value 
    ELSE                                 ! h_ML<D when ice is created 
      h_ML_n_flk = h_ML_p_flk            ! h_ML remains unchanged 
      C_T_n_flk = C_T_p_flk              ! C_T (thermocline) remains unchanged 
    END IF 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/        &
                                             (1._wp-h_ML_n_flk/depth_w)
                                         ! Update the bottom temperature 

  ELSE IF(T_bot_p_flk.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, 
                                         ! molecular heat transfer 
    h_ML_n_flk = h_ML_p_flk              ! h_ML remains unchanged 
    C_T_n_flk = C_T_p_flk                ! C_T (thermocline) remains unchanged 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/        &
                                             (1._wp-h_ML_n_flk/depth_w)
                                         ! Update the bottom temperature 

  ELSE                                   ! Ice exists and T_bot = T_r, 
                                         ! convection due to bottom heating 
    T_bot_n_flk = tpl_T_r                ! T_bot is equal to the temperature 
                                         ! of maximum density 
    IF(h_ML_p_flk.GE.c_small_flk) THEN   ! h_ML > 0 
      C_T_n_flk = C_T_p_flk              ! C_T (thermocline) remains unchanged 
      h_ML_n_flk = depth_w*(1._wp-(T_wML_n_flk-T_mnw_n_flk)/            &
                                       (T_wML_n_flk-T_bot_n_flk)/C_T_n_flk)
      h_ML_n_flk = MAX(h_ML_n_flk, 0._wp)   ! Update the mixed-layer depth  
      ! Security, h_ML is equal or very close to the lake depth
      IF(h_ML_n_flk >= (depth_w-h_ML_min_flk)) THEN
        h_ML_n_flk = 0._wp                     ! Set mixed-layer depth to zero
        T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk
                                               ! Adjust the bottom temperature 
      END IF
    ELSE                                 ! h_ML = 0 
      h_ML_n_flk = h_ML_p_flk            ! h_ML remains unchanged 
      C_T_n_flk = (T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk) 
      C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min)) 
                                         ! Update the shape factor (thermocline)
    END IF 
  END IF 

  ! Security, limit the bottom temperature by T_r 
  T_bot_n_flk = MIN(T_bot_n_flk, tpl_T_r)    

ELSE HTC_Water                                      ! Open water

! Generalized buoyancy flux scale and convective velocity scale
  flk_str_1 = flake_buoypar(T_wML_p_flk)*Q_star_flk/tpl_rho_w_r/tpl_c_w                    
  IF(flk_str_1.LT.0._wp) THEN       
    ! Convection     
    w_star_sfc_flk = (-flk_str_1*h_ML_p_flk)**(1._wp/3._wp)
  ELSE 
    ! Neutral or stable stratification
    w_star_sfc_flk = 0._wp
  END IF 

!_cdm>
! The equilibrium depth of the CBL due to surface cooling with the volumetric
! heating is not computed as a solution to the transcendental equation.
! Instead, an algebraic formula is used
! that interpolates between the two asymptotic limits.
!_cdm<

  conv_equil_h_scale = -Q_w_flk/MAX(I_w_flk, c_small_flk)
  IF(conv_equil_h_scale.GT.0._wp .AND. conv_equil_h_scale.LT.1._wp  &
    .AND. T_wML_p_flk.GT.tpl_T_r) THEN   
    ! The equilibrium CBL depth scale is only used above T_r
    conv_equil_h_scale = SQRT(6._wp*conv_equil_h_scale)                 &
              + 2._wp*conv_equil_h_scale/(1._wp-conv_equil_h_scale)
    conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water_typ)
  ELSE
    ! Set the equilibrium CBL depth to zero
    conv_equil_h_scale = 0._wp
  END IF

! Mean buoyancy frequency in the thermocline
  N_T_mean = flake_buoypar(0.5_wp*(T_wML_p_flk+T_bot_p_flk)) *      &
                                      (T_wML_p_flk-T_bot_p_flk)
  IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
    N_T_mean = SQRT(N_T_mean/(depth_w-h_ML_p_flk))  ! Compute N                   
  ELSE 
    N_T_mean = 0._wp                            ! h_ML=D, set N to zero
  END IF 

! The rate of change of C_T
  d_C_T_dt = MAX(w_star_sfc_flk, u_star_w_flk, u_star_min_flk)**2

  ! Relaxation time scale for C_T
  d_C_T_dt = N_T_mean*(depth_w-h_ML_p_flk)**2       &
           / c_relax_C/d_C_T_dt                               

  ! Rate-of-change of C_T 
  d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small_flk)     

! Compute the shape factor and the mixed-layer depth, 
! using different formulations for convection and wind mixing

  ! C_TT, using C_T at the previous time step
  C_TT_flk = C_TT_1*C_T_p_flk-C_TT_2      

  ! C_Q using C_T at the previous time step
  C_Q_flk = 2._wp*C_TT_flk/C_T_p_flk  

  Mixing_regime: IF(flk_str_1.LT.0._wp) THEN  ! Convective mixing 

    ! Update C_T, assuming dh_ML/dt>0
    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time           

    ! Limit C_T 
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))   

    ! Re-compute dC_T/dt
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time           

    IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       ! Compute dh_ML/dt
      IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    
        ! Use a reduced entrainment equation (spin-up)
        d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(w_star_sfc_flk, c_small_flk)

!_dbg>
!       IF (idbg > 10) THEN
!         PRINT *, ' FLake: reduced entrainment eq. D_time*d_h_ML_dt  = ', &
!                                                 d_h_ML_dt*del_time
!         PRINT *, '         w_*       = ', w_star_sfc_flk
!         PRINT *, '         \beta*Q_* = ', flk_str_1
!       ENDIF
!_dbg<

      ELSE
        ! Use a complete entrainment equation 
        R_H_icesnow     = depth_w/h_ML_p_flk
        R_rho_c_icesnow = R_H_icesnow-1._wp
        R_TI_icesnow    = C_T_p_flk/C_TT_flk
        R_Tstar_icesnow = (R_TI_icesnow/2._wp-1._wp)*R_rho_c_icesnow &
                          + 1._wp
        d_h_ML_dt = -Q_star_flk*(R_Tstar_icesnow*(1._wp+c_cbl_1)-1._wp) &
                    - Q_bot_flk
        ! Q_* and Q_b flux terms
        d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w

        flk_str_2 = (depth_w-h_ML_p_flk)*(T_wML_p_flk-T_bot_p_flk)*C_TT_2/   &
                                                              C_TT_flk*d_C_T_dt 

        ! Add dC_T/dt term
        d_h_ML_dt = d_h_ML_dt + flk_str_2

        flk_str_2 = I_bot_flk + (R_TI_icesnow-1._wp)*I_h_flk -           &
                                 R_TI_icesnow*I_intm_h_D_flk
        flk_str_2 = flk_str_2 + (R_TI_icesnow-2._wp)*R_rho_c_icesnow *   &
                                (I_h_flk-I_intm_0_h_flk)
        flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w

        ! Add radiation terms
        d_h_ML_dt = d_h_ML_dt + flk_str_2
        flk_str_2 = -c_cbl_2*R_Tstar_icesnow*Q_star_flk/tpl_rho_w_r/tpl_c_w/ &
                                   MAX(w_star_sfc_flk, c_small_flk)
        flk_str_2 = flk_str_2 + C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)

        ! dh_ML/dt = r.h.s.
        d_h_ML_dt = d_h_ML_dt/flk_str_2
      END IF 

!_cdm>
! Notice that dh_ML/dt may appear to be negative  
! (e.g. due to buoyancy loss to bottom sediments and/or
! the effect of volumetric radiation heating),
! although a negative generalized buoyancy flux scale indicates 
! that the equilibrium CBL depth has not yet been reached
! and convective deepening of the mixed layer should take place. Physically, 
! this situation reflects an approximate character of the lake model.
! Using the self-similar temperature profile in the thermocline, 
! there is always communication between the mixed layer, the thermocline 
! and the lake bottom. As a result, the rate of change of the CBL depth
! is always dependent on the bottom heat flux and the radiation heating of the 
! thermocline. In reality, convective mixed-layer deepening may be completely 
! decoupled from the processes underneath. In order to account for this fact,
! the rate of CBL deepening is set to a small value
! if dh_ML/dt proves to be negative.
! This is "double insurance" however, 
! as a negative dh_ML/dt is encountered very rarely.
!_cdm<

!_dbg>
!     IF (idbg > 10) THEN
!       IF(d_h_ML_dt.LT.0._wp) THEN 
!         PRINT *, 'FLake: negative d_h_ML_dt during convection, = ', d_h_ML_dt
!         PRINT *, '                d_h_ML_dt*del_time = ',        &
!                                          MAX(d_h_ML_dt, c_small_flk)*del_time
!         PRINT *, '         u_*       = ', u_star_w_flk   
!         PRINT *, '         w_*       = ', w_star_sfc_flk
!         PRINT *, '         h_CBL_eqi = ', conv_equil_h_scale
!         PRINT *, '         ZM scale  = ', ZM_h_scale
!         PRINT *, '        h_ML_p_flk = ', h_ML_p_flk
!       END IF
!       PRINT *, 'FLake: Convection, = ', d_h_ML_dt
!       PRINT *, '         Q_*       = ', Q_star_flk
!       PRINT *, '         \beta*Q_* = ', flk_str_1
!     ENDIF
!_dbg<

      d_h_ML_dt = MAX(d_h_ML_dt, c_small_flk)    

      ! Update h_ML 
      h_ML_n_flk = h_ML_p_flk + d_h_ML_dt*del_time

      ! Security, limit h_ML
      h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))
    ELSE
      ! Mixing down to the lake bottom
      h_ML_n_flk = depth_w
    END IF

  ELSE Mixing_regime                              ! Wind mixing

    ! The surface friction velocity
    d_h_ML_dt  = MAX(u_star_w_flk, u_star_min_flk)

    ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*     &
                                                       d_h_ML_dt**2
    ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
    ZM_h_scale = MAX(ZM_h_scale, c_small_flk)
    ZM_h_scale = d_h_ML_dt**3 / ZM_h_scale 

    ! The ZM96 SBL depth scale 
    ZM_h_scale = MAX(h_ML_min_flk, MIN(ZM_h_scale, h_ML_max_flk))

    ! Equilibrium mixed-layer depth 
    ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)

!_cdm> 
! In order to avoid numerical discretization problems,
! an analytical solution to the evolution equation 
! for the wind-mixed layer depth is used.
! That is, an exponential relaxation formula is applied
! over the time interval equal to the model time step.
!_cdm<

    d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
    ! Update h_ML 
    h_ML_n_flk = ZM_h_scale - (ZM_h_scale-h_ML_p_flk)*EXP(-d_h_ML_dt)

    ! Limit h_ML 
    h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))

    ! Re-compute dh_ML/dt
    d_h_ML_dt = (h_ML_n_flk-h_ML_p_flk)/del_time

    IF(h_ML_n_flk.LE.h_ML_p_flk)           &
      ! Mixed-layer retreat or stationary state, dC_T/dt<0
      d_C_T_dt = -d_C_T_dt

    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time           ! Update C_T
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))   ! Limit C_T 
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time           ! Re-compute dC_T/dt


!_dbg>
!   IF (idbg > 10) THEN
!     PRINT *, 'FLake: wind mixing: d_h_ML_dt*del_time = ', d_h_ML_dt*del_time
!     PRINT *, '              h_CBL_eqi = ', conv_equil_h_scale
!     PRINT *, '              ZM scale  = ', ZM_h_scale
!     PRINT *, '              w_*       = ', w_star_sfc_flk
!     PRINT *, '              u_*       = ', u_star_w_flk
!     PRINT *, '             h_ML_p_flk = ', h_ML_p_flk
!   ENDIF
!_dbg<

  END IF Mixing_regime

! Compute the time-rate-of-change of the the bottom temperature, 
! depending on the sign of dh_ML/dt 
! Update the bottom temperature and the mixed-layer temperature

  IF(h_ML_n_flk <= depth_w-h_ML_min_flk) THEN       
    ! Mixing did not reach the bottom 

    IF(h_ML_n_flk > h_ML_p_flk) THEN   
      ! Mixed-layer deepening 
      R_H_icesnow     = h_ML_p_flk/depth_w
      R_rho_c_icesnow = 1._wp-R_H_icesnow 
      R_TI_icesnow    = 0.5_wp*C_T_p_flk*R_rho_c_icesnow+C_TT_flk*     &
                                           (2._wp*R_H_icesnow-1._wp)
      R_Tstar_icesnow = (0.5_wp+C_TT_flk-C_Q_flk)/R_TI_icesnow
      R_TI_icesnow    = (1._wp-C_T_p_flk*R_rho_c_icesnow)/R_TI_icesnow
     
      d_T_bot_dt = (Q_w_flk-Q_bot_flk+I_w_flk-I_bot_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt - C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)*d_h_ML_dt

      ! Q+I fluxes and dh_ML/dt term
      d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w

      flk_str_2 = I_intm_h_D_flk - (1._wp-C_Q_flk)*I_h_flk - C_Q_flk * &
                                                                 I_bot_flk
      flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-h_ML_p_flk)/tpl_rho_w_r/ &
                                                                 tpl_c_w
      ! Add radiation-flux term
      d_T_bot_dt = d_T_bot_dt + flk_str_2

      flk_str_2 = (1._wp-C_TT_2*R_TI_icesnow)/C_T_p_flk
      flk_str_2 = flk_str_2*(T_wML_p_flk-T_bot_p_flk)*d_C_T_dt

      ! Add dC_T/dt term
      d_T_bot_dt = d_T_bot_dt + flk_str_2
      
    ELSE
      ! Mixed-layer retreat or stationary state
      ! dT_bot/dt=0
      d_T_bot_dt = 0._wp                                            
    END IF

    ! Update T_bot  
    T_bot_n_flk = T_bot_p_flk + d_T_bot_dt*del_time   

    ! Security, limit T_bot by the freezing point
    T_bot_n_flk = MAX(T_bot_n_flk, tpl_T_f)

    flk_str_2 = (T_bot_n_flk-tpl_T_r)*(T_mnw_n_flk-tpl_T_r)

    ! Security, avoid T_r crossover 
    IF(flk_str_2 <= 0._wp) T_bot_n_flk = tpl_T_r  

    T_wML_n_flk = C_T_n_flk*(1._wp-h_ML_n_flk/depth_w)
    T_wML_n_flk = (T_mnw_n_flk-T_bot_n_flk*T_wML_n_flk)/(1._wp-T_wML_n_flk)

    ! Security, limit T_wML by the freezing point
    T_wML_n_flk = MAX(T_wML_n_flk, tpl_T_f)

  ELSE
    ! Mixing down to the lake bottom 

    h_ML_n_flk = depth_w
    T_wML_n_flk = T_mnw_n_flk
    T_bot_n_flk = T_mnw_n_flk
    C_T_n_flk = C_T_min

  END IF

  ! In case of unstable stratification, force mixing down to the bottom
  flk_str_2 = (T_wML_n_flk-T_bot_n_flk)*(T_mnw_n_flk-tpl_T_r)
  IF(flk_str_2 < 0._wp) THEN 
    h_ML_n_flk = depth_w
    T_wML_n_flk = T_mnw_n_flk
    T_bot_n_flk = T_mnw_n_flk
    C_T_n_flk = C_T_min
  END IF

END IF HTC_Water

!------------------------------------------------------------------------------
!  Compute the depth of the upper layer of bottom sediments
!  and the temperature at that depth.
!------------------------------------------------------------------------------

Use_sediment: IF(lflk_botsed_use) THEN   ! The bottom-sediment scheme is used
  
  IF (H_B1_p_flk >= depth_bs-H_B1_min_flk) THEN  
    ! No T(z) maximum (no thermal wave) 
    H_B1_p_flk = 0._wp               ! Set H_B1_p to zero
    T_B1_p_flk = T_bot_p_flk             ! Set T_B1_p to the bottom temperature
  END IF 

  flk_str_1 = 2._wp*Phi_B1_pr0/(1._wp-C_B1)*tpl_kappa_w/tpl_rho_w_r/ &
                                                          tpl_c_w*del_time
  ! Threshold value of H_B1
  h_ice_threshold = SQRT(flk_str_1)

  ! Limit H_B1
  h_ice_threshold = MIN(0.9_wp*depth_bs, h_ice_threshold)    

  flk_str_2 = C_B2/(1._wp-C_B2)*(T_bs-T_B1_p_flk)/(depth_bs-H_B1_p_flk)

  IF (H_B1_p_flk < h_ice_threshold) THEN
    ! Use a truncated equation for H_B1(t)
    H_B1_n_flk = SQRT(H_B1_p_flk**2 +flk_str_1)  ! Advance H_B1
    d_H_B1_dt = (H_B1_n_flk-H_B1_p_flk)/del_time          ! Re-compute dH_B1/dt 
  ELSE
    ! Use a full equation for H_B1(t)
    flk_str_1 = (Q_bot_flk+I_bot_flk)/H_B1_p_flk/tpl_rho_w_r/tpl_c_w
    flk_str_1 = flk_str_1 - (1._wp-C_B1)*(T_bot_n_flk-T_bot_p_flk)/del_time
    d_H_B1_dt = (1._wp-C_B1)*(T_bot_p_flk-T_B1_p_flk)/H_B1_p_flk +    &
                                                              C_B1*flk_str_2
    d_H_B1_dt = flk_str_1/d_H_B1_dt
    H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time          ! Advance H_B1
  END IF 
  d_T_B1_dt = flk_str_2*d_H_B1_dt
  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time            ! Advance T_B1


!_dbg>
! IF (idbg > 10) THEN
!   PRINT *, 'BS module: '
!   PRINT *, '  Q_bot   = ', Q_bot_flk
!   PRINT *, '  d_H_B1_dt = ', d_H_B1_dt
!   PRINT *, '  d_T_B1_dt = ', d_T_B1_dt
!   PRINT *, '  H_B1    = ', H_B1_n_flk
!   PRINT *, '    T_bot = ', T_bot_n_flk
!   PRINT *, '  T_B1    = ', T_B1_n_flk
!   PRINT *, '    T_bs  = ',  T_bs
! ENDIF
!_dbg>

! Use a very simplistic procedure, where only the upper layer profile is used, 
! H_B1 is always set to depth_bs, and T_B1 is always set to T_bs. Then, the 
! time derivatives are zero, and the sign of the bottom heat flux depends on 
! whether T_bot is smaller or greater than T_bs.
! This is, of course, an oversimplified scheme.
!_nu  d_H_B1_dt = 0._wp
!_nu  d_T_B1_dt = 0._wp
!_nu  H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
!_nu  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1

  l_snow_exists = (H_B1_n_flk >= depth_bs-H_B1_min_flk)       &
                                             ! H_B1 reached depth_bs, or
             .OR. (H_B1_n_flk <  H_B1_min_flk)                &
                                             ! H_B1 decreased to zero, or
             .OR. ((T_bot_n_flk-T_B1_n_flk)*(T_bs-T_B1_n_flk) <= 0._wp)
                                             ! there is no T(z) maximum
  IF(l_snow_exists) THEN      
    H_B1_n_flk = depth_bs ! Set H_B1 to the depth of the thermally active layer
    T_B1_n_flk = T_bs     ! Set T_B1 to the climatological temperature 
  END IF

ELSE Use_sediment
  ! The bottom-sediment scheme is not used

  ! H_B1 is set to a reference value 
  H_B1_n_flk = rflk_depth_bs_ref   

  ! T_B1 is set to the temperature of maximum density
  T_B1_n_flk = tpl_T_r

END IF Use_sediment

!------------------------------------------------------------------------------
!  Update the surface temperature.
!------------------------------------------------------------------------------

IF     (h_snow_n_flk >= h_Snow_min_flk) THEN   
  ! Snow exists, use the snow temperature
  T_sfc_n = T_snow_n_flk
ELSE IF (h_ice_n_flk >= h_Ice_min_flk) THEN
  ! Ice exists but there is no snow, use the ice temperature
  T_sfc_n = T_ice_n_flk
ELSE 
  ! No ice-snow cover, use the mixed-layer temperature
  T_sfc_n = T_wML_n_flk
END IF

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_driver

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

REAL (KIND = wp)     FUNCTION flake_buoypar (T_water)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the buoyancy parameter,
!  using a quadratic equation of state for the fresh-water.
!  
!------------------------------------------------------------------------------

!  Input (function argument) 
REAL (KIND = wp)    , INTENT(IN) :: &
  T_water                             ! Water temperature [K]

!------------------------------------------------------------------------------
!  Start calculations
!------------------------------------------------------------------------------


  !$ACC ROUTINE SEQ


! Buoyancy parameter [m s^{-2} K^{-1}]

  flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_T_r)

!------------------------------------------------------------------------------
!  End calculations
!------------------------------------------------------------------------------

END FUNCTION flake_buoypar

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

REAL (KIND = wp)     FUNCTION flake_snowdensity (hz_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow density,
!  using an empirical approximation from Heise et al. (2003).
!  
!------------------------------------------------------------------------------

!  Input (function argument) 
REAL (KIND = wp)    , INTENT(IN) :: &
  hz_snow                              ! Snow thickness [m]

!------------------------------------------------------------------------------
!  Start calculations
!------------------------------------------------------------------------------


  !$ACC ROUTINE SEQ


! Snow density [kg m^{-3}]

! Security. Ensure that the expression in () does not become negative at a 
! very large hz_snow.
  flake_snowdensity = MAX( c_small_flk,                                 &
                          (1._wp - hz_snow*tpl_Gamma_rho_S/tpl_rho_w_r) )
  flake_snowdensity = MIN( tpl_rho_S_max, tpl_rho_S_min/flake_snowdensity )

!------------------------------------------------------------------------------
!  End calculations
!------------------------------------------------------------------------------

END FUNCTION flake_snowdensity 

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

REAL (KIND = wp)     FUNCTION flake_snowheatconduct (hz_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow heat conductivity,
!  using an empirical approximation from Heise et al. (2003).
!  
!------------------------------------------------------------------------------
 
!  Input (function argument) 
REAL (KIND = wp)    , INTENT(IN) :: &
  hz_snow                              ! Snow thickness [m]

!------------------------------------------------------------------------------
!  Start calculations
!------------------------------------------------------------------------------


  !$ACC ROUTINE SEQ


! Snow heat conductivity [J m^{-1} s^{-1} K^{-1} = kg m s^{-3} K^{-1}]

  ! Compute snow density
  flake_snowheatconduct = flake_snowdensity( hz_snow )   

  flake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min            &
              + hz_snow*tpl_Gamma_kappa_S*flake_snowheatconduct/tpl_rho_w_r )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowheatconduct

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

!---------------------------------------------------------------------------------------------------
!  End of FLake interface
!===================================================================================================

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

END MODULE sfc_flake

