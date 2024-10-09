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

! Data types and variables used by the AES physics package.
!
! This module contains
! 
!  definition of data types for organising the physical quantities in the
!    AES physics package,
!  the actual variables that are declared of these types, and
!  subroutines for (de-)allocating memory for the variables.
! 
! This module uses derived data types in order to allow for local refinement.

!NEC$ options "-O1"

MODULE mo_aes_phy_memory

  USE mo_kind,                ONLY: dp, wp
  USE mo_impl_constants,      ONLY: SUCCESS, vname_len,        &
    &                               VINTP_METHOD_PRES,         &
    &                               VINTP_METHOD_LIN,          &
    &                               VINTP_METHOD_LIN_NLEVP1
  USE mo_impl_constants,      ONLY: TASK_COMPUTE_PV
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL,    &
    &                               GRID_CELL
  USE mo_exception,           ONLY: message, finish
  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_fortran_tools,       ONLY: t_ptr_2d, t_ptr_3d
  USE mo_parallel_config,     ONLY: nproma
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_name_list_output_config,   ONLY: is_variable_in_output
  USE mtime,                  ONLY: timedelta, OPERATOR(>)
  USE mo_time_config,         ONLY: time_config
  USE mo_aes_phy_config,      ONLY: aes_phy_tc, dt_zero
  USE mo_aes_rad_config,      ONLY: aes_rad_config
  USE mo_aes_vdf_config,      ONLY: aes_vdf_config
  USE mo_aes_sfc_indices,     ONLY: nsfc_type, csfc
  USE mo_model_domain,        ONLY: t_patch
  USE mo_var_list,            ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_metadata,        ONLY: create_vert_interp_metadata, vintp_types, get_timelevel_string
  USE mo_action_types,        ONLY: ACTION_RESET, new_action, actions
  USE mo_nonhydro_state,      ONLY: p_nh_state_lists
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var, t_grib2_int_key, OPERATOR(+)
  USE mo_gribout_config,      ONLY: gribout_config
  USE mo_cdi,                 ONLY: DATATYPE_PACK16, DATATYPE_PACK24,  &
    &                               DATATYPE_FLT32,  DATATYPE_FLT64,   &
    &                               GRID_UNSTRUCTURED, GRID_LONLAT,    &
    &                               TSTEP_INSTANT, TSTEP_CONSTANT,     &
    &                               TSTEP_MIN, TSTEP_MAX,              &
    &                               cdiInqMissval
  USE mo_zaxis_type,          ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF,           &
    &                               ZA_REFERENCE_HALF_HHL,                     &
    &                               ZA_SURFACE, ZA_GENERIC_ICE, ZA_TROPOPAUSE, &
    &                               ZA_HEIGHT_2M, ZA_HEIGHT_10M, ZA_TOA,       &
    &                               ZA_ATMOSPHERE
  USE mo_sea_ice_nml,         ONLY: kice
  USE mo_run_config,          ONLY: iqv ,iqc ,iqi ,     &
    &                               iqr ,iqs ,iqg, iqh, &
    &                               io3
  USE mo_dynamics_config,     ONLY: nnew
  USE mo_advection_config,    ONLY: advection_config

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prm_field, prm_tend                         !< variables
  PUBLIC :: prm_field_list, prm_tend_list               !< variable lists
  PUBLIC :: construct_aes_phy_memory                    !< subroutine
  PUBLIC :: destruct_aes_phy_memory                     !< subroutines
  PUBLIC :: t_aes_phy_field, t_aes_phy_tend             !< derived types

  PUBLIC :: cdimissval

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_aes_phy_memory'

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------
  !>
  !! Derived data type: t_aes_phy_field
  !!
  !! This structure contains two kinds of components:
  !! <ol>
  !! <li> quantities involved in the parameterisation scheme but not in the
  !!      dynamical core, e.g., cloud cover, 10-m wind speed;
  !! <li> atmospheric state variables shared by dynamics and physics, e.g.
  !!      wind, temperature, pressure and tracer concentrations.
  !!      At each time step, the dynamical core provides initial values of
  !!      these quantites on the dynamics grid, which are then interpolated
  !!      to the physics grid and passed on to the physics package.
  !!      In the physics package, these variables may be updated once or
  !!      even more times, depending on the actual numerical schemes chosen
  !!      for the physics-dynamics coupling and the coupling between
  !!      individual parameterisation schemes;
  !! </ol>
  !!
  !! All components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma,           nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy, ntracers)
  !! </ol>
  !! Currently the physics grid has the same spatial resolution as the
  !! dynamics grid, but is unstaggered. This means
  !!
  !!    nlev_phy = nlev
  !!   nblks_phy = patch%nblks_c
  !!
  !! In the long run, the physics grid and dynamics grid may differ in
  !! horizontal and/or vertical resolution, or even in shape.

  TYPE t_aes_phy_field

    ! Metrics
    REAL(wp),POINTER ::     &
      ! horizontal grid
      & clon      (:,:)=>NULL(),    &!< [rad]   longitude at cell center
      & clat      (:,:)=>NULL(),    &!< [rad]   longitude at cell center
      & areacella (:,:)=>NULL(),    &!< [m2]    atmosphere grid-cell area
      ! vertical grid
      & zh        (:,:,:)=>NULL(),  &!< [m]     geometric height at half levels
      & zf        (:,:,:)=>NULL(),  &!< [m]     geometric height at full levels
      & dz        (:,:,:)=>NULL()    !< [m]     geometric height thickness of layer
      
    ! Meteorology and tracers
    REAL(wp),POINTER ::     &
      & ua        (:,:,:)=>NULL(),  &!< [m/s]   zonal wind
      & va        (:,:,:)=>NULL(),  &!< [m/s]   meridional wind
      & wa        (:,:,:)=>NULL(),  &!< [m/s]   vertical velocity
      & ta        (:,:,:)=>NULL(),  &!< [K]     temperature
      & tv        (:,:,:)=>NULL(),  &!< [K]     virtual temperature
      & qtrc_dyn  (:,:,:,:)=>NULL(),&!< [kg/kg] mass fraction of tracer in air
      & qtrc_phy  (:,:,:,:)=>NULL(),&!< [kg/kg] mass fraction of tracer in air
      & mtrcvi    (:,:,:)=>NULL(),  &!< [kg/m2] atmosphere mass content of tracer
      & cptgzvi   (:,:)=>NULL(),    &!< [kg/m2] dry static energy  , vertically integrated through the atmospheric column
      & udynvi    (:,:)=>NULL(),    &!< [kg/m2] vertically integrated moist internal energy -- after dynamics
      & duphyvi   (:,:)=>NULL(),    &!< [kg/m2] change of vertically integrated moist internal energy by physics
      & utmxvi    (:,:)=>NULL(),    &!< [kg/m2] vertically integrated moist internal energy -- after tmix
      & rho       (:,:,:)=>NULL(),  &!< [kg/m3] air density
      & mair      (:,:,:)=>NULL(),  &!< [kg/m2] air content
      & pv        (:,:,:)=>NULL(),  &!< [K/m2/kg/s] Ertel potential vorticity (NWP pp)
      & geoi      (:,:,:)=>NULL(),  &!< [m2/s2] geopotential above ground at half levels (vertical interfaces)
      & geom      (:,:,:)=>NULL(),  &!< [m2/s2] geopotential above ground at full levels (layer ave. or mid-point value)
      & geop      (:,:,:)=>NULL(),  &!< [m2/s2] geopotential at full levels (layer ave. or mid-point value)
      & pfull     (:,:,:)=>NULL(),  &!< [Pa]    air pressure at model levels
      & phalf     (:,:,:)=>NULL()    !< [Pa]    air pressure at model half levels

    ! surface fluxes of internal energy (positive downward)
    REAL(wp),POINTER ::             &
      & ufts      (:,  :)=>NULL(),  &!< energy flux at surface from thermal exchange [K kg/m2/s ] 
      & ufvs      (:,  :)=>NULL(),  &!< energy flux at surface from vapor exchange   [K kg/m2/s ] 
      & ufcs      (:,  :)=>NULL()    !< energy flux at surface from condensate       [K kg/m2/s ] 
    TYPE(t_ptr_2d),ALLOCATABLE :: mtrcvi_ptr(:)

    ! Radiation
    REAL(wp),POINTER ::       &
      !
      & cosmu0      (:,  :)=>NULL(),  &!< [ ]    cos of zenith angle mu0 for radiative heating  calculation
      & cosmu0_rt   (:,  :)=>NULL(),  &!< [ ]    cos of zenith angle mu0 at radiation time step
      & daylght_frc (:,  :)=>NULL(),  &!< [ ]    daylight fraction at each grid point
      & daylght_frc_rt(:,:)=>NULL(),  &!< [ ]    daylight fraction at each grid point for radiation time step
      !
      ! shortwave fluxes
      ! - through the atmosphere at the radiation time step
      & rsd_rt      (:,:,:)=>NULL(),  &!< [W/m2] downwelling shortwave radiation
      & rsu_rt      (:,:,:)=>NULL(),  &!< [W/m2] upwelling   shortwave radiation
      & rsdcs_rt    (:,:,:)=>NULL(),  &!< [W/m2] downwelling clear-sky shortwave radiation
      & rsucs_rt    (:,:,:)=>NULL(),  &!< [W/m2] upwelling   clear-sky shortwave radiation
      ! - at the top of the atmosphere at all times
      & rsdt        (:,  :)=>NULL(),  &!< [W/m2] toa incident shortwave radiation
      & rsut        (:,  :)=>NULL(),  &!< [W/m2] toa outgoing shortwave radiation
      & rsutcs      (:,  :)=>NULL(),  &!< [W/m2] toa outgoing clear-sky shortwave radiation
      & rsnt        (:,  :)=>NULL(),  &!< [W/m2] toa net shortwave radiation
      ! - at the surface at all times
      & rsds        (:,  :)=>NULL(),  &!< [W/m2] surface downwelling shortwave radiation
      & rsus        (:,  :)=>NULL(),  &!< [W/m2] surface upwelling   shortwave radiation
      & rsdscs      (:,  :)=>NULL(),  &!< [W/m2] surface downwelling clear-sky shortwave radiation
      & rsuscs      (:,  :)=>NULL(),  &!< [W/m2] surface upwelling   clear-sky shortwave radiation
      & rsns        (:,  :)=>NULL(),  &!< [W/m2] surface net shortwave radiation
      !
      ! shortwave flux components at the surface
      ! - at radiation times
      & rvds_dir_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  visible            radiation
      & rpds_dir_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  photosynth. active radiation
      & rnds_dir_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  near-infrared      radiation
      & rvds_dif_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse visible            radiation
      & rpds_dif_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse photosynth. active radiation
      & rnds_dif_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse near-infrared      radiation
      & rvus_rt     (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         visible            radiation
      & rpus_rt     (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         photosynth. active radiation
      & rnus_rt     (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         near-infrared      radiation
      ! - at all times
      & rvds_dir    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  visible            radiation
      & rpds_dir    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  photosynth. active radiation
      & rnds_dir    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  near-infrared      radiation
      & rvds_dif    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse visible            radiation
      & rpds_dif    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse photosynth. active radiation
      & rnds_dif    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse near-infrared      radiation
      & rvus        (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         visible            radiation
      & rpus        (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         photosynth. active radiation
      & rnus        (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         near-infrared      radiation
      !
      ! longwave fluxes
      ! - through the atmosphere at the radiation time step
      & rld_rt      (:,:,:)=>NULL(),  &!< [W/m2] downwelling longwave radiation
      & rlu_rt      (:,:,:)=>NULL(),  &!< [W/m2] upwelling   longwave radiation
      & rldcs_rt    (:,:,:)=>NULL(),  &!< [W/m2] downwelling clear-sky longwave radiation
      & rlucs_rt    (:,:,:)=>NULL(),  &!< [W/m2] upwelling   clear-sky longwave radiation
      ! - at the top of the atmosphere at all times
      & rlut        (:,  :)=>NULL(),  &!< [W/m2] toa outgoing longwave radiation
      & rlutcs      (:,  :)=>NULL(),  &!< [W/m2] toa outgoing clear-sky longwave radiation
      & rlnt        (:,  :)=>NULL(),  &!< [W/m2] TOA net longwave radiation
      ! - at the surface at all times
      & rlds        (:,  :)=>NULL(),  &!< [W/m2] surface downwelling longwave radiation
      & rlus        (:,  :)=>NULL(),  &!< [W/m2] surface upwelling   longwave radiation
      & rldscs      (:,  :)=>NULL(),  &!< [W/m2] surface downwelling clear-sky longwave radiation
      & rlns        (:,  :)=>NULL(),  &!< [W/m2] surface net longwave radiation
      & o3          (:,:,:)=>NULL()    !< [mol/mol] ozone volume mixing ratio
    ! effective radius of ice
    REAL(wp), POINTER ::      &
      & x_snow      (:,:,:)=>NULL(),  &!< average mass of snow flake in kg
      & x_ice       (:,:,:)=>NULL(),  &!< average mass of ice crystal in kg
      & acinc       (:,:,:)=>NULL(),  &!< cloud ice number concentration [1/m^3]
      & acsnc       (:,:,:)=>NULL(),  &!< cloud snow number concentration [1/m^3]
      & reff_ice    (:,:,:)=>NULL(),  &!< effective radius of snow in [um]
      & tau_ice     (:,:,:)=>NULL()    !< 3d optical depth of ice in clouds   
    ! effective radius of snow
    REAL(wp), POINTER ::      &
      & reff_snow   (:,:,:)=>NULL(),  &!< effective radius of snow in [um]
      & tau_snow    (:,:,:)=>NULL()    !< 3d optical depth of snow
    ! arbitrary 2d-field in radiation for output
    REAL(wp), POINTER ::      &
      & rad_2d      (:,:)=>NULL()      !< arbitrary 2d field in radiation for output     
    ! aerosol optical properties
    REAL(wp),POINTER ::      &
      & aer_aod_533 (:,:,:)=>NULL(),  &!< aerosol optical depth at 533 nm
      & aer_ssa_533 (:,:,:)=>NULL(),  &!< aerosol single scattering albedo at 533 nm
      & aer_asy_533 (:,:,:)=>NULL(),  &!< aerosol asymmetry factor at 533 nm
      & aer_aod_2325(:,:,:)=>NULL(),  &!< aerosol optical depth at 2325 nm
      & aer_ssa_2325(:,:,:)=>NULL(),  &!< aerosol single scattering albedo at 2325 nm
      & aer_asy_2325(:,:,:)=>NULL(),  &!< aerosol asymmetry factor at 2325 nm
      & aer_aod_9731(:,:,:)=>NULL()    !< effective aerosol optical depth at 9731 nm
            !< the last quantity is in the thermal wavelength ranch, 
            !< the first lie in the solar spectrum

    ! Clouds
    REAL(wp),POINTER ::     &
      & aclc      (:,:,:)=>NULL(),  &!< [m2/m2] cloud area fractional
      & aclcov    (:,  :)=>NULL(),  &!< [m2/m2] total cloud cover
      & acdnc     (:,:,:)=>NULL(),  &!< cloud droplet number concentration [1/m**3]
      & hur       (:,:,:)=>NULL()    !< relative humidity

    ! Precipitation
    REAL(wp),POINTER ::     &
      & rsfl      (:,  :)=>NULL(),  &!< sfc rain    flux, large scale [kg m-2 s-1]
      & ssfl      (:,  :)=>NULL(),  &!< sfc snow    flux, large scale [kg m-2 s-1]
!
! move to two
      & rain_gsp_rate(:,  :)=>NULL(),  &!< gridscale rain rate     [kg m-2 s-1]
      & ice_gsp_rate (:,  :)=>NULL(),  &!< gridscale ice rate      [kg m-2 s-1]
      & snow_gsp_rate(:,  :)=>NULL(),  &!< gridscale snow rate     [kg m-2 s-1]
      & graupel_gsp_rate(:,  :)=>NULL(),  &!< gridscale graupel rate     [kg m-2 s-1]
      & hail_gsp_rate(:,  :)=>NULL(),  &!< gridscale hail rate     [kg m-2 s-1]
!
!
      & pr        (:,  :)=>NULL()    !< precipitation flux         [kg m-2 s-1]

    ! Tropopause
    REAL(wp),POINTER ::     &
      & ptp       (:,  :)=>NULL()    !< tropopause air pressure [Pa]

    REAL(wp),POINTER ::     &
      & siced  (:,  :)=>NULL(),     &!< ice depth
      & alb    (:,  :)=>NULL(),     &!< surface background albedo
      & seaice (:,  :)=>NULL()       !< sea ice as read in from amip input

    ! Energy and moisture budget related diagnostic variables
    REAL(wp),POINTER ::     &
      & cvair    (:,:,:)=>NULL(),   &!< specific heat of air at constant volume   [J/kg/K]
      !
      & q_phy    (:,:,:)=>NULL(),   &!< layer heating by physics [W/m^2]
      & q_phy_vi (:,  :)=>NULL(),   &!< vertically integrated heating by physics [W/m^2]
      !
      & q_rad    (:,:,:)=>NULL(),   &!< Layer heating by LW+SW radiation
      & q_rad_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by LW+SW radiation
      & q_rlw    (:,:,:)=>NULL(),   &!< Layer heating by LW radiation
      & q_rlw_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by LW radiation
      & q_rsw    (:,:,:)=>NULL(),   &!< Layer heating by SW radiation
      & q_rsw_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by SW radiation
      & q_vdf    (:,:,:)=>NULL(),   &!< Layer heating by vertical diffusion
      & q_vdf_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by vertical diffusion
      & q_cld    (:,:,:)=>NULL(),   &!< Layer heating by cloud processes
      & q_cld_vi (:,  :)=>NULL()     !< Vertically integrated heating by cloud processes
      !
!!$      & sh_vdiff (:,  :)=>NULL(),   &!< sensible heat flux of vdiff
!!$      & qv_vdiff (:,  :)=>NULL(),   &!< qv flux of vdiff

    ! JSBACH
    REAL(wp),POINTER :: &
      & ts_rad     (:,  :)=>NULL(),  &!< [K] radiative sfc. temperature for use in radiation
      & ts_rad_rt  (:,  :)=>NULL(),  &!< [K] radiative sfc. temperature at radiation time
      & csat       (:,  :)=>NULL(),  &!<
      & cair       (:,  :)=>NULL(),  &!<
      & q_snocpymlt(:,  :)=>NULL(),  &!< [W/m2] heating used to melt snow on the canopy
      & q_rlw_impl (:,  :)=>NULL(),  &!< [W/m2] heating correction due to implicit land surface coupling
      & q_rlw_nlev (:,  :)=>NULL()    !< [W/m2] heating in the lowest layer

    ! CO2
    REAL(wp),POINTER :: &
      & co2_flux_tile   (:,:,:)=>NULL(),  &!< CO2 flux on tiles (land, ocean)
      & fco2nat         (:,  :)=>NULL()    !< Surface Carbon Mass Flux into the Atmosphere Due to Natural Sources

    TYPE(t_ptr_2d),ALLOCATABLE :: co2_flux_tile_ptr(:)

    ! Sea ice.
    ! See also sea_ice/thermodyn/mo_sea_ice_types.f90
    INTEGER, POINTER :: kice  ! Number of ice-thickness classes
    REAL(wp),POINTER     ::     &
      & Tsurf   (:,:,:)=>NULL(),        &! Ice surface temperature [degC]
      & T1      (:,:,:)=>NULL(),        &! Temperature of upper ice layer [degC]
      & T2      (:,:,:)=>NULL(),        &! Temperature of lower ice layer [degC]
      & hi      (:,:,:)=>NULL(),        &! Ice thickness [m]
      & hs      (:,:,:)=>NULL(),        &! Snow thickness on ice [m]
      & Qtop    (:,:,:)=>NULL(),        &! Energy flux available for surface melting [W/m^2]
      & Qbot    (:,:,:)=>NULL(),        &! Energy flux at ice-ocean interface [W/m^2]
      & conc    (:,:,:)=>NULL(),        &! Ice concentration [0,1]
      & albvisdir_ice(:,:,:)=>NULL(),   &! Ice surface albedo for visible range, direct
      & albvisdif_ice(:,:,:)=>NULL(),   &! Ice surface albedo for visible range, diffuse
      & albnirdir_ice(:,:,:)=>NULL(),   &! Ice surface albedo for near IR range, direct
      & albnirdif_ice(:,:,:)=>NULL()     ! Ice surface albedo for near IR range, diffuse

    ! Turbulence
    REAL(wp),POINTER ::     &
      & totte       (:,:,:)=>NULL(),  &!< total turbulent energy at step n+1
      & tottem0     (:,:,:)=>NULL(),  &!< total turbulent energy at step n
      & tottem1     (:,:,:)=>NULL()    !< total turbulent energy at step n-1

    ! need only for vdiff ++++
    REAL(wp),POINTER ::     &
      & ri_atm    (:,:,:)=>NULL(),  &!< moist Richardson number at layer interfaces
      & mixlen    (:,:,:)=>NULL()    !< mixing length at layer interfaces

    REAL(wp),POINTER ::     &
      & cptgz       (:,:,:)=>NULL()   !< dry static energy

    REAL(wp),POINTER ::      &
      & cfm     (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfm_tile(:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfh     (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfh_tile(:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfv     (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cftotte (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfthv   (:,:,:)=>NULL()       !< turbulent exchange coefficient

    TYPE(t_ptr_2d),ALLOCATABLE :: cfm_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: cfh_tile_ptr(:)

    REAL(wp),POINTER ::     &
      & coriol(:,:)=>NULL(),        &!< Coriolis parameter
      & hdtcbl (:,:)=>NULL(),       &!< height of the top of the atmospheric dry convective boundary layer
      & z0m_tile(:,:,:)=>NULL(),    &!< aerodynamic roughness length (over each surface type)
      & z0h_tile(:,:,:)=>NULL(),    &!< aerodynamic roughness length (over each surface type)
      & z0m   (:,:)=>NULL(),        &!< aerodynamic roughness length (grid box mean)
      & z0h   (:,:)=>NULL(),        &!< aerodynamic roughness length (grid box mean)
      & z0h_lnd(:,:)=>NULL(),       &!< roughness length for heat (over land)
      & ustar (:,:)=>NULL(),        &!<
      & wstar (:,:)=>NULL(),        &!< convective velocity scale
      & wstar_tile(:,:,:)=>NULL(),  &!< convective velocity scale (over each surface type)
      & kedisp(:,:)=>NULL(),        &!< vertically integrated dissipation of kinetic energy
      & ocu   (:,:)=>NULL(),        &!< eastward  velocity of ocean surface current
      & ocv   (:,:)=>NULL()          !< northward velocity of ocean surface current

      !
    REAL(wp),POINTER ::     &
      ! net fluxes at TOA and surface
      & swflxsfc_tile(:,:,:)=>NULL(),  &!< [ W/m2] shortwave net flux at surface
      & lwflxsfc_tile(:,:,:)=>NULL(),  &!< [ W/m2] longwave net flux at surface
      & dlwflxsfc_dT (:,:)  =>NULL()    !< [ W/m2/K] longwave net flux temp tend at surface

    TYPE(t_ptr_2d),ALLOCATABLE :: swflxsfc_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: lwflxsfc_tile_ptr(:)

    TYPE(t_ptr_2d),ALLOCATABLE :: z0m_tile_ptr(:), z0h_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: wstar_tile_ptr(:)

    ! need only for vdiff ----

    ! Surface variables

    REAL(wp),POINTER :: &
      & orog(:,:)=>NULL(),          &!< surface altitude [m]
      & sftlf (:,:)=>NULL(),        &!< cell area fraction occupied by land including lakes (1. = land and/or lakes only, 0. = ocean only)
      & sftgif(:,:)=>NULL(),        &!< cell area fraction occupied by land ice             (1. = land ice only, 0. = no land ice)
      & sftof (:,:)=>NULL(),        &!< cell area fraction occupied by ocean                (1. = ocean only, 0. = land and/or lakes only)
      & lsmask(:,:)=>NULL(),        &!< cell area fraction occupied by land excluding lakes (1. = land, 0. = ocean or lake only) 
      & alake (:,:)=>NULL(),        &!< cell area fraction occupied by lakes
      & glac  (:,:)=>NULL(),        &!< land area fraction that is glaciated
      & lake_ice_frc(:,:)=>NULL(),  &!< lake area fraction that is ice covered
      & icefrc(:,:)=>NULL(),        &!< ice cover given as the fraction of grid box (friac  in memory_g3b)
      & ts_tile(:,:,:)=>NULL(),     &!< surface temperature over land/water/ice
      & ts     (:,  :)=>NULL(),     &!< surface temperature, grid box mean
      & qs_sfc_tile(:,:,:)=>NULL()   !< saturation specific humidity at surface 

    TYPE(t_ptr_2d),ALLOCATABLE :: ts_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: qs_sfc_tile_ptr(:)

    ! Surface albedo
    REAL(wp),POINTER :: &
      & albvisdir_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for visible range, direct
      & albvisdif_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for visible range, diffuse
      & albnirdir_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for near-IR range, direct
      & albnirdif_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for near-IR range, diffuse
      & albedo_tile    (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles
      & albvisdir      (:,:  )=>NULL(),  &!< [ ] surface albedo for visible range, direct, grid-box mean
      & albvisdif      (:,:  )=>NULL(),  &!< [ ] surface albedo for visible range, diffuse, grid-box mean
      & albnirdir      (:,:  )=>NULL(),  &!< [ ] surface albedo for near-IR range, direct, grid-box mean
      & albnirdif      (:,:  )=>NULL(),  &!< [ ] surface albedo for near-IR range, diffuse, grid-box mean
      & albedo         (:,:  )=>NULL()    !< [ ] surface albedo, grid-box mean

    ! Surface emissivity
    REAL(wp),POINTER :: &
      & emissivity     (:,:  )=>NULL()    !< [ ] surface emissivity, grid-box mean

    TYPE(t_ptr_2d),ALLOCATABLE :: albvisdir_tile_ptr(:), albvisdif_tile_ptr(:), &
      & albnirdir_tile_ptr(:), albnirdif_tile_ptr(:), albedo_tile_ptr(:)

    REAL(wp),POINTER :: &
      & lhflx     (:,  :)=>NULL(),    &!< grid box mean latent   heat flux at surface
      & shflx     (:,  :)=>NULL(),    &!< grid box mean sensible heat flux at surface
      & evap      (:,  :)=>NULL(),    &!< grid box mean evaporation at surface
      & lhflx_tile(:,:,:)=>NULL(),    &!< latent   heat flux at surface on tiles
      & shflx_tile(:,:,:)=>NULL(),    &!< sensible heat flux at surface on tiles
      & evap_tile (:,:,:)=>NULL(),    &!< evaporation at surface on tiles
      & frac_tile (:,:,:)=>NULL()      !< surface fraction of tiles:
                                       !  - fraction of land without lakes
                                       !  - fraction of ice covered water in the grid box, for sea and lakes
                                       !  - fraction of open water in the grid box, for sea and lakes

    TYPE(t_ptr_2d),ALLOCATABLE :: lhflx_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: shflx_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: evap_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: frac_tile_ptr(:)

    REAL(wp),POINTER :: &
      & u_stress     (:,  :)=>NULL(), &!< grid box mean wind stress
      & v_stress     (:,  :)=>NULL(), &!< grid box mean wind stress
      & u_stress_tile(:,:,:)=>NULL(), &!< wind stress on tiles
      & v_stress_tile(:,:,:)=>NULL()   !< wind stress on tiles

    TYPE(t_ptr_2d),ALLOCATABLE :: u_stress_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: v_stress_tile_ptr(:)

    ! Near surface diagnostics (2m temp; 2m dew point temp; 10m wind)
    !
    REAL(wp),POINTER ::        &
      & sfcwind     (:,  :)=>NULL(),   &!< grid box mean 10 m wind
      & uas         (:,  :)=>NULL(),   &!< grid box mean 10m u-velocity
      & vas         (:,  :)=>NULL(),   &!< grid box mean 10m v-velocity
      & tas         (:,  :)=>NULL(),   &!< grid box mean 2m temperature
      & dew2        (:,  :)=>NULL(),   &!< grid box mean 2m dew point temperature
      & qv2m        (:,  :)=>NULL(),   &!< grid box mean 2m specific humidity
      & tasmax      (:,  :)=>NULL(),   &!< grid box mean maximum 2m temperature
      & tasmin      (:,  :)=>NULL(),   &!< grid box mean minimum 2m temperature
      & sfcwind_tile(:,:,:)=>NULL(),   &!< 10 m wind on tiles
      & uas_tile    (:,:,:)=>NULL(),   &!< 10m u-velocity on tiles
      & vas_tile    (:,:,:)=>NULL(),   &!< 10m v-velocity on tiles
      & tas_tile    (:,:,:)=>NULL(),   &!< 2m temperature on tiles
      & dew2_tile   (:,:,:)=>NULL(),   &!< 2m dew point temperature on tiles
      & qv2m_tile   (:,:,:)=>NULL()     !< 2m specific humidity on tiles

    TYPE(t_ptr_2d),ALLOCATABLE :: sfcwind_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: uas_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: vas_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: tas_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: dew2_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: qv2m_tile_ptr(:)

    ! global diagnostics
    REAL(wp),POINTER ::       &
      !
      & tas_gmean    (:)=>NULL(),      &!< [K] global mean 2m-temperature
      & rsdt_gmean   (:)=>NULL(),      &!< [W/m2] global mean toa incident shortwave radiation
      & rsut_gmean   (:)=>NULL(),      &!< [W/m2] global mean toa outgoing shortwave radiation
      & rlut_gmean   (:)=>NULL(),      &!< [W/m2] global mean toa outgoing longwave radiation
      & prec_gmean   (:)=>NULL(),      &!< [kg/m2/s] global mean precipitation flux
      & evap_gmean   (:)=>NULL(),      &!< [kg/m2/s] global mean evaporation flux
      & radtop_gmean (:)=>NULL(),      &!< [W/m2] global mean toa net total radiation, derived variable
      & radbot_gmean (:)=>NULL(),      &!< [W/m2] global mean surface net total radiation, derived variable
      & radbal_gmean (:)=>NULL(),      &!< [W/m2] global mean net radiative flux into atmosphere, derived variable
      & fwfoce_gmean (:)=>NULL(),      &!< [kg/m2/s] global mean freshwater flux over ocean area, derived variable
      & icefrc_gmean (:)=>NULL(),      &!< global mean ice cover given as the fraction of grid box
      & udynvi_gmean (:)=>NULL(),      &!< [kg/m2] global mean vertically integrated moist internal energy - after dynamics
      & duphyvi_gmean(:)=>NULL(),      &!< [kg/m2] global mean vertically integrated moist internal energy change by physics
      & utmxvi_gmean (:)=>NULL(),      &!< [kg/m2] global mean vertically integrated moist internal energy - after tmx
      & ufts_gmean   (:)=>NULL(),      &!< [K kg/m2/s] global mean energy flux at surface from thermal exchange
      & ufvs_gmean   (:)=>NULL(),      &!< [K kg/m2/s] global mean energy flux at surface from vapor exchange
      & ufcs_gmean   (:)=>NULL(),      &!< [K kg/m2/s] global mean energy flux at surface from condensate
      & kedisp_gmean (:)=>NULL(),      &!< [W/m2] global mean of vertically integrated dissipation of kinetic energy
      & uphybal_gmean   (:)=>NULL()     !< [W/m2] global energy balance in aes physics (positive: gain of internal energy)

  END TYPE t_aes_phy_field

  !>
  !! Data type containing the tendencies returned by the individual parameterizations
  !!
  !! The components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma, nlev_phy, nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy, ntracers)
  !! </ol>
  !!
  TYPE t_aes_phy_tend

    REAL(wp), POINTER ::   &
      !
      ! tendency due to model physics in:
      !
      &   ua_phy (:,:,:)=>NULL()  , & !< [m/s2]    u-wind
      &   va_phy (:,:,:)=>NULL()  , & !< [m/s2]    v-wind
      &   wa_phy (:,:,:)=>NULL()  , & !< [m/s2]    w-wind
      &   ta_phy (:,:,:)=>NULL()  , & !< [K/s]     temperature (for const. volume)
      & qtrc_phy (:,:,:,:)=>NULL(), & !< [kg/kg/s] mass fraction of tracer in air
      & mtrcvi_phy(:,:,  :)=>NULL(),& !< [kg/m2/s] atmosphere mass content of tracer
      & utmxvi   (:,:)=>NULL(),     & !< [J/m2/s] vertically integrated moist internal energy -- after tmix

      ! tendency due to vertical diffusion ("vdiff") in:
      !
      &   ta_vdf (:,:,:)=>NULL()  , & !< temperature (for const. pressure)
      &   ua_vdf (:,:,:)=>NULL()  , & !< u-wind 
      &   va_vdf (:,:,:)=>NULL()  , & !< v-wind
      &   wa_vdf (:,:,:)=>NULL()  , & !< w-wind
      & qtrc_vdf (:,:,:,:)=>NULL(), & !< mass fraction of tracer in air
      !
      ! tendency due to surface physics in
      !
      &   ta_sfc (:,:)=>NULL()  , & !< temperature (for const. pressure)
      !
      ! tendency due to radiation in:
      !
      &   ta_rsw (:,:,:)=>NULL()  , & !< temperature, due to shortwave radiation (for const. pressure)
      &   ta_rlw (:,:,:)=>NULL()  , & !< temperature, due to longwave  radiation (for const. pressure)
      &   ta_rad (:,:,:)=>NULL()  , & !< temperature, due to SW + LW   radiation (for const. pressure)
      &   ta_rlw_impl(:,:)=>NULL(), & !< temperature, due to implicit land surface temperature change (for const. pressure)
      !
      ! tendency due to Cariolle linearized ozone in:
      ! 
      &   o3_car (:,:,:)=>NULL()      !< mass fraction of ozone in air

    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_phy_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_vdf_ptr(:)
              
    TYPE(t_ptr_2d),ALLOCATABLE :: mtrcvi_phy_ptr(:)

  END TYPE t_aes_phy_tend

  !!--------------------------------------------------------------------------
  !!                          STATE VARIABLES 
  !!--------------------------------------------------------------------------
  !! The variable names have the prefix "prm_" in order to emphasize that they
  !! are defined for and used in parameterisations.

  TYPE(t_aes_phy_field),ALLOCATABLE,TARGET :: prm_field(:)  !< shape: (n_dom)
  TYPE(t_aes_phy_tend ),ALLOCATABLE,TARGET :: prm_tend (:)  !< shape: (n_dom)

  !!--------------------------------------------------------------------------
  !!                          VARIABLE LISTS
  !!--------------------------------------------------------------------------
  TYPE(t_var_list_ptr),ALLOCATABLE :: prm_field_list(:)  !< shape: (n_dom)
  TYPE(t_var_list_ptr),ALLOCATABLE :: prm_tend_list (:)  !< shape: (n_dom)

  REAL(dp), SAVE :: cdimissval

CONTAINS


  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the physics state
  !!
  SUBROUTINE construct_aes_phy_memory( patch_array, ntracer )

    TYPE(t_patch),INTENT(IN) :: patch_array(:)
    INTEGER,INTENT(IN) :: ntracer
    CHARACTER(len=13) :: listname_f
    CHARACTER(len=12) :: listname_t
    INTEGER :: ndomain, jg, ist, nblks, nlev
    TYPE(timedelta) :: dt_dyn

    !---

    CALL message(thismodule,'Construction of AES physics state started.')

    cdimissval = cdiInqMissval()

    ! Allocate pointer arrays prm_field and prm_tend, 
    ! as well as the corresponding list arrays.

    ndomain = SIZE(patch_array)

    ALLOCATE( prm_field(ndomain), prm_tend(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      &'allocation of prm_field/tend array failed')

    ALLOCATE( prm_field_list(ndomain), prm_tend_list(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      &'allocation of prm_field/tend list array failed')

    ! Build a field list and a tendency list for each grid level.
    ! This includes memory allocation. 

    DO jg = 1,ndomain

      nblks = patch_array(jg)%nblks_c
      nlev  = patch_array(jg)%nlev
      dt_dyn = time_config%get_model_timestep_td(patch_array(jg)%nest_level)

      WRITE(listname_f,'(a,i2.2)') 'prm_field_D',jg
      CALL new_aes_phy_field_list( jg, nproma, nlev, nblks, ntracer, nsfc_type,   &
                                   & nnew(jg),                                    &
                                   & dt_dyn, listname_f, '',                      &
                                   & p_nh_state_lists(jg)%prog_list(nnew(jg)),    &
                                   & p_nh_state_lists(jg)%diag_list,              &
                                   & p_nh_state_lists(jg)%metrics_list,           &
                                   & ext_data(jg)%atm_list,                       &
                                   & prm_field_list(jg), prm_field(jg)            )

      WRITE(listname_t,'(a,i2.2)') 'prm_tend_D',jg
      CALL new_aes_phy_tend_list( jg, nproma, nlev, nblks, ntracer,   &
                                  & dt_dyn, listname_t, 'tend_',      &
                                  & prm_tend_list(jg), prm_tend(jg)   )
    ENDDO

    CALL message(thismodule,'Construction of AES physics state finished.')

  END SUBROUTINE construct_aes_phy_memory
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_aes_phy_memory

    INTEGER :: ndomain  !< total # of grid levels/domains
    INTEGER :: jg       !< grid level/domain index
    INTEGER :: ist      !< system status code

    !---
    CALL message(thismodule,'Destruction of AES physics state started.')

    ndomain = SIZE(prm_field)

    DO jg = 1,ndomain
      CALL vlr_del(prm_field_list(jg))
      CALL vlr_del(prm_tend_list (jg))
    ENDDO

    DEALLOCATE( prm_field_list, prm_tend_list, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      & 'deallocation of prm_field/tend list array failed')

    DEALLOCATE( prm_field, prm_tend, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      & 'deallocation of prm_field/tend array failed')

    CALL message(thismodule,'Destruction of AES physics state finished.')

  END SUBROUTINE destruct_aes_phy_memory
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  SUBROUTINE new_aes_phy_field_list  ( jg, kproma, klev, kblks, ktracer, ksfc_type, &
                                     & jt,                                          &
                                     & dt_dyn, listname, prefix,                    &
                                     & prog_list,                                   &
                                     & diag_list,                                   &
                                     & metrics_list,                                &
                                     & ext_atm_list,                                &
                                     & field_list, field                            )
    INTEGER,INTENT(IN) :: jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer, ksfc_type  !< dimension sizes
    INTEGER,INTENT(IN) :: jt                                       !< index for time level

    TYPE(timedelta),         INTENT(IN)    :: dt_dyn !< Dynamics timestep.

    CHARACTER(*),            INTENT(IN)    :: listname, prefix

    TYPE(t_var_list_ptr),    INTENT(INOUT) :: prog_list
    TYPE(t_var_list_ptr),    INTENT(INOUT) :: diag_list
    TYPE(t_var_list_ptr),    INTENT(INOUT) :: metrics_list
    TYPE(t_var_list_ptr),    INTENT(INOUT) :: ext_atm_list
    TYPE(t_var_list_ptr),    INTENT(INOUT) :: field_list
    TYPE(t_aes_phy_field),   INTENT(INOUT) :: field

    ! Local variables

    CHARACTER(LEN=vname_len) :: trc_name, cfstd_name, long_name, var_name, var_suffix
    CHARACTER(len=4) :: tl_suffix
    LOGICAL :: lclrsky_lw, lclrsky_sw
    LOGICAL :: contvar_is_in_output
    LOGICAL :: use_tmx

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shapesfc(3), shapeice(3), shape3d_layer_interfaces(3), shape3d_1level(3), shape3d_wrk(3)
    INTEGER :: ibits, iextbits, ivarbits
    INTEGER :: datatype_flt
    INTEGER :: jsfc, jtrc

    use_tmx = aes_vdf_config(jg)%use_tmx

    ibits = DATATYPE_PACK16
    iextbits = DATATYPE_PACK24

    ivarbits     = MERGE(DATATYPE_PACK24, DATATYPE_PACK16, gribout_config(jg)%lgribout_24bit)
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

    shape2d  = (/kproma,       kblks/)
    shape3d  = (/kproma, klev, kblks/)
    shapesfc = (/kproma, kblks, ksfc_type/)
    shape3d_layer_interfaces = (/kproma,klev+1,kblks/)
    shape3d_1level = (/kproma,1,kblks/)

    tl_suffix = get_timelevel_string(jt)

    !$ACC ENTER DATA COPYIN(field)
    ! Register a field list and apply default settings

    CALL vlr_add(field_list, listname, patch_id=jg, lrestart=.TRUE., &
      &          model_type=get_my_process_name())

    !------------------------------
    ! Metrics
    !------------------------------

    cf_desc    = t_cf_var('cell_longitude', 'rad',                              &
                &         'cell center longitude',                              &
                &         datatype_flt)
    grib2_desc = grib2_var(0,191,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'clon', field%clon,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE.  )
    __acc_attach(field%clon)

    cf_desc    = t_cf_var('cell_latitude', 'rad',                               &
                &         'cell center latitude',                               &
                &         datatype_flt)
    grib2_desc = grib2_var(0,191,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'clat', field%clat,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE.  )
    __acc_attach(field%clat)

    cf_desc    = t_cf_var('cell_area', 'm2',                                    &
                &         'Atmosphere Grid-Cell Area',                          &
                &         datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'areacella', field%areacella,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE.  )
    __acc_attach(field%areacella)

    cf_desc    = t_cf_var('height_above_reference_ellipsoid', 'm',             &
                &         'height above reference ellipsoid, half levels',     &
                &         datatype_flt)
    grib2_desc = grib2_var(0, 3, 6, ivarbits, GRID_UNSTRUCTURED, GRID_CELL)    &
                & + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_ref( metrics_list, 'z_ifc',                                       &
                & prefix//'zhalf', field%zh,                                   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF_HHL,               &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d_layer_interfaces,                   &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),    &
                & isteptype=TSTEP_CONSTANT                                     )
    __acc_attach(field%zh)

    cf_desc    = t_cf_var('height_above_reference_ellipsoid', 'm',             &
                &         'height above reference ellipsoid, full level',      &
                &         datatype_flt)
    grib2_desc = grib2_var(0,3,6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( metrics_list, 'z_mc',                                        &
                & prefix//'zfull', field%zf,                                   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN ),           &
                & isteptype=TSTEP_CONSTANT                                     )
    __acc_attach(field%zf)

    cf_desc    = t_cf_var('layer_thickness', 'm',                              &
                &         'layer thickness',                                   &
                &         datatype_flt)
    grib2_desc = grib2_var(0,3,12, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( metrics_list, 'ddqz_z_full',                                 &
                & prefix//'dzhalf', field%dz,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN ),           &
                & isteptype=TSTEP_CONSTANT                                     )
    __acc_attach(field%dz)


    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! &       field% ua        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'eastward wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( diag_list, 'u',                                              &
                & prefix//'ua', field%ua,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & lrestart = .FALSE.,                                          &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I") )      )
    __acc_attach(field%ua)

    ! &       field% va        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'northward wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( diag_list, 'v',                                              &
                & prefix//'va', field%va,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & lrestart = .FALSE.,                                          &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I") )      )
    __acc_attach(field%va)

    ! &       field% wa     (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('upward_air_velocity', 'm s-1', 'vertical velocity in m/s', datatype_flt)
    grib2_desc = grib2_var(0,2,9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( prog_list, 'w'//tl_suffix,                                       &
                & prefix//'wa_phy', field%wa,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                       &
                & cf_desc, grib2_desc,                                             &
                & ref_idx=1, ldims=shape3d_layer_interfaces,                       &
                & lrestart = .FALSE.,                                              &
                & vert_interp=create_vert_interp_metadata(                         &
                &             vert_intp_type=vintp_types("P","Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )
    ! Note: __acc_attach(field%<var>) must not be used here. The reason is that
    ! the pointer field%<var> is dynamic, ie. generally changes every time step.
    ! Therefore field%<var> needs to be attached after the host pointer has been
    ! changed to its new target, which is done by a subsequent directive
    ! ACC DATA PRESENT(field%<var>).

    ! &       field% ta        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('air_temperature', 'K', 'air temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( diag_list, 'temp',                                           &
                & prefix//'ta', field%ta,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & lrestart = .FALSE.,                                          &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN )            )
    ! Note: __acc_attach(field%<var>) must not be used here. The reason is that
    ! the pointer field%<var> is dynamic, ie. generally changes every time step.
    ! Therefore field%<var> needs to be attached after the host pointer has been
    ! changed to its new target, which is done by a subsequent directive
    ! ACC DATA PRESENT(field%<var>).

    ! &       field% tv        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature', datatype_flt)
    grib2_desc = grib2_var(0,0,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( diag_list, 'tempv',                                          &
                & prefix//'tv', field%tv,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & lrestart = .FALSE.,                                          &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN )            )
    __acc_attach(field%tv)

    ! average mass of snow flake
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('x_snow', 'kg', 'average mass of snow flake', datatype_flt)
    CALL add_var( field_list, prefix//'x_snow', field%x_snow,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%x_snow)

    ! average mass of ice crystal
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('x_ice', 'kg', 'average mass of ice crystal', datatype_flt)
    CALL add_var( field_list, prefix//'x_ice', field%x_ice,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%x_ice)

    ! cloud snow number concentration
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('acsnc', 'm^-3', 'cloud snow number concentration', datatype_flt)
    CALL add_var( field_list, prefix//'acsnc', field%acsnc,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%acsnc)

    ! cloud ice number concentration
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('acinc', 'm^-3', 'cloud ice number concentration', datatype_flt)
    CALL add_var( field_list, prefix//'acinc', field%acinc,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%acinc)

    ! effective radius of ice
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('reff_ice', 'um', 'effective radius of cloud ice', datatype_flt)
    CALL add_var( field_list, prefix//'reff_ice', field%reff_ice,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%reff_ice)

    ! optical depth of ice
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('tau_ice', '', 'optical depth of ice, integral over all bands', datatype_flt)
    CALL add_var( field_list, prefix//'tau_ice', field%tau_ice,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%tau_ice)
    
    ! effective radius of snow
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('reff_snow', 'um', 'effective radius of snow', datatype_flt)
    CALL add_var( field_list, prefix//'reff_snow', field%reff_snow,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%reff_snow)

    ! optical depth of snow
    grib2_desc = grib2_var(0,4,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('tau_snow', '', 'optical depth of snow, integral over all bands', datatype_flt)
    CALL add_var( field_list, prefix//'tau_snow', field%tau_snow,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%tau_snow)
    
    ! OZONE 
    ! &       field% o3        (nproma,nlev  ,nblks),          &
    grib2_desc = grib2_var(0,14,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('mole_fraction_of_ozone_in_air', 'mol/mol', 'ozone volume mixing ratio', datatype_flt)
    CALL add_var( field_list, prefix//'xo3', field%o3,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lopenacc=.TRUE.)
    __acc_attach(field%o3)

    ! aerosol optical properties
    ! at 533 nm
    cf_desc    = t_cf_var('aer_aod_533','-','aerosol optical depth at 533 nm', &
                & datatype_flt)
    grib2_desc = grib2_var(0,20,102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_533', field%aer_aod_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ),                                      &
                & lopenacc=.TRUE.)
    __acc_attach(field%aer_aod_533)
    cf_desc    = t_cf_var('aer_ssa_533','-',                                   &
                & 'aerosol single scattering albedo at 533 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,103, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_533', field%aer_ssa_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ),                                      &
                & lopenacc=.TRUE.)
    __acc_attach(field%aer_ssa_533)
    cf_desc    = t_cf_var('aer_asy_533','-',                                   &
                & 'aerosol asymmetry factor at 533 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,104, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_533', field%aer_asy_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN ),                       &
                & lopenacc=.TRUE.)
    __acc_attach(field%aer_asy_533)
    ! at 2325 nm
    cf_desc    = t_cf_var('aer_aod_2325','-',                                  &
                & 'aerosol optical depth at 2325 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_2325', field%aer_aod_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ),                                      &
                & lopenacc=.TRUE.)
    __acc_attach(field%aer_aod_2325)
    cf_desc    = t_cf_var('aer_ssa_2325','-',                                  &
                & 'aerosol single scattering albedo at 2325 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,103, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_2325', field%aer_ssa_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ),                                      &
                & lopenacc=.TRUE.)
    __acc_attach(field%aer_ssa_2325)
    cf_desc    = t_cf_var('aer_asy_2325','-',                                  &
                & 'aerosol asymmetry factor at 2325 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,104, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_2325', field%aer_asy_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ),                                      &
                & lopenacc=.TRUE.)
    __acc_attach(field%aer_asy_2325)
    ! at 9731 nm
    cf_desc    = t_cf_var('aer_aod_9731','-',                                  &
                & 'effective aerosol optical depth at 9731 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_9731', field%aer_aod_9731,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ),                                      &
                & lopenacc=.TRUE.)
    __acc_attach(field%aer_aod_9731)

    !--------
    ! Tracers
    !--------

    IF (ktracer > 0) THEN
      !
      ! Reference to array for mass fractions of tracers in air, for computations in dynamics
      !
      ! &       field% qtrc_dyn  (nproma,nlev  ,nblks,ntracer),  &
      CALL add_ref( prog_list, 'tracer'//tl_suffix,                              &
                  & prefix//'qtrc_dyn', field%qtrc_dyn,                          &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & t_cf_var('mass_fraction_of_tracer_in_air', 'kg kg-1',        &
                  &          'mass fraction of tracer in air (dynamics)',        &
                  &          datatype_flt),                                      &
                  & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                  & ref_idx=1, ldims = (/kproma,klev,kblks,ktracer/),            &
                  & lrestart=.FALSE., loutput=.FALSE.)
      ! Note: __acc_attach(field%<var>) must not be used here. The reason is that
      ! the pointer field%<var> is dynamic, ie. generally changes every time step.
      ! Therefore field%<var> needs to be attached after the host pointer has been
      ! changed to its new target, which is done by a subsequent directive
      ! ACC DATA PRESENT(field%<var>).

      ! Reference to array for mass fractions of tracers in air, for computations in physics
      !
      CALL add_ref( prog_list, 'tracer'//tl_suffix,                              &
                  & prefix//'qtrc_phy', field%qtrc_phy,                          &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & t_cf_var('mass_fraction_of_tracer_in_air', 'kg kg-1',        &
                  &          'mass fraction of tracer in air (physics)',         &
                  &          datatype_flt),                                      &
                  & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                  & ref_idx=1, ldims = (/kproma,klev,kblks,ktracer/),            &
                  & lrestart=.FALSE., loutput=.FALSE.)
      ! Note: __acc_attach(field%<var>) must not be used here. The reason is that
      ! the pointer field%<var> is dynamic, ie. generally changes every time step.
      ! Therefore field%<var> needs to be attached after the host pointer has been
      ! changed to its new target, which is done by a subsequent directive
      ! ACC DATA PRESENT(field%<var>).

      ! Reference to array for tracer paths, for computations
      !
      CALL add_ref( diag_list, 'tracer_vi',                                      &
                  & prefix//'mtrcvi', field%mtrcvi,                              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('atmosphere_mass_content_of_tracers',               &
                  &          'kg m-2',                                           &
                  &          'tracer paths',                                     &
                  &          datatype_flt),                                      &
                  & grib2_var(0,20,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                  & ref_idx=1, ldims = (/kproma,kblks,ktracer/),                 &
                  & lrestart=.FALSE., loutput=.FALSE.)
      __acc_attach(field%mtrcvi)

      ! References for single tracer paths, for output
      !
      ALLOCATE(field%mtrcvi_ptr(ktracer))
      !
      DO jtrc = 1,ktracer
        !
        field%mtrcvi_ptr(jtrc)%p => NULL()
        !
        ! preset names and grib2_desc for all tracers
        !
        trc_name   = TRIM(advection_config(jg)%tracer_names(jtrc))
        cfstd_name = TRIM(advection_config(jg)% cfstd_names(jtrc))
        long_name  = TRIM(advection_config(jg)%  long_names(jtrc))
        !
        var_suffix = 'vi'
        grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        !
        ! adjust names and grib2_desc for specific tracers 
        !
        IF (jtrc == iqv) THEN
          trc_name   = 'prw'
          var_suffix = ''
          cfstd_name = 'water_vapor'
          long_name  = 'water vapor'
          grib2_desc = grib2_var(0, 1, 64, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
        ELSE IF (jtrc == iqc) THEN
          trc_name   = 'cll'
          grib2_desc = grib2_var(0, 1, 69, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
        ELSE IF (jtrc == iqi) THEN
          grib2_desc = grib2_var(0, 1, 70, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
        ELSE IF (jtrc == iqr) THEN
          grib2_desc = grib2_var(0, 1, 45, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
        ELSE IF (jtrc == iqs) THEN
          grib2_desc = grib2_var(0, 1, 46, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
        ELSE IF (jtrc == iqg) THEN
          grib2_desc = grib2_var(0, 1, 74, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
        ELSE IF (jtrc == iqh) THEN
          grib2_desc = grib2_var(0, 1, 72, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
        END IF
        !
        ! set var_name and cf_desc
        !
        var_name = prefix//TRIM(trc_name)//TRIM(var_suffix)
        cf_desc = t_cf_var( 'atmosphere_mass_content_of_'//TRIM(cfstd_name), &
                          & 'kg m-2',                                        &
                          & TRIM(long_name)//' path',                        &
                          & datatype_flt )
        !
        ! set memory references for fields which are requested for output
        !
        IF ( is_variable_in_output(var_name=TRIM(var_name)) ) THEN
          CALL add_ref( diag_list, 'tracer_vi',                   &
                      & TRIM(var_name), field%mtrcvi_ptr(jtrc)%p, &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,       &
                      & cf_desc, grib2_desc,                      &
                      & ref_idx=jtrc, ldims=(/kproma,kblks/),     &
                      & lrestart = .FALSE. )
        END IF
        !
      END DO
      !
    END IF ! (ktracer > 0)
 
    ! &       field% rho        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('air_density', 'kg m-3', 'density of air',           &
         &                datatype_flt)
    grib2_desc = grib2_var(0,3,10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( prog_list, 'rho'//tl_suffix,                                 &
         &        prefix//'rho_phy', field%rho,                                &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
         &        cf_desc, grib2_desc,                                         &
         &        ref_idx=1, ldims=shape3d,                                    &
         &        lrestart = .FALSE.,                                          &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN ) )
    ! Note: __acc_attach(field%<var>) must not be used here. The reason is that
    ! the pointer field%<var> is dynamic, ie. generally changes every time step.
    ! Therefore field%<var> needs to be attached after the host pointer has been
    ! changed to its new target, which is done by a subsequent directive
    ! ACC DATA PRESENT(field%<var>).

    ! &       field% mair        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('air_mass', 'kg m-2', 'air mass in layer', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( diag_list, 'airmass_new',                                    &
         &        prefix//'mair_phy', field%mair,                              &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
         &        cf_desc, grib2_desc,                                         &
         &        ref_idx=1, ldims=shape3d,                                    &
         &        lrestart = .FALSE., initval=1.0_wp,                          &
         &        vert_interp = create_vert_interp_metadata(                   &
         &                      vert_intp_type=vintp_types("P","Z","I"),       &
         &                      vert_intp_method=VINTP_METHOD_LIN,             &
         &                      l_loglin=.FALSE.,                              &
         &                      l_extrapol=.TRUE., l_pd_limit=.FALSE.,         &
         &                      lower_limit=0._wp ) )
    __acc_attach(field%mair)

    IF (is_variable_in_output(var_name=prefix//'pv')) THEN
       cf_desc    = t_cf_var('potential_vorticity', 'K m2 kg-1 s-1', 'potential vorticity', datatype_flt)
       grib2_desc = grib2_var(0, 2, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list,                                                  &
                   & prefix//"pv", field%pv,                                      &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                   & cf_desc, grib2_desc,                                         &
                   & ldims=shape3d,                                               &
                   & vert_interp=create_vert_interp_metadata(                     &
                   &             vert_intp_type=vintp_types("P","Z","I"),         &
                   &             vert_intp_method=VINTP_METHOD_LIN,               &
                   &             l_loglin=.FALSE.,                                &
                   &             l_extrapol=.FALSE.),                             &
                   & l_pp_scheduler_task=TASK_COMPUTE_PV, lrestart=.FALSE.,       &
                   & lopenacc=.TRUE.)
        __acc_attach(field%pv)
    END IF

    ! &       field% geom      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('geopotential_above_surface', 'm2 s-2', 'geopotential above surface', datatype_flt)
    grib2_desc = grib2_var(0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( metrics_list, 'geopot_agl',                                  &
                & prefix//'gpsm', field%geom,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN,             &
                &               l_extrapol=.TRUE., l_pd_limit=.FALSE.)         )
    __acc_attach(field%geom)

    ! &       field% geop      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential', datatype_flt)
    grib2_desc = grib2_var(0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( metrics_list, 'geopot',                                      &
                & prefix//'geop', field%geop,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN,             &
                &               l_extrapol=.TRUE., l_pd_limit=.FALSE.)         )
    __acc_attach(field%geom)

    ! &       field% pfull (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('air_pressure', 'Pa', 'air pressure on model levels', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( diag_list, 'pres',                                           &
                & prefix//'pfull', field%pfull,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d,                                    &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &             vert_intp_type=vintp_types("Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_PRES )             )
    __acc_attach(field%pfull)

    !-- Variables defined at layer interfaces --

    ! &       field% geoi      (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('geopotential above surface', 'm2 s-2', 'geopotential above surface', datatype_flt)
    grib2_desc = grib2_var(0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( metrics_list, 'geopot_agl_ifc',                              &
                & prefix//'gpsi', field%geoi,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                   &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d_layer_interfaces,                   &
                & vert_interp = create_vert_interp_metadata(                   &
                &               vert_intp_type=vintp_types("P","Z","I"),       &
                &               vert_intp_method=VINTP_METHOD_LIN_NLEVP1 )     )
    __acc_attach(field%geoi)

    ! &       field% phalf (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('air_pressure', 'Pa', 'air pressure on model half levels', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( diag_list, 'pres_ifc',                                       &
                & prefix//'phalf', field%phalf,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                   &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape3d_layer_interfaces,                   &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &             vert_intp_type=vintp_types("Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 )       )
    __acc_attach(field%phalf)

    !------------------
    ! Radiation
    !------------------
    !

    cf_desc    = t_cf_var( 'cosmu0'                                      , &
         &                 ''                                            , &
         &                 'cosine of the zenith angle for rad. heating' , &
         &                 datatype_flt                                  )
    grib2_desc = grib2_var(0,4,214, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'cosmu0' , field%cosmu0, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE        , &
         &       cf_desc , grib2_desc                       , &
         &       lrestart = .FALSE.                         , &
         &       ldims=shape2d                              , &
         &       lopenacc=.TRUE.                            )
    __acc_attach(field%cosmu0)

    cf_desc    = t_cf_var( 'cosmu0_rt'                                    , &
         &                 ''                                             , &
         &                 'cosine of the zenith angle at radiation time' , &
         &                 datatype_flt                                   )
    grib2_desc = grib2_var(192,214,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'cosmu0_rt', field%cosmu0_rt, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE             , &
         &       cf_desc , grib2_desc                            , &
         &       lrestart = .TRUE.                               , &
         &       ldims=shape2d                                   , &
         &       lopenacc=.TRUE.                                 )
    __acc_attach(field%cosmu0_rt)

    cf_desc    = t_cf_var( 'daylght_frc', &
         &                 ''           , &
         &                 ''           , &
         &                 datatype_flt )
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'daylght_frc', field%daylght_frc, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE                 , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.                                  , &
         &       ldims=shape2d                                       , &
         &       lopenacc=.TRUE.                                     )
    __acc_attach(field%daylght_frc)

    cf_desc    = t_cf_var( 'daylght_frc_rt' , &
         &                 ''               , &
         &                 ''               , &
         &                 datatype_flt     )
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'daylght_frc_rt', field%daylght_frc_rt, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE                       , &
         &       cf_desc, grib2_desc                                       , &
         &       lrestart = .TRUE.                                         , &
         &       ldims=shape2d                                             , &
         &       lopenacc=.TRUE.                                           )
    __acc_attach(field%daylght_frc_rt)

    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       ! shortwave fluxes
       !
       ! - flag for clear sky computations
       lclrsky_sw = is_variable_in_output(var_name=prefix//'rsdcs')  .OR. &
            &       is_variable_in_output(var_name=prefix//'rsucs')  .OR. &
            &       is_variable_in_output(var_name=prefix//'rsutcs') .OR. &
            &       is_variable_in_output(var_name=prefix//'rsdscs') .OR. &
            &       is_variable_in_output(var_name=prefix//'rsuscs')
       aes_rad_config(jg)%lclrsky_sw = lclrsky_sw
       !
       ! - allocate clear sky 3d radiation fields with a single level only if
       !   they are not used for computations, but still needed as arguments
       IF (lclrsky_sw) THEN
          shape3d_wrk = shape3d_layer_interfaces
       ELSE
          shape3d_wrk = shape3d_1level
       END IF
       !
       ! - through the atmosphere
       !
       cf_desc    = t_cf_var('downwelling_shortwave_flux_in_air', &
            &                'W m-2'                            , &
            &                'downwelling shortwave radiation'  , &
            &                datatype_flt                       )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsd', field%rsd_rt     , &
            &       GRID_UNSTRUCTURED_CELL   , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape3d_layer_interfaces              , &
            &       vert_interp=create_vert_interp_metadata       &
            &         (vert_intp_type=vintp_types("P","Z","I") ,  &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1), &
            &       lopenacc=.TRUE.)
       __acc_attach(field%rsd_rt     )

       cf_desc    = t_cf_var('upwelling_shortwave_flux_in_air', &
            &                'W m-2'                          , &
            &                'upwelling shortwave radiation'  , &
            &                datatype_flt                     )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsu' , field%rsu_rt    , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF  , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape3d_layer_interfaces              , &
            &       vert_interp=create_vert_interp_metadata       &
            &         (vert_intp_type=vintp_types("P","Z","I") ,  &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1), &
            &       lopenacc=.TRUE.)
       __acc_attach(field%rsu_rt    )

       cf_desc    = t_cf_var('downwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                               , &
            &                'downwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                          )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsdcs' , field%rsdcs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_wrk                            , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1),  &
            &       lopenacc=.TRUE.)
       __acc_attach(field%rsdcs_rt )

       cf_desc    = t_cf_var('upwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                             , &
            &                'upwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                        )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsucs' , field%rsucs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_wrk                            , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1),  &
            &       lopenacc=.TRUE.)
       __acc_attach(field%rsucs_rt )

       ! - at the top of the atmosphere
       !
       cf_desc    = t_cf_var('toa_incoming_shortwave_flux'     , &
            &                'W m-2'                           , &
            &                'toa incident shortwave radiation', &
            &                datatype_flt                      )
       grib2_desc = grib2_var(0,4,201, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsdt', field%rsdt, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_TOA    , &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         , &
            &       lopenacc=.TRUE.                       )
       __acc_attach(field%rsdt)

       cf_desc    = t_cf_var('toa_outgoing_shortwave_flux'     , &
            &                'W m-2'                           , &
            &                'toa outgoing shortwave radiation', &
            &                datatype_flt                      )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsut', field%rsut, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_TOA    , &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         , &
            &       lopenacc=.TRUE.                       )
       __acc_attach(field%rsut)

       cf_desc    = t_cf_var('toa_outgoing_shortwave_flux_assuming_clear_sky', &
            &                'W m-2'                                         , &
            &                'toa outgoing clear-sky shortwave radiation'    , &
            &                datatype_flt                                    )
       grib2_desc = grib2_var(0,4,208, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsutcs', field%rsutcs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_TOA      , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             , &
            &       lopenacc=.TRUE.                           )
       __acc_attach(field%rsutcs)

    END IF

    ! - at the surface (also used in update surface)
    !
    cf_desc    = t_cf_var('surface_downwelling_shortwave_flux_in_air', &
         &                'W m-2'                                    , &
         &                'surface downwelling shortwave radiation'  , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rsds', field%rsds, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rsds)

    cf_desc    = t_cf_var('surface_upwelling_shortwave_flux_in_air', &
         &                'W m-2'                                  , &
         &                'surface upwelling shortwave radiation'  , &
         &                datatype_flt                             )
    grib2_desc = grib2_var(0,4,199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rsus', field%rsus, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rsus)


    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('surface_downwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                                       , &
            &                'surface downwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                                  )
       grib2_desc = grib2_var(0,4,207, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsdscs', field%rsdscs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             , &
            &       lopenacc=.TRUE.                           )
       __acc_attach(field%rsdscs)

       cf_desc    = t_cf_var('surface_upwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                                     , &
            &                'surface upwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                                )
       grib2_desc = grib2_var(0,4,209, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsuscs', field%rsuscs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             , &
            &       lopenacc=.TRUE.                           )
       __acc_attach(field%rsuscs)

       !-----------------------------------------------------------------------------------
       ! shortwave flux components at the surface
       ! - at radiation times
       cf_desc    = t_cf_var('surface_downwelling_direct_visible_flux_in_air_at_rad_time'    , &
            &                'W m-2'                                                         , &
            &                'surface downwelling direct visible radiation at radiation time', &
            &                datatype_flt                                                    )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvds_dir_rt', field%rvds_dir_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       , &
            &       lopenacc=.TRUE.                                     )
       __acc_attach(field%rvds_dir_rt)

       cf_desc    = t_cf_var('surface_downwelling_direct_par_flux_in_air_at_rad_time'                          , &
            &                'W m-2'                                                                           , &
            &                'surface downwelling direct photosynthetically active radiation at radiation time', &
            &                datatype_flt                                                                      )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpds_dir_rt', field%rpds_dir_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       , &
            &       lopenacc=.TRUE.                                     )
       __acc_attach(field%rpds_dir_rt)

       cf_desc    = t_cf_var('surface_downwelling_direct_nearir_flux_in_air_at_rad_time'           , &
            &                'W m-2'                                                               , &
            &                'surface downwelling direct near infrared radiation at radiation time', &
            &                datatype_flt                                                          )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnds_dir_rt', field%rnds_dir_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       , &
            &       lopenacc=.TRUE.                                     )
       __acc_attach(field%rnds_dir_rt)


       cf_desc    = t_cf_var('surface_downwelling_diffuse_visible_flux_in_air_at_rad_time'    , &
            &                'W m-2'                                                          , &
            &                'surface downwelling diffuse visible radiation at radiation time', &
            &                datatype_flt                                                     )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvds_dif_rt', field%rvds_dif_rt, &
            &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE          , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       , &
            &       lopenacc=.TRUE.                                     )
       __acc_attach(field%rvds_dif_rt)

       cf_desc    = t_cf_var('surface_downwelling_diffuse_par_flux_in_air_at_rad_time'                          , &
            &                'W m-2'                                                                            , &
            &                'surface downwelling diffuse photosynthetically active radiation at radiation time', &
            &                datatype_flt                                                                       )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpds_dif_rt', field%rpds_dif_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       , &
            &       lopenacc=.TRUE.                                     )
       __acc_attach(field%rpds_dif_rt)

       cf_desc    = t_cf_var('surface_downwelling_diffuse_nearir_flux_in_air_at_rad_time'           , &
            &                'W m-2'                                                                , &
            &                'surface downwelling diffuse near infrared radiation at radiation time', &
            &                datatype_flt                                                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnds_dif_rt', field%rnds_dif_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       , &
            &       lopenacc=.TRUE.                                     )
       __acc_attach(field%rnds_dif_rt)


       cf_desc    = t_cf_var('surface_upwelling_visible_flux_in_air_at_rad_time'    , &
            &                'W m-2'                                                , &
            &                'surface upwelling visible radiation at radiation time', &
            &                datatype_flt                                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvus_rt', field%rvus_rt, &
            &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape2d                               , &
            &       lopenacc=.TRUE.                             )
       __acc_attach(field%rvus_rt)

       cf_desc    = t_cf_var('surface_upwelling_par_flux_in_air_at_rad_time'                          , &
            &                'W m-2'                                                                  , &
            &                'surface upwelling photosynthetically active radiation at radiation time', &
            &                datatype_flt                                                             )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpus_rt', field%rpus_rt, &
            &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape2d                               , &
            &       lopenacc=.TRUE.                             )
       __acc_attach(field%rpus_rt)

       cf_desc    = t_cf_var('surface_upwelling_nearir_flux_in_air_at_rad_time'           , &
            &                'W m-2'                                                      , &
            &                'surface upwelling near infrared radiation at radiation time', &
            &                datatype_flt                                                 )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnus_rt', field%rnus_rt, &
            &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape2d                               , &
            &       lopenacc=.TRUE.                             )
       __acc_attach(field%rnus_rt)

    END IF

    !-----------------------------------------------------------------------------------------
    ! shortwave flux components at the surface (also used in update_surface)
    ! - at all times
    cf_desc    = t_cf_var('surface_downwelling_direct_visible_flux_in_air', &
         &                'W m-2'                                         , &
         &                'surface downwelling direct visible radiation'  , &
         &                datatype_flt                                    )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rvds_dir', field%rvds_dir, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 , &
         &       lopenacc=.TRUE.                               )
    __acc_attach(field%rvds_dir)

    cf_desc    = t_cf_var('surface_downwelling_direct_par_flux_in_air'                    , &
         &                'W m-2'                                                         , &
         &                'surface downwelling direct photosynthetically active radiation', &
         &                datatype_flt                                                    )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rpds_dir', field%rpds_dir, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 , &
         &       lopenacc=.TRUE.                               )
    __acc_attach(field%rpds_dir)

    cf_desc    = t_cf_var('surface_downwelling_direct_nearir_flux_in_air'     , &
         &                'W m-2'                                             , &
         &                'surface downwelling direct near infrared radiation', &
         &                datatype_flt                                        )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rnds_dir', field%rnds_dir, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 , &
         &       lopenacc=.TRUE.                               )
    __acc_attach(field%rnds_dir)


    cf_desc    = t_cf_var('surface_downwelling_diffuse_visible_flux_in_air', &
         &                'W m-2'                                          , &
         &                'surface downwelling diffuse visible radiation'  , &
         &                datatype_flt                                     )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rvds_dif', field%rvds_dif, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 , &
         &       lopenacc=.TRUE.                               )
    __acc_attach(field%rvds_dif)

    cf_desc    = t_cf_var('surface_downwelling_diffuse_par_flux_in_air'                    , &
         &                'W m-2'                                                          , &
         &                'surface downwelling diffuse photosynthetically active radiation', &
         &                datatype_flt                                                     )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rpds_dif', field%rpds_dif, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 , &
         &       lopenacc=.TRUE.                               )
    __acc_attach(field%rpds_dif)

    cf_desc    = t_cf_var('surface_downwelling_diffuse_nearir_flux_in_air'     , &
         &                'W m-2'                                              , &
         &                'surface downwelling diffuse near infrared radiation', &
         &                datatype_flt                                         )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rnds_dif', field%rnds_dif, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 , &
         &       lopenacc=.TRUE.                               )
    __acc_attach(field%rnds_dif)


    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('surface_upwelling_visible_flux_in_air', &
            &                'W m-2'                                , &
            &                'surface upwelling visible radiation'  , &
            &                datatype_flt                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvus', field%rvus, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         , &
            &       lopenacc=.TRUE.                       )
       __acc_attach(field%rvus)

       cf_desc    = t_cf_var('surface_upwelling_par_flux_in_air'                    , &
            &                'W m-2'                                                , &
            &                'surface upwelling photosynthetically active radiation', &
            &                datatype_flt                                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpus', field%rpus, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         , &
            &       lopenacc=.TRUE.                       )
       __acc_attach(field%rpus)

       cf_desc    = t_cf_var('surface_upwelling_nearir_flux_in_air'     , &
            &                'W m-2'                                    , &
            &                'surface upwelling near infrared radiation', &
            &                datatype_flt                               )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnus', field%rnus, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         , &
            &       lopenacc=.TRUE.                       )
       __acc_attach(field%rnus)
       !---------------------------------------------------------

    END IF
    !
    ! net shortwave fluxes only needed for diagnostic output for destinE
    ! 176: rsns: surface net shortwave radiation flux: rsds - rsus  0-4-9-ffs1-sp1
    ! 178: rsnt: TOA net shortwave radiation flux:     rsdt - rsut  0-4-9-ffs8-sp1

    cf_desc    = t_cf_var('surface_net_shortwave_radiation_flux_in_air', &
         &                'W m-2'                                      , &
         &                'surface net shortwave radiation flux'       , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,4,9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rsns', field%rsns, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rsns)

    cf_desc    = t_cf_var('toa_net_shortwave_radiation_flux_in_air', &
         &                'W m-2'                                      , &
         &                'toa net shortwave radiation flux'       , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,4,9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rsnt', field%rsnt, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_TOA    , &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rsnt)

    !------------------
    !

    ! longwave  fluxes
    !
    ! - flag for clear sky computations
    lclrsky_lw = is_variable_in_output(var_name=prefix//'rldcs')  .OR. &
         &       is_variable_in_output(var_name=prefix//'rlucs')  .OR. &
         &       is_variable_in_output(var_name=prefix//'rlutcs') .OR. &
         &       is_variable_in_output(var_name=prefix//'rldscs')
    aes_rad_config(jg)%lclrsky_lw = lclrsky_lw
    !
    ! - allocate clear sky 3d radiation fields with a single level only if
    !   they are not used for computations, but still needed as arguments
    IF (lclrsky_lw) THEN
       shape3d_wrk = shape3d_layer_interfaces
    ELSE
       shape3d_wrk = shape3d_1level
    END IF
    !
    ! - through the atmosphere
    !
    ! (rld_rt and rlu_rt are also needed in update_surface)

    cf_desc    = t_cf_var('downwelling_longwave_flux_in_air', &
         &                'W m-2'                           , &
         &                'downwelling longwave radiation'  , &
         &                datatype_flt                       )
    grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rld' , field%rld_rt     , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .TRUE.                            , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1),  &
         &        lopenacc=.TRUE.)
    __acc_attach(field%rld_rt     )

    cf_desc    = t_cf_var('upwelling_longwave_flux_in_air', &
         &                'W m-2'                         , &
         &                'upwelling longwave radiation'  , &
         &                datatype_flt                     )
    grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlu' , field%rlu_rt     , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .TRUE.                            , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1),  &
         &        lopenacc=.TRUE.)
    __acc_attach(field%rlu_rt     )

    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('downwelling_longwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                              , &
            &                'downwelling clear-sky longwave radiation'           , &
            &                datatype_flt                                         )
       grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rldcs' , field%rldcs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_wrk                            , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1),  &
            &       lopenacc=.TRUE.)
       __acc_attach(field%rldcs_rt )

       cf_desc    = t_cf_var('upwelling_longwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                            , &
            &                'upwelling clear-sky longwave radiation'           , &
            &                datatype_flt                                       )
       grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rlucs' , field%rlucs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_wrk                            , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1),  &
            &       lopenacc=.TRUE.)
       __acc_attach(field%rlucs_rt )

       ! - at the top of the atmosphere
       !
       cf_desc    = t_cf_var('toa_outgoing_longwave_flux'     , &
            &                'W m-2'                          , &
            &                'toa outgoing longwave radiation', &
            &                datatype_flt                     )
       grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rlut', field%rlut, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_TOA    , &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         , &
            &       lopenacc=.TRUE.                       )
       __acc_attach(field%rlut)

       cf_desc    = t_cf_var('toa_outgoing_longwave_flux_assuming_clear_sky', &
            &                'W m-2'                                        , &
            &                'toa outgoing clear-sky longwave radiation'    , &
            &                datatype_flt                                   )
       grib2_desc = grib2_var(0,5,204, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rlutcs', field%rlutcs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_TOA      , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             , &
            &       lopenacc=.TRUE.                           )
       __acc_attach(field%rlutcs)

    END IF

    ! - at the surface (also used in update_surface)
    !
    cf_desc    = t_cf_var('surface_downwelling_longwave_flux_in_air', &
         &                'W m-2'                                   , &
         &                'surface downwelling longwave radiation'  , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlds', field%rlds, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rlds)

    cf_desc    = t_cf_var('surface_upwelling_longwave_flux_in_air', &
         &                'W m-2'                                 , &
         &                'surface upwelling longwave radiation'  , &
         &                datatype_flt                             )
    grib2_desc = grib2_var(0,5,199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlus', field%rlus, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rlus)

    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('surface_downwelling_longwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                                      , &
            &                'surface downwelling clear-sky longwave radiation'           , &
            &                datatype_flt                                                 )
       grib2_desc = grib2_var(0,5,203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rldscs', field%rldscs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             , &
            &       lopenacc=.TRUE.                           )
       __acc_attach(field%rldscs)
       !
    END IF

    ! net longwave fluxes only needed for diagnostic output for destinE
    !177: rlns: surface net longwave radiation flux:  rlds - rlus  0-5-5-ffs1-sp1
    !179: rlnt: TOA net longwave radiation flux:           - rlut  0-5-5-ffs8-sp1
    
    cf_desc    = t_cf_var('surface_net_longwave_radiation_flux_in_air', &
         &                'W m-2'                                     , &
         &                'surface net longwave radiation flux'       , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,5,5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlns', field%rlns, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rlns)

    cf_desc    = t_cf_var('toa_net_longwave_radiation_flux_in_air', &
         &                'W m-2'                                 , &
         &                'toa net longwave radiation flux'       , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,5,5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlnt', field%rlnt, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_TOA    , &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         , &
         &       lopenacc=.TRUE.                       )
    __acc_attach(field%rlnt)

    !
    !------------------
    !

    IF (.NOT. use_tmx) THEN
      cf_desc    = t_cf_var('drlns_dT', 'W m-2 K-1', 'longwave net flux T-derivative at surface', &
          &                datatype_flt)
      grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'drlns_dT', field%dlwflxsfc_dT, &
          &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                 &
          &        cf_desc, grib2_desc,                                &
          &        ldims=shape2d,                                      &
          &        lrestart = .FALSE.,                                 &
          &        isteptype=TSTEP_INSTANT,                            &
          &        lopenacc=.TRUE.)
      __acc_attach(field%dlwflxsfc_dT)
    END IF

    cf_desc    = t_cf_var('siced', 'm', 'sea ice thickness', datatype_flt)
    grib2_desc = grib2_var(10,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sit', field%siced,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., ldims=shape2d,                        &
                & lopenacc=.TRUE. )
    __acc_attach(field%siced)

    cf_desc    = t_cf_var('alb', '', 'surface albedo from external file', datatype_flt)
    grib2_desc = grib2_var(0,19,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'alb', field%alb,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & lopenacc=.TRUE. )
    __acc_attach(field%alb)

    cf_desc    = t_cf_var('ts_rad', 'K', 'radiative surface temperature', datatype_flt)
    grib2_desc = grib2_var(0,0,17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ts_rad', field%ts_rad,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., ldims=shape2d,                        &
                & lopenacc=.TRUE. )
    __acc_attach(field%ts_rad)

    cf_desc    = t_cf_var('ts_rad_rt', 'K', 'radiative surface temperature at rad. time', datatype_flt)
    grib2_desc = grib2_var(0,0,17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ts_rad_rt', field%ts_rad_rt,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., ldims=shape2d,                        &
                & lopenacc=.TRUE. )
    __acc_attach(field%ts_rad_rt)

    IF (aes_phy_tc(jg)%dt_vdf > dt_dyn .OR.                                &
      & is_variable_in_output(var_name=prefix//'q_snocpymlt')) THEN
       cf_desc    = t_cf_var('q_snocpymlt', 'W/m2', 'heating for snow melt on canopy', datatype_flt)
       grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'q_snocpymlt', field%q_snocpymlt,    &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & lopenacc=.TRUE. )
       __acc_attach(field%q_snocpymlt)
    END IF

    IF (is_variable_in_output(var_name=prefix//'q_rlw_impl')) THEN
       cf_desc    = t_cf_var('q_rlw_impl', 'W/m2', 'heating correction due to implicit land surface coupling', datatype_flt)
       grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'q_rlw_impl', field%q_rlw_impl,    &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                     &
                   & lopenacc=.TRUE. )
       __acc_attach(field%q_rlw_impl)
    END IF

    cf_desc    = t_cf_var('q_rlw_nlev', 'W/m2', 'LW heating in the lowest layer', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'q_rlw_nlev', field%q_rlw_nlev,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                     &
                & lopenacc=.TRUE. )
    __acc_attach(field%q_rlw_nlev)
    !
    !------------------
    !
    ! CO2

    cf_desc = t_cf_var('fco2nat', 'kg m-2 s-1',                                &
                & 'Surface Carbon Mass Flux into the Atmosphere Due to Natural Sources', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)

    CALL add_var( field_list, prefix//'fco2nat', field%fco2nat,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,     &
                & lrestart = .TRUE., initval =  0.0_wp, ldims=shape2d,         &
                & lopenacc=.TRUE. )

    __acc_attach(field%fco2nat)

    ! &       field% co2_flux_tile(nproma,nblks,nsfc_type), &
    CALL add_var( field_list, prefix//'co2_flux_tile', field%co2_flux_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('co2_flux_tile',  'kg m-2 s-1',                     &
                & 'surface_upward_mass_flux_of_carbon_dioxide', datatype_flt), &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                & ldims=shapesfc, initval=0.0_wp,                              &
                & lcontainer=.TRUE., lrestart=.FALSE.,                         &
                & lopenacc=.TRUE. )
    __acc_attach(field%co2_flux_tile)

    ALLOCATE(field%co2_flux_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'co2_flux_tile',                         &
                  & prefix//'co2_flux_'//csfc(jsfc), field%co2_flux_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('co2_flux_'//csfc(jsfc), 'kg m-2 s-1',              &
                  & 'surface_upward_mass_flux_of_carbon_dioxide', datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.TRUE., initval=0.0_wp,&
                  & lmiss=.TRUE., missval=cdimissval )

    END DO

    !
    !------------------
    !

    cf_desc    = t_cf_var('csat', '', '', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'csat', field%csat,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., initval =  1.0_wp, ldims=shape2d,     &
                & lopenacc=.TRUE. )
    __acc_attach(field%csat)

    cf_desc    = t_cf_var('cair', '', '', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'cair', field%cair,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., initval =  1.0_wp, ldims=shape2d,     &
                & lopenacc=.TRUE. )
    __acc_attach(field%cair)

    !-------------------------
    ! Sea ice
    !-------------------------

    ALLOCATE(field%kice)
    field%kice = kice ! Number of thickness classes - always 1, as of yet

    shapeice = (/kproma, field%kice, kblks/)

    CALL add_var( field_list, prefix//'ts_icecl', field%Tsurf ,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('ts_icecl', 'C','surface temperature',datatype_flt),&
      &          grib2_var(10,2,8, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.,                             &
      &          lopenacc=.TRUE.)

    __acc_attach(field%Tsurf )
    CALL add_var( field_list, prefix//'t1_icecl', field%T1 ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('t1_icecl','C','Temperature upper layer',datatype_flt), &
      &          grib2_var(10,2,8, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.,                             &
      &          lopenacc=.TRUE.)
    __acc_attach(field%T1 )
    CALL add_var( field_list, prefix//'t2_icecl', field%T2 ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('t2_icecl','C','Temperature lower layer', datatype_flt),&
      &          grib2_var(10,2,8, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.,                             &
      &          lopenacc=.TRUE.)
    __acc_attach(field%T2 )
    CALL add_var( field_list, prefix//'sit_icecl', field%hi ,                 &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('sit_icecl', 'm', 'ice thickness', datatype_flt),   &
      &          grib2_var(10,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.,                             &
      &          lopenacc=.TRUE.)
    __acc_attach(field%hi )
    CALL add_var( field_list, prefix//'hs_icecl', field%hs ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('hs_icecl', 'm', 'snow thickness', datatype_flt),   &
      &          grib2_var(10,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
      &          ldims=shapeice, lrestart=.TRUE.,                             &
      &          lopenacc=.TRUE.)
    __acc_attach(field%hs )
    CALL add_var( field_list, prefix//'qtop_icecl', field%Qtop ,              &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('qtop_icecl', 'W/m^2', 'Energy flux available for surface melting', &
      &                   datatype_flt),                                      &
      &          grib2_var(10,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
      &          ldims=shapeice, lrestart=.FALSE.,                            &
      &          lopenacc=.TRUE.)
    __acc_attach(field%Qtop )
    CALL add_var( field_list, prefix//'qbot_icecl', field%Qbot ,              &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('qbot_icecl', 'W/m^2', 'Energy flux at ice-ocean interface', datatype_flt),&
      &          grib2_var(10,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
      &          ldims=shapeice, lrestart=.FALSE.,                            &
      &          lopenacc=.TRUE.)
    __acc_attach(field%Qbot )


    CALL add_var( field_list, prefix//'sic_icecl', field%conc ,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('sic_icecl', '', 'ice concentration in each ice class', datatype_flt),&
      &          grib2_var(10,2,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.,                             &
      &          lopenacc=.TRUE.)
    __acc_attach(field%conc )

    IF (.NOT. use_tmx) THEN
     ! &       field% albvisdir_ice (nproma,field%kice,nblks),          &
     cf_desc    = t_cf_var('albvisdir_icecl', '', 'ice albedo VIS direct', datatype_flt)
     grib2_desc = grib2_var(192,128,15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( field_list, prefix//'albvisdir_icecl', field%albvisdir_ice,  &
                    & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
                    & ldims=shapeice, lrestart=.TRUE. ,                            &
                    & lopenacc=.TRUE.)
     __acc_attach(field%albvisdir_ice)

     ! &       field% albvisdif_ice (nproma,field%kice,nblks),          &
     cf_desc    = t_cf_var('albvisdif_icecl', '', 'ice albedo VIS diffuse', datatype_flt)
     grib2_desc = grib2_var(0,19,222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( field_list, prefix//'albvisdif_icecl', field%albvisdif_ice,  &
                    & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
                    & ldims=shapeice, lrestart=.TRUE. ,                            &
                    &          lopenacc=.TRUE.)
     __acc_attach(field%albvisdif_ice)

     ! &       field% albnirdir_ice (nproma,field%kice,nblks),          &
     cf_desc    = t_cf_var('albnirdir_icecl', '', 'ice albedo NIR direct', datatype_flt)
     grib2_desc = grib2_var(192,128,17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( field_list, prefix//'albnirdir_icecl', field%albnirdir_ice,  &
                    & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
                    & ldims=shapeice, lrestart=.TRUE. ,                            &
                    & lopenacc=.TRUE.)
     __acc_attach(field%albnirdir_ice)

     ! &       field% albnirdif_ice (nproma,field%kice,nblks),          &
     cf_desc    = t_cf_var('albnirdif_icecl', '', 'ice albedo NIR diffuse', datatype_flt)
     grib2_desc = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( field_list, prefix//'albnirdif_icecl', field%albnirdif_ice,  &
                    & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
                    & ldims=shapeice, lrestart=.TRUE. ,                            &
                    & lopenacc=.TRUE.)
     __acc_attach(field%albnirdif_ice)
    END IF

    !--------------------------------------
    ! Arbitrary output fields for radiation
    !--------------------------------------

    cf_desc    = t_cf_var('rad_2d', '', &
               & 'arbitrary radiation 2d-field', datatype_flt)
    grib2_desc = grib2_var(0,6,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'rad_2d', field%rad_2d,       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%rad_2d)
    
    !-------
    ! Clouds
    !-------
    cf_desc    = t_cf_var('cl', 'm2 m-2', 'cloud area fraction', datatype_flt)
    grib2_desc = grib2_var(0,6,22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'cl', field%aclc,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                                    &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                    &
                &             lower_limit=0._wp ),                                       &
                & lopenacc=.TRUE.)
    __acc_attach(field%aclc)

    cf_desc    = t_cf_var('clt', 'm2 m-2', &
               & 'total cloud cover', datatype_flt)
    grib2_desc = grib2_var(0,6,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'clt', field%aclcov,       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%aclcov)

    ! &       field% acdnc  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration', datatype_flt)
    grib2_desc = grib2_var(0,6,28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'acdnc', field%acdnc,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                                    &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,                     &
                &             lower_limit=0._wp ),                                       &
                & lopenacc=.TRUE.)
    __acc_attach(field%acdnc)

    IF (is_variable_in_output(var_name=prefix//'hur')) THEN
       cf_desc    = t_cf_var('relative_humidity', '', 'relative humidity', datatype_flt)
       grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'hur', field%hur   ,                               &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                   & lrestart = .FALSE.,                                                    &
                   & vert_interp=create_vert_interp_metadata(                               &
                   &             vert_intp_type=vintp_types("P","Z","I"),                   &
                   &             vert_intp_method=VINTP_METHOD_LIN,                         &
                   &             l_loglin=.FALSE.,                                          &
                   &             l_extrapol=.FALSE., l_pd_limit=.TRUE.,                     &
                   &             lower_limit=0._wp ),                                       &
                   & lopenacc=.TRUE.)
       __acc_attach(field%hur   )
    END IF

    cf_desc    = t_cf_var('ufts', 'W m-2',    &
               & 'energy flux at surface from thermal exchange', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ufts', field%ufts,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%ufts)
    
    cf_desc    = t_cf_var('ufvs', 'W m-2',    &
               & 'energy flux at surface from vapor exchange', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ufvs', field%ufvs,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%ufvs)
    
    cf_desc    = t_cf_var('ufcs', 'W m-2',    &
               & 'energy flux at surface from condensate', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ufcs', field%ufcs,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%ufcs)
    
    !--------------
    ! Precipitation
    !--------------
    cf_desc    = t_cf_var('prlr', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (water)', datatype_flt)
    grib2_desc = grib2_var(0,1,77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'prlr', field%rsfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%rsfl)

    cf_desc    = t_cf_var('prls', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (snow)', datatype_flt)
    grib2_desc = grib2_var(0,1,59, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'prls', field%ssfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%ssfl)

    cf_desc    = t_cf_var('rain_gsp_rate', 'kg m-2 s-1',    &
               & 'gridscale rain rate ', datatype_flt)
    grib2_desc = grib2_var(0,1,77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'rain_gsp_rate', field%rain_gsp_rate,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%rain_gsp_rate)

    cf_desc    = t_cf_var('ice_gsp_rate', 'kg m-2 s-1',    &
               & 'gridscale ice rate ', datatype_flt)
    grib2_desc = grib2_var(0,1,77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ice_gsp_rate', field%ice_gsp_rate,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%ice_gsp_rate)

    cf_desc    = t_cf_var('snow_gsp_rate', 'kg m-2 s-1',    &
               & 'gridscale snow rate ', datatype_flt)
    grib2_desc = grib2_var(0,1,56, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'snow_gsp_rate', field%snow_gsp_rate,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%snow_gsp_rate)

    cf_desc    = t_cf_var('graupel_gsp_rate', 'kg m-2 s-1',    &
               & 'gridscale graupel rate ', datatype_flt)
    grib2_desc = grib2_var(0,1,75, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'graupel_gsp_rate', field%graupel_gsp_rate,  &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%graupel_gsp_rate)

    cf_desc    = t_cf_var('hail_gsp_rate', 'kg m-2 s-1',    &
               & 'gridscale hail rate ', datatype_flt)
    grib2_desc = grib2_var(0,1,75, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'hail_gsp_rate', field%hail_gsp_rate,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%hail_gsp_rate)

    cf_desc    = t_cf_var('pr', 'kg m-2 s-1',                    &
         &                'precipitation flux',                  &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pr', field%pr,            &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%pr)

    ! &       field% totte  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('total_turbulent_energy', 'J kg-1', 'total turbulent energy', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,19,11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'totte', field%totte,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                & lrestart = .FALSE., initval = 1.e-4_wp, ldims=shape3d,   &
                & vert_interp=create_vert_interp_metadata(                 &
                &             vert_intp_type=vintp_types("P","Z","I"),     &
                &             vert_intp_method=VINTP_METHOD_LIN,           &
                &             l_loglin=.FALSE.,                            &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,       &
                &             lower_limit=0._wp ),                         &
                & lopenacc=.TRUE.)
    __acc_attach(field%totte)

    !---------------------------
    ! WMO tropopause
    !---------------------------

    ! &       field% ptp (nproma,       nblks), &
    cf_desc    = t_cf_var('ptp', 'Pa', 'tropopause air pressure', datatype_flt)
    grib2_desc = grib2_var(0,3,201, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ptp', field%ptp,          &
         &        GRID_UNSTRUCTURED_CELL, ZA_TROPOPAUSE,         &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE., initval = 20000.0_wp,       &
         &        isteptype=TSTEP_INSTANT,                       &
         &        lopenacc=.TRUE.)
    __acc_attach(field%ptp)

    !---------------------------
    ! Variables for energy diagnostic of aes physics
    !---------------------------

    CALL add_var( field_list, prefix//'cvair', field%cvair,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('cvair', 'J/kg/K',                                     &
                &          'specific heat of air at constant volume',             &
                &          datatype_flt),                                         &
                & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape3d,                                                  &
                & lrestart = .FALSE.,                                             &
                & vert_interp=create_vert_interp_metadata(                        &
                &   vert_intp_type=vintp_types("P","Z","I"),                      &
                &   vert_intp_method=VINTP_METHOD_LIN ),                          &
                & lopenacc=.TRUE.)

    __acc_attach(field%cvair)

    IF (is_variable_in_output(var_name=prefix//'q_phy')) THEN
       CALL add_var( field_list, prefix//'q_phy', field%q_phy,                       &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                   & t_cf_var('q_phy', 'W m-2',                                      &
                   &          'layer heating by physics',                            &
                   &          datatype_flt),                                         &
                   & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                   & ldims=shape3d,                                                  &
                   & lrestart = .FALSE.,                                             &
                   & vert_interp=create_vert_interp_metadata(                        &
                   &   vert_intp_type=vintp_types("P","Z","I"),                      &
                   &   vert_intp_method=VINTP_METHOD_LIN ),                          &
                   & lopenacc=.TRUE.)
       __acc_attach(field%q_phy)
    END IF

    IF (is_variable_in_output(var_name=prefix//'q_phy_vi')) THEN
       CALL add_var( field_list, prefix//'q_phy_vi', field%q_phy_vi,                 &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                   & t_cf_var('q_phy_vi', 'W m-2',                                   &
                   &          'vert. integr. heating by physics',                    &
                   &          datatype_flt),                                         &
                   & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                   & ldims=shape2d,                                                  &
                   & lrestart = .FALSE.,                                             &
                   & lopenacc=.TRUE.)
       __acc_attach(field%q_phy_vi)
    END IF

    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       IF (is_variable_in_output(var_name=prefix//'q_rad')) THEN
          CALL add_var( field_list, prefix//'q_rlw', field%q_rad,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_rad', 'W m-2',                                      &
                      &          'layer heating by LW+SW radiation',                    &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ),                          &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_rad)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'q_rad_vi')) THEN
          CALL add_var( field_list, prefix//'q_rad_vi', field%q_rad_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_rlw_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by LW+SW radiation',            &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_rad_vi)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'q_rlw')) THEN
          CALL add_var( field_list, prefix//'q_rlw', field%q_rlw,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_rlw', 'W m-2',                                      &
                      &          'layer heating by LW radiation',                       &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ),                          &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_rlw)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'q_rlw_vi')) THEN
          CALL add_var( field_list, prefix//'q_rlw_vi', field%q_rlw_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_rlw_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by LW radiation',               &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_rlw_vi)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'q_rsw')) THEN
          CALL add_var( field_list, prefix//'q_rsw', field%q_rsw,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_rsw', 'W m-2',                                      &
                      &          'layer heating by SW radiation',                       &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ),                          &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_rsw)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'q_rsw_vi')) THEN
          CALL add_var( field_list, prefix//'q_rsw_vi', field%q_rsw_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_rsw_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by SW radiation',               &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_rsw_vi)
       END IF
       !
    END IF

    IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (aes_phy_tc(jg)%dt_vdf > dt_dyn .OR.                                          &
         & is_variable_in_output(var_name=prefix//'q_vdf')) THEN
          CALL add_var( field_list, prefix//'q_vdf', field%q_vdf,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_vdf', 'W m-2',                                      &
                      &          'layer heating by vertical diffusion',                 &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ),                          &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_vdf)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'q_vdf_vi')) THEN
          CALL add_var( field_list, prefix//'q_vdf_vi', field%q_vdf_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_vdf_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by vertical diffusion',         &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & lopenacc=.TRUE.)
          __acc_attach(field%q_vdf_vi)
       END IF
       !
    END IF

!!$       cf_desc    = t_cf_var('sh_vdiff','J m-2 s-1', '', datatype_flt)
!!$       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
!!$       CALL add_var( field_list, prefix//'sh_vdiff', field%sh_vdiff,          &
!!$                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
!!$                   & lrestart = .FALSE., ldims=shape2d,                       &
!!$                   & lopenacc=.TRUE.)
!!$       __acc_attach(field%sh_vdiff)

!!$       cf_desc    = t_cf_var('qv_vdiff','kg/m^2/s', '', datatype_flt)
!!$       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
!!$       CALL add_var( field_list, prefix//'qv_vdiff', field%qv_vdiff,          &
!!$                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
!!$                   & lrestart = .FALSE., ldims=shape2d,                       &
!!$                   & lopenacc=.TRUE.)
!!$       __acc_attach(field%qv_vdiff)

    !--------------------
    ! Turbulence
    !--------------------

     ! shape2d  = (/kproma,            kblks/)
     ! shape3d  = (/kproma, klev,      kblks/)
     !shapesfc = (/kproma, ksfc_type, kblks/)
     ! shapesfc = (/kproma, kblks, ksfc_type/)

      IF (is_variable_in_output(var_name=prefix//'ri_atm')) THEN
         cf_desc    = t_cf_var('richardson_number', ' ', 'moist Richardson number', datatype_flt)
         grib2_desc = grib2_var(0, 19, 202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'ri_atm', field%ri_atm,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape3d,                         &
                   & lopenacc=.TRUE.)
         __acc_attach(field%ri_atm)
      END IF

      IF (is_variable_in_output(var_name=prefix//'mixlen')) THEN
         cf_desc    = t_cf_var('mixing_length', 'm', 'mixing length', datatype_flt)
         grib2_desc = grib2_var(0, 19, 201, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'mixlen', field%mixlen,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., initval = -999._wp, ldims=shape3d,     &
                   & lopenacc=.TRUE.)
         __acc_attach(field%mixlen)
      END IF

      ! &       field% tottem0 (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('totte', 'm2 s-2', 'TTE at step t', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'tottem0', field%tottem0,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                  & lrestart = .FALSE., initval = 1.e-4_wp, ldims=shape3d,     &
                  & lopenacc=.TRUE.)
      __acc_attach(field%tottem0)

      ! &       field% tottem1  (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('totte', 'm2 s-2', 'TTE at step t-dt', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'tottem1', field%tottem1,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                  & lrestart = .TRUE., initval = 1.e-4_wp, ldims=shape3d,      &
                  & lopenacc=.TRUE.)
      
      __acc_attach(field%tottem1)
      
      ! &       field% cptgz  (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('cptgz', 'm2 s-2', 'dry static energy', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'cptgz', field%cptgz,                  &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                  & lrestart = .TRUE., initval = 1.e-4_wp, ldims=shape3d,      &
                  & lopenacc=.TRUE.)
      __acc_attach(field%cptgz)

      ! &       field% cptgzvi     (nproma,nblks),          &
      IF (is_variable_in_output(var_name=prefix//'cptgzvi')) THEN
          cf_desc    = t_cf_var('vertically integrated dry static energy', 'm2 s-2', 'vert_int_dry_static_energy', &
               &                datatype_flt)
          grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( field_list, prefix//'cptgzvi', field%cptgzvi,              &
               &        GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE,                       &
               &        cf_desc, grib2_desc,                                         &
               &        ldims=shape2d,                                               &
               &        lrestart = .FALSE.,                                          &
               &        isteptype=TSTEP_INSTANT,                                     &
               &        lopenacc=.TRUE.)
          __acc_attach(field%cptgzvi)
      END IF

      ! &       field% udynvi     (nproma,nblks),          &
      cf_desc    = t_cf_var('u_dyn_vi','J m-2','vertically integrated moist internal energy after dynamics', &
           &                datatype_flt)
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'udynvi', field%udynvi,                      &
           &        GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE,                       &
           &        cf_desc, grib2_desc,                                         &
           &        ldims=shape2d,                                               &
           &        lrestart = .FALSE.,                                          &
           &        isteptype=TSTEP_INSTANT,                                     &
           &        lopenacc=.TRUE.)
      __acc_attach(field%udynvi)

      ! &       field% duphyvi     (nproma,nblks),          &
      cf_desc    = t_cf_var('du_phy_vi','J m-2','change of vertically integrated moist internal energy by physics', &
           &                datatype_flt)
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'duphyvi', field%duphyvi,                      &
           &        GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE,                       &
           &        cf_desc, grib2_desc,                                         &
           &        ldims=shape2d,                                               &
           &        lrestart = .FALSE.,                                          &
           &        isteptype=TSTEP_INSTANT,                                     &
           &        lopenacc=.TRUE.)
      __acc_attach(field%duphyvi)

      ! &       field% utmxvi     (nproma,nblks),          &
      IF (use_tmx .AND. (     is_variable_in_output(var_name=prefix//'utmxvi') &
                     &   .OR. is_variable_in_output(var_name=prefix//'utmxvi_gmean'))) THEN
          cf_desc    = t_cf_var('u_tmx_vi','J m-2','vertically integrated moist internal energy after tmx', &
               &                datatype_flt)
          grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( field_list, prefix//'utmxvi', field%utmxvi,                  &
               &        GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE,                       &
               &        cf_desc, grib2_desc,                                         &
               &        ldims=shape2d,                                               &
               &        lrestart = .FALSE.,                                          &
               &        isteptype=TSTEP_INSTANT,                                     &
               &        lopenacc=.TRUE.)
          __acc_attach(field%utmxvi)
      END IF

      ! REMARK: required for art emmision handling
      !IF (is_variable_in_output(var_name=prefix//'cfm')) THEN
      cf_desc    = t_cf_var('turb_exchng_coeff_momentum', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'cfm', field%cfm,                      &
           &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
           &        lrestart = .FALSE., ldims=shape3d,                         &
           &        lopenacc=.TRUE.)

      __acc_attach(field%cfm)
      !END IF
      !
      contvar_is_in_output = .FALSE.
      DO jsfc = 1,ksfc_type
         var_name=prefix//'cfm_'//csfc(jsfc)
         IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
            contvar_is_in_output = .TRUE.
         END IF
      END DO
      !
      IF (contvar_is_in_output .OR. use_tmx) THEN
         CALL add_var( field_list, prefix//'cfm_tile', field%cfm_tile,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                     & t_cf_var('turb_exchng_coeff_momentum', '', '', datatype_flt),&
                     & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                     & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                     & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,        &
                     & lopenacc=.TRUE.)
         __acc_attach(field%cfm_tile)
         ALLOCATE(field%cfm_tile_ptr(ksfc_type))
      END IF
      !
      DO jsfc = 1,ksfc_type
         var_name=prefix//'cfm_'//csfc(jsfc)
         IF (is_variable_in_output(var_name=TRIM(var_name)) .OR. use_tmx) THEN
            CALL add_ref( field_list, prefix//'cfm_tile',                              &
                        & TRIM(var_name), field%cfm_tile_ptr(jsfc)%p,                  &
                        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                        & t_cf_var('turb_exchng_coeff_momentum_'//csfc(jsfc), '', '', datatype_flt),&
                        & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                        & ref_idx=jsfc, ldims=shape2d, lrestart=use_tmx,               &
                        & lmiss=.TRUE., missval=cdimissval )
         END IF
      END DO


      ! REMARK: required for art sedimentation handling
      !IF (is_variable_in_output(var_name=prefix//'cfh')) THEN
      cf_desc    = t_cf_var('turb_exchng_coeff_heat', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'cfh', field%cfh,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                  & lrestart = .FALSE., ldims=shape3d,                         &
                  & lopenacc=.TRUE.)
      __acc_attach(field%cfh)
      !END IF
      !
      contvar_is_in_output = .FALSE.
      DO jsfc = 1,ksfc_type
         var_name=prefix//'cfh_'//csfc(jsfc)
         IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
            contvar_is_in_output = .TRUE.
         END IF
      END DO
      !
      IF (contvar_is_in_output .OR. use_tmx) THEN
         CALL add_var( field_list, prefix//'cfh_tile', field%cfh_tile,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                     & t_cf_var('turb_exchng_coeff_heat', '', '', datatype_flt),    &
                     & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                     & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                     & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,        &
                     & lopenacc=.TRUE.)
         __acc_attach(field%cfh_tile)
         ALLOCATE(field%cfh_tile_ptr(ksfc_type))
      END IF
      !
      DO jsfc = 1,ksfc_type
         var_name=prefix//'cfh_'//csfc(jsfc)
         IF (is_variable_in_output(var_name=TRIM(var_name)) .OR. use_tmx) THEN
            CALL add_ref( field_list, prefix//'cfh_tile',                              &
                        & TRIM(var_name), field%cfh_tile_ptr(jsfc)%p,                  &
                        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                        & t_cf_var('turb_exchng_coeff_heat_'//csfc(jsfc), '', '',      &
                        &          datatype_flt),                                      &
                        & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                        & ref_idx=jsfc, ldims=shape2d, lrestart = .FALSE.,             &
                        & lmiss=.TRUE., missval=cdimissval )
         END IF
      END DO


      IF (is_variable_in_output(var_name=prefix//'cfv')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_water_var', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cfv', field%cfv,                      &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d,                         &
                     & lopenacc=.TRUE.)
         __acc_attach(field%cfv)
      END IF

      IF (is_variable_in_output(var_name=prefix//'cfv')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_totte', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cftotte', field%cftotte,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d,                         &
                     & lopenacc=.TRUE.)
         __acc_attach(field%cftotte)
      END IF

      IF (is_variable_in_output(var_name=prefix//'cfthv')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_thv', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cfthv', field%cfthv,                  &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d,                         &
                     & lopenacc=.TRUE.)
         __acc_attach(field%cfthv)
      END IF

      cf_desc    = t_cf_var('Coriolis_param', 's-1', 'Coriolis parameter', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'coriol', field%coriol,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                  & lrestart = .FALSE., ldims=shape2d,                       &
                  & lopenacc=.TRUE.)
      __acc_attach(field%coriol)

      IF (is_variable_in_output(var_name=prefix//'hdtcbl')) THEN
         cf_desc    = t_cf_var('height_pbl_top', 'm', 'height of PBL top', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'hdtcbl', field%hdtcbl,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape2d,                       &
                     & lopenacc=.TRUE.)
         __acc_attach(field%hdtcbl)
      END IF

      !-----------------------------------
      ! &       field% z0m(nproma,nblks), &
      IF (is_variable_in_output(var_name=prefix//'z0m') .OR. .NOT. use_tmx) THEN
        cf_desc    = t_cf_var('z0m', '', 'aerodynamic roughness length (mom), grid box mean', &
             &                datatype_flt)
        grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( field_list, prefix//'z0m', field%z0m,                    &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                  & lrestart = .FALSE., ldims=shape2d,                       &
                  & lopenacc=.TRUE.)
        __acc_attach(field%z0m)
      END IF

      ! &       field% z0m_tile(nproma,nblks,nsfc_type), &
      contvar_is_in_output = .FALSE.
      DO jsfc = 1,ksfc_type
         var_name=prefix//'z0m_'//csfc(jsfc)
         IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
            contvar_is_in_output = .TRUE.
         END IF
      END DO
      !
      IF (contvar_is_in_output .OR. .NOT. use_tmx) THEN
        CALL add_var( field_list, prefix//'z0m_tile', field%z0m_tile,                         &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                     &
                    & t_cf_var('z0m_tile', '', 'aerodynamic roughness length (mom)', datatype_flt), &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),           &
                    & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                    & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                   &
                    & lopenacc=.TRUE.)
        __acc_attach(field%z0m_tile)
      END IF

      ALLOCATE(field%z0m_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        var_name=prefix//'z0m_'//csfc(jsfc)
        IF (is_variable_in_output(var_name=TRIM(var_name)) .OR. .NOT. use_tmx) THEN
          CALL add_ref( field_list, prefix//'z0m_tile',                              &
                      & TRIM(var_name), field%z0m_tile_ptr(jsfc)%p,                  &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                      & t_cf_var('z0m_'//csfc(jsfc), '','', datatype_flt),           &
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ref_idx=jsfc, ldims=shape2d, lrestart = .NOT. use_tmx,       &
                      & lmiss=.TRUE., missval=cdimissval )
        END IF
      END DO

      ! &       field% z0h(nproma,nblks), &
      IF (use_tmx) THEN
        IF (is_variable_in_output(var_name=prefix//'z0h')) THEN
          cf_desc    = t_cf_var('z0h', '', 'aerodynamic roughness length (heat), grid box mean', &
              &                datatype_flt)
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( field_list, prefix//'z0h', field%z0h,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                      & lrestart = .FALSE., ldims=shape2d,                       &
                      & lopenacc=.TRUE.)
          __acc_attach(field%z0h)
        END IF

        ! &       field% z0h_tile(nproma,nblks,nsfc_type), &
        contvar_is_in_output = .FALSE.
        DO jsfc = 1,ksfc_type
           var_name=prefix//'z0h_'//csfc(jsfc)
           IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
              contvar_is_in_output = .TRUE.
           END IF
        END DO
        !
        IF (contvar_is_in_output) THEN
          CALL add_var( field_list, prefix//'z0h_tile', field%z0h_tile,                         &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                     &
                    & t_cf_var('z0h_tile', '', 'aerodynamic roughness length (heat)', datatype_flt), &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),           &
                    & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                    & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                   &
                    & lopenacc=.TRUE.)
          __acc_attach(field%z0h_tile)
        END IF

        ALLOCATE(field%z0h_tile_ptr(ksfc_type))
        DO jsfc = 1,ksfc_type
          var_name=prefix//'z0h_'//csfc(jsfc)
          IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
            CALL add_ref( field_list, prefix//'z0h_tile',                            &
                      & TRIM(var_name), field%z0h_tile_ptr(jsfc)%p,                  &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                      & t_cf_var('z0h_'//csfc(jsfc), '','', datatype_flt),           &
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ref_idx=jsfc, ldims=shape2d, lrestart = .FALSE.,             &
                      & lmiss=.TRUE., missval=cdimissval )
          END IF
        END DO
      ELSE
        ! &        field% z0h_lnd(nproma, nblks), &
         cf_desc    = t_cf_var('z0h_lnd', '', 'roughness length heat, land', &
           &                datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'z0h_lnd', field%z0h_lnd,            &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                     & lrestart = .TRUE., ldims=shape2d,                        &
                     & lmiss=.TRUE., missval=cdimissval,                        &
                     & lopenacc=.TRUE.)
         __acc_attach(field%z0h_lnd)
      END IF

      !-----------------------------------

      ! &       field% ustar  (nproma,nblks),                &
      cf_desc    = t_cf_var('friction_velocity', 'm s-1', 'friction velocity', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'ustar', field%ustar,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                  & lrestart = .TRUE., initval = 1._wp, ldims=shape2d,       &
                  & lopenacc=.TRUE.)
      __acc_attach(field%ustar)

      IF (is_variable_in_output(var_name=prefix//'wstar')) THEN
         cf_desc    = t_cf_var('conv_velocity_scale', 'm s-1', 'convective velocity scale', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'wstar', field%wstar,                &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape2d,                       &
                     & lopenacc=.TRUE.)
         __acc_attach(field%wstar)
      END IF

      ! &       field% wstar_tile(nproma,nblks,nsfc_type), &
      CALL add_var( field_list, prefix//'wstar_tile', field%wstar_tile,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('wstar_tile', '', 'convective velocity scale', datatype_flt),&
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                  & ldims=shapesfc, lcontainer=.TRUE.,                           &
                  & lrestart=.FALSE., loutput=.FALSE.,                           &
                  & lopenacc=.TRUE.)
      __acc_attach(field%wstar_tile)

      ALLOCATE(field%wstar_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        CALL add_ref( field_list, prefix//'wstar_tile',                            &
                    & prefix//'wstar_'//csfc(jsfc), field%wstar_tile_ptr(jsfc)%p,  &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('wstar_'//csfc(jsfc), '','', datatype_flt),         &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                    & ref_idx=jsfc, ldims=shape2d, lrestart=.TRUE.                 )
      END DO


      IF (     is_variable_in_output(var_name=prefix//'kedisp')                     &
          .OR. is_variable_in_output(var_name=prefix//'kedisp_gmean')               &
          .OR. is_variable_in_output(var_name=prefix//'uphybal_gmean')) THEN
         cf_desc    = t_cf_var('vert_int_dissip_kin_energy', 'W/m2',                &
                     &         'vert. integr. dissip. kin. energy', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'kedisp', field%kedisp,                  &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,     &
                     & lrestart=.FALSE., ldims=shape2d,                             &
                     & lopenacc=.TRUE.)
         __acc_attach(field%kedisp)
      END IF

      ! &       field% ocu    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_u', 'm/s', 'u-component of ocean current/ice', datatype_flt)
      grib2_desc = grib2_var(10,1,2, iextbits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'ocu', field%ocu, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,   &
        &           cf_desc, grib2_desc, ldims=shape2d,   &
        &           lrestart=.TRUE., initval=0._wp,       &
        &           lopenacc=.TRUE.)
      __acc_attach(field%ocu)

      ! &       field% ocv    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_v', 'm/s', 'v-component of ocean current/ice', datatype_flt)
      grib2_desc = grib2_var(10,1,3, iextbits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'ocv', field%ocv, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,   &
        &           cf_desc, grib2_desc, ldims=shape2d,   &
        &           lrestart=.TRUE., initval=0._wp,       &
        &           lopenacc=.TRUE.)
      __acc_attach(field%ocv)

    !-----------------------
    ! Surface
    !-----------------------

    cf_desc    = t_cf_var('surface_altitude', 'm',   &
                &         'surface altitude', datatype_flt)
    grib2_desc = grib2_var(0,3,6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_ref( ext_atm_list, 'topography_c',                                &
                & prefix//'orog', field%orog,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & cf_desc, grib2_desc,                                         &
                & ref_idx=1, ldims=shape2d,                                    &
                & loutput=.TRUE.,                                              &
                & isteptype=TSTEP_CONSTANT                                     )
    __acc_attach(field%orog)

    cf_desc    = t_cf_var('land_area_fraction', 'm2/m2',   &
                &         'cell area fraction occupied by land including lakes', datatype_flt)
    grib2_desc = grib2_var(2,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sftlf', field%sftlf,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_CONSTANT,                                &
                & lopenacc=.TRUE.)
    __acc_attach(field%sftlf)

    cf_desc    = t_cf_var('land_ice_area_fraction', 'm2/m2',   &
                &         'cell area fraction occupied by land ice', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sftgif', field%sftgif,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_CONSTANT,                                &
                & lopenacc=.TRUE.)
    __acc_attach(field%sftgif)

    cf_desc    = t_cf_var('ocean_area_fraction', 'm2/m2',   &
                &         'cell area fraction occupied by ocean', datatype_flt)
    grib2_desc = grib2_var(2,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sftof', field%sftof,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_CONSTANT,                                &
                & lopenacc=.TRUE.)
    __acc_attach(field%sftof)


    ! &       field% lsmask (nproma, nblks),                 &
    cf_desc    = t_cf_var('land_cover', '', 'land cover', datatype_flt)
    grib2_desc = grib2_var(1, 2, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'land', field%lsmask,              &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.FALSE., ldims=shape2d,                         &
              & lopenacc=.TRUE.)
    __acc_attach(field%lsmask)

    ! &       field% glac   (nproma, nblks),                 &
    cf_desc    = t_cf_var('glacier_cover', '', 'fraction of land covered by glaciers', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'glac', field%glac,                &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.FALSE., ldims=shape2d,                         &
              & lopenacc=.TRUE.)
    __acc_attach(field%glac)

    ! &       field% seaice (nproma, nblks),                 &
    cf_desc    = t_cf_var('sea_ice_cover', '', 'fraction of ocean covered by sea ice', &
         &                datatype_flt)
    grib2_desc = grib2_var(10,2,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sic', field%seaice,    &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,         &
      &           cf_desc, grib2_desc, ldims=shape2d,         &
      &           lrestart=.TRUE.,                            &
      & lopenacc=.TRUE.)
    __acc_attach(field%seaice)

    ! &       field% alake (nproma, nblks),                 &
    cf_desc    = t_cf_var('alake', '', 'fraction of lakes', &
         &                datatype_flt)
    grib2_desc = grib2_var(1,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'alake', field%alake,              &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.FALSE., ldims=shape2d,                         &
              & lopenacc=.TRUE.)
    __acc_attach(field%alake)

    ! &       field% lake_ice_frc (nproma, nblks),                 &
    cf_desc    = t_cf_var('lake_ice_frc', '', 'fraction of ice on lakes', & 
         &                datatype_flt)
    grib2_desc = grib2_var(1,2,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'lake_ice_frc', field%lake_ice_frc,  &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
              & initval=0._wp, lrestart=.TRUE., ldims=shape2d,             &
              & lopenacc=.TRUE.)
    __acc_attach(field%lake_ice_frc)

    !-----------------------------------
    ! &       field% ts(nproma,nblks), &
    cf_desc    = t_cf_var('surface_temperature', 'K', 'surface temperature', datatype_flt)
    grib2_desc = grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ts', field%ts,                    &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.TRUE., ldims=shape2d,                          &
              & lopenacc=.TRUE.)
    __acc_attach(field%ts)

    ! &       field% ts_tile(nproma,nblks,nsfc_type), &
    CALL add_var( field_list, prefix//'ts_tile', field%ts_tile,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('ts_tile', 'K', 'surface temperature on tiles', datatype_flt), &
                & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,        &
                & lopenacc=.TRUE.)
    __acc_attach(field%ts_tile)

    ALLOCATE(field%ts_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      CALL add_ref( field_list, prefix//'ts_tile',                                   &
                  & prefix//'ts_'//csfc(jsfc), field%ts_tile_ptr(jsfc)%p,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                  & t_cf_var('ts_'//csfc(jsfc), 'K',                                 &
                  &          'surface temperature on '//csfc(jsfc), datatype_flt),   &
                  & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),           &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.TRUE.,                    &
                  & lmiss=.TRUE., missval=cdimissval )
    END DO
    !-----------------------------------

    contvar_is_in_output = .FALSE.
    DO jsfc = 1,ksfc_type
       var_name=prefix//'qs_sfc_'//csfc(jsfc)
       IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
          contvar_is_in_output = .TRUE.
       END IF
    END DO
    !
    IF (contvar_is_in_output) THEN
       CALL add_var( field_list, prefix//'qs_sfc_tile', field%qs_sfc_tile,        &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                   & t_cf_var('qs_sfc_tile', '', '', datatype_flt),               &
                   & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                   & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                   & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,        &
                   & lopenacc=.TRUE.)
       __acc_attach(field%qs_sfc_tile)
       ALLOCATE(field%qs_sfc_tile_ptr(ksfc_type))
    END IF
    !
    DO jsfc = 1,ksfc_type
       var_name=prefix//'qs_sfc_'//csfc(jsfc)
       IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
          CALL add_ref( field_list, prefix//'qs_sfc_tile',                           &
                      & TRIM(var_name), field%qs_sfc_tile_ptr(jsfc)%p,               &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                      & t_cf_var('qs_sfc_'//csfc(jsfc), '', '', datatype_flt),       &
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,               &
                      & lmiss=.TRUE., missval=cdimissval )
       END IF
    END DO

    !-----------------------------------
    ! &       field% albedo (nproma,nblks),          &
    cf_desc    = t_cf_var('albedo', '', 'surface albedo', datatype_flt)
    grib2_desc = grib2_var(0,19,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albedo', field%albedo,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d,                         &
                & lopenacc=.TRUE.)
    __acc_attach(field%albedo)

    ! &       field% albvisdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdir', '', 'albedo VIS direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir', field%albvisdir,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d ,                         &
                & lopenacc=.TRUE.)
    __acc_attach(field%albvisdir)

    ! &       field% albvisdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdif', '', 'albedo VIS diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif', field%albvisdif,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d ,                         &
                & lopenacc=.TRUE.)
    __acc_attach(field%albvisdif)

    ! &       field% albnirdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdir', '', 'albedo NIR direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir', field%albnirdir,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d ,                         &
                & lopenacc=.TRUE.)
    __acc_attach(field%albnirdir)

    ! &       field% albnirdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdif', '', 'albedo NIR diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif', field%albnirdif,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d ,                         &
                & lopenacc=.TRUE.)
    __acc_attach(field%albnirdif)

    ! &       field% albvisdir_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albvisdir_tile', '', 'albedo VIS direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir_tile', field%albvisdir_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.,                                    &
                & lopenacc=.TRUE.)
    __acc_attach(field%albvisdir_tile)

    ! &       field% albvisdif_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albvisdif_tile', '', 'albedo VIS diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif_tile', field%albvisdif_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.,                                    &
                & lopenacc=.TRUE.)
    __acc_attach(field%albvisdif_tile)

    ! &       field% albnirdir_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albnirdir_tile', '', 'albedo NIR direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir_tile', field%albnirdir_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.,                                    &
                & lopenacc=.TRUE.)
    __acc_attach(field%albnirdir_tile)

    ! &       field% albnirdif_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albnirdif_tile', '', 'albedo NIR diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif_tile', field%albnirdif_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.,                                    &
                & lopenacc=.TRUE.)
    __acc_attach(field%albnirdif_tile)

    ! &       field% albedo_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albedo_tile', '', 'albedo', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albedo_tile', field%albedo_tile,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.,                                    &
                & lopenacc=.TRUE.)
    __acc_attach(field%albedo_tile)

    ALLOCATE(field%albvisdir_tile_ptr(ksfc_type), field%albvisdif_tile_ptr(ksfc_type), &
             field%albnirdir_tile_ptr(ksfc_type), field%albnirdif_tile_ptr(ksfc_type), &
             field%albedo_tile_ptr(ksfc_type)                                          )

    DO jsfc = 1,ksfc_type
      CALL add_ref( field_list, prefix//'albvisdir_tile',                              &
                  & prefix//'albvisdir_'//csfc(jsfc), field%albvisdir_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albvisdir_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=use_tmx,                     &
                  & lmiss=.TRUE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albvisdif_tile',                              &
                  & prefix//'albvisdif_'//csfc(jsfc), field%albvisdif_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albvisdif_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=use_tmx,                     &
                  & lmiss=.TRUE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albnirdir_tile',                              &
                  & prefix//'albnirdir_'//csfc(jsfc), field%albnirdir_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albnirdir_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=use_tmx,                     &
                  & lmiss=.TRUE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albnirdif_tile',                              &
                  & prefix//'albnirdif_'//csfc(jsfc), field%albnirdif_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albnirdif_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=use_tmx,                     &
                  & lmiss=.FALSE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albedo_tile',                                 &
                  & prefix//'albedo_'//csfc(jsfc), field%albedo_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albedo_'//csfc(jsfc), '', '', datatype_flt),             &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                     &
                  & lmiss=.TRUE., missval=cdimissval )
    END DO

    ! &       field% emissivity (nproma,nblks),          &
    cf_desc    = t_cf_var('emissivity', '', 'longwave surface emissivity', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'emissivity', field%emissivity,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d,                         &
                & lopenacc=.TRUE.)
    __acc_attach(field%emissivity)

    !---------------------------
    ! Surface fluxes
    !---------------------------
    ! gridbox mean

    CALL add_var( field_list, prefix//'evspsbl', field%evap,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evap', 'kg m-2 s-1', 'evaporation',           &
                & datatype_flt),                                          &
                & grib2_var(0,1,6,iextbits, GRID_UNSTRUCTURED, GRID_CELL),&
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT,                                &
                & lopenacc=.TRUE.)


    __acc_attach(field%evap)

    CALL add_var( field_list, prefix//'hfls', field%lhflx,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('lhflx', 'W m-2 ', 'latent heat flux',         &
                & datatype_flt),                                          &
                & grib2_var(0,0,10, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT,                                &
                & lopenacc=.TRUE.)

    __acc_attach(field%lhflx)

    CALL add_var( field_list, prefix//'hfss', field%shflx,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('shflx', 'W m-2 ', 'sensible heat flux',       &
                & datatype_flt),                                          &
                & grib2_var(0,0,11, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT,                                &
                & lopenacc=.TRUE.)

    __acc_attach(field%shflx)

    !---------------------------------
    ! values on tiles

    CALL add_var( field_list, prefix//'rsns_tile',field%swflxsfc_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('rsns_tile', 'W m-2',                          &
                &          'shortwave net flux at surface on tiles',      &
                &          datatype_flt),                                 &
                & grib2_var(0,4,9, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,   &
                & lopenacc=.TRUE.)

    __acc_attach(field%swflxsfc_tile)

    CALL add_var( field_list, prefix//'rlns_tile',field%lwflxsfc_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('rlns_tile', 'W m-2',                          &
                &          'longwave net flux at surface on tiles',       &
                &          datatype_flt),                                 &
                & grib2_var(0,5,5, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,   &
                & lopenacc=.TRUE.)

    __acc_attach(field%lwflxsfc_tile)

    CALL add_var( field_list, prefix//'evspsbl_tile', field%evap_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evspsbl_tile', 'kg m-2 s-1',                  &
                &          'evaporation on tiles', datatype_flt),         &
                & grib2_var(0,1,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,   &
                & lopenacc=.TRUE.)

    __acc_attach(field%evap_tile)

    CALL add_var( field_list, prefix//'hfls_tile', field%lhflx_tile,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfls_tile', 'W m-2',                          &
                &          'latent heat flux on tiles', datatype_flt),    &
                & grib2_var(0,0,10, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,   &
                & lopenacc=.TRUE.)

    __acc_attach(field%lhflx_tile)

    CALL add_var( field_list, prefix//'hfss_tile', field%shflx_tile,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfss_tile', 'W m-2',                          &
                &          'sensible heat flux on tiles', datatype_flt),  &
                & grib2_var(0,0,11, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,   &
                & lopenacc=.TRUE.)

    __acc_attach(field%shflx_tile)

    CALL add_var( field_list, prefix//'frac_tile', field%frac_tile,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('frac_tile', '%',                              &
                &          'surface fraction of tiles', datatype_flt),    &
                & grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,   &
                & lopenacc=.TRUE.)

    __acc_attach(field%frac_tile)

    ALLOCATE(field%swflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%lwflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%evap_tile_ptr(ksfc_type))
    ALLOCATE(field%lhflx_tile_ptr(ksfc_type))
    ALLOCATE(field%shflx_tile_ptr(ksfc_type))
    ALLOCATE(field%frac_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'rsns_tile',                                   &
                  & prefix//'rsns_'//csfc(jsfc), field%swflxsfc_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('rsns_'//csfc(jsfc), 'W m-2',                             &
                  &          'shortwave net flux at surface on tile '//csfc(jsfc),     &
                  &          datatype_flt),                                            &
                  & grib2_var(0,4,9, ibits, GRID_UNSTRUCTURED, GRID_CELL),             &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                     &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'rlns_tile',                                   &
                  & prefix//'rlns_'//csfc(jsfc), field%lwflxsfc_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('rlns_'//csfc(jsfc), 'W m-2',                             &
                  &          'longwave net flux at surface on tile '//csfc(jsfc),      &
                  &          datatype_flt),                                            &
                  & grib2_var(0,5,5, ibits, GRID_UNSTRUCTURED, GRID_CELL),             &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                     &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'evspsbl_tile',                                &
                  & prefix//'evspsbl_'//csfc(jsfc), field%evap_tile_ptr(jsfc)%p,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('evspsbl_'//csfc(jsfc), 'kg m-2 s-1',                     &
                  &          'evaporation on tile '//csfc(jsfc), datatype_flt),        &
                  & grib2_var(0,1,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),             &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                     &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'hfls_tile',                                   &
                  & prefix//'hfls_'//csfc(jsfc), field%lhflx_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('hfls_'//csfc(jsfc), 'W m-2',                             &
                  &          'latent heat flux on tile '//csfc(jsfc), datatype_flt),   &
                  & grib2_var(0,0,10, ibits, GRID_UNSTRUCTURED, GRID_CELL),            &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                     &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'hfss_tile',                                   &
                  & prefix//'hfss_'//csfc(jsfc), field%shflx_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('hfss_'//csfc(jsfc), 'W m-2',                             &
                  &          'sensible heat flux on tile '//csfc(jsfc),datatype_flt),  &
                  & grib2_var(0,0,11, ibits, GRID_UNSTRUCTURED, GRID_CELL),            &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                     &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'frac_tile',                                   &
                  & prefix//'frac_'//csfc(jsfc), field%frac_tile_ptr(jsfc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('frac_'//csfc(jsfc), '%',                                 &
                  &          'surface fraction of tile '//csfc(jsfc),                  &
                  &          datatype_flt),                                            &
                  & grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                  & ref_idx=jsfc, ldims=shape2d, lmiss=.TRUE., missval=cdimissval )

    END DO

    !-----------------------------------------
    ! wind stress, grid box mean
    !-----------------------------------------

    CALL add_var( field_list, prefix//'tauu', field%u_stress        ,           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress','N m-2','u-momentum flux at the surface', &
                &          datatype_flt),                                       &
                & grib2_var(0,2,17, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%u_stress        )

    CALL add_var( field_list, prefix//'tauv', field%v_stress,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress','N m-2','v-momentum flux at the surface', &
                &          datatype_flt),                                       &
                & grib2_var(0,2,18, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%v_stress)

    ! wind stress, instantaneous tile values 

    CALL add_var( field_list, prefix//'tauu_tile', field%u_stress_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress_tile', 'N m-2',                            &
                &          'u-momentum flux at the surface on tiles',           &
                &          datatype_flt),                                       &
                & grib2_var(0,2,17, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,             &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
                & lopenacc=.TRUE.)

    __acc_attach(field%u_stress_tile)

    CALL add_var( field_list, prefix//'tauv_tile', field%v_stress_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress_tile', 'N m-2',                            &
                &          'v-momentum flux at the surface on tiles',           &
                &          datatype_flt),                                       &
                & grib2_var(0,2,18, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,             &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
                & lopenacc=.TRUE.)

    __acc_attach(field%v_stress_tile)

    ALLOCATE(field%u_stress_tile_ptr(ksfc_type))
    ALLOCATE(field%v_stress_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'tauu_tile',                                &
                  & prefix//'tauu_'//csfc(jsfc), field%u_stress_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('u_stress_'//csfc(jsfc), 'N m-2',                      &
                  &          'u-momentum flux at the surface on tile '//csfc(jsfc), &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,17, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                  &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'tauv_tile',                                &
                  & prefix//'tauv_'//csfc(jsfc), field%v_stress_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('v_stress_'//csfc(jsfc), 'N m-2',                      &
                  &          'v-momentum flux at the surface on tile '//csfc(jsfc), &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,18, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                  &
                  & lmiss=.TRUE., missval=cdimissval )
    END DO

    !-----------------------------------------
    ! near surface diagnostics, grid box mean
    !-----------------------------------------

    CALL add_var( field_list, prefix//'sfcwind', field%sfcwind,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                        &
                & t_cf_var('sfcwind','m s-1','10m windspeed',                   &
                &          datatype_flt),                                       &
                & grib2_var(0,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%sfcwind)

    CALL add_var( field_list, prefix//'uas', field%uas,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                        &
                & t_cf_var('uas','m s-1','zonal wind in 10m',                   &
                &          datatype_flt),                                       &
                & grib2_var(0,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%uas)

    CALL add_var( field_list, prefix//'vas', field%vas,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                        &
                & t_cf_var('vas','m s-1','meridional wind in 10m',              &
                &          datatype_flt),                                       &
                & grib2_var(0,2,3, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%vas)

    CALL add_var( field_list, prefix//'tas', field%tas,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                         &
                & t_cf_var('tas','K','temperature in 2m',                       &
                &          datatype_flt),                                       &
                & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%tas)

    CALL add_var( field_list, prefix//'qv2m', field%qv2m,                      &
                   & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                     &
                   & t_cf_var('qv2m','kg kg-1','specific humidity in 2m',      &
                   &          datatype_flt),                                   &
                   & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                   & ldims=shape2d,                                            &
                   & lrestart = .FALSE.,                                       &
                   & isteptype=TSTEP_INSTANT,                                  &
                   & lopenacc=.TRUE.)
    __acc_attach(field%qv2m)

    CALL add_var( field_list, prefix//'dew2', field%dew2,                      &
                   & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                     &
                   & t_cf_var('dew2','K','dew point temperature in 2m',        &
                   &          datatype_flt),                                   &
                   & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                   & ldims=shape2d,                                            &
                   & lrestart = .FALSE.,                                       &
                   & isteptype=TSTEP_INSTANT,                                  &
                   & lopenacc=.TRUE.)
    __acc_attach(field%dew2)

    CALL add_var( field_list, prefix//'tasmax', field%tasmax,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                         &
                & t_cf_var('tasmax','K','maximum 2m temperature',               &
                &          datatype_flt),                                       &
                & grib2_var(0,0,4, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & initval = -99._wp, resetval = -99._wp,                        &
                & isteptype=TSTEP_MAX,                                          &
                & action_list=actions(new_action(ACTION_RESET,"P1D")),          &
                & lopenacc=.TRUE.)

    __acc_attach(field%tasmax)

    CALL add_var( field_list, prefix//'tasmin', field%tasmin,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                         &
                & t_cf_var('tasmin','K','minimum 2m temperature',               &
                &          datatype_flt),                                       &
                & grib2_var(0,0,5, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & initval = 999._wp, resetval = 999._wp,                        &
                & isteptype=TSTEP_MIN,                                          &
                & action_list=actions(new_action(ACTION_RESET,"P1D")),          &
                & lopenacc=.TRUE.)

    __acc_attach(field%tasmin)

    !--------------------------------------
    ! near surface diagnostics, tile values
    !--------------------------------------

    CALL add_var( field_list, prefix//'sfcwind_tile', field%sfcwind_tile,       &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                        &
                & t_cf_var('sfcwind_tile','m s-1','10m windspeed on tiles',     &
                &          datatype_flt),                                       &
                & grib2_var(0,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%sfcwind_tile)

    CALL add_var( field_list, prefix//'uas_tile', field%uas_tile,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                        &
                & t_cf_var('uas_tile','m s-1','zonal wind in 10m on tiles',     &
                &          datatype_flt),                                       &
                & grib2_var(0,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%uas_tile)

    CALL add_var( field_list, prefix//'vas_tile', field%vas_tile,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                        &
                & t_cf_var('vas_tile','m s-1','meridional wind in 10m on tiles',&
                &          datatype_flt),                                       &
                & grib2_var(0,2,3, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%vas_tile)

    CALL add_var( field_list, prefix//'tas_tile', field%tas_tile,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                         &
                & t_cf_var('tas_tile','K','temperature in 2m on tiles',         &
                &          datatype_flt),                                       &
                & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT,                                      &
                & lopenacc=.TRUE.)

    __acc_attach(field%tas_tile)

    CALL add_var( field_list, prefix//'qv2m_tile', field%qv2m_tile,            &
                   & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                     &
                   & t_cf_var('qv2m_tile','kg kg-1','specific humidity in 2m on tiles',&
                   &          datatype_flt),                                   &
                   & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                   & ldims=shapesfc,                                           &
                   & lcontainer=.TRUE., lrestart=.FALSE.,                      &
                   & isteptype=TSTEP_INSTANT,                                  &
                   & lopenacc=.TRUE.)
    __acc_attach(field%qv2m_tile)

    CALL add_var( field_list, prefix//'dew2_tile', field%dew2_tile,            &
                   & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M,                     &
                   & t_cf_var('dew2_tile','K','dew point temperature in 2m on tiles',&
                   &          datatype_flt),                                   &
                   & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                   & ldims=shapesfc,                                           &
                   & lcontainer=.TRUE., lrestart=.FALSE.,                      &
                   & isteptype=TSTEP_INSTANT,                                  &
                   & lopenacc=.TRUE.)
    __acc_attach(field%dew2_tile)

    ALLOCATE(field%sfcwind_tile_ptr(ksfc_type))
    ALLOCATE(field%uas_tile_ptr(ksfc_type))
    ALLOCATE(field%vas_tile_ptr(ksfc_type))
    ALLOCATE(field%tas_tile_ptr(ksfc_type))
    ALLOCATE(field%qv2m_tile_ptr(ksfc_type))
    ALLOCATE(field%dew2_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'sfcwind_tile',                             &
                  & prefix//'sfcwind_'//csfc(jsfc), field%sfcwind_tile_ptr(jsfc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('sfcwind_'//csfc(jsfc), 'm s-1',                       &
                  &          '10m windspeed on tile '//csfc(jsfc),                  &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=use_tmx,                  &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'uas_tile',                                 &
                  & prefix//'uas_'//csfc(jsfc), field%uas_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('uas_'//csfc(jsfc), 'm s-1',                           &
                  &          'zonal wind in 10m on tile '//csfc(jsfc),              &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                  &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'vas_tile',                                 &
                  & prefix//'vas_'//csfc(jsfc), field%vas_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('vas_'//csfc(jsfc), 'm s-1',                           &
                  &          'meridional wind in 10m on tile '//csfc(jsfc),         &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,3, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                  &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'tas_tile',                                 &
                  & prefix//'tas_'//csfc(jsfc), field%tas_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('tas_'//csfc(jsfc), 'K',                               &
                  &          'temperature in 2m on tile '//csfc(jsfc),              &
                  &          datatype_flt),                                         &
                  & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                  &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'qv2m_tile',                              &
                & prefix//'qv2m_'//csfc(jsfc), field%qv2m_tile_ptr(jsfc)%p,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                & t_cf_var('qv2m_'//csfc(jsfc), 'kg kg-1',                        &
                &          'specific humidity in 2m on tile '//csfc(jsfc),        &
                &          datatype_flt),                                         &
                & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                  &
                & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'dew2_tile',                              &
                & prefix//'dew2_'//csfc(jsfc), field%dew2_tile_ptr(jsfc)%p,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                & t_cf_var('dew2_'//csfc(jsfc), 'K',                              &
                &          'dew point temperature in 2m on tile '//csfc(jsfc),    &
                &          datatype_flt),                                         &
                & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                & ref_idx=jsfc, ldims=shape2d, lrestart=.FALSE.,                  &
                & lmiss=.TRUE., missval=cdimissval )

    END DO

    ! global diagnostics
    cf_desc    = t_cf_var('tas_gmean', 'K', 'global mean temperature at 2m', datatype_flt,'tas_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'tas_gmean', field%tas_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%tas_gmean)

    cf_desc    = t_cf_var('rsdt_gmean', 'W m-2', 'global mean toa incident shortwave radiation', datatype_flt,'rsdt_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'rsdt_gmean', field%rsdt_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%rsdt_gmean)

    cf_desc    = t_cf_var('rsut_gmean', 'W m-2', 'global mean toa outgoing shortwave radiation', datatype_flt,'rsut_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'rsut_gmean', field%rsut_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%rsut_gmean)

    cf_desc    = t_cf_var('rlut_gmean', 'W m-2', 'global mean toa outgoing longwave radiation', datatype_flt,'rlut_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'rlut_gmean', field%rlut_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%rlut_gmean)

    cf_desc    = t_cf_var('prec_gmean', 'kg m-2 s-1', 'global mean precipitation flux', datatype_flt,'prec_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'prec_gmean', field%prec_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%prec_gmean)

    cf_desc    = t_cf_var('evap_gmean', 'kg m-2 s-1', 'global mean evaporation flux', datatype_flt,'evap_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'evap_gmean', field%evap_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%evap_gmean)

!   derived variable
    cf_desc    = t_cf_var('radtop_gmean', 'W m-2', 'global mean toa net total radiation', datatype_flt,'radtop_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'radtop_gmean', field%radtop_gmean,      &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%radtop_gmean)

    !   derived variable
    cf_desc    = t_cf_var('radbot_gmean', 'W m-2', 'global mean surface net total radiation', datatype_flt,'radbot_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'radbot_gmean', field%radbot_gmean,      &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%radbot_gmean)

    !   derived variable
    cf_desc    = t_cf_var('radbal_gmean', 'W m-2', 'global mean net radiative flux into atmosphere', datatype_flt,'radbal_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'radbal_gmean', field%radbal_gmean,      &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%radbal_gmean)

    !   derived variable
    cf_desc    = t_cf_var('fwfoce_gmean', 'kg m-2 s-1', 'mean surface freshwater flux over ocean surface', &
                & datatype_flt,'fwfoce_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'fwfoce_gmean', field%fwfoce_gmean,      &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%fwfoce_gmean)

!   derived variable
    cf_desc    = t_cf_var('udynvi_gmean', 'J m-2', 'mean vertically integrated moist internal energy after dynamics', &
                & datatype_flt,'udynvi_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'udynvi_gmean', field%udynvi_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%udynvi_gmean)
!   derived variable
    cf_desc    = t_cf_var('duphyvi_gmean', 'J m-2', 'mean vertically integrated moist internal energy change by physics', &
                & datatype_flt,'duphyvi_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'duphyvi_gmean', field%duphyvi_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%duphyvi_gmean)
    !   derived variable
    IF (use_tmx .AND. is_variable_in_output(var_name=prefix//'utmxvi_gmean')) THEN
      cf_desc    = t_cf_var('utmxvi_gmean', 'J m-2', 'mean vertically integrated moist internal energy after tmx', &
                     & datatype_flt,'utmxvi_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( field_list, prefix//'utmxvi_gmean', field%utmxvi_gmean,          &
                     & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                     & lrestart = .FALSE., ldims=(/1/),                             &
                     & lopenacc=.TRUE.)
      __acc_attach(field%utmxvi_gmean)
    END IF
!   derived variable
    cf_desc    = t_cf_var('ufts_gmean', 'W m-2', 'mean energy flux at surface from thermal exchange', &
                & datatype_flt,'ufts_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'ufts_gmean', field%ufts_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%ufts_gmean)
!   derived variable
    cf_desc    = t_cf_var('ufvs_gmean', 'W m-2', 'mean energy flux at surface from vapor exchange', &
                & datatype_flt,'ufvs_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'ufvs_gmean', field%ufvs_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%ufvs_gmean)
!   derived variable
    cf_desc    = t_cf_var('ufcs_gmean', 'W m-2', 'mean energy flux at surface from condensate', &
                & datatype_flt,'ufcs_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'ufcs_gmean', field%ufcs_gmean,          &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%ufcs_gmean)
!   derived variable
    cf_desc    = t_cf_var('kedisp_gmean', 'W m-2', 'mean vert. integr. dissip. kin. energy', &
                & datatype_flt,'kedisp_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'kedisp_gmean', field%kedisp_gmean,      &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%kedisp_gmean)

!   derived variable
    cf_desc    = t_cf_var('uphybal_gmean', 'W m-2', 'mean energy balance in aes physics', &
                & datatype_flt,'uphybal_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'uphybal_gmean', field%uphybal_gmean,      &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                & lrestart = .FALSE., ldims=(/1/),                             &
                & lopenacc=.TRUE.)
    __acc_attach(field%uphybal_gmean)

! icefrc not allocated in atmosphere
!   cf_desc    = t_cf_var('icefrc_gmean', 'frac', 'global mean ice cover of grid box', datatype_flt,'icefrc_gmean')
!   grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
!   CALL add_var( field_list, prefix//'icefrc_gmean', field%icefrc_gmean,       &
!               & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                 &
!               & lrestart = .FALSE., ldims=(/1/),                              &
!               & lopenacc=.TRUE.)

  END SUBROUTINE new_aes_phy_field_list
  !-------------
  !>
  !!
  !!
  SUBROUTINE new_aes_phy_tend_list  ( jg, kproma, klev, kblks, ktracer, &
                                    & dt_dyn, listname, prefix,         &
                                    & tend_list, tend )
    INTEGER,INTENT(IN) :: jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer  !< dimension sizes
    TYPE(timedelta), INTENT(IN) :: dt_dyn !< Dynamics timestep.
    CHARACTER(*), INTENT(IN) :: listname, prefix
    TYPE(t_var_list_ptr), INTENT(INOUT) :: tend_list
    TYPE(t_aes_phy_tend), INTENT(INOUT) :: tend
    CHARACTER(len=vname_len) :: trc_name, cfstd_name, long_name, var_name, var_suffix
    LOGICAL :: contvar_is_in_output
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc
    INTEGER :: shape2d(2), shape3d(3), shape_trc(4), shape3d_layer_interfaces(3)
    INTEGER :: ibits, jtrc
    INTEGER :: datatype_flt
    !------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    shape2d   = (/kproma, kblks/)
    shape3d   = (/kproma, klev, kblks/)
    shape_trc = (/kproma, klev, kblks, ktracer/)
    shape3d_layer_interfaces = (/kproma,klev+1,kblks/)

    !$ACC ENTER DATA COPYIN(tend)

    CALL vlr_add(tend_list, listname, patch_id=jg ,lrestart=.FALSE., &
      &          model_type=get_my_process_name())

    !------------------------------
    ! Temperature tendencies
    !------------------------------
    ! &       tend% ta_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_phy', 'K s-1',                           &
                &         'temperature tendency due to model physics (cv)',              &
                &         datatype_flt)
    grib2_desc = grib2_var(0,0,210, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_phy', tend%  ta_phy,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ),                                                &
                & lopenacc=.TRUE.)
    __acc_attach(tend%  ta_phy)

    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       IF (is_variable_in_output(var_name=prefix//'ta_rsw')) THEN
          cf_desc    = t_cf_var('temperature_tendency_rsw', 'K s-1',                           &
                      &         'temperature tendency due to shortwave radiation (cp)',        &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,205, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_rsw', tend%  ta_rsw,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ),                                                &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%  ta_rsw)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'ta_rlw')) THEN
          cf_desc    = t_cf_var('temperature_tendency_rlw', 'K s-1',                           &
                      &         'temperature tendency due to longwave radiation (cp)',         &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,204, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_rlw', tend%  ta_rlw,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ),                                                &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%  ta_rlw)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'ta_rad')) THEN
          cf_desc    = t_cf_var('temperature_tendency_rad', 'K s-1',                           &
                      &         'temperature tendency due to radiation (cp)',                  &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_rad', tend%  ta_rad,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ),                                                &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%  ta_rad)
       END IF
       !
    END IF

    IF (is_variable_in_output(var_name=prefix//'ta_rlw_impl')) THEN
       cf_desc    = t_cf_var('temperature_tendency_rlw_impl', 'K s-1',                         &
                   &         'temperature tendency due to LW rad. due to implicit land surface temperature change (cp)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( tend_list, prefix//'ta_rlw_impl', tend%  ta_rlw_impl,                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                  &
                   & ldims=(/kproma,kblks/),                                                   &
                   & lopenacc=.TRUE.)
       __acc_attach(tend%  ta_rlw_impl)
    END IF

  !  IF ( aes_phy_tc(jg)%dt_mig > dt_zero ) THEN -> See:  mo_cloud_mig/mo_cloud_mig_memory
  !  IF ( aes_phy_tc(jg)%dt_two > dt_zero ) THEN -> See:  mo_cloud_two/mo_cloud_two_memory

    IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (is_variable_in_output(var_name=prefix//'ta_vdf')) THEN
          cf_desc    = t_cf_var('temperature_tendency_turbulent', 'K s-1',                     &
                      &         'temperature tendency due to vertical diffusion (cp)',         &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_vdf', tend%  ta_vdf,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ),                                                &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%  ta_vdf)
       END IF
       !
       IF (is_variable_in_output(var_name=prefix//'ta_sfc')) THEN
          cf_desc    = t_cf_var('temperature_tendency_surface',   'K s-1',                     &
                      &         'temperature tendency due to surface porcesses (cp)',          &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_sfc', tend%  ta_sfc,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                    &
                      & cf_desc, grib2_desc, ldims=shape2d,                                    &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%  ta_sfc)
       END IF
       !
    END IF

    !------------------------------
    ! U-wind tendencies
    !------------------------------
    ! &       tend%    ua_phy (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_phy', 'm s-2',                                &
                &         'u-wind tendency due to model physics',                        &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_phy', tend%ua_phy,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ),                                                &
                & lopenacc=.TRUE.)
    __acc_attach(tend%ua_phy)

    IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (aes_phy_tc(jg)%dt_vdf > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name=prefix//'ua_vdf')) THEN
          cf_desc    = t_cf_var('u_wind_tendency_turbulent', 'm s-2',                          &
                      &         'u-wind tendency due to vertical diffusion',                   &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ua_vdf', tend%ua_vdf,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ),                                                &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%ua_vdf)
       END IF
       !
    END IF

    !------------------------------
    ! V-wind tendencies
    !------------------------------
    ! &       tend%    va_phy (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_phy', 'm s-2',                                &
                &         'v-wind tendency due to model physics',                        &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,213, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_phy', tend%va_phy,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ),                                                &
                & lopenacc=.TRUE.)
    __acc_attach(tend%va_phy)

    IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (aes_phy_tc(jg)%dt_vdf > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name=prefix//'va_vdf')) THEN
          cf_desc    = t_cf_var('v_wind_tendency_turbulent', 'm s-2',                          &
                      &         'v-wind tendency due to vertical diffusion',                   &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,212, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'va_vdf', tend%va_vdf,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ),                                                &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%va_vdf)
       END IF
       !
    END IF

    !------------------------------
    ! W-wind tendencies
    !------------------------------
    ! &       tend%    wa_phy (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('tendency_of_upward_air_velocity', 'm s-2',                    &
                &         'tendency of upward air velocity due to param. processes',     &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'wa_phy', tend%wa_phy,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                             &
                & cf_desc,grib2_desc,                                                    &
                & ldims=shape3d_layer_interfaces,                                        &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ),                                                &
                & lopenacc=.TRUE.)
    __acc_attach(tend%wa_phy)

    IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (aes_phy_tc(jg)%dt_vdf > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name=prefix//'wa_vdf')) THEN
          cf_desc    = t_cf_var('tendency_of_upward_air_velocity', 'm s-2',                    &
                      &         'tendency of upward air velocity due to turbulent diffusion',  &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'wa_vdf', tend%wa_vdf,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                             &
                      & cf_desc,grib2_desc,                                                    &
                      & ldims=shape3d_layer_interfaces,                                        &
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ),                                                &
                      & lopenacc=.TRUE.)
          __acc_attach(tend%wa_vdf)
       END IF
       !
    END IF

    !-------------------
    ! Tracer tendencies
    !-------------------

    IF (ktracer > 0) THEN
       !
       ! Tendencies of mass fraction of tracers in air due to vertical diffusion
       !
       IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
          !
          contvar_is_in_output = .FALSE.
          DO jtrc = 1,ktracer
             trc_name = advection_config(jg)%tracer_names(jtrc)
             var_name = prefix//TRIM(trc_name)//'_vdf'
             IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
                contvar_is_in_output = .TRUE.
             END IF
          END DO
          !
          IF (aes_phy_tc(jg)%dt_vdf > dt_dyn .OR. contvar_is_in_output) THEN
             !
             var_name = prefix//'qtrc_vdf'
             cf_desc = t_cf_var( 'tendency_of_mass_fraction_of_tracer_in_air_due_to_vertical_diffusion', &
                               & 'kg kg-1 s-1',                                                          &
                               & 'tendency of mass fraction of tracer in air due to vertical diffusion', &
                               & datatype_flt                                                            )
             grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
             !
             CALL add_var( tend_list, TRIM(var_name), tend%qtrc_vdf, &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,     &
                         & cf_desc, grib2_desc,                      &
                         & ldims = shape_trc,                        &
                         & lcontainer=.TRUE.,                        &
                         & lrestart=.FALSE., loutput=.FALSE.,        &
                         & lopenacc=.TRUE. )
             __acc_attach(tend%qtrc_vdf)
          END IF
          !
          ! References for tendencies of mass fraction of tracer in air due to vertical diffusion, for output
          !
          ALLOCATE(tend%qtrc_vdf_ptr(ktracer))
          !
          DO jtrc = 1,ktracer
             !
             tend%qtrc_vdf_ptr(jtrc)%p => NULL()
             !
             ! preset names and grib2_desc for all tracers
             !
             trc_name   = TRIM(advection_config(jg)%tracer_names(jtrc))
             cfstd_name = TRIM(advection_config(jg)% cfstd_names(jtrc))
             long_name  = TRIM(advection_config(jg)%  long_names(jtrc))
             !
             grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
             !
             ! adjust names and grib2_desc for specific tracers 
             !
             IF (jtrc == iqv) THEN
                grib2_desc = grib2_var(0, 1, 202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
                !
             ELSE IF (jtrc == iqc) THEN
                grib2_desc = grib2_var(0, 6, 202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
                !
             ELSE IF (jtrc == iqi) THEN
                grib2_desc = grib2_var(0, 6, 212, ibits, GRID_UNSTRUCTURED, GRID_CELL)
                !
             END IF
             !
             ! finalize var_name and cf_desc
             !
             var_name = prefix//TRIM(trc_name)//'_vdf'
             cf_desc = t_cf_var( 'tendency_of_mass_fraction_of_'//TRIM(cfstd_name)//'_due_to_vertical_diffusion', &
                               & 'kg kg-1 s-1',                                                                   &
                               & 'tendency of mass fraction of '//TRIM(long_name)//' due to vertical diffusion',  &
                               & datatype_flt )
             !
             ! set memory references for fields which are requested for output
             !
             IF (is_variable_in_output(var_name=TRIM(var_name))) THEN
                CALL add_ref( tend_list, prefix//'qtrc_vdf',                &
                     & TRIM(var_name), tend%qtrc_vdf_ptr(jtrc)%p,           &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                &
                     & cf_desc, grib2_desc,                                 &
                     & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),           &
                     & lrestart = .FALSE.,                                  &
                     & vert_interp=create_vert_interp_metadata(             &
                     &             vert_intp_type=vintp_types("P","Z","I"), &
                     &             vert_intp_method=VINTP_METHOD_LIN )      )
             END IF
             !
          END DO
          !
       END IF

       ! Tendencies of mass fraction of ozone in air due to linearized ozone chemistry (Cariolle)
       !
       IF ( (aes_phy_tc(jg)%dt_car > dt_zero) ) THEN
          !
          trc_name   = TRIM(advection_config(jg)%tracer_names(io3))
          cfstd_name = TRIM(advection_config(jg)% cfstd_names(io3))
          long_name  = TRIM(advection_config(jg)%  long_names(io3))
          !
          var_name = prefix//TRIM(trc_name)//'_car'
          cf_desc = t_cf_var( 'tendency_of_mass_fraction_of_'//TRIM(cfstd_name)//'_due_to_linearized_ozone_chemistry_(Cariolle)', &
                            & 'kg kg-1 s-1',                                                                                      &
                            & 'tendency of mass fraction of '//TRIM(long_name)//' due to linearized ozone chemistry (Cariolle)',  &
                            & datatype_flt )
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL)
          !
          IF (aes_phy_tc(jg)%dt_car > dt_dyn .OR. is_variable_in_output(var_name=TRIM(var_name))) THEN
             CALL add_var( tend_list, var_name, tend%o3_car,     &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, &
                         & cf_desc, grib2_desc,                  &
                         & ldims = shape3d,                      &
                         & lrestart=.FALSE., loutput=.FALSE.,    &
                         & lopenacc=.TRUE. )
             __acc_attach(tend%o3_car)
          END IF
          !
       END IF

       ! Tendencies of mass fraction of tracers in air due to model physics
       !
       var_name = prefix//'qtrc_phy'
       cf_desc = t_cf_var( 'tendency_of_mass_fraction_of_tracer_in_air_due_to_model_physics', &
                         & 'kg kg-1 s-1',                                                     &
                         & 'tendency of mass fraction of tracer in air due to model physics', &
                         & datatype_flt )
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       !
       CALL add_var( tend_list, var_name, tend%qtrc_phy,   &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, &
                   & cf_desc, grib2_desc,                  &
                   & ldims = shape_trc,                    &
                   & lcontainer=.TRUE.,                    &
                   & lrestart=.FALSE., loutput=.FALSE.,    &
                   & lopenacc=.TRUE.                       )
       __acc_attach(tend%qtrc_phy)
       !
       ALLOCATE(tend%qtrc_phy_ptr(ktracer))
       !
       DO jtrc = 1,ktracer
          !
          tend%qtrc_phy_ptr(jtrc)%p => NULL()
          !
          ! preset names and grib2_desc for all tracers
          !
          trc_name   = TRIM(advection_config(jg)%tracer_names(jtrc))
          cfstd_name = TRIM(advection_config(jg)% cfstd_names(jtrc))
          long_name  = TRIM(advection_config(jg)%  long_names(jtrc))
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
          ! adjust names and grib2_desc for specific tracers 
          !
          IF (jtrc == iqv) THEN
             grib2_desc = grib2_var(0, 1, 108, ibits, GRID_UNSTRUCTURED, GRID_CELL)
             !
          END IF
          !
          ! set var_name and cf_desc
          !
          var_name = prefix//TRIM(trc_name)//'_phy'
          cf_desc = t_cf_var( 'tendency_of_mass_fraction_of_'//TRIM(cfstd_name)//'_in_air_due_to_model_physics', &
                            & 'kg kg-1 s-1',                                                                     &
                            & 'tendency of mass fraction of '//TRIM(long_name)//' in air due to model physics',  &
                            & datatype_flt )
          !
          ! set memory references for fields which are requested for output
          !
          IF ( is_variable_in_output(var_name=TRIM(var_name)) ) THEN
             CALL add_ref( tend_list, prefix//'qtrc_phy',                       &
                         & TRIM(var_name), tend%qtrc_phy_ptr(jtrc)%p,           &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                &
                         & cf_desc, grib2_desc,                                 &
                         & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),           &
                         & lrestart = .FALSE.,                                  &
                         & vert_interp=create_vert_interp_metadata(             &
                         &             vert_intp_type=vintp_types("P","Z","I"), &
                         &             vert_intp_method=VINTP_METHOD_LIN )      )
          END IF
          !
       END DO

       ! Array for tendencies of tracer paths due to model physics, for computations
       !
       var_name = prefix//'mtrcvi_phy'
       cf_desc = t_cf_var( 'tendency_of_atmosphere_mass_content_of_tracer_to_model_physics', &
                         & 'kg m-2 s-1',                                                     &
                         & 'tendency of tracer path due to model physics',                   &
                         & datatype_flt )
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       !
       CALL add_var( tend_list, var_name, tend%mtrcvi_phy, &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,   &
                   & cf_desc, grib2_desc,                  &
                   & ldims = (/kproma,kblks,ktracer/),     &
                   & lcontainer=.TRUE.,                    &
                   & lrestart=.FALSE., loutput=.FALSE.,    &
                   & lopenacc=.TRUE.                       )
       __acc_attach(tend%mtrcvi_phy)

       ! References for tendencies of tracer paths due to model physics, for output
       !
       ALLOCATE(tend%mtrcvi_phy_ptr(ktracer))
       !
       DO jtrc = 1,ktracer
          !
          tend%mtrcvi_phy_ptr(jtrc)%p => NULL()
          !
          ! preset names and grib2_desc for all tracers
          !
          trc_name   = TRIM(advection_config(jg)%tracer_names(jtrc))
          cfstd_name = TRIM(advection_config(jg)% cfstd_names(jtrc))
          long_name  = TRIM(advection_config(jg)%  long_names(jtrc))
          !
          var_suffix = 'vi'
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          !
          ! adjust names and grib2_desc for specific tracers 
          !
          IF (jtrc == iqv) THEN
             trc_name   = 'prw'
             var_suffix = ''
             cfstd_name = 'water_vapor'
             long_name  = 'water vapor'
             grib2_desc = grib2_var(0, 1, 224, ibits, GRID_UNSTRUCTURED, GRID_CELL)
             !
          ELSE IF (jtrc == iqc) THEN
             trc_name   = 'cll'
             grib2_desc = grib2_var(0, 1, 229, ibits, GRID_UNSTRUCTURED, GRID_CELL)
             !
          ELSE IF (jtrc == iqi) THEN
             grib2_desc = grib2_var(0, 1, 230, ibits, GRID_UNSTRUCTURED, GRID_CELL)
             !
          END IF
          !
          ! set var_name and cf_desc
          !
          var_name = prefix//TRIM(trc_name)//TRIM(var_suffix)//'_phy'
          cf_desc = t_cf_var( 'tendency_of_atmosphere_mass_content_of_'//TRIM(cfstd_name)//'_due_to_model_physics', &
                            & 'kg kg-1 s-1',                                                                        &
                            & 'tendency of '//TRIM(long_name)//' path due to model physics',                        &
                            & datatype_flt )
          !
          ! set memory references for fields which are requested for output
          !
          IF ( is_variable_in_output(var_name=TRIM(var_name)) ) THEN
             CALL add_ref( tend_list, prefix//'mtrcvi_phy',             &
                         & TRIM(var_name), tend%mtrcvi_phy_ptr(jtrc)%p, &
                         & GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE,       &
                         & cf_desc, grib2_desc,                         &
                         & ref_idx=jtrc, ldims=(/kproma,kblks/),        &
                         & lrestart = .FALSE. )
          END IF
          !
       END DO
       !
    END IF ! (ktracer > 0)

     ! &       tend% utmxvi     (nproma,nblks),          &
     IF (aes_vdf_config(jg)%use_tmx .AND. is_variable_in_output(var_name=prefix//'utmxvi')) THEN
      CALL add_var( tend_list, prefix//'utmxvi', tend%utmxvi,                  &
            &       GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE,                       &
            &       t_cf_var('u_phy_vi','J m-2 s-1',                             &
            &                'tendency of vertically integrated moist internal energy by tmx', &
            &                datatype_flt),                                      &
            &       grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
            &       ldims=[kproma,kblks],                                        &
            &       lrestart = .FALSE.,                                          &
            &       isteptype=TSTEP_INSTANT,                                     &
            &       lopenacc=.TRUE.)
      __acc_attach(tend%utmxvi)
    END IF

  END SUBROUTINE new_aes_phy_tend_list
  !-------------

END MODULE mo_aes_phy_memory
