!> Contains the memory class for the assimilation process.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>

MODULE mo_assimi_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_assimi_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 80

  !------------------------------------------------------------------------------------------------
  !> Type definition for assimi memory
  !! 
  !------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_assimi_memory
     ! t_canopy_sum variables
    TYPE(t_jsb_var_real2d) ::     &
      gross_assimilation_ca,      & ! Gross photosynthesis at each time step
                                    ! [mol(co2) / ( m^2(canopy ground) * s) ]  (mean value)
                                    ! Note, this is the gross assimilation on the ground below canopy area.
                                    ! This area is independent of the LAI. So even when the LAI is below 1 this area remains
                                    ! the same (it is still "covered").
      gross_assimilation,         & ! ...on tile area
                                    ! quincy: canopy-integrated gross assimilation rate [µmol CO2 m-2 s-1] 
      dark_respiration_ca,        & ! Dark respiration of leaf at each time step
                                    ! [mol(co2) / (m^2(canopy ground) * s) ]  (mean value)
      dark_respiration,           & ! see above but here on tile area
      NPP_pot_rate_ca,            & ! Net primary production rate [mol(CO2)/(m^2 (canopy ground) * s)]
      seconds_year,               & ! seconds till beginning of the current year
      NPP_sum_year,               & ! Net primary production sum of the current year
      NPP_mean_5year,             & ! Net primary production running average for 5 years
      land_cover_class,           & ! preliminary, to pass lct information on pft-tiles to a process running on the veg-tile
      bclimit_min_cold_mmtemp,    & ! preliminary, to pass lct information on pft-tiles to a process running on the veg-tile
      bclimit_max_cold_mmtemp,    & ! preliminary, to pass lct information on pft-tiles to a process running on the veg-tile
      bclimit_max_warm_mmtemp,    & ! preliminary, to pass lct information on pft-tiles to a process running on the veg-tile
      bclimit_min_temprange,      & ! preliminary, to pass lct information on pft-tiles to a process running on the veg-tile
      bclimit_min_gdd,            & ! preliminary, to pass lct information on pft-tiles to a process running on the veg-tile

      tau_c_woods                   ! preliminary, to pass lct information on pft-tiles to a process running on the veg-tile


!                                         CO2_conc_leaf,                  & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)]
!                                         CO2_conc_leaf_acc,              & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)]
!                                                                           ! (mean value)
!                                         CO2_conc_leaf_unlimited_acc,    & ! CO2 concentration inside leaf for no water stress
!                                                                           ! [MOL(CO2)/MOL(AIR)] (mean value)
!
!                                         carbox_rate_max,                & ! Maximum carboxylation rate (=VCmax)
!                                                                           ! [MOL(CO2)/M^2 S]
!                                         e_transport_rate_max,           & ! Maximum electron transport rate (=Jmax)
!                                                                           ! (only C3 plants)[MOL(CO2)/M^2 S]
!                                         carbox_rate,                    & ! Actual carboxylation rate (=JC)
!                                                                           ! [MOL(CO2)/M^2 S]
!                                         e_transport_rate,               & ! Actual electron transport rate (=JE)
!                                                                           ! [MOL(CO2)/M^2 S]
!                                         net_assimilation,               & ! Canopy net assimilation assimilation
!                                                                           ! [MOL(CO2)/m^2 s] (mean value)

!
!      REAL(wp), POINTER, DIMENSION(:,:) ::                               &
!                                         CO2_flux_net_acc,               & ! net CO2 flux to the atmosphere
!                                         CO2_flux_herbivory_acc,         & ! CO2 flux to the atmosphere from grazing
!                                         CO2_emission_landcover_change_acc, & ! CO2 flux to the atmosphere from land cover changes
!                                         CO2_emission_harvest_acc,       & ! CO2 flux to the atmosphere from harvest
!                                         CO2_flux_nlcc_acc                 ! CO2 flux to the atmosphere from vegetation dynamics
!                                                                           ! (disturbance?)
!


    ! t_canopy variables
    TYPE(t_jsb_var_real3d) ::     &
      & scaling_fact_cl,          &
      & canopy_cond_cl,           & ! Canopy conductance per layer and per leaf area
                                    ! quincy: conductance for water per canopy layer [m s-1]
      & gross_assimilation_cl,    & ! Gross photosynthesis [ mol(co2) / (m^2(leaf area) * s) ]
                                    ! Note, this is the gross assimilation on the leaf area that actually exists.
                                    ! The influence of the lai on the existing leaf area is aready included.
                                    ! It is not the ground area below canopy!
                                    ! quincy: gross assimilation rate per canopy layer [µmol CO2 m-2 s-1]
      & dark_respiration_cl,      & ! Dark respiration of leaf [ mol(co2) / (m^2(leaf area) * s) ]
      & carbox_rate_max_cl,       & ! [ mol(co2)/(m^2(leaf area) * s) ]
      & e_transport_rate_max_cl,  &
      & carbox_rate_cl,           &
      & e_transport_rate_cl,      &
      & CO2_conc_leaf_cl            ! CO2 concentration inside leaf [ mol(CO2)/ mol(AIR) ]
                                    ! quincy: internal leaf CO2 concentration per canopy layer [ppm]

    TYPE(t_jsb_var_real2d) ::     &
      & canopy_cond_unlimited,    & ! Canopy conductance without water limit  [m/s] (mean value)
      & canopy_cond_limited         ! Canopy conductance with water limit  [m/s] (mean value)

    ! R: This variable is needed for process DISTURB and should be later in the process vegdyn
    !    However it has to exist on each pft, therefore I put it to this process for now!
    TYPE(t_jsb_var_real2d) :: &
      & cover_fract_pot             ! Potential Natural Land Cover Fraction


    !!---------------------------------------------------------------------------------------------
    !! quincy - additional variables
    !! 
    ! fluxes
    TYPE(t_jsb_var_real2d)    :: &
                              gross_assimilation_C13    , &  !< canopy-integrated gross assimilation rate [µmol 13CO2 m-2 s-1]
                              gross_assimilation_C14    , &  !< canopy-integrated gross assimilation rate [µmol 14CO2 m-2 s-1]
                              net_assimilation          , &  !< canopy-integrated net assimilation rate [µmol CO2 m-2 s-1] 
                              net_assimilation_boc      , &  !< net assimilation rate at the bottom of the canopy [µmol CO2 m-2 s-1] 
                              maint_respiration_leaf         !< canopy-integrated foliar respiration rate [µmol CO2 m-2 s-1]

    TYPE(t_jsb_var_real3d)    :: &
                              net_assimilation_cl       , &  !< net assimilation rate per canopy layer [µmol CO2 m-2 s-1]
                              maint_respiration_leaf_cl      !< foliar respiration rate per canopy layer [µmol CO2 m-2 s-1]

    ! other 
    TYPE(t_jsb_var_real2d)    :: &
                              beta_air                  , &  !< scaling factor to account for air humidity effects on conductance [unitless]
                              beta_soa                  , &  !< scaling factor to account for spring in coniferous forests [unitless]
                              beta_soil_ps              , &  !< scaling factor to account for soil moisture constraints on photosynthesis [unitless]
                              beta_soil_gs              , &  !< scaling factor to account for soil moisture constraints on conductance [unitless]
                              canopy_cond               , &  !< canopy-integrated conductance for water [m s-1]
                              co2_conc_leaf             , &  !< canopy-integrated internal leaf CO2 concentration [ppm]
                              t_jmax_opt                     !< optimal temperature for electron transport of photosynthesis [°C]

    ! diagnostics
    TYPE(t_jsb_var_real3d)    :: &
                              jmax_cl                   , &  !< maximum rate of electron transport [µmol CO2 m-2 s-1]
                              vcmax_cl                       !< maximum rate of carboxylation [µmol CO2 m-2 s-1]

    ! averaging
    TYPE(t_jsb_var_real2d)    :: &
                              soa_tsoa_mavg                  !< state of acclimation, used in calculation of beta_soa time-averaged

    TYPE(t_jsb_var_real2d)    :: &
                              beta_air_daytime          , &  !< scaling factor to account for air humidity effects on conductance average of the previous day [unitless]
                              beta_air_daytime_dacc     , &  !< scaling factor to account for air humidity effects on conductance average of the previous day [unitless]
                              beta_air_tfrac_mavg       , &  !< scaling factor to account for air humidity effects on conductance time-averaged [unitless]
                              beta_air_tcnl_mavg        , &  !< scaling factor to account for air humidity effects on conductance time-averaged [unitless]
                              beta_soa_tphen_mavg            !< scaling factor to account for spring in coniferous forests time-averaged [unitless]

    TYPE(t_jsb_var_real2d)    :: &
                              beta_soil_ps_daytime      , &  !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                              beta_soil_ps_daytime_dacc , &  !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                              beta_soil_ps_tfrac_mavg   , &  !< scaling factor to account for soil moisture constraints on photosynthesis time-averaged [unitless]
                              beta_soil_ps_tcnl_mavg    , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                              beta_soil_gs_daytime      , &  !< scaling factor to account for soil moisture constraints on conductance average of the previous day [unitless]
                              beta_soil_gs_daytime_dacc , &  !< scaling factor to account for soil moisture constraints on conductance average of the previous day [unitless]
                              beta_soil_gs_tphen_mavg   , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                              beta_soil_gs_tfrac_mavg   , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                              beta_soil_gs_tcnl_mavg         !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]

  CONTAINS
    PROCEDURE :: Init => Init_assimi_memory
  END TYPE t_assimi_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_memory_class'

CONTAINS

  !------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: assimilation
  !!
  !!  this routine is called for all tiles
  !------------------------------------------------------------------------------------------------
  SUBROUTINE Init_assimi_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables !, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_assimi_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_canopy                 ! Vertical grid
    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_assimi_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    canopy_layer => Get_vgrid('canopy_layer')
    vgrid_canopy => Get_vgrid('canopy_layer_q_')  !! quincy
    model   => Get_model(model_id)

    IF (.NOT. model%config%use_quincy) THEN

      CALL mem%Add_var('gross_assimilation_ca', mem%gross_assimilation_ca,           &
        & hgrid, surface,                                                            &
        & t_cf('gross_assimilation_ca', 'mol(CO2)/m^2(canopy ground) / s',           &
        &      'Gross photosynthesis on canopy ground area'),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & output_level=BASIC,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL mem%Add_var('gross_assimilation', mem%gross_assimilation,                                                &
        & hgrid, surface,                                                                                           &
        & t_cf('gross_assimilation', 'mol(CO2)/m^2 (tile area) / s', 'Gross photosynthesis on tile area'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL mem%Add_var('dark_respiration_ca', mem%dark_respiration_ca,                                              &
        & hgrid, surface,                                                                                           &
        & t_cf('dark_respiration_ca', 'mol(CO2)/m^2(canopy ground) / s', 'Dark_respiration on canopy ground area'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL mem%Add_var('dark_respiration', mem%dark_respiration,                                                    &
        & hgrid, surface,                                                                                           &
        & t_cf('dark_respiration', 'mol(CO2)/m^2 (tile area) / s', 'Dark_respiration on tile area'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL mem%Add_var('npp_pot_rate_ca', mem%NPP_pot_rate_ca,                              &
        & hgrid, surface,                                                                   &
        & t_cf('NPP_pot_rate_ca', 'mol(CO2)/(m^2 (canopy ground) / s)', 'NPP_pot_rate_ca'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('seconds_year', mem%seconds_year,                                    &
        & hgrid, surface,                                                                   &
        & t_cf('seconds_year', 'sum of seconds in year', 'time sum year'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('npp_sum_year', mem%NPP_sum_year,                                    &
        & hgrid, surface,                                                                   &
        & t_cf('NPP_sum_year', 'mol(CO2)/(m^2 (canopy ground) / s)', 'NPP_sum_year'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.TRUE., initval_r=1.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('npp_mean_5year', mem%NPP_mean_5year,                                &
        & hgrid, surface,                                                                   &
        & t_cf('NPP_mean_5year', 'mol(CO2)/(m^2 (canopy ground) / s)', 'NPP_mean_5year'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.TRUE., initval_r=-1000.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('land_cover_class', mem%land_cover_class,                            &
        & hgrid, surface,                                                                   &
        & t_cf('land_cover_class', '', 'land_cover_class'),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('bclimit_min_cold_mmtemp', mem%bclimit_min_cold_mmtemp,              &
        & hgrid, surface,                                                                   &
        & t_cf('bclimit_min_cold_mmtemp', '', 'bclimit_min_cold_mmtemp'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('bclimit_max_cold_mmtemp', mem%bclimit_max_cold_mmtemp,              &
        & hgrid, surface,                                                                   &
        & t_cf('bclimit_max_cold_mmtemp', '', 'bclimit_max_cold_mmtemp'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('bclimit_max_warm_mmtemp', mem%bclimit_max_warm_mmtemp,              &
        & hgrid, surface,                                                                   &
        & t_cf('bclimit_max_warm_mmtemp', '', 'bclimit_max_warm_mmtemp'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('bclimit_min_temprange', mem%bclimit_min_temprange,                  &
        & hgrid, surface,                                                                   &
        & t_cf('bclimit_min_temprange', '', 'bclimit_min_temprange'),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('bclimit_min_gdd', mem%bclimit_min_gdd,                              &
        & hgrid, surface,                                                                   &
        & t_cf('bclimit_min_gdd', '', 'bclimit_min_gdd'),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('tau_c_woods', mem%tau_c_woods,                                      &
        & hgrid, surface,                                                                   &
        & t_cf('tau_c_woods', '', 'tau_c_woods'),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.FALSE., initval_r=10000.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('scaling_fact_cl', mem%scaling_fact_cl,               &
        & hgrid, canopy_layer,                                               &
        & t_cf('scaling_fact_cl', '-',                                       &
        'assimilation scaling factors per canopy layer'),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL mem%Add_var('canopy_cond_cl', mem%canopy_cond_cl,                 &
        & hgrid, canopy_layer,                                               &
        & t_cf('canopy_cond_cl', '-',                                        &
        'Canopy cond. per canopy layer and per leaf area'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL mem%Add_var('gross_assimilation_cl', mem%gross_assimilation_cl,                                                    &
        & hgrid, canopy_layer,                                                                                                &
        & t_cf('gross_assimilation_cl', 'mol(CO2)/m^2(leaf area) / s', 'Gross photosynthesis per canopy layer on leaf area'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                                  &
        & prefix, suffix,                                                                                                     &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL mem%Add_var('dark_respiration_cl', mem%dark_respiration_cl,                                                         &
        & hgrid, canopy_layer,                                                                                                 &
        & t_cf('dark_respiration_cl', 'mol(CO2)/m^2(leaf area) / s', 'Dark respiration of leaf per canopy layer on leaf area'),&
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                                   &
        & prefix, suffix,                                                                                                      &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL mem%Add_var('carbox_rate_max_cl', mem%carbox_rate_max_cl,                                                &
        & hgrid, canopy_layer,                                                                                      &
        & t_cf('carbox_rate_max_cl', 'mol(CO2)/m^2(leaf area) / s', 'Maximum carboxylation rate per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL mem%Add_var('e_transport_rate_max_cl', mem%e_transport_rate_max_cl,       &
        & hgrid, canopy_layer,                                                       &
        & t_cf('e_transport_rate_max_cl', 'dimensionen',                             &
        'Maximum electron transport rate per canopy layer'),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL mem%Add_var('carbox_rate_cl', mem%carbox_rate_cl,                                                   &
        & hgrid, canopy_layer,                                                                                 &
        & t_cf('carbox_rate_cl', 'mol(CO2)/m^2(leaf area) / s', 'Actual carboxylation rate per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                   &
        & prefix, suffix,                                                                                      &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL mem%Add_var('e_transport_rate_cl', mem%e_transport_rate_cl,                                                   &
        & hgrid, canopy_layer,                                                                                           &
        & t_cf('e_transport_rate_cl', 'mol(CO2)/m^2(leaf area) / s', 'Actual electron transport rate per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )! initval taken from JSBACH3

      CALL mem%Add_var('co2_conc_leaf_cl', mem%CO2_conc_leaf_cl,                                     &
        & hgrid, canopy_layer,                                                                       &
        & t_cf('CO2_conc_leaf_cl', 'dimensionen', 'CO2 concentration inside leaf per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
        & prefix, suffix,                                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3 CO2_conc_leaf_acc

      CALL mem%Add_var('canopy_cond_unlimited', mem%canopy_cond_unlimited,       &
        & hgrid, surface,                                                        &
        & t_cf('canopy_conductance_unlimited', '', ''),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & initval_r=0.0_wp )

      CALL mem%Add_var('canopy_cond_limited', mem%canopy_cond_limited,           &
        & hgrid, surface,                                                        &
        & t_cf('canopy_conductance_limited', '', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & initval_r=0.0_wp )

      CALL mem%Add_var('cover_fract_pot', mem%cover_fract_pot,               &
        & hgrid, surface,                                                    &
        & t_cf('cover_fract_pot', 'fract', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & initval_r=0.1_wp ) ! initval choosen freely
    END IF  !  IF(.NOT. use_quincy)

    !!---------------------------------------------------------------------------------------------
    !! quincy
    IF (model%config%use_quincy) THEN
      IF ( (     One_of(LAND_TYPE, lct_ids(:)) > 0 &
        &   .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0)) THEN
        
!!! temporary fix until the PHENO & other processes are using only the quincy tasks with use_quincy == TRUE
!!!  few jsb4 variables are needed by carbon_interface:update_C_NPP_pot_allocation
        CALL mem%Add_var('npp_pot_rate_ca', mem%NPP_pot_rate_ca,                              &
          & hgrid, surface,                                                                   &
          & t_cf('NPP_pot_rate_ca', 'mol(CO2)/(m^2 (canopy ground) / s)', 'NPP_pot_rate_ca'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
          & prefix, suffix,                                                                   &
          & output_level=BASIC,                                                               &
          & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var('gross_assimilation_ca', mem%gross_assimilation_ca,           &
          & hgrid, surface,                                                            &
          & t_cf('gross_assimilation_ca', 'mol(CO2)/m^2(canopy ground) / s',           &
          &      'Gross photosynthesis on canopy ground area'),                        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
          & prefix, suffix,                                                            &
          & output_level=BASIC,                                                        &
          & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3       

          ! NOTE this jsb4 var is identical to what the canopy_cond is in quincy
          CALL mem%Add_var('canopy_cond_limited', mem%canopy_cond_limited,         &
          & hgrid, surface,                                                        &
          & t_cf('canopy_conductance_limited', '', ''),                            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & initval_r=0.0_wp )

!!! temp fix end
                
        CALL mem%Add_var('gross_assimilation', mem%gross_assimilation, &
              hgrid, surface, &
              t_cf('gross_assimilation', 'µmol CO2 m-2 s-1', 'canopy-integrated gross assimilation rate'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('gross_assimilation_C13', mem%gross_assimilation_C13, &
              hgrid, surface, &
              t_cf('gross_assimilation_C13', 'µmol 13CO2 m-2 s-1', 'canopy-integrated gross assimilation rate'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('gross_assimilation_C14', mem%gross_assimilation_C14, &
              hgrid, surface, &
              t_cf('gross_assimilation_C14', 'µmol 14CO2 m-2 s-1', 'canopy-integrated gross assimilation rate'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('net_assimilation', mem%net_assimilation, &
              hgrid, surface, &
              t_cf('net_assimilation', 'µmol CO2 m-2 s-1', 'canopy-integrated net assimilation rate'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('net_assimilation_boc', mem%net_assimilation_boc, &
              hgrid, surface, &
              t_cf('net_assimilation_boc', 'µmol CO2 m-2 s-1', 'net assimilation rate at the bottom of the canopy'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('maint_respiration_leaf', mem%maint_respiration_leaf, &
              hgrid, surface, &
              t_cf('maint_respiration_leaf', 'µmol CO2 m-2 s-1', 'canopy-integrated foliar respiration rate'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('gross_assimilation_cl', mem%gross_assimilation_cl, &
              hgrid, vgrid_canopy, &
              t_cf('gross_assimilation_cl', 'µmol CO2 m-2 s-1', 'gross assimilation rate per canopy layer'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('net_assimilation_cl', mem%net_assimilation_cl, &
              hgrid, vgrid_canopy, &
              t_cf('net_assimilation_cl', 'µmol CO2 m-2 s-1', 'net assimilation rate per canopy layer'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('maint_respiration_leaf_cl', mem%maint_respiration_leaf_cl, &
              hgrid, vgrid_canopy, &
              t_cf('maint_respiration_leaf_cl', 'µmol CO2 m-2 s-1', 'foliar respiration rate per canopy layer'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('beta_air', mem%beta_air, &
              hgrid, surface, &
              t_cf('beta_air', 'unitless', 'scaling factor to account for air humidity effects on conductance'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
   
        CALL mem%Add_var('beta_soa', mem%beta_soa, &
              hgrid, surface, &
              t_cf('beta_soa', 'unitless', 'scaling factor to account for spring in coniferous forests'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('beta_soil_ps', mem%beta_soil_ps, &
              hgrid, surface, &
              t_cf('beta_soil_ps', 'unitless', 'scaling factor to account for soil moisture constraints on photosynthesis'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_gs', mem%beta_soil_gs, &
              hgrid, surface, &
              t_cf('beta_soil_gs', 'unitless', 'scaling factor to account for soil moisture constraints on conductance'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('canopy_cond', mem%canopy_cond, &
              hgrid, surface, &
              t_cf('canopy_cond', 'm s-1', 'canopy-integrated conductance for water'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('co2_conc_leaf', mem%co2_conc_leaf, &
              hgrid, surface, &
              t_cf('co2_conc_leaf', 'ppm', 'canopy-integrated internal leaf CO2 concentration'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('t_jmax_opt', mem%t_jmax_opt, &
              hgrid, surface, &
              t_cf('t_jmax_opt', 'deg C', 'optimal temperature for electron transport of photosynthesis'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
   
        CALL mem%Add_var('canopy_cond_cl', mem%canopy_cond_cl, &
              hgrid, vgrid_canopy, &
              t_cf('canopy_cond_cl', 'm s-1', 'conductance for water per canopy layer'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
   
        CALL mem%Add_var('co2_conc_leaf_cl', mem%co2_conc_leaf_cl, &
              hgrid, vgrid_canopy, &
              t_cf('co2_conc_leaf_cl', 'ppm', 'internal leaf CO2 concentration per canopy layer'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
   
        CALL mem%Add_var('jmax_cl', mem%jmax_cl, &
              hgrid, vgrid_canopy, &
              t_cf('jmax_cl', 'µmol CO2 m-2 s-1', 'maximum rate of electron transport'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
   
        CALL mem%Add_var('vcmax_cl', mem%vcmax_cl, &
              hgrid, vgrid_canopy, &
              t_cf('vcmax_cl', 'µmol CO2 m-2 s-1', 'maximum rate of carboxylation'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.FALSE.,&
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('soa_tsoa_mavg', mem%soa_tsoa_mavg, &
              hgrid, surface, &
              t_cf('soa_tsoa_mavg', '', 'state of acclimation, used in calculation of beta_soa time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('beta_air_daytime', mem%beta_air_daytime, &
              hgrid, surface, &
              t_cf('beta_air_daytime', 'unitless', 'scaling factor beta_air average of the previous day'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_air_daytime_dacc', mem%beta_air_daytime_dacc, &
              hgrid, surface, &
              t_cf('beta_air_daytime_dacc', 'unitless', 'scaling factor beta_air average of the previous day'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_air_tfrac_mavg', mem%beta_air_tfrac_mavg, &
              hgrid, surface, &
              t_cf('beta_air_tfrac_mavg', 'unitless', 'scaling factor beta_air time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_air_tcnl_mavg', mem%beta_air_tcnl_mavg, &
              hgrid, surface, &
              t_cf('beta_air_tcnl_mavg', 'unitless', 'scaling factor beta_air time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('beta_soa_tphen_mavg', mem%beta_soa_tphen_mavg, &
              hgrid, surface, &
              t_cf('beta_soa_tphen_mavg', 'unitless', 'scaling factor account for spring in coniferous forests time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
  
        CALL mem%Add_var('beta_soil_ps_daytime', mem%beta_soil_ps_daytime, &
              hgrid, surface, &
              t_cf('beta_soil_ps_daytime', 'unitless', 'scaling factor beta_soil_ps average of the previous day'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_ps_daytime_dacc', mem%beta_soil_ps_daytime_dacc, &
              hgrid, surface, &
              t_cf('beta_soil_ps_daytime_dacc', 'unitless', 'scaling factor beta_soil_ps average of the previous day'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_ps_tfrac_mavg', mem%beta_soil_ps_tfrac_mavg, &
              hgrid, surface, &
              t_cf('beta_soil_ps_tfrac_mavg', 'unitless', 'scaling factor beta_soil_ps time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_ps_tcnl_mavg', mem%beta_soil_ps_tcnl_mavg, &
              hgrid, surface, &
              t_cf('beta_soil_ps_tcnl_mavg', 'unitless', 'scaling factor beta_soil_ps time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_gs_daytime', mem%beta_soil_gs_daytime, &
              hgrid, surface, &
              t_cf('beta_soil_gs_daytime', 'unitless', 'scaling factor beta_soil_gs average of the previous day'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_gs_daytime_dacc', mem%beta_soil_gs_daytime_dacc, &
              hgrid, surface, &
              t_cf('beta_soil_gs_daytime_dacc', 'unitless', 'scaling factor beta_soil_gs average of the previous day'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_gs_tphen_mavg', mem%beta_soil_gs_tphen_mavg, &
              hgrid, surface, &
              t_cf('beta_soil_gs_tphen_mavg', 'unitless', 'scaling factor beta_soil_gs time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_gs_tfrac_mavg', mem%beta_soil_gs_tfrac_mavg, &
              hgrid, surface, &
              t_cf('beta_soil_gs_tfrac_mavg', 'unitless', 'scaling factor beta_soil_gs time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
    
        CALL mem%Add_var('beta_soil_gs_tcnl_mavg', mem%beta_soil_gs_tcnl_mavg, &
              hgrid, surface, &
              t_cf('beta_soil_gs_tcnl_mavg', 'unitless', 'scaling factor beta_soil_gs time-averaged'), &
              t_grib1(table, 255, grib_bits), &
              t_grib2(255, 255, 255, grib_bits), &
              prefix, suffix, &
              loutput=.TRUE., lrestart=.TRUE., &
              initval_r=0.0_wp) ! initvalue
      END IF            
    END IF  ! use_quincy

  END SUBROUTINE Init_assimi_memory

#endif
END MODULE mo_assimi_memory_class
