!> Contains the memory class for the carbon process.
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
MODULE mo_carbon_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_cqt_class,          ONLY: LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes PHENO_, CARBON_

  ! Use process configurations
  dsl4jsb_Use_config(PHENO_)    ! needed for l_forestRegrowth
  dsl4jsb_Use_config(CARBON_)    ! needed for l_forestRegrowth

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_carbon_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 155

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_carbon_memory
    ! 1d reals:
    !REAL(wp) ::            &
    TYPE(t_jsb_var_real2d) ::     &
       N_pseudo_temp,       & ! Constants used for updating the pseudo-15-day mean temperature
       F_pseudo_temp,       & ! Constants used for updating the pseudo-15-day mean temperature
       N_pseudo_precip,     & ! Constants used for updating the pseudo-15-day mean precipitation
       F_pseudo_precip        ! Constants used for updating the pseudo-15-day mean precipitation

    ! cbalance_type variables
    TYPE(t_jsb_var_real2d) ::     &
      pseudo_temp,                & ! a proxy for 15-day mean air temperature, [K] needed with yasso
      pseudo_temp_yDay,           & !
      pseudo_precip,              & ! a proxy for 15-day mean preciptation rate [Kg/m^2 s] needed with yasso
      pseudo_precip_yDay            !

    TYPE(t_jsb_var_real2d) ::     &
      cconservation_calcCpools,   & ! Carbon conservation test: Deviation from conservation [mol(C)/m^2(grid box)]
      cflux_c_greenwood_2_litter, & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
                                    !    [mol(C)/m^2(canopy) s]
      c_green,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      c_woods,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      c_reserve,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive
                                    ! .. bad times[mol(C)/m^2(canopy)]
!@todo why is the c_crop_harvest realised as a density and not as a per tile variable in the first place?
! ++ only required for crop pfts, or?
      c_crop_harvest,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
     !SIZE CLASS 1
      c_acid_ag1,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_ag1,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_ag1,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_ag1,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_acid_bg1,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_bg1,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_bg1,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_bg1,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_humus_1,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(canopy)]
     !SIZE CLASS 2
      c_acid_ag2,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_ag2,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_ag2,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_ag2,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_acid_bg2,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_bg2,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_bg2,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_bg2,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_humus_2,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(canopy)]

      c_bg_sum,              & ! Sum of belowground Yasso pools  [mol(C)/m^2(canopy)]
      c_ag_sum_1,            & ! Sum of aboveground green Yasso pools  [mol(C)/m^2(canopy)]
      max_green_bio,         & ! Annual maximum of C-pool green for NLCC process [mol(C)/m^2(canopy)]
      sla,                   & ! specific leave area for NLCC process [m2(leaf)/mol(C)]

  !JN: vars for analytically equilibrating YASSO humus pools outside of jsb4 simulations (compare jsb3 svn ref 9652)
      c_decomp_humus_1_sum,      & ! annual cflux sum of humus decomposed in YASSO (leaf) [mol(C)/m^2(canopy)]
      c_decomp_humus_2_sum,      & ! annual cflux sum of humus decomposed in YASSO (wood) [mol(C)/m^2(canopy)]
      c_into_humus_1_sum,        & ! annual cflux sum into humus pool in YASSO (leaf) [mol(C)/m^2(canopy)]
      c_into_humus_2_sum,        & ! annual cflux sum into humus pool in YASSO (wood) [mol(C)/m^2(canopy)]

      !
      LAI_sum,                    & ! used to accumulate LAI over a day
      NPP_pot_sum,                    & ! used to accumulated NPP-Rate over a day
      GPP_sum,                    & ! used to accumulated GPP-Rate over a day
      !runoff_sum,                 & ! used to accumulate runoff over a day.  R: Nicht notwendig in JS4.
      soil_respiration,           & ! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(canopy)].
                                    ! Without N limitation!
      NPP_flux_correction,        & ! Daily updated flux correction from yesterdays carbon balance
                                    ! [mol(CO2)/m^2(canopy) s]
                                    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                    ! correction arises either because otherwise the reserve pool would get
                                    ! negative (positive correction).
      excess_NPP,                 & ! That part of NPP that because of structural limits could not be allocated in
                                    ! carbon pools [mol(CO2)/m^2(canopy) s]
      LAI_yyDayMean,              & ! previous previous day mean
      LAI_yDayMean,               & ! previous day mean
      NPP_pot_yDayMean,           & ! mean value of NPP-Rate yesterday (from NPP_pot_sum()) [mol(CO2)/(m^2(canopy) s)]
      NPP_act_yDayMean,           & ! mean value of actual NPP-Rate after N-limitation and excess carbon drop.
                                    ! "yDayMean" because the actual NPP depends on/is calculated from results of the
                                    ! previous day.
      GPP_yDayMean,               & !
      root_exudates,              & ! Total root exudates entering to the litter green pools [mol(C)/m^2(canopy) s]
      cflux_c_green_2_herb,            & ! Total carbon flux from C pool green to herbivory.
      cflux_herb_2_littergreen,         & ! ..part of cflux_c_green_2_herb that goes into green litter
      cflux_herb_2_atm,      & ! ..part of cflux_c_green_2_herb that goes into atmosphere
      !
      cflux_dist_greenreserve_2_soil,  & ! R: should be later renamed to C_2_GreenLitterPools
      cflux_dist_woods_2_soil,    & ! R: should be later renamed to C_2_WoodLitterPools
      cflux_fire_all_2_atm,       & ! R: should be later renamed to C_2_atmos
                                    ! Carbon immedieatly released by fire per tile [mol(C)m-2(subtile)s-1]
      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      !fract_litter_wood_new,       & ! new fraction in above ground wood litter pool

  !JN: variables required when l_forestRegrowth = true, i.e. when maxLAI is derived from biomass via allometric relationships
      current_max_green,          & ! Running maximum of the green pool [mol(C) m-2(canopy)]
                                    ! since last maxLAI calculation
      veg_carbon_at_max_green       ! Total vegetation carbon [mol(C) m-2(canopy)]
                                    ! at the time of the green pool maximum (current_max_green)

   ! Output on tile areas:
   TYPE(t_jsb_var_real2d) ::         &
      cflux_c_greenwood_2_litter_ta, & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
                                       !    [mol(C)/m^2(tile area) s]
      c_green_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      c_woods_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      c_reserve_ta,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive
                                       ! .. bad times[mol(C)/m^2(tile area)]
      c_crop_harvest_ta,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
      !CBAL c_slow_ta,                 & ! C-pool for below ground organic material (coming from the litter pools)
      !CBAL                              !     [mol(C)/m^2(tile area)]
     !SIZE CLASS 1
      c_acid_ag1_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_ag1_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_ag1_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_ag1_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_acid_bg1_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_bg1_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_bg1_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_bg1_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_humus_1_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(tile area)]
     !SIZE CLASS 2
      c_acid_ag2_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_ag2_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_ag2_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_ag2_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_acid_bg2_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_bg2_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_bg2_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_bg2_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_humus_2_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(tile area)]

  !JN: vars for analytically equilibrating YASSO humus pools outside of jsb4 simulations (compare jsb3 svn ref 9652)
      c_decomp_humus_1_sum_ta,  & ! annual cflux sum of humus decomposed in YASSO (leaf) [mol(C)/m^2(tile area)]
      c_decomp_humus_2_sum_ta,  & ! annual cflux sum of humus decomposed in YASSO (wood) [mol(C)/m^2(tile area)]
      c_into_humus_1_sum_ta,    & ! annual cflux sum into humus pool in YASSO (leaf) [mol(C)/m^2(tile area)]
      c_into_humus_2_sum_ta,    & ! annual cflux sum into humus pool in YASSO (wood) [mol(C)/m^2(tile area)]

      !
      ! diagnostic carbon sums
      c_sum_veg_ta,             & ! Sum of all vegetation carbon pools [mol(C)/m^2(tile area)]
      c_sum_litter_ag_ta,       & ! Sum of all above ground litter carbon pools [mol(C)/m^2(tile area)]
      c_sum_litter_bg_ta,       & ! Sum of all below ground litter carbon pools [mol(C)/m^2(tile area)]
      c_sum_humus_ta,           & ! Sum of all soil / humus carbon pools [mol(C)/m^2(tile area)]
      c_sum_natural_ta,         & ! Sum of all natural carbon pools [mol(C)/m^2(tile area)]
      !
      !@todo: NPP_pot_sum_ta, GPP_sum_ta, LAI_yyDayMean_ta and LAI_yDayMean_ta are never calculated
!      NPP_pot_sum_ta,                    & ! used to accumulated NPP-Rate over a day
!      GPP_sum_ta,                    & ! used to accumulated GPP-Rate over a day
      soil_respiration_ta,           & ! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(tile area)].
                                       ! Without N limitation!
      NPP_flux_correction_ta,        & ! Daily updated flux correction from yesterdays carbon balance
                                    ! [mol(CO2)/m^2(tile area) s]
                                    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                    ! correction arises either because otherwise the reserve pool would get
                                    ! negative (positive correction).
      excess_NPP_ta,                 & ! That part of NPP that because of structural limits could not be allocated in
                                    ! carbon pools [mol(CO2)/m^2(tile area) s]
!      LAI_yyDayMean_ta,              & ! previous previous day mean
!      LAI_yDayMean_ta,               & ! previous day mean
      NPP_pot_yDayMean_ta,           & ! mean value of NPP-Rate yesterday (from NPP_pot_sum()) [mol(CO2)/(m^2(tile area) s)]
      NPP_act_yDayMean_ta,           & ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation.
                                       ! Actual NPP (after N-limitation and) excess carbon drop.
                                       ! [mol(CO2)/(m^2(tile area) s)]
      GPP_yDayMean_ta,               & !
      root_exudates_ta,              & ! Total root exudates entering to the litter green pools [mol(C)/m^2(tile area) s]
      cflux_c_green_2_herb_ta,            & ! Total carbon flux from C pool green to herbivory.
      !cflux_herb_2_atm_ta,      & ! ..part of cflux_c_green_2_herb that goes into atmosphere
      cflux_dist_greenreserve_2_soil_ta,  & ! R: should be later renamed to C_2_GreenLitterPools_ta
      cflux_dist_woods_2_soil_ta,   & ! R: should be later renamed to C_2_WoodLitterPools_ta

      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      fract_litter_wood_new_ta,      & ! new fraction in above ground wood litter pool

      ! tile averages of net CO2 fluxes between various compartments and atmosphere [kg(CO2)/m^2 s]  
      ! .. positive for emission to atmosphere
      ! .. negative for absorption by biosphere
      co2flux_npp_2_atm_ta,          & ! fluxes due to NPP
      co2flux_soilresp_2_atm_ta,     & ! fluxes due to soil respiration
      co2flux_herb_2_atm_ta,         & ! fluxes due to grazing
      co2flux_fire_all_2_atm_ta,     & ! Carbon immediatly released by fire per tile
      co2flux_npp_2_atm_yday_ta        ! previous day CO2 flux from actual NPP, required for C conservation test


    ! Diagnostic global carbon sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) ::        &
      c_sum_veg_gsum,                & ! Global sum of all vegetation carbon pools [PgC]
      c_sum_litter_ag_gsum,          & ! Global sum of all above ground litter carbon pools [PgC]
      c_sum_litter_bg_gsum,          & ! Global sum of all below ground litter carbon pools [PgC]
      c_sum_humus_gsum,              & ! Global sum of all soil / humus carbon pools [PgC]
      c_sum_natural_gsum,            & ! Global sum of all natural carbon pools [PgC]
      NPP_act_yDayMean_gsum,         & ! Global actual NPP rate [PgC/yr]
      GPP_yDayMean_gsum,             & ! Global GPP rate [PgC/yr]
      soil_respiration_gsum


  CONTAINS
    PROCEDURE :: Init => Init_carbon_memory
  END TYPE t_carbon_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_memory_class'

CONTAINS

  SUBROUTINE Init_carbon_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_carbon_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(PHENO_)    ! needed for l_forestRegrowth
    dsl4jsb_Def_config(CARBON_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid     ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface   ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_memory'
    IF (lct_ids(1) > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_config(CARBON_)

! R: only if we would want to have these variables in 2d (see mo_carbon_init.f90 for their initialization):
    CALL mem%Add_var('n_pseudo_temp', mem%N_pseudo_temp,                                             &
      & hgrid, surface,                                                                              &
      & t_cf('N_pseudo_temp', 'K', 'Normalization for computing pseudo-15-day mean air temperature'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var('f_pseudo_temp', mem%F_pseudo_temp,                                             &
      & hgrid, surface,                                                                              &
      & t_cf('F_pseudo_temp', '-', 'Exponential factor used for updating the pseudo-15-day mean temperature'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var('n_pseudo_precip', mem%N_pseudo_precip,                                         &
      & hgrid, surface,                                                                              &
      & t_cf('N_pseudo_precip', 'Kg m-2 s-1', 'Normalization for computing pseudo-15-day mean precipitation rate'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('f_pseudo_precip', mem%F_pseudo_precip,                                         &
      & hgrid, surface,                                                                              &
      & t_cf('F_pseudo_precip', '-', 'Exponential factor used for updating the pseudo-15-day mean precip'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp )


    CALL mem%Add_var('pseudo_temp', mem%pseudo_temp,                                                 &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_temp', 'K', 'Pseudo-mean air temperature over previous 15 days'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & output_level=MEDIUM,                                                                         &
      & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

    ! R: In the output it should be the same as pseudo_temp, I do not understand why it was written out in JS3.
    CALL mem%Add_var('pseudo_temp_yDay', mem%pseudo_temp_yDay,                                         &
      & hgrid, surface,                                                                                &
      & t_cf('pseudo_temp_yDay', 'K', 'Yesterdays pseudo-mean air temperature over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

    CALL mem%Add_var('pseudo_precip', mem%pseudo_precip,                                             &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_precip', 'Kg m-2 s-1', 'Pseudo-mean precipitation rate over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & output_level=MEDIUM,                                                                         &
      & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

    ! R: In the output it should be the same as pseudo_precip, I do not understand why it was written out in JS3.
    CALL mem%Add_var('pseudo_precip_yDay', mem%pseudo_precip_yDay,                                   &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_precip_yDay', 'Kg m-2 s-1', 'Yesterdays pseudo-mean precipitation rate over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

    ! Carbon variables on canopy area
    CALL mem%Add_var('c_greenwood2litter', mem%cflux_c_greenwood_2_litter,                             &
      & hgrid, surface,                                                                                &
      & t_cf('c_greenwood2litter', 'mol(CO2) m-2(canopy) s-1', 'Total litter flux entering the soil pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL,                                                                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_green', mem%c_green,                                                   &
      & hgrid, surface,                                                                        &
      & t_cf('c_green', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_woods', mem%c_woods,                                                              &
      & hgrid, surface,                                                                                   &
      & t_cf('c_woods', 'mol(C) m-2(canopy)', 'C-Pool for structural material/parts of vegetation.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
      & prefix, suffix,                                                                                   &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,               &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_reserve', mem%c_reserve,                                                                           &
      & hgrid, surface,                                                                                                    &
      & t_cf('c_reserve', 'mol(C) m-2(canopy)', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,            &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_crop_harvest', mem%c_crop_harvest,                                             &
      & hgrid, surface,                                                                                &
      & t_cf('c_crop_harvest', 'mol(C) m-2(canopy)', 'C-Pool for Harvested Material from Crops.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = PRODUCT_CARBON_CQ_TYPE,         &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 1
    CALL mem%Add_var('c_acid_ag1', mem%c_acid_ag1,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_ag1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag1', mem%c_water_ag1,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_ag1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag1', mem%c_ethanol_ag1,                                                               &
      & hgrid, surface,                                                                                                &
      & t_cf('c_ethanol_ag1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag1', mem%c_nonsoluble_ag1,                                                        &
      & hgrid, surface,                                                                                               &
      & t_cf('c_nonsoluble_ag1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg1', mem%c_acid_bg1,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_bg1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg1', mem%c_water_bg1,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_bg1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                       &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg1', mem%c_ethanol_bg1,                                                              &
      & hgrid, surface,                                                                                               &
      & t_cf('c_ethanol_bg1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg1', mem%c_nonsoluble_bg1,                                                     &
      & hgrid, surface,                                                                                            &
      & t_cf('c_nonsoluble_bg1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_1', mem%c_humus_1,                                     &
      & hgrid, surface,                                                              &
      & t_cf('c_humus_1', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 2
    CALL mem%Add_var('c_acid_ag2', mem%c_acid_ag2,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_ag2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag2', mem%c_water_ag2,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_ag2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag2', mem%c_ethanol_ag2,                                                               &
      & hgrid, surface,                                                                                                &
      & t_cf('c_ethanol_ag2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                            &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag2', mem%c_nonsoluble_ag2,                                                        &
      & hgrid, surface,                                                                                               &
      & t_cf('c_nonsoluble_ag2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg2', mem%c_acid_bg2,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_bg2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg2', mem%c_water_bg2,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_bg2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg2', mem%c_ethanol_bg2,                                                              &
      & hgrid, surface,                                                                                               &
      & t_cf('c_ethanol_bg2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg2', mem%c_nonsoluble_bg2,                                                       &
      & hgrid, surface,                                                                                              &
      & t_cf('c_nonsoluble_bg2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_2', mem%c_humus_2,                                     &
      & hgrid, surface,                                                              &
      & t_cf('c_humus_2', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_bg_sum', mem%c_bg_sum,                                         &
      & hgrid, surface,                                                                &
      & t_cf('c_bg_sum', 'mol(C) m-2(canopy)', 'C-Pool of below ground C in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                                &
      & output_level=FULL,                                                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_ag_sum_1', mem%c_ag_sum_1,                                             &
      & hgrid, surface,                                                                        &
      & t_cf('c_ag_sum_1', 'mol(C) m-2(canopy)', 'C-Pool of above ground green C in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL,                                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('max_green_bio', mem%max_green_bio,                                       &
      & hgrid, surface,                                                                        &
      & t_cf('max_green_bio', 'mol(C) m-2(canopy)', 'annual maximum content of green C-Pool.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL,                                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('sla', mem%sla,                                                           &
      & hgrid, surface,                                                                        &
      & t_cf('sla', 'm2(leave)/mol(C)', 'specific leave area.'),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL,                                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('lai_sum', mem%LAI_sum,                                         &
      & hgrid, surface,                                                              &
      & t_cf('LAI_sum', '-', 'Used to accumulate LAI over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_sum', mem%NPP_pot_sum,                                         &
      & hgrid, surface,                                                                      &
      & t_cf('NPP_pot_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate NPP over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.) ! initval taken from JSBACH3

    CALL mem%Add_var('gpp_sum', mem%GPP_sum,                                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate GPP over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('soil_respiration', mem%soil_respiration,                 &
      & hgrid, surface,                                                        &
      & t_cf('soil_respiration', 'mol(C) m-2(canopy) s-1', 'Soil respiration.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_flux_correction', mem%NPP_flux_correction,                                 &
      & hgrid, surface,                                                                              &
      & t_cf('NPP_flux_correction', 'mol(C) m-2(canopy) s-1',                                        &
      &      'Amount by which the NPP rate has to be corrected to avoid a negative reserve pool.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('excess_npp', mem%excess_NPP,                                                                     &
      & hgrid, surface,                                                                                                &
      & t_cf('excess_NPP', 'mol(C) m-2 (canopy) s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('lai_yydaymean', mem%LAI_yyDayMean,                    &
      & hgrid, surface,                                                     &
      & t_cf('LAI_yyDayMean', '-', 'Previous previous day mean.'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('lai_ydaymean', mem%LAI_yDayMean,                      &
      & hgrid, surface,                                                     &
      & t_cf('LAI_yDayMean', '-', 'Previous day mean.'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_ydaymean', mem%NPP_pot_yDayMean,                                    &
      & hgrid, surface,                                                                           &
      & t_cf('NPP_pot_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & lrestart=.TRUE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would 
                                                  ! result in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('npp_act_ydaymean', mem%NPP_act_yDayMean,                                      &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean', 'mol(C) m-2(canopy) s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & lrestart=.TRUE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would 
                                                  ! result in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean', mem%GPP_yDayMean,                                        &
      & hgrid, surface,                                                                       &
      & t_cf('GPP_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
      & prefix, suffix,                                                                       &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('root_exudates', mem%root_exudates,                   &
      & hgrid, surface,                                                    &
      & t_cf('root_exudates', 'mol(C) m-2(canopy) s-1', 'Root exudates'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_c_green_2_herb', mem%cflux_c_green_2_herb,                                          &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_c_green_2_herb', 'mol(C) m-2(canopy) s-1', 'Total carbon flux from C pool green to herbivory.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_c_green_2_herb_lg', mem%cflux_herb_2_littergreen,                                   &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_herb_2_littergreen', 'mol(C) m-2(canopy) s-1',                                              &
      &      'Carbon flux from cflux_c_green_2_herb into green litter (complementary to cflux_herb_2_atm).'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_herb_2_atm', mem%cflux_herb_2_atm,                                                  &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_herb_2_atm', 'mol(C) m-2(canopy) s-1',                                                      &
      &  'Carbon flux from cflux_c_green_2_herb into atmosphere (complementary to cflux_herb_2_littergreen).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2greenlitterpools', mem%cflux_dist_greenreserve_2_soil,                          &
      & hgrid, surface,                                                                                 &
      & t_cf('cflux_dist_greenreserve_2_soil', 'mol(C) m-2(canopy) s-1',                                &
      &      'Amount of C relocated by wind damage to the green litter pools.'),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=MEDIUM,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2woodlitterpools', mem%cflux_dist_woods_2_soil,                                  &
      & hgrid, surface,                                                                                 &
      & t_cf('cflux_dist_woods_2_soil', 'mol(C) m-2(canopy) s-1',                                       &
      &      'Amount of C relocated by wind damage to the woody litter pools..'),                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=MEDIUM,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_fire_2_atmos', mem%cflux_fire_all_2_atm,                                    &
      & hgrid, surface,                                                                                           &
      & t_cf('cflux_fire_all_2_atm', 'mol(C) m-2(canopy) s-1', 'Carbon flux into atmosphere.'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM,                                                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3




    ! Carbon variables on tile area
    CALL mem%Add_var('c_greenwood2litter_ta', mem%cflux_c_greenwood_2_litter_ta,                       &
      & hgrid, surface,                                                                                &
      & t_cf('c_greenwood2litter_ta', 'mol(CO2) m-2 s-1', 'Total litter flux entering the soil pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_green_ta', mem%c_green_ta,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('c_green_ta', 'mol(C) m-2', 'C-Pool for green parts of vegetation.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_woods_ta', mem%c_woods_ta,                                                   &
      & hgrid, surface,                                                                              &
      & t_cf('c_woods_ta', 'mol(C) m-2', 'C-Pool for structural material/parts of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & loutput = .TRUE., output_level=BASIC,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_reserve_ta', mem%c_reserve_ta,                                                                &
      & hgrid, surface,                                                                                               &
      & t_cf('c_reserve_ta', 'mol(C) m-2', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_crop_harvest_ta', mem%c_crop_harvest_ta,                                  &
      & hgrid, surface,                                                                           &
      & t_cf('c_crop_harvest_ta', 'mol(C) m-2', 'C-Pool for Harvested Material from Crops.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 1
    CALL mem%Add_var('c_acid_ag1_ta', mem%c_acid_ag1_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_ag1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag1_ta', mem%c_water_ag1_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_ag1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag1_ta', mem%c_ethanol_ag1_ta,                                                    &
      & hgrid, surface,                                                                                           &
      & t_cf('c_ethanol_ag1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag1_ta', mem%c_nonsoluble_ag1_ta,                                             &
      & hgrid, surface,                                                                                          &
      & t_cf('c_nonsoluble_ag1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg1_ta', mem%c_acid_bg1_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_bg1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg1_ta', mem%c_water_bg1_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_bg1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg1_ta', mem%c_ethanol_bg1_ta,                                                   &
      & hgrid, surface,                                                                                          &
      & t_cf('c_ethanol_bg1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg1_ta', mem%c_nonsoluble_bg1_ta,                                            &
      & hgrid, surface,                                                                                         &
      & t_cf('c_nonsoluble_bg1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_1_ta', mem%c_humus_1_ta,                          &
      & hgrid, surface,                                                         &
      & t_cf('c_humus_1_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
      & prefix, suffix,                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 2
    CALL mem%Add_var('c_acid_ag2_ta', mem%c_acid_ag2_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_ag2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag2_ta', mem%c_water_ag2_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_ag2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag2_ta', mem%c_ethanol_ag2_ta,                                                    &
      & hgrid, surface,                                                                                           &
      & t_cf('c_ethanol_ag2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag2_ta', mem%c_nonsoluble_ag2_ta,                                             &
      & hgrid, surface,                                                                                          &
      & t_cf('c_nonsoluble_ag2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg2_ta', mem%c_acid_bg2_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_bg2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg2_ta', mem%c_water_bg2_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_bg2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg2_ta', mem%c_ethanol_bg2_ta,                                                   &
      & hgrid, surface,                                                                                          &
      & t_cf('c_ethanol_bg2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg2_ta', mem%c_nonsoluble_bg2_ta,                                            &
      & hgrid, surface,                                                                                         &
      & t_cf('c_nonsoluble_bg2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_2_ta', mem%c_humus_2_ta,                          &
      & hgrid, surface,                                                         &
      & t_cf('c_humus_2_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
      & prefix, suffix,                                                         &
      & loutput = .TRUE., output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_sum_veg_ta', mem%C_sum_veg_ta,                         &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_veg_ta', 'mol(C) m-2', 'Sum of all vegetation C pools.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_litter_ag_ta', mem%C_sum_litter_ag_ta,                           &
      & hgrid, surface,                                                                      &
      & t_cf('C_sum_litter_ag_ta', 'mol(C) m-2', 'Sum of all above ground litter C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & output_level=BASIC,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'c_sum_litter_bg_ta', mem%C_sum_litter_bg_ta,                          &
      & hgrid, surface,                                                                      &
      & t_cf('C_sum_litter_bg_ta', 'mol(C) m-2', 'Sum of all below ground litter C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & output_level=BASIC,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_humus_ta', mem%C_sum_humus_ta,                     &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_humus_ta', 'mol(C) m-2', 'Sum of all humus C pools.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_natural_ta', mem%C_sum_natural_ta,                 &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_natural_ta', 'mol(C) m-2', 'Sum of all natural C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('soil_respiration_ta', mem%soil_respiration_ta,        &
      & hgrid, surface,                                                     &
      & t_cf('soil_respiration_ta', 'mol(C) m-2 s-1', 'Soil respiration.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & output_level=BASIC,                                                 &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_flux_corr_ta', mem%NPP_flux_correction_ta,                                 &
      & hgrid, surface,                                                                              &
      & t_cf('NPP_flux_correction_ta', 'mol(C) m-2 s-1',                                             &
      &      'Amount by which the NPP rate has to be corrected to avoid a negative reserve pool.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('excess_npp_ta', mem%excess_NPP_ta,                                                          &
      & hgrid, surface,                                                                                           &
      & t_cf('excess_NPP_ta', 'mol(C) m-2  s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    !      CALL mem%Add_var('lai_yydaymean_ta', mem%LAI_yyDayMean_ta,              &
    !        & hgrid, surface,                                                     &
    !        & t_cf('LAI_yyDayMean_ta', '-', 'Previous previous day mean.'),          &
    !        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !        & prefix, suffix,                                                     &
    !        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    !     CALL mem%Add_var('lai_ydaymean_ta', mem%LAI_yDayMean_ta,                 &
    !        & hgrid, surface,                                                     &
    !        & t_cf('LAI_yDayMean_ta', '-', 'Previous day mean.'),                    &
    !        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !        & prefix, suffix,                                                     &
    !        & output_level=BASIC,                                                                                &
    !        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_ydaymean_ta', mem%NPP_pot_yDayMean_ta,                              &
      & hgrid, surface,                                                                           &
      & t_cf('NPP_pot_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                  ! in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('npp_act_ydaymean_ta', mem%NPP_act_yDayMean_ta,                                &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean_ta', 'mol(C) m-2 s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                  ! in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean_ta', mem%GPP_yDayMean_ta,                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


    CALL mem%Add_var('root_exudates_ta', mem%root_exudates_ta,             &
      & hgrid, surface,                                                    &
      & t_cf('root_exudates_ta', 'mol(C) m-2 s-1', 'Root exudates'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


    CALL mem%Add_var('cflux_c_green_2_herb_ta', mem%cflux_c_green_2_herb_ta,                                   &
      & hgrid, surface,                                                                                    &
      & t_cf('cflux_c_green_2_herb_ta', 'mol(C) m-2 s-1', 'Total carbon flux from C pool green to herbivory.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
      & prefix, suffix,                                                                                    &
      & output_level=MEDIUM,                                                                               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2greenlitter_ta', mem%cflux_dist_greenreserve_2_soil_ta,                 &
      & hgrid, surface,                                                                         &
      & t_cf('cflux_dist_greenreserve_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to green litter.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
      & prefix, suffix,                                                                         &
      & output_level=MEDIUM,                                                                               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2woodlitter_ta', mem%cflux_dist_woods_2_soil_ta,                        &
      & hgrid, surface,                                                                        &
      & t_cf('cflux_dist_woods_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to woody litter.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=MEDIUM,                                                                               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_fire_ta', mem%co2flux_fire_all_2_atm_ta,                                    &
      & hgrid, surface,                                                                                   &
      & t_cf('co2flux_fire_all_2_atm_ta', 'kg(CO2) m-2 s-1', 'CO2 flux from C pools to atmosphere due to fire.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
      & prefix, suffix,                                                                                   &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! R: In JSBACH3 the following variables are not written into the output.
    CALL mem%Add_var('co2_l2a_npp_ta', mem%co2flux_npp_2_atm_ta,                                                     &
      & hgrid, surface,                                                                                              &
      & t_cf('co2flux_npp_2_atm_ta', 'kg(CO2) m-2 s-1', 'Net CO2 fluxes between biosphere and atmosphere due to NPP.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & output_level=BASIC,                                                                                &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_resp_ta', mem%co2flux_soilresp_2_atm_ta,                                     &
      & hgrid, surface,                                                                                    &
      & t_cf('co2flux_soilresp_2_atm_ta', 'kg(CO2) m-2 s-1',                                               &
      &      'Net CO2 fluxes between biosphere and atmosphere due to soil respiration.'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
      & prefix, suffix,                                                                                    &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_herb_ta', mem%co2flux_herb_2_atm_ta,                                          &
      & hgrid, surface,                                                                                     &
      & t_cf('co2flux_herb_2_atm_ta', 'kg(CO2) m-2 s-1',                                                    &
      &      'Net CO2 fluxes between biosphere and atmosphere due to herbivory.'),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_npp_yday_ta', mem%co2flux_npp_2_atm_yday_ta,                                           &
      & hgrid, surface,                                                                                              &
      & t_cf('co2flux_npp_2_atm_yday_ta', 'kg(CO2) m-2 s-1', 'Net daily CO2 flux, only needed for conservation test.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('ccons_calccpools', mem%cconservation_calcCpools,            &
      & hgrid, surface,                                                           &
      & t_cf('cconservation_calcCpools', 'mol(C) m-2(tile area)',                 &
      & 'C conservation test: should be zero'),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
      & prefix, suffix,                                                           &
      & output_level=BASIC,                                                       &
      & loutput=.TRUE.,lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_decomp_humus_1_sum', mem%c_decomp_humus_1_sum, hgrid, surface,                                 &
      & t_cf('c_decomp_humus_1_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon decomposed from YASSO humus (leaf)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_decomp_humus_2_sum', mem%c_decomp_humus_2_sum, hgrid, surface,                                 &
      & t_cf('c_decomp_humus_2_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon decomposed from YASSO humus (wood)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_1_sum', mem%c_into_humus_1_sum, hgrid, surface,                                     &
      & t_cf('c_into_humus_1_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon flux to humus pool in YASSO (leaf)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_2_sum', mem%c_into_humus_2_sum, hgrid, surface,                                     &
      & t_cf('c_into_humus_2_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon flux to humus pool in YASSO (wood)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)

    CALL mem%Add_var('c_decomp_humus_1_sum_ta', mem%c_decomp_humus_1_sum_ta, hgrid, surface,                      &
      & t_cf('c_decomp_humus_1_sum_ta', 'mol(C) m-2', 'annual sum of carbon decomposed from YASSO humus (leaf)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_decomp_humus_2_sum_ta', mem%c_decomp_humus_2_sum_ta, hgrid, surface,                      &
      & t_cf('c_decomp_humus_2_sum_ta', 'mol(C) m-2', 'annual sum of carbon decomposed from YASSO humus (wood)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_1_sum_ta', mem%c_into_humus_1_sum_ta, hgrid, surface,                          &
      & t_cf('c_into_humus_1_sum_ta', 'mol(C) m-2', 'annual sum of carbon flux to humus pool in YASSO (leaf)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_2_sum_ta', mem%c_into_humus_2_sum_ta, hgrid, surface,                          &
      & t_cf('c_into_humus_2_sum_ta', 'mol(C) m-2', 'annual sum of carbon flux to humus pool in YASSO (wood)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)

    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      CALL mem%Add_var('max_green', mem%current_max_green,                   &
        & hgrid, surface,                                                    &
        & t_cf('current_max_green', 'mol(C) m-2(canopy)', 'Running maximum of the green pool'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
      CALL mem%Add_var('max_green_veg_c', mem%veg_carbon_at_max_green,       &
        & hgrid, surface,                                                    &
        & t_cf('veg_carbon_at_max_green', 'mol(C) m-2(canopy)', 'Total vegetation carbon at green pool maximum'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
    ENDIF

#ifdef __ICON__
    ! Diagnostic global carbon sums for experiment monitoring 
    !       (1d stream variables are not supported with echam)
    !
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var('c_sum_veg_gsum', mem%C_sum_veg_gsum,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_veg_gsum', 'PgC', 'Global sum of all vegetation C pools.'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_litter_ag_gsum', mem%C_sum_litter_ag_gsum,                         &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_litter_ag_gsum', 'PgC', 'Global sum of all above ground litter C pools'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_litter_bg_gsum', mem%C_sum_litter_bg_gsum,                         &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_litter_bg_gsum', 'PgC', 'Global sum of all below ground litter C pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_humus_gsum', mem%C_sum_humus_gsum,                                 &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_humus_gsum', 'PgC', 'Global sum of all humus C pools.'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_natural_gsum', mem%C_sum_natural_gsum,                             &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_natural_gsum', 'PgC', 'Global sum of all natural C pools.'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('npp_act_ydaymean_gsum', mem%NPP_act_yDayMean_gsum,                       &
        & hgrid, surface,                                                                        &
        & t_cf('NPP_act_yDayMean_gsum', 'PgC yr-1', 'Global actual NPP rate of the previous day'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('gpp_ydaymean_gsum', mem%GPP_yDayMean_gsum,                               &
        & hgrid, surface,                                                                        &
        & t_cf('GPP_yDayMean_gsum', 'PgC yr-1', 'Global GPP rate of the previous day'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('soil_respiration_gsum', mem%soil_respiration_gsum,                       &
        & hgrid, surface,                                                                        &
        & t_cf('soil_respiration_gsum', 'PgC yr-1', 'Global soil respiration rate'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
    END IF
#endif

  END SUBROUTINE Init_carbon_memory

#endif
END MODULE mo_carbon_memory_class
