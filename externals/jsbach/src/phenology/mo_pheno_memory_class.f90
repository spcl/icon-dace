!> Contains the memory class for phenology.
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
MODULE mo_pheno_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of processes in this module
  dsl4jsb_Use_processes PHENO_
  ! Use process configurations
  dsl4jsb_Use_config(PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pheno_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  ! ======================================================================================================= !
  !>Type definition for pheno memory
  !>
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_pheno_memory

    TYPE(t_jsb_var_real2d) :: &
      & lai                               !< actual leaf area index [m2 m-2]
    TYPE(t_jsb_var_real3d) :: &
      & lai_cl                            !< leaf area index per canopy layer [m2 m-2]

    ! variables NOT added to varlist
    REAL(wp), POINTER :: lai_mon(:,:,:)       => NULL()   !< Monthly leaf area index (per year or climatology) <nproma x nblks x 14> !
    REAL(wp), POINTER :: fract_fpc_mon(:,:,:) => NULL()   !< Monthly fpc (per year or climatology) <nproma x nblks x 14> !

    ! Diagnostic land mean phenology e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) ::           &
      & lai_ta_gmean,                   & !< Global Land mean LAI (box tile)
      & fract_fpc_gmean                   !< Global Land mean foliage cover (box tile)

    TYPE(t_jsb_var_real2d)           :: &
      & lai_ta,                         & !< Leaf area index relative to the tile area
      & fract_fpc,                      & !< Foliage projected cover of the pft
                                          !  (fraction of ground area occupied by the vertical projection of foliage)
      & fract_fpc_max,                  & !< Maximum fraction of tile that can be covered by leaves (assuming infinite lai)
      & fract_forest,                   & !< Forest fraction
      & growth_phase_SG,                & !< For summergreen only: values are:
                                          !<   -1.0  .. during rest phase (i.e. from autumn event to spring event)
                                          !<    0.0  .. during growth phase (i.e. some weeks after spring event)
                                          !<    1.0  .. in vegetative phase (i.e. from end of growth phase to autum event)
                                          !< Updated every day by the subroutine "update_growth_phase()".
      & growth_phase_EG,                & !< For evergreen only: values are:
                                          !<    0.0  .. during growth phase (i.e. some weeks after spring event)
                                          !<    1.0  .. in vegetative phase (i.e. from ent of growth phase to spring event)
                                          !< Updated every day by the subroutine "update_growth_phase()".
      & growth_phase_CRP,               & !< For extra-tropical crops only: values are:
                                          !<     -2.0  .. during rest phase of winter crops (i.e. from spring to autumn)
                                          !<     -1.0  .. during rest phase of summer crops (i.e. from autumn to spring)
                                          !<      0.0  .. during growth phase of summer crops (i.e. from spring to autumn)
                                          !<      1.0  .. during growth phase of second crop (i.e. from summer to autumn)
                                          !<      2.0  .. during growth phase of winter crops (i.e. from autumn to spring)
                                          !< Updated every day by the subroutine "update_growth_phase()".
      & days_since_growth_begin_EG,     & !< For summergreen and evergreen only: serves to remember the number of
                                          !< days elapsed since the growth phase began
      & days_since_growth_begin_SG,     & !< (as real numbers to save it in a stream).
      & previous_day_NPP_pot_rate_ca,   & !< The mean primary production rate [mol(CO2)/(m^2 s)] at each grid point
                                          !< at previous day, because not known for current day; updated at ..
                                          !< the beginning of each new day (stream)

                                          !! JN: variables required when l_forestRegrowth = true, i.e. when the
                                          !      maxLAI is derived from biomass via allometric relationships
      & number_of_individuals,          & !< Number of individuals in last calculation [trees ha-1]
      & biomass_per_individual,         & !< Biomass per individual according to allometric relationships [kg/tree]

      & cconservation_allom,            & !< Tester for c conservation upon rescaling due to changed max lai
      & maxLAI_allom,                   & !< dynamical maximum LAI [m^2(leaf)/m^2(ground)]:
                                          !< Current max LAI according to allometric relationships --> only changes ones a year!

      ! & laiMax_dyn,                    & ! Dynamical maximum LAI [m^2(leaf)/m^2(ground)]-- not larger than the maximum  lai
      !                                    ! coming from JSBACH, but smaller than that value when NPP is negative for more
      !                                    ! than 1 year. (stream)
      ! & year_NPP_sum,                  & ! This field sums the NPP-rate during the whole year, starting on Jan,1, in the
      !                                    ! northern hemisphere and at July,1, in the southern hemisphere. A negative value
      !                                    ! indicates that summation has not yet started.
      & day_NPP_sum,                    & !< sum of NPP since midnight (needed to compute previous_day_NPP_pot_rate_ca(:,:)) [mol(CO2)/m^2]
      & heat_sum_EG,                    & !< Heat sums --- needed for the phenology of summergreens, evergreens and crops.
      & heat_sum_SG,                    & !< Where heat summation has not yet started, "heat_sum" is set to a value < -1.
      & heat_sum_CRP,                   & !< For summer crops heat summation starts with "0".
      & heat_sum_winter,                & !< For winter crops heat summation starts with EPSILON(1.) and is positive in the
                                          !< growing period (autumn to spring) and negative during summer. (stream)
      & chill_days_SG,                  & !< Number of chill days --- needed for the phenology of summer- and
      & chill_days_EG,                  & !< evergreens (stream)
      & veg_fract_correction              !< Factor representing clumpyness of vegetation (for each tile type).
                                          !< To write it in the output file we define it as 2d even when there is no change
                                          !< between gridcells as it is in fact a 1d constant for each tile.
                                          !   JN: with l_forestRegrowth = true veg_fract_correction can change over time
                                          !       and is spatially variable
                                          !       if l_forestRegrowth = false there is still no change between gridcells
                                          !       and it is a 1d constant for each tile in this case.
                                          ! Gives the fraction of canopy area from the tile area. Between 0 and 1.
                                          ! canopy area = veg_fract_correction * tile area.
                                          ! var (tile area) = veg_fract_correction * var(canopy area).
                                          ! ACHTUNG: da wir vorerst fract_fpc_max zur Abgrenzung zwischen "Wueste" und Vegetation
                                          !          in den JS4 Berechnungen haben, gilt eigentlich dies hier:
                                          ! var (tile area) = veg_fract_correction * var(canopy area)  * fract_fpc_max.

  CONTAINS
    PROCEDURE :: Init => Init_pheno_memory
  END TYPE t_pheno_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: phenology
  !>
  SUBROUTINE Init_pheno_memory(mem, prefix, suffix, lct_ids, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_pheno_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(PHENO_)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),   POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface                      ! Vertical grid
    INTEGER                     :: table                        ! ...
    TYPE(t_grib2)               :: grib2_desc                   ! used for lai settings
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pheno_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model        => Get_model(model_id)
    table        = tables(1)
    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')

    dsl4jsb_Get_config(PHENO_)

    IF ( (     One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &  ) ) THEN

      ! IF (TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'climatology') THEN   ! TODO
        ALLOCATE(mem%lai_mon      (hgrid%nproma, hgrid%nblks, 0:13))
        ALLOCATE(mem%fract_fpc_mon(hgrid%nproma, hgrid%nblks, 0:13))
        !$ACC ENTER DATA CREATE(mem%lai_mon, mem%fract_fpc_mon)
        !$ACC ENTER DATA ATTACH(mem%lai_mon, mem%fract_fpc_mon)
      ! END IF

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(2,0,28, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'lai', mem%lai,                                        &
        & hgrid, surface,                                                      &
        & t_cf('leaf_area_index', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), grib2_desc,                          &
        & prefix, suffix,                                                      &
        & loutput=.TRUE., output_level=BASIC,                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lai_ta', mem%lai_ta,                                  &
        & hgrid, surface,                                                      &
        & t_cf('leaf_area_index', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (TRIM(mem%owner_tile_name) == 'veg') THEN
        grib2_desc = t_grib2(2,0,228, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'fract_fpc', mem%fract_fpc,                            &
        & hgrid, surface,                                                      &
        & t_cf('foliage_projected_cover', '', ''),                             &
        & t_grib1(table, 255, grib_bits), grib2_desc,                          &
        & prefix, suffix,                                                      &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'fract_fpc_max', mem%fract_fpc_max,                    &
        & hgrid, surface,                                                      &
        & t_cf('maximum_fpc', '', ''),                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'fract_forest', mem%fract_forest,                      &
        & hgrid, surface,                                                      &
        & t_cf('forest_fraction', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'veg_fract_correction', mem%veg_fract_correction,      &
        & hgrid, surface,                                                      &
        & t_cf('Correction factor to convert tile to canopy area', '', '') ,   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=FULL,                                                   &
        & lrestart=.TRUE.,                                                     &
        & initval_r=0.0_wp )

      IF (TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'logrop') THEN

        CALL mem%Add_var('growth_phase_sg', mem%growth_phase_SG,               &
          & hgrid, surface,                                                    &
          & t_cf('growth_phase_SummerGreen', '', ''),                          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & initval_r=-1.0_wp) ! Start in rest phase for summergreens

        CALL mem%Add_var('growth_phase_eg', mem%growth_phase_EG,               &
          & hgrid, surface,                                                    &
          & t_cf('growth_phase_EverGreen', '', ''),                            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & initval_r=1.0_wp ) ! Start in vegetative phase for evergreens

        CALL mem%Add_var('growth_phase_crp', mem%growth_phase_CRP,             &
          & hgrid, surface,                                                    &
          & t_cf('growth_phase_Crop', '', ''),                                 &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & initval_r=-1.0_wp )

        CALL mem%Add_var('growth_days_eg', mem%days_since_growth_begin_EG,      &
          & hgrid, surface,                                                     &
          & t_cf('days_since_growth_begin_EverGreen', '', ''),                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
          & prefix, suffix,                                                     &
          & initval_r=0.0_wp )

        CALL mem%Add_var('growth_days_sg', mem%days_since_growth_begin_SG,      &
          & hgrid, surface,                                                     &
          & t_cf('days_since_growth_begin_SummerGreen', '', ''),                &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
          & prefix, suffix,                                                     &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=0.0_wp )

        CALL mem%Add_var('yday_npp_pot_rate', mem%previous_day_NPP_pot_rate_ca,      &
          & hgrid, surface,                                                          &
          & t_cf('yday_NPP_pot_rate_ca', 'prev day NPP_pot_rate_ca', ''),            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
          & prefix, suffix,                                                          &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=0._wp )

        ! CALL mem%Add_var('laimax_dyn', mem%laiMax_dyn,                         &
        !   & hgrid, surface,                                                    &
        !   & t_cf('laiMax_dyn', '', ''),                                        &
        !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        !   & prefix, suffix,                                                    &
        !   & loutput=.FALSE., lrestart=.TRUE., initval_r=-1.0_wp)

        ! CALL mem%Add_var('year_npp_sum', mem%year_NPP_sum,                      &
        !   & hgrid, surface,                                                     &
        !   & t_cf('year_NPP_sum', '', ''),                                       &
        !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        !   & prefix, suffix,                                                     &
        !   & loutput=.FALSE., lrestart=.TRUE., initval_r=2.0_wp)

        CALL mem%Add_var('day_npp_sum', mem%day_NPP_sum,                       &
          & hgrid, surface,                                                    &
          & t_cf('day_NPP_sum', '', ''),                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.TRUE., lrestart=.TRUE., initval_r=0._wp)

        CALL mem%Add_var('heat_sum_eg', mem%heat_sum_EG,                       &
          & hgrid, surface,                                                    &
          & t_cf('heat_sum_EG', '', ''),                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=-99.0_wp)

        CALL mem%Add_var('heat_sum_sg', mem%heat_sum_SG,                       &
          & hgrid, surface,                                                    &
          & t_cf('heat_sum_SG', '', ''),                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=-99.0_wp) ! has to be lower than -1

        CALL mem%Add_var('heat_sum_crp', mem%heat_sum_CRP,                     &
          & hgrid, surface,                                                    &
          & t_cf('heat_sum_CRP', '', ''),                                      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=9999.0_wp)

        CALL mem%Add_var('heat_sum_winter', mem%heat_sum_winter,               &
          & hgrid, surface,                                                    &
          & t_cf('heat_sum_winter', '', ''),                                   &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=-EPSILON(1._wp))

        CALL mem%Add_var('chill_days_eg', mem%chill_days_EG,                   &
          & hgrid, surface,                                                    &
          & t_cf('chill_days_EG', '', ''),                                     &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=0._wp)

        CALL mem%Add_var('chill_days_sg', mem%chill_days_SG,                   &
          & hgrid, surface,                                                    &
          & t_cf('chill_days_SG', '', ''),                                     &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.FALSE., lrestart=.TRUE., initval_r=0._wp)

        IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN

          CALL mem%Add_var('maxLAI_allom', mem%maxLAI_allom,                     &
            & hgrid, surface,                                                    &
            & t_cf('maxLAI_allom', 'm^2(leaf)/m^2(ground)', 'maximum LAI according to allometric relationships'),  &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & loutput=.TRUE., lrestart=.TRUE., initval_r=0.4_wp)

          !@todo: variables are only relevant for forest pfts - how to constrain memory allocation to forest pfts?
          CALL mem%Add_var('individuals', mem%number_of_individuals,             &
            & hgrid, surface,                                                    &
            & t_cf('number_of_individuals', 'trees ha-1', 'Number of individuals in last calculation'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & loutput=.TRUE., lrestart=.TRUE., initval_r=15000.0_wp)

          CALL mem%Add_var('biomass_per_ind', mem%biomass_per_individual,        &
            & hgrid, surface,                                                    &
            & t_cf('biomass_per_individual', 'kg/tree', 'Biomass per individual according to allometric relationships'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)

          CALL mem%Add_var( 'cconservation_allom', mem%cconservation_allom,             &
            & hgrid, surface,                                                           &
            & t_cf('cconservation_allom', 'mol(C) m-2(canopy)', &
              & 'Test for c conservation on rescaling of c states for allom.'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
            & prefix, suffix,                                                           &
            & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        END IF ! l_forestRegrowth
      END IF ! logrop

#ifdef __ICON__
      ! Diagnostic global carbon sums for experiment monitoring
      !   (1d stream variables are not supported with echam)
      IF (TRIM(suffix) == 'box') THEN

        CALL mem%Add_var( 'lai_ta_gmean', mem%lai_ta_gmean,                           &
          & hgrid, surface,                                                           &
          & t_cf('lai_ta_gmean', '1', 'Global land mean leaf_area_index'),            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
          & prefix, suffix,                                                           &
          & lrestart=.FALSE., initval_r=0.0_wp )

        CALL mem%Add_var( 'fract_fpc_gmean', mem%fract_fpc_gmean,                     &
          & hgrid, surface,                                                           &
          & t_cf('fract_fpc_gmean', '1', 'Global land mean foilage projected cover'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
          & prefix, suffix,                                                           &
          & lrestart=.FALSE., initval_r=0.0_wp )

      END IF
#endif
    END IF !  VEG / LAND

  END SUBROUTINE Init_pheno_memory

#endif
END MODULE mo_pheno_memory_class
