!> Methods for JSBACH initializations
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
!>#### Initialize JSBACH instance with grid, task queue, tile structure, initial/boundary conditions etc.
!>
MODULE mo_jsb_model_init

#ifndef __NO_JSBACH__

  dsl4jsb_Use_processes SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, &
    &                   FUEL_, PPLCC_, FAGE_, ALCC_, FLCC_, WLCC_, NLCC_, A2L_, VEG_, SB_, SPQ_,  &
    &                   Q_PHENO_, Q_ASSIMI_, Q_RAD_

  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(FAGE_)
  dsl4jsb_Use_config(NLCC_)
  dsl4jsb_Use_config(HYDRO_)
#ifndef __NO_QUINCY__
  dsl4jsb_Use_config(VEG_)
#endif

  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_control,      ONLY: get_no_of_models, debug_on
  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_grid,         ONLY: Get_grid
  USE mo_jsb_grid_class,   ONLY: t_jsb_grid

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message, finish
  USE mo_jsb_parallel,     ONLY: Get_omp_thread, Is_omp_inside_serial
  USE mo_jsb_time,         ONLY: is_time_restart, get_time_start

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_setup_grid, jsbach_init, jsbach_init_after_restart

  !> Interface to set up grid of a model instance
  !
  INTERFACE jsbach_setup_grid
    MODULE PROCEDURE setup_grid_standalone
    MODULE PROCEDURE setup_grid_patch
  END INTERFACE jsbach_setup_grid

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'


CONTAINS


  ! ====================================================================================================== !
  !
  !> Set up grid for standalone model if using ECHAM infrastructure
  !
  SUBROUTINE setup_grid_standalone(nproca, nprocb, npedim, model_id, model_name)

#ifndef __ICON__
    ! >> TODO: echam
    USE mo_gaussgrid,            ONLY : inigau
    USE mo_echam_convect_tables, ONLY : init_convect_tables
    USE mo_netCDF,               ONLY: add_dim
    ! <<
    USE mo_jsb_grid_iface,       ONLY: get_lons_and_lats_from_file, get_nlon_nlat_from_file
    USE mo_jsb_domain_iface,     ONLY: t_patch, jsbalone_init_decomp
    USE mo_jsb_grid_class,       ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,             ONLY: Register_grid
#endif

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in)                    :: nproca !< Number of processors in 1st dimension
                                                     !< of 2-d decomposition
    INTEGER, INTENT(in)                    :: nprocb
        !< Number of processors in 2nd dimension of 2-d decomposition
    INTEGER, INTENT(in)                    :: npedim
        !< Size of blocks on each domain
    INTEGER, INTENT(in)                    :: model_id
        !< ID of model instance
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: model_name
        !< Name of model instance
    ! -------------------------------------------------------------------------------------------------- !

#ifndef __ICON__
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid
    TYPE(t_patch)               :: decomposition

    INTEGER :: nlon, nlat
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':setup_grid_standalone'
    ! ---------------------------------------------------------------------------------------------------!

#ifdef __ICON__
    ! avoid compiler warnings about dummy arguments not being used
    IF (.FALSE.) WRITE(*,*) nproca, nprocb, npedim, model_id, model_name
    CALL finish(TRIM(routine),'This subroutine should not be called with ICON')
#else
    model => Get_model(model_id)

    IF (PRESENT(model_name)) THEN
      IF (ADJUSTL(TRIM(model_name)) /= ADJUSTL(TRIM(model%name))) &
        CALL finish(TRIM(routine), 'Model ID and name are inconsistent with namelist config.')
    END IF

    CALL get_nlon_nlat_from_file(model%config%grid_filename, nlon, nlat)
    CALL jsbalone_init_decomp(nlon, nlat, nproca, nprocb, npedim, decomposition)

    CALL inigau               ! beside others required for get_area, i.e. required before new_grid
    CALL init_convect_tables

    CALL add_dim("lon", nlon, "longitude", "degrees_east" ) !required to write restart files
    CALL add_dim("lat", nlat, "latitude", "degrees_north")

    CALL get_lons_and_lats_from_file(model%config%grid_filename, decomposition%nproma, decomposition%ngpblks)
        ! needs to be done after the decomposition is initialised; required for get_lon and get_lat,
        ! especially for site-level sims
    grid => new_grid(decomposition, type='echam')

    CALL register_grid(grid)
    model%grid_id = grid%id

    ! CALL init_model_common(model)

    NULLIFY(model)
#endif

  END SUBROUTINE setup_grid_standalone


  ! ====================================================================================================== !
  !
  !> Set up grid provided by host model (ICON coupled/standalone or ECHAM coupled)
  !
  SUBROUTINE setup_grid_patch( model_id, patch, type)

    USE mo_jsb_domain_iface, ONLY: t_patch
    USE mo_jsb_grid_class,   ONLY: new_grid
    USE mo_jsb_grid,         ONLY: Register_grid

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in)               :: model_id !< ID of model instance
    TYPE(t_patch), TARGET, INTENT(in) :: patch    !< Grid information provided by host model
    CHARACTER(LEN=*), INTENT(in)      :: type     !< Type of host model ('icon' or 'echam')
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':setup_grid_patch'
    ! -------------------------------------------------------------------------------------------------- !

    IF (model_id > get_no_of_models()) &
      & CALL finish(routine, 'Host request to more jsbach models than are defined')

    model => Get_model(model_id)

    grid => new_grid(patch, type)
    CALL Register_grid(grid)
    model%grid_id = grid%id

    ! CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE setup_grid_patch


  ! ====================================================================================================== !
  !>
  !> Initializes an ICON-Land model instance
  !>
  !>   model_scheme: MODEL_JSBACH, MODEL_QUINCY
  !>
  SUBROUTINE jsbach_init(model_id)

    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_, AGGREGATE_
    USE mo_jsb_model_class,     ONLY: print_model, MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, &
      &                               t_jsb_aggregator, t_jsb_aggregator_weighted_by_fract
    USE mo_jsb_lct_class,       ONLY: Contains_lct
    USE mo_jsb_io_netcdf,       ONLY: jsb_netcdf_open_input
    USE mo_seb_init,            ONLY: seb_init
    USE mo_turb_init,           ONLY: turb_init
    USE mo_sse_init,            ONLY: sse_init, init_soil_properties
    USE mo_hydro_init,          ONLY: hydro_init
    USE mo_pheno_init,          ONLY: pheno_init
    USE mo_rad_init,            ONLY: rad_init
#ifndef __NO_QUINCY__
    USE mo_quincy_model_config, ONLY: QLAND, QPLANT, QCANOPY, QSOIL
    USE mo_q_rad_init,          ONLY: q_rad_init
    USE mo_q_assimi_init,       ONLY: q_assimi_init
    USE mo_q_pheno_init,        ONLY: q_pheno_init
    USE mo_spq_init,            ONLY: spq_init
    USE mo_veg_init,            ONLY: veg_init
    USE mo_sb_init,             ONLY: sb_init
#endif
    USE mo_fage_init,           ONLY: fage_init

    USE mo_fage_init_lcc,       ONLY: fage_init_lcc
    USE mo_alcc_init_lcc,       ONLY: alcc_init_lcc
    USE mo_flcc_init_lcc,       ONLY: flcc_init_lcc
    USE mo_wlcc_init_lcc,       ONLY: wlcc_init_lcc

#ifndef __NO_JSBACH_HD__
    USE mo_hd_init,             ONLY: hd_init_ic, hd_init_bc
#endif
    USE mo_jsb_io,              ONLY: cdiDefMissval, missval
    USE mo_util,                ONLY: one_of

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: model_id !< ID of the model instance to be initialized
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_grid),  POINTER          :: grid
    CLASS(t_jsb_tile_abstract), POINTER :: tile, child
    INTEGER                             :: ilct, jlct, i, no_children, lct_start, blks, blke, blk, jcs, jce
    REAL(wp), ALLOCATABLE               :: fract_tmp(:,:), tile_lct_fract_tmp(:,:)
    INTEGER                             :: no_omp_thread
    LOGICAL                             :: l_fixed_fractions
    CLASS(t_jsb_aggregator), POINTER    :: weighted_by_fract

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_config(FAGE_)
#ifndef __NO_QUINCY__
    dsl4jsb_Def_config(VEG_)
#endif

    CHARACTER(len=*), PARAMETER         :: routine = modname//':jsbach_init'
    ! -------------------------------------------------------------------------------------------------- !

    no_omp_thread = Get_omp_thread()

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    CALL cdiDefMissval(missval)

    model => Get_model(model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_config(FAGE_)

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 1. Populate task queue
    !>
    !>     @todo Should be in separate subroutine
    !>
    !>     Three Queues: MODEL_JSBACH / MODEL_QUINCY
    !>
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      ! Note: JN changed the order of the processes
      ! - PHENO processes were moved up: need to happen before C_NPP_pot_allocation for allometric
      !   considerations
      ! - processes potentially changing the area of a tile are conducted first in the timestep.
      ! - the task npp_buffer is performed even before to provide already updated multi annual average
      !   NPP to the NLCC process at the beginning of a year

      IF (model%processes(NLCC_)%p%config%active) THEN
        CALL model%Queue_task(ASSIMI_, "npp_buffer")
      END IF

      ! - area damaged by natural disturbances in this timestep are calculated before any lcc takes place
      IF (model%processes(DISTURB_)%p%config%active) THEN
        CALL model%Queue_task(DISTURB_, "damaged_area")
        CALL model%Queue_task(DISTURB_, "natural_disturbances")
      END IF

      ! ------------------------ tasks that can change area
      !
      ! Note: PPLCC_ is automatically actived in mo_jsb_base.f90:jsbach_setup_models
      ! if any process that changes area (e.g. NLCC_, ALCC_) is active.
      !
      IF (model%processes(PPLCC_)%p%config%active) THEN
        ! Prepare LCC tasks
        CALL model%Queue_task(PPLCC_, "pplcc_pre_lcc")

        !***** Add LCC processes between pplcc_pre and pplcc_post *****
        ! order of tasks of course determines area available for change for the subsequent lcc tasks!
        ! suggestion: natural changes before anthropogenic changes; wind before fire (prev. order in disturb)

        IF (model%processes(FAGE_)%p%config%active) THEN
          CALL model%Queue_task(FAGE_, "fage")
        END IF

        IF (model%processes(WLCC_)%p%config%active) THEN
          CALL model%Queue_task(WLCC_, "wlcc")
        END IF

        IF (model%processes(FLCC_)%p%config%active) THEN
          CALL model%Queue_task(FLCC_, "flcc")
        END IF

        IF (model%processes(NLCC_)%p%config%active) THEN
          CALL model%Queue_task(NLCC_, "nlcc")
        END IF

        IF (model%processes(ALCC_)%p%config%active) THEN
          CALL model%Queue_task(ALCC_, "alcc")
        END IF

        ! Do postprocessing after all LCC tasks were executed
        CALL model%Queue_task(PPLCC_, "pplcc_post_lcc")
      END IF

      ! FLCC (which runs on the veg tile) changes the co2flux_fire_all_2_atm_ta variable on the leaves
      ! ... it therefore requires a carbon aggregation task on the leaves to run subsequent to the flcc task
      IF (model%processes(FLCC_)%p%config%active .AND. model%processes(CARBON_)%p%config%active) THEN
        CALL model%Queue_task(CARBON_, "ta_vars_after_lcc")
      END IF

      ! ------------------------ tasks that do not change area
      !
      CALL model%Queue_task(RAD_,    "surface_radiation")
      CALL model%Queue_task(PHENO_,  "phenology")
        ! JN moved phenology up - for allom it needs to happen before C_NPP_pot_allocation and
        ! with area changes it should happen before assimi
        ! --- before it was between humidity_scaling and foliage_projected_cover
        ! + prev comment: can be placed between the task "surface_radiation" and "asselin_filter"
      CALL model%Queue_task(PHENO_,  "foliage_projected_cover")
      IF (model%config%use_tmx) THEN
        CALL model%Queue_task(TURB_, "surface_roughness")
        CALL model%Queue_task(TURB_, "exchange_coefficients")
      END IF
      ! The first call to task `snow_and_wet_fraction` currently needs to be before `surface_hydrology`,
      ! `soil_hydrology` and `update_albedo` because of pond scheme
      ! TODO: add cross-check
      CALL model%Queue_task(HYDRO_,  "snow_and_wet_fraction")
      CALL model%Queue_task(HYDRO_,  "soil_properties")
      CALL model%Queue_task(SEB_,    "surface_energy")
      CALL model%Queue_task(HYDRO_,  "evaporation")
      CALL model%Queue_task(SEB_,    "surface_fluxes")
      IF (model%processes(ASSIMI_)%p%config%active) THEN
        CALL model%Queue_task(RAD_,    "radiation_par")
      END IF
      CALL model%Queue_task(HYDRO_,  "surface_hydrology")
      CALL model%Queue_task(HYDRO_,  "snow_and_wet_fraction")
      CALL model%Queue_task(SSE_,    "soil_and_snow_properties")
      CALL model%Queue_task(SSE_,    "soil_and_snow_temperature")
      CALL model%Queue_task(HYDRO_,  "soil_hydrology")
#ifndef __NO_JSBACH_HD__
      CALL model%Queue_task(HD_,     "hd_local_outflow")
#endif
      CALL model%Queue_task(SEB_,    "snowmelt_correction")
      IF (model%processes(ASSIMI_)%p%config%active) THEN
        ! Normally called before update_canopy_cond_unstressed and update_assimilation
        CALL model%Queue_task(ASSIMI_, "assimilation_scaling_factors")
      END IF
      IF (model%processes(ASSIMI_)%p%config%active) THEN
        CALL model%Queue_task(ASSIMI_, "canopy_cond_unstressed")
      ELSE
        CALL model%Queue_task(HYDRO_,  "canopy_cond_unstressed")
      END IF
      CALL model%Queue_task(HYDRO_, "water_stress")
      IF (model%processes(ASSIMI_)%p%config%active) THEN
        CALL model%Queue_task(ASSIMI_, "canopy_cond_stressed")
        CALL model%Queue_task(ASSIMI_, "assimilation")
      ELSE
        CALL model%Queue_task(HYDRO_,  "canopy_cond_stressed")
      END IF
      ! CALL model%Queue_task(HYDRO_,  "snow_and_wet_fraction")
      IF (model%processes(CARBON_)%p%config%active) THEN
         CALL model%Queue_task(CARBON_, "C_NPP_pot_allocation")
      END IF
      CALL model%Queue_task(TURB_,   "humidity_scaling")
      IF (.NOT. model%config%use_tmx) THEN
        CALL model%Queue_task(TURB_,   "surface_roughness")
      END IF
      CALL model%Queue_task(RAD_,    "albedo")
          ! Albedo is updated at the end of the time step, one time step before the radiation is updated in
          ! the atmosphere. Otherwise, in between, the albedo changes would be consindered in the jsbach
          ! radiation but where not considered in the atmospheric radiation calculations. This would be
          ! inconsistent. Furthermore, atmospheric cloud changes affect radiation also only at radiation
          ! time step and jsbach should be also consistent with this.
          ! In the albedo update time step the albedo should be updated after the SEB process, so that also
          ! the SEB remains consistent within this time step.
      CALL model%Queue_task(HYDRO_,  "water_balance")
      CALL model%Queue_task(SEB_,    "asselin_filter")      ! Should always be the last task

      ! Fuel is calculated last in a timestep to be right available in the next timestep
      IF (model%processes(FUEL_)%p%config%active) THEN
        CALL model%Queue_task(FUEL_, "fuel")
      END IF

#ifndef __NO_QUINCY__
    ! ----------------------------------------------------------------------------------------------------- !
    !
    ! quincy - using the SPQ_ process for soil physics
    !
    ! QUINCY model modes represent different levels of model complexity:
    !
    ! By setting the 'quincy_model' in the namelist, 'qmodel_id' is set to QCANOPY, QPLANT, QLAND, QSOIL
    ! QCANOPY, QPLANT, QLAND share many tasks (indicated with 'ALL' below)
    !
    ! QCANOPY: CANOPY model mode - basic set of tasks ('ALL') excluding the ones specific to other model modes
    !                              (calling only tasks to calculate biogeophysical processes)
    ! QPLANT:  PLANT  model mode - with additional tasks calculating vegetation biogeochemistry and dynamics
    ! QLAND:   LAND   model mode - adding tasks to calculate soil biogeochemistry
    !                              (on top of the tasks specific for the PLANT model mode)
    !
    ! Later a separated task queue set-up might be specified for QSOIL
    ! # QSOIL:   SOIL   model mode - standalone soil biogeochemistry with forcing from vegetation
    !
    CASE (MODEL_QUINCY)
      dsl4jsb_Get_config(VEG_)

      !-- currently only QLAND, QPLANT and QCANOPY are ported from Quincy standalone (QS)
      IF(.NOT. (model%config%qmodel_id == QLAND .OR. model%config%qmodel_id == QPLANT &
                .OR. (model%config%qmodel_id == QCANOPY))) THEN
        CALL finish(TRIM(routine), &
          & 'Unexpected / not yet implemented Quincy model mode: '//TRIM(model%config%quincy_model))
      ENDIF

      !-- set flux variables to zero
      CALL model%Queue_task(Q_ASSIMI_, "reset_q_assimi_fluxes")   ! ALL
      CALL model%Queue_task(SPQ_, "reset_spq_fluxes")             ! ALL
      CALL model%Queue_task(VEG_, "reset_veg_fluxes")             ! ALL
      IF(model%config%qmodel_id == QLAND) THEN
        CALL model%Queue_task(SB_, "reset_sb_fluxes")             ! LAND ONLY
      ENDIF

      !-- science
      IF(model%config%qmodel_id == QLAND) THEN
        CALL model%Queue_task(SB_, "sb_rate_modifier")            ! LAND ONLY
      ENDIF
      CALL model%Queue_task(SPQ_, "soil_turbulence")              ! ALL
      CALL model%Queue_task(Q_RAD_, "q_radiation")                ! ALL
      CALL model%Queue_task(Q_PHENO_, "q_phenology")              ! ALL
      CALL model%Queue_task(Q_ASSIMI_, "canopy_fluxes")           ! ALL
      CALL model%Queue_task(SPQ_, "spq_physics")                  ! ALL

      IF(model%config%qmodel_id == QLAND .OR. model%config%qmodel_id == QPLANT) THEN
        IF(dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
          CALL model%Queue_task(VEG_, "products_decay")           ! LAND or PLANT
        ENDIF
        CALL model%Queue_task(VEG_, "veg_turnover")               ! LAND or PLANT
        CALL model%Queue_task(VEG_, "veg_dynamics")               ! LAND or PLANT
        IF(model%config%qmodel_id == QLAND) THEN
          CALL model%Queue_task(SB_, "sb_simple_model")           ! LAND ONLY
        ENDIF
        CALL model%Queue_task(VEG_, "veg_growth")                 ! LAND or PLANT
      ENDIF

      CALL model%Queue_task(VEG_, "veg_pools")                    ! ALL
      IF(model%config%qmodel_id == QLAND) THEN
        CALL model%Queue_task(SB_, "sb_pools")                    ! LAND ONLY
      ENDIF

      !-- averaging variables in time
      CALL model%Queue_task(Q_RAD_, "tavrg_q_radiation")          ! ALL
      CALL model%Queue_task(Q_ASSIMI_, "tavrg_assimilation")      ! ALL
      CALL model%Queue_task(VEG_, "tavrg_vegetation")             ! ALL
      IF(model%config%qmodel_id == QLAND .OR. model%config%qmodel_id == QPLANT) THEN
        CALL model%Queue_task(SB_, "tavrg_soilbiogeochemistry")   ! LAND or PLANT
      ENDIF
#endif

    END SELECT
    ! Task Queue created
    !--------------------------------------------------------------------------------------------------- !

    ! Debug message after creating task queue
    IF (debug_on()) THEN
      CALL print_model(model)
      CALL model%Print_actions()
    END IF

    ! At this point, the model configuration and tile structure is set up.

    ! We continue with initializing the memory structure and different model processes.

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 2. Iterate over tiles and call process initialization if necessary
    !>
    l_fixed_fractions = .NOT. model%Do_fractions_change()
    IF (l_fixed_fractions) THEN
      CALL message(routine, 'Tile fractions are fixed over time')
    ELSE
      CALL message(routine, 'Tile fractions vary over time')
    END IF

    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)

    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      ! ---------------------------------------------------------------------------------------------- !
      !>
      !>     2.1 Initialize each tile
      !>
      ! Initialize tile process memories and read tile fractions
      CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name),                &
        & grid_id = model%grid_id,                                              &
        & in_var_groups=one_of(TRIM(tile%name), model%config%output_tiles) > 0)

#ifndef __NO_QUINCY__
      ! setup bgc material store (if running with processes that have bgcms)
      CALL tile%Count_and_classify_bgc_materials()
      CALL tile%Collect_bgc_materials()
#endif

      ! ---------------------------------------------------------------------------------------------- !
      !>
      !>     2.2 Initialize fractions of each tile
      !>
      CALL tile%Init_fractions(                                           &
        & TRIM(model%shortname), '', TRIM(tile%name),                     &
        & l_fixed_fractions=l_fixed_fractions,                            &
        & l_rel_fractions=model%config%relative_fractions_in_file)

      ! ---------------------------------------------------------------------------------------------- !
      !>
      !>     2.3 Initialize processes on each tile
      !>
      IF (model%config%init_from_ifs) &
        & model%config%ifs_input_file = jsb_netcdf_open_input(model%config%ifs_filename, model%grid_id)

      ! TBD: Loop over processes instead of explicit

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(RAD_)) THEN
          CALL rad_init(tile)
        END IF
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        IF (tile%Is_process_active(Q_RAD_)) THEN
          CALL q_rad_init(tile)
        END IF
#endif
      END SELECT

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(SEB_)) THEN
          CALL seb_init(tile)
        END IF
        IF (tile%Is_process_active(TURB_)) THEN
          CALL turb_init(tile)
        END IF
      END SELECT

      ! Initialize PHENO_ / Q_PHENO_ before HYDRO_ and VEG_
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(PHENO_)) THEN
          CALL pheno_init(tile)
        END IF
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        IF (tile%Is_process_active(Q_PHENO_)) THEN
          CALL q_pheno_init(tile)
        END IF
#endif
      END SELECT

      ! HYDRO_ & SSE_
      ! SPQ_, SB_
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(SSE_)) THEN
          CALL sse_init(tile)
        END IF
        IF (tile%Is_process_active(HYDRO_)) THEN
          CALL hydro_init(tile)
        END IF
      ! Needs both HYDRO_ and SSE_ to continue soil initialization from IFS
        IF (tile%Is_process_active(SSE_)) THEN
          CALL init_soil_properties(tile)
        END IF

#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        IF (tile%Is_process_active(SPQ_)) THEN
          CALL spq_init(tile)
        END IF
        IF (tile%Is_process_active(VEG_)) THEN
          CALL veg_init(tile)
        END IF
        IF(model%config%qmodel_id == QLAND .OR. model%config%qmodel_id == QPLANT .OR. model%config%qmodel_id == QSOIL) THEN
          IF (tile%Is_process_active(SB_)) THEN
            CALL sb_init(tile)
          END IF
        ENDIF
#endif
      END SELECT

      SELECT CASE (model%config%model_scheme)
      ! CASE (MODEL_JSBACH)
        !IF (tile%Is_process_active(ASSIMI_)) THEN
        !  CALL assimi_init(tile)
        !END IF
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        IF (tile%Is_process_active(Q_ASSIMI_)) THEN
          CALL q_assimi_init(tile)
        END IF
#endif
      END SELECT

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(FAGE_)) THEN
          CALL fage_init(tile)
        END IF
      END SELECT

#ifndef __NO_JSBACH_HD__
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(HD_)) THEN
          CALL hd_init_bc(tile)
          CALL hd_init_ic(tile)
        END IF
      END SELECT
#endif

      CALL model%Goto_next_tile(tile)
      ! ---------------------------------------------------------------------------------------------- !
    END DO

    IF (model%config%ifs_input_file%is_open) CALL model%config%ifs_input_file%Close()

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 3. Iterate over tiles again and
    !>
    !>     3.1 set initial tile fractions for aggregator(s)
    !>
    !>     3.2 set LCT fractions of the specific tile
    !>
    ! It is important to only do this after all tile fractions have been set above
    !
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)

    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    ALLOCATE(fract_tmp(SIZE(tile%lcts(1)%fract,1), SIZE(tile%lcts(1)%fract,2)))
    ALLOCATE(tile_lct_fract_tmp(SIZE(tile%lcts(1)%fract,1), SIZE(tile%lcts(1)%fract,2)))

    ! Loop over the tiles, starting from the top tile (i.e. box)
    DO WHILE (ASSOCIATED(tile))

      ! First only handle the leaves: skip the tile in case it has children.
      IF (.NOT. tile%visited(no_omp_thread) .AND. tile%Has_children()) THEN
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      ! Set initial tile fractions for 'weighted_by_fract' aggregator
      IF (tile%Has_children()) THEN
        ! Tile has children, the tile thus had been visited before and we're on the way up
        weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
        SELECT TYPE (weighted_by_fract)
        TYPE IS (t_jsb_aggregator_weighted_by_fract)
          CALL weighted_by_fract%Set_fractions(tile)
        END SELECT
      END IF

      ! Set (relative) land cover type fractions of the current tile

      ! The rel. fraction of the primary LCT is 1: The tile is fully covered by this LCT.
      ! Todo: Are these arrays (i.e. fract_glac_glac, fract_veg_veg, ...) of any interest, or
      ! should we remove the definition here and also remove the corresponding Add_var calls?
      tile%lcts(1)%fract(:,:) = 1._wp

      ! We start with second LCT, as the tile's primary LCT fraction has just been set.
      lct_start=2
      ! There is however an exception: The primary LCT of the box tile is the LCT of the first child.
      IF (.NOT. ASSOCIATED(tile%parent_tile)) lct_start=1
      DO ilct=lct_start,SIZE(tile%lcts)

        ! Leaf tiles should have only the primary LCT and should not enter this loop.
        IF (.NOT. tile%Has_children()) CALL finish(TRIM(routine), 'Leave tiles should only have one LCT')

        ! Tile has children. The tile thus had been visited before (s.above) and we're on the way up.
        ! Loop over the children and add their contributions to the specific LCT fraction.
        tile_lct_fract_tmp(:,:) = 0._wp
        no_children = tile%Get_no_of_children()
        child => tile%Get_first_child_tile()
        DO i=1,no_children
          IF (Contains_lct(child%lcts, tile%lcts(ilct)%id)) THEN
            ! The child contains the specific land cover type
            DO jlct=1,SIZE(child%lcts)
              IF (child%lcts(jlct)%id == tile%lcts(ilct)%id) THEN
                ! Get the (absolute) child tile fraction
                CALL child%Get_fraction(fract=fract_tmp(:,:))
                ! In case the child has only one LTC, LCT fraction and tile fraction are
                ! identical and we are done. If it has several LCTs, the fraction of the
                ! current LCT still needs to be calculated:
                IF (jlct > 1) THEN
                  fract_tmp(:,:) = fract_tmp(:,:) * child%lcts(jlct)%fract(:,:)
                END IF
                ! Sum up the children's contributions
                tile_lct_fract_tmp(:,:) = tile_lct_fract_tmp(:,:) + fract_tmp(:,:)
                EXIT
              END IF
            END DO
            ! Get the (absolute) fraction of the current tile
            CALL tile%Get_fraction(fract=fract_tmp(:,:))
            ! Calculate the (relative) LCT fraction of the tile
            WHERE (tile_lct_fract_tmp(:,:) == 0._wp)
              tile%lcts(ilct)%fract(:,:) = 0._wp
            ELSE WHERE (fract_tmp(:,:) > 0._wp)
              tile%lcts(ilct)%fract(:,:) = tile_lct_fract_tmp(:,:) / fract_tmp(:,:)
            ELSE WHERE
              tile%lcts(ilct)%fract(:,:) = 1._wp
            END WHERE
          END IF
          child => child%Get_next_sibling_tile()
        END DO
        !$ACC UPDATE DEVICE(tile%lcts(ilct)%fract) ASYNC(1)
      END DO

      CALL model%Goto_next_tile(tile)

    END DO

    DEALLOCATE(fract_tmp)

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 4. Iterate over tiles again and initialize hierarchical structure for memories, variables, cq and lcc
    !>
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      ! IF (.NOT. tile%visited(no_omp_thread) .AND. tile%Has_children()) THEN
      IF (.NOT. tile%visited(no_omp_thread)) THEN
          CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      ! tile is not a leaf and we're on the way up
      CALL tile%Init_vars()
      CALL tile%Cache_GPU_pointers()

      CALL model%Goto_next_tile(tile)

    END DO

    ! Iterate over tiles, crosscheck tile structure and collect conserved quantities
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      IF (tile%Is_process_calculated(PHENO_)) THEN
        dsl4jsb_Get_config(PHENO_)
        IF (        tile%is_vegetation                                   &
            & .AND. TRIM(dsl4jsb_Config(PHENO_)%scheme) /= 'climatology' &
            & .AND. tile%lcts(1)%lib_id == 0)                            &
          & CALL finish(TRIM(routine), &
          &             'Must use phenology climatology on general, non-PFT VEG tile '//TRIM(tile%name))
      END IF

      CALL tile%Count_conserved_quantities()
      CALL tile%Collect_conserved_variables()

      CALL model%Goto_next_tile(tile)

    END DO

    ! Iterate over tiles, init lcc structures
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(FAGE_)) THEN
          CALL fage_init_lcc(tile)
        END IF

        IF (tile%Is_process_active(ALCC_)) THEN
          CALL alcc_init_lcc(tile)
        END IF

        IF (tile%Is_process_active(FLCC_)) THEN
          CALL flcc_init_lcc(tile)
        END IF

        IF (tile%Is_process_active(WLCC_)) THEN
          CALL wlcc_init_lcc(tile)
        END IF
      END SELECT

      CALL model%Goto_next_tile(tile)
    ENDDO

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> 5. Iterate over tiles again and print tile information (if debug is ON)
    !>
    IF (debug_on()) THEN
      CALL model%Reset_tiles()

      CALL model%Get_top_tile(tile)

      IF (.NOT. ASSOCIATED(tile)) &
        & CALL finish(TRIM(routine), 'Top tile not set')

      DO WHILE (ASSOCIATED(tile))

        IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
          CALL model%Goto_next_tile(tile)
          CYCLE
        END IF

        CALL tile%Print()

        CALL model%Goto_next_tile(tile)

      END DO
    END IF

! #ifdef _OPENACC
!     IF (.NOT. is_time_restart(get_time_start())) CALL model%cpu_to_gpu()
! #endif

  END SUBROUTINE jsbach_init


  ! ====================================================================================================== !
  !
  !> Do some more initialization after restart files have been read
  !
  SUBROUTINE jsbach_init_after_restart(model_id)

    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_, AGGREGATE_
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,     ONLY: MODEL_JSBACH
    USE mo_carbon_init,         ONLY: carbon_read_cpools
    USE mo_hydro_init,          ONLY: hydro_sanitize_state
    USE mo_util,                ONLY: one_of

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: model_id !< ID of model instance to be initialized
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_grid),  POINTER          :: grid
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                             :: no_omp_thread, blks, blke, blk, jcs, jce

    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(HYDRO_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_init_after_restart'
    ! -------------------------------------------------------------------------------------------------- !

    no_omp_thread = Get_omp_thread()

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    model => Get_model(model_id)
    grid  => Get_grid(model%grid_id)
    blks = grid%get_blk_start()
    blke = grid%get_blk_end()

    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      dsl4jsb_Get_config(CARBON_)
      dsl4jsb_Get_config(HYDRO_)
    END SELECT

    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        call model%Goto_next_tile(tile)
        CYCLE
      END IF

#ifndef __NO_QUINCY__
      IF (ASSOCIATED(tile%bgcm_store)) THEN
        DO blk = blks, blke
          jcs = grid%get_col_start(blk)
          jce = grid%get_col_end  (blk)
          CALL tile%bgcm_store%Store_bgcms_in_matrices(tile%name, jcs, jce, blk)
        ENDDO
      ENDIF
#endif

      ! Indicate to the weighted average aggregator that tile fractions might have changed after restart
      tile%l_fract_children_changed(:,:) = model%Do_fractions_change()

      ! TBD: Loop over processes instead of explicit
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        IF (tile%Is_process_active(CARBON_)) THEN
          IF (dsl4jsb_Config(CARBON_)%read_cpools) CALL carbon_read_cpools(tile)
        END IF
        IF (dsl4jsb_Config(HYDRO_)%sanitize_restart .AND. &
          & tile%Is_process_active(HYDRO_)) THEN
          CALL hydro_sanitize_state(tile)
        END IF
      END SELECT

      CALL model%Goto_next_tile(tile)

    END DO

#ifdef _OPENACC
    CALL model%cpu_to_gpu()
#endif

  END SUBROUTINE jsbach_init_after_restart
  ! ====================================================================================================== !

#endif
END MODULE mo_jsb_model_init
