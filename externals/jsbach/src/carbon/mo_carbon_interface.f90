!> Contains the interfaces to the carbon process
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_carbon_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  USE mo_carbon_constants,   ONLY: molarMassCO2_kg, sec_per_day,                    &
                                 & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, &
                                 & i_lctlib_nonsoluble, i_lctlib_humus

  ! Use of processes in this module (Get_carbon_memory and Get_carbon_config)
  dsl4jsb_Use_processes CARBON_, ASSIMI_, PHENO_, DISTURB_, A2L_, L2A_

  ! Use of process configurations
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(CARBON_)

  ! Use of process memories (t_carbon_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(L2A_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_carbon_tasks
  PUBLIC ::  recalc_carbon_per_tile_vars, recalc_per_tile_vars, calculate_current_c_ag_1_and_bg_sums,      &
    &        calculate_current_c_ta_state_sum, check_carbon_conservation, yday_carbon_conservation_test,   &
    &        rescale_carbon_upon_reference_area_change,                                                    &
    &        carbon_transfer_from_active_to_passive_vars_onChunk, global_carbon_diagnostics,               &
    &        rescale_carbon_vars_upon_ageing_induced_ref_area_change, merge_carbon_vars_upon_ageing

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_interface'

  !> Type definition for C_NPP_pot_allocation task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_C_NPP_pot_allocation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_C_NPP_pot_allocation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_C_NPP_pot_allocation !< Aggregates computed task variables
  END TYPE   tsk_C_NPP_pot_allocation

  !> Constructor interface for C_NPP_pot_allocation task
  INTERFACE tsk_C_NPP_pot_allocation
    PROCEDURE Create_task_C_NPP_pot_allocation        !< Constructor function for task
  END INTERFACE tsk_C_NPP_pot_allocation

  !> Type definition for ta_vars_after_lcc task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_ta_vars_after_lcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_ta_vars_after_lcc    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_ta_vars_after_lcc !< Aggregates computed task variables
  END TYPE   tsk_ta_vars_after_lcc

  !> Constructor interface for ta_vars_after_lcc task
  INTERFACE tsk_ta_vars_after_lcc
    PROCEDURE Create_task_ta_vars_after_lcc        !< Constructor function for task
  END INTERFACE tsk_ta_vars_after_lcc

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for carbon task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "carbon"
  !!
  FUNCTION Create_task_C_NPP_pot_allocation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_C_NPP_pot_allocation::return_ptr)
    CALL return_ptr%Construct(name='C_NPP_pot_allocation', process_id=CARBON_, owner_model_id=model_id)

  END FUNCTION Create_task_C_NPP_pot_allocation

  ! ====================================================================================================== !
  !
  !> Constructor for ta variables after lcc
  !
  FUNCTION Create_task_ta_vars_after_lcc(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_ta_vars_after_lcc::return_ptr)
    CALL return_ptr%Construct(name='ta_vars_after_lcc', process_id=CARBON_, owner_model_id=model_id)

  END FUNCTION Create_task_ta_vars_after_lcc

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for carbon process
  !!
  !! @param[in,out] this      Instance of carbon process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_carbon_tasks(this, model_id)

    ! in/out
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_C_NPP_pot_allocation(model_id))
    CALL this%Register_task(tsk_ta_vars_after_lcc(model_id))

  END SUBROUTINE Register_carbon_tasks

  ! ================================================================================================================================
  !>
  !> Implementation to calculate the allocation of NPP carbon to the plant and soil carbon pools
  !!        R: Corresponds to update_cbalance_bethy in JSBACH3.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_C_NPP_pot_allocation(tile, options)

    ! Use declarations
    USE mo_carbon_process,    ONLY: calc_Cpools
    USE mo_jsb_time,          ONLY: is_newday, is_newyear, timesteps_per_day, is_time_experiment_start

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare local variables
    TYPE(t_jsb_model), POINTER   :: model
    LOGICAL                      :: lstart, new_day, new_year
    INTEGER                      :: iblk, ics, ice, nc, ic
    REAL(wp)                     :: dtime
    REAL(wp), allocatable        :: precip_total(:)
    REAL(wp), DIMENSION(options%nc) :: MaxLai
    REAL(wp), DIMENSION(options%nc) :: old_c_state_sum_ta, current_fluxes
    CHARACTER(len=*), PARAMETER  :: routine = modname//':update_C_NPP_pot_allocation'

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(A2L_)

    ! Declare pointers to variables in memory
    ! R: former cbalance_type variables
    dsl4jsb_Real2D_onChunk ::  pseudo_temp
    dsl4jsb_Real2D_onChunk ::  N_pseudo_temp
    dsl4jsb_Real2D_onChunk ::  F_pseudo_temp
    dsl4jsb_Real2D_onChunk ::  pseudo_temp_yDay
    dsl4jsb_Real2D_onChunk ::  pseudo_precip
    dsl4jsb_Real2D_onChunk ::  pseudo_precip_yDay
    dsl4jsb_Real2D_onChunk ::  N_pseudo_precip
    dsl4jsb_Real2D_onChunk ::  F_pseudo_precip

    dsl4jsb_Real2D_onChunk ::  cconservation_calcCpools ! C conservation test: Deviation from conservation

    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter           ! Carbon flux from the veget. to the litter pools
                                                     ! [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  c_green           ! C-pool for leaves, fine roots, vegetative organs and
                                                     ! other green (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_reserve         ! C-pool for carbohydrate reserve (sugars, starches) that
                                                     ! allows plants to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_woods           ! C-pool for stems, thick roots and other (dead) structural
                                                     !  material of living plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest    ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]

    !SIZE CLASS 1
    dsl4jsb_Real2D_onChunk ::  c_acid_ag1        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag1       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg1       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1  ! Yasso below ground litter-pool for non-soluble litter
                                                      ! [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_1         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)/m^2(canopy)]
    !SIZE CLASS 2
    dsl4jsb_Real2D_onChunk ::  c_acid_ag2        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag2       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg2       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2  ! Yasso below ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_2         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)/m^2(canopy)]

    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_1_sum  ! Annual sum of humus decomposed in YASSO (leaf)
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_2_sum  ! Annual sum of humus decomposed in YASSO (wood)
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_into_humus_1_sum    ! Annual sum of cflux into humus in YASSO (leaf)
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_into_humus_2_sum    ! Annual sum of cflux into humus in YASSO (wood)
                                                      !  [mol(C)/m^2(canopy)]

    dsl4jsb_Real2D_onChunk ::  gross_assimilation_ca  ! Gross assimilation of all canopy layers together
    dsl4jsb_Real2D_onChunk ::  NPP_pot_rate_ca        ! The instantaneous (potential) NPP rate [mol(C)/m^2(canopy) s]
    !dsl4jsb_Real2D_onChunk ::  NPP_pot_Rate_acc      ! averaged NPP rate over one day [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  LAI_sum                ! used to accumulate LAI over a day. Sum of LAI since Midnight.
    dsl4jsb_Real2D_onChunk ::  NPP_pot_sum            ! used to accumulated NPP-Rate over a day
    dsl4jsb_Real2D_onChunk ::  GPP_sum                ! used to accumulated GPP-Rate over a day
    dsl4jsb_Real2D_onChunk ::  LAI_yDayMean           ! previous days mean
    dsl4jsb_Real2D_onChunk ::  LAI_yyDayMean          ! previous previous days mean
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean       ! mean value of NPP-Rate yesterday (from NPP_pot_sum())
                                                      ! [mol(CO2)/(m^2(canopy) s)]
                                                      ! = cbalance%NPP_pot_sum(...)/time_steps_per_day
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean       ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation
                                                      ! Actual NPP after N-limitation and excess carbon drop.
                                                      ! [mol(CO2)/(m^2(canopy) s)]
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean

    dsl4jsb_Real2D_onChunk ::  soil_respiration       ! mean daily rate of heterotrophic (soil) respiration
                                                      ! [mol(CO2)/m^2(ground)].
                                                      ! Without N limitation!
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction    ! Daily updated flux correction from yesterdays carbon balance
                                                      ! [mol(CO2)/m^2(canopy) s]
                                                      ! Amount by which the NPP rate entering the routine has to be corrected.
                                                      ! This correction arises either because otherwise the reserve pool would
                                                      ! get negative (positive correction), or the wood pool would exceed its
                                                      ! maximum value (negative correction).
    dsl4jsb_Real2D_onChunk ::  excess_NPP             ! That part of NPP that because of structural limits could not be
                                                      ! allocated in carbon pools [mol(CO2)/m^2(canopy) s]

    dsl4jsb_Real2D_onChunk ::  root_exudates          ! Total root exudates entering to the litter green pools
                                                      ! [mol(C)/m^2(canopy) s]

    dsl4jsb_Real2D_onChunk ::  c_sum_veg_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_litter_ag_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_litter_bg_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_humus_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_natural_ta

    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb     !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_littergreen !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_atm         !
    dsl4jsb_Real2D_onChunk ::  co2flux_npp_2_atm_yday_ta  ! day CO2 flux from actual NPP, required for cconservation test
    dsl4jsb_Real2D_onChunk ::  co2flux_npp_2_atm_ta       ! grid cell averages of net CO2 fluxes between
    dsl4jsb_Real2D_onChunk ::  co2flux_soilresp_2_atm_ta  ! .. biosphere (due to NPP, soil respiration and
    dsl4jsb_Real2D_onChunk ::  co2flux_herb_2_atm_ta      ! .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)]


    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !dsl4jsb_Real2D_onChunk :: fract_litter_wood_new
    !
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: LAI
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: t_air                       ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow

    ! variables needed for NLCC process
    dsl4jsb_Real2D_onChunk :: max_green_bio
    dsl4jsb_Real2D_onChunk :: sla

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(CARBON_)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(A2L_)

    ! Set process variables
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_temp)                ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_temp_yDay)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  N_pseudo_temp)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  F_pseudo_temp)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_precip)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_precip_yDay)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  N_pseudo_precip)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  F_pseudo_precip)            ! inout

    !
    ! This is on canopy area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cconservation_calcCpools )  ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter ) ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)         ! inout
    !SIZE CLASS 1
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! inout
    !SIZE CLASS 2
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)            ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum)   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum)   ! inout
    ! xx
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_ca )     ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  NPP_pot_rate_ca )           ! inout
    !dsl4jsb_Get_var2D_onChunk(ASSIMI_,  NPP_pot_Rate_acc )         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_sum )                   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_sum )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_sum )                   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_yDayMean )              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_yyDayMean )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean )              ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction)        ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  excess_NPP)                 ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates )             ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_veg_ta )               ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_litter_ag_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_litter_bg_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_humus_ta )             ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_natural_ta )           ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb )       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_littergreen  )  ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_atm )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_npp_2_atm_yday_ta )  ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_npp_2_atm_ta )       ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_soilresp_2_atm_ta )  ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_herb_2_atm_ta )      ! out

    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !dsl4jsb_Get_var2D_onChunk(CARBON_,  fract_litter_wood_new )
    !
    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )      ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   LAI )                       ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     t_air)                      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     rain)                       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     snow)                       ! in

    ! variables needed for NLCC process
    dsl4jsb_Get_var2D_onChunk(CARBON_,  max_green_bio)
    dsl4jsb_Get_var2D_onChunk(CARBON_,  sla)

    lstart   = is_time_experiment_start(options%current_datetime)
    new_day  = is_newday(options%current_datetime, dtime)
    new_year = is_newyear(options%current_datetime, dtime)

    ! ---------------------------
    ! Go

    ! Put a comment on the screen...
    IF (debug_on() .AND. iblk==1 .AND. new_day  ) CALL message(TRIM(routine), &
                                                           'update_C_NPP_pot_allocation: First time step of new day')
    IF (debug_on() .AND. iblk==1 .AND. new_year ) CALL message(TRIM(routine), &
                                                           'update_C_NPP_pot_allocation: First time step of new year')


    ! Initializations
    co2flux_npp_2_atm_yday_ta  = 0._wp
    co2flux_npp_2_atm_ta       = 0._wp
    co2flux_soilresp_2_atm_ta  = 0._wp
    co2flux_herb_2_atm_ta = 0._wp

    IF( new_year .OR. lstart) THEN ! reset annual sums if newyear or if just started exp
      c_decomp_humus_1_sum =  0._wp
      c_decomp_humus_2_sum =  0._wp
      c_into_humus_1_sum =  0._wp
      c_into_humus_2_sum =  0._wp
    ENDIF

    ! R: "IF (read_cpools .AND. (lstart .OR. lresume)) THEN ..." of JS3 is in JS4 now in SUBROUTINE carbon_init_ic(tile)!

    ! R: To keep in mind: In JS3: NPP_pot_Rate is calculated from function NPP_pot_rate_bethy. This function gets
    !    grossAssimilation from theLand%Bethy%gross_assimilation,
    !    which is GPP relative to ground area [MOL(CO2) / M^2(ground) S] (mean value), WHICH MEANS CANOPY AREA!!!
    !    -> cbalance%NPP_pot_sum = cbalance%NPP_pot_sum + cbalance%NPP_pot_rate.
    !    AND: -> cbalance_diag%NPP_pot_yDayMean = cbalance%NPP_pot_sum/time_steps_per_day. -> cbalance_diag%NPP_pot_yDayMean
    !    is used in update_Cpools to calculate the carbon storage of the pools.
    ! R: Should be arbitary now:
    ! Compute net primary production rate
    ! R: In JS4 ist fogende Zeile jetzt im assimi_process.f90: calc_NPP_pot_rate
    !    cbalance%NPP_pot_rate(kidx0:kidx1,:) = NPP_pot_rate_bethy(grossAssimilation(1:nidx,:),darkRespiration(1:nidx,:))
    ! R: NPP_pot_Rate_acc unterscheidet sich in JS3 von NPP_pot_sum darin, dass NPP_pot_rate_acc zu jedem Zeitschritt eine
    !    momentane Aufaddierung vom NPP ist aber NPP_pot_sum immer eine Tagessumme ist.
    !    NPP_pot_Rate_acc wurde nur in den Output geschrieben aber nie weiter verwendet, daher habe ich es raus genommen.
    !    Wenn man es doch drin haben will:
    !NPP_pot_Rate_ca_acc =NPP_pot_Rate_ca_acc  +   NPP_pot_rate_ca / time_steps_per_day
    ! R: The following line should be handeled differently in JS4:
    !    In JS3 the areaWeightingFactor rescales the canopy area up to the whole Gridbox. E.g. for the C pools.
    !    In JS4 this upscaling has normally to be done by the accumulation procedures. However, the canopy area
    !    still has to be scaled to the PFT tile area for JS4. As the C pools are calculated on m^2 canopy also in JS4
    !    => Cpool(on PFT tile) =Cpool(canopy) * veg_fract_correction.
    !
    ! Therefore this gets arbitary:
    ! Prepare area weighting factor to rescale from 1/[m^2(canopy)] to 1/[m^2(grid box)]
    ! areaWeightingFactor(:,:) = veg_fract_correction(:,:) * surface%cover_fract(kidx0:kidx1,:) &
    !                            * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,NCOPIES=ntiles)


    ! Update pseudo-15day-mean air temperature and precipitation for pseudo_temp_yDay and pseudo_precip_yDay
    ! These variables are needed for the "CALL calc_Cpools", where they given to yasso to calculte C decomposition.
    ! pseudo_temp is a temperature that depicts the smoth course of the soil temperature following the air temperature.
    pseudo_temp   = t_air * N_pseudo_temp  +  F_pseudo_temp   * pseudo_temp
    ! R: I would like to put this precip_total with another name into hydro memory and calculate this
    !    in HYDRO and here just load it from there...
    precip_total  = rain + snow  ! Precipitation rate [kg/(m^2 s)]
    pseudo_precip = precip_total * N_pseudo_precip  +  F_pseudo_precip * pseudo_precip

    ! All per tile variables need to be re-calculated in case they are written to output
    CALL recalc_per_tile_vars(tile, options,                                      &
      & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
      & c_crop_harvest = c_crop_harvest(:),                                       &
      & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
      & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
      & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
      & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
      & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
      & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
      & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
      & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
      & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:), root_exudates = root_exudates(:),   &
      & soil_respiration = soil_respiration(:), NPP_flux_correction = NPP_flux_correction(:),   &
      & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),                             &
      & cflux_c_green_2_herb = cflux_c_green_2_herb(:), NPP_act_yDayMean = NPP_act_yDayMean(:), &
      & NPP_pot_yDayMean = NPP_pot_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

    IF( .NOT. new_day .OR. lstart) THEN ! perform daily sums if WE ARE NOT STARTING A NEW DAY
                                            ! or if we JUST STARTED AN TOTALLY NEW EXPERIMENT

      LAI_sum          = LAI_sum          + LAI
      NPP_pot_sum      = NPP_pot_sum      + NPP_pot_rate_ca
      GPP_sum          = GPP_sum          + gross_assimilation_ca

    ELSE                                    ! A NEW DAY BEGINS and we DID NOT START AN
                                            ! TOTALLY NEW EXPERIMENT ==> perform carbon balance

      ! First save as "previous previous days" and "previous days" means
      LAI_yyDayMean  = LAI_yDayMean
      LAI_yDayMean   = LAI_sum/timesteps_per_day(dtime)

      NPP_pot_yDayMean     = NPP_pot_sum/timesteps_per_day(dtime)
      GPP_yDayMean         = GPP_sum/timesteps_per_day(dtime)

      ! Then restart summing of this days values
      LAI_sum             = LAI
      NPP_pot_sum         = NPP_pot_rate_ca
      GPP_sum             = gross_assimilation_ca
      !NPP_pot_rate_ca_acc        = NPP_pot_rate_ca/timesteps_per_day(dtime) ! s.o.

      ! Save other variables (for yasso) as "previous days"
      pseudo_temp_yDay   = pseudo_temp
      pseudo_precip_yDay = pseudo_precip

      ! In case of l_forestRegrowth = true the maximum LAI can change over time
      IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
         MaxLai(:) = dsl4jsb_var2D_onChunk(PHENO_, maxLAI_allom )
      ELSE
        ! R: Rename MaxLAI in lctlib to LAI_max in JS3 comes from theLand%Vegetation%LAI_max,
        ! which is for each vegetated tile taken from MaxLAI in the lct library.
        MaxLai(:) = dsl4jsb_Lctlib_param(MaxLAI)
      END IF

      !Calculate current sum before operation for c conservation test
      CALL recalc_carbon_per_tile_vars(tile, options)
      CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

      ! Now update the C pools
      CALL calc_Cpools( &
        &  LAI_yDayMean,                                      & ! in
        &  LAI_yyDayMean,                                     & ! in
        &  MaxLai,                                            & ! in
        &  NPP_pot_yDayMean,                                  & ! in
        !
        &  dsl4jsb_Lctlib_param(fract_npp_2_woodPool),        & ! in
        &  dsl4jsb_Lctlib_param(fract_NPP_2_reservePool),     & ! in
        &  dsl4jsb_Lctlib_param(fract_NPP_2_exudates),        & ! in
        !
        &  dsl4jsb_Lctlib_param(fract_green_2_herbivory),     & ! in
        &  dsl4jsb_Lctlib_param(tau_c_woods),                 & ! in
        &  dsl4jsb_Lctlib_param(LAI_shed_constant),           & ! in
        &  dsl4jsb_Lctlib_param(Max_C_content_woods),         & ! in
        &  dsl4jsb_Lctlib_param(specificLeafArea_C),          & ! in
        &  dsl4jsb_Lctlib_param(reserveC2leafC),              & ! in
        &  dsl4jsb_Lctlib_param(PhenologyType),               & ! in  R: is it crop? If yes then this is 5!
        !
        &  c_green,                                           & ! inout
        &  c_woods,                                           & ! inout
        &  c_reserve,                                         & ! inout
        &  c_crop_harvest,                                    & ! inout
        !
        &  soil_respiration,                                  & ! out
        &  NPP_flux_correction,                               & ! out
        &  excess_NPP,                                        & ! out
        &  root_exudates,                                     & ! out
        &  cflux_c_greenwood_2_litter,                        & ! out
        &  cflux_c_green_2_herb,                              & ! out
        &  cflux_herb_2_littergreen,                          & ! out
        &  cflux_herb_2_atm,                                  & ! out
        &  NPP_act_yDayMean,                                  & ! out
        !
        !&  fract_litter_wood_new    = fract_litter_wood_new,       & ! inout, optional. Only for spitfire
        !
        ! variables only needed with yasso:
        &  temp2_30d        =  pseudo_temp_yDay,                      & ! in
        &  precip_30d       =  pseudo_precip_yDay,                    & ! in
                             !
        &  c_acid_ag1       =  c_acid_ag1,                            & ! inout
        &  c_water_ag1      =  c_water_ag1,                           & ! inout
        &  c_ethanol_ag1    =  c_ethanol_ag1,                         & ! inout
        &  c_nonsoluble_ag1 =  c_nonsoluble_ag1,                      & ! inout
        &  c_acid_bg1       =  c_acid_bg1,                            & ! inout
        &  c_water_bg1      =  c_water_bg1,                           & ! inout
        &  c_ethanol_bg1    =  c_ethanol_bg1,                         & ! inout
        &  c_nonsoluble_bg1 =  c_nonsoluble_bg1,                      & ! inout
        &  c_humus_1        =  c_humus_1,                             & ! inout
        &  c_acid_ag2       =  c_acid_ag2,                            & ! inout
        &  c_water_ag2      =  c_water_ag2,                           & ! inout
        &  c_ethanol_ag2    =  c_ethanol_ag2,                         & ! inout
        &  c_nonsoluble_ag2 =  c_nonsoluble_ag2,                      & ! inout
        &  c_acid_bg2       =  c_acid_bg2,                            & ! inout
        &  c_water_bg2      =  c_water_bg2,                           & ! inout
        &  c_ethanol_bg2    =  c_ethanol_bg2,                         & ! inout
        &  c_nonsoluble_bg2 =  c_nonsoluble_bg2,                      & ! inout
        &  c_humus_2        =  c_humus_2,                             & ! inout
                              !
        &  LeafLit_coef_acid       =  dsl4jsb_Lctlib_param(LeafLit_coef(i_lctlib_acid)),       & ! in
        &  LeafLit_coef_water      =  dsl4jsb_Lctlib_param(LeafLit_coef(i_lctlib_water)),      & ! in
        &  LeafLit_coef_ethanol    =  dsl4jsb_Lctlib_param(LeafLit_coef(i_lctlib_ethanol)),    & ! in
        &  LeafLit_coef_nonsoluble =  dsl4jsb_Lctlib_param(LeafLit_coef(i_lctlib_nonsoluble)), & ! in
        &  LeafLit_coef_humus      =  dsl4jsb_Lctlib_param(LeafLit_coef(i_lctlib_humus)),      & ! in
        &  WoodLit_coef_acid       =  dsl4jsb_Lctlib_param(WoodLit_coef(i_lctlib_acid)),       & ! in
        &  WoodLit_coef_water      =  dsl4jsb_Lctlib_param(WoodLit_coef(i_lctlib_water)),      & ! in
        &  WoodLit_coef_ethanol    =  dsl4jsb_Lctlib_param(WoodLit_coef(i_lctlib_ethanol)),    & ! in
        &  WoodLit_coef_nonsoluble =  dsl4jsb_Lctlib_param(WoodLit_coef(i_lctlib_nonsoluble)), & ! in
        &  WoodLit_coef_humus      =  dsl4jsb_Lctlib_param(WoodLit_coef(i_lctlib_humus)),      & ! in
        &  WoodLitterSize          =  dsl4jsb_Lctlib_param(WoodLitterSize),                    & ! in

        & c_decomp_humus_1_sum  = c_decomp_humus_1_sum,               & ! out
        & c_decomp_humus_2_sum  = c_decomp_humus_2_sum,               & ! out
        & c_into_humus_1_sum    = c_into_humus_1_sum,                 & ! out
        & c_into_humus_2_sum    = c_into_humus_2_sum                  & ! out
        & )

      ! determine annual maximum content of green pool for NLCC process
      IF (new_year) THEN
         max_green_bio = c_green
      ELSEIF (new_day) THEN
         max_green_bio = MAX(max_green_bio, c_green)
      ENDIF

      ! Compute carbon contents from canopy to tile box area by weighting pools with fractions of grid box covered by vegetation
      ! ------------------------------------------------------------------------------------------------------------------------

      CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
      CALL recalc_per_tile_vars(tile, options,                                      &
        & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
        & c_crop_harvest = c_crop_harvest(:),                                       &
        & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
        & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
        & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
        & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
        & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
        & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
        & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
        & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
        & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
        & c_decomp_humus_1_sum = c_decomp_humus_1_sum(:),                           &
        & c_decomp_humus_2_sum = c_decomp_humus_2_sum(:),                           &
        & c_into_humus_1_sum = c_into_humus_1_sum(:),                               &
        & c_into_humus_2_sum = c_into_humus_2_sum(:), root_exudates = root_exudates(:),           &
        & soil_respiration = soil_respiration(:), NPP_flux_correction = NPP_flux_correction(:),   &
        & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),                             &
        & cflux_c_green_2_herb = cflux_c_green_2_herb(:), NPP_act_yDayMean = NPP_act_yDayMean(:), &
        & NPP_pot_yDayMean = NPP_pot_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

      ! cflux_herb_2_atm and soil_respiration are negative (fluxes away from land)
      ! thus + here instead of - to reduce NPP by these fluxes
      current_fluxes = (NPP_act_yDayMean + cflux_herb_2_atm + soil_respiration) &
        & * sec_per_day * veg_fract_correction * fract_fpc_max
      CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
        & current_fluxes(:), cconservation_calcCpools(:))

    END IF ! new day

    ! Compute net CO2 fluxes exchanged with atmosphere at each time step
    !-------------------------------------------------------------------
    ! Note: carbon loss of biosphere means a positive CO2 flux to atmosphere (i.e. NEP and net CO2-flux have opposite signs)

    co2flux_npp_2_atm_ta = molarMassCO2_kg *                    & ! Conversion factor from mol to kg CO2
               (- veg_fract_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                    * (NPP_pot_rate_ca                          & ! current (not actual) NPP rate
                        - (NPP_pot_yDayMean-NPP_act_yDayMean)) )  ! corrected with yesterdays actual NPP defizit

    co2flux_soilresp_2_atm_ta = molarMassCO2_kg *               & ! Conversion factor from mol to kg CO2
               (- veg_fract_correction *fract_fpc_max           & ! Minus: atmosphere gain is positive
                    * soil_respiration      )                    ! .. soil respiration

    co2flux_herb_2_atm_ta = molarMassCO2_kg *                   & ! Conversion factor from mol to kg CO2
               (- veg_fract_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                    * cflux_herb_2_atm )                          ! .. herbivory

    co2flux_npp_2_atm_yday_ta = molarMassCO2_kg *               & ! daily C conservation test cannot deal with
              (- veg_fract_correction * fract_fpc_max           & ! diurnal cycle of NPP,
                   * (NPP_act_yDayMean))                          ! therefore here additionaly the day mean

    ! Calculate diagnostic carbon sums
    CALL calculate_current_c_ta_state_sum(tile, options, c_sum_natural_ta(:),            &
      & c_sum_veg_ta(:), c_sum_litter_ag_ta(:), c_sum_litter_bg_ta(:), c_sum_humus_ta(:))

    ! write lct information on variable to make it available on pft-tiles via function collect_var for NLCC process
    sla(:) = dsl4jsb_Lctlib_param(specificLeafArea_C)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')


  END SUBROUTINE update_C_NPP_pot_allocation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "C_NPP_pot_allocation"
  !!
  !! @param[in,out] tile    Tile for which aggregation of child tiles is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_C_NPP_pot_allocation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_C_NPP_pot_allocation'

    INTEGER  :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean,          weighted_by_fract)

    !@todo Those commented out are currently not calculated
    ! non-yasso-pools
    dsl4jsb_Aggregate_onChunk(CARBON_, c_green_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_woods_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_reserve_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_crop_harvest_ta,     weighted_by_fract)
    ! ag1 and bg1
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_ag1_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_ag1_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_ag1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_ag1_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg1_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg1_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg1_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_1_ta,          weighted_by_fract)
    !ag2 and bg2
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_ag2_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_ag2_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_ag2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_ag2_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg2_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg2_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg2_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_2_ta,          weighted_by_fract)

    !JN: vars for analytically equilibrating YASSO humus pools outside of jsb4 simulations (compare jsb3 svn ref 9652)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_decomp_humus_1_sum_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_decomp_humus_2_sum_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_into_humus_1_sum_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_into_humus_2_sum_ta,          weighted_by_fract)

    ! sums
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_veg_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_litter_ag_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_litter_bg_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_humus_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_natural_ta,      weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_sum_ta,        weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_sum_ta,            weighted_by_fract)
    ! fluxes
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_greenwood_2_litter_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, soil_respiration_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_flux_correction_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, excess_NPP_ta,            weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yyDayMean_ta,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_act_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, root_exudates_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_green_2_herb_ta,  weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_littergreen_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, fract_litter_wood_new_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cconservation_calcCpools, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_npp_2_atm_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_soilresp_2_atm_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_herb_2_atm_ta,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_npp_2_atm_yday_ta, weighted_by_fract)

    ! variables for NLCC process
    dsl4jsb_Aggregate_onChunk(CARBON_, max_green_bio,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, sla,                      weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_C_NPP_pot_allocation

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "ta_vars_after_lcc"
  !>
  !> Task "ta_vars_after_lcc" updates the ta variables
  !
  SUBROUTINE update_ta_vars_after_lcc(tile, options)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !

    ! ... currently not required, but could contain something like:
    ! CALL recalc_carbon_per_tile_vars(tile, options)
    ! CALL calculate_current_c_ag_1_and_bg_sums(tile, options)

  END SUBROUTINE update_ta_vars_after_lcc

  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "ta_vars_after_lcc"
  !>
  !> Aggregate for task "ta_vars_after_lcc" is currently required to aggregate the ta for the co2flux
  !
  ! -------------------------------------------------------------------------------------------------------!
  SUBROUTINE aggregate_ta_vars_after_lcc(tile, options)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_ta_vars_after_lcc'

    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! if running with FLCC (which runs on the veg tile)
    ! the co2flux_fire_all_2_atm_ta variable is changed on the leaves
    ! ... it therefore requires a carbon aggregation task on the leaves to run subsequent to the flcc task
    dsl4jsb_Aggregate_onChunk(CARBON_,  co2flux_fire_all_2_atm_ta, weighted_by_fract)

    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_green_ta,              weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_woods_ta,              weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_reserve_ta,            weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_green_2_soil_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag1_ta,           weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag1_ta,          weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag1_ta,        weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag1_ta,     weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag2_ta,           weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag2_ta,          weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag2_ta,        weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag2_ta,     weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_1_ta,            weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_2_ta,            weighted_by_fract)


  END SUBROUTINE aggregate_ta_vars_after_lcc


!##########################################################################################################################
! other subroutines (no further tasks)
!##########################################################################################################################

  ! ================================================================================================================================
  !>
  !> scales all carbon state variables upon reference area change
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to scale from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE rescale_carbon_upon_reference_area_change(tile, options, oldRefArea, newRefArea)

    USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(wp),                   INTENT(in)            :: oldRefArea(:),  &
                                                         newRefArea(:)

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':rescale_carbon_upon_reference_area_change'
    INTEGER                     :: nc, ic, ics, ice, iblk

    dsl4jsb_Real2D_onChunk ::  c_green
    dsl4jsb_Real2D_onChunk ::  c_reserve
    dsl4jsb_Real2D_onChunk ::  c_woods
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1
    dsl4jsb_Real2D_onChunk ::  c_water_ag1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1
    dsl4jsb_Real2D_onChunk ::  c_water_bg1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk ::  c_humus_1

    dsl4jsb_Real2D_onChunk ::  c_acid_ag2
    dsl4jsb_Real2D_onChunk ::  c_water_ag2
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2
    dsl4jsb_Real2D_onChunk ::  c_water_bg2
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2
    dsl4jsb_Real2D_onChunk ::  c_humus_2

    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_1_sum
    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_2_sum
    dsl4jsb_Real2D_onChunk ::  c_into_humus_1_sum
    dsl4jsb_Real2D_onChunk ::  c_into_humus_2_sum

    ! Assertion: 0 < oldRedArea <= 1 and 0 < newRefArea <= 1
    ! TODO on GPU
#ifndef _OPENACC
    IF ( ANY(oldRefArea .LE. 0) .OR. ANY(oldRefArea .GT. 1) ) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: Reference areas (here old) need to be >0 and <=1 - min: '&
        & //real2string(MINVAL(oldRefArea)) //'; max: '//real2string(MAXVAL(oldRefArea))  )
    ENDIF
    IF ( ANY(newRefArea .LE. 0) .OR. ANY(newRefArea .GT. 1) ) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: Reference areas (here new) need to be >0 and <=1 - min: '&
        & //real2string(MINVAL(newRefArea)) //'; max: '//real2string(MAXVAL(newRefArea))  )
    ENDIF
#endif

    dsl4jsb_Get_memory(CARBON_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)         ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)            ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum)   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum)   ! inout

    c_green           = c_green           * (oldRefArea / newRefArea)
    c_woods           = c_woods           * (oldRefArea / newRefArea)
    c_reserve         = c_reserve         * (oldRefArea / newRefArea)
    c_crop_harvest    = c_crop_harvest    * (oldRefArea / newRefArea)

    c_acid_ag1       = c_acid_ag1       * (oldRefArea / newRefArea)
    c_water_ag1      = c_water_ag1      * (oldRefArea / newRefArea)
    c_ethanol_ag1    = c_ethanol_ag1    * (oldRefArea / newRefArea)
    c_nonsoluble_ag1 = c_nonsoluble_ag1 * (oldRefArea / newRefArea)
    c_acid_bg1       = c_acid_bg1       * (oldRefArea / newRefArea)
    c_water_bg1      = c_water_bg1      * (oldRefArea / newRefArea)
    c_ethanol_bg1    = c_ethanol_bg1    * (oldRefArea / newRefArea)
    c_nonsoluble_bg1 = c_nonsoluble_bg1 * (oldRefArea / newRefArea)
    c_humus_1        = c_humus_1        * (oldRefArea / newRefArea)
    c_acid_ag2       = c_acid_ag2       * (oldRefArea / newRefArea)
    c_water_ag2      = c_water_ag2      * (oldRefArea / newRefArea)
    c_ethanol_ag2    = c_ethanol_ag2    * (oldRefArea / newRefArea)
    c_nonsoluble_ag2 = c_nonsoluble_ag2 * (oldRefArea / newRefArea)
    c_acid_bg2       = c_acid_bg2       * (oldRefArea / newRefArea)
    c_water_bg2      = c_water_bg2      * (oldRefArea / newRefArea)
    c_ethanol_bg2    = c_ethanol_bg2    * (oldRefArea / newRefArea)
    c_nonsoluble_bg2 = c_nonsoluble_bg2 * (oldRefArea / newRefArea)
    c_humus_2        = c_humus_2        * (oldRefArea / newRefArea)

    c_decomp_humus_1_sum = c_decomp_humus_1_sum * (oldRefArea / newRefArea)
    c_decomp_humus_2_sum = c_decomp_humus_2_sum * (oldRefArea / newRefArea)
    c_into_humus_1_sum   = c_into_humus_1_sum   * (oldRefArea / newRefArea)
    c_into_humus_2_sum   = c_into_humus_2_sum   * (oldRefArea / newRefArea)

    CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
    CALL recalc_per_tile_vars(tile, options, &
      & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:), &
      & c_crop_harvest = c_crop_harvest(:), &
      & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:), &
      & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
      & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:), &
      & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
      & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:), &
      & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
      & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:), &
      & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
      & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
      & c_decomp_humus_1_sum  = c_decomp_humus_1_sum(:),                          &
      & c_decomp_humus_2_sum  = c_decomp_humus_2_sum(:),                          &
      & c_into_humus_1_sum    = c_into_humus_1_sum(:),                            &
      & c_into_humus_2_sum    = c_into_humus_2_sum(:))

  END SUBROUTINE rescale_carbon_upon_reference_area_change


  ! ====================================================================================================== !
  !
  !> Calculates all carbon tile variables (known to this function) from current canopy states
  !
  SUBROUTINE recalc_carbon_per_tile_vars(tile, options)

  ! -------------------------------------------------------------------------------------------------- !
  CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
  TYPE(t_jsb_task_options),   INTENT(in)            :: options
  ! -------------------------------------------------------------------------------------------------- !

  dsl4jsb_Real2D_onChunk ::  c_green
  dsl4jsb_Real2D_onChunk ::  c_woods
  dsl4jsb_Real2D_onChunk ::  c_reserve
  dsl4jsb_Real2D_onChunk ::  c_crop_harvest

  dsl4jsb_Real2D_onChunk ::  c_acid_ag1
  dsl4jsb_Real2D_onChunk ::  c_water_ag1
  dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1
  dsl4jsb_Real2D_onChunk ::  c_acid_bg1
  dsl4jsb_Real2D_onChunk ::  c_water_bg1
  dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1

  dsl4jsb_Real2D_onChunk ::  c_acid_ag2
  dsl4jsb_Real2D_onChunk ::  c_water_ag2
  dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2
  dsl4jsb_Real2D_onChunk ::  c_acid_bg2
  dsl4jsb_Real2D_onChunk ::  c_water_bg2
  dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2

  dsl4jsb_Real2D_onChunk ::  c_humus_1
  dsl4jsb_Real2D_onChunk ::  c_humus_2

  dsl4jsb_Real2D_onChunk ::  cflux_dist_green_2_soil
  dsl4jsb_Real2D_onChunk ::  cflux_dist_woods_2_soil
  dsl4jsb_Real2D_onChunk ::  co2flux_fire_all_2_atm

  dsl4jsb_Real2D_onChunk ::  soil_respiration
  dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean
  dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean
  dsl4jsb_Real2D_onChunk ::  NPP_flux_correction
  dsl4jsb_Real2D_onChunk ::  GPP_yDayMean
  dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter
  dsl4jsb_Real2D_onChunk ::  root_exudates
  dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb

  dsl4jsb_Def_memory(CARBON_)

  INTEGER                     :: ics, ice, iblk
  ! -------------------------------------------------------------------------------------------------- !
  dsl4jsb_Get_memory(CARBON_)

  iblk  = options%iblk
  ics   = options%ics
  ice   = options%ice

  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_green_2_soil)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_woods_2_soil)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_fire_all_2_atm)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb)

  CALL recalc_per_tile_vars(tile, options,                                    &
  & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
  & c_crop_harvest = c_crop_harvest(:),                                       &
  & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
  & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
  & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
  & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
  & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
  & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
  & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
  & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
  & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:), root_exudates = root_exudates(:),   &
  & soil_respiration = soil_respiration(:), NPP_flux_correction = NPP_flux_correction(:),   &
  & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),                             &
  & cflux_dist_green_2_soil = cflux_dist_green_2_soil(:),                                   &
  & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:),                                   &
  & co2flux_fire_all_2_atm = co2flux_fire_all_2_atm(:),                                     &
  & cflux_c_green_2_herb = cflux_c_green_2_herb(:), NPP_act_yDayMean = NPP_act_yDayMean(:), &
  & NPP_pot_yDayMean = NPP_pot_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

END SUBROUTINE recalc_carbon_per_tile_vars


  ! ================================================================================================================================
  !>
  !> Calculates the per tile states from current canopy states
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE recalc_per_tile_vars(tile, options,                  &
      & c_green , c_woods, c_reserve, c_crop_harvest,             &
      & c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1, &
      & c_acid_bg1, c_water_bg1, c_ethanol_bg1, c_nonsoluble_bg1, &
      & c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2, &
      & c_acid_bg2, c_water_bg2, c_ethanol_bg2, c_nonsoluble_bg2, &
      & c_humus_1, c_humus_2,                                     &
      & c_decomp_humus_1_sum, c_decomp_humus_2_sum, c_into_humus_1_sum, c_into_humus_2_sum, &
      & soil_respiration, NPP_flux_correction, NPP_pot_yDayMean, NPP_act_yDayMean,          &
      & GPP_yDayMean, root_exudates, cflux_c_green_2_herb,                                  &
      & cflux_dist_green_2_soil, cflux_dist_woods_2_soil, co2flux_fire_all_2_atm,           &
      & cflux_c_greenwood_2_litter)

    USE mo_carbon_process,    ONLY: get_per_tile

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_green(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_woods(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_reserve(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_crop_harvest(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_bg1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_bg1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg1(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_bg2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_bg2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg2(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_humus_1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_humus_2(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_decomp_humus_1_sum(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_decomp_humus_2_sum(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_into_humus_1_sum(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_into_humus_2_sum(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_dist_green_2_soil(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_dist_woods_2_soil(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: co2flux_fire_all_2_atm(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: soil_respiration(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_pot_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_act_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_flux_correction(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: GPP_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_greenwood_2_litter(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: root_exudates(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_green_2_herb(:)

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)

    ! Local variables
    CHARACTER(len=*), PARAMETER     :: routine = modname//':recalc_per_tile_vars'
    INTEGER                         :: ics, ice, iblk

    dsl4jsb_Real2D_onChunk ::  c_green_ta
    dsl4jsb_Real2D_onChunk ::  c_reserve_ta
    dsl4jsb_Real2D_onChunk ::  c_woods_ta
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_1_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_2_ta

    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_1_sum_ta
    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_2_sum_ta
    dsl4jsb_Real2D_onChunk ::  c_into_humus_1_sum_ta
    dsl4jsb_Real2D_onChunk ::  c_into_humus_2_sum_ta

    dsl4jsb_Real2D_onChunk ::  cflux_dist_green_2_soil_ta
    dsl4jsb_Real2D_onChunk ::  cflux_dist_woods_2_soil_ta
    dsl4jsb_Real2D_onChunk ::  co2flux_fire_all_2_atm_ta
    dsl4jsb_Real2D_onChunk ::  soil_respiration_ta
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction_ta
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter_ta
    dsl4jsb_Real2D_onChunk ::  root_exudates_ta
    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb_ta
    dsl4jsb_Real2D_onChunk ::  co2flux_herb_2_atm_ta

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction
    dsl4jsb_Real2D_onChunk ::  fract_fpc_max

    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green_ta )               ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods_ta )               ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest_ta)         ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1_ta)            ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1_ta)         ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1_ta)      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1_ta)            ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1_ta)         ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1_ta)      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1_ta)             ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2_ta)          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2_ta)     ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2_ta)          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2_ta)     ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2_ta)            ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum_ta)    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum_ta)    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum_ta)      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum_ta)      ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_green_2_soil_ta )   ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta )   ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_fire_all_2_atm_ta )    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean_ta )              ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter_ta )! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb_ta )      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_herb_2_atm_ta )        ! opt out

    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )       ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )              ! in

    IF (PRESENT(c_green)) CALL get_per_tile(c_green_ta(:), c_green(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_woods)) CALL get_per_tile(c_woods_ta(:), c_woods(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_reserve)) CALL get_per_tile(c_reserve_ta(:), c_reserve(:), veg_fract_correction(:), fract_fpc_max(:))

    IF (PRESENT(c_crop_harvest)) CALL get_per_tile(c_crop_harvest_ta(:), &
      & c_crop_harvest(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_ag1)) CALL get_per_tile(c_acid_ag1_ta(:), &
      & c_acid_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_ag1)) CALL get_per_tile(c_water_ag1_ta(:), &
      & c_water_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_ag1)) CALL get_per_tile(c_ethanol_ag1_ta(:), &
      & c_ethanol_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_ag1)) CALL get_per_tile(c_nonsoluble_ag1_ta(:), &
      & c_nonsoluble_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_bg1)) CALL get_per_tile(c_acid_bg1_ta(:), &
      & c_acid_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_bg1)) CALL get_per_tile(c_water_bg1_ta(:), &
      & c_water_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_bg1)) CALL get_per_tile(c_ethanol_bg1_ta(:), &
      & c_ethanol_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_bg1)) CALL get_per_tile(c_nonsoluble_bg1_ta(:), &
      & c_nonsoluble_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_humus_1)) CALL get_per_tile(c_humus_1_ta(:), &
      & c_humus_1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_ag2)) CALL get_per_tile(c_acid_ag2_ta(:), &
      & c_acid_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_ag2)) CALL get_per_tile(c_water_ag2_ta(:), &
      & c_water_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_ag2)) CALL get_per_tile(c_ethanol_ag2_ta(:), &
      & c_ethanol_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_ag2)) CALL get_per_tile(c_nonsoluble_ag2_ta(:), &
      & c_nonsoluble_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_bg2)) CALL get_per_tile(c_acid_bg2_ta(:), &
      & c_acid_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_bg2)) CALL get_per_tile(c_water_bg2_ta(:), &
      & c_water_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_bg2)) CALL get_per_tile(c_ethanol_bg2_ta(:), &
      & c_ethanol_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_bg2)) CALL get_per_tile(c_nonsoluble_bg2_ta(:), &
      & c_nonsoluble_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_humus_2)) CALL get_per_tile(c_humus_2_ta(:), &
      & c_humus_2(:), veg_fract_correction(:), fract_fpc_max(:))

    IF (PRESENT(c_decomp_humus_1_sum)) CALL get_per_tile(c_decomp_humus_1_sum_ta(:), &
      & c_decomp_humus_1_sum(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_decomp_humus_2_sum)) CALL get_per_tile(c_decomp_humus_2_sum_ta(:), &
      & c_decomp_humus_2_sum(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_into_humus_1_sum)) CALL get_per_tile(c_into_humus_1_sum_ta(:), &
      & c_into_humus_1_sum(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_into_humus_2_sum)) CALL get_per_tile(c_into_humus_2_sum_ta(:), &
      & c_into_humus_2_sum(:), veg_fract_correction(:), fract_fpc_max(:))

    ! C fluxes
    IF (PRESENT(cflux_dist_green_2_soil)) CALL get_per_tile(cflux_dist_green_2_soil_ta(:), &
      & cflux_dist_green_2_soil(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_dist_woods_2_soil)) CALL get_per_tile(cflux_dist_woods_2_soil_ta(:), &
      & cflux_dist_woods_2_soil(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(co2flux_fire_all_2_atm))CALL get_per_tile(co2flux_fire_all_2_atm_ta(:), &
      & co2flux_fire_all_2_atm(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_pot_yDayMean)) CALL get_per_tile(NPP_pot_yDayMean_ta(:), &
      & NPP_pot_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_act_yDayMean)) CALL get_per_tile(NPP_act_yDayMean_ta(:), &
      & NPP_act_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(soil_respiration)) CALL get_per_tile(soil_respiration_ta(:), &
      & soil_respiration(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_flux_correction)) CALL get_per_tile(NPP_flux_correction_ta(:), &
      & NPP_flux_correction(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(GPP_yDayMean)) CALL get_per_tile(GPP_yDayMean_ta(:), &
      & GPP_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_c_greenwood_2_litter)) CALL get_per_tile(cflux_c_greenwood_2_litter_ta(:), &
      & cflux_c_greenwood_2_litter(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(root_exudates)) CALL get_per_tile(root_exudates_ta(:), &
      & root_exudates(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_c_green_2_herb)) CALL get_per_tile(cflux_c_green_2_herb_ta(:), &
      & cflux_c_green_2_herb(:), veg_fract_correction(:), fract_fpc_max(:))

  END SUBROUTINE recalc_per_tile_vars

  ! ================================================================================================================================
  !>
  !> calculates the current sum of all carbon state variables -> on ta! (per canopy vars are not aggregated!)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE calculate_current_c_ta_state_sum(tile, options, current_c_state_sum_ta, &
    & c_sum_veg_ta, c_sum_litter_ag_ta, c_sum_litter_bg_ta, c_sum_humus_ta)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(wp),                   INTENT(out)        :: current_c_state_sum_ta(:)
    REAL(wp), OPTIONAL,         INTENT(out)        :: c_sum_veg_ta(:)
    REAL(wp), OPTIONAL,         INTENT(out)        :: c_sum_litter_ag_ta(:)
    REAL(wp), OPTIONAL,         INTENT(out)        :: c_sum_litter_bg_ta(:)
    REAL(wp), OPTIONAL,         INTENT(out)        :: c_sum_humus_ta(:)

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    INTEGER :: nc, ic, ics, ice, iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':calculate_current_c_ta_state_sum'

    dsl4jsb_Real2D_onChunk ::  c_green_ta
    dsl4jsb_Real2D_onChunk ::  c_reserve_ta
    dsl4jsb_Real2D_onChunk ::  c_woods_ta
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_1_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_2_ta

    dsl4jsb_Get_memory(CARBON_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green_ta )               ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve_ta )             ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods_ta )               ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest_ta)         ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1_ta)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1_ta)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1_ta)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1_ta)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1_ta)             ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2_ta)          ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2_ta)        ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2_ta)     ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2_ta)          ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2_ta)        ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2_ta)     ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2_ta)            ! in

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      current_c_state_sum_ta(ic) = c_green_ta(ic) + c_woods_ta(ic) + c_reserve_ta(ic) + c_crop_harvest_ta(ic) + &
        & c_acid_ag1_ta(ic) + c_water_ag1_ta(ic) + c_ethanol_ag1_ta(ic) + c_nonsoluble_ag1_ta(ic) + &
        & c_acid_ag2_ta(ic) + c_water_ag2_ta(ic) + c_ethanol_ag2_ta(ic) + c_nonsoluble_ag2_ta(ic) + &
        & c_acid_bg1_ta(ic) + c_water_bg1_ta(ic) + c_ethanol_bg1_ta(ic) + c_nonsoluble_bg1_ta(ic) + &
        & c_acid_bg2_ta(ic) + c_water_bg2_ta(ic) + c_ethanol_bg2_ta(ic) + c_nonsoluble_bg2_ta(ic) + &
        & c_humus_1_ta(ic)  + c_humus_2_ta(ic)
    END DO
    !$ACC END PARALLEL LOOP

    IF (PRESENT(c_sum_veg_ta)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        c_sum_veg_ta(ic) = c_green_ta(ic) + c_woods_ta(ic) + c_reserve_ta(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (PRESENT(c_sum_litter_ag_ta)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        c_sum_litter_ag_ta(ic) = &
          & c_acid_ag1_ta(ic) + c_water_ag1_ta(ic) + c_ethanol_ag1_ta(ic) + c_nonsoluble_ag1_ta(ic) + &
          & c_acid_ag2_ta(ic) + c_water_ag2_ta(ic) + c_ethanol_ag2_ta(ic) + c_nonsoluble_ag2_ta(ic)
      END DO
      !$ACC END PARALLEL LOOP
    ENDIF

    IF (PRESENT(c_sum_litter_bg_ta)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        c_sum_litter_bg_ta(ic) = &
          & c_acid_bg1_ta(ic) + c_water_bg1_ta(ic) + c_ethanol_bg1_ta(ic) + c_nonsoluble_bg1_ta(ic) + &
          & c_acid_bg2_ta(ic) + c_water_bg2_ta(ic) + c_ethanol_bg2_ta(ic) + c_nonsoluble_bg2_ta(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (PRESENT(c_sum_humus_ta)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        c_sum_humus_ta(ic) = c_humus_1_ta(ic)  + c_humus_2_ta(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

  END SUBROUTINE calculate_current_c_ta_state_sum


  ! ================================================================================================================================
  !>
  !> calculates the current sum of all bg carbon state variables (per canopy)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE calculate_current_c_ag_1_and_bg_sums(tile, options)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':calculate_current_c_ag_1_and_bg_sums'
    INTEGER                     :: nc, ic, ics, ice, iblk

    dsl4jsb_Real2D_onChunk ::  c_bg_sum
    dsl4jsb_Real2D_onChunk ::  c_ag_sum_1

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1
    dsl4jsb_Real2D_onChunk ::  c_water_ag1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1

    dsl4jsb_Real2D_onChunk ::  c_acid_bg1
    dsl4jsb_Real2D_onChunk ::  c_water_bg1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk ::  c_humus_1

    dsl4jsb_Real2D_onChunk ::  c_acid_bg2
    dsl4jsb_Real2D_onChunk ::  c_water_bg2
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2
    dsl4jsb_Real2D_onChunk ::  c_humus_2

    dsl4jsb_Get_memory(CARBON_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_bg_sum )             ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ag_sum_1 )           ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)             ! in

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      c_bg_sum(ic) = c_acid_bg1(ic) + c_water_bg1(ic) + c_ethanol_bg1(ic) + c_nonsoluble_bg1(ic) + &
              & c_acid_bg2(ic) + c_water_bg2(ic) + c_ethanol_bg2(ic) + c_nonsoluble_bg2(ic) + &
              & c_humus_1(ic)  + c_humus_2(ic)

      c_ag_sum_1(ic) =  c_acid_ag1(ic) + c_water_ag1(ic) + c_ethanol_ag1(ic) + c_nonsoluble_ag1(ic) + c_humus_1(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calculate_current_c_ag_1_and_bg_sums


  ! ====================================================================================================== !
  !
  !> Relocates content of (expected) own active to own passive vars that are part of lcc_relocations
  !
  SUBROUTINE carbon_transfer_from_active_to_passive_vars_onChunk(lcc_relocations, tile, i_tile, options)

    USE mo_util,                ONLY: one_of
    USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc, collect_matter_of_active_vars_onChunk
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_carbon_constants,    ONLY: fract_green_aboveGround, fract_wood_aboveGround,           &
      &                               nr_of_yasso_pools, carbon_required_passive_vars, coefficient_ind_of_passive_vars
    USE mo_carbon_process,      ONLY: add_litter_to_yasso_pool
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing matter from an active to a passive var
    CLASS(t_jsb_tile_abstract), INTENT(IN)    :: tile    !< tile for which the relocation is be conducted
    INTEGER,                    INTENT(IN)    :: i_tile  !< index of the tile in lcc structure
    TYPE(t_jsb_task_options),   INTENT(IN)    :: options !< run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_transfer_from_active_to_passive_vars_onChunk'

    TYPE(t_jsb_model), POINTER   :: model

    INTEGER :: i, k, i_pool, i_passive, i_litter_coefficient, nc, ics, ice, iblk

    REAL(wp) :: litter_coefficient, fraction, this_fract
    REAL(wp) :: this_fract_green_aboveGround, this_fract_wood_aboveGround, fire_fract_wood_2_atmos, not_burned_AG_fract
    REAL(wp), DIMENSION(options%nc) :: non_woody_litter, woody_litter, this_litter

    dsl4jsb_Def_config(CARBON_)
    ! -------------------------------------------------------------------------------------------------- !

    !JN-TODO: is it less expensive to only do it for tiles with actual area changes?
    nc = options%nc
    ics = options%ics
    ice = options%ice
    iblk = options%iblk

    IF (debug_on() .AND. options%iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(CARBON_)

    ! Collect matter from potential active variables -- here unfortunately explicit
    non_woody_litter = 0.0
    woody_litter = 0.0
    CALL collect_matter_of_active_vars_onChunk( non_woody_litter, lcc_relocations, &
      & i_tile, iblk, ics, ice, [character(len=VARNAME_LEN) :: 'c_green', 'c_reserve' ])
    CALL collect_matter_of_active_vars_onChunk( woody_litter, lcc_relocations,     &
      & i_tile, iblk, ics, ice, [character(len=VARNAME_LEN) ::  'c_woods' ])

    ! Usually green, reserve and wood are distributed to ag and bg pools according to given fractions
    this_fract_green_aboveGround = fract_green_aboveGround
    this_fract_wood_aboveGround = fract_wood_aboveGround
    IF(lcc_relocations%name .EQ. 'flcc_lcc') THEN
      ! in this case there is no transport of green or reserve to AG litter pools, since it burned
      this_fract_green_aboveGround = 0.0_wp
      ! # in this case part of the ag wood has also already been added to the atmosphere flux
      fire_fract_wood_2_atmos = dsl4jsb_Config(CARBON_)%fire_fract_wood_2_atmos
      not_burned_AG_fract = (1.0_wp - fire_fract_wood_2_atmos) * fract_wood_aboveGround
      ! (not_burned_AG_fract + (1.0_wp - fract_wood_aboveGround)) = Rest of wood after burning parts of AG wood
      this_fract_wood_aboveGround = not_burned_AG_fract / (not_burned_AG_fract + (1.0_wp - fract_wood_aboveGround))
    ENDIF

    !--- And distribute them according to the yasso scheme
    DO i_pool = 1, nr_of_yasso_pools
      i_passive = one_of(carbon_required_passive_vars(i_pool), lcc_relocations%passive_vars_names)
      i_litter_coefficient = coefficient_ind_of_passive_vars(i_pool)

      IF (INDEX(carbon_required_passive_vars(i_pool), '1') > 0) THEN
        this_litter = non_woody_litter
        this_fract = this_fract_green_aboveGround
        litter_coefficient = dsl4jsb_Lctlib_param(LeafLit_coef(i_litter_coefficient))
      ELSEIF (INDEX(carbon_required_passive_vars(i_pool), '2') > 0) THEN
        this_litter = woody_litter
        this_fract = this_fract_wood_aboveGround
        litter_coefficient = dsl4jsb_Lctlib_param(WoodLit_coef(i_litter_coefficient))
      ELSE
        CALL finish(TRIM(routine), 'Violation of assertion: ' // carbon_required_passive_vars(i_pool) &
          & // ' did not contain "1" or "2" which are key to identifiy woody vs non-woody yasso pools. Please check!')
      ENDIF

      IF (INDEX(carbon_required_passive_vars(i_pool), 'ag') > 0) THEN
        fraction = this_fract
      ELSEIF (INDEX(carbon_required_passive_vars(i_pool), 'bg') > 0) THEN
        fraction = 1.0_wp - this_fract
      ELSEIF (INDEX(carbon_required_passive_vars(i_pool), 'humus') > 0) THEN
        fraction = 1.0_wp
      ELSE
        CALL finish(TRIM(routine), 'Violation of assertion: ' // carbon_required_passive_vars(i_pool) &
          & // ' did not contain "ag" or "bg" which are key to identifiy above vs below ground yasso pools. Please check!')
      ENDIF

      !JN-TODO: better way? Reiner? DSL?
      DO i = 1,nc
        k = ics + i - 1
        CALL add_litter_to_yasso_pool( &
          & lcc_relocations%passive_vars(i_passive)%p%relocate_this(k,i_tile,iblk), &
          & this_litter(i), fraction, litter_coefficient)
      ENDDO
    ENDDO

  END SUBROUTINE carbon_transfer_from_active_to_passive_vars_onChunk


  ! ====================================================================================================== !
  !
  !> Rescales those carbon variables which need to be merged upon ageing (FAGE process), i.e.
  !> when the fraction with the maximum age of one age class is pushed into the next age class
  !
  SUBROUTINE rescale_carbon_vars_upon_ageing_induced_ref_area_change(tile, options, old_veg_fract_correction)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract),   INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),        INTENT(in)         :: options
    REAL(wp), DIMENSION(options%nc), INTENT(in)         :: old_veg_fract_correction
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk, ics, ice

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction

    dsl4jsb_Real2D_onChunk ::  GPP_sum
    dsl4jsb_Real2D_onChunk ::  NPP_pot_sum
    dsl4jsb_Real2D_onChunk ::  LAI_yDayMean
    dsl4jsb_Real2D_onChunk ::  LAI_sum
    dsl4jsb_Real2D_onChunk ::  current_max_green
    dsl4jsb_Real2D_onChunk ::  veg_carbon_at_max_green
    ! -------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(CARBON_, GPP_sum)
    dsl4jsb_Get_var2D_onChunk(CARBON_, NPP_pot_sum)
    dsl4jsb_Get_var2D_onChunk(CARBON_, LAI_yDayMean)
    dsl4jsb_Get_var2D_onChunk(CARBON_, LAI_sum)
    dsl4jsb_Get_var2D_onChunk(CARBON_, current_max_green)
    dsl4jsb_Get_var2D_onChunk(CARBON_, veg_carbon_at_max_green)

    dsl4jsb_Get_var2D_onChunk(PHENO_, veg_fract_correction)

    GPP_sum(:)                 = GPP_sum(:) * (old_veg_fract_correction(:) / veg_fract_correction(:))
    NPP_pot_sum(:)             = NPP_pot_sum(:) * (old_veg_fract_correction(:) / veg_fract_correction(:))
    LAI_yDayMean(:)            = LAI_yDayMean(:) * (old_veg_fract_correction(:) / veg_fract_correction(:))
    LAI_sum(:)                 = LAI_sum(:) * (old_veg_fract_correction(:) / veg_fract_correction(:))
    current_max_green(:)       = current_max_green(:) * (old_veg_fract_correction(:) / veg_fract_correction(:))
    veg_carbon_at_max_green(:) = veg_carbon_at_max_green(:) * (old_veg_fract_correction(:) / veg_fract_correction(:))

  END SUBROUTINE rescale_carbon_vars_upon_ageing_induced_ref_area_change


  ! ====================================================================================================== !
  !
  !> Induces merging for those carbon variables which need to be merged upon ageing (FAGE process)
  !
  SUBROUTINE merge_carbon_vars_upon_ageing(target, source, options, moved_area)

    USE mo_fage_process,     ONLY : weighted_avg_per_canopy_var_upon_area_movement
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract),      INTENT(inout) :: target     !< target age class
    CLASS(t_jsb_tile_abstract),      INTENT(inout) :: source     !< source age class
    TYPE(t_jsb_task_options),        INTENT(in)    :: options    !< Additional run-time parameters
    REAL(wp), DIMENSION(options%nc), INTENT(in)    :: moved_area !< area moved from source to target ac
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk, ics, ice, nc

    dsl4jsb_Def_memory_tile(PHENO_, target)
    dsl4jsb_Def_memory_tile(PHENO_, source)

    dsl4jsb_Def_memory_tile(CARBON_, target)
    dsl4jsb_Def_memory_tile(CARBON_, source)

    dsl4jsb_Real2D_onChunk :: veg_fract_correction_source
    dsl4jsb_Real2D_onChunk :: veg_fract_correction_target

    dsl4jsb_Real2D_onChunk :: GPP_sum_source
    dsl4jsb_Real2D_onChunk :: GPP_sum_target
    dsl4jsb_Real2D_onChunk :: NPP_pot_sum_source
    dsl4jsb_Real2D_onChunk :: NPP_pot_sum_target
    dsl4jsb_Real2D_onChunk :: LAI_sum_source
    dsl4jsb_Real2D_onChunk :: LAI_sum_target
    dsl4jsb_Real2D_onChunk :: LAI_yDayMean_source
    dsl4jsb_Real2D_onChunk :: LAI_yDayMean_target
    dsl4jsb_Real2D_onChunk :: current_max_green_source
    dsl4jsb_Real2D_onChunk :: current_max_green_target
    dsl4jsb_Real2D_onChunk :: veg_carbon_at_max_green_source
    dsl4jsb_Real2D_onChunk :: veg_carbon_at_max_green_target

    REAL(wp), DIMENSION(options%nc) :: target_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':merge_carbon_vars_upon_ageing'
    ! -------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For target '//TRIM(target%name)//' ...')

    dsl4jsb_Get_memory_tile(PHENO_, source)
    dsl4jsb_Get_memory_tile(PHENO_, target)

    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, target)

    dsl4jsb_Get_memory_tile(CARBON_, source)
    dsl4jsb_Get_memory_tile(CARBON_, target)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, GPP_sum, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, GPP_sum, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, NPP_pot_sum, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, NPP_pot_sum, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_yDayMean, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_yDayMean, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_sum, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_sum, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, current_max_green, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, current_max_green, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, veg_carbon_at_max_green, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, veg_carbon_at_max_green, target)

    CALL target%Get_fraction(ics, ice, iblk, fract=target_fract(:))

    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  GPP_sum_source(:), GPP_sum_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  NPP_pot_sum_source(:), NPP_pot_sum_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  LAI_sum_source(:), LAI_sum_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  LAI_yDayMean_source(:), LAI_yDayMean_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  current_max_green_source(:), current_max_green_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  veg_carbon_at_max_green_source(:), veg_carbon_at_max_green_target(:))

  END SUBROUTINE merge_carbon_vars_upon_ageing


  ! ================================================================================================================================
  !>
  !> conducts a carbon conservation test on this tile's current states and given old sum + fluxes
  !!        -> on ta! (per canopy vars are not aggregated!)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE check_carbon_conservation(tile, options, old_c_state_sum_ta, cflux_ta, Cconserve)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(wp),                   INTENT(in)         :: old_c_state_sum_ta(:)
    REAL(wp),                   INTENT(in)         :: cflux_ta(:)
        ! Note: cflux is not really a flux, but already multiplied with time unit!
    REAL(wp),                   INTENT(inout)      :: Cconserve(:)

    !dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    INTEGER :: ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':check_carbon_conservation'

    CALL calculate_current_c_ta_state_sum(tile, options, Cconserve)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, SIZE(cflux_ta)
      Cconserve(ic) = Cconserve(ic) - old_c_state_sum_ta(ic) - cflux_ta(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE check_carbon_conservation

  ! ================================================================================================================================
  !>
  !> conducts a box carbon conservation test for the last day, called each new day -- (jsbach_start_timestep)
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE yday_carbon_conservation_test(tile)

    USE mo_jsb_grid,           ONLY: Get_grid

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(L2A_)
    dsl4jsb_Real2D_onChunk ::  carbon_conservation_test, yday_c_state_sum

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid
    TYPE(t_jsb_task_options)        :: options

    REAL(wp), ALLOCATABLE :: C_flux_yday(:)
    INTEGER               :: ics, ice, iblk, nc

    dsl4jsb_Get_memory(L2A_)
    dsl4jsb_Get_memory(CARBON_)

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Note: options is defined here as a local variable. It is needed in calculate_current_c_ta_state_sum
    options%ics = 1
    options%ice = grid%nproma

    ics = options%ics
    ice = options%ice

    ALLOCATE(C_flux_yday(grid%nproma))

    DO iblk = 1, grid%nblks

      options%iblk = iblk
      IF (iblk == grid%nblks) THEN
        options%ice = grid%npromz
        ice = grid%npromz
        ! TODO: ics and ice should be retrieved for each iblk (are not necessarily constant for each block in ICON)
      END IF
      options%nc = ice - ics + 1
      nc = options%nc

      dsl4jsb_Get_var2D_onChunk(L2A_,  carbon_conservation_test )     ! inout
      dsl4jsb_Get_var2D_onChunk(L2A_,  yday_c_state_sum )             ! inout

      ! Do not test as long as carbon_conservation_test still carries the initialisation value
      IF (.NOT. ALL(carbon_conservation_test .EQ. -999.0_wp)) THEN

        ! CO2 flux due to npp, soil respiration and herbivory
        C_flux_yday(1:nc) = &
          & - ((  dsl4jsb_var2D_onChunk(CARBON_, co2flux_npp_2_atm_yday_ta)  &
          &     + dsl4jsb_var2D_onChunk(CARBON_, co2flux_soilresp_2_atm_ta)  &
          &     + dsl4jsb_var2D_onChunk(CARBON_, co2flux_herb_2_atm_ta)      &
          &   ) * sec_per_day / molarMassCO2_kg)

        ! CO2 flux due to fire
        IF (model%processes(DISTURB_)%p%config%active) THEN
          C_flux_yday(1:nc) =  C_flux_yday(1:nc) &
            & - dsl4jsb_var2D_onChunk(CARBON_, co2flux_fire_all_2_atm_ta) * sec_per_day / molarMassCO2_kg
        ENDIF

        CALL calculate_current_c_ta_state_sum(tile, options, carbon_conservation_test)
        carbon_conservation_test = carbon_conservation_test - yday_c_state_sum  - C_flux_yday(1:nc)

      ELSE
        carbon_conservation_test = 0._wp
      ENDIF

      CALL calculate_current_c_ta_state_sum(tile, options, yday_c_state_sum)

    END DO

    DEALLOCATE(C_flux_yday)

  END SUBROUTINE yday_carbon_conservation_test

  ! ================================================================================================================================
  !>
  !> calculations of diagnostic global sums
  !!        called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !! @param[in,out] tile    Tile for which routine is executed.
  !
  SUBROUTINE global_carbon_diagnostics(tile)

#ifndef __ICON__
    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_carbon_diagnostics'
    IF (debug_on()) CALL message(TRIM(routine), 'Global diagnostics only available with ICON')
#else

    USE mo_carbon_constants,      ONLY: molarMassC_kg, sec_per_year
    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_memory(CARBON_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_carbon_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: C_sum_veg_ta
    dsl4jsb_Real2D_onDomain :: C_sum_litter_ag_ta
    dsl4jsb_Real2D_onDomain :: C_sum_litter_bg_ta
    dsl4jsb_Real2D_onDomain :: C_sum_humus_ta
    dsl4jsb_Real2D_onDomain :: C_sum_natural_ta
    dsl4jsb_Real2D_onDomain :: NPP_act_yDayMean_ta
    dsl4jsb_Real2D_onDomain :: GPP_yDayMean_ta
    dsl4jsb_Real2D_onDomain :: soil_respiration_ta

    REAL(wp), POINTER       :: C_sum_veg_gsum(:)
    REAL(wp), POINTER       :: C_sum_litter_ag_gsum(:)
    REAL(wp), POINTER       :: C_sum_litter_bg_gsum(:)
    REAL(wp), POINTER       :: C_sum_humus_gsum(:)
    REAL(wp), POINTER       :: C_sum_natural_gsum(:)
    REAL(wp), POINTER       :: NPP_act_yDayMean_gsum(:)
    REAL(wp), POINTER       :: GPP_yDayMean_gsum(:)
    REAL(wp), POINTER       :: soil_respiration_gsum(:)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)


    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_veg_ta)                ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_litter_ag_ta)          ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_litter_bg_ta)          ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_humus_ta)              ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_natural_ta)            ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  NPP_act_yDayMean_ta)         ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  GPP_yDayMean_ta)             ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  soil_respiration_ta)         ! in

    C_sum_veg_gsum        => CARBON__mem%C_sum_veg_gsum%ptr(:)        ! out
    C_sum_litter_ag_gsum  => CARBON__mem%C_sum_litter_ag_gsum%ptr(:)  ! out
    C_sum_litter_bg_gsum  => CARBON__mem%C_sum_litter_bg_gsum%ptr(:)  ! out
    C_sum_humus_gsum      => CARBON__mem%C_sum_humus_gsum%ptr(:)      ! out
    C_sum_natural_gsum    => CARBON__mem%C_sum_natural_gsum%ptr(:)    ! out
    NPP_act_yDayMean_gsum => CARBON__mem%NPP_act_yDayMean_gsum%ptr(:) ! out
    GPP_yDayMean_gsum     => CARBON__mem%GPP_yDayMean_gsum%ptr(:)     ! out
    soil_respiration_gsum => CARBON__mem%soil_respiration_gsum%ptr(:) ! out

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea


    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Domain Mask - to mask all halo cells for global sums (otherwise these
    ! cells are counted twice)
    ALLOCATE (in_domain(grid%nproma,grid%nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ALLOCATE (scaling(grid%nproma,grid%nblks))

    ! Calculate global carbon inventories, if requested for output
    !  => Conversion from [mol(C)/m^2] to [GtC]
    !     1 mol C = molarMassC_kg kg C   => 1 mol C = molarMassC_kg * e-12 Gt C
    scaling(:,:) = molarMassC_kg * 1.e-12_wp * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (CARBON__mem%C_sum_veg_gsum%is_in_output)        &
      &  c_sum_veg_gsum        = global_sum_array(C_sum_veg_ta(:,:)        * scaling(:,:))
    IF (CARBON__mem%C_sum_litter_ag_gsum%is_in_output)  &
      &  c_sum_litter_ag_gsum  = global_sum_array(C_sum_litter_ag_ta(:,:)  * scaling(:,:))
    IF (CARBON__mem%C_sum_litter_bg_gsum%is_in_output)  &
      &  c_sum_litter_bg_gsum  = global_sum_array(C_sum_litter_bg_ta(:,:)  * scaling(:,:))
    IF (CARBON__mem%C_sum_humus_gsum%is_in_output)      &
      &  c_sum_humus_gsum      = global_sum_array(C_sum_humus_ta(:,:)      * scaling(:,:))
    IF (CARBON__mem%C_sum_natural_gsum%is_in_output)    &
      &  c_sum_natural_gsum    = global_sum_array(C_sum_natural_ta(:,:)    * scaling(:,:))
    IF (CARBON__mem%NPP_act_yDayMean_gsum%is_in_output) &
      &  NPP_act_yDayMean_gsum = global_sum_array(NPP_act_yDayMean_ta(:,:) * scaling(:,:) * sec_per_year)
    IF (CARBON__mem%GPP_yDayMean_gsum%is_in_output)     &
      &  GPP_yDayMean_gsum     = global_sum_array(GPP_yDayMean_ta(:,:)     * scaling(:,:) * sec_per_year)
    IF (CARBON__mem%soil_respiration_gsum%is_in_output) &
      &  soil_respiration_gsum = global_sum_array(soil_respiration_ta(:,:) * scaling(:,:) * sec_per_year)

    DEALLOCATE (scaling, in_domain)

#endif
  END SUBROUTINE global_carbon_diagnostics

#endif
END MODULE mo_carbon_interface
