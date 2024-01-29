!> Contains the interfaces to the disturb process
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

MODULE mo_disturb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_carbon_constants,   ONLY: molarMassCO2_kg, sec_per_day

  ! Use of processes in this module (Get_disturb_memory and Get_disturb_config)
  dsl4jsb_Use_processes DISTURB_, FUEL_, ASSIMI_, A2L_, CARBON_, SEB_

  ! Use of process configurations (t_disturb_config)
  dsl4jsb_Use_config(DISTURB_)
  dsl4jsb_Use_config(CARBON_)

  ! Use of process memories (t_disturb_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(DISTURB_)
  dsl4jsb_Use_memory(FUEL_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(SEB_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  ! R: Disturbance sollte nicht zu Dynveg, da sie auch stattfinden soll, wenn dynveg aus!
  !    Disturbance sollte nicht zu CARBON, da er nur zum Teil mit C zu tun hat.
  !    So verändert er über die disturbed fractions wenn dynveg an ist evtl. die
  !    cover fractions der pfts. Disturbance ist eigentlich ein Vegetations-Prozess neben Dynveg.
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_disturb_tasks !,t_disturb_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_interface'

  !> Type definition for natural_disturbances task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_natural_disturbances
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_natural_disturbances    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_natural_disturbances !< Aggregates computed task variables
  END TYPE tsk_natural_disturbances

  !> Constructor interface for natural_disturbances task
  INTERFACE tsk_natural_disturbances
    PROCEDURE Create_task_natural_disturbances     !< Constructor function for task
  END INTERFACE tsk_natural_disturbances

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for C_natural_disturbance task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "natural_disturbances"
  !!
  FUNCTION Create_task_natural_disturbances(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_natural_disturbances::return_ptr)
    CALL return_ptr%Construct(name='natural_disturbances', process_id=DISTURB_, owner_model_id=model_id)

  END FUNCTION Create_task_natural_disturbances

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for disturb process
  !!
  !! @param[in,out] this      Instance of disturb process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_disturb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_natural_disturbances(model_id))

  END SUBROUTINE Register_disturb_tasks


  ! ================================================================================================================================
  !>
  !> Implementation of task natural_disturbances
  !! Task "natural_disturbances" calculates the disturbance of plants
  !! due to fire and windbreak.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_natural_disturbances(tile, options)

    USE mo_disturb_process,       ONLY: burned_fract_jsbach, broken_woody_fract_jsbach, get_relative_humidity_air
    USE mo_carbon_process,        ONLY: relocate_carbon_fire, relocate_carbon_damage
    USE mo_carbon_interface,      ONLY: recalc_per_tile_vars, calculate_current_c_ag_1_and_bg_sums, &
                                      & calculate_current_c_ta_state_sum, check_carbon_conservation
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_time,              ONLY: is_time_experiment_start, timesteps_per_day, is_newday
    USE mo_disturbance_constants, ONLY: persist_rel_hum, persist_wind_10m

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model
    REAL(wp), DIMENSION(options%nc)           :: old_c_state_sum_ta, dummy_flux

    INTEGER :: iblk, nc, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_natural_disturbances'

    ! Declare process memories
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_memory(FUEL_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_config(CARBON_)

    ! Declare pointers to variables in memory
    !INTEGER :: ncanopy                 ! number of canopy layers

    !dsl4jsb_Real2D_onChunk :: C4flag ! Photosynthetic pathway (C3: 0; C4: 1)
                                      ! R: nur falls man C4flag doch als REAL haben wollte


!     dsl4jsb_Real2D_onChunk :: CO2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    dsl4jsb_Real2D_onChunk :: fuel                   ! Amount of fuel on which the fire was based.
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: wind_10m               ! 10m wind speed [m/s]
    dsl4jsb_Real2D_onChunk :: max_wind_10m_act
    dsl4jsb_Real2D_onChunk :: prev_day_max_wind_10m
    dsl4jsb_Real2D_onChunk :: max_wind_10m
    dsl4jsb_Real2D_onChunk :: q_rel_air_climbuf
    dsl4jsb_Real2D_onChunk :: q_rel_air_climbuf_yDay

    dsl4jsb_Real2D_onChunk :: cconservation_fire
    dsl4jsb_Real2D_onChunk :: cconservation_wind

    dsl4jsb_Real2D_onChunk :: burned_fract
    dsl4jsb_Real2D_onChunk :: damaged_fract
    !
    dsl4jsb_Real2D_onChunk :: c_green
    dsl4jsb_Real2D_onChunk :: c_woods
    dsl4jsb_Real2D_onChunk :: c_reserve
    dsl4jsb_Real2D_onChunk :: cflux_dist_greenreserve_2_soil
    dsl4jsb_Real2D_onChunk :: cflux_dist_woods_2_soil
    !
    dsl4jsb_Real2D_onChunk :: c_acid_ag1
    dsl4jsb_Real2D_onChunk :: c_water_ag1
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag1
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag1
    dsl4jsb_Real2D_onChunk :: c_acid_ag2
    dsl4jsb_Real2D_onChunk :: c_water_ag2
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag2
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag2
    ! below ground
    dsl4jsb_Real2D_onChunk :: c_acid_bg1
    dsl4jsb_Real2D_onChunk :: c_water_bg1
    dsl4jsb_Real2D_onChunk :: c_ethanol_bg1
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk :: c_acid_bg2
    dsl4jsb_Real2D_onChunk :: c_water_bg2
    dsl4jsb_Real2D_onChunk :: c_ethanol_bg2
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_bg2
    !
    dsl4jsb_Real2D_onChunk :: c_humus_1
    dsl4jsb_Real2D_onChunk :: c_humus_2
    dsl4jsb_Real2D_onChunk :: cflux_fire_all_2_atm
    !
    dsl4jsb_Real2D_onChunk :: cover_fract_pot
    !
    dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm_ta

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & q_rel_air

    !LOGICAL  :: pft_mask ! True for pfts to which the disturbance should be applied.
    REAL(wp) :: persist ! Factor (as fraction) determining the persistance of q_rel_air_climbuf
                        ! relative to the time-step-actual q_rel_air
    REAL(wp) :: fire_minimum, fire_tau

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(DISTURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process configurations
    dsl4jsb_Get_config(DISTURB_)
    dsl4jsb_Get_config(CARBON_)

    ! Get process memories
    dsl4jsb_Get_memory(DISTURB_)
    dsl4jsb_Get_memory(FUEL_)
    dsl4jsb_Get_memory(CARBON_)
    !dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(FUEL_,      fuel)                     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       q_air)                    ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       press_srf)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       wind_10m)                 ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,       t_unfilt)                 ! in
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   max_wind_10m_act)         ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   prev_day_max_wind_10m)    ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   max_wind_10m)             ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_air_climbuf)        ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_air_climbuf_yDay)   ! out

    dsl4jsb_Get_var2D_onChunk(DISTURB_,   burned_fract)             ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   damaged_fract)            ! out
    !
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   cconservation_wind)     ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   cconservation_fire)     ! out
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_green)                  ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_woods)                  ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_reserve)                ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_greenreserve_2_soil)! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_woods_2_soil)  ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag1)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag1)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag1)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag2)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag2)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag2)        ! inout
    ! below ground
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_bg1)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_bg1)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_bg1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_bg1)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_bg2)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_bg2)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_bg2)        ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_humus_1)               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_humus_2)               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_fire_all_2_atm)    ! inout
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,    cover_fract_pot)         ! in

    !ta variable required for the carbon conservation test
    dsl4jsb_Get_var2D_onChunk(CARBON_,    co2flux_fire_all_2_atm_ta) ! out

    ! ---------------------------
    ! Go

    q_rel_air(:) = get_relative_humidity_air(nc, q_air(:), t_unfilt(:), press_srf(:))

    ! For each time step update the climate buffer

    ! R: This may be later in a climbuf process for dynveg and disturbance.
    ! Maximum daily wind speed and relative air humidity smoothed in time
    persist = persist_rel_hum ** (1._wp/REAL(timesteps_per_day(options%dtime),wp))
    q_rel_air_climbuf(:) = q_rel_air_climbuf(:) * persist  +  MIN(q_rel_air(:),100._wp) * (1._wp - persist)

    IF (.NOT. is_newday(options%current_datetime, options%dtime) .OR. is_time_experiment_start(options%current_datetime)) THEN
      WHERE (wind_10m(:) > max_wind_10m_act(:))
         max_wind_10m_act(:) = wind_10m(:)
      END WHERE
    ELSE
      prev_day_max_wind_10m(:) = max_wind_10m_act(:)
      max_wind_10m(:) = max_wind_10m(:) * persist_wind_10m + max_wind_10m_act(:) * (1._wp - persist_wind_10m)
      max_wind_10m_act(:) = wind_10m(:)
      q_rel_air_climbuf_yDay(:) = q_rel_air_climbuf(:)
    END IF

    IF (is_newday(options%current_datetime, options%dtime)) THEN ! only once per day
      cflux_dist_greenreserve_2_soil(:) = 0._wp
      cflux_dist_woods_2_soil(:) = 0._wp

      burned_fract(:) = 0._wp
      IF (      dsl4jsb_Lctlib_param(dynamic_PFT) &
          .OR. (dsl4jsb_Config(DISTURB_)%lburn_pasture .AND. dsl4jsb_Lctlib_param(pasture_PFT)) &
         ) THEN

        SELECT CASE (dsl4jsb_Config(DISTURB_)%fire_algorithm)
          CASE (0) !! No fire algorithm
          CASE (1) !! jsbach algorithm

            ! Calculate burned area. Woody and grass types call the subroutine with different parameters
            IF (dsl4jsb_Lctlib_param(woody_PFT)) THEN  ! Woody Type
              fire_minimum = dsl4jsb_Config(DISTURB_)%fire_minimum_woody
              fire_tau     = dsl4jsb_Config(DISTURB_)%fire_tau_woody
            ELSE                             ! all other types
              fire_minimum = dsl4jsb_Config(DISTURB_)%fire_minimum_grass
              fire_tau     = dsl4jsb_Config(DISTURB_)%fire_tau_grass
            ENDIF

            CALL burned_fract_jsbach(                            &
              & dsl4jsb_Config(DISTURB_)%fire_rel_hum_threshold, &
              & dsl4jsb_Config(DISTURB_)%fire_litter_threshold,  &
              & fire_minimum,                                    & ! in
              & fire_tau,                                        & ! in
              & q_rel_air_climbuf(:),                            & ! in
              & fuel(:),                                         & ! in
              & burned_fract(:)                                  & ! inout
              & )

          CASE (2) !! Arora & Boer algorithm
            CALL finish('disturbed_fract','Arora & Boer algorithm not implemented yet.')
          CASE (3)
            CALL finish('disturbed_fract','Thornicke spit fire not implemented yet.')
          CASE (4) !! Read burned_fract
          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown fire algorithm')
        END SELECT

      END IF

      damaged_fract(:) = 0._wp
      IF (dsl4jsb_Lctlib_param(woody_PFT) .AND. dsl4jsb_Lctlib_param(dynamic_PFT)) THEN

        SELECT CASE (dsl4jsb_Config(DISTURB_)%windbreak_algorithm)
          CASE (0) !! No windbreak algorithm
          CASE (1) !! jsbach algorithm
            CALL broken_woody_fract_jsbach(                 & ! in
              & dsl4jsb_Config(DISTURB_)%wnd_threshold,     & ! in
              & dsl4jsb_Config(DISTURB_)%wnd_damage_scale,  & ! in
              & cover_fract_pot(:),                         & ! in
              & prev_day_max_wind_10m(:),                   & ! in
              & max_wind_10m(:),                            & ! in
              & damaged_fract(:)                            & ! inout
              & )
          CASE DEFAULT
            CALL finish('disturbed_frac','Unknown windbreak algorithm')
        END SELECT

      END IF

      IF ( tile%Is_process_active(CARBON_) ) THEN

        ! Calculate carbon relocation connected with damages and fires
        ! R: hier sollen nur die Variablen übergeben werden, die von dieser Prozedur hier im Folgenden benötigt werden.
        !    Andere Variablen werden erst später bereit gestellt, denn in diesem Fall rufen wir ein weiteres
        !    Interface beim CARBON Prozess auf (SUBROUTINE relocate_disturbed_carbon).
        !    Erst dort sollen dann alle Variablen an die eigentlichen Prozess-Prozeduren übergeben werden.

        SELECT CASE (dsl4jsb_Config(DISTURB_)%windbreak_algorithm)
          CASE (0) !! No windbreak algorithm
          CASE (1) !! jsbach algorithm

            !Calculate current sum before operation for c conservation test
            CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

            CALL relocate_carbon_damage(                 &
              & damaged_fract(:),                        & ! in
              & dsl4jsb_Lctlib_param(LeafLit_coef(1:5)), & ! in
              & dsl4jsb_Lctlib_param(WoodLit_coef(1:5)), & ! in
              & c_green(:),                              & ! inout
              & c_woods(:),                              & ! inout
              & c_reserve(:),                            & ! inout
              & cflux_dist_greenreserve_2_soil(:),       & ! inout
              & cflux_dist_woods_2_soil(:),              & ! inout
              & c_acid_ag1(:),                           & ! inout
              & c_water_ag1(:),                          & ! inout
              & c_ethanol_ag1(:),                        & ! inout
              & c_nonsoluble_ag1(:),                     & ! inout
              & c_acid_ag2(:),                           & ! inout
              & c_water_ag2(:),                          & ! inout
              & c_ethanol_ag2(:),                        & ! inout
              & c_nonsoluble_ag2(:),                     & ! inout
              & c_acid_bg1(:),                           & ! inout
              & c_water_bg1(:),                          & ! inout
              & c_ethanol_bg1(:),                        & ! inout
              & c_nonsoluble_bg1(:),                     & ! inout
              & c_acid_bg2(:),                           & ! inout
              & c_water_bg2(:),                          & ! inout
              & c_ethanol_bg2(:),                        & ! inout
              & c_nonsoluble_bg2(:),                     & ! inout
              & c_humus_1(:),                            & ! inout
              & c_humus_2(:)                             & ! inout
              & )

            CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
            CALL recalc_per_tile_vars(tile, options,                                      &
              & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
              & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
              & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
              & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
              & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
              & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
              & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
              & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
              & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
              & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
              & cflux_dist_greenreserve_2_soil = cflux_dist_greenreserve_2_soil(:),       &
              & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:))

            ! wind has no extra flux to atmos, only redistribution among cpools - therefore cflux = 0._wp
            dummy_flux = 0._wp
            CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
              & dummy_flux(:), cconservation_wind(:))

          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown windbreak algorithm')
        END SELECT

        cflux_fire_all_2_atm = 0._wp
        SELECT CASE (dsl4jsb_Config(DISTURB_)%fire_algorithm)
          CASE (0) !! No fire algorithm
          CASE (1) !! jsbach algorithm

           !Calculate current sum before operation for c conservation test
            CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

            CALL relocate_carbon_fire(                    &
              & burned_fract(:),                          & ! in, only for windbreak
              & dsl4jsb_Config(CARBON_)%fire_fract_wood_2_atmos , & ! in
              & dsl4jsb_Lctlib_param(LeafLit_coef(1:5)),  & ! in
              & dsl4jsb_Lctlib_param(WoodLit_coef(1:5)),  & ! in
              & c_green(:),                               & ! inout
              & c_woods(:),                               & ! inout
              & c_reserve(:),                             & ! inout
              & cflux_dist_greenreserve_2_soil(:),        & ! inout
              & cflux_dist_woods_2_soil(:),               & ! inout
              & c_acid_ag1(:),                            & ! inout
              & c_water_ag1(:),                           & ! inout
              & c_ethanol_ag1(:),                         & ! inout
              & c_nonsoluble_ag1(:),                      & ! inout
              & c_acid_ag2(:),                            & ! inout
              & c_water_ag2(:),                           & ! inout
              & c_ethanol_ag2(:),                         & ! inout
              & c_nonsoluble_ag2(:),                      & ! inout
              & c_acid_bg1(:),                            & ! inout
              & c_water_bg1(:),                           & ! inout
              & c_ethanol_bg1(:),                         & ! inout
              & c_nonsoluble_bg1(:),                      & ! inout
              & c_acid_bg2(:),                            & ! inout
              & c_water_bg2(:),                           & ! inout
              & c_ethanol_bg2(:),                         & ! inout
              & c_nonsoluble_bg2(:),                      & ! inout
              & c_humus_1(:),                             & ! inout
              & c_humus_2(:),                             & ! inout
              & cflux_fire_all_2_atm                      & ! out
              & )

            CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
            CALL recalc_per_tile_vars(tile, options,                                      &
              & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
              & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
              & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
              & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
              & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
              & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
              & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
              & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
              & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
              & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
              & cflux_dist_greenreserve_2_soil = cflux_dist_greenreserve_2_soil(:),       &
              & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:),                     &
              & cflux_fire_all_2_atm = cflux_fire_all_2_atm(:))

            ! For conservation test: negate currently positve co2flux_fire_all_2_atm_ta (flux away from land!)
            ! and convert back from CO2 flux to carbon change per day
            CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
              & -co2flux_fire_all_2_atm_ta(:) * sec_per_day / molarMassCO2_kg, cconservation_fire(:))

          CASE (2) !! Arora & Boer algorithm
            CALL finish('disturbed_fract','Arora & Boer algorithm not implemented yet.')
          CASE (3)
            CALL finish('disturbed_fract','Thornicke spit fire not implemented yet.')
          CASE (4) !! Read burned_fract
          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown fire algorithm')
        END SELECT

        ! R: This corresponds to JS3, but I did not understand it...
        !    disturbance%disturb_2_atmos corresponds to dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_fire_all_2_atm )
        !    However, I guess it is better not to take this from JS3.
        !IF ( tile%Is_process_active(CARBON_) ) THEN
        !  ! Calculate CO2 and other fluxes to the atmosphere (updated once a day)
        !  CO2_emission(:) = disturbance%disturb_2_atmos(kidx0:kidx1)   * molarMassCO2_kg/86400._wp
        !ENDIF

      ENDIF

    ENDIF ! is_newday

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_natural_disturbances

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "natural_disturbances"
  !!
  !! @param[in,out] tile    Tile for which aggregation of child tiles is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_natural_disturbances(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_natural_disturbance'

    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(DISTURB_)
    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(DISTURB_, burned_fract,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, damaged_fract,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, q_rel_air_climbuf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, q_rel_air_climbuf_yDay,  weighted_by_fract)

    ! R: I do not think these variables should be aggregated:
    !dsl4jsb_Aggregate_onChunk(DISTURB_, max_wind_10m_act,       weighted_by_fract)
    !dsl4jsb_Aggregate_onChunk(DISTURB_, prev_day_max_wind_10m,  weighted_by_fract)
    !dsl4jsb_Aggregate_onChunk(DISTURB_, max_wind_10m,           weighted_by_fract)

    IF ( tile%Is_process_active(CARBON_) ) THEN
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_green_ta,              weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_woods_ta,              weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_reserve_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_greenreserve_2_soil_ta, weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta, weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag1_ta,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag1_ta,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag1_ta,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag1_ta,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag2_ta,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag2_ta,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag2_ta,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag2_ta,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_1_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_2_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  co2flux_fire_all_2_atm_ta, weighted_by_fract)
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_natural_disturbances

#endif
END MODULE mo_disturb_interface
