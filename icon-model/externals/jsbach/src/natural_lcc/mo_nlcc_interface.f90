!> Contains the interfaces to the nlcc process
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
MODULE mo_nlcc_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_jsb_time,        ONLY: is_time_experiment_start
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes A2L_, PHENO_, ASSIMI_, CARBON_, NLCC_, DISTURB_

  ! Use of process configurations
  dsl4jsb_Use_config(NLCC_)

  ! Use of process memories
!! TR$  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(NLCC_)
  dsl4jsb_Use_memory(DISTURB_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_nlcc_tasks !, t_nlcc_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nlcc_interface'

  !> Type definition for nlcc task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_nlcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_nlcc     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_nlcc  !< Aggregates computed task variables
  END TYPE tsk_nlcc

  !> Constructor interface for nlcc task
  INTERFACE tsk_nlcc
    PROCEDURE Create_task_nlcc                        !< Constructor function for task
  END INTERFACE tsk_nlcc

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for nlcc task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "nlcc"
  !!
  FUNCTION Create_task_nlcc(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_nlcc::return_ptr)
    CALL return_ptr%Construct(name='nlcc', process_id=NLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_nlcc

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for nlcc process
  !!
  !! @param[in,out] this      Instance of nlcc process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_nlcc_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_nlcc(model_id))

  END SUBROUTINE Register_nlcc_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "nlcc"
  !! Task "nlcc" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_nlcc(tile, options)

    USE mo_nlcc_process,   ONLY: calc_climbuf, bioclim_limits, potential_tree_fpc, fpc_to_cover_fract_pot, &
                                 cover_fract_pot_to_cover_fract, desert_fraction, fpc_daily, scale_fpc
    USE mo_jsb_time,       ONLY: is_newyear, is_newmonth, is_newday
    USE mo_jsb_grid,       ONLY: Get_grid
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_tile_class, ONLY: t_jsb_tile_abstract

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER  ::  model
    TYPE(t_jsb_grid),  POINTER  ::  grid
    CLASS(t_jsb_tile_abstract),  POINTER  ::  ptr_tile
    LOGICAL  :: lstart, init_running_means, init_mmtemp20
    INTEGER  :: iblk, ics, ice, nc, no_children
    INTEGER  :: i, itile
    REAL(wp) :: dtime
    REAL(wp) :: accelerate_nlcc
    REAL(wp), ALLOCATABLE :: fuel(:) !!$ TR this should be the sum of the YASSO AWEN pools
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_nlcc'
    REAL(wp), POINTER  ::  var_collect(:,:)
    LOGICAL, ALLOCATABLE  :: is_woody(:,:)
    LOGICAL, ALLOCATABLE  :: is_grass(:,:)
    LOGICAL, ALLOCATABLE  :: is_pasture(:,:)
    LOGICAL, ALLOCATABLE  :: is_crop(:,:)
    LOGICAL, ALLOCATABLE  :: is_dynamic(:,:)
    REAL(wp), POINTER  ::  bclimit_min_cold_mmtemp(:,:)
    REAL(wp), POINTER  ::  bclimit_max_cold_mmtemp(:,:)
    REAL(wp), POINTER  ::  bclimit_max_warm_mmtemp(:,:)
    REAL(wp), POINTER  ::  bclimit_min_temprange(:,:)
    REAL(wp), POINTER  ::  bclimit_min_gdd(:,:)
    REAL(wp), POINTER  ::  tau_c_woods(:,:)
    REAL(wp), POINTER  ::  NPP_mean_5year(:,:)
    REAL(wp), POINTER  ::  burned_fract(:,:)
    REAL(wp), POINTER  ::  damaged_fract(:,:)
    REAL(wp), POINTER  ::  max_green_bio(:,:)
    REAL(wp), POINTER  ::  sla(:,:)
    REAL(wp), ALLOCATABLE  ::  cover_fract_pot_previous_year(:,:)
    REAL(wp), ALLOCATABLE  ::  cover_fract_sum(:)
    REAL(wp), ALLOCATABLE  :: fract_fpc_max(:,:)

    ! Declare process memories
    dsl4jsb_Def_config(NLCC_)
    dsl4jsb_Def_memory(NLCC_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(DISTURB_)

    ! Declare pointers to variables in nlcc memory
    dsl4jsb_Real2D_onChunk :: seconds_day
    dsl4jsb_Real2D_onChunk :: seconds_month
    dsl4jsb_Real2D_onChunk :: temp_sum_day
    dsl4jsb_Real2D_onChunk :: temp_sum_month
    dsl4jsb_Real2D_onChunk :: min_mmtemp_of_yr
    dsl4jsb_Real2D_onChunk :: max_mmtemp_of_yr
    dsl4jsb_Real2D_onChunk :: min_mmtemp20
    dsl4jsb_Real2D_onChunk :: max_mmtemp20
    dsl4jsb_Real3D_onChunk :: act_fpc
    dsl4jsb_Real2D_onChunk :: bare_fpc
    dsl4jsb_Real3D_onChunk :: pot_fpc
    dsl4jsb_Real3D_onChunk :: cover_fract_pot
    dsl4jsb_Real3D_onChunk :: cover_fract
    dsl4jsb_Real3D_onChunk :: bio_exist
    dsl4jsb_Real2D_onChunk :: gdd_sum_year
    dsl4jsb_Real2D_onChunk :: gdd_prev_year
    dsl4jsb_Real2D_onChunk :: sum_green_bio_memory
    dsl4jsb_Real2D_onChunk :: desert_fpc
    ! Declare pointers to variables in foreign memory
!!$ TR    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: t_air

    ! Local variables
    LOGICAL :: new_year, new_month, new_day

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    !X if necessary: steplen = options%steplen

    ! If process is not to be calculated on this tile, do nothing
    !IF (.NOT. tile%Is_process_calculated(NLCC_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(NLCC_)

    ! Get process memories
    dsl4jsb_Get_memory(NLCC_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(DISTURB_)

    ! import lct information from PFT tiles
    !!$ TR collect_var_real2d has to be replaced by appropriate dsl4jsb invocations!
    ! build masks for different vegetation types (is_dynamic, is_woody, etc)
    no_children = tile%Get_no_of_children()
    var_collect => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%land_cover_class)
    ALLOCATE(is_woody(nc,no_children))
    ALLOCATE(is_grass(nc,no_children))
    ALLOCATE(is_pasture(nc,no_children))
    ALLOCATE(is_crop(nc,no_children))
    ALLOCATE(is_dynamic(nc,no_children))
    is_woody(:,:) = .FALSE.
    is_grass(:,:) = .FALSE.
    is_pasture(:,:) = .FALSE.
    is_crop(:,:) = .FALSE.
    is_dynamic(:,:) = .FALSE.
    DO i = 1,nc
      DO itile = 1,no_children
        IF (NINT(var_collect(i,itile)) == 3 .OR. NINT(var_collect(i,itile)) == 5) is_woody(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 4) is_grass(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 7) is_pasture(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 6) is_crop(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 3 .OR. NINT(var_collect(i,itile)) == 4 .OR. &
            NINT(var_collect(i,itile)) == 5) is_dynamic(i,itile) = .TRUE.
      END DO
    END DO

    ! import bioclimatic limits from lct library
    bclimit_min_cold_mmtemp => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_min_cold_mmtemp)
    bclimit_max_cold_mmtemp => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_max_cold_mmtemp)
    bclimit_max_warm_mmtemp => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_max_warm_mmtemp)
    bclimit_min_temprange => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_min_temprange)
    bclimit_min_gdd => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_min_gdd)
    tau_c_woods => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%tau_c_woods)

    ! import specific leave area from lct library
    sla => carbon__mem%collect_var_real2d(ics, ice, iblk, no_children, carbon__mem%sla)

    ! import annual maximum carbon content of the green pool
    max_green_bio => carbon__mem%collect_var_real2d(ics, ice, iblk, no_children, carbon__mem%max_green_bio)

    ! import 5 year running average of potential NPP
    NPP_mean_5year => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%NPP_mean_5year)

    ! import burned fractions and damaged fractions
    burned_fract => disturb__mem%collect_var_real2d(ics, ice, iblk, no_children, disturb__mem%burned_fract)
    damaged_fract => disturb__mem%collect_var_real2d(ics, ice, iblk, no_children, disturb__mem%damaged_fract)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(NLCC_,seconds_day)      ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,seconds_month)    ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,temp_sum_day)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,temp_sum_month)   ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,min_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,max_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,min_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,max_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,gdd_sum_year)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,gdd_prev_year)    ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,act_fpc)          ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,bare_fpc)         ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,pot_fpc)          ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,cover_fract_pot)  ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,cover_fract)      ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,bio_exist)        ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,sum_green_bio_memory) ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,desert_fpc)       ! inout
!!$ TR    dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)    ! out
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)            ! in

    IF (is_time_experiment_start(options%current_datetime) .OR. dsl4jsb_Config(NLCC_)%init_nlcc) THEN
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF
    new_year = is_newyear(options%current_datetime, options%dtime)
    new_month = is_newmonth(options%current_datetime, options%dtime)
    new_day = is_newday(options%current_datetime, options%dtime)

    ! import cover fractions from grid
    ALLOCATE (cover_fract_sum(nc))
    cover_fract_sum(:) = 0._wp
    grid => Get_grid(model%grid_id)
    DO itile = 1,no_children
      ! Point to the correct pft
      IF (itile == 1) THEN
        ptr_tile => tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      ! get cover fractions
      CALL ptr_tile%Get_fraction(ics, ice, iblk, fract=cover_fract(:,itile))
      cover_fract_sum(:) = cover_fract_sum(:) + cover_fract(:,itile)
    END DO ! itile
    ! scale cover_fract to fractions relative to the veg_tile
    cover_fract(:,:) = cover_fract(:,:) / MAX(SPREAD(cover_fract_sum(:), NCOPIES=no_children, DIM=2),1.e-10_wp)
!!$    WRITE (*,*) 'sum of cover_fract: ',SUM(cover_fract(:,:),DIM=2)

    ! allocate more local variables
    ALLOCATE (cover_fract_pot_previous_year(nc,no_children))
    ALLOCATE (fuel(nc))

    ! initialize act_fpc, bare fpc, and pot_fpc
    IF (lstart) THEN
       act_fpc(:,:) = cover_fract(:,:)
       bare_fpc(:) = 0._wp
       CALL scale_fpc(routine, nc, iblk, no_children, is_dynamic, act_fpc, bare_fpc)
       pot_fpc(:,:) = act_fpc(:,:) !!$ TR to be improved (scale, non-woodies = 0, etc.)
    ENDIF

    accelerate_nlcc    = dsl4jsb_Config(NLCC_)%accelerate_nlcc
    init_running_means = dsl4jsb_Config(NLCC_)%init_running_means
    init_mmtemp20 = .TRUE.
    IF (MAXVAL(min_mmtemp20(:)) > 999._wp .OR. MINVAL(max_mmtemp20(:)) < -999._wp) init_mmtemp20 = .FALSE.
    CALL calc_climbuf(lstart, init_running_means, new_day, new_month, new_year, dtime, t_air, &
                      seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
                      min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20, &
                      gdd_sum_year, gdd_prev_year)

    IF (.NOT. init_mmtemp20 .AND. debug_on() .AND. iblk==1) THEN
       IF (MAXVAL(min_mmtemp20(:)) <= 999._wp .OR. MINVAL(max_mmtemp20(:)) >= -999._wp) &
          CALL message(TRIM(routine), 'initialisation of mmtemp20')
    ENDIF

!!$ TR    WRITE (*,*) 't_air: ', lstart, new_day, t_air(:)
!!$ TR    WRITE (*,*) 'temp_sum_day: ', lstart, new_day, temp_sum_day(:)

    ! annual calculations
    IF (new_year) THEN
       DO itile = 1,no_children
          CALL bioclim_limits (is_dynamic(:,itile), bclimit_min_cold_mmtemp(:,itile), bclimit_max_cold_mmtemp(:,itile),    &
                               bclimit_max_warm_mmtemp(:,itile), bclimit_min_temprange(:,itile), bclimit_min_gdd(:,itile), &
                               min_mmtemp20(:), max_mmtemp20(:), gdd_prev_year(:), bio_exist(:,itile))
       ENDDO
       IF (.NOT. lstart) THEN
          CALL potential_tree_fpc (nc, no_children, is_dynamic(:,:), is_woody(:,:), bio_exist(:,:), &
                                   NPP_mean_5year(:,:), act_fpc(:,:), pot_fpc(:,:))
       ENDIF
       cover_fract_pot_previous_year(:,:) = cover_fract_pot(:,:)
       CALL fpc_to_cover_fract_pot(nc, no_children, is_dynamic(:,:), is_woody(:,:), &
                                   act_fpc(:,:), bare_fpc(:), cover_fract_pot(:,:))
       IF (.NOT. lstart) THEN
          ! new cover fractions
          CALL cover_fract_pot_to_cover_fract(nc, no_children, is_dynamic(:,:), is_woody(:,:), is_pasture(:,:), is_crop(:,:), &
                                              cover_fract_pot(:,:), cover_fract_pot_previous_year(:,:), cover_fract(:,:))
          DO itile = 1,no_children
             ! Point to the correct pft
             IF (itile == 1) THEN
                ptr_tile => tile%Get_first_child_tile()
             ELSE
                ptr_tile => ptr_tile%Get_next_sibling_tile()
             END IF
             ! set cover fractions
             CALL ptr_tile%Set_fraction(ics, ice, iblk, fract=cover_fract(:,itile) * cover_fract_sum(:))
          END DO ! itile
          ! new desert fraction
          CALL desert_fraction(nc, no_children, init_running_means, accelerate_nlcc,   &
                               is_woody(:,:), is_dynamic(:,:),                         &
                               act_fpc(:,:), bare_fpc(:),                              &
                               sla(:,:), max_green_bio(:,:),                           &
                               sum_green_bio_memory(:), desert_fpc(:))
!!$ TR          fract_fpc_max(:) = MAX(1.e-10_wp,1._wp - desert_fpc(:))
          !  calculate fract_fpc_max and hand it back to child tiles (pft tiles)
          ALLOCATE(fract_fpc_max(nc,no_children))
          fract_fpc_max(:,:) = SPREAD(MAX(1.e-10_wp,1._wp - desert_fpc(:)),NCOPIES=no_children, DIM=2)
          CALL pheno__mem%HandDown_var_real2d(ics, ice, iblk, no_children, pheno__mem%fract_fpc_max, fract_fpc_max)
       ENDIF ! not lstart
    ENDIF ! new_year

    ! daily calculations
    IF (new_day) THEN
       CALL fpc_daily(nc, no_children, accelerate_nlcc, is_woody(:,:), is_dynamic(:,:), tau_c_woods(:,:),       &
                      NPP_mean_5year(:,:), bio_exist(:,:), burned_fract(:,:), damaged_fract(:,:), pot_fpc(:,:), &
                      act_fpc(:,:), bare_fpc(:), fuel(:))
    ENDIF

    DEALLOCATE(is_woody, is_grass, is_pasture, is_crop, is_dynamic)
    DEALLOCATE(bclimit_min_cold_mmtemp, bclimit_max_cold_mmtemp, &
               bclimit_max_warm_mmtemp, bclimit_min_temprange, bclimit_min_gdd, tau_c_woods)
    IF (new_year .AND. .NOT. lstart) DEALLOCATE(fract_fpc_max)
    DEALLOCATE(NPP_mean_5year)
    DEALLOCATE(burned_fract, damaged_fract)
    DEALLOCATE(max_green_bio, sla)
    DEALLOCATE(fuel)
    DEALLOCATE(cover_fract_pot_previous_year)
    DEALLOCATE(cover_fract_sum)
    DEALLOCATE(var_collect)


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_nlcc

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "nlcc"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_nlcc(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(NLCC_)

    ! Local variables
    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_nlcc'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(NLCC_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(NLCC_, min_mmtemp20,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, max_mmtemp20,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, gdd_prev_year,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, bare_fpc,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, desert_fpc,            weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(NLCC_, act_fpc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, pot_fpc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, cover_fract_pot,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, cover_fract,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, bio_exist,             weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_nlcc

#endif
END MODULE mo_nlcc_interface
