!> Contains the interfaces to the HD model
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

MODULE mo_hd_interface
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish, message, message_text
  USE mo_sync,             ONLY: SYNC_C, sync_patch_array, global_sum_array

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_grid,           ONLY: Get_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_time,           ONLY: is_time_experiment_start
  USE mo_jsb_math_constants, ONLY: deg2rad

  dsl4jsb_Use_processes HD_, HYDRO_
  dsl4jsb_Use_config(HD_)
  dsl4jsb_Use_memory(HD_)
  dsl4jsb_Use_memory(HYDRO_)

  USE mo_hd_reservoir_cascade, ONLY: reservoir_cascade
  USE mo_jsb_physical_constants, ONLY: rhoh2o

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_hd, hd_lateral_flow, hd_check_water_budget

  PUBLIC :: Register_hd_tasks

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_hd_local_outflow
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_hd
    PROCEDURE, NOPASS :: Aggregate => aggregate_hd
  END TYPE tsk_hd_local_outflow

  INTERFACE tsk_hd_local_outflow
    PROCEDURE Create_task_hd_local_outflow
  END INTERFACE tsk_hd_local_outflow

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_interface'

CONTAINS

  SUBROUTINE Register_hd_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_hd_local_outflow(model_id))

  END SUBROUTINE Register_hd_tasks

  FUNCTION Create_task_hd_local_outflow(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_hd_local_outflow::return_ptr)
    CALL return_ptr%Construct(name='hd_local_outflow', process_id=HD_, owner_model_id=model_id)

  END FUNCTION Create_task_hd_local_outflow

  !------------------------------------------------------------------------------------------------
  SUBROUTINE update_hd(tile, options)
  !------------------------------------------------------------------------------------------------
  ! the routine is called within the block-loop from jsbach_interface. There is no information on
  ! neighbouring grid cells (on other patches) availabe.
  ! We calculate here the local outflow of each grid cell.
  !------------------------------------------------------------------------------------------------

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: discharge
    dsl4jsb_Real2D_onChunk :: outflow_runoff
    dsl4jsb_Real2D_onChunk :: outflow_drainage
    dsl4jsb_Real2D_onChunk :: outflow_rivers
    dsl4jsb_Real2D_onChunk :: ret_overlflow
    dsl4jsb_Real2D_onChunk :: nres_overlflow
    dsl4jsb_Real2D_onChunk :: ret_baseflow
    dsl4jsb_Real2D_onChunk :: nres_baseflow
    dsl4jsb_Real2D_onChunk :: ret_riverflow
    dsl4jsb_Real2D_onChunk :: nres_riverflow
    dsl4jsb_Real3D_onChunk :: overlflow_res
    dsl4jsb_Real3D_onChunk :: baseflow_res
    dsl4jsb_Real3D_onChunk :: riverflow_res

    REAL(wp), POINTER :: area(:)

    REAL(wp), DIMENSION(options%nc) :: &
      & tile_fract, &
      & inflow

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_hd'

    CHARACTER(len=32) :: routing_scheme

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    model => Get_model(tile%owner_model_id)

    IF (.NOT. tile%Is_process_active(HD_)) RETURN

    dsl4jsb_Get_config(HD_)

    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    grid => Get_grid(model%grid_id)

    ! Configuration parameters from the namelist
    routing_scheme    = dsl4jsb_Config(HD_)%routing_scheme

    ! Get reference to variables in memory
    !
    dsl4jsb_Get_memory(HD_)
    !dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    area => grid %area (ics:ice, iblk)

    CALL tile%Get_fraction(ics, ice, iblk, fract=tile_fract(:))

    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, discharge)   ! in

    dsl4jsb_Get_var2D_onChunk(HD_, ret_overlflow)    ! in
    dsl4jsb_Get_var2D_onChunk(HD_, nres_overlflow)   ! in
    dsl4jsb_Get_var2D_onChunk(HD_, ret_baseflow)     ! in
    dsl4jsb_Get_var2D_onChunk(HD_, nres_baseflow)    ! in
    dsl4jsb_Get_var2D_onChunk(HD_, ret_riverflow)    ! in
    dsl4jsb_Get_var2D_onChunk(HD_, nres_riverflow)   ! in
    dsl4jsb_Get_var3D_onChunk(HD_, overlflow_res)    ! inout
    dsl4jsb_Get_var3D_onChunk(HD_, baseflow_res)     ! inout
    dsl4jsb_Get_var3D_onChunk(HD_, riverflow_res)    ! inout
    dsl4jsb_Get_var2D_onChunk(HD_, outflow_runoff)   ! out
    dsl4jsb_Get_var2D_onChunk(HD_, outflow_drainage) ! out
    dsl4jsb_Get_var2D_onChunk(HD_, outflow_rivers)   ! out

    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')

      ! unit conversions: [kg m-2 s-1] -> [m3/s]
      inflow(:) = runoff(:) * tile_fract(:) * area / rhoh2o
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_overlflow(:), nres_overlflow(:), overlflow_res(:,:), outflow_runoff(:))

      inflow(:) = drainage(:) * tile_fract(:) * area / rhoh2o
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_baseflow(:),  nres_baseflow(:),  baseflow_res(:,:),  outflow_drainage(:))

      inflow(:) = discharge(:)
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_riverflow(:), nres_riverflow(:), riverflow_res(:,:), outflow_rivers(:))

    CASE ('weighted_to_coast')
      !
      ! dummy version of HD model
      !
      ! unit conversion: [kg m-2 s-1] -> [m3/s]
      discharge(:) = tile_fract(:) * (runoff(:) + drainage(:)) * area(:) / rhoh2o

    CASE ('zero')
      !
      ! even more dummy: set discharge to zero
      !
      discharge(:) = 0._wp

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for routing_scheme: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_hd

  SUBROUTINE aggregate_hd(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    IF (options%nc > 0) CONTINUE ! only to avoid compiler warnings
    IF (tile%level > 0) CONTINUE ! only to avoid compiler warnings

    ! Nothing to do because the HD model only runs on the top tile

  END SUBROUTINE aggregate_hd
  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_lateral_flow(tile, options)
  !------------------------------------------------------------------------------------------------
  ! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  ! It is thus possible to synchronize the halos. And to use data from surounding grid cells.
  !------------------------------------------------------------------------------------------------

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    INTEGER                       :: jc, blk, n
    INTEGER                       :: jc_us, blk_us

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':hd_lateral_flow'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: &
      & hd_mask,               &
      & coast_ocean,           &
      & outflow_runoff,        &
      & outflow_drainage,      &
      & outflow_rivers,        &
      !& outflow_resid,         &
      !& outflow_count,         &
      & discharge,             &
      & discharge_ocean,       &
      & internal_drain,        &
      & nsplit

    dsl4jsb_Real3D_onDomain :: &
      & nidx_upstream,         &
      & bidx_upstream

    REAL(wp) ::                 &
      & global_discharge,       &  ! global sum of discharge to the ocean
      & global_outflow,         &
      & global_outflow_resid,   &
      & global_internal_drain,  &
      & global_discharge_ocean, &
      & sum_weights

    INTEGER  ::         &
      & nproma,         &
      & nblks,          &
      & blks, blke, jcs, jce

    REAL(wp), ALLOCATABLE :: cos_lat(:,:) ! Cosines of latitudes of cell centers

    LOGICAL, POINTER  ::  &
      & is_in_domain    (:,:) ! T: cell in domain (not halo)

    CHARACTER(len=32) :: &
      & routing_scheme        ! river routing scheme

    REAL(wp), ALLOCATABLE  :: &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell

    LOGICAL  :: debug_hd,use_bifurcated_rivers

    IF (options%nc > 0) CONTINUE ! only to avoid compiler warnings

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD_)

    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HD_)
    !dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    nproma =  grid%nproma
    nblks  =  grid%nblks

    dsl4jsb_Get_var2D_onDomain(HD_, hd_mask)          ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_runoff)   ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_drainage) ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_rivers)   ! IN
    !dsl4jsb_Get_var2D_onDomain(HD_, outflow_resid)    ! OUT
    !dsl4jsb_Get_var2D_onDomain(HD_, outflow_count)    ! OUT
    dsl4jsb_Get_var3D_onDomain(HD_, nidx_upstream)    ! IN
    dsl4jsb_Get_var3D_onDomain(HD_, bidx_upstream)    ! IN

    dsl4jsb_Get_var2D_onDomain(HD_, coast_ocean) ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, nsplit) ! IN

    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge)       ! OUT
    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge_ocean) ! OUT
    dsl4jsb_Get_var2D_onDomain(HYDRO_, internal_drain)  ! OUT

    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    ! Configuration parameters from the namelist
    routing_scheme        = dsl4jsb_Config(HD_)%routing_scheme
    debug_hd              = dsl4jsb_Config(HD_)%debug_hd
    use_bifurcated_rivers = dsl4jsb_Config(HD_)%use_bifurcated_rivers


    ! Domain Mask - to mask all halo cells for global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(nproma,nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')

      ! Synchronize halos
      CALL sync_patch_array(sync_c, grid%patch, outflow_runoff(:,:))
      CALL sync_patch_array(sync_c, grid%patch, outflow_drainage(:,:))
      CALL sync_patch_array(sync_c, grid%patch, outflow_rivers(:,:))

      discharge(:,:) = 0._wp
      discharge_ocean(:,:) = 0._wp
      internal_drain(:,:) = 0._wp
      ! Enable the next line for debugging ... but this only works on one (!) processor, i.e. no halos
      ! outflow_count(:,:) = nsplit(:,:)

      global_outflow        = global_sum_array((outflow_runoff(:,:)+outflow_drainage(:,:)+outflow_rivers(:,:))*in_domain(:,:))

      blks = grid%get_blk_start()
      blke = grid%get_blk_end()
      DO blk = blks, blke
        jcs = grid%get_col_start(blk)
        jce = grid%get_col_end(blk)
        DO jc = jcs, jce
          ! inflow from grid cells upstream: outflow from runnoff, drainage and riverflow
          DO n = 1, hd__mem%nneigh
            ! local indices of grid cells upstream
            jc_us  =  nidx_upstream(jc,n,blk)
            blk_us =  bidx_upstream(jc,n,blk)
            IF (jc_us > 0 .and. blk_us > 0) THEN
            ! Workaround for NVIDIA compiler < 22.9 internal compiler error
            ! Note that nsplit is set to 1._wp for use_bifurcated_rivers=false
#if defined(__NVCOMPILER) && (__NVCOMPILER_MAJOR__ < 22  || (__NVCOMPILER_MAJOR__ == 22 && __NVCOMPILER_MINOR__ < 9))
              discharge(jc,blk) =  discharge(jc,blk)                      &
                                 + ( ( outflow_runoff(jc_us,blk_us)       &
                                     + outflow_drainage(jc_us,blk_us)     &
                                     + outflow_rivers(jc_us,blk_us)   )   &
                                   / nsplit(jc_us,blk_us)               )
#else
              IF (use_bifurcated_rivers) THEN
                discharge(jc,blk) =  discharge(jc,blk)                      &
                                   + ( ( outflow_runoff(jc_us,blk_us)       &
                                       + outflow_drainage(jc_us,blk_us)     &
                                       + outflow_rivers(jc_us,blk_us)   )   &
                                     / nsplit(jc_us,blk_us)               )
              ELSE
                discharge(jc,blk) =  discharge(jc,blk)               &
                                   + outflow_runoff(jc_us,blk_us)    &
                                   + outflow_drainage(jc_us,blk_us)  &
                                   + outflow_rivers(jc_us,blk_us)
              END IF
#endif
              ! Enable the next six lines for debugging ... but this only works on one (!) processor, i.e. no halos
              !outflow_count(jc_us,blk_us) = outflow_count(jc_us,blk_us) - 1._wp
              !IF (outflow_count(jc_us,blk_us) == 0) THEN
              !  outflow_runoff  (jc_us,blk_us) = 0._wp
              !  outflow_drainage(jc_us,blk_us) = 0._wp
              !  outflow_rivers  (jc_us,blk_us) = 0._wp
              !END IF
            END IF
          END DO
        END DO
      END DO

      WHERE (hd_mask(:,:) == 0._wp)
        ! ocean inflow cells
        discharge_ocean(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + discharge(:,:)
        outflow_runoff(:,:) = 0._wp
        outflow_drainage(:,:) = 0._wp
        discharge(:,:) = 0._wp
      ELSE WHERE (hd_mask(:,:) == 2._wp)
        ! internal drainage cells
        internal_drain(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + discharge(:,:)
        outflow_runoff(:,:) = 0._wp
        outflow_drainage(:,:) = 0._wp
        discharge(:,:) = 0._wp
      END WHERE

      ! Enable this block for debugging ... but this only works on one (!) processor, i.e. no halos
      !outflow_resid(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + outflow_rivers(:,:)
      !DO blk = blks, blke
      !  jcs = grid%get_col_start(blk)
      !  jce = grid%get_col_end(blk)
      !  DO jc = jcs, jce
      !    IF (outflow_resid(jc,blk) /= 0._wp) print*,'AAA ', jc,blk, outflow_resid(jc,blk)
      !    IF (outflow_count(jc,blk) /= 0._wp .AND. &
      !             hd_mask(jc,blk) == 1._wp ) print*,'BBB ', jc,blk, outflow_count(jc,blk)
      !  END DO
      !END DO

      ! global integrals
      ! Switch the next two lines for debugging ... but this only works on one (!) processor, i.e. no halos
      ! global_outflow_resid  = global_sum_array(outflow_resid(:,:)*in_domain(:,:))
      global_outflow_resid   = 0._wp
      global_discharge       = global_sum_array(discharge(:,:)*in_domain(:,:))
      global_internal_drain  = global_sum_array(internal_drain(:,:)*in_domain(:,:))
      global_discharge_ocean = global_sum_array(discharge_ocean(:,:)*in_domain(:,:))

      IF (debug_hd) THEN
        ! water conservation test
        WRITE (message_text,*) 'outflow from reservoirs: ', global_outflow, ' m3/s'
        CALL message (routine, message_text)
        ! Enable the next two lines for debugging ... but this only works on one (!) processor, i.e. no halos
        ! WRITE (message_text,*) 'outflow residual: ', global_outflow_resid, ' m3/s'
        ! CALL message (routine, message_text)
        WRITE (message_text,*) 'global discharge (land): ', global_discharge, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'discharge to the ocean : ', global_discharge_ocean, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'internal drainage      : ', global_internal_drain, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'Water budget error     : ', &
          & global_outflow - global_outflow_resid - global_discharge - global_discharge_ocean - global_internal_drain, ' m3/s'
        CALL message (routine, message_text)
      END IF

      ! distribute internal drainage to ocean inflow cells (weighted with ocean inflow)
      IF (global_discharge_ocean > 0._wp) THEN
        discharge_ocean(:,:) = discharge_ocean(:,:) + &
          & (discharge_ocean(:,:) * (global_internal_drain+global_outflow_resid) / global_discharge_ocean)
      ELSE IF (global_internal_drain+global_outflow_resid > 0._wp) THEN
        WRITE (message_text,*) 'no discharge to ocean! Internal drainage of ', &
          & global_internal_drain+global_outflow_resid, ' m3/s will be lost.'
        CALL finish (routine, message_text)
      END IF

      IF (debug_hd) THEN
        ! water conservation test
        CALL message (routine, 'check re-distribution of water from internal drainage')
        global_discharge_ocean = global_sum_array(discharge_ocean(:,:)*in_domain(:,:))
        WRITE (message_text,*) 'discharge to ocean now : ', global_discharge_ocean, ' m3/s'
        CALL message (routine, message_text)
      END IF

    CASE ('weighted_to_coast')

      ! Dummy discharge: the discharge is distributed to all coastal ocean cells, weighted by cosine of latitude

      global_discharge = global_sum_array(discharge(:,:) * in_domain(:,:))

      ALLOCATE(cos_lat(SIZE(is_in_domain,1),SIZE(is_in_domain,2)))
      cos_lat(:,:) = MERGE(COS(deg2rad * grid%lat(:,:)), 0._wp, coast_ocean(:,:) > 0.5_wp)  ! coast_ocean is either 0. or 1.
      sum_weights = global_sum_array(coast_ocean(:,:) * cos_lat(:,:) * in_domain(:,:))
      WHERE (coast_ocean(:,:) > 0.5_wp)
        discharge_ocean(:,:) = global_discharge * cos_lat(:,:) / sum_weights
      ELSEWHERE
        discharge_ocean(:,:) = 0._wp
      END WHERE
      ! discharge(:,:) = 0._wp

    CASE ('zero')

      discharge_ocean(:,:) = 0._wp

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for routing_scheme: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

    DEALLOCATE (in_domain)

  END SUBROUTINE hd_lateral_flow

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_check_water_budget(tile, options)
  !------------------------------------------------------------------------------------------------
  ! calculation of the global water budget at the end of a time step, and compare it to the
  ! budget from the previous time step.
  !------------------------------------------------------------------------------------------------

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':hd_check_water_budget'

    INTEGER :: nres_o_max, nres_b_max, nres_r_max
    INTEGER :: nproma, nblks

    REAL(wp), POINTER ::       &
      & area            (:,:)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & runoff          , &
      & drainage        , &
      !& outflow_resid   , &
      !& outflow_count   , &
      & discharge       , &
      & discharge_ocean , &
      & water_budget    , &
      & water_budget_old, &
      & water_error     , &
      & water_flux

    dsl4jsb_Real3D_onDomain :: &
      & overlflow_res ,        &
      & baseflow_res  ,        &
      & riverflow_res

    LOGICAL,  POINTER ::    &
      & is_in_domain (:,:)      ! T: cell in domain (not halo)

    REAL(wp), ALLOCATABLE  ::  &
      & tile_fract (:,:), &
      & in_domain  (:,:) ! 1: cell in domain, 0: halo cell

    REAL(wp) ::                   &
      & steplen,                  &
      & global_water_budget,      &  ! water budget of current step
      & global_water_budget_old,  &  ! water budget of previous step
      & global_water_error,       &  ! water budget change within the time step
      & global_water_flux            ! water fluxes changing the budget

    LOGICAL :: debug_hd, diag_water_budget

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD_)

    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    steplen = options%steplen

    dsl4jsb_Get_memory(HD_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Configuration parameters from the namelist
    debug_hd          = dsl4jsb_Config(HD_)%debug_hd
    diag_water_budget = dsl4jsb_Config(HD_)%diag_water_budget

    ! Get reference to variables in memory

    nres_o_max = hd__mem%nres_o_max
    nres_b_max = hd__mem%nres_b_max
    nres_r_max = hd__mem%nres_r_max

    nproma          = grid%nproma
    nblks           = grid%nblks

    area            => grid%area                                (:,:)
    is_in_domain    => grid%patch%cells%decomp_info%owner_mask  (:,:)

    dsl4jsb_Get_var3D_onDomain(HD_, overlflow_res) ! IN
    dsl4jsb_Get_var3D_onDomain(HD_, baseflow_res) ! IN
    dsl4jsb_Get_var3D_onDomain(HD_, riverflow_res) ! IN

    dsl4jsb_Get_var2D_onDomain(HYDRO_, runoff)          ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, drainage)        ! IN
    !dsl4jsb_Get_var2D_onDomain(HD_   , outflow_resid)   ! IN
    !dsl4jsb_Get_var2D_onDomain(HD_   , outflow_count)   ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge)       ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge_ocean) ! IN

    dsl4jsb_Get_var2D_onDomain(HD_, water_budget)       ! out
    dsl4jsb_Get_var2D_onDomain(HD_, water_budget_old)   ! out
    dsl4jsb_Get_var2D_onDomain(HD_, water_flux)         ! out
    dsl4jsb_Get_var2D_onDomain(HD_, water_error)        ! out

    ! Mask out all halo cells for the global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(nproma, nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSE WHERE
      in_domain = 0._wp
    END WHERE

    IF (TRIM(dsl4jsb_Config(HD_)%routing_scheme) /= 'full') THEN
      IF (diag_water_budget) THEN
        WRITE (message_text,*) '  global discharge to the ocean: ', global_sum_array(discharge_ocean*in_domain), ' m3/s'
        CALL message (TRIM(routine),message_text)
      END IF
      RETURN
    END IF

    water_budget_old(:,:) = water_budget(:,:) * in_domain(:,:)
    global_water_budget_old = global_sum_array(water_budget_old(:,:))

    ! re-calculate global land water budget
    !
    water_budget(:,:) = &  ! [m3]
      &   SUM(overlflow_res(:,:,:), DIM=2) * in_domain(:,:) & !< water in overland flow reservoirs
      & + SUM(baseflow_res (:,:,:), DIM=2) * in_domain(:,:) & !< water in baseflow reservoir
      & + SUM(riverflow_res(:,:,:), DIM=2) * in_domain(:,:) & !< water in riverflow reservoirs
      & + discharge(:,:) * in_domain(:,:) * steplen                              !< river discharge [m3/s -> m3]
    global_water_budget = global_sum_array(water_budget(:,:))

    IF (debug_hd) THEN
      CALL message(routine, 'global water budget:')
      WRITE(message_text,*) '-   overlandflow reservoir: ', &
           global_sum_array(overlflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_o_max)), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   baseflow reservoir    : ', &
           global_sum_array(baseflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_b_max)), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   riverflow reservoir   : ', &
           global_sum_array(riverflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_r_max)), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   discharge             : ', &
           global_sum_array(discharge*in_domain)*steplen, ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-->          total budget : ', global_water_budget, ' m3'
      CALL message(routine, message_text)
    END IF

    ALLOCATE(tile_fract(nproma, nblks))
    CALL tile%Get_fraction(fract=tile_fract(:,:))
    water_flux(:,:) = &
      &   runoff         (:,:) * area(:,:)/rhoh2o*tile_fract(:,:)*in_domain(:,:)  & !< runoff   [kg m-2 s-1] -> [m3/s]
      & + drainage       (:,:) * area(:,:)/rhoh2o*tile_fract(:,:)*in_domain(:,:)  & !< drainage [kg m-2 s-1] -> [m3/s]
      & - discharge_ocean(:,:) * in_domain(:,:)                                     !< discharge to the ocean
    global_water_flux = global_sum_array(water_flux(:,:))

    IF (debug_hd) THEN
      CALL message(routine, 'global water fluxes during time interval:')
      WRITE(message_text,*) '-   runoff (incl. P-E on glacier/lakes) : ', &
           global_sum_array(runoff/rhoh2o*area*tile_fract*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   drainage              : ', &
           global_sum_array(drainage/rhoh2o*area*tile_fract*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   discharge_ocean (neg.): ', &
           global_sum_array(discharge_ocean*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-->      sum of the fluxes: ', global_water_flux*steplen, ' m3'
      CALL message(routine, message_text)
    END IF

    IF (diag_water_budget) THEN
      WRITE (message_text,*) '  global discharge to the ocean: ', global_sum_array(discharge_ocean*in_domain), ' m3/s'
      CALL message (TRIM(routine),message_text)
    END IF

    IF (is_time_experiment_start(options%current_datetime)) RETURN

    ! water conservation test
    !
    water_error(:,:) = water_budget_old(:,:) + water_flux(:,:) * steplen - water_budget(:,:)
    global_water_error = global_sum_array(water_error(:,:))

    IF (debug_hd) THEN
      WRITE (message_text,*) 'Water budget imbalance during time step: ', global_water_error,' m3'
      CALL message (routine,message_text)
    END IF

    IF (global_water_budget_old /= 0._wp .AND. global_water_budget /= 0._wp) THEN        !vg: zero at the very first time step
      IF (ABS(global_water_error/global_water_budget) > 10._wp*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Water conservation problem: budget imbalance: ', &
             global_water_error,' m3'
        IF (dsl4jsb_Config(HD_)%enforce_water_budget) THEN
          CALL finish (TRIM(routine), message_text)
        ELSE
          CALL message (TRIM(routine), message_text)
        END IF
      END IF
    END IF

    DEALLOCATE (in_domain)

  END SUBROUTINE hd_check_water_budget

#endif
END MODULE mo_hd_interface
