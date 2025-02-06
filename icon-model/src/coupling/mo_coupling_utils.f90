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

! Set of routines shared by various coupling related modules

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_coupling_utils

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, warning, finish
  USE mo_model_domain,    ONLY: t_patch
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
  USE mo_parallel_config, ONLY: nproma
  USE mo_run_config,      ONLY: ltimer
  USE mo_master_control,  ONLY: get_my_process_name
  USE mo_time_config,     ONLY: time_config
  USE mo_mpi,             ONLY: p_pe_work
  USE mtime,              ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_timer,           ONLY: timer_start, timer_stop, timer_coupling_put, &
    &                           timer_coupling_get, timer_coupling_very_1stget, &
    &                           timer_coupling_1stget, timer_coupling_init, &
    &                           timer_coupling_init_def_comp, &
    &                           timer_coupling_init_enddef
#ifdef YAC_coupling
  USE mo_mpi,             ONLY: p_comm_yac
  USE yac,                ONLY: yac_finit, yac_finit_comm, &
    &                           yac_fread_config_yaml, &
    &                           yac_ffinalize, yac_fdef_comp, &
    &                           yac_fdef_comps, yac_fdef_grid, &
    &                           yac_fget_version, yac_fdef_mask, &
    &                           yac_fdef_points, yac_fdef_field, &
    &                           yac_fdef_datetime, yac_fdef_field_mask, &
    &                           yac_fset_global_index, yac_fset_core_mask, &
    &                           yac_fget_action, yac_fupdate, &
    &                           yac_dble_ptr, yac_fput, yac_fget, &
    &                           yac_fget_field_collection_size, &
    &                           yac_fsync_def, yac_fenddef, &
    &                           yac_fget_grid_size, &
    &                           YAC_LOCATION_CELL, &
    &                           YAC_LOCATION_CORNER, &
    &                           YAC_LOCATION_EDGE, &
    &                           YAC_TIME_UNIT_ISO_FORMAT, &
    &                           YAC_ACTION_NONE, &
    &                           YAC_ACTION_COUPLING, &
    &                           YAC_ACTION_PUT_FOR_RESTART, &
    &                           YAC_ACTION_GET_FOR_RESTART, &
    &                           YAC_ACTION_OUT_OF_BOUND
  USE mpi
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cpl_construct
  PUBLIC :: cpl_destruct
  PUBLIC :: cpl_is_initialised
  PUBLIC :: cpl_get_instance_id
  PUBLIC :: cpl_config_file_exists
  PUBLIC :: cpl_def_main
  PUBLIC :: cpl_def_main_dummy
  PUBLIC :: cpl_def_cell_field_mask
  PUBLIC :: cpl_def_field
  PUBLIC :: cpl_get_field
  PUBLIC :: cpl_get_field_collection_size
  PUBLIC :: cpl_put_field
  PUBLIC :: cpl_sync_def
  PUBLIC :: cpl_enddef

  CHARACTER(LEN=*), PARAMETER :: yaml_filename = "coupling.yaml"

  LOGICAL :: yac_is_initialised = .FALSE.
  INTEGER :: yac_instance_id = -1

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_coupling_utils'

  ! register the main component (and optionally the output component)
  INTERFACE cpl_def_main
    MODULE PROCEDURE cpl_def_main_without_output
    MODULE PROCEDURE cpl_def_main_with_output
  END INTERFACE cpl_def_main

  ! registers a field to the coupler
  INTERFACE cpl_def_field
    MODULE PROCEDURE cpl_def_field_no_mask
    MODULE PROCEDURE cpl_def_field_mask
  END INTERFACE cpl_def_field

  ! receives cell-based field data through the coupler
  INTERFACE cpl_get_field
    MODULE PROCEDURE cpl_get_field_idx_lev_blk
    MODULE PROCEDURE cpl_get_field_idx_blk_collection
    MODULE PROCEDURE cpl_get_field_n_collection
  END INTERFACE cpl_get_field

  ! sends cell-based field data through the coupler
  INTERFACE cpl_put_field
    MODULE PROCEDURE cpl_put_field_idx_blk_collection
  END INTERFACE cpl_put_field

CONTAINS

  SUBROUTINE cpl_construct()

#if !defined NOMPI && defined YAC_coupling

    INTEGER :: global_rank, ierror

    IF (ltimer) CALL timer_start (timer_coupling_init)

    yac_is_initialised = .TRUE.

    CALL yac_finit_comm ( p_comm_yac, yac_instance_id )

    IF (ltimer) CALL timer_stop(timer_coupling_init)

#endif

  END SUBROUTINE cpl_construct

  SUBROUTINE cpl_destruct()

#if !defined NOMPI && defined YAC_coupling
    IF (yac_is_initialised) CALL yac_ffinalize( yac_instance_id )
#endif

  END SUBROUTINE cpl_destruct

  FUNCTION cpl_config_file_exists()

    LOGICAL :: cpl_config_file_exists

    LOGICAL, SAVE :: config_files_exist = .FALSE.
    LOGICAL, SAVE :: config_files_have_been_checked = .FALSE.

    LOGICAL :: yaml_exists

    IF (config_files_have_been_checked) THEN

      cpl_config_file_exists = config_files_exist

    ELSE

      INQUIRE(FILE=TRIM(ADJUSTL(yaml_filename)), EXIST=yaml_exists)

      config_files_have_been_checked = .TRUE.
      config_files_exist = yaml_exists
      cpl_config_file_exists = config_files_exist

    END IF

  END FUNCTION cpl_config_file_exists

  FUNCTION cpl_is_initialised()

    LOGICAL :: cpl_is_initialised

    cpl_is_initialised = yac_is_initialised

  END FUNCTION cpl_is_initialised

  FUNCTION cpl_get_instance_id()

    INTEGER :: cpl_get_instance_id

    CHARACTER(*), PARAMETER :: &
      routine = modname // ":cpl_get_instance_id"

    IF (.NOT. yac_is_initialised) &
      CALL finish(routine, "YAC has not been initialised")

    cpl_get_instance_id = yac_instance_id

  END FUNCTION cpl_get_instance_id

  ! registers the main component, grid and points
  SUBROUTINE def_main(                       &
    caller, p_patch, grid_name, with_output, &
    comp_id, output_comp_id, grid_id,        &
    cell_point_id, vertex_point_id,          &
    nbr_inner_cells)

    TYPE(t_patch), INTENT(IN) :: p_patch      ! basic patch
    CHARACTER(LEN=*), INTENT(IN) :: caller    ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: grid_name ! name of the grid
    LOGICAL, INTENT(IN) :: with_output        ! should the output component be registered as well
    INTEGER, INTENT(OUT) :: comp_id           ! component id
    INTEGER, INTENT(OUT) :: output_comp_id    ! component id of the output
    INTEGER, INTENT(OUT) :: grid_id           ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id     ! cell coordinate id
    INTEGER, INTENT(OUT) :: vertex_point_id   ! vertex coordinate id (only with output)
    INTEGER, INTENT(OUT) :: nbr_inner_cells   ! number of core cells

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring

    INTEGER :: jc, jv, jb, nblks, nn, comp_ids(2)
    INTEGER :: comp_comm, comp_rank, ierror

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)

    LOGICAL,  ALLOCATABLE :: is_valid(:)

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':def_main', &
      'built without coupling support.')
#else

    ! Print the YAC version
    CALL message( &
      'Running ICON ' // TRIM(get_my_process_name()) // &
      ' in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)
    CALL yac_fdef_datetime (                  &
      yac_instance_id,                        &
      start_datetime = TRIM(startdatestring), & !in
      end_datetime   = TRIM(stopdatestring)   ) !in

    ! Inform the coupler about what we are
    IF (ltimer) CALL timer_start(timer_coupling_init_def_comp)
    IF (with_output) THEN

      CALL yac_fdef_comps (                       &
        yac_instance_id,                          &
        [TRIM(get_my_process_name())//"       ",  &
         TRIM(get_my_process_name())//"_output"], & !in
         2,                                       & !in
         comp_ids )                                 !out
      comp_id = comp_ids(1)
      output_comp_id = comp_ids(2)
    ELSE
      CALL yac_fdef_comp (           &
        yac_instance_id,             &
        TRIM(get_my_process_name()), & !in
        comp_id )                      !out
      output_comp_id = -1
    END IF
    IF (ltimer) CALL timer_stop(timer_coupling_init_def_comp)

    ! root process of the component reads in the configuration file
    CALL yac_fget_comp_comm(comp_id, comp_comm)
    CALL MPI_COMM_RANK(comp_comm, comp_rank, ierror)
    IF (comp_rank == 0 .AND. cpl_config_file_exists()) &
      CALL yac_fread_config_yaml( yac_instance_id, TRIM(yaml_filename) )
    CALL MPI_Comm_free(comp_comm, ierror)

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! p_patch%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes in rad.

    nblks = MAX(p_patch%nblks_c,p_patch%nblks_v)

    ALLOCATE(                     &
      buffer_lon(nproma * nblks), &
      buffer_lat(nproma * nblks), &
      buffer_c(3, nproma * nblks))

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, jv, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, p_patch%nblks_v
      DO jv = 1, nproma
        nn = (jb-1)*nproma+jv
        buffer_lon(nn) = p_patch%verts%vertex(jv,jb)%lon
        buffer_lat(nn) = p_patch%verts%vertex(jv,jb)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, p_patch%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_c(1,nn) = (p_patch%cells%vertex_blk(jc,jb,1)-1)*nproma + &
          &               p_patch%cells%vertex_idx(jc,jb,1)
        buffer_c(2,nn) = (p_patch%cells%vertex_blk(jc,jb,2)-1)*nproma + &
          &               p_patch%cells%vertex_idx(jc,jb,2)
        buffer_c(3,nn) = (p_patch%cells%vertex_blk(jc,jb,3)-1)*nproma + &
                          p_patch%cells%vertex_idx(jc,jb,3)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! Definition of unstructured horizontal grid
    CALL yac_fdef_grid(                                &
      & grid_name             = TRIM(grid_name),       & !in
      & nbr_vertices          = p_patch%n_patch_verts, & !in
      & nbr_cells             = p_patch%n_patch_cells, & !in
      & nbr_vertices_per_cell = 3,                     & !in
      & x_vertices            = buffer_lon,            & !in
      & y_vertices            = buffer_lat,            & !in
      & cell_to_vertex        = buffer_c,              & !in
      & grid_id               = grid_id)                 !out

    IF (with_output) THEN

      ! the output component also defines fields located at the vertices
      CALL yac_fdef_points(    &
        grid_id,               & !in
        p_patch%n_patch_verts, & !in
        YAC_LOCATION_CORNER,   & !in
        buffer_lon,            & !in
        buffer_lat,            & !in
        vertex_point_id)         !out
    ELSE
      vertex_point_id = -1
    END IF

    ! Define cell center points
    !
    ! cartesian coordinates of cell centers are stored in
    ! p_patch%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, p_patch%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = p_patch%cells%center(jc,jb)%lon
        buffer_lat(nn) = p_patch%cells%center(jc,jb)%lat
      ENDDO
    ENDDO
    !ICON_OMP_END_PARALLEL_DO

    ! Definition of cell center points
    CALL yac_fdef_points (     &
      & grid_id,               & !in
      & p_patch%n_patch_cells, & !in
      & YAC_LOCATION_CELL,     & !in
      & buffer_lon,            & !in
      & buffer_lat,            & !in
      & cell_point_id )          !out

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    nblks = &
      MAX(p_patch%nblks_c, p_patch%nblks_v)
    ALLOCATE(is_valid(nproma*nblks))

    ! set global indices and core masks
    CALL set_basic_info( &
      p_patch%cells%decomp_info, YAC_LOCATION_CELL, .TRUE., nbr_inner_cells)
    CALL set_basic_info( &
      p_patch%verts%decomp_info, YAC_LOCATION_CORNER, .FALSE.)
!     CALL set_basic_info( &
!       p_patch%edges%decomp_info, YAC_LOCATION_EDGE)

    DEALLOCATE (is_valid)

  CONTAINS

    SUBROUTINE set_basic_info(decomp_info, location, set_core_mask, core_count)

      TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info
      INTEGER, INTENT(IN) :: location
      LOGICAL, INTENT(IN) :: set_core_mask
      INTEGER, OPTIONAL, INTENT(OUT) :: core_count

      INTEGER :: i, point_count

      point_count = SIZE(decomp_info%glb_index)

      CALL yac_fset_global_index ( &
        decomp_info%glb_index, location, grid_id)

      ! Generate core mask
      IF (PRESENT(core_count)) THEN
        core_count = 0
!ICON_OMP_PARALLEL_DO PRIVATE(i) REDUCTION(+:core_count) ICON_OMP_RUNTIME_SCHEDULE
        DO i = 1, point_count
          IF ( p_pe_work == decomp_info%owner_local(i) ) THEN
            is_valid(i) = .TRUE.
            core_count = core_count + 1
          ELSE
            is_valid(i) = .FALSE.
          ENDIF
        ENDDO
!ICON_OMP_END_PARALLEL_DO
      ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(i) ICON_OMP_RUNTIME_SCHEDULE
        DO i = 1, point_count
          is_valid(i) = p_pe_work == decomp_info%owner_local(i)
        ENDDO
!ICON_OMP_END_PARALLEL_DO
      END IF

      IF (set_core_mask) THEN
        ! Define core mask
        CALL yac_fset_core_mask ( is_valid, location, grid_id )
      END IF

    END SUBROUTINE set_basic_info

! YAC_coupling
#endif

  END SUBROUTINE def_main

  SUBROUTINE cpl_def_main_without_output( &
    caller, p_patch, grid_name,           &
    comp_id, grid_id, cell_point_id,      &
    nbr_inner_cells)

    TYPE(t_patch), INTENT(IN) :: p_patch      ! basic patch
    CHARACTER(LEN=*), INTENT(IN) :: caller    ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: grid_name ! name of the grid
    INTEGER, INTENT(OUT) :: comp_id           ! component id
    INTEGER, INTENT(OUT) :: grid_id           ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id     ! cell coordinate id
    INTEGER, INTENT(OUT) :: nbr_inner_cells   ! number of core cells

    INTEGER :: dummy_output_comp_id, dummy_vertex_point_id

    CALL def_main(                                        &
      caller // ':cpl_def_main_without_output', p_patch,  &
      grid_name, .FALSE., comp_id, dummy_output_comp_id,  &
      grid_id, cell_point_id, dummy_vertex_point_id,      &
      nbr_inner_cells)

  END SUBROUTINE cpl_def_main_without_output

  SUBROUTINE cpl_def_main_with_output( &
    caller, p_patch, grid_name,        &
    comp_id, output_comp_id, grid_id,  &
    cell_point_id, vertex_point_id,    &
    nbr_inner_cells)

    TYPE(t_patch), INTENT(IN) :: p_patch      ! basic patch
    CHARACTER(LEN=*), INTENT(IN) :: caller    ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: grid_name ! name of the grid
    INTEGER, INTENT(OUT) :: comp_id           ! component id
    INTEGER, INTENT(OUT) :: output_comp_id    ! component id of the output
    INTEGER, INTENT(OUT) :: grid_id           ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id     ! cell coordinate id
    INTEGER, INTENT(OUT) :: vertex_point_id   ! vertex coordinate id (only with output)
    INTEGER, INTENT(OUT) :: nbr_inner_cells   ! number of core cells

    CALL def_main(                                    &
      caller // ':cpl_def_main_with_output', p_patch, &
      grid_name, .TRUE., comp_id, output_comp_id,     &
      grid_id, cell_point_id, vertex_point_id,        &
      nbr_inner_cells)

  END SUBROUTINE cpl_def_main_with_output

  ! registers a dummy main component
  SUBROUTINE cpl_def_main_dummy(caller, comp_name)

    CHARACTER(LEN=*), INTENT(IN) :: caller    ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: comp_name ! component name

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_def_main_dummy', 'built without coupling support.')
#else

    INTEGER :: comp_id

    IF (ltimer) CALL timer_start(timer_coupling_init_def_comp)
    CALL yac_fdef_comp( &
      yac_instance_id, TRIM(comp_name), comp_id)
    IF (ltimer) CALL timer_stop(timer_coupling_init_def_comp)
#endif

  END SUBROUTINE cpl_def_main_dummy

  ! registers a field mask for a grid
  SUBROUTINE cpl_def_cell_field_mask( &
    caller, grid_id, is_valid, mask_id)

    USE, INTRINSIC :: iso_c_binding, ONLY : c_size_t, c_int


    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN)  :: grid_id        ! grid identifier
    LOGICAL, INTENT(IN)  :: is_valid(:)    ! mask values
    INTEGER, INTENT(OUT) :: mask_id        ! mask identifier

#ifdef YAC_coupling
    CALL yac_fdef_mask (                  &
      & grid_id,                          &
      & INT(yac_fget_grid_size(           &
      &     YAC_LOCATION_CELL, grid_id)), &
      & YAC_LOCATION_CELL,                &
      & is_valid,                         &
      & mask_id )
#endif

  END SUBROUTINE cpl_def_cell_field_mask

  ! registers a field to the coupler without a mask
  SUBROUTINE cpl_def_field_no_mask( &
    comp_id, cell_point_id, timestepstring, &
    field_name, collection_size, field_id)

    INTEGER, INTENT(IN) :: comp_id                 ! component id
    INTEGER, INTENT(IN) :: cell_point_id           ! cell coordinate id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring ! time step of the field
    CHARACTER(LEN=*), INTENT(IN) :: field_name     ! name of the field
    INTEGER, INTENT(IN) :: collection_size         ! number of levels/bundle size
    INTEGER, INTENT(OUT) :: field_id               ! id of the field

#ifdef YAC_coupling
    CALL yac_fdef_field (                           &
      & field_name      = TRIM(field_name),         & !in
      & component_id    = comp_id,                  & !in
      & point_ids       = (/cell_point_id/),        & !in
      & num_pointsets   = 1,                        & !in
      & collection_size = collection_size,          & !in
      & timestep        = timestepstring,           & !in
      & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
      & field_id        = field_id )                  !out
#endif

  END SUBROUTINE cpl_def_field_no_mask

  ! registers a field to the coupler with a mask
  SUBROUTINE cpl_def_field_mask( &
    comp_id, cell_point_id, cell_mask_id, timestepstring, &
    field_name, collection_size, field_id)

    INTEGER, INTENT(IN) :: comp_id                 ! component id
    INTEGER, INTENT(IN) :: cell_point_id           ! cell coordinate id
    INTEGER, INTENT(IN) :: cell_mask_id            ! cell mask id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring ! time step of the field
    CHARACTER(LEN=*), INTENT(IN) :: field_name     ! name of the field
    INTEGER, INTENT(IN) :: collection_size         ! number of levels/bundle size
    INTEGER, INTENT(OUT) :: field_id               ! id of the field

#ifdef YAC_coupling
    CALL yac_fdef_field_mask (                      &
      & field_name      = TRIM(field_name),         & !in
      & component_id    = comp_id,                  & !in
      & point_ids       = (/cell_point_id/),        & !in
      & mask_ids        = (/cell_mask_id/),         & !in
      & num_pointsets   = 1,                        & !in
      & collection_size = collection_size,          & !in
      & timestep        = timestepstring,           & !in
      & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
      & field_id        = field_id )                  !out
#endif

  END SUBROUTINE cpl_def_field_mask

  ! gets the collection size of a field
  ! (only works after the respective field has been definied and
  !  its information has been distributed among all processes either
  !  by a call to yac_fsync_def or yac_fenddef)
  FUNCTION cpl_get_field_collection_size( &
    caller, comp_name, grid_name, field_name)

    CHARACTER(LEN=*), INTENT(IN) :: comp_name  ! name of the component
    CHARACTER(LEN=*), INTENT(IN) :: grid_name  ! name of the grid
    CHARACTER(LEN=*), INTENT(IN) :: field_name ! name of the field
    CHARACTER(LEN=*), INTENT(IN) :: caller     ! name of the calling routine (for debugging)

    INTEGER :: cpl_get_field_collection_size

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_collection_size', &
      'built without coupling support.')
#else

    cpl_get_field_collection_size = &
      yac_fget_field_collection_size(yac_instance_id, &
        comp_name, grid_name, field_name)

! YAC_coupling
#endif

  END FUNCTION cpl_get_field_collection_size

#ifdef YAC_coupling

  ! basic routine for sending a field through the coupler
  SUBROUTINE put( &
    caller, field_id, field_name, field, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller
    INTEGER, INTENT(IN) :: field_id
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    TYPE(yac_dble_ptr), INTENT(IN) :: field(:, :)
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart

    INTEGER :: num_pointsets, collection_size
    INTEGER :: info, ierr

    num_pointsets = SIZE(field, 1)
    collection_size = SIZE(field, 2)

    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fget_action(field_id, info)

    IF (info == YAC_ACTION_NONE) THEN

      ! update internal clock without an actual put
      CALL yac_fupdate(field_id)

    ELSE

      CALL yac_fput( &
        field_id, num_pointsets, collection_size, field, info, ierr)

    END IF

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( info == YAC_ACTION_PUT_FOR_RESTART ) THEN
      CALL message( &
        caller // ':put', &
        'YAC says it is put for restart - ' // TRIM(field_name))
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      CALL warning( &
        caller // ':put', &
        'YAC says put called after end of run - ' // TRIM(field_name))
    ENDIF

    IF (PRESENT(write_restart)) &
      write_restart = (info == YAC_ACTION_PUT_FOR_RESTART)

  END SUBROUTINE put

! YAC_coupling
#endif

  ! sends one or more fields through the coupler
  ! remark:
  !   * field data has the dimensions (nidx,nblk)
  !     with nidx * nblk >= num_points
  !   * number of provided fields has to match the collection size
  !     associated with the provided field id
  SUBROUTINE cpl_put_field_idx_blk_collection( &
    caller, field_collection_id, field_collection_name, num_points, &
    field_1, field_2, field_3, field_4, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller                            ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_collection_id                        ! field id of the field collection
    CHARACTER(LEN=*), INTENT(IN) :: field_collection_name             ! name of the field collection (for debugging)
    INTEGER, INTENT(IN) :: num_points                                 ! number of points in the field data (e.g. number of cells)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN):: field_1(:,:)           ! field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(IN):: field_2(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(IN):: field_3(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(IN):: field_4(:,:) ! optional field data
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart                   ! .TRUE. if it was the last valid put

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_put_field_idx_blk_collection', &
      'built without coupling support.')
#else

    INTEGER :: collection_size

    TYPE(yac_dble_ptr) :: send_field_collection(1, 4)

    collection_size = 1
    send_field_collection(1,collection_size)%p(1:num_points) => &
      field_1(:,:)
    IF (PRESENT(field_2)) THEN
      collection_size = collection_size + 1
      send_field_collection(1,collection_size)%p(1:num_points) => &
        field_2(:,:)
    END IF
    IF (PRESENT(field_3)) THEN
      collection_size = collection_size + 1
      send_field_collection(1,collection_size)%p(1:num_points) => &
        field_3(:,:)
    END IF
    IF (PRESENT(field_4)) THEN
      collection_size = collection_size + 1
      send_field_collection(1,collection_size)%p(1:num_points) => &
        field_4(:,:)
    END IF

    CALL put( &
      caller // ':cpl_put_field_idx_blk_collection', field_collection_id, &
      field_collection_name, send_field_collection(:,1:collection_size), &
      write_restart)

! YAC_coupling
#endif

  END SUBROUTINE cpl_put_field_idx_blk_collection

#ifdef YAC_coupling

  ! basic routine for receiving a field through the coupler
  SUBROUTINE get( &
    caller, field_id, field_name, field, &
    first_get, received_data, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller
    INTEGER, INTENT(IN) :: field_id
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    TYPE(yac_dble_ptr), INTENT(INOUT) :: field(:)
    LOGICAL, INTENT(IN) :: first_get
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart

    INTEGER :: collection_size
    INTEGER :: get_timer, info, ierr

    ! Skip time measurement of the very first yac_fget
    ! as this will measure mainly the wait time caused
    ! by the initialisation of the model components
    ! and does not tell us much about the load balancing
    ! in subsequent calls.
    LOGICAL, SAVE :: lyac_very_1st_get = .TRUE.

    collection_size = SIZE(field, 1)

    IF (ltimer) THEN
      get_timer = &
        MERGE( &
          timer_coupling_very_1stget, &
          MERGE( &
            timer_coupling_1stget, timer_coupling_get, first_get), &
          lyac_very_1st_get)
      CALL timer_start(get_timer)
      lyac_very_1st_get = .FALSE.
    END IF

    CALL yac_fget_action(field_id, info)

    IF (info == YAC_ACTION_NONE) THEN

      ! update internal clock without an actual get
      CALL yac_fupdate(field_id)

    ELSE

      CALL yac_fget(field_id, collection_size, field, info, ierr)

    END IF

    IF (ltimer) CALL timer_stop(get_timer)

    IF ( info == YAC_ACTION_GET_FOR_RESTART ) THEN
      CALL message( &
        caller // ':get', &
        'YAC says it is get for restart - ' // TRIM(field_name))
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      CALL warning( &
        caller // ':get', &
        'YAC says get called after end of run - ' // TRIM(field_name))
    ENDIF

    IF (PRESENT(received_data)) &
      received_data = &
        (info == YAC_ACTION_COUPLING) .OR. (info == YAC_ACTION_GET_FOR_RESTART)
    IF (PRESENT(write_restart)) &
      write_restart = (info == YAC_ACTION_GET_FOR_RESTART)

  END SUBROUTINE get

! YAC_coupling
#endif

  ! receives one or more fields through the coupler
  ! remark:
  !   * field data has the dimensions (num points, collection size)
  !   * all field data is provided in a single contiguous buffer
  !   * collection size has to match the one associated with the provided
  !     field id
  !   * depending on the field and coupling timestep, no data my actually
  !     be received by this call
  SUBROUTINE cpl_get_field_n_collection( &
    caller, field_collection_id, field_collection_name, &
    field_collection, first_get, received_data, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller                ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_collection_id            ! field id of the field collection
    CHARACTER(LEN=*), INTENT(IN) :: field_collection_name ! name of the field collection (for debugging)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT):: &
      field_collection(:,:)                               ! field data
    LOGICAL, OPTIONAL, INTENT(IN) :: first_get            ! is first get of timestep
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data       ! .TRUE. if data was received by this call
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart       ! .TRUE. if it was the last valid get

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_n_collection', &
      'built without coupling support.')
#else

    INTEGER :: num_points, collection_size
    INTEGER :: i
    LOGICAL :: lfirst_get

    TYPE(yac_dble_ptr) :: recv_field_collection(SIZE(field_collection,2))

    num_points = SIZE(field_collection, 1)
    collection_size = SIZE(field_collection, 2)

    DO i = 1, collection_size
      recv_field_collection(i)%p(1:num_points) => field_collection(:,i)
    END DO

    lfirst_get = .FALSE.
    IF (PRESENT(first_get)) lfirst_get = first_get

    CALL get( &
      caller // ':cpl_get_field_n_collection', field_collection_id, &
      field_collection_name, recv_field_collection, &
      lfirst_get, received_data, write_restart)

! YAC_coupling
#endif

  END SUBROUTINE cpl_get_field_n_collection

  ! receives one or more fields through the coupler
  ! remark:
  !   * field data has the dimensions (nidx,nblk)
  !     with nidx * nblk >= num_points
  !   * number of provided fields has to match the collection size
  !     associated with the provided field id
  !   * depending on the field and coupling timestep, no data my actually
  !     be received by this call
  SUBROUTINE cpl_get_field_idx_blk_collection( &
    caller, field_collection_id, field_collection_name, num_points, &
    field_1, field_2, field_3, field_4, first_get, received_data, &
    write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller                               ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_collection_id                           ! field id of the field collection
    CHARACTER(LEN=*), INTENT(IN) :: field_collection_name                ! name of the field collection (for debugging)
    INTEGER, INTENT(IN) :: num_points                                    ! number of points in the field data (e.g. number of cells)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT):: field_1(:,:)           ! field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(INOUT):: field_2(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(INOUT):: field_3(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(INOUT):: field_4(:,:) ! optional field data
    LOGICAL, OPTIONAL, INTENT(IN) :: first_get                           ! is first get of timestep
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data                      ! .TRUE. if data was received by this call
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart                      ! .TRUE. if it was the last valid get

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_idx_blk_collection', &
      'built without coupling support.')
#else

    INTEGER :: collection_size
    LOGICAL :: lfirst_get

    TYPE(yac_dble_ptr) :: recv_field_collection(4)

    collection_size = 1
    recv_field_collection(collection_size)%p(1:num_points) => &
      field_1(:,:)
    IF (PRESENT(field_2)) THEN
      collection_size = collection_size + 1
      recv_field_collection(collection_size)%p(1:num_points) => &
        field_2(:,:)
    END IF
    IF (PRESENT(field_3)) THEN
      collection_size = collection_size + 1
      recv_field_collection(collection_size)%p(1:num_points) => &
        field_3(:,:)
    END IF
    IF (PRESENT(field_4)) THEN
      collection_size = collection_size + 1
      recv_field_collection(collection_size)%p(1:num_points) => &
        field_4(:,:)
    END IF

    lfirst_get = .FALSE.
    IF (PRESENT(first_get)) lfirst_get = first_get

    CALL get( &
      caller // ':cpl_get_field_idx_blk_collection', field_collection_id, &
      field_collection_name, recv_field_collection(1:collection_size), &
      lfirst_get, received_data, write_restart)

! YAC_coupling
#endif

  END SUBROUTINE cpl_get_field_idx_blk_collection

  ! receives multiple levels of a single field through the coupler
  ! remark:
  !   * field data has the dimensions (nidx,nlev,nblk)
  !     with nidx * nblk >= num_points
  !   * receive buffer has the dimensions (num points, nlev_)
  !     with nlev_ >= nlev
  !   * number of levels match the collection size associated with
  !     the provided field id
  !   * depending on the field and coupling timestep, no data my actually
  !     be received by this call
  SUBROUTINE cpl_get_field_idx_lev_blk( &
    caller, field_id, field_name, field, recv_buf, scale_factor, &
    first_get, received_data, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller          ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_id                 ! field id of the field
    CHARACTER(LEN=*), INTENT(IN) :: field_name      ! name of the field (for debugging)
    REAL(wp), INTENT(INOUT) :: field(:,:,:)         ! field data
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: &
      recv_buf(:,:)                                 ! contiguous temporary buffer used by this routine
    REAL(wp), OPTIONAL, INTENT(IN) :: scale_factor  ! optional: multiply whole field by this factor
                                                    ! (only if data was received)
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data ! .TRUE. if data was received by this call
    LOGICAL, OPTIONAL, INTENT(IN) :: first_get      ! is first get of timestep
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart ! .TRUE. if it was the last valid get

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_idx_lev_blk', &
      'built without coupling support.')
#else

    INTEGER :: nidx, nlev, nblk, num_points
    LOGICAL :: coupling
    INTEGER :: i, j, k
    LOGICAL :: lfirst_get

    TYPE(yac_dble_ptr) :: recv_field(SIZE(field, 2))

    nidx = SIZE(field, 1)
    nlev = SIZE(field, 2)
    nblk = SIZE(field, 3)
    num_points = SIZE(recv_buf, 1)

    IF (nlev > SIZE(recv_buf, 2)) &
      CALL finish( &
        TRIM(caller) // ':cpl_get_field_idx_lev_blk', &
        'insufficient recv_buf size')

    ! MoHa:
    !   remarks:
    !     * Since independent fields may use different field mask
    !       cells written by yac_fget may be different as well. Therefore
    !       it is easier if the caller provides the receive buffer. Otherwise
    !       we would have to reinitialise it every time.

    DO i = 1, nlev
      recv_field(i)%p => recv_buf(:,i)

    lfirst_get = .FALSE.
    IF (PRESENT(first_get)) lfirst_get = first_get
    END DO

    CALL get( &
      caller // ':cpl_get_field_idx_lev_blk', field_id, field_name, recv_field, &
      lfirst_get, coupling, write_restart)

    ! if data was received
    IF (coupling) THEN

      ! unpack data
      IF (PRESENT(scale_factor)) THEN
        DO i = 1, nblk
          DO j = 1, nlev
            DO k = 1, nidx
              IF ((i-1) * nidx + k > num_points) CYCLE
              field(k,j,i) = scale_factor * recv_buf((i-1)*nidx+k,j)
            END DO
          END DO
        END DO
      ELSE
        DO i = 1, nblk
          DO j = 1, nlev
            DO k = 1, nidx
              IF ((i-1) * nidx + k > num_points) CYCLE
              field(k,j,i) = recv_buf((i-1)*nidx+k,j)
            END DO
          END DO
        END DO
      END IF

    END IF

    IF (PRESENT(received_data)) received_data = coupling

! YAC_coupling
#endif

  END SUBROUTINE cpl_get_field_idx_lev_blk

  SUBROUTINE cpl_sync_def(caller)

    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_sync_def', 'built without coupling support.')
#else
    CALL yac_fsync_def(yac_instance_id)
#endif

  END SUBROUTINE

  SUBROUTINE cpl_enddef(caller)

    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_enddef', 'built without coupling support.')
#else
    IF (ltimer) CALL timer_start(timer_coupling_init_enddef)
    CALL yac_fenddef(yac_instance_id)
    IF (ltimer) CALL timer_stop(timer_coupling_init_enddef)
#endif

  END SUBROUTINE

END MODULE mo_coupling_utils
