# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"
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

# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/include/omp_definitions.inc" 1
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

# 16 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90" 2
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
# 61 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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

# 129 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

  END SUBROUTINE cpl_construct

  SUBROUTINE cpl_destruct()





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


    CALL finish( &
      TRIM(caller) // ':def_main', &
      'built without coupling support.')
# 415 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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


    CALL finish( &
      TRIM(caller) // ':cpl_def_main_dummy', 'built without coupling support.')
# 483 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

  END SUBROUTINE cpl_def_main_dummy

  ! registers a field mask for a grid
  SUBROUTINE cpl_def_cell_field_mask( &
    caller, grid_id, is_valid, mask_id)

    USE, INTRINSIC :: iso_c_binding, ONLY : c_size_t, c_int


    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN)  :: grid_id        ! grid identifier
    LOGICAL, INTENT(IN)  :: is_valid(:)    ! mask values
    INTEGER, INTENT(OUT) :: mask_id        ! mask identifier

# 507 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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

# 533 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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

# 561 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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


    CALL finish( &
      TRIM(caller) // ':cpl_get_field_collection_size', &
      'built without coupling support.')
# 590 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

  END FUNCTION cpl_get_field_collection_size

# 647 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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


    CALL finish( &
      TRIM(caller) // ':cpl_put_field_idx_blk_collection', &
      'built without coupling support.')
# 704 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

  END SUBROUTINE cpl_put_field_idx_blk_collection

# 781 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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


    CALL finish( &
      TRIM(caller) // ':cpl_get_field_n_collection', &
      'built without coupling support.')
# 832 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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


    CALL finish( &
      TRIM(caller) // ':cpl_get_field_idx_blk_collection', &
      'built without coupling support.')
# 900 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

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


    CALL finish( &
      TRIM(caller) // ':cpl_get_field_idx_lev_blk', &
      'built without coupling support.')
# 1000 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/coupling/mo_coupling_utils.f90"

  END SUBROUTINE cpl_get_field_idx_lev_blk

  SUBROUTINE cpl_sync_def(caller)

    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)


    CALL finish( &
      TRIM(caller) // ':cpl_sync_def', 'built without coupling support.')




  END SUBROUTINE

  SUBROUTINE cpl_enddef(caller)

    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)


    CALL finish( &
      TRIM(caller) // ':cpl_enddef', 'built without coupling support.')






  END SUBROUTINE

END MODULE mo_coupling_utils
