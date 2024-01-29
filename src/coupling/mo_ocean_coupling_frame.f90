!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ocean_coupling_frame

  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: warning, message
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_mpi,                 ONLY: p_pe_work, p_comm_work, p_sum
  USE mo_run_config,          ONLY: ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, &
       &                            timer_coupling_init
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mtime,                  ONLY: datetimeToString, MAX_DATETIME_STR_LEN, &
       &                            timedeltaToString, MAX_TIMEDELTA_STR_LEN

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_config,     ONLY: nproma
  USE mo_yac_finterface,      ONLY: yac_fdef_comp, yac_fdef_comps,          &
    &                               yac_fget_version,                       &
    &                               yac_fdef_datetime, yac_fdef_grid,       &
    &                               yac_fdef_points, yac_fset_global_index, &
    &                               yac_fset_core_mask, yac_fdef_mask,      &
    &                               yac_fdef_field_mask, yac_fenddef,       &
    &                               YAC_LOCATION_CELL, YAC_TIME_UNIT_ISO_FORMAT, &
    &                               YAC_LOCATION_CORNER
  USE mo_coupling_config,     ONLY: is_coupled_run, is_coupled_to_output
  USE mo_time_config,         ONLY: time_config
  USE mo_output_coupling,     ONLY: construct_output_coupling, winnow_field_list

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_ocean_coupling, destruct_ocean_coupling
  PUBLIC :: nbr_inner_cells, field_id

  INTEGER, PARAMETER    :: no_of_fields = 14
  INTEGER               :: field_id(no_of_fields)

  INTEGER, SAVE         :: nbr_inner_cells

CONTAINS

  !--------------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and ocean. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d

    CHARACTER(LEN=max_char_length) ::  field_name(no_of_fields)
    INTEGER                        :: collection_size(no_of_fields)
    INTEGER :: error_status

    INTEGER                :: patch_no
    TYPE(t_patch), POINTER :: patch_horz

    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    INTEGER :: comp_ids(2)
    INTEGER :: cell_point_ids(1), vertex_point_ids(1)
    INTEGER :: cell_mask_ids(2)
    INTEGER :: grid_id
    INTEGER :: nbr_vertices_per_cell

    INTEGER :: mask_checksum
    INTEGER :: nblks
    INTEGER :: blockNo, cell_index, nn

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER, ALLOCATABLE  :: buffer_c(:,:)
    LOGICAL, ALLOCATABLE  :: is_valid(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: timestepstring

    IF (.NOT. is_coupled_run()) RETURN

    IF (ltimer) CALL timer_start(timer_coupling_init)

    comp_name = TRIM(get_my_process_name())

    patch_no = 1
    patch_horz => patch_3d%p_patch_2d(patch_no)

    ! Inform the coupler about what we are
    IF( is_coupled_to_output() ) THEN
       CALL yac_fdef_comps ( [TRIM(comp_name)//"       ", TRIM(comp_name)//"_output"], 2, comp_ids )
    ELSE
       CALL yac_fdef_comp ( TRIM(comp_name), comp_ids(1) )
    ENDIF

    ! Print the YAC version
    CALL message('Running ICON ocean in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)
    CALL timedeltaToString(time_config%tc_dt_model, timestepstring)

    CALL yac_fdef_datetime ( start_datetime = TRIM(startdatestring), &
         &                   end_datetime   = TRIM(stopdatestring)   )

    ! Announce one subdomain (patch) to the coupler
    grid_name = "icon_ocean_grid"

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks))
    ALLOCATE(buffer_lat(nproma*nblks))
    ALLOCATE(buffer_c(3,nproma*nblks))

    nbr_vertices_per_cell = 3

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(blockNo, cell_index, nn) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_v
      DO cell_index = 1, nproma
        nn = (blockNo-1)*nproma+cell_index
        buffer_lon(nn) = patch_horz%verts%vertex(cell_index,blockNo)%lon
        buffer_lat(nn) = patch_horz%verts%vertex(cell_index,blockNo)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(blockNo, cell_index, nn) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      DO cell_index = 1, nproma
        nn = (blockNo-1)*nproma+cell_index
        buffer_c(1,nn) = (patch_horz%cells%vertex_blk(cell_index,blockNo,1)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(cell_index,blockNo,1)
        buffer_c(2,nn) = (patch_horz%cells%vertex_blk(cell_index,blockNo,2)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(cell_index,blockNo,2)
        buffer_c(3,nn) = (patch_horz%cells%vertex_blk(cell_index,blockNo,3)-1)*nproma + &
                          patch_horz%cells%vertex_idx(cell_index,blockNo,3)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_grid (          &
      & TRIM(grid_name),          &
      & patch_horz%n_patch_verts, &
      & patch_horz%n_patch_cells, &
      & nbr_vertices_per_cell,    &
      & buffer_lon,               &
      & buffer_lat,               &
      & buffer_c,                 &
      & grid_id )

        ! vertex points (needed for output_coupling)
    CALL yac_fdef_points (        &
         & grid_id,                  &
         & patch_horz%n_patch_verts, &
         & YAC_LOCATION_CORNER,      &
         & buffer_lon(1:patch_horz%n_patch_verts),               &
         & buffer_lat(1:patch_horz%n_patch_verts),               &
         & vertex_point_ids(1) )

    ! Can we have two fdef_point calls for the same subdomain, i.e.
    ! one single set of cells?
    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      DO cell_index = 1, nproma
        nn = (blockNo-1)*nproma+cell_index
        buffer_lon(nn) = patch_horz%cells%center(cell_index,blockNo)%lon
        buffer_lat(nn) = patch_horz%cells%center(cell_index,blockNo)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour)
    CALL yac_fdef_points (        &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & buffer_lon,               &
      & buffer_lat,               &
      & cell_point_ids(1) )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    CALL yac_fset_global_index (                &
      & patch_horz%cells%decomp_info%glb_index, &
      & YAC_LOCATION_CELL,                      &
      & grid_id )

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c))

    nbr_inner_cells = 0
!ICON_OMP_PARALLEL DO PRIVATE(cell_index) REDUCTION(+:nbr_inner_cells) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_index = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(cell_index) ) THEN
         is_valid(cell_index) =.TRUE.
         nbr_inner_cells = nbr_inner_cells + 1
       ELSE
         is_valid(cell_index) = .FALSE.
       ENDIF
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    CALL yac_fset_core_mask ( &
      & is_valid,             &
      & YAC_LOCATION_CELL,    &
      & grid_id )

    IF( is_coupled_to_output() ) THEN

      CALL construct_output_coupling ( &
        patch_3d%p_patch_2d(1:), comp_ids(2), cell_point_ids, vertex_point_ids, &
        timestepstring)

    END IF


    !
    ! mask generation : ... not yet defined ...
    !
    ! We could use the patch_horz%cells%decomp_info%owner_local information
    ! e.g. to mask out halo points. We do we get the info about what is local and what
    ! is remote.
    !
    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! This integer mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    ! The logical mask for the coupler is set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.

    mask_checksum = 0
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,cell_index) REDUCTION(+:mask_checksum) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      DO cell_index = 1, nproma
        mask_checksum = mask_checksum + ABS(patch_3d%surface_cell_sea_land_mask(cell_index, blockNo))
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    mask_checksum = p_sum(mask_checksum, comm=p_comm_work)

    IF ( mask_checksum > 0 ) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        DO cell_index = 1, nproma
          IF ( patch_3d%surface_cell_sea_land_mask(cell_index, blockNo) < 0 ) THEN
            ! ocean and ocean-coast is valid (-2, -1)
            is_valid((blockNo-1)*nproma+cell_index) = .TRUE.
          ELSE
            ! land is undef (1, 2)
            is_valid((blockNo-1)*nproma+cell_index) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ELSE

!ICON_OMP_PARALLEL_DO PRIVATE(cell_index) ICON_OMP_DEFAULT_SCHEDULE
      DO cell_index = 1, patch_horz%nblks_c * nproma
        is_valid(cell_index) = .TRUE.
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ENDIF

    CALL yac_fdef_mask (          &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & is_valid,                 &
      & cell_mask_ids(1) )

    field_name(1) = "surface_downward_eastward_stress"   ! bundled field containing two components
    collection_size(1) = 2
    field_name(2) = "surface_downward_northward_stress"  ! bundled field containing two components
    collection_size(2) = 2
    field_name(3) = "surface_fresh_water_flux"           ! bundled field containing three components
    collection_size(3) = 3
    field_name(4) = "total_heat_flux"                    ! bundled field containing four components
    collection_size(4) = 4
    field_name(5) = "atmosphere_sea_ice_bundle"          ! bundled field containing two components
    collection_size(5) = 2
    field_name(6) = "sea_surface_temperature"
    collection_size(6) = 1
    field_name(7) = "eastward_sea_water_velocity"
    collection_size(7) = 1
    field_name(8) = "northward_sea_water_velocity"
    collection_size(8) = 1
    field_name(9) = "ocean_sea_ice_bundle"               ! bundled field containing three components
    collection_size(9) = 3
    field_name(10) = "10m_wind_speed"
    collection_size(10) = 1
    field_name(11) = "river_runoff"
    collection_size(11) = 1
    field_name(12) = "co2_mixing_ratio"
    collection_size(12) = 1
    field_name(13) = "co2_flux"
    collection_size(13) = 1
    field_name(14) = "sea_level_pressure"
    collection_size(14) = 1

    ! Define the mask for all fields but the runoff

    DO cell_index = 1, no_of_fields 
      if(field_name(cell_index).ne. "river_runoff")then
      CALL yac_fdef_field_mask (        &
        & TRIM(field_name(cell_index)), &
        & comp_ids(1),                   &
        & cell_point_ids,               &
        & cell_mask_ids(1),             &
        & 1,                            &
        & collection_size(cell_index),  &
        & timestepstring,                &
        & YAC_TIME_UNIT_ISO_FORMAT,     &
        & field_id(cell_index) )
      endif
    ENDDO

    ! Define cell_mask_ids(2) for runoff: all ocean points are valid.
    !!slo! Define cell_mask_ids(2) for runoff: ocean coastal points only are valid.
    !!slo!  - todo: use same mask as for other ones: all points, better wet points only

    IF ( mask_checksum > 0 ) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        DO cell_index = 1, nproma
            ! ocean coast (-1) is valid
!         IF ( patch_3d%surface_cell_sea_land_mask(cell_index, blockNo) == -1 ) THEN
          ! all ocean points (-1, -2) are valid
          IF ( patch_3d%surface_cell_sea_land_mask(cell_index, blockNo) <= -1 ) THEN
            is_valid((blockNo-1)*nproma+cell_index) = .TRUE.
          ELSE
            ! elsewhere (land or open ocean 1, 2, -2) is undef
            is_valid((blockNo-1)*nproma+cell_index) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ELSE

!ICON_OMP_PARALLEL_DO PRIVATE(cell_index) ICON_OMP_DEFAULT_SCHEDULE
      DO cell_index = 1, patch_horz%nblks_c * nproma
        is_valid(cell_index) = .TRUE.
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ENDIF

    CALL yac_fdef_mask (          &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & is_valid,                 &
      & cell_mask_ids(1) )

    DEALLOCATE(is_valid)

    ! Define the mask for runoff
    !  - new cell_mask_ids(1) shall contain ocean coast points only for source point mapping

    CALL yac_fdef_field_mask (          &
      & TRIM("river_runoff"),           &
      & comp_ids(1),                     &
      & cell_point_ids,                 &
      & cell_mask_ids(1),               &
      & 1,                              &
      & 1,                              &
      & timestepstring,                  &
      & YAC_TIME_UNIT_ISO_FORMAT,       &
      & field_id(11) )

    CALL yac_fenddef ( )
    IF( is_coupled_to_output() ) &
         CALL winnow_field_list()

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_ocean_coupling

  !--------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE destruct_ocean_coupling()

    IF (.NOT. is_coupled_run()) RETURN

  END SUBROUTINE destruct_ocean_coupling

  !--------------------------------------------------------------------------

END MODULE mo_ocean_coupling_frame

