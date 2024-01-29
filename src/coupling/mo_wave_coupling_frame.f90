! Initialisation of wave-atmosphere coupling
!
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

MODULE mo_wave_coupling_frame

  USE mo_kind,            ONLY: wp
  USE mo_impl_constants,  ONLY: MAX_CHAR_LENGTH, SUCCESS 
  USE mo_exception,       ONLY: finish, message
  USE mo_model_domain,    ONLY: t_patch
  USE mo_parallel_config, ONLY: nproma
  USE mo_run_config,      ONLY: ltimer
  USE mo_master_control,  ONLY: get_my_process_name
  USE mo_mpi,             ONLY: p_pe_work
  USE mo_time_config,     ONLY: time_config
  USE mtime,              ONLY: datetimeToString, timedeltaToString, &
    &                           MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN
  USE mo_yac_finterface,  ONLY: yac_fget_version, yac_fdef_comp,        &
    &                           yac_fdef_datetime, yac_fdef_grid,       &
    &                           yac_fdef_points, yac_fset_global_index, &
    &                           yac_fset_core_mask, yac_fdef_mask,      &
    &                           yac_fdef_field_mask,                    &
    &                           YAC_TIME_UNIT_ISO_FORMAT, YAC_LOCATION_CELL
  USE mo_timer,           ONLY: timer_start, timer_stop, timer_coupling_init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_wave_coupling
  PUBLIC :: nbr_inner_cells
  PUBLIC :: field_id
  PUBLIC :: collection_size
  PUBLIC :: CPF_U10M, CPF_V10M, CPF_FR_SEAICE, CPF_Z0

  ! Output of module for debug
  CHARACTER(len=*), PARAMETER :: str_module = 'mo_wave_coupling_frame'

  INTEGER, PARAMETER    :: no_of_fields = 4
  INTEGER               :: field_id(no_of_fields)
  INTEGER               :: collection_size(no_of_fields)

  INTEGER, SAVE         :: nbr_inner_cells

  ! These constants are only valid AFTER calling construct_wave_coupling !
  INTEGER, PARAMETER :: CPF_U10M      = 1 !< zonal wind speed at 10m above surface
  INTEGER, PARAMETER :: CPF_V10M      = 2 !< meridional wind speed at 10m above surface
  INTEGER, PARAMETER :: CPF_FR_SEAICE = 3 !< fractional seaice cover
  INTEGER, PARAMETER :: CPF_Z0        = 4 !< surface roughness length

CONTAINS

  !>
  !! SUBROUTINE construct_wave_coupling -- the initialisation for the coupling
  !! of wave model and the atmosphere, through a coupler
  !!
  !! Note that the corresponding routine for the ATMO model construct_atmo_wave_coupling
  !! stored in src/atm_coupling/mo_atmo_waves_coupling_frame.f90
  !!
  SUBROUTINE construct_wave_coupling (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)

    CHARACTER(len=*), PARAMETER :: routine = str_module//':construct_wave_coupling'

    CHARACTER(LEN=max_char_length) :: field_name(no_of_fields)

    TYPE(t_patch), POINTER :: patch_horz

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

    CHARACTER(LEN=max_char_length) :: grid_name   ! name of component-specific horizontal grid
    CHARACTER(LEN=max_char_length) :: comp_name   ! component name

    INTEGER :: comp_id              ! component identifier
    INTEGER :: grid_id              ! grid identifier
    INTEGER :: comp_ids(1)
    INTEGER :: cell_point_ids(1)
    INTEGER :: cell_mask_ids(1)

    INTEGER :: jg
    INTEGER :: nblks
    INTEGER :: jb, jc, jv, nn
    INTEGER :: jn                  ! loop index for fields
    INTEGER :: nlen                ! block length
    INTEGER :: nbr_vertices_per_cell
    INTEGER :: ist                 ! error status

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)

    LOGICAL,  ALLOCATABLE :: is_valid(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN):: tc_dt_model_string
    

    IF (ltimer) CALL timer_start (timer_coupling_init)

    ! define name of fields which are going to be exchanged
    !
    field_name(CPF_U10M) = "zonal_wind_in_10m"
    collection_size(CPF_U10M) = 1

    field_name(CPF_V10M) = "meridional_wind_in_10m"
    collection_size(CPF_V10M) = 1

    field_name(CPF_FR_SEAICE) = "fraction_of_ocean_covered_by_sea_ice"
    collection_size(CPF_FR_SEAICE) = 1

    field_name(CPF_Z0) = "roughness_length"
    collection_size(CPF_Z0) = 1

    jg = 1
    patch_horz => p_patch(jg)

    comp_name = TRIM(get_my_process_name())

    ! Inform the coupler about what we are
    CALL yac_fdef_comp ( comp_name = TRIM(comp_name), & !in
      &                  comp_id   = comp_id )          !out
    comp_ids(1) = comp_id

    ! Print the YAC version
    CALL message('Running ICON-waves in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)

    CALL yac_fdef_datetime ( start_datetime = TRIM(startdatestring), & !in
      &                      end_datetime   = TRIM(stopdatestring)   ) !in

    ! Announce one grid (patch) to the coupler
    grid_name = "icon_waves_grid"

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes in rad.

    nblks = MAX(patch_horz%nblks_c,patch_horz%nblks_v)

    nbr_vertices_per_cell = 3

    ALLOCATE(buffer_lon(nproma*nblks), STAT=ist)
    IF (ist /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed for buffer_lon!')
    !
    ALLOCATE(buffer_lat(nproma*nblks), STAT=ist)
    IF (ist /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed for buffer_lat!')
    !
    ALLOCATE(buffer_c(nbr_vertices_per_cell,nproma*nblks),STAT=ist)
    IF (ist /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed for buffer_c!')

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, jv, nlen, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_v
      IF (jb /= patch_horz%nblks_v) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_v
      END IF
      DO jv = 1, nlen
        nn = (jb-1)*nproma+jv
        buffer_lon(nn) = patch_horz%verts%vertex(jv,jb)%lon
        buffer_lat(nn) = patch_horz%verts%vertex(jv,jb)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(jb, jc, nlen, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      IF (jb /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO jc = 1, nlen
        nn = (jb-1)*nproma+jc
        buffer_c(1,nn) = (patch_horz%cells%vertex_blk(jc,jb,1)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,1)
        buffer_c(2,nn) = (patch_horz%cells%vertex_blk(jc,jb,2)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,2)
        buffer_c(3,nn) = (patch_horz%cells%vertex_blk(jc,jb,3)-1)*nproma + &
                          patch_horz%cells%vertex_idx(jc,jb,3)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! Definition of unstructured horizontal grid
    CALL yac_fdef_grid(                                   &
      & grid_name             = TRIM(grid_name),          & !in
      & nbr_vertices          = patch_horz%n_patch_verts, & !in
      & nbr_cells             = patch_horz%n_patch_cells, & !in
      & nbr_vertices_per_cell = nbr_vertices_per_cell,    & !in
      & x_vertices            = buffer_lon,               & !in
      & y_vertices            = buffer_lat,               & !in
      & cell_to_vertex        = buffer_c,                 & !in
      & grid_id               = grid_id)                    !out

    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nlen, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      IF (jb /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO jc = 1, nlen
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = patch_horz%cells%center(jc,jb)%lon
        buffer_lat(nn) = patch_horz%cells%center(jc,jb)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour interpolation)
    CALL yac_fdef_points (                     &
      & grid_id    = grid_id,                  & !in
      & nbr_points = patch_horz%n_patch_cells, & !in
      & location   = YAC_LOCATION_CELL,        & !in
      & x_points   = buffer_lon,               & !in
      & y_points   = buffer_lat,               & !in
      & point_id   = cell_point_ids(1) )         !out

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c, STAT=ist)
    IF (ist /= SUCCESS)  THEN
      CALL finish (routine, 'DEALLOCATE failed for buffer_lon, buffer_lat, buffer_c!')
    ENDIF

    ! set global ids for grid cells
    CALL yac_fset_global_index (                               &
      & global_index = patch_horz%cells%decomp_info%glb_index, & !in
      & location     = YAC_LOCATION_CELL,                      & !in
      & grid_id      = grid_id )                                 !in

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c), STAT=ist)
    IF (ist /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed for is_valid!')

    ! TODO
    ! I am not fully sure what this scalar npr_inner_cells is good for, anth whether
    ! we need it for atmo-wave coupling.
    nbr_inner_cells = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jc) REDUCTION(+:nbr_inner_cells) ICON_OMP_RUNTIME_SCHEDULE
    DO jc = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(jc) ) THEN
         is_valid(jc) = .TRUE.
         nbr_inner_cells = nbr_inner_cells + 1
       ELSE
         is_valid(jc) = .FALSE.
       ENDIF
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! set grid core mask
    CALL yac_fset_core_mask (            &
      & is_core  = is_valid,             & !in
      & location = YAC_LOCATION_CELL,    & !in
      & grid_id  = grid_id )               !in

    ! Note that the wave grid consists of water cells, only.
    ! Hence, all cells in ICON-waves grid are valid
!ICON_OMP_PARALLEL_DO PRIVATE(jc) ICON_OMP_RUNTIME_SCHEDULE
    DO jc = 1, patch_horz%nblks_c * nproma
       is_valid(jc) = .TRUE.
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! define cell mask.
    ! Points with is_valid=.FALSE. are masked out
    CALL yac_fdef_mask (                       &
      & grid_id    = grid_id,                  & !in
      & nbr_points = patch_horz%n_patch_cells, & !in
      & location   = YAC_LOCATION_CELL,        & !in
      & is_valid   = is_valid,                 & !in
      & mask_id    = cell_mask_ids(1) )          !out


    ! convert model timestep into ISO string
    CALL timedeltaToString(time_config%tc_dt_model, tc_dt_model_string)
    
    DO jn = 1, no_of_fields
      CALL yac_fdef_field_mask (                      &
        & field_name      = TRIM(field_name(jn)),     & !in
        & component_id    = comp_id,                  & !in
        & point_ids       = cell_point_ids(1),        & !in
        & mask_ids        = cell_mask_ids(1),         & !in
        & num_pointsets   = 1,                        & !in
        & collection_size = collection_size(jn),      & !in
        & timestep        = tc_dt_model_string,       & !in
        & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
        & field_id        = field_id(jn) )              !out
    ENDDO

    DEALLOCATE (is_valid, STAT=ist)
    IF (ist /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed for is_valid!')

    ! End definition of coupling fields and search
    CALL yac_fenddef ( )

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_wave_coupling

END MODULE mo_wave_coupling_frame
