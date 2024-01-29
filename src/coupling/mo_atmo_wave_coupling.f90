! Interface between atmosphere physics and the ocean surface waves, through a coupler
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

MODULE mo_atmo_wave_coupling

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: warning, message_text, message
  USE mo_model_domain,     ONLY: t_patch
  USE mo_fortran_tools,    ONLY: assert_acc_host_only
  USE mo_run_config,       ONLY: ltimer
  USE mo_timer,            ONLY: timer_start, timer_stop, &
    &                            timer_coupling_put, timer_coupling_get, &
    &                            timer_coupling_1stget
  USE mo_sync,             ONLY: SYNC_C, sync_patch_array
  USE mo_yac_finterface,   ONLY: yac_fput, yac_fget, yac_dble_ptr, &
    &                            yac_fdef_field, &
    &                            YAC_ACTION_COUPLING, YAC_ACTION_OUT_OF_BOUND, &
    &                            YAC_TIME_UNIT_ISO_FORMAT
  USE mo_coupling,         ONLY: lyac_very_1st_get

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_atmo_wave_coupling, couple_atmo_to_wave

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_atmo_wave_coupling'

  INTEGER :: field_id_u10m
  INTEGER :: field_id_v10m
  INTEGER :: field_id_fr_seaice
  INTEGER :: field_id_z0

CONTAINS

  SUBROUTINE def_field( &
    comp_id, cell_point_id, timestepstring, &
    field_name, collection_size, field_id)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    INTEGER, INTENT(IN) :: collection_size
    INTEGER, INTENT(OUT) :: field_id

    CALL yac_fdef_field (                           &
      & field_name      = TRIM(field_name),         & !in
      & component_id    = comp_id,                  & !in
      & point_ids       = (/cell_point_id/),        & !in
      & num_pointsets   = 1,                        & !in
      & collection_size = collection_size,          & !in
      & timestep        = timestepstring,           & !in
      & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
      & field_id        = field_id )              !out

  END SUBROUTINE def_field

  !>
  !! Registers fields required for the coupling between atmosphere and wave
  !!
  !! This subroutine is called from constrcut_atmo_coupling.
  !!
  SUBROUTINE construct_atmo_wave_coupling( &
    comp_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "zonal_wind_in_10m", 1, field_id_u10m)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "meridional_wind_in_10m", 1, field_id_v10m)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "fraction_of_ocean_covered_by_sea_ice", 1, field_id_fr_seaice)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "roughness_length", 1, field_id_z0)

  END SUBROUTINE construct_atmo_wave_coupling

  !>
  !! Exchange fields between atmosphere and wave model
  !!
  !! Send fields to the wave model:
  !!   "meridional_wind_in_10m"
  !!   "zonal_wind_in_10m"
  !!   "fraction_of_ocean_covered_by_sea_ice"
  !!
  !! Receive fields from the wave model:
  !!   "roughness_length"
  !!
  !! This subroutine is called from nwp_nh_interface.
  !!
  SUBROUTINE couple_atmo_to_wave(p_patch, u10m, v10m, fr_seaice, z0_waves, lacc)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':couple_atmo_to_wave'

    TYPE(t_patch),                INTENT(IN)   :: p_patch
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)   :: u10m(:,:)      !< zonal wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)   :: v10m(:,:)      !< meridional wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)   :: fr_seaice(:,:) !< fraction_of_ocean_covered_by_sea_ice [1]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT):: z0_waves(:,:)  !< surface roughness length [m]
    LOGICAL,  OPTIONAL,           INTENT(IN)   :: lacc           ! If true, use openacc

    LOGICAL :: write_coupler_restart
    INTEGER :: info, ierror
    TYPE(yac_dble_ptr) :: yac_ptr(1,1)


    CALL assert_acc_host_only('couple_atmo_to_wave', lacc)

    write_coupler_restart = .FALSE.

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from atmosphere to wave
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !------------------------------------------------
    !  Send 10m zonal wind u10m
    !  'zonal_wind_in_10m'
    !------------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
    !
    yac_ptr(1,1)%p(1:p_patch%n_patch_cells) => u10m(:,:)
    !
    CALL yac_fput (                      &
      & field_id        = field_id_u10m, &
      & nbr_pointsets   = 1,             &
      & collection_size = 1,             &
      & send_field      = yac_ptr,       &
      & info            = info,          &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - U10'
      CALL warning(routine, message_text)
    END IF

    ! ----------------------------------------------
    !  Send 10m meridional wind v10m
    !  'meridional_wind_in_10m'
    ! ----------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
    !
    yac_ptr(1,1)%p(1:p_patch%n_patch_cells) => v10m(:,:)
    !
    CALL yac_fput (                      &
      & field_id        = field_id_v10m, &
      & nbr_pointsets   = 1,             &
      & collection_size = 1,             &
      & send_field      = yac_ptr,       &
      & info            = info,          &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - V10'
      CALL warning(routine, message_text)
    END IF

    ! ------------------------------------------------------------------
    !  Send fraction of sea ice
    !  'fraction_of_ocean_covered_by_sea_ice'
    ! ------------------------------------------------------------------
    !
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
    !
    yac_ptr(1,1)%p(1:p_patch%n_patch_cells) => fr_seaice(:,:)
    !
    CALL yac_fput (                           &
      & field_id        = field_id_fr_seaice, &
      & nbr_pointsets   = 1,                  &
      & collection_size = 1,                  &
      & send_field      = yac_ptr,            &
      & info            = info,               &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - fr_seaice'
      CALL warning(routine, message_text)
    END IF

    IF (write_coupler_restart) THEN
      WRITE(message_text, '(a,3i3,a)') 'YAC says it is put for restart - U10, V10, fr_seaice'
      CALL message(routine, message_text)
    ENDIF


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from wave to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    ! --------------------------------------------
    !  Receive roughness length z0 from the wave model
    !  'roughness_length'
    ! --------------------------------------------
    !
    IF (lyac_very_1st_get) THEN
      IF (ltimer) CALL timer_start(timer_coupling_1stget)
    ELSE
      IF (ltimer) CALL timer_start(timer_coupling_get)
    ENDIF
    !
    yac_ptr(1,1)%p(1:p_patch%n_patch_cells) => z0_waves(:,:)
    !
    CALL yac_fget(                      &
      & field_id        = field_id_z0,  &
      & collection_size = 1,            &
      & recv_field      = yac_ptr(:,1), &
      & info            = info,         &
      & ierror          = ierror )
    !
    IF (lyac_very_1st_get) THEN
      IF (ltimer) CALL timer_stop(timer_coupling_1stget)
      lyac_very_1st_get = .FALSE.
    ELSE
      IF (ltimer) CALL timer_stop(timer_coupling_get)
    ENDIF
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says it is get for restart - Z0'
      CALL message(routine, message_text)
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fget called after end of run - Z0'
      CALL warning(routine, message_text)
    ENDIF

    ! halo synchronization for fields recieved from the atmosphere
    !
    CALL sync_patch_array(SYNC_C, p_patch, z0_waves(:,:), opt_varname="z0")

  END SUBROUTINE couple_atmo_to_wave

END MODULE mo_atmo_wave_coupling
