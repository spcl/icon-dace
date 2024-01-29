! Interface between ocean surface waves and atmosphere, through a coupler
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

MODULE mo_wave_atmo_coupling

  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: warning, message_text, message
  USE mo_model_domain,             ONLY: t_patch
  USE mo_run_config,               ONLY: ltimer
  USE mo_timer,                    ONLY: timer_start, timer_stop, &
    &                                    timer_coupling_put, timer_coupling_get, &
    &                                    timer_coupling_1stget
  USE mo_sync,                     ONLY: SYNC_C, sync_patch_array
  USE mo_coupling,                 ONLY: lyac_very_1st_get
  USE mo_wave_coupling_frame,      ONLY: field_id, &
    &                                    collection_size, CPF_U10M, CPF_V10M, CPF_FR_SEAICE, &
    &                                    CPF_Z0
  USE mo_yac_finterface,           ONLY: yac_fput, yac_fget, yac_dble_ptr,       &
    &                                    YAC_ACTION_COUPLING,                    &
    &                                    YAC_ACTION_OUT_OF_BOUND

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: couple_wave_to_atmo

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_atmo_coupling'

CONTAINS

  !>
  !! Exchange fields between the wave model and the atmosphere model
  !!
  !! Send fields to atmosphere:
  !!   "roughness_length"
  !!
  !! Receive fields from atmosphere:
  !!   "meridional_wind_in_10m"
  !!   "zonal_wind_in_10m"
  !!   "fraction_of_ocean_covered_by_sea_ice"
  !!
  !! This subroutine is called from perform_wave_stepping.
  !!
  SUBROUTINE couple_wave_to_atmo(p_patch, z0, u10m, v10m, sea_ice_c)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':couple_wave_to_atmo'

    TYPE(t_patch),                INTENT(IN)    :: p_patch
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: z0(:,:)        !< surface roughness length [m]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: u10m(:,:)      !< zonal wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: v10m(:,:)      !< meridional wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: sea_ice_c(:,:) !< fraction_of_ocean_covered_by_sea_ice [1]

    INTEGER :: info, ierror

    TYPE(yac_dble_ptr):: yac_ptr(1,SIZE(field_id))


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from waves to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    ! --------------------------------------------
    !  Send roughness length z0 to the atmosphere
    !  field_name(CPF_Z0) = 'roughness_length'
    ! --------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
    !
    yac_ptr(1,CPF_Z0)%p(1:p_patch%n_patch_cells) => z0(:,:)
    !
    CALL yac_fput (                                              &
      & field_id        = field_id(CPF_Z0),                      &
      & nbr_pointsets   = SIZE(yac_ptr,1),                       &
      & collection_size = collection_size(CPF_Z0),               &
      & send_field      = yac_ptr(:,CPF_Z0:CPF_Z0),              &
      & info            = info,                                  &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says it is put for restart - id=',CPF_Z0,', Z0'
      CALL message(routine, message_text)
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - id=',CPF_Z0,', Z0'
      CALL warning(routine, message_text)
    END IF


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !   Receive fields from atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    ! -----------------------------------------
    !  Receive 10m zonal wind u10m
    !  field_id(CPF_U10M) = 'zonal_wind_in_10m'
    ! -----------------------------------------
    !
    IF (lyac_very_1st_get) THEN
      IF (ltimer) CALL timer_start(timer_coupling_1stget)
    ELSE
      IF (ltimer) CALL timer_start(timer_coupling_get)
    ENDIF
    !
    yac_ptr(1,CPF_U10M)%p(1:p_patch%n_patch_cells) => u10m(:,:)
    !
    CALL yac_fget(                                   &
      & field_id        = field_id(CPF_U10M),        &
      & collection_size = collection_size(CPF_U10M), &
      & recv_field      = yac_ptr(:,CPF_U10M),       &
      & info            = info,                      &
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
      WRITE(message_text, '(a,i3,a)') 'YAC says it is get for restart - id=',CPF_U10M,', U10'
      CALL message(routine, message_text)
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - id=',CPF_U10M,', U10'
      CALL warning(routine, message_text)
    ENDIF


    ! ----------------------------------------------
    !  Receive 10m meridional wind v10m
    !  field_id(CPF_V10M) = 'meridional_wind_in_10m'
    ! ----------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)
    !
    yac_ptr(1,CPF_V10M)%p(1:p_patch%n_patch_cells) => v10m(:,:)
    !
    CALL yac_fget(                                   &
      & field_id        = field_id(CPF_V10M),        &
      & collection_size = collection_size(CPF_V10M), &
      & recv_field      = yac_ptr(:,CPF_V10M),       &
      & info            = info,                      &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says it is get for restart - id=',CPF_V10M,', V10'
      CALL message(routine, message_text)
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - id=',CPF_V10M,', V10'
      CALL warning(routine, message_text)
    ENDIF

    ! ------------------------------------------------------------------
    !  Receive fraction of sea ice
    !  field_id(CPF_FR_SEAICE) = 'fraction_of_ocean_covered_by_sea_ice'
    ! ------------------------------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)
    !
    yac_ptr(1,CPF_FR_SEAICE)%p(1:p_patch%n_patch_cells) => sea_ice_c(:,:)
    !
    CALL yac_fget(                                        &
      & field_id        = field_id(CPF_FR_SEAICE),        &
      & collection_size = collection_size(CPF_FR_SEAICE), &
      & recv_field      = yac_ptr(:,CPF_FR_SEAICE),       &
      & info            = info,                           &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says it is get for restart - id=',CPF_FR_SEAICE,', fr_seaice'
      CALL message(routine, message_text)
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - id=',CPF_FR_SEAICE,', fr_seaice'
      CALL warning(routine, message_text)
    ENDIF


    ! halo synchronization for fields which have been received from the atmosphere
    !
    CALL sync_patch_array(SYNC_C, p_patch, u10m(:,:), opt_varname="u10m")
    CALL sync_patch_array(SYNC_C, p_patch, v10m(:,:), opt_varname="v10m")
    CALL sync_patch_array(SYNC_C, p_patch, sea_ice_c(:,:), opt_varname="sea_ice_c")

  END SUBROUTINE couple_wave_to_atmo

END MODULE mo_wave_atmo_coupling
