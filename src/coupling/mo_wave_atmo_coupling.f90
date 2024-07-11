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

  USE mo_kind,           ONLY: wp
  USE mo_model_domain,   ONLY: t_patch
  USE mo_coupling_utils, ONLY: cpl_def_field, cpl_put_field, cpl_get_field
  USE mo_sync,           ONLY: SYNC_C, sync_patch_array

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_wave_atmo_coupling, couple_wave_to_atmo

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_atmo_coupling'

  INTEGER :: field_id_u10m
  INTEGER :: field_id_v10m
  INTEGER :: field_id_fr_seaice
  INTEGER :: field_id_z0

CONTAINS

  !>
  !! Registers fields required for the coupling between wave and atmosphere
  !!
  !! This subroutine is called from construct_wave_coupling.
  !!
  SUBROUTINE construct_wave_atmo_coupling( &
    comp_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "zonal_wind_in_10m", 1, field_id_u10m)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "meridional_wind_in_10m", 1, field_id_v10m)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "fraction_of_ocean_covered_by_sea_ice", 1, field_id_fr_seaice)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "roughness_length", 1, field_id_z0)

  END SUBROUTINE construct_wave_atmo_coupling

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

    LOGICAL :: received_data

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from waves to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    ! --------------------------------------------
    !  Send roughness length z0 to the atmosphere
    !  'roughness_length'
    ! --------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_z0, 'Z0', p_patch%n_patch_cells, z0)

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !   Receive fields from atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    ! -----------------------------------------
    !  Receive 10m zonal wind u10m
    !  'zonal_wind_in_10m'
    ! -----------------------------------------
    !
    CALL cpl_get_field( &
      routine, field_id_u10m, 'u10m', p_patch%n_patch_cells, u10m, &
      first_get=.TRUE., received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, u10m, opt_varname='u10m')

    ! ----------------------------------------------
    !  Receive 10m meridional wind v10m
    !  'meridional_wind_in_10m'
    ! ----------------------------------------------
    !
    CALL cpl_get_field( &
      routine, field_id_v10m, 'v10m', p_patch%n_patch_cells, v10m, &
      received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, v10m, opt_varname='v10m')

    ! ------------------------------------------------------------------
    !  Receive fraction of sea ice
    !  'fraction_of_ocean_covered_by_sea_ice'
    ! ------------------------------------------------------------------
    !
    CALL cpl_get_field( &
      routine, field_id_fr_seaice, 'sea_ice_c', p_patch%n_patch_cells, &
      sea_ice_c, received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, sea_ice_c, opt_varname='sea_ice_c')

  END SUBROUTINE couple_wave_to_atmo

END MODULE mo_wave_atmo_coupling
