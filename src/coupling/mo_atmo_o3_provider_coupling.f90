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

! @brief Interface between atmosphere and the O3 provider, through a coupler

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_o3_provider_coupling

  USE mo_kind,            ONLY: wp
  USE mo_model_domain,    ONLY: t_patch
  USE mo_exception,       ONLY: finish
  USE mo_aes_rad_config,  ONLY: aes_rad_config
  USE mo_coupling_config, ONLY: is_coupled_to_o3
  USE mo_coupling_utils,  ONLY: cpl_def_field, cpl_get_field, &
                                cpl_get_field_collection_size
  USE mo_sync,            ONLY: SYNC_C, sync_patch_array

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_o3_provider_coupling' ! Output of module for debug

  PUBLIC :: construct_atmo_o3_provider_coupling_post_sync, &
            couple_atmo_to_o3_provider, nplev_o3_provider

  INTEGER :: field_id_o3
  INTEGER :: nplev_o3_provider

  REAL(wp), ALLOCATABLE :: recv_buf(:,:)

CONTAINS

  !>
  !! Registers fields required for the coupling between atmo and
  !! o3 provider
  !!
  !! This subroutine is called from construct_atmo_coupling.
  !!
  SUBROUTINE construct_atmo_o3_provider_coupling_post_sync( &
    comp_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    INTEGER, PARAMETER :: jg = 1

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_atmo_o3_provider_coupling_post_sync'

    IF  (.NOT. is_coupled_to_o3() .OR. &
         (aes_rad_config(jg)%irad_o3 /= 5 .AND. &
          aes_rad_config(jg)%irad_o3 /= 6)) &
      CALL finish(routine, "invalid configuration")

    nplev_o3_provider = &
      cpl_get_field_collection_size( &
        routine, "o3_provider", "o3_grid", "o3")

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "o3", nplev_o3_provider, field_id_o3)

  END SUBROUTINE construct_atmo_o3_provider_coupling_post_sync

  !>
  !! Receives fields from the o3 provider in the atmosphere model
  !!
  SUBROUTINE couple_atmo_to_o3_provider(p_patch, vmr2mmr_o3, o3_plev)

    TYPE(t_patch), INTENT(in) :: p_patch
    REAL(wp), INTENT(in) :: vmr2mmr_o3
    REAL(wp), TARGET, INTENT(inout) :: o3_plev(:,:,:,:)

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':couple_atmo_to_o3_provider'

    LOGICAL :: received_data

    IF ( .NOT. ALLOCATED(recv_buf) ) THEN
      ALLOCATE(recv_buf(p_patch%n_patch_cells, nplev_o3_provider))
      recv_buf = 0.0_wp
    END IF

    CALL cpl_get_field( &
      routine, field_id_o3, 'o3', o3_plev(:,:,:,1), recv_buf, &
      scale_factor=vmr2mmr_o3, first_get=.TRUE., received_data=received_data)
    IF (received_data) &
      CALL sync_patch_array( &
        SYNC_C, p_patch, o3_plev(:,:,:,1), opt_varname='o3')

  END SUBROUTINE couple_atmo_to_o3_provider

END MODULE mo_atmo_o3_provider_coupling
