! @brief Interface between atmosphere and the aero provider, through a coupler
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

MODULE mo_atmo_aero_provider_coupling

  USE mo_kind,            ONLY: wp
  USE mo_model_domain,    ONLY: t_patch
  USE mo_exception,       ONLY: finish
  USE mo_aes_rad_config,  ONLY: aes_rad_config
  USE mo_coupling_config, ONLY: is_coupled_to_aero
  USE mo_coupling_utils,  ONLY: cpl_def_field, cpl_get_field, &
                                cpl_get_field_collection_size

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_aero_provider_coupling' ! Output of module for debug

  PUBLIC :: construct_atmo_aero_provider_coupling_post_sync, &
            couple_atmo_to_aero_provider

  INTEGER :: field_id_aod_c_f
  INTEGER :: field_id_ssa_c_f
  INTEGER :: field_id_z_km_aer_c_mo
  INTEGER :: field_id_aod_c_s
  INTEGER :: field_id_ssa_c_s
  INTEGER :: field_id_asy_c_s
  INTEGER :: field_id_aod_f_s
  INTEGER :: field_id_ssa_f_s
  INTEGER :: field_id_asy_f_s
  INTEGER :: field_id_z_km_aer_f_mo

  INTEGER, PARAMETER   :: nb_lw = 16, nb_sw = 14, lev_clim = 40

  REAL(wp), ALLOCATABLE :: recv_buf(:,:)

CONTAINS

  !>
  !! Registers fields required for the coupling between atmo and
  !! aero provider
  !!
  !! This subroutine is called from construct_atmo_coupling.
  !!
  SUBROUTINE construct_atmo_aero_provider_coupling_post_sync( &
    comp_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    INTEGER, PARAMETER :: jg = 1

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_atmo_aero_provider_coupling_post_sync'

    IF (.NOT. is_coupled_to_aero() .OR. &
        (aes_rad_config(jg)%irad_aero /= 12 .AND. &
         aes_rad_config(jg)%irad_aero /= 13 .AND. &
         aes_rad_config(jg)%irad_aero /= 15 .AND. &
         aes_rad_config(jg)%irad_aero /= 18 .AND. &
         aes_rad_config(jg)%irad_aero /= 19)) &
      CALL finish(routine, "invalid configuration")

    IF ( nb_lw /= &
         cpl_get_field_collection_size( &
          routine, "aero_provider", "aero_grid", "aod_lw_b16_coa") ) &
      CALL finish(routine, 'inconsistent number of lw bands')
    IF ( nb_sw /= &
         cpl_get_field_collection_size( &
           routine, "aero_provider", "aero_grid", "aod_sw_b14_coa") ) &
      CALL finish(routine, 'inconsistent number of sw bands')
    IF ( lev_clim /= &
        cpl_get_field_collection_size( &
          routine, "aero_provider", "aero_grid", "aer_lw_b16_coa") ) &
      CALL finish(routine, 'inconsistent number of level')

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aod_lw_b16_coa", nb_lw, field_id_aod_c_f)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "ssa_lw_b16_coa", nb_lw, field_id_ssa_c_f)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aer_lw_b16_coa", lev_clim, field_id_z_km_aer_c_mo)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aod_sw_b14_coa", nb_sw, field_id_aod_c_s)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "ssa_sw_b14_coa", nb_sw, field_id_ssa_c_s)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "asy_sw_b14_coa", nb_sw, field_id_asy_c_s)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aod_sw_b14_fin", nb_sw, field_id_aod_f_s)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "ssa_sw_b14_fin", nb_sw, field_id_ssa_f_s)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "asy_sw_b14_fin", nb_sw, field_id_asy_f_s)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aer_sw_b14_fin", lev_clim, field_id_z_km_aer_f_mo)

  END SUBROUTINE construct_atmo_aero_provider_coupling_post_sync

  !>
  !! Receives fields from the aero provider in the atmosphere model
  !!
  SUBROUTINE couple_atmo_to_aero_provider( &
    p_patch, aod_f_s,  ssa_f_s, asy_f_s, aod_c_s, ssa_c_s, asy_c_s, &
    aod_c_f, ssa_c_f, z_km_aer_f_mo, z_km_aer_c_mo)

    TYPE(t_patch), INTENT(in) :: p_patch

    ! Fine mode SW
    REAL(wp), TARGET, INTENT(inout) :: aod_f_s(:,:,:,:)
    REAL(wp), TARGET, INTENT(inout) :: ssa_f_s(:,:,:,:)
    REAL(wp), TARGET, INTENT(inout) :: asy_f_s(:,:,:,:)
    ! Coarse mode SW
    REAL(wp), TARGET, INTENT(inout) :: aod_c_s(:,:,:,:)
    REAL(wp), TARGET, INTENT(inout) :: ssa_c_s(:,:,:,:)
    REAL(wp), TARGET, INTENT(inout) :: asy_c_s(:,:,:,:)
    ! Coarse mode LW
    REAL(wp), TARGET, INTENT(inout) :: aod_c_f(:,:,:,:)
    REAL(wp), TARGET, INTENT(inout) :: ssa_c_f(:,:,:,:)
    ! Fine mode height profiles
    REAL(wp), TARGET, INTENT(inout) :: z_km_aer_f_mo(:,:,:,:)
    ! Coarse mode height profiles
    REAL(wp), TARGET, INTENT(inout) :: z_km_aer_c_mo(:,:,:,:)

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':couple_atmo_to_aero_provider'

    INTEGER :: num_cells

    num_cells = p_patch%n_patch_cells

    IF ( .NOT. ALLOCATED(recv_buf) ) THEN
      ALLOCATE(recv_buf(num_cells, MAX(lev_clim, nb_sw, nb_lw)))
      recv_buf = 0.0_wp
    END IF

    ! aod_c_f       = aod_lw_b16_coa -> paer_tau_lw_vr ( lw band )
    ! ssa_c_f       = ssa_lw_b16_coa -> zs_i           ( lw band )
    ! z_km_aer_c_mo = aer_lw_b16_coa -> zq_aod_c       ( level )

    ! aod_lw_b16_coa -> aod_c_f ( band )
    CALL cpl_get_field( &
      routine, field_id_aod_c_f, 'aod_c_f', &
      aod_c_f(:,1:nb_lw,:,1), recv_buf, first_get=.TRUE.)

    !ssa_lw_b16_coa -> ssa_c_f ( band )
    CALL cpl_get_field( &
      routine, field_id_ssa_c_f, 'ssa_c_f', &
      ssa_c_f(:,1:nb_lw,:,1), recv_buf)

    ! aer_lw_b16_coa -> z_km_aer_c_mo ( level )
    CALL cpl_get_field( &
      routine, field_id_z_km_aer_c_mo, 'z_km_aer_c_mo', &
      z_km_aer_c_mo(:,1:lev_clim,:,1), recv_buf)

    ! aod_c_s       = aod_sw_b14_coa -> zt_c  ( sw band )
    ! ssa_c_s       = ssa_sw_b14_coa -> zs_c  ( sw band )
    ! asy_c_s       = asy_sw_b14_coa -> zg_c   (sw band )
    ! z_km_aer_c_mo = aer_sw_b14_coa -> duplicated, take from aer_lw_b16_coa

    ! aod_lw_b16_coa -> aod_c_f ( band )
    CALL cpl_get_field( &
      routine, field_id_aod_c_s, 'aod_c_f', &
      aod_c_s(:,1:nb_sw,:,1), recv_buf)

    !ssa_sw_b14_coa -> ssa_c_s ( band )
    CALL cpl_get_field( &
      routine, field_id_ssa_c_s, 'ssa_c_s', &
      ssa_c_s(:,1:nb_sw,:,1), recv_buf)

    !asy_sw_b14_coa -> asy_c_s ( band )
    CALL cpl_get_field( &
      routine, field_id_asy_c_s, 'asy_c_s', &
      asy_c_s(:,1:nb_sw,:,1), recv_buf)

    ! aod_f_s       = aod_sw_b14_fin -> zt_f     ( band )
    ! ssa_f_s       = ssa_sw_b14_fin -> zs_f     ( band )
    ! asy_f_s       = asy_sw_b14_fin -> zg_f     ( band )
    ! z_km_aer_f_mo = aer_sw_b14_fin -> zq_aod_f ( level )

    ! aod_sw_b14_fin -> aod_f_s ( band )
    CALL cpl_get_field( &
      routine, field_id_aod_f_s, 'aod_f_s', &
      aod_f_s(:,1:nb_sw,:,1), recv_buf)

    !ssa_sw_b14_fin -> ssa_f_s ( band )
    CALL cpl_get_field( &
      routine, field_id_ssa_f_s, 'ssa_f_s', &
      ssa_f_s(:,1:nb_sw,:,1), recv_buf)

    !asy_sw_b14_fin -> asy_f_s ( band )
    CALL cpl_get_field( &
      routine, field_id_asy_f_s, 'asy_f_s', &
      asy_f_s(:,1:nb_sw,:,1), recv_buf)

    ! aer_sw_b14_fin -> z_km_aer_f_mo ( level )
    CALL cpl_get_field( &
      routine, field_id_z_km_aer_f_mo, 'z_km_aer_f_mo', &
      z_km_aer_f_mo(:,1:lev_clim,:,1), recv_buf)

  END SUBROUTINE couple_atmo_to_aero_provider

END MODULE mo_atmo_aero_provider_coupling
