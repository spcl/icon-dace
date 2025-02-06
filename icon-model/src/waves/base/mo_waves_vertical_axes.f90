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

! Specification of vertical axes for the wave model

MODULE mo_waves_vertical_axes

  USE mo_kind,                              ONLY: dp
  USE mo_zaxis_type,                        ONLY: ZA_SURFACE, ZA_HEIGHT_10M, ZA_reference
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxisList
  USE mo_name_list_output_zaxes,            ONLY: single_level_axis, vertical_axis
  USE mo_level_selection_types,             ONLY: t_level_selection
  USE mo_run_config,                        ONLY: num_lev

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_zaxes_waves

CONTAINS


  SUBROUTINE setup_zaxes_waves(verticalAxisList, level_selection, log_patch_id)
    TYPE(t_verticalAxisList), INTENT(INOUT)       :: verticalAxisList
    TYPE(t_level_selection),  INTENT(IN), POINTER :: level_selection  ! in general non-associated for waves
    INTEGER,                  INTENT(IN)          :: log_patch_id

    ! local
    INTEGER :: k

    ! --------------------------------------------------------------------------------------
    ! Definitions for single levels --------------------------------------------------------
    ! --------------------------------------------------------------------------------------

    ! surface level
    CALL verticalAxisList%append(single_level_axis(ZA_surface))

    ! Specified height level above ground: 10m
    CALL verticalAxisList%append(single_level_axis(ZA_height_10m, opt_level_value=10._dp))

    ! --------------------------------------------------------------------------------------
    ! Dummy vertical axis with a single full level -----------------------------------------
    ! --------------------------------------------------------------------------------------
    ! REFERENCE
    CALL verticalAxisList%append(vertical_axis(ZA_reference, num_lev(log_patch_id),                       &
      &                           levels           = (/ ( REAL(k,dp),   k=1,num_lev(log_patch_id)+1 ) /), &
      &                           level_selection  = level_selection,                                     &
      &                           opt_set_bounds   = .TRUE. )                                             )

  END SUBROUTINE setup_zaxes_waves

END MODULE mo_waves_vertical_axes

