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

MODULE mo_aes_phy_cleanup

  USE mtime,                       ONLY: OPERATOR(>)
  USE mo_grid_config,              ONLY: n_dom
  USE mo_aes_phy_memory,           ONLY: destruct_aes_phy_memory
  USE mo_cloud_mig_memory,         ONLY: destruct_cloud_mig_memory
  USE mo_radiation_forcing_memory, ONLY: destruct_radiation_forcing_list
  USE mo_aes_phy_config,           ONLY: aes_phy_tc, dt_zero
  USE mo_turb_vdiff,               ONLY: vdiff_cleanup

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cleanup_aes_phy

CONTAINS
  !>
  !! Top-level routine for the cleanup for ECHAM6 physics.
  !! It calls a series of subroutines to deallocate parameter arrays with
  !! "allocatable" attribute.
  !!
  SUBROUTINE cleanup_aes_phy

    LOGICAL :: lany
    INTEGER :: jg

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (aes_phy_tc(jg)%dt_vdf > dt_zero)
    END DO
    IF (lany) CALL vdiff_cleanup             ! deallocate array "matrix_idx"

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (aes_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) CALL destruct_radiation_forcing_list
   
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (aes_phy_tc(jg)%dt_mig > dt_zero)
    END DO
    IF (lany) CALL destruct_cloud_mig_memory

    CALL destruct_aes_phy_memory

  END SUBROUTINE cleanup_aes_phy
  !-------------

END MODULE mo_aes_phy_cleanup

