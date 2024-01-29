! abstract type for lhs-matrix generators
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
MODULE mo_ocean_solve_lhs_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_lhs_agen

  TYPE, ABSTRACT :: t_lhs_agen
    LOGICAL :: is_const, use_shortcut
    LOGICAL :: is_init = .false.
  CONTAINS
    PROCEDURE(a_lhs_agen_wp), DEFERRED :: lhs_wp
    PROCEDURE(a_lhs_matrix_shortcut), DEFERRED :: lhs_matrix_shortcut
    GENERIC :: apply => lhs_wp
    GENERIC :: matrix_shortcut => lhs_matrix_shortcut
  END TYPE t_lhs_agen

  ABSTRACT INTERFACE
    SUBROUTINE a_lhs_agen_wp(this, x, ax, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_lhs_agen
      CLASS(t_lhs_agen), INTENT(INOUT) :: this
      REAL(KIND=wp), INTENT(IN) :: x(:,:)
      REAL(KIND=wp), INTENT(OUT) :: ax(:,:)
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_lhs_agen_wp
    SUBROUTINE a_lhs_matrix_shortcut(this, idx, blk, coeff)
      USE mo_kind, ONLY: wp
      IMPORT t_lhs_agen
      CLASS(t_lhs_agen), INTENT(INOUT) :: this
      INTEGER, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: idx, blk
      REAL(KIND=wp), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: coeff
    END SUBROUTINE a_lhs_matrix_shortcut
  END INTERFACE

END MODULE mo_ocean_solve_lhs_type
