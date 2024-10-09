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

! contains extended lhs-matrix-generator type
! provides the primal flip flop lhs for the "mass matrix inversion" - solve

MODULE mo_primal_flip_flop_lhs

  USE mo_exception, ONLY: finish
  USE mo_kind, ONLY: wp
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_solve_aux, ONLY: solve_invalid
  USE mo_model_domain, ONLY: t_patch_3d, t_patch
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_2D_per_level
  USE mo_ocean_types, ONLY: t_operator_coeff
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_primal_flip_flop_lhs

  TYPE, EXTENDS(t_lhs_agen) :: t_primal_flip_flop_lhs
    PRIVATE
    TYPE(t_patch_3d), POINTER :: patch_3d => NULL()
    TYPE(t_patch), POINTER :: patch_2d => NULL()
    TYPE(t_operator_coeff), POINTER :: op_coeffs
    INTEGER :: jk
  CONTAINS
    PROCEDURE :: lhs_wp => lhs_primal_flip_flop_wp ! deferred override
    PROCEDURE :: construct => lhs_primal_flip_flop_construct ! deferred override
    PROCEDURE :: lhs_matrix_shortcut => lhs_primal_flip_flop_matrix_shortcut
  END TYPE t_primal_flip_flop_lhs

CONTAINS

! init generator object
  SUBROUTINE lhs_primal_flip_flop_construct(this, patch_3d, op_coeffs, jk, lacc)
    CLASS(t_primal_flip_flop_lhs), INTENT(INOUT) :: this
    TYPE(t_patch_3d), TARGET, INTENT(IN) :: patch_3d
    TYPE(t_operator_coeff), TARGET, INTENT(IN) :: op_coeffs
    INTEGER, INTENT(IN) :: jk
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish("lhs_primal_flip_flop_wp()", "OpenACC version not implemented yet")
#endif

    this%is_const = .false.
    this%patch_3d => patch_3d
    this%patch_2d => patch_3d%p_patch_2d(1)
    this%op_coeffs => op_coeffs
    this%use_shortcut = .false.
    this%jk = jk
    this%is_init = .true.
  END SUBROUTINE lhs_primal_flip_flop_construct

! apply operator to vector x
  SUBROUTINE lhs_primal_flip_flop_wp(this, x, ax, lacc)
    CLASS(t_primal_flip_flop_lhs), INTENT(INOUT) :: this
    REAL(KIND=wp), INTENT(IN) :: x(:,:)
    REAL(KIND=wp), INTENT(INOUT) ::ax(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish("lhs_primal_flip_flop_wp()", "OpenACC version not implemented yet")
#endif

    IF(.NOT.this%is_init .OR. this%jk .EQ. solve_invalid) &
      CALL finish("lhs_primal_flip_flop_wp()", "not correctly initialized")
    CALL map_edges2edges_viacell_2D_per_level(this%patch_3d, &
      & x(:,:), this%op_coeffs, ax(:,:), this%jk)
  END SUBROUTINE lhs_primal_flip_flop_wp

! clear object internals

  SUBROUTINE lhs_primal_flip_flop_matrix_shortcut(this, idx, blk, coeff, lacc)
    CLASS(t_primal_flip_flop_lhs), INTENT(INOUT) :: this
    INTEGER, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: idx, blk
    REAL(KIND=wp), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: coeff
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    CALL finish("lhs_primal_flip_flop_matrix_shortcut", &
      & "not implemented -- go away!")
  END SUBROUTINE lhs_primal_flip_flop_matrix_shortcut

END MODULE mo_primal_flip_flop_lhs
