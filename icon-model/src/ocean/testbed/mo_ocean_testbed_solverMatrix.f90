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

MODULE mo_ocean_testbed_solverMatrix
  USE mo_kind,                   ONLY: wp
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range
  USE mo_exception,              ONLY: finish
  USE mo_ocean_types, ONLY: t_hydro_ocean_state, t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_mpi,                    ONLY: num_work_procs
  USE mo_parallel_config,        ONLY: nproma
  USE mo_surface_height_lhs, ONLY: t_surface_height_lhs

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: createSolverMatrix

CONTAINS

  SUBROUTINE createSolverMatrix( patch_3d, ocean_state, operators_coefficients)
    TYPE(t_patch_3d), POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in) :: ocean_state(1)
    TYPE(t_operator_coeff), TARGET, INTENT(in) :: operators_coefficients
    TYPE(t_surface_height_lhs) :: lhs
    REAL(wp) :: x(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: ax(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    INTEGER :: column,row
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    INTEGER, PARAMETER :: fileNo = 501
    TYPE(t_solverCoeff_singlePrecision), POINTER :: dummy => NULL()
    TYPE(t_operator_coeff), POINTER :: operators_coefficients_ptr
    CHARACTER(len=*), PARAMETER :: fileName = 'ocean_matrix.txt'
    CHARACTER(len=*), PARAMETER :: method_name = 'createSolverMatrix'

    operators_coefficients_ptr => operators_coefficients
    IF (nproma /= 1) &
      CALL finish(method_name, "nproma should be = 1")
    IF (num_work_procs /= 1) &
      CALL finish(method_name, "mpi processes should be = 1")
    patch_2D           => patch_3d%p_patch_2d(1)
    cells_in_domain    => patch_2D%cells%in_domain
    OPEN(fileNo, FILE=fileName, STATUS='new')
    CALL lhs%construct(patch_3d, ocean_state(1)%p_diag%thick_e, &
        & operators_coefficients_ptr, dummy)
    DO column = cells_in_domain%start_block, cells_in_domain%end_block
      x(:,:) = 0.0_wp
      ax(:,:) = 0.0_wp
      x(1,column) = 1.0_wp
      CALL lhs%apply(x, ax)
      DO row = cells_in_domain%start_block, cells_in_domain%end_block
        IF (ax(1, row) /= 0.0_wp) &
          WRITE(fileNo, *) "(", row, ",", column, ")", ax(1, row)
      ENDDO
    ENDDO
    CLOSE(fileNo)
  END SUBROUTINE createSolverMatrix

END MODULE mo_ocean_testbed_solverMatrix
