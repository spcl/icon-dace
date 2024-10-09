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

!>
!! Contains the implementation of interpolation and reconstruction.
!!
!! Contains the implementation of interpolation and reconstruction
!! routines used by the shallow water model, including the RBF
!! reconstruction routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_lib_intp_rbf

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64, &
                                           dp => real64, &
                                           sp => real32
  USE mo_lib_loopindices, ONLY: get_indices_c_lib, get_indices_e_lib, get_indices_v_lib
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rbf_vec_interpol_cell_lib, rbf_interpol_c2grad_lib,     &
           & rbf_vec_interpol_vertex_lib, rbf_vec_interpol_edge_lib

  INTERFACE rbf_vec_interpol_vertex_lib
    MODULE PROCEDURE rbf_vec_interpol_vertex_wp_lib
    MODULE PROCEDURE rbf_vec_interpol_vertex_vp_lib
  END INTERFACE

CONTAINS

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at cell center.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each cell center.
!!
  SUBROUTINE rbf_vec_interpol_cell_lib(p_vn_in, rbf_vec_idx_c, rbf_vec_blk_c, &
    &                                  rbf_vec_coeff_c, p_u_out, p_v_out, &
    &                                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                  slev, elev, nproma, lacc, acc_async)

! !INPUT PARAMETERS

    ! input normal components of (velocity) vectors at edge midpoints
    ! dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: p_vn_in(:, :, :)

    ! index array defining the stencil of surrounding edges for vector rbf interpolation at each cell center
    ! (rbf_vec_dim_c,nproma,nblks_c)
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_idx_c(:, :, :)

    ! ... dito for the blocks
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_blk_c(:, :, :)

    ! array containing the surrounding edges in the stencil for vector rbf interpolation at each cell center
    ! (nproma,nblks_c)
    REAL(wp), TARGET, INTENT(IN) :: rbf_vec_coeff_c(:, :, :, :)

    ! reconstructed x-component (u) of velocity vector, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: p_u_out(:, :, :)

    ! reconstructed y-component (v) of velocity vector, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: p_v_out(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !< async OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: acc_async

! !LOCAL VARIABLES
    INTEGER :: jc, jk, jb ! integer over cells, levels, and blocks,

    INTEGER :: i_startidx ! start index
    INTEGER :: i_endidx ! end index

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:, :, :, :), POINTER :: ptr_coeff

    LOGICAL :: lzacc
!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => rbf_vec_idx_c
    iblk => rbf_vec_blk_c

    ptr_coeff => rbf_vec_coeff_c

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), ICON_OMP_RUNTIME_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          p_u_out(jc, jk, jb) = &
            ptr_coeff(1, 1, jc, jb)*p_vn_in(iidx(1, jc, jb), jk, iblk(1, jc, jb)) + &
            ptr_coeff(2, 1, jc, jb)*p_vn_in(iidx(2, jc, jb), jk, iblk(2, jc, jb)) + &
            ptr_coeff(3, 1, jc, jb)*p_vn_in(iidx(3, jc, jb), jk, iblk(3, jc, jb)) + &
            ptr_coeff(4, 1, jc, jb)*p_vn_in(iidx(4, jc, jb), jk, iblk(4, jc, jb)) + &
            ptr_coeff(5, 1, jc, jb)*p_vn_in(iidx(5, jc, jb), jk, iblk(5, jc, jb)) + &
            ptr_coeff(6, 1, jc, jb)*p_vn_in(iidx(6, jc, jb), jk, iblk(6, jc, jb)) + &
            ptr_coeff(7, 1, jc, jb)*p_vn_in(iidx(7, jc, jb), jk, iblk(7, jc, jb)) + &
            ptr_coeff(8, 1, jc, jb)*p_vn_in(iidx(8, jc, jb), jk, iblk(8, jc, jb)) + &
            ptr_coeff(9, 1, jc, jb)*p_vn_in(iidx(9, jc, jb), jk, iblk(9, jc, jb))
          p_v_out(jc, jk, jb) = &
            ptr_coeff(1, 2, jc, jb)*p_vn_in(iidx(1, jc, jb), jk, iblk(1, jc, jb)) + &
            ptr_coeff(2, 2, jc, jb)*p_vn_in(iidx(2, jc, jb), jk, iblk(2, jc, jb)) + &
            ptr_coeff(3, 2, jc, jb)*p_vn_in(iidx(3, jc, jb), jk, iblk(3, jc, jb)) + &
            ptr_coeff(4, 2, jc, jb)*p_vn_in(iidx(4, jc, jb), jk, iblk(4, jc, jb)) + &
            ptr_coeff(5, 2, jc, jb)*p_vn_in(iidx(5, jc, jb), jk, iblk(5, jc, jb)) + &
            ptr_coeff(6, 2, jc, jb)*p_vn_in(iidx(6, jc, jb), jk, iblk(6, jc, jb)) + &
            ptr_coeff(7, 2, jc, jb)*p_vn_in(iidx(7, jc, jb), jk, iblk(7, jc, jb)) + &
            ptr_coeff(8, 2, jc, jb)*p_vn_in(iidx(8, jc, jb), jk, iblk(8, jc, jb)) + &
            ptr_coeff(9, 2, jc, jb)*p_vn_in(iidx(9, jc, jb), jk, iblk(9, jc, jb))
        END DO
      END DO
      !$ACC END PARALLEL

    END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (PRESENT(acc_async)) THEN
      IF (.NOT. acc_async) THEN
        !$ACC WAIT
      END IF
    ELSE
      !$ACC WAIT
    END IF

  END SUBROUTINE rbf_vec_interpol_cell_lib
!====================================================================================

!====================================================================================
  SUBROUTINE rbf_interpol_c2grad_lib(p_cell_in, rbf_c2grad_idx, rbf_c2grad_blk, &
    &                                rbf_c2grad_coeff, grad_x, grad_y, &
    &                                i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                slev, elev, nproma, lacc)

! !INPUT PARAMETERS
!
    ! input cell-based variable for which gradient at cell center is computed
    ! dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: p_cell_in(:, :, :)
 
    ! index array defining the stencil of surrounding cells for 2D gradient reconstruction at each cell center
    ! (rbf_c2grad_dim,nproma,nblks_c)
    INTEGER, TARGET, INTENT(IN) :: rbf_c2grad_idx(:, :, :)
 
    ! ... dito for the blocks
    INTEGER, TARGET, INTENT(IN) :: rbf_c2grad_blk(:, :, :)
 
    ! array containing the coefficients used for 2D gradient reconstruction at each cell center
    ! (rbf_c2grad_dim,2,nproma,nblks_c)
    REAL(wp), TARGET, INTENT(IN) :: rbf_c2grad_coeff(:, :, :, :)
 
    ! reconstructed zonal (x) component of gradient vector, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: grad_x(:, :, :)
 
    ! reconstructed meridional (y) component of gradient vector, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: grad_y(:, :, :)
 
    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk
 
    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk
 
    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in
 
    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in
 
    ! vertical start level
    INTEGER, INTENT(IN) :: slev
 
    ! vertical end level
    INTEGER, INTENT(IN) :: elev
 
    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma
 
    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

! !LOCAL VARIABLES
    INTEGER :: jc, jk, jb ! integer over cells, levels, and blocks,

    INTEGER :: i_startidx ! start index
    INTEGER :: i_endidx ! end index

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:, :, :, :), POINTER :: ptr_coeff

    LOGICAL :: lzacc
!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => rbf_c2grad_idx
    iblk => rbf_c2grad_blk

    ptr_coeff => rbf_c2grad_coeff

!$OMP PARALLEL

!$ACC DATA PRESENT(p_cell_in, grad_x, grad_y, iidx, iblk, ptr_coeff) IF(lzacc)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), ICON_OMP_RUNTIME_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          grad_x(jc, jk, jb) = &
            ptr_coeff(1, 1, jc, jb)*p_cell_in(jc, jk, jb) + &
            ptr_coeff(2, 1, jc, jb)*p_cell_in(iidx(2, jc, jb), jk, iblk(2, jc, jb)) + &
            ptr_coeff(3, 1, jc, jb)*p_cell_in(iidx(3, jc, jb), jk, iblk(3, jc, jb)) + &
            ptr_coeff(4, 1, jc, jb)*p_cell_in(iidx(4, jc, jb), jk, iblk(4, jc, jb)) + &
            ptr_coeff(5, 1, jc, jb)*p_cell_in(iidx(5, jc, jb), jk, iblk(5, jc, jb)) + &
            ptr_coeff(6, 1, jc, jb)*p_cell_in(iidx(6, jc, jb), jk, iblk(6, jc, jb)) + &
            ptr_coeff(7, 1, jc, jb)*p_cell_in(iidx(7, jc, jb), jk, iblk(7, jc, jb)) + &
            ptr_coeff(8, 1, jc, jb)*p_cell_in(iidx(8, jc, jb), jk, iblk(8, jc, jb)) + &
            ptr_coeff(9, 1, jc, jb)*p_cell_in(iidx(9, jc, jb), jk, iblk(9, jc, jb)) + &
            ptr_coeff(10, 1, jc, jb)*p_cell_in(iidx(10, jc, jb), jk, iblk(10, jc, jb))
          grad_y(jc, jk, jb) = &
            ptr_coeff(1, 2, jc, jb)*p_cell_in(jc, jk, jb) + &
            ptr_coeff(2, 2, jc, jb)*p_cell_in(iidx(2, jc, jb), jk, iblk(2, jc, jb)) + &
            ptr_coeff(3, 2, jc, jb)*p_cell_in(iidx(3, jc, jb), jk, iblk(3, jc, jb)) + &
            ptr_coeff(4, 2, jc, jb)*p_cell_in(iidx(4, jc, jb), jk, iblk(4, jc, jb)) + &
            ptr_coeff(5, 2, jc, jb)*p_cell_in(iidx(5, jc, jb), jk, iblk(5, jc, jb)) + &
            ptr_coeff(6, 2, jc, jb)*p_cell_in(iidx(6, jc, jb), jk, iblk(6, jc, jb)) + &
            ptr_coeff(7, 2, jc, jb)*p_cell_in(iidx(7, jc, jb), jk, iblk(7, jc, jb)) + &
            ptr_coeff(8, 2, jc, jb)*p_cell_in(iidx(8, jc, jb), jk, iblk(8, jc, jb)) + &
            ptr_coeff(9, 2, jc, jb)*p_cell_in(iidx(9, jc, jb), jk, iblk(9, jc, jb)) + &
            ptr_coeff(10, 2, jc, jb)*p_cell_in(iidx(10, jc, jb), jk, iblk(10, jc, jb))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE rbf_interpol_c2grad_lib

!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at triangle vertices.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each vertex.
!!
  SUBROUTINE rbf_vec_interpol_vertex_wp_lib(p_e_in, rbf_vec_idx_v, rbf_vec_blk_v, &
    &                                       rbf_vec_coeff_v, p_u_out, p_v_out, &
    &                                       i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                       slev, elev, nproma, lacc, acc_async)

! !INPUT PARAMETERS
!
    ! input components of velocity or horizontal vorticity vectors at edge midpoints
    ! dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: p_e_in(:, :, :)

    ! index array defining the stencil of surrounding edges for vector rbf interpolation at each triangle vertex
    ! (rbf_vec_dim_v,nproma,nblks_v)
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_idx_v(:, :, :)

    ! ... dito for the blocks
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_blk_v(:, :, :)

    ! coefficients are working precision array containing the coefficients used for vector rbf interpolation 
    ! at each tringle vertex (input is normal component), dim: (rbf_vec_dim_v,2,nproma,nblks_v)
    REAL(wp), TARGET, INTENT(IN) :: rbf_vec_coeff_v(:, :, :, :)

    ! reconstructed x-component (u) of velocity vector, dim: (nproma,nlev,nblks_v)
    REAL(wp), INTENT(INOUT) :: p_u_out(:, :, :)

    ! reconstructed y-component (v) of velocity vector, dim: (nproma,nlev,nblks_v)
    REAL(wp), INTENT(INOUT) :: p_v_out(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !< async OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: acc_async

! !LOCAL VARIABLES

    INTEGER :: jv, jk, jb ! integer over vertices, levels, and blocks,

    INTEGER :: i_startidx ! start index
    INTEGER :: i_endidx ! end index

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:, :, :, :), POINTER :: ptr_coeff

    LOGICAL :: lzacc
!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => rbf_vec_idx_v
    iblk => rbf_vec_blk_v

    ptr_coeff => rbf_vec_coeff_v

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jv), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(lzacc)
#ifdef __LOOP_EXCHANGE
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
#endif

          p_u_out(jv, jk, jb) = &
            ptr_coeff(1, 1, jv, jb)*p_e_in(iidx(1, jv, jb), jk, iblk(1, jv, jb)) + &
            ptr_coeff(2, 1, jv, jb)*p_e_in(iidx(2, jv, jb), jk, iblk(2, jv, jb)) + &
            ptr_coeff(3, 1, jv, jb)*p_e_in(iidx(3, jv, jb), jk, iblk(3, jv, jb)) + &
            ptr_coeff(4, 1, jv, jb)*p_e_in(iidx(4, jv, jb), jk, iblk(4, jv, jb)) + &
            ptr_coeff(5, 1, jv, jb)*p_e_in(iidx(5, jv, jb), jk, iblk(5, jv, jb)) + &
            ptr_coeff(6, 1, jv, jb)*p_e_in(iidx(6, jv, jb), jk, iblk(6, jv, jb))
          p_v_out(jv, jk, jb) = &
            ptr_coeff(1, 2, jv, jb)*p_e_in(iidx(1, jv, jb), jk, iblk(1, jv, jb)) + &
            ptr_coeff(2, 2, jv, jb)*p_e_in(iidx(2, jv, jb), jk, iblk(2, jv, jb)) + &
            ptr_coeff(3, 2, jv, jb)*p_e_in(iidx(3, jv, jb), jk, iblk(3, jv, jb)) + &
            ptr_coeff(4, 2, jv, jb)*p_e_in(iidx(4, jv, jb), jk, iblk(4, jv, jb)) + &
            ptr_coeff(5, 2, jv, jb)*p_e_in(iidx(5, jv, jb), jk, iblk(5, jv, jb)) + &
            ptr_coeff(6, 2, jv, jb)*p_e_in(iidx(6, jv, jb), jk, iblk(6, jv, jb))

        END DO
      END DO

    END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (PRESENT(acc_async)) THEN
      IF (.NOT. acc_async) THEN
        !$ACC WAIT
      END IF
    ELSE
      !$ACC WAIT
    END IF

  END SUBROUTINE rbf_vec_interpol_vertex_wp_lib

! Variant for mixed precision mode (output fields in single precision)
  SUBROUTINE rbf_vec_interpol_vertex_vp_lib(p_e_in, rbf_vec_idx_v, rbf_vec_blk_v, &
    &                                       rbf_vec_coeff_v, p_u_out, p_v_out, &
    &                                       i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                       slev, elev, nproma, lacc, acc_async)

    ! input components of velocity or horizontal vorticity vectors at edge midpoints
    ! dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: p_e_in(:, :, :)

    ! index array defining the stencil of surrounding edges for vector rbf interpolation at each triangle vertex
    ! dim(rbf_vec_dim_v,nproma,nblks_v)
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_idx_v(:, :, :)

    ! ... dito for the blocks
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_blk_v(:, :, :)

    ! coefficients are working precision array containing the coefficients used for vector rbf interpolation
    ! at each triangle vertex (input is normal component), dim(rbf_vec_dim_v,2,nproma,nblks_v)
    REAL(wp), TARGET, INTENT(IN) :: rbf_vec_coeff_v(:, :, :, :)

    ! reconstructed x-component (u) of velocity vector, dim: (nproma,nlev,nblks_v)
    REAL(sp), INTENT(INOUT) :: p_u_out(:, :, :)

    ! reconstructed y-component (v) of velocity vector, dim: (nproma,nlev,nblks_v)
    REAL(sp), INTENT(INOUT) :: p_v_out(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

   ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !< async OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: acc_async

! !LOCAL VARIABLES

    INTEGER :: jv, jk, jb ! integer over vertices, levels, and blocks,

    INTEGER :: i_startidx ! start index
    INTEGER :: i_endidx ! end index

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:, :, :, :), POINTER :: ptr_coeff

    LOGICAL :: lzacc
!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => rbf_vec_idx_v
    iblk => rbf_vec_blk_v

    ptr_coeff => rbf_vec_coeff_v

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jv), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
#endif

          p_u_out(jv, jk, jb) = &
            ptr_coeff(1, 1, jv, jb)*p_e_in(iidx(1, jv, jb), jk, iblk(1, jv, jb)) + &
            ptr_coeff(2, 1, jv, jb)*p_e_in(iidx(2, jv, jb), jk, iblk(2, jv, jb)) + &
            ptr_coeff(3, 1, jv, jb)*p_e_in(iidx(3, jv, jb), jk, iblk(3, jv, jb)) + &
            ptr_coeff(4, 1, jv, jb)*p_e_in(iidx(4, jv, jb), jk, iblk(4, jv, jb)) + &
            ptr_coeff(5, 1, jv, jb)*p_e_in(iidx(5, jv, jb), jk, iblk(5, jv, jb)) + &
            ptr_coeff(6, 1, jv, jb)*p_e_in(iidx(6, jv, jb), jk, iblk(6, jv, jb))
          p_v_out(jv, jk, jb) = &
            ptr_coeff(1, 2, jv, jb)*p_e_in(iidx(1, jv, jb), jk, iblk(1, jv, jb)) + &
            ptr_coeff(2, 2, jv, jb)*p_e_in(iidx(2, jv, jb), jk, iblk(2, jv, jb)) + &
            ptr_coeff(3, 2, jv, jb)*p_e_in(iidx(3, jv, jb), jk, iblk(3, jv, jb)) + &
            ptr_coeff(4, 2, jv, jb)*p_e_in(iidx(4, jv, jb), jk, iblk(4, jv, jb)) + &
            ptr_coeff(5, 2, jv, jb)*p_e_in(iidx(5, jv, jb), jk, iblk(5, jv, jb)) + &
            ptr_coeff(6, 2, jv, jb)*p_e_in(iidx(6, jv, jb), jk, iblk(6, jv, jb))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (PRESENT(acc_async)) THEN
      IF (.NOT. acc_async) THEN
        !$ACC WAIT
      END IF
    ELSE
      !$ACC WAIT
    END IF

  END SUBROUTINE rbf_vec_interpol_vertex_vp_lib

!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at edge midpoints.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each edge.
!!
  SUBROUTINE rbf_vec_interpol_edge_lib(p_vn_in, rbf_vec_idx_e, rbf_vec_blk_e, &
    &                                  rbf_vec_coeff_e, p_vt_out, &
    &                                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                  slev, elev, nproma, lacc, acc_async)

! !INPUT PARAMETERS
!

    ! input normal components of velocity vectors at edge midpoints
    ! dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) :: p_vn_in(:, :, :)

    ! index array defining the stencil of surrounding edges for vector rbf interpolation at each triangle edge
    ! dim: (rbf_vec_dim_e,nproma,nblks_e)
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_idx_e(:, :, :)

    ! ... dito for the blocks
    INTEGER, TARGET, INTENT(IN) :: rbf_vec_blk_e(:, :, :) ! ... dito for the blocks

    ! coefficients are working precision array containing the coefficients used for vector rbf interpolation
    ! of the tangential velocity component (from the surrounding normals) at each triangle edge
    ! dim: (rbf_vec_dim_e,nproma,nblks_e)
    REAL(wp), TARGET, INTENT(IN) :: rbf_vec_coeff_e(:, :, :)

    ! reconstructed tangential velocity component, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(INOUT) :: p_vt_out(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !< async OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: acc_async

    INTEGER :: je, jk, jb ! integer over edges, levels, and blocks,

    INTEGER :: i_startidx ! start index
    INTEGER :: i_endidx ! end index

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:, :, :), POINTER :: ptr_coeff

    LOGICAL :: lzacc
!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => rbf_vec_idx_e
    iblk => rbf_vec_blk_e

    ptr_coeff => rbf_vec_coeff_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif

          p_vt_out(je, jk, jb) = &
            ptr_coeff(1, je, jb)*p_vn_in(iidx(1, je, jb), jk, iblk(1, je, jb)) + &
            ptr_coeff(2, je, jb)*p_vn_in(iidx(2, je, jb), jk, iblk(2, je, jb)) + &
            ptr_coeff(3, je, jb)*p_vn_in(iidx(3, je, jb), jk, iblk(3, je, jb)) + &
            ptr_coeff(4, je, jb)*p_vn_in(iidx(4, je, jb), jk, iblk(4, je, jb))

        END DO
      END DO
      !$ACC END PARALLEL
    END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (PRESENT(acc_async)) THEN
      IF (.NOT. acc_async) THEN
        !$ACC WAIT
      END IF
    ELSE
      !$ACC WAIT
    END IF

  END SUBROUTINE rbf_vec_interpol_edge_lib

END MODULE mo_lib_intp_rbf
