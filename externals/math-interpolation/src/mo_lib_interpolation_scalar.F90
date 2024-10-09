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
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_lib_interpolation_scalar
  !-------------------------------------------------------------------------
  !
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64, &
    &                                      dp => real64, &
    &                                      sp => real32
  USE mo_exception, ONLY: finish
  USE mo_lib_loopindices, ONLY: get_indices_c_lib, get_indices_e_lib, get_indices_v_lib
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: verts2edges_scalar_lib
  PUBLIC :: cells2edges_scalar_lib
  PUBLIC :: edges2verts_scalar_lib
  PUBLIC :: edges2cells_scalar_lib
  PUBLIC :: cells2verts_scalar_lib
  PUBLIC :: cells2verts_scalar_ri_lib
  PUBLIC :: verts2cells_scalar_lib
  PUBLIC :: cell_avg_lib

#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp
#else
  INTEGER, PARAMETER :: vp = wp
#endif

  INTERFACE edges2cells_scalar_lib
    MODULE PROCEDURE edges2cells_scalar_dp_lib, edges2cells_scalar_sp_lib
  END INTERFACE edges2cells_scalar_lib

  INTERFACE cells2verts_scalar_lib
    MODULE PROCEDURE cells2verts_scalar_dp_lib, cells2verts_scalar_sp_lib
    MODULE PROCEDURE cells2verts_scalar_sp2dp_lib
  END INTERFACE cells2verts_scalar_lib

CONTAINS

!-----------------------------------------------------------------------
!
!  ! averaging and interpolation routines and
!  ! routines needed to compute the coefficients therein
!
!-----------------------------------------------------------------------
!
!>
!!  Performs  average of scalar fields from vertices to velocity points.
!!
!!  The coefficients are given by c_int.
!!
  SUBROUTINE verts2edges_scalar_lib(p_vertex_in, edge_vertex_idx, edge_vertex_blk, &
    &                               c_int, p_edge_out, i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                               slev, elev, nproma, lacc)

    ! vertex based scalar input field, dim: (nproma,nlev,nblks_v)
    REAL(wp), INTENT(IN) ::  p_vertex_in(:, :, :)

    ! line indices of edge vertices, dim: (nproma,nblks_e, 4)
    INTEGER, TARGET, INTENT(IN) :: edge_vertex_idx(:, :, :)

    ! block indices of edge vertices, dim: (nproma,nblks_e, 4)
    INTEGER, TARGET, INTENT(IN) :: edge_vertex_blk(:, :, :)

    ! interpolation field, dim: (nproma,2,nblks_e)
    REAL(wp), INTENT(IN) ::  c_int(:, :, :)

    ! edge based scalar output field, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(INOUT) :: p_edge_out(:, :, :)

    ! start_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! start index
    INTEGER :: i_startidx

    ! end index
    INTEGER :: i_endidx

    INTEGER :: je, jk, jb

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => edge_vertex_idx
    iblk => edge_vertex_blk

! loop over edges and blocks
!$ACC DATA COPYIN(p_vertex_in, c_int) PRESENT(p_edge_out, iidx, iblk) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif

          p_edge_out(je, jk, jb) =  &
            c_int(je, 1, jb)*p_vertex_in(iidx(je, jb, 1), jk, iblk(je, jb, 1)) + &
            c_int(je, 2, jb)*p_vertex_in(iidx(je, jb, 2), jk, iblk(je, jb, 2))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

  END SUBROUTINE verts2edges_scalar_lib
!------------------------------------------------------------------------
!
!
!>
!!  Computes  average of scalar fields from centers of triangular faces to.
!!
!!  Computes  average of scalar fields from centers of triangular faces to
!!  velocity points.
!!
  SUBROUTINE cells2edges_scalar_lib(p_cell_in, edge_cell_idx, edge_cell_blk, c_int, p_edge_out, &
    &                               i_startblk_in, i_endblk_in, i_startidx_in, i_endidx_in, &
    &                               slev, elev, nproma, patch_id, l_limited_area, lfill_latbc, lacc)

    ! cell based scalar input field, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: p_cell_in(:, :, :)

    ! line indices of adjacent cells, dim: (nproma,nblks_e, 2)
    INTEGER, TARGET, INTENT(IN) :: edge_cell_idx(:, :, :)

    ! block indices of adjacent cells, dim: (nproma,nblks_e, 2)
    INTEGER, TARGET, INTENT(IN) :: edge_cell_blk(:, :, :)

    ! coefficients for linear interpolation, dim: (nproma,2,nblks_e)
    REAL(wp), INTENT(IN) :: c_int(:, :, :)

    ! edge based scalar output field, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(INOUT) :: p_edge_out(:, :, :)

    ! start_block needed for get_indices_e_lib, an array with two elements as it's required
    INTEGER, DIMENSION(2), INTENT(IN) :: i_startblk_in

    ! end_block needed for get_indices_e_lib, an array with two elements as it's required
    INTEGER, DIMENSION(2), INTENT(IN) :: i_endblk_in

    ! start_index needed for get_indices_e_lib, an array with two elements as it's required
    INTEGER, DIMENSION(2), INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_e_lib, an array with two elements as it's required
    INTEGER, DIMENSION(2), INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! domain ID of current domain
    INTEGER, INTENT(IN) :: patch_id

    ! limited area setup where forcing comes in from sides
    LOGICAL, INTENT(IN)  :: l_limited_area

    ! if true, fill lateral nest boundaries
    LOGICAL, INTENT(in) :: lfill_latbc

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: je, jk, jb
    INTEGER :: i_startblk ! start block
    INTEGER :: i_endblk ! end block
    INTEGER :: i_startidx ! start index
    INTEGER :: i_endidx ! end index

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => edge_cell_idx
    iblk => edge_cell_blk

!$ACC DATA NO_CREATE(iidx, iblk, c_int, p_cell_in, p_edge_out) IF(lzacc)

!$OMP PARALLEL PRIVATE(i_startblk, i_endblk)

    IF ((l_limited_area .OR. patch_id > 1) .AND. lfill_latbc) THEN ! Fill outermost nest boundary

      i_startblk = i_startblk_in(1)
      i_endblk = i_endblk_in(1)

! DA: OpenACC needs to collapse loops
#ifndef _OPENACC
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_e_lib(i_startidx_in(1), i_endidx_in(1), nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)

        DO je = i_startidx, i_endidx
          IF (iidx(je, jb, 1) >= 1 .AND. iblk(je, jb, 1) >= 1) THEN
            DO jk = slev, elev
              p_edge_out(je, jk, jb) = p_cell_in(iidx(je, jb, 1), jk, iblk(je, jb, 1))
            END DO
          ELSE IF (iidx(je, jb, 2) >= 1 .AND. iblk(je, jb, 2) >= 1) THEN
            DO jk = slev, elev
              p_edge_out(je, jk, jb) = p_cell_in(iidx(je, jb, 2), jk, iblk(je, jb, 2))
            END DO
          ELSE
            CALL finish('mo_interpolation:cells2edges_scalar',  &
              &          'error in lateral boundary filling')
          END IF
        END DO
      END DO
!$OMP END DO

#else

      DO jb = i_startblk, i_endblk
        CALL get_indices_e_lib(i_startidx_in(1), i_endidx_in(1), nproma, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR TILE(32, 4)
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            IF (iidx(je, jb, 1) >= 1 .AND. iblk(je, jb, 1) >= 1) THEN
              p_edge_out(je, jk, jb) = p_cell_in(iidx(je, jb, 1), jk, iblk(je, jb, 1))
            ELSE IF (iidx(je, jb, 2) >= 1 .AND. iblk(je, jb, 2) >= 1) THEN
              p_edge_out(je, jk, jb) = p_cell_in(iidx(je, jb, 2), jk, iblk(je, jb, 2))
            END IF
          END DO
        END DO
        !$ACC END PARALLEL
      END DO

#endif

    END IF

! Process the remaining grid points for which a real interpolation is possible
    i_startblk = i_startblk_in(2)
    i_endblk = i_endblk_in(2)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e_lib(i_startidx_in(2), i_endidx_in(2), nproma, jb, i_startblk, i_endblk, &
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

          p_edge_out(je, jk, jb) =  &
            c_int(je, 1, jb)*p_cell_in(iidx(je, jb, 1), jk, iblk(je, jb, 1)) +  &
            c_int(je, 2, jb)*p_cell_in(iidx(je, jb, 2), jk, iblk(je, jb, 2))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

!$ACC END DATA

  END SUBROUTINE cells2edges_scalar_lib
!------------------------------------------------------------------------
!
!>
!!  Computes average of scalar fields from velocity points to.
!!
!!  Computes average of scalar fields from velocity points to
!!  centers of dual faces.
!!
  SUBROUTINE edges2verts_scalar_lib(p_edge_in, vert_edge_idx, vert_edge_blk, v_int, p_vert_out, &
    &                               i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                               slev, elev, nproma, lacc)

    ! edge based scalar input field, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) ::  p_edge_in(:, :, :)

    ! line indices of edges around a vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_edge_idx(:, :, :)

    ! block indices of edges around a vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_edge_blk(:, :, :)

    ! coefficients for (area weighted) interpolation, dim: (nproma,cell_type,nblks_v)
    REAL(wp), INTENT(IN) ::  v_int(:, :, :)

    ! vertex based scalar output field, dim: (nproma,nlev,nblks_v)
    REAL(wp), INTENT(INOUT) :: p_vert_out(:, :, :)

    ! start_block needed for get_indices_v_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_v_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_v_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_v_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jv, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => vert_edge_idx
    iblk => vert_edge_blk

!$ACC DATA PRESENT(v_int, p_edge_in, p_vert_out, iidx, iblk) IF(lzacc)

!loop over blocks and verts

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
#endif

          p_vert_out(jv, jk, jb) = &
            v_int(jv, 1, jb)*p_edge_in(iidx(jv, jb, 1), jk, iblk(jv, jb, 1)) + &
            v_int(jv, 2, jb)*p_edge_in(iidx(jv, jb, 2), jk, iblk(jv, jb, 2)) + &
            v_int(jv, 3, jb)*p_edge_in(iidx(jv, jb, 3), jk, iblk(jv, jb, 3)) + &
            v_int(jv, 4, jb)*p_edge_in(iidx(jv, jb, 4), jk, iblk(jv, jb, 4)) + &
            v_int(jv, 5, jb)*p_edge_in(iidx(jv, jb, 5), jk, iblk(jv, jb, 5)) + &
            v_int(jv, 6, jb)*p_edge_in(iidx(jv, jb, 6), jk, iblk(jv, jb, 6))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO !loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

  END SUBROUTINE edges2verts_scalar_lib

!------------------------------------------------------------------------
!
!>
!!  Computes interpolation from edges to cells
!!
!!  Computes interpolation of scalar fields from velocity points to
!!  cell centers via given interpolation weights
!!
  SUBROUTINE edges2cells_scalar_dp_lib(p_edge_in, edge_idx, edge_blk, c_int, p_cell_out, &
    &                                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                  slev, elev, nproma, lacc)

    ! edge based scalar input field, dim: (nproma,nlev,nblks_e)
    REAL(dp), INTENT(IN) ::  p_edge_in(:, :, :)

    ! line indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: edge_blk(:, :, :)

    ! coefficients for (area weighted) interpolation, dim: (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) ::  c_int(:, :, :)

    ! cell based scalar output field, dim: (nproma,nlev,nblks_c)
    REAL(dp), INTENT(INOUT) :: p_cell_out(:, :, :)

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

    INTEGER :: jc, jk, jb
    INTEGER ::  i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => edge_idx
    iblk => edge_blk

!$ACC DATA PRESENT(c_int, p_edge_in, p_cell_out, iidx, iblk) IF(lzacc)

!loop over blocks and cells

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          p_cell_out(jc, jk, jb) = &
            c_int(jc, 1, jb)*p_edge_in(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            c_int(jc, 2, jb)*p_edge_in(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            c_int(jc, 3, jb)*p_edge_in(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO !loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

  END SUBROUTINE edges2cells_scalar_dp_lib
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!>
!!  Computes interpolation from edges to cells
!!
!!  Computes interpolation of scalar fields from velocity points to
!!  cell centers via given interpolation weights
!!
  SUBROUTINE edges2cells_scalar_sp_lib(p_edge_in, cell_edge_idx, cell_edge_blk, c_int, p_cell_out, &
    &                                  i_startblk, i_endblk, i_startidx_in, i_endidx_in,  &
    &                                  slev, elev, nproma, lacc)

    ! edge based scalar input field, dim: (nproma,nlev,nblks_e)
    REAL(sp), INTENT(IN) ::  p_edge_in(:, :, :)

    ! line indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_blk(:, :, :)

    ! coefficients for (area weighted) interpolation, dim: (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) ::  c_int(:, :, :)

    ! cell based scalar output field, dim: (nproma,nlev,nblks_c)
    REAL(sp), INTENT(INOUT) :: p_cell_out(:, :, :)

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

    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_edge_idx
    iblk => cell_edge_blk

!$ACC DATA PRESENT(c_int, p_edge_in, p_cell_out, iidx, iblk) IF(lzacc)

!loop over blocks and cells

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          p_cell_out(jc, jk, jb) = REAL( &
            c_int(jc, 1, jb)*REAL(p_edge_in(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)), wp) + &
            c_int(jc, 2, jb)*REAL(p_edge_in(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)), wp) + &
            c_int(jc, 3, jb)*REAL(p_edge_in(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)), wp),  &
            sp)

        END DO
      END DO
      !$ACC END PARALLEL

    END DO !loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

  END SUBROUTINE edges2cells_scalar_sp_lib
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!>
!!  Computes  average of scalar fields from centers of cells to vertices.
!!
  SUBROUTINE cells2verts_scalar_dp_lib(p_cell_in, vert_cell_idx, vert_cell_blk, c_int, p_vert_out, &
    &                                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                  slev, elev, nproma, lacc, acc_async)

    ! cell based scalar input field, dim: (nproma,nlev,nblks_c)
    REAL(dp), INTENT(IN) :: p_cell_in(:, :, :)

    ! line indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_idx(:, :, :)

    ! block indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_blk(:, :, :)

    ! coefficients for interpolation, dim: (nproma,9-cell_type,nblks_v)
    REAL(wp), INTENT(IN) :: c_int(:, :, :)

    ! vertex based scalar output field, dim: (nproma,nlev,nblks_v)
    REAL(dp), INTENT(INOUT) :: p_vert_out(:, :, :)

    ! start_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_e_lib
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

    INTEGER :: jv, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => vert_cell_idx
    iblk => vert_cell_blk

!$ACC DATA NO_CREATE(iidx, iblk, p_cell_in, c_int, p_vert_out) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
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

          p_vert_out(jv, jk, jb) = &
            c_int(jv, 1, jb)*p_cell_in(iidx(jv, jb, 1), jk, iblk(jv, jb, 1)) + &
            c_int(jv, 2, jb)*p_cell_in(iidx(jv, jb, 2), jk, iblk(jv, jb, 2)) + &
            c_int(jv, 3, jb)*p_cell_in(iidx(jv, jb, 3), jk, iblk(jv, jb, 3)) + &
            c_int(jv, 4, jb)*p_cell_in(iidx(jv, jb, 4), jk, iblk(jv, jb, 4)) + &
            c_int(jv, 5, jb)*p_cell_in(iidx(jv, jb, 5), jk, iblk(jv, jb, 5)) + &
            c_int(jv, 6, jb)*p_cell_in(iidx(jv, jb, 6), jk, iblk(jv, jb, 6))

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

!$ACC END DATA

  END SUBROUTINE cells2verts_scalar_dp_lib
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!>
!!  Computes  average of scalar fields from centers of cells to vertices.
!!
  SUBROUTINE cells2verts_scalar_sp_lib(p_cell_in, vert_cell_idx, vert_cell_blk, c_int, p_vert_out,&
    &                                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                  slev, elev, nproma, lacc, acc_async)

    ! cell based scalar input field, dim: (nproma,nlev,nblks_c)
    REAL(sp), INTENT(IN) :: p_cell_in(:, :, :)

    ! line indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_idx(:, :, :)

    ! block indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_blk(:, :, :)

    ! coefficients for interpolation, dim: (nproma,9-cell_type,nblks_v)
    REAL(wp), INTENT(IN) :: c_int(:, :, :)

    ! vertex based scalar output field, dim: (nproma,nlev,nblks_v)
    REAL(sp), INTENT(INOUT) :: p_vert_out(:, :, :)

    ! start_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_e_lib
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

    INTEGER :: jv, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => vert_cell_idx
    iblk => vert_cell_blk

!$ACC DATA NO_CREATE(p_cell_in, c_int, p_vert_out, iidx, iblk) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=6
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
#endif

          p_vert_out(jv, jk, jb) = &
            c_int(jv, 1, jb)*p_cell_in(iidx(jv, jb, 1), jk, iblk(jv, jb, 1)) + &
            c_int(jv, 2, jb)*p_cell_in(iidx(jv, jb, 2), jk, iblk(jv, jb, 2)) + &
            c_int(jv, 3, jb)*p_cell_in(iidx(jv, jb, 3), jk, iblk(jv, jb, 3)) + &
            c_int(jv, 4, jb)*p_cell_in(iidx(jv, jb, 4), jk, iblk(jv, jb, 4)) + &
            c_int(jv, 5, jb)*p_cell_in(iidx(jv, jb, 5), jk, iblk(jv, jb, 5)) + &
            c_int(jv, 6, jb)*p_cell_in(iidx(jv, jb, 6), jk, iblk(jv, jb, 6))

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

!$ACC END DATA

  END SUBROUTINE cells2verts_scalar_sp_lib
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!>
!!  Computes  average of scalar fields from centers of cells to vertices.
!!
  SUBROUTINE cells2verts_scalar_sp2dp_lib(p_cell_in, vert_cell_idx, vert_cell_blk, c_int, p_vert_out, &
    &                                     i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                     slev, elev, nproma, lacc, acc_async)

    ! cell based scalar input field, dim: (nproma,nlev,nblks_c)
    REAL(sp), INTENT(IN) :: p_cell_in(:, :, :)

    ! line indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_idx(:, :, :)

    ! block indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_blk(:, :, :)

    ! coefficients for interpolation, dim: (nproma,9-cell_type,nblks_v)
    REAL(wp), INTENT(IN) :: c_int(:, :, :)

    ! vertex based scalar output field, dim: (nproma,nlev,nblks_v)
    REAL(dp), INTENT(INOUT) :: p_vert_out(:, :, :)

    ! start_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_e_lib
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

    INTEGER :: jv, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => vert_cell_idx
    iblk => vert_cell_blk

    !$ACC DATA NO_CREATE(p_cell_in, c_int, p_vert_out, iidx, iblk) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=6
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
#endif

          p_vert_out(jv, jk, jb) = &
            c_int(jv, 1, jb)*REAL(p_cell_in(iidx(jv, jb, 1), jk, iblk(jv, jb, 1)), dp) + &
            c_int(jv, 2, jb)*REAL(p_cell_in(iidx(jv, jb, 2), jk, iblk(jv, jb, 2)), dp) + &
            c_int(jv, 3, jb)*REAL(p_cell_in(iidx(jv, jb, 3), jk, iblk(jv, jb, 3)), dp) + &
            c_int(jv, 4, jb)*REAL(p_cell_in(iidx(jv, jb, 4), jk, iblk(jv, jb, 4)), dp) + &
            c_int(jv, 5, jb)*REAL(p_cell_in(iidx(jv, jb, 5), jk, iblk(jv, jb, 5)), dp) + &
            c_int(jv, 6, jb)*REAL(p_cell_in(iidx(jv, jb, 6), jk, iblk(jv, jb, 6)), dp)

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

    !$ACC END DATA

  END SUBROUTINE cells2verts_scalar_sp2dp_lib
!------------------------------------------------------------------------

!>
!!  Same as above, but provides output optionally in single precision and
!!  assumes reversed index order of the output field in loop exchange mode
!!
!!
  SUBROUTINE cells2verts_scalar_ri_lib(p_cell_in, vert_cell_idx, vert_cell_blk, c_int, p_vert_out, &
    &                                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                                  slev, elev, nproma, lacc, acc_async)

    ! cell based scalar input field, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: p_cell_in(:, :, :)

    ! line indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_idx(:, :, :)

    ! block indices of cells around each vertex, dim: (nproma,nblks_v, 6)
    INTEGER, TARGET, INTENT(IN) :: vert_cell_blk(:, :, :)

    ! coefficients for interpolation, dim: (nproma,9-cell_type,nblks_v)
    REAL(wp), INTENT(IN) :: c_int(:, :, :)

    ! vertex based scalar output field, dim: (nproma,nlev,nblks_v)
    REAL(vp), INTENT(INOUT) :: p_vert_out(:, :, :)

    ! start_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_e_lib
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

    INTEGER :: jv, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => vert_cell_idx
    iblk => vert_cell_blk

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jv = i_startidx, i_endidx
        DO jk = slev, elev
          p_vert_out(jk, jv, jb) = &
#else
        !$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx
          p_vert_out(jv, jk, jb) = &
#endif
            c_int(jv, 1, jb)*p_cell_in(iidx(jv, jb, 1), jk, iblk(jv, jb, 1)) + &
            c_int(jv, 2, jb)*p_cell_in(iidx(jv, jb, 2), jk, iblk(jv, jb, 2)) + &
            c_int(jv, 3, jb)*p_cell_in(iidx(jv, jb, 3), jk, iblk(jv, jb, 3)) + &
            c_int(jv, 4, jb)*p_cell_in(iidx(jv, jb, 4), jk, iblk(jv, jb, 4)) + &
            c_int(jv, 5, jb)*p_cell_in(iidx(jv, jb, 5), jk, iblk(jv, jb, 5)) + &
            c_int(jv, 6, jb)*p_cell_in(iidx(jv, jb, 6), jk, iblk(jv, jb, 6))

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

  END SUBROUTINE cells2verts_scalar_ri_lib
!------------------------------------------------------------------------

!
!
!>
!!  Computes  average of scalar fields from vertices to centers of cells.
!!
  SUBROUTINE verts2cells_scalar_lib(p_vert_in, cell_vertex_idx, cell_vertex_blk, &
    &                               c_int, p_cell_out, nblks_c, npromz_c, &
    &                               slev, elev, nproma, lacc)

    ! cell based scalar input field, dim: (nproma,nlev,nblks_v)
    REAL(wp), INTENT(IN) :: p_vert_in(:, :, :)

    ! line indices of vertices of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_vertex_idx(:, :, :)

    ! block indices of vertices of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_vertex_blk(:, :, :)

    ! coefficients for interpolation, dim: (nproma,cell_type,nblks_c)
    REAL(wp), INTENT(IN) :: c_int(:, :, :)

    ! vertex based scalar output field, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: p_cell_out(:, :, :)

    ! number of blocks for the cells: t_patch%nblks_c
    INTEGER, INTENT(IN) :: nblks_c

    ! chunck length in last block for the cells: t_patch%npromz_c
    INTEGER, INTENT(IN) :: npromz_c

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    ! vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jk, jb, jc, nlen

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_vertex_idx
    iblk => cell_vertex_blk

!$ACC DATA PRESENT(p_vert_in, c_int, p_cell_out, iidx, iblk) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      END IF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = 1, nlen
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jc = 1, nlen
#endif

          p_cell_out(jc, jk, jb) = &
            c_int(jc, 1, jb)*p_vert_in(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            c_int(jc, 2, jb)*p_vert_in(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            c_int(jc, 3, jb)*p_vert_in(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

  END SUBROUTINE verts2cells_scalar_lib
!-------------------------------------------------------------------------
!
!
!>
!! Computes the average of a cell-based variable.
!!
!! Computes the average of a cell-based variable
!! over its original location and the neighboring triangles.
!! Version with variable weighting coefficients, computed such that
!! linear horizontal gradients are not aliased into a checkerboard noise
!! input:  lives on centers of triangles
!! output: lives on centers of triangles
!!
  SUBROUTINE cell_avg_lib(psi_c, cell_neighbor_idx, cell_neighbor_blk, avg_coeff, avg_psi_c, &
    &                     i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                     slev, elev, nproma, lacc)

    ! cell based variable before averaging, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: psi_c(:, :, :)

    ! line indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_idx(:, :, :)

    ! block indices of triangles next to each cell, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_neighbor_blk(:, :, :)

    ! averaging coefficients, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(IN) :: avg_coeff(:, :, :)

    ! cell based variable after averaging, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: avg_psi_c(:, :, :)

    ! start_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_e_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) :: slev

    !vertical end level
    INTEGER, INTENT(IN) :: elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    ! if true, use OpenACC
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

    LOGICAL :: lzacc

!-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    iidx => cell_neighbor_idx
    iblk => cell_neighbor_blk

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          !  calculate the weighted average
          !
          avg_psi_c(jc, jk, jb) =  &
            psi_c(jc, jk, jb)*avg_coeff(jc, 1, jb) + &
            psi_c(iidx(jc, jb, 1), jk, iblk(jc, jb, 1))*avg_coeff(jc, 2, jb) + &
            psi_c(iidx(jc, jb, 2), jk, iblk(jc, jb, 2))*avg_coeff(jc, 3, jb) + &
            psi_c(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))*avg_coeff(jc, 4, jb)

        END DO !cell loop
      END DO !vertical levels loop
      !$ACC END PARALLEL

    END DO !block loop

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE cell_avg_lib
!-------------------------------------------------------------------------

END MODULE mo_lib_interpolation_scalar
