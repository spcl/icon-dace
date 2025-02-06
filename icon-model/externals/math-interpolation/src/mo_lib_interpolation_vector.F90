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

MODULE mo_lib_interpolation_vector
  !-------------------------------------------------------------------------
  !
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64, &
                                           sp => real32
  USE mo_lib_loopindices, ONLY: get_indices_c_lib

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: edges2cells_vector_lib

#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp
#else
  INTEGER, PARAMETER :: vp = wp
#endif

CONTAINS

!------------------------------------------------------------------------
!
!>
!!  Bilinear interpolation of normal and tangential velocity components.
!!
!!  Bilinear interpolation of normal and tangential velocity components
!!  at the edges to u and v at the cells
!!  Works only for triangles (bilinear interpolation weights are not implemented
!!  for hexagons)
!!
  SUBROUTINE edges2cells_vector_lib(p_vn_in, p_vt_in, cell_edge_idx, cell_edge_blk, &
    &                               e_bln_c_u, e_bln_c_v, p_u_out, p_v_out, &
    &                               i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                               slev, elev, nproma)

    ! normal velocity component at edges, dim: (nproma,nlev,nblks_e)
    REAL(wp), INTENT(IN) ::  p_vn_in(:, :, :)

    ! (reconstructed) tangential velocity component at edges, dim: (nproma,nlev,nblks_e)
    REAL(vp), INTENT(IN) ::  p_vt_in(:, :, :)

    ! line indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_idx(:, :, :)

    ! block indices of edges of triangles, dim: (nproma,nblks_c, 3)
    INTEGER, TARGET, INTENT(IN) :: cell_edge_blk(:, :, :)

    ! coefficient for bilinear interpolation from edges to cells for vector components
    ! (input: v_t, v_n, output: u)
    REAL(wp), INTENT(IN) :: e_bln_c_u(:, :, :)

    ! coefficient for bilinear interpolation from edges to cells for vector components
    ! (input: v_t, v_n, output: v)
    REAL(wp), INTENT(IN) :: e_bln_c_v(:, :, :)

    ! cell based output fields: u and v, dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) :: p_u_out(:, :, :), p_v_out(:, :, :)

    ! start_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startblk

    ! end_block needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endblk

    ! start_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_startidx_in

    ! end_index needed for get_indices_c_lib
    INTEGER, INTENT(IN) :: i_endidx_in

    ! vertical start level
    INTEGER, INTENT(IN) ::  slev

    ! vertical end level
    INTEGER, INTENT(IN) ::  elev

    ! inner loop length/vector length
    INTEGER, INTENT(IN) :: nproma

    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx

    INTEGER, DIMENSION(:, :, :), POINTER :: iidx, iblk

!-------------------------------------------------------------------------

    iidx => cell_edge_idx
    iblk => cell_edge_blk

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=6
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          p_u_out(jc, jk, jb) = &
            e_bln_c_u(jc, 1, jb)*p_vn_in(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            e_bln_c_u(jc, 2, jb)*p_vt_in(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            e_bln_c_u(jc, 3, jb)*p_vn_in(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            e_bln_c_u(jc, 4, jb)*p_vt_in(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            e_bln_c_u(jc, 5, jb)*p_vn_in(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) + &
            e_bln_c_u(jc, 6, jb)*p_vt_in(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

          p_v_out(jc, jk, jb) = &
            e_bln_c_v(jc, 1, jb)*p_vn_in(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            e_bln_c_v(jc, 2, jb)*p_vt_in(iidx(jc, jb, 1), jk, iblk(jc, jb, 1)) + &
            e_bln_c_v(jc, 3, jb)*p_vn_in(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            e_bln_c_v(jc, 4, jb)*p_vt_in(iidx(jc, jb, 2), jk, iblk(jc, jb, 2)) + &
            e_bln_c_v(jc, 5, jb)*p_vn_in(iidx(jc, jb, 3), jk, iblk(jc, jb, 3)) + &
            e_bln_c_v(jc, 6, jb)*p_vt_in(iidx(jc, jb, 3), jk, iblk(jc, jb, 3))

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE edges2cells_vector_lib

END MODULE mo_lib_interpolation_vector
