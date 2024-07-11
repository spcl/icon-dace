!> Contains utilities for the hydrology process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_hydro_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_amount_in_rootzone

  INTERFACE get_amount_in_rootzone
    ! MODULE PROCEDURE get_amount_in_rootzone_1d
    MODULE PROCEDURE get_amount_in_rootzone_2d
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_util'

CONTAINS

  SUBROUTINE get_amount_in_rootzone_1d(ws, dz, dz_root, ws_root)

    REAL(wp), INTENT(in) :: &
      & ws     (:),       &        !< Water content of soil layers
      & dz     (:),       &        !< Thicknesses of soil layers
      & dz_root(:)                 !< Thicknesses of root layers
    REAL(wp) :: ws_root            !< Water content in root zone

    REAL(wp) ::                       &
      & bottom    , &   !< Cumulative depth of soil down to a layer (bottom depth of layer)
      & bottom_m1  , &  !< Cumulative depth of soil down to previous layer
      & root_depth      !< Total depth of root zone

    INTEGER :: nsoil, is

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_amount_in_rootzone_1d'

    !$ACC ROUTINE SEQ

    nsoil  = SIZE(dz)       ! Number of soil layers
    bottom    = SUM(dz (:))
    root_depth = SUM(dz_root(:))
#ifndef _OPENACC
    IF (root_depth > bottom + EPSILON(1._wp)) THEN
      CALL finish(TRIM(routine), 'Root zone larger than soil depth.')
    ELSE
#endif
      bottom    = 0._wp
      bottom_m1 = 0._wp
      ws_root   = 0._wp
      DO is = 1, nsoil
        bottom    = bottom + dz(is)
        IF (root_depth >= bottom) THEN                ! Root zone extends to below current layer
          ws_root = ws_root + ws(is)
        ELSE IF (root_depth > bottom_m1) THEN         ! Root zone above current layer but below previous layer, i.e.
          ws_root = ws_root + ws(is) * (root_depth - bottom_m1) / dz(is) ! Root z. partially extends into current layer
        END IF
        bottom_m1 = bottom_m1 + dz(is)
      END DO
#ifndef _OPENACC
    END IF
#endif

  END SUBROUTINE get_amount_in_rootzone_1d



  SUBROUTINE get_amount_in_rootzone_2d(ws, dz, dz_root, ws_root)

    REAL(wp), INTENT(in) :: &
      & ws     (:,:),       &        !< Water content of soil layers
      & dz     (:,:),       &        !< Thicknesses of soil layers
      & dz_root(:,:)                 !< Thicknesses of root layers

    REAL(wp), INTENT(inout) :: ws_root(SIZE(ws,1))  !< Water content in root zone

    REAL(wp) ::                       &
      & bottom     (SIZE(ws     ,1)), & !< Cumulative depth of soil down to a layer (bottom depth of layer)
      & bottom_m1  (SIZE(ws_root,1)), & !< Cumulative depth of soil down to previous layer
      & root_depth (SIZE(ws,1))         !< Total depth of root zone

    INTEGER :: nc, ic, nsoil, is

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_amount_in_rootzone_2d'

    nc     = SIZE(ws,1)
    nsoil  = SIZE(dz,2)       ! Number of soil layers

    !$ACC DATA CREATE(bottom, bottom_m1, root_depth)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      root_depth(ic) = 0._wp
      bottom(ic)     = 0._wp
      bottom_m1(ic)  = 0._wp
      ws_root(ic)    = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        root_depth(ic) = root_depth(ic) + dz_root(ic,is)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        bottom(ic)    = bottom(ic) + dz(ic,is)
          IF (root_depth(ic) >= bottom(ic)) THEN              ! Root zone extends to below current layer
            ws_root(ic) = ws_root(ic) + ws(ic,is)
          ELSE IF (root_depth(ic) > bottom_m1(ic)) THEN       ! Root zone above current layer but below previous layer, i.e.
            ! Root zone partially extends into current layer
            ws_root(ic) = ws_root(ic) + ws(ic,is) * (root_depth(ic) - bottom_m1(ic)) / dz(ic,is)
          END IF
          bottom_m1(ic) = bottom_m1(ic) + dz(ic,is)
        END DO
      !$ACC END LOOP
    END DO
    !$ACC END PARALLEL

#ifndef _OPENACC
    DO ic = 1, SIZE(ws, 1)
      IF (root_depth(ic) > bottom(ic) + EPSILON(1._wp)) THEN
        CALL finish(TRIM(routine), 'Root zone larger than soil depth.')
      END IF
    END DO
#endif

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE get_amount_in_rootzone_2d

#endif
END MODULE mo_hydro_util
