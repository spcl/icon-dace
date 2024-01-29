!
! mo_art_clipping
! This module provides clipping routines.
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

MODULE mo_art_clipping
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_fortran_tools,                 ONLY: set_acc_host_or_device
   
  IMPLICIT NONE
 
  PRIVATE
 
  PUBLIC :: art_clip_gt
  PUBLIC :: art_clip_lt


INTERFACE art_clip_gt
  MODULE PROCEDURE art_clip_gt_r0d
  MODULE PROCEDURE art_clip_gt_r1d
  MODULE PROCEDURE art_clip_gt_r2d
  MODULE PROCEDURE art_clip_gt_r3d
  MODULE PROCEDURE art_clip_gt_r4d
  MODULE PROCEDURE art_clip_gt_r5d
END INTERFACE art_clip_gt

INTERFACE art_clip_lt
  MODULE PROCEDURE art_clip_lt_r0d
  MODULE PROCEDURE art_clip_lt_r1d
  MODULE PROCEDURE art_clip_lt_r2d
  MODULE PROCEDURE art_clip_lt_r3d
  MODULE PROCEDURE art_clip_lt_r4d
  MODULE PROCEDURE art_clip_lt_r5d
END INTERFACE art_clip_lt

    
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_gt_r0d(r0d,value)
!<
! SUBROUTINE art_clip_gt_r0d
! Clips a 0d field to a given value if larger than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r0d
  REAL(wp),INTENT(in)    :: value
    
  IF (r0d > value) THEN
    r0d = value
  END IF
END SUBROUTINE art_clip_gt_r0d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_gt_r1d(r1d,value)
!<
! SUBROUTINE art_clip_gt_r1d
! Clips a 1d field to a given value if larger than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r1d(:)
  REAL(wp),INTENT(in)    :: value

  WHERE (r1d > value)
    r1d = value
  END WHERE
END SUBROUTINE art_clip_gt_r1d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_gt_r2d(r2d,value)
!<
! SUBROUTINE art_clip_gt_r2d
! Clips a 2d field to a given value if larger than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r2d(:,:)
  REAL(wp),INTENT(in)    :: value

  WHERE (r2d > value)
    r2d = value
  END WHERE
END SUBROUTINE art_clip_gt_r2d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_gt_r3d(r3d,value)
!<
! SUBROUTINE art_clip_gt_r3d
! Clips a 3d field to a given value if larger than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r3d(:,:,:)
  REAL(wp),INTENT(in)    :: value

  WHERE (r3d > value)
    r3d = value
  END WHERE
END SUBROUTINE art_clip_gt_r3d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_gt_r4d(r4d,value)
!<
! SUBROUTINE art_clip_gt_r4d
! Clips a 4d field to a given value if larger than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r4d(:,:,:,:)
  REAL(wp),INTENT(in)    :: value

  WHERE (r4d > value)
    r4d = value
  END WHERE
END SUBROUTINE art_clip_gt_r4d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_gt_r5d(r5d,value)
!<
! SUBROUTINE art_clip_gt_r5d
! Clips a 5d field to a given value if larger than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r5d(:,:,:,:,:)
  REAL(wp),INTENT(in)    :: value

  WHERE (r5d > value)
    r5d = value
  END WHERE
END SUBROUTINE art_clip_gt_r5d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_lt_r0d(r0d,value,lacc)
!<
! SUBROUTINE art_clip_gt_r0d
! Clips a 0d field to a given value if smaller than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r0d
  REAL(wp),INTENT(in)    :: value

  LOGICAL, OPTIONAL, INTENT(in) :: lacc

#ifdef _OPENACC
  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  !$ACC SERIAL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  IF (r0d < value) THEN
    r0d = value
  END IF
  !$ACC END SERIAL
#else
  IF (r0d < value) THEN
    r0d = value
  END IF
#endif

END SUBROUTINE art_clip_lt_r0d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_lt_r1d(r1d,value,lacc)
!<
! SUBROUTINE art_clip_gt_r1d
! Clips a 1d field to a given value if smaller than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r1d(:)
  REAL(wp),INTENT(in)    :: value

  LOGICAL, OPTIONAL, INTENT(in) :: lacc

#ifdef _OPENACC
  INTEGER :: i
  INTEGER :: i_size

  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  i_size = SIZE(r1d,1)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i = 1, i_size
    IF (r1d(i) < value) THEN
      r1d(i) = value
    END IF
  END DO
  !$ACC END PARALLEL
#else
  WHERE (r1d < value)
    r1d = value
  END WHERE
#endif

END SUBROUTINE art_clip_lt_r1d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_lt_r2d(r2d,value,lacc)
!<
! SUBROUTINE art_clip_gt_r2d
! Clips a 2d field to a given value if smaller than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r2d(:,:)
  REAL(wp),INTENT(in)    :: value

  LOGICAL, OPTIONAL, INTENT(in) :: lacc

#ifdef _OPENACC
  INTEGER :: i, j
  INTEGER :: i_size, j_size

  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  i_size = SIZE(r2d,1)
  j_size = SIZE(r2d,2)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO j = 1, j_size
    DO i = 1, i_size
      IF (r2d(i,j) < value) THEN
        r2d(i,j) = value
      END IF
    END DO
  END DO
  !$ACC END PARALLEL
#else
  WHERE (r2d < value)
    r2d = value
  END WHERE
#endif

END SUBROUTINE art_clip_lt_r2d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_lt_r3d(r3d,value,lacc)
!<
! SUBROUTINE art_clip_gt_r3d
! Clips a 3d field to a given value if smaller than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r3d(:,:,:)
  REAL(wp),INTENT(in)    :: value

  LOGICAL, OPTIONAL, INTENT(in) :: lacc

#ifdef _OPENACC
  INTEGER :: i, j, k
  INTEGER :: i_size, j_size, k_size

  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  i_size = SIZE(r3d,1)
  j_size = SIZE(r3d,2)
  k_size = SIZE(r3d,3)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(3)
  DO k = 1, k_size
    DO j = 1, j_size
      DO i = 1, i_size
        IF (r3d(i,j,k) < value) THEN
          r3d(i,j,k) = value
        END IF
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
#else
  WHERE (r3d < value)
    r3d = value
  END WHERE
#endif

END SUBROUTINE art_clip_lt_r3d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_lt_r4d(r4d,value,lacc)
!<
! SUBROUTINE art_clip_gt_r4d
! Clips a 4d field to a given value if smaller than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r4d(:,:,:,:)
  REAL(wp),INTENT(in)    :: value

  LOGICAL, OPTIONAL, INTENT(in) :: lacc

#ifdef _OPENACC
  INTEGER :: i, j, k, m
  INTEGER :: i_size, j_size, k_size, m_size

  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  i_size = SIZE(r4d,1)
  j_size = SIZE(r4d,2)
  k_size = SIZE(r4d,3)
  m_size = SIZE(r4d,4)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(4)
  DO m = 1, m_size
    DO k = 1, k_size
      DO j = 1, j_size
        DO i = 1, i_size
          IF (r4d(i,j,k,m) < value) THEN
            r4d(i,j,k,m) = value
          END IF
        END DO
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
#else
  WHERE (r4d < value)
    r4d = value
  END WHERE
#endif

END SUBROUTINE art_clip_lt_r4d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clip_lt_r5d(r5d,value,lacc)
!<
! SUBROUTINE art_clip_gt_r5d
! Clips a 5d field to a given value if smaller than that value
! Based on: -
! Part of Module: mo_art_clipping
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp),INTENT(inout) :: r5d(:,:,:,:,:)
  REAL(wp),INTENT(in)    :: value

  LOGICAL, OPTIONAL, INTENT(in) :: lacc

#ifdef _OPENACC
  INTEGER :: i, j, k, m, n
  INTEGER :: i_size, j_size, k_size, m_size, n_size

  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  i_size = SIZE(r5d,1)
  j_size = SIZE(r5d,2)
  k_size = SIZE(r5d,3)
  m_size = SIZE(r5d,4)
  n_size = SIZE(r5d,5)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(5)
  DO n = 1, n_size
    DO m = 1, m_size
      DO k = 1, k_size
        DO j = 1, j_size
          DO i = 1, i_size
            IF (r5d(i,j,k,m,n) < value) THEN
              r5d(i,j,k,m,n) = value
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
#else
  WHERE (r5d < value)
    r5d = value
  END WHERE
#endif

END SUBROUTINE art_clip_lt_r5d
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_clipping
