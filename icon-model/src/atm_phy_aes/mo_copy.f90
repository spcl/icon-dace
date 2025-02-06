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

! Copy to data with given index ranges

MODULE mo_copy

  USE mo_kind, ONLY: wp


  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: copy

  INTERFACE copy
     ! integer  0d --> nd
     MODULE PROCEDURE copy_i0d_to_i0d
     MODULE PROCEDURE copy_i0d_to_i1d
     MODULE PROCEDURE copy_i0d_to_i2d
     ! real(wp) 0d --> nd
     MODULE PROCEDURE copy_r0d_to_r0d
     MODULE PROCEDURE copy_r0d_to_r1d
     MODULE PROCEDURE copy_r0d_to_r2d
     ! real(wp) nd --> nd
     MODULE PROCEDURE copy_r1d_to_r1d
     MODULE PROCEDURE copy_r2d_to_r2d
  END INTERFACE copy

CONTAINS

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_i0d_to_i0d(i0d_in,i0d_out)

    ! Arguments
    !
    INTEGER , INTENT(in)  :: i0d_in
    INTEGER , INTENT(out) :: i0d_out

    i0d_out = i0d_in

  END SUBROUTINE copy_i0d_to_i0d

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_i0d_to_i1d(j1s,j1e, i0d_in,i1d_out)

    ! Arguments
    !
    INTEGER , INTENT(in)    :: j1s,j1e
    INTEGER , INTENT(in)    :: i0d_in
    INTEGER , INTENT(inout) :: i1d_out(:)

    ! INTENT(inout) is used for i1d_out so that the content
    ! outside of the index range j1s:j1e is conserved.

    ! Local variables
    !
    INTEGER :: j1

    !$ACC DATA PRESENT(i1d_out)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(j1s, j1e, i0d_in)
    !$ACC LOOP GANG VECTOR
    DO j1 = j1s,j1e
       i1d_out(j1) = i0d_in
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE copy_i0d_to_i1d

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_i0d_to_i2d(j1s,j1e, j2s,j2e, i0d_in,i2d_out)

    ! Arguments
    !
    INTEGER , INTENT(in)    :: j1s,j1e
    INTEGER , INTENT(in)    :: j2s,j2e
    INTEGER , INTENT(in)    :: i0d_in
    INTEGER , INTENT(inout) :: i2d_out(:,:)

    ! INTENT(inout) is used for i2d_out so that the content
    ! outside of the index range j1s:j1e and j2s:j2e is conserved.

    ! Local variables
    !
    INTEGER :: j1, j2

    !$ACC DATA PRESENT(i2d_out)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(j1s, j1e, j2s, j2e, i0d_in)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO j2 = j2s,j2e
       DO j1 = j1s,j1e
          i2d_out(j1,j2) = i0d_in
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE copy_i0d_to_i2d

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_r0d_to_r0d(r0d_in,r0d_out)

    ! Arguments
    !
    REAL(wp), INTENT(in)  :: r0d_in
    REAL(wp), INTENT(out) :: r0d_out

    r0d_out = r0d_in

  END SUBROUTINE copy_r0d_to_r0d

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_r0d_to_r1d(j1s,j1e, r0d_in,r1d_out)

    ! Arguments
    !
    INTEGER , INTENT(in)    :: j1s,j1e
    REAL(wp), INTENT(in)    :: r0d_in
    REAL(wp), INTENT(inout) :: r1d_out(:)

    ! INTENT(inout) is used for r1d_out so that the content
    ! outside of the index range j1s:j1e is conserved.

    ! Local variables
    !
    INTEGER :: j1

    !$ACC DATA PRESENT(r1d_out)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(j1s, j1e, r0d_in)
    !$ACC LOOP GANG VECTOR
    DO j1 = j1s,j1e
       r1d_out(j1) = r0d_in
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE copy_r0d_to_r1d

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_r0d_to_r2d(j1s,j1e, j2s,j2e, r0d_in,r2d_out)

    ! Arguments
    !
    INTEGER , INTENT(in)    :: j1s,j1e
    INTEGER , INTENT(in)    :: j2s,j2e
    REAL(wp), INTENT(in)    :: r0d_in
    REAL(wp), INTENT(inout) :: r2d_out(:,:)

    ! INTENT(inout) is used for r2d_out so that the content
    ! outside of the index range j1s:j1e and j2s:j2e is conserved.

    ! Local variables
    !
    INTEGER :: j1,j2

    !$ACC DATA PRESENT(r2d_out)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(j1s, j1e, j2s, j2e, r0d_in)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO j2 = j2s,j2e
       DO j1 = j1s,j1e
          r2d_out(j1,j2) = r0d_in
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE copy_r0d_to_r2d

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_r1d_to_r1d(j1s,j1e, r1d_in,r1d_out)

    ! Arguments
    !
    INTEGER , INTENT(in)    :: j1s,j1e
    REAL(wp), INTENT(in)    :: r1d_in (:)
    REAL(wp), INTENT(inout) :: r1d_out(:)

    ! INTENT(inout) is used for r1d_out so that the content
    ! outside of the index range j1s:j1e is conserved.

    ! Local variables
    !
    INTEGER :: j1

    !$ACC DATA PRESENT(r1d_in, r1d_out)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(j1s, j1e)
    !$ACC LOOP GANG VECTOR
    DO j1 = j1s,j1e
       r1d_out(j1) = r1d_in(j1)
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE copy_r1d_to_r1d

  !-----------------------------------------------------------------------------

  SUBROUTINE copy_r2d_to_r2d(j1s,j1e, j2s,j2e, r2d_in,r2d_out)

    ! Arguments
    !
    INTEGER , INTENT(in)    :: j1s,j1e
    INTEGER , INTENT(in)    :: j2s,j2e
    REAL(wp), INTENT(in)    :: r2d_in (:,:)
    REAL(wp), INTENT(inout) :: r2d_out(:,:)

    ! INTENT(inout) is used for r2d_out so that the content
    ! outside of the index range j1s:j1e and j2s:j2e is conserved.

    ! Local variables
    !
    INTEGER :: j1,j2

    !$ACC DATA PRESENT(r2d_in, r2d_out)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(j1s, j1e, j2s, j2e)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO j2 = j2s,j2e
       DO j1 = j1s,j1e
          r2d_out(j1,j2) = r2d_in(j1,j2)
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE copy_r2d_to_r2d

  !-----------------------------------------------------------------------------

END MODULE mo_copy
