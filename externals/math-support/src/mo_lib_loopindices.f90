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
!! This module contains subroutines needed to determine the start and end
!! indices of do loops for a given patch and block index.
!!

MODULE mo_lib_loopindices

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_indices_c_lib, get_indices_e_lib, get_indices_v_lib

CONTAINS

!-------------------------------------------------------------------------
!
!
!
!>
!! Computes the start and end indices of do loops for cell-based variables.
!!
!!
  SUBROUTINE get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, i_blk, i_startblk, i_endblk, &
                               i_startidx_out, i_endidx_out)

    INTEGER, INTENT(IN) :: i_startidx_in ! Start index as input
    INTEGER, INTENT(IN) :: i_endidx_in ! End index as input
    INTEGER, INTENT(IN) :: nproma ! inner loop length/vector length
    INTEGER, INTENT(IN) :: i_blk ! Current block (variable jb in do loops)
    INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
    INTEGER, INTENT(IN) :: i_endblk ! End block of do loop

    INTEGER, INTENT(OUT) :: i_startidx_out, i_endidx_out ! Start and end indices (jc loop)

    IF (i_blk == i_startblk) THEN
      i_startidx_out = MAX(1, i_startidx_in)
      i_endidx_out = nproma
      IF (i_blk == i_endblk) i_endidx_out = i_endidx_in
    ELSE IF (i_blk == i_endblk) THEN
      i_startidx_out = 1
      i_endidx_out = i_endidx_in
    ELSE
      i_startidx_out = 1
      i_endidx_out = nproma
    END IF

  END SUBROUTINE get_indices_c_lib

!-------------------------------------------------------------------------
!
!
!
!>
!! Computes the start and end indices of do loops for edge-based variables.
!!
!!
  SUBROUTINE get_indices_e_lib(i_startidx_in, i_endidx_in, nproma, i_blk, i_startblk, i_endblk, &
                               i_startidx_out, i_endidx_out)

    INTEGER, INTENT(IN) :: i_startidx_in ! Start index as input
    INTEGER, INTENT(IN) :: i_endidx_in ! End index as input
    INTEGER, INTENT(IN) :: nproma ! inner loop length/vector length
    INTEGER, INTENT(IN) :: i_blk ! Current block (variable jb in do loops)
    INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
    INTEGER, INTENT(IN) :: i_endblk ! End block of do loop

    INTEGER, INTENT(OUT) :: i_startidx_out, i_endidx_out ! Start and end indices (je loop)

    i_startidx_out = MERGE(1, &
      &                MAX(1, i_startidx_in), &
      &                i_blk /= i_startblk)

    i_endidx_out = MERGE(nproma, i_endidx_in, &
      &                i_blk /= i_endblk)

  END SUBROUTINE get_indices_e_lib

!-------------------------------------------------------------------------
!
!
!
!>
!! Computes the start and end indices of do loops for cell-based variables.
!!
!!
  SUBROUTINE get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, i_blk, i_startblk, i_endblk, &
                               i_startidx_out, i_endidx_out)

    INTEGER, INTENT(IN) :: i_startidx_in ! Start index as input
    INTEGER, INTENT(IN) :: i_endidx_in ! End index as input
    INTEGER, INTENT(IN) :: nproma ! inner loop length/vector length
    INTEGER, INTENT(IN) :: i_blk ! Current block (variable jb in do loops)
    INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
    INTEGER, INTENT(IN) :: i_endblk ! End block of do loop

    INTEGER, INTENT(OUT) :: i_startidx_out, i_endidx_out ! Start and end indices (jv loop)

    IF (i_blk == i_startblk) THEN
      i_startidx_out = i_startidx_in
      i_endidx_out = nproma
      IF (i_blk == i_endblk) i_endidx_out = i_endidx_in
    ELSE IF (i_blk == i_endblk) THEN
      i_startidx_out = 1
      i_endidx_out = i_endidx_in
    ELSE
      i_startidx_out = 1
      i_endidx_out = nproma
    END IF

  END SUBROUTINE get_indices_v_lib

END MODULE mo_lib_loopindices

