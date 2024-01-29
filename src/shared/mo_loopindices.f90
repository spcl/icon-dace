! This module contains subroutines needed to determine the start and end
! indices of do loops for a given patch and block index.
!
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

MODULE mo_loopindices
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

USE mo_model_domain,    ONLY: t_patch
USE mo_impl_constants,  ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_parallel_config,  ONLY: nproma

IMPLICIT NONE

PRIVATE

PUBLIC :: get_indices_c, get_indices_e, get_indices_v

CONTAINS


!-------------------------------------------------------------------------
!
!
!! Computes the start and end indices of do loops for cell-based variables.
!!
SUBROUTINE get_indices_c(p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                         i_endidx, irl_start, opt_rl_end)


  TYPE(t_patch), INTENT(IN) :: p_patch
  INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
  INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
  INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
  INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts

  INTEGER, OPTIONAL, INTENT(IN) :: opt_rl_end ! refin_ctrl level where do loop ends

  INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (jc loop)

  ! Local variables

  INTEGER :: irl_end

  IF (PRESENT(opt_rl_end)) THEN
    irl_end = opt_rl_end
  ELSE
    irl_end = min_rlcell
  ENDIF

  IF (i_blk == i_startblk) THEN
    i_startidx = MAX(1,p_patch%cells%start_index(irl_start))
    i_endidx   = nproma
    IF (i_blk == i_endblk) i_endidx = p_patch%cells%end_index(irl_end)
  ELSE IF (i_blk == i_endblk) THEN
    i_startidx = 1
    i_endidx   = p_patch%cells%end_index(irl_end)
  ELSE
    i_startidx = 1
    i_endidx = nproma
  ENDIF

END SUBROUTINE get_indices_c

!-------------------------------------------------------------------------
!
!! Computes the start and end indices of do loops for edge-based variables.
!!
SUBROUTINE get_indices_e(p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                         i_endidx, irl_start, opt_rl_end)


  TYPE(t_patch), INTENT(IN) :: p_patch
  INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
  INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
  INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
  INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts

  INTEGER, OPTIONAL, INTENT(IN) :: opt_rl_end ! refin_ctrl level where do loop ends

  INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (je loop)


  ! Local variables

  INTEGER :: irl_end

  IF (PRESENT(opt_rl_end)) THEN
    irl_end = opt_rl_end
  ELSE
    irl_end = min_rledge
  ENDIF

  i_startidx = MERGE(1, &
    &                MAX(1,p_patch%edges%start_index(irl_start)), &
    &                i_blk /= i_startblk)
  i_endidx   = MERGE(nproma,           p_patch%edges%end_index(irl_end), &
    &                i_blk /= i_endblk)

END SUBROUTINE get_indices_e

!-------------------------------------------------------------------------
!
!! Computes the start and end indices of do loops for cell-based variables.
!!
SUBROUTINE get_indices_v(p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                         i_endidx, irl_start, opt_rl_end)


  TYPE(t_patch), INTENT(IN) :: p_patch
  INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
  INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
  INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
  INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts

  INTEGER, OPTIONAL, INTENT(IN) :: opt_rl_end ! refin_ctrl level where do loop ends

  INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (jv loop)

  ! Local variables

  INTEGER :: irl_end

  IF (PRESENT(opt_rl_end)) THEN
    irl_end = opt_rl_end
  ELSE
    irl_end = min_rlvert
  ENDIF

  IF (i_blk == i_startblk) THEN
    i_startidx = p_patch%verts%start_index(irl_start)
    i_endidx   = nproma
    IF (i_blk == i_endblk) i_endidx = p_patch%verts%end_index(irl_end)
  ELSE IF (i_blk == i_endblk) THEN
    i_startidx = 1
    i_endidx   = p_patch%verts%end_index(irl_end)
  ELSE
    i_startidx = 1
    i_endidx = nproma
  ENDIF

END SUBROUTINE get_indices_v


END MODULE mo_loopindices

