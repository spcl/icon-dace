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

! This module contains subroutines needed to determine the start and end
! indices of do loops for a given patch and block index.

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
USE mo_parallel_config, ONLY: nproma
USE mo_lib_loopindices, ONLY: get_indices_c_lib, get_indices_e_lib, get_indices_v_lib

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

  INTEGER :: irl_end, i_startidx_in, i_endidx_in

  i_startidx_in = p_patch%cells%start_index(irl_start)

  IF (PRESENT(opt_rl_end)) THEN
    irl_end = opt_rl_end
  ELSE
    irl_end = min_rlcell
  ENDIF

  i_endidx_in = p_patch%cells%end_index(irl_end)

  CALL get_indices_c_lib(i_startidx_in, i_endidx_in, nproma, i_blk, i_startblk, i_endblk, &
                         i_startidx, i_endidx)

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

  INTEGER :: irl_end, i_startidx_in, i_endidx_in

  i_startidx_in = p_patch%edges%start_index(irl_start)

  IF (PRESENT(opt_rl_end)) THEN
    irl_end = opt_rl_end
  ELSE
    irl_end = min_rledge
  ENDIF

  i_endidx_in = p_patch%edges%end_index(irl_end)

  CALL get_indices_e_lib(i_startidx_in, i_endidx_in, nproma, i_blk, i_startblk, i_endblk, &
                         i_startidx, i_endidx)

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

  INTEGER :: irl_end, i_startidx_in, i_endidx_in

  i_startidx_in = p_patch%verts%start_index(irl_start)

  IF (PRESENT(opt_rl_end)) THEN
    irl_end = opt_rl_end
  ELSE
    irl_end = min_rlvert
  ENDIF

  i_endidx_in = p_patch%verts%end_index(irl_end)

  CALL get_indices_v_lib(i_startidx_in, i_endidx_in, nproma, i_blk, i_startblk, i_endblk, &
                         i_startidx, i_endidx)

END SUBROUTINE get_indices_v

END MODULE mo_loopindices

