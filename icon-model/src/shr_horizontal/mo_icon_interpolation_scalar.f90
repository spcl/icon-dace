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

! Contains the implementation of interpolation and reconstruction
! routines used by the shallow water model, including the RBF
! reconstruction routines.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_icon_interpolation_scalar
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: dp, sp, wp, vp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: timers_level
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_intp
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device, assert_lacc_equals_i_am_accel_node
  USE mo_mpi,                 ONLY: i_am_accel_node
  USE mo_lib_interpolation_scalar, ONLY: verts2edges_scalar_lib, cells2edges_scalar_lib, &
                                         edges2verts_scalar_lib, edges2cells_scalar_lib, &
                                         cells2verts_scalar_lib, cells2verts_scalar_ri_lib, &
                                         verts2cells_scalar_lib, cell_avg_lib

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: verts2edges_scalar
  PUBLIC :: cells2edges_scalar
  PUBLIC :: edges2verts_scalar
  PUBLIC :: edges2cells_scalar
  PUBLIC :: cells2verts_scalar
  PUBLIC :: cells2verts_scalar_ri
  PUBLIC :: verts2cells_scalar
  PUBLIC :: cell_avg

  INTERFACE edges2cells_scalar
    MODULE PROCEDURE edges2cells_scalar_dp, edges2cells_scalar_sp
  END INTERFACE edges2cells_scalar

  INTERFACE cells2verts_scalar
    MODULE PROCEDURE cells2verts_scalar_dp, cells2verts_scalar_sp
    MODULE PROCEDURE cells2verts_scalar_sp2dp
  END INTERFACE cells2verts_scalar

CONTAINS


!-----------------------------------------------------------------------
!
!  ! averaging and interpolation routines and
!  ! routines needed to compute the coefficients therein
!
!-----------------------------------------------------------------------
!
!! Performs  average of scalar fields from vertices to velocity points.
!!
!! The coefficients are given by c_int.
!!
SUBROUTINE verts2edges_scalar( p_vertex_in, ptr_patch, c_int, p_edge_out, &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! vertex based scalar input field
REAL(wp), INTENT(in) ::  p_vertex_in(:,:,:)  ! dim: (nproma,nlev,nblks_v)
! interpolation field
REAL(wp), INTENT(in) ::  c_int(:,:,:)        ! dim: (nproma,2,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev   ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev   ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! edge based scalar output field
REAL(wp), INTENT(inout) :: p_edge_out(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: i_startblk     ! start block
INTEGER :: i_endblk       ! end block
INTEGER :: rl_start, rl_end, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_vertex_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

i_startblk = ptr_patch%edges%start_block(rl_start)
i_endblk   = ptr_patch%edges%end_block(rl_end)

i_startidx_in = ptr_patch%edges%start_index(rl_start)
i_endidx_in   = ptr_patch%edges%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL verts2edges_scalar_lib( p_vertex_in, ptr_patch%edges%vertex_idx,  ptr_patch%edges%vertex_blk, &
  &                      c_int, p_edge_out, i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                      slev, elev, nproma, lacc=i_am_accel_node )

IF (timers_level > 10) CALL timer_stop(timer_intp)


END SUBROUTINE verts2edges_scalar

!------------------------------------------------------------------------
!
!
!!  Computes  average of scalar fields from centers of triangular faces to
!!  velocity points.
!!
SUBROUTINE cells2edges_scalar( p_cell_in, ptr_patch, c_int, p_edge_out,                    &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend, opt_fill_latbc, &
  &                            lacc)
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for linear interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,2,nblks_e)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

LOGICAL, INTENT(in), OPTIONAL :: opt_fill_latbc  ! if true, fill lateral nest boundaries
LOGICAL, INTENT(in), OPTIONAL :: lacc  ! if true, use openACC

! edge based scalar output field
REAL(wp), INTENT(inout) :: p_edge_out(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level

INTEGER, DIMENSION(2) :: i_startblk_in                ! start block
INTEGER, DIMENSION(2) :: i_endblk_in                  ! end block
INTEGER, DIMENSION(2) :: i_startidx_in                ! start index
INTEGER, DIMENSION(2) :: i_endidx_in                  ! end index

INTEGER :: rl_start, rl_end
LOGICAL :: lfill_latbc, lzacc

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_cell_in,2) 
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  ! The calculation cannot be done for boundary edges
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_interpolation:cells2edges_scalar',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF
IF ( PRESENT(opt_fill_latbc) ) THEN
  lfill_latbc = opt_fill_latbc
ELSE
  lfill_latbc = .FALSE.
END IF

CALL set_acc_host_or_device(lzacc, lacc)
CALL assert_lacc_equals_i_am_accel_node('mo_interpolation:cells2edges_scalar', lzacc, i_am_accel_node)

i_startblk_in(1) = ptr_patch%edges%start_block(1)
i_endblk_in(1)   = ptr_patch%edges%end_block(1)

i_startblk_in(2) = ptr_patch%edges%start_block(rl_start)
i_endblk_in(2)   = ptr_patch%edges%end_block(rl_end)

i_startidx_in(1) = ptr_patch%edges%start_index(1)
i_endidx_in(1)   = ptr_patch%edges%end_index(1)

i_startidx_in(2) = ptr_patch%edges%start_index(rl_start)
i_endidx_in(2)   = ptr_patch%edges%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL cells2edges_scalar_lib( p_cell_in, ptr_patch%edges%cell_idx, ptr_patch%edges%cell_blk, c_int, p_edge_out, & 
  &                          i_startblk_in, i_endblk_in, i_startidx_in, i_endidx_in, & 
  &                          slev, elev, nproma, ptr_patch%id, l_limited_area, lfill_latbc, lzacc)

IF (timers_level > 10) CALL timer_stop(timer_intp)

END SUBROUTINE cells2edges_scalar


!------------------------------------------------------------------------
!
!!  Computes average of scalar fields from velocity points to
!!  centers of dual faces.
!!
SUBROUTINE edges2verts_scalar( p_edge_in, ptr_patch, v_int, p_vert_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! edge based scalar input field
REAL(wp), INTENT(in) ::  p_edge_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! coefficients for (area weighted) interpolation
REAL(wp), INTENT(in) ::  v_int(:,:,:)      ! dim: (nproma,cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_vert_out(:,:,:)  ! dim: (nproma,nlev,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-------------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_edge_in,2)
END IF


IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

! values for the blocking
i_startblk = ptr_patch%verts%start_block(rl_start)
i_endblk   = ptr_patch%verts%end_block(rl_end)

i_startidx_in = ptr_patch%verts%start_index(rl_start)
i_endidx_in   = ptr_patch%verts%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL edges2verts_scalar_lib( p_edge_in, ptr_patch%verts%edge_idx, ptr_patch%verts%edge_blk, v_int, p_vert_out, & 
   &                         i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
   &                         slev, elev, nproma, lacc=i_am_accel_node )

IF (timers_level > 10) CALL timer_stop(timer_intp)


END SUBROUTINE edges2verts_scalar

!------------------------------------------------------------------------
!
!!  Computes interpolation of scalar fields from velocity points to
!!  cell centers via given interpolation weights
!!
SUBROUTINE edges2cells_scalar_dp(p_edge_in, ptr_patch, c_int, p_cell_out,  &
  &                              opt_slev, opt_elev, opt_rlstart, opt_rlend)
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! edge based scalar input field
REAL(dp), INTENT(in) ::  p_edge_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! coefficients for (area weighted) interpolation
REAL(wp), INTENT(in) ::  c_int(:,:,:)      ! dim: (nproma,cell_type,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! cell based scalar output field
REAL(dp), INTENT(inout) :: p_cell_out(:,:,:)  ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-------------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_edge_in,2)
END IF


IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL edges2cells_scalar_lib(p_edge_in, ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, c_int, p_cell_out, &
  &                         i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                         slev, elev, nproma, lacc=i_am_accel_node )

IF (timers_level > 10) CALL timer_stop(timer_intp)

END SUBROUTINE edges2cells_scalar_dp
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
!!  Computes interpolation from edges to cells
!!
!!  Computes interpolation of scalar fields from velocity points to
!!  cell centers via given interpolation weights
!!
SUBROUTINE edges2cells_scalar_sp( p_edge_in, ptr_patch, c_int, p_cell_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! edge based scalar input field
REAL(sp), INTENT(in) ::  p_edge_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! coefficients for (area weighted) interpolation
REAL(wp), INTENT(in) ::  c_int(:,:,:)      ! dim: (nproma,cell_type,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! cell based scalar output field
REAL(sp), INTENT(inout) :: p_cell_out(:,:,:)  ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-------------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_edge_in,2)
END IF


IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL edges2cells_scalar_lib(p_edge_in, ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, c_int, p_cell_out, &
  &                         i_startblk, i_endblk, i_startidx_in, i_endidx_in,  &
  &                         slev, elev, nproma, lacc=i_am_accel_node )

IF (timers_level > 10) CALL timer_stop(timer_intp)

END SUBROUTINE edges2cells_scalar_sp
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!!  Computes  average of scalar fields from centers of cells to vertices.
!!
SUBROUTINE cells2verts_scalar_dp( p_cell_in, ptr_patch, c_int, p_vert_out,    &
  &                               opt_slev, opt_elev, opt_rlstart, opt_rlend, &
  &                               opt_acc_async )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(dp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,9-cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

LOGICAL, INTENT(in), OPTIONAL :: opt_acc_async

! vertex based scalar output field
REAL(dp), INTENT(inout) :: p_vert_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_cell_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

! values for the blocking
i_startblk = ptr_patch%verts%start_block(rl_start)
i_endblk   = ptr_patch%verts%end_block(rl_end)

i_startidx_in = ptr_patch%verts%start_index(rl_start)
i_endidx_in   = ptr_patch%verts%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL cells2verts_scalar_lib( p_cell_in, ptr_patch%verts%cell_idx, ptr_patch%verts%cell_blk, c_int, p_vert_out,&
  &                          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                          slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )

IF (timers_level > 10) CALL timer_stop(timer_intp)


END SUBROUTINE cells2verts_scalar_dp
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!! Computes  average of scalar fields from centers of cells to vertices.
!!
SUBROUTINE cells2verts_scalar_sp( p_cell_in, ptr_patch, c_int, p_vert_out,    &
  &                               opt_slev, opt_elev, opt_rlstart, opt_rlend, &
  &                               opt_acc_async )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(sp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,9-cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

LOGICAL, INTENT(in), OPTIONAL :: opt_acc_async

! vertex based scalar output field
REAL(sp), INTENT(inout) :: p_vert_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_cell_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

! values for the blocking
i_startblk = ptr_patch%verts%start_block(rl_start)
i_endblk   = ptr_patch%verts%end_block(rl_end)

i_startidx_in = ptr_patch%verts%start_index(rl_start)
i_endidx_in   = ptr_patch%verts%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL cells2verts_scalar_lib( p_cell_in, ptr_patch%verts%cell_idx, ptr_patch%verts%cell_blk, c_int, p_vert_out, &
  &                          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                          slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )

IF (timers_level > 10) CALL timer_stop(timer_intp)

END SUBROUTINE cells2verts_scalar_sp
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!! Computes  average of scalar fields from centers of cells to vertices.
!!
  SUBROUTINE cells2verts_scalar_sp2dp(p_cell_in, ptr_patch, c_int, p_vert_out, &
    &                                 opt_slev, opt_elev, opt_rlstart, &
    &                                 opt_rlend, opt_acc_async)
    !
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

    ! cell based scalar input field
    REAL(sp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

    ! coefficients for interpolation
    REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,9-cell_type,nblks_v)

    INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

    INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

    ! start and end values of refin_ctrl flag
    INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

    LOGICAL, INTENT(in), OPTIONAL :: opt_acc_async

    ! vertex based scalar output field
    REAL(dp), INTENT(inout) :: p_vert_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

    !-----------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = UBOUND(p_cell_in,2)
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      rl_start = opt_rlstart
    ELSE
      rl_start = 2
    END IF
    IF ( PRESENT(opt_rlend) ) THEN
      rl_end = opt_rlend
    ELSE
      rl_end = min_rlvert
    END IF

    ! values for the blocking
    i_startblk = ptr_patch%verts%start_block(rl_start)
    i_endblk   = ptr_patch%verts%end_block(rl_end)

    i_startidx_in = ptr_patch%verts%start_index(rl_start)
    i_endidx_in   = ptr_patch%verts%end_index(rl_end)

    IF (timers_level > 10) CALL timer_start(timer_intp)

    CALL cells2verts_scalar_lib( p_cell_in, ptr_patch%verts%cell_idx, ptr_patch%verts%cell_blk, &
      &                          c_int, p_vert_out, &
      &                          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
      &                          slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )

    IF (timers_level > 10) CALL timer_stop(timer_intp)

  END SUBROUTINE cells2verts_scalar_sp2dp
!------------------------------------------------------------------------

!>
!!  Same as above, but provides output optionally in single precision and
!!  assumes reversed index order of the output field in loop exchange mode
!!
!!
SUBROUTINE cells2verts_scalar_ri( p_cell_in, ptr_patch, c_int, p_vert_out,    &
  &                               opt_slev, opt_elev, opt_rlstart, opt_rlend, &
  &                               opt_acc_async )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,9-cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! vertex based scalar output field
REAL(vp), INTENT(inout) :: p_vert_out(:,:,:) ! dim: (nlev,nproma,nblks_v) or (nproma,nlev,nblks_v)

LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async   !< optional async OpenACC

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_cell_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

! values for the blocking
i_startblk = ptr_patch%verts%start_block(rl_start)
i_endblk   = ptr_patch%verts%end_block(rl_end)

i_startidx_in = ptr_patch%verts%start_index(rl_start)
i_endidx_in   = ptr_patch%verts%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL cells2verts_scalar_ri_lib( p_cell_in, ptr_patch%verts%cell_idx, ptr_patch%verts%cell_blk, &
  &                               c_int, p_vert_out, &
  &                               i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                               slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )

IF (timers_level > 10) CALL timer_stop(timer_intp)

END SUBROUTINE cells2verts_scalar_ri
!------------------------------------------------------------------------

!
!! Computes  average of scalar fields from vertices to centers of cells.
!!
SUBROUTINE verts2cells_scalar( p_vert_in, ptr_patch, c_int, p_cell_out,  &
  &                            opt_slev, opt_elev )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_vert_in(:,:,:)   ! dim: (nproma,nlev,nblks_v)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,cell_type,nblks_c)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_cell_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: nblks_c, npromz_c

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_vert_in,2)
END IF

! values for the blocking
nblks_c  = ptr_patch%nblks_c
npromz_c = ptr_patch%npromz_c

IF (timers_level > 10) CALL timer_start(timer_intp)

CALL verts2cells_scalar_lib( p_vert_in, ptr_patch%cells%vertex_idx, ptr_patch%cells%vertex_blk, &
  &                          c_int, p_cell_out, nblks_c, npromz_c, &
  &                          slev, elev, nproma, lacc=i_am_accel_node )

IF (timers_level > 10) CALL timer_stop(timer_intp)

END SUBROUTINE verts2cells_scalar

!-------------------------------------------------------------------------
!
!
!! Computes the average of a cell-based variable
!! over its original location and the neighboring triangles.
!! Version with variable weighting coefficients, computed such that
!! linear horizontal gradients are not aliased into a checkerboard noise
!! input:  lives on centers of triangles
!! output: lives on centers of triangles
!!
SUBROUTINE cell_avg( psi_c, ptr_patch, avg_coeff, avg_psi_c,     &
  &                  opt_slev, opt_elev, opt_rlstart, opt_rlend, &
  &                  lacc )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  averaging coefficients
!
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)

!
!  cell based variable before averaging
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

LOGICAL, INTENT(IN), OPTIONAL :: & 
  &  lacc    ! if true, use OpenACC

!
!   cell based variable after averaging
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  avg_psi_c(:,:,:)  ! dim: (nproma,nlev,nblks_c)


INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
LOGICAL :: lzacc ! non-optional version of lacc

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(psi_c,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_interpolation:cell_avg',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

CALL set_acc_host_or_device(lzacc, lacc)

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

!
! loop through all patch cells (and blocks)
!
IF (timers_level > 10) CALL timer_start(timer_intp)

CALL cell_avg_lib( psi_c, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, avg_coeff, avg_psi_c, &
  &                i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                slev, elev, nproma, lacc=lzacc )

  IF (timers_level > 10) CALL timer_stop(timer_intp)


END SUBROUTINE cell_avg

!-------------------------------------------------------------------------

END MODULE mo_icon_interpolation_scalar
