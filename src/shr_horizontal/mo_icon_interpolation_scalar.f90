! Contains the implementation of interpolation and reconstruction
! routines used by the shallow water model, including the RBF
! reconstruction routines.
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

!----------------------------
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
INTEGER :: i_startidx     ! start index
INTEGER :: i_endidx       ! end index
INTEGER :: rl_start, rl_end, i_nchdom
INTEGER :: je, jk, jb

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%edges%vertex_idx
iblk => ptr_patch%edges%vertex_blk

i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

! loop over edges and blocks
!$ACC DATA COPYIN(p_vertex_in, c_int) PRESENT(p_edge_out, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jk = slev, elev
    DO je = i_startidx, i_endidx

      p_edge_out(je,jk,jb) =  &
        &   c_int(je,1,jb)*p_vertex_in(iidx(je,jb,1),jk,iblk(je,jb,1))  &
        & + c_int(je,2,jb)*p_vertex_in(iidx(je,jb,2),jk,iblk(je,jb,2))

    END DO

  END DO
  !$ACC END PARALLEL

END DO
!$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

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
INTEGER :: je, jk, jb
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: rl_start, rl_end, i_nchdom
LOGICAL :: lfill_latbc, lzacc

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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
CALL assert_lacc_equals_i_am_accel_node('mo_interpolation:cells2edges_scalar', lacc)

iidx => ptr_patch%edges%cell_idx
iblk => ptr_patch%edges%cell_blk

i_nchdom   = MAX(1,ptr_patch%n_childdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

!$ACC DATA NO_CREATE(iidx, iblk, c_int, p_cell_in, p_edge_out)

!$OMP PARALLEL PRIVATE(i_startblk, i_endblk)

IF ( (l_limited_area .OR. ptr_patch%id > 1) .AND. lfill_latbc) THEN ! Fill outermost nest boundary

  i_startblk = ptr_patch%edges%start_blk(1,1)
  i_endblk   = ptr_patch%edges%end_blk(1,1)

! DA: OpenACC needs to collapse loops
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk
    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1, 1)

    DO je = i_startidx, i_endidx
      IF (iidx(je,jb,1) >= 1 .AND. iblk(je,jb,1) >= 1) THEN
        DO jk = slev, elev
          p_edge_out(je,jk,jb) =  p_cell_in(iidx(je,jb,1),jk,iblk(je,jb,1))
        END DO
      ELSE IF (iidx(je,jb,2) >= 1 .AND. iblk(je,jb,2) >= 1) THEN
        DO jk = slev, elev
          p_edge_out(je,jk,jb) =  p_cell_in(iidx(je,jb,2),jk,iblk(je,jb,2))
        END DO
      ELSE
        CALL finish ('mo_interpolation:cells2edges_scalar',  &
          &          'error in lateral boundary filling')
      ENDIF
    END DO
  END DO
!$OMP END DO


ENDIF

! Process the remaining grid points for which a real interpolation is possible
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO je = i_startidx, i_endidx

        p_edge_out(je,jk,jb) =  &
          &    c_int(je,1,jb) * p_cell_in(iidx(je,jb,1),jk,iblk(je,jb,1))  &
          &  + c_int(je,2,jb) * p_cell_in(iidx(je,jb,2),jk,iblk(je,jb,2))

      END DO
    END DO
    !$ACC END PARALLEL

  END DO
!$OMP END DO NOWAIT


!$OMP END PARALLEL

!$ACC END DATA

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
INTEGER :: jv, jk, jb
INTEGER :: nblks_v, npromz_v, nlen
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%verts%edge_idx
iblk => ptr_patch%verts%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

  !$ACC DATA PRESENT(v_int, p_edge_in, p_vert_out, iidx, iblk) IF(i_am_accel_node)

!loop over blocks and verts

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx

        p_vert_out(jv,jk,jb) =  &
          v_int(jv,1,jb) * p_edge_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
          v_int(jv,2,jb) * p_edge_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
          v_int(jv,3,jb) * p_edge_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
          v_int(jv,4,jb) * p_edge_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
          v_int(jv,5,jb) * p_edge_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
          v_int(jv,6,jb) * p_edge_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO  !loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA

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
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

!$ACC DATA PRESENT(c_int, p_edge_in, p_cell_out, iidx, iblk) IF(i_am_accel_node)

!loop over blocks and cells

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        p_cell_out(jc,jk,jb) =  &
          c_int(jc,1,jb) * p_edge_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          c_int(jc,2,jb) * p_edge_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          c_int(jc,3,jb) * p_edge_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO  !loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA

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
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

!$ACC DATA PRESENT(c_int, p_edge_in, p_cell_out, iidx, iblk) IF(i_am_accel_node)

!loop over blocks and cells

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        p_cell_out(jc,jk,jb) = REAL( &
          c_int(jc,1,jb) * REAL(p_edge_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)), wp)+&
          c_int(jc,2,jb) * REAL(p_edge_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)), wp)+&
          c_int(jc,3,jb) * REAL(p_edge_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)), wp),&
          sp)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO  !loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA

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
INTEGER :: jv, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%verts%cell_idx
iblk => ptr_patch%verts%cell_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

!$ACC DATA NO_CREATE(iidx, iblk, p_cell_in, c_int, p_vert_out)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx

         p_vert_out(jv,jk,jb) =                       &
           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
           c_int(jv,4,jb) * p_cell_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
           c_int(jv,5,jb) * p_cell_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
           c_int(jv,6,jb) * p_cell_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF ( PRESENT(opt_acc_async) ) THEN
    IF ( .NOT. opt_acc_async ) THEN
      !$ACC WAIT
    END IF
  ELSE
    !$ACC WAIT
  END IF

  !$ACC END DATA

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
INTEGER :: jv, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%verts%cell_idx
iblk => ptr_patch%verts%cell_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

!$ACC DATA NO_CREATE(p_cell_in, c_int, p_vert_out, iidx, iblk)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx

         p_vert_out(jv,jk,jb) =                       &
           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
           c_int(jv,4,jb) * p_cell_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
           c_int(jv,5,jb) * p_cell_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
           c_int(jv,6,jb) * p_cell_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF ( PRESENT(opt_acc_async) ) THEN
    IF ( .NOT. opt_acc_async ) THEN
      !$ACC WAIT
    END IF
  ELSE
    !$ACC WAIT
  END IF

  !$ACC END DATA

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
    INTEGER :: jv, jk, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

    iidx => ptr_patch%verts%cell_idx
    iblk => ptr_patch%verts%cell_blk

    ! values for the blocking
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%verts%start_blk(rl_start,1)
    i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


    IF (timers_level > 10) CALL timer_start(timer_intp)

    !$ACC DATA NO_CREATE(p_cell_in, c_int, p_vert_out, iidx, iblk)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=6
      DO jk = slev, elev
        DO jv = i_startidx, i_endidx

          p_vert_out(jv,jk,jb) =                                                    &
            c_int(jv,1,jb) * REAL(p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)), dp) + &
            c_int(jv,2,jb) * REAL(p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)), dp) + &
            c_int(jv,3,jb) * REAL(p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)), dp) + &
            c_int(jv,4,jb) * REAL(p_cell_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)), dp) + &
            c_int(jv,5,jb) * REAL(p_cell_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)), dp) + &
            c_int(jv,6,jb) * REAL(p_cell_in(iidx(jv,jb,6),jk,iblk(jv,jb,6)), dp)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF ( PRESENT(opt_acc_async) ) THEN
      IF ( .NOT. opt_acc_async ) THEN
        !$ACC WAIT
      END IF
    ELSE
      !$ACC WAIT
    END IF

    !$ACC END DATA

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
INTEGER :: jv, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%verts%cell_idx
iblk => ptr_patch%verts%cell_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
         p_vert_out(jv,jk,jb) =                                         &
           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
           c_int(jv,4,jb) * p_cell_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
           c_int(jv,5,jb) * p_cell_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
           c_int(jv,6,jb) * p_cell_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF ( PRESENT(opt_acc_async) ) THEN
    IF ( .NOT. opt_acc_async ) THEN
      !$ACC WAIT
    END IF
  ELSE
    !$ACC WAIT
  END IF

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
INTEGER :: jk, jb, jc, nlen, nblks_c, npromz_c

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%cells%vertex_idx
iblk => ptr_patch%cells%vertex_blk

! values for the blocking
nblks_c  = ptr_patch%nblks_c
npromz_c = ptr_patch%npromz_c


IF (timers_level > 10) CALL timer_start(timer_intp)

!$ACC DATA PRESENT(p_vert_in, c_int, p_cell_out, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_c

    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
    DO jk = slev, elev
      DO jc = 1, nlen

        p_cell_out(jc,jk,jb) =                                        &
          c_int(jc,1,jb)* p_vert_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          c_int(jc,2,jb)* p_vert_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          c_int(jc,3,jb)* p_vert_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      ENDDO
    ENDDO
    !$ACC END PARALLEL

  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA

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
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
LOGICAL :: lzacc ! non-optional version of lacc

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

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

iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!


IF (timers_level > 10) CALL timer_start(timer_intp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        !  calculate the weighted average
        !
        avg_psi_c(jc,jk,jb) =  &
          &    psi_c(jc,jk,jb)                       * avg_coeff(jc,1,jb) &
          &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop

    END DO !vertical levels loop
    !$ACC END PARALLEL

  END DO !block loop

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF (timers_level > 10) CALL timer_stop(timer_intp)


END SUBROUTINE cell_avg

!-------------------------------------------------------------------------


END MODULE mo_icon_interpolation_scalar
