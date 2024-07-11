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
#include "omp_definitions.inc"
!----------------------------

MODULE mo_intp_rbf
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp, sp
USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int
USE mo_model_domain,        ONLY: t_patch
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_fortran_tools,       ONLY: init
USE mo_mpi,                 ONLY: i_am_accel_node

IMPLICIT NONE


PRIVATE

PUBLIC :: rbf_vec_interpol_cell, rbf_interpol_c2grad,     &
        & rbf_vec_interpol_vertex, rbf_vec_interpol_edge

INTERFACE rbf_vec_interpol_vertex
  MODULE PROCEDURE rbf_vec_interpol_vertex_wp
  MODULE PROCEDURE rbf_vec_interpol_vertex_vp
END INTERFACE


CONTAINS

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!! Performs vector RBF reconstruction at cell center.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each cell center.
!!
SUBROUTINE rbf_vec_interpol_cell( p_vn_in, ptr_patch, ptr_int, p_u_out,  &
  &                               p_v_out, opt_slev, opt_elev, opt_rlstart, &
  &                               opt_rlend, opt_acc_async )

! !INPUT PARAMETERS
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input normal components of (velocity) vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed x-component (u) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_u_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

! reconstructed y-component (v) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

! if set, run in an asynchrounous device stream
LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

! !LOCAL VARIABLES
INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jc, jk, jb                ! integer over cells, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom, jk0, jkk


INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  elev = UBOUND(p_vn_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell_int-1
END IF

iidx => ptr_int%rbf_vec_idx_c
iblk => ptr_int%rbf_vec_blk_c

ptr_coeff => ptr_int%rbf_vec_coeff_c

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), ICON_OMP_RUNTIME_SCHEDULE

DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
  !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = slev, elev
#else
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx
#endif
      p_u_out(jc,jk,jb) =  &
        ptr_coeff(1,1,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
        ptr_coeff(2,1,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,1,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,1,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,1,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,1,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,1,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,1,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,1,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))
      p_v_out(jc,jk,jb) =  &
        ptr_coeff(1,2,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
        ptr_coeff(2,2,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,2,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,2,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,2,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,2,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,2,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,2,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,2,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))
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

END SUBROUTINE rbf_vec_interpol_cell
!====================================================================================


!====================================================================================
SUBROUTINE rbf_interpol_c2grad( p_cell_in, ptr_patch, ptr_int, grad_x,  &
  &                             grad_y, opt_slev, opt_elev, opt_rlstart, opt_rlend )

! !INPUT PARAMETERS
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input cell-based variable for which gradient at cell center is computed
REAL(wp), INTENT(IN) ::  &
  &  p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed zonal (x) component of gradient vector
REAL(wp),INTENT(INOUT) ::  &
  &  grad_x(:,:,:) ! dim: (nproma,nlev,nblks_c)

! reconstructed zonal (x) component of gradient vector
REAL(wp),INTENT(INOUT) ::  &
  &  grad_y(:,:,:) ! dim: (nproma,nlev,nblks_c)

! !LOCAL VARIABLES
INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jc, jk, jb                ! integer over cells, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom


INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  rl_end = min_rlcell_int
END IF


iidx => ptr_int%rbf_c2grad_idx
iblk => ptr_int%rbf_c2grad_blk

ptr_coeff => ptr_int%rbf_c2grad_coeff

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL

IF (ptr_patch%id > 1) THEN
#ifdef _OPENACC
  !$ACC KERNELS ASYNC(1) IF(i_am_accel_node)
  grad_x(:,:,1:i_startblk) = 0._wp
  grad_y(:,:,1:i_startblk) = 0._wp
  !$ACC END KERNELS
#else
  CALL init(grad_x(:,:,1:i_startblk), lacc=i_am_accel_node)
  CALL init(grad_y(:,:,1:i_startblk), lacc=i_am_accel_node)
!$OMP BARRIER
#endif
ENDIF

!$ACC DATA PRESENT(p_cell_in, grad_x, grad_y, iidx, iblk, ptr_coeff) IF(i_am_accel_node)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), ICON_OMP_RUNTIME_SCHEDULE

DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
  !$ACC LOOP GANG
#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    !$ACC LOOP VECTOR
    DO jk = slev, elev
#else
  DO jk = slev, elev
    !$ACC LOOP VECTOR
    DO jc = i_startidx, i_endidx
#endif

      grad_x(jc,jk,jb) =  &
        ptr_coeff(1,1,jc,jb)*p_cell_in(jc,jk,jb) + &
        ptr_coeff(2,1,jc,jb)*p_cell_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,1,jc,jb)*p_cell_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,1,jc,jb)*p_cell_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,1,jc,jb)*p_cell_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,1,jc,jb)*p_cell_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,1,jc,jb)*p_cell_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,1,jc,jb)*p_cell_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,1,jc,jb)*p_cell_in(iidx(9,jc,jb),jk,iblk(9,jc,jb)) + &
        ptr_coeff(10,1,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))
      grad_y(jc,jk,jb) =  &
        ptr_coeff(1,2,jc,jb)*p_cell_in(jc,jk,jb) + &
        ptr_coeff(2,2,jc,jb)*p_cell_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,2,jc,jb)*p_cell_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,2,jc,jb)*p_cell_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,2,jc,jb)*p_cell_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,2,jc,jb)*p_cell_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,2,jc,jb)*p_cell_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,2,jc,jb)*p_cell_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,2,jc,jb)*p_cell_in(iidx(9,jc,jb),jk,iblk(9,jc,jb)) + &
        ptr_coeff(10,2,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))

    ENDDO
  ENDDO
  !$ACC END PARALLEL

ENDDO
!$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE rbf_interpol_c2grad

!-------------------------------------------------------------------------
!
!! Performs vector RBF reconstruction at triangle vertices.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each vertex.
!!
SUBROUTINE rbf_vec_interpol_vertex_wp( p_e_in, ptr_patch, ptr_int,                 &
                                       p_u_out, p_v_out,                           &
                                       opt_slev, opt_elev, opt_rlstart, opt_rlend, &
                                       opt_acc_async )
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input components of velocity or horizontal vorticity vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_e_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed x-component (u) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_u_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

! reconstructed y-component (v) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

! !LOCAL VARIABLES

INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jv, jk, jb                ! integer over vertices, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  elev = UBOUND(p_e_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert_int-1
END IF


iidx => ptr_int%rbf_vec_idx_v
iblk => ptr_int%rbf_vec_blk_v

ptr_coeff => ptr_int%rbf_vec_coeff_v

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jv), ICON_OMP_RUNTIME_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
  DO jv = i_startidx, i_endidx
    DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
  DO jk = slev, elev
    DO jv = i_startidx, i_endidx
#endif

      p_u_out(jv,jk,jb) =  &
        ptr_coeff(1,1,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,1,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,1,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,1,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,1,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,1,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))
      p_v_out(jv,jk,jb) =  &
        ptr_coeff(1,2,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,2,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,2,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,2,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,2,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,2,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))

      ENDDO
    ENDDO

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

END SUBROUTINE rbf_vec_interpol_vertex_wp

! Variant for mixed precision mode (output fields in single precision)
SUBROUTINE rbf_vec_interpol_vertex_vp( p_e_in, ptr_patch, ptr_int, &
                                       p_u_out, p_v_out,           &
                                       opt_slev, opt_elev, opt_rlstart, opt_rlend,  &
                                       opt_acc_async )
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input components of velocity or horizontal vorticity vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_e_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed x-component (u) of velocity vector
REAL(sp),INTENT(INOUT) ::  &
  &  p_u_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

! reconstructed y-component (v) of velocity vector
REAL(sp),INTENT(INOUT) ::  &
  &  p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

! !LOCAL VARIABLES

INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jv, jk, jb                ! integer over vertices, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  elev = UBOUND(p_e_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert_int-1
END IF

iidx => ptr_int%rbf_vec_idx_v
iblk => ptr_int%rbf_vec_blk_v

ptr_coeff => ptr_int%rbf_vec_coeff_v

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jv), ICON_OMP_RUNTIME_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
  DO jv = i_startidx, i_endidx
    DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
  DO jk = slev, elev
    DO jv = i_startidx, i_endidx
#endif

      p_u_out(jv,jk,jb) =  &
        ptr_coeff(1,1,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,1,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,1,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,1,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,1,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,1,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))
      p_v_out(jv,jk,jb) =  &
        ptr_coeff(1,2,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,2,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,2,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,2,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,2,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,2,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))

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

END SUBROUTINE rbf_vec_interpol_vertex_vp

!-------------------------------------------------------------------------
!
!! Performs vector RBF reconstruction at edge midpoints.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each edge.
!!
SUBROUTINE rbf_vec_interpol_edge( p_vn_in, ptr_patch, ptr_int, p_vt_out,      &
  &                               opt_slev, opt_elev, opt_rlstart, opt_rlend, &
  &                               opt_acc_async )
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input normal components of velocity vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed tangential velocity component
REAL(wp),INTENT(INOUT) ::  &
  &  p_vt_out(:,:,:) ! dim: (nproma,nlev,nblks_e)

! if set, run in an asynchrounous device stream
LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: je, jk, jb                ! integer over edges, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: rl_start, rl_end, i_nchdom

INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff

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
  elev = UBOUND(p_vn_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge_int-2
END IF

iidx => ptr_int%rbf_vec_idx_e
iblk => ptr_int%rbf_vec_blk_e

ptr_coeff => ptr_int%rbf_vec_coeff_e

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO je = i_startidx, i_endidx
#endif

        p_vt_out(je,jk,jb) =  &
          ptr_coeff(1,je,jb)*p_vn_in(iidx(1,je,jb),jk,iblk(1,je,jb)) + &
          ptr_coeff(2,je,jb)*p_vn_in(iidx(2,je,jb),jk,iblk(2,je,jb)) + &
          ptr_coeff(3,je,jb)*p_vn_in(iidx(3,je,jb),jk,iblk(3,je,jb)) + &
          ptr_coeff(4,je,jb)*p_vn_in(iidx(4,je,jb),jk,iblk(4,je,jb))

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
  
END SUBROUTINE rbf_vec_interpol_edge


END MODULE mo_intp_rbf
