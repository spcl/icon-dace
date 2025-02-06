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
USE mo_parallel_config,     ONLY: nproma
use mo_lib_intp_rbf,        ONLY: rbf_vec_interpol_cell_lib, rbf_interpol_c2grad_lib, &
                                  rbf_vec_interpol_vertex_lib, rbf_vec_interpol_edge_lib
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

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx_in      ! start index
INTEGER :: i_endidx_in        ! end index
INTEGER :: rl_start, rl_end

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

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

CALL rbf_vec_interpol_cell_lib( p_vn_in, ptr_int%rbf_vec_idx_c, ptr_int%rbf_vec_blk_c, &
  &                             ptr_int%rbf_vec_coeff_c, p_u_out, p_v_out, &
  &                             i_startblk, i_endblk, i_startidx_in, i_endidx_in, & 
  &                             slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )

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

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx_in      ! start index
INTEGER :: i_endidx_in        ! end index
INTEGER :: rl_start, rl_end

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

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

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

!$OMP END PARALLEL

CALL rbf_interpol_c2grad_lib( p_cell_in, ptr_int%rbf_c2grad_idx, ptr_int%rbf_c2grad_blk, &
  &                           ptr_int%rbf_c2grad_coeff, grad_x, grad_y, & 
  &                           i_startblk, i_endblk, i_startidx_in, i_endidx_in, & 
  &                           slev, elev, nproma, lacc=i_am_accel_node )

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

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx_in      ! start index
INTEGER :: i_endidx_in        ! end index
INTEGER :: rl_start, rl_end

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

! values for the blocking
i_startblk = ptr_patch%verts%start_block(rl_start)
i_endblk   = ptr_patch%verts%end_block(rl_end)

i_startidx_in = ptr_patch%verts%start_index(rl_start)
i_endidx_in   = ptr_patch%verts%end_index(rl_end)

CALL rbf_vec_interpol_vertex_lib( p_e_in, ptr_int%rbf_vec_idx_v, ptr_int%rbf_vec_blk_v, &
  &                               ptr_int%rbf_vec_coeff_v, p_u_out, p_v_out, &
  &                               i_startblk, i_endblk, i_startidx_in, i_endidx_in, & 
  &                               slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )

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

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx_in      ! start index
INTEGER :: i_endidx_in        ! end index
INTEGER :: rl_start, rl_end

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

! values for the blocking
i_startblk = ptr_patch%verts%start_block(rl_start)
i_endblk   = ptr_patch%verts%end_block(rl_end)

i_startidx_in = ptr_patch%verts%start_index(rl_start)
i_endidx_in   = ptr_patch%verts%end_index(rl_end)

CALL rbf_vec_interpol_vertex_lib( p_e_in, ptr_int%rbf_vec_idx_v, ptr_int%rbf_vec_blk_v, &
  &                               ptr_int%rbf_vec_coeff_v, p_u_out, p_v_out, &
  &                               i_startblk, i_endblk, i_startidx_in, i_endidx_in, & 
  &                               slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )

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

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx_in      ! start index
INTEGER :: i_endidx_in        ! end index
INTEGER :: rl_start, rl_end

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

! values for the blocking
i_startblk = ptr_patch%edges%start_block(rl_start)
i_endblk   = ptr_patch%edges%end_block(rl_end)

i_startidx_in = ptr_patch%edges%start_index(rl_start)
i_endidx_in   = ptr_patch%edges%end_index(rl_end)

CALL rbf_vec_interpol_edge_lib( p_vn_in, ptr_int%rbf_vec_idx_e, ptr_int%rbf_vec_blk_e, & 
  &                             ptr_int%rbf_vec_coeff_e, p_vt_out, &
  &                             i_startblk, i_endblk, i_startidx_in, i_endidx_in, & 
  &                             slev, elev, nproma, lacc=i_am_accel_node, acc_async=opt_acc_async )
  
END SUBROUTINE rbf_vec_interpol_edge


END MODULE mo_intp_rbf
