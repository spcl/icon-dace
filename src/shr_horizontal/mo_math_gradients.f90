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

! Contains the implementation of the mathematical grad operators
! employed by the shallow water prototype.
!
! @par To Do
! Boundary exchange, nblks in presence of halos and dummy edge

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_math_gradients
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,               ONLY: wp, vp
USE mo_impl_constants,     ONLY: min_rlcell, min_rledge
USE mo_intp_data_strc,     ONLY: t_int_state
USE mo_intp,               ONLY: cells2edges_scalar
USE mo_model_domain,       ONLY: t_patch
USE mo_parallel_config,    ONLY: nproma
USE mo_run_config,         ONLY: timers_level
USE mo_exception,          ONLY: finish
USE mo_timer,              ONLY: timer_start, timer_stop, timer_grad
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
USE mo_fortran_tools,      ONLY: init
USE mo_mpi,                ONLY: i_am_accel_node

IMPLICIT NONE

PRIVATE

PUBLIC :: grad_fd_norm, grad_fd_tang
PUBLIC :: grad_green_gauss_cell
PUBLIC :: grad_fe_cell

INTERFACE grad_green_gauss_cell
  MODULE PROCEDURE grad_green_gauss_cell_adv
  MODULE PROCEDURE grad_green_gauss_cell_dycore
END INTERFACE

INTERFACE grad_fe_cell
  MODULE PROCEDURE grad_fe_cell_3d
  MODULE PROCEDURE grad_fe_cell_2d
END INTERFACE

CONTAINS

!-------------------------------------------------------------------------
!

!-------------------------------------------------------------------------
!
!!  Computes directional  derivative of a cell centered variable.
!!
!!  Computes directional  derivative of a cell centered variable
!!  with respect to direction normal to triangle edge.
!! input: lives on centres of triangles
!! output:  lives on edges (velocity points)
!!
SUBROUTINE grad_fd_norm( psi_c, ptr_patch, grad_norm_psi_e, &
  &                      opt_slev, opt_elev, opt_rlstart, opt_rlend )

!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  cell based variable of which normal derivative is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:)       ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which normal derivative is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  grad_norm_psi_e(:,:,:)  ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!
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
  ! rl_start=1 means edges located along a lateral boundary of a nested
  ! domain. For those, gradient computation is not possible
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:grad_fd_norm',  &
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

iidx => ptr_patch%edges%cell_idx
iblk => ptr_patch%edges%cell_blk

i_nchdom   = MAX(1,ptr_patch%n_childdom)

i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
!
!  loop through all patch edges (and blocks)
!

  IF (timers_level > 10) CALL timer_start(timer_grad)

  !$ACC DATA PRESENT(psi_c, grad_norm_psi_e, ptr_patch%edges%inv_dual_edge_length, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO je = i_startidx, i_endidx
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO je = i_startidx, i_endidx
#endif
      !
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
      !
       grad_norm_psi_e(je,jk,jb) =  &
          &  ( psi_c(iidx(je,jb,2),jk,iblk(je,jb,2)) - &
          &    psi_c(iidx(je,jb,1),jk,iblk(je,jb,1)) )  &
          &  * ptr_patch%edges%inv_dual_edge_length(je,jb)

      ENDDO
    END DO
    !$ACC END PARALLEL

  END DO
  !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA

  IF (timers_level > 10) CALL timer_stop(timer_grad)


END SUBROUTINE grad_fd_norm

!-------------------------------------------------------------------------
!
! RESTRUCT: @Marco: please adjust calls to this routine to your needs.
!!
!! Computes directional derivative of a vertex centered variable with.
!!
!! Computes directional derivative of a vertex centered variable with
!! respect to direction tanget to triangle edge. Notice that the
!! tangential direction is defined by
!!   iorient*(vertex2 - vertex1)
!! input: lives on vertices of triangles
!! output: lives on edges (velocity points)
!!
SUBROUTINE grad_fd_tang( psi_v, ptr_patch, grad_tang_psi_e,  &
  &                      opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
! vertex based variable of which tangential derivative is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_v(:,:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which tangential derivative is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  grad_tang_psi_e(:,:,:)

REAL(wp) :: iorient

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER, DIMENSION(nproma) ::  &
  &  ilv1, ibv1, ilv2, ibv2
!
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
  elev = UBOUND(psi_v,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  ! The possible domain extent depends on the reconstruction algorithm
  ! used to calculate psi_v; the following values are valid for a 6-point
  ! stencil. In the hexagon case (where prognostic variables are located
  ! at vertices), rl_start may be set to 1.
  IF ((opt_rlstart >= 1) .AND. (opt_rlstart <= 2)) THEN
    CALL finish ('mo_math_operators:grad_fd_tang',  &
          &      'opt_rlstart must not be equal to 1 or 2')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 3
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$ACC DATA PRESENT(psi_v, grad_tang_psi_e, ptr_patch) CREATE(ilv1, ibv1, ilv2, ibv2) IF(i_am_accel_node)

!
! TODO: OpenMP
!

!
!  loop through all patch edges (and blocks)
!
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO je = i_startidx, i_endidx
    !
    !  get the line and block indices of the vertices of edge je
    !
    ilv1(je) = ptr_patch%edges%vertex_idx(je,jb,1)
    ibv1(je) = ptr_patch%edges%vertex_blk(je,jb,1)
    ilv2(je) = ptr_patch%edges%vertex_idx(je,jb,2)
    ibv2(je) = ptr_patch%edges%vertex_blk(je,jb,2)
  END DO

  DO jk = slev, elev

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(iorient)
    DO je = i_startidx, i_endidx
      !
      ! compute the tangential derivative
      ! by the finite difference approximation
      iorient = ptr_patch%edges%tangent_orientation(je,jb)
      grad_tang_psi_e(je,jk,jb) = iorient  &
        &  * ( psi_v(ilv2(je),jk,ibv2(je)) - psi_v(ilv1(je),jk,ibv1(je)) )  &
        &    / ptr_patch%edges%primal_edge_length(je,jb)
    END DO

  END DO
  !$ACC END PARALLEL

END DO
!
! TODO: OpenMP
!

!$ACC WAIT(1)
!$ACC END DATA

END SUBROUTINE grad_fd_tang

!-------------------------------------------------------------------------
!
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! The triangular element is spanned by the cell circumcenters of the three
!! direct neighbours. In contrast to the Green-Gauss approach, this
!! approach does not involve the cell center value of the central triangle.
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons, Sec. 7.2, 7.6
!!
!!
SUBROUTINE grad_fe_cell_3d( p_cc, ptr_patch, ptr_int, p_grad, &
  &                      opt_slev, opt_elev, opt_rlstart,  &
  &                      opt_rlend                         )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(vp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)      ! dim:(2,nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

INTEGER, POINTER, CONTIGUOUS :: iidx(:,:,:), iblk(:,:,:)

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
  elev = UBOUND(p_cc,2)
END IF
IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF


iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

!
! 2. reconstruction of cell based geographical gradient
!

  !$ACC DATA PRESENT(p_cc, p_grad, ptr_int%gradc_bmat, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  IF (ptr_patch%id > 1) THEN
  ! Fill nest boundaries with zero to avoid trouble with MPI synchronization

#ifdef _OPENACC
    !$ACC KERNELS PRESENT(p_grad) ASYNC(1) IF(i_am_accel_node)
    p_grad(:,:,:,1:i_startblk) = 0._wp
    !$ACC END KERNELS
#else
    CALL init(p_grad(:,:,:,1:i_startblk), lacc=i_am_accel_node)
!$OMP BARRIER
#endif
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
#endif

        ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
        ! since it is extremely slow on the SX9, when combined with indirect
        ! addressing.

        ! multiply cell-based input values with precomputed grid geometry factor

        ! zonal(u)-component of Green-Gauss gradient
        p_grad(1,jc,jk,jb) = &
          &    ptr_int%gradc_bmat(jc,1,1,jb)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,1,2,jb)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,1,3,jb)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! meridional(v)-component of Green-Gauss gradient
        p_grad(2,jc,jk,jb) =  &
          &    ptr_int%gradc_bmat(jc,2,1,jb)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1))  &
          &  + ptr_int%gradc_bmat(jc,2,2,jb)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2))  &
          &  + ptr_int%gradc_bmat(jc,2,3,jb)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA


END SUBROUTINE grad_fe_cell_3d




!-------------------------------------------------------------------------
!
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! The triangular element is spanned by the cell circumcenters of the three
!! direct neighbours. In contrast to the Green-Gauss approach, this
!! approach does not involve the cell center value of the central triangle.
!!
!! 2D version, i.e. for a single vertical level
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons, Sec. 7.2, 7.6
!!
!!
SUBROUTINE grad_fe_cell_2d( p_cc, ptr_patch, ptr_int, p_grad, &
  &                         opt_rlstart, opt_rlend            )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(wp), INTENT(inout) ::  &
  &  p_grad(:,:,:)      ! dim:(2,nproma,nblks_c)

INTEGER :: jc, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

INTEGER, POINTER, CONTIGUOUS :: iidx(:,:,:), iblk(:,:,:)

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF


iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

!
! 2. reconstruction of cell based geographical gradient
!
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  IF (ptr_patch%id > 1) THEN
  ! Fill nest boundaries with zero to avoid trouble with MPI synchronization
    CALL init(p_grad(:,:,1:i_startblk), lacc=i_am_accel_node)
!$OMP BARRIER
  ENDIF

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


    DO jc = i_startidx, i_endidx

      ! We do not make use of the intrinsic function DOT_PRODUCT on purpose,
      ! since it is extremely slow on the SX9, when combined with indirect
      ! addressing.

      ! multiply cell-based input values with precomputed grid geometry factor

      ! zonal(u)-component of gradient
      p_grad(1,jc,jb) = &
        &    ptr_int%gradc_bmat(jc,1,1,jb)*p_cc(iidx(jc,jb,1),iblk(jc,jb,1))  &
        &  + ptr_int%gradc_bmat(jc,1,2,jb)*p_cc(iidx(jc,jb,2),iblk(jc,jb,2))  &
        &  + ptr_int%gradc_bmat(jc,1,3,jb)*p_cc(iidx(jc,jb,3),iblk(jc,jb,3))

      ! meridional(v)-component of gradient
      p_grad(2,jc,jb) =  &
        &    ptr_int%gradc_bmat(jc,2,1,jb)*p_cc(iidx(jc,jb,1),iblk(jc,jb,1))  &
        &  + ptr_int%gradc_bmat(jc,2,2,jb)*p_cc(iidx(jc,jb,2),iblk(jc,jb,2))  &
        &  + ptr_int%gradc_bmat(jc,2,3,jb)*p_cc(iidx(jc,jb,3),iblk(jc,jb,3))

    END DO ! end loop over cells

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL


END SUBROUTINE grad_fe_cell_2d


!-------------------------------------------------------------------------
!
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The Green-Gauss approach is used. See for example:
!! http://www.cfd-online.com/Wiki/Gradient_computation
!!
SUBROUTINE grad_green_gauss_cell_adv( p_cc, ptr_patch, ptr_int, p_grad, &
  &                                   opt_slev, opt_elev, opt_p_face,   &
  &                                   opt_rlstart, opt_rlend            )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:,:)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(vp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)      ! dim:(2,nproma,nlev,nblks_c)

! optional: calculated face values of cell centered quantity
REAL(wp), INTENT(inout), OPTIONAL ::  &
  &  opt_p_face(:,:,:)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
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
  elev = UBOUND(p_cc,2)
END IF
IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF


iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

i_nchdom = MAX(1,ptr_patch%n_childdom)


! save face values in optional output field
! (the cell-to-edge interpolation is no longer needed otherwise because
!  of using precomputed geometrical factors)
IF ( PRESENT(opt_p_face) ) THEN
  CALL cells2edges_scalar( p_cc, ptr_patch, ptr_int%c_lin_e, opt_p_face,  &
    &                      slev, elev, lacc=i_am_accel_node)
ENDIF


!
! 2. reconstruction of cell based geographical gradient
!
  !$ACC DATA PRESENT(p_cc, p_grad, ptr_int%geofac_grg, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

  IF (ptr_patch%id > 1) THEN
  ! Fill nest boundaries with zero to avoid trouble with MPI synchronization
#ifdef _OPENACC
    !$ACC KERNELS ASYNC(1) IF(i_am_accel_node)
    p_grad(:,:,:,1:i_startblk) = 0._wp
    !$ACC END KERNELS
#else
    CALL init(p_grad(:,:,:,1:i_startblk), lacc=i_am_accel_node)
!$OMP BARRIER
#endif
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
#endif

        ! multiply cell-based input values with precomputed grid geometry factor

        ! zonal(u)-component of Green-Gauss gradient
        p_grad(1,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,1)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,1)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,1)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,1)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

        ! meridional(v)-component of Green-Gauss gradient
        p_grad(2,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,2)*p_cc(jc,jk,jb)    + &
          ptr_int%geofac_grg(jc,2,jb,2)*p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
          ptr_int%geofac_grg(jc,3,jb,2)*p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
          ptr_int%geofac_grg(jc,4,jb,2)*p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA

END SUBROUTINE grad_green_gauss_cell_adv

SUBROUTINE grad_green_gauss_cell_dycore(p_ccpr, ptr_patch, ptr_int, p_grad,         &
    &                                   opt_slev, opt_elev, opt_rlstart, opt_rlend, &
    &                                   opt_acc_async)
  !
  !  patch on which computation is performed
  !
  TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
  !
  !  data structure for interpolation
  !
  TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

  !  cell centered I/O variables
  !
  REAL(vp), INTENT(in) :: p_ccpr(:,:,:,:) ! perturbation fields passed from dycore (2,nproma,nlev,nblks_c)

  INTEGER, INTENT(in), OPTIONAL :: opt_slev    ! optional vertical start level

  INTEGER, INTENT(in), OPTIONAL :: opt_elev    ! optional vertical end level

  INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

  LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
  !
  ! cell based Green-Gauss reconstructed geographical gradient vector
  !
  REAL(vp), INTENT(inout) :: p_grad(:,:,:,:)      ! dim:(4,nproma,nlev,nblks_c)

  INTEGER :: slev, elev     ! vertical start and end level
  INTEGER :: jc, jk, jb
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
#ifdef __SWAPDIM
    elev = UBOUND(p_ccpr,2)
#else
    elev = UBOUND(p_ccpr,3)
#endif
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF

  iidx => ptr_patch%cells%neighbor_idx
  iblk => ptr_patch%cells%neighbor_blk

  i_nchdom = MAX(1,ptr_patch%n_childdom)

  !
  ! 2. reconstruction of cell based geographical gradient
  !

  !$ACC DATA PRESENT(p_ccpr, p_grad, ptr_int, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else

      !$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(8)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
#ifdef __SWAPDIM
          ! zonal(u)-component of Green-Gauss gradient, field 1
          p_grad(jc,jk,jb,1) = &
            ptr_int%geofac_grg(jc,1,jb,1)*p_ccpr(jc,jk,jb,1) + &
            ptr_int%geofac_grg(jc,2,jb,1)*p_ccpr(iidx(jc,jb,1),jk,iblk(jc,jb,1),1) + &
            ptr_int%geofac_grg(jc,3,jb,1)*p_ccpr(iidx(jc,jb,2),jk,iblk(jc,jb,2),1) + &
            ptr_int%geofac_grg(jc,4,jb,1)*p_ccpr(iidx(jc,jb,3),jk,iblk(jc,jb,3),1)
          ! meridional(v)-component of Green-Gauss gradient, field 1
          p_grad(jc,jk,jb,2) = &
            ptr_int%geofac_grg(jc,1,jb,2)*p_ccpr(jc,jk,jb,1) + &
            ptr_int%geofac_grg(jc,2,jb,2)*p_ccpr(iidx(jc,jb,1),jk,iblk(jc,jb,1),1) + &
            ptr_int%geofac_grg(jc,3,jb,2)*p_ccpr(iidx(jc,jb,2),jk,iblk(jc,jb,2),1) + &
            ptr_int%geofac_grg(jc,4,jb,2)*p_ccpr(iidx(jc,jb,3),jk,iblk(jc,jb,3),1)
          ! zonal(u)-component of Green-Gauss gradient, field 2
          p_grad(jc,jk,jb,3) = &
            ptr_int%geofac_grg(jc,1,jb,1)*p_ccpr(jc,jk,jb,2) + &
            ptr_int%geofac_grg(jc,2,jb,1)*p_ccpr(iidx(jc,jb,1),jk,iblk(jc,jb,1),2) + &
            ptr_int%geofac_grg(jc,3,jb,1)*p_ccpr(iidx(jc,jb,2),jk,iblk(jc,jb,2),2) + &
            ptr_int%geofac_grg(jc,4,jb,1)*p_ccpr(iidx(jc,jb,3),jk,iblk(jc,jb,3),2)
          ! meridional(v)-component of Green-Gauss gradient, field 2
          p_grad(jc,jk,jb,4) = &
            ptr_int%geofac_grg(jc,1,jb,2)*p_ccpr(jc,jk,jb,2) + &
            ptr_int%geofac_grg(jc,2,jb,2)*p_ccpr(iidx(jc,jb,1),jk,iblk(jc,jb,1),2) + &
            ptr_int%geofac_grg(jc,3,jb,2)*p_ccpr(iidx(jc,jb,2),jk,iblk(jc,jb,2),2) + &
            ptr_int%geofac_grg(jc,4,jb,2)*p_ccpr(iidx(jc,jb,3),jk,iblk(jc,jb,3),2)
#else
          ! zonal(u)-component of Green-Gauss gradient, field 1
          p_grad(1,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,1)*p_ccpr(1,jc,jk,jb)+     &
            ptr_int%geofac_grg(jc,2,jb,1)*p_ccpr(1,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,1)*p_ccpr(1,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,1)*p_ccpr(1,iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! meridional(v)-component of Green-Gauss gradient, field 1
          p_grad(2,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,2)*p_ccpr(1,jc,jk,jb)    + &
            ptr_int%geofac_grg(jc,2,jb,2)*p_ccpr(1,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,2)*p_ccpr(1,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,2)*p_ccpr(1,iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! zonal(u)-component of Green-Gauss gradient, field 2
          p_grad(3,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,1)*p_ccpr(2,jc,jk,jb)    + &
            ptr_int%geofac_grg(jc,2,jb,1)*p_ccpr(2,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,1)*p_ccpr(2,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,1)*p_ccpr(2,iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! meridional(v)-component of Green-Gauss gradient, field 2
          p_grad(4,jc,jk,jb) = ptr_int%geofac_grg(jc,1,jb,2)*p_ccpr(2,jc,jk,jb)    + &
            ptr_int%geofac_grg(jc,2,jb,2)*p_ccpr(2,iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
            ptr_int%geofac_grg(jc,3,jb,2)*p_ccpr(2,iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
            ptr_int%geofac_grg(jc,4,jb,2)*p_ccpr(2,iidx(jc,jb,3),jk,iblk(jc,jb,3))
#endif
        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    END DO ! end loop over blocks

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
  END SUBROUTINE grad_green_gauss_cell_dycore

END MODULE mo_math_gradients
