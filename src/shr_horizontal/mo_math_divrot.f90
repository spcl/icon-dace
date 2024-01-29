! Contains the implementation of the div,rot,recon mathematical
! operators employed by the shallow water prototype.
!
! @par To Do
! Boundary exchange, nblks in presence of halos and dummy edge
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

MODULE mo_math_divrot
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp, vp
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_intp_data_strc,      ONLY: t_int_state, t_lsq
USE mo_interpol_config,     ONLY: lsq_high_set
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: l_limited_area
USE mo_parallel_config,     ONLY: nproma
USE mo_exception,           ONLY: finish
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_fortran_tools,       ONLY: init
#ifdef _OPENACC
USE mo_mpi,                 ONLY: i_am_accel_node
#endif

! USE mo_timer,              ONLY: timer_start, timer_stop, timer_div

IMPLICIT NONE

PRIVATE


PUBLIC :: recon_lsq_cell_l, recon_lsq_cell_l_svd
PUBLIC :: recon_lsq_cell_q, recon_lsq_cell_q_svd
PUBLIC :: recon_lsq_cell_c, recon_lsq_cell_c_svd
PUBLIC :: div, div_avg
PUBLIC :: rot_vertex, rot_vertex_ri
PUBLIC :: rot_vertex_atmos

INTERFACE rot_vertex

  MODULE PROCEDURE rot_vertex_atmos

END INTERFACE


INTERFACE div

  MODULE PROCEDURE div3d
  MODULE PROCEDURE div3d_2field
  MODULE PROCEDURE div4d

END INTERFACE


CONTAINS


!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Rx = Q^T d.
!! R: upper triangular matrix (2 x 2)
!! Q: orthogonal matrix (3 x 2)
!! d: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! works only on triangular grid yet
!!
SUBROUTINE recon_lsq_cell_l( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                          opt_slev, opt_elev, opt_rlstart,       &
  &                          opt_rlend, opt_lconsv, opt_acc_async )

  TYPE(t_patch), TARGET, INTENT(IN) :: &  !< patch on which computation 
    &  ptr_patch                          !<is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN)          ::  &   !<  cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< if true, conservative reconstruction is used
    &  opt_lconsv

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< optional async OpenACC
    &  opt_acc_async 
#ifdef __LOOP_EXCHANGE
  REAL(wp)  ::   &               !< weights * difference of scalars i j
    &  z_d(3,nproma,ptr_patch%nlev)
#else
  REAL(wp)  :: z_d(3)
#endif
  REAL(wp)  ::   &               !< matrix product of transposed Q matrix and d
    &  z_qt_times_d(2)

  INTEGER, POINTER ::   &             !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)       !< required stencil
  INTEGER :: slev, elev               !< vertical start and end level
  INTEGER :: jc, jk, jb               !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
  LOGICAL :: l_consv

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
  IF ( PRESENT(opt_lconsv) ) THEN
    l_consv = opt_lconsv
  ELSE
    l_consv = .FALSE.
  END IF


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  ! pointers to line and block indices of stencil
  iidx => ptr_patch%cells%neighbor_idx
  iblk => ptr_patch%cells%neighbor_blk



  !
  ! 1. reconstruction of cell based gradient (geographical components)
  !
  !$ACC DATA PRESENT(p_coeff, p_cc, iidx, iblk, ptr_int_lsq)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    !$ACC DATA CREATE(z_d)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev

        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the multiplication with
        ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
        ! Q-matrix.
        z_d(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_qt_times_d)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! matrix multiplication Q^T d (partitioned into 2 dot products)
        z_qt_times_d(1) = ptr_int_lsq%lsq_qtmat_c(jc,1,1,jb) * z_d(1,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,2,jb) * z_d(2,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,3,jb) * z_d(3,jc,jk)
        z_qt_times_d(2) = ptr_int_lsq%lsq_qtmat_c(jc,2,1,jb) * z_d(1,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,2,jb) * z_d(2,jc,jk)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,3,jb) * z_d(3,jc,jk)


        ! Solve linear system by backward substitution
        ! Gradient in zonal and meridional direction
        !
        ! meridional
        p_coeff(3,jc,jk,jb) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,2,jb) * z_qt_times_d(2)

        ! zonal
        p_coeff(2,jc,jk,jb) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,1,jb)                  &
          & * (z_qt_times_d(1) - ptr_int_lsq%lsq_rmat_utri_c(jc,1,jb)                &
          & * p_coeff(3,jc,jk,jb))

        ! constant
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)
        
      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

#else
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_d, z_qt_times_d)
!$NEC outerloop_unroll(4)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the multiplication with
        ! lsq_weights_c(jc,js,jb) has been shifted into the transposed
        ! Q-matrix.
        z_d(1) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

        ! matrix multiplication Q^T d (partitioned into 2 dot products)
        z_qt_times_d(1) = ptr_int_lsq%lsq_qtmat_c(jc,1,1,jb) * z_d(1)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,2,jb) * z_d(2)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,1,3,jb) * z_d(3)
        z_qt_times_d(2) = ptr_int_lsq%lsq_qtmat_c(jc,2,1,jb) * z_d(1)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,2,jb) * z_d(2)  &
          &             + ptr_int_lsq%lsq_qtmat_c(jc,2,3,jb) * z_d(3)


        ! Solve linear system by backward substitution
        ! Gradient in zonal and meridional direction
        !
        ! meridional
        p_coeff(3,jc,jk,jb) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,2,jb) * z_qt_times_d(2)

        ! zonal
        p_coeff(2,jc,jk,jb) = ptr_int_lsq%lsq_rmat_rdiag_c(jc,1,jb)                  &
          & * (z_qt_times_d(1) - ptr_int_lsq%lsq_rmat_utri_c(jc,1,jb)                &
          & * p_coeff(3,jc,jk,jb))

        ! constant
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)
        
      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

#endif


    IF (l_consv) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          ! constant
          p_coeff(1,jc,jk,jb) = p_coeff(1,jc,jk,jb)                                    &
            &                 - p_coeff(2,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,1) &
            &                 - p_coeff(3,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,2)


        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    ENDIF


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

END SUBROUTINE recon_lsq_cell_l



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Ax = b via Singular 
!! Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape:
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (2 x 3)
!! b: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! only works on triangular grid yet
!!
SUBROUTINE recon_lsq_cell_l_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff,      &
  &                              opt_slev, opt_elev, opt_rlstart, opt_rlend, &
  &                              opt_lconsv, opt_acc_async )

  TYPE(t_patch), TARGET, INTENT(IN) :: &  !< patch on which computation 
    &  ptr_patch                          !< is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN)          ::  &   !<  cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< if true, conservative reconstruction is used
    &  opt_lconsv

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< optional async OpenACC
    &  opt_acc_async 
#ifdef __LOOP_EXCHANGE
  REAL(wp)  ::   &               !< weights * difference of scalars i j
    &  z_b(3,nproma,ptr_patch%nlev)
#else
  REAL(wp)  ::  z_b(3)
#endif
  INTEGER, POINTER ::   &            !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)      !< required stencil
  INTEGER :: slev, elev              !< vertical start and end level
  INTEGER :: jc, jk, jb              !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
  LOGICAL :: l_consv

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
  IF ( PRESENT(opt_lconsv) ) THEN
    l_consv = opt_lconsv
  ELSE
    l_consv = .FALSE.
  END IF


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  ! pointers to line and block indices of stencil
  iidx => ptr_patch%cells%neighbor_idx
  iblk => ptr_patch%cells%neighbor_blk


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    !$ACC DATA CREATE(z_b)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev


        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the weights have been shifted 
        ! into the pseudoinverse.
        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

    !
    ! 2. compute cell based coefficients for linear reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! meridional
        p_coeff(3,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,2,1,jb) * z_b(1,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,2,jb) * z_b(2,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,3,jb) * z_b(3,jc,jk)

        ! zonal
        p_coeff(2,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,1,1,jb) * z_b(1,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,2,jb) * z_b(2,jc,jk)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,3,jb) * z_b(3,jc,jk)

        ! constant
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

#else
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_b)
!$NEC outerloop_unroll(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! note that the multiplication with lsq_weights_c(jc,js,jb) at
        ! runtime is now avoided. Instead, the weights have been shifted 
        ! into the pseudoinverse.
        z_b(1) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)


        ! meridional
        p_coeff(3,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,2,1,jb) * z_b(1)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,2,jb) * z_b(2)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,2,3,jb) * z_b(3)

        ! zonal
        p_coeff(2,jc,jk,jb) = ptr_int_lsq%lsq_pseudoinv(jc,1,1,jb) * z_b(1)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,2,jb) * z_b(2)  &
          &                 + ptr_int_lsq%lsq_pseudoinv(jc,1,3,jb) * z_b(3)

        ! constant
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)

      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL
#endif

    IF (l_consv) THEN

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! In the case of a conservative reconstruction, 
          ! the coefficient c0 is derived from the linear constraint
          !
          p_coeff(1,jc,jk,jb) = p_coeff(1,jc,jk,jb)                                    &
            &                 - p_coeff(2,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,1) &
            &                 - p_coeff(3,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,2)

        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      !$ACC END PARALLEL

    ENDIF

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

END SUBROUTINE recon_lsq_cell_l_svd


!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes the coefficients (derivatives) for a quadratic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (5 x 5)
!! Q  : orthogonal matrix (9 x 5)
!! Q^T: transposed of Q (5 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_q( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,  &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) ::   & !< patch on which computation
    &  ptr_patch                   !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< derivatives, one mixed derivative and a constant
                                  !< coefficient for zonal and meridional direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_d(lsq_high_set%dim_c,nproma,ptr_patch%nlev)
  REAL(wp)  ::           &        !< matrix-vector product of transposed
    &  z_qt_times_d(5)            !< Q matrix and d

  REAL(wp), POINTER ::   &        !< Pointer to reciprocal diagonal R-matrix-elements
    &  ptr_rrdiag(:,:,:)
  REAL(wp), POINTER ::   &        !< Pointer to upper triangular R-matrix-elements
    &  ptr_rutri(:,:,:)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

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


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c

  ! pointer to reciprocal diagonal R-elements
  ptr_rrdiag => ptr_int_lsq%lsq_rmat_rdiag_c(:,:,:)

  ! pointer to upper triangular R-elements
  ptr_rutri => ptr_int_lsq%lsq_rmat_utri_c(:,:,:)


  !$ACC DATA PRESENT(p_cc, p_coeff, ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_qtmat_c, iidx, iblk) &
  !$ACC   CREATE(z_d) IF(i_am_accel_node)
!$OMP PARALLEL

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,1:6,1:i_startblk))
!$OMP BARRIER
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !
    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
!$NEC outerloop_unroll(4)
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
#endif

        z_d(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_d(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_d(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_d(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_d(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_d(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_d(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG
    DO jk = slev, elev

      !$ACC LOOP VECTOR PRIVATE(z_qt_times_d)
      DO jc = i_startidx, i_endidx

        ! calculate matrix vector product Q^T d (transposed of Q times LHS)
        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
        z_qt_times_d(1) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,1,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(2) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,2,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(3) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,3,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(4) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,4,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(5) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,5,1:9,jb),z_d(1:9,jc,jk))

        !
        ! Solve linear system Rx=Q^T d by back substitution
        !
        p_coeff(6,jc,jk,jb) = ptr_rrdiag(jc,5,jb) * z_qt_times_d(5)
        p_coeff(5,jc,jk,jb) = ptr_rrdiag(jc,4,jb)                                         &
          &                 * ( z_qt_times_d(4) - ptr_rutri(jc,1,jb)*p_coeff(6,jc,jk,jb) )
        p_coeff(4,jc,jk,jb) = ptr_rrdiag(jc,3,jb)                                         &
          &                 * ( z_qt_times_d(3) - ptr_rutri(jc,2,jb)*p_coeff(5,jc,jk,jb)  &
          &                 - ptr_rutri(jc,3,jb) * p_coeff(6,jc,jk,jb) )
        p_coeff(3,jc,jk,jb) = ptr_rrdiag(jc,2,jb)                                         &
          &                 * ( z_qt_times_d(2) - ptr_rutri(jc,4,jb)*p_coeff(4,jc,jk,jb)  &
          &                 - ptr_rutri(jc,5,jb) * p_coeff(5,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,6,jb) * p_coeff(6,jc,jk,jb) )
        p_coeff(2,jc,jk,jb) = ptr_rrdiag(jc,1,jb)                                         &
          &                 * ( z_qt_times_d(1) - ptr_rutri(jc,7,jb)*p_coeff(3,jc,jk,jb)  &
          &                 - ptr_rutri(jc,8,jb) * p_coeff(4,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,9,jb) * p_coeff(5,jc,jk,jb)                    &
          &                 - ptr_rutri(jc,10,jb)* p_coeff(6,jc,jk,jb) )

        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb)                                              &
          &                  - p_coeff(2,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  - p_coeff(3,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  - p_coeff(4,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  - p_coeff(5,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  - p_coeff(6,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,5)


      END DO ! end loop over cells

    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC WAIT(1)
  !$ACC END DATA


END SUBROUTINE recon_lsq_cell_q



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a quadratic polynomial,
!! using the least-squares method. The coefficients are provided at cell 
!! centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (5 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients:
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_q_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,      &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) ::   & !< patch on which computation
    &  ptr_patch                   !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< derivatives, one mixed derivative and a constant
                                  !< coefficient for zonal and meridional direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_b(lsq_high_set%dim_c,nproma,ptr_patch%nlev)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

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


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c



  !$ACC DATA PRESENT(p_cc, p_coeff, ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_pseudoinv, iidx, iblk) &
  !$ACC   CREATE(z_b) IF(i_am_accel_node)
!$OMP PARALLEL

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,1:6,1:i_startblk))
!$OMP BARRIER
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! 1. compute right hand side of linear system
    !

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
!$NEC outerloop_unroll(4)
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
#endif

        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_b(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_b(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_b(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_b(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_b(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_b(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG
    DO jk = slev, elev

      !$ACC LOOP VECTOR
!$NEC ivdep
      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is used.
        p_coeff(6,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(5,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(4,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(3,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(2,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))

        ! At the end, the coefficient c0 is derived from the linear constraint
        !
        p_coeff(1,jc,jk,jb) = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(2:6,jc,jk,jb), &
          &                   ptr_int_lsq%lsq_moments(jc,jb,1:5))


      END DO ! end loop over cells

    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
  !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA

END SUBROUTINE recon_lsq_cell_q_svd



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes the coefficients (derivatives) for a cubic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (9 x 9)
!! Q  : orthogonal matrix (9 x 9)
!! Q^T: transposed of Q (9 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_c( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,  &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_d(lsq_high_set%dim_c,nproma,ptr_patch%nlev)
  REAL(wp)  ::           &        !< matrix-vector product of transposed
    &  z_qt_times_d(9)            !< Q matrix and d

  REAL(wp), POINTER ::   &        !< Pointer to reciprocal diagonal R-matrix-elements
    &  ptr_rrdiag(:,:,:)
  REAL(wp), POINTER ::   &        !< Pointer to upper triangular R-matrix-elements
    &  ptr_rutri(:,:,:)

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

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


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c

  ! pointer to reciprocal diagonal R-elements
  ptr_rrdiag => ptr_int_lsq%lsq_rmat_rdiag_c(:,:,:)

  ! pointer to upper triangular R-elements
  ptr_rutri => ptr_int_lsq%lsq_rmat_utri_c(:,:,:)


  !$ACC DATA PRESENT(p_cc, p_coeff, ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_qtmat_c, iidx, iblk) &
  !$ACC   CREATE(z_d) IF(i_am_accel_node)
!$OMP PARALLEL

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    CALL init(p_coeff(:,:,1:10,1:i_startblk))
!$OMP BARRIER
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_d,z_qt_times_d), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !
    ! 1. compute right hand side of linear system
    !
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
      !$ACC LOOP VECTOR
      DO jk = slev, elev
#else
    !$ACC LOOP GANG
!$NEC outerloop_unroll(4)
    DO jk = slev, elev
      !$ACC LOOP VECTOR
!NEC$ ivdep
      DO jc = i_startidx, i_endidx
#endif

        z_d(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_d(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_d(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_d(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_d(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_d(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_d(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_d(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_d(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !
    ! 2. compute cell based coefficients for quadratic reconstruction
    !
    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG
    DO jk = slev, elev

      !$ACC LOOP VECTOR PRIVATE(z_qt_times_d)
!NEC$ ivdep
      DO jc = i_startidx, i_endidx

        ! calculate matrix vector product Q^T d (transposed of Q times LHS)
        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied
!TODO:  these should be nine scalars, since they should reside in registers
        z_qt_times_d(1) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,1,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(2) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,2,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(3) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,3,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(4) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,4,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(5) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,5,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(6) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,6,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(7) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,7,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(8) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,8,1:9,jb),z_d(1:9,jc,jk))
        z_qt_times_d(9) = DOT_PRODUCT(ptr_int_lsq%lsq_qtmat_c(jc,9,1:9,jb),z_d(1:9,jc,jk))


        !
        ! Solve linear system Rx=Q^T d by back substitution
        !
        p_coeff(10,jc,jk,jb) = ptr_rrdiag(jc,9,jb) * z_qt_times_d(9)
        p_coeff(9,jc,jk,jb)  = ptr_rrdiag(jc,8,jb)                                         &
          &                  * ( z_qt_times_d(8) - ptr_rutri(jc,1,jb)*p_coeff(10,jc,jk,jb) )
        p_coeff(8,jc,jk,jb)  = ptr_rrdiag(jc,7,jb)                                         &
          &                  * ( z_qt_times_d(7) - (ptr_rutri(jc,2,jb)*p_coeff(9,jc,jk,jb) &
          &                  + ptr_rutri(jc,3,jb) * p_coeff(10,jc,jk,jb)) )
        p_coeff(7,jc,jk,jb)  = ptr_rrdiag(jc,6,jb)                                         &
          &                  * ( z_qt_times_d(6) - (ptr_rutri(jc,4,jb)*p_coeff(8,jc,jk,jb) &
          &                  + ptr_rutri(jc,5,jb) * p_coeff(9,jc,jk,jb)                    &
          &                  + ptr_rutri(jc,6,jb) * p_coeff(10,jc,jk,jb)) )
        p_coeff(6,jc,jk,jb)  = ptr_rrdiag(jc,5,jb)                                         &
          &                  * ( z_qt_times_d(5) - (ptr_rutri(jc,7,jb)*p_coeff(7,jc,jk,jb) &
          &                  + ptr_rutri(jc,8,jb) * p_coeff(8,jc,jk,jb)                    &
          &                  + ptr_rutri(jc,9,jb) * p_coeff(9,jc,jk,jb)                    &
          &                  + ptr_rutri(jc,10,jb)* p_coeff(10,jc,jk,jb)) )
        p_coeff(5,jc,jk,jb)  = ptr_rrdiag(jc,4,jb)                                         &
          &                  * ( z_qt_times_d(4) - (ptr_rutri(jc,11,jb)*p_coeff(6,jc,jk,jb)&
          &                  + ptr_rutri(jc,12,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,13,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,14,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,15,jb) * p_coeff(10,jc,jk,jb)) )
        p_coeff(4,jc,jk,jb)  = ptr_rrdiag(jc,3,jb)                                         &
          &                  * ( z_qt_times_d(3) - (ptr_rutri(jc,16,jb)*p_coeff(5,jc,jk,jb)&
          &                  + ptr_rutri(jc,17,jb) * p_coeff(6,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,18,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,19,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,20,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,21,jb) * p_coeff(10,jc,jk,jb)) )
        p_coeff(3,jc,jk,jb)  = ptr_rrdiag(jc,2,jb)                                         &
          &                  * ( z_qt_times_d(2) - (ptr_rutri(jc,22,jb)*p_coeff(4,jc,jk,jb)&
          &                  + ptr_rutri(jc,23,jb) * p_coeff(5,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,24,jb) * p_coeff(6,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,25,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,26,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,27,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,28,jb) * p_coeff(10,jc,jk,jb)) )
        p_coeff(2,jc,jk,jb)  = ptr_rrdiag(jc,1,jb)                                         &
          &                  * ( z_qt_times_d(1) - (ptr_rutri(jc,29,jb)*p_coeff(3,jc,jk,jb)&
          &                  + ptr_rutri(jc,30,jb) * p_coeff(4,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,31,jb) * p_coeff(5,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,32,jb) * p_coeff(6,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,33,jb) * p_coeff(7,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,34,jb) * p_coeff(8,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,35,jb) * p_coeff(9,jc,jk,jb)                   &
          &                  + ptr_rutri(jc,36,jb) * p_coeff(10,jc,jk,jb)) )


        p_coeff(1,jc,jk,jb)  = p_cc(jc,jk,jb) - (                                          &
          &                    p_coeff(2,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,1)     &
          &                  + p_coeff(3,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,2)     &
          &                  + p_coeff(4,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,3)     &
          &                  + p_coeff(5,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,4)     &
          &                  + p_coeff(6,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,5)     &
          &                  + p_coeff(7,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,6)     &
          &                  + p_coeff(8,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,7)     &
          &                  + p_coeff(9,jc,jk,jb)  * ptr_int_lsq%lsq_moments(jc,jb,8)     &
          &                  + p_coeff(10,jc,jk,jb) * ptr_int_lsq%lsq_moments(jc,jb,9))

      END DO ! end loop over cells

    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

  END DO ! end loop over blocks
  !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA


END SUBROUTINE recon_lsq_cell_c




!--------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a cubic polynomial,
!! using the least-squares method. The coefficients are provided at 
!! cell centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (9 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_c_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                               opt_slev, opt_elev, opt_rlstart,      &
  &                               opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed
                                        
  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

#ifdef __LOOP_EXCHANGE
  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_b(lsq_high_set%dim_c,nproma,ptr_patch%nlev)
#else
  REAL(wp)  ::           &        !< difference of scalars i j
    &  z_b(9)
#endif

  INTEGER, POINTER  ::   &        !< Pointer to line and block indices of
    &  iidx(:,:,:), iblk(:,:,:)   !< required stencil
  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: jc, jk, jb           !< index of cell, vertical level and block
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: rl_start_init, rl_end_init, i_startblk_init, i_endblk_init
  INTEGER :: ji

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


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  ! pointers to line and block indices of required stencil
  iidx => ptr_int_lsq%lsq_idx_c
  iblk => ptr_int_lsq%lsq_blk_c

  !$ACC DATA PRESENT(p_cc, p_coeff, ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_pseudoinv, iidx, iblk) &
  !$ACC   IF(i_am_accel_node)
!$OMP PARALLEL

  IF (ptr_patch%id > 1 .OR. l_limited_area) THEN
    ! Only zero-init the lateral boundary points

    ! values for the blocking
    rl_start_init   = 1
    rl_end_init     = MAX(1,rl_start-1)
    i_startblk_init = ptr_patch%cells%start_block(rl_start_init)
    i_endblk_init   = ptr_patch%cells%end_block(rl_end_init)

!$OMP DO PRIVATE(jb,jc,jk,ji,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk_init, i_endblk_init

      CALL get_indices_c(ptr_patch, jb, i_startblk_init, i_endblk_init, &
                         i_startidx, i_endidx, rl_start_init, rl_end_init)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
!NEC$ forced_collapse
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          DO ji = 1, 10
            p_coeff(ji,jc,jk,jb) = 0._wp
          ENDDO
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
!$OMP BARRIER
  ENDIF

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_b), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)
    !
    ! 1. compute right hand side of linear system
    !

#ifdef __LOOP_EXCHANGE
    !$ACC DATA CREATE(z_b)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev

        z_b(1,jc,jk) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2,jc,jk) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3,jc,jk) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_b(4,jc,jk) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_b(5,jc,jk) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_b(6,jc,jk) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_b(7,jc,jk) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_b(8,jc,jk) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_b(9,jc,jk) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !
    ! 2. compute cell based coefficients for cubic reconstruction
    !    calculate matrix vector product PINV(A) * b
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied

        p_coeff(10,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,9,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(9, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,8,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(8, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,7,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(7, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,6,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(6, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(5, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(4, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(3, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9,jc,jk))
        p_coeff(2, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9,jc,jk))


        p_coeff(1,jc,jk,jb)  = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(2:10,jc,jk,jb), &
          &                    ptr_int_lsq%lsq_moments(jc,jb,1:9))


      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

#else
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR TILE(32, 4) PRIVATE(z_b)
    DO jk = slev, elev
!$NEC ivdep
      DO jc = i_startidx, i_endidx

        !
        ! 1. compute right hand side of linear system
        !
        z_b(1) = p_cc(iidx(jc,jb,1),jk,iblk(jc,jb,1)) - p_cc(jc,jk,jb)
        z_b(2) = p_cc(iidx(jc,jb,2),jk,iblk(jc,jb,2)) - p_cc(jc,jk,jb)
        z_b(3) = p_cc(iidx(jc,jb,3),jk,iblk(jc,jb,3)) - p_cc(jc,jk,jb)
        z_b(4) = p_cc(iidx(jc,jb,4),jk,iblk(jc,jb,4)) - p_cc(jc,jk,jb)
        z_b(5) = p_cc(iidx(jc,jb,5),jk,iblk(jc,jb,5)) - p_cc(jc,jk,jb)
        z_b(6) = p_cc(iidx(jc,jb,6),jk,iblk(jc,jb,6)) - p_cc(jc,jk,jb)
        z_b(7) = p_cc(iidx(jc,jb,7),jk,iblk(jc,jb,7)) - p_cc(jc,jk,jb)
        z_b(8) = p_cc(iidx(jc,jb,8),jk,iblk(jc,jb,8)) - p_cc(jc,jk,jb)
        z_b(9) = p_cc(iidx(jc,jb,9),jk,iblk(jc,jb,9)) - p_cc(jc,jk,jb)

        !
        ! 2. compute cell based coefficients for cubic reconstruction
        !    calculate matrix vector product PINV(A) * b
        !

        ! (intrinsic function matmul not applied, due to massive
        ! performance penalty on the NEC. Instead the intrinsic dot product
        ! function is applied

        p_coeff(10,jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,9,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(9, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,8,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(8, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,7,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(7, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,6,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(6, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,5,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(5, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,4,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(4, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,3,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(3, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,2,1:9,jb), &
          &                               z_b(1:9))
        p_coeff(2, jc,jk,jb) = DOT_PRODUCT(ptr_int_lsq%lsq_pseudoinv(jc,1,1:9,jb), &
          &                               z_b(1:9))


        p_coeff(1,jc,jk,jb)  = p_cc(jc,jk,jb) - DOT_PRODUCT(p_coeff(2:10,jc,jk,jb), &
          &                    ptr_int_lsq%lsq_moments(jc,jb,1:9))


      END DO ! end loop over cells
    END DO ! end loop over vertical levels
    !$ACC END PARALLEL

#endif

  END DO ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC WAIT
  !$ACC END DATA

END SUBROUTINE recon_lsq_cell_c_svd


!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges.
!! The midpoint rule is used for quadrature.
!! input:  lives on edges (velocity points)
!! output: lives on centers of triangles
!!
SUBROUTINE div3d( vec_e, ptr_patch, ptr_int, div_vec_c, &
  &               opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

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
  elev = UBOUND(vec_e,2)
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

! loop through all patch cells (and blocks)
!

!IF(ltimer) CALL timer_start(timer_div)


  !$ACC DATA PRESENT(vec_e, div_vec_c, ptr_int%geofac_div, iidx, iblk) &
  !$ACC   IF(i_am_accel_node)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    ! original comment for divergence computation;
    ! everything that follows in this explanation has been combined into geofac_div

    ! compute the discrete divergence for cell jc by finite volume
    ! approximation (see Bonaventura and Ringler MWR 2005);
    ! multiplication of the normal vector component vec_e at the edges
    ! by the appropriate cell based edge_orientation is required to
    ! obtain the correct value for the application of Gauss theorem
    ! (which requires the scalar product of the vector field with the
    ! OUTWARD pointing unit vector with respect to cell jc; since the
    ! positive direction for the vector components is not necessarily
    ! the outward pointing one with respect to cell jc, a correction
    ! coefficient (equal to +-1) is necessary, given by
    ! ptr_patch%grid%cells%edge_orientation)

      !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO
      !$ACC END PARALLEL

  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div3d

SUBROUTINE div3d_2field( vec_e, ptr_patch, ptr_int, div_vec_c, &
  &                      opt_slev, opt_elev, in2, out2, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! second input field for more efficient processing in NH core
REAL(wp), INTENT(in) ::  &
  &  in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

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
  elev = UBOUND(vec_e,2)
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

! loop through all patch cells (and blocks)
!

!IF(ltimer) CALL timer_start(timer_div)


  !$ACC DATA PRESENT(vec_e, in2, div_vec_c, out2, ptr_int%geofac_div, iidx, iblk) &
  !$ACC   IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    ! original comment for divergence computation;
    ! everything that follows in this explanation has been combined into geofac_div

    ! compute the discrete divergence for cell jc by finite volume
    ! approximation (see Bonaventura and Ringler MWR 2005);
    ! multiplication of the normal vector component vec_e at the edges
    ! by the appropriate cell based edge_orientation is required to
    ! obtain the correct value for the application of Gauss theorem
    ! (which requires the scalar product of the vector field with the
    ! OUTWARD pointing unit vector with respect to cell jc; since the
    ! positive direction for the vector components is not necessarily
    ! the outward pointing one with respect to cell jc, a correction
    ! coefficient (equal to +-1) is necessary, given by
    ! ptr_patch%grid%cells%edge_orientation)

      !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO jc = i_startidx, i_endidx
        !$ACC LOOP VECTOR
        DO jk = slev, elev
#else
      !$ACC LOOP GANG
      DO jk = slev, elev
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

          out2(jc,jk,jb) =  &
            in2(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            in2(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            in2(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO
      !$ACC END PARALLEL

  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div3d_2field

!-------------------------------------------------------------------------
!
!
!>
!! Special version of div that processes 4D fields in one step
!!
!! See standard routine (div3d) for further description
!!
SUBROUTINE div4d( ptr_patch, ptr_int, f4din, f4dout, dim4d, &
  &              opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based 4D input field of which divergence is computed
!
REAL(wp), INTENT(in) ::  &
  &  f4din(:,:,:,:) ! dim: (nproma,nlev,nblks_e,dim4d)

INTEGER, INTENT(in) :: dim4d ! Last dimension of the input/output fields

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev(dim4d)    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev(dim4d)    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! 4D cell based variable in which divergence is stored
!
REAL(vp), INTENT(inout) ::  &
  &  f4dout(:,:,:,:) ! dim: (nproma,nlev,nblks_c,dim4d)

INTEGER :: slev(dim4d), elev(dim4d)     ! vertical start and end level
INTEGER :: jc, jk, jb, ji
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
  elev = UBOUND(f4din,2)
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

! loop through all patch cells (and blocks)
!

!IF(ltimer) CALL timer_start(timer_div)


  !$ACC DATA PRESENT(f4din, f4dout, ptr_int%geofac_div, iidx, iblk) &
  !$ACC   IF(i_am_accel_node)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,ji) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG
    DO jc = i_startidx, i_endidx
      DO ji = 1, dim4d
      !$ACC LOOP VECTOR
        DO jk = slev(ji), elev(ji)
#else
    !$ACC LOOP GANG
    DO ji = 1, dim4d
      DO jk = slev(ji), elev(ji)
        !$ACC LOOP VECTOR
        DO jc = i_startidx, i_endidx
#endif

          f4dout(jc,jk,jb,ji) =  &
            f4din(iidx(jc,jb,1),jk,iblk(jc,jb,1),ji) * ptr_int%geofac_div(jc,1,jb) + &
            f4din(iidx(jc,jb,2),jk,iblk(jc,jb,2),ji) * ptr_int%geofac_div(jc,2,jb) + &
            f4din(iidx(jc,jb,3),jk,iblk(jc,jb,3),ji) * ptr_int%geofac_div(jc,3,jb)

        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div4d

!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges,
!! followed by bilinear averaging to remove checkerboard noise
!! (Combines div_midpoint and cell_avg_varwgt to increase computing efficiency)
!!
SUBROUTINE div_avg( vec_e, ptr_patch, ptr_int, avg_coeff, div_vec_c,    &
  &                 opt_in2, opt_out2, opt_slev, opt_elev, opt_rlstart, &
  &                 opt_rlend )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  averaging coefficients
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)
!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! optional second input field for more efficient processing in NH core
REAL(wp), OPTIONAL, INTENT(in) ::  &
  &  opt_in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! optional second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  opt_out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: jc, jk, jb
INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end, rl_start_l2, rl_end_l1
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

REAL(wp), DIMENSION (nproma,ptr_patch%nlev,ptr_patch%nblks_c) :: aux_c, aux_c2

INTEGER,  DIMENSION(:,:,:),   POINTER :: inidx, inblk, ieidx, ieblk
LOGICAL :: l2fields

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
  elev = UBOUND(vec_e,2)
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
IF ( PRESENT(opt_in2) .AND. PRESENT(opt_out2)) THEN
  l2fields = .TRUE.
ELSE
  l2fields = .FALSE.
ENDIF

rl_start_l2 = rl_start + 1

IF ( PRESENT(opt_rlend) .AND. rl_end < 0 .AND. rl_end > min_rlcell ) THEN
  rl_end_l1 = rl_end - 1
ELSE
  rl_end_l1 = rl_end
END IF

inidx => ptr_patch%cells%neighbor_idx
inblk => ptr_patch%cells%neighbor_blk
ieidx => ptr_patch%cells%edge_idx
ieblk => ptr_patch%cells%edge_blk


! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)

! First compute divergence
!
!$ACC DATA PRESENT(vec_e, avg_coeff, div_vec_c) CREATE(aux_c) &
!$ACC   PRESENT(ptr_int, ieidx, ieblk, inidx, inblk) IF(i_am_accel_node)
!$ACC DATA PRESENT(opt_in2, opt_out2) CREATE(aux_c2) IF(i_am_accel_node .AND. l2fields)
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end_l1,i_nchdom)

IF (l2fields) THEN

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end_l1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)

#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        aux_c(jc,jk,jb) =  &
          vec_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          vec_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          vec_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        aux_c2(jc,jk,jb) =  &
          opt_in2(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          opt_in2(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          opt_in2(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

      END DO
    END DO
    !$ACC END PARALLEL

  END DO
!$OMP END DO

ELSE

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end_l1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif
        aux_c(jc,jk,jb) =  &
          vec_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          vec_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          vec_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

      END DO
    END DO
    !$ACC END PARALLEL

  END DO
!$OMP END DO

ENDIF

IF (l_limited_area .OR. ptr_patch%id > 1) THEN
  ! Fill div_vec_c along the lateral boundaries of nests

  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)

  !$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk ! like copy(aux_c, div_vec_c)

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_start_l2)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
        div_vec_c(jc,jk,jb) = aux_c(jc,jk,jb)
      END DO
    END DO
    !$ACC END PARALLEL
  END DO
  !$OMP END DO

  IF (l2fields) THEN
    !$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk ! like copy(aux_c2, opt_out2)

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_start_l2)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          opt_out2(jc,jk,jb) = aux_c2(jc,jk,jb)
        END DO
      END DO
    !$ACC END PARALLEL
    END DO
  !$OMP END DO
  ENDIF
ENDIF

!
! Now do averaging with weights given by avg_coeff

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!

IF (l2fields) THEN

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start_l2, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        !  calculate the weighted average
        div_vec_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

        opt_out2(jc,jk,jb) =  &
          &    aux_c2(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c2(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c2(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c2(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop
    END DO !vertical levels loop
    !$ACC END PARALLEL

  END DO !block loop

!$OMP END DO NOWAIT

ELSE

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start_l2, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jc = i_startidx, i_endidx
      DO jk = slev, elev
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
#endif

        !  calculate the weighted average
        div_vec_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop
    END DO !vertical levels loop
    !$ACC END PARALLEL

  END DO !block loop

!$OMP END DO NOWAIT

ENDIF

!$OMP END PARALLEL
!$ACC WAIT
!$ACC END DATA
!$ACC END DATA

END SUBROUTINE div_avg


!-------------------------------------------------------------------------
!
!>
!! Computes discrete rotation.
!!
!! Computes discrete rotation at
!! vertices of triangle cells (centers of dual hexagon cells)
!! from a vector field given by its components in the directions normal
!! to triangle edges.
!! input:  lives on edges (velocity points)
!! output: lives on dual of cells (vertices for triangular grid)
!!
SUBROUTINE rot_vertex_atmos( vec_e, ptr_patch, ptr_int, rot_vec, &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which rotation is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  vertex based variable in which rotation is stored
!
REAL(wp), INTENT(inout) ::  &
  &  rot_vec(:,:,:) ! dim: (nproma,nlev,nblks_v or nblks_e)

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
  elev = UBOUND(vec_e,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:rot_vertex_atmos',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

!
!  loop through over all patch vertices (and blocks)
!
! The special treatment of 2D fields is essential for efficiency on the NEC

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%verts%start_blk(rl_start,1)
  i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)

  !$ACC DATA PRESENT(vec_e, rot_vec, ptr_int%geofac_rot, iidx, iblk) &
  !$ACC   IF(i_am_accel_node)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! compute the discrete rotation for vertex jv by
    ! finite volume approximation
    ! (see Bonaventura and Ringler MWR 2005);
    ! multiplication of the vector component vec_e by
    ! the appropriate dual cell based verts%edge_orientation
    ! is required to obtain the correct value for the
    ! application of Stokes theorem (which requires the scalar
    ! product of the vector field with the tangent unit vectors
    ! going around dual cell jv COUNTERCLOKWISE;
    ! since the positive direction for the vec_e components is
    ! not necessarily the one yelding counterclockwise rotation
    ! around dual cell jv, a correction coefficient (equal to +-1)
    ! is necessary, given by g%verts%edge_orientation
    !

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif
        !
        ! calculate rotation, i.e.
        ! add individual edge contributions to rotation
        ! (remark: for pentagon points the 6th weighting is 0)
        !

        rot_vec(jv,jk,jb) =   &
          vec_e(iidx(jv,jb,1),jk,iblk(jv,jb,1)) * ptr_int%geofac_rot(jv,1,jb) + &
          vec_e(iidx(jv,jb,2),jk,iblk(jv,jb,2)) * ptr_int%geofac_rot(jv,2,jb) + &
          vec_e(iidx(jv,jb,3),jk,iblk(jv,jb,3)) * ptr_int%geofac_rot(jv,3,jb) + &
          vec_e(iidx(jv,jb,4),jk,iblk(jv,jb,4)) * ptr_int%geofac_rot(jv,4,jb) + &
          vec_e(iidx(jv,jb,5),jk,iblk(jv,jb,5)) * ptr_int%geofac_rot(jv,5,jb) + &
          vec_e(iidx(jv,jb,6),jk,iblk(jv,jb,6)) * ptr_int%geofac_rot(jv,6,jb)

      END DO

    END DO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA

END SUBROUTINE rot_vertex_atmos

!>
!! Same as above routine, but expects reversed index order (vertical first)
!! of the output field if __LOOP_EXCHANGE is specified. In addition, the 
!! output field (vorticity) has single precision if __MIXED_PRECISION is specified
!!
!!
SUBROUTINE rot_vertex_ri( vec_e, ptr_patch, ptr_int, rot_vec, &
  &                       opt_slev, opt_elev, opt_rlend, opt_acc_async )
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which rotation is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlend   ! end value of refin_ctrl flag

LOGICAL, INTENT(IN), OPTIONAL ::  &   
  &  opt_acc_async ! optional async OpenACC

!
!  vertex based variable in which rotation is stored
!
REAL(vp), INTENT(inout) ::  &
  &  rot_vec(:,:,:) ! dim: (nproma,nlev,nblks_v) or (nlev,nproma,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb

INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

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
  elev = UBOUND(vec_e,2)
END IF

rl_start = 2

IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF


!
!  loop through over all patch vertices (and blocks)
!

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  i_startblk = ptr_patch%verts%start_block(rl_start)
  i_endblk   = ptr_patch%verts%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), ICON_OMP_RUNTIME_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)

    ! calculate rotation, i.e.
    ! add individual edge contributions to rotation
    !
#ifdef __LOOP_EXCHANGE
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
        rot_vec(jk,jv,jb) =   &
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
        rot_vec(jv,jk,jb) =   &
#endif
          vec_e(iidx(jv,jb,1),jk,iblk(jv,jb,1)) * ptr_int%geofac_rot(jv,1,jb) + &
          vec_e(iidx(jv,jb,2),jk,iblk(jv,jb,2)) * ptr_int%geofac_rot(jv,2,jb) + &
          vec_e(iidx(jv,jb,3),jk,iblk(jv,jb,3)) * ptr_int%geofac_rot(jv,3,jb) + &
          vec_e(iidx(jv,jb,4),jk,iblk(jv,jb,4)) * ptr_int%geofac_rot(jv,4,jb) + &
          vec_e(iidx(jv,jb,5),jk,iblk(jv,jb,5)) * ptr_int%geofac_rot(jv,5,jb) + &
          vec_e(iidx(jv,jb,6),jk,iblk(jv,jb,6)) * ptr_int%geofac_rot(jv,6,jb)

      END DO
    END DO
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

END SUBROUTINE rot_vertex_ri

END MODULE mo_math_divrot
