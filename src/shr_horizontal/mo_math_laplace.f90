!   Contains the implementation of the nabla mathematical operators.
!
!   Contains the implementation of the mathematical operators
!   employed by the shallow water prototype.
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

MODULE mo_math_laplace
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_intp,                ONLY: edges2verts_scalar, verts2edges_scalar
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: l_limited_area
USE mo_parallel_config,     ONLY: nproma, p_test_run
USE mo_exception,           ONLY: finish
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
USE mo_math_gradients,      ONLY: grad_fd_norm
USE mo_math_divrot,         ONLY: div, rot_vertex
USE mo_fortran_tools,       ONLY: copy


IMPLICIT NONE

PRIVATE


PUBLIC :: nabla2_vec
PUBLIC :: nabla2_scalar, nabla2_scalar_avg
PUBLIC :: nabla4_vec
PUBLIC :: nabla4_scalar

INTERFACE nabla2_vec

  MODULE PROCEDURE nabla2_vec_atmos

END INTERFACE


CONTAINS


!-------------------------------------------------------------------------
!
!>
!!  Computes  laplacian of a vector field.
!!
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
SUBROUTINE nabla2_vec_atmos( vec_e, ptr_patch, ptr_int, nabla2_vec_e, &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which laplacian is computed
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
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_c, rl_end_c, rl_start_v, rl_end_v
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlen, nblks_e, npromz_e

REAL(wp) ::  &
  &  z_div_c(nproma,ptr_patch%nlev,ptr_patch%nblks_c),  &
  &  z_rot_v(nproma,ptr_patch%nlev,ptr_patch%nblks_v)

INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk

!-----------------------------------------------------------------------
IF (p_test_run) THEN
  z_div_c(:,:,:)=0.0_wp
  z_rot_v(:,:,:)=0.0_wp
ENDIF

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
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 2)) THEN
    CALL finish ('mo_math_operators:nabla2_vec_atmos',  &
          &      'opt_rlstart must not be between 0 and 2')
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

rl_start_c = rl_start/2

IF (rl_start > 0) THEN
  rl_start_v = (rl_start+1)/2
ELSE
  rl_start_v = (rl_start-1)/2
ENDIF

IF (rl_end > 0) THEN
  rl_end_c = (rl_end+1)/2
  rl_end_v = rl_end/2+1
ELSE
  rl_end_c = (rl_end-1)/2
  rl_end_v = rl_end/2-1
ENDIF

rl_end_c = MAX(min_rlcell,rl_end_c)
rl_end_v = MAX(min_rlvert,rl_end_v)

icidx => ptr_patch%edges%cell_idx
icblk => ptr_patch%edges%cell_blk
ividx => ptr_patch%edges%vertex_idx
ivblk => ptr_patch%edges%vertex_blk

i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$ACC DATA CREATE(z_div_c, z_rot_v) PRESENT(vec_e) PRESENT(nabla2_vec_e) &
!$ACC   PRESENT(ptr_patch, ptr_int) IF(i_am_accel_node)

! Initialization of unused elements of nabla2_vec_e
! DO jb = 1, i_startblk
!   nabla2_vec_e(:,:,jb) = 0._wp
! ENDDO
! DO jb = i_endblk, ptr_patch%nblks_e
!   nabla2_vec_e(:,:,jb) = 0._wp
! ENDDO

! compute divergence of vector field
CALL div( vec_e, ptr_patch, ptr_int, z_div_c, slev, elev, &
          opt_rlstart=rl_start_c, opt_rlend=rl_end_c )

!
!  loop through over all patch edges (and blocks)
!

! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (ptr_patch%geometry_info%cell_type)

CASE (3) ! (cell_type == 3)

  ! compute rotation of vector field
  CALL rot_vertex( vec_e, ptr_patch, ptr_int, z_rot_v, slev, elev, &
                   opt_rlstart=rl_start_v, opt_rlend=rl_end_v)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
!CDIR UNROLL=3
   !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO je = i_startidx, i_endidx

        nabla2_vec_e(je,jk,jb) =  &
          &   ptr_patch%edges%tangent_orientation(je,jb) *  &
          &   ( z_rot_v(ividx(je,jb,2),jk,ivblk(je,jb,2))  &
          &   - z_rot_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) )  &
          &   * ptr_patch%edges%inv_primal_edge_length(je,jb)  &
          & + ( z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))    &
          &   - z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1)) )  &
          &   * ptr_patch%edges%inv_dual_edge_length(je,jb)

      END DO
    END DO
    !$ACC END PARALLEL

  END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SELECT
!$ACC WAIT(1)

!$ACC END DATA

END SUBROUTINE nabla2_vec_atmos


!-------------------------------------------------------------------------
!

!>
!! Computes biharmonic laplacian @f$\nabla ^4@f$ of a vector field without boundaries as used in atmospheric model.
!!
!! Computes biharmonic laplacian @f$\nabla ^4@f$ of a vector field without boundaries as used in atmospheric model.
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
SUBROUTINE nabla4_vec( vec_e, ptr_patch, ptr_int, nabla4_vec_e, &
  &                    opt_nabla2, opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which laplacian is computed
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
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla4_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Optional argument for passing nabla2 to the calling program
! (to avoid double computation for Smagorinsky diffusion and nest boundary diffusion)
REAL(wp), INTENT(inout), TARGET, OPTIONAL  ::  &
  &  opt_nabla2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_s1, rl_end_s1

REAL(wp), ALLOCATABLE, TARGET :: z_nabla2_vec_e(:,:,:) ! dim: (nproma,nlev,ptr_patch%nblks_e)
REAL(wp), POINTER :: p_nabla2(:,:,:)

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
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 4)) THEN
    CALL finish ('mo_math_operators:nabla4_vec',  &
          &      'opt_rlstart must not be between 0 and 4')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 5
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

IF (rl_start > 0) THEN
  rl_start_s1 = rl_start - 2
ELSE
  rl_start_s1 = rl_start + 2
ENDIF

IF (rl_end > 0) THEN
  rl_end_s1 = rl_end + 2
ELSE
  rl_end_s1 = rl_end - 2
ENDIF

rl_end_s1 = MAX(min_rledge,rl_end_s1)

IF (PRESENT(opt_nabla2) ) THEN
  p_nabla2 => opt_nabla2
ELSE
  ALLOCATE (z_nabla2_vec_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e))

  p_nabla2 => z_nabla2_vec_e
ENDIF

!$ACC DATA CREATE(z_nabla2_vec_e) PRESENT(vec_e) PRESENT(nabla4_vec_e) &
!$ACC   PRESENT(ptr_patch, ptr_int) IF(i_am_accel_node)

!
! apply second order Laplacian twice
!
IF (p_test_run) THEN
  p_nabla2(:,:,:) = 0.0_wp
!   rl_start_s1 = 1
!   rl_end_s1 = min_rledge
ENDIF

CALL nabla2_vec( vec_e, ptr_patch, ptr_int, p_nabla2,  &
  &              slev, elev, opt_rlstart=rl_start_s1, opt_rlend=rl_end_s1 )

CALL sync_patch_array(SYNC_E, ptr_patch, p_nabla2)

CALL nabla2_vec( p_nabla2, ptr_patch, ptr_int, nabla4_vec_e,  &
  &              slev, elev, opt_rlstart=rl_start, opt_rlend=rl_end )

IF (.NOT. PRESENT(opt_nabla2) ) THEN
  DEALLOCATE (z_nabla2_vec_e)
ENDIF

!$ACC END DATA

END SUBROUTINE nabla4_vec
!-----------------------------------------------------------------------

!>
!!  Computes laplacian @f$\nabla ^2 @f$ of a scalar field.
!!
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
SUBROUTINE nabla2_scalar( psi_c, ptr_patch, ptr_int, nabla2_psi_c, &
  &                       slev, elev, rl_start, rl_end )

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch         !< patch on which computation is performed
TYPE(t_int_state),     INTENT(in) :: ptr_int           !< interpolation state

REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) !< cells based variable of which biharmonic laplacian is computed, dim: (nproma,nlev,nblks_c)

INTEGER,               INTENT(in) :: slev              !< vertical start level
INTEGER,               INTENT(in) :: elev              !< vertical end level
INTEGER,               INTENT(in) :: rl_start,rl_end   !< start and end values of refin_ctrl flag

! cell based variable in which biharmonic laplacian is stored
REAL(wp), INTENT(inout) ::  &
  &  nabla2_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)
INTEGER :: jb, jc, jk, i_startblk, i_endblk, i_startidx, i_endidx
REAL(wp) ::  &
  &  z_grad_fd_norm_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e)
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

!$ACC DATA CREATE(z_grad_fd_norm_e) PRESENT(psi_c) PRESENT(nabla2_psi_c) &
!$ACC   PRESENT(ptr_patch, ptr_int) IF(i_am_accel_node)

SELECT CASE (ptr_patch%geometry_info%cell_type)

CASE (3) ! (cell_type == 3)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)







    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        !
        !  calculate div(grad) in one step
        !
        nabla2_psi_c(jc,jk,jb) =  &
          &    psi_c(jc,jk,jb)                       * ptr_int%geofac_n2s(jc,1,jb) &
          &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_n2s(jc,2,jb) &
          &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_n2s(jc,3,jb) &
          &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_n2s(jc,4,jb)
      END DO !cell loop
    END DO !vertical levels loop
    !$ACC END PARALLEL

  END DO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

CASE (6) ! (cell_type == 6) THEN ! Use unoptimized version for the time being

  ! compute finite difference gradient in normal direction
  CALL grad_fd_norm( psi_c, ptr_patch, z_grad_fd_norm_e, slev, elev)

  ! compute divergence of resulting vector field
  CALL div( z_grad_fd_norm_e, ptr_patch, ptr_int, nabla2_psi_c, slev, elev)

END SELECT
!$ACC WAIT(1)

!$ACC END DATA

END SUBROUTINE nabla2_scalar

!-------------------------------------------------------------------------
!

!>
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging.
!!
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging
!!  with the neighboring cells to increase computing efficiency.
!!  NOTE: This optimized routine works for triangular grids only.
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
SUBROUTINE nabla2_scalar_avg( psi_c, ptr_patch, ptr_int, avg_coeff, nabla2_psi_c, &
  &                           opt_slev, opt_elev )
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  averaging coefficients
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)

!
!  cells based variable of which biharmonic laplacian is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level


!
!  cell based variable in which biharmonic laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end, rl_start_l2
INTEGER :: jb, jc, jk, i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

REAL(wp), DIMENSION (nproma,ptr_patch%nlev,ptr_patch%nblks_c) :: aux_c

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

iidx => ptr_patch%cells%neighbor_idx
iblk => ptr_patch%cells%neighbor_blk

rl_start = 2
rl_start_l2 = rl_start + 1
rl_end = min_rlcell
i_nchdom   = MAX(1,ptr_patch%n_childdom)

! The special treatment of 2D fields is essential for efficiency on the NEC

!$ACC DATA CREATE(aux_c) PRESENT(avg_coeff, psi_c) PRESENT(nabla2_psi_c) &
!$ACC   PRESENT(ptr_patch, ptr_int) IF(i_am_accel_node)

SELECT CASE (ptr_patch%geometry_info%cell_type)

CASE (3) ! (cell_type == 3)

IF (slev == elev) THEN
  jk = slev

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    IF (jb == i_startblk) THEN
      i_startidx = ptr_patch%cells%start_idx(rl_start,1)
      i_endidx   = nproma
      IF (jb == i_endblk) i_endidx = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE IF (jb == i_endblk) THEN
      i_startidx = 1
      i_endidx   = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE
      i_startidx = 1
      i_endidx   = nproma
    ENDIF

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx

      !
      !  calculate div(grad) in one step
      !
      aux_c(jc,jk,jb) =  &
        &    psi_c(jc,jk,jb)                       * ptr_int%geofac_n2s(jc,1,jb) &
        &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_n2s(jc,2,jb) &
        &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_n2s(jc,3,jb) &
        &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_n2s(jc,4,jb)

    END DO
    !$ACC END PARALLEL

  END DO
!$OMP END DO

  IF (l_limited_area .OR. ptr_patch%id > 1) THEN
    ! Fill nabla2_psi_c along the lateral boundaries of nests

    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)

    CALL copy(aux_c(:,jk,i_startblk:i_endblk), &
         nabla2_psi_c(:,jk,i_startblk:i_endblk))
!$OMP BARRIER
  ENDIF

!
! Now do averaging with weights given by avg_coeff

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    IF (jb == i_startblk) THEN
      i_startidx = ptr_patch%cells%start_idx(rl_start_l2,1)
      i_endidx   = nproma
      IF (jb == i_endblk) i_endidx = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE IF (jb == i_endblk) THEN
      i_startidx = 1
      i_endidx   = ptr_patch%cells%end_idx(rl_end,i_nchdom)
    ELSE
      i_startidx = 1
      i_endidx   = nproma
    ENDIF

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP VECTOR
    DO jc = i_startidx, i_endidx
      !
      !  calculate the weighted average
      !
      nabla2_psi_c(jc,jk,jb) =  &
        &    aux_c(jc,jk,jb)                       * avg_coeff(jc,1,jb) &
        &  + aux_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
        &  + aux_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
        &  + aux_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * avg_coeff(jc,4,jb)

    END DO !cell loop
    !$ACC END PARALLEL

  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

ELSE

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
        !
        !  calculate div(grad) in one step
        !
        aux_c(jc,jk,jb) =  &
          &    psi_c(jc,jk,jb)                       * ptr_int%geofac_n2s(jc,1,jb) &
          &  + psi_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_n2s(jc,2,jb) &
          &  + psi_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_n2s(jc,3,jb) &
          &  + psi_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_n2s(jc,4,jb)
      END DO !cell loop
    END DO !vertical levels loop
    !$ACC END PARALLEL

  END DO
!$OMP END DO

  IF (l_limited_area .OR. ptr_patch%id > 1) THEN
    ! Fill nabla2_psi_c along the lateral boundaries of nests

    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)

    CALL copy(aux_c(:,:,i_startblk:i_endblk), &
         nabla2_psi_c(:,:,i_startblk:i_endblk))
!$OMP BARRIER
  ENDIF

!
! Now do averaging with weights given by avg_coeff

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start_l2, rl_end)

    !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG
    DO jk = slev, elev
      !$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx
        !
        !  calculate the weighted average
        !
        nabla2_psi_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                       * avg_coeff(jc,1,jb) &
          &  + aux_c(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop
    END DO !vertical levels loop
    !$ACC END PARALLEL

  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL
ENDIF

END SELECT
!$ACC WAIT(1)

!$ACC END DATA

END SUBROUTINE nabla2_scalar_avg

!-------------------------------------------------------------------------
!
!

!>
!!  Computes biharmonic laplacian @f$\nabla ^4 @f$ of a scalar field.
!!
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
SUBROUTINE nabla4_scalar( psi_c, ptr_patch, ptr_int, nabla4_psi_c, &
  &                       slev, elev, rl_start, rl_end, p_nabla2  )

TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch           !< patch on which computation is performed
TYPE(t_int_state),     INTENT(in)    :: ptr_int             !< interpolation state

REAL(wp),              INTENT(in)    ::  &
  &  psi_c(:,:,:) !< cells based variable of which biharmonic laplacian is computed, dim: (nproma,nlev,nblks_c)

INTEGER,               INTENT(in)    ::  slev               !< vertical start level
INTEGER,               INTENT(in)    ::  elev               !< vertical end level
INTEGER,               INTENT(in)    ::  rl_start, rl_end   !< start and end values of refin_ctrl flag

! cell based variable in which biharmonic laplacian is stored
REAL(wp), INTENT(inout) ::  nabla4_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)
! argument for passing nabla2 to the calling program
! (to avoid double computation for Smagorinsky diffusion and nest boundary diffusion)
REAL(wp), INTENT(inout)  ::  &
  &  p_nabla2(:,:,:) ! dim: (nproma,nlev,nblks_e)
INTEGER :: rl_start_s1, rl_end_s1

!-----------------------------------------------------------------------

IF (rl_start > 0) THEN
  rl_start_s1 = rl_start - 1
ELSE
  rl_start_s1 = rl_start + 1
ENDIF
IF (rl_end > 0) THEN
  rl_end_s1 = rl_end + 1
ELSE
  rl_end_s1 = rl_end - 1
ENDIF

rl_end_s1 = MAX(min_rlcell,rl_end_s1)

!$ACC DATA PRESENT(psi_c) PRESENT(nabla4_psi_c) &
!$ACC   PRESENT(ptr_patch, ptr_int) IF(i_am_accel_node)

! apply second order Laplacian twice
IF (p_test_run) p_nabla2(:,:,:) = 0.0_wp

CALL nabla2_scalar( psi_c, ptr_patch, ptr_int, p_nabla2, &
                    slev, elev, rl_start=rl_start_s1, rl_end=rl_end_s1 )

CALL sync_patch_array(SYNC_C, ptr_patch, p_nabla2)

CALL nabla2_scalar( p_nabla2, ptr_patch, ptr_int, nabla4_psi_c, &
                    slev, elev, rl_start=rl_start, rl_end=rl_end )

!$ACC END DATA

END SUBROUTINE nabla4_scalar


END MODULE mo_math_laplace
