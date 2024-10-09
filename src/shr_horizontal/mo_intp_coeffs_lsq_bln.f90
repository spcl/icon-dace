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

#ifdef __xlC__
! @PROCESS nosmp
! @PROCESS NOOPTimize
! @PROCESS smp=noopt
@PROCESS noopt
#endif
#ifdef __PGI
!pgi$g opt=1
#endif

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_intp_coeffs_lsq_bln
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
USE mo_math_constants,      ONLY: pi2
USE mo_exception,           ONLY: message, message_text, finish
USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell, min_rlcell_int
USE mo_master_control,      ONLY: get_my_process_type, wave_process
USE mo_model_domain,        ONLY: t_patch
USE mo_math_types,          ONLY: t_cartesian_coordinates
USE mo_math_utilities,      ONLY: gnomonic_proj, rotate_latlon, &
                                  plane_torus_closest_coordinates
USE mo_math_utility_solvers, ONLY: qrdec
USE mo_parallel_config,     ONLY: nproma
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_advection_config,    ONLY: advection_config
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_idx
USE mo_grid_config,         ONLY: grid_sphere_radius
USE mo_lib_grid_geometry_info,  ONLY: planar_torus_geometry, sphere_geometry
USE mo_intp_data_strc,      ONLY: t_lsq, t_int_state
USE mo_fortran_tools,       ONLY: copy

IMPLICIT NONE

PRIVATE

!> module name string
CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_coeffs_lsq_bln'

PUBLIC :: lsq_stencil_create
PUBLIC :: lsq_compute_coeff_cell
PUBLIC :: scalar_int_coeff
PUBLIC :: bln_int_coeff_e2c

CONTAINS

  !-------------------------------------------------------------------------
  !
  !! This routine initializes the indices used to define the stencil
  !! of the lsq reconstruction. The stencil is cell based and includes
  !! a variable number of cells (lsq_dim_c) around each control volume
  !! (currently 3, 9, or 12)
  !!
  SUBROUTINE lsq_stencil_create( ptr_patch, ptr_int_lsq, lsq_dim_c)
    !
    TYPE(t_patch), INTENT(INOUT) :: ptr_patch
    TYPE(t_lsq),   INTENT(INOUT) :: ptr_int_lsq
    INTEGER,       INTENT(IN)    :: lsq_dim_c    !< specifies size of the lsq stencil

    ! local
    INTEGER :: cnt                      ! counter

    CHARACTER(len=*), PARAMETER :: routine = modname//':lsq_stencil_create'

    !--------------------------------------------------------------------

    CALL message(routine, '')


    SELECT CASE(lsq_dim_c)
    !
    CASE (3) ! lsq_dim_c == 3
      IF (get_my_process_type() /= wave_process) THEN
        !
        ! 3-point stencil
        !
        CALL create_stencil_c3(p_patch     = ptr_patch,                 & !in
          &                    idx_c       = ptr_int_lsq%lsq_idx_c,     & !out
          &                    blk_c       = ptr_int_lsq%lsq_blk_c,     & !out
          &                    dim_stencil = ptr_int_lsq%lsq_dim_stencil) !out
      ELSE
        !
        ! 3-point stencil with special treatment of lateral boundaries
        !
        CALL create_stencil_c3_bnd(p_patch     = ptr_patch,                 & !in
          &                        idx_c       = ptr_int_lsq%lsq_idx_c,     & !out
          &                        blk_c       = ptr_int_lsq%lsq_blk_c,     & !out
          &                        dim_stencil = ptr_int_lsq%lsq_dim_stencil) !out
      ENDIF

    CASE (9) ! lsq_dim_c == 9
      !
      ! 9-point stencil
      !
      CALL create_stencil_c9(p_patch     = ptr_patch,                 & !in
        &                    idx_c       = ptr_int_lsq%lsq_idx_c,     & !out
        &                    blk_c       = ptr_int_lsq%lsq_blk_c,     & !out
        &                    dim_stencil = ptr_int_lsq%lsq_dim_stencil) !out

    CASE (12) ! lsq_dim_c == 12
      !
      ! 12-point stencil
      !
      CALL create_stencil_c12(p_patch     = ptr_patch,                 & !in
        &                     idx_c       = ptr_int_lsq%lsq_idx_c,     & !out
        &                     blk_c       = ptr_int_lsq%lsq_blk_c,     & !out
        &                     dim_stencil = ptr_int_lsq%lsq_dim_stencil) !out

    CASE DEFAULT
      WRITE(message_text,'(a,i2)') 'Could not create lsq stencil; invalid stencil size lsq_dim_c=', lsq_dim_c
      CALL finish(routine, message_text)
    END SELECT


    DO cnt = 1, lsq_dim_c
      CALL sync_idx(SYNC_C, SYNC_C, ptr_patch, ptr_int_lsq%lsq_idx_c(:,:,cnt), &
        &                                      ptr_int_lsq%lsq_blk_c(:,:,cnt))
    ENDDO
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_int_lsq%lsq_dim_stencil)

  END SUBROUTINE lsq_stencil_create


  !>
  !! For each cell create a 3-point stencil which consists of the 3 cells
  !! surrounding the control volume, i.e. the direct neighbors.
  !!
  !!
  SUBROUTINE create_stencil_c3 (p_patch, idx_c, blk_c, dim_stencil)

    TYPE(t_patch), INTENT(IN   ) :: p_patch

    INTEGER,       INTENT(INOUT) :: idx_c(:,:,:)      !< cell indizes
    INTEGER,       INTENT(INOUT) :: blk_c(:,:,:)      !< block indices
    INTEGER,       INTENT(INOUT) :: dim_stencil(:,:)  !< stencil size
    CHARACTER(len=*), PARAMETER :: routine = modname//':create_stencil_c3'

    CALL message(routine, 'create 3-point stencil')

    ! sanity check
    IF (ANY((/SIZE(idx_c,3),SIZE(blk_c,3)/) /= 3)) THEN
      CALL finish(routine, "Invalid size of output fields")
    ENDIF

    ! the cell and block indices are copied from p_patch%cells%neighbor_idx
    ! and p_patch%cells%neighbor_blk

!$OMP PARALLEL
    CALL copy(p_patch%cells%neighbor_idx(:,:,:), idx_c(:,:,:), lacc=.FALSE.)
    CALL copy(p_patch%cells%neighbor_blk(:,:,:), blk_c(:,:,:), lacc=.FALSE.)
    CALL copy(p_patch%cells%num_edges(:,:),      dim_stencil(:,:), lacc=.FALSE.)
!$OMP END PARALLEL

  END SUBROUTINE create_stencil_c3


  !>
  !! For each cell create a 3-point stencil which consists of the 3 cells
  !! surrounding the control volume, i.e. the direct neighbors.
  !!
  !! If a cell has less than 3 neighbors (i.e. if it is a boundary cell),
  !! the algorithm tries to create an asymmetric 3-point stencil, by searching
  !! for neighbors of the neighbors.
  !!
  !!
  SUBROUTINE create_stencil_c3_bnd (p_patch, idx_c, blk_c, dim_stencil)

    TYPE(t_patch), INTENT(IN   ) :: p_patch

    INTEGER,       INTENT(INOUT) :: idx_c(:,:,:)      !< cell indizes
    INTEGER,       INTENT(INOUT) :: blk_c(:,:,:)      !< block indices
    INTEGER,       INTENT(INOUT) :: dim_stencil(:,:)  !< stencil size

    INTEGER :: jb, jc                   !< loop index blocks/cells
    INTEGER :: jec, js                  !< loop index
    INTEGER :: cnt                      !< counter
    INTEGER :: i_rlstart, i_rlend       !< cell start and end row
    INTEGER :: i_startblk, i_endblk     !< start/end block
    INTEGER :: i_startidx, i_endidx     !< start/end index
    INTEGER :: ilc, ibc                 !< line and block index of neighbors
    INTEGER :: ilc_n, ibc_n             !< line and block index for neighbors of
                                        !  direct neighbors

    CHARACTER(len=*), PARAMETER :: routine = modname//':create_stencil_c3_bnd'

    CALL message(routine, 'create 3-point stencil with special boundary treatment')

    ! sanity check
    IF (ANY((/SIZE(idx_c,3),SIZE(blk_c,3)/) /= 3)) THEN
      CALL finish(routine, "Invalid size of output fields")
    ENDIF

    ! the cell and block indices are copied from p_patch%cells%neighbor_idx
    ! and p_patch%cells%neighbor_blk

!$OMP PARALLEL
    CALL copy(p_patch%cells%neighbor_idx(:,:,:), idx_c(:,:,:), lacc=.FALSE.)
    CALL copy(p_patch%cells%neighbor_blk(:,:,:), blk_c(:,:,:), lacc=.FALSE.)
    CALL copy(p_patch%cells%num_edges(:,:),      dim_stencil(:,:), lacc=.FALSE.)
!$OMP BARRIER

    ! special treatment of boundary cell row defined by refin_c_ctrl==1:
    ! construct asymmetric stencil for boundary cells
    !
    i_rlstart = 1
    i_rlend   = 1
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jc,jec,js,i_startidx,i_endidx,cnt,ilc,ibc,ilc_n,ibc_n) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,     &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF(.NOT. p_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        cnt = 0

        ! get valid direct neighbours
        DO jec = 1, 3
          ilc = p_patch%cells%neighbor_idx(jc,jb,jec)
          ibc = p_patch%cells%neighbor_blk(jc,jb,jec)

          IF (ilc > 0) THEN
            cnt = cnt + 1
            !
            idx_c(jc,jb,cnt) = ilc
            blk_c(jc,jb,cnt) = ibc
          ENDIF
        ENDDO


        ! check for incomplete stencil
        !
        IF (cnt==2) THEN
          ! search for valid neighbors of neighbors
          jsloop: DO js=1,2
            ilc = idx_c(jc,jb,js)
            ibc = blk_c(jc,jb,js)

            jecloop: DO jec=1,3
              ilc_n = p_patch%cells%neighbor_idx(ilc,ibc,jec)
              ibc_n = p_patch%cells%neighbor_blk(ilc,ibc,jec)

              ! if this is a valid cell and if it is not identical to
              ! the 'source' cell, add it to the stencil
              IF ( (ilc_n > 0) .AND. (ilc_n /= jc .OR. ibc_n /= jb)) THEN
                cnt = cnt + 1
                idx_c(jc,jb,cnt) = ilc_n
                blk_c(jc,jb,cnt) = ibc_n
                EXIT jsloop    ! we have found a valid candidate
              ENDIF
            ENDDO jecloop
          ENDDO jsloop

        ELSE IF (cnt==1) THEN

          ! search for valid neighbors of neighbors
          ! The current code assumes that we will find 2 more valid cells,
          ! which means that all neighbors of the neighbors must be valid cells.
          ! This is not necessarily the case for arbitrary grid resolutions!!
          !
          ilc = idx_c(jc,jb,1)
          ibc = blk_c(jc,jb,1)

          DO jec=1,3
            ilc_n = p_patch%cells%neighbor_idx(ilc,ibc,jec)
            ibc_n = p_patch%cells%neighbor_blk(ilc,ibc,jec)

            ! if this is a valid cell and if it is not identical to
            ! the 'source' cell, add it to the stencil
            IF ( (ilc_n > 0) .AND. (ilc_n /= jc .OR. ibc_n /= jb)) THEN
              cnt = cnt + 1
              idx_c(jc,jb,cnt) = ilc_n
              blk_c(jc,jb,cnt) = ibc_n
            ENDIF
          ENDDO

        ELSE if (cnt==0) THEN
          WRITE(message_text,'(a,i2,a)') 'insufficient stencil size. Unable to deal with an isolated cell.'
          CALL finish(routine, message_text)
        ENDIF


        ! sanity check
        IF (cnt < 3) THEN
          WRITE(message_text,'(a,i2,a)') 'Insufficient stencil size ', cnt, ' < 3'
          CALL finish(routine, message_text)
        ENDIF

        dim_stencil(jc,jb) = cnt

      ENDDO ! loop over cells
    ENDDO ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE create_stencil_c3_bnd



  !>
  !! For each cell create a 9-point stencil which consists of the
  !! 3 cells surrounding the control volume, i.e. the direct neighbors,
  !! and the neighbors of the neighbors.
  !!
  !!
  SUBROUTINE create_stencil_c9 (p_patch, idx_c, blk_c, dim_stencil)

    TYPE(t_patch), INTENT(IN   ) :: p_patch

    INTEGER,       INTENT(INOUT) :: idx_c(:,:,:)      !< cell indizes
    INTEGER,       INTENT(INOUT) :: blk_c(:,:,:)      !< block indices
    INTEGER,       INTENT(INOUT) :: dim_stencil(:,:)  !< stencil size

    INTEGER :: jb, jc                   !< loop index blocks/cells
    INTEGER :: jec, jj                  !< loop index
    INTEGER :: cnt                      !< counter
    INTEGER :: i_rlstart, i_rlend       !< cell start and end row
    INTEGER :: i_startblk, i_endblk     !< start/end block
    INTEGER :: i_startidx, i_endidx     !< start/end index
    INTEGER :: ilc, ibc                 !< line and block index of neighbors
    INTEGER :: ilc_n(3), ibc_n(3)       !< line and block index for neighbors of
                                        !  direct neighbors
    CHARACTER(len=*), PARAMETER :: routine = modname//':create_stencil_c9'

    CALL message(routine, 'create 9-point stencil')

    ! sanity check
    IF (ANY((/SIZE(idx_c,3),SIZE(blk_c,3)/) /= 9)) THEN
      CALL finish(routine, "Invalid size of output fields")
    ENDIF

    ! The start block depends on the width of the stencil
    i_rlstart = 2
    i_rlend   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jec,jj,i_startidx,i_endidx,cnt,ilc,ibc,ilc_n,ibc_n) &
!$OMP ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,     &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF(.NOT. p_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        cnt = 0

        DO jec = 1, 3

          ilc = p_patch%cells%neighbor_idx(jc,jb,jec)
          ibc = p_patch%cells%neighbor_blk(jc,jb,jec)

          ! direct neighbors
          cnt = cnt + 1
          idx_c(jc,jb,cnt) = ilc
          blk_c(jc,jb,cnt) = ibc

          ! neighbors of direct neighbors
          DO jj = 1,3

            ilc_n(jj) = p_patch%cells%neighbor_idx(ilc,ibc,jj)
            ibc_n(jj) = p_patch%cells%neighbor_blk(ilc,ibc,jj)

            IF (ilc_n(jj) /= jc .OR. ibc_n(jj) /= jb) THEN
              cnt = cnt + 1
              idx_c(jc,jb,cnt) = ilc_n(jj)
              blk_c(jc,jb,cnt) = ibc_n(jj)
            ENDIF
          ENDDO

        ENDDO ! jec loop

        dim_stencil(jc,jb) = cnt

      ENDDO ! loop over cells
    ENDDO ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE create_stencil_c9


  !>
  !! For each cell create a 12-point stencil which consists of all
  !! vertex neighbors.
  !! It is similar to the 9-point stencil, but more isotropic.
  !!
  !! Note: At pentagon points the stencil size reduces to 11.
  !!
  SUBROUTINE create_stencil_c12 (p_patch, idx_c, blk_c, dim_stencil)

    TYPE(t_patch), INTENT(IN   ) :: p_patch

    INTEGER,       INTENT(INOUT) :: idx_c(:,:,:)      !< cell indizes
    INTEGER,       INTENT(INOUT) :: blk_c(:,:,:)      !< block indices
    INTEGER,       INTENT(INOUT) :: dim_stencil(:,:)  !< stencil size

    INTEGER :: jb, jc                   !< loop index blocks/cells
    INTEGER :: jec, jj, jtri            !< loop index
    INTEGER :: cnt                      !< counter
    INTEGER :: i_rlstart, i_rlend       !< cell start and end row
    INTEGER :: i_startblk, i_endblk     !< start/end block
    INTEGER :: i_startidx, i_endidx     !< start/end index
    INTEGER :: ilc_n(3), ibc_n(3)       !< line and block index for neighbors of
                                        !  direct neighbors
    INTEGER :: ilv(3), ibv(3)           !< vertex line and block indices
    INTEGER :: ilc_v(3,6), ibc_v(3,6)   !< cell line and block indices
                                        !  around each of the three vertices
    CHARACTER(len=*), PARAMETER :: routine = modname//':create_stencil_c12'

    CALL message(routine, 'create 12-point stencil')

    ! sanity check
    IF (ANY((/SIZE(idx_c,3),SIZE(blk_c,3)/) /= 12)) THEN
      CALL finish(routine, "Invalid size of output fields")
    ENDIF

    ! The start block depends on the width of the stencil
    i_rlstart = 2
    i_rlend   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jec,jj,jtri,i_startidx,i_endidx,cnt,ilv,ibv, &
!$OMP            ilc_v,ibc_v,ilc_n,ibc_n) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,     &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF(.NOT. p_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        cnt = 0

        ! get line and block indices of cell vertices
        ilv(1:3) = p_patch%cells%vertex_idx(jc,jb,1:3)
        ibv(1:3) = p_patch%cells%vertex_blk(jc,jb,1:3)

        ! for each vertex: get all cells which share this vertex
        DO jj = 1,3
          ilc_v(jj,:)=p_patch%verts%cell_idx(ilv(jj),ibv(jj),:)
          ibc_v(jj,:)=p_patch%verts%cell_blk(ilv(jj),ibv(jj),:)
        ENDDO

        !
        ! 1. add the 3 direct neighbors to the stencil
        !
        DO jec = 1,3
          ! get line and block indices of direct neighbors
          ilc_n(jec) = p_patch%cells%neighbor_idx(jc,jb,jec)
          ibc_n(jec) = p_patch%cells%neighbor_blk(jc,jb,jec)

          cnt = cnt + 1

          idx_c(jc,jb,cnt) = ilc_n(jec)
          blk_c(jc,jb,cnt) = ibc_n(jec)
        ENDDO

        !
        ! 2. loop over the vertices and add all cells
        !    that are no direct neighbors and not our CV.
        !
        DO jj = 1,3   ! loop over vertices
          DO jtri=1,6 ! loop over cells around each vertex

            IF (.NOT.( (ilc_v(jj,jtri) == ilc_n(1) .AND. ibc_v(jj,jtri) == ibc_n(1))  &
              &  .OR.  (ilc_v(jj,jtri) == ilc_n(2) .AND. ibc_v(jj,jtri) == ibc_n(2))  &
              &  .OR.  (ilc_v(jj,jtri) == ilc_n(3) .AND. ibc_v(jj,jtri) == ibc_n(3))  &
              &  .OR.  (ilc_v(jj,jtri) == jc       .AND. ibc_v(jj,jtri) == jb)        &
              &  .OR.  (ilc_v(jj,jtri) == 0        .AND. ibc_v(jj,jtri) == 0 ) ) ) THEN

              cnt = cnt + 1

              idx_c(jc,jb,cnt) = ilc_v(jj,jtri)
              blk_c(jc,jb,cnt) = ibc_v(jj,jtri)
            ENDIF
          ENDDO
        ENDDO

        dim_stencil(jc,jb) = cnt

      ENDDO ! loop over cells

    ENDDO ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE create_stencil_c12


!-------------------------------------------------------------------------
!
!
!>
!! This routine bifurcates into lsq_compute_coeff_cell based on geometry type
!!
!! Anurag Dipankar, MPIM (2013-04)
!!
SUBROUTINE lsq_compute_coeff_cell( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
  &                                lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )
!

!
TYPE(t_patch), INTENT(INOUT) ::  ptr_patch

TYPE(t_lsq), TARGET, INTENT(INOUT) ::  ptr_int_lsq

LOGICAL, INTENT(IN) ::   &  ! flag determining whether the least
  &  llsq_rec_consv         ! squares reconstruction should be conservative

INTEGER, INTENT(IN)  ::  &  ! parameter determining the size of the lsq stencil
  &  lsq_dim_c

INTEGER, INTENT(IN)  ::  &  ! parameter determining the dimension of the solution
  &  lsq_dim_unk

INTEGER, INTENT(IN)  ::  &  ! least squares weighting exponent
  &  lsq_wgt_exp

 !Local variable
 CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs_lsq_bln:lsq_compute_coeff_cell'

    !
    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)

     CALL lsq_compute_coeff_cell_torus( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
                                        lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )

    CASE (sphere_geometry)

     CALL lsq_compute_coeff_cell_sphere( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
                                         lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )

    CASE DEFAULT

      CALL finish(method_name, "Undefined geometry type")

    END SELECT


END SUBROUTINE lsq_compute_coeff_cell


!-------------------------------------------------------------------------
!
!! AD: This routine has been just renamed with affix "_sphere"
!!
!! This routine computes the coefficients needed for a weighted least-squares.
!!
!! This routine computes the coefficients needed for a weighted least-squares
!! reconstruction at cell centers. Optionally, the reconstruction can be
!! enforced to be conservative in the sense that, when integrated over the
!! control volume, it recovers the area average stored at the mass point.
!! Works for triangular and hexagonal control volumes.
!!
SUBROUTINE lsq_compute_coeff_cell_sphere( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
  &                                       lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )
!

!
TYPE(t_patch), INTENT(INOUT) ::  ptr_patch

TYPE(t_lsq), TARGET, INTENT(INOUT) ::  ptr_int_lsq

LOGICAL, INTENT(IN) ::   &  ! flag determining whether the least
  &  llsq_rec_consv         ! squares reconstruction should be conservative

INTEGER, INTENT(IN)  ::  &  ! parameter determining the size of the lsq stencil
  &  lsq_dim_c

INTEGER, INTENT(IN)  ::  &  ! parameter determining the dimension of the solution
  &  lsq_dim_unk

INTEGER, INTENT(IN)  ::  &  ! least squares weighting exponent
  &  lsq_wgt_exp

REAL(wp), DIMENSION(lsq_dim_c,2) ::  &      ! geographical coordinates of all cell centers
  & xytemp_c                                ! in the stencil

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type,2)   ::  &  ! geogr. coordinates of vertices of the
  & xytemp_v                                ! control volume

REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c,lsq_dim_c,2) ::  &
  & z_dist_g                                ! for each cell:
                                            ! distance vectors from control volume cell center 
                                            ! to cell centers of all cells in the stencil

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type,2)   ::  &  ! lat/lon distance vector edge midpoint -> cvertex
  & distxy_v

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) :: dely, delx ! difference in latitude and longitude between
                                               ! vertices

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk) ::  & ! lsq matrix
  & z_lsq_mat_c

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk)   ::  &  ! Q matrix of QR-factorization
  & z_qmat

REAL(wp), DIMENSION(nproma,lsq_dim_unk,lsq_dim_unk) ::  &  ! R matrix of QR-factorization
  & z_rmat

REAL(wp) :: z_rcarea                   ! reciprocal of cell area

REAL(wp) :: xloc, yloc                 ! geographical coordinates of
                                       ! point under consideration

REAL(wp) :: z_norm                     ! vector length (distance between control volume
                                       ! center and cell centers in the stencil on tangent
                                       ! plane) (also used for normalization)

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) ::  & ! integrand for each edge
  & fx, fy, fxx, fyy, fxy,           & ! for analytical calculation of moments
  & fxxx, fyyy, fxxy, fxyy

INTEGER, POINTER  ::       &           ! pointer to stencil size (cell dependent)
  & ptr_ncells(:,:)

INTEGER, DIMENSION(lsq_dim_c) ::  &    ! line and block indices of cells in the stencil
  & ilc_s, ibc_s
INTEGER, DIMENSION(ptr_patch%geometry_info%cell_type) :: jlv, jbv      ! line and block indices of vertex
INTEGER :: cnt                         ! counter
INTEGER :: jrow                        ! matrix row-identifier
INTEGER :: nel                         ! number of matrix elements
INTEGER :: nblks_c
INTEGER :: pid                         ! patch ID
INTEGER :: jb                          ! index of current block
INTEGER :: jc                          ! index of current cell
INTEGER :: js                          ! index of current control volume in the stencil
INTEGER :: ju                          ! loop index for column of lsq matrix
INTEGER :: jec                         ! loop index for cell's edge
INTEGER :: i_startblk                  ! start block
INTEGER :: i_startidx                  ! start index
INTEGER :: i_endidx                    ! end index
INTEGER :: i_rcstartlev                ! refinement control start level
INTEGER :: ist, icheck                 ! status
INTEGER :: nverts
INTEGER :: jecp
INTEGER :: jja, jjb, jjk               ! loop indices for Moore-Penrose inverse

REAL(wp) ::   &                        ! singular values of lsq design matrix A
  &  zs(lsq_dim_unk,nproma)            ! min(lsq_dim_c,lsq_dim_unk)

REAL(wp) ::   &                        ! U matrix of SVD. Columns of U are the left
  &  zu  (lsq_dim_c,lsq_dim_c,nproma)  ! singular vectors of A

REAL(wp) ::   &                        ! TRANSPOSE of V matrix of SVD. Columns of V are
  &  zv_t(lsq_dim_unk,lsq_dim_unk,nproma) ! the right singular vectors of A.


INTEGER, PARAMETER  :: &     ! size of work array for SVD lapack routine
  &  lwork=10000
REAL(wp) ::   &              ! work array for SVD lapack routine
  &  zwork(lwork)
INTEGER  ::   &              ! work array for SVD lapack routine
  & ziwork(8*min(lsq_dim_c,lsq_dim_unk))


!DR for DEBUG purposes
! #ifdef DEBUG_COEFF LL it's used in openmp directives,
REAL(wp) :: za_debug(nproma,lsq_dim_c,lsq_dim_unk)
! #endif
!--------------------------------------------------------------------


  CALL message('mo_interpolation:lsq_compute_coeff_cell_sphere', '')

  i_rcstartlev = 2

  ! get patch id
  pid = ptr_patch%id

  ! stencil size
  ptr_ncells => ptr_int_lsq%lsq_dim_stencil(:,:)

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)


!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jc,js,jec,i_startidx,i_endidx,jlv,jbv,ilc_s,ibc_s, &
!!$OMP            xloc,yloc,xytemp_c,xytemp_v,z_norm,distxy_v,z_rcarea, &
!!$OMP            delx,dely,fx,fy,fxx,fyy,fxy,jecp,nverts) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,     &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each cell, calculate weights, moments, matrix coefficients
    ! and QR decomposition of normal equation matrix
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      IF (ptr_patch%geometry_info%cell_type == 3 )THEN
        nverts  = 3
      ELSE
        nverts = ptr_patch%cells%num_edges(jc,jb)
      ENDIF
    !
    ! Gather some information about the control volume and
    ! all the cells in the stencil.
    !
    ! get line and block indices of edge vertices
      jlv(1:nverts) = ptr_patch%cells%vertex_idx(jc,jb,1:nverts)
      jbv(1:nverts) = ptr_patch%cells%vertex_blk(jc,jb,1:nverts)

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    !
    ! 1. Get geographical coordinates of the control volume center
    !    and all cell centers in the stencil. In addition, get geographical
    !    coordinates of the vertices of the control volume.
    !
    ! get geographical coordinates of control volume center
      xloc = ptr_patch%cells%center(jc,jb)%lon
      yloc = ptr_patch%cells%center(jc,jb)%lat

    ! get geogr. coordinates of all cell centers in the stencil
      DO js = 1, ptr_ncells(jc,jb)
        xytemp_c(js,1) = ptr_patch%cells%center(ilc_s(js),ibc_s(js))%lon
        xytemp_c(js,2) = ptr_patch%cells%center(ilc_s(js),ibc_s(js))%lat
      ENDDO

    ! get geogr. coordinates of edge-vertices (only for control volume)
      DO jec = 1, nverts
        xytemp_v(jec,1) = ptr_patch%verts%vertex(jlv(jec),jbv(jec))%lon
        xytemp_v(jec,2) = ptr_patch%verts%vertex(jlv(jec),jbv(jec))%lat
      ENDDO


    !
    ! 2. Now, all points (centers, vertices) are projected onto a tangent
    !    plane having its center at the cell center of the control volume.
    !    At the same time distance vectors (r_x,r_y) are calculated between
    !    the control volume center and all other mass-points in the stencil.
    !    From these the lsq weights can be deduced.

    !
    ! a: Project centers and calculate distance vectors between control volume
    !    center and all cell centers in the stencil in geographical coordinates.
    !    Since the control volume center has the local coordinates (x,y)=(0,0),
    !    the coordinates of the other cell centers equal the distance
    !    vectors.
!$NEC novector
      DO js = 1, ptr_ncells(jc,jb)
        CALL gnomonic_proj( xloc, yloc, xytemp_c(js,1), xytemp_c(js,2),  &! in
         &                  z_dist_g(jc,jb,js,1), z_dist_g(jc,jb,js,2) )  ! out

      ENDDO
      ! multiply with earth radius and store
      z_dist_g(jc,jb,1:ptr_ncells(jc,jb),:) = grid_sphere_radius * &
        & z_dist_g(jc,jb,1:ptr_ncells(jc,jb),:)



    !
    ! b: compute normalized weights for weighted least squares system
    !    The closest cell circumcenter is assigned weight of w=1.
    !
      DO js = 1, ptr_ncells(jc,jb)
        z_norm = SQRT(DOT_PRODUCT(z_dist_g(jc,jb,js,1:2),z_dist_g(jc,jb,js,1:2)))

        !
        ! weights for weighted least squares system
        !
        ptr_int_lsq%lsq_weights_c(jc,js,jb)= 1._wp/(z_norm**lsq_wgt_exp)

      ENDDO
      !
      ! Normalization
      !
      ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)=                &
        &     ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)         &
        &   / MAXVAL(ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb))



    ! 3. the following part (including calculation of moments) will only
    !    be called, if a conservative least squares reconstruction
    !    is chosen. Otherwise all moments will be equal to zero and the
    !    reconstruction simplifies to the standard non-conservative
    !    reconstruction.
      IF (llsq_rec_consv) THEN
      !
      ! a: Project control volume vertices and calculate distance vectors
      !    between cell center and vertices
      !
!$NEC novector
        DO jec=1,nverts
          CALL gnomonic_proj( xloc, yloc, xytemp_v(jec,1), xytemp_v(jec,2),  &
           &                 distxy_v(jec,1), distxy_v(jec,2) )

        ENDDO
      ! multiply with earth radius
        distxy_v(1:nverts,1:2) = grid_sphere_radius * distxy_v(1:nverts,1:2)


      !
      ! b: calculate moments for given cell
      !    (calculated analytically; see Lauritzen CSLAM 09 for first 5 moments)
      !
      ! !DR: Those moments have been re-rechecked, using an alternative,
      !      quadrature-based formulation. Results have been identical up to
      !      roundoff-errors. Similarly the hat-moments have been checked. The
      !      inconsistency caused by the different projections involved do not
      !      seem to negatively effect the results.

      ! Storage docu for x^ny^m:
      ! lsq_moments(:,:,1) : x^1y^0
      ! lsq_moments(:,:,2) : x^0y^1
      ! lsq_moments(:,:,3) : x^2y^0
      ! lsq_moments(:,:,4) : x^0y^2
      ! lsq_moments(:,:,5) : x^1y^1
      ! lsq_moments(:,:,6) : x^3y^0
      ! lsq_moments(:,:,7) : x^0y^3
      ! lsq_moments(:,:,8) : x^2y^1
      ! lsq_moments(:,:,9) : x^1y^2
      !

        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! note that the distance vector distxy_v between each vertex and
          ! the center of the tangent plane are identical to the coordinates of
          ! each vertex on the tangent plane. Thus the distances between the
          ! vertices in x and y direction can be derived as follows:
          !
          ! longitudinal-distance between vertices on tangent plane
          delx(jec) = distxy_v(jecp,1) - distxy_v(jec,1)

          ! latitudinal-distance between vertices on tangent plane
          dely(jec) = distxy_v(jecp,2) - distxy_v(jec,2)


          !
          ! analytic moment calculation
          !
          ! 0: control volume area (reciprocal value)
          fx(jec) = distxy_v(jecp,1) + distxy_v(jec,1)

        ENDDO

        z_rcarea = 2._wp/DOT_PRODUCT(fx(1:nverts),dely(1:nverts))


        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! I. x^1y^0
          fx(jec) = distxy_v(jecp,1)**2              &
          &       + distxy_v(jecp,1)*distxy_v(jec,1) &
          &       + distxy_v(jec ,1)**2

          ! II. x^0y^1
          fy(jec) = distxy_v(jecp,2)**2              &
          &       + distxy_v(jecp,2)*distxy_v(jec,2) &
          &       + distxy_v(jec,2)**2

          IF ( lsq_dim_unk > 2 ) THEN

            ! III. x^2y^0
            fxx(jec) = (distxy_v(jecp,1)    + distxy_v(jec,1)       ) &
            &        * (distxy_v(jecp,1)**2 + distxy_v(jec,1)**2)

            ! IV. x^0y^2
            fyy(jec) = (distxy_v(jecp,2)    + distxy_v(jec,2)       ) &
            &        * (distxy_v(jecp,2)**2 + distxy_v(jec,2)**2)

            ! V. x^1y^1
            fxy(jec) = distxy_v(jecp,2) * (3._wp*distxy_v(jecp,1)**2                         &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + distxy_v(jec,1)**2 )     &
            &        + distxy_v(jec,2) * ( distxy_v(jecp,1)**2                               &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + 3._wp*distxy_v(jec,1)**2 )

          ENDIF ! lsq_dim_unk > 2

          IF ( lsq_dim_unk > 5 ) THEN

            ! VI.  x^3y^0
            fxxx(jec) = 5._wp*distxy_v(jec,1)**4                   &
              &       + 10._wp*distxy_v(jec,1)**3 * delx(jec)      &
              &       + 10._wp*distxy_v(jec,1)**2 * delx(jec)**2   &
              &       + 5._wp *distxy_v(jec,1)    * delx(jec)**3   &
              &       + delx(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fxxx(jec) = distxy_v(jecp,1)**4                       &
            !DR    &       + distxy_v(jec,1) * distxy_v(jecp,1)**3     &
            !DR    &       + distxy_v(jec,1)**2 * distxy_v(jecp,1)**2  &
            !DR    &       + distxy_v(jec,1)**3 * distxy_v(jecp,1)     &
            !DR    &       + distxy_v(jec,1)**4


            ! VII. x^0y^3
            fyyy(jec) = 5._wp*distxy_v(jec,2)**4                   &
              &       + 10._wp*distxy_v(jec,2)**3 * dely(jec)      &
              &       + 10._wp*distxy_v(jec,2)**2 * dely(jec)**2   &
              &       + 5._wp *distxy_v(jec,2)    * dely(jec)**3   &
              &       + dely(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fyyy(jec) = distxy_v(jecp,2)**4                       &
            !DR    &       + distxy_v(jec,2)    * distxy_v(jecp,2)**3  &
            !DR    &       + distxy_v(jec,2)**2 * distxy_v(jecp,2)**2  &
            !DR    &       + distxy_v(jec,2)**3 * distxy_v(jecp,2)     &
            !DR    &       + distxy_v(jec,2)**4


            ! VIII. x^2y^1
            fxxy(jec) = 4._wp*distxy_v(jecp,1)**3 *distxy_v(jecp,2)                    &
              &       + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jecp,2)   &
              &       + 2._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jecp,2) &
              &       +       distxy_v(jec,1)**3 * distxy_v(jecp,2)                    &
              &       +       distxy_v(jecp,1)**3 * distxy_v(jec,2)                    &
              &       + 2._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jec,2)    &
              &       + 3._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jec,2)  &
              &       + 4._wp*distxy_v(jec,1)**3 * distxy_v(jec,2)

            ! IX. x^1y^2
            fxyy(jec) = 6._wp*distxy_v(jecp,1)**2 * distxy_v(jecp,2)**2                     &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jecp,2)**2              &
              &   +       distxy_v(jec,1)**2 * distxy_v(jecp,2)**2                          &
              &   + 3._wp*distxy_v(jecp,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)            &
              &   + 4._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)*distxy_v(jecp,2) &
              &   + 3._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)             &
              &   +       distxy_v(jecp,1)**2 * distxy_v(jec,2)**2                          &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)**2               &
              &   + 6._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)**2

          ENDIF ! lsq_dim_unk > 5


        ENDDO ! loop over nverts

        ptr_int_lsq%lsq_moments(jc,jb,1)= z_rcarea/6._wp*DOT_PRODUCT(fx(1:nverts),dely(1:nverts))
        ptr_int_lsq%lsq_moments(jc,jb,2)=-z_rcarea/6._wp*DOT_PRODUCT(fy(1:nverts),delx(1:nverts))

        IF ( lsq_dim_unk > 2 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,3) =  z_rcarea/12._wp * &
            & DOT_PRODUCT(fxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,4) = -z_rcarea/12._wp * &
            & DOT_PRODUCT(fyy(1:nverts),delx(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,5) =  z_rcarea/24._wp * &
            & DOT_PRODUCT(fxy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 2


        IF ( lsq_dim_unk > 5 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,6) =  z_rcarea/20._wp * &
            & DOT_PRODUCT(fxxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,7) = -z_rcarea/20._wp * &
            & DOT_PRODUCT(fyyy(1:nverts),delx(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,8) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxxy(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,9) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxyy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 5

      END IF  ! llsq_rec_consv

    END DO  ! loop over cells

  END DO  ! loop over blocks
!!$OMP END DO NOWAIT
!! For unknown reasons, closing the parallel section here is needed to get the above
!! loop parallelized.
!!$OMP END PARALLEL

  DO jb = 1, lsq_dim_c
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_weights_c(:,jb,:))
  ENDDO
  DO jb = 1, lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments(:,:,jb))
  ENDDO


!$OMP PARALLEL PRIVATE(jb,jc,js,ju,jja,jjb,jjk,i_startidx,i_endidx,ilc_s,ibc_s, &
!$OMP            z_lsq_mat_c,zs,zu,zv_t,zwork,ziwork,ist,icheck,za_debug, &
!$OMP            z_qmat,z_rmat,cnt,jrow,nel)
!$OMP DO ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! 4. for each cell, calculate LSQ design matrix A
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
        ! Take care that z_lsq_mat_c isn't singular
        z_lsq_mat_c(jc,:,:) = 0.0_wp
        DO js = 1, MIN(lsq_dim_unk, lsq_dim_c)
          z_lsq_mat_c(jc,js,js) = 1.0_wp
        ENDDO
        CYCLE
      ENDIF

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    ! Calculate full moments lsq_moments_hat(ilc_s(js),ibc_s(js),ju)
    !
    ! Storage docu for x^ny^m:
    ! lsq_moments_hat(:,:,:,1) : \hat{x^1y^0}
    ! lsq_moments_hat(:,:,:,2) : \hat{x^0y^1}
    ! lsq_moments_hat(:,:,:,3) : \hat{x^2y^0}
    ! lsq_moments_hat(:,:,:,4) : \hat{x^0y^2}
    ! lsq_moments_hat(:,:,:,5) : \hat{x^1y^1}
    ! lsq_moments_hat(:,:,:,6) : \hat{x^3y^0}
    ! lsq_moments_hat(:,:,:,7) : \hat{x^0y^3}
    ! lsq_moments_hat(:,:,:,8) : \hat{x^2y^1}
    ! lsq_moments_hat(:,:,:,9) : \hat{x^1y^2}
    !
      DO js = 1, ptr_ncells(jc,jb)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,1) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1) + z_dist_g(jc,jb,js,1)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,2) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2) + z_dist_g(jc,jb,js,2)

        IF (lsq_dim_unk > 2) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,3) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
           &    + z_dist_g(jc,jb,js,1)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,4) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2) &
           &    + z_dist_g(jc,jb,js,2)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,5) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)                              &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)        &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)        &
           &    + z_dist_g(jc,jb,js,1) * z_dist_g(jc,jb,js,2)
        ENDIF

        IF ( lsq_dim_unk > 5 ) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,6) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),6)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,1) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)                       &
            &    * z_dist_g(jc,jb,js,1)**2                                                    &
            &    + z_dist_g(jc,jb,js,1)**3

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,7) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),7)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,2) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)                       &
            &    * z_dist_g(jc,jb,js,2)**2                                                    &
            &    + z_dist_g(jc,jb,js,2)**3

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,8) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),8)                             &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,2)       &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,1) &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
            &    * z_dist_g(jc,jb,js,2)                                                       &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)**2    &
            &    + z_dist_g(jc,jb,js,1)**2 * z_dist_g(jc,jb,js,2)


          ptr_int_lsq%lsq_moments_hat(jc,jb,js,9) =                                          &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),9)                            &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,1)      &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,2)&
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2)&
            &    * z_dist_g(jc,jb,js,1)                                                      &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)**2   &
            &    + z_dist_g(jc,jb,js,2)**2 * z_dist_g(jc,jb,js,1)
        ENDIF
      ENDDO


      ! loop over rows of lsq design matrix (all cells in the stencil)
      DO js = 1, ptr_ncells(jc,jb)
        ! loop over columns of lsq design matrix (number of unknowns)
        DO ju = 1,lsq_dim_unk

          z_lsq_mat_c(jc,js,ju) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                         &
           &    * (ptr_int_lsq%lsq_moments_hat(jc,jb,js,ju)-ptr_int_lsq%lsq_moments(jc,jb,ju))

        END DO
       END DO
       IF(ptr_ncells(jc,jb) < lsq_dim_c) THEN
         z_lsq_mat_c(jc,lsq_dim_c,:) = 0.0_wp
       ENDIF

    ENDDO   ! loop over cells



    !
    ! compute QR decomposition and Singular Value Decomposition (SVD)
    ! of least squares design matrix A. For the time being both methods are
    ! retained.

    !
    ! 5a. QR-factorization of design matrix A
    !
    IF (.NOT. advection_config(pid)%llsq_svd) THEN
!CDIR NOIEXPAND
    CALL qrdec(lsq_dim_c, lsq_dim_unk, i_startidx, & ! in
     &         i_endidx, z_lsq_mat_c,              & ! in
     &         z_qmat, z_rmat)                       ! out


    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! 7. Save transposed Q-Matrix
      ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,1:lsq_dim_c,jb)  =  &
       &      TRANSPOSE(z_qmat(jc,1:lsq_dim_c,1:lsq_dim_unk))

      ! 8. Save R-Matrix
      !
      ! a. Save reciprocal values of the diagonal elements
      !
      DO ju = 1,lsq_dim_unk
        ptr_int_lsq%lsq_rmat_rdiag_c(jc,ju,jb) = 1._wp/z_rmat(jc,ju,ju)
      ENDDO

      !
      ! b. Save upper triangular elements without diagonal elements in a 1D-array
      !    (starting from the bottom right)
      !
      cnt = 1
      DO jrow = lsq_dim_unk-1,1,-1
        ! number of elements to store
        nel = lsq_dim_unk - jrow
        ptr_int_lsq%lsq_rmat_utri_c(jc,cnt:cnt+nel-1,jb) = z_rmat(jc,jrow,jrow+1:lsq_dim_unk)
        cnt = cnt + nel
      ENDDO


      ! Multiply ith column of the transposed Q-matrix (corresponds to the
      ! different members of the stencil) with the ith weight. This avoids
      ! multiplication of the RHS of the LSQ-System with this weight during
      ! runtime.
      DO js = 1,lsq_dim_c
        ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  =                  &
          &                ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  &
          &              * ptr_int_lsq%lsq_weights_c(jc,js,jb)
      ENDDO

    END DO  ! loop over cells

    ELSE   ! llsq_svd=.TRUE.

    !
    ! 5b. Singular value decomposition of lsq design matrix A
    ! !!! does not vectorize, unfortunately !!!
    ist = 0
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! A = U * SIGMA * transpose(V)
      !
      ! z_lsq_mat_c : M x N least squares design matrix A            (IN)
      ! zu          : M x M orthogonal matrix U                      (OUT)
      ! zv_t        : N x N orthogonal matrix V                      (OUT)
      ! zs          : min(M,N) Singular values of A                  (OUT)
      ! zwork       : workspace(1,LWORK)                             (OUT)
      ! lwork       : 3*min(M,N)                                     (IN)
      !              + max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N))  (IN)
      ! iwork       : workspace(8*min(M,N))                          (IN)

      CALL DGESDD('A',                 & !in
        &         lsq_dim_c,           & !in
        &         lsq_dim_unk,         & !in
        &         z_lsq_mat_c(jc,:,:), & !inout Note: destroyed on output
        &         lsq_dim_c,           & !in
        &         zs(:,jc),            & !out
        &         zu(:,:,jc),          & !out
        &         lsq_dim_c,           & !in
        &         zv_t(:,:,jc),        & !out
        &         lsq_dim_unk,         & !in
        &         zwork,               & !out
        &         lwork,               & !in
        &         ziwork,              & !inout
        &         icheck               ) !out
      ist = ist + icheck
    ENDDO
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:lsq_compute_coeff_cell_sphere',   &
        &             'singular value decomposition failed')
    ENDIF

    ! compute Moore-Penrose inverse
    ! INVERSE(A):: V * INVERSE(SIGMA) * TRANSPOSE(U) and store
    ! note that the ith column is multiplied with the ith weight
    ! in order to avoid the weighting of the r.h.s. during runtime.
    DO jja = 1, lsq_dim_unk
      DO jjb = 1, lsq_dim_c
        DO jjk = 1, lsq_dim_unk
          DO jc = i_startidx, i_endidx
            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
            ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb) =            &
              &  ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb)         &
              &  + zv_t(jjk,jja,jc) /zs(jjk,jc) * zu(jjb,jjk,jc)  &
              &  * ptr_int_lsq%lsq_weights_c(jc,jjb,jb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!!!!!DEBUG !!!
#ifdef DEBUG_COEFF
    za_debug(:,:,:)  = 0._wp
    ! re-COMPUTE A:: U * SIGMA * TRANSPOSE(V)  !!! Funktioniert
    DO jja = 1, lsq_dim_c
      DO jjb = 1, lsq_dim_unk
        DO jjk = 1, lsq_dim_unk  !lsq_dim_c
          DO jc = i_startidx, i_endidx
            za_debug(jc,jja,jjb) = za_debug(jc,jja,jjb)  &
              &  + zu(jja,jjk,jc) * zs(jjk,jc) * zv_t(jjk,jjb,jc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    write(0,*) "za_debug(10,:,:)- z_lsq_mat_c(10,:,:)",za_debug(10,:,:)- z_lsq_mat_c(10,:,:)
    write(0,*) "za_debug(55,:,:)- z_lsq_mat_c(55,:,:)",za_debug(55,:,:)- z_lsq_mat_c(55,:,:)
    write(0,*) "zs(:,10) ",zs(:,10)
    write(0,*) "zs(:,55) ",zs(:,55)
    write(0,*) "ptr_int_lsq%lsq_weights_c(10,:,jb)",ptr_int_lsq%lsq_weights_c(10,:,jb)
    write(0,*) "ptr_int_lsq%lsq_weights_c(55,:,jb)",ptr_int_lsq%lsq_weights_c(55,:,jb)
#endif
!!!! END DEBUG !!!

    ENDIF  ! llsq_svd

  END DO  ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  DO ju = 1, lsq_dim_unk
    DO jc = 1, lsq_dim_c
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments_hat(:,:,jc,ju))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_pseudoinv(:,ju,jc,:))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_qtmat_c(:,ju,jc,:))
    ENDDO
  ENDDO

  DO jc = 1, UBOUND(ptr_int_lsq%lsq_rmat_utri_c, 2)
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_utri_c(:,jc,:))
  ENDDO
  DO ju = 1,lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_rdiag_c(:,ju,:))
  ENDDO

END SUBROUTINE lsq_compute_coeff_cell_sphere


!-------------------------------------------------------------------------
!
!! This is same routine as lsq_compute_coeff_cell_sphere just modified for
!! flat geometry
!!
SUBROUTINE lsq_compute_coeff_cell_torus( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
  &                                      lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )
!

!
TYPE(t_patch), INTENT(INOUT) ::  ptr_patch

TYPE(t_lsq), TARGET, INTENT(INOUT) ::  ptr_int_lsq

LOGICAL, INTENT(IN) ::   &  ! flag determining whether the least
  &  llsq_rec_consv         ! squares reconstruction should be conservative

INTEGER, INTENT(IN)  ::  &  ! parameter determining the size of the lsq stencil
  &  lsq_dim_c

INTEGER, INTENT(IN)  ::  &  ! parameter determining the dimension of the solution
  &  lsq_dim_unk

INTEGER, INTENT(IN)  ::  &  ! least squares weighting exponent
  &  lsq_wgt_exp


!CC of points in the stencil
TYPE(t_cartesian_coordinates) :: cc_cv, cc_cell(lsq_dim_c), cc_vert(ptr_patch%geometry_info%cell_type)

REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c,lsq_dim_c,2) ::  &
  & z_dist_g                                ! for each cell:
                                            ! distance vectors from control volume cell center 
                                            ! to cell centers of all cells in the stencil

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type,2)   ::  &  ! lat/lon distance vector edge midpoint -> cvertex
  & distxy_v

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) :: dely, delx ! difference in latitude and longitude between
                                               ! vertices

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk) ::  & ! lsq matrix
  & z_lsq_mat_c

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk)   ::  &  ! Q matrix of QR-factorization
  & z_qmat

REAL(wp), DIMENSION(nproma,lsq_dim_unk,lsq_dim_unk) ::  &  ! R matrix of QR-factorization
  & z_rmat

REAL(wp) :: z_rcarea                   ! reciprocal of cell area

REAL(wp) :: z_norm                     ! vector length (distance between control volume
                                       ! center and cell centers in the stencil on tangent
                                       ! plane) (also used for normalization)

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) ::  & ! integrand for each edge
  & fx, fy, fxx, fyy, fxy,           & ! for analytical calculation of moments
  & fxxx, fyyy, fxxy, fxyy

INTEGER, POINTER  ::       &           ! pointer to stencil size (cell dependent)
  & ptr_ncells(:,:)

INTEGER, DIMENSION(lsq_dim_c) ::  &    ! line and block indices of cells in the stencil
  & ilc_s, ibc_s
INTEGER, DIMENSION(ptr_patch%geometry_info%cell_type) :: jlv, jbv      ! line and block indices of vertex
INTEGER :: cnt                         ! counter
INTEGER :: jrow                        ! matrix row-identifier
INTEGER :: nel                         ! number of matrix elements
INTEGER :: nblks_c
INTEGER :: pid                         ! patch ID
INTEGER :: jb                          ! index of current block
INTEGER :: jc                          ! index of current cell
INTEGER :: js                          ! index of current control volume in the stencil
INTEGER :: ju                          ! loop index for column of lsq matrix
INTEGER :: jec                         ! loop index for cell's edge
INTEGER :: i_startblk                  ! start block
INTEGER :: i_startidx                  ! start index
INTEGER :: i_endidx                    ! end index
INTEGER :: i_rcstartlev                ! refinement control start level
INTEGER :: ist, icheck                 ! status
INTEGER :: nverts
INTEGER :: jecp
INTEGER :: jja, jjb, jjk               ! loop indices for Moore-Penrose inverse

REAL(wp) ::   &                        ! singular values of lsq design matrix A
  &  zs(lsq_dim_unk,nproma)            ! min(lsq_dim_c,lsq_dim_unk)

REAL(wp) ::   &                        ! U matrix of SVD. Columns of U are the left
  &  zu  (lsq_dim_c,lsq_dim_c,nproma)  ! singular vectors of A

REAL(wp) ::   &                        ! TRANSPOSE of V matrix of SVD. Columns of V are
  &  zv_t(lsq_dim_unk,lsq_dim_unk,nproma) ! the right singular vectors of A.


INTEGER, PARAMETER  :: &     ! size of work array for SVD lapack routine
  &  lwork=10000
REAL(wp) ::   &              ! work array for SVD lapack routine
  &  zwork(lwork)
INTEGER  ::   &              ! work array for SVD lapack routine
  & ziwork(8*min(lsq_dim_c,lsq_dim_unk))


!DR for DEBUG purposes
! #ifdef DEBUG_COEFF LL it's used in openmp directives,
REAL(wp) :: za_debug(nproma,lsq_dim_c,lsq_dim_unk)
! #endif
!--------------------------------------------------------------------


  CALL message('mo_interpolation:lsq_compute_coeff_cell_torus', '')

  i_rcstartlev = 2

  ! get patch id
  pid = ptr_patch%id

  ! stencil size
  ptr_ncells => ptr_int_lsq%lsq_dim_stencil(:,:)

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)


!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jc,js,jec,i_startidx,i_endidx,jlv,jbv,ilc_s,ibc_s, &
!!$OMP            cc_cv,cc_cell,cc_vert,z_norm,distxy_v,z_rcarea, &
!!$OMP            delx,dely,fx,fy,fxx,fyy,fxy,jecp,nverts) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,     &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each cell, calculate weights, moments, matrix coefficients
    ! and QR decomposition of normal equation matrix
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      IF (ptr_patch%geometry_info%cell_type == 3 )THEN
        nverts  = 3
      ELSE
        nverts = ptr_patch%cells%num_edges(jc,jb)
      ENDIF
    !
    ! Gather some information about the control volume and
    ! all the cells in the stencil.
    !
    ! get line and block indices of edge vertices
      jlv(1:nverts) = ptr_patch%cells%vertex_idx(jc,jb,1:nverts)
      jbv(1:nverts) = ptr_patch%cells%vertex_blk(jc,jb,1:nverts)

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    !
    ! 1. Get CC of the control volume center and all cell centers in the stencil.
    !    In addition, get cartesian coordinates of the vertices of the control volume.
    !

    ! get CC of control volume center
      cc_cv = ptr_patch%cells%cartesian_center(jc,jb)

    ! get CC of all cell centers in the stencil
      DO js = 1, ptr_ncells(jc,jb)
        cc_cell(js) = ptr_patch%cells%cartesian_center(ilc_s(js),ibc_s(js))
      ENDDO

    ! get cartesian coordinates of edge-vertices (only for control volume)
      DO jec = 1, nverts
        cc_vert(jec) = ptr_patch%verts%cartesian(jlv(jec),jbv(jec))
      ENDDO


    !
    ! 2. Now, all points (centers, vertices) are projected onto a tangent
    !    plane having its center at the cell center of the control volume.
    !    At the same time distance vectors (r_x,r_y) are calculated between
    !    the control volume center and all other mass-points in the stencil.
    !    From these the lsq weights can be deduced.

    !
    ! a: Calculate distance vectors between control volume center and all cell
    !    centers in the stencil.
      DO js = 1, ptr_ncells(jc,jb)
        !Get the actual location of the cell w.r.t the cc_cv
        cc_cell(js) = plane_torus_closest_coordinates(cc_cv%x, cc_cell(js)%x, &
                                                      ptr_patch%geometry_info)

        !the distance vector: z coord is 0
        z_dist_g(jc,jb,js,1) = cc_cell(js)%x(1) - cc_cv%x(1)
        z_dist_g(jc,jb,js,2) = cc_cell(js)%x(2) - cc_cv%x(2)
      ENDDO


    !
    ! b: compute normalized weights for weighted least squares system
    !    The closest cell circumcenter is assigned weight of w=1.
    !
      DO js = 1, ptr_ncells(jc,jb)
        z_norm = SQRT(DOT_PRODUCT(z_dist_g(jc,jb,js,1:2),z_dist_g(jc,jb,js,1:2)))

        !
        ! weights for weighted least squares system
        !
        ptr_int_lsq%lsq_weights_c(jc,js,jb)= 1._wp/(z_norm**lsq_wgt_exp)

      ENDDO
      !
      ! Normalization
      !
      ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)=                &
        &     ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)         &
        &   / MAXVAL(ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb))



    ! 3. the following part (including calculation of moments) will only
    !    be called, if a conservative least squares reconstruction
    !    is chosen. Otherwise all moments will be equal to zero and the
    !    reconstruction simplifies to the standard non-conservative
    !    reconstruction.
      IF (llsq_rec_consv) THEN
      !
      ! a: Project control volume vertices and calculate distance vectors
      !    between cell center and vertices
      !
        DO jec=1,nverts

          !Get the actual location of the cell w.r.t the cc_cv
          cc_vert(jec) = plane_torus_closest_coordinates(cc_cv%x, cc_vert(jec)%x, &
                                                         ptr_patch%geometry_info)

          !the distance vector: z coord is 0
          distxy_v(jec,1) = cc_vert(jec)%x(1) - cc_cv%x(1)
          distxy_v(jec,2) = cc_vert(jec)%x(2) - cc_cv%x(2)

        ENDDO

      !AD: Remaining part of the code is same as the spherical part

      !
      ! b: calculate moments for given cell
      !    (calculated analytically; see Lauritzen CSLAM 09 for first 5 moments)
      !
      ! !DR: Those moments have been re-rechecked, using an alternative,
      !      quadrature-based formulation. Results have been identical up to
      !      roundoff-errors. Similarly the hat-moments have been checked. The
      !      inconsistency caused by the different projections involved do not
      !      seem to negatively effect the results.

      ! Storage docu for x^ny^m:
      ! lsq_moments(:,:,1) : x^1y^0
      ! lsq_moments(:,:,2) : x^0y^1
      ! lsq_moments(:,:,3) : x^2y^0
      ! lsq_moments(:,:,4) : x^0y^2
      ! lsq_moments(:,:,5) : x^1y^1
      ! lsq_moments(:,:,6) : x^3y^0
      ! lsq_moments(:,:,7) : x^0y^3
      ! lsq_moments(:,:,8) : x^2y^1
      ! lsq_moments(:,:,9) : x^1y^2
      !

        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! note that the distance vector distxy_v between each vertex and
          ! the center of the tangent plane are identical to the coordinates of
          ! each vertex on the tangent plane. Thus the distances between the
          ! vertices in x and y direction can be derived as follows:
          !
          ! longitudinal-distance between vertices on tangent plane
          delx(jec) = distxy_v(jecp,1) - distxy_v(jec,1)

          ! latitudinal-distance between vertices on tangent plane
          dely(jec) = distxy_v(jecp,2) - distxy_v(jec,2)


          !
          ! analytic moment calculation
          !
          ! 0: control volume area (reciprocal value)
          fx(jec) = distxy_v(jecp,1) + distxy_v(jec,1)

        ENDDO

        z_rcarea = 2._wp/DOT_PRODUCT(fx(1:nverts),dely(1:nverts))


        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! I. x^1y^0
          fx(jec) = distxy_v(jecp,1)**2              &
          &       + distxy_v(jecp,1)*distxy_v(jec,1) &
          &       + distxy_v(jec ,1)**2

          ! II. x^0y^1
          fy(jec) = distxy_v(jecp,2)**2              &
          &       + distxy_v(jecp,2)*distxy_v(jec,2) &
          &       + distxy_v(jec,2)**2

          IF ( lsq_dim_unk > 2 ) THEN

            ! III. x^2y^0
            fxx(jec) = (distxy_v(jecp,1)    + distxy_v(jec,1)       ) &
            &        * (distxy_v(jecp,1)**2 + distxy_v(jec,1)**2)

            ! IV. x^0y^2
            fyy(jec) = (distxy_v(jecp,2)    + distxy_v(jec,2)       ) &
            &        * (distxy_v(jecp,2)**2 + distxy_v(jec,2)**2)

            ! V. x^1y^1
            fxy(jec) = distxy_v(jecp,2) * (3._wp*distxy_v(jecp,1)**2                         &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + distxy_v(jec,1)**2 )     &
            &        + distxy_v(jec,2) * ( distxy_v(jecp,1)**2                               &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + 3._wp*distxy_v(jec,1)**2 )

          ENDIF ! lsq_dim_unk > 2

          IF ( lsq_dim_unk > 5 ) THEN

            ! VI.  x^3y^0
            fxxx(jec) = 5._wp*distxy_v(jec,1)**4                   &
              &       + 10._wp*distxy_v(jec,1)**3 * delx(jec)      &
              &       + 10._wp*distxy_v(jec,1)**2 * delx(jec)**2   &
              &       + 5._wp *distxy_v(jec,1)    * delx(jec)**3   &
              &       + delx(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fxxx(jec) = distxy_v(jecp,1)**4                       &
            !DR    &       + distxy_v(jec,1) * distxy_v(jecp,1)**3     &
            !DR    &       + distxy_v(jec,1)**2 * distxy_v(jecp,1)**2  &
            !DR    &       + distxy_v(jec,1)**3 * distxy_v(jecp,1)     &
            !DR    &       + distxy_v(jec,1)**4


            ! VII. x^0y^3
            fyyy(jec) = 5._wp*distxy_v(jec,2)**4                   &
              &       + 10._wp*distxy_v(jec,2)**3 * dely(jec)      &
              &       + 10._wp*distxy_v(jec,2)**2 * dely(jec)**2   &
              &       + 5._wp *distxy_v(jec,2)    * dely(jec)**3   &
              &       + dely(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fyyy(jec) = distxy_v(jecp,2)**4                       &
            !DR    &       + distxy_v(jec,2)    * distxy_v(jecp,2)**3  &
            !DR    &       + distxy_v(jec,2)**2 * distxy_v(jecp,2)**2  &
            !DR    &       + distxy_v(jec,2)**3 * distxy_v(jecp,2)     &
            !DR    &       + distxy_v(jec,2)**4


            ! VIII. x^2y^1
            fxxy(jec) = 4._wp*distxy_v(jecp,1)**3 *distxy_v(jecp,2)                    &
              &       + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jecp,2)   &
              &       + 2._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jecp,2) &
              &       +       distxy_v(jec,1)**3 * distxy_v(jecp,2)                    &
              &       +       distxy_v(jecp,1)**3 * distxy_v(jec,2)                    &
              &       + 2._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jec,2)    &
              &       + 3._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jec,2)  &
              &       + 4._wp*distxy_v(jec,1)**3 * distxy_v(jec,2)

            ! IX. x^1y^2
            fxyy(jec) = 6._wp*distxy_v(jecp,1)**2 * distxy_v(jecp,2)**2                     &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jecp,2)**2              &
              &   +       distxy_v(jec,1)**2 * distxy_v(jecp,2)**2                          &
              &   + 3._wp*distxy_v(jecp,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)            &
              &   + 4._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)*distxy_v(jecp,2) &
              &   + 3._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)             &
              &   +       distxy_v(jecp,1)**2 * distxy_v(jec,2)**2                          &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)**2               &
              &   + 6._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)**2

          ENDIF ! lsq_dim_unk > 5


        ENDDO ! loop over nverts

        ptr_int_lsq%lsq_moments(jc,jb,1)= z_rcarea/6._wp*DOT_PRODUCT(fx(1:nverts),dely(1:nverts))
        ptr_int_lsq%lsq_moments(jc,jb,2)=-z_rcarea/6._wp*DOT_PRODUCT(fy(1:nverts),delx(1:nverts))

        IF ( lsq_dim_unk > 2 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,3) =  z_rcarea/12._wp * &
            & DOT_PRODUCT(fxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,4) = -z_rcarea/12._wp * &
            & DOT_PRODUCT(fyy(1:nverts),delx(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,5) =  z_rcarea/24._wp * &
            & DOT_PRODUCT(fxy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 2


        IF ( lsq_dim_unk > 5 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,6) =  z_rcarea/20._wp * &
            & DOT_PRODUCT(fxxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,7) = -z_rcarea/20._wp * &
            & DOT_PRODUCT(fyyy(1:nverts),delx(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,8) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxxy(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,9) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxyy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 5

      END IF  ! llsq_rec_consv

    END DO  ! loop over cells

  END DO  ! loop over blocks
!!$OMP END DO NOWAIT
!! For unknown reasons, closing the parallel section here is needed to get the above
!! loop parallelized.
!!$OMP END PARALLEL

  DO jb = 1, lsq_dim_c
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_weights_c(:,jb,:))
  ENDDO
  DO jb = 1, lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments(:,:,jb))
  ENDDO


!$OMP PARALLEL PRIVATE(jb,jc,js,ju,jja,jjb,jjk,i_startidx,i_endidx,ilc_s,ibc_s, &
!$OMP            z_lsq_mat_c,zs,zu,zv_t,zwork,ziwork,ist,icheck,za_debug, &
!$OMP            z_qmat,z_rmat,cnt,jrow,nel)
!$OMP DO ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! 4. for each cell, calculate LSQ design matrix A
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
        ! Take care that z_lsq_mat_c isn't singular
        z_lsq_mat_c(jc,:,:) = 0.0_wp
        DO js = 1, MIN(lsq_dim_unk, lsq_dim_c)
          z_lsq_mat_c(jc,js,js) = 1.0_wp
        ENDDO
        CYCLE
      ENDIF

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    ! Calculate full moments lsq_moments_hat(ilc_s(js),ibc_s(js),ju)
    !
    ! Storage docu for x^ny^m:
    ! lsq_moments_hat(:,:,:,1) : \hat{x^1y^0}
    ! lsq_moments_hat(:,:,:,2) : \hat{x^0y^1}
    ! lsq_moments_hat(:,:,:,3) : \hat{x^2y^0}
    ! lsq_moments_hat(:,:,:,4) : \hat{x^0y^2}
    ! lsq_moments_hat(:,:,:,5) : \hat{x^1y^1}
    ! lsq_moments_hat(:,:,:,6) : \hat{x^3y^0}
    ! lsq_moments_hat(:,:,:,7) : \hat{x^0y^3}
    ! lsq_moments_hat(:,:,:,8) : \hat{x^2y^1}
    ! lsq_moments_hat(:,:,:,9) : \hat{x^1y^2}
    !
      DO js = 1, ptr_ncells(jc,jb)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,1) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1) + z_dist_g(jc,jb,js,1)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,2) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2) + z_dist_g(jc,jb,js,2)

        IF (lsq_dim_unk > 2) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,3) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
           &    + z_dist_g(jc,jb,js,1)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,4) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2) &
           &    + z_dist_g(jc,jb,js,2)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,5) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)                              &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)        &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)        &
           &    + z_dist_g(jc,jb,js,1) * z_dist_g(jc,jb,js,2)
        ENDIF

        IF ( lsq_dim_unk > 5 ) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,6) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),6)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,1) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)                       &
            &    * z_dist_g(jc,jb,js,1)**2                                                    &
            &    + z_dist_g(jc,jb,js,1)**3

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,7) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),7)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,2) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)                       &
            &    * z_dist_g(jc,jb,js,2)**2                                                    &
            &    + z_dist_g(jc,jb,js,2)**3

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,8) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),8)                             &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,2)       &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,1) &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
            &    * z_dist_g(jc,jb,js,2)                                                       &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)**2    &
            &    + z_dist_g(jc,jb,js,1)**2 * z_dist_g(jc,jb,js,2)


          ptr_int_lsq%lsq_moments_hat(jc,jb,js,9) =                                          &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),9)                            &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,1)      &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,2)&
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2)&
            &    * z_dist_g(jc,jb,js,1)                                                      &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)**2   &
            &    + z_dist_g(jc,jb,js,2)**2 * z_dist_g(jc,jb,js,1)
        ENDIF
      ENDDO


      ! loop over rows of lsq design matrix (all cells in the stencil)
      DO js = 1, ptr_ncells(jc,jb)
        ! loop over columns of lsq design matrix (number of unknowns)
        DO ju = 1,lsq_dim_unk

          z_lsq_mat_c(jc,js,ju) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                         &
           &    * (ptr_int_lsq%lsq_moments_hat(jc,jb,js,ju)-ptr_int_lsq%lsq_moments(jc,jb,ju))

        END DO
       END DO
       IF(ptr_ncells(jc,jb) < lsq_dim_c) THEN
         z_lsq_mat_c(jc,lsq_dim_c,:) = 0.0_wp
       ENDIF

    ENDDO   ! loop over cells



    !
    ! compute QR decomposition and Singular Value Decomposition (SVD)
    ! of least squares design matrix A. For the time being both methods are
    ! retained.

    !
    ! 5a. QR-factorization of design matrix A
    !
    IF (.NOT. advection_config(pid)%llsq_svd) THEN
!CDIR NOIEXPAND
    CALL qrdec(lsq_dim_c, lsq_dim_unk, i_startidx, & ! in
     &         i_endidx, z_lsq_mat_c,              & ! in
     &         z_qmat, z_rmat)                       ! out


    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! 7. Save transposed Q-Matrix
      ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,1:lsq_dim_c,jb)  =  &
       &      TRANSPOSE(z_qmat(jc,1:lsq_dim_c,1:lsq_dim_unk))

      ! 8. Save R-Matrix
      !
      ! a. Save reciprocal values of the diagonal elements
      !
      DO ju = 1,lsq_dim_unk
        ptr_int_lsq%lsq_rmat_rdiag_c(jc,ju,jb) = 1._wp/z_rmat(jc,ju,ju)
      ENDDO

      !
      ! b. Save upper triangular elements without diagonal elements in a 1D-array
      !    (starting from the bottom right)
      !
      cnt = 1
      DO jrow = lsq_dim_unk-1,1,-1
        ! number of elements to store
        nel = lsq_dim_unk - jrow
        ptr_int_lsq%lsq_rmat_utri_c(jc,cnt:cnt+nel-1,jb) = z_rmat(jc,jrow,jrow+1:lsq_dim_unk)
        cnt = cnt + nel
      ENDDO


      ! Multiply ith column of the transposed Q-matrix (corresponds to the
      ! different members of the stencil) with the ith weight. This avoids
      ! multiplication of the RHS of the LSQ-System with this weight during
      ! runtime.
      DO js = 1,lsq_dim_c
        ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  =                  &
          &                ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  &
          &              * ptr_int_lsq%lsq_weights_c(jc,js,jb)
      ENDDO

    END DO  ! loop over cells

    ELSE   ! llsq_svd=.TRUE.

    !
    ! 5b. Singular value decomposition of lsq design matrix A
    ! !!! does not vectorize, unfortunately !!!
    ist = 0
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! A = U * SIGMA * transpose(V)
      !
      ! z_lsq_mat_c : M x N least squares design matrix A            (IN)
      ! zu          : M x M orthogonal matrix U                      (OUT)
      ! zv_t        : N x N orthogonal matrix V                      (OUT)
      ! zs          : min(M,N) Singular values of A                  (OUT)
      ! zwork       : workspace(1,LWORK)                             (OUT)
      ! lwork       : 3*min(M,N)                                     (IN)
      !              + max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N))  (IN)
      ! iwork       : workspace(8*min(M,N))                          (IN)

      CALL DGESDD('A',                 & !in
        &         lsq_dim_c,           & !in
        &         lsq_dim_unk,         & !in
        &         z_lsq_mat_c(jc,:,:), & !inout Note: destroyed on output
        &         lsq_dim_c,           & !in
        &         zs(:,jc),            & !out
        &         zu(:,:,jc),          & !out
        &         lsq_dim_c,           & !in
        &         zv_t(:,:,jc),        & !out
        &         lsq_dim_unk,         & !in
        &         zwork,               & !out
        &         lwork,               & !in
        &         ziwork,              & !inout
        &         icheck               ) !out
      ist = ist + icheck
    ENDDO
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:lsq_compute_coeff_cell_torus',   &
        &             'singular value decomposition failed')
    ENDIF

    ! compute Moore-Penrose inverse
    ! INVERSE(A):: V * INVERSE(SIGMA) * TRANSPOSE(U) and store
    ! note that the ith column is multiplied with the ith weight
    ! in order to avoid the weighting of the r.h.s. during runtime.
    DO jja = 1, lsq_dim_unk
      DO jjb = 1, lsq_dim_c
        DO jjk = 1, lsq_dim_unk
          DO jc = i_startidx, i_endidx
            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
            ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb) =            &
              &  ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb)         &
              &  + zv_t(jjk,jja,jc) /zs(jjk,jc) * zu(jjb,jjk,jc)  &
              &  * ptr_int_lsq%lsq_weights_c(jc,jjb,jb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!!!!!DEBUG !!!
#ifdef DEBUG_COEFF
    za_debug(:,:,:)  = 0._wp
    ! re-COMPUTE A:: U * SIGMA * TRANSPOSE(V)  !!! Funktioniert
    DO jja = 1, lsq_dim_c
      DO jjb = 1, lsq_dim_unk
        DO jjk = 1, lsq_dim_unk  !lsq_dim_c
          DO jc = i_startidx, i_endidx
            za_debug(jc,jja,jjb) = za_debug(jc,jja,jjb)  &
              &  + zu(jja,jjk,jc) * zs(jjk,jc) * zv_t(jjk,jjb,jc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    write(0,*) "za_debug(10,:,:)- z_lsq_mat_c(10,:,:)",za_debug(10,:,:)- z_lsq_mat_c(10,:,:)
    write(0,*) "za_debug(55,:,:)- z_lsq_mat_c(55,:,:)",za_debug(55,:,:)- z_lsq_mat_c(55,:,:)
    write(0,*) "zs(:,10) ",zs(:,10)
    write(0,*) "zs(:,55) ",zs(:,55)
    write(0,*) "ptr_int_lsq%lsq_weights_c(10,:,jb)",ptr_int_lsq%lsq_weights_c(10,:,jb)
    write(0,*) "ptr_int_lsq%lsq_weights_c(55,:,jb)",ptr_int_lsq%lsq_weights_c(55,:,jb)
#endif
!!!! END DEBUG !!!

    ENDIF  ! llsq_svd

  END DO  ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  DO ju = 1, lsq_dim_unk
    DO jc = 1, lsq_dim_c
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments_hat(:,:,jc,ju))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_pseudoinv(:,ju,jc,:))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_qtmat_c(:,ju,jc,:))
    ENDDO
  ENDDO

  DO jc = 1, UBOUND(ptr_int_lsq%lsq_rmat_utri_c, 2)
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_utri_c(:,jc,:))
  ENDDO
  DO ju = 1,lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_rdiag_c(:,ju,:))
  ENDDO

END SUBROUTINE lsq_compute_coeff_cell_torus



!-------------------------------------------------------------------------
!!
!! This routine initializes the coefficients used
!! for interpolations needed for scalars. The original routines were aw_int_coeff
!! and cell2edge_lin_int_coeff
!!
!! ptr_int_state% c_lin_e
!!                verts_aw_cells
!!                cells_aw_verts
!!                e_inn_c
!!
SUBROUTINE scalar_int_coeff( ptr_patch, ptr_int_state )
!

TYPE(t_patch), INTENT(in) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int_state

INTEGER :: nlen, nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, npromz_v
INTEGER :: jc, je, jb, jv  ! integer over edges, blocks and levels
INTEGER :: ile, ibe, &
           ilc, ibc, ilc1, ilc2, ibc1, ibc2, idx_ve,&
           ilv, ibv, ilv1, ilv2, ibv1, ibv2, idx_ce
INTEGER :: i_startblk                ! start block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c
  npromz_c = ptr_patch%npromz_c
  nblks_e  = ptr_patch%nblks_e
  npromz_e = ptr_patch%npromz_e
  nblks_v  = ptr_patch%nblks_v
  npromz_v = ptr_patch%npromz_v


!$OMP PARALLEL PRIVATE(i_startblk)

  ! a) cell to edge averages
  !-------------------------
  ! The calculation cannot be done for boundary edges
  i_startblk = ptr_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, 2)

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

      ! inverse distance averaging (former subroutine cell2edge_lin_int_coeff)
      ! For the hexagonal grid, this is also the direct distance averaging
      ! because the edge is exactly half way between the cells
      ! (needed for proper bracket formalism)
      ptr_int_state%c_lin_e(je,1,jb) = ptr_patch%edges%edge_cell_length(je,jb,2)/&
                                           ptr_patch%edges%dual_edge_length(je,jb)
      ptr_int_state%c_lin_e(je,2,jb) = 1._wp - ptr_int_state%c_lin_e(je,1,jb)

    ENDDO
  ENDDO
!$OMP END DO



  ! b) vert to cell averagings, edge to cell inner product
  !-------------------------------------------------------
  ! loop over all blocks and cells

!$OMP DO PRIVATE(jb,jc,je,jv,nlen,ile,ibe,idx_ce,ilv1,ilv2,ibv1,ibv2,&
!$OMP ilv,ibv) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

    DO jc = 1, nlen

       IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

       ptr_int_state%verts_aw_cells(jc,:,jb) = 0.0_wp

       DO je = 1, ptr_patch%cells%num_edges(jc,jb)

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)
          IF ( ptr_patch%edges%cell_idx(ile,ibe,1) == jc .AND. &
               ptr_patch%edges%cell_blk(ile,ibe,1) == jb ) THEN
               idx_ce = 1
          ELSE
               idx_ce = 2
          ENDIF

          ptr_int_state%e_inn_c(jc,je,jb) = &
              ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
              ptr_patch%edges%primal_edge_length(ile,ibe)/&
              ptr_patch%cells%area(jc,jb)

          ilv1 = ptr_patch%edges%vertex_idx(ile,ibe,1)
          ibv1 = ptr_patch%edges%vertex_blk(ile,ibe,1)
          ilv2 = ptr_patch%edges%vertex_idx(ile,ibe,2)
          ibv2 = ptr_patch%edges%vertex_blk(ile,ibe,2)

          DO jv = 1, ptr_patch%cells%num_edges(jc,jb)
            ilv = ptr_patch%cells%vertex_idx(jc,jb,jv)
            ibv = ptr_patch%cells%vertex_blk(jc,jb,jv)

            IF (ilv == ilv1 .AND. ibv == ibv1) THEN
              ptr_int_state%verts_aw_cells(jc,jv,jb) =   &
                ptr_int_state%verts_aw_cells(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,1)
            ENDIF
            IF (ilv == ilv2 .AND. ibv == ibv2) THEN
              ptr_int_state%verts_aw_cells(jc,jv,jb)  =  &
                ptr_int_state%verts_aw_cells(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,2)
            ENDIF

          ENDDO

       ENDDO

    ENDDO !loop over all cells

  ENDDO   !loop over all blocks
!$OMP END DO



  ! c) cells to verts averagings, edge to verts averagings
  !-------------------------------------------------------
  ! loop over all blocks and verts

  i_startblk = ptr_patch%verts%start_blk(2,1)
!$OMP DO PRIVATE(jb,jc,je,jv,i_startidx,i_endidx,ile,ibe,idx_ve,ilc,ibc, &
!$OMP            ilc1,ilc2,ibc1,ibc2 ) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, 2)

    DO jv = i_startidx, i_endidx

       IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

       ptr_int_state%cells_aw_verts(jv,:,jb) = 0.0_wp

       DO je = 1, ptr_patch%verts%num_edges(jv,jb)

          ile = ptr_patch%verts%edge_idx(jv,jb,je)
          ibe = ptr_patch%verts%edge_blk(jv,jb,je)
          IF ( ptr_patch%edges%vertex_idx(ile,ibe,1) == jv .AND. &
               ptr_patch%edges%vertex_blk(ile,ibe,1) == jb ) THEN
               idx_ve = 1
          ELSE
               idx_ve = 2
          ENDIF

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          DO jc = 1, ptr_patch%verts%num_edges(jv,jb)
            ilc = ptr_patch%verts%cell_idx(jv,jb,jc)
            ibc = ptr_patch%verts%cell_blk(jv,jb,jc)

            IF (ilc == ilc1 .AND. ibc == ibc1) THEN
              ptr_int_state%cells_aw_verts(jv,jc,jb) =   &
                ptr_int_state%cells_aw_verts(jv,jc,jb) + &
                0.5_wp/ptr_patch%verts%dual_area(jv,jb) *             &
                ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve)*&
                ptr_patch%edges%edge_cell_length(ile,ibe,1)
            ENDIF
            IF (ilc == ilc2 .AND. ibc == ibc2) THEN
              ptr_int_state%cells_aw_verts(jv,jc,jb)  =  &
                ptr_int_state%cells_aw_verts(jv,jc,jb) + &
                0.5_wp/ptr_patch%verts%dual_area(jv,jb) *             &
                ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve)*&
                ptr_patch%edges%edge_cell_length(ile,ibe,2)
            ENDIF

          ENDDO

       ENDDO

    ENDDO !loop over all cells

  ENDDO   !loop over all blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int_state%c_lin_e)
  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%verts_aw_cells)
  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_inn_c)
  CALL sync_patch_array(SYNC_V,ptr_patch,ptr_int_state%cells_aw_verts)

END SUBROUTINE scalar_int_coeff
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!>
!! Calls routines to calculate the weighting coefficients for bilinear
!! edge-to-cell interpolation depending on grid geometry.
!!
SUBROUTINE bln_int_coeff_e2c( ptr_patch, ptr_int_state )
  TYPE(t_patch),     INTENT(inout) :: ptr_patch
  TYPE(t_int_state), INTENT(inout) :: ptr_int_state

  CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs_lsq_bln:bln_int_coeff_e2c'

  !
  SELECT CASE(ptr_patch%geometry_info%geometry_type)

  CASE (planar_torus_geometry)
    CALL flat_scalar_coeffs( ptr_patch, ptr_int_state )
    CALL vector_coeffs(ptr_patch, ptr_int_state)

  CASE (sphere_geometry)
    CALL spherical_scalar_coeffs( ptr_patch, ptr_int_state )
    CALL vector_coeffs(ptr_patch, ptr_int_state)

  CASE DEFAULT
    CALL finish(method_name, "Undefined geometry type")

  END SELECT

END SUBROUTINE bln_int_coeff_e2c
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!
!! Computes the weighting coefficients for bilinear edge-to-cell interpolation.
!!
!! Results are stored in ptr_int_state\\%e_bln_c_s
!!
SUBROUTINE spherical_scalar_coeffs ( ptr_patch, ptr_int_state )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int_state
!

INTEGER :: jc, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3

REAL(wp) :: xtemp,ytemp,wgt(3),xloc,yloc,x(3),y(3), &
            pollat,pollon

!-----------------------------------------------------------------------

rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,yloc,xloc,pollat,pollon,ile1,ibe1,&
!$OMP            ile2,ibe2,ile3,ibe3,xtemp,ytemp,wgt,x,y) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    yloc = ptr_patch%cells%center(jc,jb)%lat
    xloc = ptr_patch%cells%center(jc,jb)%lon

    ! Rotate local point into the equator for better accuracy of bilinear weights
    IF (yloc >= 0._wp) THEN
      pollat = yloc - pi2/4._wp
    ELSE
      pollat = yloc + pi2/4._wp
    ENDIF
    pollon = xloc

    CALL rotate_latlon( yloc, xloc, pollat, pollon )

    !  get the line and block indices of the cell edges

    ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
    ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
    ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
    ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
    ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
    ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

    ! x and y are the zonal and meridional distances from the local
    ! cell point (ignoring the earth's radius, which drops out anyway)

    xtemp = ptr_patch%edges%center(ile1,ibe1)%lon
    ytemp = ptr_patch%edges%center(ile1,ibe1)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(1)  = ytemp-yloc
    x(1)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(1) >  3.5_wp) x(1) = x(1) - pi2
    IF (x(1) < -3.5_wp) x(1) = x(1) + pi2

    xtemp = ptr_patch%edges%center(ile2,ibe2)%lon
    ytemp = ptr_patch%edges%center(ile2,ibe2)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(2)  = ytemp-yloc
    x(2)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(2) >  3.5_wp) x(2) = x(2) - pi2
    IF (x(2) < -3.5_wp) x(2) = x(2) + pi2

    xtemp = ptr_patch%edges%center(ile3,ibe3)%lon
    ytemp = ptr_patch%edges%center(ile3,ibe3)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(3)  = ytemp-yloc
    x(3)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(3) >  3.5_wp) x(3) = x(3) - pi2
    IF (x(3) < -3.5_wp) x(3) = x(3) + pi2

    ! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
    ! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
    ! are not aliased into a checkerboard pattern between upward- and downward
    ! directed cells. The third condition is sum(w(i)) = 1. Analytical elimination yields...

    IF (ABS(x(2)-x(1)) > 1.e-11_wp .AND. ABS(y(3)-y(1)) > 1.e-11_wp ) THEN
      wgt(3) = 1.0_wp/( (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))/(x(2)-x(1)) ) * &
                  ( -y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) )
      wgt(2) = (-x(1) - wgt(3)*(x(3)-x(1)))/(x(2)-x(1))
      wgt(1) = 1.0_wp - wgt(2) - wgt(3)
    ELSE
      wgt(2) = 1.0_wp/( (y(2)-y(1)) - (x(2)-x(1))*(y(3)-y(1))/(x(3)-x(1)) ) * &
                  ( -y(1) + x(1)*(y(3)-y(1))/(x(3)-x(1)) )
      wgt(3) = (-x(1) - wgt(2)*(x(2)-x(1)))/(x(3)-x(1))
      wgt(1) = 1.0_wp - wgt(2) - wgt(3)
    ENDIF

    ! Store results in ptr_int_state%e_bln_c_s
    ptr_int_state%e_bln_c_s(jc,1,jb) = wgt(1)
    ptr_int_state%e_bln_c_s(jc,2,jb) = wgt(2)
    ptr_int_state%e_bln_c_s(jc,3,jb) = wgt(3)

  ENDDO !cell loop

END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_s)


END SUBROUTINE spherical_scalar_coeffs
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!
!! Computes the vector weighting coefficients using the scalar part computed earlier
!!
SUBROUTINE vector_coeffs ( ptr_patch, ptr_int_state )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int_state
!

INTEGER :: jc, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3

!-----------------------------------------------------------------------

! Now compute vector interpolation weights: These take the normal and tangential
! wind components at the edges and reconstruct u and v at the cell midpoints

rl_start = 2
rl_end = min_rlcell

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ile1,ibe1,ile2,ibe2,&
!$OMP ile3,ibe3) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    !  get the line and block indices of the cell edges

    ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
    ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
    ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
    ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
    ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
    ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

    IF (ptr_patch%edges%cell_idx(ile1,ibe1,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile1,ibe1,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v1
      ptr_int_state%e_bln_c_u(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,1)%v1

      ptr_int_state%e_bln_c_v(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v2
      ptr_int_state%e_bln_c_v(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile1,ibe1,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile1,ibe1,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v1
      ptr_int_state%e_bln_c_u(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,2)%v1

      ptr_int_state%e_bln_c_v(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v2
      ptr_int_state%e_bln_c_v(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,2)%v2

    ENDIF

    IF (ptr_patch%edges%cell_idx(ile2,ibe2,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile2,ibe2,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,1)%v1
      ptr_int_state%e_bln_c_u(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,1)%v1

      ptr_int_state%e_bln_c_v(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,1)%v2
      ptr_int_state%e_bln_c_v(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile2,ibe2,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile2,ibe2,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,2)%v1
      ptr_int_state%e_bln_c_u(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,2)%v1

      ptr_int_state%e_bln_c_v(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,2)%v2
      ptr_int_state%e_bln_c_v(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,2)%v2

    ENDIF

    IF (ptr_patch%edges%cell_idx(ile3,ibe3,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile3,ibe3,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,1)%v1
      ptr_int_state%e_bln_c_u(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,1)%v1

      ptr_int_state%e_bln_c_v(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,1)%v2
      ptr_int_state%e_bln_c_v(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile3,ibe3,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile3,ibe3,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,2)%v1
      ptr_int_state%e_bln_c_u(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,2)%v1

      ptr_int_state%e_bln_c_v(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,2)%v2
      ptr_int_state%e_bln_c_v(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,2)%v2

    ENDIF
  ENDDO !cell loop

END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_u)
CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_v)

END SUBROUTINE vector_coeffs
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!! Computes the weighting coefficients for bilinear edge-to-cell interpolation for
!! flat geometry with equilateral triangular cells
!!
SUBROUTINE flat_scalar_coeffs ( ptr_patch, ptr_int_state )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int_state
!

INTEGER  :: jc, jb
INTEGER  :: rl_start, rl_end
INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
REAL(wp) :: wgt

!-----------------------------------------------------------------------

rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

wgt = 1._wp/3._wp

!
! loop through all patch cells (and blocks)
!
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)
    DO jc = i_startidx, i_endidx

      ! Simple for torus
      ptr_int_state%e_bln_c_s(jc,1:3,jb) = wgt

    ENDDO !cell loop
  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_s)


END SUBROUTINE flat_scalar_coeffs
!-------------------------------------------------------------------------


END MODULE mo_intp_coeffs_lsq_bln
