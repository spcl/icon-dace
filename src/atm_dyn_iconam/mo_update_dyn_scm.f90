!
! Updates dynamical fields with slow physics tendencies in the special case
! that dynamics are switched off (ldynamics=F). An option is added where Coriolis force
! should be used (lcoriolis=T).
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

MODULE mo_update_dyn_scm

  USE mo_kind,               ONLY: wp
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_types,     ONLY: t_nh_state
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_physical_constants, ONLY: p0ref, rd, cvd_o_rd
  USE mo_impl_constants,     ONLY: min_rlcell_int, min_rledge_int
  USE mo_sync,               ONLY: SYNC_E, SYNC_C, sync_patch_array
  USE mo_dynamics_config,    ONLY: lcoriolis
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e




  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: add_slowphys_scm, rbf_coeff_scm


CONTAINS

!-------------------------------------------------------------------------
!>
!! Updates dynamical fields with slow physics tendencies
!!
!!
  SUBROUTINE add_slowphys_scm(p_nh, p_patch, p_int, nnow, nnew, dtime)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state), TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),     TARGET, INTENT(INOUT) :: p_patch

    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime

    INTEGER :: nlev                  ! number of vertical (full) levels

    INTEGER :: jc, je, jk, jb        ! loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk  ! start and end block
    INTEGER :: i_startidx, i_endidx  ! start and end indices
    INTEGER :: i_nchdom
    REAL(wp) :: p_diag_vt, ddt_corio_vn
    ! Pointers
    INTEGER, POINTER :: iqidx(:,:,:), iqblk(:,:,:) ! to quad edge indices

!-----------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk
    i_nchdom = MAX(1,p_patch%n_childdom)

    !$ACC DATA PRESENT(p_nh) IF(i_am_accel_node)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    rl_start = grf_bdywidth_c+1 
    rl_end   = min_rlcell_int 

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! exner (update)
          p_nh%prog(nnew)%exner(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)  &
             &                            + dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)

          ! diagnose theta_v from updated exner
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = (p0ref/(rd*p_nh%prog(nnew)%rho(jc,jk,jb)))  &
             &                              * p_nh%prog(nnew)%exner(jc,jk,jb)**cvd_o_rd

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP ENDDO NOWAIT

    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk  (rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,p_diag_vt,ddt_corio_vn)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          ! add Coriolis force
          IF ( lcoriolis ) THEN
            ! RBF reconstruction of tangential wind component
            p_diag_vt = &
              p_int%rbf_vec_coeff_e(1,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) + &
              p_int%rbf_vec_coeff_e(2,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) + &
              p_int%rbf_vec_coeff_e(3,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) + &
              p_int%rbf_vec_coeff_e(4,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            !Coriolis tendency
            ddt_corio_vn  = - p_diag_vt * p_patch%edges%f_e(je,jb)
          ELSE
            ddt_corio_vn = 0._wp
          END IF

          ! vn (update)
          p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)  &
            &                + dtime * (ddt_corio_vn+p_nh%diag%ddt_vn_phy(je,jk,jb))

        ENDDO  ! je
      ENDDO  ! jk
      !$ACC END PARALLEL

    ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL

    ! Synchronize updated prognostic variables
    CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%exner)
    CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn)

    !$ACC END DATA

  END SUBROUTINE add_slowphys_scm


!-------------------------------------------------------------------------
!>
!! SCM: initialize RBF coefficients as well as blocks and indices on a single central point
!!
!! The interpolation of u/v to vn and back each time step results in differences in the 32
!! points of the SCM.  Therefore the coefficients are selected such that only a single point is choosen
!! for the RBF interpolation.  As a random choice the first point is taken (i_startidx,i_startblk).
!!

  SUBROUTINE rbf_coeff_scm( ptr_patch, ptr_int )

  TYPE(t_patch)    , INTENT(inout) :: ptr_patch

  TYPE(t_int_state), INTENT(inout) :: ptr_int

  INTEGER :: jc, jb                ! loop indices
  INTEGER :: i_startblk, i_endblk  ! start and end block
  INTEGER :: i_startidx, i_endidx  ! start and end indices
  INTEGER :: i_rcstartlev          ! refinement control start level

!-----------------------------------------------------------------------

  i_rcstartlev = 2

  ! values for the blocking
  i_endblk  = ptr_patch%nblks_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, i_rcstartlev)

    DO jc = i_startidx, i_endidx
      ptr_int%rbf_vec_idx_c  (:,  jc,jb) = ptr_int%rbf_vec_idx_c  (:,  i_startidx,i_startblk)
      ptr_int%rbf_vec_blk_c  (:,  jc,jb) = ptr_int%rbf_vec_blk_c  (:,  i_startidx,i_startblk)
      ptr_int%rbf_vec_coeff_c(:,:,jc,jb) = ptr_int%rbf_vec_coeff_c(:,:,i_startidx,i_startblk)
    END DO

  END DO

  END SUBROUTINE rbf_coeff_scm

END MODULE mo_update_dyn_scm



