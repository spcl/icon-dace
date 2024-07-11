! This module contains subroutines that update the prep_adv state
! with meaningful values for standalone advection runs.
! Note that in real case runs, the prep_adv state is updated by the
! dynamical core (nh_solve).
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

MODULE mo_prepadv_util

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,     ONLY: min_rledge_int, min_rlcell_int, min_rlcell
  USE mo_timer,              ONLY: timers_level, timer_start, timer_stop, timer_prep_tracer
  USE mo_fortran_tools,      ONLY: init, set_acc_host_or_device


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prepare_tracer

CONTAINS

  !--------------------------------------------------------------------------
  !>
  !! Diagnose mass fluxes and trajectory-velocities for the 
  !! standalone tracer transport scheme.
  !!
  !! Note that upper boundary fluxes (when using horizontal/vertical nesting)
  !! are set in mo_nh_nest_utilities/boundary_interpolation.
  !!
  SUBROUTINE prepare_tracer( p_patch, p_now, p_new, p_metrics, p_nh_diag, &
    &                        p_vn_traj, p_mass_flx_me, p_mass_flx_ic, lacc)

    TYPE(t_patch),     INTENT(IN) :: p_patch
    TYPE(t_nh_prog),   INTENT(IN) :: p_now, p_new
    TYPE(t_nh_metrics),INTENT(IN) :: p_metrics
    TYPE(t_nh_diag),   INTENT(IN) :: p_nh_diag

    REAL(wp),          INTENT(INOUT) :: p_vn_traj(:,:,:)      ! (nproma,  nlev,p_patch%nblks_e)
    REAL(wp),          INTENT(INOUT) :: p_mass_flx_me(:,:,:)  ! (nproma,  nlev,p_patch%nblks_e)
    REAL(wp),          INTENT(INOUT) :: p_mass_flx_ic(:,:,:)  ! (nproma,nlevp1,p_patch%nblks_c)
    LOGICAL, OPTIONAL, INTENT(IN   ) :: lacc                  ! flag to run on GPU

    ! local variables
    REAL(wp):: w_tavg               !< contravariant vertical velocity at n+\alpha

    INTEGER :: je, jc, jk, jb       !< loop indices and domain ID
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart_e, i_rlend_e, i_rlstart_c, i_rlend_c
    INTEGER :: nlev, nlevp1         !< number of full and half levels
    LOGICAL :: lzacc                ! non-optional version of lacc

   !--------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (timers_level > 5) CALL timer_start(timer_prep_tracer)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! refinement control start/end level for cells
    i_rlstart_c = 1
    i_rlend_c   = min_rlcell_int

    ! refinement control start/end level for edges
    i_rlstart_e = 1
    i_rlend_e   = min_rledge_int - 2

    !$ACC DATA PRESENT(p_vn_traj, p_mass_flx_me, p_mass_flx_ic) IF(lzacc)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk  = p_patch%edges%start_block(i_rlstart_e)
    i_endblk    = p_patch%edges%end_block(i_rlend_e)

    !
    ! contravariant normal velocites at n+1/2
    ! necessary for computation of backward trajectories.
    !
    ! horizontal (contravariant) mass flux at full level edges
    !
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart_e, i_rlend_e )

      ! reset mass fluxes and trajectory-velocities to start new integration sweep
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      p_vn_traj    (:,:,jb) = 0._wp
      p_mass_flx_me(:,:,jb) = 0._wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev
!DIR$ IVDEP
        DO je = i_startidx, i_endidx

          ! trajectory-velocity
          p_vn_traj(je,jk,jb) = 0.5_wp * ( p_now%vn(je,jk,jb) + p_new%vn(je,jk,jb) )

          ! mass flux
          p_mass_flx_me(je,jk,jb) = p_nh_diag%mass_fl_e(je,jk,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO NOWAIT

    IF (p_test_run) THEN ! Reset also halo points to zero
      CALL init(p_vn_traj    (:,:,i_endblk+1:p_patch%nblks_e), lacc=lzacc)
      CALL init(p_mass_flx_me(:,:,i_endblk+1:p_patch%nblks_e), lacc=lzacc)
!$OMP BARRIER
    ENDIF


    i_startblk   = p_patch%cells%start_block(i_rlstart_c)
    i_endblk     = p_patch%cells%end_block(i_rlend_c)

    !
    ! vertical mass flux at half level centers
    !
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,w_tavg) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

      ! reset mass fluxes and trajectory-velocities to start new integration sweep
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      p_mass_flx_ic(:,:,jb) = 0._wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(w_tavg)
      DO jk = 1, nlevp1
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx

          ! Note(DR): This is somewhat inconsistent since for horizontal trajectories
          ! v_n at n+1/2 is used.
          ! w_concorr_c at TL n+1
          w_tavg = p_metrics%vwind_expl_wgt(jc,jb)*p_now%w(jc,jk,jb)  &
            &    + p_metrics%vwind_impl_wgt(jc,jb)*p_new%w(jc,jk,jb)  &
            &    - p_nh_diag%w_concorr_c(jc,jk,jb)

          p_mass_flx_ic(jc,jk,jb) = p_nh_diag%rho_ic(jc,jk,jb) * w_tavg

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO

    IF (p_test_run) THEN ! Reset also halo points to zero
      CALL init(p_mass_flx_ic(:,:,i_endblk+1:p_patch%nblks_c), lacc=lzacc)
!$OMP BARRIER
    ENDIF
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

    IF (timers_level > 5) CALL timer_stop(timer_prep_tracer)

  END SUBROUTINE prepare_tracer

END MODULE mo_prepadv_util

