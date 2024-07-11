! Contains correction term for moist thermodynamics (moisture-dependence of heat capacities)
!
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

MODULE mo_nh_moist_thdyn

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_prepadv_types,       ONLY: t_prepare_adv
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_run_config,          ONLY: ntracer, iqv, iqc, iqi, iqr, iqs, iqg, iqh, iqgl, iqhl
  USE mo_impl_constants,      ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_physical_constants,  ONLY: rd, cvv, cvd, cpv, cpd, cvd_o_rd, clw, ci
  USE mo_parallel_config,     ONLY: nproma

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: thermo_src_term

  CONTAINS


  !>
  !! Apply correction term for moisture-dependence of heat capacities
  !!
  SUBROUTINE thermo_src_term(p_patch, p_int, p_nh, prep_adv, dtime, nnow_rcf, nnew)

    TYPE(t_patch), TARGET, INTENT(INOUT)    :: p_patch
    TYPE(t_nh_state), TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state), TARGET,INTENT(IN)    :: p_int
    TYPE(t_prepare_adv), INTENT(INOUT)      :: prep_adv

    REAL(wp), INTENT(IN) :: dtime
    INTEGER , INTENT(IN) :: nnow_rcf, nnew

    INTEGER :: i_startidx, i_endidx, jc, jk, jt, jb, nlev
    INTEGER :: i_rlstart , i_rlend , i_startblk, i_endblk
    INTEGER :: jg

    REAL(wp), POINTER :: p_csum(:,:)
    INTEGER, POINTER, CONTIGUOUS :: ieidx(:,:,:), ieblk(:,:,:)
    TYPE(t_nh_prog), POINTER :: p_prog_rcf, p_prog

    REAL(wp), DIMENSION(nproma,p_patch%nlev), TARGET :: qsum_liq, qsum_ice, chi_q, v_flxdiv
    REAL(wp) :: z_a, z_b, exner_sv, beta_q


    jg = p_patch%id

    ! number of vertical levels
    nlev = p_patch%nlev

    ! Set pointers to neighbor edges
    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    p_prog_rcf => p_nh%prog(nnow_rcf)
    p_prog     => p_nh%prog(nnew)

    !$ACC DATA CREATE(qsum_liq, qsum_ice, chi_q, v_flxdiv) COPYIN(kstart_moist)

!$OMP PARALLEL PRIVATE(i_rlstart, i_rlend, i_startblk, i_endblk)

    ! boundary zone and halo points are not needed
    i_rlstart  = grf_bdywidth_c + 1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jc,jk,jt,jb,i_startidx,i_endidx,p_csum,z_a,z_b,qsum_liq,qsum_ice,chi_q,v_flxdiv, &
!$OMP            exner_sv,beta_q) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
      &                   i_startidx, i_endidx, i_rlstart, i_rlend )


      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
!$NEC outerloop_unroll(8)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          qsum_liq(jc,jk) = 0.0_wp
          qsum_ice(jc,jk) = 0.0_wp

          v_flxdiv(jc,jk) = p_nh%metrics%deepatmo_divh_mc(jk) * (                              &
            p_nh%metrics%ddqz_z_full_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))* &
              prep_adv%vn_traj(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) + &
            p_nh%metrics%ddqz_z_full_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))* &
              prep_adv%vn_traj(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) + &
            p_nh%metrics%ddqz_z_full_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))* &
              prep_adv%vn_traj(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb)   )
        END DO
      END DO
      !$ACC END PARALLEL

      DO jt = 1, ntracer

        IF (ANY((/iqc,iqr,iqhl,iqgl/)==jt)) THEN
          p_csum => qsum_liq
        ELSE IF (ANY((/iqi,iqs,iqg,iqh/)==jt)) THEN
          p_csum => qsum_ice
        ELSE
          CYCLE
        END IF

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
         DO jk = kstart_moist(jg),nlev
           DO jc = i_startidx, i_endidx

             p_csum(jc,jk) = p_csum(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,jt)

           END DO ! jc
         END DO ! jk
        !$ACC END PARALLEL

      END DO ! jt

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_a, z_b)
      DO jk = kstart_moist(jg),nlev
        DO jc = i_startidx, i_endidx

          z_a = 1._wp - ( p_prog_rcf%tracer(jc,jk,jb,iqv) + qsum_liq(jc,jk) + qsum_ice(jc,jk) )

          z_b = clw * qsum_liq(jc,jk) + ci * qsum_ice(jc,jk)

          chi_q(jc,jk) = 1._wp - (                                             &
              ( z_a + ( z_b  + cpv * p_prog_rcf%tracer(jc,jk,jb,iqv) ) / cpd ) &
            / ( z_a + ( z_b  + cvv * p_prog_rcf%tracer(jc,jk,jb,iqv) ) / cvd ) )

        END DO ! jc
      END DO ! jk
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_a)
      DO jk = 2, kstart_moist(jg)-1
        DO jc = i_startidx, i_endidx

          z_a = 1._wp - p_prog_rcf%tracer(jc,jk,jb,iqv) 

          chi_q(jc,jk) = 1._wp - (                                   &
              ( z_a + cpv * p_prog_rcf%tracer(jc,jk,jb,iqv)  / cpd ) &
            / ( z_a + cvv * p_prog_rcf%tracer(jc,jk,jb,iqv)  / cvd ) )

        END DO ! jc
      END DO ! jk
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(beta_q, exner_sv)
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          exner_sv = p_prog%exner(jc,jk,jb)

          beta_q = dtime*rd*chi_q(jc,jk)*exner_sv*p_nh%metrics%inv_ddqz_z_full(jc,jk,jb) / cvd

          p_prog%exner(jc,jk,jb) = exner_sv + beta_q * ( v_flxdiv(jc,jk)           &
            + prep_adv%vol_flx_ic(jc,jk,jb)   * p_nh%metrics%deepatmo_divzU_mc(jk) &
            - prep_adv%vol_flx_ic(jc,jk+1,jb) * p_nh%metrics%deepatmo_divzL_mc(jk) ) 

          ! recompute theta
          p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) * &
            ( (p_prog%exner(jc,jk,jb)/exner_sv-1.0_wp) * cvd_o_rd+1.0_wp )
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    END DO ! jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE thermo_src_term

END MODULE mo_nh_moist_thdyn

