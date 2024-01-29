!
! Module containing subroutine computing the vertical integral of tracer mass
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

MODULE mo_diagnose_qvi

  USE mo_kind            ,ONLY: wp

  USE mo_run_config      ,ONLY: ntracer, num_lev
  USE mo_dynamics_config ,ONLY: nnew_rcf

  USE mo_nonhydro_state  ,ONLY: p_nh_state
  USE mo_aes_phy_memory  ,ONLY: prm_field, prm_tend

  USE mo_timer           ,ONLY: ltimer, timer_start, timer_stop, timer_qvi

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: diagnose_qvi

CONTAINS

  !-------------------------------------------------------------------

  SUBROUTINE diagnose_qvi(jg, jb, jcs, jce)

    INTEGER , INTENT(in)    :: jg, jb, jcs, jce

    INTEGER                 :: jc, jk, jt, jtl_trc

    IF (ltimer) call timer_start(timer_qvi)

    jtl_trc = nnew_rcf(jg)

    DO jt = 1,ntracer
      !
      IF (ASSOCIATED(prm_field(jg)%mtrcvi_ptr(jt)%p)) THEN
        ! 
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = jcs, jce
          prm_field(jg)%mtrcvi(jc,jb,jt) = 0.0_wp
        END DO ! jc
        !
        !$ACC LOOP SEQ
        DO jk = 1,num_lev(jg)
        !$ACC LOOP GANG VECTOR
          DO jc = jcs, jce
            !
            ! tracer path
            prm_field(jg)%mtrcvi(jc,jb,jt) =   prm_field(jg)%mtrcvi(jc,jb,jt) &
                  & + p_nh_state(jg)%diag%airmass_new(jc,jk,jb)*p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,jt)
            !
          END DO ! jc
        END DO ! jk
        !$ACC END PARALLEL
        !
      END IF
      !
      IF (ASSOCIATED(prm_tend(jg)%mtrcvi_phy_ptr(jt)%p)) THEN
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = jcs, jce
          prm_tend(jg)%mtrcvi_phy(jc,jb,jt) = 0.0_wp
        END DO ! jc
        !
        !$ACC LOOP SEQ
        DO jk = 1,num_lev(jg)
          !$ACC LOOP GANG VECTOR
          DO jc = jcs, jce
            !
            ! tendency of tracer path
            prm_tend(jg)%mtrcvi_phy(jc,jb,jt) = prm_tend(jg)%mtrcvi_phy(jc,jb,jt) &
                  & + p_nh_state(jg)%diag%airmass_new(jc,jk,jb)*prm_tend(jg)%qtrc_phy(jc,jk,jb,jt)
            !
          END DO ! jc
        END DO ! jk
        !$ACC END PARALLEL
        !
      END IF
      !
    END DO ! jt

    IF (ltimer) call timer_stop(timer_qvi)

  END SUBROUTINE  diagnose_qvi

  !-------------------------------------------------------------------

END MODULE mo_diagnose_qvi
