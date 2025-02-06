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

! Module containing subroutine computing the vertical integral of the internal energy

MODULE mo_diagnose_uvi

  USE mo_kind            ,ONLY: wp

  USE mo_run_config      ,ONLY: num_lev, iqv, iqc, iqi, iqr, iqs, iqg
  USE mo_dynamics_config ,ONLY: nnew, nnew_rcf

  USE mo_nonhydro_state  ,ONLY: p_nh_state
  USE mo_aes_phy_memory  ,ONLY: prm_field

  USE mo_timer           ,ONLY: ltimer, timer_start, timer_stop, timer_uvi

  USE mo_aes_thermo      ,ONLY: internal_energy

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: diagnose_uvd, diagnose_uvp

CONTAINS

  !-------------------------------------------------------------------

  SUBROUTINE diagnose_uvd(jg, jb, jcs, jce)

    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    CALL u_vertical_integral  (jg, jb, jcs, jce, prm_field(jg)%udynvi(:,jb))

  END SUBROUTINE diagnose_uvd

  !-------------------------------------------------------------------

  SUBROUTINE diagnose_uvp(jg, jb, jcs, jce)

    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    REAL(wp) :: uphyvi(SIZE(prm_field(jg)%udynvi,1))
    INTEGER :: jc

    !$ACC DATA CREATE(uphyvi)

    CALL u_vertical_integral  (jg, jb, jcs, jce, uphyvi(:))

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jc = jcs, jce
      prm_field(jg)%duphyvi(jc,jb) = uphyvi(jc) - prm_field(jg)%udynvi(jc,jb)
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE diagnose_uvp

  !-------------------------------------------------------------------

  SUBROUTINE u_vertical_integral(jg, jb, jcs, jce, uvi)

    INTEGER , INTENT(in)    :: jg, jb, jcs, jce
    REAL(wp), INTENT(out)   :: uvi(:)

    INTEGER                 :: jc, jk, jtl_dyn, jtl_trc
    REAL(wp)                :: qliquid, qfrozen

    IF (ltimer) call timer_start(timer_uvi)

    jtl_dyn = nnew(jg)
    jtl_trc = nnew_rcf(jg)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = jcs, jce
      uvi(jc) = 0.0_wp
    END DO

    !$ACC LOOP SEQ
    DO jk = 1,num_lev(jg)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = jcs, jce

        qliquid =   p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqc) &
                & + p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqr) 

        qfrozen =   p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqi) &
                & + p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqs) &
                & + p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqg)

        uvi(jc) =   uvi(jc)                                                            &
             &    + internal_energy(p_nh_state(jg)%diag%temp(jc,jk,jb),                & ! temperature
             &                      p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqv), & ! qv
             &                      qliquid,                                           & ! sum of liquid phases
             &                      qfrozen,                                           & ! sum of frozen phases
             &                      p_nh_state(jg)%prog(jtl_dyn)%rho(jc,jk,jb),        & ! air density
             &                      p_nh_state(jg)%metrics%ddqz_z_full(jc,jk,jb)     )   ! layer thickness

      END DO !jc
    END DO !jk
    !$ACC END PARALLEL

    IF (ltimer) call timer_stop(timer_uvi)

  END SUBROUTINE  u_vertical_integral

  !-------------------------------------------------------------------

END MODULE mo_diagnose_uvi
