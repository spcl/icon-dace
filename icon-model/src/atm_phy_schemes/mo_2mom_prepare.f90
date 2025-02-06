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

MODULE mo_2mom_prepare

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_2mom_mcrph_main, ONLY: particle, particle_lwf, atmosphere
#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PUBLIC :: prepare_twomoment, post_twomoment

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_prepare'

CONTAINS

  SUBROUTINE prepare_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
       rho, rhocorr, rhocld, pres, w, tk, hhl, tke, &
       nccn, ninpot, ninact, &
       qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl, &
       lprogccn, lprogin, lprogmelt, its, ite, kts, kte)

    TYPE(atmosphere), INTENT(inout)   :: atmo
    CLASS(particle),  INTENT(inout)   :: cloud, rain, ice, snow
    CLASS(particle),  INTENT(inout)   :: graupel, hail
    REAL(wp), TARGET, DIMENSION(:, :), INTENT(in) :: &
         rho, rhocorr, rhocld, pres, w, tk, hhl
    REAL(wp), POINTER, DIMENSION(:, :), INTENT(in) :: tke
    REAL(wp), DIMENSION(:,:), INTENT(inout) , TARGET :: &
         &               qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh
    LOGICAL, INTENT(in) :: lprogccn, lprogin, lprogmelt
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               nccn, ninpot, ninact
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               qgl,qhl
    INTEGER, INTENT(in) :: its, ite, kts, kte
    INTEGER :: ii, kk

    ! ... Transformation of microphysics variables to densities
    !$ACC DATA PRESENT(qv, qc, qr, qi, qs, qg, qh, qnc, qnr, qni, qns, qng, qnh) &
    !$ACC   PRESENT(ninact, ninpot, nccn, rho, cloud) &
    !$ACC   PRESENT(rain, ice, snow, graupel, hail)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(kts, kte, its, ite, lprogccn, lprogin)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO kk = kts, kte
      DO ii = its, ite

        ! ... concentrations --> number densities
        qnc(ii,kk) = rho(ii,kk) * qnc(ii,kk)
        qnr(ii,kk) = rho(ii,kk) * qnr(ii,kk)
        qni(ii,kk) = rho(ii,kk) * qni(ii,kk)
        qns(ii,kk) = rho(ii,kk) * qns(ii,kk)
        qng(ii,kk) = rho(ii,kk) * qng(ii,kk)
        qnh(ii,kk) = rho(ii,kk) * qnh(ii,kk)

        ! ... mixing ratios -> mass densities
        qv(ii,kk) = rho(ii,kk) * qv(ii,kk)
        qc(ii,kk) = rho(ii,kk) * qc(ii,kk)
        qr(ii,kk) = rho(ii,kk) * qr(ii,kk)
        qi(ii,kk) = rho(ii,kk) * qi(ii,kk)
        qs(ii,kk) = rho(ii,kk) * qs(ii,kk)
        qg(ii,kk) = rho(ii,kk) * qg(ii,kk)
        qh(ii,kk) = rho(ii,kk) * qh(ii,kk)

        ninact(ii,kk)  = rho(ii,kk) * ninact(ii,kk)

        IF (lprogccn) THEN
          nccn(ii,kk) = rho(ii,kk) * nccn(ii,kk)
        END IF
        IF (lprogin) THEN
          ninpot(ii,kk)  = rho(ii,kk) * ninpot(ii,kk)
        END IF
#ifndef _OPENACC
        IF (lprogmelt) THEN
          qgl(ii,kk)  = rho(ii,kk) * qgl(ii,kk)
          qhl(ii,kk)  = rho(ii,kk) * qhl(ii,kk)
        END IF
#endif

      END DO
    END DO
    !$ACC END PARALLEL

    IF (lprogmelt.AND.(.not.PRESENT(qgl).or..not.PRESENT(qhl))) THEN
      CALL finish(TRIM(routine),'Error in prepare_twomoment, something wrong with qgl or qhl')
    END IF

    ! set pointers
    atmo%w   => w
    atmo%T   => tk
    atmo%p   => pres
    atmo%qv  => qv
    atmo%rho => rho
    atmo%zh  => hhl

    IF (ASSOCIATED(tke)) THEN
      atmo%tke => tke
    ELSE
      atmo%tke=>NULL()
    END IF

    __acc_attach(atmo%w)
    __acc_attach(atmo%T)
    __acc_attach(atmo%p)
    __acc_attach(atmo%qv)
    __acc_attach(atmo%rho)
    __acc_attach(atmo%zh)
    __acc_attach(atmo%tke)

    cloud%rho_v   => rhocld
    rain%rho_v    => rhocorr
    ice%rho_v     => rhocorr
    graupel%rho_v => rhocorr
    snow%rho_v    => rhocorr
    hail%rho_v    => rhocorr

    cloud%q   => qc
    cloud%n   => qnc
    rain%q    => qr
    rain%n    => qnr
    ice%q     => qi
    ice%n     => qni
    snow%q    => qs
    snow%n    => qns
    graupel%q => qg
    graupel%n => qng
    hail%q    => qh
    hail%n    => qnh

    SELECT TYPE (graupel)
    CLASS IS (particle_lwf)
       graupel%l => qgl
    END SELECT

    SELECT TYPE (hail)
    CLASS IS (particle_lwf)
       hail%l    => qhl
    END SELECT

    __acc_attach(ice%rho_v)
    __acc_attach(graupel%rho_v)
    __acc_attach(snow%rho_v)
    __acc_attach(hail%rho_v)
    __acc_attach(cloud%rho_v)
    __acc_attach(rain%rho_v)
    __acc_attach(ice%q)
    __acc_attach(ice%n)
    __acc_attach(snow%q)
    __acc_attach(snow%n)
    __acc_attach(graupel%q)
    __acc_attach(graupel%n)
    __acc_attach(hail%q)
    __acc_attach(hail%n)
    __acc_attach(cloud%q)
    __acc_attach(cloud%n)
    __acc_attach(rain%q)
    __acc_attach(rain%n)

    ! enforce upper and lower bounds for number concentrations
    ! (may not be necessary or only at initial time)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(kts, kte, its, ite)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO kk=kts,kte
      DO ii=its,ite
        rain%n(ii,kk) = MIN(rain%n(ii,kk), rain%q(ii,kk)/rain%x_min)
        rain%n(ii,kk) = MAX(rain%n(ii,kk), rain%q(ii,kk)/rain%x_max)
        ice%n(ii,kk) = MIN(ice%n(ii,kk), ice%q(ii,kk)/ice%x_min)
        ice%n(ii,kk) = MAX(ice%n(ii,kk), ice%q(ii,kk)/ice%x_max)
        snow%n(ii,kk) = MIN(snow%n(ii,kk), snow%q(ii,kk)/snow%x_min)
        snow%n(ii,kk) = MAX(snow%n(ii,kk), snow%q(ii,kk)/snow%x_max)
        graupel%n(ii,kk) = MIN(graupel%n(ii,kk), graupel%q(ii,kk)/graupel%x_min)
        graupel%n(ii,kk) = MAX(graupel%n(ii,kk), graupel%q(ii,kk)/graupel%x_max)
        hail%n(ii,kk) = MIN(hail%n(ii,kk), hail%q(ii,kk)/hail%x_min)
        hail%n(ii,kk) = MAX(hail%n(ii,kk), hail%q(ii,kk)/hail%x_max)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(kts, kte, its, ite)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO kk=kts,kte
      DO ii=its,ite
        IF(cloud%q(ii,kk) <= 1.0e-12) cloud%n(ii,kk) = 0.0_wp
        IF(rain%q(ii,kk) <= 1.0e-12) rain%n(ii,kk) = 0.0_wp
        IF(ice%q(ii,kk) <= 1.0e-12) ice%n(ii,kk) = 0.0_wp
        IF(snow%q(ii,kk) <= 1.0e-12) snow%n(ii,kk) = 0.0_wp
        IF(graupel%q(ii,kk) <= 1.0e-12) graupel%n(ii,kk) = 0.0_wp
        IF(hail%q(ii,kk) <= 1.0e-12) hail%n(ii,kk) = 0.0_wp
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA ! DATA PRESENT

  END SUBROUTINE prepare_twomoment

  SUBROUTINE post_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
       rho_r, qnc, nccn, ninpot, ninact, &
       qv, qc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl,  &
       lprogccn, lprogin, lprogmelt, its, ite, kts, kte)

    TYPE(atmosphere), INTENT(inout)   :: atmo
    CLASS(particle), INTENT(inout)    :: cloud, rain, ice, snow
    CLASS(particle), INTENT(inout)    :: graupel, hail
    REAL(wp), INTENT(in) :: rho_r(:, :)
    REAL(wp), DIMENSION(:,:), INTENT(inout) :: &
         &           qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &           nccn, ninpot, ninact
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               qgl,qhl
    LOGICAL, INTENT(in) :: lprogccn, lprogin, lprogmelt
    INTEGER, INTENT(in) :: its, ite, kts, kte
    INTEGER :: ii, kk
    REAL(wp) :: hlp

    IF (lprogmelt.AND.(.not.PRESENT(qgl).or..not.PRESENT(qhl))) THEN
      CALL finish(TRIM(routine),'Error in post_twomoment, something wrong with qgl or qhl')
    END IF

    ! nullify pointers
    atmo%w   => NULL()
    atmo%T   => NULL()
    atmo%p   => NULL()
    atmo%qv  => NULL()
    atmo%rho => NULL()
    atmo%zh  => NULL()
    atmo%tke => NULL()
    
    cloud%rho_v   => NULL()
    rain%rho_v    => NULL()
    ice%rho_v     => NULL()
    graupel%rho_v => NULL()
    snow%rho_v    => NULL()
    hail%rho_v    => NULL()
    
    cloud%q   => NULL()
    cloud%n   => NULL()
    rain%q    => NULL()
    rain%n    => NULL()
    ice%q     => NULL()
    ice%n     => NULL()
    snow%q    => NULL()
    snow%n    => NULL()
    graupel%q => NULL()
    graupel%n => NULL()
    hail%q    => NULL()
    hail%n    => NULL()

    SELECT TYPE (graupel)
    CLASS IS (particle_lwf) 
      graupel%l => NULL()
    END SELECT

    SELECT TYPE (hail)
    CLASS IS (particle_lwf) 
      hail%l    => NULL()
    END SELECT

    ! ... Transformation of variables back to ICON standard variables
    !$ACC DATA PRESENT(qv, qc, qr, qi, qs, qg, qh, qnc, qnr, qni, qns, qng, qnh) &
    !$ACC   PRESENT(ninact, ninpot, nccn, rho_r)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(kts, kte, its, ite, lprogccn, lprogin)
    !$ACC LOOP GANG VECTOR PRIVATE(hlp) COLLAPSE(2)
    DO kk = kts, kte
      DO ii = its, ite

        hlp = rho_r(ii,kk)

        ! ... from mass densities back to mixing ratios
        qv(ii,kk) = hlp * qv(ii,kk)
        qc(ii,kk) = hlp * qc(ii,kk)
        qr(ii,kk) = hlp * qr(ii,kk)
        qi(ii,kk) = hlp * qi(ii,kk)
        qs(ii,kk) = hlp * qs(ii,kk)
        qg(ii,kk) = hlp * qg(ii,kk)
        qh(ii,kk) = hlp * qh(ii,kk)

        ! ... number concentrations
        qnc(ii,kk) = hlp * qnc(ii,kk)
        qnr(ii,kk) = hlp * qnr(ii,kk)
        qni(ii,kk) = hlp * qni(ii,kk)
        qns(ii,kk) = hlp * qns(ii,kk)
        qng(ii,kk) = hlp * qng(ii,kk)
        qnh(ii,kk) = hlp * qnh(ii,kk)

        ninact(ii,kk)  = hlp * ninact(ii,kk)

        IF (lprogccn) THEN
          nccn(ii,kk) = hlp * nccn(ii,kk)
        END IF
        IF (lprogin) THEN
          ninpot(ii,kk)  = hlp * ninpot(ii,kk)
        END IF
#ifndef _OPENACC
        IF (lprogmelt) THEN
          qgl(ii,kk)  = hlp * qgl(ii,kk)
          qhl(ii,kk)  = hlp * qhl(ii,kk)
        END IF
#endif

      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE post_twomoment

END MODULE mo_2mom_prepare
