!
! Subroutine cloud_mig calls the saturation adjustment and cloud scheme
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

MODULE mo_cloud_mig

  USE mo_kind                ,ONLY: wp

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, &
       &                            timer_sat, timer_grp

  USE mo_aes_thermo          ,ONLY: saturation_adjustment
  USE mo_aes_graupel         ,ONLY: graupel_run

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cloud_mig

  REAL(wp), PARAMETER ::  cloud_num =  50.00e+06_wp   ! cloud droplet number concentration

CONTAINS

  SUBROUTINE cloud_mig     ( jcs, jce   ,&
       &                     pdtime     ,&
       &                     dz         ,&
       &                     rho        ,&
       &                     pf         ,&
       &                     ta         ,&
       &                     qv         ,&
       &                     qc         ,&
       &                     qi         ,&
       &                     qr         ,&
       &                     qs         ,&
       &                     qg         ,&
       &                     tend_ta    ,&
       &                     tend_qv    ,&
       &                     tend_qc    ,&
       &                     tend_qi    ,&
       &                     tend_qr    ,&
       &                     tend_qs    ,&
       &                     tend_qg    ,&
       &                     pr_rain    ,&
       &                     pr_ice     ,&
       &                     pr_snow    ,&
       &                     pr_grpl    ,&
       &                     pr_eflx    )

    ! Arguments
    !
    INTEGER , INTENT(in)  :: jcs, jce      !< column index range
    REAL(wp), INTENT(in)  :: pdtime        !< timestep
    !
    REAL(wp), INTENT(in)  :: dz      (:,:) !< vertical layer thickness
    REAL(wp), INTENT(in)  :: rho     (:,:) !< density
    REAL(wp), INTENT(in)  :: pf      (:,:) !< pressure
    !
    REAL(wp), INTENT(in)  :: ta      (:,:) !< temperature
    REAL(wp), INTENT(in)  :: qv      (:,:) !< sp humidity
    REAL(wp), INTENT(in)  :: qc      (:,:) !< cloud water
    REAL(wp), INTENT(in)  :: qi      (:,:) !< ice
    REAL(wp), INTENT(in)  :: qr      (:,:) !< rain
    REAL(wp), INTENT(in)  :: qs      (:,:) !< snow
    REAL(wp), INTENT(in)  :: qg      (:,:) !< graupel
    !
    REAL(wp), INTENT(out) :: tend_ta (:,:) !< tendency of temperature
    REAL(wp), INTENT(out) :: tend_qv (:,:) !< tendency of water vapor
    REAL(wp), INTENT(out) :: tend_qc (:,:) !< tendency of cloud water
    REAL(wp), INTENT(out) :: tend_qi (:,:) !< tendency of cloud ice
    REAL(wp), INTENT(out) :: tend_qr (:,:) !< tendency of rain
    REAL(wp), INTENT(out) :: tend_qs (:,:) !< tendency of snow
    REAL(wp), INTENT(out) :: tend_qg (:,:) !< tendency of graupel

    REAL(wp), INTENT(out) :: pr_rain (:)   !< precip rate rain
    REAL(wp), INTENT(out) :: pr_ice  (:)   !< precip rate ice
    REAL(wp), INTENT(out) :: pr_snow (:)   !< precip rate snow
    REAL(wp), INTENT(out) :: pr_grpl (:)   !< precip rate graupel
    REAL(wp), INTENT(out) :: pr_eflx (:)   !< energy flux from precip

    ! Local variables
    !
    INTEGER  :: nproma
    INTEGER  :: jc
    INTEGER  :: jk, jks, jke
    !
    REAL(wp) :: zqnc(SIZE(dz,1))
    !
    REAL(wp) :: zta(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqv(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqc(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqi(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqr(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqs(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqg(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: total_ice(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqrsflux (SIZE(dz,1),SIZE(dz,2))
    !
    REAL(wp) :: zdtr ! reciprocal of timestep
    REAL(wp), PARAMETER :: cloud_num = 200.00e+06_wp

    !$ACC DATA PRESENT(dz, rho, pf, ta, qv, qc, qi, qr, qs, qg) &
    !$ACC   PRESENT(tend_ta, tend_qv, tend_qc, tend_qi, tend_qr, tend_qs, tend_qg) &
    !$ACC   PRESENT(pr_ice, pr_rain, pr_snow, pr_grpl, pr_eflx) &
    !$ACC   CREATE(zqnc, zta, zqv, zqc, zqi, zqr, zqs, zqg, total_ice, zqrsflux)

    nproma = SIZE(dz,1)

    jks = 1
    jke = (SIZE(dz,2))

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,jce
       zqnc(jc) = cloud_num
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = jks,jke
       DO jc = jcs,jce
          zta(jc,jk) = ta(jc,jk)
          zqv(jc,jk) = qv(jc,jk)
          zqc(jc,jk) = qc(jc,jk)
          zqi(jc,jk) = qi(jc,jk)
          zqr(jc,jk) = qr(jc,jk)
          zqs(jc,jk) = qs(jc,jk)
          zqg(jc,jk) = qg(jc,jk)
          total_ice(jc,jk) = qg(jc,jk)+qs(jc,jk)+qi(jc,jk)
       END DO
    END DO
    !$ACC END PARALLEL

    zdtr = 1._wp/pdtime

    ! Initial saturation adjustment
    !
    IF (ltimer) call timer_start(timer_sat)
    !
    CALL saturation_adjustment(                  &
         &           ilo      = jcs             ,& !> in
         &           iup      = jce             ,& !> in
         &           klo      = jks             ,& !> in
         &           kup      = jke             ,& !> in
         &           te       = zta       (:,:) ,& !> inout
         &           qve      = zqv       (:,:) ,& !> inout
         &           qce      = zqc       (:,:) ,& !> inout
         &           qre      = zqr       (:,:) ,& !> in
         &           qti      = total_ice (:,:) ,& !> in
         &           rho      = rho       (:,:) )  !> in

    IF (ltimer) call timer_stop(timer_sat)

    ! Single moment cloud microphyiscs for water vapor,
    ! cloud water, cloud ice, rain, snow and graupel
    !
    IF (ltimer) call timer_start(timer_grp)
    !
    CALL graupel_run( nvec= nproma        ,& !< in
         &        ke      = jke           ,& !< in
         &        ivstart = jcs           ,& !< in
         &        ivend   = jce           ,& !< in
         &        kstart  = jks           ,& !< in
         &        dt      = pdtime        ,& !< in   : timestep
         &        qnc     = zqnc    (:)   ,& !< in
         &        dz      = dz      (:,:) ,& !< in   : vertical layer thickness
         &        rho     = rho     (:,:) ,& !< in   : density
         &        p       = pf      (:,:) ,& !< in   : pressure
         !
         &        t       = zta     (:,:) ,& !< inout: temp
         &        qv      = zqv     (:,:) ,& !< inout: sp humidity
         &        qc      = zqc     (:,:) ,& !< inout: cloud water
         &        qi      = zqi     (:,:) ,& !< inout: ice
         &        qr      = zqr     (:,:) ,& !< inout: rain
         &        qs      = zqs     (:,:) ,& !< inout: snow
         &        qg      = zqg     (:,:) ,& !< inout: graupel

         &        pflx    = zqrsflux(:,:) ,& !<   out: precip flux in atmosphere
         &        prr_gsp = pr_rain (:)   ,& !<   out: precip rate rain
         &        pri_gsp = pr_ice  (:)   ,& !<   out: precip rate cloud ice
         &        prs_gsp = pr_snow (:)   ,& !<   out: precip rate snow
         &        prg_gsp = pr_grpl (:)   ,& !<   out: precip rate graupel
         &        pre_gsp = pr_eflx (:)   )  !<   out: precip energy flux
    !
    IF (ltimer) call timer_stop(timer_grp)

    ! Final saturation adjustment
    !
    IF (ltimer) call timer_start(timer_sat)
    !
    CALL saturation_adjustment(                  &
         &           ilo      = jcs             ,& !> in
         &           iup      = jce             ,& !> in
         &           klo      = jks             ,& !> in
         &           kup      = jke             ,& !> in
         &           te       = zta       (:,:) ,& !> inout
         &           qve      = zqv       (:,:) ,& !> inout
         &           qce      = zqc       (:,:) ,& !> inout
         &           qre      = zqr       (:,:) ,& !> in
         &           qti      = total_ice (:,:) ,& !> in
         &           rho      = rho       (:,:) )  !> in
    !
    IF (ltimer) call timer_stop(timer_sat)

    ! Calculate tendencies and convert temperature tendency, as computed
    ! in satad/graupel for constant volume to constant pressure
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = jks,jke
       DO jc = jcs,jce
          tend_ta(jc,jk) =     (zta(jc,jk)-ta(jc,jk))*zdtr
          tend_qv(jc,jk) = MAX((zqv(jc,jk)-qv(jc,jk))*zdtr,-qv(jc,jk)*zdtr)
          tend_qc(jc,jk) = MAX((zqc(jc,jk)-qc(jc,jk))*zdtr,-qc(jc,jk)*zdtr)
          tend_qi(jc,jk) = MAX((zqi(jc,jk)-qi(jc,jk))*zdtr,-qi(jc,jk)*zdtr)
          tend_qr(jc,jk) = MAX((zqr(jc,jk)-qr(jc,jk))*zdtr,-qr(jc,jk)*zdtr)
          tend_qs(jc,jk) = MAX((zqs(jc,jk)-qs(jc,jk))*zdtr,-qs(jc,jk)*zdtr)
          tend_qg(jc,jk) = MAX((zqg(jc,jk)-qg(jc,jk))*zdtr,-qg(jc,jk)*zdtr)
       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE cloud_mig

  !-------------------------
END MODULE mo_cloud_mig
