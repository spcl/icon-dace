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

! Subroutine cloud_two calls the saturation adjustment
! and the two-moment bulk microphysics by Seifert and Beheng (2006)
!                  with prognostic cloud droplet number

MODULE mo_cloud_two

  USE mo_kind                ,ONLY: wp

!  USE mo_cloud_two_config    ,ONLY: cloud_two_config
  USE mo_physical_constants  ,ONLY: cvd
  USE mo_satad               ,ONLY: satad_v_3d
  USE mo_2mom_mcrph_driver   ,ONLY: two_moment_mcrph

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cloud_two

CONTAINS

  SUBROUTINE cloud_two     ( nproma, nlev  ,&
       &                     jcs, jce      ,&
       &                     jks, jke      ,&
       &                     msg_level     ,&
       &                     pdtime        ,& ! in: time step
       &                     dz            ,& ! in: vertical layer thickness
       &                     zh            ,& ! in: height of half levels
       &                     rho           ,& ! in: density
       &                     pf            ,& ! in: pressure
       &                     tke           ,& ! in: TKE (on half levels, size nlev+1)
       &                     ta            ,& ! inout: temp
       &                     qv            ,& ! inout: specific humidity
       &                     qc, qnc       ,& ! inout: cloud water
       &                     qi, qni       ,& ! inout: ice
       &                     qr, qnr       ,& ! inout: rain
       &                     qs, qns       ,& ! inout: snow
       &                     qg, qng       ,& ! inout: graupel
       &                     qh, qnh       ,& ! inout: hail
       &                     ninact        ,& ! inout: activated ice nuclei
       &                     wa            ,& ! inout: w (on half levels, size nlev+1)
       &                     tend_ta       ,& ! out: tendency of temperature
       &                     tend_qv       ,& ! out: tendency of specific humidity
       &                     tend_qc       ,& ! out: tendency of cloud water
       &                     tend_qnc      ,& ! out: tendency of cloud water droplets
       &                     tend_qi       ,& ! out: tendency of cloud ice
       &                     tend_qni      ,& ! out: tendency of cloud ice droplets
       &                     tend_qr       ,& ! out: tendency of rain
       &                     tend_qnr      ,& ! out: tendency of rain droplets
       &                     tend_qs       ,& ! out: tendency of snow
       &                     tend_qns      ,& ! out: tendency of snow droplets
       &                     tend_qg       ,& ! out: tendency of graupel
       &                     tend_qng      ,& ! out: tendency of graupel droplets
       &                     tend_qh       ,& ! out: tendency of hail
       &                     tend_qnh      ,& ! out: tendency of hail droplets
       &                     tend_ninact   ,& ! out: tendency of activated ice nuclei
       &                     pr_rain       ,& ! inout: precip rate rain
       &                     pr_ice        ,& ! inout: precip rate ice
       &                     pr_snow       ,& ! inout: precip rate snow
       &                     pr_grpl       ,& ! inout: precip rate graupel
       &                     pr_hail       ,& ! inout: precip rate hail
       &                     zqrsflux      )  ! inout: unused 3d total prec. rate

    ! Arguments
    !
    INTEGER , INTENT(in)  :: nproma        !< grid index
    INTEGER , INTENT(in)  :: nlev          !< grid index
    INTEGER , INTENT(in)  :: jcs, jce      !< column index range
    INTEGER , INTENT(in)  :: jks, jke      !< column index range
    INTEGER , INTENT(in)  :: msg_level     !< message level
    REAL(wp), INTENT(in)  :: pdtime        !< timestep
    !
    REAL(wp), DIMENSION(:,:), INTENT(in)  :: dz       !< vertical layer thickness
    REAL(wp), DIMENSION(:,:), INTENT(in)  :: zh       !< height of half levels
    REAL(wp), DIMENSION(:,:), INTENT(in)  :: rho      !< density
    REAL(wp), DIMENSION(:,:), INTENT(in)  :: pf       !< pressure
    REAL(wp), DIMENSION(:,:), INTENT(in), POINTER  :: tke      !< TKE
    !
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: ta       !< temperature
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qv       !< sp humidity
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qc       !< cloud water
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qnc      !< cloud water number
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qi       !< ice
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qni      !< ice number
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qr       !< rain
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qnr      !< rain number
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qs       !< snow
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qns      !< snow number
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qg       !< graupel
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qng      !< graupel number
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qh       !< hail
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: qnh      !< hail number
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: ninact   !< activated ice nuclei
    REAL(wp), DIMENSION(:,:), INTENT(inout)  :: wa       !< vertical velocity
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_ta      !< tendency of temperature
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qv      !< tendency of water vapor
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qc      !< tendency of cloud water
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qi      !< tendency of cloud ice
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qr      !< tendency of rain
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qs      !< tendency of snow
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qg      !< tendency of graupel
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qh      !< tendency of hail
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qnc     !< tendency of cloud water droplets
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qni     !< tendency of cloud ice droplets
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qnr     !< tendency of rain droplets
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qns     !< tendency of snow droplets
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qng     !< tendency of graupel droplets
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_qnh     !< tendency of hail droplets
    REAL(wp), DIMENSION(:,:), INTENT(out) :: tend_ninact  !< tendency of activated ice nuclei

    REAL(wp), DIMENSION(:), INTENT(inout) :: pr_rain    !< precip rate rain
    REAL(wp), DIMENSION(:), INTENT(inout) :: pr_ice     !< precip rate ice
    REAL(wp), DIMENSION(:), INTENT(inout) :: pr_snow    !< precip rate snow
    REAL(wp), DIMENSION(:), INTENT(inout) :: pr_grpl    !< precip rate graupel
    REAL(wp), DIMENSION(:), INTENT(inout) :: pr_hail    !< precip rate hail

    REAL(wp), DIMENSION(:,:), INTENT(inout) :: zqrsflux   !< precip rate hail

    ! Local variables
    !
    INTEGER  :: jc, jk
    !
    !
    REAL(wp) :: zta(nproma,nlev)
    REAL(wp) :: zqv(nproma,nlev)
    REAL(wp) :: zqc(nproma,nlev)
    REAL(wp) :: zqi(nproma,nlev)
    REAL(wp) :: zqr(nproma,nlev)
    REAL(wp) :: zqs(nproma,nlev)
    REAL(wp) :: zqg(nproma,nlev)
    REAL(wp) :: zqh(nproma,nlev)
    REAL(wp) :: zqnc(nproma,nlev)
    REAL(wp) :: zqni(nproma,nlev)
    REAL(wp) :: zqnr(nproma,nlev)
    REAL(wp) :: zqns(nproma,nlev)
    REAL(wp) :: zqng(nproma,nlev)
    REAL(wp) :: zqnh(nproma,nlev)
    REAL(wp) :: zninact(nproma,nlev)
    REAL(wp) :: zhhl(nproma,nlev+1)
    !
    REAL(wp) :: zdtr ! reciprocal of timestep


    zta(:,:)     = ta(:,:)
    zqv(:,:)     = qv(:,:)
    zqc(:,:)     = qc(:,:)
    zqi(:,:)     = qi(:,:)
    zqr(:,:)     = qr(:,:)
    zqs(:,:)     = qs(:,:)
    zqg(:,:)     = qg(:,:)
    zqh(:,:)     = qh(:,:)
    zqnc(:,:)    = qnc(:,:)
    zqni(:,:)    = qni(:,:)
    zqnr(:,:)    = qnr(:,:)
    zqns(:,:)    = qns(:,:)
    zqng(:,:)    = qng(:,:)
    zqnh(:,:)    = qnh(:,:)
    zninact(:,:) = ninact(:,:)
    zhhl  (:,:)  = zh  (:,:)

    zdtr = 1._wp/pdtime

    ! Initial saturation adjustment
    !
       CALL satad_v_3d( maxiter  = 10              ,& !> in
            &           idim     = nproma          ,& !> in
            &           kdim     = jke             ,& !> in
            &           ilo      = jcs             ,& !> in
            &           iup      = jce             ,& !> in
            &           klo      = jks             ,& !> in
            &           kup      = jke             ,& !> in
            &           tol      = 1.e-3_wp        ,& !> in
            &           te       = zta       (:,:) ,& !> inout
            &           qve      = zqv       (:,:) ,& !> inout
            &           qce      = zqc       (:,:) ,& !> inout
            &           rhotot   = rho       (:,:) )
          !
    ! Single moment cloud microphyiscs for water vapor,
    ! cloud water, cloud ice, rain, snow and graupel
    !
    !
    CALL two_moment_mcrph(                 &
         &        isize   = nproma        ,& !<    in: array size
         &        ke      = jke           ,& !<    in: end level/array size
         &        is      = jcs           ,& !<    in: start index
         &        ie      = jce           ,& !<    in: end index
         &        ks      = jks           ,& !<    in: start level
         &        dt      = pdtime        ,& !<    in: time step
         &        dz      = dz      (:,:) ,& !<    in: vertical layer thickness
         &        hhl     = zhhl    (:,:) ,& !<    in: height of half levels
         &        rho     = rho     (:,:) ,& !<    in: density
         &        pres    = pf      (:,:) ,& !<    in: pressure
         &        tke     = tke     (:,:) ,& !<    in: TKE (on half levels, size nlev+1)
         &        qv      = zqv     (:,:) ,& !< inout: sp humidity
         &        qc      = zqc     (:,:) ,& !< inout: cloud water
         &        qnc     = zqnc    (:,:) ,& !< inout: cloud droplet number 
         &        qr      = zqr     (:,:) ,& !< inout: rain
         &        qnr     = zqnr    (:,:) ,& !< inout: rain droplet number 
         &        qi      = zqi     (:,:) ,& !< inout: ice
         &        qni     = zqni    (:,:) ,& !< inout: cloud ice number
         &        qs      = zqs     (:,:) ,& !< inout: snow
         &        qns     = zqns    (:,:) ,& !< inout: snow number
         &        qg      = zqg     (:,:) ,& !< inout: graupel
         &        qng     = zqng    (:,:) ,& !< inout: graupel number
         &        qh      = zqh     (:,:) ,& !< inout: hail 
         &        qnh     = zqnh    (:,:) ,& !< inout: hail number
         &        ninact  = zninact (:,:) ,& !< inout: IN number
         &        tk      = zta     (:,:) ,& !< inout: temp
         &        w       = wa      (:,:) ,& !< inout: w (on half levels, size nlev+1)
         &        prec_r  = pr_rain (:)   ,& !< inout: precip rate rain
         &        prec_i  = pr_ice  (:)   ,& !< inout: precip rate ice
         &        prec_s  = pr_snow (:)   ,& !< inout: precip rate snow
         &        prec_g  = pr_grpl (:)   ,& !< inout: precip rate graupel
         &        prec_h  = pr_hail (:)   ,& !< inout: precip rate hail
         &        qrsflux = zqrsflux(:,:) ,& !< inout: 3D total precipitation rate (unused)
         &        msg_level = msg_level   ,& !< in   : message leve
         &        l_cv    = .TRUE.     )     !< in   : if temp. is changed for const. volumes

    !

    ! Final saturation adjustment
    !
    !
       CALL satad_v_3d( maxiter  = 10              ,& !> in
            &           idim     = nproma          ,& !> in
            &           kdim     = jke             ,& !> in
            &           ilo      = jcs             ,& !> in
            &           iup      = jce             ,& !> in
            &           klo      = jks             ,& !> in
            &           kup      = jke             ,& !> in
            &           tol      = 1.e-3_wp        ,& !> in
            &           te       = zta       (:,:) ,& !> inout
            &           qve      = zqv       (:,:) ,& !> inout
            &           qce      = zqc       (:,:) ,& !> inout
            &           rhotot   = rho       (:,:) )
    !
    ! Calculate tendencies and convert temperature tendency, as computed
    ! in satad/graupel for constant volume to constant pressure
    !
    DO jk = jks,jke
       DO jc = jcs,jce
          tend_ta (jc,jk) =     (zta (jc,jk)-ta (jc,jk))*zdtr
          tend_qv (jc,jk) = MAX((zqv (jc,jk)-qv (jc,jk))*zdtr,-qv (jc,jk)*zdtr)
          tend_qc (jc,jk) = MAX((zqc (jc,jk)-qc (jc,jk))*zdtr,-qc (jc,jk)*zdtr)
          tend_qnc(jc,jk) = MAX((zqnc(jc,jk)-qnc(jc,jk))*zdtr,-qnc(jc,jk)*zdtr)
          tend_qi (jc,jk) = MAX((zqi (jc,jk)-qi (jc,jk))*zdtr,-qi (jc,jk)*zdtr)
          tend_qni(jc,jk) = MAX((zqni(jc,jk)-qni(jc,jk))*zdtr,-qni(jc,jk)*zdtr)
          tend_qr (jc,jk) = MAX((zqr (jc,jk)-qr (jc,jk))*zdtr,-qr (jc,jk)*zdtr)
          tend_qnr(jc,jk) = MAX((zqnr(jc,jk)-qnr(jc,jk))*zdtr,-qnr(jc,jk)*zdtr)
          tend_qs (jc,jk) = MAX((zqs (jc,jk)-qs (jc,jk))*zdtr,-qs (jc,jk)*zdtr)
          tend_qns(jc,jk) = MAX((zqns(jc,jk)-qns(jc,jk))*zdtr,-qns(jc,jk)*zdtr)
          tend_qg (jc,jk) = MAX((zqg (jc,jk)-qg (jc,jk))*zdtr,-qg (jc,jk)*zdtr)
          tend_qng(jc,jk) = MAX((zqng(jc,jk)-qng(jc,jk))*zdtr,-qng(jc,jk)*zdtr)
          tend_qh (jc,jk) = MAX((zqh (jc,jk)-qh (jc,jk))*zdtr,-qh (jc,jk)*zdtr)
          tend_qnh(jc,jk) = MAX((zqnh(jc,jk)-qnh(jc,jk))*zdtr,-qnh(jc,jk)*zdtr)
          tend_ninact(jc,jk) = MAX((zninact(jc,jk)-ninact(jc,jk))*zdtr,-ninact(jc,jk)*zdtr)
       END DO
    END DO

  END SUBROUTINE cloud_two

END MODULE mo_cloud_two
