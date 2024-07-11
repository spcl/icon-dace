!
! Subroutine interface_cloud_two calls NWP two-moment bulk microphysics
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

MODULE mo_interface_cloud_two

  USE mo_kind                ,ONLY: wp
  USE mo_copy                ,ONLY: copy

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqh,  &
       &                            iqnc, iqni, iqnr, iqns, iqng, iqnh, &
       &                            ininact, msg_level

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field, &
       &                            t_aes_phy_tend,  prm_tend

  USE mo_cloud_two_types     ,ONLY: t_cloud_two_input, t_cloud_two_output
  USE mo_cloud_two_memory    ,ONLY:   cloud_two_input,   cloud_two_output
  USE mo_cloud_two           ,ONLY:   cloud_two

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_two

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_cloud_two

CONTAINS

  SUBROUTINE interface_cloud_two(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)       :: jg, jb, jcs, jce

    ! Pointers
    !
    ! to aes_phy_memory
    TYPE(t_aes_phy_field),    POINTER :: field
    TYPE(t_aes_phy_tend),     POINTER :: tend
    !
    ! to cloud_two_memory
    TYPE(t_cloud_two_input ), POINTER :: input
    TYPE(t_cloud_two_output), POINTER :: output

    ! Dummy pointer to tke
    REAL(wp)                ,POINTER    :: ptr_tke(:,:)

    ! Local variables
    !
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: fc_two
    !
    REAL(wp) :: tend_ta_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qv_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qc_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qnc_two   (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qi_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qni_two   (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qr_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qnr_two   (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qs_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qns_two   (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qg_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qng_two   (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qh_two    (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qnh_two   (aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_ninact_two(aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    !
    ! Dummy variable for 2mom scheme. Used only when
    !  assimilation of radar data using latent heat nudging
    !  is switched on (ldass_lhn=true). 
    REAL(wp) :: zqrsflux(aes_phy_dims(jg)%nproma, aes_phy_dims(jg)%nlev)
    !
    INTEGER  :: jk, jks, jke, jl
    INTEGER  :: jc

    nlev    = aes_phy_dims(jg)%nlev
    nproma  = aes_phy_dims(jg)%nproma

    pdtime               =  aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_art
    is_active            =  aes_phy_tc(jg)%is_active_art

    ! associate pointers
    !
    ! to aes_phy_memory
    field  => prm_field(jg)
    tend   => prm_tend (jg)
    !
    ! to cloud_two memory
    input  => cloud_two_input (jg)
    output => cloud_two_output(jg)

    ! associate pointers
    jks       =  aes_phy_config(jg)%jks_cloudy
    fc_two    =  aes_phy_config(jg)%fc_two
    
    ptr_tke => NULL()

    jke       =  nlev

    ! security
    tend_ta_two    (:,:)   = 0.0_wp
    tend_qv_two    (:,:)   = 0.0_wp
    tend_qc_two    (:,:)   = 0.0_wp
    tend_qnc_two   (:,:)   = 0.0_wp
    tend_qi_two    (:,:)   = 0.0_wp
    tend_qni_two   (:,:)   = 0.0_wp
    tend_qr_two    (:,:)   = 0.0_wp
    tend_qnr_two   (:,:)   = 0.0_wp
    tend_qs_two    (:,:)   = 0.0_wp
    tend_qns_two   (:,:)   = 0.0_wp
    tend_qg_two    (:,:)   = 0.0_wp
    tend_qng_two   (:,:)   = 0.0_wp
    tend_qh_two    (:,:)   = 0.0_wp
    tend_qnh_two   (:,:)   = 0.0_wp
    tend_ninact_two(:,:)   = 0.0_wp

    zqrsflux(:,:)  = 0.0_wp

    ! store input in memory for output
    !
    ! input parameters
    !
    IF (ASSOCIATED(input% jcs       )) CALL copy(jcs,jce, jcs      , input% jcs (:,jb))
    IF (ASSOCIATED(input% jce       )) CALL copy(jcs,jce, jce      , input% jce (:,jb))
    IF (ASSOCIATED(input% msg_level )) CALL copy(jcs,jce, msg_level, input% msg_level (:,jb))
    IF (ASSOCIATED(input% pdtime    )) CALL copy(jcs,jce, pdtime   , input% pdtime    (:,jb))
    !
    ! input fields, in
    !
    IF (ASSOCIATED(input% dz    )) CALL copy(jcs,jce, jks,jke, field% dz        (:,:,jb)    , input% dz    (:,:,jb))
    IF (ASSOCIATED(input% zh    )) CALL copy(jcs,jce, jks,jke, field% zh        (:,:,jb)    , input% zh    (:,:,jb))
    IF (ASSOCIATED(input% rho   )) CALL copy(jcs,jce, jks,jke, field% rho       (:,:,jb)    , input% rho   (:,:,jb))
    IF (ASSOCIATED(input% pf    )) CALL copy(jcs,jce, jks,jke, field% pfull     (:,:,jb)    , input% pf    (:,:,jb))
    !
    ! input fields, inout
    !
    IF (ASSOCIATED(input% ta    )) CALL copy(jcs,jce, jks,jke, field% ta   (:,:,jb)        , input% ta    (:,:,jb))
    IF (ASSOCIATED(input% qv    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqv)    , input% qv    (:,:,jb))
    IF (ASSOCIATED(input% qc    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqc)    , input% qc    (:,:,jb))
    IF (ASSOCIATED(input% qnc   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqnc)   , input% qnc   (:,:,jb))
    IF (ASSOCIATED(input% qi    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqi)    , input% qi    (:,:,jb))
    IF (ASSOCIATED(input% qni   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqni)   , input% qni   (:,:,jb))
    IF (ASSOCIATED(input% qr    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqr)    , input% qr    (:,:,jb))
    IF (ASSOCIATED(input% qnr   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqnr)   , input% qnr   (:,:,jb))
    IF (ASSOCIATED(input% qs    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqs)    , input% qs    (:,:,jb))
    IF (ASSOCIATED(input% qns   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqns)   , input% qns   (:,:,jb))
    IF (ASSOCIATED(input% qg    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqg)    , input% qg    (:,:,jb))
    IF (ASSOCIATED(input% qng   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqng)   , input% qng   (:,:,jb))
    IF (ASSOCIATED(input% qh    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqh)    , input% qh    (:,:,jb))
    IF (ASSOCIATED(input% qnh   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqnh)   , input% qnh   (:,:,jb))
    IF (ASSOCIATED(input% ninact)) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,ininact), input% ninact(:,:,jb))
    IF (ASSOCIATED(input% w     )) CALL copy(jcs,jce, jks,jke+1, field% wa(:,:,jb)           , input% w     (:,:,jb))
    !
    IF (ASSOCIATED(input% pr_rain)) CALL copy(jcs,jce, field% rain_gsp_rate (:,jb), input% pr_rain   (:,jb))
    IF (ASSOCIATED(input% pr_ice )) CALL copy(jcs,jce, field%  ice_gsp_rate (:,jb), input% pr_ice    (:,jb))
    IF (ASSOCIATED(input% pr_snow)) CALL copy(jcs,jce, field% snow_gsp_rate (:,jb), input% pr_snow   (:,jb))
    IF (ASSOCIATED(input% pr_grpl)) CALL copy(jcs,jce, field% graupel_gsp_rate (:,jb), input% pr_grpl   (:,jb))
    IF (ASSOCIATED(input% pr_hail)) CALL copy(jcs,jce, field% hail_gsp_rate (:,jb), input% pr_hail   (:,jb))

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN

          IF (ltimer) CALL timer_start(timer_two)
          !

          CALL cloud_two( nproma, nlev                          ,& !< in : grid index
               &          jcs, jce                              ,& !< in : column index range
               &          jks, jke                              ,& !< in : column index range
               &          msg_level                             ,& !< in : message level 
               &          pdtime                                ,& !< in : timestep
               &          field% dz        (:,:,jb)       ,& !< in : vertical layer thickness
               &          field% zh        (:,:,jb)       ,& !< in : height of half levels
               &          field% rho       (:,:,jb)       ,& !< in : density
               &          field% pfull     (:,:,jb)       ,& !< in : pressure
               &          ptr_tke                         ,& !< in : TKE  (dummy because of no TKE in AES)
               &          field% ta        (:,:,jb)       ,& !< inout : temperature
               &          field% qtrc_phy  (:,:,jb,iqv)   ,& !< inout : sp humidity
               &          field% qtrc_phy  (:,:,jb,iqc)   ,& !< inout : cloud water
               &          field% qtrc_phy  (:,:,jb,iqnc)  ,& !< inout : cloud water number
               &          field% qtrc_phy  (:,:,jb,iqi)   ,& !< inout : ice
               &          field% qtrc_phy  (:,:,jb,iqni)  ,& !< inout : ice number
               &          field% qtrc_phy  (:,:,jb,iqr)   ,& !< inout : rain
               &          field% qtrc_phy  (:,:,jb,iqnr)  ,& !< inout : rain number
               &          field% qtrc_phy  (:,:,jb,iqs)   ,& !< inout : snow
               &          field% qtrc_phy  (:,:,jb,iqns)  ,& !< inout : snow number
               &          field% qtrc_phy  (:,:,jb,iqg)   ,& !< inout : graupel
               &          field% qtrc_phy  (:,:,jb,iqng)  ,& !< inout : graupel number
               &          field% qtrc_phy  (:,:,jb,iqh)   ,& !< inout : hail
               &          field% qtrc_phy  (:,:,jb,iqnh)  ,& !< inout : hail number
               &          field% qtrc_phy  (:,:,jb,ininact) ,& !< inout : activated ice nuclei
               &          field% wa        (:,:,jb)     ,& !< in : vertical velocity
               &          tend_ta_two      (:,:)          ,& !< out: tendency of temperature
               &          tend_qv_two      (:,:)          ,& !< out: tendency of water vapor
               &          tend_qc_two      (:,:)          ,& !< out: tendency of cloud water
               &          tend_qnc_two     (:,:)          ,& !< out: tendency of cloud water
               &          tend_qi_two      (:,:)          ,& !< out: tendency of cloud ice
               &          tend_qni_two     (:,:)          ,& !< out: tendency of cloud ice
               &          tend_qr_two      (:,:)          ,& !< out: tendency of rain
               &          tend_qnr_two     (:,:)          ,& !< out: tendency of rain droplets
               &          tend_qs_two      (:,:)          ,& !< out: tendency of snow
               &          tend_qns_two     (:,:)          ,& !< out: tendency of snow droplets
               &          tend_qg_two      (:,:)          ,& !< out: tendency of graupel 
               &          tend_qng_two     (:,:)          ,& !< out: tendency of graupel droplets
               &          tend_qh_two      (:,:)          ,& !< out: tendency of hail
               &          tend_qnh_two     (:,:)          ,& !< out: tendency of hail droplets
               &          tend_ninact_two  (:,:)          ,& !< out: tendency of activated ice
               &          field% rain_gsp_rate (:,jb)           ,& !& inout: precip rate rain
               &          field%  ice_gsp_rate (:,jb)           ,& !& inout: precip rate ice
               &          field% snow_gsp_rate (:,jb)           ,& !& inout: precip rate snow
               &          field% graupel_gsp_rate (:,jb)        ,& !& inout: precip rate graupel
               &          field% hail_gsp_rate (:,jb)           ,& !& inout: precip rate hail
               &          zqrsflux         (:,:)          )  !< inout: unused 3d total prec. rate
          !
          IF (ltimer) CALL timer_stop(timer_two)

          !
          ! Calculate rain and snow
          !
          DO jc = jcs, jce
            field% rsfl(jc,jb) = field%    rain_gsp_rate (jc,jb)    ! = liquid precip rate
            field% ssfl(jc,jb) = field%    snow_gsp_rate (jc,jb) &  ! = frozen precip rate
                 &             + field%     ice_gsp_rate (jc,jb) &
                 &             + field% graupel_gsp_rate (jc,jb) &
                 &             + field%    hail_gsp_rate (jc,jb)
            field% pr  (jc,jb) = field% rsfl(jc,jb) &               ! = total  precip rate
                 &             + field% ssfl(jc,jb)   
          END DO
          !
          !
          ! store output in memory for output or recycling
          !
          IF (ASSOCIATED(output% tend_ta_two   )) CALL copy(jcs,jce, jks,jke, tend_ta_two   (:,:), output% tend_ta_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qv_two   )) CALL copy(jcs,jce, jks,jke, tend_qv_two   (:,:), output% tend_qv_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qc_two   )) CALL copy(jcs,jce, jks,jke, tend_qc_two   (:,:), output% tend_qc_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qnc_two  )) CALL copy(jcs,jce, jks,jke, tend_qnc_two  (:,:), output% tend_qnc_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qi_two   )) CALL copy(jcs,jce, jks,jke, tend_qi_two   (:,:), output% tend_qi_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qni_two  )) CALL copy(jcs,jce, jks,jke, tend_qni_two  (:,:), output% tend_qni_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qr_two   )) CALL copy(jcs,jce, jks,jke, tend_qr_two   (:,:), output% tend_qr_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qnr_two  )) CALL copy(jcs,jce, jks,jke, tend_qnr_two  (:,:), output% tend_qnr_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qs_two   )) CALL copy(jcs,jce, jks,jke, tend_qs_two   (:,:), output% tend_qs_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qns_two  )) CALL copy(jcs,jce, jks,jke, tend_qns_two  (:,:), output% tend_qns_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qg_two   )) CALL copy(jcs,jce, jks,jke, tend_qg_two   (:,:), output% tend_qg_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qng_two  )) CALL copy(jcs,jce, jks,jke, tend_qng_two  (:,:), output% tend_qng_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qh_two   )) CALL copy(jcs,jce, jks,jke, tend_qh_two   (:,:), output% tend_qh_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qnh_two  )) CALL copy(jcs,jce, jks,jke, tend_qnh_two  (:,:), output% tend_qnh_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_ninact_two)) CALL copy(jcs,jce, jks,jke, tend_ninact_two(:,:), output% tend_ninact_two(:,:,jb))
          !
          IF (ASSOCIATED(output% pr_rain     )) CALL copy(jcs,jce, field% rain_gsp_rate (:,jb), output% pr_rain (:,jb))
          IF (ASSOCIATED(output% pr_ice      )) CALL copy(jcs,jce, field%  ice_gsp_rate (:,jb), output% pr_ice  (:,jb))
          IF (ASSOCIATED(output% pr_snow     )) CALL copy(jcs,jce, field% snow_gsp_rate (:,jb), output% pr_snow (:,jb))
          IF (ASSOCIATED(output% pr_grpl     )) CALL copy(jcs,jce, field% graupel_gsp_rate (:,jb), output% pr_grpl (:,jb))
          IF (ASSOCIATED(output% pr_hail     )) CALL copy(jcs,jce, field% hail_gsp_rate (:,jb), output% pr_hail (:,jb))
          !
       ELSE    ! is_active
          !
          ! retrieve output from memory for recycling
          !
          IF (ASSOCIATED(output% tend_ta_two   )) CALL copy(jcs,jce, jks,jke, output% tend_ta_two   (:,:,jb), tend_ta_two   (:,:))
          IF (ASSOCIATED(output% tend_qv_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qv_two   (:,:,jb), tend_qv_two   (:,:))
          IF (ASSOCIATED(output% tend_qc_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qc_two   (:,:,jb), tend_qc_two   (:,:))
          IF (ASSOCIATED(output% tend_qnc_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qnc_two  (:,:,jb), tend_qnc_two  (:,:))
          IF (ASSOCIATED(output% tend_qi_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qi_two   (:,:,jb), tend_qi_two   (:,:))
          IF (ASSOCIATED(output% tend_qni_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qni_two  (:,:,jb), tend_qni_two  (:,:))
          IF (ASSOCIATED(output% tend_qr_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qr_two   (:,:,jb), tend_qr_two   (:,:))
          IF (ASSOCIATED(output% tend_qnr_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qnr_two  (:,:,jb), tend_qnr_two  (:,:))
          IF (ASSOCIATED(output% tend_qs_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qs_two   (:,:,jb), tend_qs_two   (:,:))
          IF (ASSOCIATED(output% tend_qns_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qns_two  (:,:,jb), tend_qns_two  (:,:))
          IF (ASSOCIATED(output% tend_qg_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qg_two   (:,:,jb), tend_qg_two   (:,:))
          IF (ASSOCIATED(output% tend_qng_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qng_two  (:,:,jb), tend_qng_two  (:,:))
          IF (ASSOCIATED(output% tend_qh_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qh_two   (:,:,jb), tend_qh_two   (:,:))
          IF (ASSOCIATED(output% tend_qnh_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qnh_two  (:,:,jb), tend_qnh_two  (:,:))
          IF (ASSOCIATED(output% tend_ninact_two)) CALL copy(jcs,jce, jks,jke, output% tend_ninact_two(:,:,jb), tend_ninact_two(:,:))
          !
       END IF  ! is_active

       !

       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_two)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          DO jk = jks, jke
            DO jl = jcs, jce
              ! use tendency to update the model state
              tend%   ta_phy(jl,jk,jb)        = tend%   ta_phy(jl,jk,jb)         + tend_ta_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqv)    = tend% qtrc_phy(jl,jk,jb,iqv)     + tend_qv_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqc)    = tend% qtrc_phy(jl,jk,jb,iqc)     + tend_qc_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqnc)   = tend% qtrc_phy(jl,jk,jb,iqnc)    + tend_qnc_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqi)    = tend% qtrc_phy(jl,jk,jb,iqi)     + tend_qi_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqni)   = tend% qtrc_phy(jl,jk,jb,iqni)    + tend_qni_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqr)    = tend% qtrc_phy(jl,jk,jb,iqr)     + tend_qr_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqnr)   = tend% qtrc_phy(jl,jk,jb,iqnr)    + tend_qnr_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqs)    = tend% qtrc_phy(jl,jk,jb,iqs)     + tend_qs_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqns)   = tend% qtrc_phy(jl,jk,jb,iqns)    + tend_qns_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqg)    = tend% qtrc_phy(jl,jk,jb,iqg)     + tend_qg_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqng)   = tend% qtrc_phy(jl,jk,jb,iqng)    + tend_qng_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqh)    = tend% qtrc_phy(jl,jk,jb,iqh)     + tend_qh_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqnh)   = tend% qtrc_phy(jl,jk,jb,iqnh)    + tend_qnh_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,ininact)= tend% qtrc_phy(jl,jk,jb,ininact) + tend_ninact_two(jl,jk)
            END DO
          END DO
       END SELECT
       !
       ! update physics state for input to the next physics process
       SELECT CASE(fc_two)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          DO jk = jks, jke
            DO jl = jcs, jce
              field%       ta(jl,jk,jb)        = field%       ta(jl,jk,jb)         + tend_ta_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqv)    = field% qtrc_phy(jl,jk,jb,iqv)     + tend_qv_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqc)    = field% qtrc_phy(jl,jk,jb,iqc)     + tend_qc_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqnc)   = field% qtrc_phy(jl,jk,jb,iqnc)    + tend_qnc_two   (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqi)    = field% qtrc_phy(jl,jk,jb,iqi)     + tend_qi_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqni)   = field% qtrc_phy(jl,jk,jb,iqni)    + tend_qni_two   (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqr)    = field% qtrc_phy(jl,jk,jb,iqr)     + tend_qr_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqnr)   = field% qtrc_phy(jl,jk,jb,iqnr)    + tend_qnr_two   (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqs)    = field% qtrc_phy(jl,jk,jb,iqs)     + tend_qs_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqns)   = field% qtrc_phy(jl,jk,jb,iqns)    + tend_qns_two   (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqg)    = field% qtrc_phy(jl,jk,jb,iqg)     + tend_qg_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqng)   = field% qtrc_phy(jl,jk,jb,iqng)    + tend_qng_two   (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqh)    = field% qtrc_phy(jl,jk,jb,iqh)     + tend_qh_two    (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,iqnh)   = field% qtrc_phy(jl,jk,jb,iqnh)    + tend_qnh_two   (jl,jk)*pdtime
              field% qtrc_phy(jl,jk,jb,ininact)= field% qtrc_phy(jl,jk,jb,ininact) + tend_ninact_two(jl,jk)*pdtime
            END DO
          END DO
          !
          ! output fields, inout
          !
          IF (ASSOCIATED(output% ta    )) CALL copy(jcs,jce, jks,jke, field% ta       (:,:,jb)        , output% ta    (:,:,jb))
          IF (ASSOCIATED(output% qv    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqv)    , output% qv    (:,:,jb))
          IF (ASSOCIATED(output% qc    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqc)    , output% qc    (:,:,jb))
          IF (ASSOCIATED(output% qnc   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqnc)   , output% qnc   (:,:,jb))
          IF (ASSOCIATED(output% qi    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqi)    , output% qi    (:,:,jb))
          IF (ASSOCIATED(output% qni   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqni)   , output% qni   (:,:,jb))
          IF (ASSOCIATED(output% qr    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqr)    , output% qr    (:,:,jb))
          IF (ASSOCIATED(output% qnr   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqnr)   , output% qnr   (:,:,jb))
          IF (ASSOCIATED(output% qs    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqs)    , output% qs    (:,:,jb))
          IF (ASSOCIATED(output% qns   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqns)   , output% qns   (:,:,jb))
          IF (ASSOCIATED(output% qg    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqg)    , output% qg    (:,:,jb))
          IF (ASSOCIATED(output% qng   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqng)   , output% qng   (:,:,jb))
          IF (ASSOCIATED(output% qh    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqh)    , output% qh    (:,:,jb))
          IF (ASSOCIATED(output% qnh   )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,iqnh)   , output% qnh   (:,:,jb))
          IF (ASSOCIATED(output% ninact)) CALL copy(jcs,jce, jks,jke, field% qtrc_phy (:,:,jb,ininact), output% ninact(:,:,jb))
          IF (ASSOCIATED(output% w     )) CALL copy(jcs,jce, jks,jke+1, field% wa       (:,:,jb)        , output% w     (:,:,jb))
          !
       END SELECT
       !
    ELSE       ! is_in_sd_ed_interval : bypass area
       !
       IF (ASSOCIATED(output% tend_ta_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_ta_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qv_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qv_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qc_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qc_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qnc_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qnc_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qi_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qi_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qni_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qni_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qr_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qr_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qnr_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qnr_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qs_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qs_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qns_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qns_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qg_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qg_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qng_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qng_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qh_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qh_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qnh_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qnh_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_ninact_two )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_ninact_two (:,:,jb))
       !
       IF (ASSOCIATED(output% pr_rain )) CALL copy(jcs,jce, 0._wp, output% pr_rain (:,jb))
       IF (ASSOCIATED(output% pr_ice  )) CALL copy(jcs,jce, 0._wp, output% pr_ice  (:,jb))
       IF (ASSOCIATED(output% pr_snow )) CALL copy(jcs,jce, 0._wp, output% pr_snow (:,jb))
       IF (ASSOCIATED(output% pr_grpl )) CALL copy(jcs,jce, 0._wp, output% pr_grpl (:,jb))
       IF (ASSOCIATED(output% pr_hail )) CALL copy(jcs,jce, 0._wp, output% pr_hail (:,jb))
       !
    END IF     ! is_in_sd_ed_interval

    ! disassociate pointers

    NULLIFY(field)
    NULLIFY(tend)

  END SUBROUTINE interface_cloud_two

END MODULE mo_interface_cloud_two
