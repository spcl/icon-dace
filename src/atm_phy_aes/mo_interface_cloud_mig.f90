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

! Subroutine interface_cloud_mig calls NWP graupel scheme

MODULE mo_interface_cloud_mig

  USE mo_kind                ,ONLY: wp
  USE mo_copy                ,ONLY: copy

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, t_aes_phy_tend, &
       &                                    prm_field,         prm_tend

  USE mo_cloud_mig_types     ,ONLY: t_cloud_mig_input, t_cloud_mig_output
  USE mo_cloud_mig_memory    ,ONLY:   cloud_mig_input,   cloud_mig_output
  USE mo_cloud_mig           ,ONLY:   cloud_mig

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, &
       &                            timer_mig, timer_cld_mig

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_cloud_mig

CONTAINS

  SUBROUTINE interface_cloud_mig(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)       :: jg, jb, jcs, jce

    ! Pointers
    !
    ! to aes_phy_memory
    TYPE(t_aes_phy_field),    POINTER :: field
    TYPE(t_aes_phy_tend),     POINTER :: tend
    !
    ! to cloud_mig_memory
    TYPE(t_cloud_mig_input ), POINTER :: input
    TYPE(t_cloud_mig_output), POINTER :: output

    ! Local variables
    !
    INTEGER  :: nlev
    !
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: fc_mig
    !
    REAL(wp) :: tend_ta_mig(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qv_mig(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qc_mig(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qi_mig(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qr_mig(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qs_mig(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: tend_qg_mig(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    !
    REAL(wp) :: pr_rain(aes_phy_dims(jg)%nproma)
    REAL(wp) :: pr_ice(aes_phy_dims(jg)%nproma)
    REAL(wp) :: pr_snow(aes_phy_dims(jg)%nproma)
    REAL(wp) :: pr_grpl(aes_phy_dims(jg)%nproma)
    REAL(wp) :: pr_eflx(aes_phy_dims(jg)%nproma)
    !
    INTEGER  :: jk, jks, jke
    INTEGER  :: jc

    IF (ltimer) call timer_start(timer_mig)

    nlev    = aes_phy_dims(jg)%nlev

    pdtime               = aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval = aes_phy_tc(jg)%is_in_sd_ed_interval_mig
    is_active            = aes_phy_tc(jg)%is_active_mig

    ! associate pointers
    !
    ! to aes_phy_memory
    field  => prm_field(jg)
    tend   => prm_tend (jg)
    !
    ! to cloud_mig memory
    input  => cloud_mig_input (jg)
    output => cloud_mig_output(jg)

    !$ACC DATA PRESENT(field, tend, input, output) &
    !$ACC   PRESENT(field%qtrc_phy) & ! ACCWA (nvhpc on levante): to prevent illegal address during kernel execution
    !$ACC   CREATE(tend_ta_mig) &
    !$ACC   CREATE(tend_qv_mig, tend_qc_mig, tend_qi_mig) &
    !$ACC   CREATE(tend_qr_mig, tend_qs_mig, tend_qg_mig) &
    !$ACC   CREATE(pr_rain, pr_ice, pr_snow, pr_grpl, pr_eflx)

    fc_mig = aes_phy_config(jg)% fc_mig

    jks    = aes_phy_config(jg)% jks_cloudy
    jke    = nlev

    ! store input in memory for output
    !
    ! input parameters
    !
    IF (ASSOCIATED(input% jcs       )) CALL copy(jcs,jce, jcs      , input% jcs       (:,jb))
    IF (ASSOCIATED(input% jce       )) CALL copy(jcs,jce, jce      , input% jce       (:,jb))
    IF (ASSOCIATED(input% pdtime    )) CALL copy(jcs,jce, pdtime   , input% pdtime    (:,jb))
    !
    ! input fields
    !
    IF (ASSOCIATED(input% dz    )) CALL copy(jcs,jce, jks,jke, field% dz        (:,:,jb)    , input% dz    (:,:,jb))
    IF (ASSOCIATED(input% rho   )) CALL copy(jcs,jce, jks,jke, field% rho       (:,:,jb)    , input% rho   (:,:,jb))
    IF (ASSOCIATED(input% pf    )) CALL copy(jcs,jce, jks,jke, field% pfull     (:,:,jb)    , input% pf    (:,:,jb))
    IF (ASSOCIATED(input% ta    )) CALL copy(jcs,jce, jks,jke, field% ta        (:,:,jb)    , input% ta    (:,:,jb))
    IF (ASSOCIATED(input% qv    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy  (:,:,jb,iqv), input% qv    (:,:,jb))
    IF (ASSOCIATED(input% qc    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy  (:,:,jb,iqc), input% qc    (:,:,jb))
    IF (ASSOCIATED(input% qi    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy  (:,:,jb,iqi), input% qi    (:,:,jb))
    IF (ASSOCIATED(input% qr    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy  (:,:,jb,iqr), input% qr    (:,:,jb))
    IF (ASSOCIATED(input% qs    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy  (:,:,jb,iqs), input% qs    (:,:,jb))
    IF (ASSOCIATED(input% qg    )) CALL copy(jcs,jce, jks,jke, field% qtrc_phy  (:,:,jb,iqg), input% qg    (:,:,jb))

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN

          IF (ltimer) CALL timer_start(timer_cld_mig)
          !
          CALL cloud_mig( jcs, jce                            ,& !< in : column index range
               &          pdtime                              ,& !< in : timestep
               &          field% dz        (:,jks:jke,jb)     ,& !< in : vertical layer thickness
               &          field% rho       (:,jks:jke,jb)     ,& !< in : density
               &          field% pfull     (:,jks:jke,jb)     ,& !< in : pressure
               &          field% ta        (:,jks:jke,jb)     ,& !< in : temperature
               &          field% qtrc_phy  (:,jks:jke,jb,iqv) ,& !< in : sp humidity
               &          field% qtrc_phy  (:,jks:jke,jb,iqc) ,& !< in : cloud water
               &          field% qtrc_phy  (:,jks:jke,jb,iqi) ,& !< in : ice
               &          field% qtrc_phy  (:,jks:jke,jb,iqr) ,& !< in : rain
               &          field% qtrc_phy  (:,jks:jke,jb,iqs) ,& !< in : snow
               &          field% qtrc_phy  (:,jks:jke,jb,iqg) ,& !< in : graupel
               &          tend_ta_mig      (:,jks:jke)        ,& !< out: tendency of temperature
               &          tend_qv_mig      (:,jks:jke)        ,& !< out: tendency of water vapor
               &          tend_qc_mig      (:,jks:jke)        ,& !< out: tendency of cloud water
               &          tend_qi_mig      (:,jks:jke)        ,& !< out: tendency of cloud ice
               &          tend_qr_mig      (:,jks:jke)        ,& !< out: tendency of rain
               &          tend_qs_mig      (:,jks:jke)        ,& !< out: tendency of snow
               &          tend_qg_mig      (:,jks:jke)        ,& !< out: tendency of graupel 
               &          pr_rain          (:)                ,& !& out: precip rate rain
               &          pr_ice           (:)                ,& !& out: precip rate rain
               &          pr_snow          (:)                ,& !& out: precip rate snow
               &          pr_grpl          (:)                ,& !& out: precip rate graupel
               &          pr_eflx          (:)                )  !& out: precip rate graupel
          !
          IF (ltimer) CALL timer_stop(timer_cld_mig)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO jc = jcs,jce
             field% ufcs(jc,jb) = pr_eflx(jc)          ! = energy flux from precip
             field% rsfl(jc,jb) = pr_rain(jc)          ! = liquid precip rate
             field% ssfl(jc,jb) = pr_ice(jc)  &        ! = frozen precip rate
                  &             + pr_snow(jc) &
                  &             + pr_grpl(jc)
             field% pr  (jc,jb) = field%rsfl(jc,jb) &  ! = total  precip rate
                  &             + field%ssfl(jc,jb)
          END DO
          !$ACC END PARALLEL

          ! store output in memory for output or recycling
          !
          IF (ASSOCIATED(output% tend_ta_mig )) CALL copy(jcs,jce, jks,jke, tend_ta_mig (:,:), output% tend_ta_mig (:,:,jb))
          IF (ASSOCIATED(output% tend_qv_mig )) CALL copy(jcs,jce, jks,jke, tend_qv_mig (:,:), output% tend_qv_mig (:,:,jb))
          IF (ASSOCIATED(output% tend_qc_mig )) CALL copy(jcs,jce, jks,jke, tend_qc_mig (:,:), output% tend_qc_mig (:,:,jb))
          IF (ASSOCIATED(output% tend_qi_mig )) CALL copy(jcs,jce, jks,jke, tend_qi_mig (:,:), output% tend_qi_mig (:,:,jb))
          IF (ASSOCIATED(output% tend_qr_mig )) CALL copy(jcs,jce, jks,jke, tend_qr_mig (:,:), output% tend_qr_mig (:,:,jb))
          IF (ASSOCIATED(output% tend_qs_mig )) CALL copy(jcs,jce, jks,jke, tend_qs_mig (:,:), output% tend_qs_mig (:,:,jb))
          IF (ASSOCIATED(output% tend_qg_mig )) CALL copy(jcs,jce, jks,jke, tend_qg_mig (:,:), output% tend_qg_mig (:,:,jb))
          !
          IF (ASSOCIATED(output% pr_rain     )) CALL copy(jcs,jce, pr_rain (:), output% pr_rain (:,jb))
          IF (ASSOCIATED(output% pr_ice      )) CALL copy(jcs,jce, pr_ice  (:), output% pr_ice  (:,jb))
          IF (ASSOCIATED(output% pr_snow     )) CALL copy(jcs,jce, pr_snow (:), output% pr_snow (:,jb))
          IF (ASSOCIATED(output% pr_grpl     )) CALL copy(jcs,jce, pr_grpl (:), output% pr_grpl (:,jb))
          IF (ASSOCIATED(output% pr_eflx     )) CALL copy(jcs,jce, pr_eflx (:), output% pr_eflx (:,jb))
          !
       ELSE    ! is_active
          !
          ! retrieve output from memory for recycling
          !
          IF (ASSOCIATED(output% tend_ta_mig )) CALL copy(jcs,jce, jks,jke, output% tend_ta_mig (:,:,jb), tend_ta_mig (:,:))
          IF (ASSOCIATED(output% tend_qv_mig )) CALL copy(jcs,jce, jks,jke, output% tend_qv_mig (:,:,jb), tend_qv_mig (:,:))
          IF (ASSOCIATED(output% tend_qc_mig )) CALL copy(jcs,jce, jks,jke, output% tend_qc_mig (:,:,jb), tend_qc_mig (:,:))
          IF (ASSOCIATED(output% tend_qi_mig )) CALL copy(jcs,jce, jks,jke, output% tend_qi_mig (:,:,jb), tend_qi_mig (:,:))
          IF (ASSOCIATED(output% tend_qr_mig )) CALL copy(jcs,jce, jks,jke, output% tend_qr_mig (:,:,jb), tend_qr_mig (:,:))
          IF (ASSOCIATED(output% tend_qs_mig )) CALL copy(jcs,jce, jks,jke, output% tend_qs_mig (:,:,jb), tend_qs_mig (:,:))
          IF (ASSOCIATED(output% tend_qg_mig )) CALL copy(jcs,jce, jks,jke, output% tend_qg_mig (:,:,jb), tend_qg_mig (:,:))
          !
       END IF  ! is_active
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_mig)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = jks,jke
            DO jc = jcs,jce
              ! use tendency to update the model state
              tend%   ta_phy(jc,jk,jb)      = tend%   ta_phy(jc,jk,jb)     + tend_ta_mig(jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqv)  = tend% qtrc_phy(jc,jk,jb,iqv) + tend_qv_mig(jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqc)  = tend% qtrc_phy(jc,jk,jb,iqc) + tend_qc_mig(jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqi)  = tend% qtrc_phy(jc,jk,jb,iqi) + tend_qi_mig(jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqr)  = tend% qtrc_phy(jc,jk,jb,iqr) + tend_qr_mig(jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqs)  = tend% qtrc_phy(jc,jk,jb,iqs) + tend_qs_mig(jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqg)  = tend% qtrc_phy(jc,jk,jb,iqg) + tend_qg_mig(jc,jk)
            END DO
          END DO
          !$ACC END PARALLEL
       END SELECT
       !
       ! update physics state for input to the next physics process
       SELECT CASE(fc_mig)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = jks,jke
            DO jc = jcs,jce
              field%       ta(jc,jk,jb)      = field%       ta(jc,jk,jb)      + tend_ta_mig(jc,jk)*pdtime
              field% qtrc_phy(jc,jk,jb,iqv)  = field% qtrc_phy(jc,jk,jb,iqv)  + tend_qv_mig(jc,jk)*pdtime
              field% qtrc_phy(jc,jk,jb,iqc)  = field% qtrc_phy(jc,jk,jb,iqc)  + tend_qc_mig(jc,jk)*pdtime
              field% qtrc_phy(jc,jk,jb,iqi)  = field% qtrc_phy(jc,jk,jb,iqi)  + tend_qi_mig(jc,jk)*pdtime
              field% qtrc_phy(jc,jk,jb,iqr)  = field% qtrc_phy(jc,jk,jb,iqr)  + tend_qr_mig(jc,jk)*pdtime
              field% qtrc_phy(jc,jk,jb,iqs)  = field% qtrc_phy(jc,jk,jb,iqs)  + tend_qs_mig(jc,jk)*pdtime
              field% qtrc_phy(jc,jk,jb,iqg)  = field% qtrc_phy(jc,jk,jb,iqg)  + tend_qg_mig(jc,jk)*pdtime
            END DO
          END DO
          !$ACC END PARALLEL
       END SELECT
       !
    ELSE       ! is_in_sd_ed_interval
       !
       ! initialize
       !
       IF (ASSOCIATED(output% tend_ta_mig )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_ta_mig (:,:,jb))
       IF (ASSOCIATED(output% tend_qv_mig )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qv_mig (:,:,jb))
       IF (ASSOCIATED(output% tend_qc_mig )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qc_mig (:,:,jb))
       IF (ASSOCIATED(output% tend_qi_mig )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qi_mig (:,:,jb))
       IF (ASSOCIATED(output% tend_qr_mig )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qr_mig (:,:,jb))
       IF (ASSOCIATED(output% tend_qs_mig )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qs_mig (:,:,jb))
       IF (ASSOCIATED(output% tend_qg_mig )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qg_mig (:,:,jb))
       !
       IF (ASSOCIATED(output% pr_rain )) CALL copy(jcs,jce, 0._wp, output% pr_rain (:,jb))
       IF (ASSOCIATED(output% pr_ice  )) CALL copy(jcs,jce, 0._wp, output% pr_ice  (:,jb))
       IF (ASSOCIATED(output% pr_snow )) CALL copy(jcs,jce, 0._wp, output% pr_snow (:,jb))
       IF (ASSOCIATED(output% pr_grpl )) CALL copy(jcs,jce, 0._wp, output% pr_grpl (:,jb))
       IF (ASSOCIATED(output% pr_eflx )) CALL copy(jcs,jce, 0._wp, output% pr_eflx (:,jb))
       !
    END IF     ! is_in_sd_ed_interval

    !$ACC WAIT(1)
    !$ACC END DATA

    ! disassociate pointers
    NULLIFY(field, tend)
    NULLIFY(input, output)

    IF (ltimer) CALL timer_stop(timer_mig)

  END SUBROUTINE interface_cloud_mig

END MODULE mo_interface_cloud_mig
