!
! Subroutine interface_aes_rht calls the radiative heating scheme.
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

MODULE mo_interface_aes_rht

  USE mo_kind,                       ONLY: wp
  USE mo_exception,                  ONLY: finish

  USE mo_aes_phy_dims,               ONLY: aes_phy_dims
  USE mo_aes_phy_config,             ONLY: aes_phy_config, aes_phy_tc
  USE mo_aes_phy_memory,             ONLY: t_aes_phy_field, prm_field, &
    &                                      t_aes_phy_tend,  prm_tend

  USE mo_radheating,                 ONLY: radheating
  USE mo_radiation_solar_data,       ONLY: psctm

  USE mo_timer,                      ONLY: ltimer, timer_start, timer_stop, timer_rht

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_rht

CONTAINS

  SUBROUTINE interface_aes_rht(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Pointers
    !
    TYPE(t_aes_phy_field)   ,POINTER    :: field
    TYPE(t_aes_phy_tend)    ,POINTER    :: tend

    ! Shortcuts
    !
    INTEGER  :: fc_rht

    ! Local variables
    !
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: nlevp1, jc, jk
    !
    REAL(wp) :: q_rad(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: q_rlw(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    REAL(wp) :: q_rsw(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    !
    REAL(wp) :: tend_ta_rad(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)

    IF (ltimer) CALL timer_start(timer_rht)

    nlev   = aes_phy_dims(jg)%nlev
    nproma = aes_phy_dims(jg)%nproma

    pdtime               = aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval = aes_phy_tc(jg)%is_in_sd_ed_interval_rad
    is_active            = aes_phy_tc(jg)%is_active_rad

    ! associate pointers
    field  => prm_field(jg)
    tend   => prm_tend (jg)

    fc_rht = aes_phy_config(jg)%fc_rht

    IF ( is_in_sd_ed_interval ) THEN
       !$ACC DATA PRESENT(field, tend) &
       !$ACC   CREATE(q_rad, q_rsw, q_rlw, tend_ta_rad)
       !
       IF (is_active) THEN
          !
          nlevp1 = nlev+1
          !
          CALL radheating (                                   &
               !
               ! input
               ! -----
               !
               & jcs        = jcs                            ,&! loop start index
               & jce        = jce                            ,&! loop end index
               & kbdim      = nproma                         ,&! dimension size
               & klev       = nlev                           ,&! vertical dimension size
               & klevp1     = nlevp1                         ,&! vertical dimension size
               !
               & rsdt0      = psctm(jg)                      ,&! toa incident shortwave radiation for sun in zenith
               & cosmu0     = field%cosmu0    (:,jb)         ,&! solar zenith angle at current time
               & daylght_frc= field%daylght_frc(:,jb)        ,&! daylight fraction
               !
               & emiss      = field%emissivity (:,jb)        ,&! lw sfc emissivity
               & tsr        = field%ts_rad (:,jb)            ,&! radiative surface temperature at current   time [K]
               & tsr_rt     = field%ts_rad_rt(:,jb)          ,&! radiative surface temperature at radiation time [K]
               !
               & rsd_rt     = field%rsd_rt           (:,:,jb),&! all-sky   shortwave downward flux at radiation time [W/m2]
               & rsu_rt     = field%rsu_rt           (:,:,jb),&! all-sky   shortwave upward   flux at radiation time [W/m2]
               !
               & rsdcs_rt   = field%rsdcs_rt         (:,:,jb),&! clear-sky shortwave downward flux at radiation time [W/m2]
               & rsucs_rt   = field%rsucs_rt         (:,:,jb),&! clear-sky shortwave upward   flux at radiation time [W/m2]
               !
               & rld_rt     = field%rld_rt           (:,:,jb),&! all-sky   longwave  downward flux at radiation time [W/m2]
               & rlu_rt     = field%rlu_rt           (:,:,jb),&! all-sky   longwave  upward   flux at radiation time [W/m2]
               !
               & rldcs_rt   = field%rldcs_rt         (:,:,jb),&! clear-sky longwave  downward flux at radiation time [W/m2]
               & rlucs_rt   = field%rlucs_rt         (:,:,jb),&! clear-sky longwave  upward   flux at radiation time [W/m2]
               !
               & rvds_dir_rt= field%rvds_dir_rt        (:,jb),&!< out  all-sky downward direct visible radiation at surface
               & rpds_dir_rt= field%rpds_dir_rt        (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
               & rnds_dir_rt= field%rnds_dir_rt        (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
               & rvds_dif_rt= field%rvds_dif_rt        (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
               & rpds_dif_rt= field%rpds_dif_rt        (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
               & rnds_dif_rt= field%rnds_dif_rt        (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
               & rvus_rt    = field%rvus_rt            (:,jb),&!< out  all-sky upward visible radiation at surface
               & rpus_rt    = field%rpus_rt            (:,jb),&!< out  all-sky upward PAR     radiation at surfac
               & rnus_rt    = field%rnus_rt            (:,jb),&!< out  all-sky upward near-IR radiation at surface
               !
               ! output
               ! ------
               !
               & rsdt       = field%rsdt               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
               & rsut       = field%rsut               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
               & rsds       = field%rsds               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
               & rsus       = field%rsus               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
               !
               & rsutcs     = field%rsutcs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
               & rsdscs     = field%rsdscs             (:,jb),&! clear-sky shortwave downward flux at current   time [W/m2]
               & rsuscs     = field%rsuscs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
               !
               & rvds_dir   = field%rvds_dir           (:,jb),&!< out  all-sky downward direct visible radiation at surface
               & rpds_dir   = field%rpds_dir           (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
               & rnds_dir   = field%rnds_dir           (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
               & rvds_dif   = field%rvds_dif           (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
               & rpds_dif   = field%rpds_dif           (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
               & rnds_dif   = field%rnds_dif           (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
               & rvus       = field%rvus               (:,jb),&!< out  all-sky upward visible radiation at surface
               & rpus       = field%rpus               (:,jb),&!< out  all-sky upward PAR     radiation at surfac
               & rnus       = field%rnus               (:,jb),&!< out  all-sky upward near-IR radiation at surface
               !
               & rlut       = field%rlut               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
               & rlds       = field%rlds               (:,jb),&! all-sky   longwave  downward flux at current   time [W/m2]
               & rlus       = field%rlus               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
               !
               & rlutcs     = field%rlutcs             (:,jb),&! clear-sky longwave  upward   flux at current   time [W/m2]
               & rldscs     = field%rldscs             (:,jb),&! clear-sky longwave  downward flux at current   time [W/m2]
               !
               & q_rsw      = q_rsw                    (:,:) ,&! rad. heating by SW           [W/m2]
               & q_rlw      = q_rlw                    (:,:)  )! rad. heating by LW           [W/m2]
          !
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = jcs, jce
              q_rad(jc,jk) = q_rsw(jc,jk)+q_rlw(jc,jk) ! rad. heating by SW+LW        [W/m2]
            END DO
          END DO
          !
          ! for output: SW+LW heating
          IF (ASSOCIATED(field% q_rad)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                field% q_rad(jc,jk,jb) = q_rad(jc,jk)
              END DO
            END DO
          END IF
          !
          ! for output: SW heating
          IF (ASSOCIATED(field% q_rsw)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                field% q_rsw(jc,jk,jb) = q_rsw(jc,jk)
              END DO
            END DO
          END IF
          !
          ! for output: LW heating
          IF (ASSOCIATED(field% q_rlw)) THEN 
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                field% q_rlw(jc,jk,jb) = q_rlw(jc,jk)
              END DO
            END DO
          END IF
          !$ACC END PARALLEL
          !
          !
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          IF (ASSOCIATED(field% q_rad_vi)) THEN
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_rad_vi(jc, jb) = SUM(q_rad(jc,:))
            END DO
          END IF
          !
          IF (ASSOCIATED(field% q_rsw_vi)) THEN
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_rsw_vi(jc,jb) = SUM(q_rsw(jc,:))
            END DO
          END IF
          !
          IF (ASSOCIATED(field% q_rlw_vi)) THEN
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_rlw_vi(jc,jb) = SUM(q_rlw(jc,:))
            END DO
          END IF
          !
          ! store LW heating in lowermost layer separately,
          ! which is needed for computing q_rlw_impl
          !$ACC LOOP GANG VECTOR
          DO jc = jcs, jce
            field% q_rlw_nlev(jc,jb) = q_rlw(jc,nlev)
          END DO
          !$ACC END PARALLEL
          !
       ELSE
          CALL finish('mo_interface_aes_rht','interface_aes_rht must not be called with is_active=.FALSE.')
       END IF
       !
       ! convert    heating
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
       DO jk = 1, nlev
         DO jc = jcs, jce
           tend_ta_rad(jc,jk) = q_rad(jc,jk) / field% cvair(jc,jk,jb)
         END DO
       END DO
       !
       IF (ASSOCIATED(tend% ta_rad)) THEN
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             tend% ta_rad(jc,jk,jb) = tend_ta_rad(jc,jk)
           END DO
         END DO
       END IF
       IF (ASSOCIATED(tend% ta_rsw)) THEN
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             tend% ta_rsw(jc,jk,jb) = q_rsw(jc,jk) / field% cvair(jc,jk,jb)
           END DO
         END DO
       END IF
       IF (ASSOCIATED(tend% ta_rlw)) THEN
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             tend% ta_rlw(jc,jk,jb) = q_rlw(jc,jk) / field% cvair(jc,jk,jb)
           END DO
         END DO
       END IF
       !
       ! accumulate heating
       IF (ASSOCIATED(field% q_phy   )) THEN
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             field% q_phy(jc,jk,jb) = field% q_phy(jc,jk,jb) + q_rad(jc,jk)
           END DO
         END DO
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_rht)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             tend% ta_phy(jc,jk,jb) = tend% ta_phy(jc,jk,jb) + tend_ta_rad (jc,jk)
           END DO
         END DO
       END SELECT
       !
       ! update physics state for input to the next physics process
       SELECT CASE(fc_rht)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the physics state
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = jcs, jce
              field% ta(jc,jk,jb) = field% ta(jc,jk,jb) + tend_ta_rad(jc,jk)*pdtime
            END DO
          END DO
       END SELECT
       !$ACC END PARALLEL

       IF (ASSOCIATED(field% q_phy_vi)) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = jcs, jce
          field% q_phy_vi(jc,jb) = field% q_phy_vi(jc,jb) + SUM(q_rad(jc,:))
        END DO
        !$ACC END PARALLEL
       END IF
       !
       !$ACC WAIT(1)
       !$ACC END DATA
       !
    ELSE
       !
       !$ACC DATA PRESENT(field, tend)
       !
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       IF (ASSOCIATED(field% q_rad)) THEN
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             field% q_rad(jc,jk,jb) = 0.0_wp
           END DO
         END DO
       END IF
       !
       IF (ASSOCIATED(field% q_rsw)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = jcs, jce
            field% q_rsw(jc,jk,jb) = 0.0_wp
          END DO
        END DO
       END IF
       !
       IF (ASSOCIATED(field% q_rlw)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = jcs, jce
            field% q_rlw(jc,jk,jb) = 0.0_wp
          END DO
        END DO
       END IF
       !
       IF (ASSOCIATED(field% q_rlw)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = jcs, jce
            field% q_rlw(jc,jk,jb) = 0.0_wp
          END DO
        END DO
       END IF
       !
       IF (ASSOCIATED(tend% ta_rad)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = jcs, jce
            tend% ta_rad(jc,jk,jb) = 0.0_wp
          END DO
        END DO
       END IF
       !
       IF (ASSOCIATED(tend% ta_rsw)) THEN
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             tend% ta_rsw(jc,jk,jb) = 0.0_wp
           END DO
         END DO
       END IF
       !
       IF (ASSOCIATED(tend% ta_rlw)) THEN
         !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
         DO jk = 1, nlev
           DO jc = jcs, jce
             tend% ta_rlw(jc,jk,jb) = 0.0_wp
           END DO
         END DO
       END IF
       !$ACC END PARALLEL
       !
       !
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       IF (ASSOCIATED(field% q_rad_vi)) THEN
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_rad_vi(jc, jb) = 0.0_wp
         END DO
       END IF
       !
       IF (ASSOCIATED(field% q_rsw_vi)) THEN
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_rsw_vi(jc,jb) = 0.0_wp
         END DO
       END IF
       !
       IF (ASSOCIATED(field% q_rlw_vi)) THEN
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_rlw_vi(jc,jb) = 0.0_wp
         END DO
       END IF
       !
       !$ACC LOOP GANG VECTOR
       DO jc = jcs, jce
         field% q_rlw_nlev(jc,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
    END IF

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_rht)

  END SUBROUTINE interface_aes_rht

END MODULE mo_interface_aes_rht
