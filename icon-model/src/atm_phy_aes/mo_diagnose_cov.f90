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

! Module containing subroutine computing the 0/1 cloud cover.

MODULE mo_diagnose_cov

  USE mo_kind                ,ONLY: wp

  USE mo_run_config          ,ONLY: iqc, iqi
  USE mo_aes_phy_config      ,ONLY: aes_phy_config
  USE mo_aes_cov_config      ,ONLY: aes_cov_config

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_memory      ,ONLY: prm_field

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cov

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: diagnose_cov

CONTAINS

  SUBROUTINE diagnose_cov(jg, jb, jcs, jce) 

    ! Arguments
    !
    INTEGER, INTENT(in)  :: jg, jb, jcs, jce

    ! Local variables
    !
    INTEGER  :: jc, jk, jks, jke
    REAL(wp) :: cqx, zqx

    IF (ltimer) CALL timer_start(timer_cov)

    jks = aes_phy_config(jg)%jks_cloudy
    jke = aes_phy_dims(jg)%nlev

    cqx = aes_cov_config(jg)%cqx

    ! Calculate the binary, 0 or 1, cloud cover in a cell based on the  mass fraction of
    ! cloud water + cloud ice in air, in levels jks_cloudy to klev. There cloud cover is 1
    ! if the cloud mass fraction in air exceeds the critical value cqx, and 0 otherwise.
    ! Cloud cover in the levels 1 to jks-1 is 0.
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,jks-1
       DO jc = jcs,jce
          !
          prm_field(jg)%aclc(jc,jk,jb) = 0.0_wp
          !
       END DO  !jc
    END DO   !jk
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = jks,jke
       DO jc = jcs,jce
          !
          zqx =  prm_field(jg)%qtrc_phy(jc,jk,jb,iqc) + prm_field(jg)%qtrc_phy(jc,jk,jb,iqi)
          !
          prm_field(jg)%aclc(jc,jk,jb) = MERGE(1.0_wp, 0.0_wp, zqx >= cqx)
          !
       END DO  !jc
    END DO   !jk
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    ! Calculate the total cloud cover as the maximum cloud cover in the column based on the
    ! binary cloud cover in the cells of the column. Thus total cloud cover is also 0 or 1.
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,jce
       !
       prm_field(jg)%aclcov(jc,jb) = MAXVAL(prm_field(jg)%aclc(jc,jks:jke,jb))
       !
    END DO  !jc
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    IF (ltimer) CALL timer_stop(timer_cov)

  END SUBROUTINE diagnose_cov

END MODULE mo_diagnose_cov
