# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/mo_interface_aes_rad.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/mo_interface_aes_rad.f90"
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

! Subroutine interface_aes_rad calls the radiative transfer scheme.

!----------------------------

# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/include/omp_definitions.inc" 1
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

# 16 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/mo_interface_aes_rad.f90" 2
!----------------------------

MODULE mo_interface_aes_rad

  USE mo_kind,                 ONLY: wp
  USE mtime,                   ONLY: t_datetime => datetime, OPERATOR(<=), OPERATOR(>)

  USE mo_aes_phy_dims,         ONLY: aes_phy_dims
  USE mo_aes_phy_config,       ONLY: aes_phy_tc
  USE mo_aes_rad_config,       ONLY: aes_rad_config
  USE mo_aes_phy_memory,       ONLY: t_aes_phy_field, prm_field
  USE mo_snow_ice_reff,        ONLY: ice_reff_moss, snow_reff_funeedles, snow_x, ice_x
  USE mo_run_config,           ONLY: iqs, iqi
  USE mo_aes_graupel,          ONLY: snow_number, snow_lambda, ice_number





  USE mo_timer,                ONLY: ltimer, timer_start, timer_stop, timer_rad

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_aes_rad

CONTAINS

  SUBROUTINE interface_aes_rad(jg, jb, jcs, jce)

    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Pointers
    !
    TYPE(t_aes_phy_field), POINTER :: field

    ! Local variables
    !
    INTEGER  :: ntracer
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    TYPE(t_datetime), POINTER :: datetime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER :: i1, i2, i3
    LOGICAL :: loland(aes_phy_dims(jg)%nproma)
    LOGICAL :: loglac(aes_phy_dims(jg)%nproma)
    !
    ! temp variable since non-contiguous slicing is not supported in OpenACC
    REAL(wp):: qtrc_phy(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev,aes_phy_dims(jg)%ntracer)
    !
    LOGICAL :: lclrsky_lw, lclrsky_sw

    IF (ltimer) CALL timer_start(timer_rad)

    ntracer = aes_phy_dims(jg)%ntracer
    nlev    = aes_phy_dims(jg)%nlev
    nproma  = aes_phy_dims(jg)%nproma

    datetime             => aes_phy_tc(jg)%datetime
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_rad
    is_active            =  aes_phy_tc(jg)%is_active_rad

    lclrsky_lw = aes_rad_config(jg)%lclrsky_lw
    lclrsky_sw = aes_rad_config(jg)%lclrsky_sw

    ! associate pointers
    field => prm_field(jg)

    IF ( is_in_sd_ed_interval ) THEN
        !
# 205 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/mo_interface_aes_rad.f90"
          !
       ELSE
          !
          ! LW
          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
          field%rld_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          field%rlu_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          !
          ! SW all
          field%rsd_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          field%rsu_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          !
          ! SW vis, par and nir
          field%rvds_dir_rt(:,:) = 0.0_wp !< out  all-sky downward direct visible radiation at surface
          field%rpds_dir_rt(:,:) = 0.0_wp !< all-sky downward direct PAR     radiation at surface
          field%rnds_dir_rt(:,:) = 0.0_wp !< all-sky downward direct near-IR radiation at surface
          field%rvds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse visible radiation at surface
          field%rpds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse PAR     radiation at surface
          field%rnds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse near-IR radiation at surface
          field%rvus_rt    (:,:) = 0.0_wp !< all-sky upward visible radiation at surface
          field%rpus_rt    (:,:) = 0.0_wp !< all-sky upward PAR     radiation at surfac
          field%rnus_rt    (:,:) = 0.0_wp !< all-sky upward near-IR radiation at surface
          !
          ! total cloud cover diagnostics
          field%aclcov(:,:)      = 0.0_wp !< out  total cloud cover
          ! cloud ice and snow optical depth integrated over all bands
          field%tau_ice(:,:,:) = 0.0_wp
          field%tau_snow(:,:,:) = 0.0_wp
          !$ACC END KERNELS
          !
          ! LW clear sky
          IF (lclrsky_lw) THEN
             !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
             field%rldcs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net longwave  at all levels
             field%rlucs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net longwave  at all levels
             !$ACC END KERNELS
          END IF
          !
          ! SW clear sky
          IF (lclrsky_sw) THEN
             !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
             field%rsdcs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net shortwave at all levels
             field%rsucs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net shortwave at all levels
             !$ACC END KERNELS
          END IF
          !
       END IF

     IF (ltimer) CALL timer_stop(timer_rad)

  END SUBROUTINE interface_aes_rad

END MODULE mo_interface_aes_rad
