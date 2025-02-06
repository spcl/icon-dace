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

!> @brief
!! Module containing diagnostics for atmospheric energy
!!
!! Author: Marco Giorgetta, MPI-M, 2024

MODULE mo_atm_energy_diag

  USE mo_kind               ,ONLY: wp                                   !< working precision

  USE mo_parallel_config    ,ONLY: nproma                               !< size of cell dimension
  USE mo_run_config         ,ONLY: num_lev,                           & !< size of vertical dimension
       &                           iqv, iqc, iqi, iqr, iqs, iqg,      & !< tracer indices
       &                           dtime                                !< time step
  USE mo_dynamics_config    ,ONLY: nnow, nnow_rcf, nnew, nnew_rcf       !< time level indices

  USE mo_model_domain       ,ONLY: p => p_patch                         !< memory for patch  variables
  USE mo_nonhydro_state     ,ONLY: s => p_nh_state                      !< memory for state  variables in dynamics
  USE mo_aes_phy_memory     ,ONLY: f => prm_field                       !< memory for state  variables in physics
  USE mo_atm_energy_memory  ,ONLY: e => atm_energy                      !< memory for energy variables

  USE mo_aes_thermo         ,ONLY: internal_energy                      !< function for internal energy
  USE mo_statistics         ,ONLY: horizontal_sum                       !< subroutine for global integral

  USE mo_copy               ,ONLY: copy                                 !< copy variables
  USE mo_fortran_tools      ,ONLY: init                                 !< initialize variables

  USE mo_timer              ,ONLY: ltimer, timer_start, timer_stop,   & !< timer mechanism and
       &                           timer_atm_energy_diag,             & !< timer for energy dignostics
       &                           timer_atm_energy_hint,             & !< timer for global integral of energy
       &                           timer_atm_energy_vint                !< timer for vertical integral of energy

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: atm_energy_diag_d1, atm_energy_diag_p1, atm_energy_hint_1, &
       &    atm_energy_diag_d2, atm_energy_diag_p2, atm_energy_hint_2, &
       &    atm_energy_copy_2_1_3d_vi, atm_energy_copy_2_1_hi_ti, &
       &    atm_energy_copy_2_3_3d_vi, atm_energy_copy_2_3_hi_ti, &
       &    atm_energy_tend_dyn_3d_vi, atm_energy_tend_dyn_hi_ti, &
       &    atm_energy_tend_phy_3d_vi, atm_energy_tend_phy_hi_ti, &
       &    atm_energy_tend_cld_3d_vi, atm_energy_tend_cld_hi_ti, &
       &    atm_energy_tend_rad_3d_vi, atm_energy_tend_rad_hi_ti, &
       &    atm_energy_tend_tmx_3d_vi, atm_energy_tend_tmx_hi_ti

CONTAINS

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_diag_d1(jg, jb, jcs, jce)

    INTEGER, INTENT(in)      :: jg, jb, jcs, jce                        !< horizontal domain indices

    INTEGER                  :: jtl_dyn, jtl_trc                        !< time level indices
    INTEGER                  :: jks, jke                                !< vertical domain indices

    LOGICAL                  :: lein,   lekh,   lekv,   legp,   leto    !< logicals for computation of energy densities
    LOGICAL                  :: leinvi, lekhvi, lekvvi, legpvi, letovi  !< logicals for computation of energy contents

    REAL(wp), DIMENSION(nproma, num_lev(jg)) :: ein, ekh, ekv, egp, eto !< local fields for energy densities
    REAL(wp), DIMENSION(nproma) :: einvi, ekhvi, ekvvi, egpvi, etovi    !< local fields for energy contents

    !$ACC DATA CREATE(ein, ekh, ekv, egp, eto, einvi, ekhvi, ekvvi, egpvi, etovi)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)

    !> Compute energy density (J/m3) and energy content (J/m2)
    !> - Input : from the "now" state of the dynamics memory
    !> - Output: to the "1" state of energy memory
    !
    jtl_dyn = nnow(jg)
    jtl_trc = nnow_rcf(jg)
    !
    jks = 1
    jke = num_lev(jg)
    !
    lein   = ASSOCIATED(e(jg)%ein1)
    lekh   = ASSOCIATED(e(jg)%ekh1)
    lekv   = ASSOCIATED(e(jg)%ekv1)
    legp   = ASSOCIATED(e(jg)%egp1)
    leto   = ASSOCIATED(e(jg)%eto1)
    !
    leinvi = ASSOCIATED(e(jg)%ein1vi)
    lekhvi = ASSOCIATED(e(jg)%ekh1vi)
    lekvvi = ASSOCIATED(e(jg)%ekv1vi)
    legpvi = ASSOCIATED(e(jg)%egp1vi)
    letovi = ASSOCIATED(e(jg)%eto1vi)
    !
    CALL atm_energy_diag(jcs, jce,                                    & !< domain indices
         !
         &               s(jg)%metrics%ddqz_z_full(:,:,jb),           & !< layer thickness
         &               s(jg)%prog(jtl_dyn)%rho(:,:,jb),             & !< mass density
         &               s(jg)%diag%temp(:,:,jb),                     & !< air temperature
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqv),      & !< mass fraction in air of water vapor
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqc),      & !< mass fraction in air of cloud liquid
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqi),      & !< mass fraction in air of cloud ice
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqr),      & !< mass fraction in air of rain
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqs),      & !< mass fraction in air of snow
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqg),      & !< mass fraction in air of graupel
         &               s(jg)%diag%u(:,:,jb),                        & !< zonal velocity
         &               s(jg)%diag%v(:,:,jb),                        & !< meridional velocity
         &               s(jg)%prog(jtl_dyn)%w(:,:,jb),               & !< vertical velocity
         &               s(jg)%metrics%geopot(:,:,jb),                & !< geopotential
         !
         &               lein, ein(:,:),                              & !< internal energy density
         &               lekh, ekh(:,:),                              & !< horizontal kinetic energy density
         &               lekv, ekv(:,:),                              & !< vertical kinetic energy density
         &               legp, egp(:,:),                              & !< geopotential energy density
         &               leto, eto(:,:),                              & !< total energy density
         !
         &               leinvi, einvi(:),                            & !< internal energy content
         &               lekhvi, ekhvi(:),                            & !< horizontal kinetic energy content
         &               lekvvi, ekvvi(:),                            & !< vertical kinetic energy content
         &               legpvi, egpvi(:),                            & !< geopotential energy content
         &               letovi, etovi(:))                              !< total energy content
    !
    !> Store local fields in global fields, where the latter exist
    !
    IF (lein)   CALL copy(jcs, jce, jks, jke, ein(:,:), e(jg)%ein1(:,:,jb))
    IF (lekh)   CALL copy(jcs, jce, jks, jke, ekh(:,:), e(jg)%ekh1(:,:,jb))
    IF (lekv)   CALL copy(jcs, jce, jks, jke, ekv(:,:), e(jg)%ekv1(:,:,jb))
    IF (legp)   CALL copy(jcs, jce, jks, jke, egp(:,:), e(jg)%egp1(:,:,jb))
    IF (leto)   CALL copy(jcs, jce, jks, jke, eto(:,:), e(jg)%eto1(:,:,jb))
    !
    IF (leinvi) CALL copy(jcs, jce, einvi(:), e(jg)%ein1vi(:,jb))
    IF (lekhvi) CALL copy(jcs, jce, ekhvi(:), e(jg)%ekh1vi(:,jb))
    IF (lekvvi) CALL copy(jcs, jce, ekvvi(:), e(jg)%ekv1vi(:,jb))
    IF (legpvi) CALL copy(jcs, jce, egpvi(:), e(jg)%egp1vi(:,jb))
    IF (letovi) CALL copy(jcs, jce, etovi(:), e(jg)%eto1vi(:,jb))

    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

    !$ACC END DATA

  END SUBROUTINE atm_energy_diag_d1

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_diag_d2(jg, jb, jcs, jce)

    INTEGER, INTENT(in)      :: jg, jb, jcs, jce                        !< horizontal domain indices

    INTEGER                  :: jtl_dyn, jtl_trc                        !< time level indices
    INTEGER                  :: jks, jke                                !< vertical domain indices

    LOGICAL                  :: lein,   lekh,   lekv,   legp,   leto    !< logicals for computation of energy densities
    LOGICAL                  :: leinvi, lekhvi, lekvvi, legpvi, letovi  !< logicals for computation of energy contents

    REAL(wp), DIMENSION(nproma, num_lev(jg)) :: ein, ekh, ekv, egp, eto !< local fields for energy densities
    REAL(wp), DIMENSION(nproma) :: einvi, ekhvi, ekvvi, egpvi, etovi    !< local fields for energy contents

    !$ACC DATA CREATE(ein, ekh, ekv, egp, eto, einvi, ekhvi, ekvvi, egpvi, etovi)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)

    !> Compute energy density (J/m3) and energy content (J/m2)
    !> - Input : from the "new" state of the dynamics memory
    !> - Output: to the "2" state of energy memory
    !
    jtl_dyn = nnew(jg)
    jtl_trc = nnew_rcf(jg)
    !
    jks = 1
    jke = num_lev(jg)
    !
    lein   = ASSOCIATED(e(jg)%ein2)
    lekh   = ASSOCIATED(e(jg)%ekh2)
    lekv   = ASSOCIATED(e(jg)%ekv2)
    legp   = ASSOCIATED(e(jg)%egp2)
    leto   = ASSOCIATED(e(jg)%eto2)
    !
    leinvi = ASSOCIATED(e(jg)%ein2vi)
    lekhvi = ASSOCIATED(e(jg)%ekh2vi)
    lekvvi = ASSOCIATED(e(jg)%ekv2vi)
    legpvi = ASSOCIATED(e(jg)%egp2vi)
    letovi = ASSOCIATED(e(jg)%eto2vi)
    !
    CALL atm_energy_diag(jcs, jce,                                    & !< domain indices
         !
         &               s(jg)%metrics%ddqz_z_full(:,:,jb),           & !< layer thickness
         &               s(jg)%prog(jtl_dyn)%rho(:,:,jb),             & !< mass density
         &               s(jg)%diag%temp(:,:,jb),                     & !< air temperature
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqv),      & !< mass fraction in air of water vapor
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqc),      & !< mass fraction in air of cloud liquid
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqi),      & !< mass fraction in air of cloud ice
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqr),      & !< mass fraction in air of rain
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqs),      & !< mass fraction in air of snow
         &               s(jg)%prog(jtl_trc)%tracer(:,:,jb,iqg),      & !< mass fraction in air of graupel
         &               s(jg)%diag%u(:,:,jb),                        & !< zonal velocity
         &               s(jg)%diag%v(:,:,jb),                        & !< meridional velocity
         &               s(jg)%prog(jtl_dyn)%w(:,:,jb),               & !< vertical velocity
         &               s(jg)%metrics%geopot(:,:,jb),                & !< geopotential
         !
         &               lein, ein(:,:),                              & !< internal energy density
         &               lekh, ekh(:,:),                              & !< horizontal kinetic energy density
         &               lekv, ekv(:,:),                              & !< vertical kinetic energy density
         &               legp, egp(:,:),                              & !< geopotential energy density
         &               leto, eto(:,:),                              & !< total energy density
         !
         &               leinvi, einvi(:),                            & !< internal energy content
         &               lekhvi, ekhvi(:),                            & !< horizontal kinetic energy content
         &               lekvvi, ekvvi(:),                            & !< vertical kinetic energy content
         &               legpvi, egpvi(:),                            & !< geopotential energy content
         &               letovi, etovi(:))                              !< total energy content

    !> Store local fields in global fields, where the latter exist
    !
    IF (lein)   CALL copy(jcs, jce, jks, jke, ein(:,:), e(jg)%ein2(:,:,jb))
    IF (lekh)   CALL copy(jcs, jce, jks, jke, ekh(:,:), e(jg)%ekh2(:,:,jb))
    IF (lekv)   CALL copy(jcs, jce, jks, jke, ekv(:,:), e(jg)%ekv2(:,:,jb))
    IF (legp)   CALL copy(jcs, jce, jks, jke, egp(:,:), e(jg)%egp2(:,:,jb))
    IF (leto)   CALL copy(jcs, jce, jks, jke, eto(:,:), e(jg)%eto2(:,:,jb))
    !
    IF (leinvi) CALL copy(jcs, jce, einvi(:), e(jg)%ein2vi(:,jb))
    IF (lekhvi) CALL copy(jcs, jce, ekhvi(:), e(jg)%ekh2vi(:,jb))
    IF (lekvvi) CALL copy(jcs, jce, ekvvi(:), e(jg)%ekv2vi(:,jb))
    IF (legpvi) CALL copy(jcs, jce, egpvi(:), e(jg)%egp2vi(:,jb))
    IF (letovi) CALL copy(jcs, jce, etovi(:), e(jg)%eto2vi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

    !$ACC END DATA

  END SUBROUTINE atm_energy_diag_d2

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_diag_p1(jg, jb, jcs, jce)

    INTEGER, INTENT(in)      :: jg, jb, jcs, jce                        !< horizontal domain indices

    INTEGER                  :: jks, jke                                !< vertical domain indices

    LOGICAL                  :: lein,   lekh,   lekv,   legp,   leto    !< logicals for computation of energy densities
    LOGICAL                  :: leinvi, lekhvi, lekvvi, legpvi, letovi  !< logicals for computation of energy contents

    REAL(wp), DIMENSION(nproma, num_lev(jg)) :: ein, ekh, ekv, egp, eto !< local fields for energy densities
    REAL(wp), DIMENSION(nproma) :: einvi, ekhvi, ekvvi, egpvi, etovi    !< local fields for energy contents

    !$ACC DATA CREATE(ein, ekh, ekv, egp, eto, einvi, ekhvi, ekvvi, egpvi, etovi)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)

    !> Compute energy density (J/m3) and energy content (J/m2)
    !> - Input : from the physics memory
    !> - Output: to the "1" state of energy memory
    !
    jks = 1
    jke = num_lev(jg)
    !
    lein   = ASSOCIATED(e(jg)%ein1)
    lekh   = ASSOCIATED(e(jg)%ekh1)
    lekv   = ASSOCIATED(e(jg)%ekv1)
    legp   = ASSOCIATED(e(jg)%egp1)
    leto   = ASSOCIATED(e(jg)%eto1)
    !
    leinvi = ASSOCIATED(e(jg)%ein1vi)
    lekhvi = ASSOCIATED(e(jg)%ekh1vi)
    lekvvi = ASSOCIATED(e(jg)%ekv1vi)
    legpvi = ASSOCIATED(e(jg)%egp1vi)
    letovi = ASSOCIATED(e(jg)%eto1vi)
    !
    CALL atm_energy_diag(jcs, jce,                                    & !< domain indices
         !
         &               s(jg)%metrics%ddqz_z_full(:,:,jb),           & !< layer thickness
         &               f(jg)%rho(:,:,jb),                           & !< mass density
         &               f(jg)%ta(:,:,jb),                            & !< air temperature
         &               f(jg)%qtrc_phy(:,:,jb,iqv),                  & !< mass fraction in air of water vapor
         &               f(jg)%qtrc_phy(:,:,jb,iqc),                  & !< mass fraction in air of cloud liquid
         &               f(jg)%qtrc_phy(:,:,jb,iqi),                  & !< mass fraction in air of cloud ice
         &               f(jg)%qtrc_phy(:,:,jb,iqr),                  & !< mass fraction in air of rain
         &               f(jg)%qtrc_phy(:,:,jb,iqs),                  & !< mass fraction in air of snow
         &               f(jg)%qtrc_phy(:,:,jb,iqg),                  & !< mass fraction in air of graupel
         &               f(jg)%ua(:,:,jb),                            & !< zonal velocity
         &               f(jg)%va(:,:,jb),                            & !< meridional velocity
         &               f(jg)%wa(:,:,jb),                            & !< vertical velocity
         &               s(jg)%metrics%geopot(:,:,jb),                & !< geopotential
         !
         &               lein, ein(:,:),                              & !< internal energy density
         &               lekh, ekh(:,:),                              & !< horizontal kinetic energy density
         &               lekv, ekv(:,:),                              & !< vertical kinetic energy density
         &               legp, egp(:,:),                              & !< geopotential energy density
         &               leto, eto(:,:),                              & !< total energy density
         !
         &               leinvi, einvi(:),                            & !< internal energy content
         &               lekhvi, ekhvi(:),                            & !< horizontal kinetic energy content
         &               lekvvi, ekvvi(:),                            & !< vertical kinetic energy content
         &               legpvi, egpvi(:),                            & !< geopotential energy content
         &               letovi, etovi(:))                              !< total energy content
    !
    !> Store local fields in global fields, where the latter exist
    !
    IF (lein)   CALL copy(jcs, jce, jks, jke, ein(:,:), e(jg)%ein1(:,:,jb))
    IF (lekh)   CALL copy(jcs, jce, jks, jke, ekh(:,:), e(jg)%ekh1(:,:,jb))
    IF (lekv)   CALL copy(jcs, jce, jks, jke, ekv(:,:), e(jg)%ekv1(:,:,jb))
    IF (legp)   CALL copy(jcs, jce, jks, jke, egp(:,:), e(jg)%egp1(:,:,jb))
    IF (leto)   CALL copy(jcs, jce, jks, jke, eto(:,:), e(jg)%eto1(:,:,jb))
    !
    IF (leinvi) CALL copy(jcs, jce, einvi(:), e(jg)%ein1vi(:,jb))
    IF (lekhvi) CALL copy(jcs, jce, ekhvi(:), e(jg)%ekh1vi(:,jb))
    IF (lekvvi) CALL copy(jcs, jce, ekvvi(:), e(jg)%ekv1vi(:,jb))
    IF (legpvi) CALL copy(jcs, jce, egpvi(:), e(jg)%egp1vi(:,jb))
    IF (letovi) CALL copy(jcs, jce, etovi(:), e(jg)%eto1vi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

    !$ACC END DATA

  END SUBROUTINE atm_energy_diag_p1

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_diag_p2(jg, jb, jcs, jce)

    INTEGER, INTENT(in)      :: jg, jb, jcs, jce                        !< horizontal domain indices

    INTEGER                  :: jks, jke                                !< vertical domain indices

    LOGICAL                  :: lein,   lekh,   lekv,   legp,   leto    !< logicals for computation of energy densities
    LOGICAL                  :: leinvi, lekhvi, lekvvi, legpvi, letovi  !< logicals for computation of energy contents

    REAL(wp), DIMENSION(nproma, num_lev(jg)) :: ein, ekh, ekv, egp, eto !< local fields for energy densities
    REAL(wp), DIMENSION(nproma) :: einvi, ekhvi, ekvvi, egpvi, etovi    !< local fields for energy contents

    !$ACC DATA CREATE(ein, ekh, ekv, egp, eto, einvi, ekhvi, ekvvi, egpvi, etovi)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)

    !> Compute energy density (J/m3) and energy content (J/m2)
    !> - Input : from the physics memory
    !> - Output: to the "2" state of energy memory
    !
    jks = 1
    jke = num_lev(jg)
    !
    lein   = ASSOCIATED(e(jg)%ein2)
    lekh   = ASSOCIATED(e(jg)%ekh2)
    lekv   = ASSOCIATED(e(jg)%ekv2)
    legp   = ASSOCIATED(e(jg)%egp2)
    leto   = ASSOCIATED(e(jg)%eto2)
    !
    leinvi = ASSOCIATED(e(jg)%ein2vi)
    lekhvi = ASSOCIATED(e(jg)%ekh2vi)
    lekvvi = ASSOCIATED(e(jg)%ekv2vi)
    legpvi = ASSOCIATED(e(jg)%egp2vi)
    letovi = ASSOCIATED(e(jg)%eto2vi)
    !
    CALL atm_energy_diag(jcs, jce,                                    & !< domain indices
         !
         &               s(jg)%metrics%ddqz_z_full(:,:,jb),           & !< layer thickness
         &               f(jg)%rho(:,:,jb),                           & !< mass density
         &               f(jg)%ta(:,:,jb),                            & !< air temperature
         &               f(jg)%qtrc_phy(:,:,jb,iqv),                  & !< mass fraction in air of water vapor
         &               f(jg)%qtrc_phy(:,:,jb,iqc),                  & !< mass fraction in air of cloud liquid
         &               f(jg)%qtrc_phy(:,:,jb,iqi),                  & !< mass fraction in air of cloud ice
         &               f(jg)%qtrc_phy(:,:,jb,iqr),                  & !< mass fraction in air of rain
         &               f(jg)%qtrc_phy(:,:,jb,iqs),                  & !< mass fraction in air of snow
         &               f(jg)%qtrc_phy(:,:,jb,iqg),                  & !< mass fraction in air of graupel
         &               f(jg)%ua(:,:,jb),                            & !< zonal velocity
         &               f(jg)%va(:,:,jb),                            & !< meridional velocity
         &               f(jg)%wa(:,:,jb),                            & !< vertical velocity
         &               s(jg)%metrics%geopot(:,:,jb),                & !< geopotential
         !
         &               lein, ein(:,:),                              & !< internal energy density
         &               lekh, ekh(:,:),                              & !< horizontal kinetic energy density
         &               lekv, ekv(:,:),                              & !< vertical kinetic energy density
         &               legp, egp(:,:),                              & !< geopotential energy density
         &               leto, eto(:,:),                              & !< total energy density
         !
         &               leinvi, einvi(:),                            & !< internal energy content
         &               lekhvi, ekhvi(:),                            & !< horizontal kinetic energy content
         &               lekvvi, ekvvi(:),                            & !< vertical kinetic energy content
         &               legpvi, egpvi(:),                            & !< geopotential energy content
         &               letovi, etovi(:))                              !< total energy content

    !> Store local fields in global fields, where the latter exist
    !
    IF (lein)   CALL copy(jcs, jce, jks, jke, ein(:,:), e(jg)%ein2(:,:,jb))
    IF (lekh)   CALL copy(jcs, jce, jks, jke, ekh(:,:), e(jg)%ekh2(:,:,jb))
    IF (lekv)   CALL copy(jcs, jce, jks, jke, ekv(:,:), e(jg)%ekv2(:,:,jb))
    IF (legp)   CALL copy(jcs, jce, jks, jke, egp(:,:), e(jg)%egp2(:,:,jb))
    IF (leto)   CALL copy(jcs, jce, jks, jke, eto(:,:), e(jg)%eto2(:,:,jb))
    !
    IF (leinvi) CALL copy(jcs, jce, einvi(:), e(jg)%ein2vi(:,jb))
    IF (lekhvi) CALL copy(jcs, jce, ekhvi(:), e(jg)%ekh2vi(:,jb))
    IF (lekvvi) CALL copy(jcs, jce, ekvvi(:), e(jg)%ekv2vi(:,jb))
    IF (legpvi) CALL copy(jcs, jce, egpvi(:), e(jg)%egp2vi(:,jb))
    IF (letovi) CALL copy(jcs, jce, etovi(:), e(jg)%eto2vi(:,jb))

    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

    !$ACC END DATA

  END SUBROUTINE atm_energy_diag_p2

  !-------------------------------------------------------------------

  PURE SUBROUTINE atm_energy_diag(jcs, jce,                           & !< horizontal domain indices
       !
       &                          dz,                                 & !< layer thickness
       &                          rho,                                & !< mass density
       &                          temp,                               & !< air temperature
       &                          qv,                                 & !< mass fraction in air of water vapor
       &                          qc,                                 & !< mass fraction in air of cloud liquid
       &                          qi,                                 & !< mass fraction in air of cloud ice
       &                          qr,                                 & !< mass fraction in air of rain
       &                          qs,                                 & !< mass fraction in air of snow
       &                          qg,                                 & !< mass fraction in air of graupel
       &                          u,                                  & !< zonal velocity
       &                          v,                                  & !< meridional velocity
       &                          w,                                  & !< vertical velocity
       &                          geopot,                             & !< geopotential
       !
       &                          lein, ein,                          & !< internal energy density
       &                          lekh, ekh,                          & !< horizontal kinetic energy density
       &                          lekv, ekv,                          & !< vertical kinetic energy density
       &                          legp, egp,                          & !< geopotential energy density
       &                          leto, eto,                          & !< total energy density
       !
       &                          leinvi, einvi,                      & !< internal energy content
       &                          lekhvi, ekhvi,                      & !< horizontal kinetic energy content
       &                          lekvvi, ekvvi,                      & !< vertical kinetic energy content
       &                          legpvi, egpvi,                      & !< geopotential energy content
       &                          letovi, etovi)                        !< total energy content


    INTEGER,  INTENT(in)  :: jcs, jce

    REAL(wp), INTENT(in)  ::     dz(:,:)
    REAL(wp), INTENT(in)  ::    rho(:,:)
    REAL(wp), INTENT(in)  ::   temp(:,:)
    REAL(wp), INTENT(in)  ::     qv(:,:)
    REAL(wp), INTENT(in)  ::     qc(:,:)
    REAL(wp), INTENT(in)  ::     qi(:,:)
    REAL(wp), INTENT(in)  ::     qr(:,:)
    REAL(wp), INTENT(in)  ::     qs(:,:)
    REAL(wp), INTENT(in)  ::     qg(:,:)
    REAL(wp), INTENT(in)  ::      u(:,:)
    REAL(wp), INTENT(in)  ::      v(:,:)
    REAL(wp), INTENT(in)  ::      w(:,:)
    REAL(wp), INTENT(in)  :: geopot(:,:)

    LOGICAL,  INTENT(in)  ::    lein,   lekh,   lekv,   legp,   leto
    LOGICAL,  INTENT(in)  ::    leinvi, lekhvi, lekvvi, legpvi, letovi

    REAL(wp), INTENT(out) ::    ein(:,:)
    REAL(wp), INTENT(out) ::    ekh(:,:)
    REAL(wp), INTENT(out) ::    ekv(:,:)
    REAL(wp), INTENT(out) ::    egp(:,:)
    REAL(wp), INTENT(out) ::    eto(:,:)

    REAL(wp), INTENT(out) ::    einvi(:)
    REAL(wp), INTENT(out) ::    ekhvi(:)
    REAL(wp), INTENT(out) ::    ekvvi(:)
    REAL(wp), INTENT(out) ::    egpvi(:)
    REAL(wp), INTENT(out) ::    etovi(:)

    INTEGER               :: jks, jke                                   !< vertical domain indices

    jks = 1
    jke = SIZE(rho,2)

    !> Compute internal, horizontal and vertical kinetic, geopotential and total energy density (J/m3)
    !
    IF (lein .OR. leinvi .OR. leto .OR. letovi) THEN
       ein = atm_energy_diag_ein(                                        & !< internal energy density
            &                    jcs, jce, jks, jke,                     & !< domain indices
            &                    rho,                                    & !< mass density
            &                    temp,                                   & !< air temperature
            &                    qv,                                     & !< mass fraction in air of water vapor
            &                    qc,                                     & !< mass fraction in air of cloud liquid
            &                    qi,                                     & !< mass fraction in air of cloud ice
            &                    qr,                                     & !< mass fraction in air of rain
            &                    qs,                                     & !< mass fraction in air of snow
            &                    qg)                                       !< mass fraction in air of graupel
    END IF
    !
    IF (lekh .OR. lekhvi .OR. leto .OR. letovi) THEN
       ekh = atm_energy_diag_ekh(                                        & !< horizontal kinetic energy density
            &                    jcs, jce, jks, jke,                     & !< domain indices
            &                    rho,                                    & !< mass density
            &                    u,                                      & !< zonal velocity
            &                    v)                                        !< meridional velocity
    END IF
    !
    IF (lekv .OR. lekvvi .OR. leto .OR. letovi) THEN
       ekv = atm_energy_diag_ekv(                                        & !< vertical kinetic energy density
            &                    jcs, jce, jks, jke,                     & !< domain indices
            &                    rho,                                    & !< mass density
            &                    w)                                        !< vertical velocity
    END IF
    !
    IF (legp .OR. legpvi .OR. leto .OR. letovi) THEN
       egp = atm_energy_diag_egp(                                        & !< geopotential energy density
            &                    jcs, jce, jks, jke,                     & !< domain indices
            &                    rho,                                    & !< mass density
            &                    geopot)                                   !< geopotential
    END IF
    !
    IF (leto .OR. letovi) THEN
       eto = atm_energy_diag_eto(                                        & !< total energy density
            &                    jcs, jce, jks, jke,                     & !< domain indices
            &                    ein,                                    & !< internal energy density
            &                    ekh,                                    & !< horizontal kinetic energy density
            &                    ekv,                                    & !< vertical kinetic energy density
            &                    egp)                                      !< geopotential energy density
    END IF
    
    !> Compute vertical integrals: energy density (J/m3) -> energy content (J/m2)
    !
    IF (leinvi) THEN
       einvi = atm_energy_vint(                                          & !< internal energy content
            &                  jcs, jce, jks, jke,                       & !< domain indices
            &                  dz,                                       & !< layer thickness
            &                  ein)                                        !< internal energy density
    END IF
    !
    IF (lekhvi) THEN
       ekhvi = atm_energy_vint(                                          & !< horizontal kinetic energy content
            &                  jcs, jce, jks, jke,                       & !< domain indices
            &                  dz,                                       & !< layer thickness
            &                  ekh)                                        !< horizontal kinetic energy density
    END IF
    !
    IF (lekvvi) THEN
       ekvvi = atm_energy_vint(                                          & !< vertical kinetic energy content
            &                  jcs, jce, jks, jke,                       & !< domain indices
            &                  dz,                                       & !< layer thickness
            &                  ekv)                                        !< vertical kinetic energy density
    END IF
    !
    IF (legpvi) THEN
       egpvi = atm_energy_vint(                                          & !< geopotential energy content
            &                  jcs, jce, jks, jke,                       & !< domain indices
            &                  dz,                                       & !< layer thickness
            &                  egp)                                        !< geopotential energy density
    END IF
    !
    IF (letovi) THEN
       etovi = atm_energy_vint(                                          & !< total energy content
            &                  jcs, jce, jks, jke,                       & !< domain indices
            &                  dz,                                       & !< layer thickness
            &                  eto)                                        !< total energy density
    END IF

  END SUBROUTINE atm_energy_diag

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_dyn_3d_vi(jg, jb, jcs, jce)
    
    INTEGER, INTENT(in)             :: jg, jb, jcs, jce                 !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy density
    !
    IF (ASSOCIATED(e(jg)%eindyn)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ein1(:,:,jb), e(jg)%ein2(:,:,jb), e(jg)%eindyn(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekhdyn)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ekh1(:,:,jb), e(jg)%ekh2(:,:,jb), e(jg)%ekhdyn(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekvdyn)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ekv1(:,:,jb), e(jg)%ekv2(:,:,jb), e(jg)%ekvdyn(:,:,jb))
    IF (ASSOCIATED(e(jg)%egpdyn)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%egp1(:,:,jb), e(jg)%egp2(:,:,jb), e(jg)%egpdyn(:,:,jb))
    IF (ASSOCIATED(e(jg)%etodyn)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%eto1(:,:,jb), e(jg)%eto2(:,:,jb), e(jg)%etodyn(:,:,jb))
    !
    !> Tendency of energy content
    !
    IF (ASSOCIATED(e(jg)%eindynvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ein1vi(:,jb), e(jg)%ein2vi(:,jb), e(jg)%eindynvi(:,jb))
    IF (ASSOCIATED(e(jg)%ekhdynvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ekh1vi(:,jb), e(jg)%ekh2vi(:,jb), e(jg)%ekhdynvi(:,jb))
    IF (ASSOCIATED(e(jg)%ekvdynvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ekv1vi(:,jb), e(jg)%ekv2vi(:,jb), e(jg)%ekvdynvi(:,jb))
    IF (ASSOCIATED(e(jg)%egpdynvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%egp1vi(:,jb), e(jg)%egp2vi(:,jb), e(jg)%egpdynvi(:,jb))
    IF (ASSOCIATED(e(jg)%etodynvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%eto1vi(:,jb), e(jg)%eto2vi(:,jb), e(jg)%etodynvi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_dyn_3d_vi

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_dyn_hi_ti(jg)
    
    INTEGER, INTENT(in)             :: jg                               !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy profile
    !
    IF (ASSOCIATED(e(jg)%eindynhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ein1hi(:), e(jg)%ein2hi(:), e(jg)%eindynhi(:))
    IF (ASSOCIATED(e(jg)%ekhdynhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ekh1hi(:), e(jg)%ekh2hi(:), e(jg)%ekhdynhi(:))
    IF (ASSOCIATED(e(jg)%ekvdynhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ekv1hi(:), e(jg)%ekv2hi(:), e(jg)%ekvdynhi(:))
    IF (ASSOCIATED(e(jg)%egpdynhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%egp1hi(:), e(jg)%egp2hi(:), e(jg)%egpdynhi(:))
    IF (ASSOCIATED(e(jg)%etodynhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%eto1hi(:), e(jg)%eto2hi(:), e(jg)%etodynhi(:))
    !
    !> Tendency of energy
    !
    IF (ASSOCIATED(e(jg)%eindynti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ein1ti(:), e(jg)%ein2ti(:), e(jg)%eindynti(:))
    IF (ASSOCIATED(e(jg)%ekhdynti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ekh1ti(:), e(jg)%ekh2ti(:), e(jg)%ekhdynti(:))
    IF (ASSOCIATED(e(jg)%ekvdynti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ekv1ti(:), e(jg)%ekv2ti(:), e(jg)%ekvdynti(:))
    IF (ASSOCIATED(e(jg)%egpdynti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%egp1ti(:), e(jg)%egp2ti(:), e(jg)%egpdynti(:))
    IF (ASSOCIATED(e(jg)%etodynti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%eto1ti(:), e(jg)%eto2ti(:), e(jg)%etodynti(:))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_dyn_hi_ti

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_phy_3d_vi(jg, jb, jcs, jce)
    
    INTEGER, INTENT(in)             :: jg, jb, jcs, jce                 !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy density
    !
    IF (ASSOCIATED(e(jg)%einphy)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ein3(:,:,jb), e(jg)%ein2(:,:,jb), e(jg)%einphy(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekhphy)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ekh3(:,:,jb), e(jg)%ekh2(:,:,jb), e(jg)%ekhphy(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekvphy)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ekv3(:,:,jb), e(jg)%ekv2(:,:,jb), e(jg)%ekvphy(:,:,jb))
    IF (ASSOCIATED(e(jg)%egpphy)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%egp3(:,:,jb), e(jg)%egp2(:,:,jb), e(jg)%egpphy(:,:,jb))
    IF (ASSOCIATED(e(jg)%etophy)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%eto3(:,:,jb), e(jg)%eto2(:,:,jb), e(jg)%etophy(:,:,jb))
    !
    !> Tendency of energy content
    !
    IF (ASSOCIATED(e(jg)%einphyvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ein3vi(:,jb), e(jg)%ein2vi(:,jb), e(jg)%einphyvi(:,jb))
    IF (ASSOCIATED(e(jg)%ekhphyvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ekh3vi(:,jb), e(jg)%ekh2vi(:,jb), e(jg)%ekhphyvi(:,jb))
    IF (ASSOCIATED(e(jg)%ekvphyvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ekv3vi(:,jb), e(jg)%ekv2vi(:,jb), e(jg)%ekvphyvi(:,jb))
    IF (ASSOCIATED(e(jg)%egpphyvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%egp3vi(:,jb), e(jg)%egp2vi(:,jb), e(jg)%egpphyvi(:,jb))
    IF (ASSOCIATED(e(jg)%etophyvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%eto3vi(:,jb), e(jg)%eto2vi(:,jb), e(jg)%etophyvi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_phy_3d_vi

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_phy_hi_ti(jg)

    INTEGER, INTENT(in)             :: jg                               !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy profile
    !
    IF (ASSOCIATED(e(jg)%einphyhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ein3hi(:), e(jg)%ein2hi(:), e(jg)%einphyhi(:))
    IF (ASSOCIATED(e(jg)%ekhphyhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ekh3hi(:), e(jg)%ekh2hi(:), e(jg)%ekhphyhi(:))
    IF (ASSOCIATED(e(jg)%ekvphyhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ekv3hi(:), e(jg)%ekv2hi(:), e(jg)%ekvphyhi(:))
    IF (ASSOCIATED(e(jg)%egpphyhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%egp3hi(:), e(jg)%egp2hi(:), e(jg)%egpphyhi(:))
    IF (ASSOCIATED(e(jg)%etophyhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%eto3hi(:), e(jg)%eto2hi(:), e(jg)%etophyhi(:))
    !
    !> Tendency of energy
    !
    IF (ASSOCIATED(e(jg)%einphyti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ein3ti(:), e(jg)%ein2ti(:), e(jg)%einphyti(:))
    IF (ASSOCIATED(e(jg)%ekhphyti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ekh3ti(:), e(jg)%ekh2ti(:), e(jg)%ekhphyti(:))
    IF (ASSOCIATED(e(jg)%ekvphyti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ekv3ti(:), e(jg)%ekv2ti(:), e(jg)%ekvphyti(:))
    IF (ASSOCIATED(e(jg)%egpphyti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%egp3ti(:), e(jg)%egp2ti(:), e(jg)%egpphyti(:))
    IF (ASSOCIATED(e(jg)%etophyti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%eto3ti(:), e(jg)%eto2ti(:), e(jg)%etophyti(:))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_phy_hi_ti

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_cld_3d_vi(jg, jb, jcs, jce)
    
    INTEGER, INTENT(in)             :: jg, jb, jcs, jce                !< horizontal domain indices

    INTEGER                         :: jks, jke                        !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy density
    !
    IF (ASSOCIATED(e(jg)%eincld)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ein1(:,:,jb), e(jg)%ein2(:,:,jb), e(jg)%eincld(:,:,jb))
    !
    !> Tendency of energy content
    !
    IF (ASSOCIATED(e(jg)%eincldvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ein1vi(:,jb), e(jg)%ein2vi(:,jb), e(jg)%eincldvi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_cld_3d_vi

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_cld_hi_ti(jg)
    
    INTEGER, INTENT(in)             :: jg                               !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy profile
    !
    IF (ASSOCIATED(e(jg)%eincldhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ein1hi(:), e(jg)%ein2hi(:), e(jg)%eincldhi(:))
    !
    !> Tendency of energy
    !
    IF (ASSOCIATED(e(jg)%eincldti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ein1ti(:), e(jg)%ein2ti(:), e(jg)%eincldti(:))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_cld_hi_ti

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_rad_3d_vi(jg, jb, jcs, jce)
    
    INTEGER, INTENT(in)             :: jg, jb, jcs, jce                 !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy density
    !
    IF (ASSOCIATED(e(jg)%einrad)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ein1(:,:,jb), e(jg)%ein2(:,:,jb), e(jg)%einrad(:,:,jb))
    !
    !> Tendency of energy content
    !
    IF (ASSOCIATED(e(jg)%einradvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ein1vi(:,jb), e(jg)%ein2vi(:,jb), e(jg)%einradvi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_rad_3d_vi

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_rad_hi_ti(jg)
    
    INTEGER, INTENT(in)             :: jg                               !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy profile
    !
    IF (ASSOCIATED(e(jg)%einradhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ein1hi(:), e(jg)%ein2hi(:), e(jg)%einradhi(:))
    !
    !> Tendency of energy
    !
    IF (ASSOCIATED(e(jg)%einradti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ein1ti(:), e(jg)%ein2ti(:), e(jg)%einradti(:))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_rad_hi_ti

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_tmx_3d_vi(jg, jb, jcs, jce)
    
    INTEGER, INTENT(in)             :: jg, jb, jcs, jce                 !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy density
    !
    IF (ASSOCIATED(e(jg)%eintmx)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ein1(:,:,jb), e(jg)%ein2(:,:,jb), e(jg)%eintmx(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekhtmx)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ekh1(:,:,jb), e(jg)%ekh2(:,:,jb), e(jg)%ekhtmx(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekvtmx)) CALL atm_energy_tend_2d(jcs, jce, jks, jke, dtime, e(jg)%ekv1(:,:,jb), e(jg)%ekv2(:,:,jb), e(jg)%ekvtmx(:,:,jb))
    !
    !> Tendency of energy content
    !
    IF (ASSOCIATED(e(jg)%eintmxvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ein1vi(:,jb), e(jg)%ein2vi(:,jb), e(jg)%eintmxvi(:,jb))
    IF (ASSOCIATED(e(jg)%ekhtmxvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ekh1vi(:,jb), e(jg)%ekh2vi(:,jb), e(jg)%ekhtmxvi(:,jb))
    IF (ASSOCIATED(e(jg)%ekvtmxvi)) CALL atm_energy_tend_1d(jcs, jce, dtime, e(jg)%ekv1vi(:,jb), e(jg)%ekv2vi(:,jb), e(jg)%ekvtmxvi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_tmx_3d_vi

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_tend_tmx_hi_ti(jg)
    
    INTEGER, INTENT(in)             :: jg                               !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Tendency of energy profile
    !
    IF (ASSOCIATED(e(jg)%eintmxhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ein1hi(:), e(jg)%ein2hi(:), e(jg)%eintmxhi(:))
    IF (ASSOCIATED(e(jg)%ekhtmxhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ekh1hi(:), e(jg)%ekh2hi(:), e(jg)%ekhtmxhi(:))
    IF (ASSOCIATED(e(jg)%ekvtmxhi)) CALL atm_energy_tend_1d(jks, jke, dtime, e(jg)%ekv1hi(:), e(jg)%ekv2hi(:), e(jg)%ekvtmxhi(:))
    !
    !> Tendency of energy
    !
    IF (ASSOCIATED(e(jg)%eintmxti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ein1ti(:), e(jg)%ein2ti(:), e(jg)%eintmxti(:))
    IF (ASSOCIATED(e(jg)%ekhtmxti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ekh1ti(:), e(jg)%ekh2ti(:), e(jg)%ekhtmxti(:))
    IF (ASSOCIATED(e(jg)%ekvtmxti)) CALL atm_energy_tend_1d(  1,   1, dtime, e(jg)%ekv1ti(:), e(jg)%ekv2ti(:), e(jg)%ekvtmxti(:))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_tend_tmx_hi_ti

  !-------------------------------------------------------------------

  PURE SUBROUTINE atm_energy_tend_2d(js, je, ks, ke, &
       &                             dtime,          &
       &                             x1,             &
       &                             x2,             &
       &                             tend_x_2d)

    INTEGER , INTENT(in)          :: js, je, ks, ke                     !< domain indices
    REAL(wp), INTENT(in)          :: dtime                              !< (s)     time step
    REAL(wp), INTENT(in)          :: x1(:,:)                            !< (<x>)   x at stage 1
    REAL(wp), INTENT(in)          :: x2(:,:)                            !< (<x>)   x at stage 2
    REAL(wp), INTENT(out)         :: tend_x_2d(SIZE(x1,1),SIZE(x1,2))   !< (<x>/s) tendency of x

    INTEGER                       :: j, k                               !< loop indices
    REAL(wp)                      :: rdtime                             !< 1/dtime

    rdtime = 1.0_wp/dtime

    !$ACC DATA PRESENT(x1, x2, tend_x_2d)
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(js, je, ks, ke, rdtime) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO k = ks, ke
       DO j = js, je

          tend_x_2d(j,k) = (x2(j,k) - x1(j,k))*rdtime                   !< (<x>/s) tendency of x

       END DO
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE atm_energy_tend_2d

  !-------------------------------------------------------------------

  PURE SUBROUTINE atm_energy_tend_1d(js, je, &
       &                             dtime,  &
       &                             x1,     &
       &                             x2,     &
       &                             tend_x_1d)

    INTEGER , INTENT(in)          :: js, je                             !< domain indices
    REAL(wp), INTENT(in)          :: dtime                              !< (s)        time step
    REAL(wp), INTENT(in)          :: x1(:)                              !< (<x>)   x at stage 1
    REAL(wp), INTENT(in)          :: x2(:)                              !< (<x>)   x at stage 2
    REAL(wp), INTENT(out)         :: tend_x_1d(SIZE(x1,1))              !< (<x>/s) tendency of x

    INTEGER                       :: j                                  !< loop indices
    REAL(wp)                      :: rdtime                             !< 1/dtime

    rdtime = 1.0_wp/dtime

    !$ACC DATA PRESENT(x1, x2, tend_x_1d)
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(js, je, rdtime) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO j = js, je

       tend_x_1d(j) = (x2(j) - x1(j))*rdtime                            !< (<x>/s)  tendency of x

    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE atm_energy_tend_1d

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_copy_2_1_3d_vi(jg, jb, jcs, jce)

    INTEGER, INTENT(in)             :: jg, jb, jcs, jce                 !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Copy energy density
    !
    IF (ASSOCIATED(e(jg)%ein1))   CALL copy(jcs, jce, jks, jke, e(jg)%ein2(:,:,jb), e(jg)%ein1(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekh1))   CALL copy(jcs, jce, jks, jke, e(jg)%ekh2(:,:,jb), e(jg)%ekh1(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekv1))   CALL copy(jcs, jce, jks, jke, e(jg)%ekv2(:,:,jb), e(jg)%ekv1(:,:,jb))
    IF (ASSOCIATED(e(jg)%egp1))   CALL copy(jcs, jce, jks, jke, e(jg)%egp2(:,:,jb), e(jg)%egp1(:,:,jb))
    IF (ASSOCIATED(e(jg)%eto1))   CALL copy(jcs, jce, jks, jke, e(jg)%eto2(:,:,jb), e(jg)%eto1(:,:,jb))
    !
    !> Copy energy content
    !
    IF (ASSOCIATED(e(jg)%ein1vi)) CALL copy(jcs, jce, e(jg)%ein2vi(:,jb), e(jg)%ein1vi(:,jb))
    IF (ASSOCIATED(e(jg)%ekh1vi)) CALL copy(jcs, jce, e(jg)%ekh2vi(:,jb), e(jg)%ekh1vi(:,jb))
    IF (ASSOCIATED(e(jg)%ekv1vi)) CALL copy(jcs, jce, e(jg)%ekv2vi(:,jb), e(jg)%ekv1vi(:,jb))
    IF (ASSOCIATED(e(jg)%egp1vi)) CALL copy(jcs, jce, e(jg)%egp2vi(:,jb), e(jg)%egp1vi(:,jb))
    IF (ASSOCIATED(e(jg)%eto1vi)) CALL copy(jcs, jce, e(jg)%eto2vi(:,jb), e(jg)%eto1vi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_copy_2_1_3d_vi

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_copy_2_1_hi_ti(jg)
    
    INTEGER, INTENT(in)             :: jg                               !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Copy energy profile
    !
    IF (ASSOCIATED(e(jg)%ein1hi)) CALL copy(jks, jke, e(jg)%ein2hi(:), e(jg)%ein1hi(:))
    IF (ASSOCIATED(e(jg)%ekh1hi)) CALL copy(jks, jke, e(jg)%ekh2hi(:), e(jg)%ekh1hi(:))
    IF (ASSOCIATED(e(jg)%ekv1hi)) CALL copy(jks, jke, e(jg)%ekv2hi(:), e(jg)%ekv1hi(:))
    IF (ASSOCIATED(e(jg)%egp1hi)) CALL copy(jks, jke, e(jg)%egp2hi(:), e(jg)%egp1hi(:))
    IF (ASSOCIATED(e(jg)%eto1hi)) CALL copy(jks, jke, e(jg)%eto2hi(:), e(jg)%eto1hi(:))
    !
    !> Copy energy
    !
    IF (ASSOCIATED(e(jg)%ein1ti)) CALL copy(  1,   1, e(jg)%ein2ti(:), e(jg)%ein1ti(:))
    IF (ASSOCIATED(e(jg)%ekh1ti)) CALL copy(  1,   1, e(jg)%ekh2ti(:), e(jg)%ekh1ti(:))
    IF (ASSOCIATED(e(jg)%ekv1ti)) CALL copy(  1,   1, e(jg)%ekv2ti(:), e(jg)%ekv1ti(:))
    IF (ASSOCIATED(e(jg)%egp1ti)) CALL copy(  1,   1, e(jg)%egp2ti(:), e(jg)%egp1ti(:))
    IF (ASSOCIATED(e(jg)%eto1ti)) CALL copy(  1,   1, e(jg)%eto2ti(:), e(jg)%eto1ti(:))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_copy_2_1_hi_ti

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_copy_2_3_3d_vi(jg, jb, jcs, jce)

    INTEGER, INTENT(in)             :: jg, jb, jcs, jce                 !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Copy energy density
    !
    IF (ASSOCIATED(e(jg)%ein3))   CALL copy(jcs, jce, jks, jke, e(jg)%ein2(:,:,jb), e(jg)%ein3(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekh3))   CALL copy(jcs, jce, jks, jke, e(jg)%ekh2(:,:,jb), e(jg)%ekh3(:,:,jb))
    IF (ASSOCIATED(e(jg)%ekv3))   CALL copy(jcs, jce, jks, jke, e(jg)%ekv2(:,:,jb), e(jg)%ekv3(:,:,jb))
    IF (ASSOCIATED(e(jg)%egp3))   CALL copy(jcs, jce, jks, jke, e(jg)%egp2(:,:,jb), e(jg)%egp3(:,:,jb))
    IF (ASSOCIATED(e(jg)%eto3))   CALL copy(jcs, jce, jks, jke, e(jg)%eto2(:,:,jb), e(jg)%eto3(:,:,jb))
    !
    !> Copy energy content
    !
    IF (ASSOCIATED(e(jg)%ein3vi)) CALL copy(jcs, jce, e(jg)%ein2vi(:,jb), e(jg)%ein3vi(:,jb))
    IF (ASSOCIATED(e(jg)%ekh3vi)) CALL copy(jcs, jce, e(jg)%ekh2vi(:,jb), e(jg)%ekh3vi(:,jb))
    IF (ASSOCIATED(e(jg)%ekv3vi)) CALL copy(jcs, jce, e(jg)%ekv2vi(:,jb), e(jg)%ekv3vi(:,jb))
    IF (ASSOCIATED(e(jg)%egp3vi)) CALL copy(jcs, jce, e(jg)%egp2vi(:,jb), e(jg)%egp3vi(:,jb))
    IF (ASSOCIATED(e(jg)%eto3vi)) CALL copy(jcs, jce, e(jg)%eto2vi(:,jb), e(jg)%eto3vi(:,jb))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_copy_2_3_3d_vi

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_copy_2_3_hi_ti(jg)
    
    INTEGER, INTENT(in)             :: jg                               !< horizontal domain indices

    INTEGER                         :: jks, jke                         !< vertical domain indices

    jks = 1
    jke = num_lev(jg)

    IF (ltimer) CALL timer_start(timer_atm_energy_diag)
    !
    !> Copy energy profile
    !
    IF (ASSOCIATED(e(jg)%ein3hi)) CALL copy(jks, jke, e(jg)%ein2hi(:), e(jg)%ein3hi(:))
    IF (ASSOCIATED(e(jg)%ekh3hi)) CALL copy(jks, jke, e(jg)%ekh2hi(:), e(jg)%ekh3hi(:))
    IF (ASSOCIATED(e(jg)%ekv3hi)) CALL copy(jks, jke, e(jg)%ekv2hi(:), e(jg)%ekv3hi(:))
    IF (ASSOCIATED(e(jg)%egp3hi)) CALL copy(jks, jke, e(jg)%egp2hi(:), e(jg)%egp3hi(:))
    IF (ASSOCIATED(e(jg)%eto3hi)) CALL copy(jks, jke, e(jg)%eto2hi(:), e(jg)%eto3hi(:))
    !
    !> Copy energy
    !
    IF (ASSOCIATED(e(jg)%ein3ti)) CALL copy(  1,   1, e(jg)%ein2ti(:), e(jg)%ein3ti(:))
    IF (ASSOCIATED(e(jg)%ekh3ti)) CALL copy(  1,   1, e(jg)%ekh2ti(:), e(jg)%ekh3ti(:))
    IF (ASSOCIATED(e(jg)%ekv3ti)) CALL copy(  1,   1, e(jg)%ekv2ti(:), e(jg)%ekv3ti(:))
    IF (ASSOCIATED(e(jg)%egp3ti)) CALL copy(  1,   1, e(jg)%egp2ti(:), e(jg)%egp3ti(:))
    IF (ASSOCIATED(e(jg)%eto3ti)) CALL copy(  1,   1, e(jg)%eto2ti(:), e(jg)%eto3ti(:))
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_diag)

  END SUBROUTINE atm_energy_copy_2_3_hi_ti

  !-------------------------------------------------------------------

  PURE FUNCTION atm_energy_diag_ein(jcs, jce, jks, jke, &
       &                            rho,                &
       &                            ta,                 &
       &                            qv,                 &
       &                            qc,                 &
       &                            qi,                 &
       &                            qr,                 &
       &                            qs,                 &
       &                            qg)                 &
       &                     RESULT(ein)

    INTEGER , INTENT(in)         :: jcs, jce, jks, jke                  !< domain indices

    REAL(wp), INTENT(in)         :: rho(:,:)                            !< (kg/m3)  mass density
    REAL(wp), INTENT(in)         :: ta(:,:)                             !< (K)      temperature of air
    REAL(wp), INTENT(in)         :: qv(:,:)                             !< (kg/kg)  mass fraction in air of water vapor
    REAL(wp), INTENT(in)         :: qc(:,:)                             !< (kg/kg)  mass fraction in air of cloud liquid
    REAL(wp), INTENT(in)         :: qi(:,:)                             !< (kg/kg)  mass fraction in air of cloud ice
    REAL(wp), INTENT(in)         :: qr(:,:)                             !< (kg/kg)  mass fraction in air of rain
    REAL(wp), INTENT(in)         :: qs(:,:)                             !< (kg/kg)  mass fraction in air of snow
    REAL(wp), INTENT(in)         :: qg(:,:)                             !< (kg/kg)  mass fraction in air of graupel

    REAL(wp)                     :: ein(SIZE(rho,1),SIZE(rho,2))        !< (J/m3)  internal  energy density

    INTEGER                      :: jc, jk                              !< loop indices

    !$ACC DATA PRESENT(rho, ta, qv, qc, qi, qr, qs, qg, ein)

    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, jks, jke) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jk = jks, jke
       DO jc = jcs, jce

          ein(jc,jk) =                                                & !< (J/m3)  internal energy density
               &       rho(jc,jk)                                     & !< (kg/m3) mass density of air
               &      *internal_energy(                               & !< (J/kg)  specific internal energy
               &                       ta(jc,jk),                     & !< (K)     temperature
               &                       qv(jc,jk),                     & !< (kg/kg) specific humidity
               &                       qc(jc,jk)+qr(jc,jk),           & !< (kg/kg) mass fraction of liquid phases in air
               &                       qi(jc,jk)+qs(jc,jk)+qg(jc,jk), & !< (kg/kg) mass fraction of frozen phases in air
               &                       1.0_wp   ,                     & !< ()      unused linear scaling factor
               &                       1.0_wp   )                       !< ()      unused linear scaling factor

       END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END FUNCTION atm_energy_diag_ein

  !-------------------------------------------------------------------

  PURE FUNCTION atm_energy_diag_ekh(jcs, jce, jks, jke, &
       &                            rho,                &
       &                            ua,                 &
       &                            va)                 &
       &                     RESULT(ekh)

    INTEGER , INTENT(in)         :: jcs, jce, jks, jke                  !< domain indices

    REAL(wp), INTENT(in)         :: rho(:,:)                            !< (kg/m3) mass density
    REAL(wp), INTENT(in)         :: ua(:,:)                             !< (m/s)   zonal velocity
    REAL(wp), INTENT(in)         :: va(:,:)                             !< (m/s)   meridional velocity

    REAL(wp)                     :: ekh(SIZE(rho,1),SIZE(rho,2))        !< (J/m3)  horizontal kinetic energy density

    INTEGER                      :: jc, jk                              !< loop indices

    !$ACC DATA PRESENT(rho, ua, va, ekh)

    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, jks, jke) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jk = jks, jke
       DO jc = jcs, jce

          !> Here horizontal kinetic energy is computed on full levels in the
          !! circumcenters of the cells, where the mass density is defined.
          !! The horizontal velocity components in the mass points are diagnosed
          !! from the normal wind components on the same full level in the mid
          !! points of the edges of a cell.
          !! Further it is assumed that the components (u,v,w) are locally orthogonal
          !! to each other.
          ! 
          ekh(jc,jk) =                                                & !< (J/m3)  horizontal kinetic energy density
               &        rho(jc,jk)                                    & !< (kg/m3) mass density of air
               &       *0.5_wp*(ua(jc,jk)**2 + va(jc,jk)**2)            !< (J/kg)  specific horizontal kinetic energy

       END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END FUNCTION atm_energy_diag_ekh

  !-------------------------------------------------------------------

  PURE FUNCTION atm_energy_diag_ekv(jcs, jce, jks, jke, &
       &                            rho,                &
       &                            wa)                 &
       &                     RESULT(ekv)

    INTEGER , INTENT(in)         :: jcs, jce, jks, jke                  !< domain indices

    REAL(wp), INTENT(in)         :: rho(:,:)                            !< (kg/m3) mass density
    REAL(wp), INTENT(in)         :: wa(:,:)                             !< (m/s)   vertical velocity, at half level

    REAL(wp)                     :: ekv(SIZE(rho,1),SIZE(rho,2))        !< (J/m3)  vertical kinetic energy density

    INTEGER                      :: jc, jk                              !< loop indices

    !$ACC DATA PRESENT(rho, wa, ekv)

    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, jks, jke) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jk = jks, jke
       DO jc = jcs, jce

          !> Here vertical kinetic energy is computed on full levels in the
          !! circumcenters of the cells, where the mass density is defined.
          !! The vertical velocity in the mass point is interpolated from the
          !! circumcenters on the upper and lower interfaces of a cell to the
          !! center on the mid level of a cell.
          !! Further it is assumed that the components (u,v,w) are locally orthogonal
          !! to each other.
          ! 
          ekv(jc,jk) =                                                & !< (J/m3)  vertical kinetic energy density
               &        rho(jc,jk)                                    & !< (kg/m3) mass density of air
               &       *0.5_wp*((0.5_wp*(wa(jc,jk) + wa(jc,jk+1)))**2)  !< (J/kg)  specific vertical kinetic energy

       END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END FUNCTION atm_energy_diag_ekv

  !-------------------------------------------------------------------

  PURE FUNCTION atm_energy_diag_egp(jcs, jce, jks, jke, &
       &                            rho,                &
       &                            geopot)             &
       &                     RESULT(egp)

    INTEGER , INTENT(in)         :: jcs, jce, jks, jke                  !< domain indices

    REAL(wp), INTENT(in)         :: rho(:,:)                            !< (kg/m3) mass density
    REAL(wp), INTENT(in)         :: geopot(:,:)                         !< (J/kg)  geopotential

    REAL(wp)                     :: egp(SIZE(rho,1),SIZE(rho,2))        !< (J/m3)  geopotential energy density

    INTEGER                      :: jc, jk                              !< loop indices

    !$ACC DATA PRESENT(rho, geopot, egp)

    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, jks, jke) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jk = jks, jke
       DO jc = jcs, jce

          egp(jc,jk) =                                                & !< (J/m3)  geopotential energy density
               &       rho(jc,jk)                                     & !< (kg/m3) mass density
               &      *geopot(jc,jk)                                    !< (J/kg)  spec. geopot. energy = geopotential

       END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END FUNCTION atm_energy_diag_egp

  !-------------------------------------------------------------------

  PURE FUNCTION atm_energy_diag_eto(jcs, jce, jks, jke, &
       &                            ein,                &
       &                            ekh,                &
       &                            ekv,                &
       &                            egp)                &
       &                     RESULT(eto)

    INTEGER , INTENT(in)         :: jcs, jce, jks, jke                  !< domain indices

    REAL(wp), INTENT(in)         :: ein(:,:)                            !< (J/m3)  internal     energy density
    REAL(wp), INTENT(in)         :: ekh(:,:)                            !< (J/m3)  hor. kinetic energy density
    REAL(wp), INTENT(in)         :: ekv(:,:)                            !< (J/m3)  vert.kinetic energy density
    REAL(wp), INTENT(in)         :: egp(:,:)                            !< (J/m3)  geopotential energy density

    REAL(wp)                     :: eto(SIZE(ein,1),SIZE(ein,2))        !< (J/m3)  total        energy density

    INTEGER                      :: jc, jk                              !< loop indices

    !$ACC DATA PRESENT(ein, ekh, ekv, egp, eto)

    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, jks, jke) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jk = jks, jke
       DO jc = jcs, jce

          eto(jc,jk) = ein(jc,jk) + ekh(jc,jk) + ekv(jc,jk) + egp(jc,jk)!< (J/m3)  total energy density

       END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END FUNCTION atm_energy_diag_eto

  !-------------------------------------------------------------------

  PURE FUNCTION atm_energy_vint(jcs, jce, jks, jke, &
       &                        dz,                 &
       &                        x)                  &
       &                 RESULT(xvi)

    INTEGER , INTENT(in)     :: jcs, jce, jks, jke                      !< domain indices

    REAL(wp), INTENT(in)     :: dz(:,:)                                 !< (m)    layer thickness
    REAL(wp), INTENT(in)     :: x(:,:)                                  !< (<x>)  vertically resolved field

    REAL(wp)                 :: xvi(SIZE(x,1))                          !< (<x>m) vertically integrated field

    INTEGER                  :: jc, jk                                  !< loop indices

    !$ACC DATA PRESENT(dz, x, xvi)

    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, jks, jke) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = jcs, jce
       xvi(jc) = 0.0_wp
    END DO

    !$ACC LOOP SEQ
    DO jk = jks, jke
       !$ACC LOOP GANG(STATIC: 1) VECTOR
       DO jc = jcs, jce

          xvi(jc) = xvi(jc) + x(jc,jk)*dz(jc,jk)

       END DO !jc
    END DO !jk
    !$ACC END PARALLEL

    !$ACC END DATA

  END FUNCTION atm_energy_vint

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_hint_1(jg)

    INTEGER, INTENT(in)     :: jg                                       !< grid index

    IF (ltimer) CALL timer_start(timer_atm_energy_hint)
    !
    !> Here the atmosphere energies are horizontally integrated
    !  for stage 1: ex1 -> ex1hi and ex1vi -> ex1ti
    !
    !> vertically resolved (vr) energy fields
    IF (ASSOCIATED(e(jg)%ein1hi)) CALL atm_energy_hint_vr(jg, e(jg)%ein1, e(jg)%ein1hi)
    IF (ASSOCIATED(e(jg)%ekh1hi)) CALL atm_energy_hint_vr(jg, e(jg)%ekh1, e(jg)%ekh1hi)
    IF (ASSOCIATED(e(jg)%ekv1hi)) CALL atm_energy_hint_vr(jg, e(jg)%ekv1, e(jg)%ekv1hi)
    IF (ASSOCIATED(e(jg)%egp1hi)) CALL atm_energy_hint_vr(jg, e(jg)%egp1, e(jg)%egp1hi)
    IF (ASSOCIATED(e(jg)%eto1hi)) CALL atm_energy_hint_vr(jg, e(jg)%eto1, e(jg)%eto1hi)
    !
    !> vertically integrated (vi) energy fields
    IF (ASSOCIATED(e(jg)%ein1ti)) CALL atm_energy_hint_vi(jg, e(jg)%ein1vi, e(jg)%ein1ti)
    IF (ASSOCIATED(e(jg)%ekh1ti)) CALL atm_energy_hint_vi(jg, e(jg)%ekh1vi, e(jg)%ekh1ti)
    IF (ASSOCIATED(e(jg)%ekv1ti)) CALL atm_energy_hint_vi(jg, e(jg)%ekv1vi, e(jg)%ekv1ti)
    IF (ASSOCIATED(e(jg)%egp1ti)) CALL atm_energy_hint_vi(jg, e(jg)%egp1vi, e(jg)%egp1ti)
    IF (ASSOCIATED(e(jg)%eto1ti)) CALL atm_energy_hint_vi(jg, e(jg)%eto1vi, e(jg)%eto1ti)
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_hint)

  END SUBROUTINE atm_energy_hint_1

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_hint_2(jg)

    INTEGER, INTENT(in)     :: jg                                       !< grid index

    IF (ltimer) CALL timer_start(timer_atm_energy_hint)
    !
    !> Here the atmosphere energies are horizontally integrated
    !  for stage 2: ex2 -> ex2hi and ex2vi -> ex2ti
    !
    !> vertically resolved (vr) energy fields
    IF (ASSOCIATED(e(jg)%ein2hi)) CALL atm_energy_hint_vr(jg, e(jg)%ein2, e(jg)%ein2hi)
    IF (ASSOCIATED(e(jg)%ekh2hi)) CALL atm_energy_hint_vr(jg, e(jg)%ekh2, e(jg)%ekh2hi)
    IF (ASSOCIATED(e(jg)%ekv2hi)) CALL atm_energy_hint_vr(jg, e(jg)%ekv2, e(jg)%ekv2hi)
    IF (ASSOCIATED(e(jg)%egp2hi)) CALL atm_energy_hint_vr(jg, e(jg)%egp2, e(jg)%egp2hi)
    IF (ASSOCIATED(e(jg)%eto2hi)) CALL atm_energy_hint_vr(jg, e(jg)%eto2, e(jg)%eto2hi)
    !
    !> vertically integrated (vi) energy fields
    IF (ASSOCIATED(e(jg)%ein2ti)) CALL atm_energy_hint_vi(jg, e(jg)%ein2vi, e(jg)%ein2ti)
    IF (ASSOCIATED(e(jg)%ekh2ti)) CALL atm_energy_hint_vi(jg, e(jg)%ekh2vi, e(jg)%ekh2ti)
    IF (ASSOCIATED(e(jg)%ekv2ti)) CALL atm_energy_hint_vi(jg, e(jg)%ekv2vi, e(jg)%ekv2ti)
    IF (ASSOCIATED(e(jg)%egp2ti)) CALL atm_energy_hint_vi(jg, e(jg)%egp2vi, e(jg)%egp2ti)
    IF (ASSOCIATED(e(jg)%eto2ti)) CALL atm_energy_hint_vi(jg, e(jg)%eto2vi, e(jg)%eto2ti)
    !
    IF (ltimer) CALL timer_stop(timer_atm_energy_hint)

  END SUBROUTINE atm_energy_hint_2

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_hint_vr(jg, x,  xgi)

    INTEGER, INTENT(in)      :: jg                                      !< grid index

    REAL(wp), INTENT(in)     :: x(:,:,:)                                !< globally resolved field

    REAL(wp), INTENT(out)    :: xgi(:)                                  !< globally integrated field

!$OMP PARALLEL
    CALL init(xgi)
!$OMP END PARALLEL

    CALL horizontal_sum(x, p(jg)%cells%area, p(jg)%cells%owned, xgi, lopenacc=.TRUE.)

  END SUBROUTINE atm_energy_hint_vr

  !-------------------------------------------------------------------

  SUBROUTINE atm_energy_hint_vi(jg, x, xgi)

    INTEGER, INTENT(in)      :: jg                                      !< grid index

    REAL(wp), INTENT(in)     :: x(:,:)                                  !< globally resolved field

    REAL(wp), INTENT(out)    :: xgi(:)                                  !< globally integrated field

    ! -----------------------------------------------------------------
    ! WORKAROUND START
    !
    ! Using 'horizontal_sum' of mo_statistics with a 2d input field on GPUs
    ! currently does not work: The resulting horizontal weighted sum is zero.
    ! 'horizontal_sum' is the public interface to private subroutines, which
    ! differ in the dimensionality of the arguments:
    ! - 2d input arrays and  2d weights :      'HorizontalSum_2D_InRange_2Dweights'
    ! - 3d input arrays with 2d weights : 'LevelHorizontalSum_3D_InRange_2Dweights'
    !
    ! As the subroutine for 3d input fields works on GPUs, this workaround simply
    ! uses the 3d variant to process the 2d field as a 3d field with 1 level.
    !
    ! => create a 3d field with 1 level
    !    copy the 2d field to the provisional 3d field
    !    compute the weighted sum using the 3d variant of 'horizontal_sum'

    REAL(wp) :: x3d(SIZE(x,1), 1, SIZE(x,2))                            !< 3d variable for x

    INTEGER  :: jb, jc                                                  !< loop indices

    !$ACC DATA CREATE(x3d) PRESENT(x)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jb = 1, SIZE(x,2)
       DO jc = 1, SIZE(x,1)

          x3d(jc,1,jb) = x(jc,jb)                                       !< copy over

       END DO
    END DO
    !$ACC END PARALLEL

    CALL horizontal_sum(x3d, p(jg)%cells%area, p(jg)%cells%owned, xgi, lopenacc=.TRUE.)

    !$ACC END DATA

    ! WORKAROUND END
    ! -----------------------------------------------------------------
    ! ORIGINAL START

!!$    CALL horizontal_sum(x, p(jg)%cells%area, p(jg)%cells%owned, xgi(1), lopenacc=.TRUE.)

    ! ORIGINAL END
    ! -----------------------------------------------------------------

  END SUBROUTINE atm_energy_hint_vi

  !-------------------------------------------------------------------

END MODULE mo_atm_energy_diag
