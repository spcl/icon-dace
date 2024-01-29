!> Contains some functions for physcial schemes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_phy_schemes
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_convect_tables_iface, ONLY: tlucua, jptlucu1, jptlucu2
  USE mo_jsb_thermo_iface, ONLY: specific_humidity, sat_pres_water, sat_pres_ice, potential_temperature
  USE mo_jsb_surface_exchange_iface, ONLY: sfc_exchange_coefficients

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: specific_humidity, qsat_water, qsat_ice, sat_pres_water, qsat, &
    & tlucua, jptlucu1, jptlucu2, &
    & q_effective, surface_dry_static_energy, heat_transfer_coef, &
    & thermal_radiation, lwnet_from_lwdown, potential_temperature, &
    & sfc_exchange_coefficients
    ! Note: nvhpc doesn't support function pointers
    ! & register_exchange_coefficients_procedure, exchange_coefficients, registered_exchange_coefficients_procedure

  ! ABSTRACT INTERFACE
  !   PURE SUBROUTINE i_exchange_coefficients_procedure(           &
  !     & dz,                                                 &
  !     & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
  !     & km, kh, km_neutral, kh_neutral                      &
  !     & )
  !     IMPORT :: wp

  !     REAL(wp), INTENT(in) :: &
  !       dz,        &
  !       thetam1,   &
  !       pqm1 ,     &
  !       mwind,     &
  !       rough_m,   &
  !       theta_sfc, &
  !       qsat_sfc
  !       !
  !     REAL(wp), INTENT(out) :: &
  !       km,         &
  !       kh,         &
  !       km_neutral, &
  !       kh_neutral

  !   END SUBROUTINE
  ! END INTERFACE

  ! PROCEDURE(i_exchange_coefficients_procedure), POINTER :: registered_exchange_coefficients_procedure => NULL()

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes'

CONTAINS

  ! SUBROUTINE register_exchange_coefficients_procedure(exchange_coefficients_procedure)

  !   PROCEDURE(i_exchange_coefficients_procedure) :: exchange_coefficients_procedure

  !   registered_exchange_coefficients_procedure => exchange_coefficients_procedure

  ! END SUBROUTINE register_exchange_coefficients_procedure

  ! SUBROUTINE exchange_coefficients(                       &
  !   & dz,                                                 &
  !   & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
  !   & km, kh, km_neutral, kh_neutral                      &
  !   & )

  !   REAL(wp), DIMENSION(:), INTENT(in) :: &
  !     dz,        &
  !     thetam1,   &
  !     pqm1 ,     &
  !     mwind,     &
  !     rough_m,   &
  !     theta_sfc, &
  !     qsat_sfc
  !     !
  !   REAL(wp), DIMENSION(:), INTENT(out) :: &
  !     km,         &
  !     kh,         &
  !     km_neutral, &
  !     kh_neutral

  !   INTEGER :: ic, nc

  !   nc = SIZE(dz)

  !   !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
  !   DO ic = 1, nc
  !     CALL registered_exchange_coefficients_procedure(                                &
  !       & dz(ic),                                                                     &
  !       & pqm1(ic), thetam1(ic), mwind(ic), rough_m(ic), theta_sfc(ic), qsat_sfc(ic), &
  !       & km(ic), kh(ic), km_neutral(ic), kh_neutral(ic)                              &
  !       & )
  !   END DO
  !   !$ACC END PARALLEL LOOP

  ! END SUBROUTINE exchange_coefficients

  REAL(wp) FUNCTION qsat_water(temperature, pressure, use_convect_tables)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]
    LOGICAL,  INTENT(in) :: use_convect_tables

    IF (use_convect_tables) THEN
      qsat_water = qsat(temperature, pressure)
    ELSE
      qsat_water = specific_humidity(sat_pres_water(temperature), pressure)
    END IF

  END FUNCTION qsat_water

  REAL(wp) FUNCTION qsat_ice(temperature, pressure, use_convect_tables)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]
    LOGICAL,  INTENT(in) :: use_convect_tables

    IF (use_convect_tables) THEN
      qsat_ice = qsat(temperature, pressure)
    ELSE
      qsat_ice = specific_humidity(sat_pres_ice(temperature), pressure)
    END IF

  END FUNCTION qsat_ice

  ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
  ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
  ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
  !
  ! Note: temporary for backward compatibility
  REAL(wp) FUNCTION qsat(temperature, pressure)

    !$ACC ROUTINE SEQ

    USE mo_jsb_convect_tables_iface, ONLY: tlucua, &        ! Table for water vapor pressure e_s multiplied by R_d/R_v
      &                                    jptlucu1, jptlucu2
    USE mo_jsb_physical_constants, ONLY: rvd1                    ! = R_v/R_d - 1

    REAL(wp), INTENT(in) :: temperature      ! Air temperature [K]
    REAL(wp), INTENT(in) :: pressure         ! Pressure [Pa]

    INTEGER :: it
    REAL(wp) :: tluc

    it = NINT(temperature*1000._wp)
!!$    it = MIN(MAX(NINT(temperature*1000._wp), 50000),400000)    ! Temporary fix
    tluc = tlucua(it) !<- TODO catch seg fault for unplausible temperatures? (e.g. see mo_jsb4_forcing)

    IF (it >= jptlucu1 .AND. it <= jptlucu2) THEN
      qsat = tluc / (pressure - rvd1*tluc)
    ELSE
      qsat = HUGE(1._wp)
    END IF

  END FUNCTION qsat

  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION q_effective
  !! 
  !! out: q_effective
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) FUNCTION q_effective(qsat, qair, fsat, fair)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(in) :: &
      & qsat,               & !< Surface saturation specific humidity
      & qair,               & !< Air specific humidity
      & fsat, fair            !< Weighing factors for qsat and qair accounting for only partially wet surface.

    q_effective = fsat * qsat + (1._wp -fair) * qair

  END FUNCTION q_effective


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION surface_dry_static_energy
  !! 
  !! out: surface_dry_static_energy
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) FUNCTION surface_dry_static_energy(t_srf, qsat_srf, cpd_or_cvd)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: &
      & cpvd1        !< cpv/cpd-1

    REAL(wp), INTENT(in) :: &
      & t_srf,                      & !< Surface temperature
      & qsat_srf,                   & !< Surface saturation specific humidity
      & cpd_or_cvd

    ! surface_dry_static_energy = t_srf * cpd * ( 1._wp + cpvd1 * qsat_srf)
    surface_dry_static_energy = t_srf * cpd_or_cvd

  END FUNCTION surface_dry_static_energy


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION thermal_radiation
  !! 
  !! out: thermal_radiation
  !-----------------------------------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(wp) FUNCTION thermal_radiation(t_srf)

    USE mo_jsb_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: t_srf

    thermal_radiation = stbo * zemiss_def * t_srf**4._wp

  END FUNCTION thermal_radiation


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION lwnet_from_lwdown
  !! 
  !! out: lwnet_from_lwdown
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) FUNCTION lwnet_from_lwdown(lwdown, t_srf)

    !$ACC ROUTINE SEQ

    USE mo_physical_constants, ONLY: &
      stbo,       &  !< Stefan-Boltzmann constant
      zemiss_def     !< Surface emissivity

    REAL(wp), INTENT(in) :: lwdown   ! downward longwave radiation
    REAL(wp), INTENT(in) :: t_srf    ! surface temperature

    lwnet_from_lwdown = zemiss_def * (lwdown - stbo * t_srf**4._wp)

  END FUNCTION lwnet_from_lwdown


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION heat_transfer_coef
  !! 
  !! out: heat_transfer_coef
  !-----------------------------------------------------------------------------------------------------
  REAL(wp) FUNCTION heat_transfer_coef(drag, steplen)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: grav, cvdifts

    REAL(wp), INTENT(in) :: &
      & drag,      &
      & steplen

    heat_transfer_coef = drag / (cvdifts * grav * steplen)

  END FUNCTION heat_transfer_coef

#endif
END MODULE mo_phy_schemes
