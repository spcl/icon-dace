!> calculate vegetation respiration (QUINCY)
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### calculate vegetation maintenance-respiration considering, e.g., temperature and acclimation
!>
MODULE mo_q_veg_respiration
#ifndef __NO_QUINCY__

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_maintenance_respiration, temperature_response_respiration

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_respiration'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Calculate temperature response of maintenance respiration
  !! following the temperature function of Lloyd and Taylor, 1996, Functional Ecology, but
  !! assuming that the base rate of respiration acclimates to growth temperature according to
  !! Atkin et al, XXXX New Phytologist
  !!
  !! Input: temperature, growth temperature
  !!
  !! Output: temperature rate modifier
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION temperature_response_respiration(temp, temp_mavg) RESULT(temp_response)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_veg_constants,           ONLY: t_maint_k1, t_maint_k2, t_maint_k3, f_resp_acclim, t_acclim_zero

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: temp                      !< instantaneous temperature (K)
    REAL(wp), INTENT(in) :: temp_mavg                 !< average temperature for acclimation (K)
    REAL(wp)             :: temp_response             !< temperature response of respiration (1 at 10 degC)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)             :: hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':temperature_response_respiration'


    hlp1 = t_maint_k3 + 1.0_wp / (1._wp / t_maint_k2 - log(eps8) / t_maint_k1)

    IF(temp >= hlp1) THEN
      ! IF flag_t_resp_acclimation is false: temp_mavg == t_acclim_zero
      temp_response = (10._wp**(f_resp_acclim * (temp_mavg - t_acclim_zero))) * &         ! this term is 1 if flag_t_resp_acclimation is false
                       exp(t_maint_k1 * (1._wp / t_maint_k2 - 1._wp / (temp - t_maint_k3)))
    ELSE
      temp_response = 0.0_wp
    ENDIF

  END FUNCTION temperature_response_respiration


  !-----------------------------------------------------------------------------------------------------
  !> Calculate maintenance respiration
  !! following the base rate as defined by Sprugel, 1996, modified by Temperature according to
  !!   Lloyd and Taylor, 1996.
  !!
  !! Input: tissue N and temperature (leaf respiration calculated in photosynthesis and NOT modified
  !!        here)
  !!
  !! Output: maintenance respiration
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_maintenance_respiration( leaf_carbon           , &
                                                        t_wood                , &
                                                        t_wood_mavg           , &
                                                        t_root                , &
                                                        t_root_mavg           , &
                                                        sap_wood_nitrogen     , &
                                                        fine_root_nitrogen    , &
                                                        coarse_root_nitrogen)   &
                                                        RESULT(maint_respiration)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_veg_constants,           ONLY: fmaint_rate_base, fmaint_rate_wood, k_maint_resp_suppress

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: leaf_carbon              !< carbon content of leafs (veg_pool)
    REAL(wp), INTENT(in) :: t_wood                   !< temperature of the stem (K)
    REAL(wp), INTENT(in) :: t_wood_mavg              !< average temperature for acclimation (K)
    REAL(wp), INTENT(in) :: t_root                   !< temperature of the roots (K)
    REAL(wp), INTENT(in) :: t_root_mavg              !< average temperature for acclimation (K)
    REAL(wp), INTENT(in) :: sap_wood_nitrogen        !< nitrogen content of sapwood (mol N / m2)
    REAL(wp), INTENT(in) :: coarse_root_nitrogen     !< nitrogen content of sapwood (mol N / m2)
    REAL(wp), INTENT(in) :: fine_root_nitrogen       !< nitrogen content of sapwood (mol N / m2)
                                                     ! (calculated in photosynthesis; micro-mol C / m2 / s)
    REAL(wp)             :: maint_respiration        !< maintenance respiration (micro-mol C / m2 / s)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_maintenance_respiration'

    !> 1.0 calculate maintenance respiration
    maint_respiration = temperature_response_respiration(t_wood, t_wood_mavg) * &
                        fmaint_rate_wood * sap_wood_nitrogen + &
                        temperature_response_respiration(t_root, t_root_mavg) * &
                        fmaint_rate_base * ( coarse_root_nitrogen + fine_root_nitrogen )

    !> 2.0 reduce maintenance respiration if there are no leafs
    !!
    IF (leaf_carbon < eps8) THEN
      maint_respiration = maint_respiration * k_maint_resp_suppress
    ENDIF

  END FUNCTION calc_maintenance_respiration

#endif
END MODULE mo_q_veg_respiration
