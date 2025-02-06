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

MODULE mo_rte_rrtmgp_setup

  USE mo_kind, ONLY: wp
  USE mo_gas_optics_rrtmgp, ONLY: ty_gas_optics_rrtmgp
  USE mo_gas_concentrations, ONLY: ty_gas_concs
  USE mo_cloud_optics,      ONLY: ty_cloud_optics
  USE mo_load_coefficients, ONLY: load_and_init
  USE mo_load_cloud_coefficients, ONLY: load_cld_lutcoeff, load_cld_padecoeff
  USE mo_cloud_gas_profiles, ONLY: init_gas_profiles
  USE mo_rte_kind,   ONLY: wl
  USE mo_rte_config, ONLY: rte_config_checks
  USE mo_snow_ice_reff, ONLY: init_reff_funeedles

  IMPLICIT NONE

  TYPE(ty_gas_optics_rrtmgp) :: k_dist_lw, k_dist_sw
  TYPE(ty_cloud_optics)      :: cloud_optics_lw, cloud_optics_sw
  REAL(wp)                   :: inhoml, inhomi, inhoms

  PUBLIC :: rte_rrtmgp_basic_setup
  PUBLIC :: k_dist_lw, k_dist_sw
  PUBLIC :: cloud_optics_lw, cloud_optics_sw
  PUBLIC :: stop_on_err
  PUBLIC :: inhoml, inhomi, inhoms

CONTAINS

  SUBROUTINE stop_on_err(msg)
    USE iso_fortran_env, ONLY : error_unit
    CHARACTER(len=*), INTENT(IN) :: msg

    IF(msg /= "") THEN
      WRITE(error_unit, *) msg
      STOP
    END IF
  END SUBROUTINE

  SUBROUTINE rte_rrtmgp_basic_setup(nproma, nlev, &
    pressure_scale_, droplet_scale_, &
    zinhoml_, zinhomi_, zinhoms_)

!!    USE mo_radiation_cloud_optics, ONLY: setup_cloud_optics

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nproma, nlev
    REAL(wp), INTENT(IN) :: pressure_scale_, droplet_scale_, &
      zinhoml_, zinhomi_, zinhoms_
    INTEGER :: i
    REAL(wp) :: log_factor
    INTEGER, PARAMETER :: n_gas_names = 8
    LOGICAL, PARAMETER :: luse_luts = .TRUE.
    CHARACTER(len=5), PARAMETER :: gas_names(n_gas_names) = (/ &
       'h2o  ', 'co2  ', 'ch4  ', 'o2   ', 'o3   ', 'n2o  ','cfc11', 'cfc12'/)
    TYPE(ty_gas_concs) :: gas_concs

    !pressure_scale = pressure_scale_
    !droplet_scale = droplet_scale_

!!    CALL setup_cloud_optics(droplet_scale_, &
!!      zinhoml1_, zinhoml2_, zinhoml3_, zinhomi_)
    inhoml = zinhoml_
    inhomi = zinhomi_
    inhoms = zinhoms_

    ! Turn off error checking
    CALL rte_config_checks(LOGICAL(.FALSE.,wl))  ! Dressing to please gfortran

    !
    ! Initialize RRTMGP with the names of the gases that will be available
    !
    CALL stop_on_err(gas_concs%init(gas_names))
    DO i = 1, n_gas_names
      CALL stop_on_err(gas_concs%set_vmr(gas_names(i), 0._wp))
    ENDDO
    CALL load_and_init(k_dist_lw, 'coefficients_lw.nc', gas_concs)
    CALL load_and_init(k_dist_sw, 'coefficients_sw.nc', gas_concs)
!!    CALL load_and_init
!!    CALL load_and_init(cloud_optics_sw, 'rrtmgp-cloud-optics-coeffs-sw.nc')
    IF (luse_luts) THEN
      CALL load_cld_lutcoeff(cloud_optics_lw, 'rrtmgp-cloud-optics-coeffs-lw.nc')
      CALL load_cld_lutcoeff(cloud_optics_sw, 'rrtmgp-cloud-optics-coeffs-sw.nc')
    ELSE
      CALL load_cld_padecoeff(cloud_optics_lw, 'rrtmgp-cloud-optics-coeffs-lw.nc')
      CALL load_cld_padecoeff(cloud_optics_sw, 'rrtmgp-cloud-optics-coeffs-sw.nc')
    ENDIF 
    CALL stop_on_err(cloud_optics_lw%set_ice_roughness(2))
    CALL stop_on_err(cloud_optics_sw%set_ice_roughness(2))

    CALL init_gas_profiles
    CALL init_reff_funeedles
  END SUBROUTINE rte_rrtmgp_basic_setup


END MODULE mo_rte_rrtmgp_setup
