!> Contains structures and methods for seb config
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
MODULE mo_seb_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_config

  TYPE, EXTENDS(t_jsb_config) :: t_seb_config
     CHARACTER(len=10)           :: scheme_albedo !< Scheme to use for computation of albedo (echam5 only for now)
     REAL(wp)                    :: lake_mixed_layer_depth !< Depth of lake mixed layer                             [m]
     REAL(wp)                    :: lake_min_ice_depth     !< Minimum ice thickness to start ice formation on lakes [m]
     LOGICAL                     :: l_ice_on_lakes         !< Whether to use ice on lakes
     INTEGER                     :: niter_tmx              !< For tmx: number of iterations in seb solver
   CONTAINS
     PROCEDURE :: Init => Init_seb_config
  END type t_seb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_config_class'

CONTAINS

  SUBROUTINE Init_seb_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_seb_config), INTENT(inout) :: config

    LOGICAL  :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(len=10)           :: scheme_albedo
    REAL(wp)                    :: lake_mixed_layer_depth, lake_min_ice_depth
    LOGICAL                     :: l_ice_on_lakes
    INTEGER                     :: niter_tmx

    NAMELIST /jsb_seb_nml/         &
      active,                      &
      ic_filename, bc_filename,    &
      scheme_albedo,               &
      lake_mixed_layer_depth,      &
      lake_min_ice_depth,          &
      l_ice_on_lakes,              &
      niter_tmx

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_seb_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting seb configuration')

    ! Set defaults
    active          = .TRUE.
    ic_filename     = 'ic_land_srf.nc'
    bc_filename     = 'bc_land_srf.nc'
    scheme_albedo   = 'echam5'
    lake_mixed_layer_depth = 10._wp
    lake_min_ice_depth     = 0.05_wp
    l_ice_on_lakes  = .TRUE.
    niter_tmx       = 0

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_seb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_seb_nml)

    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%scheme_albedo             = scheme_albedo
    config%lake_mixed_layer_depth    = lake_mixed_layer_depth
    config%lake_min_ice_depth        = lake_min_ice_depth
    config%l_ice_on_lakes            = l_ice_on_lakes
    config%niter_tmx                 = niter_tmx

    CALL close_nml(nml_handler)

  END SUBROUTINE Init_seb_config

#endif
END MODULE mo_seb_config_class
