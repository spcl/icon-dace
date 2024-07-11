!> QUINCY phenology process config
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
!>#### define phenology config structure, read phenology namelist and init configuration parameters
!>
MODULE mo_q_pheno_config_class
#ifndef __NO_QUINCY__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_pheno_config

  ! ======================================================================================================= !
  !> Q_PHENO_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_q_pheno_config
    REAL(wp)            :: lai_max              !< site specific lai, used to init lai_max in Q_PHENO_
    CONTAINS
      PROCEDURE :: Init => Init_q_pheno_config
  END type t_q_pheno_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_pheno_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read Q_PHENO_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_q_pheno_config(config)
#ifdef __QUINCY_STANDALONE__
    USE mo_namelist,              ONLY: open_nml, POSITIONED, position_nml, close_nml
#else
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_pheno_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    REAL(wp)                    :: lai_site_specific      !< site specific lai, used to init lai_max in Q_PHENO_
#ifdef __QUINCY_STANDALONE__
    NAMELIST /phenology_ctl/  &
#else
    NAMELIST /q_pheno_nml/  &
#endif
      & active,               &
      & ic_filename,          &
      & bc_filename,          &
      & lai_site_specific
    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_pheno_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting phenology configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_phenology.nc'
    ic_filename         = 'ic_land_phenology.nc'
    lai_site_specific   = 6.0_wp
    !first_call_phenology_this_day = .FALSE.
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
#ifdef __QUINCY_STANDALONE__
    nml_unit = position_nml('phenology_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, phenology_ctl)
#else
    nml_unit = position_nml('q_pheno_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, q_pheno_nml)
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                     = active
    config%ic_filename                = ic_filename
    config%bc_filename                = bc_filename
    config%lai_max                    = lai_site_specific
    CALL close_nml(nml_handler)
  END SUBROUTINE Init_q_pheno_config

#endif
END MODULE mo_q_pheno_config_class
