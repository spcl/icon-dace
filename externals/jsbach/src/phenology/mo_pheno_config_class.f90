!> Contains structures and methods for phenology config
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
MODULE mo_pheno_config_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pheno_config

  ! ======================================================================================================= !
  !> PHENO_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_pheno_config
    LOGICAL             :: l_forestRegrowth           !< If TRUE, the maxLAI is derived from biomass via allometric relationships
    CHARACTER(len=15)   :: scheme                     !< Phenology scheme to use ('logrop', 'climatology')
    !LOGICAL :: first_call_phenology_this_day         !< TRUE marks if a new day has just started and
    !                                                 !< update_phenology is called for the top tile
    CONTAINS
      PROCEDURE :: Init => Init_pheno_config
  END type t_pheno_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read PHENO_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_pheno_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_pheno_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    LOGICAL                     :: l_forestRegrowth
    CHARACTER(len=15)           :: scheme
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    NAMELIST /jsb_pheno_nml/  &
      & active,               &
      & l_forestRegrowth,     &
      & scheme,               &
      & ic_filename,          &
      & bc_filename
    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pheno_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting phenology configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active              = .TRUE.
    l_forestRegrowth    = .FALSE.
    scheme              = 'logrop'
    bc_filename         = 'bc_land_phenology.nc'
    ic_filename         = 'ic_land_phenology.nc'
    !first_call_phenology_this_day = .FALSE.
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('jsb_pheno_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_pheno_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                     = active
    config%l_forestRegrowth           = l_forestRegrowth
    config%scheme                     = TRIM(scheme)
    config%ic_filename                = ic_filename
    config%bc_filename                = bc_filename
    CALL close_nml(nml_handler)
  END SUBROUTINE Init_pheno_config

#endif
END MODULE mo_pheno_config_class
