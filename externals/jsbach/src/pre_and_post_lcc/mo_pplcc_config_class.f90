!> pplcc (pre- and post lcc) config
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
!>#### Could contain namelist info for pplcc, but currently of no use beyond infrastructural requirements
!>
!> currently of no use beyond being required by the infrastructure
!>
MODULE mo_pplcc_config_class
#ifndef __NO_JSBACH__
  
  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_pplcc_config

  TYPE, EXTENDS(t_jsb_config) :: t_pplcc_config

   CONTAINS
     PROCEDURE :: Init => Init_pplcc_config
  END type t_pplcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pplcc_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize pplcc process
  !
  SUBROUTINE Init_pplcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_pplcc_config), INTENT(inout) :: config !<Configuration type of process (t_pplcc_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    ! -------------------------------------------------------------------------------------------------- !
    
    NAMELIST /jsb_pplcc_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pplcc_config'

    CALL message(TRIM(routine), 'Starting pplcc configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_pplcc.nc'
    ic_filename         = 'ic_land_pplcc.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_pplcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_pplcc_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    
  END SUBROUTINE Init_pplcc_config

#endif
END MODULE mo_pplcc_config_class
