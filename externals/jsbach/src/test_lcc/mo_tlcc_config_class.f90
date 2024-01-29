!> tlcc config
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
!>#### Could contain namelist info for the test lcc proc
!>
!> currently of no use beyond being required by the infrastructure
!>
MODULE mo_tlcc_config_class
#ifndef __NO_JSBACH__
  
  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_tlcc_config

  TYPE, EXTENDS(t_jsb_config) :: t_tlcc_config

   CONTAINS
     PROCEDURE :: Init => Init_tlcc_config
  END type t_tlcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tlcc_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize tlcc process
  !
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize tlcc process
  !!
  !! @param[inout]     config     Configuration type of process (t_tlcc_config)
  !!  
  SUBROUTINE Init_tlcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_tlcc_config), INTENT(inout) :: config !<Configuration type of process (t_tlcc_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    ! -------------------------------------------------------------------------------------------------- !

    NAMELIST /jsb_tlcc_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_tlcc_config'

    CALL message(TRIM(routine), 'Starting tlcc configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_tlcc.nc'
    ic_filename         = 'ic_land_tlcc.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_tlcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_tlcc_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    
  END SUBROUTINE Init_tlcc_config

#endif
END MODULE mo_tlcc_config_class
