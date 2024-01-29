! ===============================================================================================================================
! === THIS IS A TEMPLATE. PLEASE REPLACE ...                                                                                  ===
! === ... <FIRST NAME LAST NAME>    with your name e.g. "Martin Mustermann"                                                   ===
! === ... <DATE>                    with the present date in YYYY-MM-DD format e.g. "2016-04-30"                              ===
! === ... <PROCESS_NAME_LOWER_CASE> with your process name e.g. "seb" for surface energy balance                              ===
! ===                                                                                                                         ===
! === Then go through the code line by line to adapt it to your needs:                                                        ===
! === !X marks lines with examples (e.g:), templates (if necessary:) and implementation points (Implementation:) that you     ===
! === have to adapt to your process.                                                                                          ===
! ===                                                                                                                         ===
! === Finally delete this header.                                                                                             ===
! ===============================================================================================================================
!> Contains structures and methods for <PROCESS_NAME_LOWER_CASE> config
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
MODULE mo_<PROCESS_NAME_LOWER_CASE>_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_<PROCESS_NAME_LOWER_CASE>_config

  TYPE, EXTENDS(t_jsb_config) :: t_<PROCESS_NAME_LOWER_CASE>_config
     !X e.g: LOGICAL                     :: use_alb_veg_simple

   CONTAINS
     PROCEDURE :: Init => Init_<PROCESS_NAME_LOWER_CASE>_config
  END type t_<PROCESS_NAME_LOWER_CASE>_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_<PROCESS_NAME_LOWER_CASE>_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize <PROCESS_NAME_LOWER_CASE> process
  !!
  !! @param[inout]     config     Configuration type of process (t_<PROCESS_NAME_LOWER_CASE>_config)
  !!  
  SUBROUTINE Init_<PROCESS_NAME_LOWER_CASE>_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_<PROCESS_NAME_LOWER_CASE>_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    !X e.g: LOGICAL                     :: use_alb_veg_simple

    NAMELIST /jsb_<PROCESS_NAME_LOWER_CASE>_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename,                 &
         !X e.g: use_<PROCESS_NAME_LOWER_CASE>_veg_simple

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_<PROCESS_NAME_LOWER_CASE>_config'

    CALL message(TRIM(routine), 'Starting <PROCESS_NAME_LOWER_CASE> configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_<PROCESS_NAME_LOWER_CASE>.nc'
    ic_filename         = 'ic_land_<PROCESS_NAME_LOWER_CASE>.nc'
    !X e.g: use_alb_veg_simple        = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_<PROCESS_NAME_LOWER_CASE>_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_<PROCESS_NAME_LOWER_CASE>_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    !X e.g: config%use_alb_veg_simple        = use_alb_veg_simple
    
  END SUBROUTINE Init_<PROCESS_NAME_LOWER_CASE>_config

#endif
END MODULE mo_<PROCESS_NAME_LOWER_CASE>_config_class
