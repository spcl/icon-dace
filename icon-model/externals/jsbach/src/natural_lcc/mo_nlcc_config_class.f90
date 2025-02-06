!> Contains structures and methods for nlcc config
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
MODULE mo_nlcc_config_class
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
  PUBLIC :: t_nlcc_config

  TYPE, EXTENDS(t_jsb_config) :: t_nlcc_config
    LOGICAL             :: &
      & init_running_means, &
      & init_nlcc
    REAL(wp)            :: &
      & accelerate_nlcc

   CONTAINS
     PROCEDURE :: Init => Init_nlcc_config
  END type t_nlcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nlcc_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize nlcc process
  !!
  !! @param[inout]     config     Configuration type of process (t_nlcc_config)
  !!
  SUBROUTINE Init_nlcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    CLASS(t_nlcc_config), INTENT(inout) :: config
    TYPE(t_jsb_vgrid), POINTER  :: pft_layer
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: init_running_means
    LOGICAL                     :: init_nlcc
    REAL(wp)                    :: accelerate_nlcc

    NAMELIST /jsb_nlcc_nml/  &
      &  active,             &
      &  ic_filename,        &
      &  bc_filename,        &
      &  init_running_means, &
      &  init_nlcc,          &
      &  accelerate_nlcc

    INTEGER :: nml_handler, nml_unit, istat, npft, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_nlcc_config'

    CALL message(TRIM(routine), 'Starting nlcc configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_nlcc.nc'
    ic_filename         = 'ic_land_nlcc.nc'
    init_running_means  = .FALSE.
    init_nlcc           = .FALSE.
    accelerate_nlcc     = 1._wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_nlcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_nlcc_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%init_running_means  = init_running_means
    IF (init_running_means) CALL message(TRIM(routine), 'Initializing climate buffer')
    config%init_nlcc         = init_nlcc
    IF (init_nlcc) CALL message(TRIM(routine), 'Initializing NLCC')
    config%accelerate_nlcc   = accelerate_nlcc
    IF (accelerate_nlcc /= 1._wp) THEN
      WRITE(message_text, *) 'NLCC calculations accelerated by factor ', accelerate_nlcc
      CALL message(TRIM(routine), message_text)
    END IF

    ! Create vertical axis
    npft=11 ! TODO here we can not grab it from tile%no_of_children or mem%no_of_children
    pft_layer  => new_vgrid('pfts', ZAXIS_GENERIC, npft, &
      & levels=(/ (REAL(i,kind=wp),i=1,npft) /), units='')
    CALL Register_vgrid(pft_layer)

  END SUBROUTINE Init_nlcc_config

#endif
END MODULE mo_nlcc_config_class
