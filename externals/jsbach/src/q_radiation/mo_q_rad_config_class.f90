!> QUINCY radiation process config
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
!>#### define radiation config structure, read radiation namelist and init configuration parameters
!>
MODULE mo_q_rad_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_rad_config

  ! ======================================================================================================= !
  !> Q_RAD_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_q_rad_config
    ! ...
    ! ...
  CONTAINS
    PROCEDURE :: Init => Init_q_rad_config
  END type t_q_rad_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_rad_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read Q_RAD_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_q_rad_config(config)
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC
#ifdef __QUINCY_STANDALONE__
    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml
#else
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_rad_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    INTEGER                     :: nspec
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_spectral

#ifdef __QUINCY_STANDALONE__
    NAMELIST /radiation_ctl/ &
#else
    NAMELIST /q_rad_nml/ &
#endif
      & active,                      &
      & ic_filename,                 &
      & bc_filename

    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_rad_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting radiation configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                    = .TRUE.
    bc_filename               = 'bc_land_q_rad.nc'  !< input file not yet available
    ic_filename               = 'ic_land_q_rad.nc'  !< input file not yet available
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
#ifdef __QUINCY_STANDALONE__
    nml_unit = position_nml('radiation_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, radiation_ctl)
#else
    nml_unit = position_nml('q_rad_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, q_rad_nml)
#endif
    CALL close_nml(nml_handler)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename

    ! number of spectal bands (3rd dimension of 'vgrid_spectral')
    nspec = 2162

    !> Create vertical axis - spectral bands
    !>   spectral bands of the light for Q_RAD_ variables
    !>
#ifdef __QUINCY_STANDALONE__
    vgrid_spectral  => new_vgrid('spectral_bands', ZAXIS_GENERIC, nspec, &
                      longname='Number of spectral bands the light is divided into, from Q_RAD_ QUINCY') !, &
                      ! units=''  , &
                      ! levels=   , &
                      ! lbounds=  , &
                      ! ubounds=  , &
                      ! dz=       )
    CALL register_vgrid(vgrid_spectral)
#endif
  END SUBROUTINE Init_q_rad_config

#endif
END MODULE mo_q_rad_config_class
