!> vegetation process config (QUINCY)
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
!>#### define vegetation config structure, read vegetation namelist and init configuration parameters
!>
MODULE mo_veg_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the vegetation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_veg_config
    ! INTEGER          :: pft_id
    ! REAL(wp)         :: cohort_harvest_interval    !< interval of harvesting plants [yr], used with veg_dynamics_scheme="cohort"
    ! CHARACTER(15)    :: bnf_scheme                 !< select scheme to simulate BNF (biological nitrogen fixation)
    !                                                !! "none fixed standard optimal resistance unlimited"   
    ! CHARACTER(15)    :: veg_dynamics_scheme        !< select scheme to calculate within-tile vegetation dynamics + start from bareground
    !                                                !! none: constant mortality
    !                                                !! "none population cohort"   
    ! CHARACTER(15)    :: biomass_alloc_scheme       !< select scheme for veg biomass allocation: fixed dynamic optimal
    ! CHARACTER(15)    :: leaf_stoichom_scheme       !< select scheme for C/N leaf stoichometry:  fixed dynamic optimal
   CONTAINS
     PROCEDURE :: Init => Init_veg_config
  END type t_veg_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_config_class'

CONTAINS
  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_veg_config
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_veg_config), INTENT(inout) :: config    !< config type for veg
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename,   &
                                   bc_filename
    ! INTEGER                     :: plant_functional_type_id
    ! REAL(wp)                    :: cohort_harvest_interval
    ! CHARACTER(15)               :: bnf_scheme       
    ! CHARACTER(15)               :: veg_dynamics_scheme
    ! CHARACTER(15)               :: biomass_alloc_scheme
    ! CHARACTER(15)               :: leaf_stoichom_scheme
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_config'

    NAMELIST /jsb_veg_nml/       & 
      active,                       &
      bc_filename,                  &
      ic_filename!,                  &
      ! plant_functional_type_id,     &
      ! cohort_harvest_interval,      &
      ! bnf_scheme,                   & 
      ! veg_dynamics_scheme,          &
      ! biomass_alloc_scheme,         &
      ! leaf_stoichom_scheme

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting veg configuration')

    ! Set defaults
    active                         = .TRUE.
    bc_filename                    = 'bc_land_vegetation.nc'
    ic_filename                    = 'ic_land_vegetation.nc'
    ! plant_functional_type_id       = 1
    ! cohort_harvest_interval        = 80.0_wp
    ! bnf_scheme                     = "standard"
    ! veg_dynamics_scheme            = "population"
    ! biomass_alloc_scheme           = "fixed"
    ! leaf_stoichom_scheme           = "fixed"

    ! read the namelist
    ! &jsb_veg_nml
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_veg_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_veg_nml)

    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active                        = active
    config%bc_filename                   = bc_filename
    config%ic_filename                   = ic_filename
    ! config%pft_id                        = plant_functional_type_id
    ! config%cohort_harvest_interval       = cohort_harvest_interval
    ! config%bnf_scheme                    = TRIM(bnf_scheme)
    ! config%veg_dynamics_scheme           = TRIM(veg_dynamics_scheme)
    ! config%biomass_alloc_scheme          = biomass_alloc_scheme
    ! config%leaf_stoichom_scheme          = leaf_stoichom_scheme

  END SUBROUTINE Init_veg_config

#endif
END MODULE mo_veg_config_class
