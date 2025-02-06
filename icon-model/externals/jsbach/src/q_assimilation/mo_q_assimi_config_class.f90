!> QUINCY assimilation process config
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
!>#### define assimilation config structure, read assimilation namelist and init configuration parameters
!>
MODULE mo_q_assimi_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_math_constants,ONLY: eps4
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_q_assimi_config

  ! ======================================================================================================= !
  !> Q_ASSIMI_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_q_assimi_config
    INTEGER                 :: ncanopy            !< number of canopy layers
    LOGICAL                 :: flag_optimal_Nfraction     !< optimise leaf internal N allocation
    LOGICAL                 :: flag_t_resp_acclimation    !< whether or not to use the respiration temperature acclimation factor
    LOGICAL                 :: flag_t_jmax_acclimation    !< acclimation of optimum temperature for Jmax
    CHARACTER(15)           :: canopy_layer_scheme        !< select scheme to calculate canopy layer thicknesses [standard|fapar]
    CHARACTER(15)           :: canopy_conductance_scheme  !< select scheme to calculate canopy conductance [medlyn|ballberry]
  CONTAINS
    PROCEDURE :: Init => Init_q_assimi_config
  END TYPE t_q_assimi_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_assimi_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read Q_ASSIMI_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_q_assimi_config(config)
#ifdef __QUINCY_STANDALONE__
    USE mo_namelist,              ONLY: open_nml, POSITIONED, position_nml, close_nml
#else
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
#endif
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_io,                ONLY: ZAXIS_GENERIC
    USE mo_jsb_math_constants,    ONLY: eps4
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_assimi_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    INTEGER                     :: ncanopy
    LOGICAL                     :: flag_optimal_Nfraction
    LOGICAL                     :: flag_t_resp_acclimation
    LOGICAL                     :: flag_t_jmax_acclimation
    CHARACTER(15)               :: canopy_layer_scheme
    CHARACTER(15)               :: canopy_conductance_scheme
    REAL(wp), ALLOCATABLE       :: canopy_layer_thickness_profile(:)  !< use LAI per layer [m2 m-2], indirectly defining the thickness (i.e., vertical height) of each single canopy layer
    REAL(wp)                    :: dAPAR
    REAL(wp)                    :: k_canopy_layer                     !< extinction coefficient to calculate canopy layer depth
    INTEGER                     :: icanopy
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_canopy_q_assimi

#ifdef __QUINCY_STANDALONE__
    NAMELIST /assimilation_ctl/  &
#else
    NAMELIST /q_assimi_nml/  &
#endif
      & active,                   &
      & ic_filename,              &
      & bc_filename,              &
      & ncanopy,                  &
      & flag_optimal_Nfraction,   &
      & flag_t_resp_acclimation,  &
      & flag_t_jmax_acclimation,  &
      & canopy_layer_scheme,      &
      & canopy_conductance_scheme

    INTEGER                     :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_assimi_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting assimi configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                        = .TRUE.
    bc_filename                   = 'bc_land_assimi.nc'
    ic_filename                   = 'ic_land_assimi.nc'
    ncanopy                       = 10           ! (default for FAPAR based canopy layers would be 10)
    flag_optimal_Nfraction        = .FALSE.      ! .TRUE. .FALSE.
    flag_t_resp_acclimation       = .FALSE.      ! .TRUE. .FALSE.
    flag_t_jmax_acclimation       = .FALSE.      ! .TRUE. .FALSE.
    canopy_layer_scheme           = "fapar"      !
    canopy_conductance_scheme     = "medlyn"
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
#ifdef __QUINCY_STANDALONE__
    nml_unit    = position_nml('assimilation_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, assimilation_ctl)
#else
    nml_unit    = position_nml('q_assimi_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, q_assimi_nml)
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                       = active
    config%ic_filename                  = ic_filename
    config%bc_filename                  = bc_filename
    config%ncanopy                      = ncanopy
    config%flag_optimal_Nfraction       = flag_optimal_Nfraction
    config%flag_t_resp_acclimation      = flag_t_resp_acclimation
    config%flag_t_jmax_acclimation      = flag_t_jmax_acclimation
    config%canopy_layer_scheme          = TRIM(canopy_layer_scheme)
    config%canopy_conductance_scheme    = TRIM(canopy_conductance_scheme)
    CALL close_nml(nml_handler)

    ! extinction coefficient for canopy layer depth calculation
    k_canopy_layer = 0.7_wp

    !> canopy layer depth calculation
    !>
    !! here the thickness of each canopy layer is calculated (saved in the vgrid as dz(:)) \n
    !! this thickness is not directly given (by e.g. using the unit meter), but represented by the leaf mass in this layer \n
    !! this leaf mass is measured as LAI [m2 m-2], that is, the actual thickness (in m) of canopy layer 2 depends on how much
    !!   vertical space (hight) needs to be covered to get the LAI specific to canopy layer 2
    !!
    !! layer 1 is the top canopy layer (reversed order of vgrid axis) \n
    !! and the value of the upper bound of a layer is larger than the value of the lower bound of the same layer
    !!
    ALLOCATE(canopy_layer_thickness_profile(config%ncanopy))
    SELECT CASE(TRIM(config%canopy_layer_scheme))
    CASE ("standard")
      canopy_layer_thickness_profile(:) = 0.5_wp
    CASE ("fapar")
      ! alternative canopy layer depth calculation, approximately proportional to FAPAR
      ! to reduce impact of actual canopy layering on PS calculation
      dAPAR = 1.0_wp / REAL(config%ncanopy, KIND = wp)
      DO icanopy = 1,config%ncanopy
        canopy_layer_thickness_profile(icanopy) = -1.0_wp / k_canopy_layer * LOG(1._wp - (dAPAR * icanopy) + eps4)
      END DO
    CASE DEFAULT
      WRITE(message_text,'(2a)') 'invalid canopy_layer_scheme', config%canopy_layer_scheme
      CALL finish(routine, message_text)
    END SELECT

    !> check validity of canopy-conductance scheme
    !>
    SELECT CASE(TRIM(config%canopy_conductance_scheme))
    CASE ("medlyn")
      CALL message(TRIM(routine), 'Using the canopy_conductance_scheme: '//TRIM(config%canopy_conductance_scheme)// &
        &                         ' the local lctlib_g1 variable is set to lctlib%g1_medlyn')
    CASE ("ballberry")
      CALL message(TRIM(routine), 'Using the canopy_conductance_scheme: '//TRIM(config%canopy_conductance_scheme)// &
        &                         ' the local lctlib_g1 variable is set to lctlib%g1_bberry')
    CASE DEFAULT
      WRITE(message_text,'(2a)') 'invalid canopy_conductance_scheme: ', TRIM(config%canopy_conductance_scheme)
      CALL finish(routine, message_text)
    END SELECT

    !> Create vertical axis - canopy layers
    !>
    !! these are "infrastructure" values, constant across PFTs and sites, and they do not interfere with tree height \n
    !! vgrid_canopy_q_assimi% lbounds, ubounds, and levels is calculated by 'new_vgrid()'
    !!
    !! note: the axis of the q_canopy_layer vgrid is used as "reversed", i.e., layer 1 is the top canopy layer \n
    !! levels is defined as the depth at the center of the layer: levels(:) = 0.5_wp * (lbounds(:) + ubounds(:))
    vgrid_canopy_q_assimi  => new_vgrid('q_canopy_layer', ZAXIS_GENERIC, config%ncanopy, &
      & longname = 'Canopy layers (leaf mass per layer as LAI) from Q_ASSIMI_ QUINCY', &
      & units = 'm2 m-2 LAI', &
      ! levels=   , &
      ! lbounds=  , &
      ! ubounds=  , &
      & dz = canopy_layer_thickness_profile(:))
    CALL register_vgrid(vgrid_canopy_q_assimi)

  END SUBROUTINE Init_q_assimi_config

#endif
END MODULE mo_q_assimi_config_class
