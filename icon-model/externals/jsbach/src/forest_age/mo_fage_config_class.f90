!> fage (forest age) config
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
!>#### Contains namelist info for the fage (forest age) proc
!>
MODULE mo_fage_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_util,               ONLY: int2string
  USE mo_io_units,           ONLY: filename_max
  USE mo_kind,               ONLY: wp
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN, MEDIUM_NAME_LEN

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_fage_config

  TYPE, EXTENDS(t_jsb_config) :: t_fage_config
    INTEGER              :: max_age                       !< Maximum age tracked
    INTEGER              :: nacs                          !< number of age classes
    CHARACTER(len=MEDIUM_NAME_LEN)    :: age_class_scheme
        !< scheme how age is distributed over classes (currently: shorterYoung or equal)
    CHARACTER(len=SHORT_NAME_LEN)     :: init_ac_fracts_scheme
        !< scheme how to initialise the age class fractions (so far only: allFirst)
    CHARACTER(len=MEDIUM_NAME_LEN)    :: disturbance_fracts_scheme
        !< scheme how to distribute disturbances over the age fractions in an age classes
        !< (so far only: prop -- i.e. proportionally)
    CHARACTER(len=MEDIUM_NAME_LEN)    :: alcc_age_fract_scheme
        !< scheme how to distribute alcc losses over the age classes and the age fractions
        !< within an age class (so far only: prop -- i.e. proportionally)
  CONTAINS
    PROCEDURE            :: Init => Init_fage_config
  END type t_fage_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize fage process
  !
  ! -------------------------------------------------------------------------------------------------------
  SUBROUTINE Init_fage_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_fage_config), INTENT(inout) :: config !< Configuration type of process (t_fage_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                       :: active
    CHARACTER(len=filename_max)   :: ic_filename, bc_filename

    INTEGER            :: nacs, max_age
    CHARACTER(len=MEDIUM_NAME_LEN)  :: age_class_scheme
    CHARACTER(len=SHORT_NAME_LEN)   :: init_ac_fracts_scheme
    CHARACTER(len=MEDIUM_NAME_LEN)  :: disturbance_fracts_scheme
    CHARACTER(len=MEDIUM_NAME_LEN)  :: alcc_age_fract_scheme

    NAMELIST /jsb_fage_nml/         &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & max_age,                    &
      & nacs,                       &
      & age_class_scheme,           &
      & init_ac_fracts_scheme,      &
      & disturbance_fracts_scheme,  &
      & alcc_age_fract_scheme

    INTEGER :: nml_handler, nml_unit, istat
    INTEGER :: i
    REAL(wp), ALLOCATABLE :: fage_grid_levels(:)
    REAL(wp), ALLOCATABLE :: age_class_max_age(:), age_class_mean_age(:)
    TYPE(t_jsb_vgrid), POINTER :: forest_age, age_classes_grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_fage_config'
    ! -------------------------------------------------------------------------------------------------- !

    CALL message(TRIM(routine), 'Starting fage configuration')

    ! Set defaults
    active            = .FALSE.
    bc_filename       = 'bc_land_fage.nc'  !JN-TODO do I need it at all?
    ic_filename       = 'ic_land_fage.nc'  !JN-TODO
    max_age           = 150
    nacs              = 11
    age_class_scheme  = "shorterYoung"     ! alternative: "equal"
    init_ac_fracts_scheme     = "allFirst"
    disturbance_fracts_scheme = "prop"
    alcc_age_fract_scheme     = "prop"

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_fage_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_fage_nml)
    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%max_age             = max_age
    config%nacs                = nacs
    config%age_class_scheme          = TRIM(age_class_scheme)
    config%init_ac_fracts_scheme     = TRIM(init_ac_fracts_scheme)
    config%disturbance_fracts_scheme = TRIM(disturbance_fracts_scheme)
    config%alcc_age_fract_scheme     = TRIM(alcc_age_fract_scheme)

    IF (.NOT. active) RETURN

    ALLOCATE(fage_grid_levels(max_age))
    ALLOCATE(age_class_max_age(nacs+1))
    ALLOCATE(age_class_mean_age(nacs))

    ! Assertion:
    IF ( max_age .LT. nacs) CALL finish(TRIM(routine), 'Assertion: number of age classes needs to be <= max age!')

    !> specifies vertical grid with age-class dimension
    ! since the boundaries of the age-class dimensions are static and the same all over the globe, they are now
    ! defined in this grid (previous implementations used the fage memory variable 'age_class_max_age' for this)
    CALL derive_age_class_boundaries(age_class_scheme, max_age, nacs, age_class_max_age, age_class_mean_age)

    age_classes_grid  => new_vgrid('age_classes_grid', ZAXIS_GENERIC, nacs, &
      & levels  = age_class_mean_age,          &
      & lbounds = age_class_max_age(1:nacs),   &
      & ubounds = age_class_max_age(2:nacs+1), &
      & units='y')
    CALL register_vgrid(age_classes_grid)

    !> specifies vertical grid with tracked age dimension
    fage_grid_levels(1:max_age) = (/ (REAL(i,kind=wp),i=1,max_age) /)
    forest_age  => new_vgrid('forest_age', ZAXIS_GENERIC, max_age, &
      & levels  = fage_grid_levels (1:max_age) - 0.5_wp, &
      & lbounds = fage_grid_levels (1:max_age) - 1._wp,  &
      & ubounds = fage_grid_levels (1:max_age),          &
      & units='y')
    CALL register_vgrid(forest_age)

    CALL message(TRIM(routine), 'Tracked forest ages [y]: '//TRIM(int2string(max_age)))
    CALL message(TRIM(routine), 'Run with n age classes: '//TRIM(int2string(nacs)))
    WRITE(message_text, *) 'Age classes in fage - upper boundary [y]: ', int(age_classes_grid%ubounds)
    CALL message(TRIM(routine), message_text)

    DEALLOCATE(fage_grid_levels, age_class_max_age, age_class_mean_age)

  END SUBROUTINE Init_fage_config

  ! ====================================================================================================== !
  !
  !> Derive the boundaries of the age-classes
  !
  SUBROUTINE derive_age_class_boundaries(age_class_scheme, max_age, nacs, age_class_max_age, age_class_mean_age)
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), INTENT(IN) ::  age_class_scheme
      !< the age class scheme to determine age class boundaries: currently either 'equal' or 'shorterYoung'
    INTEGER, INTENT(IN)          ::  max_age               !<
    INTEGER, INTENT(IN)          ::  nacs                  !<
    REAL(wp), INTENT(INOUT)      ::  age_class_max_age(:)  !<
    REAL(wp), INTENT(INOUT)      ::  age_class_mean_age(:) !<
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: i
    REAL(wp):: age_interval
    CHARACTER(len=*), PARAMETER :: routine = modname//':derive_age_class_boundaries'
    ! -------------------------------------------------------------------------------------------------- !

    ! Assertion: check expected vector dimensions
    IF (SIZE(age_class_max_age) .NE. (nacs + 1)) THEN
      CALL finish(TRIM(routine), 'Assertion: max age vector is assumed to have number of age classes + 1 as dimension!')
    END IF
    IF (SIZE(age_class_mean_age) .NE. nacs) THEN
      CALL finish(TRIM(routine), 'Assertion: mean age vector is assumed to have number of age classes as dimension!')
    END IF

    ! First class always contains the first year only
    age_class_max_age(1) = 0._wp
    age_class_max_age(2) = 1._wp

    IF (age_class_scheme == 'equal' ) THEN
      ! Equal intervals distributed equally over all but the first age class
      age_interval = (max_age - 1) / (nacs - 1)
      DO i = 2,nacs
        age_class_max_age(i + 1) = age_class_max_age (i) + int(age_interval)
      END DO
    ELSE IF (age_class_scheme == 'shorterYoung' ) THEN
      ! Equal intervals distributed unequally such that younger age classes have shorter steps
      ! (1, 1+n, 1+n+2*n, 1+n+2*n+3*n, ...)
      age_interval = 0
      DO i = 1,nacs-1
        age_interval = age_interval + i
      END DO
      age_interval = max_age / age_interval
      DO i = 2,nacs-1
        age_class_max_age(i + 1) = age_class_max_age(i) + int(i * age_interval)
      END DO
    ELSE
      CALL finish(TRIM(routine), 'Currently only "shorterYoung" and "equal" implemented as forest ageing schemes')
    END IF

    ! Last age class is not really restricted to a max age, but set to max_age for compliance with the fract_per_age field
    age_class_max_age(nacs+1) = max_age

    ! Assertion (currently this is of course the case - just in case that others want to define distributions themselves)
    IF (age_class_max_age(2) .NE. 1._wp) THEN
      CALL finish(TRIM(routine), 'Assertion: upper boundary of lowest age class always needs to have max_age of 1!')
    END IF
    IF (ANY(age_class_max_age(:) .GT. max_age)) THEN
      CALL finish(TRIM(routine), 'Assertion: no age class should have a max age larger than the total max_age!')
    END IF

    DO i = 1,nacs
      age_class_mean_age(i) = (age_class_max_age(i) + age_class_max_age(i+1)) * 0.5_wp
    END DO

  END SUBROUTINE derive_age_class_boundaries

#endif
END MODULE mo_fage_config_class
