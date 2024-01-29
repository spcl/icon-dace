! This module contains types for polymorphic tracer metadata. Additionally,
! a type bound procedure (TBP) is given as a constructor for the base type
! t_tracer_meta. This procedure is private and may only be adressed as TBP.
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_tracer_metadata_types

  USE mo_kind,            ONLY: wp
  USE mo_key_value_store, ONLY: t_key_value_store

  IMPLICIT NONE

  PRIVATE

  ! maximum string length for variable names
  INTEGER, PARAMETER :: VARNAME_LEN = 32

  ! Tracer metadata
  ! 
  ! Polymorphic type that contains tracer metadata according to the type 
  ! of tracer (aerosol, chemical, hydrometeor). 
  TYPE t_tracer_meta
    !
    LOGICAL :: lis_tracer         ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=VARNAME_LEN) :: name         ! Name of tracer
    ! One-way vs. Two-way nesting
    LOGICAL :: lfeedback          ! feedback from child- to parent domain (TRUE/FALSE)
    ! Advection
    INTEGER :: ihadv_tracer       ! Method for horizontal transport
    INTEGER :: ivadv_tracer       ! Method for vertical transport
    ! Other processes covered by ICON
    LOGICAL :: lturb_tracer       ! Turbulent transport (TRUE/FALSE)
    LOGICAL :: lconv_tracer       ! Convection  (TRUE/FALSE)
    ! Processes not covered by ICON (requires ART extension)
    TYPE(t_key_value_store) :: opt_meta   ! Storage container for optional metadata
    
    CONTAINS
      procedure :: construct_base => construct_t_tracer_meta
  END TYPE t_tracer_meta

  ! Aerosol-specific metadata
  TYPE, extends(t_tracer_meta) :: t_aero_meta
    ! Non-optional metadata
    INTEGER  :: moment            ! moment of distribution (e.g. 0=number, 3=proportional to mass)
    CHARACTER(LEN=VARNAME_LEN) :: &
      &         mode              ! name of mode the tracer is contained in
    CHARACTER(LEN=VARNAME_LEN) :: &
      &         substance         ! to which substance this tracer belongs
    REAL(wp) :: rho               ! Density [kg m-3]
    REAL(wp) :: mol_weight        ! Molar mass [g mol-1]
    LOGICAL  :: linsol            ! TRUE for insoluble tracers

    INTEGER :: ised_tracer        ! Sedimentation
                                  !   0 = No sedimentation
                                  !   1 = Monodisperse aerosol
                                  !   2 = As part of an according aerosol mode
                                  !   ... (hydrometeors to be added)
    LOGICAL :: ldep_tracer        ! Dry deposition  (TRUE/FALSE)
    INTEGER :: iwash_tracer       ! Washout
                                  !   0 = No washout
                                  !   1 = Monodisperse aerosol
                                  !   2 = As part of an according aerosol mode
  END TYPE

  ! Chemical tracer metadata
  TYPE, extends(t_tracer_meta) :: t_chem_meta
    ! Non-optional metadata
    INTEGER :: ised_tracer        ! Sedimentation
                                  !   0 = No sedimentation
                                  !   1 = Monodisperse aerosol
                                  !   2 = As part of an according aerosol mode
    LOGICAL :: ldep_tracer        ! Dry deposition  (TRUE/FALSE)
    INTEGER :: iwash_tracer       ! Washout
                                  !   0 = No washout
                                  !   1 = Monodisperse aerosol
                                  !   2 = As part of an according aerosol mode

    CONTAINS
        PROCEDURE, PASS(this) :: set_tracer_meta => create_tracer_metadata_chem
  END TYPE

  ! Hydrometeor metadata
  TYPE, extends(t_tracer_meta) :: t_hydro_meta
    ! 
  END TYPE

  PUBLIC :: t_tracer_meta
  PUBLIC :: t_aero_meta
  PUBLIC :: t_chem_meta
  PUBLIC :: t_hydro_meta

CONTAINS

  SUBROUTINE construct_t_tracer_meta(meta, lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
    &                                lturb_tracer, lconv_tracer)

    CLASS(t_tracer_meta),INTENT(OUT)     :: meta          ! Base meta container to be filled
    LOGICAL, INTENT(IN), OPTIONAL        :: lis_tracer    ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name          ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL        :: lfeedback     ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL        :: ihadv_tracer  ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL        :: ivadv_tracer  ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL        :: lturb_tracer  ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL        :: lconv_tracer  ! Switch for convection

    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      meta%lis_tracer = lis_tracer
    ELSE
      meta%lis_tracer = .FALSE.
    ENDIF

    ! name
    IF ( PRESENT(name) ) THEN
      meta%name = name
    ELSE
      meta%name = "unnamed"
    ENDIF

    ! lfeedback
    IF ( PRESENT(lfeedback) ) THEN
      meta%lfeedback = lfeedback
    ELSE
      meta%lfeedback = .FALSE.
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      meta%ihadv_tracer = ihadv_tracer
    ELSE
      meta%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      meta%ivadv_tracer = ivadv_tracer
    ELSE
      meta%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer
    IF ( PRESENT(lturb_tracer) ) THEN
      meta%lturb_tracer = lturb_tracer
    ELSE
      meta%lturb_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      meta%lconv_tracer = lconv_tracer
    ELSE
      meta%lconv_tracer = .FALSE.
    ENDIF

  END SUBROUTINE construct_t_tracer_meta


  SUBROUTINE create_tracer_metadata_chem(this, lis_tracer, name, lfeedback, ihadv_tracer,     &
                      &                  ivadv_tracer, lturb_tracer, lconv_tracer,            &
                      &                  ised_tracer, ldep_tracer, iwash_tracer)
    ! Base type (t_tracer_meta) content
    CLASS(t_chem_meta), INTENT(INOUT)  :: this
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name       ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer      ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer      ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer     ! Method for washout

    ! Fill the metadata of the base type
    CALL this%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                             lturb_tracer, lconv_tracer)

    ! Fill the meta of the extended type (t_chem_meta)
    ! ised_tracer
    IF ( PRESENT(ised_tracer) ) THEN
      this%ised_tracer = ised_tracer
    ELSE
      this%ised_tracer = 0
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      this%ldep_tracer = ldep_tracer
    ELSE
      this%ldep_tracer = .FALSE.
    ENDIF

    ! iwash_tracer
    IF ( PRESENT(iwash_tracer) ) THEN
      this%iwash_tracer = iwash_tracer
    ELSE
      this%iwash_tracer = 0
    ENDIF


END SUBROUTINE create_tracer_metadata_chem

END MODULE mo_tracer_metadata_types
