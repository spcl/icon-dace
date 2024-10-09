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

! This module contains the constructors for polymorphic tracer metadata.
! For each extended object, a constructor for the base type is called
! via a type bound procedure.

MODULE mo_tracer_metadata

  USE mo_kind,                  ONLY: wp
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta, t_aero_meta, &
                                  &   t_chem_meta, t_hydro_meta

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: create_tracer_metadata
  PUBLIC  :: create_tracer_metadata_aero
  PUBLIC  :: create_tracer_metadata_hydro

  !------------------------------------------------------------------------------------------------
  !>create_tracer_metadata
  ! Quasi-constructors for tracer metadata
  ! (public routine. Can be used for two things:
  ! 1.) default settings: If used without any argument, the routines return a variable
  !     of the according extended type, containing the default settings.
  ! 2.) Setting of metadata: If used with arguments, the routines return a variable
  !     of the according extended type, containing the default settings except for those components
  !     which are given in the argument list.

CONTAINS

  TYPE(t_tracer_meta) FUNCTION create_tracer_metadata(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
                      &                               lturb_tracer, lconv_tracer)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name       ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection

    ! Fill the metadata of the base type
    CALL create_tracer_metadata%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                        lturb_tracer, lconv_tracer)

  END FUNCTION create_tracer_metadata



  TYPE(t_aero_meta) FUNCTION create_tracer_metadata_aero(lis_tracer, name, lfeedback, ihadv_tracer,      &
                      &                                  ivadv_tracer, lturb_tracer, lconv_tracer, ised_tracer, &
                      &                                  ldep_tracer, iwash_tracer, moment, mode, substance,   &
                      &                                  rho, mol_weight, linsol)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer      ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback       ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer    ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer    ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer    ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer    ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer     ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer     ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer    ! Method for washout
    ! Extended type (t_aero_meta) content
    INTEGER, INTENT(IN), OPTIONAL  :: moment          ! moment of distribution (e.g. 0=number, 3=proportional to mass)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: &
      &                               mode             ! Name of mode the tracer is contained in
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: substance  ! substance the tracer represents
    REAL(wp), INTENT(IN), OPTIONAL :: rho              ! Density [kg m-3]
    REAL(wp), INTENT(IN), OPTIONAL :: mol_weight       ! Molar mass [g mol-1]
    LOGICAL,  INTENT(IN), OPTIONAL :: linsol           ! TRUE for insoluble tracer

    ! Fill the metadata of the base type
    CALL create_tracer_metadata_aero%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                             lturb_tracer, lconv_tracer)


    ! Fill the meta of the extended type (t_aero_meta)
    IF(PRESENT(moment)) THEN
      create_tracer_metadata_aero%moment = moment
    ELSE
      create_tracer_metadata_aero%moment = -1
    ENDIF

    ! Fill the meta of the extended type (t_aero_meta)
    IF(PRESENT(mode)) THEN
      create_tracer_metadata_aero%mode = TRIM(mode)
    ELSE
      create_tracer_metadata_aero%mode = 'no_mode'
    ENDIF

    ! Fill the meta of the extended type (t_aero_meta)
    IF(PRESENT(substance)) THEN
      create_tracer_metadata_aero%substance = TRIM(substance)
    ELSE
      create_tracer_metadata_aero%substance = 'no_substance'
    ENDIF

    IF(PRESENT(rho)) THEN
      create_tracer_metadata_aero%rho = rho
    ELSE
      create_tracer_metadata_aero%rho = -1._wp
    ENDIF

    IF(PRESENT(mol_weight)) THEN
      create_tracer_metadata_aero%mol_weight = mol_weight
    ELSE
      create_tracer_metadata_aero%mol_weight = -1._wp
    ENDIF

    IF(PRESENT(linsol)) THEN
      create_tracer_metadata_aero%linsol = linsol
    ELSE
      create_tracer_metadata_aero%linsol = .FALSE.
    ENDIF

    ! ised_tracer
    IF ( PRESENT(ised_tracer) ) THEN
      create_tracer_metadata_aero%ised_tracer = ised_tracer
    ELSE
      create_tracer_metadata_aero%ised_tracer = 0
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      create_tracer_metadata_aero%ldep_tracer = ldep_tracer
    ELSE
      create_tracer_metadata_aero%ldep_tracer = .FALSE.
    ENDIF

    ! iwash_tracer
    IF ( PRESENT(iwash_tracer) ) THEN
      create_tracer_metadata_aero%iwash_tracer = iwash_tracer
    ELSE
      create_tracer_metadata_aero%iwash_tracer = 0
    ENDIF

  END FUNCTION create_tracer_metadata_aero


  TYPE(t_hydro_meta) FUNCTION create_tracer_metadata_hydro(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
                      &                                  lturb_tracer, lconv_tracer)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    ! Extended type (t_hydro_meta) content
    ! ...

    ! Fill the metadata of the base type
    CALL create_tracer_metadata_hydro%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                              lturb_tracer, lconv_tracer)

    ! Fill the meta of the extended type (t_hydro_meta)
    ! ...

  END FUNCTION create_tracer_metadata_hydro

END MODULE mo_tracer_metadata

