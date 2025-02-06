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

! @brief configuration setup for Grib output
!
! configuration setup for tracer transport

MODULE mo_gribout_config

  USE mo_impl_constants,     ONLY: max_phys_dom
  USE mo_exception,          ONLY: message
  USE mo_grib2_tile,         ONLY: t_grib2_template_tile
  USE mo_cdi,                ONLY: gribapiLibraryVersion

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: t_gribout_config
  PUBLIC :: gribout_config
  PUBLIC :: configure_gribout
  PUBLIC :: gribout_crosscheck

  PUBLIC :: GRIB_UNDEFVAL
  PUBLIC :: GRIB_NOTINUSEVAL
  PUBLIC :: GRIB_LIB_COMPAT_ECC_2_31_0

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_gribout_config'

  INTEGER, PARAMETER :: GRIB_UNDEFVAL    = -1
  INTEGER, PARAMETER :: GRIB_NOTINUSEVAL = 0

  INTEGER, PARAMETER :: GRIB_LIB_COMPAT_ECC_2_31_0 = 1

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for grib output
  !!--------------------------------------------------------------------------
  TYPE :: t_gribout_config

    ! namelist variables

    INTEGER :: tablesVersion              ! Main switch for table version
    INTEGER :: localTablesVersion         ! Switch for centre local table version

    INTEGER :: &                          ! Table 1.2
      & significanceOfReferenceTime       ! 0: Analysis
                                          ! 1: Start of forecast
                                          ! 2: Verifying time of forecast
                                          ! 4: ...

    INTEGER :: &                          ! Table 1.3
      & productionStatusOfProcessedData   ! 0: Oper. products
                                          ! 1: Oper. test products
                                          ! 2: Research products
                                          ! 3: ...

    INTEGER :: &                          ! Table 1.4
      & typeOfProcessedData               ! 0: Analysis products
                                          ! 1: Forecast products
                                          ! 2: Analysis and forecast products
                                          ! 3: ...

    INTEGER :: &                          ! Table 4.3
      & typeOfGeneratingProcess           ! 0: Analysis
                                          ! 1: Initialization
                                          ! 2: Forecast
                                          ! 3: ...
                                          ! 196: invariant data

    INTEGER :: &                          ! Table: backgroundProcess
      & backgroundProcess                 ! 0: main run
                                          ! 1: pre-assimilation
                                          ! 2: assimilation
                                          ! 3: ...

    INTEGER :: &                          ! Table: generatingProcessIdentifier
      & generatingProcessIdentifier       ! 1: icogl
                                          ! 2: icrgl
                                          ! 3: icoeu
                                          ! 4: ...

    INTEGER :: &                          ! Table: local.78.254.def
      & localDefinitionNumber             ! 252: Ensemble system incl. postprocessing
                                          ! 253: Ensemble system
                                          ! 254: Deterministic system

    INTEGER :: &                          ! Table: local.78.254.def
      & localNumberOfExperiment           !


    INTEGER :: &                          ! Output generating center
      & generatingCenter                  ! DWD  : 78
                                          ! ECMWF: 98


    INTEGER :: &                          ! Output generating subcenter
      & generatingSubcenter               ! DWD  : 255
                                          ! ECMWF: 0

    LOGICAL :: lspecialdate_invar         ! .TRUE.: use special date 10101 for encoding
                                          ! invariant and climatological fields
                                          ! .FALSE.: no special treatment of invariant
                                          ! and climatological fields.

    LOGICAL :: ldate_grib_act             ! add Creation date to GRIB file
                                          ! .TRUE. : activated
                                          ! .FALSE.: deactivated (use dummy date/time)

    LOGICAL :: lgribout_24bit             ! write thermodynamic fields rho, theta_v, T, p
                                          ! with 24bit precision

    LOGICAL :: lgribout_compress_ccsds    ! enable CCSDS second level compression

    INTEGER :: typeOfEnsembleForecast,        &
      &        localTypeOfEnsembleForecast,   &
      &        numberOfForecastsInEnsemble,   &
      &        perturbationNumber

    CHARACTER(len=32) ::  &               !< type of GRIB2 templates used for surface tile fields
      &  typeOfGrib2TileTemplate          !  'wmo': official WMO templates 55, 59, ...
                                          !  'dwd': local 'DWD' templates 40455, 40456, ...

    INTEGER :: grib_lib_compat            !< Type of GRIB library backward compatibility adjustment:
                                          !< 0: none
                                          !< 1: for ecCodes versions >= 2.32.0

    ! derived variables
    !
    TYPE(t_grib2_template_tile)::  &      !< defines set of employed GRIB2 tile templates for writing
      &  grib2_template_tile              !< contains tile template numbers as well as the corresponding
                                          !< set GRIB2 tile keys.
                                          !< variant 1: templates 40455, 40456, etc
                                          !<            DWD's local tile templates
                                          !< variant 2: 55, 59, etc
                                          !< official WMO tile templates
                                          !< the variant in use is determined by the namelist parameter itype_tiletemplate


  END TYPE t_gribout_config

  !>
  !!
  TYPE(t_gribout_config), TARGET :: gribout_config(1:max_phys_dom)


CONTAINS


  !! potentially modify generatingCenter and generatingSubcenter
  !!
  !! If generatingCenter and generatingSubcenter are not set via namelist,
  !! they are filled with values read from the grid file.
  !!
  SUBROUTINE configure_gribout(grid_generatingCenter, grid_generatingSubcenter, &
    &                          n_dom)
  !
    INTEGER,       INTENT(IN)  :: grid_generatingCenter(0:)
    INTEGER,       INTENT(IN)  :: grid_generatingSubcenter(0:)
    INTEGER,       INTENT(IN)  :: n_dom

    ! local fields
    INTEGER  :: jg

    !-----------------------------------------------------------------------

    DO jg = 1, n_dom
      !
      ! check, whether generatingCenter was set in gribout_nml
      !
      IF ( gribout_config(jg)%generatingCenter == GRIB_UNDEFVAL ) THEN
        ! If not, then fill with grid generating center
        gribout_config(jg)%generatingCenter = grid_generatingCenter(jg)
      ENDIF

      ! check, whether generatingSubcenter was set in gribout_nml
      !
      IF ( gribout_config(jg)%generatingSubcenter == GRIB_UNDEFVAL ) THEN
        ! If not, then fill with grid generating subcenter
        gribout_config(jg)%generatingSubcenter = grid_generatingSubcenter(jg)
      ENDIF

    ENDDO  ! jg

  END SUBROUTINE configure_gribout

  !> Crosscheck with settings beyond gribout_nml
  !!
  SUBROUTINE gribout_crosscheck(n_dom, verbose)

    ! Arguments
    INTEGER, INTENT(IN) :: n_dom
    LOGICAL, INTENT(IN) :: verbose

    ! Local variables
    INTEGER :: eccodes_version(3)
    INTEGER :: jg
    LOGICAL :: is_ecc_vers_ge_2_32_0

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::gribout_crosscheck'

    !-----------------------------------------------------------------------

    ! Inquire ecCodes version
    CALL gribapiLibraryVersion(eccodes_version(1), eccodes_version(2), eccodes_version(3))

    ! For GRIB library backward compatibility: gribout_nml / grib_lib_compat:

    ! The compatibility adjustment is only necessary for ecCodes versions >= 2.32.0
    is_ecc_vers_ge_2_32_0 = ( version_compare(eccodes_version, [2, 32, 0]) >= 0 )

    IF (.NOT. is_ecc_vers_ge_2_32_0) THEN
      DO jg = 1, n_dom
        gribout_config(jg)%grib_lib_compat = GRIB_NOTINUSEVAL
      END DO
      IF (verbose) CALL message(routine, "ecCodes version < 2.32.0: Compatiblity adjustment switched off.")
    ENDIF

  END SUBROUTINE gribout_crosscheck

  !> Auxiliary function for (three-part) version number comparison
  !!
  INTEGER FUNCTION version_compare(a, b)

    INTEGER, INTENT(IN) :: a(3)
    INTEGER, INTENT(IN) :: b(3)

    INTEGER :: i

    !---------------------------------

    DO i = 1, 3
      IF (a(i) /= b(i)) THEN
        version_compare = a(i) - b(i)
        RETURN
      ENDIF
    END DO

    version_compare = 0

  END FUNCTION version_compare

END MODULE mo_gribout_config
