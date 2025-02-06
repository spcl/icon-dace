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

! Namelist for Grib output
!
! These subroutines are called by  read_atmo_namelists and do the transport
! setup.

MODULE mo_gribout_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_gribout_config,      ONLY: gribout_config, GRIB_UNDEFVAL, GRIB_NOTINUSEVAL, &
    &                               GRIB_LIB_COMPAT_ECC_2_31_0
  USE mo_grib2_tile,          ONLY: grib2_keys_tile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_util_string,         ONLY: int2string, tolower, one_of
  USE mo_exception,           ONLY: finish, message


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_gribout_namelist


  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_gribout_nml'

  !----------------------------------!
  ! gribout_nml namelist variables   !
  !----------------------------------!

  ! Namelist parameter "preset": main switch (string), possible options are
  !  - "none"
  !  - "deterministic"
  !  - "ensemble"
  !
  ! Setting this different to "none" enables a couple of defaults for
  ! the other gribout_nml namelist parameters. If, additionally, the
  ! user tries to set any of these other parameters to a conflicting
  ! value, an error message is thrown.
  CHARACTER(LEN=32) :: preset

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
    & typeOfGeneratingProcess           ! 0  : Analysis
                                        ! 1  : Initialization
                                        ! 2  : Forecast
                                        ! 3  : ...
                                        ! 196: invariant data

  INTEGER :: &                          ! Table: backgroundProcess
    & backgroundProcess                 ! 0: main run
                                        ! 1: pre-assimilation
                                        ! 2: assimilation
                                        ! 3: ...

  INTEGER :: &                          ! Table: generatingProcessIdentifier
    & generatingProcessIdentifier(1:max_dom) ! 1: icogl
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
    & generatingCenter                  !


  INTEGER :: &                          ! Output generating subcenter
    & generatingSubcenter               !


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

  CHARACTER(LEN=22) :: &                ! Type of GRIB library backward compatibility adjustment:
    & grib_lib_compat                   ! 'current'
                                        ! 'eccodes:2.31.0'

  NAMELIST/gribout_nml/  &
    &                    preset, tablesVersion,           &
    &                    localTablesVersion,              &
    &                    significanceOfReferenceTime,     &
    &                    productionStatusOfProcessedData, &
    &                    typeOfProcessedData,             &
    &                    typeOfGeneratingProcess,         &
    &                    backgroundProcess,               &
    &                    generatingProcessIdentifier,     &
    &                    localDefinitionNumber,           &
    &                    localNumberOfExperiment,         &
    &                    generatingCenter,                &
    &                    generatingSubcenter,             &
    &                    lspecialdate_invar,              &
    &                    ldate_grib_act,                  &
    &                    typeOfEnsembleForecast,          &
    &                    localTypeOfEnsembleForecast,     &
    &                    numberOfForecastsInEnsemble,     &
    &                    perturbationNumber,              &
    &                    lgribout_24bit,                  &
    &                    lgribout_compress_ccsds,         &
    &                    typeOfGrib2TileTemplate,         &
    &                    grib_lib_compat


CONTAINS


  !-------------------------------------------------------------------------
  !
  !! Read Namelist for gribout.
  !!
  !! This subroutine
  !! - reads the Namelist for gribout
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_gribout_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit
    INTEGER :: grib_lib_compat_int

    CHARACTER(len=*), PARAMETER :: routine = modname//"::read_gribout_nml"

    !-----------------------
    ! 1. default settings
    !-----------------------
    preset                               = "deterministic"

    tablesVersion                        = 15
    localTablesVersion                   = 1

    significanceOfReferenceTime          = 1   ! 1: Start of forecast
    productionStatusOfProcessedData      = 1   ! 1: Oper. test products
    backgroundProcess                    = 0   ! 0: main run
    generatingProcessIdentifier(:)       = 1   ! 1: icogl
    localNumberOfExperiment              = 1
    lspecialdate_invar                   = .FALSE.  ! no special date for invar fields
    ldate_grib_act                       = .TRUE.
    lgribout_24bit                       = .FALSE.  ! use 16bit precision for all fields
    lgribout_compress_ccsds              = .FALSE.  ! do not use second level compression by default

    typeOfGrib2TileTemplate              = 'dwd'    ! use DWD templates 40455, etc

    typeOfProcessedData                  = GRIB_UNDEFVAL
    typeOfGeneratingProcess              = GRIB_UNDEFVAL
    localDefinitionNumber                = GRIB_UNDEFVAL
    generatingCenter                     = GRIB_UNDEFVAL  ! output generating center
    generatingSubcenter                  = GRIB_UNDEFVAL  ! output generating subcenter
    typeOfEnsembleForecast               = GRIB_UNDEFVAL  ! (undefined, will not be set if unchanged)
    localTypeOfEnsembleForecast          = GRIB_UNDEFVAL  ! (undefined, will not be set if unchanged)
    numberOfForecastsInEnsemble          = GRIB_UNDEFVAL  ! (undefined, will not be set if unchanged)
    perturbationNumber                   = GRIB_UNDEFVAL  ! (undefined, will not be set if unchanged)

    grib_lib_compat                      = 'current'      ! i.e. switched off

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('gribout_nml')
      READ(funit,NML=gribout_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('gribout_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, gribout_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, gribout_nml)                                      ! overwrite default settings
      ! Preset values, when main switch is provided.
      CALL preset_namelist()
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, gribout_nml)  ! write settings to temporary text file
      END IF
    CASE DEFAULT
      CALL preset_namelist()
    END SELECT
    CALL close_nml


    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    IF ( one_of(tolower(TRIM(typeOfGrib2TileTemplate)), (/'wmo','dwd'/)) < 0 ) THEN
      CALL finish(routine, "Illegal Namelist setting for typeOfGrib2TileTemplate: "//TRIM(typeOfGrib2TileTemplate))
    ENDIF

    ! Currently, the following types of GRIB library backward compatibility adjustment are available:
    ! - 'current': Switched off.
    ! - 'eccodes:2.31.0': For ecCodes versions >= 2.32.0.
    !    CDI uses the ecCodes sample file "GRIB2.tmpl" as a starting file for ecCodes.
    !    The SecondFixedSurface keys are assigned the following values in GRIB2.tmpl:
    !
    !    - typeOfSecondFixedSurface        = 255 (Missing)
    !    - scaleFactorOfSecondFixedSurface = MISSING (= 255)
    !    - scaledValueOfSecondFixedSurface = MISSING (= 2,147,483,647)
    !
    !    Now, if typeOfSecondFixedSurface is set to a value >= 10, 102 say,
    !    but scaleFactorOfSecondFixedSurface and scaledValueOfSecondFixedSurface
    !    are not explicitly set, the result is as follows:
    !
    !    For ecCodes version < 2.32.0:
    !    -----------------------------
    !    - typeOfSecondFixedSurface        = 102 (Specific altitude above mean sea level)
    !    - scaleFactorOfSecondFixedSurface = 0
    !    - scaledValueOfSecondFixedSurface = 0
    !
    !    For ecCodes version >= 2.32.0:
    !    ------------------------------
    !    - typeOfSecondFixedSurface        = 102 (Specific altitude above mean sea level)
    !    - scaleFactorOfSecondFixedSurface = MISSING
    !    - scaledValueOfSecondFixedSurface = MISSING
    !
    !    With grib_lib_compat = 'eccodes:2.31.0', we try to overwrite the behavior of
    !    ecCodes versions >= 2.32.0 with the behavior of versions < 2.32.0, in this respect.
    !
    grib_lib_compat_int = GRIB_NOTINUSEVAL
    IF (LEN_TRIM(grib_lib_compat) > 0) THEN
      ! lowercase namelist entry
      grib_lib_compat = tolower(grib_lib_compat)
      IF (TRIM(grib_lib_compat) == 'current') THEN
        grib_lib_compat_int = GRIB_NOTINUSEVAL
      ELSEIF (TRIM(grib_lib_compat) == 'eccodes:2.31.0') THEN
        grib_lib_compat_int = GRIB_LIB_COMPAT_ECC_2_31_0
      ELSE
        ! invalid namelist entry
        CALL finish(routine, "Invalid namelist setting for grib_lib_compat: "//TRIM(grib_lib_compat))
      ENDIF
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 1,max_dom
      gribout_config(jg)%tablesVersion                     = &
        &                tablesVersion
      gribout_config(jg)%localTablesVersion                = &
        &                localTablesVersion
      gribout_config(jg)%significanceOfReferenceTime       = &
        &                significanceOfReferenceTime
      gribout_config(jg)%productionStatusOfProcessedData   = &
        &                productionStatusOfProcessedData
      gribout_config(jg)%typeOfProcessedData               = &
        &                typeOfProcessedData
      gribout_config(jg)%typeOfGeneratingProcess           = &
        &                typeOfGeneratingProcess
      gribout_config(jg)%backgroundProcess                 = &
        &                backgroundProcess
      gribout_config(jg)%generatingProcessIdentifier       = &
        &                generatingProcessIdentifier(jg)
      gribout_config(jg)%localDefinitionNumber             = &
        &                localDefinitionNumber
      gribout_config(jg)%localNumberOfExperiment           = &
        &                localNumberOfExperiment
      gribout_config(jg)%generatingCenter                  = &
        &                generatingCenter
      gribout_config(jg)%generatingSubcenter               = &
        &                generatingSubcenter
      gribout_config(jg)%lspecialdate_invar                = &
        &                lspecialdate_invar
      gribout_config(jg)%ldate_grib_act                    = &
        &                ldate_grib_act
      gribout_config(jg)%typeOfEnsembleForecast            = &
        &                typeOfEnsembleForecast
      gribout_config(jg)%localTypeOfEnsembleForecast       = &
        &                localTypeOfEnsembleForecast
      gribout_config(jg)%numberOfForecastsInEnsemble       = &
        &                numberOfForecastsInEnsemble
      gribout_config(jg)%perturbationNumber                = &
        &                perturbationNumber
      gribout_config(jg)%lgribout_24bit                    = &
        &                lgribout_24bit
      gribout_config(jg)%lgribout_compress_ccsds           = &
        &                lgribout_compress_ccsds
      gribout_config(jg)%typeOfGrib2TileTemplate           = &
        &                tolower(typeOfGrib2TileTemplate)
      gribout_config(jg)%grib_lib_compat                   = &
        &                grib_lib_compat_int
    ENDDO



    !----------------------------------------------------------------
    ! 5b. Define the set of employed GRIB2 tile templates for writing
    !----------------------------------------------------------------
    !
    ! By placing it here we make sure that this information is also
    ! known to the output PEs, which are detached right after.
    !
    DO jg = 1, max_dom
      IF (TRIM(gribout_config(jg)%typeOfGrib2TileTemplate)=="wmo") THEN
        !
        ! allowed set of templates to choose from
        !
        gribout_config(jg)%grib2_template_tile%tpl_inst     = 55
        gribout_config(jg)%grib2_template_tile%tpl_acc      = -999  ! 62 not yet available (validation by WMO pending)
        gribout_config(jg)%grib2_template_tile%tpl_inst_ens = 59
        gribout_config(jg)%grib2_template_tile%tpl_acc_ens  = -999  ! 63 not yet available (validation by WMO pending)
        !
        ! keynames for official WMO templates 55, 59, ...
        !
        gribout_config(jg)%grib2_template_tile%keys = &
          &                grib2_keys_tile(str_tileClassification              = "tileClassification",              &
          &                                str_totalNumberOfTileAttributePairs = "totalNumberOfTileAttributePairs", &
          &                                str_numberOfUsedSpatialTiles        = "numberOfUsedSpatialTiles",        &
          &                                str_tileIndex                       = "tileIndex",                       &
          &                                str_numberOfUsedTileAttributes      = "numberOfUsedTileAttributes",      &
          &                                str_attributeOfTile                 = "attributeOfTile"                  )
      ELSE
        !
        ! allowed set of templates to choose from
        !
        gribout_config(jg)%grib2_template_tile%tpl_inst     = 40455
        gribout_config(jg)%grib2_template_tile%tpl_acc      = -999  ! not available (never)
        gribout_config(jg)%grib2_template_tile%tpl_inst_ens = 40456
        gribout_config(jg)%grib2_template_tile%tpl_acc_ens  = -999  ! not available (never)
        !
        ! keynames for local 'DWD' templates 40455, 40456, ...
        !
        gribout_config(jg)%grib2_template_tile%keys =  &
          &                grib2_keys_tile(str_tileClassification              = "tileClassification",              &
          &                                str_totalNumberOfTileAttributePairs = "totalNumberOfTileAttributePairs", &
          &                                str_numberOfUsedSpatialTiles        = "numberOfTiles",                   &
          &                                str_tileIndex                       = "tileIndex",                       &
          &                                str_numberOfUsedTileAttributes      = "numberOfTileAttributes",          &
          &                                str_attributeOfTile                 = "tileAttribute"                    )
      ENDIF
    ENDDO


    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=gribout_nml)
      CALL store_and_close_namelist(funit, 'gribout_nml')
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=gribout_nml)


  END SUBROUTINE read_gribout_namelist


  !> Preset values, when main switch is provided.
  !
  SUBROUTINE preset_namelist()
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::preset_namelist"

    IF (TRIM(preset) == "deterministic") THEN
      !
      ! deterministic forecast
      !
      ! 2: Forecast
      CALL preset_value("typeOfGeneratingProcess",typeOfGeneratingProcess,             2 , quiet=.TRUE.)
      ! 254: Deterministic system
      CALL preset_value("localDefinitionNumber",  localDefinitionNumber,             254 , quiet=.TRUE. )
      ! 1: Forecast products
      CALL preset_value("typeOfProcessedData",    typeOfProcessedData,                 1 , quiet=.TRUE. )

    ELSE IF (TRIM(preset) == "ensemble") THEN
      !
      ! ensemble forecast
      !
      ! values are preset according to
      !  GME   : pp_makepdt.f90
      !  COSMO : io_metadata.f90
      !
      ! The user is expected to set only
      !
      ! - perturbationNumber          (COSMO: "iepsmem")
      ! - numberOfForecastsInEnsemble (COSMO: "iepstot")
      ! - localTypeOfEnsembleForecast (COSMO: "iepstyp")
      !
      CALL preset_value("typeOfGeneratingProcess",         typeOfGeneratingProcess,             4 , quiet=.FALSE. )
      ! 253: ensemble
      CALL preset_value("localDefinitionNumber",           localDefinitionNumber,             253 , quiet=.FALSE. )
      ! 5  : "control and perturbed forecast products"
      CALL preset_value("typeOfProcessedData",             typeOfProcessedData,                 5 , quiet=.FALSE. )
      ! Note: atmospheric chemical constituents -> 41
      !       statistically processed data      -> 11
      CALL preset_value("typeOfEnsembleForecast",          typeOfEnsembleForecast,            192 , quiet=.FALSE. )

    ELSE IF (LEN_TRIM(preset) > 0) THEN
      !
      ! invalid namelist entry for preset
      !
      CALL finish(routine, "Invalid namelist setting for preset: "//TRIM(preset))

    END IF
  END SUBROUTINE preset_namelist


  !> Sets a variable to a given value, but only if current value was
  !> GRIB_UNDEFVAL = -1
  !
  SUBROUTINE preset_value(name, ival, preset_val, quiet)
    CHARACTER(LEN=*), INTENT(IN)    :: name              !< name of this parameter (for screen output)
    INTEGER,          INTENT(INOUT) :: ival              !< value to be altered
    INTEGER,          INTENT(IN)    :: preset_val        !< value that shall be preset
    LOGICAL,          INTENT(IN)    :: quiet             !< LOGICAL: if .FALSE. we print some screen output
    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//"::preset_value"
    CHARACTER(len=20) :: cval, cpre

    IF (ival == GRIB_UNDEFVAL) THEN
      ival = preset_val
      IF (.NOT. quiet) THEN
        CALL message(routine, "presetting namelist parameter '"//TRIM(name)//"' as "//TRIM(int2string(preset_val)))
      END IF
    ELSE
      IF (ival /= preset_val) THEN
        ! obviously, the user tried to set both: the main switch for
        ! presetting values and the local value
        WRITE (cval,'(i0)') ival
        WRITE (cpre,'(i0)') preset_val
        CALL finish(routine, "Namelist setting "//TRIM(name)//"="//TRIM(cval)// &
                             " contradicts preset value "//TRIM(cpre))
      END IF
    END IF

  END SUBROUTINE preset_value

END MODULE mo_gribout_nml
