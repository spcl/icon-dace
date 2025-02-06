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

! A list of requests for input.

! The use case for which this class has been developed is this:
!  1. A list is created with `myList = t_InputRequestList()`, and the names (only the names) of the requested fields are added with `request()`.
!  2. Files are read with readFile(). This reads all DATA ASSOCIATED with the requested variable names into memory (already distributing it to the worker PEs to keep memory footprint down).
!  3. User code inspects AND retrieves the fetched DATA.

MODULE mo_input_request_list
    USE ISO_C_BINDING, ONLY: C_CHAR, C_SIGNED_CHAR, C_INT, C_DOUBLE, C_NULL_PTR, C_NULL_CHAR, C_ASSOCIATED

    USE mo_cdi, ONLY: t_CdiIterator, cdiIterator_new, cdiIterator_nextField, cdiIterator_delete, cdiIterator_inqVTime, &
                    & cdiIterator_inqLevelType, cdiIterator_inqLevel, cdiIterator_inqGridId, cdiIterator_inqVariableName, &
                    & gridInqType, gridInqUuid, CDI_UNDEFID, ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_HYBRID, &
                    & ZAXIS_HYBRID_HALF, &
                    & ZAXIS_PRESSURE, ZAXIS_HEIGHT, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_ISENTROPIC, &
                    & ZAXIS_TRAJECTORY, ZAXIS_ALTITUDE, ZAXIS_SIGMA, ZAXIS_MEANSEA, ZAXIS_TOA, ZAXIS_SEA_BOTTOM, &
                    & ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP, ZAXIS_ISOTHERM_ZERO, ZAXIS_SNOW, ZAXIS_LAKE_BOTTOM, &
                    & ZAXIS_SEDIMENT_BOTTOM, ZAXIS_SEDIMENT_BOTTOM_TA, ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_MIX_LAYER, &
                    & ZAXIS_REFERENCE, cdiIterator_inqTile, CDI_NOERR, CDI_EINVAL, GRID_UNSTRUCTURED, t_CdiParam, &
                    & cdiIterator_inqParamParts, gridInqNumber, gridInqPosition, cdiGribIterator_inqLongValue, t_CdiGribIterator, &
                    & cdiGribIterator_clone, cdiGribIterator_delete, cdiIterator_inqRTime, CDI_UUID_SIZE, &
                    & cdiIterator_inqFiletype, FILETYPE_GRB, FILETYPE_GRB2, institutInq, institutInqNamePtr
    USE mo_dictionary, ONLY: t_dictionary
    USE mo_exception, ONLY: message, finish
    USE mo_grid_config, ONLY: n_dom
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_initicon_config, ONLY: timeshift, lconsistency_checks
    USE mo_parallel_config, ONLY: use_omp_input
    USE mo_initicon_utils, ONLY: initicon_inverse_post_op
    USE mo_input_container, ONLY: t_InputContainer, inputContainer_make
    USE mo_kind, ONLY: wp, dp
    USE mo_lnd_nwp_config, ONLY: tile_list
    USE mo_nwp_sfc_tiles, ONLY: t_tile_att, t_tileinfo_icon, t_tileinfo_grb2, trivial_tile_att
    USE mo_math_types, ONLY: t_Statistics
    USE mo_model_domain, ONLY: t_patch
    USE mo_mpi, ONLY: my_process_is_mpi_workroot, p_bcast, process_mpi_root_id, p_comm_work, &
                    & p_pe, p_isEqual, p_mpi_wtime
    USE mo_run_config, ONLY: msg_level
    USE mo_time_config, ONLY: time_config
    USE mo_util_string, ONLY: real2string, int2string, &
         toCharArray, toCharacter, c2f_char, &
         charArray_equal, charArray_toLower, tolower, &
                            & charArray_dup, one_of
    USE mo_util_table, ONLY: t_table, initialize_table, add_table_column, set_table_entry, print_table, finalize_table
    USE mo_util_uuid_types, ONLY: t_uuid, uuid_string_length
    USE mo_util_uuid, ONLY: uuid_unparse, OPERATOR(==)
    USE mtime, ONLY: datetime, timedelta, newDatetime, datetimeToString, newTimedelta, timedeltaToString, deallocateDatetime, &
                   & deallocateTimedelta, max_timedelta_str_len, max_datetime_str_len, OPERATOR(-), OPERATOR(+), OPERATOR(==)

    IMPLICIT NONE

PUBLIC :: t_InputRequestList, InputRequestList_create

    TYPE :: t_InputRequestList
        PRIVATE
        TYPE(t_ListEntry), POINTER :: list(:)
        INTEGER :: variableCount

    CONTAINS
        PROCEDURE :: request => InputRequestList_request    !< Require that a variable be read.
        PROCEDURE :: requestMultiple => InputRequestList_requestMultiple    !< Require that a list of variables. Unlike request() this will request the 
                                                                            ! trimmed strings (because it's impossible to pass an array of strings of different LEN).
        PROCEDURE :: readFile => InputRequestList_readFile  !< Scan a file for input data to satisfy the requests.

        PROCEDURE :: getLevels => InputRequestList_getLevels    !< Get the count AND height values (elevation/presure/whatever) of all the levels encountered IN the file.

        !> The `fetchXXX()` methods simply RETURN FALSE IF the DATA could NOT be fetched entirely.
        !> This IS an atomic operation: Either the entire output array IS overwritten OR it IS NOT touched at all.
        PROCEDURE :: fetch2d => InputRequestList_fetch2d
        PROCEDURE :: fetch3d => InputRequestList_fetch3d
        PROCEDURE :: fetchSurface => InputRequestList_fetchSurface  !No level given, fail IF there are several levels.
        PROCEDURE :: fetchTiled2d => InputRequestList_fetchTiled2d
        PROCEDURE :: fetchTiled3d => InputRequestList_fetchTiled3d
        PROCEDURE :: fetchTiledSurface => InputRequestList_fetchTiledSurface  !No level given, fail IF there are several levels.

        !> The `fetchRequiredXXX()` methods will CALL `finish()` IF there are holes IN the DATA that was READ.
        PROCEDURE :: fetchRequired2d => InputRequestList_fetchRequired2d
        PROCEDURE :: fetchRequired3d => InputRequestList_fetchRequired3d
        PROCEDURE :: fetchRequiredSurface => InputRequestList_fetchRequiredSurface
        PROCEDURE :: fetchRequiredTiled2d => InputRequestList_fetchRequiredTiled2d
        PROCEDURE :: fetchRequiredTiled3d => InputRequestList_fetchRequiredTiled3d
        PROCEDURE :: fetchRequiredTiledSurface => InputRequestList_fetchRequiredTiledSurface

        PROCEDURE :: printInventory => InputRequestList_printInventory
        PROCEDURE :: checkRuntypeAndUuids => InputRequestList_checkRuntypeAndUuids

        PROCEDURE :: findIconName => InputRequestList_findIconName    !< Retrieve a t_ListEntry for the given ICON variable name if it exists already.

        PROCEDURE :: destruct => InputRequestList_destruct  !< Destructor.

        PROCEDURE, PRIVATE :: checkRequests => InputRequestList_checkRequests   !< Check that all processes IN the communicator have the same view on which variables are needed.
        PROCEDURE, PRIVATE :: translateNames => InputRequestList_translateNames !< Recalculates the translatedNames of all list entries using the given dictionary.
        PROCEDURE, PRIVATE :: findTranslatedName => InputRequestList_findTranslatedName    !< As findIconName, but uses the translatedVarName.
        PROCEDURE, PRIVATE :: sendStopMessage => InputRequestList_sendStopMessage
        PROCEDURE, PRIVATE :: sendFieldMetadata => InputRequestList_sendFieldMetadata
        PROCEDURE, PRIVATE :: receiveFieldMetadata => InputRequestList_receiveFieldMetadata
        PROCEDURE, PRIVATE :: isRecordValid => InputRequestList_isRecordValid
        PROCEDURE, PRIVATE :: nextField => InputRequestList_nextField
    END TYPE

PRIVATE

    ! These objects are created via findDomainData(, , opt_lcreate = .TRUE.), which will already instanciate an empty container.
    ! On the I/O PE, it IS the job of InputRequestList_isRecordValid() to immediately add a MetadataCache, so that ANY DomainData 
    ! object returned by findDomainData() CONTAINS both a valid InputContainer AND a valid MetadataCache.
    TYPE :: t_DomainData
        INTEGER :: jg
        TYPE(t_DomainData), POINTER :: next

        CLASS(t_InputContainer), POINTER :: container
        TYPE(t_MetadataCache), POINTER :: metadata  !< Some metadata connected with the variable, which IS only used for consistency checking AND printing of the inventory table.
        TYPE(t_Statistics) :: statistics
    END TYPE

    TYPE :: t_ListEntry
      !> The name as it has been requested.
      CHARACTER(len=:), ALLOCATABLE :: iconVarName
      !> The name as it is matched against the stored name. This has
      !! dictionary translation, trimming, and case canonicalization
      !! applied.
      CHARACTER(len=:), ALLOCATABLE :: translatedVarName
        TYPE(t_DomainData), POINTER :: domainData   !< A linked list of an InputContainer AND a MetadataCache for each domain. Only accessed via findDomainData().
    END TYPE

    TYPE :: t_MetadataCache
        CHARACTER(KIND = C_CHAR), POINTER :: rtime(:), vtime(:)
        INTEGER :: levelType, gridNumber, gridPosition, runClass, experimentId, generatingProcessType
        TYPE(t_CdiParam) :: param
        TYPE(t_uuid) :: gridUuid

    CONTAINS
        PROCEDURE :: equalTo => MetadataCache_equalTo
        PROCEDURE :: destruct => MetadataCache_destruct
    END TYPE

    CHARACTER(*), PARAMETER :: modname = "mo_input_request_list"
    LOGICAL, PARAMETER :: debugModule = .FALSE.

CONTAINS

    !Can't use a type constructor interface t_InputRequestList() since the cray compiler looses the list pointer while returning + assigning the function result.
    FUNCTION InputRequestList_create() RESULT(resultVar)
        TYPE(t_InputRequestList), POINTER :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_create"
        INTEGER :: error, i

        ALLOCATE(resultVar, STAT = error)
        if(error /= SUCCESS) CALL finish(routine, "error allocating memory")

        resultVar%variableCount = 0
        ALLOCATE(resultVar%list(8), STAT = error)
        if(error /= SUCCESS) CALL finish(routine, "error allocating memory")
        DO i = 1, SIZE(resultVar%list, 1)
            resultVar%list(i)%domainData => NULL()
        END DO
    END FUNCTION InputRequestList_create

    FUNCTION findDomainData(listEntry, jg, opt_lcreate) RESULT(resultVar)
        TYPE(t_ListEntry), POINTER, INTENT(INOUT) :: listEntry
        INTEGER, INTENT(IN) :: jg
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lcreate
        TYPE(t_DomainData), POINTER :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":findDomainData"
        INTEGER :: error

        IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, "assertion failed, listEntry IS NOT ASSOCIATED")

        ! Try to find a preexisting DomainData object.
        resultVar => listEntry%domainData
        DO WHILE (ASSOCIATED(resultVar))
            IF(resultVar%jg == jg) RETURN
            resultVar => resultVar%next
        END DO

        ! Nothing preexisting found, should we create a new one?
        if(PRESENT(opt_lcreate)) THEN
            IF(opt_lcreate) THEN
                ALLOCATE(resultVar, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "error allocating memory")
                resultVar%jg = jg
                resultVar%next => listEntry%domainData
                resultVar%container => InputContainer_make()
                resultVar%metadata => NULL()
                CALL resultVar%statistics%reset()
                listEntry%domainData => resultVar
            END IF
        END IF
    END FUNCTION findDomainData

    SUBROUTINE InputRequestList_translateNames(me, opt_dict)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_dictionary), OPTIONAL, INTENT(IN) :: opt_dict

        INTEGER :: i

        DO i = 1, me%variableCount
          IF(PRESENT(opt_dict)) THEN
            me%list(i)%translatedVarName &
                 = tolower(TRIM(opt_dict%get(me%list(i)%iconVarName, &
                 &                           me%list(i)%iconVarName)))
          ELSE
            me%list(i)%translatedVarName = tolower(me%list(i)%iconVarName)
          END IF
        END DO
    END SUBROUTINE InputRequestList_translateNames

    FUNCTION InputRequestList_findIconName(me, fieldName, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: fieldName
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug
        TYPE(t_ListEntry), POINTER :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_findIconName"
        INTEGER :: i
        LOGICAL :: debugInfo
        CHARACTER(:), POINTER :: tempName

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        resultVar => NULL()
        DO i = 1, me%variableCount
            tempName => me%list(i)%iconVarName
            IF(fieldName == tempName) THEN
                IF(debugInfo) CALL message(routine, fieldName//" == "//tempName)
                resultVar => me%list(i)
                RETURN
            ELSE
                IF(debugInfo) CALL message(routine, fieldName//" /= "//tempName)
            END IF
        END DO
    END FUNCTION InputRequestList_findIconName

    FUNCTION InputRequestList_findTranslatedName(me, fieldName) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(INOUT) :: fieldName
        TYPE(t_ListEntry), POINTER :: resultVar

        INTEGER :: i

        fieldname = toLower(fieldName)
        resultVar => NULL()
        DO i = 1, me%variableCount
          IF(fieldName == me%list(i)%translatedVarName) THEN
            resultVar => me%list(i)
            RETURN
          END IF
        END DO
    END FUNCTION InputRequestList_findTranslatedName

    !XXX: This also ensures that the requests have been given in the same order on all processes. Not technically necessary, but easier to 
    !     implement and I guess, if the order is not the same, that's a hint that there is a bug somewhere else.
    SUBROUTINE InputRequestList_checkRequests(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_checkRequests"
        INTEGER :: i, j, error, concatenatedSize, curSize, accumulatedSize
        CHARACTER(KIND = C_CHAR), ALLOCATABLE :: concatenatedNames(:)

        !compute the concatenation of all requested variables
        concatenatedSize = 0
        DO i = 1, me%variableCount
          concatenatedSize = concatenatedSize + LEN(me%list(i)%iconVarName) + 1
        END DO
        ALLOCATE(concatenatedNames(concatenatedSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        accumulatedSize = 0
        DO i = 1, me%variableCount
          curSize = LEN(me%list(i)%iconVarName)
            DO j = 1, curSize
                concatenatedNames(accumulatedSize + j) = me%list(i)%iconVarName(j:j)
            END DO
            concatenatedNames(accumulatedSize + curSize + 1) = C_NULL_CHAR
            accumulatedSize = accumulatedSize + curSize + 1
        END DO

        !check that all processes have the same concatenatedNames string
        IF(.NOT. p_isEqual(concatenatedNames, p_comm_work)) THEN
            print*, "process ", p_pe, " has the variable list: ", concatenatedNames
            CALL finish(routine, "not all processes have the same requests in their t_InputRequestList")
        END IF
    END SUBROUTINE InputRequestList_checkRequests

    SUBROUTINE InputRequestList_request(me, fieldName)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: fieldName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_request"
        INTEGER :: i, listSize, error
        TYPE(t_ListEntry), POINTER :: tempList(:), newEntry

        ! don't add a name twice
        IF(ASSOCIATED(me%findIconName(fieldName))) RETURN

        IF(debugModule .AND. my_process_is_mpi_workroot()) print*, 'Adding request for variable "', fieldName, '"'

        ! ensure space for the new container
        listSize = SIZE(me%list, 1)
        IF(me%variableCount == listSize) THEN
            ALLOCATE(tempList(2*listSize), STAT = error)
            if(error /= SUCCESS) CALL finish(routine, "error allocating memory")
            DO i = 1, listSize
                tempList(i) = me%list(i)
            END DO
            DO i = listSize + 1, 2*listSize
                tempList(i)%domainData => NULL()
            END DO
            DEALLOCATE(me%list)
            me%list => tempList
        END IF

        ! add the entry to our list
        me%variableCount = me%variableCount + 1
        newEntry => me%list(me%variableCount)

        newEntry%iconVarName = fieldName
    END SUBROUTINE InputRequestList_request

    SUBROUTINE InputRequestList_requestMultiple(me, fieldNames)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: fieldNames(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_requestMultiple"
        INTEGER :: i

        DO i = 1, SIZE(fieldNames, 1)
            CALL me%request(TRIM(fieldNames(i)))
        END DO
    END SUBROUTINE InputRequestList_requestMultiple

    SUBROUTINE fail(message, variableName, resultVar)
        CHARACTER(LEN = *), INTENT(IN) :: message
        CHARACTER(LEN = *), INTENT(in) :: variableName
        LOGICAL, INTENT(inout) :: resultVar

        IF(msg_level >= 1) print*, 'invalid record for variable "', variableName, '" encountered: '//message
        resultVar = .FALSE.
    END SUBROUTINE fail

    LOGICAL FUNCTION InputRequestList_isRecordValid(me, iterator, p_patch, level, tileId, variableName, lIsFg) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_CdiIterator) :: iterator
        TYPE(t_patch), INTENT(IN) :: p_patch
        REAL(dp), INTENT(OUT) :: level
        INTEGER, INTENT(OUT) :: tileId
        CHARACTER(:), ALLOCATABLE, INTENT(out) :: variableName
        LOGICAL, INTENT(IN) :: lIsFg

        INTEGER(KIND = C_INT) :: error, gridId, gridType, tileIndex, tileAttribute
        REAL(KIND = C_DOUBLE) :: levelValue
        TYPE(t_MetadataCache), POINTER :: metadata
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        TYPE(t_CdiGribIterator) :: gribIterator
        TYPE(datetime), POINTER :: tempTime, iniTime, startTime
        CHARACTER(:), POINTER :: vtimeString
        CHARACTER(max_datetime_str_len) :: debugDatetimeString
        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_isRecordValid"
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute
        TYPE(t_tileinfo_grb2) :: tileinfo_grb2
        TYPE(t_tileinfo_icon) :: tileinfo_icon
        CHARACTER(KIND = C_CHAR), POINTER :: instName(:)  ! institute name
        INTEGER :: instId                                 ! institute ID
        INTEGER :: generatingCenter, generatingSubCenter
        CHARACTER(KIND=C_CHAR), POINTER :: variableName_(:)

        resultVar = .TRUE.
        variableName_ => cdiIterator_inqVariableName(iterator)
        CALL c2f_char(variableName, variableName_)
        DEALLOCATE(variableName_)
        metadata => MetadataCache_create()
        CALL cdiIterator_inqParamParts(iterator, metadata%param%discipline, metadata%param%category, metadata%param%number)

        !Check the time.
        metadata%vtime => cdiIterator_inqVTime(iterator)
        metadata%rtime => cdiIterator_inqRTime(iterator)
        IF (.NOT. ASSOCIATED(metadata%vtime) .OR. .NOT. ASSOCIATED(metadata%rtime)) THEN
          CALL finish(routine, "Internal error!")
        END IF

        IF(lconsistency_checks) THEN
            vtimeString => toCharacter(metadata%vtime)
            tempTime => newDatetime(vtimeString)

            ALLOCATE(iniTime, STAT = error)
            IF(error /= SUCCESS) CALL fail("memory allocation failure", variableName, resultVar)
            iniTime = time_config%tc_startdate
            IF(lIsFg) THEN
                ! add timeshift to INI-datetime to get true starting time
                ALLOCATE(startTime, STAT = error)
                IF(error /= SUCCESS) CALL fail("memory allocation failure", variableName, resultVar)
                startTime = iniTime + timeshift%mtime_shift
                IF(.NOT.(tempTime == startTime)) THEN
                    CALL datetimeToString(startTime, debugDatetimeString)
                    CALL fail("vtime of first-guess field ("//vtimeString//") does not match model start time (" &
                             &//TRIM(debugDatetimeString)//")", variableName, resultVar)
                END IF
                DEALLOCATE(startTime)
            ELSE
                IF(.NOT.(tempTime == iniTime)) THEN
                    CALL datetimeToString(iniTime, debugDatetimeString)
                    CALL fail("vtime of analysis field ("//vtimeString//") does not match model initialization time (" &
                             &//TRIM(debugDatetimeString)//")", variableName, resultVar)
                END IF
            END IF
            DEALLOCATE(iniTime)
            DEALLOCATE(vtimeString)
            CALL deallocateDatetime(tempTime)
        END IF

        !We only check the primary (top) level (selector = 0). Usually that's the only one, but GRIB does allow a secondary lower boundary level.
        metadata%levelType = cdiIterator_inqLevelType(iterator, 0)
        SELECT CASE(metadata%levelType)
            !the level types that translate to a single height VALUE
            CASE(ZAXIS_SURFACE, ZAXIS_PRESSURE, ZAXIS_HEIGHT, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_ALTITUDE, &
                &ZAXIS_REFERENCE, ZAXIS_SNOW)

                levelValue = 0.0
                error = cdiIterator_inqLevel(iterator, 1, outValue1 = levelValue)
                level = REAL(levelValue, dp)
                IF(error /= 0) CALL fail("cdiIterator_inqLevel() failed", variableName, resultVar)
                !TODO: check the zaxis UUID

            !the level types for special levels
            CASE(ZAXIS_TOA, ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP, ZAXIS_ISOTHERM_ZERO, &
                &ZAXIS_MEANSEA, ZAXIS_SEA_BOTTOM, ZAXIS_LAKE_BOTTOM, ZAXIS_SEDIMENT_BOTTOM, ZAXIS_SEDIMENT_BOTTOM_TA, &
                &ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_MIX_LAYER)

                level = REAL(-metadata%levelType, dp)

            !the known z-axis types that are NOT handled by this code
            CASE(ZAXIS_GENERIC)
                CALL fail("z-axis type ZAXIS_GENERIC is not implemented", variableName, resultVar)
            CASE(ZAXIS_HYBRID)
                CALL fail("z-axis type ZAXIS_HYBRID is not implemented", variableName, resultVar)
            CASE(ZAXIS_HYBRID_HALF)
                CALL fail("z-axis type ZAXIS_HYBRID_HALF is not implemented", variableName, resultVar)
            CASE(ZAXIS_ISENTROPIC)
                CALL fail("z-axis type ZAXIS_ISENTROPIC is not implemented", variableName, resultVar)
            CASE(ZAXIS_TRAJECTORY)
                CALL fail("z-axis type ZAXIS_TRAJECTORY is not implemented", variableName, resultVar)
            CASE(ZAXIS_SIGMA)
                CALL fail("z-axis type ZAXIS_SIGMA is not implemented", variableName, resultVar)

            !fallback to catch future expansions of the list of available z-axis types
            CASE DEFAULT
                CALL fail("unknown z-axis TYPE ("//int2string(metadata%levelType)//")", variableName, resultVar)
        END SELECT

        !Check the grid.
        gridId = cdiIterator_inqGridId(iterator)
        IF(gridId /= CDI_UNDEFID) THEN
            !XXX: I believe, it's enough sanity checking if we check the type and the size of the grid,
            !     I don't want to go into checking the lon/lat for all its vertices here...
            !     A test for the correct grid SIZE IS IMPLICIT IN the t_InputContainer when it selects the scatter pattern to USE.
            gridType = gridInqType(gridId)
            IF(gridType /= GRID_UNSTRUCTURED) THEN
                CALL fail("support for this gridtype is not implemented (CDI grid type = "//TRIM(int2string(gridType))//")", variableName, resultVar)
            ELSE
                CALL gridInqUuid(gridId, metadata%gridUuid%DATA)
                metadata%gridNumber = gridInqNumber(gridID)
                metadata%gridPosition = gridInqPosition(gridID)
            END IF
        ELSE
            CALL fail("couldn't inquire grid ID", variableName, resultVar)
        END IF

        error = cdiIterator_inqTile(iterator, tileIndex, tileAttribute);
        SELECT CASE(error)
            CASE(CDI_NOERR)
              tileinfo_grb2%idx = tileIndex
              tileinfo_grb2%att = tileAttribute
              this_att => tile_list%getTileAtt(tileinfo_grb2)
              tileinfo_icon = this_att%getTileinfo_icon()
              tileId = tileinfo_icon%idx

            CASE(CDI_EINVAL)
              !There IS no tile information connected to this field, so we USE the trivial tileId.
              tileinfo_icon = trivial_tile_att%getTileinfo_icon()
              tileId = tileinfo_icon%idx

            CASE DEFAULT
                CALL finish(routine, "unexpected error while reading tile information")
        END SELECT

        !Fetch some additional metadata.
        metadata%runClass = -1
        metadata%experimentId = -1
        metadata%generatingProcessType = -1
        IF(resultVar) THEN
            gribIterator = cdiGribIterator_clone(iterator)
            IF(C_ASSOCIATED(gribIterator%ptr)) THEN
                metadata%runClass = INT(cdiGribIterator_inqLongValue(gribIterator, "backgroundProcess"))
                metadata%generatingProcessType = INT(cdiGribIterator_inqLongValue(gribIterator, "typeOfGeneratingProcess"))
                !
                ! fetching GRIB2 key "localNumberOfExperiment" is restricted to input data generated by DWD
                generatingCenter = INT(cdiGribIterator_inqLongValue(gribIterator, "centre"))
                generatingSubCenter = INT(cdiGribIterator_inqLongValue(gribIterator, "subCentre"))
                instId = institutInq(generatingCenter, generatingSubcenter, '', '')
                instName => institutInqNamePtr(instId) 
                !
                ! Check if instName is associated as CDI returns a NULL pointer if the center is not found within
                ! a CDI internal list (instituteDefaultEntries)
                IF (ASSOCIATED(instName)) THEN
                  IF (TRIM(toCharacter(instName))=="DWD") THEN 
                    metadata%experimentId = INT(cdiGribIterator_inqLongValue(gribIterator, "localNumberOfExperiment"))
                  ENDIF
                ENDIF
                CALL cdiGribIterator_delete(gribIterator)
            END IF
        END IF

        !Check whether the metadata of this record IS consistent with the metadata we've already seen for this variable.
        listEntry => me%findTranslatedName(variableName)
        IF(.NOT.ASSOCIATED(listEntry)) THEN
            resultVar = .FALSE.    !We are NOT interested IN this variable.
            CALL metadata%destruct()
            DEALLOCATE(metadata)
            RETURN
        END IF
        IF(resultVar) domainData => findDomainData(listEntry, p_patch%id)
        IF(resultVar .AND. ASSOCIATED(domainData)) resultVar = metadata%equalTo(domainData%metadata)

        !Commit AND cleanup.
        IF(resultVar) THEN
            IF(ASSOCIATED(domainData)) THEN
                CALL metadata%destruct()
                DEALLOCATE(metadata)
            ELSE
                domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
                domainData%metadata => metadata  !We don't have a metadata cache yet, so we just remember this one.
            END IF
        ELSE
            !The record was not valid.
            DEALLOCATE(variableName)
            CALL metadata%destruct()
            DEALLOCATE(metadata)
        END IF

    END FUNCTION InputRequestList_isRecordValid

    ! Format of the messages which are broadcasted by the following three FUNCTION:
    !
    ! REAL(dp) :: message(3)
    ! message(1) = length of variable NAME, zero IF this is a stop message
    ! message(2) = level
    ! message(3) = tileId
    !
    ! In the case that the variable NAME length is nonzero, this is followed by another message containing the NAME itself.
    ! Note: The broadcasts IN receiveFieldMetadata() are matched with the broadcasts within sendFieldMetadata() and sendStopMessage().
    SUBROUTINE InputRequestList_sendStopMessage(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_sendStopMessage"
        REAL(dp) :: message(3)

        message(1) = 0.0_dp
        message(2) = 0.0_dp
        message(3) = 0.0_dp
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)
    END SUBROUTINE InputRequestList_sendStopMessage

    SUBROUTINE InputRequestList_sendFieldMetadata(me, level, tileId, variableName)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        REAL(dp), INTENT(in) :: level
        INTEGER, INTENT(in) :: tileId
        CHARACTER(*), INTENT(IN) :: variableName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_sendFieldMetadata"
        REAL(dp) :: message(3)
        CHARACTER(:), POINTER :: tempName
        INTEGER :: error

        message(1) = REAL(LEN(variableName), dp)
        message(2) = level
        message(3) = REAL(tileId, dp)
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)

        IF(debugModule) print*, 'Reading field for variable "'//variableName//'"'
        CALL p_bcast(variableName, process_mpi_root_id, p_comm_work)
    END SUBROUTINE InputRequestList_sendFieldMetadata

    LOGICAL FUNCTION InputRequestList_receiveFieldMetadata(me, level, tileId, variableName) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        REAL(dp), INTENT(INOUT) :: level
        INTEGER, INTENT(INOUT) :: tileId
        CHARACTER(:), ALLOCATABLE, INTENT(out) :: variableName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_receiveFieldMetadata"
        REAL(dp) :: message(3)
        CHARACTER(:), ALLOCATABLE :: tempName
        INTEGER :: error

        message(1) = 0.0_dp
        message(2) = 0.0_dp
        message(3) = 0.0_dp
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)
        level = message(2)
        tileId = INT(message(3))
        resultVar = message(1) /= 0.0_dp
        IF (resultVar) THEN
          ALLOCATE(CHARACTER(LEN = INT(message(1))) :: variableName, STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "error allocating memory")
          CALL p_bcast(variableName, process_mpi_root_id, p_comm_work)
        END IF
    END FUNCTION InputRequestList_receiveFieldMetadata

    ! Find the next field that we are interested IN.
    ! This FUNCTION is collective: either all processes RETURN .TRUE. or all RETURN .FALSE. .
    ! When this FUNCTION returns FALSE, THEN there is no further field IN the file.
    !
    ! ignoredRecords IS NOT reset by this FUNCTION, it IS ONLY incremented
    LOGICAL FUNCTION InputRequestList_nextField(me, iterator, p_patch, level, tileId, variableName, ignoredRecords, lIsFg) &
    &RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_CdiIterator), INTENT(IN) :: iterator
        TYPE(t_patch), INTENT(IN) :: p_patch
        REAL(dp), INTENT(OUT) :: level
        INTEGER, INTENT(OUT) :: tileId
        CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: variableName
        INTEGER, INTENT(INOUT) :: ignoredRecords
        LOGICAL, INTENT(IN) :: lIsFg

        resultVar = .FALSE.
        IF(my_process_is_mpi_workroot()) THEN
            ! Scan the file until we find a field that we are interested in.
            DO
                IF(cdiIterator_nextField(iterator) /= 0) THEN
                  IF (.NOT. use_omp_input) THEN
                    CALL me%sendStopMessage()
                  ENDIF
                    RETURN
                ELSE
                    IF(me%isRecordValid(iterator, p_patch, level, tileId, variableName, lIsFg)) THEN
                        IF(ASSOCIATED(me%findTranslatedName(variableName))) THEN
                          IF (.NOT. use_omp_input) THEN
                            !NEC: skip communcation here in VH_OMP case, but do in read routine
                            CALL me%sendFieldMetadata(level, tileId, variableName)
                          ENDIF
                          resultVar = .TRUE.
                          RETURN
                        END IF
                    ELSE
                        ignoredRecords = ignoredRecords + 1
                    END IF
                END IF
            END DO
        ELSE
            resultVar = me%receiveFieldMetadata(level, tileId, variableName)
        END IF
    END FUNCTION InputRequestList_nextField

    SUBROUTINE InputRequestList_readFile(me, p_patch, path, lIsFg, opt_dict)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: p_patch
        CHARACTER(LEN = *, KIND = C_CHAR), INTENT(IN) :: path
        TYPE(t_dictionary), OPTIONAL, INTENT(IN) :: opt_dict
        LOGICAL, INTENT(IN) :: lIsFg

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_readFile"
        TYPE(t_CdiIterator) :: iterator
        REAL(dp) :: level
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER :: vtime
        CHARACTER(len = :), ALLOCATABLE :: variableName, variableName_prev
        INTEGER :: i, tileId, recordsRead, recordsIgnored
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp) :: timer(5), savetime
        LOGICAL  :: ret, l_exist

        INTEGER :: iread
        iread = 0

        recordsRead = 0
        recordsIgnored = 0
        timer(1) = p_mpi_wtime()
        timer(2:5) = 0._dp

        CALL me%checkRequests() !sanity checks
        CALL me%translateNames(opt_dict)

        iterator%ptr = C_NULL_PTR
        IF(my_process_is_mpi_workroot()) THEN

            INQUIRE (FILE=path, EXIST=l_exist)
            IF (.NOT.l_exist) THEN
              CALL finish(TRIM(routine),'File is not found: '//TRIM(path))
            ENDIF

            iterator = cdiIterator_new(path)
            IF(.NOT. C_ASSOCIATED(iterator%ptr)) THEN
              CALL finish(routine, "can't open file "//'"'//path//'" for reading')
            END IF
            ! Check whether a Map file for translating fileInputName<=>internalName is required
            ! - a Map File is mandatory, if input is read in GRIB2-Format
            IF (((cdiIterator_inqFiletype(iterator) == FILETYPE_GRB)  .OR.   &
              &  (cdiIterator_inqFiletype(iterator) == FILETYPE_GRB2)) .AND. &
              &  .NOT. PRESENT(opt_dict)) THEN
              CALL finish( routine,                         &
                &  'dictionary missing. It is required when trying to read data in GRIB format.')
            END IF
        END IF
        IF (use_omp_input .AND. my_process_is_mpi_workroot()) THEN
           !NEC_RP: if masterprocess: use readField_omp routine that OMP parallelizes read, statistics and distribution
           DO 
               savetime = p_mpi_wtime()
               ret = me%nextField(iterator, p_patch, level, tileId, variableName, recordsIgnored, lIsFg)
               timer(2) = timer(2) + p_mpi_wtime() - savetime
               IF (.NOT. ret) EXIT
               IF (ALLOCATED(variableName_prev)) THEN
               IF (variableName /= variableName_prev) THEN
                  CALL domainData%container%readField_omp(variableName_prev, level, tileId, timer, &
                       p_patch%id, iterator, domainData%statistics, -1)
                  iread = 0
               END IF
               END IF
               recordsRead = recordsRead + 1
               ! We have now found the next field that we are interested in.
               listEntry => me%findTranslatedName(variableName)
               IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, &
                 "Assertion failed: Processes have different input request lists!")
               domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
               iread = iread + 1
               CALL domainData%container%readField_omp(variableName, level, tileId, timer, &
                  p_patch%id, iterator, domainData%statistics, iread)
               CALL MOVE_ALLOC(variableName, variableName_prev)
           END DO
           IF (ALLOCATED(variableName_prev)) THEN
              CALL domainData%container%readField_omp(variableName_prev, level, tileId, timer, &
                p_patch%id, iterator, domainData%statistics, -2)
              DEALLOCATE(variableName_prev)
           END IF
           CALL me%sendStopMessage()
        ELSE
          !NEC_RP: all other processes use original code
          DO 
            savetime = p_mpi_wtime()
            ret = me%nextField(iterator, p_patch, level, tileId, variableName, recordsIgnored, lIsFg)
            timer(2) = timer(2) + p_mpi_wtime() - savetime
            IF (.NOT. ret) EXIT
            recordsRead = recordsRead + 1
            ! We have now found the next field that we are interested IN.
            listEntry => me%findTranslatedName(variableName)
            IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, "Assertion failed: Processes have different input request lists!")
            domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
            CALL domainData%container%readField(variableName, level, tileId, timer, p_patch%id, iterator, domainData%statistics)
          END DO
        END IF

        timer(1) = p_mpi_wtime() - timer(1)
        IF(my_process_is_mpi_workroot()) THEN
            IF(msg_level > 4) THEN
              WRITE(0, *) routine//": READ "//TRIM(int2string(recordsRead))//" records from file '"//path//"', &
                         &ignoring "//TRIM(int2string(recordsIgnored))//" records"
              WRITE(0, '(3(a,f10.5),a)') ' Timer report: Total ', timer(1), ' s, Read metadata ', timer(2), &
                                       & ' s, Read data ', timer(3), ' s'
              WRITE(0, '(2(a,f10.5),a)') '               Compute statistics ', timer(4), ' s, Distribute data ', timer(5), 's'
            ENDIF
            CALL cdiIterator_delete(iterator)
        END IF
    END SUBROUTINE InputRequestList_readFile

    FUNCTION InputRequestList_getLevels(me, varName, jg, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug
        REAL(dp), POINTER :: resultVar(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_getLevels"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch level data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar => NULL()
        IF(ASSOCIATED(domainData)) resultVar => domainData%container%getLevels()
    END FUNCTION InputRequestList_getLevels

    LOGICAL FUNCTION InputRequestList_fetch2d(me, varName, level, tile, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(IN) :: level
        INTEGER, INTENT(IN) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetch2d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER   :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetch2d(level, tile, outData, opt_lDebug)
        IF(resultVar) THEN
          this_att => tile_list%getTileAtt(t_tileinfo_icon(tile))
          CALL initicon_inverse_post_op( &
            &   TRIM(varName//TRIM(this_att%getTileSuffix())), &
            &   outData)
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetch2d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetch2d

    LOGICAL FUNCTION InputRequestList_fetch3d(me, varName, tile, jg, outData, optLevelDimension, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetch3d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        TYPE(t_tile_att), POINTER   :: this_att  ! pointer to attribute
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetch3d(tile, outData, optLevelDimension, opt_lDebug)
        IF(resultVar .AND. varName /= 'smi' .AND. varName /= 'SMI') THEN   !SMI IS NOT IN the ICON variable lists, so we need to skip inverse postprocessing for it manually.
          this_att => tile_list%getTileAtt(t_tileinfo_icon(tile))
          CALL initicon_inverse_post_op( &
            &   TRIM(varName//TRIM(this_att%getTileSuffix())), &
            &   outData)
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetch3d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetch3d

    LOGICAL FUNCTION InputRequestList_fetchSurface(me, varName, tile, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchSurface"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp), POINTER :: levels(:)
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) THEN
            levels => domainData%container%getLevels()
            SELECT CASE(SIZE(levels, 1))
                CASE(0)
                    resultVar = .FALSE.
                    IF(debugInfo) CALL message(routine, "no levels found")
                CASE(1)
                    resultVar = domainData%container%fetch2d(levels(1), tile, outData, opt_lDebug)
                    IF(debugInfo .AND. .NOT. resultVar) CALL message(routine, "InputContainer_fetch2d() returned an error")
                CASE DEFAULT
                    CALL finish(routine, "trying to read '"//varName//"' as a surface variable, but the file contains several &
                                         &levels of this variable")
            END SELECT
        END IF
        IF(resultVar) THEN
            this_att => tile_list%getTileAtt(t_tileinfo_icon(tile))
            CALL initicon_inverse_post_op( &
            &   TRIM(varName//TRIM(this_att%getTileSuffix())), &
            &   outData)
        END IF
    END FUNCTION InputRequestList_fetchSurface

    LOGICAL FUNCTION InputRequestList_fetchTiled2d(me, varName, level, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(IN) :: level
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiled2d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        INTEGER :: i
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetchTiled2d(level, outData, opt_lDebug)
        IF(resultVar) THEN
            DO i = 1, SIZE(outData, 3)
                this_att => tile_list%getTileAtt(t_tileinfo_icon(i))
                CALL initicon_inverse_post_op(TRIM(varName//TRIM(this_att%getTileSuffix())), &
                  &                           outData(:,:,i))
            END DO
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetchTiled2d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetchTiled2d

    LOGICAL FUNCTION InputRequestList_fetchTiled3d(me, varName, jg, outData, optLevelDimension, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiled3d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        INTEGER :: i
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetchTiled3d(outData, optLevelDimension, opt_lDebug)
        IF(resultVar .AND. varName /= 'smi' .AND. varName /= 'SMI') THEN   !SMI IS NOT IN the ICON variable lists, so we need to skip inverse postprocessing for it manually.
            DO i = 1, SIZE(outData, 4)
                this_att => tile_list%getTileAtt(t_tileinfo_icon(i))
                CALL initicon_inverse_post_op(TRIM(varName//TRIM(this_att%getTileSuffix())), &
                  &                           outData(:,:,:,i))
            END DO
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetchTiled3d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetchTiled3d

    LOGICAL FUNCTION InputRequestList_fetchTiledSurface(me, varName, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiledSurface"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp), POINTER :: levels(:)
        INTEGER :: i
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) THEN
            levels => domainData%container%getLevels()
            SELECT CASE(SIZE(levels, 1))
                CASE(0)
                    resultVar = .FALSE.
                    IF(debugInfo) CALL message(routine, "no levels found")
                CASE(1)
                    resultVar = domainData%container%fetchTiled2d(levels(1), outData, opt_lDebug)
                    IF(debugInfo .AND. .NOT. resultVar) CALL message(routine, "InputContainer_fetch2d() returned an error")
                CASE DEFAULT
                    CALL finish(routine, "trying to read '"//varName//"' as a surface variable, but the file contains several &
                                         &levels of this variable")
            END SELECT
        END IF
        IF(resultVar) THEN
            DO i = 1, SIZE(outData, 3)
                this_att => tile_list%getTileAtt(t_tileinfo_icon(i))
                CALL initicon_inverse_post_op(TRIM(varName//TRIM(this_att%getTileSuffix())), &
                  &                           outData(:,:,i))
            END DO
        END IF
    END FUNCTION InputRequestList_fetchTiledSurface

    SUBROUTINE InputRequestList_fetchRequired2d(me, varName, level, tile, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(in) :: level
        INTEGER, INTENT(in) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequired2d"

        IF(.NOT. me%fetch2d(varName, level, tile, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequired2d

    SUBROUTINE InputRequestList_fetchRequired3d(me, varName, tile, jg, outData, optLevelDimension)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(in) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequired3d"

        IF(.NOT. me%fetch3d(varName, tile, jg, outData, optLevelDimension)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequired3d

    SUBROUTINE InputRequestList_fetchRequiredSurface(me, varName, tile, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(in) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredSurface"

        IF(.NOT. me%fetchSurface(varName, tile, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredSurface

    SUBROUTINE InputRequestList_fetchRequiredTiled2d(me, varName, level, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(in) :: level
        INTEGER, INTENT(in) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiled2d"

        IF(.NOT. me%fetchTiled2d(varName, level, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiled2d

    SUBROUTINE InputRequestList_fetchRequiredTiled3d(me, varName, jg, outData, optLevelDimension)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiled3d"

        IF(.NOT. me%fetchTiled3d(varName, jg, outData, optLevelDimension)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiled3d

    SUBROUTINE InputRequestList_fetchRequiredTiledSurface(me, varName, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiledSurface"

        IF(.NOT. me%fetchTiledSurface(varName, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiledSurface

    SUBROUTINE InputRequestList_checkRuntypeAndUuids(me, incrementVariables, gridUuids, lIsFg, lHardCheckUuids)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: incrementVariables(:)
        TYPE(t_uuid), INTENT(IN) :: gridUuids(:)    !< gridUuids(n_dom)
        LOGICAL, INTENT(IN) :: lIsFg, lHardCheckUuids

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_checkRuntypeAndUuids"
        INTEGER :: i, jg, expectedRuntype
        TYPE(t_ListEntry), POINTER :: curVar
        TYPE(t_DomainData), POINTER :: curDomain
        CHARACTER(:), POINTER :: varnameString
        CHARACTER(LEN = uuid_string_length) :: expectedUuid, foundUuid

        IF(.NOT.my_process_is_mpi_workroot()) CALL finish(routine, "assertion failed")
        IF(SIZE(gridUuids, 1) /= n_dom) CALL finish(routine, "assertion failed")
        DO jg = 1, n_dom
            DO i = 1, me%variableCount
                curVar => me%list(i)
                curDomain => findDomainData(curVar, jg)
                IF(.NOT.ASSOCIATED(curDomain)) CYCLE
                varnameString => curVar%iconVarName

                ! first check the TYPE of the generating process of the DATA
                IF(.NOT.lIsFg) THEN
                    IF(one_of(varnameString, incrementVariables) > 0) THEN
                        expectedRuntype = 201   ! analysis increment variables
                    ELSE
                        expectedRuntype = 0 ! analysis full variables
                    END IF
                    IF(expectedRuntype /= curDomain%metadata%generatingProcessType) THEN
                        CALL finish(routine, "detected wrong type of generating process on variable '"//varnameString//"', &
                                             &expected "//TRIM(int2string(expectedRuntype))//", &
                                             &found "//TRIM(int2string(curDomain%metadata%generatingProcessType)))
                    END IF
                END IF

                ! second check the UUID of the grid
                IF(.NOT.(gridUuids(jg) == curDomain%metadata%gridUuid)) THEN
                    CALL uuid_unparse(gridUuids(jg), expectedUuid)
                    CALL uuid_unparse(curDomain%metadata%gridUuid, foundUuid)
                    IF(lHardCheckUuids) THEN
                        CALL finish(routine, "detected wrong UUID of grid for variable '"//varnameString//"', &
                                             &expected "//expectedUuid//", &
                                             &found "//foundUuid)
                    ELSE
                        CALL message(routine, "warning: unexpected UUID of grid for variable '"//varnameString//"', &
                                             &expected "//expectedUuid//", &
                                             &found "//foundUuid)
                    END IF
                END IF
            END DO
        END DO
    END SUBROUTINE InputRequestList_checkRuntypeAndUuids

    SUBROUTINE InputRequestList_printInventory(me)
        CLASS(t_InputRequestList), INTENT(IN) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_printInventory"
        INTEGER :: i, jg, curRow, levelCount, tileCount
        LOGICAL :: lUntiledData
        TYPE(t_ListEntry), POINTER :: curVar
        TYPE(t_DomainData), POINTER :: curDomain
        CHARACTER(:), POINTER :: rtimeString, vtimeString
        TYPE(t_table) :: table
        CHARACTER(*), PARAMETER :: domainCol = "jg", &
                                 & variableCol = "variable", &
                                 & tripleCol = "triple", &
                                 & vtimeCol = "validity time", &
                                 & levelTypeCol = "levTyp", &
                                 & levelCountCol = "nlev", &
                                 & tileCountCol = "tileCnt", &
                                 & untiledCol = "untiled", &
                                 & runtypeCol = "runtype", &
                                 & vvmmCol = "vvmm", &
                                 & clasCol = "clas", &
                                 & expidCol = "expid", &
                                 & gridCol = "grid", &
                                 & rgridCol = "rgrid", &
                                 & minCol = "min", &
                                 & meanCol = "mean", &
                                 & maxCol = "max"
        CHARACTER(LEN = 3*3+2) :: parameterString
        TYPE(datetime), POINTER :: rtime, vtime
        TYPE(timedelta), POINTER :: forecastTime
        CHARACTER(len=max_timedelta_str_len) :: forecastTimeString

        CALL initialize_table(table)

        CALL add_table_column(table, domainCol)
        CALL add_table_column(table, variableCol)
        CALL add_table_column(table, tripleCol)
        CALL add_table_column(table, vtimeCol)
        CALL add_table_column(table, vvmmCol)
        CALL add_table_column(table, levelTypeCol)
        CALL add_table_column(table, levelCountCol)
        CALL add_table_column(table, tileCountCol)
        CALL add_table_column(table, untiledCol)
        CALL add_table_column(table, runtypeCol)
        CALL add_table_column(table, clasCol)
        CALL add_table_column(table, expidCol)
        CALL add_table_column(table, gridCol)
        CALL add_table_column(table, rgridCol)
        CALL add_table_column(table, minCol)
        CALL add_table_column(table, meanCol)
        CALL add_table_column(table, maxCol)

        IF(.NOT.my_process_is_mpi_workroot()) CALL finish(routine, "assertion failed")
        curRow = 1  !we can have zero to n_dom rows for each variable, so we can't USE the loop counter for the rows
        DO jg = 1, n_dom
            DO i = 1, me%variableCount
                curVar => me%list(i)
                curDomain => findDomainData(curVar, jg)
                IF(.NOT.ASSOCIATED(curDomain)) CYCLE
                CALL curDomain%container%getCounts(levelCount, tileCount, lUntiledData)

                !domain, NAME, AND triple columns
                CALL set_table_entry(table, curRow, domainCol, TRIM(int2string(curDomain%jg)))
                CALL set_table_entry(table, curRow, variableCol, curVar%iconVarName)
                WRITE(parameterString, '(3(I3,:,"."))') curDomain%metadata%param%discipline, curDomain%metadata%param%category, &
                &                                    curDomain%metadata%param%number
                CALL set_table_entry(table, curRow, tripleCol, parameterString)


                !date AND forecast time columns
                rtimeString => toCharacter(curDomain%metadata%rtime)
                vtimeString => toCharacter(curDomain%metadata%vtime)
                CALL set_table_entry(table, curRow, vtimeCol, vtimeString)

                rtime => newDatetime(rtimeString)
                vtime => newDatetime(vtimeString)

                forecastTime => newTimedelta("PT00H")  ! this 'initialization' IS necessary, IN order to correctly deal with timedelta=0.
                forecastTime = vtime - rtime
                CALL timedeltaToString(forecastTime, forecastTimeString)
                CALL set_table_entry(table, curRow, vvmmCol, TRIM(forecastTimeString))

                CALL deallocateDatetime(rtime)
                CALL deallocateDatetime(vtime)
                CALL deallocateTimedelta(forecastTime)
                DEALLOCATE(rtimeString)
                DEALLOCATE(vtimeString)


                !the simpler columns
                CALL set_table_entry(table, curRow, levelTypeCol, TRIM(int2string(curDomain%metadata%levelType)))
                CALL set_table_entry(table, curRow, levelCountCol, TRIM(int2string(levelCount)))
                IF(tileCount /= 0) CALL set_table_entry(table, curRow, tileCountCol, TRIM(int2string(tileCount)))
                IF(lUntiledData) THEN
                    CALL set_table_entry(table, curRow, untiledCol, "yes")
                ELSE
                    CALL set_table_entry(table, curRow, untiledCol, "no")
                END IF
                CALL set_table_entry(table, curRow, clasCol, TRIM(int2string(curDomain%metadata%runClass)))
                CALL set_table_entry(table, curRow, expidCol, TRIM(int2string(curDomain%metadata%experimentId)))
                IF(curDomain%metadata%generatingProcessType /= -1) THEN
                    CALL set_table_entry(table, curRow, runtypeCol, TRIM(int2string(curDomain%metadata%generatingProcessType)))
                END IF
                IF(curDomain%metadata%gridNumber /= -1) THEN
                    CALL set_table_entry(table, curRow, gridCol, TRIM(int2string(curDomain%metadata%gridNumber)))
                END IF
                IF(curDomain%metadata%gridPosition /= -1) THEN
                    CALL set_table_entry(table, curRow, rgridCol, TRIM(int2string(curDomain%metadata%gridPosition)))
                END IF
                CALL set_table_entry(table, curRow, minCol, TRIM(real2string(curDomain%statistics%MIN)))
                CALL set_table_entry(table, curRow, meanCol, TRIM(real2string(curDomain%statistics%mean)))
                CALL set_table_entry(table, curRow, MAXCol, TRIM(real2string(curDomain%statistics%MAX)))


                !next row
                curRow = curRow + 1
            END DO
        END DO

        CALL print_table(table, opt_delimiter = " | ")
        CALL finalize_table(table)

    END SUBROUTINE InputRequestList_printInventory


    SUBROUTINE InputRequestList_destruct(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        INTEGER :: i
        TYPE(t_ListEntry), POINTER :: currentEntry
        TYPE(t_DomainData), POINTER :: domainData, domainDataTemp

        DO i = 1, me%variableCount
            currentEntry => me%list(i)
            domainData => currentEntry%domainData
            DO WHILE (ASSOCIATED(domainData))
                IF(ASSOCIATED(domainData%container)) THEN
                    CALL domainData%container%destruct()
                    DEALLOCATE(domainData%container)
                END IF
                IF(ASSOCIATED(domainData%metadata)) THEN
                    CALL domainData%metadata%destruct()
                    DEALLOCATE(domainData%metadata)
                END IF
                domainDataTemp => domainData%next
                DEALLOCATE(domainData)
                domainData => domainDataTemp
            END DO
        END DO
        DEALLOCATE(me%list)
    END SUBROUTINE InputRequestList_destruct

    FUNCTION MetadataCache_create() RESULT(resultVar)
        TYPE(t_MetadataCache), POINTER :: resultVar

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":MetadataCache_create"
        INTEGER :: error

        ALLOCATE(resultVar, STAT = error)
        IF(error /= success) CALL finish(routine, "memory allocation error")
        resultVar%vtime => NULL()
        resultVar%rtime => NULL()
    END FUNCTION MetadataCache_create

    LOGICAL FUNCTION MetadataCache_equalTo(me, other) RESULT(resultVar)
        CLASS(t_MetadataCache), INTENT(IN) :: me, other

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":MetadataCache_create"

        INTEGER :: gridNumber, gridPosition, runClass, experimentId, generatingProcessType

        resultVar = .FALSE.

        !compare the time strings
        IF(.NOT.ASSOCIATED(me%rtime).OR..NOT.ASSOCIATED(other%rtime)) THEN
            CALL finish(routine, "internal error, please report this bug")
        END IF
        IF(SIZE(me%rtime) /= SIZE(other%rtime) .OR. ANY(me%rtime /= other%rtime)) THEN
            CALL message(routine, "inconsistent rtime detected")
            RETURN
        END IF

        IF(.NOT.ASSOCIATED(me%vtime).OR..NOT.ASSOCIATED(other%vtime)) THEN
            CALL finish(routine, "internal error, please report this bug")
        END IF
        IF(SIZE(me%vtime) /= SIZE(other%vtime) .OR. ANY(me%vtime /= other%vtime)) THEN
            CALL message(routine, "inconsistent vtime detected")
            RETURN
        END IF

        !compare the parameters
        IF(me%param%discipline /= other%param%discipline) THEN
            CALL message(routine, "inconsistent discipline detected")
            RETURN
        END IF
        IF(me%param%category /= other%param%category) THEN
            CALL message(routine, "inconsistent category detected")
            RETURN
        END IF
        IF(me%param%number /= other%param%number) THEN
            CALL message(routine, "inconsistent number detected")
            RETURN
        END IF

        !compare the other fields
        IF(me%levelType /= other%levelType) THEN
            CALL message(routine, "inconsistent level type detected")
            RETURN
        END IF

        IF(me%gridNumber /= other%gridNumber) THEN
            CALL message(routine, "inconsistent number of grids detected")
            RETURN
        END IF

        IF(me%gridPosition /= other%gridPosition) THEN
            CALL message(routine, "inconsistent grid index detected")
            RETURN
        END IF

        IF(me%runClass /= other%runClass) THEN
            CALL message(routine, "inconsistent run CLASS detected")
            RETURN
        END IF

        IF(me%experimentId /= other%experimentId) THEN
            CALL message(routine, "inconsistent experiment ID detected")
            RETURN
        END IF

        IF(me%generatingProcessType /= other%generatingProcessType) THEN
            CALL message(routine, "inconsistent type of generating process detected")
            RETURN
        END IF

        IF(.NOT.(me%gridUuid == other%gridUuid)) THEN
            CALL message(routine, "inconsistent UUID of grid detected")
            RETURN
        END IF

        resultVar = .TRUE.
    END FUNCTION MetadataCache_equalTo

    SUBROUTINE MetadataCache_destruct(me)
        CLASS(t_MetadataCache), INTENT(INOUT) :: me

        IF(ASSOCIATED(me%vtime)) DEALLOCATE(me%vtime)
        IF(ASSOCIATED(me%rtime)) DEALLOCATE(me%rtime)
    END SUBROUTINE MetadataCache_destruct

END MODULE mo_input_request_list
