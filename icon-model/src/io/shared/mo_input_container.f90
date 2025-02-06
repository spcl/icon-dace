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

MODULE mo_input_container
    USE ISO_C_BINDING, ONLY: C_INT32_T, C_INT64_T, C_DOUBLE, C_FLOAT, C_ASSOCIATED
    USE mo_cdi, ONLY: t_CdiIterator, gridInqSize, cdiIterator_inqGridId, cdiIterator_readField, cdiIterator_readFieldF, &
                    & cdiIterator_inqDatatype, CDI_DATATYPE_PACK23, CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32
    USE mo_communication, ONLY: t_ScatterPattern
    USE mo_exception, ONLY: message, finish, message_text
    USE mo_hash_table, ONLY: t_HashTable, hashTable_make
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: wp, dp
    USE mo_math_types, ONLY: t_Statistics
    USE mo_mpi, ONLY: p_bcast, p_comm_work, p_mpi_wtime, process_mpi_root_id
    USE mo_parallel_config, ONLY: blk_no, nproma
    USE mo_scatter_pattern_base, ONLY: lookupScatterPattern
    USE mo_nwp_sfc_tiles, ONLY: trivial_tile_att, t_tileinfo_icon
    USE mo_util_string, ONLY: int2string, real2string
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

PUBLIC :: t_InputContainer, InputContainer_make

    TYPE :: t_ValueList
        PRIVATE
        REAL(dp), ALLOCATABLE :: values(:)
        INTEGER :: valueCount
    CONTAINS
        PROCEDURE :: init => ValueList_init
        PROCEDURE :: addValue => ValueList_addValue
        PROCEDURE :: haveValue => ValueList_haveValue
        PROCEDURE :: ensureSpace => ValueList_ensureSpace
    END TYPE

    ! A t_InputContainer is used for file driven I/O to cache input DATA IN a distributed fashion
    ! until we have enough information to ALLOCATE the arrays to hold this DATA.
    ! This is especially needed for the fields READ by request_dwdfg_atm_ii() AND fetch_dwdfg_atm_ii()
    ! for which we have no prior knowledge about the count of levels PRESENT IN the input file.
    !
    ! t_InputContainers are always created via a t_InputRequestList which distributes the DATA to the different t_InputContainer objects that it created.
    !
    ! TODO (optional): make it possible to use a different communicator than p_comm_work.
    TYPE :: t_InputContainer
        PRIVATE
        !the NAME of the variable is NOT a part of the t_InputContainer since it's ONLY needed by the t_InputRequestList
        TYPE(t_HashTable), POINTER :: fields    !this is a collection of all the different 2D levels that have been READ (tile, level)
        INTEGER :: fieldCount
        TYPE(t_ValueList) :: tiles, levels

    CONTAINS
        PROCEDURE :: getLevels => InputContainer_getLevels
        PROCEDURE :: getCounts => InputContainer_getCounts
        PROCEDURE :: destruct => InputContainer_destruct

        !the fetch operations are atomical, the output array IS NOT touched upon failure
        PROCEDURE :: fetch2D => InputContainer_fetch2D  !returns .TRUE. on success
        PROCEDURE :: fetch3D => InputContainer_fetch3D  !returns .TRUE. on success
        PROCEDURE :: fetchTiled2D => InputContainer_fetchTiled2D  !returns .TRUE. on success
        PROCEDURE :: fetchTiled3D => InputContainer_fetchTiled3D  !returns .TRUE. on success

        PROCEDURE :: readField => InputContainer_readField
        !NEC: new routine only for this one communication loop
        PROCEDURE :: readField_omp => InputContainer_readField_omp
        PROCEDURE, PRIVATE :: dataAvailable => InputContainer_dataAvailable
    END TYPE

PRIVATE

    TYPE :: t_LevelKey
      PRIVATE
      REAL(dp) :: levelValue
      INTEGER :: tileId   !< the ICON-internal tile id
    END TYPE

    TYPE :: t_LevelPointer
      PRIVATE
      REAL(wp), ALLOCATABLE :: ptr(:,:)   !dimensions: (nproma, nblks_*)
    END TYPE

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_input_container"
    REAL(dp), PARAMETER :: k_levelComparisonInterval = 1e-5

CONTAINS

    INTEGER(C_INT32_T) FUNCTION InputContainer_hashKey(key) RESULT(resultVar)
        CLASS(*), POINTER, INTENT(IN) :: key

        ! Just some large primes to produce some good pseudorandom bits in the hashes.
        INTEGER(C_INT64_T), PARAMETER :: prime1 = 729326603, prime2 = 941095657
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":InputContainer_hashKey"

        INTEGER(C_INT64_T) :: temp
        INTEGER(KIND = C_INT64_T), PARAMETER :: mask = 2_C_INT64_T**31 - 1_C_INT64_T    ! a bitmask for the 31 low order bits

        SELECT TYPE(key)
        TYPE IS(t_LevelKey)
          temp = prime1*INT(key%levelValue, C_INT64_T)
          temp = temp + prime2*INT(key%tileId, C_INT64_T)
          resultVar = INT(IAND(temp, mask), C_INT32_T)
        CLASS DEFAULT
          CALL finish(routine, "illegal argument type")
        END SELECT
    END FUNCTION InputContainer_hashKey

    LOGICAL FUNCTION InputContainer_equalKeysFunction(keyA, keyB) RESULT(resultVar)
        CLASS(*), POINTER, INTENT(IN) :: keyA, keyB

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":InputContainer_equalKeysFunction"

        SELECT TYPE(keyA)
        TYPE IS(t_LevelKey)
          SELECT TYPE(keyB)
          TYPE IS(t_LevelKey)
            resultVar = keyA%levelValue > keyB%levelValue - k_levelComparisonInterval .AND. &
                 & keyA%levelValue < keyB%levelValue + k_levelComparisonInterval
            resultVar = resultVar .AND. keyA%tileId == keyB%tileId
          CLASS DEFAULT
            CALL finish(routine, "illegal argument type")
          END SELECT
        CLASS DEFAULT
          CALL finish(routine, "illegal argument type")
        END SELECT
    END FUNCTION InputContainer_equalKeysFunction

    FUNCTION InputContainer_make() RESULT(resultVar)
        CLASS(t_InputContainer), POINTER :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_make"
        INTEGER :: error

        ALLOCATE(resultVar, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        resultVar%fields => hashTable_make(InputContainer_hashKey, InputContainer_equalKeysFunction)
        resultVar%fieldCount = 0
        CALL resultVar%tiles%init()
        CALL resultVar%levels%init()
    END FUNCTION InputContainer_make

    FUNCTION InputContainer_getLevels(me) RESULT(resultVar)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), POINTER :: resultVar(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_getLevels"
        INTEGER :: i, error

        ALLOCATE(resultVar(me%levels%valueCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        DO i = 1, me%levels%valueCount
            resultVar(i) = me%levels%values(i)
        END DO
    END FUNCTION InputContainer_getLevels

    SUBROUTINE InputContainer_getCounts(me, levelCount, tileCount, lHaveUntiledData)
        CLASS(t_InputContainer), INTENT(IN) :: me
        INTEGER, INTENT(OUT) :: levelCount, tileCount
        LOGICAL, INTENT(OUT) :: lHaveUntiledData

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_getCounts"
        TYPE(t_tileinfo_icon) :: trivial_tileinfo_icon

        levelCount = me%levels%valueCount
        tileCount = me%tiles%valueCount
        trivial_tileinfo_icon = trivial_tile_att%getTileinfo_icon()
        lHaveUntiledData = me%tiles%haveValue(REAL(trivial_tileinfo_icon%idx, dp))
        IF(lHaveUntiledData) tileCount = tileCount - 1
    END SUBROUTINE InputContainer_getCounts

    LOGICAL FUNCTION InputContainer_dataAvailable(me, opt_level, opt_tile, opt_jg) RESULT(resultVar)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), OPTIONAL, INTENT(IN) :: opt_level
        INTEGER, OPTIONAL, INTENT(IN) :: opt_tile, opt_jg

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_dataAvailable"
        TYPE(t_levelkey), TARGET :: key
        CLASS(*), POINTER :: polymorphicKey, val
        INTEGER :: levelCount, tileCount, levelIndex, tileIndex, error

        resultVar = .FALSE.
        IF(me%levels%valueCount <= 0 .OR. me%tiles%valueCount <= 0) RETURN

        levelCount = me%levels%valueCount
        IF(PRESENT(opt_level)) levelCount = 1
        tileCount = me%tiles%valueCount
        IF(PRESENT(opt_tile)) tileCount = 1

        polymorphicKey => key
        resultVar = .TRUE.
        outerLoop: DO levelIndex = 1, levelCount
          key%levelValue = me%levels%values(levelIndex)
          IF(PRESENT(opt_level)) key%levelValue = opt_level
          DO tileIndex = 1, tileCount
            key%tileId = INT(me%tiles%values(tileIndex))
            IF(PRESENT(opt_tile)) key%tileId = opt_tile
            val => me%fields%getEntry(polymorphicKey)
            IF(.NOT.ASSOCIATED(val)) THEN
              resultVar = .FALSE.
              EXIT outerLoop
            END IF
          END DO
        END DO outerLoop
    END FUNCTION InputContainer_dataAvailable

    LOGICAL FUNCTION InputContainer_fetch2D(me, level, tile, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), INTENT(IN) :: level
        INTEGER, INTENT(IN) :: tile
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetch2D"
        TYPE(t_levelkey), TARGET :: key
        CLASS(*), POINTER :: polymorphicKey, val
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug
        IF(debugInfo) CALL message(routine, "level = "//TRIM(real2string(level))//", tile = "//TRIM(int2string(tile)))
        resultVar = me%dataAvailable(opt_level = level, opt_tile = tile)
        IF(.NOT.resultVar) RETURN

        key%levelValue = level
        key%tileId = tile
        polymorphicKey => key
        val => me%fields%getEntry(polymorphicKey)

        IF(.NOT.ASSOCIATED(val)) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but getEntry() failed")
        SELECT TYPE(val)
            TYPE IS(t_LevelPointer)
              IF (ANY(SHAPE(outData) /= SHAPE(val%ptr))) THEN
                   CALL finish(routine, "dimensions of output array do not match the dimensions of the data read from file")
                END IF
                outData(:,:) = val%ptr(:,:)
            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT

    END FUNCTION InputContainer_fetch2D

    LOGICAL FUNCTION InputContainer_fetch3D(me, tile, outData, optLevelDimension, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputContainer), INTENT(IN) :: me
        INTEGER, VALUE :: tile
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, INTENT(IN), OPTIONAL :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetch3D"
        INTEGER :: levelDimension, i
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        levelDimension = 2
        IF (PRESENT(optLevelDimension)) levelDimension = optLevelDimension
        IF(SIZE(outData, levelDimension) /= me%levels%valueCount) THEN
            resultVar = .FALSE.
            IF(debugInfo) CALL message(routine, "wrong level count")
            RETURN
        END IF

        resultVar = me%dataAvailable(opt_tile = tile)
        IF(.NOT.resultVar) RETURN

        DO i = 1, me%levels%valueCount
            SELECT CASE(levelDimension)
                CASE(2)
                    resultVar = me%fetch2D(me%levels%values(i), tile, outData(:, i, :))
                CASE(3)
                    resultVar = me%fetch2D(me%levels%values(i), tile, outData(:, :, i))
                CASE DEFAULT
                    CALL finish(routine, "illegal value in optleveldimension")
            END SELECT
            IF(.NOT. resultVar) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but a fetch2D() CALL failed")
        END DO
    END FUNCTION InputContainer_fetch3D

    LOGICAL FUNCTION InputContainer_fetchTiled2D(me, level, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), INTENT(IN) :: level
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetchTiled2D"
        INTEGER :: tileId
        LOGICAL :: debugInfo, useUntiledData
        TYPE(t_tileinfo_icon) :: trivial_tileinfo_icon

        useUntiledData = .FALSE.
        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug
        resultVar = me%dataAvailable(opt_level = level)
        IF(.NOT.resultVar) RETURN

        trivial_tileinfo_icon = trivial_tile_att%getTileinfo_icon()

        !check whether we have the right count of tiles, ignoring untiled fields of the same variable that may be PRESENT as well
        IF(me%tiles%haveValue(REAL(trivial_tileinfo_icon%idx, dp))) THEN
            resultVar = SIZE(outData, 3) == me%tiles%valueCount - 1
            IF(.NOT.resultVar.AND.SIZE(outData, 3) == 1) THEN
                CALL message(routine, "Warning: No tiled input data is present. Reading single tile from untiled data instead.")
                useUntiledData = .TRUE.
                resultVar = .TRUE.
            END IF
            IF(debugInfo .AND. .NOT. resultVar) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount - 1))//" &
                                      &(+ one untiled field), &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 3))))
            END IF
        ELSE
            resultVar = SIZE(outData, 3) == me%tiles%valueCount
            IF(debugInfo .AND. .NOT. resultVar) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount))//" &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 3))))
            END IF
        END IF
        IF(.NOT. resultVar) RETURN

        DO tileId = 1, SIZE(outData, 3)
            IF(useUntiledData) THEN
              trivial_tileinfo_icon = trivial_tile_att%getTileinfo_icon()
              resultVar = me%fetch2D(level, trivial_tileinfo_icon%idx, outData(:, :, tileId))
            ELSE
              resultVar = me%fetch2D(level, tileId, outData(:, :, tileId))
            END IF
            IF(.NOT. resultVar) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but a fetch2D() CALL failed")
        END DO
    END FUNCTION InputContainer_fetchTiled2D

    LOGICAL FUNCTION InputContainer_fetchTiled3D(me, outData, optLevelDimension, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, INTENT(IN), OPTIONAL :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetchTiled3D"
        INTEGER :: tileId
        LOGICAL :: debugInfo, useUntiledData
        TYPE(t_tileinfo_icon) :: trivial_tileinfo_icon

        useUntiledData = .FALSE.
        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        !XXX: This IS more strict than strictly necessary:
        !     It will check that untiled DATA that may also be PRESENT has the same amount of levels as the tiled DATA.
        !     I guess, this additional strictness shouldn't be a problem, but I may be wrong.
        resultVar = me%dataAvailable()
        IF(.NOT.resultVar) RETURN

        !check whether we have the right count of tiles, ignoring
        !untiled fields of the same variable that may be PRESENT as
        !well
        trivial_tileinfo_icon = trivial_tile_att%getTileinfo_icon()

        IF(me%tiles%haveValue(REAL(trivial_tileinfo_icon%idx, dp))) THEN
            resultVar = SIZE(outData, 4) == me%tiles%valueCount - 1
            IF(.NOT.resultVar.AND.SIZE(outData, 4) == 1) THEN
                CALL message(routine, "Warning: No tiled input data is present. Reading single tile from untiled data instead.")
                useUntiledData = .TRUE.
                resultVar = .TRUE.
            END IF
            IF(debugInfo .AND. .NOT. resultVar) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount - 1))//" &
                                      &(+ one untiled field), &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 4))))
            END IF
        ELSE
            resultVar = SIZE(outData, 4) == me%tiles%valueCount
            IF(debugInfo .AND. .NOT. resultVar) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount))//" &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 4))))
            END IF
        END IF
        IF(.NOT. resultVar) RETURN

        DO tileId = 1, SIZE(outData, 4)
            IF(useUntiledData) THEN
                resultVar = me%fetch3D(trivial_tileinfo_icon%idx, outData(:, :, :, tileId), optLevelDimension)
            ELSE
                resultVar = me%fetch3D(tileId, outData(:, :, :, tileId), optLevelDimension)
            END IF
            IF(.NOT. resultVar) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but a fetch3D() CALL failed")
        END DO
    END FUNCTION InputContainer_fetchTiled3D

    SUBROUTINE InputContainer_destruct(me)
        CLASS(t_InputContainer) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_destruct"

        CALL me%fields%destruct()
        DEALLOCATE(me%fields)
    END SUBROUTINE InputContainer_destruct

    !NEC_RP: new routine that OMP parallelizes read, statistics and distribution
    SUBROUTINE InputContainer_readField_omp(me, variableName, level, tile, timer, jg, iterator, statistics, iread)
        CLASS(t_InputContainer), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: variableName
        REAL(dp), INTENT(IN) :: level
        INTEGER, INTENT(IN) :: tile, jg
        REAL(dp), INTENT(INOUT) :: timer(:)
        TYPE(t_CdiIterator), VALUE :: iterator
        TYPE(t_Statistics), INTENT(INOUT) :: statistics ! This gets the statistics of the READ field added, but ONLY on the master process.
        INTEGER, INTENT(IN) :: iread

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_readField_omp"
        CLASS(*), POINTER :: val, key
        REAL(C_DOUBLE), SAVE, POINTER :: bufferD(:)
        REAL(C_FLOAT), SAVE, POINTER :: bufferS(:)
        INTEGER :: gridSize, datatype, packedMessage(2), error    !packedMessage(1) == gridSize, packedMessage(2) == datatype
        CLASS(t_ScatterPattern), POINTER :: distribution
        REAL(dp) :: message(3)
        ! Variables for deferred sending
        CHARACTER(LEN=:), ALLOCATABLE, SAVE :: variableName_prev
        CLASS(*), POINTER :: key_prev
        REAL(C_DOUBLE), SAVE, POINTER :: bufferD_prev(:)
        REAL(C_FLOAT), SAVE, POINTER :: bufferS_prev(:)
        REAL(C_DOUBLE), POINTER :: tmpDP(:)
        REAL(C_FLOAT), POINTER :: tmpSP(:)
        INTEGER, SAVE :: packedMessage_prev(2)
        REAL(dp), SAVE :: level_prev
        INTEGER, SAVE :: jg_prev
        INTEGER, SAVE :: tile_prev, datatype_prev

        !NEC_RP: use three cases: only read, only distribute, both
        SELECT CASE(iread)
        CASE(1)  ! First record: only read

           CALL read_data(firstcall = .TRUE.)
   
           ! prepare reading the rest of the data
           jg_prev = jg
           level_prev = level
           tile_prev = tile
           datatype_prev = datatype
           packedMessage_prev = packedMessage
           IF (ALLOCATED(variableName_prev)) DEALLOCATE(variableName_prev)
           ALLOCATE(CHARACTER(LEN=LEN(variableName)) :: variableName_prev)
           variableName_prev = variableName

        CASE (-1)
           !NEC_RP: Last record of previous variable

!$OMP PARALLEL SECTIONS NUM_THREADS(2)
!$OMP SECTION
            ! first section for distribution of previous buffer
            CALL distribute_data()
!$OMP SECTION
            ! second section for statistics on previous buffer
            CALL compute_statistics()
!$OMP END PARALLEL SECTIONS

        CASE (-2)
           !NEC_RP: Last record: Only data distribution and statistics; deallocate

!$OMP PARALLEL SECTIONS NUM_THREADS(2)
!$OMP SECTION
            ! first section for distribution of previous buffer
            CALL distribute_data()
!$OMP SECTION
            ! second section for statistics on previous buffer
            CALL compute_statistics()
!$OMP END PARALLEL SECTIONS

            IF (ASSOCIATED(bufferD_prev)) THEN
              DEALLOCATE(bufferD_prev)
              NULLIFY(bufferD_prev)
            END IF
            IF (ASSOCIATED(bufferD)) THEN 
              DEALLOCATE(bufferD)
              NULLIFY(bufferD)
            END IF  
            IF (ASSOCIATED(bufferS_prev)) THEN
               DEALLOCATE(bufferS_prev)
               NULLIFY(bufferS_prev)
            END IF
            IF (ASSOCIATED(bufferS)) THEN 
              DEALLOCATE(bufferS)
              NULLIFY(bufferS)
            END IF  
            
        CASE DEFAULT
           !NEC_RP: Read data, distribute them, and compute statistics
!$OMP PARALLEL SECTIONS NUM_THREADS(3)
!$OMP SECTION
           ! first section for read of current buffer
           CALL read_data(firstcall = .FALSE.)
!$OMP SECTION
           ! second and third sections for statistics on and distribution of previous buffer
            CALL compute_statistics()
!$OMP SECTION
            CALL distribute_data()
!$OMP END PARALLEL SECTIONS

            ! Prepare data for next iteration
            ! first flip pointer assignments (faster than copy)
            IF (ASSOCIATED(bufferD)) THEN
               IF (SIZE(bufferD) /= SIZE(bufferD_prev)) THEN
                  IF (ASSOCIATED(bufferD_prev)) DEALLOCATE(bufferD_prev)
                  ALLOCATE(bufferD_prev(SIZE(bufferD)))
               END IF
               tmpDP => bufferD_prev
               bufferD_prev => bufferD
               bufferD => tmpDP
            END IF
            IF (ASSOCIATED(bufferS)) THEN
               IF (SIZE(bufferS) /= SIZE(bufferS_prev)) THEN
                  IF (ASSOCIATED(bufferS_prev)) DEALLOCATE(bufferS_prev)
                  ALLOCATE(bufferS_prev(SIZE(bufferS)))
               END IF
               tmpSP => bufferS_prev
               bufferS_prev => bufferS
               bufferS => tmpSP
            END IF

            jg_prev = jg
            level_prev = level
            tile_prev = tile
            datatype_prev = datatype
            packedMessage_prev = packedMessage
            IF (LEN(variableName_prev) /= LEN(variableName)) THEN
               DEALLOCATE(variableName_prev)
               ALLOCATE(CHARACTER(LEN=LEN(variableName)) :: variableName_prev)
            END IF
            variableName_prev = variableName
        END SELECT

        CONTAINS

        SUBROUTINE read_data(firstcall)

           LOGICAL, INTENT(in) :: firstcall
           REAL(dp) :: savetime

           ALLOCATE(t_LevelKey :: key, STAT = error)
           IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
           SELECT TYPE(key)
           TYPE IS(t_LevelKey)
                key%levelValue = level
                key%tileId = tile
             CLASS DEFAULT
                CALL finish(routine, "assertion failed")
           END SELECT

           !Inquire buffer SIZE information AND broadcast it.
           IF(C_ASSOCIATED(iterator%ptr)) THEN
               gridSize = gridInqSize(cdiIterator_inqGridId(iterator))
               datatype = cdiIterator_inqDatatype(iterator)
               packedMessage(1) = gridSize
               packedMessage(2) = datatype
           END IF
   
           SELECT CASE(datatype)
           CASE(CDI_DATATYPE_PACK23:CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32)
               !ALLOCATE the global buffer
               IF(C_ASSOCIATED(iterator%ptr)) THEN
                   IF (firstcall) ALLOCATE(bufferD_prev(gridSize), STAT = error)
                   ! Note: SIZE(buffer) may deliver the right value even if "buffer" has been deallocated before
                   IF (.NOT. ASSOCIATED(bufferD)) ALLOCATE(bufferD(gridSize), STAT = error)
                   IF (SIZE(bufferD) /= gridSize) THEN
                      DEALLOCATE(bufferD)
                      ALLOCATE(bufferD(gridSize), STAT = error)
                   END IF
               ELSE
                   IF (firstcall) THEN
                     ALLOCATE(bufferD_prev(1), STAT = error)
                     ALLOCATE(bufferD(1), STAT = error)
                   ENDIF
                   IF (SIZE(bufferD) /= 1) THEN
                      DEALLOCATE(bufferD)
                      ALLOCATE(bufferD(1), STAT = error)
                   END IF
               END IF
               IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
   
               !READ the DATA
               IF(C_ASSOCIATED(iterator%ptr)) THEN
                   savetime = p_mpi_wtime()
                   IF (firstcall) THEN
                     CALL cdiIterator_readField(iterator, bufferD_prev)
                   ELSE
                     CALL cdiIterator_readField(iterator, bufferD)
                   ENDIF
                   timer(3) = timer(3) + p_mpi_wtime() - savetime
               END IF
           CASE DEFAULT
               !ALLOCATE the global buffer
               IF(C_ASSOCIATED(iterator%ptr)) THEN
                   IF (firstcall) ALLOCATE(bufferS_prev(gridSize), STAT = error)
                   ! Note: SIZE(buffer) may deliver the right value even if "buffer" has been deallocated before
                   IF (.NOT. ASSOCIATED(bufferS)) ALLOCATE(bufferS(gridSize), STAT = error)
                   IF (SIZE(bufferS) /= gridSize) THEN
                      DEALLOCATE(bufferS)
                      ALLOCATE(bufferS(gridSize), STAT = error)
                   END IF
               ELSE
                   IF (firstcall) THEN
                     ALLOCATE(bufferS_prev(1), STAT = error)
                     ALLOCATE(bufferS(1), STAT = error)
                   ENDIF
                   IF (SIZE(bufferS) /= 1) THEN
                      DEALLOCATE(bufferS)
                      ALLOCATE(bufferS(1), STAT = error)
                   END IF
               END IF
               IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
   
               !READ the DATA
               IF(C_ASSOCIATED(iterator%ptr)) THEN
                   savetime = p_mpi_wtime()
                   IF (firstcall) THEN
                     CALL cdiIterator_readFieldF(iterator, bufferS_prev)
                   ELSE
                     CALL cdiIterator_readFieldF(iterator, bufferS)
                   ENDIF
                   timer(3) = timer(3) + p_mpi_wtime() - savetime
               END IF
           END SELECT

        END SUBROUTINE read_data

        SUBROUTINE distribute_data

           REAL(dp) :: savetime

           message(1) = REAL(LEN(variableName_prev, 1), dp)
           message(2) = level_prev
           message(3) = REAL(tile_prev, dp)
           CALL p_bcast(message, process_mpi_root_id, p_comm_work)
           CALL p_bcast(variableName_prev, process_mpi_root_id, p_comm_work)

           CALL p_bcast(packedMessage_prev, 0, p_comm_work)

           !Get the corresponding ScatterPattern AND initialize our hash table entry.
           distribution => lookupScatterPattern(jg_prev, packedMessage_prev(1))
           IF(.NOT. ASSOCIATED(distribution)) CALL finish(routine, "could not find scatter pattern to distribute input field")
           ALLOCATE(t_LevelPointer :: val, STAT = error)
           IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
           SELECT TYPE(val)
             TYPE IS(t_LevelPointer)

             !ALLOCATE the local buffer
             ALLOCATE(val%ptr(nproma, blk_no(distribution%localSize())), STAT = error)
             IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
             val%ptr = 0.0_wp     !Avoid nondeterministic values IN the unused points for checksumming.

             SELECT CASE(packedMessage_prev(2))
             CASE(CDI_DATATYPE_PACK23:CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32)
               savetime = p_mpi_wtime()
               CALL distribution%distribute(bufferD_prev, val%ptr, .FALSE.)
               timer(5) = timer(5) + p_mpi_wtime() - savetime
             CASE DEFAULT
               savetime = p_mpi_wtime()
               CALL distribution%distribute(bufferS_prev, val%ptr, .FALSE.)
               timer(5) = timer(5) + p_mpi_wtime() - savetime
             END SELECT
           END SELECT

           ALLOCATE(t_LevelKey :: key_prev, STAT = error)
           IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
           SELECT TYPE(key_prev)
             TYPE IS(t_LevelKey)
                key_prev%levelValue = level_prev
                key_prev%tileId = tile_prev
             CLASS DEFAULT
                CALL finish(routine, "assertion failed")
           END SELECT
           
           !store the DATA IN our hash table
           CALL me%fields%setEntry(key_prev, val) !this will DEALLOCATE both the val AND the key eventually
           me%fieldCount = me%fieldCount + 1
           CALL me%tiles%addValue(REAL(tile_prev, dp))
           CALL me%levels%addValue(level_prev)

        END SUBROUTINE distribute_data


        SUBROUTINE compute_statistics

           REAL(dp) :: savetime

           SELECT CASE(packedMessage_prev(2))
           CASE(CDI_DATATYPE_PACK23:CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32)
               savetime = p_mpi_wtime()
               CALL statistics%add(bufferD_prev)
               timer(4) = timer(4) + p_mpi_wtime() - savetime
           CASE DEFAULT
               savetime = p_mpi_wtime()
               CALL statistics%add(bufferS_prev)
               timer(4) = timer(4) + p_mpi_wtime() - savetime
           END SELECT

        END SUBROUTINE compute_statistics

    END SUBROUTINE InputContainer_readField_omp


    SUBROUTINE InputContainer_readField(me, variableName, level, tile, timer, jg, iterator, statistics)
        CLASS(t_InputContainer), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: variableName
        REAL(dp), INTENT(IN) :: level
        INTEGER, INTENT(IN) :: tile, jg
        REAL(dp), INTENT(INOUT) :: timer(:)
        TYPE(t_CdiIterator), VALUE :: iterator
        TYPE(t_Statistics), INTENT(INOUT) :: statistics ! This gets the statistics of the READ field added, but ONLY on the master process.

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_readField"
        CLASS(*), POINTER :: key, val
        REAL(C_DOUBLE), POINTER :: bufferD(:)
        REAL(C_FLOAT), POINTER :: bufferS(:)
        INTEGER :: gridSize, datatype, packedMessage(2), error    !packedMessage(1) == gridSize, packedMessage(2) == datatype
        CLASS(t_ScatterPattern), POINTER :: distribution
        REAL(dp) :: savetime

        !sanity check: fail IF this field has already been READ
        ALLOCATE(t_LevelKey :: key, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        SELECT TYPE(key)
            TYPE IS(t_LevelKey)
                key%levelValue = level
                key%tileId = tile
            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT
        IF(ASSOCIATED(me%fields%getEntry(key))) THEN
            WRITE(message_text, '(a,g24.15e3,a,i2,3a)') "double definition of level-tile tuple (", &
              level,",",tile,") in variable '",variableName,"' in an input file"
            CALL finish(routine, message_text)
        END IF

        !Inquire buffer SIZE information AND broadcast it.
        IF(C_ASSOCIATED(iterator%ptr)) THEN
            gridSize = gridInqSize(cdiIterator_inqGridId(iterator))
            datatype = cdiIterator_inqDatatype(iterator)
            packedMessage(1) = gridSize
            packedMessage(2) = datatype
        END IF
        CALL p_bcast(packedMessage, 0, p_comm_work)
        gridSize = packedMessage(1)
        datatype = packedMessage(2)

        !Get the corresponding ScatterPattern AND initialize our hash table entry.
        distribution => lookupScatterPattern(jg, gridSize)
        IF(.NOT. ASSOCIATED(distribution)) CALL finish(routine, "could not find scatter pattern to distribute input field")
        ALLOCATE(t_LevelPointer :: val, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        SELECT TYPE(val)
            TYPE IS(t_LevelPointer)

                !ALLOCATE the local buffer
                ALLOCATE(val%ptr(nproma, blk_no(distribution%localSize())), STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
                val%ptr(:,:) = 0.0_wp     !Avoid nondeterministic values IN the unused points for checksumming.
                SELECT CASE(datatype)
                    CASE(CDI_DATATYPE_PACK23:CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32)
                        !ALLOCATE the global buffer
                        IF(C_ASSOCIATED(iterator%ptr)) THEN
                            ALLOCATE(bufferD(gridSize), STAT = error)
                        ELSE
                            ALLOCATE(bufferD(1), STAT = error)
                        END IF
                        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")

                        !READ the DATA
                        IF(C_ASSOCIATED(iterator%ptr)) THEN
                            savetime = p_mpi_wtime()
                            CALL cdiIterator_readField(iterator, bufferD)
                            timer(3) = timer(3) + p_mpi_wtime() - savetime
                            savetime = p_mpi_wtime()
                            CALL statistics%add(bufferD)
                            timer(4) = timer(4) + p_mpi_wtime() - savetime
                        END IF
                        savetime = p_mpi_wtime()
                        CALL distribution%distribute(bufferD, val%ptr(:, :), .FALSE.)
                        timer(5) = timer(5) + p_mpi_wtime() - savetime

                        DEALLOCATE(bufferD)

                    CASE DEFAULT
                        !ALLOCATE the global buffer
                        IF(C_ASSOCIATED(iterator%ptr)) THEN
                            ALLOCATE(bufferS(gridSize), STAT = error)
                        ELSE
                            ALLOCATE(bufferS(1), STAT = error)
                        END IF
                        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")

                        !READ the DATA
                        IF(C_ASSOCIATED(iterator%ptr)) THEN
                            savetime = p_mpi_wtime()
                            CALL cdiIterator_readFieldF(iterator, bufferS)
                            timer(3) = timer(3) + p_mpi_wtime() - savetime
                            savetime = p_mpi_wtime()
                            CALL statistics%add(bufferS)
                            timer(4) = timer(4) + p_mpi_wtime() - savetime
                        END IF
                        savetime = p_mpi_wtime()
                        CALL distribution%distribute(bufferS, val%ptr(:, :), .FALSE.)
                        timer(5) = timer(5) + p_mpi_wtime() - savetime

                        DEALLOCATE(bufferS)

                END SELECT

            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT

        !store the DATA IN our hash table
        CALL me%fields%setEntry(key, val) !this will DEALLOCATE both the val AND the key eventually
        me%fieldCount = me%fieldCount + 1
        CALL me%tiles%addValue(REAL(tile, dp))
        CALL me%levels%addValue(level)
    END SUBROUTINE InputContainer_readField

    SUBROUTINE ValueList_init(me)
        CLASS(t_ValueList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_init"
        INTEGER :: i, error

        IF(ALLOCATED(me%values)) DEALLOCATE(me%values)
        ALLOCATE(me%values(8), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        me%valueCount = 0
        me%values(:) = 0.0
    END SUBROUTINE ValueList_init

    SUBROUTINE ValueList_addValue(me, val)
        CLASS(t_ValueList), INTENT(INOUT) :: me
        REAL(dp), INTENT(IN) :: val

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_addValue"
        INTEGER :: i, j

        !check whether the val is already IN the list
        DO i = 1, me%valueCount
            IF(me%values(i) == val) RETURN    !nothing to DO, we already have this val
            IF(me%values(i) > val) EXIT
        END DO

        !insert the val at the point that i is pointing to
        CALL me%ensureSpace()
        DO j = me%valueCount, i, -1
            me%values(j+1) = me%values(j)
        END DO
        me%values(i) = val
        me%valueCount = me%valueCount + 1
    END SUBROUTINE ValueList_addValue

    FUNCTION ValueList_haveValue(me, val) RESULT(resultVar)
        CLASS(t_ValueList), INTENT(IN) :: me
        REAL(dp), VALUE :: val
        LOGICAL :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_haveValue"
        INTEGER :: i

        resultVar = .TRUE.
        DO i = 1, me%valueCount
            IF(me%values(i) == val) RETURN    !found it, RETURN success
        END DO
        resultVar = .FALSE.
    END FUNCTION ValueList_haveValue

    SUBROUTINE ValueList_ensureSpace(me)
        CLASS(t_ValueList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_ensureSpace"
        INTEGER :: i, error
        REAL(dp), ALLOCATABLE :: temp(:)

        IF(me%valueCount == SIZE(me%values, 1)) THEN
            ALLOCATE(temp(2*SIZE(me%values, 1)), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            DO i = 1, me%valueCount
                temp(i) = me%values(i)
            END DO
            temp(me%valueCount + 1:) = 0.0
            CALL move_ALLOC(temp, me%values)
        END IF
    END SUBROUTINE ValueList_ensureSpace

END MODULE mo_input_container
