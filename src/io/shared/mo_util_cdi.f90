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

! Contains utility routines for reading NetCDF and GRIB2 files (using
! the CDI library) and communicating fields in parallel.

MODULE mo_util_cdi

  USE ISO_C_BINDING,         ONLY: C_INT, C_CHAR
  USE mo_kind,               ONLY: wp, sp, dp, i8
  USE mo_exception,          ONLY: finish
  USE mo_communication,      ONLY: t_scatterPattern
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_run_config,         ONLY: msg_level
  USE mo_io_config,          ONLY: config_lmask_boundary => lmask_boundary
  USE mo_mpi,                ONLY: p_bcast, p_io, my_process_is_stdio, p_mpi_wtime,  &
    &                              my_process_is_mpi_workroot, p_barrier
  USE mo_util_string,        ONLY: tolower, int2string
  USE mo_dictionary,         ONLY: t_dictionary, DICT_MAX_STRLEN
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_gribout_config,     ONLY: t_gribout_config, GRIB_LIB_COMPAT_ECC_2_31_0
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_aes_phy_config,     ONLY: aes_phy_config
  USE mo_grib2_util,         ONLY: set_GRIB2_additional_keys, set_GRIB2_tile_keys, &
    &                              set_GRIB2_ensemble_keys, set_GRIB2_local_keys,  &
    &                              set_GRIB2_synsat_keys, set_GRIB2_chem_keys,     &
    &                              grib_lib_compatibility
  USE mo_nwp_sfc_tiles,      ONLY: t_tileinfo_icon, t_tileinfo_grb2, trivial_tile_att

  USE mo_cdi,                ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, streamInqVlist, vlistNvars, vlistInqVarDatatype, &
                                 & vlistInqVarIntKey, vlistInqVarZaxis, &
                                 & vlistInqVarGrid, gridInqSize, zaxisInqSize, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32, &
                                 & streamInqTimestep, &
                                 & vlistDefVarIntKey, &
                                 & streamReadVarSliceF, streamReadVarSlice, vlistInqVarName, &
                                 & vlistInqVarSubtype, subtypeInqSize, &
                                 & cdi_max_name, &
                                 & subtypeDefActiveIndex, CDI_DATATYPE_PACK23, CDI_DATATYPE_PACK32, cdiStringError, &
                                 & FILETYPE_GRB2, vlistDefVar, cdiEncodeParam, &
                                 & vlistDefVarName, vlistDefVarLongname, vlistDefVarTsteptype, &
    &                              vlistDefVarStdname, vlistDefVarUnits, vlistDefVarParam, vlistDefVarMissval, &
    &                              vlistDefVarDatatype, vlistDefVarIntKey, vlistDefVarDblKey, &
                                 & TIME_CONSTANT, TIME_VARYING, TSTEP_CONSTANT

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: has_filetype_netcdf
  PUBLIC :: read_cdi_2d, read_cdi_3d
  PUBLIC :: test_cdi_varID
  PUBLIC :: get_cdi_varID
  PUBLIC :: get_cdi_NlevRef
  PUBLIC :: t_inputParameters, makeInputParameters, deleteInputParameters
  PUBLIC :: cdiGetStringError
  PUBLIC :: create_cdi_variable

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_cdi'

  INTERFACE read_cdi_2d
    MODULE PROCEDURE read_cdi_2d_real
    MODULE PROCEDURE read_cdi_2d_int
    MODULE PROCEDURE read_cdi_2d_time
    MODULE PROCEDURE read_cdi_2d_real_tiles
    MODULE PROCEDURE read_cdi_2d_int_tiles
    MODULE PROCEDURE read_cdi_2d_time_tiles
    MODULE PROCEDURE read_cdi_2d_lbc
  END INTERFACE

  INTERFACE read_cdi_3d
    MODULE PROCEDURE read_cdi_3d_real
    MODULE PROCEDURE read_cdi_3d_real_tiles
    MODULE PROCEDURE read_cdi_3d_lbc
  END INTERFACE


  TYPE t_tileinfo
    TYPE(t_tileinfo_grb2), ALLOCATABLE :: tile(:)  !< variable specific
    !tile indices > CDI internal, one-dimensional index for the list
    !of (idx/att) > pairs:
    INTEGER,              ALLOCATABLE :: tile_index(:)
  END TYPE t_tileinfo



  ! This is a small type that serves two functions:
  ! 1. It encapsulates three parameters to the read functions into one, significantly reducing the hassle to call them.
  ! 2. It allows some information about the file to be distributed to the non-I/O PEs in a single broadcast that would otherwise
  !    require a broadcast in every call.
  !
  ! Create this by calling makeInputParameters(), destroy with deleteInputParameters().
  TYPE :: t_inputParameters
    !PRIVATE! Don't use outside of this file!
    INTEGER                                     :: streamId     !< CDI stream ID
    INTEGER                                     :: glb_arr_len  !< global array length
    TYPE (t_dictionary)                         :: dict     !< a dictionary that is used to translate variable names
    CLASS(t_scatterPattern), POINTER            :: distribution    !< a t_scatterPattern to distribute the data
    LOGICAL                                     :: have_dict    !< whether `dict` is used or not
    CHARACTER(LEN=MAX_CHAR_LENGTH), ALLOCATABLE :: variableNames(:) !< the names of the variables
    INTEGER, ALLOCATABLE                        :: variableDatatype(:) !< the datatypes that are used by the file to store the variables
    TYPE(t_tileinfo), ALLOCATABLE               :: variableTileinfo(:) !< variable specific tile meta info

    REAL(dp) :: readDuration    !< statistic on how long we took to read the data
    INTEGER(i8) :: readBytes    !< statistic on how much data we read from the file

  CONTAINS
    PROCEDURE :: findVarId => inputParametersFindVarId  !< determine the ID of an named variable
    PROCEDURE :: inqVarId => inputParametersInqVarId
    PROCEDURE :: inqVarIdTiles => inputParametersInqVarIdTiles
  END TYPE

CONTAINS

  !> Provides a common interface to all NetCDF flavors used within
  !> ICON.
  !
  FUNCTION has_filetype_netcdf(cdi_filetype)
    LOGICAL :: has_filetype_netcdf
    INTEGER, INTENT(IN) :: cdi_filetype

    has_filetype_netcdf = ((cdi_filetype == FILETYPE_NC)  .OR. &
      &                    (cdi_filetype == FILETYPE_NC2) .OR. &
      &                    (cdi_filetype == FILETYPE_NC4))
  END FUNCTION has_filetype_netcdf


  !---------------------------------------------------------------------------------------------------------------------------------
  !> Creates a t_inputParameters object for the given arguments.
  !---------------------------------------------------------------------------------------------------------------------------------
  FUNCTION makeInputParameters(streamId, glb_arr_len, distribution, opt_dict) RESULT(me)
    IMPLICIT NONE
    TYPE(t_inputParameters) :: me
    INTEGER, INTENT(IN) :: streamID, glb_arr_len
    CLASS(t_scatterPattern), POINTER, INTENT(IN) :: distribution
    TYPE(t_dictionary), INTENT(IN), OPTIONAL :: opt_dict

    CHARACTER(len=*), PARAMETER :: routine = modname//':makeInputParameters'
    INTEGER :: vlistId, variableCount, i, j, ierrstat
    INTEGER :: subtypeID
    INTEGER :: ientry
    INTEGER :: cnt

    INTEGER, ALLOCATABLE :: tileIdxAttTidTemp(:,:)
    INTEGER, ALLOCATABLE :: subtypeSize(:)

    INTEGER :: idx, att, tile_index
    TYPE(t_tileinfo_icon) :: tileinfo_icon
    LOGICAL :: is_workroot

    is_workroot = my_process_is_mpi_workroot()


    !first forward the arguments to the object we are building
    me%streamId = streamId
    me%glb_arr_len = glb_arr_len

    me%have_dict = .FALSE.
    IF(PRESENT(opt_dict)) THEN
      CALL me%dict%init(.FALSE.)
      CALL opt_dict%copy(me%dict)
      me%have_dict = .TRUE.
    END IF
    me%distribution => distribution
    CALL me%distribution%resetStatistics()

    !now the interesting part: introspect the file and broadcast the variable info needed to avoid broadcasting it later.
    IF (is_workroot) THEN
        vlistId = streamInqVlist(streamId)
        variableCount = vlistNvars(vlistId)
    END IF
    CALL p_bcast(variableCount, p_io, distribution%communicator)

    ALLOCATE(me%variableNames(variableCount), &
      &      me%variableDatatype(variableCount), &
      &      me%variableTileinfo(variableCount), &
      &      subtypeSize(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    subtypeSize(1:variableCount) = 0

    IF (is_workroot) THEN
        do i = 1, variableCount
            CALL vlistInqVarName(vlistId, i-1, me%variableNames(i))
            me%variableDatatype(i) = vlistInqVarDatatype(vlistId, i-1)

            subtypeID = vlistInqVarSubtype(vlistID,i-1)
            subtypeSize(i) = subtypeInqSize(subtypeID)
            IF(subtypeSize(i) < 1) subtypeSize(i) = 1   !We need this to ensure that the trivial_tileinfo is actually stored somewhere IN the CASE of non-tile variables. Otherwise, non-tile variables can never be found.

            ALLOCATE(me%variableTileinfo(i)%tile(subtypeSize(i)), &
              &      me%variableTileinfo(i)%tile_index(subtypeSize(i)), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

            IF (.not. ANY(aes_phy_config(:)%ljsb)) THEN
              IF (vlistInqVarIntKey(vlistId, i-1, "totalNumberOfTileAttributePairs") <= 0) THEN
                ! not a tile variable
                me%variableTileinfo(i)%tile(:) = trivial_tile_att%getTileinfo_grb2()
                tileinfo_icon = trivial_tile_att%getTileinfo_icon()

                me%variableTileinfo(i)%tile_index(:) = tileinfo_icon%idx
              ELSE
                ! tile
                DO ientry=1, subtypeSize(i)
                  CALL subtypeDefActiveIndex(subtypeID,ientry-1)  ! starts with 0

                  idx        = vlistInqVarIntKey(vlistId, i-1, "tileIndex")
                  att        = vlistInqVarIntKey(vlistId, i-1, "tileAttribute")
                  tile_index = ientry-1

                  me%variableTileinfo(i)%tile(ientry)       = t_tileinfo_grb2( idx, att )
                  me%variableTileinfo(i)%tile_index(ientry) = tile_index
                END DO
                ! reset active index
                CALL subtypeDefActiveIndex(subtypeID,0)
              ENDIF  ! totalNumberOfTileAttributePairs <= 0
            ENDIF
          ENDDO
    ENDIF
    
    CALL p_bcast(subtypeSize, p_io, distribution%communicator)

    ! put tile info into local 1D arrays, for broadcasting
    ! currently direct bradcast of variable variableTileinfo not possible
    ALLOCATE(tileIdxAttTidTemp(SUM(subtypeSize(1:variableCount)),3), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    IF (is_workroot) THEN
      cnt = 0
      DO i=1, variableCount
        DO j = 1, subtypeSize(i)
          tileIdxAttTidTemp(cnt+j, 1) = me%variableTileinfo(i)%tile(j)%idx
          tileIdxAttTidTemp(cnt+j, 2) = me%variableTileinfo(i)%tile(j)%att
          tileIdxAttTidTemp(cnt+j, 3) = me%variableTileinfo(i)%tile_index(j)
        END DO
        cnt = cnt + subtypeSize(i)
      END DO
    ENDIF

    CALL p_bcast(me%variableNames   , p_io, distribution%communicator)
    CALL p_bcast(me%variableDatatype, p_io, distribution%communicator)
    CALL p_bcast(tileIdxAttTidTemp  , p_io, distribution%communicator)

    ! read tile info from broadcasted 1D array and store in array variableTileinfo of TYPE t_tileinfo
    IF (.NOT. is_workroot) THEN
      cnt = 0
      DO i=1, variableCount
        IF (subtypeSize(i) > 0) THEN
          ALLOCATE(me%variableTileinfo(i)%tile(subtypeSize(i)), &
            &      me%variableTileinfo(i)%tile_index(subtypeSize(i)), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
          DO j = 1, subtypeSize(i)
            me%variableTileinfo(i)%tile(j)%idx   = tileIdxAttTidTemp(cnt+j, 1)
            me%variableTileinfo(i)%tile(j)%att   = tileIdxAttTidTemp(cnt+j, 2)
            me%variableTileinfo(i)%tile_index(j) = tileIdxAttTidTemp(cnt+j, 3)
          END DO
          cnt = cnt + subtypeSize(i)
        ENDIF
      ENDDO
    ENDIF

    me%readDuration = 0.0_dp
    me%readBytes = 0_i8

    ! cleanup
    DEALLOCATE(tileIdxAttTidTemp, subtypeSize, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END FUNCTION makeInputParameters

  ! This function is only a workaround for a compiler bug on the
  ! blizzard. May be reintegrated into inputParametersFindVarId() once
  ! we are not concerned about xlf anymore.
  PURE FUNCTION compareTiledVars(name1, idx1, att1, name2, idx2, att2) &
       RESULT(is_equal)
    CHARACTER(LEN = *), INTENT(IN) :: name1, name2
    INTEGER, INTENT(in) :: idx1, att1, idx2, att2
    LOGICAL :: is_equal
    is_equal = idx1 == idx2 .AND. att1 == att2 .AND. tolower(name1) == name2
  END FUNCTION compareTiledVars

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Try to find the ID of the given tile of the given variable in the input file.
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE inputParametersInqVarIdTiles(me, name, tileinfo, varID, tile_index)
    CLASS(t_inputParameters), INTENT(IN)  :: me !< Input parameters.
    CHARACTER(len=*),         INTENT(IN)  :: name !< Name of the variable.
    TYPE(t_tileinfo_grb2),    INTENT(IN)  :: tileinfo !< GRIB tile information.
    INTEGER,                  INTENT(OUT) :: varID !< CDI variable ID.
    INTEGER,                  INTENT(OUT) :: tile_index !< CDI tile index.

    CHARACTER(len=DICT_MAX_STRLEN) :: mapped_name
    INTEGER :: i, j, tlen

    mapped_name = name
    IF (me%have_dict) &
      & mapped_name = me%dict%get(TRIM(name), DEFAULT=name)
    mapped_name = tolower(mapped_name)
    tlen = LEN_TRIM(mapped_name)

    varID      = -1
    tile_index = -1

    DO i = 1, SIZE(me%variableNames)
      DO j = 1, SIZE(me%variableTileinfo(i)%tile)
        IF (compareTiledVars(me%variableNames(i), me%variableTileinfo(i)%tile(j)%idx, me%variableTileinfo(i)%tile(j)%att, &
          &                  mapped_name(1:tlen), tileinfo%idx, tileinfo%att)) THEN
          varID      = i-1
          tile_index = me%variableTileinfo(i)%tile_index(j)
          RETURN
        END IF
      END DO
    END DO

  END SUBROUTINE inputParametersInqVarIdTiles

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Try to find the ID of the given variable in the input file.
  !---------------------------------------------------------------------------------------------------------------------------------
  FUNCTION inputParametersInqVarId(me, name) RESULT(varID)
    CLASS(t_inputParameters), INTENT(IN)  :: me !< Input parameters.
    CHARACTER(len=*),         INTENT(IN)  :: name !< Name of the variable.
    INTEGER                               :: varID !< CDI variable ID.

    INTEGER :: tile_index

    CALL me%inqVarIdTiles(name, trivial_tile_att%getTileinfo_grb2(), varID, tile_index)

  END FUNCTION inputParametersInqVarId

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Determine the datatype of the given variable in the input file
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE inputParametersFindVarId(me, name, tileinfo, varID, tile_index)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(IN)  :: me
    CHARACTER(len=*),         INTENT(IN)  :: name
    TYPE(t_tileinfo_grb2),    INTENT(IN)  :: tileinfo
    INTEGER,                  INTENT(OUT) :: varID, tile_index

    CHARACTER(len=*), PARAMETER :: routine = modname//':inputParametersFindVarId'
    CHARACTER(len=DICT_MAX_STRLEN) :: mapped_name
    INTEGER :: i, j, tlen

    CALL me%inqVarIdTiles(name, tileinfo, varID, tile_index)

    ! insanity check
    IF(varID < 0) THEN
      IF(my_process_is_stdio()) THEN
        mapped_name = name
        IF (me%have_dict) &
          & mapped_name = me%dict%get(TRIM(name), DEFAULT=name)
        mapped_name = tolower(mapped_name)
        tlen = LEN_TRIM(mapped_name)

        PRINT '(5a)', routine, ": mapped_name = '", mapped_name(1:tlen), "'", "", &
             routine, ": tile idx = ", tileinfo%idx, ", att = ", tileinfo%att, &
             routine, ": list of variables:"
        DO i = 1, SIZE(me%variableNames)
          ! the following "tolower" function call should be enclosed
          ! in a TRIM() intrinsic; however, this throws an (erroneous)
          ! compiler error on some Cray compiler versions.
          PRINT '(4a)', routine, ":     '", tolower(me%variableNames(i)), "'"
          DO j = 1, SIZE(me%variableTileinfo(i)%tile)
            PRINT '(5a)', routine, ":         tile idx = ", &
                 me%variableTileinfo(i)%tile(j)%idx, &
                 ", att = ", me%variableTileinfo(i)%tile(j)%att
          END DO
        END DO
      END IF
      CALL finish(routine, "Variable "//name(1:tlen)//" not found!")
    END IF
  END SUBROUTINE inputParametersFindVarId

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Destroys a t_inputParameters object
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE deleteInputParameters(me)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(INOUT) :: me

    CHARACTER(len=*), PARAMETER :: routine = modname//':deleteInputParameters'

    IF(msg_level >= 10 .and. my_process_is_stdio()) then
        WRITE(0,*) routine, ": Total read statistics for stream ID ", me%streamId
        WRITE(0,'(8X,A,I19,A)')   "amount:    ", me%readBytes, " bytes"
        WRITE(0,'(8X,A,F19.3,A)') "duration:  ", me%readDuration, " seconds"
        IF (me%readDuration > 0._wp) THEN
          WRITE(0,'(8X,A,F19.3,A)') "bandwidth: ", REAL(me%readBytes, dp)/(1048576.0_dp*me%readDuration), " MiB/s"
        ENDIF
    END IF
    CALL me%distribution%printStatistics()

    IF(me%have_dict) CALL me%dict%finalize()
    DEALLOCATE(me%variableNames)
    DEALLOCATE(me%variableDatatype)
    DEALLOCATE(me%variableTileinfo)
  END SUBROUTINE deleteInputParameters

  !-------------------------------------------------------------------------
  !> @return CDI variable ID if CDI stream contains a variable of the
  !> given name, -1 otherwise
  !
  !  Uses cdilib for file access.
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  FUNCTION test_cdi_varID(streamID, name, opt_dict) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: streamID            !< link to file
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    CHARACTER(len=CDI_MAX_NAME) :: zname
    INTEGER                         :: nvars, varID, vlistID
    CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name

    IF (PRESENT(opt_dict)) THEN
      ! Search name mapping for name in NetCDF/GRIB2 file
      mapped_name = tolower(opt_dict%get(name, DEFAULT=name))
    ELSE
      mapped_name = tolower(name)
    END IF

    zname   = ""
    vlistID = streamInqVlist(streamID)

    result_varID = -1
    ! total number of available fields:
    nvars = vlistNvars(vlistID)
    ! loop over vlist, find the corresponding varID
    LOOP : DO varID=0,(nvars-1)

      CALL vlistInqVarName(vlistID, varID, zname)
      IF (tolower(zname) == mapped_name) THEN

        result_varID = varID
        EXIT LOOP
      END IF
    END DO LOOP
  END FUNCTION test_cdi_varID


  !-------------------------------------------------------------------------
  !> @return vlist variable ID for a given variable name
  !
  !  Uses cdilib for file access.
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  FUNCTION get_cdi_varID(streamID, name, opt_dict) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: streamID            !< link to file
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_cdi_varID'

    result_varID = test_cdi_varID(streamID, name, opt_dict)
    IF (result_varID== -1) THEN
      CALL finish(routine, "Variable "//TRIM(name)//" not found!")
    END IF
  END FUNCTION get_cdi_varID



  !-------------------------------------------------------------------------
  !> @return the number of levels for a given variable name
  !
  !  Uses cdilib for file access.
  !  Initial revision by D. Reinert, DWD (2014-10-24)
  !
  FUNCTION get_cdi_NlevRef(streamID, name, opt_dict) RESULT(result_NlevRef)
    INTEGER                                 :: result_NlevRef
    INTEGER,           INTENT(IN)           :: streamID            !< link to file
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    INTEGER :: varID                                              ! variable ID
    INTEGER :: vlistID
    INTEGER :: zaxisID
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_cdi_NlevRef'

    varID = test_cdi_varID(streamID, name, opt_dict)
    IF (varID== -1) THEN
      CALL finish(routine, "Variable "//TRIM(name)//" not found!")
    END IF
    vlistID = streamInqVlist(streamID)
    zaxisID = vlistInqVarZaxis(vlistID,varID)
    ! number of half levels of the generalized Z-axis
    result_NlevRef = zaxisInqSize(zaxisID)
  END FUNCTION get_cdi_NlevRef


  !---------------------------------------------------------------------------------------------------------------------------------
  !> small helper message to output the speed of an individual input operation
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE writeSpeedMessage(routineName, bytes, duration)
    CHARACTER(len=*), INTENT(IN) :: routineName
    INTEGER(i8), VALUE :: bytes
    REAL(dp), VALUE :: duration

    IF(msg_level >= 15) THEN
        IF(duration /= 0.0) THEN
            WRITE(0,*) routineName, ": Read ", bytes, " bytes in ", duration, " seconds (", &
                &      REAL(bytes, dp)/(1048576.0_dp*duration), " MiB/s)"
        ELSE
            WRITE(0,*) routineName, ": Read ", bytes, " bytes in no time. ", &
                &      "Failure to measure time may be due to compilation without MPI."
        END IF
    END IF
  END SUBROUTINE


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for streamReadVarSlice() that measures the time
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE timeStreamReadVarSlice(parameters, varID, level, buffer, nmiss)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, level
    REAL(dp), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(OUT) :: nmiss

    CHARACTER(len=*), PARAMETER :: routine = modname//':timeStreamReadVarSlice'
    REAL(dp) :: startTime, duration
    INTEGER(i8) :: bytes

    startTime = p_mpi_wtime()
    CALL streamReadVarSlice(parameters%streamID, varID, level, buffer, nmiss)
    IF(msg_level >= 5 .and. my_process_is_stdio()) THEN
        duration = p_mpi_wtime() - startTime
        bytes = INT(SIZE(buffer, 1), i8) * 8_i8
        parameters%readDuration = parameters%readDuration + duration
        parameters%readBytes = parameters%readBytes + bytes
        CALL writeSpeedMessage(routine, bytes, duration)
    END IF
  END SUBROUTINE timeStreamReadVarSlice


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for streamReadVarSliceF() that measures the time
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE timeStreamReadVarSliceF(parameters, varID, level, buffer, nmiss)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, level
    REAL(sp), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(OUT) :: nmiss

    CHARACTER(len=*), PARAMETER :: routine = modname//':timeStreamReadVarSliceF'
    REAL(dp) :: startTime, duration
    INTEGER(i8) :: bytes

    startTime = p_mpi_wtime()
    CALL streamReadVarSliceF(parameters%streamID, varID, level, buffer, nmiss)
    IF(msg_level >= 5 .and. my_process_is_stdio()) THEN
        duration = p_mpi_wtime() - startTime
        bytes = INT(SIZE(buffer, 1), i8) * 8_i8
        parameters%readDuration = parameters%readDuration + duration
        parameters%readBytes = parameters%readBytes + bytes
        CALL writeSpeedMessage(routine, bytes, duration)
    END IF
  END SUBROUTINE timeStreamReadVarSliceF


  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 3D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_3d_wp(parameters, varID, nlevs, levelDimension, var_out, lvalue_add)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, nlevs, levelDimension
    REAL(wp), TARGET, INTENT(INOUT) :: var_out(:,:,:) !< output field
    LOGICAL, INTENT(IN) :: lvalue_add         !< If .TRUE., add values to given field

    REAL(wp), POINTER :: var_out_lev(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_wp'
    INTEGER :: jk, ierrstat, nmiss
    REAL(wp), ALLOCATABLE :: tmp_buf(:) ! temporary local array
    LOGICAL :: is_workroot

    is_workroot = my_process_is_mpi_workroot()

    ! allocate a buffer for one vertical level
    ALLOCATE(tmp_buf(MERGE(parameters%glb_arr_len, 1, is_workroot)), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! initialize output field:
    !XXX: It stands to reason whether we really need to do this. It's only effective iff there are points in var_out that are not
    !     covered by blk_no/idx_no. But the unconditional version was definitely wrong as it entirely negated the effect of
    !     lvalue_add.
    IF(.not. lvalue_add) var_out(:,:,:) = 0._wp

    IF (leveldimension < 2 .OR. leveldimension > 3) &
      CALL finish(routine, "Internal error!")
    !FIXME: This code is most likely latency bound, not throughput bound. Needs some asynchronicity to hide the latencies.
    DO jk=1,nlevs
      IF (is_workroot) THEN
        ! read record as 1D field
        CALL timeStreamReadVarSlice(parameters, varID, jk-1, tmp_buf(:), nmiss)
      END IF

      SELECT CASE(levelDimension)
      CASE(2)
        var_out_lev => var_out(:, jk, :)
      CASE(3)
        var_out_lev => var_out(:, :, jk)
      END SELECT
      CALL parameters%distribution%distribute(tmp_buf(:), var_out_lev, lvalue_add)
    END DO ! jk=1,nlevs

    ! clean up
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE read_cdi_3d_wp


  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 3D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_3d_sp(parameters, varID, nlevs, levelDimension, var_out, lvalue_add)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, nlevs, levelDimension
    REAL(wp), TARGET, INTENT(INOUT) :: var_out(:,:,:) !< output field
    LOGICAL, INTENT(IN) :: lvalue_add         !< If .TRUE., add values to given field

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_sp'
    INTEGER :: jk, ierrstat, nmiss
    REAL(wp), POINTER :: var_out_lev(:,:)
    REAL(sp), ALLOCATABLE :: tmp_buf(:) ! temporary local array
    LOGICAL :: is_workroot

    is_workroot = my_process_is_mpi_workroot()

    ! allocate a buffer for one vertical level
    ALLOCATE(tmp_buf(MERGE(parameters%glb_arr_len, 1, is_workroot)), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! initialize output field:
    !XXX: It stands to reason whether we really need to do this. It's only effective iff there are points in var_out that are not
    !     covered by blk_no/idx_no. But the unconditional version was definitely wrong as it entirely negated the effect of
    !     lvalue_add.
    IF(.not. lvalue_add) var_out(:,:,:) = 0._wp

    IF (leveldimension < 2 .OR. leveldimension > 3) &
      CALL finish(routine, "Internal error!")
    !FIXME: This code is most likely latency bound, not throughput bound. Needs some asynchronicity to hide the latencies.
    DO jk=1,nlevs
      IF (is_workroot) THEN
        ! read record as 1D field
        CALL timeStreamReadVarSliceF(parameters, varID, jk-1, tmp_buf(:), nmiss)
      END IF

      SELECT CASE(levelDimension)
      CASE(2)
        var_out_lev => var_out(:, jk, :)
      CASE(3)
        var_out_lev => var_out(:, :, jk)
      END SELECT
      CALL parameters%distribution%distribute(tmp_buf(:), var_out_lev, lvalue_add)
    END DO ! jk=1,nlevs

    ! clean up
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE read_cdi_3d_sp

  !-------------------------------------------------------------------------
  !> Read 3D dataset from file.
  !
  !  Note: This implementation uses a 2D buffer.
  !
  SUBROUTINE read_cdi_3d_real_tiles(parameters, varname, nlevs, var_out, tileinfo, opt_lvalue_add, opt_lev_dim)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(wp),                INTENT(INOUT) :: var_out(:,:,:) !< output field
    TYPE(t_tileinfo_grb2),   INTENT(IN)    :: tileinfo
    LOGICAL,             INTENT(IN), OPTIONAL :: opt_lvalue_add       !< If .TRUE., add values to given field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_lev_dim          !< array dimension (of the levels)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_real_tiles'
    INTEGER :: vlistId, varId, zaxisId, gridId, levelDimension, subtypeID, tile_index
    LOGICAL :: lvalue_add

    levelDimension = 2
    IF (PRESENT(opt_lev_dim)) levelDimension = opt_lev_dim
    lvalue_add = .FALSE.
    IF (PRESENT(opt_lvalue_add)) lvalue_add = opt_lvalue_add

    CALL parameters%findVarId(varname, tileinfo, varID, tile_index)

    IF ((tile_index < 0) .AND. (tileinfo%idx > 0)) THEN
      CALL finish(routine, "Requested tile not found!")
    END IF

    IF(my_process_is_mpi_workroot()) THEN
      ! set active tile index, if this is a tile-based variable
      vlistId = streamInqVlist(parameters%streamId)
      subtypeID  = vlistInqVarSubtype(vlistID,varID)
      IF (tile_index > 0)  CALL subtypeDefActiveIndex(subtypeID, tile_index)
      ! sanity check of the variable dimensions
      zaxisId = vlistInqVarZaxis(vlistId, varId)
      gridId = vlistInqVarGrid(vlistId, varId)
      IF ((gridInqSize(gridId) /= parameters%glb_arr_len) .OR.  &
        & (zaxisInqSize(zaxisId) /= nlevs)) THEN
        CALL finish(routine, "Incompatible dimensions! Grid size = "//trim(int2string(gridInqSize(gridId)))//&
        &                    " (expected "//trim(int2string(parameters%glb_arr_len))//"),"//&
        &                    " level count = "//trim(int2string(zaxisInqSize(zaxisId)))//&
        &                    " (expected "//trim(int2string(nlevs))//")")
      END IF
    END IF

    SELECT CASE(parameters%variableDatatype(varId+1))
        CASE(CDI_DATATYPE_PACK23:CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32)
            ! int32 is treated as double precision because single precision floats would cut off up to seven bits from the integer
            CALL read_cdi_3d_wp(parameters, varId, nlevs, levelDimension, var_out, lvalue_add)
        CASE DEFAULT
            ! XXX: Broadcasting CDI_DATATYPE_PACK1..CDI_DATATYPE_PACK22 DATA as single precision may actually change their values, but this
            !      error will always be smaller than the error made by storing the DATA as CDI_DATATYPE_PACK1..CDI_DATATYPE_PACK22 IN the
            !      first place.
            CALL read_cdi_3d_sp(parameters, varId, nlevs, levelDimension, var_out, lvalue_add)
    END SELECT

    IF(my_process_is_mpi_workroot()) THEN
      ! reset tile index
      IF (tile_index > 0)  CALL subtypeDefActiveIndex(subtypeID, 0)
    END IF
  END SUBROUTINE read_cdi_3d_real_tiles


  !-------------------------------------------------------------------------
  !> Read 3D dataset from file.
  !
  !  Note: This implementation uses a 2D buffer.
  !
  SUBROUTINE read_cdi_3d_real(parameters, varname, nlevs, var_out, opt_lvalue_add, opt_lev_dim)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(wp),                INTENT(INOUT) :: var_out(:,:,:) !< output field
    LOGICAL,             INTENT(IN), OPTIONAL :: opt_lvalue_add       !< If .TRUE., add values to given field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_lev_dim          !< array dimension (of the levels)

    CALL read_cdi_3d_real_tiles(parameters, varname, nlevs, var_out, &
      &                         trivial_tile_att%getTileinfo_grb2(), opt_lvalue_add, opt_lev_dim)
  END SUBROUTINE read_cdi_3d_real


  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 3D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_3d_lbc(parameters, varname, nlevs, var_out, npoints, idx_list)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),  INTENT(IN) :: varname
    INTEGER,           INTENT(IN) :: nlevs, npoints, idx_list(:)
    REAL(wp),       INTENT(INOUT) :: var_out(:,:,:) !< output field


    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_lbc'
    INTEGER :: jk, ierrstat, nmiss, i
    INTEGER :: vlistId, varId, zaxisId, gridId, tile_index, grid_size
    REAL(sp), ALLOCATABLE :: tmp_buf(:), map_buf(:) ! temporary local array
    LOGICAL :: lmap_buf, is_workroot
    CHARACTER(len=DICT_MAX_STRLEN) :: mapped_name
    INTEGER :: tlen

    is_workroot = my_process_is_mpi_workroot()

    IF (ANY(aes_phy_config(:)%ljsb)) THEN
      mapped_name = varname
      IF (parameters%have_dict) &
      & mapped_name = parameters%dict%get(TRIM(varname), DEFAULT=varname)
      mapped_name = tolower(mapped_name)
      tlen = LEN_TRIM(mapped_name)
      varID = -1
      DO i = 1, SIZE(parameters%variableNames)
        IF (tolower(parameters%variableNames(i))==mapped_name(1:tlen)) THEN
            varID = i-1
        END IF
      END DO
    ELSE
      CALL parameters%findVarId(varname, trivial_tile_att%getTileinfo_grb2(), varID, tile_index)
    ENDIF
    lmap_buf = .FALSE.

    IF (is_workroot) THEN
      vlistId = streamInqVlist(parameters%streamId)
      ! sanity check of the variable dimensions
      zaxisId = vlistInqVarZaxis(vlistId, varId)
      IF (zaxisInqSize(zaxisId) /= nlevs) THEN
        CALL finish(routine, "Incompatible dimensions! level count = "//trim(int2string(zaxisInqSize(zaxisId)))//&
        &                    " (expected "//trim(int2string(nlevs))//")")
      END IF

      gridId = vlistInqVarGrid(vlistId, varId)
      grid_size = gridInqSize(gridId)
      IF (grid_size /= parameters%glb_arr_len) THEN
        IF (npoints == grid_size) THEN
          lmap_buf = .TRUE.
        ELSE
          CALL finish(routine, "Incompatible dimensions! Grid size = "//trim(int2string(grid_size))//&
          &                    " (expected "//trim(int2string(npoints))//")")
        ENDIF
      END IF
    END IF


    ! allocate a buffer for one vertical level
    IF (is_workroot) THEN
      IF (lmap_buf) THEN
        ALLOCATE(tmp_buf(npoints), map_buf(parameters%glb_arr_len), STAT=ierrstat)
        map_buf(:) = 0._sp
      ELSE
        ALLOCATE(tmp_buf(1),map_buf(parameters%glb_arr_len), STAT=ierrstat)
      ENDIF
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ELSE
        ALLOCATE(tmp_buf(1), map_buf(1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF



    !FIXME: This code is most likely latency bound, not throughput bound. Needs some asynchronicity to hide the latencies.
    DO jk=1,nlevs
      IF (is_workroot) THEN
        ! read record as 1D field
        IF (lmap_buf) THEN
          CALL timeStreamReadVarSliceF(parameters, varID, jk-1, tmp_buf(:), nmiss)
!$OMP PARALLEL DO PRIVATE(i)
          DO i = 1, npoints
            map_buf(idx_list(i)) = tmp_buf(i)
          ENDDO
!$OMP END PARALLEL DO
        ELSE
          CALL timeStreamReadVarSliceF(parameters, varID, jk-1, map_buf(:), nmiss)
        ENDIF
      END IF

      CALL parameters%distribution%distribute(map_buf(:), var_out(:, jk, :), .FALSE.)
 
    END DO ! jk=1,nlevs

    ! clean up
    DEALLOCATE(tmp_buf, map_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE read_cdi_3d_lbc

  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 2D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_2d_lbc(parameters, varname, var_out, npoints, idx_list)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),  INTENT(IN) :: varname
    INTEGER,           INTENT(IN) :: npoints, idx_list(:)
    REAL(wp),       INTENT(INOUT) :: var_out(:,:) !< output field


    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_lbc'
    INTEGER :: jk, ierrstat, nmiss, i
    INTEGER :: vlistId, varId, zaxisId, gridId, tile_index, grid_size
    REAL(sp), ALLOCATABLE :: tmp_buf(:), map_buf(:) ! temporary local array
    LOGICAL :: lmap_buf, is_workroot
    CHARACTER(len=DICT_MAX_STRLEN) :: mapped_name
    INTEGER :: tlen

    is_workroot = my_process_is_mpi_workroot()
    IF (ANY(aes_phy_config(:)%ljsb)) THEN !If using JSB, bypassing TERRA-init routines
      mapped_name = varname
      IF (parameters%have_dict) &
      & mapped_name = parameters%dict%get(TRIM(varname), DEFAULT=varname)
      mapped_name = tolower(mapped_name)
      tlen = LEN_TRIM(mapped_name)
      varID = -1
      DO i = 1, SIZE(parameters%variableNames)
        IF (tolower(parameters%variableNames(i))==mapped_name(1:tlen)) THEN
            varID = i-1
        END IF
      END DO
    ELSE
      CALL parameters%findVarId(varname, trivial_tile_att%getTileinfo_grb2(), varID, tile_index)
    ENDIF
    lmap_buf = .FALSE.

    IF (is_workroot) THEN
      vlistId = streamInqVlist(parameters%streamId)

      gridId = vlistInqVarGrid(vlistId, varId)
      grid_size = gridInqSize(gridId)
      IF (grid_size /= parameters%glb_arr_len) THEN
        IF (npoints == grid_size) THEN
          lmap_buf = .TRUE.
        ELSE
          CALL finish(routine, "Incompatible dimensions! Grid size = "//trim(int2string(grid_size))//&
          &                    " (expected "//trim(int2string(npoints))//")")
        ENDIF
      END IF
    END IF


    ! allocate a buffer for one vertical level
    IF (is_workroot) THEN
      IF (lmap_buf) THEN
        ALLOCATE(tmp_buf(npoints), map_buf(parameters%glb_arr_len), STAT=ierrstat)
        map_buf(:) = 0._sp
      ELSE
        ALLOCATE(tmp_buf(1),map_buf(parameters%glb_arr_len), STAT=ierrstat)
      ENDIF
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ELSE
        ALLOCATE(tmp_buf(1), map_buf(1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF



    !FIXME: This code is most likely latency bound, not throughput bound. Needs some asynchronicity to hide the latencies.
      IF (is_workroot) THEN
        ! read record as 1D field
        IF (lmap_buf) THEN
          CALL timeStreamReadVarSliceF(parameters, varID, 0, tmp_buf(:), nmiss)
!$OMP PARALLEL DO PRIVATE(i)
          DO i = 1, npoints
            map_buf(idx_list(i)) = tmp_buf(i)
          ENDDO
!$OMP END PARALLEL DO
        ELSE
          CALL timeStreamReadVarSliceF(parameters, varID, 0, map_buf(:), nmiss)
        ENDIF
      END IF

      CALL parameters%distribution%distribute(map_buf(:), var_out(:, :), .FALSE.)
 

    ! clean up
    DEALLOCATE(tmp_buf, map_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE read_cdi_2d_lbc

  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 2D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_2d_sp(parameters, varID, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID
    REAL(wp), INTENT(INOUT) :: var_out(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_sp'
    INTEGER :: nmiss, ierrstat
    REAL(sp), ALLOCATABLE :: tmp_buf(:)

    IF (my_process_is_mpi_workroot()) THEN
        ! read record as 1D field
        ALLOCATE(tmp_buf(parameters%glb_arr_len), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        CALL timeStreamReadVarSliceF(parameters, varID, 0, tmp_buf(:), nmiss)
    ELSE
        ALLOCATE(tmp_buf(1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF

    CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :), .false.)

    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE read_cdi_2d_sp

  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 2D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_2d_wp(parameters, varID, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID
    REAL(wp), INTENT(INOUT) :: var_out(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_wp'
    INTEGER :: nmiss, ierrstat
    REAL(wp), ALLOCATABLE :: tmp_buf(:)

    IF (my_process_is_mpi_workroot()) THEN
        ! read record as 1D field
        ALLOCATE(tmp_buf(parameters%glb_arr_len), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        CALL timeStreamReadVarSlice(parameters, varID, 0, tmp_buf(:), nmiss)
    ELSE
        ALLOCATE(tmp_buf(1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF

    CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :), .false.)

    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE read_cdi_2d_wp


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  SUBROUTINE read_cdi_2d_real_tiles (parameters, varname, var_out, tileinfo)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT) :: var_out(:,:)   !< output field
    TYPE(t_tileinfo_grb2),   INTENT(IN)    :: tileinfo

    ! local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_real_tiles'
    INTEGER       :: varId, vlistId, gridId, subtypeID, tile_index

    CALL parameters%findVarId(varname, tileinfo, varID, tile_index)

    IF ((tile_index < 0) .AND. (tileinfo%idx > 0)) THEN
      CALL finish(routine, "Requested tile not found!")
    END IF

    IF (my_process_is_mpi_workroot()) THEN
      ! set active tile index, if this is a tile-based variable
      vlistId   = streamInqVlist(parameters%streamId)
      subtypeID  = vlistInqVarSubtype(vlistID,varID)
      IF (tile_index > 0)  CALL subtypeDefActiveIndex(subtypeID, tile_index)
      !sanity check on the variable dimensions
      gridId    = vlistInqVarGrid(vlistId, varId)
      IF (gridInqSize(gridId) /= parameters%glb_arr_len) CALL finish(routine, "Incompatible dimensions!"//&
      & " Grid size = "//trim(int2string(gridInqSize(gridId)))//" (expected "//trim(int2string(parameters%glb_arr_len))//")")
    END IF
    SELECT CASE(parameters%variableDatatype(varId+1))
        CASE(CDI_DATATYPE_PACK23:CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32)
            ! int32 is treated as double precision because single precision floats would cut off up to seven bits from the integer
            CALL read_cdi_2d_wp(parameters, varId, var_out)
        CASE DEFAULT
            ! XXX: Broadcasting CDI_DATATYPE_PACK1..CDI_DATATYPE_PACK22 DATA as single precision may actually change their values, but this
            !      error will always be smaller than the error made by storing the DATA as CDI_DATATYPE_PACK1..CDI_DATATYPE_PACK22 IN the
            !      first place.
            CALL read_cdi_2d_sp(parameters, varId, var_out)
    END SELECT

    IF(my_process_is_mpi_workroot()) THEN
      ! reset tile index
      IF (tile_index > 0)  CALL subtypeDefActiveIndex(subtypeID, 0)
    END IF
  END SUBROUTINE read_cdi_2d_real_tiles


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  SUBROUTINE read_cdi_2d_real (parameters, varname, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT) :: var_out(:,:)   !< output field

    CALL read_cdi_2d_real_tiles (parameters, varname, var_out, trivial_tile_att%getTileinfo_grb2())
  END SUBROUTINE read_cdi_2d_real


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for INTEGER fields
  !

  SUBROUTINE read_cdi_2d_int_tiles(parameters, varname, var_out, tileinfo)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(INOUT) :: var_out(:,:)   !< output field
    TYPE(t_tileinfo_grb2),   INTENT(IN)    :: tileinfo

    ! local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_int'
    REAL(wp), ALLOCATABLE :: var_tmp(:,:)
    INTEGER               :: ierrstat

    ! allocate a temporary array:
    ALLOCATE(var_tmp(SIZE(var_out,1), SIZE(var_out,2)), STAT=ierrstat)
    ! Initialization is needed in order to avoid errors in subsequent conversion to INTEGER
    var_tmp(:,:) = 0._wp
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! read the field as a REAL-valued field:
    CALL read_cdi_2d_real_tiles (parameters, varname, var_tmp, tileinfo )
    ! perform number conversion
    var_out(:,:) = NINT(var_tmp)
    ! clean up
    DEALLOCATE(var_tmp, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE read_cdi_2d_int_tiles


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for INTEGER fields
  !

  SUBROUTINE read_cdi_2d_int(parameters, varname, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(INOUT) :: var_out(:,:)   !< output field

    CALL read_cdi_2d_int_tiles(parameters, varname, var_out, trivial_tile_att%getTileinfo_grb2())
  END SUBROUTINE read_cdi_2d_int


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  SUBROUTINE read_cdi_2d_time_tiles (parameters, ntime, varname, var_out, tileinfo)
    TYPE(t_inputParameters), INTENT(INOUT)  :: parameters
    INTEGER,                 INTENT(IN)     :: ntime          !< time levels of file
    CHARACTER(len=*),        INTENT(IN)     :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT)  :: var_out(:,:,:) !< output field
    TYPE(t_tileinfo_grb2),   INTENT(IN)     :: tileinfo

    ! local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_time_tiles'
    INTEGER :: jt, nrecs, varId, subtypeID, tile_index, vlistID
    LOGICAL :: is_workroot

    is_workroot = my_process_is_mpi_workroot()

    ! Get var ID
    CALL parameters%findVarId(varname, tileinfo, varID, tile_index)

    IF ((tile_index < 0) .AND. (tileinfo%idx > 0)) THEN
      CALL finish(routine, "Requested tile not found!")
    END IF

    DO jt = 1, ntime
      IF (is_workroot) THEN
        ! set active tile index, if this is a tile-based variable
        vlistId    = streamInqVlist(parameters%streamId)
        subtypeID  = vlistInqVarSubtype(vlistID, varID)
        IF (tile_index > 0)  CALL subtypeDefActiveIndex(subtypeID, tile_index)
        nrecs = streamInqTimestep(parameters%streamId, (jt-1))
      END IF
      SELECT CASE(parameters%variableDatatype(varId+1))
        CASE(CDI_DATATYPE_PACK23:CDI_DATATYPE_PACK32, CDI_DATATYPE_FLT64, CDI_DATATYPE_INT32)
            ! int32 is treated as double precision because single precision floats would cut off up to seven bits from the integer
            CALL read_cdi_2d_wp(parameters, varId, var_out(:,:,jt))
        CASE DEFAULT
            ! XXX: Broadcasting CDI_DATATYPE_PACK1..CDI_DATATYPE_PACK22 DATA as single precision may actually change their values, but this
            !      error will always be smaller than the error made by storing the DATA as CDI_DATATYPE_PACK1..CDI_DATATYPE_PACK22 IN the
            !      first place.
            CALL read_cdi_2d_sp(parameters, varId, var_out(:,:,jt))
      END SELECT
    END DO
    IF (is_workroot) THEN
      ! reset tile index
      IF (tile_index > 0)  CALL subtypeDefActiveIndex(subtypeID, 0)
    END IF
  END SUBROUTINE read_cdi_2d_time_tiles


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  SUBROUTINE read_cdi_2d_time (parameters, ntime, varname, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER,                 INTENT(IN)    :: ntime          !< time levels of file
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT) :: var_out(:,:,:) !< output field

    CALL read_cdi_2d_time_tiles (parameters, ntime, varname, var_out, trivial_tile_att%getTileinfo_grb2())
  END SUBROUTINE read_cdi_2d_time

  SUBROUTINE cdiGetStringError(errorId, outErrorString)
    INTEGER(C_INT), VALUE :: errorId
    CHARACTER(KIND = C_CHAR, LEN=*), INTENT(INOUT) :: outErrorString
    CHARACTER(KIND = C_CHAR), dimension(:), POINTER :: cString
    INTEGER :: i, n

    cString => cdiStringError(errorId)
    outErrorString = ""
    IF(ASSOCIATED(cString)) THEN
      n = MIN(LEN(outErrorString), SIZE(cString, 1))
      DO i = 1, n
        outErrorString(i:i) = cString(i)
      END DO
    END IF
  END SUBROUTINE cdiGetStringError


  !------------------------------------------------------------------------------------------------
  !> Define new CDI variable based on given meta-data.
  !
  FUNCTION create_cdi_variable(vlistID, gridID, zaxisID,    &
    &                          info, missval, output_type,  &
    &                          gribout_config, i_lctype,    &
    &                          out_varnames_dict) RESULT(varID)
    INTEGER :: varID
    INTEGER,                INTENT(IN)         :: vlistID, gridID, zaxisID
    TYPE(t_var_metadata),   INTENT(IN), TARGET :: info
    REAL(dp),               INTENT(IN)         :: missval
    INTEGER,                INTENT(IN)         :: output_type
    TYPE(t_gribout_config), INTENT(IN)         :: gribout_config
    INTEGER,                INTENT(IN)         :: i_lctype
    TYPE(t_dictionary),     INTENT(IN)         :: out_varnames_dict
    ! local variables
    CHARACTER(LEN=DICT_MAX_STRLEN) :: mapped_name
    TYPE(t_cf_var), POINTER        :: this_cf
    INTEGER                        :: i

    ! Search name mapping for name in NetCDF file
    IF (info%cf%short_name /= '') THEN
      mapped_name = out_varnames_dict%get(info%cf%short_name, default=info%cf%short_name)
    ELSE
      mapped_name = out_varnames_dict%get(info%name, default=info%name)
    END IF


    varID = vlistDefVar(vlistID, gridID, zaxisID, MERGE(TIME_CONSTANT, TIME_VARYING, info%isteptype == TSTEP_CONSTANT))
    CALL vlistDefVarTsteptype(vlistID, varID, info%isteptype)

    CALL vlistDefVarName(vlistID, varID, TRIM(mapped_name))
    IF (info%post_op%lnew_cf) THEN
      this_cf => info%post_op%new_cf
    ELSE
      this_cf => info%cf
    END IF

    IF (this_cf%long_name /= '')     CALL vlistDefVarLongname(vlistID, varID, TRIM(this_cf%long_name))
    IF (this_cf%standard_name /= '') CALL vlistDefVarStdname(vlistID, varID, TRIM(this_cf%standard_name))
    IF (this_cf%units /= '')         CALL vlistDefVarUnits(vlistID, varID, TRIM(this_cf%units))

    IF (info%lmiss .OR.                                                  &
      &   (info%lmask_boundary .AND. ANY(config_lmask_boundary(:)) .AND. &
      &   (info%hgrid == GRID_UNSTRUCTURED_CELL))) THEN
      CALL vlistDefVarMissval(vlistID, varID, missval)
    END IF

    ! Set GRIB2 Triplet
    IF (info%post_op%lnew_grib2) THEN
      CALL vlistDefVarParam(vlistID, varID,                   &
        &  cdiEncodeParam(info%post_op%new_grib2%number,      &
        &                 info%post_op%new_grib2%category,    &
        &                 info%post_op%new_grib2%discipline) )
      IF ( output_type == FILETYPE_GRB2 ) THEN
        CALL vlistDefVarDatatype(vlistID, varID, info%post_op%new_grib2%bits)
      END IF
    ELSE
      CALL vlistDefVarParam(vlistID, varID,                   &
        &  cdiEncodeParam(info%grib2%number,                  &
        &                 info%grib2%category,                &
        &                 info%grib2%discipline) )
      IF ( output_type == FILETYPE_GRB2 ) THEN
        CALL vlistDefVarDatatype(vlistID, varID, info%grib2%bits)
      END IF
    END IF

    IF (output_type == FILETYPE_GRB2) THEN
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          ATTENTION                    !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Note that re-setting of the surface types must come AFTER (re)setting
      ! "productDefinitionTemplateNumber" in set_additional_GRIB2_keys. It was observed
      ! (i.e. for Ensemble output), that the surface-type information is lost again, if
      ! these settings are performed prior to "productDefinitionTemplateNumber"

      ! GRIB2 Quick hack: Set additional GRIB2 keys
      CALL set_GRIB2_additional_keys(vlistID, varID, gribout_config)

      ! Set ensemble keys in SECTION 4 (if applicable)
      CALL set_GRIB2_ensemble_keys(vlistID, varID, gribout_config)

      ! Set synsat keys (if applicable)
      CALL set_GRIB2_synsat_keys(vlistID, varID, info)

      ! Set keys for atmospheric chemical constituents, if applicable
      CALL set_GRIB2_chem_keys(vlistID, varID, info)

      ! Set local use SECTION 2
      CALL set_GRIB2_local_keys(vlistID, varID, gribout_config)

#ifndef __NO_ICON_ATMO__
      ! Set tile-specific GRIB2 keys (if applicable)
      CALL set_GRIB2_tile_keys(vlistID, varID, info, i_lctype, gribout_config%grib2_template_tile)
#endif

      ! Set further additional integer keys
      DO i=1,info%grib2%additional_keys%nint_keys
        CALL vlistDefVarIntKey(vlistID, varID, TRIM(info%grib2%additional_keys%int_key(i)%key), &
          &                    info%grib2%additional_keys%int_key(i)%val)
      END DO

      ! Set further additional double keys
      DO i=1,info%grib2%additional_keys%ndbl_keys
        CALL vlistDefVarDblKey(vlistID, varID, TRIM(info%grib2%additional_keys%dbl_key(i)%key), &
          &                    info%grib2%additional_keys%dbl_key(i)%val)
      END DO

      ! Adjustments for compensating non-backward-compatible
      ! metadata changes after version updates of the GRIB library.
      ! (For now, any adjustment presupposes additional integer GRIB keys being set.)
      IF ((gribout_config%grib_lib_compat == GRIB_LIB_COMPAT_ECC_2_31_0) &
        & .AND. (info%grib2%additional_keys%nint_keys > 0)) THEN
        CALL grib_lib_compatibility(grib_lib_compat      = gribout_config%grib_lib_compat, & ! in
          &                         vlistID              = vlistID,                        & ! in
          &                         varID                = varID,                          & ! in
          &                         zaxisID              = zaxisID,                        & ! in
          &                         additional_grib_keys = info%grib2%additional_keys      ) ! in
      ENDIF

    ELSE ! NetCDF
      CALL vlistDefVarDatatype(vlistID, varID, this_cf%datatype)
    ENDIF
   
  END FUNCTION create_cdi_variable


END MODULE mo_util_cdi
