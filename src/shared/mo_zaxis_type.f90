! Global registry for vertical axis types.
!
! Please note: The purpose of this module is to register vertical axis types
! (ie. categories), which can be used by the "add_var" mechanism in ICON.
!
! The meta-data for these axis types is *not* defined in this place, but in
! the module "src/io/shared/mo_name_list_output_zaxes.f90" .
! The reason for this separation is that the meta-data
! is only required for writing variables to disk.
!
! To make a long story short: After adding a new vertical axis type here,
! you probably need to add the meta-data definition for this type in
! "mo_name_list_output_zaxes.f90", otherwise you won't be able
! to read/write the variable.
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

MODULE mo_zaxis_type

  USE ISO_C_BINDING,    ONLY: C_INT32_T
  USE mo_exception,     ONLY: finish
  USE mo_hash_table,    ONLY: t_HashTable, hashTable_make
  USE mo_cdi, ONLY: ZAXIS_ALTITUDE, ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP,        &
    &               ZAXIS_DEPTH_BELOW_LAND, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_GENERIC, ZAXIS_HEIGHT, &
    &               ZAXIS_HYBRID, ZAXIS_HYBRID_HALF, ZAXIS_ISENTROPIC, ZAXIS_ISOTHERM_ZERO,     &
    &               ZAXIS_LAKE_BOTTOM, ZAXIS_MEANSEA, ZAXIS_MIX_LAYER, ZAXIS_PRESSURE,          &
    &               ZAXIS_REFERENCE, ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_SURFACE, ZAXIS_TOA,        &
    &               CDI_UNDEFID, zaxisCreate, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE,      &
    &               cdiDefKeyInt



  IMPLICIT NONE

  PUBLIC :: t_zaxisType
  PUBLIC :: t_zaxisTypeList
  PUBLIC :: zaxisTypeList


  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_zaxis_type'


  ! Convenience variables for those ICON-internal zaxis types which
  ! are known a priori. The actual values of these variables is
  ! irrelevant - they are set by the "register" routine of
  ! t_zaxisTypeList.
  !
  ! See "t_zaxisTypeList" for details.
  !
  INTEGER, PUBLIC :: &
    &   ZA_SURFACE, ZA_CLOUD_BASE, ZA_CLOUD_TOP , ZA_ISOTHERM_ZERO, ZA_REFERENCE,       &
    &   ZA_REFERENCE_HALF, ZA_REFERENCE_HALF_HHL, ZA_HYBRID, ZA_HYBRID_HALF,            &
    &   ZA_HYBRID_HALF_HHL, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, ZA_SNOW,       &
    &   ZA_SNOW_HALF, ZA_PRESSURE, ZA_HEIGHT, ZA_HEIGHT_2M, ZA_HEIGHT_10M,              &
    &   ZA_ALTITUDE, ZA_MEANSEA, ZA_ISENTROPIC, ZA_TOA, ZA_PRESSURE_800,                &
    &   ZA_PRESSURE_400, ZA_PRESSURE_0, ZA_DEPTH_RUNOFF_S, ZA_DEPTH_RUNOFF_G,           &
    &   ZA_LAKE_BOTTOM, ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER ,                             &
    &   ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF,        &
    &   ZA_GENERIC_ICE, ZA_OCEAN_SEDIMENT, ZA_PRES_FL_BOT_TOP, ZA_ATMOSPHERE,           &
    &   ZA_HEIGHT_2M_LAYER, ZA_ECHOTOP, ZA_TROPOPAUSE, ZA_wshear,                       &
    &   ZA_PRESSURE_LAPSERATE, ZA_SRH

  !> Derived type holding a the ICON-internal key for a single
  !  vertical axis type. See "t_zaxisTypeList" for details.
  !
  TYPE :: t_zaxisKey
    INTEGER :: icon_zaxis_type
  END TYPE t_zaxisKey


  !> Derived type holding a single vertical axis type. See
  !  "t_zaxisTypeList" for details.
  !
  TYPE :: t_zaxisType
    INTEGER :: icon_zaxis_type
    INTEGER :: cdi_zaxis_type
    LOGICAL :: is_2d
  CONTAINS
    GENERIC, PUBLIC :: OPERATOR(==) => eqv
    PROCEDURE :: eqv => t_zaxisType_eqv
  END TYPE t_zaxisType


  !> This derived data type manages a list of available vertical axis
  !  types.  Every ICON variable, regardless whether it is intended
  !  for output or not, requires a vertical axis type (even 2D fields
  !  require a "trivial" axis).
  !
  !  In ICON "add_var"/"add_ref" makes use of these parameters instead
  !  of the ZAXIS types defined in cdi.inc. This makes it much easier
  !  to distinguish vertical axes of the same type
  !  (e.g. ZAXIS_HEIGHT), which may have different level heights or
  !  numbers.
  !
  !  Note that this derived data type is only a kind of registry. If
  !  all axis types were known a priori, then it could be replaced by
  !  a long list of integer constants. However, there exist axis types
  !  contributed by JSBACH, e.g., which make a dynamic list necessary.
  !
  !  Meta-data attached to the vertical axis is not handled by this
  !  derived data type. See "t_verticalAxis" for this task.
  !
  TYPE t_zaxisTypeList
    TYPE(t_HashTable) :: list
    INTEGER           :: max_icon_zaxis_type

  CONTAINS
    PROCEDURE :: finalize       => t_zaxisTypeList_finalize        !< destructor
    PROCEDURE :: is_2d          => t_zaxisTypeList_is_2d           !< return .TRUE. if this is a 2D field axis
    PROCEDURE :: register       => t_zaxisTypeList_register        !< register a new vertical axis
    PROCEDURE :: getEntry       => t_zaxisTypeList_getEntry        !< retrieve an entry from hash table
    PROCEDURE :: cdi_zaxis_type => t_zaxisTypeList_cdi_zaxis_type  !< returns associated CDI zaxis type
    PROCEDURE :: za_count       => t_zaxisTypeList_za_count        !< return max. axis key
  END TYPE t_zaxisTypeList

  !> Constructor. Registers all zaxis types known a priori.
  INTERFACE t_zaxisTypeList
    MODULE PROCEDURE new_zaxisTypeList
  END INTERFACE


  ! -----------------------------------------------------------------
  !> Global registry for zaxis types
  ! -----------------------------------------------------------------

  TYPE(t_zaxisTypeList) :: zaxisTypeList

  LOGICAL :: l_add_addition_zaxis_types = .FALSE.
  INTEGER, PUBLIC :: ZAXIS_TROPOPAUSE = CDI_UNDEFID

CONTAINS

  !> Comparison operator for two zaxis types.
  !
  LOGICAL FUNCTION t_zaxisType_eqv(zaxis1, zaxis2)
    CLASS(t_zaxisType), INTENT(IN) :: zaxis1, zaxis2

    t_zaxisType_eqv = .TRUE.
    IF (zaxis1%icon_zaxis_type /=     zaxis2%icon_zaxis_type)  t_zaxisType_eqv=.FALSE.
    IF (zaxis1%cdi_zaxis_type  /=     zaxis2%cdi_zaxis_type)   t_zaxisType_eqv=.FALSE.
    IF (zaxis1%is_2d           .NEQV. zaxis2%is_2d)            t_zaxisType_eqv=.FALSE.
  END FUNCTION t_zaxisType_eqv


  !> Constructor for t_zaxisTypeList. Registers all zaxis types known
  !  a priori.
  !
  FUNCTION new_zaxisTypeList()  RESULT(za_list)
    TYPE(t_zaxisTypeList) :: za_list

    INTEGER :: res

    IF (l_add_addition_zaxis_types) THEN
      ZAXIS_TROPOPAUSE = zaxisCreate(ZAXIS_GENERIC, 1)
      ! WMO GRIB2 code table 4.5: level type TOPOPAUSE = 7
      res = cdiDefKeyInt(ZAXIS_TROPOPAUSE, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, 7)
      l_add_addition_zaxis_types = .FALSE.
    ENDIF

    za_list%max_icon_zaxis_type = 0
    za_list%list = hashTable_make(list_hashKey, list_equalKeys)

    ZA_SURFACE                 = za_list%register(cdi_zaxis_type=ZAXIS_SURFACE            , is_2D=.TRUE.)
    ! Atmosphere
    ZA_CLOUD_BASE              = za_list%register(cdi_zaxis_type=ZAXIS_CLOUD_BASE         , is_2D=.TRUE.)
    ZA_CLOUD_TOP               = za_list%register(cdi_zaxis_type=ZAXIS_CLOUD_TOP          , is_2D=.TRUE.)
    ZA_ISOTHERM_ZERO           = za_list%register(cdi_zaxis_type=ZAXIS_ISOTHERM_ZERO      , is_2D=.TRUE.)
    ZA_REFERENCE               = za_list%register(cdi_zaxis_type=ZAXIS_REFERENCE          , is_2D=.FALSE.)
    ZA_REFERENCE_HALF          = za_list%register(cdi_zaxis_type=ZAXIS_REFERENCE          , is_2D=.FALSE.)
    ZA_REFERENCE_HALF_HHL      = za_list%register(cdi_zaxis_type=ZAXIS_REFERENCE          , is_2D=.FALSE.)
    ZA_ECHOTOP                 = za_list%register(cdi_zaxis_type=ZAXIS_GENERIC            , is_2D=.FALSE.)

    !DR *********** FIXME *************
    ! Re-set
    ! ZA_HYBRID       -> ZA_REFERENCE
    ! ZA_HYBRID_HALF  -> ZA_REFERENCE_HALF
    ! as long as ZA_hybrid/ZA_hybrid_half is used throughout the code.
    ! Should be replaced by ZA_reference/ZA_reference_half for the
    ! nonhydrostatic model.

    ! ZA_HYBRID                  = za_list%register(cdi_zaxis_type=ZAXIS_HYBRID             , is_2D=.FALSE.)
    ! ZA_HYBRID_HALF             = za_list%register(cdi_zaxis_type=ZAXIS_HYBRID_HALF        , is_2D=.FALSE.)
    ! ZA_HYBRID_HALF_HHL         = za_list%register(cdi_zaxis_type=ZAXIS_HYBRID_HALF        , is_2D=.FALSE.)
    ZA_HYBRID                  = ZA_REFERENCE
    ZA_HYBRID_HALF             = ZA_REFERENCE_HALF
    ZA_HYBRID_HALF_HHL         = ZA_REFERENCE_HALF_HHL
    !DR*********WILL BE REMOVED SOON**********

    ZA_DEPTH_BELOW_LAND        = za_list%register(cdi_zaxis_type=ZAXIS_DEPTH_BELOW_LAND   , is_2D=.FALSE.)
    ZA_DEPTH_BELOW_LAND_P1     = za_list%register(cdi_zaxis_type=ZAXIS_DEPTH_BELOW_LAND   , is_2D=.FALSE.)
    ZA_SNOW                    = za_list%register(cdi_zaxis_type=ZAXIS_GENERIC            , is_2D=.FALSE.)
    ZA_SNOW_HALF               = za_list%register(cdi_zaxis_type=ZAXIS_GENERIC            , is_2D=.FALSE.)
    ZA_PRESSURE                = za_list%register(cdi_zaxis_type=ZAXIS_PRESSURE           , is_2D=.FALSE.)
    ZA_PRESSURE_LAPSERATE      = za_list%register(cdi_zaxis_type=ZAXIS_PRESSURE           , is_2D=.TRUE.)
    ZA_HEIGHT                  = za_list%register(cdi_zaxis_type=ZAXIS_HEIGHT             , is_2D=.FALSE.)
    ZA_WSHEAR                  = za_list%register(cdi_zaxis_type=ZAXIS_HEIGHT             , is_2D=.FALSE.)
    ZA_SRH                     = za_list%register(cdi_zaxis_type=ZAXIS_HEIGHT             , is_2D=.FALSE.)
    ZA_HEIGHT_2M               = za_list%register(cdi_zaxis_type=ZAXIS_HEIGHT             , is_2D=.TRUE.)
    ZA_HEIGHT_10M              = za_list%register(cdi_zaxis_type=ZAXIS_HEIGHT             , is_2D=.TRUE.)
    ZA_ALTITUDE                = za_list%register(cdi_zaxis_type=ZAXIS_ALTITUDE           , is_2D=.FALSE.)
    ZA_MEANSEA                 = za_list%register(cdi_zaxis_type=ZAXIS_MEANSEA            , is_2D=.TRUE.)
    ZA_ISENTROPIC              = za_list%register(cdi_zaxis_type=ZAXIS_ISENTROPIC         , is_2D=.FALSE.)
    ZA_TOA                     = za_list%register(cdi_zaxis_type=ZAXIS_TOA                , is_2D=.TRUE.)
    ZA_PRESSURE_800            = za_list%register(cdi_zaxis_type=ZAXIS_PRESSURE           , is_2D=.TRUE.)
    ZA_PRESSURE_400            = za_list%register(cdi_zaxis_type=ZAXIS_PRESSURE           , is_2D=.TRUE.)
    ZA_PRESSURE_0              = za_list%register(cdi_zaxis_type=ZAXIS_PRESSURE           , is_2D=.TRUE.)
    ZA_DEPTH_RUNOFF_S          = za_list%register(cdi_zaxis_type=ZAXIS_DEPTH_BELOW_LAND   , is_2D=.TRUE.)
    ZA_DEPTH_RUNOFF_G          = za_list%register(cdi_zaxis_type=ZAXIS_DEPTH_BELOW_LAND   , is_2D=.TRUE.)
    ! Lake
    ZA_LAKE_BOTTOM             = za_list%register(cdi_zaxis_type=ZAXIS_LAKE_BOTTOM        , is_2D=.TRUE.)
    ZA_LAKE_BOTTOM_HALF        = za_list%register(cdi_zaxis_type=ZAXIS_LAKE_BOTTOM        , is_2D=.TRUE.)
    ZA_MIX_LAYER               = za_list%register(cdi_zaxis_type=ZAXIS_MIX_LAYER          , is_2D=.TRUE.)
    ZA_SEDIMENT_BOTTOM_TW_HALF = za_list%register(cdi_zaxis_type=ZAXIS_SEDIMENT_BOTTOM_TW , is_2D=.TRUE.)
    ! Ocean
    ZA_DEPTH_BELOW_SEA         = za_list%register(cdi_zaxis_type=ZAXIS_DEPTH_BELOW_SEA    , is_2D=.FALSE.)
    ZA_DEPTH_BELOW_SEA_HALF    = za_list%register(cdi_zaxis_type=ZAXIS_DEPTH_BELOW_SEA    , is_2D=.FALSE.)
    ZA_GENERIC_ICE             = za_list%register(cdi_zaxis_type=ZAXIS_GENERIC            , is_2D=.FALSE.)
    ! HAMOCC sediment
    ZA_OCEAN_SEDIMENT          = za_list%register(cdi_zaxis_type=ZAXIS_GENERIC            , is_2D=.FALSE.)
    ! needed for ICON-ART
    ZA_PRES_FL_BOT_TOP         = za_list%register(cdi_zaxis_type=ZAXIS_PRESSURE           , is_2D=.TRUE.)
    ZA_ATMOSPHERE              = za_list%register(cdi_zaxis_type=ZAXIS_ATMOSPHERE         , is_2D=.TRUE.)
    ZA_HEIGHT_2M_LAYER         = za_list%register(cdi_zaxis_type=ZAXIS_HEIGHT             , is_2D=.TRUE.)
    ZA_TROPOPAUSE              = za_list%register(cdi_zaxis_type=ZAXIS_GENERIC            , is_2D=.TRUE.)    

  END FUNCTION new_zaxisTypeList


  SUBROUTINE t_zaxisTypeList_finalize(zaxisTypeList)
    CLASS(t_zaxisTypeList), INTENT(INOUT) :: zaxisTypeList
    CALL zaxisTypeList%list%destruct()
  END SUBROUTINE t_zaxisTypeList_finalize


  !> Auxiliary function for the internal hash table: compute hash key.
  !
  INTEGER(C_INT32_T) FUNCTION list_hashKey(key) RESULT(RESULT)
    CLASS(*), POINTER, INTENT(in) :: key
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":list_hashKey"

    SELECT TYPE(key)
    TYPE IS(t_zaxisKey)
      RESULT = key%icon_zaxis_type
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for key.")
    END SELECT
  END FUNCTION list_hashKey


  !> Auxiliary function for the internal hash table: compare keys.
  !
  LOGICAL FUNCTION list_equalKeys(keyA, keyB) RESULT(RESULT)
    CLASS(*), POINTER, INTENT(in) :: keyA, keyB
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":list_equalKeys"

    SELECT TYPE(keyA)
    TYPE IS(t_zaxisKey)
      SELECT TYPE(keyB)
      TYPE IS(t_zaxisKey)
        RESULT = (keyA%icon_zaxis_type == keyB%icon_zaxis_type)
        CLASS DEFAULT
        CALL finish(routine, "Unknown type for keyB.")
      END SELECT
      CLASS DEFAULT
      CALL finish(routine, "Unknown type for keyA.")
    END SELECT
  END FUNCTION list_equalKeys


  !> Add a new zaxisType entry to the list of zaxis types.
  !
  FUNCTION t_zaxisTypeList_register(zaxisTypeList, cdi_zaxis_type, is_2d, icon_zaxis_type)  RESULT(icon_type)

    INTEGER :: icon_type
    CLASS(t_zaxisTypeList), INTENT(INOUT) :: zaxisTypeList
    INTEGER,                INTENT(IN)    :: cdi_zaxis_type
    LOGICAL,                INTENT(IN)    :: is_2d
    INTEGER, OPTIONAL,      INTENT(IN)    :: icon_zaxis_type

    CHARACTER(LEN=*), PARAMETER     :: routine = modname//'::t_zaxisTypeList_register'

    CLASS(*),POINTER :: p_key
    CLASS(*),POINTER :: p_value

    IF (PRESENT(icon_zaxis_type)) THEN
      ! check if icon_zaxis_type already set
      ALLOCATE(t_zaxisKey :: p_key)
      SELECT TYPE(p_key)
      TYPE IS(t_zaxisKey)
        p_key%icon_zaxis_type = icon_zaxis_type
      END SELECT
      p_value => zaxisTypeList%list%getEntry(p_key)
      IF (ASSOCIATED(p_value)) THEN
        CALL finish(routine, "Key already registered.")        
      ELSE
        icon_type = icon_zaxis_type
      END IF
      DEALLOCATE(p_key)
    ELSE
      icon_type = zaxisTypeList%max_icon_zaxis_type + 1
    END IF

    ALLOCATE(t_zaxisKey  :: p_key  )
    ALLOCATE(t_zaxisType :: p_value)
    SELECT TYPE(p_key)
    TYPE IS(t_zaxisKey)
      p_key%icon_zaxis_type = icon_type
    END SELECT
    SELECT TYPE(p_value)
    TYPE IS(t_zaxisType)
      p_value%icon_zaxis_type = icon_type
      p_value%cdi_zaxis_type  = cdi_zaxis_type
      p_value%is_2d           = is_2d
    END SELECT
    CALL zaxisTypeList%list%setEntry(p_key, p_value)
    zaxisTypeList%max_icon_zaxis_type = MAX(zaxisTypeList%max_icon_zaxis_type, icon_type)
  END FUNCTION t_zaxisTypeList_register


  !> Retrieve a zaxisType entry from the list of zaxis types.
  !
  FUNCTION t_zaxisTypeList_getEntry(zaxisTypeList, icon_zaxis_type, ierror)  RESULT(zaxisType)
    TYPE(t_zaxisType) :: zaxisType
    CLASS(t_zaxisTypeList), INTENT(IN)  :: zaxisTypeList
    INTEGER,                INTENT(IN)  :: icon_zaxis_type
    INTEGER, OPTIONAL,      INTENT(OUT) :: ierror

    CHARACTER(LEN=*), PARAMETER        :: routine = modname//'::t_zaxisTypeList_getEntry'
    CLASS(*), POINTER     :: p_key, p_value

    ! retrieve entry from hash table:
    IF (PRESENT(ierror))  ierror = 0
    ALLOCATE(t_zaxisKey :: p_key)
    SELECT TYPE(p_key)
    TYPE IS(t_zaxisKey)
      p_key%icon_zaxis_type = icon_zaxis_type
    END SELECT
    p_value => zaxisTypeList%list%getEntry(p_key)
    IF (ASSOCIATED(p_value)) THEN
      SELECT TYPE (p_value)
      TYPE IS(t_zaxisType)
        zaxisType = p_value
      CLASS DEFAULT
      IF (PRESENT(ierror)) THEN
        ierror = 1
      ELSE
        CALL finish(routine, "Wrong return type.")
      END IF
      END SELECT
    ELSE
      IF (PRESENT(ierror)) THEN
        ierror = 1
      ELSE
        CALL finish(routine, "Index not found!")
      END IF
    ENDIF
    DEALLOCATE(p_key)
  END FUNCTION t_zaxisTypeList_getEntry


  !> Return whether the zaxis "icon_zaxis_type" is a 2D field.
  !
  FUNCTION t_zaxisTypeList_is_2d(zaxisTypeList, icon_zaxis_type)  RESULT(is_2d)
    LOGICAL                            :: is_2d
    CLASS(t_zaxisTypeList), INTENT(IN) :: zaxisTypeList
    INTEGER,                INTENT(IN) :: icon_zaxis_type
    TYPE(t_zaxisType) :: zaxisType
    zaxisType = zaxisTypeList%getEntry(icon_zaxis_type)
    is_2d = zaxisType%is_2d
  END FUNCTION t_zaxisTypeList_is_2d


  !> Return the CDI axis type associated with "icon_zaxis_type".
  !
  FUNCTION t_zaxisTypeList_cdi_zaxis_type(zaxisTypeList, icon_zaxis_type)  RESULT(cdi_zaxis_type)
    INTEGER                            :: cdi_zaxis_type
    CLASS(t_zaxisTypeList), INTENT(IN) :: zaxisTypeList
    INTEGER,                INTENT(IN) :: icon_zaxis_type
    TYPE(t_zaxisType) :: zaxisType
    INTEGER           :: ierror
    zaxisType = zaxisTypeList%getEntry(icon_zaxis_type, ierror)
    IF (ierror == 0) THEN
      cdi_zaxis_type = zaxisType%cdi_zaxis_type
    ELSE
      cdi_zaxis_type = CDI_UNDEFID
    END IF
  END FUNCTION t_zaxisTypeList_cdi_zaxis_type


  !> Return the largest ICON zaxis ID (sometimes needed for array allocation).
  !
  FUNCTION t_zaxisTypeList_za_count(zaxisTypeList)  RESULT(za_count)
    INTEGER :: za_count
    CLASS(t_zaxisTypeList), INTENT(IN) :: zaxisTypeList
    za_count = zaxisTypeList%max_icon_zaxis_type
  END FUNCTION t_zaxisTypeList_za_count

END MODULE mo_zaxis_type
