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

! Module handling the meta-data types for vertical axes for the output
! module.
!
! See module mo_name_list_output_zaxes for details of the
! implementation.

MODULE mo_name_list_output_zaxes_types

  USE ISO_C_BINDING,                        ONLY: C_SIGNED_CHAR
  USE mo_cdi,                               ONLY: CDI_UNDEFID, zaxisCreate, zaxisDefNumber, zaxisDefUUID,        &
    &                                             zaxisDefLevels, zaxisDefLbounds, zaxisDefUbounds, zaxisDefVct, &
    &                                             zaxisDefUnits, zaxisDefNlevRef, zaxisDefName, zaxisDestroy,    &
    &                                             zaxisDefLongname, cdiDefKeyInt, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &
    &                                             CDI_GLOBAL
  USE mo_zaxis_type,                        ONLY: t_zaxisType
  USE mo_kind,                              ONLY: wp, dp
  USE mo_exception,                         ONLY: finish
  USE mo_packed_message,                    ONLY: t_PackedMessage
  USE mo_mpi,                               ONLY: p_get_bcast_role

  IMPLICIT NONE

  PRIVATE

  ! subroutines
  PUBLIC :: t_verticalAxis, t_verticalAxisList


  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_zaxes_types'

  !> axis ID which signals the end of an MPI send/receive process
  INTEGER, PARAMETER :: MSG_SENTINEL = -9999


  !> Derived type which contains all meta-data of a vertical axis.
  !
  TYPE t_verticalAxis

    TYPE(t_zaxisType)                    :: zaxisType       !< CDI axis type        
    INTEGER                              :: zaxisNlev       !< no. of levels

    INTEGER                              :: cdi_id = CDI_UNDEFID       !< CDI object ID (or undefined)

    REAL(dp),                    ALLOCATABLE :: zaxisLevels(:)  !< level list
    REAL(dp),                    ALLOCATABLE :: zaxisLbounds(:) !< lower bounds
    REAL(dp),                    ALLOCATABLE :: zaxisUbounds(:) !< upper bounds

    CHARACTER(LEN=:),            ALLOCATABLE :: zaxisName       !< name string
    CHARACTER(LEN=:),            ALLOCATABLE :: zaxisLongname   !< extended name string
    CHARACTER(LEN=:),            ALLOCATABLE :: zaxisUnits      !< axis units

    INTEGER,                     ALLOCATABLE :: zaxisNumber     !< reference no. for unstructured grid
    INTEGER,                     ALLOCATABLE :: zaxisNlevRef    !< no. of half levels of generalized Z-axis

    INTEGER(KIND=C_SIGNED_CHAR), ALLOCATABLE :: zaxisUUID(:)    !< UUID of vertical grid
    REAL(dp),                    ALLOCATABLE :: zaxisVct(:)     !< vertical coordinate table

    INTEGER,                     ALLOCATABLE :: zaxisDefLtype   !< additional GRIB2 level type

  CONTAINS
    PROCEDURE :: finalize       => t_verticalAxis_finalize        !< destructor
    PROCEDURE :: cdiZaxisCreate => t_verticalAxis_cdiZaxisCreate  !< call to CDI library
    PROCEDURE :: set            => t_verticalAxis_set             !< set meta-data

    PROCEDURE :: packMPI        => t_verticalAxis_packMPI
    PROCEDURE :: unpackMPI      => t_verticalAxis_unpackMPI
    PROCEDURE :: sendMPI        => t_verticalAxis_sendMPI
    PROCEDURE :: recvMPI        => t_verticalAxis_recvMPI
    PROCEDURE :: bcastMPI       => t_verticalAxis_bcastMPI

    GENERIC, PUBLIC :: OPERATOR(==) => eqv
    PROCEDURE :: eqv            => t_verticalAxis_eqv
  END TYPE t_verticalAxis

  INTERFACE t_verticalAxis
    MODULE PROCEDURE new_verticalAxis
  END INTERFACE



  !> Linked list of vertical axis objects
  !
  TYPE t_verticalAxisList
    
    TYPE(t_verticalAxis),     POINTER :: axis => NULL()
    TYPE(t_verticalAxisList), POINTER :: next => NULL()

  CONTAINS
    PROCEDURE :: finalize       => t_verticalAxisList_finalize    !< destructor
    PROCEDURE :: getEntry       => t_verticalAxisList_getEntry    !< search for specific axis in list

    PROCEDURE :: sendMPI        => t_verticalAxisList_sendMPI
    PROCEDURE :: recvMPI        => t_verticalAxisList_recvMPI
    PROCEDURE :: bcastMPI       => t_verticalAxisList_bcastMPI

    GENERIC, PUBLIC :: append   => append_axis, append_list
    PROCEDURE :: append_axis    => t_verticalAxisList_append      !< append axis object to list
    PROCEDURE :: append_list    => t_verticalAxisList_append_list !< append axis list to this list
  END TYPE t_verticalAxisList

  INTERFACE t_verticalAxisList
    MODULE PROCEDURE new_verticalAxisList
  END INTERFACE


CONTAINS

  ! --------------------------------------------------------------------------------------
  !> Destructor.
  !
  SUBROUTINE t_verticalAxis_finalize(axis)
    CLASS(t_verticalAxis), INTENT(INOUT) :: axis

    axis%zaxisNlev = 0
    
    IF (ALLOCATED(axis%zaxisLevels))    DEALLOCATE(axis%zaxisLevels)
    IF (ALLOCATED(axis%zaxisLbounds))   DEALLOCATE(axis%zaxisLbounds)
    IF (ALLOCATED(axis%zaxisUbounds))   DEALLOCATE(axis%zaxisUbounds)
    IF (ALLOCATED(axis%zaxisName))      DEALLOCATE(axis%zaxisName)
    IF (ALLOCATED(axis%zaxisLongname))  DEALLOCATE(axis%zaxisLongname)
    IF (ALLOCATED(axis%zaxisUnits))     DEALLOCATE(axis%zaxisUnits)
    IF (ALLOCATED(axis%zaxisNumber))    DEALLOCATE(axis%zaxisNumber)
    IF (ALLOCATED(axis%zaxisNlevRef))   DEALLOCATE(axis%zaxisNlevRef)
    IF (ALLOCATED(axis%zaxisUUID))      DEALLOCATE(axis%zaxisUUID)
    IF (ALLOCATED(axis%zaxisVct))       DEALLOCATE(axis%zaxisVct)
    IF (ALLOCATED(axis%zaxisDefLtype))  DEALLOCATE(axis%zaxisDefLtype)

    IF (axis%cdi_id /= CDI_UNDEFID) THEN
      CALL zaxisDestroy(axis%cdi_id)
      axis%cdi_id = CDI_UNDEFID
    END IF
  END SUBROUTINE t_verticalAxis_finalize


  ! --------------------------------------------------------------------------------------
  !> Based on the object data, create a CDI axis and return the ID.
  !
  SUBROUTINE t_verticalAxis_cdiZaxisCreate(axis) 
    CLASS(t_verticalAxis), INTENT(INOUT) :: axis
    INTEGER :: cdiID
    INTEGER :: res
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::t_verticalAxis_cdiZaxisCreate'

    cdiID = zaxisCreate(axis%zaxisType%cdi_zaxis_type, axis%zaxisNlev)
    IF ((cdiID == CDI_UNDEFID) .OR. (axis%zaxisNlev <= 0)) THEN
      CALL finish(routine, "CDI zaxis creation failed!")
    END IF
    axis%cdi_id = cdiID

    IF (ALLOCATED(axis%zaxisLevels))   &
      CALL zaxisDefLevels(cdiID, axis%zaxisLevels)
    IF (ALLOCATED(axis%zaxisLbounds) ) &
      CALL zaxisDefLbounds(cdiID, axis%zaxisLbounds)
    IF (ALLOCATED(axis%zaxisUbounds))  &
      CALL zaxisDefUbounds(cdiID, axis%zaxisUbounds)
    IF (ALLOCATED(axis%zaxisName))     &
      CALL zaxisDefName(cdiID, TRIM(axis%zaxisName))
    IF (ALLOCATED(axis%zaxisLongname)) &
      CALL zaxisDefLongname(cdiID, TRIM(axis%zaxisLongname))
    IF (ALLOCATED(axis%zaxisUnits))    &
      CALL zaxisDefUnits(cdiID, TRIM(axis%zaxisUnits))
    IF (ALLOCATED(axis%zaxisNumber))   &
      CALL zaxisDefNumber(cdiID, axis%zaxisNumber)
    IF (ALLOCATED(axis%zaxisNlevRef))  &
      CALL zaxisDefNlevRef(cdiID, axis%zaxisNlevRef)
    IF (ALLOCATED(axis%zaxisUUID))     &
      CALL zaxisDefUUID (cdiID, axis%zaxisUUID)
    IF (ALLOCATED(axis%zaxisVct))      &
      CALL zaxisDefVct(cdiID, SIZE(axis%zaxisVct), axis%zaxisVct)
    IF (ALLOCATED(axis%zaxisDefLtype)) &
      res = cdiDefKeyInt(cdiID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, axis%zaxisDefLtype)

  END SUBROUTINE t_verticalAxis_cdiZaxisCreate


  ! --------------------------------------------------------------------------------------
  !> (Partially) sets meta-data of the vertical axis.
  !
  SUBROUTINE t_verticalAxis_set(axis, &
    &                           zaxisLevels,  zaxisLbounds, zaxisUbounds, zaxisName, zaxisLongname, &
    &                           zaxisUnits, zaxisNumber, zaxisNlevRef, zaxisUUID, zaxisVct, zaxisDefLtype)

    CLASS(t_verticalAxis),       INTENT(INOUT) :: axis
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisLevels(:)
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisLbounds(:)
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisUbounds(:)
    CHARACTER(LEN=*),            INTENT(IN), OPTIONAL :: zaxisName
    CHARACTER(LEN=*),            INTENT(IN), OPTIONAL :: zaxisLongname
    CHARACTER(LEN=*),            INTENT(IN), OPTIONAL :: zaxisUnits
    INTEGER,                     INTENT(IN), OPTIONAL :: zaxisNumber
    INTEGER,                     INTENT(IN), OPTIONAL :: zaxisNlevRef
    INTEGER(KIND=C_SIGNED_CHAR), INTENT(IN), OPTIONAL :: zaxisUUID(:)
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisVct(:)
    INTEGER,                     INTENT(IN), OPTIONAL :: zaxisDefLtype

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::t_verticalAxis_set'

    IF (PRESENT(zaxisLevels)) THEN
      IF (ALLOCATED(axis%zaxisLevels))   CALL finish(routine, "'zaxisLevels' already initialized!")
      ALLOCATE(axis%zaxisLevels(SIZE(zaxisLevels)))
      axis%zaxisLevels = zaxisLevels
    END IF
    IF (PRESENT(zaxisLbounds)) THEN
      IF (ALLOCATED(axis%zaxisLbounds))  CALL finish(routine, "'zaxisLbounds' already initialized!")
      ALLOCATE(axis%zaxisLbounds(SIZE(zaxisLbounds)))
      axis%zaxisLbounds = zaxisLbounds
    END IF
    IF (PRESENT(zaxisUbounds)) THEN
      IF (ALLOCATED(axis%zaxisUbounds))  CALL finish(routine, "'zaxisUbounds' already initialized!")
      ALLOCATE(axis%zaxisUbounds(SIZE(zaxisUbounds)))
      axis%zaxisUbounds = zaxisUbounds
    END IF
    IF (PRESENT(zaxisName)) THEN
      IF (ALLOCATED(axis%zaxisName))     CALL finish(routine, "'zaxisName' already initialized!")
      ALLOCATE(CHARACTER(LEN_TRIM(zaxisName)) :: axis%zaxisName)
      axis%zaxisName = zaxisName
    END IF
    IF (PRESENT(zaxisLongname)) THEN
      IF (ALLOCATED(axis%zaxisLongname)) CALL finish(routine, "'zaxisLongname' already initialized!")
      ALLOCATE(CHARACTER(LEN_TRIM(zaxisLongname)) :: axis%zaxisLongname)
      axis%zaxisLongname = zaxisLongname
    END IF
    IF (PRESENT(zaxisUnits)) THEN
      IF (ALLOCATED(axis%zaxisUnits))    CALL finish(routine, "'zaxisUnits' already initialized!")
      ALLOCATE(CHARACTER(LEN_TRIM(zaxisUnits)) :: axis%zaxisUnits)
      axis%zaxisUnits = zaxisUnits
    END IF
    IF (PRESENT(zaxisNumber)) THEN
      IF (ALLOCATED(axis%zaxisNumber))   CALL finish(routine, "'zaxisNumber' already initialized!")
      ALLOCATE(axis%zaxisNumber)
      axis%zaxisNumber = zaxisNumber
    END IF
    IF (PRESENT(zaxisNlevRef)) THEN
      IF (ALLOCATED(axis%zaxisNlevRef))   CALL finish(routine, "'zaxisNlevRef' already initialized!")
      ALLOCATE(axis%zaxisNlevRef)
      axis%zaxisNlevRef = zaxisNlevRef
    END IF
    IF (PRESENT(zaxisUUID)) THEN
      IF (ALLOCATED(axis%zaxisUUID))     CALL finish(routine, "'zaxisUUID' already initialized!")
      ALLOCATE(axis%zaxisUUID(SIZE(zaxisUUID)))
      axis%zaxisUUID = zaxisUUID
    END IF
    IF (PRESENT(zaxisVct)) THEN
      IF (ALLOCATED(axis%zaxisVct))      CALL finish(routine, "'zaxisVct' already initialized!")
      ALLOCATE(axis%zaxisVct(SIZE(zaxisVct)))
      axis%zaxisVct = zaxisVct
    END IF
    IF (PRESENT(zaxisDefLtype)) THEN
      IF (ALLOCATED(axis%zaxisDefLtype))  CALL finish(routine, "'zaxisDefLtype' already initialized!")
      ALLOCATE(axis%zaxisDefLtype)
      axis%zaxisDefLtype = zaxisDefLtype
    END IF
  END SUBROUTINE t_verticalAxis_set


  ! --------------------------------------------------------------------------------------
  !> Send object via MPI point-to-point communication.
  !
  SUBROUTINE t_verticalAxis_packMPI(axis, packedMessage)
    CLASS(t_verticalAxis),       INTENT(IN), TARGET    :: axis
    TYPE(t_PackedMessage),       INTENT(INOUT) :: packedMessage
    LOGICAL :: is_zaxisName, is_zaxisLongname, is_zaxisUnits, &
      & is_zaxisNumber, is_zaxisNlevRef

    is_zaxisName     = ALLOCATED(axis%zaxisName)
    is_zaxisLongname = ALLOCATED(axis%zaxisLongname)
    is_zaxisUnits    = ALLOCATED(axis%zaxisUnits)
    is_zaxisNumber   = ALLOCATED(axis%zaxisNumber)
    is_zaxisNlevRef  = ALLOCATED(axis%zaxisNlevRef)

    CALL packedMessage%pack(axis%zaxisNlev)
    CALL packedMessage%pack(axis%zaxisType%icon_zaxis_type)
    CALL packedMessage%pack(axis%zaxisType%cdi_zaxis_type)
    CALL packedMessage%pack(axis%zaxisType%is_2d)
    CALL packedMessage%pack(axis%cdi_id)
    CALL packedMessage%pack(axis%zaxisLevels)
    CALL packedMessage%pack(axis%zaxisLbounds)
    CALL packedMessage%pack(axis%zaxisUbounds)
    CALL packedMessage%pack(is_zaxisName)
    IF (is_zaxisName) CALL pack_string(axis%zaxisName)
    CALL packedMessage%pack(is_zaxisLongName)
    IF (is_zaxisLongName) CALL pack_string(axis%zaxisLongname)
    CALL packedMessage%pack(is_zaxisUnits)
    IF (is_zaxisUnits) CALL pack_string(axis%zaxisUnits)
    CALL packedMessage%pack(is_zaxisNumber)
    IF (is_zaxisNumber) CALL packedMessage%pack(axis%zaxisNumber)
    CALL packedMessage%pack(is_zaxisNlevRef)
    IF (is_zaxisNlevRef) CALL packedMessage%pack(axis%zaxisNlevRef)
    CALL packedMessage%pack(axis%zaxisUUID)
    CALL packedMessage%pack(axis%zaxisVct)
    CALL packedMessage%pack(axis%zaxisDefLType)
  CONTAINS

    SUBROUTINE pack_string(string)
      CHARACTER(*), INTENT(IN) :: string

      CALL packedMessage%pack(LEN(string))
      CALL packedMessage%pack(string)
    END SUBROUTINE pack_string
  END SUBROUTINE t_verticalAxis_packMPI


  ! --------------------------------------------------------------------------------------
  !> Receive object via MPI point-to-point communication.
  !
  SUBROUTINE t_verticalAxis_unpackMPI(axis, packedMessage)
    CLASS(t_verticalAxis),       INTENT(INOUT), TARGET :: axis
    TYPE(t_PackedMessage),       INTENT(INOUT) :: packedMessage
    LOGICAL :: is_zaxisName, is_zaxisLongname, is_zaxisUnits, &
      & is_zaxisNumber, is_zaxisNlevRef
    INTEGER :: strlen

    CALL packedMessage%unpack(axis%zaxisNlev)
    CALL packedMessage%unpack(axis%zaxisType%icon_zaxis_type)
    CALL packedMessage%unpack(axis%zaxisType%cdi_zaxis_type)
    CALL packedMessage%unpack(axis%zaxisType%is_2D)
    CALL packedMessage%unpack(axis%cdi_id)
    CALL packedMessage%unpack(axis%zaxisLevels)
    CALL packedMessage%unpack(axis%zaxisLbounds)
    CALL packedMessage%unpack(axis%zaxisUbounds)
    CALL packedMessage%unpack(is_zaxisName)
    IF (is_zaxisName) THEN
      CALL packedMessage%unpack(strlen)
      ALLOCATE(CHARACTER(strlen) :: axis%zaxisName)
      CALL unpack_string(axis%zaxisName)
    END IF
    CALL packedMessage%unpack(is_zaxisLongName)
    IF (is_zaxisLongName) THEN
      CALL packedMessage%unpack(strlen)
      ALLOCATE(CHARACTER(strlen) :: axis%zaxisLongname)
      CALL unpack_string(axis%zaxisLongname)
    END IF
    CALL packedMessage%unpack(is_zaxisUnits)
    IF (is_zaxisUnits) THEN
      CALL packedMessage%unpack(strlen)
      ALLOCATE(CHARACTER(strlen) :: axis%zaxisUnits)
      CALL unpack_string(axis%zaxisUnits)
    END IF
    CALL packedMessage%unpack(is_zaxisNumber)
    IF (is_zaxisNumber) THEN
      ALLOCATE(axis%zaxisNumber)
      CALL packedMessage%unpack(axis%zaxisNumber)
    END IF
    CALL packedMessage%unpack(is_zaxisNlevRef)
    IF (is_zaxisNlevRef) THEN
      ALLOCATE(axis%zaxisNlevRef)
      CALL packedMessage%unpack(axis%zaxisNlevRef)
    END IF
    CALL packedMessage%unpack(axis%zaxisUUID)
    CALL packedMessage%unpack(axis%zaxisVct)
    CALL packedMessage%unpack(axis%zaxisDefLtype)
  CONTAINS

    SUBROUTINE unpack_string(string)
      CHARACTER(*), INTENT(OUT) :: string

      CALL packedMessage%unpack(string)
    END SUBROUTINE unpack_string
  END SUBROUTINE t_verticalAxis_unpackMPI


  ! --------------------------------------------------------------------------------------
  !> Send object via MPI point-to-point communication.
  !
  SUBROUTINE t_verticalAxis_sendMPI(axis, dest, comm)
    CLASS(t_verticalAxis),       INTENT(IN) :: axis
    INTEGER,                     INTENT(IN) :: dest  !< destination PE
    INTEGER,                     INTENT(IN) :: comm  !< MPI communicator
    TYPE(t_PackedMessage) :: packedMessage
    CALL axis%packMPI(packedMessage)
    CALL packedMessage%send(dest, 0, comm)
  END SUBROUTINE t_verticalAxis_sendMPI


  ! --------------------------------------------------------------------------------------
  !> Receive object via MPI point-to-point communication.
  !
  SUBROUTINE t_verticalAxis_recvMPI(axis, src, comm)
    CLASS(t_verticalAxis),       INTENT(INOUT) :: axis
    INTEGER,                     INTENT(IN) :: src   !< source PE
    INTEGER,                     INTENT(IN) :: comm  !< MPI communicator
    TYPE(t_PackedMessage) :: packedMessage

    CALL packedMessage%recv(src, 0, comm)
    CALL axis%finalize()  ! clear receive object
    CALL axis%unpackMPI(packedMessage)
  END SUBROUTINE t_verticalAxis_recvMPI


  ! --------------------------------------------------------------------------------------
  !> Broadcast object via MPI communication.
  !
  SUBROUTINE t_verticalAxis_bcastMPI(axis, root, comm)
    CLASS(t_verticalAxis),   INTENT(INOUT) :: axis
    INTEGER,                 INTENT(IN)    :: root  !< root PE
    INTEGER,                 INTENT(IN)    :: comm  !< MPI communicator
    TYPE(t_PackedMessage) :: packedMessage
    LOGICAL               :: isRoot, isReceiver

    CALL p_get_bcast_role(root, comm, isRoot, isReceiver)
    IF (isRoot)  CALL axis%packMPI(packedMessage)
    CALL packedMessage%bcast(root, comm)
    IF (isReceiver)  CALL axis%unpackMPI(packedMessage)
  END SUBROUTINE t_verticalAxis_bcastMPI


  ! --------------------------------------------------------------------------------------
  !> Auxiliary function: Check two REAL(dp) arrays for equality.
  !
  LOGICAL FUNCTION array_eqv_dp(array1, array2)
    REAL(wp), ALLOCATABLE, INTENT(IN) :: array1(:), array2(:)
    
    array_eqv_dp = .TRUE.
    ! check pointer association status
    IF (.NOT. (ALLOCATED(array1) .EQV. ALLOCATED(array2))) THEN
      array_eqv_dp=.FALSE.
      RETURN
    END IF
    ! check array size
    IF (ALLOCATED(array1)) THEN
      IF (.NOT. (SIZE(array1) == SIZE(array2))) THEN
        array_eqv_dp=.FALSE.
      ELSE
        IF (.NOT. ALL(array1(:) == array2(:)))  array_eqv_dp=.FALSE.
      END IF
    END IF
  END FUNCTION array_eqv_dp


  ! --------------------------------------------------------------------------------------
  !> Auxiliary function: Check two REAL(dp) arrays for equality.
  !
  LOGICAL FUNCTION array_eqv_intc(array1, array2)
    INTEGER(KIND=C_SIGNED_CHAR), ALLOCATABLE, INTENT(IN) :: array1(:), array2(:)
    
    array_eqv_intc = .TRUE.
    ! check pointer association status
    IF (.NOT. (ALLOCATED(array1) .EQV. ALLOCATED(array2))) THEN
      array_eqv_intc=.FALSE.
      RETURN
    END IF
    ! check array size
    IF (ALLOCATED(array1)) THEN
      IF (.NOT. (SIZE(array1) == SIZE(array2))) THEN
        array_eqv_intc=.FALSE.
      ELSE
        IF (.NOT. ALL(array1(:) == array2(:)))  array_eqv_intc=.FALSE.
      END IF
    END IF
  END FUNCTION array_eqv_intc


  ! --------------------------------------------------------------------------------------
  !> Check two t_verticalAxis objects for equality.
  !
  LOGICAL FUNCTION t_verticalAxis_eqv(axis1, axis2)
    CLASS(t_verticalAxis), INTENT(IN) :: axis1, axis2

    t_verticalAxis_eqv =  .TRUE. 
    IF (.NOT. (axis1%zaxisType == axis2%zaxisType))  t_verticalAxis_eqv=.FALSE.
    IF (.NOT. (axis1%zaxisNlev == axis2%zaxisNlev))  t_verticalAxis_eqv=.FALSE.
    ! (note that we do not check the cdi_id!)

    ! check pointer association status 
    IF (.NOT. (ALLOCATED(axis1%zaxisName    ) .EQV. ALLOCATED(axis2%zaxisName    )))  t_verticalAxis_eqv=.FALSE.
    IF (.NOT. (ALLOCATED(axis1%zaxisLongname) .EQV. ALLOCATED(axis2%zaxisLongname)))  t_verticalAxis_eqv=.FALSE.
    IF (.NOT. (ALLOCATED(axis1%zaxisUnits   ) .EQV. ALLOCATED(axis2%zaxisUnits   )))  t_verticalAxis_eqv=.FALSE.
    IF (.NOT. (ALLOCATED(axis1%zaxisNumber  ) .EQV. ALLOCATED(axis2%zaxisNumber  )))  t_verticalAxis_eqv=.FALSE.
    IF (.NOT. (ALLOCATED(axis1%zaxisNlevRef ) .EQV. ALLOCATED(axis2%zaxisNlevRef )))  t_verticalAxis_eqv=.FALSE.

    ! check data members (if available):
    IF (t_verticalAxis_eqv) THEN
      IF (ALLOCATED(axis1%zaxisName)) THEN
        IF (axis1%zaxisName /= axis2%zaxisName)          t_verticalAxis_eqv=.FALSE.
      END IF
      IF (ALLOCATED(axis1%zaxisLongname)) THEN
        IF (axis1%zaxisLongname /= axis2%zaxisLongname)  t_verticalAxis_eqv=.FALSE.
      END IF
      IF (ALLOCATED(axis1%zaxisUnits)) THEN
        IF (axis1%zaxisUnits /= axis2%zaxisUnits)        t_verticalAxis_eqv=.FALSE.
      END IF
      IF (ALLOCATED(axis1%zaxisNumber)) THEN
        IF (axis1%zaxisNumber /= axis2%zaxisNumber)      t_verticalAxis_eqv=.FALSE.
      END IF
      IF (ALLOCATED(axis1%zaxisNlevRef)) THEN
        IF (axis1%zaxisNlevRef /= axis2%zaxisNlevRef)    t_verticalAxis_eqv=.FALSE.
      END IF
      IF (ALLOCATED(axis1%zaxisDefLtype)) THEN
        IF (axis1%zaxisDefLtype /= axis2%zaxisDefLtype)  t_verticalAxis_eqv=.FALSE.
      END IF
   
      ! check array contents (if available):
      IF (.NOT. array_eqv_dp(axis1%zaxisLevels,  axis2%zaxisLevels))   t_verticalAxis_eqv=.FALSE.
      IF (.NOT. array_eqv_dp(axis1%zaxisLbounds, axis2%zaxisLbounds))  t_verticalAxis_eqv=.FALSE.
      IF (.NOT. array_eqv_dp(axis1%zaxisUbounds, axis2%zaxisUbounds))  t_verticalAxis_eqv=.FALSE.
      IF (.NOT. array_eqv_dp(axis1%zaxisVct,     axis2%zaxisVct))      t_verticalAxis_eqv=.FALSE.
      IF (.NOT. array_eqv_intc(axis1%zaxisUUID,  axis2%zaxisUUID))     t_verticalAxis_eqv=.FALSE.
    END IF
  END FUNCTION t_verticalAxis_eqv


  ! --------------------------------------------------------------------------------------
  !> Constructor for t_verticalAxis type
  !
  FUNCTION new_verticalAxis(zaxisType, zaxisNlev,                                 &
    &                       zaxisLevels,  zaxisLbounds, zaxisUbounds, zaxisName,  &
    &                       zaxisLongname, zaxisUnits, zaxisNumber, zaxisNlevRef, &
    &                       zaxisUUID, zaxisVct, zaxisDefLtype)

    TYPE(t_verticalAxis) :: new_verticalAxis
    TYPE(t_zaxisType),           INTENT(IN)           :: zaxisType
    INTEGER,                     INTENT(IN), OPTIONAL :: zaxisNlev
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisLevels(:)
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisLbounds(:)
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisUbounds(:)
    CHARACTER(LEN=*),            INTENT(IN), OPTIONAL :: zaxisName
    CHARACTER(LEN=*),            INTENT(IN), OPTIONAL :: zaxisLongname
    CHARACTER(LEN=*),            INTENT(IN), OPTIONAL :: zaxisUnits
    INTEGER,                     INTENT(IN), OPTIONAL :: zaxisNumber
    INTEGER,                     INTENT(IN), OPTIONAL :: zaxisNlevRef
    INTEGER(KIND=C_SIGNED_CHAR), INTENT(IN), OPTIONAL :: zaxisUUID(:)
    REAL(dp),                    INTENT(IN), OPTIONAL :: zaxisVct(:)
    INTEGER,                     INTENT(IN), OPTIONAL :: zaxisDefLtype
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::new_verticalAxis'

    new_verticalAxis%zaxisType = zaxisType
    new_verticalAxis%zaxisNlev = 1
    IF (PRESENT(zaxisNlev)) new_verticalAxis%zaxisNlev = zaxisNlev
    IF (new_verticalAxis%zaxisNlev <= 0) THEN
      CALL finish(routine, "Invalid no. of vertical axis levels!")
    END IF

    CALL new_verticalAxis%set(zaxisLevels,  zaxisLbounds, zaxisUbounds, zaxisName,  &
      &                       zaxisLongname, zaxisUnits, zaxisNumber, zaxisNlevRef, &
      &                       zaxisUUID, zaxisVct, zaxisDefLtype)
  END FUNCTION new_verticalAxis


  ! --------------------------------------------------------------------------------------
  !> Destructor.
  !
  RECURSIVE SUBROUTINE t_verticalAxisList_finalize(axisList)
    CLASS(t_verticalAxisList), INTENT(INOUT) :: axisList
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::t_verticalAxisList_finalize'

    !CALL message(routine, "Deallocate vertical axes meta-data.")
    IF (ASSOCIATED(axisList%axis)) THEN
      CALL axisList%axis%finalize()
      DEALLOCATE(axisList%axis)
    END IF
    IF (ASSOCIATED(axisList%next)) THEN
      CALL axisList%next%finalize()
      DEALLOCATE(axisList%next)
    END IF
  END SUBROUTINE t_verticalAxisList_finalize


  ! --------------------------------------------------------------------------------------
  !> Append vertical axis to linked list.
  !  Duplicates are avoided.
  !
  SUBROUTINE t_verticalAxisList_append(axisList, axis)
    CLASS(t_verticalAxisList), TARGET, INTENT(INOUT) :: axisList
    TYPE(t_verticalAxis),              INTENT(IN)    :: axis
    TYPE(t_verticalAxisList), POINTER :: it

    ! if no identical zaxis object has been found, create a new one
    IF (.NOT. ASSOCIATED(axisList%axis)) THEN
      ALLOCATE(axisList%axis)
      axisList%axis = axis
      axisList%next => NULL()
    ELSE
      it => axisList
      DO
        ! do not add axis if it is already available
        IF (it%axis == axis)  EXIT

        IF (.NOT. ASSOCIATED(it%next)) THEN
          ALLOCATE(it%next)
          ALLOCATE(it%next%axis)
          it%next%axis = axis
          it%next%next => NULL()
          EXIT
        END IF
        it => it%next
      END DO
    END IF
  END SUBROUTINE t_verticalAxisList_append


  ! --------------------------------------------------------------------------------------
  !> Append vertical axis list to this list.
  !  Duplicates are avoided.
  !
  SUBROUTINE t_verticalAxisList_append_list(axisList1, axisList2)
    CLASS(t_verticalAxisList), TARGET, INTENT(INOUT) :: axisList1
    TYPE(t_verticalAxisList),  TARGET, INTENT(IN)    :: axisList2
    TYPE(t_verticalAxisList), POINTER :: it

    it => axisList2
    DO
      IF (.NOT. ASSOCIATED(it))  EXIT
      CALL axisList1%append(it%axis)
      it => it%next
    END DO
  END SUBROUTINE t_verticalAxisList_append_list


  ! --------------------------------------------------------------------------------------
  !> Search linked list for a specific axis.
  !
  !  The search result must match all optional arguments present
  !  (name, type, ...).
  !
  FUNCTION t_verticalAxisList_getEntry(axisList, &
    &                                  icon_zaxis_type, cdi_zaxis_type, &
    &                                  zaxisName, zaxisNlev)  RESULT(axis)

    TYPE(t_verticalAxis), POINTER :: axis
    CLASS(t_verticalAxisList), TARGET, INTENT(INOUT) :: axisList
    INTEGER,          OPTIONAL,        INTENT(IN)    :: icon_zaxis_type
    INTEGER,          OPTIONAL,        INTENT(IN)    :: cdi_zaxis_type
    CHARACTER(LEN=*), OPTIONAL,        INTENT(IN)    :: zaxisName
    INTEGER,          OPTIONAL,        INTENT(IN)    :: zaxisNlev
    TYPE(t_verticalAxisList), POINTER :: it
    LOGICAL :: match

    it   => axisList
    axis => NULL()
    DO
      IF (.NOT. ASSOCIATED(it))      EXIT
      IF (.NOT. ASSOCIATED(it%axis)) EXIT

      match = .TRUE.
      IF (PRESENT(icon_zaxis_type)) THEN
        IF (it%axis%zaxisType%icon_zaxis_type /= icon_zaxis_type)  match=.FALSE.
      END IF
      IF (PRESENT(cdi_zaxis_type)) THEN
        IF (it%axis%zaxisType%cdi_zaxis_type /= cdi_zaxis_type)  match=.FALSE.
      END IF
      IF (PRESENT(zaxisName)) THEN
        IF (it%axis%zaxisName /= zaxisName)  match=.FALSE.
      END IF
      IF (PRESENT(zaxisNlev)) THEN
        IF (it%axis%zaxisNlev /= zaxisNlev)  match=.FALSE.
      END IF
      IF (match) THEN
        axis => it%axis ; EXIT
      END IF
      it => it%next
    END DO
  END FUNCTION t_verticalAxisList_getEntry


  ! --------------------------------------------------------------------------------------
  !> Send object via MPI point-to-point communication.
  !
  SUBROUTINE t_verticalAxisList_sendMPI(axisList, dest, comm)
    CLASS(t_verticalAxisList), TARGET, INTENT(IN) :: axisList
    INTEGER,                           INTENT(IN) :: dest  !< destination PE
    INTEGER,                           INTENT(IN) :: comm  !< MPI communicator
    
    TYPE(t_verticalAxisList), POINTER :: it
    TYPE(t_verticalAxis) :: sentinel

    it => axisList
    DO
      IF (.NOT. ASSOCIATED(it%next))  EXIT
      CALL it%next%axis%sendMPI(dest,comm)
      it => it%next
    END DO
    ! the last item of the list is marked by a subsequent "sentinel
    ! message":
    sentinel%zaxisNlev = MSG_SENTINEL
    CALL sentinel%sendMPI(dest,comm)
  END SUBROUTINE t_verticalAxisList_sendMPI


  ! --------------------------------------------------------------------------------------
  !> Receive object via MPI point-to-point communication.
  !
  SUBROUTINE t_verticalAxisList_recvMPI(axisList, src, comm)
    CLASS(t_verticalAxisList),   INTENT(INOUT) :: axisList
    INTEGER,                     INTENT(IN) :: src   !< source PE
    INTEGER,                     INTENT(IN) :: comm  !< MPI communicator

    TYPE(t_verticalAxis),     POINTER :: axis
    DO
      ALLOCATE(axis)
      CALL axis%recvMPI(src,comm)
      IF (axis%zaxisNlev == MSG_SENTINEL) THEN
        DEALLOCATE(axis)
        EXIT
      END IF
      CALL axisList%append(axis)
    END DO
  END SUBROUTINE t_verticalAxisList_recvMPI


  ! --------------------------------------------------------------------------------------
  !> Receive object via MPI point-to-point communication.
  !
  SUBROUTINE t_verticalAxisList_bcastMPI(axisList, root, comm)
    CLASS(t_verticalAxisList), INTENT(INOUT), TARGET :: axisList
    INTEGER,                   INTENT(IN) :: root  !< root PE
    INTEGER,                   INTENT(IN) :: comm  !< MPI communicator

    TYPE(t_verticalAxisList), POINTER :: it
    TYPE(t_verticalAxis),     POINTER :: axis
    TYPE(t_verticalAxis)              :: sentinel
    LOGICAL                           :: isRoot, isReceiver

    CALL p_get_bcast_role(root, comm, isRoot, isReceiver)
    IF (isRoot) THEN
      it => axisList
      DO
        IF (.NOT. ASSOCIATED(it%next))  EXIT
        CALL it%next%axis%bcastMPI(root,comm)
        it => it%next
      END DO
      ! the last item of the list is marked by a subsequent "sentinel
      ! message":
      sentinel%zaxisNlev = MSG_SENTINEL
      CALL sentinel%bcastMPI(root,comm)
    END IF
    IF (isReceiver) THEN
      DO
        ALLOCATE(axis)
        CALL axis%bcastMPI(root,comm)
        IF (axis%zaxisNlev == MSG_SENTINEL) THEN
          DEALLOCATE(axis)
          EXIT
        END IF
        CALL axisList%append(axis)
      END DO
    END IF
  END SUBROUTINE t_verticalAxisList_bcastMPI


  ! --------------------------------------------------------------------------------------
  !> Constructor for t_verticalAxisList type
  FUNCTION new_verticalAxisList()
    TYPE(t_verticalAxisList) :: new_verticalAxisList

    new_verticalAxisList%axis => NULL()
    new_verticalAxisList%next => NULL()
  END FUNCTION new_verticalAxisList

END MODULE mo_name_list_output_zaxes_types
