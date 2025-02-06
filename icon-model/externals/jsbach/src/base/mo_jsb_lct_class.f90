!> Contains definitions and methods for landcover types.
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
!>
!>### Landcover types
!>
!> This module provides possible ids for landcover types (lcts), a constructor for new lct instances and
!> further methods.
!>
!> Each lct has a named id. Furthermore, it can have a name, a longname, a fraction,
!> and in case that it is a PFT tile an id for the lctlib.
!>
!> lcts are used in usecases upon tile generation.
!> Each tile has a primary lct and, in addition, contains all lcts of its descendant tiles.
!> Leaf tiles only have one lct.
!> The lcts of a tile are used to define properties of the tile.
!> They are used to determine the logical flow in update routines;
!> which variables need to be allocated in process memories (for primary ltc but also for aggregation of
!> contained lcts); and they might provide an lctlib id as well as fractional shares of a tile.
!>
!> @todo currently the constructor does not catch cases where the id is not contained in the enumerator
!>
!> @todo clarify language by adding another term and renaming members of t_jsb_lct:
!>       -> issue: the id of the lct is not a unique id of the created t_jsb_lct, at least not in case of PFTs
!>                 (all PFTs share the same id (VEG_TYPE) and only differ in fract and lib_id?!)
!>       => suggestion: call it type_id
!>
!> @todo change default id to -1 instead of 0?
!>
MODULE mo_jsb_lct_class
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message_text, message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LAND_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE, LAKE_TYPE
  PUBLIC :: max_no_of_lct, t_jsb_lct, Contains_lct

  INTEGER, PARAMETER :: max_no_of_lct = 10

  ENUM, BIND(C)
    ENUMERATOR ::      &
      & LAND_TYPE = 1, & !< Land type to specify tile that computes the implicit surface energy balance
      & VEG_TYPE,      & !< Vegetated
      & BARE_TYPE,     & !< Bare soil -- currently not used
      & GLACIER_TYPE,  & !< Glacier
      & LAKE_TYPE        !< Lake
  END ENUM

  TYPE t_jsb_lct
    INTEGER           :: id = 0               !< Id for lct, one of enumeration above
    CHARACTER(LEN=5)  :: name = ""            !< Short name for lct
    CHARACTER(LEN=30) :: long_name = ""       !< Long name for lct
    REAL(wp), POINTER :: fract(:,:) => NULL() !< Cover fraction of this lct (relative to tile that owns this lct)
    INTEGER           :: lib_id               !< Id of lct in lctlib
  CONTAINS
    PROCEDURE :: Get_name     => Get_lct_type_name
    PROCEDURE :: Get_longname => Get_lct_type_longname
!!$    PROCEDURE :: Contains
  END TYPE t_jsb_lct

  INTERFACE t_jsb_lct
    PROCEDURE Construct_jsb_lct
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lct_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for t_jsb_lct
  !
  FUNCTION Construct_jsb_lct(id, name, long_name, lib_id) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in)                    :: id            !< the lct type id
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: name          !< the name of this lct
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: long_name     !< a long name for this lct
    INTEGER,          INTENT(in), OPTIONAL :: lib_id        !< in case of PFTs the id in the lctlib file
    TYPE(t_jsb_lct)                        :: return_value  !< the new instance
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Construct_jsb_lct'
    ! -------------------------------------------------------------------------------------------------- !

    return_value%id = id
    IF (PRESENT(name)) THEN
      return_value%name = TRIM(name)
    ELSE
      return_value%name = TRIM(return_value%Get_name())
    END IF

    IF (PRESENT(long_name)) THEN
      return_value%long_name = TRIM(long_name)
    ELSE
      return_value%long_name = TRIM(return_value%Get_longname())
    END IF

    IF (PRESENT(lib_id)) THEN
      return_value%lib_id = lib_id
    ELSE
      return_value%lib_id = 0
    END IF

  END FUNCTION Construct_jsb_lct

!!$  FUNCTION Contains(this, lct_type) RESULT(return_value)
!!$
!!$    CLASS(t_jsb_lct), INTENT(in) :: this
!!$    INTEGER,         INTENT(in) :: lct_type
!!$    LOGICAL                     :: return_value
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = modname//':Contains'
!!$
!!$    return_value = Contains_lct(this%id, lct_type)
!!$
!!$  END FUNCTION Contains
!!$

  ! ====================================================================================================== !
  !
  !> Checks if one of the given lcts has the given type
  !
  FUNCTION Contains_lct(lcts, lct_id) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lct), INTENT(in) :: lcts(:)
    INTEGER,         INTENT(in) :: lct_id
    LOGICAL                     :: return_value
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: ilct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Contains_lct'
    ! -------------------------------------------------------------------------------------------------- !

    return_value = .FALSE.
    DO ilct=1,SIZE(lcts)
      IF (lcts(ilct)%id == lct_id) THEN
        return_value = .TRUE.
        EXIT
      END IF
    END DO

  END FUNCTION Contains_lct

  ! ====================================================================================================== !
  !
  !> If the lct that is passed is named, the name is returned, else the default name for the type id
  !
  FUNCTION Get_lct_type_name(this) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_lct),  INTENT(in)  :: this
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_lct_type_name'
    ! -------------------------------------------------------------------------------------------------- !
    IF (this%name /= '') THEN
      return_value = TRIM(this%name)
    ELSE
      SELECT CASE(this%id)
      CASE (LAND_TYPE)
        return_value = 'land'
      CASE (VEG_TYPE)
        return_value = 'veg'
      CASE (BARE_TYPE)
        return_value = 'bare'
      CASE (GLACIER_TYPE)
        return_value = 'glac'
      CASE (LAKE_TYPE)
        return_value = 'lake'
      END SELECT
    END IF

  END FUNCTION Get_lct_type_name

  ! ====================================================================================================== !
  !
  !> Retuns the long name of the passed lct, or - if not specified - the default name for the type id
  !
  FUNCTION Get_lct_type_longname(this) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_lct),  INTENT(in)  :: this
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_lct_type_longname'
    ! -------------------------------------------------------------------------------------------------- !
    IF (this%long_name /= '') THEN
      return_value = TRIM(this%long_name)
    ELSE
      SELECT CASE(this%id)
      CASE (LAND_TYPE)
        return_value = 'Land (veg/soil/glac)'
      CASE (VEG_TYPE)
        return_value = 'Vegetated'
      CASE (BARE_TYPE)
        return_value = 'Bare soil'
      CASE (GLACIER_TYPE)
        return_value = 'Glacier'
      CASE (LAKE_TYPE)
        return_value = 'Lake'
      END SELECT
    END IF

  END FUNCTION Get_lct_type_longname

#endif
END MODULE mo_jsb_lct_class
