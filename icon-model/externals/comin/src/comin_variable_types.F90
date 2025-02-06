!> @file comin_variable_types.F90
!! @brief Data types for variable definition
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_variable_types

  USE iso_c_binding,           ONLY: c_int, c_char, c_bool, c_ptr, c_null_ptr
  USE comin_setup_constants,   ONLY: wp, COMIN_ZAXIS_3D, &
    &                                COMIN_HGRID_UNSTRUCTURED_CELL
  USE comin_metadata_types,    ONLY: t_comin_var_metadata
  IMPLICIT NONE

  PUBLIC

#include "comin_global.inc"

  ! ------------------------------------
  ! data types for variable definition
  ! ------------------------------------

  !> Variable descriptor.
  !> identifies (uniquely) a variable. Do not confuse with meta-data
  !! @ingroup common
  TYPE :: t_comin_var_descriptor
    CHARACTER(LEN=:), ALLOCATABLE :: name
    ! domain id
    INTEGER                       :: id
  END TYPE t_comin_var_descriptor

  TYPE, BIND(C) :: t_comin_var_descriptor_c
    CHARACTER(KIND=c_char) :: name(MAX_LEN_VAR_NAME+1)
    ! domain id
    INTEGER(kind=c_int)    :: id
  END TYPE t_comin_var_descriptor_c

  !> Variable pointer. Basically this wraps a 5-dimensional REAL(wp)
  !> pointer together with some interpretation of the different array
  !> dimensions.
  !! @ingroup common
  TYPE :: t_comin_var_ptr
    !> the var_descriptor
    TYPE(t_comin_var_descriptor) :: descriptor

    !> the (current) pointer to the data
    REAL(wp), POINTER :: ptr(:,:,:,:,:) => NULL()

    !> the (current) device ptr to the data (if any)
    TYPE(c_ptr) :: device_ptr = c_null_ptr

    ! index positions in the 5D array.
    INTEGER :: pos_jc = -1, pos_jk = -1, pos_jb = -1, pos_jn = -1

    !> if (tracer==.TRUE.) and (ncontained > 0), then the variable
    !  pointer refers to an array slice pointer
    !  ptr(:,:,:,:,ncontained)
    INTEGER :: ncontained = 0
    !> LOGICAL flag. TRUE, if this is a container (contains variables)
    LOGICAL(kind=c_bool) :: lcontainer = .FALSE.
  END TYPE t_comin_var_ptr

  !> Variable item
  TYPE :: t_comin_var_item
    TYPE(t_comin_var_metadata)           :: metadata
    TYPE(t_comin_var_ptr),   POINTER     :: p => NULL()
  END TYPE t_comin_var_item

  !> Variable list for context access
  TYPE, EXTENDS(t_comin_var_item) :: t_comin_var_context_item
    INTEGER :: access_flag
  END TYPE t_comin_var_context_item

  !> Information on requested variables
  TYPE :: t_comin_request_item
    TYPE(t_comin_var_descriptor)    :: descriptor
    TYPE(t_comin_var_metadata)      :: metadata
    INTEGER, ALLOCATABLE            :: moduleID(:)

    !> LOGICAL flag. TRUE, if this variable is intended to be used
    !> exclusively by a particular 3rd party plugin:
    LOGICAL(kind=c_bool) :: lmodexclusive = .FALSE.
  END TYPE t_comin_request_item

  INTERFACE
    SUBROUTINE comin_var_sync_device_mem_fct(var_ptr, direction)
      import t_comin_var_ptr
      TYPE(t_comin_var_ptr), POINTER, INTENT(IN) :: var_ptr
      LOGICAL, INTENT(IN) :: direction
    END SUBROUTINE comin_var_sync_device_mem_fct
  END INTERFACE

  ! ------------------------------------
  ! lists of exposed variables
  ! ------------------------------------

#include "comin_var_linked_list_header.inc"

  !> Array of variable lists (array of pointer lists) each entry
  !  stores the lists of variables registered for the context
  !  (dimension of array) points to the first element of the variable
  !  list
  !  - contains TYPE(t_comin_var_item)
  TYPE :: t_comin_var_list_context
    TYPE(t_var_context_list) :: var_list
  END TYPE t_comin_var_list_context

CONTAINS

#include "comin_var_linked_list_body.inc"

  !> compare two variable descriptors.
  FUNCTION comin_var_descr_match(var_descriptor1, var_descriptor2)
    TYPE(t_comin_var_descriptor), INTENT(IN) :: var_descriptor1, var_descriptor2
    LOGICAL :: comin_var_descr_match
    ! local
    LOGICAL :: l_name, l_domain

    l_name = (TRIM(ADJUSTL(var_descriptor1%name)) == TRIM(ADJUSTL(var_descriptor2%name)))
    l_domain = (var_descriptor1%id == var_descriptor2%id)

    comin_var_descr_match = l_name .AND. l_domain
  END FUNCTION comin_var_descr_match

END MODULE comin_variable_types
