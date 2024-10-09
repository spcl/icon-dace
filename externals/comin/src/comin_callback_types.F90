!> @file comin_callback_types.F90
!! @brief Type definitions for handling third party plugin callbacks.
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_callback_types
  USE comin_plugin_types, ONLY: t_comin_plugin_info
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_callback_routine
  PUBLIC :: t_comin_callback_element
  PUBLIC :: t_comin_callback_context
  PUBLIC :: t_callback_list
  PUBLIC :: t_callback_list_item

  ABSTRACT INTERFACE
    SUBROUTINE comin_callback_routine() BIND(C)
    END SUBROUTINE comin_callback_routine
  END INTERFACE

  !> information about each entry point/callback
  TYPE :: t_comin_callback_element
    INTEGER :: entry_point_id
    TYPE(t_comin_plugin_info), POINTER :: plugin_info
    PROCEDURE(comin_callback_routine), POINTER, NOPASS:: comin_callback => NULL()
  END TYPE t_comin_callback_element

  !> linked list of all callbacks, items have TYPE(t_comin_callback_element)
#include "comin_callback_linked_list_header.inc"

  !> Array of variable lists (array of pointer lists)
  !  each entry stores the lists of callbacks registered for the
  !  context (dimension of array)
  TYPE :: t_comin_callback_context
    TYPE(t_comin_callback_element), POINTER :: vl => NULL()
  END TYPE t_comin_callback_context

CONTAINS

#include "comin_callback_linked_list_body.inc"

END MODULE comin_callback_types
