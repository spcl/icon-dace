!> @file comin_plugin_types.F90
!! @brief Data type Definition and variables describing the dynamic libs (plugins).
!
!  @authors 06/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_plugin_types

  USE ISO_C_BINDING, ONLY : C_INT, C_CHAR
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_comin_plugin_description
  PUBLIC :: t_comin_plugin_info

#include "comin_global.inc"

  !> Data type, describing the dynamic libraries
  !! @ingroup host_interface
  TYPE :: t_comin_plugin_description
    ! name of the plugin - currently only used for messages
    CHARACTER(LEN=MAX_LEN_PLUGIN_NAME) :: name = ""

    ! full name of plugin shared library (including `.so` file
    ! extension) or "icon" for static linking.
    CHARACTER(LEN=MAX_LEN_PLUGIN_LIBRARY) :: plugin_library = ""

    ! name of primary constructor.
    CHARACTER(LEN=MAX_LEN_PRIMARY_CONSTRUCTOR)  :: primary_constructor  = "comin_main"

    ! options string: offers the possibility to pass a character
    ! string (e.g. a python script filename) to the plugin.
    CHARACTER(LEN=MAX_LEN_OPTIONS)  :: options  = ""

    ! name of MPI communicator. left as an empty string if the
    ! application does not require a communicator for this plugin.
    CHARACTER(LEN=MAX_LEN_COMM)  :: comm        = ""

  END TYPE t_comin_plugin_description

  !> The elements of this derived data type describe a 3rd party plugin.
  !! @ingroup plugin_interface
  TYPE :: t_comin_plugin_info
    INTEGER                       :: id
    CHARACTER(LEN=:), ALLOCATABLE :: name
    CHARACTER(LEN=:), ALLOCATABLE :: options
    CHARACTER(LEN=:), ALLOCATABLE :: comm
    LOGICAL                       :: errors_return = .FALSE.
  END TYPE t_comin_plugin_info

END MODULE comin_plugin_types
