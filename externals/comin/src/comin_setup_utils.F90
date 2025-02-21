!> @file comin_setup_utils.F90
!! @brief ComIn setup utilities, containing version compatibility checks, ComIn setup, Third party plugin setup
!
!  @authors 01/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_setup_utils

  USE iso_c_binding,      ONLY: c_int
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_comin_setup_version_info, comin_setup_get_version
  PUBLIC :: comin_setup_version_compatible






  !> The elements of this derived data type describe the current
  !> community interface.
  !! @ingroup constants
  TYPE, BIND(C) :: t_comin_setup_version_info
    !> ComIn versioning (major and minor version info)
     INTEGER(kind=c_int) :: version_no_major, version_no_minor, &
          & version_no_patch
  END TYPE t_comin_setup_version_info


CONTAINS

  !  Convention: The major version has to agree between ComIn and
  !              ICON/3rd party plugins and minor version should be
  !              backward compatible.
  FUNCTION comin_setup_version_compatible(setupA, setupB)  RESULT(lret)
    LOGICAL :: lret
    TYPE(t_comin_setup_version_info), INTENT(IN) :: setupA, setupB

    lret = ( setupA%version_no_major == setupB%version_no_major )
  END FUNCTION comin_setup_version_compatible

  !> Returns version info.
  !! @ingroup constants
  !!
  !! `BIND(C)` is needed to have it as symbol
  !! `comin_setup_get_version` in the shared library
  FUNCTION comin_setup_get_version() BIND(C)
    TYPE(t_comin_setup_version_info) :: comin_setup_get_version
    comin_setup_get_version%version_no_major = 0
    comin_setup_get_version%version_no_minor = 1
    comin_setup_get_version%version_no_patch = 0
  END FUNCTION comin_setup_get_version


END MODULE comin_setup_utils
