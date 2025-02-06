!> @file comin_errhandler_types.F90
!! @brief Common data structures for error handling.
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_errhandler_types

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_host_errhandler_fct

  INTERFACE
    !> In order to be compatible with ICON, the interface contains
    !> OPTIONAL arguments which are probably not C-compliant.
    SUBROUTINE comin_host_errhandler_fct(routine, text)
      CHARACTER(LEN=*), INTENT(IN) :: routine
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE comin_host_errhandler_fct
  END INTERFACE

END MODULE comin_errhandler_types
