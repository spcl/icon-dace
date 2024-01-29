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

!> Wrapper module containing the Fortran-C-Interface for the
!  Fortran namelist scanner.

MODULE mo_util_nml

  USE, INTRINSIC ::  ISO_C_BINDING, ONLY: c_int, c_char, c_null_char

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    FUNCTION private_annotate_nml(in_filename, out_filename) RESULT(iret) &
      &      BIND(C, NAME='util_annotate_nml')
      IMPORT :: c_int, c_char
      INTEGER(c_int) :: iret
      CHARACTER(c_char), DIMENSION(*), INTENT(in) :: in_filename, out_filename
    END FUNCTION private_annotate_nml
  END INTERFACE

  PUBLIC :: util_annotate_nml

CONTAINS

  FUNCTION util_annotate_nml(in_filename, out_filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: in_filename, out_filename
    iret = private_annotate_nml(TRIM(in_filename)//c_null_char, &
      &                         TRIM(out_filename)//c_null_char)
  END FUNCTION util_annotate_nml

END MODULE mo_util_nml
