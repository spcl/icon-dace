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

MODULE mo_util_uuid_types
  
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIGNED_CHAR
  
  IMPLICIT NONE 
  
  PRIVATE
  
  PUBLIC :: t_uuid
  PUBLIC :: UUID_STRING_LENGTH
  PUBLIC :: UUID_DATA_LENGTH
  
  INTEGER, PARAMETER :: UUID_STRING_LENGTH = 36
  INTEGER, PARAMETER :: UUID_DATA_LENGTH   = 16
  
  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_uuid_types'
  
  TYPE, BIND(C) :: t_uuid
    INTEGER(C_SIGNED_CHAR) :: DATA(16)
  END TYPE t_uuid

END MODULE mo_util_uuid_types
