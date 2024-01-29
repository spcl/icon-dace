! This module provides fortran interface for OpenACC routines,
! because the OpenACC standard only defines C-Interfaces.
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

! (GZ, 2013-08-30): So far, the Cray compiler is the only one for which an OpenMP parallelization
! of copying data into / back from the MPI-buffer seems to give a benefit. Further compilers may
! be added here once the OpenMP implementation is sufficiently efficient

!----------------------------
#include "icon_definitions.inc"
!----------------------------
MODULE mo_openacc
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
#ifdef _OPENACC
USE mo_kind,                 ONLY: dp, sp
USE mo_mpi,                  ONLY: p_real_dp_byte, p_int_byte, p_real_sp_byte
USE iso_c_binding,           ONLY: c_ptr, c_size_t, c_loc

IMPLICIT NONE

PRIVATE

PUBLIC :: acc_deviceptr, acc_hostptr, &
          acc_map_data, acc_unmap_data, &
          acc_malloc, acc_free

!
!------------------------------------------------------------------------------------------------
!

CHARACTER(*), PARAMETER :: modname = "mo_openacc"

   INTERFACE
      SUBROUTINE c_acc_map_data(host, device, bytes) &
          bind(C, name='acc_map_data')
        IMPORT :: c_ptr, c_size_t
        TYPE(c_ptr), VALUE, INTENT(IN) :: host
        TYPE(c_ptr), VALUE, INTENT(IN) :: device
        INTEGER(c_size_t), VALUE, INTENT(IN) :: bytes
      END SUBROUTINE c_acc_map_data
   END INTERFACE
   INTERFACE
      SUBROUTINE c_acc_unmap_data(host) &
          bind(C, name='acc_unmap_data')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(IN) :: host
      END SUBROUTINE c_acc_unmap_data
   END INTERFACE
   INTERFACE
      FUNCTION c_acc_deviceptr(host) &
          bind(C, name='acc_deviceptr') RESULT(device)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(IN) :: host
        TYPE(c_ptr) :: device
      END FUNCTION c_acc_deviceptr
   END INTERFACE
   INTERFACE
      FUNCTION c_acc_hostptr(device) &
          bind(C, name='acc_hostptr') RESULT(host)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(IN) :: device
        TYPE(c_ptr) :: host
      END FUNCTION c_acc_hostptr
   END INTERFACE

   INTERFACE acc_deviceptr
      MODULE PROCEDURE acc_deviceptr_c_ptr
      MODULE PROCEDURE acc_deviceptr_1d_dp
      MODULE PROCEDURE acc_deviceptr_1d_sp
      MODULE PROCEDURE acc_deviceptr_1d_int
      MODULE PROCEDURE acc_deviceptr_2d_dp
      MODULE PROCEDURE acc_deviceptr_2d_sp
      MODULE PROCEDURE acc_deviceptr_2d_int
      MODULE PROCEDURE acc_deviceptr_3d_dp
      MODULE PROCEDURE acc_deviceptr_3d_sp
      MODULE PROCEDURE acc_deviceptr_3d_int
   END INTERFACE acc_deviceptr

   INTERFACE acc_hostptr
      MODULE PROCEDURE acc_hostptr_c_ptr
      MODULE PROCEDURE acc_hostptr_1d_dp
      MODULE PROCEDURE acc_hostptr_1d_sp
      MODULE PROCEDURE acc_hostptr_1d_int
      MODULE PROCEDURE acc_hostptr_2d_dp
      MODULE PROCEDURE acc_hostptr_2d_sp
      MODULE PROCEDURE acc_hostptr_2d_int
      MODULE PROCEDURE acc_hostptr_3d_dp
      MODULE PROCEDURE acc_hostptr_3d_sp
      MODULE PROCEDURE acc_hostptr_3d_int
   END INTERFACE acc_hostptr

   INTERFACE acc_map_data
      MODULE PROCEDURE acc_map_data_c_ptr
      MODULE PROCEDURE acc_map_data_1d_dp
      MODULE PROCEDURE acc_map_data_1d_sp
      MODULE PROCEDURE acc_map_data_1d_int
      MODULE PROCEDURE acc_map_data_2d_dp
      MODULE PROCEDURE acc_map_data_2d_sp
      MODULE PROCEDURE acc_map_data_2d_int
      MODULE PROCEDURE acc_map_data_3d_dp
      MODULE PROCEDURE acc_map_data_3d_sp
      MODULE PROCEDURE acc_map_data_3d_int
   END INTERFACE acc_map_data

   INTERFACE acc_unmap_data
      MODULE PROCEDURE acc_unmap_data_c_ptr
      MODULE PROCEDURE acc_unmap_data_1d_dp
      MODULE PROCEDURE acc_unmap_data_1d_sp
      MODULE PROCEDURE acc_unmap_data_1d_int
      MODULE PROCEDURE acc_unmap_data_2d_dp
      MODULE PROCEDURE acc_unmap_data_2d_sp
      MODULE PROCEDURE acc_unmap_data_2d_int
      MODULE PROCEDURE acc_unmap_data_3d_dp
      MODULE PROCEDURE acc_unmap_data_3d_sp
      MODULE PROCEDURE acc_unmap_data_3d_int
   END INTERFACE acc_unmap_data

   INTERFACE
      FUNCTION acc_malloc(bytes) bind(C, name='acc_malloc') RESULT(device)
        IMPORT :: c_ptr, c_size_t
        INTEGER(c_size_t), VALUE, INTENT(IN) :: bytes
        TYPE(c_ptr) :: device
      END FUNCTION acc_malloc
   END INTERFACE

   INTERFACE
      SUBROUTINE acc_free(device) bind(C, name='acc_free')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(IN) :: device
      END SUBROUTINE acc_free
   END INTERFACE

!-------------------------------------------------------------------------

CONTAINS

  FUNCTION acc_deviceptr_c_ptr(host) RESULT(device)
    TYPE(c_ptr), INTENT(IN), TARGET :: host
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(host)
  END FUNCTION acc_deviceptr_c_ptr

  FUNCTION acc_deviceptr_1d_dp(host) RESULT(device)
    REAL(dp), INTENT(IN), TARGET :: host(*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_1d_dp

  FUNCTION acc_deviceptr_1d_sp(host) RESULT(device)
    REAL(sp), INTENT(IN), TARGET :: host(*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_1d_sp

  FUNCTION acc_deviceptr_1d_int(host) RESULT(device)
    INTEGER, INTENT(IN), TARGET :: host(*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_1d_int

  FUNCTION acc_deviceptr_2d_dp(host) RESULT(device)
    REAL(dp), INTENT(IN), TARGET :: host(1,*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_2d_dp

  FUNCTION acc_deviceptr_2d_sp(host) RESULT(device)
    REAL(sp), INTENT(IN), TARGET :: host(1,*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_2d_sp

  FUNCTION acc_deviceptr_2d_int(host) RESULT(device)
    INTEGER, INTENT(IN), TARGET :: host(1,*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_2d_int

  FUNCTION acc_deviceptr_3d_dp(host) RESULT(device)
    REAL(dp), INTENT(IN), TARGET :: host(1,1,*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_3d_dp

  FUNCTION acc_deviceptr_3d_sp(host) RESULT(device)
    REAL(sp), INTENT(IN), TARGET :: host(1,1,*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_3d_sp

  FUNCTION acc_deviceptr_3d_int(host) RESULT(device)
    INTEGER, INTENT(IN), TARGET :: host(1,1,*)
    TYPE(c_ptr) :: device
    device = c_acc_deviceptr(c_loc(host))
  END FUNCTION acc_deviceptr_3d_int

  FUNCTION acc_hostptr_c_ptr(device) RESULT(host)
    TYPE(c_ptr), INTENT(IN), TARGET :: device
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(device)
  END FUNCTION acc_hostptr_c_ptr

  FUNCTION acc_hostptr_1d_dp(device) RESULT(host)
    REAL(dp), INTENT(IN), TARGET :: device(*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_1d_dp

  FUNCTION acc_hostptr_1d_sp(device) RESULT(host)
    REAL(sp), INTENT(IN), TARGET :: device(*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_1d_sp

  FUNCTION acc_hostptr_1d_int(device) RESULT(host)
    INTEGER, INTENT(IN), TARGET :: device(*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_1d_int

  FUNCTION acc_hostptr_2d_dp(device) RESULT(host)
    REAL(dp), INTENT(IN), TARGET :: device(1,*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_2d_dp

  FUNCTION acc_hostptr_2d_sp(device) RESULT(host)
    REAL(sp), INTENT(IN), TARGET :: device(1,*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_2d_sp

  FUNCTION acc_hostptr_2d_int(device) RESULT(host)
    INTEGER, INTENT(IN), TARGET :: device(1,*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_2d_int

  FUNCTION acc_hostptr_3d_dp(device) RESULT(host)
    REAL(dp), INTENT(IN), TARGET :: device(1,1,*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_3d_dp

  FUNCTION acc_hostptr_3d_sp(device) RESULT(host)
    REAL(sp), INTENT(IN), TARGET :: device(1,1,*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_3d_sp

  FUNCTION acc_hostptr_3d_int(device) RESULT(host)
    INTEGER, INTENT(IN), TARGET :: device(1,1,*)
    TYPE(c_ptr) :: host
    host = c_acc_hostptr(c_loc(device))
  END FUNCTION acc_hostptr_3d_int

  SUBROUTINE acc_map_data_c_ptr(host, device, data_size)
    TYPE(c_ptr), INTENT(IN) :: host
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data(host, device, INT(data_size, c_size_t))
  END SUBROUTINE acc_map_data_c_ptr

  SUBROUTINE acc_map_data_1d_dp(host, device, data_size)
    REAL(dp), INTENT(IN), TARGET :: host(*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_real_dp_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_1d_dp

  SUBROUTINE acc_map_data_1d_sp(host, device, data_size)
    REAL(sp), INTENT(IN), TARGET :: host(*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_real_sp_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_1d_sp

  SUBROUTINE acc_map_data_1d_int(host, device, data_size)
    INTEGER, INTENT(IN), TARGET :: host(*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_int_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_1d_int

  SUBROUTINE acc_map_data_2d_dp(host, device, data_size)
    REAL(dp), INTENT(IN), TARGET :: host(1,*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_real_dp_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_2d_dp

  SUBROUTINE acc_map_data_2d_sp(host, device, data_size)
    REAL(sp), INTENT(IN), TARGET :: host(1,*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_real_sp_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_2d_sp

  SUBROUTINE acc_map_data_2d_int(host, device, data_size)
    INTEGER, INTENT(IN), TARGET :: host(1,*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_int_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_2d_int

  SUBROUTINE acc_map_data_3d_dp(host, device, data_size)
    REAL(dp), INTENT(IN), TARGET :: host(1,1,*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_real_dp_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_3d_dp

  SUBROUTINE acc_map_data_3d_sp(host, device, data_size)
    REAL(sp), INTENT(IN), TARGET :: host(1,1,*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_real_sp_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_3d_sp

  SUBROUTINE acc_map_data_3d_int(host, device, data_size)
    INTEGER, INTENT(IN), TARGET :: host(1,1,*)
    TYPE(c_ptr), INTENT(IN) :: device
    INTEGER, INTENT(IN) :: data_size
    CALL c_acc_map_data( &
      c_loc(host), device, INT(p_int_byte * data_size, c_size_t))
  END SUBROUTINE acc_map_data_3d_int

  SUBROUTINE acc_unmap_data_c_ptr(host)
    TYPE(c_ptr), INTENT(IN) :: host
    CALL c_acc_unmap_data(host)
  END SUBROUTINE acc_unmap_data_c_ptr

  SUBROUTINE acc_unmap_data_1d_dp(host)
    REAL(dp), INTENT(IN), TARGET :: host(*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_1d_dp

  SUBROUTINE acc_unmap_data_1d_sp(host)
    REAL(sp), INTENT(IN), TARGET :: host(*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_1d_sp

  SUBROUTINE acc_unmap_data_1d_int(host)
    INTEGER, INTENT(IN), TARGET :: host(*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_1d_int

  SUBROUTINE acc_unmap_data_2d_dp(host)
    REAL(dp), INTENT(IN), TARGET :: host(1,*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_2d_dp

  SUBROUTINE acc_unmap_data_2d_sp(host)
    REAL(sp), INTENT(IN), TARGET :: host(1,*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_2d_sp

  SUBROUTINE acc_unmap_data_2d_int(host)
    INTEGER, INTENT(IN), TARGET :: host(1,*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_2d_int

  SUBROUTINE acc_unmap_data_3d_dp(host)
    REAL(dp), INTENT(IN), TARGET :: host(1,1,*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_3d_dp

  SUBROUTINE acc_unmap_data_3d_sp(host)
    REAL(sp), INTENT(IN), TARGET :: host(1,1,*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_3d_sp

  SUBROUTINE acc_unmap_data_3d_int(host)
    INTEGER, INTENT(IN), TARGET :: host(1,1,*)
    CALL c_acc_unmap_data(c_loc(host))
  END SUBROUTINE acc_unmap_data_3d_int
#endif

END MODULE mo_openacc
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
