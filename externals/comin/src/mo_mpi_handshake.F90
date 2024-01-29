!> @file mo_mpi_handshake.F90
!! @brief Interfaces for MPI handshake routines (written in C).
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE mo_mpi_handshake


  PUBLIC MAX_GROUPNAME_LEN, mpi_handshake, mpi_handshake_dummy
  INTEGER, PARAMETER :: MAX_GROUPNAME_LEN = 256

CONTAINS

  !> Procedure for the communicator splitting ("MPI handshake") that
  !> has been harmonized with the respective algorithm of the YAC
  !> coupler software
  !! @ingroup host_interface
  SUBROUTINE mpi_handshake ( comm, group_names, group_comms )
    USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_char, c_null_char, c_loc

    implicit none

    interface
       SUBROUTINE mpi_handshake_c2f (n, group_names, group_comms, comm) &
          bind ( c, name='mpi_handshake_c2f')
         use, intrinsic :: iso_c_binding, only : c_ptr, c_int
         implicit none
         INTEGER(KIND=c_int), INTENT(in), VALUE :: n
         TYPE (c_ptr) , INTENT(in) :: group_names(n)
         INTEGER(KIND=c_int), INTENT(out) :: group_comms(n)
         INTEGER(KIND=c_int), INTENT(in), VALUE :: comm
       end subroutine mpi_handshake_c2f
    end interface

    integer, intent(in) :: comm
    character(len=MAX_GROUPNAME_LEN), intent(in) :: group_names(:)
    integer, intent(inout) :: group_comms(SIZE(group_names))

    CHARACTER (kind=c_char, len=MAX_GROUPNAME_LEN), TARGET :: group_names_cpy(SIZE(group_names))
    type( c_ptr ) :: group_name_ptr(SIZE(group_names))
    integer :: i
    DO i=1,SIZE(group_names)
       group_names_cpy(i) = TRIM(group_names(i)) // c_null_char
       group_name_ptr(i) = c_loc(group_names_cpy(i))
    END DO

    CALL mpi_handshake_c2f(SIZE(group_names), group_name_ptr, group_comms, comm)

  END SUBROUTINE mpi_handshake

    !! @ingroup host_interface
  SUBROUTINE mpi_handshake_dummy(comm)
    INTEGER, INTENT(IN) :: comm
    INTEGER :: empty_int_array(0)
    CHARACTER(LEN=MAX_GROUPNAME_LEN) :: empty_char_array(0)

    CALL mpi_handshake(comm, empty_char_array, empty_int_array)
  END SUBROUTINE mpi_handshake_dummy

END MODULE mo_mpi_handshake

