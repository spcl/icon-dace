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

MODULE mo_util_file

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char, c_null_char, c_long
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: i8 => int64
  USE mo_exception, ONLY: finish
  USE mo_io_units, ONLY: filename_max

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    FUNCTION private_symlink(file, link) RESULT(iret) BIND(C, NAME='symlink')
      IMPORT :: c_int, c_char
      INTEGER(c_int) :: iret
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: file
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: link
    END FUNCTION private_symlink
  END INTERFACE

  INTERFACE
    FUNCTION private_unlink(filename) RESULT(iret) BIND(C, NAME='unlink')
      IMPORT :: c_int, c_char
      INTEGER(c_int) :: iret
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: filename
    END FUNCTION private_unlink
  END INTERFACE

  INTERFACE
    FUNCTION private_islink(filename) RESULT(iret) BIND(C, NAME='util_islink')
      IMPORT :: c_int, c_char
      INTEGER(c_int) :: iret
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: filename
    END FUNCTION private_islink
  END INTERFACE

  INTERFACE
    FUNCTION private_rename(old_filename, new_filename) RESULT(iret) BIND(C, NAME='rename')
      IMPORT :: c_int, c_char
      INTEGER(c_int) :: iret
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: old_filename
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: new_filename
    END FUNCTION private_rename
  END INTERFACE

  INTERFACE
    FUNCTION private_create_tmpfile(filename, max_len) RESULT(flen) BIND(C, NAME='util_create_tmpfile')
      IMPORT :: c_char, c_int
      INTEGER(c_int) :: flen
      CHARACTER(c_char), DIMENSION(*), INTENT(INOUT) :: filename
      INTEGER(c_int), VALUE, INTENT(IN) :: max_len
    END FUNCTION private_create_tmpfile
  END INTERFACE

  INTERFACE
    FUNCTION private_filesize(filename) RESULT(flen) BIND(C, NAME='util_filesize')
      IMPORT :: c_long, c_char
      INTEGER(c_long) :: flen
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: filename
    END FUNCTION private_filesize
  END INTERFACE

  INTERFACE
    FUNCTION private_file_is_writable(filename) RESULT(iwritable) BIND(C, NAME='util_file_is_writable')
      IMPORT :: c_int, c_char
      INTEGER(c_int) :: iwritable
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: filename
    END FUNCTION private_file_is_writable
  END INTERFACE

  PUBLIC :: util_symlink
  PUBLIC :: util_unlink
  PUBLIC :: util_islink
  PUBLIC :: util_rename
  PUBLIC :: util_tmpnam
  PUBLIC :: util_filesize
  PUBLIC :: util_file_is_writable
  PUBLIC :: createSymlink
  PUBLIC :: get_filename
  PUBLIC :: get_filename_noext
  PUBLIC :: get_path

  CHARACTER(*), PARAMETER :: modname = "mo_util_file"

CONTAINS

  FUNCTION util_symlink(file, link) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(IN) :: file
    CHARACTER(len=*), INTENT(IN) :: link
    iret = private_symlink(TRIM(file)//c_null_char, TRIM(link)//c_null_char)
  END FUNCTION util_symlink

  FUNCTION util_unlink(filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(IN) :: filename
    iret = private_unlink(TRIM(filename)//c_null_char)
  END FUNCTION util_unlink

  FUNCTION util_islink(filename) RESULT(islink)
    LOGICAL :: islink
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER :: iret
    iret = private_islink(TRIM(filename)//c_null_char)
    islink = .FALSE.
    IF (iret == 1) islink = .TRUE.
  END FUNCTION util_islink

  FUNCTION util_rename(old_filename, new_filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(IN) :: old_filename
    CHARACTER(len=*), INTENT(IN) :: new_filename
    iret = private_rename(TRIM(old_filename)//c_null_char, TRIM(new_filename)//c_null_char)
  END FUNCTION util_rename

  FUNCTION create_tmpfile(filename, max_len) RESULT(flen)
    INTEGER :: flen
    CHARACTER(len=*), INTENT(OUT) :: filename
    INTEGER, INTENT(IN)  :: max_len
    ! local variables
    INTEGER :: i
    !
    CHARACTER(c_char), ALLOCATABLE :: c_filename(:)
    !
    ALLOCATE (c_filename(max_len + 1))
    flen = private_create_tmpfile(c_filename, max_len + 1) - 1
    DO i = 1, flen
      filename(i:i) = c_filename(i)
    END DO
    DEALLOCATE (c_filename)
  END FUNCTION create_tmpfile

  FUNCTION util_tmpnam(filename, klen) RESULT(flen)
    INTEGER :: flen
    CHARACTER(len=*), INTENT(OUT) :: filename
    INTEGER, INTENT(IN)  :: klen

    flen = create_tmpfile(filename, klen)
    IF (flen < 0) THEN
      CALL finish('mo_util_file::util_tmpnam', 'Failed to find a tmp filename!')
    END IF
  END FUNCTION util_tmpnam

  FUNCTION util_filesize(filename) RESULT(flen)
    INTEGER(KIND=i8) :: flen
    CHARACTER(len=*), INTENT(IN) :: filename
    flen = private_filesize(TRIM(filename)//c_null_char)
  END FUNCTION util_filesize

  FUNCTION util_file_is_writable(filename) RESULT(lwritable)
    LOGICAL :: lwritable
    CHARACTER(len=*), INTENT(IN) :: filename
    lwritable = (private_file_is_writable(TRIM(filename)//c_null_char) == 1)
  END FUNCTION util_file_is_writable

  INTEGER FUNCTION createSymlink(targetPath, linkName) RESULT(error)
    CHARACTER(*), INTENT(IN) :: targetPath, linkName
    INTEGER :: i
    CHARACTER(KIND=c_char) :: linkNameCopy(LEN(linkName) + 1), targetPathCopy(LEN(targetPath) + 1)
    INTERFACE
      INTEGER(c_int) FUNCTION c_createSymlink(c_targetPath, c_linkName) BIND(C, NAME="createSymlink")
        IMPORT c_int, c_char
        CHARACTER(KIND=c_char) :: c_targetPath(*), c_linkName(*)
      END FUNCTION c_createSymlink
    END INTERFACE

    DO i = 1, LEN(targetPath)
      targetPathCopy(i) = targetPath(i:i)
    END DO
    targetPathCopy(i) = c_null_char
    DO i = 1, LEN(linkName)
      linkNameCopy(i) = linkName(i:i)
    END DO
    linkNameCopy(i) = c_null_char
    error = c_createSymlink(targetPathCopy, linkNameCopy)
  END FUNCTION createSymlink

  FUNCTION get_filename(in_str) RESULT(out_str)
    ! Parameters
    CHARACTER(len=*), INTENT(IN) :: in_str
    CHARACTER(len=LEN(in_str))   :: out_str
    ! Local parameters
    INTEGER :: pos

    ! start after last occurrence of '/':
    pos = INDEX(in_str, '/', .TRUE.)
    pos = pos + 1
    out_str = in_str(pos:LEN(in_str))
  END FUNCTION get_filename

  FUNCTION get_filename_noext(name) RESULT(filename)
    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(LEN=filename_max) :: filename

    INTEGER :: end_name

    end_name = INDEX(name, '.', back=.TRUE.)
    IF (end_name > 0) THEN
      filename = name(1:end_name - 1)
    ELSE
      filename = TRIM(name)
    END IF

  END FUNCTION get_filename_noext

  FUNCTION get_path(in_str) RESULT(out_str)
    ! Parameters
    CHARACTER(len=*), INTENT(IN) :: in_str
    CHARACTER(len=LEN(in_str))   :: out_str
    ! Local parameters
    INTEGER :: pos

    ! start after last occurrence of '/':
    out_str = ""
    pos = INDEX(in_str, '/', .TRUE.)
    IF (pos > 0) THEN
      out_str = in_str(1:pos)
    END IF
  END FUNCTION get_path

END MODULE mo_util_file
