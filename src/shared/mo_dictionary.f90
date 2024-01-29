! Module for associating pairs of strings,
! e.g. for translating between variable names.
!
! @note This implementation is aimed at *small* dictionaries,
!       e.g. translation between variable names. The algorithms
!       contained in this module are therefore simple and will
!       not scale with large data bases.
!
! @note When loading key-value pairs from an external text file, the
!       current implementation is restricted to strings that do not
!       contain spaces.
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

MODULE mo_dictionary

  USE mo_exception,      ONLY: finish
  USE mo_util_string,    ONLY: tocompact
  USE mo_io_units,       ONLY: find_next_free_unit
  USE mo_key_value_store, ONLY: t_key_value_store

  IMPLICIT NONE
  PRIVATE

  CHARACTER(*), PARAMETER :: modname = 'mo_dictionary'
  INTEGER, PARAMETER :: DICT_MAX_STRLEN = 132    !< maximum string length

  TYPE t_dictionary
    TYPE(t_key_value_store), PRIVATE :: dic, idic
  CONTAINS
    PROCEDURE :: init => dict_init
    PROCEDURE :: finalize => dict_finalize
    PROCEDURE :: get => dict_get
    PROCEDURE :: loadfile => dict_loadfile
    PROCEDURE :: copy => dict_copy
    PROCEDURE :: bcast => dict_bcast
  END TYPE t_dictionary

  PUBLIC :: t_dictionary, DICT_MAX_STRLEN
 
CONTAINS

  SUBROUTINE dict_bcast(dict, root, comm)
    CLASS(t_dictionary), INTENT(INOUT) :: dict
    INTEGER, INTENT(IN) :: root, comm

    CALL dict%dic%bcast(root, comm)
    CALL dict%dic%output(inverse=dict%idic)
  END SUBROUTINE dict_bcast

  SUBROUTINE dict_init(dict, lcase_sensitive)
    CLASS(t_dictionary), INTENT(INOUT) :: dict            !< dictionary data structure
    LOGICAL,            INTENT(IN)    :: lcase_sensitive !< Flag. If .TRUE. key strings are handled case-sensitively

    CALL dict%dic%init(lcase_sensitive)
    CALL dict%idic%init(lcase_sensitive)
  END SUBROUTINE dict_init

  SUBROUTINE dict_finalize(dict)
    CLASS(t_dictionary), INTENT(INOUT) :: dict !< dictionary data structure

    CALL dict%dic%destruct()
    CALL dict%idic%destruct()
  END SUBROUTINE dict_finalize

  SUBROUTINE dict_copy(src_dict, dst_dict)
    CLASS(t_dictionary), INTENT(IN)    :: src_dict !< source dictionary data structure
    TYPE(t_dictionary), INTENT(INOUT) :: dst_dict !< destination dictionary data structure

    CALL src_dict%dic%output(copy=dst_dict%dic)
    CALL src_dict%idic%output(copy=dst_dict%idic)
  END SUBROUTINE dict_copy

  !--------------------------------------------------------------------------
  !> Retrieve a value from the given dictionary which is assciated
  !  with @p key.
  !
  !  @param[in] linverse   optional flag: If .TRUE., the dictionary is queried in inverse order.
  !
  !  @return If the key is not found, this function returns the
  !          "default" value. If this optional argument has not been
  !          provided, this function throws an error message.
  !
  FUNCTION dict_get(dict, key, default, linverse) RESULT(val)
    CHARACTER(LEN=DICT_MAX_STRLEN) :: val !< return value.
    CLASS(t_dictionary), INTENT(IN), TARGET :: dict !< dictionary data structure
    CHARACTER(*), INTENT(IN) :: key  !< new search key
    CHARACTER(*), INTENT(IN), OPTIONAL :: default  !< default value
    LOGICAL, INTENT(IN), OPTIONAL :: linverse !< Flag. If .TRUE., the dictionary is queried in inverse order.
    CHARACTER(*), PARAMETER :: routine = modname//':dict_get'
    CHARACTER(:), ALLOCATABLE :: tmp
    INTEGER :: opt_err
    LOGICAL :: linv
    TYPE(t_key_value_store), POINTER :: dic_ptr

    linv = .FALSE.
    IF (PRESENT(linverse)) linv = linverse
    dic_ptr => dict%dic
    IF (linv) dic_ptr => dict%idic
    CALL dic_ptr%get(key, tmp, opt_err=opt_err)
    IF (opt_err .EQ. 0) THEN
      val = tmp
    ELSE IF (PRESENT(default)) THEN
      val = TRIM(default)
    ELSE
      CALL finish(routine, 'Requested dictionary key ' // TRIM(key) // ' not found!')
    END IF
  END FUNCTION dict_get

  !--------------------------------------------------------------------------
  !> Load the contents of a text file into the given dictionary data
  !  structure.
  ! 
  !  @note The file format is assumed to be as follows: Each line of
  !        the text file contains a key-value pair, where both strings
  !        are separated by one or more spaces. Comment lines
  !        (beginning with "#") are ignored.
  !
  SUBROUTINE dict_loadfile(dict, filename, linverse)
    CLASS(t_dictionary), INTENT(INOUT) :: dict      !< dictionary data structure
    CHARACTER(*),   INTENT(IN)    :: filename  !< text file name
    LOGICAL, OPTIONAL,  INTENT(IN)    :: linverse  !< Flag. If .TRUE., dictionary columns are read in inverse order.
    CHARACTER(*), PARAMETER :: routine = modname//':dict_loadfile'
    INTEGER :: iunit, ist, llen, klen
    CHARACTER(LEN=256), TARGET :: line
! Cray Fortran Compiler < 8.7.11 (tested with 8.7.7) does not interpret the strings put into the dictionary correctly.
! Hence, for older cray compiler versions, the following code circumvents this bug. 
! This can be removed when recent cray compiler versions are available at ECMWF's Cray
#if (_CRAYFTN == 1) && ( _RELEASE_MAJOR <= 8)
    CHARACTER(LEN=DICT_MAX_STRLEN) :: key, val
#else
    CHARACTER(:), POINTER :: key, val
#endif
    LOGICAL :: lread_inverse

    lread_inverse = .FALSE.
    IF (PRESENT(linverse)) lread_inverse = linverse
    iunit = find_next_free_unit(10,99)
    OPEN (unit=iunit,file=filename,access='SEQUENTIAL', &
      &  form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)
    IF(ist/=0) CALL finish (routine, 'Open of dictionary file '//TRIM(filename)//' failed')
    DO
      line = ''
      READ(iunit,'(a)',iostat=ist) line
      IF(ist < 0) EXIT ! No more lines
      IF(ist > 0) CALL finish(routine, 'Read error in dictionary file '//TRIM(filename))
      CALL tocompact(line)
      llen = LEN_TRIM(line)
      IF(llen .EQ. 0) CYCLE ! blank line
      IF(line(1:1) == '#') CYCLE ! comment
      klen = INDEX(line, ' ')
      IF (llen .LE. klen) CALL finish(routine, "Illegal line in name_map.")
      IF (klen-1 .GT. DICT_MAX_STRLEN) CALL finish(routine, "Key too long")
      IF (llen - klen - 1 .GT. DICT_MAX_STRLEN) CALL finish(routine, "Value too long.")
! Cray Fortran Compiler < 8.7.11 (tested with 8.7.7) does not interpret the strings put into the dictionary correctly.
! Hence, for older cray compiler versions, the following code circumvents this bug.
! This can be removed when recent cray compiler versions are available at ECMWF's Cray
#if (_CRAYFTN == 1) && ( _RELEASE_MAJOR <= 8)
      key = line(1:klen-1)
      val = line(klen+1:llen)
      IF (.NOT. lread_inverse) THEN
        CALL dict%dic%put(TRIM(key), TRIM(val))
        CALL dict%idic%put(TRIM(val), TRIM(key))
      ELSE
        CALL dict%dic%put(TRIM(val), TRIM(key))
        CALL dict%idic%put(TRIM(key), TRIM(val))
      END IF
#else
      key => line(1:klen-1)
      val => line(klen+1:llen)
      IF (.NOT. lread_inverse) THEN
        CALL dict%dic%put(key, val)
        CALL dict%idic%put(val, key)
      ELSE
        CALL dict%dic%put(val, key)
        CALL dict%idic%put(key, val)
      END IF
#endif
    ENDDO
    CLOSE(unit=iunit)
  END SUBROUTINE dict_loadfile

END MODULE mo_dictionary
