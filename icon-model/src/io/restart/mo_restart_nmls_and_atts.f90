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

MODULE mo_restart_nml_and_att

  USE mo_exception,             ONLY: finish, message
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_master_config,         ONLY: isRestart
  USE mo_mpi, ONLY: p_comm_rank
  USE mo_util_file, ONLY: util_tmpnam, util_filesize, util_unlink
  USE mo_util_string, ONLY: tocompact
  USE mo_io_units, ONLY: filename_max, find_next_free_unit
  USE mo_key_value_store, ONLY: t_key_value_store
  USE mo_kind, ONLY: wp
  USE mo_hash_table, ONLY: t_HashIterator
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf

  IMPLICIT NONE
  PRIVATE

  ! used throughout the icon code
  PUBLIC :: open_tmpfile, store_and_close_namelist
  PUBLIC :: open_and_restore_namelist, close_tmpfile

  PUBLIC :: restartAttributeList_make, restartAttributeList_write_to_ncdf
  PUBLIC :: getAttributesForRestarting, restartAttributeList_read

  PUBLIC :: bcastNamelistStore

  TYPE(t_key_value_store), TARGET :: gAttributeStore
  TYPE(t_key_value_store) :: gNamelistStore

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_nml_and_att"
  LOGICAL, SAVE, PUBLIC :: ocean_initFromRestart_OVERRIDE = .FALSE. 

CONTAINS

  INTEGER FUNCTION open_tmpfile() RESULT(funit)
    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename

    flen = util_tmpnam(filename, filename_max)
    funit = find_next_free_unit(10,100)
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='sequential', DELIM='apostrophe')
  END FUNCTION open_tmpfile

  SUBROUTINE store_and_close_namelist(funit, nmlname)
    INTEGER, INTENT(IN) :: funit
    CHARACTER(*), INTENT(in) :: nmlname
    CHARACTER(filename_max) :: filename
    INTEGER :: nmllen, ierr
    CHARACTER(:), ALLOCATABLE :: nmlbuf
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":store_and_close_namelist"

    INQUIRE(funit, NAME=filename)
    CLOSE(funit)
    nmllen = INT(util_filesize(filename))
    IF (nmllen == 0) CALL message(routine, 'empty nml "'//TRIM(nmlname)//'" : not saved')
    ALLOCATE(CHARACTER(LEN = nmllen) :: nmlbuf, STAT=ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation failed")
    nmlbuf(1:nmllen) = ''
    OPEN(UNIT=funit, FILE=TRIM(filename), ACTION='read', ACCESS='stream', FORM='unformatted')
    READ(funit) nmlbuf(1:nmllen)
    CLOSE(funit)
    CALL tocompact(nmlbuf)
    IF(.NOT.gNamelistStore%is_init) CALL gNamelistStore%init(.true.)
    CALL gNamelistStore%put('nml_'//nmlname, TRIM(nmlbuf))
    ierr = util_unlink(TRIM(filename))
  END SUBROUTINE store_and_close_namelist

  INTEGER FUNCTION open_and_restore_namelist(nmlname) RESULT(funit)
    CHARACTER(len=*), INTENT(in) :: nmlname
    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename
    CHARACTER(:), ALLOCATABLE :: nmlbuf

    IF (.NOT.gNamelistStore%is_init) &
      & CALL finish(modname//':open_and_restore_namelist', "gNamelistStore not initialized")
    CALL gNamelistStore%get('nml_'//nmlname, nmlbuf)
    flen = util_tmpnam(filename, filename_max)
    funit = find_next_free_unit(10, 100)
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='stream', FORM='unformatted')
    WRITE(funit) TRIM(nmlbuf)
    CLOSE(funit)
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='read', ACCESS='sequential', RECL=65535, DELIM='apostrophe')
  END FUNCTION open_and_restore_namelist

  SUBROUTINE close_tmpfile(funit)
    INTEGER, INTENT(in) :: funit
    CHARACTER(len=filename_max) :: filename
    INTEGER :: ierr
    CHARACTER(*), PARAMETER :: routine = modname//":close_tmpfile"

    INQUIRE(funit,NAME=filename)
    CLOSE(funit)
    ierr = util_unlink(TRIM(filename))
    IF(ierr /= SUCCESS) CALL message(routine, "warning: error unlinking '"//TRIM(filename)//"'")
  END SUBROUTINE close_tmpfile

  SUBROUTINE getAttributesForRestarting(ptr)
    TYPE(t_key_value_store), POINTER, INTENT(OUT) :: ptr
    CHARACTER(*), PARAMETER :: routine = modname//":getAttributesForRestarting"

    IF (.NOT. gAttributeStore%is_init .AND. &
      & (isRestart() .OR. ocean_initFromRestart_OVERRIDE)) &
      & CALL finish(routine, "restart attributes not yet set up")
    ptr => gAttributeStore
  END SUBROUTINE getAttributesForRestarting

  SUBROUTINE restartAttributeList_read(root_pe, comm, ncid)
    INTEGER, INTENT(IN) :: root_pe, comm, ncid
    CHARACTER(*), PARAMETER :: routine = modname//":restartAttributeList_read"
    INTEGER :: natt, i, aType, aLen, namLen, iTmp(1)
    CHARACTER(NF90_MAX_NAME) :: aName
    REAL(wp) :: dTmp(1)
    CHARACTER(:), ALLOCATABLE :: aTxt

    IF (.NOT.isRestart() .AND. .NOT.ocean_initFromRestart_OVERRIDE) &
      & CALL finish(routine, "not a restart run")
    IF (gAttributeStore%is_init) &
      & CALL finish(routine, "no second assignment of gAttributeStore allowed")
    IF (p_comm_rank(comm) == root_pe) THEN
      CALL gAttributeStore%init(.true.)
      ALLOCATE(CHARACTER(LEN=4096) :: aTxt)
      CALL nf(nf90_inquire(ncid, nAttributes = natt), routine)
      DO i = 1, natt
        CALL nf(nf90_inq_attname(ncid, NF90_GLOBAL, i, aName), routine)
        namLen = LEN_TRIM(aName)
        CALL nf(nf90_inquire_attribute(ncid, NF90_GLOBAL, aName(1:namLen), xtype = aType, len = aLen), routine)
        SELECT CASE(aType)
        CASE(NF90_DOUBLE)
          CALL nf(nf90_get_att(ncid, NF90_GLOBAL, aName(1:namLen), dTmp), routine)
          CALL gAttributeStore%put(aName(1:namLen), dTmp(1))
        CASE(NF90_INT)
          CALL nf(nf90_get_att(ncid, NF90_GLOBAL, aName(1:namLen), iTmp), routine)
          IF ('bool_' == aName(1:MIN(5,namLen))) THEN
            CALL gAttributeStore%put(aName(6:namLen), iTmp(1) .NE. 0)
          ELSE
            CALL gAttributeStore%put(aName(1:namLen), iTmp(1))
          END IF
        CASE(NF90_CHAR)
          IF (aLen .GT. LEN(aTxt)) DEALLOCATE(aTxt)
          IF (.NOT.ALLOCATED(aTxt)) ALLOCATE(CHARACTER(LEN=aLen) :: aTxt)
          CALL nf(nf90_get_att(ncid, NF90_GLOBAL, aName(1:namLen), aTxt(1:aLen)), routine)
          CALL gAttributeStore%put(aName(1:namLen), aTxt(1:aLen))
        END SELECT
      ENDDO
    END IF
    CALL gAttributeStore%bcast(root_pe, comm)
    IF (.NOT.ocean_initFromRestart_OVERRIDE) &
      CALL gAttributeStore%output(archive=gNamelistStore)
  END SUBROUTINE restartAttributeList_read

  SUBROUTINE restartAttributeList_write_to_ncdf(kvs, ncid)
#ifdef __PGI
    USE mo_util_texthash, ONLY: t_char_workaround
#endif
    TYPE(t_key_value_store), INTENT(IN) :: kvs
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), PARAMETER :: routine = modname//"write_to_ncdf"
    TYPE(t_HashIterator), POINTER :: iterator
    CLASS(*), POINTER :: curKey, curVal
    CHARACTER(:), POINTER :: ccKey
    CHARACTER(:), ALLOCATABLE :: ccTmp

    iterator => kvs%getIterator()
    DO WHILE (iterator%nextEntry(curKey, curVal))
      SELECT TYPE(curKey)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        ccKey => curKey%c
#endif
      TYPE IS(CHARACTER(*))
        ccKey => curKey
      CLASS DEFAULT
        CALL finish(routine, "key: invalid type")
      END SELECT
      SELECT TYPE(curVal)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        CALL nf(nf90_put_att(ncid, NF90_GLOBAL, ccKey, curVal%c), routine)
#else
      TYPE IS(CHARACTER(*))
        CALL nf(nf90_put_att(ncid, NF90_GLOBAL, ccKey, curVal), routine)
#endif
      TYPE IS(REAL(wp))
        CALL nf(nf90_put_att(ncid, NF90_GLOBAL, ccKey, curVal), routine)
      TYPE IS(INTEGER)
        CALL nf(nf90_put_att(ncid, NF90_GLOBAL, ccKey, curVal), routine)
      TYPE IS(LOGICAL)
        ccTmp = 'bool_'//ccKey
        CALL nf(nf90x_put_att_converted(ncid, NF90_GLOBAL, ccTmp, curVal), routine)
      CLASS DEFAULT
        CALL finish(routine, "val: invalid type")
      END SELECT
    END DO
    DEALLOCATE(iterator)
  END SUBROUTINE restartAttributeList_write_to_ncdf

  SUBROUTINE restartAttributeList_make(ptr)
    TYPE(t_key_value_store), ALLOCATABLE, INTENT(OUT) :: ptr

    ALLOCATE(ptr)
    CALL gNamelistStore%output(copy=ptr)
  END SUBROUTINE restartAttributeList_make

  SUBROUTINE bcastNamelistStore(root, comm)
    INTEGER, INTENT(IN) :: root, comm
    
    CALL gNamelistStore%bcast(root, comm)
  END SUBROUTINE bcastNamelistStore

END MODULE mo_restart_nml_and_att
