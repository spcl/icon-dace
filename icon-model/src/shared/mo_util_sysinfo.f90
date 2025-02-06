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

MODULE mo_util_sysinfo

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_LONG
  USE mo_mpi,            ONLY: p_mpi_wtime
  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: find_next_free_unit
  USE mo_exception,      ONLY: finish

  IMPLICIT NONE

  PRIVATE

  ! module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_sysinfo'


  ! log stream errors:
  ENUM, BIND(C)
    ENUMERATOR ::                 &
    SUCCESS                 =  0, &
    ERROR_FILE_NOT_FOUND    = -2, &
    ERROR_FILE_NOT_READABLE = -3 
  END ENUM


  INTERFACE
    SUBROUTINE private_util_user_name(name, name_len) BIND(C,NAME='util_user_name') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_user_name
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_os_system(name, name_len) BIND(C,NAME='util_os_system') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_os_system
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_node_name(name, name_len) BIND(C,NAME='util_node_name') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_node_name
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_get_maxrss(maxrss) BIND(C,NAME='util_get_maxrss') 
      IMPORT :: C_INT
      INTEGER(C_INT), INTENT(out) :: maxrss
    END SUBROUTINE private_util_get_maxrss
  END INTERFACE


  INTERFACE
    SUBROUTINE private_util_compiler_release(release_str, rstr_len) BIND(C,NAME='util_compiler_release') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: release_str
      INTEGER(C_INT), INTENT(out) :: rstr_len
    END SUBROUTINE private_util_compiler_release
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_c_getpid(pid) BIND(C,NAME='util_c_getpid') 
      IMPORT :: C_LONG
      INTEGER(C_LONG), INTENT(out) :: pid
    END SUBROUTINE private_util_c_getpid
  END INTERFACE


  PUBLIC :: util_user_name
  PUBLIC :: util_os_system
  PUBLIC :: util_node_name
  PUBLIC :: util_get_maxrss
  PUBLIC :: util_compiler_release
  PUBLIC :: util_c_getpid
  PUBLIC :: get_smaps_sum
  PUBLIC :: check_file_exists

CONTAINS

  SUBROUTINE util_user_name(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_user_name(name, name_len)
  END SUBROUTINE util_user_name

  SUBROUTINE util_os_system(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_os_system(name, name_len)
  END SUBROUTINE util_os_system

  SUBROUTINE util_node_name(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_node_name(name, name_len)
  END SUBROUTINE util_node_name

  SUBROUTINE util_get_maxrss(maxrss)
    INTEGER, INTENT(out) :: maxrss
    CALL private_util_get_maxrss(maxrss)
  END SUBROUTINE util_get_maxrss

  SUBROUTINE util_compiler_release(release_str, rstr_len)
    CHARACTER(len=*), INTENT(out) :: release_str
    INTEGER, INTENT(out) :: rstr_len
    CALL private_util_compiler_release(release_str, rstr_len)
  END SUBROUTINE util_compiler_release

  SUBROUTINE util_c_getpid(pid)
    INTEGER(c_long), INTENT(out) :: pid
    CALL private_util_c_getpid(pid)
  END SUBROUTINE util_c_getpid


  ! auxiliary function to measure huge page consumption on XCE.
  !
  ! This should do the same as D. Sternkopf's suggestion:
  !
  ! grep -B 11 'KernelPageSize:     2048 kB' /proc/$pid/smaps |  \
  !       grep "^Size:" | awk 'BEGIN{sum=0}{sum+=$2}END{print sum/1024}'  
  !
  INTEGER FUNCTION get_smaps_sum(filename, opt_ierr)
    CHARACTER(LEN=*),   INTENT(IN)    :: filename   !< source file name.
    INTEGER, INTENT(OUT), OPTIONAL    :: opt_ierr   !< error code (0=SUCCESS).
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_smaps_sum"
#ifndef _CRAYFTN
    CALL finish("get_smaps_sum", "Not implemented!")
#else
    ! local variables
    INTEGER, PARAMETER               :: BUFFER_LEN = 512
    CHARACTER(LEN=*), PARAMETER      :: KEYWORD1   = "Size:"
    CHARACTER(LEN=*), PARAMETER      :: KEYWORD2   = "KernelPageSize"
    INTEGER,          PARAMETER      :: KEYVAL2    = 2048

    CHARACTER(len=BUFFER_LEN)        :: buffer, str1, str2
    LOGICAL :: l_exist
    INTEGER :: io_error, iunit, ierr, ival, sum, ival2

    ierr = SUCCESS
    INQUIRE (FILE=filename, EXIST=l_exist)
    IF (.NOT.l_exist) THEN
      ierr = ERROR_FILE_NOT_FOUND
      IF (PRESENT(opt_ierr))  opt_ierr = ierr
      RETURN
    END IF

    ! print input file's contents:
    iunit = find_next_free_unit(10,100)
    OPEN(unit=iunit, file=filename, status='old',action='read', iostat=io_error) 
    IF ( io_error /= 0) THEN
      ierr = ERROR_FILE_NOT_READABLE
      IF (PRESENT(opt_ierr))  opt_ierr = ierr
      RETURN
    END IF

    sum = 0
    DO
      READ (iunit, '(A)', iostat=io_error) buffer
      IF (io_error /= 0)  EXIT

      ! parse lines with trigger word KEYWORD1
      IF (buffer(1:LEN(KEYWORD1)) == KEYWORD1) THEN
        READ (buffer, *) str1, ival, str2
      END IF

      ! add last value if trigger word KEYWORD2 occurs:
      IF (buffer(1:LEN(KEYWORD2)) == KEYWORD2) THEN
        READ (buffer, *) str1, ival2, str2
        IF (ival2 == KEYVAL2) THEN
          sum = sum + ival
        END IF
      END IF
    END DO

    CLOSE(iunit) 
    IF (PRESENT(opt_ierr))  opt_ierr = ierr
    get_smaps_sum = sum/1024
#endif
  END FUNCTION get_smaps_sum


  !> Idle wait for @p wait_sec seconds. Non-collective.
  !
  SUBROUTINE wait_idle(wait_sec)
    INTEGER, INTENT(IN) :: wait_sec
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::wait_idle"
#ifdef NOMPI
    CALL finish(routine, "Only implemented for parallel setups!")
#else
    ! local variables
    REAL(wp) starttime, elapsed_sec
    starttime = p_mpi_wtime()
    LOOP: DO
      ! Polling loop: Do not use in regions of spacetime exhibiting
      ! strong gravitational effects, right?
      elapsed_sec = p_mpi_wtime() - starttime
      IF (INT(elapsed_sec) > wait_sec)  EXIT LOOP
    END DO LOOP
#endif
  END SUBROUTINE wait_idle


  !> Checks with INQUIRE, if a file with the (full path) file name @p
  !  filename exists. Optionally: If the file does not exist, perform
  !  an idle-wait-and-retry loop.
  !
  !  @return .TRUE. if file exists.
  !
  LOGICAL FUNCTION check_file_exists(filename, nretries, retry_wait_sec)
    CHARACTER(LEN=*),  INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: nretries
    INTEGER, INTENT(IN), OPTIONAL :: retry_wait_sec
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::check_file_exists"
    INTEGER :: i
    LOGICAL :: l_exist

    ! consistency check
    IF (PRESENT(retry_wait_sec) .NEQV. PRESENT(nretries)) THEN
      CALL finish(routine, "Internal error!")
    END IF

    INQUIRE(FILE=TRIM(ADJUSTL(filename)), EXIST=l_exist)
    IF (PRESENT(nretries) .AND. .NOT. l_exist) THEN
      LOOP : DO i=1,nretries
        WRITE (0,'(4a,i0,a)') routine, ' :: file not found: "', &
          &                   TRIM(ADJUSTL(filename)), '", retry after ', retry_wait_sec, ' seconds.'
        CALL wait_idle(retry_wait_sec)
        INQUIRE(FILE=TRIM(ADJUSTL(filename)), EXIST=l_exist)
        IF (l_exist)  EXIT LOOP
      END DO LOOP
    END IF
    check_file_exists = l_exist
  END FUNCTION check_file_exists

END MODULE mo_util_sysinfo
