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

! Auxiliary module for status information of the version control system (VCS).
!
! Most of the routines which are defined within this module are ISO-C bindings
! to a related C program "version.c", which is automatically generated at the
! build time.

MODULE mo_util_vcs

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    C_ASSOCIATED, C_F_POINTER, C_LOC, c_char, c_int, c_ptr, c_size_t
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: COMPILER_VERSION
  USE mo_util_sysinfo, ONLY: util_node_name, util_os_system, util_user_name
  USE mo_cdi, ONLY: cdiLibraryVersion
  USE mo_cf_convention, ONLY: set_cf_global
  USE mo_exception, ONLY: message
  USE mo_mpi, ONLY: my_process_is_global_root

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_icon_version
  PUBLIC :: get_revision
  PUBLIC :: get_remote_url
  PUBLIC :: get_local_branch
  PUBLIC :: show_version

  CHARACTER(*), PARAMETER :: unknown_value = 'unknown'

CONTAINS

  FUNCTION get_icon_version(fallback)
    CHARACTER(:), ALLOCATABLE :: get_icon_version
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      SUBROUTINE pvcs_get_icon_version(ver, ver_len) BIND(c)
        IMPORT c_ptr, c_size_t
        TYPE(c_ptr), INTENT(out) :: ver
        INTEGER(c_size_t), INTENT(out) :: ver_len
      END SUBROUTINE pvcs_get_icon_version
    END INTERFACE
    TYPE(c_ptr) :: c_ver
    INTEGER(c_size_t) :: c_ver_len
    CALL pvcs_get_icon_version(c_ver, c_ver_len)
    IF (PRESENT(fallback)) THEN
      get_icon_version = c_ptr_to_str(c_ver, c_ver_len, fallback)
    ELSE
      get_icon_version = c_ptr_to_str(c_ver, c_ver_len, unknown_value)
    END IF
  END FUNCTION get_icon_version

  FUNCTION get_revision(key, fallback)
    CHARACTER(:), ALLOCATABLE :: get_revision
    CHARACTER(*, c_char), INTENT(in) :: key
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      INTEGER(c_int) FUNCTION pvcs_get_revision(key, key_len, val, val_len) BIND(c)
        IMPORT c_int, c_ptr, c_size_t
        TYPE(c_ptr), VALUE, INTENT(in) :: key
        INTEGER(c_size_t), VALUE, INTENT(in) :: key_len
        TYPE(c_ptr), INTENT(out) :: val
        INTEGER(c_size_t), INTENT(out) :: val_len
      END FUNCTION pvcs_get_revision
    END INTERFACE
    TYPE(c_ptr) :: c_val
    INTEGER(c_size_t) :: c_val_len
    IF (pvcs_get_revision(str_to_c_ptr(key), &
                          INT(LEN(key), kind=c_size_t), &
                          c_val, &
                          c_val_len) == 0_c_int) THEN
      IF (PRESENT(fallback)) THEN
        get_revision = c_ptr_to_str(c_val, c_val_len, fallback)
      ELSE
        get_revision = c_ptr_to_str(c_val, c_val_len, unknown_value)
      END IF
    ELSE IF (PRESENT(fallback)) THEN
      get_revision = fallback
    ELSE
      get_revision = 'external'
    END IF
  END FUNCTION get_revision

  FUNCTION get_remote_url(key, fallback)
    CHARACTER(:), ALLOCATABLE :: get_remote_url
    CHARACTER(*, c_char), INTENT(in) :: key
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      INTEGER(c_int) FUNCTION pvcs_get_remote_url(key, key_len, val, val_len) BIND(c)
        IMPORT c_int, c_ptr, c_size_t
        TYPE(c_ptr), VALUE, INTENT(in) :: key
        INTEGER(c_size_t), VALUE, INTENT(in) :: key_len
        TYPE(c_ptr), INTENT(out) :: val
        INTEGER(c_size_t), INTENT(out) :: val_len
      END FUNCTION pvcs_get_remote_url
    END INTERFACE
    TYPE(c_ptr) :: c_val
    INTEGER(c_size_t) :: c_val_len
    IF (pvcs_get_remote_url(str_to_c_ptr(key), &
                            INT(LEN(key), kind=c_size_t), &
                            c_val, &
                            c_val_len) == 0_c_int) THEN
      IF (PRESENT(fallback)) THEN
        get_remote_url = c_ptr_to_str(c_val, c_val_len, fallback)
      ELSE
        get_remote_url = c_ptr_to_str(c_val, c_val_len, unknown_value)
      END IF
    ELSE IF (PRESENT(fallback)) THEN
      get_remote_url = fallback
    ELSE
      get_remote_url = 'external'
    END IF
  END FUNCTION get_remote_url

  FUNCTION get_local_branch(key, fallback)
    CHARACTER(:), ALLOCATABLE :: get_local_branch
    CHARACTER(*, c_char), INTENT(in) :: key
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      INTEGER(c_int) FUNCTION pvcs_get_local_branch(key, key_len, val, val_len) BIND(c)
        IMPORT c_int, c_ptr, c_size_t
        TYPE(c_ptr), VALUE, INTENT(in) :: key
        INTEGER(c_size_t), VALUE, INTENT(in) :: key_len
        TYPE(c_ptr), INTENT(out) :: val
        INTEGER(c_size_t), INTENT(out) :: val_len
      END FUNCTION pvcs_get_local_branch
    END INTERFACE
    TYPE(c_ptr) :: c_val
    INTEGER(c_size_t) :: c_val_len
    IF (pvcs_get_local_branch(str_to_c_ptr(key), &
                              INT(LEN(key), kind=c_size_t), &
                              c_val, &
                              c_val_len) == 0_c_int) THEN
      IF (PRESENT(fallback)) THEN
        get_local_branch = c_ptr_to_str(c_val, c_val_len, fallback)
      ELSE
        get_local_branch = c_ptr_to_str(c_val, c_val_len, unknown_value)
      END IF
    ELSE IF (PRESENT(fallback)) THEN
      get_local_branch = fallback
    ELSE
      get_local_branch = 'external'
    END IF
  END FUNCTION get_local_branch

  SUBROUTINE show_version()
    CHARACTER(:), ALLOCATABLE :: &
      icon_version, icon_revision, icon_remote_url, icon_local_branch, &
      executable, date_string, time_string, user_name, host_name, os_name, &
      tmp_string1, tmp_string2
    CHARACTER(256) :: buf1, buf2
    INTEGER :: buf_length

    icon_version = get_icon_version()
    icon_revision = get_revision('icon')
    icon_remote_url = get_remote_url('icon')
    icon_local_branch = get_local_branch('icon')

    buf1 = ''
    CALL GET_COMMAND_ARGUMENT(0, buf1, buf_length)
    executable = buf1(1:MIN(buf_length, LEN_TRIM(buf1)))

    buf1 = ''
    buf2 = ''
    CALL DATE_AND_TIME(buf1, buf2)
    IF (LEN_TRIM(buf1) > 0) THEN
      date_string = buf1(1:8)
    ELSE
      date_string = ''
    END IF
    IF (LEN_TRIM(buf2) > 0) THEN
      time_string = buf2(1:6)
    ELSE
      time_string = ''
    END IF

    buf1 = ''
    CALL util_user_name(buf1, buf_length)
    user_name = buf1(1:MIN(buf_length, LEN_TRIM(buf1)))

    buf1 = ''
    CALL util_node_name(buf1, buf_length)
    host_name = buf1(1:MIN(buf_length, LEN_TRIM(buf1)))

    buf1 = ''
    CALL util_os_system(buf1, buf_length)
    os_name = buf1(1:MIN(buf_length, LEN_TRIM(buf1)))

    IF (my_process_is_global_root()) THEN

      tmp_string1 = unknown_value
      IF (LEN(executable) > 0) tmp_string1 = executable
      CALL message('', 'executable: '//tmp_string1)

      tmp_string1 = unknown_value
      IF (LEN(date_string) > 0) tmp_string1 = date_string
      CALL message('', 'date: '//tmp_string1)

      tmp_string1 = unknown_value
      IF (LEN(time_string) > 0) tmp_string1 = time_string
      CALL message('', 'time: '//tmp_string1)

      ! The following is currently redundant because util_user_name returns its
      ! own 'unknown' value when it fails to identify the username:
      tmp_string1 = unknown_value
      IF (LEN(user_name) > 0) tmp_string1 = user_name
      CALL message('', 'user: '//tmp_string1)

      tmp_string1 = 'host: '
      ! The following is currently redundant because util_node_name returns its
      ! own 'unknown' value when it fails to identify the hostname:
      tmp_string2 = unknown_value
      IF (LEN(host_name) > 0) tmp_string2 = host_name
      IF (LEN(os_name) > 0) THEN
        CALL message('', tmp_string1//tmp_string2//' ('//os_name//')')
      ELSE
        ! The following is currently redundant because util_os_system returns
        ! its own 'unknown' value when it fails to identify the OS name:
        CALL message('', tmp_string1//tmp_string2)
      END IF
      CALL message('', 'version: '//icon_version)
      CALL message('', 'revision: '//icon_revision)
      CALL message('', 'repository: '//icon_remote_url)
      CALL message('', 'local branch: '//icon_local_branch)

      CALL message('', 'model components:')
      CALL message('', '  JSBACH: '//get_revision('jsbach'))

      CALL message('', 'application libraries:')
      CALL message('', '  RTE-RRTMGP: '//get_revision('rte-rrtmgp'))

      CALL message('', 'infrastructure and support libraries:')
      CALL message('', '  MATH-INTERPOLATION: '//get_revision('math-interpolation'))
      CALL message('', '  MATH-SUPPORT: '//get_revision('math-support'))
      CALL message('', '  FORTRAN-SUPPORT: '//get_revision('fortran-support'))
      CALL message('', '  MTIME: '//get_revision('mtime'))
      CALL message('', '  CDI:')
      CALL message('', '    version: '//f_ptr_to_str(cdiLibraryVersion()))
      CALL message('', '    revision: '//get_revision('cdi'))

      CALL message('', 'other libraries:')
      tmp_string1 = get_eccodes_version()
      IF (LEN(tmp_string1) > 0) THEN
        CALL message('', '  ECCODES: '//tmp_string1)
      END IF
      CALL message('', '  NetCDF-C: '//get_netcdf_c_version())
      CALL message('', 'compilers:')
      tmp_string1 = get_fortran_compiler_name(fallback='')
      tmp_string2 = get_first_line(COMPILER_VERSION())
      IF (LEN(tmp_string2) > 0) tmp_string2 = ' ('//tmp_string2//')'
      IF (LEN(tmp_string1) > 0) THEN
        CALL message('', '  Fortran: '//tmp_string1//' '//get_fortran_compiler_version()//tmp_string2)
      ELSE
        CALL message('', '  Fortran: '//unknown_value//tmp_string2)
      END IF
      tmp_string1 = get_c_compiler_name(fallback='')
      IF (LEN(tmp_string1) > 0) THEN
        tmp_string2 = get_c_compiler_secondary_name(fallback='')
        IF (LEN(tmp_string2) > 0) THEN
          CALL message( &
            '', &
            '  C: '//tmp_string1//' '//get_c_compiler_version()// &
              & ' ('//tmp_string2//' '//get_c_compiler_secondary_version()//')')
        ELSE
          CALL message('', '  C: '//tmp_string1//' '//get_c_compiler_version())
        END IF
      ELSE
        CALL message('', '  C: '//unknown_value)
      END IF
      ! Probtest fails if it does not find the following lines in the log:
      CALL message('', 'probtest metadata:')
      CALL message('', '  Revision: '//icon_revision)
      CALL message('', '  Branch: '//icon_local_branch)
      CALL message('', '')
    END IF

    IF (LEN(executable) == 0) executable = 'unknown executable'
    IF (LEN(date_string) == 0) date_string = 'unknown date'
    IF (LEN(time_string) == 0) time_string = 'unknown time'
    ! The following is currently redundant because util_user_name returns its
    ! own 'unknown' value when it fails to identify the username:
    IF (LEN(user_name) == 0) user_name = 'unknown user'
    ! The following is currently redundant because util_node_name returns its
    ! own 'unknown' value when it fails to identify the hostname:
    IF (LEN(host_name) == 0) host_name = 'unknown host'
    ! The following is currently redundant because util_os_system returns its
    ! own 'unknown' value when it fails to identify the OS name:
    IF (LEN(os_name) == 0) os_name = 'unknown OS'
    CALL set_cf_global( &
      title='ICON simulation', &
      institution='Max Planck Institute for Meteorology/Deutscher Wetterdienst', &
      source='version: '//icon_version//'; revision: '//icon_revision//'; URL: '//icon_remote_url, &
      history=executable//' at '//date_string//' '//time_string, &
      references='see MPIM/DWD publications', &
      comment=user_name//' on '//host_name//' ('//os_name//')')

  END SUBROUTINE show_version




  FUNCTION get_eccodes_version()
    USE mo_cdi, ONLY: gribapiLibraryVersion
    CHARACTER(:), ALLOCATABLE :: get_eccodes_version
    CHARACTER(32) :: buf = ''
    INTEGER(c_int) :: major, minor, patch
    major = 0; minor = 0; patch = 0
    CALL gribapiLibraryVersion(major, minor, patch)
    IF (major /= 0 .OR. minor /= 0 .OR. patch /= 0) THEN
      WRITE (buf, '(i0,2(a,i0))') major, '.', minor, '.', patch
      get_eccodes_version = buf(1:LEN_TRIM(buf))
    ELSE
      ! We might be simply not using ECCODES:
      get_eccodes_version = ''
    END IF
  END FUNCTION get_eccodes_version

  FUNCTION get_netcdf_c_version()
    USE netcdf, ONLY: nf90_inq_libvers
    CHARACTER(:), ALLOCATABLE :: get_netcdf_c_version
    CHARACTER(80) :: buf = ''
    INTEGER :: delim_idx
    buf = nf90_inq_libvers()
    delim_idx = INDEX(buf, ' of ')
    IF (delim_idx > 1) THEN
      get_netcdf_c_version = buf(1:delim_idx - 1)
    ELSE
      get_netcdf_c_version = unknown_value
    END IF
  END FUNCTION get_netcdf_c_version




! Intel

  FUNCTION get_fortran_compiler_name(fallback)
    CHARACTER(:), ALLOCATABLE :: get_fortran_compiler_name
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    ! Avoid the 'unused dummy variable' compiler warning:
    IF (PRESENT(fallback)) THEN
      get_fortran_compiler_name = "GNU"
    ELSE
      get_fortran_compiler_name = "GNU"
    END IF
  END FUNCTION get_fortran_compiler_name

  FUNCTION get_fortran_compiler_version(fallback)
    CHARACTER(:), ALLOCATABLE :: get_fortran_compiler_version
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    CHARACTER(32) :: buf = ''
    WRITE (buf, '(i0,2(a,i0))') 14, &
      '.', 2, '.', 0
    ! Avoid the 'unused dummy variable' compiler warning:
    IF (PRESENT(fallback)) THEN
      get_fortran_compiler_version = buf(1:LEN_TRIM(buf))
    ELSE
      get_fortran_compiler_version = buf(1:LEN_TRIM(buf))
    END IF
  END FUNCTION get_fortran_compiler_version

  FUNCTION get_c_compiler_name(fallback)
    CHARACTER(:), ALLOCATABLE :: get_c_compiler_name
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      SUBROUTINE pvcs_get_c_compiler_name(ver, ver_len) BIND(c)
        IMPORT c_ptr, c_size_t
        TYPE(c_ptr), INTENT(out) :: ver
        INTEGER(c_size_t), INTENT(out) :: ver_len
      END SUBROUTINE pvcs_get_c_compiler_name
    END INTERFACE
    TYPE(c_ptr) :: c_ver
    INTEGER(c_size_t) :: c_ver_len
    CALL pvcs_get_c_compiler_name(c_ver, c_ver_len)
    IF (PRESENT(fallback)) THEN
      get_c_compiler_name = c_ptr_to_str(c_ver, c_ver_len, fallback)
    ELSE
      get_c_compiler_name = c_ptr_to_str(c_ver, c_ver_len, unknown_value)
    END IF
  END FUNCTION get_c_compiler_name

  FUNCTION get_c_compiler_secondary_name(fallback)
    CHARACTER(:), ALLOCATABLE :: get_c_compiler_secondary_name
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      SUBROUTINE pvcs_get_c_compiler_secondary_name(ver, ver_len) BIND(c)
        IMPORT c_ptr, c_size_t
        TYPE(c_ptr), INTENT(out) :: ver
        INTEGER(c_size_t), INTENT(out) :: ver_len
      END SUBROUTINE pvcs_get_c_compiler_secondary_name
    END INTERFACE
    TYPE(c_ptr) :: c_ver
    INTEGER(c_size_t) :: c_ver_len
    CALL pvcs_get_c_compiler_secondary_name(c_ver, c_ver_len)
    IF (PRESENT(fallback)) THEN
      get_c_compiler_secondary_name = c_ptr_to_str(c_ver, c_ver_len, fallback)
    ELSE
      get_c_compiler_secondary_name = c_ptr_to_str(c_ver, c_ver_len, unknown_value)
    END IF
  END FUNCTION get_c_compiler_secondary_name

  FUNCTION get_c_compiler_version(fallback)
    CHARACTER(:), ALLOCATABLE :: get_c_compiler_version
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      SUBROUTINE pvcs_get_c_compiler_version(ver, ver_len) BIND(c)
        IMPORT c_ptr, c_size_t
        TYPE(c_ptr), INTENT(out) :: ver
        INTEGER(c_size_t), INTENT(out) :: ver_len
      END SUBROUTINE pvcs_get_c_compiler_version
    END INTERFACE
    TYPE(c_ptr) :: c_ver
    INTEGER(c_size_t) :: c_ver_len
    CALL pvcs_get_c_compiler_version(c_ver, c_ver_len)
    IF (PRESENT(fallback)) THEN
      get_c_compiler_version = c_ptr_to_str(c_ver, c_ver_len, fallback)
    ELSE
      get_c_compiler_version = c_ptr_to_str(c_ver, c_ver_len, unknown_value)
    END IF
  END FUNCTION get_c_compiler_version

  FUNCTION get_c_compiler_secondary_version(fallback)
    CHARACTER(:), ALLOCATABLE :: get_c_compiler_secondary_version
    CHARACTER(*), OPTIONAL, INTENT(in) :: fallback
    INTERFACE
      SUBROUTINE pvcs_get_c_compiler_secondary_version(ver, ver_len) BIND(c)
        IMPORT c_ptr, c_size_t
        TYPE(c_ptr), INTENT(out) :: ver
        INTEGER(c_size_t), INTENT(out) :: ver_len
      END SUBROUTINE pvcs_get_c_compiler_secondary_version
    END INTERFACE
    TYPE(c_ptr) :: c_ver
    INTEGER(c_size_t) :: c_ver_len
    CALL pvcs_get_c_compiler_secondary_version(c_ver, c_ver_len)
    IF (PRESENT(fallback)) THEN
      get_c_compiler_secondary_version = c_ptr_to_str(c_ver, c_ver_len, fallback)
    ELSE
      get_c_compiler_secondary_version = c_ptr_to_str(c_ver, c_ver_len, unknown_value)
    END IF
  END FUNCTION get_c_compiler_secondary_version

  FUNCTION get_first_line(str, str_len_lim)
    CHARACTER(:), ALLOCATABLE :: get_first_line
    CHARACTER(*), INTENT(in) :: str
    INTEGER, OPTIONAL, INTENT(in) :: str_len_lim
    CHARACTER(*), PARAMETER :: ellipsis = '...'
    INTEGER :: str_len, new_str_len
    str_len = LEN_TRIM(str)
    new_str_len = INDEX(str(1:str_len), NEW_LINE(str))
    IF (new_str_len > 1) THEN
      get_first_line = str(1:new_str_len - 1)//ellipsis
      new_str_len = new_str_len + 2
    ELSE
      get_first_line = str(1:str_len)
      new_str_len = str_len
    END IF
    IF (PRESENT(str_len_lim)) THEN
      IF (new_str_len > str_len_lim) THEN
        IF (str_len_lim > LEN(ellipsis)*2 + 1) THEN
          get_first_line = get_first_line(1:str_len_lim - LEN(ellipsis))//ellipsis
        ELSE
          get_first_line = get_first_line(1:str_len_lim)
        END IF
      END IF
    END IF
  END FUNCTION

  FUNCTION f_ptr_to_str(ptr)
    CHARACTER(:), ALLOCATABLE :: f_ptr_to_str
    CHARACTER(c_char), POINTER, INTENT(in) :: ptr(:)
    INTEGER(c_size_t) :: ii
    ALLOCATE (CHARACTER(SIZE(ptr)) :: f_ptr_to_str)
    DO ii = 1, SIZE(ptr)
      f_ptr_to_str(ii:ii) = ptr(ii)
    END DO
  END FUNCTION f_ptr_to_str

  FUNCTION c_ptr_to_str(ptr, str_len, fallback)
    CHARACTER(:), ALLOCATABLE :: c_ptr_to_str
    TYPE(c_ptr), INTENT(in) :: ptr
    INTEGER(c_size_t), INTENT(in) :: str_len
    CHARACTER(*), INTENT(in) :: fallback
    CHARACTER(c_char), POINTER :: f_ptr(:)
    IF (C_ASSOCIATED(ptr)) THEN
      CALL C_F_POINTER(ptr, f_ptr, [str_len])
      c_ptr_to_str = f_ptr_to_str(f_ptr)
    ELSE
      c_ptr_to_str = fallback
    END IF
  END FUNCTION c_ptr_to_str

  TYPE(c_ptr) FUNCTION str_to_c_ptr(str)
    CHARACTER(*, c_char), TARGET, INTENT(in) :: str
    str_to_c_ptr = C_LOC(str)
  END FUNCTION str_to_c_ptr

END MODULE mo_util_vcs
