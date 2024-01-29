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

MODULE mo_mtrace
  USE iso_c_binding, ONLY: c_char, c_null_char, c_int
  USE mpi
  USE mo_exception, ONLY: em_warn, message
  PRIVATE
  PUBLIC :: start_memory_tracing
  INTERFACE
    FUNCTION setenv(name, VALUE, overwrite) RESULT(ierror) &
         BIND(c, name='setenv')
      IMPORT :: c_char, c_int
      INTEGER(c_int) :: ierror
      CHARACTER(len=1, kind=c_char), INTENT(in) :: name(*), VALUE(*)
      INTEGER(c_int), VALUE :: overwrite
    END FUNCTION setenv
    SUBROUTINE mtrace() BIND(c, name='mtrace')
    END SUBROUTINE mtrace
  END INTERFACE
  CHARACTER(len=*), PARAMETER :: modname = 'mo_mtrace::'
CONTAINS
  SUBROUTINE start_memory_tracing
    CHARACTER(len=1, kind=c_char), PARAMETER :: &
         mtrace_fn_env(13) = (/ 'M', 'A', 'L', 'L', 'O', 'C', '_', &
         &                      'T', 'R', 'A', 'C', 'E', c_null_char /)
    CHARACTER(len=*), PARAMETER :: fn_pfx = 'ICON_MTRACE_'
    INTEGER, PARAMETER :: pfx_len = LEN(fn_pfx)
    CHARACTER(len=10) :: mt_fn_sfx, sfx_fmt
    CHARACTER(len=1) :: mt_fn(128)
    INTEGER :: i, tlen, ndig, nmax, world_size, world_rank
    INTEGER(c_int) :: se_rc, ow
    LOGICAL :: got_rank
    CHARACTER(len=*), PARAMETER :: routine = modname//'start_memory_tracing'

    got_rank = get_world_rank_slurm(world_size, world_rank)
    IF (.NOT. got_rank) got_rank = get_world_rank_mpi(world_size, world_rank)
    IF (.NOT. got_rank) THEN
      CALL message(routine, 'Cannot determine output file for memory tracing', &
           all_print=.TRUE., level=em_warn)
      RETURN
    END IF
    ndig = 1
    nmax = 10
    DO WHILE (world_size >= nmax)
      ndig = ndig + 1
      nmax = nmax * 10
    END DO
    WRITE (sfx_fmt, '(a,i0,a,i0,a)') '(i', ndig, '.', ndig, ')'
    WRITE (mt_fn_sfx, sfx_fmt) world_rank
    tlen = LEN_TRIM(mt_fn_sfx)
    DO i = 1, pfx_len
      mt_fn(i) = fn_pfx(i:i)
    END DO
    DO i = pfx_len+1, pfx_len+tlen
      mt_fn(i) = mt_fn_sfx(i-pfx_len:i-pfx_len)
    END DO
    mt_fn(pfx_len+tlen+1) = c_null_char
    ow = 1_c_int
    se_rc =  setenv(mtrace_fn_env, mt_fn, ow)
    IF (se_rc /= 0) THEN
      CALL message(routine, 'Failed to set environment variable!', &
           level=em_warn, all_print=.TRUE.)
      RETURN
    END IF
    CALL mtrace
  END SUBROUTINE start_memory_tracing

  FUNCTION get_world_rank_mpi(world_size, world_rank) RESULT(succeeded)
    INTEGER, INTENT(out) :: world_size, world_rank

    LOGICAL :: is_mpi_initialized, succeeded
    INTEGER :: ierror
    CHARACTER(len=*), PARAMETER :: routine = modname//'get_world_rank_mpi'

    succeeded = .FALSE.
    CALL mpi_initialized(is_mpi_initialized, ierror)
    IF (ierror /= mpi_success) THEN
      CALL message(routine, 'Error in mpi_initialized!', level=em_warn, &
        &          all_print=.TRUE.)
      RETURN
    END IF
    IF (.NOT. is_mpi_initialized) RETURN
    CALL mpi_comm_size(mpi_comm_world, world_size, ierror)
    IF (ierror /= mpi_success) THEN
      CALL message(routine, 'Failed to retrieve world size!', level=em_warn, &
        &          all_print=.TRUE.)
      RETURN
    END IF
    CALL mpi_comm_rank(mpi_comm_world, world_rank, ierror)
    IF (ierror /= mpi_success) THEN
      CALL message(routine, 'Failed to retrieve world rank!', level=em_warn, &
        &          all_print=.TRUE.)
      RETURN
    END IF
    succeeded = .TRUE.
  END FUNCTION get_world_rank_mpi

  FUNCTION get_world_rank_slurm(world_size, world_rank) RESULT(succeeded)
    INTEGER, INTENT(out) :: world_size, world_rank
    LOGICAL :: succeeded

    CALL get_int_from_env('SLURM_PROCID', world_rank, succeeded)
    IF (succeeded) THEN
      CALL get_int_from_env('SLURM_NTASKS', world_size, succeeded)
      IF (.NOT. succeeded) &
        &  CALL get_int_from_env('SLURM_NPROCS', world_size, succeeded)
    END IF
  END FUNCTION get_world_rank_slurm

  SUBROUTINE get_int_from_env(ename, ival, succeeded)
    CHARACTER(len=*), INTENT(in) :: ename
    INTEGER, INTENT(out) :: ival
    LOGICAL, INTENT(out) :: succeeded

    INTEGER :: ierror, e_len
    CHARACTER(LEN=18) :: evar

    succeeded = .FALSE.
    CALL get_environment_VARIABLE(ename, evar, e_len, ierror)
    IF (ierror /= 0) RETURN
    READ (evar, '(i18)', iostat=ierror) ival
    succeeded = ierror == 0
  END SUBROUTINE get_int_from_env
END MODULE mo_mtrace
