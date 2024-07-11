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

!> Providing a query tool for the maximum resident size.

MODULE mo_util_rusage

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_long

  IMPLICIT NONE

  PRIVATE

  TYPE, BIND(c) :: timeval
    INTEGER(c_long) :: tv_sec ! seconds
#ifdef __linux__
    INTEGER(c_long)  :: tv_usec ! and microseconds
#else
    INTEGER(c_int)  :: tv_usec ! and microseconds
#endif
  END TYPE timeval

  TYPE, BIND(c) :: rusage
    TYPE(timeval)   :: ru_utime ! user time used
    TYPE(timeval)   :: ru_stime ! system time used
    INTEGER(c_long) :: ru_maxrss ! max resident set size
    INTEGER(c_long) :: ru_ixrss ! integral shared text memory size
    INTEGER(c_long) :: ru_idrss ! integral unshared data size
    INTEGER(c_long) :: ru_isrss ! integral unshared stack size
    INTEGER(c_long) :: ru_minflt ! page reclaims
    INTEGER(c_long) :: ru_majflt ! page faults
    INTEGER(c_long) :: ru_nswap ! swaps
    INTEGER(c_long) :: ru_inblock ! block input operations
    INTEGER(c_long) :: ru_oublock ! block output operations
    INTEGER(c_long) :: ru_msgsnd ! messages sent
    INTEGER(c_long) :: ru_msgrcv ! messages received
    INTEGER(c_long) :: ru_nsignals ! signals received
    INTEGER(c_long) :: ru_nvcsw ! voluntary context switches
    INTEGER(c_long) :: ru_nivcsw ! involuntary context switches
  END TYPE rusage

  INTEGER(c_int), PARAMETER :: RUSAGE_SELF = 0
  INTEGER(c_int), PARAMETER :: RUSAGE_CHILDREN = -1

  INTERFACE
    FUNCTION getrusage(who, r_usage) RESULT(r) BIND(c, name='getrusage')
      IMPORT :: c_int, c_long, rusage
      INTEGER(c_int) :: r
      INTEGER(c_int), VALUE :: who
      TYPE(rusage), INTENT(INOUT) :: r_usage
    END FUNCTION getrusage
  END INTERFACE

  TYPE rss
    TYPE(rusage) :: used_rss
    INTEGER :: idx
  END TYPE rss

  TYPE rss_list
    CHARACTER(len=32) :: name
    INTEGER :: idx
    INTEGER :: used_rss
    TYPE(rss), ALLOCATABLE :: rss_usage(:)
    CHARACTER(len=256) :: filename
    INTEGER :: fileunit
  END TYPE rss_list

  INTEGER, SAVE :: used_rss_lists = 0
  TYPE(rss_list), ALLOCATABLE :: rss_lists(:)

  INTEGER, PARAMETER :: max_lists = 32
  INTEGER, PARAMETER :: max_list_size = 409600
  INTEGER, PARAMETER :: line_length = 1024

  PUBLIC :: add_rss_list
  PUBLIC :: add_rss_usage
  PUBLIC :: print_rss_usage
  PUBLIC :: close_rss_lists

CONTAINS

  FUNCTION add_rss_list(name, tag) RESULT(idx)
    INTEGER :: idx
    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=10), INTENT(INOUT), OPTIONAL :: tag

    TYPE(rss_list), ALLOCATABLE :: tmp_rss_lists(:)
    INTEGER :: iostat

    IF (.NOT. ALLOCATED(rss_lists)) THEN
      ALLOCATE (rss_lists(max_lists))
      used_rss_lists = 0
    END IF

    IF (used_rss_lists == SIZE(rss_lists)) THEN
      ALLOCATE (tmp_rss_lists(2*used_rss_lists))
      tmp_rss_lists(1:SIZE(rss_lists)) = rss_lists(:)
      CALL MOVE_ALLOC(to=rss_lists, from=tmp_rss_lists)
    END IF

    used_rss_lists = used_rss_lists + 1

    idx = used_rss_lists
    rss_lists(idx)%idx = used_rss_lists
    rss_lists(idx)%name = name
    ALLOCATE (rss_lists(idx)%rss_usage(max_list_size))
    rss_lists(idx)%rss_usage(:)%idx = 0
    rss_lists(idx)%used_rss = 0

    IF (PRESENT(tag)) THEN
      rss_lists(idx)%filename = TRIM(name)//'_'//TRIM(tag)//'.log'
      rss_lists(idx)%fileunit = find_next_free_unit(10, 999)
      OPEN (UNIT=rss_lists(idx)%fileunit, FILE=rss_lists(idx)%filename, IOSTAT=iostat, Recl=line_length)
      WRITE (rss_lists(idx)%fileunit, '(1x,a)') 'idx        maxrss      majflt      minflt       nvcsw      nivcsw'
    ELSE
      rss_lists(idx)%filename = ''
      rss_lists(idx)%fileunit = -1
    END IF

  END FUNCTION add_rss_list

  SUBROUTINE add_rss_usage(list_idx)
    INTEGER, INTENT(IN) :: list_idx

    TYPE(rusage) :: ru
    INTEGER :: max_size, idx, ret
    TYPE(rss), ALLOCATABLE :: tmp_rss_usage(:)
    CHARACTER(len=line_length)     :: line

    IF ((list_idx < 1) .OR. (list_idx > used_rss_lists)) THEN
      PRINT *, 'ERROR: list does not exist ...'
    END IF

    max_size = SIZE(rss_lists(list_idx)%rss_usage)
    idx = rss_lists(list_idx)%used_rss

    IF (idx == max_list_size) THEN
      ALLOCATE (tmp_rss_usage(2*max_size))
      tmp_rss_usage(1:max_size) = rss_lists(list_idx)%rss_usage(:)
      CALL MOVE_ALLOC(to=rss_lists(list_idx)%rss_usage, from=tmp_rss_usage)
    END IF

    rss_lists(list_idx)%used_rss = rss_lists(list_idx)%used_rss + 1

    idx = rss_lists(list_idx)%used_rss
    rss_lists(list_idx)%rss_usage%idx = idx
    ret = getrusage(RUSAGE_SELF, ru)
    rss_lists(list_idx)%rss_usage(idx)%used_rss = ru

    IF (-1 .NE. rss_lists(list_idx)%fileunit) THEN
      line = ''
      WRITE (line, '(1x,i5,5i12)') &
          &          idx, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_maxrss, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_majflt, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_minflt, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_nvcsw, &
          &          rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_nivcsw
      !write(line,'(i6,a,i12)')idx,' ',rss_lists(list_idx)%rss_usage(idx)%used_rss%ru_maxrss
      WRITE (rss_lists(list_idx)%fileunit, '(a)') TRIM(line)
      FLUSH (rss_lists(list_idx)%fileunit)
    END IF

  END SUBROUTINE add_rss_usage

  SUBROUTINE print_rss_usage()
    INTEGER                        :: il, idx

    DO il = 1, used_rss_lists
      PRINT *, 'List: ', TRIM(rss_lists(il)%name), ' index: ', rss_lists(il)%idx
      PRINT '(1x,a)', '            maxrss      majflt      minflt       nvcsw      nivcsw'
      DO idx = 1, rss_lists(il)%used_rss
        PRINT '(1x,i5,5i12)', &
          &          idx, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_maxrss, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_majflt, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_minflt, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_nvcsw, &
          &          rss_lists(il)%rss_usage(idx)%used_rss%ru_nivcsw
      END DO
    END DO

    PRINT *, ''
#ifdef __linux__
    PRINT *, ' maxrss - Maximum resident set size (kbytes)'
#else
    PRINT *, ' maxrss - Maximum resident set size (bytes)'
#endif
    PRINT *, ' majflt - Major (requiring I/O) page faults'
    PRINT *, ' minflt - Minor (reclaiming a frame) page faults'
    PRINT *, ' nvcsw  - Voluntary context switches'
    PRINT *, ' nivcsw - Involuntary context switches'

  END SUBROUTINE print_rss_usage

  SUBROUTINE close_rss_lists(verbose)
    LOGICAL, OPTIONAL :: verbose

    LOGICAL :: my_verbose
    INTEGER :: il

    my_verbose = .FALSE.
    IF (PRESENT(verbose)) my_verbose = verbose

    DO il = 1, used_rss_lists
      IF (-1 .NE. rss_lists(il)%fileunit) THEN
        IF (my_verbose) WRITE (0, *) 'CLOSE: ', TRIM(rss_lists(il)%filename), ' index: ', rss_lists(il)%idx
        CLOSE (unit=rss_lists(il)%fileunit)
      END IF
    END DO
  END SUBROUTINE close_rss_lists

  FUNCTION find_next_free_unit(istart, istop) RESULT(iunit)
    INTEGER :: iunit
    INTEGER, INTENT(IN) :: istart, istop
    !
    INTEGER :: kstart, kstop
    LOGICAL :: lfound, lopened
    INTEGER :: i
    !
    lfound = .FALSE.
    !
    kstart = istart
    kstop = istop
    IF (kstart < 10) kstart = 10
    IF (kstop <= kstart) kstop = kstart + 10
    !
    DO i = kstart, kstop
      INQUIRE (unit=i, opened=lopened)
      IF (.NOT. lopened) THEN
        iunit = i
        lfound = .TRUE.
        EXIT
      END IF
    END DO
    !
    IF (.NOT. lfound) THEN
      iunit = -1
    END IF
    !
  END FUNCTION find_next_free_unit
END MODULE mo_util_rusage
