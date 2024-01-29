! Routines for handling of regular output steps and ready file events
! on multiple I/O PEs.
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
!
!
! Note that this module contains the abstract implementation of the
! output event handling, while the actual control of the trigger
! mechanism is located in the module "mo_output_event_control".
!
! Detailed description:
!
! 1. The writing of regular output is triggered at so-called "output
!    event steps". Additional tasks for such an output event step may
!    be the opening of a new output file or the closing of an output
!    file. There may be several output events, each with a different
!    number of output event steps.
!
! 2. Output event steps happen at regular intervals. These are given
!    by an interval size and the time stamps for begin and end.  One
!    may think of each output event as being handled by its own
!    PE. The completion of an output event step is communicated via
!    non-blocking MPI messages to the root PE, which keeps track of
!    the overall event status.
!
! 3. Different output events/PEs may be joined together to form a
!    single new event. Then not every PE necessarily participates in
!    every event step. The purpose of joining/grouping events is that
!    certain actions can be triggered for the group, e.g. writing
!    "ready files".
!
!    A "ready file" is a technique for handling dependencies between
!    the NWP processes at DWD: When a program - parallel or
!    sequential, shell script or binary - produces some output which
!    is necessary for other running applications, then the completion
!    of the write process signals this by creating a small file
!    (size: a few bytes). Only when this file exists, the second
!    program starts reading its input data. Implicity, this assumes
!    that a file system creates (and closes) files in the same order
!    as they are written by the program.
!
! @note This event handling is *static*, i.e. all event occurrences
!       are pre-defined during the initialization. This is necessary
!       for generating a complete list of "ready" files and output
!       files at the beginning of the simulation.
!
! MPI ROLES
! ---------
!
! The output event handling requires communication between the I/O
! PEs and a root I/O MPI rank "ROOT_OUTEVENT" - usually chosen as
! rank 0 from the I/O PEs.
!
! Root PE: Parallel communication is necessary
!
! - During the setup phase:
!
!   The root I/O MPI rank asks all participating I/O PEs for their
!   output event info. This information about these events is
!   forwarded (broadcasted) to the worker PEs. Afterwards, each of
!   the I/O and worker PEs generates a unified output event,
!   indicating which PE performs a write process at which step. Thus,
!   every PE with the exception of the asynchronous restart PEs has
!   the same view on the sequence of output events.
!
!   see FUNCTION union_of_all_events
!
! - During the loop of output steps:
!
!   The root I/O rank places MPI_IRECV calls, thereby asking all
!   participating PEs for acknowledging the completion of the current
!   output step.
!
!   see  SUBROUTINE trigger_output_step_irecv
!        FUNCTION is_output_step_complete
!
! All I/O PEs: Parallel communication is performed
!
! - During setup:
!
!   Output events are generated and the necessary meta-data for
!   replication on the root I/O PE is send via MPI to the root
!   rank. An event name which is different from the string
!   DEFAULT_EVENT_NAME ( = "default") leads to the creation of ready
!   files during the output loop.
!
!   see FUNCTION new_parallel_output_event
!
! - During the loop of output steps:
!
!   The I/O PEs acknowledge the completion of their respective output
!   step by sending non-blocking MPI messages.  The non-blocking
!   messages have unique MPI tags that are computed as i_tag =
!   SENDRECV_TAG_OUTEVENT + this_pe + (local_event_no-1)*nranks
!
!   see SUBROUTINE pass_output_step

MODULE mo_output_event_handler

#ifndef NOMPI
  USE mpi
#endif

  USE mo_kind,                   ONLY: i8, wp
  USE mo_impl_constants,         ONLY: SUCCESS, MAX_TIME_INTERVALS, &
    &                                  pio_type_cdipio
  USE mo_exception,              ONLY: finish, message_text
  USE mo_io_units,               ONLY: FILENAME_MAX, find_next_free_unit
  USE mo_util_string,            ONLY: int2string
  USE mo_mpi,                    ONLY: p_int, p_int_i4, p_int_i8, p_real,                   &
    &                                  get_my_global_mpi_id,                                &
    &                                  p_pe,                                                &
    &                                  my_process_is_mpi_workroot,                          &
    &                                  process_mpi_all_comm,                                &
    &                                  p_comm_rank, p_comm_size,                            &
    &                                  p_work_pe0, p_io_pe0,                                &
    &                                  p_send, p_recv,                                      &
    &                                  p_gather, p_gatherv, mpi_comm_null
  USE mtime,                     ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN,         &
    &                                  datetime, timedelta,  newTimedelta,                  &
    &                                  deallocateDatetime, datetimeToString,                &
    &                                  newDatetime, OPERATOR(>=), OPERATOR(==),&
    &                                  OPERATOR(<), OPERATOR(<=), &
    &                                  OPERATOR(>), OPERATOR(+), OPERATOR(/=),              &
    &                                  deallocateTimedelta, JulianDay,                      &
    &                                  getJulianDayFromDatetime, getDatetimeFromJulianDay
  USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data,                  &
    &                                  t_event_step, t_output_event, t_par_output_event,    &
    &                                  MAX_FILENAME_STR_LEN, MAX_EVENT_NAME_STR_LEN,        &
    &                                  DEFAULT_EVENT_NAME
  USE mo_name_list_output_types, ONLY: t_fname_metadata, t_event_data_local
  USE mo_util_table,             ONLY: initialize_table, finalize_table, add_table_column,  &
    &                                  set_table_entry, print_table, t_table
  USE mo_name_list_output_config,ONLY: use_async_name_list_io
  USE mo_parallel_config,        ONLY: pio_type
  USE iso_c_binding,             ONLY: c_size_t, c_loc
  USE mo_util_libc,              ONLY: memcmp_f, memset_f

  IMPLICIT NONE

  ! public subroutines + functions:
  PRIVATE

  ! initialization + destruction
  PUBLIC :: new_parallel_output_event
  PUBLIC :: union_of_all_events
  PUBLIC :: deallocate_output_event
  PUBLIC :: blocking_wait_for_irecvs
  ! inquiry functions
  PUBLIC :: get_current_date
  PUBLIC :: get_current_step
  PUBLIC :: get_current_filename
  PUBLIC :: get_current_jfile
  PUBLIC :: check_write_readyfile
  PUBLIC :: check_open_file
  PUBLIC :: check_close_file
  PUBLIC :: is_event_root_pe
  PUBLIC :: is_output_step
  PUBLIC :: is_output_step_complete
  PUBLIC :: is_output_event_finished
  ! handshake routines (called after steps)
  PUBLIC :: pass_output_step
#ifndef NOMPI
  PUBLIC :: trigger_output_step_irecv
#endif
  ! auxiliary functions
  PUBLIC :: print_output_event
  PUBLIC :: set_event_to_simstep

  INTERFACE deallocate_output_event
    MODULE PROCEDURE deallocate_output_event
    MODULE PROCEDURE deallocate_par_output_event
  END INTERFACE

  INTERFACE is_output_event_finished
    MODULE PROCEDURE is_output_event_finished
    MODULE PROCEDURE is_par_output_event_finished
  END INTERFACE

  INTERFACE is_output_step
    MODULE PROCEDURE is_output_step
    MODULE PROCEDURE is_par_output_step
  END INTERFACE

  INTERFACE check_close_file
    MODULE PROCEDURE check_close_file
    MODULE PROCEDURE check_close_file_step
  END INTERFACE check_close_file

  INTERFACE get_current_date
    MODULE PROCEDURE get_current_date
    MODULE PROCEDURE get_current_date_par
  END INTERFACE

  INTERFACE get_current_step
    MODULE PROCEDURE get_current_step
    MODULE PROCEDURE get_current_step_par
  END INTERFACE

  INTERFACE print_output_event
    MODULE PROCEDURE print_output_event
    MODULE PROCEDURE print_par_output_event
  END INTERFACE

  INTERFACE set_event_to_simstep
    MODULE PROCEDURE set_event_to_simstep
    MODULE PROCEDURE set_event_to_simstep_par
  END INTERFACE set_event_to_simstep

  INTERFACE p_gatherv
    MODULE PROCEDURE p_gatherv_event_data_1d1d
  END INTERFACE p_gatherv

#ifndef NOMPI
  INTERFACE p_send
    MODULE PROCEDURE p_send_event_data_1d
  END INTERFACE p_send

  INTERFACE p_recv
    MODULE PROCEDURE p_recv_event_data_1d
  END INTERFACE p_recv
#endif

  !---------------------------------------------------------------
  ! constants

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output_event_handler'

  !> maximum buffer size for sending event meta-data (MPI_PACK)
  INTEGER, PARAMETER :: MAX_BUF_SIZE =    (    MAX_EVENT_NAME_STR_LEN    &
    &                                      + 9*(MAX_DATETIME_STR_LEN+1)  &
    &                                      + MAX_TIMEDELTA_STR_LEN       &
    &                                      + 3*FILENAME_MAX              &
    &                                      + 1024 )

  !> MPI message tag for setup of output events
  INTEGER, PARAMETER :: SENDRECV_TAG_SETUP    = 1002

  !> MPI message tag for output event handshake
  INTEGER, PARAMETER :: SENDRECV_TAG_OUTEVENT = 1001

  !> Internal switch for debugging output
  LOGICAL, PARAMETER :: ldebug                = .FALSE.

  !> Max. no. of steps printed out to stderr
  INTEGER, PARAMETER :: MAX_PRINTOUT          = 2000


#ifndef NOMPI
  INTEGER :: event_data_dt = mpi_datatype_null
#endif
CONTAINS

  !---------------------------------------------------------------
  ! ROUTINES FOR DEALLOCATION
  !---------------------------------------------------------------

  !> Deallocate t_event_step data structure.
  !
  SUBROUTINE deallocate_event_step(event_step)
    TYPE(t_event_step), INTENT(INOUT) :: event_step
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_event_step"
    INTEGER :: ierrstat
    DEALLOCATE(event_step%event_step_data, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE deallocate_event_step


  !> Deallocate t_output_event data structure.
  !
  SUBROUTINE deallocate_output_event(event)
    TYPE(t_output_event), ALLOCATABLE, INTENT(inout) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_output_event"
    INTEGER :: ierrstat,i
    IF (ALLOCATED(event)) THEN
      DO i=1,event%n_event_steps
        CALL deallocate_event_step(event%event_step(i))
      END DO
      event%n_event_steps = 0
      DEALLOCATE(event%event_step, STAT=ierrstat)
      IF (ierrstat == SUCCESS) DEALLOCATE(event, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
  END SUBROUTINE deallocate_output_event


  !> Deallocate t_par_output_event data structure.
  !
  !  @note This subroutine recursively deallocates the complete
  !        singly-linked list rooted at @p event.
  !
  SUBROUTINE deallocate_par_output_event(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_par_output_event"
    INTEGER :: ierrstat

    TYPE(t_par_output_event), POINTER :: ev, next

    ierrstat = success
    ev => event
    DO WHILE (ASSOCIATED(ev))
      next => ev%next
      CALL deallocate_output_event(ev%output_event)
      IF (ALLOCATED(ev%irecv_req)) THEN
        DEALLOCATE(ev%output_event, ev%irecv_req, ev%irecv_buf, STAT=ierrstat)
      END IF
      IF (ierrstat == success) DEALLOCATE(ev, STAT=ierrstat)
      IF (ierrstat /= success) CALL finish (routine, 'DEALLOCATE failed.')
      ev => next
    END DO
    NULLIFY(event)
  END SUBROUTINE deallocate_par_output_event


  !---------------------------------------------------------------
  ! ROUTINES FOR PRINTING EVENT DESCRIPTIONS
  !---------------------------------------------------------------

  !> Screen print-out of an output event.
  !
  !  @note If an (optional) destination file is provided, we
  !        implicitly assume that this file has already been opened by
  !        the called and will be closed by the caller.
  !
  SUBROUTINE print_output_event(event, opt_dstfile)
    TYPE(t_output_event),          INTENT(IN) :: event           !< output event data structure
    INTEGER, OPTIONAL,             INTENT(IN) :: opt_dstfile     !< optional destination ASCII file unit
    ! local variables
    INTEGER                     :: i, irow, dst, n
    TYPE(t_table)               :: table

    dst = 0
    IF (PRESENT(opt_dstfile)) dst = opt_dstfile

    IF (ldebug) THEN
      WRITE (0,*) "PE ", get_my_global_mpi_id(), ": event%n_event_steps = ", event%n_event_steps
    END IF

    IF (check_write_readyfile(event)) THEN
      WRITE (dst,'(3a)') 'output "', TRIM(event%event_data%name), '", writes ready files:'
    ELSE
      WRITE (dst,'(3a)') 'output "', TRIM(event%event_data%name), '", does not write ready files:'
    END IF

    ! do not print, if the table is excessively long
    IF (event%n_event_steps > MAX_PRINTOUT) THEN
      WRITE (dst,*) "detailed print-out of output steps has been omitted, cf. 'output_schedule.txt'."
      RETURN
    END IF

    ! table-based output
    CALL initialize_table(table)
    CALL add_table_column(table, "model step")
    CALL add_table_column(table, "model date")
    CALL add_table_column(table, "filename")
    CALL add_table_column(table, "I/O PE")
    CALL add_table_column(table, "output date")

    CALL add_table_column(table, "#")

    CALL add_table_column(table, "open")
    CALL add_table_column(table, "close")
    irow = 0
    n = event%n_event_steps
    DO i = 1, n
      CALL tabulate_event_step(table, irow, dst, event%event_step(i), i, n)
    END DO
    CALL print_table(table, opt_delimiter='   ', opt_dstfile=dst)
    CALL finalize_table(table)
  END SUBROUTINE print_output_event

  SUBROUTINE tabulate_event_step(table, irow, dst, event_step, &
    &                            step_idx, num_steps)
    TYPE(t_table), INTENT(inout) :: table
    INTEGER, INTENT(inout) :: irow
    INTEGER, INTENT(in) :: dst, step_idx, num_steps
    TYPE(t_event_step), INTENT(in) :: event_step
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: event_datetime_string
    CHARACTER(LEN=1+9+9) :: cfilepart

    INTEGER :: j, tlen
    LOGICAL :: lclose, lopen
    DO j=1,event_step%n_pes
      irow = irow + 1
      IF (j==1) THEN
        CALL set_table_entry(table,irow,"model step", int2string(event_step%i_sim_step))
        CALL datetimetostring(event_step%exact_date, event_datetime_string)
        CALL set_table_entry(table,irow,"model date", TRIM(event_datetime_string))
      ELSE
        CALL set_table_entry(table,irow,"model step", " ")
        CALL set_table_entry(table,irow,"model date", " ")
      END IF
      tlen = LEN_TRIM(event_step%event_step_data(j)%filename_string)

      CALL set_table_entry(table,irow,"filename",    event_step%event_step_data(j)%filename_string(1:tlen))

      CALL set_table_entry(table,irow,"I/O PE",      int2string(event_step%event_step_data(j)%i_pe))

      CALL datetimetostring(event_step%event_step_data(j)%datetime, &
        &                   event_datetime_string)
      CALL set_table_entry(table,irow,"output date", TRIM(event_datetime_string))

      WRITE (cfilepart, '(i0,a,i0)') event_step%event_step_data(j)%jfile, ".", &
           event_step%event_step_data(j)%jpart
      CALL set_table_entry(table,irow,"#", cfilepart)

      ! append "+ open" or "+ close" according to event step data:
      lopen = event_step%event_step_data(j)%l_open_file
      CALL set_table_entry(table,irow,"open", MERGE("x", " ", lopen))
      lclose = step_idx > 1 &
        .AND. (     step_idx == num_steps &
        &      .OR. event_step%event_step_data(j)%l_open_file)
      CALL set_table_entry(table,irow,"close", MERGE("x", " ", lclose))
    END DO
  END SUBROUTINE tabulate_event_step


  !> Screen print-out of a parallel output event.
  !
  !  @note This subroutine recursively prints the complete
  !        singly-linked list rooted at @p event.
  !
  SUBROUTINE print_par_output_event(event, opt_filename, opt_dstfile)
    TYPE(t_par_output_event), TARGET,  INTENT(IN) :: event
    CHARACTER(LEN=*), OPTIONAL,        INTENT(IN) :: opt_filename    !< name of ASCII file (optional)
    INTEGER,          OPTIONAL,        INTENT(IN) :: opt_dstfile     !< optional destination ASCII file unit
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::print_par_output_event"
    INTEGER :: dst, ierrstat
    TYPE(t_par_output_event), POINTER :: ev

    ! consistency check:
    IF (PRESENT(opt_dstfile) .AND. PRESENT(opt_filename)) THEN
      CALL finish(routine, "Routine was called with both a file unit and a filename to open!")
    END IF

    ! do not print, if the table is excessively long
    dst = 0
    IF (event%output_event%n_event_steps > MAX_PRINTOUT) THEN
      WRITE (dst,*) "detailed print-out of output steps has been omitted, cf. 'output_schedule.txt'."
      RETURN
    END IF

    ! open ASCII output file (if necessary):
    IF (PRESENT(opt_dstfile)) THEN
      dst = opt_dstfile
    ELSE IF (PRESENT(opt_filename)) THEN
      dst = find_next_free_unit(10,100)
      OPEN (dst, file=TRIM(opt_filename), status='REPLACE', iostat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'OPEN failed.')
    END IF

    ev => event
    DO WHILE (ASSOCIATED(ev))
      WRITE (dst,'(a)') "" ! newline
      CALL print_output_event(ev%output_event, opt_dstfile=dst)
      ev => ev%next
    END DO
    WRITE (dst,'(a)') "" ! newline

    ! close ASCII output file (if necessary):
    IF (PRESENT(opt_filename)) THEN
      CLOSE (dst, iostat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'CLOSE failed.')
    END IF
  END SUBROUTINE print_par_output_event


  !---------------------------------------------------------------
  ! ROUTINES FOR CREATING NEW OUTPUT EVENTS
  !---------------------------------------------------------------

  !> Create a simple output event, happening at regular intervals.
  !
  !  These are given by an interval size @p intvl_str and the time
  !  stamps for begin and end, @p begin_str and @p end_str.
  !
  !  The begin and end time stamps may contain "modifier symbols",
  !  e.g. ">2014-06-01T00:00:00.000", where ">" means that output
  !  should not include the start date.
  !
  !  Note that this subroutine generates the output event steps but
  !  does not map these time stamps onto the corresponding simulation
  !  steps. This simulation-specific task is performed by a different
  !  subroutine, which is given as a function parameter
  !  "fct_time2simstep".
  !
  SUBROUTINE init_output_event(p_event, evd, i_pe, fct_time2simstep,    &
       &                       fct_generate_filenames)
    TYPE(t_output_event), ALLOCATABLE, INTENT(OUT) :: p_event
    TYPE(t_event_data_local),          INTENT(IN)  :: evd
    !> rank of participating PE
    INTEGER,                           INTENT(IN)  :: i_pe

    INTERFACE
      !> As an argument of this function, the user must provide a
      !!  conversion "time stamp -> simulation step"
      SUBROUTINE fct_time2simstep(num_dates, dates, sim_step_info, &
        &                         result_steps, result_exactdate)
        IMPORT :: t_sim_step_info, datetime
        !> no. of dates to convert
        INTEGER,               INTENT(IN)    :: num_dates
        !> array of mtime datetime objects
        TYPE(datetime),        INTENT(IN)    :: dates(:)
        TYPE(t_sim_step_info), INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,               INTENT(OUT)   :: result_steps(:)      !< resulting step indices
        TYPE(datetime),        INTENT(OUT)   :: result_exactdate(:)  !< resulting (exact) time step
      END SUBROUTINE fct_time2simstep

      !> As an argument of this function, the user must provide a
      !! function for generating output file names
      SUBROUTINE fct_generate_filenames(num_dates, dates, sim_steps, &
        &                             sim_step_info, fname_metadata, &
        &                             skipped_dates, result_fnames)
        IMPORT :: t_sim_step_info, t_event_step_data, t_fname_metadata, &
             datetime

        !> no. of dates to convert
        INTEGER,                 INTENT(IN)    :: num_dates
        !> array of mtime datetime objects
        TYPE(datetime), TARGET,  INTENT(IN)    :: dates(:)
        INTEGER,                 INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),   INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata),  INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        INTEGER,                 INTENT(IN)    :: skipped_dates
        TYPE(t_event_step_data), INTENT(out)   :: result_fnames(SIZE(dates))
      END SUBROUTINE fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_output_event"
    !> Max. no. of event steps (used for local array sizes)
    INTEGER, PARAMETER :: INITIAL_NEVENT_STEPS = 256 ! tested to be a good compromise
    TYPE(datetime) :: mtime_date
    TYPE(datetime),  POINTER :: mtime_begin, mtime_end
    TYPE(datetime) :: sim_end, run_start, mtime_restart, &
      &               mtime_dom_start, mtime_dom_end
    TYPE(datetime), TARGET, ALLOCATABLE :: mtime_dates(:), tmp_dates(:)
    TYPE(timedelta), POINTER :: delta
    INTEGER                  :: ierrstat, i, j, n_event_steps, &
      &                         nintvls, iintvl, skipped_dates, &
      &                         old_size, new_size
    LOGICAL                  :: l_active, l_append_step
    CHARACTER(len=MAX_DATETIME_STR_LEN) :: dtime_string
    INTEGER,                             ALLOCATABLE :: mtime_sim_steps(:)
    TYPE(datetime), ALLOCATABLE :: mtime_exactdate(:)
    TYPE(t_event_step_data),             ALLOCATABLE :: filename_metadata(:)
    TYPE(t_event_step_data),             ALLOCATABLE :: step_data(:)
    CHARACTER(len=MAX_DATETIME_STR_LEN+1)            :: dt_string
    CHARACTER(len=MAX_DATETIME_STR_LEN)              :: begin_str2(MAX_TIME_INTERVALS)
    CHARACTER(len=MAX_DATETIME_STR_LEN)              :: end_str2(MAX_TIME_INTERVALS)
    LOGICAL                                          :: incl_begin(MAX_TIME_INTERVALS)
    LOGICAL                                          :: incl_end(MAX_TIME_INTERVALS)
    CHARACTER                                        :: char
    TYPE(julianday), POINTER :: tmp_jd

    TYPE(julianday), ALLOCATABLE, TARGET :: mtime_date_container_a(:), &
      &                                     mtime_date_container_b(:)
    TYPE(julianday), ALLOCATABLE :: tmp(:)
    TYPE(julianday), ALLOCATABLE :: mtime_date_uniq(:)
    
    TYPE(julianday), POINTER :: mtime_date_container(:)
    
    INTEGER, ALLOCATABLE :: indices_to_use(:)
    INTEGER :: remaining_intvls, iselected_intvl
    
    INTEGER :: n_event_steps_a, n_event_steps_b, remaining_event_steps
    
    ! allocate event data structure
    ALLOCATE(p_event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    n_event_steps                    = 0
    p_event%i_event_step             = 1
    p_event%event_data%name          = evd%name
    p_event%event_data%sim_start     = evd%sim_step_info%sim_start

    ! count the number of different time intervals for this event (usually 1)
    DO nintvls = 0, max_time_intervals-1
      IF (LEN_TRIM(evd%begin_str(nintvls+1)) == 0) EXIT
    END DO

    ! status output
    IF (ldebug) THEN
      WRITE (0,*) "PE ",get_my_global_mpi_id(), ":"
      WRITE (0,'(7a)') 'Defining output event "', &
           TRIM(evd%begin_str(1)), '", "', &
           TRIM(evd%end_str(1)), '", "', &
           TRIM(evd%intvl_str(1)), '"'
      DO i=2,nintvls
        WRITE (0,'(7a)') ' +  "', &
             TRIM(evd%begin_str(i)), '", "', &
             TRIM(evd%end_str(i)), '", "',   &
             TRIM(evd%intvl_str(i)), '"'
      END DO
      CALL datetimeToString(evd%sim_step_info%sim_start, begin_str2(1))
      CALL datetimeToString(evd%sim_step_info%sim_end, begin_str2(2))
      CALL datetimeToString(evd%sim_step_info%restart_time, begin_str2(3))

      WRITE (0,*) 'Simulation bounds:    "'//TRIM(dtime_string)//'", "'//TRIM(begin_str2(2))//'"'
      WRITE (0,*) 'restart bound: "'//TRIM(begin_str2(3))
    END IF

    ! set some dates used later:
    sim_end     = evd%sim_step_info%sim_end
    run_start   = evd%sim_step_info%run_start

    ! Domains (and their output) can be activated and deactivated
    ! during the simulation. This is determined by the parameters
    ! "dom_start_time" and "dom_end_time". Therefore, we must create
    ! a corresponding event.
    mtime_dom_start = evd%sim_step_info%dom_start_time
    mtime_dom_end   = evd%sim_step_info%dom_end_time ! this carries the domain-specific namelist value "end_time"

    ! To avoid further case discriminations, sim_end is set to the minimum of the simulation end time
    ! and the time at which a nested domain is turned off
    IF (sim_end > mtime_dom_end) THEN
      sim_end = evd%sim_step_info%dom_end_time
    ENDIF

    ! Compute the end time wrt. "dt_restart": It might be that the
    ! simulation end is limited by this parameter
    mtime_restart   = evd%sim_step_info%restart_time

    ! loop over the event occurrences
    ALLOCATE(mtime_date_container_a(256), &
      &      mtime_date_container_b(256), stat=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! The begin and end time stamps may contain "modifier symbols",
    ! e.g. ">2014-06-01T00:00:00.000", where ">" means that output
    ! should not include the start date.
    !
    ! At this point, we strip these modifiers from the begin and end
    ! time stamps and store this information in logical flags.
    incl_begin(:) = .TRUE.
    incl_end(:)   = .TRUE.
    DO iintvl=1,nintvls
      ! begin time stamp
      dt_string = evd%begin_str(iintvl)
      char      = dt_string(1:1)
      SELECT CASE(char)
      CASE ('>')
        incl_begin(iintvl) = .FALSE.
        dt_string = dt_string(2:)
      END SELECT
      begin_str2(iintvl) = dt_string

      ! end time stamp
      dt_string = evd%end_str(iintvl)
      char      = dt_string(1:1)
      SELECT CASE(char)
      CASE ('<')
        incl_end(iintvl) = .FALSE.
        dt_string = dt_string(2:)
      END SELECT
      end_str2(iintvl) = dt_string
    END DO

    allocate(indices_to_use(nintvls))
    
    CALL remove_duplicate_intervals(begin_str2, end_str2, evd%intvl_str, nintvls, &
         &                          indices_to_use, remaining_intvls)
    
    ! there may be multiple starts/ends/intervals (usually only one):

    mtime_date_container => mtime_date_container_a
    n_event_steps = 0
    skipped_dates = 0

    INTERVAL_LOOP: DO iselected_intvl = 1, remaining_intvls    
      iintvl = indices_to_use(iselected_intvl)
      
      mtime_begin => newDatetime(TRIM(begin_str2(iintvl)))
      mtime_end   => newDatetime(TRIM(end_str2(iintvl)))
      IF (.NOT. ASSOCIATED(mtime_begin) .OR. .NOT. ASSOCIATED(mtime_end))  THEN
        CALL finish(routine, "date-time conversion error: "&
             //TRIM(begin_str2(iintvl))//" "//TRIM(end_str2(iintvl)))
      END IF

      mtime_date  = mtime_begin
      delta       => newTimedelta(TRIM(evd%intvl_str(iintvl))) ! create a time delta
      IF (mtime_end >= mtime_begin) THEN
        EVENT_LOOP: DO
          IF ((mtime_date    >= run_start)                         .AND. &
            & (sim_end       >= mtime_date)                        .AND. &
            & (mtime_restart >= mtime_date)                        .AND. &
            & (incl_begin(iintvl) .OR. (mtime_date > mtime_begin)) .AND. &
            & (incl_end(iintvl)   .OR. (mtime_end  > mtime_date)) )  THEN

            IF  (mtime_date >= mtime_dom_start) THEN

              n_event_steps = n_event_steps + 1

              old_size = SIZE(mtime_date_container)
              IF (n_event_steps > old_size) THEN
                ALLOCATE(tmp(2*old_size), stat=ierrstat)
                IF (ierrstat /= 0) STOP 'allocate failed'
                tmp(1:old_size) = mtime_date_container(:)
                IF (ASSOCIATED(mtime_date_container, mtime_date_container_a)) THEN
                  CALL MOVE_ALLOC(tmp, mtime_date_container_a)
                  mtime_date_container => mtime_date_container_a
                ELSE
                  CALL MOVE_ALLOC(tmp, mtime_date_container_b)
                  mtime_date_container => mtime_date_container_b
                ENDIF
              ENDIF

              tmp_jd => mtime_date_container(n_event_steps)
              CALL getJulianDayFromDatetime(mtime_date, tmp_jd)

            ELSE
              ! we skip an output date when the domain is not yet
              ! "active" - this leads to a file count offset.
              skipped_dates = skipped_dates+1
            END IF
          END IF

          mtime_date = mtime_date + delta

          l_active = mtime_date <= mtime_end .AND.   &
            &        mtime_date <= sim_end   .AND.   &
            &        mtime_date <= mtime_restart
          IF (.NOT. l_active) EXIT EVENT_LOOP
        END DO EVENT_LOOP
      END IF

      CALL deallocateTimedelta(delta)
      IF (iselected_intvl /= remaining_intvls) CALL deallocateDatetime(mtime_end)
      CALL deallocateDatetime(mtime_begin)

      IF (iintvl == 1) THEN
        n_event_steps_a = n_event_steps
        mtime_date_container => mtime_date_container_b
        n_event_steps = 0
        CYCLE interval_loop
      ELSE
        n_event_steps_b = n_event_steps
        n_event_steps = 0      
        
        IF (n_event_steps_b > 0) THEN
          CALL merge2SortedAndRemoveDuplicates(mtime_date_container_a, n_event_steps_a, &
            &                               mtime_date_container_b, n_event_steps_b, &
            &                               mtime_date_uniq, remaining_event_steps)

          IF (remaining_event_steps > SIZE(mtime_date_container_a)) THEN
            ALLOCATE(tmp(remaining_event_steps), stat=ierrstat)
            IF (ierrstat /= 0) STOP 'allocate failed'
            tmp(1:remaining_event_steps) = mtime_date_uniq(1:remaining_event_steps)
            ! FIXME: this is probably a bug and the deallocation of
            ! mtime_date_container_a needs to be framed with
            ! 1.
            ! l_assoc = ASSOCIATED(mtime_date_container, mtime_date_container_a)
            ! ...
            CALL MOVE_ALLOC(tmp, mtime_date_container_a)
            ! ... and
            ! 2.
            ! mtime_date_container => mtime_date_container_a
          ELSE
            mtime_date_container_a(1:remaining_event_steps) = mtime_date_uniq(1:remaining_event_steps)
          ENDIF
          
          n_event_steps_a = remaining_event_steps
          
          DEALLOCATE(mtime_date_uniq)
        END IF
      ENDIF
      
    END DO INTERVAL_LOOP

    ! copy back results into original data structures
    n_event_steps = n_event_steps_a
    ! to prevent a potential reallocation in next step add 1 element
    ALLOCATE(mtime_dates(n_event_steps+1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! copy mtime_date_container_a to mtime_date_string 
    DO i = 1, n_event_steps
      CALL getDatetimeFromJulianDay(mtime_date_container_a(i), mtime_dates(i))
    END DO

    ! Optional: Append the last event time step
    IF (n_event_steps > 0) THEN
      IF (evd%l_output_last .AND. mtime_dates(n_event_steps) < sim_end .AND. mtime_end >= sim_end) THEN
        ! check, that we do not duplicate the last time step:
        l_append_step = mtime_dates(n_event_steps) /= sim_end
        IF (l_append_step) THEN
          n_event_steps = n_event_steps + 1
          old_size = SIZE(mtime_dates)
          IF (n_event_steps > old_size) THEN
            ! resize buffer
            new_size = old_size + INITIAL_NEVENT_STEPS
            ALLOCATE(tmp_dates(new_size), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            tmp_dates(1:old_size) = mtime_dates
            CALL MOVE_ALLOC(tmp_dates, mtime_dates)
          END IF
          mtime_dates(n_event_steps) = sim_end
          IF (ldebug) THEN
            CALL datetimeToString(sim_end, dtime_string)
            WRITE (0,*) get_my_global_mpi_id(), ": ", &
                 &      n_event_steps, ": output event '", dtime_string, "'"
          END IF
        END IF
      END IF
    END IF
    
    ALLOCATE(mtime_sim_steps(SIZE(mtime_dates)),   &
         &   mtime_exactdate(SIZE(mtime_dates)),   &
         &   filename_metadata(SIZE(mtime_dates)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    CALL fct_time2simstep(n_event_steps, mtime_dates, &
      &                   evd%sim_step_info, mtime_sim_steps, mtime_exactdate)

    ! remove all those event steps which have no corresponding simulation
    ! step (mtime_sim_steps(i) < 0):
    DO i=1,n_event_steps
      IF (mtime_sim_steps(i) < 0)  EXIT
    END DO
    n_event_steps = (i-1)
    
    IF (n_event_steps > 0) THEN
      CALL fct_generate_filenames(n_event_steps, mtime_dates,       &
        &                   mtime_sim_steps, evd%sim_step_info, evd%fname_metadata, &
        &                   skipped_dates, filename_metadata)
    END IF
    
    ! from this list of time stamp strings: generate the event steps
    ! for this event
    p_event%n_event_steps = n_event_steps
    ALLOCATE(p_event%event_step(n_event_steps), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    DO i=1,n_event_steps
      p_event%event_step(i)%exact_date         = mtime_exactdate(i)
      p_event%event_step(i)%i_sim_step         = mtime_sim_steps(i)
      p_event%event_step(i)%n_pes              = 1
      ALLOCATE(step_data(1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      step_data(1)%i_pe            = i_pe
      step_data(1)%datetime = mtime_dates(i)
      step_data(1)%i_tag           = evd%i_tag
      step_data(1)%filename_string = filename_metadata(i)%filename_string
      step_data(1)%jfile           = filename_metadata(i)%jfile
      step_data(1)%jpart           = filename_metadata(i)%jpart
      step_data(1)%l_open_file     = filename_metadata(i)%l_open_file
      CALL MOVE_ALLOC(step_data, p_event%event_step(i)%event_step_data)
    END DO
    IF (ldebug .AND. n_event_steps > 0) THEN
      WRITE (0,'(3a,2(a,i0))') routine, ": defined event ", &
           & TRIM(p_event%event_step(1)%event_step_data(1)%filename_string), &
           & "; tag = ", p_event%event_step(1)%event_step_data(1)%i_tag, &
           & "; pe = ", i_pe
    END IF

    ! clean up
    CALL deallocateDatetime(mtime_end)

    DEALLOCATE(mtime_sim_steps, &
      &        mtime_exactdate, filename_metadata, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  CONTAINS

    SUBROUTINE merge2SortedAndRemoveDuplicates(InputArray1, nsize_IA1, &
         &                                     InputArray2, nsize_IA2, &
         &                                     OutputArray, nsize_OA)
      TYPE(julianday), INTENT(in) :: InputArray1(:)
      TYPE(julianday), INTENT(in) :: InputArray2(:)
      TYPE(julianday), ALLOCATABLE, INTENT(out) :: OutputArray(:)
      INTEGER, INTENT(in) :: nsize_IA1
      INTEGER, INTENT(in) :: nsize_IA2
      INTEGER, INTENT(out) :: nsize_OA
      
      INTEGER(i8) :: diff
      
      INTEGER :: n, na, nb
      INTEGER :: i, j, k


      IF (nsize_IA1 == 0 .AND. nsize_IA2 == 0) THEN
        CALL finish('merge2SortedAndRemoveDuplicates', 'Invalid size of both input arrays.')
      ENDIF

      IF (nsize_IA1 == 0) THEN
        nsize_OA = nsize_IA2
        ALLOCATE(OutputArray(nsize_IA2))
        OutputArray(:) = InputArray2(:nsize_IA2)
        RETURN
      ENDIF
      IF (nsize_IA2 == 0) THEN
        nsize_OA = nsize_IA1
        ALLOCATE(OutputArray(nsize_IA1))
        OutputArray(:) = InputArray1(:nsize_IA1)
        RETURN
      ENDIF
      
      na = nsize_IA1
      nb = nsize_IA2 
      n = na + nb
      
      ALLOCATE(OutputArray(n))
      
      i = 1
      j = 1

      ! handle special case for k == 1 to be Fortran conforming in start-up step

      diff = 86400000_i8 * (InputArray1(i)%day - InputArray2(j)%day) + InputArray1(i)%ms - InputArray2(j)%ms
      IF (diff < 0_i8) THEN
        OutputArray(1) = InputArray1(1)
        i = 2
        j = 1
      ELSE IF (diff > 0_i8) THEN
        OutputArray(1) = InputArray2(1)
        i = 1
        j = 2
      ELSE
        OutputArray(1) = InputArray1(1)
        i = 2
        j = 2
      ENDIF
      k = 2

      DO WHILE (i <= na .AND. j <= nb)
        diff = 86400000_i8 * (InputArray1(i)%day - InputArray2(j)%day) + InputArray1(i)%ms - InputArray2(j)%ms
        IF (diff < 0_i8) THEN
          IF ((InputArray1(i)%day /= OutputArray(k-1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k-1)%ms)) THEN
            OutputArray(k) = InputArray1(i)
            k = k+1
          ENDIF
          i = i+1
        ELSE IF (diff > 0_i8) THEN
          IF ((InputArray2(j)%day /= OutputArray(k-1)%day) .OR. (InputArray2(j)%ms /= OutputArray(k-1)%ms)) THEN
            OutputArray(k) = InputArray2(j)
            k = k+1
          ENDIF
          j = j+1
        ELSE
          IF ((InputArray1(i)%day /= OutputArray(k-1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k-1)%ms)) THEN
            OutputArray(k) = InputArray1(i)
            k = k+1
          ENDIF
          i = i+1
          j = j+1
        ENDIF
      ENDDO
      
      DO WHILE (i <= na)
        IF ((InputArray1(i)%day /= OutputArray(k-1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k-1)%ms)) THEN
          OutputArray(k) = InputArray1(i)
          k = k+1
          i = i+1
        ELSE
          i = i+1
        ENDIF
      ENDDO
      
      DO WHILE (j <= nb)
        IF ((InputArray2(j)%day /= OutputArray(k-1)%day) .OR. (InputArray2(j)%ms /= OutputArray(k-1)%ms)) THEN
          OutputArray(k) = InputArray2(j)
          k = k+1
          j = j+1
        ELSE
          j = j+1
        ENDIF
      ENDDO
      
      nsize_OA = k-1
      
    END SUBROUTINE merge2SortedAndRemoveDuplicates

    SUBROUTINE remove_duplicate_intervals(starts, ends, intvls, n, indices_to_use, remaining)
      CHARACTER(len=*), INTENT(in) :: starts(:)
      CHARACTER(len=*), INTENT(in) :: ends(:)
      CHARACTER(len=*), INTENT(in) :: intvls(:)
      INTEGER, INTENT(in) :: n         
      INTEGER, INTENT(out) :: indices_to_use(n)
      INTEGER, INTENT(out) :: remaining
      
      INTEGER :: i, j
      
      remaining = 1
      indices_to_use(1) = 1
      
      OUTER_LOOP: DO i = 2, n
        DO j = 1, remaining
          IF (TRIM(starts(j)) == TRIM(starts(i))) THEN
            IF (TRIM(ends(j)) == TRIM(ends(i))) THEN
              IF (TRIM(intvls(j)) == TRIM(intvls(i))) THEN
                ! found a match so start looking again
                CYCLE OUTER_LOOP
              END IF
            END IF
          END IF
        END DO
        ! no match found so add it to the output
        remaining = remaining + 1
        indices_to_use(remaining) = i
      END DO OUTER_LOOP
      
    END SUBROUTINE remove_duplicate_intervals
    
  END SUBROUTINE init_output_event


  !> Create a simple *parallel* output event, happening at regular intervals.
  !
  !  This subroutine calls the local version "new_output_event", adds
  !  data structures for parallel communication and launches a
  !  non-blocking MPI send to the root I/O PE.
  !
  !  The begin and end time stamps may contain "modifier symbols",
  !  e.g. ">2014-06-01T00:00:00.000", where ">" means that output
  !  should not include the start date.
  !
  FUNCTION new_parallel_output_event(name, begin_str, end_str, intvl_str,                               &
    &                                l_output_last, sim_step_info, fname_metadata, fct_time2simstep,    &
    &                                fct_generate_filenames, local_event_no, icomm,                     &
    &                                event_list_local, ievent_list_local) RESULT(p_event)
    TYPE(t_par_output_event), POINTER :: p_event
    CHARACTER(LEN=*),                      INTENT(IN)  :: name                 !< output event name
    !> start time stamp + modifier
    CHARACTER(len=MAX_DATETIME_STR_LEN+1), INTENT(IN)  :: begin_str(MAX_TIME_INTERVALS)
    !> start time stamp + modifier
    CHARACTER(len=MAX_DATETIME_STR_LEN+1), INTENT(IN)  :: end_str(MAX_TIME_INTERVALS)
    CHARACTER(len=MAX_DATETIME_STR_LEN),   INTENT(IN)  :: intvl_str(MAX_TIME_INTERVALS)
    !> Flag. If .TRUE. the last step is always written
    LOGICAL,                               INTENT(IN)  :: l_output_last
    !> definitions for conversion "time stamp -> simulation step"
    TYPE(t_sim_step_info),                 INTENT(IN)  :: sim_step_info
    !> additional meta-data for generating output filename
    TYPE(t_fname_metadata),                INTENT(IN)  :: fname_metadata
    INTEGER,                               INTENT(IN)  :: local_event_no       !< local index of this event on local PE
    INTEGER,                               INTENT(IN)  :: icomm                !< MPI communicator
    !> local list of output events
    TYPE(t_event_data_local),              INTENT(INOUT), CONTIGUOUS :: event_list_local(:)
    !> length of local list of output events
    INTEGER,                               INTENT(INOUT) :: ievent_list_local


    INTERFACE
      !> As an argument of this function, the user must provide a
      !! conversion "time stamp -> simulation step"
      SUBROUTINE fct_time2simstep(num_dates, dates, sim_step_info, &
        &                         result_steps, result_exactdate)
        IMPORT :: t_sim_step_info, datetime
        !> no. of dates to convert
        INTEGER,                   INTENT(IN)    :: num_dates
        !> array of mtime datetime objects
        TYPE(datetime),            INTENT(IN)    :: dates(:)
        TYPE(t_sim_step_info),     INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,                   INTENT(OUT)   :: result_steps(:)      !< resulting step indices
        TYPE(datetime),            INTENT(OUT)   :: result_exactdate(:)  !< resulting (exact) time step
      END SUBROUTINE fct_time2simstep

      !> As an argument of this function, the user must provide a
      !! function for generating output file names
      SUBROUTINE fct_generate_filenames(num_dates, dates, sim_steps, &
        &                               sim_step_info, fname_metadata, &
        &                               skipped_dates, result_fnames)
        IMPORT :: t_sim_step_info, t_event_step_data, t_fname_metadata, &
             datetime

        !> no. of dates to convert
        INTEGER,                   INTENT(IN)    :: num_dates
        !> array of mtime datetime objects
        TYPE(datetime), TARGET,    INTENT(IN)    :: dates(:)
        INTEGER,                   INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),     INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata),    INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        INTEGER,                   INTENT(IN)    :: skipped_dates
        TYPE(t_event_step_data),   INTENT(OUT)   :: result_fnames(SIZE(dates))
      END SUBROUTINE fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::new_parallel_output_event"
    INTEGER :: ierrstat, this_pe, i_tag, nranks, nintvls

    ! determine this PE's MPI rank wrt. the given MPI communicator:
    IF (icomm /= mpi_comm_null) THEN
      this_pe = p_comm_rank(icomm)
      nranks = p_comm_size(icomm)
      IF (ldebug) THEN
        WRITE (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", this_pe, "; icomm has size ", nranks
      END IF
    ELSE
      this_pe = 0
      nranks  = 1
    END IF

    ! compute i_tag ID st. it stays unique even for multiple events
    ! running on the same I/O PE:
    i_tag = SENDRECV_TAG_OUTEVENT + this_pe + (local_event_no-1)*nranks

    ! allocate parallel event data structure
    ALLOCATE(p_event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    NULLIFY(p_event%next)

    ! count the number of different time intervals for this event (usually 1)
    DO nintvls = 0, MAX_TIME_INTERVALS-1
      IF (LEN_TRIM(begin_str(nintvls+1)) == 0) EXIT
    END DO

    ievent_list_local = ievent_list_local + 1
    event_list_local(ievent_list_local)%name = name
    event_list_local(ievent_list_local)%begin_str = begin_str
    event_list_local(ievent_list_local)%end_str = end_str
    event_list_local(ievent_list_local)%intvl_str = intvl_str
    event_list_local(ievent_list_local)%l_output_last = l_output_last
    event_list_local(ievent_list_local)%sim_step_info = sim_step_info
    event_list_local(ievent_list_local)%fname_metadata = fname_metadata
    event_list_local(ievent_list_local)%i_tag = i_tag
    ! create the local (non-parallel) output event data structure:
    CALL init_output_event(p_event%output_event, &
      &                    event_list_local(ievent_list_local), this_pe, &
      &                    fct_time2simstep, fct_generate_filenames)

    IF (ldebug) THEN
      WRITE (0,*) "PE ", get_my_global_mpi_id(), ": created event with ", &
        &      p_event%output_event%n_event_steps, " steps."
    END IF

    ! set the other MPI-related data fields:
    p_event%icomm         = icomm
    p_event%iroot         = 0
    p_event%irecv_nreq    = 0
    p_event%isend_req     = 0
#ifndef NOMPI
    p_event%isend_req     = MPI_REQUEST_NULL
#endif
  END FUNCTION new_parallel_output_event


  !---------------------------------------------------------------
  ! ROUTINES FOR JOINING OUTPUT EVENTS
  !---------------------------------------------------------------

  !> Receives event meta-data from all PEs within the given MPI
  !> communicator; creates the union of these events.
  !
  FUNCTION union_of_all_events(fct_time2simstep, fct_generate_filenames, &
    &                          icomm, root_outevent, &
    &                          event_list_local, ievent_list_local)
    TYPE(t_par_output_event), POINTER :: union_of_all_events
    !> MPI communicator for intra-I/O communication
    INTEGER,                INTENT(IN)  :: icomm, root_outevent
    TYPE(t_event_data_local), INTENT(INOUT), CONTIGUOUS :: event_list_local(:)
    !> length of local list of output events
    INTEGER, INTENT(in) :: ievent_list_local

    INTERFACE
      !> As an argument of this function, the user must provide a
      !! conversion "time stamp -> simulation step"
      SUBROUTINE fct_time2simstep(num_dates, dates, sim_step_info, &
        &                         result_steps, result_exactdate)
        IMPORT :: t_sim_step_info, datetime
        !> no. of dates to convert
        INTEGER,                  INTENT(IN)    :: num_dates
        !> array of mtime datetime objects
        TYPE(datetime),           INTENT(IN)    :: dates(:)
        TYPE(t_sim_step_info),    INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,                  INTENT(OUT)   :: result_steps(:)      !< resulting step indices
        TYPE(datetime),           INTENT(OUT)   :: result_exactdate(:)  !< resulting (exact) time step
      END SUBROUTINE fct_time2simstep

      !> As an argument of this function, the user must provide a
      !! function for generating output file names
      SUBROUTINE fct_generate_filenames(num_dates, dates, sim_steps, &
        &                               sim_step_info, fname_metadata, &
        &                               skipped_dates, result_fnames)
        IMPORT :: t_sim_step_info, t_event_step_data, t_fname_metadata, &
             datetime

        !> no. of dates to convert
        INTEGER,                   INTENT(IN)    :: num_dates
        !> array of mtime datetime objects
        TYPE(datetime), TARGET,    INTENT(IN)    :: dates(:)
        INTEGER,                   INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),     INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata),    INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        INTEGER,                   INTENT(IN)    :: skipped_dates
        TYPE(t_event_step_data),  INTENT(OUT)    :: result_fnames(SIZE(dates))
      END SUBROUTINE fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::union_of_all_events"
    TYPE(t_par_output_event), POINTER     :: par_event, last_node
    TYPE(t_output_event), ALLOCATABLE     :: ev2
    INTEGER                               :: i_pe, nranks, ierror, &
      &                                      this_pe, i, acc,   &
      &                                      num_all_events, cnt
    INTEGER, ALLOCATABLE                  :: evd_counts(:), evd_rank(:)
    TYPE(t_event_data_local), ALLOCATABLE :: evd_all(:)

    IF (ldebug)  WRITE (0,*) routine, " enter."

    NULLIFY(union_of_all_events)

    IF (icomm /= mpi_comm_null) THEN
      this_pe = p_comm_rank(icomm)
      nranks = p_comm_size(icomm)
      IF (ldebug) THEN
        WRITE (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", &
             this_pe, "; icomm has size ", nranks
      END IF
#ifndef NOMPI
      IF (pio_type == pio_type_cdipio) THEN
        IF (this_pe == root_outevent) THEN
          CALL p_recv(num_all_events, 0, p_tag=156, comm=icomm)
          ALLOCATE(evd_all(num_all_events), evd_rank(num_all_events), STAT=ierror)
          IF (ierror /= 0) CALL finish (routine, 'ALLOCATE failed.')
          CALL p_recv(evd_all, 0, p_tag=156, comm=icomm)
          evd_rank = 0
        ELSE IF (p_pe == p_work_pe0) THEN
          CALL p_send(ievent_list_local, p_io_pe0, p_tag=156, comm=icomm)
          CALL p_send(event_list_local, p_io_pe0, p_tag=156, comm=icomm)
          num_all_events = 0
        ELSE
          num_all_events = 0
        END IF
      ELSE
#endif
        CALL p_gatherv(event_list_local(1:ievent_list_local), evd_all, &
          &            root_outevent, counts=evd_counts, comm=icomm)
        IF (this_pe == root_outevent) THEN
          num_all_events = SIZE(evd_all)
          ALLOCATE(evd_rank(num_all_events), STAT=ierror)
          IF (ierror /= 0) CALL finish (routine, 'ALLOCATE failed.')
          acc = 0
          DO i = 1, nranks
            cnt = evd_counts(i)
            evd_rank(acc+1:acc+cnt) = i-1
            acc = acc + cnt
          END DO
        ELSE
          num_all_events = 0
        END IF
#ifndef NOMPI
      END IF
#endif

      IF (num_all_events > 0) THEN
        ! create the event steps from the received meta-data:
        ALLOCATE(union_of_all_events, STAT=ierror)
        IF (ierror /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
        par_event => union_of_all_events
        CALL init_output_event(par_event%output_event, &
          &                    evd_all(1), evd_rank(1), &
          &                    fct_time2simstep, fct_generate_filenames)
        par_event%icomm         = icomm
        par_event%iroot         = root_outevent
        par_event%irecv_nreq    = 0
        NULLIFY(par_event%next)
        IF (ldebug) THEN
          WRITE (0,*) "PE ", get_my_global_mpi_id(), ": n_event_steps = ", &
            & union_of_all_events%output_event%n_event_steps
        END IF
        event_append_loop: DO i = 2, num_all_events
          i_pe = evd_rank(i)
          ! create the event steps from the received meta-data:
          CALL init_output_event(ev2, evd_all(i), i_pe, &
            &                    fct_time2simstep, fct_generate_filenames)
          ! find the parallel output event in the linked list that
          ! matches the event name
          NULLIFY(last_node)
          par_event => union_of_all_events
          DO WHILE (ASSOCIATED(par_event))
            IF (par_event%output_event%event_data%name == ev2%event_data%name) THEN
              CALL merge_events(par_event%output_event, ev2)
              CYCLE event_append_loop
            END IF

            last_node => par_event
            par_event => par_event%next
          END DO
          ! if there is no such parallel output event, then create one
          ! at the end of the linked list:
          ALLOCATE(last_node%next, STAT=ierror)
          IF (ierror /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
          par_event => last_node%next
          ! set the MPI-related data fields:
          par_event%icomm         = icomm
          par_event%iroot         = root_outevent
          par_event%irecv_nreq    = 0
          NULLIFY(par_event%next)
          CALL MOVE_ALLOC(ev2, par_event%output_event)
        END DO event_append_loop
      END IF
    END IF

    IF (ldebug)  WRITE (0,*) routine, " done."
  END FUNCTION union_of_all_events


  !> Merge compatible event into another.
  !
  !  Different output events/PEs may be joined together to form a
  !  single new event. Then not every PE necessarily participates in
  !  every event step. The purpose of joining/grouping events is that
  !  certain actions can be triggered for the group, e.g. writing
  !  "ready files".
  !
  SUBROUTINE merge_events(event1, event2)
    TYPE(t_output_event), INTENT(INOUT) :: event1
    TYPE(t_output_event), INTENT(IN) :: event2
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::merge_event"
    INTEGER                             :: i1, i2, ierrstat, i_sim_step1, i_sim_step2, &
      &                                    max_sim_step, j1, j2, nsteps
    LOGICAL :: copy_i1
    CHARACTER(LEN=MAX_FILENAME_STR_LEN) :: filename_string1
    TYPE(t_event_step), ALLOCATABLE :: new_steps(:)
    ! allocate event data structure
    event1%i_event_step             = 1
    IF (event1%event_data%name /= event2%event_data%name) THEN
      i1 = LEN_TRIM(event1%event_data%name)
      event1%event_data%name(i1:) = ","//event2%event_data%name
    END IF
    IF (event1%event_data%sim_start /= event2%event_data%sim_start) THEN
      CALL finish(routine, "Simulation start dates do not match!")
    END IF

    ! consistency check: test, if any of the filenames in event1
    ! occurs in event2 (avoid duplicate names)
    DO i1=1,event1%n_event_steps
      DO j1=1,event1%event_step(i1)%n_pes
        IF (.NOT. event1%event_step(i1)%event_step_data(j1)%l_open_file) CYCLE
        filename_string1 = event1%event_step(i1)%event_step_data(j1)%filename_string
        DO i2=1,event2%n_event_steps
          DO j2=1,event2%event_step(i2)%n_pes
            IF (.NOT. event2%event_step(i2)%event_step_data(j2)%l_open_file) CYCLE
            IF (event2%event_step(i2)%event_step_data(j2)%filename_string == filename_string1) THEN
              ! found a duplicate filename:
              CALL finish(routine, "Error! Ambiguous output file name: '"//TRIM(filename_string1)//"'")
            END IF
          END DO
        END DO
      END DO
    END DO

    ! choose a maximum simulation step number as abort criterion:
    max_sim_step = 0
    IF (event1%n_event_steps > 0)  max_sim_step = MAX(max_sim_step, event1%event_step(event1%n_event_steps)%i_sim_step)
    IF (event2%n_event_steps > 0)  max_sim_step = MAX(max_sim_step, event2%event_step(event2%n_event_steps)%i_sim_step)
    max_sim_step = max_sim_step + 1

    nsteps = 0
    i1 = 1
    i2 = 1
    DO WHILE (i1 <= event1%n_event_steps .OR. i2 <= event2%n_event_steps)
      IF (i1 > event1%n_event_steps) THEN
        i_sim_step1 = max_sim_step
      ELSE
        i_sim_step1 = event1%event_step(i1)%i_sim_step
      END IF
      IF (i2 > event2%n_event_steps) THEN
        i_sim_step2 = max_sim_step
      ELSE
        i_sim_step2 = event2%event_step(i2)%i_sim_step
      END IF
      nsteps = nsteps+1
      ! copy event step i1 from event1:
      i1 = i1 + MERGE(1, 0, i_sim_step2 >= i_sim_step1)
      ! copy event step i2 from event2:
      i2 = i2 + MERGE(1, 0, i_sim_step2 <= i_sim_step1)
    END DO

    ! now create the new event:
    ALLOCATE(new_steps(nsteps), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    nsteps = 0
    i1 = 1
    i2 = 1
    DO WHILE (i1 <= event1%n_event_steps .OR. i2 <= event2%n_event_steps)
      IF (i1 > event1%n_event_steps) THEN
        i_sim_step1 = max_sim_step
      ELSE
        i_sim_step1 = event1%event_step(i1)%i_sim_step
      END IF
      IF (i2 > event2%n_event_steps) THEN
        i_sim_step2 = max_sim_step
      ELSE
        i_sim_step2 = event2%event_step(i2)%i_sim_step
      END IF
      nsteps = nsteps+1
      copy_i1 = i_sim_step1 <= i_sim_step2
      IF (copy_i1) THEN
        ! copy event step i1 from event1:
        new_steps(nsteps)%i_sim_step = event1%event_step(i1)%i_sim_step
        new_steps(nsteps)%exact_date = event1%event_step(i1)%exact_date
        new_steps(nsteps)%n_pes = event1%event_step(i1)%n_pes
        CALL MOVE_ALLOC(from=event1%event_step(i1)%event_step_data, &
          &             to=new_steps(nsteps)%event_step_data)
        i1 = i1 + 1
      END IF
      IF (i_sim_step2 <= i_sim_step1) THEN
        ! copy event step i2 from event2:
        CALL append_event_step(new_steps(nsteps), event2%event_step(i2), l_create=(.NOT. copy_i1))
        i2 = i2 + 1
      END IF
    END DO
    CALL MOVE_ALLOC(from=new_steps, to=event1%event_step)
    event1%n_event_steps = nsteps
  END SUBROUTINE merge_events


  !> Appends the data of an event step to the data of another event step.
  !
  SUBROUTINE append_event_step(dst_event_step, src_event_step, l_create)
    TYPE(t_event_step), INTENT(INOUT)  :: dst_event_step  !< source event step
    TYPE(t_event_step), INTENT(IN)     :: src_event_step  !< destination event step
    LOGICAL,            INTENT(IN)     :: l_create        !< Flag. If .TRUE., the the destination event step is created
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::append_event_step"
    TYPE(t_event_step_data), ALLOCATABLE  :: tmp_event_step_data(:)
    INTEGER :: iold, ierrstat, i1, i2, inew, iadd

    IF (l_create) THEN
      dst_event_step%i_sim_step = src_event_step%i_sim_step
      dst_event_step%exact_date = src_event_step%exact_date
    ELSE
      IF (src_event_step%i_sim_step /= dst_event_step%i_sim_step) &
        &   CALL finish(routine, "sim_steps do not match!")
      IF (src_event_step%exact_date /= dst_event_step%exact_date) &
        &   CALL finish(routine, "exact dates do not match!")
    END IF
    ! event_step%event_step_data
    IF (l_create) THEN
      dst_event_step%n_pes = src_event_step%n_pes
      ALLOCATE(dst_event_step%event_step_data(dst_event_step%n_pes), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      dst_event_step%event_step_data(1:dst_event_step%n_pes) = &
        &   src_event_step%event_step_data(1:src_event_step%n_pes)
    ELSE
      iold = dst_event_step%n_pes
      iadd = src_event_step%n_pes
      inew = iold + iadd
      ALLOCATE(tmp_event_step_data(inew), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      tmp_event_step_data(1:iold) = dst_event_step%event_step_data(1:iold)
      tmp_event_step_data(iold+1:inew) = &
           & src_event_step%event_step_data(1:iadd)
      CALL MOVE_ALLOC(from=tmp_event_step_data, &
        &             to=dst_event_step%event_step_data)
      dst_event_step%n_pes = inew

      ! consistency check: test, if any of the filenames in group1
      ! occurs in group2 (avoid duplicate names)
      DO i1=1,iold
        DO i2=(iold+1),(iold+src_event_step%n_pes)
          IF (dst_event_step%event_step_data(i1)%filename_string &
            == dst_event_step%event_step_data(i2)%filename_string) THEN
            ! found a duplicate filename:
            CALL finish(routine, "Ambiguous output file name: '"//TRIM(dst_event_step%event_step_data(i1)%filename_string)//"'")
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE append_event_step


  !---------------------------------------------------------------
  ! ROUTINES FOR TESTING / INQURING OUTPUT EVENTS
  !---------------------------------------------------------------

  !> @return .TRUE. if the output event has exceeded its final step.
  !
  FUNCTION is_output_event_finished(event)
    LOGICAL :: is_output_event_finished
    TYPE(t_output_event), INTENT(IN) :: event
    is_output_event_finished = (event%i_event_step > event%n_event_steps)
  END FUNCTION is_output_event_finished


  !> @return .TRUE. if the output event has exceeded its final step.
  !
  FUNCTION is_par_output_event_finished(event)
    LOGICAL :: is_par_output_event_finished
    TYPE(t_par_output_event), INTENT(in)             :: event

    is_par_output_event_finished = is_output_event_finished(event%output_event)
  END FUNCTION is_par_output_event_finished


  !> @return .TRUE. if the given step index matches the step index of
  !          the current event step, ie. if the current event step is
  !          active.
  !
  FUNCTION is_output_step(event, jstep)
    LOGICAL :: is_output_step
    TYPE(t_output_event), INTENT(IN) :: event  !< output event
    INTEGER,              INTENT(IN) :: jstep  !< given step index
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::is_output_step"
    INTEGER :: istep

    IF (is_output_event_finished(event)) THEN
      is_output_step = .FALSE.
    ELSE
      istep = event%i_event_step
      is_output_step = (event%event_step(istep)%i_sim_step == jstep)
    END IF
  END FUNCTION is_output_step


  !> @return .TRUE. if the given step index matches the step index of
  !          the current event step, ie. if the current event step is
  !          active.
  !
  !  @note This function recursively checks the complete
  !        singly-linked list rooted at @p event.
  !
  FUNCTION is_par_output_step(event, jstep) RESULT(p)
    TYPE(t_par_output_event), POINTER, INTENT(IN) :: event  !< output event
    INTEGER,                           INTENT(IN) :: jstep  !< given step index
    TYPE(t_par_output_event), POINTER :: ev
    ! local variables
    LOGICAL :: p

    p = .FALSE.
    ev => event
    DO WHILE (.NOT. p .AND. ASSOCIATED(ev))
      p = p .OR. is_output_step(ev%output_event, jstep)
      ev => ev%next
    END DO
  END FUNCTION is_par_output_step


  !> @return .TRUE. if all participants of the parallel output event
  !>         have acknowledged the completion of the current step.
  !
  !  @note We do not use the wrapper routines from "mo_mpi" here,
  !        since we need direct control over the non-blocking P2P
  !        requests.
  !
  FUNCTION is_output_step_complete(event) RESULT(ret)
    LOGICAL :: ret
    TYPE(t_par_output_event), INTENT(inout) :: event
    ! local variables
#ifndef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::is_output_step_complete"
    INTEGER              :: ierrstat
    ! INTEGER, ALLOCATABLE :: irecv_status(:,:)

    ret = .TRUE.
    IF (event%irecv_nreq /= 0) THEN
      ! ALLOCATE(irecv_status(MPI_STATUS_SIZE,event%irecv_nreq), STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      CALL MPI_TESTALL(event%irecv_nreq, event%irecv_req, ret, &
        &              mpi_statuses_ignore, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_TESTALL.')
      ! DEALLOCATE(irecv_status, STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
#else
    ret = .TRUE.
#endif
  END FUNCTION is_output_step_complete


  !> @return current date-time stamp string.
  !
  FUNCTION get_current_date(event) RESULT(current_date)
    TYPE(datetime) :: current_date
    TYPE(t_output_event), INTENT(in) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_current_date"
    INTEGER :: istep

    istep = MIN(event%i_event_step, event%n_event_steps)
    ! For events that are the union of multiple other events we return
    ! the *model date-time string*:
    IF (event%event_step(istep)%n_pes > 1) THEN
      current_date = event%event_step(istep)%exact_date
    ELSE
      current_date = event%event_step(istep)%event_step_data(1)%datetime
    END IF
  END FUNCTION get_current_date


  !> @return current date-time stamp string.
  !
  FUNCTION get_current_date_par(event) RESULT(current_date)
    TYPE(datetime) :: current_date
    TYPE(t_par_output_event), INTENT(in) :: event
    current_date =  get_current_date(event%output_event)
  END FUNCTION get_current_date_par


  !> @return current output step.
  !
  FUNCTION get_current_step(event)
    INTEGER :: get_current_step
    TYPE(t_output_event), INTENT(in) :: event
    get_current_step = MIN(event%i_event_step, event%n_event_steps)
  END FUNCTION get_current_step


  !> @return current output step.
  !
  FUNCTION get_current_step_par(event)
    INTEGER :: get_current_step_par
    TYPE(t_par_output_event), POINTER :: event
    IF (.NOT. ASSOCIATED(event)) THEN
      get_current_step_par = -1
    ELSE
      get_current_step_par =  get_current_step(event%output_event)
    END IF
  END FUNCTION get_current_step_par


  !> @return current filename string.
  !
  FUNCTION get_current_filename(event)
    CHARACTER(LEN=MAX_FILENAME_STR_LEN) :: get_current_filename
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_current_filename"
    INTEGER :: istep
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: event_datetime_string

    istep = event%output_event%i_event_step
    ! We cannot check events that are the union of multiple other
    ! events:
    IF (event%output_event%event_step(istep)%n_pes > 1) THEN
      CALL datetimetostring(event%output_event%event_step(istep)%exact_date, &
        &                   event_datetime_string)
      WRITE (0,*) "Event step ", istep,  "(", TRIM(event_datetime_string), ")", &
        &         ", shared by ", event%output_event%event_step(istep)%n_pes, " PEs."
      CALL finish(routine, "Error! Multi-part event step!")
    END IF
    get_current_filename = event%output_event%event_step(istep)%event_step_data(1)%filename_string
  END FUNCTION get_current_filename


  !> @return current file number.
  !
  FUNCTION get_current_jfile(event)
    INTEGER :: get_current_jfile
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_current_jfile"
    INTEGER :: istep
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: event_datetime_string

    istep = MIN(event%output_event%i_event_step, event%output_event%n_event_steps)
    IF (istep == 0) THEN
      get_current_jfile = 0
      RETURN
    END IF
    ! We cannot check events that are the union of multiple other
    ! events:
    IF (event%output_event%event_step(istep)%n_pes > 1) THEN
      CALL datetimetostring(event%output_event%event_step(istep)%exact_date, &
        &                   event_datetime_string)
      WRITE (0,*) "Event step ", istep,  "(", TRIM(event_datetime_string), ")", &
        &         ", shared by ", event%output_event%event_step(istep)%n_pes, " PEs."
      CALL finish(routine, "Error! Multi-part event step!")
    END IF
    get_current_jfile = event%output_event%event_step(istep)%event_step_data(1)%jfile
  END FUNCTION get_current_jfile


  !> @return .TRUE. if this PE should write a ready file for the given
  !          event.
  !
  FUNCTION check_write_readyfile(event) RESULT(do_write)
    LOGICAL :: do_write
    TYPE(t_output_event), INTENT(in) :: event

    do_write =       (event%event_data%name /= DEFAULT_EVENT_NAME) &
      &        .AND. (event%n_event_steps > 0)
  END FUNCTION check_write_readyfile


  !> @return .TRUE. if current event step has the "open file" flag
  !          enabled.
  !
  FUNCTION check_open_file(event)
    LOGICAL :: check_open_file
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::check_open_file"
    INTEGER :: istep
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: event_datetime_string

    istep = event%output_event%i_event_step
    ! We cannot check events that are the union of multiple other
    ! events:
    IF (event%output_event%event_step(istep)%n_pes > 1) THEN
      CALL datetimetostring(event%output_event%event_step(istep)%exact_date, &
        &                   event_datetime_string)
      WRITE (0,*) "Event step ", istep,  "(", TRIM(event_datetime_string), ")",&
        &         ", shared by ", event%output_event%event_step(istep)%n_pes,  &
        &         " PEs."
      CALL finish(routine, "Error! Multi-part event step!")
    END IF
    check_open_file = event%output_event%event_step(istep)%event_step_data(1)%l_open_file
  END FUNCTION check_open_file


  !> @return .TRUE. if current event step has the "close file" flag
  !          enabled.
  !
  FUNCTION check_close_file(event)
    LOGICAL :: check_close_file
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::check_close_file"

    check_close_file = check_close_file_step(event, event%output_event%i_event_step)
  END FUNCTION check_close_file

  !> @return .TRUE. if current event step has the "close file" flag
  !          enabled.
  !
  FUNCTION check_close_file_step(event, istep) RESULT(check_close_file)
    LOGICAL :: check_close_file
    TYPE(t_par_output_event), INTENT(IN) :: event
    INTEGER, INTENT(in) :: istep
    ! local variables
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: event_datetime_string
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::check_close_file"

    ! We cannot check events that are the union of multiple other
    ! events:
    IF (istep <= SIZE(event%output_event%event_step)) THEN
      IF (event%output_event%event_step(istep)%n_pes > 1) THEN
        CALL datetimetostring(event%output_event%event_step(istep)%exact_date, &
          &                   event_datetime_string)
        WRITE (0,*) "Event step ", istep,  "(", TRIM(event_datetime_string), ")", &
        &         ", shared by ", event%output_event%event_step(istep)%n_pes, " PEs."
        CALL finish(routine, "Error! Multi-part event step!")
      END IF
      check_close_file &
        = event%output_event%event_step(istep)%event_step_data(1)%l_open_file
    ELSE
      check_close_file = .TRUE.
    END IF
    check_close_file = check_close_file .AND. istep > 1
  END FUNCTION check_close_file_step


  !> @return .TRUE. if this PE is the event's root PE.  If no event
  !          data structure has been provided as an argument to this
  !          function, we test for ROOT_OUTEVENT.
  !
  FUNCTION is_event_root_pe(opt_event, opt_icomm)
    LOGICAL :: is_event_root_pe
    TYPE(t_par_output_event), INTENT(inout), OPTIONAL :: opt_event
    INTEGER, INTENT(IN),               OPTIONAL :: opt_icomm
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::is_event_root_pe"
#ifndef NOMPI
    INTEGER :: this_pe, p_comm
#endif

    ! at least one of the two arguments must be provided:
    IF (.NOT. PRESENT(opt_event) .AND. .NOT. PRESENT(opt_icomm)) THEN
      CALL finish(routine, "Internal error!")
    END IF

#ifndef NOMPI
    IF (PRESENT(opt_event)) THEN
      p_comm = opt_event%icomm
    ELSE
      p_comm = opt_icomm
    END IF
    ! determine this PE's MPI rank wrt. the given MPI communicator and
    ! return if this PE is not root PE for the given event:
    this_pe = p_comm_rank(p_comm)
    IF (PRESENT(opt_event)) THEN
      is_event_root_pe = (this_pe == opt_event%iroot)
    ELSE
      is_event_root_pe = (this_pe == 0)
    END IF
#else
    is_event_root_pe = .TRUE.
#endif
  END FUNCTION is_event_root_pe


  !---------------------------------------------------------------
  ! ROUTINES PERFORMING DATA TRANSFER TO ROOT PE DURING SETUP
  !---------------------------------------------------------------
  SUBROUTINE p_gatherv_event_data_1d1d(sbuf, rbuf, p_dest, counts, comm)
    TYPE(t_event_data_local), INTENT(in), CONTIGUOUS :: sbuf(:)
    TYPE(t_event_data_local), ALLOCATABLE, INTENT(inout) :: rbuf(:)
    INTEGER, ALLOCATABLE, INTENT(inout) :: counts(:)
    INTEGER, INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

    INTEGER :: sum_counts
#ifndef NOMPI
    CHARACTER(len=*), PARAMETER :: &
         routine = modname//"::p_gatherv_event_data_1d1d"
    INTEGER, ALLOCATABLE :: displs(:)
    INTEGER :: comm_rank, comm_size, p_comm, acc, i, ierror
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF
    comm_size = p_comm_size(p_comm)
    comm_rank = p_comm_rank(p_comm)
    IF (comm_rank == p_dest) THEN
      IF (ALLOCATED(counts)) THEN
        IF (SIZE(counts) /= comm_size) DEALLOCATE(counts)
      END IF
      IF (.NOT. ALLOCATED(counts)) THEN
        ALLOCATE(counts(comm_size))
        counts(1) = -1
      END IF
      ALLOCATE(displs(comm_size))
    ELSE
      IF (ALLOCATED(counts)) THEN
        IF (SIZE(counts) < 1) DEALLOCATE(counts)
      END IF
      IF (.NOT. ALLOCATED(counts)) THEN
        ALLOCATE(counts(1))
        counts(1) = -1
      END IF
      ALLOCATE(displs(1))
    END IF
    IF (counts(1) < 0) CALL p_gather(SIZE(sbuf), counts, p_dest, p_comm)
    IF (comm_rank == p_dest) THEN
      acc = 0
      DO i = 1, comm_size
        displs(i) = acc
        acc = acc + counts(i)
      END DO
      sum_counts = acc
    ELSE
      sum_counts = 1
    END IF
    IF (ALLOCATED(rbuf)) THEN
      IF (SIZE(rbuf) < sum_counts) DEALLOCATE(rbuf)
    END IF
    IF (.NOT. ALLOCATED(rbuf)) ALLOCATE(rbuf(sum_counts))
    IF (event_data_dt == mpi_datatype_null) CALL create_event_data_dt
    CALL mpi_gatherv(sbuf, SIZE(sbuf), event_data_dt, &
      &              rbuf, counts, displs, event_data_dt, &
      &              p_dest, p_comm, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_gatherv error')
#else
    sum_counts = SIZE(sbuf)
    IF (ALLOCATED(counts)) THEN
      IF (SIZE(counts) /= 1) DEALLOCATE(counts)
    END IF
    IF (.NOT. ALLOCATED(counts)) ALLOCATE(counts(1))
    counts(1) = sum_counts
    IF (ALLOCATED(rbuf)) THEN
      IF (SIZE(rbuf) < sum_counts) DEALLOCATE(rbuf)
    END IF
    IF (.NOT. ALLOCATED(rbuf)) ALLOCATE(rbuf(sum_counts))
    rbuf(1:sum_counts) = sbuf
#endif
  END SUBROUTINE p_gatherv_event_data_1d1d

#ifndef NOMPI
  SUBROUTINE p_send_event_data_1d(t_buffer, p_destination, p_tag, p_count, comm)
    TYPE(t_event_data_local), INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER :: p_comm, icount, ierror
    CHARACTER(len=*), PARAMETER :: routine = modname//'::p_send_event_data_1d'

    IF (PRESENT(comm)) THEN
      p_comm = comm
    ELSE
      p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
      icount = p_count
    ELSE
      icount = SIZE(t_buffer)
    END IF

    IF (event_data_dt == mpi_datatype_null) CALL create_event_data_dt
    CALL mpi_send(t_buffer, icount, event_data_dt, p_destination, p_tag, &
            p_comm, ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (message_text,'(2(a,i4),a,i6,2a,i0)') 'mpi_send from ', p_pe, &
         ' to ', p_destination, ' for tag ', p_tag, ' failed.', &
         ' error code=', ierror
       CALL finish(routine, message_text)
    END IF
  END SUBROUTINE p_send_event_data_1d

  SUBROUTINE p_recv_event_data_1d(t_buffer, p_source, p_tag, p_count, comm)
    TYPE(t_event_data_local), INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER :: p_comm, icount, ierror
    CHARACTER(len=*), PARAMETER :: routine = modname//'::p_recv_event_data_1d'

    IF (PRESENT(comm)) THEN
      p_comm = comm
    ELSE
      p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
      icount = p_count
    ELSE
      icount = SIZE(t_buffer)
    END IF

    IF (event_data_dt == mpi_datatype_null) CALL create_event_data_dt
    CALL mpi_recv(t_buffer, icount, event_data_dt, p_source, p_tag, &
            p_comm, mpi_status_ignore, ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (message_text,'(2(a,i4),a,i6,2a,i0)') 'mpi_recv from ', p_source, &
         ' to ', p_pe, ' for tag ', p_tag, ' failed.', &
         ' error code=', ierror
       CALL finish(routine, message_text)
    END IF
  END SUBROUTINE p_recv_event_data_1d


  !> create mpi datatype for variables of type t_event_data_local
  SUBROUTINE create_event_data_dt
    INTEGER, PARAMETER :: num_dt_elem = 8
    INTEGER(mpi_address_kind) :: base, displs(num_dt_elem), ext, stride
    INTEGER :: elem_dt(num_dt_elem), i, ierror, resized_dt
    CHARACTER(len=*), PARAMETER :: routine = modname//"::create_event_data_dt"
    INTEGER, PARAMETER :: blocklens(num_dt_elem) &
      = (/ 1, max_time_intervals, max_time_intervals, max_time_intervals, &
      &    1, 1, 1, 1 /)
    TYPE(t_event_data_local), TARGET :: dummy(2)
    TYPE(datetime), POINTER :: temp_datetime
    CALL mpi_type_contiguous(max_event_name_str_len, mpi_character, &
      &                      elem_dt(1), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_contiguous error')
    CALL mpi_type_contiguous(max_datetime_str_len, mpi_character, &
      &                      elem_dt(2), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_contiguous error')
    elem_dt(3) = elem_dt(2)
    CALL mpi_type_contiguous(max_timedelta_str_len, mpi_character, &
      &                      elem_dt(4), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_contiguous error')
    elem_dt(5) = mpi_integer
    elem_dt(6) = mpi_logical
    CALL create_sim_step_info_dt(elem_dt(7))
    CALL create_fname_metadata_dt(elem_dt(8))
    CALL mpi_get_address(dummy(1), base, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%name, displs(1), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%begin_str, displs(2), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%end_str, displs(3), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%intvl_str, displs(4), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%i_tag, displs(5), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%l_output_last, displs(6), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%sim_step_info, displs(7), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%fname_metadata, displs(8), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    displs = displs - base
    CALL mpi_type_create_struct(num_dt_elem, blocklens, displs, elem_dt, &
      &                         event_data_dt, ierror)
    IF (ierror /= mpi_success) &
         & CALL finish(routine, 'mpi_type_create_struct error')
    CALL mpi_type_get_extent(event_data_dt, stride, ext, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_get_extent error')
    CALL mpi_get_address(dummy(2), stride, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    stride = stride - base
    IF (stride < ext) THEN
      CALL finish(routine, 'fatal extent mismatch')
    ELSE IF (stride > ext) THEN
      CALL mpi_type_create_resized(event_data_dt, 0_mpi_address_kind, &
        &                          stride, resized_dt, ierror)
      IF (ierror /= mpi_success) &
        CALL finish(routine, 'mpi_type_create_resized error')
      CALL mpi_type_free(event_data_dt, ierror)
      IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_free error')
      event_data_dt = resized_dt
    END IF
    CALL mpi_type_commit(event_data_dt, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_commit error')
    DO i = 1, num_dt_elem
      IF (i /= 3 .AND. i < 5 .OR. i > 6) THEN
        CALL mpi_type_free(elem_dt(i), ierror)
        IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_free error')
      END IF
    END DO
    CALL memset_f(C_LOC(dummy(1)), 0, INT(2*stride, c_size_t))
    dummy(1)%name = 'test_name'
    dummy(1)%begin_str(1) = '1970-01-01'
    dummy(1)%begin_str(2:) = ''
    dummy(1)%end_str(1) = '1970-12-31'
    dummy(1)%end_str(2:) = ''
    dummy(1)%intvl_str(1) = '1month'
    dummy(1)%intvl_str(2:) = ''
    dummy(1)%i_tag = 700
    dummy(1)%l_output_last = .TRUE.
    temp_datetime => newDatetime('1969-01-01')
    dummy(1)%sim_step_info%sim_start = temp_datetime
    CALL deallocateDatetime(temp_datetime)
    temp_datetime => newDatetime('1975-12-06')
    dummy(1)%sim_step_info%sim_end = temp_datetime
    CALL deallocateDatetime(temp_datetime)
    temp_datetime => newDatetime('1970-01-01')
    dummy(1)%sim_step_info%run_start = temp_datetime
    CALL deallocateDatetime(temp_datetime)
    dummy(1)%sim_step_info%restart_time = dummy(1)%sim_step_info%sim_end
    dummy(1)%sim_step_info%dtime = 0.5_wp
    temp_datetime => newDatetime('1970-01-02')
    dummy(1)%sim_step_info%dom_start_time = temp_datetime
    CALL deallocateDatetime(temp_datetime)
    temp_datetime => newDatetime('1970-11-23')
    dummy(1)%sim_step_info%dom_end_time = temp_datetime
    CALL deallocateDatetime(temp_datetime)
    dummy(1)%sim_step_info%jstep0 = 1
    dummy(1)%fname_metadata%steps_per_file = 15
    dummy(1)%fname_metadata%steps_per_file_inclfirst = .FALSE.
    dummy(1)%fname_metadata%file_interval = '5days'
    dummy(1)%fname_metadata%phys_patch_id = -25
    dummy(1)%fname_metadata%ilev_type = 123456
    dummy(1)%fname_metadata%filename_format = 'dummy.nc'
    dummy(1)%fname_metadata%filename_pref = 'dummy-lalalal.nc'
    dummy(1)%fname_metadata%extn = 'nc'
    dummy(1)%fname_metadata%jfile_offset = 178
    dummy(1)%fname_metadata%npartitions = 2
    dummy(1)%fname_metadata%ifile_partition = 2
    CALL mpi_sendrecv(dummy(1), 1, event_data_dt, 0, 178, &
      &               dummy(2), 1, event_data_dt, 0, 178, &
      &               mpi_comm_self, mpi_status_ignore, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'transfer test error')
    IF (memcmp_f(C_LOC(dummy(1)), C_LOC(dummy(2)), INT(stride, c_size_t))) &
      CALL finish(routine, 'transfer test error')
  END SUBROUTINE create_event_data_dt

  SUBROUTINE create_sim_step_info_dt(dt)
    INTEGER, INTENT(out) :: dt
    INTEGER, PARAMETER :: num_dt_elem = 8
    INTEGER(mpi_address_kind) :: base, displs(num_dt_elem)
    INTEGER :: elem_dt(num_dt_elem), ierror
    CHARACTER(len=*), PARAMETER :: routine &
      = modname//"::create_sim_step_info_dt"
    INTEGER, PARAMETER :: blocklens(num_dt_elem) = 1
    TYPE(t_sim_step_info) :: dummy(1)
    CALL create_datetime_dt(elem_dt(1))
    elem_dt(2) = elem_dt(1)
    elem_dt(3) = elem_dt(1)
    elem_dt(4) = elem_dt(1)
    elem_dt(5) = elem_dt(1)
    elem_dt(6) = elem_dt(1)
    elem_dt(7) = p_real
    elem_dt(8) = mpi_integer
    CALL mpi_get_address(dummy(1), base, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%sim_start, displs(1), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%sim_end, displs(2), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%run_start, displs(3), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%restart_time, displs(4), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%dom_start_time, displs(5), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%dom_end_time, displs(6), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%dtime, displs(7), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%jstep0, displs(8), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    displs = displs - base
    CALL mpi_type_create_struct(num_dt_elem, blocklens, displs, elem_dt, &
      &                         dt, ierror)
    IF (ierror /= mpi_success) &
      & CALL finish(routine, 'mpi_type_create_struct error')
    CALL mpi_type_free(elem_dt(1), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_type_free error')
  END SUBROUTINE create_sim_step_info_dt

  SUBROUTINE create_datetime_dt(dt)
    INTEGER, INTENT(out) :: dt
    TYPE(datetime) :: dummy(1)
    INTEGER, PARAMETER :: num_dt_elem = 2
    INTEGER(mpi_address_kind) :: base, displs(num_dt_elem)
    INTEGER :: elem_dt(num_dt_elem), ierror
    CHARACTER(len=*), PARAMETER :: routine &
      = modname//"::create_datetime_dt"
    INTEGER, PARAMETER :: blocklens(num_dt_elem) = (/ 1, 6 /)
    elem_dt(1) = p_int_i8
    elem_dt(2) = p_int_i4
    CALL mpi_get_address(dummy(1), base, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%date%year, displs(1), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%date%month, displs(2), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    displs = displs - base
    CALL mpi_type_create_struct(num_dt_elem, blocklens, displs, elem_dt, &
      &                         dt, ierror)
    IF (ierror /= mpi_success) &
      & CALL finish(routine, 'mpi_type_create_struct error')
  END SUBROUTINE create_datetime_dt

  SUBROUTINE create_fname_metadata_dt(dt)
    INTEGER, INTENT(out) :: dt
    INTEGER, PARAMETER :: num_dt_elem = 11
    INTEGER(mpi_address_kind) :: base, displs(num_dt_elem)
    INTEGER :: elem_dt(num_dt_elem), ierror
    CHARACTER(len=*), PARAMETER :: routine &
      = modname//"::create_fname_metadata_dt"
    INTEGER, PARAMETER :: blocklens(num_dt_elem) &
      = (/ 1, 1, MAX_TIMEDELTA_STR_LEN, 1, 1, FILENAME_MAX, FILENAME_MAX, 16, &
      &    1, 1, 1 /)
    TYPE(t_fname_metadata) :: dummy(2)
    elem_dt(1) = mpi_integer
    elem_dt(2) = mpi_logical
    elem_dt(3) = mpi_character
    elem_dt(4) = mpi_integer
    elem_dt(5) = mpi_integer
    elem_dt(6) = mpi_character
    elem_dt(7) = mpi_character
    elem_dt(8) = mpi_character
    elem_dt(9) = mpi_integer
    elem_dt(10) = mpi_integer
    elem_dt(11) = mpi_integer
    CALL mpi_get_address(dummy(1), base, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%steps_per_file, displs(1), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%steps_per_file_inclfirst, displs(2), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%file_interval, displs(3), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%phys_patch_id, displs(4), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%ilev_type, displs(5), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%filename_format, displs(6), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%filename_pref, displs(7), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%extn, displs(8), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%jfile_offset, displs(9), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%npartitions, displs(10), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    CALL mpi_get_address(dummy(1)%ifile_partition, displs(11), ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'mpi_get_address error')
    displs = displs - base
    CALL mpi_type_create_struct(num_dt_elem, blocklens, displs, elem_dt, &
      &                         dt, ierror)
    IF (ierror /= mpi_success) &
      & CALL finish(routine, 'mpi_type_create_struct error')
  END SUBROUTINE create_fname_metadata_dt
#endif



  !---------------------------------------------------------------
  ! ROUTINES PERFORMING HANDSHAKE DURING EVENT STEPS
  !---------------------------------------------------------------

  !> Declare output step as being completed by this participating PE.
  !
  !  @note We do not use the wrapper routines from "mo_mpi" here,
  !        since we need direct control over the non-blocking P2P
  !        requests.
  !
  SUBROUTINE pass_output_step(event)
    TYPE(t_par_output_event), INTENT(inout) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::pass_output_step"
    LOGICAL :: step_is_not_finished
#ifndef NOMPI
    INTEGER :: ierrstat, mstat(MPI_STATUS_SIZE), istep, i_tag
#endif

    step_is_not_finished = .NOT. is_output_event_finished(event)

#ifndef NOMPI
    IF (use_async_name_list_io .AND. step_is_not_finished &
         .AND. (event%icomm /= MPI_COMM_NULL)) THEN
      ! wait for the last ISEND to be processed:
      IF (ldebug) WRITE (0,'(i0,a)') p_pe, ": waiting for request handle."
      CALL MPI_WAIT(event%isend_req, mstat, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_WAIT.')
      istep = event%output_event%i_event_step
      IF (ldebug .AND. istep > 1) THEN
        WRITE (0,'(i0,2(a,i0))') p_pe, ": completed request for rank ", &
             event%iroot, ', tag ', &
             event%output_event%event_step(istep-1)%event_step_data(1)%i_tag
      END IF
      ! launch a new non-blocking send:
      i_tag = event%output_event%event_step(istep)%event_step_data(1)%i_tag
      event%isend_buf = istep
      event%isend_buf = event%output_event%event_step(istep)%i_sim_step
      IF (ldebug) THEN
        WRITE (0,'(a,4(a,i0))') routine, ": ", p_pe, &
          &   ": sending message ", i_tag, " from ", get_my_global_mpi_id(), &
          &   " to ", event%iroot
      END IF
      CALL MPI_IBSEND(event%isend_buf, 1, p_int, event%iroot, i_tag, &
        &            event%icomm, event%isend_req, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_ISEND.')
      IF (ldebug) THEN
        WRITE (0,'(4(i0,a))') p_pe, ": pass ", event%output_event%i_event_step,&
          &         " (",  event%output_event%event_step(istep)%event_step_data(1)%jfile, &
          &         " / ", event%output_event%event_step(istep)%event_step_data(1)%jpart, &
          &         " )"
      END IF
    END IF
#endif

    IF (step_is_not_finished) THEN
      ! increment step counter
      event%output_event%i_event_step = event%output_event%i_event_step + 1
    END IF
  END SUBROUTINE pass_output_step


#ifndef NOMPI
  !> Place MPI_IRECV calls, thereby asking all participating PEs for
  !> acknowledging the completion of the output step.
  !
  !  @note We do not use the wrapper routines from "mo_mpi" here,
  !        since we need direct control over the non-blocking P2P
  !        requests.
  !
  SUBROUTINE trigger_output_step_irecv(event)
    TYPE(t_par_output_event), INTENT(inout) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::trigger_output_step_irecv"
    INTEGER                     :: ierrstat

    ! determine this PE's MPI rank wrt. the given MPI communicator and
    ! return if this PE is not root PE for the given event:
    IF (.NOT. is_event_root_pe(event))  RETURN

    ! wait for the last IRECVs to be processed:
    CALL wait_for_pending_irecvs(event)

    IF (ALLOCATED(event%irecv_req)) THEN
      DEALLOCATE(event%irecv_req, event%irecv_buf, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      event%irecv_nreq = 0
    END IF
    IF (.NOT. is_output_event_finished(event)) THEN
      CALL start_event_step_irecvs(event%icomm, &
           event%output_event%event_step(event%output_event%i_event_step), &
           event%irecv_req, event%irecv_buf, event%irecv_nreq)
    END IF
  END SUBROUTINE trigger_output_step_irecv
#endif

#ifndef NOMPI
  SUBROUTINE start_event_step_irecvs(icomm, event_step, &
       irecv_req, irecv_buf, irecv_nreq)
    INTEGER, INTENT(in) :: icomm
    TYPE(t_event_step),    INTENT(in) :: event_step
    INTEGER, ALLOCATABLE, INTENT(out) :: irecv_req(:), irecv_buf(:)
    INTEGER, INTENT(out)              :: irecv_nreq

    CHARACTER(LEN=*), PARAMETER :: routine = &
         modname//"::start_event_step_irecvs"
    INTEGER :: i, ierror, nreq, i_pe, i_tag
    ! launch a couple of new non-blocking receives:
    nreq = event_step%n_pes
    irecv_nreq = nreq

    ALLOCATE(irecv_req(nreq), irecv_buf(nreq), STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    irecv_req = MPI_REQUEST_NULL
    DO i=1, nreq
      i_pe  = event_step%event_step_data(i)%i_pe
      i_tag = event_step%event_step_data(i)%i_tag
      IF (ldebug) THEN
        WRITE (0,'(a,3(a,i0))') routine, ": launching IRECV ", i_tag, " on ", get_my_global_mpi_id(), &
          &         " from ", i_pe
      END IF
      CALL mpi_irecv(irecv_buf(i), 1, p_int, i_pe, &
        &            i_tag, icomm, irecv_req(i), ierror)
      IF (ierror /= mpi_success) CALL finish (routine, 'Error in MPI_IRECV.')
    END DO
  END SUBROUTINE start_event_step_irecvs


  !> Utility routine: blocking wait for pending "output completed"
  !> messages.
  !
  SUBROUTINE wait_for_pending_irecvs(event)
    TYPE(t_par_output_event), INTENT(inout) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::wait_for_pending_irecs"
    INTEGER                     :: ierrstat
    ! INTEGER, ALLOCATABLE        :: irecv_status(:,:)

    ! wait for the last IRECVs to be processed:
    IF (ALLOCATED(event%irecv_req)) THEN
      ! ALLOCATE(irecv_status(MPI_STATUS_SIZE,event%irecv_nreq), STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      CALL MPI_WAITALL(event%irecv_nreq, event%irecv_req, mpi_statuses_ignore, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish(routine, 'Error in MPI_WAITALL.')
      ! DEALLOCATE(irecv_status, STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      ! increment event step counter:
      event%output_event%i_event_step = event%output_event%i_event_step + 1
    END IF
  END SUBROUTINE wait_for_pending_irecvs
#endif


  !> Utility routine: blocking wait for pending "output completed"
  !> messages.
  !
  SUBROUTINE blocking_wait_for_irecvs(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
#ifndef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::blocking_wait_for_irecvs"
    INTEGER                           :: ierrstat, nreq, ireq
    INTEGER, ALLOCATABLE              :: irecv_req(:)
    TYPE(t_par_output_event), POINTER :: ev

    ! count the number of request handles
    ev   =>  event
    nreq = 0
    DO WHILE (ASSOCIATED(ev))
      IF (ALLOCATED(ev%irecv_req)) nreq = nreq + ev%irecv_nreq
      ev => ev%next
    END DO
    IF (ldebug) WRITE (0,*) "Total ", nreq, " IRECV request handles."
    IF (nreq > 0) THEN
      ! collect the request handles
      ALLOCATE(irecv_req(nreq), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      ev   =>  event
      ireq = 1
      DO WHILE (ASSOCIATED(ev))
        IF (ALLOCATED(ev%irecv_req)) &
          irecv_req(ireq:(ireq+ev%irecv_nreq-1)) = ev%irecv_req(1:ev%irecv_nreq)
        ireq = ireq + ev%irecv_nreq
        ev => ev%next
      END DO

      ! wait for the last IRECVs to be processed:
      CALL MPI_WAITALL(nreq, irecv_req, mpi_statuses_ignore, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_WAITALL.')
      DEALLOCATE(irecv_req, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

      ! clear the request handles
      ev => event
      DO WHILE (ASSOCIATED(ev))
        IF (ALLOCATED(ev%irecv_req)) &
          ev%irecv_req(:) = MPI_REQUEST_NULL
        ev => ev%next
      END DO
    END IF
#endif
  END SUBROUTINE blocking_wait_for_irecvs


  !> Spool event state fast-forward to a given event step.
  !
  !  This functionality is, e.g., required for resume after restart.
  !
  SUBROUTINE set_event_to_simstep(event, jstep, l_isrestart, lrecover_open_file)
    TYPE(t_output_event),   INTENT(INOUT), TARGET :: event              !< output event data structure
    INTEGER,                INTENT(IN)            :: jstep              !< simulation step
    LOGICAL,                INTENT(IN)            :: l_isrestart        !< .TRUE. if this is a restart run
    !> Flag. If true, we test for an existing file from previous runs
    LOGICAL,                INTENT(IN)            :: lrecover_open_file
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_event_to_simstep"
    CHARACTER(LEN=10) :: jpart_str
    INTEGER :: ev_step, istep, n_pes, i_pe, tlen
    TYPE(t_event_step_data), POINTER :: event_step_data

    ev_step = 0
    DO istep=1,event%n_event_steps
      IF (event%event_step(istep)%i_sim_step < jstep)  ev_step = istep
    END DO
    event%i_event_step = ev_step + 1

    IF (event%i_event_step <= event%n_event_steps) THEN
      istep = event%i_event_step
      n_pes = event%event_step(istep)%n_pes

      DO i_pe=1,n_pes
        event_step_data => event%event_step(istep)%event_step_data(i_pe)

        ! Spooling forward the event status to a given step means that
        ! we have to set the "open" flag.
        event_step_data%l_open_file = .TRUE.

        ! Test if the opening of the file at step "jstep" would
        ! destroy a file of the same name that has been created in a
        ! previous run (i.e. we are in a restart run).
        !
        ! This happens in the following two situations:
        !
        ! - Each file contains more than one output step and "jstep"
        !   corresponds to an output step "in the middle of the file".
        !
        ! - We are in a restart run and we are dealing with the first
        !   output file, which also contains the initial state.

        IF ((event_step_data%jpart > 1) .OR.         &
          & (l_isrestart .AND. (event_step_data%jfile == 1))) THEN
          ! Resuming after a restart means that we have to open the
          ! file for output though this has not been planned
          ! initially. We must find a unique suffix then for this new
          ! file (otherwise we would overwrite the file from the last
          ! step) and we must inform the user about this incident.
          IF (.NOT. lrecover_open_file) THEN
            ! simply throw an error message
            CALL finish(routine, "Attempt to overwrite existing file after restart!")
          ELSE
            ! otherwise: modify file name s.t. the new, resumed file
            ! is clearly distinguishable: We append "_<part>+"
            jpart_str = int2string(event_step_data%jpart)
            tlen = LEN_TRIM(event_step_data%filename_string)
            IF (my_process_is_mpi_workroot()) THEN
              WRITE (0,*) "Modify filename ", event_step_data%filename_string(1:tlen), " to ", &
                &      event_step_data%filename_string(1:tlen)//"_part_"//TRIM(jpart_str)//"+",  &
                &      " after restart."
            END IF
            CALL modify_filename(event, event_step_data%filename_string(1:tlen), &
              &       event_step_data%filename_string(1:tlen)//"_part_"//TRIM(jpart_str)//"+", &
              &       start_step=istep)
          END IF
        END IF

      END DO
    END IF
  END SUBROUTINE set_event_to_simstep


  !> Spool event state fast-forward to a given event step.
  !  Implementation for parallel events.
  !
  SUBROUTINE set_event_to_simstep_par(event, jstep, l_isrestart, lrecover_open_file)
    TYPE(t_par_output_event), TARGET    :: event              !< output event data structure
    INTEGER,                  INTENT(IN) :: jstep              !< simulation step
    LOGICAL,                  INTENT(IN) :: l_isrestart        !< .TRUE. if this is a restart run
    LOGICAL,                  INTENT(IN) :: lrecover_open_file !< Flag. If true, we test for an existing file from previous runs
    TYPE(t_par_output_event), POINTER :: ev
    ev => event
    DO WHILE (ASSOCIATED(ev))
      CALL set_event_to_simstep(ev%output_event, jstep, l_isrestart, lrecover_open_file)
      ev => ev%next
    END DO
  END SUBROUTINE set_event_to_simstep_par


  !> Utility routine: Loop over all event steps of a given event and
  !> modify all matching filenames.
  !
  SUBROUTINE modify_filename(event, old_name, new_name, start_step)
    TYPE(t_output_event), INTENT(INOUT), TARGET :: event              !< output event data structure    
    CHARACTER(LEN=*),     INTENT(IN)            :: old_name, new_name !< old file name string and replacement
    INTEGER,              INTENT(IN)            :: start_step         !< event step where to start searching
    ! local variables
    INTEGER :: istep, n_pes, i_pe, ipart
    TYPE(t_event_step_data), POINTER :: event_step_data

    ipart = 0
    DO istep=start_step,event%n_event_steps
      n_pes = event%event_step(istep)%n_pes
      DO i_pe=1,n_pes
        event_step_data => event%event_step(istep)%event_step_data(i_pe)        
        IF (event_step_data%filename_string == old_name) THEN
          event_step_data%filename_string = new_name
          ipart = ipart + 1
          event_step_data%jpart = ipart
        END IF
      END DO
    END DO
  END SUBROUTINE modify_filename

END MODULE mo_output_event_handler
