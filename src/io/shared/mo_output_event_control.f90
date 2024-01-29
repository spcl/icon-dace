! Routines that describe when regular output steps and ready file
! events shall be triggered. The definitions in this module control
! the behaviour of "mo_output_event_handler".
!
! See "mo_output_event_handler" for a detailed description.
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

MODULE mo_output_event_control

  USE mo_mpi,                ONLY: my_process_is_mpi_test, my_process_is_mpi_workroot
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,          ONLY: finish, message_text
  USE mo_kind,               ONLY: wp, i4, i8
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN,          &
    &                              MAX_TIMEDELTA_STR_LEN, datetime,                     &
    &                              deallocateDatetime, datetimeToString,                &
    &                              newDatetime, OPERATOR(>=), OPERATOR(*),              &
    &                              OPERATOR(+), OPERATOR(-), timedelta, newTimedelta,   &
    &                              deallocateTimedelta, OPERATOR(<=), OPERATOR(>),      &
    &                              OPERATOR(<), OPERATOR(==),                           &
    &                              divisionquotienttimespan,                            &
    &                              dividedatetimedifferenceinseconds,                   &
    &                              getPTStringFromMS, getPTStringFromSeconds,           &
    &                              timedeltaToString
  USE mo_var,                ONLY: lev_type_str
  USE mo_output_event_types, ONLY: t_sim_step_info, t_event_step_data
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords,    &
    &                              int2string, tolower
  USE mo_util_mtime,         ONLY: mtime_utils, FMT_DDHHMMSS, FMT_DDDHHMMSS, &
       mtime_timedelta_from_fseconds
  USE mo_name_list_output_types, ONLY: t_fname_metadata


  IMPLICIT NONE

  ! public subroutines + functions:
  PUBLIC :: compute_matching_sim_steps
  PUBLIC :: generate_output_filenames

  !---------------------------------------------------------------
  ! constants
  
  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output_event_control'

  !> Internal switch for debugging output
  LOGICAL,          PARAMETER :: ldebug  = .FALSE.

CONTAINS

  ! --------------------------------------------------------------------------------------------------
  !> Compute the matching simulation step for each date-time stamp.
  !
  !  The list of date-time stamp string is given in @p date_string. We
  !  assume that the date_string list is ordered in time. We choose
  !  the advection time step matching or following the time stamp
  !  string. The "true" date string is also returned.
  !
  !  @note The mtime library offers no arithmetic subroutines to
  !        compute dates. Therefore we must actually step through all
  !        dynamic timesteps.
  !
  ! --------------------------------------------------------------------------------------------------
  SUBROUTINE compute_matching_sim_steps(num_dates, dates, sim_step_info, &
    &                                   result_steps, result_exactdate)

    INTEGER,              INTENT(IN)    :: num_dates                             !< no. of string to convert
    !> array of mtime datetime objects
    TYPE(datetime),       INTENT(IN)    :: dates(:)
    TYPE(t_sim_step_info),INTENT(IN)    :: sim_step_info                        !< definitions: time step size, etc.
    INTEGER,              INTENT(OUT) :: result_steps(:)                      !< resulting step indices (+last sim step)
    !> resulting (exact) time step strings  (+last sim step)
    TYPE(datetime),        INTENT(OUT)   :: result_exactdate(:)

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_matching_sim_steps"
    INTEGER                  :: ilist
    TYPE(datetime) :: mtime_begin, mtime_end, mtime_dom_start, mtime_dom_end
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datetime_string

    ! Domains (and their output) can be activated and deactivated
    ! during the simulation. This is determined by the parameters
    ! "dom_start_time" and "dom_end_time". Therefore, we must create
    ! a corresponding event.
    mtime_dom_start = sim_step_info%dom_start_time
    mtime_dom_end   = sim_step_info%dom_end_time
    mtime_begin     = sim_step_info%run_start
    mtime_end       = sim_step_info%sim_end

    DO ilist = 1,num_dates
      ! check if domain is inactive:
      IF (((dates(ilist) >= mtime_dom_start) .AND. (dates(ilist) >= mtime_begin)) .OR.  &
        & ((dates(ilist) <= mtime_end) .AND. (dates(ilist) <= mtime_dom_end))) THEN 
        CALL compute_step(dates(ilist), mtime_begin, mtime_end,                &
          &               sim_step_info%dtime, sim_step_info%jstep0,           &
          &               result_steps(ilist), result_exactdate(ilist))
        IF (ldebug) THEN
          CALL datetimetostring(result_exactdate(ilist), datetime_string)
          WRITE (0,*) ilist, ": ", result_steps(ilist), " -> ", &
            &         TRIM(datetime_string)
        END IF
      ELSE
        result_steps(ilist)     = -1
      END IF
    END DO
    result_steps(num_dates+1:)     = -1

  END SUBROUTINE compute_matching_sim_steps


  ! --------------------------------------------------------------------------------------------------
  !> Utility function: Compute "(date1-sim_start)/(dtime)"
  ! --------------------------------------------------------------------------------------------------
  SUBROUTINE compute_step(mtime_current, mtime_begin, mtime_end, dtime,  &
    &                     step_offset, step, exact_date)
    TYPE(datetime),  INTENT(in)                      :: mtime_current       !< input date to translated into step
    TYPE(datetime),  INTENT(in)                      :: mtime_begin         !< begin of run (note: restart cases!)
    TYPE(datetime),  INTENT(in)                      :: mtime_end           !< end of run
    REAL(wp),                            INTENT(IN)  :: dtime               !< [s] length of a time step
    INTEGER,                             INTENT(IN)  :: step_offset
    INTEGER,                             INTENT(OUT) :: step                !< result: corresponding simulations step
    !> result: corresponding simulation date
    TYPE(datetime),                      INTENT(OUT) :: exact_date
    ! local variables
    REAL                                 :: intvlmillisec
    TYPE(datetime)                       :: mtime_step
    TYPE(divisionquotienttimespan)       :: tq
    TYPE(timedelta) :: vlsec

    ! first, we compute the dynamic time step which is equal or larger than
    ! the desired date "mtime_current"
    ! intvlsec    = REAL(dtime)
    ! step        = CEILING(datetimedividebyseconds(mtime_begin, mtime_date1, intvlsec))
    CALL mtime_timedelta_from_fseconds(dtime, mtime_begin, vlsec)

    CALL divideDatetimeDifferenceInSeconds(mtime_current, mtime_begin, vlsec, tq)

    step = INT(tq%quotient,i4)

    IF (step >= 0) exact_date = mtime_begin + step * vlsec

    ! then we add the offset "jstep0" (nonzero for restart cases):
    step        = step + step_offset

  END SUBROUTINE compute_step


  ! --------------------------------------------------------------------------------------------------
  !> Function for generating output file names.
  !
  !  This routine also prescribes, in which steps the output file is
  !  opened or closed.
  ! --------------------------------------------------------------------------------------------------
  SUBROUTINE generate_output_filenames(num_dates, dates, sim_steps, &
    &                                  sim_step_info, fname_metadata, &
    &                                  skipped_dates, result_fnames)
    !> no. of dates to convert
    INTEGER,                INTENT(IN)    :: num_dates
    TYPE(datetime), TARGET, INTENT(IN)    :: dates(:)     !< array of ISO 8601 time stamp strings
    INTEGER,                INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
    TYPE(t_sim_step_info),  INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
    TYPE(t_fname_metadata), INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
    INTEGER,                INTENT(IN)    :: skipped_dates
    TYPE(t_event_step_data), INTENT(out)  :: result_fnames(SIZE(dates))
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::generate_output_filenames"
    INTEGER                             :: i, j, ifile, ipart, total_index, this_jfile, errno
    CHARACTER(len=MAX_CHAR_LENGTH)      :: cfilename
    TYPE (t_keyword_list), POINTER      :: keywords
    TYPE(datetime) :: run_start, sim_start
    TYPE(datetime), POINTER :: step_date, mtime_begin, mtime_first, mtime_date
    TYPE(datetime) :: file_end
    TYPE(timedelta), POINTER            :: delta, forecast_delta
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string, forecast_delta_str

    run_start  = sim_step_info%run_start
    sim_start  = sim_step_info%sim_start

    ! ---------------------------------------------------
    ! prescribe, which file is used in which output step.
    ! ---------------------------------------------------
    ! a) user has specified "steps_per_file"
    ! b) user has specified "file_interval"
    IF (LEN_TRIM(fname_metadata%file_interval) == 0) THEN
      ! case a): steps_per_file
      mtime_first  => dates(1)
      ! special treatment for the initial time step written at the
      ! begin of the simulation: The first file gets one extra entry
      IF ( (run_start == mtime_first)              .AND.  &
        &  fname_metadata%steps_per_file_inclfirst   .AND.  &
        &  (fname_metadata%jfile_offset == 0) ) THEN
        result_fnames(1)%jfile = 1
        result_fnames(1)%jpart = 1
        DO i=2,num_dates
          result_fnames(i)%jfile = (i-2)/fname_metadata%steps_per_file + 1 ! INTEGER division!
          result_fnames(i)%jpart = i - 1 - (result_fnames(i)%jfile-1)*fname_metadata%steps_per_file
          IF (result_fnames(i)%jfile == 1)  result_fnames(i)%jpart = result_fnames(i)%jpart + 1
        END DO
      ELSE
        DO i=1,num_dates
          result_fnames(i)%jfile = (i-1)/fname_metadata%steps_per_file + 1 ! INTEGER division!
          result_fnames(i)%jpart = i - (result_fnames(i)%jfile-1)*fname_metadata%steps_per_file
        END DO
      END IF
    ELSE
      ! case b): file_interval
      ifile      =  1
      ipart      =  0
      delta     =>  newTimedelta(TRIM(fname_metadata%file_interval),errno)
      IF (errno /= SUCCESS) CALL finish(routine,"Wrong file_interval")
      file_end   =  dates(1) + delta

      OUTSTEP_LOOP : DO i=1,num_dates
        step_date => dates(i)
        IF (ldebug) THEN
          CALL datetimeToString(file_end, dtime_string)
          WRITE (0,'(a,i0,2a)') "i=", i, "file_end = ", dtime_string
        END IF
        IF (step_date >= file_end) THEN
          ifile = ifile + 1
          ipart = 0
          DO
            file_end =  file_end + delta
            IF (file_end > step_date) EXIT
          END DO
        END IF
        ipart = ipart + 1
        result_fnames(i)%jfile = ifile
        result_fnames(i)%jpart = ipart
      END DO OUTSTEP_LOOP
      CALL deallocateTimedelta(delta)
    END IF
    ! add offset to file numbers:
    IF (fname_metadata%jfile_offset /= 0) THEN
      result_fnames(1:num_dates)%jfile = result_fnames(1:num_dates)%jfile + fname_metadata%jfile_offset - 1
    END IF

    ! --------------------------------------------------------------
    ! prescribe, in which steps the output file is opened or closed:
    ! --------------------------------------------------------------

    IF (num_dates > 0) result_fnames(1)%l_open_file = .TRUE.
    DO i=2,num_dates
      result_fnames(i)%l_open_file &
        = result_fnames(i-1)%jfile /= result_fnames(i)%jfile
    END DO

    NULLIFY(keywords)
    ! ----------------------------------------------
    ! Set actual output file name (insert keywords):
    ! ----------------------------------------------
    forecast_delta => newTimedelta("P01D")
    DO i=1,num_dates
      IF (.NOT. result_fnames(i)%l_open_file) THEN
        ! if no file is opened: filename is identical to last step
        result_fnames(i)%filename_string = result_fnames(i-1)%filename_string
        CYCLE
      END IF
      ! otherwise, generate filename
      !
      ! define keywords:
      CALL associate_keyword("<path>",            TRIM(getModelBaseDir()),                                  keywords)
      CALL associate_keyword("<output_filename>", TRIM(fname_metadata%filename_pref),                       keywords)
      CALL associate_keyword("<physdom>",         TRIM(int2string(fname_metadata%phys_patch_id, "(i2.2)")), keywords)
      CALL associate_keyword("<levtype>",         TRIM(lev_type_str(fname_metadata%ilev_type)),             keywords)
      CALL associate_keyword("<levtype_l>",       TRIM(tolower(lev_type_str(fname_metadata%ilev_type))),    keywords)
      CALL associate_keyword("<npartitions>",     TRIM(int2string(fname_metadata%npartitions)),             keywords)
      CALL associate_keyword("<ifile_partition>", TRIM(int2string(fname_metadata%ifile_partition)),         keywords)
      CALL datetimeToString(dates(i), dtime_string)
      CALL associate_keyword("<datetime>",        TRIM(dtime_string),                                       keywords)
      ! keywords: compute current forecast time (delta):
      mtime_date => dates(i)
      CALL associate_keyword("<ddhhmmss>",                                                     &
        &                    TRIM(mtime_utils%ddhhmmss(sim_start, mtime_date, FMT_DDHHMMSS)),  &
        &                    keywords)
      CALL associate_keyword("<dddhhmmss>",                                                    &
        &                    TRIM(mtime_utils%ddhhmmss(sim_start, mtime_date, FMT_DDDHHMMSS)), &
        &                    keywords)

      forecast_delta     = mtime_date - sim_start
      forecast_delta_str = ""
      WRITE (forecast_delta_str,'(i3.3,2(i2.2))') forecast_delta%day*24 + forecast_delta%hour, &
        &                                         forecast_delta%minute, forecast_delta%second 
      CALL associate_keyword("<hhhmmss>",        TRIM(forecast_delta_str),                                 keywords)

      ! keywords: compose other variants of the absolute date-time
      !
      ! "YYYYMMDDThhmmssZ"     for the basic format of ISO8601 without the
      !                        separators "-" and ":" of the extended date-time format
      WRITE (dtime_string,'(i4.4,2(i2.2),a,3(i2.2),a)')                            &
        & mtime_date%date%year, mtime_date%date%month, mtime_date%date%day, 'T',   &
        & mtime_date%time%hour, mtime_date%time%minute, mtime_date%time%second, 'Z'
      CALL associate_keyword("<datetime2>",       TRIM(dtime_string),                                       keywords)

      ! "YYYYMMDDThhmmss.sssZ" for the basic format of ISO8601 with 3-digit 
      !                        fractions of seconds
      WRITE (dtime_string,'(i4.4,2(i2.2),a,3(i2.2),a,i3.3,a)')                                                 &
        &                      mtime_date%date%year, mtime_date%date%month, mtime_date%date%day, 'T',          &
        &                      mtime_date%time%hour, mtime_date%time%minute, mtime_date%time%second, '.',      &
        &                      mtime_date%time%ms, 'Z'
      CALL associate_keyword("<datetime3>",       TRIM(dtime_string),                                       keywords)

      ! "total_index": If we have split one namelist into concurrent,
      !                alternating files, compute the file index,
      !                which an "unsplit" namelist would have
      !                produced:
      IF (fname_metadata%npartitions > 1) THEN
        total_index = fname_metadata%npartitions*(result_fnames(i)%jfile+skipped_dates-1) + &
          &           fname_metadata%ifile_partition 
        this_jfile  = total_index
      ELSE
        total_index = result_fnames(i)%jfile
        this_jfile  = result_fnames(i)%jfile
      END IF
      CALL associate_keyword("<total_index>", TRIM(int2string(total_index,"(i4.4)")),  keywords)
      CALL associate_keyword("<jfile>",       TRIM(int2string(this_jfile, "(i4.4)")), keywords)

      cfilename = with_keywords(keywords, fname_metadata%filename_format)
      IF(my_process_is_mpi_test()) THEN
         ! (a) use filename with extension "_TEST" (then any post-processing needs to be adapted)
         WRITE(result_fnames(i)%filename_string,'(a,"_TEST",a)') TRIM(cfilename),TRIM(fname_metadata%extn)
!        ! (b) use standard filename
!        WRITE(result_fnames(i)%filename_string,'(a,a)')         TRIM(cfilename),TRIM(fname_metadata%extn)
      ELSE
        WRITE(result_fnames(i)%filename_string,'(2a)') TRIM(cfilename),TRIM(fname_metadata%extn)
      ENDIF

      ! consistency check: test if the user has accidentally chosen a
      ! file name syntax which yields duplicate file names:
      DO j=1,i-1
        IF (result_fnames(j)%l_open_file &
          & .AND.    result_fnames(i)%filename_string &
          &       == result_fnames(j)%filename_string) THEN
          ! found a duplicate filename:
          WRITE (message_text, '(3a)') &
               "Ambiguous output file name: '", &
               result_fnames(j)%filename_string, "'"
          CALL finish(routine, message_text)
        END IF
      END DO
    END DO
    CALL deallocateTimedelta(forecast_delta)


  END SUBROUTINE generate_output_filenames

END MODULE mo_output_event_control
