!> Contains time control for JSBACH
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_time
#ifndef __NO_JSBACH__

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_jsb_time_iface, ONLY: t_datetime, deallocateDatetime, get_time_next, &
                               get_date_components, get_time_start, &
                               get_time_dt, is_time_experiment_start, get_time_experiment_start, is_time_restart, &
                               get_year_length, get_month_length, get_day_length, get_year_day,   &
                               is_time_ltrig_rad_m1, get_time_interpolation_weights, get_asselin_coef

  IMPLICIT NONE
  PRIVATE

!!$  PUBLIC :: init_jsb_time
  PUBLIC :: init_time
  PUBLIC :: t_datetime, deallocateDatetime, get_time_next
  PUBLIC :: get_date_components, get_time_start, get_time_experiment_start, get_time_dt, &
            is_time_experiment_start, is_time_restart, is_time_ltrig_rad_m1, &
            get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: is_newday, is_newmonth, is_newyear, get_year_at_model_start, get_year_at_experiment_start, &
            get_year, get_month, get_day, get_secs_of_day
  PUBLIC :: get_time_interpolation_weights
  PUBLIC :: get_asselin_coef
  PUBLIC :: timestep_in_days, timesteps_per_day

  TYPE(t_datetime), POINTER, SAVE :: start_date
  TYPE(t_datetime), POINTER, SAVE :: final_date

  LOGICAL, SAVE :: time_initialized = .FALSE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time'

CONTAINS

  SUBROUTINE init_time(namelist_filename, model_shortname, restart_jsbach)

    USE mo_jsb_time_iface, ONLY: read_time_namelist, & !, configure_time
                                 configure_time_and_events, &
                                 get_time_start, get_time_stop, get_year_day

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: namelist_filename
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname
    LOGICAL, OPTIONAL, INTENT(IN) :: restart_jsbach

    REAL(wp) :: eps

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_time'

    !TODO? currently the time is initalised from the namelist of the first processed model
    !      but then used globally for all models!
    !      If different models should have different time information this would need to change
    !      -> Attention: it is not possible to specify this global information before knowing the names
    !         of the models: knowing the models is a precondition for restart with echam, because echam restart files
    !        contain the short name of a model. The time-information of a restart file in turn needs
    !        to be processed before times and events can be initialised.
    IF (.NOT. time_initialized) THEN
      CALL message(TRIM(routine), 'Initializing JSBACH time')

      ! TODO: remove #ifdef
#ifdef __ICON__
      IF (PRESENT(restart_jsbach)) CONTINUE ! avoid compiler warnings about dummy arguments not being used
      IF (PRESENT(namelist_filename)) CONTINUE ! avoid compiler warnings about dummy arguments not being used
      ! IF (PRESENT(namelist_filename)) THEN                   ! present only for standalone model
      !   CALL read_time_namelist(TRIM(namelist_filename))     ! already done from mo_util_jsbach.f90 in icon
      ! END IF
#else
      IF (PRESENT(namelist_filename) .AND. PRESENT(restart_jsbach)) THEN  ! present only for standalone model
        CALL read_time_namelist(TRIM(namelist_filename), restart_jsbach)
      END IF
#endif

      IF (PRESENT(model_shortname)) THEN
        CALL configure_time_and_events(model_shortname)
      ELSE
        CALL configure_time_and_events()
      ENDIF

      start_date => get_time_start()
      final_date => get_time_stop()

      ! just to avoid compiler warnings about variable set but never referenced
      IF (.FALSE.) THEN
        PRINT*, get_year_day(start_date)
        PRINT*, get_year_day(final_date)
      END IF

      ! Check if Asselin time filter is used and print message
      eps = get_asselin_coef()
      IF (eps > 0._wp) THEN
        WRITE(message_text,*) 'Using Asselin time filter for surface temperature: eps =', eps
      ELSE
        WRITE(message_text,*) 'Not using Asselin time filter for surface temperature'
      END IF
      CALL message(TRIM(routine), TRIM(message_text))

      time_initialized = .TRUE.
    ENDIF

  END SUBROUTINE init_time

  LOGICAL FUNCTION is_newday(current, dt)

    USE mo_jsb_time_iface, ONLY: get_time_previous

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_day, previous_day

    previous => get_time_previous(current, dt)

    ! CALL get_date_components(current, day=current_day)
    ! CALL get_date_components(previous, day=previous_day)
    current_day = INT(get_year_day(current))
    previous_day = INT(get_year_day(previous))

    is_newday = current_day /= previous_day

    CALL deallocateDatetime(previous)

  END FUNCTION is_newday

  LOGICAL FUNCTION is_newmonth(current, dt)

    USE mo_jsb_time_iface, ONLY: get_time_previous

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_month, previous_month

    previous => get_time_previous(current, dt)

    current_month  = INT(get_month(current))
    previous_month = INT(get_month(previous))

    is_newmonth = current_month /= previous_month

    CALL deallocateDatetime(previous)

  END FUNCTION is_newmonth

  LOGICAL FUNCTION is_newyear(current, dt)

    USE mo_jsb_time_iface, ONLY: get_time_previous

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_year, previous_year

    previous => get_time_previous(current, dt)

    ! CALL get_date_components(current, year=current_year)
    ! CALL get_date_components(previous, year=previous_year)
    current_year  = get_year(current)
    previous_year = get_year(previous)

    is_newyear = current_year /= previous_year

    CALL deallocateDatetime(previous)

  END FUNCTION is_newyear

  INTEGER FUNCTION get_year(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: year

    CALL get_date_components(this_datetime, year=year)

    get_year = year

  END FUNCTION get_year

  INTEGER FUNCTION get_month(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: month

    CALL get_date_components(this_datetime, month=month)

    get_month = month

  END FUNCTION get_month

  INTEGER FUNCTION get_day(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: day

    CALL get_date_components(this_datetime, day=day)

    get_day = day

  END FUNCTION get_day

  INTEGER FUNCTION get_secs_of_day(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: hour, minute, second

    CALL get_date_components(this_datetime, hour=hour, minute=minute, second=second)

    get_secs_of_day = NINT(3600._wp * REAL(hour,wp) + 60._wp * REAL(minute,wp) + REAL(second,wp))

  END FUNCTION get_secs_of_day

  INTEGER FUNCTION get_year_at_model_start()

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_year_at_model_start'

    IF (.NOT. time_initialized) &
      CALL finish(TRIM(routine), 'Time not initialized yet.')

    get_year_at_model_start = get_year(start_date)

  END FUNCTION get_year_at_model_start

  INTEGER FUNCTION get_year_at_experiment_start()

    TYPE(t_datetime), POINTER :: temp_date

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_year_at_experiment_start'

    IF (.NOT. time_initialized) &
      CALL finish(TRIM(routine), 'Time not initialized yet.')

    temp_date => get_time_experiment_start()
    get_year_at_experiment_start = get_year(temp_date)
    CALL deallocateDatetime(temp_date)

  END FUNCTION get_year_at_experiment_start

  REAL(wp) FUNCTION timestep_in_days(dt, model_id) ! EXPRESS the time step in days (as fraction of day)

    REAL(wp), OPTIONAL, INTENT(in) :: dt
    INTEGER,  OPTIONAL, INTENT(in) :: model_id

    CHARACTER(len=*), PARAMETER :: routine = modname//':timestep_in_days'

    IF (PRESENT(dt)) THEN
      timestep_in_days = dt / 86400._wp
    ELSE IF (PRESENT(model_id)) THEN
      timestep_in_days = get_time_dt(model_id) / 86400._wp
    ELSE
      CALL finish(routine, 'Must provide either dt oder model_id')
    END IF

  END FUNCTION timestep_in_days

  INTEGER FUNCTION timesteps_per_day(dt, model_id) ! NUMBER of time steps per day

    REAL(wp), OPTIONAL, INTENT(in) :: dt
    INTEGER,  OPTIONAL, INTENT(in) :: model_id

    CHARACTER(len=*), PARAMETER :: routine = modname//':timesteps_per_day'

    IF (PRESENT(dt)) THEN
      IF( mod(86400,INT(dt)) /= 0) THEN ! For computing the mean day temperature it is assumed that each day
                                        ! is computed with a fixed number of time steps. Therefore the program is
                                        ! stopped wheen day length is not an integer multiple of the time step.
         CALL finish("timesteps_per_day", 'Day length is not an integer multiple of the time step!')
      ELSE
         timesteps_per_day = 86400/INT(dt)
      END IF
    ELSE IF (PRESENT(model_id)) THEN
      timesteps_per_day = 86400 / INT(get_time_dt(model_id))
    ELSE
      CALL finish(routine, 'Must provide either dt oder model_id')
    END IF

  END FUNCTION timesteps_per_day

#endif
END MODULE mo_jsb_time
