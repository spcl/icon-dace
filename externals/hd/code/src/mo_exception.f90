! mo_exception.f90 - Message routines for information, errors, warning and debugging
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_exception

!  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
!  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
!  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
!  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
!  doi: 10.1029/2018MS001400.

  USE mo_io_units, ONLY: nerr, nlog 
  USE mo_mpi,      ONLY: p_abort, p_parallel, p_parallel_io, p_pe 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, finish
  PUBLIC :: em_none, em_info, em_warn, em_error, em_param, em_debug
  PUBLIC :: open_log, close_log
  PUBLIC :: debug_messages_on, debug_messages_off
  PUBLIC :: number_of_warnings, number_of_errors

  INTEGER, PARAMETER :: em_none  = 0   ! normal message
  INTEGER, PARAMETER :: em_info  = 1   ! informational message
  INTEGER, PARAMETER :: em_warn  = 2   ! warning message: number of warnings counted
  INTEGER, PARAMETER :: em_error = 3   ! error message: number of errors counted
  INTEGER, PARAMETER :: em_param = 4   ! report parameter value
  INTEGER, PARAMETER :: em_debug = 5   ! debugging message

  CHARACTER(len=256) :: message_text = ''         !++mgs

  LOGICAL :: l_debug = .FALSE.
  LOGICAL :: l_log   = .FALSE.

  INTEGER :: number_of_warnings  = 0
  INTEGER :: number_of_errors    = 0

CONTAINS

  SUBROUTINE debug_messages_on
    l_debug = .TRUE.
  END SUBROUTINE debug_messages_on

  SUBROUTINE debug_messages_off
    l_debug = .FALSE.
  END SUBROUTINE debug_messages_off

  SUBROUTINE finish (name, text, exit_no)

#ifdef __INTEL_COMPILER
    USE ifcore
#endif

    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: text
    INTEGER,          INTENT(in), OPTIONAL :: exit_no

    INTEGER           :: iexit

#ifndef __STANDALONE
!OSBUTL    EXTERNAL :: util_exit 
#if ! (defined (__SX__) || defined (__INTEL_COMPILER) || defined (__xlC__))
    EXTERNAL :: util_backtrace
#endif
#endif

    WRITE (nerr,'(/,80("="),/)')
    IF (l_log) WRITE (nlog,'(/,80("="),/)')

    IF (PRESENT(exit_no)) THEN
       iexit = exit_no
    ELSE
       iexit = 1     ! POSIX defines this as EXIT_FAILURE
    END IF

    IF (PRESENT(text)) THEN
      IF (iexit == 1) THEN
        WRITE (nerr,'(a,a,a,a)') 'FATAL ERROR in ', TRIM(name), ': ', TRIM(text)
      ELSE
        WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
      ENDIF
      IF (l_log) WRITE (nlog,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
      IF (iexit == 1) THEN
        WRITE (nerr,'(a,a)') 'FATAL ERROR in ', TRIM(name)
      ELSE
        WRITE (nerr,'(1x,a)') TRIM(name)
      ENDIF
      IF (l_log) WRITE (nlog,'(a,a)') TRIM(name)
    ENDIF

    IF (p_parallel .AND. iexit == 1) THEN
      WRITE (nerr,'(1x,a,i0)') 'FINISH called from PE: ', p_pe
      IF (l_log) WRITE (nlog,'(1x,a,i0)') 'FINISH called from PE: ', p_pe
    ENDIF

    ! generate traceback in error cases
    IF (iexit == 1) THEN
      WRITE (nerr,'(/,80("-"),/,/)')
      IF (l_log) WRITE (nlog,'(/,80("-"),/,/)')
#ifdef __INTEL_COMPILER
      CALL tracebackqq(STRING="Error tracing with tracebackqq:",USER_EXIT_CODE=1)

#else
#ifdef __xlC__
      CALL xl__trbk
#else
#ifndef __STANDALONE
#if ! defined (__SX__)
      CALL util_backtrace
#endif
#endif
#endif
#endif
    ENDIF

    WRITE (nerr,'(/,80("="),/)')
    IF (l_log) WRITE (nlog,'(/,80("="),/)')
    
    IF (p_parallel) THEN 
      CALL p_abort
    ELSE
#ifndef __STANDALONE
!OSBUTL       CALL util_exit(iexit)
       STOP 'util_exit'
#else
       STOP 'mo_exception: finish ..'
#endif
    END IF

  END SUBROUTINE finish

  SUBROUTINE message (name, text, out, level, all_print, adjust_right)

    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text
    INTEGER,           INTENT(in), OPTIONAL :: out
    INTEGER,           INTENT(in), OPTIONAL :: level
    LOGICAL,           INTENT(in), OPTIONAL :: all_print
    LOGICAL,           INTENT(in), OPTIONAL :: adjust_right

    INTEGER :: iout
    INTEGER :: ilevel
    LOGICAL :: lprint
    LOGICAL :: ladjustr     !++mgs renamed from ladjust to ladjustr

    CHARACTER(len=32) :: prefix

    CHARACTER(len=LEN(message_text)) :: write_text

    IF (PRESENT(all_print)) THEN
      lprint = all_print
    ELSE
      lprint = .FALSE.
    ENDIF

    IF (PRESENT(adjust_right)) THEN
      ladjustr = adjust_right 
    ELSE
      ladjustr = .FALSE.
    ENDIF

    IF (PRESENT(out)) THEN
      iout = out
    ELSE
      iout = nerr
    END IF

    IF (PRESENT(level)) THEN
      ilevel = level
    ELSE
      ilevel = em_none
    END IF

    SELECT CASE (ilevel)
    CASE (em_none)  ; prefix = '        '
    CASE (em_info)  ; prefix = 'INFO   :'
    CASE (em_warn)  ; prefix = 'WARNING:' ; number_of_warnings  = number_of_warnings+1
    CASE (em_error) ; prefix = 'ERROR  :' ; number_of_errors    = number_of_errors+1
    CASE (em_param) ; prefix = '---     '
    CASE (em_debug) ; prefix = 'DEBUG  :'
    END SELECT

    IF (.NOT. ladjustr) THEN
      message_text = TRIM(ADJUSTL(text))
    ENDIF
    IF (name /= '')  THEN
      message_text = TRIM(name) // ': ' // TRIM(message_text)
    ENDIF
    IF (ilevel > em_none) THEN
      message_text = TRIM(prefix) // ' ' // TRIM(message_text)
    ENDIF

    IF (p_parallel .AND. (l_debug .OR. ilevel == em_error)) THEN
     WRITE(write_text,'(1x,a,i6,a,a)') 'PE ', p_pe, ' ', TRIM(message_text)
     lprint = .TRUE.
   ELSE
     write_text = message_text
   END IF

   IF (p_parallel_io .OR. lprint) THEN
     WRITE(iout,'(1x,a)') TRIM(write_text)
     IF (l_log) WRITE(nlog,'(1x,a)') TRIM(write_text)
   END IF
     
  END SUBROUTINE message

  SUBROUTINE open_log (logfile_name)

    CHARACTER(len=*), INTENT(in) :: logfile_name
    LOGICAL                      :: l_opened 

    INQUIRE (UNIT=nlog,OPENED=l_opened)

    IF (l_opened) THEN
      WRITE (message_text,'(a)') 'log file unit has been used already.'
      CALL message ('open_log', message_text)
      WRITE (message_text,'(a)') 'Close unit and reopen for log file.'
      CALL message ('open_log', message_text, level=em_warn)
      CLOSE (nlog)
    ENDIF

    OPEN (nlog,file=TRIM(logfile_name))
  
    l_log = .TRUE.

  END SUBROUTINE open_log

  SUBROUTINE close_log
    LOGICAL :: l_opened
   
    INQUIRE (UNIT=nlog,OPENED=l_opened)
    IF (l_opened) THEN
      CLOSE (nlog)
    ENDIF

    l_log = .FALSE.

  END SUBROUTINE close_log

END MODULE mo_exception
