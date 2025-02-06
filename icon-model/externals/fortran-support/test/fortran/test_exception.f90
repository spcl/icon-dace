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

MODULE test_mo_exception
  USE FORTUTF
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64, &
                                                                              i8 => int64
  USE mo_exception
  USE mo_io_units, ONLY: find_next_free_unit
  USE helpers, ONLY: custom_exit, open_logfile, open_new_logfile, custom_exit_dummy

CONTAINS

  SUBROUTINE TEST_message
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CALL TAG_TEST("TEST_message")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(0, .TRUE., nerr, callback_abort=custom_exit)

    CALL message('Test', 'Message')
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    READ (nerr, '(A)') log_in_file

    CALL STRING_CONTAINS(log_in_file, 'Test: Message')
    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_message_own_unit
    INTEGER :: nerr_1, nerr_2
    CHARACTER(len=100) :: log_in_file, logfile_1, logfile_2
    logfile_1 = 'logger_output_primary.txt'
    logfile_2 = 'logger_output_secondary.txt'

    CALL open_new_logfile(nerr_1, TRIM(logfile_1))

    CALL init_logger(0, .TRUE., nerr_1, callback_abort=custom_exit)

    CALL message('Test', 'Message')
    CLOSE (nerr_1)

    CALL open_new_logfile(nerr_2, TRIM(logfile_2))
    CALL message_to_own_unit('Ocean', 'Component', nerr_2)
    CLOSE (nerr_2)

    ! check log in primary-file
    CALL TAG_TEST("TEST_message")
    CALL open_logfile(nerr_1, TRIM(logfile_1))
    READ (nerr_1, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'Test: Message')
    CLOSE (nerr_1)

    ! check log in secondary_file
    CALL TAG_TEST("TEST_message_own_unit")
    CALL open_logfile(nerr_2, TRIM(logfile_2))
    READ (nerr_1, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'Ocean: Component')
    CLOSE (nerr_2)

  END SUBROUTINE

  SUBROUTINE TEST_message_all_print
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CALL TAG_TEST("TEST_all_print")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(0, .FALSE., nerr, callback_abort=custom_exit)

    ! not written to logfile
    CALL message('Test', 'Message')
    CALL message('Test', 'ALL PRINT', all_print=.TRUE.)
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    READ (nerr, '(A)') log_in_file

    CALL STRING_CONTAINS(log_in_file, 'Test: ALL PRINT')
    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_debug
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CALL TAG_TEST("TEST_debug")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(0, .FALSE., nerr, callback_abort=custom_exit)

    ! is not written to logfile
    CALL message('Info', 'Message')

    CALL debug_on()
    ! is written to logfile
    CALL message('Debug', 'Message')
    CALL debug_off()

    ! is not written to logfile
    CALL message('Info', 'Message')
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    READ (nerr, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'DEBUG PE:     0 Debug: Message')
    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_warning
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CALL TAG_TEST("TEST_warning")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(3, .FALSE., nerr, callback_abort=custom_exit)

    ! is not written to logfile
    CALL message('Info', 'Message')

    ! is written to logfile
    CALL warning('Warning', 'Message')

    ! is not written to logfile
    CALL message('Info', 'Message')
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    READ (nerr, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'WARNING PE:     3 Warning: Message')
    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_param
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(3, .TRUE., nerr, callback_abort=custom_exit)

    ! print one parameter per data-type
    CALL print_value('-', .TRUE.)
    CALL print_value('-', 6)
    CALL print_value('-', 6_i8)
    CALL print_value('-', 6.0_wp)

    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))

    READ (nerr, '(A)') log_in_file
    CALL TAG_TEST("TEST_param_bool")
    CALL STRING_CONTAINS(log_in_file, '- : TRUE')

    READ (nerr, '(A)') log_in_file
    CALL TAG_TEST("TEST_param_integer_i4")
    CALL STRING_CONTAINS(log_in_file, '- :         6')

    READ (nerr, '(A)') log_in_file
    CALL TAG_TEST("TEST_param_integer_i8")
    CALL STRING_CONTAINS(log_in_file, '- :         6')

    READ (nerr, '(A)') log_in_file
    CALL TAG_TEST("TEST_param_real")
    CALL STRING_CONTAINS(log_in_file, '- :  6.0000')

    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_error
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CALL TAG_TEST("TEST_error")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(5, .FALSE., nerr, callback_abort=custom_exit_dummy)

    CALL finish('Testmodule', 'Problem')
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    READ (nerr, '(A)') log_in_file

    CALL STRING_CONTAINS(log_in_file, 'FINISH PE:     5 Testmodule: Problem')
    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_timestamp
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CHARACTER(len=10) :: ctime
    CHARACTER(len=8)  :: cdate
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(0, .TRUE., nerr, callback_abort=custom_exit)
    CALL set_msg_timestamp(.TRUE.)
    CALL message('Test', 'Message')
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    CALL set_msg_timestamp(.FALSE.)
    READ (nerr, '(A)') log_in_file

    CALL DATE_AND_TIME(date=cdate, time=ctime)

    CALL TAG_TEST("TEST_timestamp_date")
    CALL STRING_CONTAINS(log_in_file, cdate)

    CALL TAG_TEST("TEST_timestamp_time")
    ! check only up to minutes
    CALL STRING_CONTAINS(log_in_file, ctime(1:4))

    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_extra_output
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CALL TAG_TEST("TEST_extra_output")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(234, .FALSE., nerr, callback_abort=custom_exit, l_extra_output=.TRUE., &
                     extra_info_prefix='Extra Prefix for Test')

    CALL message('Test', 'Message')
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    READ (nerr, '(A)') log_in_file

    CALL STRING_CONTAINS(log_in_file, 'Extra Prefix for Test PE:   234 Test: Message')
    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE TEST_enable_logging
    INTEGER :: nerr
    CHARACTER(len=100) :: log_in_file, logfile
    CALL TAG_TEST("TEST_enable_logging")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(0, .TRUE., nerr, callback_abort=custom_exit_dummy)

    CALL enable_logging(.FALSE.)
    CALL message('Test', 'No write to logfile')
    CALL enable_logging(.TRUE.)
    CALL message('Test', 'Should be in logfile')
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    ! first line of file is message from abort
    READ (nerr, '(A)') log_in_file

    CALL STRING_CONTAINS(log_in_file, 'Should be in logfile')
    CLOSE (nerr)
  END SUBROUTINE

END MODULE test_mo_exception
