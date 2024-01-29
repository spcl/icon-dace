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

MODULE test_mo_util_table
  USE FORTUTF
  USE mo_util_table, ONLY: t_table, initialize_table, add_table_column, set_table_entry, print_table
  USE helpers, ONLY: open_logfile, open_new_logfile

CONTAINS

  SUBROUTINE TEST_table()
    TYPE(t_table) :: table
    INTEGER :: nerr
    CHARACTER(LEN=100) :: log_in_file, logfile

    logfile = 'table_dump.txt'
    CALL open_new_logfile(nerr, TRIM(logfile))

    ! Initialize the table
    CALL initialize_table(table)

    ! Add columns to the table
    CALL add_table_column(table, "Name")
    CALL add_table_column(table, "Age")
    CALL add_table_column(table, "City")

    ! Set entries in the table
    CALL set_table_entry(table, 1, "Name", "John Doe")
    CALL set_table_entry(table, 1, "Age", "25")
    CALL set_table_entry(table, 1, "City", "New York")

    CALL set_table_entry(table, 2, "Name", "Jane Smith")
    CALL set_table_entry(table, 2, "Age", "30")
    CALL set_table_entry(table, 2, "City", "Los Angeles")

    CALL set_table_entry(table, 3, "Name", "Mike Johnson")
    CALL set_table_entry(table, 3, "Age", "40")
    CALL set_table_entry(table, 3, "City", "Chicago")

    ! Print the table to file
    CALL print_table(table, opt_dstfile=nerr)
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))

    ! first line is empty
    READ (nerr, '(A)') log_in_file

    CALL TAG_TEST("TEST_table_header")
    READ (nerr, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'Name         | Age | City        |')

    CALL TAG_TEST("TEST_table_delimiter")
    READ (nerr, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, '------------ | --- | ----------- |')

    ! skip empty line
    READ (nerr, '(A)') log_in_file

    CALL TAG_TEST("TEST_table_1st_entry")
    READ (nerr, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'John Doe     | 25  | New York    |')

    CALL TAG_TEST("TEST_table_2nd_entry")
    READ (nerr, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'Jane Smith   | 30  | Los Angeles |')

    CALL TAG_TEST("TEST_table_3rd_entry")
    READ (nerr, '(A)') log_in_file
    CALL STRING_CONTAINS(log_in_file, 'Mike Johnson | 40  | Chicago     |')

    CLOSE (nerr)
  END SUBROUTINE

END MODULE test_mo_util_table
