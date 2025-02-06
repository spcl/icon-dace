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

! Contains utility routines for a TRUE/FALSE table.
!
! The table has the following layout: Each table row is determined by
! a row name. Then, in each column the table entry is marked by an
! "x" (or a user-defined symbol) where the column contains this row
! name.
!
! Code Example:
!
!   CALL init_bool_table(table)
!   CALL add_column(table, "cookies", (/ 'round ', 'edible' /)          )
!   CALL add_column(table, "apples",  (/ 'round ', 'juicy ', 'edible' /))
!   CALL add_column(table, "spiders", (/ 'juicy ', 'edible' /)          )
!   CALL print_bool_table(table)
!
! gives the following output:
!
!         | cookies | apples | spiders |
!
!  round  |  x      |  x     |         |
!  edible |  x      |  x     |  x      |
!  juicy  |         |  x     |  x      |
!
! @todo Implement a destruction routine.

MODULE mo_util_bool_table

  USE mo_impl_constants,      ONLY : SUCCESS
  USE mo_var_metadata_types,  ONLY : VARNAME_LEN
  USE mo_exception,           ONLY : finish
  USE mo_util_string,         ONLY : remove_duplicates, tolower, int2string


  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: init_bool_table
  PUBLIC  :: add_column
  PUBLIC  :: print_bool_table
  PUBLIC  :: t_column
  PUBLIC  :: t_bool_table

  ! ---------------------------------------------------------------------
  ! CONSTANTS
  ! ---------------------------------------------------------------------

  INTEGER, PARAMETER :: MAX_TABLE_COLUMNS =  100          !< max num. of table columns
  INTEGER, PARAMETER :: MAX_TABLE_ROWS    = 1000          !< max num. of table rows
  INTEGER, PARAMETER :: MAX_TITLE_LEN     =   64          !< max length of column title
  INTEGER, PARAMETER :: MAX_COLUMN_WIDTH  =   16          !< max width of column entry
  INTEGER, PARAMETER :: MAX_ROWNAME_LEN   = VARNAME_LEN   !< max length of row name

  CHARACTER(LEN=*), PARAMETER :: DELIMITER             = ' | '    !< vertical line
  CHARACTER(LEN=*), PARAMETER :: DEFAULT_CHAR_TRUE     = ' x '    !< TRUE  sign (default)
  CHARACTER(LEN=*), PARAMETER :: DEFAULT_CHAR_FALSE    = '   '    !< FALSE sign (default)

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_util_bool_table'

  ! ---------------------------------------------------------------------
  ! TYPE DEFINITIONS
  ! ---------------------------------------------------------------------

  !> Type definition for a single table column
  TYPE t_column
    CHARACTER(LEN=MAX_TITLE_LEN) :: title
    INTEGER                      :: width
    CHARACTER(LEN=MAX_ROWNAME_LEN),   ALLOCATABLE :: true_entries(:)
    CHARACTER(LEN=MAX_COLUMN_WIDTH),  ALLOCATABLE :: markers(:)
  END TYPE t_column

  !> Type definition for a complete table
  TYPE t_bool_table
    INTEGER                         :: n_columns
    TYPE (t_column)                 :: column(MAX_TABLE_COLUMNS) 
    CHARACTER(LEN=MAX_ROWNAME_LEN)  :: rowname(MAX_TABLE_ROWS) 
    INTEGER                         :: n_rows
    INTEGER                         :: rowname_len                !< max. length of a row name
  END TYPE t_bool_table


CONTAINS

  !> Initialize the table data structure.
  !
  SUBROUTINE init_bool_table(table)
    TYPE (t_bool_table), INTENT(INOUT) :: table

    table%n_columns   = 0
    table%n_rows      = 0
    table%rowname_len = 0
  END SUBROUTINE init_bool_table


  !> Add a column to the table.
  !
  SUBROUTINE add_column(table, colname, str_list, opt_nitems, opt_markers)
    TYPE (t_bool_table),   INTENT(INOUT)      :: table              !< true/false table object
    CHARACTER(LEN=*),      INTENT(IN)         :: colname            !< column name
    CHARACTER(LEN=*),      INTENT(IN)         :: str_list(:)        !< list of "true" entries
    INTEGER,          INTENT(IN), OPTIONAL    :: opt_nitems         !< optional: length of given list
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL    :: opt_markers(:)     !< list of entry markers (default: 'x')
    ! local variables    
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::add_column'
    INTEGER :: ierrstat, i, rowname_len, nitems

    nitems = SIZE(str_list)
    IF (PRESENT(opt_nitems))  nitems = opt_nitems

    table%n_columns = table%n_columns + 1
    ALLOCATE(table%column(table%n_columns)%true_entries(nitems), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! set column title string
    IF (LEN_TRIM(colname) > MAX_TITLE_LEN) THEN
      CALL finish(routine, "title name length too long!")
    ELSE
      table%column(table%n_columns)%title = TRIM(colname)
      table%column(table%n_columns)%width = LEN_TRIM(colname)
    END IF
    ! fill list of row entries
    DO i=1,nitems
      rowname_len = LEN_TRIM(str_list(i))
      IF (rowname_len > MAX_ROWNAME_LEN) THEN
        CALL finish(routine, "row name length too long!")
      ELSE
        table%rowname_len = MAX(table%rowname_len, rowname_len)
        table%column(table%n_columns)%true_entries(i) = TRIM(tolower(TRIM(str_list(i))))

        ! add entry also to the list of rows (duplicates will be removed later)
        table%n_rows = table%n_rows + 1
        table%rowname(table%n_rows) = TRIM(tolower(TRIM(str_list(i))))
      END IF
    END DO
    ! throw out duplicates:
    CALL remove_duplicates(table%rowname, table%n_rows)
    ! Optional: store entry markers
    IF (PRESENT(opt_markers)) THEN
      IF (SIZE(opt_markers) < nitems) &
        & CALL finish(routine, "Provided marker list is too short!")
      ALLOCATE(table%column(table%n_columns)%markers(nitems), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      DO i=1,nitems
        IF (LEN_TRIM(opt_markers(i)) > MAX_COLUMN_WIDTH) THEN
          CALL finish(routine, "Marker entry is too long!")
        ELSE
          table%column(table%n_columns)%markers(i) = TRIM(opt_markers(i))
          table%column(table%n_columns)%width = MAX(table%column(table%n_columns)%width, LEN_TRIM(opt_markers(i)))
        END IF
      END DO
    END IF
  END SUBROUTINE add_column


  !> @return string corresponding to TRUE, if @p rowname is set in
  !          column @p column, otherwise return the string
  !          corresponding to FALSE.
  !
  FUNCTION get_bool_table_entry(rowname, column)
    CHARACTER(LEN=MAX_COLUMN_WIDTH) :: get_bool_table_entry
    CHARACTER(len=*), INTENT(IN) :: rowname
    TYPE (t_column),  INTENT(IN) :: column
    ! local variables
    INTEGER :: i
    
    get_bool_table_entry = DEFAULT_CHAR_FALSE
    FIND_LOOP : DO i=1,SIZE(column%true_entries)
      IF (TRIM(tolower(TRIM(rowname))) == TRIM(tolower(TRIM(column%true_entries(i))))) THEN
        IF (ALLOCATED(column%markers)) THEN
          get_bool_table_entry = column%markers(i)
        ELSE
          get_bool_table_entry = DEFAULT_CHAR_TRUE
        END IF
        EXIT FIND_LOOP
      END IF
    END DO FIND_LOOP
  END FUNCTION get_bool_table_entry


  !> Print the TRUE/FALSE table
  !
  SUBROUTINE print_bool_table(table)
    TYPE (t_bool_table),        INTENT(IN) :: table       !< true/false table object
    ! local variables
    INTEGER :: dst, i, width, icol
    CHARACTER(LEN=MAX_TITLE_LEN)   :: title_str
    CHARACTER(LEN=MAX_ROWNAME_LEN) :: entry_str
    CHARACTER(LEN=64)              :: format_str

    dst = 0 ! output file / output stream

    ! construct and print the table header
    WRITE (dst,*) "" ! new line
    ! first column:
    DO i=1,table%rowname_len
      WRITE (dst,'(a)', advance='no') " "
    END DO
    WRITE (dst,'(a)', advance='no') DELIMITER
    ! other, TRUE/FALSE columns:
    DO icol = 1, table%n_columns
      title_str = table%column(icol)%title
      width     = table%column(icol)%width
      format_str = "(a"//TRIM(int2string(width))//',a)'
      WRITE (dst,TRIM(format_str), advance='no') ADJUSTL(title_str(1:width)), DELIMITER
    END DO
    WRITE (dst,*) "" ! new line
    WRITE (dst,*) "" ! new line

    ! write the table contents
    LOOP : DO i=1,table%n_rows
      ! first column:
      format_str = "(a"//TRIM(int2string(table%rowname_len))//',a)'
      title_str  = table%rowname(i)
      WRITE (dst,TRIM(format_str), advance='no') ADJUSTL(title_str(1:table%rowname_len)), DELIMITER
      ! other, TRUE/FALSE columns:
      DO icol = 1, table%n_columns
        entry_str = get_bool_table_entry(table%rowname(i), table%column(icol))
        format_str = "(a"//TRIM(int2string(table%column(icol)%width))//',a'//")"
        WRITE (dst,TRIM(format_str), advance='no') entry_str, DELIMITER
      END DO
      WRITE (dst,*) "" ! new line
    END DO LOOP
    WRITE (dst,*) "" ! new line
  END SUBROUTINE print_bool_table

END MODULE mo_util_bool_table
