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

MODULE TEST_STRING
  USE FORTUTF
  USE mo_util_string
  USE ISO_C_BINDING, ONLY: c_char
  USE helpers, ONLY: open_logfile, open_new_logfile

CONTAINS
  SUBROUTINE TEST_string_to_lower
    CHARACTER(len=10) :: uppercase, lowercase
    CALL TAG_TEST("TEST_to_lower")
    uppercase = 'ALLCAPITAL'
    lowercase = tolower(uppercase)
    CALL STRING_CONTAINS('allcapital', lowercase)
  END SUBROUTINE

  SUBROUTINE TEST_string_low_case
    CHARACTER(len=10) :: testcase
    CALL TAG_TEST("TEST_low_case")
    testcase = 'ALLCAPITAL'
    CALL lowcase(testcase)
    CALL STRING_CONTAINS('allcapital', testcase)
  END SUBROUTINE

  SUBROUTINE TEST_tocompact
    CHARACTER(len=30) :: testcase
    CALL TAG_TEST("TEST_tocompact")
    testcase = ' remove    extra   spaces '
    CALL tocompact(testcase)
    CALL STRING_CONTAINS('remove extra spaces', testcase)
  END SUBROUTINE

  SUBROUTINE TEST_int2string
    CHARACTER(len=6) :: testcase
    INTEGER :: n = 1234
    CALL TAG_TEST("TEST_int2string")
    testcase = int2string(n)
    CALL STRING_CONTAINS('1234', testcase)
  END SUBROUTINE

  SUBROUTINE TEST_float2string
    CHARACTER(len=32) :: testcase
    REAL :: n = 1234.5678
    CALL TAG_TEST("TEST_float2string")
    testcase = real2string(n)
    CALL STRING_CONTAINS('1234.6', testcase)
  END SUBROUTINE

  SUBROUTINE TEST_double2string
    CHARACTER(len=32) :: testcase
    DOUBLE PRECISION :: n = 1234.5678d0
    CALL TAG_TEST("TEST_double2string")
    testcase = real2string(n)
    CALL STRING_CONTAINS('1234.5678', testcase)
  END SUBROUTINE

  SUBROUTINE TEST_logical2string
    CHARACTER(len=10) :: testcase
    LOGICAL :: n = .TRUE.
    CALL TAG_TEST("TEST_logical2string_true")
    testcase = logical2string(n)
    CALL STRING_CONTAINS('T', testcase)

    CALL TAG_TEST("TEST_logical2string_false")
    n = .FALSE.
    testcase = logical2string(n)
    CALL STRING_CONTAINS('F', testcase)
  END SUBROUTINE

  SUBROUTINE TEST_str_replace
    CHARACTER(len=30) :: testcase, keyword, subst, RESULT
    CALL TAG_TEST("TEST_str_replace")
    testcase = 'Hello, world!'
    keyword = 'world'
    subst = 'universe'
    RESULT = str_replace(testcase, keyword, subst)
    CALL STRING_CONTAINS('Hello, universe!', RESULT)
  END SUBROUTINE TEST_str_replace

  SUBROUTINE TEST_remove_duplicates
    CHARACTER(len=30), DIMENSION(5) :: str_list
    INTEGER :: nitems
    CALL TAG_TEST("TEST_remove_duplicates")
    str_list = ['word1', 'word2', 'word1', 'word3', 'word2']
    nitems = SIZE(str_list)
    CALL remove_duplicates(str_list, nitems)
    CALL ASSERT_EQUAL((nitems == 3 .AND. str_list(1) == 'word1' .AND. str_list(2) == 'word2' .AND. str_list(3) == 'word3'), .TRUE.)
  END SUBROUTINE TEST_remove_duplicates

  SUBROUTINE TEST_difference
    CHARACTER(len=30), DIMENSION(5) :: str_list1
    CHARACTER(len=30), DIMENSION(2) :: str_list2
    INTEGER :: nitems1, nitems2
    CALL TAG_TEST("TEST_difference")
    str_list1 = ['word1', 'word2', 'word3', 'word4', 'word5']
    str_list2 = ['word2', 'word4']
    nitems1 = SIZE(str_list1)
    nitems2 = SIZE(str_list2)
    CALL difference(str_list1, nitems1, str_list2, nitems2)
    CALL ASSERT_EQUAL((nitems1 == 3 .AND. str_list1(1) == 'word1' .AND. str_list1(2) == 'word3' &
                       .AND. str_list1(3) == 'word5'), .TRUE.)
  END SUBROUTINE TEST_difference

  SUBROUTINE TEST_new_list
    CHARACTER(len=30), ALLOCATABLE :: str_list(:)
    INTEGER :: nitems
    CALL TAG_TEST("TEST_new_list")
    CALL new_list(str_list, nitems)
    CALL ASSERT_EQUAL(nitems == 0 .AND. ALLOCATED(str_list), .TRUE.)
  END SUBROUTINE TEST_new_list

  SUBROUTINE TEST_add_to_list
    CHARACTER(len=30), ALLOCATABLE, DIMENSION(:) :: str_list1
    CHARACTER(len=30), DIMENSION(2) :: str_list2
    INTEGER :: nitems1, nitems2
    CALL add_to_list(str_list1, nitems1, ['word1', 'word2', 'word3']) ! Implicit allocation
    CALL TAG_TEST("TEST_add_to_list (1)")
    CALL ASSERT_EQUAL(nitems1 == 3 .AND. ALLOCATED(str_list1) .AND. SIZE(str_list1) >= nitems1, .TRUE.)
    str_list2 = ['word2', 'word4']
    nitems2 = SIZE(str_list2)
    CALL add_to_list(str_list1, nitems1, str_list2, nitems2)
    CALL TAG_TEST("TEST_add_to_list (2)")
    CALL ASSERT_EQUAL((nitems1 == 4 .AND. str_list1(1) == 'word1' .AND. &
                       str_list1(2) == 'word2' .AND. str_list1(3) == 'word3' &
                       .AND. str_list1(4) == 'word4'), .TRUE.)
    CALL add_to_list(str_list1, nitems1, ['5', '6', '7', '8', '9'])
    CALL TAG_TEST("TEST_add_to_list (3)")
    CALL ASSERT_EQUAL(nitems1 == 9 .AND. ALLOCATED(str_list1) .AND. SIZE(str_list1) >= nitems1 .AND. &
                      ALL(str_list1(1:nitems1) == [ &
                          'word1', 'word2', 'word3', 'word4', '5    ', '6    ', '7    ', '8    ', '9    ']), .TRUE.)
  END SUBROUTINE TEST_add_to_list

  SUBROUTINE TEST_add_to_list_1
    CHARACTER(len=30), ALLOCATABLE, DIMENSION(:) :: str_list
    INTEGER :: nitems
    CALL add_to_list(str_list, nitems, 'word1') ! Implicit allocation
    CALL TAG_TEST("TEST_add_to_list_1 (1)")
    CALL ASSERT_EQUAL(nitems == 1 .AND. ALLOCATED(str_list) .AND. SIZE(str_list) >= nitems, .TRUE.)
    CALL add_to_list(str_list, nitems, 'word2')
    CALL add_to_list(str_list, nitems, 'word3')
    CALL add_to_list(str_list, nitems, 'word4')
    CALL add_to_list(str_list, nitems, 'word2') ! Try to add duplicate
    CALL TAG_TEST("TEST_add_to_list_1 (2)")
    CALL ASSERT_EQUAL((nitems == 4 .AND. str_list(1) == 'word1' &
                       .AND. str_list(2) == 'word2' .AND. str_list(3) == 'word3' &
                       .AND. str_list(4) == 'word4'), .TRUE.)
    CALL add_to_list(str_list, nitems, ['5', '6', '7', '8']) ! Fill to allocated size
    CALL add_to_list(str_list, nitems, '9') ! Expand list
    CALL TAG_TEST("TEST_add_to_list_1 (3)")
    CALL ASSERT_EQUAL(nitems == 9 .AND. ALLOCATED(str_list) .AND. SIZE(str_list) >= nitems .AND. &
                      ALL(str_list(1:nitems) == [ &
                          'word1', 'word2', 'word3', 'word4', '5    ', '6    ', '7    ', '8    ', '9    ']), .TRUE.)
  END SUBROUTINE TEST_add_to_list_1

  SUBROUTINE TEST_remove_whitespace
    CHARACTER(len=30) :: str
    CALL TAG_TEST("TEST_remove_whitespace")
    str = "  h e l l o  "
    CALL STRING_CONTAINS(remove_whitespace(str), "hello")
  END SUBROUTINE TEST_remove_whitespace

  SUBROUTINE TEST_c2f_char
    CHARACTER(LEN=:), ALLOCATABLE :: c
    CHARACTER(KIND=c_char), DIMENSION(5) :: s
    CALL TAG_TEST("TEST_c2f_char")
    s = ['h', 'e', 'l', 'l', 'o']
    CALL c2f_char(c, s)
    CALL STRING_CONTAINS(c, "hello")
  END SUBROUTINE TEST_c2f_char

  SUBROUTINE TEST_charArray_equal_array
    CHARACTER(KIND=c_char), DIMENSION(5) :: stringA, stringB
    CALL TAG_TEST("TEST_charArray_equal_array")
    stringA = ['h', 'e', 'l', 'l', 'o']
    stringB = ['h', 'e', 'l', 'l', 'o']
    CALL ASSERT_EQUAL((charArray_equal(stringA, stringB)), .TRUE.)
  END SUBROUTINE TEST_charArray_equal_array

  SUBROUTINE TEST_charArray_equal_char
    CHARACTER(KIND=c_char), DIMENSION(5) :: stringA
    CHARACTER(len=30) :: stringB
    LOGICAL:: res
    CALL TAG_TEST("TEST_charArray_equal_char")
    stringA = ['h', 'e', 'l', 'l', 'o']
    stringB = "hello"
    CALL ASSERT_EQUAL(charArray_equal(stringA, stringB), .TRUE.)
  END SUBROUTINE TEST_charArray_equal_char

  SUBROUTINE TEST_charArray_toLower
    CHARACTER(KIND=c_char), DIMENSION(5) :: string
    CALL TAG_TEST("TEST_charArray_toLower")
    string = ['H', 'E', 'L', 'L', 'O']
    CALL charArray_toLower(string)
    CALL ASSERT_EQUAL((ALL(string == ['h', 'e', 'l', 'l', 'o'])), .TRUE.)
  END SUBROUTINE TEST_charArray_toLower

  SUBROUTINE TEST_pretty_print_string_list
    CHARACTER(len=30), DIMENSION(5) :: list
    CHARACTER(len=200) :: logfile, log_in_file
    CALL TAG_TEST("TEST_pretty_print_string_list")
    list = ['word1', 'word2', 'word3', 'word4', 'word5']
    logfile = 'print_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL pretty_print_string_list(list, 20, nerr, "Prefix: ")
    CLOSE (nerr)
    CALL open_logfile(nerr, TRIM(logfile))
    ! first line of file is message from abort
    READ (nerr, '(A)') log_in_file

    CALL STRING_CONTAINS(log_in_file, 'Prefix:')
    CLOSE (nerr)
  END SUBROUTINE TEST_pretty_print_string_list

  SUBROUTINE TEST_find_trailing_number
    CHARACTER(len=30) :: str
    INTEGER :: pos
    CALL TAG_TEST("TEST_find_trailing_number")
    str = "word123"
    pos = find_trailing_number(str)
    CALL ASSERT_EQUAL(pos, 5)
  END SUBROUTINE TEST_find_trailing_number

  SUBROUTINE TEST_tohex
    CHARACTER(len=30) :: str, RESULT
    CALL TAG_TEST('TEST_tohex')
    str = "Hello, world!"
    RESULT = tohex(str)
    CALL STRING_CONTAINS(RESULT, "48 65 6c 6c 6f 2c 20 77 6f 72")
  END SUBROUTINE TEST_tohex

  SUBROUTINE TEST_sort_and_compress_list
    INTEGER, DIMENSION(6) :: list
    CHARACTER(len=30) :: RESULT
    CALL TAG_TEST('TEST_sort_and_compress_list')
    list = [3, 1, 2, 5, 7, 6]
    CALL sort_and_compress_list(list, RESULT)
    CALL STRING_CONTAINS(RESULT, "1-3,5-7")
  END SUBROUTINE TEST_sort_and_compress_list

  SUBROUTINE TEST_insert_group
    CHARACTER(len=7), DIMENSION(5) :: varlist
    INTEGER :: nused, ninserted
    CHARACTER(len=7) :: group_name
    CHARACTER(len=7), DIMENSION(3) :: group_list

    varlist(1) = "u      "
    varlist(2) = "v      "
    varlist(3) = "tracers"
    varlist(4) = "       "
    varlist(5) = "       "

    nused = 3
    group_name = "tracers"
    group_list = ["Q1", "Q2", "Q3"]

    CALL insert_group(varlist, nused, group_name, group_list, ninserted)

    CALL TAG_TEST('TEST_insert_group_nused')
    CALL ASSERT_EQUAL(nused, 5)
    CALL TAG_TEST('TEST_insert_group_ninserted')
    CALL ASSERT_EQUAL(ninserted, 3)
  END SUBROUTINE TEST_insert_group

  SUBROUTINE TEST_associate_keyword_with_keywords
    CHARACTER(len=100), PARAMETER :: filename = "<path>/bin/<prefix>grid.nc"
    TYPE(t_keyword_list), POINTER :: keywords => NULL()
    CHARACTER(len=100) :: result_str

    CALL TAG_TEST('TEST_associate_keyword_with_keywords')
    CALL associate_keyword("<path>", "/usr/local", keywords)
    CALL associate_keyword("<prefix>", "exp01_", keywords)

    result_str = with_keywords(keywords, filename)

    CALL STRING_CONTAINS(result_str, "/usr/local/bin/exp01_grid.nc")
  END SUBROUTINE TEST_associate_keyword_with_keywords

END MODULE TEST_STRING
