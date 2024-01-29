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

MODULE TEST_mo_hash_table
  USE FORTUTF
  USE mo_util_texthash, ONLY: text_hash_c

  TYPE :: t_Key
    INTEGER :: val
  END TYPE

  TYPE :: t_Value
    CHARACTER(5) :: val
  END TYPE

CONTAINS

  SUBROUTINE TEST_hash_table_int
    USE mo_hash_table
    USE mo_util_texthash, ONLY: text_hash, text_isEqual

    TYPE(t_HashTable), POINTER :: hashTable
    TYPE(t_HashIterator) :: iterator
    CLASS(*), POINTER :: key, val
    CHARACTER(8), TARGET, ALLOCATABLE :: keyword(:)
    INTEGER, TARGET, ALLOCATABLE :: int_value(:)
    LOGICAL :: success
    INTEGER :: i
    INTEGER :: table_size = 20

    ALLOCATE (hashTable)
    hashTable => hashTable_make(text_hash, text_isEqual)

    ALLOCATE (keyword(table_size))
    ALLOCATE (int_value(table_size))

    DO i = 1, table_size
      WRITE (keyword(i), '(A,I2)') 'nr:', i
      key => keyword(i)
      int_value(i) = i
      val => int_value(i)
      CALL hashTable%setEntry(key, val)
    END DO

    CALL TAG_TEST("TEST_hash_table_int")
    CALL ASSERT_EQUAL(hashTable%getEntryCount(), table_size)

  END SUBROUTINE TEST_hash_table_int

END MODULE TEST_mo_hash_table
