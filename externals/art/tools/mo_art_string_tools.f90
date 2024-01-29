!
! mo_art_string_tools
! This module provides routines for parsing arrays and production elements from
! XML strings
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

MODULE mo_art_string_tools
  ! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_key_value_store,               ONLY: t_key_value_store
  ! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_get_no_elements_in_string
  PUBLIC :: art_split_string_in_array
  PUBLIC :: art_parse_production
  PUBLIC :: key_value_storage_as_string

CONTAINS
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_get_no_elements_in_string(num_elem,string,sep,positions_sep)
!<
! SUBROUTINE art_get_no_elements_in_string
! Get number of elements if string is separated via sep
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: &
    &  string,                    &
    &  sep
  INTEGER, INTENT(out) ::             &
    &  num_elem
  INTEGER, INTENT(out), ALLOCATABLE, OPTIONAL :: &
    &  positions_sep(:)

  ! local variables
  INTEGER ::       &
    & number_of_seps_in_str, &
    & i,                     &
    & pos(50)
  

  IF (LEN(sep) > 1) THEN
    CALL finish('mo_art_string_tools:art_get_number_of_elements',  &
           &    'parameter sep is too long (must have length 1).')
  END IF

  number_of_seps_in_str = 0

  DO i = 1,LEN_TRIM(string)
    IF (string(i:i) == sep) THEN
      number_of_seps_in_str = number_of_seps_in_str + 1
      IF (number_of_seps_in_str <= SIZE(pos)) THEN
        pos(number_of_seps_in_str) = i 
      END IF
    END IF
  END DO

 
  num_elem = number_of_seps_in_str + 1

  IF (PRESENT(positions_sep)) THEN
    IF (number_of_seps_in_str > 50) THEN
      CALL finish('mo_art_string_tools:art_get_number_of_elements', &
             &    'number of '//sep//' in string must not exceed 50.')
    ELSE IF (number_of_seps_in_str >= 1) THEN
      ALLOCATE(positions_sep(number_of_seps_in_str))
      positions_sep(:) = pos(1:number_of_seps_in_str)
    END IF
  END IF


END SUBROUTINE art_get_no_elements_in_string
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_split_string_in_array(str_arr,string,sep)
!<
! SUBROUTINE art_split_string_in_array
! Returns a string array of string separated via sep
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  CHARACTER(LEN=*), ALLOCATABLE, INTENT(out) :: &
    &  str_arr(:)
  CHARACTER(LEN=*), INTENT(in) :: &
    &  string,                    &
    &  sep

  ! local variables
  INTEGER :: &
    &  num_elem, i
  INTEGER, ALLOCATABLE :: &
    &  positions_sep(:)

  CALL art_get_no_elements_in_string(num_elem,TRIM(string),sep,positions_sep)


  ALLOCATE(str_arr(num_elem))

  IF (num_elem == 1) THEN
    str_arr(1) = TRIM(ADJUSTL(string))
  ELSE
    DO i = 1,num_elem
      IF (i == 1) THEN
        str_arr(i) = string(1:positions_sep(i)-1)
      ELSE IF (i < num_elem) THEN
        str_arr(i) = string(positions_sep(i-1)+1:positions_sep(i)-1)
      ELSE
        str_arr(i) = string(positions_sep(i-1)+1:)
      END IF
    END DO
  END IF
END SUBROUTINE art_split_string_in_array
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_parse_production(factors,tracer_names,string)
!<
! SUBROUTINE art_parse_production
! Returns an array of factors and tracer names from a string of format:
! 0.007*TRSO2;0.993*TROCS
!
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  REAL(wp), ALLOCATABLE, INTENT(out) :: &
    &  factors(:)
  CHARACTER(LEN=*), ALLOCATABLE, INTENT(out) :: &
    &  tracer_names(:)
  CHARACTER(LEN=*), INTENT(in) :: &
    &  string

  ! local variables
  CHARACTER(LEN = 30), ALLOCATABLE :: &
    &  str_arr(:)
  CHARACTER(LEN=30) :: &
    &  prod_str
  CHARACTER(LEN = IART_VARNAMELEN) :: &
    &  tracer_name
  INTEGER ::      &
    &  ioerr,     &
    &  num_prods, &
    &  i,         &
    &  index_star


  CALL art_split_string_in_array(str_arr,TRIM(string),';')
 
  num_prods = SIZE(str_arr)
  ALLOCATE(factors(num_prods))
  ALLOCATE(tracer_names(num_prods))

  DO i = 1,SIZE(str_arr)
    prod_str = str_arr(i)
    index_star = INDEX(prod_str,'*')

    IF (index_star > 0) THEN
      READ(prod_str(1:index_star-1),*,iostat=ioerr) factors(i)
  
      IF (ioerr /= 0) THEN
        CALL finish('mo_art_string_tools:art_parse_production',  &
               &    'Could not read factor from '//prod_str(1:index_star-1)  &
               &  //' which is part of '//TRIM(string)//'.')
      END IF
    ELSE
      factors(i) = 1
    END IF

    tracer_names(i) = TRIM(ADJUSTL(prod_str(index_star+1:)))
    tracer_name = tracer_names(i)
    IF (tracer_name(1:2) /= 'TR') THEN
      CALL finish('mo_art_string_tools:art_parse_production',  &
             &    'parsed tracer name '//TRIM(tracer_names(i))  &
             &  //' does not begin with TR.')
    END IF
  END DO

END SUBROUTINE art_parse_production

SUBROUTINE key_value_storage_as_string(key_value_store, key, val, ierror)
  IMPLICIT NONE
  TYPE(t_key_value_store), INTENT(in) :: &
    &  key_value_store
  CHARACTER(LEN = *), INTENT(in) :: &
    &  key
  CHARACTER(:), ALLOCATABLE, INTENT(out) :: &   
    &  val
  INTEGER, INTENT(out), OPTIONAL :: &
    &  ierror

  IF (PRESENT(ierror)) THEN
    CALL key_value_store%get(key,val,ierror)
  ELSE
    CALL key_value_store%get(key,val)
  END IF
  
END SUBROUTINE key_value_storage_as_string
!!
!!-------------------------------------------------------------------
!!
END MODULE mo_art_string_tools
