!
! mo_art_tagging
! This module provides utility functions and tools in order to
! create 'tagged' tracers. These are duplicates of the original tracers
! with distinct changes (i.e. different starting time of a volcanic eruption)
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

MODULE mo_art_tagging
! ICON
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_impl_constants,                ONLY: SUCCESS

! ART
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_tagging'

  PUBLIC  :: get_number_tagged_tracer

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
INTEGER FUNCTION get_number_tagged_tracer(key_value_store) RESULT(ntag)
  TYPE(t_key_value_store),INTENT(in) :: &
    &  key_value_store                    !< Storage container with metadata
! Local variables
  INTEGER                    :: &
    &  jt, jtmax,               & !< Counter for tags, maximum number of tags
    &  ierror                     !< Error value from key_value_store%get
  CHARACTER(LEN=3)           :: &
    &  jt_string                  !< jt as character
  CHARACTER(:),ALLOCATABLE   :: &
    &  c_tmp

  jtmax = 999
  ntag  = 0

  DO jt = 1, jtmax
    WRITE(jt_string,'(I3)') jt
    IF (jt < 100) jt_string = '0'//TRIM(ADJUSTL(jt_string))
    IF (jt < 10)  jt_string = '0'//TRIM(ADJUSTL(jt_string))
    CALL key_value_storage_as_string(key_value_store,TRIM('tag'//jt_string),c_tmp,ierror)

    IF(ierror == SUCCESS) THEN
      ntag = jt
    ELSE
      EXIT
    ENDIF
  ENDDO

END FUNCTION get_number_tagged_tracer
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tagging
