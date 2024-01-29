! @brief Routines for the initialisation of YAC on IO processes
!
! The purpose of routines construct_io_coupling and destruct_io_coupling is
! to initialise YAC on asynchronous IO processes. The YAC component definition
! and search are collective operations in the MPI sense. Thus all MPI processes
! must call yac_fdef_comp(s) and yac_fenddef.
! Currently, the IO processes do not receive any fields via YAC for output.
! Therefore, yac_fenddef is called with an empty list of fields.
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

MODULE mo_io_coupling_frame

  USE mo_coupling_config,           ONLY: is_coupled_run
  USE mo_impl_constants,            ONLY: MAX_CHAR_LENGTH
  USE mo_exception,                 ONLY: message, message_text
  USE mo_yac_finterface,            ONLY: yac_fdef_comp, yac_fenddef

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_io_coupling
  PUBLIC :: destruct_io_coupling
 
CONTAINS

  SUBROUTINE construct_io_coupling ( comp_name )

    CHARACTER(LEN=*), INTENT(IN) :: comp_name
    CHARACTER(*), PARAMETER :: routine = "mo_io_coupling_frame:construct_io_coupling"

    INTEGER :: comp_ids(1)

    IF ( is_coupled_run() ) THEN

      WRITE(message_text,*) "YAC initialisation for asynchronous I/O processes of type ", &
        &                   TRIM(comp_name)
      CALL message(routine, message_text)

      ! Inform YAC about what we are
      CALL yac_fdef_comp ( TRIM(comp_name), comp_ids(1) )

      ! IO processes need to participate in the YAC search
      CALL yac_fenddef ( )
    ENDIF

  END SUBROUTINE construct_io_coupling

  SUBROUTINE destruct_io_coupling ( comp_name )

    CHARACTER(LEN=*), INTENT(IN) :: comp_name
 
    CHARACTER(*), PARAMETER :: routine = "mo_io_coupling_frame:construct_io_coupling"

    IF ( is_coupled_run() ) THEN

      WRITE(message_text,*) "YAC termination of I/O process of type ", &
        &                   TRIM(comp_name)
      CALL message(routine, message_text)

    ENDIF

  END SUBROUTINE destruct_io_coupling

END MODULE mo_io_coupling_frame
