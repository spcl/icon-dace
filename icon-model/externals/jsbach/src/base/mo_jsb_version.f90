!> Contains methods to define and print JSBACH version number and label a simulation.
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

#ifdef __xlC__
@PROCESS STRICT
#endif

MODULE mo_jsb_version
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: jsbach_version, jsbach_init_version, jsbach_label_run

  CHARACTER(len=  6), PARAMETER :: jsbach_version = '4.20p9' !< JSBACH version number

  CHARACTER(len=256) :: executable                         !< Name of executable
  CHARACTER(len= 80) :: label(4)                           !< Space for storing model information
  LOGICAL            :: version_printed = .FALSE.

CONTAINS

  !>
  !! Define JSBACH version and labels
  !!
  !! Define JSBACH version number and get some information about a simulation
  !! at model initialization.
  !!
  SUBROUTINE jsbach_init_version(standalone)

!!$    USE mo_netcdf, ONLY: global_att, put_att

    LOGICAL, INTENT(in), OPTIONAL :: standalone !< true: offline JSBACH, false: ECHAM/JSBACH

    !LOGICAL             :: l_standalone !< true: offline JSBACH, false: ECHAM/JSBACH
    !CHARACTER (len=  8) :: ydate        !< Date of simulation
    !CHARACTER (len= 10) :: ytime        !< Time of simulation
    !CHARACTER (len=256) :: os_name      !< Name of OS for simulation
    !CHARACTER (len=256) :: user_name    !< User that runs simulation
    !CHARACTER (len=256) :: host_name    !< Host name for simulation

    !INTEGER :: i_lena, i_lenb, i_lenc
    INTEGER :: i

    !  External subroutines
!!$    EXTERNAL :: util_os_system, util_user_name, util_node_name

    !  Executable statements

    !l_standalone = .FALSE.
    !IF (PRESENT(standalone)) l_standalone = standalone
    IF (PRESENT(standalone)) CONTINUE ! to avoid compiler warnings about dummy arguments not being used

!!$    IF (l_standalone) THEN
    IF (.NOT. version_printed) THEN
      ! Print version
      CALL message ('','')
      CALL message ('','')
      message_text = '==========================================================='
      CALL message ('',TRIM(message_text))
      CALL message ('','')

      message_text = '  JSBACH - Version '//jsbach_version
      CALL message ('',TRIM(message_text))
      message_text = '  Copyright by Max-Planck-Institute for Meteorology, 2015'
      CALL message ('',TRIM(message_text))
      IF (INDEX(executable, 'icon') /= 0) THEN
        message_text = '  License subject to the DWD and MPI-M-Software-License-Agreement'
        CALL message ('',TRIM(message_text))
      END IF

      CALL message ('','')
      message_text = '==========================================================='
      CALL message ('',TRIM(message_text))
      CALL message ('','')
      CALL message ('','')
      version_printed = .TRUE.
    END IF


!!$    IF (l_standalone) THEN
!!$      os_name   = 'unknown'
!!$      user_name = 'unknown'
!!$      host_name = 'unknown'
!!$
!!$      CALL util_os_system (os_name,   i_lena)
!!$      CALL util_user_name (user_name, i_lenb)
!!$      CALL util_node_name (host_name, i_lenc)
!!$
!!$      CALL DATE_AND_TIME(ydate, ytime)
!!$    END IF

    DO i = 1, SIZE(label)
      label(i) = ' '
    ENDDO

    WRITE (label(1), '(a)') ' Surface model version: JSBACH '//jsbach_version
!!$    IF (l_standalone) THEN
!!$      WRITE (label(2), '(a)') ' Date - '//ydate(1:8)//' Time - '//ytime(1:6)
!!$      WRITE (label(3), '(a)') ' '//user_name(1:i_lenb)//' on '//host_name(1:i_lenc)
!!$      WRITE (label(4), '(a)') ' '//os_name(1:i_lena)
!!$    END IF

    ! set global attributes to be written to NetCDF file
!!$    CALL put_att (global_att,'jsbach_version',jsbach_version)
!!$    IF (l_standalone) THEN
!!$      CALL put_att (global_att,'date_time',ydate(1:8)//' '//ytime(1:6))
!!$      CALL put_att (global_att,'operating_system',os_name(1:i_lena))
!!$      CALL put_att (global_att,'user_name',user_name(1:i_lenb))
!!$      CALL put_att (global_att,'host_name',host_name(1:i_lenc))
!!$   END IF

  END SUBROUTINE jsbach_init_version

  !>
  !! Label a simulation
  !!
  !! Write out details of a forecast run after the setup is complete, just before
  !! computing the first timestep
  !!
  SUBROUTINE jsbach_label_run(standalone)

    USE mo_jsb_parallel, ONLY: my_process_is_stdio, p_io, p_bcast, mpi_comm
    USE mo_util_string,  ONLY: separator

    LOGICAL, INTENT(in), OPTIONAL :: standalone !< true: offline JSBACH, false: ECHAM/JSBACH

    LOGICAL :: l_standalone !< true: offline JSBACH, false: ECHAM/JSBACH
    INTEGER :: i

    ! Executable statements

    l_standalone = .FALSE.
    IF (PRESENT(standalone)) l_standalone = standalone

    ! Get model name from environment
    IF (my_process_is_stdio()) THEN
      CALL get_command_argument(0, executable)
    END IF

    CALL p_bcast (executable, p_io, comm=mpi_comm)

    IF (l_standalone) THEN
      CALL message ('','')
      CALL message('',separator)
      CALL message ('',' Model: '//TRIM(executable))
      CALL message('',separator)
    END IF

    ! Type of run
    IF (l_standalone) CALL message('',separator)
    DO i = 1, SIZE(label)
      IF (label(i) /= ' ') THEN
        WRITE (message_text,'(a)') label(i)
        CALL message ('',message_text)
      END IF
    ENDDO

    IF (l_standalone) THEN
      CALL message('',separator)
      CALL message('','')
    END IF

  END SUBROUTINE jsbach_label_run

#endif
END MODULE mo_jsb_version
