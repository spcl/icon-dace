! Routines to access basic coupler functionality
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

!----------------------------
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!



!----------------------------
MODULE mo_coupling

  USE mo_io_units, ONLY: nerr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_coupler
  PUBLIC :: finalize_coupler
  PUBLIC :: coupler_config_files_exist

  PUBLIC :: lyac_very_1st_get

  CHARACTER(LEN=*), PARAMETER :: yaml_filename = "coupling.yaml"

  LOGICAL :: config_files_have_been_checked = .FALSE.
  LOGICAL :: config_files_exist = .FALSE.
  LOGICAL :: yac_is_initialised = .FALSE.

  LOGICAL, SAVE :: lyac_very_1st_get

  CHARACTER(*), PARAMETER :: modname = "mo_coupling"

CONTAINS

  LOGICAL FUNCTION coupler_config_files_exist()

    LOGICAL :: yaml_exists

    IF (config_files_have_been_checked) THEN

      coupler_config_files_exist = config_files_exist

    ELSE

      INQUIRE(FILE=TRIM(ADJUSTL(yaml_filename)), EXIST=yaml_exists)

      config_files_have_been_checked = .TRUE.
      config_files_exist = yaml_exists
      coupler_config_files_exist = config_files_exist

    END IF

  END FUNCTION

  SUBROUTINE init_coupler(opt_yac_comm, opt_world_comm)

    INTEGER, INTENT(IN), OPTIONAL :: opt_yac_comm
    INTEGER, INTENT(INOUT), OPTIONAL :: opt_world_comm



  END SUBROUTINE init_coupler

  SUBROUTINE finalize_coupler


  END SUBROUTINE finalize_coupler

  SUBROUTINE print_info_stderr (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    WRITE (nerr,'(4a)') " ", TRIM(name), ": ", TRIM(text)

  END SUBROUTINE print_info_stderr

END MODULE mo_coupling
