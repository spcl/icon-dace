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

MODULE mo_comin_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,   &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_master_config,       ONLY: isRestart
  USE mo_comin_config,        ONLY: comin_config






  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_comin_namelist









CONTAINS
  !> Read namelist settings for ICON ComIn
  !
  !  Sets defaults for various `comin_nml` user settings, then opens
  !  and reads the namelist file `filename`, overwriting some
  !  defaults.
  !
  SUBROUTINE read_comin_namelist( filename )
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !
  END SUBROUTINE read_comin_namelist

END MODULE mo_comin_nml
