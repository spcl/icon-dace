! Read namelists, make sanity checks specific to each namelist and make
! a cross check once all namelists of a component are available.
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

MODULE mo_wave_read_namelists

  USE mo_mpi,                   ONLY: my_process_is_stdio
  USE mo_namelist,              ONLY: open_nml_output, close_nml_output
  USE mo_nml_annotate,          ONLY: log_nml_settings
  USE mo_time_nml,              ONLY: read_time_namelist
  USE mo_parallel_nml,          ONLY: read_parallel_namelist
  USE mo_run_nml,               ONLY: read_run_namelist
  USE mo_gribout_nml,           ONLY: read_gribout_namelist
  USE mo_io_nml,                ONLY: read_io_namelist
  USE mo_name_list_output_init, ONLY: read_name_list_output_namelists
  USE mo_grid_nml,              ONLY: read_grid_namelist
  USE mo_grid_config,           ONLY: init_grid_configuration
  USE mo_coupling_nml,          ONLY: read_coupling_namelist
  USE mo_extpar_nml,            ONLY: read_extpar_namelist
  USE mo_wave_nml,              ONLY: read_wave_namelist
  USE mo_interpol_nml,          ONLY: read_interpol_namelist
  USE mo_energy_propagation_nml,ONLY: read_energy_propagation_nml
  USE mo_advection_nml,         ONLY: read_transport_namelist

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_wave_namelists

CONTAINS

  !---------------------------------------------------------------------
  !>
  SUBROUTINE read_wave_namelists(wave_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: wave_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    INTEGER :: tlen
    LOGICAL :: is_stdio

    is_stdio = my_process_is_stdio()
    IF (is_stdio) CALL open_nml_output('NAMELIST_ICON_output_wave')

    CALL read_time_namelist(TRIM(shr_namelist_filename))

    tlen = LEN_TRIM(wave_namelist_filename)
    CALL read_parallel_namelist       (wave_namelist_filename(1:tlen))

    CALL read_run_namelist            (wave_namelist_filename(1:tlen))

    CALL read_io_namelist             (wave_namelist_filename(1:tlen))

    CALL read_name_list_output_namelists (wave_namelist_filename(1:tlen))

    CALL read_grid_namelist           (wave_namelist_filename(1:tlen))
    CALL read_interpol_namelist       (wave_namelist_filename(1:tlen))

    CALL init_grid_configuration()

    CALL read_energy_propagation_nml  (wave_namelist_filename(1:tlen))
    ! temporary hack as long as llsq_svd is used directly from advecton_config
    CALL read_transport_namelist      (wave_namelist_filename(1:tlen))

    CALL read_wave_namelist           (wave_namelist_filename(1:tlen))

    CALL read_extpar_namelist         (wave_namelist_filename(1:tlen))

    CALL read_gribout_namelist        (wave_namelist_filename(1:tlen))

    CALL read_coupling_namelist       (wave_namelist_filename(1:tlen))

    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------

    IF (is_stdio) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (is_stdio) CALL log_nml_settings("nml.wave.log")



  END SUBROUTINE read_wave_namelists

END MODULE mo_wave_read_namelists
