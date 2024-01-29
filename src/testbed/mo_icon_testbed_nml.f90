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

MODULE mo_icon_testbed_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: filename_max
  USE mo_icon_testbed_config, ONLY: &
    & config_testbed_model        => testbed_model,       &
    & config_testbed_iterations   => testbed_iterations,  &
    & config_calculate_iterations => calculate_iterations,&
    & config_no_of_blocks         => no_of_blocks,        &
    & config_no_of_layers         => no_of_layers,        &
    & config_testfile_3D_time          => testfile_3D_time,         &
    & config_testfile_2D_time         => testfile_2D_time,        &
    & null_model
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_icon_testbed_namelist

       
  CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_icon_testbed_namelist( filename )

    ! ------------------------------------------------------------------------
    INTEGER :: testbed_model
    INTEGER :: testbed_iterations
    INTEGER :: calculate_iterations
    INTEGER :: no_of_blocks, no_of_layers
    INTEGER :: iunit

   CHARACTER(LEN=filename_max) :: testfile_3D_time(2), testfile_2D_time(2)
    
    NAMELIST /testbed_nml/ testbed_model, testbed_iterations, calculate_iterations, &
      & no_of_blocks, no_of_layers, testfile_3D_time, testfile_2D_time

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat
    !0!CHARACTER(len=*), PARAMETER ::   &
    !0!        &  method_name = 'mo_parallel_nml:read_parallel_namelist'

    !--------------------------------------------
    ! set default values
    !--------------------------------------------
    testbed_model         = null_model
    testbed_iterations   = 10
    calculate_iterations = 10
    no_of_blocks         = 16
    no_of_layers         = 80
    testfile_3D_time     = ""
    testfile_2D_time     = ""

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes) 
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('testbed_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, testbed_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, testbed_nml, iostat=istat)                          ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, testbed_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml
    
    !-----------------------------------------------------
    ! fill_config_testbed
    config_testbed_model        = testbed_model
    config_testbed_iterations   = testbed_iterations
    config_calculate_iterations = calculate_iterations
    config_no_of_blocks         = no_of_blocks
    config_no_of_layers         = no_of_layers
    config_testfile_3D_time     = testfile_3D_time
    config_testfile_2D_time     = testfile_2D_time
    
  END SUBROUTINE read_icon_testbed_namelist
  !-------------------------------------------------------------------------

END MODULE mo_icon_testbed_nml
