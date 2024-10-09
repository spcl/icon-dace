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

MODULE mo_icon_testbed_config

  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: filename_max

  IMPLICIT NONE

  PRIVATE
  
  ! Exported parameters
  PUBLIC :: null_model, test_coupler_model, test_jitter_model, test_halo_communication, &
    & test_netcdf_read_model, testbed_ocean_model, test_gather_communication, &
    & test_exchange_communication, test_bench_exchange_data_mult
    
  ! Exported variables
  PUBLIC :: testbed_model
  PUBLIC :: testbed_iterations
  PUBLIC :: calculate_iterations
  PUBLIC :: no_of_blocks, no_of_layers
  PUBLIC :: testfile_3D_time, testfile_2D_time
  
  ! ---------------
  ! the tesbed modes
  INTEGER, PARAMETER :: null_model                   = 0  ! does nothing
  INTEGER, PARAMETER :: test_coupler_model           = 1  ! test the coupler
  INTEGER, PARAMETER :: test_jitter_model            = 3  ! test the jitter
  INTEGER, PARAMETER :: test_halo_communication      = 4  ! test the mpi communication
  INTEGER, PARAMETER :: test_netcdf_read_model       = 6
  INTEGER, PARAMETER :: testbed_ocean_model          = 7
  INTEGER, PARAMETER :: test_gather_communication = 8
  INTEGER, PARAMETER :: test_exchange_communication = 9
  INTEGER, PARAMETER :: test_bench_exchange_data_mult = 10
  

  INTEGER  :: testbed_model
  ! ---------------
  INTEGER  :: testbed_iterations
  INTEGER  :: calculate_iterations
  INTEGER  :: no_of_blocks, no_of_layers

  CHARACTER(LEN=filename_max) :: testfile_3D_time(2), testfile_2D_time(2)


CONTAINS
  
  !-------------------------------------------------------------------------
  SUBROUTINE check_icon_testbed_config()
 
  END SUBROUTINE check_icon_testbed_config
  !-------------------------------------------------------------------------
    
END MODULE mo_icon_testbed_config
