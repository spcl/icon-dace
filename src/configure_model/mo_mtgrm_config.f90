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

! Setup of meteogram output.

MODULE mo_meteogram_config

  USE mo_impl_constants, ONLY: max_dom
  USE mo_math_types,     ONLY: t_geographical_coordinates
  USE mo_exception,      ONLY: message

  IMPLICIT NONE

!orig  INTEGER, PARAMETER :: MAX_NAME_LENGTH      =   48  !<  max. name string length   
!pa  
  INTEGER, PARAMETER :: MAX_NAME_LENGTH      =   256  !<  max. name string length   
!pa
  INTEGER, PARAMETER :: MAX_NUM_STATIONS     =   75  !<  max. number of meteogram locations (global)
                                                     ! NOTE: MAX_NUM_STATIONS=80 did not work with Intel and GCC.
                                                     ! Reading the meteogram namelist from the restart file
                                                     ! resulted in END OF FILE error in mo_mtgrm_nml.
  INTEGER, PARAMETER :: FTYPE_NETCDF         =    1
  INTEGER, PARAMETER :: MAX_NVARS            =  120  ! max. no. of meteogram variables


  !--------------------------------------------------------------------------
  ! Derived type
  !--------------------------------------------------------------------------

  !>
  !! Input data specifying a meteogram station.
  !!
  TYPE t_station_list
    CHARACTER(len=MAX_NAME_LENGTH)   :: zname    !< station name
    TYPE(t_geographical_coordinates) :: location !< geographical position
  END TYPE t_station_list

  !>
  !! TODO[FP] : Make some fixed-size arrays ALLOCATABLE?
  !!
  TYPE t_meteogram_output_config
    LOGICAL                        :: lenabled     !< Flag. True for output.
    CHARACTER(len=MAX_NAME_LENGTH) :: zprefix      !< file prefix string
    INTEGER                        :: ftype        !< file type (NetCDF, ...)
    LOGICAL                        :: ldistributed !< Flag. Separate files for each PE

    ! Initial step and step interval for meteogram output
    INTEGER                        :: n0_mtgrm, ninc_mtgrm

    ! activate output of some selected tile specific fields
    LOGICAL                        :: loutput_tiles
  
    !> maximal number of steps retained until flush
    INTEGER                        :: max_time_stamps
    !> warn user about (synchronizing) flush?
    LOGICAL                        :: silent_flush

    !> mpi_comm_world rank that writes the data for this meteogram file
    INTEGER                        :: io_proc_id

    !> instead of clobbering existing file, append instead
    LOGICAL                        :: append_if_exists
    ! Specification of meteogram stations.
    ! Note: This info is the same for all patches.
    INTEGER                           :: nstations
    TYPE(t_station_list), ALLOCATABLE :: station_list(:)

    ! Positive-list of variables (optional). Only variables contained
    ! in this list are included in this meteogram. If the default list
    ! is not changed by user input, then all available variables are
    ! added to the meteogram
    CHARACTER(len=MAX_NAME_LENGTH)    :: var_list(MAX_NVARS)
  END TYPE t_meteogram_output_config


  ! parameters specifying meteogram output
  TYPE(t_meteogram_output_config), TARGET   :: meteogram_output_config(1:max_dom)

  PUBLIC :: t_meteogram_output_config, t_station_list
  PUBLIC :: meteogram_output_config
  PUBLIC :: FTYPE_NETCDF, MAX_NAME_LENGTH, MAX_NUM_STATIONS, MAX_NVARS
  PUBLIC :: check_meteogram_configuration

  !-------------------------------------------------------------------------

CONTAINS

  SUBROUTINE check_meteogram_configuration(num_io_procs)
    INTEGER, INTENT(IN) :: num_io_procs
    ! local variables
    CHARACTER(*), PARAMETER :: routine &
         = "mo_mtgrm_config:check_meteogram_configuration"
    INTEGER :: idom

    ! Asynchronous output does not work with distributed meteogram
    ! file output!
    IF (num_io_procs > 0) THEN
      DO idom=1,max_dom
        IF ( meteogram_output_config(idom)%lenabled     .AND. &
          &  meteogram_output_config(idom)%ldistributed ) THEN
          CALL message(routine, 'Flag "ldistributed" collides with async IO. Resetting.')
          meteogram_output_config(idom)%ldistributed = .FALSE.
        END IF
      END DO
    END IF

  END SUBROUTINE check_meteogram_configuration

END MODULE mo_meteogram_config
