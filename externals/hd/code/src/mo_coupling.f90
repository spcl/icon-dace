! mo_coupling.f90 - Register HD to the coupler and access basic coupler functionality
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Author: Moritz Hanke (DKRZ)
! Contact: <stefan.hagemann@hereon.de>, <moritz.hanke@dkrz.de>
!_________________________________________
!----------------------------------
!>
!! Routines to access basic coupler functionality and routines for registering
!! the HD model to the coupler
!!
!! @par Revision History
!! First version by Moritz Hanke,  DKRZ, May 2022.
!!
!----------------------------

#ifdef NOMPI
#undef COUP_OAS
#undef COUP_YAC
#endif

#if defined(COUP_OAS) && defined(COUP_YAC)
#error "HD was configured with OASIS and YAC support. You can only use one of them."
#endif

MODULE mo_coupling

  USE mo_io_units, ONLY: nerr
  USE mo_kind

#ifdef COUP_OAS
  USE oas_hd_ini, only : kl_comm, oas_hd_init, ncomp_id, oas_hd_finalize
  USE mod_oasis ! OASIS3-MCT v3 module
#endif

#ifdef COUP_YAC
  USE mo_yac_finterface
  USE mpi
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_coupler
  PUBLIC :: finalize_coupler
  PUBLIC :: is_coupled_run

  PUBLIC :: set_grid_dimensions
  PUBLIC :: get_grid_dimensions
  PUBLIC :: get_grid_size
  PUBLIC :: set_local_partition
  PUBLIC :: get_local_partition

  PUBLIC :: set_coupling_type
  PUBLIC :: get_coupling_type

  INTEGER :: coupling_type   ! 0=no (def.)
                             ! 1=no interpolation
                             ! 2=interpolation in HD
  LOGICAL, PUBLIC :: lcoupling_atm  ! Switch for coupling to atmosphere  
  LOGICAL, PUBLIC :: lcoupling_oce  ! Switch for coupling to ocean (default: .False.)
  INTEGER, PUBLIC :: icpl_sinks     ! Redistribution of water in sinks for ocean coupling 
	                            ! 0=None (Def.), 1=to mouth boxes, 2=to all ocean boxes
  INTEGER, PUBLIC :: icpl_mask_tohd ! Switch for the mask on which HD may receive input for atmosphere coupling 
	                            ! 0=HD land boxes (Def.), 1=all HD boxes
 
  INTEGER,  ALLOCATABLE, PUBLIC :: fdir_hd(:,:) ! river directions with river mouth = 0
  REAL(dp), ALLOCATABLE, PUBLIC :: runoff_s(:,:)    ! surface water runoff [m/s]
  REAL(dp), ALLOCATABLE, PUBLIC :: runoff_dr(:,:)    ! soil water runoff [m/s]

  ! global grid information
  INTEGER :: nl = -1 ! number of longitudes
  INTEGER :: nb = -1 ! number of latitudes
  REAL(dp) :: florg ! north-west corner longitude of gridbox(1,1)
  REAL(dp) :: fborg ! north-west corner latitude of gridbox(1,1)
  REAL(dp) :: fscal ! grid resolution

  ! local decomposition information
  INTEGER :: local_start(2) = (/-1,-1/)
  INTEGER :: local_extent(2) = (/-1,-1/)

#ifdef COUP_YAC
  CHARACTER(LEN=*), PARAMETER :: yaml_filename = "coupling.yaml"

  LOGICAL :: config_file_have_been_checked = .FALSE.
  LOGICAL :: config_file_exist = .FALSE.
  LOGICAL :: yac_is_initialised = .FALSE.
#endif

  CHARACTER(*), PARAMETER :: modname = "mo_coupling"

CONTAINS

#ifdef COUP_YAC
  LOGICAL FUNCTION coupler_config_file_exist()

    IF (config_file_have_been_checked) THEN

      coupler_config_file_exist = config_file_exist

    ELSE

      INQUIRE(FILE=TRIM(ADJUSTL(yaml_filename)), EXIST=config_file_exist)

      config_file_have_been_checked = .TRUE.
      coupler_config_file_exist = config_file_exist

    END IF

  END FUNCTION
#endif

  SUBROUTINE init_coupler(world_communicator, global_name)

    INTEGER, INTENT(OUT) :: world_communicator
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: global_name

#ifdef COUP_OAS

    CALL oas_hd_init()
    world_communicator = kl_comm

#endif

#ifdef COUP_YAC

    CHARACTER(*), PARAMETER :: routine = modname//":init_coupler"

    INTEGER :: yac_comm, hd_comm, group_comms(2)
    CHARACTER(len=YAC_MAX_CHARLEN) :: group_names(2)
    INTEGER :: hd_group_rank
    INTEGER :: ierror

    CALL MPI_INIT (ierror)
    IF (ierror /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a)') ' MPI_INIT failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      STOP
    END IF

    IF (coupler_config_file_exist()) THEN

      yac_is_initialised = .TRUE.

      ! HD requires its own communicator, therefore we generate
      ! the YAC and HD communicator before initialising YAC
      group_names(1) = "yac"
      IF (PRESENT(global_name)) THEN
         group_names(2) = TRIM(global_name)
      ELSE
        group_names(2) = TRIM('hd')
      END IF
      CALL yac_fmpi_handshake( MPI_COMM_WORLD, group_names, group_comms)
      yac_comm = group_comms(1)
      hd_comm = group_comms(2)

      ! initialise YAC
      CALL yac_finit_comm(yac_comm)
      CALL MPI_COMM_FREE(yac_comm, ierror)

      ! read in configuration file at root rank (could also be
      ! done at every process)
      CALL MPI_COMM_RANK(hd_comm, hd_group_rank, ierror)
      IF ( hd_group_rank == 0 ) &
        CALL yac_fread_config_yaml(TRIM(yaml_filename))

    ELSE

      CALL mpi_comm_dup(MPI_COMM_WORLD, hd_comm, ierror)
      IF (ierror /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_COMM_DUP failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', ierror
        STOP
      END IF

    END IF

    world_communicator = hd_comm
#endif

  END SUBROUTINE init_coupler

  SUBROUTINE finalize_coupler
  INTEGER :: ierror

#ifdef COUP_OAS
   CALL oas_hd_finalize()
#endif

#ifdef COUP_YAC
    IF (yac_is_initialised) CALL yac_ffinalize

    CALL MPI_FINALIZE (ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', ierror
       STOP
    END IF

#endif

  END SUBROUTINE finalize_coupler

  SUBROUTINE print_info_stderr (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    WRITE (nerr,'(4a)') " ", TRIM(name), ": ", TRIM(text)

  END SUBROUTINE print_info_stderr

  LOGICAL FUNCTION is_coupled_run()

#if defined(COUP_OAS) || defined(COUP_YAC)
    is_coupled_run = .TRUE.
#else
    is_coupled_run = .FALSE.
#endif
  END FUNCTION is_coupled_run

  SUBROUTINE coupling_abort(routine, message)

    CHARACTER(*), INTENT(IN) :: routine
    CHARACTER(*), INTENT(IN) :: message

#ifdef COUP_OAS
    CALL oasis_abort( ncomp_id, routine, message )
#endif

#ifdef COUP_YAC
    CALL yac_abort_message( routine // " : " // message, __FILE__, __LINE__)
#endif

  END SUBROUTINE coupling_abort

  SUBROUTINE set_grid_dimensions(nlon, nlat, lon_org, lat_org, res)

    CHARACTER(*), PARAMETER :: routine = modname//":set_grid_dimensions"

    INTEGER, INTENT(IN) :: nlon
    INTEGER, INTENT(IN) :: nlat
    REAL(dp), INTENT(IN) :: lon_org ! north-west corner longitude of gridbox(1,1)
    REAL(dp), INTENT(IN) :: lat_org ! north-west corner latitude of gridbox(1,1)
    REAL(dp), INTENT(IN) :: res     ! grid resolution

    IF ((nlon < 1) .OR. (nlat < 1)) THEN
      CALL coupling_abort( routine, 'invalid grid dimensions' )
    END IF

    IF ((lat_org > 90.0) .OR. (lat_org <= -90.0)) THEN
      CALL coupling_abort( &
        routine, 'invalid grid dimensions (latitude of north-west corner' )
    END IF

    nl = nlon
    nb = nlat
    florg = lon_org
    fborg = lat_org
    fscal = res

  END SUBROUTINE set_grid_dimensions

  SUBROUTINE get_grid_dimensions(nlon, nlat, lon_org, lat_org, res)

    CHARACTER(*), PARAMETER :: routine = modname//":get_grid_dimensions"

    INTEGER, INTENT(OUT) :: nlon
    INTEGER, INTENT(OUT) :: nlat
    REAL(dp), INTENT(OUT) :: lon_org ! north-west corner longitude of gridbox(1,1)
    REAL(dp), INTENT(OUT) :: lat_org ! north-west corner latitude of gridbox(1,1)
    REAL(dp), INTENT(OUT) :: res     ! grid resolution

    IF ((nl == -1) .OR. (nb == -1)) THEN
      CALL coupling_abort(routine, 'grid dimensions have not yet been set')
    END IF

    nlon = nl
    nlat = nb
    lon_org = florg
    lat_org = fborg
    res = fscal

  END SUBROUTINE get_grid_dimensions

  SUBROUTINE get_grid_size(nlon, nlat)

    CHARACTER(*), PARAMETER :: routine = modname//":get_grid_size"

    INTEGER, INTENT(OUT) :: nlon ! number of longitudes
    INTEGER, INTENT(OUT) :: nlat ! number of latitudes

    IF ((nl == -1) .OR. (nb == -1)) THEN
      CALL coupling_abort(routine, 'grid dimensions have not yet been set')
    END IF

    nlon = nl
    nlat = nb

  END SUBROUTINE get_grid_size

  SUBROUTINE set_local_partition( &
    local_extent_lon, local_extent_lat, local_start_lon, local_start_lat)

    CHARACTER(*), PARAMETER :: routine = modname//":set_local_partition"

    INTEGER, INTENT(IN) :: local_extent_lon
    INTEGER, INTENT(IN) :: local_extent_lat
    INTEGER, INTENT(IN) :: local_start_lon
    INTEGER, INTENT(IN) :: local_start_lat

    IF ((nl == -1) .OR. (nb == -1)) THEN
      CALL coupling_abort(routine, 'grid dimensions have not yet been set')
    END IF

    IF ((local_extent_lon > nl) .OR. &
        (local_extent_lat > nb) .OR. &
        (local_start_lon < 1) .OR. &
        (local_start_lon > nl) .OR. &
        (local_start_lat < 1) .OR. &
        (local_start_lat > nb) .OR. &
        (local_extent_lon + local_start_lon - 1 > nl) .OR. &
        (local_extent_lat + local_start_lat - 1 > nb)) THEN
      CALL coupling_abort(routine, 'invalid local partition information')
    END IF

    local_start(1) = local_start_lon
    local_start(2) = local_start_lat
    local_extent(1) = local_extent_lon
    local_extent(2) = local_extent_lat

  END SUBROUTINE set_local_partition

  SUBROUTINE get_local_partition( &
    local_extent_lon, local_extent_lat, local_start_lon, local_start_lat)

    CHARACTER(*), PARAMETER :: routine = modname//":get_local_partition"

    INTEGER, INTENT(OUT) :: local_extent_lon
    INTEGER, INTENT(OUT) :: local_extent_lat
    INTEGER, INTENT(OUT) :: local_start_lon
    INTEGER, INTENT(OUT) :: local_start_lat

    IF ((local_start(1) == -1) .OR. (local_start(2) == -1) .OR. &
        (local_extent(1) == -1) .OR. (local_extent(2) == -1)) THEN
      CALL coupling_abort(routine, 'local partition information has not yet been set')
    END IF

    local_start_lon = local_start(1)
    local_start_lat = local_start(2)
    local_extent_lon = local_extent(1)
    local_extent_lat = local_extent(2)

  END SUBROUTINE get_local_partition

  SUBROUTINE set_coupling_type(coupling_type_)

    CHARACTER(*), PARAMETER :: routine = modname//":set_coupling_type"

    INTEGER, INTENT(IN) :: coupling_type_

    IF ((coupling_type_ < 0) .OR. (coupling_type_ > 2)) THEN
      CALL coupling_abort(routine, 'invalid value for coupling_type')
    END IF

    coupling_type = coupling_type_

  END SUBROUTINE set_coupling_type

  FUNCTION get_coupling_type()

    INTEGER :: get_coupling_type

    get_coupling_type = coupling_type

  END FUNCTION get_coupling_type

END MODULE mo_coupling
