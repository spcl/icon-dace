! mo_mpi.f90 - Settings/Routines for using MPI or OpenMPI 
! 
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!
! This file has been modified for the use in HD-model.
!_________________________________________


MODULE mo_mpi
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019).
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11,
  !  doi: 10.1029/2018MS001400.
  !
  !  - April 2018, Ha Ho-Hagemann, HZG: Modifications for coupling via OASIS
  !
  !
  ! Comment: Please use basic WRITE to nerr for messaging in the whole
  !          MPI package to achieve promer output.

  ! actual method (MPI-2)
#ifndef NOMPI
  USE mpi
#endif

#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num, &
                     omp_get_num_threads, omp_set_num_threads
#endif

  USE mo_coupling, ONLY: is_coupled_run, init_coupler, finalize_coupler
  USE mo_kind
  USE mo_io_units,    ONLY: nerr

  IMPLICIT NONE

  PRIVATE                          ! all declarations are private

  TYPE p_io_server_info
    INTEGER :: comp_source
    INTEGER :: io_target
  END TYPE p_io_server_info

  ! subroutines defined, overloaded depending on argument type

  PUBLIC :: p_start, p_stop, p_abort
  PUBLIC :: p_send, p_recv, p_sendrecv, p_bcast, p_barrier
  PUBLIC :: p_isend, p_irecv, p_wait, p_wait_any
  PUBLIC :: p_gather, p_max, p_min, p_sum, p_global_sum, p_field_sum
  PUBLIC :: p_set_communicator, p_init_communicators, p_is_iope
  PUBLIC :: p_probe
  PUBLIC :: p_inc_request
#ifndef NOMPI
  PUBLIC :: p_request, p_irequest
#endif

  ! data type matching Fortran and MPI implementation

  PUBLIC :: p_real
  PUBLIC :: p_real_sp
  PUBLIC :: p_real_dp
  PUBLIC :: p_int
  PUBLIC :: p_int_i4
  PUBLIC :: p_int_i8

  PUBLIC :: p_address_kind

  PUBLIC :: p_any_source

  ! logical switches

  PUBLIC :: p_parallel, p_parallel_io, p_ioservers

  ! PE identifier

  PUBLIC :: p_pe, p_io, p_nprocs
  PUBLIC :: p_iopes, p_gio, p_num_iopes
  PUBLIC :: p_io_server_info, p_io_info_pair

  ! communicator

  PUBLIC :: p_global_comm
  PUBLIC :: p_all_comm
  PUBLIC :: p_ioserver_comm
  PUBLIC :: p_io_comm
  PUBLIC :: p_communicator_a, p_communicator_b, p_communicator_d

  PUBLIC :: p_info_null

  ! old fashioned method (MPI-1)

!!$#ifndef NOMPI
!!$  INCLUDE 'mpif.h'
!!$#endif

  ! general run time information

#ifndef NOMPI
  INTEGER :: version, subversion   ! MPI version
#endif

  ! MPI call inherent variables

#ifndef NOMPI
  INTEGER :: p_error                     ! MPI error number
  INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV
  PUBLIC :: p_error
#endif

  ! public parallel run information

  LOGICAL :: p_parallel    = .FALSE.
  LOGICAL :: p_parallel_io = .TRUE.
  LOGICAL :: p_ioservers   = .FALSE.

  INTEGER :: p_pe        = 0     ! this is the PE number of this task
  INTEGER :: p_io        = 0     ! PE number of PE handling IO
  INTEGER :: p_nprocs    = 1     ! number of available PEs (processors)
  INTEGER :: p_num_iopes = 0     ! number of I/O server PEs

  INTEGER, ALLOCATABLE :: p_iopes(:) ! I/O server PEs internal counting
  INTEGER, ALLOCATABLE :: p_gio(:)   ! I/O server PE mapping from global ranks
  ! I/O server PEs mapping to compute PEs for information handling
  TYPE(p_io_server_info), ALLOCATABLE :: p_io_info_pair(:)

  ! communicator sets

  INTEGER :: p_global_comm             ! replaces MPI_COMM_WORLD in one application
  INTEGER :: p_all_comm                ! communicator of application without I/O servers
  INTEGER :: p_ioserver_comm           ! communicator of application for I/O servers only
  INTEGER, ALLOCATABLE :: p_io_comm(:) ! communicator of application to I/O servers
  INTEGER :: p_communicator_a          ! for Set A
  INTEGER :: p_communicator_b          ! for Set B
  INTEGER :: p_communicator_d          ! for debug node

  ! Info

#ifndef NOMPI
  INTEGER :: p_info_null = MPI_INFO_NULL
#else
  INTEGER :: p_info_null = 0
#endif
  ! non blocking calls

#ifndef NOMPI
  INTEGER, ALLOCATABLE, SAVE :: p_request(:) ! request values for non blocking calls
  INTEGER :: p_irequest ! the first p_irequest values of p_request are in use
  INTEGER :: p_mrequest ! actual size of p_request
  INTEGER, PARAMETER :: p_request_alloc_size = 4096
#endif

  ! module intrinsic names

  INTEGER :: mype                  ! this is the PE number of this task
!!#ifndef NOMPI
!!LK  INTEGER :: iope                  ! PE able to do IO
!!#endif
  INTEGER :: npes                  ! number of available PEs

  INTEGER :: nbcast = 0            ! counter for broadcasts for debugging

  ! MPI native transfer types

  INTEGER :: p_int     = 0
  INTEGER :: p_real    = 0
  INTEGER :: p_bool    = 0
  INTEGER :: p_char    = 0

  ! MPI transfer types

  INTEGER :: p_real_dp = 0
  INTEGER :: p_real_sp = 0
  INTEGER :: p_int_i4  = 0
  INTEGER :: p_int_i8  = 0

  ! MPI native transfer types byte size

  INTEGER :: p_int_byte     = 0
  INTEGER :: p_real_byte    = 0
!  INTEGER :: p_bool_byte    = 0
!  INTEGER :: p_char_byte    = 0

  ! MPI transfer types byte size

  INTEGER :: p_real_dp_byte = 0
  INTEGER :: p_real_sp_byte = 0
  INTEGER :: p_int_i4_byte  = 0
  INTEGER :: p_int_i8_byte  = 0

#ifndef NOMPI
  INTEGER, PARAMETER :: p_address_kind = MPI_ADDRESS_KIND
  INTEGER            :: p_any_source   = MPI_ANY_SOURCE
#else
  INTEGER, PARAMETER :: p_address_kind = i8    ! should not get touched at all
  INTEGER            :: p_any_source = -1
#endif

  ! define generic interfaces to allow proper compiling with picky compilers
  ! like NAG f95 for clean argument checking and shortening the call sequence.

  INTERFACE p_send
     MODULE PROCEDURE p_send_real
     MODULE PROCEDURE p_send_int
     MODULE PROCEDURE p_send_bool
     MODULE PROCEDURE p_send_real_1d
     MODULE PROCEDURE p_send_int_1d
     MODULE PROCEDURE p_send_bool_1d
     MODULE PROCEDURE p_send_real_2d
     MODULE PROCEDURE p_send_int_2d
     MODULE PROCEDURE p_send_bool_2d
     MODULE PROCEDURE p_send_real_3d
     MODULE PROCEDURE p_send_int_3d
     MODULE PROCEDURE p_send_bool_3d
     MODULE PROCEDURE p_send_real_4d
     MODULE PROCEDURE p_send_int_4d
     MODULE PROCEDURE p_send_bool_4d
     MODULE PROCEDURE p_send_char
     MODULE PROCEDURE p_send_real_5d
  END INTERFACE

  INTERFACE p_isend
     MODULE PROCEDURE p_isend_real
     MODULE PROCEDURE p_isend_int
     MODULE PROCEDURE p_isend_bool
     MODULE PROCEDURE p_isend_real_1d
     MODULE PROCEDURE p_isend_int_1d
     MODULE PROCEDURE p_isend_bool_1d
     MODULE PROCEDURE p_isend_real_2d
     MODULE PROCEDURE p_isend_int_2d
     MODULE PROCEDURE p_isend_bool_2d
     MODULE PROCEDURE p_isend_real_3d
     MODULE PROCEDURE p_isend_int_3d
     MODULE PROCEDURE p_isend_bool_3d
     MODULE PROCEDURE p_isend_real_4d
     MODULE PROCEDURE p_isend_int_4d
     MODULE PROCEDURE p_isend_bool_4d
     MODULE PROCEDURE p_isend_char
     MODULE PROCEDURE p_isend_real_5d
  END INTERFACE

  INTERFACE p_recv
     MODULE PROCEDURE p_recv_real
     MODULE PROCEDURE p_recv_int
     MODULE PROCEDURE p_recv_bool
     MODULE PROCEDURE p_recv_real_1d
     MODULE PROCEDURE p_recv_int_1d
     MODULE PROCEDURE p_recv_bool_1d
     MODULE PROCEDURE p_recv_real_2d
     MODULE PROCEDURE p_recv_int_2d
     MODULE PROCEDURE p_recv_bool_2d
     MODULE PROCEDURE p_recv_real_3d
     MODULE PROCEDURE p_recv_int_3d
     MODULE PROCEDURE p_recv_bool_3d
     MODULE PROCEDURE p_recv_real_4d
     MODULE PROCEDURE p_recv_int_4d
     MODULE PROCEDURE p_recv_bool_4d
     MODULE PROCEDURE p_recv_char
     MODULE PROCEDURE p_recv_real_5d
  END INTERFACE

  INTERFACE p_irecv
     MODULE PROCEDURE p_irecv_real
     MODULE PROCEDURE p_irecv_real_1d
     MODULE PROCEDURE p_irecv_real_2d
     MODULE PROCEDURE p_irecv_real_3d
     MODULE PROCEDURE p_irecv_real_4d
  END INTERFACE

  INTERFACE p_sendrecv
     MODULE PROCEDURE p_sendrecv_real_1d
     MODULE PROCEDURE p_sendrecv_real_2d
     MODULE PROCEDURE p_sendrecv_real_3d
     MODULE PROCEDURE p_sendrecv_real_4d
  END INTERFACE

  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_real
     MODULE PROCEDURE p_bcast_int_i4
     MODULE PROCEDURE p_bcast_int_i8
     MODULE PROCEDURE p_bcast_bool
     MODULE PROCEDURE p_bcast_real_1d
     MODULE PROCEDURE p_bcast_int_1d
     MODULE PROCEDURE p_bcast_bool_1d
     MODULE PROCEDURE p_bcast_real_2d
     MODULE PROCEDURE p_bcast_int_2d
     MODULE PROCEDURE p_bcast_bool_2d
     MODULE PROCEDURE p_bcast_real_3d
     MODULE PROCEDURE p_bcast_int_3d
     MODULE PROCEDURE p_bcast_bool_3d
     MODULE PROCEDURE p_bcast_real_4d
     MODULE PROCEDURE p_bcast_int_4d
     MODULE PROCEDURE p_bcast_bool_4d
     MODULE PROCEDURE p_bcast_real_5d
     MODULE PROCEDURE p_bcast_char
     MODULE PROCEDURE p_bcast_char_1d
  END INTERFACE

  INTERFACE p_gather
     MODULE PROCEDURE p_gather_real_1d2d
     MODULE PROCEDURE p_gather_real_5d6d
  END INTERFACE

  INTERFACE p_max
     MODULE PROCEDURE p_max_0d
     MODULE PROCEDURE p_max_1d
     MODULE PROCEDURE p_max_2d
     MODULE PROCEDURE p_max_3d
  END INTERFACE

  INTERFACE p_min
     MODULE PROCEDURE p_min_0d
     MODULE PROCEDURE p_min_1d
     MODULE PROCEDURE p_min_2d
     MODULE PROCEDURE p_min_3d
  END INTERFACE

  INTERFACE p_sum
     MODULE PROCEDURE p_sum_0d
     MODULE PROCEDURE p_sum_1d
  END INTERFACE

  INTERFACE p_global_sum
     MODULE PROCEDURE p_global_sum_1d
  END INTERFACE

  INTERFACE p_field_sum
     MODULE PROCEDURE p_field_sum_1d
  END INTERFACE

CONTAINS

  SUBROUTINE p_start(model_name)

    USE mo_util_string, ONLY: toupper

    CHARACTER(len=*), INTENT(in) :: model_name

    ! variables are required for determing I/O size in bytes of the defined
    ! KIND types for assigning the right MPI data types with the used kinds

    INTEGER     :: iig = 0
    INTEGER(i4) :: ii4 = 0_i4
    INTEGER(i8) :: ii8 = 0_i8
    REAL        :: rrg = 0.0
    REAL(sp)    :: rsp = 0.0_sp
    REAL(dp)    :: rdp = 0.0_dp
!    LOGICAL     :: llg = .FALSE.
!    CHARACTER   :: yyg = 'A'

    ! variables used for determing the OpenMP threads
    ! suitable as well for coupled models


    CHARACTER(len=32) :: env_name
    CHARACTER(len=32) :: thread_num
    INTEGER :: env_threads
#if defined _OPENMP
    INTEGER ::  thread
    INTEGER, ALLOCATABLE :: threads(:)
#endif

    ! status variable
    INTEGER :: istat
#if defined _OPENMP && !defined NOMPI
    ! loop index
    INTEGER :: jp
#endif

#ifndef NOMPI
    ! start MPI
    IF (is_coupled_run()) THEN

      ! init MPI and coupler and replace p_global_comm with
      ! HD-global-communicator, if HD is coupled to other components
      CALL init_coupler(p_global_comm, model_name)

    ELSE

      CALL MPI_INIT (p_error)
      IF (p_error /= MPI_SUCCESS) THEN
         WRITE (nerr,'(a)') ' MPI_INIT failed.'
         WRITE (nerr,'(a,i4)') ' Error =  ', p_error
         STOP
      END IF
      CALL MPI_COMM_DUP (MPI_COMM_WORLD, p_global_comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
         WRITE (nerr,'(a)') ' MPI_COMM_DUP failed.'
         WRITE (nerr,'(a,i4)') ' Error =  ', p_error
         CALL p_abort
      END IF

    END IF

    CALL MPI_ERRHANDLER_SET(p_global_comm, MPI_ERRORS_RETURN, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ERRHANDLER_SET failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF

    ! get local PE identification
    CALL MPI_COMM_RANK (p_global_comm, mype, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! get number of available PEs
    CALL MPI_COMM_SIZE (p_global_comm, npes, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! if HD was build without MPI
#else
    mype = 0
    npes = 1
#endif

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a)') ' PE ', mype, ' started.'
#endif

    ! for non blocking calls

#ifndef NOMPI
    p_mrequest = p_request_alloc_size
    ALLOCATE(p_request(p_mrequest))
    p_irequest = 0
#endif

! set the I/O pe statically to 0

    p_io = 0

    ! Information ...

    IF (mype == 0) THEN
       WRITE (nerr,'(/,a,a,a)') ' ', &
            TRIM(model_name), ' MPI interface runtime information:'
    END IF

#ifndef NOMPI
    ! lets check the available MPI version

    CALL MPI_GET_VERSION (version, subversion, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a)') ' MPI_GET_VERSION failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      CALL p_abort
    END IF

    IF (mype == 0) THEN
      WRITE (nerr,'(a,i0,a1,i0)') &
           '  Used MPI version: ', version, '.', subversion
    END IF
#endif

    IF (npes < 2) THEN
      p_parallel = .FALSE.
      p_parallel_io = .TRUE.   ! can always do I/O
      IF (mype == 0) THEN
        WRITE (nerr,'(a)') '  Single processor run.'
      END IF
      p_pe = 0
      p_nprocs = 1
    ELSE
      p_parallel = .TRUE.
      IF (mype == p_io) THEN
        p_parallel_io = .TRUE.
      ELSE
        p_parallel_io = .FALSE.
      END IF
      IF (mype == 0) THEN
        WRITE (nerr,'(a,i0,a)') '  Run on ', npes, ' processors.'
      END IF
      p_pe = mype
      p_nprocs = npes
    END IF

    ! Have to do the retrieval of the environment variable outside of the
    ! OpenMP accessible section to allow the MPI broadcast for the testing
    ! cases.

    IF (mype == 0) THEN

print*,'MOMPI: ',model_name
      env_name = TRIM(toupper(model_name)) // '_THREADS'
print*,'MOMPI: ', env_name
      CALL get_environment_variable(name=TRIM(env_name), value=thread_num, &
           status=istat)

      IF (istat == 0) THEN
        READ(thread_num,*) env_threads
      ELSE
        WRITE (nerr,'(1x,a)') ' Ha Ho-Hagemann: As there is no sub. get_environment_variable in src/:'
        WRITE (nerr,'(1x,a)') ' Number of OpenMP threads not given!'
        WRITE (nerr,'(1x,a,a,a,/,1x,a,a,a)') &
	     ' Environment variable ', TRIM(env_name), ' either not set,', &
             ' or not available to ', TRIM(model_name), ' root PE. Set to 1!'
        env_threads = 1
      ENDIF

   ENDIF

#ifndef NOMPI
    ! Make number of threads from environment available to all model PEs

    CALL MPI_BCAST (env_threads, 1, MPI_INTEGER, 0, p_global_comm, p_error)
#endif

#ifdef _OPENMP

    ! Expect that PE 0 did got the information of OMP_NUM_THREADS.
    ! That might be wrong in the coupled case when the model is
    ! started via MPI dynamic process creation. So we have to check
    ! the environment variable too.

    IF (mype == 0) THEN
      IF (p_parallel) THEN
        WRITE (nerr,'(/,a)') ' Running hybrid OpenMP-MPI mode.'
      ELSE
        WRITE (nerr,'(/,a)') ' Running OpenMP mode.'
      ENDIF
    ENDIF

    ! Inform on OpenMP thread usage

    CALL OMP_SET_DYNAMIC(.FALSE.)
    CALL OMP_SET_NUM_THREADS(env_threads)

    ! Check information for each team thread memeber

    ALLOCATE(threads(0:env_threads-1))

!$OMP PARALLEL PRIVATE(thread)
    thread = OMP_GET_THREAD_NUM()
    threads(thread) = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

    IF (ANY(threads /= threads(0))) THEN
      WRITE (nerr,'(a)') ' Setting of number of team threads failed.'
      CALL p_abort
    ENDIF

    IF (mype == 0) THEN

       WRITE (nerr,*)

       ! write out thread usage of master PE

       WRITE (nerr,'(1x,a,i2,a)') ' Use ', env_threads, ' OpenMP threads.'

    ENDIF
#endif

#ifndef NOMPI

    ! due to a possible circular dependency with mo_machine and other
    ! modules, we determine here locally the I/O size of the different
    ! kind types (assume 8 bit/byte. This is than used for determing
    ! the right MPI send/receive type parameters.

    CALL MPI_SIZEOF(iig, p_int_byte, p_error)
    CALL MPI_SIZEOF(ii4, p_int_i4_byte, p_error)
    CALL MPI_SIZEOF(ii8, p_int_i8_byte, p_error)
    CALL MPI_SIZEOF(rrg, p_real_byte, p_error)
    CALL MPI_SIZEOF(rsp, p_real_sp_byte, p_error)
    CALL MPI_SIZEOF(rdp, p_real_dp_byte, p_error)
!!$       CALL MPI_SIZEOF(llg, p_bool_byte, p_error)
!!$       CALL MPI_SIZEOF(yyg, p_char_byte, p_error)

    p_int     = MPI_INTEGER
    p_real    = MPI_REAL
    p_bool    = MPI_LOGICAL
    p_char    = MPI_CHARACTER

    CALL MPI_TYPE_CREATE_F90_REAL(ps, rs, p_real_sp, p_error)
    CALL MPI_TYPE_CREATE_F90_REAL(pd, rd, p_real_dp, p_error)

    CALL MPI_TYPE_CREATE_F90_INTEGER(pi4, p_int_i4, p_error)
    CALL MPI_TYPE_CREATE_F90_INTEGER(pi8, p_int_i8, p_error)

    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, p_real_sp_byte, p_real_sp, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, p_real_dp_byte, p_real_dp, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, p_int_i4_byte, p_int_i4, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, p_int_i8_byte, p_int_i8, p_error)

#endif

#ifdef DEBUG
    WRITE (nerr,'(/,a)')    ' MPI transfer sizes [bytes]:'
    WRITE (nerr,'(a,i4)') '  INTEGER generic:', p_int_byte
    WRITE (nerr,'(a,i4)') '  INTEGER 4 byte :', p_int_i4_byte
    WRITE (nerr,'(a,i4)') '  INTEGER 8 byte :', p_int_i8_byte
    WRITE (nerr,'(a,i4)') '  REAL generic   :', p_real_byte
    WRITE (nerr,'(a,i4)') '  REAL single    :', p_real_sp_byte
    WRITE (nerr,'(a,i4)') '  REAL double    :', p_real_dp_byte
#endif

    IF (mype == 0) THEN
      WRITE (nerr,*)
    ENDIF

  END SUBROUTINE p_start

  SUBROUTINE p_stop

    ! finish MPI and clean up all PEs
#ifndef NOMPI
    CALL MPI_COMM_FREE(p_global_comm, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_FREE failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF

    IF (is_coupled_run()) THEN

      CALL finalize_coupler()

    ELSE

     CALL MPI_FINALIZE (p_error)
     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
     END IF

    END IF

    p_parallel = .FALSE.
    DEALLOCATE(p_request)
#endif

  END SUBROUTINE p_stop

  SUBROUTINE p_abort

    ! this routine should be used instead of abort, util_abort() or STOP
    ! in all routines for proper clean up of all PEs

!OSBUTL    EXTERNAL :: util_exit

#ifndef NOMPI
    CALL MPI_ABORT (MPI_COMM_WORLD, 1, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#else
    WRITE (nerr,'(a)') 'mo_mpi: p_abort ...'
!OSBUTL    CALL util_exit(1)
    STOP 'util_exit'
#endif

  END SUBROUTINE p_abort

  ! Create computation-only communicator and I/O communicators

  SUBROUTINE p_init_communicators(nproca, nprocb, nprocio)

    INTEGER, INTENT(in) :: nproca, nprocb, nprocio

#ifndef NOMPI
    INTEGER :: jp, jio, group_world, grpid
    INTEGER, ALLOCATABLE :: pe_ranks(:)
#endif

    p_num_iopes = nprocio
    IF (p_num_iopes < 1) p_num_iopes = 0

    IF (p_num_iopes > 0) p_ioservers = .TRUE.

    IF (p_num_iopes > nproca*nprocb) THEN
        WRITE (nerr,'(a)') &
             'Cannot use more I/O servers than compute PEs.'
        CALL p_abort
    ENDIF

    IF (nproca*nprocb+nprocio > p_nprocs) THEN
        WRITE (nerr,'(a)') &
             'Not enough processors available for selected decomposition.'
        CALL p_abort
    ENDIF

#ifndef NOMPI
    IF (p_num_iopes == 0) THEN

      CALL MPI_COMM_DUP (p_global_comm, p_all_comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_COMM_DUP failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

    ELSE

      CALL MPI_COMM_GROUP(p_global_comm, group_world, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_COMM_GROUP failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

      ! build compute PEs communicator, adjust number of compute PEs
      ! so that we do not need to change the models code
      p_nprocs = p_nprocs - p_num_iopes
      ALLOCATE(pe_ranks(0:p_nprocs-1))
      pe_ranks(0:p_nprocs-1) = (/ (jp, jp = 0, p_nprocs-1) /)

      CALL MPI_GROUP_INCL(group_world, p_nprocs, pe_ranks, grpid, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_GROUP_INCL failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

      p_all_comm = MPI_COMM_NULL

      CALL MPI_COMM_CREATE(p_global_comm, grpid, p_all_comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_COMM_CREATE failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

      DEALLOCATE(pe_ranks)

      CALL MPI_GROUP_FREE(grpid, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_GROUP_FREE failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

      ! do we need an I/O PEs only communicator?  yes

      ALLOCATE(pe_ranks(0:p_num_iopes-1))
      pe_ranks(0:p_num_iopes-1) = (/ (jp, jp = p_nprocs, p_nprocs+p_num_iopes-1) /)

      CALL MPI_GROUP_INCL(group_world, p_num_iopes, pe_ranks, grpid, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_COMM_INCL failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

      p_ioserver_comm = MPI_COMM_NULL

      CALL MPI_COMM_CREATE(p_global_comm, grpid, p_ioserver_comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_COMM_CREATE failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

      DEALLOCATE(pe_ranks)

      CALL MPI_GROUP_FREE(grpid, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a)') ' MPI_GROUP_FREE failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL p_abort
      END IF

      ! build a communicator between all compute PEs and each I/O PE

      ALLOCATE(p_io_comm(0:p_num_iopes-1))
      p_io_comm(:) = MPI_COMM_NULL
      ALLOCATE(pe_ranks(0:p_nprocs))
      pe_ranks(0:p_nprocs-1) = (/ (jp, jp = 0, p_nprocs-1) /)
      jio = 0
      DO jp = p_nprocs, p_nprocs+p_num_iopes-1
        pe_ranks(p_nprocs) = jp

        CALL MPI_GROUP_INCL(group_world, p_nprocs+1, pe_ranks, grpid, p_error)
        IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
        END IF

        CALL MPI_COMM_CREATE(p_global_comm, grpid, p_io_comm(jio), p_error)
        IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
        END IF

        CALL MPI_GROUP_FREE(grpid, p_error)
        IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_GROUP_FREE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
        END IF

        jio = jio+1
      ENDDO
      DEALLOCATE(pe_ranks)

   ALLOCATE(p_iopes(0:p_num_iopes-1))
    p_iopes(0:p_num_iopes-1) = (/ (jp, jp = p_nprocs, p_nprocs+p_num_iopes-1) /)

    ALLOCATE(p_gio(p_nprocs:p_nprocs+p_num_iopes-1))
    p_gio(p_nprocs:p_nprocs+p_num_iopes-1) = (/ (jp, jp = 0, p_num_iopes-1) /)

    ALLOCATE(p_io_info_pair(0:p_nprocs-1))

    DO jp = 0, p_num_iopes-1
      p_io_info_pair(jp)%comp_source = jp
    ENDDO
    IF (p_nprocs >= p_num_iopes) THEN
      DO jp = p_num_iopes, p_nprocs-1
        p_io_info_pair(jp)%comp_source = MPI_PROC_NULL
      ENDDO
    ENDIF

    DO jp = 0, p_num_iopes-1
      p_io_info_pair(jp)%io_target   = p_iopes(jp)
    ENDDO
    IF (p_nprocs >= p_num_iopes) THEN
      DO jp = p_num_iopes, p_nprocs-1
        p_io_info_pair(jp)%io_target = MPI_PROC_NULL
      ENDDO
    ENDIF
  ENDIF

#ifdef DEBUG
  ! diagnostics of setup:

  WRITE (0,*) 'All parallel PE mappings:'
  WRITE (0,*) 'Number of total PEs:      ', nproca*nprocb+nprocio
  WRITE (0,*) 'Number of compute PEs:    ', p_nprocs
  WRITE (0,*) 'Number of I/O server PEs: ', p_num_iopes
  WRITE (0,*) 'Ranks of compute PEs in p_global_comm:'
  ALLOCATE(pe_ranks(0:p_nprocs-1))
  pe_ranks(0:p_nprocs-1) = (/ (jp, jp = 0, p_nprocs-1) /)
  WRITE (0,*) pe_ranks
  DEALLOCATE(pe_ranks)
  WRITE (0,*) 'Ranks of I/O server PEs in p_global_comm:'
  WRITE (0,*) p_iopes
  WRITE (0,*) 'Check global to local mapping for I/O server:'
  DO jp = 0, p_num_iopes-1
    WRITE (0,*) '    global index: ', p_iopes(jp), ' local index: ', p_gio(p_iopes(jp))
  ENDDO
  WRITE (0,*) 'Compute - I/O server pair:'
  DO jp = 0, SIZE(p_io_info_pair)-1
    WRITE (0,*) &
         '   PE: ',    jp,                             &
         ' source: ',  p_io_info_pair(jp)%comp_source, &
         ' target: ',  p_io_info_pair(jp)%io_target
  ENDDO
#endif

#endif

  END SUBROUTINE p_init_communicators

  ! communicator set up

  SUBROUTINE p_set_communicator (nproca, nprocb, mapmesh, debug_parallel)

    INTEGER, INTENT(in) :: nproca, nprocb
    INTEGER, INTENT(in) :: mapmesh(0:,0:)
    INTEGER, INTENT(in) :: debug_parallel

#ifndef NOMPI
    INTEGER :: all_debug_pes(SIZE(mapmesh))

    INTEGER :: group_world, group_a, group_b, group_d
    INTEGER :: p_communicator_tmp

    INTEGER :: n, members

    INTEGER :: ranks(1) = 0

    ! first set global group

    CALL MPI_COMM_GROUP (p_all_comm, group_world, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_GROUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! communicator is p_all_comm

    IF (debug_parallel >= 0 ) THEN

       CALL MPI_GROUP_INCL (group_world, 1, ranks, group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype == 0) p_communicator_d = p_communicator_tmp

       DO n = 1, SIZE(mapmesh)
          all_debug_pes(n) = n
       END DO

       CALL MPI_GROUP_INCL (group_world, SIZE(mapmesh), all_debug_pes, &
            group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype /= 0) p_communicator_d = p_communicator_tmp

    ELSE
       p_communicator_d = p_all_comm
    END IF

    DO n = 0, nproca-1
       members = nprocb
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(:,n), group_a, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_a, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_a = p_communicator_tmp

    END DO

    ! create groups for set Bs

    DO n = 0, nprocb-1
       members = nproca
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(n,:), group_b, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_b, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_b = p_communicator_tmp

    END DO

    CALL MPI_BARRIER (p_all_comm, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_BARRIER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    IF (debug_parallel >= 0 .AND. mype == 0) THEN
      p_communicator_a = p_communicator_d
      p_communicator_b = p_communicator_d
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,3i8)') &
         'p_set_communicator on PE ', mype, ': ', &
         p_communicator_d, &
         p_communicator_a, &
         p_communicator_b
#endif
#endif
  END SUBROUTINE p_set_communicator

  FUNCTION p_is_iope() RESULT(is_iope)
    ! .TRUE. if dedicated I/O PE
    logical :: is_iope
    is_iope = (p_pe >= p_nprocs)
  END FUNCTION p_is_iope

!=========================================================================

  ! send implementation

  SUBROUTINE p_send_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real

  SUBROUTINE p_send_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d

  SUBROUTINE p_send_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_2d

  SUBROUTINE p_send_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_3d

  SUBROUTINE p_send_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_4d

  SUBROUTINE p_send_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_5d

  SUBROUTINE p_send_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_1d

  SUBROUTINE p_send_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_2d

  SUBROUTINE p_send_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_3d

  SUBROUTINE p_send_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_4d


  SUBROUTINE p_send_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool

  SUBROUTINE p_send_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_1d

  SUBROUTINE p_send_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_2d

  SUBROUTINE p_send_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_3d

  SUBROUTINE p_send_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_4d

  SUBROUTINE p_send_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in) :: buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, LEN(buffer), p_char, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_char

! non-blocking sends

  SUBROUTINE p_inc_request
    INTEGER, ALLOCATABLE :: tmp(:)

#ifndef NOMPI
    p_irequest = p_irequest + 1
    IF (p_irequest > p_mrequest) THEN
      ALLOCATE(tmp(p_mrequest))
      tmp(:) = p_request(:)
      DEALLOCATE(p_request)
      ALLOCATE(p_request(p_mrequest+p_request_alloc_size))
      p_request(1:p_mrequest) = tmp(:)
      p_mrequest = p_mrequest+p_request_alloc_size
      DEALLOCATE(tmp)
    ENDIF
#endif

  END SUBROUTINE p_inc_request

  SUBROUTINE p_isend_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real

  SUBROUTINE p_isend_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d

  SUBROUTINE p_isend_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_2d

  SUBROUTINE p_isend_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_3d

  SUBROUTINE p_isend_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_4d

  SUBROUTINE p_isend_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_5d

  SUBROUTINE p_isend_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int

  SUBROUTINE p_isend_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_1d

  SUBROUTINE p_isend_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_2d

  SUBROUTINE p_isend_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_3d

  SUBROUTINE p_isend_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_4d


  SUBROUTINE p_isend_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool

  SUBROUTINE p_isend_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_1d

  SUBROUTINE p_isend_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_2d

  SUBROUTINE p_isend_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_3d

  SUBROUTINE p_isend_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_4d

  SUBROUTINE p_isend_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(inout) :: buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, LEN(buffer), p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_char

  ! recv implementation

  SUBROUTINE p_recv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real

  SUBROUTINE p_recv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d

  SUBROUTINE p_recv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_2d

  SUBROUTINE p_recv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_3d

  SUBROUTINE p_recv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_4d

  SUBROUTINE p_recv_real_5d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_5d

  SUBROUTINE p_recv_int (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_int_1d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_1d

  SUBROUTINE p_recv_int_2d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_2d

  SUBROUTINE p_recv_int_3d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_3d

  SUBROUTINE p_recv_int_4d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_4d


  SUBROUTINE p_recv_bool (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool

  SUBROUTINE p_recv_bool_1d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_1d

  SUBROUTINE p_recv_bool_2d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_2d

  SUBROUTINE p_recv_bool_3d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_3d

  SUBROUTINE p_recv_bool_4d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_4d

  SUBROUTINE p_recv_char (buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(out) :: buffer
    INTEGER,           INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, LEN(buffer), p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_char

! non-blocking receives

  SUBROUTINE p_irecv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real

  SUBROUTINE p_irecv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d

  SUBROUTINE p_irecv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_2d

  SUBROUTINE p_irecv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_3d

  SUBROUTINE p_irecv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_4d

  ! sendrecv implementation

  SUBROUTINE p_sendrecv_real_1d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_1d

  SUBROUTINE p_sendrecv_real_2d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_2d

  SUBROUTINE p_sendrecv_real_3d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_3d

  SUBROUTINE p_sendrecv_real_4d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_4d

  ! bcast implementation

  SUBROUTINE p_bcast_real (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real

  SUBROUTINE p_bcast_real_1d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_1d

  SUBROUTINE p_bcast_real_2d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_2d

  SUBROUTINE p_bcast_real_3d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_3d

  SUBROUTINE p_bcast_real_4d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_4d

  SUBROUTINE p_bcast_real_5d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_5d

  SUBROUTINE p_bcast_int_i4 (buffer, p_source, comm)

    INTEGER (i4), INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_int_i4, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i4

  SUBROUTINE p_bcast_int_i8 (buffer, p_source, comm)

    INTEGER (i8), INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_int_i8, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8

  SUBROUTINE p_bcast_int_1d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_bcast_int_2d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_2d

  SUBROUTINE p_bcast_int_3d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_3d

  SUBROUTINE p_bcast_int_4d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_4d


  SUBROUTINE p_bcast_bool (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool

  SUBROUTINE p_bcast_bool_1d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in) :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_1d

  SUBROUTINE p_bcast_bool_2d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_2d

  SUBROUTINE p_bcast_bool_3d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_3d

  SUBROUTINE p_bcast_bool_4d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_4d

  SUBROUTINE p_bcast_char (buffer, p_source, comm)

    CHARACTER(len=*),  INTENT(inout) :: buffer
    INTEGER,           INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, LEN(buffer), p_char, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_char

  SUBROUTINE p_bcast_char_1d (buffer, p_source, comm)

    CHARACTER (*), INTENT(inout) :: buffer(:)
    INTEGER,       INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                      :: lexlength, flength

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       lexlength=LEN(buffer(1))
       flength=SIZE(buffer)
       lexlength=lexlength*flength
       CALL MPI_BCAST (buffer, lexlength, p_char, p_source, p_comm, p_error)
#ifdef DEBUG
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
               ' failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
          CALL p_abort
       END IF
#endif
    ENDIF
#endif

  END SUBROUTINE p_bcast_char_1d

  ! probe implementation

  SUBROUTINE p_probe (p_tagcount, p_tagtable, p_source, &
       p_tag, p_count, comm)

    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i = 1, p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
               flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_real_dp, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO

#endif

  END SUBROUTINE p_probe

  SUBROUTINE p_wait
#ifndef NOMPI
    INTEGER :: p_status_wait(MPI_STATUS_SIZE,p_irequest)

    CALL MPI_WAITALL(p_irequest, p_request, p_status_wait, p_error)
    p_irequest = 0
#endif
  END SUBROUTINE p_wait

  SUBROUTINE p_wait_any(return_pe)

    INTEGER, INTENT(out) :: return_pe
#ifndef NOMPI
    INTEGER :: i

    CALL MPI_WAITANY(p_irequest, p_request, i, p_status, p_error)
    IF (i == MPI_UNDEFINED) THEN
      p_irequest = 0
      return_pe = -1
    ELSE
      return_pe = p_status(MPI_SOURCE)
    ENDIF
#else
    return_pe = 0
#endif

  END SUBROUTINE p_wait_any

  SUBROUTINE p_barrier (comm)
  INTEGER ,INTENT(IN) ,OPTIONAL :: comm
#ifndef NOMPI
    INTEGER :: com
    com = p_all_comm; IF(PRESENT(comm)) com = comm
    CALL MPI_BARRIER (com, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BARRIER on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
!#endif
#endif

  END SUBROUTINE p_barrier

  FUNCTION p_sum_0d (zfield, comm) RESULT (p_sum)

    REAL(dp)                      :: p_sum
    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_0d

  FUNCTION p_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_1d

  FUNCTION p_global_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum
    REAL(dp)                      :: pe_sums(SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_REDUCE (zfield, pe_sums, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_io, p_comm, p_error)
       p_sum = SUM(pe_sums)
    ELSE
       p_sum = SUM(zfield)
    END IF
#else
    p_sum = SUM(zfield)
#endif

  END FUNCTION p_global_sum_1d

  FUNCTION p_field_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_REDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_io, p_comm, p_error)
       IF (.NOT. p_parallel_io) p_sum = 0.0_dp
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_field_sum_1d

  FUNCTION p_max_0d (zfield, comm) RESULT (p_max)

    REAL(dp)                      :: p_max
    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_0d

  FUNCTION p_max_1d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_1d

  FUNCTION p_max_2d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_2d

  FUNCTION p_max_3d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_3d

  FUNCTION p_min_0d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_0d

  FUNCTION p_min_1d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_1d

  FUNCTION p_min_2d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_2d

  FUNCTION p_min_3d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_3d

  SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:,LBOUND(recvbuf,2)) = sendbuf(:)
#endif
   END SUBROUTINE p_gather_real_1d2d

   SUBROUTINE p_gather_real_5d6d (sendbuf, recvbuf, p_dest, comm)

     REAL(dp)                      :: sendbuf(:,:,:,:,:), recvbuf(:,:,:,:,:,:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = p_all_comm
     ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' p_gather_real_5d6d failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       STOP
     END IF

#else
     recvbuf(:,:,:,:,:,LBOUND(recvbuf,6)) = sendbuf(:,:,:,:,:)
#endif
   END SUBROUTINE p_gather_real_5d6d

END MODULE mo_mpi
