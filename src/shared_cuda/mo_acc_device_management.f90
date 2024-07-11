! Module acc_device_management
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

MODULE mo_acc_device_management

!------------------------------------------------------------------------------
!
! Description : 
!  This module provides routines to set up the GPU device as well wrapper
!  for cuda function calls. The module is empty if preprocessor variable 
!  _OPENACC is not set.
!  Note: routines using STOP in case of error instead of model_abort
!  as they may be called before mpi init
!  
! Routines contained:
!  - cudaGetErrorString     : get cuda error
!  - cudaGetDeviceCount     : get number of devices
!  - cudaSetDevice          : set device
!  - cudaDeviceReset        : reset device
!  - cudaMemGetInfo         : get memory usage
!  - cudaDeviceSynchronize  : halts CPU/Host thread until GPU has processed all 
!                             tasks
!  - my_cudaErrorCheck      : check and print cuda error if any
!  - initAccDevice          : init device
!  - checkAccDevice         : check that device is correctly initialized
!  - setAccDevice           : set device to current thread
!  - runSmallAccKernel      : run small test OpenACC kernel
!  - finalizeAccDevice      : finalize device
!  - printGPUMem            : print current GPU usage
!
!------------------------------------------------------------------------------

#ifdef _OPENACC

! Modules used:
!------------------------------------------------------------------------------

USE, INTRINSIC :: iso_c_binding
USE mpi
USE mo_kind,      ONLY: wp
USE mo_mpi,       ONLY: p_comm_work, p_comm_rank, num_work_procs, p_max
USE mo_exception, ONLY: message, message_text

#ifdef ICON_USE_CUDA_GRAPH
#if (__NVCOMPILER_MAJOR__ == 23 && __NVCOMPILER_MINOR__ <= 9) || (__NVCOMPILER_MAJOR__ == 99 && __NVCOMPILER_PATCHLEVEL__ <= 239500)
#define USE_OLD_GRAPH_INTERFACE
USE openacc, ONLY: accGraph => accgraph,    accx_begin_capture,       accx_end_capture,       accx_graph_exec,   accx_graph_delete
#else
USE openacc, ONLY: accGraph => acc_graph_t, accx_async_begin_capture, accx_async_end_capture, accx_graph_launch, accx_graph_delete
#endif
#endif

!------------------------------------------------------------------------------

IMPLICIT NONE


!------------------------------------------------------------------------------


PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:

PUBLIC :: initAccDevice, checkAccDevice, finalizeAccDevice, printGPUMem
#ifdef ICON_USE_CUDA_GRAPH
PUBLIC :: accGraph, accBeginCapture, accEndCapture, accGraphLaunch, accGraphDelete
#endif


!==============================================================================
! Parameters

LOGICAL, PARAMETER :: debug = .FALSE.  ! set to True, to get extensive debug output


!==============================================================================
! Module variables

INTEGER :: numid=0             ! MPI number of ranks
INTEGER :: myid = -2           ! MPI rank
INTEGER :: mylocid = -2           ! LOCAL node rank
INTEGER :: numnodev=0          ! number of MPI ranks without device

INTEGER :: numdev=0            ! number of CUDA devices on node
INTEGER :: mydev=-2            ! CUDA device number assigned to this rank

INTEGER :: numnode=0           ! number of nodes

LOGICAL :: init_done = .FALSE. ! has the device been initialized

ENUM, BIND(C) !:: cudaError
  ENUMERATOR :: cudaSuccess=0
  ENUMERATOR :: cudaErrorMissingConfiguration=1
  ENUMERATOR :: cudaErrorMemoryAllocation=2
  ENUMERATOR :: cudaErrorInitializationError=3
  ENUMERATOR :: cudaErrorLaunchFailure=4
  ENUMERATOR :: cudaErrorPriorLaunchFailure=5
  ENUMERATOR :: cudaErrorLaunchTimeout=6
  ENUMERATOR :: cudaErrorLaunchOutOfResources=7
  ENUMERATOR :: cudaErrorInvalidDeviceFunction=8
  ENUMERATOR :: cudaErrorInvalidConfiguration=9
  ENUMERATOR :: cudaErrorInvalidDevice=10
  ENUMERATOR :: cudaErrorInvalidValue=11
  ENUMERATOR :: cudaErrorInvalidPitchValue=12
  ENUMERATOR :: cudaErrorInvalidSymbol=13
  ENUMERATOR :: cudaErrorMapBufferObjectFailed=14
  ENUMERATOR :: cudaErrorUnmapBufferObjectFailed=15
  ENUMERATOR :: cudaErrorInvalidHostPointer=16
  ENUMERATOR :: cudaErrorInvalidDevicePointer=17
  ENUMERATOR :: cudaErrorInvalidTexture=18
  ENUMERATOR :: cudaErrorInvalidTextureBinding=19
  ENUMERATOR :: cudaErrorInvalidChannelDescriptor=20
  ENUMERATOR :: cudaErrorInvalidMemcpyDirection=21
  ENUMERATOR :: cudaErrorAddressOfConstant=22
  ENUMERATOR :: cudaErrorTextureFetchFailed=23
  ENUMERATOR :: cudaErrorTextureNotBound=24
  ENUMERATOR :: cudaErrorSynchronizationError=25
  ENUMERATOR :: cudaErrorInvalidFilterSetting=26
  ENUMERATOR :: cudaErrorInvalidNormSetting=27
  ENUMERATOR :: cudaErrorMixedDeviceExecution=28
  ENUMERATOR :: cudaErrorCudartUnloading=29
  ENUMERATOR :: cudaErrorUnknown=30
  ENUMERATOR :: cudaErrorNotYetImplemented=31
  ENUMERATOR :: cudaErrorMemoryValueTooLarge=32
  ENUMERATOR :: cudaErrorInvalidResourceHandle=33
  ENUMERATOR :: cudaErrorNotReady=34
  ENUMERATOR :: cudaErrorInsufficientDriver=35
  ENUMERATOR :: cudaErrorSetOnActiveProcess=36
  ENUMERATOR :: cudaErrorInvalidSurface=37
  ENUMERATOR :: cudaErrorNoDevice=38
  ENUMERATOR :: cudaErrorECCUncorrectable=39
  ENUMERATOR :: cudaErrorSharedObjectSymbolNotFound=40
  ENUMERATOR :: cudaErrorSharedObjectInitFailed=41
  ENUMERATOR :: cudaErrorUnsupportedLimit=42
  ENUMERATOR :: cudaErrorDuplicateVariableName=43
  ENUMERATOR :: cudaErrorDuplicateTextureName=44
  ENUMERATOR :: cudaErrorDuplicateSurfaceName=45
  ENUMERATOR :: cudaErrorDevicesUnavailable=46
  ENUMERATOR :: cudaErrorStartupFailure=127
  ENUMERATOR :: cudaErrorApiFailureBase=10000
END ENUM ! cudaError

!==============================================================================
! Interface definition

#ifndef __HIP__

INTERFACE ! [['const', 'char', '*'], 'cudaGetErrorString', [['cudaError_t', None, 'error']]]
  FUNCTION cudaGetErrorString(error) RESULT( res ) BIND(C, name="cudaGetErrorString")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (KIND(cudaSuccess)), VALUE :: error
    TYPE(c_ptr) :: res
  END FUNCTION cudaGetErrorString
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaGetDeviceCount', [['int', '*', 'count']]]
  FUNCTION cudaGetDeviceCount(count) RESULT( res ) BIND(C, name="cudaGetDeviceCount")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER(c_int) :: count
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaGetDeviceCount
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaSetDevice', [['int', None, 'device']]]
  FUNCTION cudaSetDevice(device) RESULT( res ) BIND(C, name="cudaSetDevice")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER(c_int), VALUE :: device
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaSetDevice
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaDeviceReset', ['void']]
  FUNCTION cudaDeviceReset() RESULT( res ) BIND(C, name="cudaDeviceReset")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaDeviceReset
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaMemGetInfo', [['int', '*', 'free'],['int', '*', 'total']]]
  FUNCTION cudaMemGetInfo(free, total) RESULT( res ) BIND(C, name="cudaMemGetInfo")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (c_size_t) :: free, total
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaMemGetInfo
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaDeviceSynchronize', ['void']]
  FUNCTION cudaDeviceSynchronize() RESULT( res ) BIND(C,name="cudaDeviceSynchronize")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaDeviceSynchronize
END INTERFACE

#else ! HIP versions

INTERFACE ! [['const', 'char', '*'], 'cudaGetErrorString', [['cudaError_t', None, 'error']]]
  FUNCTION cudaGetErrorString(error) RESULT( res ) BIND(C, name="hipGetErrorString")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (KIND(cudaSuccess)), VALUE :: error
    TYPE(c_ptr) :: res
  END FUNCTION cudaGetErrorString
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaGetDeviceCount', [['int', '*', 'count']]]
  FUNCTION cudaGetDeviceCount(count) RESULT( res ) BIND(C, name="hipGetDeviceCount")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER(c_int) :: count
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaGetDeviceCount
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaSetDevice', [['int', None, 'device']]]
  FUNCTION cudaSetDevice(device) RESULT( res ) BIND(C, name="hipSetDevice")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER(c_int), VALUE :: device
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaSetDevice
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaDeviceReset', ['void']]
  FUNCTION cudaDeviceReset() RESULT( res ) BIND(C, name="hipDeviceReset")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaDeviceReset
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaMemGetInfo', [['int', '*', 'free'],['int', '*', 'total']]]
  FUNCTION cudaMemGetInfo(free, total) RESULT( res ) BIND(C, name="hipMemGetInfo")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (c_size_t) :: free, total
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaMemGetInfo
END INTERFACE

INTERFACE ! [['cudaError_t', None], 'cudaDeviceSynchronize', ['void']]
  FUNCTION cudaDeviceSynchronize() RESULT( res ) BIND(C,name="hipDeviceSynchronize")
    USE, INTRINSIC :: ISO_C_BINDING
    IMPORT cudaSuccess
    IMPLICIT NONE
    INTEGER (KIND(cudaSuccess)) :: res
  END FUNCTION cudaDeviceSynchronize
END INTERFACE

#endif

!==============================================================================
! Module procedures in acc_device_management
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure "my_cudaErrorCheck" in "acc_device_management"
!------------------------------------------------------------------------------

SUBROUTINE my_cudaErrorCheck( err )

!------------------------------------------------------------------------------
!
! Description:
!   Check cuda error code
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER (KIND(cudaSuccess)), INTENT(in) :: err

! Locals:
! --------------------
  INTEGER :: i, j, ip, ipen
  CHARACTER, POINTER :: string_err(:)

!------------------------------------------------------------------------------
! Begin Subroutine my_cudaErrorCheck
!------------------------------------------------------------------------------


  IF ( err /= cudaSuccess ) THEN
     CALL C_F_POINTER( cudaGetErrorString( err ), string_err, [1] )
     i = 1
     DO WHILE ( string_err( i ) /= C_NULL_CHAR )
        i = i + 1
     ENDDO
     WRITE(0,'(a)') 'ERROR: CUDA Error Detected'
     WRITE(0,*) string_err( 1:i )
     STOP
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure my_cudaErrorCheck
!------------------------------------------------------------------------------

END SUBROUTINE my_cudaErrorCheck

!==============================================================================

!==============================================================================
!+ Module procedure "initAccDevice" in "acc_device_management"
!------------------------------------------------------------------------------

SUBROUTINE initAccDevice()

!------------------------------------------------------------------------------
!
! Description:
!
! Initialize GPUs. In case of multi-thread run associates one GPU to each 
! thread. The last COSMO_NPROC_NODEVICE ranks are not associated to any GPU.
! COSMO_NPROC_NODEVICE is an environement variable which should be set before
! executing the model
!
!------------------------------------------------------------------------------

! Local parameters:
! ----------------

  IMPLICIT NONE


  CHARACTER(len=256) :: snumid
  CHARACTER(len=256) :: smyid
  CHARACTER(len=256) :: smylocid

!------------------------------------------------------------------------------
! Begin Subroutine initAccDevice
!------------------------------------------------------------------------------


! get rank of this MPI task from environment since
! CUDA setup needs to be done _before_ MPI init
  CALL getenv("SLURM_NPROCS", snumid)
  CALL getenv("SLURM_PROCID", smyid)
  CALL getenv("SLURM_LOCALID", smylocid)

  IF (LEN_TRIM(snumid) == 0 .OR. LEN_TRIM(smyid) == 0 ) THEN
    ! Fall back to mvapich env for backward compatibility
    CALL getenv("MV2_COMM_WORLD_SIZE", snumid)
    CALL getenv("MV2_COMM_WORLD_RANK", smyid)
    CALL getenv("MV2_COMM_WORLD_LOCAL_RANK", smylocid)
    IF (LEN_TRIM(snumid) == 0 .OR. LEN_TRIM(smylocid) == 0) THEN
      WRITE(0,'(a)') 'INFO:    Env variables SLURM_NPROCS and SLURM_LOCALID'
      WRITE(0,'(a)') '         are not set'
      WRITE(0,'(a)') '   -->   no cudaSetDevice() called before MPI_Init()'
      WRITE(0,'(a)') '   This could be and issue on system with multi-GPUs per nodes'
      RETURN
    ENDIF
  ENDIF
  READ(snumid,'(i8)') numid
  READ(smyid,'(i8)') myid
  READ(smylocid,'(i8)') mylocid

  IF (debug) WRITE(*,*) 'DBG: ', myid, ' NPROCS = ',numid
  IF (debug) WRITE(*,*) 'DBG: ', myid, ' MYID = ',myid
  IF (debug) WRITE(*,*) 'DBG: ', myid, ' MYLOCID = ',mylocid

  
  ! get number of procesors not to connect to a device
  CALL getenv("COSMO_NPROC_NODEVICE", snumid)
  IF (LEN_TRIM(snumid) == 0) THEN
    numnodev = 0
  ELSE
    READ(snumid,'(i8)') numnodev
  ENDIF
  IF (debug) WRITE(*,*) 'DBG: ', myid, ' COSMO_NPROC_NODEVICE = ',numnodev

  ! check if rank should init device
  IF (myid < numid - numnodev) THEN
    CALL setAccDevice()
    IF (debug) WRITE(*,*) 'DBG: ', myid, ' setAccDevice ok on ',  myid
  ELSE
    IF (debug) WRITE(*,*) 'DBG: ', myid, ' excluding rank ', myid
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure initAccDevice
!------------------------------------------------------------------------------

END SUBROUTINE initAccDevice

!==============================================================================
  
!==============================================================================
!+ Module procedure "checkAccDevice" in "acc_device_management"
!------------------------------------------------------------------------------

SUBROUTINE checkAccDevice(mycomm, mynumranks, myrank)

!------------------------------------------------------------------------------
!
! Description:
!
! Check that the GPU device is correctly set - init the device otherwise
! NOTE: this subroutine should only be called from MPI ranks which actually
! should have a accelerator device associeted (thus excluding I/O PEs)
!
!------------------------------------------------------------------------------

! Modules:
! --------
  USE iso_c_binding

  IMPLICIT NONE

! externals
! --------
  INTERFACE
    FUNCTION gethostid() BIND(C)
      USE iso_c_binding
      INTEGER (C_INT) :: gethostid
    END FUNCTION gethostid
  END INTERFACE

! Subroutine arguments:
! --------------------

  INTEGER, INTENT(in) :: &
       mycomm, &           ! mpi communicator
       mynumranks, &       ! number of mpi ranks 
       myrank              ! mpi ranks

! locals:
!--------
  INTEGER :: hostid, ierr, i, j, numlocal
  INTEGER :: hostids(mynumranks), localprocs(mynumranks)
  INTEGER :: uniquehostid(mynumranks)

!------------------------------------------------------------------------------
! Begin Subroutine checkAccDevice
!------------------------------------------------------------------------------

  ! if init has already been done, check consistency, otherwise init
  IF (init_done) THEN

    IF (numdev <= 0) THEN
      WRITE(0,'(a,i)') 'ERROR: inconsistent number of devices encountered', numdev
      STOP
    ENDIF
    IF (mydev < 0) THEN
      WRITE(0,'(a,i)') 'ERROR: inconsistent device number encountered', mydev
      STOP
    ENDIF
    IF (myrank /= myid) THEN
      WRITE(0,'(a,i,i)') 'ERROR: inconsistent MPI ranks encountered', myrank, myid
      STOP
    ENDIF

!!$  if (numid - numnodev /= mynumranks) then
!!$    write(0,'(a,i,i,i)') 'ERROR: inconsistent number of ranks on communicator and passed as argument', numid, numnodev, mynumranks
!!$    stop
!!$  endif
    CALL mpi_comm_size(mycomm, numid, ierr)
    IF (ierr /= 0) THEN
      WRITE(0,'(a,i)') 'ERROR: problem with mpi_comm_size', ierr
      STOP
    ENDIF
    IF (numid /= mynumranks) THEN
      WRITE(0,'(a,i,i)') 'ERROR: inconsistent number of ranks on communicator and passed as argument', numid, mynumranks
      STOP
    ENDIF

  ELSE

    ! setup defaults
    myid = myrank
    numnodev = 0

    ! get number of ranks
    CALL mpi_comm_size(mycomm, numid, ierr)
    IF (ierr /= 0) THEN
      WRITE(0,'(a,i)') 'ERROR: problem with mpi_comm_size', ierr
      STOP
    ENDIF
    IF (numid /= mynumranks) THEN
      WRITE(0,'(a,i,i)') 'ERROR: inconsistent number of ranks on communicator and passed as argument', numid, mynumranks
      STOP
    ENDIF

    ! and now set device
    CALL setAccDevice()

  ENDIF
  
  ! get the host ID for this rank
  hostid = gethostid()

  ! special treatment for just 1 rank (since MPI might not be in use)
  IF (numid == 1) THEN

    numlocal = 1
    localprocs(1) = 0
    numnode = 1
    uniquehostid(1) = hostid

  ELSE

    ! get the hostids so we can determine what other processes are on this node
    CALL mpi_allgather( hostid, 1, mpi_integer, hostids, 1, mpi_integer, mycomm, ierr )
    IF (ierr /= 0) THEN
      WRITE(0,'(a,i)') 'ERROR: problem with mpi_allgather', ierr
      STOP
    ENDIF
    IF (debug) WRITE(*,*) 'DBG: ', myid, ' gathered hostids'

    ! determine who's on this node
    numlocal = 0
    localprocs = 0
    DO i = 1, numid
      IF (hostid .EQ. hostids(i)) THEN
        localprocs(i) = numlocal
        numlocal = numlocal+1
      ENDIF
    ENDDO

    ! determine number of nodes
    numnode = 1
    uniquehostid(1) = hostids(1)
    DO i = 2, numid
      DO j = 1, numnode
        IF (hostids(i) == uniquehostid(j)) EXIT
      ENDDO
      IF (j > numnode) THEN
        numnode = numnode + 1
        uniquehostid(numnode) = hostids(i)
      ENDIF
    ENDDO
    IF (debug) WRITE(*,*) 'DBG: ', myid, ' number of nodes = ',numnode

  ENDIF

  ! print a warning if the number of devices is less then the number
  ! of processes on this node. Some NVIDIA devices do not support multiple context on a
  ! single device. Hence, having multiple processes share devices is not recommended.
  ! It may work, but can also fail.
  !if (numdev < numlocal) then
  !  write(0,'(a,i,i)') 'ERROR: number of ranks per node is larger than the number of devices', numdev, numlocal
  !  stop
  !endif
  !if (mydev /= localprocs(myid+1)) then
  !  write(0,'(a,i,i)') 'ERROR: mapping between MPI rank and accelerator device is not linear', mydev, localprocs(myid+1)
  !  stop
  !endif
  
  ! inform user
  IF (myid == 0) THEN
    WRITE(*,'(a,1x,i3,1x,a,1x,i3,1x,a,1x,i3,1x,a)') &
         'INFO: Running',numid,'MPI tasks on',numnode,'nodes with',numdev,'devices'
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure checkAccDevice
!------------------------------------------------------------------------------

END SUBROUTINE checkAccDevice

!==============================================================================
  
!==============================================================================
!+ Module procedure "setAccDevice" in "acc_device_manangement"
!------------------------------------------------------------------------------

SUBROUTINE setAccDevice()
!------------------------------------------------------------------------------
!
! Description:
!
! Set device for this MPI rank both for cuda driver and OpenACC runtime
!
!------------------------------------------------------------------------------

! Modules:
! --------
  USE openacc

  IMPLICIT NONE

! Locals:
! -------
  INTEGER :: ierr

!------------------------------------------------------------------------------
! Begin Subroutine setAccDevice
!------------------------------------------------------------------------------
  
  
  ! setup reasonable defaults
  numdev = 1
  mydev = 0

  ! get number of available devices on this node
  ierr = cudaGetDeviceCount(numdev)
  CALL my_cudaErrorCheck(ierr)
  IF (numdev < 1) THEN
    WRITE(0,'(a,i)') 'ERROR: No CUDA devices found on this node', numdev
    STOP
  ENDIF
  IF (debug) WRITE(*,*) 'DBG: ', myid, ' rank/#dev = ',myid,numdev

  ! set device for this MPI rank
  ! NOTE: if a node is oversubscribed, two MPI tasks will get the same
  !       CUDA device. This will be checked later after MPI has been
  !       initialized
  ! NOTE: we need to set the device both for the cuda driver as well
  !       as for the OpenACC runtime
  mydev = MOD(mylocid, numdev)
  ierr = cudaSetDevice(mydev)
  CALL my_cudaErrorCheck(ierr)
  ierr = cudaDeviceSynchronize()
  CALL my_cudaErrorCheck(ierr)
  CALL acc_set_device_num(mydev, acc_device_nvidia)
  CALL acc_init(acc_device_nvidia)
  IF (debug) WRITE(*,*) 'DBG: ', myid, ' rank/dev = ',myid,mydev

  ! run small example kernel
  CALL runSmallAccKernel()
  IF (debug) WRITE(*,*) 'DBG: ', myid, ' ran small init kernel'

  ! finally, this is done
  init_done = .TRUE.

  IF (myid == 0) THEN
    WRITE(*,'(a)') 'INFO: Accelerator devices initialized and set'
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure setAccDevice
!------------------------------------------------------------------------------

END SUBROUTINE setAccDevice

!==============================================================================
  
!==============================================================================
!+ Module procedure "runSmallAccKernel" in "acc_device_manangement"
!------------------------------------------------------------------------------

SUBROUTINE runSmallAccKernel()

!------------------------------------------------------------------------------
!
! Description:
!
! Run small test OpenACC kernel. This will initialize the device.
! This is necessary to get accurate timing
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Locals:
! -------

  INTEGER :: dummy(16), idummy

!------------------------------------------------------------------------------
! Begin Subroutine runSmallAccKernel
!------------------------------------------------------------------------------
  
  !$ACC PARALLEL ASYNC(1)
  !$ACC LOOP
  DO idummy = 1,16
    dummy(idummy) = idummy
  ENDDO
  !$ACC END PARALLEL
  !$ACC WAIT(1)
  IF (myid == 0) THEN
#ifdef _OPENACC
    WRITE(*,'(a)') 'INFO: Running with OpenAcc directives'
#else
    WRITE(*,'(a)') 'INFO: Running without OpenAcc directives'
#endif
  ENDIF
  IF (myid == 12342134) THEN
    ! this is needed to avoid the compiler to optimize away the dummy
    ! array, but due to the IF above, this is unlikely to be reached.
    ! only, the compiler can not know this 
    ! (avoid dead code elimination)
    PRINT*, 'Init print: ', SUM(dummy)
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure runSmallAccKernel
!------------------------------------------------------------------------------

END SUBROUTINE runSmallAccKernel 

!==============================================================================
  
!==============================================================================
!+ Module procedure "finalizeAccDevice" in "acc_device_manangement"
!------------------------------------------------------------------------------

SUBROUTINE finalizeAccDevice()
IMPLICIT NONE

INTEGER :: ierr

IF (mydev > 0) THEN
  ierr = cudaDeviceReset()
  CALL my_cudaErrorCheck(ierr)
ENDIF

END SUBROUTINE finalizeAccDevice

!==============================================================================
  
!==============================================================================
!+ Module procedure "printGPUMem" in "acc_device_manangement"
!------------------------------------------------------------------------------

SUBROUTINE printGPUMem(mtag)
!------------------------------------------------------------------------------
!
! Description:
!
! Print current GPU memory usage
!
!------------------------------------------------------------------------------

  USE, INTRINSIC :: ISO_C_BINDING

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN)    ::   mtag      ! call-site information

! Locals:
! -------
  INTEGER(c_size_t) :: ifree, itotal
  INTEGER :: izerror, icudaerr
  REAL(wp) :: input(2), output(2)
  CHARACTER (LEN=80)  :: yzerrmsg

!------------------------------------------------------------------------------
! Begin Subroutine printGPUMem
!------------------------------------------------------------------------------


  ifree  = 0
  itotal = 0

  icudaerr = cudaMemGetInfo(ifree,itotal)
  input(1) = REAL(ifree,wp)/1073741824._wp
  input(2) = REAL(itotal,wp)/1073741824._wp

  IF (num_work_procs > 1) THEN
    output = p_max( input, comm=p_comm_work )
  ELSE
    output = input
  END IF
  IF (p_comm_rank(p_comm_work) == 0) THEN
    WRITE(message_text,'(a,E12.6,a,E12.6)') 'used GB ', output(2)-output(1), ' free GB ', output(1)
    CALL message(TRIM(mtag), message_text)
  ENDIF
!------------------------------------------------------------------------------
! End of module procedure printGPUMem
!------------------------------------------------------------------------------

END SUBROUTINE printGPUMem

!==============================================================================
  
!==============================================================================
!+ Module procedures for the graphs. Needed due to interface change in
!+ NVHPC 23.11
!------------------------------------------------------------------------------
#ifdef ICON_USE_CUDA_GRAPH

SUBROUTINE accBeginCapture(queue)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: queue

#ifndef USE_OLD_GRAPH_INTERFACE
  CALL accx_async_begin_capture(queue)
#else
  CALL accx_begin_capture(queue)
#endif

END SUBROUTINE accBeginCapture


SUBROUTINE accEndCapture(queue, graph)
  IMPLICIT NONE

  INTEGER,        INTENT(IN)  :: queue
  TYPE(accGraph), INTENT(OUT) :: graph

#ifndef USE_OLD_GRAPH_INTERFACE
  CALL accx_async_end_capture(queue, graph)
#else
  CALL accx_end_capture(graph, queue)
#endif

END SUBROUTINE accEndCapture


SUBROUTINE accGraphLaunch(graph, queue)
  IMPLICIT NONE

  TYPE(accGraph), INTENT(IN) :: graph
  INTEGER,        INTENT(IN) :: queue

#ifndef USE_OLD_GRAPH_INTERFACE
  CALL accx_graph_launch(graph, queue)
#else
  CALL accx_graph_exec(graph, queue)
#endif

END SUBROUTINE accGraphLaunch


SUBROUTINE accGraphDelete(graph)
  IMPLICIT NONE

  TYPE(accGraph), INTENT(IN) :: graph

  CALL accx_graph_delete(graph)

END SUBROUTINE accGraphDelete

! ICON_USE_CUDA_GRAPH
#endif

! _OPENACC
#endif

END MODULE mo_acc_device_management
