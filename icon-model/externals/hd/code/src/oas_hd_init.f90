! oas_hd_init.f90 - Namcouple and routines for coupling using OASIS3-MCT
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Author: Ha Ho-Hagemann
! Contact: <ha.hagemann@hereon.de>
!_________________________________________

! Ha Ho-Hagemann {
#ifdef COUP_OAS
MODULE oas_hd_ini

!     AUTHOR: Ha Ho-Hagemann, HZG-IfK
!     -------------------------------
!     created : 2016, original source
!     updated : 04.04.2018
!

USE mod_oasis_namcouple          ! OASIS3MCT namcouple variables: e.g. coupling time step
!#ifdef OASIS3MCT
  USE mpi
  USE mod_oasis                    ! OASIS3-MCT v3 module
!#else
!  USE mod_prism_proto              ! OASIS3 prism module
!  USE mod_prism_def_partition_proto! OASIS3 prism module for partitioning
!  USE mod_prism_grids_writing      ! OASIS3 prism module for writing grid files
!  USE mod_prism_put_proto          ! OASIS3 prism module for snding
!  USE mod_prism_get_proto          ! OASIS3 prism module for receiving
!#endif

use mo_kind,        ONLY: i4, wp   ! wp = dp
   
IMPLICIT NONE

SAVE

! Debug level of OASIS
!     0 : Minimum debugging
!     1 : Debugging
!     2 : Perfs measurement
!     3 : OASIS restart production

INTEGER(kind=i4) :: ncomp_id  ! id returned by oasis_init_comp
INTEGER(kind=i4) :: kl_comm   ! Local communicator
INTEGER :: kl_mype
INTEGER :: kl_npes
INTEGER(kind=i4) :: nerror    ! return error code
CHARACTER(LEN=6)       :: modname   ! Name of the model

INTEGER(KIND=i4) :: &
  OASIS_Rcv  = 1,         & ! return code if received field
  OASIS_idle = 0,         & ! return code if nothing was done by OASIS
  OASIS_Success = 0         ! return code if no error in OASIS

REAL(KIND=wp), ALLOCATABLE :: &
     exfld (:,:),             & ! Temporary buffer for receiving
     exfld_nemo (:,:),        & ! Temporary buffer for receiving
     frcv  (:,:,:)              ! all fields received from coupled model

CONTAINS

!================================================================================
SUBROUTINE oas_hd_init 
!
IMPLICIT NONE

INTEGER :: p_error

modname = 'hdmd.x'  ! Name of the model

!------------------------------------------------------------------
! 1) Initialize the OASIS system for the component
!------------------------------------------------------------------

 CALL oasis_init_comp( ncomp_id, modname, nerror )
 IF( nerror .ne. OASIS_Success) THEN
  CALL oasis_abort( ncomp_id, 'oas_hd_init', 'Failure in oasis_init_comp' )
 ENDIF

!------------------------------------------------------------------
! 2) Get an MPI communicator for HD local communication
!------------------------------------------------------------------

 CALL oasis_get_localcomm( kl_comm, nerror )
 IF( nerror .ne. OASIS_Success ) THEN
  CALL oasis_abort( ncomp_id, 'oas_hd_init', 'Failure in oasis_get_localcomm' )
 ENDIF

  ! get local PE identification
  CALL MPI_COMM_RANK (kl_comm, kl_mype, p_error)
  IF (p_error /= MPI_SUCCESS) THEN
    CALL oasis_abort( ncomp_id, 'oas_hd_init', 'Failure in MPI_COMM_RANK' )
  END IF

  ! get number of available PEs
  CALL MPI_COMM_SIZE (kl_comm, kl_npes, p_error)
  IF (p_error /= MPI_SUCCESS) THEN
    CALL oasis_abort( ncomp_id, 'oas_hd_init', 'Failure in MPI_COMM_SIZE' )
  END IF

END SUBROUTINE oas_hd_init

!================================================================================
SUBROUTINE oas_hd_finalize

!!---------------------------------------------------------------------
!!              ***  ROUTINE oas_hd_finalize  ***
!!
!! ** Purpose : - Finalizes the coupling. If MPI_init has not been
!!      called explicitly before oas_hd_init it will also close
!!      MPI communication.
!----------------------------------------------------------------------

IMPLICIT NONE

write(0,*) 'Call oasis_terminate in HD'
IF (ALLOCATED (exfld))  DEALLOCATE(exfld)
IF (ALLOCATED (exfld_nemo))  DEALLOCATE(exfld_nemo)
IF (ALLOCATED (frcv))  DEALLOCATE(frcv)
CALL oasis_terminate ( nerror )

END SUBROUTINE oas_hd_finalize

END MODULE oas_hd_ini
#endif
! Ha Ho-Hagemann }
