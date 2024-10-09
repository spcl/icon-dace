! oas_hd.f90 - Controls, definitions and variables for communication of HD with OASIS
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
MODULE oas_hd

!**** *oas_hd* - Controls, definitions and variables
!                   for OASIS communications
!
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
   
USE mo_kind,        ONLY: wp, i4    ! wp = dp
USE mo_mpi,         ONLY: p_nprocs
USE mo_control,     ONLY: nproca, nprocb
USE mo_time_control,ONLY: delta_time, get_time_step
USE mo_couple_to_ocean,  ONLY: nxocean, nyocean, discharge_on_ocean
USE oas_hd_ini,          ONLY: kl_mype, kl_npes
USE mo_coupling, ONLY: get_grid_size, get_local_partition, &
                       fdir_hd, runoff_s, runoff_dr, &
                       get_coupling_type

IMPLICIT NONE

SAVE

! Debug level of OASIS
!     0 : Minimum debugging
!     1 : Debugging
!     2 : Perfs measurement
!     3 : OASIS restart production

  
INTEGER(KIND=i4) :: &
  debug_oasis,             & ! Debug level of OASIS
  ncomp_id,                & ! id returned by oasis_init_comp   
  kl_comm,                 & ! Local communicator
  nerror                     ! return error code

INTEGER(KIND=i4) :: &
  hdshape(4),hdshape_nemo(4) ! shape of arrays passed to OASIS3

INTEGER(KIND=i4) :: &
  OASIS_Rcv  = 1,          & ! return code if received field
  OASIS_idle = 0,          & ! return code if nothing was done by OASIS
  OASIS_Success = 0          ! return code if no error in OASIS

INTEGER(KIND=i4), PARAMETER :: OASIS_Box_Params = 5

INTEGER(KIND=i4) :: &
  nulout=6                   ! unit number for joblog file
    
TYPE :: FLD_CPL                      ! Type for coupling field information
  LOGICAL               :: laction   ! To be coupled or not
  CHARACTER(LEN = 16)   :: clname    ! Name of the coupling field
  CHARACTER(LEN = 1)    :: clgrid    ! Grid type
  INTEGER(KIND=i4):: nid      ! Id of the field
END TYPE FLD_CPL

TYPE(FLD_CPL), ALLOCATABLE :: &
     srcv(:),                 & ! All fields to be received
     ssnd(:)                    ! All fields to be sent

REAL(KIND=wp), ALLOCATABLE :: &
     exfld (:,:),             & ! Temporary buffer for receiving
     exfld_nemo (:,:),        & ! Temporary buffer for receiving
     frcv  (:,:,:)              ! all fields received from coupled model

LOGICAL                 :: lpe_cpl = .TRUE.
CHARACTER(LEN=6)        :: modname   ! Name of the model

INTEGER(KIND=i4) :: &
  nfld,                    & !
  nfld_snd_tot,            & ! total number of fields to be sent
  nfld_rcv_tot,            & ! total number of fields to be received
  nfld_snd_oce,            & !
  nfld_rcv_oce

INTEGER(KIND=i4), ALLOCATABLE :: &
  nlev_snd_oce(:),         & ! number of levels per field to be received for coupling
  nlev_rcv_oce(:)            !

CHARACTER(LEN=16), ALLOCATABLE :: &
  nam_snd_oce(:),          & ! names of fields to be sent for coupling
  nam_rcv_oce(:)             !
  
INTEGER(KIND=i4)   :: istep

CONTAINS

!================================================================================
SUBROUTINE oas_hd_define
!
IMPLICIT NONE

INTEGER(KIND=i4) :: hd_part, hd_part_nemo
INTEGER(KIND=i4) :: il_paral(OASIS_Box_Params)
INTEGER(KIND=i4) :: ji
INTEGER          :: global_extent_lon, global_extent_lat
INTEGER          :: local_extent_lon, local_extent_lat
INTEGER          :: local_start_lon, local_start_lat
!!
!!--------------------------------------------------------------------
if (kl_mype == 0) then
 WRITE(0,'(A)')'HD: Setting up OASIS3-MCT interface'
endif

! -----------------------------------------------------------------
! ... Define the partition
! -----------------------------------------------------------------

!---- HD grid 'box' decomp
call get_grid_size(global_extent_lon, global_extent_lat)
call get_local_partition( &
  local_extent_lon, local_extent_lat, local_start_lon, local_start_lat)
call decomp_def( &
  il_paral, global_extent_lon, global_extent_lat, &
  local_extent_lon, local_extent_lat, local_start_lon, local_start_lat)
hdshape(:) = 1
hdshape(2) = local_extent_lon
hdshape(4) = local_extent_lat
!print*, "+++ HD1: il_paral=", il_paral
call OASIS_def_partition  (hd_part, il_paral, nerror, name='parthd')
IF( nerror .ne. OASIS_Success )THEN
    CALL oasis_abort (ncomp_id, 'oas_hd_define', 'Failure in oasis_def_partition of parthd' )
ENDIF

!---- NEMO grid 'box' decomp
! MoHa: here I assume that HD runs on a single process, if this is not the case
!       anymore, the decomposition information has to be adjusted accordingly
call decomp_def(il_paral, nxocean, nyocean, nxocean, nyocean, 1, 1)
hdshape_nemo(:) = 1
hdshape_nemo(2) = nxocean
hdshape_nemo(4) = nyocean
!print*, "+++ HD2: il_paral=", il_paral
call OASIS_def_partition  (hd_part_nemo, il_paral, nerror, name='parthd_nemo')
IF( nerror .ne. OASIS_Success )THEN
    CALL oasis_abort (ncomp_id, 'oas_hd_define', 'Failure in oasis_def_partition of parthd_nemo' )
ENDIF

!if (kl_mype == 0 ) then
! write(0,*) '++++HD: end def_partition'
!endif

! -----------------------------------------------------------------
! ... Initialize some variables
! ----------------------------------------------------------------

nfld_snd_oce = 0
nfld_rcv_oce = 0

nfld_snd_tot = 0
nfld_rcv_tot = 0

! -----------------------------------------------------------------
! ... Define list of SENT variable names per coupling
! ----------------------------------------------------------------

nfld = 0

    nfld_snd_oce = 1
    ALLOCATE ( nlev_snd_oce (nfld_snd_oce), STAT=nerror )
    ALLOCATE ( nam_snd_oce  (nfld_snd_oce), STAT=nerror )
    nlev_snd_oce = (/1 /)
    nam_snd_oce = (/'RDC2NEMO'/)
    DO ji = 1, nfld_snd_oce
      nfld = nfld + nlev_snd_oce(ji)
    ENDDO

! total number of fields to be sent:
nfld_snd_tot = nfld
  
! -----------------------------------------------------------------
! ... Define list of RECEIVED variable names per coupling
! ----------------------------------------------------------------

nfld = 0
  
    nfld_rcv_oce = 2
    ALLOCATE ( nlev_rcv_oce (nfld_rcv_oce), STAT=nerror )
    ALLOCATE ( nam_rcv_oce  (nfld_rcv_oce), STAT=nerror )
    nlev_rcv_oce = (/1,1 /)
    nam_rcv_oce = (/'RUNOFF_S','RUNOFF_G' /)
    DO ji = 1, nfld_rcv_oce
      nfld = nfld + nlev_rcv_oce(ji)
    ENDDO

! total number of fields to be received:
nfld_rcv_tot = nfld

! -----------------------------------------------------------------
! ... write variable names for all active couplings to one structure
! ----------------------------------------------------------------

! Allocate memory for data exchange:
ALLOCATE( ssnd(0:nfld_snd_tot), stat = nerror )
ALLOCATE( srcv(0:nfld_rcv_tot), stat = nerror )
ALLOCATE( exfld(hdshape(2),hdshape(4)), stat = nerror )
IF ( nerror > 0 ) THEN
  CALL oasis_abort( ncomp_id, 'oas_hd_define', 'Failure in allocating exfld, ssnd or srcv' )
  RETURN
ENDIF
ALLOCATE( exfld_nemo(hdshape_nemo(2),hdshape_nemo(4)), stat = nerror )
IF ( nerror > 0 ) THEN
  CALL oasis_abort( ncomp_id, 'oas_hd_define', 'Failure in allocating exfld_nemo, ssnd or srcv' )
  RETURN
ENDIF


! fill ssnd with names of fields to be sent
nfld = 0
DO ji = 1, nfld_snd_oce
   nfld = nfld + 1
   ssnd(nfld)%clname = TRIM(nam_snd_oce(ji))
ENDDO

nfld = 0
DO ji = 1, nfld_rcv_oce
   nfld = nfld + 1
   srcv(nfld)%clname = TRIM(nam_rcv_oce(ji))
ENDDO

! ----------------------------------------------------------------------------
! ... Variable selection
! ----------------------------------------------------------------------------

! This is still preliminary: set laction for all (this way no subsets of
! variables can be selected; this has to be changed yet)
ssnd(:)%laction = .TRUE.
srcv(:)%laction = .TRUE.
      
! -----------------------------------------------------------------
! ... Announce variables to be sent and to be received
! -----------------------------------------------------------------

! Announce variables to be sent:
DO ji = 1, nfld_snd_tot
  IF (get_coupling_type() .EQ. 2) THEN
    CALL oasis_def_var( ssnd(ji)%nid, ssnd(ji)%clname, hd_part_nemo, &
      (/ 2, 0/), OASIS_Out, hdshape_nemo, OASIS_REAL, nerror )
  ELSE
    CALL oasis_def_var( ssnd(ji)%nid, ssnd(ji)%clname, hd_part, &
      (/ 2, 0/), OASIS_Out, hdshape, OASIS_REAL, nerror )
  ENDIF
    IF ( nerror /= OASIS_Success ) CALL oasis_abort( ssnd(ji)%nid, &
      'oas_hd_define', 'Failure in oasis_def_var for '//TRIM(ssnd(ji)%clname) )
ENDDO
      
! Announce variables to be received:
DO ji = 1, nfld_rcv_tot
  CALL oasis_def_var( srcv(ji)%nid, srcv(ji)%clname, hd_part, &
    (/ 2, 0/), OASIS_In, hdshape, OASIS_REAL, nerror )
  IF ( nerror /= OASIS_Success ) CALL oasis_abort( srcv(ji)%nid, &
    'oas_hd_define', 'Failure in oasis_def_var for '//TRIM(srcv(ji)%clname) )
ENDDO

! Allocate array to store received fields between two coupling steps
ALLOCATE( frcv(hdshape(2),hdshape(4), nfld_rcv_tot), stat = nerror ); frcv(:,:,:) = 0._wp
IF ( nerror > 0 ) THEN
  CALL oasis_abort( ncomp_id, 'oas_hd_define', 'Failure in allocating frcv' )
  RETURN
ENDIF
      

! ------------------------------------------------------------------
! ... End of definition phase (must be called by all processes including 
!     the PEs not involved in the coupling)
! ------------------------------------------------------------------

CALL oasis_enddef( nerror )

IF ( nerror /= OASIS_Success ) CALL oasis_abort ( ncomp_id, &
    'oas_hd_define', 'Failure in oasis_enddef')

if (kl_mype == 0 ) then
 write(0,*) '++++HD: after oasis_enddef, end of oas_hd_define'
endif

! -----------------------------------------------------------------
! ... Deallocate arrays that are not needed anymore
! -----------------------------------------------------------------

IF (ALLOCATED (nlev_snd_oce)) THEN
  DEALLOCATE ( nlev_snd_oce, STAT=nerror )
  DEALLOCATE ( nam_snd_oce , STAT=nerror )
ENDIF
IF (ALLOCATED (nlev_rcv_oce)) THEN
  DEALLOCATE ( nlev_rcv_oce, STAT=nerror )
  DEALLOCATE ( nam_rcv_oce , STAT=nerror )
ENDIF

if (p_nprocs .gt. 1)then
 CALL MPI_Barrier(kl_comm, nerror)
endif

END SUBROUTINE oas_hd_define

!================================================================================
SUBROUTINE hd_receive_fld
!
!!---------------------------------------------------------------------
!!              ***  ROUTINE hd_receive_fld  ***
!!
!! ** Purpose : Prepare and receive coupling fields
!!      from OASIS
!!----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=i4) ::   isec, info, nerror  	 ! temporary integer
INTEGER(KIND=i4) ::   i,j,jn        	  	 ! index
INTEGER(KIND=i4) , DIMENSION(nfld_rcv_tot) :: &
                             nrcvinfo                    ! OASIS info argument
REAL(KIND=wp)           ::   ztmp1(hdshape(2),hdshape(4))  ! temporary array

istep = get_time_step()
isec  = nint(istep * delta_time)

!if ( kl_mype == 0 ) write(0,*) '++++HD: hd_receive_fld: isec = ',isec

!------------------------------------------------------------------------------
! 2. Receive all coupling fields (independent of the specific coupling) 
!------------------------------------------------------------------------------
DO jn = 1, nfld_rcv_tot
  IF( srcv(jn)%laction ) THEN
   CALL oas_hd_rcv( jn, isec, ztmp1(1:hdshape(2), 1:hdshape(4)), nrcvinfo(jn) )
   IF( nrcvinfo(jn) == OASIS_Rcv ) frcv(1:hdshape(2), 1:hdshape(4),jn)=ztmp1(1:hdshape(2), 1:hdshape(4))
  ENDIF
ENDDO

jn = 0
! 1: RUNOFF_S (kg/m2): surface water runoff; sum over forecast
jn = jn + 1
IF( srcv(jn)%laction .AND. nrcvinfo(jn) == OASIS_Rcv ) THEN
  runoff_s(1:hdshape(2), 1:hdshape(4)) = frcv(1:hdshape(2), 1:hdshape(4),jn)
!  DO j = 1, hdshape(4)
!   runoff_s(1:hdshape(2),j) = frcv(1:hdshape(2),hdshape(4)-j+1,jn)
!  ENDDO
ENDIF
! 2: RUNOFF_G (kg/m2): soil water runoff; sum over forecast
jn = jn + 1
IF( srcv(jn)%laction .AND. nrcvinfo(jn) == OASIS_Rcv ) THEN
  runoff_dr(1:hdshape(2), 1:hdshape(4)) = frcv(1:hdshape(2), 1:hdshape(4),jn)
!  DO j = 1, hdshape(4)
!   runoff_dr(1:hdshape(2),j) = frcv(1:hdshape(2),hdshape(4)-j+1,jn)
!  ENDDO
ENDIF

END SUBROUTINE hd_receive_fld

!================================================================================
SUBROUTINE oas_hd_rcv( kid, kstep, pdata, kinfo )

!!---------------------------------------------------------------------
!!              ***  ROUTINE oas_hd_rcv  ***
!!
!! ** Purpose : - At each coupling time-step, this routine call fields
!!                from the coupler.
!!----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(kind=i4), INTENT( IN    )   :: kid    ! variable index
INTEGER(kind=i4), INTENT( IN    )   :: kstep  ! time-step in seconds
REAL(kind=wp), DIMENSION(hdshape(2),hdshape(4)), INTENT( OUT )   :: pdata ! return data
INTEGER(kind=i4), INTENT(   OUT )   :: kinfo  ! OASIS info argument
LOGICAL                   :: llaction                ! local action flag

!!--------------------------------------------------------------------
!
! receive local data from OASIS 
CALL oasis_get ( srcv(kid)%nid, kstep, exfld(1:hdshape(2), 1:hdshape(4)), kinfo )
IF ( nerror .ne. OASIS_Success ) CALL oasis_abort( srcv(kid)%nid, modname,   &
         &                'Failure in oasis_get for '//TRIM(srcv(kid)%clname) )

      llaction = .false.

      IF( kinfo == OASIS_Recvd   .OR. kinfo == OASIS_FromRest .OR.   &
          kinfo == OASIS_RecvOut .OR. kinfo == OASIS_FromRestOut .OR. kinfo == OASIS_Input)   llaction = .TRUE.
       
      IF ( llaction ) THEN

         ! Declare to the calling routine that the coupling field is received from OASIS
         kinfo = OASIS_Rcv

         ! Update the output array which contains the coupling field
         pdata(1:hdshape(2), 1:hdshape(4)) = exfld(1:hdshape(2), 1:hdshape(4))

      ELSE
         ! Declare to the calling routine that the coupling field is NOT received from OASIS
         kinfo = OASIS_idle     
      ENDIF

END SUBROUTINE oas_hd_rcv

!================================================================================
SUBROUTINE hd_send_fld(friv_hd)

!!---------------------------------------------------------------------
!!              ***  ROUTINE hd_send_fld  ***
!!
!! ** Purpose : Prepare and send coupling fields
!!      to OASIS
!!----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(wp), INTENT(IN) :: friv_hd(:,:) ! river flow

  INTEGER(KIND=i4) :: isec, kinfo, jn, ig, i,j	 ! temporary integer
  REAL(KIND=wp)    :: ztmp1(hdshape(2),hdshape(4))	 ! temporary array
  REAL(KIND=wp)    :: rivmouth(hdshape(2),hdshape(4))	 ! temporary array

  ztmp1=0._wp

  istep = get_time_step()  
  isec  = nint(istep * delta_time)

  WHERE (fdir_hd == 0)
   rivmouth(1:hdshape(2), 1:hdshape(4)) = 1.
  ELSEWHERE
   rivmouth(1:hdshape(2), 1:hdshape(4)) = 0.
  END WHERE

  !  DO i = 1, hdshape(2)
  !   DO j = 1, hdshape(4)
  !    IF (rivmouth (i,j) .GT. 0. .and. friv_hd(i,j) .GT. 0. ) THEN
  !     print*, "++++ HD: oas_hd.f90: step=", istep, " rivmouth=",rivmouth(i,j), friv_hd(i,j)
  !    ENDIF
  !   ENDDO
  !  ENDDO


  jn = 0

  ! RDC2NEMO [m3/s]: river discharge to NEMO
  jn = jn + 1
  IF( ssnd(jn)%laction ) THEN
  !   ztmp1(1:hdshape(2), 1:hdshape(4)) = friv_hd(1:hdshape(2), 1:hdshape(4))*rivmouth(1:hdshape(2), 1:hdshape(4))
    DO j = 1, hdshape(4)
     ztmp1(1:hdshape(2),j) = friv_hd(1:hdshape(2),hdshape(4)-j+1)*rivmouth(1:hdshape(2),hdshape(4)-j+1)
    ENDDO
    CALL oas_hd_snd (jn, isec, ztmp1, kinfo )
  ENDIF

END SUBROUTINE hd_send_fld

!================================================================================
SUBROUTINE oas_hd_snd( kid, kstep, pdata, kinfo )

!!---------------------------------------------------------------------
!!              ***  ROUTINE oas_hd_snd  ***
!!
!! ** Purpose : - At each coupling time-step, this routine sends fields
!!                to the coupler.
!!----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(kind=i4), INTENT( IN    )   :: kid    ! variable index
INTEGER(kind=i4), INTENT(   OUT )   :: kinfo  ! OASIS info argument
INTEGER(kind=i4), INTENT( IN    )   :: kstep  ! time-step in seconds
REAL(kind=wp), DIMENSION(hdshape(2),hdshape(4)), INTENT( IN )   :: pdata ! return data

!!
!!--------------------------------------------------------------------
!
! prepare array for OASIS (only computational domain, without nboundlines)
!
exfld(1:hdshape(2), 1:hdshape(4)) = pdata(1:hdshape(2), 1:hdshape(4))

! Call OASIS at each time step,
! but the fields are sent to the other models only at coupling time steps
! (accumulation etc. is controlled in the SMIOC configuration file)
!
CALL oasis_put( ssnd(kid)%nid, kstep, exfld(1:hdshape(2), 1:hdshape(4)), kinfo)

IF ( nerror .ne. OASIS_Success ) CALL oasis_abort( ssnd(kid)%nid, modname,   &
         &      'Failure in oasis_put for '//TRIM(ssnd(kid)%clname) )

END SUBROUTINE oas_hd_snd

!================================================================================
SUBROUTINE hd_send_fld_nemo

!!---------------------------------------------------------------------
!!              ***  ROUTINE hd_send_fld  ***
!!
!! ** Purpose : Prepare and send coupling fields
!!      to OASIS
!!----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=i4) :: isec, kinfo, jn, ig, i,j	 ! temporary integer
REAL(KIND=wp)    :: ztmp1(hdshape_nemo(2),hdshape_nemo(4))	 ! temporary array

ztmp1=0._wp

istep = get_time_step()  
isec  = nint(istep * delta_time)

jn = 0

! RDC2NEMO [m3/s]: river discharge to NEMO
jn = jn + 1
IF( ssnd(jn)%laction ) THEN
 ztmp1(1:hdshape_nemo(2), 1:hdshape_nemo(4)) = discharge_on_ocean(1:hdshape_nemo(2), 1:hdshape_nemo(4))
 CALL oas_hd_snd_nemo (jn, isec, ztmp1, kinfo )
ENDIF

END SUBROUTINE hd_send_fld_nemo

!================================================================================
SUBROUTINE oas_hd_snd_nemo( kid, kstep, pdata, kinfo )

!!---------------------------------------------------------------------
!!              ***  ROUTINE oas_hd_snd  ***
!!
!! ** Purpose : - At each coupling time-step, this routine sends fields
!!                to the coupler.
!!----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(kind=i4), INTENT( IN    )   :: kid    ! variable index
INTEGER(kind=i4), INTENT(   OUT )   :: kinfo  ! OASIS info argument
INTEGER(kind=i4), INTENT( IN    )   :: kstep  ! time-step in seconds
REAL(kind=wp), DIMENSION(hdshape_nemo(2),hdshape_nemo(4)), INTENT( IN )   :: pdata ! return data

!!
!!--------------------------------------------------------------------
!
! prepare array for OASIS (only computational domain, without nboundlines)
!
exfld_nemo(1:hdshape_nemo(2), 1:hdshape_nemo(4)) = pdata(1:hdshape_nemo(2), 1:hdshape_nemo(4))

! Call OASIS at each time step,
! but the fields are sent to the other models only at coupling time steps
! (accumulation etc. is controlled in the SMIOC configuration file)
!
CALL oasis_put( ssnd(kid)%nid, kstep, exfld_nemo(1:hdshape_nemo(2), 1:hdshape_nemo(4)), kinfo)

IF ( nerror .ne. OASIS_Success ) CALL oasis_abort( ssnd(kid)%nid, modname,   &
         &      'Failure in oasis_put for '//TRIM(ssnd(kid)%clname) )

END SUBROUTINE oas_hd_snd_nemo

!================================================================================

SUBROUTINE decomp_def( &
  id_paral, global_extent_lon, global_extent_lat, &
  local_extent_lon, local_extent_lat, local_start_lon, local_start_lat)

  IMPLICIT NONE
  INTEGER, INTENT(out) :: id_paral(OASIS_Box_Params)
  INTEGER, INTENT(IN)  :: global_extent_lon
  INTEGER, INTENT(IN)  :: global_extent_lat
  INTEGER, INTENT(IN)  :: local_extent_lon
  INTEGER, INTENT(IN)  :: local_extent_lat
  INTEGER, INTENT(IN)  :: local_start_lon
  INTEGER, INTENT(IN)  :: local_start_lat

  ! define parameter for the 'box' partition of OASIS

  ! this is a 'box' partition
  id_paral(1) = 2
  ! the upper left corner global offset
  id_paral(2) = (local_start_lon - 1) + &
                (local_start_lat - 1) * global_extent_lon
  ! the local extent in x
  id_paral(3) = local_extent_lon
  ! the local extent in y
  id_paral(4) = local_extent_lat
  ! the global extent in x
  id_paral(5) = global_extent_lon

END SUBROUTINE decomp_def

!================================================================================

END MODULE oas_hd
#endif
! Ha Ho-Hagemann }
