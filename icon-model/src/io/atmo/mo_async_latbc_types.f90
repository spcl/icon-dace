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

! This Module is a description for listings and indices
! used for input prefetching routine.

MODULE mo_async_latbc_types

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, C_F_POINTER
  USE mo_kind,                     ONLY: wp, sp, i8
  USE mo_dictionary,               ONLY: DICT_MAX_STRLEN
  USE mtime,                       ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, &
    &                                    event, datetime, timedelta, OPERATOR(-), &
    &                                    getTotalSecondsTimedelta, deallocateTimedelta, &
    &                                    deallocateEvent, datetimeToString, timedeltaToString
  USE mo_util_mtime,               ONLY: mtime_timedelta_to_seconds
  USE mo_initicon_types,           ONLY: t_init_state, t_init_state_const, t_pi_atm, t_init_state_finalize
  USE mo_impl_constants,           ONLY: SUCCESS, max_ntracer, vname_len
  USE mo_exception,                ONLY: finish, message, message_text
  USE mo_run_config,               ONLY: msg_level
  USE mo_reorder_info,             ONLY: t_reorder_info, release_reorder_info
  USE mo_mpi,                      ONLY: p_comm_work_pref, p_barrier
  USE mo_cdi,                      ONLY: cdi_undefid, streamclose
#ifndef NOMPI
  USE mpi
#endif
  IMPLICIT NONE

  PRIVATE

  ! module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc_types'

  ! derived data types:
  PUBLIC :: t_latbc_state
  PUBLIC :: t_latbc_data
  PUBLIC :: t_patch_data
  PUBLIC :: t_mem_win
  PUBLIC :: t_buffer
  PUBLIC :: t_size
  PUBLIC :: t_glb_indices


  !------------------------------------------------------------------------------------------------
  ! DERIVED DATA TYPES
  !------------------------------------------------------------------------------------------------


#ifdef NOMPI
  INTEGER, PARAMETER :: mpi_win_null = 0
#endif

  !> Data structure containing variables for MPI memory window
  !
  TYPE t_mem_win
     ! Currently, we use only 1 MPI window for all input prefetching
     ! Used for async prefetch only
#ifndef NOMPI
     INTEGER                   :: mpi_win = mpi_win_null
     INTEGER(mpi_address_kind) :: f_mem_ptr
#endif
     REAL(sp), POINTER  :: mem_ptr_sp(:) => NULL() !< Pointer to memory window (REAL*4)
  END TYPE t_mem_win


  TYPE t_size
     REAL(sp), POINTER :: buffer(:,:,:) => NULL()
  END TYPE t_size


  !> Data structure containing the lateral boundary forcing data.
  !> Compared to t_init_state it allows to store forcing data 
  !< for additional child domains, if required. 
  !< This is required in case of upper boundary nudging.
  !
  TYPE, EXTENDS (t_init_state) :: t_latbc_state

    TYPE (t_pi_atm), ALLOCATABLE  :: atm_child(:)  ! atmospheric state for child domains
                                                   ! (excluding the base domain)
   CONTAINS
     PROCEDURE :: finalize => t_latbc_state_finalize   !< destructor
  END TYPE t_latbc_state


  TYPE t_buffer
     INTEGER                                     :: ngrp_vars          ! Number of variables for prefetching
     CHARACTER(LEN=DICT_MAX_STRLEN), ALLOCATABLE :: mapped_name(:)     ! name of mapped dictionary variables
     CHARACTER(LEN=DICT_MAX_STRLEN), ALLOCATABLE :: internal_name(:)   ! corresponding internal name of variables
     INTEGER,                        ALLOCATABLE :: nlev(:)            ! Size of variables for prefetching
     TYPE(t_size),                   ALLOCATABLE :: vars(:)
     INTEGER,                        ALLOCATABLE :: varID(:)           ! ID for variable to be read from file
     INTEGER,                        ALLOCATABLE :: hgrid(:)           ! CDI horizontal grid type (cell/edge)
     LOGICAL                                     :: lread_qr, lread_qs ! are qr, qs provided as input?

     ! for additional tracer variables, e.g. ART tracers
     LOGICAL                                     :: lread_tracer(max_ntracer) ! provided as input?
     CHARACTER(LEN=vname_len)                  :: name_tracer(max_ntracer)  ! names
     INTEGER                                     :: idx_tracer(max_ntracer)   ! indices in tracer container

     LOGICAL                                     :: lread_vn           ! is vn provided as input?
     LOGICAL                                     :: lread_u_v          ! is u,v provided as input?

     ! If .FALSE., input levels are computed from sfc geopotential:
     LOGICAL                                     :: lread_hhl

     ! are prognostic thermodynamic variables (= rho and theta_v) present in the input file?
     LOGICAL                                     :: lread_theta_rho      

     ! .TRUE., if pressure is read from input
     LOGICAL                                     :: lread_pres

     ! .TRUE., if temperature is read from input
     LOGICAL                                     :: lread_temp

     ! If .TRUE., surface pressure and geopotential are available in
     ! the input file
     LOGICAL                                     :: lread_ps_geop

     ! .TRUE., if vertical wind or omega is available in input data
     LOGICAL                                     :: lread_w

     ! .FALSE., if vertical component of velocity (W) is provided as input
     LOGICAL                                     :: lconvert_omega2w

     ! .TRUE., if heights are computed (hydrostatic model input):
     LOGICAL                                     :: lcompute_hhl_pres

     CHARACTER(LEN=10)                           :: geop_ml_var        ! model level surface geopotential
     CHARACTER(LEN=10)                           :: hhl_var

     ! input data validity DateTime
     TYPE(datetime)                              :: vDateTime          

   CONTAINS
     PROCEDURE :: finalize => t_buffer_finalize   !< destructor
  END TYPE t_buffer


  ! TYPE p_patch_info contains the ordering info for cells, edges and verts
  TYPE t_patch_data
     TYPE(t_reorder_info) :: cells
     TYPE(t_reorder_info) :: edges
     LOGICAL, ALLOCATABLE :: cell_mask(:,:), edge_mask(:,:)

     ! used for async prefetching only
     TYPE(t_mem_win) :: mem_win  !< data structure containing variables for MPI memory window

     ! number of full and half levels
     INTEGER :: nlev
     INTEGER :: nlevp1

     INTEGER :: level     ! patch level (e.g. xx in R03Bxx)
     INTEGER :: num_vars  ! no of input prefetch variables

     ! number of cells and edges in the local patch
     INTEGER :: n_patch_cells
     INTEGER :: n_patch_edges

     ! number of cells and edges in the global patch
     INTEGER :: n_patch_cells_g
     INTEGER :: n_patch_edges_g

     ! number of blocks
     INTEGER :: nblks_c, nblks_e

   CONTAINS
     PROCEDURE :: finalize => t_patch_data_finalize   !< destructor
  END TYPE t_patch_data


  !> Sparse latbc mode: index data for boundary rows:
  !  Derived type specifying a local-to-global index mapping for
  !  extracted subgrids.
  !
  TYPE t_glb_indices
    INTEGER, ALLOCATABLE :: cells(:), edges(:)      !< (1...local) global indices for cells and edges
    INTEGER              :: n_patch_cells_g         !< total no. of global cells
    INTEGER              :: n_patch_edges_g         !< total no. of global edges
  CONTAINS
    PROCEDURE :: finalize => t_glb_indices_finalize
  END TYPE t_glb_indices


  !> Data type containing the time control, the variable buffer, and
  !> the necessary index arrays for reordering, to read in lateral
  !> boundary data into ICON.
  !
  TYPE t_latbc_data

    !> full path of currently open file, unallocated if not open
    CHARACTER(:), ALLOCATABLE :: open_filepath
    !> CDI handle of currently open stream, CDI_UNDEFID if not open
    INTEGER :: open_cdi_stream_handle = cdi_undefid
    TYPE(datetime) :: mtime_last_read
    TYPE(event),     POINTER :: prefetchEvent   => NULL()
    TYPE(timedelta), POINTER :: delta_dtime     => NULL()

    ! time level indices for  latbc_data. can be 1 or 2.
    INTEGER :: new_latbc_tlev

    ! storage for time-constant height level data
    TYPE(t_init_state_const), POINTER :: latbc_data_const => NULL()

    ! storage for two time-level boundary data
    TYPE(t_latbc_state) :: latbc_data(2)

    ! raw buffer
    TYPE(t_buffer) :: buffer

    !< indices for async latbc prefetching
    TYPE(t_patch_data) :: patch_data

    ! for sparse latbc mode: index data for boundary rows:
    TYPE(t_glb_indices) :: global_index

    ! linear interpolation weights for time interpolation between 
    ! consecutive boundary data
    REAL(wp):: lc1, lc2
  CONTAINS
    PROCEDURE :: finalize        => t_latbc_data_finalize
    PROCEDURE :: prev_latbc_tlev => prev_latbc_tlev
    PROCEDURE :: update_intp_wgt => t_latbc_data_update_intp_wgt

  END TYPE t_latbc_data


CONTAINS

  !> Destructor for "t_glb_indices" class.
  !
  SUBROUTINE t_glb_indices_finalize(this)
    CLASS(t_glb_indices) :: this
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::t_glb_indices_finalize'
    INTEGER :: ierrstat=0

    !CALL message(routine, 't_glb_indices_finalize')

    IF (ALLOCATED(this%cells))  DEALLOCATE(this%cells, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    IF (ALLOCATED(this%edges))  DEALLOCATE(this%edges, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE t_glb_indices_finalize


  SUBROUTINE t_patch_data_finalize(patch_data)
    CLASS(t_patch_data), INTENT(INOUT) :: patch_data

#ifndef NOMPI
    CHARACTER(len=*), PARAMETER :: routine = modname//'::t_patch_data_finalize'
    LOGICAL, PARAMETER :: lprint_dbg = .FALSE.
    INTEGER            :: ierror
    TYPE(c_ptr)        :: c_mem_ptr
    INTEGER, POINTER   :: baseptr

    CALL message(routine, "")
    IF (lprint_dbg) CALL p_barrier(comm=p_comm_work_pref) ! make sure all are here
#endif
    CALL release_reorder_info(patch_data%cells)
    IF (ALLOCATED(patch_data%cell_mask)) DEALLOCATE(patch_data%cell_mask)
    CALL release_reorder_info(patch_data%edges)
    IF (ALLOCATED(patch_data%edge_mask)) DEALLOCATE(patch_data%edge_mask)
#ifndef NOMPI
    ! note: we do not touch the MPI window pointer here:
    !
    IF (lprint_dbg) CALL p_barrier(comm=p_comm_work_pref) ! make sure all are here
    IF  ((msg_level >= 15) .OR. lprint_dbg)  CALL message(routine, "Free MPI window")
    IF (patch_data%mem_win%mpi_win /= mpi_win_null) THEN
      CALL mpi_win_free(patch_data%mem_win%mpi_win, ierror)
      IF (ierror /= 0) CALL finish(routine, "mpi_win_free failed!")
    END IF
    IF (lprint_dbg) CALL p_barrier(comm=p_comm_work_pref) ! make sure all are here
    IF  ((msg_level >= 15) .OR. lprint_dbg)  CALL message(routine, "Nullify MPI window pointer")
    IF (ASSOCIATED(patch_data%mem_win%mem_ptr_sp)) THEN
      NULLIFY(patch_data%mem_win%mem_ptr_sp)
    END IF
    IF (lprint_dbg) CALL p_barrier(comm=p_comm_work_pref) ! make sure all are here
    IF  ((msg_level >= 15) .OR. lprint_dbg)  CALL message(routine, "done.")
#endif
  END SUBROUTINE t_patch_data_finalize


  SUBROUTINE t_buffer_finalize(buffer)
    CLASS(t_buffer), INTENT(INOUT) :: buffer
    INTEGER :: i

    !CALL message("", 't_buffer_finalize')

    IF (ALLOCATED(buffer%mapped_name))    DEALLOCATE(buffer%mapped_name)
    IF (ALLOCATED(buffer%internal_name))  DEALLOCATE(buffer%internal_name)
    IF (ALLOCATED(buffer%nlev))           DEALLOCATE(buffer%nlev)
    IF (ALLOCATED(buffer%vars)) THEN
      DO i=1,SIZE(buffer%vars)
        IF (ASSOCIATED(buffer%vars(i)%buffer))  DEALLOCATE(buffer%vars(i)%buffer)      
      END DO
      DEALLOCATE(buffer%vars)
    END IF
    IF (ALLOCATED(buffer%varID))          DEALLOCATE(buffer%varID)
    IF (ALLOCATED(buffer%hgrid))          DEALLOCATE(buffer%hgrid)
  END SUBROUTINE t_buffer_finalize


  SUBROUTINE t_latbc_state_finalize(init_data)
    CLASS(t_latbc_state), INTENT(INOUT) :: init_data
    INTEGER :: jg 

    CALL t_init_state_finalize(init_data)
    !
    IF (ALLOCATED(init_data%atm_child)) THEN
      DO jg=1, SIZE(init_data%atm_child)
        CALL init_data%atm_child(jg)%finalize
      ENDDO
      DEALLOCATE(init_data%atm_child)
    ENDIF

  END SUBROUTINE t_latbc_state_finalize

  !-------------------------------------------------------------------------
  !!
  SUBROUTINE t_latbc_data_finalize(latbc)
    CLASS(t_latbc_data), INTENT(INOUT) :: latbc

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::t_latbc_data_finalize"
    INTEGER :: tlev, ierror

    CALL message(routine, 'deallocating latbc data')

    ! deallocate boundary data memory
    DO tlev = 1, 2
      CALL latbc%latbc_data(tlev)%finalize()
    END DO
    
    CALL latbc%patch_data%finalize()          ! deallocate patch data
    CALL latbc%buffer%finalize()              ! deallocate intermediate storage latbc%buffer
    CALL latbc%global_index%finalize()        ! clean up global indices data structure.

    IF (ASSOCIATED(latbc%prefetchEvent)) THEN
      CALL deallocateEvent(latbc%prefetchEvent) ! deallocate prefetch input event
    END IF

    ! deallocating date and time data structures
    IF (ASSOCIATED(latbc%delta_dtime)) &
      CALL deallocateTimedelta(latbc%delta_dtime)
    IF (ASSOCIATED(latbc%latbc_data_const)) THEN
      DEALLOCATE(latbc%latbc_data_const, stat=ierror)
      IF (ierror /= SUCCESS) CALL finish(routine, "deallocate failed!")
    END IF
    ! close input file if open
    IF (latbc%open_cdi_stream_handle /= cdi_undefid) THEN
      CALL streamclose(latbc%open_cdi_stream_handle)
      latbc%open_cdi_stream_handle = cdi_undefid
    END IF
    IF  (msg_level >= 15)  CALL message(routine, 'done.')
  END SUBROUTINE t_latbc_data_finalize


  INTEGER FUNCTION prev_latbc_tlev(latbc)
    CLASS(t_latbc_data), INTENT(IN) :: latbc

    prev_latbc_tlev = 3 - latbc%new_latbc_tlev
  END FUNCTION prev_latbc_tlev


    !-------------------------------------------------------------------------
    !  Update linear interpolation weights (lc1, lc2) for the time interpolation 
    !  between two consecutive boundary forcing time slices.
    !
    !
    !      <------------------------------------->  dtime_latbc
    !                               <------------>  delta_tstep
    !                          datetime_current   
    !   ---*------------------------X------------*----------------> time axis
    !    vDate_prev                            vDate_cur
    !
    SUBROUTINE t_latbc_data_update_intp_wgt( latbc, datetime_current )
      CLASS(t_latbc_data), INTENT(INOUT)  :: latbc
      TYPE(datetime),      INTENT(IN)     :: datetime_current ! datetime for which interpolation 
                                                              ! weights shall be computed
#ifndef NOMPI
      TYPE(timedelta)         :: delta_tstep        ! time delta between current datetime and 
                                                    ! validity time of current boundary forcing time slice.
      INTEGER(i8)             :: delta_tstep_in_sec ! delta_tstep converted to seconds

      REAL(wp)                :: dtime_latbc        ! time delta between two consecutive 
                                                    ! boundary forcing time slices [s]
      TYPE(timedelta)         :: td                 ! same in mtime format
      LOGICAL                 :: failure
      INTEGER                 :: prv_tlev, cur_tlev

      CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::t_latbc_data_update_intp_wgt"
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: vDateTime_str_cur, vDateTime_str_prv
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN):: delta_tstep_str


      ! compute boundary update timedelta
      cur_tlev = latbc%new_latbc_tlev
      prv_tlev = latbc%prev_latbc_tlev()
      td = latbc%latbc_data(cur_tlev)%vDateTime - latbc%latbc_data(prv_tlev)%vDateTime
      dtime_latbc =  REAL(getTotalSecondsTimedelta(td, latbc%latbc_data(cur_tlev)%vDateTime))

      delta_tstep = latbc%latbc_data(cur_tlev)%vDateTime - datetime_current

      failure = delta_tstep%month /= 0
      IF (msg_level >= 15 .OR. failure) THEN
        CALL message(routine, "", all_print=failure)
        CALL datetimeToString(latbc%latbc_data(prv_tlev)%vDateTime, vDateTime_str_prv)
        CALL datetimeToString(latbc%latbc_data(cur_tlev)%vDateTime, vDateTime_str_cur)
        CALL timedeltaToString(delta_tstep, delta_tstep_str)
        WRITE (message_text, '(a,a)')  "lbc vdate current : ", vDateTime_str_cur
        CALL message("", message_text, all_print=failure)
        WRITE (message_text, '(a,a)')  "lbc vdate previous: ", vDateTime_str_prv
        CALL message("", message_text, all_print=failure)
        WRITE (message_text, '(a,a)')  "delta_tstep_str: ", delta_tstep_str
        CALL message("", message_text, all_print=failure)
        IF (failure) &
          CALL finish(routine, "time difference for reading boundary&
          & data must not be more than a month.")
      ENDIF

      ! compute the number of "dtime_latbc" intervals fitting into the time difference "delta_tstep":

      CALL mtime_timedelta_to_seconds(delta_tstep, delta_tstep_in_sec)
      !
      latbc%lc1 = REAL(delta_tstep_in_sec,wp) / dtime_latbc
      latbc%lc2 = 1._wp - latbc%lc1

#endif

    END SUBROUTINE t_latbc_data_update_intp_wgt


END MODULE mo_async_latbc_types

!------------------------------------------------------------------------------------------------
