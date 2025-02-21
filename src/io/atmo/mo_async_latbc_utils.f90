! This module contains the asynchronous I/O routine for lateral boundary nudging
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

!----------------------------
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

!----------------------------

  MODULE mo_async_latbc_utils




    USE mo_async_latbc_types,   ONLY: t_latbc_state, t_latbc_data, t_buffer
    USE mo_reorder_info,        ONLY: t_reorder_info
    USE mo_kind,                ONLY: wp, i8
    USE mo_util_string,         ONLY: int2string
    USE mo_parallel_config,     ONLY: nproma, proc0_offloading
    USE mo_model_domain,        ONLY: t_patch
    USE mo_grid_config,         ONLY: nroot, n_dom
    USE mo_exception,           ONLY: message, finish, message_text
    USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS, min_rlcell
    USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
    USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
    USE mo_io_units,            ONLY: filename_max
    USE mo_nonhydro_types,      ONLY: t_nh_state
    USE mo_intp_data_strc,      ONLY: t_int_state
    USE mo_nh_vert_interp,      ONLY: vert_interp
    USE mo_nh_nest_utilities,   ONLY: intp_nestubc_nudging
    USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1
    USE mo_nh_init_utils,       ONLY: convert_omega2w, compute_input_pressure_and_height
    USE mo_sync,                ONLY: sync_patch_array, SYNC_E
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
    USE mtime,                  ONLY: timedelta, newTimedelta, deallocateTimedelta, &
         &                            newEvent, datetime, newDatetime,             &
         &                            isCurrentEventActive, deallocateDatetime,    &
         &                            MAX_DATETIME_STR_LEN,                        &
         &                            getTotalSecondsTimedelta,                    &
         &                            datetimeToString,                            &
         &                            OPERATOR(>=), OPERATOR(-), OPERATOR(>),      &
         &                            OPERATOR(/=), OPERATOR(+), OPERATOR(*),      &
         &                            OPERATOR(<), OPERATOR(==),                   &
         &                            getTriggerNextEventAtDateTime
    USE mo_util_mtime,          ONLY: dummyDateTime
    USE mo_time_config,         ONLY: time_config
    USE mo_limarea_config,      ONLY: latbc_config, generate_filename
    USE mo_nudging_config,      ONLY: nudging_config, indg_type, ithermdyn_type
    USE mo_initicon_config,     ONLY: timeshift
    USE mo_ext_data_types,      ONLY: t_external_data
    USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, ltransport, msg_level, ntracer
    USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
    USE mo_cdi,                 ONLY: streamOpenRead, streamClose, streamInqVlist, vlistInqTaxis, &
      &                               taxisInqVDate, taxisInqVTime, &
      &                               cdiDecodeTime, cdiDecodeDate, &
      &                               cdi_undefid
    USE mo_util_cdi,            ONLY: cdiGetStringError, read_cdi_2d, read_cdi_3d, t_inputParameters,  &
    &                                 makeInputParameters, deleteInputParameters
    USE mo_util_file,           ONLY: util_filesize
    USE mo_master_config,       ONLY: isRestart
    USE mo_fortran_tools,       ONLY: copy, init
    USE mo_util_string,         ONLY: tolower
    USE mo_util_sysinfo,        ONLY: check_file_exists
    USE mo_dictionary,          ONLY: t_dictionary




    IMPLICIT NONE
    PRIVATE

    ! handshake subroutines
    PUBLIC :: async_pref_send_handshake
    PUBLIC :: compute_wait_for_async_pref
    PUBLIC :: async_pref_wait_for_start
    PUBLIC :: compute_shutdown_async_pref
    PUBLIC :: allocate_pref_latbc_data





    PUBLIC ::  async_init_latbc_data, prefetch_latbc_data,     &
         &     recv_latbc_data


    INTERFACE fetch_from_buffer
      MODULE PROCEDURE fetch_from_buffer_2D
      MODULE PROCEDURE fetch_from_buffer_3D_cells 
      MODULE PROCEDURE fetch_from_buffer_3D_generic
    END INTERFACE

    INTERFACE get_data
      MODULE PROCEDURE get_data_2D
      MODULE PROCEDURE get_data_3D 
    END INTERFACE

    TYPE t_read_params
      TYPE(t_inputParameters) :: cdi_params
      INTEGER                 :: npoints = 0
      INTEGER                 :: imode_asy
      INTEGER, POINTER        :: idx_ptr(:) => NULL()
    END TYPE t_read_params


    !------------------------------------------------------------------------------------------------
    ! CONSTANTS
    !------------------------------------------------------------------------------------------------

    ! module name string
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc_utils'

    ! Tags for communication between compute PEs and prefetching PEs
    INTEGER, PARAMETER :: msg_pref_start    = 56984
    INTEGER, PARAMETER :: msg_pref_done     = 26884
    INTEGER, PARAMETER :: msg_pref_shutdown = 48965

    INTEGER, PARAMETER :: TAG_PREFETCH2WORK = 2001
    INTEGER, PARAMETER :: TAG_WORK2PREFETCH = 2002
    INTEGER, PARAMETER :: TAG_VDATETIME     = 2003

    INTEGER, PARAMETER :: icell             = 1
    INTEGER, PARAMETER :: iedge             = 2

  CONTAINS


    !-------------------------------------------------------------------------
    !!
    SUBROUTINE allocate_pref_latbc_data(latbc, nlev_in, p_nh_state, ext_data, p_patch)

      TYPE(t_latbc_data), TARGET, INTENT(INOUT) :: latbc       !< variable buffer for latbc data
      INTEGER,                    INTENT(IN)    :: nlev_in     !< no. of vertical input levels
      TYPE(t_nh_state),           INTENT(INOUT) :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_external_data),      INTENT(IN)    :: ext_data    !< external data on the global domain
      TYPE(t_patch),              INTENT(IN)    :: p_patch(:)

    END SUBROUTINE allocate_pref_latbc_data




    !-------------------------------------------------------------------------
    !! Read interpolated lateral boundary conditions
    !!
    !! This subroutine is called by compute processors.
    !! Depending on the parameter read_params%imode_asy, reading is done either synchronously
    !! be PE0 (for the initial lateral boundary data), or asynchronously by the prefetch PE.
    !! In the latter case, the 'get_data' routine copies the data from the memory buffer
    !!
    !! In the final step, the data are interpolated vertically from intermediate
    !! remapicon grid to ICON grid, followed by computing/completing the prognostic NH variable set
    !!
    !! NOTE: This subroutine is MPI-collective, since
    !!       it contains several synchronization calls. It must
    !!       therefore be passed by all worker PEs. However, there is
    !!       the common situation where vertical interpolation shall
    !!       performed on a subset of PEs only, while no valid data is
    !!       available on the remaining PE. For this situation we need
    !!       the optional "opt_lmask" parameter.
    !!
    SUBROUTINE read_latbc_data(latbc, p_patch, p_nh_state, p_int, tlev, read_params, latbc_dict)
      TYPE(t_latbc_data),     INTENT(INOUT), TARGET :: latbc  !< variable buffer for latbc data
      TYPE(t_patch),          INTENT(INOUT)         :: p_patch
      TYPE(t_nh_state),       INTENT(IN)            :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      INTENT(IN)            :: p_int
      INTEGER,                INTENT(IN)            :: tlev
      TYPE(t_read_params),    INTENT(INOUT)         :: read_params(:)
      TYPE (t_dictionary),    INTENT(IN), OPTIONAL  :: latbc_dict

    END SUBROUTINE read_latbc_data

    !-------------------------------------------------------------------------
    !!
    !! ** This subroutine is only called by the prefetching PEs. **
    !!
    SUBROUTINE async_init_latbc_data(latbc)
      TYPE(t_latbc_data), INTENT(INOUT) :: latbc

    END SUBROUTINE async_init_latbc_data


    !-------------------------------------------------------------------------
    !! Read horizontally interpolated atmospheric boundary data.
    !!
    !! The subroutine reads atmospheric boundary data and projects on
    !! the ICON global grid. 
    !!
    !! The following steps are performed:
    !! - Read atmospheric input data,
    !! - Write input data to memory window buffer. The offset for data
    !!   is set such that each of dataset belongs to the respective
    !!   compute processor,
    !!
    !! ** This subroutine is only called by the prefetching PE. **
    !!

    SUBROUTINE prefetch_latbc_data(latbc, latbc_read_datetime)

      TYPE(t_latbc_data), TARGET, INTENT(INOUT)    :: latbc

      ! datetime of the next input time level
      TYPE(datetime), INTENT(IN)         :: latbc_read_datetime

    END SUBROUTINE prefetch_latbc_data

    !-------------------------------------------------------------------------
    !>
    !! Consistency check: Make sure that the requested date is actually contained in the file.
    !!                    Print the file name of the boundary file which will be read.
    !!
    SUBROUTINE check_validity_date_and_print_filename(latbc, latbc_read_datetime,  &
      &                                               mtime_vdate)

      TYPE(t_latbc_data), INTENT(in) :: latbc !< latbc state
      TYPE(datetime), INTENT(IN) :: &
        &  latbc_read_datetime     !< Requested datetime of LatBC file
      TYPE(datetime), POINTER, INTENT(INOUT),OPTIONAL :: &
        &  mtime_vdate             !< LatBC file validity date as mtime object
      ! Local variables
      TYPE(datetime), POINTER :: &
        &  mtime_vdate_loc         !< LatBC file validity date as mtime object
      INTEGER                 :: &
        &  vlistID, taxisID,     & !< CDI identifiers for variables list and time axis
        &  idate, iyear,         & !< Integer value of validity date: total date and year
        &  imonth, iday,         & !< Integer value of validity date: month and day
        &  itime, ihour,         & !< Integer value of validity time: total time and hour
        &  iminute, isecond        !< Integer value of validity time: minute and second
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: &
        &  dstringA,             & !< Requested date (used to generate filename) as string
        &  dstringB                !< File validity date as string
      CHARACTER(LEN=*),PARAMETER :: &
        &  routine = modname//"::check_validity_date_and_print_filename"

      vlistID = streamInqVlist(latbc%open_cdi_stream_handle)
      taxisID = vlistInqTaxis(vlistID)
      idate   = taxisInqVDate(taxisID)
      itime   = taxisInqVTime(taxisID)
      CALL cdiDecodeDate(idate, iyear, imonth, iday)
      CALL cdiDecodeTime(itime, ihour, iminute, isecond)
      mtime_vdate_loc => newDatetime(iyear, imonth, iday, ihour, iminute, isecond, 0)
      IF (msg_level >= 10) THEN
        CALL datetimeToString(latbc_read_datetime, dstringA)
        WRITE (message_text, '(5 A)') routine, &
          ":: reading boundary data from file ", latbc%open_filepath, &
          " for date: ", TRIM(dstringA)
        WRITE (0,*) TRIM(message_text)
      END IF
      IF (mtime_vdate_loc /= latbc_read_datetime) THEN
        CALL datetimeToString(latbc_read_datetime, dstringA)
        CALL datetimeToString(mtime_vdate_loc, dstringB)
        WRITE (message_text, '(6 A)')  "File validity date ", TRIM(dstringB), &
          " of file ", latbc%open_filepath, " does not match requested date ", &
          TRIM(dstringA)
        CALL finish(routine, message_text)
      END IF

      IF (PRESENT(mtime_vdate)) THEN
        mtime_vdate => mtime_vdate_loc
      ELSE
        CALL deallocateDatetime(mtime_vdate_loc)
      END IF
    END SUBROUTINE check_validity_date_and_print_filename

    !-------------------------------------------------------------------------
    !! Receive horizontally interpolated atmospheric boundary data
    !! from the prefetching PE. 
    !!
    !! ** This subroutine is only called by the worker PEs. **
    !!
    SUBROUTINE recv_latbc_data(latbc, p_patch, p_nh_state, p_int, cur_datetime, &
      &                        latbc_read_datetime, lcheck_read, tlev)

      TYPE(t_latbc_data),     INTENT(INOUT) :: latbc
      TYPE(t_patch),          INTENT(INOUT) :: p_patch(:)
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state   !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      INTENT(IN)    :: p_int
      TYPE(datetime), POINTER, INTENT(IN)   :: cur_datetime !< current time
      INTEGER,                INTENT(INOUT) :: tlev

      ! datetime of the next input time level
      TYPE(datetime),         INTENT(IN)    :: latbc_read_datetime

      ! flag to set wether compute processor need to read boundary
      ! data or need not and than they will return
      LOGICAL,      INTENT(IN)    :: lcheck_read


    END SUBROUTINE recv_latbc_data


    ! Wrapper routine for copying prognostic variables from initial state to the 
    ! first time level of the lateral boundary data
    !
    SUBROUTINE copy_fg_to_latbc(latbc_data, p_nh, tlev, idx_tracer)
      TYPE(t_latbc_state),    INTENT(INOUT) :: latbc_data(:)
      TYPE(t_nh_state),       INTENT(IN)    :: p_nh
      INTEGER,                INTENT(IN)    :: tlev
      INTEGER,                INTENT(IN)    :: idx_tracer(:)

      INTEGER, PARAMETER :: jg = 1
      INTEGER            :: idx

!$OMP PARALLEL
      CALL copy(p_nh%prog(nnow(jg))%vn,      latbc_data(tlev)%atm%vn)
      CALL copy(p_nh%prog(nnow(jg))%w,       latbc_data(tlev)%atm%w)
      CALL copy(p_nh%prog(nnow(jg))%rho,     latbc_data(tlev)%atm%rho)
      CALL copy(p_nh%prog(nnow(jg))%theta_v, latbc_data(tlev)%atm%theta_v)
      CALL copy(p_nh%diag%pres,              latbc_data(tlev)%atm%pres)
      CALL copy(p_nh%diag%temp,              latbc_data(tlev)%atm%temp)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqv), latbc_data(tlev)%atm%qv)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqc), latbc_data(tlev)%atm%qc)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqi), latbc_data(tlev)%atm%qi)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqr), latbc_data(tlev)%atm%qr)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqs), latbc_data(tlev)%atm%qs)

      ! Copy additional tracer variables
      DO idx=1, ntracer
        IF (ASSOCIATED(latbc_data(tlev)%atm%tracer(idx)%field)) THEN
          CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,idx_tracer(idx)), &
            &       latbc_data(tlev)%atm%tracer(idx)%field)
        ENDIF
      END DO
!$OMP END PARALLEL

    END SUBROUTINE copy_fg_to_latbc


    !-------------------------------------------------------------------------
    !!
    SUBROUTINE compute_boundary_tendencies ( latbc_data, p_patch, p_nh, tlev, idx_tracer )
      TYPE(t_latbc_state),     INTENT(IN)   :: latbc_data(:)
      TYPE(t_patch),          INTENT(IN)    :: p_patch
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh
      INTEGER,                INTENT(IN)    :: tlev
      INTEGER,                INTENT(IN)    :: idx_tracer(:)


    END SUBROUTINE compute_boundary_tendencies

    !-------------------------------------------------------------------------------------------------
    !> Send a message to the work PEs that the input prefetching PEs is ready. The
    !  counterpart on the work PEs side is compute_wait_for_async_pref
    !
    SUBROUTINE async_pref_send_handshake()
    END SUBROUTINE async_pref_send_handshake

    !-------------------------------------------------------------------------------------------------
    !> compute_wait_for_async_pref: Wait for a message that the prefetch PE is ready.
    !  The counterpart on the input Prefetching PEs side is async_pref_send_handshake
    !
    SUBROUTINE compute_wait_for_async_pref()
    END SUBROUTINE compute_wait_for_async_pref


    !-------------------------------------------------------------------------------------------------
    !> async_pref_wait_for_start: Wait for a message from compute PE that we should start
    !  tranferring the prefetch data or finish. The counterpart on the compute side is
    !  compute_start_async_pref/compute_shutdown_async_pref
    !
    SUBROUTINE async_pref_wait_for_start(done)
      LOGICAL, INTENT(INOUT) :: done ! flag if we should shut down

    END SUBROUTINE async_pref_wait_for_start

    !-------------------------------------------------------------------------------------------------
    !> compute_start_async_pref: Send a message to prefetching PEs that they should start
    !  prefetching input. The counterpart on the prefetching side is async_pref_wait_for_start
    !
    SUBROUTINE compute_start_async_pref()
    END SUBROUTINE compute_start_async_pref

    !-------------------------------------------------------------------------------------------------
    !> compute_shutdown_async_pref: Send a message to prefetching PEs that they should shut down
    !  The counterpart on the prefetching side is async_pref_wait_for_start
    !
    SUBROUTINE compute_shutdown_async_pref
    END SUBROUTINE compute_shutdown_async_pref

    !-------------------------------------------------------------------------
    !>
    ! Return the index for a given variable in mapped variable list
    !
    FUNCTION get_field_index(buffer,name) RESULT(result_varID)
      TYPE(t_buffer), INTENT(IN) :: buffer
      CHARACTER (LEN=*),   INTENT(IN) :: name !< variable name
      ! local variables
      LOGICAL, PARAMETER :: ldebug = .FALSE.
      INTEGER :: result_varID, varID, nvars
      CHARACTER(len=len(name)) :: name_lc

      result_varID = -1
      nvars = buffer%ngrp_vars
      if (nvars >= 1) name_lc = tolower(name)
      ! looping over variable list in internal name
      LOOP : DO varID=1, nvars
        IF (name_lc == tolower(buffer%internal_name(varID))) THEN
          result_varID = varID
          EXIT LOOP
        END IF
      END DO LOOP
    END FUNCTION get_field_index


    ! Wrapper routines for reading data either synchronously by PE0 via read_cdi or asynchronously
    ! by fetching them from the prefetch PE
    !
    SUBROUTINE get_data_3D(latbc, name, arr3d, read_params, opt_latbc_dict)

      TYPE(t_latbc_data),  INTENT(IN)    :: latbc
      CHARACTER(LEN=*),    INTENT(IN)    :: name
      REAL(wp),            INTENT(INOUT) :: arr3d(:,:,:)
      TYPE(t_read_params), INTENT(INOUT) :: read_params

      TYPE (t_dictionary), INTENT(IN),    OPTIONAL :: opt_latbc_dict

      ! local variables
      CHARACTER(LEN=MAX_CHAR_LENGTH) :: mapped_name
      INTEGER                        :: nlev

      IF (PRESENT(opt_latbc_dict)) THEN
        mapped_name = opt_latbc_dict%get(name, default=name)
      ELSE
        mapped_name = name
      ENDIF
      nlev = SIZE(arr3d,2)
      
      IF (read_params%imode_asy == 0) THEN
        CALL read_cdi_3d(read_params%cdi_params, TRIM(mapped_name), nlev, arr3d, read_params%npoints, read_params%idx_ptr)
      ELSE IF (read_params%imode_asy == iedge) THEN
        CALL fetch_from_buffer(latbc, TRIM(name), arr3d, latbc%patch_data%edges)
      ELSE
        CALL fetch_from_buffer(latbc, TRIM(name), arr3d)
      ENDIF

    END SUBROUTINE get_data_3D


    SUBROUTINE get_data_2D(latbc, name, arr2d, read_params, opt_latbc_dict)

      TYPE(t_latbc_data),  INTENT(IN)    :: latbc
      CHARACTER(LEN=*),    INTENT(IN)    :: name
      REAL(wp),            INTENT(INOUT) :: arr2d(:,:)
      TYPE(t_read_params), INTENT(INOUT) :: read_params

      TYPE (t_dictionary), INTENT(IN),    OPTIONAL :: opt_latbc_dict

      ! local variables
      CHARACTER(LEN=MAX_CHAR_LENGTH) :: mapped_name

      IF (PRESENT(opt_latbc_dict)) THEN
        mapped_name = opt_latbc_dict%get(name, default=name)
      ELSE
        mapped_name = name
      ENDIF

      IF (read_params%imode_asy == 0) THEN
        CALL read_cdi_2d(read_params%cdi_params, TRIM(mapped_name), arr2d, read_params%npoints, read_params%idx_ptr)
      ELSE
        CALL fetch_from_buffer(latbc, TRIM(name), arr2d)
      ENDIF

    END SUBROUTINE get_data_2D

    ! ----------------------------------------------------------------------
    ! Auxiliary routine: fetches data from latbc buffer.
    !
    SUBROUTINE fetch_from_buffer_2D(latbc, name, target_buf, opt_p_ri)
      TYPE(t_latbc_data), INTENT(IN), TARGET :: latbc             !< contains read buffer
      CHARACTER(LEN=*),   INTENT(IN)         :: name              !< variable name
      REAL(wp),           INTENT(INOUT)      :: target_buf(:,:)   !< target buffer
      TYPE(t_reorder_info), POINTER, OPTIONAL, INTENT(IN) :: opt_p_ri         !< patch indices (cells, edges)
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::fetch_from_buffer_2D"
      INTEGER :: jm, j, jb, jl
      TYPE(t_reorder_info), POINTER :: p_ri !< patch indices (cells, edges)

      p_ri => latbc%patch_data%cells
      IF (PRESENT(opt_p_ri))  p_ri => opt_p_ri

      jm = get_field_index(latbc%buffer, TRIM(name))
      IF (jm <= 0)  CALL finish(routine//"_"//TRIM(name), "Internal error, invalid field index!")

!$OMP PARALLEL DO PRIVATE (j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
        jb = p_ri%own_blk(j) ! Block index in distributed patch
        jl = p_ri%own_idx(j) ! Line  index in distributed patch
        target_buf(jl,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,1,jb), wp)
      ENDDO
!$OMP END PARALLEL DO
    END SUBROUTINE fetch_from_buffer_2D


    SUBROUTINE fetch_from_buffer_3D_cells(latbc, name, target_buf)
      TYPE(t_latbc_data), INTENT(IN), TARGET :: latbc             !< contains read buffer
      CHARACTER(LEN=*),   INTENT(IN)         :: name              !< variable name
      REAL(wp),           INTENT(INOUT)      :: target_buf(:,:,:) !< target buffer

      CALL fetch_from_buffer_3D_generic(latbc, name, target_buf, latbc%patch_data%cells)
    END SUBROUTINE fetch_from_buffer_3D_cells


    SUBROUTINE fetch_from_buffer_3D_generic(latbc, name, target_buf, p_ri)
      TYPE(t_latbc_data),   INTENT(IN), TARGET :: latbc             !< contains read buffer
      CHARACTER(LEN=*),     INTENT(IN)         :: name              !< variable name
      REAL(wp),             INTENT(INOUT)      :: target_buf(:,:,:) !< target buffer
      TYPE(t_reorder_info), INTENT(IN)         :: p_ri              !< patch indices (cells, edges)
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::fetch_from_buffer_3D_generic"
      INTEGER :: jm, jk, j, jb, jl

      jm = get_field_index(latbc%buffer, TRIM(name))
      ! buffer%internal_name is constructed based on the inverted latbc_varnames_map_file dictionary. 
      ! A wrong name in the left column of latbc_varnames_map_file (internal name) will trigger the following error.
      IF (jm <= 0)  CALL finish(routine//"_"//TRIM(name), &
        &  "Internal error, invalid field index! Is "//TRIM(name)//" listed in latbc dict?")

      ! consistency check
      IF ((SIZE(target_buf,2) < latbc%buffer%nlev(jm)) .OR.  &
        & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
        CALL finish(routine//"_"//TRIM(name), "Internal error!")
      END IF

!$OMP PARALLEL DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jm)
        DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
          jb = p_ri%own_blk(j) ! Block index in distributed patch
          jl = p_ri%own_idx(j) ! Line  index in distributed patch
          target_buf(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    END SUBROUTINE fetch_from_buffer_3D_generic

    !-------------------------------------------------------------------------

  END MODULE mo_async_latbc_utils
