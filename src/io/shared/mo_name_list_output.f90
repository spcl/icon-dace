! Module handling synchronous and asynchronous output; supporting
! multiple I/O PEs and horizontal interpolation.
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
! 
!
! @todo In asynchronous I/O mode, windows are created but not freed
!
! @todo Several fields are allocated but not freed at the end of the
!       simulation. A pseudo-destructor should be implemented!
! @note: The spelling "name_list" (with underscore) is intended to make
!        clear that this does not pertain to a FORTRAN namelist but rather
!        to a list of names of output variables
!
! -------------------------------------------------------------------------
!
! The "namelist_output" module was originally written by Rainer
! Johanni. Some data structures used therein are duplicates of those
! created in the other parts of the model: In general, variable
! fields are introduced in ICON through the "add_var" mechanism in
! the module "shared/mo_var_list". This mechanism allocates "r_ptr"
! POINTERs for REAL(wp) variable fields, see the data structure in
! "t_var_list_element" (mo_var_list_element.f90). The "p_nh_state"
! variables, for example, then point to the same location. In the
! output, however, there exists a data structure "t_var_desc"
! (variable descriptor) which also contains an "r_ptr" POINTER. This
! also points to the original "r_ptr" location in memory.
!
! Exceptions and caveats for this described mechanism:
!
! - INTEGER fields are stored in "i_ptr" POINTERs.
! - After gathering the output data, so-called "post-ops" are
!   performed which modify the copied data (for example scaling from/to
!   percent values).
! - In asynchronous output mode, the "r_ptr" POINTERs are meaningless
!   on those PEs which are dedicated for output. These are NULL
!   pointers then.
!
! MPI roles in asynchronous communication:
!
! - Compute PEs create local memory windows, buffering all variables
!   for all output files (for the local horizontal grid partition).
!
! - Asynchronous I/O servers create trivial local memory windows of
!   size 1.
!
! - Additionally, when writing, the asynchronous I/O servers allocate
!   a 3D buffer for a single variable. This temporary field serves as
!   a target buffer for the one-sided MPI_GET operation.
!
! Transfer of meta-info:
!
!  Since parts of the variable's "info-field" TYPE(t_var_metadata) may change
!  during simulation, the following mechanism updates the metadata on the
!  asynchronous output PEs:
!  For each output file, a separate MPI window is created on work PE#0, where
!  the work root stores the current variable meta-info. This is then retrieved
!  via an additional MPI_GET by the I/O PEs.

MODULE mo_name_list_output

  ! constants
  USE mo_kind,                      ONLY: wp, i4, i8, dp, sp
  USE mo_impl_constants,            ONLY: max_dom, SUCCESS, MAX_TIME_LEVELS,       &
    &                                     ihs_ocean, BOUNDARY_MISSVAL, nlat_moc
  USE mo_cdi_constants,             ONLY: GRID_REGULAR_LONLAT, GRID_UNSTRUCTURED_VERT,              &
    &                                     GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_ZONAL
  USE mo_impl_constants_grf,        ONLY: grf_bdywidth_c
  USE mo_dynamics_config,           ONLY: iequations
  USE mo_cdi,                       ONLY: streamOpenWrite, FILETYPE_GRB2, streamDefTimestep, cdiEncodeTime, cdiEncodeDate, &
      &                                   CDI_UNDEFID, TSTEP_CONSTANT, FILETYPE_GRB, taxisDestroy, gridDestroy, &
      &                                   vlistDestroy, streamClose, streamWriteVarSlice, streamWriteVarSliceF, streamDefVlist, &
      &                                   streamSync, taxisDefVdate, taxisDefVtime, GRID_LONLAT, &
      &                                   streamDefCompType, CDI_COMPRESS_SZIP, &
      &                                   streamOpenAppend, streamInqVlist, vlistInqTaxis, vlistNtsteps, &
      &                                   vlistDuplicate, taxisDuplicate, &
      &                                   cdi_datatype_flt32, cdi_datatype_flt64
  USE mo_util_cdi,                  ONLY: cdiGetStringError
  ! utility functions
  USE mo_io_units,                  ONLY: FILENAME_MAX, find_next_free_unit
  USE mo_exception,                 ONLY: finish, message, message_text
  USE mo_util_string,               ONLY: t_keyword_list, associate_keyword, with_keywords,         &
  &                                       int2string
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_output, ltimer,      &
    &                                     timer_wait_for_async_io, print_timer
  USE mo_level_selection_types,     ONLY: t_level_selection
  USE mo_name_list_output_gridinfo, ONLY: write_grid_info_grb2, GRID_INFO_NONE
  USE mo_util_file,                 ONLY: util_rename, get_filename, get_path
  ! config
  USE mo_master_config,             ONLY: getModelBaseDir
  USE mo_grid_config,               ONLY: n_dom, l_limited_area
  USE mo_run_config,                ONLY: msg_level
  USE mo_io_config,                 ONLY: lkeep_in_sync,                   &
    &                                     config_lmask_boundary => lmask_boundary




  USE mo_gribout_config,            ONLY: gribout_config
  USE mo_parallel_config,           ONLY: p_test_run, use_dp_mpi2io, &
       num_io_procs, io_proc_chunk_size, nproma, pio_type
  USE mo_name_list_output_config,   ONLY: use_async_name_list_io
  ! data types
  USE mo_var_metadata_types,        ONLY: t_var_metadata, POST_OP_SCALE, POST_OP_LUC, &
    &                                     POST_OP_LIN2DBZ, var_metadata_get_size
  USE mo_reorder_info,              ONLY: t_reorder_info, ri_cpy_part2whole
  USE mo_name_list_output_types,    ONLY: t_output_file, icell, iedge, ivert, &
    &                                     msg_io_start, msg_io_done, &
    &                                     msg_io_meteogram_flush, &
    &                                     msg_io_shutdown, all_events, &
    &                                     t_var_desc, t_output_name_list
  USE mo_output_event_types,        ONLY: t_sim_step_info, t_par_output_event
  ! parallelization
  USE mo_communication,             ONLY: exchange_data, t_comm_gather_pattern,&
       idx_no, blk_no
  USE mo_mpi,                       ONLY: p_send, p_recv, p_barrier, stop_mpi,                      &
    &                                     p_mpi_wtime, p_irecv, p_wait, p_test, p_isend,            &
    &                                     p_comm_work, p_real_dp, p_real_sp, p_int,                 &
    &                                     my_process_is_stdio, my_process_is_mpi_test,              &
    &                                     my_process_is_mpi_workroot, my_process_is_work,           &
    &                                     my_process_is_io, my_process_is_mpi_ioroot,               &
    &                                     process_mpi_all_test_id, process_mpi_all_workroot_id,     &
    &                                     num_work_procs, p_pe, p_pe_work,                          &
    &                                     p_max, p_comm_work_2_io, mpi_request_null
  ! calendar operations
  USE mtime,                        ONLY: datetime, newDatetime, deallocateDatetime, OPERATOR(-),   &
    &                                     timedelta, max_datetime_str_len, &
    &                                     datetimeToString
  ! output scheduling
  USE mo_output_event_handler,      ONLY: is_output_step, check_open_file, check_close_file,        &
    &                                     pass_output_step, get_current_filename,                   &
    &                                     get_current_date,                                         &
    &                                     is_output_step_complete, is_output_event_finished,        &
    &                                     check_write_readyfile, blocking_wait_for_irecvs
  USE mo_name_list_output_stats,    ONLY: set_reference_time, interval_start, interval_end,         &
    &                                     interval_write_psfile
  ! output initialization
  USE mo_name_list_output_init,     ONLY: init_name_list_output, setup_output_vlist,                &
    &                                     varnames_dict, out_varnames_dict,                         &
    &                                     output_file, patch_info, lonlat_info,                     &
    &                                     collect_requested_ipz_levels, &
    &                                     create_vertical_axes, nlevs_of_var, zonal_ri, profile_ri
  USE mo_name_list_output_metadata, ONLY: metainfo_write_to_memwin, metainfo_get_from_buffer,       &
    &                                     metainfo_get_timelevel
  USE mo_level_selection,           ONLY: create_mipz_level_selections
  USE mo_grib2_util,                ONLY: set_GRIB2_timedep_keys, set_GRIB2_timedep_local_keys
  ! model domain
  USE mo_model_domain,              ONLY: t_patch, p_patch
  USE mo_loopindices,               ONLY: get_indices_c
  ! post-ops

  USE mo_post_op,                   ONLY: perform_post_op
  USE mo_meteogram_output,          ONLY: meteogram_init, meteogram_finalize, &
       meteogram_flush_file
  USE mo_meteogram_config,          ONLY: meteogram_output_config
  USE mo_intp_lonlat_types,         ONLY: lonlat_grids
  USE mo_fortran_tools, ONLY: insert_dimension, init, set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_name_list_output
  PUBLIC :: close_name_list_output
  PUBLIC :: istime4name_list_output
  PUBLIC :: istime4name_list_output_dom
  PUBLIC :: name_list_io_main_proc
  PUBLIC :: write_ready_files_cdipio

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output'

  !> constant for better readability
  INTEGER, PARAMETER :: WAIT_UNTIL_FINISHED = -1

  !> Internal switch for debugging output
  LOGICAL, PARAMETER :: ldebug  = .FALSE.

  INTEGER,          PARAMETER                 :: iUNKNOWN = 0
  INTEGER,          PARAMETER                 :: iINTEGER = 1
  INTEGER,          PARAMETER                 :: iREAL    = 2
  INTEGER,          PARAMETER                 :: iREAL_sp = 3


CONTAINS


  !------------------------------------------------------------------------------------------------
  !> open_output_file:
  !  Opens a output file and sets its vlist
  !
  !  Please note that this routine is only executed on one processor
  !  (for a specific file) and thus all calls to message get the
  !  all_print=.TRUE. argument so that the messages really appear in
  !  the log.
  !
  SUBROUTINE open_output_file(of, all_print)

    TYPE(t_output_file), INTENT(INOUT) :: of
    LOGICAL, INTENT(in)                :: all_print
    ! local variables:
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::open_output_file"
    CHARACTER(LEN=filename_max)    :: filename
    INTEGER                        :: name_len, part_idx
    LOGICAL                        :: lexist, lappend

    ! open/append file: as this is a preliminary solution only, I do not try to
    ! streamline the conditionals
    filename = get_current_filename(of%out_event)
    lappend  = .FALSE.

    ! check and reset filename, if data should be appended
    part_idx = INDEX(filename, '_part_')
    IF (part_idx > 0) THEN
      ! does the file to append to exist
      name_len = part_idx-1
      INQUIRE(file=filename(1:name_len), exist=lexist)
      IF (lexist) THEN
        ! open for append
        of%cdiFileID       = streamOpenAppend(filename(1:name_len))

        lappend            = .TRUE.
      ELSE
        ! file to append to does not exist that means we can use the name without part trailer
        of%cdiFileID       = streamOpenWrite(filename(1:name_len), of%output_type)
      ENDIF
    ELSE
      name_len = LEN_TRIM(filename)
      of%cdiFileID       = streamOpenWrite(filename(1:name_len), of%output_type)
      IF (gribout_config(of%phys_patch_id)%lgribout_compress_ccsds) THEN
        CALL streamDefCompType(of%cdiFileID, CDI_COMPRESS_SZIP)
      ENDIF
    ENDIF

    IF (of%cdiFileID < 0) THEN
      CALL cdiGetStringError(of%cdiFileID, message_text)
      CALL message(routine, message_text, all_print=.TRUE.)
      CALL finish (routine, 'open failed on '//filename(1:name_len))
    ELSE IF (msg_level >= 8) THEN
      IF (lappend) THEN
        CALL message (routine, 'to add more data, reopened '//filename(1:name_len),all_print=all_print)
      ELSE
        CALL message (routine, 'opened '//filename(1:name_len),all_print=all_print)
      END IF
    ENDIF

    IF (lappend) THEN
      ! get the already stored number of time steps
      of%cdiTimeIndex = vlistNtsteps(streamInqVlist(of%cdiFileID))
    ELSE
      ! assign the vlist (which must have ben set before)
        CALL streamDefVlist(of%cdiFileID, of%cdiVlistID)
      ! set cdi internal time index to 0 for writing time slices in netCDF
      of%cdiTimeIndex = 0
    ENDIF

  END SUBROUTINE open_output_file


  !------------------------------------------------------------------------------------------------
  !> Close all name_list files
  !
  SUBROUTINE close_name_list_output()
    ! local variables
    INTEGER :: i, ierror, prev_cdi_namespace

      !-- asynchronous I/O PEs (receiver):
      DO i = 1, SIZE(output_file)
        IF (output_file(i)%cdiFileID >= 0) THEN
          ! clean up level selection (if there is one):
          IF (ASSOCIATED(output_file(i)%level_selection)) THEN
            CALL output_file(i)%level_selection%finalize()
            DEALLOCATE(output_file(i)%level_selection)
            output_file(i)%level_selection => NULL()
          END IF
          CALL close_output_file(output_file(i))
          CALL destroy_output_vlist(output_file(i))
        END IF

        ! destroy vertical axes meta-data:
        CALL output_file(i)%verticalAxisList%finalize()
      ENDDO

    DEALLOCATE(output_file)

    ! destroy variable name dictionaries:
    CALL varnames_dict%finalize()
    CALL out_varnames_dict%finalize()

  END SUBROUTINE close_name_list_output


  !------------------------------------------------------------------------------------------------
  !> Close output stream and the associated file,
  !  destroy all vlist related data for this file
  !
  SUBROUTINE close_output_file(of)
    TYPE (t_output_file), INTENT(INOUT) :: of
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::close_output_file"

    IF(of%cdiFileID /= CDI_UNDEFID) CALL streamClose(of%cdiFileID)
    of%cdiFileID = CDI_UNDEFID

  END SUBROUTINE close_output_file


  !------------------------------------------------------------------------------------------------
  !> Close output stream and the associated file,
  !  destroy all vlist related data for this file
  !
  SUBROUTINE destroy_output_vlist(of)
    TYPE (t_output_file), INTENT(INOUT) :: of
    ! local variables
    INTEGER :: j

    IF(of%cdiVlistID /= CDI_UNDEFID) THEN
      IF(of%cdiCellGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiCellGridID)
      IF(of%cdiEdgeGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiEdgeGridID)
      IF(of%cdiVertGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiVertGridID)
      IF(of%cdiLonLatGridID /= CDI_UNDEFID) CALL gridDestroy(of%cdiLonLatGridID)
      CALL taxisDestroy(vlistInqTaxis(of%cdiVlistID))
      CALL vlistDestroy(of%cdiVlistID)
    ENDIF

    of%cdiCellGridID   = CDI_UNDEFID
    of%cdiEdgeGridID   = CDI_UNDEFID
    of%cdiVertGridID   = CDI_UNDEFID
    of%cdiLonLatGridID = CDI_UNDEFID
    of%cdiVlistID      = CDI_UNDEFID

  END SUBROUTINE destroy_output_vlist


  !------------------------------------------------------------------------------------------------
  !> Loop over all output_name_list's, write the ones for which output is due
  !  This routine also cares about opening the output files the first time
  !  and reopening the files after a certain number of steps.
  !
  SUBROUTINE write_name_list_output(jstep, opt_lhas_output, lacc)
    INTEGER,           INTENT(IN)   :: jstep             !< model step
    !> (Optional) Flag: .TRUE. if this async I/O PE has written during this step:
    LOGICAL, OPTIONAL, INTENT(OUT)  :: opt_lhas_output
    LOGICAL, OPTIONAL, INTENT(IN)   :: lacc
    ! local variables
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::write_name_list_output"
    INTEGER                           :: i, idate, itime, iret
    TYPE(datetime) :: io_datetime
    CHARACTER(LEN=filename_max+100)   :: text
    TYPE(t_par_output_event), POINTER :: ev
    INTEGER                           :: noutput_pe_list, io_proc_id
    INTEGER                           :: output_pe_list(MAX(1,num_io_procs))
    INTEGER :: prev_cdi_namespace
    INTEGER :: taxisID
    LOGICAL :: is_io, is_test
    LOGICAL :: lhas_output, all_print, do_sync
    LOGICAL :: ofile_is_active(SIZE(output_file)), &
         ofile_has_first_write(SIZE(output_file)), &
         ofile_is_assigned_here(SIZE(output_file))
    CHARACTER(len=MAX_DATETIME_STR_LEN) :: current_date_string
    LOGICAL :: lzacc

    IF (ltimer) CALL timer_start(timer_write_output)

    CALL set_acc_host_or_device(lzacc, lacc)

    is_io = my_process_is_io()
    is_test = my_process_is_mpi_test()
    all_print = .TRUE.


    lhas_output = .FALSE.

    ! during the following loop, we collect a list of all I/O PEs for
    ! which output is performed:
    output_pe_list(:) = -1
    noutput_pe_list   =  0

    OUTFILE_OPEN_CLOSE_LOOP : DO i=1,SIZE(output_file)
      ofile_is_active(i) = is_output_step(output_file(i)%out_event, jstep)
      IF (ofile_is_active(i)) THEN
        io_proc_id = output_file(i)%io_proc_id
        ofile_is_assigned_here(i) = &
          &      (is_io .AND. io_proc_id == p_pe_work) &
          & .OR. (.NOT. use_async_name_list_io .AND. .NOT. is_test &
          &       .AND. p_pe_work == 0)
        ofile_has_first_write(i) = check_open_file(output_file(i)%out_event)
        IF (ofile_is_assigned_here(i)) THEN
          ! -------------------------------------------------
          ! Check if files have to be closed
          ! -------------------------------------------------
          IF (check_close_file(output_file(i)%out_event)) THEN
            CALL close_output_file(output_file(i))
            IF (msg_level >= 8) THEN
              CALL message (routine, 'closed '//TRIM(get_current_filename(output_file(i)%out_event)),all_print=all_print)
            END IF
          END IF
          ! -------------------------------------------------
          ! Check if files have to be (re)opened
          ! -------------------------------------------------
          IF (ofile_has_first_write(i)) THEN
            IF (output_file(i)%cdiVlistID == CDI_UNDEFID)  &
                 &  CALL setup_output_vlist(output_file(i))
            CALL open_output_file(output_file(i), all_print)
          END IF
        END IF
      ELSE
        ofile_is_assigned_here(i) = .FALSE.
        ofile_has_first_write(i) = .FALSE.
      END IF
    END DO OUTFILE_OPEN_CLOSE_LOOP

    ! Go over all output files
    OUTFILE_WRITE_LOOP : DO i=1,SIZE(output_file)

      ! Skip this output file if it is not due for output!
      IF (.NOT. ofile_is_active(i)) CYCLE OUTFILE_WRITE_LOOP
      io_proc_id = output_file(i)%io_proc_id
      lhas_output = lhas_output .OR. ofile_is_assigned_here(i)

      IF (ofile_is_assigned_here(i)) THEN
        ! -------------------------------------------------
        ! Do the output
        ! -------------------------------------------------

        ! Notify user
        IF (msg_level >= 8) THEN
          CALL datetimeToString(get_current_date(output_file(i)%out_event), &
            &                   current_date_string)
            WRITE(text,'(5a,i0)')                                            &
              & 'Output to ',                                                &
              & TRIM(get_current_filename(output_file(i)%out_event)),        &
              & ' at simulation time ',                                      &
              & TRIM(current_date_string), &
              & ' by PE ', p_pe
          CALL message(routine, text,all_print=all_print)
        END IF

        ! convert time stamp string into
        ! year/month/day/hour/minute/second values using the mtime
        ! library:
        io_datetime = get_current_date(output_file(i)%out_event)
        idate = cdiEncodeDate(INT(io_datetime%date%year),   &
          &                   INT(io_datetime%date%month),  &
          &                   INT(io_datetime%date%day))
        itime = cdiEncodeTime(INT(io_datetime%time%hour),   &
          &                   INT(io_datetime%time%minute), &
          &                   INT(io_datetime%time%second))
        taxisID = vlistInqTaxis(streamInqVlist(output_file(i)%cdiFileID))
        CALL taxisDefVdate(taxisID, idate)
        CALL taxisDefVtime(taxisID, itime)
        iret = streamDefTimestep(output_file(i)%cdiFileId, output_file(i)%cdiTimeIndex)
        output_file(i)%cdiTimeIndex = output_file(i)%cdiTimeIndex + 1
      END IF

      IF(is_io) THEN
      ELSE
        CALL write_name_list(output_file(i), ofile_has_first_write(i), i, lzacc)
        do_sync = lkeep_in_sync .AND. ofile_is_assigned_here(i)
      ENDIF

      ! -------------------------------------------------
      ! add I/O PE of output file to the "output_list"
      ! -------------------------------------------------
      IF (ALL(output_pe_list(1:noutput_pe_list) /= io_proc_id)) THEN
        noutput_pe_list = noutput_pe_list + 1
        output_pe_list(noutput_pe_list) = io_proc_id
      END IF

      ! GRB2 format: define geographical longitude, latitude as special
      ! variables "RLON", "RLAT"
      IF (ofile_is_assigned_here(i) &
        .AND. patch_info(output_file(i)%phys_patch_id)%grid_info_mode         &
        &     /= GRID_INFO_NONE                                               &
        .AND. output_file(i)%name_list%output_grid                            &
        .AND. output_file(i)%name_list%filetype == FILETYPE_GRB2              &
        .AND. check_close_file(output_file(i)%out_event,                      &
        &           output_file(i)%out_event%output_event%i_event_step+1)) THEN
        CALL write_grid_info_grb2(output_file(i), patch_info)
      END IF

      ! -------------------------------------------------
      ! hand-shake protocol: step finished!
      ! -------------------------------------------------
      CALL pass_output_step(output_file(i)%out_event)
    ENDDO OUTFILE_WRITE_LOOP

    ! If asynchronous I/O is enabled, the compute PEs can now start
    ! the I/O PEs
    ! Handle incoming "output step completed" messages: After all
    ! participating I/O PE's have acknowledged the completion of their
    ! write processes, we trigger a "ready file" on the first I/O PE.
    IF (.NOT. is_test) THEN
       IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
            & (.NOT. use_async_name_list_io .AND. my_process_is_mpi_workroot())) THEN
          ev => all_events
          HANDLE_COMPLETE_STEPS : DO WHILE (ASSOCIATED(ev))
            IF (is_output_step_complete(ev) .AND.  &
              & .NOT. is_output_event_finished(ev)) THEN
              !--- write ready file
              !
              ! FIXME: for CDI-PIO this needs to be moved to the I/O
              ! processes which are aware what has been written when
              IF (check_write_readyfile(ev%output_event)) CALL write_ready_file(ev)
              ! launch a non-blocking request to all participating PEs to
              ! acknowledge the completion of the next output event
                ev%output_event%i_event_step = ev%output_event%i_event_step + 1
              ev => ev%next
            ELSE
              ev => ev%next
            END IF
          END DO HANDLE_COMPLETE_STEPS
       END IF
    END IF
    IF (PRESENT(opt_lhas_output)) opt_lhas_output = lhas_output
    IF (ltimer) CALL timer_stop(timer_write_output)
    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": write_name_list_output done."
  END SUBROUTINE write_name_list_output

  SUBROUTINE write_ready_files_cdipio
    TYPE(t_par_output_event), POINTER :: ev
    !fixme: this needs a mechanism to enforce disk flushes via streamsync
    IF (p_pe_work == 0) THEN
      ev => all_events
      DO WHILE (ASSOCIATED(ev))
        IF (.NOT. is_output_event_finished(ev)) THEN
          !--- write ready file
          !
          IF (check_write_readyfile(ev%output_event)) CALL write_ready_file(ev)
          ! launch a non-blocking request to all participating PEs to
          ! acknowledge the completion of the next output event
          ev%output_event%i_event_step = ev%output_event%i_event_step + 1
        END IF
        ev => ev%next
      END DO
    END IF
  END SUBROUTINE write_ready_files_cdipio


  !------------------------------------------------------------------------------------------------
  !> Create a "ready file"
  !
  !  A "ready file" is a technique for handling dependencies between
  !  the NWP processes at DWD: When a program - parallel or
  !  sequential, shell script or binary - produces some output which
  !  is necessary for other running applications, then the completion
  !  of the write process signals this by creating a small file (size:
  !  a few bytes). Only when this file exists, the second program
  !  starts reading its input data. Implicity, this assumes that a
  !  file system creates (and closes) files in the same order as they
  !  are written by the program.
  !
  SUBROUTINE write_ready_file(ev)
    TYPE(t_par_output_event), INTENT(IN) :: ev
    ! local variables
    CHARACTER(LEN=*), PARAMETER         :: routine = modname//"::write_ready_file"
    CHARACTER(LEN=*), PARAMETER         :: tmp_prefix = ".."
    CHARACTER(LEN=FILENAME_MAX)         :: rdy_filename, tmp_filename
    CHARACTER(LEN=8)                    :: forecast_delta_str
    TYPE(datetime)                      :: mtime_begin, mtime_date, current_date
    TYPE(timedelta)                     :: forecast_delta
    INTEGER                             :: iunit, tlen, iret
    TYPE (t_keyword_list), POINTER      :: keywords
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string, current_date_string

    current_date = get_current_date(ev)
    ! compute current forecast time (delta):
    mtime_date     = current_date
    mtime_begin    = ev%output_event%event_data%sim_start
    forecast_delta = mtime_date - mtime_begin

    WRITE (forecast_delta_str,'(4(i2.2))') forecast_delta%day, forecast_delta%hour, &
      &                                    forecast_delta%minute, forecast_delta%second
    WRITE (dtime_string,'(i4.4,2(i2.2),a,3(i2.2),a)')                                                 &
      &                      mtime_date%date%year, mtime_date%date%month, mtime_date%date%day, 'T',   &
      &                      mtime_date%time%hour, mtime_date%time%minute, mtime_date%time%second, 'Z'

    NULLIFY(keywords)
    ! substitute tokens in ready file name
    CALL associate_keyword("<path>",      TRIM(getModelBaseDir()),    keywords)
    CALL datetimeToString(current_date, current_date_string)
    CALL associate_keyword("<datetime>",  TRIM(current_date_string),  keywords)
    CALL associate_keyword("<ddhhmmss>",  forecast_delta_str,         keywords)
    CALL associate_keyword("<datetime2>", TRIM(dtime_string),         keywords)
    rdy_filename = with_keywords(keywords, ev%output_event%event_data%name)
    tlen = LEN_TRIM(rdy_filename)
    IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
      & (.NOT. use_async_name_list_io .AND. my_process_is_stdio())) THEN
      WRITE (0,*) 'Write ready file "', rdy_filename(1:tlen), '"'
    END IF

    ! Actually create ready file.
    !
    ! This procedure is carried out in two steps: First, a file with
    ! the name "tmp_prefix+rdy_filename" is created. After the file
    ! has been closed, it is then renamed into "rdy_filename" in a
    ! second step.
    ! This detour is necessary when another process polls the output
    ! directory and relies on a "complete" ready file.
    tmp_filename = TRIM(get_path(rdy_filename(1:tlen)))//tmp_prefix//TRIM(get_filename(rdy_filename(1:tlen)))
    iunit = find_next_free_unit(10,20)
    OPEN (iunit, file=TRIM(tmp_filename), form='formatted')
    WRITE(iunit, '(A)') 'ready'
    CLOSE(iunit)

    iret = util_rename(TRIM(tmp_filename), rdy_filename(1:tlen))
  END SUBROUTINE write_ready_file


  !------------------------------------------------------------------------------------------------
  !> Write an output name list. Called by non-IO PEs.
  !
  SUBROUTINE write_name_list(of, is_first_write, file_idx, lacc)


    TYPE (t_output_file), INTENT(INOUT), TARGET :: of
    LOGICAL,              INTENT(IN)            :: is_first_write
    INTEGER,              INTENT(IN)            :: file_idx ! File index in output_file(:) array
    LOGICAL, OPTIONAL,    INTENT(IN)            :: lacc
    ! local variables:
    CHARACTER(LEN=*), PARAMETER                 :: routine = modname//"::write_name_list"
    INTEGER                                     :: tl, i_dom, i_log_dom, iv, jk, &
      &                                            nlevs, lonlat_id,          &
      &                                            idata_type
    INTEGER                                     :: ioff
    TYPE (t_var_metadata), POINTER              :: info
    TYPE (t_reorder_info), POINTER              :: p_ri

    REAL(wp), ALLOCATABLE, TARGET :: r_ptr_m(:,:,:)
    REAL(sp), ALLOCATABLE, TARGET :: s_ptr_m(:,:,:)
    INTEGER, ALLOCATABLE, TARGET :: i_ptr_m(:,:,:)

    REAL(wp), POINTER :: r_ptr(:,:,:)
    REAL(sp), POINTER :: s_ptr(:,:,:)
    INTEGER, POINTER :: i_ptr(:,:,:)
    TYPE(t_comm_gather_pattern), POINTER        :: p_pat
    LOGICAL                                     :: var_ignore_level_selection
    INTEGER                                     :: last_bdry_index
    INTEGER :: info_nlevs

    INTEGER :: ipost_op_type, alloc_shape(3), alloc_shape_op(3)
    LOGICAL :: post_op_apply

    LOGICAL :: is_test, is_stdio
    LOGICAL, PARAMETER :: participate_in_async_io = .FALSE.
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)
    ! Offset in memory window for async I/O
    ioff = 0

    i_dom = of%phys_patch_id
    i_log_dom = of%log_patch_id

    tl = 0 ! to prevent warning

    is_test = my_process_is_mpi_test()
    is_stdio = my_process_is_stdio()

    ! "lmask_boundary": Some of the output fields are not updated with
    ! meaningful values in the vicinity of the lateral domain
    ! boundary. To avoid spurious data on these triangle cells (which
    ! could also spoil the GRIB range compression), the user may
    ! choose to set them to a "missing value". Implementation details:
    ! In the "synchronous" output mode, the implementation exploits
    ! the fact that all (global) indices for the lateral boundary
    ! region are ordered to the start of the data array. Therefore,
    ! only the computation of a limit index "last_bdry_index" is
    ! required to mask out the lateral points. In the asynchronous
    ! output mode, on the other hand, the compute processes possess
    ! only a portion of the output field and therefore need to loop
    ! over the lateral triangles block- and line-wise. This feature
    ! can be (de-)activated for specific variables through the
    ! "info%lmask_boundary" metadata flag. It also depends on a global
    ! namelist switch "io_nml/lmask_boundary" (LOGICAL, default:
    ! false).
    !
    ! Only for synchronous output mode: communicate the largest global
    ! index of the lateral boundary cells, if required:

    IF ((.NOT. participate_in_async_io) .AND. config_lmask_boundary(i_log_dom))  THEN
      last_bdry_index = get_last_bdry_index(i_log_dom)
    ELSE
      last_bdry_index = 0
    END IF

    ! ----------------------------------------------------
    ! Go over all name list variables for this output file
    ! ----------------------------------------------------
    DO iv = 1, of%num_vars

      info => of%var_desc(iv)%info

      ! inspect time-constant variables only if we are writing the
      ! first step in this file:
      IF ((info%isteptype == TSTEP_CONSTANT) .AND. .NOT. is_first_write) CYCLE

      ! Check if first dimension of array is nproma.
      ! Otherwise we got an array which is not suitable for this output scheme.
    ! IF(info%used_dimensions(1) /= nproma) &
    !   CALL finish(routine,'1st dim is not nproma: '//TRIM(info%name))

      idata_type = iUNKNOWN

      ! For time level dependent elements: set time level and check if
      ! time level is present:
        ! set a default time level (which is not used anyway, but must
        ! be a valid array subscript):
      tl = 1

      IF (.NOT. ASSOCIATED(of%var_desc(iv)%r_ptr)  .AND. &
        & .NOT. ASSOCIATED(of%var_desc(iv)%s_ptr)  .AND. &
        & .NOT. ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
        tl = metainfo_get_timelevel(info,i_log_dom)
        IF(tl<=0 .OR. tl>max_time_levels) &
          CALL finish(routine, 'Illegal time level in nnow()/nnow_rcf()')
        ! Check if present
        IF (.NOT. ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)   .AND.   &
          & .NOT. ASSOCIATED(of%var_desc(iv)%tlev_sptr(tl)%p)   .AND.   &
          & .NOT. ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
          CALL finish(routine,'Actual timelevel not in '//TRIM(info%name))
        END IF
      ENDIF


      ! determine, if this is a REAL or an INTEGER variable:
      IF (ASSOCIATED(of%var_desc(iv)%r_ptr) .OR.  &
        & ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)) THEN
        idata_type = iREAL
      ELSE IF (ASSOCIATED(of%var_desc(iv)%s_ptr) .OR.  &
        & ASSOCIATED(of%var_desc(iv)%tlev_sptr(tl)%p)) THEN
        idata_type = iREAL_sp
      ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr) .OR.  &
        & ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
        idata_type = iINTEGER
      END IF

      CALL get_ptr_to_var_data(i_ptr, r_ptr, s_ptr, &
        &                      tl, of%var_desc(iv), info)

      ! --------------------------------------------------------
      ! Perform post-ops (small arithmetic operations on fields)
      ! --------------------------------------------------------

      ipost_op_type = info%post_op%ipost_op_type
      post_op_apply &
        = ipost_op_type == post_op_scale .OR. ipost_op_type == post_op_luc .OR. ipost_op_type == post_op_lin2dbz
      IF ( post_op_apply ) THEN
        IF (idata_type == iREAL) THEN
          alloc_shape = SHAPE(r_ptr)
          IF (ALLOCATED(r_ptr_m)) THEN
            alloc_shape_op = SHAPE(r_ptr_m)
            IF (ANY(alloc_shape_op /= alloc_shape)) THEN
              DEALLOCATE(r_ptr_m)
              ALLOCATE(r_ptr_m(alloc_shape(1), alloc_shape(2), alloc_shape(3)))
            END IF
          ELSE
            ALLOCATE(r_ptr_m(alloc_shape(1), alloc_shape(2), alloc_shape(3)))
          END IF
          r_ptr_m = r_ptr
          r_ptr => r_ptr_m
          CALL perform_post_op(info%post_op, r_ptr, lacc=lzacc)
        ELSE IF (idata_type == iREAL_sp) THEN
          alloc_shape = SHAPE(s_ptr)
          IF (ALLOCATED(s_ptr_m)) THEN
            alloc_shape_op = SHAPE(s_ptr_m)
            IF (ANY(alloc_shape_op /= alloc_shape)) THEN
              DEALLOCATE(s_ptr_m)
              ALLOCATE(s_ptr_m(alloc_shape(1), alloc_shape(2), alloc_shape(3)))
            END IF
          ELSE
            ALLOCATE(s_ptr_m(alloc_shape(1), alloc_shape(2), alloc_shape(3)))
          END IF
          s_ptr_m = s_ptr
          s_ptr => s_ptr_m
          CALL perform_post_op(info%post_op, s_ptr, lacc=lzacc)
        ELSE IF (idata_type == iINTEGER) THEN
          alloc_shape = SHAPE(i_ptr)
          IF (ALLOCATED(i_ptr_m)) THEN
            alloc_shape_op = SHAPE(i_ptr_m)
            IF (ANY(alloc_shape_op /= alloc_shape)) THEN
              DEALLOCATE(i_ptr_m)
              ALLOCATE(i_ptr_m(alloc_shape(1), alloc_shape(2), alloc_shape(3)))
            END IF
          ELSE
            ALLOCATE(i_ptr_m(alloc_shape(1), alloc_shape(2), alloc_shape(3)))
          END IF
          i_ptr_m = i_ptr
          i_ptr => i_ptr_m
          CALL perform_post_op(info%post_op, i_ptr, lacc=lzacc)
        ENDIF
      END IF

      nlevs = nlevs_of_var(info, of%level_selection, var_ignore_level_selection)
      IF (var_ignore_level_selection .AND. is_stdio .AND. msg_level >= 15) &
          &   WRITE (0,'(2a)') &
          &         "warning: ignoring level selection for variable ", &
          &         TRIM(info%name)

      ! Get pointer to appropriate reorder_info
      nullify(p_ri, p_pat)
      SELECT CASE (info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        p_ri     => patch_info(i_dom)%ri(icell)
        p_pat    => patch_info(i_dom)%p_pat_c
      CASE (GRID_LONLAT)
        p_ri => profile_ri
        NULLIFY(p_pat)
      CASE (GRID_ZONAL)
        p_ri => zonal_ri
        NULLIFY(p_pat)
      CASE (GRID_UNSTRUCTURED_EDGE)
        p_ri     => patch_info(i_dom)%ri(iedge)
        p_pat    => patch_info(i_dom)%p_pat_e
      CASE (GRID_UNSTRUCTURED_VERT)
        p_ri     => patch_info(i_dom)%ri(ivert)
        p_pat    => patch_info(i_dom)%p_pat_v
      CASE (GRID_REGULAR_LONLAT)
        lonlat_id = info%hor_interp%lonlat_id
        p_ri     => lonlat_info(lonlat_id, i_log_dom)%ri
        p_pat    => lonlat_grids%list(lonlat_id)%p_pat(i_log_dom)
      CASE default
        CALL finish(routine,'unknown grid type')
      END SELECT

        IF (.NOT.use_async_name_list_io .OR. is_test) THEN
          CALL gather_on_workroot_and_write(of, idata_type, r_ptr, s_ptr, &
            i_ptr, p_ri%n_glb, iv, last_bdry_index, &
            nlevs, var_ignore_level_selection, p_pat, info)
        END IF

    ENDDO


  END SUBROUTINE write_name_list

  FUNCTION get_last_bdry_index(i_log_dom) RESULT(last_bdry_index)
    INTEGER, INTENT(in) :: i_log_dom
    INTEGER :: last_bdry_index

    TYPE(t_patch), POINTER                      :: ptr_patch
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: jl_start, jl_end, jl
    INTEGER :: max_glb_idx, tmp_dummy
    INTEGER, POINTER, CONTIGUOUS :: glb_index(:)
    CALL finish(modname//":get_last_bdry_index", "Caution: bug ahead! Synchronous (no IO proc) is deprecated so this bug won't be fixed.")
    rl_start   = 1
    rl_end     = grf_bdywidth_c
    CALL get_bdry_blk_idx(i_log_dom, &
      &                   i_startidx, i_endidx, i_startblk, i_endblk)
    max_glb_idx = 0
    ptr_patch => p_patch(i_log_dom)
    glb_index => ptr_patch%cells%decomp_info%glb_index
    CALL get_indices_c(ptr_patch, i_startblk, i_startblk, i_endblk, &
         i_startidx, tmp_dummy, rl_start, rl_end)
    CALL get_indices_c(ptr_patch, i_endblk, i_startblk, i_endblk, &
         tmp_dummy, i_endidx, rl_start, rl_end)
    jl_start = (i_startblk-1)*nproma + i_startidx
    jl_end = (i_endblk-1)*nproma + i_endidx
    DO jl=jl_start,jl_end
      max_glb_idx = MAX(max_glb_idx, glb_index(jl))
    END DO
    last_bdry_index = p_max(max_glb_idx, p_comm_work)
  END FUNCTION get_last_bdry_index

  SUBROUTINE get_ptr_to_var_data(i_ptr, r_ptr, s_ptr, tl, var_desc, info)
    TYPE (t_var_metadata), INTENT(in) :: info
    REAL(wp), POINTER, INTENT(out) :: r_ptr(:,:,:)
    REAL(sp), POINTER, INTENT(out) :: s_ptr(:,:,:)
    INTEGER, POINTER, INTENT(out) :: i_ptr(:,:,:)
    INTEGER, INTENT(in) :: tl
    TYPE(t_var_desc), TARGET, INTENT(in) :: var_desc

    REAL(wp), SAVE, TARGET :: r_dummy(1,1,1)
    REAL(wp), POINTER :: r_ptr_t(:,:,:,:,:,:), r_ptr_5d(:,:,:,:,:)
    REAL(sp), SAVE, TARGET :: s_dummy(1,1,1)
    REAL(sp), POINTER :: s_ptr_t(:,:,:,:,:,:), s_ptr_5d(:,:,:,:,:)
    INTEGER, SAVE, TARGET :: i_dummy(1,1,1)
    INTEGER, POINTER :: i_ptr_t(:,:,:,:,:,:), i_ptr_5d(:,:,:,:,:)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_ptr_to_var_data"
    INTEGER :: var_ref_pos, nindex

    r_ptr => r_dummy
    s_ptr => s_dummy
    i_ptr => i_dummy

    NULLIFY(r_ptr_5d, s_ptr_5d, i_ptr_5d)
    IF      (     ASSOCIATED(var_desc%r_ptr) &
      &      .OR. ASSOCIATED(var_desc%tlev_rptr(tl)%p)) THEN
      IF (ASSOCIATED(var_desc%r_ptr)) THEN
        r_ptr_5d => var_desc%r_ptr
      ELSE
        r_ptr_5d => var_desc%tlev_rptr(tl)%p
      END IF
    ELSE IF (     ASSOCIATED(var_desc%s_ptr) &
      &      .OR. ASSOCIATED(var_desc%tlev_sptr(tl)%p)) THEN
      IF (ASSOCIATED(var_desc%s_ptr)) THEN
        s_ptr_5d => var_desc%s_ptr
      ELSE
        s_ptr_5d => var_desc%tlev_sptr(tl)%p
      END IF
    ELSE IF (     ASSOCIATED(var_desc%i_ptr) &
      &      .OR. ASSOCIATED(var_desc%tlev_iptr(tl)%p)) THEN
      IF (ASSOCIATED(var_desc%i_ptr)) THEN
        i_ptr_5d => var_desc%i_ptr
      ELSE
        i_ptr_5d => var_desc%tlev_iptr(tl)%p
      END IF
    ELSE
      CALL finish(routine, "Internal error!")
    END IF

    nindex = MERGE(info%ncontained, 1, info%lcontained)

    SELECT CASE (info%ndims)
    CASE (1)
      IF (info%lcontained .AND. (info%var_ref_pos /= -1))  &
           & CALL finish(routine, "Internal error (ndims=1, lcontained)")
      IF (ASSOCIATED(r_ptr_5d)) THEN
        IF (info%hgrid == grid_lonlat) THEN
          CALL insert_dimension(r_ptr_t, r_ptr_5d, 1)
          r_ptr => r_ptr_t(:,:,1:1,1,1,1)
        ELSE
          r_ptr => r_ptr_5d(:,1:1,1:1,1,1)
        END IF
      ELSE IF (ASSOCIATED(s_ptr_5d)) THEN
        s_ptr => s_ptr_5d(:,1:1,1:1,1,1)
      ELSE IF (ASSOCIATED(i_ptr_5d)) THEN
        i_ptr => i_ptr_5d(:,1:1,1:1,1,1)
      ELSE
        CALL finish(routine, "Internal error (not found vardata pointer)")
      ENDIF

    CASE (2)
      var_ref_pos = 3
      IF (info%lcontained)  var_ref_pos = info%var_ref_pos
      IF (var_ref_pos < 1 .OR. var_ref_pos > 3) THEN
        WRITE (message_text, '(2a,i0)') TRIM(info%name), &
             ": internal error! var_ref_pos=", var_ref_pos
        GO TO 999
      END IF
      IF      (ASSOCIATED(r_ptr_5d)) THEN
        SELECT CASE(var_ref_pos)
        CASE (1)
          CALL insert_dimension(r_ptr_t, r_ptr_5d, 3)
          r_ptr => r_ptr_t(nindex,:,1:1,:,1,1)
        CASE (2)
          r_ptr => r_ptr_5d(:,nindex:nindex,:,1,1)
        CASE (3)
          IF (info%hgrid == grid_zonal) THEN
            CALL insert_dimension(r_ptr, r_ptr_5d(:,:,nindex,1,1), 1)
          ELSE
            CALL insert_dimension(r_ptr_t, r_ptr_5d, 2)
            r_ptr => r_ptr_t(:,1:1,:,nindex,1,1)
          END IF
        END SELECT
      ELSE IF (ASSOCIATED(s_ptr_5d)) THEN
        SELECT CASE(var_ref_pos)
        CASE (1)
          CALL insert_dimension(s_ptr_t, s_ptr_5d, 3)
          s_ptr => s_ptr_t(nindex,:,1:1,:,1,1)
        CASE (2)
          s_ptr => s_ptr_5d(:,nindex:nindex,:,1,1)
        CASE (3)
          CALL insert_dimension(s_ptr_t, s_ptr_5d, 2)
          IF (info%hgrid == GRID_ZONAL) THEN
            CALL insert_dimension(s_ptr, s_ptr_t(:,1,:,nindex,1,1), 2)
          ELSE
            s_ptr => s_ptr_t(:,1:1,:,nindex,1,1)
          END IF
        END SELECT
      ELSE IF (ASSOCIATED(i_ptr_5d)) THEN
        SELECT CASE(var_ref_pos)
        CASE (1)
          CALL insert_dimension(i_ptr_t, i_ptr_5d, 3)
          i_ptr => i_ptr_t(nindex,:,1:1,:,1,1)
        CASE (2)
          i_ptr => i_ptr_5d(:,nindex:nindex,:,1,1)
        CASE (3)
          CALL insert_dimension(i_ptr_t, i_ptr_5d, 2)
          IF (info%hgrid == GRID_ZONAL) THEN
            CALL insert_dimension(i_ptr, i_ptr_t(:,1,:,nindex,1,1), 2)
          ELSE
            i_ptr => i_ptr_t(:,1:1,:,nindex,1,1)
          END IF
        END SELECT
      ENDIF
    CASE (3)

      var_ref_pos = 4
      IF (info%lcontained)  var_ref_pos = info%var_ref_pos
      IF (var_ref_pos < 1 .OR. var_ref_pos > 4) THEN
        WRITE (message_text, '(2a,i0)') TRIM(info%name), &
             ": internal error! var_ref_pos=", var_ref_pos
        GO TO 999
      END IF

      ! 3D fields: Here we could just set a pointer to the
      ! array... if there were no post-ops
      IF      (ASSOCIATED(r_ptr_5d)) THEN
        SELECT CASE(var_ref_pos)
        CASE (1)
          r_ptr => r_ptr_5d(nindex,:,:,:,1)
        CASE (2)
          r_ptr => r_ptr_5d(:,nindex,:,:,1)
        CASE (3)
          r_ptr => r_ptr_5d(:,:,nindex,:,1)
        CASE (4)
          r_ptr => r_ptr_5d(:,:,:,nindex,1)
        END SELECT
      ELSE IF (ASSOCIATED(s_ptr_5d)) THEN
        SELECT CASE(var_ref_pos)
        CASE (1)
          s_ptr => s_ptr_5d(nindex,:,:,:,1)
        CASE (2)
          s_ptr => s_ptr_5d(:,nindex,:,:,1)
        CASE (3)
          s_ptr => s_ptr_5d(:,:,nindex,:,1)
        CASE (4)
          s_ptr => s_ptr_5d(:,:,:,nindex,1)
        END SELECT
      ELSE IF (ASSOCIATED(i_ptr_5d)) THEN
        SELECT CASE(var_ref_pos)
        CASE (1)
          i_ptr => i_ptr_5d(nindex,:,:,:,1)
        CASE (2)
          i_ptr => i_ptr_5d(:,nindex,:,:,1)
        CASE (3)
          i_ptr => i_ptr_5d(:,:,nindex,:,1)
        CASE (4)
          i_ptr => i_ptr_5d(:,:,:,nindex,1)
        END SELECT
      END IF
    CASE DEFAULT
      WRITE (message_text, '(2a,i0)') TRIM(info%name), &
           ": internal error! unhandled info%ndims=", info%ndims
      GO TO 999
    END SELECT

    

    IF      (ASSOCIATED(r_ptr_5d)) THEN
      !$ACC UPDATE HOST(r_ptr) ASYNC(1) IF(i_am_accel_node .AND. acc_is_present(r_ptr))
    ELSE IF (ASSOCIATED(s_ptr_5d)) THEN
      !$ACC UPDATE HOST(s_ptr) ASYNC(1) IF(i_am_accel_node .AND. acc_is_present(s_ptr))
    ELSE IF (ASSOCIATED(i_ptr_5d)) THEN
      !$ACC UPDATE HOST(i_ptr) ASYNC(1) IF(i_am_accel_node .AND. acc_is_present(i_ptr))
    ENDIF
    !$ACC WAIT(1) IF(i_am_accel_node)

    RETURN
999 CALL finish(routine,message_text)

  END SUBROUTINE get_ptr_to_var_data

  SUBROUTINE gather_on_workroot_and_write(of, idata_type, r_ptr, s_ptr, &
       i_ptr, n_glb, iv, last_bdry_index, &
       nlevs, var_ignore_level_selection, pat, info)
    TYPE (t_output_file), INTENT(IN) :: of
    INTEGER, INTENT(in) :: idata_type, iv, nlevs, last_bdry_index
    LOGICAL, INTENT(in) :: var_ignore_level_selection
    REAL(dp), INTENT(in) :: r_ptr(:,:,:)
    REAL(sp), INTENT(in) :: s_ptr(:,:,:)
    INTEGER, INTENT(in)  :: i_ptr(:,:,:), n_glb

    REAL(dp), ALLOCATABLE :: r_out_dp(:)
    INTEGER, ALLOCATABLE :: r_out_int(:)
    REAL(sp), ALLOCATABLE :: r_out_sp(:)
    REAL(wp), ALLOCATABLE :: r_out_recv(:)
    REAL(wp), PARAMETER :: SYNC_ERROR_PRINT_TOL = 1e-13_wp
    REAL(wp) :: missval
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::gather_on_workroot_and_write"

    TYPE(t_comm_gather_pattern), INTENT(in), POINTER :: pat
    TYPE (t_var_metadata), INTENT(in) :: info


    INTEGER :: lev, lev_idx, i, nmiss
    LOGICAL :: l_error, have_grib, lwrite_single_precision
    LOGICAL :: is_test, is_mpi_workroot
    LOGICAL :: make_level_selection

    is_mpi_workroot = my_process_is_mpi_workroot()

    is_test = my_process_is_mpi_test()

    have_GRIB =      of%output_type == FILETYPE_GRB  &
      &         .OR. of%output_type == FILETYPE_GRB2
    lwrite_single_precision =   (.NOT. use_dp_mpi2io) .AND. (.NOT. have_GRIB) &
      &                       .OR. idata_type == iREAL_sp

    IF (idata_type == iREAL) THEN
      ALLOCATE(r_out_dp(MERGE(n_glb, 0, is_mpi_workroot)))
    END IF
    IF ((idata_type == iREAL_sp) .OR. lwrite_single_precision) THEN
      ALLOCATE(r_out_sp(MERGE(n_glb, 0, is_mpi_workroot)))
    END IF
    IF (idata_type == iINTEGER) THEN
      IF ( .NOT. ALLOCATED(r_out_sp) ) ALLOCATE(r_out_sp(MERGE(n_glb, 0, is_mpi_workroot)))
      IF ( .NOT. ALLOCATED(r_out_dp) ) ALLOCATE(r_out_dp(MERGE(n_glb, 0, is_mpi_workroot)))
      ALLOCATE(r_out_int(MERGE(n_glb, 0, is_mpi_workroot)))
    END IF

    IF(is_mpi_workroot) THEN

      IF (is_test) THEN

        IF (p_test_run .AND. use_dp_mpi2io) ALLOCATE(r_out_recv(n_glb))

      ELSE

        CALL set_time_varying_metadata(of, info, of%var_desc(iv)%info_ptr)
      END IF ! is_test
    END IF ! is_mpi_workroot


    ! set missval flag, if applicable
    !
    ! Layerwise missing value masks are available in GRIB output format
    ! only. A missing value might be set by the user (info%lmiss) or
    ! automatically on nest boundary regions.
    nmiss = MERGE(1, 0, (info%lmiss &
      &                  .OR. (info%lmask_boundary &
      &                        .AND. ANY(config_lmask_boundary(:))) ) &
      &                 .AND. last_bdry_index > 0)
    IF (.NOT. have_GRIB .AND. nmiss /= 0) THEN
      ! this is the wrong place to set nmiss for NetCDF
      CALL finish(routine, "Caution! Bug ahead. Synchronous (no IO proc) is deprecated so this bug won't be fixed.")
    END IF

    make_level_selection = ASSOCIATED(of%level_selection) &
      &              .AND. (.NOT. var_ignore_level_selection) &
      &              .AND. (info%ndims > 2)
    ! For all levels (this needs to be done level-wise in order to reduce
    !                 memory consumption)
    DO lev = 1, nlevs
      ! -------------------
      ! No asynchronous I/O
      ! -------------------
      !
      ! gather the array on stdio PE and write it out there
      IF ( info%hgrid == GRID_LONLAT ) THEN
        IF (is_mpi_workroot) THEN
          IF      (idata_type == iREAL ) THEN
            r_out_dp(:)  = r_ptr(:,1,1)
          ELSE IF (idata_type == iREAL_sp ) THEN
            r_out_sp(:)  = s_ptr(:,1,1)
          ELSE IF (idata_type == iINTEGER) THEN
            r_out_int(:) = i_ptr(:,1,1)
          END IF
        END IF
      ELSE IF ( info%hgrid == GRID_ZONAL ) THEN ! 1deg zonal grid
        lev_idx = lev
        IF (is_mpi_workroot) THEN
          IF      (idata_type == iREAL ) THEN
            r_out_dp(:)  = r_ptr(1,lev_idx,:)
          ELSE IF (idata_type == iREAL_sp ) THEN
            r_out_sp(:)  = s_ptr(1,lev_idx,:)
          ELSE IF (idata_type == iINTEGER) THEN
            r_out_int(:) = i_ptr(1,lev_idx,:)
          END IF
        END IF
      ELSE
        IF (idata_type == iREAL) THEN
          r_out_dp(:)  = 0._wp

          lev_idx = lev
          ! handle the case that a few levels have been selected out of
          ! the total number of levels:
          IF (make_level_selection) THEN
            lev_idx = of%level_selection%global_idx(lev_idx)
          END IF
          CALL exchange_data(in_array=r_ptr(:,lev_idx,:),                 &
            &                out_array=r_out_dp(:), gather_pattern=pat,   &
            &                fill_value = BOUNDARY_MISSVAL)

        ELSE IF (idata_type == iREAL_sp) THEN
          r_out_sp(:)  = 0._wp

          lev_idx = lev
          ! handle the case that a few levels have been selected out of
          ! the total number of levels:
          IF (make_level_selection) THEN
            lev_idx = of%level_selection%global_idx(lev_idx)
          END IF
          CALL exchange_data(in_array=s_ptr(:,lev_idx,:),                 &
            &                out_array=r_out_sp(:), gather_pattern=pat)
          ! FIXME: Implement and use fill_value!
        ELSE IF (idata_type == iINTEGER) THEN
          r_out_int(:) = 0

          lev_idx = lev
          ! handle the case that a few levels have been selected out of
          ! the total number of levels:
          IF (make_level_selection) THEN
            lev_idx = of%level_selection%global_idx(lev_idx)
          END IF
          CALL exchange_data(in_array=i_ptr(:,lev_idx,:),                  &
            &                out_array=r_out_int(:), gather_pattern=pat)
          ! FIXME: Implement and use fill_value!
        END IF
      END IF ! n_glb

      IF (is_mpi_workroot) THEN

        SELECT CASE(idata_type)
        CASE(iREAL)
          !
          ! "r_out_dp" contains double precision data. If single precision
          ! output is desired, we need to perform a type conversion:
          IF ( lwrite_single_precision ) THEN
            r_out_sp(:) = REAL(r_out_dp(:), sp)
          ENDIF

        CASE(iREAL_sp)
          !
          IF ( .NOT. lwrite_single_precision ) THEN
            r_out_dp(:) = REAL(r_out_sp(:), dp)
          ENDIF

        CASE(iINTEGER)
          !
          IF ( lwrite_single_precision ) THEN
            r_out_sp(:) = REAL(r_out_int(:), sp)
          ELSE
            r_out_dp(:) = REAL(r_out_int(:), dp)
          ENDIF
        END SELECT

        ! If required, set lateral boundary points to missing
        ! value. Note that this modifies only the output buffer!
        IF ( info%lmask_boundary                   .AND. &
             & (info%hgrid == GRID_UNSTRUCTURED_CELL) .AND. &
             & ANY(config_lmask_boundary(:)) ) THEN
          missval = BOUNDARY_MISSVAL
          IF (info%lmiss)  missval = info%missval%rval

          IF ( lwrite_single_precision ) THEN
            r_out_sp(1:last_bdry_index) = missval
          ELSE
            r_out_dp(1:last_bdry_index) = missval
          END IF
        END IF

        ! ------------------
        ! case of a test run
        ! ------------------
        !
        ! compare results on worker PEs and test PE
        IF (p_test_run  .AND.  use_dp_mpi2io) THEN
          ! Currently we don't do the check for REAL*4, we would need
          ! p_send/p_recv for this type
          IF (.NOT. is_test) THEN
            ! Send to test PE
            CALL p_send(r_out_dp, process_mpi_all_test_id, 1)
          ELSE IF (p_pe == process_mpi_all_test_id) THEN
            ! Receive result from parallel worker PEs
            CALL p_recv(r_out_recv, process_mpi_all_workroot_id, 1)
            ! check for correctness
            l_error = .FALSE.
            DO i = 1, n_glb
              IF (r_out_recv(i) /= r_out_dp(i)) THEN
                ! do detailed print-out only for "large" errors:
                IF (ABS(r_out_recv(i) - r_out_dp(i)) > SYNC_ERROR_PRINT_TOL) THEN
                  WRITE (0,*) 'Sync error test PE/worker PEs for ', TRIM(info%name)
                  WRITE (0,*) "global pos (", idx_no(i), ",", blk_no(i),")"
                  WRITE (0,*) "level", lev, "/", nlevs
                  WRITE (0,*) "vals: ", r_out_recv(i), r_out_dp(i)
                  l_error = .TRUE.
                END IF
              END IF
            ENDDO
            IF (l_error)   CALL finish(routine,"Sync error!")
          END IF
        END IF

        ! ----------
        ! write data
        ! ----------
        IF (.NOT. is_test) THEN
          IF (.NOT. lwrite_single_precision) THEN
            ! Note for NetCDF: We have already enabled/disabled missing values via vlistDefVarMissVal, since
            !       it is impossible to introduce a FillValue here with nmiss here.
            CALL streamWriteVarSlice (of%cdiFileID, info%cdiVarID, lev-1, r_out_dp(:), nmiss)
          ELSE
            CALL streamWriteVarSliceF(of%cdiFileID, info%cdiVarID, lev-1, r_out_sp(:), nmiss)
          END IF
        END IF

      END IF ! is_mpi_workroot
    END DO ! lev = 1, nlevs

  END SUBROUTINE gather_on_workroot_and_write

  ! Set some GRIB2 keys that may have changed during simulation.
  ! Note that (for synchronous output mode) we provide the
  ! pointer "info_ptr" to the variable's info data object and
  ! not the modified copy "info".
  SUBROUTINE set_time_varying_metadata(of, info, updated_info)
    TYPE (t_output_file), INTENT(IN) :: of
    TYPE(t_var_metadata), INTENT(in) :: info, updated_info

    IF  (of%output_type == FILETYPE_GRB2) THEN
      CALL set_GRIB2_timedep_keys( &
           & of%cdiFileID, info%cdiVarID, updated_info, &
           & of%out_event%output_event%event_data%sim_start,        &
           & get_current_date(of%out_event))
      CALL set_GRIB2_timedep_local_keys(of%cdiFileID, info%cdiVarID, &
           & gribout_config(of%phys_patch_id) )
    END IF
  END SUBROUTINE set_time_varying_metadata


  SUBROUTINE get_bdry_blk_idx(i_log_dom, &
       i_startblk, i_endblk, i_startidx, i_endidx)
    INTEGER, INTENT(in) :: i_log_dom
    INTEGER, INTENT(out) :: i_startblk, i_endblk, i_startidx, i_endidx

    INTEGER  :: rl_start, rl_end, i_nchdom
    TYPE(t_patch), POINTER :: ptr_patch

    ptr_patch => p_patch(i_log_dom)
    rl_start   = 1
    rl_end     = grf_bdywidth_c
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    CALL get_indices_c(ptr_patch, i_endblk, i_startblk, i_endblk, &
      &                i_startidx, i_endidx, rl_start, rl_end)
  END SUBROUTINE get_bdry_blk_idx


  !------------------------------------------------------------------------------------------------
  !> Returns if it is time for the next output step
  !  Please note:
  !  This function returns .TRUE. whenever the next output time of any name list
  !  is reached at the simulation step @p jstep.
  !
  FUNCTION istime4name_list_output(jstep)
    LOGICAL :: istime4name_list_output
    INTEGER, INTENT(IN)   :: jstep            ! simulation time step
    ! local variables
    INTEGER :: i
    LOGICAL :: ret

    ret = .FALSE.
    IF (ALLOCATED(output_file)) THEN
       ! note: there may be cases where no output namelist has been
       ! defined. thus we must check if "output_file" has been
       ! allocated.
       DO i = 1, SIZE(output_file)
         ret = ret .OR. is_output_step(output_file(i)%out_event, jstep)
         IF (ret) EXIT
       END DO
    END IF
    istime4name_list_output = ret
  END FUNCTION istime4name_list_output


 !------------------------------------------------------------------------------------------------
  !> Returns if it is time for output of a particular variable at a particular output step on a particular domain
  !  Please note:
  !  This function returns .TRUE. whenever the variable is due for output in any name list
  !  at the simulation step @p jstep for the domain @p jg.
  !
  FUNCTION istime4name_list_output_dom(jg, jstep)
    LOGICAL                      :: istime4name_list_output_dom
    INTEGER, INTENT(IN)          :: jg         !< domain index
    INTEGER, INTENT(IN)          :: jstep      !< simulation time step
    ! local variables
    INTEGER :: i
    LOGICAL :: ret
    TYPE(t_output_name_list), POINTER :: p_onl      !< output name list

    ret = .FALSE.
    IF (ALLOCATED(output_file)) THEN
       ! note: there may be cases where no output namelist has been
       ! defined. thus we must check if "output_file" has been
       ! allocated.
       DO i = 1, SIZE(output_file)
         p_onl => output_file(i)%name_list
         ret = ret .OR. &
           ( is_output_step(output_file(i)%out_event, jstep) .AND. &
             ( p_onl%dom == jg .OR. p_onl%dom == -1 )              &
           )
         IF (ret) EXIT
       END DO
    END IF
    istime4name_list_output_dom = ret
  END FUNCTION istime4name_list_output_dom


  !------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------
  ! The following routines are only needed for asynchronous IO
  !------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------

  ! Just define the entry point of name_list_io_main_proc, it will never be called

  SUBROUTINE name_list_io_main_proc(sim_step_info)
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
  END SUBROUTINE name_list_io_main_proc


END MODULE mo_name_list_output
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
