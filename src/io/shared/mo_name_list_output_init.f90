! Module handling the initialization of synchronous and asynchronous output.
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
!NEC$ options "-fno-loop-unroll"
MODULE mo_name_list_output_init

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_f_pointer, c_int64_t, c_double

  ! constants and global settings
  USE mo_cdi,                               ONLY: FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB2, gridCreate, cdiEncodeDate,         &
                                                & cdiEncodeTime, institutInq, vlistCreate, TUNIT_MINUTE, CDI_UNDEFID,           &
                                                & TAXIS_RELATIVE, taxisCreate, TAXIS_ABSOLUTE, GRID_UNSTRUCTURED, GRID_LONLAT,  &
                                                & gridDefPosition, gridDefXsize, gridDefXname, gridDefXunits, gridDefYsize,     &
                                                & gridDefYname, gridDefYunits, gridDefNumber, gridDefUUID, vlistDefInstitut,    &
                                                & gridDefNvertex, gridDefXvals, cdiDefAttTxt, CDI_GLOBAL, gridDefParamRLL,      &
                                                & gridDefYvals, gridDefXlongname, gridDefYlongname, gridDefReference,           &
                                                & taxisDefTunit, taxisDefCalendar, taxisDefRdate, taxisDefRtime, vlistDefTaxis, &
                                                & gridDefProj, GRID_PROJECTION, GRID_CURVILINEAR, institutDef
  USE mo_cdi_constants,                     ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE, &
                                                & GRID_REGULAR_LONLAT, GRID_VERTEX, GRID_EDGE, GRID_CELL, GRID_ZONAL
  USE mo_dynamics_config, ONLY: nnow, nnew, nold
  USE mo_kind,                              ONLY: wp, i8, dp, sp
  USE mo_impl_constants,                    ONLY: max_phys_dom, max_dom, SUCCESS, vname_len,         &
    &                                             max_var_ml, max_var_pl, max_var_hl, max_var_il,    &
    &                                             MAX_CHAR_LENGTH, MAX_NUM_IO_PROCS, nlat_moc, INWP, &
    &                                             MAX_TIME_INTERVALS, MAX_NPLEVS, pio_type_cdipio,   &
    &                                             MAX_NZLEVS, MAX_NILEVS, BOUNDARY_MISSVAL,          &
    &                                             dtime_proleptic_gregorian => proleptic_gregorian,  &
    &                                             dtime_cly360              => cly360
  USE mo_cdi_constants,                     ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT,            &
    &                                             GRID_UNSTRUCTURED_EDGE, GRID_REGULAR_LONLAT, GRID_VERTEX,  &
    &                                             GRID_EDGE, GRID_CELL, LONLAT_PREFIX
  USE mo_io_units,                          ONLY: filename_max, nnml, nnml_output
  USE mo_master_config,                     ONLY: getModelBaseDir, isRestart
  USE mo_master_control,                    ONLY: my_process_is_atmo, my_process_is_oceanic, my_process_is_jsbach, &
    &                                             my_process_is_waves
  ! basic utility modules
  USE mo_exception,                         ONLY: finish, message, message_text
  USE mo_dictionary,                        ONLY: t_dictionary
  USE mo_io_util,                           ONLY: get_file_extension
  USE mo_util_cdi,                          ONLY: create_cdi_variable
  USE mo_util_string,                       ONLY: t_keyword_list, associate_keyword,              &
    &                                             with_keywords, insert_group,                    &
    &                                             tolower, int2string, difference,                &
    &                                             sort_and_compress_list,                         &
    &                                             real2string, remove_whitespace,                 &
    &                                             lowcase, remove_duplicates
  USE mo_util_texthash,                         ONLY: text_hash_c
  USE mo_cf_convention,                     ONLY: t_cf_var, cf_global_info
  USE mo_restart_nml_and_att,               ONLY: getAttributesForRestarting
  USE mo_key_value_store,                   ONLY: t_key_value_store
  USE mo_model_domain,                      ONLY: p_patch, p_phys_patch
  USE mo_math_utilities,                    ONLY: merge_values_into_set, t_value_set
  USE mo_math_constants,                    ONLY: rad2deg
  ! config modules
  USE mo_parallel_config,                   ONLY: nproma, p_test_run, use_dp_mpi2io, pio_type
  USE mo_run_config,                        ONLY: dtime, msg_level, output_mode,                  &
    &                                             ICON_grid_file_uri, number_of_grid_used, iforcing
  USE mo_grid_config,                       ONLY: n_dom, n_phys_dom, start_time, end_time,        &
    &                                             DEFAULT_ENDTIME
  USE mo_io_config,                         ONLY: netcdf_dict, output_nml_dict, linvert_dict,     &
    &                                             config_lmask_boundary => lmask_boundary
  USE mo_name_list_output_config,           ONLY: use_async_name_list_io,                         &
    &                                             first_output_name_list
  USE mo_time_config,                       ONLY: time_config
  USE mo_gribout_config,                    ONLY: gribout_config


  USE mo_nh_pzlev_config,                   ONLY: nh_pzlev_config
  USE mo_extpar_config,                     ONLY: i_lctype
  USE mo_lnd_nwp_config,                    ONLY: ntiles_lnd, ntiles_water, ntiles_total, tile_list, &
    &                                             isub_water, isub_lake, isub_seaice, lsnowtile
  USE mo_nwp_sfc_tiles,                     ONLY: setup_tile_list
  USE mo_meteogram_config,                  ONLY: meteogram_output_config

  ! MPI Communication routines
  USE mo_mpi,                               ONLY: p_bcast, p_comm_work, p_comm_work_2_io,         &
    &                                             p_comm_io, p_comm_work_io,                      &
    &                                             MPI_COMM_NULL, MPI_COMM_SELF,                   &
    &                                             p_send, p_recv,                                 &
    &                                             p_real_dp, p_real_sp,                           &
    &                                             my_process_is_stdio, my_process_is_mpi_test,    &
    &                                             my_process_is_mpi_workroot,                     &
    &                                             my_process_is_io, my_process_is_work,           &
    &                                             my_process_is_mpi_ioroot,                       &
    &                                             process_work_io0, p_allgather,        &
    &                                             process_mpi_io_size, p_n_work,  &
    &                                             p_pe_work, p_io_pe0, p_work_pe0, p_pe



  USE mo_communication,                     ONLY: idx_no, blk_no
  ! namelist handling
  USE mo_namelist,                          ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_nml_annotate,                      ONLY: temp_defaults, temp_settings
  ! variable lists
  USE mo_var_groups,                        ONLY: var_groups_dyn, MAX_GROUPS
  USE mo_var_metadata_types,                ONLY: t_var_metadata
  USE mo_var_list_register_utils,           ONLY: vlr_group, vlr_replicate
  USE mo_var_list_register,                 ONLY: t_vl_register_iter
  USE mo_var,                               ONLY: t_var
  USE mo_var_metadata,                      ONLY: get_var_timelevel, get_var_name
  USE mo_var, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  ! lon-lat interpolation
  USE mo_lonlat_grid,                       ONLY: t_lon_lat_grid, compute_lonlat_blocking,        &
    &                                             compute_lonlat_specs, threshold_delta_or_intvls,&
    &                                             rotate_latlon_grid
  USE mo_intp_lonlat_types,                 ONLY: t_lon_lat_intp, t_lon_lat_data, lonlat_grids
  ! output events
  USE mtime,                                ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN,    &
    &                                             timedelta, newTimedelta, deallocateTimedelta,   &
    &                                             OPERATOR(<), newDatetime, deallocateDatetime,   &
    &                                             getTotalMilliSecondsTimeDelta, datetime,        &
    &                                             OPERATOR(+), datetimeToString, OPERATOR(>),     &
    &                                             timedeltaToString, calendarType,                &
    &                                             getPTStringFromSeconds, OPERATOR(/=),           &
    &                                             mtime_proleptic_gregorian => proleptic_gregorian, &
    &                                             mtime_year_of_360_days => year_of_360_days
  USE mo_util_mtime,                        ONLY: mtime_timedelta_from_fseconds
  USE mo_output_event_types,                ONLY: t_sim_step_info, MAX_EVENT_NAME_STR_LEN,        &
    &                                             DEFAULT_EVENT_NAME, t_par_output_event,         &
    &                                             max_filename_str_len
  USE mo_output_event_control,              ONLY: compute_matching_sim_steps,                     &
    &                                             generate_output_filenames
  USE mo_output_event_handler,              ONLY: new_parallel_output_event,                      &
    &                                             union_of_all_events,      &
    &                                             print_output_event,                             &
    &                                             set_event_to_simstep



  ! name list output
  USE mo_reorder_info,                      ONLY: t_reorder_info, &
    mask2reorder_info, transfer_reorder_info
  USE mo_name_list_output_types,            ONLY: l_output_phys_patch, t_output_name_list,        &
    &                                             t_output_file, t_var_desc, t_grid_info,         &
    &                                             t_patch_info, icell, iedge, ivert,              &
    &                                             REMAP_NONE, REMAP_REGULAR_LATLON,               &
    &                                             GRP_PREFIX, TILE_PREFIX,                        &
    &                                             t_fname_metadata, all_events, t_patch_info_ll,  &
    &                                             is_grid_info_var, GRB2_GRID_INFO_NAME,          &
    &                                             t_event_data_local
  USE mo_name_list_output_gridinfo,         ONLY: set_grid_info_grb2, set_grid_info_netcdf,       &
    &                                             collect_all_grid_info, copy_grid_info,          &
    &                                             allgather_grid_info, deallocate_all_grid_info,  &
    &                                             GRID_INFO_NONE, GRID_INFO_FILE, GRID_INFO_BCAST
  USE mo_name_list_output_metadata,         ONLY: metainfo_allocate_memory_window
  USE mo_name_list_output_zaxes,            ONLY: setup_ml_axes_atmo, setup_pl_axis_atmo,         &
    &                                             setup_hl_axis_atmo, setup_il_axis_atmo,         &
    &                                             setup_zaxes_oce

  USE mo_waves_vertical_axes,               ONLY: setup_zaxes_waves

  USE mo_level_selection_types,             ONLY: t_level_selection

  USE mo_atm_phy_nwp_config,                ONLY: atm_phy_nwp_config
  USE mo_impl_constants,                    ONLY: LSS_JSBACH
  USE mo_aes_phy_config,                    ONLY: aes_phy_config
  USE mo_jsb_vertical_axes,                 ONLY: setup_zaxes_jsbach

  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxisList, t_verticalAxis
  USE mo_name_list_output_printvars,        ONLY: print_var_list
  USE mo_util_vgrid_types,                  ONLY: vgrid_buffer
  USE mo_derived_variable_handling,         ONLY: init_statistics


  USE mo_vertical_coord_table,              ONLY: vct



  USE mo_coupling_config,                   ONLY: is_coupled_run
  USE mo_master_control,                    ONLY: get_my_process_name
  IMPLICIT NONE

  PRIVATE

  ! variables and data types
  PUBLIC :: out_varnames_dict
  PUBLIC :: varnames_dict
  PUBLIC :: output_file
  PUBLIC :: patch_info
  PUBLIC :: lonlat_info
  PUBLIC :: zonal_ri, profile_ri
  ! subroutines
  PUBLIC :: read_name_list_output_namelists
  PUBLIC :: parse_variable_groups
  PUBLIC :: init_name_list_output
  PUBLIC :: setup_output_vlist
  PUBLIC :: collect_requested_ipz_levels
  PUBLIC :: create_vertical_axes
  PUBLIC :: isRegistered
  PUBLIC :: nlevs_of_var


  !------------------------------------------------------------------------------------------------

  TYPE(t_output_file),   ALLOCATABLE, TARGET :: output_file(:)
  TYPE(t_patch_info),    ALLOCATABLE, TARGET :: patch_info (:)
  TYPE(t_patch_info_ll), ALLOCATABLE, TARGET :: lonlat_info(:,:)
  TYPE(t_key_value_store) :: outputRegister
  TYPE(t_reorder_info), TARGET :: zonal_ri, profile_ri

  ! Number of output domains. This depends on l_output_phys_patch and is either the number
  ! of physical or the number of logical domains.
  INTEGER :: n_dom_out

  ! Broadcast root for intercommunicator broadcasts form compute PEs to IO PEs using p_comm_work_2_io
  INTEGER :: bcast_root

  !NEC_RP: implementation for hybrid-MPI
  INTEGER :: compute_dt
  PUBLIC  :: compute_dt

  !------------------------------------------------------------------------------------------------
  ! dictionaries for variable names:
  TYPE (t_dictionary) ::     &
    &   varnames_dict      , & !< maps variable names onto the internal ICON names.
    &   out_varnames_dict      !< maps internal variable names onto names in output file (NetCDF only).
  !------------------------------------------------------------------------------------------------

  CHARACTER(*), PARAMETER :: modname = 'mo_name_list_output_init'

CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Read configuration for namelist controlled output module.
  !
  !  The name is a bit strange, but if follows the convention to read
  !  namelists with a routine called read_XXX_namelist (plural is
  !  used here since several namelists are read).
  !
  !  Please note the difference between name_list and namelist!
  !
  SUBROUTINE read_name_list_output_namelists( filename )
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    ! local variables
    CHARACTER(LEN=*), PARAMETER       :: routine = modname//'::read_name_list_output_namelists'

    INTEGER                               :: istat, idom
    TYPE(t_output_name_list), POINTER     :: p_onl
    LOGICAL                               :: lrewind

    ! Local variables corresponding to members of output_name_list
    INTEGER                               :: filetype
    INTEGER                               :: mode
    INTEGER                               :: taxis_tunit
    INTEGER                               :: dom(max_phys_dom)
    INTEGER                               :: steps_per_file
    LOGICAL                               :: steps_per_file_inclfirst         !< Flag. Do not count first step in files count
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: file_interval                    !< length of a file (ISO8601 duration)
    LOGICAL                               :: include_last
    LOGICAL                               :: output_grid
    CHARACTER(LEN=filename_max)           :: output_filename
    CHARACTER(LEN=filename_max)           :: filename_format
    CHARACTER(LEN=filename_max)           :: filename_extn                    !< user-specified filename extension (or "default")
    CHARACTER(LEN=vname_len)              :: ml_varlist(max_var_ml)
    CHARACTER(LEN=vname_len)              :: pl_varlist(max_var_pl)
    CHARACTER(LEN=vname_len)              :: hl_varlist(max_var_hl)
    CHARACTER(LEN=vname_len)              :: il_varlist(max_var_il)
    CHARACTER(len=MAX_CHAR_LENGTH)        :: m_levels                         !< level selection: model levels
    REAL(wp)                              :: p_levels(MAX_NPLEVS)             !< pressure levels
    REAL(wp)                              :: h_levels(MAX_NZLEVS)             !< height levels
    REAL(wp)                              :: i_levels(MAX_NILEVS)             !< isentropic levels
    INTEGER                               :: remap
    CHARACTER(LEN=MAX_CHAR_LENGTH)        :: operation
    REAL(wp)                              :: reg_lon_def(3)
    REAL(wp)                              :: reg_lat_def(3)
    INTEGER                               :: reg_def_mode
    REAL(wp)                              :: north_pole(2)

    REAL(wp)                              :: output_bounds(3*MAX_TIME_INTERVALS)
    INTEGER                               :: output_time_unit
    CHARACTER(LEN=MAX_DATETIME_STR_LEN+1) :: output_start(MAX_TIME_INTERVALS), &
      &                                      output_end(MAX_TIME_INTERVALS)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: output_interval(MAX_TIME_INTERVALS)
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN) :: ready_file  !< ready filename prefix (=output event name)

    TYPE(t_lon_lat_data),  POINTER        :: lonlat
    TYPE (t_keyword_list), POINTER        :: keywords
    CHARACTER(len=MAX_CHAR_LENGTH)        :: cfilename
    INTEGER                               :: iunit, lonlat_id, jg
    TYPE (t_lon_lat_grid)                 :: new_grid

    !> "stream_partitions": Split one namelist into concurrent,
    !> alternating files:
    INTEGER                               :: stream_partitions_ml, &
      &                                      stream_partitions_pl, &
      &                                      stream_partitions_hl, &
      &                                      stream_partitions_il

    !> MPI ranks which were explicitly specified by the user:
    INTEGER                               :: pe_placement_ml(MAX_NUM_IO_PROCS), &
      &                                      pe_placement_pl(MAX_NUM_IO_PROCS), &
      &                                      pe_placement_hl(MAX_NUM_IO_PROCS), &
      &                                      pe_placement_il(MAX_NUM_IO_PROCS)

    !> RBF shape parameter.
    REAL(wp)                              :: rbf_scale
    LOGICAL :: is_stdio

    ! The namelist containing all variables above
    NAMELIST /output_nml/ &
      mode, taxis_tunit, dom,                                &
      filetype, filename_format, output_filename,            &
      steps_per_file, steps_per_file_inclfirst,              &
      include_last, output_grid,                             &
      remap, reg_lon_def, reg_lat_def, north_pole,           &
      ml_varlist, pl_varlist, hl_varlist, il_varlist,        &
      m_levels, p_levels, h_levels, i_levels,                &
      output_start, output_end, output_interval,             &
      output_bounds, output_time_unit,                       &
      ready_file, file_interval, reg_def_mode,               &
      stream_partitions_ml, stream_partitions_pl,            &
      stream_partitions_hl, stream_partitions_il,            &
      pe_placement_ml, pe_placement_pl,                      &
      pe_placement_hl, pe_placement_il,                      &
      filename_extn, rbf_scale, operation

    ! Before we start: prepare the levels set objects for the vertical
    ! interpolation.
    !
    ! level ordering: zlevels, plevels, ilevels must be ordered from
    ! TOA to bottom:
    !
    DO jg = 1, max_dom
      ! pressure levels:
      nh_pzlev_config(jg)%plevels%sort_smallest_first = .TRUE.
      ! height levels
      nh_pzlev_config(jg)%zlevels%sort_smallest_first = .FALSE.
      ! isentropic levels
      nh_pzlev_config(jg)%ilevels%sort_smallest_first = .FALSE.
    END DO

    ! create variable of registering output variables. should be used later for
    ! triggering computation only in case of output request
    CALL outputRegister%init(.false.)

    ! -- Open input file and position to first namelist 'output_nml'

    CALL open_nml(TRIM(filename))

    ! As in COSMO, there may exist several output_nml namelists in the input file
    ! Loop until EOF is reached

    p_onl                  => NULL()
    first_output_name_list => NULL()
    lrewind                = .TRUE.

    is_stdio = my_process_is_stdio()

    IF (.NOT. output_mode%l_nml) RETURN ! do not read output namelists if main switch is set to false

    DO
      CALL position_nml ('output_nml', lrewind=lrewind, status=istat)
      IF(istat /= POSITIONED) THEN

        ! if no "output_nml" has been found at all, we still cannot
        ! disable this mode: There might be meteograms which have to
        ! be handled by the output PEs.  If meteograms are disabled,
        ! too, then the user's namelist settings were inconsistent.
        !
        IF (.NOT. ASSOCIATED(first_output_name_list)) THEN
          IF (.NOT. ANY(meteogram_output_config(:)%lenabled)) THEN
            CALL message(routine, "No output definition found; disabling output.")
            output_mode%l_nml = .FALSE.
          ELSE
            CALL message(routine, "No output definition found; meteogram output only.")
          END IF
        END IF

        CALL close_nml
        RETURN
      ENDIF
      lrewind = .FALSE.

      ! -- Set all variables in output_nml to their default values

      filetype                 = FILETYPE_NC2 ! NetCDF
      mode                     = 2
      taxis_tunit              = TUNIT_MINUTE
      dom(:)                   = -1
      steps_per_file           = -1
      steps_per_file_inclfirst = .TRUE.
      file_interval            = ' '
      include_last             = .TRUE.
      output_grid              = .FALSE.
      output_filename          = ' '
      filename_format          = "<output_filename>_DOM<physdom>_<levtype>_<jfile>"
      filename_extn            = "default"
      ml_varlist(:)            = ' '
      pl_varlist(:)            = ' '
      hl_varlist(:)            = ' '
      il_varlist(:)            = ' '
      m_levels                 = " "
      p_levels(:)              = -1._wp
      h_levels(:)              = -1._wp
      i_levels(:)              = -1._wp
      remap                    = REMAP_NONE
      operation                = ''
      reg_lon_def(:)           = 0._wp
      reg_lat_def(:)           = 0._wp
      reg_def_mode             = 0
      north_pole(:)            = (/ 0._wp, 90._wp /)
      output_start(:)          = ' '
      output_end(:)            = ' '
      output_interval(:)       = ' '
      output_bounds(:)         = -1._wp
      output_time_unit         = 1
      ready_file               = DEFAULT_EVENT_NAME
      lonlat_id                = 0
      stream_partitions_ml     = 1
      stream_partitions_pl     = 1
      stream_partitions_hl     = 1
      stream_partitions_il     = 1
      pe_placement_ml(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      pe_placement_pl(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      pe_placement_hl(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      pe_placement_il(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      rbf_scale                = -1._wp

      ! -- Read output_nml

      IF (is_stdio)  THEN
        iunit = temp_defaults()
        WRITE(iunit, output_nml)                                     ! write defaults to temporary text file
      END IF
      READ (nnml, output_nml)                          ! overwrite default settings
      IF (is_stdio)  THEN
        iunit = temp_settings()
        WRITE(iunit, output_nml)                                     ! write settings to temporary text file
      END IF

      ! -- Consistency checks:

      IF ((steps_per_file == -1) .AND. (LEN_TRIM(file_interval) == 0)) THEN
        CALL finish(routine, "Please specify either <steps_per_file> or <file_interval>!")
      END IF
      IF ((steps_per_file /= -1) .AND. (LEN_TRIM(file_interval) /= 0)) THEN
        CALL finish(routine, "User has specified conflicting parameters <steps_per_file>, <file_interval>!")
      END IF
      IF(remap/=REMAP_NONE .AND. remap/=REMAP_REGULAR_LATLON) THEN
        CALL finish(routine,'Unsupported value for remap')
      END IF
      IF ((remap==REMAP_REGULAR_LATLON) .AND. ALL(reg_lon_def(:) == 0.) .AND. ALL(reg_lat_def(:) == 0.)) THEN
        CALL finish(routine,'Lon-lat output: Grid not specified in namelist!')
      END IF
      IF ((reg_lon_def(3) >  reg_lon_def(1)) .AND. &
        & (reg_lon_def(2) <= 0._wp)) THEN
        CALL finish(routine,'Illegal LON increment')
      END IF
      IF ((reg_lat_def(3) /= reg_lat_def(1)) .AND. &
        & (reg_lat_def(2) == 0._wp)) THEN
        CALL finish(routine,'Illegal LAT increment')
      END IF
      IF(reg_lon_def(3)<reg_lon_def(1)) CALL finish(routine,'end lon < start lon')

      ! -- Scale output bounds

      IF (output_bounds(1) > -1._wp) THEN
        ! output_bounds is always in seconds - the question is what to do with months or years
        ! output_time_unit: 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
        SELECT CASE(output_time_unit)
        CASE(1); output_bounds(:) = output_bounds(:)
        ! Note: output_bounds == -1 is used later on to check for valid entries
        CASE(2); WHERE(output_bounds(:) >= 0._wp) output_bounds(:) = output_bounds(:)*60._wp
        CASE(3); WHERE(output_bounds(:) >= 0._wp) output_bounds(:) = output_bounds(:)*3600._wp
        CASE(4); WHERE(output_bounds(:) >= 0._wp) output_bounds(:) = output_bounds(:)*86400._wp
        CASE(5); WHERE(output_bounds(:) >= 0._wp) output_bounds(:) = output_bounds(:)*86400._wp*30._wp  ! Not a real calender month
        CASE(6); WHERE(output_bounds(:) >= 0._wp) output_bounds(:) = output_bounds(:)*86400._wp*365._wp ! Not a real calender year
        CASE DEFAULT
          CALL finish(routine,'Illegal output_time_unit')
        END SELECT
      END IF

      ! -- Read the map files into dictionary data structures

      CALL varnames_dict%init(.FALSE.)
      CALL out_varnames_dict%init(.FALSE.)

      NULLIFY(keywords)
      CALL associate_keyword("<path>", TRIM(getModelBaseDir()), keywords)
      IF(output_nml_dict     /= ' ') THEN
        cfilename = with_keywords(keywords, output_nml_dict)
        CALL message(routine, "load dictionary file.")
        CALL varnames_dict%loadfile(cfilename, linverse=linvert_dict)
      END IF
      IF(netcdf_dict /= ' ') THEN
        cfilename = with_keywords(keywords, netcdf_dict)
        CALL message(routine, "load dictionary file (output names).")
        CALL out_varnames_dict%loadfile(cfilename, linverse=.TRUE.)
      END IF

      ! -- If "remap=1": lon-lat interpolation requested

      IF (remap == REMAP_REGULAR_LATLON) THEN
        ! define lon-lat grid specification
        new_grid%reg_lon_def(:) = reg_lon_def(:)
        new_grid%reg_lat_def(:) = reg_lat_def(:)
        new_grid%north_pole(:)  = north_pole(:)
        new_grid%reg_def_mode   = reg_def_mode
        ! compute some additional entries of lon-lat grid specification:
        CALL compute_lonlat_specs(new_grid)
        CALL compute_lonlat_blocking(new_grid, nproma)
        ! If the user has explicitly specified an interpolation
        ! parameter, then we always register this as a new lon-lat
        ! grid. Otherwise we might share the lon-lat coefficients with
        ! other output namelists.
        IF (rbf_scale > 0._wp) THEN
          lonlat_id             =  lonlat_grids%add_new_grid()
          lonlat                => lonlat_grids%list(lonlat_id)
          CALL lonlat%init(new_grid, rbf_scale)
        ELSE
          ! check, if lon-lat grids has already been registered
          lonlat_id = lonlat_grids%get_ID(new_grid)
          IF (lonlat_id == -1) THEN
            ! Register a lon-lat grid data structure in global list
            lonlat_id             =  lonlat_grids%add_new_grid()
            lonlat                => lonlat_grids%list(lonlat_id)
            CALL lonlat%init(new_grid, rbf_scale)
          ELSE
            lonlat => lonlat_grids%list(lonlat_id)
          END IF
        END IF

        ! Flag those domains, which are used for this lon-lat grid:
        !     If dom(:) was not specified in namelist input, it is set
        !     completely to -1.  In this case all domains are wanted in
        !     the output
        IF (dom(1) < 0)  THEN
          lonlat%l_dom(:) = .TRUE.
        ELSE
          DOM_LOOP : DO idom = 1, max_dom
            IF (dom(idom) < 0) EXIT DOM_LOOP
            lonlat%l_dom( dom(idom) ) = .TRUE.
          ENDDO DOM_LOOP
        END IF
      ENDIF
! #ifndef __NO_ICON_ATMO__

      ! loop over domains and generate a separare copy of the namelist
      ! settings for each chosen domain.
      !
      DOM_LOOP2 : DO idom = 1, max_dom
        IF ((dom(idom) < 0) .AND. (dom(1) >= 0)) EXIT DOM_LOOP2

        ! -- Allocate next output_name_list

        IF(.NOT.ASSOCIATED(first_output_name_list)) THEN
          ! Allocate first name_list
          ALLOCATE(first_output_name_list)
          p_onl => first_output_name_list
        ELSE
          ! This is not the first one, p_onl points to the last one which was created
          ALLOCATE(p_onl%next)
          p_onl => p_onl%next
        ENDIF

        ! -- Set next output_name_list from values read

        p_onl%filetype                 = filetype
        p_onl%mode                     = mode
        p_onl%taxis_tunit              = taxis_tunit

        ! If dom(:) was not specified in namelist input, it is set
        ! completely to -1.  In this case all domains are wanted in
        ! the output.
        IF (dom(idom) < 0) THEN
          p_onl%dom                    = idom
        ELSE
          p_onl%dom                    = dom(idom)
        END IF

        p_onl%steps_per_file           = steps_per_file
        ! conditional default: steps_per_file_inclfirst=.FALSE. for GRIB output
        p_onl%steps_per_file_inclfirst = steps_per_file_inclfirst .AND. (filetype /= FILETYPE_GRB2)
        p_onl%file_interval            = file_interval
        p_onl%include_last             = include_last
        p_onl%output_grid              = output_grid
        p_onl%output_filename          = output_filename
        p_onl%filename_format          = filename_format
        p_onl%filename_extn            = filename_extn
        p_onl%ml_varlist(:)            = ml_varlist(:)
        p_onl%pl_varlist(:)            = pl_varlist(:)
        p_onl%hl_varlist(:)            = hl_varlist(:)
        p_onl%il_varlist(:)            = il_varlist(:)
        p_onl%m_levels                 = m_levels
        p_onl%p_levels                 = p_levels
        p_onl%z_levels                 = h_levels
        p_onl%i_levels                 = i_levels
        p_onl%remap                    = remap
        p_onl%operation                = operation
        p_onl%output_start(:)          = output_start(:)
        p_onl%output_end(:)            = output_end
        p_onl%output_interval(:)       = output_interval
        p_onl%output_bounds(:)         = output_bounds(:)
        p_onl%ready_file               = ready_file
        p_onl%lonlat_id                = lonlat_id
        p_onl%stream_partitions_ml     = stream_partitions_ml
        p_onl%stream_partitions_pl     = stream_partitions_pl
        p_onl%stream_partitions_hl     = stream_partitions_hl
        p_onl%stream_partitions_il     = stream_partitions_il
        p_onl%pe_placement_ml(:)       = pe_placement_ml(:)
        p_onl%pe_placement_pl(:)       = pe_placement_pl(:)
        p_onl%pe_placement_hl(:)       = pe_placement_hl(:)
        p_onl%pe_placement_il(:)       = pe_placement_il(:)

        ! -- translate variables names according to variable name dictionary:
        ! allow case-insensitive variable names:
        ! -- if the namelist switch "output_grid" has been enabled: add
        !    "clon, "clat", "elon", "elat", etc. to the list of variables:
        IF (p_onl%output_grid) THEN
          CALL lookup_lc_gridout(p_onl%ml_varlist, p_onl%remap)
          CALL lookup_lc_gridout(p_onl%pl_varlist, p_onl%remap)
          CALL lookup_lc_gridout(p_onl%hl_varlist, p_onl%remap)
          CALL lookup_lc_gridout(p_onl%il_varlist, p_onl%remap)
        ELSE
          CALL lookup_lc_gridout(p_onl%ml_varlist)
          CALL lookup_lc_gridout(p_onl%pl_varlist)
          CALL lookup_lc_gridout(p_onl%hl_varlist)
          CALL lookup_lc_gridout(p_onl%il_varlist)
        END IF
        p_onl%next => NULL()
      END DO DOM_LOOP2

      ! -- write the contents of the namelist to an ASCII file
      IF (is_stdio) WRITE(nnml_output,nml=output_nml)

    ENDDO

    CALL close_nml
  CONTAINS

    SUBROUTINE lookup_lc_gridout(varlist, remap)
      CHARACTER(LEN=vname_len), INTENT(INOUT) :: varlist(:)
      INTEGER, INTENT(IN), OPTIONAL :: remap
      INTEGER :: k

      DO k = 1, SIZE(varlist)
        IF (' ' == varlist(k)) EXIT ! since read from nml-file array is filled bottom to top...
        varlist(k) = varnames_dict%get(varlist(k), default=varlist(k))
        varlist(k) = tolower(varlist(k))
      END DO
      IF (PRESENT(remap)) THEN 
        IF (k .GT. 1) THEN
          SELECT CASE(remap)
          CASE (REMAP_NONE)
            IF (k+5 .GT. SIZE(varlist)) &
              & CALL finish(routine, "Insufficient array size!")
            varlist(k  ) = tolower(GRB2_GRID_INFO_NAME(1,1))
            varlist(k+1) = tolower(GRB2_GRID_INFO_NAME(1,2))
            varlist(k+2) = tolower(GRB2_GRID_INFO_NAME(2,1))
            varlist(k+3) = tolower(GRB2_GRID_INFO_NAME(2,2))
            varlist(k+4) = tolower(GRB2_GRID_INFO_NAME(3,1))
            varlist(k+5) = tolower(GRB2_GRID_INFO_NAME(3,2))
          CASE (REMAP_REGULAR_LATLON)
            IF (k+1 .GT. SIZE(varlist)) &
              & CALL finish(routine, "Insufficient array size!")
            varlist(k  ) = tolower(GRB2_GRID_INFO_NAME(0,1))
            varlist(k+1) = tolower(GRB2_GRID_INFO_NAME(0,2))
          END SELECT
        END IF
      END IF
    END SUBROUTINE lookup_lc_gridout
  END SUBROUTINE read_name_list_output_namelists

  !------------------------------------------------------------------------------------------------
  !> Appends the chosen p-levels, z-levels, i-levels to the levels
  !  sets for the corresponding domains (note that we do the vertical
  !  interpolation for the union of all chosen levels and only do a
  !  selection for each output namelist):
  !
  SUBROUTINE collect_requested_ipz_levels()
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::collect_requested_ipz_levels"
    TYPE(t_output_name_list), POINTER  :: onl
    INTEGER :: n_dom_out, jp, id

    ! Loop over the output namelists and create a union set of all requested vertical levels (per domain):
    onl => first_output_name_list
    n_dom_out = MERGE(n_phys_dom, n_dom, l_output_phys_patch)
    DO WHILE (ASSOCIATED(onl))
      DO jp = 1, n_dom_out
        ! append the chosen p-levels, z-levels, i-levels
        id = MERGE(p_phys_patch(jp)%logical_id, jp, l_output_phys_patch)
        IF (onl%dom .NE. id) CYCLE
        CALL merge_set(onl%p_levels, "pressure",   onl%pl_varlist(1), nh_pzlev_config(id)%plevels)
        CALL merge_set(onl%z_levels, "height",     onl%hl_varlist(1), nh_pzlev_config(id)%zlevels)
        CALL merge_set(onl%i_levels, "isentropic", onl%il_varlist(1), nh_pzlev_config(id)%ilevels)
      END DO
      onl => onl%next
    END DO ! onl
  CONTAINS

  SUBROUTINE merge_set(levels, intp_name, vl1, mergeset)
    REAL(wp), INTENT(IN) :: levels(:)
    CHARACTER(*), INTENT(IN) :: intp_name, vl1
    TYPE(t_value_set), INTENT(INOUT) :: mergeset
    INTEGER :: nlevs

    DO nlevs=1,SIZE(levels)
      IF (levels(nlevs) < 0._wp) EXIT
    END DO
    nlevs = nlevs - 1
    IF ((nlevs == 0) .AND. (vl1 /= ' ')) &
      & CALL finish(routine, "Input from output_nml: requested " // &
        & intp_name // "interpolation without specifying levels!")
    IF (nlevs > 0)  CALL merge_values_into_set(nlevs, levels, mergeset)
  END SUBROUTINE merge_set
  END SUBROUTINE collect_requested_ipz_levels


  !------------------------------------------------------------------------------------------------
  !> Looks for variable groups ("group:xyz", "tiles:xyz") and replaces them
  !
  !  @note This subroutine cannot be called directly from
  !         "read_name_list_output_namelists" because the latter
  !         routine is processed earlier, before all variable lists
  !         have been registered through "add_vars".
  !
  !  @note In more detail, this subroutine looks for variable groups
  !        ("group:xyz","tiles:xyz") and replaces them by all variables
  !        belonging to the group. Afterwards, variables can be REMOVED
  !        from this union set with the syntax "-varname". Note that typos
  !        are not detected but that the corresponding variable is
  !        simply not removed!
  !
  SUBROUTINE parse_variable_groups()
    CHARACTER(LEN=vname_len), ALLOCATABLE :: varlist(:), grp_vars(:)
    CHARACTER(LEN=vname_len)              :: vname, grp_name
    INTEGER :: nvars, ngrp_vars, i_typ, ierrstat, ivar, ntotal_vars, &
      & jvar, i, nsubtract_vars, tlen, ninserted
    CHARACTER(LEN=vname_len),  POINTER      :: in_varlist(:)
    TYPE(t_output_name_list), POINTER      :: p_onl
    TYPE(t_vl_register_iter) :: vl_iter
    CHARACTER(*), PARAMETER :: routine = modname//"::parse_variable_groups"

    ntotal_vars = 0
    DO WHILE(vl_iter%next())
      DO i = 1, vl_iter%cur%p%nvars
        IF (.NOT.vl_iter%cur%p%vl(i)%p%info%lcontainer) &
          ntotal_vars = ntotal_vars + 1
      END DO
    END DO
    ! temporary variables needed for variable group parsing
    ALLOCATE(varlist(ntotal_vars), grp_vars(ntotal_vars), stat=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! -- loop over all output namelists
    p_onl => first_output_name_list
    DO WHILE (ASSOCIATED(p_onl))
      ! process i_typ=ml_varlist, pl_varlist, hl_varlist, il_varlist:
      DO i_typ = 1, 4

        SELECT CASE(i_typ)
        CASE (level_type_ml)
          in_varlist => p_onl%ml_varlist
        CASE (level_type_pl)
          in_varlist => p_onl%pl_varlist
        CASE (level_type_hl)
          in_varlist => p_onl%hl_varlist
        CASE (level_type_il)
          in_varlist => p_onl%il_varlist
        END SELECT

        ! Get the number of variables in varlist
        nvars = 0
        DO WHILE (nvars < SIZE(in_varlist))
          IF (in_varlist(nvars+1) == ' ') EXIT
          nvars = nvars + 1
        END DO
        CALL remove_duplicates(in_varlist, nvars)
!         write(0,*) "nvars=", nvars, "  ntotal_vars=", ntotal_vars
        IF (nvars>ntotal_vars)  CALL finish(routine, "Internal error: nvars > ntotal_vars")

        ! look for variable groups ("tiles:xyz" and "group:xyz") and replace them:
        ivar = 1
        DO WHILE (ivar <= nvars)
          vname = in_varlist(ivar)

          ! FIXME: this is probably a bug: the INDEX function also
          ! gives a match if tile_prefix is found after the first
          ! position, i.e. one would need to check for index returning
          ! exactly 1 or compare vname(1:len(tile_prefix)) ==
          ! tile_prefix
          IF (INDEX(vname, TILE_PREFIX) > 0) THEN
            ! this is a tile group identifier
            ! translate group name from GRIB2 to internal nomenclature, if necessary
            grp_name = varnames_dict%get(vname(LEN(TILE_PREFIX)+1:), &
              &                          vname(LEN(TILE_PREFIX)+1:))

            tlen = len_trim(grp_name)
            grp_name(tlen+1:tlen+2) ="_t"
            tlen = tlen + 2
          ELSE IF (INDEX(vname, GRP_PREFIX) > 0) THEN
            ! this is a group identifier
            grp_name = vname(LEN(GRP_PREFIX)+1:)
            tlen = LEN_TRIM(grp_name)
          ELSE
            ! do not perform insertion if nothing matched
            ivar = ivar + 1
            CYCLE
          END IF
          ! loop over all variables and collects the variables names
          ! corresponding to the group "grp_name"
          CALL vlr_group(grp_name, grp_vars, ngrp_vars, &
            &               loutputvars_only = .TRUE., &
            &               lremap_lonlat    = (p_onl%remap == REMAP_REGULAR_LATLON), &
            &               opt_vlevel_type  = i_typ, &
            &               opt_dom_id       = p_onl%dom)
          CALL lowcase(grp_vars(1:ngrp_vars))

          ! generate varlist where "grp_name" has been replaced;
          ! duplicates are removed
          CALL insert_group(in_varlist, nvars, vname, grp_vars(1:ngrp_vars), ninserted)

          ! status output
          IF (msg_level >= 12) THEN
            CALL message(routine, "Activating group of variables: "//TRIM(grp_name))
            DO jvar=1,ngrp_vars
              CALL message(routine, "   "//TRIM(grp_vars(jvar)))
            END DO
          END IF
          ivar = ivar + ninserted
        END DO

        ! second step: look for "subtraction" of variables groups ("-varname"):
        nsubtract_vars = 0
        DO ivar = 1, nvars
          vname = in_varlist(ivar)
          IF (vname(1:1) == "-") THEN
            nsubtract_vars = nsubtract_vars + 1
            varlist(nsubtract_vars) = ADJUSTL(vname)
            nsubtract_vars = nsubtract_vars + 1
            varlist(nsubtract_vars) = ADJUSTL(vname(2:))
          END IF
        END DO
        ! remove variables
        CALL difference(in_varlist, nvars, varlist, nsubtract_vars)
      END DO ! i_typ = 1,4
      p_onl => p_onl%next

    END DO ! p_onl

    DEALLOCATE(varlist, grp_vars, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE parse_variable_groups


  !------------------------------------------------------------------------------------------------
  !> Initialize data structures for output.
  !
  !  This routine is called after reading the namelists AND setting up
  !  the domains and variables.
  !
  SUBROUTINE init_name_list_output(sim_step_info, &
    &                              opt_lprintlist, opt_l_is_ocean)


    !> Data structure containing all necessary data for mapping an
    !  output time stamp onto a corresponding simulation step index.
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info

    LOGICAL, OPTIONAL, INTENT(in) :: opt_lprintlist
    LOGICAL, OPTIONAL, INTENT(in) :: opt_l_is_ocean

    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_name_list_output"
    INTEGER,          PARAMETER :: print_patch_id = 1

    LOGICAL                              :: l_print_list ! Flag. Enables  a list of all variables
    INTEGER                              :: nfiles, &
      &                                     jp, idom, idom_log,                           &
      &                                     grid_info_mode, ierrstat,                     &
      &                                     errno
    TYPE (t_output_name_list), POINTER   :: p_onl
    TYPE(t_par_output_event),  POINTER   :: ev
    TYPE(timedelta),           POINTER   :: mtime_output_interval, &
      &                                     mtime_td, mtime_day
    TYPE(datetime),            POINTER   :: mtime_datetime_start,              &
         &                                  mtime_datetime_end
    TYPE(timedelta) :: mtime_td1, mtime_td2, mtime_td3
    TYPE(datetime) :: mtime_date1, mtime_date2
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: lower_bound_str
    INTEGER                              :: idx, istart
    LOGICAL                              :: is_io, is_stdio
    LOGICAL                              :: is_mpi_test
    INTEGER                              :: this_i_lctype

    l_print_list = .FALSE.
    IF (PRESENT(opt_lprintlist)) l_print_list = opt_lprintlist
    is_mpi_test = my_process_is_mpi_test()

    is_io = my_process_is_io()
    is_stdio = my_process_is_stdio()

    IF (.NOT. p_test_run .AND. is_stdio .AND. &
      & (l_print_list .OR. (msg_level >= 15))) THEN

      this_i_lctype = 0
      this_i_lctype = i_lctype(print_patch_id)

      CALL print_var_list(out_varnames_dict,              &
        &                 print_patch_id,                 &
        &                 gribout_config(print_patch_id), &
        &                 this_i_lctype)
    END IF

    ! -- preliminary checks:
    !
    ! We need dtime
    IF(dtime<=0._wp) CALL finish(routine, 'dtime must be set before reading output namelists')

    ! bcast_root is not used in this case
    bcast_root = 0
! 1

    ! ---------------------------------------------------------------------------


    ! Replicate physical domain setup, only the number of domains and
    ! the logical ID is needed
    IF (use_async_name_list_io .AND. .NOT. is_mpi_test) THEN
      CALL p_bcast(n_phys_dom, bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_phys_patch(1:n_phys_dom)%logical_id, bcast_root, p_comm_work_2_io)
    END IF

    ! Set the number of output domains depending on
    ! l_output_phys_patch
    n_dom_out = MERGE(n_phys_dom, n_dom, l_output_phys_patch)

    ! allocate patch info data structure for unstructured and regular
    ! grids:
    ALLOCATE(patch_info(n_dom_out), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ALLOCATE(lonlat_info(lonlat_grids%ngrids, n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! ---------------------------------------------------------------------------

    ! Set number of global cells/edges/verts and logical patch ID
    DO jp = 1, n_dom_out
      IF(l_output_phys_patch) THEN
        patch_info(jp)%log_patch_id = p_phys_patch(jp)%logical_id
        IF (.NOT. is_io) THEN
          patch_info(jp)%p_pat_c    => p_phys_patch(jp)%comm_pat_gather_c
          patch_info(jp)%nblks_glb_c = (p_phys_patch(jp)%n_patch_cells-1)/nproma + 1
          patch_info(jp)%p_pat_e    => p_phys_patch(jp)%comm_pat_gather_e
          patch_info(jp)%nblks_glb_e = (p_phys_patch(jp)%n_patch_edges-1)/nproma + 1
          patch_info(jp)%p_pat_v    => p_phys_patch(jp)%comm_pat_gather_v
          patch_info(jp)%nblks_glb_v = (p_phys_patch(jp)%n_patch_verts-1)/nproma + 1
          patch_info(jp)%max_cell_connectivity   = p_patch(patch_info(jp)%log_patch_id)%cells%max_connectivity
          patch_info(jp)%max_vertex_connectivity = p_patch(patch_info(jp)%log_patch_id)%verts%max_connectivity
        END IF
      ELSE
        patch_info(jp)%log_patch_id = jp
        IF (.NOT. is_io) THEN
          patch_info(jp)%p_pat_c    => p_patch(jp)%comm_pat_gather_c
          patch_info(jp)%nblks_glb_c = (p_patch(jp)%n_patch_cells_g-1)/nproma + 1
          patch_info(jp)%p_pat_e    => p_patch(jp)%comm_pat_gather_e
          patch_info(jp)%nblks_glb_e = (p_patch(jp)%n_patch_edges_g-1)/nproma + 1
          patch_info(jp)%p_pat_v    => p_patch(jp)%comm_pat_gather_v
          patch_info(jp)%nblks_glb_v = (p_patch(jp)%n_patch_verts_g-1)/nproma + 1
          patch_info(jp)%max_cell_connectivity   = p_patch(patch_info(jp)%log_patch_id)%cells%max_connectivity
          patch_info(jp)%max_vertex_connectivity = p_patch(patch_info(jp)%log_patch_id)%verts%max_connectivity
        END IF
      ENDIF
    ENDDO ! jp

    ! ---------------------------------------------------------------------------
    !--- determine the mode how to collect grid information for output:

    patch_info(:)%grid_info_mode    = GRID_INFO_NONE
    lonlat_info(:,:)%grid_info_mode = GRID_INFO_NONE

    p_onl => first_output_name_list
    ! Loop over all "output_nml" namelists:
    DO WHILE (ASSOCIATED(p_onl))

      idom = p_onl%dom ! domain for which this name list should be used

      ! non-existent domains are simply ignored:
      IF (p_onl%output_grid .AND. (idom <= n_dom_out)) THEN
        grid_info_mode = GRID_INFO_BCAST
        ! For hexagons, we still copy grid info from file; for
        ! triangular grids we have a faster method without file access
        ! IF (max_cell_connectivity == 6) grid_info_mode =
        ! GRID_INFO_FILE
        IF (PRESENT(opt_l_is_ocean)) THEN
          IF (opt_l_is_ocean) grid_info_mode = GRID_INFO_BCAST
        ENDIF
        IF (p_onl%remap==REMAP_REGULAR_LATLON) THEN
          lonlat_info(p_onl%lonlat_id,patch_info(idom)%log_patch_id)%grid_info_mode = grid_info_mode
        ELSE
          patch_info(idom)%grid_info_mode = grid_info_mode
        END IF
      END IF

      p_onl => p_onl%next
    ENDDO

    ! replicate grid_info_mode on I/O PEs:
    IF (use_async_name_list_io .AND. .NOT. is_mpi_test) THEN
      ! Go over all output domains
      CALL p_bcast(patch_info(1:n_dom_out)%grid_info_mode, &
        &          bcast_root, p_comm_work_2_io)
      ! A similar process as above - for the lon-lat grids
      IF (lonlat_grids%ngrids > 0) THEN ! may cause segfaults otherwise
        CALL p_bcast(lonlat_info(1:lonlat_grids%ngrids,1:n_dom)%grid_info_mode, &
          &          bcast_root, p_comm_work_2_io)
      ENDIF
    END IF

    ! Set the number of domains in output and the patch reorder information
    CALL set_patch_info

    ! Date-time computations:
    !
    ! There are two alternative implementations for setting the output
    ! intervals, "output_bounds" and "output_start" / "output_end" /
    ! "output_interval". The former defines the output events relative
    ! to the simulation start (in seconds) and the latter define the
    ! output events by setting ISO8601-conforming date-time strings.
    !
    ! If the user has set values for "output_bounds", we compute the
    ! "output_start" / "output_end" / "output_interval" from this
    ! info:
    mtime_day => newTimedelta("P1D")
    p_onl => first_output_name_list
    DO WHILE (ASSOCIATED(p_onl))
      ! there may be multiple "output_bounds" intervals, consider all:
      DO idx=1,MAX_TIME_INTERVALS
        istart = (idx-1)*3
        IF (p_onl%output_bounds(istart+1) == -1._wp) CYCLE

        CALL mtime_timedelta_from_fseconds(p_onl%output_bounds(istart+1), &
             sim_step_info%sim_start, mtime_td1)
        CALL mtime_timedelta_from_fseconds(p_onl%output_bounds(istart+2), &
             sim_step_info%sim_start, mtime_td2)
        CALL mtime_timedelta_from_fseconds(p_onl%output_bounds(istart+3), &
             sim_step_info%sim_start, mtime_td3)
        CALL timedeltaToString(mtime_td3, p_onl%output_interval(idx))

        mtime_date1 = sim_step_info%sim_start
        mtime_date1 = mtime_date1 + mtime_td1
        CALL datetimeToString(mtime_date1, p_onl%output_start(idx))
        mtime_date2 = sim_step_info%sim_start
        mtime_date2 = mtime_date2 + mtime_td2
        CALL datetimeToString(mtime_date2, p_onl%output_end(idx))

        IF (is_stdio) THEN
          WRITE (0,*) "setting output bounds as ", TRIM(p_onl%output_start(idx)), " / ", &
            &                                      TRIM(p_onl%output_end(idx)),   " / ", &
            &                                      TRIM(p_onl%output_interval(idx))
        END IF

      END DO

      !--- consistency check: do not allow output intervals < dtime:

      ! there may be multiple "output_bounds" intervals, consider all:
      INTVL_LOOP : DO idx=1,MAX_TIME_INTERVALS

        IF (p_onl%output_start(idx) /= '') THEN

          mtime_datetime_start => newDatetime(TRIM(strip_from_modifiers(p_onl%output_start(idx))))
          mtime_datetime_end   => newDatetime(TRIM(strip_from_modifiers(p_onl%output_end(idx))))

          mtime_output_interval => newTimedelta(TRIM(p_onl%output_interval(idx)),errno=errno)
          IF (errno /= SUCCESS) CALL finish(routine,"Wrong output interval")

          mtime_td => newTimedelta("PT"//TRIM(real2string(sim_step_info%dtime, '(f20.3)'))//"S")
          CALL timedeltaToString(mtime_td, lower_bound_str)
          IF (mtime_td > mtime_day)  THEN
            CALL finish(routine, "Internal error: dtime > 1 day!")
          END IF
          IF  ((mtime_output_interval < mtime_td) .and. &
            &  (mtime_datetime_start /= mtime_datetime_end)) THEN
            CALL finish(routine, "Output interval "//TRIM(p_onl%output_interval(idx))//" < dtime !")
          END IF
          CALL deallocateTimedelta(mtime_output_interval)
          CALL deallocateTimeDelta(mtime_td)

          CALL deallocateDatetime(mtime_datetime_start)
          CALL deallocateDatetime(mtime_datetime_end)
        END IF

      END DO INTVL_LOOP
      p_onl => p_onl%next
    ENDDO
    CALL deallocateTimeDelta(mtime_day)

    ! Get the number of output files needed (by counting the domains per name list)
    p_onl => first_output_name_list
    nfiles = 0
    DO WHILE(ASSOCIATED(p_onl))
      ! non-existent domains are simply ignored:
      IF(p_onl%dom <= n_dom_out) THEN

        ! Check if name_list has variables of corresponding type,
        ! then increase file counter.
        IF (p_onl%ml_varlist(1) /= ' ') nfiles = nfiles + p_onl%stream_partitions_ml
        IF (p_onl%pl_varlist(1) /= ' ') nfiles = nfiles + p_onl%stream_partitions_pl
        IF (p_onl%hl_varlist(1) /= ' ') nfiles = nfiles + p_onl%stream_partitions_hl
        IF (p_onl%il_varlist(1) /= ' ') nfiles = nfiles + p_onl%stream_partitions_il
      END IF

      p_onl => p_onl%next

    ENDDO
    WRITE(message_text,'(a,i4)') 'Number of name list output files: ',nfiles
    CALL message(routine,message_text)

    ! Allocate output_file struct for all output files

    ALLOCATE(output_file(nfiles))

    ! Init temporal accumulation fields etc.
    p_onl => first_output_name_list
    DO WHILE(ASSOCIATED(p_onl))
      IF (p_onl%dom .LE. n_dom_out) THEN
        IF (my_process_is_work()) THEN ! p_patch is not defined on io-procs
          CALL init_statistics(p_onl, use_async_name_list_io .AND. .NOT.is_mpi_test, &
            & bcast_root, p_comm_work_2_io, p_patch(patch_info(p_onl%dom)%log_patch_id))
        ELSE
          CALL init_statistics(p_onl, use_async_name_list_io .AND. .NOT.is_mpi_test, &
            & bcast_root, p_comm_work_2_io)
        END IF
      END IF
      p_onl => p_onl%next
    END DO

    ! ---------------------------------------------------------------------------
    ! If async IO is used, replicate data (mainly the variable lists) on IO procs
! 1

    output_file(:)%cdiFileID  = CDI_UNDEFID ! i.e. not opened
    output_file(:)%cdiVlistID = CDI_UNDEFID ! i.e. not defined


    ! --------------------------------------------------------------------------------------
    ! Loop over all output namelists, set up the output_file struct for all associated files
    ! --------------------------------------------------------------------------------------

    CALL output_name_lists_to_files()

    CALL assign_output_task(output_file%io_proc_id, output_file%pe_placement)

    ! ---------------------------------------------------------------------------
    ! handle grid coordinate information
    IF (.NOT. is_mpi_test) THEN

      ! Prepare the output of grid information: For each
      ! physical/logical patch we must collect the geographical
      ! locations of cells, edges, and vertices

      ! Only needed if no async name list io is used
      IF (.NOT. use_async_name_list_io) THEN
        ! Go over all output domains
        DO idom = 1, n_dom_out
          IF (patch_info(idom)%grid_info_mode == GRID_INFO_BCAST) THEN
            ! logical domain ID
            idom_log = patch_info(idom)%log_patch_id
            CALL collect_all_grid_info(p_patch(idom_log), patch_info(idom))
          END IF
        END DO
      END IF

! 1
    END IF

    CALL create_event_data(sim_step_info)

    ! If async IO is used, initialize the memory window for communication

    CALL message(routine,'Done')

  END SUBROUTINE init_name_list_output

  SUBROUTINE output_name_lists_to_files()
    TYPE (t_output_name_list), POINTER   :: p_onl
    TYPE (t_output_file),      POINTER   :: p_of
    CHARACTER(len=vname_len), POINTER :: varlist_ptr(:)
    INTEGER, POINTER                     :: pe_placement(:)
    INTEGER :: ifile, ifile_partition, npartitions, i_typ, idom, log_patch_id
    CHARACTER(*), PARAMETER :: routine = modname//"::output_name_lists_to_files"

    p_onl => first_output_name_list
    ifile = 0
    LOOP_NML : DO WHILE (ASSOCIATED(p_onl))
      idom = p_onl%dom ! domain for which this name list should be used
      ! non-existent domains are simply ignored:
      IF(idom > n_dom_out) THEN
        p_onl => p_onl%next
        CYCLE
      END IF
      log_patch_id = patch_info(idom)%log_patch_id

      ! Loop over model/pressure/height levels
      DO i_typ = 1, 4

        ! Check if name_list has variables of corresponding type
        SELECT CASE(i_typ)
        CASE (level_type_ml)
          npartitions  =  p_onl%stream_partitions_ml
          pe_placement => p_onl%pe_placement_ml(:)
          varlist_ptr  => p_onl%ml_varlist
        CASE (level_type_pl)
          npartitions  =  p_onl%stream_partitions_pl
          pe_placement => p_onl%pe_placement_pl(:)
          varlist_ptr  => p_onl%pl_varlist
        CASE (level_type_hl)
          npartitions  =  p_onl%stream_partitions_hl
          pe_placement => p_onl%pe_placement_hl(:)
          varlist_ptr  => p_onl%hl_varlist
        CASE (level_type_il)
          npartitions  =  p_onl%stream_partitions_il
          pe_placement => p_onl%pe_placement_il(:)
          varlist_ptr  => p_onl%il_varlist
        END SELECT
        IF (varlist_ptr(1) == ' ') CYCLE

        IF (npartitions > 1) THEN
          WRITE(message_text,'(a,i4,a)') "Fork file into: ", npartitions, " concurrent parts."
          CALL message(routine, message_text)
        END IF

        ! Split one namelist into concurrent, alternating files:
        DO ifile_partition = 1,npartitions

          ifile = ifile+1
          p_of => output_file(ifile)
          p_of%ilev_type = i_typ

          ! Fill data members of "t_output_file" data structures
          p_of%filename_pref   = p_onl%output_filename
          p_of%phys_patch_id   = idom
          p_of%log_patch_id    = log_patch_id
          p_of%output_type     = p_onl%filetype
          p_of%name_list       => p_onl
          p_of%remap           = p_onl%remap
          p_of%cdiCellGridID   = CDI_UNDEFID
          p_of%cdiEdgeGridID   = CDI_UNDEFID
          p_of%cdiVertGridID   = CDI_UNDEFID
          p_of%cdiLonLatGridID = CDI_UNDEFID
          p_of%cdiZonal1DegID  = CDI_UNDEFID
          p_of%cdiVlistID      = CDI_UNDEFID

          p_of%npartitions     = npartitions
          p_of%ifile_partition = ifile_partition

          ! (optional:) explicitly specified I/O rank
          p_of%io_proc_id      = -1 ! undefined MPI rank
          p_of%pe_placement    = pe_placement(ifile_partition)

          CALL add_varlist_to_output_file(p_of, varlist_ptr)
        END DO ! ifile_partition
      ENDDO ! i_typ
      p_onl => p_onl%next
    ENDDO LOOP_NML
  END SUBROUTINE output_name_lists_to_files

  SUBROUTINE assign_output_task(io_proc_id, pe_placement)
    INTEGER, INTENT(out) :: io_proc_id(:)
    INTEGER, INTENT(in) :: pe_placement(:)
    INTEGER :: i, j, nfiles
    INTEGER :: nremaining_io_procs !< no. of non-placed I/O ranks
    LOGICAL :: is_stdio
    CHARACTER(len=MAX_CHAR_LENGTH) :: proc_list_str !< string (unoccupied I/O ranks)
    INTEGER :: remaining_io_procs(process_mpi_io_size) !< non-placed I/O ranks
    LOGICAL :: occupied_pes(process_mpi_io_size) !< explicitly placed I/O ranks
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::assign_output_task"
    ! ------------------------------------------------------
    ! Set ID of process doing I/O
    ! ------------------------------------------------------

    ! --- First, set those MPI ranks which were explicitly specified
    !     by the user:
    !
    nfiles = SIZE(io_proc_id)
    occupied_pes(:) = .FALSE.
    is_stdio = my_process_is_stdio()
    IF(use_async_name_list_io) THEN
      IF (process_mpi_io_size == 0) &
        CALL finish(routine, "Asynchronous I/O but no IO procs!")
      IF (ANY(pe_placement /= -1 .AND. &
        &     (pe_placement < 0 .OR. pe_placement > process_mpi_io_size))) &
        CALL finish(routine, "Invalid explicit placement of IO rank!")
      DO i = 1, nfiles
        ! Asynchronous I/O
        !
        ! MPI ranks "p_io_pe0 ... (p_io_pe0+process_mpi_io_size-1)" are available.

        IF (pe_placement(i) /= -1) THEN
          io_proc_id(i) = pe_placement(i)
          occupied_pes(pe_placement(i)+1) = .TRUE.
        END IF
      END DO
    ELSE
      IF (ANY(pe_placement /= -1 .AND. pe_placement /= 0)) &
        &  CALL finish(routine, "Invalid explicit placement of IO rank!")

      ! Normal I/O done by the standard I/O processor
      !
      ! Only MPI rank "process_mpi_stdio_id" is available.
      io_proc_id(1:nfiles) = -2
    END IF

    ! --- Build a list of MPI ranks that have not yet been occupied:
    !
    nremaining_io_procs = 0
    DO i=1,process_mpi_io_size
      IF (.NOT. occupied_pes(i)) THEN
        nremaining_io_procs = nremaining_io_procs + 1
        remaining_io_procs(nremaining_io_procs) = (i-1)
      END IF
    END DO
    ! status output, if some ranks were explicitly specified
    IF ((process_mpi_io_size /= nremaining_io_procs) .AND. is_stdio) THEN
      WRITE (0,'(a)') " ", "I/O : Explicit placement of I/O ranks:"
      DO i = 1, nfiles
        IF (pe_placement(i) /= -1) THEN
          WRITE (0,'(2(a,i0))') "    file #", i, " placed on rank #", io_proc_id(i)
        END IF
      END DO
      IF (nremaining_io_procs > 0) THEN
        CALL sort_and_compress_list(remaining_io_procs(1:nremaining_io_procs), proc_list_str)
        WRITE (0,'(2a)') "I/O : Remaining I/O ranks: # ", TRIM(proc_list_str)
      END IF
    END IF

    ! --- Then, set MPI ranks in a Round-Robin fashion for those
    !     namelists which had no explicitly specified "pe_placement":
    !
    ! status print-out only when some PEs were explicitly set.
    IF ((process_mpi_io_size /= nremaining_io_procs) .AND. &
      & is_stdio                        .AND. &
      & (nremaining_io_procs > 0)                    .AND. &
      & ANY(pe_placement == -1)) THEN
      WRITE (0,'(a)') " ", "I/O : Round-Robin placement of I/O ranks:"
    END IF
    IF (use_async_name_list_io) THEN
      IF (ANY(pe_placement == -1) .AND. nremaining_io_procs == 0) THEN
        CALL finish(routine, "No I/O proc left after explicit placement!")
      END IF
      j = 0
      DO i = 1, nfiles
        IF (pe_placement(i) == -1) THEN
          ! Asynchronous I/O
          j = j + 1
          io_proc_id(i) = remaining_io_procs(MOD(j-1,nremaining_io_procs) + 1)
          IF ((process_mpi_io_size /= nremaining_io_procs) .AND. is_stdio) THEN
            WRITE (0,'(2(a,i0))') "    file #", i, " placed on rank #", io_proc_id(i)
          END IF
        END IF
      END DO
    END IF
    IF ((process_mpi_io_size /= nremaining_io_procs) .AND. is_stdio) THEN
      WRITE (0,'(a)') ""
    END IF

  END SUBROUTINE assign_output_task

  SUBROUTINE create_event_data(sim_step_info)
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info

    TYPE(t_event_data_local), ALLOCATABLE :: event_list_local(:)
    !> length of local list of output events
    INTEGER :: ievent_list_local, local_i, i, nfiles, num_local_events
    INTEGER :: dom_sim_step_info_jstep0
    INTEGER :: ev_tx_comm, ev_tx_root
    LOGICAL :: is_io, is_mpi_test

    !
    !  Regular output is triggered at
    ! so-called "output event steps". The completion of an output
    ! event step is communicated via non-blocking MPI messages to the
    ! root I/O PE, which keeps track of the overall event status. This
    ! event handling is static, i.e. all event occurrences are
    ! pre-defined during the initialization.
    !
    is_io = my_process_is_io()
    is_mpi_test = my_process_is_mpi_test()
    nfiles = SIZE(output_file)

    IF (.NOT. is_io) THEN
      num_local_events = nfiles
    ELSE
      num_local_events = 0
      DO i = 1, nfiles
        num_local_events &
          = num_local_events + MERGE(1,0,output_file(i)%io_proc_id == p_pe_work)
      END DO
    END IF

    ALLOCATE(event_list_local(num_local_events))
    ievent_list_local = 0
    local_i = 0
    DO i = 1, nfiles
      IF (.NOT. is_io .OR. output_file(i)%io_proc_id == p_pe_work) THEN
        local_i = local_i + 1
        output_file(i)%out_event => add_out_event(output_file(i), i, local_i, &
          &                           sim_step_info, dom_sim_step_info_jstep0, &
          &                           event_list_local, ievent_list_local)
      ELSE
        NULLIFY(output_file(i)%out_event)
      END IF
    END DO

    IF (use_async_name_list_io) THEN
      ! The root I/O MPI rank asks all participating I/O PEs for their
      ! output event info and generates a unified output event,
      ! indicating which PE performs a write process at which step.
      ev_tx_comm = p_comm_io
      ev_tx_root = 0
    ELSE IF (pio_type == pio_type_cdipio) THEN
      ! Every work rank has a full set, the root work rank transfers it to the
      ! root I/O task
      IF (.NOT. is_mpi_test .AND. p_pe_work == 0) THEN
        IF (p_pe == p_work_pe0) &
          CALL p_send(dom_sim_step_info_jstep0, p_io_pe0, p_tag=156, &
          &           comm=p_comm_work_io)
        ev_tx_comm = p_comm_work_io
        ev_tx_root = process_work_io0
      ELSE
        ev_tx_comm = mpi_comm_null
        ev_tx_root = -1
      END IF
    ELSE
      ! work rank 0 creates the output events from its own list
      ev_tx_comm = MERGE(mpi_comm_self, mpi_comm_null, p_pe_work == 0)
      ev_tx_root = 0
    END IF
    all_events => union_of_all_events(compute_matching_sim_steps, &
      &                               generate_output_filenames, &
      &                               ev_tx_comm, ev_tx_root, &
      &                               event_list_local, ievent_list_local)

    IF (dom_sim_step_info_jstep0 > 0 .AND. ASSOCIATED(all_events)) &
      &  CALL set_event_to_simstep(all_events, dom_sim_step_info_jstep0 + 1, &
      &                            isRestart(), lrecover_open_file=.TRUE.)
    ! print a table with all output events
    IF (.NOT. is_mpi_test) THEN
      IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
        & (.NOT. use_async_name_list_io .AND. my_process_is_mpi_workroot() &
        &  )) THEN
        CALL print_output_event_table(dom_sim_step_info_jstep0)
      END IF
    END IF
  END SUBROUTINE create_event_data

  SUBROUTINE print_output_event_table(dom_sim_step_info_jstep0)
    INTEGER, INTENT(in) :: dom_sim_step_info_jstep0
    CHARACTER(len=max_filename_str_len) :: osched_fname
    CHARACTER(LEN=max_char_length)       :: comp_name

    IF (ASSOCIATED(all_events)) THEN
      CALL print_output_event(all_events)                                       ! screen output
      IF (dom_sim_step_info_jstep0 > 0) THEN
        IF ( is_coupled_run() ) THEN
          comp_name = get_my_process_name()
          WRITE (osched_fname, '(3a,i0,a)') "output_schedule_", &
            TRIM(comp_name), "_steps_", dom_sim_step_info_jstep0, "+.txt"
        ELSE
          WRITE (osched_fname, '(a,i0,a)') "output_schedule_steps_", &
            dom_sim_step_info_jstep0, "+.txt"
        ENDIF
      ELSE
        IF ( is_coupled_run() ) THEN
          comp_name = get_my_process_name()
          WRITE (osched_fname, '(3a)') "output_schedule_", &
            &                          TRIM(comp_name), ".txt"
        ELSE
          osched_fname = "output_schedule.txt"
        ENDIF
      END IF
      CALL print_output_event(all_events, &
        ! ASCII file output:
        & opt_filename=TRIM(osched_fname))
    END IF
  END SUBROUTINE print_output_event_table


  FUNCTION add_out_event(of, i, local_i, sim_step_info, &
    &                    dom_sim_step_info_jstep0, &
    &                    event_list_local, ievent_list_local) &
       RESULT(out_event)
    TYPE(t_par_output_event), POINTER :: out_event
    TYPE (t_output_file), INTENT(in) :: of
    INTEGER, INTENT(in) :: i, local_i
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    INTEGER, INTENT(out) :: dom_sim_step_info_jstep0
    TYPE(t_event_data_local), INTENT(INOUT), CONTIGUOUS :: event_list_local(:)

    !> length of local list of output events
    INTEGER, INTENT(inout) :: ievent_list_local

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::add_out_event"

    TYPE(t_sim_step_info) :: dom_sim_step_info
    TYPE(t_fname_metadata) :: fname_metadata
    TYPE(timedelta), POINTER :: mtime_interval, mtime_td
    TYPE(datetime), POINTER :: mtime_datetime
    INTEGER(c_int64_t) :: total_ms
    INTEGER :: tlen
    INTEGER :: iintvl, nintvls, ifile, opt_err
    INTEGER :: errno

    TYPE(t_key_value_store), POINTER :: restartAttributes
    LOGICAL :: include_last
    CHARACTER(LEN=max_timedelta_str_len) :: time_offset_str
    CHARACTER(LEN=max_datetime_str_len) :: output_interval(max_time_intervals)
    CHARACTER(LEN=max_datetime_str_len+1) :: output_start(max_time_intervals)
    CHARACTER(len=max_char_length) :: attname !< attribute name

    ! pack file-name meta-data into a derived type to pass them
    ! to "new_parallel_output_event":
    fname_metadata%steps_per_file             = of%name_list%steps_per_file
    fname_metadata%steps_per_file_inclfirst   = of%name_list%steps_per_file_inclfirst
    fname_metadata%file_interval              = of%name_list%file_interval
    fname_metadata%phys_patch_id              = of%phys_patch_id
    fname_metadata%ilev_type                  = of%ilev_type
    fname_metadata%filename_format            = of%name_list%filename_format
    fname_metadata%filename_pref              = of%filename_pref
    fname_metadata%npartitions                = of%npartitions
    fname_metadata%ifile_partition            = of%ifile_partition
    ! set user-specified filename extension or use the default
    ! extension:
    tlen = LEN_TRIM(of%name_list%filename_extn)
    IF (of%name_list%filename_extn(1:tlen) == "default") THEN
      fname_metadata%extn = get_file_extension(of%name_list%filetype)
    ELSE
      fname_metadata%extn = of%name_list%filename_extn(1:tlen)
    END IF

    CALL getAttributesForRestarting(restartAttributes)
    IF (restartAttributes%is_init) THEN
      ! Restart case: Get starting index of ouput from restart file
      !               (if there is such an attribute available).
      WRITE(attname,'(a,i2.2)') 'output_jfile_',i
      CALL restartAttributes%get(attname, fname_metadata%jfile_offset, opt_err=opt_err)
      fname_metadata%jfile_offset = MERGE(fname_metadata%jfile_offset, 0, opt_err .EQ. 0)
    ELSE
      fname_metadata%jfile_offset = 0
    END IF

    ! set model domain start/end time
    dom_sim_step_info = sim_step_info
    CALL getPTStringFromSeconds(NINT(start_time(of%log_patch_id),i8), time_offset_str)
    mtime_td   => newTimedelta(time_offset_str)
    dom_sim_step_info%dom_start_time = time_config%tc_exp_startdate + mtime_td
    CALL deallocateTimedelta(mtime_td)

    IF (end_time(of%log_patch_id) < DEFAULT_ENDTIME) THEN
      CALL getPTStringFromSeconds(NINT(end_time(of%log_patch_id),i8), time_offset_str)
      mtime_td   => newTimedelta(time_offset_str)
      dom_sim_step_info%dom_end_time = time_config%tc_exp_startdate + mtime_td
      CALL deallocateTimedelta(mtime_td)
    ELSE
      dom_sim_step_info%dom_end_time = dom_sim_step_info%sim_end
    END IF

    include_last    = of%name_list%include_last
    output_interval = of%name_list%output_interval
    output_start    = of%name_list%output_start

    ! Handle the case that one namelist has been split into
    ! concurrent, alternating files ("streams"):
    !
    IF (of%npartitions > 1) THEN
      ! count the number of different time intervals for this event (usually 1)
      nintvls = 0
      DO WHILE (LEN_TRIM(output_start(nintvls+1)) > 0)
        nintvls = nintvls + 1
        IF (nintvls == MAX_TIME_INTERVALS) EXIT
      END DO

      DO iintvl=1,nintvls
        mtime_interval => newTimedelta(output_interval(iintvl),errno=errno)
        IF (errno /= SUCCESS) CALL finish(routine,"Wrong output interval")
        mtime_datetime => newDatetime(strip_from_modifiers(output_start(iintvl)))
        !
        ! - The start_date gets an offset of
        !         "(ifile_partition - 1) * output_interval"
        DO ifile=1,(of%ifile_partition-1)
          mtime_datetime = mtime_datetime + mtime_interval
        END DO
        CALL datetimeToString(mtime_datetime, output_start(iintvl))
        ! - The output_interval is replaced by "
        !         "npartitions * output_interval"
        total_ms = getTotalMilliSecondsTimeDelta(mtime_interval, mtime_datetime)
        total_ms = total_ms * of%npartitions

        mtime_td => newTimedelta("PT"//TRIM(int2string(INT(total_ms/1000), '(i0)'))//"S")
        CALL timedeltaToString(mtime_td, output_interval(iintvl))
        CALL deallocateTimedelta(mtime_td)
        IF (of%ifile_partition == 1) THEN
          WRITE(message_text,'(2a)') "File stream partitioning: &
               &total output interval = ", output_interval(iintvl)
          CALL message(routine, message_text)
        END IF
        ! - The "include_last" flag is set to .FALSE.
        include_last                  = .FALSE.
        ! - The "steps_per_file" counter is set to 1
        fname_metadata%steps_per_file = 1
        ! - The "steps_per_file_inclfirst" flag is set to .FALSE.
        fname_metadata%steps_per_file_inclfirst = .FALSE.
        !
        CALL deallocateTimedelta(mtime_interval)
        CALL deallocateDatetime(mtime_datetime)
      END DO ! iintvl
    END IF

    ! ------------------------------------------------------------------------------------------
    ! --- I/O PEs communicate their event data, the other PEs create
    ! --- the event data only locally for their own event control:
    out_event => new_parallel_output_event(of%name_list%ready_file,        &
         &             output_start, of%name_list%output_end, &
         &             output_interval, include_last, dom_sim_step_info,   &
         &             fname_metadata, compute_matching_sim_steps,         &
         &             generate_output_filenames, local_i, p_comm_io,      &
         &             event_list_local, ievent_list_local)
    ! ------------------------------------------------------------------------------------------
    IF (dom_sim_step_info%jstep0 > 0) &
         &  CALL set_event_to_simstep(out_event, dom_sim_step_info%jstep0 + 1, &
         &                            isRestart(), lrecover_open_file=.TRUE.)

    dom_sim_step_info_jstep0 = dom_sim_step_info%jstep0
  END FUNCTION add_out_event

  !------------------------------------------------------------------------------------------------
  !> Create meta-data for vertical axes.
  !
  SUBROUTINE create_vertical_axes(output_file)
    TYPE(t_output_file), TARGET, INTENT(INOUT) :: output_file(:)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::create_vertical_axes"
    INTEGER :: i
    TYPE (t_output_file),      POINTER   :: p_of

    DO i = 1, SIZE(output_file)
      p_of  => output_file(i)
      ! Since commit "Use intercomm for namelist output client-/server-communication."
      ! io_proc_id is -2 for all output in the synchronous case
      IF ((p_of%io_proc_id > -1) .and. (p_of%io_proc_id /= p_pe_work)) CYCLE

      p_of%verticalAxisList = t_verticalAxisList()


      IF (my_process_is_oceanic()) THEN
        CALL setup_zaxes_oce(p_of%verticalAxisList,p_of%level_selection)

      ELSE IF (my_process_is_atmo()) THEN
        SELECT CASE(p_of%ilev_type)
        CASE (level_type_ml)
          CALL setup_ml_axes_atmo(p_of%verticalAxisList, p_of%level_selection, p_of%log_patch_id)
          IF (ANY(aes_phy_config(1:n_dom)%ljsb .OR. ANY(atm_phy_nwp_config(1:n_dom)%inwp_surface == LSS_JSBACH))) &
            &  CALL setup_zaxes_jsbach(p_of%verticalAxisList)
        CASE (level_type_pl)
          CALL setup_pl_axis_atmo(p_of%verticalAxisList, nh_pzlev_config(p_of%log_patch_id)%plevels, &
            &                     p_of%level_selection)
        CASE (level_type_hl)
          CALL setup_hl_axis_atmo(p_of%verticalAxisList, nh_pzlev_config(p_of%log_patch_id)%zlevels, &
            &                     p_of%level_selection)
        CASE (level_type_il)
          CALL setup_il_axis_atmo(p_of%verticalAxisList, nh_pzlev_config(p_of%log_patch_id)%ilevels, &
            &                     p_of%level_selection)
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT

      ELSE IF (my_process_is_jsbach()) THEN
        CALL setup_zaxes_jsbach(p_of%verticalAxisList)
      ELSE IF (my_process_is_waves()) THEN
        CALL setup_zaxes_waves(p_of%verticalAxisList, p_of%level_selection, p_of%log_patch_id)
      ENDIF

    END DO
  END SUBROUTINE create_vertical_axes

  FUNCTION nlevs_of_var(info, level_selection, var_ignore_level_selection) &
       RESULT(nlevs)
    TYPE(t_var_metadata), INTENT(IN) :: info
    TYPE(t_level_selection), POINTER, INTENT(IN) :: level_selection
    LOGICAL, OPTIONAL :: var_ignore_level_selection
    INTEGER :: nlevs, info_nlevs, jk
    LOGICAL :: var_ignore_level_selection_, is_reduction_var

    var_ignore_level_selection_ = .FALSE.
    ! zonal grids are 2-dim but with a vertical axis
    is_reduction_var = info%hgrid == grid_zonal .OR. info%hgrid == grid_lonlat
    IF (info%ndims < 3 .AND. .NOT. is_reduction_var) THEN
      ! other 2-dim. var are supposed to be horizontal only
      nlevs = 1
    ELSE
      ! handle the case that a few levels have been selected out of
      ! the total number of levels:
      info_nlevs = info%used_dimensions(MERGE(1, 2, is_reduction_var))
      IF (ASSOCIATED(level_selection)) THEN
        nlevs = 0
        ! Sometimes the user mixes level-selected variables with
        ! other fields on other z-axes (e.g. soil fields) in the
        ! output namelist. We try to catch this "wrong" user input
        ! here and handle it in the following way: if the current
        ! variable does not have one (or more) of the requested
        ! levels, then we completely ignore the level selection for
        ! this variable.
        !
        ! (... but note that we accept (nlevs+1) for an nlevs variable.)
        CHECK_LOOP : DO jk=1,MIN(level_selection%n_selected, info_nlevs)
          IF (level_selection%global_idx(jk) < 1 .OR.  &
            & level_selection%global_idx(jk) > info_nlevs+1) THEN
            var_ignore_level_selection_ = .TRUE.
            nlevs = info_nlevs
            EXIT CHECK_LOOP
          ELSE
            nlevs = nlevs &
             & + MERGE(1, 0, &
             &         level_selection%global_idx(jk) >= 1 &
             &   .AND. level_selection%global_idx(jk) <= info_nlevs)
          END IF
        END DO CHECK_LOOP
      ELSE
        nlevs = info_nlevs
      END IF
    ENDIF
    IF (PRESENT(var_ignore_level_selection)) &
      var_ignore_level_selection = var_ignore_level_selection_
  END FUNCTION nlevs_of_var

  !------------------------------------------------------------------------------------------------
  SUBROUTINE add_varlist_to_output_file(of, varlist)
    TYPE(t_output_file), INTENT(INOUT), TARGET :: of
    CHARACTER(*), INTENT(IN) :: varlist(:)
    CHARACTER(*), PARAMETER :: routine = modname//"::add_varlist_to_output_file"
    INTEGER :: iv, nv, cv, tl, ivl, key_notl, svl
    CHARACTER(:), ALLOCATABLE :: vname
    LOGICAL :: found
    TYPE(t_var), POINTER :: elem
    TYPE(t_var_desc), POINTER :: var_desc   !< variable descriptor
    TYPE(t_cf_var), POINTER :: this_cf
    TYPE(t_vl_register_iter) :: vl_iter

    ! Get the number of variables in varlist
    nv = 0
    cv = 0
    svl = SIZE(varlist)
    DO WHILE (nv .LT. svl .AND. varlist(nv+1) /= ' ')
      nv = nv + 1
      IF (.NOT.is_grid_info_var(varlist(nv))) &
        & cv = cv + 1
    ENDDO
    ! Allocate a list of variable descriptors:
    of%max_vars = cv
    of%num_vars = cv ! we know this already...
    ALLOCATE(of%var_desc(cv))
    ! Allocate array of variable descriptions
    cv = 0
    DO iv = 1, nv
      IF (is_grid_info_var(varlist(iv))) CYCLE
      cv = cv + 1
      IF (of%name_list%remap .EQ. REMAP_REGULAR_LATLON) THEN
        vname = LONLAT_PREFIX//tolower(varlist(iv))
      ELSE
        vname = tolower(varlist(iv))
      END IF
      key_notl = text_hash_c(vname)
      var_desc => of%var_desc(cv)
      found = .FALSE.
      ! Loop over all var_lists listed in vl_list to find the variable
      ! Please note that there may be several variables with different time levels,
      ! we just add unconditionally all with the name varlist(ivar).
      ! Remark: The different time levels may appear in different lists
      ! or in the same list, the code will accept both
      DO WHILE(vl_iter%next())
        IF (.NOT.vl_iter%cur%p%loutput) CYCLE
        ! patch_id in var_lists always corresponds to the LOGICAL domain
        IF(vl_iter%cur%p%patch_id /= of%log_patch_id) CYCLE
        IF(vl_iter%cur%p%vlevel_type /= of%ilev_type) CYCLE
        DO ivl = 1, vl_iter%cur%p%nvars
          elem => vl_iter%cur%p%vl(ivl)%p
          ! Do not inspect element if output is disabled
          IF (.NOT.elem%info%loutput) CYCLE
          IF (elem%info%lcontainer) CYCLE
          IF (key_notl .NE. vl_iter%cur%p%key_notl(ivl)) CYCLE
          IF (of%name_list%remap .EQ. REMAP_REGULAR_LATLON) THEN
            ! If lon-lat variable is requested, skip variable if it
            ! does not correspond to the same lon-lat grid:
            IF (elem%info%hgrid .NE. GRID_REGULAR_LONLAT .OR. &
              & of%name_list%lonlat_id .NE. elem%info%hor_interp%lonlat_id) &
              & CYCLE
          ELSE
            ! On the other hand: If no lon-lat interpolation is
            ! requested for this output file, skip all variables of
            ! this kind:
            IF (elem%info%hgrid .EQ. GRID_REGULAR_LONLAT) CYCLE
          END IF ! (remap/=REMAP_REGULAR_LATLON)
          IF (vname /= tolower(get_var_name(elem%info))) CYCLE
            ! register variable
          CALL registerOutputVariable(varlist(iv))
          ! register shortnames, which are used by mvstream and possibly by the users
          IF ("" /= elem%info%cf%short_name) CALL registerOutputVariable(elem%info%cf%short_name)
          ! get time level
          tl = get_var_timelevel(elem%info%name)
          ! Found it, add it to the variable list of output file
          IF(tl == -1) THEN
            ! Not time level dependent
            IF (found) CALL finish(routine,'Duplicate var name: '//TRIM(varlist(iv)))
            var_desc%r_ptr    => elem%r_ptr
            var_desc%s_ptr    => elem%s_ptr
            var_desc%i_ptr    => elem%i_ptr
            var_desc%info     =  elem%info
            var_desc%info_ptr => elem%info
          ELSE
            IF(found) THEN
              ! We have already the info field, make some plausibility checks:
              IF (ANY(var_desc%info%used_dimensions(:) /= elem%info%used_dimensions(:))) THEN
                CALL message(routine, "Var "//TRIM(elem%info%name))
                CALL finish(routine,'Dimension mismatch TL variable: '//TRIM(varlist(iv)))
              END IF
              ! There must not be a TL independent variable with the same name
              IF (     ASSOCIATED(var_desc%r_ptr) &
                & .OR. ASSOCIATED(var_desc%s_ptr) &
                & .OR. ASSOCIATED(var_desc%i_ptr)) &
                   CALL finish(routine, 'Duplicate var name: '//TRIM(varlist(iv)))
              ! Maybe some more members of info should be tested ...
            ELSE
              ! Variable encountered the first time, set info field ...
              var_desc%info = elem%info
              ! ... and set name without .TL# suffix
              var_desc%info%name = TRIM(get_var_name(elem%info))
            ENDIF

            IF (     ASSOCIATED(var_desc%tlev_rptr(tl)%p) &
              & .OR. ASSOCIATED(var_desc%tlev_sptr(tl)%p) &
              & .OR. ASSOCIATED(var_desc%tlev_iptr(tl)%p)) &
              CALL finish(routine, 'Duplicate time level for '//TRIM(elem%info%name))
            var_desc%tlev_rptr(tl)%p => elem%r_ptr
            var_desc%tlev_sptr(tl)%p => elem%s_ptr
            var_desc%tlev_iptr(tl)%p => elem%i_ptr
            var_desc%info_ptr        => elem%info
          ENDIF
          IF (of%name_list%remap .EQ. REMAP_REGULAR_LATLON) THEN
            var_desc%info%name = elem%info%name(LEN(LONLAT_PREFIX)+1:)
            IF (tl .NE. -1) &
              & var_desc%info%name = TRIM(get_var_name(var_desc%info))
          END IF
          found = .TRUE.
        ENDDO
      ENDDO ! i = 1, SIZE(vl_list)
      ! Check that at least one element with this name has been found
      IF (found) CYCLE
      ! error reporting ...
      IF (my_process_is_stdio()) THEN
        DO WHILE(vl_iter%next())
          WRITE(message_text,'(3a, i2)') &
               'Variable list name: ',TRIM(vl_iter%cur%p%vlname), &
               ' Patch: ',vl_iter%cur%p%patch_id
          CALL message('',message_text)
          DO ivl = 1, vl_iter%cur%p%nvars
            elem => vl_iter%cur%p%vl(ivl)%p
            IF (elem%info%post_op%lnew_cf) THEN
              this_cf => elem%info%post_op%new_cf
            ELSE
              this_cf => elem%info%cf
            END IF

            WRITE (message_text,'(a,a,l1,a,a)') &
                 &     '    ',elem%info%name,              &
                 &            elem%info%loutput, '  ',     &
                 &            this_cf%long_name
            CALL message('',message_text)
          ENDDO
        ENDDO
      ENDIF
      CALL finish(routine,'Output name list variable not found: '//TRIM(vname)//&
        &", patch "//int2string(of%log_patch_id,'(i0)'))
    ENDDO ! ivar = 1,nvars
  END SUBROUTINE add_varlist_to_output_file

  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_patch_info()
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_patch_info"
    INTEGER :: jp, jl, jg, i
    LOGICAL :: is_mpi_test, is_io

    is_mpi_test = my_process_is_mpi_test()
    is_io = my_process_is_io()

    DO jp = 1, n_dom_out

      jl = patch_info(jp)%log_patch_id

      IF (.NOT. is_io) THEN
        ! Set reorder_info on work and test PE
        CALL set_reorder_info(jp, p_patch(jl)%n_patch_cells_g, p_patch(jl)%n_patch_cells,             &
          &                   p_patch(jl)%cells%decomp_info%owner_mask, p_patch(jl)%cells%phys_id,    &
          &                   p_patch(jl)%cells%decomp_info%glb_index,                                &
          &                   patch_info(jp)%ri(icell), patch_info(jp)%grid_info(icell))

        CALL set_reorder_info(jp, p_patch(jl)%n_patch_edges_g, p_patch(jl)%n_patch_edges,             &
          &                   p_patch(jl)%edges%decomp_info%owner_mask, p_patch(jl)%edges%phys_id,    &
          &                   p_patch(jl)%edges%decomp_info%glb_index,                                &
          &                   patch_info(jp)%ri(iedge), patch_info(jp)%grid_info(iedge))

        CALL set_reorder_info(jp, p_patch(jl)%n_patch_verts_g, p_patch(jl)%n_patch_verts,             &
          &                   p_patch(jl)%verts%decomp_info%owner_mask, p_patch(jl)%verts%phys_id,    &
          &                   p_patch(jl)%verts%decomp_info%glb_index,                                &
          &                   patch_info(jp)%ri(ivert), patch_info(jp)%grid_info(ivert))
        ! Set grid_filename on work and test PE
        patch_info(jp)%grid_filename = p_patch(jl)%grid_filename
        ! Set UUID on work and test PE
        patch_info(jp)%grid_uuid = p_patch(jl)%grid_uuid
        ! Set information about numberOfGridUsed on work and test PE
        patch_info(jp)%number_of_grid_used = number_of_grid_used(jl)
        patch_info(jp)%ICON_grid_file_uri  = ICON_grid_file_uri(jl)

        patch_info(jp)%max_cell_connectivity = p_patch(jl)%cells%max_connectivity
        patch_info(jp)%max_vertex_connectivity = p_patch(jl)%verts%max_connectivity

      ENDIF
! 1

    ENDDO ! jp

    ! A similar process as above - for the lon-lat grids
    DO jl = 1,lonlat_grids%ngrids
      DO jg = 1,n_dom
        IF (.NOT. lonlat_grids%list(jl)%l_dom(jg)) CYCLE
        IF (.NOT. is_io) THEN
          ! Set reorder_info on work and test PE
          CALL set_reorder_info_lonlat(lonlat_grids%list(jl)%grid,      &
            &                          lonlat_grids%list(jl)%intp(jg),  &
            &                          lonlat_info(jl,jg))
        ENDIF
! 1
      END DO ! jg
    ENDDO ! jl
! #ifndef __NO_ICON_ATMO__
    IF (.NOT. is_mpi_test) THEN
      CALL create_rank0only_ri(zonal_ri, nlat_moc, is_io)
      CALL create_rank0only_ri(profile_ri, 1, is_io)
    END IF
  END SUBROUTINE set_patch_info

  SUBROUTINE create_rank0only_ri(ri, n_pnt, is_io)
    TYPE(t_reorder_info), INTENT(out) :: ri
    INTEGER, INTENT(in) :: n_pnt
    LOGICAL, INTENT(in) :: is_io
    INTEGER :: jl
    IF (.NOT. is_io) THEN
      ri%n_glb = n_pnt
      ri%n_own = MERGE(n_pnt, 0, p_pe_work == 0)
      ALLOCATE(ri%own_idx(ri%n_own), &
           ri%own_blk(ri%n_own), &
           ri%reorder_index_own(ri%n_own), &
           ri%pe_own(0:p_n_work-1), &
           ri%pe_off(0:p_n_work-1))
      IF (p_pe_work == 0) THEN
        ! hack ahead: note that zonal data is not blocked, we add an
        ! artificial blocking of nproma=1 in
        ! mo_name_list_output::get_ptr_to_var_data!
        DO jl = 1, n_pnt
          ri%own_idx(jl) = 1
          ri%own_blk(jl) = jl
          ri%reorder_index_own(jl) = jl
        END DO
      END IF
      ri%pe_off(0) = 0
      ri%pe_own(0) = n_pnt
      DO jl = 1, p_n_work-1
        ri%pe_off(jl) = n_pnt
        ri%pe_own(jl) = 0
      END DO
    END IF
    IF (use_async_name_list_io) THEN
      CALL transfer_reorder_info(ri, &
           &    is_io, bcast_root, p_comm_work_2_io)
      ! CALL transfer_grid_info(lonlat_info(jl,jg)%grid_info, n_pnt, lonlat_info(jl,jg)%grid_info_mode)
    END IF
  END SUBROUTINE create_rank0only_ri

  !------------------------------------------------------------------------------------------------
  !> Sets the reorder_info for cells/edges/verts
  !  ATTENTION: This routine must only be called on work and test PE (i.e. not on IO PEs)
  !             The arguments don't make sense on the IO PEs anyways
  !
  SUBROUTINE set_reorder_info(phys_patch_id, n_points_g, n_points, owner_mask, phys_id, &
                              glb_index, p_ri, grid_info)

    INTEGER, INTENT(IN) :: phys_patch_id   ! Physical patch ID
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
    INTEGER, INTENT(IN) :: phys_id(:,:)    ! phys_id for logical patch
    INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch
    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info
    TYPE(t_grid_info), INTENT(INOUT) :: grid_info
    ! local variables
    INTEGER :: i, n, il, ib
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER(i8), ALLOCATABLE :: occupation_mask(:)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_info"

    ! Set the physical patch owner mask
    ALLOCATE(phys_owner_mask(n_points))
    n = 0
    DO i = 1, n_points
      il = idx_no(i)
      ib = blk_no(i)
      phys_owner_mask(i) =       owner_mask(il,ib) &
        &                  .AND. (     .NOT. l_output_phys_patch &
        &                         .OR. phys_id(il,ib) == phys_patch_id)
    ENDDO

    CALL mask2reorder_info(p_ri, phys_owner_mask, n_points_g, glb_index, &
         p_comm_work, occupation_mask, pio_type == pio_type_cdipio)
    DEALLOCATE(phys_owner_mask)


    grid_info%n_log = n_points_g ! Total points in logical domain

    IF (patch_info(phys_patch_id)%grid_info_mode == GRID_INFO_FILE) THEN
      CALL bitmask2start_count_blks(occupation_mask, INT(n_points_g, i8), &
           grid_info%log_dom_starts, grid_info%log_dom_counts)
      ! Safety check
      n = SUM(grid_info%log_dom_counts)
      IF (n /= p_ri%n_glb) THEN
        WRITE(message_text, '(2(a,i0))')  'Reordering failed, n=', n, &
             ', p_ri%n_glb = ', p_ri%n_glb
        CALL finish(routine, message_text)
      END IF
    END IF
  END SUBROUTINE set_reorder_info

  SUBROUTINE bitmask2start_count_blks(mask, nb, starts, counts)
    INTEGER(i8), INTENT(in), CONTIGUOUS :: mask(0:)
    INTEGER(i8), INTENT(in) :: nb
    INTEGER, ALLOCATABLE, INTENT(out) :: starts(:), counts(:)

    INTEGER(i8) :: i
    INTEGER(i8), PARAMETER :: nbits_i8 = BIT_SIZE(i)
    INTEGER :: num_cblk, n
    LOGICAL :: prev_is_set, current_is_set
    INTEGER(i8) :: apos, bmask

    num_cblk = 0
    prev_is_set = .FALSE.
    bmask = 1_i8
    DO i = 0_i8, nb-1_i8
      ! extract logical from single bit
      apos = i/nbits_i8
      current_is_set = IAND(mask(apos), bmask) /= 0_i8
      num_cblk = num_cblk + MERGE(1, 0, (.NOT.prev_is_set) .AND. current_is_set)
      bmask = ISHFTC(bmask, 1_i8)
      prev_is_set = current_is_set
    ENDDO
    ALLOCATE(starts(num_cblk), counts(num_cblk))
    bmask = 1_i8
    DO i = 0_i8, nb-1_i8
      apos = i/nbits_i8
      current_is_set = IAND(mask(apos), bmask) /= 0_i8
      bmask = ISHFTC(bmask, 1_i8)
      IF (current_is_set) EXIT
    END DO
    n = 0
    DO WHILE (i < nb)
      ! at this point bit i is set because of the previous loop
      n = n + 1
      ! i is zero-based but starts are Fortran indices
      starts(n) = INT(i+1_i8)
      DO WHILE (current_is_set)
        i = i + 1_i8
        IF (i >= nb) EXIT
        apos = i/nbits_i8
        current_is_set = IAND(mask(apos), bmask) /= 0_i8
        bmask = ISHFTC(bmask, 1_i8)
      END DO
      counts(n) = INT(i+1_i8) - starts(n)
      DO WHILE (.NOT. current_is_set)
        i = i + 1_i8
        IF (i >= nb) EXIT
        apos = i/nbits_i8
        current_is_set = IAND(mask(apos), bmask) /= 0_i8
        bmask = ISHFTC(bmask, 1_i8)
      END DO
    END DO
  END SUBROUTINE bitmask2start_count_blks

  !------------------------------------------------------------------------------------------------
  !> Sets the reorder_info for lon-lat-grids
  !
  SUBROUTINE set_reorder_info_lonlat(grid, intp, patch_info_ll)
    TYPE(t_lon_lat_grid),  INTENT(IN)    :: grid
    TYPE(t_lon_lat_intp),  INTENT(IN)    :: intp
    TYPE(t_patch_info_ll), INTENT(INOUT) :: patch_info_ll      ! Result: reorder info

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_info_lonlat"
    INTEGER :: ierrstat, i, n_own, n

    ! Just for safety
    IF(my_process_is_io()) CALL finish(routine, 'Must not be called on IO PEs')

    n_own                  = intp%nthis_local_pts        ! No. of own points
    patch_info_ll%ri%n_glb = grid%lon_dim * grid%lat_dim ! Total points in lon-lat grid
    patch_info_ll%ri%n_own = n_own
    ! Set index arrays to own cells/edges/verts
    ALLOCATE(patch_info_ll%ri%own_idx(n_own),           &
      &      patch_info_ll%ri%own_blk(n_own),           &
      &      patch_info_ll%ri%pe_own(0:p_n_work-1),     &
      &      patch_info_ll%ri%pe_off(0:p_n_work-1),     &
      &      patch_info_ll%ri%reorder_index_own(n_own), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    DO i=1,n_own
      patch_info_ll%ri%own_idx(i) = MOD(i - 1, nproma) + 1
      patch_info_ll%ri%own_blk(i) =    (i - 1)/nproma  + 1
    END DO ! i

    ! Gather the number of own points for every PE into p_ri%pe_own
    CALL p_allgather(n_own, patch_info_ll%ri%pe_own, comm=p_comm_work)

    ! Get offset within result array
    n = 0
    DO i = 0, p_n_work-1
      patch_info_ll%ri%pe_off(i) = n
      n = n + patch_info_ll%ri%pe_own(i)
    ENDDO

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during I/O
    DO i=1,n_own
      patch_info_ll%ri%reorder_index_own(i) = intp%global_idx(i)
    END DO


    IF (patch_info_ll%grid_info_mode == GRID_INFO_FILE) THEN
      ! mapping between logical and physical patch is trivial for
      ! lon-lat grids:
      ALLOCATE(patch_info_ll%grid_info%log_dom_starts(1), &
        &      patch_info_ll%grid_info%log_dom_counts(1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      patch_info_ll%grid_info%log_dom_starts(1) = 1
      patch_info_ll%grid_info%log_dom_counts(1) = patch_info_ll%ri%n_glb
    END IF
  END SUBROUTINE set_reorder_info_lonlat

  !------------------------------------------------------------------------------------------------
  !> Sets up the vlist for a t_output_file structure
  !
  SUBROUTINE setup_output_vlist(of)
    TYPE(t_output_file), INTENT(INOUT), TARGET :: of
    ! local constants
    CHARACTER(LEN=*), PARAMETER       :: routine = modname//"::setup_output_vlist"
    REAL(wp),         PARAMETER       :: ZERO_TOL = 1.e-15_wp
    ! local variables
    INTEGER                           :: k, i_dom, gridtype, idate, &
      &                                  itime, iret, tlen, ll_dim1, ll_dim2
    TYPE(t_lon_lat_data), POINTER     :: lonlat
    REAL(wp), PARAMETER               :: pi_180 = ATAN(1._wp)/45._wp
    INTEGER                           :: max_cell_connectivity, max_vertex_connectivity, &
      &                                  cdiInstID
    INTEGER                           :: i, cdi_grid_ids(3), nvert, errstat, &
      &                                  cdiLonLatGridID, curvilinearGridID
    REAL(wp), ALLOCATABLE             :: p_lonlat(:)
    TYPE(t_verticalAxisList), POINTER :: it
    CHARACTER(len=128)                :: comment
    LOGICAL                           :: lrotated
    REAL(wp), ALLOCATABLE             :: rotated_pts(:,:,:)
    TYPE (t_lon_lat_grid), POINTER :: grid
    INTEGER                           :: taxisID


    gridtype = GRID_UNSTRUCTURED

    i_dom = of%phys_patch_id
    max_cell_connectivity   = patch_info(i_dom)%max_cell_connectivity
    max_vertex_connectivity = patch_info(i_dom)%max_vertex_connectivity

    !
    ! The following sections add the file global properties collected in init_name_list_output
    !
    ! 1. create cdi vlist
    !
    of%cdiVlistID = vlistCreate()
    !
    ! 2. add global attributes for netCDF
    !
    ! define output generating institute
    !
    ! inquire the Institute ID from (center/subcenter)
    !
    cdiInstID = institutInq(gribout_config(i_dom)%generatingCenter,          &
      &                     gribout_config(i_dom)%generatingSubcenter, '', '')

    IF (cdiInstID == CDI_UNDEFID) &
         &   cdiInstID = institutDef(gribout_config(i_dom)%generatingCenter, &
         &                           gribout_config(i_dom)%generatingSubcenter, &
         &                           "MPIMET",    "Max-Planck-Institute for Meteorology")

    ! define Institute
    CALL vlistDefInstitut(of%cdiVlistID,cdiInstID)

    tlen = LEN_TRIM(cf_global_info%title)
    iret = cdiDefAttTxt(of%cdiVlistID, CDI_GLOBAL, 'title',       &
         &                tlen, cf_global_info%title(1:tlen))
    tlen = LEN_TRIM(cf_global_info%institution)
    iret = cdiDefAttTxt(of%cdiVlistID, CDI_GLOBAL, 'institution', &
         &                tlen, cf_global_info%institution(1:tlen))
    tlen = LEN_TRIM(cf_global_info%source)
    iret = cdiDefAttTxt(of%cdiVlistID, CDI_GLOBAL, 'source',      &
         &                tlen, cf_global_info%source(1:tlen))
    tlen = LEN_TRIM(cf_global_info%history)
    iret = cdiDefAttTxt(of%cdiVlistID, CDI_GLOBAL, 'history',     &
         &                tlen, cf_global_info%history(1:tlen))
    tlen = LEN_TRIM(cf_global_info%references)
    iret = cdiDefAttTxt(of%cdiVlistID, CDI_GLOBAL, 'references',  &
         &                tlen, cf_global_info%references(1:tlen))
    comment = cf_global_info%comment
    tlen = LEN_TRIM(comment)
    iret = cdiDefAttTxt(of%cdiVlistID, CDI_GLOBAL, 'comment',     &
      &                 tlen, comment(1:tlen))

    ! 3. add horizontal grid descriptions

    IF(of%name_list%remap == REMAP_REGULAR_LATLON) THEN

      ! Lon/Lat Interpolation requested

      of%cdiCellGridID = CDI_UNDEFID
      of%cdiEdgeGridID = CDI_UNDEFID
      of%cdiVertGridID = CDI_UNDEFID

      lonlat => lonlat_grids%list(of%name_list%lonlat_id)
      ll_dim1 = lonlat%grid%lon_dim
      ll_dim2 = lonlat%grid%lat_dim

      lrotated = ( ABS(90._wp - lonlat%grid%north_pole(2)) > ZERO_TOL .OR.  &
      &    ABS( 0._wp - lonlat%grid%north_pole(1)) > ZERO_TOL )

      IF (.NOT. lrotated) THEN
        cdiLonLatGridID = gridCreate(GRID_LONLAT, ll_dim1*ll_dim2)
        of%cdiLonLatGridID = cdiLonLatGridID
      ELSE
        cdiLonLatGridID   = gridCreate(GRID_PROJECTION, ll_dim1*ll_dim2)
        curvilinearGridID = gridCreate(GRID_CURVILINEAR, ll_dim1*ll_dim2)
        of%cdiLonLatGridID = curvilinearGridID

        CALL gridDefParamRLL(cdiLonLatGridID, lonlat%grid%north_pole(1), &
          &                  lonlat%grid%north_pole(2), 0._c_double)
      END IF

      CALL gridDefXsize(cdiLonLatGridID, ll_dim1)
      CALL gridDefXname(cdiLonLatGridID, 'lon')
      CALL gridDefXunits(cdiLonLatGridID, 'degrees_east')

      CALL gridDefYsize(cdiLonLatGridID, ll_dim2)
      CALL gridDefYname(cdiLonLatGridID, 'lat')
      CALL gridDefYunits(cdiLonLatGridID, 'degrees_north')

      ALLOCATE(p_lonlat(ll_dim1))
      IF (lonlat%grid%reg_lon_def(2) <= threshold_delta_or_intvls) THEN
        DO k=1,ll_dim1
          p_lonlat(k) = lonlat%grid%reg_lon_def(1) + REAL(k-1,wp)*lonlat%grid%reg_lon_def(2)
        END DO
      ELSE
        DO k=1,ll_dim1
          p_lonlat(k) = (lonlat%grid%start_corner(1) + REAL(k-1,wp)*lonlat%grid%delta(1)) / pi_180
        END DO
      END IF
      CALL gridDefXvals(cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)

      ALLOCATE(p_lonlat(ll_dim2))
      IF (lonlat%grid%reg_lat_def(2) <= threshold_delta_or_intvls) THEN
        DO k=1,ll_dim2
          p_lonlat(k) = lonlat%grid%reg_lat_def(1) + REAL(k-1,wp)*lonlat%grid%reg_lat_def(2)
        END DO
      ELSE
        DO k=1,ll_dim2
          p_lonlat(k) = (lonlat%grid%start_corner(2) + REAL(k-1,wp)*lonlat%grid%delta(2)) / pi_180
        END DO
      END IF
      CALL gridDefYvals(cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)

      IF (lrotated) THEN

        grid => lonlat%grid
        ! compute some entries of lon-lat grid specification:
        CALL compute_lonlat_specs(grid)
        ALLOCATE(rotated_pts(grid%lon_dim, grid%lat_dim, 2), stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed!')
        ! compute grid points of rotated lon/lat grid
        CALL rotate_latlon_grid(grid, rotated_pts)

        CALL gridDefXsize(curvilinearGridID, ll_dim1)
        CALL gridDefYsize(curvilinearGridID, ll_dim2)
        CALL gridDefXvals(curvilinearGridID, rotated_pts(:,:,1) * rad2deg)
        CALL gridDefYvals(curvilinearGridID, rotated_pts(:,:,2) * rad2deg)
        CALL gridDefXname(curvilinearGridID, 'lon')
        CALL gridDefYname(curvilinearGridID, 'lat')
        CALL gridDefProj(curvilinearGridID, cdiLonLatGridID)

        ! clean up
        DEALLOCATE(rotated_pts, stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed!')

      END IF

! #ifndef __NO_ICON_ATMO__
    ELSE

      ! Cells
        of%cdiCellGridID = gridCreate(gridtype, patch_info(i_dom)%ri(icell)%n_glb)
        CALL gridDefNvertex(of%cdiCellGridID, max_cell_connectivity)
      !
      CALL gridDefXname(of%cdiCellGridID, 'clon')
      CALL gridDefXlongname(of%cdiCellGridID, 'center longitude')
      CALL gridDefXunits(of%cdiCellGridID, 'radian')
      !
      CALL gridDefYname(of%cdiCellGridID, 'clat')
      CALL gridDefYlongname(of%cdiCellGridID, 'center latitude')
      CALL gridDefYunits(of%cdiCellGridID, 'radian')
      !
      CALL gridDefUUID(of%cdiCellGridID, patch_info(i_dom)%grid_uuid%DATA)
      CALL gridDefReference(of%cdiCellGridID,TRIM(patch_info(i_dom)%ICON_grid_file_uri))
      !
      CALL gridDefNumber(of%cdiCellGridID, patch_info(i_dom)%number_of_grid_used)

      !
      ! not clear whether meta-info GRID_CELL or GRID_UNSTRUCTURED_CELL should be used
      CALL gridDefPosition(of%cdiCellGridID, GRID_CELL)

      ! Single point grid for monitoring
      of%cdiSingleGridID = gridCreate(GRID_LONLAT, 1)
      !
      CALL griddefxsize(of%cdiSingleGridID, 1)
      CALL griddefysize(of%cdiSingleGridID, 1)
      CALL griddefxvals(of%cdiSingleGridID, (/0.0_wp/))
      CALL griddefyvals(of%cdiSingleGridID, (/0.0_wp/))

      ! Zonal 1 degree grid
      of%cdiZonal1DegID  = gridCreate(GRID_LONLAT,nlat_moc)
      CALL griddefxsize(of%cdiZonal1DegID, 1)
      CALL griddefxvals(of%cdiZonal1DegID, (/0.0_wp/))
      CALL griddefysize(of%cdiZonal1DegID, nlat_moc)
      ALLOCATE(p_lonlat(nlat_moc))
      DO k=1,nlat_moc
        p_lonlat(k) = -90.0_wp-90._wp/REAL(nlat_moc,wp) &
          &           + REAL(k*180,wp)/REAL(nlat_moc, wp)
      END DO
      CALL griddefyvals(of%cdiZonal1DegID, p_lonlat)
      DEALLOCATE(p_lonlat)

      ! Verts
      nvert = MERGE(max_vertex_connectivity, 9-max_cell_connectivity, &
           my_process_is_oceanic())
        of%cdiVertGridID = gridCreate(gridtype, patch_info(i_dom)%ri(ivert)%n_glb)
        CALL gridDefNvertex(of%cdiVertGridID, nvert)
      !
      CALL gridDefXname(of%cdiVertGridID, 'vlon')
      CALL gridDefXlongname(of%cdiVertGridID, 'vertex longitude')
      CALL gridDefXunits(of%cdiVertGridID, 'radian')
      !
      CALL gridDefYname(of%cdiVertGridID, 'vlat')
      CALL gridDefYlongname(of%cdiVertGridID, 'vertex latitude')
      CALL gridDefYunits(of%cdiVertGridID, 'radian')
      !
      CALL gridDefUUID(of%cdiVertGridID, patch_info(i_dom)%grid_uuid%DATA)
      CALL gridDefReference(of%cdiCellGridID,TRIM(patch_info(i_dom)%ICON_grid_file_uri))
      !
      CALL gridDefNumber(of%cdiVertGridID, patch_info(i_dom)%number_of_grid_used)

      !
      ! not clear whether meta-info GRID_VERTEX or GRID_UNSTRUCTURED_VERTEX should be used
      CALL gridDefPosition(of%cdiVertGridID, GRID_VERTEX)

      ! Edges

        of%cdiEdgeGridID = gridCreate(gridtype, patch_info(i_dom)%ri(iedge)%n_glb)
        CALL gridDefNvertex(of%cdiEdgeGridID, 4)
      !
      CALL gridDefXname(of%cdiEdgeGridID, 'elon')
      CALL gridDefXlongname(of%cdiEdgeGridID, 'edge midpoint longitude')
      CALL gridDefXunits(of%cdiEdgeGridID, 'radian')
      !
      CALL gridDefYname(of%cdiEdgeGridID, 'elat')
      CALL gridDefYlongname(of%cdiEdgeGridID, 'edge midpoint latitude')
      CALL gridDefYunits(of%cdiEdgeGridID, 'radian')
      !
      CALL gridDefUUID(of%cdiEdgeGridID, patch_info(i_dom)%grid_uuid%DATA)
      CALL gridDefReference(of%cdiCellGridID,TRIM(patch_info(i_dom)%ICON_grid_file_uri))
      !
      CALL gridDefNumber(of%cdiEdgeGridID, patch_info(i_dom)%number_of_grid_used)

      !
      ! not clear whether meta-info GRID_EDGE or GRID_UNSTRUCTURED_EDGE should be used
      CALL gridDefPosition(of%cdiEdgeGridID, GRID_EDGE)

      of%cdiLonLatGridID = CDI_UNDEFID

      ! If wanted, set grid info
      IF(of%name_list%output_grid) THEN
        SELECT CASE(of%name_list%filetype)
        CASE (FILETYPE_NC2, FILETYPE_NC4)
          ! encode grid info in NetCDF format:
          SELECT CASE(patch_info(i_dom)%grid_info_mode)
          CASE (GRID_INFO_FILE)
            CALL copy_grid_info(of, patch_info(i_dom))
          CASE (GRID_INFO_BCAST)
            IF (.NOT. my_process_is_mpi_test()) THEN
              cdi_grid_ids(icell) = of%cdiCellGridID
              cdi_grid_ids(iedge) = of%cdiEdgeGridID
              cdi_grid_ids(ivert) = of%cdiVertGridID
              DO i = 1, 3 ! icell, iedge, ivert
                CALL set_grid_info_netcdf(cdi_grid_ids(i), &
                     patch_info(i_dom)%grid_info(i))
              END DO
            END IF
          END SELECT
        CASE (FILETYPE_GRB2)
          ! handled later...
        CASE DEFAULT
          CALL finish(routine, "Unknown grid type")
        END SELECT
      END IF

    ENDIF

    !
    ! 4. add vertical grid descriptions

    ! generate the CDI IDs for the vertical axes in the list
    it => of%verticalAxisList
    DO WHILE (ASSOCIATED(it))
      IF (.NOT. ASSOCIATED(it%axis)) CALL finish(routine, "Internal error!")
      CALL it%axis%cdiZaxisCreate()
      it => it%next
    END DO

    !
    ! 5. output does contain absolute time
    !
    SELECT CASE (of%name_list%mode)
    CASE (1)  ! forecast mode
     taxisID = taxisCreate(TAXIS_RELATIVE)

     IF (of%name_list%taxis_tunit > 10 .OR. of%name_list%taxis_tunit < 1 ) THEN
       of%name_list%taxis_tunit=TUNIT_MINUTE
       CALL message('','invalid taxis_tunit, reset to TUNIT_MINUTE')
     END IF
     CALL taxisDefTunit (taxisID, of%name_list%taxis_tunit)

     SELECT CASE(calendarType())
     CASE (mtime_proleptic_gregorian)
       CALL taxisDefCalendar (taxisID, dtime_proleptic_gregorian)
     CASE (mtime_year_of_360_days)
       CALL taxisDefCalendar (taxisID, dtime_cly360)
     CASE default
       CALL finish(routine, "Unsupported calendar!")
     END SELECT
     idate = cdiEncodeDate(INT(time_config%tc_exp_startdate%date%year),  &
       &                   INT(time_config%tc_exp_startdate%date%month), &
       &                   INT(time_config%tc_exp_startdate%date%day))
     itime = cdiEncodeTime(time_config%tc_exp_startdate%time%hour, time_config%tc_exp_startdate%time%minute, &
                           INT(time_config%tc_exp_startdate%time%second))

     CALL taxisDefRdate (taxisID, idate)
     CALL taxisDefRtime (taxisID, itime)

     !WRITE(6,'(a,i,a,i)')'idate ',idate,' ',taxisInqRdate(taxisID)
     !WRITE(6,'(a,i,a,i)')'itime ',itime,' ',taxisInqRtime(taxisID)
    CASE (2)  ! climate mode
     taxisID = taxisCreate(TAXIS_ABSOLUTE)
    CASE DEFAULT
     taxisID = taxisCreate(TAXIS_ABSOLUTE)
    END SELECT

    !
    CALL vlistDefTaxis(of%cdiVlistID, taxisID)

    !
    ! add variables
    !
    CALL add_variables_to_vlist(of)

    ! GRB2 format: define geographical longitude, latitude as special
    ! variables "RLON", "RLAT". Note that the grid information may be
    ! contained in the patch_info data structure (for a different
    ! output_file) but still we may do not want to include it into
    ! this output file's (of) data file.
    IF ((of%name_list%output_grid)                                      .AND. &
      & (patch_info(of%phys_patch_id)%grid_info_mode /= GRID_INFO_NONE) .AND. &
      & (of%name_list%filetype == FILETYPE_GRB2)) THEN
      CALL set_grid_info_grb2(of)
    END IF

  END SUBROUTINE setup_output_vlist




  !------------------------------------------------------------------------------------------------
  !> define variables and attributes
  !
  SUBROUTINE add_variables_to_vlist(of)
    TYPE (t_output_file), INTENT(INOUT), TARGET :: of
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::add_variables_to_vlist"
    TYPE (t_var_metadata), POINTER :: info
    INTEGER                        :: iv, vlistID, gridID, zaxisID, this_i_lctype
    TYPE(t_verticalAxis), POINTER  :: zaxis
    REAL(wp)                       :: missval
    LOGICAL                        :: is_mpi_test

    is_mpi_test = my_process_is_mpi_test()
    vlistID = of%cdiVlistID

    DO iv = 1, of%num_vars
      !
      info => of%var_desc(iv)%info
      !
      ! set grid ID
      !
      SELECT CASE (info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        info%cdiGridID = of%cdiCellGridID
      CASE(GRID_LONLAT)
        info%cdiGridID = of%cdiSingleGridID
      CASE(GRID_ZONAL)
        info%cdiGridID = of%cdiZonal1DegID
      CASE(GRID_UNSTRUCTURED_VERT)
        info%cdiGridID = of%cdiVertGridID
      CASE(GRID_UNSTRUCTURED_EDGE)
        info%cdiGridID = of%cdiEdgeGridID
      CASE (GRID_REGULAR_LONLAT)
        info%cdiGridID = of%cdiLonLatGridID
      CASE DEFAULT
        CALL finish(routine, 'GRID definition missing for '//TRIM(info%name))
      END SELECT

      gridID = info%cdiGridID

      ! set z axis ID
      !
      zaxis => of%verticalAxisList%getEntry(icon_zaxis_type=info%vgrid)
      IF (.NOT. ASSOCIATED(zaxis)) THEN
        WRITE (message_text,'(a,i0,a)') 'Zaxis no. ', info%vgrid,' undefined.'
        CALL finish(routine, message_text)
      END IF
      zaxisID = zaxis%cdi_id
      IF (info%lmask_boundary .AND. ANY(config_lmask_boundary(:)) .AND. &
        &      (info%hgrid == GRID_UNSTRUCTURED_CELL)) THEN
        missval = BOUNDARY_MISSVAL
      END IF
      IF (info%lmiss) THEN
        ! Set the missing value. Currently only real valued variables
        ! are allowed, so we can always use info%missval%rval
        IF ((.NOT.use_async_name_list_io .OR. is_mpi_test) .OR. use_dp_mpi2io) THEN
          missval = info%missval%rval
        ELSE
          ! In cases, where we use asynchronous output and the data is
          ! transferred using a single-precision buffer, we need to
          ! transfer the missing value to single-precision as well.
          ! Otherwise, in pathological cases, the missing value and
          ! the masked data in the buffer might be different values.
          missval = REAL(REAL(info%missval%rval,sp),dp)
        END IF
      END IF

      this_i_lctype = 0
      this_i_lctype = i_lctype(of%phys_patch_id)

      info%cdiVarID = create_cdi_variable(vlistID, gridID, zaxisID,         &
        &                                 info, missval, of%output_type,    &
        &                                 gribout_config(of%phys_patch_id), &
        &                                 this_i_lctype,                    &
        &                                 out_varnames_dict)

    ENDDO
    !
  END SUBROUTINE add_variables_to_vlist


  SUBROUTINE transfer_grid_info(grid_info, n_glb, grid_info_mode)
    TYPE(t_grid_info), INTENT(INOUT)    :: grid_info
    INTEGER,              INTENT(IN)    :: n_glb, grid_info_mode
    INTEGER :: nblk, temp(2)
    LOGICAL :: is_io

    IF (grid_info_mode == GRID_INFO_FILE) THEN
      is_io = my_process_is_io()
      IF (.NOT. is_io) THEN
        temp(1) = grid_info%n_log
        temp(2) = SIZE(grid_info%log_dom_starts)
      END IF
      CALL p_bcast(temp, bcast_root, p_comm_work_2_io)
      IF (is_io) THEN
        grid_info%n_log = temp(1)
        nblk = temp(2)
        ALLOCATE(grid_info%log_dom_starts(nblk))
        ALLOCATE(grid_info%log_dom_counts(nblk))
      END IF
      CALL p_bcast(grid_info%log_dom_starts, bcast_root, p_comm_work_2_io)
      CALL p_bcast(grid_info%log_dom_counts, bcast_root, p_comm_work_2_io)
    END IF

  END SUBROUTINE transfer_grid_info



  SUBROUTINE registerOutputVariable(vname)
    CHARACTER(*), INTENT(IN) :: vname

    CALL outputRegister%put(vname, 1)
  END SUBROUTINE registerOutputVariable

  LOGICAL FUNCTION isRegistered(vname)
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER :: dummy, err

    CALL outputRegister%get(vname, dummy, opt_err=err)
    isRegistered = err .EQ. 0
  END FUNCTION


! 1


  !> Utility routine: Strip date-time stamp (string) from modifiers,
  !  e.g. ">", "<".
  !
  FUNCTION strip_from_modifiers(dt_string)
    CHARACTER(LEN=*), INTENT(IN) :: dt_string
    CHARACTER(LEN=LEN_TRIM(dt_string)) :: strip_from_modifiers
    ! local variables
    CHARACTER :: char

    strip_from_modifiers = remove_whitespace(dt_string)
    char = strip_from_modifiers(1:1)
    SELECT CASE(char)
    CASE ('>','<')
      strip_from_modifiers = strip_from_modifiers(2:)
    END SELECT
  END FUNCTION strip_from_modifiers

END MODULE mo_name_list_output_init
