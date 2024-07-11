! Data structures and subroutines for meteogram output.
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
! The sampling intervals for meteogram data are independent
! from global output steps. Values are buffered in memory until
! the next field output is invoked.
! Before each write operation, data is gathered from all working
! PEs by the IO PE and written to a NetCDF file.
!
! How to add new variables for sampling:
! --------------------------------------
! In SR "meteogram_setup_variables", insert
!
! a) for volume variables:
!      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
!                        "myvarname", "myvarunit", "long name", &
!                        var_info, <state_var>)
!
!      where "VAR_GROUP_ATMO_ML" (or "VAR_GROUP_ATMO_HL" or
!      "VAR_GROUP_SOIL_HL", ...)  determines the level heights of this
!      variable.  The argument <state_var> denotes a 2D, 3D, or 4D
!      pointer to the corresponding data.
!
! b) for surface variables:
!      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SFC, &
!                       "myvarname", "myvarunit", "long name", &
!                       sfc_var_info, <state_var>)
!
! How to add additional diagnostic quantities
! -------------------------------------------
! In SR "meteogram_setup_variables", insert
!
!      CALL add_atmo/sfc_var (meteogram_config, var_list, IBSET(VAR_GROUP_XX, FLAG_DIAG), &
!         &                   "myvarname", "myvarunit", "long name", &
!         &                   var_info, <var>)
! where <var> denotes a variable of _equal_size_.
!
! The computation of the diagnostic quantities should be placed inside
!  SR compute_diagnostics()
! based on sampled values. Here, one may use the utility functions get_var/get_sfcvar for
! convencience.
!
! Roles in MPI communication:
! ---------------------------
! Depending on the namelist parameter "ldistributed" and the number
! of output PEs (i.e. asynchronous or synchronous I/O mode) the MPI
! tasks have the following functions:
!
! 1) use_async_name_list_io == .FALSE. (synchronous I/O)
!
!    1a) ldistributed == .TRUE.
!
!        All MPI tasks are writing their own files from their
!        individual out_buf.
!        Thus, all PEs have the flag "l_is_writer" enabled.
!
!    1b) ldistributed == .FALSE.
!
!        All PEs are sampling meteogram data and send it to a single
!        writing PE (all PEs have the flag "l_is_sender" enabled).
!        One of the working PEs is collecting all meteogram data in its
!        out_buf structure, opens, writes, and
!        closes the NetCDF file. The MPI rank of this PE is
!        "process_mpi_all_workroot_id", this PE has the flag
!        "l_is_collecting_pe" enabled.
! 2) use_async_name_list_io == .TRUE. => num_io_procs > 0 (asynchronous I/O without CDI-PIO)
!
!    2a) ldistributed == .TRUE.
!
!        Invalid case, caught by namelist cross checks
!
!    2b) ldistributed == .FALSE.
!
!        The last I/O PE collects data from working PEs and writes
!        the NetCDF output. Thus, this PE has "l_is_collecting_pe"
!        and "l_is_writer" enabled.  Since this output PE has no
!        information on variable (levels) and patches, it has also
!        the flag "is_pure_io_pe" enabled and receives this setup
!        from a dedicated working PE (highest worker rank with data
!        or highest worker rank if no PE has any station
!        assigned). The latter has "l_is_varlist_sender" enabled.
!
! Known limitations:
! ------------------
!
! - So far, only NetCDF file format is supported.
!
! - ASCII output data (similar to COSMO) must be generated in a
!   post-processing step.
!
! - In case of an application crash, latest meteogram data, which has
!   not yet been written to hard disk, will be lost. The recommended
!   work-around is to choose the max_time_stamps namelist
!   configuration variable to match the restart frequency of the model.
!
! TODO[FP] : use the same GNAT data structure as for the RBF
!            coefficient computation!
! TODO[FP] : use cdi functionality instead of direct NetCDF access.
! TODO[FP] : MPI communication of height levels and header info is
!            necessary only once at the beginning.

MODULE mo_meteogram_output

  USE mo_kind,                  ONLY: wp, i8
  USE mtime,                    ONLY: datetime, datetimeToPosixString,    &
    &                              MAX_DATETIME_STR_LEN,                  &
    &                              max_timedelta_str_len, getPTStringFromMS
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_mpi,                   ONLY: p_n_work, p_allreduce_max,          &
    &                                 get_my_mpi_all_id, p_wait,          &
    &                                 p_send, p_isend, p_irecv, p_pe_work,&
    &                                 p_send_packed, p_irecv_packed,      &
    &                                 p_int_byte, mpi_any_source,         &
    &                                 p_pack_int,    p_pack_real,         &
    &                                 p_pack_int_1d, p_pack_real_1d,      &
    &                                 p_pack_string, p_pack_real_2d,      &
    &                                 p_pack_size_int, p_pack_size_real_dp,&
    &                                 p_unpack_int,    p_unpack_real,     &
    &                                 p_unpack_int_1d, p_unpack_real_1d,  &
    &                                 p_unpack_string, p_unpack_real_2d,  &
    &                                 my_process_is_mpi_workroot,         &
    &                                 my_process_is_io,                   &
    &                                 my_process_is_work,                 &
    &                                 my_process_is_mpi_test,             &
    &                                 p_real_dp_byte, p_int,              &
    &                                 p_comm_work, p_comm_work_2_io,      &
    &                                 process_mpi_io_size,                &
    &                                 p_comm_remote_size
#ifndef NOMPI
  USE mpi
#endif
  USE mo_model_domain,          ONLY: t_patch, t_grid_cells
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_impl_constants,        ONLY: inwp, max_dom, SUCCESS, REAL_T
  USE mo_math_constants,        ONLY: pi, pi_180
  USE mo_communication,         ONLY: idx_1d, blk_no, idx_no
  USE mo_ext_data_types,        ONLY: t_external_data, t_external_atmos
  USE mo_nonhydro_types,        ONLY: t_nh_state, t_nh_prog, t_nh_diag, &
    &                                 t_nh_metrics
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_lnd_types,         ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, &
    &                                 t_wtr_prog

  USE mo_cf_convention,         ONLY: t_cf_var, t_cf_global, cf_global_info
  USE mo_util_string,           ONLY: int2string, one_of
  USE mo_util_uuid_types,       ONLY: t_uuid, uuid_string_length
  USE mo_util_uuid,             ONLY: uuid_unparse
  USE mo_netcdf_errhandler,     ONLY: nf
  ! TODO[FP] : When using an already built GNAT, not all of the
  ! following USEs will be necessary:
  USE mo_gnat_gridsearch,       ONLY: gnat_init_grid, gnat_destroy, t_gnat_tree, &
    &                                 gnat_query_containing_triangles,           &
    &                                 gnat_merge_distributed_queries, gk
  USE mo_dynamics_config,       ONLY: nnow
  USE mo_io_config,             ONLY: inextra_2d, inextra_3d, var_in_output, &
    &                                 celltracks_interval, gust_interval, echotop_meta
  USE mo_lnd_nwp_config,        ONLY: tile_list, ntiles_total, ntiles_water, zml_soil
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqr, iqs,               &
    &                                 iqm_max, iqni,                         &
    &                                 iqns, iqng, iqnh, iqnr, iqnc, ininact, &
                                      iqg, iqh
  USE mo_meteogram_config,      ONLY: t_meteogram_output_config, t_station_list, &
    &                                 FTYPE_NETCDF, MAX_NAME_LENGTH, MAX_NUM_STATIONS
  USE mo_name_list_output_config, ONLY: use_async_name_list_io
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config, t_atm_phy_nwp_config
  USE mo_name_list_output_types,ONLY: msg_io_meteogram_flush
  USE mo_util_phys,             ONLY: rel_hum, swdir_s
  USE mo_grid_config,           ONLY: grid_sphere_radius, is_plane_torus
  USE mo_aes_phy_memory,        ONLY: prm_field
  USE mo_fortran_tools,         ONLY: assert_acc_device_only, set_acc_host_or_device
  ! generalized meteogram output
  USE mo_var_list_register,     ONLY: t_vl_register_iter
  USE mo_var,                   ONLY: t_var_ptr, level_type_ml
  USE mo_var_metadata,          ONLY: get_var_timelevel
  USE mo_var_metadata_types,    ONLY: t_var_metadata
  USE mo_var_groups,            ONLY: var_groups_dyn
  USE mo_util_string,           ONLY: toupper
  USE mo_zaxis_type,            ONLY: ZA_SURFACE, ZA_ATMOSPHERE, &
    &                                 ZA_REFERENCE, ZA_REFERENCE_HALF
#ifdef _OPENACC
  USE openacc,                  ONLY: acc_is_present
#endif

  USE mo_netcdf

  IMPLICIT NONE

  PRIVATE
  INTEGER,                        PARAMETER :: dbg_level     = 0
  CHARACTER(LEN=*), PARAMETER :: modname       = 'mo_meteogram_output'

  ! IO routines.
  ! called collectively, though non-IO PEs are occupied
  ! only for the case of distributed write mode.
  PUBLIC ::  meteogram_init
  PUBLIC ::  meteogram_is_sample_step
  PUBLIC ::  meteogram_sample_vars
  PUBLIC ::  meteogram_finalize
  PUBLIC ::  meteogram_flush_file

  INTEGER, PARAMETER :: MAX_NVARS            =  300  !< max. number of sampled 3d vars
  INTEGER, PARAMETER :: MAX_NSFCVARS         =  300  !< max. number of sampled surface vars
  INTEGER, PARAMETER :: MAX_DESCR_LENGTH     =  128  !< length of info strings (see cf_convention)
  INTEGER, PARAMETER :: MAX_DATE_LEN         =   16  !< length of iso8601 date strings
  ! arbitrarily chosen value for buffer size (somewhat large for safety reasons)
  INTEGER, PARAMETER :: mtgrm_pack_header_ints =  2  !< * p_int_byte

  INTEGER, PARAMETER :: TAG_VARLIST          =   99  !< MPI tag for communication of variable info
  INTEGER, PARAMETER :: TAG_MTGRM_MSG        = 77777 !< MPI tag (base) for communication of meteogram data
  !> separating of tag IDs for different domains
  INTEGER, PARAMETER :: tag_domain_shift     = max_num_stations
  ! flags for communication of variable info
  INTEGER, PARAMETER :: FLAG_VARLIST_ATMO    =    0
  INTEGER, PARAMETER :: FLAG_VARLIST_SFC     =    1
  INTEGER, PARAMETER :: FLAG_VARLIST_END     =   -1

  !! Groups of variables; using this we can distinguish different height axes
  INTEGER, PARAMETER :: VAR_GROUP_ATMO_ML    =    1  !< variables defined on model levels
  INTEGER, PARAMETER :: VAR_GROUP_ATMO_HL    =    2  !< variables defined on half levels
  INTEGER, PARAMETER :: VAR_GROUP_SURFACE    =    3  !< surface variables
  INTEGER, PARAMETER :: VAR_GROUP_SOIL_ML    =    4  !< variables defined on soil half levels
  INTEGER, PARAMETER :: VAR_GROUP_SOIL_MLp2  =    5  !< height levels [0m, soil half levels, -14.58m]
  INTEGER, PARAMETER :: FLAG_DIAG            =    4  !< Flag bit: if set then this variable is a diagnostic

  !>
  !! Generic interface for adding atmospheric vars to list (required
  !! to cope with 4d vars, e.g. with tiles):
  INTERFACE add_atmo_var
    MODULE PROCEDURE add_atmo_var_3d
    MODULE PROCEDURE add_atmo_var_4d
  END INTERFACE

  !>
  !! Generic interface for adding surface vars to list (required
  !! to cope with 3d vars, e.g. with tiles):
  INTERFACE add_sfc_var
    MODULE PROCEDURE add_sfc_var_2d
    MODULE PROCEDURE add_sfc_var_3d
  END INTERFACE

  !>
  !! Storage for information on a single variable
  !!
  TYPE t_var_info
    TYPE(t_cf_var)        :: cf              !< variable name, unit
    INTEGER               :: igroup_id       !< variable group (surface vars, soil temperatures, ...)
    INTEGER               :: nlevs           !< number of levels for this variable
    REAL(wp), POINTER     :: p_source(:,:,:) !< pointer to source array  (nproma, nlev, nblk)
  END TYPE t_var_info

  !>
  !! Storage for information on a single surface variable
  !!
  TYPE t_sfc_var_info
    TYPE(t_cf_var)        :: cf              !< variable name, unit
    INTEGER               :: igroup_id       !< variable group (surface vars, soil temperatures, ...)
    REAL(wp), POINTER     :: p_source(:,:)   !< pointer to source array
  END TYPE t_sfc_var_info

  !> number of time- and variable-invariant items per station
  INTEGER, PARAMETER :: num_time_inv = 6

  TYPE t_a_2d
    REAL(wp), ALLOCATABLE :: a(:,:)
  END TYPE t_a_2d

  TYPE t_a_3d
    REAL(wp), ALLOCATABLE :: a(:,:,:)
  END TYPE t_a_3d

  TYPE t_mtgrm_out_buffer
    !> each element contains data for (nstations,nlevs,ntimestamps) of variable
    TYPE(t_a_3d), ALLOCATABLE :: atmo_vars(:)
    !> each element contains data for (nstations,ntimestamps) of variable
    TYPE(t_a_2d), ALLOCATABLE :: sfc_vars(:)
    !> global index of station specification
    INTEGER, ALLOCATABLE :: station_idx(:)
  END TYPE t_mtgrm_out_buffer

  TYPE t_mtgrm_invariants
    !> height buffer for each variable
    TYPE(t_a_2d), ALLOCATABLE :: heights(:)
    ! Tile info
    !> tile fractions for each station and tile
    REAL(wp), ALLOCATABLE :: tile_frac(:,:)
    !> tile specific landuse classes
    INTEGER, ALLOCATABLE :: tile_luclass(:,:)
    !> soil type
    INTEGER, ALLOCATABLE :: soiltype(:)
    !> surface height
    REAL(wp), ALLOCATABLE :: hsurf(:)
    !> fraction of land
    REAL(wp), ALLOCATABLE :: frland(:)
    !> Coriolis parameter
    REAL(wp), ALLOCATABLE :: fc(:)
    !> triangle index (global idx,block)
    INTEGER, ALLOCATABLE :: tri_idx(:,:)
    !> notal number of tiles (ntiles_total+ntiles_water)
    !! if NWP tiles are set up, 1 otherwise
    INTEGER :: ntiles_mtgrm
    !> maximum no. of levels for variables
    INTEGER :: max_nlevs
  END TYPE t_mtgrm_invariants

  !>
  !! Data structure specifying NetCDF IDs
  !!
  TYPE t_ncid
    INTEGER  :: nstations, nvars, ntiles, charid, station_name, station_lat, station_lon, &
      &         station_idx, station_blk, station_hsurf, station_frland, station_fc,      &
      &         station_soiltype, station_tile_frac, station_tile_luclass,                &
      &         nsfcvars, var_name, var_unit, sfcvar_name, sfcvar_unit,                   &
      &         var_group_id, sfcvar_group_id, var_nlevs, max_nlevs, timeid,  &
      &         time_step, dateid, var_values, sfcvar_values, var_heights, var_longname,  &
      &         sfcvar_longname
    INTEGER                          :: file_id       !< meteogram file ID
  END TYPE t_ncid

  !>
  !! Data structure specifying output file for meteogram data.
  !!
  TYPE t_meteogram_file
    INTEGER                          :: ftype         !< file type (NetCDF, ...)
    LOGICAL                          :: ldistributed  !< Flag. Separate files for each PE
    CHARACTER(len=MAX_NAME_LENGTH)   :: zname         !< file name string
    CHARACTER(len=uuid_string_length):: uuid_string   !< unparsed grid UUID
    INTEGER                          :: number_of_grid_used  !< as it says
    TYPE(t_cf_global)                :: cf            !< meta info
    !> NetCDF dimension IDs
    TYPE(t_ncid)                     :: ncid
  END TYPE t_meteogram_file

  !>
  !! Data structure containing internal indices for variables
  !!
  TYPE t_var
    INTEGER :: no_atmo_vars       !< number of atmo variables declared so far
    INTEGER :: no_sfc_vars        !< number of surface variables declared so far
  END TYPE t_var

  ! -------------------------------------------------------------------------------------------

  !> Holds indices into meteogram_(local|global)_data%station%(sfc_)var
  TYPE meteogram_diag_var_indices
    ! several variable indices, stored for convenience (when computing additional diagnostics)
    INTEGER                 :: i_T        = -1,  &
      &                        i_REL_HUM  = -1,  &
      &                        i_QV       = -1,  &
      &                        i_PEXNER   = -1,  &
      &                        i_SWDIR_S  = -1,  &
      &                        i_ALB      = -1,  &
      &                        i_SWDIFD_S = -1,  &
      &                        i_SOBS     = -1
  END TYPE meteogram_diag_var_indices

  !>
  !! Data structure containing meteogram buffers and other data.
  !!
  TYPE t_buffer_state
    !> triangle index (nstations,2)
    INTEGER, ALLOCATABLE    :: tri_idx_local(:,:)
    !> meteogram file handle etc.
    TYPE(t_meteogram_file)  :: meteogram_file_info
    !> info on sample times (1:icurrent)
    !! iteration step of model
    INTEGER, ALLOCATABLE :: istep(:)
    !> date and time of point sample (iso8601)
    CHARACTER(len=MAX_DATE_LEN), ALLOCATABLE :: zdate(:)

    TYPE(t_mtgrm_out_buffer) :: out_buf

    !! -- data for distributed meteogram sampling (MPI) --
    INTEGER                 :: max_buf_size                     !< max buffer size for MPI messages

    ! time stamp info:
    !> current time stamp index
    INTEGER :: icurrent
    !> maximum number of time stamps stored before flush
    INTEGER :: max_time_stamps
    !> flush silently when time stamp buffer is exhausted or warn user?
    LOGICAL :: silent_flush

    ! variable info:
    !> info for each variable (1:nvars)
    TYPE(t_var_info), ALLOCATABLE :: var_info(:)
    !> info for each surface variable
    TYPE(t_sfc_var_info), ALLOCATABLE :: sfc_var_info(:)

    ! different roles in communication:
    LOGICAL                 :: l_is_sender, l_is_writer,         &
      &                        l_is_collecting_pe
    !> communicator to use in collection of data
    INTEGER                 :: io_collect_comm
    !> rank of PE which gathers data
    INTEGER                 :: io_collector_rank
    !> rank of PE which sends invariants (either identical for all
    !! stations or all time stamps)
    INTEGER                 :: io_invar_send_rank
    !> global rank of each station
    INTEGER                 :: global_idx(MAX_NUM_STATIONS)

    TYPE(meteogram_diag_var_indices) :: diag_var_indices
    !> number of stations for which data is available on this MPI rank
    INTEGER :: nstations_local
    !> time dimension offset to use when writing
    INTEGER :: time_offset
    !> "owner" PE for each station
    INTEGER :: pstation(MAX_NUM_STATIONS)
  END TYPE t_buffer_state

  TYPE mtgrm_pack_buf
    !> MPI buffer for variable info
    CHARACTER, ALLOCATABLE  :: msg_varlist(:)
    !> current position when buffering
    INTEGER :: pos
    LOGICAL :: l_is_varlist_sender
  END TYPE mtgrm_pack_buf

  !! -- module data: --
  TYPE(t_buffer_state), SAVE :: mtgrm(1:max_dom)


CONTAINS

  !! Set up list of variables for sampling.
  !!
  SUBROUTINE meteogram_setup_variables(meteogram_config, jg, iforcing, &
    &                                  ext_data, prog, diag, &
    &                                  prm_diag, lnd_prog, lnd_diag, prog_wtr, &
    &                                  prm_nwp_tend, &
    &                                  atm_phy_nwp_config, var_info, &
    &                                  sfc_var_info, diag_var_indices, &
    &                                  var_list, pack_buf)
    ! station data from namelist
    TYPE(t_meteogram_output_config),     INTENT(IN) :: meteogram_config
    !> patch index
    INTEGER,                             INTENT(IN) :: jg
    INTEGER,                             INTENT(IN) :: iforcing
    ! atmosphere external data
    TYPE(t_external_data),               INTENT(IN) :: ext_data
    ! physical model state and other auxiliary variables
    TYPE(t_nwp_phy_diag),                INTENT(IN) :: prm_diag
    !> model state for the NWP land physics, prognostic variables
    TYPE(t_nh_prog),                     INTENT(in) :: prog
    !> model state for the NWP land physics, diagnostic variables
    TYPE(t_nh_diag),                     INTENT(in) :: diag
    TYPE(t_lnd_prog),                   INTENT(IN) :: lnd_prog
    TYPE(t_lnd_diag),                   INTENT(IN) :: lnd_diag
    TYPE(t_wtr_prog),                   INTENT(IN) :: prog_wtr
    ! model state of physics tendencies
    TYPE(t_nwp_phy_tend),                INTENT(IN) :: prm_nwp_tend
    TYPE(t_atm_phy_nwp_config), INTENT(in) :: atm_phy_nwp_config
    !> list of atmospheric variables to setup
    TYPE(t_var_info), INTENT(inout) :: var_info(:)
    !> list of surface variables to setup
    TYPE(t_sfc_var_info), INTENT(inout) :: sfc_var_info(:)
    !> and corresponding packed description for remote receivers
    TYPE(mtgrm_pack_buf), INTENT(inout) :: pack_buf
    !> sizes of (surface) variable lists
    TYPE(t_var), INTENT(out) :: var_list
    !> indices at which to find variables for compute_diagnostics
    TYPE(meteogram_diag_var_indices), INTENT(out) :: diag_var_indices

    CHARACTER(len=max_timedelta_str_len) :: c_time_int
    CHARACTER(len=14) :: c_thresh_int
    INTEGER           :: i

    ! gereralized meteogram output
    INTEGER                       :: idx_group_mtgrm
    TYPE(t_vl_register_iter)      :: vl_iter
    TYPE(t_var_ptr), POINTER      :: elem
    TYPE(t_var_metadata), POINTER :: info
    INTEGER                       :: iv
    CHARACTER(len=128)            :: var_name_mtgrm

    var_list%no_atmo_vars = 0
    var_list%no_sfc_vars = 0

    ! -- atmosphere

    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "P", "Pa", "Pressure", &
      &               var_info, diag%pres(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "T", "K", "Temperature", &
      &               var_info, diag%temp(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "PEXNER", "-", "Exner pressure", &
      &               var_info, prog%exner(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "RHO", "kg/m^3", "Density", &
      &               var_info, prog%rho(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "THETAV", "K", "virtual potential temperature", &
      &               var_info, prog%theta_v(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "U", "m/s", "zonal wind", &
      &               var_info, diag%u(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "V", "m/s", "meridional wind", &
      &               var_info, diag%v(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
      &               "W", "m/s", "orthogonal vertical wind", &
      &               var_info, prog%w(:,:,:))
    ! add some output for turbulence diagnostic
    IF ( iforcing == inwp ) THEN
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
        &               "TKE", "m^2/s^2", "turbulent kinetic energy", &
        &               var_info, prog%tke(:,:,:))
#ifndef __NO_ICON_LES__
      IF ( .NOT. atm_phy_nwp_config%is_les_phy ) THEN
#endif
        CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
          &               "ddt_tke_hsh", "m^2/s^3", &
          &               "TKE tendency horizonzal shear production", &
          &               var_info, prm_nwp_tend%ddt_tke_hsh)
        CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
          &               "ddt_tke_pconv", "m^2/s^3", &
          &               "TKE tendency due to subgrid-scale convection", &
          &               var_info, prm_nwp_tend%ddt_tke_pconv)
#ifndef __NO_ICON_LES__
      END IF
#endif
    ENDIF ! iforcing == inwp
    ! For dry test cases: do not sample variables defined below this line:
    ! (but allow for TORUS moist runs; see call in mo_atmo_nonhydrostatic.F90)
    !IF (ltestcase .AND. les_config(jg)%is_dry_cbl) RETURN

    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "QV", "kg kg-1", "specific humidity", &
      &               var_info, prog%tracer_ptr(iqv)%p_3d(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "QC", "kg kg-1", "specific cloud water content", &
      &               var_info, prog%tracer_ptr(iqc)%p_3d(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "QI", "kg kg-1", "specific cloud ice content", &
      &               var_info, prog%tracer_ptr(iqi)%p_3d(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "QR", "kg kg-1", "rain_mixing_ratio", &
      &               var_info, prog%tracer_ptr(iqr)%p_3d(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
      &               "QS", "kg kg-1", "snow_mixing_ratio", &
      &               var_info, prog%tracer_ptr(iqs)%p_3d(:,:,:))
    CALL add_atmo_var(meteogram_config, var_list, &
      &               IBSET(VAR_GROUP_ATMO_ML, FLAG_DIAG), &
      &               "REL_HUM", "%", "relative humidity", &
      &               var_info, prog%tracer_ptr(iqv)%p_3d(:,:,:))
    IF(atm_phy_nwp_config%lhave_graupel) THEN
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QG", "kg kg-1", "graupel_mixing_ratio", &
        &               var_info, prog%tracer_ptr(iqg)%p_3d(:,:,:))
    END IF
    If (atm_phy_nwp_config%l2moment) THEN
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QH", "kg kg-1", "hail_mixing_ratio", &
        &               var_info, prog%tracer_ptr(iqh)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QNI", "kg-1", "number concentration ice", &
        &               var_info, prog%tracer_ptr(iqni)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QNS", "kg-1", "number concentration snow", &
        &               var_info, prog%tracer_ptr(iqns)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QNR", "kg-1", "number concentration rain droplet", &
        &               var_info, prog%tracer_ptr(iqnr)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QNG", "kg-1", "number concentration graupel", &
        &               var_info, prog%tracer_ptr(iqng)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QNH", "kg-1", "number concentration hail", &
        &               var_info, prog%tracer_ptr(iqnh)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QNC", "kg-1", "number concentration cloud water", &
        &               var_info, prog%tracer_ptr(iqnc)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "NIACT", "kg-1", &
        &               "number concentration activated ice nuclei", &
        &               var_info, prog%tracer_ptr(ininact)%p_3d(:,:,:))
    END IF

    IF ( iforcing == inwp ) THEN
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QV_DIA", "kg kg-1", &
        &               "total specific humidity (diagnostic)", &
        &               var_info, prm_diag%tot_ptr(iqv)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QC_DIA", "kg kg-1", &
        &               "total specific cloud water content (diagnostic)", &
        &               var_info, prm_diag%tot_ptr(iqc)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "QI_DIA", "kg kg-1", &
        &               "total specific cloud ice content (diagnostic)", &
        &               var_info, prm_diag%tot_ptr(iqi)%p_3d(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "CLC", "-", "cloud cover", &
        &               var_info, prm_diag%clc(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
        &               "TKVM", "m**2/s", &
        &               "turbulent diffusion coefficients for momentum", &
        &               var_info, prm_diag%tkvm(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
        &               "TKVH", "m**2/s", &
        &               "turbulent diffusion coefficients for heat", &
        &               var_info, prm_diag%tkvh(:,:,:))
    ENDIF ! iforcing == inwp

    CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
      &               "PHALF", "Pa", "Pressure on the half levels", &
      &               var_info, diag%pres_ifc(:,:,:))


    ! generalized meteogram output
    idx_group_mtgrm = var_groups_dyn%group_id("METEOGRAM")
    DO WHILE(vl_iter%next())
      ! skip e.g. p_tracer_list defined for advection
      IF (.NOT. vl_iter%cur%p%loutput) CYCLE
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE
      IF (vl_iter%cur%p%patch_id /= jg)  CYCLE
      LOOPVAR : DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)
        info => elem%p%info
        ! for time-level dependent variables: output only once
        IF (get_var_timelevel(info%name) > 1)  CYCLE
        IF ( info%in_group(idx_group_mtgrm) ) THEN
          IF (ASSOCIATED(elem%p%r_ptr))  THEN
            SELECT CASE(info%data_type)
            CASE(REAL_T)
              IF (info%grib2%category /= 18) THEN
                ! change case of variable name to upper
                var_name_mtgrm = toupper(info%cf%standard_name)
              ELSE
                ! for nuclear variables (parameterCategory = 18) do not change case
                var_name_mtgrm = info%cf%standard_name
              END IF
              IF ( info%vgrid == ZA_REFERENCE .AND. info%ndims == 3 ) THEN
                CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
                  &               var_name_mtgrm, info%cf%units, &
                  &               info%cf%long_name, var_info, elem%p%r_ptr(:,:,:,1,1))
              ELSE IF ( info%vgrid == ZA_REFERENCE_HALF .AND. info%ndims == 3 ) THEN
                CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_HL, &
                  &               var_name_mtgrm, info%cf%units, &
                  &               info%cf%long_name, var_info, elem%p%r_ptr(:,:,:,1,1))
              ELSE IF ( ANY((/ZA_SURFACE, ZA_ATMOSPHERE/) == info%vgrid) .AND. info%ndims == 2 ) THEN
                CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
                  &              var_name_mtgrm, info%cf%units, &
                  &              info%cf%long_name, sfc_var_info, elem%p%r_ptr(:,:,1,1,1))
              END IF
            END SELECT
          END IF
        END IF
      ENDDO LOOPVAR
    ENDDO

    ! -- soil related
    IF (  atm_phy_nwp_config%inwp_surface == 1 ) THEN
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_SOIL_MLp2, &
        &               "T_SO", "K", "soil temperature", &
        &               var_info, lnd_diag%t_so(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_SOIL_ML, &
        &               "W_SO", "m H2O", &
        &               "total water content (ice + liquid water)", &
        &               var_info, lnd_diag%w_so(:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_SOIL_ML, &
        &               "W_SO_ICE", "m H2O", "ice content", &
        &               var_info, lnd_diag%w_so_ice(:,:,:))

      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "PL_COV", "-", "ground fraction covered by plants", &
        &              sfc_var_info, ext_data%atm%plcov(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "LA_IND", "-", "leaf area index (vegetation period)", &
        &              sfc_var_info, ext_data%atm%lai(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RO_DEPT", "m", "root depth", &
        &              sfc_var_info, ext_data%atm%rootdp(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "Z0", "m", "roughness length*g", &
        &              sfc_var_info, prm_diag%gz0(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "QV_S", "kg/kg", "specific humidity at the surface", &
        &              sfc_var_info, lnd_diag%qv_s(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "W_I", "m H2O", "water content of interception water", &
        &              sfc_var_info, lnd_diag%w_i(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "W_SNOW", "m H2O", "water content of snow", &
        &              sfc_var_info, lnd_diag%w_snow(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RUNOFF_S", "kg/m2", &
        &              "surface water runoff; sum over forecast", &
        &              sfc_var_info, lnd_diag%runoff_s(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RUNOFF_G", "kg/m2", &
        &              "soil water runoff; sum over forecast", &
        &              sfc_var_info, lnd_diag%runoff_g(:,:))
      IF (var_in_output(jg)%res_soilwatb) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "RESID_WSO", "kg/m2", &
          &              "residuum of the soil water budget", &
          &              sfc_var_info, lnd_diag%resid_wso(:,:))
      ENDIF
      IF (var_in_output(jg)%snow_melt) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "SNOW_MELT", "kg m-2", &
          &              "snow melt amount", &
          &              sfc_var_info, lnd_diag%snow_melt(:,:))
      ENDIF
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "T_SNOW", "K", "temperature of the snow-surface", &
        &              sfc_var_info, lnd_diag%t_snow(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "T_S", "K", "temperature of the ground surface", &
        &              sfc_var_info, lnd_diag%t_s(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "T_G", "K", "weighted surface temperature", &
        &              sfc_var_info, lnd_prog%t_g(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "FRESHSNW", "-", &
        &              "indicator for age of snow in top of snow layer", &
        &              sfc_var_info, lnd_diag%freshsnow(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RHO_SNOW", "kg/m**3", "snow density", &
        &              sfc_var_info, lnd_diag%rho_snow(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "H_SNOW", "m", "snow height", &
        &              sfc_var_info, lnd_diag%h_snow(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "FR_SEAICE", "-", "fraction of sea ice", &
        &              sfc_var_info, lnd_diag%fr_seaice(:,:))
    ENDIF

    ! -- single level variables

    CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
      &              "P_SFC", "Pa", "surface pressure", &
      &              sfc_var_info, diag%pres_sfc(:,:))

    IF ( iforcing == inwp ) THEN
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TCM", "-", &
        &              "turbulent transfer coefficients for momentum", &
        &              sfc_var_info, prm_diag%tcm(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TCH", "-", "turbulent transfer coefficients for heat", &
        &              sfc_var_info, prm_diag%tch(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SHFL", "W/m2", "sensible heat flux (surface)", &
        &              sfc_var_info, prm_diag%shfl_s(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "LHFL", "W/m2", "latent heat flux (surface)", &
        &              sfc_var_info, prm_diag%lhfl_s(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "QHFL_S", "kg/m2", "evapotranspiration flux (surface)", &
        &              sfc_var_info, prm_diag%qhfl_s(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "VIO3", "Pa O3", "vertically integrated ozone amount", &
        &              sfc_var_info, prm_diag%vio3(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "HMO3", "Pa", "height of O3 maximum", &
        &              sfc_var_info, prm_diag%hmo3(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "T2M", "K", "temperature in 2m", &
        &              sfc_var_info, prm_diag%t_2m(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TD2M", "K", "dew-point temperature in 2m", &
        &              sfc_var_info, prm_diag%td_2m(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "U10M", "m/s", "zonal wind in 10m", &
        &              sfc_var_info, prm_diag%u_10m(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "V10M", "m/s", "meridional wind in 10m", &
        &              sfc_var_info, prm_diag%v_10m(:,:))
      CALL getPTStringFromMS(NINT(1000_wp*gust_interval(jg), i8), c_time_int)
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "VBMAX10M", "m/s", "gust in 10m since end of previous full "//&
        &              TRIM(ADJUSTL(c_time_int(3:)))//" since model start", &
        &              sfc_var_info, prm_diag%gust10(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "dyn_gust", "m/s", "dynamical gust", &
        &              sfc_var_info, prm_diag%dyn_gust(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "con_gust", "m/s", "convective gust", &
        &              sfc_var_info, prm_diag%con_gust(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "cape_ml", "J/kg", "cape of mean surface layer parcel", &
        &              sfc_var_info, prm_diag%cape_ml(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SOBT", "W m-2", "shortwave net flux at toa", &
        &              sfc_var_info, prm_diag%swflxtoa(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "THBT", "W m-2", "longwave net flux at toa", &
        &              sfc_var_info, prm_diag%lwflxall(:,1,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SOBS", "W m-2", "shortwave net flux at surface", &
        &              sfc_var_info, prm_diag%swflxsfc(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "THBS", "W m-2", "longwave net flux at surface", &
        &              sfc_var_info, prm_diag%lwflxsfc(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "ALB", "-", "surface shortwave albedo, diffuse", &
        &              sfc_var_info, prm_diag%albdif(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RAIN_GSP", "kg/m2", &
        &              "accumulated grid-scale surface rain", &
        &              sfc_var_info, prm_diag%rain_gsp(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SNOW_GSP", "kg/m2", &
        &              "accumulated grid-scale surface snow", &
        &              sfc_var_info, prm_diag%snow_gsp(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RAIN_CON", "kg/m2", &
        &              "accumulated convective surface rain", &
        &              sfc_var_info, prm_diag%rain_con(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SNOW_CON", "kg/m2", &
        &              "accumulated convective surface snow", &
        &              sfc_var_info, prm_diag%snow_con(:,:))
      IF(atm_phy_nwp_config%lhave_graupel) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "GRAUPEL_GSP", "kg/m2", &
        &              "accumulated grid-scale surface graupel", &
        &              sfc_var_info, prm_diag%graupel_gsp(:,:))
      ENDIF
      IF(atm_phy_nwp_config%l2moment) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "HAIL_GSP", "kg/m2", &
        &              "accumulated grid-scale surface hail", &
        &              sfc_var_info, prm_diag%hail_gsp(:,:))
      ENDIF
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "PREC_CON", "kg/m2", &
        &              "accumulated convective surface precipitation", &
        &              sfc_var_info, prm_diag%prec_con(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "PREC_GSP", "kg/m2", &
        &              "accumulated grid-scale surface precipitation", &
        &              sfc_var_info, prm_diag%prec_gsp(:,:))

      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "H_ICE", "m", "sea ice depth", &
        &              sfc_var_info, prog_wtr%h_ice(:,:))

      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "CLCT", "-", "total cloud cover", &
        &              sfc_var_info, prm_diag%clct(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "CLCL", "-", "low level cloud cover", &
        &              sfc_var_info, prm_diag%clcl(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "CLCM", "-", "mid level cloud cover", &
        &              sfc_var_info, prm_diag%clcm(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "CLCH", "-", "high level cloud cover", &
        &              sfc_var_info, prm_diag%clch(:,:))

      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "hbas_con", "m", "height of convective cloud base",&
        &              sfc_var_info, prm_diag%hbas_con(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "htop_con", "m", "height of convective cloud top",&
        &              sfc_var_info, prm_diag%htop_con(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "UMFL_S", "N m-2", "u-momentum flux at the surface", &
        &              sfc_var_info, prm_diag%umfl_s(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "VMFL_S", "N m-2", "v-momentum flux at the surface", &
        &              sfc_var_info, prm_diag%vmfl_s(:,:))

      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SWDIFU_S", "W m-2", "shortwave upward flux at surface", &
        &              sfc_var_info, prm_diag%swflx_up_sfc(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SWDIFD_S", "W m-2", &
        &              "shortwave diffuse downward flux at surface", &
        &              sfc_var_info, prm_diag%swflx_dn_sfc_diff(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "PAB_S", "W m-2", &
        &       "photosynthetically active shortwave downward flux at surface", &
        &              sfc_var_info, prm_diag%swflx_par_sfc(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SNOWFRAC", "%", &
        &              "snow-cover fraction", &
        &              sfc_var_info, lnd_diag%snowfrac(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SNOWFRAC_LC", "%", &
        &              "snow-cover fraction (related to landuse class", &
        &              sfc_var_info, lnd_diag%snowfrac_lc(:,:))
      ! There is no 2D-field for SWDIR_S available, so we provide swflx_dn_sfc_diff instead
      ! The actual values will be diagnosed at a later point
      CALL add_sfc_var(meteogram_config, var_list, &
        &              IBSET(VAR_GROUP_SURFACE, FLAG_DIAG), &
        &              "SWDIR_S", "W m-2", &
        &              "shortwave direct downward flux at surface", &
        &              sfc_var_info, prm_diag%swflx_dn_sfc_diff(:,:))
    ELSE ! ie .NOT. iforcing == inwp
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RSDS", "W m-2", &
        &              "surface downwelling shortwave radiation", &
        &              sfc_var_info, prm_field(jg)%rsds(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RVDS_DIF", "W m-2", &
        &              "surface downwelling diffuse visible radiation", &
        &              sfc_var_info, prm_field(jg)%rvds_dif(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RNDS_DIF", "W m-2", &
        &              "surface downwelling diffuse near infrared radiation", &
        &              sfc_var_info, prm_field(jg)%rnds_dif(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RVDS_DIR", "W m-2", &
        &              "surface downwelling direct visible radiation", &
        &              sfc_var_info, prm_field(jg)%rvds_dir(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RNDS_DIR", "W m-2", &
        &              "surface downwelling direct near infrared radiation", &
        &              sfc_var_info, prm_field(jg)%rnds_dir(:,:))
    ENDIF ! iforcing == inwp

    ! -- tiled surface fields
    IF (meteogram_config%loutput_tiles) THEN     ! write some selected tile specific fields
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_SOIL_ML, &
        &               "W_SO_T", "m H2O", "soil water content (water + ice)", &
        &               var_info, lnd_prog%w_so_t(:,:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_SOIL_ML, &
        &               "W_SO_ICE_T", "m H2O", "ice content", &
        &               var_info, lnd_prog%w_so_ice_t(:,:,:,:))
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_SOIL_MLp2, &
        &               "T_SO_T", "K", "soil temperature", &
        &               var_info, lnd_prog%t_so_t(:,:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "T_G_T", "K", "surface temperature", &
        &              sfc_var_info, lnd_prog%t_g_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SHFL_T", "W/m2", "sensible heat flux (surface)", &
        &              sfc_var_info, prm_diag%shfl_s_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "LHFL_T", "W/m2", "latent heat flux (surface)", &
        &              sfc_var_info, prm_diag%lhfl_s_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "QHFL_S_T", "kg/m2", "evapotranspiration flux (surface)", &
        &              sfc_var_info, prm_diag%qhfl_s_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "H_SNOW_T", "m", "snow height", &
        &              sfc_var_info, lnd_diag%h_snow_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "W_SNOW_T", "m H2O", "water content of snow", &
        &              sfc_var_info, lnd_prog%w_snow_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "SOBS_T", "W m-2", "shortwave net flux (surface)", &
        &              sfc_var_info, prm_diag%swflxsfc_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "THBS_T", "W m-2", "longwave net flux (surface)", &
        &              sfc_var_info, prm_diag%lwflxsfc_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "FRAC_T", "-", "tile fractions (time dependent)", &
        &              sfc_var_info, ext_data%atm%frac_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "snowfrac_t", "%", &
        &              "local tile-based snow-cover fraction", &
        &              sfc_var_info, lnd_diag%snowfrac_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "snowfrac_lc_t", "%", &
        &              "snow-cover fraction per land-cover class (reduced by melting)", &
        &              sfc_var_info, lnd_diag%snowfrac_lc_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "snowfrac_lcu_t", "%", &
        &              "snow-cover fraction per land-cover class (melting unmodified)", &
        &              sfc_var_info, lnd_diag%snowfrac_lcu_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
       &              "W_I_T", "m H2O", & 
       &              "water content of interception water", &
       &              sfc_var_info, lnd_prog%w_i_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RUNOFF_S_T", "kg/m2", &
        &              "surface water runoff; sum over forecast", &
        &              sfc_var_info, lnd_diag%runoff_s_t(:,:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "RUNOFF_G_T", "kg/m2", &
        &              "soil water runoff; sum over forecast", &
        &              sfc_var_info, lnd_diag%runoff_g_t(:,:,:))
      IF (var_in_output(jg)%res_soilwatb) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &            "RESID_WSO_T", "kg/m2", &
          &            "residuum of the soil water budget", &
          &            sfc_var_info, lnd_diag%resid_wso_t(:,:,:))
      ENDIF
      IF (var_in_output(jg)%snow_melt) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &            "SNOW_MELT_FLUX_T", "kg m-2 s-1", &
          &            "snow melt flux tile", &
          &            sfc_var_info, lnd_diag%snow_melt_flux_t(:,:,:))
      ENDIF
    ENDIF

    ! -- vertical integrals

    CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
      &              "TQV", "kg m-2", "column integrated water vapour", &
      &              sfc_var_info, diag%tracer_vi_ptr(iqv)%p_2d(:,:))
    CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
      &              "TQC", "kg m-2", "column integrated cloud water", &
      &              sfc_var_info, diag%tracer_vi_ptr(iqc)%p_2d(:,:))
    CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
      &              "TQI", "kg m-2", "column integrated cloud ice", &
      &              sfc_var_info, diag%tracer_vi_ptr(iqi)%p_2d(:,:))
    IF ( iqm_max >= 4) THEN
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TQR", "kg m-2", "column integrated rain", &
        &              sfc_var_info, diag%tracer_vi_ptr(iqr)%p_2d(:,:))
    ENDIF
    IF ( iqm_max >= 5) THEN
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TQS", "kg m-2", "column integrated snow", &
        &              sfc_var_info, diag%tracer_vi_ptr(iqs)%p_2d(:,:))
    END IF

    IF ( iforcing == inwp ) THEN
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TQV_DIA", "kg m-2", &
        &              "total column integrated water vapour (diagnostic)", &
        &              sfc_var_info, prm_diag%tci_ptr(iqv)%p_2d(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TQC_DIA", "kg m-2", &
        &              "total column integrated cloud water (diagnostic)", &
        &              sfc_var_info, prm_diag%tci_ptr(iqc)%p_2d(:,:))
      CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &              "TQI_DIA", "kg m-2", &
        &              "total column integrated cloud ice (diagnostic)", &
        &              sfc_var_info, prm_diag%tci_ptr(iqi)%p_2d(:,:))
      IF (var_in_output(jg)%dbzlmx_low) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "DBZLMX_LOW", "dBZ", &
          &              "max radar reflectivity in layer 1000 - 2000 m AGL", &
          &              sfc_var_info, prm_diag%dbzlmx_low(:,:))
      END IF
      IF (var_in_output(jg)%dbz850) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "DBZ_850", "dBZ", &
          &              "radar reflectivity in 1500 m AGL", &
          &              sfc_var_info, prm_diag%dbz_850(:,:))
      END IF
      IF (var_in_output(jg)%dbzcmax) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "DBZ_CMAX", "dBZ", &
          &              "column max radar reflectivity", &
          &              sfc_var_info, prm_diag%dbz_cmax(:,:))
      END IF
      IF (var_in_output(jg)%dbzctmax) THEN
        CALL getPTStringFromMS(NINT(1000_wp*celltracks_interval(jg), i8), c_time_int)
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "DBZ_CTMAX", "dBZ", &
          &              "column and time max radar reflectivity since end of previous full "// &
          &              TRIM(ADJUSTL(c_time_int(3:)))//" since model start", &
          &              sfc_var_info, prm_diag%dbz_ctmax(:,:))
      END IF
      IF (var_in_output(jg)%lpi_max) THEN
        CALL getPTStringFromMS(NINT(1000_wp*celltracks_interval(jg), i8), c_time_int)
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "LPI_MAX", "J/kg", &
          &              "time max lightning potential index since end of previous full "// &
          &              TRIM(ADJUSTL(c_time_int(3:)))//" since model start", &
          &              sfc_var_info, prm_diag%lpi_max(:,:))
      END IF
      IF (var_in_output(jg)%ceiling) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "CEILING", "m", &
          &              "ceiling height", &
          &              sfc_var_info, prm_diag%ceiling_height(:,:))
      END IF
      IF (var_in_output(jg)%echotopinm) THEN
        DO i=1, echotop_meta(jg)%nechotop
          CALL getPTStringFromMS(NINT(1000_wp*echotop_meta(jg)%time_interval, i8), c_time_int)
          WRITE (c_thresh_int, '(i10," dBZ")') NINT(echotop_meta(jg)%dbzthresh(i))
          CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
            &              "ECHOTOPinM_"//TRIM(ADJUSTL(c_thresh_int(1:10))), "m", &
            &              "maximum height of exceeding radar reflectivity threshold "// &
            &              TRIM(ADJUSTL(c_thresh_int))//" since end of previous full "// &
            &              TRIM(ADJUSTL(c_time_int(3:)))//" since model start", &
            &              sfc_var_info, prm_diag%echotopinm(:,i,:))
        END DO
      END IF
      IF (var_in_output(jg)%vis) THEN
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "VIS", "m", &
          &              "near surface visibility", &
          &              sfc_var_info, prm_diag%vis(:,:))
      END IF
      IF (var_in_output(jg)%tot_pr_max) THEN
        CALL getPTStringFromMS(NINT(1000_wp*celltracks_interval(jg), i8), c_time_int)
        CALL add_sfc_var(meteogram_config, var_list, VAR_GROUP_SURFACE, &
          &              "TOT_PR_MAX", "kg/m2/s", &
          &              "time max precipitation rate since end of previous full "// &
          &              TRIM(ADJUSTL(c_time_int(3:)))//" since model start", &
          &              sfc_var_info, prm_diag%tot_pr_max(:,:))
      END IF
    ENDIF ! iforcing == nwp

    IF (inextra_2d > 0) THEN
      ! Variable: Extra 2D
      CALL add_sfc_var (meteogram_config, var_list, VAR_GROUP_SURFACE, &
        &               "EXTRA2D","","-", &
        &               sfc_var_info, diag%extra_2d(:,:,1:inextra_2d))
    ENDIF
    IF (inextra_3d > 0) THEN
      ! Variable: Extra 3D
      CALL add_atmo_var(meteogram_config, var_list, VAR_GROUP_ATMO_ML, &
        &               "EXTRA3D","","-", &
        &               var_info, diag%extra_3d(:,:,:,1:inextra_3d))
    END IF

#ifndef NOMPI
    ! collect variable info for pure I/O PEs
    IF (pack_buf%l_is_varlist_sender) &
      CALL pack_varlists(pack_buf, var_list, var_info, sfc_var_info)
#endif

    ! several variable indices, stored for convenience (when computing
    ! additional diagnostics):
    CALL setup_diag_var_indices(diag_var_indices, &
      var_info(1:var_list%no_atmo_vars)%cf, &
      sfc_var_info(1:var_list%no_sfc_vars)%cf)

  END SUBROUTINE meteogram_setup_variables

#ifndef NOMPI
  SUBROUTINE pack_varlists(pack_buf, var_list, var_info, sfc_var_info)
    TYPE(mtgrm_pack_buf), INTENT(inout) :: pack_buf
    TYPE(t_var), INTENT(in) :: var_list
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//":pack_varlists"
    INTEGER :: var_counts(2), ivar
    var_counts(1) = var_list%no_atmo_vars
    var_counts(2) = var_list%no_sfc_vars
    CALL p_pack_int_1d(var_counts, 2, pack_buf%msg_varlist, pack_buf%pos)

    IF (dbg_level > 0) &
      CALL message(routine, "collect variable info for pure I/O PEs")
    DO ivar = 1, var_list%no_atmo_vars

      CALL p_pack_string(var_info(ivar)%cf%standard_name, pack_buf%msg_varlist, pack_buf%pos)
      CALL p_pack_string(var_info(ivar)%cf%long_name,     pack_buf%msg_varlist, pack_buf%pos)
      CALL p_pack_string(var_info(ivar)%cf%units,         pack_buf%msg_varlist, pack_buf%pos)
      CALL p_pack_int   (var_info(ivar)%igroup_id,        pack_buf%msg_varlist, pack_buf%pos)
      CALL p_pack_int   (var_info(ivar)%nlevs,            pack_buf%msg_varlist, pack_buf%pos)
    END DO
    IF (dbg_level > 0) &
      CALL message(routine, "collect surface variable info for pure I/O PEs")
    DO ivar = 1, var_list%no_sfc_vars
      CALL p_pack_string(sfc_var_info(ivar)%cf%standard_name, pack_buf%msg_varlist, pack_buf%pos)
      CALL p_pack_string(sfc_var_info(ivar)%cf%long_name,     pack_buf%msg_varlist, pack_buf%pos)
      CALL p_pack_string(sfc_var_info(ivar)%cf%units,         pack_buf%msg_varlist, pack_buf%pos)
      CALL p_pack_int   (sfc_var_info(ivar)%igroup_id,        pack_buf%msg_varlist, pack_buf%pos)
    END DO
  END SUBROUTINE pack_varlists
#endif

  SUBROUTINE setup_diag_var_indices(diag_var_indices, cf_atmo, cf_lnd)
    !> indices at which to find variables for compute_diagnostics
    TYPE(meteogram_diag_var_indices), INTENT(out) :: diag_var_indices
    TYPE(t_cf_var), INTENT(in) :: cf_atmo(:), cf_lnd(:)

    diag_var_indices%i_T        = get_var("T"       , cf_atmo)
    diag_var_indices%i_QV       = get_var("QV"      , cf_atmo)
    diag_var_indices%i_REL_HUM  = get_var("REL_HUM" , cf_atmo)
    diag_var_indices%i_PEXNER   = get_var("PEXNER"  , cf_atmo)
    diag_var_indices%i_SWDIR_S  = get_var("SWDIR_S" , cf_lnd)
    diag_var_indices%i_ALB      = get_var("ALB"     , cf_lnd)
    diag_var_indices%i_SWDIFD_S = get_var("SWDIFD_S", cf_lnd)
    diag_var_indices%i_SOBS     = get_var("SOBS"    , cf_lnd)
  END SUBROUTINE setup_diag_var_indices

  !! @return Index of (3d) variable with given name.
  !!
  FUNCTION get_var(zname, cf)
    INTEGER :: get_var
    CHARACTER(LEN=*),  INTENT(IN) :: zname
    TYPE(t_cf_var), INTENT(in) :: cf(:)

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":get_var"
    INTEGER :: ivar, nvar

    get_var = -1 ! invalid result
    nvar = SIZE(cf)
    VAR_LOOP : DO ivar=1,nvar
      IF (cf(ivar)%standard_name == zname) THEN
        get_var = ivar
        EXIT VAR_LOOP
      END IF
    END DO VAR_LOOP
    ! the following consistency check is disabled (since the user may
    ! use the namelist parameter "var_list"):
    !
    ! IF (get_var == -1)  CALL finish (routine, 'Invalid name: '//TRIM(zname))
  END FUNCTION get_var

  !! Computation of additional diagnostic quantities for meteogram.
  !!
  SUBROUTINE compute_diagnostics(diag_var_indices, i_tstep, &
    ithis_nlocal_pts, out_buf, buf_idx)
    TYPE(meteogram_diag_var_indices), INTENT(in) :: diag_var_indices
    INTEGER, INTENT(IN) :: i_tstep   ! time step index
    TYPE(t_mtgrm_out_buffer), INTENT(inout) :: out_buf
    INTEGER, INTENT(in) :: ithis_nlocal_pts
    INTEGER, INTENT(in) :: buf_idx(ithis_nlocal_pts)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":compute_diagnostics"
    INTEGER                         :: ilev, nlevs, istation, istation_buf
    INTEGER                         :: i_REL_HUM, i_T, i_QV, i_PEXNER, &
      &                                i_SWDIR_S, i_ALB, i_SWDIFD_S, i_SOBS
    REAL(wp)                        :: temp, qv, p_ex
    REAL(wp)                        :: albedo, swdifd_s, sobs

    IF (diag_var_indices%i_REL_HUM == 0) RETURN

    ! TODO[FP] : In some cases, values (slightly) greater than 100%
    !            are computed for relative humidity.

    i_REL_HUM = diag_var_indices%i_REL_HUM
    IF (i_REL_HUM /= -1) THEN
      i_T = diag_var_indices%i_T
      i_QV = diag_var_indices%i_QV
      i_PEXNER = diag_var_indices%i_PEXNER

      IF (i_T /= -1 .AND. i_QV /= -1 .AND. i_PEXNER  /= -1) THEN
        nlevs = SIZE(out_buf%atmo_vars(i_T)%a, 2)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP COLLAPSE(2) PRIVATE(istation_buf, temp, qv, p_ex)
        DO ilev=1,nlevs
          DO istation = 1, ithis_nlocal_pts
            istation_buf = buf_idx(istation)
            ! get values for temperature, etc.:
            temp = out_buf%atmo_vars(i_T)%a(istation_buf, ilev, i_tstep)
            qv   = out_buf%atmo_vars(i_QV)%a(istation_buf, ilev, i_tstep)
            p_ex = out_buf%atmo_vars(i_PEXNER)%a(istation_buf, ilev, i_tstep)
            !-- compute relative humidity as r = e/e_s:
!CDIR NEXPAND
            out_buf%atmo_vars(i_REL_HUM)%a(istation_buf, ilev, i_tstep) &
              = rel_hum(temp, qv, p_ex)
          END DO
        END DO
        !$ACC END PARALLEL
      ELSE
        CALL message(routine, ">>> meteogram: REL_HUM could not be computed&
          & (T, QV, and/or PEXNER missing)")
      END IF
    END IF

    ! compute shortwave direct downward flux at surface
    i_SWDIR_S = diag_var_indices%i_SWDIR_S
    IF (i_SWDIR_S /= -1) THEN
      i_ALB = diag_var_indices%i_ALB
      i_SWDIFD_S = diag_var_indices%i_SWDIFD_S
      i_SOBS = diag_var_indices%i_SOBS
      IF (i_ALB /= -1 .AND. i_SWDIFD_S /= -1 .AND. i_SOBS /= -1) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(istation_buf, albedo, swdifd_s, sobs)
        DO istation = 1, ithis_nlocal_pts
          istation_buf = buf_idx(istation)
          albedo   = out_buf%sfc_vars(i_ALB)%a(istation_buf, i_tstep)
          swdifd_s = out_buf%sfc_vars(i_SWDIFD_S)%a(istation_buf, i_tstep)
          sobs     = out_buf%sfc_vars(i_SOBS)%a(istation_buf, i_tstep)
          out_buf%sfc_vars(i_SWDIR_S)%a(istation_buf, i_tstep) &
            = swdir_s(albedo, swdifd_s, sobs)
        END DO
        !$ACC END PARALLEL
      ELSE
        CALL message(routine, ">>> meteogram: SWDIR_S could not be computed&
          & (ALB, SWDIFD_S, and/or SOBS missing)")
      END IF
    END IF

  END SUBROUTINE compute_diagnostics


  !! Initialize meteogram data buffer, allocating storage.
  !! This is a collective operation.
  !!
  !! Note: Patch information, model state, etc. are optional
  !!       parameters here since this SR may also be called by pure
  !!       I/O PEs in asynchronous output mode.
  !!
  SUBROUTINE meteogram_init(meteogram_output_config, jg,     &
    &                       ptr_patch, ext_data, p_nh_state, &
    &                       prm_diag, p_lnd_state, prm_nwp_tend, iforcing, &
    &                       grid_uuid, number_of_grid_used)
    !> station data from namelist
    TYPE(t_meteogram_output_config), INTENT(INOUT) :: meteogram_output_config
    !> patch index
    INTEGER,                   INTENT(IN) :: jg
    !> data structure containing grid info:
    TYPE(t_patch),             INTENT(IN), OPTIONAL :: ptr_patch
    !> atmosphere external data
    TYPE(t_external_data),     INTENT(IN), OPTIONAL :: ext_data
    !> nonhydrostatic state
    TYPE(t_nh_state),          INTENT(IN), OPTIONAL :: p_nh_state
    !> physical model state and other auxiliary variables
    TYPE(t_nwp_phy_diag),      INTENT(IN), OPTIONAL :: prm_diag
    !> model state for the NWP land physics
    TYPE(t_lnd_state),         INTENT(IN), OPTIONAL :: p_lnd_state
    ! model state for physics tendencies
    TYPE(t_nwp_phy_tend),      INTENT(IN), OPTIONAL :: prm_nwp_tend
    !> parameterized forcing (right hand side) of dynamics, affects
    !! topography specification, see "mo_extpar_config"
    INTEGER,                   INTENT(IN), OPTIONAL :: iforcing

    TYPE(t_uuid),              INTENT(IN), OPTIONAL :: grid_uuid
    !> number of grid used
    INTEGER,                   INTENT(IN), OPTIONAL :: number_of_grid_used

    ! local variables:
    CHARACTER(*), PARAMETER :: routine = modname//":meteogram_init"

    INTEGER      :: ithis_nlocal_pts, nblks,          &
      &             nstations, ierrstat,              &
      &             jb, jc, istation, istation_buf, nstations_buf
    ! list of triangles containing lon-lat grid points (first dim: index and block)
    TYPE(mtgrm_pack_buf) :: pack_buf
    TYPE(t_mtgrm_invariants) :: invariants
    !> buffer size for var list
    INTEGER :: max_varlist_buf_size

    TYPE(t_var) :: var_list
    INTEGER      :: tri_idx(2,nproma,(meteogram_output_config%nstations+nproma-1)/nproma)
    INTEGER      :: max_time_stamps
    INTEGER                            :: io_collector_rank, iowner, &
      io_invar_send_rank
    LOGICAL :: is_io, is_mpi_workroot, is_mpi_test, is_pure_io_pe

    max_varlist_buf_size &
      = 2*p_int_byte + max_nvars*(3*max_descr_length+5*p_int_byte)

    is_io = my_process_is_io()
    is_mpi_workroot = my_process_is_mpi_workroot()
    is_mpi_test = my_process_is_mpi_test()
    !-- define the different roles in the MPI communication inside
    !-- this module

    ! Flag. True, if this PE is a pure I/O PE without own patch data:
    is_pure_io_pe = use_async_name_list_io .AND. is_io

    io_collector_rank = -1
    IF (use_async_name_list_io) THEN
      io_collector_rank &
        = process_mpi_io_size - 1 - MOD(jg - 1, process_mpi_io_size)
    END IF
    meteogram_output_config%io_proc_id = io_collector_rank

    mtgrm(jg)%time_offset = 0
    ! Flag. True, if this PE collects data from (other) working PEs
    mtgrm(jg)%l_is_collecting_pe                                             &
      &    =       (.NOT. meteogram_output_config%ldistributed)              &
      &      .AND. (.NOT. is_mpi_test)                                       &
      &      .AND. (     ((.NOT. use_async_name_list_io)                     &
      &                   .AND. is_mpi_workroot                              &
      &                   .AND. (p_n_work > 1) )                             &
      &             .OR. (is_pure_io_pe                             &
      &                   .AND. (p_pe_work == io_collector_rank)))
    IF (.NOT. meteogram_output_config%ldistributed) THEN
      IF (.NOT. use_async_name_list_io) THEN
        mtgrm(jg)%io_collect_comm = p_comm_work
        mtgrm(jg)%io_collector_rank = 0
      ELSE
        mtgrm(jg)%io_collect_comm = p_comm_work_2_io
        mtgrm(jg)%io_collector_rank = io_collector_rank
      END IF
    END IF

    ! Consistency check I: If this is NOT a pure I/O PE, then patch data
    ! must be available:
    IF (      .NOT. is_pure_io_pe                                        &
      & .AND. .NOT. (      PRESENT(ptr_patch)    .AND. PRESENT(ext_data)   &
      &              .AND. PRESENT(p_nh_state)   .AND. PRESENT(p_lnd_state)&
      &              .AND. PRESENT(prm_nwp_tend) .AND. PRESENT(iforcing)) ) THEN
      CALL finish (routine, 'Missing argument(s)!')
    END IF

    ! Consistency check II: If this is a pure I/O PE, then number_of_grid_used
    ! and grid_uuid must be available.
    IF (is_pure_io_pe &
      &.AND. .NOT. (PRESENT(number_of_grid_used) .AND. PRESENT(grid_uuid))) THEN
      CALL finish (routine, 'I/O PE Missing argument(s)!')
    ENDIF


    IF (ALLOCATED(tile_list%tile)) THEN
      invariants%ntiles_mtgrm = ntiles_total + ntiles_water
    ELSE
      invariants%ntiles_mtgrm = 1
    ENDIF

    max_time_stamps = meteogram_output_config%max_time_stamps
    mtgrm(jg)%max_time_stamps = max_time_stamps
    mtgrm(jg)%silent_flush = meteogram_output_config%silent_flush

    ! set meta data (appears in NetCDF output file)
    mtgrm(jg)%meteogram_file_info%cf       = cf_global_info
    mtgrm(jg)%meteogram_file_info%cf%title = 'ICON Meteogram File'

    mtgrm(jg)%meteogram_file_info%ldistributed = meteogram_output_config%ldistributed
    CALL uuid_unparse(grid_uuid, mtgrm(jg)%meteogram_file_info%uuid_string)
    mtgrm(jg)%meteogram_file_info%number_of_grid_used = number_of_grid_used

    ! ------------------------------------------------------------
    ! Distribute stations, determine number of stations located on
    ! this PE:
    ! ------------------------------------------------------------

    nstations = meteogram_output_config%nstations

    IF (.NOT. is_pure_io_pe) THEN
      nblks = (nstations+nproma-1)/nproma
      CALL mtgrm_station_owners(mtgrm(jg)%pstation, &
        mtgrm(jg)%nstations_local, meteogram_output_config, &
        nblks, ptr_patch, mtgrm(jg)%io_collector_rank, &
        mtgrm(jg)%io_collect_comm, tri_idx, mtgrm(jg)%global_idx, &
        mtgrm(jg)%io_invar_send_rank)
      pack_buf%l_is_varlist_sender = mtgrm(jg)%io_invar_send_rank == p_pe_work
      ithis_nlocal_pts = mtgrm(jg)%nstations_local
    ELSE
      ithis_nlocal_pts = 0
      pack_buf%l_is_varlist_sender = .FALSE.
      mtgrm(jg)%l_is_sender = .FALSE.
      IF (is_io .AND. io_collector_rank == p_pe_work) THEN
        CALL p_irecv(mtgrm(jg)%pstation, MPI_ANY_SOURCE, tag_mtgrm_msg+jg-1, &
          &          comm=mtgrm(jg)%io_collect_comm)
      ELSE
        RETURN
      END IF
    END IF

    ! Pure I/O PEs must receive all variable info from elsewhere.
    ! Here, they get it from working PE#0 which has collected it in
    ! "msg_varlist_buffer" during the add_xxx_var calls
    IF (     pack_buf%l_is_varlist_sender &
      & .OR. (is_pure_io_pe .AND. mtgrm(jg)%l_is_collecting_pe)) THEN
      ALLOCATE(pack_buf%msg_varlist(max_varlist_buf_size), stat=ierrstat)
      IF (ierrstat /= SUCCESS) &
        CALL finish (routine, 'ALLOCATE of MPI buffer failed.')
      pack_buf%pos = 0

      IF (is_pure_io_pe) THEN
        io_invar_send_rank = p_comm_remote_size(mtgrm(jg)%io_collect_comm) - 1
        CALL p_wait()
        DO istation = nstations, 1, -1
          iowner = mtgrm(jg)%pstation(istation)
          IF (iowner /= -1) THEN
            io_invar_send_rank = iowner
            EXIT
          END IF
        END DO
        mtgrm(jg)%io_invar_send_rank = io_invar_send_rank
        ! launch message receive call
        CALL p_irecv_packed(pack_buf%msg_varlist, mtgrm(jg)%io_invar_send_rank,&
          &                 TAG_VARLIST, max_varlist_buf_size, &
          &                 comm=mtgrm(jg)%io_collect_comm)
      END IF
    END IF

    ! ------------------------------------------------------------
    ! Initialize local data structure, fill header
    ! ------------------------------------------------------------

    mtgrm(jg)%icurrent  = 0 ! reset current sample index
    invariants%max_nlevs = 1

    ! set up list of variables:
    var_list%no_atmo_vars = 0
    var_list%no_sfc_vars  = 0
    IF (     ithis_nlocal_pts > 0 &
      & .OR. mtgrm(jg)%io_invar_send_rank == p_pe_work &
      & .OR. (.NOT. use_async_name_list_io .AND. is_mpi_workroot)) THEN
      ALLOCATE(mtgrm(jg)%sfc_var_info(MAX_NSFCVARS),  &
        &      mtgrm(jg)%var_info(MAX_NVARS),         &
        &      stat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, &
        'ALLOCATE of meteogram data structures failed (part 1)')
      CALL meteogram_setup_variables(meteogram_output_config, jg, iforcing, &
        &                            ext_data, &
        &                            p_nh_state%prog(nnow(jg)), &
        &                            p_nh_state%diag, prm_diag, &
        &                            p_lnd_state%prog_lnd(nnow(jg)), &
        &                            p_lnd_state%diag_lnd, &
        &                            p_lnd_state%prog_wtr(nnow(jg)), &
        &                            prm_nwp_tend, &
        &                            atm_phy_nwp_config(jg), &
        &                            mtgrm(jg)%var_info, &
        &                            mtgrm(jg)%sfc_var_info, &
        &                            mtgrm(jg)%diag_var_indices, &
        &                            var_list, pack_buf)
      CALL resize_var_info(var_list, mtgrm(jg)%var_info, mtgrm(jg)%sfc_var_info)

      CALL acc_copyin_var_info(mtgrm(jg)%var_info, mtgrm(jg)%sfc_var_info)

      IF (pack_buf%l_is_varlist_sender) THEN
        CALL p_send_packed(pack_buf%msg_varlist, mtgrm(jg)%io_collector_rank, &
          &                TAG_VARLIST, pack_buf%pos, &
          &                comm=mtgrm(jg)%io_collect_comm)
      END IF
    ELSE IF (mtgrm(jg)%l_is_collecting_pe) THEN
      CALL receive_var_info(mtgrm(jg)%var_info, mtgrm(jg)%sfc_var_info, &
        var_list, pack_buf)
    ELSE
      RETURN
    END IF

    IF (     pack_buf%l_is_varlist_sender &
      & .OR. (is_pure_io_pe .AND. mtgrm(jg)%l_is_collecting_pe)) THEN
      ! deallocate buffer
      DEALLOCATE(pack_buf%msg_varlist, stat=ierrstat)
      IF (ierrstat /= SUCCESS) &
        CALL finish (routine, 'DEALLOCATE of MPI buffer failed.')
    END IF

    nstations_buf = MERGE(nstations, &
      mtgrm(jg)%nstations_local, mtgrm(jg)%l_is_collecting_pe)
    CALL allocate_out_buf(mtgrm(jg)%out_buf, &
      mtgrm(jg)%var_info, mtgrm(jg)%sfc_var_info, &
      meteogram_output_config%max_time_stamps, nstations_buf)
    CALL allocate_heights(invariants%heights, mtgrm(jg)%var_info, nstations_buf)
    ALLOCATE(mtgrm(jg)%istep(max_time_stamps), &
      invariants%tile_frac(invariants%ntiles_mtgrm, nstations_buf), &
      invariants%tile_luclass(invariants%ntiles_mtgrm, nstations_buf), &
      invariants%soiltype(nstations_buf), &
      invariants%fc(nstations_buf), &
      invariants%frland(nstations_buf), &
      invariants%hsurf(nstations_buf), &
      invariants%tri_idx(nstations_buf,2), &
      mtgrm(jg)%zdate(max_time_stamps), stat=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, &
      'ALLOCATE of meteogram time stamp data structure failed')

    invariants%max_nlevs = &
      & MAX(0, MAXVAL(mtgrm(jg)%var_info(1:var_list%no_atmo_vars)%nlevs))

    ! set up list of local stations:
    IF (.NOT. is_pure_io_pe) THEN

      ALLOCATE(mtgrm(jg)%tri_idx_local(ithis_nlocal_pts, 2), &
        stat=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish (routine, 'ALLOCATE of meteogram data structures failed (part 3)')
      ENDIF
      DO istation=1,ithis_nlocal_pts
        jb = (istation-1)/nproma + 1
        jc = MOD(istation-1, nproma)+1
        istation_buf = MERGE(mtgrm(jg)%global_idx(istation), istation, &
          mtgrm(jg)%l_is_collecting_pe)
        mtgrm(jg)%tri_idx_local(istation, :) = tri_idx(:,jc,jb)
        CALL sample_station_init(istation_buf, tri_idx(:,jc,jb), &
          ptr_patch%cells, ext_data%atm, iforcing, p_nh_state%metrics, &
          mtgrm(jg)%var_info, mtgrm(jg)%sfc_var_info, invariants)
      END DO
      IF (      .NOT. mtgrm(jg)%l_is_collecting_pe &
        & .AND. .NOT. meteogram_output_config%ldistributed) &
        CALL send_time_invariants(mtgrm(jg)%var_info, invariants, &
        mtgrm(jg)%io_collector_rank, mtgrm(jg)%io_collect_comm, &
        mtgrm(jg)%global_idx(1:ithis_nlocal_pts))
      IF (.NOT. mtgrm(jg)%l_is_collecting_pe) THEN
        DO istation=1,ithis_nlocal_pts
          mtgrm(jg)%out_buf%station_idx(istation) &
            = mtgrm(jg)%global_idx(istation)
        END DO
      END IF
    END IF

    !$ACC ENTER DATA COPYIN(mtgrm(jg)%tri_idx_local) ASYNC(1)

    ! ------------------------------------------------------------
    ! If this is the IO PE: initialize global data structure
    ! ------------------------------------------------------------

    IO_PE : IF (mtgrm(jg)%l_is_collecting_pe) THEN
      DO istation = 1, nstations
        mtgrm(jg)%out_buf%station_idx(istation) = istation
      END DO
      CALL recv_time_invariants(mtgrm(jg)%var_info, invariants, &
        mtgrm(jg)%pstation, mtgrm(jg)%io_collect_comm, is_pure_io_pe)
    END IF IO_PE

    ! ------------------------------------------------------------
    ! initialize MPI buffer
    ! ------------------------------------------------------------

    ! Flag. True, if this PE sends data to a collector via MPI
    mtgrm(jg)%l_is_sender   = .NOT. meteogram_output_config%ldistributed &
      &                 .AND. .NOT. mtgrm(jg)%l_is_collecting_pe         &
      &                 .AND. .NOT. is_mpi_test                          &
      &                 .AND. (ANY(mtgrm(jg)%pstation == p_pe_work)      &
      &                        .OR. pack_buf%l_is_varlist_sender)

    ! Flag. True, if this PE writes data to file
    mtgrm(jg)%l_is_writer &
      & =      (      meteogram_output_config%ldistributed &
      &         .AND. mtgrm(jg)%nstations_local > 0) &
      &   .OR. mtgrm(jg)%l_is_collecting_pe

    IF (.NOT. meteogram_output_config%ldistributed) THEN
      ! compute maximum buffer size for MPI messages:
      mtgrm(jg)%max_buf_size                                                 &
        = station_pack_size(max_time_stamps, mtgrm(jg)%var_info,             &
        &                   mtgrm(jg)%sfc_var_info, mtgrm(jg)%io_collect_comm)
    END IF

    ! ------------------------------------------------------------
    ! If this is the IO PE: open NetCDF file
    ! ------------------------------------------------------------
    CALL meteogram_open_file(meteogram_output_config, mtgrm(jg), jg, invariants)

  END SUBROUTINE meteogram_init

  SUBROUTINE mtgrm_station_owners(pstation, ithis_nlocal_pts, output_config, &
    nblks, ptr_patch, io_collector_rank, io_collect_comm, tri_idx, &
    global_idx, io_invar_send_rank)
    INTEGER, INTENT(out) :: pstation(:), ithis_nlocal_pts
    !> station data from namelist
    TYPE(t_meteogram_output_config), INTENT(in) :: output_config
    INTEGER, INTENT(in) :: nblks, io_collector_rank, io_collect_comm
    INTEGER, INTENT(out) :: tri_idx(2,nproma,nblks), global_idx(:), &
      io_invar_send_rank
    !> data structure containing grid info:
    TYPE(t_patch), INTENT(in) :: ptr_patch

    INTEGER :: istation, nstations, iowner
    INTEGER :: jc, jb, npromz
    !> minimal distance
    REAL(gk) :: min_dist(nproma,nblks)
    !> geographical locations
    REAL(gk) :: in_points(nproma,nblks,2)
    REAL(wp) :: grid_sphere_radius_mtg
    TYPE(t_gnat_tree) :: gnat

    ! build an array of geographical coordinates from station list:
    ! in_points(...)
    nstations = output_config%nstations
    DO istation=1,nstations
      jc = MOD(istation-1,nproma)+1
      jb = (istation+nproma-1)/nproma
      in_points(jc,jb,1) &
        = output_config%station_list(istation)%location%lon * pi_180
      in_points(jc,jb,2) &
        = output_config%station_list(istation)%location%lat * pi_180
    END DO

    ! build GNAT data structure
    CALL gnat_init_grid(gnat, ptr_patch)
    ! perform proximity query

    IF (is_plane_torus) THEN
      grid_sphere_radius_mtg = ptr_patch%geometry_info%domain_length / (2*pi)
    ELSE
      grid_sphere_radius_mtg = grid_sphere_radius
    END IF

    npromz = MOD(nstations-1,nproma)+1
    CALL gnat_query_containing_triangles(gnat, ptr_patch, in_points(:,:,:),             &
      &                                  nproma, nblks, npromz, grid_sphere_radius_mtg, &
      &                                  p_test_run, tri_idx(:,:,:), min_dist(:,:))

    CALL gnat_merge_distributed_queries(ptr_patch, nstations, nproma, nblks, min_dist,  &
      &                                 tri_idx(:,:,:), in_points(:,:,:),               &
      &                                 global_idx(:), ithis_nlocal_pts)

    ! clean up
    CALL gnat_destroy(gnat)

    ! build a list of "owner" PEs for the stations:
    pstation(:) = -1
    DO istation = 1, ithis_nlocal_pts
      pstation(global_idx(istation)) = p_pe_work
    END DO
    ! All-to-all communicate the located stations to find out if
    ! some stations are "owned" by no PE at all (in the case of
    ! regional nests):
    CALL p_allreduce_max(pstation, comm=p_comm_work)
    io_invar_send_rank = -1
    IF (use_async_name_list_io) THEN
      io_invar_send_rank = -1
      DO istation = nstations, 1, -1
        iowner = pstation(istation)
        IF (iowner /= -1) THEN
          ! PE collecting variable info to send it to pure I/O PEs.
          ! (only relevant if pure I/O PEs exist)
          io_invar_send_rank = iowner
          EXIT
        END IF
      END DO
      IF (io_invar_send_rank == -1) THEN
        io_invar_send_rank = p_n_work - 1
        IF (io_invar_send_rank == p_pe_work) &
          WRITE (0, *) &
          'warning: potential problem in meteogram output of patch ', &
          ptr_patch%id, ': no station found to be in grid'
      END IF
      IF (io_invar_send_rank == p_pe_work) THEN
        CALL p_send(pstation, io_collector_rank, tag_mtgrm_msg+ptr_patch%id-1, &
          &         comm=io_collect_comm)
      END IF
    END IF

  END SUBROUTINE mtgrm_station_owners

  SUBROUTINE sample_station_init(istation_buf, tri_idx, &
    cells, atm, iforcing, metrics, var_info, sfc_var_info, invariants)
    INTEGER, INTENT(in) :: istation_buf
    INTEGER, INTENT(in) :: tri_idx(2)
    TYPE(t_grid_cells), INTENT(in) :: cells
    TYPE(t_external_atmos), INTENT(in) :: atm
    TYPE(t_mtgrm_invariants), TARGET, INTENT(inout) :: invariants
    !> parameterized forcing (right hand side) of dynamics, affects
    !! topography specification, see "mo_extpar_config"
    INTEGER, INTENT(IN) :: iforcing
    ! nonhydrostatic state
    TYPE(t_nh_metrics), INTENT(IN) :: metrics
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)

    REAL(wp), POINTER :: station_var_heights(:)
    INTEGER :: tri_idx1, tri_idx2, glb_index, ivar, nvars, nlevs
    CHARACTER(len=*), PARAMETER :: routine = modname//"::sample_station_init"

    ! set local triangle index, block:
    tri_idx1 = tri_idx(1)
    tri_idx2 = tri_idx(2)
    ! translate local index to global index:
    glb_index = cells%decomp_info%glb_index(idx_1d(tri_idx1, tri_idx2))
    invariants%tri_idx(istation_buf, 1) = idx_no(glb_index)
    invariants%tri_idx(istation_buf, 2) = blk_no(glb_index)
    ! set Coriolis parameter for station
    invariants%fc(istation_buf) = cells%f_c(tri_idx1, tri_idx2)

    !
    ! set station information on height, soil type etc.:
    SELECT CASE ( iforcing )
    CASE ( inwp ) ! NWP physics
      invariants%hsurf(istation_buf) &
        &        =  atm%topography_c(tri_idx1, tri_idx2)
      invariants%frland(istation_buf) &
        &        =  atm%fr_land(tri_idx1, tri_idx2)
      invariants%soiltype(istation_buf) =  atm%soiltyp(tri_idx1, tri_idx2)
      !
      invariants%tile_frac(:, istation_buf) &
        &        = atm%lc_frac_t(tri_idx1, tri_idx2, 1:invariants%ntiles_mtgrm)
      invariants%tile_luclass(:, istation_buf) &
        &        = atm%lc_class_t(tri_idx1, tri_idx2, 1:invariants%ntiles_mtgrm)

    CASE DEFAULT
      !!! invariants%hsurf(istation_buf)    =  0._wp
      invariants%hsurf(istation_buf) &
        &        =  atm%topography_c(tri_idx1, tri_idx2)
      invariants%frland(istation_buf)   =  0._wp
      invariants%soiltype(istation_buf) =  0
      !
      invariants%tile_frac(:, istation_buf) = 0._wp
      invariants%tile_luclass(:, istation_buf) = 0
    END SELECT

    ! initialize value buffer and set level heights:
    nvars = SIZE(var_info)
    DO ivar = 1, nvars
      nlevs = var_info(ivar)%nlevs
      station_var_heights => invariants%heights(ivar)%a(istation_buf, :)
      ! initialize level heights:
      SELECT CASE(IBCLR(var_info(ivar)%igroup_id, FLAG_DIAG))
      CASE(VAR_GROUP_ATMO_ML)
        ! model level heights
        station_var_heights = metrics%z_mc(tri_idx1, 1:nlevs, tri_idx2)
      CASE(VAR_GROUP_ATMO_HL)
        ! half level heights
        station_var_heights = metrics%z_ifc(tri_idx1, 1:nlevs, tri_idx2)
      CASE(VAR_GROUP_SOIL_ML)
        ! soil half level heights
        station_var_heights = zml_soil(1:nlevs)
      CASE(VAR_GROUP_SOIL_MLp2)
        ! soil half level heights PLUS surface level
        station_var_heights(1) = 0._wp
        station_var_heights(2:nlevs) = zml_soil(1:nlevs-1)
      CASE DEFAULT
        CALL finish (routine, 'Invalid group ID.')
      END SELECT
    END DO

  END SUBROUTINE sample_station_init

  SUBROUTINE allocate_out_buf(out_buf, var_info, sfc_var_info, &
    max_time_stamps, nstations)
    TYPE(t_mtgrm_out_buffer), INTENT(inout) :: out_buf
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)
    INTEGER, INTENT(in) :: max_time_stamps, nstations

    INTEGER :: ivar, nvars_atmo, nvars_sfc, nlevs, ierror
    CHARACTER(len=*), PARAMETER :: &
      routine = modname//"::allocate_out_buf"

    nvars_atmo = SIZE(var_info)
    nvars_sfc = SIZE(sfc_var_info)
    ALLOCATE(out_buf%atmo_vars(nvars_atmo), out_buf%sfc_vars(nvars_sfc), &
      out_buf%station_idx(nstations), stat=ierror)
    IF (ierror == success) THEN
      !$ACC ENTER DATA CREATE(out_buf) ASYNC(1)
      !$ACC ENTER DATA CREATE(out_buf%atmo_vars, out_buf%sfc_vars) ASYNC(1)
      DO ivar = 1, nvars_atmo
        nlevs = var_info(ivar)%nlevs
        ALLOCATE(out_buf%atmo_vars(ivar)%a(nstations,nlevs,max_time_stamps))
        IF (ierror /= success) EXIT
        !$ACC ENTER DATA CREATE(out_buf%atmo_vars(ivar)%a) ASYNC(1)
      END DO
    END IF
    IF (ierror == success) THEN
      DO ivar = 1, nvars_sfc
        ALLOCATE(out_buf%sfc_vars(ivar)%a(nstations,max_time_stamps))
        IF (ierror /= success) EXIT
        !$ACC ENTER DATA CREATE(out_buf%sfc_vars(ivar)%a) ASYNC(1)
      END DO
    END IF
    IF (ierror /= SUCCESS) CALL finish(routine, &
      'ALLOCATE of meteogram data structures failed')
  END SUBROUTINE allocate_out_buf

  SUBROUTINE allocate_heights(heights, var_info, nstations)
    TYPE(t_a_2d), ALLOCATABLE, INTENT(inout) :: heights(:)
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    INTEGER, INTENT(in) :: nstations
    INTEGER :: ivar, nvars, nlevs
    nvars = SIZE(var_info)
    ALLOCATE(heights(nvars))
    DO ivar = 1, nvars
      nlevs = var_info(ivar)%nlevs
      ALLOCATE(heights(ivar)%a(nstations,nlevs))
    END DO
  END SUBROUTINE allocate_heights

  !! @return .TRUE. if meteogram data will be recorded for this step.
  !!
  FUNCTION meteogram_is_sample_step(meteogram_output_config, cur_step)
    LOGICAL :: meteogram_is_sample_step
    ! station data from namelist
    TYPE(t_meteogram_output_config), INTENT(IN) :: meteogram_output_config
    !> current model iteration step
    INTEGER,          INTENT(IN)  :: cur_step

    meteogram_is_sample_step = &
      &  meteogram_output_config%lenabled               .AND. &
      &  (cur_step >= meteogram_output_config%n0_mtgrm) .AND. &
      &  (MOD((cur_step - meteogram_output_config%n0_mtgrm),  &
      &       meteogram_output_config%ninc_mtgrm) == 0)

  END FUNCTION meteogram_is_sample_step


  !! Adds values for current model time to buffer. Also flushes the
  !! buffer when full and adds time stamp
  !!
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation.
  !!
  SUBROUTINE meteogram_sample_vars(jg, cur_step, cur_datetime, lacc)
    !> patch index
    INTEGER,          INTENT(IN)  :: jg
    !> current model iteration step
    INTEGER,          INTENT(IN)  :: cur_step
    !> date and time of point sample
    TYPE(datetime), INTENT(IN) :: cur_datetime
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":meteogram_sample_vars"
    INTEGER :: istation, ithis_nlocal_pts, i_tstep
    CHARACTER(len=MAX_DATETIME_STR_LEN) :: zdate
    INTEGER :: msg(2)
    INTEGER :: buf_idx(mtgrm(jg)%nstations_local)

    CALL assert_acc_device_only(routine, lacc)

    IF (dbg_level > 0) THEN
      WRITE(message_text,*) "Sampling at step=", cur_step
      CALL message(routine, TRIM(message_text))
    END IF

    ! increase time step counter
    IF (mtgrm(jg)%icurrent >= mtgrm(jg)%max_time_stamps) THEN
      ! buffer full
      IF (.NOT. mtgrm(jg)%silent_flush) THEN
        CALL message(routine, 'WARNING: Intermediate meteogram flush. &
             &Is the sampling buffer too small?')
      ENDIF
      IF (use_async_name_list_io .AND. p_pe_work == 0) THEN
        msg(1) = msg_io_meteogram_flush
        msg(2) = jg
        CALL p_send(msg, mtgrm(jg)%io_collector_rank, 0, &
          comm=mtgrm(jg)%io_collect_comm)
      END IF

      CALL meteogram_flush_file(jg, lacc=.TRUE.)
    END IF
    i_tstep = mtgrm(jg)%icurrent + 1
    mtgrm(jg)%icurrent = i_tstep

    ithis_nlocal_pts = mtgrm(jg)%nstations_local
    IF (mtgrm(jg)%l_is_writer &
      .OR. (mtgrm(jg)%l_is_sender &
      &     .AND. p_pe_work == mtgrm(jg)%io_invar_send_rank)) THEN
      CALL datetimeToPosixString(cur_datetime, zdate, "%Y%m%dT%H%M%SZ")
      mtgrm(jg)%zdate(i_tstep) = zdate
      mtgrm(jg)%istep(i_tstep) = cur_step
    END IF
    IF (ithis_nlocal_pts > 0 .OR. mtgrm(jg)%l_is_collecting_pe) THEN
      IF (mtgrm(jg)%l_is_collecting_pe) THEN
        buf_idx = mtgrm(jg)%global_idx(1:ithis_nlocal_pts)
      ELSE
        DO istation = 1, ithis_nlocal_pts
          buf_idx(istation) = istation
        END DO
      END IF
      !$ACC DATA COPYIN(buf_idx)
      ! fill time step with values
      !$ACC WAIT
      CALL sample_station_vars(&
        mtgrm(jg)%tri_idx_local, &
        mtgrm(jg)%var_info, mtgrm(jg)%sfc_var_info, &
        mtgrm(jg)%diag_var_indices, i_tstep, ithis_nlocal_pts, &
        mtgrm(jg)%out_buf, buf_idx)
      !$ACC WAIT
      !$ACC END DATA
    END IF
  END SUBROUTINE meteogram_sample_vars

  !> transfers data for the current time step for each station X
  !! variable + plus diagnostics into output buffer
  SUBROUTINE sample_station_vars(tri_idx_local, var_info, sfc_var_info, &
    diag_var_indices, i_tstep, ithis_nlocal_pts, out_buf, buf_idx)
    INTEGER, INTENT(in) :: tri_idx_local(:,:)
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)
    TYPE(meteogram_diag_var_indices), INTENT(in) :: diag_var_indices
    INTEGER, INTENT(in) :: i_tstep, ithis_nlocal_pts
    TYPE(t_mtgrm_out_buffer), INTENT(inout) :: out_buf
    INTEGER, INTENT(in) :: buf_idx(ithis_nlocal_pts)

    INTEGER :: iidx, iblk, ivar, nvars, istation, istation_buf
    CHARACTER(len=*), PARAMETER :: routine &
      = modname//'::sample_station_vars'

    REAL(wp), PARAMETER :: z10olog10 = 10.0_wp / LOG(10.0_wp)
    REAL(wp), PARAMETER :: eps_dbz   = 1.0e-15_wp

    ! sample 3D variables:
    nvars = SIZE(var_info)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG
    VAR_LOOP : DO ivar=1,nvars
      IF (.NOT. BTEST(var_info(ivar)%igroup_id, FLAG_DIAG)) THEN
        !$ACC LOOP VECTOR PRIVATE(istation_buf, iidx, iblk)
        DO istation = 1, ithis_nlocal_pts
          istation_buf = buf_idx(istation)
          iidx  = tri_idx_local(istation, 1)
          iblk  = tri_idx_local(istation, 2)
          out_buf%atmo_vars(ivar)%a(istation_buf, :, i_tstep) &
            = var_info(ivar)%p_source(iidx, :, iblk)
        END DO
      END IF
    END DO VAR_LOOP
    !$ACC END PARALLEL

    ! sample surface variables:
    nvars = SIZE(sfc_var_info)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG
    SFCVAR_LOOP : DO ivar=1,nvars
      IF (.NOT. BTEST(sfc_var_info(ivar)%igroup_id, FLAG_DIAG)) THEN
        !$ACC LOOP VECTOR PRIVATE(istation_buf, iidx, iblk)
        DO istation = 1, ithis_nlocal_pts
          istation_buf = buf_idx(istation)
          iidx  = tri_idx_local(istation, 1)
          iblk  = tri_idx_local(istation, 2)
          out_buf%sfc_vars(ivar)%a(istation_buf, i_tstep) &
            = sfc_var_info(ivar)%p_source(iidx, iblk)
        END DO
        ! convert some units of variables to the desired output units:
!!! THE FOLLOWING WOULD BE THE CONVENIENT FORTRAN SOLUTION,
!!! BUT TRIM() IN THIS CONTEXT SEEMS NOT TO BE SUPPORTED ON DAINT_GPU
!        SELECT CASE (TRIM(sfc_var_info(ivar)%cf%standard_name))
!        CASE ('DBZLMX_LOW','DBZ_850','DBZ_CMAX','DBZ_CTMAX')
!          !$ACC LOOP VECTOR PRIVATE(istation_buf)
!          DO istation = 1, ithis_nlocal_pts
!            istation_buf = buf_idx(istation)
!            out_buf%sfc_vars(ivar)%a(istation_buf, i_tstep) &
!              =  z10olog10 * LOG( MAX(out_buf%sfc_vars(ivar)%a(istation_buf, i_tstep), eps_dbz) )
!          END DO
!        END SELECT
!!! SO WE TRY THIS HACK INSTEAD:
        IF (sfc_var_info(ivar)%cf%standard_name(1:3) == 'DBZ') THEN
          !$ACC LOOP VECTOR PRIVATE(istation_buf)
          DO istation = 1, ithis_nlocal_pts
            istation_buf = buf_idx(istation)
            out_buf%sfc_vars(ivar)%a(istation_buf, i_tstep) &
              =  z10olog10 * LOG( MAX(out_buf%sfc_vars(ivar)%a(istation_buf, i_tstep), eps_dbz) )
          END DO
        END IF
      END IF
    END DO SFCVAR_LOOP
    !$ACC END PARALLEL

    ! compute additional diagnostic quantities:
    CALL compute_diagnostics(diag_var_indices, i_tstep, &
      ithis_nlocal_pts, out_buf, buf_idx)

  END SUBROUTINE sample_station_vars

  !! Destroy meteogram data structure.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation.
  !!
  SUBROUTINE meteogram_finalize(jg, lacc)
    INTEGER, INTENT(IN)         :: jg    !< patch index
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
    ! local variables:
    CHARACTER(*), PARAMETER     :: routine = modname//":meteogram_finalize"
    INTEGER                     :: ierror
    LOGICAL :: is_mpi_workroot

    ! ------------------------------------------------------------
    ! flush, and if this is the IO PE, close NetCDF file
    ! ------------------------------------------------------------
    CALL meteogram_close_file(jg, lacc=lacc)

    is_mpi_workroot = my_process_is_mpi_workroot()
    IF (     mtgrm(jg)%nstations_local > 0 &
      & .OR. (.NOT. use_async_name_list_io .AND. is_mpi_workroot) &
      & .OR. mtgrm(jg)%l_is_collecting_pe) THEN
      DEALLOCATE(mtgrm(jg)%istep, mtgrm(jg)%zdate, mtgrm(jg)%var_info, &
        mtgrm(jg)%sfc_var_info, stat=ierror)
      IF (ierror /= SUCCESS) &
        CALL finish (routine, 'DEALLOCATE of meteogram data structures failed')
    END IF
  END SUBROUTINE meteogram_finalize


  !! IO PE gathers all buffer information from working PEs and copies
  !! the contents to the global meteogram buffer.
  !! Afterwards, all working PEs flush their local buffers.
  !!
  !! For gathered NetCDF output, this is a collective operation.
  !!
  SUBROUTINE meteogram_collect_buffers(mtgrm, jg)
    ! patch buffer
    TYPE(t_buffer_state), INTENT(inout) :: mtgrm
    INTEGER, INTENT(in) :: jg


#ifndef NOMPI
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":meteogram_collect_buffers"
    INTEGER     :: station_idx, position, icurrent,   &
      &            icurrent_recv(max_num_stations), &
      &            istation, nstations, &
      &            iowner
    !> MPI buffer for station data
    CHARACTER, ALLOCATABLE :: msg_buffer(:,:)

    INTEGER :: ierror, req(2+max_num_stations), &
      stati(mpi_status_size, 2+max_num_stations)
    LOGICAL :: is_pure_io_pe
    INTEGER :: max_time_stamps

    IF (dbg_level > 5)  WRITE (*,*) routine, " Enter (collecting PE=", mtgrm%l_is_collecting_pe, ")"

    is_pure_io_pe = use_async_name_list_io .AND. my_process_is_io()
    IF (mtgrm%l_is_collecting_pe .NEQV. mtgrm%l_is_sender) THEN
      max_time_stamps = mtgrm%max_time_stamps
      ALLOCATE(msg_buffer(mtgrm%max_buf_size, &
        MERGE(max_num_stations, 1, mtgrm%l_is_collecting_pe)), &
        stat=ierror)
      IF (ierror /= SUCCESS) THEN
        WRITE (0,'(a,i0)') "jg = ", jg, &
          " : message buffer: mtgrm%max_buf_size = ", &
          & mtgrm%max_buf_size, "; MAX_NUM_STATIONS = ", MAX_NUM_STATIONS, &
          "mtgrm_pack_header_ints  = ", mtgrm_pack_header_ints, &
          "p_real_dp_byte          = ", p_real_dp_byte, &
          "p_int_byte              = ", p_int_byte, &
          "max_time_stamps         = ", max_time_stamps, &
          "MAX_DATE_LEN            = ", MAX_DATE_LEN, &
          "max_station_pack_size   = ", &
          station_pack_size(mtgrm%max_time_stamps, mtgrm%var_info, &
          &                 mtgrm%sfc_var_info, mtgrm%io_collect_comm)
        WRITE (message_text, '(3a)') &
          'ALLOCATE of meteogram message buffer failed (', &
          MERGE('collector', 'sender   ', mtgrm%l_is_collecting_pe), ')'
        CALL finish(routine, message_text)
      END IF
      msg_buffer(:,:) = ''
    END IF

    ! -- RECEIVER CODE --
    RECEIVER : IF (mtgrm%l_is_collecting_pe) THEN
      ! launch MPI message requests for station data on foreign PEs
      nstations = SIZE(mtgrm%out_buf%station_idx)
      IF (is_pure_io_pe) THEN
        CALL p_irecv(mtgrm%istep, mtgrm%io_invar_send_rank, &
          &          tag_mtgrm_msg+max_dom+2*jg-2, comm=mtgrm%io_collect_comm, &
          &          request=req(1))
        CALL p_irecv(mtgrm%zdate, mtgrm%io_invar_send_rank, &
          &          tag_mtgrm_msg+max_dom+2*jg-1, comm=mtgrm%io_collect_comm, &
          &          request=req(2))
      ELSE
        req(1) = mpi_request_null
        req(2) = mpi_request_null
      END IF
      DO istation=1,nstations
        iowner = mtgrm%pstation(istation)
        IF ((is_pure_io_pe .OR. iowner /= p_pe_work) .AND. (iowner >= 0)) THEN
          CALL p_irecv_packed(msg_buffer(:,istation), iowner, &
            &    tag_mtgrm_msg+3*max_dom + (jg-1)*tag_domain_shift + istation, &
            &    mtgrm%max_buf_size, comm=mtgrm%io_collect_comm, &
            &    request=req(2+istation))
        ELSE
          req(2+istation) = mpi_request_null
        END IF
      END DO

      ! wait for messages to arrive:
      IF (dbg_level > 5)  WRITE (*,*) routine, " :: call p_wait"
      CALL p_wait(req(:2+nstations), stati)
      IF (dbg_level > 5)  WRITE (*,*) routine, " :: p_wait call done."

      ! unpack received messages:
      DO istation=1,nstations
        IF (dbg_level > 5) WRITE (*,*) "Receiver side: Station ", istation
        iowner = mtgrm%pstation(istation)
        IF (iowner >= 0) THEN
          IF (iowner /= p_pe_work .OR. is_pure_io_pe) THEN
            CALL unpack_station_sample(mtgrm%var_info, mtgrm%sfc_var_info, &
              mtgrm%out_buf, msg_buffer(:,istation), istation, nstations, &
              icurrent_recv(istation))
          ELSE
            icurrent_recv(istation) = mtgrm%icurrent
          END IF
        ELSE
          icurrent_recv(istation) = -2
          IF (dbg_level > 5) WRITE (*,*) "skipping station!"
        END IF
      END DO

      IF (is_pure_io_pe) THEN
        CALL mpi_get_count(stati(:,1), p_int, icurrent, ierror)
      ELSE
        icurrent = mtgrm%icurrent
      END IF

      ! consistency check
      ! Note: We only check the number of received time stamps, not the
      !       exact sample dates themselves
      IF (ANY(icurrent_recv(1:nstations) /= icurrent &
        &     .AND. mtgrm%pstation(1:nstations) /= -1)) &
        CALL finish(routine, "Received inconsistent time slice data!")
      mtgrm%icurrent = icurrent

    END IF RECEIVER

    ! -- SENDER CODE --
    SENDER : IF (mtgrm%l_is_sender) THEN
      IF (p_pe_work == mtgrm%io_invar_send_rank) THEN
        CALL p_isend(mtgrm%istep, mtgrm%io_collector_rank, &
          &          tag_mtgrm_msg+max_dom+2*jg-2, &
          &          p_count=mtgrm%icurrent, &
          &          comm=mtgrm%io_collect_comm)
        CALL p_isend(mtgrm%zdate, mtgrm%io_collector_rank, &
          &          tag_mtgrm_msg+max_dom+2*jg-1, &
          &          p_count=mtgrm%icurrent, &
          &          comm=mtgrm%io_collect_comm)
      END IF
      ! pack station into buffer; send it
      DO istation=1,mtgrm%nstations_local
        station_idx = mtgrm%global_idx(istation)
        CALL pack_station_sample(msg_buffer(:,1), position, &
          mtgrm%icurrent, &
          mtgrm%out_buf, istation, mtgrm%var_info, &
          mtgrm%global_idx(istation), mtgrm%io_collect_comm)
        ! (blocking) send of packed station data to IO PE:
        CALL p_send_packed(msg_buffer, mtgrm%io_collector_rank, &
          &    tag_mtgrm_msg+3*max_dom + (jg-1)*tag_domain_shift + station_idx,&
          &    position, comm=mtgrm%io_collect_comm)
        IF (dbg_level > 0) &
          WRITE (*,*) "Sending ", icurrent, " time slices, station ", &
          station_idx

      END DO
      IF (p_pe_work == mtgrm%io_invar_send_rank) CALL p_wait()
    END IF SENDER

    IF (dbg_level > 5)  WRITE (*,*) routine, " Leave (collecting PE=", mtgrm%l_is_collecting_pe, ")"

#endif

  END SUBROUTINE meteogram_collect_buffers

  !> unpacks all data transferred for a single station into gathering ranks output buffer
  SUBROUTINE unpack_station_sample(var_info, sfc_var_info, out_buf, &
      sttn_buffer, istation, nstations, icurrent)
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)
    TYPE(t_mtgrm_out_buffer), INTENT(inout) :: out_buf
    CHARACTER, INTENT(in) :: sttn_buffer(:)
    INTEGER, INTENT(in) :: istation, nstations
    INTEGER, INTENT(out) :: icurrent

    INTEGER :: position, ivar, nvars, nvals, station_idx
    CHARACTER(len=*), PARAMETER :: routine = modname//"::unpack_station_sample"

    position = 0

    !-- unpack global time stamp index
    CALL p_unpack_int(sttn_buffer, position, icurrent)
    IF (dbg_level > 0) &
      WRITE (*,'(3(a,i0))') "Receiving ", icurrent, &
      & " time slices from station ", istation, "/", nstations

    !-- unpack station header information
    CALL p_unpack_int(sttn_buffer, position, station_idx)
    IF (out_buf%station_idx(istation) /= station_idx) &
      CALL finish(routine, "non-matching global indices")

    !-- unpack meteogram data:
    nvars = SIZE(var_info)
    DO ivar = 1, nvars
      nvals = var_info(ivar)%nlevs * icurrent
      CALL p_unpack_real_2d(sttn_buffer, position, &
        &   out_buf%atmo_vars(ivar)%a(istation, :, :), nvals)
    END DO
    nvars = SIZE(sfc_var_info)
    nvals = icurrent
    DO ivar = 1, nvars
      CALL p_unpack_real_1d(sttn_buffer, position, &
        &  out_buf%sfc_vars(ivar)%a(istation, :), icurrent)
    END DO
  END SUBROUTINE unpack_station_sample


  !> fills mpi packed buffer for sending to gathering rank
  SUBROUTINE pack_station_sample(sttn_buffer, pos, icurrent, out_buf, istation,&
    var_info, global_idx, io_collect_comm)
    CHARACTER, INTENT(out) :: sttn_buffer(:)
    INTEGER, INTENT(in) :: icurrent, global_idx, istation, io_collect_comm
    INTEGER, INTENT(out) :: pos
    TYPE(t_mtgrm_out_buffer), INTENT(inout) :: out_buf
    TYPE(t_var_info), INTENT(in) :: var_info(:)

    INTEGER :: ivar, nvars, nlevs
    pos = 0

    !-- pack global time stamp index
    CALL p_pack_int(icurrent, sttn_buffer, pos, io_collect_comm)

    !-- pack meteogram header (information on location, ...)
    CALL p_pack_int(global_idx, sttn_buffer, pos, io_collect_comm)

    !-- pack meteogram data:
    nvars = SIZE(out_buf%atmo_vars)
    DO ivar = 1, nvars
      nlevs = var_info(ivar)%nlevs
      CALL p_pack_real_2d(out_buf%atmo_vars(ivar)%a(istation,:,1:icurrent), &
        &                 nlevs*icurrent, sttn_buffer, pos, io_collect_comm)
    END DO
    nvars = SIZE(out_buf%sfc_vars)
    DO ivar = 1, nvars
      CALL p_pack_real_1d(out_buf%sfc_vars(ivar)%a(istation,1:icurrent), &
        &                 icurrent, sttn_buffer, pos, io_collect_comm)
    END DO

  END SUBROUTINE pack_station_sample

  FUNCTION station_pack_size(max_time_stamps, var_info, &
    sfc_var_info, io_collect_comm) RESULT(pack_size)
    INTEGER, INTENT(in) :: max_time_stamps, io_collect_comm
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)
    INTEGER :: pack_size, ivar
    pack_size = p_pack_size_int(1, io_collect_comm) * mtgrm_pack_header_ints
    DO ivar = 1, SIZE(var_info)
      pack_size = pack_size &
        + p_pack_size_real_dp(var_info(ivar)%nlevs * max_time_stamps, &
        &                     io_collect_comm)
    END DO
    pack_size = pack_size + SIZE(sfc_var_info) &
      * p_pack_size_real_dp(max_time_stamps, io_collect_comm)
  END FUNCTION station_pack_size

  !! The IO PE creates and opens a disk file for output.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  SUBROUTINE meteogram_open_file(meteogram_output_config, mtgrm, jg, invariants)
    !> station data from namelist
    TYPE(t_meteogram_output_config), INTENT(IN) :: meteogram_output_config
    !> patch index
    INTEGER,                             INTENT(in) :: jg
    !> patch buffer
    TYPE(t_buffer_state), INTENT(inout) :: mtgrm
    TYPE(t_mtgrm_invariants), INTENT(in) :: invariants

    ! local variables:
    CHARACTER(len=*), PARAMETER :: &
      &  routine = "mo_meteogram_output:meteogram_open_file"
    INTEGER :: old_mode
    LOGICAL :: mtgrm_file_exists

    IF (meteogram_output_config%ftype /= FTYPE_NETCDF) &
      CALL finish(routine, "Output format not yet implemented.")

    ! In "non-distributed" mode, station data is gathered by PE #0
    ! which writes a single file.

    ! skip routine, if this PE has nothing to do...
    IF  (.NOT. mtgrm%l_is_writer) RETURN

    IF (dbg_level > 5)  WRITE (*,*) routine, " Enter"

    ! create a file name for this PE:
    CALL meteogram_create_filename(meteogram_output_config, jg)

    INQUIRE(file=TRIM(mtgrm%meteogram_file_info%zname), &
      exist=mtgrm_file_exists)
    IF (.NOT. mtgrm_file_exists .OR. &
      .NOT. meteogram_output_config%append_if_exists) THEN
      CALL meteogram_create_file(meteogram_output_config, &
        mtgrm%meteogram_file_info, mtgrm%pstation, mtgrm%out_buf, &
        mtgrm%var_info, mtgrm%sfc_var_info, invariants)
    ELSE
      CALL meteogram_append_file(mtgrm%meteogram_file_info, mtgrm%sfc_var_info)
      ! inquire about current number of records in file:
      CALL nf(nf90_inquire_dimension(mtgrm%meteogram_file_info%ncid%file_id, &
        mtgrm%meteogram_file_info%ncid%timeid, len = mtgrm%time_offset), routine)

    END IF
    IF (dbg_level > 5)  WRITE (*,*) routine, " Leave"

  END SUBROUTINE meteogram_open_file

  !> create meteogram netcdf dataset and call create_netcdf_meta for
  !! corresponding variables and meta-data, afterwards add
  !! time-invariant data via put_station_invariants
  SUBROUTINE meteogram_create_file(output_config, file_info, pstation, &
    out_buf, var_info, sfc_var_info, invariants)
    TYPE(t_meteogram_output_config), INTENT(IN) :: output_config
    TYPE(t_meteogram_file), INTENT(inout) :: file_info
    TYPE(t_mtgrm_out_buffer), INTENT(in) :: out_buf
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)
    INTEGER, INTENT(in) :: pstation(:)
    TYPE(t_mtgrm_invariants), INTENT(in) :: invariants

    INTEGER :: iowner, station_idx
    INTEGER :: ncfile, nstations, istation
    CHARACTER(len=*), PARAMETER :: routine = modname//":meteogram_create_file"

    ! create NetCDF file:
    CALL nf(nf90_create(TRIM(file_info%zname), &
      &                 IOR(nf90_clobber, nf90_64bit_offset), &
      &                 file_info%ncid%file_id), routine)
    CALL create_netcdf_meta(file_info%ncid, file_info%cf, &
      &                     file_info%uuid_string, &
      &                     file_info%number_of_grid_used, &
      &                     var_info, sfc_var_info, invariants)

    nstations = SIZE(invariants%fc)
    DO istation=1,nstations
      IF (dbg_level > 5)  WRITE (*,*) "station ", istation

      iowner = pstation(istation)
      IF (iowner >= 0) THEN
        station_idx = out_buf%station_idx(istation)
        CALL put_station_invariants(file_info%ncid, istation, &
          output_config%station_list(station_idx))
      ELSE IF (dbg_level > 5) THEN
        WRITE (*,*) "skipping station!"
      END IF
    END DO
    CALL put_invariants(file_info%ncid, var_info, invariants)
  END SUBROUTINE meteogram_create_file

  !> add meteogram netcdf dataset variables and meta-data, most
  !! importantly records the per-variable id into argument ncid
  SUBROUTINE create_netcdf_meta(ncid, cf, uuid_string, number_of_grid_used, &
    var_info, sfc_var_info, invariants)
    TYPE(t_ncid), INTENT(inout) :: ncid
    TYPE(t_cf_global), INTENT(in) :: cf
    CHARACTER(len=*), INTENT(in) :: uuid_string
    INTEGER, INTENT(in) :: number_of_grid_used
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)
    TYPE(t_mtgrm_invariants), INTENT(in) :: invariants

    INTEGER :: station_name_dims(2), var_name_dims(2), &
      &        time_string_dims(2), tile_dims(2),      &
      &        var_dims(4), sfcvar_dims(3),            &
      &        height_level_dims(3),                   &
      &        tlen, istart(2), icount(2)
    INTEGER :: old_mode, ncfile, nstations, ivar, nvars, nsfcvars

    CHARACTER(len=*), PARAMETER :: routine = modname//'::create_netcdf_meta'

    nvars    = SIZE(var_info)
    nsfcvars = SIZE(sfc_var_info)

    ncfile = ncid%file_id
    CALL nf(nf90_set_fill(ncfile, nf90_nofill, old_mode), routine)
    CALL put_global_txt_att('title', TRIM(cf%title))
    CALL put_global_txt_att('history', TRIM(cf%history))
    CALL put_global_txt_att('institution', TRIM(cf%institution))
    CALL put_global_txt_att('source', TRIM(cf%source))
    CALL put_global_txt_att('comment', TRIM(cf%comment))
    CALL put_global_txt_att('references', TRIM(cf%references))
    CALL put_global_txt_att('uuidOfHGrid', uuid_string)
    CALL nf(nf90_put_att(ncfile, NF90_GLOBAL, 'numberOfGridUsed', &
      &                  number_of_grid_used), routine)


    ! for the definition of a character-string variable define
    ! character-position dimension for strings
    CALL nf(nf90_def_dim(ncfile, "stringlen", MAX_DESCR_LENGTH, &
      &                ncid%charid), routine)
    ! station header:
    nstations = SIZE(invariants%fc)
    CALL nf(nf90_def_dim(ncfile, 'nstations', nstations, &
      &                ncid%nstations), routine)
    ! write variables:
    CALL nf(nf90_def_dim(ncfile, 'nvars', nvars, ncid%nvars), routine)
    CALL nf(nf90_def_dim(ncfile, 'ntiles', invariants%ntiles_mtgrm, &
      &                ncid%ntiles), routine)
    IF (nsfcvars > 0) &
      CALL nf(nf90_def_dim(ncfile, 'nsfcvars', nsfcvars, &
      &                  ncid%nsfcvars), routine)
    CALL nf(nf90_def_dim(ncfile, 'max_nlevs',  invariants%max_nlevs, &
      &     ncid%max_nlevs), routine)
    ! create time dimension:
    CALL nf(nf90_def_dim(ncfile, 'time', NF90_UNLIMITED, ncid%timeid), routine)

    ! create station variables:
    station_name_dims(1) = ncid%charid
    station_name_dims(2) = ncid%nstations
    CALL nf(nf90_def_var(ncfile, "station_name", NF90_CHAR, station_name_dims, &
      &                  ncid%station_name), routine)
    CALL nf_add_descr("Station name (character string)", ncfile, &
      &               ncid%station_name)
    CALL nf(nf90_def_var(ncfile, "station_lon", NF90_DOUBLE, &
      &                  ncid%nstations, ncid%station_lon), routine)
    CALL nf_add_descr("Longitude of meteogram station", ncfile, &
      &               ncid%station_lon)
    CALL nf(nf90_def_var(ncfile, "station_lat", NF90_DOUBLE, &
      &                  ncid%nstations, ncid%station_lat), routine)
    CALL nf_add_descr("Latitude of meteogram station", ncfile, ncid%station_lat)
    CALL nf(nf90_def_var(ncfile, "station_idx", NF90_INT, &
      &                  ncid%nstations, ncid%station_idx), routine)
    CALL nf_add_descr("Global triangle adjacent to meteogram station (index)", &
      &               ncfile, ncid%station_idx)
    CALL nf(nf90_def_var(ncfile, "station_blk", NF90_INT, &
      &                  ncid%nstations, ncid%station_blk), routine)
    CALL nf_add_descr("Global triangle adjacent to meteogram station (block)", &
      &               ncfile, ncid%station_blk)
    CALL nf(nf90_def_var(ncfile, "station_hsurf", NF90_DOUBLE, &
      &                  ncid%nstations, ncid%station_hsurf), routine)
    CALL nf_add_descr("Meteogram station surface height", ncfile, &
      &               ncid%station_hsurf)
    CALL nf(nf90_def_var(ncfile, "station_frland", NF90_DOUBLE, &
      &                  ncid%nstations, ncid%station_frland), routine)
    CALL nf_add_descr("Meteogram station land fraction", ncfile, &
      &               ncid%station_frland)
    CALL nf(nf90_def_var(ncfile, "station_fc", NF90_DOUBLE, &
      &                  ncid%nstations, ncid%station_fc), routine)
    CALL nf_add_descr("Meteogram station Coriolis parameter", ncfile, &
      &               ncid%station_fc)
    CALL nf(nf90_def_var(ncfile, "station_soiltype", NF90_INT, &
      &                  ncid%nstations, ncid%station_soiltype), routine)
    CALL nf_add_descr("Meteogram station soil type", ncfile, &
      &     ncid%station_soiltype)

    tile_dims(1) = ncid%ntiles
    tile_dims(2) = ncid%nstations
    CALL nf(nf90_def_var(ncfile, "station_tile_frac", NF90_DOUBLE, tile_dims, &
      &                  ncid%station_tile_frac), routine)
    CALL nf_add_descr("Meteogram station tile fractions", ncfile, &
      &               ncid%station_tile_frac)
    CALL nf(nf90_def_var(ncfile, "station_tile_luclass", NF90_INT, tile_dims, &
      &                  ncid%station_tile_luclass), routine)
    CALL nf_add_descr("Meteogram station tile specific land-use classes", &
      &               ncfile, ncid%station_tile_luclass)


    ! create variable info fields:
    ! volume variables
    var_name_dims(1) = ncid%charid
    var_name_dims(2) = ncid%nvars
    CALL nf(nf90_def_var(ncfile, "var_name", NF90_CHAR, var_name_dims, &
      &                  ncid%var_name), routine)
    CALL nf_add_descr("Variable name (character string)", ncfile, ncid%var_name)
    CALL nf(nf90_def_var(ncfile, "var_long_name", NF90_CHAR, var_name_dims, &
      &                  ncid%var_longname), routine)
    CALL nf_add_descr("Variable name (long, character string)", ncfile, &
      &               ncid%var_longname)
    CALL nf(nf90_def_var(ncfile, "var_unit", NF90_CHAR, var_name_dims, &
      &                  ncid%var_unit), routine)
    CALL nf_add_descr("Variable unit (character string)", ncfile, ncid%var_unit)
    CALL nf(nf90_def_var(ncfile, "var_group_id", NF90_INT, &
      &                  ncid%nvars, ncid%var_group_id), routine)
    CALL nf_add_descr("Variable group ID", ncfile, ncid%var_group_id)
    CALL nf(nf90_def_var(ncfile, "var_nlevs", NF90_INT, &
      &                  ncid%nvars, ncid%var_nlevs), routine)
    CALL nf_add_descr("No. of levels for volume variable", ncfile, &
      &               ncid%var_nlevs)
    ! surface variables:
    IF (nsfcvars > 0) THEN
      var_name_dims(1) = ncid%charid
      var_name_dims(2) = ncid%nsfcvars
      CALL nf(nf90_def_var(ncfile, "sfcvar_name", NF90_CHAR, var_name_dims, &
        &                  ncid%sfcvar_name), routine)
      CALL nf_add_descr("Surface variable name (character string)", ncfile, &
        &               ncid%sfcvar_name)
      CALL nf(nf90_def_var(ncfile, "sfcvar_long_name", NF90_CHAR, var_name_dims, &
        &                  ncid%sfcvar_longname), routine)
      CALL nf_add_descr("Surface variable name (long, character string)", &
        &               ncfile, ncid%sfcvar_longname)
      CALL nf(nf90_def_var(ncfile, "sfcvar_unit", NF90_CHAR, var_name_dims, &
        &                  ncid%sfcvar_unit), routine)
      CALL nf_add_descr("Surface variable unit (character string)", ncfile, &
        &               ncid%sfcvar_unit)
      CALL nf(nf90_def_var(ncfile, "sfcvar_group_id", NF90_INT, &
        &                  ncid%nsfcvars, ncid%sfcvar_group_id), routine)
      CALL nf_add_descr("Surface variable group ID", ncfile, &
      &     ncid%sfcvar_group_id)
    END IF

    ! create variables for time slice info:
    CALL nf(nf90_def_var(ncfile, "time_step", NF90_INT, &
      &                  ncid%timeid, ncid%time_step), routine)
    CALL nf_add_descr("Time step indices", ncfile, &
      &     ncid%time_step)
    time_string_dims(1) = ncid%charid
    time_string_dims(2) = ncid%timeid
    CALL nf(nf90_def_var(ncfile, "date", NF90_CHAR, time_string_dims, &
      &                  ncid%dateid), routine)
    CALL nf_add_descr("Sample dates (character string)", ncfile, &
      &     ncid%dateid)

    ! height levels
    height_level_dims(1) = ncid%nstations
    height_level_dims(2) = ncid%nvars
    height_level_dims(3) = ncid%max_nlevs
    CALL nf(nf90_def_var(ncfile, "heights", NF90_DOUBLE, height_level_dims, &
      &                  ncid%var_heights), routine)
    CALL nf_add_descr("level heights for volume variables", ncfile, &
      &               ncid%var_heights)

    ! add value buffer for volume variables:
    var_dims(1) = ncid%nstations
    var_dims(2) = ncid%nvars
    var_dims(3) = ncid%max_nlevs
    var_dims(4) = ncid%timeid
    CALL nf(nf90_def_var(ncfile, "values", NF90_DOUBLE, var_dims, &
      &                  ncid%var_values), routine)
    CALL nf_add_descr("value buffer for volume variables", ncfile, &
      &               ncid%var_values)
    ! add value buffer for surface variables:
    IF (nsfcvars > 0) THEN
      sfcvar_dims(1) = ncid%nstations
      sfcvar_dims(2) = ncid%nsfcvars
      sfcvar_dims(3) = ncid%timeid
      CALL nf(nf90_def_var(ncfile, "sfcvalues", NF90_DOUBLE, sfcvar_dims, &
        &                  ncid%sfcvar_values), routine)
      CALL nf_add_descr("value buffer for surface variables", ncfile, &
      &     ncid%sfcvar_values)
    END IF

    ! ----------------------
    ! End of definition mode
    CALL nf(nf90_enddef(ncfile), routine)
    IF (dbg_level > 7)  WRITE (*,*) routine, " : End of definition mode"

    istart(1) = 1
    icount(2) = 1
    DO ivar=1,nvars
      istart(2) = ivar
      tlen = LEN_TRIM(var_info(ivar)%cf%standard_name)
      icount(1) = tlen
      CALL nf(nf90_put_var(ncfile, ncid%var_name, &
        &                  var_info(ivar)%cf%standard_name(1:tlen), &
        &                  istart, icount), &
        &     routine)
      tlen = LEN_TRIM(var_info(ivar)%cf%long_name)
      icount(1) = tlen
      CALL nf(nf90_put_var(ncfile, ncid%var_longname, &
        &                  var_info(ivar)%cf%long_name(1:tlen), &
        &                  istart, icount), &
        &     routine)
      tlen = LEN_TRIM(var_info(ivar)%cf%units)
      icount(1) = tlen
      CALL nf(nf90_put_var(ncfile, ncid%var_unit, &
        &                  var_info(ivar)%cf%units(1:tlen), &
        &                  istart, icount), &
        &     routine)
      CALL nf(nf90_put_var(ncfile, ncid%var_group_id, &
        &     var_info(ivar)%igroup_id, [ivar]), routine)
      CALL nf(nf90_put_var(ncfile, ncid%var_nlevs, &
        &     var_info(ivar)%nlevs, [ivar]), routine)
    END DO

    DO ivar=1,nsfcvars
      istart(2) = ivar
      tlen = LEN_TRIM(sfc_var_info(ivar)%cf%standard_name)
      icount(1) = tlen
      CALL nf(nf90_put_var(ncfile, ncid%sfcvar_name, &
        &                  sfc_var_info(ivar)%cf%standard_name(1:tlen), &
        &                  istart, icount), &
        &     routine)
      tlen = LEN_TRIM(sfc_var_info(ivar)%cf%long_name)
      icount(1) = tlen
      CALL nf(nf90_put_var(ncfile, ncid%sfcvar_longname, &
        &                  sfc_var_info(ivar)%cf%long_name(1:tlen), &
        &                  istart, icount), &
        &     routine)
      tlen = LEN_TRIM(sfc_var_info(ivar)%cf%units)
      icount(1) = tlen
      CALL nf(nf90_put_var(ncfile, ncid%sfcvar_unit, &
        &                  sfc_var_info(ivar)%cf%units(1:tlen), &
        &                  istart, icount), &
        &     routine)
      CALL nf(nf90_put_var(ncfile, ncid%sfcvar_group_id, &
        &                  sfc_var_info(ivar)%igroup_id, [ivar]), &
        &     routine)
    END DO

  CONTAINS
    SUBROUTINE put_global_txt_att(attname, attval)
      CHARACTER(len=*), INTENT(in) :: attname, attval
      CHARACTER(len=*), PARAMETER :: routine = modname//":put_global_txt_att"
      CALL nf(nf90_put_att(ncfile, NF90_GLOBAL, attname, attval), &
        &     routine)
    END SUBROUTINE put_global_txt_att
  END SUBROUTINE create_netcdf_meta

  !> write time-invariant meteogram data for single station to netcdf
  !! dataset
  SUBROUTINE put_station_invariants(ncid, istation, station_cfg)
    TYPE(t_ncid), INTENT(in) :: ncid
    INTEGER, INTENT(in) :: istation
    TYPE(t_station_list), INTENT(in) :: station_cfg

    INTEGER :: tlen, ncfile
    INTEGER :: istart(2), icount(2)
    CHARACTER(len=*), PARAMETER :: routine = modname//"::put_station_invariants"

    ncfile = ncid%file_id
    tlen = LEN_TRIM(station_cfg%zname)
    istart(1) = 1
    istart(2) = istation
    icount(1) = tlen
    icount(2) = 1
    CALL nf(nf90_put_var(ncfile, ncid%station_name, &
      &                  station_cfg%zname(1:tlen), &
      &                  istart(1:2), icount(1:2)), &
      &     routine)
    CALL nf(nf90_put_var(ncfile, ncid%station_lon, station_cfg%location%lon, &
      &                  [istation]), routine)
    CALL nf(nf90_put_var(ncfile, ncid%station_lat, station_cfg%location%lat, &
      &                  [istation]), routine)
  END SUBROUTINE put_station_invariants

  !> write time-invariant and station-independent meteogram data to
  !! netcdf dataset
  SUBROUTINE put_invariants(ncid, var_info, invariants)
    TYPE(t_ncid), INTENT(in) :: ncid
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_mtgrm_invariants), INTENT(in) :: invariants

    INTEGER :: ivar, nvars
    INTEGER :: istart(3), icount(3)
    CHARACTER(len=*), PARAMETER :: routine = modname//"::put_invariants"

    istart = 1
    nvars = SIZE(var_info)
    ! nstations
    icount(1) = SIZE(invariants%fc)
    ! model level heights
    icount(2) = 1
    DO ivar=1,nvars
      istart(2) = ivar
      icount(3) = var_info(ivar)%nlevs
      CALL nf(nf90_put_var(ncid%file_id, ncid%var_heights, &
        &                  invariants%heights(ivar)%a, &
        &     istart, icount), &
        &     routine)
    END DO
    istart(2) = 1
    icount(2) = icount(1)
    icount(1) = invariants%ntiles_mtgrm
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_tile_frac, &
      &                  invariants%tile_frac, &
      &                  istart(1:2), icount(1:2)),&
      &     routine)
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_tile_luclass, &
      &                  invariants%tile_luclass, &
      &                  istart(1:2), icount(1:2)), &
      &     routine)
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_soiltype, &
      &                  invariants%soiltype, &
      &                  istart(2:2), icount(2:2)), &
      &     routine)
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_hsurf, &
      &                  invariants%hsurf, &
      &                  istart(2:2), icount(2:2)), &
      &     routine)
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_frland, &
      &                  invariants%frland, &
      &                  istart(2:2), icount(2:2)), &
      &    routine)
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_fc, &
      &                  invariants%fc, &
      &                  istart(2:2), icount(2:2)), &
      &     routine)
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_idx, &
      &                  invariants%tri_idx(:,1), &
      &                  istart(2:2), icount(2:2)), &
      &     routine)
    CALL nf(nf90_put_var(ncid%file_id, ncid%station_blk, &
      &                  invariants%tri_idx(:,2), &
      &                  istart(2:2), icount(2:2)), &
      &     routine)
  END SUBROUTINE put_invariants

  !> instead of creating a dataset, open an existing meteogram file
  !! and query variable IDs for appending.
  SUBROUTINE meteogram_append_file(file_info, sfc_var_info)
    TYPE(t_meteogram_file), INTENT(inout) :: file_info
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//":meteogram_append_file"
    INTEGER :: old_mode, ncfile, nsfcvars

    nsfcvars = SIZE(sfc_var_info)
    CALL nf(nf90_open(TRIM(file_info%zname), NF90_WRITE, &
      &               file_info%ncid%file_id), routine)
    ncfile = file_info%ncid%file_id
    CALL nf(nf90_set_fill(ncfile, nf90_nofill, old_mode), routine)
    CALL nf(nf90_inq_dimid(ncfile, "stringlen", file_info%ncid%charid), &
      &     routine)
    CALL nf(nf90_inq_dimid(ncfile, 'nstations', file_info%ncid%nstations), &
      &     routine)
    CALL nf(nf90_inq_dimid(ncfile, 'nvars', file_info%ncid%nvars), &
      &     routine)
    IF (nsfcvars > 0) &
      CALL nf(nf90_inq_dimid(ncfile, 'nsfcvars', file_info%ncid%nsfcvars), &
      &       routine)
    CALL nf(nf90_inq_dimid(ncfile, 'max_nlevs', file_info%ncid%max_nlevs), &
      &     routine)
    CALL nf(nf90_inq_dimid(ncfile, 'time', file_info%ncid%timeid), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "station_name", file_info%ncid%station_name), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "station_lon", file_info%ncid%station_lon), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "station_lat", file_info%ncid%station_lat), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "station_idx", file_info%ncid%station_idx), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "station_blk", file_info%ncid%station_blk), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "station_hsurf", file_info%ncid%station_hsurf), &
      &      routine)
    CALL nf(nf90_inq_varid(ncfile, "station_frland", file_info%ncid%station_frland), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "station_fc", file_info%ncid%station_fc), &
      &      routine)
    CALL nf(nf90_inq_varid(ncfile, "station_soiltype", file_info%ncid%station_soiltype), &
      &     routine)

    ! inquire variable info fields:
    ! volume variables
    CALL nf(nf90_inq_varid(ncfile, "var_name", file_info%ncid%var_name), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "var_long_name", file_info%ncid%var_longname), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "var_unit", file_info%ncid%var_unit), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "var_group_id", file_info%ncid%var_group_id), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "var_nlevs", file_info%ncid%var_nlevs), &
      &     routine)
    ! surface variables:
    IF (nsfcvars > 0) THEN
      CALL nf(nf90_inq_varid(ncfile, "sfcvar_name", file_info%ncid%sfcvar_name), &
        &     routine)
      CALL nf(nf90_inq_varid(ncfile, "sfcvar_long_name", file_info%ncid%sfcvar_longname), &
        &     routine)
      CALL nf(nf90_inq_varid(ncfile, "sfcvar_unit", file_info%ncid%sfcvar_unit), &
        &     routine)
      CALL nf(nf90_inq_varid(ncfile, "sfcvar_group_id", file_info%ncid%sfcvar_group_id), &
        &     routine)
    END IF

    ! create variables for time slice info:
    CALL nf(nf90_inq_varid(ncfile, "time_step", file_info%ncid%time_step), &
      &     routine)
    CALL nf(nf90_inq_varid(ncfile, "date", file_info%ncid%dateid), &
      &     routine)

    ! height levels
    CALL nf(nf90_inq_varid(ncfile, "heights", file_info%ncid%var_heights), &
      &     routine)

    ! add value buffer for volume variables:
    CALL nf(nf90_inq_varid(ncfile, "values", file_info%ncid%var_values), &
      &     routine)
    ! add value buffer for surface variables:
    IF (nsfcvars > 0) THEN
      CALL nf(nf90_inq_varid(ncfile, "sfcvalues", file_info%ncid%sfcvar_values), &
        &     routine)
    END IF

  END SUBROUTINE meteogram_append_file

  !! The IO PE writes the global meteogram buffer to the output
  !! file. Afterwards, the global meteogram buffer is cleared.
  !!
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  SUBROUTINE meteogram_flush_file(jg, lacc)
    INTEGER, INTENT(IN)         :: jg       !< patch index
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
    ! local variables:
    INTEGER :: ivar
    CHARACTER(len=*), PARAMETER :: routine = modname//":meteogram_flush_file"
    LOGICAL :: lzacc ! non-optional version of lacc

    IF (dbg_level > 5)  WRITE (*,*) routine, " Enter"

    CALL set_acc_host_or_device(lzacc, lacc)
    DO ivar=1,SIZE(mtgrm(jg)%var_info)
      !$ACC UPDATE HOST(mtgrm(jg)%out_buf%atmo_vars(ivar)%a) ASYNC(1) IF(lzacc)
    END DO
    DO ivar=1,SIZE(mtgrm(jg)%sfc_var_info)
      !$ACC UPDATE HOST(mtgrm(jg)%out_buf%sfc_vars(ivar)%a) ASYNC(1) IF(lzacc)
    END DO
    !$ACC WAIT(1)

    ! In "non-distributed" mode, station data is gathered by PE #0
    ! which writes a single file:
    IF (.NOT. mtgrm(jg)%meteogram_file_info%ldistributed) THEN
      CALL meteogram_collect_buffers(mtgrm(jg), jg)
    ELSE
#ifdef _OPENACC
      CALL finish(routine, "ldistributed has not been tested with OpenACC.")
      ! MJ assumes that the code should work. Please test it and if it works,
      ! please add a ldistributed test to the testing and remove this warning.
#endif
    ENDIF
    IF (mtgrm(jg)%l_is_writer) THEN
      CALL disk_flush(mtgrm(jg)%var_info, mtgrm(jg)%sfc_var_info, &
      &               mtgrm(jg)%out_buf, mtgrm(jg)%time_offset, &
      &               mtgrm(jg)%icurrent, mtgrm(jg)%istep, mtgrm(jg)%zdate, &
      &               mtgrm(jg)%meteogram_file_info%ncid)
      mtgrm(jg)%time_offset = mtgrm(jg)%time_offset + mtgrm(jg)%icurrent
      mtgrm(jg)%icurrent = 0
    END IF

    ! finally, reset buffer counter for new data
    mtgrm(jg)%icurrent = 0

    IF (dbg_level > 5)  WRITE (*,*) routine, " Leave"

  END SUBROUTINE meteogram_flush_file

  !> writer part of intermediate disk flushing
  SUBROUTINE disk_flush(var_info, sfc_var_info, &
    out_buf, time_offset, icurrent, istep, zdate, ncid)
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_sfc_var_info), INTENT(in) :: sfc_var_info(:)
    TYPE(t_mtgrm_out_buffer), INTENT(in) :: out_buf
    INTEGER, INTENT(in) :: icurrent, istep(:), time_offset
    CHARACTER(len=MAX_DATE_LEN), INTENT(in) :: zdate(:)

    TYPE(t_ncid), INTENT(in) :: ncid

    INTEGER :: itime, nstations, ivar, nlevs, nvars, nsfcvars
    INTEGER :: istart(4), icount(4), ncfile
    CHARACTER(len=*), PARAMETER :: routine = modname//"::disk_flush"

    nvars    = SIZE(var_info)
    nsfcvars = SIZE(sfc_var_info)
    ncfile   = ncid%file_id

    IF (dbg_level > 0) THEN
      WRITE(message_text,*) "Meteogram"
      CALL message(routine, TRIM(message_text))
    END IF

    IF (dbg_level > 0) &
      WRITE (*,'(a,i0,a)') "Writing ", icurrent, " time slices to disk."

    ! write time stamp info:
    istart(1) = time_offset+1
    icount(1) = icurrent
    CALL nf(nf90_put_var(ncfile, ncid%time_step, istep, &
      &                  istart(1:1), icount(1:1)), &
      &     routine)
    ! volume variables:
    IF (nvars > 0) THEN
      nstations = SIZE(out_buf%atmo_vars(1)%a, 1)
      istart(1) = 1
      icount(1) = nstations
      icount(2) = 1 ! 1 variable at a time
      istart(3) = 1 ! always start at ilev=1
      istart(4) = time_offset+1
      icount(4) = icurrent
      DO ivar=1,nvars
        nlevs = var_info(ivar)%nlevs
        istart(2) = ivar
        icount(3) = nlevs
        CALL nf(nf90_put_var(ncfile, ncid%var_values,   &
          &                  out_buf%atmo_vars(ivar)%a, &
          &                  istart, icount), &
          &     routine)
      END DO
    END IF
    ! surface variables:
    IF (nsfcvars > 0) THEN
      nstations = SIZE(out_buf%sfc_vars(1)%a, 1)
      istart(1) = 1
      icount(1) = nstations
      icount(2) = 1 ! 1 variable at a time
      istart(3) = time_offset+1
      icount(3) = icurrent
      DO ivar=1,nsfcvars
        istart(2) = ivar
        CALL nf(nf90_put_var(ncfile, ncid%sfcvar_values, &
          &                  out_buf%sfc_vars(ivar)%a, &
          &                  istart(1:3), icount(1:3)), &
          &     routine)
      END DO
    END IF

    icount = 1
    DO itime=1,icurrent
      istart(3) = 1
      icount(3) = LEN_TRIM(zdate(itime))
      istart(4) = time_offset+itime
      icount(4) = 1
      CALL nf(nf90_put_var(ncfile, ncid%dateid, &
        &                  zdate(itime), &
        &                  istart(3:4), icount(3:4)), &
        &     routine)
    END DO
    CALL nf(nf90_sync(ncfile), routine)

  END SUBROUTINE disk_flush

  !! The IO PE closes the meteogram output file.
  !! For gathered NetCDF output, this is a collective operation,
  !! otherwise this is a local operation for the IO PE.
  !!
  SUBROUTINE meteogram_close_file(jg, lacc)
    INTEGER, INTENT(IN)  :: jg    !< patch index
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    CHARACTER(len=*), PARAMETER :: routine=modname//"::meteogram_close_file"

    ! write remaining buffers:
    CALL meteogram_flush_file(jg, lacc=lacc)

    ! Close NetCDF file
    ! skip routine, if this PE has nothing to do...
    IF (mtgrm(jg)%l_is_writer) THEN
      CALL nf(nf90_close(mtgrm(jg)%meteogram_file_info%ncid%file_id), routine)
    END IF
  END SUBROUTINE meteogram_close_file

  !! @return file name of meteogram file.
  !! This is a local operation.
  !!
  SUBROUTINE meteogram_create_filename (meteogram_output_config, jg)

    ! station data from namelist
    TYPE(t_meteogram_output_config), INTENT(IN) :: meteogram_output_config
    ! patch index
    INTEGER, INTENT(IN) :: jg
    ! Local variables
    INTEGER :: my_id, dist_prefix_len
    CHARACTER(len=3+10) :: dist_prefix

#ifndef NOMPI
    IF (meteogram_output_config%ldistributed) THEN
      my_id = get_my_mpi_all_id()
      WRITE(dist_prefix, '(a,i3.3,a)') "PE", my_id, "_"
      dist_prefix_len = LEN_TRIM(dist_prefix)
    ELSE
#endif
      dist_prefix = ''
      dist_prefix_len = 0
#ifndef NOMPI
    END IF
#endif

    SELECT CASE (meteogram_output_config%ftype)
    CASE (FTYPE_NETCDF)
      WRITE (mtgrm(jg)%meteogram_file_info%zname,'(3a,i3.3,a)') &
        TRIM(meteogram_output_config%zprefix), &
        dist_prefix(1:dist_prefix_len), "patch", jg, ".nc"
    END SELECT
  END SUBROUTINE meteogram_create_filename



  !>
  !!  Help functions for NetCDF I/O
  !!
  !!  Adds a string attribute containing variable description.
  SUBROUTINE nf_add_descr(description_str, ncfile, var_id)
    CHARACTER(LEN=*), INTENT(in) :: description_str
    INTEGER         , INTENT(in) :: ncfile, var_id

    CHARACTER(LEN=*), PARAMETER  :: descr_label = "description"
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//":nf_add_descr"
    INTEGER :: desc_tlen

    desc_tlen = LEN_TRIM(description_str)
    CALL nf(nf90_put_att(ncfile, var_id, descr_label, &
      &     description_str(1:desc_tlen)), routine)

  END SUBROUTINE nf_add_descr

  SUBROUTINE set_cf_info(cf, zname, zlong_name, zunit)
    TYPE(t_cf_var), INTENT(inout) :: cf              !< variable name, unit
    CHARACTER(LEN=*), INTENT(in) :: zname, zunit, zlong_name

    cf%standard_name = zname
    cf%long_name = zlong_name
    cf%units = zunit
  END SUBROUTINE set_cf_info

  !>
  !!  Utility function (3d formulation of generic interface).
  !!
  !!  Registers a new atmospheric (volume) variable.
  SUBROUTINE add_atmo_var_3d(meteogram_config, var_list, igroup_id, &
       zname, zunit, zlong_name, var_info, source)
    TYPE(t_meteogram_output_config), INTENT(IN) :: meteogram_config
    TYPE(t_var), INTENT(inout) :: var_list
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    INTEGER,           INTENT(IN)    :: igroup_id
    TYPE(t_var_info), INTENT(inout) :: var_info(:)
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:,:)   !< source array
    ! Local variables
    CHARACTER(*), PARAMETER :: routine = modname//":add_atmo_var_3d"
    INTEGER                         :: nlev, ivar

    IF (dbg_level > 0) &
      &  WRITE ( * , * ) " add_atmo_var_3d: calling for ", TRIM(zname)

    IF (LEN_TRIM(meteogram_config%var_list(1)) /= 0) THEN
      ! If the user has specified a list of variable names to be
      ! included in the meteogram, check if this variable is contained
      ! in the list:
      IF (one_of(zname, meteogram_config%var_list) == -1) RETURN
    END IF

    IF (dbg_level > 0) &
      CALL message(routine, "add atmo var "//zname)

    nlev = SIZE(source, 2) ! get level no from array dimensions

    ! create new variable index
    ivar = var_list%no_atmo_vars + 1
    var_list%no_atmo_vars = ivar

    ! create meteogram data structure
    CALL set_cf_info(var_info(ivar)%cf, zname, zlong_name, zunit)
    var_info(ivar)%igroup_id        = igroup_id
    var_info(ivar)%nlevs            = nlev
    var_info(ivar)%p_source => source

    IF (.NOT. ASSOCIATED(var_info(ivar)%p_source)) THEN
      WRITE (message_text, '(3a)') 'Source array ', &
        TRIM(var_info(ivar)%cf%standard_name), ' not associated!'
      CALL finish(routine, message_text)
    END IF

#ifdef _OPENACC
    IF (.NOT. acc_is_present(var_info(ivar)%p_source)) THEN
      WRITE (message_text, '(3a)') 'Source array ', &
        TRIM(var_info(ivar)%cf%standard_name), ' not present on accelerator!'
      CALL finish(routine, message_text)
    END IF
    !$ACC ENTER DATA ATTACH(var_info(ivar)%p_source)
#endif

  END SUBROUTINE add_atmo_var_3d


  !>
  !!  Utility function (4d formulation of generic interface).
  !!  Adds the 4d var as separate 3d var slices.
  !!
  !!  Registers a new atmospheric (volume) variable.
  SUBROUTINE add_atmo_var_4d(meteogram_config, var_list, igroup_id, &
       zname, zunit, zlong_name, var_info, source, iidx)
    TYPE(t_meteogram_output_config), INTENT(IN) :: meteogram_config
    TYPE(t_var), INTENT(inout) :: var_list
    INTEGER,           INTENT(IN)    :: igroup_id
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    TYPE(t_var_info), INTENT(inout) :: var_info(:)
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:,:,:)   !< source array
    INTEGER,           INTENT(IN), OPTIONAL :: iidx
    ! Local variables
    INTEGER                          :: isource_idx, nidx

    IF (dbg_level > 0) &
      &  WRITE ( * , * ) " add_atmo_var_4d: calling for ", TRIM(zname)

    IF (PRESENT(iidx)) THEN
      CALL add_atmo_var_3d(meteogram_config, var_list, igroup_id, zname, &
        &                  zunit, zlong_name, var_info, source(:,:,:,iidx))
    ELSE
      nidx = SIZE(source, 4) ! get number of 3d var indices (e.g. tile number)
      DO isource_idx=1,nidx
        CALL add_atmo_var_3d(meteogram_config, var_list, igroup_id, &
          &                  zname//"_"//int2string(isource_idx), &
          &                  zunit, zlong_name, var_info, &
          &                  source(:,:,:,isource_idx))
      END DO
    END IF
  END SUBROUTINE add_atmo_var_4d


  !>
  !!  Utility function.
  !!
  !!  Registers a new surface variable.
  SUBROUTINE add_sfc_var_2d(meteogram_config, var_list, igroup_id, &
       zname, zunit, zlong_name, sfc_var_info, source)
    TYPE(t_meteogram_output_config), INTENT(IN) :: meteogram_config
    TYPE(t_var), INTENT(inout) :: var_list
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    INTEGER,           INTENT(IN)    :: igroup_id
    TYPE(t_sfc_var_info), INTENT(inout) :: sfc_var_info(:)
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:)   !< source array
    ! Local variables
    CHARACTER(*), PARAMETER :: routine = modname//":add_sfc_var_2d"
    INTEGER                         :: ivar

    IF (dbg_level > 0) &
      &  WRITE ( * , * ) " add_sfc_var_2d: calling for ", TRIM(zname)

    IF (LEN_TRIM(meteogram_config%var_list(1)) /= 0) THEN
      ! If the user has specified a list of variable names to be
      ! included in the meteogram, check if this variable is contained
      ! in the list:
      IF (one_of(TRIM(zname), meteogram_config%var_list) == -1) RETURN
    END IF

    IF (dbg_level > 0) &
      CALL message(routine, "add surface var "//zname)

    ! create new variable index
    ivar = var_list%no_sfc_vars + 1
    var_list%no_sfc_vars = ivar
    ! create meteogram data structure
    CALL set_cf_info(sfc_var_info(ivar)%cf, zname, zlong_name, zunit)
    sfc_var_info(ivar)%igroup_id        = igroup_id
    sfc_var_info(ivar)%p_source  => source

#ifdef _OPENACC
    IF (.NOT. acc_is_present(sfc_var_info(ivar)%p_source)) THEN
      WRITE (message_text, '(3a)') 'Source array ', &
        TRIM(sfc_var_info(ivar)%cf%standard_name), ' not present on accelerator!'
      CALL finish(routine, message_text)
    END IF
    !$ACC ENTER DATA ATTACH(sfc_var_info(ivar)%p_source)
#endif
  END SUBROUTINE add_sfc_var_2d


  !>
  !!  Utility function.
  !!  Adds the 3d var as separate 2d var slices.
  !!
  !!  Registers a new surface variable.
  SUBROUTINE add_sfc_var_3d(meteogram_config, var_list, igroup_id,  &
      &                     zname, zunit, zlong_name, sfc_var_info, &
      &                     source, iidx)
    TYPE(t_meteogram_output_config), INTENT(IN) :: meteogram_config
    TYPE(t_var), INTENT(inout) :: var_list
    CHARACTER(LEN=*),  INTENT(IN)    :: zname, zunit, zlong_name
    INTEGER,           INTENT(IN)    :: igroup_id
    TYPE(t_sfc_var_info), INTENT(inout) :: sfc_var_info(:)
    REAL(wp), TARGET,  INTENT(IN)    :: source(:,:,:)   !< source array
    INTEGER,           INTENT(IN), OPTIONAL :: iidx     !< third dim of variable (e.g. tile-nr)
    ! Local variables
    INTEGER                 :: isource_idx, nidx

    IF (dbg_level > 0) &
      &  WRITE ( * , * ) " add_sfc_var_3d: calling for ", TRIM(zname)

    IF (PRESENT(iidx)) THEN
      CALL add_sfc_var_2d(meteogram_config, var_list, igroup_id, zname, &
      &                  zunit, zlong_name, sfc_var_info, source(:,:,iidx))
    ELSE 
      nidx = SIZE(source, 3) ! get number of 2d var indices (e.g. tile number)
      DO isource_idx=1,nidx
        CALL add_sfc_var_2d(meteogram_config, var_list, igroup_id, &
          &                 zname//"_"//int2string(isource_idx), &
          &                 zunit, zlong_name, sfc_var_info, &
          &                 source(:,:,isource_idx))
      END DO
    ENDIF
  END SUBROUTINE add_sfc_var_3d

  !> after adding the variables to maximally sized buffers, shrink
  !! wrap the buffers to match the number of variables actually
  !! defined
  SUBROUTINE resize_var_info(var_list, var_info, sfc_var_info)
    TYPE(t_var), INTENT(in) :: var_list
    TYPE(t_var_info), ALLOCATABLE, INTENT(inout) :: var_info(:)
    TYPE(t_sfc_var_info), ALLOCATABLE, INTENT(inout) :: sfc_var_info(:)

    TYPE(t_var_info), ALLOCATABLE :: tmp_var_info(:)
    TYPE(t_sfc_var_info), ALLOCATABLE ::tmp_sfc_var_info(:)
    INTEGER :: nvars

    nvars = var_list%no_atmo_vars
    ALLOCATE(tmp_var_info(nvars))
    tmp_var_info = var_info(1:nvars)
    CALL MOVE_ALLOC(tmp_var_info, var_info)

    nvars = var_list%no_sfc_vars
    ALLOCATE(tmp_sfc_var_info(nvars))
    tmp_sfc_var_info = sfc_var_info(1:nvars)
    CALL MOVE_ALLOC(tmp_sfc_var_info, sfc_var_info)
  END SUBROUTINE resize_var_info

  !> copy var list and buffers on OpenACC device and attach pointers
  SUBROUTINE acc_copyin_var_info(var_info, sfc_var_info)
    TYPE(t_var_info), ALLOCATABLE, INTENT(inout) :: var_info(:)
    TYPE(t_sfc_var_info), ALLOCATABLE, INTENT(inout) :: sfc_var_info(:)

    INTEGER :: ivar

    !$ACC ENTER DATA COPYIN(var_info, sfc_var_info) ASYNC(1)
    DO ivar = 1, SIZE(var_info)
      !$ACC ENTER DATA ATTACH(var_info(ivar)%p_source) ASYNC(1)
    ENDDO
    DO ivar = 1, SIZE(sfc_var_info)
      !$ACC ENTER DATA ATTACH(sfc_var_info(ivar)%p_source) ASYNC(1)
    ENDDO
  END SUBROUTINE acc_copyin_var_info

  !>
  !!  Utility function.
  !!  Receive meteogram var list via MPI communication.
  !!
  SUBROUTINE receive_var_info(var_info, sfc_var_info, var_list, pack_buf)
    TYPE(t_var_info), ALLOCATABLE, TARGET, INTENT(out) :: var_info(:)
    TYPE(t_sfc_var_info), ALLOCATABLE, TARGET, INTENT(out) :: sfc_var_info(:)
    TYPE(t_var), INTENT(inout) :: var_list
    TYPE(mtgrm_pack_buf), INTENT(inout) :: pack_buf

#ifndef NOMPI
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":receive_var_info"
    INTEGER                 :: ivar, var_counts(2), ierror

    ! wait for messages to arrive:
    CALL p_wait()

    CALL p_unpack_int_1d(pack_buf%msg_varlist, pack_buf%pos, var_counts, 2)
    ALLOCATE(var_info(var_counts(1)), &
      &      sfc_var_info(var_counts(2)), stat=ierror)
    IF (ierror /= SUCCESS) &
      CALL finish (routine, 'ALLOCATE of var_info arrays failed.')


    ! from the received message, unpack the atmosphere/surface
    ! variables one by one:
    DO ivar = 1, var_counts(1)
      CALL unpack_cf(var_info(ivar)%cf, pack_buf)
      CALL p_unpack_int(pack_buf%msg_varlist, pack_buf%pos, &
        &               var_info(ivar)%igroup_id)
      CALL p_unpack_int(pack_buf%msg_varlist, pack_buf%pos, &
        &               var_info(ivar)%nlevs)
      NULLIFY(var_info(ivar)%p_source)
      IF (dbg_level > 0) &
        WRITE (*,*) "Added variable ", var_info(ivar)%cf%standard_name
    END DO
    DO ivar = 1, var_counts(2)
      CALL unpack_cf(sfc_var_info(ivar)%cf, pack_buf)
      CALL p_unpack_int(pack_buf%msg_varlist, pack_buf%pos, &
        &               sfc_var_info(ivar)%igroup_id)
      NULLIFY(sfc_var_info(ivar)%p_source)
      IF (dbg_level > 0) &
        WRITE (*,*) "Added variable ", sfc_var_info(ivar)%cf%standard_name
    END DO
    var_list%no_atmo_vars = var_counts(1)
    var_list%no_sfc_vars = var_counts(2)
#endif
  END SUBROUTINE receive_var_info

  !> fill t_cf_var struct for meteogram variable description from pack buffer
  SUBROUTINE unpack_cf(cf, pack_buf)
    TYPE(t_cf_var), INTENT(out) :: cf
    TYPE(mtgrm_pack_buf), INTENT(inout) :: pack_buf
    CALL p_unpack_string(pack_buf%msg_varlist, pack_buf%pos, cf%standard_name)
    CALL p_unpack_string(pack_buf%msg_varlist, pack_buf%pos, cf%long_name)
    CALL p_unpack_string(pack_buf%msg_varlist, pack_buf%pos, cf%units)
  END SUBROUTINE unpack_cf

  !> send time invariants in case output happens through a rank dedicated to asynchronous I/O
  SUBROUTINE send_time_invariants(var_info, invariants, &
    io_collector_rank, io_collect_comm, global_idx)
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_mtgrm_invariants), INTENT(in) :: invariants
    INTEGER, INTENT(in) :: io_collector_rank, io_collect_comm, global_idx(:)

    REAL(wp), ALLOCATABLE :: buf(:,:)

    INTEGER :: ivar, nvars, nlevs, pos, istation, nstations, ntiles

    ntiles = invariants%ntiles_mtgrm
    nstations = SIZE(global_idx)
    ALLOCATE(buf(num_time_inv + 2*ntiles + SUM(var_info%nlevs),nstations))
    nvars = SIZE(var_info)
    DO istation = 1, nstations
      pos = num_time_inv + 2 * ntiles
      buf(1,istation) = invariants%hsurf(istation)
      buf(2,istation) = invariants%frland(istation)
      buf(3,istation) = invariants%fc(istation)
      buf(4,istation) = REAL(invariants%soiltype(istation), wp)
      buf(5:6,istation) = REAL(invariants%tri_idx(istation,:), wp)
      buf(7:6+ntiles,istation) = invariants%tile_frac(:,istation)
      buf(7+ntiles:6+2*ntiles,istation) &
        = REAL(invariants%tile_luclass(:,istation), wp)
      DO ivar = 1, nvars
        nlevs = var_info(ivar)%nlevs
        buf(pos+1:pos+nlevs,istation) = invariants%heights(ivar)%a(istation,:)
        pos = pos + nlevs
      END DO
      CALL p_isend(buf(:,istation), io_collector_rank, &
        TAG_VARLIST+global_idx(istation), comm=io_collect_comm)
    END DO
    CALL p_wait
  END SUBROUTINE send_time_invariants

  !> recv and fill time invariant fields on dedicated async I/O rank from pack buffer
  SUBROUTINE recv_time_invariants(var_info, invariants, pstation, &
    io_collect_comm, is_pure_io_pe)
    TYPE(t_var_info), INTENT(in) :: var_info(:)
    TYPE(t_mtgrm_invariants), INTENT(inout) :: invariants
    INTEGER, INTENT(in) :: io_collect_comm, pstation(:)
    LOGICAL, INTENT(in) :: is_pure_io_pe

    REAL(wp), ALLOCATABLE :: buf(:,:)
    INTEGER :: ivar, nvars, nlevs, pos, istation, nstations, iowner, ntiles

    ntiles = invariants%ntiles_mtgrm
    nstations = SIZE(invariants%hsurf)
    ALLOCATE(buf(num_time_inv + 2*ntiles + SUM(var_info%nlevs),nstations))
    nvars = SIZE(var_info)
    DO istation = 1, nstations
      iowner = pstation(istation)
      IF ((is_pure_io_pe .OR. iowner /= p_pe_work) .AND. iowner >= 0) THEN
        CALL p_irecv(buf(:,istation), pstation(istation), &
          TAG_VARLIST+istation, comm=io_collect_comm)
      END IF
    END DO
    CALL p_wait()
    DO istation = 1, nstations
      iowner = pstation(istation)
      IF ((is_pure_io_pe .OR. iowner /= p_pe_work) .AND. iowner >= 0) THEN
        pos = num_time_inv + 2*ntiles
        invariants%hsurf(istation) = buf(1,istation)
        invariants%frland(istation) = buf(2,istation)
        invariants%fc(istation) = buf(3,istation)
        invariants%soiltype(istation) = INT(buf(4,istation))
        invariants%tri_idx(istation,:) = INT(buf(5:6,istation))
        invariants%tile_frac(:,istation) = buf(7:6+ntiles,istation)
        invariants%tile_luclass(:,istation) &
          = INT(buf(7+ntiles:6+2*ntiles,istation))
        DO ivar = 1, nvars
          nlevs = var_info(ivar)%nlevs
          invariants%heights(ivar)%a(istation,:) = buf(pos+1:pos+nlevs,istation)
          pos = pos + nlevs
        END DO
      END IF
    END DO
  END SUBROUTINE recv_time_invariants

END MODULE mo_meteogram_output
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
