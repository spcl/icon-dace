! Interface DACE observation operators from ICON
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

MODULE mo_icon2dace

!-----------------------------------------------------------------------------
!
! Description:
!   Interface DACE observation operators from ICON
!
!-----------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------

  !-------------
  ! ICON modules
  !-------------
  USE mo_kind,           ONLY: wp                    ! working precision kind
  USE mo_exception,      ONLY: finish, message
  USE mo_namelist,       ONLY: open_nml,            &! open namelist file
                               close_nml,           &! close namelist file
                               position_nml,        &! position to nml group
                               nnml,                &! Fortran namelist unit
                               POSITIONED            ! position_nml: OK return flag
  USE mo_util_vcs,       ONLY: util_repository_url, &!
                               util_branch_name,    &!
                               util_revision_key     !
  USE mo_grid_config,    ONLY: nroot, start_lev
  USE mo_gribout_config, ONLY: gribout_config
  USE mo_parallel_config,ONLY: nproma, idx_1d
  USE mo_model_domain,   ONLY: p_patch,             &!
               t_patch_icon => t_patch
  USE mo_communication,  ONLY: t_comm_pattern, t_comm_gather_pattern
  USE mo_loopindices,    ONLY: get_indices_c, get_indices_v
  USE mo_impl_constants, ONLY: min_rlcell, min_rlvert, min_rledge
  USE mo_sync,           ONLY: SYNC_C, SYNC_V, sync_patch_array

  USE mo_time_config,    ONLY: time_config
  USE mtime,             ONLY: &
                         MD => MAX_DATETIME_STR_LEN
  USE mo_util_mtime,     ONLY: getElapsedSimTimeInSeconds

  USE mo_ext_data_state, ONLY: ext_data       ! fr_land etc.
  USE mo_nonhydro_state, ONLY: p_nh_state
  USE mo_nonhydro_types, ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_state,  ONLY: prm_diag
  USE mo_nwp_phy_types,  ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_state,  ONLY: p_lnd_state
  USE mo_nwp_lnd_types,  ONLY: t_lnd_diag, t_lnd_prog
  USE mo_run_config,     ONLY: iqv, iqc, iqi, iqg, iqr, iqs
  USE mo_dynamics_config,ONLY: nnow, nnow_rcf
  USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info
  USE mo_assimilation_config,ONLY: assimilation_config
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_cell
  USE mo_intp_data_strc,        ONLY: p_int_state
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config 
  !-------------------
  ! ICON event control
  !-------------------
  USE mtime,             ONLY: datetime,  newDatetime,  deallocateDatetime, &
       &                       timedelta, newTimedelta, deallocateTimedelta,&
       &                       datetimetostring, timedeltaToString,         &
       &                       MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, &
       &                       MAX_MTIME_ERROR_STR_LEN, mtime_strerror,     &
       &                       OPERATOR(+), OPERATOR(-), OPERATOR(>=),      &
       &                       ASSIGNMENT(=), event, eventGroup, newEvent,  &
       &                       addEventToEventGroup, isCurrentEventActive
  USE mo_event_manager,  ONLY: addEventGroup, getEventGroup, printEventGroup

#ifdef __DACE__
  !-----------------------------
  ! DACE general purpose modules
  !-----------------------------
#ifndef NOMPI
  use mpi,            only: MPI_Comm_rank
#endif
  use mo_mpi_dace,    only: dace,            &! DACE communicator
                            set_dace_comm,   &! set the DACE MPI communicator
#ifdef __PGI
  ! PGI and NVIDIA compilers cannot choose the correct implementention of the
  ! generic procedure p_bcast.
                            p_bcast_alias => p_bcast, &
# define p_bcast p_bcast_alias
#else
                            p_bcast,         &! generic MPI bcast routine
#endif
                            p_ibcast,        &! generic MPI bcast routine (non-blocking)
                            p_sum,           &! generic MPI sum
                            p_allgather,     &! generic MPI allgather routine
                            p_waitall         ! wait for MPI requests to complete
  use mo_run_params,  only: nml_run_flags,   &! routine to read namelist /RUN/
                            nproc1,          &! partitioning parameter
                            nproc2,          &! partitioning parameter
                            path_file,       &! concatenate path+file name
                            output,          &! output files path
                            aux               ! auxiliary output path
! use mo_cpu_time,    only: stop_time         ! determine cpu and wall time
  use mo_time,        only: t_time,          &!
                            operator(+),     &! add times
                            operator(*),     &! multiply
                            time_cyyyymmddhhmmss
  use mo_dace_string, only: byte2hex          ! Convert to 'hexadecimal' CHAR(2)
  use mo_p_output,    only: flush_buf         ! write buffer to stdout
  use mo_fortran_units,only:get_unit_number, &!request free unit number
                            return_unit_number!return unit number
  use info_3dvar,     only: INFO_RevisionTag  ! Technical: enforce dependency
  !---------------------
  ! DACE grid handling,
  ! domain decomposition
  !---------------------
  use mo_atm_decomp,  only: t_atm_dec,       &! parallel decomposition data type
                            setup_parallel,  &! setup multi  processor run
                            setup_decomp,    &! setup domain decomposition
                            print_decomp,    &! print domain decomposition
                            setup_com_icon,  &! setup communication buffers
                            destruct          ! clean up t_atm_dec
  use mo_atm_grid,    only: t_grid, print,   &! grid information
                            construct,       &! constructor routine
                            destruct,        &! destructor  routine
                            allocate,        &! alloc. pointer components
                            deallocate,      &! dealloc. pointer components
                            set_ptopf,       &! set topmost level pressure
                            set_plev_indices,&! set some level indices
                            set_zlev_ref,    &! estimate unperturbed levels
                            set_empty_grid,  &! set default values for empty grid
                            set_grid_fields, &! mnemonics, bounds for specific fields
                            set_pointers,    &! update pointers
                            set_geoid,       &! set geoid
                            VCT_Z_GEN,       &! generalised z coordinate
                            MO_ICON,         &! ICON  model
                            setup_global_coord! cartesian coordinates of the mass points
  use mo_atm_state,   only: t_atm,           &! atmosphere derived type
                            construct,       &! constructor routine
                            destruct,        &! destructor  routine
                            print,           &! state information
                            allocate,        &! alloc. pointer components
                            assignment(=),   &! t_atm = ...
                            set_pointers,    &! update pointers
                            set_geo,         &! set geopotential
                            set_rh,          &! set relative humidity from q, t, p
                            deallocate        ! deallocate pointers
  use mo_memory,      only: t_m,             &! data type holding 3D fields
                            update            ! update field bounds, alloc.status
  use mo_wmo_tables,  only: DWD6_ICON,       &!
                            WMO3_GENV         !
  USE mo_t_icon,      ONLY: t_grid_cells,    &! Derived type for grid cells
                            t_grid_edges,    &! Derived type for grid edges
                            t_grid_vertices, &! Derived type for grid vertices
            t_patch_dace => t_patch,         &! Derived type for ICON patch
             nproma_dace => nproma,          &! Vector length for blocking
                            idx_no,          &! Block index from global index
                            blk_no            ! Line  index within block
  use mo_icon_grid,   only: t_grid_icon,     &! ICON grid metadata
                            check_neighbours,&! check neighbor refs for correct order
                            dist_to_bound,   &! determine distance to boundary
                            set_search_grid, &! auxiliary grid for interpolation
                            search_icon_global
  !use mo_letkf,            lcheck_ptopf      ! discard obs for plev < ptopf?
  !-----------------------
  ! observation processing
  !-----------------------
  use mo_t_obs,       only: read_obs_nml,    &! read namelist /observation/s
                            read_cdfin,      &! flag to read COSMO-CDFIN files
                            t_obs,           &! observation data type
                            t_spot,          &! component of t_obs
                            obs_local,       &! keep obs on gridpoint PE
                            fdb_files,       &! fof-files to read
                            obs_files,       &! CDFIN-Files (not COSMO) to read
                            add_source,      &! add source-file to list
                            source,          &! list of source-files
                            n_source,        &! number of source files
                            m_source,        &! max. number of source files
                            reorder_obs,     &! order observ. according to PEs
                            scatter_obs,     &! scatter observ. into boxes
                            release_mem,     &! release unused memory in t_obs
                            construct,       &! constructor routine for t_obs
                            destruct,        &! destructor  routine for t_obs
                            unique_spot,     &! ensure unique obs id
                            FT_FEEDBACK,     &! flag for feedback file
                            FT_MISSING,      &! flag for missing  file
                            TSK_SET_CHR,     &!
                            TSK_SETUP_OP,    &!
                            TSK_SETUP_COLS,  &!
                            TSK_SETUP_FULL,  &!
                            TSK_INIT,        &!
                            TSK_SHRINK,      &!
                            TSK_READ,        &!
                            TSK_R,           &!
                            TSK_Y,           &! flag to process_obs
                            n_dace_op,     &! lenght of non-empty entries in dace_op
                            COSMO             ! observation module numbers
  use mo_obs_set,     only: t_obs_set,       &! observation data type set
                            destruct,        &! t_obs_set destructor routine
                            t_obs_block       ! observation data type box
  use mo_obs_tables,  only: read_nml_report, &! set table rept_use from namelist
                            rept_use,        &!     table rept_use
                            rept_char,       &!     table rept_char
                            derive_rept_stat,&! re-derive obs.statistics
                            gather_rept_stat,&! gather statistics on I/O PE
                            print_rept_stat   ! print observ.type statistics
  use mo_obs,         only: process_obs,     &! general purpose routine
                            set_rules_cosmo, &! set default /RULES/ for COSMO
                            set_time_slice    ! set metadata for time interpol.
  use mo_obs_sndrcv,  only: p_alltoall,      &! MPI_ALLTOALL (t_obs)
                            p_bcast           ! MPI_BCAST    (t_obs)
  use mo_obs_rules,   only: read_nml_rules    ! read namelists /RULES/
  use mo_test_obs,    only: test_obs          ! test integrity of obs.data
  use mo_t_use,       only: STAT_PAS_REJ,    &! status flag: passive-rejected
                            STAT_REJECTED,   &! status flag: rejected
                            CHK_DOMAIN        ! out of domain check
  use mo_set_matrix,  only: PB_NONE,         &! no calculation of Pb flag
                            set_flags         ! set flags in t_obs
  use mo_std,         only: read_std_nml_dace ! read namelist /STD_OBS/
  use mo_boxes,       only: set_input_boxes, &! distribute files over PEs
                            set_veri_boxes    ! associate obs. with boxes
  use mo_thinning,    only: check_domain,    &! check for out of domain
                            thinning,        &! thinning routine
                            read_nml_thin     ! read namelist /THINNING/ 
  use mo_blacklist,   only: read_blacklists   ! read blacklist file
  use mo_fg_checks,   only: check_obs,       &! check for valid report
                            check_cons,      &! consistency check
                            check_rule,      &! check for rules
                            check_black,     &! apply blacklist
                            check_gross,     &! gross error check
                            check_suff,      &! check sufficient data in report
!                           check_operator,  &! check inapplicable obs.operator
                            check_obstype     ! check if obstype is in list
  use mo_t_enkf,      only: apply_H           ! apply obsv. operator
  use mo_cosmo_conv,  only: lqc               ! apply fg-check as in COSMO
  use mo_t_veri,      only: setup_veri_obs,  &! set up obs.interpol. operators
                            add_veri,        &! add entry to feedback file
                            prefix_out        ! feedback output file prefix
  use mo_tovs,        only: read_tovs_nml,   &! read TOVS_* namelists
                            use_reff,        &! use ICON effective radii in RTTOV
                            prep_H_btcs,     &! prepare t_vector
                            thin_superob_tovs !
  use mo_rad,         only: rad_set,         &! Options for radiance datasets
                            n_set
  !------------------------
  ! feedback file interface
  !------------------------
  use mo_fdbk_tables, only: init_fdbk_tables,&! initialise tables
                            OT_RAD,          &! obstype RAD
                            OF_BT_CLEAR_SKY   ! operator-flag for clear-sky TBs
  use mo_fdbk_3dvar,  only: t_fdbk_3dv,      &! 3dvar feedback file derived type
                            destruct,        &! t_fdbk_3dv destructor routine
                            write_fdbk_3dv    ! write (NetCDF) feedback file
  !---------------
  ! linear algebra
  !---------------
  use mo_dec_matrix,  only: t_vector,        &! vector derived type
                            construct,       &! t_vector constructor routine
                            destruct,        &! deallocate t_vector components
                            assignment(=)     ! assignment operator
  !---------------------
  ! model state handling
  !---------------------
  use mo_cntrlvar,    only: disable_gh        ! disable gen. humidity transform
  use mo_physics,     only: gacc,            &! gravity acceleration
                      Rd => R,               &! gas constant of dry air [J/(kg*K)]
                            r2d,             &! factor radians -> degree
                            d2r               ! factor degree -> radians
#endif /* __DACE__ */

  IMPLICIT NONE
  !----------------
  ! Public entities
  !----------------
  PRIVATE
  PUBLIC :: init_dace      ! initialise DACE as a subsystem
  PUBLIC :: init_dace_op   ! initialise DACE observation operators
  PUBLIC :: run_dace_op    ! run the DACE observation operators for a time step
  PUBLIC :: finish_dace    ! write fof-Files, clean up DACE
  PUBLIC :: mec_event      ! event handle for invoking the MEC
  PUBLIC :: dace_op_init   ! DACE operators initialized?

  !-----------------
  ! Module variables
  !-----------------
  logical, save        :: dace_op_init = .false. ! DACE operators initialized?
  PROTECTED            :: dace_op_init
  !----------------------
  ! Event control (timer)
  !----------------------
  TYPE(event),    POINTER :: mec_Event   => NULL()
  TYPE(datetime), POINTER :: mec_RefDate => NULL()
  !----------------------------------------------------------------------------
#ifdef __DACE__
  !----------------------------------------------------------------------------
  type(t_obs_set)      :: obs              ! observations and meta-data
  type(t_obs), target  :: obs_in(1)        ! observation input buffer
  type(t_obs_block)    :: obs_b            ! obs.input buffer container
  type(t_fdbk_3dv)     :: mon              ! feedback file meta data
  type(t_vector)       :: H_det            ! model equivalents
  type(t_atm), pointer :: atm(:) => NULL() ! model states
  !----------------------------------------------------------------------------
  ! DACE diagnostic output control
  !-------------------------------
  integer              :: dace_output_unit = 0  ! old version (backward compatible)
! integer              :: dace_output_unit = 6  ! new version (stdout redirection)
  character(255)       :: dace_output_file = "" ! filename
  !----------------------------------------------------------------------------
  ! ICON model information
  !-----------------------
  character(256) :: icon_repo     = ""
  character(256) :: icon_branch   = ""
  character(256) :: icon_revision = ""
  integer        :: grid_root     = -1
  integer        :: grid_level    = -1
  character      :: grid_uuid(16) = ACHAR (0)
  integer        :: exp_id        = -1  ! localNumberOfExperiment
  integer        :: ens_id        = -1  ! localTypeOfEnsembleForecast
  integer        :: ens_mem       = -1  ! perturbationNumber
  integer        :: num_mem       = -1  ! numberOfForecastsInEnsemble
  !------------------
  ! DACE time control
  !------------------
  integer                   :: dace_time_ctrl(3) = 0  ! Start,end,step [s]
  integer                   :: mec_calls         = 0
  integer                   :: mec_slot1         = 0
  integer                   :: mec_slotn         = 0
  integer                   :: mec_sloti         = -1
  type(t_time),save         :: exp_start
  type(t_time),save         :: exp_end
  type(t_time),save         :: tinc           ! time increments of model fields
  !----------------------------------------------------------------------------
  ! ICON grid information
  ! domain decomposition
  !----------------------
  integer,  allocatable :: n_cells(:)    ! 0..npe-1
  integer,  allocatable :: n_verts(:)    ! 0..npe-1
  integer,  allocatable :: n_edges(:)    ! 0..npe-1
  !----------------------------------------------------------------------------
  type(t_atm_dec)       :: icon_dc       ! ICON domain decomposition (owners)
  type(t_atm_dec)       :: icon_dc_d     ! ICON domain decomposition, dual grid
  type(t_grid), pointer :: grid => NULL()
  integer     , pointer :: marr_c(:,:)   ! (3,i) cells index field: (pe,i,j)
  integer     , pointer :: marr_v(:,:)   ! (3,i) verts index field: (pe,i,j)
  !----------------------------------------------------------------------------
  integer,  allocatable :: cell_glb_idx(:,:,:) ! cells    global indices
  integer,  allocatable :: vert_glb_idx(:,:,:) ! vertices global indices
  integer,  allocatable :: ext_glb_idx (:)     ! external cells global indices
  !----------------------------------------------------------------------------
  integer               :: dbg_level = 1 ! 0=no debug, 1=minimal, 2=full
  !----------------------------------------------------------------------------
  !-------------------
  ! Namelist /mec_obs/
  !-------------------
                      ! obstypes to handle (tested so far):
  character(len=512) :: obstypes     = 'TEMP PILOT SYNOP DRIBU AIREP SATOB'
  integer            :: interpolation= -1    ! >0: use nth slot, 0:nn -1:interval
  integer            :: fg_check     = -1    ! switch how to apply quality checks

  namelist /mec_obs/ obstypes, prefix_out, interpolation, fg_check

#endif /* __DACE__ */
  !============================================================================
contains
  !============================================================================
  subroutine init_dace (comm, p_io, ldetached)
    !----------------------------------------
    ! set the DACE MPI communicator from ICON
    ! initialize domain decomposition, grid
    !----------------------------------------
    integer ,intent(in) :: comm  ! communicator to use
    integer ,intent(in) :: p_io  ! PE to use for I/O
    logical, intent(in) :: ldetached ! indicates that current PE is a detached IO-PE
    !----------------
    ! Local variables
    !----------------
    integer       :: n
    integer       :: ios     ! I/O status
    integer       :: iu      ! temp. unit number
    character(MD) :: cdate   ! yyyy-mm-ddThh:mm:ss.000
    character(14) :: adate   ! yyyymmddhhmmss
    character(14) :: refdate ! yyyymmddhhmmss

    call message ("","")
    call message ("icon2dace","initializing DACE coupling")

#ifndef __DACE__
    CALL finish ("init_dace","DACE coupling requested but not compiled in")
#else

    dbg_level = max (dbg_level, assimilation_config(1)% dace_debug)

    if (assimilation_config(1)% dace_output_file /= "") then
       dace_output_unit = 6         ! new version (stdout redirection)
       dace_output_file = adjustl (assimilation_config(1)% dace_output_file)
    else
       dace_output_unit = 0         ! old version (backward compatible)
    end if

    IF (.not. ldetached) then
      call set_dace_comm (comm, p_io)
    ENDIF

    !--------------------------
    ! Set up time slots for MEC
    !--------------------------
    call set_dace_timer () ! the mtime event needs to be created on a detached IO-PE as well;
    IF (ldetached) RETURN  ! afterwards, this process is excluded from DACE coupling

    !---------------------------
    ! Get ICON model information
    !---------------------------
    n = len (icon_repo    ); call util_repository_url (icon_repo,     n)
    n = len (icon_branch  ); call util_branch_name    (icon_branch,   n)
    n = len (icon_revision); call util_revision_key   (icon_revision, n)

    !-----------------
    ! Grid information
    !-----------------
    grid_root  = nroot
    grid_level = start_lev
    grid_uuid  = transfer (p_patch(1)% grid_uuid, grid_uuid)

    if (dace% lpio .and. dbg_level > 1) then
       write(0,*) "root,level =", grid_root, grid_level
       write(0,*) "grid_uuid  = '", byte2hex (grid_uuid), "'"
       write (1000+dace% pe,*) "#pe      =",dace% pe
       write (1000+dace% pe,*) "#cells_g =",p_patch(1)% n_patch_cells_g
       write (1000+dace% pe,*) "#cells   =",p_patch(1)% n_patch_cells
       write (1000+dace% pe,*) "#nblks_c =",p_patch(1)% nblks_c
       write (1000+dace% pe,*) "#npromz_c=",p_patch(1)% npromz_c
       write (1000+dace% pe,*) "#nlev    =",p_patch(1)% nlev
       write (1000+dace% pe,*) "#shape   =",shape  (p_patch(1)% cells% center)
       write (1000+dace% pe,*) "#shape (owner_mask):", &
            shape (p_patch(1)% cells% decomp_info% owner_mask)
       write (1000+dace% pe,*) "#fr_land =",shape  (ext_data(1)% atm% fr_land)
       write (1000+dace% pe,*) "#hhl     =",shape  (p_nh_state(1)% metrics% z_ifc)
    end if

    call dace_open_output (first=.true.)

    !---------------------------
    ! Start date (analysis date)
    !---------------------------
    call datetimetostring (time_config%tc_exp_startdate, cdate)
    adate = cdate(1:4) // cdate(6:7) // cdate(9:10)    &
         // cdate(12:13) // cdate(15:16) // cdate(18:19)
    exp_start = time_cyyyymmddhhmmss (adate)

    call datetimetostring (time_config%tc_exp_stopdate, cdate)
    refdate = cdate(1:4) // cdate(6:7) // cdate(9:10)    &
         // cdate(12:13) // cdate(15:16) // cdate(18:19)
    exp_end = time_cyyyymmddhhmmss (refdate)

    if (dace% lpio .and. dbg_level > 1) then
       write(0,*) "cdate = '", trim (cdate), "'"
       write(0,*) "adate = '", trim (adate), "'"
    end if

    exp_id  = gribout_config(1)% localNumberOfExperiment
    ens_id  = gribout_config(1)% localTypeOfEnsembleForecast
    ens_mem = gribout_config(1)% perturbationNumber
    num_mem = gribout_config(1)% numberOfForecastsInEnsemble
    if (dace% lpio) then
      write(0,*) "icon2dace: ens_id  = ", ens_id
      write(0,*) "icon2dace: ens_mem = ", ens_mem
      write(0,*) "icon2dace: num_mem = ", num_mem
    end if

    call message ("icon2dace","set custom domain decomposition")
    call setup_dace_decomp ()

    call message ("icon2dace","reading DACE namelist /RUN/")
    !=================================================
    ! open namelist file (unit: nnml)
    ! read namelists:
    !   /RUN/ for I/O paths (obsinput, output,...)
    !=================================================
    if (dace% lpio) call open_nml (assimilation_config(1)%dace_namelist_file)
    call nml_run_flags ()
    flush (6)
    call message ("","")
    !----------------------------------------------------------------
    ! Check whether given directories 'output' and 'aux' are writable
    !----------------------------------------------------------------
    if (dace% lpio) then
       iu = get_unit_number ()
       open (iu, file=path_file (output,"o.test"), status='new', iostat=ios)
       if (ios /= 0) call finish ("init_dace_op",                          &
                                  "directory not writable: "//trim (output))
       close (iu, status='delete')
       if (aux /= output) then
          open (iu, file=path_file (aux,"a.test"), status='new', iostat=ios)
          if (ios /= 0) call finish ("init_dace_op",                       &
                                     "directory not writable: "//trim (aux))
          close (iu, status='delete')
       end if
       call return_unit_number (iu)
    end if
    if (dace% lpio) call close_nml ()

    call message ("icon2dace","grid setup")
    call message ("","")

    !-------------------------------------------
    ! Set up ICON grid metadata required by DACE
    !-------------------------------------------
    allocate (grid% icongrid)
    call set_global_indices (grid)
    call icongrid_from_icon (grid% icongrid, p_patch(1),comm)
    call grid_from_icon     (grid)
    deallocate              (cell_glb_idx, vert_glb_idx)

    call message ("icon2dace","grid setup completed")

    call message ("icon2dace","init_dace done.")

    call dace_close_output ()

#endif /* __DACE__ */

  end subroutine init_dace
  !============================================================================
#ifdef __DACE__

  subroutine setup_dace_decomp ()
    !----------------------------------------------------------
    ! Set up domain decomposition info for "linear" addressing.
    ! We keep the basic domain decomposition of ICON, but do
    ! not use NPROMA blocking.
    !----------------------------------------------------------
    integer               :: n, i
    integer               :: lb  (4)    ! lower bounds
    integer               :: ub  (4)    ! upper bounds
    integer               :: lb_d(4)    ! lower bounds (dual grid)
    integer               :: ub_d(4)    ! upper bounds (dual grid)
    integer               :: n_cells_g  ! global number of cells
    integer               :: n_verts_g  ! global number of vertices
    type(t_grid), pointer :: g          ! Grid template
    integer, allocatable  :: ilim1  (:) ! domain decomposition
    integer, allocatable  :: ilim2  (:) ! domain decomposition
    integer, allocatable  :: ilim1_d(:) ! domain decomposition
    integer, allocatable  :: ilim2_d(:) ! domain decomposition

    !-----------------------------------
    ! Get domain decomposition from ICON
    !-----------------------------------
    allocate (n_cells(0:dace% npe-1))
    allocate (n_verts(0:dace% npe-1))
    allocate (n_edges(0:dace% npe-1))

    n = count (p_patch(1)% cells% decomp_info% owner_mask)
    call p_allgather (n, n_cells)
    n = count (p_patch(1)% verts% decomp_info% owner_mask)
    call p_allgather (n, n_verts)
    n = count (p_patch(1)% edges% decomp_info% owner_mask)
    call p_allgather (n, n_edges)

    n_cells_g = p_patch(1)% n_patch_cells_g
    n_verts_g = p_patch(1)% n_patch_verts_g
    if (sum (n_cells) /= n_cells_g) then
       write (0,*) "#cells_g =",n_cells_g, "/=", sum (n_cells)
       call finish ("setup_dace_decomp","cells mismatch")
    end if
    if (sum (n_verts) /= n_verts_g) then
       write (0,*) "#verts_g =",n_verts_g, "/=", sum (n_verts)
       call finish ("setup_dace_decomp","vertices mismatch")
    end if

    if (dace% lpio .and. dbg_level > 1) then
       write (1000+dace% pe,*)
       write (1000+dace% pe,*) "#owner_mask:"
       write (1000+dace% pe,*) "#cells_g =",p_patch(1)% n_patch_cells_g, sum (n_cells)
       write (1000+dace% pe,*) "#n_cells =", n_cells
       write (1000+dace% pe,*) "#verts_g =",p_patch(1)% n_patch_verts_g, sum (n_verts)
       write (1000+dace% pe,*) "#n_verts =", n_verts
       write (1000+dace% pe,*) "#edges_g =",p_patch(1)% n_patch_edges_g, sum (n_edges)
       write (1000+dace% pe,*) "#n_edges =", n_edges
    end if

    !--------------------------------------------------------
    ! Set up grid template (c.f. mo_atm_grid::construct_grid)
    !--------------------------------------------------------
    allocate (g)
    call set_empty_grid (g)
    g% gridtype    = DWD6_ICON
    g% nxny        = n_cells_g
    g% nx          = n_cells_g
    g% ny          = n_cells_g
    g% nz          = p_patch(1)% nlev
    g% nd          = 1
    g% lbg         = 1
    g% ubg         = [ g% nxny   , 1 , g% nz , 1 ]
    g% lb          = g% lbg
    g% ub          = g% ubg
    g% lb(1)       = sum (n_cells(0:dace% pe-1)) + 1
    g% ub(1)       = sum (n_cells(0:dace% pe  ))
    g% shape       = g% ub - g% lb + 1
    g% size        = g% nxny * max (1, g% nz)

    g% nxny_d      = n_verts_g
    g% nx_d        = n_verts_g
    g% ny_d        = n_verts_g
    g% lbg_d       = 1
    g% ubg_d       = [ g% nxny_d , 1 , g% nz , 1 ]
    g% lb_d        = g% lbg_d
    g% ub_d        = g% ubg_d
    g% lb_d(1)     = sum (n_verts(0:dace% pe-1)) + 1
    g% ub_d(1)     = sum (n_verts(0:dace% pe  ))

    g% model       = MO_ICON
    g% levtyp      = WMO3_GENV
    g% vct         = VCT_Z_GEN  ! DACE internal vertical coord.type
    g% vc% ivctype = VCT_Z_GEN

    !------------------------------------------------------
    ! Set up custom decomposition matching DACE conventions
    !------------------------------------------------------
    nproc1 = dace% npe
    nproc2 = 1
    allocate (ilim1  (0:nproc1))
    allocate (ilim1_d(0:nproc1))
    allocate (ilim2  (0:nproc2))
    allocate (ilim2_d(0:nproc2))
    ilim1  (0) = 1
    ilim1_d(0) = 1
    do i = 1, nproc1
       ilim1  (i) = ilim1  (i-1) + n_cells(i-1)
       ilim1_d(i) = ilim1_d(i-1) + n_verts(i-1)
    end do
    ilim2  (0) = 1
    ilim2  (1) = 2
    ilim2_d(0) = 1
    ilim2_d(1) = 2

    call setup_parallel (g% dc,     nproc1, nproc2, dace% comm)
    call setup_parallel (g% dc_d,   nproc1, nproc2, dace% comm)
    call setup_parallel (icon_dc,   nproc1, nproc2, dace% comm)
    call setup_parallel (icon_dc_d, nproc1, nproc2, dace% comm)
    lb    = g% lbg
    ub    = g% ubg
    call setup_decomp   (icon_dc,   lb, ub, ilim1, ilim2)
    lb    = g% lbg
    ub    = g% ubg
    call setup_decomp   (g% dc,     lb, ub, ilim1, ilim2)
    CALL setup_com_icon (g% dc,     g% marr, g% lbg, g% ubg)
!   call print_decomp   (g% dc, unit=1000+dace% pe)

    lb_d  = g% lbg_d
    ub_d  = g% ubg_d
    call setup_decomp   (icon_dc_d, lb_d, ub_d, ilim1_d, ilim2_d)
    lb_d  = g% lbg_d
    ub_d  = g% ubg_d
    call setup_decomp   (g% dc_d,   lb_d, ub_d, ilim1_d, ilim2_d)

    !-------------------
    ! Consistency checks
    !-------------------
    if (g% lb(1) /= g% dc% ilim1(dace% pe)       .or. &
        g% ub(1) /= g% dc% ilim1(dace% pe+1) - 1      ) then
       write(unit=1000+dace% pe, fmt=*) "###(ilim1):", &
            g% lb(1), g% dc% ilim1(dace% pe),    &
            g% ub(1), g% dc% ilim1(dace% pe+1) - 1
       call finish ("setup_dace_decomp","ilim1")
    end if
    if (g% lb(2) /= g% dc% ilim2(0)       .or. &
        g% ub(2) /= g% dc% ilim2(1) - 1      ) then
       write(unit=1000+dace% pe, fmt=*) "###(ilim2):", &
            g% lb(2), g% dc% ilim2(0),  &
            g% ub(2), g% dc% ilim2(1) - 1
       call finish ("setup_dace_decomp","ilim2")
    end if

    !---------------------------------------------------------------
    ! reference to 'original' indices and processor, ICON convention
    !---------------------------------------------------------------
    call setup_com_custom (icon_dc,   marr_c, g% lbg,   g% ubg,   &
                           p_patch(1)% cells% decomp_info         )
    call setup_com_custom (icon_dc_d, marr_v, g% lbg_d, g% ubg_d, &
                           p_patch(1)% verts% decomp_info         )

    grid => g

  end subroutine setup_dace_decomp
  !----------------------------------------------------------------------------
  subroutine set_global_indices (grid)
    type(t_grid), intent(in) :: grid

    integer               :: lb  (4)    ! lower bounds
    integer               :: ub  (4)    ! upper bounds
    integer               :: lb_d(4)    ! lower bounds (dual grid)
    integer               :: ub_d(4)    ! upper bounds (dual grid)
    integer               :: i
    !----------------------
    ! Derive global indices
    !----------------------
    lb  (1) = minval (marr_c(2,:))
    ub  (1) = maxval (marr_c(2,:))
    lb  (2) = minval (marr_c(3,:))
    ub  (2) = maxval (marr_c(3,:))
    lb_d(1) = minval (marr_v(2,:))
    ub_d(1) = maxval (marr_v(2,:))
    lb_d(2) = minval (marr_v(3,:))
    ub_d(2) = maxval (marr_v(3,:))
    lb  (3) = 0
    ub  (3) = dace% npe - 1
    allocate (cell_glb_idx(lb  (1):ub  (1),lb  (2):ub  (2),lb(3):ub(3)))
    allocate (vert_glb_idx(lb_d(1):ub_d(1),lb_d(2):ub_d(2),lb(3):ub(3)))
    cell_glb_idx = -1
    vert_glb_idx = -1

    do i = grid% lbg  (1), grid% ubg  (1)
       cell_glb_idx(marr_c(2,i),marr_c(3,i),marr_c(1,i)) = i
    end do
    do i = grid% lbg_d(1), grid% ubg_d(1)
       vert_glb_idx(marr_v(2,i),marr_v(3,i),marr_v(1,i)) = i
    end do

    if (dace% lpio .and. dbg_level > 1) then
       write (1000+dace% pe,*)
       write (1000+dace% pe,*) "### ub  (1:2)", lb(1:2), ub(1:2)
       write (1000+dace% pe,*) "### ub_d(1:2)", lb_d(1:2), ub_d(1:2)
       write (1000+dace% pe,*) "### count (cell_glb_idx > 0) =",  count (cell_glb_idx > 0)
       write (1000+dace% pe,*) "### count (vert_glb_idx > 0) =",  count (vert_glb_idx > 0)
       write (1000+dace% pe,*) "### count (cell_glb_idx < 1) =",  count (cell_glb_idx < 1)
       write (1000+dace% pe,*) "### count (vert_glb_idx < 1) =",  count (vert_glb_idx < 1)
       write (1000+dace% pe,*) "### count (cell_glb_idx > 0,pe=0) =", count (cell_glb_idx(:,:,0) > 0)
       write (1000+dace% pe,*) "### count (vert_glb_idx > 0,pe=0) =", count (vert_glb_idx(:,:,0) > 0)
       write (1000+dace% pe,*) "### count (cell_glb_idx < 1,pe=0) =", count (cell_glb_idx(:,:,0) < 1)
       write (1000+dace% pe,*) "### count (vert_glb_idx < 1,pe=0) =", count (vert_glb_idx(:,:,0) < 1)
    end if
  end subroutine set_global_indices
  !----------------------------------------------------------------------------
  SUBROUTINE setup_com_custom (dc, marr, lbg, ubg, decomp_info)
    type(t_atm_dec)                ,INTENT(in) :: dc
    INTEGER                        ,POINTER    :: marr(:,:)    ! (3,i)
    INTEGER                        ,INTENT(in) :: lbg (4)
    INTEGER                        ,INTENT(in) :: ubg (4)
    TYPE(t_grid_domain_decomp_info),INTENT(in) :: decomp_info
    !---------------------------------------------------
    ! derive communication info for custom "linear" grid
    !---------------------------------------------------
    integer              :: j, j1, j2      ! Loop indices
    integer              :: pe             ! Processor index
    integer              :: nb, nl, ib, il ! Block/line indices
    integer, allocatable :: blk(:), idx(:) ! Block/line index lists

    if (ubg(2) /= lbg(2) .or. dc% nproc1 /= dace% npe .or. dc% nproc2 /= 1)  &
         call finish('setup_com_custom','unimplemented/invalid decomposition')

    !-------------------------------
    ! Determine locally owned points
    !-------------------------------
    allocate (idx(lbg(1):ubg(1)))
    allocate (blk(lbg(1):ubg(1)))
    idx = -1
    blk = -1
    nl = size (decomp_info% owner_mask, dim=1)
    nb = size (decomp_info% owner_mask, dim=2)
    j  = dc% ilim1(dace% pe)
    do ib = 1, nb
       do il = 1, nl
          if (decomp_info% owner_mask(il,ib)) then
             idx(j) = il
             blk(j) = ib
             j = j + 1
          end if
       end do
    end do
    if (j /= dc% ilim1(dace% pe+1)) &
         call finish ('setup_com_custom','owner mask?')

    !----------------------------
    ! Exchange decomposition info
    !----------------------------
    do pe = 0, dc% nproc1-1
       j1 = dc% ilim1(pe)
       j2 = dc% ilim1(pe+1) - 1
!      call p_bcast (idx(j1:j2), pe)
!      call p_bcast (blk(j1:j2), pe)
       call p_ibcast (idx(j1:j2), pe)
       call p_ibcast (blk(j1:j2), pe)
    end do
    call p_waitall ()
    if (any (idx < 0)) call finish ("setup_com_custom", "bad idx")
    if (any (blk < 0)) call finish ("setup_com_custom", "bad blk")

    allocate (marr(3,lbg(1):ubg(1)))
    marr(:,:) = -1

    do pe = 0, dc% nproc1-1
       j1 = dc% ilim1(pe)
       j2 = dc% ilim1(pe+1) - 1
       do j = j1, j2
          marr(1,j) = pe
          marr(2,j) = idx(j)
          marr(3,j) = blk(j)
       end do
    end do
    if (any (marr < 0))      call finish ("setup_com_custom", "bad marr")
    if (any (marr > ubg(1))) call finish ("setup_com_custom", "marr ???")

  end SUBROUTINE setup_com_custom
  !----------------------------------------------------------------------------
  subroutine icongrid_from_icon (icongrid, patch, comm, verbose)
    type(t_grid_icon),  intent(out)   :: icongrid ! DACE: metadata of ICON grid
    type(t_patch_icon), intent(inout) :: patch    ! ICON native patch
    integer,            intent(in)    :: comm     ! communicator handle
    logical, optional,  intent(in)    :: verbose  ! Enable debugging
    !-------------------------------------------
    ! Gather ICON grid metadata required by DACE
    !-------------------------------------------
    type(t_patch_dace), pointer :: p              ! Patch data type, DACE
    logical                     :: global         ! global grid?
    integer                     :: euler          ! Euler characteristic
    integer                     :: j, k, pe       ! Loop indices
    integer                     :: n_cells, n_edges, n_verts
    integer                     :: nproma_c, nproma_e, nproma_v
    integer                     :: j1, j2, k1, k2, ne
    integer                     :: jc, jb, kc, kb, nc, nb, p_pe, ierr
    INTEGER                     :: i_startblk, i_endblk
    INTEGER                     :: i_startidx, i_endidx
    integer,  allocatable       :: i1(:), i2(:), i3(:)
    real(wp), allocatable       :: c1(:), c2(:), d1(:), d2(:)
    integer,  allocatable       :: owner_c(:,:,:)
    integer,  allocatable       :: owner_v(:,:,:)
    integer                     :: ou

    ou = dace_output_unit

    ! This is needed for consistency with set_dace_comm
#ifndef NOMPI
    call MPI_COMM_RANK (comm, p_pe, ierr)
#else
    p_pe = 0
#endif

    !-----------------------------------------------------
    ! Reconstruct owner information for cells and vertices
    !-----------------------------------------------------
    allocate (owner_c(nproma,3,patch% nblks_c))
    owner_c    = -1
    i_startblk = 1
    i_endblk   = patch% nblks_c
    do jb = i_startblk,i_endblk
       call get_indices_c (patch, jb, i_startblk, i_endblk,   &
                           i_startidx, i_endidx, 1, min_rlcell)
       do jc = i_startidx,i_endidx
          owner_c(jc,1,jb) = p_pe
          owner_c(jc,2,jb) = jc
          owner_c(jc,3,jb) = jb
       end do
    end do
    call sync_patch_array (SYNC_C, patch, owner_c)
    ! Kontrolle
    do jb = i_startblk, i_endblk
       call get_indices_c (patch, jb, i_startblk, i_endblk,   &
                           i_startidx, i_endidx, 1, min_rlcell)
       do jc = i_startidx,i_endidx
          if (.not. patch% cells% decomp_info% owner_mask(jc,jb) .and. p_pe==owner_c(jc,1,jb)) then
             write(dace%pe + 1000,*) p_pe,"%%% pe,jc,jb=", owner_c(jc,1:3,jb)
          endif
          if (any (owner_c(jc,1:3,jb) < 0)) then
             write(dace%pe + 1000,*) p_pe,"#%# pe,jc,jb=", owner_c(jc,1:3,jb)
          endif
       enddo
    enddo

    allocate (owner_v(nproma,3,patch% nblks_c))
    owner_v    = -1
    i_startblk = 1
    i_endblk   = patch% nblks_v
    do jb = i_startblk,i_endblk
       call get_indices_v (patch, jb, i_startblk, i_endblk,   &
                           i_startidx, i_endidx, 1, min_rlvert)
       do jc = i_startidx,i_endidx
          owner_v(jc,1,jb) = p_pe
          owner_v(jc,2,jb) = jc
          owner_v(jc,3,jb) = jb
       end do
    end do
    call sync_patch_array (SYNC_V, patch, owner_v)
    ! Kontrolle
    do jb = i_startblk, i_endblk
       call get_indices_v (patch, jb, i_startblk, i_endblk,   &
                           i_startidx, i_endidx, 1, min_rlvert)
       do jc = i_startidx,i_endidx
          if (.not. patch% verts% decomp_info% owner_mask(jc,jb) .and. p_pe==owner_v(jc,1,jb)) then
             write(dace%pe + 1000,*) p_pe,"%%% pe,jc,jb=", owner_v(jc,1:3,jb)
          endif
          if (any (owner_v(jc,1:3,jb) < 0)) then
             write(dace%pe + 1000,*) p_pe,"#%# pe,jc,jb=", owner_v(jc,1:3,jb)
          endif
       enddo
    enddo

    n_cells = patch% n_patch_cells_g
    n_verts = patch% n_patch_verts_g
    n_edges = patch% n_patch_edges_g

    icongrid% grid_root  = nroot
    icongrid% grid_level = start_lev
    icongrid% uuid       = transfer (patch% grid_uuid, icongrid% uuid)

    allocate (icongrid% patch)
    p => icongrid% patch
    p% cell_type       = 3
    p% n_patch_cells_g = n_cells
    p% n_patch_verts_g = n_verts
    p% n_patch_edges_g = n_edges
    p% n_patch_cells   = n_cells
    p% n_patch_verts   = n_verts
    p% n_patch_edges   = n_edges

    euler = p% n_patch_verts_g - p% n_patch_edges_g + p% n_patch_cells_g
    select case (euler)
    case (2)
       global = .true.
    case (1)
       global = .false.
    case default
       call finish ("icongrid_from_icon", "invalid grid")
    end select
    icongrid% global = global
    if (dace% lpio) then
      if (dbg_level > 0) then
        write(ou,*) "# patch is global     = ", icongrid% global
      end if
    end if

    !-----------------------------------------
    ! compute the no. of blocks (currently 1):
    !-----------------------------------------
    p%nblks_c  = blk_no (p%n_patch_cells)
    p%npromz_c = idx_no (p%n_patch_cells)
    p%nblks_e  = blk_no (p%n_patch_edges)
    p%npromz_e = idx_no (p%n_patch_edges)
    p%nblks_v  = blk_no (p%n_patch_verts)
    p%npromz_v = idx_no (p%n_patch_verts)

    ! Handle case where nproma is ridiculously large (DACE)
    nproma_c   = min (nproma_dace, n_cells)
    nproma_e   = min (nproma_dace, n_edges)
    nproma_v   = min (nproma_dace, n_verts)

    ALLOCATE(p%verts%vertex       (nproma_v,p%nblks_v),             &! mand.
      &      p%verts%num_edges    (nproma_v,p%nblks_v),             &! mand.
      &      p%verts%neighbor_idx (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%verts%neighbor_blk (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%verts%cell_idx     (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%verts%cell_blk     (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%cells%center       (nproma_c,p%nblks_c),             &! mand.
      &      p%cells%neighbor_idx (nproma_c,p%nblks_c,p%cell_type), &! mand.
      &      p%cells%neighbor_blk (nproma_c,p%nblks_c,p%cell_type), &! mand.
      &      p%cells%c_ctrl       (nproma_c,p%nblks_c)              )!

    !------------------------------------------------------------------
    ! Derive global index corresponding to external data representation
    ! (cells only so far)
    !------------------------------------------------------------------
    allocate (ext_glb_idx(n_cells))
    ext_glb_idx = -1
    j1 = icon_dc% ilim1(dace% pe)
    j2 = icon_dc% ilim1(dace% pe+1) - 1
    do j = j1, j2
       jc = marr_c(2,j)
       jb = marr_c(3,j)
       ext_glb_idx(j) = patch% cells% decomp_info% glb_index(idx_1d(jc,jb))
    end do
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc% ilim1(pe)
       j2 = icon_dc% ilim1(pe+1) - 1
       call p_ibcast (ext_glb_idx(j1:j2), pe)
    end do
    call p_waitall ()

    if (dace% lpio .and. dbg_level > 1) then
       write(dace%pe + 1000,*) "### ext_glb_idx =", ext_glb_idx(:)
    end if

    if (any (ext_glb_idx < 1)) then
       if (dace% lpio) then
          do j = 1, n_cells
             if (ext_glb_idx(j) < 1) &
                  write(ou,*) "### BAD:", j,":", ext_glb_idx(j)
          end do
       end if
       call finish ("icongrid_from_icon","bad global index")
    end if

    !------
    ! Cells
    !------
    allocate (c1(n_cells), c2(n_cells))

    !------------------------
    ! Cell center coordinates
    !------------------------
    c1 = -HUGE(0._wp)
    c2 = -HUGE(0._wp)
    j1 = icon_dc% ilim1(dace% pe)
    j2 = icon_dc% ilim1(dace% pe+1) - 1
    do j = j1, j2
       jc = marr_c(2,j)
       jb = marr_c(3,j)
       c1(j) = patch% cells% center(jc,jb)% lon  ! [rad]
       c2(j) = patch% cells% center(jc,jb)% lat  ! [rad]
    end do
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc% ilim1(pe)
       j2 = icon_dc% ilim1(pe+1) - 1
       call p_ibcast (c1(j1:j2), pe)
       call p_ibcast (c2(j1:j2), pe)
    end do
    call p_waitall ()
    do j = 1, n_cells
       p% cells% center(j,1)% lon = c1(j)
       p% cells% center(j,1)% lat = c2(j)
    end do
    deallocate (c1, c2)

    !---------------
    ! Neighbor cells
    !---------------
    p% cells% neighbor_idx(:,:,:) = -1
    p% cells% neighbor_blk(:,:,:) = -1

    j1 = icon_dc% ilim1(dace% pe)
    j2 = icon_dc% ilim1(dace% pe+1) - 1
    do j = j1, j2
       jc = marr_c(2,j)
       jb = marr_c(3,j)
       ne = 0
       if (.not. patch% cells% decomp_info% owner_mask(jc,jb)) then
          write(1000+dace% pe,*) "##",j,"something really bad happened for",jc,jb,&
               "real owner:", owner_c(jc,1,jb)
       end if
       do k = 1, p%cell_type
          nc = patch% cells% neighbor_idx(jc,jb,k)
          nb = patch% cells% neighbor_blk(jc,jb,k)
          if (nc <= 0) then
             cycle  ! k loop
          end if
          ne = ne + 1
          pe = dace% pe
          if (.not. patch% cells% decomp_info% owner_mask(nc,nb)) then
             pe = owner_c(nc,1,nb)
             kc = owner_c(nc,2,nb)
             kb = owner_c(nc,3,nb)
             if (pe == dace% pe) then
                write(1000+dace% pe,*) "##",j,k,"jc,jb=",jc,jb,"bad cell_neighbor?",&
                     patch% cells% decomp_info% owner_mask(jc,jb),nc,nb, &
                     patch% cells% decomp_info% owner_mask(nc,nb),pe,kc,kb
             else
                nc = kc
                nb = kb
             end if
          end if
          if (nc > 0 .and. nb > 0) then
             p% cells% neighbor_idx(j,1,k) = cell_glb_idx(nc,nb,pe)
             p% cells% neighbor_blk(j,1,k) = 1
          end if
          if (p% cells% neighbor_idx(j,1,k) < 1) then
             kc = patch% cells% neighbor_idx(jc,jb,k)
             kb = patch% cells% neighbor_blk(jc,jb,k)
             write(1000+dace% pe,*) "##",j,k,":cell_neighbor:",nc,nb,pe,&
                  p% cells% neighbor_idx(j,1,k),jc,jb,owner_c(kc,1:3,kb)
          end if
       end do
    end do
    !---------------------------------
    ! Communicate p% cells% neighbor_*
    !---------------------------------
    allocate (i1(n_cells*p%cell_type), i2(n_cells*p%cell_type))
    i1 = -1
    i2 = -1
    j1 = icon_dc% ilim1(dace% pe)
    j2 = icon_dc% ilim1(dace% pe+1) - 1
    k1 = (j1-1)*p%cell_type+1
    k2 =  j2   *p%cell_type
    i1(k1:k2) = reshape (p% cells% neighbor_idx(j1:j2, 1, 1:p%cell_type), &
                                                   [k2-k1+1]              )
    i2(k1:k2) = reshape (p% cells% neighbor_blk(j1:j2, 1, 1:p%cell_type), &
                                                   [k2-k1+1]              )
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc% ilim1(pe)
       j2 = icon_dc% ilim1(pe+1) - 1
       k1 = (j1-1)*p%cell_type+1
       k2 =  j2   *p%cell_type
       call p_ibcast (i1(k1:k2), pe)
       call p_ibcast (i2(k1:k2), pe)
    end do
    call p_waitall ()
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc% ilim1(pe)
       j2 = icon_dc% ilim1(pe+1) - 1
       k1 = (j1-1)*p%cell_type+1
       k2 =  j2   *p%cell_type
       p% cells% neighbor_idx(j1:j2, 1, 1:p%cell_type) = &
            reshape (i1(k1:k2), [j2-j1+1, p%cell_type])
       p% cells% neighbor_blk(j1:j2, 1, 1:p%cell_type) = &
            reshape (i2(k1:k2), [j2-j1+1, p%cell_type])
    end do
    deallocate (i1,i2)

    !---------
    ! Vertices
    !---------
    allocate (d1(n_verts), d2(n_verts))
    d1 = -HUGE(0._wp)
    d2 = -HUGE(0._wp)

    !-------------------
    ! Vertex coordinates
    !-------------------
    j1 = icon_dc_d% ilim1(dace% pe)
    j2 = icon_dc_d% ilim1(dace% pe+1) - 1
    do j = j1, j2
       jc = marr_v(2,j)
       jb = marr_v(3,j)
       d1(j) = patch% verts% vertex(jc,jb)% lon  ! [rad]
       d2(j) = patch% verts% vertex(jc,jb)% lat  ! [rad]
    end do
    do pe = 0, icon_dc_d% nproc1-1
       j1 = icon_dc_d% ilim1(pe)
       j2 = icon_dc_d% ilim1(pe+1) - 1
       call p_ibcast (d1(j1:j2), pe)
       call p_ibcast (d2(j1:j2), pe)
    end do
    call p_waitall ()
    do j = 1, n_verts
       p% verts% vertex(j,1)% lon = d1(j)
       p% verts% vertex(j,1)% lat = d2(j)
    end do
    deallocate (d1, d2)

    !------------------
    ! Neighbor vertices
    !------------------
    p% verts% num_edges   (:,:)   = -1
    p% verts% neighbor_idx(:,:,:) = -1
    p% verts% neighbor_blk(:,:,:) = -1
    j1 = icon_dc_d% ilim1(dace% pe)
    j2 = icon_dc_d% ilim1(dace% pe+1) - 1
    do j = j1, j2
       jc = marr_v(2,j)
       jb = marr_v(3,j)
       if (.not. patch% verts% decomp_info% owner_mask(jc,jb)) then
          write(1000+dace% pe,*) "##",j,"something really bad happened for",jc,jb,&
               "real owner:", owner_v(jc,1,jb)
       end if
       ne = patch% verts% num_edges(jc,jb)
       if (ne < 6) then
          !--------------------------------
          ! Remove duplicate neighbor cells
          !--------------------------------
          do k = 2, ne
             if (any (      (patch% verts% cell_idx(jc,jb,  k  ) == &
                             patch% verts% cell_idx(jc,jb,1:k-1)    ) &
                      .and. (patch% verts% cell_blk(jc,jb,  k  ) == &
                             patch% verts% cell_blk(jc,jb,1:k-1))   ) ) then
                if (dbg_level > 1) then
                  write(ou,*) "### 'Duplicate' neighbor cells: ne,jc,jb,k=", ne, jc, jb, k
                  write(ou,*) "# patch% verts% cell_idx=", patch% verts% cell_idx(jc,jb,1:k)
                  write(ou,*) "# patch% verts% cell_blk=", patch% verts% cell_blk(jc,jb,1:k)
                end if
                !--------------------------------------------------------------------
                ! Adjust number of edges for boundary points with duplicate neighbors
                !--------------------------------------------------------------------
                ne = k - 1
                exit
             end if
          end do
       end if
       p% verts% num_edges(j,1) = ne
       do k = 1, ne
          nc = patch% verts% neighbor_idx(jc,jb,k)
          nb = patch% verts% neighbor_blk(jc,jb,k)
          if (nc <= 0) then
             cycle  ! k loop
          end if
          pe = dace% pe
          if (.not. patch% verts% decomp_info% owner_mask(nc,nb)) then
             pe = owner_v(nc,1,nb)
             kc = owner_v(nc,2,nb)
             kb = owner_v(nc,3,nb)
             if (pe == dace% pe) then
                write(1000+dace% pe,*) "##",j,k,"jc,jb=",jc,jb,"bad vert_neighbor?",&
                     patch% cells% decomp_info% owner_mask(jc,jb),nc,nb, &
                     patch% cells% decomp_info% owner_mask(nc,nb),pe,kc,kb
             else
                nc = kc
                nb = kb
             end if
          end if
          if (nc > 0 .and. nb > 0) then
             p% verts% neighbor_idx(j,1,k) = vert_glb_idx(nc,nb,pe)
             p% verts% neighbor_blk(j,1,k) = 1
          end if
          if (p% verts% neighbor_idx(j,1,k) < 1) then
             kc = patch% verts% neighbor_idx(jc,jb,k)
             kb = patch% verts% neighbor_blk(jc,jb,k)
             write(1000+dace% pe,*) "##",j,k,":vert_neighbor:",nc,nb,pe,&
                  p% verts% neighbor_idx(j,1,k),jc,jb,owner_v(kc,1:3,kb)
          end if
       end do
    end do
    !--------------------------------
    ! Communicate p% verts% num_edges
    !--------------------------------
    allocate (i3(n_verts))
    i3 = -1
    j1 = icon_dc_d% ilim1(dace% pe)
    j2 = icon_dc_d% ilim1(dace% pe+1) - 1
    i3(j1:j2) = p% verts% num_edges(j1:j2,1)
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc_d% ilim1(pe)
       j2 = icon_dc_d% ilim1(pe+1) - 1
       call p_ibcast (i3(j1:j2), pe)
    end do
    call p_waitall ()
    p% verts% num_edges(:,1) = i3(:)
    deallocate (i3)
    !---------------------------------
    ! Communicate p% verts% neighbor_*
    !---------------------------------
    allocate (i1(n_verts*6), i2(n_verts*6))
    i1 = -1
    i2 = -1
    j1 = icon_dc_d% ilim1(dace% pe)
    j2 = icon_dc_d% ilim1(dace% pe+1) - 1
    k1 = (j1-1)*6+1
    k2 =  j2   *6
    i1(k1:k2) = reshape (p% verts% neighbor_idx(j1:j2, 1, 1:6), [k2-k1+1])
    i2(k1:k2) = reshape (p% verts% neighbor_blk(j1:j2, 1, 1:6), [k2-k1+1])
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc_d% ilim1(pe)
       j2 = icon_dc_d% ilim1(pe+1) - 1
       k1 = (j1-1)*6+1
       k2 =  j2   *6
       call p_ibcast (i1(k1:k2), pe)
       call p_ibcast (i2(k1:k2), pe)
    end do
    call p_waitall ()
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc_d% ilim1(pe)
       j2 = icon_dc_d% ilim1(pe+1) - 1
       k1 = (j1-1)*6+1
       k2 =  j2   *6
       p% verts% neighbor_idx(j1:j2, 1, 1:6) = reshape (i1(k1:k2), [j2-j1+1, 6])
       p% verts% neighbor_blk(j1:j2, 1, 1:6) = reshape (i2(k1:k2), [j2-j1+1, 6])
    end do
!   deallocate (i1,i2)

    !---------------------------
    ! Neighbor cells of vertices
    !---------------------------
    p% verts% cell_idx(:,:,:) = -1
    p% verts% cell_blk(:,:,:) = -1
    j1 = icon_dc_d% ilim1(dace% pe)
    j2 = icon_dc_d% ilim1(dace% pe+1) - 1
    do j = j1, j2
       jc = marr_v(2,j)
       jb = marr_v(3,j)
!       ne = patch% verts% num_edges(jc,jb)
       ne = p% verts% num_edges(j,1)            ! "Adjusted" number of edges
       do k = 1, ne
          nc = patch% verts% cell_idx(jc,jb,k)
          nb = patch% verts% cell_blk(jc,jb,k)
          if (nc <= 0) then
             cycle  ! k loop
          end if
          pe = dace% pe
          if (.not. patch% cells% decomp_info% owner_mask(nc,nb)) then
             pe = owner_c(nc,1,nb)
             kc = owner_c(nc,2,nb)
             kb = owner_c(nc,3,nb)
             if (pe == dace% pe) then
                write(1000+dace% pe,*) "##",j,k,"jc,jb=",jc,jb,"bad vert_cell?",&
                     patch% verts% decomp_info% owner_mask(jc,jb),nc,nb, &
                     patch% cells% decomp_info% owner_mask(nc,nb),pe,kc,kb
             else
                nc = kc
                nb = kb
             end if
          end if
          if (nc > 0 .and. nb > 0) then
             p% verts% cell_idx(j,1,k) = cell_glb_idx(nc,nb,pe)
             p% verts% cell_blk(j,1,k) = 1
          end if
          if (p% verts% cell_idx(j,1,k) < 1) then
             kc = patch% verts% cell_idx(jc,jb,k)
             kb = patch% verts% cell_blk(jc,jb,k)
             write(1000+dace% pe,*) "##",j,k,":vert_cell:",nc,nb,pe,&
                  p% verts% cell_idx(j,1,k),jc,jb,owner_c(kc,1:3,kb)
          end if
       end do
    end do
    !-----------------------------
    ! Communicate p% verts% cell_*
    !-----------------------------
!   allocate (i1(n_verts*6), i2(n_verts*6))
    i1 = -1
    i2 = -1
    j1 = icon_dc_d% ilim1(dace% pe)
    j2 = icon_dc_d% ilim1(dace% pe+1) - 1
    k1 = (j1-1)*6+1
    k2 =  j2   *6
    i1(k1:k2) = reshape (p% verts% cell_idx(j1:j2, 1, 1:6), [k2-k1+1])
    i2(k1:k2) = reshape (p% verts% cell_blk(j1:j2, 1, 1:6), [k2-k1+1])
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc_d% ilim1(pe)
       j2 = icon_dc_d% ilim1(pe+1) - 1
       k1 = (j1-1)*6+1
       k2 =  j2   *6
       call p_ibcast (i1(k1:k2), pe)
       call p_ibcast (i2(k1:k2), pe)
    end do
    call p_waitall ()
    do pe = 0, icon_dc% nproc1-1
       j1 = icon_dc_d% ilim1(pe)
       j2 = icon_dc_d% ilim1(pe+1) - 1
       k1 = (j1-1)*6+1
       k2 =  j2   *6
       p% verts% cell_idx(j1:j2, 1, 1:6) = reshape (i1(k1:k2), [j2-j1+1, 6])
       p% verts% cell_blk(j1:j2, 1, 1:6) = reshape (i2(k1:k2), [j2-j1+1, 6])
    end do
    deallocate (i1,i2)

    !--------------------------------------------------------------
    ! Set up coarser auxiliary grid for searching and interpolation
    !--------------------------------------------------------------
    nullify  (icongrid% uvec_cell, icongrid% vert,  &
              icongrid% uvec_vert, icongrid% aux_grd)
    allocate (icongrid% g)
    call set_search_grid (icongrid, global)

!   if (dace%lpio .and. dbg_level > 0) then
!      write(*,*) "### auxgrid:nx,ny,global=", icongrid% g% nx, icongrid% g% ny, icongrid% g% global
!      write(*,*) "### auxgrid:min/max(lon)=", icongrid% g% lon
!      write(*,*) "### auxgrid:min/max(lat)=", icongrid% g% lat
!      write(*,*) "### auxgrid:lo1,dxi=", icongrid% g% lo1, icongrid% g% dxi
!      write(*,*) "### auxgrid:la1,dyi=", icongrid% g% la1, icongrid% g% dyi
!      write(*,*)
!   end if

    if (dace% lpio) write(ou,*) 'checking neighbor relationships in ICON grid'
    call check_neighbours (icongrid)

    if (dace% lpio) write(ou,*) 'deriving gridpoint distance to grid boundary'
    call dist_to_bound (p, 50)

    !-----------------------------------------------------
    ! Test with model grid points (cell centers, vertices)
    !-----------------------------------------------------
    if (dbg_level > 1) then
#if defined(__PGI) && __PGIC__ <= 19
      ! The compiler does not support BLOCK constructs
#else
       block
          integer  :: idx(3,n_cells), blk(3,n_cells), ngp(n_cells)
          real(wp) :: clon(n_cells), clat(n_cells), w(3,n_cells)
          real(wp) :: vlon(n_verts), vlat(n_verts)
          real(wp) :: t1, t2, d, dmax
          integer  :: j0

          do j = 1, n_verts
             vlon(j) = p% verts% vertex(j,1)% lon * r2d
             vlat(j) = p% verts% vertex(j,1)% lat * r2d
          end do
          if (dace% lpio) then
            write(ou,*)
            write(ou,*) "checking DACE interpolation routines on vertices"
            write(ou,*) "vert: min/max(lon) =", minval (vlon), maxval (vlon)
            write(ou,*) "vert: min/max(lat) =", minval (vlat), maxval (vlat)
          end if
          idx = 0
          blk = 0
          w   = 0
          ngp = 0
          call cpu_time (t1)
          call search_icon_global (icongrid, vlon, vlat, idx, blk, w, ngp)
          call cpu_time (t2)
          if (dace% lpio) then
            write(ou,*) "min/max(ngp) =", minval (ngp(1:n_verts)), maxval (ngp(1:n_verts))
            write(ou,*) "min/max(w)   =", minval (w(:,1:n_verts)), maxval (w(:,1:n_verts))
            write(ou,*) "Time for search [s]:", real (t2-t1)
            write(ou,*) "Time for search [microseconds/point]:", real ((t2-t1)/n_verts*1.e6_wp)
            write(ou,*)
          end if
          !---------------------------------------
          ! Check validity of interpolation points
          !---------------------------------------
          dmax = 0._wp
          j0   = 0
          do j = 1, n_verts
             if (ngp(j) == 3) then
                d = minval ((p% cells% center(idx(:,j),1)% lon * r2d - vlon(j))**2 &
                           +(p% cells% center(idx(:,j),1)% lat * r2d - vlat(j))**2 )
                d = sqrt (d)
                if (d > dmax) then
                   dmax = d
                   j0   = j
                end if
             end if
          end do
          if (dace% lpio) then
            write(ou,*) "Maximum distance to closest cell center", real (dmax), &
                 "for index", j0, ": vlon,vlat=", real([vlon(j0), vlat(j0)])
            write(ou,*)
          end if
          ! Derive points from cell center coordinates with small perturbations
          call random_number (clon)
          call random_number (clat)
          do j = 1, n_cells
             clon(j) = (2*clon(j)-1)*0.005_wp + p% cells% center(j,1)% lon * r2d
             clat(j) = (2*clat(j)-1)*0.005_wp + p% cells% center(j,1)% lat * r2d
          end do
          if (dace% lpio) then
            write(ou,*) "checking DACE interpolation routines near cell centers"
            write(ou,*) "cell: min/max(lon) =", minval (clon), maxval (clon)
            write(ou,*) "cell: min/max(lat) =", minval (clat), maxval (clat)
          end if
          idx = 0
          blk = 0
          w   = 0
          ngp = 0
          call cpu_time (t1)
          call search_icon_global (icongrid, clon, clat, idx, blk, w, ngp)
          call cpu_time (t2)
          if (dace% lpio) then
            write(ou,*) "min/max(ngp) =", minval (ngp(1:n_cells)), maxval (ngp(1:n_cells))
            write(ou,*) "min/max(w)   =", minval (w(:,1:n_cells)), maxval (w(:,1:n_cells))
            write(ou,*) "Time for search [s]:", real (t2-t1)
            write(ou,*) "Time for search [microseconds/point]:", real ((t2-t1)/n_cells*1.e6_wp)
            write(ou,*)
          end if
          !---------------------------------------
          ! Check validity of interpolation points
          !---------------------------------------
          dmax = 0._wp
          j0   = 0
          do j = 1, n_cells
             if (ngp(j) == 3) then
                d = minval ((p% cells% center(idx(:,j),1)% lon * r2d - clon(j))**2 &
                           +(p% cells% center(idx(:,j),1)% lat * r2d - clat(j))**2 )
                d = sqrt (d)
                if (d > dmax) then
                   dmax = d
                   j0   = j
                end if
             end if
          end do
          if (dace% lpio) then
            write(ou,*) "Maximum distance to closest cell center", real (dmax), &
                 "for index", j0, ": vlon,vlat=", real([clon(j0), clat(j0)])
            write(ou,*)
          end if
       end block
#endif
    end if

  end subroutine icongrid_from_icon
  !----------------------------------------------------------------------------
  subroutine grid_from_icon (g)
    type(t_grid), intent(inout) :: g

    type(t_patch_dace), pointer :: p            ! Patch data type, DACE
    integer                     :: j, k, pe     ! Loop indices
    integer                     :: j1, j2       ! Temporary bounds
    integer                     :: idx, blk     ! Temporary indices
    integer                     :: nz           ! Number of full levels
    integer                     :: ierr
    real(wp)                    :: frmin, frmax ! minval, maxval of fr_land

    if (g% icongrid% grid_root > 0 .and. g% icongrid% grid_level > 0) then
       g% ni = g% icongrid% grid_root * 2 ** g% icongrid% grid_level
    end if
    if (g% ni > 0) then
       g% d_km  = 7054._wp / (g% ni * sqrt(2._wp))
       g% d_deg = r2d * g% d_km / (0.001_wp * g% a)
    end if

    g% global =  g% icongrid% global
    p         => g% icongrid% patch

    call set_grid_fields (g)
    call update (g%m)
    call set_pointers (g)
    !-----------------------------------------------------
    ! Set up arrays xnglob, rlon, rlat for given ICON grid
    !-----------------------------------------------------
    allocate (g% dlon (g% nx))
    allocate (g% dlat (g% ny))
    g% dlon(:)  =  p% cells% center(:,1)% lon * r2d
    g% dlat(:)  =  p% cells% center(:,1)% lat * r2d
    call setup_global_coord (g% xnglob, g% rlon, g% rlat,            &
                             reshape (g% dlon * d2r,(/g% nxny,1,1/)),&
                             reshape (g% dlat * d2r,(/g% nxny,1,1/)),&
                             g% lbg, g% ubg, ierr                    )

    nz = g% nz
    !------------------------
    ! Get HHL, hsurf, fr_land
    !------------------------
    call allocate (g, 'hsurf')    ! non-decomposed
    call allocate (g, 'lsm')
    call allocate (g, 'soiltyp')
    call allocate (g, 'fr_lake')
    call allocate (g, 'depth_lk')
    call allocate (g, 'hhl')      ! decomposed
    call allocate (g, 'sso_stdh')
    j1 = icon_dc% ilim1(dace% pe)
    j2 = icon_dc% ilim1(dace% pe+1) - 1
    do j = j1, j2
       if (marr_c(1,j) /= dace% pe) call finish ("grid_from_icon","bad marr_c")
       idx = marr_c(2,j)
       blk = marr_c(3,j)
       do k = 1, nz+1
          g% hhl  (j,1,k,1) = p_nh_state(1)% metrics% z_ifc(idx,k,blk)
       end do
       g% hsurf   (j,1,1,1) = g% hhl(j,1,nz+1,1)
       g% lsm     (j,1,1,1) = ext_data(1)% atm% fr_land (idx,blk)
       g% soiltyp (j,1,1,1) = ext_data(1)% atm% soiltyp (idx,blk)
       g% fr_lake (j,1,1,1) = ext_data(1)% atm% fr_lake (idx,blk)
       g% depth_lk(j,1,1,1) = ext_data(1)% atm% depth_lk(idx,blk)
       g% sso_stdh(j,1,1,1) = ext_data(1)% atm% sso_stdh(idx,blk)
    end do
    !-------------------------------------------
    ! Communicate components to be held globally
    !-------------------------------------------
    do pe = 0, g% dc% nproc1-1
       j1 = icon_dc% ilim1(pe)
       j2 = icon_dc% ilim1(pe+1) - 1
       call p_ibcast (g% hsurf   (j1:j2,1,1,1), pe)
       call p_ibcast (g% lsm     (j1:j2,1,1,1), pe)
       call p_ibcast (g% soiltyp (j1:j2,1,1,1), pe)
       call p_ibcast (g% fr_lake (j1:j2,1,1,1), pe)
       call p_ibcast (g% depth_lk(j1:j2,1,1,1), pe)
    end do
    call p_waitall ()

    frmin = minval (g% lsm)
    frmax = maxval (g% lsm)
    !---------------------------------------------
    ! workaround: values of fr_land > 1.0 occur,
    ! restrict to [0.0:1.0]
    !---------------------------------------------
    if (frmin < 0._wp) then
       if (dace% lpio) write(0,*) "WARNING: ICON bug: MINVAL(fr_land) =", frmin
       g% lsm = max (g% lsm, 0.0_wp)
    end if
    if (frmax > 1._wp) then
       if (dace% lpio) write(0,*) "WARNING: ICON bug: MAXVAL(fr_land) =", frmax
       g% lsm = min (g% lsm, 1.0_wp)
    end if
    flush (0)

    !-------------------
    ! Derived quantities
    !-------------------
    call allocate (g, 'geosp')
    g% geosp = g% hsurf * gacc

    !-------------------------------
    ! geoid correction for GNSS data
    !-------------------------------
    call set_geoid (g)
    !-------------------------------------
    ! Estimate topmost full level pressure
    ! Set topmost height (gpm)
    !-------------------------------------
    call set_ptopf (g)
    !-----------------------------------------
    ! estimate unperturbed model-level heights
    !-----------------------------------------
    call set_zlev_ref (g)
    !--------------------------------
    ! set some pressure level indices
    !--------------------------------
    call set_plev_indices (g)

    if (dbg_level > 0) then
       call message ("grid_from_icon","ICON grid:")
       call print   (g, iunit=dace_output_unit, verbose=dbg_level > 1)
       call message ("","")
    end if

  end subroutine grid_from_icon
  !============================================================================
  subroutine atm_from_icon (state)
    type(t_atm), intent(inout), target :: state

    type(t_grid),         pointer :: g
    integer                       :: j, k     ! Loop indices
    integer                       :: j1, j2   ! Temporary bounds
    integer                       :: idx, blk ! Temporary indices
    integer                       :: nz       ! Number of full levels
    logical                       :: lqr, lqs ! Rain, snow
    logical                       :: lqg      ! Graupel
    type(t_nh_prog),      pointer :: atm_p    ! atmosphere, prognostic vars.
    type(t_nh_prog),      pointer :: atm_r    ! dto., reduced calling frequency
    type(t_nh_diag),      pointer :: atm_d    ! atmosphere, diagnostic vars.
    type(t_nwp_phy_diag), pointer :: phy_d    ! physical model diagnostic vars.
    type(t_lnd_prog),     pointer :: lnd_p    ! land model, prognostic vars.
    type(t_lnd_diag),     pointer :: lnd_d    ! land model, diagnostic vars.
    logical                       :: l_rad_cld
    ! character(*), parameter :: fields = &
    !      "ps pf ph t u v den q qcl qci qv_s z0 qv_dia qc_dia qi_dia& 
    !      & t2m td2m rh2m u_10m v_10m clct clcl clcm clch clc&
    !      & tsurf h_snow fr_ice" ! "t_so" currently not used
    character(*), parameter :: fields_default = &
         "ps pf ph t u v den q qcl qci qv_s z0& 
         & t2m td2m rh2m u_10m v_10m clct clcl clcm clch&
         & tsurf h_snow fr_ice" ! "t_so" currently not used
    character(*), parameter :: fields_rad_cld      = &
                                                  "qv_dia qc_dia qi_dia clc"
    character(*), parameter :: fields_rad_reff  = "reff_qc reff_qi" 
    character(512)          :: fields
    integer                 :: ou

    ! Pointer to the acutally used variant of clc:
    REAL(wp), DIMENSION(:,:,:), POINTER ::  ptr_clc

    g  => state% grid
    nz =  g% nz
    atm_d => p_nh_state(1)% diag
    atm_p => p_nh_state(1)% prog(nnow    (1))
    atm_r => p_nh_state(1)% prog(nnow_rcf(1))
    phy_d => prm_diag  (1)
    lnd_d => p_lnd_state(1)% diag_lnd
    lnd_p => p_lnd_state(1)% prog_lnd(nnow_rcf(1))

    l_rad_cld = .false.
    do j = 1, n_set
      if (any(rad_set(j)%iopts(1:rad_set(j)%n_instr)%cloud_mode > 0)) then
        l_rad_cld = .true.
      end if
    end do

    fields = fields_default
    if (l_rad_cld) fields = trim(fields)//' '//trim(fields_rad_cld)
    if (l_rad_cld .and. use_reff .and. atm_phy_nwp_config(1)% icalc_reff .gt. 0) &
       fields = trim(fields)//' '//trim(fields_rad_reff)

    IF (l_rad_cld) THEN
      NULLIFY(ptr_clc)
      ! Decide which field for cloud cover has to be used:
      IF (atm_phy_nwp_config(1)%luse_clc_rad) THEN
        ptr_clc => phy_d % clc_rad(:,:,:)
      ELSE
        ptr_clc => phy_d % clc(:,:,:)
      END IF
    END IF
    
    call allocate (state, fields)
    
    ! Ensure that diagnostic fields are up-to-date (HR, CW)
    call rbf_vec_interpol_cell (atm_p%vn, p_patch(1), p_int_state(1), atm_d%u, atm_d%v)

    call diagnose_pres_temp ( p_nh_state(1)%metrics, atm_p, atm_r, atm_d, p_patch(1), &
         &                    opt_calc_temp=.TRUE.,                                   &
         &                    opt_calc_pres=.TRUE. )

    if (iqv < 1) call finish ("atm_from_icon","no qv!")
    if (iqc < 1) call finish ("atm_from_icon","no qc!")
    if (iqi < 1) call finish ("atm_from_icon","no qi!")
    lqr = iqr > 0
    lqs = iqs > 0
    lqg = iqg > 0
    if (lqr) call allocate (state, "qr")
    if (lqs) call allocate (state, "qs")
    if (lqg) call allocate (state, "qg")

    if (dbg_level > 1) then
       if (dace% lpio) write(0,*) "iqv,iqc,iqi,iqr,iqs,iqg=",iqv,iqc,iqi,iqr,iqs,iqg
    end if

    j1 = icon_dc% ilim1(dace% pe)
    j2 = icon_dc% ilim1(dace% pe+1) - 1

    do j = j1, j2
       if (marr_c(1,j) /= dace% pe) call finish ("atm_from_icon","bad marr_c")
       idx = marr_c(2,j)
       blk = marr_c(3,j)
       state% ps    (j,1,1,1) = atm_d% pres_sfc(idx,blk)
       do k = 1, nz
          state% pf (j,1,k,1) = atm_d% pres     (idx,k,blk)
          state% t  (j,1,k,1) = atm_d% temp     (idx,k,blk)
          state% u  (j,1,k,1) = atm_d% u        (idx,k,blk)
          state% v  (j,1,k,1) = atm_d% v        (idx,k,blk)
       end do
       do k = 1, nz+1
          state% ph (j,1,k,1) = atm_d% pres_ifc (idx,k,blk)
       end do

       do k = 1, nz
          state% den  (j,1,k,1) = atm_p% rho    (idx,k,blk)
       end do
       
       do k = 1, nz
          state%          q      (j,1,k,1) = atm_r% tracer(idx,k,blk,iqv)
          state%          qcl    (j,1,k,1) = atm_r% tracer(idx,k,blk,iqc)
          state%          qci    (j,1,k,1) = atm_r% tracer(idx,k,blk,iqi)
          if (l_rad_cld) then
            state%          clc    (j,1,k,1) = ptr_clc         (idx,k,blk)*100._wp !convert to percent
            state%          qv_dia (j,1,k,1) = phy_d% tot_cld  (idx,k,blk,iqv)
            state%          qc_dia (j,1,k,1) = phy_d% tot_cld  (idx,k,blk,iqc)
            state%          qi_dia (j,1,k,1) = phy_d% tot_cld  (idx,k,blk,iqi)

            if (atm_phy_nwp_config(1)% icalc_reff .gt. 0 .and. use_reff ) then 
               state%       reff_qc(j,1,k,1) = phy_d% reff_qc (idx,k,blk)    
               state%       reff_qi(j,1,k,1) = phy_d% reff_qi (idx,k,blk) 
            end if
          end if
          if (lqr) state% qr (j,1,k,1) = atm_r% tracer(idx,k,blk,iqr)
          if (lqs) state% qs (j,1,k,1) = atm_r% tracer(idx,k,blk,iqs)
          if (lqg) state% qg (j,1,k,1) = atm_r% tracer(idx,k,blk,iqg)
       end do

       state% z0     (j,1,1,1) = phy_d% gz0      (idx,blk) / gacc
       state% t2m    (j,1,1,1) = phy_d% t_2m     (idx,blk)
       state% td2m   (j,1,1,1) = phy_d% td_2m    (idx,blk)
       state% rh2m   (j,1,1,1) = phy_d% rh_2m    (idx,blk)
       state% u_10m  (j,1,1,1) = phy_d% u_10m    (idx,blk)
       state% v_10m  (j,1,1,1) = phy_d% v_10m    (idx,blk)
       state% clct   (j,1,1,1) = phy_d% clct     (idx,blk)*100._wp !convert to percent
       state% clcl   (j,1,1,1) = phy_d% clcl     (idx,blk)*100._wp !convert to percent
       state% clcm   (j,1,1,1) = phy_d% clcm     (idx,blk)*100._wp !convert to percent
       state% clch   (j,1,1,1) = phy_d% clch     (idx,blk)*100._wp !convert to percent
       state% tsurf  (j,1,1,1) = lnd_p% t_g      (idx,blk)
       state% h_snow (j,1,1,1) = lnd_d% h_snow   (idx,blk)
       state% fr_ice (j,1,1,1) = lnd_d% fr_seaice(idx,blk)
       state% qv_s   (j,1,1,1) = lnd_d% qv_s     (idx,blk)
!      state% t_so   (j,1,:,1) = lnd_p% t_so   (idx,:,blk) ! needs checking
    end do

    call set_geo (state, geof=.true.)
    call set_rh  (state)

    if (dbg_level > 0) then
       ou = dace_output_unit
       call print (state, grid=.false., comment="ICON state", &
                          iunit=ou,     verbose= dbg_level > 1)
       flush (ou)
    end if

  end subroutine atm_from_icon
  !============================================================================
  subroutine set_refatm (atm, time)
    type(t_atm),  intent(inout) :: atm
    type(t_time), intent(in)    :: time

    integer             :: lb  (4)    ! lower bounds
    integer             :: ub  (4)    ! upper bounds
    integer             :: i,j,k,d    ! loop indices
!   integer             :: nz
    real(wp)            :: z, zs
    real(wp), parameter :: P_00 = 101325._wp
    real(wp), parameter :: T_00 =    213.15_wp ! -60 degree Celsius
    real(wp), parameter :: DT   =     75._wp
    real(wp), parameter :: H_0  =  10000._wp   ! 10 km

    call allocate (atm, "ps pf t")
    atm = -HUGE (0._wp)

    call set_geo  (atm, geof=.true.)

    lb = atm% grid% lb
    ub = atm% grid% ub
!   nz = ub(3)
    do       d = lb(4), ub(4)
       do    j = lb(2), ub(2)
          do i = lb(1), ub(1)
             zs = atm% grid% hsurf(i,j,1,d)
             atm% ps   (i,j,1,d) = pref (zs)
             do k = lb(3), ub(3)
                z = atm% geof (i,j,k,d) / gacc
                atm% t (i,j,k,d) = tref (z)
                atm% pf(i,j,k,d) = pref (z)
             end do
          end do
       end do
    end do

    atm% ref_time = time
    atm% time     = time

  contains

    elemental function tref (z)
      real(wp), intent(in) :: z
      real(wp)             :: tref
      tref = T_00 + DT * exp (-z/H_0)
    end function tref

    elemental function pref (z)
      real(wp), intent(in) :: z
      real(wp)             :: pref
      real(wp), parameter  :: gamma = gacc*H_0/(Rd*T_00)
      pref = P_00 * exp (-gamma * log ( (T_00*exp (z/H_0) + DT) / (T_00+DT) ) )
    end function pref

  end subroutine set_refatm

#endif /* __DACE__ */
  !============================================================================
  subroutine init_dace_op ()
    !-----------------------------------------------------------------------
    ! initialise the DACE observation operators from ICON
    ! 1) set grid according to the specifications by ICON
    ! 2) read relevant namelists in DACE
    ! 3) read observation CDFIN files, set up operators
    ! 4) distribute observations over PEs, set up interpolation coefficients
    ! 5) return information on required time steps to ICON
    !-----------------------------------------------------------------------
#ifndef __DACE__
    CALL finish ("init_dace_op","DACE coupling requested but not compiled in")
#else
    type(t_atm), pointer :: refatm => NULL()
    type(t_obs)          :: sdobs          ! send buffer
    type(t_obs)          :: rcobs          ! recv buffer
    integer              :: i              ! loop index
    integer              :: ntri           ! number of boxes
    integer     ,pointer :: pes (:)        ! processors / box
    type(t_time)         :: last_time      ! time of last slot
    ! Testing:
    integer :: n_slot

    flush (6)
    call message ("init_dace_op","")
    if (.not. associated (grid)) call finish ("init_dace_op","grid not set!")
    if (grid% dc% myproc /= dace% pe) then
       write(0,*) "pe, myproc =", dace% pe, grid% dc% myproc
       call finish ("init_dace_op","bad myproc")
    end if

    call dace_open_output ()

    !----------------------------------
    ! "Reference" atmosphere for MEC:
    ! get current model state from ICON
    !----------------------------------
    call message ("init_dace_op","setting refatm from ICON")
    allocate (refatm)
    call construct  (refatm, grid=grid)
    refatm% runtype     = "forecast"
    refatm% runclass    = 2               ! 0=haupt, 2=ass
    refatm% expid       = exp_id
    refatm% member      = ens_mem
    refatm% members     = num_mem         ! Ensemble size
    refatm% ensemble_id = ens_id
    refatm% ref_time    = exp_start       ! as needed by read_cosmo_obs
    refatm% time        = exp_end         ! as needed by read_cosmo_obs
    call atm_from_icon (state=refatm)
    call message ("init_dace_op","get refatm from ICON done.")

    call message ("init_dace_op","reading DACE remaining namelists")
    !=================================================
    ! open namelist file (unit: nnml)
    ! read namelists:
    !   /OBSERVATIONS/
    !   /.../
    !=================================================
    if (dace% lpio) call open_nml (assimilation_config(1)%dace_namelist_file)

    call read_nml_mec_obs () ! read namelist /MEC_OBS/
    call read_obs_nml     () ! read namelist /observations/
    call read_std_nml_dace() ! read namelist /STD_OBS/
    call init_fdbk_tables () ! initialise tables
    call disable_gh       () ! disable generalized humidity transformation
    call read_nml_report     ! set defaults in table 'rept_use'
    call read_tovs_nml       ! read namelists /TOVS_OBS/ and /TOVS_OBS_CHAN_NML/
    call read_nml_thin       ! read namelist /THINNING/
    if (n_dace_op > 0) then
       call read_blacklists  ! read namelist /BLACKLIST/ (and blacklist file)
    end if
    flush (6)


    if  (use_reff .and. atm_phy_nwp_config(1)% icalc_reff <= 0) &
         call finish('init_dace_op', 'use_reff (DACE namelist) requires &
         &icalc_reff > 0 (ICON namelist)')

    !===================================
    ! read observations from CDFIN files
    ! (see read_veri_obs for details).
    !===================================
    read_cdfin = .true.
    !-------------------------------
    ! initialize observation modules
    !-------------------------------
    call message ('init_dace_op','initialize observation modules')
    call set_rules_cosmo
    call read_nml_rules
    call construct (obs_in(1))
    obs_b% o => obs_in(1)
    call process_obs (TSK_INIT, obs_b, refatm)
    flush (6)
!   call stop_time ('scan observations')
    call message ('init_dace_op','scan observations')
    !------------------
    ! read observations
    !------------------
    call set_input_boxes
    flush (6)
!   call stop_time ('read observations')
    call message ('init_dace_op','read observations')
    call process_obs (TSK_READ, obs_b, refatm)
    call flush_buf
    call unique_spot (obs_b% o)
    !-----------------------------------
    ! check for wanted observation types
    ! release memory not used
    !-----------------------------------
    call check_obstype(obs_in, obstypes) ! check if obstype is in list
    call release_mem (obs_b% o, keep = obs_b% o% spot(:)%o%n > 0)
    !------------------------------------------------------
    ! check for out of domain observations
    ! + for COSMO operators do reject,
    ! + dismiss does not work yet
    !------------------------------------------------------
    ! check if problems with global conventional operators!
    ! HA: lecheck_ptopf noetig (wie in dace)?
    ! Use status rejected only for COSMO operators; not for all observation types
    where(rept_char(:)% mod == COSMO) rept_use(:)% use(CHK_DOMAIN) = STAT_REJECTED
         
    call check_domain (obs_in(1), grid, horizontal=.true.)
    !lcheck_ptopf=lcheck_ptopf)

    obs_in% pe = dace% pe
    !call stop_time ('thinning routine')
    call message ('init_dace_op','thinning routine')
    call thinning (obs_in, grid, 0)
    !------------------------------------------
    ! Set up boxes to run observation operators
    !------------------------------------------
!   call stop_time ('set up boxes')
    call message ('init_dace_op','set up boxes')
    !obs_in% pe = dace% pe
    tinc       = t_time(0,dace_time_ctrl(3))
    mec_slot1  = 1
    mec_slotn  = max (mec_calls, 1)
    n_slot     = mec_slotn -  mec_slot1 + 1
    last_time  = exp_start + (mec_slotn - 1) * tinc
    call set_time_slice (obs_in, last_time, tinc, n_slot, interpolation)
    call set_veri_boxes (obs_in, pes, n_slot, interpolation, obs_local, refatm)
    ntri = size (pes)

    !---------------------------------
    ! distribute observations over PEs
    !---------------------------------
!   call stop_time ('broadcast observations to PEs')
    call message ('init_dace_op','broadcast observations to PEs')
    call process_obs (TSK_SHRINK, obs_b, refatm, state=STAT_PAS_REJ)
    call test_obs    (obs_in,'TSK_SHRINK',0)
    allocate         (obs% o (ntri))
    call test_obs    (obs_in(1),0,'obs_in',0)
    call reorder_obs (obs_in(1), sdobs, obs_in(1)% spot% fgbx)
    call test_obs    (sdobs,0,'sdobs',0)
    call p_alltoall  (sdobs, rcobs, sdobs% spot% fgbx)
    call destruct    (sdobs)
    call scatter_obs (rcobs, obs% o, rcobs% spot% fgbx)
    call release_mem (obs% o)
    call destruct    (rcobs)
    obs% o% pe = pes
    do i=1,ntri
      call p_bcast (obs% o(i)% n_spot, obs% o(i)% pe)
      call p_bcast (obs% o(i)% n_obs,  obs% o(i)% pe)
      call p_bcast (obs% o(i)% n_int,  obs% o(i)% pe)
      call p_bcast (obs% o(i)% n_par,  obs% o(i)% pe)
    end do
    deallocate    (pes)
    call destruct (obs_in(1))

    !------------------------------
    ! write report usage statistics
    !------------------------------
    call derive_rept_stat (obs% o)
    call gather_rept_stat
    call print_rept_stat  ('after observation input', unit=dace_output_unit)

    !------------------------------
    ! apply some checks on the data
    !------------------------------
!   call   check_cons         ! consistency check
!   call   check_gross        ! generic gross check (to be called from MEC)
!   call   check_black        ! check for blacklisting
!   call   check_rule         ! check for specific rules
    call   check_obs (obs% o) ! check for valid report
    call   check_suff(obs% o) ! final check for sufficient data in report

    !------------------------------
    ! write report usage statistics
    !------------------------------
    call derive_rept_stat (obs% o)
    call gather_rept_stat
    call print_rept_stat  ('after additional checks', unit=dace_output_unit)

    !--------
    ! Cleanup
    !--------
    if (dace% lpio) call close_nml ()

    call destruct (refatm)
    deallocate    (refatm)

    dace_op_init = .true.
    call message ("init_dace_op","done.")
    CALL message ('','')
    call dace_close_output ()

#endif /* __DACE__ */

  end subroutine init_dace_op
  !============================================================================
#ifdef __DACE__

  subroutine set_dace_timer ()
    integer                   :: i
    integer                   :: secs           ! Seconds since exp_start
    integer                   :: ierr
    logical                   :: lret
    INTEGER                   :: mec_Events
    type(t_time), allocatable :: mec_time(:)
    TYPE(eventGroup), POINTER :: mec_EventGroup => NULL()
    TYPE(datetime),   POINTER :: mec_StartDate  => NULL()
    TYPE(datetime),   POINTER :: mec_StopDate   => NULL()
    TYPE(timedelta),  POINTER :: mec_Start      => NULL()
    TYPE(timedelta),  POINTER :: mec_Stop       => NULL()
    TYPE(timedelta),  POINTER :: mec_Interval   => NULL()
    TYPE(timedelta),  POINTER :: time_step      => NULL()
    CHARACTER(MAX_TIMEDELTA_STR_LEN)   :: td_string
    CHARACTER(MAX_DATETIME_STR_LEN)    :: dt_string
    CHARACTER(MAX_MTIME_ERROR_STR_LEN) :: errstring

    time_step  => time_config%tc_dt_model
    flush (6)
    CALL message ('','')
    CALL message ("init_dace","initializing timer for MEC")
    dace_time_ctrl =  assimilation_config(1)%  dace_time_ctrl
    mec_Start      => newTimeDelta(secs2string(dace_time_ctrl(1)))
    mec_Stop       => newTimeDelta(secs2string(dace_time_ctrl(2)))
!   mec_Stop       => newTimeDelta(secs2string(dace_time_ctrl(2) + 11)) ! not OK!
!   mec_Stop       => newTimeDelta(secs2string(dace_time_ctrl(2) + 12)) ! OK
    mec_Interval   => newTimeDelta(secs2string(dace_time_ctrl(3)))

    if (dace_time_ctrl(1) /= 0) then
       call finish ("set_dace_timer","dace_time_ctrl(1) /= 0 not supported")
    end if

    if (dace_time_ctrl(1) > dace_time_ctrl(2)) then
       if (dace% lpio) write(0,*) "Invalid dace_time_ctrl:", dace_time_ctrl
       call finish ("set_dace_timer","invalid dace_time_ctrl")
    end if
    if (dace_time_ctrl(3) <= 0) then
       if (dace% lpio) write(0,*) "Invalid time step:", dace_time_ctrl(3)
       call finish ("set_dace_timer","invalid time step for MEC")
    end if
    mec_calls = (dace_time_ctrl(2) - dace_time_ctrl(1)) / dace_time_ctrl(3) + 1

    allocate (mec_time(1:mec_calls))
    do i = 1, mec_calls
       secs        =  dace_time_ctrl(1) + (i-1)*dace_time_ctrl(3)
       mec_time(i) =  exp_start + t_time(0,secs)
    end do

    mec_RefDate    => time_config%tc_exp_startdate
    mec_StartDate  => newDatetime (mec_RefDate)
    mec_StopDate   => newDatetime (mec_RefDate)
    mec_StartDate  =  mec_StartDate + mec_Start
    mec_StopDate   =  mec_StopDate  + mec_Stop

    ! Workaround for mtime bug with plus_slack at end of interval:
    CALL message ("init_dace","enabling workaround for mtime bug with slack")
    if (dace% lpio) then
       CALL datetimeToString (mec_StopDate, dt_string)
       WRITE(0,*) "  stop date before: ",   dt_string
       CALL timedeltaToString (time_step,   td_string)
       WRITE(0,*) "  shift           : ",   td_string
    end if
    mec_StopDate   =  mec_StopDate + time_step
    if (dace% lpio) then
       CALL datetimeToString (mec_StopDate, dt_string)
       WRITE(0,*) "  stop date after : ",   dt_string
    end if
    flush (0)

    mec_Events     =  addEventGroup('mecEventGroup')
    mec_EventGroup => getEventGroup(mec_Events)
    mec_Event      => newEvent('mec', mec_RefDate,  mec_StartDate,          &
                                      mec_StopDate, mec_Interval, errno=ierr)

    if (dace% lpio .and. (dbg_level > 0 .or. ierr /= 0)) then
       CALL timedeltaToString (mec_Start,    td_string)
       WRITE (0,*) "MEC start time    : ",   td_string
       CALL timedeltaToString (mec_Stop,     td_string)
       WRITE (0,*) "MEC stop time     : ",   td_string
       CALL timedeltaToString (mec_Interval, td_string)
       WRITE (0,*) "MEC interval      : ",   td_string
       CALL datetimeToString (mec_RefDate,   dt_string)
       WRITE (0,*) "MEC reference date: ",   dt_string
       CALL datetimeToString (mec_StartDate, dt_string)
       WRITE (0,*) "MEC start date    : ",   dt_string
       CALL datetimeToString (mec_StopDate,  dt_string)
       WRITE (0,*) "MEC stop date     : ",   dt_string
    end if
    if (ierr /= 0) then
       CALL mtime_strerror(ierr, errstring)
       CALL finish ("init_dace","initializing MEC events: "//errstring)
    end if

    lret = addEventToEventGroup(mec_Event, mec_EventGroup)
    if (dace% lpio .and. (.not. lret .or. dbg_level > 0)) then
       WRITE (0,*) "addEventToEventGroup returns:", lret
    end if

    if (dbg_level > 0) call check_dace_timer ()

    CALL message ('','')
    CALL printEventGroup (mec_Events)
    CALL message ('','')

    call deallocateTimedelta (mec_Start)
    call deallocateTimedelta (mec_Stop)
    call deallocateTimedelta (mec_Interval)
    call deallocateDatetime  (mec_StartDate)
    call deallocateDatetime  (mec_StopDate)
  contains
    subroutine check_dace_timer ()
      TYPE(datetime),  POINTER :: mtime
      TYPE(event),     POINTER :: next_Event => NULL()
      integer                  :: i
      integer                  :: secs

      mtime      => newDatetime (mec_RefDate)
      next_Event => newEvent('prep_mec', mec_RefDate,  mec_StartDate,          &
                                         mec_StopDate, mec_Interval, errno=ierr)
      lret = addEventToEventGroup(next_Event, mec_EventGroup)
      if (dace% lpio .and. (.not. lret)) then
         WRITE (0,*) "check_dace_timer: addEventToEventGroup returns:", lret
      end if

      CALL message ('','')
      CALL message ('check_dace_timer','')
      i = 0
      mtime = mec_RefDate
      do
         if (isCurrentEventActive (next_Event, mtime, plus_slack=time_step)) then
            secs = nint (getElapsedSimTimeInSeconds (mtime))
            i    = i + 1
            if (dace% lpio) then
               CALL datetimeToString (mtime, dt_string)
               WRITE (0,*) "MEC will be called on: ", trim (dt_string), &
                    "  (elapsed:", secs, "s)"
            end if
         end if
         IF (mtime >= time_config%tc_stopdate) THEN
            exit
         end IF
         mtime = mtime + time_step
      end do
      if (dace% lpio) write(0,*) "check_dace_timer: total MEC calls:", i, &
           "(expected:", mec_calls, ")"
      if (i /= mec_calls) call finish ("check_dace_timer","internal error")
      CALL message ('','')

      call deallocateDatetime  (mtime)
    end subroutine check_dace_timer
  end subroutine set_dace_timer
  !============================================================================
  function secs2string (secs)
    integer, intent(in) :: secs
    character(len=12)   :: secs2string
    secs2string = ""
    write (secs2string,'("PT",i0,"S")') secs
  end function secs2string
  !============================================================================
  function string2secs (str) result (secs)
    character(*), intent(in) :: str
    integer                  :: secs
    !-----------------------------------------
    ! Convert string from timedelta to seconds
    !-----------------------------------------
    integer             :: i, n, scale, ios
    character(len(str)) :: s
    s = str
    n = len_trim (s)
    if (s(1:2) /= "PT") call finish ("string2secs","bad td_string: "//trim (s))
    select case (s(n:n))
    case ("S")
       scale = 1
    case ("M")
       scale = 60
    case ("H")
       scale = 3600
    case default
       call finish ("string2secs","unsupported: "//trim (s))
    end select
    do i = 3, n-1
       if (s(i:i) == ".") s(i:i) = " "
    end do
    secs = 0
    read (s(3:n-1),*,iostat=ios) secs
    if (ios /= 0) then
       write(0,*) "string2secs: s(3:n-1) = '", s(3:n-1), "'"
       call finish ("string2secs","bad input: "//trim (str))
    end if
    secs = secs * scale
  end function string2secs

#endif /* __DACE__ */
  !============================================================================
  subroutine run_dace_op (mtime_current)
    TYPE(datetime), POINTER :: mtime_current    !< current datetime (mtime)
    !---------------------------------------------------------
    ! run the dace observation operators for a given time step
    !---------------------------------------------------------
#ifndef __DACE__
    CALL finish ("run_dace_op","DACE coupling requested but not compiled in")
#else
    CHARACTER(MAX_TIMEDELTA_STR_LEN) :: td_string
    CHARACTER(MAX_DATETIME_STR_LEN)  :: dt_string
    TYPE(timedelta),  POINTER        :: mec_time  => NULL()
    type(t_atm), pointer  :: x
    integer               :: secs
    logical               :: first = .true.

    if (.not. dace_op_init) then
       call finish ("run_dace_op","DACE operators not initialized.")
    end if

    call dace_open_output ()

    mec_time => newtimedelta('PT0S')
    mec_time =  mtime_current - mec_RefDate
    CALL datetimeToString  (mtime_current, dt_string)
    CALL timedeltaToString (mec_time,      td_string)
    secs = nint (getElapsedSimTimeInSeconds (mtime_current))

    if (dace% lpio .and. dbg_level > 0) then
       WRITE (0,*) "MEC current date : ",  trim (dt_string)
       WRITE (0,*) "MEC relative time: ",  trim (td_string), "  (", secs, "s)"
    end if

    call deallocateTimedelta (mec_time)

    if (mec_sloti == -1) then
       allocate (atm(mec_slot1:mec_slotn))
       mec_sloti = mec_slot1
    else
       mec_sloti = mec_sloti + 1
    end if

    if (mec_sloti > mec_slotn) then
       call finish ("run_dace_op","FATAL: mec_sloti > mec_slotn")
    end if

    !----------------------------------
    ! Get current model state from ICON
    !----------------------------------
    x              => atm(mec_sloti)
    call construct (x, grid=grid)
    x% runtype     = "forecast"
    x% runclass    = 2               ! 0=haupt, 2=ass
    x% expid       = exp_id
    x% member      = ens_mem
    x% members     = num_mem         ! Ensemble size
    x% ensemble_id = ens_id
    x% ref_time    = exp_start
    x% time        = exp_start + t_time(0,secs)
    call atm_from_icon (state=x)

    !--------------------------
    ! run observation operators
    !--------------------------
    !   info:   apply_H in mo_t_enkf.f90
    !     calls process_obs1 -> process_obs0 -> call_specific in mo_obs.f90
    !   which calls process_cosmo_conv -> run_operator in mo_cosmo_conv.f90

    if (mec_sloti == mec_slotn) then
       call setup_veri_obs (obs, x)
       call check_obs      (obs% o)
       lqc = (mod(fg_check,4) == 1)
       if (interpolation <= 0) then
          !----------------------------
          ! use all time slices at once
          !----------------------------
          call apply_H  (H_det , obs, atm,                      bg=first)
       else
          !------------------------
          ! use one time slice only
          !------------------------
          call apply_H  (H_det , obs, atm(mec_slotn:mec_slotn), bg=first)
       end if

       call thin_superob_tovs(obs, H_det)
       call process_obs (TSK_R,    obs)


       if (first) then
          !-------------------------
          ! apply consistency checks
          !-------------------------
          call check_cons (obs% o)
          !---------------------------------
          ! apply observation specific rules
          !---------------------------------
          call check_rule (obs% o) ! check for specific rules
          if (n_dace_op > 0) then
             !call read_blacklists
             call check_black (obs% o) ! check for blacklisting
          end if
          if (fg_check >= 4) then
             !----------------------------------------------------
             ! apply quality checks as done in LETKF
             ! but no fg-check (no reasonable e_o, e_fg available)
             !----------------------------------------------------
             call check_black (obs% o)        ! check for blacklisting
             call check_gross (obs% o, H_det) ! generic gross check
          endif
          call check_suff (obs% o) ! check for sufficient data in report
          first = .false.
       end if
       

    end if

    call dace_close_output ()

#endif /* __DACE__ */

  end subroutine run_dace_op
  !============================================================================
  subroutine finish_dace ()
    !-------------------
    ! 1) write fof-Files
    ! 2) clean up DACE
    !-------------------
#ifndef __DACE__
    CALL finish ("finish_dace","DACE coupling requested but not compiled in")
#else
    character(22) :: comment
    type(t_vector) ::  H_btcs              !clear-sky Tbs
    integer, allocatable :: op_flag(:)  ! operator flag for each fdbk file category

    call dace_open_output ()

    if (ens_mem > 1) then
       write (comment,'("ensemble member",i4)') ens_mem
    else
       comment =        "deterministic run"
    end if
    !------------------------------------------
    ! create initial feedback files if required
    !------------------------------------------
    call fix_gp_indices ()
    call write_fdbk_3dv (mon, obs, atm(mec_slot1)% time, grid, step=4 ,  &
                              name='', prefix=prefix_out, comment=comment)
    call destruct       (mon)

    !------------------------------
    ! write report usage statistics
    !------------------------------
    call derive_rept_stat (obs% o)
    call gather_rept_stat
    call print_rept_stat  (unit=dace_output_unit)

    !------------------------
    ! append to feedback file
    !------------------------
    call message ("","")
    call message ("icon2dace","writing feedback files")
    call add_veri (H_det, obs, atm(mec_slotn), ensm=ens_mem)
    !--------------------------------------------------
    ! clear-sky model equivalent brightness temperatures
    !--------------------------------------------------
    call construct (H_btcs, obs% oi,'H_btcs' )  !clear-sky TB
    allocate(op_flag(n_source))
    call prep_H_btcs(H_btcs, obs,op_flag)
    if (any(op_flag == OF_BT_CLEAR_SKY)) then
      call add_veri (H_btcs,obs, atm(mec_slotn), op_flag=op_flag)
    end if
    call destruct(H_btcs)
    !-------------------------------
    ! clean up observation operators
    !-------------------------------
    call destruct (H_det)

    !---------
    ! clean up
    !---------
    call destruct (obs)

    if (associated (atm)) then
       call destruct (atm)
       deallocate    (atm)
       nullify       (atm)
    end if

    if (associated (grid)) then
       call destruct (grid)
       nullify (grid)
    end if

    call message ("icon2dace","Done.")
    call message ("","")
    call dace_close_output ()

#endif /* __DACE__ */

  end subroutine finish_dace
  !============================================================================
#ifdef __DACE__

  subroutine fix_gp_indices ()
    !----------------------------------------------------------------------
    ! Convert gridpoint indices (i,j,d) from icon2dace internal conventions
    ! to unique external representation
    !----------------------------------------------------------------------
    integer :: nb                       ! number of observation 'boxes'
    integer :: ib, is                   ! loop indices
    integer :: ijd(3)                   ! gridpoint (i,j,d)
    integer :: n                        ! counter
    type(t_obs) ,pointer :: o           ! pointer to current observation box
    type(t_spot),pointer :: spt         ! pointer to report

    n  = 0
    nb = size (obs% o)                          ! number of 'boxes'
    do ib = 1, nb                               ! loop over 'boxes'
       if (obs% o(ib)% pe /= dace% pe) cycle    ! handled on this PE only
       if (obs% o(ib)% n_spot == 0)    cycle
       o => obs%  o(ib)
       do is = 1, o% n_spot                     ! loop over reports
          spt => o% spot(is)
          ijd(1:3) = spt% col% h% ijdp(1:3)
          if (all (ijd > 0)) then
             spt% col% h% ijdp(1) = ext_glb_idx(ijd(1))
             spt% col% h% ijdp(2) = 1
             n = n + 1
          end if
       end do
    end do

    if (dbg_level > 1) then
       n = p_sum (n)
       if (dace% lpio) write (dace_output_unit,*) "fix_gp_indices: ", n, "reports"
    end if
  end subroutine fix_gp_indices

  !============================================================================

  subroutine read_nml_mec_obs ()
    !------------------------
    ! read namelist /MEC_OBS/
    !------------------------
    integer              :: ierr
    !-------------
    ! set defaults
    !-------------
    obstypes     = 'TEMP PILOT SYNOP DRIBU AIREP SATOB' ! obstypes to process
    prefix_out   = 'fof' ! feedback output file prefix
    interpolation= -1    ! >0: use nth slot, 0:nn -1:interval
    fg_check     = -1    ! switch how to apply quality checks

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('MEC_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=MEC_OBS, iostat=ierr)
        if (ierr/=0) call finish ('read_nml_mec_obs','ERROR in namelist /MEC_OBS/')
#else
        read (nnml ,nml=MEC_OBS)
#endif
      end select
      !---------
      ! printout
      !---------
      flush (6)
      write (6,'(a,/ )')   repeat('-',79)
      write (6,'(a,/ )')   ' namelist /MEC_OBS/'
      write (6,'(a,a )')   '  obstypes      = ',trim(obstypes)
      write (6,'(a,a )')   '  prefix_out    = ',     prefix_out
      write (6,'(a,i0)')   '  fg_check      = ',fg_check
      select case (interpolation)
      case default
        write (6,'(a,i0)') '  interpolation = ',interpolation
      case (-1)
        write (6,'(a)')    '  interpolation = -1 : using time interpolation'
      case (0)
        write (6,'(a)')    '  interpolation =  0 : using nearest time index'
      case (1)
        write (6,'(a)')    '  interpolation = +1 : NO time interpolation'
      end select
      write (6,'()')
      flush (6)
    end if
    !---------------------------
    ! broadcast namelist entries
    !---------------------------
    call p_bcast (obstypes      ,dace% pio)
    call p_bcast (prefix_out    ,dace% pio)
    call p_bcast (interpolation ,dace% pio)
    call p_bcast (fg_check      ,dace% pio)
    if (interpolation < -1) call finish ("read_nml_mec_obs",    &
                                         "invalid interpolation")
    if (interpolation >  1) call finish ("read_nml_mec_obs",              &
                                         "interpolation > 1 not supported")

  end subroutine read_nml_mec_obs

  !============================================================================

  subroutine dace_open_output (first)
    logical, intent(in), optional :: first  ! Check if writable on first open
    !-------------------------------------------------------------
    ! (Re-)Open DACE output_unit on I/O processor (in append mode)
    !-------------------------------------------------------------
    integer :: iu, ios

    if (.not. dace% lpio .or. dace_output_unit == 0 &
                         .or. dace_output_file == "") return
    if (present (first)) then
       if (first) then
          iu = get_unit_number ()
          open (iu, file=dace_output_file, status='new', iostat=ios)
          if (ios /= 0) call finish ("dace_open_output",                      &
                                     "cannot write: "//trim (dace_output_file))
          close (iu, status='delete')
          call return_unit_number (iu)
       end if
    end if
    open (unit=dace_output_unit, file=dace_output_file, position="APPEND")
  end subroutine dace_open_output

  !============================================================================

  subroutine dace_close_output ()
    !----------------------------------------
    ! Close DACE output_unit on I/O processor
    !----------------------------------------
    if (.not. dace% lpio .or. dace_output_unit == 0 &
                         .or. dace_output_file == "") return
    close (unit=dace_output_unit)
  end subroutine dace_close_output

#endif /* __DACE__ */
  !============================================================================

end module mo_icon2dace
