! This module contains the I/O routines for lateral boundary nudging
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
!
!
! ------------------------------------------------------------------------
! Which fields are read from the lateral boundary conditions file?
! ------------------------------------------------------------------------
!
! This question is answered independently from the "init_icon"
! namelist parameter of the initial state setup!
! - If "VN" is available, then it is read from file, otherwise "U","V".
! - "W" is optional and read if available in the input data (note that "W" may in fact contain OMEGA).
! - "QV", "QC", "QI" are always read
! - "QR", "QS" are read if available
!
! The other fields for the lateral boundary conditions are read
! from file, based on the following decision tree.  Note that the
! basic distinction between input from non-hydrostatic and
! hydrostatic models is made on the availability of the HHL field.
!
!                            +------------------+
!                            |  HHL available?  |
!                            +------------------+
!                                     |
!                     ______yes_______|________no________
!                     |                                  |
!         +--------------------------+            +------------------------+
!         | RHO & THETA_V available? |            | PS,GEOP & T available? |
!         +--------------------------+            +------------------------+
!                     |                                        |
!           ____yes___|____no______                    ___yes__|________no___________
!           |                      |                  |                              |
!           |               +----------------+        |                          +--------+
! * read HHL,RHO,THETA_V    | P,T available? |     * read in PS,GEOP,T           | ERROR! |
! * read W if available     |                |     * read OMEGA if available     |        |
! * ignore PS,GEOP          +----------------+     * compute P,HHL               +--------+
! * diagnose P,T                   |               * CALL OMEGA -> W
!                          ___yes__|___no____
!                         |                  |
!                         |               +--------+
!                    * read HHL,P,T       | ERROR! |
!                    * ignore PS,GEOP     +--------+
!                    * read W if available
!
! Afterwards, we
! - re-compute the virtual temperature (inside the vertical
!   interpolation subroutine)
! - (re-)compute RHO (inside the vertical interpolation subroutine)
! - perform vertical interpolation
!
! ------------------------------------------------------------------------
! Read-in of lateral boundary data: General Overview of the Implementation
! ------------------------------------------------------------------------
!
! Note: This short documentation focuses on the "asynchronous
! prefetching" mode only, the old (possibly deprecated) synchronous
! read-in of the boundary data, "mo_sync_latbc.f90", is not covered.
!
! Read-in of boundary data is invoked via
!   CALL recv_latbc_data
! in the time loop (module "mo_nh_stepping").
!
! Modules related to (asynchronous) boundary data read-in:
! --------------------------------------------------------
!
! src/io/atmo/mo_async_latbc.f90              : Setup of the boundary data read-in functionality
! src/io/atmo/mo_async_utils.f90              : * Initialisation, allocation
!                                               * Top level routines: "read-in" (e.g. "recv_latbc_data")
!                                               * Top level routines: "fetch" (see explanation below)
! src/io/atmo/mo_async_latbc_types.f90        : Declaration of data types.
! src/io/atmo/mo_latbc_read_recv.f90          : Low level routines:
!                                               read-in and sending of field data via MPI
!
! src/atm_dyn_iconam/mo_initicon_types.f90    : Type declarations for the final destination buffers
! src/atm_dyn_iconam/mo_nh_nest_utilities.f90 : Actual usage of the boundary data: buffers -> tendencies
!
! src/namelists/mo_limarea_nml.f90            : Namelist "limarea_nml"
! src/configure_model/mo_limarea_config.f90   : Configuration state (filled by namelist)
!
! Important notes for understanding the read-in process:
! ------------------------------------------------------
!
! * General switch: "latbc_config%itype_latbc" (INTEGER)
!
!   This is set to the value LATBC_TYPE_EXT when the field data is
!   read from an external file (default situation).
!
! * General switch: "latbc_config%lsparse_latbc" (LOGICAL)
!
!   Lateral boundary data can be provided either as a boundary strip
!   or as the complete field, including the (unused) interior points.
!
! * Variables which are considered for boundary data read-in are
!   contained in the variable group LATBC_PREFETCH_VARS.
!   Buffers are set up and allocated for group members.
!
! * Most important data types: "t_latbc_data", and contained therein: "t_buffer"
!
!   Defined in src/io/atmo/mo_async_latbc_types.f90 These derived
!   data types contain the raw data, that has been read from file,
!   together with a number of configurations settings
!   (LOGICAL-Flags), e.g. if specific optional fields are requested.
!
! * "Decision tree":
!
!   During the setup phase, the boundary data module opens the first
!   boundary data file and analyzes its contents (subroutine
!   "check_variables" in "mo_async_latbc", contains numerous calls of
!   "test_cdi_varID").
!
!   According to the available data, several LOGICAL flags
!   "lread_XXX" are set in the intermediate buffer "latbc%buffer".
!
! * Distinction between "read-in" and "fetch"
!
!   (In the module "mo_async_latbc_utils":)
!   "read-in": A list of variables is read from file into a buffer.
!              * executed on the read-in process.
!              * reads all fields of the group LATBC_TYPE_EXT
!              * the intermediate buffer may have only the size of
!              * the boundary strip.
!              * most important subroutine in this context:
!                "prefetch_latbc_data" and, therein, "CALL
!                prefetch_cdi_3d".
!
!   "Fetch": Copies the data from the intermediate buffer into
!            ICON-internal variables.
!              * executed on the compute processes.
!              * copy/processing of the fields depends on the
!                "lread_XXX" flags (see above), followed by vertical
!                interpolation.
!              * most important subroutine in this context:
!                compute_latbc_intp_data (mo_async_utils) and,
!                therein, "CALL fetch_from_buffer".
!              * Note that the "ICON-internal variables" are not the
!                prognostic fields, but again intermediate buffers of
!                the data type "t_init_state" (Modul
!                "mo_initicon_types").
!                Example: "latbc%latbc_data(tlev)%atm_in%u", which is
!                then later processed into tendencies (in
!                "mo_nh_nest_utilities").
!
! * Recipe: How to implement additional boundary data
!
!   1) Modify the "add_var" calls of the new fields in the ICON code,
!      such that these variables become members of the group
!      "LATBC_TYPE_EXT".
!   2) Extend the data type "t_buffer" by a LOGICAL flag for the new
!      field "lread_XXX".
!   3) Subroutine "check_variables" in "mo_async_latbc": Additional
!      test, if the new field is available in the boundary data file;
!      set the flag "lread_XXX" accordingly.
!   4) Extend the "t_latbc_data" data structure by a buffer for the
!      new field, similar to the contents of
!      "latbc%latbc_data(tlev)%atm_in".
!   5) Subroutine "compute_latbc_icon_data": Additional call to
!      "fetch_from_buffer", contained in an IF-condition "IF
!      (lread_XXX) ...".
!   6) Implement the vertical interpolation for the new field.
!   7) Use the new data in the application code.
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

MODULE mo_async_latbc

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_f_pointer, c_int
  USE mo_kind,                      ONLY: sp
  USE mo_exception,                 ONLY: finish, message, message_text
  USE mo_mpi,                       ONLY: stop_mpi, my_process_is_pref, &
       &                                  my_process_is_mpi_test, p_real_sp, &
       &                                  p_reduce, mpi_sum, p_allgather, &
       &                                  p_allgatherv, p_bcast
  USE mo_parallel_config,           ONLY: nproma, num_prefetch_proc
  USE mo_decomposition_tools,       ONLY: t_grid_domain_decomp_info
  USE mo_model_domain,              ONLY: p_patch, t_patch
  USE mo_mpi, ONLY: p_comm_work, p_comm_work_pref, p_comm_work_2_pref, my_process_is_work, &
    & num_work_procs, p_pe_work, p_work_pe0, p_comm_work_pref_compute_pe0
  USE mo_time_config,               ONLY: time_config
  USE mo_reorder_info,              ONLY: t_reorder_info
  USE mo_async_latbc_types,         ONLY: t_patch_data, t_buffer, &
                                          t_latbc_data, t_mem_win
  USE mo_async_latbc_utils,         ONLY: prefetch_latbc_data, async_init_latbc_data,&
       &                                  compute_wait_for_async_pref, compute_shutdown_async_pref, &
       &                                  async_pref_send_handshake,  async_pref_wait_for_start, &
       &                                  allocate_pref_latbc_data
  USE mo_impl_constants,            ONLY: SUCCESS, TIMELEVEL_SUFFIX, vname_len
  USE mo_cdi_constants,             ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
  USE mo_communication,             ONLY: idx_no, blk_no
  USE mo_nonhydro_state,            ONLY: p_nh_state, p_nh_state_lists
  USE mo_intp_data_strc,            ONLY: p_int_state
  USE mo_ext_data_state,            ONLY: ext_data
  USE mo_var_metadata_types,        ONLY: t_var_metadata_ptr
  USE mo_var_list_register_utils,   ONLY: vlr_group, vlr_replicate
  USE mo_var_list_register,         ONLY: t_vl_register_iter
  USE mo_var_metadata,              ONLY: get_var_name
  USE mo_var,                       ONLY: t_var
  USE mo_limarea_config,            ONLY: latbc_config
  USE mo_dictionary,                ONLY: t_dictionary
  USE mo_util_string,               ONLY: add_to_list, tolower
  USE mo_run_config,                ONLY: iqs
  USE mo_util_sort,                 ONLY: quicksort
  USE mo_time_config,               ONLY: time_config
  USE mtime,                        ONLY: datetime, OPERATOR(+)
  USE mo_cdi,                       ONLY: vlistInqVarZaxis, streamInqVlist, &
       &                                  vlistNvars, zaxisInqSize, vlistInqVarName,         &
       &                                  streamInqFiletype,                                 &
       &                                  FILETYPE_NC2, FILETYPE_NC4, FILETYPE_NC,           &
       &                                  cdi_max_name
  USE mo_io_util,                   ONLY: read_netcdf_int_1d, t_netcdf_att_int
  USE mo_util_cdi,                  ONLY: test_cdi_varID

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'



    PUBLIC :: init_prefetch
    PUBLIC :: close_prefetch

  !------------------------------------------------------------------------------------------------
  ! CONSTANTS
  !------------------------------------------------------------------------------------------------

  ! common constant strings
  CHARACTER(*), PARAMETER :: modname = 'mo_async_latbc'
  ! variables in this group are ICON, COSMO or IFS data which are
  ! read by the prefetch PE:
  CHARACTER(*), PARAMETER :: LATBC_PREFETCH_VARS = 'LATBC_PREFETCH_VARS'
  INTEGER,          PARAMETER :: MAX_NUM_GRPVARS = 200

  INTERFACE
  FUNCTION streaminqfileid(streamid) BIND(c, name='streamInqFileID') &
    RESULT(fileid)
    IMPORT :: c_int
    INTEGER(c_int), VALUE :: streamid
    INTEGER(c_int) :: fileid
  END FUNCTION streaminqfileid
  END INTERFACE

CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Close all name_list files and deallocate variables
  SUBROUTINE close_prefetch()







    END SUBROUTINE close_prefetch
  !--------------------------------------------------------------------------
  !
  !> Initialize data structures for prefetching boundary data
  !
  !  This routine is called after reading the namelists AND setting up
  !  the domains and variables.
  !

  ! ------------------------------------------------------------------------
  ! replicate data on prefetch proc
  ! ------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  !> Replicates data (mainly the variable lists) needed for async prefetching
  !  on the prefetching procs.
  SUBROUTINE init_prefetch(latbc)
    TYPE(t_latbc_data), INTENT(INOUT), TARGET :: latbc
  END SUBROUTINE init_prefetch


  !-------------------------------------------------------------------------------------------------
  !> open files containing first variable list and analysis
  !
  SUBROUTINE read_init_file(latbc, StrLowCasegrp, latbc_varnames_dict, p_patch)
    TYPE (t_latbc_data),        INTENT(INOUT) :: latbc
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: StrLowCasegrp(:) !< grp name in lower case letter
    TYPE (t_dictionary),        INTENT(IN)    :: latbc_varnames_dict
    TYPE(t_patch), OPTIONAL,    INTENT(IN)    :: p_patch(:)

  END SUBROUTINE read_init_file



  !-------------------------------------------------------------------------------------------------
  !> Replicates data needed for async prefetch on the prefetch proc.
  !  ATTENTION: The data is not completely replicated, only as far as needed for prefetch.
  !
  !  This routine has to be called by all PEs (work and prefetch)
  !

END MODULE mo_async_latbc
