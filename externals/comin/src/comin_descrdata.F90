!> @file comin_descrdata.F90
!! @brief Accessor functions for ComIn descriptive data structures.
!
!  @authors 10/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_descrdata

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_ptr, c_loc, c_f_pointer, c_null_ptr, c_double, c_bool
  USE comin_setup_constants, ONLY: DOMAIN_UNDEFINED, wp
  USE comin_state,           ONLY: state
  USE comin_descrdata_types, ONLY: t_comin_descrdata_global,       &
    &                              t_comin_descrdata_domain,       &
    &                              t_comin_descrdata_simulation_interval,   &
    &                              t_comin_descrdata_domain_cells, &
    &                              t_comin_descrdata_domain_edges, &
    &                              t_comin_descrdata_domain_verts
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS
  USE comin_errhandler,      ONLY: comin_plugin_finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_descrdata_set_global, comin_descrdata_get_global
  PUBLIC :: comin_descrdata_get_domain, comin_descrdata_set_domain
  PUBLIC :: comin_descrdata_set_simulation_interval, comin_descrdata_get_simulation_interval
  PUBLIC :: comin_descrdata_finalize
  PUBLIC :: comin_descrdata_get_index, comin_descrdata_get_block, comin_descrdata_get_cell_indices
  PUBLIC :: comin_descrdata_get_cell_npromz, comin_descrdata_get_edge_npromz, comin_descrdata_get_vert_npromz
  PUBLIC :: comin_descrdata_index_lookup_glb2loc_cell
  PUBLIC :: comin_current_get_datetime
  PUBLIC :: comin_current_set_datetime
  PUBLIC :: comin_descrdata_set_timesteplength, comin_descrdata_get_timesteplength



CONTAINS

  !> Fill global data.
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_set_global(comin_global_info)
    TYPE(t_comin_descrdata_global), INTENT(IN)  :: comin_global_info
    state%comin_descrdata_global = comin_global_info
  END SUBROUTINE comin_descrdata_set_global

  !> Set up data type for grid data.
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_set_domain(comin_domain_info)
    TYPE(t_comin_descrdata_domain), INTENT(IN)  :: comin_domain_info(:)
    state%comin_descrdata_domain = comin_domain_info
  END SUBROUTINE comin_descrdata_set_domain

  !> Fill time stamp info.
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_set_simulation_interval(comin_time_info)
    TYPE(t_comin_descrdata_simulation_interval), INTENT(IN)  :: comin_time_info
    state%comin_descrdata_simulation_interval = comin_time_info
  END SUBROUTINE comin_descrdata_set_simulation_interval

  !> Retrieve time stamp info, current time information
  !! @ingroup plugin_interface
  SUBROUTINE comin_current_get_datetime(sim_time_current)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: sim_time_current
    sim_time_current = state%current_datetime
  END SUBROUTINE comin_current_get_datetime

  !> request simulation interval information, C interface
  SUBROUTINE comin_current_get_datetime_c(val, len) &
    & BIND(C, NAME="comin_current_get_datetime")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%current_datetime)
    len = LEN_TRIM(state%current_datetime)
  END SUBROUTINE comin_current_get_datetime_c

  !> Update time stamp info, current time information.
  !! @ingroup host_interface
  SUBROUTINE comin_current_set_datetime(sim_time_current)
    CHARACTER(LEN=*),  INTENT(IN) :: sim_time_current
    state%current_datetime = sim_time_current
  END SUBROUTINE comin_current_set_datetime

  !> request a pointer to the global data type
  !! @ingroup plugin_interface
  FUNCTION comin_descrdata_get_global()
    TYPE(t_comin_descrdata_global), POINTER :: comin_descrdata_get_global
    ! local

    comin_descrdata_get_global => NULL()
    comin_descrdata_get_global => state%comin_descrdata_global
    IF(.NOT. ASSOCIATED(comin_descrdata_get_global)) THEN
      CALL comin_plugin_finish("comin_descrdata_get_global", " ERROR: Pointer not associated.")
    END IF
  END FUNCTION comin_descrdata_get_global

  !> request a pointer to the grid data type for a specific computational domain
  !! @ingroup plugin_interface
  FUNCTION comin_descrdata_get_domain(jg)
    INTEGER,                       INTENT(IN)  :: jg
    ! local
    TYPE(t_comin_descrdata_domain), POINTER    :: comin_descrdata_get_domain

    comin_descrdata_get_domain => NULL()
    comin_descrdata_get_domain => state%comin_descrdata_domain(jg)
    IF(.NOT. ASSOCIATED(comin_descrdata_get_domain)) THEN
      CALL comin_plugin_finish("comin_descrdata_get_domain", " ERROR: Pointer not associated.")
    END IF

  END FUNCTION comin_descrdata_get_domain

  !> request a pointer to simulation status
  !! @ingroup plugin_interface
  FUNCTION comin_descrdata_get_simulation_interval()
    TYPE(t_comin_descrdata_simulation_interval), POINTER :: comin_descrdata_get_simulation_interval

    comin_descrdata_get_simulation_interval => state%comin_descrdata_simulation_interval
    IF(.NOT. ASSOCIATED(comin_descrdata_get_simulation_interval)) THEN
      CALL comin_plugin_finish("comin_descrdata_get_simulation_interval", " ERROR: Pointer not associated.")
    END IF
  END FUNCTION comin_descrdata_get_simulation_interval

  !> request simulation interval information, C interface
  SUBROUTINE comin_descrdata_get_simulation_interval_exp_start(val, len) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_exp_start")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%comin_descrdata_simulation_interval%exp_start)
    len = LEN_TRIM(state%comin_descrdata_simulation_interval%exp_start)
  END SUBROUTINE comin_descrdata_get_simulation_interval_exp_start

  !> request simulation interval information, C interface
  SUBROUTINE comin_descrdata_get_simulation_interval_exp_stop(val, len) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_exp_stop")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%comin_descrdata_simulation_interval%exp_stop)
    len = LEN_TRIM(state%comin_descrdata_simulation_interval%exp_stop)
  END SUBROUTINE comin_descrdata_get_simulation_interval_exp_stop

  !> request simulation interval information, C interface
  SUBROUTINE comin_descrdata_get_simulation_interval_run_start(val, len) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_run_start")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%comin_descrdata_simulation_interval%run_start)
    len = LEN_TRIM(state%comin_descrdata_simulation_interval%run_start)
  END SUBROUTINE comin_descrdata_get_simulation_interval_run_start

  !> request simulation interval information, C interface
  SUBROUTINE comin_descrdata_get_simulation_interval_run_stop(val, len) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_run_stop")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%comin_descrdata_simulation_interval%run_stop)
    len = LEN_TRIM(state%comin_descrdata_simulation_interval%run_stop)
  END SUBROUTINE comin_descrdata_get_simulation_interval_run_stop

  !> Receive pointer on array storing timestep information for all domains
  !! @ingroup plugin_interface
  FUNCTION comin_descrdata_get_timesteplength(jg) BIND(C)
    REAL(wp) :: comin_descrdata_get_timesteplength
    INTEGER(c_int), INTENT(IN), VALUE :: jg

    comin_descrdata_get_timesteplength = state%comin_descrdata_timesteplength(jg)
  END FUNCTION comin_descrdata_get_timesteplength

  !> Fill array with timestep.
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_set_timesteplength(jg, dt_current)
    INTEGER,               INTENT(IN)  :: jg
    REAL(wp),              INTENT(IN)  :: dt_current

    IF (.NOT. ALLOCATED(state%comin_descrdata_timesteplength)) THEN
      ALLOCATE(state%comin_descrdata_timesteplength(state%comin_descrdata_global%n_dom+4))
    END IF
    state%comin_descrdata_timesteplength(jg) = dt_current
  END SUBROUTINE comin_descrdata_set_timesteplength

  !> Clean descriptive data structure in ComIn
  !> currently no content but keep for future use
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_finalize()

  END SUBROUTINE comin_descrdata_finalize

  !!
  !> auxiliary functions taken from ICON, version 2.6.5
  !!

  !> from mo_parallel_config
  !> names of routines in ICON are: blk_no, idx_no
  !-------------------------------------------------------------------------
  ! The following two functions are for conversion of 1D to 2D indices and vice versa
  !
  ! Treatment of 0 (important for empty domains) and negative numbers:
  !
  ! Converting 1D => 2D:
  !
  ! 0 always is mapped to blk_no = 1, idx_no = 0
  ! negative numbers: Convert usings ABS(j) and negate idx_no
  !
  ! Thus: blk_no >= 1 always!
  !       idx_no > 0  for j > 0
  !       idx_no = 0  for j = 0
  !       idx_no < 0  for j < 0
  !
  ! This mimics mostly the behaviour of reshape_idx in mo_model_domimp_patches
  ! with a difference for nproma=1 and j=0 (where reshape_idx returns blk_no=0, idx_no=1)
  !
  ! The consistent treatment of 0 in the above way is very important for empty domains
  ! where start_index=1, end_index=0
  !
  ! Converting 2D => 1D:
  ! Trying to invert the above and catching cases with blk_no < 1
  !-------------------------------------------------------------------------

  !> Auxiliary function: conversion of 1D to 2D indices.
  !! @ingroup plugin_interface
  INTEGER(c_int) FUNCTION comin_descrdata_get_block(idx1D) BIND(C, name="comin_descrdata_get_block")
    INTEGER(c_int), INTENT(IN), VALUE :: idx1D
    comin_descrdata_get_block = MAX((ABS(idx1D)-1)/state%comin_descrdata_global%nproma + 1, 1) ! i.e. also 1 for idx1D=0, nproma=1
  END FUNCTION comin_descrdata_get_block

  !> Auxiliary function: conversion of 1D to 2D indices.
  !! @ingroup plugin_interface
  INTEGER(c_int) FUNCTION comin_descrdata_get_index(idx1D) BIND(C, name="comin_descrdata_get_index")
    INTEGER(c_int), INTENT(IN), VALUE :: idx1D
    IF(idx1D==0) THEN
      comin_descrdata_get_index = 0
    ELSE
      comin_descrdata_get_index = SIGN(MOD(ABS(idx1D)-1,state%comin_descrdata_global%nproma)+1, idx1D)
    ENDIF
  END FUNCTION comin_descrdata_get_index

  !> Computes the start and end indices of do loops for cell-based variables.
  !! @ingroup plugin_interface
  !!
  !! From ICON's `mo_loopindices`; name of corresponding ICON routine: `get_indices_c`.
  !!
  SUBROUTINE comin_descrdata_get_cell_indices(jg, i_blk, i_startblk, i_endblk, i_startidx, &
       i_endidx, irl_start, irl_end) &
    &  BIND(C, NAME="comin_descrdata_get_cell_indices")

    INTEGER(c_int), INTENT(IN), VALUE :: jg         ! Patch index for comin_domain
    INTEGER(c_int), INTENT(IN), VALUE :: i_blk      ! Current block (variable jb in do loops)
    INTEGER(c_int), INTENT(IN), VALUE :: i_startblk ! Start block of do loop
    INTEGER(c_int), INTENT(IN), VALUE :: i_endblk   ! End block of do loop
    INTEGER(c_int), INTENT(IN), VALUE :: irl_start  ! refin_ctrl level where do loop starts
    INTEGER(c_int), INTENT(IN), VALUE :: irl_end    ! refin_ctrl level where do loop ends

    INTEGER(c_int), INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (jc loop)

    IF (i_blk == i_startblk) THEN
      i_startidx = MAX(1,state%comin_descrdata_domain(jg)%cells%start_index(irl_start))
      i_endidx   = state%comin_descrdata_global%nproma
      IF (i_blk == i_endblk) i_endidx = state%comin_descrdata_domain(jg)%cells%end_index(irl_end)
    ELSE IF (i_blk == i_endblk) THEN
      i_startidx = 1
      i_endidx   = state%comin_descrdata_domain(jg)%cells%end_index(irl_end)
    ELSE
      i_startidx = 1
      i_endidx = state%comin_descrdata_global%nproma
    ENDIF

  END SUBROUTINE comin_descrdata_get_cell_indices

  !> Calculate `npromz` value for the blocking, needed for patch allocation.
  !> ... for the cells
  !! @ingroup plugin_interface
  !!
  !! NB: Avoid the case nblks=0 for empty patches, this might cause troubles
  !! if a empty patch is used somewhere (and npromz gets wrong in the formulas below).
  !!
  !! from `mo_setup_subdivision`; name of ICON routine: `npromz_c`.
  INTEGER(c_int) FUNCTION comin_descrdata_get_cell_npromz(jg) BIND(C)
    INTEGER(c_int),  INTENT(IN), VALUE  :: jg     ! domain index for comin_domain

    comin_descrdata_get_cell_npromz = state%comin_descrdata_domain(jg)%cells%ncells - &
             &  (state%comin_descrdata_domain(jg)%cells%nblks-1)*state%comin_descrdata_global%nproma
  END FUNCTION comin_descrdata_get_cell_npromz

  !> Calculate `npromz` value for the blocking, needed for patch allocation.
  !> ... for the edges
  !! @ingroup plugin_interface
  !!
  !! NB: Avoid the case nblks=0 for empty patches, this might cause troubles
  !! if a empty patch is used somewhere (and npromz gets wrong in the formulas below).
  !!
  !! from `mo_setup_subdivision`; name of ICON routine: `npromz_e`.
  INTEGER(c_int) FUNCTION comin_descrdata_get_edge_npromz(jg) BIND(C)
    INTEGER(c_int),  INTENT(IN), VALUE  :: jg     ! domain index for comin_domain

    comin_descrdata_get_edge_npromz = state%comin_descrdata_domain(jg)%edges%nedges - &
             &  (state%comin_descrdata_domain(jg)%edges%nblks-1)*state%comin_descrdata_global%nproma
  END FUNCTION comin_descrdata_get_edge_npromz

  !> Calculate `npromz` value for the blocking, needed for patch allocation.
  !> ... for the vertices
  !! @ingroup plugin_interface
  !!
  !! NB: Avoid the case nblks=0 for empty patches, this might cause troubles
  !! if a empty patch is used somewhere (and npromz gets wrong in the formulas below).
  !!
  !! from `mo_setup_subdivision`; name of ICON routine: `npromz_v`.
  INTEGER(c_int) FUNCTION comin_descrdata_get_vert_npromz(jg) BIND(C)
    INTEGER(c_int),  INTENT(IN), VALUE  :: jg     ! domain index for comin_domain

    comin_descrdata_get_vert_npromz = state%comin_descrdata_domain(jg)%verts%nverts - &
             &  (state%comin_descrdata_domain(jg)%verts%nblks-1)*state%comin_descrdata_global%nproma
  END FUNCTION comin_descrdata_get_vert_npromz

  !> Conversion of global cell index to MPI-process local index.
  !! @ingroup plugin_interface
  INTEGER(C_INT) FUNCTION comin_descrdata_index_lookup_glb2loc_cell(jg, global_idx) &
    & RESULT(loc) BIND(C)
    INTEGER(kind=C_INT), INTENT(IN), VALUE :: jg          !< domain index
    INTEGER(kind=C_INT), INTENT(IN), VALUE :: global_idx  !< global cell index
    loc = state%comin_descrdata_fct_glb2loc_cell(jg, INT(global_idx))
  END FUNCTION comin_descrdata_index_lookup_glb2loc_cell

  ! Query topo data routines generated by python script (comin_descrdata_get_domain.F90.py) in ../utils. !
!  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

! *** DO NOT EDIT MANUALLY!  Generated by python script in utils/. DO NOT EDIT MANUALLY! *** !

  SUBROUTINE comin_descrdata_get_domain_grid_filename(jg, grid_filename, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_grid_filename")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: grid_filename
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()
    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ASSOCIATED(p%grid_filename)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      &  "ERROR: Pointer of grid_filename not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%grid_filename)
    grid_filename = C_LOC(p%grid_filename)
  END SUBROUTINE comin_descrdata_get_domain_grid_filename

  SUBROUTINE comin_descrdata_get_domain_grid_uuid(jg, grid_uuid, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_grid_uuid")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: grid_uuid
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()
    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ASSOCIATED(p%grid_uuid)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      &  "ERROR: Pointer of grid_uuid not associated.")
    END IF
    arr_size = SHAPE(p%grid_uuid)
    grid_uuid = C_LOC(p%grid_uuid)
  END SUBROUTINE comin_descrdata_get_domain_grid_uuid

  SUBROUTINE comin_descrdata_get_domain_number_of_grid_used(jg, number_of_grid_used, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_number_of_grid_used")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: number_of_grid_used
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()
    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ALLOCATED(p%number_of_grid_used)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      &  "ERROR: Pointer of number_of_grid_used not associated.")
    END IF
    arr_size = SHAPE(p%number_of_grid_used)
    number_of_grid_used = C_LOC(p%number_of_grid_used)
  END SUBROUTINE comin_descrdata_get_domain_number_of_grid_used

  FUNCTION comin_descrdata_get_domain_id(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_id") &
      &  RESULT(id)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: id
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()

    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ASSOCIATED(p%id)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      & "ERROR: Pointer of id not associated.")
    END IF
    id  = p%id
  END FUNCTION comin_descrdata_get_domain_id

  FUNCTION comin_descrdata_get_domain_n_childdom(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_n_childdom") &
      &  RESULT(n_childdom)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: n_childdom
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()

    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ASSOCIATED(p%n_childdom)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      & "ERROR: Pointer of n_childdom not associated.")
    END IF
    n_childdom  = p%n_childdom
  END FUNCTION comin_descrdata_get_domain_n_childdom

  FUNCTION comin_descrdata_get_domain_dom_start(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_dom_start") &
      &  RESULT(dom_start)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    REAL(wp)                      :: dom_start
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()

    p => comin_descrdata_get_domain(jg)
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      & "ERROR: Pointer of dom_start not associated.")
    END IF
    dom_start  = p%dom_start
  END FUNCTION comin_descrdata_get_domain_dom_start

  FUNCTION comin_descrdata_get_domain_dom_end(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_dom_end") &
      &  RESULT(dom_end)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    REAL(wp)                      :: dom_end
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()

    p => comin_descrdata_get_domain(jg)
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      & "ERROR: Pointer of dom_end not associated.")
    END IF
    dom_end  = p%dom_end
  END FUNCTION comin_descrdata_get_domain_dom_end

  FUNCTION comin_descrdata_get_domain_nlev(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_nlev") &
      &  RESULT(nlev)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nlev
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()

    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ASSOCIATED(p%nlev)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      & "ERROR: Pointer of nlev not associated.")
    END IF
    nlev  = p%nlev
  END FUNCTION comin_descrdata_get_domain_nlev

  FUNCTION comin_descrdata_get_domain_nshift(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_nshift") &
      &  RESULT(nshift)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nshift
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()

    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ASSOCIATED(p%nshift)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      & "ERROR: Pointer of nshift not associated.")
    END IF
    nshift  = p%nshift
  END FUNCTION comin_descrdata_get_domain_nshift

  FUNCTION comin_descrdata_get_domain_nshift_total(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_nshift_total") &
      &  RESULT(nshift_total)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nshift_total
    !
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()

    p => comin_descrdata_get_domain(jg)
    IF (.NOT. ASSOCIATED(p%nshift_total)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain", &
      & "ERROR: Pointer of nshift_total not associated.")
    END IF
    nshift_total  = p%nshift_total
  END FUNCTION comin_descrdata_get_domain_nshift_total

  FUNCTION comin_descrdata_get_domain_cells(jg)
    USE iso_c_binding, ONLY: C_INT
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()
    TYPE(t_comin_descrdata_domain_cells), POINTER :: comin_descrdata_get_domain_cells

    p => comin_descrdata_get_domain(jg)
    comin_descrdata_get_domain_cells => p%cells
  END FUNCTION comin_descrdata_get_domain_cells

  FUNCTION comin_descrdata_get_domain_cells_ncells(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_ncells") &
      &  RESULT(ncells)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: ncells
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%ncells)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      & "ERROR: Pointer of ncells not associated.")
    END IF
    ncells  = p%ncells
  END FUNCTION comin_descrdata_get_domain_cells_ncells

  FUNCTION comin_descrdata_get_domain_cells_ncells_global(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_ncells_global") &
      &  RESULT(ncells_global)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: ncells_global
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%ncells_global)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      & "ERROR: Pointer of ncells_global not associated.")
    END IF
    ncells_global  = p%ncells_global
  END FUNCTION comin_descrdata_get_domain_cells_ncells_global

  FUNCTION comin_descrdata_get_domain_cells_nblks(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_nblks") &
      &  RESULT(nblks)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nblks
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%nblks)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      & "ERROR: Pointer of nblks not associated.")
    END IF
    nblks  = p%nblks
  END FUNCTION comin_descrdata_get_domain_cells_nblks

  FUNCTION comin_descrdata_get_domain_cells_max_connectivity(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_max_connectivity") &
      &  RESULT(max_connectivity)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: max_connectivity
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%max_connectivity)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      & "ERROR: Pointer of max_connectivity not associated.")
    END IF
    max_connectivity  = p%max_connectivity
  END FUNCTION comin_descrdata_get_domain_cells_max_connectivity

  SUBROUTINE comin_descrdata_get_domain_cells_num_edges(jg, num_edges, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_num_edges")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: num_edges
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%num_edges)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of num_edges not associated.")
    END IF
    arr_size = SHAPE(p%num_edges)
    num_edges = C_LOC(p%num_edges)
  END SUBROUTINE comin_descrdata_get_domain_cells_num_edges

  SUBROUTINE comin_descrdata_get_domain_cells_refin_ctrl(jg, refin_ctrl, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_refin_ctrl")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: refin_ctrl
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%refin_ctrl)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of refin_ctrl not associated.")
    END IF
    arr_size = SHAPE(p%refin_ctrl)
    refin_ctrl = C_LOC(p%refin_ctrl)
  END SUBROUTINE comin_descrdata_get_domain_cells_refin_ctrl

  SUBROUTINE comin_descrdata_get_domain_cells_start_index(jg, start_index, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_start_index")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: start_index
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%start_index)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of start_index not associated.")
    END IF
    arr_size = SHAPE(p%start_index)
    start_index = C_LOC(p%start_index)
  END SUBROUTINE comin_descrdata_get_domain_cells_start_index

  SUBROUTINE comin_descrdata_get_domain_cells_end_index(jg, end_index, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_end_index")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: end_index
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%end_index)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of end_index not associated.")
    END IF
    arr_size = SHAPE(p%end_index)
    end_index = C_LOC(p%end_index)
  END SUBROUTINE comin_descrdata_get_domain_cells_end_index

  SUBROUTINE comin_descrdata_get_domain_cells_start_block(jg, start_block, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_start_block")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: start_block
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%start_block)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of start_block not associated.")
    END IF
    arr_size = SHAPE(p%start_block)
    start_block = C_LOC(p%start_block)
  END SUBROUTINE comin_descrdata_get_domain_cells_start_block

  SUBROUTINE comin_descrdata_get_domain_cells_end_block(jg, end_block, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_end_block")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: end_block
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%end_block)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of end_block not associated.")
    END IF
    arr_size = SHAPE(p%end_block)
    end_block = C_LOC(p%end_block)
  END SUBROUTINE comin_descrdata_get_domain_cells_end_block

  SUBROUTINE comin_descrdata_get_domain_cells_child_id(jg, child_id, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_child_id")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: child_id
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%child_id)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of child_id not associated.")
    END IF
    arr_size = SHAPE(p%child_id)
    child_id = C_LOC(p%child_id)
  END SUBROUTINE comin_descrdata_get_domain_cells_child_id

  SUBROUTINE comin_descrdata_get_domain_cells_child_idx(jg, child_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_child_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: child_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%child_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of child_idx not associated.")
    END IF
    arr_size = SHAPE(p%child_idx)
    child_idx = C_LOC(p%child_idx)
  END SUBROUTINE comin_descrdata_get_domain_cells_child_idx

  SUBROUTINE comin_descrdata_get_domain_cells_child_blk(jg, child_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_child_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: child_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%child_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of child_blk not associated.")
    END IF
    arr_size = SHAPE(p%child_blk)
    child_blk = C_LOC(p%child_blk)
  END SUBROUTINE comin_descrdata_get_domain_cells_child_blk

  SUBROUTINE comin_descrdata_get_domain_cells_parent_glb_idx(jg, parent_glb_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_parent_glb_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: parent_glb_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%parent_glb_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of parent_glb_idx not associated.")
    END IF
    arr_size = SHAPE(p%parent_glb_idx)
    parent_glb_idx = C_LOC(p%parent_glb_idx)
  END SUBROUTINE comin_descrdata_get_domain_cells_parent_glb_idx

  SUBROUTINE comin_descrdata_get_domain_cells_parent_glb_blk(jg, parent_glb_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_parent_glb_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: parent_glb_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%parent_glb_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of parent_glb_blk not associated.")
    END IF
    arr_size = SHAPE(p%parent_glb_blk)
    parent_glb_blk = C_LOC(p%parent_glb_blk)
  END SUBROUTINE comin_descrdata_get_domain_cells_parent_glb_blk

  SUBROUTINE comin_descrdata_get_domain_cells_vertex_idx(jg, vertex_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_vertex_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: vertex_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%vertex_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of vertex_idx not associated.")
    END IF
    arr_size = SHAPE(p%vertex_idx)
    vertex_idx = C_LOC(p%vertex_idx)
  END SUBROUTINE comin_descrdata_get_domain_cells_vertex_idx

  SUBROUTINE comin_descrdata_get_domain_cells_vertex_blk(jg, vertex_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_vertex_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: vertex_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%vertex_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of vertex_blk not associated.")
    END IF
    arr_size = SHAPE(p%vertex_blk)
    vertex_blk = C_LOC(p%vertex_blk)
  END SUBROUTINE comin_descrdata_get_domain_cells_vertex_blk

  SUBROUTINE comin_descrdata_get_domain_cells_neighbor_blk(jg, neighbor_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_neighbor_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: neighbor_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%neighbor_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of neighbor_blk not associated.")
    END IF
    arr_size = SHAPE(p%neighbor_blk)
    neighbor_blk = C_LOC(p%neighbor_blk)
  END SUBROUTINE comin_descrdata_get_domain_cells_neighbor_blk

  SUBROUTINE comin_descrdata_get_domain_cells_neighbor_idx(jg, neighbor_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_neighbor_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: neighbor_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%neighbor_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of neighbor_idx not associated.")
    END IF
    arr_size = SHAPE(p%neighbor_idx)
    neighbor_idx = C_LOC(p%neighbor_idx)
  END SUBROUTINE comin_descrdata_get_domain_cells_neighbor_idx

  SUBROUTINE comin_descrdata_get_domain_cells_edge_idx(jg, edge_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_edge_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: edge_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%edge_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of edge_idx not associated.")
    END IF
    arr_size = SHAPE(p%edge_idx)
    edge_idx = C_LOC(p%edge_idx)
  END SUBROUTINE comin_descrdata_get_domain_cells_edge_idx

  SUBROUTINE comin_descrdata_get_domain_cells_edge_blk(jg, edge_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_edge_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: edge_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%edge_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of edge_blk not associated.")
    END IF
    arr_size = SHAPE(p%edge_blk)
    edge_blk = C_LOC(p%edge_blk)
  END SUBROUTINE comin_descrdata_get_domain_cells_edge_blk

  SUBROUTINE comin_descrdata_get_domain_cells_clon(jg, clon, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_clon")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: clon
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ALLOCATED(p%clon)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of clon not associated.")
    END IF
    arr_size = SHAPE(p%clon)
    clon = C_LOC(p%clon)
  END SUBROUTINE comin_descrdata_get_domain_cells_clon

  SUBROUTINE comin_descrdata_get_domain_cells_clat(jg, clat, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_clat")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: clat
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ALLOCATED(p%clat)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of clat not associated.")
    END IF
    arr_size = SHAPE(p%clat)
    clat = C_LOC(p%clat)
  END SUBROUTINE comin_descrdata_get_domain_cells_clat

  SUBROUTINE comin_descrdata_get_domain_cells_area(jg, area, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_area")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: area
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%area)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of area not associated.")
    END IF
    arr_size = SHAPE(p%area)
    area = C_LOC(p%area)
  END SUBROUTINE comin_descrdata_get_domain_cells_area

  SUBROUTINE comin_descrdata_get_domain_cells_hhl(jg, hhl, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_hhl")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: hhl
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ALLOCATED(p%hhl)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of hhl not associated.")
    END IF
    arr_size = SHAPE(p%hhl)
    hhl = C_LOC(p%hhl)
  END SUBROUTINE comin_descrdata_get_domain_cells_hhl

  SUBROUTINE comin_descrdata_get_domain_cells_glb_index(jg, glb_index, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_glb_index")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: glb_index
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%glb_index)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of glb_index not associated.")
    END IF
    arr_size = SHAPE(p%glb_index)
    glb_index = C_LOC(p%glb_index)
  END SUBROUTINE comin_descrdata_get_domain_cells_glb_index

  SUBROUTINE comin_descrdata_get_domain_cells_decomp_domain(jg, decomp_domain, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_cells_decomp_domain")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: decomp_domain
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_cells), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_cells(jg)
    IF (.NOT. ASSOCIATED(p%decomp_domain)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_cells", &
      &  "ERROR: Pointer of decomp_domain not associated.")
    END IF
    arr_size = SHAPE(p%decomp_domain)
    decomp_domain = C_LOC(p%decomp_domain)
  END SUBROUTINE comin_descrdata_get_domain_cells_decomp_domain

  FUNCTION comin_descrdata_get_domain_verts(jg)
    USE iso_c_binding, ONLY: C_INT
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()
    TYPE(t_comin_descrdata_domain_verts), POINTER :: comin_descrdata_get_domain_verts

    p => comin_descrdata_get_domain(jg)
    comin_descrdata_get_domain_verts => p%verts
  END FUNCTION comin_descrdata_get_domain_verts

  FUNCTION comin_descrdata_get_domain_verts_nverts(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_nverts") &
      &  RESULT(nverts)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nverts
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%nverts)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      & "ERROR: Pointer of nverts not associated.")
    END IF
    nverts  = p%nverts
  END FUNCTION comin_descrdata_get_domain_verts_nverts

  FUNCTION comin_descrdata_get_domain_verts_nverts_global(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_nverts_global") &
      &  RESULT(nverts_global)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nverts_global
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%nverts_global)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      & "ERROR: Pointer of nverts_global not associated.")
    END IF
    nverts_global  = p%nverts_global
  END FUNCTION comin_descrdata_get_domain_verts_nverts_global

  FUNCTION comin_descrdata_get_domain_verts_nblks(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_nblks") &
      &  RESULT(nblks)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nblks
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%nblks)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      & "ERROR: Pointer of nblks not associated.")
    END IF
    nblks  = p%nblks
  END FUNCTION comin_descrdata_get_domain_verts_nblks

  SUBROUTINE comin_descrdata_get_domain_verts_num_edges(jg, num_edges, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_num_edges")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: num_edges
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%num_edges)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of num_edges not associated.")
    END IF
    arr_size = SHAPE(p%num_edges)
    num_edges = C_LOC(p%num_edges)
  END SUBROUTINE comin_descrdata_get_domain_verts_num_edges

  SUBROUTINE comin_descrdata_get_domain_verts_refin_ctrl(jg, refin_ctrl, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_refin_ctrl")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: refin_ctrl
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%refin_ctrl)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of refin_ctrl not associated.")
    END IF
    arr_size = SHAPE(p%refin_ctrl)
    refin_ctrl = C_LOC(p%refin_ctrl)
  END SUBROUTINE comin_descrdata_get_domain_verts_refin_ctrl

  SUBROUTINE comin_descrdata_get_domain_verts_start_index(jg, start_index, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_start_index")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: start_index
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%start_index)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of start_index not associated.")
    END IF
    arr_size = SHAPE(p%start_index)
    start_index = C_LOC(p%start_index)
  END SUBROUTINE comin_descrdata_get_domain_verts_start_index

  SUBROUTINE comin_descrdata_get_domain_verts_end_index(jg, end_index, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_end_index")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: end_index
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%end_index)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of end_index not associated.")
    END IF
    arr_size = SHAPE(p%end_index)
    end_index = C_LOC(p%end_index)
  END SUBROUTINE comin_descrdata_get_domain_verts_end_index

  SUBROUTINE comin_descrdata_get_domain_verts_start_block(jg, start_block, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_start_block")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: start_block
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%start_block)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of start_block not associated.")
    END IF
    arr_size = SHAPE(p%start_block)
    start_block = C_LOC(p%start_block)
  END SUBROUTINE comin_descrdata_get_domain_verts_start_block

  SUBROUTINE comin_descrdata_get_domain_verts_end_block(jg, end_block, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_end_block")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: end_block
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%end_block)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of end_block not associated.")
    END IF
    arr_size = SHAPE(p%end_block)
    end_block = C_LOC(p%end_block)
  END SUBROUTINE comin_descrdata_get_domain_verts_end_block

  SUBROUTINE comin_descrdata_get_domain_verts_neighbor_blk(jg, neighbor_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_neighbor_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: neighbor_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%neighbor_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of neighbor_blk not associated.")
    END IF
    arr_size = SHAPE(p%neighbor_blk)
    neighbor_blk = C_LOC(p%neighbor_blk)
  END SUBROUTINE comin_descrdata_get_domain_verts_neighbor_blk

  SUBROUTINE comin_descrdata_get_domain_verts_neighbor_idx(jg, neighbor_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_neighbor_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: neighbor_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%neighbor_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of neighbor_idx not associated.")
    END IF
    arr_size = SHAPE(p%neighbor_idx)
    neighbor_idx = C_LOC(p%neighbor_idx)
  END SUBROUTINE comin_descrdata_get_domain_verts_neighbor_idx

  SUBROUTINE comin_descrdata_get_domain_verts_cell_idx(jg, cell_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_cell_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: cell_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%cell_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of cell_idx not associated.")
    END IF
    arr_size = SHAPE(p%cell_idx)
    cell_idx = C_LOC(p%cell_idx)
  END SUBROUTINE comin_descrdata_get_domain_verts_cell_idx

  SUBROUTINE comin_descrdata_get_domain_verts_cell_blk(jg, cell_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_cell_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: cell_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%cell_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of cell_blk not associated.")
    END IF
    arr_size = SHAPE(p%cell_blk)
    cell_blk = C_LOC(p%cell_blk)
  END SUBROUTINE comin_descrdata_get_domain_verts_cell_blk

  SUBROUTINE comin_descrdata_get_domain_verts_edge_idx(jg, edge_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_edge_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: edge_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%edge_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of edge_idx not associated.")
    END IF
    arr_size = SHAPE(p%edge_idx)
    edge_idx = C_LOC(p%edge_idx)
  END SUBROUTINE comin_descrdata_get_domain_verts_edge_idx

  SUBROUTINE comin_descrdata_get_domain_verts_edge_blk(jg, edge_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_edge_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: edge_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ASSOCIATED(p%edge_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of edge_blk not associated.")
    END IF
    arr_size = SHAPE(p%edge_blk)
    edge_blk = C_LOC(p%edge_blk)
  END SUBROUTINE comin_descrdata_get_domain_verts_edge_blk

  SUBROUTINE comin_descrdata_get_domain_verts_vlon(jg, vlon, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_vlon")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: vlon
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ALLOCATED(p%vlon)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of vlon not associated.")
    END IF
    arr_size = SHAPE(p%vlon)
    vlon = C_LOC(p%vlon)
  END SUBROUTINE comin_descrdata_get_domain_verts_vlon

  SUBROUTINE comin_descrdata_get_domain_verts_vlat(jg, vlat, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_verts_vlat")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: vlat
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_verts), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_verts(jg)
    IF (.NOT. ALLOCATED(p%vlat)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_verts", &
      &  "ERROR: Pointer of vlat not associated.")
    END IF
    arr_size = SHAPE(p%vlat)
    vlat = C_LOC(p%vlat)
  END SUBROUTINE comin_descrdata_get_domain_verts_vlat

  FUNCTION comin_descrdata_get_domain_edges(jg)
    USE iso_c_binding, ONLY: C_INT
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(t_comin_descrdata_domain), POINTER :: p => NULL()
    TYPE(t_comin_descrdata_domain_edges), POINTER :: comin_descrdata_get_domain_edges

    p => comin_descrdata_get_domain(jg)
    comin_descrdata_get_domain_edges => p%edges
  END FUNCTION comin_descrdata_get_domain_edges

  FUNCTION comin_descrdata_get_domain_edges_nedges(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_nedges") &
      &  RESULT(nedges)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nedges
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%nedges)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      & "ERROR: Pointer of nedges not associated.")
    END IF
    nedges  = p%nedges
  END FUNCTION comin_descrdata_get_domain_edges_nedges

  FUNCTION comin_descrdata_get_domain_edges_nedges_global(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_nedges_global") &
      &  RESULT(nedges_global)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nedges_global
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%nedges_global)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      & "ERROR: Pointer of nedges_global not associated.")
    END IF
    nedges_global  = p%nedges_global
  END FUNCTION comin_descrdata_get_domain_edges_nedges_global

  FUNCTION comin_descrdata_get_domain_edges_nblks(jg) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_nblks") &
      &  RESULT(nblks)
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    INTEGER(C_INT)                      :: nblks
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()

    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%nblks)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      & "ERROR: Pointer of nblks not associated.")
    END IF
    nblks  = p%nblks
  END FUNCTION comin_descrdata_get_domain_edges_nblks

  SUBROUTINE comin_descrdata_get_domain_edges_refin_ctrl(jg, refin_ctrl, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_refin_ctrl")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: refin_ctrl
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%refin_ctrl)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of refin_ctrl not associated.")
    END IF
    arr_size = SHAPE(p%refin_ctrl)
    refin_ctrl = C_LOC(p%refin_ctrl)
  END SUBROUTINE comin_descrdata_get_domain_edges_refin_ctrl

  SUBROUTINE comin_descrdata_get_domain_edges_start_index(jg, start_index, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_start_index")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: start_index
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%start_index)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of start_index not associated.")
    END IF
    arr_size = SHAPE(p%start_index)
    start_index = C_LOC(p%start_index)
  END SUBROUTINE comin_descrdata_get_domain_edges_start_index

  SUBROUTINE comin_descrdata_get_domain_edges_end_index(jg, end_index, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_end_index")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: end_index
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%end_index)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of end_index not associated.")
    END IF
    arr_size = SHAPE(p%end_index)
    end_index = C_LOC(p%end_index)
  END SUBROUTINE comin_descrdata_get_domain_edges_end_index

  SUBROUTINE comin_descrdata_get_domain_edges_start_block(jg, start_block, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_start_block")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: start_block
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%start_block)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of start_block not associated.")
    END IF
    arr_size = SHAPE(p%start_block)
    start_block = C_LOC(p%start_block)
  END SUBROUTINE comin_descrdata_get_domain_edges_start_block

  SUBROUTINE comin_descrdata_get_domain_edges_end_block(jg, end_block, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_end_block")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: end_block
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%end_block)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of end_block not associated.")
    END IF
    arr_size = SHAPE(p%end_block)
    end_block = C_LOC(p%end_block)
  END SUBROUTINE comin_descrdata_get_domain_edges_end_block

  SUBROUTINE comin_descrdata_get_domain_edges_child_id(jg, child_id, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_child_id")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: child_id
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%child_id)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of child_id not associated.")
    END IF
    arr_size = SHAPE(p%child_id)
    child_id = C_LOC(p%child_id)
  END SUBROUTINE comin_descrdata_get_domain_edges_child_id

  SUBROUTINE comin_descrdata_get_domain_edges_child_idx(jg, child_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_child_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: child_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%child_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of child_idx not associated.")
    END IF
    arr_size = SHAPE(p%child_idx)
    child_idx = C_LOC(p%child_idx)
  END SUBROUTINE comin_descrdata_get_domain_edges_child_idx

  SUBROUTINE comin_descrdata_get_domain_edges_child_blk(jg, child_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_child_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: child_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%child_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of child_blk not associated.")
    END IF
    arr_size = SHAPE(p%child_blk)
    child_blk = C_LOC(p%child_blk)
  END SUBROUTINE comin_descrdata_get_domain_edges_child_blk

  SUBROUTINE comin_descrdata_get_domain_edges_parent_glb_idx(jg, parent_glb_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_parent_glb_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: parent_glb_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%parent_glb_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of parent_glb_idx not associated.")
    END IF
    arr_size = SHAPE(p%parent_glb_idx)
    parent_glb_idx = C_LOC(p%parent_glb_idx)
  END SUBROUTINE comin_descrdata_get_domain_edges_parent_glb_idx

  SUBROUTINE comin_descrdata_get_domain_edges_parent_glb_blk(jg, parent_glb_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_parent_glb_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: parent_glb_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%parent_glb_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of parent_glb_blk not associated.")
    END IF
    arr_size = SHAPE(p%parent_glb_blk)
    parent_glb_blk = C_LOC(p%parent_glb_blk)
  END SUBROUTINE comin_descrdata_get_domain_edges_parent_glb_blk

  SUBROUTINE comin_descrdata_get_domain_edges_cell_idx(jg, cell_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_cell_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: cell_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%cell_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of cell_idx not associated.")
    END IF
    arr_size = SHAPE(p%cell_idx)
    cell_idx = C_LOC(p%cell_idx)
  END SUBROUTINE comin_descrdata_get_domain_edges_cell_idx

  SUBROUTINE comin_descrdata_get_domain_edges_cell_blk(jg, cell_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_cell_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: cell_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%cell_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of cell_blk not associated.")
    END IF
    arr_size = SHAPE(p%cell_blk)
    cell_blk = C_LOC(p%cell_blk)
  END SUBROUTINE comin_descrdata_get_domain_edges_cell_blk

  SUBROUTINE comin_descrdata_get_domain_edges_vertex_idx(jg, vertex_idx, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_vertex_idx")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: vertex_idx
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%vertex_idx)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of vertex_idx not associated.")
    END IF
    arr_size = SHAPE(p%vertex_idx)
    vertex_idx = C_LOC(p%vertex_idx)
  END SUBROUTINE comin_descrdata_get_domain_edges_vertex_idx

  SUBROUTINE comin_descrdata_get_domain_edges_vertex_blk(jg, vertex_blk, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_vertex_blk")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: vertex_blk
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(3)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ASSOCIATED(p%vertex_blk)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of vertex_blk not associated.")
    END IF
    arr_size = SHAPE(p%vertex_blk)
    vertex_blk = C_LOC(p%vertex_blk)
  END SUBROUTINE comin_descrdata_get_domain_edges_vertex_blk

  SUBROUTINE comin_descrdata_get_domain_edges_elon(jg, elon, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_elon")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: elon
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ALLOCATED(p%elon)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of elon not associated.")
    END IF
    arr_size = SHAPE(p%elon)
    elon = C_LOC(p%elon)
  END SUBROUTINE comin_descrdata_get_domain_edges_elon

  SUBROUTINE comin_descrdata_get_domain_edges_elat(jg, elat, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_domain_edges_elat")
    INTEGER(C_INT), INTENT(IN), VALUE  :: jg
    TYPE(C_PTR),    INTENT(OUT) :: elat
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(2)
    !
    TYPE(t_comin_descrdata_domain_edges), POINTER :: p => NULL()
    p => comin_descrdata_get_domain_edges(jg)
    IF (.NOT. ALLOCATED(p%elat)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_domain_edges", &
      &  "ERROR: Pointer of elat not associated.")
    END IF
    arr_size = SHAPE(p%elat)
    elat = C_LOC(p%elat)
  END SUBROUTINE comin_descrdata_get_domain_edges_elat

  ! Query global data routines generated by python script (comin_descrdata_get_global.F90.py) in ../utils. !
!  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

! *** DO NOT EDIT MANUALLY!  Generated by python script in utils/. DO NOT EDIT MANUALLY! *** !

  FUNCTION comin_descrdata_get_global_n_dom() &
      &  BIND(C, NAME="comin_descrdata_get_global_n_dom") &
      &  RESULT(n_dom)

    INTEGER(C_INT)                      :: n_dom
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of n_dom not associated.")
    END IF
    n_dom  = p%n_dom
  END FUNCTION comin_descrdata_get_global_n_dom

  FUNCTION comin_descrdata_get_global_max_dom() &
      &  BIND(C, NAME="comin_descrdata_get_global_max_dom") &
      &  RESULT(max_dom)

    INTEGER(C_INT)                      :: max_dom
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of max_dom not associated.")
    END IF
    max_dom  = p%max_dom
  END FUNCTION comin_descrdata_get_global_max_dom

  FUNCTION comin_descrdata_get_global_nproma() &
      &  BIND(C, NAME="comin_descrdata_get_global_nproma") &
      &  RESULT(nproma)

    INTEGER(C_INT)                      :: nproma
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of nproma not associated.")
    END IF
    nproma  = p%nproma
  END FUNCTION comin_descrdata_get_global_nproma

  FUNCTION comin_descrdata_get_global_wp() &
      &  BIND(C, NAME="comin_descrdata_get_global_wp") &
      &  RESULT(wp)

    INTEGER(C_INT)                      :: wp
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of wp not associated.")
    END IF
    wp  = p%wp
  END FUNCTION comin_descrdata_get_global_wp

  FUNCTION comin_descrdata_get_global_min_rlcell_int() &
      &  BIND(C, NAME="comin_descrdata_get_global_min_rlcell_int") &
      &  RESULT(min_rlcell_int)

    INTEGER(C_INT)                      :: min_rlcell_int
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of min_rlcell_int not associated.")
    END IF
    min_rlcell_int  = p%min_rlcell_int
  END FUNCTION comin_descrdata_get_global_min_rlcell_int

  FUNCTION comin_descrdata_get_global_min_rlcell() &
      &  BIND(C, NAME="comin_descrdata_get_global_min_rlcell") &
      &  RESULT(min_rlcell)

    INTEGER(C_INT)                      :: min_rlcell
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of min_rlcell not associated.")
    END IF
    min_rlcell  = p%min_rlcell
  END FUNCTION comin_descrdata_get_global_min_rlcell

  FUNCTION comin_descrdata_get_global_max_rlcell() &
      &  BIND(C, NAME="comin_descrdata_get_global_max_rlcell") &
      &  RESULT(max_rlcell)

    INTEGER(C_INT)                      :: max_rlcell
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of max_rlcell not associated.")
    END IF
    max_rlcell  = p%max_rlcell
  END FUNCTION comin_descrdata_get_global_max_rlcell

  FUNCTION comin_descrdata_get_global_min_rlvert_int() &
      &  BIND(C, NAME="comin_descrdata_get_global_min_rlvert_int") &
      &  RESULT(min_rlvert_int)

    INTEGER(C_INT)                      :: min_rlvert_int
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of min_rlvert_int not associated.")
    END IF
    min_rlvert_int  = p%min_rlvert_int
  END FUNCTION comin_descrdata_get_global_min_rlvert_int

  FUNCTION comin_descrdata_get_global_min_rlvert() &
      &  BIND(C, NAME="comin_descrdata_get_global_min_rlvert") &
      &  RESULT(min_rlvert)

    INTEGER(C_INT)                      :: min_rlvert
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of min_rlvert not associated.")
    END IF
    min_rlvert  = p%min_rlvert
  END FUNCTION comin_descrdata_get_global_min_rlvert

  FUNCTION comin_descrdata_get_global_max_rlvert() &
      &  BIND(C, NAME="comin_descrdata_get_global_max_rlvert") &
      &  RESULT(max_rlvert)

    INTEGER(C_INT)                      :: max_rlvert
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of max_rlvert not associated.")
    END IF
    max_rlvert  = p%max_rlvert
  END FUNCTION comin_descrdata_get_global_max_rlvert

  FUNCTION comin_descrdata_get_global_min_rledge_int() &
      &  BIND(C, NAME="comin_descrdata_get_global_min_rledge_int") &
      &  RESULT(min_rledge_int)

    INTEGER(C_INT)                      :: min_rledge_int
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of min_rledge_int not associated.")
    END IF
    min_rledge_int  = p%min_rledge_int
  END FUNCTION comin_descrdata_get_global_min_rledge_int

  FUNCTION comin_descrdata_get_global_min_rledge() &
      &  BIND(C, NAME="comin_descrdata_get_global_min_rledge") &
      &  RESULT(min_rledge)

    INTEGER(C_INT)                      :: min_rledge
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of min_rledge not associated.")
    END IF
    min_rledge  = p%min_rledge
  END FUNCTION comin_descrdata_get_global_min_rledge

  FUNCTION comin_descrdata_get_global_max_rledge() &
      &  BIND(C, NAME="comin_descrdata_get_global_max_rledge") &
      &  RESULT(max_rledge)

    INTEGER(C_INT)                      :: max_rledge
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of max_rledge not associated.")
    END IF
    max_rledge  = p%max_rledge
  END FUNCTION comin_descrdata_get_global_max_rledge

  FUNCTION comin_descrdata_get_global_grf_bdywidth_c() &
      &  BIND(C, NAME="comin_descrdata_get_global_grf_bdywidth_c") &
      &  RESULT(grf_bdywidth_c)

    INTEGER(C_INT)                      :: grf_bdywidth_c
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of grf_bdywidth_c not associated.")
    END IF
    grf_bdywidth_c  = p%grf_bdywidth_c
  END FUNCTION comin_descrdata_get_global_grf_bdywidth_c

  FUNCTION comin_descrdata_get_global_grf_bdywidth_e() &
      &  BIND(C, NAME="comin_descrdata_get_global_grf_bdywidth_e") &
      &  RESULT(grf_bdywidth_e)

    INTEGER(C_INT)                      :: grf_bdywidth_e
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of grf_bdywidth_e not associated.")
    END IF
    grf_bdywidth_e  = p%grf_bdywidth_e
  END FUNCTION comin_descrdata_get_global_grf_bdywidth_e

  FUNCTION comin_descrdata_get_global_lrestartrun() &
      &  BIND(C, NAME="comin_descrdata_get_global_lrestartrun") &
      &  RESULT(lrestartrun)

    LOGICAL(C_BOOL)                      :: lrestartrun
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of lrestartrun not associated.")
    END IF
    lrestartrun  = p%lrestartrun
  END FUNCTION comin_descrdata_get_global_lrestartrun

  SUBROUTINE comin_descrdata_get_global_vct_a( vct_a, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_vct_a")

    TYPE(C_PTR),    INTENT(OUT) :: vct_a
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%vct_a)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of vct_a not associated.")
    END IF
    arr_size = SHAPE(p%vct_a)
    vct_a = C_LOC(p%vct_a)
  END SUBROUTINE comin_descrdata_get_global_vct_a

  FUNCTION comin_descrdata_get_global_yac_instance_id() &
      &  BIND(C, NAME="comin_descrdata_get_global_yac_instance_id") &
      &  RESULT(yac_instance_id)

    INTEGER(C_INT)                      :: yac_instance_id
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of yac_instance_id not associated.")
    END IF
    yac_instance_id  = p%yac_instance_id
  END FUNCTION comin_descrdata_get_global_yac_instance_id

  SUBROUTINE comin_descrdata_get_global_host_git_remote_url( host_git_remote_url, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_host_git_remote_url")

    TYPE(C_PTR),    INTENT(OUT) :: host_git_remote_url
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%host_git_remote_url)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of host_git_remote_url not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%host_git_remote_url)
    host_git_remote_url = C_LOC(p%host_git_remote_url)
  END SUBROUTINE comin_descrdata_get_global_host_git_remote_url

  SUBROUTINE comin_descrdata_get_global_host_git_branch( host_git_branch, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_host_git_branch")

    TYPE(C_PTR),    INTENT(OUT) :: host_git_branch
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%host_git_branch)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of host_git_branch not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%host_git_branch)
    host_git_branch = C_LOC(p%host_git_branch)
  END SUBROUTINE comin_descrdata_get_global_host_git_branch

  SUBROUTINE comin_descrdata_get_global_host_git_tag( host_git_tag, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_host_git_tag")

    TYPE(C_PTR),    INTENT(OUT) :: host_git_tag
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%host_git_tag)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of host_git_tag not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%host_git_tag)
    host_git_tag = C_LOC(p%host_git_tag)
  END SUBROUTINE comin_descrdata_get_global_host_git_tag

  SUBROUTINE comin_descrdata_get_global_host_revision( host_revision, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_host_revision")

    TYPE(C_PTR),    INTENT(OUT) :: host_revision
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%host_revision)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of host_revision not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%host_revision)
    host_revision = C_LOC(p%host_revision)
  END SUBROUTINE comin_descrdata_get_global_host_revision

  FUNCTION comin_descrdata_get_global_has_device() &
      &  BIND(C, NAME="comin_descrdata_get_global_has_device") &
      &  RESULT(has_device)

    LOGICAL(C_BOOL)                      :: has_device
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()

    p => comin_descrdata_get_global()
    IF (.NOT. .TRUE.) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      & "ERROR: Pointer of has_device not associated.")
    END IF
    has_device  = p%has_device
  END FUNCTION comin_descrdata_get_global_has_device

  SUBROUTINE comin_descrdata_get_global_device_name( device_name, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_device_name")

    TYPE(C_PTR),    INTENT(OUT) :: device_name
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%device_name)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of device_name not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%device_name)
    device_name = C_LOC(p%device_name)
  END SUBROUTINE comin_descrdata_get_global_device_name

  SUBROUTINE comin_descrdata_get_global_device_vendor( device_vendor, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_device_vendor")

    TYPE(C_PTR),    INTENT(OUT) :: device_vendor
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%device_vendor)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of device_vendor not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%device_vendor)
    device_vendor = C_LOC(p%device_vendor)
  END SUBROUTINE comin_descrdata_get_global_device_vendor

  SUBROUTINE comin_descrdata_get_global_device_driver( device_driver, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_global_device_driver")

    TYPE(C_PTR),    INTENT(OUT) :: device_driver
    INTEGER(C_INT), INTENT(INOUT) :: arr_size(1)
    !
    TYPE(t_comin_descrdata_global), POINTER :: p => NULL()
    p => comin_descrdata_get_global()
    IF (.NOT. ALLOCATED(p%device_driver)) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_global", &
      &  "ERROR: Pointer of device_driver not associated.")
    END IF
    arr_size(1) = LEN_TRIM(p%device_driver)
    device_driver = C_LOC(p%device_driver)
  END SUBROUTINE comin_descrdata_get_global_device_driver

END MODULE comin_descrdata
