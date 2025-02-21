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
  SUBROUTINE comin_current_get_datetime_c(val, len, ierr) &
    & BIND(C, NAME="comin_current_get_datetime")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
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
  SUBROUTINE comin_descrdata_get_simulation_interval_exp_start(val, len, ierr) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_exp_start")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
    val = C_LOC(state%comin_descrdata_simulation_interval%exp_start)
    len = LEN_TRIM(state%comin_descrdata_simulation_interval%exp_start)
  END SUBROUTINE comin_descrdata_get_simulation_interval_exp_start


  !> request simulation interval information, C interface
  SUBROUTINE comin_descrdata_get_simulation_interval_exp_stop(val, len, ierr) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_exp_stop")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
    val = C_LOC(state%comin_descrdata_simulation_interval%exp_stop)
    len = LEN_TRIM(state%comin_descrdata_simulation_interval%exp_stop)
  END SUBROUTINE comin_descrdata_get_simulation_interval_exp_stop


  !> request simulation interval information, C interface
  SUBROUTINE comin_descrdata_get_simulation_interval_run_start(val, len, ierr) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_run_start")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
    val = C_LOC(state%comin_descrdata_simulation_interval%run_start)
    len = LEN_TRIM(state%comin_descrdata_simulation_interval%run_start)
  END SUBROUTINE comin_descrdata_get_simulation_interval_run_start


  !> request simulation interval information, C interface
  SUBROUTINE comin_descrdata_get_simulation_interval_run_stop(val, len, ierr) &
    & BIND(C, NAME="comin_descrdata_get_simulation_interval_run_stop")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< datetime string (ISO 8601)
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
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
  SUBROUTINE comin_descrdata_set_timesteplength(jg, dt_current, ierr) BIND(C)
    INTEGER(C_INT), VALUE, INTENT(IN)  :: jg
    REAL(wp),              INTENT(IN)  :: dt_current
    INTEGER(C_INT),        INTENT(OUT) :: ierr

    ierr = COMIN_SUCCESS ! currently unused
    IF (.NOT. ALLOCATED(state%comin_descrdata_timesteplength)) THEN
      ALLOCATE(state%comin_descrdata_timesteplength(state%comin_descrdata_global%n_dom+4))
    END IF
    state%comin_descrdata_timesteplength(jg) = dt_current
  END SUBROUTINE comin_descrdata_set_timesteplength

  !> Clean descriptive data structure in ComIn
  !> currently no content but keep for future use
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_finalize(ierr)
    INTEGER, INTENT(OUT) :: ierr !< error code

    ierr = COMIN_SUCCESS ! currently unused: COMIN_ERROR_DESCRDATA_FINALIZE
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
  INCLUDE 'comin_descrdata_query_domain.F90'

  ! Query global data routines generated by python script (comin_descrdata_get_global.F90.py) in ../utils. !
  INCLUDE 'comin_descrdata_query_global.F90'

END MODULE comin_descrdata
