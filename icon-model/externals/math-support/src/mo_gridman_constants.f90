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

!>
!!  Module contains some constants relevant for implementational issues.
!!

MODULE mo_gridman_constants
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64

  IMPLICIT NONE

  PUBLIC

  ! define model name and version
  CHARACTER(len=*), PARAMETER :: modelname = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  INTEGER, PARAMETER :: MAX_CHAR_LENGTH = 1024

  INTEGER, PARAMETER :: SUCCESS = 0
  INTEGER, PARAMETER :: CELLS = 123
  INTEGER, PARAMETER :: EDGES = 345
  INTEGER, PARAMETER :: VERTS = 678

  INTEGER, PARAMETER :: ON_CELLS = 1
  INTEGER, PARAMETER :: ON_EDGES = 2
  INTEGER, PARAMETER :: ON_VERTICES = 3
  INTEGER, PARAMETER :: HALO_LEVELS_CEILING = 256 ! should be greater than the max level
  ! of halo levels
!-------------------------------------------------------------------------------
! Comments by Hui:
! According to Luis' explanation, the declarations above related to the blocking
! are not correct. There should be a single NPROMA, and different NBLKS values
! for edges, cells and vertices.
!
! Take the triangular cells for example. Considering that the number of cells is
! different from patch to patch, it may be a good idea NOT to declare "nblks_c" as
! a 2D array HERE, but as a component of the "patch" type. Since we want to do
! calculations ONLY for the internal cells, we also need to know how many blocks
! (e.g. "nblks_i_c") we have for the internal cells. Furthermore, it is also
! necessary to define a a parameter "npromz_i_c" for each patch, which stores the
! number of valid items in the last block of the internal cells.
!   At the very beginning of the execution of the program, "nproma" can be read
! from the namelist file. Then, after reading in the number of cells and external halo
! cells of a certain patch, we calculate the value of "nblks_c", "nblks_i_c"
! and "npromz_i_c" for that patch by:
!
!   n_patch_cell_all     = ptr_patch%ncells + ptr_patch%n_e_halo_cells
!   ptr_patch%nblks_c    = ( n_patch_cell_all - 1 )/nproma + 1
!
!   ptr_patch%nblks_i_c  = ( ptr_patch%ncell - 1 )/nproma + 1
!   ptr_patch%npromz_i_c = ptr_patch%ncell - (ptr_patch%nblks_i_c - 1)*nproma
!
!   (cf. echam5, mo_decomposition.f90,
!        the last few lines of SUBROUTINE decompose_grid)
!
! These calculations can be done in {\it mo\_model\_domain\_import} .
!-------------------------------------------------------------------------------

  ! Width of lateral boundary zones (as seen from the child domain) for which
  ! tendencies are interpolated from the parent domain
  ! These two come from mo_impl_constants_grf.f90; all others come from mo_impl_constants.f90
  INTEGER, PARAMETER :: grf_bdywidth_c = 4
  INTEGER, PARAMETER :: grf_bdywidth_e = 9

  ! dimensions of index list fields and index ranges for which grid points are reordered
  ! according to the refin_ctrl values
  ! Specifically:
  ! max_rl* indicates the number of cell/edge/vertex rows along the lateral boundary of nested
  !   domains for which grid points are reordered, i.e. moved to the beginning of the index lists;
  !   the number of cell rows for which the refin_ctrl flag is set is determined by the variable
  !   bdy_indexing_depth in prepare_gridref; it is in general larger than max_rlcell
  !  (the refin_ctrl flag here counts the distance from the lateral boundary in units of cell rows)
  ! ABS(min_rl*_int)-1 indicates the number of cell/edge/vertex rows overlapping with the lateral boundary
  !   of a nested domain for which grid points are reordered, i.e. moved to the end of the index lists;
  !   min_rl*_int refers to grid points overlapping with interior points of nested domains
  ! Finally, the indices between min_rl*_int-1 and min_rl*_int are reserved for halo points emerging
  !   from the MPI domain decomposition; these parts of the index lists are empty on exit of the
  !   grid generator and are filled only on mo_subdivision. However, the index list fields are always
  !   dimensioned with (min_rl*:max_rl*). The values set below are sufficient for a halo
  !   width of two full cell rows; normally we use one, but stencils for high-order schemes may
  !   sometime require a halo width of two full rows
!
!   ------------------------------------------
!   LL: copied from the icon_flowcontrol as described by Guenther:
!
!   The ordering of the halo points is as follows:
!   min_rlcell_int- 1: halo cells having a prognostic cell as neighbor
!   min_rlcell_int- 2: halo cells in the first cell row having no prognostic cell as neighbor
!   and analogously for the second halo cell row if present. For n_ghost_rows = 1, the index segments
!   corresponding to min_rlcell_int - 3 and min_rlcell_int - 4 (= min_rlcell) are empty.
!
!   For edges and vertices, one needs to be aware of the fact that outer boundary edges/vertices of a prognostic
!   cell may not be owned by the current PE because the PE of the neighboring cell has the ownership (otherwise
!   there would be double-counting). There are, however, operations for which even such edges/vertices can be
!   excluded from prognostic computation because a halo synchronization follows immediately afterwards (and
!   has to be there anyway). Thus, the following ordering is applied:
!   min_rledge_int - 1: outer boundary edges of a prognostic cell not owned by the current PE\\
!   min_rledge_int - 2: edges connecting halo cells of the first row
!   min_rledge_int - 3: outer boundary edges of the first halo cells row, or edges connecting cells
!   of the first halo cell row with cells of the second halo cell row.
!   For n_ghost_rows = 2, an analogous setting applies to min_rledge_int - 4 and
!   min_rledge_int - 5 (= min_rledge). For vertices, we have
!   min_rlvert_int - 1: outer boundary vertices of a prognostic cell not owned by the current PE
!   min_rlvert_int - 2: outer boundary vertices of the first halo cells row, or vertices connecting cells
!   of the first halo cell row with cells of the second halo cell row.
!   For n_ghost_rows = 2, an analogous setting applies to min_rlvert_int - 3  (= min_rlvert).
  !---------------------------------------------
  !
  ! Ordering Scheme:
  !
  ! Following is the order of the grid entities (an all the associated variables)
  ! in ascending order.
  !
  ! A. The indices from 1 to max_rl
  !    Mark the lateral boundaries of the patch
  !    start_idx(1) = start of the first boundary level. It is always 1
  !    end_idx(1)   = end of the first boundary level end_idx(1)
  !    start_idx(2) = start of the second boundary level, it is always end_idx(1)+1
  !    ..... etc until end_idx(max_rl) which the the end of the lateral boundary levels
  !
  ! B. The index 0
  !    Marks the internal entities, that do not overlap with child patches
  !    start_idx(0) = end_idx(max_rl) + 1
  !    end_idx(0)   = start_idx(-1) -1
  !
  ! C. The indices from -1 to min_rl_int
  !    Mark the internal entities that overlap with child patches
  !    (they are defined for each child patch)
  !    start_idx(-1) = start of the internal entities overlapping with the first (two) levels
  !                    of the lateral boundaries of the child patch
  !    end_idx(-1)   = end of the internal entities overlapping with the first (two) levels
  !                    of the lateral boundaries of the child patch
  !    start_idx(-2) = start of the internal entities overlapping with the next (two) levels
  !                    of the lateral boundaries of the child patch = end_idx(-1) + 1
  !    ... etc
  !    end_idx(minrl_int) = end of all the internal entities overlapping with the the child patch
  !
  ! D. The indices from min_rl_int-1 to min_rl
  !    Mark the halo entities, when they do not overlap with a child patch
  !    Note: See above
  !
  !---------------------------------------------
  !
  ! Examples:
  !  A. Get all entities in the grid:    start_idx(1) -- end_idx(min_rl)
  !     This is the default range for most operators
  !  B. Get all owned entities: start_idx(1) -- end_idx(min_rl_int)
  !     Note that this may still contain halo entities if they overlap with child patches
  !  C. Get all entities, except halos: for cells: start_idx(1) -- end_idx(min_rl_int)
  !        For verts/edges: start_idx(1) -- end_idx(min_rl_int - 1)
  !
  !---------------------------------------------

  INTEGER, PARAMETER :: max_hw = 2 ! maximum halo width (n_ghost_rows)
  !
  INTEGER, PARAMETER :: min_rlcell_int = -4 ! previously -6
  INTEGER, PARAMETER :: min_rlcell = min_rlcell_int - 2*max_hw ! = -8
  INTEGER, PARAMETER :: max_rlcell = 5 ! previously 8
  INTEGER, PARAMETER :: min_rlvert_int = min_rlcell_int
  INTEGER, PARAMETER :: min_rlvert = min_rlvert_int - (max_hw + 1)
  INTEGER, PARAMETER :: max_rlvert = max_rlcell
  INTEGER, PARAMETER :: min_rledge_int = 2*min_rlcell_int ! -8
  INTEGER, PARAMETER :: min_rledge = min_rledge_int - (2*max_hw + 1) ! -13
  INTEGER, PARAMETER :: max_rledge = 2*max_rlcell ! 10

  ! Aliases for loop control
  !
  ! start index level for prognostic cells (excluding a possible nest boundary interpolation zone)
  INTEGER, PARAMETER :: start_prog_cells = grf_bdywidth_c + 1
  ! start index level for prognostic edges (excluding a possible nest boundary interpolation zone)
  INTEGER, PARAMETER :: start_prog_edges = grf_bdywidth_e + 1
  ! remark: the corresponding parameter for vertices would be unused
  !
  ! end index level for computations excluding all halo cells
  INTEGER, PARAMETER :: end_prog_cells = min_rlcell_int
  ! end index level for computations including halo level 1 cells (direct neighbors, i.e. halo cells
  ! sharing an edge with a prognostic cell)
  INTEGER, PARAMETER :: end_halo_lev1_cells = min_rlcell_int - 1
  ! end index for computations including all halo cells (this adds indirect neighbors, i.e.
  ! halo cells sharing only a vertex with a prognostic cell)
  ! remark: n_ghost_rows=2, i.e. using a full second row of halo cells, is currently not foreseen in the code,
  ! so an end index level of min_rlcell_int - 2 is equivalent to min_rlcell
  INTEGER, PARAMETER :: end_all_cells = min_rlcell
  !
  ! end index level for computations excluding all halo edges
  INTEGER, PARAMETER :: end_prog_edges = min_rledge_int
  ! end index level for edges of prognostic cells (including those not owned by the current PE)
  INTEGER, PARAMETER :: end_edges_of_prog_cells = min_rledge_int - 1
  ! end index level for computations including all edges of halo level 1 cells
  INTEGER, PARAMETER :: end_halo_lev1_edges = min_rledge_int - 2
  ! end index for computations including all halo edges (remark: consistent with what has been said above,
  ! min_rledge is equivalent to min_rledge_int - 3)
  INTEGER, PARAMETER :: end_all_edges = min_rledge
  !
  ! end index level for computations excluding all halo verices
  INTEGER, PARAMETER :: end_prog_verts = min_rlvert_int
  ! end index level for vertices of prognostic cells (including those not owned by the current PE)
  INTEGER, PARAMETER :: end_verts_of_prog_cells = min_rlvert_int - 1
  ! end index for computations including all halo vertices (remark: consistent with what has been said above,
  ! min_rledge is equivalent to min_rlvert_int - 2)
  INTEGER, PARAMETER :: end_all_verts = min_rlvert

  ! maximum allowed number of model domains (10 should be enough for the time being)
  INTEGER, PARAMETER :: max_dom = 10
  ! maximum number of decimal digits to print domain id, roughly(log10(max_dom))
  INTEGER, PARAMETER :: max_dom_dig10 = 2
  ! Maximum allowed number of physical model domains
  INTEGER, PARAMETER :: max_phys_dom = 30

  ! maximum allowed number of tracers (20 should be enough for the time being)
  ! DRIEG: For ART, more than 20 tracers are needed
  ! For ICON-waves the minimum value is 900
  INTEGER, PARAMETER :: max_ntracer = 1600

  ! maximum allowed number of echotop levels:
  INTEGER, PARAMETER :: max_echotop = 10

  ! maximum allowed number of wshear levels:
  INTEGER, PARAMETER :: max_wshear = 10

  ! maximum allowed number of srh levels:
  INTEGER, PARAMETER :: max_srh = 10

  ! The lon-lat parameterization of the torus is
  !    (lon,lat) = [0, 2*pi] x [-max_lat, max_lat]
  ! where max_lat := pi/180 = 10 degrees
  ! (hard-coded in the torus grid generator)
  REAL(wp), PARAMETER :: TORUS_MAX_LAT = 4._wp/18._wp*ATAN(1._wp)

!--------------------------------------------------------------------
END MODULE mo_gridman_constants
