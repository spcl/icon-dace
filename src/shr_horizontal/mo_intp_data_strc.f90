! Contains the the interpolation data structures for the triangular grid.
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

MODULE mo_intp_data_strc
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, max_dom, SUCCESS,     &
    &                               HINTP_TYPE_NONE, HINTP_TYPE_LONLAT_RBF,                   &
    &                               HINTP_TYPE_LONLAT_NNB, HINTP_TYPE_LONLAT_BCTR
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_lonlat_grid,         ONLY: t_lon_lat_grid, OPERATOR(==)
  USE mo_communication,       ONLY: t_comm_gather_pattern
  USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_dim_c2l, l_mono_c2l

  IMPLICIT NONE


  ! NOTE: The variables will be use along the mo_interpolation sub-modules
  !       They are declared to be public
  PUBLIC
  
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_data_strc'
  


  TYPE t_lsq
    ! fields related to weighted least squares polynomial reconstruction
    !---------------------------------------------------------------------
    INTEGER, ALLOCATABLE  :: lsq_dim_stencil(:,:)  ! stencil size as a function of jc and jb
                                                   ! necessary in order to account for pentagon
                                                   ! points.
    INTEGER, ALLOCATABLE  :: lsq_idx_c(:,:,:)      ! index array defining the stencil for
                                                   ! lsq reconstruction (nproma,nblks_c,lsq_dim_c)
    INTEGER, ALLOCATABLE  :: lsq_blk_c(:,:,:)      ! dito for the blocks
    REAL(wp), ALLOCATABLE :: lsq_weights_c(:,:,:)  ! weights for lsq reconstruction
                                                   ! (nproma,lsq_dim_c,nblks_c)
    REAL(wp), ALLOCATABLE :: lsq_qtmat_c(:,:,:,:)  ! transposed Q of QR-factorization for
                                                   ! lsq reconstruction
                                                   ! (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), ALLOCATABLE :: lsq_rmat_rdiag_c(:,:,:)! reciprocal diagonal elements of R-matrix
                                                    ! resulting from QR-decomposition
                                                    ! (nproma,lsq_dim_unk,nblks_c)
    REAL(wp), ALLOCATABLE :: lsq_rmat_utri_c(:,:,:)! upper triangular elements without diagonal
                                                   ! elements of R-matrix (starting from the bottom
                                                   ! right)
                                                   ! (nproma,(lsq_dim_unk^2-lsq_dim_unk)/2,nblks_c)
    REAL(wp), ALLOCATABLE :: lsq_pseudoinv(:,:,:,:)! pseudo (or Moore-Penrose) inverse of lsq
                                                   ! design matrix A
                                                   ! (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
    REAL(wp), ALLOCATABLE :: lsq_moments(:,:,:)      ! Moments (x^ny^m)_{i} for control volume
                                                     ! (nproma,nblks_c,lsq_dim_unk)
    REAL(wp), ALLOCATABLE :: lsq_moments_hat(:,:,:,:)! Moments (\hat{x^ny^m}_{ij}) for control volume
                                                     ! (nproma,nblks_c,lsq_dim_c,lsq_dim_unk)
  END TYPE t_lsq


  TYPE t_gauss_quad
    !
    ! quadrature points for intergration over triangular element
    !
    TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< linear (nproma, nblks_c)
      &  qpts_tri_l(:,:)
    TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< quadratic (nproma, nblks_c,3)
      &  qpts_tri_q(:,:,:)
    TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< cubic (nproma, nblks_c,4)
      &  qpts_tri_c(:,:,:)
    REAL(wp), ALLOCATABLE :: weights_tri_q(:)      !< quadratic weight (3)
    REAL(wp), ALLOCATABLE :: weights_tri_c(:)      !< cubic weights (4)
  END TYPE t_gauss_quad


  TYPE t_cell_environ
    !
    ! This derived type stores an index list of cells lying within a certain radius
    ! around a center cell (ic,ib) 
    ! and additionally limited by the number of halo lines.
    ! The related variable is defined in subr. 'gen_index_list_radius'.
    !

    LOGICAL   :: is_used = .FALSE.      ! currently: set .TRUE. if SDI or LPI shall be computed

    INTEGER   :: nmbr_nghbr_cells_alloc = 13  ! max. number of cells around *any* center cell for allocation
                                ! This value must be increased for larger numbers of max_nmbr_iter or radius
                                ! (see subr. gen_index_list_radius)

    INTEGER   :: max_nmbr_nghbr_cells   ! max. number of cells around *any* center cell
                                        ! (including the center cell)

    INTEGER, ALLOCATABLE :: nmbr_nghbr_cells(:,:)  ! number of cells around the center cell (ic,ib)
                                                   ! (including the center cell).
         ! Recommendation: for efficient vectorization, build loops using max_nmbr_nghbr_cells 
         ! instead of this field.

    INTEGER,  ALLOCATABLE :: idx(:,:,:)             ! Index jc of the l-th neighbour cell (ic, ib, l)

    INTEGER,  ALLOCATABLE :: blk(:,:,:)             ! Block jb of the l-th neighbour cell (ic, ib, l)

    REAL(wp), ALLOCATABLE :: area_norm(:,:,:)       ! area of the l-th neighbour cell (ic, ib, l)
                                                    ! normalized by the area sum over all neighbour cells.
                                                    ! Note: area_norm=0 if l>nmbr_nghbr_cells(ic,ib)

    ! the following variables are only for validation purposes in any calling subroutine
    REAL(wp)             :: radius                 ! the limiting radius
    INTEGER              :: max_nmbr_iter          ! maximum number of iterations used

  END TYPE t_cell_environ
 


  TYPE t_int_state
  
    ! a) weights which are inconsistent with the Hamiltonian viewpoint
    !-----------------------------------------------------------------
  
    REAL(wp), ALLOCATABLE :: c_lin_e(:,:,:)   ! coefficient for interpolation
                                              ! from adjacent cells onto edge
                                              ! (nproma,2,nblks_e)
  
    REAL(wp), ALLOCATABLE :: e_bln_c_s(:,:,:) ! coefficient for bilinear
                                              ! interpolation from edges to cells
                                              ! for scalar quantities
  
    REAL(wp), ALLOCATABLE :: e_bln_c_u(:,:,:) ! coefficient for bilinear interpolation
                                              ! from edges to cells for vector components
                                              ! (input: v_t, v_n, output: u)
  
    REAL(wp), ALLOCATABLE :: e_bln_c_v(:,:,:) ! coefficient for bilinear interpolation
                                              ! from edges to cells for vector components
                                              ! (input: v_t, v_n, output: v)
  
    REAL(wp), ALLOCATABLE :: c_bln_avg(:,:,:) ! coefficients for bilinear divergence
                                              ! averaging (nproma,4,nblks_c)
  
    REAL(wp), ALLOCATABLE :: e_flx_avg(:,:,:) ! coefficients for related velocity or mass flux
                                              ! averaging (nproma,5,nblks_e)
  
    REAL(wp), ALLOCATABLE :: v_1o2_e(:,:,:)   ! coefficient for interpolation
                                              ! from vertices onto edges by 1/2
                                              ! weighting (nproma,2,nblks_e),

    REAL(wp), ALLOCATABLE :: gradc_bmat(:,:,:,:) ! Bmatrix for cell centered shape function based
                                              ! gradient (nproma,2,3,nblks_c)
  
  
    ! b) weights which are consistent with the Hamiltonian viewpoint
    !---------------------------------------------------------------
    ! The following weights are needed for the mass and theta brackets
  
    REAL(wp), ALLOCATABLE :: e_inn_c(:,:,:)   ! coefficient for inner product
                                              ! of 2 vector components
                                              ! from edges to cells

    REAL(wp), ALLOCATABLE :: verts_aw_cells(:,:,:)! coefficient for interpolation
                                              ! from vertices to cells by
                                              ! area weighting
  
    REAL(wp), ALLOCATABLE :: cells_aw_verts(:,:,:)! coefficient for interpolation
                                              ! from cells to verts by
                                              ! area weighting
  
    ! c) RBF related fields
    !----------------------
    INTEGER, ALLOCATABLE  :: rbf_vec_idx_c(:,:,:)  ! index array defining the
                                              ! stencil of surrounding edges for
                                              ! vector rbf interpolation at each
                                              ! cell center
                                              ! (rbf_vec_dim_c,nproma,nblks_c)
    INTEGER, ALLOCATABLE  :: rbf_vec_blk_c(:,:,:)  ! ... dito for the blocks
  
    INTEGER, ALLOCATABLE  :: rbf_vec_stencil_c(:,:) ! array defining number of
                                              ! surrounding edges in the stencil
                                              ! for vector rbf interpolation at
                                              ! each cell center
                                              ! (nproma,nblks_c)
    REAL(wp), ALLOCATABLE :: rbf_vec_coeff_c(:,:,:,:) ! array containing the
                                              ! coefficients used for
                                              ! vector rbf interpolation
                                              ! at each cell center
                                              ! (rbf_vec_dim_c,2,nproma,nblks_c)
  
    INTEGER, ALLOCATABLE  :: rbf_c2grad_idx(:,:,:)  ! index array defining the
                                              ! stencil of surrounding cells for
                                              ! 2D gradient reconstruction at each
                                              ! cell center
                                              ! (rbf_c2grad_dim,nproma,nblks_c)
    INTEGER, ALLOCATABLE  :: rbf_c2grad_blk(:,:,:)  ! ... dito for the blocks
  
    REAL(wp), ALLOCATABLE :: rbf_c2grad_coeff(:,:,:,:) ! array containing the
                                              ! coefficients used for
                                              ! 2D gradient reconstruction
                                              ! at each cell center
                                              ! (rbf_c2grad_dim,2,nproma,nblks_c)
  
    INTEGER, ALLOCATABLE  :: rbf_vec_idx_v(:,:,:) ! index array defining the
                                              ! stencil of surrounding edges for
                                              ! vector rbf interpolation at each
                                              ! triangle vertex
                                              ! (rbf_vec_dim_v,nproma,nblks_v)
    INTEGER, ALLOCATABLE  :: rbf_vec_blk_v(:,:,:) ! ... dito for the blocks
  
    INTEGER, ALLOCATABLE  :: rbf_vec_stencil_v(:,:) ! array defining number of
                                              ! surrounding edges in the stencil
                                              ! for vector rbf interpolation at
                                              ! each triangle vertex
                                              ! (nproma,nblks_v)
  
    REAL(wp), ALLOCATABLE :: rbf_vec_coeff_v(:,:,:,:) ! array containing the
                                              ! coefficients used for vector rbf
                                              ! interpolation at each tringle
                                              ! vertex (input is normal component)
                                              ! (rbf_vec_dim_v,2,nproma,nblks_v)
  
    INTEGER, ALLOCATABLE  :: rbf_vec_idx_e(:,:,:) ! index array defining the
                                              ! stencil of surrounding edges for
                                              ! vector rbf interpolation at each
                                              ! triangle edge
                                              ! (rbf_vec_dim_e,nproma,nblks_e)
    INTEGER, ALLOCATABLE  :: rbf_vec_blk_e(:,:,:) ! ... dito for the blocks
  
    INTEGER, ALLOCATABLE  :: rbf_vec_stencil_e(:,:) ! array defining number of
                                              ! surrounding edges in the stencil
                                              ! for vector rbf interpolation at
                                              ! each triangle edge
                                              ! (nproma,nblks_e)
  
    REAL(wp), ALLOCATABLE :: rbf_vec_coeff_e(:,:,:) ! array containing the
                                              ! coefficients used for rbf inter-
                                              ! polation of the tangential velo-
                                              ! city component (from the
                                              ! surrounding normals) at each
                                              ! triangle edge
                                              ! (rbf_vec_dim_e,nproma,nblks_e)

    ! d) precomputed geometrical factors for mathematical operators (for efficiency)
    !------------------------------------------------------------------------------
    REAL(wp), ALLOCATABLE :: geofac_div(:,:,:)    ! factor for divergence (nproma,cell_type,nblks_c)
    REAL(wp), ALLOCATABLE :: geofac_grdiv(:,:,:)  ! factor for gradient of divergence (nproma,5,nblks_e)
    REAL(wp), ALLOCATABLE :: geofac_rot(:,:,:)    ! factor for divergence (nproma,9-cell_type,nblks_v)
    REAL(wp), ALLOCATABLE :: geofac_n2s(:,:,:)    ! factor for nabla2-scalar (nproma,cell_type+1,nblks_c)
    REAL(wp), ALLOCATABLE :: geofac_grg(:,:,:,:)  ! factor for Green-Gauss gradient (nproma,4,nblks_c,2)
  
    ! e) precomputed Cartesian orientation and location vectors of edge midpoints
    !    and location of cell centers(for efficiency) : it is now computed in grid generator stored
    !    in p_patch
    !------------------------------------------------------------------------------
  
    ! f) patch elements restored from edges to cells to reduce frequency of indirect addressing
    !------------------------------------------------------------------------------
    REAL(wp), ALLOCATABLE :: primal_normal_ec(:,:,:,:) ! p_patch%edges%primal_normal_cell stored on
                                                       ! the cell data type (nproma,nblks_c,3,2)
    REAL(wp), ALLOCATABLE :: edge_cell_length(:,:,:)   ! p_patch%edges%edge_cell_length stored on
                                                       ! the cell data type (nproma,nblks_c,3)
  
    ! g) distance from cells to vertices on local cartesian grid with origin at the cell center
    !    (used for gradient limiter)
    !------------------------------------------------------------------------------
    REAL(wp), ALLOCATABLE :: cell_vert_dist(:,:,:,:)   ! (nproma,3,2,nblks_c)
  
    ! h) fields related to calculation of backward trajectories on local plane
    !    tangential to the edge midpoint
    !------------------------------------------------------------------------------
    REAL(wp), ALLOCATABLE :: pos_on_tplane_e(:,:,:,:)  ! positions of various points on local plane
                                                       ! tangential to the edge midpoint.
                                                       ! (nproma,4,2,nblks_e)

    TYPE(t_geographical_coordinates), ALLOCATABLE ::  &! positions of vertices and butterfly
      &  pos_on_tplane_c_edge(:,:,:,:)                 ! neighbors on local plane tangential to the
                                                       ! cell circumcenter.
                                                       ! stored in an edge-based data structure
                                                       ! (nproma,nblks_e,2,5)
  
  
    ! i) fields related to weighted least squares polynomial reconstruction
    !------------------------------------------------------------------------------
    TYPE(t_lsq) :: lsq_lin,  &  ! coefficients for linear lsq-reconstruction
      &            lsq_high     ! coefficients for higher order lsq-reconstruction


    ! j) Nudging coefficients used for 1-way nesting and limited-area mode (defined here
    !    rather than in grf_state because the limited-area mode may be used without nesting)
    !------------------------------------------------------------------------------
    REAL(wp), ALLOCATABLE :: nudgecoeff_c(:,:)  !< Nudging coefficient for cells
    REAL(wp), ALLOCATABLE :: nudgecoeff_e(:,:)  !< Nudging coefficient for cells
  
  
    ! k) Quadrature points and weights for integration over triangular element
    !--------------------------------------------------------------------------
    TYPE(t_gauss_quad) ::gquad

    ! index list for neighbouring cells within a certain radius
    TYPE(t_cell_environ) :: cell_environ
 
  END TYPE t_int_state


  ! MODULE VARIABLES --------------------------------------------------------------

  TYPE(t_int_state),TARGET,ALLOCATABLE :: p_int_state(:), p_int_state_local_parent(:)  

END MODULE mo_intp_data_strc
