! Contains the implementation of interpolation onto regular grids.
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









  MODULE mo_intp_lonlat
    !-------------------------------------------------------------------------
    !
    !    ProTeX FORTRAN source: Style 2
    !    modified for ICON project, DWD/MPI-M 2006
    !
    !-------------------------------------------------------------------------
    !



    USE mo_kind,                ONLY: wp
    USE mo_exception,           ONLY: message, message_text, finish
    USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int,                                &
      &                               HINTP_TYPE_LONLAT_RBF,                                  &
      &                               HINTP_TYPE_LONLAT_NNB, HINTP_TYPE_LONLAT_BCTR,          &
      &                               SCALE_MODE_TABLE, SCALE_MODE_AUTO, SCALE_MODE_PRESET
    USE mo_cdi_constants,       ONLY: GRID_REGULAR_LONLAT, GRID_CELL
    USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
    USE mo_model_domain,        ONLY: t_patch
    USE mo_run_config,          ONLY: timers_level
    USE mo_grid_config,         ONLY: n_dom, grid_sphere_radius, is_plane_torus, l_limited_area
    USE mo_timer,               ONLY: timer_start, timer_stop, timer_lonlat_setup
    USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
    USE mo_math_utilities,      ONLY: gc2cc, gvec2cvec, arc_length_v
    USE mo_math_utility_solvers, ONLY: solve_chol_v, choldec_v
    USE mo_lonlat_grid,         ONLY: t_lon_lat_grid, rotate_latlon_grid
    USE mo_parallel_config,     ONLY: nproma, p_test_run
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
    USE mo_intp_data_strc,      ONLY: t_int_state
    USE mo_intp_lonlat_types,   ONLY: t_lon_lat_intp, lonlat_grids, t_intp_scalar_coeff,      &
      &                               dbg_level
    USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_c2grad_dim, rbf_vec_kern_ll,         &
      &                               rbf_vec_scale_ll, rbf_dim_c2l,                          &
      &                               rbf_scale_mode_ll, support_baryctr_intp,    &
      &                               lreduced_nestbdry_stencil
    USE mo_gnat_gridsearch,     ONLY: gnat_init_grid, gnat_destroy, t_gnat_tree,              &
      &                               gnat_query_containing_triangles,                        &
      &                               gnat_merge_distributed_queries, gk, SKIP_NODE,          &
      &                               INVALID_NODE, gnat_recursive_proximity_query
    USE mo_mpi,                 ONLY: get_my_mpi_work_id,                                     &
      &                               p_min, p_bcast, p_comm_work, p_io, p_pe_work,              &
      &                               my_process_is_mpi_test, p_max, p_send,                  &
      &                               p_recv, process_mpi_all_test_id,                        &
      &                               process_mpi_all_workroot_id, my_process_is_stdio
    USE mo_communication,       ONLY: blk_no, idx_no, setup_comm_gather_pattern, idx_1d
    USE mo_rbf_errana,          ONLY: estimate_rbf_parameter
    USE mo_delaunay_types,      ONLY: t_point_list, OPERATOR(/), t_triangulation, OPERATOR(+)
    USE mo_intp_lonlat_baryctr, ONLY: setup_barycentric_intp_lonlat,                          &
      &                               setup_barycentric_intp_lonlat_repartition,              &
      &                               compute_auxiliary_triangulation, try_triangulation_readin
    USE mo_util_string,         ONLY: int2string

    IMPLICIT NONE

    !> module name
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_lonlat'

    PRIVATE

    PUBLIC :: compute_lonlat_intp_coeffs

  CONTAINS

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

!-------------------------------------------------------------------------   
!>
!! Gaussian kernel for RBF interpolation.
!! 
!! @f$\phi(r)=e^{-r^2}@f$
!! 
FUNCTION gaussi (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = -1._wp * p_rbf_val * p_rbf_val
p_rbf_val = EXP(p_rbf_val)

END FUNCTION gaussi



!-------------------------------------------------------------------------
!
!  
!>
!! Multiquadric kernel for RBF interpolation.
!! 
!! @f$\phi(r)=(1+r^2)^{\frac{1}{2}}@f$
!! 
FUNCTION multiq (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = p_rbf_val * p_rbf_val
p_rbf_val = SQRT(1._wp + p_rbf_val)

END FUNCTION multiq

!-------------------------------------------------------------------------
!
!  
!>
!! Inverse multiquadric kernel for RBF interpolation.
!! 
!! @f$\phi(r)=(1+r^2)^{-\frac{1}{2}}@f$
!! 
FUNCTION inv_multiq (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = p_rbf_val * p_rbf_val
p_rbf_val = SQRT(1._wp + p_rbf_val)
p_rbf_val = 1._wp / p_rbf_val

END FUNCTION inv_multiq

!-------------------------------------------------------------------------
!
!  
!>
!! Weighting function for IDW interpolation.
!! 
!! @f$\phi(r)=\frac{1}{(r)^p_{exp}}@f$
!! 
FUNCTION invdwgt (p_x, p_exp)  RESULT (wgtfac)
REAL(wp) , INTENT(in) :: p_x            ! radial distance
REAL(wp) , INTENT(in) :: p_exp          ! exponent

REAL(wp)              :: wgtfac       ! weighting factor

!-----------------------------------------------------------------------  

  wgtfac = 1._wp/max(1.e-10_wp , p_x**p_exp)

END FUNCTION invdwgt


!-------------------------------------------------------------------------  
!
!  
!>
!! Modified inverse multiquadric kernel for RBF interpolation.
!! 
!! @f$\phi(r)=(1+r^2)^{-1}@f$
!! 
FUNCTION inv_multiq2 (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = p_rbf_val * p_rbf_val
p_rbf_val = 1._wp / (1._wp + p_rbf_val)

END FUNCTION inv_multiq2



    !---------------------------------------------------------------
    !> Setup of lon-lat interpolation
    !
    SUBROUTINE compute_lonlat_intp_coeffs(p_patch, p_int_state)
      TYPE(t_patch),        INTENT(INOUT) :: p_patch(:)
      TYPE(t_int_state),    INTENT(INOUT) :: p_int_state(:)
      ! local variables
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_lonlat_intp_coeffs"
      LOGICAL,      PARAMETER :: lrepartition = .TRUE.

      INTEGER                :: i, jg, n_points, my_id, nthis_local_pts, &
        &                       triangulation_read_from_file, ist
      TYPE(t_gnat_tree)      :: gnat
      TYPE (t_triangulation) :: tri_global
      TYPE (t_point_list)    :: p_global

      ! list of triangles containing lon-lat grid points (first dim: index and block)
      INTEGER, ALLOCATABLE  :: tri_idx(:,:,:) ! 2, nproma, nblks_lonlat


      ! -----------------------------------------------------------

      IF (dbg_level > 5) CALL message(routine, "Enter")
      ! allocate and compute coefficients needed for lon-lat
      ! interpolation:
      DO jg=1,n_dom
        IF (.NOT. ANY(lonlat_grids%list(1:lonlat_grids%ngrids)%l_dom(jg))) CYCLE

        ! build GNAT data structure
        IF (dbg_level > 1) CALL message(routine, "build GNAT data structure")
        CALL gnat_init_grid(gnat, p_patch(jg))

        DO i=1, lonlat_grids%ngrids
          IF (lonlat_grids%list(i)%l_dom(jg)) THEN
            IF (timers_level > 3) CALL timer_start(timer_lonlat_setup)

            ! allocate global arrays for distributed computation:
            ALLOCATE(tri_idx(2, nproma, lonlat_grids%list(i)%grid%nblks),                           &
              &      lonlat_grids%list(i)%intp(jg)%global_idx(lonlat_grids%list(i)%grid%total_dim), &
              &      STAT=ist )
            IF (ist /= SUCCESS) CALL finish (routine, 'allocation for lon-lat point distribution failed')

            ! for each PE determine lon-lat points 
            CALL distribute_lonlat_points(lonlat_grids%list(i)%grid,     &
              &                           gnat, p_patch(jg), tri_idx,   &
              &                           lonlat_grids%list(i)%intp(jg) )

            ! After the lon-lat points have been distributed, we know
            ! exactly for each PE the number of points
            CALL lonlat_grids%list(i)%intp(jg)%init(nproma)

            ! compute weights for Radial Basis Function (RBF)
            ! interpolation:
            CALL rbf_setup_interpol_lonlat_grid(p_patch(jg), lonlat_grids%list(i)%intp(jg), &
              &                                 tri_idx, p_int_state(jg) )

            ! compute "weights" for Nearest-Neighbor interpolation:
            CALL nnb_setup_interpol_lonlat_grid(tri_idx, lonlat_grids%list(i)%intp(jg))

            ! compute weights for barycentric interpolation:
            IF (support_baryctr_intp) THEN
              
              ! --- try to read auxiliary triangulation from file:
              triangulation_read_from_file = try_triangulation_readin(p_patch(jg), tri_global, p_global)

              ! --- compute barycentric weights
              IF (((jg == 1) .AND. lrepartition) .OR. &
                &  (triangulation_read_from_file == SUCCESS)) THEN

                CALL setup_barycentric_intp_lonlat_repartition(p_patch(jg), tri_global, p_global, &
                  &                                            triangulation_read_from_file,      &
                  &                                            lonlat_grids%list(i)%intp(jg))

              ELSE

                CALL compute_auxiliary_triangulation(p_patch(jg), tri_global, p_global)
                CALL setup_barycentric_intp_lonlat(tri_global, p_global, lonlat_grids%list(i)%intp(jg))

              END IF

              CALL p_global%destructor()
              CALL tri_global%destructor()

              ! The inexact "point-inside-triangle" test which is used
              ! by the barycentric interpolation algorithm may cause
              ! (very small) weights outside [0,1]. We therefore apply
              ! a cut-off and renormalization.
              CALL normalize_weights(lonlat_grids%list(i)%intp(jg), &
                &                    lonlat_grids%list(i)%intp(jg)%baryctr)
            END IF

            ! --- for debugging purposes: visualize stencil for a
            !   - given lon-lat coordinate (writes in *.asy format)
            ! CALL visualize_stencil(lonlat_grids%list(i)%intp(jg), p_patch(jg), jg, 51.5_wp, 40.5_wp)

            ! mask out nest boundary points
            !
            ! Tentative fix: "lreduced_nestbdry_stencil" controls
            ! only barycentric interpolation.
            !
            !            IF ((jg > 1) .AND. lreduced_nestbdry_stencil) THEN 
            IF ((jg > 1) .OR. l_limited_area) THEN 
              CALL mask_out_boundary( p_patch(jg), lonlat_grids%list(i)%intp(jg) )
            END IF

            IF (timers_level > 3) CALL timer_stop(timer_lonlat_setup)
            lonlat_grids%list(i)%intp(jg)%l_initialized = .TRUE.

            ! clean-ip:
            DEALLOCATE(tri_idx, STAT=ist )
            IF (ist /= SUCCESS) CALL finish (routine, "Deallocate failed!")
          END IF
        END DO ! i

        ! clean up
        CALL gnat_destroy(gnat)
      END DO ! jg

      ! --- set communication pattern for gathering lon-lat grids:

      DO i=1, lonlat_grids%ngrids
        DO jg=1,n_dom
          IF (lonlat_grids%list(i)%l_dom(jg)) THEN
            ! n_points     : total number of points in the RECEIVER array
            n_points = lonlat_grids%list(i)%grid%lon_dim * &
              &        lonlat_grids%list(i)%grid%lat_dim

            my_id = get_my_mpi_work_id()
            nthis_local_pts = lonlat_grids%list(i)%intp(jg)%nthis_local_pts

            ! resize global data structures, after the setup only
            ! local information is needed:
            CALL lonlat_grids%list(i)%intp(jg)%contract()
            
            CALL setup_comm_gather_pattern(n_points,                       &
              (/(my_id, i = 1, nthis_local_pts)/),                         &
              lonlat_grids%list(i)%intp(jg)%global_idx(1:nthis_local_pts), &
              lonlat_grids%list(i)%p_pat(jg),.true.)

          END IF
        END DO ! jg

      END DO ! i
      IF (dbg_level > 5) CALL message(routine, "Done")

    END SUBROUTINE compute_lonlat_intp_coeffs


    !> For debugging purposes: visualize stencil at a
    !! given lon-lat coordinate (writes in *.asy format).
    !!
    SUBROUTINE visualize_stencil(ptr_int_lonlat, ptr_patch, jg, lon,lat)
      TYPE(t_lon_lat_intp), INTENT(inout) :: ptr_int_lonlat
      TYPE(t_patch),        INTENT(IN)    :: ptr_patch
      INTEGER,              INTENT(IN)    :: jg
      REAL(wp),             INTENT(IN)    :: lon,lat
      ! local variables
      REAL(wp), PARAMETER              :: ZERO_TOL = 1.e-5_wp
      TYPE(t_geographical_coordinates) :: test_coord
      INTEGER                          :: test_jc, test_jb, i_startidx, i_endidx, &
        &                                 nblks_lonlat, npromz_lonlat
      REAL(wp)                         :: pi_180, query_lon, query_lat

      pi_180 = ATAN(1._wp)/45._wp
      query_lon = lon*pi_180
      query_lat = lat*pi_180

      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)
      DO test_jb = 1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (test_jb == nblks_lonlat) i_endidx = npromz_lonlat
        
        DO test_jc=i_startidx,i_endidx
          test_coord = ptr_int_lonlat%ll_coord(test_jc,test_jb)
          IF ((ABS(test_coord%lon - query_lon) < ZERO_TOL) .AND.  &
            & (ABS(test_coord%lat - query_lat) < ZERO_TOL)) THEN
            CALL ptr_int_lonlat%baryctr%visualize_stencil(&
              &               test_coord, test_jc, test_jb, ptr_patch, &
              &               "test_baryctr_"//TRIM(int2string(jg,"(i0)"))//".asy")
            CALL ptr_int_lonlat%rbf_c2l%visualize_stencil(&
              &               test_coord, test_jc, test_jb, ptr_patch, &
              &               "test_rbf_"//TRIM(int2string(jg,"(i0)"))//".asy")
          END IF
        END DO
      END DO
    END SUBROUTINE visualize_stencil


    !>
    !! This routine creates the index list needed for RBF reconstruction
    !! of scalar values at cell centers to lon-lat points.
    !!
    SUBROUTINE rbf_c2l_index( ptr_patch, tri_idx, ptr_int, ptr_int_lonlat )

      ! patch on which computation is performed
      TYPE(t_patch),        INTENT(INOUT)  :: ptr_patch
      INTEGER,              INTENT(IN)     :: tri_idx(:,:,:) ! 2, nproma, nblks_lonlat
      TYPE (t_int_state),   INTENT(IN)     :: ptr_int
      TYPE(t_lon_lat_intp), INTENT(inout)  :: ptr_int_lonlat

      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_c2l_index"
      INTEGER :: ilc_n(3), ibc_n(3)       ! line and block index for neighbors of
      ! direct neighbors
      INTEGER :: ilv(3), ibv(3)           ! vertex line and block indices
      INTEGER :: ilc_v(3,6), ibc_v(3,6)   ! cell line and block indices
      ! around each of the three vertices
      INTEGER :: jb                       ! loop index blocks
      INTEGER :: jc                       ! loop index cells
      INTEGER :: jj                       ! loop index
      INTEGER :: jec                      ! loop index
      INTEGER :: jtri                     ! loop index
      INTEGER :: jb_lonlat,jc_lonlat
      INTEGER :: cnt                      ! counter
      INTEGER :: i_startidx               ! start index
      INTEGER :: i_endidx                 ! end index
      INTEGER :: nblks_lonlat, npromz_lonlat

      !--------------------------------------------------------------------

      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

      ! stencil size: 4  = nearest neighbor,
      !               13 = vertex stencil,
      !               rbf_c2grad_dim = edge stencil

      IF ((rbf_dim_c2l == 4) .OR. (rbf_dim_c2l == 13)) THEN

        ptr_int_lonlat%rbf_c2l%idx(:,:,:) = 0

!$OMP PARALLEL

        ! The stencil consists of 12 cells surrounding the control
        ! volume that share a vertex with the control volume, plus the
        ! control volume itself.
        !
        ! Note: At pentagon points the size of the stencil reduces to 12.

!$OMP DO PRIVATE(jb,jc,jec,jj,jtri,i_startidx,i_endidx,cnt,ilv,ibv, &
!$OMP            ilc_v,ibc_v,ilc_n,ibc_n,jb_lonlat,jc_lonlat)
        DO jb_lonlat = 1,nblks_lonlat
          
          i_startidx = 1
          i_endidx   = nproma
          IF (jb_lonlat == nblks_lonlat) i_endidx = npromz_lonlat
          
          DO jc_lonlat = i_startidx, i_endidx
            jc = tri_idx(1,jc_lonlat, jb_lonlat)
            jb = tri_idx(2,jc_lonlat, jb_lonlat)
          
            cnt = 1

            ! get get line and block indices of cell vertices
            ilv(1:3) = ptr_patch%cells%vertex_idx(jc,jb,1:3)
            ibv(1:3) = ptr_patch%cells%vertex_blk(jc,jb,1:3)

            ! for each vertex: get all the cells which share this vertex
            DO jj = 1,3
              ilc_v(jj,:)=ptr_patch%verts%cell_idx(ilv(jj),ibv(jj),:)
              ibc_v(jj,:)=ptr_patch%verts%cell_blk(ilv(jj),ibv(jj),:)
            ENDDO

            ! 1. add the 3 direct neighbors to the stencil
            DO jec = 1, 3
              ! get line and block indices of direct neighbors
              ilc_n(jec) = ptr_patch%cells%neighbor_idx(jc,jb,jec)
              ibc_n(jec) = ptr_patch%cells%neighbor_blk(jc,jb,jec)

              ptr_int_lonlat%rbf_c2l%idx(cnt,jc_lonlat,jb_lonlat) = ilc_n(jec)
              ptr_int_lonlat%rbf_c2l%blk(cnt,jc_lonlat,jb_lonlat) = ibc_n(jec)

              cnt = cnt + 1
            ENDDO

            ! 2. loop over the vertices and add all the cells
            !    that are no direct neighbors and not our CV.

            IF (rbf_dim_c2l == 13) THEN
              DO jj = 1,3   ! loop over vertices
                DO jtri=1,6 ! loop over cells around each vertex

                  IF (.NOT.( (ilc_v(jj,jtri) == ilc_n(1) .AND. ibc_v(jj,jtri) == ibc_n(1))  &
                    &  .OR.  (ilc_v(jj,jtri) == ilc_n(2) .AND. ibc_v(jj,jtri) == ibc_n(2))  &
                    &  .OR.  (ilc_v(jj,jtri) == ilc_n(3) .AND. ibc_v(jj,jtri) == ibc_n(3))  &
                    &  .OR.  (ilc_v(jj,jtri) == jc       .AND. ibc_v(jj,jtri) == jb)        &
                    &  .OR.  (ilc_v(jj,jtri) == 0        .AND. ibc_v(jj,jtri) == 0 ) ) ) THEN

                    ptr_int_lonlat%rbf_c2l%idx(cnt,jc_lonlat,jb_lonlat) = ilc_v(jj,jtri)
                    ptr_int_lonlat%rbf_c2l%blk(cnt,jc_lonlat,jb_lonlat) = ibc_v(jj,jtri)

                    cnt = cnt + 1
                  ENDIF
                ENDDO
              ENDDO
            END IF

            ! 3. Finally, add the control volume itself
            ptr_int_lonlat%rbf_c2l%idx(cnt,jc_lonlat,jb_lonlat) = jc
            ptr_int_lonlat%rbf_c2l%blk(cnt,jc_lonlat,jb_lonlat) = jb

            ptr_int_lonlat%rbf_c2l%stencil(jc_lonlat,jb_lonlat) = cnt

          ENDDO ! loop over cells
        ENDDO ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ELSE IF (rbf_dim_c2l == rbf_c2grad_dim) THEN

        DO jb_lonlat = 1,nblks_lonlat
          
          i_startidx = 1
          i_endidx   = nproma
          IF (jb_lonlat == nblks_lonlat) i_endidx = npromz_lonlat
          
          DO jc_lonlat = i_startidx, i_endidx
            jc = tri_idx(1,jc_lonlat, jb_lonlat)
            jb = tri_idx(2,jc_lonlat, jb_lonlat)

            ! copy stencil info from "ptr_int%rbf_c2grad_idx":
            ptr_int_lonlat%rbf_c2l%idx(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_c2grad_idx(:,jc,jb)
            ptr_int_lonlat%rbf_c2l%blk(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_c2grad_blk(:,jc,jb)
            ptr_int_lonlat%rbf_c2l%stencil(jc_lonlat,jb_lonlat) = rbf_dim_c2l
          END DO
        END DO

      ELSE

        CALL finish(routine, "Unknown stencil!")

     END IF

    END SUBROUTINE rbf_c2l_index


    !-------------------------------------------------------------------------
    !> Computes the RBF coefficients for lon-lat grid points.
    !-------------------------------------------------------------------------
    !
    ! This routine is based on mo_intp_rbf_coeffs::rbf_vec_compute_coeff_cell()
    !
    SUBROUTINE rbf_compute_coeff_vec( ptr_patch, ptr_int_lonlat, rbf_shape_param )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
      REAL(wp),                      INTENT(IN)    :: rbf_shape_param
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_compute_coeff_vec"
      REAL(wp)                         :: cc_e1(3), cc_e2(3), cc_c(nproma,3)  ! coordinates of edge midpoints
      TYPE(t_cartesian_coordinates)    :: cc_center                  ! cartes. coordinates of lon-lat points
      REAL(wp), DIMENSION (nproma,3)   :: z_nx1, z_nx2, z_nx3        ! 3d  normal velocity vectors at edge midpoints
      INTEGER                          :: nblks_lonlat, npromz_lonlat ! blocking info
      REAL(wp)                         :: z_lon, z_lat,            & ! longitude and latitude
        &                                 z_norm,                  & ! norm of velocity vectors
        &                                 z_dist,                  & ! distance between data points
        &                                 z_nxprod                   ! scalar prod of normal veloc.
      REAL(wp),ALLOCATABLE             :: z_rbfmat(:,:,:),         & ! RBF interpolation matrix
        &                                 z_diag(:,:),             & ! diagonal of Cholesky decomp.
        &                                 z_rbfval(:,:),           & ! RBF function value
        &                                 z_rhs1(:,:),             & ! right hand side of linear
        &                                 z_rhs2(:,:)                ! interpolation system in 2d
      INTEGER                          :: &
        &                                 jc, jb,                  & ! integer over lon-lat grid points
        &                                 je1, je2,                & ! integer over edges
        &                                 ile1, ibe1, ile2, ibe2,  & ! edge indices
        &                                 ist,                     & ! return value of array allocation
        &                                 i_startidx, i_endidx,    & ! start/end index
        &                                 istencil(nproma)           ! actual number of stencil points
      REAL(wp)                         :: checksum_u,checksum_v      ! to check if sum of interpolation coefficients is correct
      TYPE(t_geographical_coordinates) :: grid_point

      !--------------------------------------------------------------------

      CALL message(routine, '')
      IF (.NOT. ptr_patch%domain_is_owned) RETURN;

      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist,   &
!$OMP                   grid_point),                                   &
!$OMP          SHARED  (nproma, rbf_vec_dim_c, nblks_lonlat,           &
!$OMP                   npromz_lonlat, ptr_int_lonlat, ptr_patch,      &
!$OMP                   rbf_vec_kern_ll)

      ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_c,rbf_vec_dim_c),  &
        z_diag(nproma,rbf_vec_dim_c),                  &
        z_rbfval(nproma,rbf_vec_dim_c),                &
        z_rhs1(nproma,rbf_vec_dim_c),                  &
        z_rhs2(nproma,rbf_vec_dim_c), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, 'allocation for working arrays failed')
      ENDIF

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil,       &
!$OMP             ist,ile1,ibe1,cc_e1,z_lon,z_lat,z_norm,           &
!$OMP             z_nx1,ile2,ibe2,cc_e2,cc_c,z_nx2,z_nxprod,z_dist, &
!$OMP             cc_center,z_nx3,checksum_u,checksum_v)
      BLOCKS: DO jb = 1, nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        if (jb == nblks_lonlat) i_endidx = npromz_lonlat

        ! for each cell, build the vector RBF interpolation matrix
        DO je1 = 1, rbf_vec_dim_c
          DO je2 = 1, je1
!$NEC ivdep
            DO jc = i_startidx, i_endidx

              ! Get actual number of stencil points
              istencil(jc) = ptr_int_lonlat%rbf_vec%stencil(jc,jb)

              ! paranoia:
              IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE

              ! line and block indices for each edge je1 and je2 of RBF stencil
              ile1 = ptr_int_lonlat%rbf_vec%idx(je1,jc,jb)
              ibe1 = ptr_int_lonlat%rbf_vec%blk(je1,jc,jb)
              ile2 = ptr_int_lonlat%rbf_vec%idx(je2,jc,jb)
              ibe2 = ptr_int_lonlat%rbf_vec%blk(je2,jc,jb)

              ! get Cartesian coordinates and orientation vectors
              cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
              cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
              !
              z_nx1(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)
              z_nx2(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

              ! compute dot product of normal vectors and distance between edge midpoints
              z_nxprod = DOT_PRODUCT(z_nx1(jc,:),z_nx2(jc,:))
              z_dist   = arc_length_v(cc_e1,cc_e2,ptr_patch%geometry_info)

              ! set up interpolation matrix
              IF      (rbf_vec_kern_ll == 1) THEN
                z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_shape_param)
              ELSE IF (rbf_vec_kern_ll == 3) THEN
                z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_shape_param)
              ENDIF

              IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)

            END DO
          END DO
        END DO

        ! apply Cholesky decomposition to matrix
        !
        CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
!$NEC ivdep
        DO jc = i_startidx, i_endidx

          !
          ! Solve immediately for coefficients
          !
          ! convert coordinates to cartesian vector
          !
          grid_point = ptr_int_lonlat%ll_coord(jc,jb)
          cc_center = gc2cc(grid_point,ptr_patch%geometry_info)
          cc_c(jc,1:3) = cc_center%x(1:3)

          z_lon = grid_point%lon
          z_lat = grid_point%lat

          ! Zonal wind component
          CALL gvec2cvec(1._wp,0._wp, z_lon,z_lat, &
            &            z_nx1(jc,1),z_nx1(jc,2),z_nx1(jc,3),ptr_patch%geometry_info)

          z_norm = SQRT( DOT_PRODUCT(z_nx1(jc,:),z_nx1(jc,:)) )
          z_nx1(jc,:)  = 1._wp/z_norm * z_nx1(jc,:)

          ! Meridional wind component
          CALL gvec2cvec(0._wp,1._wp, z_lon,z_lat, &
            &            z_nx2(jc,1),z_nx2(jc,2),z_nx2(jc,3),ptr_patch%geometry_info)

          z_norm = SQRT( DOT_PRODUCT(z_nx2(jc,:),z_nx2(jc,:)) )
          z_nx2(jc,:)  = 1._wp/z_norm * z_nx2(jc,:)

        END DO

        !
        ! set up right hand side for interpolation system
        !
        DO je2 = 1, rbf_vec_dim_c
!$NEC ivdep
          DO jc = i_startidx, i_endidx

            IF (je2 > istencil(jc)) CYCLE

            ! get indices and coordinates of edge midpoints and compute distance
            ! to cell center
            ile2   = ptr_int_lonlat%rbf_vec%idx(je2,jc,jb)
            ibe2   = ptr_int_lonlat%rbf_vec%blk(je2,jc,jb)

            cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)

            z_dist = arc_length_v(cc_c(jc,:), cc_e2,ptr_patch%geometry_info)

            ! get Cartesian orientation vector
            z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

            IF (rbf_vec_kern_ll == 1) THEN
              z_rbfval(jc,je2) = gaussi(z_dist,rbf_shape_param)
            ELSE IF (rbf_vec_kern_ll == 3) THEN
              z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_shape_param)
            ENDIF

            ! compute projection on target vector orientation
            z_rhs1(jc,je2) = z_rbfval(jc,je2) * &
              &                     DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
            z_rhs2(jc,je2) = z_rbfval(jc,je2) * &
              &                     DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))

          END DO
        END DO

        ! compute vector coefficients
        CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
          &               z_diag, z_rhs1, ptr_int_lonlat%rbf_vec%coeff(:,1,:,jb))
        CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
          &               z_diag, z_rhs2, ptr_int_lonlat%rbf_vec%coeff(:,2,:,jb))

        DO jc = i_startidx, i_endidx

          ! Ensure that sum of interpolation coefficients is correct
          checksum_u = 0._wp
          checksum_v = 0._wp

          DO je1 = 1, istencil(jc)
            ile1   = ptr_int_lonlat%rbf_vec%idx(je1,jc,jb)
            ibe1   = ptr_int_lonlat%rbf_vec%blk(je1,jc,jb)

            ! get Cartesian orientation vector
            z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)

            checksum_u = checksum_u + ptr_int_lonlat%rbf_vec%coeff(je1,1,jc,jb)* &
              DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
            checksum_v = checksum_v + ptr_int_lonlat%rbf_vec%coeff(je1,2,jc,jb)* &
              DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))
          ENDDO

          DO je1 = 1, istencil(jc)
            ptr_int_lonlat%rbf_vec%coeff(je1,1,jc,jb) = &
              ptr_int_lonlat%rbf_vec%coeff(je1,1,jc,jb) / checksum_u
          END DO
          DO je1 = 1, istencil(jc)
            ptr_int_lonlat%rbf_vec%coeff(je1,2,jc,jb) = &
              ptr_int_lonlat%rbf_vec%coeff(je1,2,jc,jb) / checksum_v
          END DO

        END DO ! jc

      END DO BLOCKS
!$OMP END DO

      DEALLOCATE( z_rbfmat, z_diag, z_rbfval, z_rhs1, z_rhs2,  STAT=ist )
      IF (ist /= SUCCESS)  CALL finish (routine, 'deallocation for working arrays failed')
!$OMP END PARALLEL

    END SUBROUTINE rbf_compute_coeff_vec


    !-------------------------------------------------------------------------
    !> Computes the RBF coefficients for lon-lat grid points.
    !
    !  We directly interpolate from cell centers to lon-lat points,
    !  otherwise we do gradient interpolation and reconstruction.
    !  -------------------------------------------------------------------------
    !
    SUBROUTINE rbf_compute_coeff_c2l( ptr_patch, ptr_int_lonlat, rbf_shape_param )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
      REAL(wp),                      INTENT(IN)    :: rbf_shape_param
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_compute_coeff_c2l"
      REAL(wp)                         :: cc_c(nproma,3)             ! coordinates of cell centers
      TYPE(t_cartesian_coordinates)    :: cc_center, cc_1, cc_2      ! temporary variables
      REAL(wp)                         :: z_dist                     ! distance between data points
      REAL(wp),ALLOCATABLE             :: z_rbfmat(:,:,:),         & ! RBF interpolation matrix
        &                                 z_diag(:,:),             & ! diagonal of Cholesky decomp.
        &                                 z_rbfval(:,:)              ! RBF function value
      INTEGER                          :: jc, jb,                  & ! integer over lon-lat grid points
        &                                 je1, je2,                & ! integer over edges
        &                                 ilc1, ibc1, ilc2, ibc2,  & ! cell indices
        &                                 ist,                     & ! return value of array allocation
        &                                 i_startidx, i_endidx,    & ! start/end index
        &                                 istencil(nproma)           ! actual number of stencil points
      REAL(wp)                         :: checksum                   ! to check if sum of interpolation coefficients is correct
      INTEGER                          :: nblks_lonlat, npromz_lonlat
      TYPE(t_geographical_coordinates) :: grid_point

      !--------------------------------------------------------------------

      CALL message(routine, '')
      IF (.NOT. ptr_patch%domain_is_owned) RETURN;

      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval, ist)

      ALLOCATE( z_rbfmat(nproma,rbf_dim_c2l,rbf_dim_c2l),    &
        &       z_diag(nproma,rbf_dim_c2l),                  &
        &       z_rbfval(nproma,rbf_dim_c2l),                &
        &       STAT=ist )
      IF (ist /= SUCCESS) &
        & CALL finish (routine, 'allocation for working arrays failed')

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil, grid_point, &
!$OMP             ist,ilc1,ibc1, ilc2,ibc2,cc_1, cc_2,cc_c,z_dist,        &
!$OMP             cc_center, checksum )
      BLOCKS: DO jb = 1, nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        if (jb == nblks_lonlat) i_endidx = npromz_lonlat

        ! for each cell, build the vector RBF interpolation matrix
        z_rbfmat(:,:,:) = 0._wp
        DO je1 = 1, rbf_dim_c2l
          DO je2 = 1, je1
!$NEC ivdep
            DO jc = i_startidx, i_endidx

              ! Get actual number of stencil points
              istencil(jc) = ptr_int_lonlat%rbf_c2l%stencil(jc,jb)
              ! paranoia:
              IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE

              ! line and block indices for each cell of RBF stencil
              ilc1 = ptr_int_lonlat%rbf_c2l%idx(je1,jc,jb)
              ibc1 = ptr_int_lonlat%rbf_c2l%blk(je1,jc,jb)
              ilc2 = ptr_int_lonlat%rbf_c2l%idx(je2,jc,jb)
              ibc2 = ptr_int_lonlat%rbf_c2l%blk(je2,jc,jb)

              ! get Cartesian coordinates and orientation vectors
              cc_1 = gc2cc(ptr_patch%cells%center(ilc1,ibc1),ptr_patch%geometry_info)
              cc_2 = gc2cc(ptr_patch%cells%center(ilc2,ibc2),ptr_patch%geometry_info)
              ! compute distance between cell centers
              z_dist = arc_length_v(cc_1%x(:),cc_2%x(:),ptr_patch%geometry_info)

              ! set up interpolation matrix
              IF      (rbf_vec_kern_ll == 1) THEN
                z_rbfmat(jc,je1,je2) = gaussi(z_dist,rbf_shape_param)
              ELSE IF (rbf_vec_kern_ll == 3) THEN
                z_rbfmat(jc,je1,je2) = inv_multiq(z_dist,rbf_shape_param)
              ENDIF

              IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)

            END DO
          END DO
        END DO

        ! apply Cholesky decomposition to matrix
        !
        CALL choldec_v(i_startidx,i_endidx,istencil,            z_rbfmat,z_diag)

        ! compute RHS for coefficient computation
!$NEC ivdep
        DO jc = i_startidx, i_endidx

          grid_point = ptr_int_lonlat%ll_coord(jc,jb)
          ! convert coordinates to cartesian vector
          cc_center = gc2cc(grid_point,ptr_patch%geometry_info)
          cc_c(jc,1:3) = cc_center%x(1:3)

        END DO

        !
        ! set up right hand side for interpolation system
        !
        DO je2 = 1, rbf_dim_c2l
!$NEC ivdep
          DO jc = i_startidx, i_endidx

            IF (je2 > istencil(jc)) CYCLE

            ! get indices and coordinates of cell centers and compute
            ! distance to lon-lat point:
            ilc2   = ptr_int_lonlat%rbf_c2l%idx(je2,jc,jb)
            ibc2   = ptr_int_lonlat%rbf_c2l%blk(je2,jc,jb)
            cc_1 = gc2cc(ptr_patch%cells%center(ilc2,ibc2),ptr_patch%geometry_info)

            z_dist = arc_length_v(cc_c(jc,:), cc_1%x(:),ptr_patch%geometry_info)

            IF (rbf_vec_kern_ll == 1) THEN
              z_rbfval(jc,je2) = gaussi(z_dist,rbf_shape_param)
            ELSE IF (rbf_vec_kern_ll == 3) THEN
              z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_shape_param)
            ENDIF

          END DO
        END DO

        ! compute vector coefficients
        CALL solve_chol_v(i_startidx, i_endidx, istencil,              z_rbfmat,  &
          &               z_diag, z_rbfval, ptr_int_lonlat%rbf_c2l%coeff(:,:,jb))

        DO jc = i_startidx, i_endidx

          ! Ensure that sum of interpolation coefficients is correct
          checksum = 0._wp
          DO je1 = 1, istencil(jc)
            checksum = checksum + ptr_int_lonlat%rbf_c2l%coeff(je1,jc,jb)
          ENDDO
          DO je1 = 1, istencil(jc)
            ptr_int_lonlat%rbf_c2l%coeff(je1,jc,jb) = &
              ptr_int_lonlat%rbf_c2l%coeff(je1,jc,jb) / checksum
          END DO
        END DO ! jc

      END DO BLOCKS
!$OMP END DO

      DEALLOCATE( z_rbfmat, z_diag, z_rbfval, STAT=ist )
      IF (ist /= SUCCESS)  CALL finish (routine, 'deallocation for working arrays failed')
!$OMP END PARALLEL

    END SUBROUTINE rbf_compute_coeff_c2l


    !-------------------------------------------------------------------------
    !> Mask out nest boundary zone.
    !  -------------------------------------------------------------------------
    !
    SUBROUTINE mask_out_boundary( ptr_patch, ptr_int_lonlat )

      ! Input parameters
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      ! Local parameters
      CHARACTER(*), PARAMETER :: routine      = modname//"::mask_out_boundary"
      INTEGER,      PARAMETER :: intp_types(3) = (/ HINTP_TYPE_LONLAT_RBF, &
        &                                           HINTP_TYPE_LONLAT_NNB, &
        &                                           HINTP_TYPE_LONLAT_BCTR /)

      INTEGER  :: last_bdry_index, rl_start, rl_end, i_nchdom, i_startblk, i_endblk, &
        &         i_startidx, i_endidx, local_idx, glb_idx, jc, jb, je1,             &
        &         istencil, istencil2, itype, hintp_type, nblks_lonlat,  &
        &         npromz_lonlat, jc_src, jb_src, jc_nb, jb_nb
      REAL(wp) :: area
      TYPE(t_intp_scalar_coeff), POINTER :: ptr_intp

      !--------------------------------------------------------------------

      IF (dbg_level > 1)  CALL message(routine,'')

      IF (.NOT. ptr_patch%domain_is_owned) RETURN

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

      ! communicate the largest global index of the lateral boundary cells.
      ! 
      ! note that this implicitly assumes that the indices of the nest
      ! boundary region are stored at the beginning of the solution
      ! vector in a contiguous fashion.
      last_bdry_index = 0

      rl_start   = 1
      rl_end     = grf_bdywidth_c
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      glb_idx    = 0
      DO jb=i_startblk,i_endblk
        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)
        DO jc=i_startidx, i_endidx
          local_idx = idx_1d(jc,jb)
          glb_idx   = MAX(glb_idx, ptr_patch%cells%decomp_info%glb_index(local_idx))
        END DO
      END DO
      last_bdry_index = p_max(glb_idx, p_comm_work)

      DO itype=1,SIZE(intp_types)
        hintp_type = intp_types(itype)
        IF (dbg_level > 1)  CALL message(routine, "interpolation type: "//int2string(hintp_type))

        ! depending on the interpolation type: set a pointer to the
        ! coefficient vector
        NULLIFY(ptr_intp)
        SELECT CASE(hintp_type)
        CASE (HINTP_TYPE_LONLAT_RBF)
          ptr_intp  => ptr_int_lonlat%rbf_c2l
        CASE (HINTP_TYPE_LONLAT_NNB)
          ptr_intp  => ptr_int_lonlat%nnb
        CASE (HINTP_TYPE_LONLAT_BCTR)
          IF (.NOT. support_baryctr_intp)  CYCLE

          ! Tentative fix: "lreduced_nestbdry_stencil" controls
          ! only barycentric interpolation.
          IF (.NOT. lreduced_nestbdry_stencil) CYCLE

          ptr_intp  => ptr_int_lonlat%baryctr
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT

        ! Treat the masking differently for nearest-neighbor
        ! interpolation (where masking would result in a missing
        ! value). The general case is handled below.
        !
        IF (hintp_type == HINTP_TYPE_LONLAT_NNB) THEN

!$OMP PARALLEL DO PRIVATE(jb, i_startidx, i_endidx, jc, je1, local_idx, glb_idx, &
!$OMP                     jc_src, jb_src, jc_nb, jb_nb)
          DO jb = 1, nblks_lonlat

            i_startidx = 1
            i_endidx   = nproma
            IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

            ! loop over all lonlat points and their input stencils
            CELL_LOOP : DO jc = i_startidx, i_endidx

              ! if stencil point is from the nest boundary region,
              ! then replace it by one of its neighbors.
              jc_src = ptr_intp%idx(1,jc,jb)
              jb_src = ptr_intp%blk(1,jc,jb)
              local_idx = idx_1d(jc_src, jb_src)
              glb_idx   = ptr_patch%cells%decomp_info%glb_index(local_idx)
              IF (glb_idx <= last_bdry_index) THEN
                DO je1 = 1,3
                  ! get line and block indices of direct neighbors
                  jc_nb   = ptr_patch%cells%neighbor_idx(jc_src,jb_src,je1)
                  jb_nb   = ptr_patch%cells%neighbor_blk(jc_src,jb_src,je1)
                  local_idx = idx_1d(jc_nb, jb_nb)
                  glb_idx   = ptr_patch%cells%decomp_info%glb_index(local_idx)
                  IF (glb_idx > last_bdry_index) THEN
                    ptr_intp%idx(1,jc,jb) = jc_nb
                    ptr_intp%blk(1,jc,jb) = jb_nb
                    CYCLE CELL_LOOP
                  END IF
                END DO
              END IF
            END DO CELL_LOOP ! jc

          END DO
!$OMP END PARALLEL DO

        ELSE

!$OMP PARALLEL DO PRIVATE(jb, i_startidx, i_endidx, jc, istencil, &
!$OMP                     istencil2, area, je1, local_idx, glb_idx)
          BLOCKS: DO jb = 1, nblks_lonlat

            i_startidx = 1
            i_endidx   = nproma
            IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

            ! loop over all lonlat points and their input stencils
            DO jc = i_startidx, i_endidx

              ! Get actual number of stencil points
              istencil  = ptr_intp%stencil(jc,jb)
              istencil2 = istencil
              ! first, multiply interpolation weights by total "area" (the
              ! sum of weights in the stencil). we perform this step
              ! though we may assume normalized coefficients (area=1).
              area = 0._wp
              DO je1 = 1, istencil
                area = area + ptr_intp%coeff(je1,jc,jb)
              ENDDO
              DO je1 = 1, istencil
                ptr_intp%coeff(je1,jc,jb) = ptr_intp%coeff(je1,jc,jb) * area
              END DO
              ! remove stencil points from the nest boundary region
              ! together with their coefficients
              DO je1 = 1, istencil
                local_idx = idx_1d(ptr_intp%idx(je1,jc,jb), ptr_intp%blk(je1,jc,jb))
                IF (local_idx < 1) THEN
                  WRITE (0,*) "interpolation type: hintp_type=", hintp_type
                  WRITE (0,*) "PE ", get_my_mpi_work_id(), ": je1, jc,jb = ", je1, jc,jb, &
                    & ", indices: ", ptr_intp%idx(je1,jc,jb), ptr_intp%blk(je1,jc,jb)
                  CALL finish(routine, "Internal error!")
                END IF
                glb_idx = ptr_patch%cells%decomp_info%glb_index(local_idx)
                IF (glb_idx <= last_bdry_index) THEN
                  area = area - ptr_intp%coeff(je1,jc,jb)
                  ptr_intp%coeff(je1,jc,jb) = 0._wp
                  istencil2 = istencil2 - 1
                  IF (istencil2 == 0)  area = 1._wp ! avoid zero divide
                END IF
              END DO
              ! re-normalize weights
              DO je1 = 1, istencil
                ptr_intp%coeff(je1,jc,jb) = ptr_intp%coeff(je1,jc,jb) / area
              END DO
            END DO ! jc

          END DO BLOCKS
!$OMP END PARALLEL DO

        END IF

      END DO  
      IF (dbg_level > 1)  CALL message(routine, "done.")
      
    END SUBROUTINE mask_out_boundary


    ! -----------------------------------------------------------------------------------
    !> Utility function: distance computation
    !
    PURE FUNCTION dist_p(p1, p2)
      REAL(gk)              :: dist_p
      REAL(gk), INTENT(IN)  :: p1(2), p2(2)
      ! local variables
      REAL(gk) :: val
      
      ! spherical distance:
      val = SIN(p1(2))*SIN(p2(2)) + COS(p1(2))*COS(p2(2))*COS(p1(1)-p2(1))
      dist_p = ACOS( MIN(1._gk, MAX(-1._gk, val)) )
    END FUNCTION dist_p


    ! -----------------------------------------------------------------------------------
    !> Utility function: Divide-and-Conquer algorithm to determine
    !  points in lon-lat grid which are "far off" the current patch of
    !  the triangular grid.
    !
    RECURSIVE SUBROUTINE flag_ll_points(rotated_pts, s_lon, e_lon, s_lat, e_lat, &
      &                                 pts_flags, recursion_depth, max_recursion, gnat, min_radius)
      REAL(wp), INTENT(IN)    :: rotated_pts(:,:,:)  ! dim (lon,lat,1:2)
      INTEGER,  INTENT(IN)    :: s_lon, e_lon, s_lat, e_lat
      INTEGER,  INTENT(INOUT) :: pts_flags(:,:)
      INTEGER,  INTENT(IN)    :: recursion_depth, max_recursion
      TYPE(t_gnat_tree), INTENT(IN) :: gnat
      REAL(wp), INTENT(IN)    :: min_radius
      ! local variables
      INTEGER  :: m_lon, m_lat, c_lon(4), c_lat(4), v_lon(4,4), v_lat(4,4), &
        &         s_lon2(4), s_lat2(4), e_lon2(4), e_lat2(4), icirc, ivert, &
        &         i1, i2
      REAL(gk) :: p1(4,2), p2(2), radius(4)
      LOGICAL  :: is_pt_inside(4)

      IF (recursion_depth == max_recursion) RETURN
      if ((s_lon == e_lon) .or. (s_lat == e_lat)) return

      ! determine mid longitude and latitude in current range, thus
      ! dividing the regular grid into 4 subgrids:
      m_lon = (e_lon + s_lon)/2
      m_lat = (e_lat + s_lat)/2

      ! nomenclature
      !
      ! *-------*-------*   - e_lat
      ! |       |       |
      ! |   3   |   4   |
      ! |       |       |
      ! *-------*-------*   - m_lat
      ! |       |       |
      ! |   1   |   2   |
      ! |       |       |
      ! *-------*-------*   - s_lat
      !
      ! |       |       |
      !  s_lon   m_lon   e_lon

      s_lon2(1:4) = (/ s_lon, m_lon, s_lon, m_lon /)
      e_lon2(1:4) = (/ m_lon, e_lon, m_lon, e_lon /)
      s_lat2(1:4) = (/ s_lat, s_lat, m_lat, m_lat /)
      e_lat2(1:4) = (/ m_lat, m_lat, e_lat, e_lat /)
      ! for each of the four subblocks, determine centre point and
      ! vertex indices:
      DO icirc=1,4
        c_lon(icirc) = (s_lon2(icirc) + e_lon2(icirc)) / 2
        c_lat(icirc) = (s_lat2(icirc) + e_lat2(icirc)) / 2
        ! vertex numbering counter-clockwise
        v_lon(icirc,1:4) = (/ s_lon2(icirc), e_lon2(icirc), e_lon2(icirc), s_lon2(icirc) /)
        v_lat(icirc,1:4) = (/ s_lat2(icirc), s_lat2(icirc), e_lat2(icirc), e_lat2(icirc) /)
      END DO
      ! determine radii of circles overlapping each subblock:
      radius(:) = 0._gk
      DO icirc=1,4
        p1(icirc,1:2) = rotated_pts(c_lon(icirc),c_lat(icirc),1:2)
        DO ivert=1,4
          p2(1:2) = rotated_pts(v_lon(icirc,ivert),v_lat(icirc,ivert),1:2)
          radius(icirc) = MAX(radius(icirc), dist_p(p1(icirc,1:2), p2))
        END DO
      END DO
      radius(1:4) = MAX(min_radius, 1.25_gk * radius(1:4)) ! safety margin
      ! test, if there exists a triangle center inside the circles:
      DO icirc=1,4
        is_pt_inside(icirc) = .FALSE.
        CALL gnat_recursive_proximity_query(gnat, gnat%gnat_tree, p1(icirc,1:2), radius(icirc), is_pt_inside(icirc))
      END DO
      ! for each of the four subblocks:
      DO icirc=1,4
        IF (.NOT. is_pt_inside(icirc)) THEN
          ! if patch of triangular grid does not intersect circle:
          ! mark corresponding points in regular grid with "SKIP_NODE"
          DO i2=s_lat2(icirc),e_lat2(icirc)
            DO i1=s_lon2(icirc),e_lon2(icirc)
              pts_flags(i1,i2) = SKIP_NODE
            END DO
          END DO
        ELSE
          ! if patch of triangular grid does intersect circle:
          ! recursion and further subdivision:
          CALL flag_ll_points(rotated_pts, s_lon2(icirc), e_lon2(icirc), s_lat2(icirc), e_lat2(icirc), &
            &                 pts_flags, (recursion_depth+1), max_recursion, gnat, min_radius)          
        END IF
      END DO
    END SUBROUTINE flag_ll_points


    !-------------------------------------------------------------------------
    !> Setup routine for lon-lat grid points for an arbitrary grid.
    !
    ! @note   This subroutine assumes that the GNAT data structure for
    !         this patch has already been initialized.
    !
    ! Note:   Two variables are created globally: "tri_idx" and "global_idx"
    !         This is necessary for the splitting of the lon-lat grid
    !         points over the PEs and for creating the corresponding
    !         communication patterns. However, after this initial phase,
    !         these data structures are resized to the local sizes.
    !
    SUBROUTINE distribute_lonlat_points(grid, gnat, ptr_patch, tri_idx, ptr_int_lonlat)

      TYPE (t_lon_lat_grid), INTENT(INOUT)         :: grid
      TYPE (t_gnat_tree),    INTENT(INOUT)         :: gnat
      ! data structure containing grid info:
      TYPE(t_patch),         INTENT(INOUT)         :: ptr_patch
      INTEGER,               INTENT(INOUT)         :: tri_idx(:,:,:) ! 2, nproma, nblks_lonlat
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), INTENT(INOUT)         :: ptr_int_lonlat
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::distribute_lonlat_points"
      ! Flag: .TRUE., if we want raw performance measurements
      LOGICAL,      PARAMETER :: l_measure_time         = .FALSE.
      ! Zero threshold constant for equality test
      REAL(wp),     PARAMETER :: ZERO_TOL               = 1.e-15_wp

      REAL(wp), ALLOCATABLE            :: rotated_pts(:,:,:)
      REAL(gk), ALLOCATABLE            :: in_points(:,:,:)
      INTEGER,  ALLOCATABLE            :: pts_flags(:,:)
      REAL(gk), ALLOCATABLE            :: min_dist(:,:)       ! minimal distance
      INTEGER                          :: jb, jc, i_startidx, i_endidx,           &
        &                                 nblks_lonlat, npromz_lonlat, errstat,   &
        &                                 rl_start, rl_end, i_startblk,           &
        &                                 i, j, row_idx(2), ibeg_idx, ibeg_blk,   &
        &                                 iend_idx, iend_blk, tri_idx_pole,       &
        &                                 ibeg_glb, ist
      REAL(gk)                         :: min_dist_pole
      REAL(wp)                         :: point(2),                               &
        &                                 max_dist, start_radius
      LOGICAL                          :: l_grid_is_unrotated, l_grid_contains_poles, &
           l_my_process_is_mpi_test


      !-----------------------------------------------------------------------

      IF (dbg_level > 1) THEN
        WRITE(message_text,*) "SETUP : interpol_lonlat"
        CALL message(routine, message_text)
      END IF
      IF (ptr_patch%geometry_info%cell_type == 6) &
        &   CALL finish(routine, "Lon-lat interpolation not yet implemented for cell_type == 6!")

      nblks_lonlat  = grid%nblks
      npromz_lonlat = grid%npromz

      ! ------------------------------------------
      ! check, if we have a global, unrotated grid
      ! ------------------------------------------

      l_grid_is_unrotated   = ((ABS(90._wp - grid%north_pole(2)) < ZERO_TOL)       .AND.  &
        &                       ABS( 0._wp - grid%north_pole(1)) < ZERO_TOL)
      l_grid_contains_poles = ((ABS(90._wp - ABS(grid%reg_lat_def(1))) < ZERO_TOL) .AND.  &
        &                      (ABS(90._wp - ABS(grid%reg_lat_def(3))) < ZERO_TOL))

      ! ---------------------------------------------------------
      ! preliminaries: make a sensible guess for maximum distance
      ! ---------------------------------------------------------

      ! make a sensible guess for maximum distance (just taking a
      ! "randomly chosen" edge length)
      rl_start   = 2
      rl_end     = min_rlcell_int
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      CALL get_indices_e(ptr_patch, i_startblk, i_startblk, i_startblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      max_dist = 3._gk * &
        & REAL(ptr_patch%edges%primal_edge_length(i_startidx,i_startblk)/grid_sphere_radius, gk)
      ! for MPI-independent behaviour: determine global max.
      max_dist = p_max(max_dist, comm=p_comm_work)
      IF (p_test_run) THEN
        l_my_process_is_mpi_test = my_process_is_mpi_test()
        IF (p_pe_work == p_io) THEN
          IF(.NOT. l_my_process_is_mpi_test) THEN
            ! Send to test PE
            CALL p_send(max_dist, process_mpi_all_test_id, 1)
          ELSE
            ! Receive result from parallel worker PEs
            CALL p_recv(max_dist, process_mpi_all_workroot_id, 1)
          END IF
        END IF
        IF(l_my_process_is_mpi_test) THEN
          CALL p_bcast(max_dist, p_io, comm=p_comm_work)
        END IF
      END IF

      ! -------------------------------------------
      ! compute grid points of rotated lon/lat grid
      ! -------------------------------------------

      ALLOCATE(in_points(nproma, nblks_lonlat, 2),            &
        &      rotated_pts(grid%lon_dim, grid%lat_dim, 2),    &
        &      pts_flags(grid%lon_dim, grid%lat_dim),         &
        &      min_dist(nproma, nblks_lonlat),                &
        &      stat=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'allocation for working arrays failed')

      IF (dbg_level > 1) CALL message(routine, "rotate lon-lat grid points")
      ! compute grid points (in radians): "rotated_pts(lon,lat,1:2)"
      ! contains the longitude and latitude of the point with grid
      ! index lon/lat.
      CALL rotate_latlon_grid( grid, rotated_pts )
      in_points(:,:,1) = RESHAPE(REAL(rotated_pts(:,:,1), gk), (/ nproma, nblks_lonlat /), (/ 0._gk /) )
      in_points(:,:,2) = RESHAPE(REAL(rotated_pts(:,:,2), gk), (/ nproma, nblks_lonlat /), (/ 0._gk /) )

      ! --------------------------------------------------------------
      ! discard lon/lat grid points which are "far off"
      ! --------------------------------------------------------------

      pts_flags(:,:) = INVALID_NODE
      CALL flag_ll_points(rotated_pts, 1,grid%lon_dim, 1,grid%lat_dim, pts_flags, 1,5, gnat, max_dist)

      !-- if we have a global, unrotated grid: skip all points of the
      !   first and last latitude row (except one):
      IF (l_grid_is_unrotated .AND. l_grid_contains_poles) THEN
        pts_flags(2:,           1) = SKIP_NODE
        pts_flags(2:,grid%lat_dim) = SKIP_NODE
      END IF

      ! --------------------------------------------------------------
      ! proximity query, build a list of cell indices that contain the
      ! lon/lat grid points
      ! --------------------------------------------------------------

      tri_idx(1,:,:) = RESHAPE(pts_flags(:,:), (/ nproma, nblks_lonlat /), (/ INVALID_NODE /) )
      tri_idx(2,:,:) = 0


      ! Perform query. Note that for distributed patches we receive a
      ! local list of "in_points" actually located on this portion of the
      ! domain.
      IF (dbg_level > 1) CALL message(routine, "proximity query")
      start_radius = grid_sphere_radius
      ! work-around: increased search radius on torus
      if (is_plane_torus)  start_radius = 10._wp
      CALL gnat_query_containing_triangles(gnat, ptr_patch, in_points(:,:,:),      &
        &                 nproma, nblks_lonlat, npromz_lonlat, start_radius,       &
        &                 p_test_run, tri_idx(:,:,:), min_dist(:,:))


      !-- if we have a global, unrotated grid: copy the distance
      !   result for all points of the first and last latitude row:
      IF (l_grid_is_unrotated .AND. l_grid_contains_poles) THEN
        row_idx = (/ 1, grid%lat_dim /)
        DO i=1,2
          ! compute indices for pole row
          ibeg_glb = (row_idx(i)-1)*grid%lon_dim
          ibeg_idx = idx_no(ibeg_glb + 1)
          ibeg_blk = blk_no(ibeg_glb + 1)
          iend_idx = idx_no(ibeg_glb + grid%lon_dim)
          iend_blk = blk_no(ibeg_glb + grid%lon_dim)
          ! copy indices and distances:
          DO j=1,2
            tri_idx_pole = tri_idx(j,ibeg_idx,ibeg_blk)
            IF (iend_blk > ibeg_blk) THEN
              tri_idx(j, (ibeg_idx+1):nproma,                  ibeg_blk) = tri_idx_pole
              tri_idx(j,                   :, (ibeg_blk+1):(iend_blk-1)) = tri_idx_pole
              tri_idx(j,          1:iend_idx,                  iend_blk) = tri_idx_pole
            ELSE
              tri_idx(j, (ibeg_idx+1):iend_idx, ibeg_blk) = tri_idx_pole
            END IF
          END DO
          min_dist_pole = min_dist(ibeg_idx,ibeg_blk)
          IF (iend_blk > ibeg_blk) THEN
            min_dist((ibeg_idx+1):nproma,                  ibeg_blk) = min_dist_pole
            min_dist(                  :, (ibeg_blk+1):(iend_blk-1)) = min_dist_pole
            min_dist(         1:iend_idx,                  iend_blk) = min_dist_pole
          ELSE
            min_dist((ibeg_idx+1):iend_idx, ibeg_blk) = min_dist_pole
          END IF
        END DO
      END IF

      CALL gnat_merge_distributed_queries(ptr_patch, grid%total_dim, nproma, grid%nblks, &
        &                 min_dist, tri_idx(:,:,:), in_points(:,:,:),                    &
        &                 ptr_int_lonlat%global_idx(:), ptr_int_lonlat%nthis_local_pts)


      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

      ! store the geographical coordinates for the lon-lat points
      !
      ALLOCATE(ptr_int_lonlat%ll_coord(nproma, nblks_lonlat), stat=ist)
      IF (ist /= SUCCESS)  CALL finish (routine, 'Allocation of array with lon-lat coordinates!')

!$OMP PARALLEL DO PRIVATE(jb,jc,i_startidx,i_endidx,point)
      DO jb = 1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc = i_startidx, i_endidx
          point(1:2) = REAL(in_points(jc,jb,1:2), wp)
          ptr_int_lonlat%ll_coord(jc,jb)%lon = point(1)
          ptr_int_lonlat%ll_coord(jc,jb)%lat = point(2)
        END DO
      END DO
!$OMP END PARALLEL DO

      ! clean up
      DEALLOCATE(in_points, pts_flags, rotated_pts, min_dist, stat=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE of working arrays failed')

    END SUBROUTINE distribute_lonlat_points


    !-------------------------------------------------------------------------
    !> Setup routine for RBF interpolation on lon-lat grid points.
    !
    SUBROUTINE nnb_setup_interpol_lonlat_grid(tri_idx, ptr_int_lonlat)
      INTEGER,               INTENT(IN)    :: tri_idx(:,:,:) ! 2, nproma, nblks_lonlat
      TYPE (t_lon_lat_intp), INTENT(INOUT) :: ptr_int_lonlat
      ! local variables
      CHARACTER(*), PARAMETER :: routine = modname//"::nnb_setup_interpol_lonlat_grid"
      INTEGER :: jc_lonlat, jb_lonlat, jc, jb, nblks_lonlat, &
        &        npromz_lonlat, i_startidx, i_endidx

      IF (dbg_level > 5) CALL message(routine, "Enter")

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

!$OMP PARALLEL DO PRIVATE(jb_lonlat,jc_lonlat,i_startidx,i_endidx,jc,jb)
      DO jb_lonlat = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb_lonlat == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc_lonlat = i_startidx, i_endidx

          jc = tri_idx(1,jc_lonlat, jb_lonlat)
          jb = tri_idx(2,jc_lonlat, jb_lonlat)

          ptr_int_lonlat%nnb%stencil(jc_lonlat,jb_lonlat) = 1
          ptr_int_lonlat%nnb%coeff(1,jc_lonlat,jb_lonlat) = 1._wp
          ptr_int_lonlat%nnb%idx(1,jc_lonlat,jb_lonlat)   = jc
          ptr_int_lonlat%nnb%blk(1,jc_lonlat,jb_lonlat)   = jb

        END DO
      END DO
!$OMP END PARALLEL DO

    END SUBROUTINE nnb_setup_interpol_lonlat_grid


    !-------------------------------------------------------------------------
    !> Setup routine for RBF interpolation on lon-lat grid points.
    !
    SUBROUTINE rbf_setup_interpol_lonlat_grid(ptr_patch, ptr_int_lonlat, tri_idx, ptr_int)

      TYPE(t_patch),         INTENT(INOUT) :: ptr_patch       !< data structure containing grid info
      TYPE (t_lon_lat_intp), INTENT(INOUT) :: ptr_int_lonlat  !< Indices of source points and interpolation coefficients
      INTEGER,               INTENT(IN)    :: tri_idx(:,:,:)  !< 2, nproma, nblks_lonlat
      TYPE (t_int_state),    INTENT(IN)    :: ptr_int
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::rbf_setup_interpol_lonlat_grid"

      INTEGER                          :: jb, jc, i_startidx, i_endidx,           &
        &                                 nblks_lonlat, npromz_lonlat,            &
        &                                 jb_lonlat, jc_lonlat
      REAL(wp)                         :: rbf_shape_param

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

      !-- edge indices of stencil
      ! are available through previous calls to SR rbf_vec_index_cell()

      ! -----------------------------------------------------
      ! copy neighbor indices from standard RBF interpolation
      ! -----------------------------------------------------

      IF (dbg_level > 1) CALL message(routine, "copy neighbor indices from standard RBF interpolation")

!$OMP PARALLEL SHARED(nblks_lonlat, nproma, npromz_lonlat,   &
!$OMP                 ptr_int_lonlat, ptr_int, tri_idx),     &
!$OMP          DEFAULT (NONE)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jc_lonlat,jb_lonlat), SCHEDULE(runtime)
      DO jb_lonlat = 1,nblks_lonlat

        i_startidx = 1
        i_endidx   = nproma
        IF (jb_lonlat == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc_lonlat = i_startidx, i_endidx

          jc = tri_idx(1,jc_lonlat, jb_lonlat)
          jb = tri_idx(2,jc_lonlat, jb_lonlat)

          ptr_int_lonlat%rbf_vec%stencil(jc_lonlat,jb_lonlat) = ptr_int%rbf_vec_stencil_c(jc,jb)
          ptr_int_lonlat%rbf_vec%idx(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_vec_idx_c(:,jc,jb)
          ptr_int_lonlat%rbf_vec%blk(:,jc_lonlat,jb_lonlat)   = ptr_int%rbf_vec_blk_c(:,jc,jb)
        ENDDO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      IF (dbg_level > 1) CALL message(routine, "compute lon-lat interpolation coefficients")

      ! -------------------------------------------------------------------------
      ! compute interpolation coefficients for RBF interpolation of vector fields
      ! -------------------------------------------------------------------------

      SELECT CASE (rbf_scale_mode_ll)
      CASE (SCALE_MODE_TABLE) 
        rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      CASE (SCALE_MODE_AUTO)
        ! if no shape parameter has been set: compute an estimate 
        CALL estimate_rbf_parameter(nblks_lonlat, npromz_lonlat, ptr_patch%edges%center,              &
          &                         ptr_int_lonlat%rbf_vec%idx, ptr_int_lonlat%rbf_vec%blk,           &
          &                         ptr_int_lonlat%rbf_vec%stencil, rbf_vec_dim_c,                    &
          &                         ptr_int_lonlat%global_idx, rbf_shape_param)
        rbf_shape_param          = p_min(rbf_shape_param, comm=p_comm_work)
        ptr_int_lonlat%rbf_scale = rbf_shape_param
        IF (my_process_is_stdio()) THEN
          WRITE(0,*) routine, ": auto-estimated shape_param = ", rbf_shape_param
        END IF
      CASE (SCALE_MODE_PRESET)
        IF (ptr_int_lonlat%rbf_scale <= 0._wp) THEN
          CALL finish(routine, "Explicitly presetting RBF shape parameter... invalid value!")
        END IF
        rbf_shape_param = ptr_int_lonlat%rbf_scale
      CASE DEFAULT
        CALL finish(routine, "Unknown value for rbf_scale_mode_ll!")
      END SELECT
      ! fallback solution:
      IF (rbf_shape_param <= 0._wp) &
        &  rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      WRITE(message_text,*) ": chosen rbf_shape_param = ", rbf_shape_param
      CALL message(routine, TRIM(message_text))

      CALL rbf_compute_coeff_vec( ptr_patch, ptr_int_lonlat, rbf_shape_param )

      ! -------------------------------------------------------------------------
      ! compute interpolation coefficients for RBF interpolation of scalar values
      ! -------------------------------------------------------------------------

      CALL rbf_c2l_index( ptr_patch, tri_idx, ptr_int, ptr_int_lonlat )

      SELECT CASE (rbf_scale_mode_ll)
      CASE (SCALE_MODE_TABLE) 
        rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      CASE (SCALE_MODE_AUTO)
        ! if no shape parameter has been set: compute an estimate 
        CALL estimate_rbf_parameter(nblks_lonlat, npromz_lonlat, ptr_patch%cells%center,       &
          &                         ptr_int_lonlat%rbf_c2l%idx, ptr_int_lonlat%rbf_c2l%blk,    &
          &                         ptr_int_lonlat%rbf_c2l%stencil, rbf_dim_c2l,               &
          &                         ptr_int_lonlat%global_idx, rbf_shape_param)
        rbf_shape_param = p_min(rbf_shape_param, comm=p_comm_work)
        ptr_int_lonlat%rbf_scale = rbf_shape_param
        IF (my_process_is_stdio()) THEN
          WRITE(0,*) routine, ": auto-estimated shape_param = ", rbf_shape_param
        END IF
      CASE (SCALE_MODE_PRESET)
        rbf_shape_param = ptr_int_lonlat%rbf_scale
      CASE DEFAULT
        CALL finish(routine, "Unknown value for rbf_scale_mode_ll!")
      END SELECT
      ! fallback solution:
      IF (rbf_shape_param <= 0._wp) &
        &  rbf_shape_param = rbf_vec_scale_ll(MAX(ptr_patch%id,1))
      WRITE(message_text,*) ": chosen rbf_shape_param = ", rbf_shape_param
      CALL message(routine, TRIM(message_text))
      ! compute coefficients:
      CALL rbf_compute_coeff_c2l( ptr_patch, ptr_int_lonlat, rbf_shape_param )

    END SUBROUTINE rbf_setup_interpol_lonlat_grid


    ! Interpolation weight algorithms may cause small weights outside
    ! [0,1]. This subroutine applies a cut-off and renormalization of
    ! cell-based interpolation coefficients.
    !
    SUBROUTINE normalize_weights(ptr_int_lonlat, intp_coeff)
      TYPE (t_lon_lat_intp),     INTENT(INOUT) :: ptr_int_lonlat
      TYPE(t_intp_scalar_coeff), INTENT(INOUT) :: intp_coeff
      ! local variables
      CHARACTER(*), PARAMETER :: routine = modname//"::normalize_weights"
      LOGICAL,      PARAMETER :: ldebug  = .FALSE.

      INTEGER  :: jc, jb, nblks_lonlat, npromz_lonlat, i_startidx, i_endidx, &
        &         je1, istencil, icount = 0
      REAL(wp) :: checksum, mmin, mmax, &
        &         rmin = 1.0_wp, rmax = 0.0_wp

      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

      icount = 0
      rmin   = 1.0_wp
      rmax   = 0.0_wp

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx, je1, checksum, istencil, mmin, mmax)  &
!$OMP    REDUCTION(+:icount), REDUCTION(min:rmin), REDUCTION(max:rmax)
      BLOCKS: DO jb = 1, nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        if (jb == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc = i_startidx, i_endidx
          istencil = intp_coeff%stencil(jc,jb) ! actual number of stencil points

          ! statistics
          mmin   = MINVAL(intp_coeff%coeff(1:istencil,jc,jb))
          rmin   = MIN(rmin, mmin)
          mmax   = MAXVAL(intp_coeff%coeff(1:istencil,jc,jb))
          rmax   = MAX(rmax, mmax)
          icount = icount + COUNT(intp_coeff%coeff(1:istencil,jc,jb) < 0._wp)
          icount = icount + COUNT(intp_coeff%coeff(1:istencil,jc,jb) > 1._wp)

          ! Apply cut-off, st. weights in [0,1]
          DO je1 = 1, istencil
            intp_coeff%coeff(je1,jc,jb) = MAX(MIN(1.0_wp, intp_coeff%coeff(je1,jc,jb)), 0.0_wp)
          END DO

          ! Ensure that sum of interpolation coefficients = 1.0
          checksum = 0._wp
          DO je1 = 1, istencil
            checksum = checksum + intp_coeff%coeff(je1,jc,jb)
          ENDDO
          DO je1 = 1, istencil
            intp_coeff%coeff(je1,jc,jb) = intp_coeff%coeff(je1,jc,jb) / checksum
          END DO
        END DO ! jc
      END DO BLOCKS
!$OMP END DO

      IF (ldebug .AND. (icount > 0)) THEN
        WRITE (0,*) routine, ": PE ", get_my_mpi_work_id(), &
          &         " handled ", icount, " cases; min/max = ", rmin, ", ", rmax
      END IF
    END SUBROUTINE normalize_weights

  END MODULE mo_intp_lonlat
