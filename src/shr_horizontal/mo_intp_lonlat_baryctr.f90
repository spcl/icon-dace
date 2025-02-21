! Contains the implementation of interpolation onto regular grids
! with barycentric interpolation.
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









  MODULE mo_intp_lonlat_baryctr
    !-------------------------------------------------------------------------
    !
    !    ProTeX FORTRAN source: Style 2
    !    modified for ICON project, DWD/MPI-M 2006
    !
    !-------------------------------------------------------------------------
    !



    USE mo_kind,                ONLY: wp
    USE mo_exception,           ONLY: message, finish
    USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, min_rlcell
    USE mo_model_domain,        ONLY: t_patch
    USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
    USE mo_math_utilities,      ONLY: gc2cc, cc2gc, cc_dot_product,           &
      &                               project_point_to_plane
    USE mo_math_constants,      ONLY: pi, pi_2
    USE mo_io_units,            ONLY: filename_max 
    USE mo_parallel_config,     ONLY: nproma
    USE mo_grid_config,         ONLY: set_patches_grid_filename, n_dom_start, n_dom
    USE mo_loopindices,         ONLY: get_indices_c
    USE mo_intp_lonlat_types,   ONLY: t_lon_lat_intp, dbg_level
    USE mo_mpi,                 ONLY: get_my_mpi_work_id, p_io, p_comm_work,                &
      &                               my_process_is_stdio, p_n_work, p_bcast, p_barrier
    USE mo_communication,       ONLY: idx_1d, blk_no, idx_no
    USE mo_delaunay_types,      ONLY: t_point_list, point_list, point, t_spherical_cap,     &
      &                               spherical_cap, OPERATOR(/), OPERATOR(+), OPERATOR(*), &
      &                               ccw_spherical, t_point, t_triangulation,              &
      &                               triangulation, t_sphcap_list, sagitta_on_unit_sphere
    USE mo_delaunay,            ONLY: point_cloud_diam, triangulate,                        &
      &                               triangulate_mthreaded, create_thin_covering
    USE mo_util_string,         ONLY: int2string
    USE mo_octree,              ONLY: t_range_octree, octree_init, octree_count_point,      &
      &                               octree_query_point, octree_finalize

    IMPLICIT NONE

    !> module name
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_lonlat_baryctr'

    PRIVATE

    PUBLIC :: setup_barycentric_intp_lonlat
    PUBLIC :: setup_barycentric_intp_lonlat_repartition
    PUBLIC :: compute_auxiliary_triangulation
    PUBLIC :: try_triangulation_readin

    INTERFACE compute_auxiliary_triangulation
      MODULE PROCEDURE compute_triangulation_local_partition
      MODULE PROCEDURE compute_triangulation_repartition
    END INTERFACE

  CONTAINS

    !> Generate a sequence of points that are evenly distributed
    !> (approximately) on the sphere.
    !
    !  This spiral approximation formula originates from
    !
    !  E.A Rakhmanov, E.B Saff, Y.M Zhou: Electrons on the
    !  sphere. Series in Approximations and Decompositions, 1994
    !
    SUBROUTINE compute_point_distribution(n, points)
      INTEGER,             INTENT(IN)    :: n         !< number of points to distribute
      TYPE (t_point_list), INTENT(INOUT) :: points
      ! local variables
      TYPE(t_geographical_coordinates) :: points_gc(n)
      TYPE(t_cartesian_coordinates)    :: cc
      TYPE(t_point)                    :: pt
      INTEGER  :: k
      REAL(wp) :: h

      CALL points%initialize()
      DO k=1,n
        h = MIN(MAX(-1._wp + 2._wp*(k-1)/(N-1), -1._wp), 1._wp)
        points_gc(k)%lat = ACOS(h) 
        IF ((k==1) .OR. (k==n)) THEN
          points_gc(k)%lon = 0._wp
        ELSE
          points_gc(k)%lon = MOD( (points_gc(k-1)%lon + 3.6_wp/SQRT(n*(1._wp-h*h))) , 2._wp*pi ) 
        END IF
        points_gc(k)%lon = points_gc(k)%lon - pi
        points_gc(k)%lat = points_gc(k)%lat - pi_2
        cc = gc2cc(points_gc(k))
        pt = point(cc%x(1),cc%x(2),cc%x(3), -1, -1._wp)
        CALL points%push_back(pt)
      END DO
    END SUBROUTINE compute_point_distribution


    !-------------------------------------------------------------------------

    !> compute barycentric coordinates u(1...3) for the point pt
    !> inside the triangle pidx(1...3). The vertex positions are
    !> provided through the array @v vertex.
    !
    SUBROUTINE compute_barycentric_coords(pt, v, u)
      TYPE(t_geographical_coordinates), INTENT(IN)  :: pt                  !< query point (longitude/latitude)  
      !> triangle vertices, first index is x,y,z,
      !! second index is 1st,2nd or 3rd vertex
      REAL(wp),                         INTENT(IN)  :: v(3,3)
      REAL(wp),                         INTENT(OUT) :: u(3)                !< barycentric coordinates (dim: 3)
      ! local variables
      TYPE(t_cartesian_coordinates) :: a,b,c,p
      REAL(wp) :: B1(2), B2(2), r(2), det
      INTEGER :: i

      ! Cartesian coordinates of triangle vertices
      a%x(:) = v(:, 1)
      b%x(:) = v(:, 2)
      c%x(:) = v(:, 3)
      ! query point and project in onto the plane (p1,p2,p3):
      p = project_point_to_plane (gc2cc(pt), a,b,c)
      ! solve linear system for the barycentric coordinates
      !
      ! define matrix and right hand side
      DO i = 1, 3
        a%x(i) = a%x(i) - c%x(i)
        b%x(i) = b%x(i) - c%x(i)
        p%x(i) = p%x(i) - c%x(i)
      END DO
      B1(1:2) = (/ cc_dot_product(a,a) , cc_dot_product(a,b) /)   ! (a-c)^2     ,  (a-c)*(b-c)
      B2(1:2) = (/ B1(2)               , cc_dot_product(b,b) /)   ! (a-c)*(b-c) ,  (b-c)^2
      r(1:2)  = (/ cc_dot_product(p,a) , cc_dot_product(p,b) /)   ! [ (p-c)*(a-c), (p-c)*(b-c) ]
      ! solve using Cramer's rule:
      det     = B2(2)*B1(1) - B1(2)*B2(1) 
      u(1)    = ( r(1)*B2(2) - r(2)*B1(2) )/det
      u(2)    = ( r(2)*B1(1) - r(1)*B2(1) )/det
      u(3)    = 1._wp - (u(1) + u(2))
    END SUBROUTINE compute_barycentric_coords


    !-------------------------------------------------------------------------
    !> Simple test if a point on the unit sphere is inside a triangle.
    !
    FUNCTION inside_triangle(v, v1,v2,v3)
      INTEGER :: inside_triangle
      REAL(wp),       INTENT(IN)     :: v(3)          !< query point
      TYPE(t_point),  INTENT(IN)     :: v1,v2,v3      !< vertex longitudes/latitudes
      ! local variables
      INTEGER       :: c1,c2,c3
      TYPE(t_point) :: p

      p%x = v(1)
      p%y = v(2)
      p%z = v(3)

      c1  = ccw_spherical(v1,v2,p)
      c2  = ccw_spherical(v2,v3,p)
      c3  = ccw_spherical(v3,v1,p)
      IF ((c1 + c2 + c3) == 3 .OR. (c1 + c2 + c3) == -3) THEN
        inside_triangle = 1 ! numerically clearly inside
      ELSEIF (c1 == 0 .OR. c2 == 0 .OR. c3 == 0) THEN
        inside_triangle = 0 ! numerically indefinite
      ELSE
        inside_triangle = -1 ! numerically clearly outside
      ENDIF
    END FUNCTION inside_triangle

    FUNCTION inside_scaled_triangle(v, v1, v2, v3, rscale)
      ! Test if point v is inside a triangle around (v1, v2, v3), that is scaled in size by given rscale.
      ! The area of (v1, v2, v3) is approximately increased to (1+rscale)**2
      INTEGER :: inside_scaled_triangle
      REAL(wp),       INTENT(IN)     :: v(3)          !< query point
      REAL(wp),       INTENT(IN)     :: rscale        !< factor by which to increase the triangle size (>0 for growth)
      TYPE(t_point),  INTENT(IN)     :: v1,v2,v3      !< vertex longitudes/latitudes
      ! local variables
      TYPE(t_point) :: center, g1, g2, g3
      REAL(wp) :: rscale1

      center = (v1 + v2 + v3)/3._wp
      
      ! grow the v-triangle from the center by rscale
      ! \vec{g_i} = \vec{v_i} + scale * (\vec{v_i} - \vec{center}) = (scale+1) * \vec{v_i} - scale * \vec{center}
      rscale1 = 1._wp + rscale
      g1 = rscale1 * v1 + (-rscale) * center
      g2 = rscale1 * v2 + (-rscale) * center
      g3 = rscale1 * v3 + (-rscale) * center

      inside_scaled_triangle = inside_triangle(v, g1, g2, g3)
    END FUNCTION inside_scaled_triangle

    !-------------------------------------------------------------------------
    !> Build a global list of cell circumcenters of "ptr_patch".
    !
    SUBROUTINE create_global_pointlist(ptr_patch, p_global, ldisturb)
      TYPE(t_patch),          INTENT(IN)           :: ptr_patch          !< data structure containing grid info:
      TYPE (t_point_list),    INTENT(INOUT)        :: p_global           !< resulting point set
      LOGICAL,                INTENT(IN)           :: ldisturb
      ! local variables
      TYPE (t_point_list)               :: p_local
      INTEGER                           :: idx, jb, jc, i_startblk, i_endblk,  &
        &                                  rl_start, rl_end, i_startidx,       &
        &                                  i_endidx, i_nchdom, i
      TYPE(t_cartesian_coordinates)     :: p_x, cc     ! Cart. coordinates
      TYPE(t_geographical_coordinates)  :: gc          ! geo. coordinates
      REAL(wp)                          :: perturbation

      ! --- create an array-like data structure containing the local
      ! --- mass points
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# create an array-like data structure containing the local mass points"
      END IF
      CALL p_local%initialize()
      CALL p_local%reserve(ptr_patch%n_patch_cells)

      rl_start   = 1
      rl_end     = min_rlcell_int
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb,        &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
          p_x = gc2cc(ptr_patch%cells%center(jc,jb))
          i = i + 1
          idx = idx_1d(jc,jb)
          CALL p_local%push_back(point(p_x%x(1),p_x%x(2),p_x%x(3), &
            &                          iindex=ptr_patch%cells%decomp_info%glb_index(idx)-1))
        END DO
      END DO
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# min ICON grid coord: ", MINVAL(p_local%a(0:(p_local%nentries-1))%z)
        WRITE (0,*) "# max ICON grid coord: ", MAXVAL(p_local%a(0:(p_local%nentries-1))%z)
      END IF

      ! --- create a global copy of all points

      IF (dbg_level > 10)  WRITE (0,*) "# create a global copy of all points"

      p_global = point_list(p_local)
      CALL p_local%destructor()
      CALL p_global%sync()

      IF (ldisturb) THEN
        ! slightly disturb symmetric coordinates; this should make the
        ! Delaunay triangulation unique, cf. [Lawson1984]
        !
        ! we disturb the triangulation by adding
        !  f(lon,lat) = 1.e-10 * cos(lon) * sin(lat)
        
        IF (dbg_level > 10)  WRITE (0,*) "# slightly disturb symmetric coordinates"
        
!$OMP PARALLEL DO PRIVATE(perturbation,gc,cc)
        DO i=0,(p_global%nentries-1)
          cc%x(1:3)       = (/ p_global%a(i)%x, p_global%a(i)%y, p_global%a(i)%z /)
          gc              = cc2gc(cc)
          perturbation    = COS(gc%lon) * SIN(gc%lat) * 1.e-10_wp
          gc%lon          = gc%lon + perturbation
          cc              = gc2cc(gc)
          p_global%a(i)%x = cc%x(1)
          p_global%a(i)%y = cc%x(2)
          p_global%a(i)%z = cc%x(3)
        END DO
!$OMP END PARALLEL DO
      END IF

    END SUBROUTINE create_global_pointlist


    !-------------------------------------------------------------------------
    !> Build a Delaunay triangulation connecting the cell circumcenters
    !  of "ptr_patch".
    !
    SUBROUTINE compute_triangulation_repartition (ptr_patch, tri_global, p_global, g2l_index)
      TYPE(t_patch),          INTENT(IN)           :: ptr_patch          !< data structure containing grid info:
      TYPE (t_triangulation), INTENT(INOUT)        :: tri_global         !< resulting triangulation
      TYPE (t_point_list),    INTENT(INOUT)        :: p_global           !< resulting point set
      INTEGER, ALLOCATABLE,   INTENT(INOUT)        :: g2l_index(:)       !< point index mapping: global->local
      ! Local parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_triangulation_repartition"
      INTEGER, ALLOCATABLE              :: permutation(:)
      TYPE(t_point)                     :: centroid
      TYPE (t_point_list)               :: pivot_points
      TYPE (t_sphcap_list)              :: subset_list
      INTEGER                           :: errstat, i
      TYPE (t_spherical_cap)            :: subset




      ! --- create an array-like data structure containing the mass points
      CALL create_global_pointlist(ptr_patch, p_global, ldisturb=.TRUE.)

      IF (dbg_level > 10) THEN
        WRITE (0,*) "# total no. of points to triangulate: ", p_global%nentries
      END IF




      CALL p_global%quicksort()

      ALLOCATE(permutation(0:(p_global%nentries-1)), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')

      DO i=0,(p_global%nentries-1)
        IF (p_global%a(i)%gindex /= -1) THEN
          permutation(i) = g2l_index(p_global%a(i)%gindex+1)
        ELSE
          permutation(i) = -1
        END IF
        p_global%a(i)%gindex = i
      END DO
      DEALLOCATE(g2l_index, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')
      CALL MOVE_ALLOC(from=permutation, to=g2l_index)

      ! --- create an auxiliary triangulation, where the vertices are
      ! --- the mass points
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# create an auxiliary triangulation, where the vertices are the mass points"
      END IF
      CALL tri_global%initialize()

      IF (p_n_work < 3) THEN
        ! note: for regular grids the position of the spherical cap
        ! should be one of the poles (meridional convergence of grid
        ! lines)
        centroid = point(0._wp, 0._wp, 1._wp)
        subset = spherical_cap(centroid, -1._wp)
        CALL triangulate_mthreaded(p_global, tri_global, subset, ignore_completeness=.FALSE.)
      ELSE




        ! generate list of spherical cap centers (one for each MPI task)
        IF (dbg_level >= 10) THEN
          WRITE (0,*) "# triangulate pivot points"
        END IF
        CALL compute_point_distribution(p_n_work, pivot_points)
        pivot_points%a(get_my_mpi_work_id())%gindex = get_my_mpi_work_id()
        CALL subset_list%initialize()
        CALL create_thin_covering(pivot_points, subset_list, get_my_mpi_work_id())
        IF (subset_list%nentries > 1)  CALL finish(routine, "Internal error!")
        IF (dbg_level >= 10) THEN
          WRITE (0,*) "# triangulate subset with radius ", subset_list%a(0)%radius
        END IF
        CALL triangulate_mthreaded(p_global, tri_global, subset_list%a(0), ignore_completeness=.FALSE.)

        IF (dbg_level >= 10) THEN
          WRITE (0,'(a,i0,a,a,i0)') "# done. triangulation: ", tri_global%nentries, " triangles.", &
            &                       "; pts = ", p_global%nentries
        END IF
        CALL tri_global%sync()
        ! clean up
        CALL pivot_points%destructor()
        CALL subset_list%destructor()
      END IF
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# no. of cells in auxiliary triangulation: ", tri_global%nentries
      END IF

      ! --- plotting for debugging purposes:
      !
      IF (dbg_level > 20) THEN
        IF (my_process_is_stdio()) THEN
          CALL tri_global%write_vtk("tri_global.vtk", p_global, .FALSE.)
        END IF
        CALL p_barrier(p_comm_work)
      END IF

    END SUBROUTINE compute_triangulation_repartition


    !-------------------------------------------------------------------------
    ! Compute Delaunay triangulation of mass points.
    !
    SUBROUTINE compute_triangulation_local_partition (ptr_patch, tri, p_global)
      ! data structure containing grid info:
      TYPE(t_patch),          INTENT(IN)           :: ptr_patch
      TYPE (t_triangulation), INTENT(INOUT)        :: tri
      TYPE (t_point_list),    INTENT(INOUT)        :: p_global
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_triangulation_local_partition"
      ! enlarge the local triangulation area by this factor
      REAL(wp),     PARAMETER :: RADIUS_FACTOR = 1.1_wp

      INTEGER                         :: jb, jc, i_nchdom,                            &
        &                                i_startidx, i_endidx, i_startblk, i_endblk,  &
        &                                rl_start, rl_end, i, errstat,                &
        &                                idx, nthreads, dim, ierrstat
      TYPE (t_point_list)             :: p_local
      TYPE (t_spherical_cap)          :: subset
      TYPE(t_point)                   :: centroid
      TYPE(t_cartesian_coordinates)   :: p_x
      INTEGER, ALLOCATABLE            :: g2l_index(:)
      TYPE (t_triangulation)          :: tri_global
      INTEGER, ALLOCATABLE            :: permutation(:) ! point index permutation: sorted -> ICON ordering

      !-----------------------------------------------------------------------

      CALL message(routine, '')

      ! --- create an array-like data structure containing the local
      ! --- mass points

      CALL p_local%initialize()
      CALL p_local%reserve(ptr_patch%n_patch_cells)
      rl_start   = 1
      rl_end     = min_rlcell_int
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb,        &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
          p_x = gc2cc(ptr_patch%cells%center(jc,jb))
          i = i + 1
          idx = idx_1d(jc,jb)
          CALL p_local%push_back(point(p_x%x(1),p_x%x(2),p_x%x(3), &
            &                          iindex=ptr_patch%cells%decomp_info%glb_index(idx)))
        END DO
      END DO

      IF (dbg_level > 1) THEN
        WRITE (0,*) "min ICON grid coord: ", MINVAL(p_local%a(0:(p_local%nentries-1))%z)
        WRITE (0,*) "max ICON grid coord: ", MAXVAL(p_local%a(0:(p_local%nentries-1))%z)
      END IF

      ! --- create a translation table global 1D index -> local index
      ALLOCATE(g2l_index(ptr_patch%n_patch_cells_g), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      g2l_index(:) = -1

      rl_start   = 2
      rl_end     = min_rlcell ! note that we include halo cells!
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
!$OMP PARALLEL DO PRIVATE(jb,jc,i_startidx,i_endidx,idx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb,        &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          idx = idx_1d(jc,jb)
          g2l_index(ptr_patch%cells%decomp_info%glb_index(idx)) = idx
        END DO
      END DO
!$OMP END PARALLEL DO

      ! --- create a global copy of all points
      
      p_global = point_list(p_local)
      CALL p_global%sync()

      ALLOCATE(permutation(0:(p_global%nentries-1)), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')

      CALL p_global%quicksort()

      ! slightly disturb symmetric coordinates; this should make the
      ! Delaunay triangulation unique, cf. [Lawson1984]
      dim = 0
      DO i=0,(p_global%nentries-1)
        SELECT CASE(dim)
        CASE (0)
          p_global%a(i)%x = p_global%a(i)%x + 1.e-10_wp
        CASE (1)
          p_global%a(i)%y = p_global%a(i)%y + 1.e-10_wp
        CASE (2)
          p_global%a(i)%z = p_global%a(i)%z + 1.e-10_wp
        END SELECT
        dim = MOD(dim+1,3)
      END DO

      DO i=0,(p_global%nentries-1)
        IF (p_global%a(i)%gindex /= -1) THEN
          permutation(i) = g2l_index(p_global%a(i)%gindex)
        ELSE
          permutation(i) = -1
        END IF
        p_global%a(i)%gindex = i
      END DO
      DEALLOCATE(g2l_index, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

      ! --- create an auxiliary triangulation, where the vertices are
      ! --- the mass points

      ! create a spherical cap around the centroid of local mass
      ! points:
      ! NEC_RP: prevent processes without local data from accessing on allocated data
      IF (p_local%nentries > 0) THEN
        centroid = point(0._wp, 0._wp, 0._wp)
        DO i=0,(p_local%nentries-1)
          centroid = centroid + p_local%a(i)
        END DO
        centroid = centroid/REAL(p_local%nentries,wp)
        subset = spherical_cap(centroid, MIN(RADIUS_FACTOR*point_cloud_diam(p_local, centroid), pi - 1.e-12_wp))      
        IF (dbg_level > 1) THEN
          WRITE (0,*) "spherical cap around ", p_local%a(0)%x, p_local%a(0)%y, p_local%a(0)%z, "; radius ", subset%radius
        END IF
      END IF

      CALL tri%initialize()
      nthreads = 1

      ! for local domains we do not force complete Delaunay
      ! triangulations, since these domains contain pathological
      ! triangles near the boundary which would lead to a
      ! time-consuming triangulation process.

      ! NEC_RP: prevent processes without local data from accessing on allocated data
      IF (p_local%nentries > 0) THEN
        IF (nthreads > 1) THEN
          CALL triangulate_mthreaded(p_global, tri, subset, &
            &                        ignore_completeness = (ptr_patch%id > 1))
        ELSE
          CALL triangulate(p_global, tri, subset, ignore_completeness = (ptr_patch%id > 1))
        END IF
      END IF

      CALL p_local%destructor()

      IF (dbg_level > 1) THEN
        WRITE (0,*) "no. of cells in auxiliary triangulation: ", tri%nentries
      END IF

      DO i=0,(p_global%nentries-1)
        p_global%a(i)%gindex = permutation(i)
      END DO
      DEALLOCATE(permutation, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')


      ! --- plotting for debugging purposes:
      !
      ! --- write a plot of the local triangulation
      IF (dbg_level > 2) THEN
        WRITE (0,'(a,i0,a)') "# formed ", tri%nentries, " triangles."
        IF (dbg_level > 10) THEN
          CALL tri%write_vtk("test"//TRIM(int2string(get_my_mpi_work_id()))//".vtk", p_global, .FALSE.)
        END IF
      END IF

      ! --- write a plot of the global triangulation
      CALL tri%quicksort()
      tri_global=triangulation(tri)
      CALL tri_global%sync()
      IF (my_process_is_stdio() .AND. (dbg_level > 10)) THEN
        CALL tri_global%write_vtk("tri_global.vtk", p_global, .FALSE.)
      END IF

    END SUBROUTINE compute_triangulation_local_partition


    !-------------------------------------------------------------------------
    !> Build an octree out of triangles bounding boxes. Triangles
    !  outside a given range are ommitted.
    !
    SUBROUTINE compute_triangle_bboxes(p_global, tri_global, octree, opt_minrange, opt_maxrange)
      TYPE (t_point_list), TARGET, INTENT(IN) :: p_global           !< set of global points
      TYPE (t_triangulation), INTENT(IN)    :: tri_global           !< global auxiliary triangulation
      REAL(wp), OPTIONAL,     INTENT(IN)    :: opt_minrange(0:2)    !< corner of minimum coords (Cart. coords)
      REAL(wp), OPTIONAL,     INTENT(IN)    :: opt_maxrange(0:2)    !< corner of maximum coords (Cart. coords)
      TYPE (t_range_octree),  INTENT(INOUT) :: octree               !< octree data structure

      ! local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_triangle_bboxes"
      ! local variables
      INTEGER                         :: nlocal_triangles, l, i, j, k, errstat
      REAL(wp)                        :: pmin0(0:2), pmax0(0:2)
      REAL(wp), ALLOCATABLE           :: pmin(:,:), pmax(:,:)
      INTEGER,  ALLOCATABLE           :: glb_index_tri(:)
      REAL(wp)                        :: pp(0:2)
      REAL(wp)                        :: max_sagitta
      REAL(wp)                        :: brange(2,3)          !< box range (min/max, dim=1,2,3)
      TYPE (t_point), POINTER         :: previous_p, current_p
      LOGICAL                         :: llocal_partition

      ! Flag: .TRUE. if only a limited area of the triangulation is
      ! processed:
      llocal_partition = PRESENT(opt_minrange) .AND. PRESENT(opt_maxrange)

      IF (llocal_partition) THEN

        ! --- count the no. of triangles that are not far-off
        IF (dbg_level > 10)  WRITE (0,*) "# count the no. of triangles that are not far-off"

        
        nlocal_triangles = 0 
        ! TODO: OpenMP parallelization
        DO l=0,(tri_global%nentries-1)
          pmin0(:) =  99._wp
          pmax0(:) = -99._wp
          max_sagitta = 0._wp
          current_p => p_global%a(tri_global%a(l)%p(2)) ! initialize with the last point
          DO j=0,2
            previous_p => current_p
            current_p => p_global%a(tri_global%a(l)%p(j))
            pp(0) = current_p%x
            pp(1) = current_p%y
            pp(2) = current_p%z

            DO k=0,2
              pmin0(k) = MIN(pmin0(k), pp(k))
              pmax0(k) = MAX(pmax0(k), pp(k))
            END DO

            ! Determine sagitta of the longest triangle edge
            max_sagitta = MAX(max_sagitta, sagitta_on_unit_sphere(current_p, previous_p))
          END DO

          ! Enlarge the triangle bounding to include the bulge of the triangle.
          ! See comment and descriptions of "bulge height" and "sagitta" below
          pmin0(:) = pmin0(:) - 2._wp*max_sagitta
          pmax0(:) = pmax0(:) + 2._wp*max_sagitta
          
          IF (.NOT. ((pmax0(0) < opt_minrange(0)) .OR. (pmin0(0) > opt_maxrange(0)) .OR.  &
            &        (pmax0(1) < opt_minrange(1)) .OR. (pmin0(1) > opt_maxrange(1)) .OR.  &
            &        (pmax0(2) < opt_minrange(2)) .OR. (pmin0(2) > opt_maxrange(2)))) THEN
            nlocal_triangles = nlocal_triangles + 1
          END IF
        END DO
        IF (dbg_level > 10) THEN
          WRITE (0,*) "# ", nlocal_triangles, " triangles are local to this PE."
        END IF


      ELSE

        nlocal_triangles = tri_global%nentries

      END IF

      ! --- build a list of triangle bounding boxes s.t. we can find
      ! --- the triangles containing our lon-lat points

      IF (dbg_level > 10)  WRITE (0,*) "# build a list of triangle bounding boxes"


      ALLOCATE(pmin(nlocal_triangles,3), pmax(nlocal_triangles,3), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      IF (llocal_partition) THEN
        ALLOCATE(glb_index_tri(nlocal_triangles), STAT=errstat)
        IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      END IF

      i = 1
      DO l=0,(tri_global%nentries-1)
        pmin0(:) =  99._wp
        pmax0(:) = -99._wp
        max_sagitta = 0._wp
        current_p => p_global%a(tri_global%a(l)%p(2)) ! initialize with the last point
        DO j=0,2
          previous_p => current_p
          current_p => p_global%a(tri_global%a(l)%p(j))
          pp(0) = current_p%x
          pp(1) = current_p%y
          pp(2) = current_p%z

          DO k=0,2
            pmin0(k) = MIN(pmin0(k), pp(k))
            pmax0(k) = MAX(pmax0(k), pp(k))
          END DO

          ! Determine sagitta of the longest triangle edge
          max_sagitta = MAX(max_sagitta, sagitta_on_unit_sphere(current_p, previous_p))
        END DO
        ! Enlarge the triangle bounding by the maximum sagitta of all edges to include the bulge 
        ! height of the triangle.
        ! The "bulge height" describes the maximum distance of any point on the spherical triangle 
        ! from the plane trough all vertices.  The bulge height of a spherical triangle might be 
        ! more than the sagitta of its longest edge, but is less then _two times_ this sagitta (3D 
        ! effect).  The sagitta of an edge is distance from the center of the arc (i.e. the edge) to
        ! the center of its base (i.e. n chord). The longest edge however has the smallest dot 
        ! product between a pair of vertices, i.e. min_dot_product.

        pmin0(:) = pmin0(:) - 2._wp*max_sagitta
        pmax0(:) = pmax0(:) + 2._wp*max_sagitta

        IF (.NOT. llocal_partition) THEN
          pmin(i,1:3) = pmin0(0:2)
          pmax(i,1:3) = pmax0(0:2)
          i = i + 1
        ELSE
          IF (.NOT. ((pmax0(0) < opt_minrange(0)) .OR. (pmin0(0) > opt_maxrange(0)) .OR.  &
            &        (pmax0(1) < opt_minrange(1)) .OR. (pmin0(1) > opt_maxrange(1)) .OR.  &
            &        (pmax0(2) < opt_minrange(2)) .OR. (pmin0(2) > opt_maxrange(2)))) THEN
            ! create a translation table "local -> global"
            glb_index_tri(i) = l
            pmin(i,1:3) = pmin0(0:2)
            pmax(i,1:3) = pmax0(0:2)
            i = i + 1
          END IF
        END IF
      END DO

      ! --- insert local triangles into a tree-like data structure

      IF (dbg_level > 10)  WRITE (0,*) "# insert local triangles into a tree-like data structure"
      brange(1,:) = (/ -1._wp, -1._wp, -1._wp /)
      brange(2,:) = (/  1._wp,  1._wp,  1._wp /)

      !NEC_RP: exclude processes without local data from octree_init
      IF (nlocal_triangles > 0) THEN
        IF (llocal_partition) THEN
          CALL octree_init(octree, brange, pmin, pmax, opt_index=glb_index_tri)
          DEALLOCATE(glb_index_tri, STAT=errstat)
          IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')
        ELSE
          CALL octree_init(octree, brange, pmin, pmax)
        END IF
      END IF

      DEALLOCATE(pmin, pmax, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

    END SUBROUTINE compute_triangle_bboxes

    
    !-------------------------------------------------------------------------
    !> Compute barycentric coordinates for a set of points, based on a
    !  given auxiliary triangulation.
    !
    SUBROUTINE compute_barycentric_coordinates(p_global, tri_global, octree, &
      &                                        g2l_index, lcheck_locality, ptr_int_lonlat)
      TYPE (t_point_list),    INTENT(IN)    :: p_global             !< set of global points
      TYPE (t_triangulation), INTENT(IN)    :: tri_global           !< global auxiliary triangulation
      TYPE (t_range_octree),  INTENT(IN)    :: octree               !< octree data structure
      INTEGER,                INTENT(IN)    :: g2l_index(:)         !< point index mapping: global->local
      LOGICAL,                INTENT(IN)    :: lcheck_locality      !< Flag. .FALSE. for relaxed consistency check
      TYPE (t_lon_lat_intp),  INTENT(INOUT) :: ptr_int_lonlat       !< Indices of source points and interpolation coefficients 

      ! local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_barycentric_coordinates"
      ! we use the barycentric coords for the "point in triangle
      ! test"; this is the threshold for this test
      REAL(wp),     PARAMETER :: INSIDETEST_TOL = 1.e-4_wp
      ! the factor by which the interpolation triangles are scaled-up, if no enclosing triangle was found initially
      REAL(wp),     PARAMETER :: UP_SCALE_FACTOR = 1.e6_wp * EPSILON(1._wp) ! eps is 2.2e-16 for dp
      ! max. iterations of UP_SCALE_FACTOR doubling for triangle up-scaling
      INTEGER,      PARAMETER :: NMAX_UP_SCALE = 13
      ! max. no. of triangles (bounding boxes) containing a single lat-lon point.
      INTEGER,      PARAMETER :: NMAX_HITS = 99
      ! local variables
      INTEGER                               :: jb, jc, start_idx, end_idx, nobjects, &
        &                                      i, j, k, idx0, idx1(3), nblks_lonlat, &
        &                                      npromz_lonlat, i_scale
      INTEGER                               :: obj_list(NMAX_HITS)  !< query result (triangle search)
      TYPE(t_cartesian_coordinates)         :: ll_point_c           !< cartes. coordinates of lon-lat points
      REAL(wp)                              :: v(3,0:2)
      LOGICAL                               :: inside_test
      INTEGER                               :: last_idx1(3)

      IF (dbg_level > 10)  WRITE (0,*) "# compute barycentric coordinates"

      ! make sure that the interpolation data structure for the
      ! barycentric interpolation has been allocated:
      IF (.NOT. ALLOCATED(ptr_int_lonlat%baryctr%coeff)) THEN
        CALL finish(routine, "Data structure for the barycentric interpolation not allocated!")
      END IF

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

      ! set dummy indices
      ptr_int_lonlat%baryctr%stencil(:,:)   = 0
      ptr_int_lonlat%baryctr%idx    (:,:,:) = 1
      ptr_int_lonlat%baryctr%blk    (:,:,:) = 1
      ptr_int_lonlat%baryctr%coeff  (:,:,:) = 0._wp


!$OMP PARALLEL DO PRIVATE(jb,jc,start_idx,end_idx,ll_point_c,nobjects,obj_list,      &
!$OMP                     idx0, idx1, v,i,j,k, inside_test, last_idx1)
      DO jb=1,nblks_lonlat
        start_idx = 1
        end_idx   = nproma
        IF (jb == nblks_lonlat) end_idx = npromz_lonlat

        DO jc=start_idx,end_idx

          ! --- determine the triangle in our auxiliary triangulation
          ! --- which contains the lon-lat grid point:

          ! lon-lat point in Cartesian coordinates
          ll_point_c = gc2cc( ptr_int_lonlat%ll_coord(jc,jb) )
          ! query triangles whose bounding boxes contain this point:
          nobjects = octree_count_point(octree, ll_point_c%x(1:3))
          IF (nobjects > NMAX_HITS) THEN
            WRITE (0,*) "point ", ll_point_c%x, " hits ", nobjects, " objects."
            CALL finish(routine, "Internal error!")
          ELSE
            CALL octree_query_point(octree, ll_point_c%x(1:3), obj_list)
          END IF

          ! now test which of the triangles in "obj_list" actually
          ! contains "ll_point_c":
          idx0         = -1
          last_idx1(:) = -1
          ! If no triangle is found in the first round of LOOP_SCALING, increase the size of the test triangles 
          ! by scaling them up. The up-scaling is necessary to solve numerically unclear situations when the
          ! test point is on the edge of a triangle.
          LOOP_SCALING: DO i_scale=0,NMAX_UP_SCALE
            LOOP: DO i=1,nobjects
              j = obj_list(i)
              inside_test = .FALSE.

              IF (i_scale == 0) THEN
                DO k=0,2
                  v(1,k) = p_global%a(tri_global%a(j)%p(k))%x
                  v(2,k) = p_global%a(tri_global%a(j)%p(k))%y
                  v(3,k) = p_global%a(tri_global%a(j)%p(k))%z
                END DO

                ! --- compute the barycentric interpolation weights for
                ! --- this triangle

                CALL compute_barycentric_coords(ptr_int_lonlat%ll_coord(jc,jb), v,               &
                  &                             ptr_int_lonlat%baryctr%coeff(1:3,jc,jb))

                ! test if either the barycentric interpolation weights
                ! indicate that "ll_point_c" lies inside the triangle or
                ! if the test by dot-product succeeds:
                inside_test = ( ALL(ptr_int_lonlat%baryctr%coeff(1:3,jc,jb) >= -1._wp*INSIDETEST_TOL)  .AND. &
                  &             ALL(ptr_int_lonlat%baryctr%coeff(1:3,jc,jb) <=  1._wp+INSIDETEST_TOL))
                IF (.NOT. inside_test) THEN
                  obj_list(i) = -1 ! this candidate must not be considered in the second round of LOOP_SCALING
                  CYCLE LOOP
                END IF

                SELECT CASE (inside_triangle(ll_point_c%x, &
                  &                           p_global%a(tri_global%a(j)%p(0)), &
                  &                           p_global%a(tri_global%a(j)%p(1)), &
                  &                           p_global%a(tri_global%a(j)%p(2))) )
                CASE (1) ! clearly inside
                  inside_test = .TRUE.
                CASE (0) ! unsure. Do nothing and retry with scaling if no other candidate can be found.
                  CYCLE LOOP
                CASE (-1) ! clearly outside
                  obj_list(i) = -1 ! this candidate must not be considered in the second round of LOOP_SCALING
                  CYCLE LOOP
                END SELECT

              ELSE ! i_scale > 0
                IF (j == -1) CYCLE LOOP

                SELECT CASE (inside_scaled_triangle(ll_point_c%x, &
                  &                           p_global%a(tri_global%a(j)%p(0)), &
                  &                           p_global%a(tri_global%a(j)%p(1)), &
                  &                           p_global%a(tri_global%a(j)%p(2)), &
                  &                           2**i_scale * UP_SCALE_FACTOR) )
                CASE (1) ! clearly inside
                  inside_test = .TRUE.
                CASE (0) ! unsure.
                  IF (i_scale == NMAX_UP_SCALE) THEN
                    inside_test = .TRUE. ! accept the numerical undecidable case in the last scaling iteration.
                  ELSE
                    CYCLE LOOP ! Do nothing and retry with more scaling if no other candidate can be found.
                  ENDIF
                CASE (-1) ! clearly outside
                  obj_list(i) = -1 ! this candidate must not be considered in the second round of LOOP_SCALING
                  CYCLE LOOP
                END SELECT

                IF (inside_test) THEN
                  ! If the inside test passed after scaling up the triangle, we have to retrieve the coefficients.
                  DO k=0,2
                    v(1,k) = p_global%a(tri_global%a(j)%p(k))%x
                    v(2,k) = p_global%a(tri_global%a(j)%p(k))%y
                    v(3,k) = p_global%a(tri_global%a(j)%p(k))%z
                  END DO
                  ! --- compute the barycentric interpolation weights for
                  ! --- this enclosing triangle

                  CALL compute_barycentric_coords(ptr_int_lonlat%ll_coord(jc,jb), v,               &
                    &                             ptr_int_lonlat%baryctr%coeff(1:3,jc,jb))

                END IF
              END IF

              IF (inside_test) THEN
                idx0    = j
                ! get global index of cell circumcenters:
                idx1(:) = p_global%a(tri_global%a(idx0)%p(0:2))%gindex + 1
                ! a query point may lie exactly on the edge between two
                ! triangles. We need to make sure that the result is
                ! processor-independent by choosing the triangles with
                ! larger indices.
                IF  (.NOT.  ((idx1(1) >  last_idx1(1))                                    .OR. &
                  &         ((idx1(1) == last_idx1(1)) .AND. (idx1(2) >  last_idx1(2)))   .OR. &
                  &         ((idx1(1) == last_idx1(1)) .AND. (idx1(2) == last_idx1(2)) &
                  &                                    .AND. (idx1(3) >  last_idx1(3))))) THEN
                  CYCLE LOOP
                END IF
                ! convert global to local indices
                idx1(:) = g2l_index(idx1(:))
                last_idx1(:) = idx1(:)

                IF (lcheck_locality) THEN
                  IF (ANY(idx1(:) == -1)) THEN
                    ! the containing triangle is not local for this PE?
                    WRITE (0,*) "indices: ", idx1
                    WRITE (0,*) "baryctr coeffs: ", ptr_int_lonlat%baryctr%coeff(1:3,jc,jb)
                    WRITE (0,*) "lon-lat point: ", ll_point_c%x
                    CALL finish(routine, "Internal error: The containing triangle is not local for this PE!")
                  END IF
                END IF

                EXIT LOOP
              END IF

            END DO LOOP
            ! a matching triangle was found?
            IF (last_idx1(1) /= -1) EXIT LOOP_SCALING
          END DO LOOP_SCALING

          ! no triangle was matching?
          IF (last_idx1(1) == -1) THEN

            IF (lcheck_locality) THEN
              WRITE (0,*) "point ", ll_point_c%x, " hits ", nobjects, " objects."
              CALL finish(routine, "Internal error: no triangle was matching!")
            END IF
            ptr_int_lonlat%baryctr%coeff(1:3,jc,jb) = 0._wp

          ELSE
            CALL compute_barycentric_coords(ptr_int_lonlat%ll_coord(jc,jb), v, &
              &                             ptr_int_lonlat%baryctr%coeff(1:3,jc,jb))

            IF (dbg_level > 5) THEN
              if (i_scale > 0) call message (routine, &
                "Triangle up-scaling was necessary "//TRIM(int2string(i_scale))// &
                " time(s) to achieve a numerical stable result.")
              ptr_int_lonlat%baryctr%v(:,1,jc,jb) = v(:,0)
              ptr_int_lonlat%baryctr%v(:,2,jc,jb) = v(:,1)
              ptr_int_lonlat%baryctr%v(:,3,jc,jb) = v(:,2)
            END IF

            IF (ALL(last_idx1(1:3) >= 1)) THEN
              ! get indices of the containing triangle
              ptr_int_lonlat%baryctr%stencil(jc,jb) = 3            
              ptr_int_lonlat%baryctr%idx(1:3,jc,jb) = idx_no(last_idx1(1:3))
              ptr_int_lonlat%baryctr%blk(1:3,jc,jb) = blk_no(last_idx1(1:3))
            ELSE
              IF (lcheck_locality) THEN                
                WRITE (0,*) "g2l_index(tri_global%a(idx0)%p(0:2)) = ", last_idx1(:)
                CALL finish(routine, "Internal error!")
              END IF
            END IF

          END IF
          
        END DO
      END DO
!$OMP END PARALLEL DO


    END SUBROUTINE compute_barycentric_coordinates


    !-------------------------------------------------------------------------
    !> Try to read the auxiliary triangulation of the cell mass points
    !  from the grid file.
    !
    !  Note: This is an MPI-collective operation!
    !
    INTEGER FUNCTION try_triangulation_readin(ptr_patch, tri_global, p_global)
      TYPE(t_patch),          INTENT(IN)    :: ptr_patch    !< data structure containing grid info
      TYPE (t_triangulation), INTENT(INOUT) :: tri_global
      TYPE (t_point_list),    INTENT(INOUT) :: p_global
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::try_triangulation_readin"
      CHARACTER(LEN=filename_max)       :: grid_filename(n_dom_start:n_dom)
      CHARACTER(LEN=filename_max)       :: grid_filename_grfinfo(n_dom_start:n_dom)
      INTEGER                           :: triangulation_read_from_file

      ! --- try to read triangulation of the mass points from file.
      IF (my_process_is_stdio()) THEN
        CALL tri_global%initialize()
        CALL set_patches_grid_filename(grid_filename(n_dom_start:n_dom), &
          &                            grid_filename_grfinfo(n_dom_start:n_dom))
        triangulation_read_from_file = tri_global%read_netcdf(grid_filename(ptr_patch%id))
      END IF
      CALL p_bcast(triangulation_read_from_file, p_io, p_comm_work)

      IF (triangulation_read_from_file == SUCCESS) THEN

        IF (my_process_is_stdio()) THEN
          WRITE (0,*) routine, ": Auxiliary triangulation read from file '", &
            &         TRIM(grid_filename(ptr_patch%id)), "'."
        ELSE
          CALL tri_global%initialize()
        END IF
        CALL tri_global%bcast(0,p_comm_work)

        ! --- create an array-like data structure containing the mass points
        CALL create_global_pointlist(ptr_patch, p_global, ldisturb=.FALSE.)
        p_global%a(:)%ps = REAL(p_global%a(:)%gindex, wp)
        CALL p_global%quicksort()  ! order point list by their global indices
        IF (dbg_level > 10) THEN
          WRITE (0,*) "synchronization done, ", tri_global%nentries, " triangles."
        END IF

      ELSE

        IF (my_process_is_stdio()) THEN
          WRITE (0,*) routine
        END IF

      END IF
      try_triangulation_readin = triangulation_read_from_file
    END FUNCTION try_triangulation_readin


    !-------------------------------------------------------------------------
    ! Setup routine for barycentric interpolation at lon-lat grid
    ! points for an arbitrary grid.
    !
    SUBROUTINE setup_barycentric_intp_lonlat_repartition(ptr_patch, tri_global, p_global,  &
      &                                                  triangulation_read_from_file,     &
      &                                                  ptr_int_lonlat)

      ! data structure containing grid info:
      TYPE(t_patch),                 INTENT(IN)    :: ptr_patch
      ! data structure for global triangulation of mass points
      TYPE (t_triangulation),        INTENT(INOUT) :: tri_global
      ! data structure for global set of mass points
      TYPE (t_point_list),           INTENT(INOUT) :: p_global
      ! Flag. .TRUE. if triangulation is already initialized
      INTEGER,                       INTENT(IN)    :: triangulation_read_from_file
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat

      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::setup_barycentric_intp_lonlat_repartition"

      INTEGER                           :: jb, jc, start_idx, end_idx,       &
        &                                  i, errstat, idx,                  &
        &                                  i_startblk, i_endblk,             &
        &                                  rl_start, rl_end, i_startidx,     &
        &                                  i_endidx, i_nchdom, nblks_lonlat, &
        &                                  npromz_lonlat, pidx
      INTEGER, ALLOCATABLE              :: g2l_index(:)
      TYPE (t_range_octree)             :: octree               !< octree data structure
      TYPE(t_cartesian_coordinates)     :: ll_point_c, p_c      !< cartes. coordinates of lon-lat points
      REAL(wp)                          :: minrange(0:2), maxrange(0:2),diff(3), maxval
      LOGICAL                           :: lcheck_locality

      CALL message(routine, '')

      ! --- create a translation table global 1D index -> local index

      IF (dbg_level > 10)  WRITE (0,*) "# create a translation table global 1D index -> local index"
      ALLOCATE(g2l_index(ptr_patch%n_patch_cells_g), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      g2l_index(:) = -1

      rl_start   = 2
      rl_end     = min_rlcell ! note that we include halo cells
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          i = i + 1
          idx = idx_1d(jc,jb)
          g2l_index(ptr_patch%cells%decomp_info%glb_index(idx)) = idx
        END DO
      END DO

      IF (triangulation_read_from_file /= SUCCESS) THEN

        ! --- compute a triangulation of the mass points; this modifies
        !     the "g2l_index" array since the point ordering is altered.
        CALL compute_auxiliary_triangulation(ptr_patch, tri_global, p_global, g2l_index)

        ! --- since the point set for the auxiliary triangulation had
        !     been distorted, we need to restore the original triangle
        !     circumcenters:

        maxval = 0._wp        !
!$OMP PARALLEL DO PRIVATE(i,pidx,p_c,diff), REDUCTION(MAX:maxval)
        DO i=0,(p_global%nentries-1)
          pidx = g2l_index(p_global%a(i)%gindex) 
          IF (pidx > 0) THEN
            p_c  = gc2cc(ptr_patch%cells%center(idx_no(pidx),blk_no(pidx)))
            ! consistency check: compute distance between old and new:
            diff(:) = (/ p_global%a(i)%x,p_global%a(i)%y,p_global%a(i)%z /) - p_c%x(:)
            maxval = MAX(maxval, SQRT(DOT_PRODUCT(diff,diff)))
            p_global%a(i)%x = p_c%x(1)
            p_global%a(i)%y = p_c%x(2)
            p_global%a(i)%z = p_c%x(3)
          END IF
        END DO
!$OMP END PARALLEL DO
        IF (maxval > 1.e-5_wp) THEN
          WRITE (0,*) "maxval = ", maxval
          CALL finish(routine, "Internal error!")
        END IF
      END IF

      ! --- determine bounding box of query points

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

      ! TODO: OpenMP parallelization
      minrange(:) =  99._wp
      maxrange(:) = -99._wp
      DO jb=1,nblks_lonlat
        start_idx = 1
        end_idx   = nproma
        IF (jb == nblks_lonlat) end_idx = npromz_lonlat

        DO jc=start_idx,end_idx
          ! lon-lat point in Cartesian coordinates
          ll_point_c = gc2cc( ptr_int_lonlat%ll_coord(jc,jb) )

          DO i=0,2
            minrange(i) = MIN(minrange(i), ll_point_c%x(i+1))
            maxrange(i) = MAX(maxrange(i), ll_point_c%x(i+1))
          END DO
        END DO
      END DO

      ! --- determine bounding box of triangles
      CALL compute_triangle_bboxes(p_global, tri_global, octree, minrange, maxrange)

      ! --- the containing triangle is not local for this PE;
      !     this may happen for nested regions; we therefore do not stop.
      lcheck_locality = (ptr_patch%id <= 1)

      ! --- compute barycentric coordinates
      CALL compute_barycentric_coordinates(p_global, tri_global, octree,  &
        &                                  g2l_index, lcheck_locality, ptr_int_lonlat)

      ! clean up
      !NEC_RP: only call finalization when octree is allocated 
      IF (ALLOCATED(octree%box)) THEN
        CALL octree_finalize(octree)
      END IF

      DEALLOCATE(g2l_index, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

    END SUBROUTINE setup_barycentric_intp_lonlat_repartition


    !-------------------------------------------------------------------------
    ! Setup routine for barycentric interpolation at lon-lat grid
    ! points for an arbitrary grid.
    !
    SUBROUTINE setup_barycentric_intp_lonlat(tri, p_global, ptr_int_lonlat)
      ! triangulation of mass points.
      TYPE (t_triangulation),        INTENT(IN)    :: tri
      TYPE (t_point_list),           INTENT(IN)    :: p_global
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp),         INTENT(INOUT) :: ptr_int_lonlat
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::setup_barycentric_intp_lonlat"
      ! max. no. of triangle (bounding boxes) containing a single lat-lon point.
      INTEGER,      PARAMETER :: NMAX_HITS = 99
      ! we use the barycentric coords for the "point in triangle
      ! test"; this is the threshold for this test
      REAL(wp),     PARAMETER :: INSIDETEST_TOL = 1.e-6

      INTEGER                         :: nblks_lonlat, npromz_lonlat, jb, jc,    &
        &                                i_startidx, i_endidx, i, j, iv,        &
        &                                nobjects, idx0
      REAL(wp)                        :: v(3,3), baryctr_coeff(3)
      TYPE (t_range_octree)           :: octree               !< octree data structure
      INTEGER                         :: obj_list(NMAX_HITS)  !< query result (triangle search)
      TYPE(t_cartesian_coordinates)   :: ll_point_c           !< cartes. coordinates of lon-lat points
      LOGICAL                         :: inside_test1, inside_test2

      !-----------------------------------------------------------------------

      CALL message(routine, '')

      ! make sure that the interpolation data structure for the
      ! barycentric interpolation has been allocated:
      IF (.NOT. ALLOCATED(ptr_int_lonlat%baryctr%coeff)) THEN
        CALL finish(routine, "Data structure for the barycentric interpolation not allocated!")
      END IF

      ! set dummy indices
      ptr_int_lonlat%baryctr%stencil(:,:)   = 0
      ptr_int_lonlat%baryctr%idx    (:,:,:) = 1
      ptr_int_lonlat%baryctr%blk    (:,:,:) = 1
      ptr_int_lonlat%baryctr%coeff  (:,:,:) = 0._wp

      ! --- determine bounding box of triangles
      CALL compute_triangle_bboxes(p_global, tri, octree)

      ! --- compute barycentric coordinates

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = ptr_int_lonlat%nblks_lonlat(nproma)
      npromz_lonlat = ptr_int_lonlat%npromz_lonlat(nproma)

!$OMP PARALLEL DO PRIVATE(jb,jc,i_startidx,i_endidx,ll_point_c,nobjects,obj_list, &
!$OMP                     idx0, v,i,j,inside_test1,inside_test2, baryctr_coeff )
      DO jb=1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc=i_startidx,i_endidx

          ! --- determine the triangle in our auxiliary triangulation
          ! --- which contains the lon-lat grid point:

          ! lon-lat point in Cartesian coordinates
          ll_point_c = gc2cc( ptr_int_lonlat%ll_coord(jc,jb) )
          ! query triangles whose bounding boxes contain this point:
          nobjects = octree_count_point(octree, ll_point_c%x(1:3))
          IF (nobjects > NMAX_HITS) THEN
            WRITE (0,*) "point ", ll_point_c%x, " hits ", nobjects, " objects."
            CALL finish(routine, "Internal error!")
          ELSE
            CALL octree_query_point(octree, ll_point_c%x(1:3), obj_list)
          END IF

          ! now test which of the triangles in "obj_list" actually
          ! contains "ll_point_c":
          idx0 = -1
          LOOP: DO i=1,nobjects
            j = obj_list(i) - 1
            DO iv = 1, 3
              v(1,iv) = p_global%a(tri%a(j)%p(iv-1))%x
              v(2,iv) = p_global%a(tri%a(j)%p(iv-1))%y
              v(3,iv) = p_global%a(tri%a(j)%p(iv-1))%z
            END DO
            ! --- compute the barycentric interpolation weights for
            ! --- this triangle
            baryctr_coeff(:) = 0._wp
            CALL compute_barycentric_coords(ptr_int_lonlat%ll_coord(jc,jb), v, baryctr_coeff)

            IF (dbg_level > 5) ptr_int_lonlat%baryctr%v(:,:,jc,jb) = v

            ! test if either the barycentric interpolation weights
            ! indicate that "ll_point_c" lies inside the triangle or
            ! if the test by dot-product succeeds:
            inside_test1 = ( ALL(baryctr_coeff >= -1._wp*INSIDETEST_TOL)  .AND. &
              &              ALL(baryctr_coeff <=  1._wp+INSIDETEST_TOL))
            inside_test2 = (1 == inside_triangle(ll_point_c%x, p_global%a(tri%a(j)%p(0)), p_global%a(tri%a(j)%p(1)), &
              &                            p_global%a(tri%a(j)%p(2))) )

            IF (inside_test1 .OR. inside_test2) THEN
              idx0 = j

              IF (ALL(p_global%a(tri%a(idx0)%p(0:2))%gindex >= 1)) THEN
                ! get indices of the containing triangle
                ptr_int_lonlat%baryctr%stencil(jc,jb)   = 3
                ptr_int_lonlat%baryctr%coeff(1:3,jc,jb) = baryctr_coeff(:)
                ptr_int_lonlat%baryctr%idx(1:3,jc,jb)   = idx_no(p_global%a(tri%a(idx0)%p(0:2))%gindex)
                ptr_int_lonlat%baryctr%blk(1:3,jc,jb)   = blk_no(p_global%a(tri%a(idx0)%p(0:2))%gindex)

                IF (ANY(p_global%a(tri%a(idx0)%p(0:2))%gindex < 1)) THEN
                  WRITE (0,*) "permutation(tri%a(idx0)%p(0:2)) = ", p_global%a(tri%a(idx0)%p(0:2))%gindex
                  CALL finish(routine, "Internal error!")
                END IF
              ELSE
                ! the containing triangle is not local for this PE;
                ! this may happen for nested regions; we therefore do not stop.
              END IF

              EXIT LOOP
            END IF
          END DO LOOP

        END DO
      END DO
!$OMP END PARALLEL DO

      ! clean up
      !NEC_RP: only call finalization when octree is allocated 
      IF (ALLOCATED(octree%box)) THEN
        CALL octree_finalize(octree)
      END IF

    END SUBROUTINE setup_barycentric_intp_lonlat

  END MODULE mo_intp_lonlat_baryctr
