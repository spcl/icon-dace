! Type and object declarations, operators and utility routines for
! the construction of a Delaunay
! triangulation on the sphere.
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

MODULE mo_delaunay_types

#ifndef NOMPI
  USE mpi
#endif

#ifdef _OPENMP
  USE OMP_LIB
#endif

  USE mo_util_file,         ONLY: util_file_is_writable
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf
  USE mo_exception,         ONLY: finish
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_kind,              ONLY: wp
  USE mo_mpi,               ONLY: p_comm_work, p_real_dp
#ifdef __SX__
  USE mo_util_sort,         ONLY: radixsort
#endif
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_edge, t_point, t_triangle
  PUBLIC :: t_spherical_cap
  PUBLIC :: t_triangulation
  PUBLIC :: t_point_list, t_sphcap_list
  PUBLIC :: point, triangle, point_list, triangulation, triangulation_ptr, spherical_cap
  PUBLIC :: OPERATOR(<), OPERATOR(/), OPERATOR(*), OPERATOR(==), OPERATOR(+)
  PUBLIC :: circum_circle_spherical, ccw_spherical, sagitta_on_unit_sphere

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_delaunay_types'

  ! quadruple precision, needed for some determinant computations
#if ( defined __PGI )
#define _NO_QUAD_PRECISION
#elif ( ! defined NAGFOR && ! defined _SX )
  INTEGER, PARAMETER :: QR_K = SELECTED_REAL_KIND (32)
#else
  INTEGER, PARAMETER :: QR_K = SELECTED_REAL_KIND (2*precision(1.0_wp))
#endif


  ! --------------------------------------------------------------------
  ! OBJECT DECLARATIONS
  ! --------------------------------------------------------------------

  !> type declaration: single triangle edge
  TYPE t_edge
    INTEGER        :: p1,p2                       !< edge vertex indices
  END TYPE t_edge

  !> type declaration: single point in R^3
  TYPE t_point
    REAL(wp)       :: x,y,z, ps
    INTEGER        :: gindex                      !< global index
  CONTAINS
    PROCEDURE :: norm2          => point_norm2
    PROCEDURE :: spherical_dist => point_spherical_dist
  END TYPE t_point

  !> type declaration: single triangle
  TYPE t_triangle
    INTEGER        :: p(0:2)                      !< vertex indices
    INTEGER        :: oedge                       !< local index of edge opposite to ghost point (if present)
    INTEGER        :: complete                    !< 0: flag not set; 1: triangle must not be visited again; 2: triangle invalid
    TYPE(t_point)  :: cc                          !< circumcircle center
    REAL(wp)       :: r                           !< circumcircle radius
    REAL(wp)       :: cap_distance                !< min distance from subset center to circumcircle
    REAL(wp)       :: rdiscard                    !< max distance from subset center, projected onto sort direction
  CONTAINS
    PROCEDURE, PUBLIC :: edge                 => triangle_edge
    PROCEDURE, PUBLIC :: compute_circumcenter => triangle_compute_circumcenter
    PROCEDURE, PUBLIC :: normalize_indices    => triangle_normalize_indices
  END TYPE t_triangle

  TYPE :: t_mpi_point
    SEQUENCE
    REAL(wp) :: x,y,z,ps
    INTEGER  :: gindex, padding
  END TYPE t_mpi_point

  TYPE :: t_mpi_triangle
    SEQUENCE
    INTEGER  :: p(0:2)
  END TYPE t_mpi_triangle

  !> type declaration: spherical cap
  TYPE t_spherical_cap
    TYPE (t_point) :: point, sorting_direction
    REAL(wp)       :: radius                      !< cap radius (radius<0 means "global")
  END TYPE t_spherical_cap

  !> type declaration: list of abstract objects
  TYPE, ABSTRACT :: t_object_list
    INTEGER :: nentries
  CONTAINS
    PRIVATE
    PROCEDURE, PUBLIC :: initialize => initialize_object_list
    PROCEDURE(routine1),  DEFERRED, PASS :: destructor !< destructor... not all compilers support "FINAL"
  END TYPE t_object_list

  !> type declaration: list of abstract objects
  TYPE, ABSTRACT, EXTENDS(t_object_list) :: t_sortable_list
  CONTAINS
    PRIVATE
    PROCEDURE(arrayfct),    DEFERRED, PASS :: swap       !< swap operator for generic sort algorithm
    PROCEDURE(binaryop),    DEFERRED, PASS :: cmp        !< "<" operator for generic sort algorithm
    PROCEDURE(routine3),    DEFERRED, PASS :: reserve    !< re-allocate data structure
    PROCEDURE, PUBLIC :: quicksort  => quicksort_sortable_list
  END TYPE t_sortable_list

  !> type declaration: list of points
  TYPE, EXTENDS(t_sortable_list) :: t_point_list
    TYPE(t_point), ALLOCATABLE :: a(:)
  CONTAINS
    PROCEDURE :: swap       => swap_point_list
    PROCEDURE :: cmp        => cmp_point_list
    PROCEDURE :: destructor => destructor_point_list
    PROCEDURE :: reserve    => reserve_point_list
    PROCEDURE :: resize     => resize_point_list
    PROCEDURE :: push_back  => push_back_point_list
    PROCEDURE :: sync       => sync_point_list
    PROCEDURE, PUBLIC :: quicksort  => quicksort_point_list
#ifdef __SX__
    PROCEDURE, PUBLIC :: radixsort => radixsort_point_list
#endif
  END TYPE t_point_list

  !> type declaration: list of triangles
  TYPE, EXTENDS(t_sortable_list) :: t_triangulation
    TYPE(t_triangle), ALLOCATABLE :: a(:)
  CONTAINS
    PROCEDURE :: swap         => swap_triangulation
    PROCEDURE :: cmp          => cmp_triangulation
    PROCEDURE :: destructor   => destructor_triangulation
    PROCEDURE :: reserve      => reserve_triangulation
    PROCEDURE :: resize       => resize_triangulation
    PROCEDURE :: push_back    => push_back_triangulation
    PROCEDURE :: sync         => sync_triangulation
    PROCEDURE :: bcast        => bcast_triangulation

    PROCEDURE :: write_vtk    => write_vtk_triangulation
    PROCEDURE :: write_netcdf => write_netcdf_triangulation
    PROCEDURE :: read_netcdf  => read_netcdf_triangulation
  END TYPE t_triangulation

  !> type declaration: list of triangles
  TYPE, EXTENDS(t_object_list) :: t_sphcap_list
    TYPE(t_spherical_cap), ALLOCATABLE :: a(:)
  CONTAINS
    PROCEDURE :: destructor => destructor_sphcap_list
    PROCEDURE :: reserve    => reserve_sphcap_list
    PROCEDURE :: resize     => resize_sphcap_list
    PROCEDURE :: push_back  => push_back_sphcap_list
  END TYPE t_sphcap_list

  TYPE t_triptr
    TYPE(t_triangulation), POINTER :: ptr
  END TYPE t_triptr


  ! --------------------------------------------------------------------
  ! DATA TYPES FOR K-WAY MERGE ALGORITHM
  ! --------------------------------------------------------------------

  TYPE t_min_heap_elt
    TYPE(t_triangle) :: p
  END TYPE t_min_heap_elt

  ! min heap node
  TYPE t_min_heap_node
    TYPE (t_min_heap_elt) :: elt ! element to be stored
    INTEGER               :: i   ! array from which the element is taken
    INTEGER               :: j   ! index of the next element to be picked 
  END TYPE t_min_heap_node

  TYPE t_min_heap
    TYPE (t_min_heap_node), ALLOCATABLE :: harr(:)  ! array of elements in heap
    INTEGER                             :: isize    ! size of min heap    
  END TYPE t_min_heap



  ! --------------------------------------------------------------------
  ! DECLARATION OF OPERATORS AND TYPE-BOUND PROCEDURES
  ! --------------------------------------------------------------------

  ABSTRACT INTERFACE
    SUBROUTINE routine1(this)
      IMPORT :: t_object_list
      CLASS(t_object_list), INTENT(INOUT) :: this
    END SUBROUTINE routine1
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE routine3(this, new_size)
      IMPORT :: t_sortable_list
      CLASS(t_sortable_list), INTENT(INOUT) :: this
      INTEGER,                INTENT(IN)    :: new_size
    END SUBROUTINE routine3
  END INTERFACE

  ABSTRACT INTERFACE
    LOGICAL PURE FUNCTION binaryop(this, i,j)
      IMPORT :: t_sortable_list
      CLASS(t_sortable_list), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: i,j
    END FUNCTION binaryop
  END INTERFACE

  ABSTRACT INTERFACE
    PURE SUBROUTINE arrayfct(this, i,j)
      IMPORT :: t_sortable_list
      CLASS(t_sortable_list), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: i,j
    END SUBROUTINE arrayfct
  END INTERFACE

  !> operator: comparison of two edges
  INTERFACE OPERATOR ( == )
    MODULE PROCEDURE t_edge_equal
    MODULE PROCEDURE t_point_equal
    MODULE PROCEDURE t_triangle_equal
  END INTERFACE

  !> operator: subtraction of two points
  INTERFACE OPERATOR ( - )
    MODULE PROCEDURE t_point_minus
  END INTERFACE

  !> operator: subtraction of two points
  INTERFACE OPERATOR ( + )
    MODULE PROCEDURE t_point_plus
  END INTERFACE

  !> operator: dot product and product with scalar
  INTERFACE OPERATOR ( * )
    MODULE PROCEDURE t_point_mul
    MODULE PROCEDURE t_point_mul_scalar
  END INTERFACE

  !> operator: division of point by scalar
  INTERFACE OPERATOR ( / )
    MODULE PROCEDURE t_point_div_scalar
  END INTERFACE

  !> operator: comparison of points and triangles (wrt. sorting direction)
  INTERFACE OPERATOR ( < )
    MODULE PROCEDURE t_point_cmp
    MODULE PROCEDURE t_triangle_cmp
  END INTERFACE

CONTAINS

  ! --------------------------------------------------------------------
  ! GENERAL GEOMETRIC SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> Return TRUE if a point (xp,yp) is inside the circumcircle made up
  !  of the points ip[1,2,3].
  ! 
  !  See, e.g., 
  !  Renka, R. J. Interpolation of Data on the Surface of a Sphere
  !               ACM Trans. Math. Softw., ACM, 1984, 10, 417-436
  !  Renka's STRIPACK algorithm (http://www.netlib.org/toms/772)
#ifndef _NO_QUAD_PRECISION
  PURE &
#endif
  FUNCTION circum_circle_spherical(p, pxyz, ip)
    LOGICAL :: circum_circle_spherical
    TYPE (t_point),      INTENT(IN)   :: p
    TYPE (t_point_list), INTENT(IN)   :: pxyz
    INTEGER,             INTENT(IN)   :: ip(0:2)
    ! local variables
    REAL(wp) :: d1_x, d1_y, d1_z,d2_x, d2_y, d2_z,d3_x, d3_y, d3_z, d1, d2

    ! p lies above the plane of (p1,p3,p2) iff p2 lies above the plane
    ! of (p3,p1,p) iff Det(p2-p,p3-p,p1-p) = (p2-p,p3-p X p1-p) > 0.
    !
    ! det =      (pxyz%a(ip(1))%x - p%x)*(d3%y*d1%z - d1%y*d3%z) &
    !   &     -  (pxyz%a(ip(1))%y - p%y)*(d3%x*d1%z - d1%x*d3%z) &
    !   &     +  (pxyz%a(ip(1))%z - p%z)*(d3%x*d1%y - d1%x*d3%y)
    d1_x = pxyz%a(ip(0))%x - p%x
    d1_y = pxyz%a(ip(0))%y - p%y
    d1_z = pxyz%a(ip(0))%z - p%z
    d2_x = pxyz%a(ip(1))%x - p%x
    d2_y = pxyz%a(ip(1))%y - p%y
    d2_z = pxyz%a(ip(1))%z - p%z
    d3_x = pxyz%a(ip(2))%x - p%x
    d3_y = pxyz%a(ip(2))%y - p%y
    d3_z = pxyz%a(ip(2))%z - p%z

    d1 = d2_x*d3_y*d1_z  +  d2_y*d1_x*d3_z + d2_z*d3_x*d1_y
    d2 = d2_y*d3_x*d1_z + d2_x*d1_y*d3_z + d2_z*d1_x*d3_y 
    IF (ABS(d1-d2) < 1.e-12_wp) THEN
      circum_circle_spherical = circum_circle_spherical_q128(p, pxyz, ip)
    ELSE
      circum_circle_spherical = (d1 > d2)
    END IF
  END FUNCTION circum_circle_spherical


  ! --------------------------------------------------------------------
  !> Return TRUE if a point (xp,yp) is inside the circumcircle made up
  !  of the points ip[1,2,3].
  ! 
  !  See, e.g., 
  !  Renka, R. J. Interpolation of Data on the Surface of a Sphere
  !               ACM Trans. Math. Softw., ACM, 1984, 10, 417-436
  !  Renka's STRIPACK algorithm (http://www.netlib.org/toms/772)
#ifndef _NO_QUAD_PRECISION
  PURE FUNCTION circum_circle_spherical_q128(p, pxyz, ip)
    LOGICAL :: circum_circle_spherical_q128
    TYPE (t_point),      INTENT(IN)   :: p
    TYPE (t_point_list), INTENT(IN)   :: pxyz
    INTEGER,             INTENT(IN)   :: ip(0:2)
    ! local variables
    REAL(QR_K) :: d1_x, d1_y, d1_z,d2_x, d2_y, d2_z,d3_x, d3_y, d3_z

    ! p lies above the plane of (p1,p3,p2) iff p2 lies above the plane
    ! of (p3,p1,p) iff Det(p2-p,p3-p,p1-p) = (p2-p,p3-p X p1-p) > 0.
    !
    ! det =      (pxyz%a(ip(1))%x - p%x)*(d3%y*d1%z - d1%y*d3%z) &
    !   &     -  (pxyz%a(ip(1))%y - p%y)*(d3%x*d1%z - d1%x*d3%z) &
    !   &     +  (pxyz%a(ip(1))%z - p%z)*(d3%x*d1%y - d1%x*d3%y)
    d1_x = pxyz%a(ip(0))%x - p%x
    d1_y = pxyz%a(ip(0))%y - p%y
    d1_z = pxyz%a(ip(0))%z - p%z
    d2_x = pxyz%a(ip(1))%x - p%x
    d2_y = pxyz%a(ip(1))%y - p%y
    d2_z = pxyz%a(ip(1))%z - p%z
    d3_x = pxyz%a(ip(2))%x - p%x
    d3_y = pxyz%a(ip(2))%y - p%y
    d3_z = pxyz%a(ip(2))%z - p%z

    IF ( d2_x*d3_y*d1_z  +  d2_y*d1_x*d3_z + d2_z*d3_x*d1_y &
      > d2_y*d3_x*d1_z + d2_x*d1_y*d3_z + d2_z*d1_x*d3_y ) THEN 
      circum_circle_spherical_q128 = .TRUE.
    ELSE
      circum_circle_spherical_q128 = .FALSE.
    END IF
  END FUNCTION circum_circle_spherical_q128
#else
  FUNCTION circum_circle_spherical_q128(p, pxyz, ip)
    LOGICAL :: circum_circle_spherical_q128
    TYPE (t_point),      INTENT(IN)   :: p
    TYPE (t_point_list), INTENT(IN)   :: pxyz
    INTEGER,             INTENT(IN)   :: ip(0:2)
    CALL finish("circum_circle_spherical_q128", &
      "This function cannot be used without quadruple precision support. " // &
      "Please consider applying `icondeluny` from the DWD ICON Tools on " // &
      "the grid files manually to pre-generate the variable `cc_delaunay`.")
    circum_circle_spherical_q128 = .FALSE.
  END FUNCTION circum_circle_spherical_q128
#endif


  ! --------------------------------------------------------------------
  !> Locates a point relative to a directed arc. Let v1, v2, and v3 be
  !  distinct points on the sphere, and denote by v1 -> v2 the
  !  geodesic connecting v1 and v2 and directed toward v2. This test
  !  determines which of the two hemispheres defined by v1 -> v2
  !  contains v3, or, in other words, if we "turn left" when going
  !  from v1 to v2 to v3.
  PURE FUNCTION ccw_spherical(v1,v2,v3)
    INTEGER :: ccw_spherical ! -1: clearly left, 0: numerically uncertain, 1: clearly right
    TYPE (t_point), INTENT(IN)  :: v1,v2,v3
    REAL(wp), PARAMETER :: FLOAT_ERROR_MARGIN = 1.1e-14_wp
    REAL(wp) :: ccw

    ! det(v1,v2,v3) = <v1 x v2, v3> = | v1 x v2 | cos(a) 
    !  
    ! where a is the angle between v3 and the normal to the plane
    ! defined by v1 and v2.
    
    ccw =           v3%x*(v1%y*v2%z - v2%y*v1%z) &
      &        -    v3%y*(v1%x*v2%z - v2%x*v1%z) &
      &        +    v3%z*(v1%x*v2%y - v2%x*v1%y) 


    ! we apply a static error of
    !   | e - e'| <= 3*2^-48
    ! to decide if a floating-point evaluation e' of an expression e
    ! has the correct sign, see Section 2.2 of
    !
    ! Burnikel, C.; Funke, S. & Seel, M. 
    ! "Exact geometric computation using Cascading"
    ! International Journal of Computational Geometry & Applications, 
    ! World Scientific, 2001, 11, 245-266
    IF (ccw < -FLOAT_ERROR_MARGIN) THEN
      ccw_spherical = -1
    ELSE IF (ccw > FLOAT_ERROR_MARGIN) THEN
      ccw_spherical = 1
    ELSE
      ccw_spherical = 0
    END IF
  END FUNCTION ccw_spherical



  ! --------------------------------------------------------------------
  !> Computes sagitta of great circle arc, which connects points a1 and a2 on a unit sphere
  REAL(wp) PURE FUNCTION sagitta_on_unit_sphere(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1, a2
    ! The sagitta of an arc is distance from the center of the arc to the center of its base.
    ! The base is the chord that connects a1 and a2.
    ! 
    ! This function assumes 3D points on a unit sphere.
    !
    !       a1          Let 2*q be the chord length between a1 and a2.
    !     , '\          Let s be sagitta that we are looking for and h the height of the 
    !   ,   | \         triangle (center, a1, a2).
    !  ,   q|  \ 1      Then s = 1 - h
    ! ,     |   \         len(a1) = len(a2) = 1
    ! ,_____|___ \cen     h**2 + q**2 = 1**2 ! Pythagoras
    ! , s   | h  /ter     (2*q)**2 = len(a1)**2 + len(a2)**2 - 2 * len(a1) * len(a1) &
    !  ,   q|   /                  & * cos(angle(a1,a2)) ! Law of cosines
    !   ,   |  /1         cos(angle(a1,a2)) = abs(dot(a1, a2)) / (len(a1) * len(a2)) &
    !     , | /                    & ! definition of scalar product
    !       '/          Combine and simplify: q**2 = 0.5 - 0.5 * abs(dot(a1, a2))
    !       a2          Thus sagitta s = 1 - SQRT(0.5 * (1 + ABS(dot(a1, a2))))

    sagitta_on_unit_sphere = 1._wp - SQRT( 0.5_wp * (1._wp + ABS(a1 * a2)) )
  END FUNCTION sagitta_on_unit_sphere


  ! --------------------------------------------------------------------
  ! DEFINITION OF OVERLOADED OPERATORS
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> operator: comparison of two edges
  LOGICAL PURE FUNCTION t_edge_equal(e1,e2)
    TYPE(t_edge), INTENT(IN) :: e1,e2
    ! note: we test opposite directions first, since this is the
    ! situation for neighboring triangles
    t_edge_equal =                                            &
      &    ( (e1%p1 == e2%p2) .AND. (e1%p2 == e2%p1) )  .OR.  &  
      &    ( (e1%p1 == e2%p1) .AND. (e1%p2 == e2%p2) )
  END FUNCTION t_edge_equal


  ! --------------------------------------------------------------------
  !> operator: comparison of two triangles
  LOGICAL PURE FUNCTION t_triangle_equal(t1,t2)
    TYPE(t_triangle), INTENT(IN) :: t1,t2
    t_triangle_equal =                                        &
      &    ( (t1%p(0) == t2%p(0)) .AND. (t1%p(1) == t2%p(1)) .AND. (t1%p(2) == t2%p(2)) )
  END FUNCTION t_triangle_equal


  ! --------------------------------------------------------------------
  !> operator: subtraction of two points
  TYPE(t_point) PURE FUNCTION t_point_minus(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_minus = t_point(x=a1%x-a2%x, y=a1%y-a2%y, z=a1%z-a2%z, &
      &                     ps=-2._wp, gindex=0)
  END FUNCTION t_point_minus


  ! --------------------------------------------------------------------
  !> operator: addition of two points
  TYPE(t_point) PURE FUNCTION t_point_plus(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_plus = t_point(x=a1%x+a2%x, y=a1%y+a2%y, z=a1%z+a2%z, &
      &                    ps=-2._wp, gindex=0)
  END FUNCTION t_point_plus


  ! --------------------------------------------------------------------
  !> operator: dot product
  REAL(wp) PURE FUNCTION t_point_mul(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_mul = a1%x*a2%x + a1%y*a2%y + a1%z*a2%z
  END FUNCTION t_point_mul


  ! --------------------------------------------------------------------
  !> operator: product of scalar and point
  TYPE(t_point) PURE FUNCTION t_point_mul_scalar(v,t)
    TYPE(t_point), INTENT(IN) :: t
    REAL(wp),      INTENT(IN) :: v
    t_point_mul_scalar = t_point(x=v*t%x, y=v*t%y, z=v*t%z, &
      &                   ps=t%ps, gindex=t%gindex)
  END FUNCTION t_point_mul_scalar


  ! --------------------------------------------------------------------
  !> operator: division of point by scalar
  TYPE(t_point) PURE FUNCTION t_point_div_scalar(t,v)
    TYPE(t_point), INTENT(IN) :: t
    REAL(wp),      INTENT(IN) :: v
    t_point_div_scalar = t_point(x=t%x/v, y=t%y/v, z=t%z/v, &
      &                   ps=t%ps, gindex=t%gindex)
  END FUNCTION t_point_div_scalar


  ! --------------------------------------------------------------------
  ! DEFINITION OF "t_spherical_cap" SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> constructor for t_spherical_cap type
  PURE FUNCTION spherical_cap(pp, rr)
    TYPE(t_spherical_cap) :: spherical_cap
    TYPE(t_point), INTENT(IN)  :: pp
    REAL(wp),      INTENT(IN)  :: rr
    TYPE(t_point) :: point
    point = pp/pp%norm2()
    spherical_cap = t_spherical_cap(point=point, sorting_direction=(point/(-1._wp)), radius=rr)
  END FUNCTION spherical_cap


  ! --------------------------------------------------------------------
  ! DEFINITION OF "t_point" SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> constructor for t_point type
  PURE FUNCTION point(ix,iy,iz,iindex,ips)
    TYPE(t_point) :: point
    REAL(wp), INTENT(IN) :: ix, iy, iz
    INTEGER,  INTENT(IN), OPTIONAL :: iindex
    REAL(wp), INTENT(IN), OPTIONAL :: ips
    INTEGER  :: index
    REAL(wp) :: ps
    index =  0
    ps    = -2._wp
    IF (PRESENT(iindex))  index=iindex
    IF (PRESENT(ips))     ps=ips
    point = t_point(x=ix,y=iy,z=iz, ps=ps,gindex=index)
  END FUNCTION point

  ! --------------------------------------------------------------------
  !> operator: comparison of points (wrt. sorting direction)
  LOGICAL PURE FUNCTION t_point_cmp(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_cmp = (a2%ps > a1%ps)
  END FUNCTION t_point_cmp


  ! --------------------------------------------------------------------
  !> operator: comparison of points (wrt. sorting direction)
  LOGICAL PURE FUNCTION t_point_equal(a1,a2)
    TYPE(t_point), INTENT(IN) :: a1,a2
    t_point_equal = (a2%ps == a1%ps)
  END FUNCTION t_point_equal


  ! --------------------------------------------------------------------
  !> compute Euclidean norm in R^3
  REAL(wp) PURE FUNCTION point_norm2(this)
    CLASS(t_point), INTENT(IN) :: this
    point_norm2 = SQRT(this%x*this%x + this%y*this%y + this%z*this%z)
  END FUNCTION point_norm2


  ! --------------------------------------------------------------------
  !> compute spherical distance
  REAL(wp) PURE FUNCTION point_spherical_dist(a1,a2)
    CLASS(t_point), INTENT(IN) :: a1
    TYPE(t_point),  INTENT(IN) :: a2
    point_spherical_dist = ACOS(a1%x*a2%x + a1%y*a2%y + a1%z*a2%z)
  END FUNCTION point_spherical_dist


  ! --------------------------------------------------------------------
  ! DEFINITION OF "t_triangle" SUBROUTINES
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> constructor for t_triangle type
  PURE FUNCTION triangle(ip1, ip2, ip3, ioedge, icomplete)
    TYPE(t_triangle) :: triangle 
    INTEGER, INTENT(IN)           :: ip1, ip2, ip3
    INTEGER, INTENT(IN), OPTIONAL :: ioedge
    INTEGER, INTENT(IN), OPTIONAL :: icomplete
    INTEGER :: oedge, complete

    oedge    = -1
    complete =  0
    IF (PRESENT(ioedge))     oedge=ioedge
    IF (PRESENT(icomplete))  complete=icomplete
    triangle = t_triangle(p=(/ip1,ip2,ip3/), oedge=oedge, complete=complete,  &
      &                   cc=point(-2._wp,-2._wp,-2._wp), r=0._wp, cap_distance=2._wp, rdiscard=-1._wp)
  END FUNCTION triangle


  ! --------------------------------------------------------------------
  !> operator: comparison of triangles (wrt. point indices)
  LOGICAL PURE FUNCTION t_triangle_cmp(a1,a2)
    TYPE(t_triangle), INTENT(IN) :: a1,a2
    t_triangle_cmp =  (a1%p(0) < a2%p(0))                        .OR.                       &
      &        ((a1%p(0) == a2%p(0)) .AND. (a1%p(1) < a2%p(1)))  .OR.                       &
      &        ((a1%p(0) == a2%p(0)) .AND. (a1%p(1) == a2%p(1)) .AND. (a1%p(2) < a2%p(2)))
  END FUNCTION t_triangle_cmp


  ! --------------------------------------------------------------------
  PURE FUNCTION triangle_edge(this, i)
    TYPE(t_edge) :: triangle_edge
    CLASS(t_triangle),  INTENT(IN) :: this
    INTEGER,            INTENT(IN) :: i
    triangle_edge = t_edge( p1=this%p(i), p2=this%p(MOD(i+1,3)) )
  END FUNCTION triangle_edge

  
  
  ! --------------------------------------------------------------------
  PURE SUBROUTINE triangle_compute_circumcenter(this, pxyz, subset)
    CLASS(t_triangle),     INTENT(INOUT) :: this
    TYPE(t_point_list),    INTENT(IN)    :: pxyz
    TYPE(t_spherical_cap), INTENT(IN)    :: subset
    ! local variables
    TYPE (t_point) :: a,b,cc0

    ASSOCIATE (p1 => pxyz%a(this%p(0)), p2 => pxyz%a(this%p(1)), p3 => pxyz%a(this%p(2)))
      ! compute circumcenter as (p2-p3) x (p1-p3) normalized:
      a = p2 - p3
      b = p1 - p3
      cc0 = t_point( x=(a%y*b%z - b%y*a%z), y=(a%z*b%x - b%z*a%x), z=(a%x*b%y - a%y*b%x), &
        &            ps=0._wp, gindex=0) 
      this%cc = cc0/cc0%norm2()
      ! compute square radius
      this%r     = this%cc*p1                   ! dot product
      this%r     = SQRT(1._wp - this%r*this%r)    ! Pythagoras
      this%cc%ps = this%cc * subset%sorting_direction
      this%cap_distance = subset%point%spherical_dist(this%cc) - this%r
    END ASSOCIATE
  END SUBROUTINE triangle_compute_circumcenter


  ! --------------------------------------------------------------------
  !> cycle triangle indices st. smallest point index is first.
  PURE SUBROUTINE triangle_normalize_indices(this, pts)
    CLASS(t_triangle),     INTENT(INOUT) :: this
    TYPE(t_point_list),    INTENT(IN)    :: pts
    ! local variables
    INTEGER :: tp(0:2), ishift, j

    tp(0:2) = pts%a(this%p(0:2))%gindex
    ishift = 0
    IF      ((tp(1) <= tp(0)) .AND. (tp(1) <= tp(2))) THEN
      ishift = 1
    ELSE IF ((tp(2) <= tp(0)) .AND. (tp(2) <= tp(1))) THEN
      ishift = 2
    END IF
    DO j=0,2
      this%p(j) = tp(MOD(j + ishift, 3))
    END DO
  END SUBROUTINE triangle_normalize_indices


  ! --------------------------------------------------------------------
  ! SUBROUTINES FOR POINT/TRIANGLE LISTS (SORTABLE)
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> exchanges elements a(i), a(j)
  PURE SUBROUTINE initialize_object_list(this)
    CLASS(t_object_list), INTENT(INOUT) :: this
    this%nentries = 0
  END SUBROUTINE initialize_object_list

  ! --------------------------------------------------------------------
  !> constructor for t_point_list type
  FUNCTION point_list(list)
    TYPE(t_point_list) :: point_list
    TYPE(t_point_list), INTENT(IN)  :: list
    CALL point_list%resize(list%nentries)
    point_list%a(0:(list%nentries-1)) = list%a(0:(list%nentries-1))
  END FUNCTION point_list

  ! --------------------------------------------------------------------
  !> constructor for t_triangulation type
  FUNCTION triangulation(list)
    TYPE(t_triangulation) :: triangulation
    TYPE(t_triangulation), INTENT(IN)  :: list
    CALL triangulation%resize(list%nentries)
    triangulation%a(0:(list%nentries-1)) = list%a(0:(list%nentries-1))
  END FUNCTION triangulation

  ! --------------------------------------------------------------------
  !> constructor for t_triangulation type
  FUNCTION triangulation_ptr()
    TYPE(t_triptr) :: triangulation_ptr
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":triangulation_ptr"
    INTEGER :: ierrstat
    ALLOCATE(triangulation_ptr%ptr, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    CALL triangulation_ptr%ptr%initialize()
  END FUNCTION triangulation_ptr


  ! --------------------------------------------------------------------
  !> exchanges elements a(i), a(j)
  PURE SUBROUTINE swap_point_list(this, i,j)
    CLASS(t_point_list), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: i,j
    TYPE(t_point) :: t
    t         = this%a(i)
    this%a(i) = this%a(j)
    this%a(j) = t
  END SUBROUTINE swap_point_list


  ! --------------------------------------------------------------------
  !> @return .TRUE. if a(i) < a(j)
  LOGICAL PURE FUNCTION cmp_point_list(this, i,j)
    CLASS(t_point_list), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i,j
    cmp_point_list = (this%a(i) < this%a(j))
  END FUNCTION cmp_point_list


  ! --------------------------------------------------------------------
  !> finalization subroutine
  SUBROUTINE destructor_point_list(this)
    CLASS(t_point_list), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":destructor_point_list"
    INTEGER :: ierrstat
    IF (ALLOCATED(this%a)) THEN
      DEALLOCATE(this%a, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    this%nentries = 0
  END SUBROUTINE destructor_point_list


  ! --------------------------------------------------------------------
  !> Reserve subroutine.
  !
  !  Requests that the vector capacity be at least enough to contain
  !  @p new_size elements. If new_size is greater than the current
  !  vector capacity, the function causes the array to reallocate its
  !  storage increasing its capacity to new_size.  In all other cases,
  !  the function call does not cause a reallocation and the capacity
  !  is not affected.
  SUBROUTINE reserve_point_list(this, new_size)
    CLASS(t_point_list), INTENT(INOUT) :: this
    INTEGER,             INTENT(IN)    :: new_size
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":reserve_point_list"
    INTEGER :: ierrstat
    TYPE(t_point), ALLOCATABLE :: tmp(:)

    ! allocate temporary storage
    IF (ALLOCATED(this%a)) THEN
      IF (SIZE(this%a) < new_size) THEN
        ALLOCATE(tmp(0:(this%nentries-1)), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        ! copy existing data and resize
        tmp(0:(this%nentries-1)) = this%a(0:(this%nentries-1))
        DEALLOCATE(this%a, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      END IF
    END IF
    IF (.NOT. ALLOCATED(this%a)) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ELSE IF (SIZE(this%a) < new_size) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF
    IF (ALLOCATED(tmp)) this%a(0:(this%nentries-1)) = tmp(0:(this%nentries-1))
    ! clean up
    IF (ALLOCATED(tmp)) THEN
      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE reserve_point_list


  ! --------------------------------------------------------------------
  !> Resize subroutine.
  !
  !  Resizes the container so that it contains n elements.  If n is
  !  smaller than the current container size, the content is reduced
  !  to its first n elements.  If n is greater than the current
  !  container size, the content is expanded by inserting at the end
  !  as many elements as needed to reach a size of n.  If n is also
  !  greater than the current container capacity, an automatic
  !  reallocation of the allocated storage space takes place.
  !
  !  This subroutine does not initialize new elements!
  SUBROUTINE resize_point_list(this, new_size)
    CLASS(t_point_list) :: this
    INTEGER,   INTENT(IN) :: new_size
    CALL this%reserve(new_size)
    this%nentries = new_size
  END SUBROUTINE resize_point_list


  ! --------------------------------------------------------------------
  !> "push_back" subroutine, appends element to end of list.
  SUBROUTINE push_back_point_list(this, element)
    CLASS(t_point_list) :: this
    TYPE(t_point), INTENT(IN) :: element
    IF (.NOT. ALLOCATED(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    ELSE IF (this%nentries == SIZE(this%a)) THEN  
      CALL this%reserve(this%nentries + 1)
    END IF
    this%a(this%nentries) = element
    this%nentries         = this%nentries + 1
  END SUBROUTINE push_back_point_list


  ! --------------------------------------------------------------------
  !> exchanges elements a(i), a(j)
  PURE SUBROUTINE swap_triangulation(this, i,j)
    CLASS(t_triangulation), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: i,j
    TYPE(t_triangle) :: t
    t         = this%a(i)
    this%a(i) = this%a(j)
    this%a(j) = t
  END SUBROUTINE swap_triangulation


  ! --------------------------------------------------------------------
  !> @return .TRUE. if a(i) < a(j)
  LOGICAL PURE FUNCTION cmp_triangulation(this, i,j)
    CLASS(t_triangulation), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i,j
    cmp_triangulation = (this%a(i) < this%a(j))
  END FUNCTION cmp_triangulation


  ! --------------------------------------------------------------------
  !> finalization subroutine
  SUBROUTINE destructor_triangulation(this)
    CLASS(t_triangulation), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":destructor_triangulation"
    INTEGER :: ierrstat
    IF (ALLOCATED(this%a)) THEN
      DEALLOCATE(this%a, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    this%nentries = 0
  END SUBROUTINE destructor_triangulation


  ! --------------------------------------------------------------------
  !> Reserve subroutine.
  !
  !  Requests that the vector capacity be at least enough to contain
  !  @p new_size elements. If new_size is greater than the current
  !  vector capacity, the function causes the array to reallocate its
  !  storage increasing its capacity to new_size.  In all other cases,
  !  the function call does not cause a reallocation and the capacity
  !  is not affected.
  SUBROUTINE reserve_triangulation(this, new_size)
    CLASS(t_triangulation), INTENT(INOUT) :: this
    INTEGER,                INTENT(IN)    :: new_size
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":reserve_triangulation"
    INTEGER :: ierrstat
    TYPE(t_triangle), ALLOCATABLE :: tmp(:)

    ! allocate temporary storage
    IF (ALLOCATED(this%a)) THEN
      IF (SIZE(this%a) < new_size) THEN
        ALLOCATE(tmp(0:(this%nentries-1)), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        ! copy existing data and resize
        tmp(0:(this%nentries-1)) = this%a(0:(this%nentries-1))
        DEALLOCATE(this%a, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      END IF
    END IF
    IF (.NOT. ALLOCATED(this%a)) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ELSE IF (SIZE(this%a) < new_size) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF
    IF (ALLOCATED(tmp)) this%a(0:(this%nentries-1)) = tmp(0:(this%nentries-1))
    ! clean up
    IF (ALLOCATED(tmp)) THEN
      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE reserve_triangulation


  ! --------------------------------------------------------------------
  !> Resize subroutine.
  !
  !  Resizes the container so that it contains n elements.  If n is
  !  smaller than the current container size, the content is reduced
  !  to its first n elements.  If n is greater than the current
  !  container size, the content is expanded by inserting at the end
  !  as many elements as needed to reach a size of n.  If n is also
  !  greater than the current container capacity, an automatic
  !  reallocation of the allocated storage space takes place.
  !
  !  This subroutine does not initialize new elements!
  SUBROUTINE resize_triangulation(this, new_size)
    CLASS(t_triangulation) :: this
    INTEGER,   INTENT(IN) :: new_size
    CALL this%reserve(new_size)
    this%nentries = new_size
  END SUBROUTINE resize_triangulation


  ! --------------------------------------------------------------------
  !> "push_back" subroutine, appends element to end of list.
  SUBROUTINE push_back_triangulation(this, element)
    CLASS(t_triangulation) :: this
    TYPE(t_triangle), INTENT(IN) :: element
    IF (.NOT. ALLOCATED(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    ELSE IF (this%nentries == SIZE(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    END IF

    this%a(this%nentries) = element
    this%nentries         = this%nentries + 1
  END SUBROUTINE push_back_triangulation


  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  ! 
  !  Ordering after the sorting process: smallest...largest.
  !
  !  The user provides two callback functions: "cmp"("<") and "swap"
  !  which emulate a "Fortran template".
  !
  RECURSIVE SUBROUTINE quicksort_sortable_list(this, in_l, in_r)
    CLASS(t_sortable_list) :: this
    INTEGER, INTENT(IN), OPTIONAL :: in_l,in_r     !< left, right partition indices
    ! local variables
    INTEGER :: i,j,m,l,r

    l = 0
    r = this%nentries-1
    IF (PRESENT(in_l))  l=in_l
    IF (PRESENT(in_r))  r=in_r

    IF (r>l) THEN
      ! median-of-three selection of partitioning element
      IF ((r-l) > 3) THEN 
        m = (l+r)/2
        IF (this%cmp(m,l)) CALL this%swap(l,m)
        IF (this%cmp(r,l)) THEN
          CALL this%swap(l,r)
        ELSE IF (this%cmp(m,r)) THEN
          CALL this%swap(r,m)
        END IF
      END IF
      i = l-1
      j = r
      LOOP : DO
        CNTLOOP1 : DO
          i = i+1
          IF (.NOT. this%cmp(i,r)) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2 : DO
          j = j-1
          IF (.NOT. this%cmp(r,j) .OR. (j==0)) EXIT CNTLOOP2
        END DO CNTLOOP2
        CALL this%swap(i,j)
        IF (j <= i) EXIT LOOP
      END DO LOOP
      CALL this%swap(i,j) ! undo extra exchange
      CALL this%swap(i,r) ! put partitioning element into position
      CALL this%quicksort(l,i-1)
      CALL this%quicksort(i+1,r)
    END IF
  END SUBROUTINE quicksort_sortable_list

#ifdef __SX__
  ! Wrapper for NEC's vecorized radix sort algorithm
  SUBROUTINE radixsort_point_list(this)
    CLASS(t_point_list), INTENT(INOUT) :: this
    TYPE(t_point) :: tmp(size(this%a))
    REAL(wp) :: ps_arr(size(this%a))
    INTEGER :: co_arr(size(this%a))
    INTEGER :: n, i

    n = size(this%a)

    DO i=1,n
      ps_arr(i) = this%a(i-1)%ps
      tmp(i) = this%a(i-1)
      co_arr(i) = i
    END DO

    CALL radixsort(ps_arr,co_arr)

    DO i=1,n
      this%a(i-1) = tmp(co_arr(i))
    END DO

  END SUBROUTINE radixsort_point_list
#endif

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  ! 
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_point_list(this, in_l, in_r)
    CLASS(t_point_list) :: this
    INTEGER, INTENT(IN), OPTIONAL :: in_l,in_r     !< left, right partition indices
    ! local variables
    INTEGER :: i,j,m,l,r

    l = 0
    r = this%nentries-1
    IF (PRESENT(in_l))  l=in_l
    IF (PRESENT(in_r))  r=in_r

    IF (r>l) THEN
      ! median-of-three selection of partitioning element
      IF ((r-l) > 3) THEN 
        m = (l+r)/2

        IF (this%a(m) < this%a(l)) CALL swap_point_list(this,l,m)
        IF (this%a(r) < this%a(l)) THEN
          CALL swap_point_list(this,l,r)
        ELSE IF (this%a(m) < this%a(r)) THEN
          CALL swap_point_list(this,r,m)
        END IF
      END IF
      i = l-1
      j = r
      LOOP : DO
        CNTLOOP1 : DO
          i = i+1
          IF (.NOT. (this%a(i) < this%a(r))) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2 : DO
          j = j-1
          IF (.NOT. (this%a(r) < this%a(j)) .OR. (j==0)) EXIT CNTLOOP2
        END DO CNTLOOP2
        CALL swap_point_list(this,i,j)
        IF (j <= i) EXIT LOOP
      END DO LOOP
      CALL swap_point_list(this,i,j) ! undo extra exchange
      CALL swap_point_list(this,i,r) ! put partitioning element into position
      CALL this%quicksort(l,i-1)
      CALL this%quicksort(i+1,r)
    END IF
  END SUBROUTINE quicksort_point_list


  ! --------------------------------------------------------------------
  !> Synchronizes a point list between different processors with MPI.
  !
  !  Duplicate points are NOT removed.
  !
  SUBROUTINE sync_point_list(this, sync_gindex)
    CLASS(t_point_list) :: this
    LOGICAL, INTENT(IN), OPTIONAL :: sync_gindex
    ! local variables
#if (!defined(NOMPI))
    CHARACTER(*), PARAMETER :: routine = modname//":sync_point_list"
    INTEGER                        :: ierr, mpi_t_point, local_nentries, global_nentries, mpi_comm,    &
      &                               oldtypes(2), blockcounts(2), ierrstat, nranks, i, irank, &
      &                               this_start, this_end
    INTEGER(MPI_ADDRESS_KIND)      :: offsets(2), typeLB, extent
    INTEGER, ALLOCATABLE           :: recv_count(:), recv_displs(:)
    TYPE(t_mpi_point), ALLOCATABLE :: tmp(:), recv_tmp(:)
    LOGICAL                        :: lsync_gindex

    lsync_gindex = .TRUE.
    IF (PRESENT(sync_gindex))  lsync_gindex = sync_gindex

    mpi_comm = p_comm_work

    ! get no. of available MPI ranks:
    CALL MPI_COMM_SIZE (mpi_comm, nranks, ierr)
    CALL MPI_COMM_RANK (mpi_comm, irank,  ierr)

    ! communicate max. number of points
    local_nentries = this%nentries
    ALLOCATE(recv_count(0:(nranks-1)), recv_displs(0:(nranks-1)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    CALL MPI_ALLGATHER(local_nentries, 1, MPI_INTEGER, recv_count, 1, MPI_INTEGER, mpi_comm, ierr)
    global_nentries = SUM(recv_count(:))

    ! create an MPI-sendable copy of the point list
    ALLOCATE(tmp(0:(local_nentries-1)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    tmp(0:(local_nentries-1))%x      = this%a(0:(local_nentries-1))%x
    tmp(0:(local_nentries-1))%y      = this%a(0:(local_nentries-1))%y
    tmp(0:(local_nentries-1))%z      = this%a(0:(local_nentries-1))%z
    tmp(0:(local_nentries-1))%ps     = this%a(0:(local_nentries-1))%ps
    tmp(0:(local_nentries-1))%gindex = this%a(0:(local_nentries-1))%gindex

    ! create an MPI type for a single point
    !
    ! setup description of the 4 REAL(wp) fields: x, y, z, ps
    offsets(1)     = 0_MPI_ADDRESS_KIND
    oldtypes(1)    = p_real_dp
    blockcounts(1) = 4 
    CALL MPI_TYPE_GET_EXTENT(p_real_dp, typeLB, extent, ierr) 
    offsets(2)     = 4*extent
    oldtypes(2)    = MPI_INTEGER
    blockcounts(2) = 2
    CALL MPI_TYPE_CREATE_STRUCT(2, blockcounts, offsets, oldtypes, mpi_t_point, ierr) 
    CALL MPI_TYPE_COMMIT(mpi_t_point, ierr) 

    ! perform gather operation
    ALLOCATE(recv_tmp(0:(global_nentries-1)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    recv_displs(0) = 0
    DO i=1,(nranks-1)
      recv_displs(i) = recv_displs(i-1) + recv_count(i-1)
    END DO
    CALL MPI_ALLGATHERV(tmp, local_nentries, mpi_t_point, recv_tmp, recv_count, recv_displs, &
      &                 mpi_t_point, mpi_comm, ierr)

    CALL MPI_TYPE_FREE(mpi_t_point, ierr) 

    CALL this%resize(global_nentries)
    this%a(0:(this%nentries-1))%x      = recv_tmp(0:(this%nentries-1))%x 
    this%a(0:(this%nentries-1))%y      = recv_tmp(0:(this%nentries-1))%y 
    this%a(0:(this%nentries-1))%z      = recv_tmp(0:(this%nentries-1))%z 
    this%a(0:(this%nentries-1))%ps     = recv_tmp(0:(this%nentries-1))%ps
    IF (lsync_gindex) THEN
      this%a(0:(this%nentries-1))%gindex = recv_tmp(0:(this%nentries-1))%gindex
    ELSE
      this%a(0:(this%nentries-1))%gindex = -1
      this_start = recv_displs(irank)
      this_end   = this_start+recv_count(irank)-1
      this%a(this_start:this_end)%gindex = recv_tmp(this_start:this_end)%gindex
    END IF

    ! clean up
    DEALLOCATE(tmp, recv_tmp, recv_count, recv_displs, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#endif
  END SUBROUTINE sync_point_list


  ! utility function to swap two elements
  ELEMENTAL SUBROUTINE heap_swap(heap, i, j)
    TYPE(t_min_heap),      INTENT(INOUT)   :: heap
    INTEGER,               INTENT(IN)      :: i,j
    TYPE(t_min_heap_node) :: tmp
    tmp = heap%harr(i)
    heap%harr(i) = heap%harr(j)
    heap%harr(j) = tmp
  END SUBROUTINE heap_swap


  ! to get index of left child of node at index i
  ELEMENTAL INTEGER FUNCTION heap_left(i) 
    INTEGER, INTENT(IN) :: i
    heap_left = (2*i + 1)
  END FUNCTION heap_left

 
  ! build a heap from a given array a[] of given size
  SUBROUTINE construct_heap(heap, a)
    TYPE(t_min_heap),      INTENT(INOUT)  :: heap
    TYPE(t_min_heap_node), INTENT(IN)     :: a(0:)
    INTEGER :: i
    
    heap%isize = SIZE(a)
    ALLOCATE(heap%harr(0:(heap%isize-1)))
    heap%harr(:) = a(:)
    DO i=(heap%isize - 1)/2,0,-1
      CALL heap_heapify(heap, i)
    END DO
  END SUBROUTINE construct_heap


  ! Recursive method to heapify a subtree with root at given index.
  ! This method assumes that the subtrees are already heapified.
  RECURSIVE SUBROUTINE heap_heapify(heap, i)
    TYPE(t_min_heap), INTENT(INOUT)   :: heap
    INTEGER,          INTENT(IN)      :: i
    INTEGER :: l,r,smallest

    l        = heap_left(i)
    smallest = i
    IF (l < heap%isize) THEN
      IF (heap%harr(l)%elt%p < heap%harr(i)%elt%p)           smallest = l
      r = l+1
      IF (r < heap%isize) THEN
        IF (heap%harr(r)%elt%p < heap%harr(smallest)%elt%p)  smallest = r
      END IF
    END IF
    IF (smallest /= i) THEN
      CALL heap_swap(heap, i, smallest)
      CALL heap_heapify(heap, smallest)
    END IF
  END SUBROUTINE heap_heapify


  ! This function takes an array of arrays as an argument where all
  ! arrays are assumed to be sorted and merges them together.
  !
  ! Time complexity:
  !  O(nk Logk) ,  where n := SIZE(arr) , k := size(isize)
  !
  ! Based on C++ program
  !
  ! http://www.geeksforgeeks.org/merge-k-sorted-arrays
  !
  SUBROUTINE kway_merge(arr, isize, output, count)
    TYPE(t_mpi_triangle),  INTENT(IN)     :: arr(0:)
    INTEGER,               INTENT(IN)     :: isize(0:)
    TYPE (t_min_heap_elt), INTENT(INOUT)  :: output(0:)
    INTEGER,               INTENT(OUT)    :: count
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::kway_merge"
    INTEGER                 :: istart(0:(SIZE(isize)-1))
    INTEGER                 :: k, i, j
    TYPE(t_min_heap_node)   :: heap_init(0:(SIZE(isize)-1)), root
    TYPE (t_min_heap)       :: heap
    TYPE (t_min_heap_elt)   :: last_elt

    k     = SIZE(isize)
    count = 0

    ! Create a min heap with k heap nodes. Every heap node has first
    ! element of an array.
    istart(0) =  0
    j         = -1
    DO i=0,(k-1)
      IF (isize(i)>0) THEN
        j = j + 1
        heap_init(j)%elt%p%p = arr(istart(i))%p  ! Store the first element
        heap_init(j)%j = 1                       ! index of next element to be stored from array
        heap_init(j)%i = i                       ! index of array
      END IF
      IF (i /= (k-1)) THEN
        istart(i+1) = istart(i) + isize(i)
      END IF
    END DO

    IF (j == -1)  RETURN ! skip routine if no input
    CALL construct_heap(heap, heap_init(0:j)) ! create the heap

    ! Now one by one get the minimum element from min heap and replace
    ! it with next element of its array
    IF (SIZE(arr) > 0) THEN
      DO j=0,(SIZE(arr)-1)
        ! get the minimum element and store it in output
        root          = heap%harr(0)
        IF ((j == 0) .OR. (.NOT. (root%elt%p == last_elt%p))) THEN
          output(count)%p = root%elt%p
          last_elt        = root%elt
          count = count + 1
        END IF
        
        ! Find the next element that will replace current root of heap.
        ! The next element belongs to same array as the current root.
        IF (root%j < isize(root%i)) THEN
          root%elt%p%p = arr(istart(root%i)+root%j)%p
          root%j = root%j + 1
        ELSE
          ! If root was the last element of its array
          root%elt%p%p(0) =  HUGE(1)
        END IF
        ! Replace root with next element of array
        heap%harr(0) = root
        CALL heap_heapify(heap, 0); 
      END DO
    END IF
    ! clean up
    DEALLOCATE(heap%harr)
  END SUBROUTINE kway_merge


  ! --------------------------------------------------------------------
  !> Synchronizes a triangulation between different processors with MPI.
  !
  !  Duplicate triangles are removed.
  !
  !  First, this subroutine first calls itself recursively with a
  !  split MPI communicator, followed by a second call with the
  !  communicator that contains all worker PEs. This way, the gather
  !  process is processed in two stages which reqires smaller MPI
  !  buffers.
  !
  RECURSIVE SUBROUTINE sync_triangulation(this, opt_mpi_comm)
    CLASS(t_triangulation) :: this
    INTEGER, INTENT(IN), OPTIONAL :: opt_mpi_comm
#if (!defined(NOMPI))
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":sync_triangulation"
    INTEGER,      PARAMETER :: root_pe      = 0
    LOGICAL,      PARAMETER :: print_timers = .FALSE.

    INTEGER                           :: ierr, mpi_t_triangle, local_nentries,   &
      &                                  global_nentries, mpi_comm, oldtypes(2), &
      &                                  blockcounts(2), ierrstat, nranks, i,    &
      &                                  count, ngroups, comm1, comm2,           &
      &                                  irank, irank2, color1, color2,          &
      &                                  alloc_size
    INTEGER(MPI_ADDRESS_KIND)         :: offsets(2)
    INTEGER, ALLOCATABLE              :: recv_count(:), recv_displs(:)
    TYPE(t_mpi_triangle), ALLOCATABLE :: tmp(:), recv_tmp(:)
    TYPE(t_min_heap_elt), ALLOCATABLE :: kway_merge_array_out(:)
#ifdef _OPENMP
    DOUBLE PRECISION                  :: time_s, toc
#endif

    mpi_comm = p_comm_work
    IF (PRESENT(opt_mpi_comm))  mpi_comm = opt_mpi_comm

    IF (mpi_comm /= MPI_COMM_NULL) THEN
      ! create an MPI type for a single triangle
      offsets(1)     = 0_MPI_ADDRESS_KIND
      oldtypes(1)    = MPI_INTEGER
      blockcounts(1) = 3
      CALL MPI_TYPE_CREATE_STRUCT(1, blockcounts, offsets, oldtypes, mpi_t_triangle, ierr) 
      CALL MPI_TYPE_COMMIT(mpi_t_triangle, ierr) 
    ELSE
      RETURN
    END IF

    ! get no. of available MPI ranks
    CALL MPI_COMM_SIZE (mpi_comm, nranks, ierr)
    ! determine this worker's rank
    CALL MPI_COMM_RANK (mpi_comm, irank, ierr)

    IF (.NOT. PRESENT(opt_mpi_comm)) THEN
      ! we form SQRT(nranks) different groups
      ngroups = INT( SQRT(REAL(nranks)) )
      ! split MPI communicator: communicator for gather stage 1
      color1 = INT(irank/ngroups)
      CALL MPI_COMM_SPLIT(mpi_comm, color1, irank, comm1, ierr)
      ! split MPI communicator: communicator for gather stage 2
      CALL MPI_COMM_RANK (comm1, irank2, ierr)
      color2 = MPI_UNDEFINED
      IF (MOD(irank2, ngroups) == root_pe)  color2 = 1
      CALL MPI_COMM_SPLIT(mpi_comm, color2, irank, comm2, ierr)
#ifdef _OPENMP
      time_s = omp_get_wtime()
#endif

      ! local sort on each PE
      IF (this%nentries > 0)  CALL this%quicksort()
#ifdef _OPENMP
      toc = omp_get_wtime() - time_s
      IF (print_timers .AND. (irank == 0))  WRITE (0,*) "sorting: ", toc
      time_s = omp_get_wtime()
#endif
      ! recursive call, gather stage 1
      CALL sync_triangulation(this, comm1)
#ifdef _OPENMP
      toc = omp_get_wtime() - time_s

      IF (print_timers .AND. (irank == 0))  WRITE (0,*) "sync, stage 1: ", toc
      time_s = omp_get_wtime()
#endif
      ! recursive call, gather stage 2
      CALL sync_triangulation(this, comm2)
#ifdef _OPENMP
      toc = omp_get_wtime() - time_s

      IF (print_timers .AND. (irank == 0))  WRITE (0,*) "sync, stage 2: ", toc
      time_s = omp_get_wtime()
#endif
      IF (irank /= 0) THEN
        local_nentries = 0
      ELSE
        local_nentries = this%nentries
      END IF
      alloc_size = local_nentries
      ! broadcast buffer size from PE "irank == 0"
      CALL MPI_BCAST(alloc_size, 1, MPI_INTEGER, 0, mpi_comm, ierr)
#ifdef _OPENMP
      toc = omp_get_wtime() - time_s
      IF (print_timers .AND. (irank == 0))  WRITE (0,*) "epilogue 1: ", toc
#endif
    ELSE
      local_nentries = this%nentries
      alloc_size     = local_nentries
    END IF

    ! create an MPI-sendable copy of the triangle list
    IF (alloc_size > 0) THEN
      ALLOCATE(tmp(0:(alloc_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      IF (local_nentries > 0) THEN
        tmp(0:(local_nentries-1))%p(0) = this%a(0:(local_nentries-1))%p(0)
        tmp(0:(local_nentries-1))%p(1) = this%a(0:(local_nentries-1))%p(1)
        tmp(0:(local_nentries-1))%p(2) = this%a(0:(local_nentries-1))%p(2)
      END IF
    ELSE
      ALLOCATE(tmp(0:0), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF

    IF (.NOT. PRESENT(opt_mpi_comm)) THEN
      ! broadcast result from PE "irank == 0"
#ifdef _OPENMP
      time_s = omp_get_wtime()
#endif
      CALL MPI_BCAST(tmp, alloc_size, mpi_t_triangle, 0, mpi_comm, ierr)
#ifdef _OPENMP
      toc = omp_get_wtime() - time_s
      IF (print_timers .AND. (irank == 0))  WRITE (0,*) "bcast: ", toc
      time_s = omp_get_wtime()
#endif
      IF (irank /= 0) THEN
        CALL this%resize(alloc_size)
        IF (alloc_size > 0) THEN
          this%a(0:(alloc_size-1))%p(0) = tmp(0:(alloc_size-1))%p(0)
          this%a(0:(alloc_size-1))%p(1) = tmp(0:(alloc_size-1))%p(1)
          this%a(0:(alloc_size-1))%p(2) = tmp(0:(alloc_size-1))%p(2)
        END IF
      END IF

      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#ifdef _OPENMP
      toc = omp_get_wtime() - time_s
      IF (print_timers .AND. (irank == 0))  WRITE (0,*) "epilogue 2: ", toc
#endif
    ELSE
      ! communicate max. number of points
      ALLOCATE(recv_count(0:(nranks-1)), recv_displs(0:(nranks-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      CALL MPI_GATHER(local_nentries, 1, MPI_INTEGER, recv_count, 1, MPI_INTEGER, root_pe, mpi_comm, ierr)
      IF (irank == root_pe) THEN
        global_nentries = SUM(recv_count(:)) ! gather a total of "global_nentries" entries
      
        ! perform gather operation
        IF (global_nentries > 0) THEN
          ALLOCATE(recv_tmp(0:(global_nentries-1)), STAT=ierrstat)
        ELSE
          ALLOCATE(recv_tmp(0:0), STAT=ierrstat)
        END IF
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        recv_displs(0) = 0
        DO i=1,(nranks-1)
          recv_displs(i) = recv_displs(i-1) + recv_count(i-1)
        END DO
      ELSE
        ALLOCATE(recv_tmp(0:1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      END IF
      CALL MPI_GATHERV(tmp, local_nentries, mpi_t_triangle, recv_tmp, recv_count, recv_displs, &
        &              mpi_t_triangle, root_pe, mpi_comm, ierr)

      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
               
      IF (irank == root_pe) THEN
        ! NEC_RP exclude processes without local data from kway_merge
        IF (local_nentries > 0) THEN
          IF (global_nentries > 0) THEN
            ALLOCATE(kway_merge_array_out(0:global_nentries-1), STAT=ierrstat)
          ELSE
            ALLOCATE(kway_merge_array_out(0:0), STAT=ierrstat)
          END IF
          IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      
          ! perform a k-way merging of the sorted arrays from the k MPI
          ! processes, merge this to "count" entries.
          CALL kway_merge(recv_tmp, recv_count, kway_merge_array_out, count)

          DEALLOCATE(recv_tmp, recv_count, recv_displs, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

          CALL this%resize(count)
          IF (this%nentries > 0) THEN
            this%a(0:(this%nentries-1))%p(0) = kway_merge_array_out(0:(this%nentries-1))%p%p(0)
            this%a(0:(this%nentries-1))%p(1) = kway_merge_array_out(0:(this%nentries-1))%p%p(1)
            this%a(0:(this%nentries-1))%p(2) = kway_merge_array_out(0:(this%nentries-1))%p%p(2)
          END IF
          DEALLOCATE(kway_merge_array_out, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
        END IF
      ELSE
        DEALLOCATE(recv_tmp, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      END IF
    END IF

    ! clean up
    CALL MPI_TYPE_FREE(mpi_t_triangle, ierr) 
#endif
  END SUBROUTINE sync_triangulation


  ! --------------------------------------------------------------------
  !> Synchronizes a triangulation between different processors with MPI.
  !
  !  Duplicate triangles are removed.
  !
  !  First, this subroutine first calls itself recursively with a
  !  split MPI communicator, followed by a second call with the
  !  communicator that contains all worker PEs. This way, the gather
  !  process is processed in two stages which reqires smaller MPI
  !  buffers.
  !
  SUBROUTINE bcast_triangulation(this, root, mpi_comm)
    CLASS(t_triangulation) :: this
    INTEGER, INTENT(IN) :: root, mpi_comm
#if (!defined(NOMPI))
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":bcast_triangulation"
    INTEGER                           :: ierr, irank, mpi_t_triangle, oldtypes(2), &
      &                                  blockcounts(2), alloc_size
    INTEGER(MPI_ADDRESS_KIND)         :: offsets(2)
    TYPE(t_mpi_triangle), ALLOCATABLE :: tmp(:)

    IF (mpi_comm == MPI_COMM_NULL)  RETURN

    ! determine this worker's rank
    CALL MPI_COMM_RANK (mpi_comm, irank, ierr)
    IF ((irank /= root) .AND. (this%nentries > 0)) THEN
      CALL finish(routine, "Triangulation broadcast failed, non-empty receive buffer!")
    END IF

    ! broadcast buffer size from root PE
    alloc_size = this%nentries
    CALL MPI_BCAST(alloc_size, 1, MPI_INTEGER, root, mpi_comm, ierr)

    ! create an MPI type for a single triangle
    offsets(1)     = 0_MPI_ADDRESS_KIND
    oldtypes(1)    = MPI_INTEGER
    blockcounts(1) = 3
    CALL MPI_TYPE_CREATE_STRUCT(1, blockcounts, offsets, oldtypes, mpi_t_triangle, ierr) 
    CALL MPI_TYPE_COMMIT(mpi_t_triangle, ierr) 

    ! create an MPI-sendable copy of the triangle list
    ALLOCATE(tmp(0:(alloc_size-1)), STAT=ierr)
    IF (ierr /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    IF ((irank == root) .AND. (alloc_size > 0)) THEN
      tmp(0:(alloc_size-1))%p(0) = this%a(0:(alloc_size-1))%p(0)
      tmp(0:(alloc_size-1))%p(1) = this%a(0:(alloc_size-1))%p(1)
      tmp(0:(alloc_size-1))%p(2) = this%a(0:(alloc_size-1))%p(2)
    END IF

    ! broadcast data
    CALL MPI_BCAST(tmp, alloc_size, mpi_t_triangle, root, mpi_comm, ierr)

    ! copy-back of received triangle data
    IF (irank /= root) THEN
      CALL this%resize(alloc_size)
      IF (this%nentries > 0) THEN
        this%a(0:(this%nentries-1))%p(0) = tmp(0:(this%nentries-1))%p(0)
        this%a(0:(this%nentries-1))%p(1) = tmp(0:(this%nentries-1))%p(1)
        this%a(0:(this%nentries-1))%p(2) = tmp(0:(this%nentries-1))%p(2)
      END IF
    END IF

    ! clean up
    CALL MPI_TYPE_FREE(mpi_t_triangle, ierr) 
    DEALLOCATE(tmp, STAT=ierr)
    IF (ierr /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#endif
  END SUBROUTINE bcast_triangulation


  ! --------------------------------------------------------------------
  !> Write out triangulation in VTK ASCII format.
  !
  SUBROUTINE write_vtk_triangulation(tri, filename, p, ldata)
    CLASS(t_triangulation), INTENT(IN) :: tri
    CHARACTER(LEN=*),       INTENT(IN) :: filename
    TYPE(t_point_list),     INTENT(IN) :: p
    LOGICAL,                INTENT(IN) :: ldata
    ! local variables
    INTEGER :: out_unit, i
    
    ! write triangles in vtk format
    out_unit=20
    OPEN (unit=out_unit,file=TRIM(filename),action="write",status="replace")
    WRITE (out_unit,'(a)') "# vtk DataFile Version 2.0"
    WRITE (out_unit,'(a)') "# delaunay example"
    WRITE (out_unit,'(a)') "ASCII"
    WRITE (out_unit,'(a)') "DATASET UNSTRUCTURED_GRID"
    WRITE (out_unit,'(a,i0,a)') "POINTS ", p%nentries, " float"
    DO i=0,(p%nentries-1)
      WRITE (out_unit,'(3f12.8)') p%a(i)%x, p%a(i)%y, p%a(i)%z
    END DO
    WRITE (out_unit,'(a,i0,a,i0)') "CELLS ", tri%nentries," ", 4*tri%nentries
    DO i=0,(tri%nentries-1)
      WRITE (out_unit,'(3(a,i0))') "3 ", tri%a(i)%p(0)," ",tri%a(i)%p(1)," ",tri%a(i)%p(2)
    END DO
    WRITE (out_unit,'(a,i0)') "CELL_TYPES ", tri%nentries
    DO i=0,(tri%nentries-1)
      WRITE (out_unit,'(a)') "5"
    END DO
    IF (ldata) THEN
      WRITE (out_unit,'(a,i0)') "CELL_DATA ", tri%nentries
      WRITE (out_unit,'(a)')    "SCALARS rdiscard float 1"
      WRITE (out_unit,'(a)')    "LOOKUP_TABLE default"
      DO i=0,(tri%nentries-1)
        WRITE (out_unit,'(f7.3)') tri%a(i)%rdiscard
      END DO
    END IF
    CLOSE(out_unit)
  END SUBROUTINE write_vtk_triangulation


  ! --------------------------------------------------------------------
  !> Write out triangulation in NetCDF format.
  !
  !  Note: If the output file already exists, then the data is
  !        appended.
  !
  SUBROUTINE write_netcdf_triangulation(tri, filename)
    CLASS(t_triangulation), INTENT(IN) :: tri
    CHARACTER(LEN=*),       INTENT(IN) :: filename
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":write_netcdf_triangulation"
    CHARACTER(LEN=*), PARAMETER :: dimname1 = "cell_delaunay"
    CHARACTER(LEN=*), PARAMETER :: dimname2 = "vert_delaunay"
    CHARACTER(LEN=*), PARAMETER :: varname  = "cc_delaunay"
    LOGICAL :: lexist, lreadonly
    INTEGER :: ncfileid, has_id, has_var, dimID1, dimID2, varID, i, ierrstat
    INTEGER, ALLOCATABLE :: tri_data(:,:)

    INQUIRE (FILE=TRIM(filename), EXIST=lexist)
    IF (lexist) THEN
      ! check if we can write to the file
      lreadonly = .NOT. util_file_is_writable(filename)
      IF (lreadonly)  CALL finish(routine, "NetCDF file exists and is read-only!")
      CALL nf(nf90_open(TRIM(filename), NF90_WRITE, ncfileid), routine)
      ! check if the dimension and variable already exist:
      has_id = nf90_inq_dimid(ncfileid, dimname1, dimID1)
      IF (has_id == nf90_noerr)  CALL finish(routine, "Dimension 1 seems to exist!")
      has_id = nf90_inq_dimid(ncfileid, dimname2, dimID2)
      IF (has_id == nf90_noerr)  CALL finish(routine, "Dimension 2 seems to exist!")
      has_var = nf90_inq_varid(ncfileid, varname, varID)
      IF (has_var == nf90_noerr)  CALL finish(routine, "Variable seems to exist!")
    ELSE
      ! create file from scratch
      CALL nf(nf90_create(TRIM(filename), NF90_CLOBBER, ncfileid), routine)
    END IF

    ! define necessary dimension
    CALL nf(nf90_redef(ncfileid), routine)
    CALL nf(nf90_def_dim(ncfileid, dimname1, tri%nentries, dimID1), routine)    
    CALL nf(nf90_def_dim(ncfileid, dimname2,            3, dimID2), routine)    
    CALL nf(nf90_def_var(ncfileid, varname, NF90_INT, (/ dimID1, dimID2 /), varID), routine)    
    ! leave define mode
    CALL nf(nf90_enddef(ncfileid), routine)
    
    ALLOCATE(tri_data(0:(tri%nentries-1), 3), stat=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish(routine, "ALLOCATE failed!")
    DO i=0,(tri%nentries-1)
      tri_data(i,:) = tri%a(i)%p(:)
    END DO
    CALL nf(nf90_put_var(ncfileid, varid, tri_data), routine)
    DEALLOCATE(tri_data, stat=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish(routine, "DEALLOCATE failed!")

    CALL nf(nf90_close(ncfileid), routine)  ! close NetCDF file
  END SUBROUTINE write_netcdf_triangulation


  ! --------------------------------------------------------------------
  !> Read triangulation from NetCDF format.
  !
  INTEGER FUNCTION read_netcdf_triangulation(tri, filename)
    CLASS(t_triangulation), INTENT(INOUT) :: tri
    CHARACTER(LEN=*),       INTENT(IN)    :: filename
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":read_netcdf_triangulation"
    CHARACTER(LEN=*), PARAMETER :: dimname1 = "cell_delaunay"
    CHARACTER(LEN=*), PARAMETER :: dimname2 = "vert_delaunay"
    CHARACTER(LEN=*), PARAMETER :: varname  = "cc_delaunay"
    INTEGER :: ncfileid, has_id1, has_id2, has_var, dimID1, dimID2, varID, dimlen, ierrstat, i
    LOGICAL :: lexist
    INTEGER, ALLOCATABLE :: tri_data(:,:)

    read_netcdf_triangulation = SUCCESS

    INQUIRE (FILE=TRIM(filename), EXIST=lexist)
    IF (.NOT. lexist)  CALL finish(routine, "File '"//TRIM(filename)//"' not found!")

    ! open file in read-only mode:
    CALL nf(nf90_open(TRIM(filename), NF90_NOWRITE, ncfileid), routine)

    ! check if the dimension and variable exist:
    has_id1 = nf90_inq_dimid(ncfileid, dimname1, dimID1)
    has_id2 = nf90_inq_dimid(ncfileid, dimname2, dimID2)
    has_var = nf90_inq_varid(ncfileid, varname, varID)
    IF ((has_id1 /= nf90_noerr) .OR. (has_id2 /= nf90_noerr) .OR. (has_var /= nf90_noerr)) THEN

      read_netcdf_triangulation = -1

    ELSE

      CALL nf(nf90_inquire_dimension(ncfileid, dimID1, len = dimlen), routine)
      CALL tri%resize(dimlen)

      ALLOCATE(tri_data(0:(tri%nentries-1), 3), stat=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish(routine, "ALLOCATE failed!")
      CALL nf(nf90_get_var(ncfileid, varid, tri_data), routine)
      DO i=0,(tri%nentries-1)
        ! note that we shift the point indices st. they are suitable for
        ! 0-based arrays.
        tri%a(i)%p(:) = tri_data(i,:) - 1
      END DO
      DEALLOCATE(tri_data, stat=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish(routine, "DEALLOCATE failed!")

    END IF
    CALL nf(nf90_close(ncfileid), routine)  ! close NetCDF file
  END FUNCTION read_netcdf_triangulation


  ! --------------------------------------------------------------------
  ! SUBROUTINES FOR SPHERICAL CAP LISTS
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !> finalization subroutine
  SUBROUTINE destructor_sphcap_list(this)
    CLASS(t_sphcap_list), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":destructor_sphcap_list"
    INTEGER :: ierrstat
    IF (ALLOCATED(this%a)) THEN
      DEALLOCATE(this%a, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    this%nentries = 0
  END SUBROUTINE destructor_sphcap_list


  ! --------------------------------------------------------------------
  !> Reserve subroutine.
  !
  !  Requests that the vector capacity be at least enough to contain
  !  @p new_size elements. If new_size is greater than the current
  !  vector capacity, the function causes the array to reallocate its
  !  storage increasing its capacity to new_size.  In all other cases,
  !  the function call does not cause a reallocation and the capacity
  !  is not affected.
  SUBROUTINE reserve_sphcap_list(this, new_size)
    CLASS(t_sphcap_list) :: this
    INTEGER,   INTENT(IN) :: new_size
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":reserve_sphcap_list"
    INTEGER :: ierrstat
    TYPE(t_spherical_cap), ALLOCATABLE :: tmp(:)

    ! allocate temporary storage
    IF (ALLOCATED(this%a)) THEN
      IF (SIZE(this%a) < new_size) THEN
        ALLOCATE(tmp(0:(this%nentries-1)), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        ! copy existing data and resize
        tmp(0:(this%nentries-1)) = this%a(0:(this%nentries-1))
        DEALLOCATE(this%a, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      END IF
    END IF
    IF (.NOT. ALLOCATED(this%a)) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ELSE IF (SIZE(this%a) < new_size) THEN
      ALLOCATE(this%a(0:(new_size-1)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF
    IF (ALLOCATED(tmp)) this%a(0:(this%nentries-1)) = tmp(0:(this%nentries-1))
    ! clean up
    IF (ALLOCATED(tmp)) THEN
      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE reserve_sphcap_list


  ! --------------------------------------------------------------------
  !> Resize subroutine.
  !
  !  Resizes the container so that it contains n elements.  If n is
  !  smaller than the current container size, the content is reduced
  !  to its first n elements.  If n is greater than the current
  !  container size, the content is expanded by inserting at the end
  !  as many elements as needed to reach a size of n.  If n is also
  !  greater than the current container capacity, an automatic
  !  reallocation of the allocated storage space takes place.
  !
  !  This subroutine does not initialize new elements!
  SUBROUTINE resize_sphcap_list(this, new_size)
    CLASS(t_sphcap_list) :: this
    INTEGER,   INTENT(IN) :: new_size
    CALL this%reserve(new_size)
    this%nentries = new_size
  END SUBROUTINE resize_sphcap_list


  ! --------------------------------------------------------------------
  !> "push_back" subroutine, appends element to end of list.
  SUBROUTINE push_back_sphcap_list(this, element)
    CLASS(t_sphcap_list) :: this
    TYPE(t_spherical_cap), INTENT(IN) :: element
    IF (.NOT. ALLOCATED(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    ELSE IF (this%nentries == SIZE(this%a)) THEN
      CALL this%reserve(this%nentries + 1)
    END IF
    this%a(this%nentries) = element
    this%nentries         = this%nentries + 1
  END SUBROUTINE push_back_sphcap_list

END MODULE mo_delaunay_types
