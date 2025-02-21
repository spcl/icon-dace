! Contains the implementation of various mathematical algorithms
! used by the shallow water model.
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

! LL: xlc has trouble optimizing routines with implicit shaped parameters
! #ifdef __xlC__
! @process nohot
! ! @PROCESS NOOPTIMIZE
! #endif

MODULE mo_math_utilities
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi, pi_2, dbl_eps
  USE mo_impl_constants,      ONLY: SUCCESS, TORUS_MAX_LAT
  USE mo_exception,           ONLY: finish
  USE mo_grid_geometry_info,  ONLY: t_grid_geometry_info, planar_torus_geometry, sphere_geometry
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates, &
    &                               t_line, t_tangent_vectors
  USE mo_util_sort,           ONLY: quicksort



  IMPLICIT NONE

  PRIVATE


  PUBLIC :: cc2gc
  PUBLIC :: gc2cc
  PUBLIC :: cc2tv
  PUBLIC :: gvec2cvec
  PUBLIC :: cvec2gvec
  PUBLIC :: arc_length
  PUBLIC :: arc_length_v
  PUBLIC :: cartesian_to_geographical, geographical_to_cartesian
  PUBLIC :: arc_length_on_unitsphere
  PUBLIC :: cos_arc_length
  PUBLIC :: plane_torus_distance, plane_torus_closest_coordinates
  PUBLIC :: plane_angle_of_three_points
  PUBLIC :: sphere_cartesian_midpoint
  PUBLIC :: norma, normalize
  PUBLIC :: angle_of_vectors, sin_cc
  PUBLIC :: triangle_area
  PUBLIC :: circum_center
  PUBLIC :: bary_center
  PUBLIC :: middle
  PUBLIC :: sphere_tangent_coordinates
  PUBLIC :: normal_vector
  PUBLIC :: spherical_intersection, spherical_intersection2, spherical_intersection2_v
  PUBLIC :: rotate_z,rotate_x,rotate_y
  PUBLIC :: norma_of_vector_product
  PUBLIC :: project_point_to_plane

  PUBLIC :: disp_new
  PUBLIC :: disp_new_vect

  PUBLIC :: cc_dot_product
  PUBLIC :: cc_norm
  PUBLIC :: vector_product
  PUBLIC :: integral_over_triangle
  PUBLIC :: rotate_latlon
  PUBLIC :: rotate_latlon_vec

  PUBLIC :: gnomonic_proj
  PUBLIC :: orthogr_proj
  PUBLIC :: az_eqdist_proj
!   PUBLIC :: sphere_cell_mean_char_length
  PUBLIC :: ccw
  PUBLIC :: line_intersect
  PUBLIC :: lintersect
  PUBLIC :: tdma_solver
  PUBLIC :: tdma_solver_vec
  PUBLIC :: check_orientation

  !  vertical coordinates routines
  PUBLIC :: set_zlev

  ! subroutines for value sets
  PUBLIC :: t_value_set
  PUBLIC :: merge_values_into_set
  PUBLIC :: find_values_in_set
  PUBLIC :: deallocate_set

  PUBLIC :: OPERATOR(+)
  PUBLIC :: OPERATOR(-)
  PUBLIC :: OPERATOR(*)

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cartesian_coordinates_plus
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cartesian_coordinates_minus
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE cartesian_coordinates_prod
  END INTERFACE

  INTERFACE arc_length
    MODULE PROCEDURE arc_length_sphere
    MODULE PROCEDURE arc_length_generic
  END INTERFACE
  INTERFACE arc_length_v
    MODULE PROCEDURE arc_length_v_sphere
    MODULE PROCEDURE arc_length_v_generic
  END INTERFACE
  INTERFACE gc2cc
    MODULE PROCEDURE gc2cc_sphere
    MODULE PROCEDURE gc2cc_generic
  END INTERFACE
  INTERFACE cc2gc
    MODULE PROCEDURE cc2gc_sphere
    MODULE PROCEDURE cc2gc_generic
  END INTERFACE


  !> Derived type specifying a set of REAL(wp) values.
  !
  TYPE t_value_set
    INTEGER               :: nvalues                         !< no. of values in set
    LOGICAL               :: sort_smallest_first             !< Flag. If true, set is sorted smallest-to-largest
    REAL(wp), ALLOCATABLE :: values(:)                       !< sorted set of values.

    !> Everytime when new values are inserted into the set, the total
    !  list of values is re-sorted.
  END TYPE t_value_set

  REAL(wp), PARAMETER :: fcirc = 2.0_wp * pi
  REAL(wp), PARAMETER :: hcirc = fcirc * 0.5_wp, qcirc = fcirc * 0.25_wp
  ! constants for fixpoint<->float conversion of latitude/longitude angles
  REAL(wp), PARAMETER :: lat_spread = 0.25_wp * REAL(HUGE(1), wp)/qcirc, &
       lat_contract = 4.0_wp * qcirc/REAL(HUGE(1), wp), &
       lon_spread = 0.25_wp * REAL(HUGE(1), wp)/hcirc, &
       lon_contract = 4.0_wp * hcirc/REAL(HUGE(1), wp)

  PUBLIC :: normalized_lat, fxp_lat, flp_lat
  PUBLIC :: normalized_lon, fxp_lon, flp_lon

!-----------------------------------------------------------------------
! Basic geometry functions definitions








  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_math_utilities'




CONTAINS

  !-------------------------------------------------------------------------
  INTEGER FUNCTION check_orientation (lonc, lon, lat, n)
    INTEGER, INTENT(in) :: n
    REAL(wp), INTENT(in) :: lonc
    REAL(wp), INTENT(in) :: lon(n), lat(n)

    REAL(wp) :: lonl(n), latl(n)

    REAL(wp) :: area

    INTEGER :: i,j

    lonl(:) = lon(:)
    latl(:) = lat(:)

    DO i = 1, n
      lonl(i) = lonl(i) - lonc
      IF (lonl(i) < -pi) THEN
        lonl(i) =  pi+MOD(lonl(i), pi)
      ENDIF
      IF (lonl(i) >  pi) THEN
        lonl(i) = -pi+MOD(lonl(i), pi)
      ENDIF
    ENDDO

    area = 0.0_wp
    DO i = 1, n
      j = MOD(i,n)+1
      area = area+lonl(i)*latl(j)
      area = area-latl(i)*lonl(j)
    ENDDO

    IF (area >= 0.0_wp) THEN
      check_orientation = +1
    ELSE
      check_orientation = -1
    END IF

  END FUNCTION check_orientation
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian.
  !!
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian
  !! ones @f$(p\_cu, p\_cv, p\_cw)@f$ using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & -\sin p\_lat \cdot \cos p\_long \\  !! \cos p\_long & -\sin p\_lat \cdot \sin p\_long \\  !! 0 & -\cos p\_lat


  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix}
  !! @f}
  !!
  ELEMENTAL SUBROUTINE gvec2cvec (p_gu, p_gv, p_long, p_lat, p_cu, p_cv, p_cw, geometry_info)
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(in)  :: p_gu, p_gv     ! zonal and meridional vec. component
    REAL(wp), INTENT(in)  :: p_long, p_lat  ! geo. coord. of data point
    TYPE(t_grid_geometry_info), INTENT(in), OPTIONAL :: geometry_info

    REAL(wp), INTENT(out) :: p_cu, p_cv, p_cw            ! Cart. vector

    REAL(wp)              :: z_cln, z_sln, z_clt, z_slt  ! sin and cos of
    INTEGER               :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_math_utils:gvec2cvec'
    !-------------------------------------------------------------------------

    geometry_type = sphere_geometry
    IF (PRESENT(geometry_info)) &
      geometry_type = geometry_info%geometry_type

    SELECT CASE(geometry_type)

    CASE (planar_torus_geometry)
      p_cu  = p_gu
      p_cv  = p_gv
      p_cw  = 0._wp
    CASE (sphere_geometry)
      z_sln = SIN(p_long)
      z_cln = COS(p_long)
      z_slt = SIN(p_lat)
      z_clt = COS(p_lat)

      p_cu = z_sln * p_gu + z_slt * z_cln * p_gv
      p_cu = -1._wp * p_cu
      p_cv = z_cln * p_gu - z_slt * z_sln * p_gv
      p_cw = z_clt * p_gv
   CASE DEFAULT
  ! Commented by GZ because this destroys vectorization
  !    CALL finish(method_name, "Undefined geometry type")
      p_cu  = p_gu
      p_cv  = p_gv
      p_cw  = 0._wp
   END SELECT

  END SUBROUTINE gvec2cvec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$.
  !!
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$
  !! into zonal @f$p\_gu@f$ and meridional @f$g\_v@f$ velocity components
  !! using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & \cos p\_long & 0 \\  !! -\sin p\_lat \cdot \cos p\_long & -\sin p\_lat \cdot \sin p\_long & \cos p\_lat

  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix}
  !! @f}
  !!
  ELEMENTAL SUBROUTINE cvec2gvec (p_cu, p_cv, p_cw, p_long, p_lat, p_gu, p_gv, geometry_info)
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(in)  :: p_cu, p_cv, p_cw  ! Cart. vector
    REAL(wp), INTENT(in)  :: p_long, p_lat     ! geo. coord. of data point
    TYPE(t_grid_geometry_info), INTENT(in), OPTIONAL :: geometry_info

    REAL(wp), INTENT(out) :: p_gu, p_gv        ! zonal and meridional vec. comp.

    REAL(wp)              :: z_cln, z_clt, z_sln, z_slt  ! sin and cos of
    INTEGER               :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_math_utils:cvec2gvec'
    !-------------------------------------------------------------------------

    geometry_type = sphere_geometry
    IF (PRESENT(geometry_info)) &
      geometry_type = geometry_info%geometry_type

    SELECT CASE(geometry_type)

    CASE (planar_torus_geometry)
       p_gu = p_cu
       p_gv = p_cv
    CASE (sphere_geometry)
      ! p_long and p_lat
      z_sln = SIN(p_long)
      z_cln = COS(p_long)
      z_slt = SIN(p_lat)
      z_clt = COS(p_lat)

      p_gu = z_cln * p_cv - z_sln * p_cu
      p_gv = z_cln * p_cu + z_sln * p_cv
      p_gv = z_slt * p_gv
      p_gv = z_clt * p_cw - p_gv
    CASE DEFAULT
  ! Commented by GZ because this destroys vectorization
  !    CALL finish(method_name, "Undefined geometry type")
       p_gu = p_cu
       p_gv = p_cv
    END SELECT

  END SUBROUTINE cvec2gvec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts  vectors (in)  cartesian coordinate representation
  !! to tangent vectors.
  FUNCTION cc2tv(xx,position) result (tt)

    TYPE(t_cartesian_coordinates), INTENT(in) :: xx
    TYPE(t_geographical_coordinates),INTENT(in)  :: position
    TYPE(t_tangent_vectors):: tt

    REAL(wp) :: z_sinlo,z_coslo,z_sinlacoslo,z_sinlasinlo

    z_sinlo=SIN(position%lon)
    z_coslo=COS(position%lon)
    z_sinlacoslo=SIN(position%lat)*z_coslo
    z_sinlasinlo=SIN(position%lat)*z_sinlo
    tt%v1= -z_sinlo*xx%x(1)+z_coslo*xx%x(2)
    tt%v2= -z_sinlacoslo*xx%x(1)-z_sinlasinlo*xx%x(2) + COS(position%lat)*xx%x(3)

  END FUNCTION cc2tv
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE cartesian_to_geographical(cartesian, no_of_points, geo_coordinates, &
    & start_point, end_point )
    TYPE(t_cartesian_coordinates), POINTER :: cartesian(:)
    TYPE(t_geographical_coordinates), POINTER :: geo_coordinates(:)
    INTEGER, INTENT(in) :: no_of_points
    INTEGER, INTENT(in), OPTIONAL :: start_point, end_point

    INTEGER :: i, start_p, end_p

    IF (PRESENT(start_point)) THEN
      start_p = start_point
      end_p   = end_point
    ELSE
      start_p = 1
      end_p   = no_of_points
    ENDIF
    IF (.NOT. ASSOCIATED(geo_coordinates)) THEN
      ALLOCATE (geo_coordinates(start_p:end_p), stat=i)
      IF (i >0) THEN
        CALL finish ('cartesian_to_geographical', 'Problem in allocating local arrays')
      ENDIF
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(i)
    DO i=start_p, end_p
        geo_coordinates(i) = cc2gc(cartesian(i))
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE cartesian_to_geographical
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE geographical_to_cartesian(geo_coordinates, no_of_points, cartesian)
    TYPE(t_geographical_coordinates), POINTER :: geo_coordinates(:)
    TYPE(t_cartesian_coordinates), POINTER :: cartesian(:)
    INTEGER, INTENT(in) :: no_of_points

    INTEGER :: i
    REAL(wp) :: cos_lat
    REAL(wp) :: check

    IF (.NOT. ASSOCIATED(cartesian)) THEN
      ALLOCATE (cartesian(no_of_points), stat=i)
      IF (i >0) THEN
        CALL finish ('geographical_to_cartesian', 'Problem in allocating local arrays')
      ENDIF
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(i, cos_lat, check)
    DO i=1,no_of_points
      cos_lat           = COS(geo_coordinates(i)%lat)
      cartesian(i)%x(1) = COS(geo_coordinates(i)%lon) * cos_lat
      cartesian(i)%x(2) = SIN(geo_coordinates(i)%lon) * cos_lat
      cartesian(i)%x(3) = SIN(geo_coordinates(i)%lat)

      check = SQRT(DOT_PRODUCT(cartesian(i)%x,cartesian(i)%x))
      IF (ABS(check-1.0_wp) > 1.0e-8_wp) &
        & CALL finish ('geographicalToCartesian', 'Error in calculation')

      !    WRITE(*,*) i, " LonLat=", geoCoordinates(i)
      !    WRITE(*,*) i, " Cart=", cartesian(i)

    ENDDO ! i=1,noOfPoints
!$OMP END DO
!$OMP END PARALLEL

    RETURN

  END SUBROUTINE geographical_to_cartesian
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION norma (v) result(length)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v
    REAL(wp) :: length

    length = SQRT(DOT_PRODUCT(v%x,v%x))

  END FUNCTION norma
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  PURE FUNCTION normalize(v) result(n)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v
    TYPE(t_cartesian_coordinates) :: n

    n%x = v%x / SQRT(DOT_PRODUCT(v%x,v%x))

  END FUNCTION normalize
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION angle_of_vectors(v1,v2) result(angle)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v1,v2
    REAL(wp) :: angle

    REAL(wp) :: prod
   ! TYPE(cartesian_coordinates) :: v3

    prod = DOT_PRODUCT(v1%x,v2%x)/ &
      & (SQRT(DOT_PRODUCT(v1%x,v1%x)) * SQRT(DOT_PRODUCT(v2%x,v2%x)))
    angle = ATAN2(SQRT(1.0_wp - prod*prod), prod)

  END FUNCTION angle_of_vectors
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns a coordinate system (x,y),
  !! for the tanget plane at the spehere on the given point (lon,lat)
  PURE SUBROUTINE sphere_tangent_coordinates (geo_point,x,y)
    TYPE(t_geographical_coordinates), INTENT(in)  :: geo_point
    TYPE(t_cartesian_coordinates), INTENT(out) :: x,y

    REAL(wp) :: lon,lat
    !-------------------------------------------------------------------------
    lon = geo_point%lon
    lat = geo_point%lat
    x%x(1) = -SIN(lon)
    x%x(2) =  COS(lon)
    x%x(3) =  0._wp
    y%x(1) = -COS(lon) * SIN(lat)
    y%x(2) = -SIN(lon) * SIN(lat)
    y%x(3) =  COS(lat)

  END SUBROUTINE sphere_tangent_coordinates
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Computes normal vector to plane determined by the vectors x0,x1.
  ELEMENTAL FUNCTION normal_vector (x0, x1) result(x2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

!    x2%x = x2%x/ SQRT(x2%x(1)**2+x2%x(2)**2+x2%x(3)**2)
    x2%x = x2%x/ SQRT(DOT_PRODUCT(x2%x,x2%x))

  END FUNCTION normal_vector
  !-------------------------------------------------------------------------

  !--------------------------------------------------------------------
  !>
  ELEMENTAL FUNCTION middle(p1,p2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p1,p2
    TYPE(t_cartesian_coordinates) :: middle

    middle%x = (p1%x + p2%x) / 2.0_wp

  END FUNCTION middle
  !--------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  ELEMENTAL FUNCTION arc_length_sphere (p_x, p_y)  result (p_arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x, p_y  ! endpoints

    REAL(wp)            :: p_arc          ! length of geodesic arc

    REAL(wp)            :: z_lx,  z_ly    ! length of vector p_x and p_y
    REAL(wp)            :: z_cc           ! cos of angle between endpoints

    !-----------------------------------------------------------------------

    z_lx = SQRT(DOT_PRODUCT(p_x%x,p_x%x))
    z_ly = SQRT(DOT_PRODUCT(p_y%x,p_y%x))

    z_cc = DOT_PRODUCT(p_x%x, p_y%x)/(z_lx*z_ly)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    p_arc = ACOS(z_cc)

  END FUNCTION arc_length_sphere

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  FUNCTION arc_length_generic (p_x, p_y, geometry_info)  result (p_arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x, p_y  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    REAL(wp)            :: p_arc          ! length of geodesic arc

    INTEGER             :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_math_utils:arc_length'
    !-----------------------------------------------------------------------

    geometry_type = geometry_info%geometry_type

    SELECT CASE(geometry_type)

    CASE (planar_torus_geometry)
      !Assuming that the flat geometry is nothing but a small arc
      !over sphere. This assumption doesn't really affect any calculation
      !AD (20 Sept 2013) Now we use the planar distance instead. RBF scale
      !has been adjusted accordingly
      p_arc = plane_torus_distance(p_x%x,p_y%x,geometry_info)
    CASE (sphere_geometry)
      !
      p_arc = arc_length_sphere (p_x, p_y)
    CASE DEFAULT
      !
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END FUNCTION arc_length_generic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints x0,x1.
  ELEMENTAL FUNCTION arc_length_on_unitsphere (x0, x1) result(arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: arc
    REAL(wp) :: z_cc


    z_cc = DOT_PRODUCT(x0%x, x1%x)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    arc = ACOS(z_cc)

  END FUNCTION arc_length_on_unitsphere
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  !! Vectorizable version
  !!
  PURE FUNCTION arc_length_v_sphere (p_x, p_y)  result (p_arc)
    REAL(wp), INTENT(in) :: p_x(3), p_y(3)  ! endpoints

    REAL(wp)            :: p_arc          ! length of geodesic arc

    REAL(wp)            :: z_lx,  z_ly    ! length of vector p_x and p_y
    REAL(wp)            :: z_cc           ! cos of angle between endpoints

    !-----------------------------------------------------------------------

    z_lx = SQRT(DOT_PRODUCT(p_x,p_x))
    z_ly = SQRT(DOT_PRODUCT(p_y,p_y))

    z_cc = DOT_PRODUCT(p_x, p_y)/(z_lx*z_ly)

    ! in case we get numerically incorrect solutions

    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    p_arc = ACOS(z_cc)

  END FUNCTION arc_length_v_sphere

  !-------------------------------------------------------------------------
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  !! Vectorizable version
  FUNCTION arc_length_v_generic (p_x, p_y, geometry_info)  result (p_arc)
    REAL(wp), INTENT(in) :: p_x(3), p_y(3)  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    REAL(wp)            :: p_arc          ! length of geodesic arc

    INTEGER             :: geometry_type
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_math_utils:arc_length'
    !-----------------------------------------------------------------------

    geometry_type = geometry_info%geometry_type

    SELECT CASE(geometry_type)

    CASE (planar_torus_geometry)
      !Assuming that the flat geometry is nothing but a small arc
      !over sphere. This assumption doesn't really affect any calculation
      !AD (20 Sept 2013) Now we use the planar distance in these calculations
      !for TORUS. RBF scale has been adjusted accordingly
      p_arc = plane_torus_distance(p_x,p_y,geometry_info) !/ &
              !geometry_info%sphere_radius
    CASE (sphere_geometry)
      !
      p_arc = arc_length_v_sphere (p_x, p_y)
    CASE DEFAULT
      !
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END FUNCTION arc_length_v_generic
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !! Computes the cosine of the length of geodesic arc with endpoints x0,x1.
  ELEMENTAL FUNCTION cos_arc_length (x0, x1) result(carc)
   TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: carc
    REAL(wp) :: z_d0,  z_d1,  z_cc

    z_d0 = SQRT(DOT_PRODUCT(x0%x,x0%x))
    z_d1 = SQRT(DOT_PRODUCT(x1%x,x1%x))

    z_cc = DOT_PRODUCT(x0%x, x1%x)/(z_d0*z_d1)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    carc = z_cc

  END FUNCTION cos_arc_length
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! returns the torus modulo coordinates of v1 that are
  ! closest to v0
  FUNCTION plane_torus_closest_coordinates(v0, v1, &
    & geometry_info) result(new_v1_coord)

    REAL(wp), INTENT(in) :: v0(3), v1(3)
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    REAL(wp) :: length_of_torus, height_of_torus
    TYPE(t_cartesian_coordinates) :: new_v1_coord

    length_of_torus = geometry_info%domain_length
    height_of_torus = geometry_info%domain_height
    ! check the x coordinate
    IF ( ABS(v0(1) - v1(1)) >  length_of_torus * 0.5_wp) THEN
      ! we will wrap around + or -  length_of_torus
      IF (v0(1) > v1(1)) THEN
        new_v1_coord%x(1) = v1(1) + length_of_torus
      ELSE
        new_v1_coord%x(1) = v1(1) - length_of_torus
      ENDIF
    ELSE
      ! keep the same
      new_v1_coord%x(1) = v1(1)
    ENDIF

    ! check the y coordinate
    IF ( ABS(v0(2) - v1(2)) >  height_of_torus * 0.5_wp) THEN
      ! we will wrap around + or -  length_of_torus
      IF (v0(2) > v1(2)) THEN
        new_v1_coord%x(2) = v1(2) + height_of_torus
      ELSE
        new_v1_coord%x(2) = v1(2) - height_of_torus
      ENDIF
    ELSE
      ! keep the same
      new_v1_coord%x(2) = v1(2)
    ENDIF

    ! this should be zero
    new_v1_coord%x(3) = v1(3)

  END FUNCTION plane_torus_closest_coordinates
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  REAL(wp) FUNCTION plane_torus_distance(v0, v1, &
    & geometry_info)
    REAL(wp), INTENT(in) :: v0(3), v1(3)
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    TYPE(t_cartesian_coordinates) :: dv

    dv = plane_torus_closest_coordinates(v0, v1, geometry_info)
    dv%x = dv%x - v0
    plane_torus_distance = SQRT(DOT_PRODUCT(dv%x,dv%x))

  END FUNCTION plane_torus_distance
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION sphere_cartesian_midpoint (p1, p2) result(mid_point)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p1, p2
    TYPE(t_cartesian_coordinates) :: mid_point

      mid_point%x = (p1%x + p2%x)
      mid_point%x = mid_point%x / &
        & SQRT(DOT_PRODUCT(mid_point%x,mid_point%x))

  END FUNCTION sphere_cartesian_midpoint
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !! Computes the sin between two vectors x0, x1
  ELEMENTAL FUNCTION sin_cc (x0, x1)
    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: sin_cc

    TYPE(t_cartesian_coordinates) :: r

    r = vector_product(x0, x1)
    sin_cc = SQRT(DOT_PRODUCT(r%x,r%x)) / &
      & (SQRT(DOT_PRODUCT(x0%x,x0%x)) * SQRT(DOT_PRODUCT(x1%x,x1%x)))

  END FUNCTION sin_cc
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes spherical area of a triangle
  ELEMENTAL FUNCTION triangle_area (x0, x1, x2) result(area)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1, x2

    REAL(wp) :: area
    REAL(wp) :: z_s12, z_s23, z_s31, z_ca1, z_ca2, z_ca3, z_a1, z_a2, z_a3

    TYPE(t_cartesian_coordinates) :: u12, u23, u31

    ! This variant to calculate the area of a spherical triangle
    ! is more precise.

    !  Compute cross products Uij = Vi x Vj.
    u12 = vector_product (x0, x1)
    u23 = vector_product (x1, x2)
    u31 = vector_product (x2, x0)

    !  Normalize Uij to unit vectors.
    z_s12 = DOT_PRODUCT ( u12%x(1:3), u12%x(1:3) )
    z_s23 = DOT_PRODUCT ( u23%x(1:3), u23%x(1:3) )
    z_s31 = DOT_PRODUCT ( u31%x(1:3), u31%x(1:3) )

    !  Test for a degenerate triangle associated with collinear vertices.
    IF ( z_s12 == 0.0_wp .or. z_s23 == 0.0_wp  .or. z_s31 == 0.0_wp ) THEN
      area = 0.0_wp
      RETURN
    END IF

    z_s12 = SQRT(z_s12)
    z_s23 = SQRT(z_s23)
    z_s31 = SQRT(z_s31)

    u12%x(1:3) = u12%x(1:3)/z_s12
    u23%x(1:3) = u23%x(1:3)/z_s23
    u31%x(1:3) = u31%x(1:3)/z_s31

    !  Compute interior angles Ai as the dihedral angles between planes:
    !  CA1 = cos(A1) = -<U12,U31>
    !  CA2 = cos(A2) = -<U23,U12>
    !  CA3 = cos(A3) = -<U31,U23>
    z_ca1 = -u12%x(1)*u31%x(1)-u12%x(2)*u31%x(2)-u12%x(3)*u31%x(3)
    z_ca2 = -u23%x(1)*u12%x(1)-u23%x(2)*u12%x(2)-u23%x(3)*u12%x(3)
    z_ca3 = -u31%x(1)*u23%x(1)-u31%x(2)*u23%x(2)-u31%x(3)*u23%x(3)

    IF (z_ca1 < -1.0_wp) z_ca1 = -1.0_wp
    IF (z_ca1 >  1.0_wp) z_ca1 =  1.0_wp
    IF (z_ca2 < -1.0_wp) z_ca2 = -1.0_wp
    IF (z_ca2 >  1.0_wp) z_ca2 =  1.0_wp
    IF (z_ca3 < -1.0_wp) z_ca3 = -1.0_wp
    IF (z_ca3 >  1.0_wp) z_ca3 =  1.0_wp

    z_a1 = ACOS(z_ca1)
    z_a2 = ACOS(z_ca2)
    z_a3 = ACOS(z_ca3)

    !  Compute areas = z_a1 + z_a2 + z_a3 - pi.
    area = z_a1+z_a2+z_a3-pi

    IF ( area < 0.0_wp ) area = 0.0_wp

  END FUNCTION triangle_area
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Determines the circum_center of triangle with vertices v0,v1,v2.
  ELEMENTAL FUNCTION circum_center (v0, v1, v2) result(center)
    !> the coordinates of the three triangle vertices (unit vectors)
    !! in counter clockwise order.
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1, v2
    !> the coordinates of the circumcenter unless co-linear
    TYPE(t_cartesian_coordinates) :: center

    !  Local variables:
    TYPE(t_cartesian_coordinates) :: e1, e2 ! edges of the underlying planar triangle:
    ! v1-v0 ands v2-v0, respectively
    TYPE(t_cartesian_coordinates) :: cu     ! vector product of center:  e1 x e2
    REAL(wp) :: z_cnorm                     ! norm of cu

    !-----------------------------------------------------------------------

    e1%x = v1%x - v0%x
    e2%x = v2%x - v0%x

    ! compute cu = e1 x e2 and cnorm**2.
    cu = vector_product (e1, e2)

    IF (DOT_PRODUCT(cu%x, v0%x) < 0.0_wp) THEN
      cu%x = -cu%x
    END IF

    z_cnorm = SQRT(DOT_PRODUCT(cu%x,cu%x))
    center%x = cu%x/z_cnorm

  END FUNCTION circum_center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Determines the bary_center of triangle with vertices v0,v1,v2.
  !!
  !! @par Remarks
  !! currently not used in the code, just for testing purposes
  ELEMENTAL FUNCTION bary_center (v0, v1, v2) result(center)
    !> the coordinates of the three triangle vertices (unit vectors)
    !! in counter clockwise order.
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1, v2
    !> the coordinates of the bary_center, unless co-linear
    TYPE(t_cartesian_coordinates) :: center

    !  Local variables:
    TYPE(t_cartesian_coordinates) :: bc ! barycenter of the underlying planar
    ! triangle:
    REAL(wp) :: z_cnorm               ! norm of e1

    !-----------------------------------------------------------------------

    ! get barycenter of planar triangle

    bc%x = v0%x + v1%x + v2%x
    bc%x = bc%x/3._wp

    ! norm of barycenter
    z_cnorm = SQRT(DOT_PRODUCT(bc%x,bc%x))

    ! new center
    center%x = bc%x/z_cnorm

  END FUNCTION bary_center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Finds the intersection of two great circles.
  !!
  !! Cannot be made elemental, as long as the
  !! subroutine calls to message and finish are maintained.
  !!
  ELEMENTAL FUNCTION spherical_intersection (p0, p1, v0, v1) result(p)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0, p1
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1
    TYPE(t_cartesian_coordinates) :: p

    TYPE(t_cartesian_coordinates) :: en, cn

    REAL(wp) :: z_pn

    ! normal to plane where edge lies

    en = vector_product (v0, v1)

    ! normal tp plane where triangle centers lies

    cn = vector_product (p0, p1)

    p = vector_product (en, cn)

    ! rescale on sphere surface

    z_pn = SQRT(DOT_PRODUCT(p%x,p%x))

    p%x = p%x/z_pn

    IF (DOT_PRODUCT(p0%x, p%x) < 0.0_wp) THEN
      p%x = -p%x
    END IF

  END FUNCTION spherical_intersection
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! Projects a point to a plane defined by three points:
  ! to_v0, to_v1, to_v2
  ELEMENTAL FUNCTION project_point_to_plane (in_point, to_v0, to_v1, to_v2) result(projected_point)
    TYPE(t_cartesian_coordinates), INTENT(in) :: in_point, to_v0, to_v1, to_v2
    TYPE(t_cartesian_coordinates) :: projected_point

    TYPE(t_cartesian_coordinates) :: in_p, v1, v2, norm_to_plane

    ! translate everybody to the center
    v1%x = to_v1%x - to_v0%x
    v2%x = to_v2%x - to_v0%x
    in_p%x = in_point%x - to_v0%x

    norm_to_plane = vector_product (v1, v2)
    norm_to_plane%x=norm_to_plane%x/SQRT(DOT_PRODUCT(norm_to_plane%x,norm_to_plane%x))
    norm_to_plane%x = norm_to_plane%x * &
      & DOT_PRODUCT(norm_to_plane%x, (v1%x - in_p%x))

    projected_point%x = in_p%x + norm_to_plane%x + to_v0%x

  END FUNCTION project_point_to_plane
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION plane_angle_of_three_points(p1, p2, p3) result(angle)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p1,p2,p3

    REAL(wp) :: prod,angle
    TYPE(t_cartesian_coordinates) :: v1,v2

    v1%x = p1%x - p2%x
    v2%x = p3%x - p2%x
    prod = DOT_PRODUCT(v1%x,v2%x)/ (norma(v1) * norma(v2))
    angle = ACOS(prod)

  END FUNCTION plane_angle_of_three_points
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Finds the intersection of two great circles.
  !!
  !! Cannot be made elemental, as long as the
  !! subroutine calls to message and finish are maintained.
  !!
  ELEMENTAL FUNCTION spherical_intersection2 (p0, p1, v0, v1) result(p)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0, p1
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1
    TYPE(t_cartesian_coordinates) :: p

    TYPE(t_geographical_coordinates) :: g1,g2
    TYPE(t_cartesian_coordinates) :: en, cn

    !!$    REAL(wp), PARAMETER :: voronoi_tolerance = 1.0e-9_wp

    REAL(wp) :: z_pn

    ! normal to plane where edge lies
    !  en = vector_product (v0, v1)

    g1=cc2gc(v0)
    g2=cc2gc(v1)

    en = vecpro (g1, g2)

    ! normal tp plane where triangle centers lies
    ! cn = vector_product (p0, p1)


    g1=cc2gc(p0)
    g2=cc2gc(p1)

    cn = vecpro (g1, g2)

    ! check Voronoi/Delaunay properties
    !!$    IF (ABS(DOT_PRODUCT(en%x, cn%x)) > voronoi_tolerance) THEN
    !!$      WRITE (message_text,'(a,e12.6)') &
    !!$           'Orthogonality = ',                &
    !!$           ABS(DOT_PRODUCT(en%x, cn%x))
    !!$      CALL message ('inter_sect', TRIM(message_text))
    !!$      CALL finish ('inter_sect','Voronoi orthogonality violated.')
    !!$    END IF

    g1=cc2gc(en)
    g2=cc2gc(cn)

    !  p = vector_product (en, cn)
    p = vecpro (g1, g2)

    ! rescale on sphere surface
    z_pn = SQRT(DOT_PRODUCT(p%x,p%x))

    p%x = p%x/z_pn

    IF (DOT_PRODUCT(p0%x, p%x) < 0.0_wp) THEN
      p%x = -p%x
    END IF

  END FUNCTION spherical_intersection2

  !>
  !! Vectorized version of spherical_intersection2. Originally developed for the NEC,
  !! but also used for the Cray compiler because the original code (above) generates a segfault.
  !!
  SUBROUTINE spherical_intersection2_v (p0, p1, v0, v1, p, nlen)
    INTEGER, INTENT(in)                       :: nlen
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0(0:nlen-1), p1(0:nlen-1)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0(0:nlen-1), v1(0:nlen-1)
    TYPE(t_cartesian_coordinates), INTENT(out):: p(0:nlen-1)

    TYPE(t_geographical_coordinates) :: g1,g2
    TYPE(t_cartesian_coordinates) :: en, cn

    REAL(wp) :: z_pn
    INTEGER  :: ji

!$OMP PARALLEL DO PRIVATE(g1,g2,en,cn,z_pn)
    DO ji = 0, nlen-1
      g1=cc2gc(v0(ji))
      g2=cc2gc(v1(ji))

      en = vecpro (g1, g2)

      g1=cc2gc(p0(ji))
      g2=cc2gc(p1(ji))

      cn = vecpro (g1, g2)

      g1=cc2gc(en)
      g2=cc2gc(cn)

      !  p = vector_product (en, cn)
      p(ji) = vecpro (g1, g2)

      ! rescale on sphere surface
      z_pn = SQRT(DOT_PRODUCT(p(ji)%x,p(ji)%x))

      p(ji)%x = p(ji)%x/z_pn

      IF (DOT_PRODUCT(p0(ji)%x, p(ji)%x) < 0.0_wp) THEN
        p(ji)%x = -p(ji)%x
      END IF
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE spherical_intersection2_v

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the vector product of of 2 unit vectors in lon,lat
  ELEMENTAL FUNCTION vecpro(geog1,geog2) result (xx)

    TYPE(t_geographical_coordinates), INTENT(in) :: geog2, geog1

    TYPE(t_cartesian_coordinates)  :: xx
    REAL (wp)  :: long1, lat1,long2,lat2
    REAL (wp) :: sn1,sn2,sn3,sn4,cs1,cs2

    long1=geog1%lon
    lat1=geog1%lat
    long2=geog2%lon
    lat2=geog2%lat

    sn1=SIN(lat1-lat2)
    sn2=SIN(lat1+lat2)
    sn3=SIN(0.5_wp*(long1-long2))
    sn4=SIN(0.5_wp*(long1+long2))
    cs1=COS(0.5_wp*(long1-long2))
    cs2=COS(0.5_wp*(long1+long2))

    xx%x(1)=sn1*sn4*cs1-sn2*cs2*sn3
    xx%x(2)=sn1*cs2*cs1+sn2*sn4*sn3
    xx%x(3)=COS(lat1)*cos(lat2)*sin(long1-long2)

  END FUNCTION vecpro
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Computes the norma of the vector product of x0,x1.
  ELEMENTAL FUNCTION norma_of_vector_product (x0, x1) result(r)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: r

    REAL(wp) :: y1, y2, y3

    y1 = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    y2 = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    y3 = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

    r = SQRT(y1*y1 + y2*y2 + y3*y3)

  END FUNCTION norma_of_vector_product
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  Calculates coordinates of a grid point (longin, latin) in relation to
  !!  a rotated coordinate system with origin in (reflongin, reflatin)
  !!  with rotation matrix and Euler angles. <br>
  !!  The rotation is done by two basic rotations. First around unit vector in
  !!  noth pole direction (z-direction), second around new (y-direction).
  !!  The basic roatation matrices are:
  !!  @f{align*}{
  !!  R_z (\lambda) &=
  !!  \begin{pmatrix}
  !!   \cos\lambda & \sin\lambda & 0 \\  !!   -\sin\lambda & \cos\lambda & 0 \\  !!   0 & 0 & 1


  !!  \end{pmatrix}\\  !!  R_y (\varphi) &=

  !!  \begin{pmatrix}
  !!   \cos\varphi & 0 & \sin\varphi \\  !!   0 & 1 & 0 \\  !!   -\sin\varphi & 0 &\cos\varphi


  !!  \end{pmatrix}
  !!  @f}
  !!  Thus, the rotation matrix @f$R@f$ is
  !!  @f{align*}{
  !!  R(\lambda,\varphi) ~=~  R_y (\varphi) \cdot R_z (\lambda) &=
  !!  \begin{pmatrix}
  !!   \cos\lambda \cos\varphi & \sin\lambda \cos\varphi & \sin\varphi\\  !!   -\sin\lambda & \cos\lambda & 0\\  !!   -\cos\lambda \sin\varphi & -\sin\lambda \sin\varphi & \cos\varphi


  !!  \end{pmatrix}
  !!  @f}
  !!  Given the Cartesian coordinates of point
  !!  @f$(\lambda, \varphi) \mapsto (x,y,z)^T@f$
  !!  and applying the rotation matrix
  !!  @f{align*}{
  !!  R(\lambda,\varphi) \cdot (x,y,z)^T &= (1,0,0)^T
  !!  @f}
  !!  this point is the origin of the rotated coordinate system. All other points
  !!  are then rotated to this new system in the similar way by applying the
  !!  rotation matrix to their Cartesian coordinates.<br>
  !!  Algorithm:
  !!  <ol>
  !!  <li>  Given the geographical coordinates of a
  !!   point @f$(\lambda_{in}, \varphi_{in})@f$ and the origin of the rotated system
  !!  @f$(\lambda_{in}^{ref}, \varphi_{in}^{ref})@f$
  !!  <li>  Get Cartesian coordinates of @f$(\lambda_{in}, \varphi_{in})@f$ by
  !!  <i>gc2cc</i>
  !!  <li>  Apply roatation matrix @f$R(\lambda_{in}^{ref}, \varphi_{in}^{ref})@f$
  !!  <li>  Get rotated geographical coordinates @f$(\lambda_{dis}, \varphi_{dis})@f$
  !!  of new Cartesian ones by <i>cc2gc</i>
  !!  </ol>
  !!
  PURE SUBROUTINE disp_new(longin, latin, reflongin, reflatin, dislong, dislat)

    REAL(wp), INTENT(in)   :: longin, latin, reflongin, reflatin
    REAL(wp), INTENT(out)  :: dislong, dislat

    TYPE(t_geographical_coordinates) :: gpoint, gpoint_rot
    TYPE(t_cartesian_coordinates)    :: cpoint, row, cpoint_rot

    REAL(wp)               :: z_long, z_lat, z_reflong, z_reflat
    REAL(wp)               :: z_cos_a, z_cos_b, z_sin_a, z_sin_b
    REAL(wp)               :: z_twopi

    z_twopi = 2._wp * pi

    ! check if all angles are in the wright range
    z_long = longin
    IF (z_long < 0._wp) z_long = z_long + z_twopi

    z_reflong = reflongin
    IF (z_reflong < 0._wp) z_reflong = z_reflong + z_twopi

    z_lat    = latin
    z_reflat = reflatin

    ! get cartesian coordinates of (z_long, z_lat)

    gpoint%lon = z_long
    gpoint%lat = z_lat

    cpoint =  gc2cc(gpoint)

    ! calculate rotation matrix
    z_cos_a = COS(z_reflong)
    z_sin_a = SIN(z_reflong)
    z_cos_b = COS(z_reflat)
    z_sin_b = SIN(z_reflat)

    ! first row stored in a cartesian vector (just for using cc_dot_product)
    row%x(1) = z_cos_a * z_cos_b
    row%x(2) = z_sin_a * z_cos_b
    row%x(3) = z_sin_b

    ! first component of rotated cpoint
    cpoint_rot%x(1) = cc_dot_product (row, cpoint)

    ! second row stored in a cartesian vector (just for using cc_dot_product)
    row%x(1) = -1._wp * z_sin_a
    row%x(2) = z_cos_a
    row%x(3) = 0._wp

    ! second component of rotated cpoint
    cpoint_rot%x(2) = cc_dot_product (row, cpoint)

    ! third row stored in a cartesian vector (just for using cc_dot_product)
    row%x(1) = -1._wp * z_cos_a * z_sin_b
    row%x(2) = -1._wp * z_sin_a * z_sin_b
    row%x(3) = z_cos_b

    ! third component of rotated cpoint
    cpoint_rot%x(3) = cc_dot_product (row, cpoint)

    ! get geographic coordinates of cpoint_rot
    gpoint_rot = cc2gc(cpoint_rot)

    dislong = gpoint_rot%lon
    dislat  = gpoint_rot%lat

  END SUBROUTINE disp_new
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Compute the zonal and meridional components of a vector with respect.
  !!
  !! Compute the zonal and meridional components of a vector with respect
  !! to a rotated coordinate system. This function uses the same
  !! conventions and formalism as "disp_new". The algortithm is as
  !! follows:
  !! 1) convert the vector into Cartesian coordinates
  !! 2) apply the rotation to the vector in cartesian components
  !! 3) transform the rotated vector in geographical coordinates, using
  !! the new values of latitude and longitude.
  !!
  SUBROUTINE disp_new_vect(p_gu,p_gv,lon,lat,barlon,barlat, &
    & new_lon,new_lat,new_p_gu,new_p_gv)

    REAL(wp), INTENT(in) :: &
      & p_gu, p_gv,       & ! original lon-lat components
      & lon, lat,         & ! original lon and lat
      & barlon, barlat,   & ! new origin
      & new_lon, new_lat    ! point (lon,lat) in the new system (this
    ! value should be computed with disp_new)
    !
    REAL(wp), INTENT(out) :: &
      & new_p_gu, new_p_gv  ! new geographical components
    !
    REAL(wp) :: p_c(3), r1(3,3), r2(3,3), rot_p_c(3)
    !
    !--------------------------------------------------------------------

    ! 1) convert to Cartesian coordinates
    CALL gvec2cvec(p_gu,p_gv,lon,lat,p_c(1),p_c(2),p_c(3))

    ! 2) apply the rotation
    r1(1,:) = (/  COS(barlon) , SIN(barlon) , 0.0_wp /)
    r1(2,:) = (/ -SIN(barlon) , COS(barlon) , 0.0_wp /)
    r1(3,:) = (/     0.0_wp   ,    0.0_wp   , 1.0_wp /)

    r2(1,:) = (/  COS(barlat) , 0.0_wp , SIN(barlat) /)
    r2(2,:) = (/     0.0_wp   , 1.0_wp ,    0.0_wp   /)
    r2(3,:) = (/ -SIN(barlat) , 0.0_wp , COS(barlat) /)

    rot_p_c = MATMUL(r2,MATMUL(r1,p_c))

    ! 3) back to geographical coordinates
    CALL cvec2gvec(rot_p_c(1),rot_p_c(2),rot_p_c(3), &
      & new_lon,new_lat,new_p_gu,new_p_gv)

  END SUBROUTINE disp_new_vect
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts cartesian coordinates to geographical.
  !!
  ELEMENTAL FUNCTION cc2gc_sphere(p_x) result (p_pos)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x          ! Cart. coordinates
    TYPE(t_geographical_coordinates)          :: p_pos        ! geo. coordinates

    REAL(wp)                                :: z_x, z_y, z_z, z_r
    !-----------------------------------------------------------------------

    z_x = p_x%x(1)
    z_y = p_x%x(2)
    z_z = p_x%x(3)

    z_r = z_x * z_x + z_y * z_y
    z_r = SQRT(z_r)

    IF (ABS(z_r) < dbl_eps) THEN    ! one of the poles

      IF (z_z > 0.0_wp) THEN
        p_pos%lat = pi_2
      ELSE
        p_pos%lat = -1._wp * pi_2
      END IF
      p_pos%lon = 0._wp

    ELSE

      p_pos%lat = ATAN2 ( z_z, z_r)

      IF (ABS(z_x) < dbl_eps) THEN    ! z_x == 0 ?

        IF (z_y >= 0.0_wp) THEN
          p_pos%lon = pi_2
        ELSE
          p_pos%lon = -1._wp * pi_2
        END IF

      ELSE
        p_pos%lon = ATAN2( z_y, z_x)
      END IF

    END IF

  END FUNCTION cc2gc_sphere
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts cartesian coordinates to geographical.
  !!
  !! Implementation for plane torus geometry.
  !!
  !! The lon-lat parameterization of the torus is 
  !!    (lon,lat) = [0, 2*pi] x [-max_lat, max_lat]
  !! where max_lat := pi/18 = 10 degrees  (hard-coded in the torus grid generator)
  !!
  !! The Cartesian coordinates of the torus grid are: v = (x,y,0)
  !! where
  !!    (x,y) = [0,domain_length] x [0,domain_height]
  !! and the lengths are stored as global attributes in the grid file.
  !!
  FUNCTION cc2gc_plane_torus(p_x, geometry_info) RESULT (p_pos)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x          ! Cart. coordinates
    TYPE(t_grid_geometry_info),    INTENT(in) :: geometry_info
    TYPE(t_geographical_coordinates)          :: p_pos        ! geo. coordinates

    p_pos%lon = p_x%x(1) * 2._wp * pi          / geometry_info%domain_length
    p_pos%lat = p_x%x(2) * 2._wp*TORUS_MAX_LAT / geometry_info%domain_height - TORUS_MAX_LAT
  END FUNCTION cc2gc_plane_torus
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Converts cartesian coordinates to geographical.
  !!
  !! Generic implementation that handles sphere and other geometries.
  !!
  FUNCTION cc2gc_generic(p_x, geometry_info) RESULT (p_pos)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x          ! Cart. coordinates
    TYPE(t_grid_geometry_info),    INTENT(in) :: geometry_info
    TYPE(t_geographical_coordinates)          :: p_pos        ! geo. coordinates

    CHARACTER(LEN=*), PARAMETER :: method_name = modname//':cc2gc_generic'

    SELECT CASE(geometry_info%geometry_type)
    CASE (planar_torus_geometry)
      p_pos = cc2gc_plane_torus (p_x,geometry_info)
      !
    CASE (sphere_geometry)
      p_pos = cc2gc_sphere (p_x)
      !
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT
  END FUNCTION cc2gc_generic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Converts geographical to cartesian coordinates.
  !!
  ELEMENTAL FUNCTION gc2cc_sphere (p_pos)  result(p_x)
    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos     ! geo. coordinates
    TYPE(t_cartesian_coordinates)                :: p_x       ! Cart. coordinates
    REAL(wp) :: z_cln, z_sln, z_clt, z_slt

    z_sln = SIN(p_pos%lon)
    z_cln = COS(p_pos%lon)
    z_slt = SIN(p_pos%lat)
    z_clt = COS(p_pos%lat)

    p_x%x(1) = z_cln*z_clt
    p_x%x(2) = z_sln*z_clt
    p_x%x(3) = z_slt

  END FUNCTION gc2cc_sphere
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts geographical to cartesian coordinates.
  !!
  !! Implementation for plane torus geometry.
  !!
  !! The lon-lat parameterization of the torus is 
  !!    (lon,lat) = [0, 2*pi] x [-max_lat, max_lat]
  !! where max_lat := pi/18 = 10 degrees  (hard-coded in the torus grid generator)
  !!
  !! The Cartesian coordinates of the torus grid are: v = (x,y,0)
  !! where
  !!    (x,y) = [0,domain_length] x [0,domain_height]
  !! and the lengths are stored as global attributes in the grid file.
  !!
  FUNCTION gc2cc_plane_torus (p_pos, geometry_info)  RESULT(p_x)

    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos          ! geo. coordinates
    TYPE(t_grid_geometry_info),       INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates)                :: p_x            ! Cart. coordinates

    p_x%x(:) = (/ p_pos%lon * geometry_info%domain_length/(2._wp*pi),                              &
      &           (p_pos%lat + TORUS_MAX_LAT) * geometry_info%domain_height/(2._wp*TORUS_MAX_LAT), &
      &           0._wp /)
  END FUNCTION gc2cc_plane_torus
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts geographical to cartesian coordinates.
  !!
  !! Generic implementation that handles sphere and other geometries.
  !!
  FUNCTION gc2cc_generic (p_pos, geometry_info)  RESULT(p_x)

    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos          ! geo. coordinates
    TYPE(t_grid_geometry_info),       INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates)                :: p_x            ! Cart. coordinates

    CHARACTER(LEN=*), PARAMETER :: method_name = modname//':gc2cc_generic'

    SELECT CASE(geometry_info%geometry_type)
    CASE (planar_torus_geometry)
      p_x = gc2cc_plane_torus (p_pos,geometry_info)
      !
    CASE (sphere_geometry)
      p_x = gc2cc_sphere (p_pos)
      !
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT
  END FUNCTION gc2cc_generic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!  Calculates dot product of to cartesian coordinates.
  !!
  ELEMENTAL FUNCTION cc_dot_product(cc_x1, cc_x2)  result (p_prod)
    !$ACC ROUTINE SEQ
    TYPE(t_cartesian_coordinates), INTENT(in):: cc_x1, cc_x2 ! cart. coordinates

    REAL(wp) :: p_prod    ! scalar product of cart.coordinates

    p_prod = DOT_PRODUCT (cc_x1%x, cc_x2%x)

  END FUNCTION cc_dot_product
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!  Calculates 2 norm for vectors in cartesian coordinates.
  !!
  ELEMENTAL FUNCTION cc_norm(cc_x1)  result (norm)
    !$ACC ROUTINE SEQ
    TYPE(t_cartesian_coordinates), INTENT(in) :: cc_x1 ! cart. coordinates

    REAL(wp) :: norm

    norm = SQRT(cc_dot_product(cc_x1,cc_x1))

  END FUNCTION cc_norm
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Computes vector product of x0,x1.
  !!
  ELEMENTAL FUNCTION vector_product (x0, x1) result(x2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1

    TYPE(t_cartesian_coordinates) :: x2

    !-----------------------------------------------------------------------
    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

  END FUNCTION vector_product
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  PURE FUNCTION cartesian_coordinates_plus(x,y) result(z)

    TYPE(t_cartesian_coordinates) :: z
    TYPE(t_cartesian_coordinates), INTENT(in) :: x,y

    z%x = x%x + y%x

  END FUNCTION cartesian_coordinates_plus
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  PURE FUNCTION cartesian_coordinates_minus(x,y) result(z)

    TYPE(t_cartesian_coordinates) :: z
    TYPE(t_cartesian_coordinates), INTENT(in) :: x,y

    z%x = x%x - y%x

  END FUNCTION cartesian_coordinates_minus
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  PURE FUNCTION cartesian_coordinates_prod(a,x) result(y)

    TYPE(t_cartesian_coordinates) :: y
    REAL(wp), INTENT(in) :: a
    TYPE(t_cartesian_coordinates), INTENT(in) :: x

    y%x = a * x%x

  END FUNCTION cartesian_coordinates_prod
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !! Calculates integral @f$\int_{\tau} \Psi@f$ over spherical triangle @f$\tau@f$.
  !!
  !! Calculates integral @f$\int_{\tau} \Psi@f$ over spherical triangle @f$\tau@f$
  !! using second order quadrature rule
  !! @f{align*}{
  !!  I(\Psi, \tau) &= \frac{1}{6} \sum_{i=1}^3 \omega_i \Psi_i
  !! @f}
  !! See Boal and Sayas (2004), page 66 for details.
  !!
  FUNCTION integral_over_triangle (weight, values) result(p_integ)

    REAL(wp), INTENT(in) :: weight(3), values(3)

    REAL(wp)             :: p_integ

    !-----------------------------------------------------------------------

    p_integ = DOT_PRODUCT( weight, values)
    p_integ = p_integ / 6._wp

  END FUNCTION integral_over_triangle
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Rotate counterclockwise (seen from positive axis direction)
  !! vector x1 around axis z for phirot radians
  ELEMENTAL FUNCTION rotate_z (p_phirot,x1) result(x2)

    REAL(wp), INTENT(in) :: p_phirot
    TYPE(t_cartesian_coordinates), INTENT(in) :: x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(1) = COS(p_phirot)*x1%x(1) - SIN(p_phirot)*x1%x(2)
    x2%x(2) = SIN(p_phirot)*x1%x(1) + COS(p_phirot)*x1%x(2)
    x2%x(3) = x1%x(3)

  END FUNCTION rotate_z
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Rotate counterclockwise (seen from positive axis direction)
  !! vector x1 around axis y for phirot radians
  !!
  ELEMENTAL FUNCTION rotate_y (p_phirot,x1) result(x2)

    !-----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: p_phirot
    TYPE(t_cartesian_coordinates), INTENT(in) :: x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(1) =  COS(p_phirot)*x1%x(1) + SIN(p_phirot)*x1%x(3)
    x2%x(3) = -SIN(p_phirot)*x1%x(1) + COS(p_phirot)*x1%x(3)
    x2%x(2) =  x1%x(2)

  END FUNCTION rotate_y
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Rotate counterclockwise (seen from positive axis direction)
  !! vector x1 around axis x for phirot radians
  !!
  ELEMENTAL FUNCTION rotate_x (p_phirot,x1) result(x2)

    !-----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: p_phirot
    TYPE(t_cartesian_coordinates), INTENT(in) :: x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(2) =  COS(p_phirot)*x1%x(2) + SIN(p_phirot)*x1%x(3)
    x2%x(3) = -SIN(p_phirot)*x1%x(2) + COS(p_phirot)*x1%x(3)
    x2%x(1) =  x1%x(1)

  END FUNCTION rotate_x
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Rotates latitude and longitude for more accurate computation.
  !!
  !! Rotates latitude and longitude for more accurate computation
  !! of bilinear interpolation
  !!
  !! See the routine "mo_lonlat_grid::rotate_latlon_grid" for a
  !! detailed description of the transformation process.
  !!
  SUBROUTINE rotate_latlon( lat, lon, pollat, pollon )

    REAL(wp), INTENT(inout) :: lat, lon
    REAL(wp), INTENT(in)    :: pollat, pollon

    REAL(wp) :: rotlat, rotlon

    !-----------------------------------------------------------------------
    rotlat = ASIN(SIN(lat)*SIN(pollat) + COS(lat)*COS(pollat)*COS(lon-pollon))
    rotlon = ATAN2( COS(lat)*SIN(lon-pollon) , &
      & (COS(lat)*SIN(pollat)*COS(lon-pollon)- SIN(lat)*COS(pollat)) )

    lat = rotlat
    lon = rotlon

  END SUBROUTINE rotate_latlon
  !-------------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !! Provides rotation angle between coordinate systems
  !!
  !! Provides sin(delta) and cos(delta), the entries of the rotation matrix
  !! for vectors from the geographical to the rotated coordinate system,
  !! of which the pole is known.
  !!
  SUBROUTINE rotate_latlon_vec( lon, lat, pollon, pollat, sin_d, cos_d )
    !
    !
    REAL(wp), INTENT(in) :: lat, lon
    REAL(wp), INTENT(in) :: pollat, pollon
    REAL(wp), INTENT(out):: sin_d, cos_d

    REAL(wp) :: z_lamdiff, z_a, z_b, z_sq
    !-----------------------------------------------------------------------

    z_lamdiff = pollon - lon
    z_a       = COS(pollat)*SIN(z_lamdiff)
    z_b       = COS(lat)*SIN(pollat)-SIN(lat)*COS(pollat)*COS(z_lamdiff)
    z_sq      = SQRT(z_a*z_a + z_b*z_b)
    sin_d     = z_a/z_sq
    cos_d     = z_b/z_sq

  END SUBROUTINE rotate_latlon_vec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! gnomonic projection for unit sphere.
  !!
  !! Projects a point
  !! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
  !! The results are the local cartesian coordinates (x,y).
  !!
  !!
  !! Features:
  !! - azimuthal
  !! - neither conformal nor equal-area
  !! - no distortion at the center only
  !! - directions from the center are true
  !! - radial scale factor increases as distances increase from the center
  !! - scale factor in a direction perpendicular to a line radiating from
  !!   the center increases as distances increase from the center
  !! - all great circles are shown as straight lines
  !!
  SUBROUTINE gnomonic_proj( lon_c, lat_c, lon, lat, x, y )
    !
    ! !LITERATURE:
    ! Map Projections: A Working Manual, Snyder, 1987, p. 165
    !
    !
    REAL(wp), INTENT(in) :: lat_c, lon_c  ! center on tangent plane
    REAL(wp), INTENT(in) :: lat, lon      ! point to be projected
    REAL(wp), INTENT(out):: x, y          ! coordinates of projected point

    REAL(wp) :: zk   ! scale factor perpendicular to the radius from the
    ! center of the map
    REAL(wp) :: cosc ! cosine of the angular distance of the given point
    ! (lat,lon) from the center of projection

    !-----------------------------------------------------------------------

    cosc = SIN(lat_c)*SIN(lat) + COS(lat_c)*COS(lat)*COS(lon-lon_c)
    zk   = 1._wp/cosc

    x    = zk * COS(lat)*SIN(lon-lon_c)
    y    = zk * ( COS(lat_c)*SIN(lat) - SIN(lat_c)*COS(lat)*COS(lon-lon_c) )

  END SUBROUTINE gnomonic_proj
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! azimuthal equidistant projection for unit sphere.
  !!
  !! Projects a point
  !! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
  !! The results are the local cartesian coordinates (x,y).
  !!
  !!
  !! Features:
  !! - azimuthal
  !! - neither conformal nor equal-area
  !! - no distortion at the center only
  !! - directions from the center are true
  !! - radial scale is true
  !! - scale factor in a direction perpendicular to a line radiating from
  !!   the center linearly increases as distances increase from the center
  !!
  SUBROUTINE az_eqdist_proj( lon_c, lat_c, lon, lat, x, y )
    !
    ! !LITERATURE:
    ! Map Projections: A Working Manual, Snyder, 1987, p. 200
    !
    !
    REAL(wp), INTENT(in) :: lat_c, lon_c  ! center on tangent plane
    REAL(wp), INTENT(in) :: lat, lon      ! point to be projected
    REAL(wp), INTENT(out):: x, y          ! coordinates of projected point

    REAL(wp) :: zk   ! scale factor perpendicular to the radius from the
    ! center of the map
    REAL(wp) :: c    ! angular distance of the given point (lat,lon)
    ! from the center of projection

    !-----------------------------------------------------------------------

    c  = ACOS(SIN(lat_c)*SIN(lat) + COS(lat_c)*COS(lat)*COS(lon-lon_c))
    zk = c/SIN(c)

    x  = zk * COS(lat)*SIN(lon-lon_c)
    y  = zk * ( COS(lat_c)*SIN(lat) - SIN(lat_c)*COS(lat)*COS(lon-lon_c) )

  END SUBROUTINE az_eqdist_proj
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! orthographic projection for unit sphere.
  !!
  !! Projects a point
  !! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
  !! The results are the local cartesian coordinates (x,y).
  !!
  !! Features:
  !! - azimuthal
  !! - neither conformal nor equal-area
  !! - no distortion at the center only
  !! - directions from the center are true
  !! - radial scale factor decreases as distances increase from the center
  !! - scale in a direction perpendicular to a line radiating from the center is true
  !!
  SUBROUTINE orthogr_proj( lon_c, lat_c, lon, lat, x, y )
    !
    ! !LITERATURE:
    ! Map Projections: A Working Manual, Snyder, 1987,  p. 149
    !
    !
    REAL(wp), INTENT(in) :: lat_c, lon_c  ! center on tangent plane
    REAL(wp), INTENT(in) :: lat, lon      ! point to be projected
    REAL(wp), INTENT(out):: x, y          ! coordinates of projected point

    !-----------------------------------------------------------------------

    x = COS(lat)*SIN(lon-lon_c)
    y = COS(lat_c)*SIN(lat) - SIN(lat_c)*COS(lat)*COS(lon-lon_c)

  END SUBROUTINE orthogr_proj
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Calculates the domain mean of the characteristical
  !! length scale of the area. Needed for physics.
  !! Not used since the charecteristic lenntgh os now part of the grid_geometry_info. LL, MPI-M (2012-12)
!   SUBROUTINE sphere_cell_mean_char_length( total_number_of_cells, mean_charlen ) ! output
!
!     INTEGER , INTENT(in)  :: total_number_of_cells
!     REAL(wp), INTENT(out) :: mean_charlen
!
!     mean_charlen = SQRT (4._wp*pi*grid_sphere_radius**2 /REAL(total_number_of_cells,wp))
!
!   END SUBROUTINE sphere_cell_mean_char_length
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Checking turn when travelling along three points
  !!
  !! Given three points in a cartesian system, we want to know
  !! whether, in travelling from the first to the second to the third
  !! we turn counterclockwise or clockwise.
  !! Can be used to check whether two line segments intersect, or not.
  !!
  !! @par LITERATURE
  !! Sedgewick, R. (1988): Algorithms, 2nd edition, pp. 350
  !!
  ELEMENTAL FUNCTION ccw( p0, p1, p2 )
    !$ACC ROUTINE SEQ
    !

    IMPLICIT NONE

    TYPE(t_geographical_coordinates), INTENT(in) :: p0
    TYPE(t_geographical_coordinates), INTENT(in) :: p1
    TYPE(t_geographical_coordinates), INTENT(in) :: p2

    INTEGER :: ccw

    REAL(wp) :: dx1, dx2, dy1, dy2  ! segment lengths in x and y direction

    REAL(wp) :: dx1dy2, dy1dx2

    LOGICAL :: lccw

    !-----------------------------------------------------------------------
    ! segment lengths in x and y direction P0-->P1
    dx1 = p1%lon - p0%lon
    dy1 = p1%lat - p0%lat

    ! segment lengths in x and y direction P0-->P2
    dx2 = p2%lon - p0%lon
    dy2 = p2%lat - p0%lat

    dx1dy2 = dx1 * dy2
    dy1dx2 = dy1 * dx2


    ! in this case we turn counterclockwise
    ! dy2/dx2 > dy1/dx1
    lccw = dx1dy2 > dy1dx2

    ! we set ccw=1 when turning counterclockwise and -1 when
    ! turning clockwise.
    ! The case dy2/dx2 = dy1/dx1 is neglected. In this case
    ! we set ccw = -1
    ccw = MERGE(1, -1, lccw)

  END FUNCTION ccw
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Checks, whether two lines intersect
  !!
  !! Checks whether two lines intersect (2D cartesian geometry)
  !!
  !! @par LITERATURE
  !! Sedgewick, R. (1988): Algorithms, 2nd edition, pp. 351
  !!
  ELEMENTAL FUNCTION lintersect( line1, line2 )
    !$ACC ROUTINE SEQ

    TYPE(t_line), INTENT(in) :: line1
    TYPE(t_line), INTENT(in) :: line2

    INTEGER :: intersect1, intersect2

    LOGICAL :: lintersect

    !-----------------------------------------------------------------------
    intersect1 = ccw(line1%p1,line1%p2,line2%p1)      &
      & * ccw(line1%p1,line1%p2,line2%p2)

    intersect2 = ccw(line2%p1,line2%p2,line1%p1)      &
      & * ccw(line2%p1,line2%p2,line1%p2)

    lintersect = (intersect1 + intersect2) == -2

  END FUNCTION lintersect
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes intersection point of two lines in 2D
  !!
  !! Computes intersection point of two lines in 2D (cartesian geometry).
  !! Note that this function does not check, whether the two lines
  !! do intersect at all. So it is up to the user to make sure that
  !! the two lines indeed do intersect (e.g. by making use of the
  !! ccw-function)
  !!
  !! The two lines are given in the form:
  !! f1(x) = y_a1 + m1*(x-x_a1)
  !! f2(x) = y_b1 + m2*(x-x_b1)
  !! with
  !! m = (y_2-y_1)/(x_2-x_1)
  !!
  FUNCTION line_intersect( line1, line2 ) result(intersect)
    !$ACC ROUTINE SEQ

    TYPE(t_line), INTENT(in) :: line1
    TYPE(t_line), INTENT(in) :: line2

    REAL(wp) :: m1, m2          !< slopes
    REAL(wp) :: intersect(2)    ! coordinates of intersection point

    !-----------------------------------------------------------------------

    ! determine slopes of the two lines
    m1 = (line1%p2%lat - line1%p1%lat)/(line1%p2%lon - line1%p1%lon)
    m2 = (line2%p2%lat - line2%p1%lat)/(line2%p2%lon - line2%p1%lon)

    intersect(1) = (line2%p1%lat - line1%p1%lat + m1*line1%p1%lon - m2*line2%p1%lon) &
      & / (m1 - m2)

    intersect(2) = line1%p1%lat + m1*(intersect(1) - line1%p1%lon)

  END FUNCTION line_intersect


  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! TDMA tridiagonal matrix solver for a_i*x_(i-1) + b_i*x_i + c_i*x_(i+1) = d_i
  !!
  !!       a - sub-diagonal (means it is the diagonal below the main diagonal)
  !!       b - the main diagonal
  !!       c - sup-diagonal (means it is the diagonal above the main diagonal)
  !!       d - right part
  !!       x - the answer
  !!       n - number of equations
  SUBROUTINE tdma_solver(a,b,c,d,n,varout)

       !$ACC ROUTINE SEQ

       INTEGER, INTENT(in) :: n
       REAL(wp),DIMENSION(n),INTENT(in)  :: a,b,c,d
       REAL(wp),DIMENSION(n),INTENT(out) :: varout

       REAL(wp):: m, cp(n), dp(n)
       INTEGER :: i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize varout
         varout(n) = dp(n)
! solve for varout from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          varout(i) = dp(i)-cp(i)*varout(i+1)
        end do

  END SUBROUTINE tdma_solver

  !-------------------------------------------------------------------------
  !>
  !! TDMA tridiagonal matrix solver for a_i*x_(i-1) + b_i*x_i + c_i*x_(i+1) = d_i
  !!
  !!       a - sub-diagonal (means it is the diagonal below the main diagonal)
  !!       b - the main diagonal
  !!       c - sup-diagonal (means it is the diagonal above the main diagonal)
  !!       d - right part
  !!  varout - the answer (identical to x in description above)
  !!    slev - start level (top)
  !!    elev - end level (bottom)
  SUBROUTINE tdma_solver_vec(a,b,c,d,slev,elev,startidx,endidx,varout)
    INTEGER,   INTENT (IN) :: slev, elev
    INTEGER,   INTENT (IN) :: startidx, endidx
    REAL(wp),  INTENT (IN) :: a(:,:),b(:,:),c(:,:),d(:,:)
    REAL(wp),  INTENT(OUT) :: varout(:,:)
    !
    ! local
    REAL(wp):: m, cp(SIZE(a,1),SIZE(a,2)), dp(SIZE(a,1),SIZE(a,2))
    INTEGER :: i
    INTEGER :: jc

    !$ACC DATA CREATE(cp, dp) PRESENT(a, b, c, d) &
    !$ACC   PRESENT(varout) IF(i_am_accel_node .AND. acc_on)

    ! initialize c-prime and d-prime
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(slev, elev, startidx, endidx) ASYNC(1) IF(i_am_accel_node .AND. acc_on)
    !$ACC LOOP GANG VECTOR
    DO jc=startidx, endidx
      cp(jc,slev) = c(jc,slev)/b(jc,slev)
      dp(jc,slev) = d(jc,slev)/b(jc,slev)
    ENDDO
    !$ACC END PARALLEL
    ! solve for vectors c-prime and d-prime
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(slev, elev, startidx, endidx) ASYNC(1) IF(i_am_accel_node .AND. acc_on)
    !$ACC LOOP SEQ
!NEC$ outerloop_unroll(4)
    DO i = slev+1,elev
      !$ACC LOOP GANG VECTOR PRIVATE(m)
      DO jc=startidx, endidx
        m = 1._wp/(b(jc,i)-cp(jc,i-1)*a(jc,i))
        cp(jc,i) = c(jc,i) * m
        dp(jc,i) = (d(jc,i)-dp(jc,i-1)*a(jc,i)) * m
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    ! initialize varout
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node .AND. acc_on)
    varout(startidx:endidx,elev) = dp(startidx:endidx,elev)
    !$ACC END KERNELS
    ! solve for varout from the vectors c-prime and d-prime
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(slev, elev, startidx, endidx) ASYNC(1) IF(i_am_accel_node .AND. acc_on)
    !$ACC LOOP SEQ
!NEC$ outerloop_unroll(4)
    DO i = elev-1, slev, -1
      !$ACC LOOP GANG VECTOR
      DO jc=startidx, endidx
        varout(jc,i) = dp(jc,i)-cp(jc,i)*varout(jc,i+1)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA
  END SUBROUTINE tdma_solver_vec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !!! Helper functions for computing the vertical layer structure
  SUBROUTINE set_zlev(zlev_i, zlev_m, n_zlev, dzlev_m)
    INTEGER , INTENT(IN)  :: n_zlev
    REAL(wp), INTENT(OUT) :: zlev_i(n_zlev+1), zlev_m(n_zlev)
    REAL(wp), INTENT(IN)  :: dzlev_m(:)  ! namelist input of layer thickness

    INTEGER :: jk
    REAL(wp) :: accum


    accum = 0.0_wp
    ! zlev_i    : upper border surface of vertical cells
    ! zlev_m    : position of coordinate surfaces in meters below zero surface.
    DO jk = 1, n_zlev
      zlev_i(jk) = accum
      zlev_m(jk) = 0.5_wp * ( 2.0_wp * accum + dzlev_m(jk)  )
      accum = accum + dzlev_m(jk)
    END DO
    zlev_i(n_zlev+1) = accum

  END SUBROUTINE set_zlev

  !> normalize latitude to range [-pi/2,pi/2]
  ELEMENTAL FUNCTION normalized_lat(lat) RESULT(normalized)
    REAL(wp), INTENT(in) :: lat
    REAL(wp) :: normalized
    REAL(wp) :: norm_hcirc, qcirc_div, offset, sign_mult
    offset = SIGN(qcirc, lat)
    qcirc_div = AINT((lat + offset)/qcirc)
    norm_hcirc = MOD(lat + offset, hcirc) - offset
    sign_mult = -2.0_wp * ABS(MOD(AINT(qcirc_div*.5_wp), 2.0_wp)) + 1.0_wp
    normalized = sign_mult * norm_hcirc
  END FUNCTION normalized_lat

  !> convert latitude to normalized fixed-point integer
  ELEMENTAL FUNCTION fxp_lat(lat)
    REAL(wp), INTENT(in) :: lat
    INTEGER :: fxp_lat
    fxp_lat = NINT(lat_spread * normalized_lat(lat))
  END FUNCTION fxp_lat

  !> convert fixed-point latitude to floating-point value in range [-pi/2,pi/2]
  ELEMENTAL FUNCTION flp_lat(lat)
    INTEGER, INTENT(in) :: lat
    REAL(wp) :: flp_lat
    flp_lat = lat_contract * REAL(lat, wp)
  END FUNCTION flp_lat

  !> normalize latitude to range [-pi,pi]
  ELEMENTAL FUNCTION normalized_lon(lon) RESULT(normalized)
    REAL(wp), INTENT(in) :: lon
    REAL(wp) :: normalized,offset
    offset = SIGN(hcirc, lon)
    normalized = MOD(lon + offset, fcirc) - offset
  END FUNCTION normalized_lon

  !> convert longitude to normalized fixed-point integer
  ELEMENTAL FUNCTION fxp_lon(lon)
    REAL(wp), INTENT(in) :: lon
    INTEGER :: fxp_lon
    fxp_lon = NINT(lon_spread * normalized_lon(lon))
  END FUNCTION fxp_lon

  !> convert fixed-point longitude to floating-point value in range [-pi,pi]
  ELEMENTAL FUNCTION flp_lon(lon)
    INTEGER, INTENT(in) :: lon
    REAL(wp) :: flp_lon
    flp_lon = lon_contract * REAL(lon, wp)
  END FUNCTION flp_lon


  !-------------------------------------------------------------------------
  !> Merges a list of REAL(wp) numbers into another list, yielding the
  !  union set (without duplicates).
  !
  !  Initial implementation: F. Prill, DWD (2014-08-18)
  !
  SUBROUTINE merge_values_into_set(nin_values, in_values, value_set, opt_tol)
    INTEGER,           INTENT(IN)           :: nin_values     !< no. of input values
    REAL(wp),          INTENT(IN)           :: in_values(:)   !< list of input values
    TYPE(t_value_set), INTENT(INOUT)        :: value_set      !< result values (union set)
    REAL(wp),          INTENT(IN), OPTIONAL :: opt_tol        !< (optional:) tolerance for floating-point equality
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::merge_values_into_set'
    REAL(wp),         PARAMETER :: default_tol = 1.e-12_wp
    REAL(wp) :: tol
    INTEGER  :: i, k, n0, nn0, ierrstat
    INTEGER, ALLOCATABLE :: idx(:)

    ! set tolerance for floating-point equality test
    tol = default_tol
    IF (PRESENT(opt_tol)) tol = opt_tol
    ! make sure that we have enough space in our "t_value_set" object:
    CALL resize_set(value_set, nin_values)

    n0  = value_set%nvalues
    nn0 = n0+nin_values ! size of new set (with duplicates)
    ! first step: simply append the value to the end of the list:
    value_set%values      ((n0+1):nn0) = in_values(1:nin_values)
    ! second step: sort the whole set (including duplicates):
    CALL quicksort(value_set%values(1:nn0))
    ! revert levels largest-to-smallest if required
    IF (.NOT. value_set%sort_smallest_first) THEN
      value_set%values(1:nn0) = value_set%values(nn0:1:-1)
    END IF
    ! allocate temporary index buffer
    ALLOCATE(idx(nn0), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE of temporary failed!")
    ! third step: loop over the sorted list and mark duplicates
    k = 0
    DO i=1,nn0
      IF (k == 0) THEN
        k = k + 1
      ELSE
        IF (ABS(value_set%values(i-1) - value_set%values(i)) > tol)  k = k + 1
      END IF
      idx(i) = k
    END DO
    value_set%nvalues = k
    ! update value and index list (remove duplicates):
    DO i=1,nn0
      value_set%values(idx(i))  = value_set%values(i)
    END DO
    ! clean up
    DEALLOCATE(idx, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE of temporary failed!")
  END SUBROUTINE merge_values_into_set


  !-------------------------------------------------------------------------
  !> Compares a list of REAL(wp) numbers with the contents of a
  !  "t_value_set" data object. For a "value_set" with N values, this
  !  subroutine returns a LOGICAL array "out_flag" of size 1:N, where
  !  "out_flag(i) == .TRUE.", if the value "value_set%values(i)"
  !  matches one of the input numbers.
  !
  !  @todo For the sake of simplicity, this algorithm has quadratic
  !        complexity, though we could exploit the fact the the value
  !        set is sorted.
  !
  !  Initial implementation: F. Prill, DWD (2014-08-18)
  !
  SUBROUTINE find_values_in_set(nin_values, in_values, value_set, out_flag, opt_tol)
    INTEGER,           INTENT(IN)           :: nin_values     !< no. of input values
    REAL(wp),          INTENT(IN)           :: in_values(:)   !< list of input values
    TYPE(t_value_set), INTENT(IN)           :: value_set      !< set of values
    LOGICAL,           INTENT(OUT)          :: out_flag(:)    !<
    REAL(wp),          INTENT(IN), OPTIONAL :: opt_tol        !< (optional:) tolerance for floating-point equality
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::find_values_in_set'
    REAL(wp),         PARAMETER :: default_tol = 1.e-12_wp
    REAL(wp) :: tol
    INTEGER  :: i,j

    ! set tolerance for floating-point equality test
    tol = default_tol
    IF (PRESENT(opt_tol)) tol = opt_tol
    ! make sure that our result array has the required size:
    IF (SIZE(out_flag) < value_set%nvalues) THEN
      CALL finish(routine, "Wrong dimension of output argument!")
    END IF
    out_flag(:) = .FALSE.
    DO i=1,value_set%nvalues
      inner_loop: DO j=1,nin_values
        IF (ABS(value_set%values(i) - in_values(j)) <= tol) THEN
          out_flag(i) = .TRUE.
          EXIT inner_loop
        END IF
      END DO inner_loop
    END DO

  END SUBROUTINE find_values_in_set


  !> Resizes the data structures of a "t_value_set" derived type object.
  !
  !  Initial implementation: F. Prill, DWD (2014-08-18)
  !
  SUBROUTINE resize_set(value_set, nadditional_values)
    TYPE(t_value_set), INTENT(INOUT) :: value_set       !< result values (union set)
    INTEGER,           INTENT(IN)    :: nadditional_values
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::deallocate_set'
    INTEGER               :: ierrstat, oldsize, newsize
    REAL(wp), ALLOCATABLE :: tmp_rbuf(:)

    ! increase size of REAL(wp) buffer
    oldsize = 0
    IF (ALLOCATED(value_set%values)) THEN
      oldsize = SIZE(value_set%values)
      ALLOCATE(tmp_rbuf(oldsize), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      tmp_rbuf(1:oldsize) = value_set%values(1:oldsize)
      DEALLOCATE(value_set%values, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    ELSE
      value_set%nvalues = 0
    END IF
    newsize = oldsize + nadditional_values
    ALLOCATE(value_set%values(newsize), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    IF (ALLOCATED(tmp_rbuf)) THEN
      value_set%values(1:oldsize) = tmp_rbuf(1:oldsize)
      ! clean up
      DEALLOCATE(tmp_rbuf, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE resize_set


  !> Frees a "t_value_set" data type.
  !
  !  Initial implementation: F. Prill, DWD (2014-08-18)
  !
  SUBROUTINE deallocate_set(value_set)
    TYPE(t_value_set), INTENT(INOUT) :: value_set   !< result values (union set)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::deallocate_set'
    INTEGER :: ierrstat

    value_set%nvalues = 0
    DEALLOCATE(value_set%values, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE deallocate_set

END MODULE mo_math_utilities
