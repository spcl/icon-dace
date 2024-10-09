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

! Type definition and utility routines of lon-lat grids used in the model.

MODULE mo_lonlat_grid

  USE mo_kind,              ONLY: wp
  USE mo_mpi,               ONLY: p_bcast
  USE mo_exception,         ONLY: finish

  IMPLICIT NONE

  ! subroutines + functions:
  PUBLIC :: rotate_latlon_grid
  PUBLIC :: latlon_compute_area_weights
  PUBLIC :: compute_lonlat_specs
  PUBLIC :: compute_lonlat_blocking
  PUBLIC :: bcast_lonlat_specs
  PUBLIC :: OPERATOR(==)
  ! types and variables:
  PUBLIC :: t_lon_lat_grid
  PUBLIC :: threshold_delta_or_intvls

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_lonlat_grid'


  !---------------------------------------------------------------
  ! constants

  ! Threshold constant: If a "delta" value is provided, which is
  ! larger than this constant, then it is interpreted as the number of
  ! intervals for the lon-lat grid.
  REAL(wp), PARAMETER :: threshold_delta_or_intvls = 5._wp



  !> specification of (rotated) lon-lat grid
  !
  !  @note Variables of this type contain redundant information. The
  !        user must take care if all entries have been consistently
  !        initialized.
  !
  TYPE t_lon_lat_grid

    ! description of regular grid (in degrees) as in namelist
    INTEGER  :: reg_def_mode              ! 0: switch automatically between increment and no. of grid points
                                          ! 1: reg_lon/lat_def(2) specifies increment
                                          ! 2: reg_lon/lat_def(2) specifies no. of grid points

    REAL(wp) :: reg_lon_def(3)            ! start, increment OR grid points, end longitude in degrees
    REAL(wp) :: reg_lat_def(3)            ! start, increment OR grid points, end latitude in degrees
    REAL(wp) :: north_pole (2)            ! position of north pole (lon,lat) in degrees
                                           
    INTEGER  :: lon_dim                   ! Number of points in lon direction
    INTEGER  :: lat_dim                   ! Number of points in lat direction

    ! Alternative description of regular grid in rad,
    !   start_corner(lon,lat) + delta(lon,lat)*[0,1,2, ..., dim(lon,lat)-1]
    REAL(wp) ::                     &
      &  delta       (2),           &     ! lon-lat grid resolution,                unit:rad
      &  start_corner(2)                  ! south western corner of area (lon/lat), unit:rad

    ! Blocking information, computed from above values:
    INTEGER  :: total_dim                 ! total number of grid points 
    INTEGER  :: nblks, npromz             ! blocking info
    
  END TYPE t_lon_lat_grid

  INTERFACE OPERATOR (==)
    MODULE PROCEDURE lonlat_grid_compare
  END INTERFACE OPERATOR(==)

  
CONTAINS

  !---------------------------------------------------------------
  !> Computes rad-unit entries in lon-lat-grid definition.
  !
  !  This subroutine takes care of some special conventions:
  !
  !  1. For longitude values the last grid point is omitted if the end
  !     point matches the start point, e.g. for 0 and 360 degrees.
  !
  !  2. Instead of defining an increment it is also possible to
  !     prescribe the number of grid points, where it is expected that
  !     this value is larger than 5.0.
  !
  SUBROUTINE compute_lonlat_specs(lonlat_grid)
    TYPE (t_lon_lat_grid), INTENT(INOUT) :: lonlat_grid
    ! local variables
    REAL(wp),     PARAMETER :: ZERO_TOL = 1.e-15_wp
    CHARACTER(*), PARAMETER :: routine =  modname//"compute_lonlat_specs"
    REAL(wp) :: pi_180, lon_s, lon_e
    INTEGER  :: nintvls
    LOGICAL  :: lnpts_given(2), lskip_last_lon

    pi_180 = ATAN(1._wp)/45._wp
    lonlat_grid%delta(1)        = lonlat_grid%reg_lon_def(2) * pi_180
    lonlat_grid%delta(2)        = lonlat_grid%reg_lat_def(2) * pi_180
    lonlat_grid%start_corner(1) = lonlat_grid%reg_lon_def(1) * pi_180
    lonlat_grid%start_corner(2) = lonlat_grid%reg_lat_def(1) * pi_180

    ! Check for some special cases:
    !
    ! 1. --- check, if the longitude end value equals the longitude
    !        start value; skip the last point then
    lon_s          = lonlat_grid%reg_lon_def(1)
    lon_e          = lonlat_grid%reg_lon_def(3)
    IF (lon_s >= 360._wp) lon_s = lon_s - 360._wp
    IF (lon_e >= 360._wp) lon_e = lon_e - 360._wp
    lskip_last_lon = (ABS(lon_s - lon_e) < ZERO_TOL)
    ! 2. --- check if the "delta" value prescribes an interval size or
    !        the total *number* of intervals:
    SELECT CASE(lonlat_grid%reg_def_mode)
    CASE(0)
      ! 0: switch automatically between increment and no. of grid points
      lnpts_given(1) = (lonlat_grid%reg_lon_def(2) > threshold_delta_or_intvls)
      lnpts_given(2) = (lonlat_grid%reg_lat_def(2) > threshold_delta_or_intvls)
    CASE(1)
      ! 1: reg_lon/lat_def(2) specifies increment
      lnpts_given(1) = .FALSE.
      lnpts_given(2) = .FALSE.
    CASE(2)
      ! 2: reg_lon/lat_def(2) specifies no. of grid points
      lnpts_given(1) = .TRUE.
      lnpts_given(2) = .TRUE.
    END SELECT

    ! set longitude dimension
    IF (lnpts_given(1)) THEN
      ! no. of points was given:
      lonlat_grid%lon_dim  = NINT(lonlat_grid%reg_lon_def(2))
      IF (lskip_last_lon) THEN
        nintvls = lonlat_grid%lon_dim
      ELSE
        nintvls = lonlat_grid%lon_dim - 1
      END IF
      lonlat_grid%delta(1) = (lonlat_grid%reg_lon_def(3)-lonlat_grid%reg_lon_def(1))/nintvls * pi_180
    ELSE
      ! increment was given:
      IF (lonlat_grid%reg_lon_def(2) == 0._wp) THEN
        CALL finish(routine, "Invalid setting for reg_lon_def increment!")
      END IF
      lonlat_grid%lon_dim  = INT( (lonlat_grid%reg_lon_def(3)-lonlat_grid%reg_lon_def(1))/lonlat_grid%reg_lon_def(2) ) + 1
      IF (lskip_last_lon .AND. (lonlat_grid%lon_dim > 0)) THEN
         lonlat_grid%lon_dim = lonlat_grid%lon_dim - 1
      END IF
    END IF
    ! set latitude dimension
    IF (lnpts_given(2)) THEN
      ! no. of points was given:
      lonlat_grid%lat_dim  = NINT(lonlat_grid%reg_lat_def(2))
      IF (lonlat_grid%lat_dim == 1) THEN
         lonlat_grid%delta(2) = 0.0_wp
      ELSE
         lonlat_grid%delta(2) = (lonlat_grid%reg_lat_def(3)-lonlat_grid%reg_lat_def(1))/(lonlat_grid%lat_dim-1) * pi_180
      END IF
    ELSE
      ! increment was given:
      IF (lonlat_grid%reg_lat_def(2) == 0._wp) THEN
        CALL finish(routine, "Invalid setting for reg_lat_def increment!")
      END IF
      lonlat_grid%lat_dim  = INT( (lonlat_grid%reg_lat_def(3)-lonlat_grid%reg_lat_def(1))/lonlat_grid%reg_lat_def(2) ) + 1
    END IF
  END SUBROUTINE compute_lonlat_specs


  SUBROUTINE bcast_lonlat_specs(lonlat_grid, source, comm)
    TYPE (t_lon_lat_grid), INTENT(INOUT) :: lonlat_grid
    INTEGER,               INTENT(IN)    :: source, comm
    ! local variables
    REAL(wp) :: buf(12)
    INTEGER  :: int_buf(2)

    !-- broadcast REAL data:    
    buf = (/ lonlat_grid%delta(1),        lonlat_grid%delta(2),                                  &
      &      lonlat_grid%start_corner(1), lonlat_grid%start_corner(2),                           &
      &      lonlat_grid%reg_lon_def(1), lonlat_grid%reg_lon_def(2), lonlat_grid%reg_lon_def(3), &
      &      lonlat_grid%reg_lat_def(1), lonlat_grid%reg_lat_def(2), lonlat_grid%reg_lat_def(3), & 
      &      lonlat_grid%north_pole(1), lonlat_grid%north_pole(2) /)
    CALL p_bcast(buf, source, comm)
    lonlat_grid%delta(1)        = buf( 1)
    lonlat_grid%delta(2)        = buf( 2)
    lonlat_grid%start_corner(1) = buf( 3)
    lonlat_grid%start_corner(2) = buf( 4)
    lonlat_grid%reg_lon_def(1)  = buf( 5)
    lonlat_grid%reg_lon_def(2)  = buf( 6)
    lonlat_grid%reg_lon_def(3)  = buf( 7)
    lonlat_grid%reg_lat_def(1)  = buf( 8)
    lonlat_grid%reg_lat_def(2)  = buf( 9)
    lonlat_grid%reg_lat_def(3)  = buf(10)
    lonlat_grid%north_pole(1)   = buf(11)
    lonlat_grid%north_pole(2)   = buf(12)  
    !-- broadcast INTEGER data:
    int_buf = (/ lonlat_grid%lon_dim, lonlat_grid%lat_dim /)
    CALL p_bcast(int_buf, source, comm)
    lonlat_grid%lon_dim = int_buf(1)
    lonlat_grid%lat_dim = int_buf(2)
  END SUBROUTINE bcast_lonlat_specs


  !---------------------------------------------------------------
  !> Setup of nproma-blocking for lon-lat-grid
  !
  SUBROUTINE compute_lonlat_blocking(lonlat_grid, inproma)
    TYPE (t_lon_lat_grid), INTENT(INOUT) :: lonlat_grid
    INTEGER,               INTENT(IN)    :: inproma

    lonlat_grid%total_dim = lonlat_grid%lon_dim * lonlat_grid%lat_dim
    lonlat_grid%nblks     = (lonlat_grid%total_dim - 1)/inproma + 1
    lonlat_grid%npromz    = lonlat_grid%total_dim - (lonlat_grid%nblks-1)*inproma
  END SUBROUTINE compute_lonlat_blocking


  !-------------------------------------------------------------------------
  !> Rotates lon-lat grid
  !
  !  Rotates latitude and longitude for all grid points to standard grid
  !  coordinates.
  !
  !  Description of the rotated spherical coordinates system:
  !  [ cf. COSMO User Guide, Part I - Dynamis and Numerics, p.21 ]
  ! 
  !  "The origin of this new system is also located at the earth's
  !  centre, but the \tilde{Z}-axis is tilted against the Z-axis. By
  !  defining the \tilde{Z}-axis to point from the centre to a point
  !  P_N = (\lambda_g^N, \phi_g^N) in which \lambda_g^N is
  !  geographical longitude and \phi_g^N is geographical latitude of
  !  the point, the transformation is uniquely specified. P_N defines
  !  the north pole of the rotated coordinate system."
  !
  !  Transformation relations:
  !  [ cf. COSMO User Guide, Part I - Dynamis and Numerics, p.25 ] 
  !
  !  To transform the geographical longitude/latitude (\lambda_g, \phi_g)
  !  to the rotated horizontal coordinates (\lambda,\phi):
  !
  !  (1) :=  \phi_g    = arcsin( sin(\phi) \sin(\phi_g^N) + \cos(\phi) \cos(\lambda) \cos(\phi_g^N) )
  !
  !  (2) :=  \lambda_g = arctan( \frac{\cos(\phi) \sin(\phi)}
  !                        {\SIN(\phi_g^N) \COS(\phi) \COS(\lambda) - \sin(\phi) \cos(\phi_g^N)}
  !                      )  +  \lambda_g^N
  !
  !  Since we have atan(a/b) + c = atan((a + b*tan(c))/(b-a*tan(c)))
  !                              = atan((cos(c)*a + sin(c)*b)/(cos(c)*b-sin(c)*a))
  !  equation (2) is equivalent to the formulation below, which is also
  !  given in the COSMO database description, appendices A.1, A.2.
  !
  SUBROUTINE rotate_latlon_grid( lon_lat_grid, rotated_pts )
    
    TYPE (t_lon_lat_grid), INTENT(in)    :: lon_lat_grid
    REAL(wp),              INTENT(inout) :: rotated_pts(:,:,:)
    
    ! Local parameters
    REAL(wp), PARAMETER :: ZERO_TOL = 1.e-15_wp

    REAL(wp) :: sincos_pole(2,2),                   &  ! (lon/lat, sin/cos)
      &         sincos_lon(lon_lat_grid%lon_dim,2), &
      &         sincos_lat(lon_lat_grid%lat_dim,2), &
      &         pi_180, npole_rad(2), arg1, arg2,   &
      &         rlon(lon_lat_grid%lon_dim),         &
      &         rlat(lon_lat_grid%lat_dim)
    INTEGER  :: k, j
    LOGICAL  :: ltrivial_rotation

    pi_180 = ATAN(1._wp)/45._wp   

    ! check for "trivial rotation" (no.pole at +90,0):
    ltrivial_rotation  = ((ABS(90._wp - lon_lat_grid%north_pole(2)) < ZERO_TOL) .AND.  &
      &                    ABS( 0._wp - lon_lat_grid%north_pole(1)) < ZERO_TOL)

    ! compute the non-rotated lon/lat values
    DO k=1,lon_lat_grid%lon_dim
      rlon(k)         = lon_lat_grid%start_corner(1) + REAL(k-1,wp)*lon_lat_grid%delta(1)
      sincos_lon(k,1) = SIN(rlon(k))
      sincos_lon(k,2) = COS(rlon(k))
    END DO
    DO k=1,lon_lat_grid%lat_dim
      rlat(k)         = lon_lat_grid%start_corner(2) + REAL(k-1,wp)*lon_lat_grid%delta(2)
      sincos_lat(k,1) = SIN(rlat(k))
      sincos_lat(k,2) = COS(rlat(k))
    END DO

    ! special treatment for "trivial rotation" (no.pole at +90,0):
    IF (ltrivial_rotation) THEN

      DO j = 1, lon_lat_grid%lat_dim
        DO k = 1, lon_lat_grid%lon_dim
          rotated_pts(k,j,1) = rlon(k)
          rotated_pts(k,j,2) = rlat(j)
        END DO
      END DO

    ELSE

      ! convert north pole: degree -> rad:
      npole_rad(:) = lon_lat_grid%north_pole(:)*pi_180

      sincos_pole(:,1) = SIN(npole_rad(:))
      sincos_pole(:,2) = COS(npole_rad(:))  

      DO j = 1, lon_lat_grid%lat_dim
        DO k = 1, lon_lat_grid%lon_dim

          ! ATAN2(COS(phi)*SIN(lambda), SIN(poleY)*COS(phi)*COS(lambda) - SIN(phi)*COS(poleY)) + poleX

          arg1 = sincos_lat(j,2)*sincos_lon(k,1)
          arg2 = sincos_pole(2,1)*sincos_lat(j,2)*sincos_lon(k,2) - sincos_lat(j,1)*sincos_pole(2,2)

          IF ((ABS(arg1) > ZERO_TOL) .OR. (ABS(arg2) > ZERO_TOL)) THEN
            ! rotated_pts(k,j,1) = ATAN2( arg1, arg2 ) + npole_rad(1)
            rotated_pts(k,j,1) = ATAN2( -1._wp*sincos_pole(1,1)*arg2 - sincos_pole(1,2)*arg1 ,&
              &                         -1._wp*sincos_pole(1,2)*arg2 + sincos_pole(1,1)*arg1 )
          ELSE
            rotated_pts(k,j,1) = 0.0_wp ! ATAN2(0,0) is undefined, so we just have to set something
          ENDIF

          ! ASIN( SIN(phi)*SIN(poleY) + COS(phi)*COS(lambda)*COS(poleY) )
          rotated_pts(k,j,2) = &
            & ASIN( sincos_lat(j,1)*sincos_pole(2,1) + &
            & sincos_lat(j,2)*sincos_lon(k,2)*sincos_pole(2,2) )

        ENDDO
      ENDDO
    END IF
  END SUBROUTINE rotate_latlon_grid


  !-------------------------------------------------------------------------
  !> Compute normalized area weights for lon-lat grid
  !!
  !! @return 1D array (index: latitude)
  !!
  SUBROUTINE latlon_compute_area_weights( grid, earth_radius, area )
    
    TYPE (t_lon_lat_grid), INTENT(in)    :: grid
    REAL(wp),              INTENT(IN)    :: earth_radius
    REAL(wp),              INTENT(inout) :: area(:)
    ! local variables
    REAL(wp) :: start_lat, delta_lon, delta_lat, delta_lat_2,  &
      & pi_180, radius, rr_dlon, tot_area
    REAL(wp) :: latitude(grid%lat_dim)
    INTEGER :: k, pole1, pole2
    
    radius = earth_radius ! earth's radius (average)
    pi_180 = ATAN(1._wp)/45._wp
    start_lat   = grid%reg_lat_def(1) * pi_180
    delta_lon   = grid%reg_lon_def(2) * pi_180
    delta_lat   = grid%reg_lat_def(2) * pi_180
    delta_lat_2 = delta_lat / 2._wp
    rr_dlon     = radius*radius * delta_lon
    
    ! for each latitude, compute area of a grid box with a lon-lat
    ! point at its center
    DO k = 1, grid%lat_dim
      latitude(k)  = start_lat + REAL(k-1,wp)*delta_lat
      area(k)      = 2._wp*rr_dlon * SIN(delta_lat_2) * COS(latitude(k))
    END DO
    ! treat the special case of the poles (compute area of "triangle"
    ! with one vertex at the pole and the opposite side with constant
    ! latitude)
    pole1 = INT(( 90._wp - grid%reg_lat_def(1))/grid%reg_lat_def(2)) + 1
    pole2 = INT((-90._wp - grid%reg_lat_def(1))/grid%reg_lat_def(2)) + 1
    IF ((pole1 > 0) .AND. (pole1 <= grid%lat_dim)) THEN
      area(pole1) = rr_dlon*(1._wp-SIN(latitude(pole1)-delta_lat_2))
    END IF
    IF ((pole2 > 0) .AND. (pole2 <= grid%lat_dim)) THEN
      area(pole2) = rr_dlon*(1._wp-ABS(SIN(latitude(pole2)+delta_lat_2)))
    END IF
    ! normalize
    tot_area = REAL(grid%lon_dim,wp) * SUM(area)
    DO k = 1, grid%lat_dim
      area(k) = area(k)/tot_area
    END DO
  END SUBROUTINE latlon_compute_area_weights


  !---------------------------------------------------------------
  ! Test two floating point arrays for equality.
  !
  FUNCTION float_cmp(arr1, arr2, zero_tol)
    LOGICAL :: float_cmp
    REAL(wp), INTENT(IN) :: arr1(:), arr2(:)
    REAL(wp), INTENT(IN) :: zero_tol
    ! local variables
    INTEGER :: i

    float_cmp = .TRUE.
    IF (SIZE(arr1) /= SIZE(arr2)) THEN
      float_cmp = .FALSE.
      RETURN
    END IF
    DO i=1,SIZE(arr1)
      IF (ABS(arr1(i) - arr2(i)) > zero_tol) THEN
        float_cmp = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION float_cmp


  !---------------------------------------------------------------
  ! Test two lon-lat grid specifications for equality.
  !
  FUNCTION lonlat_grid_compare(grid1, grid2)
    LOGICAL :: lonlat_grid_compare
    TYPE(t_lon_lat_grid), INTENT(IN) :: grid1, grid2
    ! local variables
    REAL(wp), PARAMETER :: ZERO_TOL = 1.e-15_wp

    lonlat_grid_compare = .TRUE.
    IF  ( (grid1%total_dim    /= grid2%total_dim)                            .OR.  &
      &   (grid1%nblks        /= grid2%nblks    )                            .OR.  &
      &   (grid1%npromz       /= grid2%npromz   )                            .OR.  &
      &   (grid1%lon_dim      /= grid2%lon_dim  )                            .OR.  &
      &   (grid1%lat_dim      /= grid2%lat_dim  )                            .OR.  &
      &   (grid1%reg_def_mode /= grid2%reg_def_mode  )                       .OR.  &
      &   (.NOT. float_cmp(grid1%reg_lon_def, grid2%reg_lon_def , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%reg_lat_def, grid2%reg_lat_def , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%north_pole,  grid2%north_pole  , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%delta,       grid2%delta       , ZERO_TOL)) .OR.  &
      &   (.NOT. float_cmp(grid1%start_corner,grid2%start_corner, ZERO_TOL)) )  THEN
      lonlat_grid_compare = .FALSE.
    END IF
  END FUNCTION lonlat_grid_compare

END MODULE mo_lonlat_grid
