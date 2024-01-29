!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"
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

MODULE mo_emvorado_warmbubbles

  USE mo_kind,                      ONLY: wp
  USE mo_emvorado_warmbubbles_type, ONLY: t_bubblecontainer, t_warmbubble
  USE mo_parallel_config,           ONLY: nproma
  USE mo_impl_constants,            ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices,               ONLY: get_indices_c
  USE mo_model_domain,              ONLY: t_patch
  USE mo_nonhydro_types,            ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_run_config,                ONLY: iqv
  USE mo_physical_constants,        ONLY: earth_radius
  USE mo_math_types,                ONLY: t_geographical_coordinates
  USE mo_math_constants,            ONLY: pi, deg2rad
  USE mo_math_utilities,            ONLY: plane_torus_distance
  USE mo_grid_geometry_info,        ONLY: planar_torus_geometry, sphere_geometry
  USE mo_satad,                     ONLY: sat_pres_water
  USE mo_exception,                 ONLY: message, finish

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: set_artif_heatrate_dist
  
CONTAINS

  !=======================================================================================
  !
  ! "Warm bubbles": isochoric heating disturbances, optionally constant RH during heating
  !
  ! Method: - warm bubbles are requested on input by the list of bubbles in bublist
  !         - T-disturbances are added to the model temperature in p_diag%temp
  !         - optionally, QV is increased accordingly to keep RH constant at constant
  !           total density
  !         - The shape of the disturbances is 3D elliptic in terrain-following coordinates
  !
  !=======================================================================================
 
  SUBROUTINE set_artif_heatrate_dist(jg, sim_time, bublist, dt, p_patch, p_metrics, p_prog_rcf, p_diag)

    INTEGER                , intent(in)    :: jg         ! Domain index to decorate debug output
    REAL(wp)               , INTENT(in)    :: sim_time   ! simulation time in sec. since tc_exp_startdate
    TYPE(t_bubblecontainer), INTENT(INOUT) :: bublist    ! INOUT due to possible testmode=.TRUE.
    REAL(wp)               , INTENT(in)    :: dt         ! physics time step for this domain
    TYPE(t_patch)          , INTENT(IN)    :: p_patch
    TYPE(t_nh_metrics)     , INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog)        , INTENT(INOUT) :: p_prog_rcf ! contains tracers
    TYPE(t_nh_diag)        , INTENT(INOUT) :: p_diag     ! contains temperature

    CHARACTER(len=*), PARAMETER :: routine="mo_emvorado_warmbubbles::set_artif_heatrate_dist"
    CHARACTER(len=250)          :: message_text

    INTEGER  :: i_rlstart, i_rlend, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i, jb, jk, jc
    REAL(wp) :: f_xyz_3D(nproma, p_patch%nlev, p_patch%nblks_c), tnew, told, tinc
    TYPE(t_warmbubble), POINTER    :: p_bub
    LOGICAL :: bub_is_active

    LOGICAL, PARAMETER :: testmode = .FALSE.

    ! For testing, define 30 test bubbles scattered over central Germany and the first 1.5 h simulation time:
    IF (testmode .AND. .NOT.ASSOCIATED(bublist%bubs)) THEN
      bublist%num_bubs = 30
      ALLOCATE(bublist%bubs(bublist%num_bubs))
      bublist%bubs(1:30)%timecounter       = 0.0                           ! Actual duration during which the bubble was already active [s]
      bublist%bubs(1:30)%ltempdist         = .TRUE.                        ! Switch to set temperature disturbance to ACTIVE

      ! 3 Prototypes, which are cloned to 30 bubbles below:
      bublist%bubs(1:3)%ctype_tempdist     = (/'cos-instant ','cos-instant ','cos-hrd     '/) ! Type of perturbation 'cos-hrd', 'cos-instant'
      bublist%bubs(1:3)%ladd_bubblenoise_t = (/.FALSE.,.FALSE.,.FALSE./)   ! Switch to overlay random noise on the disturbance (not yet implemented)
      bublist%bubs(1:3)%lbub_rhconst       = (/.TRUE.,.TRUE.,.TRUE./)      ! Switch to activate a moisture increment such that rel. hum. stays constant during heating 
      bublist%bubs(1:3)%htempdist          = (/120.,280.,420./)            ! Time for beginning of temperature disturbance since model start time [s]
      bublist%bubs(1:3)%centlon            = (/9.55,9.8,10.35/)            ! Center (lon) of temperature disturbance [deg]
      bublist%bubs(1:3)%centlat            = (/49.57,50.34,50.58/)         ! Center (lat) of temperature disturbance [deg]
      bublist%bubs(1:3)%centz              = (/1500.0,2500.0,1700.0/)      ! Center (Z) of temperature disturbance [m AGL]
      bublist%bubs(1:3)%timespan           = (/0.0,0.0,340.0/)             ! Total duration for release of temperature disturbance [s]
      bublist%bubs(1:3)%radx               = (/10000.0,20000.,25000.0/)    ! Horizontal main axis (radius) in X (longitudinal) direction of temperature disturbances [m]
      bublist%bubs(1:3)%rady               = (/20000.0,10000.0,15000.0/)   ! Horizontal main axis (radius) in Y (latitudinal) direction of temperature disturbances [m]
      bublist%bubs(1:3)%radz               = (/1500.0,2500.0,2500.0/)      ! Vertical main axis (radius) of temperature disturbances [m]
      bublist%bubs(1:3)%rotangle           = (/0.,45.0,-45.0/)             ! Rotation angle of main axes of temperature disturbances [degrees]
      bublist%bubs(1:3)%heatingrate        = (/0.0,0.0,.0238/)             ! Constant heating rate for 'cos-hrd' bubble [K/s]
      bublist%bubs(1:3)%dT                 = (/5.0,4.0,3.0/)               ! Temperature increment of 'cos-instant' bubble [K]
      bublist%bubs(1:3)%dT_bubblenoise     = (/0.1,0.2,0.3/)               ! In case of ladd_bubblenoise_t=.true., relative noise level, such that
                                                                           !   dT          = dT          * (1 + dT_bubblenoise * random_noise[-1,1])   ('cos-instant')
                                                                           !   heatingrate = heatingrate * (1 + dT_bubblenoise * random_noise[-1,1])   ('cos-hrd')
      ! Clone the 3 prototypes and shift in space and time:
      DO i=4,28,3
        bublist%bubs(i:i+2)           = bublist%bubs(1:3)
        bublist%bubs(i:i+2)%htempdist = bublist%bubs(1:3)%htempdist + (i-3)*200.0
        bublist%bubs(i:i+2)%centlon   = bublist%bubs(1:3)%centlon + (i-20)*0.18
        bublist%bubs(i:i+2)%centlat   = bublist%bubs(1:3)%centlat + (20-i)*0.2
      END DO
    END IF

    !======================================================
    !
    ! update p_prog_rcf%tracer(iqv) and p_diag%temp:
    !
    !======================================================
    
    IF (bublist%num_bubs > 0) THEN

      DO i=1, bublist%num_bubs

        p_bub => bublist%bubs(i)

        ! bubble is active, if the actual simulation time is within [p_bub%htempdist; p_bub%htempdist + p_bub%timespan - dt]
        !  so that the heating rate is applied exactly for the time span p_bub%timespan.
        bub_is_active = ( TRIM(p_bub%ctype_tempdist) == 'cos-hrd' .AND. &
             sim_time >= p_bub%htempdist-0.5_wp*dt .AND. sim_time < p_bub%htempdist+(CEILING(p_bub%timespan/dt)+0.5_wp)*dt ) &
                        .OR. &
                        ( TRIM(p_bub%ctype_tempdist) == 'cos-instant' .AND. &
                            sim_time >= p_bub%htempdist-0.5_wp*dt .AND. sim_time < p_bub%htempdist+0.5_wp*dt )
        
        IF (p_bub%ltempdist .AND. bub_is_active) THEN

          SELECT CASE (TRIM(p_bub%ctype_tempdist))

          CASE ('cos-instant','cos-hrd')
          
            ! Compute spatial pattern of the bubble in the range [0,1]:
            CALL calc_f_xyz_cos(p_bub, p_patch, p_metrics, f_xyz_3D)

            message_text(:) = ' '
            ! Max. temperature increment for the bubble:
            IF (TRIM(p_bub%ctype_tempdist) == 'cos-instant') THEN
              tinc = p_bub%dT
              p_bub%ltempdist = .FALSE.  ! deactivate bubble for the next timesteps
              WRITE (message_text, '(a,i3,2(a,f0.5),a,f0.1,a,f0.2,a,f0.1)') &
                   'automatic warm bubble on domain=', jg, ' type='//TRIM(p_bub%ctype_tempdist)//' lon=', p_bub%centlon, &
                   ' lat=', p_bub%centlat, ' height[mAGL]=', p_bub%centz, ' dT[K]=', p_bub%dT, &
                   ' time[s]=', sim_time
            ELSE
              p_bub%timecounter = p_bub%timecounter + dt
              IF (p_bub%timecounter > p_bub%timespan) THEN
                tinc = p_bub%heatingrate * (dt - MODULO(p_bub%timecounter, p_bub%timespan))
                p_bub%timecounter = p_bub%timespan
              ELSE
                tinc = p_bub%heatingrate * dt
              END IF
              WRITE (message_text, '(a,i3,2(a,f0.5),a,f0.1,a,f0.2,a,f0.1,a,f0.1)') &
                   'automatic warm bubble on domain=', jg, ' type='//TRIM(p_bub%ctype_tempdist)//' lon=', p_bub%centlon, &
                   ' lat=', p_bub%centlat, ' height[mAGL]=', p_bub%centz, ' heatingrate[K/s]=', p_bub%heatingrate, &
                   ' time[s]=', sim_time,' heattime[s]=', p_bub%timecounter 
            END IF
            CALL message(TRIM(routine), TRIM(message_text))
      
            ! without halo or boundary  points:
            i_rlstart = grf_bdywidth_c + 1
            i_rlend   = min_rlcell_int

            i_startblk = p_patch%cells%start_block( i_rlstart )
            i_endblk   = p_patch%cells%end_block  ( i_rlend   )
            
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,told,tnew)
            DO jb = i_startblk, i_endblk

              CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
                                  i_startidx, i_endidx, i_rlstart, i_rlend)

              DO jk = 1, p_patch%nlev
                DO jc = i_startidx, i_endidx

                  IF (f_xyz_3D(jc,jk,jb) > 1E-20_wp) THEN
                    ! T-increment:
                    told = p_diag%temp(jc,jk,jb)
                    p_diag%temp(jc,jk,jb) = told + tinc*f_xyz_3D(jc,jk,jb)
                    IF (p_bub%lbub_rhconst) THEN
                      tnew = p_diag%temp(jc,jk,jb)
                      ! qv-increment assuming constant RH at constant total density:
                      p_prog_rcf%tracer(jc,jk,jb,iqv) = p_prog_rcf%tracer(jc,jk,jb,iqv) * &
                           told * sat_pres_water(tnew) / ( tnew * sat_pres_water(told) )
                    END IF
                  END IF
                  
                END DO
              END DO
              
            END DO
!$OMP END DO
!$OMP END PARALLEL

          CASE default
            CALL finish(TRIM(routine), "Undefined bubble type for automatic warm bubbles! "// &
                        "Possible are ''cos-instant'' or ''cos-hrd''")
          END SELECT


        END IF

      END DO
      
    END IF
    
  END SUBROUTINE set_artif_heatrate_dist


  SUBROUTINE calc_f_xyz_cos(bub, p_patch, p_metrics, f_xyz_3D)

    TYPE(t_warmbubble)     , INTENT(IN)    :: bub
    TYPE(t_patch)          , INTENT(IN)    :: p_patch
    TYPE(t_nh_metrics)     , INTENT(IN)    :: p_metrics
    REAL(wp)               , INTENT(OUT)   :: f_xyz_3D(:,:,:)

    INTEGER  :: i_rlstart, i_rlend, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: jb, jk, jc
    REAL(wp) :: hsurf(nproma), zdx_rot(nproma), zdy_rot(nproma)
    REAL(wp) :: x_loc_1, x_loc_2, x_loc_3, x_c(3), x_bubble(3), dist, bub_hor_width

    CHARACTER(len=*), PARAMETER :: routine="mo_emvorado_warmbubbles::calc_f_xyz_3d"

    TYPE(t_geographical_coordinates) :: geo_bub


    f_xyz_3D(:,:,:) = 0.0_wp
    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = p_patch%cells%start_block( i_rlstart )
    i_endblk   = p_patch%cells%end_block  ( i_rlend   )

    
    SELECT CASE(p_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)

!!!!!!!!!!!! NOT YET TESTED! CHECK BUBBLE GEOMETRY AND COORDINATES !!!!!!!!!!!!

      CALL finish(TRIM(routine), "Untested grid geometry type for automatic warm bubbles: planar_torus_geometry! STOP!")

      ! Bubble center : valid on the torus domain
      x_bubble = (/ bub%centlon*deg2rad*earth_radius, bub%centlat*deg2rad*earth_radius, bub%centz /)
      
      ! Non-dimensionalize the bubble center
      bub_hor_width = SQRT(bub%radx*bub%rady)  ! no ellipse, just an equal-area circle
      x_c(1) = x_bubble(1) / bub_hor_width
      x_c(2) = x_bubble(2) / bub_hor_width
      x_c(3) = x_bubble(3) / bub%radz

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,hsurf,dist,x_loc_1,x_loc_2,x_loc_3)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
             i_startidx, i_endidx, i_rlstart, i_rlend)
        
        DO jc = i_startidx, i_endidx
          hsurf(jc) = p_metrics%z_ifc(jc,p_patch%nlev+1,jb)
        END DO

        DO jk = 1, p_patch%nlev
          DO jc = i_startidx, i_endidx

            x_loc_1 = p_patch%cells%cartesian_center(jc,jb)%x(1) / bub_hor_width
            x_loc_2 = p_patch%cells%cartesian_center(jc,jb)%x(2) / bub_hor_width
            x_loc_3 = (p_metrics%z_mc(jc,jk,jb) - hsurf(jc)) / bub%radz
            dist = plane_torus_distance( (/ x_loc_1, x_loc_2, x_loc_3 /), x_c, p_patch%geometry_info)
            IF(dist < 1.0_wp)THEN
              f_xyz_3D(jc,jk,jb) = COS(0.5_wp*pi*dist)**2
            END IF
          
          END DO
      END DO
      
    END DO
!$OMP END DO
!$OMP END PARALLEL

      
    CASE (sphere_geometry)

      geo_bub%lon = bub%centlon * deg2rad
      geo_bub%lat = bub%centlat * deg2rad

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,hsurf,dist,zdx_rot,zdy_rot)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

!NEC$ ivdep
      DO jc = i_startidx, i_endidx
        hsurf(jc) = p_metrics%z_ifc(jc,p_patch%nlev+1,jb)
        CALL hill_rot_coords( p_patch%cells%center(jc,jb), geo_bub, &
                              bub%rotangle*deg2rad, 0.0_wp, zdx_rot(jc), zdy_rot(jc))
      END DO
      
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx

          ! Terrain-following ellsoidal bubble:
          dist  = SQRT( (zdx_rot(jc)/bub%radx)**2 + (zdy_rot(jc)/bub%rady)**2 + &
               (( p_metrics%z_mc(jc,jk,jb) - hsurf(jc) - bub%centz )/bub%radz)**2 ) ! terrain-following bubble
          IF (dist < 1.0_wp) THEN
            f_xyz_3d(jc,jk,jb) = COS(0.5_wp*pi*dist)**2
          END IF
          
        END DO
      END DO
        
    END DO
!$OMP END DO
!$OMP END PARALLEL

      
    CASE DEFAULT
      CALL finish(TRIM(routine), "Undefined grid geometry type for automatic warm bubbles! "// &
                                 "Possible are ''planar_torus_geometry'' or ''sphere_geometry''")
    END SELECT

  END SUBROUTINE calc_f_xyz_cos

!=================================================================================
!
! Functions for coordinate operations in spherical (geographical) coordinates:
!
!=================================================================================

  ! Distances relative to hill/bubble elliptic main axes (these can be rotated by <rotangle> 
  ! relative to the rotated North direction) for coordinates
  ! given at the mass points. The coordinates must be geographic coordinates in radians.
  ! 
  ! The main axes are along great circles and the distances are 
  ! also measured along great circles.
  ! 
  SUBROUTINE hill_rot_coords(geo_coord, bub_coord, rotangle, height, rx, ry)

    IMPLICIT NONE

    !.. Input/Output parameters:
    !---------------------------
    !   global index of model grid point:
    TYPE(t_geographical_coordinates), INTENT(in) :: geo_coord
    !   global index of center of bubble/hill:
    TYPE(t_geographical_coordinates), INTENT(in) :: bub_coord
    !   rotation angle of main hill/bubble y-axis clockwise relative to north in rad:
    REAL(KIND=wp),     INTENT(in)       :: rotangle
    !   height level where the arc lengths are referenced to in m:
    REAL(KIND=wp),     INTENT(in)       :: height
    !   arc lengths along great circles perpendicular to the 
    !   hill/bubble main axes (also great circles):
    REAL(KIND=wp),     INTENT(out)      :: rx, ry


    !.. Local variables:
    !---------------------------
    !   angles and arc lengths on the unit sphere in rad:
    REAL(KIND=wp)     ::  d, delta, &
         zlon, zlat, zlon_c, zlat_c, tmp, cos_arg

    zlon_c = bub_coord%lon
    zlat_c = bub_coord%lat
    zlon   = geo_coord%lon
    zlat   = geo_coord%lat    

    d     =  geo_distance ( bub_coord, geo_coord, 0.0_wp ) / earth_radius
    delta =  geo_course   ( bub_coord, geo_coord )

    !.. take rotated bubble orientation into account:
    delta = delta - rotangle

    rx = ASIN(SIN(d)*SIN(delta))
    tmp = 1.0_wp - SIN(delta)*SIN(rx)*SIN(d)
    IF (ABS(tmp) < 1.0E-20_wp) THEN
      ry = 0.0_wp
    ELSE
      cos_arg = COS(rx)*COS(d) / tmp
      cos_arg = MAX(MIN(cos_arg, 1.0_wp), -1.0_wp)
      ry = ACOS(cos_arg)
      IF (zlat_c > zlat) ry = -ry
    END IF
    rx = (earth_radius + height) * rx
    ry = (earth_radius + height) * ry

  END SUBROUTINE hill_rot_coords

  FUNCTION geo_distance(gc1, gc2, height) RESULT (dist)

    IMPLICIT NONE

    TYPE(t_geographical_coordinates), INTENT(in) :: gc1  ! in radians
    TYPE(t_geographical_coordinates), INTENT(in) :: gc2  ! in radians
    !   height level where the arc lengths are referenced to in m:
    REAL(KIND=wp),     INTENT(in)       :: height

    !   great circle distance between the two points in m:
    REAL(KIND=wp)                       :: dist

    dist = SIN(gc2%lat)*SIN(gc1%lat)+COS(gc2%lat)*COS(gc1%lat)*COS(gc2%lon-gc1%lon)
    dist = MAX(dist,-1.0_wp)
    dist = MIN(dist, 1.0_wp)
    dist = (earth_radius+height)*ACOS( dist)

  END FUNCTION geo_distance

  FUNCTION geo_course(gc1, gc2) RESULT (truecourse)

    IMPLICIT NONE

    TYPE(t_geographical_coordinates), INTENT(in) :: gc1  ! in radians
    TYPE(t_geographical_coordinates), INTENT(in) :: gc2  ! in radians
    
    !   geogr. direction ("course") from the start point gc1 to the target point gc2
    !   at the start point in rad:
    REAL(KIND=wp)     :: truecourse

    REAL(KIND=wp)     :: cangle, cos_cangle, cos_arg

    cos_cangle =  SIN(gc2%lat)*SIN(gc1%lat)+COS(gc2%lat)*COS(gc1%lat)*COS(gc2%lon-gc1%lon) 
    cangle = ACOS(cos_cangle)
    IF (ABS(cangle) < 1.0E-20_wp) cangle = 1.0E-20_wp

    cos_arg = (SIN(gc2%lat) - SIN(gc1%lat)*cos_cangle) / (COS(gc1%lat)*SIN(cangle))
    cos_arg = MAX(MIN(cos_arg, 1.0_wp), -1.0_wp)

    truecourse = ACOS(cos_arg)
    IF (gc1%lon > gc2%lon) truecourse = 2.0_wp*pi - truecourse

  END FUNCTION geo_course

END MODULE mo_emvorado_warmbubbles
