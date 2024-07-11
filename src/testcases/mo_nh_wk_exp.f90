! Subroutine to initialized the Weisman Klemp test case
! for the NH-Core in limited area mode
!
! Literature
! -M. L. Weisman and J. B. Klemp, 1982
!  The Dependence of Numerically Simulated Convective Storms on
!  Vertical Wind Shear and Buoyancy.
!  Monthly Weather Review, 110,504-520
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

MODULE mo_nh_wk_exp

!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!
  
   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd_o_cpd, p0ref, grav, &
     &                                cvd_o_rd, cpd ,     &
     &                                vtmpc1 , rdv,  rd,             &
     &                                cp_d => cpd, &
     &                                earth_radius
   USE mo_math_types,           ONLY: t_cartesian_coordinates, t_geographical_coordinates
   USE mo_math_constants,       ONLY: pi, deg2rad
   USE mo_math_utilities,       ONLY: gc2cc, arc_length, plane_torus_distance
   USE mo_grid_geometry_info,   ONLY: planar_torus_geometry, sphere_geometry
   USE mo_model_domain,         ONLY: t_patch
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_run_config,           ONLY: iqv
   USE mo_impl_constants,       ONLY: min_rlcell_int
   USE mo_parallel_config,      ONLY: nproma
   USE mo_satad,                ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
     &                                sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
     &                                spec_humi          !! Specific humidity
   USE mo_exception,            ONLY: message, finish, message_text
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_loopindices,          ONLY: get_indices_e
   USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
   USE mo_extpar_config,        ONLY: itopo
   USE mo_sync,                 ONLY: sync_patch_array, SYNC_C
   USE mo_nh_init_utils,        ONLY: init_w
   USE mo_hydro_adjust,         ONLY: hydro_adjust
   USE mo_vertical_coord_table, ONLY: vct_a

   IMPLICIT NONE

   PUBLIC  :: init_nh_topo_wk, init_nh_env_wk, init_nh_buble_wk
  
   PRIVATE

   CHARACTER(len=*), PARAMETER :: modname = 'mo_nh_wk_exp'
   
   REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd
   REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd


! !DEFINED PARAMETERS for Weisman Klemp:
   REAL(wp), PARAMETER :: hmin_wk         = 0.0_wp         ! [m], base height of the profile
   REAL(wp), PARAMETER :: p_base_wk       = 1.0e5_wp       ! [Pa], pressure at height hmin_wk

   REAL(wp), PARAMETER :: h_tropo_wk      = 12000._wp      ! [m], height of the tropopause
   REAL(wp), PARAMETER :: theta_0_wk      = 300._wp        ! [K], pot. temp. at z=0
   REAL(wp), PARAMETER :: theta_tropo_wk  = 343._wp        ! [K], pot. temp. at z=h_tropo_wk
   REAL(wp), PARAMETER :: expo_theta_wk   = 1.25_wp
   REAL(wp), PARAMETER :: expo_relhum_wk  = 1.25_wp
   REAL(wp), PARAMETER :: t_tropo_wk      = 213.0_wp


  ! Data for the initial profile of relative humidity:
   REAL(wp), PARAMETER :: rh_min_wk   = 0.10_wp        ! [%] ! rel. hum. above the tropopause level
   REAL(wp), PARAMETER :: rh_max_wk   = 0.95_wp        ! [%]  

   INTEGER,  PARAMETER :: niter=20
   REAL(wp), PUBLIC :: qv_max_wk   != 0.012_wp  - 0.016_wp  !NAMELIST PARAMETER
  ! Data for the wind profile
   REAL(wp), PUBLIC :: u_infty_wk  != 0 - 45  ms-1          !NAMELIST PARAMETER
   REAL(wp), PARAMETER :: href_wk  = 3000._wp ! Scaling height (height of 70 % windspeed) 
                                    ! for the Weisman-Klemp wind profile [m]
  ! Data for the thermal perturbation

   REAL(wp), PUBLIC :: bubctr_lat, bubctr_lon ! buble position in degrees   !NAMELIST PARAMETER
   REAL(wp), PUBLIC :: bubctr_z               ! buble position in meters    !NAMELIST PARAMETER
   REAL(wp), PUBLIC :: bub_hor_width, bub_ver_width !buble width in meters  !NAMELIST PARAMETER
   REAL(wp), PUBLIC :: bub_amp                !buble amplitude in Kelvin    !NAMELIST PARAMETER
  
!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of topography for the Weisman Klemp test case
  !! The topography is set to 0, but the possibility  to have 
  !!  a different topography is open
  !!
  SUBROUTINE init_nh_topo_wk( topo_c, nblks_c, npromz_c )

    INTEGER,  INTENT (IN) ::  nblks_c, npromz_c
    REAL(wp), INTENT(INOUT) :: topo_c    (nproma,nblks_c)

    ! local variables

  INTEGER        :: jc, jv, jb, nlen
!--------------------------------------------------------------------

      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  npromz_c
        ENDIF
        DO jc = 1, nlen
          IF ( itopo==0 ) THEN
           topo_c(jc,jb) = 0._wp
          END IF
        ENDDO
      ENDDO


  END SUBROUTINE init_nh_topo_wk
!------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the Weisman Klemp  test case 
  !! 
  !!
  SUBROUTINE init_nh_env_wk( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                p_metrics, p_int, l_hydro_adjust )

    TYPE(t_patch), TARGET, INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    TYPE(t_int_state), INTENT(IN)       :: p_int
    LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition

    INTEGER        ::  jc, jb, jk, je,   &
                       nlen, nblks_e,  nblks_c, npromz_c
    INTEGER        :: i_startidx, i_endidx, i_startblk
    INTEGER        :: nlev        !< number of full levels
    INTEGER        :: k_tropo

    REAL(wp)       :: z_u, exner_tropo, e_tropo, pres_tropo, qv_tropo, theta_v_tropo,    &
                      exner_aux, temp_aux, e_aux, ew_aux, ei_aux ,pres_aux, qv_aux, theta_v_aux, qv_extrap

    REAL(wp), DIMENSION(ptr_patch%nlev) :: z_full, theta, exner, pres, qv, theta_v, rh, &
                                           temp, rhwcheck, rhicheck

!--------------------------------------------------------------------
!

    nblks_c   = ptr_patch%nblks_c
    npromz_c  = ptr_patch%npromz_c
    nblks_e   = ptr_patch%nblks_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev

    ! height of main levels    
    DO jk = 1, nlev
      z_full(jk) = 0.5_wp*(vct_a(jk)+vct_a(jk+1))
    ENDDO

    ! Determine level index right above the tropopause
    DO jk = 1, nlev-1
      IF (z_full(jk) >= h_tropo_wk .AND. z_full(jk+1) < h_tropo_wk) THEN
        k_tropo = jk
        EXIT
      ENDIF
      IF (jk == nlev-1) CALL finish ('WK initialization', &
         'model top must be higher than 12 km')
    ENDDO

    ! Tropopause parameters
    exner_tropo   = t_tropo_wk/theta_tropo_wk
! UB    e_tropo       = rh_min_wk*sat_pres_ice(t_tropo_wk)
    e_tropo       = rh_min_wk*sat_pres_water(t_tropo_wk)
    pres_tropo    = p0ref*(exner_tropo**cpd_o_rd)
    qv_tropo      = spec_humi(e_tropo,pres_tropo)
    theta_v_tropo = theta_tropo_wk*(1._wp+vtmpc1*qv_tropo)

    ! profiles above the tropopause: note that T = const here (and therefore e = e_tropo)
    DO jk = 1, k_tropo
      theta(jk) = theta_tropo_wk*EXP(grav_o_cpd/t_tropo_wk*(z_full(jk)-h_tropo_wk))
      exner(jk) = exner_tropo*EXP(-grav_o_cpd/t_tropo_wk*(z_full(jk)-h_tropo_wk))
      pres(jk)  = p0ref*(exner(jk)**cpd_o_rd)
      qv(jk)    = spec_humi(e_tropo,pres(jk))
      theta_v(jk) = theta(jk)*(1._wp+vtmpc1*qv(jk))
      rh(jk)    = rh_min_wk
      temp(jk)  = t_tropo_wk
    ENDDO

    ! known quantities below the tropopause
    DO jk = k_tropo+1, nlev
      theta(jk) = theta_0_wk+(theta_tropo_wk-theta_0_wk)*(z_full(jk)/h_tropo_wk)**expo_theta_wk
      rh(jk)    = 1._wp-0.75_wp*(z_full(jk)/h_tropo_wk)**expo_relhum_wk
      rh(jk)    = min(rh(jk),rh_max_wk)
    ENDDO

    ! Piecewise vertical integration from the TP towards the bottom
    jk = k_tropo+1

    ! 1st step: preliminary estimate
    theta_v_aux = theta(jk)*(1._wp+vtmpc1*qv_tropo)
    exner_aux   = exner_tropo-grav_o_cpd*(z_full(jk)-h_tropo_wk)/(theta_v_aux-theta_v_tropo)*&
                  LOG(theta_v_aux/theta_v_tropo)
    temp_aux    = theta(jk)*exner_aux
! UB    IF (temp_aux > tmelt) THEN
      e_aux     = rh(jk)*sat_pres_water(temp_aux)
! UB    ELSE
! UB      e_aux     = rh(jk)*sat_pres_ice(temp_aux)
! UB    ENDIF
    pres_aux    = p0ref*(exner_aux**cpd_o_rd)
    qv_aux      = spec_humi(e_aux,pres_aux)
    theta_v_aux = theta(jk)*(1._wp+vtmpc1*qv_aux) 

    ! 2nd step: final computation  
    exner(jk)   = exner_tropo-grav_o_cpd*(z_full(jk)-h_tropo_wk)/(theta_v_aux-theta_v_tropo)*&
                  LOG(theta_v_aux/theta_v_tropo)
    temp(jk)    = theta(jk)*exner(jk)
! UB    IF (temp(jk) > tmelt) THEN
      e_aux     = rh(jk)*sat_pres_water(temp(jk))
! UB    ELSE
! UB      e_aux     = rh(jk)*sat_pres_ice(temp(jk))
! UB    ENDIF
    pres(jk)    = p0ref*(exner(jk)**cpd_o_rd)
    qv(jk)      = spec_humi(e_aux,pres(jk))
    theta_v(jk) = theta(jk)*(1._wp+vtmpc1*qv(jk)) 

    ! Remaining model layers
    DO jk = k_tropo+2, nlev

      ! 1st step: preliminary estimate
!>FR
      qv_extrap   = MIN(qv_max_wk,qv(jk-1)+(qv(jk-2)-qv(jk-1))/          &
      (z_full(jk-2)-z_full(jk-1))*(z_full(jk)-z_full(jk-1)) )
!      qv_extrap   = qv(jk-1)+(qv(jk-2)-qv(jk-1))/          &
!        (z_full(jk-2)-z_full(jk-1))*(z_full(jk)-z_full(jk-1)) 
!<FR
      theta_v_aux = theta(jk)*(1._wp+vtmpc1*qv_extrap)
      exner_aux   = exner(jk-1)-grav_o_cpd*(z_full(jk)-z_full(jk-1))/ &
                   (theta_v_aux-theta_v(jk-1))*LOG(theta_v_aux/theta_v(jk-1))
      temp_aux    = theta(jk)*exner_aux
! UB      IF (temp_aux > tmelt) THEN
        e_aux     = rh(jk)*sat_pres_water(temp_aux)
! UB      ELSE
! UB        e_aux     = rh(jk)*sat_pres_ice(temp_aux)
! UB      ENDIF
      pres_aux    = p0ref*(exner_aux**cpd_o_rd)
!>FR
      qv_aux      = MIN(qv_max_wk,spec_humi(e_aux,pres_aux))
!      qv_aux      = spec_humi(e_aux,pres_aux)
!<FR
      theta_v_aux = theta(jk)*(1._wp+vtmpc1*qv_aux) 

      ! 2nd step: final computation  
      exner(jk)   = exner(jk-1)-grav_o_cpd*(z_full(jk)-z_full(jk-1))/ &
                   (theta_v_aux-theta_v(jk-1))*LOG(theta_v_aux/theta_v(jk-1))
      temp(jk)    = theta(jk)*exner(jk)
! UB      IF (temp(jk) > tmelt) THEN
        e_aux     = rh(jk)*sat_pres_water(temp(jk))
! UB      ELSE
! UB       e_aux     = rh(jk)*sat_pres_ice(temp(jk))
! UB      ENDIF
      pres(jk)    = p0ref*(exner(jk)**cpd_o_rd)
!>FR
      qv(jk)      = MIN(qv_max_wk,spec_humi(e_aux,pres(jk)))
!      qv(jk)      = spec_humi(e_aux,pres(jk))
!<FR
      theta_v(jk) = theta(jk)*(1._wp+vtmpc1*qv(jk)) 


    ENDDO

    ! check on rh
    WRITE(message_text,'(a)') ' Check on initial WK humidity profile '
    CALL message('', TRIM(message_text))
    DO jk = 1, nlev
      ! recalculate RH from final atmospheric profile
      ew_aux = sat_pres_water(temp(jk))
      ei_aux = sat_pres_ice(temp(jk))
      rhwcheck(jk) = qv(jk)*pres(jk)/(rdv - (1._wp-rdv)*qv(jk))/ew_aux
      rhicheck(jk) = qv(jk)*pres(jk)/(rdv - (1._wp-rdv)*qv(jk))/ei_aux
      WRITE(message_text,'(a,i4,3(a,f7.2))') 'level ',jk,': analytical RH(%)',rh(jk)*100._wp, &
        ', actual RH(%) ',rhwcheck(jk)*100._wp,', actual RHi(%) ',rhicheck(jk)*100._wp
      CALL message('', TRIM(message_text))
    ENDDO

   ! Copy to prognostic model fields
!$OMP PARALLEL PRIVATE(i_startblk)
!$OMP DO PRIVATE(jb,jk,jc,nlen)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
        DO jc = 1, nlen  
          ptr_nh_prog%theta_v(jc,jk,jb)    = theta_v(jk)
          ptr_nh_prog%exner(jc,jk,jb)      = exner(jk)
          ptr_nh_prog%tracer(jc,jk,jb,iqv) = qv(jk)
          ptr_nh_prog%rho(jc,jk,jb)  = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref &
                                       /rd/ptr_nh_prog%theta_v(jc,jk,jb)
        ENDDO !jc
      ENDDO !jk     
    ENDDO !jb
!$OMP END DO

! initialize horizontal velocities
    i_startblk = ptr_patch%edges%start_blk(2,1)

    ! horizontal normal components of the velocity
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_u)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_u = u_infty_wk * ( TANH((z_full(jk)-hmin_wk)/(href_wk-hmin_wk)) - 0.45_wp) 
            ptr_nh_prog%vn(je,jk,jb) = &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO !je
        ENDDO !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    CALL diagnose_pres_temp (p_metrics, ptr_nh_prog,ptr_nh_prog, ptr_nh_diag,     &
                             ptr_patch, opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)

! initialize vertical velocity
   CALL init_w(ptr_patch, p_int, ptr_nh_prog%vn, p_metrics%z_ifc, ptr_nh_prog%w)
   CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)


  IF (l_hydro_adjust) THEN

   CALL hydro_adjust ( ptr_patch, p_metrics, ptr_nh_prog%rho,  &
                     & ptr_nh_prog%exner, ptr_nh_prog%theta_v  )

  END IF

  END SUBROUTINE init_nh_env_wk
!--------------------------------------------------------------------
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of thermal buble  for the Weisman Klemp  test case 
  !! 
  !!
  SUBROUTINE init_nh_buble_wk( ptr_patch, p_metrics, ptr_nh_prog, ptr_nh_diag )

  
    TYPE(t_patch), TARGET, INTENT(INOUT):: &  !< patch on which computation is performed
         &  ptr_patch

    TYPE(t_nh_metrics), INTENT(IN)     :: p_metrics !< NH metrics state
    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
         &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
         &  ptr_nh_diag
 ! local variables  

    INTEGER        :: jc, jk, jb, nlen, nlev
    INTEGER        :: nblks_c, npromz_c

    REAL(wp) :: x_loc(3), x_c(3), x_bubble(3), dis, dis_hor, dis_vert
    TYPE(t_geographical_coordinates) :: geo_bub
    TYPE(t_cartesian_coordinates)  :: cart_loc, cart_bub

    CHARACTER(len=250) :: message_text
!--------------------------------------------------------------------

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nblks_c   = ptr_patch%nblks_c
    npromz_c  = ptr_patch%npromz_c
    
    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)

      message_text(:) = ' '
      WRITE (message_text, '(2(a,f0.5),a,f0.1,a,f0.2)') 'wk82: bubble on torus domain lon=', bubctr_lon, &
           ' lat=', bubctr_lat, ' height=', bubctr_z, ' ampl=', bub_amp
      CALL message(modname//': init_nh_buble_wk', TRIM(message_text))

      ! Bubble center : valid on the torus domain
      x_bubble = (/ bubctr_lon, bubctr_lat, bubctr_z /)

      ! Non-dimensionalize the bubble center
      x_c(1) = x_bubble(1) / bub_hor_width
      x_c(2) = x_bubble(2) / bub_hor_width
      x_c(3) = x_bubble(3) / bub_ver_width

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,x_loc,dis)
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        
        DO jk = 1 , nlev
          DO jc = 1 , nlen
            x_loc(1) = ptr_patch%cells%cartesian_center(jc,jb)%x(1)/bub_hor_width
            x_loc(2) = ptr_patch%cells%cartesian_center(jc,jb)%x(2)/bub_hor_width
            x_loc(3) = p_metrics%z_mc(jc,jk,jb)/bub_ver_width
            dis = plane_torus_distance(x_loc,x_c,ptr_patch%geometry_info)
            IF(dis < 1._wp)THEN
              ptr_nh_prog%theta_v(jc,jk,jb) = ptr_nh_prog%theta_v(jc,jk,jb) + &
                   bub_amp * COS(dis*pi/2._wp )**2  * (1.0_wp + vtmpc1*ptr_nh_prog%tracer(jc,jk,jb,iqv))
              ptr_nh_prog%rho(jc,jk,jb) = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd * &
                   p0ref / (rd * ptr_nh_prog%theta_v(jc,jk,jb) )
            END IF
          END DO !jk
        END DO !jc
      ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    CASE (sphere_geometry)

      message_text(:) = ' '
      WRITE (message_text, '(2(a,f0.5),a,f0.1,a,f0.2)') 'wk82: bubble on spherical domain lon=', bubctr_lon, &
           ' lat=', bubctr_lat, ' height=', bubctr_z, ' ampl=', bub_amp
      CALL message(modname//': init_nh_buble_wk', TRIM(message_text))

      geo_bub % lon = bubctr_lon * deg2rad
      geo_bub % lat = bubctr_lat * deg2rad
      cart_bub = gc2cc( geo_bub )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,x_loc,dis)
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        
        DO jk = 1 , nlev
          DO jc = 1 , nlen
            cart_loc = gc2cc( ptr_patch%cells%center(jc,jb) )
            dis_hor  = arc_length(cart_bub, cart_loc, ptr_patch%geometry_info) * earth_radius / bub_hor_width
            dis_vert = ( p_metrics%z_mc(jc,jk,jb) - bubctr_z ) / bub_ver_width
            dis = SQRT( dis_hor*dis_hor + dis_vert*dis_vert )
            IF(dis < 1._wp)THEN
              ptr_nh_prog%theta_v(jc,jk,jb) = ptr_nh_prog%theta_v(jc,jk,jb) +     &
                   bub_amp * COS(dis*pi/2._wp )**2  * (1.0_wp + vtmpc1*ptr_nh_prog%tracer(jc,jk,jb,iqv))
              ptr_nh_prog%rho(jc,jk,jb) = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd * &
                   p0ref / (rd * ptr_nh_prog%theta_v(jc,jk,jb))
            END IF
          END DO !jk
        END DO !jc
      ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL
      
    CASE DEFAULT
      CALL finish(modname//': init_nh_buble_wk', "Undefined geometry type for wk82 bubbles")
    END SELECT
      
    CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog,     &
                              ptr_nh_prog, ptr_nh_diag,   &
                              ptr_patch,                  &
                              opt_calc_temp=.TRUE.,       &
                              opt_calc_pres=.FALSE.,       &
                              opt_rlend=min_rlcell_int )
     
  END SUBROUTINE init_nh_buble_wk
!--------------------------------------------------------------------
!!$! Function taken from COSMO
!!$
!!$  ! Specific humidity as function of T, p, and relhum 
!!$  !   (and qcrs = sum of all hydrometeor contents):
!!$  ! NOTE: on input, relhum has to be smaller than p / E(T)!
!!$  REAL(wp) FUNCTION qv_Tprelhum(p, temp, relhum, qcrs)
!!$    IMPLICIT NONE
!!$    REAL(wp), INTENT(IN) :: p, temp, relhum, qcrs
!!$    REAL(wp) :: zesat_w, coeff
!!$    !zesat_w  =b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
!!$    zesat_w  = sat_pres_water(temp)
!!$    coeff = relhum * zesat_w * rdv / p
!!$    qv_Tprelhum = coeff * (1.0_wp + qcrs) / (1.0_wp - rvd_m_o*coeff)
!!$  END FUNCTION qv_Tprelhum
!--------------------------------------------------------------------
  ! Gas constant of moist air containing hydrometeors (qcrs is the sum of
  ! the specific hydrometeor contents):
  REAL(wp) FUNCTION rd_moist(qv,qcrs)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: qv, qcrs

    rd_moist = rd * (1.0_wp + vtmpc1*qv - qcrs)

  END FUNCTION rd_moist
!--------------------------------------------------------------------
  ! COSMO-APPROXIMATION: CP IS APPROXIMATED TO BE THAT OF DRY AIR
  REAL(wp) FUNCTION cp_moist_cosmo(qv,ql,qi)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: qv, ql, qi
    REAL(wp)             :: qd

    ! Dry air content:
    qd = 1.0_wp - qv - ql -qi

    ! cp:
    cp_moist_cosmo = qd*cp_d + qv*cp_d + ql*cp_d + qi*cp_d

  END FUNCTION cp_moist_cosmo
!--------------------------------------------------------------------
  END MODULE mo_nh_wk_exp

